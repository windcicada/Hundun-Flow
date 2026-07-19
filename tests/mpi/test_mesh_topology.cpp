// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::mesh::BoundaryPatch;
using hundun::mesh::EntityOwnership;
using hundun::mesh::FaceAxis;
using hundun::mesh::GlobalCellId;
using hundun::mesh::GlobalFaceId;
using hundun::mesh::LocalCellId;
using hundun::mesh::LocalFaceId;
using hundun::mesh::LogicalFace;
using hundun::mesh::MeshTopology;
using hundun::mesh::PatchPairingKind;
using hundun::runtime::Box3;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kExtent{7, 5, 3};
constexpr std::array<bool, 3> kPeriodic{true, false, false};

static_assert(!std::is_default_constructible_v<BoundaryPatch>);
static_assert(!std::is_copy_constructible_v<BoundaryPatch>);

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool contains(Box3 box, Int3 cell) {
  return cell.x >= box.begin.x && cell.x < box.end.x && cell.y >= box.begin.y &&
         cell.y < box.end.y && cell.z >= box.begin.z && cell.z < box.end.z;
}

Int3 process_grid_for(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {2, 1, 1};
    case 4:
      return {2, 2, 1};
    default:
      throw hundun::runtime::Error("unsupported mesh-topology test size");
  }
}

Int3 extent_one_grid_for(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {1, 2, 1};
    case 4:
      return {1, 2, 2};
    default:
      throw hundun::runtime::Error("unsupported extent-one test size");
  }
}

std::uint64_t cell_count(Int3 extent) {
  return static_cast<std::uint64_t>(extent.x) *
         static_cast<std::uint64_t>(extent.y) *
         static_cast<std::uint64_t>(extent.z);
}

GlobalCellId cell_id(Int3 extent, Int3 cell) {
  return (static_cast<std::uint64_t>(cell.z) *
              static_cast<std::uint64_t>(extent.y) +
          static_cast<std::uint64_t>(cell.y)) *
             static_cast<std::uint64_t>(extent.x) +
         static_cast<std::uint64_t>(cell.x);
}

std::uint64_t face_count(Int3 extent, FaceAxis axis) {
  switch (axis) {
    case FaceAxis::x:
      return static_cast<std::uint64_t>(extent.x + 1) *
             static_cast<std::uint64_t>(extent.y) *
             static_cast<std::uint64_t>(extent.z);
    case FaceAxis::y:
      return static_cast<std::uint64_t>(extent.x) *
             static_cast<std::uint64_t>(extent.y + 1) *
             static_cast<std::uint64_t>(extent.z);
    case FaceAxis::z:
      return static_cast<std::uint64_t>(extent.x) *
             static_cast<std::uint64_t>(extent.y) *
             static_cast<std::uint64_t>(extent.z + 1);
  }
  throw hundun::runtime::Error("invalid oracle face axis");
}

GlobalFaceId face_id(Int3 extent, LogicalFace face) {
  const std::uint64_t x_count = face_count(extent, FaceAxis::x);
  const std::uint64_t y_count = face_count(extent, FaceAxis::y);
  const auto i = static_cast<std::uint64_t>(face.coordinate.x);
  const auto j = static_cast<std::uint64_t>(face.coordinate.y);
  const auto k = static_cast<std::uint64_t>(face.coordinate.z);
  switch (face.axis) {
    case FaceAxis::x:
      return (k * static_cast<std::uint64_t>(extent.y) + j) *
                 static_cast<std::uint64_t>(extent.x + 1) +
             i;
    case FaceAxis::y:
      return x_count +
             (k * static_cast<std::uint64_t>(extent.y + 1) + j) *
                 static_cast<std::uint64_t>(extent.x) +
             i;
    case FaceAxis::z:
      return x_count + y_count +
             (k * static_cast<std::uint64_t>(extent.y) + j) *
                 static_cast<std::uint64_t>(extent.x) +
             i;
  }
  throw hundun::runtime::Error("invalid oracle face axis");
}

std::array<LogicalFace, 6> cell_faces(Int3 cell) {
  return {LogicalFace{FaceAxis::x, Int3{cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::x, Int3{cell.x + 1, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, Int3{cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, Int3{cell.x, cell.y + 1, cell.z}},
          LogicalFace{FaceAxis::z, Int3{cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::z, Int3{cell.x, cell.y, cell.z + 1}}};
}

struct Connectivity {
  Int3 owner;
  std::optional<Int3> neighbour;
};

Connectivity connectivity(Int3 extent, std::array<bool, 3> periodic,
                          LogicalFace face) {
  Int3 lower = face.coordinate;
  Int3 upper = face.coordinate;
  int plane = 0;
  int cells = 0;
  bool wraps = false;
  switch (face.axis) {
    case FaceAxis::x:
      plane = face.coordinate.x;
      cells = extent.x;
      lower.x = plane - 1;
      upper.x = plane;
      wraps = periodic[0];
      break;
    case FaceAxis::y:
      plane = face.coordinate.y;
      cells = extent.y;
      lower.y = plane - 1;
      upper.y = plane;
      wraps = periodic[1];
      break;
    case FaceAxis::z:
      plane = face.coordinate.z;
      cells = extent.z;
      lower.z = plane - 1;
      upper.z = plane;
      wraps = periodic[2];
      break;
  }
  if (plane > 0 && plane < cells) {
    return {lower, upper};
  }

  Int3 owner = plane == 0 ? upper : lower;
  if (!wraps) {
    return {owner, std::nullopt};
  }
  Int3 opposite = owner;
  switch (face.axis) {
    case FaceAxis::x:
      opposite.x = plane == 0 ? extent.x - 1 : 0;
      break;
    case FaceAxis::y:
      opposite.y = plane == 0 ? extent.y - 1 : 0;
      break;
    case FaceAxis::z:
      opposite.z = plane == 0 ? extent.z - 1 : 0;
      break;
  }
  return {owner, opposite};
}

std::optional<std::uint32_t> expected_patch(Int3 extent, LogicalFace face) {
  switch (face.axis) {
    case FaceAxis::x:
      if (face.coordinate.x == 0)
        return 0U;
      if (face.coordinate.x == extent.x)
        return 1U;
      return std::nullopt;
    case FaceAxis::y:
      if (face.coordinate.y == 0)
        return 2U;
      if (face.coordinate.y == extent.y)
        return 3U;
      return std::nullopt;
    case FaceAxis::z:
      if (face.coordinate.z == 0)
        return 4U;
      if (face.coordinate.z == extent.z)
        return 5U;
      return std::nullopt;
  }
  throw hundun::runtime::Error("invalid oracle face axis");
}

LogicalFace opposite_periodic_face(Int3 extent, LogicalFace face) {
  switch (face.axis) {
    case FaceAxis::x:
      face.coordinate.x = face.coordinate.x == 0 ? extent.x : 0;
      break;
    case FaceAxis::y:
      face.coordinate.y = face.coordinate.y == 0 ? extent.y : 0;
      break;
    case FaceAxis::z:
      face.coordinate.z = face.coordinate.z == 0 ? extent.z : 0;
      break;
  }
  return face;
}

template <class Function> void expect_error(Function&& function) {
  bool rejected = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    rejected = !std::string(error.what()).empty();
  }
  HUNDUN_CHECK(rejected);
}

std::vector<std::uint64_t> gather_ids(const MpiContext& mpi,
                                      const std::vector<std::uint64_t>& local) {
  HUNDUN_CHECK(local.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> displacements(static_cast<std::size_t>(mpi.size()));
  int total = 0;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    displacements[static_cast<std::size_t>(rank)] = total;
    HUNDUN_CHECK(counts[static_cast<std::size_t>(rank)] >= 0);
    HUNDUN_CHECK(total <= std::numeric_limits<int>::max() -
                              counts[static_cast<std::size_t>(rank)]);
    total += counts[static_cast<std::size_t>(rank)];
  }
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(total));
  HUNDUN_CHECK(MPI_Allgatherv(local.data(), local_count, MPI_UINT64_T,
                              gathered.data(), counts.data(),
                              displacements.data(), MPI_UINT64_T,
                              mpi.comm()) == MPI_SUCCESS);
  return gathered;
}

std::vector<GlobalCellId> expected_ghost_cells(Box3 owned, Int3 extent,
                                               std::array<bool, 3> periodic) {
  std::vector<GlobalCellId> ghosts;
  constexpr std::array<Int3, 6> offsets{Int3{-1, 0, 0}, Int3{1, 0, 0},
                                        Int3{0, -1, 0}, Int3{0, 1, 0},
                                        Int3{0, 0, -1}, Int3{0, 0, 1}};
  for (int k = owned.begin.z; k < owned.end.z; ++k) {
    for (int j = owned.begin.y; j < owned.end.y; ++j) {
      for (int i = owned.begin.x; i < owned.end.x; ++i) {
        for (const Int3 offset : offsets) {
          Int3 adjacent{i + offset.x, j + offset.y, k + offset.z};
          bool present = true;
          int* coordinates[] = {&adjacent.x, &adjacent.y, &adjacent.z};
          const int cells[] = {extent.x, extent.y, extent.z};
          for (std::size_t axis = 0; axis < 3U; ++axis) {
            if (*coordinates[axis] < 0 || *coordinates[axis] >= cells[axis]) {
              if (!periodic[axis]) {
                present = false;
              } else {
                *coordinates[axis] =
                    *coordinates[axis] < 0 ? cells[axis] - 1 : 0;
              }
            }
          }
          if (present && !contains(owned, adjacent)) {
            ghosts.push_back(cell_id(extent, adjacent));
          }
        }
      }
    }
  }
  std::sort(ghosts.begin(), ghosts.end());
  ghosts.erase(std::unique(ghosts.begin(), ghosts.end()), ghosts.end());
  return ghosts;
}

std::vector<GlobalFaceId> expected_local_faces(Box3 owned, Int3 extent) {
  std::vector<GlobalFaceId> faces;
  for (int k = owned.begin.z; k < owned.end.z; ++k) {
    for (int j = owned.begin.y; j < owned.end.y; ++j) {
      for (int i = owned.begin.x; i < owned.end.x; ++i) {
        for (const LogicalFace face : cell_faces(Int3{i, j, k})) {
          faces.push_back(face_id(extent, face));
        }
      }
    }
  }
  std::sort(faces.begin(), faces.end());
  faces.erase(std::unique(faces.begin(), faces.end()), faces.end());
  return faces;
}

void test_counts_ids_and_local_layout(const MeshTopology& topology) {
  const Box3 owned = topology.owned_global_box();
  HUNDUN_CHECK(same(topology.global_extent(), kExtent));
  HUNDUN_CHECK(topology.global_cell_count() == cell_count(kExtent));
  std::uint64_t total_faces = 0;
  for (const FaceAxis axis : {FaceAxis::x, FaceAxis::y, FaceAxis::z}) {
    HUNDUN_CHECK(topology.global_face_count(axis) == face_count(kExtent, axis));
    total_faces += face_count(kExtent, axis);
  }
  HUNDUN_CHECK(topology.global_face_count() == total_faces);

  LocalCellId local_cell = 0;
  for (int k = owned.begin.z; k < owned.end.z; ++k) {
    for (int j = owned.begin.y; j < owned.end.y; ++j) {
      for (int i = owned.begin.x; i < owned.end.x; ++i) {
        const Int3 global{i, j, k};
        HUNDUN_CHECK(topology.cell_ownership(local_cell) ==
                     EntityOwnership::owned);
        HUNDUN_CHECK(same(topology.global_cell(local_cell), global));
        HUNDUN_CHECK(topology.global_cell_id(local_cell) ==
                     cell_id(kExtent, global));
        HUNDUN_CHECK(topology.global_cell_id(global) ==
                     cell_id(kExtent, global));
        HUNDUN_CHECK(topology.find_local_cell(cell_id(kExtent, global)) ==
                     local_cell);
        ++local_cell;
      }
    }
  }
  HUNDUN_CHECK(local_cell == topology.owned_cell_count());

  const std::vector<GlobalCellId> ghosts =
      expected_ghost_cells(owned, kExtent, kPeriodic);
  HUNDUN_CHECK(ghosts.size() == topology.ghost_cell_count());
  for (std::size_t index = 0; index < ghosts.size(); ++index) {
    const LocalCellId local = topology.owned_cell_count() + index;
    HUNDUN_CHECK(topology.cell_ownership(local) == EntityOwnership::ghost);
    HUNDUN_CHECK(topology.global_cell_id(local) == ghosts[index]);
    HUNDUN_CHECK(topology.find_local_cell(ghosts[index]) == local);
  }
  HUNDUN_CHECK(topology.local_cell_count() ==
               topology.owned_cell_count() + topology.ghost_cell_count());

  const std::vector<GlobalFaceId> expected_faces =
      expected_local_faces(owned, kExtent);
  HUNDUN_CHECK(expected_faces.size() == topology.local_face_count());
  std::size_t owned_faces = 0;
  for (LocalFaceId local = 0; local < topology.local_face_count(); ++local) {
    const GlobalFaceId expected_id = expected_faces[local];
    const LogicalFace face = topology.logical_face(local);
    HUNDUN_CHECK(topology.global_face_id(face) == expected_id);
    HUNDUN_CHECK(topology.global_face_id(local) == expected_id);
    HUNDUN_CHECK(topology.find_local_face(expected_id) == local);

    const Connectivity expected = connectivity(kExtent, kPeriodic, face);
    const LocalCellId owner = topology.owner(local);
    HUNDUN_CHECK(topology.global_cell_id(owner) ==
                 cell_id(kExtent, expected.owner));
    HUNDUN_CHECK(topology.cell_ownership(owner) ==
                 (contains(owned, expected.owner) ? EntityOwnership::owned
                                                  : EntityOwnership::ghost));
    const bool should_be_owned = contains(owned, expected.owner);
    HUNDUN_CHECK(
        topology.face_ownership(local) ==
        (should_be_owned ? EntityOwnership::owned : EntityOwnership::ghost));
    owned_faces += should_be_owned ? 1U : 0U;

    const auto neighbour = topology.neighbour(local);
    HUNDUN_CHECK(neighbour.has_value() == expected.neighbour.has_value());
    if (expected.neighbour.has_value()) {
      HUNDUN_CHECK(topology.global_cell_id(*neighbour) ==
                   cell_id(kExtent, *expected.neighbour));
    }
    HUNDUN_CHECK(topology.patch_id(local) == expected_patch(kExtent, face));
  }
  HUNDUN_CHECK(owned_faces == topology.owned_face_count());
  HUNDUN_CHECK(topology.local_face_count() ==
               topology.owned_face_count() + topology.ghost_face_count());
}

void test_patches_and_periodic_pairs(const MeshTopology& topology) {
  constexpr std::array<std::string_view, 6> names{"x_min", "x_max", "y_min",
                                                  "y_max", "z_min", "z_max"};
  constexpr std::array<std::uint32_t, 6> pairs{1U, 0U, 3U, 2U, 5U, 4U};
  HUNDUN_CHECK(topology.patches().size() == 6U);
  for (std::uint32_t id = 0; id < 6U; ++id) {
    const BoundaryPatch& patch = topology.patch(id);
    HUNDUN_CHECK(&patch == &topology.patches()[id]);
    HUNDUN_CHECK(patch.stable_id() == id);
    HUNDUN_CHECK(patch.name() == names[id]);
    const bool periodic = id < 2U;
    HUNDUN_CHECK(patch.pairing_kind() == (periodic ? PatchPairingKind::periodic
                                                   : PatchPairingKind::none));
    HUNDUN_CHECK(
        patch.paired_patch_id() ==
        (periodic ? std::optional<std::uint32_t>{pairs[id]} : std::nullopt));
    HUNDUN_CHECK(std::is_sorted(
        patch.local_faces().begin(), patch.local_faces().end(),
        [&](LocalFaceId lhs, LocalFaceId rhs) {
          return topology.global_face_id(lhs) < topology.global_face_id(rhs);
        }));
    for (const LocalFaceId face : patch.local_faces()) {
      HUNDUN_CHECK(patch.contains(face));
      HUNDUN_CHECK(topology.patch_id(face) == id);
      HUNDUN_CHECK(topology.face_ownership(face) == EntityOwnership::owned);
      const auto pair = topology.periodic_pair(face);
      HUNDUN_CHECK(pair.has_value() == periodic);
      if (pair.has_value()) {
        const LogicalFace logical = topology.logical_face(face);
        const GlobalFaceId expected =
            face_id(kExtent, opposite_periodic_face(kExtent, logical));
        HUNDUN_CHECK(*pair == expected);
        const LogicalFace opposite = opposite_periodic_face(kExtent, logical);
        HUNDUN_CHECK(
            face_id(kExtent, opposite_periodic_face(kExtent, opposite)) ==
            topology.global_face_id(face));
      }
    }
  }
}

void test_invalid_queries(const MeshTopology& topology) {
  expect_error([&] {
    static_cast<void>(topology.global_cell_id(Int3{-1, 0, 0}));
  });
  expect_error([&] {
    static_cast<void>(topology.global_cell_id(Int3{kExtent.x, 0, 0}));
  });
  expect_error([&] {
    static_cast<void>(topology.global_cell(topology.local_cell_count()));
  });
  expect_error([&] {
    static_cast<void>(topology.cell_ownership(topology.local_cell_count()));
  });
  expect_error([&] {
    static_cast<void>(topology.find_local_cell(topology.global_cell_count()));
  });
  expect_error([&] {
    static_cast<void>(topology.global_face_count(static_cast<FaceAxis>(255)));
  });
  expect_error([&] {
    static_cast<void>(
        topology.global_face_id(LogicalFace{FaceAxis::x, Int3{-1, 0, 0}}));
  });
  expect_error([&] {
    static_cast<void>(topology.global_face_id(
        LogicalFace{FaceAxis::x, Int3{kExtent.x + 1, 0, 0}}));
  });
  expect_error([&] {
    static_cast<void>(topology.global_face_id(
        LogicalFace{FaceAxis::y, Int3{0, kExtent.y + 1, 0}}));
  });
  expect_error([&] {
    static_cast<void>(topology.global_face_id(
        LogicalFace{FaceAxis::z, Int3{0, 0, kExtent.z + 1}}));
  });
  expect_error([&] {
    static_cast<void>(topology.global_face_id(
        LogicalFace{static_cast<FaceAxis>(255), Int3{0, 0, 0}}));
  });
  expect_error([&] {
    static_cast<void>(topology.logical_face(topology.local_face_count()));
  });
  expect_error(
      [&] { static_cast<void>(topology.owner(topology.local_face_count())); });
  expect_error([&] {
    static_cast<void>(topology.find_local_face(topology.global_face_count()));
  });
  expect_error([&] { static_cast<void>(topology.patch(6U)); });

  if (topology.local_cell_count() < topology.global_cell_count()) {
    GlobalCellId absent = 0;
    while (topology.find_local_cell(absent).has_value()) {
      ++absent;
    }
    HUNDUN_CHECK(absent < topology.global_cell_count());
    HUNDUN_CHECK(!topology.find_local_cell(absent).has_value());
  }
  if (topology.local_face_count() < topology.global_face_count()) {
    GlobalFaceId absent = 0;
    while (topology.find_local_face(absent).has_value()) {
      ++absent;
    }
    HUNDUN_CHECK(absent < topology.global_face_count());
    HUNDUN_CHECK(!topology.find_local_face(absent).has_value());
  }
}

void test_global_inventories(const MpiContext& mpi,
                             const MeshTopology& topology) {
  std::vector<std::uint64_t> owned_cells;
  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    owned_cells.push_back(topology.global_cell_id(local));
  }
  owned_cells = gather_ids(mpi, owned_cells);
  std::sort(owned_cells.begin(), owned_cells.end());
  HUNDUN_CHECK(owned_cells.size() ==
               static_cast<std::size_t>(topology.global_cell_count()));
  for (std::uint64_t id = 0; id < topology.global_cell_count(); ++id) {
    HUNDUN_CHECK(owned_cells[static_cast<std::size_t>(id)] == id);
  }

  std::vector<std::uint64_t> owned_faces;
  for (LocalFaceId local = 0; local < topology.local_face_count(); ++local) {
    if (topology.face_ownership(local) == EntityOwnership::owned) {
      owned_faces.push_back(topology.global_face_id(local));
    }
  }
  owned_faces = gather_ids(mpi, owned_faces);
  std::sort(owned_faces.begin(), owned_faces.end());
  HUNDUN_CHECK(owned_faces.size() ==
               static_cast<std::size_t>(topology.global_face_count()));
  for (std::uint64_t id = 0; id < topology.global_face_count(); ++id) {
    HUNDUN_CHECK(owned_faces[static_cast<std::size_t>(id)] == id);
  }

  for (std::uint32_t patch_id = 0; patch_id < 6U; ++patch_id) {
    std::vector<std::uint64_t> local;
    for (const LocalFaceId face : topology.patch(patch_id).local_faces()) {
      local.push_back(topology.global_face_id(face));
    }
    std::vector<std::uint64_t> gathered = gather_ids(mpi, local);
    std::sort(gathered.begin(), gathered.end());
    std::vector<std::uint64_t> expected;
    const bool x_axis = patch_id < 2U;
    const bool y_axis = patch_id >= 2U && patch_id < 4U;
    const int outer = x_axis ? kExtent.z : (y_axis ? kExtent.z : kExtent.y);
    const int inner = x_axis ? kExtent.y : kExtent.x;
    for (int a = 0; a < outer; ++a) {
      for (int b = 0; b < inner; ++b) {
        LogicalFace face{};
        if (x_axis) {
          face = {FaceAxis::x, Int3{patch_id == 0U ? 0 : kExtent.x, b, a}};
        } else if (y_axis) {
          face = {FaceAxis::y, Int3{b, patch_id == 2U ? 0 : kExtent.y, a}};
        } else {
          face = {FaceAxis::z, Int3{b, a, patch_id == 4U ? 0 : kExtent.z}};
        }
        expected.push_back(face_id(kExtent, face));
      }
    }
    std::sort(expected.begin(), expected.end());
    HUNDUN_CHECK(gathered == expected);
  }
}

void test_partition_face(const MpiContext& mpi, const MeshTopology& topology) {
  if (mpi.size() == 1) {
    return;
  }
  const int quotient = kExtent.x / 2;
  const int remainder = kExtent.x % 2;
  const int plane = quotient + (remainder > 0 ? 1 : 0);
  const LogicalFace logical{FaceAxis::x, Int3{plane, 0, 0}};
  const GlobalFaceId target = face_id(kExtent, logical);
  constexpr std::uint64_t absent = std::numeric_limits<std::uint64_t>::max();
  std::array<std::uint64_t, 5> local{absent, absent, absent, absent, absent};
  const auto found = topology.find_local_face(target);
  if (found.has_value()) {
    local[0] =
        topology.face_ownership(*found) == EntityOwnership::owned ? 1U : 0U;
    local[1] = topology.global_cell_id(topology.owner(*found));
    local[2] = topology.global_cell_id(*topology.neighbour(*found));
    local[3] = topology.patch_id(*found).has_value()
                   ? *topology.patch_id(*found)
                   : absent;
    local[4] = topology.global_face_id(*found);
  }
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(mpi.size()) *
                                      local.size());
  HUNDUN_CHECK(MPI_Allgather(local.data(), static_cast<int>(local.size()),
                             MPI_UINT64_T, gathered.data(),
                             static_cast<int>(local.size()), MPI_UINT64_T,
                             mpi.comm()) == MPI_SUCCESS);
  int copies = 0;
  int owned_copies = 0;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const std::size_t base = static_cast<std::size_t>(rank) * local.size();
    if (gathered[base] == absent) {
      continue;
    }
    ++copies;
    owned_copies += gathered[base] == 1U ? 1 : 0;
    HUNDUN_CHECK(gathered[base + 1U] ==
                 cell_id(kExtent, Int3{plane - 1, 0, 0}));
    HUNDUN_CHECK(gathered[base + 2U] == cell_id(kExtent, Int3{plane, 0, 0}));
    HUNDUN_CHECK(gathered[base + 3U] == absent);
    HUNDUN_CHECK(gathered[base + 4U] == target);
  }
  HUNDUN_CHECK(copies == 2);
  HUNDUN_CHECK(owned_copies == 1);
}

void test_extent_one_periodic(const MpiContext& mpi) {
  constexpr Int3 extent{1, 5, 3};
  constexpr std::array<bool, 3> periodic{true, false, false};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, periodic,
      DecompositionOptions{extent_one_grid_for(mpi.size())});
  const MeshTopology topology(decomposition);
  for (const std::uint32_t patch_id : {0U, 1U}) {
    for (const LocalFaceId face : topology.patch(patch_id).local_faces()) {
      const LocalCellId owner = topology.owner(face);
      const auto neighbour = topology.neighbour(face);
      HUNDUN_CHECK(neighbour.has_value());
      HUNDUN_CHECK(topology.global_cell_id(owner) ==
                   topology.global_cell_id(*neighbour));
      HUNDUN_CHECK(owner == *neighbour);
      const auto pair = topology.periodic_pair(face);
      HUNDUN_CHECK(pair.has_value());
      HUNDUN_CHECK(*pair != topology.global_face_id(face));
    }
  }
  mpi.barrier();
}

void test_overflow_rejection(const MpiContext& mpi) {
  constexpr Int3 near_limit{INT_MAX, INT_MAX, 2};
  auto decomposition = StructuredDecomposition::create(
      mpi, near_limit, kPeriodic,
      DecompositionOptions{process_grid_for(mpi.size())});
  expect_error([&] { static_cast<void>(MeshTopology(decomposition)); });
  mpi.barrier();
}

void run_tests(const MpiContext& mpi, const MeshTopology& topology) {
  test_counts_ids_and_local_layout(topology);
  test_patches_and_periodic_pairs(topology);
  test_invalid_queries(topology);
  test_global_inventories(mpi, topology);
  test_partition_face(mpi, topology);
  test_extent_one_periodic(mpi);
  test_overflow_rejection(mpi);
}

}  // namespace

int main(int argc, char** argv) {
  std::optional<MeshTopology> detached;
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    result = hundun::test::run([&] {
      {
        auto decomposition = StructuredDecomposition::create(
            mpi, kExtent, kPeriodic,
            DecompositionOptions{process_grid_for(mpi.size())});
        detached.emplace(decomposition);
      }
      HUNDUN_CHECK(detached->global_cell_count() == cell_count(kExtent));
      run_tests(mpi, *detached);
      mpi.barrier();
    });
  }

  const int post_finalize = hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);
    HUNDUN_CHECK(detached->global_cell_count() == cell_count(kExtent));
    HUNDUN_CHECK(detached->patch(0U).name() == "x_min");
    detached.reset();
  });
  return result == EXIT_SUCCESS ? post_finalize : result;
}
