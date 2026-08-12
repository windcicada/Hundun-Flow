// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "hundun/rt_types.hpp"
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
#include <utility>
#include <vector>

namespace {

using hundun::runtime::Box3;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kGlobalExtent{17, 11, 7};
constexpr std::array<bool, 3> kPeriodic{true, false, true};

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

std::array<int, 3> as_array(Int3 value) {
  return {value.x, value.y, value.z};
}

std::array<int, 2> expected_axis_bounds(int global_cells,
                                        int process_count,
                                        int process_coordinate) {
  const int quotient = global_cells / process_count;
  const int remainder = global_cells % process_count;
  const int earlier_remainder_cells =
      process_coordinate < remainder ? process_coordinate : remainder;
  const int begin = process_coordinate * quotient + earlier_remainder_cells;
  const int end = begin + quotient +
                  (process_coordinate < remainder ? 1 : 0);
  return {begin, end};
}

template <class Function>
void expect_runtime_error(Function&& function) {
  bool threw = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()).empty() == false);
  }
  HUNDUN_CHECK(threw);
}

template <class Function>
void expect_runtime_error_with_text(Function&& function,
                                    std::string_view expected_text) {
  bool threw = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string_view{error.what()}.find(expected_text) !=
                 std::string_view::npos);
  }
  HUNDUN_CHECK(threw);
}

template <class Function>
void expect_collective_error(const MpiContext& mpi, Function&& function) {
  bool threw = false;
  std::string message;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    message = error.what();
  }

  const int local_threw = threw ? 1 : 0;
  int total_threw = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_threw, &total_threw, 1, MPI_INT, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);

  HUNDUN_CHECK(message.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  const int local_length = static_cast<int>(message.size());
  std::vector<int> lengths(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_length, 1, MPI_INT, lengths.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> displacements(static_cast<std::size_t>(mpi.size()));
  int total_length = 0;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    displacements[static_cast<std::size_t>(rank)] = total_length;
    HUNDUN_CHECK(lengths[static_cast<std::size_t>(rank)] >= 0);
    HUNDUN_CHECK(total_length <=
                 std::numeric_limits<int>::max() -
                     lengths[static_cast<std::size_t>(rank)]);
    total_length += lengths[static_cast<std::size_t>(rank)];
  }
  std::vector<char> messages(static_cast<std::size_t>(total_length));
  HUNDUN_CHECK(MPI_Allgatherv(message.data(), local_length, MPI_BYTE,
                              messages.data(), lengths.data(),
                              displacements.data(), MPI_BYTE,
                              mpi.comm()) == MPI_SUCCESS);

  HUNDUN_CHECK(total_threw == mpi.size());
  HUNDUN_CHECK(!message.empty());
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const std::size_t begin = static_cast<std::size_t>(
        displacements[static_cast<std::size_t>(rank)]);
    const std::size_t length =
        static_cast<std::size_t>(lengths[static_cast<std::size_t>(rank)]);
    HUNDUN_CHECK(std::string(messages.data() + begin, length) == message);
  }

  mpi.barrier();
}

bool boxes_overlap(const Box3& lhs, const Box3& rhs) {
  return std::max(lhs.begin.x, rhs.begin.x) <
             std::min(lhs.end.x, rhs.end.x) &&
         std::max(lhs.begin.y, rhs.begin.y) <
             std::min(lhs.end.y, rhs.end.y) &&
         std::max(lhs.begin.z, rhs.begin.z) <
             std::min(lhs.end.z, rhs.end.z);
}

Int3 expected_invariant_grid(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {1, 2, 1};
    case 4:
      return {2, 2, 1};
    case 12:
      return {4, 3, 1};
    default:
      throw hundun::runtime::Error("unsupported decomposition test size");
  }
}

void test_cartesian_topology(const StructuredDecomposition& decomposition,
                             const MpiContext& mpi,
                             Int3 expected_grid) {
  HUNDUN_CHECK(same(decomposition.global_extent(), kGlobalExtent));

  int topology = MPI_UNDEFINED;
  HUNDUN_CHECK(MPI_Topo_test(decomposition.comm(), &topology) == MPI_SUCCESS);
  HUNDUN_CHECK(topology == MPI_CART);

  int comparison = MPI_UNEQUAL;
  HUNDUN_CHECK(MPI_Comm_compare(mpi.comm(), decomposition.comm(),
                                &comparison) == MPI_SUCCESS);
  HUNDUN_CHECK(comparison == MPI_CONGRUENT);

  int cart_rank = -1;
  HUNDUN_CHECK(MPI_Comm_rank(decomposition.comm(), &cart_rank) == MPI_SUCCESS);
  HUNDUN_CHECK(cart_rank == mpi.rank());

  std::array<int, 3> dims{};
  std::array<int, 3> periods{};
  std::array<int, 3> coordinates{};
  HUNDUN_CHECK(MPI_Cart_get(decomposition.comm(), 3, dims.data(),
                            periods.data(), coordinates.data()) == MPI_SUCCESS);
  HUNDUN_CHECK(dims == as_array(expected_grid));
  HUNDUN_CHECK(dims == as_array(decomposition.process_grid()));
  HUNDUN_CHECK(coordinates ==
               as_array(decomposition.process_coordinates()));
  HUNDUN_CHECK(periods[0] != 0);
  HUNDUN_CHECK(periods[1] == 0);
  HUNDUN_CHECK(periods[2] != 0);
}

std::vector<Box3> gather_boxes(const StructuredDecomposition& decomposition,
                               const MpiContext& mpi) {
  const Box3 owned = decomposition.owned_box();
  const std::array<int, 6> local{owned.begin.x, owned.begin.y, owned.begin.z,
                                 owned.end.x, owned.end.y, owned.end.z};
  std::vector<int> packed(static_cast<std::size_t>(mpi.size()) * local.size());
  HUNDUN_CHECK(MPI_Allgather(local.data(), static_cast<int>(local.size()),
                             MPI_INT, packed.data(),
                             static_cast<int>(local.size()), MPI_INT,
                             decomposition.comm()) == MPI_SUCCESS);

  std::vector<Box3> boxes;
  boxes.reserve(static_cast<std::size_t>(mpi.size()));
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const std::size_t base = static_cast<std::size_t>(rank) * local.size();
    boxes.push_back(Box3{Int3{packed[base], packed[base + 1U],
                                  packed[base + 2U]},
                         Int3{packed[base + 3U], packed[base + 4U],
                                  packed[base + 5U]}});
  }
  return boxes;
}

void test_boxes_and_coverage(const StructuredDecomposition& decomposition,
                             const MpiContext& mpi) {
  const Int3 global_extent = decomposition.global_extent();
  const Int3 process_grid = decomposition.process_grid();
  const Int3 process_coordinates = decomposition.process_coordinates();
  const Box3 owned = decomposition.owned_box();
  const auto expected_x = expected_axis_bounds(
      global_extent.x, process_grid.x, process_coordinates.x);
  const auto expected_y = expected_axis_bounds(
      global_extent.y, process_grid.y, process_coordinates.y);
  const auto expected_z = expected_axis_bounds(
      global_extent.z, process_grid.z, process_coordinates.z);
  HUNDUN_CHECK(owned.begin.x == expected_x[0]);
  HUNDUN_CHECK(owned.end.x == expected_x[1]);
  HUNDUN_CHECK(owned.begin.y == expected_y[0]);
  HUNDUN_CHECK(owned.end.y == expected_y[1]);
  HUNDUN_CHECK(owned.begin.z == expected_z[0]);
  HUNDUN_CHECK(owned.end.z == expected_z[1]);

  const Int3 local_extent = decomposition.local_extent();
  HUNDUN_CHECK(local_extent.x > 0);
  HUNDUN_CHECK(local_extent.y > 0);
  HUNDUN_CHECK(local_extent.z > 0);
  HUNDUN_CHECK(hundun::runtime::volume(local_extent) > 0);

  const std::vector<Box3> boxes = gather_boxes(decomposition, mpi);
  std::int64_t summed_volume = 0;
  for (std::size_t index = 0; index < boxes.size(); ++index) {
    const Box3& box = boxes[index];
    HUNDUN_CHECK(box.begin.x >= 0 && box.end.x <= kGlobalExtent.x);
    HUNDUN_CHECK(box.begin.y >= 0 && box.end.y <= kGlobalExtent.y);
    HUNDUN_CHECK(box.begin.z >= 0 && box.end.z <= kGlobalExtent.z);
    const Int3 extent{box.end.x - box.begin.x, box.end.y - box.begin.y,
                      box.end.z - box.begin.z};
    HUNDUN_CHECK(extent.x > 0 && extent.y > 0 && extent.z > 0);
    summed_volume += hundun::runtime::volume(extent);
    for (std::size_t other = index + 1U; other < boxes.size(); ++other) {
      HUNDUN_CHECK(boxes_overlap(box, boxes[other]) == false);
    }
  }

  constexpr std::int64_t total =
      static_cast<std::int64_t>(kGlobalExtent.x) * kGlobalExtent.y *
      kGlobalExtent.z;
  HUNDUN_CHECK(summed_volume == total);

  std::vector<int> coverage(static_cast<std::size_t>(total), 0);
  for (const Box3& box : boxes) {
    for (int k = box.begin.z; k < box.end.z; ++k) {
      for (int j = box.begin.y; j < box.end.y; ++j) {
        for (int i = box.begin.x; i < box.end.x; ++i) {
          const std::int64_t id =
              (static_cast<std::int64_t>(k) * kGlobalExtent.y + j) *
                  kGlobalExtent.x +
              i;
          ++coverage[static_cast<std::size_t>(id)];
        }
      }
    }
  }
  HUNDUN_CHECK(std::all_of(coverage.begin(), coverage.end(),
                           [](int visits) { return visits == 1; }));
}

void test_global_cells_and_ids(const StructuredDecomposition& decomposition,
                               const MpiContext& mpi) {
  const Int3 extent = decomposition.local_extent();
  const Box3 owned = decomposition.owned_box();
  const int local_count =
      static_cast<int>(hundun::runtime::volume(extent));
  std::vector<std::uint64_t> local_ids;
  local_ids.reserve(static_cast<std::size_t>(local_count));

  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        const Int3 local{i, j, k};
        const Int3 global = decomposition.global_cell(local);
        HUNDUN_CHECK(same(global,
                          Int3{owned.begin.x + i, owned.begin.y + j,
                               owned.begin.z + k}));
        const std::uint64_t expected =
            (static_cast<std::uint64_t>(global.z) *
                 static_cast<std::uint64_t>(kGlobalExtent.y) +
             static_cast<std::uint64_t>(global.y)) *
                static_cast<std::uint64_t>(kGlobalExtent.x) +
            static_cast<std::uint64_t>(global.x);
        const std::uint64_t id = decomposition.global_cell_id(local);
        HUNDUN_CHECK(id == expected);
        local_ids.push_back(id);
      }
    }
  }

  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1,
                             MPI_INT, decomposition.comm()) == MPI_SUCCESS);
  std::vector<int> displacements(static_cast<std::size_t>(mpi.size()));
  int gathered_count = 0;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    displacements[static_cast<std::size_t>(rank)] = gathered_count;
    gathered_count += counts[static_cast<std::size_t>(rank)];
  }
  std::vector<std::uint64_t> gathered_ids(
      static_cast<std::size_t>(gathered_count));
  HUNDUN_CHECK(MPI_Allgatherv(
                   local_ids.data(), local_count, MPI_UINT64_T,
                   gathered_ids.data(), counts.data(), displacements.data(),
                   MPI_UINT64_T, decomposition.comm()) == MPI_SUCCESS);
  std::sort(gathered_ids.begin(), gathered_ids.end());

  constexpr std::uint64_t total =
      static_cast<std::uint64_t>(kGlobalExtent.x) * kGlobalExtent.y *
      kGlobalExtent.z;
  HUNDUN_CHECK(gathered_ids.size() == static_cast<std::size_t>(total));
  for (std::uint64_t id = 0; id < total; ++id) {
    HUNDUN_CHECK(gathered_ids[static_cast<std::size_t>(id)] == id);
  }

  for (Int3 invalid : {Int3{-1, 0, 0}, Int3{0, -1, 0}, Int3{0, 0, -1},
                       Int3{extent.x, 0, 0}, Int3{0, extent.y, 0},
                       Int3{0, 0, extent.z}}) {
    expect_runtime_error(
        [&] { static_cast<void>(decomposition.global_cell(invalid)); });
    expect_runtime_error(
        [&] { static_cast<void>(decomposition.global_cell_id(invalid)); });
  }
}

int rank_with_coordinates(const std::vector<int>& packed_coordinates,
                          Int3 expected) {
  for (std::size_t rank = 0; rank * 3U < packed_coordinates.size(); ++rank) {
    const std::size_t base = rank * 3U;
    if (packed_coordinates[base] == expected.x &&
        packed_coordinates[base + 1U] == expected.y &&
        packed_coordinates[base + 2U] == expected.z) {
      return static_cast<int>(rank);
    }
  }
  return MPI_PROC_NULL;
}

void test_neighbors(const StructuredDecomposition& decomposition,
                    const MpiContext& mpi) {
  const Int3 coordinates = decomposition.process_coordinates();
  const Int3 grid = decomposition.process_grid();
  const std::array<int, 3> local_coordinates = as_array(coordinates);
  std::vector<int> packed_coordinates(
      static_cast<std::size_t>(mpi.size()) * local_coordinates.size());
  HUNDUN_CHECK(MPI_Allgather(
                   local_coordinates.data(),
                   static_cast<int>(local_coordinates.size()), MPI_INT,
                   packed_coordinates.data(),
                   static_cast<int>(local_coordinates.size()), MPI_INT,
                   decomposition.comm()) == MPI_SUCCESS);

  int visited = 0;
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        if (x == 0 && y == 0 && z == 0) {
          continue;
        }
        ++visited;
        Int3 target{coordinates.x + x, coordinates.y + y,
                    coordinates.z + z};
        bool physical_boundary = false;
        if (target.x < 0) {
          target.x = grid.x - 1;
        } else if (target.x >= grid.x) {
          target.x = 0;
        }
        if (target.y < 0 || target.y >= grid.y) {
          physical_boundary = true;
        }
        if (target.z < 0) {
          target.z = grid.z - 1;
        } else if (target.z >= grid.z) {
          target.z = 0;
        }

        const int expected =
            physical_boundary
                ? MPI_PROC_NULL
                : rank_with_coordinates(packed_coordinates, target);
        HUNDUN_CHECK(decomposition.neighbor_rank(Int3{x, y, z}) == expected);
      }
    }
  }
  HUNDUN_CHECK(visited == 26);

  for (Int3 invalid : {Int3{0, 0, 0}, Int3{-2, 0, 0}, Int3{2, 0, 0},
                       Int3{0, -2, 0}, Int3{0, 2, 0}, Int3{0, 0, -2},
                       Int3{0, 0, 2}}) {
    expect_runtime_error(
        [&] { static_cast<void>(decomposition.neighbor_rank(invalid)); });
  }
}

Int3 legal_explicit_grid(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {2, 1, 1};
    case 4:
      return {4, 1, 1};
    case 12:
      return {3, 4, 1};
    default:
      throw hundun::runtime::Error("unsupported explicit-grid test size");
  }
}

Int3 alternate_explicit_grid(int ranks) {
  switch (ranks) {
    case 2:
      return {1, 2, 1};
    case 4:
      return {2, 2, 1};
    case 12:
      return {4, 3, 1};
    default:
      throw hundun::runtime::Error("unsupported alternate-grid test size");
  }
}

void test_frozen_shape_selection(const MpiContext& mpi) {
  if (mpi.size() == 4) {
    auto decomposition = StructuredDecomposition::create(
        mpi, Int3{1, 2, 2}, {false, false, false});
    HUNDUN_CHECK(same(decomposition.process_grid(), Int3{1, 2, 2}));
  }
  if (mpi.size() == 12) {
    auto decomposition = StructuredDecomposition::create(
        mpi, Int3{1, 3, 4}, {false, false, false});
    HUNDUN_CHECK(same(decomposition.process_grid(), Int3{1, 3, 4}));
  }
  mpi.barrier();
}

void test_legal_explicit_grid(const MpiContext& mpi) {
  const Int3 expected = legal_explicit_grid(mpi.size());
  auto decomposition = StructuredDecomposition::create(
      mpi, kGlobalExtent, kPeriodic, DecompositionOptions{expected});
  HUNDUN_CHECK(same(decomposition.process_grid(), expected));
  mpi.barrier();
}

void test_collective_input_agreement(const MpiContext& mpi) {
  if (mpi.size() == 1) {
    return;
  }

  expect_collective_error(mpi, [&] {
    const Int3 extent = mpi.rank() == 0 ? Int3{18, 11, 7} : kGlobalExtent;
    static_cast<void>(
        StructuredDecomposition::create(mpi, extent, kPeriodic));
  });

  expect_collective_error(mpi, [&] {
    auto periodic = kPeriodic;
    periodic[1] = mpi.rank() == 0;
    static_cast<void>(
        StructuredDecomposition::create(mpi, kGlobalExtent, periodic));
  });

  expect_collective_error(mpi, [&] {
    DecompositionOptions options{};
    if (mpi.rank() == 0) {
      options.process_grid = legal_explicit_grid(mpi.size());
    }
    static_cast<void>(StructuredDecomposition::create(
        mpi, kGlobalExtent, kPeriodic, options));
  });

  expect_collective_error(mpi, [&] {
    const Int3 grid = mpi.rank() == 0
                          ? legal_explicit_grid(mpi.size())
                          : alternate_explicit_grid(mpi.size());
    static_cast<void>(StructuredDecomposition::create(
        mpi, kGlobalExtent, kPeriodic, DecompositionOptions{grid}));
  });
}

void test_collective_selection_failures(const MpiContext& mpi) {
  for (const Int3 grid : {Int3{0, 1, 1}, Int3{-1, 1, 1},
                          Int3{INT_MAX, INT_MAX, INT_MAX}}) {
    expect_collective_error(mpi, [&] {
      static_cast<void>(StructuredDecomposition::create(
          mpi, kGlobalExtent, kPeriodic, DecompositionOptions{grid}));
    });
  }

  const Int3 wrong_product =
      mpi.size() == 1 ? Int3{2, 1, 1} : Int3{1, 1, 1};
  expect_collective_error(mpi, [&] {
    static_cast<void>(StructuredDecomposition::create(
        mpi, kGlobalExtent, kPeriodic,
        DecompositionOptions{wrong_product}));
  });

  if (mpi.size() > 1) {
    const Int3 cells{mpi.size() - 1, 2, 2};
    const Int3 axis_too_large{mpi.size(), 1, 1};
    expect_collective_error(mpi, [&] {
      static_cast<void>(StructuredDecomposition::create(
          mpi, cells, kPeriodic, DecompositionOptions{axis_too_large}));
    });
  }

  if (mpi.size() == 12) {
    MPI_Comm subgroup = MPI_COMM_NULL;
    const int color = mpi.rank() < 6 ? 0 : 1;
    HUNDUN_CHECK(MPI_Comm_split(mpi.comm(), color, mpi.rank(), &subgroup) ==
                 MPI_SUCCESS);
    {
      MpiContext six_rank_context = MpiContext::duplicate(subgroup);
      HUNDUN_CHECK(MPI_Comm_free(&subgroup) == MPI_SUCCESS);
      expect_collective_error(six_rank_context, [&] {
        static_cast<void>(StructuredDecomposition::create(
            six_rank_context, Int3{2, 2, 2}, {false, false, false}));
      });
    }
    mpi.barrier();
  }
}

void test_move_ownership(const MpiContext& mpi) {
  std::optional<StructuredDecomposition> moved;
  MPI_Comm cartesian = MPI_COMM_NULL;
  {
    auto source = StructuredDecomposition::create(
        mpi, kGlobalExtent, kPeriodic);
    cartesian = source.comm();
    moved.emplace(std::move(source));
  }

  HUNDUN_CHECK(moved.has_value());
  HUNDUN_CHECK(moved->comm() == cartesian);
  int moved_rank = -1;
  HUNDUN_CHECK(MPI_Comm_rank(moved->comm(), &moved_rank) == MPI_SUCCESS);
  HUNDUN_CHECK(moved_rank == mpi.rank());
}

void test_extent_validation(const MpiContext& mpi) {
  for (Int3 invalid : {Int3{0, 11, 7}, Int3{17, 0, 7}, Int3{17, 11, 0},
                       Int3{-1, 11, 7}, Int3{17, -1, 7},
                       Int3{17, 11, -1}}) {
    expect_runtime_error([&] {
      static_cast<void>(StructuredDecomposition::create(
          mpi, invalid, kPeriodic));
    });
  }

  expect_runtime_error([&] {
    static_cast<void>(StructuredDecomposition::create(
        mpi, Int3{INT_MAX, INT_MAX, INT_MAX}, kPeriodic));
  });

  const Int3 near_limit{INT_MAX, INT_MAX, 2};
  auto accepted = StructuredDecomposition::create(
      mpi, near_limit, kPeriodic);
  const Int3 last_local{accepted.local_extent().x - 1,
                        accepted.local_extent().y - 1,
                        accepted.local_extent().z - 1};
  const Int3 global = accepted.global_cell(last_local);
  const std::uint64_t expected =
      (static_cast<std::uint64_t>(global.z) *
           static_cast<std::uint64_t>(near_limit.y) +
       static_cast<std::uint64_t>(global.y)) *
          static_cast<std::uint64_t>(near_limit.x) +
      static_cast<std::uint64_t>(global.x);
  HUNDUN_CHECK(accepted.global_cell_id(last_local) == expected);
  HUNDUN_CHECK(expected <=
               static_cast<std::uint64_t>(
                   std::numeric_limits<std::int64_t>::max()));

  if (mpi.size() > 1) {
    expect_runtime_error([&] {
      static_cast<void>(StructuredDecomposition::create(
          mpi, Int3{1, 1, 1}, kPeriodic));
    });
  }
}

void test_intercommunicator_rejection(const MpiContext& mpi) {
  if (mpi.size() != 2) {
    return;
  }

  MPI_Comm local = MPI_COMM_NULL;
  HUNDUN_CHECK(MPI_Comm_split(mpi.comm(), mpi.rank(), 0, &local) ==
               MPI_SUCCESS);
  MPI_Comm intercommunicator = MPI_COMM_NULL;
  HUNDUN_CHECK(MPI_Intercomm_create(local, 0, mpi.comm(), 1 - mpi.rank(), 29,
                                    &intercommunicator) == MPI_SUCCESS);
  int is_intercommunicator = 0;
  HUNDUN_CHECK(MPI_Comm_test_inter(intercommunicator,
                                   &is_intercommunicator) == MPI_SUCCESS);
  HUNDUN_CHECK(is_intercommunicator != 0);
  expect_runtime_error([&] {
    static_cast<void>(MpiContext::duplicate(intercommunicator));
  });
  HUNDUN_CHECK(MPI_Comm_free(&intercommunicator) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_free(&local) == MPI_SUCCESS);
}

void run_decomposition_tests(const MpiContext& mpi,
                             const StructuredDecomposition& decomposition) {
  HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4 ||
               mpi.size() == 12);
  test_cartesian_topology(decomposition, mpi,
                          expected_invariant_grid(mpi.size()));
  test_boxes_and_coverage(decomposition, mpi);
  test_global_cells_and_ids(decomposition, mpi);
  test_neighbors(decomposition, mpi);
  test_frozen_shape_selection(mpi);
  test_legal_explicit_grid(mpi);
  test_collective_input_agreement(mpi);
  test_collective_selection_failures(mpi);
  test_move_ownership(mpi);
  test_extent_validation(mpi);
  test_intercommunicator_rejection(mpi);
  mpi.barrier();
}

}  // namespace

int main(int argc, char** argv) {
  std::optional<MpiContext> mpi;
  std::optional<StructuredDecomposition> decomposition;
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    mpi.emplace(MpiContext::duplicate(MPI_COMM_WORLD));
    decomposition.emplace(StructuredDecomposition::create(
        *mpi, kGlobalExtent, kPeriodic));
    const std::string_view mode = argc > 1 ? argv[1] : "full";
    result = hundun::test::run([&] {
      if (mode == "full") {
        run_decomposition_tests(*mpi, *decomposition);
      } else if (mode == "mismatch") {
        test_collective_input_agreement(*mpi);
      } else if (mode == "shape") {
        test_frozen_shape_selection(*mpi);
      } else {
        throw hundun::runtime::Error("unknown decomposition test mode");
      }
    });
  }

  const int finalized_result = hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);
    expect_runtime_error_with_text(
        [&] {
          static_cast<void>(StructuredDecomposition::create(
              *mpi, kGlobalExtent, kPeriodic));
        },
        "after MPI_Finalize");
    expect_runtime_error_with_text(
        [&] {
          static_cast<void>(decomposition->neighbor_rank(Int3{1, 0, 0}));
        },
        "after MPI_Finalize");
    decomposition.reset();
    mpi.reset();
  });

  return result == EXIT_SUCCESS ? finalized_result : result;
}
