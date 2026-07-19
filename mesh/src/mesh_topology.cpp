// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/mesh_topology.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace hundun::mesh {
namespace {

struct FaceCounts {
  std::uint64_t x{};
  std::uint64_t y{};
  std::uint64_t z{};
  std::uint64_t total{};
};

struct CellRecord {
  GlobalCellId global_id{};
  runtime::Int3 global_cell{};
  EntityOwnership ownership{EntityOwnership::ghost};
};

struct FaceRecord {
  GlobalFaceId global_id{};
  LogicalFace logical{};
  EntityOwnership ownership{EntityOwnership::ghost};
  LocalCellId owner{};
  std::optional<LocalCellId> neighbour{};
  std::optional<std::uint32_t> patch{};
  std::optional<GlobalFaceId> periodic_pair{};
};

struct FaceConnectivity {
  runtime::Int3 owner{};
  std::optional<runtime::Int3> neighbour{};
};

struct BuildData {
  runtime::Int3 global_extent{};
  runtime::Box3 owned_box{};
  std::array<bool, 3> periodic{};
  std::uint64_t global_cell_count{};
  FaceCounts face_counts{};
  std::size_t owned_cell_count{};
  std::size_t owned_face_count{};
  std::vector<CellRecord> cells{};
  std::vector<std::pair<GlobalCellId, LocalCellId>> cell_lookup{};
  std::vector<FaceRecord> faces{};
  std::array<std::vector<LocalFaceId>, 6> patch_faces{};
};

std::uint64_t checked_add(std::uint64_t lhs, std::uint64_t rhs,
                          const char* subject) {
  constexpr std::uint64_t limit = std::numeric_limits<std::uint64_t>::max();
  if (lhs > limit - rhs) {
    throw runtime::Error(std::string(subject) + " addition overflow");
  }
  return lhs + rhs;
}

std::uint64_t checked_multiply(std::uint64_t lhs, std::uint64_t rhs,
                               const char* subject) {
  constexpr std::uint64_t limit = std::numeric_limits<std::uint64_t>::max();
  if (rhs != 0U && lhs > limit / rhs) {
    throw runtime::Error(std::string(subject) + " multiplication overflow");
  }
  return lhs * rhs;
}

std::size_t checked_size(std::uint64_t value, const char* subject) {
  if (value >
      static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw runtime::Error(std::string(subject) +
                         " cannot be represented as a local size");
  }
  return static_cast<std::size_t>(value);
}

void validate_extent(runtime::Int3 extent) {
  if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0) {
    throw runtime::Error("mesh topology requires a positive global extent");
  }
}

bool valid_box(runtime::Box3 box, runtime::Int3 extent) {
  return box.begin.x >= 0 && box.begin.y >= 0 && box.begin.z >= 0 &&
         box.end.x <= extent.x && box.end.y <= extent.y &&
         box.end.z <= extent.z && box.begin.x < box.end.x &&
         box.begin.y < box.end.y && box.begin.z < box.end.z;
}

bool contains(runtime::Box3 box, runtime::Int3 cell) noexcept {
  return cell.x >= box.begin.x && cell.x < box.end.x && cell.y >= box.begin.y &&
         cell.y < box.end.y && cell.z >= box.begin.z && cell.z < box.end.z;
}

std::uint64_t global_cell_count(runtime::Int3 extent) {
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto nz = static_cast<std::uint64_t>(extent.z);
  return checked_multiply(checked_multiply(nx, ny, "global cell count"), nz,
                          "global cell count");
}

FaceCounts global_face_counts(runtime::Int3 extent) {
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto nz = static_cast<std::uint64_t>(extent.z);
  const auto nx_plus_one = checked_add(nx, 1U, "x-face extent");
  const auto ny_plus_one = checked_add(ny, 1U, "y-face extent");
  const auto nz_plus_one = checked_add(nz, 1U, "z-face extent");

  FaceCounts counts{};
  counts.x = checked_multiply(checked_multiply(nx_plus_one, ny, "x-face count"),
                              nz, "x-face count");
  counts.y = checked_multiply(checked_multiply(nx, ny_plus_one, "y-face count"),
                              nz, "y-face count");
  counts.z = checked_multiply(checked_multiply(nx, ny, "z-face count"),
                              nz_plus_one, "z-face count");
  counts.total =
      checked_add(checked_add(counts.x, counts.y, "global face count"),
                  counts.z, "global face count");
  return counts;
}

void validate_global_cell(runtime::Int3 extent, runtime::Int3 cell) {
  if (cell.x < 0 || cell.x >= extent.x || cell.y < 0 || cell.y >= extent.y ||
      cell.z < 0 || cell.z >= extent.z) {
    throw runtime::Error("global cell coordinate is outside the mesh");
  }
}

GlobalCellId calculate_cell_id(runtime::Int3 extent, runtime::Int3 cell) {
  validate_global_cell(extent, cell);
  const auto i = static_cast<std::uint64_t>(cell.x);
  const auto j = static_cast<std::uint64_t>(cell.y);
  const auto k = static_cast<std::uint64_t>(cell.z);
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const std::uint64_t row = checked_add(
      checked_multiply(k, ny, "global cell ID"), j, "global cell ID");
  return checked_add(checked_multiply(row, nx, "global cell ID"), i,
                     "global cell ID");
}

void validate_logical_face(runtime::Int3 extent, LogicalFace face) {
  const runtime::Int3 coordinate = face.coordinate;
  bool valid = false;
  switch (face.axis) {
    case FaceAxis::x:
      valid = coordinate.x >= 0 && coordinate.x <= extent.x &&
              coordinate.y >= 0 && coordinate.y < extent.y &&
              coordinate.z >= 0 && coordinate.z < extent.z;
      break;
    case FaceAxis::y:
      valid = coordinate.x >= 0 && coordinate.x < extent.x &&
              coordinate.y >= 0 && coordinate.y <= extent.y &&
              coordinate.z >= 0 && coordinate.z < extent.z;
      break;
    case FaceAxis::z:
      valid = coordinate.x >= 0 && coordinate.x < extent.x &&
              coordinate.y >= 0 && coordinate.y < extent.y &&
              coordinate.z >= 0 && coordinate.z <= extent.z;
      break;
    default:
      throw runtime::Error("invalid face axis");
  }
  if (!valid) {
    throw runtime::Error("logical face coordinate is outside its lattice");
  }
}

GlobalFaceId calculate_face_id(runtime::Int3 extent, const FaceCounts& counts,
                               LogicalFace face) {
  validate_logical_face(extent, face);
  const auto i = static_cast<std::uint64_t>(face.coordinate.x);
  const auto j = static_cast<std::uint64_t>(face.coordinate.y);
  const auto k = static_cast<std::uint64_t>(face.coordinate.z);
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  switch (face.axis) {
    case FaceAxis::x: {
      const std::uint64_t row =
          checked_add(checked_multiply(k, ny, "x-face ID"), j, "x-face ID");
      return checked_add(checked_multiply(row,
                                          checked_add(nx, 1U, "x-face extent"),
                                          "x-face ID"),
                         i, "x-face ID");
    }
    case FaceAxis::y: {
      const std::uint64_t row =
          checked_add(checked_multiply(k, checked_add(ny, 1U, "y-face extent"),
                                       "y-face ID"),
                      j, "y-face ID");
      const std::uint64_t local =
          checked_add(checked_multiply(row, nx, "y-face ID"), i, "y-face ID");
      return checked_add(counts.x, local, "y-face ID");
    }
    case FaceAxis::z: {
      const std::uint64_t row =
          checked_add(checked_multiply(k, ny, "z-face ID"), j, "z-face ID");
      const std::uint64_t local =
          checked_add(checked_multiply(row, nx, "z-face ID"), i, "z-face ID");
      return checked_add(checked_add(counts.x, counts.y, "z-face offset"),
                         local, "z-face ID");
    }
    default:
      throw runtime::Error("invalid face axis");
  }
}

std::array<LogicalFace, 6> cell_faces(runtime::Int3 cell) {
  return {LogicalFace{FaceAxis::x, runtime::Int3{cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::x, runtime::Int3{cell.x + 1, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, runtime::Int3{cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, runtime::Int3{cell.x, cell.y + 1, cell.z}},
          LogicalFace{FaceAxis::z, runtime::Int3{cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::z, runtime::Int3{cell.x, cell.y, cell.z + 1}}};
}

FaceConnectivity face_connectivity(runtime::Int3 extent,
                                   std::array<bool, 3> periodic,
                                   LogicalFace face) {
  validate_logical_face(extent, face);
  runtime::Int3 lower = face.coordinate;
  runtime::Int3 upper = face.coordinate;
  int plane = 0;
  int cells = 0;
  bool wraps = false;
  switch (face.axis) {
    case FaceAxis::x:
      plane = face.coordinate.x;
      cells = extent.x;
      wraps = periodic[0];
      --lower.x;
      break;
    case FaceAxis::y:
      plane = face.coordinate.y;
      cells = extent.y;
      wraps = periodic[1];
      --lower.y;
      break;
    case FaceAxis::z:
      plane = face.coordinate.z;
      cells = extent.z;
      wraps = periodic[2];
      --lower.z;
      break;
    default:
      throw runtime::Error("invalid face axis");
  }

  if (plane > 0 && plane < cells) {
    return FaceConnectivity{lower, upper};
  }
  runtime::Int3 owner = plane == 0 ? upper : lower;
  if (!wraps) {
    return FaceConnectivity{owner, std::nullopt};
  }

  runtime::Int3 opposite = owner;
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
    default:
      throw runtime::Error("invalid face axis");
  }
  return FaceConnectivity{owner, opposite};
}

std::optional<std::uint32_t> boundary_patch(runtime::Int3 extent,
                                            LogicalFace face) {
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
    default:
      throw runtime::Error("invalid face axis");
  }
}

std::optional<GlobalFaceId> periodic_pair(runtime::Int3 extent,
                                          const FaceCounts& counts,
                                          std::array<bool, 3> periodic,
                                          LogicalFace face) {
  const auto patch = boundary_patch(extent, face);
  if (!patch.has_value()) {
    return std::nullopt;
  }
  const std::size_t axis = static_cast<std::size_t>(*patch / 2U);
  if (!periodic[axis]) {
    return std::nullopt;
  }
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
    default:
      throw runtime::Error("invalid face axis");
  }
  return calculate_face_id(extent, counts, face);
}

std::optional<LocalCellId>
find_cell(const std::vector<std::pair<GlobalCellId, LocalCellId>>& lookup,
          GlobalCellId global_id) {
  const auto found = std::lower_bound(
      lookup.begin(), lookup.end(), global_id,
      [](const auto& entry, GlobalCellId id) { return entry.first < id; });
  if (found == lookup.end() || found->first != global_id) {
    return std::nullopt;
  }
  return found->second;
}

std::size_t local_volume(runtime::Box3 box) {
  const auto nx = static_cast<std::uint64_t>(box.end.x - box.begin.x);
  const auto ny = static_cast<std::uint64_t>(box.end.y - box.begin.y);
  const auto nz = static_cast<std::uint64_t>(box.end.z - box.begin.z);
  return checked_size(
      checked_multiply(checked_multiply(nx, ny, "owned cell count"), nz,
                       "owned cell count"),
      "owned cell count");
}

BuildData build_data(const runtime::StructuredDecomposition& decomposition) {
  BuildData data{};
  data.global_extent = decomposition.global_extent();
  data.owned_box = decomposition.owned_box();
  data.periodic = decomposition.periodic();
  validate_extent(data.global_extent);
  if (!valid_box(data.owned_box, data.global_extent)) {
    throw runtime::Error("mesh topology received an invalid owned box");
  }

  data.global_cell_count = global_cell_count(data.global_extent);
  data.face_counts = global_face_counts(data.global_extent);
  static_cast<void>(checked_size(data.global_cell_count, "global cell count"));
  static_cast<void>(checked_size(data.face_counts.total, "global face count"));
  data.owned_cell_count = local_volume(data.owned_box);

  data.cells.reserve(data.owned_cell_count);
  for (int k = data.owned_box.begin.z; k < data.owned_box.end.z; ++k) {
    for (int j = data.owned_box.begin.y; j < data.owned_box.end.y; ++j) {
      for (int i = data.owned_box.begin.x; i < data.owned_box.end.x; ++i) {
        const runtime::Int3 cell{i, j, k};
        data.cells.push_back(
            CellRecord{calculate_cell_id(data.global_extent, cell), cell,
                       EntityOwnership::owned});
      }
    }
  }

  const std::uint64_t candidate_count =
      checked_multiply(static_cast<std::uint64_t>(data.owned_cell_count), 6U,
                       "local face candidate count");
  std::vector<LogicalFace> face_candidates;
  face_candidates.reserve(
      checked_size(candidate_count, "local face candidate count"));
  for (std::size_t local = 0; local < data.owned_cell_count; ++local) {
    const auto faces = cell_faces(data.cells[local].global_cell);
    face_candidates.insert(face_candidates.end(), faces.begin(), faces.end());
  }
  const auto face_less = [&](LogicalFace lhs, LogicalFace rhs) {
    return calculate_face_id(data.global_extent, data.face_counts, lhs) <
           calculate_face_id(data.global_extent, data.face_counts, rhs);
  };
  std::sort(face_candidates.begin(), face_candidates.end(), face_less);
  face_candidates.erase(
      std::unique(face_candidates.begin(), face_candidates.end(),
                  [&](LogicalFace lhs, LogicalFace rhs) {
                    return calculate_face_id(data.global_extent,
                                             data.face_counts, lhs) ==
                           calculate_face_id(data.global_extent,
                                             data.face_counts, rhs);
                  }),
      face_candidates.end());

  std::vector<runtime::Int3> ghost_candidates;
  ghost_candidates.reserve(face_candidates.size());
  for (const LogicalFace face : face_candidates) {
    const FaceConnectivity connection =
        face_connectivity(data.global_extent, data.periodic, face);
    if (!contains(data.owned_box, connection.owner)) {
      ghost_candidates.push_back(connection.owner);
    }
    if (connection.neighbour.has_value() &&
        !contains(data.owned_box, *connection.neighbour)) {
      ghost_candidates.push_back(*connection.neighbour);
    }
  }
  const auto cell_less = [&](runtime::Int3 lhs, runtime::Int3 rhs) {
    return calculate_cell_id(data.global_extent, lhs) <
           calculate_cell_id(data.global_extent, rhs);
  };
  std::sort(ghost_candidates.begin(), ghost_candidates.end(), cell_less);
  ghost_candidates.erase(
      std::unique(ghost_candidates.begin(), ghost_candidates.end(),
                  [&](runtime::Int3 lhs, runtime::Int3 rhs) {
                    return calculate_cell_id(data.global_extent, lhs) ==
                           calculate_cell_id(data.global_extent, rhs);
                  }),
      ghost_candidates.end());

  if (ghost_candidates.size() >
      std::numeric_limits<std::size_t>::max() - data.cells.size()) {
    throw runtime::Error("local cell count addition overflow");
  }
  data.cells.reserve(data.cells.size() + ghost_candidates.size());
  for (const runtime::Int3 cell : ghost_candidates) {
    data.cells.push_back(CellRecord{calculate_cell_id(data.global_extent, cell),
                                    cell, EntityOwnership::ghost});
  }

  data.cell_lookup.reserve(data.cells.size());
  for (LocalCellId local = 0; local < data.cells.size(); ++local) {
    data.cell_lookup.emplace_back(data.cells[local].global_id, local);
  }
  std::sort(
      data.cell_lookup.begin(), data.cell_lookup.end(),
      [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

  data.faces.reserve(face_candidates.size());
  for (const LogicalFace face : face_candidates) {
    const FaceConnectivity connection =
        face_connectivity(data.global_extent, data.periodic, face);
    const auto owner =
        find_cell(data.cell_lookup,
                  calculate_cell_id(data.global_extent, connection.owner));
    if (!owner.has_value()) {
      throw runtime::Error("local face owner is missing from compact cells");
    }
    std::optional<LocalCellId> neighbour;
    if (connection.neighbour.has_value()) {
      neighbour =
          find_cell(data.cell_lookup, calculate_cell_id(data.global_extent,
                                                        *connection.neighbour));
      if (!neighbour.has_value()) {
        throw runtime::Error(
            "local face neighbour is missing from compact cells");
      }
    }
    const EntityOwnership ownership = contains(data.owned_box, connection.owner)
                                          ? EntityOwnership::owned
                                          : EntityOwnership::ghost;
    const std::optional<std::uint32_t> patch =
        boundary_patch(data.global_extent, face);
    const LocalFaceId local_face = data.faces.size();
    data.faces.push_back(FaceRecord{
        calculate_face_id(data.global_extent, data.face_counts, face), face,
        ownership, *owner, neighbour, patch,
        periodic_pair(data.global_extent, data.face_counts, data.periodic,
                      face)});
    if (ownership == EntityOwnership::owned) {
      ++data.owned_face_count;
      if (patch.has_value()) {
        data.patch_faces[*patch].push_back(local_face);
      }
    }
  }
  return data;
}

}  // namespace

struct MeshTopology::Impl {
  Impl(BuildData data, std::array<BoundaryPatch, 6> built_patches)
      : global_extent(data.global_extent), owned_box(data.owned_box),
        periodic(data.periodic), global_cell_count(data.global_cell_count),
        face_counts(data.face_counts), owned_cell_count(data.owned_cell_count),
        owned_face_count(data.owned_face_count), cells(std::move(data.cells)),
        cell_lookup(std::move(data.cell_lookup)), faces(std::move(data.faces)),
        patches(std::move(built_patches)) {}

  runtime::Int3 global_extent{};
  runtime::Box3 owned_box{};
  std::array<bool, 3> periodic{};
  std::uint64_t global_cell_count{};
  FaceCounts face_counts{};
  std::size_t owned_cell_count{};
  std::size_t owned_face_count{};
  std::vector<CellRecord> cells{};
  std::vector<std::pair<GlobalCellId, LocalCellId>> cell_lookup{};
  std::vector<FaceRecord> faces{};
  std::array<BoundaryPatch, 6> patches;
};

BoundaryPatch::BoundaryPatch(std::uint32_t stable_id, std::string_view name,
                             PatchPairingKind pairing_kind,
                             std::optional<std::uint32_t> paired_patch_id,
                             std::vector<LocalFaceId> local_faces)
    : stable_id_(stable_id), name_(name), pairing_kind_(pairing_kind),
      paired_patch_id_(paired_patch_id), local_faces_(std::move(local_faces)) {}

std::uint32_t BoundaryPatch::stable_id() const noexcept { return stable_id_; }

std::string_view BoundaryPatch::name() const noexcept { return name_; }

PatchPairingKind BoundaryPatch::pairing_kind() const noexcept {
  return pairing_kind_;
}

std::optional<std::uint32_t> BoundaryPatch::paired_patch_id() const noexcept {
  return paired_patch_id_;
}

const std::vector<LocalFaceId>& BoundaryPatch::local_faces() const noexcept {
  return local_faces_;
}

bool BoundaryPatch::contains(LocalFaceId local_face) const noexcept {
  return std::binary_search(local_faces_.begin(), local_faces_.end(),
                            local_face);
}

std::unique_ptr<MeshTopology::Impl> MeshTopology::create_impl(
    const runtime::StructuredDecomposition& decomposition) {
  BuildData data = build_data(decomposition);
  const auto pairing_kind = [&](std::size_t axis) {
    return data.periodic[axis] ? PatchPairingKind::periodic
                               : PatchPairingKind::none;
  };
  const auto paired_id = [&](std::size_t axis,
                             std::uint32_t id) -> std::optional<std::uint32_t> {
    if (!data.periodic[axis]) {
      return std::nullopt;
    }
    return id % 2U == 0U ? id + 1U : id - 1U;
  };

  std::array<BoundaryPatch, 6> patches{
      BoundaryPatch{0U, "x_min", pairing_kind(0U), paired_id(0U, 0U),
                    std::move(data.patch_faces[0])},
      BoundaryPatch{1U, "x_max", pairing_kind(0U), paired_id(0U, 1U),
                    std::move(data.patch_faces[1])},
      BoundaryPatch{2U, "y_min", pairing_kind(1U), paired_id(1U, 2U),
                    std::move(data.patch_faces[2])},
      BoundaryPatch{3U, "y_max", pairing_kind(1U), paired_id(1U, 3U),
                    std::move(data.patch_faces[3])},
      BoundaryPatch{4U, "z_min", pairing_kind(2U), paired_id(2U, 4U),
                    std::move(data.patch_faces[4])},
      BoundaryPatch{5U, "z_max", pairing_kind(2U), paired_id(2U, 5U),
                    std::move(data.patch_faces[5])}};
  return std::make_unique<Impl>(std::move(data), std::move(patches));
}

MeshTopology::MeshTopology(
    const runtime::StructuredDecomposition& decomposition) try
    : impl_(create_impl(decomposition)) {
} catch (const runtime::Error&) {
  throw;
} catch (const std::bad_alloc&) {
  throw runtime::Error("mesh topology allocation failed");
} catch (const std::length_error&) {
  throw runtime::Error("mesh topology allocation size is unsupported");
}

MeshTopology::~MeshTopology() = default;
MeshTopology::MeshTopology(MeshTopology&&) noexcept = default;
MeshTopology& MeshTopology::operator=(MeshTopology&&) noexcept = default;

runtime::Int3 MeshTopology::global_extent() const noexcept {
  return impl_->global_extent;
}

runtime::Box3 MeshTopology::owned_global_box() const noexcept {
  return impl_->owned_box;
}

std::uint64_t MeshTopology::global_cell_count() const noexcept {
  return impl_->global_cell_count;
}

std::uint64_t MeshTopology::global_face_count() const noexcept {
  return impl_->face_counts.total;
}

std::uint64_t MeshTopology::global_face_count(FaceAxis axis) const {
  switch (axis) {
    case FaceAxis::x:
      return impl_->face_counts.x;
    case FaceAxis::y:
      return impl_->face_counts.y;
    case FaceAxis::z:
      return impl_->face_counts.z;
    default:
      throw runtime::Error("invalid face axis");
  }
}

std::size_t MeshTopology::owned_cell_count() const noexcept {
  return impl_->owned_cell_count;
}

std::size_t MeshTopology::ghost_cell_count() const noexcept {
  return impl_->cells.size() - impl_->owned_cell_count;
}

std::size_t MeshTopology::local_cell_count() const noexcept {
  return impl_->cells.size();
}

std::size_t MeshTopology::owned_face_count() const noexcept {
  return impl_->owned_face_count;
}

std::size_t MeshTopology::ghost_face_count() const noexcept {
  return impl_->faces.size() - impl_->owned_face_count;
}

std::size_t MeshTopology::local_face_count() const noexcept {
  return impl_->faces.size();
}

GlobalCellId MeshTopology::global_cell_id(runtime::Int3 global_cell) const {
  return calculate_cell_id(impl_->global_extent, global_cell);
}

template <class ImplType>
const CellRecord& require_cell(const ImplType& impl, LocalCellId local_cell) {
  if (local_cell >= impl.cells.size()) {
    throw runtime::Error("local cell index is outside compact topology");
  }
  return impl.cells[local_cell];
}

GlobalCellId MeshTopology::global_cell_id(LocalCellId local_cell) const {
  return require_cell(*impl_, local_cell).global_id;
}

runtime::Int3 MeshTopology::global_cell(LocalCellId local_cell) const {
  return require_cell(*impl_, local_cell).global_cell;
}

EntityOwnership MeshTopology::cell_ownership(LocalCellId local_cell) const {
  return require_cell(*impl_, local_cell).ownership;
}

std::optional<LocalCellId>
MeshTopology::find_local_cell(GlobalCellId global_id) const {
  if (global_id >= impl_->global_cell_count) {
    throw runtime::Error("global cell ID is outside the mesh");
  }
  return find_cell(impl_->cell_lookup, global_id);
}

GlobalFaceId MeshTopology::global_face_id(LogicalFace face) const {
  return calculate_face_id(impl_->global_extent, impl_->face_counts, face);
}

template <class ImplType>
const FaceRecord& require_face(const ImplType& impl, LocalFaceId local_face) {
  if (local_face >= impl.faces.size()) {
    throw runtime::Error("local face index is outside compact topology");
  }
  return impl.faces[local_face];
}

GlobalFaceId MeshTopology::global_face_id(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).global_id;
}

LogicalFace MeshTopology::logical_face(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).logical;
}

EntityOwnership MeshTopology::face_ownership(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).ownership;
}

LocalCellId MeshTopology::owner(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).owner;
}

std::optional<LocalCellId>
MeshTopology::neighbour(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).neighbour;
}

std::optional<LocalFaceId>
MeshTopology::find_local_face(GlobalFaceId global_id) const {
  if (global_id >= impl_->face_counts.total) {
    throw runtime::Error("global face ID is outside the mesh");
  }
  const auto found =
      std::lower_bound(impl_->faces.begin(), impl_->faces.end(), global_id,
                       [](const FaceRecord& face, GlobalFaceId id) {
                         return face.global_id < id;
                       });
  if (found == impl_->faces.end() || found->global_id != global_id) {
    return std::nullopt;
  }
  return static_cast<LocalFaceId>(found - impl_->faces.begin());
}

std::optional<std::uint32_t>
MeshTopology::patch_id(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).patch;
}

std::optional<GlobalFaceId>
MeshTopology::periodic_pair(LocalFaceId local_face) const {
  return require_face(*impl_, local_face).periodic_pair;
}

const std::array<BoundaryPatch, 6>& MeshTopology::patches() const noexcept {
  return impl_->patches;
}

const BoundaryPatch& MeshTopology::patch(std::uint32_t stable_id) const {
  if (stable_id >= impl_->patches.size()) {
    throw runtime::Error("boundary patch stable ID is outside [0, 6)");
  }
  return impl_->patches[stable_id];
}

}  // namespace hundun::mesh
