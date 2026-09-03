// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_domain.hpp"

#include "ib_periodic_surface_window_detail.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace hundun::immersed::detail {

struct ActiveCellLayoutStorage final {
  std::size_t local_cell_count{};
  std::size_t owned_active_count{};
  std::vector<mesh::GlobalCellId> ordered_global_ids;
  std::vector<std::size_t> active_indices;
  std::uint64_t fingerprint{};
};

struct ActiveBoundaryLayoutStorage final {
  std::array<std::vector<mesh::GlobalFaceId>, 6> patch_faces;
  bool open_domain{};
  bool has_pressure_reference{};
  std::uint64_t fingerprint{};
};

struct ImmersedDomainStorage final {
  ImmersedDomainStorage(std::vector<CellRegion> regions_value,
                        std::vector<ImmersedLink> links_value,
                        ActiveCellLayout active_cells_value,
                        ActiveBoundaryLayout active_boundaries_value,
                        std::uint64_t classification_fingerprint_value,
                        std::uint64_t surface_coverage_fingerprint_value,
                        std::uint64_t classified_cells_value)
      : regions(std::move(regions_value)), links(std::move(links_value)),
        active_cells(std::move(active_cells_value)),
        active_boundaries(std::move(active_boundaries_value)),
        classification_fingerprint(classification_fingerprint_value),
        surface_coverage_fingerprint(surface_coverage_fingerprint_value) {
    performance.init_classification_cells = classified_cells_value;
  }

  std::vector<CellRegion> regions;
  std::vector<ImmersedLink> links;
  ActiveCellLayout active_cells;
  ActiveBoundaryLayout active_boundaries;
  std::uint64_t classification_fingerprint{};
  std::uint64_t surface_coverage_fingerprint{};
  diagnostics::Stage3PerformanceCounters performance;
};

} // namespace hundun::immersed::detail

namespace hundun::immersed {
namespace {

constexpr std::uint64_t fnv_offset = UINT64_C(14695981039346656037);
constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);
constexpr std::size_t no_active_index = std::numeric_limits<std::size_t>::max();
constexpr std::size_t link_record_width = 22U;
constexpr std::size_t boundary_record_width = 3U;
constexpr std::size_t boundary_signature_patch_width = 9U;
using BoundarySignature =
    std::array<std::uint64_t, 3U + 6U * boundary_signature_patch_width>;

struct LinkKey final {
  mesh::GlobalCellId first{};
  mesh::GlobalCellId second{};
  TriangleId triangle{};
  std::array<std::uint64_t, 3> intercept_bits{};
};

struct LinkRecord final {
  mesh::GlobalFaceId face{};
  bool owned{};
  LinkKey key{};
  mesh::GlobalCellId fluid_cell{};
  mesh::GlobalCellId solid_cell{};
  runtime::Real3 intercept{};
  runtime::Real3 normal{};
  double fraction{};
  runtime::Real3 fluid_center{};
  runtime::Real3 solid_center{};
  double fluid_h{};
  double solid_h{};
};

struct CanonicalLink final {
  ImmersedLinkId id{};
  LinkRecord record{};
};

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    hash ^= (value >> (8U * byte)) & UINT64_C(0xff);
    hash *= fnv_prime;
  }
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t result = 0U;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double double_from_bits(std::uint64_t value) noexcept {
  double result = 0.0;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

runtime::Real3 subtract(runtime::Real3 first, runtime::Real3 second) noexcept {
  return {first.x - second.x, first.y - second.y, first.z - second.z};
}

runtime::Real3 multiply(runtime::Real3 value, double factor) noexcept {
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(runtime::Real3 first, runtime::Real3 second) noexcept {
  return first.x * second.x + first.y * second.y + first.z * second.z;
}

double norm(runtime::Real3 value) noexcept {
  return std::hypot(value.x, value.y, value.z);
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

std::string mpi_error_message(int code, std::string_view operation) {
  std::array<char, MPI_MAX_ERROR_STRING> buffer{};
  int length = 0;
  if (MPI_Error_string(code, buffer.data(), &length) != MPI_SUCCESS ||
      length < 0 || static_cast<std::size_t>(length) > buffer.size()) {
    return std::string(operation) + " failed";
  }
  return std::string(operation) + " failed: " +
         std::string(buffer.data(), static_cast<std::size_t>(length));
}

void check_mpi(int code, std::string_view operation) {
  if (code != MPI_SUCCESS) {
    throw runtime::Error(mpi_error_message(code, operation));
  }
}

void require_collective(const runtime::MpiContext &mpi, bool local_ok,
                        std::string_view invariant) {
  const runtime::CollectiveStatus status =
      runtime::collective_status(mpi, local_ok, invariant);
  if (!status.ok) {
    throw runtime::Error(status.message + " (lowest failing rank " +
                         std::to_string(status.failing_rank) + ")");
  }
}

std::vector<std::uint64_t>
gather_u64_chunked(const std::vector<std::uint64_t> &local,
                   const runtime::MpiContext &mpi) {
  require_collective(
      mpi,
      local.size() <=
          static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max()),
      "immersed domain: invariant local_record_size failed");
  bool sizes_allocation_ok = true;
  std::vector<std::uint64_t> sizes;
  try {
    sizes.resize(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    sizes_allocation_ok = false;
  }
  require_collective(
      mpi, sizes_allocation_ok,
      "immersed domain: invariant collective_size_allocation failed");

  for (int root = 0; root < mpi.size(); ++root) {
    std::uint64_t size =
        mpi.rank() == root ? static_cast<std::uint64_t>(local.size()) : 0U;
    check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, root, mpi.comm()),
              "MPI_Bcast immersed-domain record count");
    sizes[static_cast<std::size_t>(root)] = size;
  }

  std::size_t total = 0U;
  bool representable = true;
  for (const std::uint64_t size : sizes) {
    if (size > std::numeric_limits<std::size_t>::max() ||
        static_cast<std::size_t>(size) >
            std::numeric_limits<std::size_t>::max() - total) {
      representable = false;
      break;
    }
    total += static_cast<std::size_t>(size);
  }
  representable =
      representable && total <= static_cast<std::size_t>(
                                    std::numeric_limits<std::ptrdiff_t>::max());
  require_collective(
      mpi, representable,
      "immersed domain: invariant collective_record_size failed");

  bool allocation_ok = true;
  std::vector<std::uint64_t> gathered;
  try {
    gathered.resize(total);
  } catch (...) {
    allocation_ok = false;
  }
  require_collective(
      mpi, allocation_ok,
      "immersed domain: invariant collective_record_allocation failed");

  std::size_t offset = 0U;
  constexpr std::size_t maximum_chunk =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  for (int root = 0; root < mpi.size(); ++root) {
    const std::size_t root_size =
        static_cast<std::size_t>(sizes[static_cast<std::size_t>(root)]);
    if (mpi.rank() == root) {
      std::copy(local.begin(), local.end(),
                gathered.begin() + static_cast<std::ptrdiff_t>(offset));
    }
    std::size_t sent = 0U;
    while (sent < root_size) {
      const std::size_t chunk = std::min(maximum_chunk, root_size - sent);
      check_mpi(MPI_Bcast(gathered.data() + offset + sent,
                          static_cast<int>(chunk), MPI_UINT64_T, root,
                          mpi.comm()),
                "MPI_Bcast immersed-domain records");
      sent += chunk;
    }
    offset += root_size;
  }
  return gathered;
}

bool less_key(const LinkKey &first, const LinkKey &second) noexcept {
  return std::tie(first.first, first.second, first.triangle,
                  first.intercept_bits) < std::tie(second.first, second.second,
                                                   second.triangle,
                                                   second.intercept_bits);
}

bool same_key(const LinkKey &first, const LinkKey &second) noexcept {
  return !less_key(first, second) && !less_key(second, first);
}

bool same_record(const LinkRecord &first, const LinkRecord &second) noexcept {
  return first.face == second.face && same_key(first.key, second.key) &&
         first.fluid_cell == second.fluid_cell &&
         first.solid_cell == second.solid_cell &&
         double_bits(first.intercept.x) == double_bits(second.intercept.x) &&
         double_bits(first.intercept.y) == double_bits(second.intercept.y) &&
         double_bits(first.intercept.z) == double_bits(second.intercept.z) &&
         double_bits(first.normal.x) == double_bits(second.normal.x) &&
         double_bits(first.normal.y) == double_bits(second.normal.y) &&
         double_bits(first.normal.z) == double_bits(second.normal.z) &&
         double_bits(first.fraction) == double_bits(second.fraction) &&
         double_bits(first.fluid_center.x) ==
             double_bits(second.fluid_center.x) &&
         double_bits(first.fluid_center.y) ==
             double_bits(second.fluid_center.y) &&
         double_bits(first.fluid_center.z) ==
             double_bits(second.fluid_center.z) &&
         double_bits(first.solid_center.x) ==
             double_bits(second.solid_center.x) &&
         double_bits(first.solid_center.y) ==
             double_bits(second.solid_center.y) &&
         double_bits(first.solid_center.z) ==
             double_bits(second.solid_center.z) &&
         double_bits(first.fluid_h) == double_bits(second.fluid_h) &&
         double_bits(first.solid_h) == double_bits(second.solid_h);
}

void append_link_record(std::vector<std::uint64_t> &encoded,
                        const LinkRecord &record) {
  encoded.insert(encoded.end(), {record.face,
                                 record.owned ? 1U : 0U,
                                 record.key.first,
                                 record.key.second,
                                 record.fluid_cell,
                                 record.solid_cell,
                                 record.key.triangle,
                                 double_bits(record.intercept.x),
                                 double_bits(record.intercept.y),
                                 double_bits(record.intercept.z),
                                 double_bits(record.normal.x),
                                 double_bits(record.normal.y),
                                 double_bits(record.normal.z),
                                 double_bits(record.fraction),
                                 double_bits(record.fluid_center.x),
                                 double_bits(record.fluid_center.y),
                                 double_bits(record.fluid_center.z),
                                 double_bits(record.solid_center.x),
                                 double_bits(record.solid_center.y),
                                 double_bits(record.solid_center.z),
                                 double_bits(record.fluid_h),
                                 double_bits(record.solid_h)});
}

LinkRecord decode_link_record(const std::uint64_t *record) {
  LinkRecord result;
  result.face = record[0];
  result.owned = record[1] != 0U;
  result.key = {
      record[2], record[3], record[6], {record[7], record[8], record[9]}};
  result.fluid_cell = record[4];
  result.solid_cell = record[5];
  result.intercept = {double_from_bits(record[7]), double_from_bits(record[8]),
                      double_from_bits(record[9])};
  result.normal = {double_from_bits(record[10]), double_from_bits(record[11]),
                   double_from_bits(record[12])};
  result.fraction = double_from_bits(record[13]);
  result.fluid_center = {double_from_bits(record[14]),
                         double_from_bits(record[15]),
                         double_from_bits(record[16])};
  result.solid_center = {double_from_bits(record[17]),
                         double_from_bits(record[18]),
                         double_from_bits(record[19])};
  result.fluid_h = double_from_bits(record[20]);
  result.solid_h = double_from_bits(record[21]);
  return result;
}

std::uint64_t compute_classification_fingerprint(
    const ImmersedSurface &surface, const SurfaceQuery &query,
    config::ImmersedFluidSide fluid_side,
    const std::vector<CellRegion> &global_regions) noexcept {
  std::uint64_t hash = fnv_offset;
  hash_u64(hash, UINT64_C(0x48554e44434c5333));
  hash_u64(hash, surface.fingerprint());
  hash_u64(hash, query.fingerprint());
  hash_u64(hash, static_cast<std::uint64_t>(fluid_side));
  hash_u64(hash, global_regions.size());
  for (std::size_t id = 0U; id < global_regions.size(); ++id) {
    hash_u64(hash, id);
    hash_u64(hash, static_cast<std::uint64_t>(global_regions[id]));
  }
  return hash;
}

std::shared_ptr<const detail::ActiveCellLayoutStorage>
build_active_cells(const mesh::MeshTopology &topology,
                   const std::vector<CellRegion> &regions) {
  auto storage = std::make_shared<detail::ActiveCellLayoutStorage>();
  storage->local_cell_count = topology.local_cell_count();
  storage->active_indices.assign(storage->local_cell_count, no_active_index);

  std::vector<std::pair<mesh::GlobalCellId, mesh::LocalCellId>> owned;
  std::vector<std::pair<mesh::GlobalCellId, mesh::LocalCellId>> ghosts;
  for (mesh::LocalCellId local = 0U; local < topology.local_cell_count();
       ++local) {
    if (regions[local] != CellRegion::fluid) {
      continue;
    }
    auto &destination =
        topology.cell_ownership(local) == mesh::EntityOwnership::owned ? owned
                                                                       : ghosts;
    destination.push_back({topology.global_cell_id(local), local});
  }
  std::sort(owned.begin(), owned.end());
  std::sort(ghosts.begin(), ghosts.end());
  storage->owned_active_count = owned.size();
  storage->ordered_global_ids.reserve(owned.size() + ghosts.size());
  const auto append = [&](const auto &source) {
    for (const auto &[global, local] : source) {
      const std::size_t active = storage->ordered_global_ids.size();
      storage->ordered_global_ids.push_back(global);
      storage->active_indices[local] = active;
    }
  };
  append(owned);
  append(ghosts);

  std::uint64_t hash = fnv_offset;
  hash_u64(hash, UINT64_C(0x48554e4441435433));
  hash_u64(hash, storage->local_cell_count);
  hash_u64(hash, storage->owned_active_count);
  for (const mesh::GlobalCellId id : storage->ordered_global_ids) {
    hash_u64(hash, id);
  }
  storage->fingerprint = hash;
  return storage;
}

bool make_boundary_signature(const mesh::MeshTopology &topology,
                             const boundary::BoundaryRegistry &boundaries,
                             BoundarySignature &signature) noexcept {
  constexpr std::uint64_t absent = std::numeric_limits<std::uint64_t>::max();
  bool compatible = true;
  try {
    std::optional<std::uint32_t> inlet;
    std::optional<std::uint32_t> outlet;
    for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
      const boundary::BoundaryDescriptor &descriptor = boundaries.patch(patch);
      const auto topology_pair = topology.patch(patch).paired_patch_id();
      const std::size_t offset =
          3U + static_cast<std::size_t>(patch) * boundary_signature_patch_width;
      signature[offset] = descriptor.stable_id();
      signature[offset + 1U] = static_cast<std::uint64_t>(descriptor.kind());
      signature[offset + 2U] =
          static_cast<std::uint64_t>(descriptor.velocity_rule());
      signature[offset + 3U] =
          static_cast<std::uint64_t>(descriptor.pressure_rule());
      signature[offset + 4U] =
          static_cast<std::uint64_t>(descriptor.density_rule());
      signature[offset + 5U] =
          static_cast<std::uint64_t>(descriptor.enthalpy_rule());
      signature[offset + 6U] =
          static_cast<std::uint64_t>(descriptor.scalar_rule());
      signature[offset + 7U] =
          static_cast<std::uint64_t>(descriptor.mass_flux_rule());
      signature[offset + 8U] = descriptor.paired_patch_id().has_value()
                                   ? *descriptor.paired_patch_id()
                                   : absent;
      compatible = compatible && descriptor.stable_id() == patch &&
                   descriptor.paired_patch_id() == topology_pair &&
                   ((descriptor.kind() == boundary::BoundaryKind::periodic) ==
                    topology_pair.has_value());
      if (descriptor.kind() == boundary::BoundaryKind::velocity_inlet) {
        compatible = compatible && !inlet.has_value();
        inlet = patch;
      } else if (descriptor.kind() == boundary::BoundaryKind::pressure_outlet) {
        compatible = compatible && !outlet.has_value();
        outlet = patch;
      }
    }
    compatible = compatible && inlet.has_value() == outlet.has_value();
    signature[0] = inlet.has_value() ? 1U : 0U;
    signature[1] = inlet.has_value() ? *inlet : absent;
    signature[2] = outlet.has_value() ? *outlet : absent;
  } catch (...) {
    compatible = false;
  }
  return compatible;
}

std::shared_ptr<const detail::ActiveBoundaryLayoutStorage>
build_active_boundaries(const mesh::MeshTopology &topology,
                        const boundary::BoundaryRegistry &boundaries,
                        const std::vector<CellRegion> &regions,
                        const runtime::MpiContext &mpi) {
  BoundarySignature boundary_signature{};
  require_collective(
      mpi, make_boundary_signature(topology, boundaries, boundary_signature),
      "immersed domain: invariant boundary_registry_compatible failed");
  BoundarySignature root_signature = boundary_signature;
  check_mpi(MPI_Bcast(root_signature.data(),
                      static_cast<int>(root_signature.size()), MPI_UINT64_T, 0,
                      mpi.comm()),
            "MPI_Bcast immersed-domain boundary signature");
  require_collective(
      mpi, boundary_signature == root_signature,
      "immersed domain: invariant boundary_registry_consistent failed");

  std::vector<std::uint64_t> local;
  bool local_records_ok = true;
  try {
    for (std::uint32_t patch_id = 0U; patch_id < 6U; ++patch_id) {
      for (const mesh::LocalFaceId face :
           topology.patch(patch_id).local_faces()) {
        if (regions[topology.owner(face)] != CellRegion::fluid) {
          continue;
        }
        const auto pair = topology.periodic_pair(face);
        local.insert(
            local.end(),
            {patch_id, topology.global_face_id(face),
             pair.value_or(std::numeric_limits<std::uint64_t>::max())});
      }
    }
  } catch (...) {
    local_records_ok = false;
  }
  require_collective(
      mpi, local_records_ok,
      "immersed domain: invariant active_boundary_record_allocation failed");
  const auto gathered = gather_u64_chunked(local, mpi);
  require_collective(
      mpi, gathered.size() % boundary_record_width == 0U,
      "immersed domain: invariant active_boundary_record_size failed");

  struct BoundaryRecord final {
    std::uint32_t patch{};
    mesh::GlobalFaceId face{};
    std::optional<mesh::GlobalFaceId> pair;
  };
  std::vector<BoundaryRecord> records;
  bool unique = true;
  bool patch_ids_ok = true;
  bool materialization_ok = true;
  std::shared_ptr<detail::ActiveBoundaryLayoutStorage> storage;
  try {
    records.reserve(gathered.size() / boundary_record_width);
    for (std::size_t offset = 0U; offset < gathered.size();
         offset += boundary_record_width) {
      const std::uint64_t raw_patch = gathered[offset];
      const std::uint64_t raw_pair = gathered[offset + 2U];
      if (raw_patch >= 6U) {
        patch_ids_ok = false;
        break;
      }
      records.push_back({static_cast<std::uint32_t>(raw_patch),
                         gathered[offset + 1U],
                         raw_pair == std::numeric_limits<std::uint64_t>::max()
                             ? std::nullopt
                             : std::optional<mesh::GlobalFaceId>{raw_pair}});
    }
    std::sort(records.begin(), records.end(),
              [](const BoundaryRecord &first, const BoundaryRecord &second) {
                return std::tie(first.patch, first.face) <
                       std::tie(second.patch, second.face);
              });
    for (std::size_t index = 1U; index < records.size(); ++index) {
      if (records[index - 1U].patch == records[index].patch &&
          records[index - 1U].face == records[index].face) {
        unique = false;
      }
    }
    storage = std::make_shared<detail::ActiveBoundaryLayoutStorage>();
    for (const BoundaryRecord &record : records) {
      storage->patch_faces[record.patch].push_back(record.face);
    }
  } catch (...) {
    materialization_ok = false;
  }
  require_collective(
      mpi, materialization_ok,
      "immersed domain: invariant active_boundary_materialization failed");
  require_collective(
      mpi, patch_ids_ok,
      "immersed domain: invariant active_boundary_patch_id failed");
  require_collective(
      mpi, unique, "immersed domain: invariant active_boundary_unique failed");

  std::optional<std::uint32_t> inlet;
  std::optional<std::uint32_t> outlet;
  bool kinds_ok = true;
  std::string_view kind_invariant =
      "immersed domain: invariant active_boundary_kind failed";
  for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
    const auto kind = boundaries.patch(patch).kind();
    if (kind == boundary::BoundaryKind::velocity_inlet) {
      inlet = patch;
      if (storage->patch_faces[patch].empty()) {
        kinds_ok = false;
        kind_invariant =
            "immersed domain: invariant active_velocity_inlet_nonzero failed";
        break;
      }
    } else if (kind == boundary::BoundaryKind::pressure_outlet) {
      outlet = patch;
      if (storage->patch_faces[patch].empty()) {
        kinds_ok = false;
        kind_invariant =
            "immersed domain: invariant active_pressure_outlet_nonzero failed";
        break;
      }
    }
  }
  require_collective(mpi, kinds_ok, kind_invariant);

  bool periodic_ok = true;
  for (const BoundaryRecord &record : records) {
    const auto &descriptor = boundaries.patch(record.patch);
    if (descriptor.kind() != boundary::BoundaryKind::periodic) {
      if (record.pair.has_value()) {
        periodic_ok = false;
        break;
      }
      continue;
    }
    const auto paired_patch = descriptor.paired_patch_id();
    if (!paired_patch.has_value() || !record.pair.has_value() ||
        *paired_patch >= 6U) {
      periodic_ok = false;
      break;
    }
    const auto paired_record =
        std::lower_bound(records.begin(), records.end(),
                         std::pair<std::uint32_t, mesh::GlobalFaceId>{
                             *paired_patch, *record.pair},
                         [](const BoundaryRecord &candidate, const auto &key) {
                           return std::tie(candidate.patch, candidate.face) <
                                  std::tie(key.first, key.second);
                         });
    if (paired_record == records.end() ||
        paired_record->patch != *paired_patch ||
        paired_record->face != *record.pair ||
        !paired_record->pair.has_value() ||
        *paired_record->pair != record.face) {
      periodic_ok = false;
      break;
    }
  }
  require_collective(
      mpi, periodic_ok,
      "immersed domain: invariant active_periodic_pairing failed");

  storage->open_domain = inlet.has_value() || outlet.has_value();
  storage->has_pressure_reference = outlet.has_value();
  std::uint64_t hash = fnv_offset;
  hash_u64(hash, UINT64_C(0x48554e44424e4433));
  for (const std::uint64_t value : boundary_signature) {
    hash_u64(hash, value);
  }
  hash_u64(hash, records.size());
  for (const BoundaryRecord &record : records) {
    hash_u64(hash, record.patch);
    hash_u64(hash, record.face);
    hash_u64(hash, record.pair.has_value()
                       ? *record.pair
                       : std::numeric_limits<std::uint64_t>::max());
  }
  hash_u64(hash, storage->open_domain ? 1U : 0U);
  hash_u64(hash, storage->has_pressure_reference ? 1U : 0U);
  for (std::uint32_t patch = 0U; patch < 6U; ++patch) {
    hash_u64(hash, patch);
    hash_u64(hash, storage->patch_faces[patch].size());
    for (const mesh::GlobalFaceId face : storage->patch_faces[patch]) {
      hash_u64(hash, face);
    }
  }
  storage->fingerprint = hash;
  return storage;
}

std::uint64_t
coverage_fingerprint(const ImmersedSurface &surface,
                     const std::vector<CanonicalLink> &canonical_links,
                     const std::vector<CellRegion> &global_regions,
                     const detail::PeriodicSurfaceWindow &surface_window,
                     const runtime::MpiContext &mpi) {
  constexpr std::array<std::array<double, 3>, 3> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  constexpr double witness_factor = 2.0 * 1.7320508075688772935;

  std::uint64_t hash = fnv_offset;
  hash_u64(hash, UINT64_C(0x48554e4443565233));
  hash_u64(hash, surface.fingerprint());
  hash_u64(hash, surface.triangle_count());
  std::size_t active_triangle_count = 0U;
  double maximum_distance = 0.0;
  for (TriangleId triangle_id = 0U; triangle_id < surface.triangle_count();
       ++triangle_id) {
    const SurfaceTriangle &triangle = surface.triangle(triangle_id);
    if (!surface_window.triangle_is_active(triangle)) {
      continue;
    }
    ++active_triangle_count;
    for (std::size_t point_index = 0U; point_index < barycentric.size();
         ++point_index) {
      runtime::Real3 point{};
      for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
        point.x +=
            barycentric[point_index][vertex] * triangle.vertices_m[vertex].x;
        point.y +=
            barycentric[point_index][vertex] * triangle.vertices_m[vertex].y;
        point.z +=
            barycentric[point_index][vertex] * triangle.vertices_m[vertex].z;
      }

      const CanonicalLink *selected = nullptr;
      double selected_distance = std::numeric_limits<double>::infinity();
      for (const CanonicalLink &link : canonical_links) {
        const double distance = norm(subtract(point, link.record.intercept));
        if (distance <= 2.0 * link.record.fluid_h &&
            (distance < selected_distance ||
             (distance == selected_distance &&
              (selected == nullptr || link.id < selected->id)))) {
          selected = &link;
          selected_distance = distance;
        }
      }
      require_collective(
          mpi, selected != nullptr,
          "immersed domain: invariant coverage_point_link_distance failed");
      const bool active_row =
          selected->record.fluid_cell < global_regions.size() &&
          global_regions[static_cast<std::size_t>(
              selected->record.fluid_cell)] == CellRegion::fluid;
      require_collective(
          mpi, active_row,
          "immersed domain: invariant coverage_active_row failed");
      const double fluid_distance =
          norm(subtract(point, selected->record.fluid_center));
      const double solid_distance =
          norm(subtract(point, selected->record.solid_center));
      require_collective(
          mpi, fluid_distance <= witness_factor * selected->record.fluid_h,
          "immersed domain: invariant coverage_fluid_witness failed");
      require_collective(
          mpi, solid_distance <= witness_factor * selected->record.solid_h,
          "immersed domain: invariant coverage_solid_witness failed");
      maximum_distance = std::max(maximum_distance, selected_distance);
      hash_u64(hash, triangle_id);
      hash_u64(hash, point_index);
      hash_u64(hash, selected->id);
      hash_u64(hash, selected->record.fluid_cell);
      hash_u64(hash, selected->record.solid_cell);
      hash_u64(hash, double_bits(selected_distance));
      hash_u64(hash, double_bits(fluid_distance));
      hash_u64(hash, double_bits(solid_distance));
    }
  }
  require_collective(
      mpi, active_triangle_count != 0U,
      "immersed domain: invariant active_surface_triangle_nonzero failed");
  hash_u64(hash, double_bits(maximum_distance));
  return hash;
}

} // namespace

ImmersedDomain ImmersedDomain::create(
    const ImmersedSurface &surface, const SurfaceQuery &query,
    config::ImmersedFluidSide fluid_side, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi) {
  if (mpi.comm() == MPI_COMM_NULL) {
    throw runtime::Error("immersed domain: MPI context is invalid");
  }

  const bool side_ok = fluid_side == config::ImmersedFluidSide::outside ||
                       fluid_side == config::ImmersedFluidSide::inside;
  require_collective(mpi, side_ok,
                     "immersed domain: invariant fluid_side_valid failed");
  require_collective(mpi, geometry.compatible(topology),
                     "immersed domain: invariant geometry_compatible failed");
  bool query_surface_ok = true;
  try {
    query_surface_ok =
        query.fingerprint() == SurfaceQuery::create(surface).fingerprint();
  } catch (...) {
    query_surface_ok = false;
  }
  require_collective(
      mpi, query_surface_ok,
      "immersed domain: invariant query_surface_compatible failed");

  bool classification_ok = true;
  std::vector<CellRegion> regions;
  try {
    regions.resize(topology.local_cell_count());
    for (mesh::LocalCellId local = 0U; local < topology.local_cell_count();
         ++local) {
      regions[local] =
          query.classify(geometry.cell_center_m(local), fluid_side);
    }
  } catch (...) {
    classification_ok = false;
  }
  require_collective(mpi, classification_ok,
                     "immersed domain: invariant cell_classification failed");

  std::vector<std::uint64_t> local_classification;
  bool classification_record_ok = true;
  try {
    if (topology.owned_cell_count() >
        std::numeric_limits<std::size_t>::max() / 2U) {
      classification_record_ok = false;
    } else {
      local_classification.reserve(2U * topology.owned_cell_count());
      for (mesh::LocalCellId local = 0U; local < topology.owned_cell_count();
           ++local) {
        local_classification.push_back(topology.global_cell_id(local));
        local_classification.push_back(
            static_cast<std::uint64_t>(regions[local]));
      }
    }
  } catch (...) {
    classification_record_ok = false;
  }
  require_collective(
      mpi, classification_record_ok,
      "immersed domain: invariant classification_record_allocation failed");
  const auto gathered_classification =
      gather_u64_chunked(local_classification, mpi);
  require_collective(
      mpi, gathered_classification.size() % 2U == 0U,
      "immersed domain: invariant classification_record_size failed");
  std::vector<std::pair<mesh::GlobalCellId, CellRegion>> global_classification;
  bool global_classification_allocation_ok = true;
  try {
    global_classification.reserve(gathered_classification.size() / 2U);
    for (std::size_t offset = 0U; offset < gathered_classification.size();
         offset += 2U) {
      const std::uint64_t raw_region = gathered_classification[offset + 1U];
      if (raw_region > static_cast<std::uint64_t>(CellRegion::solid)) {
        classification_record_ok = false;
        break;
      }
      global_classification.push_back({gathered_classification[offset],
                                       static_cast<CellRegion>(raw_region)});
    }
  } catch (...) {
    global_classification_allocation_ok = false;
  }
  require_collective(
      mpi, global_classification_allocation_ok && classification_record_ok,
      "immersed domain: invariant classification_record_decode failed");
  bool exact_global_cells = true;
  std::vector<CellRegion> global_regions;
  std::size_t fluid_count = 0U;
  std::size_t solid_count = 0U;
  bool global_regions_ok = true;
  try {
    std::sort(global_classification.begin(), global_classification.end());
    exact_global_cells =
        global_classification.size() == topology.global_cell_count();
    for (std::size_t index = 0U;
         exact_global_cells && index < global_classification.size(); ++index) {
      exact_global_cells = global_classification[index].first == index;
    }
    global_regions.reserve(global_classification.size());
    for (const auto &[id, region] : global_classification) {
      static_cast<void>(id);
      global_regions.push_back(region);
      fluid_count += region == CellRegion::fluid ? 1U : 0U;
      solid_count += region == CellRegion::solid ? 1U : 0U;
    }
  } catch (...) {
    global_regions_ok = false;
  }
  require_collective(
      mpi, global_regions_ok && exact_global_cells,
      "immersed domain: invariant global_classification_partition failed");
  bool local_classification_agrees = true;
  for (mesh::LocalCellId local = 0U; local < topology.local_cell_count();
       ++local) {
    const mesh::GlobalCellId global = topology.global_cell_id(local);
    if (global >= global_regions.size() ||
        regions[local] != global_regions[static_cast<std::size_t>(global)]) {
      local_classification_agrees = false;
      break;
    }
  }
  require_collective(mpi, local_classification_agrees,
                     "immersed domain: invariant "
                     "local_global_classification_agreement failed");

  require_collective(mpi, fluid_count != 0U,
                     "immersed domain: invariant fluid_cell_nonzero failed");
  require_collective(mpi, solid_count != 0U,
                     "immersed domain: invariant solid_cell_nonzero failed");
  const std::size_t enclosed_count =
      fluid_side == config::ImmersedFluidSide::inside ? fluid_count
                                                      : solid_count;
  require_collective(
      mpi, enclosed_count > 1U,
      "immersed domain: invariant surface_resolved_by_multiple_centres failed");

  bool links_ok = true;
  std::vector<LinkRecord> local_shadow_links;
  try {
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value() || topology.patch_id(face).has_value()) {
        continue;
      }
      const mesh::LocalCellId owner = topology.owner(face);
      if (regions[owner] == regions[*neighbour]) {
        continue;
      }
      const runtime::Real3 owner_center = geometry.cell_center_m(owner);
      const runtime::Real3 neighbour_center =
          geometry.cell_center_m(*neighbour);
      const auto intersections =
          query.segment_intersections(owner_center, neighbour_center);
      if (intersections.size() != 1U) {
        links_ok = false;
        break;
      }
      const SegmentIntersection &intersection = intersections.front();
      const bool owner_fluid = regions[owner] == CellRegion::fluid;
      const mesh::LocalCellId fluid_local = owner_fluid ? owner : *neighbour;
      const mesh::LocalCellId solid_local = owner_fluid ? *neighbour : owner;
      const mesh::GlobalCellId owner_global = topology.global_cell_id(owner);
      const mesh::GlobalCellId neighbour_global =
          topology.global_cell_id(*neighbour);
      const SurfaceTriangle &triangle = surface.triangle(intersection.triangle);
      const double normal_sign =
          fluid_side == config::ImmersedFluidSide::outside ? 1.0 : -1.0;
      LinkRecord record;
      record.face = topology.global_face_id(face);
      record.owned =
          topology.face_ownership(face) == mesh::EntityOwnership::owned;
      record.key = {std::min(owner_global, neighbour_global),
                    std::max(owner_global, neighbour_global),
                    intersection.triangle,
                    {double_bits(intersection.point_m.x),
                     double_bits(intersection.point_m.y),
                     double_bits(intersection.point_m.z)}};
      record.fluid_cell = topology.global_cell_id(fluid_local);
      record.solid_cell = topology.global_cell_id(solid_local);
      record.intercept = intersection.point_m;
      record.normal = multiply(triangle.geometric_outward_normal, normal_sign);
      record.fraction = owner_fluid ? intersection.segment_fraction
                                    : 1.0 - intersection.segment_fraction;
      record.fluid_center = geometry.cell_center_m(fluid_local);
      record.solid_center = geometry.cell_center_m(solid_local);
      record.fluid_h = std::cbrt(geometry.cell_volume_m3(fluid_local));
      record.solid_h = std::cbrt(geometry.cell_volume_m3(solid_local));
      const runtime::Real3 toward_fluid =
          subtract(record.fluid_center, record.intercept);
      if (!finite(record.intercept) || !finite(record.normal) ||
          !std::isfinite(record.fraction) || record.fraction <= 0.0 ||
          record.fraction >= 1.0 || !std::isfinite(record.fluid_h) ||
          record.fluid_h <= 0.0 || !std::isfinite(record.solid_h) ||
          record.solid_h <= 0.0 || dot(record.normal, toward_fluid) <= 0.0) {
        links_ok = false;
        break;
      }
      local_shadow_links.push_back(record);
    }
  } catch (...) {
    links_ok = false;
  }
  require_collective(
      mpi, links_ok,
      "immersed domain: invariant link_single_intersection failed");

  std::vector<std::uint64_t> encoded_links;
  bool link_encode_ok = true;
  try {
    if (local_shadow_links.size() >
        std::numeric_limits<std::size_t>::max() / link_record_width) {
      link_encode_ok = false;
    } else {
      encoded_links.reserve(local_shadow_links.size() * link_record_width);
      for (const LinkRecord &record : local_shadow_links) {
        append_link_record(encoded_links, record);
      }
    }
  } catch (...) {
    link_encode_ok = false;
  }
  require_collective(
      mpi, link_encode_ok,
      "immersed domain: invariant link_record_allocation failed");
  const auto gathered_links = gather_u64_chunked(encoded_links, mpi);
  require_collective(mpi, gathered_links.size() % link_record_width == 0U,
                     "immersed domain: invariant link_record_size failed");

  std::vector<LinkRecord> shadows;
  bool link_decode_ok = true;
  try {
    shadows.reserve(gathered_links.size() / link_record_width);
    for (std::size_t offset = 0U; offset < gathered_links.size();
         offset += link_record_width) {
      shadows.push_back(decode_link_record(gathered_links.data() + offset));
    }
  } catch (...) {
    link_decode_ok = false;
  }
  require_collective(mpi, link_decode_ok,
                     "immersed domain: invariant link_record_decode failed");

  std::vector<LinkRecord> unique_links;
  bool shadow_agreement = true;
  bool link_agreement_materialization_ok = true;
  try {
    std::sort(shadows.begin(), shadows.end(),
              [](const LinkRecord &first, const LinkRecord &second) {
                return std::tie(first.face, first.owned) <
                       std::tie(second.face, second.owned);
              });
    for (std::size_t begin = 0U; begin < shadows.size();) {
      std::size_t end = begin + 1U;
      std::size_t owned_count = shadows[begin].owned ? 1U : 0U;
      while (end < shadows.size() && shadows[end].face == shadows[begin].face) {
        if (!same_record(shadows[begin], shadows[end])) {
          shadow_agreement = false;
        }
        owned_count += shadows[end].owned ? 1U : 0U;
        ++end;
      }
      if (owned_count != 1U) {
        shadow_agreement = false;
      }
      unique_links.push_back(shadows[begin]);
      begin = end;
    }
  } catch (...) {
    link_agreement_materialization_ok = false;
  }
  require_collective(
      mpi, link_agreement_materialization_ok,
      "immersed domain: invariant link_agreement_materialization failed");
  require_collective(
      mpi, shadow_agreement,
      "immersed domain: invariant partition_link_agreement failed");
  require_collective(
      mpi, !unique_links.empty(),
      "immersed domain: invariant interface_link_nonzero failed");

  bool unique_keys = true;
  bool canonical_materialization_ok = true;
  std::vector<CanonicalLink> canonical_links;
  try {
    std::sort(unique_links.begin(), unique_links.end(),
              [](const LinkRecord &first, const LinkRecord &second) {
                return less_key(first.key, second.key);
              });
    for (std::size_t index = 1U; index < unique_links.size(); ++index) {
      if (same_key(unique_links[index - 1U].key, unique_links[index].key)) {
        unique_keys = false;
      }
    }
    canonical_links.reserve(unique_links.size());
    for (std::size_t id = 0U; id < unique_links.size(); ++id) {
      canonical_links.push_back({id, unique_links[id]});
    }
  } catch (...) {
    canonical_materialization_ok = false;
  }
  require_collective(
      mpi, canonical_materialization_ok,
      "immersed domain: invariant canonical_link_materialization failed");
  require_collective(
      mpi, unique_keys,
      "immersed domain: invariant interface_link_key_unique failed");

  std::vector<ImmersedLink> local_links;
  bool local_link_materialization_ok = true;
  bool local_link_lookup_ok = true;
  try {
    for (const LinkRecord &local : local_shadow_links) {
      if (!local.owned) {
        continue;
      }
      const auto found = std::lower_bound(
          canonical_links.begin(), canonical_links.end(), local.key,
          [](const CanonicalLink &candidate, const LinkKey &key) {
            return less_key(candidate.record.key, key);
          });
      if (found == canonical_links.end() ||
          !same_key(found->record.key, local.key)) {
        local_link_lookup_ok = false;
        break;
      }
      local_links.push_back({found->id, local.fluid_cell, local.solid_cell,
                             local.key.triangle, local.intercept, local.normal,
                             local.fraction});
    }
    std::sort(local_links.begin(), local_links.end(),
              [](const ImmersedLink &first, const ImmersedLink &second) {
                return first.id < second.id;
              });
  } catch (...) {
    local_link_materialization_ok = false;
  }
  require_collective(
      mpi, local_link_materialization_ok,
      "immersed domain: invariant local_link_materialization failed");
  require_collective(mpi, local_link_lookup_ok,
                     "immersed domain: invariant local_link_lookup failed");

  std::shared_ptr<const detail::ActiveCellLayoutStorage> active_storage;
  bool active_cell_materialization_ok = true;
  try {
    active_storage = build_active_cells(topology, regions);
  } catch (...) {
    active_cell_materialization_ok = false;
  }
  require_collective(
      mpi, active_cell_materialization_ok,
      "immersed domain: invariant active_cell_materialization failed");
  const ActiveCellLayout active_cells(active_storage);
  const auto boundary_storage =
      build_active_boundaries(topology, boundaries, regions, mpi);
  const ActiveBoundaryLayout active_boundaries(boundary_storage);

  double local_h_max = 0.0;
  bool local_h_ok = true;
  for (mesh::LocalCellId local = 0U; local < topology.owned_cell_count();
       ++local) {
    if (regions[local] != CellRegion::fluid) {
      continue;
    }
    const double h = std::cbrt(geometry.cell_volume_m3(local));
    if (!std::isfinite(h) || h <= 0.0) {
      local_h_ok = false;
      break;
    }
    local_h_max = std::max(local_h_max, h);
  }
  require_collective(mpi, local_h_ok,
                     "immersed domain: invariant active_cell_size failed");
  mpi.allreduce_fp64_in_place(&local_h_max, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  const runtime::Real3 minimum = surface.bounding_box_min_m();
  const runtime::Real3 maximum = surface.bounding_box_max_m();
  const runtime::Real3 origin = geometry.origin_m();
  const runtime::Real3 extent = {origin.x + geometry.length_m().x,
                                 origin.y + geometry.length_m().y,
                                 origin.z + geometry.length_m().z};
  const detail::PeriodicSurfaceWindow surface_window(
      origin, extent, topology.periodicity());
  const bool separated = surface_window.bounding_box_is_admissible(
      minimum, maximum, local_h_max);
  require_collective(
      mpi, separated,
      "immersed domain: invariant surface_domain_separation failed");
  bool periodic_split_ok = true;
  for (TriangleId triangle_id = 0U;
       periodic_split_ok && triangle_id < surface.triangle_count();
       ++triangle_id) {
    periodic_split_ok = surface_window.triangle_is_split_at_periodic_planes(
        surface.triangle(triangle_id));
  }
  require_collective(
      mpi, periodic_split_ok,
      "immersed domain: invariant periodic_surface_plane_split failed");

  const std::uint64_t classification = compute_classification_fingerprint(
      surface, query, fluid_side, global_regions);
  const std::uint64_t coverage = coverage_fingerprint(
      surface, canonical_links, global_regions, surface_window, mpi);
  std::shared_ptr<detail::ImmersedDomainStorage> storage;
  bool publication_ok = true;
  try {
    storage = std::make_shared<detail::ImmersedDomainStorage>(
        std::move(regions), std::move(local_links), active_cells,
        active_boundaries, classification, coverage,
        static_cast<std::uint64_t>(topology.local_cell_count()));
  } catch (...) {
    publication_ok = false;
  }
  require_collective(
      mpi, publication_ok,
      "immersed domain: invariant publication_materialization failed");
  return ImmersedDomain(std::move(storage));
}

std::size_t ActiveCellLayout::owned_active_count() const noexcept {
  return storage_->owned_active_count;
}

std::size_t ActiveCellLayout::local_active_count() const noexcept {
  return storage_->ordered_global_ids.size();
}

bool ActiveCellLayout::active(mesh::LocalCellId local_cell) const {
  if (local_cell >= storage_->local_cell_count) {
    throw runtime::Error(
        "immersed active layout: local cell index is outside layout");
  }
  return storage_->active_indices[local_cell] != no_active_index;
}

std::optional<std::size_t>
ActiveCellLayout::active_index(mesh::LocalCellId local_cell) const {
  if (local_cell >= storage_->local_cell_count) {
    throw runtime::Error(
        "immersed active layout: local cell index is outside layout");
  }
  const std::size_t index = storage_->active_indices[local_cell];
  return index == no_active_index ? std::nullopt
                                  : std::optional<std::size_t>{index};
}

const std::vector<mesh::GlobalCellId> &
ActiveCellLayout::ordered_global_ids() const noexcept {
  return storage_->ordered_global_ids;
}

std::uint64_t ActiveCellLayout::fingerprint() const noexcept {
  return storage_->fingerprint;
}

const std::vector<mesh::GlobalFaceId> &
ActiveBoundaryLayout::patch_faces(std::uint32_t stable_patch_id) const {
  if (stable_patch_id >= storage_->patch_faces.size()) {
    throw runtime::Error(
        "active boundary layout: stable patch ID is outside layout");
  }
  return storage_->patch_faces[stable_patch_id];
}

bool ActiveBoundaryLayout::open_domain() const noexcept {
  return storage_->open_domain;
}

bool ActiveBoundaryLayout::has_pressure_reference() const noexcept {
  return storage_->has_pressure_reference;
}

std::uint64_t ActiveBoundaryLayout::fingerprint() const noexcept {
  return storage_->fingerprint;
}

CellRegion ImmersedDomain::region(mesh::LocalCellId local_cell) const {
  if (local_cell >= storage_->regions.size()) {
    throw runtime::Error(
        "immersed domain: local cell index is outside classification");
  }
  return storage_->regions[local_cell];
}

const std::vector<ImmersedLink> &ImmersedDomain::links() const noexcept {
  return storage_->links;
}

const ActiveCellLayout &ImmersedDomain::active_cells() const noexcept {
  return storage_->active_cells;
}

const ActiveBoundaryLayout &ImmersedDomain::active_boundaries() const noexcept {
  return storage_->active_boundaries;
}

std::uint64_t ImmersedDomain::classification_fingerprint() const noexcept {
  return storage_->classification_fingerprint;
}

std::uint64_t ImmersedDomain::surface_coverage_fingerprint() const noexcept {
  return storage_->surface_coverage_fingerprint;
}

diagnostics::Stage3PerformanceCounters
ImmersedDomain::performance_counters() const noexcept {
  return storage_->performance;
}

} // namespace hundun::immersed
