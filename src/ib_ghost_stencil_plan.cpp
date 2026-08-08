// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_ghost_stencil_plan.hpp"

#include "ib_deterministic_qr_detail.hpp"
#include "ib_ghost_stencil_plan_detail.hpp"
#include "ib_quadratic_reconstruction_detail.hpp"

#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace hundun::immersed {
namespace {

constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);
constexpr std::size_t kCellRecordWidth = 12U;
constexpr std::size_t kLinkRecordWidth = 14U;
constexpr std::size_t kMinimumDonors = 14U;
constexpr std::size_t kMaximumDonors = 32U;
constexpr std::size_t kParallelWorkerBudget = 96U;
constexpr int kMaximumReach = 4;

struct CellRecord final {
  mesh::GlobalCellId id{};
  runtime::Int3 logical{};
  runtime::Real3 center_m{};
  double volume_m3{};
  CellRegion region{CellRegion::solid};
  int owner_rank{};
};

struct GlobalLink final {
  ImmersedLink record{};
};

struct Frame final {
  runtime::Real3 normal{};
  runtime::Real3 tangent1{};
  runtime::Real3 tangent2{};
};

struct CandidateReconstruction final {
  QuadraticReconstruction reconstruction;
  std::vector<runtime::Int3> donor_cells;
};

struct GhostFunctionalTarget final {
  runtime::Real3 point_m{};
  runtime::Int3 logical{};
};

enum class ReconstructionSelection {
  widest_accepted_nearest_prefix,
  widest_accepted_directional_prefix,
  minimum_wall_value_amplification,
  minimum_wall_normal_gradient_amplification,
  minimum_ghost_pressure_neumann_amplification,
  minimum_directional_ghost_pressure_neumann_amplification,
};

class ReconstructionSelectionError final : public runtime::Error {
public:
  using runtime::Error::Error;
};

double weight_l1(const std::vector<WeightedDonor> &weights) noexcept {
  double result = 0.0;
  for (const auto &weight : weights)
    result += std::abs(weight.weight);
  return result;
}

double difference_weight_l1(const std::vector<WeightedDonor> &first,
                            const std::vector<WeightedDonor> &second) {
  if (first.size() != second.size())
    throw runtime::Error("ghost plan functional donor layouts disagree");
  double result = 0.0;
  for (std::size_t index = 0U; index < first.size(); ++index) {
    if (first[index].global_cell != second[index].global_cell)
      throw runtime::Error("ghost plan functional donor IDs disagree");
    const double value = first[index].weight - second[index].weight;
    if (!std::isfinite(value))
      throw runtime::Error("ghost plan functional weight is non-finite");
    result += std::abs(value);
  }
  return result;
}

double combined_weight_l1(const std::vector<WeightedDonor> &first,
                          double first_scale,
                          const std::vector<WeightedDonor> &second,
                          double second_scale) {
  if (first.size() != second.size())
    throw runtime::Error("ghost plan functional donor layouts disagree");
  double result = 0.0;
  for (std::size_t index = 0U; index < first.size(); ++index) {
    if (first[index].global_cell != second[index].global_cell)
      throw runtime::Error("ghost plan functional donor IDs disagree");
    const double value =
        first_scale * first[index].weight + second_scale * second[index].weight;
    if (!std::isfinite(value))
      throw runtime::Error("ghost plan functional weight is non-finite");
    result += std::abs(value);
  }
  return result;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

runtime::Real3 add(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

runtime::Real3 subtract(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

runtime::Real3 multiply(runtime::Real3 value, double factor) noexcept {
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

runtime::Real3 cross(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double norm(runtime::Real3 value) noexcept {
  return std::hypot(value.x, value.y, value.z);
}

double squared_distance(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  const auto difference = subtract(lhs, rhs);
  return dot(difference, difference);
}

std::uint32_t logical_reach(runtime::Int3 lhs, runtime::Int3 rhs) noexcept {
  const auto distance = [](int first, int second) {
    const auto difference =
        static_cast<std::int64_t>(first) - static_cast<std::int64_t>(second);
    return static_cast<std::uint32_t>(difference < 0 ? -difference
                                                     : difference);
  };
  return std::max(
      {distance(lhs.x, rhs.x), distance(lhs.y, rhs.y), distance(lhs.z, rhs.z)});
}

std::uint32_t
donor_cloud_reach(const std::vector<runtime::Int3> &donors) noexcept {
  std::uint32_t result = 0U;
  for (const auto first : donors)
    for (const auto second : donors)
      result = std::max(result, logical_reach(first, second));
  return result;
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

std::uint64_t signed_bits(int value) noexcept {
  const auto wide = static_cast<std::int64_t>(value);
  std::uint64_t result = 0U;
  static_assert(sizeof(result) == sizeof(wide));
  std::memcpy(&result, &wide, sizeof(result));
  return result;
}

int int_from_signed_bits(std::uint64_t value) {
  std::int64_t wide = 0;
  static_assert(sizeof(value) == sizeof(wide));
  std::memcpy(&wide, &value, sizeof(wide));
  if (wide < std::numeric_limits<int>::min() ||
      wide > std::numeric_limits<int>::max()) {
    throw runtime::Error("ghost plan signed integer record is out of range");
  }
  return static_cast<int>(wide);
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    hash ^= (value >> (8U * byte)) & UINT64_C(0xff);
    hash *= kFnvPrime;
  }
}

void hash_real3(std::uint64_t &hash, runtime::Real3 value) noexcept {
  hash_u64(hash, double_bits(value.x));
  hash_u64(hash, double_bits(value.y));
  hash_u64(hash, double_bits(value.z));
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
  const auto status = runtime::collective_status(mpi, local_ok, invariant);
  if (!status.ok) {
    throw runtime::Error(status.message + " (lowest failing rank " +
                         std::to_string(status.failing_rank) + ")");
  }
}

std::vector<std::uint64_t>
gather_u64_chunked(const std::vector<std::uint64_t> &local,
                   const runtime::MpiContext &mpi) {
  const bool local_size_ok =
      local.size() <= std::numeric_limits<std::uint64_t>::max();
  require_collective(mpi, local_size_ok,
                     "ghost plan: local record size is not representable");

  std::vector<std::uint64_t> sizes;
  bool sizes_ok = true;
  try {
    sizes.resize(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    sizes_ok = false;
  }
  require_collective(mpi, sizes_ok,
                     "ghost plan: collective size allocation failed");
  for (int root = 0; root < mpi.size(); ++root) {
    std::uint64_t size =
        mpi.rank() == root ? static_cast<std::uint64_t>(local.size()) : 0U;
    check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, root, mpi.comm()),
              "MPI_Bcast ghost-plan record count");
    sizes[static_cast<std::size_t>(root)] = size;
  }

  std::size_t total = 0U;
  bool representable = true;
  for (const auto size : sizes) {
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
  require_collective(mpi, representable,
                     "ghost plan: collective record size is not representable");

  std::vector<std::uint64_t> gathered;
  bool allocation_ok = true;
  try {
    gathered.resize(total);
  } catch (...) {
    allocation_ok = false;
  }
  require_collective(mpi, allocation_ok,
                     "ghost plan: collective record allocation failed");

  constexpr std::size_t kMaximumChunk =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t offset = 0U;
  for (int root = 0; root < mpi.size(); ++root) {
    const auto root_size =
        static_cast<std::size_t>(sizes[static_cast<std::size_t>(root)]);
    if (mpi.rank() == root) {
      std::copy(local.begin(), local.end(),
                gathered.begin() + static_cast<std::ptrdiff_t>(offset));
    }
    std::size_t sent = 0U;
    while (sent < root_size) {
      const std::size_t chunk = std::min(kMaximumChunk, root_size - sent);
      check_mpi(MPI_Bcast(gathered.data() + offset + sent,
                          static_cast<int>(chunk), MPI_UINT64_T, root,
                          mpi.comm()),
                "MPI_Bcast ghost-plan records");
      sent += chunk;
    }
    offset += root_size;
  }
  return gathered;
}

void validate_common_inputs(const ImmersedSurface &surface,
                            const SurfaceQuery &query,
                            const ImmersedDomain &domain,
                            const mesh::MeshTopology &topology,
                            const mesh::MeshGeometry &geometry,
                            const runtime::MpiContext &mpi) {
  const bool local_ok =
      mpi.comm() != MPI_COMM_NULL && surface.triangle_count() > 0U &&
      surface.fingerprint() != 0U && query.fingerprint() != 0U &&
      domain.classification_fingerprint() != 0U &&
      domain.surface_coverage_fingerprint() != 0U &&
      geometry.compatible(topology);
  require_collective(mpi, local_ok,
                     "ghost plan: common input validation failed");

  const std::array<std::uint64_t, 6> local_signature{
      surface.fingerprint(),
      query.fingerprint(),
      domain.classification_fingerprint(),
      domain.surface_coverage_fingerprint(),
      static_cast<std::uint64_t>(topology.global_cell_count()),
      static_cast<std::uint64_t>(surface.triangle_count())};
  auto reference = local_signature;
  check_mpi(MPI_Bcast(reference.data(), static_cast<int>(reference.size()),
                      MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast ghost-plan signature");
  require_collective(mpi, local_signature == reference,
                     "ghost plan: rank input signature disagrees");
}

void validate_decomposition(
    const mesh::MeshTopology &topology,
    const runtime::StructuredDecomposition &decomposition,
    const runtime::MpiContext &mpi) {
  int relation = MPI_UNEQUAL;
  check_mpi(MPI_Comm_compare(decomposition.comm(), mpi.comm(), &relation),
            "MPI_Comm_compare ghost-plan decomposition");
  const bool communicator_ok =
      relation == MPI_IDENT || relation == MPI_CONGRUENT;
  const auto topology_box = topology.owned_global_box();
  const auto decomposition_box = decomposition.owned_box();
  const bool layout_ok =
      topology.global_extent().x == decomposition.global_extent().x &&
      topology.global_extent().y == decomposition.global_extent().y &&
      topology.global_extent().z == decomposition.global_extent().z &&
      topology_box.begin.x == decomposition_box.begin.x &&
      topology_box.begin.y == decomposition_box.begin.y &&
      topology_box.begin.z == decomposition_box.begin.z &&
      topology_box.end.x == decomposition_box.end.x &&
      topology_box.end.y == decomposition_box.end.y &&
      topology_box.end.z == decomposition_box.end.z;
  require_collective(mpi, communicator_ok && layout_ok,
                     "ghost plan: decomposition is incompatible");
}

std::vector<CellRecord> build_global_cells(const ImmersedDomain &domain,
                                           const mesh::MeshTopology &topology,
                                           const mesh::MeshGeometry &geometry,
                                           const runtime::MpiContext &mpi) {
  std::vector<std::uint64_t> encoded;
  bool local_ok = true;
  try {
    if (topology.owned_cell_count() >
        std::numeric_limits<std::size_t>::max() / kCellRecordWidth) {
      local_ok = false;
    } else {
      encoded.reserve(topology.owned_cell_count() * kCellRecordWidth);
      for (mesh::LocalCellId local = 0U; local < topology.owned_cell_count();
           ++local) {
        const auto logical = topology.global_cell(local);
        const auto center = geometry.cell_center_m(local);
        const double volume = geometry.cell_volume_m3(local);
        const auto region = domain.region(local);
        if (!finite(center) || !std::isfinite(volume) || volume <= 0.0) {
          local_ok = false;
          break;
        }
        encoded.insert(encoded.end(),
                       {topology.global_cell_id(local), signed_bits(logical.x),
                        signed_bits(logical.y), signed_bits(logical.z),
                        double_bits(center.x), double_bits(center.y),
                        double_bits(center.z), double_bits(volume),
                        static_cast<std::uint64_t>(region),
                        signed_bits(mpi.rank()), 0U, 0U});
      }
    }
  } catch (...) {
    local_ok = false;
  }
  require_collective(mpi, local_ok,
                     "ghost plan: local cell catalog construction failed");
  const auto gathered = gather_u64_chunked(encoded, mpi);
  require_collective(mpi, gathered.size() % kCellRecordWidth == 0U,
                     "ghost plan: cell catalog record width is invalid");

  std::vector<CellRecord> records;
  bool decode_ok = true;
  try {
    records.reserve(gathered.size() / kCellRecordWidth);
    for (std::size_t offset = 0U; offset < gathered.size();
         offset += kCellRecordWidth) {
      const auto raw_region = gathered[offset + 8U];
      if (raw_region > static_cast<std::uint64_t>(CellRegion::solid)) {
        decode_ok = false;
        break;
      }
      records.push_back({gathered[offset],
                         {int_from_signed_bits(gathered[offset + 1U]),
                          int_from_signed_bits(gathered[offset + 2U]),
                          int_from_signed_bits(gathered[offset + 3U])},
                         {double_from_bits(gathered[offset + 4U]),
                          double_from_bits(gathered[offset + 5U]),
                          double_from_bits(gathered[offset + 6U])},
                         double_from_bits(gathered[offset + 7U]),
                         static_cast<CellRegion>(raw_region),
                         int_from_signed_bits(gathered[offset + 9U])});
    }
  } catch (...) {
    decode_ok = false;
  }
  require_collective(mpi, decode_ok, "ghost plan: cell catalog decode failed");
  std::sort(records.begin(), records.end(),
            [](const auto &lhs, const auto &rhs) { return lhs.id < rhs.id; });
  bool complete = records.size() == topology.global_cell_count();
  const auto extent = topology.global_extent();
  for (std::size_t index = 0U; complete && index < records.size(); ++index) {
    const auto &record = records[index];
    complete = record.id == index && record.logical.x >= 0 &&
               record.logical.x < extent.x && record.logical.y >= 0 &&
               record.logical.y < extent.y && record.logical.z >= 0 &&
               record.logical.z < extent.z && finite(record.center_m) &&
               std::isfinite(record.volume_m3) && record.volume_m3 > 0.0 &&
               record.owner_rank >= 0 && record.owner_rank < mpi.size() &&
               topology.global_cell_id(record.logical) == record.id;
  }
  require_collective(mpi, complete,
                     "ghost plan: global cell catalog is incomplete");
  return records;
}

void append_link(std::vector<std::uint64_t> &encoded,
                 const ImmersedLink &link) {
  encoded.insert(encoded.end(),
                 {link.id, link.fluid_cell, link.solid_cell, link.triangle,
                  double_bits(link.wall_intercept_m.x),
                  double_bits(link.wall_intercept_m.y),
                  double_bits(link.wall_intercept_m.z),
                  double_bits(link.solid_to_fluid_normal.x),
                  double_bits(link.solid_to_fluid_normal.y),
                  double_bits(link.solid_to_fluid_normal.z),
                  double_bits(link.fluid_to_wall_fraction), 0U, 0U, 0U});
}

ImmersedLink decode_link(const std::uint64_t *record) {
  return {record[0],
          record[1],
          record[2],
          record[3],
          {double_from_bits(record[4]), double_from_bits(record[5]),
           double_from_bits(record[6])},
          {double_from_bits(record[7]), double_from_bits(record[8]),
           double_from_bits(record[9])},
          double_from_bits(record[10])};
}

bool same_link(const ImmersedLink &lhs, const ImmersedLink &rhs) noexcept {
  return lhs.id == rhs.id && lhs.fluid_cell == rhs.fluid_cell &&
         lhs.solid_cell == rhs.solid_cell && lhs.triangle == rhs.triangle &&
         double_bits(lhs.wall_intercept_m.x) ==
             double_bits(rhs.wall_intercept_m.x) &&
         double_bits(lhs.wall_intercept_m.y) ==
             double_bits(rhs.wall_intercept_m.y) &&
         double_bits(lhs.wall_intercept_m.z) ==
             double_bits(rhs.wall_intercept_m.z) &&
         double_bits(lhs.solid_to_fluid_normal.x) ==
             double_bits(rhs.solid_to_fluid_normal.x) &&
         double_bits(lhs.solid_to_fluid_normal.y) ==
             double_bits(rhs.solid_to_fluid_normal.y) &&
         double_bits(lhs.solid_to_fluid_normal.z) ==
             double_bits(rhs.solid_to_fluid_normal.z) &&
         double_bits(lhs.fluid_to_wall_fraction) ==
             double_bits(rhs.fluid_to_wall_fraction);
}

std::vector<GlobalLink> build_global_links(const ImmersedDomain &domain,
                                           const std::vector<CellRecord> &cells,
                                           const ImmersedSurface &surface,
                                           const runtime::MpiContext &mpi) {
  std::vector<std::uint64_t> encoded;
  bool local_ok = true;
  try {
    if (domain.links().size() >
        std::numeric_limits<std::size_t>::max() / kLinkRecordWidth) {
      local_ok = false;
    } else {
      encoded.reserve(domain.links().size() * kLinkRecordWidth);
      for (const auto &link : domain.links()) {
        append_link(encoded, link);
      }
    }
  } catch (...) {
    local_ok = false;
  }
  require_collective(mpi, local_ok,
                     "ghost plan: local link catalog construction failed");
  const auto gathered = gather_u64_chunked(encoded, mpi);
  require_collective(mpi, gathered.size() % kLinkRecordWidth == 0U,
                     "ghost plan: link catalog record width is invalid");
  std::vector<ImmersedLink> decoded;
  bool decode_ok = true;
  try {
    decoded.reserve(gathered.size() / kLinkRecordWidth);
    for (std::size_t offset = 0U; offset < gathered.size();
         offset += kLinkRecordWidth) {
      decoded.push_back(decode_link(gathered.data() + offset));
    }
  } catch (...) {
    decode_ok = false;
  }
  require_collective(mpi, decode_ok, "ghost plan: link catalog decode failed");
  std::sort(decoded.begin(), decoded.end(),
            [](const auto &lhs, const auto &rhs) { return lhs.id < rhs.id; });

  std::vector<GlobalLink> links;
  bool canonical_ok = true;
  try {
    for (std::size_t index = 0U; index < decoded.size();) {
      const auto canonical = decoded[index];
      std::size_t next = index + 1U;
      while (next < decoded.size() && decoded[next].id == canonical.id) {
        if (!same_link(canonical, decoded[next])) {
          canonical_ok = false;
          break;
        }
        ++next;
      }
      if (!canonical_ok) {
        break;
      }
      links.push_back({canonical});
      index = next;
    }
  } catch (...) {
    canonical_ok = false;
  }
  if (links.empty()) {
    canonical_ok = false;
  }
  for (std::size_t index = 0U; canonical_ok && index < links.size(); ++index) {
    const auto &link = links[index].record;
    canonical_ok =
        link.id == index && link.fluid_cell < cells.size() &&
        link.solid_cell < cells.size() &&
        link.triangle < surface.triangle_count() &&
        cells[static_cast<std::size_t>(link.fluid_cell)].region ==
            CellRegion::fluid &&
        cells[static_cast<std::size_t>(link.solid_cell)].region ==
            CellRegion::solid &&
        finite(link.wall_intercept_m) && finite(link.solid_to_fluid_normal) &&
        std::isfinite(link.fluid_to_wall_fraction) &&
        link.fluid_to_wall_fraction > 0.0 && link.fluid_to_wall_fraction < 1.0;
  }
  require_collective(mpi, canonical_ok,
                     "ghost plan: global link catalog is inconsistent");
  return links;
}

// Per-link true body-surface measure vector and patch centroid, aggregated
// with the same triangle-point-to-link association the WallQuadraturePlan
// uses (nearest link within two local cells, deterministic tie-break). This
// is the surface partition the force-authority operator row integrates; the
// sum over the closed body of the measure vectors is zero.
struct SurfaceMeasureByLink final {
  std::vector<runtime::Real3> measures;
  std::vector<runtime::Real3> centroids;
};

SurfaceMeasureByLink build_surface_measure_by_link(
    const ImmersedSurface &surface, const std::vector<GlobalLink> &links,
    const std::vector<CellRecord> &cells) {
  using SpatialBin = std::array<std::int64_t, 3>;
  double maximum_link_h = 0.0;
  for (const auto &global_link : links) {
    const auto fluid_cell = global_link.record.fluid_cell;
    if (fluid_cell >= cells.size())
      throw runtime::Error(
          "ghost plan surface link fluid cell is out of range");
    const double h =
        std::cbrt(cells[static_cast<std::size_t>(fluid_cell)].volume_m3);
    if (!(h > 0.0) || !std::isfinite(h))
      throw runtime::Error("ghost plan surface link scale is invalid");
    maximum_link_h = std::max(maximum_link_h, h);
  }
  if (!(maximum_link_h > 0.0) || !std::isfinite(maximum_link_h))
    throw runtime::Error("ghost plan surface maximum link scale is invalid");
  const double bin_width = 2.0 * maximum_link_h;
  const auto bin_coordinate = [bin_width](double value) {
    const double coordinate = std::floor(value / bin_width);
    if (!std::isfinite(coordinate) ||
        coordinate <
            static_cast<double>(std::numeric_limits<std::int64_t>::min()) ||
        coordinate >
            static_cast<double>(std::numeric_limits<std::int64_t>::max()))
      throw runtime::Error(
          "ghost plan surface spatial-bin coordinate is invalid");
    return static_cast<std::int64_t>(coordinate);
  };
  std::map<SpatialBin, std::vector<std::size_t>> link_bins;
  for (std::size_t link_index = 0U; link_index < links.size(); ++link_index) {
    const auto &link = links[link_index].record;
    const double h =
        std::cbrt(cells[static_cast<std::size_t>(link.fluid_cell)].volume_m3);
    const double radius = 2.0 * h;
    const SpatialBin first{bin_coordinate(link.wall_intercept_m.x - radius),
                           bin_coordinate(link.wall_intercept_m.y - radius),
                           bin_coordinate(link.wall_intercept_m.z - radius)};
    const SpatialBin last{bin_coordinate(link.wall_intercept_m.x + radius),
                          bin_coordinate(link.wall_intercept_m.y + radius),
                          bin_coordinate(link.wall_intercept_m.z + radius)};
    std::array<int, 3> span{};
    for (std::size_t axis = 0U; axis < span.size(); ++axis) {
      const auto difference = last[axis] - first[axis];
      if (difference < 0 || difference > 3)
        throw runtime::Error(
            "ghost plan surface spatial-bin span is invalid");
      span[axis] = static_cast<int>(difference);
    }
    for (int k = 0; k <= span[2]; ++k)
      for (int j = 0; j <= span[1]; ++j)
        for (int i = 0; i <= span[0]; ++i)
          link_bins[{first[0] + i, first[1] + j, first[2] + k}].push_back(
              link_index);
  }
  struct SurfacePatch final {
    runtime::Real3 measure{0.0, 0.0, 0.0};
    runtime::Real3 weighted_point{0.0, 0.0, 0.0};
    double weight{};
  };
  std::vector<SurfacePatch> surface_patches(links.size());
  constexpr std::array<std::array<double, 3>, 3> kBarycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  struct CanonicalTriangle final {
    std::uint64_t id{};
    std::array<std::uint64_t, 13> key{};
  };
  std::vector<CanonicalTriangle> canonical_triangles;
  canonical_triangles.reserve(
      static_cast<std::size_t>(surface.triangle_count()));
  for (std::uint64_t triangle_id = 0U;
       triangle_id < surface.triangle_count(); ++triangle_id) {
    const auto &triangle = surface.triangle(triangle_id);
    CanonicalTriangle canonical;
    canonical.id = triangle_id;
    std::size_t key_index = 0U;
    for (const auto vertex : triangle.vertices_m) {
      canonical.key[key_index++] = double_bits(vertex.x);
      canonical.key[key_index++] = double_bits(vertex.y);
      canonical.key[key_index++] = double_bits(vertex.z);
    }
    canonical.key[key_index++] =
        double_bits(triangle.geometric_outward_normal.x);
    canonical.key[key_index++] =
        double_bits(triangle.geometric_outward_normal.y);
    canonical.key[key_index++] =
        double_bits(triangle.geometric_outward_normal.z);
    canonical.key[key_index] = double_bits(triangle.area_m2);
    canonical_triangles.push_back(canonical);
  }
  std::sort(canonical_triangles.begin(), canonical_triangles.end(),
            [](const CanonicalTriangle &first,
               const CanonicalTriangle &second) {
              return first.key < second.key;
            });
  for (const auto &canonical_triangle : canonical_triangles) {
    const auto triangle_id = canonical_triangle.id;
    const auto &triangle = surface.triangle(triangle_id);
    if (!finite(triangle.geometric_outward_normal) ||
        !std::isfinite(triangle.area_m2) || triangle.area_m2 <= 0.0)
      throw runtime::Error("ghost plan surface triangle geometry is invalid");
    const double point_weight = triangle.area_m2 / 3.0;
    for (std::size_t point_index = 0U; point_index < 3U; ++point_index) {
      runtime::Real3 point{};
      for (std::size_t vertex = 0U; vertex < 3U; ++vertex)
        point = add(point, multiply(triangle.vertices_m[vertex],
                                    kBarycentric[point_index][vertex]));
      const ImmersedLink *selected = nullptr;
      std::size_t selected_index = 0U;
      double selected_distance = std::numeric_limits<double>::infinity();
      const SpatialBin point_bin{bin_coordinate(point.x),
                                 bin_coordinate(point.y),
                                 bin_coordinate(point.z)};
      const auto candidates = link_bins.find(point_bin);
      if (candidates == link_bins.end())
        throw runtime::Error(
            "ghost plan surface point has no spatial-link candidates");
      for (const auto link_index : candidates->second) {
        const auto &link = links[link_index].record;
        const double distance =
            std::sqrt(squared_distance(point, link.wall_intercept_m));
        const double h =
            std::cbrt(cells[static_cast<std::size_t>(link.fluid_cell)]
                          .volume_m3);
        const double normal_alignment = dot(triangle.geometric_outward_normal,
                                            link.solid_to_fluid_normal);
        const double alignment_tolerance =
            512.0 * std::numeric_limits<double>::epsilon();
        if (std::isfinite(distance) && std::isfinite(h) && h > 0.0 &&
            std::isfinite(normal_alignment) &&
            std::abs(normal_alignment) > alignment_tolerance &&
            distance <= 2.0 * h &&
            (distance < selected_distance ||
             (distance == selected_distance &&
              (selected == nullptr || link.id < selected->id)))) {
          selected = &link;
          selected_index = link_index;
          selected_distance = distance;
        }
      }
      if (selected == nullptr)
        throw runtime::Error(
            "ghost plan surface point has no associated interface link");
      const double alignment = dot(triangle.geometric_outward_normal,
                                   selected->solid_to_fluid_normal);
      if (!std::isfinite(alignment) || alignment == 0.0)
        throw runtime::Error(
            "ghost plan surface normal orientation is ambiguous");
      const auto normal = multiply(triangle.geometric_outward_normal,
                                   alignment > 0.0 ? 1.0 : -1.0);
      if (dot(normal, selected->solid_to_fluid_normal) <= 0.0)
        throw runtime::Error(
            "ghost plan surface associated normal is inconsistent");
      auto &patch = surface_patches[selected_index];
      patch.measure = add(patch.measure, multiply(normal, point_weight));
      patch.weighted_point =
          add(patch.weighted_point, multiply(point, point_weight));
      patch.weight += point_weight;
    }
  }
  SurfaceMeasureByLink result;
  result.measures.resize(links.size());
  result.centroids.resize(links.size());
  for (std::size_t link_index = 0U; link_index < links.size(); ++link_index) {
    const auto &patch = surface_patches[link_index];
    if (!(patch.weight > 0.0) || !std::isfinite(patch.weight)) {
      const auto &link = links[link_index].record;
      throw runtime::Error(
          "ghost plan surface patch is empty for link " +
          std::to_string(link.id) + " at intercept (" +
          std::to_string(link.wall_intercept_m.x) + ", " +
          std::to_string(link.wall_intercept_m.y) + ", " +
          std::to_string(link.wall_intercept_m.z) + ")");
    }
    result.measures[link_index] = patch.measure;
    result.centroids[link_index] =
        multiply(patch.weighted_point, 1.0 / patch.weight);
    if (!finite(result.measures[link_index]) ||
        !finite(result.centroids[link_index]))
      throw runtime::Error("ghost plan surface patch result is non-finite");
  }
  return result;
}

Frame deterministic_frame(runtime::Real3 normal) {
  const double length = norm(normal);
  if (!finite(normal) || !std::isfinite(length) || length <= 0.0) {
    throw runtime::Error("ghost plan normal must be finite and nonzero");
  }
  normal = multiply(normal, 1.0 / length);
  const std::array<double, 3> magnitude{std::abs(normal.x), std::abs(normal.y),
                                        std::abs(normal.z)};
  std::size_t selected = 0U;
  for (std::size_t axis = 1U; axis < magnitude.size(); ++axis) {
    if (magnitude[axis] < magnitude[selected]) {
      selected = axis;
    }
  }
  const std::array<runtime::Real3, 3> axes{runtime::Real3{1.0, 0.0, 0.0},
                                           runtime::Real3{0.0, 1.0, 0.0},
                                           runtime::Real3{0.0, 0.0, 1.0}};
  auto tangent1 = cross(axes[selected], normal);
  const double tangent_length = norm(tangent1);
  if (!std::isfinite(tangent_length) || tangent_length <= 0.0) {
    throw runtime::Error("ghost plan tangent construction failed");
  }
  tangent1 = multiply(tangent1, 1.0 / tangent_length);
  const auto tangent2 = cross(normal, tangent1);
  if (!finite(tangent2)) {
    throw runtime::Error("ghost plan tangent frame is non-finite");
  }
  return {normal, tangent1, tangent2};
}

CandidateReconstruction build_reconstruction(
    runtime::Real3 point_m, runtime::Real3 normal, mesh::GlobalCellId anchor_id,
    const std::vector<CellRecord> &cells, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry, ReconstructionSelection selection,
    std::optional<GhostFunctionalTarget> ghost_target = std::nullopt) {
  if (anchor_id >= cells.size()) {
    throw runtime::Error("ghost plan reconstruction anchor is missing");
  }
  const auto &anchor = cells[static_cast<std::size_t>(anchor_id)];
  if (anchor.region != CellRegion::fluid) {
    throw runtime::Error("ghost plan reconstruction anchor is not fluid");
  }
  const Frame frame = deterministic_frame(normal);
  const double scale = std::cbrt(anchor.volume_m3);
  if (!std::isfinite(scale) || scale <= 0.0) {
    throw runtime::Error("ghost plan reconstruction scale is invalid");
  }
  struct Candidate final {
    double distance_squared{};
    mesh::GlobalCellId id{};
    runtime::Int3 logical{};
  };
  std::vector<Candidate> candidates;
  candidates.reserve(729U);
  const auto extent = topology.global_extent();
  const runtime::Int3 first{std::max(0, anchor.logical.x - kMaximumReach),
                            std::max(0, anchor.logical.y - kMaximumReach),
                            std::max(0, anchor.logical.z - kMaximumReach)};
  const runtime::Int3 last{
      std::min(extent.x - 1, anchor.logical.x + kMaximumReach),
      std::min(extent.y - 1, anchor.logical.y + kMaximumReach),
      std::min(extent.z - 1, anchor.logical.z + kMaximumReach)};
  for (int k = first.z; k <= last.z; ++k)
    for (int j = first.y; j <= last.y; ++j)
      for (int i = first.x; i <= last.x; ++i) {
        const runtime::Int3 logical{i, j, k};
        const auto id = topology.global_cell_id(logical);
        if (id >= cells.size())
          throw runtime::Error(
              "ghost plan donor catalog lookup is out of range");
        const auto &cell = cells[static_cast<std::size_t>(id)];
        if (cell.id != id || cell.region != CellRegion::fluid)
          continue;
        const auto offset = subtract(cell.center_m, point_m);
        const double normal_coordinate = dot(offset, frame.normal);
        const double distance = squared_distance(cell.center_m, point_m);
        if (!std::isfinite(normal_coordinate) || normal_coordinate <= 0.0 ||
            !std::isfinite(distance))
          continue;
        candidates.push_back({distance, cell.id, cell.logical});
      }
  std::sort(candidates.begin(), candidates.end(),
            [](const auto &lhs, const auto &rhs) {
              return std::tie(lhs.distance_squared, lhs.id) <
                     std::tie(rhs.distance_squared, rhs.id);
            });
  if (candidates.size() < kMinimumDonors) {
    throw runtime::Error("ghost plan has fewer than fourteen fluid donors");
  }
  const std::size_t maximum = std::min(kMaximumDonors, candidates.size());
  std::vector<Candidate> selected;
  selected.reserve(maximum);
  const bool nearest_order =
      selection == ReconstructionSelection::widest_accepted_nearest_prefix ||
      selection ==
          ReconstructionSelection::minimum_ghost_pressure_neumann_amplification;
  if (nearest_order) {
    selected.insert(selected.end(), candidates.begin(),
                    candidates.begin() + static_cast<std::ptrdiff_t>(maximum));
  } else {
    double maximum_normal_coordinate = 0.0;
    for (const auto &candidate : candidates) {
      const auto &cell = cells[static_cast<std::size_t>(candidate.id)];
      maximum_normal_coordinate =
          std::max(maximum_normal_coordinate,
                   dot(subtract(cell.center_m, point_m), frame.normal));
    }
    if (!(maximum_normal_coordinate > 0.0) ||
        !std::isfinite(maximum_normal_coordinate))
      throw runtime::Error("ghost plan directional donor extent is invalid");
    std::array<std::size_t, 12> directional{};
    directional.fill(candidates.size());
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
      const auto &cell = cells[static_cast<std::size_t>(candidates[index].id)];
      const auto offset = subtract(cell.center_m, point_m);
      const double normal_coordinate = dot(offset, frame.normal);
      const double tangent1_coordinate = dot(offset, frame.tangent1);
      const double tangent2_coordinate = dot(offset, frame.tangent2);
      const std::size_t layer = std::min<std::size_t>(
          2U, static_cast<std::size_t>(3.0 * normal_coordinate /
                                       maximum_normal_coordinate));
      const std::size_t quadrant = (tangent1_coordinate < 0.0 ? 2U : 0U) +
                                   (tangent2_coordinate < 0.0 ? 1U : 0U);
      const std::size_t slot = layer * 4U + quadrant;
      if (directional[slot] == candidates.size() ||
          std::tie(candidates[index].distance_squared, candidates[index].id) <
              std::tie(candidates[directional[slot]].distance_squared,
                       candidates[directional[slot]].id))
        directional[slot] = index;
    }
    std::vector<bool> used(candidates.size(), false);
    std::vector<std::size_t> representatives;
    representatives.reserve(directional.size());
    for (const auto index : directional)
      if (index != candidates.size())
        representatives.push_back(index);
    std::sort(representatives.begin(), representatives.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                return std::tie(candidates[lhs].distance_squared,
                                candidates[lhs].id) <
                       std::tie(candidates[rhs].distance_squared,
                                candidates[rhs].id);
              });
    for (const auto index : representatives) {
      if (!used[index] && selected.size() < maximum) {
        used[index] = true;
        selected.push_back(candidates[index]);
      }
    }
    for (std::size_t index = 0U;
         index < candidates.size() && selected.size() < maximum; ++index) {
      if (!used[index]) {
        used[index] = true;
        selected.push_back(candidates[index]);
      }
    }
    if (selected.size() != maximum)
      throw runtime::Error("ghost plan directional donor selection failed");
  }
  std::string last_error = "ghost plan has no accepted donor prefix";
  std::optional<CandidateReconstruction> best;
  std::optional<detail::ReconstructionFunctionalScore> best_score;
  std::vector<std::size_t> prefix_counts;
  prefix_counts.reserve(maximum - kMinimumDonors + 1U);
  if (selection == ReconstructionSelection::widest_accepted_nearest_prefix ||
      selection ==
          ReconstructionSelection::widest_accepted_directional_prefix) {
    for (std::size_t count = maximum;; --count) {
      prefix_counts.push_back(count);
      if (count == kMinimumDonors)
        break;
    }
  } else {
    for (std::size_t count = kMinimumDonors; count <= maximum; ++count)
      prefix_counts.push_back(count);
  }
  for (const auto count : prefix_counts) {
    std::vector<runtime::Int3> donors;
    donors.reserve(count);
    for (std::size_t index = 0U; index < count; ++index) {
      donors.push_back(selected[index].logical);
    }
    try {
      auto reconstruction = QuadraticReconstruction::create(
          point_m, frame.normal, frame.tangent1, frame.tangent2, scale,
          anchor.logical, donors, topology, geometry);
      if (selection ==
              ReconstructionSelection::widest_accepted_nearest_prefix ||
          selection ==
              ReconstructionSelection::widest_accepted_directional_prefix) {
        return {std::move(reconstruction), std::move(donors)};
      }
      double amplification = 0.0;
      if (selection ==
          ReconstructionSelection::minimum_wall_value_amplification) {
        amplification =
            weight_l1(detail::QuadraticReconstructionWeights::value_weights(
                reconstruction, point_m));
      } else if (selection == ReconstructionSelection::
                                  minimum_wall_normal_gradient_amplification) {
        amplification =
            scale *
            weight_l1(detail::QuadraticReconstructionWeights::
                          origin_constrained_directional_gradient_weights(
                              reconstruction, point_m, frame.normal));
      } else {
        if (!ghost_target.has_value())
          throw runtime::Error("ghost functional target is unavailable");
        const auto value_at_ghost =
            detail::QuadraticReconstructionWeights::value_weights(
                reconstruction, ghost_target->point_m);
        const auto value_at_wall =
            detail::QuadraticReconstructionWeights::value_weights(
                reconstruction, point_m);
        const auto cell_average_at_ghost =
            detail::QuadraticReconstructionWeights::cell_average_weights(
                reconstruction, ghost_target->logical, topology, geometry);
        const auto normal_gradient = detail::QuadraticReconstructionWeights::
            directional_gradient_weights(reconstruction, point_m, frame.normal);
        const double distance =
            dot(subtract(ghost_target->point_m, point_m), frame.normal);
        if (!std::isfinite(distance) || !(distance < 0.0))
          throw runtime::Error("ghost plan functional distance is invalid");
        double row_constraint_coefficient = 0.0;
        for (std::size_t index = 0U; index < cell_average_at_ghost.size();
             ++index) {
          if (cell_average_at_ghost[index].global_cell !=
              normal_gradient[index].global_cell)
            throw runtime::Error("ghost plan functional donor IDs disagree");
          const double coefficient = cell_average_at_ghost[index].weight -
                                     distance * normal_gradient[index].weight;
          if (!std::isfinite(coefficient))
            throw runtime::Error(
                "ghost plan pressure row coefficient is non-finite");
          if (cell_average_at_ghost[index].global_cell != anchor_id)
            row_constraint_coefficient += coefficient;
        }
        if (!(row_constraint_coefficient > 0.0) ||
            !std::isfinite(row_constraint_coefficient))
          throw runtime::Error(
              "ghost plan pressure row coefficient is not positive");
        amplification =
            std::max({weight_l1(value_at_ghost),
                      difference_weight_l1(value_at_ghost, value_at_wall),
                      combined_weight_l1(value_at_ghost, 1.0, normal_gradient,
                                         -distance),
                      combined_weight_l1(cell_average_at_ghost, 1.0,
                                         normal_gradient, -distance),
                      scale * weight_l1(normal_gradient)});
      }
      if (!std::isfinite(amplification)) {
        throw runtime::Error(
            "reconstruction functional amplification is non-finite");
      }
      const detail::ReconstructionFunctionalScore score{
          amplification, reconstruction.quality().condition_estimate, count,
          reconstruction.quality().pivot_fingerprint};
      bool replace = !best.has_value();
      if (best_score.has_value())
        replace = detail::functional_score_less(score, *best_score);
      if (replace) {
        best_score = score;
        best.emplace(CandidateReconstruction{std::move(reconstruction),
                                             std::move(donors)});
      }
    } catch (const runtime::Error &error) {
      last_error = error.what();
    }
  }
  if (best.has_value()) {
    return std::move(*best);
  }
  throw ReconstructionSelectionError(
      "ghost plan donor selection failed at anchor (" +
      std::to_string(anchor.logical.x) + "," +
      std::to_string(anchor.logical.y) + "," +
      std::to_string(anchor.logical.z) + ") with " + std::to_string(maximum) +
      " candidates: " + last_error);
}

CandidateReconstruction build_pressure_optimized_reconstruction(
    const ImmersedLink &link, const std::vector<CellRecord> &cells,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  if (link.solid_cell >= cells.size())
    throw runtime::Error(
        "pressure authority solid cell is outside the cell catalog");
  const auto &solid = cells[static_cast<std::size_t>(link.solid_cell)];
  const GhostFunctionalTarget target{solid.center_m, solid.logical};
  try {
    return build_reconstruction(
        link.wall_intercept_m, link.solid_to_fluid_normal, link.fluid_cell,
        cells, topology, geometry,
        ReconstructionSelection::minimum_ghost_pressure_neumann_amplification,
        target);
  } catch (const ReconstructionSelectionError &) {
    return build_reconstruction(
        link.wall_intercept_m, link.solid_to_fluid_normal, link.fluid_cell,
        cells, topology, geometry,
        ReconstructionSelection::
            minimum_directional_ghost_pressure_neumann_amplification,
        target);
  }
}

CandidateReconstruction rebuild_authority_at_point(
    runtime::Real3 point_m, runtime::Real3 normal, const ImmersedLink &link,
    const std::vector<runtime::Int3> &authoritative_donors,
    const std::vector<CellRecord> &cells, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    bool allow_tangent_plane_crossing = false) {
  if (link.fluid_cell >= cells.size())
    throw runtime::Error("boundary authority anchor is outside the catalog");
  const auto &anchor = cells[static_cast<std::size_t>(link.fluid_cell)];
  const auto frame = deterministic_frame(normal);
  const double scale = std::cbrt(anchor.volume_m3);
  std::optional<detail::BoundaryAuthorityCoverageScope> coverage_scope;
  if (allow_tangent_plane_crossing)
    coverage_scope.emplace();
  auto reconstruction = QuadraticReconstruction::create(
      point_m, frame.normal, frame.tangent1, frame.tangent2, scale,
      anchor.logical, authoritative_donors, topology, geometry);
  return {std::move(reconstruction), authoritative_donors};
}

double ghost_pressure_authority_score(
    const QuadraticReconstruction &reconstruction, const ImmersedLink &link,
    const std::vector<CellRecord> &cells, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry) {
  if (link.solid_cell >= cells.size() || link.fluid_cell >= cells.size())
    throw runtime::Error("boundary authority link is outside the catalog");
  const auto &solid = cells[static_cast<std::size_t>(link.solid_cell)];
  const auto value_at_ghost =
      detail::QuadraticReconstructionWeights::value_weights(reconstruction,
                                                            solid.center_m);
  const auto value_at_wall =
      detail::QuadraticReconstructionWeights::value_weights(
          reconstruction, link.wall_intercept_m);
  const auto cell_average_at_ghost =
      detail::QuadraticReconstructionWeights::cell_average_weights(
          reconstruction, solid.logical, topology, geometry);
  const auto normal_gradient =
      detail::QuadraticReconstructionWeights::directional_gradient_weights(
          reconstruction, link.wall_intercept_m, link.solid_to_fluid_normal);
  const double distance = dot(subtract(solid.center_m, link.wall_intercept_m),
                              link.solid_to_fluid_normal);
  if (!std::isfinite(distance) || !(distance < 0.0))
    throw runtime::Error("boundary authority pressure distance is invalid");
  double row_constraint_coefficient = 0.0;
  for (std::size_t index = 0U; index < cell_average_at_ghost.size(); ++index) {
    if (cell_average_at_ghost[index].global_cell !=
        normal_gradient[index].global_cell)
      throw runtime::Error("boundary authority donor IDs disagree");
    const double coefficient = cell_average_at_ghost[index].weight -
                               distance * normal_gradient[index].weight;
    if (cell_average_at_ghost[index].global_cell != link.fluid_cell)
      row_constraint_coefficient += coefficient;
  }
  if (!(row_constraint_coefficient > 0.0) ||
      !std::isfinite(row_constraint_coefficient))
    throw runtime::Error(
        "boundary authority pressure row coefficient is not positive");
  const double scale =
      std::cbrt(cells[static_cast<std::size_t>(link.fluid_cell)].volume_m3);
  return std::max(
      {weight_l1(value_at_ghost),
       difference_weight_l1(value_at_ghost, value_at_wall),
       combined_weight_l1(value_at_ghost, 1.0, normal_gradient, -distance),
       combined_weight_l1(cell_average_at_ghost, 1.0, normal_gradient,
                          -distance),
       scale * weight_l1(normal_gradient)});
}

std::vector<std::vector<runtime::Int3>>
build_pressure_authority_donor_candidates(
    const ImmersedLink &link, const std::vector<CellRecord> &cells,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  std::vector<std::vector<runtime::Int3>> donor_candidates;
  const auto same_donors = [](const auto &left, const auto &right) {
    if (left.size() != right.size())
      return false;
    for (std::size_t index = 0U; index < left.size(); ++index)
      if (left[index].x != right[index].x || left[index].y != right[index].y ||
          left[index].z != right[index].z)
        return false;
    return true;
  };
  const auto append = [&](CandidateReconstruction candidate) {
    auto donors = std::move(candidate.donor_cells);
    std::sort(donors.begin(), donors.end(), [&](auto left, auto right) {
      return topology.global_cell_id(left) < topology.global_cell_id(right);
    });
    if (std::find_if(donor_candidates.begin(), donor_candidates.end(),
                     [&](const auto &existing) {
                       return same_donors(existing, donors);
                     }) == donor_candidates.end())
      donor_candidates.push_back(std::move(donors));
  };
  append(
      build_pressure_optimized_reconstruction(link, cells, topology, geometry));
  append(build_reconstruction(
      link.wall_intercept_m, link.solid_to_fluid_normal, link.fluid_cell, cells,
      topology, geometry,
      ReconstructionSelection::widest_accepted_directional_prefix));
  append(build_reconstruction(
      link.wall_intercept_m, link.solid_to_fluid_normal, link.fluid_cell, cells,
      topology, geometry,
      ReconstructionSelection::minimum_wall_normal_gradient_amplification));
  return donor_candidates;
}

CandidateReconstruction build_pressure_authority_reconstruction(
    const ImmersedLink &link, const std::vector<CellRecord> &cells,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  const auto donor_candidates = build_pressure_authority_donor_candidates(
      link, cells, topology, geometry);

  std::optional<CandidateReconstruction> best;
  std::optional<detail::ReconstructionFunctionalScore> best_score;
  std::string last_candidate_error = "no donor candidate was evaluated";
  const double scale =
      std::cbrt(cells[static_cast<std::size_t>(link.fluid_cell)].volume_m3);
  for (const auto &donors : donor_candidates) {
    try {
      auto link_candidate = rebuild_authority_at_point(
          link.wall_intercept_m, link.solid_to_fluid_normal, link, donors,
          cells, topology, geometry);
      double amplification = ghost_pressure_authority_score(
          link_candidate.reconstruction, link, cells, topology, geometry);
      amplification = std::max(
          amplification,
          scale * weight_l1(detail::QuadraticReconstructionWeights::
                                origin_constrained_directional_gradient_weights(
                                    link_candidate.reconstruction,
                                    link.wall_intercept_m,
                                    link.solid_to_fluid_normal)));
      const detail::ReconstructionFunctionalScore score{
          amplification,
          link_candidate.reconstruction.quality().condition_estimate,
          donors.size(),
          link_candidate.reconstruction.quality().pivot_fingerprint};
      if (!best_score.has_value() ||
          detail::functional_score_less(score, *best_score)) {
        best_score = score;
        best = std::move(link_candidate);
      }
    } catch (const runtime::Error &error) {
      last_candidate_error = error.what();
    }
  }
  if (!best.has_value())
    throw ReconstructionSelectionError(
        "boundary authority has no common pressure/velocity donor set: " +
        last_candidate_error);
  return std::move(*best);
}

std::vector<WeightedDonor>
combine_weights(const std::vector<WeightedDonor> &first, double first_scale,
                const std::vector<WeightedDonor> &second, double second_scale) {
  if (first.size() != second.size()) {
    throw runtime::Error("ghost plan functional donor layouts disagree");
  }
  std::vector<WeightedDonor> result(first.size());
  for (std::size_t index = 0U; index < result.size(); ++index) {
    if (first[index].global_cell != second[index].global_cell) {
      throw runtime::Error("ghost plan functional donor IDs disagree");
    }
    const double weight =
        first_scale * first[index].weight + second_scale * second[index].weight;
    if (!std::isfinite(weight)) {
      throw runtime::Error("ghost plan functional weight is non-finite");
    }
    result[index] = {first[index].global_cell, weight};
  }
  return result;
}

void hash_donors(std::uint64_t &hash,
                 const std::vector<WeightedDonor> &donors) noexcept {
  hash_u64(hash, donors.size());
  for (const auto &donor : donors) {
    hash_u64(hash, donor.global_cell);
    hash_u64(hash, double_bits(donor.weight));
  }
}

} // namespace

namespace detail {

struct GhostEntry final {
  ImmersedLinkId link{};
  ImmersedLink record;
  QuadraticReconstruction reconstruction;
  std::array<AffineGhostConstraint, 3> velocity;
  AffineGhostConstraint zero_normal;
  FluidExtrapolation density;
  runtime::Real3 surface_measure_vector_m2{};
  runtime::Real3 surface_patch_centroid_m{};
};

struct GhostStencilPlanStorage final {
  std::vector<GhostEntry> entries;
  std::uint32_t maximum_halo_reach{};
  std::uint64_t fingerprint{};
};

struct WallQuadraturePlanStorage final {
  std::vector<WallQuadraturePoint> local_points;
  std::uint64_t fingerprint{};
};

} // namespace detail

namespace {

[[noreturn]] void
throw_collective_failure(const runtime::CollectiveStatus &status) {
  throw runtime::Error(status.message + " (lowest failing rank " +
                       std::to_string(status.failing_rank) + ")");
}

const detail::GhostEntry &
find_entry(const detail::GhostStencilPlanStorage &storage,
           ImmersedLinkId link) {
  const auto iterator =
      std::lower_bound(storage.entries.begin(), storage.entries.end(), link,
                       [](const auto &entry, ImmersedLinkId value) {
                         return entry.link < value;
                       });
  if (iterator == storage.entries.end() || iterator->link != link) {
    throw runtime::Error("ghost plan link ID is outside the plan");
  }
  return *iterator;
}

std::vector<runtime::Box3> owner_boxes(const std::vector<CellRecord> &cells,
                                       const runtime::MpiContext &mpi) {
  std::vector<runtime::Box3> boxes(static_cast<std::size_t>(mpi.size()));
  std::vector<std::size_t> counts(static_cast<std::size_t>(mpi.size()), 0U);
  for (auto &box : boxes) {
    box.begin = {std::numeric_limits<int>::max(),
                 std::numeric_limits<int>::max(),
                 std::numeric_limits<int>::max()};
    box.end = {std::numeric_limits<int>::min(), std::numeric_limits<int>::min(),
               std::numeric_limits<int>::min()};
  }
  for (const auto &cell : cells) {
    auto &box = boxes[static_cast<std::size_t>(cell.owner_rank)];
    box.begin.x = std::min(box.begin.x, cell.logical.x);
    box.begin.y = std::min(box.begin.y, cell.logical.y);
    box.begin.z = std::min(box.begin.z, cell.logical.z);
    box.end.x = std::max(box.end.x, cell.logical.x + 1);
    box.end.y = std::max(box.end.y, cell.logical.y + 1);
    box.end.z = std::max(box.end.z, cell.logical.z + 1);
    ++counts[static_cast<std::size_t>(cell.owner_rank)];
  }
  for (std::size_t rank = 0U; rank < boxes.size(); ++rank) {
    if (counts[rank] == 0U || boxes[rank].begin.x >= boxes[rank].end.x ||
        boxes[rank].begin.y >= boxes[rank].end.y ||
        boxes[rank].begin.z >= boxes[rank].end.z) {
      throw runtime::Error("wall quadrature owner box is empty");
    }
  }
  return boxes;
}

void hash_quality(std::uint64_t &hash,
                  const ReconstructionQuality &quality) noexcept {
  hash_u64(hash, quality.rank);
  hash_u64(hash, double_bits(quality.condition_estimate));
  hash_u64(hash, quality.halo_reach);
  hash_u64(hash, quality.pivot_fingerprint);
}

} // namespace

GhostStencilPlan GhostStencilPlan::create(
    const ImmersedSurface &surface, const SurfaceQuery &query,
    const ImmersedDomain &domain, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const runtime::StructuredDecomposition &decomposition,
    const runtime::MpiContext &mpi) {
  validate_common_inputs(surface, query, domain, topology, geometry, mpi);
  validate_decomposition(topology, decomposition, mpi);
  const auto cells = build_global_cells(domain, topology, geometry, mpi);
  const auto links = build_global_links(domain, cells, surface, mpi);
  const auto surface_measure_by_link =
      build_surface_measure_by_link(surface, links, cells);

  std::shared_ptr<detail::GhostStencilPlanStorage> storage;
  bool local_ok = true;
  std::string local_message;
  try {
    storage = std::make_shared<detail::GhostStencilPlanStorage>();
    storage->entries.reserve(links.size());
    std::uint64_t fingerprint = kFnvOffset;
    hash_u64(fingerprint, UINT64_C(0x48554e4447535033));
    hash_u64(fingerprint, surface.fingerprint());
    hash_u64(fingerprint, query.fingerprint());
    hash_u64(fingerprint, domain.classification_fingerprint());
    hash_u64(fingerprint, domain.surface_coverage_fingerprint());
    hash_u64(fingerprint, links.size());
    struct PlannedGhost final {
      detail::GhostEntry entry;
      std::uint32_t required_reach{};
    };
    const auto build_link = [&](std::size_t link_index) {
    const auto &link = links[link_index].record;
    const auto &solid = cells[static_cast<std::size_t>(link.solid_cell)];
    const auto surface_measure =
        surface_measure_by_link.measures[link_index];
    const auto surface_centroid =
        surface_measure_by_link.centroids[link_index];
    if (!finite(surface_measure))
      throw runtime::Error("ghost plan link surface measure is non-finite");
    if (!finite(surface_centroid))
      throw runtime::Error("ghost plan link surface centroid is non-finite");
      auto candidate = build_pressure_authority_reconstruction(
          link, cells, topology, geometry);
      const auto value_at_ghost =
          detail::QuadraticReconstructionWeights::value_weights(
              candidate.reconstruction, solid.center_m);
      const auto velocity_value_at_wall =
          detail::QuadraticReconstructionWeights::value_weights(
              candidate.reconstruction, link.wall_intercept_m);
      const auto normal_gradient =
          detail::QuadraticReconstructionWeights::directional_gradient_weights(
              candidate.reconstruction, link.wall_intercept_m,
              link.solid_to_fluid_normal);
      const double distance =
          dot(subtract(solid.center_m, link.wall_intercept_m),
              link.solid_to_fluid_normal);
      if (!std::isfinite(distance) || distance >= 0.0) {
        throw runtime::Error(
            "ghost plan solid-neighbour centre is not on the solid side");
      }
      AffineGhostConstraint velocity{
          link.id,
          combine_weights(value_at_ghost, 1.0, velocity_value_at_wall, -1.0),
          1.0, 0.0};
      AffineGhostConstraint zero_normal{
          link.id,
          combine_weights(value_at_ghost, 1.0, normal_gradient, -distance), 0.0,
          distance};
      FluidExtrapolation density{link.id, value_at_ghost};
      std::array<AffineGhostConstraint, 3> velocity_components{
          velocity, velocity, std::move(velocity)};
      const auto quality = candidate.reconstruction.quality();
      const auto cloud_reach = donor_cloud_reach(candidate.donor_cells);
      const auto interface_reach =
          cloud_reach <= 4U ? cloud_reach : quality.halo_reach;
      const auto required_reach = std::max(quality.halo_reach, interface_reach);
      return PlannedGhost{{link.id, link, std::move(candidate.reconstruction),
                           std::move(velocity_components),
                           std::move(zero_normal), std::move(density),
                           surface_measure, surface_centroid},
                          required_reach};
    };
    const auto consume_link = [&](PlannedGhost &planned) {
      const auto &entry = planned.entry;
      const auto &link = entry.record;
      const auto quality = entry.reconstruction.quality();
      storage->maximum_halo_reach =
          std::max(storage->maximum_halo_reach, planned.required_reach);

      hash_u64(fingerprint, link.id);
      hash_u64(fingerprint, link.fluid_cell);
      hash_u64(fingerprint, link.solid_cell);
      hash_real3(fingerprint, link.wall_intercept_m);
      hash_real3(fingerprint, link.solid_to_fluid_normal);
      hash_real3(fingerprint, entry.surface_measure_vector_m2);
      hash_real3(fingerprint, entry.surface_patch_centroid_m);
      hash_quality(fingerprint, quality);
      hash_donors(fingerprint, entry.velocity[0].donors);
      hash_u64(fingerprint, double_bits(entry.velocity[0].wall_value_weight));
      hash_u64(fingerprint,
               double_bits(entry.velocity[0].wall_normal_gradient_weight_m));
      hash_donors(fingerprint, entry.zero_normal.donors);
      hash_u64(fingerprint, double_bits(entry.zero_normal.wall_value_weight));
      hash_u64(fingerprint,
               double_bits(entry.zero_normal.wall_normal_gradient_weight_m));
      hash_donors(fingerprint, entry.density.donors);

      storage->entries.push_back(std::move(planned.entry));
    };
    const auto hardware_threads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::size_t rank_thread_budget = std::max<std::size_t>(
        1U, kParallelWorkerBudget / static_cast<std::size_t>(mpi.size()));
    const std::size_t worker_count = std::min<std::size_t>(
        std::min<std::size_t>(hardware_threads, rank_thread_budget),
        links.size());
    std::vector<std::optional<PlannedGhost>> planned(links.size());
    std::vector<std::exception_ptr> errors(links.size());
    std::atomic<std::size_t> next{0U};
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    try {
      for (std::size_t worker = 0U; worker < worker_count; ++worker) {
        workers.emplace_back([&] {
          for (;;) {
            const std::size_t link_index =
                next.fetch_add(1U, std::memory_order_relaxed);
            if (link_index >= links.size())
              return;
            try {
              detail::ValidatedGeometryScope validated_geometry(topology,
                                                                geometry);
              planned[link_index].emplace(build_link(link_index));
            } catch (...) {
              errors[link_index] = std::current_exception();
            }
          }
        });
      }
    } catch (...) {
      next.store(links.size(), std::memory_order_relaxed);
      for (auto &worker : workers)
        worker.join();
      throw;
    }
    for (auto &worker : workers)
      worker.join();
    for (std::size_t link_index = 0U; link_index < links.size(); ++link_index) {
      if (errors[link_index])
        std::rethrow_exception(errors[link_index]);
      if (!planned[link_index].has_value())
        throw runtime::Error("ghost plan parallel link result is missing");
      consume_link(*planned[link_index]);
    }
    hash_u64(fingerprint, storage->maximum_halo_reach);
    storage->fingerprint = fingerprint;
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "ghost plan construction failed";
  }
  const auto status = runtime::collective_status(mpi, local_ok, local_message);
  if (!status.ok) {
    throw_collective_failure(status);
  }

  std::uint64_t minimum = storage->fingerprint;
  std::uint64_t maximum = storage->fingerprint;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                          mpi.comm()),
            "MPI_Allreduce ghost-plan minimum fingerprint");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce ghost-plan maximum fingerprint");
  require_collective(mpi, minimum == maximum,
                     "ghost plan: rank fingerprints disagree");
  return GhostStencilPlan(std::move(storage));
}

const AffineGhostConstraint &
GhostStencilPlan::velocity_constraint(ImmersedLinkId link,
                                      std::size_t component) const {
  if (component >= 3U) {
    throw runtime::Error("ghost plan velocity component is out of bounds");
  }
  return find_entry(*storage_, link).velocity[component];
}

const AffineGhostConstraint &
GhostStencilPlan::zero_normal_constraint(ImmersedLinkId link) const {
  return find_entry(*storage_, link).zero_normal;
}

const FluidExtrapolation &
GhostStencilPlan::density_extrapolation(ImmersedLinkId link) const {
  return find_entry(*storage_, link).density;
}

const QuadraticReconstruction &
GhostStencilPlan::reconstruction(ImmersedLinkId link) const {
  return find_entry(*storage_, link).reconstruction;
}

std::uint32_t GhostStencilPlan::maximum_halo_reach() const noexcept {
  return storage_->maximum_halo_reach;
}

std::uint64_t GhostStencilPlan::fingerprint() const noexcept {
  return storage_->fingerprint;
}

std::size_t GhostStencilPlan::immersed_operator_link_count() const noexcept {
  return storage_->entries.size();
}

const ImmersedLink &
GhostStencilPlan::link_for_immersed_operator(ImmersedLinkId link) const {
  return find_entry(*storage_, link).record;
}

runtime::Real3 GhostStencilPlan::
    surface_measure_vector_m2_for_immersed_operator(ImmersedLinkId link) const {
  return find_entry(*storage_, link).surface_measure_vector_m2;
}

runtime::Real3
GhostStencilPlan::surface_patch_centroid_m_for_immersed_operator(
    ImmersedLinkId link) const {
  return find_entry(*storage_, link).surface_patch_centroid_m;
}

WallQuadraturePlan WallQuadraturePlan::create(
    const ImmersedSurface &surface, const SurfaceQuery &query,
    const ImmersedDomain &domain, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry, const runtime::MpiContext &mpi) {
  validate_common_inputs(surface, query, domain, topology, geometry, mpi);
  const auto cells = build_global_cells(domain, topology, geometry, mpi);
  const auto links = build_global_links(domain, cells, surface, mpi);

  std::shared_ptr<detail::WallQuadraturePlanStorage> storage;
  bool local_ok = true;
  std::string local_message;
  try {
    const auto boxes = owner_boxes(cells, mpi);
    using SpatialBin = std::array<std::int64_t, 3>;
    double maximum_link_h = 0.0;
    for (const auto &global_link : links) {
      const auto fluid_cell = global_link.record.fluid_cell;
      if (fluid_cell >= cells.size())
        throw runtime::Error("wall quadrature link fluid cell is out of range");
      const double h =
          std::cbrt(cells[static_cast<std::size_t>(fluid_cell)].volume_m3);
      if (!(h > 0.0) || !std::isfinite(h))
        throw runtime::Error("wall quadrature link scale is invalid");
      maximum_link_h = std::max(maximum_link_h, h);
    }
    if (!(maximum_link_h > 0.0) || !std::isfinite(maximum_link_h))
      throw runtime::Error("wall quadrature maximum link scale is invalid");
    const double bin_width = 2.0 * maximum_link_h;
    const auto bin_coordinate = [bin_width](double value) {
      const double coordinate = std::floor(value / bin_width);
      if (!std::isfinite(coordinate) ||
          coordinate <
              static_cast<double>(std::numeric_limits<std::int64_t>::min()) ||
          coordinate >
              static_cast<double>(std::numeric_limits<std::int64_t>::max()))
        throw runtime::Error(
            "wall quadrature spatial-bin coordinate is invalid");
      return static_cast<std::int64_t>(coordinate);
    };
    std::map<SpatialBin, std::vector<std::size_t>> link_bins;
    for (std::size_t link_index = 0U; link_index < links.size(); ++link_index) {
      const auto &link = links[link_index].record;
      const double h =
          std::cbrt(cells[static_cast<std::size_t>(link.fluid_cell)].volume_m3);
      const double radius = 2.0 * h;
      const SpatialBin first{bin_coordinate(link.wall_intercept_m.x - radius),
                             bin_coordinate(link.wall_intercept_m.y - radius),
                             bin_coordinate(link.wall_intercept_m.z - radius)};
      const SpatialBin last{bin_coordinate(link.wall_intercept_m.x + radius),
                            bin_coordinate(link.wall_intercept_m.y + radius),
                            bin_coordinate(link.wall_intercept_m.z + radius)};
      std::array<int, 3> span{};
      for (std::size_t axis = 0U; axis < span.size(); ++axis) {
        const auto difference = last[axis] - first[axis];
        if (difference < 0 || difference > 3)
          throw runtime::Error("wall quadrature spatial-bin span is invalid");
        span[axis] = static_cast<int>(difference);
      }
      for (int k = 0; k <= span[2]; ++k)
        for (int j = 0; j <= span[1]; ++j)
          for (int i = 0; i <= span[0]; ++i)
            link_bins[{first[0] + i, first[1] + j, first[2] + k}].push_back(
                link_index);
    }
    storage = std::make_shared<detail::WallQuadraturePlanStorage>();
    auto mutable_authority_catalog =
        std::make_shared<std::vector<detail::BoundaryAuthorityOwner>>();
    mutable_authority_catalog->reserve(links.size());
    for (const auto &global_link : links) {
      mutable_authority_catalog->push_back(
          {global_link.record.id,
           cells[static_cast<std::size_t>(global_link.record.fluid_cell)]
               .owner_rank});
    }
    const std::shared_ptr<const std::vector<detail::BoundaryAuthorityOwner>>
        authority_catalog = std::move(mutable_authority_catalog);
    std::vector<std::optional<CandidateReconstruction>>
        pressure_authority_by_link(links.size());
    std::vector<std::exception_ptr> pressure_authority_errors(links.size());
    {
      const auto hardware_threads =
          std::max(1U, std::thread::hardware_concurrency());
      const std::size_t rank_thread_budget = std::max<std::size_t>(
          1U, kParallelWorkerBudget / static_cast<std::size_t>(mpi.size()));
      const std::size_t worker_count = std::min<std::size_t>(
          std::min<std::size_t>(hardware_threads, rank_thread_budget),
          links.size());
      std::atomic<std::size_t> next_link{0U};
      std::vector<std::thread> workers;
      workers.reserve(worker_count);
      try {
        for (std::size_t worker = 0U; worker < worker_count; ++worker) {
          workers.emplace_back([&] {
            for (;;) {
              const std::size_t link_index =
                  next_link.fetch_add(1U, std::memory_order_relaxed);
              if (link_index >= links.size())
                return;
              try {
                detail::ValidatedGeometryScope validated_geometry(topology,
                                                                  geometry);
                pressure_authority_by_link[link_index].emplace(
                    build_pressure_authority_reconstruction(
                        links[link_index].record, cells, topology, geometry));
              } catch (...) {
                pressure_authority_errors[link_index] =
                    std::current_exception();
              }
            }
          });
        }
      } catch (...) {
        next_link.store(links.size(), std::memory_order_relaxed);
        for (auto &worker : workers)
          worker.join();
        throw;
      }
      for (auto &worker : workers)
        worker.join();
    }
    for (std::size_t link_index = 0U; link_index < links.size(); ++link_index) {
      if (pressure_authority_errors[link_index])
        std::rethrow_exception(pressure_authority_errors[link_index]);
      if (!pressure_authority_by_link[link_index].has_value())
        throw runtime::Error(
            "wall quadrature boundary authority result is missing");
    }
    std::map<mesh::GlobalCellId, std::vector<std::size_t>> links_by_fluid_row;
    for (std::size_t link_index = 0U; link_index < links.size(); ++link_index)
      links_by_fluid_row[links[link_index].record.fluid_cell].push_back(
          link_index);
    const auto global_extent = topology.global_extent();
    constexpr std::array<runtime::Int3, 6> neighbour_offsets{
        runtime::Int3{-1, 0, 0}, runtime::Int3{1, 0, 0},
        runtime::Int3{0, -1, 0}, runtime::Int3{0, 1, 0},
        runtime::Int3{0, 0, -1}, runtime::Int3{0, 0, 1}};
    for (const auto &[fluid_row, row_links] : links_by_fluid_row) {
      if (fluid_row >= cells.size() || row_links.empty())
        throw runtime::Error(
            "wall quadrature boundary row association is invalid");
      std::vector<std::pair<mesh::GlobalCellId, runtime::Int3>> ordered_donors;
      for (const auto link_index : row_links)
        for (const auto donor :
             pressure_authority_by_link[link_index]->donor_cells)
          ordered_donors.push_back({topology.global_cell_id(donor), donor});
      const auto row_logical = cells[static_cast<std::size_t>(fluid_row)].logical;
      for (const auto offset : neighbour_offsets) {
        const runtime::Int3 neighbour{row_logical.x + offset.x,
                                      row_logical.y + offset.y,
                                      row_logical.z + offset.z};
        if (neighbour.x < 0 || neighbour.x >= global_extent.x ||
            neighbour.y < 0 || neighbour.y >= global_extent.y ||
            neighbour.z < 0 || neighbour.z >= global_extent.z)
          continue;
        const auto neighbour_id = topology.global_cell_id(neighbour);
        if (neighbour_id < cells.size() &&
            cells[static_cast<std::size_t>(neighbour_id)].region ==
                CellRegion::fluid)
          ordered_donors.push_back({neighbour_id, neighbour});
      }
      std::sort(ordered_donors.begin(), ordered_donors.end(),
                [](const auto &left, const auto &right) {
                  return left.first < right.first;
                });
      ordered_donors.erase(
          std::unique(ordered_donors.begin(), ordered_donors.end(),
                      [](const auto &left, const auto &right) {
                        return left.first == right.first;
                      }),
          ordered_donors.end());
      std::vector<runtime::Int3> shared_donors;
      shared_donors.reserve(ordered_donors.size());
      for (const auto &[id, donor] : ordered_donors) {
        static_cast<void>(id);
        shared_donors.push_back(donor);
      }
      for (const auto link_index : row_links) {
        const auto &link = links[link_index].record;
        pressure_authority_by_link[link_index].emplace(
            rebuild_authority_at_point(
                link.wall_intercept_m, link.solid_to_fluid_normal, link,
                shared_donors, cells, topology, geometry, true));
      }
    }
    std::uint64_t fingerprint = kFnvOffset;
    hash_u64(fingerprint, UINT64_C(0x48554e4457515033));
    hash_u64(fingerprint, surface.fingerprint());
    hash_u64(fingerprint, query.fingerprint());
    hash_u64(fingerprint, domain.classification_fingerprint());
    hash_u64(fingerprint, domain.surface_coverage_fingerprint());
    hash_u64(fingerprint, surface.triangle_count());

    constexpr std::array<std::array<double, 3>, 3> kBarycentric{{
        {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
        {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
        {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
    }};
    struct PointAssociation final {
      runtime::Real3 point{};
      const ImmersedLink *link{};
    };
    struct PlannedPoint final {
      runtime::Real3 point{};
      runtime::Real3 normal{};
      runtime::Real3 pressure_authority_normal{};
      double weight{};
      int owner{};
      int authority_owner{};
      ImmersedLinkId link{};
      QuadraticReconstruction point_reconstruction;
      QuadraticReconstruction pressure_authority_reconstruction;
    };
    struct PlannedTriangle final {
      std::array<std::optional<PlannedPoint>, 3> points;
    };
    const auto build_triangle = [&](TriangleId triangle_id) {
      PlannedTriangle planned;
      const auto &triangle = surface.triangle(triangle_id);
      if (!finite(triangle.geometric_outward_normal) ||
          !std::isfinite(triangle.area_m2) || triangle.area_m2 <= 0.0) {
        throw runtime::Error("wall quadrature triangle geometry is invalid");
      }
      std::array<PointAssociation, 3> associations{};
      for (std::size_t point_index = 0U; point_index < associations.size();
           ++point_index) {
        runtime::Real3 point{};
        for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
          point = add(point, multiply(triangle.vertices_m[vertex],
                                      kBarycentric[point_index][vertex]));
        }
        const ImmersedLink *selected = nullptr;
        double selected_distance = std::numeric_limits<double>::infinity();
        const SpatialBin point_bin{bin_coordinate(point.x),
                                   bin_coordinate(point.y),
                                   bin_coordinate(point.z)};
        const auto candidates = link_bins.find(point_bin);
        if (candidates == link_bins.end())
          throw runtime::Error(
              "wall quadrature point has no spatial-link candidates");
        for (const auto link_index : candidates->second) {
          const auto &link = links[link_index].record;
          const auto &fluid = cells[static_cast<std::size_t>(link.fluid_cell)];
          const double distance =
              std::sqrt(squared_distance(point, link.wall_intercept_m));
          const double h = std::cbrt(fluid.volume_m3);
          const double normal_alignment = dot(triangle.geometric_outward_normal,
                                              link.solid_to_fluid_normal);
          const double alignment_tolerance =
              512.0 * std::numeric_limits<double>::epsilon();
          if (std::isfinite(distance) && std::isfinite(h) && h > 0.0 &&
              std::isfinite(normal_alignment) &&
              std::abs(normal_alignment) > alignment_tolerance &&
              distance <= 2.0 * h &&
              (distance < selected_distance ||
               (distance == selected_distance &&
                (selected == nullptr || link.id < selected->id)))) {
            selected = &link;
            selected_distance = distance;
          }
        }
        if (selected == nullptr) {
          throw runtime::Error(
              "wall quadrature point has no associated interface link");
        }
        associations[point_index] = {point, selected};
      }

      double alignment = dot(triangle.geometric_outward_normal,
                             associations.front().link->solid_to_fluid_normal);
      if (!std::isfinite(alignment) || alignment == 0.0) {
        throw runtime::Error("wall quadrature normal orientation is ambiguous");
      }
      const auto normal = multiply(triangle.geometric_outward_normal,
                                   alignment > 0.0 ? 1.0 : -1.0);
      for (const auto &association : associations) {
        if (dot(normal, association.link->solid_to_fluid_normal) <= 0.0) {
          throw runtime::Error(
              "wall quadrature associated normal is inconsistent");
        }
      }

      std::vector<runtime::Int3> point_donors;
      std::vector<runtime::Int3> pressure_authority_donors;
      const auto append_logical_donors = [&](
                                               const auto &reconstruction,
                                               auto &destination) {
        const auto &donor_ids =
            detail::QuadraticReconstructionWeights::donor_global_ids(
                reconstruction);
        destination.reserve(destination.size() + donor_ids.size());
        for (const auto donor : donor_ids) {
          if (donor >= cells.size())
            throw runtime::Error(
                "wall quadrature donor exceeds owner Halo reach");
          destination.push_back(
              cells[static_cast<std::size_t>(donor)].logical);
        }
      };
      for (std::size_t point_index = 0U; point_index < associations.size();
           ++point_index) {
        const auto &association = associations[point_index];
        const auto &pressure_authority_candidate =
            *pressure_authority_by_link[association.link->id];
        auto point_authority_candidate = rebuild_authority_at_point(
            association.point, normal, *association.link,
            pressure_authority_candidate.donor_cells, cells, topology, geometry,
            true);
        append_logical_donors(point_authority_candidate.reconstruction,
                              point_donors);
        append_logical_donors(pressure_authority_candidate.reconstruction,
                              pressure_authority_donors);
        const double weight = triangle.area_m2 / 3.0;
        if (!finite(association.point) || !finite(normal) ||
            !std::isfinite(weight) || weight <= 0.0) {
          throw runtime::Error("wall quadrature point is invalid");
        }
        planned.points[point_index].emplace(PlannedPoint{
            association.point, normal, association.link->solid_to_fluid_normal,
            weight, mpi.size(),
            cells[static_cast<std::size_t>(association.link->fluid_cell)]
                .owner_rank,
            association.link->id, point_authority_candidate.reconstruction,
            pressure_authority_candidate.reconstruction});
      }
      const int owner = detail::select_wall_quadrature_execution_owner(
          boxes, point_donors, pressure_authority_donors, kMaximumReach);
      for (auto &point : planned.points) {
        if (!point.has_value())
          throw runtime::Error("wall quadrature planned point is missing");
        point->owner = owner;
      }
      return planned;
    };
    const auto consume_triangle = [&](TriangleId triangle_id,
                                      PlannedTriangle &planned) {
      for (std::size_t point_index = 0U; point_index < planned.points.size();
           ++point_index) {
        if (!planned.points[point_index].has_value())
          throw runtime::Error("wall quadrature planned point is missing");
        auto &point = *planned.points[point_index];
        hash_u64(fingerprint, triangle_id);
        hash_u64(fingerprint, point_index);
        hash_real3(fingerprint, point.point);
        hash_real3(fingerprint, point.normal);
        hash_real3(fingerprint, point.pressure_authority_normal);
        hash_u64(fingerprint, double_bits(point.weight));
        hash_u64(fingerprint, point.link);
        hash_quality(fingerprint, point.point_reconstruction.quality());
        hash_donors(fingerprint,
                    detail::QuadraticReconstructionWeights::value_weights(
                        point.point_reconstruction, point.point));
        hash_quality(fingerprint, point.point_reconstruction.quality());
        hash_donors(fingerprint,
                    detail::QuadraticReconstructionWeights::value_weights(
                        point.point_reconstruction, point.point));
        hash_quality(fingerprint,
                     point.pressure_authority_reconstruction.quality());
        hash_donors(fingerprint,
                    detail::QuadraticReconstructionWeights::value_weights(
                        point.pressure_authority_reconstruction, point.point));
        hash_quality(fingerprint, point.point_reconstruction.quality());
        for (const auto direction :
             {runtime::Real3{1.0, 0.0, 0.0}, runtime::Real3{0.0, 1.0, 0.0},
              runtime::Real3{0.0, 0.0, 1.0}}) {
          hash_donors(
              fingerprint,
              detail::QuadraticReconstructionWeights::
                  origin_constrained_directional_gradient_weights(
                      point.point_reconstruction, point.point, direction));
        }
        if (mpi.rank() == point.owner) {
          auto reconstruction = detail::with_boundary_authority(
              point.point_reconstruction, point.link, point.authority_owner,
              point.pressure_authority_reconstruction, authority_catalog);
          storage->local_points.push_back(
              {triangle_id, static_cast<std::uint32_t>(point_index),
               point.point, point.normal, point.weight, point.owner,
               std::move(reconstruction)});
        }
      }
    };
    const auto hardware_threads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::size_t rank_thread_budget = std::max<std::size_t>(
        1U, kParallelWorkerBudget / static_cast<std::size_t>(mpi.size()));
    const std::size_t worker_count = std::min<std::size_t>(
        std::min<std::size_t>(hardware_threads, rank_thread_budget),
        surface.triangle_count());
    std::vector<std::optional<PlannedTriangle>> planned(
        surface.triangle_count());
    std::vector<std::exception_ptr> errors(surface.triangle_count());
    std::atomic<std::size_t> next{0U};
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    try {
      for (std::size_t worker = 0U; worker < worker_count; ++worker) {
        workers.emplace_back([&] {
          for (;;) {
            const std::size_t triangle_index =
                next.fetch_add(1U, std::memory_order_relaxed);
            if (triangle_index >= surface.triangle_count())
              return;
            try {
              detail::ValidatedGeometryScope validated_geometry(topology,
                                                                geometry);
              planned[triangle_index].emplace(build_triangle(triangle_index));
            } catch (...) {
              errors[triangle_index] = std::current_exception();
            }
          }
        });
      }
    } catch (...) {
      next.store(surface.triangle_count(), std::memory_order_relaxed);
      for (auto &worker : workers)
        worker.join();
      throw;
    }
    for (auto &worker : workers)
      worker.join();
    for (TriangleId triangle = 0U; triangle < surface.triangle_count();
         ++triangle) {
      if (errors[triangle])
        std::rethrow_exception(errors[triangle]);
      if (!planned[triangle].has_value())
        throw runtime::Error(
            "wall quadrature parallel triangle result is missing");
      consume_triangle(triangle, *planned[triangle]);
    }
    storage->fingerprint = fingerprint;
  } catch (const std::exception &error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "wall quadrature plan construction failed";
  }
  const auto status = runtime::collective_status(mpi, local_ok, local_message);
  if (!status.ok) {
    throw_collective_failure(status);
  }

  std::uint64_t minimum = storage->fingerprint;
  std::uint64_t maximum = storage->fingerprint;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                          mpi.comm()),
            "MPI_Allreduce wall-plan minimum fingerprint");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce wall-plan maximum fingerprint");
  require_collective(mpi, minimum == maximum,
                     "wall quadrature plan: rank fingerprints disagree");
  return WallQuadraturePlan(std::move(storage));
}

const std::vector<WallQuadraturePoint> &
WallQuadraturePlan::local_points() const noexcept {
  return storage_->local_points;
}

std::uint32_t WallQuadraturePlan::maximum_halo_reach() const noexcept {
  return static_cast<std::uint32_t>(kMaximumReach);
}

int detail::select_wall_quadrature_execution_owner(
    const std::vector<runtime::Box3> &boxes,
    const std::vector<runtime::Int3> &point_donors,
    const std::vector<runtime::Int3> &pressure_authority_donors, int reach) {
  const auto donor_is_covered = [reach](runtime::Int3 donor,
                                        runtime::Box3 box) noexcept {
    return reach >= 0 && donor.x >= box.begin.x - reach &&
           donor.x < box.end.x + reach &&
           donor.y >= box.begin.y - reach && donor.y < box.end.y + reach &&
           donor.z >= box.begin.z - reach && donor.z < box.end.z + reach;
  };
  for (std::size_t rank = 0U; rank < boxes.size(); ++rank) {
    const auto covers = [&](const std::vector<runtime::Int3> &donors) {
      return std::all_of(donors.begin(), donors.end(), [&](const auto donor) {
        return donor_is_covered(donor, boxes[rank]);
      });
    };
    if (covers(point_donors) && covers(pressure_authority_donors))
      return static_cast<int>(rank);
  }
  throw runtime::Error("wall quadrature donor exceeds owner Halo reach");
}

std::size_t detail::select_minimum_functional_score(
    const std::vector<ReconstructionFunctionalScore> &scores) {
  if (scores.empty())
    throw runtime::Error("reconstruction functional score list is empty");
  std::size_t selected = 0U;
  for (std::size_t index = 0U; index < scores.size(); ++index) {
    const auto &score = scores[index];
    if (!(score.amplification >= 0.0) || !std::isfinite(score.amplification) ||
        !(score.condition_estimate >= 0.0) ||
        !std::isfinite(score.condition_estimate) || score.donor_count == 0U)
      throw runtime::Error("reconstruction functional score is invalid");
    if (functional_score_less(score, scores[selected]))
      selected = index;
  }
  return selected;
}

bool detail::functional_score_less(const ReconstructionFunctionalScore &left,
                                   const ReconstructionFunctionalScore &right) {
  return std::tuple{left.amplification, left.condition_estimate,
                    left.donor_count, left.pivot_fingerprint} <
         std::tuple{right.amplification, right.condition_estimate,
                    right.donor_count, right.pivot_fingerprint};
}

std::uint64_t WallQuadraturePlan::fingerprint() const noexcept {
  return storage_->fingerprint;
}

} // namespace hundun::immersed
