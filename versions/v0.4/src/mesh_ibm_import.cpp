// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kImportInput = 13801U;
constexpr std::uint32_t kImportCollective = 13802U;
constexpr std::uint32_t kImportRegion = 13803U;
constexpr std::uint32_t kImportLinks = 13804U;
constexpr std::uint32_t kImportDonors = 13805U;
constexpr std::uint32_t kImportMetric = 13806U;
constexpr std::uint32_t kImportAllocation = 13807U;
constexpr std::uint32_t kQuadraticCoverage = 1304U;
constexpr std::uint32_t kImportSchemaRevision = 1U;
constexpr std::uint8_t kRegionHalo = 4U;
constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

constexpr std::array<Int3, 6U> kDirectionDelta{{
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
    {0, 1, 0},  {0, 0, -1}, {0, 0, 1},
}};

static_assert(std::is_nothrow_move_assignable_v<EBTopology>);
static_assert(std::is_nothrow_move_assignable_v<BoundaryStencilPlan>);
static_assert(std::is_nothrow_move_assignable_v<SurfaceQuadraturePlan>);

class Hash64 {
 public:
  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    const Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t shift = 0U; shift < sizeof(bits) * 8U; shift += 8U) {
      value_ ^= static_cast<std::uint8_t>((bits >> shift) & Unsigned{0xffU});
      value_ *= kFnvPrime;
    }
  }

  void real(double value) noexcept {
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  PlanFingerprint finish() const noexcept { return value_ == 0U ? 1U : value_; }

 private:
  std::uint64_t value_{kFnvOffset};
};

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& out) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) return false;
  out = left + right;
  return true;
}

bool checked_multiply(std::uint64_t left, std::uint64_t right,
                      std::uint64_t& out) noexcept {
  if (left != 0U &&
      right > std::numeric_limits<std::uint64_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

bool checked_product(Int3 cells, std::uint64_t& out) noexcept {
  std::uint64_t xy{};
  return cells.x > 0 && cells.y > 0 && cells.z > 0 &&
         checked_multiply(static_cast<std::uint64_t>(cells.x),
                          static_cast<std::uint64_t>(cells.y), xy) &&
         checked_multiply(xy, static_cast<std::uint64_t>(cells.z), out);
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

void split_axis(std::int32_t global, std::int32_t partitions,
                std::int32_t coordinate, std::int32_t& begin,
                std::int32_t& cells) noexcept {
  const std::int32_t quotient = global / partitions;
  const std::int32_t remainder = global % partitions;
  cells = quotient + (coordinate < remainder ? 1 : 0);
  begin = coordinate * quotient + std::min(coordinate, remainder);
}

bool valid_patch(Int3 global, const MeshPatch& patch, int rank,
                 int size) noexcept {
  if (size <= 0 || patch.process_grid.x <= 0 ||
      patch.process_grid.y <= 0 || patch.process_grid.z <= 0 ||
      patch.process_grid.x > global.x || patch.process_grid.y > global.y ||
      patch.process_grid.z > global.z || patch.process_coord.x < 0 ||
      patch.process_coord.y < 0 || patch.process_coord.z < 0 ||
      patch.process_coord.x >= patch.process_grid.x ||
      patch.process_coord.y >= patch.process_grid.y ||
      patch.process_coord.z >= patch.process_grid.z) {
    return false;
  }
  std::uint64_t process_count{};
  if (!checked_product(patch.process_grid, process_count) ||
      process_count != static_cast<std::uint64_t>(size)) {
    return false;
  }
  const std::int64_t expected_rank =
      patch.process_coord.x +
      static_cast<std::int64_t>(patch.process_grid.x) *
          (patch.process_coord.y +
           static_cast<std::int64_t>(patch.process_grid.y) *
               patch.process_coord.z);
  if (expected_rank != rank) return false;
  Int3 expected_begin{};
  Int3 expected_cells{};
  split_axis(global.x, patch.process_grid.x, patch.process_coord.x,
             expected_begin.x, expected_cells.x);
  split_axis(global.y, patch.process_grid.y, patch.process_coord.y,
             expected_begin.y, expected_cells.y);
  split_axis(global.z, patch.process_grid.z, patch.process_coord.z,
             expected_begin.z, expected_cells.z);
  return same(patch.begin, expected_begin) && same(patch.cells, expected_cells);
}

bool valid_stencil_limits(QuadraticStencilLimits limits) noexcept {
  return limits.minimum_donors >= 14U &&
         limits.minimum_donors <= limits.maximum_donors &&
         limits.maximum_donors <= 32U && limits.maximum_reach > 0U &&
         limits.maximum_reach <= kRegionHalo &&
         limits.minimum_normal_bands >= 3U &&
         std::isfinite(limits.condition_limit) &&
         limits.condition_limit >= 1.0 && limits.condition_limit <= 1.0e8 &&
         limits.policy <= IbmReconstructionPolicy::adaptive_order &&
         limits.minimum_linear_donors >= 4U &&
         limits.minimum_linear_donors <= limits.maximum_donors &&
         limits.minimum_linear_normal_bands >= 2U &&
         limits.minimum_linear_normal_bands <= limits.minimum_normal_bands &&
         limits.standard_reach > 0U &&
         limits.standard_reach <= limits.maximum_reach;
}

bool inside(Int3 cell, Int3 global) noexcept {
  return cell.x >= 0 && cell.y >= 0 && cell.z >= 0 && cell.x < global.x &&
         cell.y < global.y && cell.z < global.z;
}

std::size_t flat(Int3 shape, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(cell.z));
}

GlobalCellId global_id(Int3 cell, Int3 global) noexcept {
  return static_cast<GlobalCellId>(cell.x) +
         static_cast<GlobalCellId>(global.x) *
             (static_cast<GlobalCellId>(cell.y) +
              static_cast<GlobalCellId>(global.y) *
                  static_cast<GlobalCellId>(cell.z));
}

std::uint64_t mix_cell(GlobalCellId id, std::uint8_t region) noexcept {
  std::uint64_t value = id ^ (static_cast<std::uint64_t>(region) << 63U) ^
                        UINT64_C(0x9e3779b97f4a7c15);
  value ^= value >> 30U;
  value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27U;
  value *= UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

double component(Real3 value, std::size_t axis) noexcept {
  return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
}

double centre(const CartesianGeometryPlan& geometry, Int3 cell,
              std::size_t axis) noexcept {
  const std::int32_t index =
      axis == 0U ? cell.x : (axis == 1U ? cell.y : cell.z);
  return geometry.axis(static_cast<CartesianAxis>(axis))
      .centres()
      .data[static_cast<std::size_t>(index)];
}

double width(const CartesianGeometryPlan& geometry, Int3 cell,
             std::size_t axis) noexcept {
  const std::int32_t index =
      axis == 0U ? cell.x : (axis == 1U ? cell.y : cell.z);
  return geometry.axis(static_cast<CartesianAxis>(axis))
      .widths()
      .data[static_cast<std::size_t>(index)];
}

Real3 cell_centre(const CartesianGeometryPlan& geometry, Int3 cell) noexcept {
  return {centre(geometry, cell, 0U), centre(geometry, cell, 1U),
          centre(geometry, cell, 2U)};
}

Real3 subtract(Real3 left, Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

double dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Real3 cross(Real3 left, Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double norm_squared(Real3 value) noexcept { return dot(value, value); }

Real3 tangent1_for(Real3 normal) noexcept {
  const Real3 reference = std::abs(normal.x) <= std::abs(normal.y) &&
                                  std::abs(normal.x) <= std::abs(normal.z)
                              ? Real3{1.0, 0.0, 0.0}
                          : std::abs(normal.y) <= std::abs(normal.z)
                              ? Real3{0.0, 1.0, 0.0}
                              : Real3{0.0, 0.0, 1.0};
  Real3 tangent = cross(normal, reference);
  const double inverse = 1.0 / std::sqrt(norm_squared(tangent));
  tangent.x *= inverse;
  tangent.y *= inverse;
  tangent.z *= inverse;
  return tangent;
}

Real3 link_normal(std::size_t direction) noexcept {
  const Int3 delta = kDirectionDelta[direction];
  return {-static_cast<double>(delta.x), -static_cast<double>(delta.y),
          -static_cast<double>(delta.z)};
}

Real3 wall_point(const CartesianGeometryPlan& geometry, Int3 fluid,
                 std::size_t direction) noexcept {
  Real3 point = cell_centre(geometry, fluid);
  const std::size_t axis = direction / 2U;
  const bool positive = (direction & 1U) != 0U;
  const std::int32_t index =
      axis == 0U ? fluid.x : (axis == 1U ? fluid.y : fluid.z);
  const double face = geometry.axis(static_cast<CartesianAxis>(axis))
                          .faces()
                          .data[static_cast<std::size_t>(index +
                                                         (positive ? 1 : 0))];
  if (axis == 0U)
    point.x = face;
  else if (axis == 1U)
    point.y = face;
  else
    point.z = face;
  return point;
}

double face_area(const CartesianGeometryPlan& geometry, Int3 fluid,
                 std::size_t direction) noexcept {
  const std::size_t axis = direction / 2U;
  if (axis == 0U) return width(geometry, fluid, 1U) * width(geometry, fluid, 2U);
  if (axis == 1U) return width(geometry, fluid, 0U) * width(geometry, fluid, 2U);
  return width(geometry, fluid, 0U) * width(geometry, fluid, 1U);
}

Int3 anchor_cell(const CartesianGeometryPlan& geometry, Real3 point) noexcept {
  const auto select = [](Span<const double> values, double coordinate) {
    const double* const found =
        std::lower_bound(values.data, values.data + values.size, coordinate);
    const std::ptrdiff_t raw = found - values.data;
    return static_cast<std::int32_t>(std::max<std::ptrdiff_t>(
        0, std::min<std::ptrdiff_t>(
               raw, static_cast<std::ptrdiff_t>(values.size - 1U))));
  };
  return {select(geometry.x().centres(), point.x),
          select(geometry.y().centres(), point.y),
          select(geometry.z().centres(), point.z)};
}

Status mpi_context(MPI_Comm communicator, int& rank, int& size) noexcept {
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0 ||
      MPI_Finalized(&finalized) != MPI_SUCCESS || finalized != 0 ||
      communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || rank < 0 ||
      size <= 0) {
    rank = -1;
    size = 0;
    return {StatusCode::mpi_failure, kImportCollective};
  }
  return {};
}

Status consensus(MPI_Comm communicator, int rank, int size, Status local,
                 int& lowest) noexcept {
  const int candidate = local ? size : rank;
  int selected = size;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kImportCollective};
  }
  if (selected == size) {
    lowest = -1;
    return {};
  }
  std::array<std::uint32_t, 2U> wire{};
  if (rank == selected) {
    wire = {static_cast<std::uint32_t>(local.code), local.detail};
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT32_T,
                selected, communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kImportCollective};
  }
  lowest = selected;
  return {static_cast<StatusCode>(wire[0]), wire[1]};
}

struct CandidateDonor {
  QuadraticDonorCell donor{};
  double distance_squared{};
  double normal_coordinate{};
  double tangent1_coordinate{};
  double tangent2_coordinate{};
};

std::uint8_t donor_quadrant(const CandidateDonor& donor) noexcept {
  const std::uint8_t high_t1 = donor.tangent1_coordinate >= 0.0 ? 1U : 0U;
  const std::uint8_t high_t2 = donor.tangent2_coordinate >= 0.0 ? 1U : 0U;
  return static_cast<std::uint8_t>(high_t1 * 2U + high_t2);
}

bool donor_less(const CandidateDonor& left,
                const CandidateDonor& right) noexcept {
  return left.distance_squared != right.distance_squared
             ? left.distance_squared < right.distance_squared
             : left.donor.global_cell < right.donor.global_cell;
}

bool same_coordinate_band(double left, double right) noexcept {
  return std::abs(left - right) <=
         512.0 * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(left), std::abs(right)});
}

std::uint8_t choose_donors(const std::vector<CandidateDonor>& candidates,
                           std::size_t maximum,
                           std::vector<QuadraticDonorCell>& out) {
  std::vector<std::uint8_t> selected(candidates.size(), 0U);
  out.clear();
  out.reserve(maximum);
  const auto take = [&](std::size_t index) {
    if (out.size() < maximum && selected[index] == 0U) {
      selected[index] = 1U;
      out.push_back(candidates[index].donor);
    }
  };
  for (std::uint8_t quadrant = 0U; quadrant < 4U; ++quadrant) {
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
      if (donor_quadrant(candidates[index]) == quadrant) {
        take(index);
        break;
      }
    }
  }
  std::vector<std::size_t> coordinate_order(candidates.size());
  // On stretched meshes, distance-only truncation can fill all 32 slots from
  // one coordinate plane even though valid off-plane donors are in reach.
  // Preserve the local polynomial dimensions before filling by proximity.
  const auto take_bands = [&](const auto& coordinate) {
    for (std::size_t index = 0U; index < coordinate_order.size(); ++index)
      coordinate_order[index] = index;
    std::sort(coordinate_order.begin(), coordinate_order.end(),
              [&](std::size_t left, std::size_t right) {
                const double a = coordinate(candidates[left]);
                const double b = coordinate(candidates[right]);
                if (std::abs(a) != std::abs(b))
                  return std::abs(a) < std::abs(b);
                if (a != b) return a < b;
                return donor_less(candidates[left], candidates[right]);
              });
    std::array<double, 3U> values{};
    std::uint8_t bands = 0U;
    for (const std::size_t index : coordinate_order) {
      const double value = coordinate(candidates[index]);
      bool distinct = true;
      for (std::size_t band = 0U; band < bands; ++band) {
        distinct &= !same_coordinate_band(values[band], value);
      }
      if (!distinct) continue;
      take(index);
      values[bands++] = value;
      if (bands == values.size()) break;
    }
    return bands;
  };
  const std::uint8_t bands = take_bands(
      [](const CandidateDonor& value) { return value.normal_coordinate; });
  take_bands(
      [](const CandidateDonor& value) { return value.tangent1_coordinate; });
  take_bands(
      [](const CandidateDonor& value) { return value.tangent2_coordinate; });
  for (std::size_t index = 0U;
       index < candidates.size() && out.size() < maximum; ++index) {
    take(index);
  }
  std::sort(out.begin(), out.end(), [](const QuadraticDonorCell& left,
                                       const QuadraticDonorCell& right) {
    return left.global_cell < right.global_cell;
  });
  return bands;
}

Status collect_donors(const CartesianGeometryPlan& geometry,
                      const MeshPatch& patch,
                      Span<const std::uint8_t> marker, Real3 origin,
                      Real3 normal, QuadraticStencilLimits limits,
                      QuadraticFrame& frame,
                      std::vector<CandidateDonor>& candidates,
                      std::vector<QuadraticDonorCell>& out,
                      std::uint8_t& normal_bands) {
  frame.origin = origin;
  frame.normal = normal;
  frame.tangent1 = tangent1_for(normal);
  frame.tangent2 = cross(normal, frame.tangent1);
  frame.anchor_global_cell = anchor_cell(geometry, origin);
  frame.scale = std::cbrt(width(geometry, frame.anchor_global_cell, 0U) *
                          width(geometry, frame.anchor_global_cell, 1U) *
                          width(geometry, frame.anchor_global_cell, 2U));
  const Int3 global = geometry.global_cells();
  const int reach = static_cast<int>(limits.maximum_reach);
  candidates.clear();
  std::uint8_t available_quadrants = 0U;
  for (int dz = -reach; dz <= reach; ++dz) {
    for (int dy = -reach; dy <= reach; ++dy) {
      for (int dx = -reach; dx <= reach; ++dx) {
        const Int3 index{frame.anchor_global_cell.x + dx,
                         frame.anchor_global_cell.y + dy,
                         frame.anchor_global_cell.z + dz};
        if (!inside(index, global)) continue;
        const GlobalCellId id = global_id(index, global);
        if (marker.data[static_cast<std::size_t>(id)] !=
            static_cast<std::uint8_t>(RegionFlag::fluid)) {
          continue;
        }
        const Real3 position = cell_centre(geometry, index);
        const Real3 delta = subtract(position, origin);
        const double normal_coordinate = dot(delta, normal) / frame.scale;
        if (!(normal_coordinate >
              64.0 * std::numeric_limits<double>::epsilon())) {
          continue;
        }
        QuadraticDonorCell donor;
        donor.global_cell = id;
        donor.global_index = index;
        donor.local_index = {index.x - patch.begin.x,
                             index.y - patch.begin.y,
                             index.z - patch.begin.z};
        donor.centre = position;
        donor.widths = {width(geometry, index, 0U),
                        width(geometry, index, 1U),
                        width(geometry, index, 2U)};
        CandidateDonor candidate{
            donor, norm_squared(delta), normal_coordinate,
            dot(delta, frame.tangent1) / frame.scale,
            dot(delta, frame.tangent2) / frame.scale};
        available_quadrants = static_cast<std::uint8_t>(
            available_quadrants |
            (UINT8_C(1) << donor_quadrant(candidate)));
        candidates.push_back(candidate);
      }
    }
  }
  if (available_quadrants == 0U) {
    return {StatusCode::invalid_plan, kImportDonors};
  }
  frame.required_quadrant_mask = available_quadrants;
  std::sort(candidates.begin(), candidates.end(), donor_less);
  const std::size_t minimum =
      limits.policy == IbmReconstructionPolicy::adaptive_order
          ? std::min(limits.minimum_donors, limits.minimum_linear_donors)
          : limits.minimum_donors;
  if (candidates.size() < minimum) {
    return {StatusCode::invalid_plan, kImportDonors};
  }
  normal_bands = choose_donors(
      candidates,
      std::min<std::size_t>(candidates.size(), limits.maximum_donors), out);
  return {};
}

PlanFingerprint marker_fingerprint(const CartesianGeometryPlan& geometry,
                                   PlanFingerprint marker_source,
                                   Span<const std::uint8_t> marker) noexcept {
  Hash64 hash;
  hash.integer(kImportSchemaRevision);
  hash.integer(geometry.fingerprint());
  hash.integer(marker_source);
  hash.integer(static_cast<std::uint64_t>(marker.size));
  for (std::size_t index = 0U; index < marker.size; ++index)
    hash.integer(marker.data[index]);
  return hash.finish();
}

PlanFingerprint compile_contract(const CartesianGeometryPlan& geometry,
                                 const MeshPatch& patch,
                                 PlanFingerprint surface_fingerprint,
                                 ImmersedPlanLimits limits, int size) noexcept {
  Hash64 hash;
  hash.integer(kImportSchemaRevision);
  hash.integer(geometry.fingerprint());
  hash.integer(geometry.topology_revision());
  hash.integer(surface_fingerprint);
  hash.integer(patch.process_grid.x);
  hash.integer(patch.process_grid.y);
  hash.integer(patch.process_grid.z);
  hash.integer(size);
  hash.integer(limits.stencil.minimum_donors);
  hash.integer(limits.stencil.maximum_donors);
  hash.integer(limits.stencil.maximum_reach);
  hash.integer(limits.stencil.minimum_normal_bands);
  hash.real(limits.stencil.condition_limit);
  hash.integer(static_cast<std::uint8_t>(limits.stencil.policy));
  hash.integer(limits.stencil.minimum_linear_donors);
  hash.integer(limits.stencil.minimum_linear_normal_bands);
  hash.integer(limits.stencil.standard_reach);
  hash.integer(limits.maximum_persistent_bytes_per_rank);
  hash.integer(limits.maximum_peak_bytes_per_rank);
  hash.integer(limits.maximum_local_links);
  hash.integer(limits.maximum_local_quadrature_points);
  return hash.finish();
}

}  // namespace

Status ImportedIbmCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, Span<const std::uint8_t> global_marker,
    PlanFingerprint marker_source, ImmersedPlanLimits limits,
    EBTopology& topology, BoundaryStencilPlan& boundary,
    SurfaceQuadraturePlan& quadrature) noexcept {
  int rank = -1;
  int size = 0;
  const Status context = mpi_context(communicator, rank, size);
  if (!context) return context;

  const Int3 global = geometry.global_cells();
  std::uint64_t global_count{};
  std::uint64_t local_count{};
  Status local{};
  if (geometry.fingerprint() == 0U || geometry.topology_revision() == 0U ||
      marker_source == 0U || !checked_product(global, global_count) ||
      !checked_product(patch.cells, local_count) ||
      global_count > std::numeric_limits<std::size_t>::max() ||
      global_count > (UINT64_MAX - 5U) / 6U ||
      global_marker.data == nullptr || global_marker.size != global_count ||
      local_count > UINT32_MAX || !valid_patch(global, patch, rank, size) ||
      !valid_stencil_limits(limits.stencil) ||
      limits.maximum_persistent_bytes_per_rank == 0U ||
      limits.maximum_peak_bytes_per_rank <
          limits.maximum_persistent_bytes_per_rank ||
      limits.maximum_local_links == 0U ||
      limits.maximum_local_links > UINT32_MAX ||
      limits.maximum_local_quadrature_points == 0U) {
    local = {StatusCode::invalid_plan, kImportInput};
  }
  if (local) {
    for (std::size_t index = 0U; index < global_marker.size; ++index) {
      if (global_marker.data[index] !=
              static_cast<std::uint8_t>(RegionFlag::solid) &&
          global_marker.data[index] !=
              static_cast<std::uint8_t>(RegionFlag::fluid)) {
        local = {StatusCode::invalid_plan, kImportRegion};
        break;
      }
    }
  }
  int lowest = -1;
  Status agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  const PlanFingerprint imported_surface =
      marker_fingerprint(geometry, marker_source, global_marker);
  const PlanFingerprint contract =
      compile_contract(geometry, patch, imported_surface, limits, size);
  PlanFingerprint root_contract = contract;
  if (MPI_Bcast(&root_contract, 1, MPI_UINT64_T, 0, communicator) !=
      MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kImportCollective};
  } else if (root_contract != contract) {
    local = {StatusCode::invalid_plan, kImportInput};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  EBTopology topology_candidate;
  BoundaryStencilPlan boundary_candidate;
  SurfaceQuadraturePlan quadrature_candidate;
  const std::int64_t halo_x =
      static_cast<std::int64_t>(patch.cells.x) + 2 * kRegionHalo;
  const std::int64_t halo_y =
      static_cast<std::int64_t>(patch.cells.y) + 2 * kRegionHalo;
  const std::int64_t halo_z =
      static_cast<std::int64_t>(patch.cells.z) + 2 * kRegionHalo;
  Int3 halo_shape{};
  std::uint64_t halo_count{};
  if (halo_x > std::numeric_limits<std::int32_t>::max() ||
      halo_y > std::numeric_limits<std::int32_t>::max() ||
      halo_z > std::numeric_limits<std::int32_t>::max()) {
    local = {StatusCode::invalid_plan, kImportRegion};
  } else {
    halo_shape = {static_cast<std::int32_t>(halo_x),
                  static_cast<std::int32_t>(halo_y),
                  static_cast<std::int32_t>(halo_z)};
    if (!checked_product(halo_shape, halo_count) ||
        halo_count > std::numeric_limits<std::size_t>::max()) {
      local = {StatusCode::invalid_plan, kImportRegion};
    }
  }
  std::uint64_t base_persistent{};
  if (local &&
      (!checked_add(local_count, halo_count, base_persistent) ||
       base_persistent > limits.maximum_persistent_bytes_per_rank)) {
    local = {StatusCode::invalid_plan, kImportRegion};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  std::uint64_t local_xor = 0U;
  std::uint64_t local_sum = 0U;
  try {
    topology_candidate.region_.resize(static_cast<std::size_t>(local_count));
    topology_candidate.halo_region_.assign(
        static_cast<std::size_t>(halo_count),
        static_cast<std::uint8_t>(RegionFlag::solid));
    topology_candidate.halo_stride_y_ = static_cast<std::size_t>(halo_shape.x);
    topology_candidate.halo_stride_z_ =
        topology_candidate.halo_stride_y_ *
        static_cast<std::size_t>(halo_shape.y);
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        for (std::int32_t x = 0; x < patch.cells.x; ++x) {
          const Int3 owned{x, y, z};
          const Int3 global_cell{patch.begin.x + x, patch.begin.y + y,
                                 patch.begin.z + z};
          const std::uint8_t region = global_marker.data[static_cast<std::size_t>(
              global_id(global_cell, global))];
          topology_candidate.region_[flat(patch.cells, owned)] = region;
          local_xor ^= mix_cell(global_id(global_cell, global), region);
          local_sum += mix_cell(global_id(global_cell, global), region);
        }
      }
    }
    for (std::int32_t z = 0; z < halo_shape.z; ++z) {
      for (std::int32_t y = 0; y < halo_shape.y; ++y) {
        for (std::int32_t x = 0; x < halo_shape.x; ++x) {
          const Int3 global_cell{patch.begin.x + x - kRegionHalo,
                                 patch.begin.y + y - kRegionHalo,
                                 patch.begin.z + z - kRegionHalo};
          if (inside(global_cell, global)) {
            topology_candidate.halo_region_[flat(halo_shape, {x, y, z})] =
                global_marker.data[static_cast<std::size_t>(
                    global_id(global_cell, global))];
          }
        }
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kImportAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kImportRegion};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  std::uint64_t local_link_count = 0U;
  std::uint64_t local_interface_count = 0U;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 owned{x, y, z};
        if (topology_candidate.region_[flat(patch.cells, owned)] !=
            static_cast<std::uint8_t>(RegionFlag::fluid)) {
          continue;
        }
        const Int3 fluid{patch.begin.x + x, patch.begin.y + y,
                         patch.begin.z + z};
        bool interface_cell = false;
        for (const Int3 delta : kDirectionDelta) {
          const Int3 solid{fluid.x + delta.x, fluid.y + delta.y,
                           fluid.z + delta.z};
          if (inside(solid, global) &&
              global_marker.data[static_cast<std::size_t>(
                  global_id(solid, global))] ==
                  static_cast<std::uint8_t>(RegionFlag::solid)) {
            interface_cell = true;
            ++local_link_count;
          }
        }
        if (interface_cell) ++local_interface_count;
      }
    }
  }
  std::uint64_t global_link_count = 0U;
  if (MPI_Allreduce(&local_link_count, &global_link_count, 1, MPI_UINT64_T,
                    MPI_SUM, communicator) != MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kImportCollective};
  }
  std::uint64_t link_bytes{};
  std::uint64_t interface_bytes{};
  std::uint64_t metric_bytes{};
  std::uint64_t topology_persistent{};
  if (local &&
      (global_link_count == 0U ||
       local_link_count > limits.maximum_local_links ||
       local_link_count > limits.maximum_local_quadrature_points ||
       local_link_count > UINT32_MAX || local_interface_count > UINT32_MAX ||
       !checked_multiply(local_link_count, sizeof(ImmersedLink), link_bytes) ||
       !checked_multiply(local_interface_count, sizeof(std::uint32_t),
                         interface_bytes) ||
       !checked_multiply(local_link_count, sizeof(IbmInterfaceLinkMetric),
                         metric_bytes) ||
       !checked_add(base_persistent, link_bytes, topology_persistent) ||
       !checked_add(topology_persistent, interface_bytes,
                    topology_persistent) ||
       !checked_add(topology_persistent, metric_bytes,
                    topology_persistent) ||
       topology_persistent > limits.maximum_persistent_bytes_per_rank)) {
    local = {StatusCode::invalid_plan, kImportLinks};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  try {
    topology_candidate.links_.reserve(static_cast<std::size_t>(local_link_count));
    topology_candidate.interface_cells_.reserve(
        static_cast<std::size_t>(local_interface_count));
    topology_candidate.interface_metric_.links_.reserve(
        static_cast<std::size_t>(local_link_count));
    boundary_candidate.links_.reserve(static_cast<std::size_t>(local_link_count));
    quadrature_candidate.local_points_.reserve(
        static_cast<std::size_t>(local_link_count));
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        for (std::int32_t x = 0; x < patch.cells.x; ++x) {
          const Int3 fluid_local{x, y, z};
          if (topology_candidate.region_[flat(patch.cells, fluid_local)] !=
              static_cast<std::uint8_t>(RegionFlag::fluid)) {
            continue;
          }
          const Int3 fluid{patch.begin.x + x, patch.begin.y + y,
                           patch.begin.z + z};
          bool interface_cell = false;
          for (std::size_t direction = 0U; direction < kDirectionDelta.size();
               ++direction) {
            const Int3 delta = kDirectionDelta[direction];
            const Int3 solid{fluid.x + delta.x, fluid.y + delta.y,
                             fluid.z + delta.z};
            if (!inside(solid, global) ||
                global_marker.data[static_cast<std::size_t>(
                    global_id(solid, global))] !=
                    static_cast<std::uint8_t>(RegionFlag::solid)) {
              continue;
            }
            interface_cell = true;
            ImmersedLink link;
            link.global_link = 6U * global_id(fluid, global) + direction;
            link.fluid_cell = global_id(fluid, global);
            link.solid_cell = global_id(solid, global);
            link.fluid_global_index = fluid;
            link.solid_global_index = solid;
            link.fluid_local_index = fluid_local;
            link.solid_local_index = {solid.x - patch.begin.x,
                                      solid.y - patch.begin.y,
                                      solid.z - patch.begin.z};
            link.direction = static_cast<ImmersedFaceDirection>(direction);
            link.triangle = link.global_link;
            link.wall_point = wall_point(geometry, fluid, direction);
            link.solid_to_fluid_normal = link_normal(direction);
            link.cartesian_control_face_area =
                face_area(geometry, fluid, direction);
            link.surface_patch_centroid = link.wall_point;
            topology_candidate.links_.push_back(link);

            IbmInterfaceLinkMetric metric;
            metric.global_link = link.global_link;
            metric.source_triangle = link.triangle;
            metric.physical_quadrature_area =
                link.cartesian_control_face_area;
            metric.physical_area_vector =
                {link.cartesian_control_face_area * link.solid_to_fluid_normal.x,
                 link.cartesian_control_face_area * link.solid_to_fluid_normal.y,
                 link.cartesian_control_face_area * link.solid_to_fluid_normal.z};
            metric.physical_first_moment =
                {link.cartesian_control_face_area * link.wall_point.x,
                 link.cartesian_control_face_area * link.wall_point.y,
                 link.cartesian_control_face_area * link.wall_point.z};
            for (std::size_t row = 0U; row < 3U; ++row) {
              for (std::size_t column = 0U; column < 3U; ++column) {
                metric.normal_first_moment[3U * row + column] =
                    link.cartesian_control_face_area *
                    component(link.wall_point, row) *
                    component(link.solid_to_fluid_normal, column);
                metric.normal_second_moment[3U * row + column] =
                    link.cartesian_control_face_area *
                    component(link.solid_to_fluid_normal, row) *
                    component(link.solid_to_fluid_normal, column);
              }
            }
            topology_candidate.interface_metric_.links_.push_back(metric);
          }
          if (interface_cell) {
            topology_candidate.interface_cells_.push_back(
                static_cast<std::uint32_t>(flat(patch.cells, fluid_local)));
          }
        }
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kImportAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kImportLinks};
  }
  if (local &&
      (topology_candidate.links_.size() != local_link_count ||
       topology_candidate.interface_cells_.size() != local_interface_count ||
       topology_candidate.interface_metric_.links_.size() !=
           local_link_count)) {
    local = {StatusCode::invalid_plan, kImportLinks};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  std::uint64_t global_xor{};
  std::uint64_t global_sum{};
  if (MPI_Allreduce(&local_xor, &global_xor, 1, MPI_UINT64_T, MPI_BXOR,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_UINT64_T, MPI_SUM,
                    communicator) != MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kImportCollective};
  }

  constexpr std::size_t kMetricChannels = 26U;
  std::array<double, kMetricChannels> local_metric{};
  for (std::size_t index = 0U; index < topology_candidate.links_.size();
       ++index) {
    const ImmersedLink& link = topology_candidate.links_[index];
    const IbmInterfaceLinkMetric& metric =
        topology_candidate.interface_metric_.links_[index];
    local_metric[0U] += link.cartesian_control_face_area;
    local_metric[1U] += metric.physical_quadrature_area;
    local_metric[2U] += metric.physical_area_vector.x;
    local_metric[3U] += metric.physical_area_vector.y;
    local_metric[4U] += metric.physical_area_vector.z;
    local_metric[5U] += metric.physical_first_moment.x;
    local_metric[6U] += metric.physical_first_moment.y;
    local_metric[7U] += metric.physical_first_moment.z;
    for (std::size_t entry = 0U; entry < 9U; ++entry) {
      local_metric[8U + entry] += metric.normal_first_moment[entry];
      local_metric[17U + entry] += metric.normal_second_moment[entry];
    }
  }
  std::array<double, kMetricChannels> global_metric{};
  if (MPI_Allreduce(local_metric.data(), global_metric.data(),
                    static_cast<int>(global_metric.size()), MPI_DOUBLE,
                    MPI_SUM, communicator) != MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kImportCollective};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  topology_candidate.global_cells_ = global;
  topology_candidate.patch_ = patch;
  topology_candidate.fluid_side_ = ImmersedFluidSide::outside;
  topology_candidate.region_halo_width_ = kRegionHalo;
  topology_candidate.geometry_revision_ = geometry.topology_revision();
  topology_candidate.geometry_fingerprint_ = geometry.fingerprint();
  topology_candidate.surface_fingerprint_ = imported_surface;
  Hash64 topology_hash;
  topology_hash.integer(kImportSchemaRevision);
  topology_hash.integer(geometry.fingerprint());
  topology_hash.integer(imported_surface);
  topology_hash.integer(global_count);
  topology_hash.integer(global_xor);
  topology_hash.integer(global_sum);
  topology_candidate.fingerprint_ = topology_hash.finish();

  IbmInterfaceMetricPlan& metric_plan = topology_candidate.interface_metric_;
  metric_plan.conservation_.cartesian_control_area = global_metric[0U];
  metric_plan.conservation_.physical_quadrature_area = global_metric[1U];
  metric_plan.conservation_.physical_area_vector =
      {global_metric[2U], global_metric[3U], global_metric[4U]};
  metric_plan.conservation_.physical_first_moment =
      {global_metric[5U], global_metric[6U], global_metric[7U]};
  for (std::size_t entry = 0U; entry < 9U; ++entry) {
    metric_plan.conservation_.normal_first_moment[entry] =
        global_metric[8U + entry];
    metric_plan.conservation_.normal_second_moment[entry] =
        global_metric[17U + entry];
  }
  metric_plan.resources_.persistent_bytes_per_rank = metric_bytes;
  metric_plan.resources_.peak_bytes_per_rank = metric_bytes;
  metric_plan.resources_.collective_doubles_per_rank = kMetricChannels;
  metric_plan.geometry_revision_ = geometry.topology_revision();
  metric_plan.geometry_fingerprint_ = geometry.fingerprint();
  metric_plan.surface_fingerprint_ = imported_surface;
  Hash64 physical_hash;
  physical_hash.integer(kImportSchemaRevision);
  physical_hash.integer(geometry.fingerprint());
  physical_hash.integer(imported_surface);
  physical_hash.integer(global_link_count);
  metric_plan.physical_fingerprint_ = physical_hash.finish();
  Hash64 metric_layout;
  metric_layout.integer(metric_plan.physical_fingerprint_);
  metric_layout.integer(rank);
  metric_layout.integer(static_cast<std::uint64_t>(metric_plan.links_.size()));
  for (const IbmInterfaceLinkMetric& metric : metric_plan.links_) {
    metric_layout.integer(metric.global_link);
    metric_layout.integer(metric.source_triangle);
    metric_layout.real(metric.physical_quadrature_area);
    metric_layout.real(metric.physical_area_vector.x);
    metric_layout.real(metric.physical_area_vector.y);
    metric_layout.real(metric.physical_area_vector.z);
  }
  metric_plan.fingerprint_ = metric_layout.finish();
  Hash64 sealed_topology;
  sealed_topology.integer(topology_candidate.fingerprint_);
  sealed_topology.integer(metric_plan.physical_fingerprint_);
  topology_candidate.fingerprint_ = sealed_topology.finish();

  try {
    struct Prepared {
      QuadraticFrame frame{};
      std::size_t donor_begin{};
      std::size_t donor_count{};
      std::size_t row_begin{};
      std::uint32_t regular_group{kInvalidIbmIndex};
      bool single_normal_band{};
    };
    std::vector<Prepared> prepared;
    std::vector<QuadraticDonorCell> donor_storage;
    std::vector<QuadraticFunctionalRequest> row_storage;
    std::vector<CandidateDonor> candidates;
    std::vector<QuadraticDonorCell> donors;
    prepared.reserve(topology_candidate.links_.size());
    donor_storage.reserve(topology_candidate.links_.size() *
                          limits.stencil.maximum_donors);
    row_storage.reserve(topology_candidate.links_.size() * 4U);
    for (const ImmersedLink& link : topology_candidate.links_) {
      QuadraticFrame frame;
      std::uint8_t normal_bands = 0U;
      local = collect_donors(geometry, patch, global_marker, link.wall_point,
                             link.solid_to_fluid_normal, limits.stencil, frame,
                             candidates, donors, normal_bands);
      if (!local) break;
      const Real3 ghost = cell_centre(geometry, link.solid_global_index);
      const std::array<QuadraticFunctionalRequest, 4U> rows{{
          {QuadraticFunctionalKind::value,
           QuadraticConstraint::origin_value, ghost, {}},
          {QuadraticFunctionalKind::value,
           QuadraticConstraint::origin_normal_gradient, ghost, {}},
          {QuadraticFunctionalKind::value, QuadraticConstraint::none,
           link.wall_point, {}},
          {QuadraticFunctionalKind::directional_derivative,
           QuadraticConstraint::origin_value, link.wall_point,
           link.solid_to_fluid_normal},
      }};
      const std::size_t donor_begin = donor_storage.size();
      const std::size_t row_begin = row_storage.size();
      donor_storage.insert(donor_storage.end(), donors.begin(), donors.end());
      row_storage.insert(row_storage.end(), rows.begin(), rows.end());
      prepared.push_back({frame, donor_begin, donors.size(), row_begin,
                          kInvalidIbmIndex, normal_bands == 1U});
    }
    std::vector<QuadraticStencilRequest> regular_requests;
    QuadraticStencilPlan regular_plan;
    std::uint64_t thin_link_count = 0U;
    if (local) regular_requests.reserve(prepared.size());
    for (std::size_t index = 0U; local && index < prepared.size(); ++index) {
      Prepared& item = prepared[index];
      const QuadraticStencilRequest request{
          item.frame,
          {donor_storage.data() + item.donor_begin, item.donor_count},
          {row_storage.data() + item.row_begin, 4U}};
      if (item.single_normal_band &&
          limits.stencil.policy == IbmReconstructionPolicy::adaptive_order) {
        // A one-cell-thick Cartesian fluid sheet has no second positive-normal
        // sample at any configured reach.  First prove that the ordinary
        // compiler rejects exactly this mathematical coverage defect; no
        // other failure is eligible for the explicit two-point closure.
        QuadraticStencilPlan rejected;
        const Status ordinary = QuadraticStencilCompiler::compile(
            {&request, 1U}, limits.stencil, rejected);
        if (ordinary || ordinary.code != StatusCode::invalid_plan ||
            ordinary.detail != kQuadraticCoverage) {
          local = ordinary ? Status{StatusCode::invalid_plan, kImportDonors}
                           : ordinary;
          break;
        }
        ++thin_link_count;
        continue;
      }
      if (regular_requests.size() >= UINT32_MAX) {
        local = {StatusCode::invalid_plan, kImportDonors};
        break;
      }
      item.regular_group =
          static_cast<std::uint32_t>(regular_requests.size());
      regular_requests.push_back(request);
    }
    if (local && !regular_requests.empty()) {
      local = QuadraticStencilCompiler::compile(
          {regular_requests.data(), regular_requests.size()}, limits.stencil,
          regular_plan);
    }
    if (local) {
      QuadraticStencilPlan& merged = boundary_candidate.reconstruction_;
      const auto append_regular = [&](std::uint32_t source_group) -> Status {
        if (source_group >= regular_plan.groups_.size())
          return {StatusCode::invalid_plan, kImportDonors};
        const QuadraticStencilGroup& source =
            regular_plan.groups_[source_group];
        const std::size_t donor_count = source.quality.donor_count;
        if (source.row_count != 4U ||
            source.donor_begin > regular_plan.donor_global_cells_.size() ||
            donor_count > regular_plan.donor_global_cells_.size() -
                              source.donor_begin ||
            source.donor_begin > regular_plan.donor_local_indices_.size() ||
            donor_count > regular_plan.donor_local_indices_.size() -
                              source.donor_begin ||
            source.row_begin > regular_plan.rows_.size() ||
            source.row_count > regular_plan.rows_.size() - source.row_begin ||
            merged.groups_.size() > UINT32_MAX ||
            merged.rows_.size() > UINT32_MAX - 4U ||
            merged.donor_global_cells_.size() > UINT32_MAX - donor_count) {
          return {StatusCode::invalid_plan, kImportDonors};
        }
        QuadraticStencilGroup group = source;
        group.donor_begin = static_cast<std::uint32_t>(
            merged.donor_global_cells_.size());
        group.row_begin = static_cast<std::uint32_t>(merged.rows_.size());
        merged.donor_global_cells_.insert(
            merged.donor_global_cells_.end(),
            regular_plan.donor_global_cells_.begin() + source.donor_begin,
            regular_plan.donor_global_cells_.begin() + source.donor_begin +
                donor_count);
        merged.donor_local_indices_.insert(
            merged.donor_local_indices_.end(),
            regular_plan.donor_local_indices_.begin() + source.donor_begin,
            regular_plan.donor_local_indices_.begin() + source.donor_begin +
                donor_count);
        for (std::size_t row_offset = 0U; row_offset < 4U; ++row_offset) {
          const QuadraticAffineRow& source_row =
              regular_plan.rows_[source.row_begin + row_offset];
          if (source_row.weight_begin > regular_plan.weights_.size() ||
              donor_count >
                  regular_plan.weights_.size() - source_row.weight_begin ||
              merged.weights_.size() > UINT32_MAX - donor_count) {
            return {StatusCode::invalid_plan, kImportDonors};
          }
          QuadraticAffineRow row = source_row;
          row.group = static_cast<std::uint32_t>(merged.groups_.size());
          row.weight_begin =
              static_cast<std::uint32_t>(merged.weights_.size());
          merged.weights_.insert(
              merged.weights_.end(),
              regular_plan.weights_.begin() + source_row.weight_begin,
              regular_plan.weights_.begin() + source_row.weight_begin +
                  donor_count);
          merged.rows_.push_back(row);
        }
        merged.groups_.push_back(group);
        merged.maximum_halo_reach_ =
            std::max(merged.maximum_halo_reach_, source.quality.reach);
        return {};
      };
      const auto append_two_point = [&](std::size_t index) -> Status {
        const Prepared& item = prepared[index];
        const ImmersedLink& link = topology_candidate.links_[index];
        const Real3 fluid = cell_centre(geometry, link.fluid_global_index);
        const Real3 solid = cell_centre(geometry, link.solid_global_index);
        const double fluid_distance =
            std::sqrt(norm_squared(subtract(fluid, link.wall_point)));
        const double solid_distance =
            std::sqrt(norm_squared(subtract(solid, link.wall_point)));
        if (!std::isfinite(fluid_distance) || !(fluid_distance > 0.0) ||
            !std::isfinite(solid_distance) || !(solid_distance > 0.0) ||
            merged.groups_.size() > UINT32_MAX ||
            merged.rows_.size() > UINT32_MAX - 4U ||
            merged.donor_global_cells_.size() > UINT32_MAX - 1U ||
            merged.weights_.size() > UINT32_MAX - 4U) {
          return {StatusCode::invalid_plan, kImportDonors};
        }
        const double ratio = solid_distance / fluid_distance;
        const std::array<double, 4U> weights{{
            -ratio, 1.0, 1.0, 1.0 / fluid_distance,
        }};
        const std::uint32_t group_index =
            static_cast<std::uint32_t>(merged.groups_.size());
        const std::uint32_t donor_begin =
            static_cast<std::uint32_t>(merged.donor_global_cells_.size());
        const std::uint32_t row_begin =
            static_cast<std::uint32_t>(merged.rows_.size());
        merged.donor_global_cells_.push_back(link.fluid_cell);
        merged.donor_local_indices_.push_back(link.fluid_local_index);
        merged.rows_.push_back(
            {group_index, static_cast<std::uint32_t>(merged.weights_.size()),
             1.0 + ratio, 0.0});
        merged.weights_.push_back(weights[0U]);
        merged.rows_.push_back(
            {group_index, static_cast<std::uint32_t>(merged.weights_.size()),
             0.0, -(fluid_distance + solid_distance)});
        merged.weights_.push_back(weights[1U]);
        merged.rows_.push_back(
            {group_index, static_cast<std::uint32_t>(merged.weights_.size()),
             0.0, 0.0});
        merged.weights_.push_back(weights[2U]);
        merged.rows_.push_back(
            {group_index, static_cast<std::uint32_t>(merged.weights_.size()),
             -1.0 / fluid_distance, 0.0});
        merged.weights_.push_back(weights[3U]);

        const auto absolute_row_sum = [](double donor, double wall,
                                         double gradient) noexcept {
          return std::abs(donor) + std::abs(wall) + std::abs(gradient);
        };
        const double functional_l1 =
            std::max({absolute_row_sum(weights[0U], 1.0 + ratio, 0.0),
                      absolute_row_sum(
                          weights[1U], 0.0,
                          -(fluid_distance + solid_distance)),
                      absolute_row_sum(weights[2U], 0.0, 0.0),
                      absolute_row_sum(weights[3U],
                                       -1.0 / fluid_distance, 0.0)});
        Hash64 pivot;
        pivot.integer(UINT64_C(0x74776f706f696e74));
        pivot.integer(link.global_link);
        pivot.integer(link.fluid_cell);
        pivot.real(fluid_distance);
        pivot.real(solid_distance);
        QuadraticStencilGroup group;
        group.donor_begin = donor_begin;
        group.row_begin = row_begin;
        group.row_count = 4U;
        group.quality.donor_count = 1U;
        group.quality.normal_band_count = 1U;
        group.quality.quadrant_mask = 0x08U;
        group.quality.required_quadrant_mask = 0x08U;
        const std::uint64_t dx = static_cast<std::uint64_t>(std::abs(
            link.fluid_global_index.x - item.frame.anchor_global_cell.x));
        const std::uint64_t dy = static_cast<std::uint64_t>(std::abs(
            link.fluid_global_index.y - item.frame.anchor_global_cell.y));
        const std::uint64_t dz = static_cast<std::uint64_t>(std::abs(
            link.fluid_global_index.z - item.frame.anchor_global_cell.z));
        const std::uint64_t reach = std::max({dx, dy, dz});
        if (reach > UINT8_MAX) {
          return {StatusCode::invalid_plan, kImportDonors};
        }
        group.quality.reach = static_cast<std::uint8_t>(reach);
        group.quality.rank = 1U;
        group.quality.condition_estimate = 1.0;
        group.quality.functional_l1 = functional_l1;
        group.quality.pivot_fingerprint = pivot.finish();
        group.quality.order = IbmReconstructionOrder::linear;
        group.quality.fallback_reason =
            IbmReconstructionFallbackReason::quadratic_coverage;
        Hash64 fingerprint;
        fingerprint.integer(UINT64_C(0x696d706f72743270));
        fingerprint.integer(link.global_link);
        fingerprint.integer(link.fluid_cell);
        fingerprint.integer(group.quality.pivot_fingerprint);
        fingerprint.real(fluid_distance);
        fingerprint.real(solid_distance);
        for (const double weight : weights) fingerprint.real(weight);
        group.fingerprint = fingerprint.finish();
        merged.groups_.push_back(group);
        merged.maximum_halo_reach_ = std::max(
            merged.maximum_halo_reach_, group.quality.reach);
        merged.audit_.maximum_condition_estimate = std::max(
            merged.audit_.maximum_condition_estimate, 1.0);
        merged.audit_.maximum_functional_l1 = std::max(
            merged.audit_.maximum_functional_l1, functional_l1);
        return {};
      };

      merged.audit_ = regular_plan.audit_;
      merged.audit_.valid = true;
      merged.audit_.policy = limits.stencil.policy;
      merged.audit_.standard_reach = limits.stencil.standard_reach;
      merged.audit_.linear_groups += thin_link_count;
      merged.audit_.coverage_fallback_groups += thin_link_count;
      for (std::size_t index = 0U; local && index < prepared.size(); ++index) {
        if (prepared[index].regular_group != kInvalidIbmIndex) {
          local = append_regular(prepared[index].regular_group);
        } else {
          local = append_two_point(index);
        }
      }
      merged.audit_.group_count = merged.groups_.size();
      if (local && (merged.groups_.size() != prepared.size() ||
                    merged.audit_.group_count !=
                        merged.audit_.quadratic_groups +
                            merged.audit_.linear_groups)) {
        local = {StatusCode::invalid_plan, kImportDonors};
      }
      if (local) merged.refresh_fingerprint();
    }
    if (local) {
      const Span<const QuadraticStencilGroup> groups =
          boundary_candidate.reconstruction_.groups();
      if (groups.size != topology_candidate.links_.size()) {
        local = {StatusCode::invalid_plan, kImportDonors};
      } else {
        for (std::size_t index = 0U; index < groups.size; ++index) {
          const QuadraticStencilGroup& group = groups.data[index];
          if (group.row_count != 4U || group.row_begin > UINT32_MAX - 3U) {
            local = {StatusCode::invalid_plan, kImportDonors};
            break;
          }
          boundary_candidate.links_.push_back(
              {static_cast<std::uint32_t>(index),
               static_cast<std::uint32_t>(index), group.row_begin,
               group.row_begin + 1U, group.row_begin + 2U,
               group.row_begin + 3U});
        }
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kImportAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kImportDonors};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  Hash64 boundary_hash;
  boundary_hash.integer(topology_candidate.fingerprint_);
  boundary_hash.integer(imported_surface);
  boundary_hash.integer(boundary_candidate.reconstruction_.fingerprint());
  boundary_hash.integer(
      static_cast<std::uint64_t>(boundary_candidate.links_.size()));
  boundary_candidate.fingerprint_ = boundary_hash.finish();
  boundary_candidate.lowest_failing_rank_ = -1;

  try {
    QuadraticStencilPlan& target = quadrature_candidate.reconstruction_;
    const QuadraticStencilPlan& source = boundary_candidate.reconstruction_;
    target.groups_ = source.groups_;
    target.rows_ = source.rows_;
    target.donor_global_cells_ = source.donor_global_cells_;
    target.donor_local_indices_ = source.donor_local_indices_;
    target.weights_ = source.weights_;
    target.maximum_halo_reach_ = source.maximum_halo_reach_;
    target.periodic_axes_ = source.periodic_axes_;
    target.audit_ = source.audit_;
    target.refresh_fingerprint();
    if (target.fingerprint_ != source.fingerprint_) {
      local = {StatusCode::invalid_plan, kImportDonors};
    }
    if (local) {
      quadrature_candidate.global_point_count_ = global_link_count;
      for (std::size_t index = 0U;
           index < topology_candidate.links_.size(); ++index) {
        const ImmersedLink& link = topology_candidate.links_[index];
        const QuadraticStencilGroup& group =
            quadrature_candidate.reconstruction_.groups_.at(index);
        quadrature_candidate.local_points_.push_back(
            {link.triangle, 0U, link.wall_point, link.solid_to_fluid_normal,
             link.cartesian_control_face_area, rank, link.fluid_cell,
             static_cast<std::uint32_t>(index), group.row_begin + 2U,
             group.row_begin + 3U});
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kImportAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kImportDonors};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) return agreed;

  quadrature_candidate.physical_fingerprint_ = metric_plan.physical_fingerprint_;
  Hash64 quadrature_layout;
  quadrature_layout.integer(quadrature_candidate.physical_fingerprint_);
  quadrature_layout.integer(topology_candidate.fingerprint_);
  quadrature_layout.integer(
      quadrature_candidate.reconstruction_.fingerprint());
  quadrature_layout.integer(rank);
  quadrature_layout.integer(
      static_cast<std::uint64_t>(quadrature_candidate.local_points_.size()));
  for (const SurfaceQuadraturePoint& point :
       quadrature_candidate.local_points_) {
    quadrature_layout.integer(point.triangle);
    quadrature_layout.integer(point.owner_cell);
  }
  quadrature_candidate.local_layout_fingerprint_ =
      quadrature_layout.finish();
  quadrature_candidate.lowest_failing_rank_ = -1;
  topology_candidate.lowest_failing_rank_ = -1;

  topology = std::move(topology_candidate);
  boundary = std::move(boundary_candidate);
  quadrature = std::move(quadrature_candidate);
  return {};
}

}  // namespace hundun::v04
