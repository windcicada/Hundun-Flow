// SPDX-License-Identifier: Apache-2.0

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
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kPlanInput = 13401U;
constexpr std::uint32_t kPlanCollective = 13402U;
constexpr std::uint32_t kPlanGeometry = 13403U;
constexpr std::uint32_t kPlanDonors = 13404U;
constexpr std::uint32_t kPlanQuadrature = 13405U;

constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

class Hash64 {
 public:
  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    const Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t shift = 0U; shift < sizeof(bits) * 8U; shift += 8U) {
      byte(static_cast<std::uint8_t>((bits >> shift) & Unsigned{0xffU}));
    }
  }
  void real(double value) noexcept {
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }
  PlanFingerprint finish() const noexcept { return value_ == 0U ? 1U : value_; }

 private:
  void byte(std::uint8_t value) noexcept {
    value_ ^= value;
    value_ *= kFnvPrime;
  }
  std::uint64_t value_{kFnvOffset};
};

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Real3 subtract(Real3 left, Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Real3 cross(Real3 left, Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double norm_squared(Real3 value) noexcept { return dot(value, value); }

std::int32_t component(Int3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

double centre(const CartesianGeometryPlan& geometry, Int3 global,
              int axis) noexcept {
  const Span<const double> values = geometry.axis(
      static_cast<CartesianAxis>(axis)).centres();
  return values.data[static_cast<std::size_t>(component(global, axis))];
}

double width(const CartesianGeometryPlan& geometry, Int3 global,
             int axis) noexcept {
  const Span<const double> values = geometry.axis(
      static_cast<CartesianAxis>(axis)).widths();
  return values.data[static_cast<std::size_t>(component(global, axis))];
}

Int3 anchor_cell(const CartesianGeometryPlan& geometry,
                 Real3 point) noexcept {
  const auto select = [](Span<const double> centres,
                         double coordinate) noexcept {
    const double* const begin = centres.data;
    const double* const end = centres.data + centres.size;
    const double* const selected = std::lower_bound(begin, end, coordinate);
    const std::ptrdiff_t raw = selected - begin;
    return static_cast<std::int32_t>(
        std::max<std::ptrdiff_t>(
            0, std::min<std::ptrdiff_t>(raw,
                                        static_cast<std::ptrdiff_t>(
                                            centres.size - 1U))));
  };
  return {select(geometry.x().centres(), point.x),
          select(geometry.y().centres(), point.y),
          select(geometry.z().centres(), point.z)};
}

bool inside(Int3 index, Int3 cells) noexcept {
  return index.x >= 0 && index.y >= 0 && index.z >= 0 &&
         index.x < cells.x && index.y < cells.y && index.z < cells.z;
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(const MeshPatch& left, const MeshPatch& right) noexcept {
  return same(left.begin, right.begin) && same(left.cells, right.cells) &&
         same(left.process_grid, right.process_grid) &&
         same(left.process_coord, right.process_coord);
}

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& out) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    return false;
  }
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
  const std::uint64_t grid_x =
      static_cast<std::uint64_t>(patch.process_grid.x);
  const std::uint64_t grid_y =
      static_cast<std::uint64_t>(patch.process_grid.y);
  const std::uint64_t grid_z =
      static_cast<std::uint64_t>(patch.process_grid.z);
  std::uint64_t process_xy = 0U;
  std::uint64_t process_count = 0U;
  std::uint64_t z_offset = 0U;
  std::uint64_t yz_coordinate = 0U;
  std::uint64_t rank_offset = 0U;
  std::uint64_t expected_rank = 0U;
  if (!checked_multiply(grid_x, grid_y, process_xy) ||
      !checked_multiply(process_xy, grid_z, process_count) ||
      !checked_multiply(
          grid_y, static_cast<std::uint64_t>(patch.process_coord.z),
          z_offset) ||
      !checked_add(static_cast<std::uint64_t>(patch.process_coord.y),
                   z_offset, yz_coordinate) ||
      !checked_multiply(grid_x, yz_coordinate, rank_offset) ||
      !checked_add(static_cast<std::uint64_t>(patch.process_coord.x),
                   rank_offset, expected_rank) ||
      process_count != static_cast<std::uint64_t>(size) ||
      expected_rank != static_cast<std::uint64_t>(rank)) {
    return false;
  }
  Int3 expected_begin{};
  Int3 expected_cells{};
  split_axis(global.x, patch.process_grid.x, patch.process_coord.x,
             expected_begin.x, expected_cells.x);
  split_axis(global.y, patch.process_grid.y, patch.process_coord.y,
             expected_begin.y, expected_cells.y);
  split_axis(global.z, patch.process_grid.z, patch.process_coord.z,
             expected_begin.z, expected_cells.z);
  return same(patch.begin, expected_begin) &&
         same(patch.cells, expected_cells);
}

bool valid_stencil_limits(QuadraticStencilLimits limits) noexcept {
  return limits.minimum_donors >= 14U &&
         limits.minimum_donors <= limits.maximum_donors &&
         limits.maximum_donors <= 32U && limits.maximum_reach > 0U &&
         limits.maximum_reach <= 4U &&
         limits.minimum_normal_bands >= 3U &&
         std::isfinite(limits.condition_limit) &&
         limits.condition_limit >= 1.0 && limits.condition_limit <= 1.0e8;
}

GlobalCellId global_id(Int3 index, Int3 cells) noexcept {
  return static_cast<GlobalCellId>(index.x) +
         static_cast<GlobalCellId>(cells.x) *
             (static_cast<GlobalCellId>(index.y) +
              static_cast<GlobalCellId>(cells.y) *
                  static_cast<GlobalCellId>(index.z));
}

int owner_coordinate(std::int32_t global, std::int32_t partitions,
                     std::int32_t index) noexcept {
  const std::int32_t quotient = global / partitions;
  const std::int32_t remainder = global % partitions;
  const std::int32_t long_end = (quotient + 1) * remainder;
  if (index < long_end) {
    return index / (quotient + 1);
  }
  return remainder + (index - long_end) / quotient;
}

int owner_rank(Int3 index, Int3 global, Int3 process_grid) noexcept {
  const int x = owner_coordinate(global.x, process_grid.x, index.x);
  const int y = owner_coordinate(global.y, process_grid.y, index.y);
  const int z = owner_coordinate(global.z, process_grid.z, index.z);
  return x + process_grid.x * (y + process_grid.y * z);
}

Real3 quadrature_position(const SurfaceTriangle& triangle,
                          std::size_t point_index) noexcept {
  constexpr std::array<std::array<double, 3U>, 3U> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  const std::array<double, 3U>& weights = barycentric[point_index];
  return {weights[0] * triangle.vertices[0].x +
              weights[1] * triangle.vertices[1].x +
              weights[2] * triangle.vertices[2].x,
          weights[0] * triangle.vertices[0].y +
              weights[1] * triangle.vertices[1].y +
              weights[2] * triangle.vertices[2].y,
          weights[0] * triangle.vertices[0].z +
              weights[1] * triangle.vertices[1].z +
              weights[2] * triangle.vertices[2].z};
}

Status consensus(MPI_Comm communicator, int rank, Status local,
                 int& lowest) noexcept {
  const int local_failure = local ? std::numeric_limits<int>::max() : rank;
  int first = std::numeric_limits<int>::max();
  if (MPI_Allreduce(&local_failure, &first, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kPlanCollective};
  }
  lowest = first == std::numeric_limits<int>::max() ? -1 : first;
  std::array<std::uint32_t, 2U> wire{};
  if (rank == lowest) {
    wire = {static_cast<std::uint32_t>(local.code), local.detail};
  }
  if (lowest >= 0 &&
      MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT32_T,
                lowest, communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kPlanCollective};
  }
  return lowest < 0 ? Status{}
                    : Status{static_cast<StatusCode>(wire[0]), wire[1]};
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
    return {StatusCode::mpi_failure, kPlanCollective};
  }
  return {};
}

PlanFingerprint shared_compile_contract(
    const CartesianGeometryPlan& geometry,
    const ImmersedSurfacePlan& surface, const EBTopology& topology,
    const MeshPatch& patch, ImmersedPlanLimits limits, int communicator_size,
    std::uint8_t compiler_kind) noexcept {
  Hash64 hash;
  hash.integer(compiler_kind);
  hash.integer(geometry.fingerprint());
  hash.integer(geometry.topology_revision());
  hash.integer(surface.fingerprint());
  hash.integer(topology.fingerprint());
  hash.integer(topology.geometry_fingerprint());
  hash.integer(topology.surface_fingerprint());
  hash.integer(static_cast<std::uint8_t>(topology.fluid_side()));
  hash.integer(patch.process_grid.x);
  hash.integer(patch.process_grid.y);
  hash.integer(patch.process_grid.z);
  hash.integer(communicator_size);
  hash.integer(limits.stencil.minimum_donors);
  hash.integer(limits.stencil.maximum_donors);
  hash.integer(limits.stencil.maximum_reach);
  hash.integer(limits.stencil.minimum_normal_bands);
  hash.real(limits.stencil.condition_limit);
  hash.integer(limits.maximum_persistent_bytes_per_rank);
  hash.integer(limits.maximum_peak_bytes_per_rank);
  hash.integer(limits.maximum_local_links);
  hash.integer(limits.maximum_local_quadrature_points);
  return hash.finish();
}

Status agree_compile_contract(MPI_Comm communicator, int rank,
                              PlanFingerprint local_contract,
                              int& lowest) noexcept {
  PlanFingerprint root_contract = local_contract;
  Status local{};
  if (MPI_Bcast(&root_contract, 1, MPI_UINT64_T, 0, communicator) !=
      MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kPlanCollective};
  } else if (root_contract != local_contract) {
    local = {StatusCode::invalid_plan, kPlanInput};
  }
  return consensus(communicator, rank, local, lowest);
}

struct CandidateDonor {
  QuadraticDonorCell donor{};
  double distance_squared{};
  double normal_coordinate{};
  double tangent1_coordinate{};
  double tangent2_coordinate{};
};

struct DonorSearchScratch {
  std::vector<CandidateDonor> candidates;
  std::vector<std::uint8_t> selected;
  std::vector<std::size_t> normal_order;
};

struct PreparedQuadraticStencil {
  QuadraticFrame frame{};
  std::size_t donor_begin{};
  std::size_t donor_count{};
  std::size_t functional_begin{};
  std::size_t functional_count{};
};

bool plan_memory_within_limits(std::uint64_t group_count,
                               std::uint64_t row_count,
                               std::uint64_t donor_count,
                               std::uint64_t weight_count,
                               std::uint64_t mapping_count,
                               std::uint64_t mapping_bytes,
                               ImmersedPlanLimits limits) noexcept {
  std::uint64_t persistent = 0U;
  std::uint64_t bytes = 0U;
  const auto add_array = [&](std::uint64_t count,
                             std::uint64_t element_bytes) noexcept {
    return checked_multiply(count, element_bytes, bytes) &&
           checked_add(persistent, bytes, persistent);
  };
  if (!add_array(group_count, sizeof(QuadraticStencilGroup)) ||
      !add_array(row_count, sizeof(QuadraticAffineRow)) ||
      !add_array(donor_count, sizeof(GlobalCellId)) ||
      !add_array(donor_count, sizeof(Int3)) ||
      !add_array(weight_count, sizeof(double)) ||
      !add_array(mapping_count, mapping_bytes) ||
      persistent > limits.maximum_persistent_bytes_per_rank) {
    return false;
  }

  std::uint64_t peak = persistent;
  const auto add_scratch = [&](std::uint64_t count,
                               std::uint64_t element_bytes) noexcept {
    return checked_multiply(count, element_bytes, bytes) &&
           checked_add(peak, bytes, peak);
  };
  const std::uint64_t reach = limits.stencil.maximum_reach;
  const std::uint64_t diameter = 2U * reach + 1U;
  std::uint64_t search_count = 0U;
  if (!checked_multiply(diameter, diameter, search_count) ||
      !checked_multiply(search_count, diameter, search_count) ||
      !add_scratch(group_count, sizeof(PreparedQuadraticStencil)) ||
      !add_scratch(group_count, sizeof(QuadraticStencilRequest)) ||
      !add_scratch(donor_count, sizeof(QuadraticDonorCell)) ||
      !add_scratch(row_count, sizeof(QuadraticFunctionalRequest)) ||
      !add_scratch(search_count, sizeof(CandidateDonor)) ||
      !add_scratch(search_count, sizeof(std::uint8_t)) ||
      !add_scratch(search_count, sizeof(std::size_t)) ||
      !add_scratch(limits.stencil.maximum_donors,
                   sizeof(QuadraticDonorCell))) {
    return false;
  }
  return peak <= limits.maximum_peak_bytes_per_rank;
}

bool worst_case_plan_counts(std::uint64_t group_count,
                            std::uint64_t rows_per_group,
                            QuadraticStencilLimits limits,
                            std::uint64_t& row_count,
                            std::uint64_t& donor_count,
                            std::uint64_t& weight_count) noexcept {
  return checked_multiply(group_count, rows_per_group, row_count) &&
         checked_multiply(group_count, limits.maximum_donors, donor_count) &&
         checked_multiply(donor_count, rows_per_group, weight_count) &&
         group_count <= UINT32_MAX && row_count <= UINT32_MAX &&
         donor_count <= UINT32_MAX && weight_count <= UINT32_MAX;
}

bool donor_less(const CandidateDonor& left,
                const CandidateDonor& right) noexcept {
  if (left.distance_squared != right.distance_squared) {
    return left.distance_squared < right.distance_squared;
  }
  return left.donor.global_cell < right.donor.global_cell;
}

std::uint8_t donor_quadrant(const CandidateDonor& donor) noexcept {
  const std::uint8_t high_t1 = donor.tangent1_coordinate >= 0.0 ? 1U : 0U;
  const std::uint8_t high_t2 = donor.tangent2_coordinate >= 0.0 ? 1U : 0U;
  return static_cast<std::uint8_t>(high_t1 * 2U + high_t2);
}

bool same_normal_band(double left, double right) noexcept {
  return std::abs(left - right) <=
         256.0 * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(left), std::abs(right)});
}

void choose_coverage_complete(
    const std::vector<CandidateDonor>& candidates, std::size_t maximum,
    DonorSearchScratch& scratch, std::vector<QuadraticDonorCell>& out) {
  scratch.selected.assign(candidates.size(), 0U);
  out.clear();
  out.reserve(maximum);
  const auto take = [&](std::size_t index) {
    if (out.size() < maximum && scratch.selected[index] == 0U) {
      scratch.selected[index] = 1U;
      out.push_back(candidates[index].donor);
    }
  };

  // Seed every tangential quadrant with its nearest donor.  This prevents a
  // curved or oblique wall from filling the fixed-width row with only the
  // geometrically closest half-plane before the strict coverage gate runs.
  for (std::uint8_t quadrant = 0U; quadrant < 4U; ++quadrant) {
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
      if (donor_quadrant(candidates[index]) == quadrant) {
        take(index);
        break;
      }
    }
  }

  // Seed three distinct positive-normal bands, preferring the closest donor
  // in each band.  Bands are mathematical coordinate bands, not case names or
  // Cartesian direction heuristics.
  scratch.normal_order.resize(candidates.size());
  for (std::size_t index = 0U; index < candidates.size(); ++index) {
    scratch.normal_order[index] = index;
  }
  std::sort(scratch.normal_order.begin(), scratch.normal_order.end(),
            [&](std::size_t left, std::size_t right) {
              if (candidates[left].normal_coordinate !=
                  candidates[right].normal_coordinate) {
                return candidates[left].normal_coordinate <
                       candidates[right].normal_coordinate;
              }
              return donor_less(candidates[left], candidates[right]);
            });
  double previous = 0.0;
  bool have_previous = false;
  std::uint8_t bands = 0U;
  for (const std::size_t index : scratch.normal_order) {
    if (!have_previous ||
        !same_normal_band(previous, candidates[index].normal_coordinate)) {
      take(index);
      previous = candidates[index].normal_coordinate;
      have_previous = true;
      if (++bands == 3U) {
        break;
      }
    }
  }
  for (std::size_t index = 0U;
       index < candidates.size() && out.size() < maximum; ++index) {
    take(index);
  }
}

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

Status collect_donors(const CartesianGeometryPlan& geometry,
                      const MeshPatch& patch, const EBTopology& topology,
                      Real3 origin, Real3 normal,
                      QuadraticStencilLimits limits,
                      DonorSearchScratch& scratch, QuadraticFrame& frame,
                      std::vector<QuadraticDonorCell>& out) {
  if (!finite(origin) || !finite(normal) ||
      std::abs(norm_squared(normal) - 1.0) > 1.0e-8) {
    return {StatusCode::invalid_plan, kPlanGeometry};
  }
  frame.origin = origin;
  frame.normal = normal;
  frame.tangent1 = tangent1_for(normal);
  frame.tangent2 = cross(normal, frame.tangent1);
  std::vector<CandidateDonor>& candidates = scratch.candidates;
  candidates.clear();
  const int reach = static_cast<int>(limits.maximum_reach);
  const Int3 global = geometry.global_cells();
  frame.anchor_global_cell = anchor_cell(geometry, origin);
  frame.scale = std::cbrt(
      width(geometry, frame.anchor_global_cell, 0) *
      width(geometry, frame.anchor_global_cell, 1) *
      width(geometry, frame.anchor_global_cell, 2));
  for (int dz = -reach; dz <= reach; ++dz) {
    for (int dy = -reach; dy <= reach; ++dy) {
      for (int dx = -reach; dx <= reach; ++dx) {
        const Int3 index{frame.anchor_global_cell.x + dx,
                         frame.anchor_global_cell.y + dy,
                         frame.anchor_global_cell.z + dz};
        if (!inside(index, global) || !topology.is_fluid_global(index)) {
          continue;
        }
        const Real3 position{centre(geometry, index, 0),
                             centre(geometry, index, 1),
                             centre(geometry, index, 2)};
        const Real3 delta = subtract(position, origin);
        const double wall_normal = dot(delta, normal);
        if (!(wall_normal >
              64.0 * std::numeric_limits<double>::epsilon() * frame.scale)) {
          continue;
        }
        QuadraticDonorCell donor;
        donor.global_cell = global_id(index, global);
        donor.global_index = index;
        donor.local_index = {index.x - patch.begin.x,
                             index.y - patch.begin.y,
                             index.z - patch.begin.z};
        donor.centre = position;
        donor.widths = {width(geometry, index, 0), width(geometry, index, 1),
                        width(geometry, index, 2)};
        candidates.push_back(
            {donor, norm_squared(delta), wall_normal / frame.scale,
             dot(delta, frame.tangent1) / frame.scale,
             dot(delta, frame.tangent2) / frame.scale});
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(), donor_less);
  if (candidates.size() < limits.minimum_donors) {
    return {StatusCode::invalid_plan, kPlanDonors};
  }
  const std::size_t count =
      std::min<std::size_t>(candidates.size(), limits.maximum_donors);
  choose_coverage_complete(candidates, count, scratch, out);
  std::sort(out.begin(), out.end(), [](const QuadraticDonorCell& left,
                                       const QuadraticDonorCell& right) {
    return left.global_cell < right.global_cell;
  });
  return {};
}

}  // namespace

void QuadraticStencilPlan::refresh_fingerprint() noexcept {
  Hash64 hash;
  hash.integer(static_cast<std::uint64_t>(groups_.size()));
  hash.integer(static_cast<std::uint64_t>(rows_.size()));
  hash.integer(static_cast<std::uint64_t>(donor_global_cells_.size()));
  hash.integer(static_cast<std::uint64_t>(weights_.size()));
  hash.integer(maximum_halo_reach_);
  for (const QuadraticStencilGroup& group : groups_) {
    hash.integer(group.fingerprint);
  }
  for (const QuadraticAffineRow& row : rows_) {
    hash.integer(row.group);
    hash.integer(row.weight_begin);
    hash.real(row.wall_value_weight);
    hash.real(row.wall_normal_gradient_weight_m);
  }
  for (const GlobalCellId donor : donor_global_cells_) {
    hash.integer(donor);
  }
  for (const Int3 index : donor_local_indices_) {
    hash.integer(index.x);
    hash.integer(index.y);
    hash.integer(index.z);
  }
  for (const double weight : weights_) {
    hash.real(weight);
  }
  fingerprint_ = hash.finish();
}

Status BoundaryStencilCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology, ImmersedPlanLimits limits,
    BoundaryStencilPlan& out) noexcept {
  int rank = -1;
  int size = 0;
  const Status context = mpi_context(communicator, rank, size);
  if (!context) {
    return context;
  }
  Status local{};
  BoundaryStencilPlan candidate;
  try {
    if (geometry.fingerprint() == 0U || topology.fingerprint() == 0U ||
        surface.fingerprint() == 0U ||
        topology.geometry_revision() != geometry.topology_revision() ||
        topology.geometry_fingerprint() != geometry.fingerprint() ||
        topology.surface_fingerprint() != surface.fingerprint() ||
        !same(patch, topology.patch_) ||
        !valid_patch(geometry.global_cells(), patch, rank, size) ||
        !valid_stencil_limits(limits.stencil) ||
        limits.maximum_persistent_bytes_per_rank == 0U ||
        limits.maximum_peak_bytes_per_rank <
            limits.maximum_persistent_bytes_per_rank) {
      local = {StatusCode::invalid_plan, kPlanInput};
    }
    int lowest = -1;
    Status agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    agreed = agree_compile_contract(
        communicator, rank,
        shared_compile_contract(geometry, surface, topology, patch, limits,
                                size, 1U),
        lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    if (local) {
      const std::size_t link_count = topology.links_.size();
      std::uint64_t row_count = 0U;
      std::uint64_t donor_count = 0U;
      std::uint64_t weight_count = 0U;
      if (!worst_case_plan_counts(link_count, 4U, limits.stencil, row_count,
                                  donor_count, weight_count) ||
          !plan_memory_within_limits(
              link_count, row_count, donor_count, weight_count, link_count,
              sizeof(BoundaryStencilLink), limits)) {
        local = {StatusCode::invalid_plan, kPlanDonors};
      }
      std::vector<PreparedQuadraticStencil> prepared;
      std::vector<QuadraticDonorCell> donor_storage;
      std::vector<QuadraticFunctionalRequest> functional_storage;
      std::vector<QuadraticDonorCell> donors;
      DonorSearchScratch donor_scratch;
      std::array<QuadraticFunctionalRequest, 4U> rows{};
      if (local) {
        prepared.reserve(link_count);
        functional_storage.reserve(link_count * rows.size());
        if (link_count > donor_storage.max_size() /
                             limits.stencil.maximum_donors) {
          local = {StatusCode::invalid_plan, kPlanDonors};
        } else {
          donor_storage.reserve(link_count * limits.stencil.maximum_donors);
        }
      }
      for (std::size_t link_index = 0U;
           local && link_index < link_count; ++link_index) {
        const ImmersedLink& link = topology.links_[link_index];
        QuadraticFrame frame;
        local = collect_donors(geometry, patch, topology, link.wall_point,
                               link.solid_to_fluid_normal, limits.stencil,
                               donor_scratch, frame, donors);
        if (!local) {
          break;
        }
        const Real3 ghost{centre(geometry, link.solid_global_index, 0),
                          centre(geometry, link.solid_global_index, 1),
                          centre(geometry, link.solid_global_index, 2)};
        rows[0] = {QuadraticFunctionalKind::value,
                   QuadraticConstraint::origin_value, ghost, {}};
        rows[1] = {QuadraticFunctionalKind::value,
                   QuadraticConstraint::origin_normal_gradient, ghost, {}};
        rows[2] = {QuadraticFunctionalKind::value,
                   QuadraticConstraint::none, link.wall_point, {}};
        rows[3] = {QuadraticFunctionalKind::directional_derivative,
                   QuadraticConstraint::none, link.wall_point,
                   link.solid_to_fluid_normal};
        const std::size_t donor_begin = donor_storage.size();
        const std::size_t functional_begin = functional_storage.size();
        donor_storage.insert(donor_storage.end(), donors.begin(),
                             donors.end());
        functional_storage.insert(functional_storage.end(), rows.begin(),
                                  rows.end());
        prepared.push_back({frame, donor_begin, donors.size(),
                            functional_begin, rows.size()});
      }
      std::vector<QuadraticStencilRequest> requests;
      if (local && !prepared.empty()) {
        // Construct non-owning spans only after every owning vector has
        // finished growing.  This keeps all request pointers stable during
        // the single batch factorization pass.
        requests.reserve(prepared.size());
        for (const PreparedQuadraticStencil& stencil : prepared) {
          requests.push_back(
              {stencil.frame,
               {donor_storage.data() + stencil.donor_begin,
                stencil.donor_count},
               {functional_storage.data() + stencil.functional_begin,
                stencil.functional_count}});
        }
        local = QuadraticStencilCompiler::compile(
            {requests.data(), requests.size()}, limits.stencil,
            candidate.reconstruction_);
      }
      if (local) {
        const Span<const QuadraticStencilGroup> groups =
            candidate.reconstruction_.groups();
        if (groups.size != prepared.size()) {
          local = {StatusCode::invalid_plan, kPlanDonors};
        } else {
          candidate.links_.reserve(link_count);
          for (std::size_t link_index = 0U; link_index < link_count;
               ++link_index) {
            const QuadraticStencilGroup& group = groups.data[link_index];
            if (group.row_count != 4U || group.row_begin > UINT32_MAX - 3U) {
              local = {StatusCode::invalid_plan, kPlanDonors};
              break;
            }
            candidate.links_.push_back(
                {static_cast<std::uint32_t>(link_index),
                 static_cast<std::uint32_t>(link_index), group.row_begin,
                 group.row_begin + 1U, group.row_begin + 2U,
                 group.row_begin + 3U});
          }
        }
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kPlanDonors};
  } catch (...) {
    local = {StatusCode::invalid_plan, kPlanDonors};
  }
  int lowest = -1;
  const Status agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  candidate.lowest_failing_rank_ = -1;
  candidate.reconstruction_.refresh_fingerprint();
  Hash64 hash;
  hash.integer(topology.fingerprint());
  hash.integer(surface.fingerprint());
  hash.integer(candidate.reconstruction_.fingerprint());
  hash.integer(static_cast<std::uint64_t>(candidate.links_.size()));
  candidate.fingerprint_ = hash.finish();
  out = std::move(candidate);
  return {};
}

Status SurfaceQuadratureCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology, ImmersedPlanLimits limits,
    SurfaceQuadraturePlan& out) noexcept {
  int rank = -1;
  int size = 0;
  const Status context = mpi_context(communicator, rank, size);
  if (!context) {
    return context;
  }
  Status local{};
  SurfaceQuadraturePlan candidate;
  try {
    if (geometry.fingerprint() == 0U || topology.fingerprint() == 0U ||
        surface.fingerprint() == 0U ||
        topology.geometry_revision() != geometry.topology_revision() ||
        topology.geometry_fingerprint() != geometry.fingerprint() ||
        topology.surface_fingerprint() != surface.fingerprint() ||
        !same(patch, topology.patch_) ||
        !valid_patch(geometry.global_cells(), patch, rank, size) ||
        !valid_stencil_limits(limits.stencil) ||
        limits.maximum_local_quadrature_points == 0U ||
        limits.maximum_persistent_bytes_per_rank == 0U ||
        limits.maximum_peak_bytes_per_rank <
            limits.maximum_persistent_bytes_per_rank) {
      local = {StatusCode::invalid_plan, kPlanInput};
    }
    int lowest = -1;
    Status agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    agreed = agree_compile_contract(
        communicator, rank,
        shared_compile_contract(geometry, surface, topology, patch, limits,
                                size, 2U),
        lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    const Span<const SurfaceTriangle> surface_triangles = surface.triangles();
    const std::uint64_t triangle_count = surface_triangles.size;
    if (local &&
        !checked_multiply(triangle_count, 3U,
                          candidate.global_point_count_)) {
      local = {StatusCode::invalid_plan, kPlanQuadrature};
    }
    std::uint64_t local_point_count = 0U;
    if (local) {
      // Ownership is a pure Cartesian lookup.  Count it before allocating so
      // every rank sizes the cold plan from its exact local work, not from the
      // replicated STL point count.
      for (std::size_t triangle_index = 0U;
           local && triangle_index < surface_triangles.size;
           ++triangle_index) {
        const SurfaceTriangle& triangle =
            surface_triangles.data[triangle_index];
        for (std::size_t point_index = 0U; point_index < 3U; ++point_index) {
          const Real3 point = quadrature_position(triangle, point_index);
          const Int3 execution_anchor = anchor_cell(geometry, point);
          const int selected_rank = owner_rank(execution_anchor,
                                               geometry.global_cells(),
                                               patch.process_grid);
          if (selected_rank == rank &&
              !checked_add(local_point_count, 1U, local_point_count)) {
            local = {StatusCode::invalid_plan, kPlanQuadrature};
            break;
          }
        }
      }
    }
    if (local) {
      std::uint64_t row_count = 0U;
      std::uint64_t donor_count = 0U;
      std::uint64_t weight_count = 0U;
      if (local_point_count > limits.maximum_local_quadrature_points ||
          !worst_case_plan_counts(local_point_count, 2U, limits.stencil,
                                  row_count, donor_count, weight_count) ||
          !plan_memory_within_limits(
              local_point_count, row_count, donor_count, weight_count,
              local_point_count, sizeof(SurfaceQuadraturePoint), limits)) {
        local = {StatusCode::invalid_plan, kPlanQuadrature};
      }
    }
    // Do not let one rank enter donor construction while another has already
    // rejected its exact local count or memory budget.
    agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    if (local) {
      std::vector<PreparedQuadraticStencil> prepared;
      std::vector<QuadraticDonorCell> donor_storage;
      std::vector<QuadraticFunctionalRequest> functional_storage;
      std::vector<QuadraticDonorCell> donors;
      DonorSearchScratch donor_scratch;
      std::array<QuadraticFunctionalRequest, 2U> rows{};
      const std::size_t exact_local_count =
          static_cast<std::size_t>(local_point_count);
      const std::size_t exact_row_count = exact_local_count * rows.size();
      const std::size_t exact_donor_capacity =
          exact_local_count * limits.stencil.maximum_donors;
      prepared.reserve(exact_local_count);
      candidate.local_points_.reserve(exact_local_count);
      functional_storage.reserve(exact_row_count);
      donor_storage.reserve(exact_donor_capacity);
      for (std::size_t triangle_index = 0U;
           triangle_index < surface_triangles.size; ++triangle_index) {
        const SurfaceTriangle& triangle =
            surface_triangles.data[triangle_index];
        const Real3 solid_to_fluid =
            topology.fluid_side_ == ImmersedFluidSide::outside
                ? triangle.geometric_outward_normal
                : Real3{-triangle.geometric_outward_normal.x,
                        -triangle.geometric_outward_normal.y,
                        -triangle.geometric_outward_normal.z};
        for (std::size_t point_index = 0U; point_index < 3U; ++point_index) {
          const Real3 point = quadrature_position(triangle, point_index);

          const Int3 execution_anchor = anchor_cell(geometry, point);
          const int selected_rank = owner_rank(execution_anchor,
                                               geometry.global_cells(),
                                               patch.process_grid);
          if (selected_rank != rank) {
            continue;
          }
          if (candidate.local_points_.size() >= exact_local_count) {
            local = {StatusCode::invalid_plan, kPlanQuadrature};
            break;
          }

          QuadraticFrame frame;
          local = collect_donors(geometry, patch, topology, point,
                                 solid_to_fluid, limits.stencil, donor_scratch,
                                 frame, donors);
          if (!local) {
            break;
          }
          GlobalCellId owner_cell = kInvalidSurfaceTriangle;
          double owner_distance = std::numeric_limits<double>::infinity();
          for (const QuadraticDonorCell& donor : donors) {
            const double distance = norm_squared(subtract(donor.centre, point));
            if (distance < owner_distance ||
                (distance == owner_distance &&
                 donor.global_cell < owner_cell)) {
              owner_distance = distance;
              owner_cell = donor.global_cell;
            }
          }
          // Pressure and other unknown wall values are reconstructed from the
          // fluid-side donors.  They must not depend on a caller-supplied wall
          // value (which is only meaningful for prescribed wall fields).
          rows[0] = {QuadraticFunctionalKind::value,
                     QuadraticConstraint::none, point, {}};
          rows[1] = {QuadraticFunctionalKind::directional_derivative,
                     QuadraticConstraint::origin_value, point,
                     solid_to_fluid};
          const std::size_t donor_begin = donor_storage.size();
          const std::size_t functional_begin = functional_storage.size();
          donor_storage.insert(donor_storage.end(), donors.begin(),
                               donors.end());
          functional_storage.insert(functional_storage.end(), rows.begin(),
                                    rows.end());
          prepared.push_back({frame, donor_begin, donors.size(),
                              functional_begin, rows.size()});
          candidate.local_points_.push_back(
              {triangle.id, static_cast<std::uint8_t>(point_index), point,
               solid_to_fluid, triangle.area / 3.0, rank, owner_cell,
               kInvalidIbmIndex, kInvalidIbmIndex, kInvalidIbmIndex});
        }
        if (!local) {
          break;
        }
      }
      if (local && candidate.local_points_.size() != exact_local_count) {
        local = {StatusCode::invalid_plan, kPlanQuadrature};
      }
      std::vector<QuadraticStencilRequest> requests;
      if (local && !prepared.empty()) {
        requests.reserve(prepared.size());
        for (const PreparedQuadraticStencil& stencil : prepared) {
          requests.push_back(
              {stencil.frame,
               {donor_storage.data() + stencil.donor_begin,
                stencil.donor_count},
               {functional_storage.data() + stencil.functional_begin,
                stencil.functional_count}});
        }
        local = QuadraticStencilCompiler::compile(
            {requests.data(), requests.size()}, limits.stencil,
            candidate.reconstruction_);
      }
      if (local && !prepared.empty()) {
        const Span<const QuadraticStencilGroup> groups =
            candidate.reconstruction_.groups();
        if (groups.size != prepared.size() ||
            groups.size != candidate.local_points_.size()) {
          local = {StatusCode::invalid_plan, kPlanQuadrature};
        } else {
          for (std::size_t point_index = 0U; point_index < groups.size;
               ++point_index) {
            const QuadraticStencilGroup& group = groups.data[point_index];
            if (point_index > UINT32_MAX || group.row_count != 2U ||
                group.row_begin == UINT32_MAX) {
              local = {StatusCode::invalid_plan, kPlanQuadrature};
              break;
            }
            SurfaceQuadraturePoint& point =
                candidate.local_points_[point_index];
            point.reconstruction_group =
                static_cast<std::uint32_t>(point_index);
            point.wall_value_row = group.row_begin;
            point.wall_normal_gradient_row = group.row_begin + 1U;
          }
        }
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kPlanQuadrature};
  } catch (...) {
    local = {StatusCode::invalid_plan, kPlanQuadrature};
  }
  int lowest = -1;
  const Status agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  candidate.reconstruction_.refresh_fingerprint();
  Hash64 physical;
  physical.integer(surface.fingerprint());
  physical.integer(topology.fluid_side_ == ImmersedFluidSide::outside ? 0U
                                                                       : 1U);
  physical.integer(candidate.global_point_count_);
  candidate.physical_fingerprint_ = physical.finish();
  Hash64 layout;
  layout.integer(candidate.physical_fingerprint_);
  layout.integer(geometry.fingerprint());
  layout.integer(topology.fingerprint());
  layout.integer(candidate.reconstruction_.fingerprint());
  layout.integer(rank);
  layout.integer(static_cast<std::uint64_t>(candidate.local_points_.size()));
  for (const SurfaceQuadraturePoint& point : candidate.local_points_) {
    layout.integer(point.triangle);
    layout.integer(point.point_index);
    layout.integer(point.owner_cell);
  }
  candidate.local_layout_fingerprint_ = layout.finish();
  candidate.lowest_failing_rank_ = -1;
  out = std::move(candidate);
  return {};
}

}  // namespace hundun::v04
