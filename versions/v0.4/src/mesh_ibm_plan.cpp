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
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kPlanInput = 13401U;
constexpr std::uint32_t kPlanCollective = 13402U;
constexpr std::uint32_t kPlanGeometry = 13403U;
constexpr std::uint32_t kPlanDonors = 13404U;
constexpr std::uint32_t kPlanQuadrature = 13405U;
constexpr std::uint32_t kPlanBoundaryIntersection = 13406U;
constexpr std::uint32_t kPlanInterfaceMetric = 13407U;
constexpr std::uint32_t kInterfaceMetricSchemaRevision = 1U;

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

int face_for(int axis, bool upper) noexcept {
  return 2 * axis + (upper ? 1 : 0);
}

bool periodic_pair(const ImmersedDomainBoundaryPolicy& policy,
                   int axis) noexcept {
  return policy.allow_periodic_images[face_for(axis, false)] &&
         policy.allow_periodic_images[face_for(axis, true)];
}

bool periodic_policy_valid(const ImmersedDomainBoundaryPolicy& policy,
                           Int3 global, int reach) noexcept {
  const std::int32_t extents[3]{global.x, global.y, global.z};
  for (int axis = 0; axis < 3; ++axis) {
    const bool lower =
        policy.allow_periodic_images[face_for(axis, false)];
    const bool upper = policy.allow_periodic_images[face_for(axis, true)];
    if (lower != upper) return false;
    if ((lower || upper) && extents[axis] <= 2 * reach) return false;
  }
  return true;
}

bool wrap_periodic_index(std::int32_t raw, std::int32_t extent, int axis,
                         const ImmersedDomainBoundaryPolicy& policy,
                         int reach, std::int32_t& canonical,
                         double& image_shift, double domain_lower,
                         double domain_upper) noexcept {
  canonical = raw;
  image_shift = 0.0;
  if (raw >= 0 && raw < extent) return true;
  if (!periodic_pair(policy, axis) || extent <= 2 * reach) return false;
  if (raw < 0) {
    if (raw < -reach) return false;
    canonical = raw + extent;
    // The lower logical image is represented by the upper canonical cell,
    // translated back by one domain length.
    image_shift = domain_lower - domain_upper;
    return canonical >= 0 && canonical < extent;
  }
  if (raw >= extent) {
    if (raw >= extent + reach) return false;
    canonical = raw - extent;
    // The upper logical image is represented by the lower canonical cell,
    // translated forward by one domain length.
    image_shift = domain_upper - domain_lower;
    return canonical >= 0 && canonical < extent;
  }
  return false;
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
         limits.condition_limit >= 1.0 && limits.condition_limit <= 1.0e8 &&
         limits.policy <= IbmReconstructionPolicy::adaptive_order &&
         limits.minimum_linear_donors >= 4U &&
         limits.minimum_linear_donors <= limits.maximum_donors &&
         limits.minimum_linear_normal_bands >= 2U &&
         limits.minimum_linear_normal_bands <=
             limits.minimum_normal_bands &&
         limits.standard_reach > 0U &&
         limits.standard_reach <= limits.maximum_reach;
}

GlobalCellId global_id(Int3 index, Int3 cells) noexcept {
  return static_cast<GlobalCellId>(index.x) +
         static_cast<GlobalCellId>(cells.x) *
             (static_cast<GlobalCellId>(index.y) +
              static_cast<GlobalCellId>(cells.y) *
                  static_cast<GlobalCellId>(index.z));
}

Int3 decode_global_id(GlobalCellId id, Int3 cells, bool& valid) noexcept {
  valid = false;
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) return {};
  const std::uint64_t gx = static_cast<std::uint64_t>(cells.x);
  const std::uint64_t gy = static_cast<std::uint64_t>(cells.y);
  const std::uint64_t gz = static_cast<std::uint64_t>(cells.z);
  if (gx > std::numeric_limits<std::uint64_t>::max() / gy) return {};
  const std::uint64_t plane = gx * gy;
  if (plane > std::numeric_limits<std::uint64_t>::max() / gz) return {};
  const std::uint64_t total = plane * gz;
  if (id >= total) return {};
  const std::uint64_t z = id / plane;
  const std::uint64_t remainder = id - z * plane;
  const std::uint64_t y = remainder / gx;
  const std::uint64_t x = remainder - y * gx;
  valid = true;
  return {static_cast<std::int32_t>(x), static_cast<std::int32_t>(y),
          static_cast<std::int32_t>(z)};
}

bool offset_index(std::int32_t base, int offset,
                  std::int32_t& result) noexcept {
  const std::int64_t value = static_cast<std::int64_t>(base) + offset;
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max()) {
    return false;
  }
  result = static_cast<std::int32_t>(value);
  return true;
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

double real_component(Real3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

void set_real_component(Real3& value, int axis, double component_value) noexcept {
  if (axis == 0) {
    value.x = component_value;
  } else if (axis == 1) {
    value.y = component_value;
  } else {
    value.z = component_value;
  }
}

bool same_point(Real3 left, Real3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

struct ClippedPolygon {
  std::array<Real3, 16U> vertices{};
  std::size_t size{};
};

bool append_clipped_vertex(ClippedPolygon& polygon, Real3 point) noexcept {
  if (polygon.size != 0U &&
      same_point(polygon.vertices[polygon.size - 1U], point)) {
    return true;
  }
  if (polygon.size >= polygon.vertices.size()) return false;
  polygon.vertices[polygon.size++] = point;
  return true;
}

bool clip_halfspace(const ClippedPolygon& input, int axis, double bound,
                    bool keep_greater, ClippedPolygon& output) noexcept {
  output.size = 0U;
  if (input.size == 0U) return true;
  Real3 previous = input.vertices[input.size - 1U];
  double previous_value = real_component(previous, axis);
  bool previous_inside =
      keep_greater ? previous_value >= bound : previous_value <= bound;
  for (std::size_t index = 0U; index < input.size; ++index) {
    const Real3 current = input.vertices[index];
    const double current_value = real_component(current, axis);
    const bool current_inside =
        keep_greater ? current_value >= bound : current_value <= bound;
    if (current_inside != previous_inside) {
      const double denominator = current_value - previous_value;
      if (denominator == 0.0 || !std::isfinite(denominator)) return false;
      const double fraction = (bound - previous_value) / denominator;
      if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
        return false;
      Real3 intersection{
          previous.x + fraction * (current.x - previous.x),
          previous.y + fraction * (current.y - previous.y),
          previous.z + fraction * (current.z - previous.z)};
      set_real_component(intersection, axis, bound);
      if (!append_clipped_vertex(output, intersection)) return false;
    }
    if (current_inside && !append_clipped_vertex(output, current)) return false;
    previous = current;
    previous_value = current_value;
    previous_inside = current_inside;
  }
  if (output.size > 1U &&
      same_point(output.vertices.front(), output.vertices[output.size - 1U])) {
    --output.size;
  }
  return true;
}

bool lies_on_domain_face(const std::array<Real3, 3U>& triangle,
                         Real3 lower, Real3 upper) noexcept {
  for (int axis = 0; axis < 3; ++axis) {
    const double lo = real_component(lower, axis);
    const double hi = real_component(upper, axis);
    bool all_lower = true;
    bool all_upper = true;
    for (const Real3 point : triangle) {
      all_lower &= real_component(point, axis) == lo;
      all_upper &= real_component(point, axis) == hi;
    }
    if (all_lower || all_upper) return true;
  }
  return false;
}

template <class Visitor>
Status visit_domain_clipped_quadrature(
    const CartesianGeometryPlan& geometry,
    Span<const SurfaceTriangle> surface_triangles,
    Visitor&& visitor) noexcept {
  const Real3 lower = geometry.lower();
  const Real3 upper = geometry.upper();
  for (std::size_t triangle_index = 0U;
       triangle_index < surface_triangles.size; ++triangle_index) {
    const SurfaceTriangle& source = surface_triangles.data[triangle_index];
    ClippedPolygon first;
    first.size = 3U;
    first.vertices[0U] = source.vertices[0U];
    first.vertices[1U] = source.vertices[1U];
    first.vertices[2U] = source.vertices[2U];
    ClippedPolygon second;
    ClippedPolygon* input = &first;
    ClippedPolygon* output = &second;
    for (int axis = 0; axis < 3; ++axis) {
      if (!clip_halfspace(*input, axis, real_component(lower, axis), true,
                          *output))
        return {StatusCode::invalid_plan, kPlanQuadrature};
      std::swap(input, output);
      if (!clip_halfspace(*input, axis, real_component(upper, axis), false,
                          *output))
        return {StatusCode::invalid_plan, kPlanQuadrature};
      std::swap(input, output);
    }
    if (input->size < 3U) continue;
    std::size_t point_index = 0U;
    for (std::size_t vertex = 1U; vertex + 1U < input->size; ++vertex) {
      const std::array<Real3, 3U> clipped{
          input->vertices[0U], input->vertices[vertex],
          input->vertices[vertex + 1U]};
      if (lies_on_domain_face(clipped, lower, upper)) continue;
      const Real3 edge1 = subtract(clipped[1U], clipped[0U]);
      const Real3 edge2 = subtract(clipped[2U], clipped[0U]);
      const Real3 area_vector = cross(edge1, edge2);
      const double twice_area = std::sqrt(norm_squared(area_vector));
      const double scale = std::max(
          {norm_squared(edge1), norm_squared(edge2), 1.0});
      if (!(twice_area > 256.0 * std::numeric_limits<double>::epsilon() *
                              scale))
        continue;
      if (!std::isfinite(twice_area) ||
          !(dot(area_vector, source.geometric_outward_normal) > 0.0))
        return {StatusCode::invalid_plan, kPlanQuadrature};
      for (std::size_t local = 0U; local < 3U; ++local) {
        if (point_index > UINT8_MAX)
          return {StatusCode::invalid_plan, kPlanQuadrature};
        constexpr std::array<std::array<double, 3U>, 3U> barycentric{{
            {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
            {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
            {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
        }};
        const auto weights = barycentric[local];
        const Real3 point{
            weights[0U] * clipped[0U].x + weights[1U] * clipped[1U].x +
                weights[2U] * clipped[2U].x,
            weights[0U] * clipped[0U].y + weights[1U] * clipped[1U].y +
                weights[2U] * clipped[2U].y,
            weights[0U] * clipped[0U].z + weights[1U] * clipped[1U].z +
                weights[2U] * clipped[2U].z};
        const Status visited = visitor(
            source, static_cast<std::uint8_t>(point_index++), point,
            0.5 * twice_area / 3.0);
        if (!visited) return visited;
      }
    }
  }
  return {};
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
    const MeshPatch& patch, ImmersedDomainBoundaryPolicy boundary_policy,
    ImmersedPlanLimits limits, int communicator_size,
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
  for (const bool allowed : boundary_policy.allow_one_sided_quadratic) {
    hash.integer(allowed ? 1U : 0U);
  }
  for (const bool allowed : boundary_policy.allow_periodic_images) {
    hash.integer(allowed ? 1U : 0U);
  }
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

struct SparseTopologyFluidValue {
  GlobalCellId id{};
  std::uint8_t fluid{};
};

struct SparseTopologyFluidOracle {
  std::vector<SparseTopologyFluidValue> values;

  bool lookup(GlobalCellId id, std::uint8_t& fluid) const noexcept {
    const auto found = std::lower_bound(
        values.begin(), values.end(), id,
        [](const SparseTopologyFluidValue& value, GlobalCellId key) {
          return value.id < key;
        });
    if (found == values.end() || found->id != id) return false;
    fluid = found->fluid;
    return true;
  }
};

bool topology_covers(Int3 global_index, const MeshPatch& patch,
                     const EBTopology& topology) noexcept {
  const std::int64_t halo = topology.region_halo_width();
  const std::int64_t x = static_cast<std::int64_t>(global_index.x) -
                         patch.begin.x;
  const std::int64_t y = static_cast<std::int64_t>(global_index.y) -
                         patch.begin.y;
  const std::int64_t z = static_cast<std::int64_t>(global_index.z) -
                         patch.begin.z;
  return x >= -halo && y >= -halo && z >= -halo &&
         x < static_cast<std::int64_t>(patch.cells.x) + halo &&
         y < static_cast<std::int64_t>(patch.cells.y) + halo &&
         z < static_cast<std::int64_t>(patch.cells.z) + halo;
}

Status append_query_ids_for_point(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    const EBTopology& topology, Real3 point,
    ImmersedDomainBoundaryPolicy boundary_policy,
    QuadraticStencilLimits limits, std::uint64_t maximum_peak_bytes,
    std::vector<GlobalCellId>& ids) {
  const int reach = static_cast<int>(limits.maximum_reach);
  const Int3 global = geometry.global_cells();
  if (maximum_peak_bytes < sizeof(GlobalCellId)) {
    return {StatusCode::invalid_plan, kPlanDonors};
  }
  const std::uint64_t maximum_ids_u =
      maximum_peak_bytes / sizeof(GlobalCellId);
  const std::size_t maximum_ids =
      maximum_ids_u > ids.max_size()
          ? ids.max_size()
          : static_cast<std::size_t>(maximum_ids_u);
  if (ids.capacity() > maximum_ids) {
    return {StatusCode::invalid_plan, kPlanDonors};
  }
  const Int3 anchor = anchor_cell(geometry, point);
  const double domain_lower[3]{geometry.lower().x, geometry.lower().y,
                               geometry.lower().z};
  const double domain_upper[3]{geometry.upper().x, geometry.upper().y,
                               geometry.upper().z};
  const std::int32_t extents[3]{global.x, global.y, global.z};
  for (int dz = -reach; dz <= reach; ++dz) {
    for (int dy = -reach; dy <= reach; ++dy) {
      for (int dx = -reach; dx <= reach; ++dx) {
        Int3 raw{};
        if (!offset_index(anchor.x, dx, raw.x) ||
            !offset_index(anchor.y, dy, raw.y) ||
            !offset_index(anchor.z, dz, raw.z)) {
          return {StatusCode::invalid_plan, kPlanDonors};
        }
        const std::int32_t raw_values[3]{raw.x, raw.y, raw.z};
        Int3 canonical{};
        bool mapped = true;
        for (int axis = 0; axis < 3; ++axis) {
          std::int32_t value = 0;
          double ignored_shift = 0.0;
          if (!wrap_periodic_index(
                  raw_values[axis], extents[axis], axis, boundary_policy,
                  reach, value, ignored_shift, domain_lower[axis],
                  domain_upper[axis])) {
            mapped = false;
            break;
          }
          if (axis == 0)
            canonical.x = value;
          else if (axis == 1)
            canonical.y = value;
          else
            canonical.z = value;
        }
        if (!mapped || !inside(canonical, global) ||
            topology_covers(canonical, patch, topology)) {
          continue;
        }
        if (ids.size() >= maximum_ids) {
          return {StatusCode::invalid_plan, kPlanDonors};
        }
        if (ids.size() == ids.capacity()) {
          const std::size_t current = ids.capacity();
          const std::size_t next =
              current == 0U
                  ? 1U
                  : current > maximum_ids - current ? maximum_ids
                                                     : current * 2U;
          if (next <= current || next > maximum_ids) {
            return {StatusCode::invalid_plan, kPlanDonors};
          }
          ids.reserve(next);
        }
        ids.push_back(global_id(canonical, global));
        if (ids.capacity() > maximum_ids) {
          return {StatusCode::invalid_plan, kPlanDonors};
        }
      }
    }
  }
  return {};
}

struct SparseTopologyQueryRequest {
  int owner{};
  GlobalCellId id{};
};

bool query_request_less(const SparseTopologyQueryRequest& left,
                        const SparseTopologyQueryRequest& right) noexcept {
  return left.owner < right.owner ||
         (left.owner == right.owner && left.id < right.id);
}

bool topology_region_shape_valid(const EBTopology& topology,
                                 Int3 cells) noexcept {
  std::uint64_t expected = 0U;
  if (!checked_multiply(static_cast<std::uint64_t>(cells.x),
                        static_cast<std::uint64_t>(cells.y), expected) ||
      !checked_multiply(expected, static_cast<std::uint64_t>(cells.z),
                        expected)) {
    return false;
  }
  return expected == topology.region().size;
}

std::size_t local_region_offset(Int3 local, Int3 cells) noexcept {
  return static_cast<std::size_t>(local.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(local.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(local.z));
}

Status query_topology_fluid(
    MPI_Comm communicator, int rank, int size, Int3 global,
    const MeshPatch& patch, const EBTopology& topology,
    const std::vector<GlobalCellId>& query_ids,
    std::uint64_t maximum_peak_bytes,
    SparseTopologyFluidOracle& out) noexcept {
  const int local_has_queries = query_ids.empty() ? 0 : 1;
  int global_has_queries = 0;
  if (MPI_Allreduce(&local_has_queries, &global_has_queries, 1, MPI_INT,
                    MPI_MAX, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPlanCollective};
  }
  if (global_has_queries == 0) {
    out.values.clear();
    return {};
  }
  std::vector<SparseTopologyQueryRequest> requests;
  std::vector<GlobalCellId> request_values;
  std::vector<GlobalCellId> supplied_values;
  std::vector<std::uint8_t> reply_values;
  std::vector<std::uint8_t> received_values;
  std::vector<int> request_counts;
  std::vector<int> request_displacements;
  std::vector<int> supply_counts;
  std::vector<int> supply_displacements;
  Status local{};
  std::uint64_t initial_bytes = 0U;
  try {
    if (query_ids.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        maximum_peak_bytes == 0U ||
        !topology_region_shape_valid(topology, patch.cells)) {
      local = {StatusCode::invalid_plan, kPlanDonors};
    }
    for (std::size_t index = 1U; local && index < query_ids.size(); ++index) {
      if (query_ids[index - 1U] >= query_ids[index]) {
        local = {StatusCode::invalid_plan, kPlanDonors};
      }
    }
    std::uint64_t bytes = 0U;
    std::uint64_t count_bytes = 0U;
    std::uint64_t capacity_bytes = 0U;
    bool memory_ok =
        checked_multiply(static_cast<std::uint64_t>(size), sizeof(int),
                         count_bytes) &&
        checked_multiply(count_bytes, 4U, count_bytes) &&
        checked_multiply(static_cast<std::uint64_t>(query_ids.capacity()),
                         sizeof(GlobalCellId), capacity_bytes) &&
        checked_add(bytes, capacity_bytes, bytes) &&
        checked_multiply(static_cast<std::uint64_t>(query_ids.size()),
                         sizeof(SparseTopologyQueryRequest), capacity_bytes) &&
        checked_add(bytes, capacity_bytes, bytes) &&
        checked_multiply(static_cast<std::uint64_t>(query_ids.size()),
                         sizeof(GlobalCellId), capacity_bytes) &&
        checked_add(bytes, capacity_bytes, bytes) &&
        checked_multiply(static_cast<std::uint64_t>(query_ids.size()),
                         sizeof(SparseTopologyFluidValue), capacity_bytes) &&
        checked_add(bytes, capacity_bytes, bytes) &&
        checked_multiply(static_cast<std::uint64_t>(query_ids.size()),
                         sizeof(std::uint8_t), capacity_bytes) &&
        checked_add(bytes, capacity_bytes, bytes) &&
        checked_add(bytes, count_bytes, bytes);
    if (local && (!memory_ok || bytes > maximum_peak_bytes)) {
      local = {StatusCode::invalid_plan, kPlanDonors};
    } else if (local) {
      initial_bytes = bytes;
    }
    if (local) {
      requests.reserve(query_ids.size());
      request_values.resize(query_ids.size());
      received_values.resize(query_ids.size());
      out.values.resize(query_ids.size());
      request_counts.assign(static_cast<std::size_t>(size), 0);
      request_displacements.assign(static_cast<std::size_t>(size), 0);
      supply_counts.assign(static_cast<std::size_t>(size), 0);
      supply_displacements.assign(static_cast<std::size_t>(size), 0);
      for (const GlobalCellId id : query_ids) {
        bool decoded = false;
        const Int3 index = decode_global_id(id, global, decoded);
        if (!decoded || global_id(index, global) != id) {
          local = {StatusCode::invalid_plan, kPlanDonors};
          break;
        }
        const int owner = owner_rank(index, global, patch.process_grid);
        if (owner < 0 || owner >= size) {
          local = {StatusCode::invalid_plan, kPlanDonors};
          break;
        }
        requests.push_back({owner, id});
      }
      std::sort(requests.begin(), requests.end(), query_request_less);
      for (std::size_t index = 0U; local && index < requests.size(); ++index) {
        if (index != 0U && requests[index - 1U].id == requests[index].id) {
          local = {StatusCode::invalid_plan, kPlanDonors};
          break;
        }
        int& count = request_counts[static_cast<std::size_t>(requests[index].owner)];
        if (count == std::numeric_limits<int>::max()) {
          local = {StatusCode::invalid_plan, kPlanDonors};
          break;
        }
        ++count;
        request_values[index] = requests[index].id;
      }
      int request_total = 0;
      for (int peer = 0; local && peer < size; ++peer) {
        request_displacements[static_cast<std::size_t>(peer)] = request_total;
        const int count = request_counts[static_cast<std::size_t>(peer)];
        if (count < 0 || count > std::numeric_limits<int>::max() - request_total) {
          local = {StatusCode::invalid_plan, kPlanDonors};
          break;
        }
        request_total += count;
      }
      if (local && request_total != static_cast<int>(requests.size())) {
        local = {StatusCode::invalid_plan, kPlanDonors};
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kPlanDonors};
  } catch (...) {
    local = {StatusCode::invalid_plan, kPlanDonors};
  }
  int lowest = -1;
  Status agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  if (MPI_Alltoall(request_counts.data(), 1, MPI_INT, supply_counts.data(), 1,
                   MPI_INT, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPlanCollective};
  }
  local = {};
  int supply_total = 0;
  for (int peer = 0; peer < size; ++peer) {
    supply_displacements[static_cast<std::size_t>(peer)] = supply_total;
    const int count = supply_counts[static_cast<std::size_t>(peer)];
    if (count < 0 || count > std::numeric_limits<int>::max() - supply_total) {
      local = {StatusCode::invalid_plan, kPlanDonors};
      break;
    }
    supply_total += count;
  }
  try {
    std::uint64_t bytes = 0U;
    std::uint64_t inbound = 0U;
    if (!local ||
        !checked_multiply(static_cast<std::uint64_t>(supply_total),
                          sizeof(GlobalCellId), inbound) ||
        !checked_multiply(static_cast<std::uint64_t>(supply_total),
                          sizeof(std::uint8_t), bytes) ||
        !checked_add(bytes, inbound, bytes) ||
        !checked_add(initial_bytes, bytes, inbound) ||
        inbound > maximum_peak_bytes) {
      local = {StatusCode::invalid_plan, kPlanDonors};
    } else {
      supplied_values.resize(static_cast<std::size_t>(supply_total));
      reply_values.resize(static_cast<std::size_t>(supply_total));
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kPlanDonors};
  } catch (...) {
    local = {StatusCode::invalid_plan, kPlanDonors};
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  GlobalCellId dummy_id = 0U;
  if (MPI_Alltoallv(
          request_values.empty() ? &dummy_id : request_values.data(),
          request_counts.data(), request_displacements.data(), MPI_UINT64_T,
          supplied_values.empty() ? &dummy_id : supplied_values.data(),
          supply_counts.data(), supply_displacements.data(), MPI_UINT64_T,
          communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPlanCollective};
  }

  const Span<const std::uint8_t> region = topology.region();
  local = {};
  for (int peer = 0; peer < size; ++peer) {
    const int begin = supply_displacements[static_cast<std::size_t>(peer)];
    const int count = supply_counts[static_cast<std::size_t>(peer)];
    for (int offset = 0; offset < count; ++offset) {
      const std::size_t position =
          static_cast<std::size_t>(begin + offset);
      const GlobalCellId id = supplied_values[position];
      bool decoded = false;
      const Int3 index = decode_global_id(id, global, decoded);
      const Int3 local_index{index.x - patch.begin.x,
                             index.y - patch.begin.y,
                             index.z - patch.begin.z};
      bool valid = decoded && global_id(index, global) == id &&
                   inside(local_index, patch.cells) &&
                   owner_rank(index, global, patch.process_grid) == rank;
      std::uint8_t fluid = 0U;
      if (valid) {
        const std::size_t region_index =
            local_region_offset(local_index, patch.cells);
        if (region_index >= region.size ||
            (region.data[region_index] !=
                 static_cast<std::uint8_t>(RegionFlag::solid) &&
             region.data[region_index] !=
                 static_cast<std::uint8_t>(RegionFlag::fluid))) {
          valid = false;
        } else {
          fluid = region.data[region_index] ==
                          static_cast<std::uint8_t>(RegionFlag::fluid)
                      ? 1U
                      : 0U;
        }
      }
      if (!valid) local = {StatusCode::invalid_plan, kPlanDonors};
      reply_values[position] = fluid;
    }
    for (int offset = 1; offset < count; ++offset) {
      const std::size_t position =
          static_cast<std::size_t>(begin + offset);
      const std::size_t prior =
          static_cast<std::size_t>(begin + offset - 1);
      if (supplied_values[position] == supplied_values[prior]) {
        local = {StatusCode::invalid_plan, kPlanDonors};
      }
    }
  }

  std::uint8_t dummy_byte = 0U;
  if (MPI_Alltoallv(
          reply_values.empty() ? &dummy_byte : reply_values.data(),
          supply_counts.data(), supply_displacements.data(), MPI_UINT8_T,
          received_values.empty() ? &dummy_byte : received_values.data(),
          request_counts.data(), request_displacements.data(), MPI_UINT8_T,
          communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kPlanCollective};
  }
  for (std::size_t index = 0U; index < requests.size(); ++index) {
    if (received_values[index] > 1U) {
      local = {StatusCode::invalid_plan, kPlanDonors};
      break;
    }
    out.values[index] = {requests[index].id, received_values[index]};
  }
  std::sort(out.values.begin(), out.values.end(),
            [](const SparseTopologyFluidValue& left,
               const SparseTopologyFluidValue& right) {
              return left.id < right.id;
            });
  for (std::size_t index = 1U; local && index < out.values.size(); ++index) {
    if (out.values[index - 1U].id == out.values[index].id) {
      local = {StatusCode::invalid_plan, kPlanDonors};
    }
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  return {};
}

Status prepare_topology_fluid_oracle_from_links(
    MPI_Comm communicator, int rank, int size,
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    const EBTopology& topology,
    ImmersedDomainBoundaryPolicy boundary_policy,
    QuadraticStencilLimits stencil_limits, std::uint64_t maximum_peak_bytes,
    SparseTopologyFluidOracle& out) noexcept {
  std::vector<GlobalCellId> query_ids;
  Status local{};
  try {
    const Span<const ImmersedLink> links = topology.links();
    for (std::size_t index = 0U; index < links.size && local; ++index) {
      local = append_query_ids_for_point(
          geometry, patch, topology, links.data[index].wall_point,
          boundary_policy, stencil_limits, maximum_peak_bytes, query_ids);
    }
    std::sort(query_ids.begin(), query_ids.end());
    query_ids.erase(std::unique(query_ids.begin(), query_ids.end()),
                    query_ids.end());
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kPlanDonors};
  } catch (...) {
    local = {StatusCode::invalid_plan, kPlanDonors};
  }
  int lowest = -1;
  Status agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  return query_topology_fluid(communicator, rank, size, geometry.global_cells(),
                              patch, topology, query_ids,
                              maximum_peak_bytes, out);
}

bool plan_memory_within_limits(std::uint64_t group_count,
                               std::uint64_t row_count,
                               std::uint64_t donor_count,
                               std::uint64_t weight_count,
                               std::uint64_t mapping_count,
                               std::uint64_t mapping_bytes,
                               ImmersedPlanLimits limits,
                               std::uint64_t* peak_bytes_out = nullptr) noexcept {
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
  if (peak_bytes_out != nullptr) *peak_bytes_out = peak;
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
    bool domain_truncated, DonorSearchScratch& scratch,
    std::vector<QuadraticDonorCell>& out) {
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
  if (domain_truncated) {
    // A one-sided physical-domain intersection needs tangential breadth, not
    // another cluster of nearest cells.  Deterministic farthest-point fill
    // improves the full quadratic row geometry while the downstream pivoted
    // QR remains the rank/condition authority.
    while (out.size() < maximum) {
      std::size_t best = candidates.size();
      double best_separation = -1.0;
      for (std::size_t index = 0U; index < candidates.size(); ++index) {
        if (scratch.selected[index] != 0U) continue;
        double minimum_separation = std::numeric_limits<double>::infinity();
        for (std::size_t selected = 0U; selected < candidates.size();
             ++selected) {
          if (scratch.selected[selected] == 0U) continue;
          const double dn = candidates[index].normal_coordinate -
                            candidates[selected].normal_coordinate;
          const double dt1 = candidates[index].tangent1_coordinate -
                             candidates[selected].tangent1_coordinate;
          const double dt2 = candidates[index].tangent2_coordinate -
                             candidates[selected].tangent2_coordinate;
          minimum_separation = std::min(
              minimum_separation, dn * dn + dt1 * dt1 + dt2 * dt2);
        }
        if (best == candidates.size() ||
            minimum_separation > best_separation ||
            (minimum_separation == best_separation &&
             donor_less(candidates[index], candidates[best]))) {
          best = index;
          best_separation = minimum_separation;
        }
      }
      if (best == candidates.size()) break;
      take(best);
    }
  } else {
    for (std::size_t index = 0U;
         index < candidates.size() && out.size() < maximum; ++index) {
      take(index);
    }
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
                      ImmersedDomainBoundaryPolicy boundary_policy,
                      QuadraticStencilLimits limits,
                      const SparseTopologyFluidOracle& topology_oracle,
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
  if (!periodic_policy_valid(boundary_policy, global, reach)) {
    return {StatusCode::invalid_plan, kPlanBoundaryIntersection};
  }
  std::uint8_t available_quadrants = 0U;
  for (int dz = -reach; dz <= reach; ++dz) {
    for (int dy = -reach; dy <= reach; ++dy) {
      for (int dx = -reach; dx <= reach; ++dx) {
        const Int3 raw_index{frame.anchor_global_cell.x + dx,
                             frame.anchor_global_cell.y + dy,
                             frame.anchor_global_cell.z + dz};
        Int3 canonical_index{};
        const std::int32_t raw_components[3]{raw_index.x, raw_index.y,
                                             raw_index.z};
        const std::int32_t extents[3]{global.x, global.y, global.z};
        const double domain_lower[3]{geometry.lower().x, geometry.lower().y,
                                     geometry.lower().z};
        const double domain_upper[3]{geometry.upper().x, geometry.upper().y,
                                     geometry.upper().z};
        double image_shift[3]{};
        bool mapped = true;
        for (int axis = 0; axis < 3; ++axis) {
          std::int32_t canonical = 0;
          double shift = 0.0;
          if (!wrap_periodic_index(raw_components[axis], extents[axis], axis,
                                   boundary_policy, reach, canonical, shift,
                                   domain_lower[axis], domain_upper[axis])) {
            mapped = false;
            break;
          }
          if (axis == 0)
            canonical_index.x = canonical;
          else if (axis == 1)
            canonical_index.y = canonical;
          else
            canonical_index.z = canonical;
          image_shift[axis] = shift;
        }
        if (!mapped || !inside(canonical_index, global)) {
          continue;
        }
        Real3 position{centre(geometry, canonical_index, 0),
                       centre(geometry, canonical_index, 1),
                       centre(geometry, canonical_index, 2)};
        position.x += image_shift[0];
        position.y += image_shift[1];
        position.z += image_shift[2];
        const Real3 delta = subtract(position, origin);
        const double wall_normal = dot(delta, normal);
        if (!(wall_normal >
              64.0 * std::numeric_limits<double>::epsilon() * frame.scale)) {
          continue;
        }
        const CandidateDonor geometric{
            {}, norm_squared(delta), wall_normal / frame.scale,
            dot(delta, frame.tangent1) / frame.scale,
            dot(delta, frame.tangent2) / frame.scale};
        const GlobalCellId canonical_id = global_id(canonical_index, global);
        std::uint8_t sparse_fluid = 0U;
        bool fluid = false;
        if (topology_covers(canonical_index, patch, topology)) {
          fluid = topology.is_fluid_global(canonical_index);
        } else if (topology_oracle.lookup(canonical_id, sparse_fluid)) {
          fluid = sparse_fluid != 0U;
        } else {
          // Every remote canonical donor must have been resolved by the cold
          // sparse topology exchange.  Never guess at an uncovered cell.
          return {StatusCode::invalid_plan, kPlanDonors};
        }
        if (!fluid) {
          continue;
        }
        // Coverage is a donor-side contract, not a Cartesian-search-box
        // contract.  Solid/material-opposite cells cannot certify a
        // tangential quadrant that the reconstruction is unable to use.
        available_quadrants = static_cast<std::uint8_t>(
            available_quadrants |
            (UINT8_C(1) << donor_quadrant(geometric)));
        QuadraticDonorCell donor;
        donor.global_cell = canonical_id;
        // Keep the unwrapped logical search index for the reconstruction
        // reach contract.  The canonical global_cell identifies the unique
        // storage donor, while local_index remains the raw alias coordinate
        // consumed by field/remote exchange; centre is the translated image.
        donor.global_index = raw_index;
        donor.local_index = {raw_index.x - patch.begin.x,
                             raw_index.y - patch.begin.y,
                             raw_index.z - patch.begin.z};
        donor.centre = position;
        donor.widths = {width(geometry, canonical_index, 0),
                        width(geometry, canonical_index, 1),
                        width(geometry, canonical_index, 2)};
        candidates.push_back(
            {donor, norm_squared(delta), wall_normal / frame.scale,
             dot(delta, frame.tangent1) / frame.scale,
             dot(delta, frame.tangent2) / frame.scale});
      }
    }
  }
  if (available_quadrants == 0U) {
    return {StatusCode::invalid_plan, kPlanDonors};
  }
  if (available_quadrants != 0x0fU) {
    const std::array<bool, 6U> near_boundary{
        frame.anchor_global_cell.x < reach,
        frame.anchor_global_cell.x >= global.x - reach,
        frame.anchor_global_cell.y < reach,
        frame.anchor_global_cell.y >= global.y - reach,
        frame.anchor_global_cell.z < reach,
        frame.anchor_global_cell.z >= global.z - reach};
    for (std::size_t face = 0U; face < near_boundary.size(); ++face) {
      if (near_boundary[face] &&
          !boundary_policy.allow_one_sided_quadratic[face] &&
          !boundary_policy.allow_periodic_images[face]) {
        return {StatusCode::invalid_plan, kPlanBoundaryIntersection};
      }
    }
  }
  std::vector<GlobalCellId> canonical_ids;
  canonical_ids.reserve(candidates.size());
  for (const CandidateDonor& candidate : candidates)
    canonical_ids.push_back(candidate.donor.global_cell);
  std::sort(canonical_ids.begin(), canonical_ids.end());
  if (std::adjacent_find(canonical_ids.begin(), canonical_ids.end()) !=
      canonical_ids.end()) {
    return {StatusCode::invalid_plan, kPlanDonors};
  }
  frame.required_quadrant_mask = available_quadrants;
  std::sort(candidates.begin(), candidates.end(), donor_less);
  const std::size_t minimum_donors =
      limits.policy == IbmReconstructionPolicy::adaptive_order
          ? std::min(limits.minimum_donors, limits.minimum_linear_donors)
          : limits.minimum_donors;
  if (candidates.size() < minimum_donors) {
    return {StatusCode::invalid_plan, kPlanDonors};
  }
  const std::size_t count =
      std::min<std::size_t>(candidates.size(), limits.maximum_donors);
  choose_coverage_complete(candidates, count,
                           frame.required_quadrant_mask != 0x0fU, scratch,
                           out);
  std::sort(out.begin(), out.end(), [](const QuadraticDonorCell& left,
                                       const QuadraticDonorCell& right) {
    return left.global_cell < right.global_cell;
  });
  return {};
}

constexpr std::size_t kMetricChannels = 25U;
constexpr std::size_t kMetricArea = 0U;
constexpr std::size_t kMetricAreaVector = 1U;
constexpr std::size_t kMetricFirstMoment = 4U;
constexpr std::size_t kMetricNormalFirstMoment = 7U;
constexpr std::size_t kMetricNormalSecondMoment = 16U;

double component(Real3 value, std::size_t index) noexcept {
  return index == 0U ? value.x : (index == 1U ? value.y : value.z);
}

void add_metric(double* values, Real3 position, Real3 normal,
                double weight) noexcept {
  values[kMetricArea] += weight;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    values[kMetricAreaVector + axis] += weight * component(normal, axis);
    values[kMetricFirstMoment + axis] += weight * component(position, axis);
    for (std::size_t normal_axis = 0U; normal_axis < 3U; ++normal_axis) {
      values[kMetricNormalFirstMoment + 3U * axis + normal_axis] +=
          weight * component(position, axis) *
          component(normal, normal_axis);
      values[kMetricNormalSecondMoment + 3U * axis + normal_axis] +=
          weight * component(normal, axis) *
          component(normal, normal_axis);
    }
  }
}

Status allreduce_sum_doubles(MPI_Comm communicator,
                             std::vector<double>& values) noexcept {
  std::size_t offset = 0U;
  while (offset < values.size()) {
    const std::size_t remaining = values.size() - offset;
    const int count = static_cast<int>(std::min<std::size_t>(
        remaining, static_cast<std::size_t>(std::numeric_limits<int>::max())));
    if (MPI_Allreduce(MPI_IN_PLACE, values.data() + offset, count, MPI_DOUBLE,
                      MPI_SUM, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kPlanCollective};
    }
    offset += static_cast<std::size_t>(count);
  }
  return {};
}

double link_control_area(const CartesianGeometryPlan& geometry,
                         const ImmersedLink& link) noexcept {
  const Int3 cell = link.fluid_global_index;
  const std::size_t direction = static_cast<std::size_t>(link.direction);
  const int axis = static_cast<int>(direction / 2U);
  if (axis == 0)
    return width(geometry, cell, 1) * width(geometry, cell, 2);
  if (axis == 1)
    return width(geometry, cell, 0) * width(geometry, cell, 2);
  return width(geometry, cell, 0) * width(geometry, cell, 1);
}

bool metric_close(double value, double reference,
                  std::uint64_t terms) noexcept {
  const double scale = std::max({1.0, std::abs(value), std::abs(reference)});
  const double tolerance =
      4096.0 * std::numeric_limits<double>::epsilon() *
      static_cast<double>(std::max<std::uint64_t>(terms, 1U)) * scale;
  return std::isfinite(value) && std::isfinite(reference) &&
         std::abs(value - reference) <= tolerance;
}

}  // namespace

void QuadraticStencilPlan::refresh_fingerprint() noexcept {
  Hash64 hash;
  hash.integer(static_cast<std::uint64_t>(groups_.size()));
  hash.integer(static_cast<std::uint64_t>(rows_.size()));
  hash.integer(static_cast<std::uint64_t>(donor_global_cells_.size()));
  hash.integer(static_cast<std::uint64_t>(weights_.size()));
  hash.integer(maximum_halo_reach_);
  hash.integer(audit_.valid ? 1U : 0U);
  hash.integer(static_cast<std::uint8_t>(audit_.policy));
  hash.integer(audit_.standard_reach);
  hash.integer(audit_.group_count);
  hash.integer(audit_.quadratic_groups);
  hash.integer(audit_.linear_groups);
  hash.integer(audit_.expanded_search_groups);
  hash.integer(audit_.rank_fallback_groups);
  hash.integer(audit_.condition_fallback_groups);
  hash.integer(audit_.coverage_fallback_groups);
  hash.integer(audit_.donor_fallback_groups);
  hash.real(audit_.maximum_condition_estimate);
  hash.real(audit_.maximum_functional_l1);
  for (const bool periodic : periodic_axes_) {
    hash.integer(periodic ? 1U : 0U);
  }
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

Status IbmInterfaceMetricCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology, ImmersedPlanLimits limits,
    IbmInterfaceMetricPlan& out) noexcept {
  return compile_with_resident_storage(communicator, geometry, patch, surface,
                                       topology, limits, 0U, out);
}

Status IbmInterfaceMetricCompiler::compile_with_resident_storage(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology, ImmersedPlanLimits limits,
    std::uint64_t resident_persistent_bytes,
    IbmInterfaceMetricPlan& out) noexcept {
  int rank = -1;
  int size = 0;
  const Status context = mpi_context(communicator, rank, size);
  if (!context) return context;

  const Span<const SurfaceTriangle> triangles = surface.triangles();
  const Span<const std::array<std::uint32_t, 3U>> neighbours =
      surface.triangle_neighbours();
  const Span<const ImmersedLink> links = topology.links();
  Status local{};
  if (geometry.fingerprint() == 0U || surface.fingerprint() == 0U ||
      topology.fingerprint() == 0U ||
      topology.geometry_revision() != geometry.topology_revision() ||
      topology.geometry_fingerprint() != geometry.fingerprint() ||
      topology.surface_fingerprint() != surface.fingerprint() ||
      !same(patch, topology.patch_) ||
      !valid_patch(geometry.global_cells(), patch, rank, size) ||
      triangles.size == 0U || neighbours.size != triangles.size ||
      triangles.size > static_cast<std::size_t>(UINT32_MAX) ||
      links.size > limits.maximum_local_links ||
      limits.maximum_persistent_bytes_per_rank == 0U ||
      limits.maximum_peak_bytes_per_rank <
          limits.maximum_persistent_bytes_per_rank ||
      resident_persistent_bytes >
          limits.maximum_persistent_bytes_per_rank ||
      resident_persistent_bytes > limits.maximum_peak_bytes_per_rank) {
    local = {StatusCode::invalid_plan, kPlanInput};
  }
  for (std::size_t triangle = 0U; local && triangle < neighbours.size;
       ++triangle) {
    for (const std::uint32_t adjacent : neighbours.data[triangle]) {
      if (adjacent >= triangles.size || adjacent == triangle) {
        local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
        break;
      }
    }
  }
  int lowest = -1;
  Status agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  agreed = agree_compile_contract(
      communicator, rank,
      shared_compile_contract(geometry, surface, topology, patch,
                              ImmersedDomainBoundaryPolicy{}, limits, size,
                              3U),
      lowest);
  if (!agreed) return agreed;

  const std::uint64_t triangle_count =
      static_cast<std::uint64_t>(triangles.size);
  const std::uint64_t local_link_count =
      static_cast<std::uint64_t>(links.size);
  std::uint64_t metric_values = 0U;
  std::uint64_t persistent_bytes = 0U;
  std::uint64_t control_bytes = 0U;
  std::uint64_t metric_bytes = 0U;
  std::uint64_t paired_metric_bytes = 0U;
  std::uint64_t index_bytes = 0U;
  std::uint64_t audit_bytes =
      static_cast<std::uint64_t>((kMetricChannels + 1U) * sizeof(double));
  std::uint64_t temporary_bytes = 0U;
  std::uint64_t peak_bytes = 0U;
  std::uint64_t combined_persistent_bytes = 0U;
  std::uint64_t combined_peak_bytes = 0U;
  if (!checked_multiply(triangle_count, kMetricChannels, metric_values) ||
      !checked_multiply(local_link_count, sizeof(IbmInterfaceLinkMetric),
                        persistent_bytes) ||
      !checked_multiply(triangle_count, sizeof(double), control_bytes) ||
      !checked_multiply(metric_values, sizeof(double), metric_bytes) ||
      !checked_multiply(metric_bytes, 2U, paired_metric_bytes) ||
      !checked_multiply(triangle_count, 3U * sizeof(std::uint32_t),
                        index_bytes) ||
      !checked_add(control_bytes, paired_metric_bytes, temporary_bytes) ||
      !checked_add(temporary_bytes, index_bytes, temporary_bytes) ||
      !checked_add(temporary_bytes, audit_bytes, temporary_bytes) ||
      !checked_add(persistent_bytes, temporary_bytes, peak_bytes) ||
      !checked_add(resident_persistent_bytes, persistent_bytes,
                   combined_persistent_bytes) ||
      !checked_add(resident_persistent_bytes, peak_bytes,
                   combined_peak_bytes) ||
      combined_persistent_bytes >
          limits.maximum_persistent_bytes_per_rank ||
      combined_peak_bytes > limits.maximum_peak_bytes_per_rank ||
      metric_values > std::numeric_limits<std::size_t>::max()) {
    local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  IbmInterfaceMetricPlan candidate;
  std::vector<double> control;
  std::vector<double> triangle_metric;
  std::vector<double> source_metric;
  std::vector<std::uint32_t> distance;
  std::vector<std::uint32_t> source;
  std::vector<std::uint32_t> queue;
  std::vector<double> audit;
  try {
    candidate.links_.reserve(links.size);
    control.assign(triangles.size, 0.0);
    triangle_metric.assign(static_cast<std::size_t>(metric_values), 0.0);
    source_metric.assign(static_cast<std::size_t>(metric_values), 0.0);
    distance.assign(triangles.size, UINT32_MAX);
    source.assign(triangles.size, kInvalidIbmIndex);
    queue.reserve(triangles.size);
    audit.assign(kMetricChannels + 1U, 0.0);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kPlanInterfaceMetric};
  } catch (...) {
    local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  const auto vector_bytes = [](std::size_t capacity, std::size_t element,
                               std::uint64_t& bytes) noexcept {
    return checked_multiply(static_cast<std::uint64_t>(capacity),
                            static_cast<std::uint64_t>(element), bytes);
  };
  std::array<std::uint64_t, 7U> scratch_bytes{};
  temporary_bytes = 0U;
  if (!vector_bytes(candidate.links_.capacity(),
                    sizeof(IbmInterfaceLinkMetric), persistent_bytes) ||
      !vector_bytes(control.capacity(), sizeof(double), scratch_bytes[0U]) ||
      !vector_bytes(triangle_metric.capacity(), sizeof(double),
                    scratch_bytes[1U]) ||
      !vector_bytes(source_metric.capacity(), sizeof(double),
                    scratch_bytes[2U]) ||
      !vector_bytes(distance.capacity(), sizeof(std::uint32_t),
                    scratch_bytes[3U]) ||
      !vector_bytes(source.capacity(), sizeof(std::uint32_t),
                    scratch_bytes[4U]) ||
      !vector_bytes(queue.capacity(), sizeof(std::uint32_t),
                    scratch_bytes[5U]) ||
      !vector_bytes(audit.capacity(), sizeof(double), scratch_bytes[6U])) {
    local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  }
  for (const std::uint64_t bytes : scratch_bytes) {
    if (local && !checked_add(temporary_bytes, bytes, temporary_bytes))
      local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  }
  if (local &&
      (!checked_add(persistent_bytes, temporary_bytes, peak_bytes) ||
       !checked_add(resident_persistent_bytes, persistent_bytes,
                    combined_persistent_bytes) ||
       !checked_add(resident_persistent_bytes, peak_bytes,
                    combined_peak_bytes) ||
       combined_persistent_bytes >
           limits.maximum_persistent_bytes_per_rank ||
       combined_peak_bytes > limits.maximum_peak_bytes_per_rank)) {
    local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  std::uint64_t global_link_count = local_link_count;
  if (MPI_Allreduce(MPI_IN_PLACE, &global_link_count, 1, MPI_UINT64_T, MPI_SUM,
                    communicator) != MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kPlanCollective};
  }
  for (std::size_t index = 0U; local && index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    const double expected = link_control_area(geometry, link);
    if (link.triangle >= triangles.size ||
        !std::isfinite(link.cartesian_control_face_area) ||
        !(link.cartesian_control_face_area > 0.0) ||
        link.cartesian_control_face_area != expected) {
      local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
      break;
    }
    control.data()[static_cast<std::size_t>(link.triangle)] +=
        link.cartesian_control_face_area;
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  local = allreduce_sum_doubles(communicator, control);
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  local = visit_domain_clipped_quadrature(
      geometry, triangles,
      [&](const SurfaceTriangle& triangle, std::uint8_t, Real3 point,
          double weight) noexcept -> Status {
        const Int3 execution_anchor = anchor_cell(geometry, point);
        if (owner_rank(execution_anchor, geometry.global_cells(),
                       patch.process_grid) != rank) {
          return {};
        }
        const Real3 normal =
            topology.fluid_side_ == ImmersedFluidSide::outside
                ? triangle.geometric_outward_normal
                : Real3{-triangle.geometric_outward_normal.x,
                        -triangle.geometric_outward_normal.y,
                        -triangle.geometric_outward_normal.z};
        add_metric(triangle_metric.data() +
                       static_cast<std::size_t>(triangle.id) *
                           kMetricChannels,
                   point, normal, weight);
        return {};
      });
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  local = allreduce_sum_doubles(communicator, triangle_metric);
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  for (std::uint32_t triangle = 0U; triangle < triangles.size; ++triangle) {
    if (control[triangle] > 0.0) {
      distance[triangle] = 0U;
      source[triangle] = triangle;
      queue.push_back(triangle);
    }
  }
  for (std::size_t head = 0U; head < queue.size(); ++head) {
    const std::uint32_t triangle = queue[head];
    for (const std::uint32_t adjacent : neighbours.data[triangle]) {
      if (distance[adjacent] == UINT32_MAX) {
        distance[adjacent] = distance[triangle] + 1U;
        queue.push_back(adjacent);
      }
    }
  }
  // The first pass finds only graph distance.  This second pass selects the
  // globally deterministic minimum source among all shortest paths.
  for (const std::uint32_t triangle : queue) {
    if (distance[triangle] == 0U) continue;
    std::uint32_t selected = kInvalidIbmIndex;
    for (const std::uint32_t adjacent : neighbours.data[triangle]) {
      if (distance[adjacent] + 1U == distance[triangle])
        selected = std::min(selected, source[adjacent]);
    }
    if (selected == kInvalidIbmIndex) {
      local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
      break;
    }
    source[triangle] = selected;
  }
  for (std::size_t triangle = 0U; local && triangle < triangles.size;
       ++triangle) {
    const double area =
        triangle_metric[triangle * kMetricChannels + kMetricArea];
    if (source[triangle] == kInvalidIbmIndex) {
      if (area != 0.0)
        local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
      continue;
    }
    double* const target =
        source_metric.data() +
        static_cast<std::size_t>(source[triangle]) * kMetricChannels;
    const double* const contribution =
        triangle_metric.data() + triangle * kMetricChannels;
    for (std::size_t channel = 0U; channel < kMetricChannels; ++channel)
      target[channel] += contribution[channel];
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  const auto real3_from = [](const double* values) noexcept {
    return Real3{values[0U], values[1U], values[2U]};
  };
  for (std::size_t index = 0U; local && index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    const std::size_t seed = static_cast<std::size_t>(link.triangle);
    const double denominator = control[seed];
    const double ratio = link.cartesian_control_face_area / denominator;
    const double* const assigned =
        source_metric.data() + seed * kMetricChannels;
    if (!std::isfinite(denominator) || !(denominator > 0.0) ||
        !std::isfinite(ratio) || !(ratio > 0.0) ||
        !std::isfinite(assigned[kMetricArea]) ||
        !(assigned[kMetricArea] > 0.0)) {
      local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
      break;
    }
    IbmInterfaceLinkMetric metric;
    metric.global_link = link.global_link;
    metric.source_triangle = link.triangle;
    metric.physical_quadrature_area = ratio * assigned[kMetricArea];
    std::array<double, 3U> area_vector{};
    std::array<double, 3U> first_moment{};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      area_vector[axis] = ratio * assigned[kMetricAreaVector + axis];
      first_moment[axis] = ratio * assigned[kMetricFirstMoment + axis];
    }
    metric.physical_area_vector = real3_from(area_vector.data());
    metric.physical_first_moment = real3_from(first_moment.data());
    for (std::size_t entry = 0U; entry < 9U; ++entry) {
      metric.normal_first_moment[entry] =
          ratio * assigned[kMetricNormalFirstMoment + entry];
      metric.normal_second_moment[entry] =
          ratio * assigned[kMetricNormalSecondMoment + entry];
    }
    candidate.links_.push_back(metric);
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  for (std::size_t index = 0U; index < candidate.links_.size(); ++index) {
    const IbmInterfaceLinkMetric& metric = candidate.links_[index];
    audit[0U] += links.data[index].cartesian_control_face_area;
    audit[1U + kMetricArea] += metric.physical_quadrature_area;
    audit[1U + kMetricAreaVector + 0U] += metric.physical_area_vector.x;
    audit[1U + kMetricAreaVector + 1U] += metric.physical_area_vector.y;
    audit[1U + kMetricAreaVector + 2U] += metric.physical_area_vector.z;
    audit[1U + kMetricFirstMoment + 0U] += metric.physical_first_moment.x;
    audit[1U + kMetricFirstMoment + 1U] += metric.physical_first_moment.y;
    audit[1U + kMetricFirstMoment + 2U] += metric.physical_first_moment.z;
    for (std::size_t entry = 0U; entry < 9U; ++entry) {
      audit[1U + kMetricNormalFirstMoment + entry] +=
          metric.normal_first_moment[entry];
      audit[1U + kMetricNormalSecondMoment + entry] +=
          metric.normal_second_moment[entry];
    }
  }
  local = allreduce_sum_doubles(communicator, audit);
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;
  std::array<double, kMetricChannels> total{};
  for (std::size_t triangle = 0U; triangle < triangles.size; ++triangle)
    for (std::size_t channel = 0U; channel < kMetricChannels; ++channel)
      total[channel] +=
          triangle_metric[triangle * kMetricChannels + channel];
  double total_control = 0.0;
  for (const double value : control) total_control += value;
  const std::uint64_t audit_terms =
      std::max(global_link_count, triangle_count);
  if (!metric_close(audit[0U], total_control, audit_terms))
    local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  for (std::size_t channel = 0U; local && channel < kMetricChannels;
       ++channel) {
    if (!metric_close(audit[1U + channel], total[channel], audit_terms))
      local = {StatusCode::invalid_plan, kPlanInterfaceMetric};
  }
  agreed = consensus(communicator, rank, local, lowest);
  if (!agreed) return agreed;

  candidate.conservation_.cartesian_control_area = total_control;
  candidate.conservation_.physical_quadrature_area = total[kMetricArea];
  candidate.conservation_.physical_area_vector =
      real3_from(total.data() + kMetricAreaVector);
  candidate.conservation_.physical_first_moment =
      real3_from(total.data() + kMetricFirstMoment);
  for (std::size_t entry = 0U; entry < 9U; ++entry) {
    candidate.conservation_.normal_first_moment[entry] =
        total[kMetricNormalFirstMoment + entry];
    candidate.conservation_.normal_second_moment[entry] =
        total[kMetricNormalSecondMoment + entry];
  }
  candidate.resources_.persistent_bytes_per_rank = persistent_bytes;
  candidate.resources_.peak_bytes_per_rank = peak_bytes;
  candidate.resources_.collective_doubles_per_rank =
      (1U + kMetricChannels) * triangle_count + kMetricChannels + 1U;
  candidate.geometry_revision_ = geometry.topology_revision();
  candidate.geometry_fingerprint_ = geometry.fingerprint();
  candidate.surface_fingerprint_ = surface.fingerprint();
  Hash64 physical_hash;
  physical_hash.integer(kInterfaceMetricSchemaRevision);
  physical_hash.integer(geometry.fingerprint());
  physical_hash.integer(surface.fingerprint());
  physical_hash.integer(static_cast<std::uint8_t>(topology.fluid_side_));
  physical_hash.integer(global_link_count);
  candidate.physical_fingerprint_ = physical_hash.finish();
  Hash64 layout_hash;
  layout_hash.integer(candidate.physical_fingerprint_);
  layout_hash.integer(rank);
  layout_hash.integer(static_cast<std::uint64_t>(candidate.links_.size()));
  for (const IbmInterfaceLinkMetric& metric : candidate.links_) {
    layout_hash.integer(metric.global_link);
    layout_hash.integer(metric.source_triangle);
    layout_hash.real(metric.physical_quadrature_area);
    layout_hash.real(metric.physical_area_vector.x);
    layout_hash.real(metric.physical_area_vector.y);
    layout_hash.real(metric.physical_area_vector.z);
    layout_hash.real(metric.physical_first_moment.x);
    layout_hash.real(metric.physical_first_moment.y);
    layout_hash.real(metric.physical_first_moment.z);
    for (const double value : metric.normal_first_moment)
      layout_hash.real(value);
    for (const double value : metric.normal_second_moment)
      layout_hash.real(value);
  }
  candidate.fingerprint_ = layout_hash.finish();
  out = std::move(candidate);
  return {};
}

Status BoundaryStencilCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology, ImmersedPlanLimits limits,
    BoundaryStencilPlan& out) noexcept {
  return compile(communicator, geometry, patch, surface, topology,
                 ImmersedDomainBoundaryPolicy{}, limits, out);
}

Status BoundaryStencilCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology,
    ImmersedDomainBoundaryPolicy boundary_policy, ImmersedPlanLimits limits,
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
        !periodic_policy_valid(
            boundary_policy, geometry.global_cells(),
            static_cast<int>(limits.stencil.maximum_reach)) ||
        limits.maximum_persistent_bytes_per_rank == 0U ||
        limits.maximum_peak_bytes_per_rank <
            limits.maximum_persistent_bytes_per_rank) {
      local = {StatusCode::invalid_plan, kPlanInput};
    }
    for (std::size_t face = 0U; face < boundary_policy.allow_periodic_images.size();
         ++face) {
      if (boundary_policy.allow_periodic_images[face] &&
          boundary_policy.allow_one_sided_quadratic[face]) {
        local = {StatusCode::invalid_plan, kPlanInput};
      }
    }
    int lowest = -1;
    Status agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    agreed = agree_compile_contract(
        communicator, rank,
        shared_compile_contract(geometry, surface, topology, patch,
                                boundary_policy, limits, size, 1U),
        lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    // Resolve every donor cell that lies outside this rank's topology halo
    // before fitting any stencil.  This is a cold compile-time exchange;
    // no sparse topology traffic is introduced into the timestep path.
    SparseTopologyFluidOracle topology_oracle;
    std::size_t link_count = 0U;
    std::uint64_t planned_peak_bytes = 0U;
    if (local) {
      link_count = topology.links_.size();
      std::uint64_t row_count = 0U;
      std::uint64_t donor_count = 0U;
      std::uint64_t weight_count = 0U;
      if (!worst_case_plan_counts(link_count, 4U, limits.stencil, row_count,
                                  donor_count, weight_count) ||
          !plan_memory_within_limits(
              link_count, row_count, donor_count, weight_count, link_count,
              sizeof(BoundaryStencilLink), limits, &planned_peak_bytes)) {
        local = {StatusCode::invalid_plan, kPlanDonors};
      }
    }
    agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    const std::uint64_t query_budget =
        planned_peak_bytes < limits.maximum_peak_bytes_per_rank
            ? limits.maximum_peak_bytes_per_rank - planned_peak_bytes
            : 0U;
    local = prepare_topology_fluid_oracle_from_links(
        communicator, rank, size, geometry, patch, topology,
        boundary_policy, limits.stencil, query_budget, topology_oracle);
    if (local) {
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
                               link.solid_to_fluid_normal, boundary_policy,
                               limits.stencil, topology_oracle, donor_scratch,
                               frame, donors);
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
                   QuadraticConstraint::origin_value, link.wall_point,
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
      if (local && prepared.empty()) {
        candidate.reconstruction_.audit_.valid = true;
        candidate.reconstruction_.audit_.policy = limits.stencil.policy;
        candidate.reconstruction_.audit_.standard_reach =
            limits.stencil.standard_reach;
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
      for (int axis = 0; axis < 3; ++axis)
        candidate.reconstruction_.periodic_axes_[static_cast<std::size_t>(axis)] =
            periodic_pair(boundary_policy, axis);
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
  return compile(communicator, geometry, patch, surface, topology,
                 ImmersedDomainBoundaryPolicy{}, limits, out);
}

Status SurfaceQuadratureCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const ImmersedSurfacePlan& surface,
    const EBTopology& topology,
    ImmersedDomainBoundaryPolicy boundary_policy, ImmersedPlanLimits limits,
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
        !periodic_policy_valid(
            boundary_policy, geometry.global_cells(),
            static_cast<int>(limits.stencil.maximum_reach)) ||
        limits.maximum_local_quadrature_points == 0U ||
        limits.maximum_persistent_bytes_per_rank == 0U ||
        limits.maximum_peak_bytes_per_rank <
            limits.maximum_persistent_bytes_per_rank) {
      local = {StatusCode::invalid_plan, kPlanInput};
    }
    for (std::size_t face = 0U; face < boundary_policy.allow_periodic_images.size();
         ++face) {
      if (boundary_policy.allow_periodic_images[face] &&
          boundary_policy.allow_one_sided_quadratic[face]) {
        local = {StatusCode::invalid_plan, kPlanInput};
      }
    }
    int lowest = -1;
    Status agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    agreed = agree_compile_contract(
        communicator, rank,
        shared_compile_contract(geometry, surface, topology, patch,
                                boundary_policy, limits, size, 2U),
        lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    const Span<const SurfaceTriangle> surface_triangles = surface.triangles();
    std::uint64_t local_point_count = 0U;
    std::uint64_t planned_peak_bytes = 0U;
    std::vector<GlobalCellId> topology_query_ids;
    SparseTopologyFluidOracle topology_oracle;
    if (local) {
      // Clip the replicated immutable STL to the open Cartesian domain and
      // count ownership before allocating.  Domain-exterior closure facets
      // are not immersed traction surface.
      local = visit_domain_clipped_quadrature(
          geometry, surface_triangles,
          [&](const SurfaceTriangle&, std::uint8_t, Real3 point,
              double) noexcept -> Status {
            if (!checked_add(candidate.global_point_count_, 1U,
                             candidate.global_point_count_))
              return {StatusCode::invalid_plan, kPlanQuadrature};
            const Int3 execution_anchor = anchor_cell(geometry, point);
            const int selected_rank = owner_rank(execution_anchor,
                                                 geometry.global_cells(),
                                                 patch.process_grid);
            if (selected_rank == rank) {
              if (!checked_add(local_point_count, 1U, local_point_count))
                return {StatusCode::invalid_plan, kPlanQuadrature};
            }
            return {};
          });
      if (local && candidate.global_point_count_ == 0U)
        local = {StatusCode::invalid_plan, kPlanQuadrature};
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
              local_point_count, sizeof(SurfaceQuadraturePoint), limits,
              &planned_peak_bytes)) {
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
    const std::uint64_t query_budget =
        planned_peak_bytes < limits.maximum_peak_bytes_per_rank
            ? limits.maximum_peak_bytes_per_rank - planned_peak_bytes
            : 0U;
    if (local) {
      // Enumerate the exact local quadrature points only after the normal
      // plan-memory gate has reserved its peak budget.  The sparse topology
      // request/reply/lookup storage then consumes the remaining budget.
      local = visit_domain_clipped_quadrature(
          geometry, surface_triangles,
          [&](const SurfaceTriangle&, std::uint8_t, Real3 point,
              double) noexcept -> Status {
            const Int3 execution_anchor = anchor_cell(geometry, point);
            const int selected_rank = owner_rank(execution_anchor,
                                                 geometry.global_cells(),
                                                 patch.process_grid);
            if (selected_rank != rank) return {};
            try {
              return append_query_ids_for_point(
                  geometry, patch, topology, point, boundary_policy,
                  limits.stencil, query_budget, topology_query_ids);
            } catch (const std::bad_alloc&) {
              return {StatusCode::allocation_failure, kPlanDonors};
            } catch (...) {
              return {StatusCode::invalid_plan, kPlanDonors};
            }
          });
    }
    agreed = consensus(communicator, rank, local, lowest);
    if (!agreed) {
      out.lowest_failing_rank_ = lowest;
      return agreed;
    }
    if (local) {
      std::sort(topology_query_ids.begin(), topology_query_ids.end());
      topology_query_ids.erase(
          std::unique(topology_query_ids.begin(), topology_query_ids.end()),
          topology_query_ids.end());
      local = query_topology_fluid(
          communicator, rank, size, geometry.global_cells(), patch, topology,
          topology_query_ids, query_budget, topology_oracle);
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
      local = visit_domain_clipped_quadrature(
          geometry, surface_triangles,
          [&](const SurfaceTriangle& triangle, std::uint8_t point_index,
              Real3 point, double point_weight) -> Status {
            const Real3 solid_to_fluid =
                topology.fluid_side_ == ImmersedFluidSide::outside
                    ? triangle.geometric_outward_normal
                    : Real3{-triangle.geometric_outward_normal.x,
                            -triangle.geometric_outward_normal.y,
                            -triangle.geometric_outward_normal.z};
            const Int3 execution_anchor = anchor_cell(geometry, point);
            const int selected_rank = owner_rank(execution_anchor,
                                                 geometry.global_cells(),
                                                 patch.process_grid);
            if (selected_rank != rank) {
              return {};
            }
            if (candidate.local_points_.size() >= exact_local_count) {
              return {StatusCode::invalid_plan, kPlanQuadrature};
            }

            QuadraticFrame frame;
            const Status donor_status = collect_donors(
                geometry, patch, topology, point, solid_to_fluid,
                boundary_policy, limits.stencil, topology_oracle,
                donor_scratch, frame, donors);
            if (!donor_status) return donor_status;
            GlobalCellId owner_cell = kInvalidSurfaceTriangle;
            double owner_distance = std::numeric_limits<double>::infinity();
            for (const QuadraticDonorCell& donor : donors) {
              const double distance =
                  norm_squared(subtract(donor.centre, point));
              if (distance < owner_distance ||
                  (distance == owner_distance &&
                   donor.global_cell < owner_cell)) {
                owner_distance = distance;
                owner_cell = donor.global_cell;
              }
            }
            // Pressure and other unknown wall values are reconstructed from
            // fluid-side donors, never a caller-supplied wall value.
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
                {triangle.id, point_index, point, solid_to_fluid, point_weight,
                 rank, owner_cell, kInvalidIbmIndex, kInvalidIbmIndex,
                 kInvalidIbmIndex});
            return {};
          });
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
      if (local && prepared.empty()) {
        candidate.reconstruction_.audit_.valid = true;
        candidate.reconstruction_.audit_.policy = limits.stencil.policy;
        candidate.reconstruction_.audit_.standard_reach =
            limits.stencil.standard_reach;
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
      for (int axis = 0; axis < 3; ++axis)
        candidate.reconstruction_.periodic_axes_[static_cast<std::size_t>(axis)] =
            periodic_pair(boundary_policy, axis);
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
  physical.integer(geometry.fingerprint());
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
