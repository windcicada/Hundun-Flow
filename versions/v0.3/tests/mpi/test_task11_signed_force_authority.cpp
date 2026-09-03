// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_immersed_operator.hpp"

#include "tests/support/fvm_immersed_operator_test_access.hpp"
#include "tests/support/flow_immersed_test_access.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "src/ib_wall_force_detail.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using namespace hundun;

constexpr runtime::Int3 kExtent{12, 12, 12};
constexpr runtime::PhaseId kPhase = 111U;
constexpr runtime::ActorId kActor = 111U;
constexpr double kPressure = 2.75;
constexpr double kViscosity = 0.35;
constexpr double kArithmeticFactor = 262144.0;
constexpr int kWallReach = 4;
constexpr immersed::ImmersedLinkId kMixedAuthorityLink = 77U;

runtime::Real3 add(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

runtime::Real3 subtract(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

runtime::Real3 multiply(double scale, runtime::Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

runtime::Real3 cross(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double norm(runtime::Real3 value) noexcept {
  return std::sqrt(dot(value, value));
}

double max_abs(runtime::Real3 value) noexcept {
  return std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

runtime::Real3 normalized(runtime::Real3 value) {
  const double magnitude = norm(value);
  if (!(magnitude > 0.0) || !std::isfinite(magnitude))
    throw runtime::Error("Task 11 RED-S1 normalization is invalid");
  return multiply(1.0 / magnitude, value);
}

double arithmetic_bound(runtime::Real3 actual,
                        runtime::Real3 expected) noexcept {
  return kArithmeticFactor * std::numeric_limits<double>::epsilon() *
         std::max({1.0, max_abs(actual), max_abs(expected)});
}

bool near(runtime::Real3 actual, runtime::Real3 expected) noexcept {
  const double tolerance = arithmetic_bound(actual, expected);
  return std::abs(actual.x - expected.x) <= tolerance &&
         std::abs(actual.y - expected.y) <= tolerance &&
         std::abs(actual.z - expected.z) <= tolerance;
}

bool near(double actual, double expected) noexcept {
  const double tolerance =
      kArithmeticFactor * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(actual), std::abs(expected)});
  return std::abs(actual - expected) <= tolerance;
}

immersed::ForceComponents force(runtime::Real3 pressure,
                                runtime::Real3 viscous) noexcept {
  return {pressure, viscous, add(pressure, viscous)};
}

immersed::ForceComponents add_force(const immersed::ForceComponents &left,
                                    const immersed::ForceComponents &right) {
  return force(add(left.pressure_N, right.pressure_N),
               add(left.viscous_N, right.viscous_N));
}

immersed::ForceComponents subtract_force(
    const immersed::ForceComponents &left,
    const immersed::ForceComponents &right) {
  return force(subtract(left.pressure_N, right.pressure_N),
               subtract(left.viscous_N, right.viscous_N));
}

immersed::ForceComponents negate_force(
    const immersed::ForceComponents &value) {
  return force(multiply(-1.0, value.pressure_N),
               multiply(-1.0, value.viscous_N));
}

bool near_force(const immersed::ForceComponents &actual,
                const immersed::ForceComponents &expected) noexcept {
  return near(actual.pressure_N, expected.pressure_N) &&
         near(actual.viscous_N, expected.viscous_N) &&
         near(actual.total_N, expected.total_N);
}

void check_mpi(int code, const char *operation) {
  if (code != MPI_SUCCESS)
    throw runtime::Error(std::string(operation) + " failed");
}

template <class Value>
void allreduce_sum(std::vector<Value> &values, MPI_Datatype datatype,
                   const runtime::MpiContext &mpi, const char *operation) {
  if (values.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error(std::string(operation) + " count overflows");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, values.data(),
                          static_cast<int>(values.size()), datatype, MPI_SUM,
                          mpi.comm()),
            operation);
}

bool all_true(bool local, const runtime::MpiContext &mpi,
              const char *operation) {
  int value = local ? 1 : 0;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_MIN,
                          mpi.comm()),
            operation);
  return value == 1;
}

std::array<double, 9> outer(runtime::Real3 left,
                            runtime::Real3 right) noexcept {
  return {left.x * right.x, left.x * right.y, left.x * right.z,
          left.y * right.x, left.y * right.y, left.y * right.z,
          left.z * right.x, left.z * right.y, left.z * right.z};
}

runtime::Real3 multiply(const std::array<double, 9> &matrix,
                        runtime::Real3 vector) noexcept {
  return {matrix[0] * vector.x + matrix[1] * vector.y +
              matrix[2] * vector.z,
          matrix[3] * vector.x + matrix[4] * vector.y +
              matrix[5] * vector.z,
          matrix[6] * vector.x + matrix[7] * vector.y +
              matrix[8] * vector.z};
}

struct AnalyticStress final {
  runtime::Real3 tangent{};
  runtime::Real3 a{};
  std::array<double, 9> gradient{};
  std::array<double, 9> tau{};
};

AnalyticStress analytic_stress(runtime::Real3 normal) {
  const std::array<double, 3> absolute{std::abs(normal.x), std::abs(normal.y),
                                      std::abs(normal.z)};
  std::size_t least = 0U;
  for (std::size_t axis = 1U; axis < absolute.size(); ++axis)
    if (absolute[axis] < absolute[least])
      least = axis;
  const runtime::Real3 basis = least == 0U   ? runtime::Real3{1.0, 0.0, 0.0}
                               : least == 1U ? runtime::Real3{0.0, 1.0, 0.0}
                                             : runtime::Real3{0.0, 0.0, 1.0};
  AnalyticStress result;
  result.tangent = normalized(cross(basis, normal));
  result.a = add(multiply(0.70, result.tangent), multiply(0.40, normal));
  result.gradient = outer(result.a, normal);
  const double divergence =
      result.gradient[0] + result.gradient[4] + result.gradient[8];
  for (std::size_t row = 0U; row < 3U; ++row)
    for (std::size_t column = 0U; column < 3U; ++column)
      result.tau[row * 3U + column] =
          kViscosity *
          (result.gradient[row * 3U + column] +
           result.gradient[column * 3U + row] -
           (row == column ? (2.0 / 3.0) * divergence : 0.0));
  return result;
}

test::StlFixtureTriangle triangle_from_vertices(
    std::array<runtime::Real3, 3> vertices) {
  const auto area_vector =
      cross(subtract(vertices[1], vertices[0]),
            subtract(vertices[2], vertices[0]));
  const auto normal = normalized(area_vector);
  return {normal, std::move(vertices)};
}

runtime::Real3 midpoint(runtime::Real3 left, runtime::Real3 right) noexcept {
  return multiply(0.5, add(left, right));
}

std::vector<test::StlFixtureTriangle> signed_force_intermediate_triangles() {
  auto triangles = test::outward_cube();
  for (int level = 0; level < 2; ++level) {
    std::vector<test::StlFixtureTriangle> refined;
    refined.reserve(4U * triangles.size());
    for (const auto &triangle : triangles) {
      const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
      const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
      const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
      refined.push_back(triangle_from_vertices(
          {triangle.vertices[0], ab, ca}));
      refined.push_back(triangle_from_vertices(
          {ab, triangle.vertices[1], bc}));
      refined.push_back(triangle_from_vertices(
          {ca, bc, triangle.vertices[2]}));
      refined.push_back(triangle_from_vertices({ab, bc, ca}));
    }
    triangles = std::move(refined);
  }
  constexpr runtime::Real3 q1{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
  constexpr runtime::Real3 q2{2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0};
  constexpr runtime::Real3 q3{-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0};
  for (auto &triangle : triangles) {
    std::array<runtime::Real3, 3> transformed{};
    for (std::size_t vertex = 0U; vertex < transformed.size(); ++vertex) {
      const auto centered =
          subtract(triangle.vertices[vertex], {0.5, 0.5, 0.5});
      const auto rotated =
          add(add(multiply(centered.x, q1), multiply(centered.y, q2)),
              multiply(centered.z, q3));
      transformed[vertex] =
          add({0.5, 0.5, 0.5}, multiply(0.36, rotated));
    }
    triangle = triangle_from_vertices(transformed);
  }
  return triangles;
}

double coordinate(runtime::Real3 point, std::size_t axis) noexcept {
  if (axis == 0U)
    return point.x;
  if (axis == 1U)
    return point.y;
  return point.z;
}

runtime::Int3 logical_cell(mesh::GlobalCellId id) {
  const auto nx = static_cast<std::uint64_t>(kExtent.x);
  const auto ny = static_cast<std::uint64_t>(kExtent.y);
  const auto plane = nx * ny;
  if (id >= plane * static_cast<std::uint64_t>(kExtent.z))
    throw runtime::Error("Task 11 RED-S1 donor global cell is invalid");
  return {static_cast<int>(id % nx), static_cast<int>((id / nx) % ny),
          static_cast<int>(id / plane)};
}

bool inside_expanded(runtime::Box3 box, runtime::Int3 cell,
                     int reach) noexcept {
  return cell.x >= box.begin.x - reach && cell.x < box.end.x + reach &&
         cell.y >= box.begin.y - reach && cell.y < box.end.y + reach &&
         cell.z >= box.begin.z - reach && cell.z < box.end.z + reach;
}

void set_coordinate(runtime::Real3 &point, std::size_t axis,
                    double value) noexcept {
  if (axis == 0U)
    point.x = value;
  else if (axis == 1U)
    point.y = value;
  else
    point.z = value;
}

using VertexBits = std::array<std::uint64_t, 3>;

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t result{};
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

VertexBits vertex_bits(runtime::Real3 point) noexcept {
  return {double_bits(point.x), double_bits(point.y), double_bits(point.z)};
}

bool bitwise_equal(runtime::Real3 left, runtime::Real3 right) noexcept {
  return vertex_bits(left) == vertex_bits(right);
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

runtime::Real3 oriented_area_vector(
    const test::StlFixtureTriangle &triangle) noexcept {
  return multiply(0.5,
                  cross(subtract(triangle.vertices[1], triangle.vertices[0]),
                        subtract(triangle.vertices[2], triangle.vertices[0])));
}

double scalar_area(const test::StlFixtureTriangle &triangle) noexcept {
  return norm(oriented_area_vector(triangle));
}

runtime::Real3 canonical_plane_intersection(runtime::Real3 first,
                                            runtime::Real3 second,
                                            std::size_t axis) {
  constexpr double plane = 0.5;
  const double first_coordinate = coordinate(first, axis);
  const double second_coordinate = coordinate(second, axis);
  if (first_coordinate == plane) {
    set_coordinate(first, axis, plane);
    return first;
  }
  if (second_coordinate == plane) {
    set_coordinate(second, axis, plane);
    return second;
  }
  const auto lower = first_coordinate < second_coordinate ? first : second;
  const auto higher = first_coordinate < second_coordinate ? second : first;
  const double lower_coordinate = coordinate(lower, axis);
  const double higher_coordinate = coordinate(higher, axis);
  if (!(lower_coordinate < plane && higher_coordinate > plane))
    throw runtime::Error("Task 11 RED-S1 split edge does not cross plane");
  const double fraction =
      (plane - lower_coordinate) / (higher_coordinate - lower_coordinate);
  auto intersection =
      add(lower, multiply(fraction, subtract(higher, lower)));
  set_coordinate(intersection, axis, plane);
  return intersection;
}

std::vector<runtime::Real3>
remove_adjacent_duplicate_vertices(std::vector<runtime::Real3> polygon) {
  std::vector<runtime::Real3> result;
  result.reserve(polygon.size());
  for (const auto vertex : polygon)
    if (result.empty() || !bitwise_equal(result.back(), vertex))
      result.push_back(vertex);
  if (result.size() > 1U && bitwise_equal(result.front(), result.back()))
    result.pop_back();
  return result;
}

std::vector<runtime::Real3>
clip_polygon(const std::array<runtime::Real3, 3> &vertices,
             std::size_t axis, bool keep_lower) {
  constexpr double plane = 0.5;
  std::vector<runtime::Real3> result;
  result.reserve(4U);
  auto previous = vertices.back();
  bool previous_inside = keep_lower ? coordinate(previous, axis) <= plane
                                    : coordinate(previous, axis) >= plane;
  for (const auto current : vertices) {
    const bool current_inside = keep_lower
                                    ? coordinate(current, axis) <= plane
                                    : coordinate(current, axis) >= plane;
    if (current_inside != previous_inside)
      result.push_back(
          canonical_plane_intersection(previous, current, axis));
    if (current_inside)
      result.push_back(current);
    previous = current;
    previous_inside = current_inside;
  }
  return remove_adjacent_duplicate_vertices(std::move(result));
}

std::vector<test::StlFixtureTriangle>
split_triangle(const test::StlFixtureTriangle &triangle, std::size_t axis) {
  constexpr double plane = 0.5;
  double minimum = coordinate(triangle.vertices[0], axis);
  double maximum = minimum;
  for (const auto vertex : triangle.vertices) {
    minimum = std::min(minimum, coordinate(vertex, axis));
    maximum = std::max(maximum, coordinate(vertex, axis));
  }
  if (!(minimum < plane && maximum > plane))
    return {triangle};
  std::vector<test::StlFixtureTriangle> result;
  for (const bool keep_lower : {true, false}) {
    const auto polygon = clip_polygon(triangle.vertices, axis, keep_lower);
    if (polygon.size() < 3U)
      continue;
    for (std::size_t vertex = 1U; vertex + 1U < polygon.size(); ++vertex)
      result.push_back(test::StlFixtureTriangle{
          triangle.file_normal,
          {polygon[0], polygon[vertex], polygon[vertex + 1U]}});
  }
  return result;
}

bool strictly_crosses_plane(const test::StlFixtureTriangle &triangle,
                            std::size_t axis) noexcept {
  double minimum = coordinate(triangle.vertices[0], axis);
  double maximum = minimum;
  for (const auto vertex : triangle.vertices) {
    minimum = std::min(minimum, coordinate(vertex, axis));
    maximum = std::max(maximum, coordinate(vertex, axis));
  }
  return minimum < 0.5 && maximum > 0.5;
}

struct Bounds final {
  runtime::Real3 minimum{};
  runtime::Real3 maximum{};
};

Bounds bounds(const std::vector<test::StlFixtureTriangle> &triangles) {
  const double infinity = std::numeric_limits<double>::infinity();
  Bounds result{{infinity, infinity, infinity},
                {-infinity, -infinity, -infinity}};
  for (const auto &triangle : triangles)
    for (const auto vertex : triangle.vertices) {
      result.minimum.x = std::min(result.minimum.x, vertex.x);
      result.minimum.y = std::min(result.minimum.y, vertex.y);
      result.minimum.z = std::min(result.minimum.z, vertex.z);
      result.maximum.x = std::max(result.maximum.x, vertex.x);
      result.maximum.y = std::max(result.maximum.y, vertex.y);
      result.maximum.z = std::max(result.maximum.z, vertex.z);
    }
  return result;
}

struct SignedForceFixture final {
  std::vector<test::StlFixtureTriangle> intermediate;
  std::vector<test::StlFixtureTriangle> emitted;
  bool children_finite_positive{true};
  bool child_normal_orientation{true};
  bool no_strict_process_plane_crossing{true};
  bool per_source_scalar_area_preserved{true};
  bool per_source_oriented_area_preserved{true};
  bool total_scalar_area_preserved{};
  bool total_oriented_area_preserved{};
  bool aabb_unchanged{};
  bool bitwise_edge_incidence_closed{};

  bool valid() const noexcept {
    return children_finite_positive && child_normal_orientation &&
           no_strict_process_plane_crossing &&
           per_source_scalar_area_preserved &&
           per_source_oriented_area_preserved &&
           total_scalar_area_preserved && total_oriented_area_preserved &&
           aabb_unchanged && bitwise_edge_incidence_closed;
  }
};

SignedForceFixture signed_force_fixture() {
  SignedForceFixture result;
  result.intermediate = signed_force_intermediate_triangles();
  double intermediate_area = 0.0;
  runtime::Real3 intermediate_oriented_area{};
  for (const auto &source : result.intermediate) {
    intermediate_area += scalar_area(source);
    intermediate_oriented_area =
        add(intermediate_oriented_area, oriented_area_vector(source));
    std::vector<test::StlFixtureTriangle> children{source};
    for (const std::size_t axis : {0U, 1U}) {
      std::vector<test::StlFixtureTriangle> next;
      for (const auto &child : children) {
        auto split = split_triangle(child, axis);
        next.insert(next.end(), split.begin(), split.end());
      }
      children = std::move(next);
    }
    double child_area = 0.0;
    runtime::Real3 child_oriented_area{};
    for (const auto &child : children) {
      const auto area_vector = oriented_area_vector(child);
      const double area = norm(area_vector);
      child_area += area;
      child_oriented_area = add(child_oriented_area, area_vector);
      result.children_finite_positive =
          result.children_finite_positive && finite(child.file_normal) &&
          finite(child.vertices[0]) && finite(child.vertices[1]) &&
          finite(child.vertices[2]) && std::isfinite(area) && area > 0.0;
      result.child_normal_orientation =
          result.child_normal_orientation &&
          dot(area_vector, child.file_normal) > 0.0;
      result.no_strict_process_plane_crossing =
          result.no_strict_process_plane_crossing &&
          !strictly_crosses_plane(child, 0U) &&
          !strictly_crosses_plane(child, 1U);
    }
    result.per_source_scalar_area_preserved =
        result.per_source_scalar_area_preserved &&
        near(child_area, scalar_area(source));
    result.per_source_oriented_area_preserved =
        result.per_source_oriented_area_preserved &&
        near(child_oriented_area, oriented_area_vector(source));
    result.emitted.insert(result.emitted.end(), children.begin(),
                          children.end());
  }

  double emitted_area = 0.0;
  runtime::Real3 emitted_oriented_area{};
  using Edge = std::pair<VertexBits, VertexBits>;
  std::map<Edge, std::array<int, 2>> edges;
  for (const auto &triangle : result.emitted) {
    emitted_area += scalar_area(triangle);
    emitted_oriented_area =
        add(emitted_oriented_area, oriented_area_vector(triangle));
    for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
      const auto first = vertex_bits(triangle.vertices[vertex]);
      const auto second = vertex_bits(triangle.vertices[(vertex + 1U) % 3U]);
      if (first == second) {
        result.children_finite_positive = false;
        continue;
      }
      const bool forward = first < second;
      const Edge edge = forward ? Edge{first, second} : Edge{second, first};
      auto &incidence = edges[edge];
      ++incidence[0];
      incidence[1] += forward ? 1 : -1;
    }
  }
  result.total_scalar_area_preserved = near(emitted_area, intermediate_area);
  result.total_oriented_area_preserved =
      near(emitted_oriented_area, intermediate_oriented_area);
  const auto intermediate_bounds = bounds(result.intermediate);
  const auto emitted_bounds = bounds(result.emitted);
  result.aabb_unchanged = near(emitted_bounds.minimum,
                               intermediate_bounds.minimum) &&
                          near(emitted_bounds.maximum,
                               intermediate_bounds.maximum);
  result.bitwise_edge_incidence_closed = !edges.empty();
  for (const auto &[edge, incidence] : edges) {
    static_cast<void>(edge);
    result.bitwise_edge_incidence_closed =
        result.bitwise_edge_incidence_closed && incidence[0] == 2 &&
        incidence[1] == 0;
  }
  return result;
}

struct PlanFixtureObservations final {
  bool loaded_triangle_count_matches_emitted{};
  bool local_point_owner_matches_rank{};
  bool every_triangle_has_three_points{};
  bool every_triangle_point_index_mask_exact{};
  bool every_triangle_has_one_owner_rank{};
  bool triangle_zero_oracle_applicable{};
  bool triangle_zero_point_donors_present{};
  bool triangle_zero_pressure_authority_donors_present{};
  bool triangle_zero_authority_owners_exact{};
  bool triangle_zero_donor_x_bounds_exact{};
  bool triangle_zero_rank0_interval_exact{};
  bool triangle_zero_rank1_interval_exact{};
  bool triangle_zero_rank0_infeasible{};
  bool triangle_zero_rank1_feasible{};
  bool triangle_zero_execution_owner_one{};
  std::array<int, 3> triangle_zero_authority_owners{{-1, -1, -1}};
  int triangle_zero_donor_x_min{kExtent.x};
  int triangle_zero_donor_x_max{-1};
  int triangle_zero_rank0_interval_begin{};
  int triangle_zero_rank0_interval_end{};
  int triangle_zero_rank1_interval_begin{};
  int triangle_zero_rank1_interval_end{};
  int triangle_zero_execution_owner{-1};
  std::uint64_t triangle_zero_point_donor_count{};
  std::uint64_t triangle_zero_pressure_authority_donor_count{};

  bool valid() const noexcept {
    const bool common = loaded_triangle_count_matches_emitted &&
                        local_point_owner_matches_rank &&
                        every_triangle_has_three_points &&
                        every_triangle_point_index_mask_exact &&
                        every_triangle_has_one_owner_rank;
    const bool triangle_zero =
        !triangle_zero_oracle_applicable ||
        (triangle_zero_point_donors_present &&
         triangle_zero_pressure_authority_donors_present &&
         triangle_zero_authority_owners_exact &&
         triangle_zero_donor_x_bounds_exact &&
         triangle_zero_rank0_interval_exact &&
         triangle_zero_rank1_interval_exact && triangle_zero_rank0_infeasible &&
         triangle_zero_rank1_feasible && triangle_zero_execution_owner_one);
    return common && triangle_zero;
  }
};

PlanFixtureObservations
plan_fixture_observations(const runtime::MpiContext &mpi,
                          const runtime::StructuredDecomposition &decomposition,
                          const immersed::ImmersedSurface &surface,
                          const immersed::WallQuadraturePlan &wall_plan,
                          std::size_t emitted_triangle_count) {
  if (emitted_triangle_count >
      static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw runtime::Error("Task 11 RED-S1 emitted triangle count overflows");
  std::vector<int> point_count(emitted_triangle_count);
  std::vector<int> point_index_mask(emitted_triangle_count);
  std::vector<int> owner_min(emitted_triangle_count, mpi.size());
  std::vector<int> owner_max(emitted_triangle_count, -1);
  const std::size_t global_cell_count = static_cast<std::size_t>(kExtent.x) *
                                        static_cast<std::size_t>(kExtent.y) *
                                        static_cast<std::size_t>(kExtent.z);
  std::vector<int> triangle_zero_point_donor_mask(global_cell_count);
  std::vector<int> triangle_zero_pressure_donor_mask(global_cell_count);
  std::array<int, 3> triangle_zero_authority_owner_min{mpi.size(), mpi.size(),
                                                       mpi.size()};
  std::array<int, 3> triangle_zero_authority_owner_max{{-1, -1, -1}};
  std::uint64_t triangle_zero_point_donor_count = 0U;
  std::uint64_t triangle_zero_pressure_donor_count = 0U;
  bool local_identity_valid = true;
  bool local_owner_matches = true;
  for (const auto &point : wall_plan.local_points()) {
    local_owner_matches =
        local_owner_matches && point.owner_rank == mpi.rank();
    if (point.triangle >= emitted_triangle_count || point.point_index >= 3U ||
        point.owner_rank < 0 || point.owner_rank >= mpi.size()) {
      local_identity_valid = false;
      continue;
    }
    const auto triangle = static_cast<std::size_t>(point.triangle);
    ++point_count[triangle];
    point_index_mask[triangle] |= 1 << point.point_index;
    owner_min[triangle] = std::min(owner_min[triangle], point.owner_rank);
    owner_max[triangle] = std::max(owner_max[triangle], point.owner_rank);
    if (point.triangle != 0U)
      continue;
    const auto point_index = static_cast<std::size_t>(point.point_index);
    const int authority_owner =
        immersed::detail::boundary_authority_owner_rank(point.reconstruction);
    triangle_zero_authority_owner_min[point_index] = std::min(
        triangle_zero_authority_owner_min[point_index], authority_owner);
    triangle_zero_authority_owner_max[point_index] = std::max(
        triangle_zero_authority_owner_max[point_index], authority_owner);
    const auto mark_donors = [&](const immersed::QuadraticReconstruction
                                     &reconstruction,
                                 std::vector<int> &mask, std::uint64_t &count) {
      const auto &donors =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              reconstruction);
      count += donors.size();
      for (const auto donor : donors) {
        if (donor >= global_cell_count) {
          local_identity_valid = false;
          continue;
        }
        mask[static_cast<std::size_t>(donor)] = 1;
      }
    };
    mark_donors(point.reconstruction, triangle_zero_point_donor_mask,
                triangle_zero_point_donor_count);
    mark_donors(immersed::detail::boundary_authority_reconstruction(
                    point.reconstruction),
                triangle_zero_pressure_donor_mask,
                triangle_zero_pressure_donor_count);
  }
  allreduce_sum(point_count, MPI_INT, mpi,
                "MPI_Allreduce fixture triangle point count");
  if (!point_index_mask.empty())
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, point_index_mask.data(),
                            static_cast<int>(point_index_mask.size()), MPI_INT,
                            MPI_BOR, mpi.comm()),
              "MPI_Allreduce fixture triangle point masks");
  if (!owner_min.empty()) {
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, owner_min.data(),
                            static_cast<int>(owner_min.size()), MPI_INT,
                            MPI_MIN, mpi.comm()),
              "MPI_Allreduce fixture triangle owner minimum");
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, owner_max.data(),
                            static_cast<int>(owner_max.size()), MPI_INT,
                            MPI_MAX, mpi.comm()),
              "MPI_Allreduce fixture triangle owner maximum");
  }
  check_mpi(
      MPI_Allreduce(MPI_IN_PLACE, triangle_zero_authority_owner_min.data(),
                    static_cast<int>(triangle_zero_authority_owner_min.size()),
                    MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce triangle zero authority owner minimum");
  check_mpi(
      MPI_Allreduce(MPI_IN_PLACE, triangle_zero_authority_owner_max.data(),
                    static_cast<int>(triangle_zero_authority_owner_max.size()),
                    MPI_INT, MPI_MAX, mpi.comm()),
      "MPI_Allreduce triangle zero authority owner maximum");
  for (auto *mask :
       {&triangle_zero_point_donor_mask, &triangle_zero_pressure_donor_mask})
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, mask->data(),
                            static_cast<int>(mask->size()), MPI_INT, MPI_BOR,
                            mpi.comm()),
              "MPI_Allreduce triangle zero donor mask");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &triangle_zero_point_donor_count, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce triangle zero point donor count");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &triangle_zero_pressure_donor_count, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce triangle zero pressure donor count");
  const auto local_box = decomposition.owned_box();
  const std::array<int, 6> local_box_wire{local_box.begin.x, local_box.begin.y,
                                          local_box.begin.z, local_box.end.x,
                                          local_box.end.y,   local_box.end.z};
  std::vector<int> owner_box_wire(static_cast<std::size_t>(mpi.size()) * 6U);
  check_mpi(MPI_Allgather(
                local_box_wire.data(), static_cast<int>(local_box_wire.size()),
                MPI_INT, owner_box_wire.data(),
                static_cast<int>(local_box_wire.size()), MPI_INT, mpi.comm()),
            "MPI_Allgather triangle zero owner boxes");

  PlanFixtureObservations result;
  result.loaded_triangle_count_matches_emitted = all_true(
      surface.triangle_count() == emitted_triangle_count, mpi,
      "MPI_Allreduce loaded fixture triangle count");
  result.local_point_owner_matches_rank =
      all_true(local_identity_valid && local_owner_matches, mpi,
               "MPI_Allreduce fixture local point ownership");
  result.every_triangle_has_three_points = emitted_triangle_count > 0U;
  result.every_triangle_point_index_mask_exact = emitted_triangle_count > 0U;
  result.every_triangle_has_one_owner_rank = emitted_triangle_count > 0U;
  for (std::size_t triangle = 0U; triangle < emitted_triangle_count;
       ++triangle) {
    result.every_triangle_has_three_points =
        result.every_triangle_has_three_points && point_count[triangle] == 3;
    result.every_triangle_point_index_mask_exact =
        result.every_triangle_point_index_mask_exact &&
        point_index_mask[triangle] == 0b111;
    result.every_triangle_has_one_owner_rank =
        result.every_triangle_has_one_owner_rank && owner_min[triangle] >= 0 &&
        owner_min[triangle] < mpi.size() &&
        owner_min[triangle] == owner_max[triangle];
  }
  result.triangle_zero_oracle_applicable = mpi.size() == 2;
  if (result.triangle_zero_oracle_applicable) {
    result.triangle_zero_authority_owners = triangle_zero_authority_owner_min;
    result.triangle_zero_point_donor_count = triangle_zero_point_donor_count;
    result.triangle_zero_pressure_authority_donor_count =
        triangle_zero_pressure_donor_count;
    result.triangle_zero_point_donors_present =
        triangle_zero_point_donor_count > 0U;
    result.triangle_zero_pressure_authority_donors_present =
        triangle_zero_pressure_donor_count > 0U;
    result.triangle_zero_authority_owners_exact =
        triangle_zero_authority_owner_min == std::array<int, 3>{0, 0, 1} &&
        triangle_zero_authority_owner_min == triangle_zero_authority_owner_max;
    const runtime::Box3 rank0{
        {owner_box_wire[0], owner_box_wire[1], owner_box_wire[2]},
        {owner_box_wire[3], owner_box_wire[4], owner_box_wire[5]}};
    const runtime::Box3 rank1{
        {owner_box_wire[6], owner_box_wire[7], owner_box_wire[8]},
        {owner_box_wire[9], owner_box_wire[10], owner_box_wire[11]}};
    bool donor_present = false;
    bool rank0_feasible = true;
    bool rank1_feasible = true;
    for (std::size_t donor = 0U; donor < global_cell_count; ++donor) {
      if (triangle_zero_point_donor_mask[donor] == 0 &&
          triangle_zero_pressure_donor_mask[donor] == 0)
        continue;
      donor_present = true;
      const auto logical = logical_cell(donor);
      result.triangle_zero_donor_x_min =
          std::min(result.triangle_zero_donor_x_min, logical.x);
      result.triangle_zero_donor_x_max =
          std::max(result.triangle_zero_donor_x_max, logical.x);
      rank0_feasible =
          rank0_feasible && inside_expanded(rank0, logical, kWallReach);
      rank1_feasible =
          rank1_feasible && inside_expanded(rank1, logical, kWallReach);
    }
    result.triangle_zero_donor_x_bounds_exact =
        donor_present && result.triangle_zero_donor_x_min == 4 &&
        result.triangle_zero_donor_x_max == 10;
    result.triangle_zero_rank0_interval_begin = rank0.begin.x - kWallReach;
    result.triangle_zero_rank0_interval_end = rank0.end.x + kWallReach;
    result.triangle_zero_rank1_interval_begin = rank1.begin.x - kWallReach;
    result.triangle_zero_rank1_interval_end = rank1.end.x + kWallReach;
    result.triangle_zero_rank0_interval_exact =
        result.triangle_zero_rank0_interval_begin == -4 &&
        result.triangle_zero_rank0_interval_end == 10;
    result.triangle_zero_rank1_interval_exact =
        result.triangle_zero_rank1_interval_begin == 2 &&
        result.triangle_zero_rank1_interval_end == 16;
    result.triangle_zero_rank0_infeasible = !rank0_feasible;
    result.triangle_zero_rank1_feasible = rank1_feasible;
    if (!owner_min.empty() && owner_min[0] == owner_max[0])
      result.triangle_zero_execution_owner = owner_min[0];
    result.triangle_zero_execution_owner_one =
        result.triangle_zero_execution_owner == 1;
  }
  return result;
}

void print_pre_product_fixture_observations(
    const SignedForceFixture &fixture, const runtime::MpiContext &mpi) {
  if (mpi.rank() != 0)
    return;
  std::cout << "fixture.intermediate_triangle_count="
            << fixture.intermediate.size() << '\n'
            << "fixture.emitted_triangle_count=" << fixture.emitted.size()
            << '\n'
            << "fixture.children_finite_positive=" << std::boolalpha
            << fixture.children_finite_positive << '\n'
            << "fixture.child_normal_orientation="
            << fixture.child_normal_orientation << '\n'
            << "fixture.no_strict_process_plane_crossing="
            << fixture.no_strict_process_plane_crossing << '\n'
            << "fixture.per_source_scalar_area_preserved="
            << fixture.per_source_scalar_area_preserved << '\n'
            << "fixture.per_source_oriented_area_preserved="
            << fixture.per_source_oriented_area_preserved << '\n'
            << "fixture.total_scalar_area_preserved="
            << fixture.total_scalar_area_preserved << '\n'
            << "fixture.total_oriented_area_preserved="
            << fixture.total_oriented_area_preserved << '\n'
            << "fixture.aabb_unchanged=" << fixture.aabb_unchanged << '\n'
            << "fixture.bitwise_edge_incidence_closed="
            << fixture.bitwise_edge_incidence_closed << '\n';
}

void print_plan_fixture_observations(
    const PlanFixtureObservations &observations,
    const runtime::MpiContext &mpi) {
  if (mpi.rank() != 0)
    return;
  std::cout << "fixture.loaded_triangle_count_matches_emitted="
            << std::boolalpha
            << observations.loaded_triangle_count_matches_emitted << '\n'
            << "fixture.local_point_owner_matches_rank="
            << observations.local_point_owner_matches_rank << '\n'
            << "fixture.every_triangle_has_three_points="
            << observations.every_triangle_has_three_points << '\n'
            << "fixture.every_triangle_point_index_mask_exact="
            << observations.every_triangle_point_index_mask_exact << '\n'
            << "fixture.every_triangle_has_one_owner_rank="
            << observations.every_triangle_has_one_owner_rank << '\n'
            << "fixture.triangle0_owner_oracle_applicable="
            << observations.triangle_zero_oracle_applicable << '\n';
  if (observations.triangle_zero_oracle_applicable)
    std::cout << "fixture.triangle0_point_authority_owners=["
              << observations.triangle_zero_authority_owners[0] << ','
              << observations.triangle_zero_authority_owners[1] << ','
              << observations.triangle_zero_authority_owners[2] << "]\n"
              << "fixture.triangle0_authority_owners_exact="
              << observations.triangle_zero_authority_owners_exact << '\n'
              << "fixture.triangle0_point_donor_count="
              << observations.triangle_zero_point_donor_count << '\n'
              << "fixture.triangle0_pressure_authority_donor_count="
              << observations.triangle_zero_pressure_authority_donor_count
              << '\n'
              << "fixture.triangle0_complete_donor_x_bounds=["
              << observations.triangle_zero_donor_x_min << ','
              << observations.triangle_zero_donor_x_max << "]\n"
              << "fixture.triangle0_donor_x_bounds_exact="
              << observations.triangle_zero_donor_x_bounds_exact << '\n'
              << "fixture.triangle0_rank0_reach_interval=["
              << observations.triangle_zero_rank0_interval_begin << ','
              << observations.triangle_zero_rank0_interval_end << ")\n"
              << "fixture.triangle0_rank1_reach_interval=["
              << observations.triangle_zero_rank1_interval_begin << ','
              << observations.triangle_zero_rank1_interval_end << ")\n"
              << "fixture.triangle0_rank0_infeasible="
              << observations.triangle_zero_rank0_infeasible << '\n'
              << "fixture.triangle0_rank1_feasible="
              << observations.triangle_zero_rank1_feasible << '\n'
              << "fixture.triangle0_execution_owner="
              << observations.triangle_zero_execution_owner << '\n'
              << "fixture.triangle0_execution_owner_one="
              << observations.triangle_zero_execution_owner_one << '\n';
  std::cout.flush();
}

class FixtureFile final {
public:
  FixtureFile(const runtime::MpiContext &mpi,
              const std::vector<test::StlFixtureTriangle> &triangles)
      : mpi_(&mpi) {
    std::string path_text;
    if (mpi.rank() == 0) {
      path_ = std::filesystem::temp_directory_path() /
              ("hundun-task11-red-s1-" +
               std::to_string(std::chrono::steady_clock::now()
                                  .time_since_epoch()
                                  .count()) +
               ".stl");
      test::write_text(path_,
                       test::ascii_stl(triangles, "task11-red-s1-box"));
      path_text = path_.string();
    }
    std::uint64_t length = path_text.size();
    check_mpi(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, mpi.comm()),
              "MPI_Bcast fixture path length");
    if (length > static_cast<std::uint64_t>(
                     std::numeric_limits<int>::max()))
      throw runtime::Error("Task 11 RED-S1 fixture path is too long");
    path_text.resize(static_cast<std::size_t>(length));
    check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                        MPI_BYTE, 0, mpi.comm()),
              "MPI_Bcast fixture path");
    path_ = path_text;
    check_mpi(MPI_Barrier(mpi.comm()), "MPI_Barrier fixture creation");
  }

  ~FixtureFile() {
    MPI_Barrier(mpi_->comm());
    if (mpi_->rank() == 0) {
      std::error_code error;
      std::filesystem::remove(path_, error);
    }
  }

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  const runtime::MpiContext *mpi_{};
  std::filesystem::path path_;
};

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("Task 11 RED-S1 requires 1, 2, or 4 MPI ranks");
}

config::FlowCaseConfig case_config(int ranks, runtime::Int3 grid) {
  config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task11-red-s1-signed-force";
  config.simulation_type = config::SimulationType::variable_density_flow;
  config.density_model = config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = grid;
  config.mesh.cells = kExtent;
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = config::MeshMapping::uniform_box;
  config.time.mode = config::TimeMode::fixed;
  config.time.steps = 1;
  config.time.initial_dt_s = 0.01;
  config.time.min_dt_s = 0.01;
  config.time.max_dt_s = 0.01;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = kViscosity;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = config::BoundaryType::no_slip_wall;
  }
  return config;
}

struct IndependentPoint final {
  runtime::Real3 position{};
  runtime::Real3 normal{};
  double weight{};
};

IndependentPoint independent_point(
    const std::vector<test::StlFixtureTriangle> &triangles,
    immersed::TriangleId triangle_id, std::uint32_t point_index) {
  if (triangle_id >= triangles.size() || point_index >= 3U)
    throw runtime::Error("Task 11 RED-S1 quadrature identity is invalid");
  constexpr std::array<std::array<double, 3>, 3> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  const auto &triangle = triangles[static_cast<std::size_t>(triangle_id)];
  IndependentPoint result;
  for (std::size_t vertex = 0U; vertex < 3U; ++vertex)
    result.position =
        add(result.position,
            multiply(barycentric[point_index][vertex],
                     triangle.vertices[vertex]));
  const auto area_vector =
      cross(subtract(triangle.vertices[1], triangle.vertices[0]),
            subtract(triangle.vertices[2], triangle.vertices[0]));
  const double twice_area = norm(area_vector);
  result.normal = multiply(1.0 / twice_area, area_vector);
  result.weight = 0.5 * twice_area / 3.0;
  return result;
}

struct ExpectedForces final {
  immersed::ForceComponents operator_force;
  immersed::ForceComponents surface_force;
};

ExpectedForces expected_forces(runtime::Real3 normal,
                               runtime::Real3 background_area,
                               runtime::Real3 quadrature_normal,
                               double surface_area) {
  const auto stress = analytic_stress(normal);
  const auto operator_pressure = multiply(kPressure, background_area);
  const auto operator_viscous =
      multiply(-1.0, multiply(stress.tau, background_area));
  const auto surface_pressure =
      multiply(-kPressure * surface_area, quadrature_normal);
  const auto surface_viscous =
      multiply(surface_area, multiply(stress.tau, quadrature_normal));
  return {force(operator_pressure, operator_viscous),
          force(surface_pressure, surface_viscous)};
}

struct Preflight final {
  bool valid{};
  immersed::ImmersedLinkId link{};
  mesh::GlobalCellId active_cell{};
  mesh::GlobalCellId fluid_cell{};
  mesh::GlobalCellId solid_cell{};
  immersed::TriangleId triangle{};
  int domain_owner{-1};
  int row_owner{-1};
  int point_count{};
  runtime::Real3 wall_intercept{};
  runtime::Real3 solid_to_fluid_normal{};
  runtime::Real3 background_area{};
  runtime::Real3 quadrature_normal{};
  double surface_area{};
  double signed_background_measure{};
  bool selected_one_link_row{};
  bool oblique_normal{};
  bool surface_points_present{};
  bool planar_link_patch{};
  bool domain_normal_equals_triangle_normal{};
  bool wall_intercept_on_triangle_plane{};
  bool quadrature_geometry_equals_fixture{};
  bool operator_measure_nonzero{};
  bool surface_measure_nonzero{};
  bool operator_surface_measures_distinct{};
  ExpectedForces expected{};
};

Preflight select_preflight(
    const runtime::MpiContext &mpi, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain,
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const immersed::WallQuadraturePlan &wall_plan,
    const std::vector<test::StlFixtureTriangle> &triangles) {
  const auto rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
  std::uint64_t local_max = 0U;
  int local_any = 0;
  const auto observe_id = [&](std::uint64_t id) {
    local_any = 1;
    local_max = std::max(local_max, id);
  };
  for (const auto &link : domain.links())
    observe_id(link.id);
  for (const auto &row : rows)
    for (const auto &link : row.links)
      observe_id(link.id);
  for (const auto &point : wall_plan.local_points())
    observe_id(immersed::detail::boundary_authority_link(
        point.reconstruction));
  int global_any = 0;
  std::uint64_t global_max = 0U;
  check_mpi(MPI_Allreduce(&local_any, &global_any, 1, MPI_INT, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce preflight link presence");
  check_mpi(MPI_Allreduce(&local_max, &global_max, 1, MPI_UINT64_T, MPI_MAX,
                          mpi.comm()),
            "MPI_Allreduce preflight maximum link");
  if (global_any == 0 ||
      global_max >= static_cast<std::uint64_t>(
                        std::numeric_limits<int>::max()))
    return {};
  const std::size_t count = static_cast<std::size_t>(global_max + 1U);
  std::vector<int> domain_count(count), row_membership_count(count),
      single_row_count(count), face_count(count), point_count(count),
      geometry_bad(count), domain_owner_plus_one(count),
      row_owner_plus_one(count);
  std::vector<std::uint64_t> fluid_cell(count), solid_cell(count),
      triangle_id(count), active_cell(count);
  std::vector<double> wall_x(count), wall_y(count), wall_z(count), normal_x(count),
      normal_y(count), normal_z(count), area_x(count), area_y(count),
      area_z(count), point_weight(count), point_normal_x(count),
      point_normal_y(count), point_normal_z(count);

  for (const auto &link : domain.links()) {
    const auto index = static_cast<std::size_t>(link.id);
    ++domain_count[index];
    domain_owner_plus_one[index] += mpi.rank() + 1;
    fluid_cell[index] += link.fluid_cell;
    solid_cell[index] += link.solid_cell;
    triangle_id[index] += link.triangle;
    wall_x[index] += link.wall_intercept_m.x;
    wall_y[index] += link.wall_intercept_m.y;
    wall_z[index] += link.wall_intercept_m.z;
    normal_x[index] += link.solid_to_fluid_normal.x;
    normal_y[index] += link.solid_to_fluid_normal.y;
    normal_z[index] += link.solid_to_fluid_normal.z;
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value())
        continue;
      const auto owner_global = topology.global_cell_id(topology.owner(face));
      const auto neighbour_global = topology.global_cell_id(*neighbour);
      const bool owner_is_fluid = owner_global == link.fluid_cell &&
                                  neighbour_global == link.solid_cell;
      const bool neighbour_is_fluid = neighbour_global == link.fluid_cell &&
                                      owner_global == link.solid_cell;
      if (!owner_is_fluid && !neighbour_is_fluid)
        continue;
      ++face_count[index];
      const auto area = geometry.face_area_vector_m2(
          face, owner_is_fluid ? mesh::FaceSide::owner
                               : mesh::FaceSide::neighbour);
      area_x[index] += area.x;
      area_y[index] += area.y;
      area_z[index] += area.z;
    }
  }
  for (const auto &row : rows) {
    for (const auto &link : row.links)
      ++row_membership_count[static_cast<std::size_t>(link.id)];
    if (row.links.size() == 1U) {
      const auto index = static_cast<std::size_t>(row.links.front().id);
      ++single_row_count[index];
      row_owner_plus_one[index] += mpi.rank() + 1;
      active_cell[index] += row.active_cell;
    }
  }
  const std::size_t key_count_size = triangles.size() * 3U;
  std::vector<int> membership_key_count(key_count_size);
  std::vector<std::uint64_t> membership_key_link_plus_one(key_count_size);
  for (const auto &point : wall_plan.local_points()) {
    const auto link =
        immersed::detail::boundary_authority_link(point.reconstruction);
    const auto index = static_cast<std::size_t>(link);
    ++point_count[index];
    const auto expected =
        independent_point(triangles, point.triangle, point.point_index);
    point_weight[index] += expected.weight;
    point_normal_x[index] += expected.weight * expected.normal.x;
    point_normal_y[index] += expected.weight * expected.normal.y;
    point_normal_z[index] += expected.weight * expected.normal.z;
    if (!near(point.position_m, expected.position) ||
        !near(point.solid_to_fluid_normal, expected.normal) ||
        !near(point.weight_m2, expected.weight))
      ++geometry_bad[index];
    const auto key = static_cast<std::size_t>(point.triangle) * 3U +
                     static_cast<std::size_t>(point.point_index);
    if (key >= membership_key_count.size())
      throw runtime::Error("Task 11 RED-S1 membership key is invalid");
    ++membership_key_count[key];
    membership_key_link_plus_one[key] += link + 1U;
  }

  for (auto *values : {&domain_count, &row_membership_count,
                       &single_row_count, &face_count, &point_count,
                       &geometry_bad, &domain_owner_plus_one,
                       &row_owner_plus_one, &membership_key_count})
    allreduce_sum(*values, MPI_INT, mpi,
                  "MPI_Allreduce preflight integer table");
  for (auto *values : {&fluid_cell, &solid_cell, &triangle_id, &active_cell,
                       &membership_key_link_plus_one})
    allreduce_sum(*values, MPI_UINT64_T, mpi,
                  "MPI_Allreduce preflight identity table");
  for (auto *values : {&wall_x, &wall_y, &wall_z, &normal_x, &normal_y,
                       &normal_z, &area_x, &area_y, &area_z, &point_weight,
                       &point_normal_x, &point_normal_y, &point_normal_z})
    allreduce_sum(*values, MPI_DOUBLE, mpi,
                  "MPI_Allreduce preflight geometry table");

  for (std::size_t index = 0U; index < count; ++index) {
    if (domain_count[index] != 1 || row_membership_count[index] != 1 ||
        single_row_count[index] != 1 || face_count[index] != 1 ||
        point_count[index] < 3 || geometry_bad[index] != 0)
      continue;
    if (triangle_id[index] >= triangles.size())
      continue;
    const runtime::Real3 normal{normal_x[index], normal_y[index],
                                normal_z[index]};
    const runtime::Real3 background_area{area_x[index], area_y[index],
                                         area_z[index]};
    const runtime::Real3 weighted_normal{point_normal_x[index],
                                         point_normal_y[index],
                                         point_normal_z[index]};
    const double surface_area = point_weight[index];
    if (!(surface_area > 0.0) || !(norm(weighted_normal) > 0.0))
      continue;
    const auto quadrature_normal = normalized(weighted_normal);
    const auto &triangle = triangles[static_cast<std::size_t>(
        triangle_id[index])];
    const auto fixture_normal = triangle.file_normal;
    const runtime::Real3 wall_intercept{wall_x[index], wall_y[index],
                                        wall_z[index]};
    const bool oblique =
        std::min({std::abs(normal.x), std::abs(normal.y),
                  std::abs(normal.z)}) >= 0.2;
    const double planar_tolerance =
        4096.0 * std::numeric_limits<double>::epsilon();
    const bool planar =
        std::abs(norm(weighted_normal) / surface_area - 1.0) <=
            planar_tolerance &&
        dot(quadrature_normal, fixture_normal) >= 1.0 - planar_tolerance;
    const bool domain_normal = near(normal, fixture_normal);
    const bool wall_plane =
        std::abs(dot(subtract(wall_intercept, triangle.vertices[0]),
                     fixture_normal)) <=
        kArithmeticFactor * std::numeric_limits<double>::epsilon();
    const double signed_measure = -dot(background_area, normal);
    const auto expected = expected_forces(normal, background_area,
                                          quadrature_normal, surface_area);
    const bool operator_nonparallel =
        norm(cross(expected.operator_force.pressure_N,
                   expected.operator_force.viscous_N)) >
        kArithmeticFactor * std::numeric_limits<double>::epsilon();
    const bool surface_nonparallel =
        norm(cross(expected.surface_force.pressure_N,
                   expected.surface_force.viscous_N)) >
        kArithmeticFactor * std::numeric_limits<double>::epsilon();
    const bool operator_nonzero =
        max_abs(expected.operator_force.total_N) >
        kArithmeticFactor * std::numeric_limits<double>::epsilon();
    const bool surface_nonzero =
        max_abs(expected.surface_force.total_N) >
        kArithmeticFactor * std::numeric_limits<double>::epsilon();
    const auto operator_measure = background_area;
    const auto surface_measure = multiply(-surface_area, quadrature_normal);
    const bool distinct =
        max_abs(subtract(operator_measure, surface_measure)) >
        arithmetic_bound(operator_measure, surface_measure);
    int unique_memberships = 0;
    for (std::size_t key = 0U; key < membership_key_count.size(); ++key)
      if (membership_key_count[key] == 1 &&
          membership_key_link_plus_one[key] == index + 1U)
        ++unique_memberships;
    const bool membership_unique = unique_memberships == point_count[index];
    if (!oblique || !planar || !domain_normal || !wall_plane ||
        !(signed_measure > 0.0) || !operator_nonparallel ||
        !surface_nonparallel || !operator_nonzero || !surface_nonzero ||
        !distinct || !membership_unique)
      continue;
    Preflight result;
    result.valid = true;
    result.link = static_cast<immersed::ImmersedLinkId>(index);
    result.active_cell = active_cell[index];
    result.fluid_cell = fluid_cell[index];
    result.solid_cell = solid_cell[index];
    result.triangle = triangle_id[index];
    result.domain_owner = domain_owner_plus_one[index] - 1;
    result.row_owner = row_owner_plus_one[index] - 1;
    result.point_count = point_count[index];
    result.wall_intercept = wall_intercept;
    result.solid_to_fluid_normal = normal;
    result.background_area = background_area;
    result.quadrature_normal = quadrature_normal;
    result.surface_area = surface_area;
    result.signed_background_measure = signed_measure;
    result.selected_one_link_row = true;
    result.oblique_normal = oblique;
    result.surface_points_present = point_count[index] >= 3;
    result.planar_link_patch = planar && membership_unique;
    result.domain_normal_equals_triangle_normal = domain_normal;
    result.wall_intercept_on_triangle_plane = wall_plane;
    result.quadrature_geometry_equals_fixture = geometry_bad[index] == 0;
    result.operator_measure_nonzero = operator_nonzero;
    result.surface_measure_nonzero = surface_nonzero;
    result.operator_surface_measures_distinct = distinct;
    result.expected = expected;
    return result;
  }
  return {};
}

void print_vector(std::ostream &stream, runtime::Real3 value) {
  stream << '(' << value.x << ',' << value.y << ',' << value.z << ')';
}

void print_preflight(const Preflight &preflight,
                     const runtime::MpiContext &mpi) {
  for (int rank = 0; rank < mpi.size(); ++rank) {
    check_mpi(MPI_Barrier(mpi.comm()), "MPI_Barrier preflight print");
    if (mpi.rank() != rank)
      continue;
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "fixture_preflight.rank=" << rank
              << " fixture_preflight.size=" << mpi.size()
              << " selected_link=" << preflight.link
              << " selected_row=" << preflight.active_cell
              << " selected_triangle=" << preflight.triangle
              << " domain_owner=" << preflight.domain_owner
              << " row_owner=" << preflight.row_owner
              << " point_count=" << preflight.point_count << " n_s=";
    print_vector(std::cout, preflight.solid_to_fluid_normal);
    std::cout << " x_wall=";
    print_vector(std::cout, preflight.wall_intercept);
    std::cout << " A_background=";
    print_vector(std::cout, preflight.background_area);
    std::cout << " n_q=";
    print_vector(std::cout, preflight.quadrature_normal);
    std::cout << " A_surface=" << preflight.surface_area
              << " A_signed=" << preflight.signed_background_measure
              << " preflight.ok=" << std::boolalpha << preflight.valid
              << '\n';
  }
  check_mpi(MPI_Barrier(mpi.comm()), "MPI_Barrier preflight print complete");
}

runtime::FieldDescriptor cell_descriptor(std::string name,
                                         std::uint32_t components,
                                         int ghost) {
  return {std::move(name),
          "1",
          "task11-red-s1",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_descriptor(std::string name,
                                         std::uint32_t components) {
  return {std::move(name),
          "1",
          "task11-red-s1",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

struct FieldIds final {
  runtime::FieldId pressure{};
  runtime::FieldId velocity{};
  runtime::FieldId gradient{};
  runtime::FieldId cell_viscosity{};
  runtime::FieldId face_velocity{};
  runtime::FieldId face_viscosity{};
  runtime::FieldId mass_flux{};
  runtime::FieldId residual{};
};

FieldIds declare_fields(runtime::FieldRegistry &registry, int ghost) {
  return {registry.declare_field(cell_descriptor("pressure", 1U, ghost)),
          registry.declare_field(cell_descriptor("velocity", 3U, ghost)),
          registry.declare_field(cell_descriptor("gradient", 9U, ghost)),
          registry.declare_field(
              cell_descriptor("cell_viscosity", 1U, ghost)),
          registry.declare_field(face_descriptor("face_velocity", 3U)),
          registry.declare_field(face_descriptor("face_viscosity", 1U)),
          finite_volume::declare_face_mass_flux(registry),
          registry.declare_field(cell_descriptor("residual", 3U, 0))};
}

struct CaseInput final {
  double pressure_offset{};
  runtime::Real3 pressure_gradient{};
  std::array<double, 9> velocity_gradient{};
  double viscosity{};
};

bool inside(runtime::Int3 global) noexcept {
  return global.x >= 0 && global.y >= 0 && global.z >= 0 &&
         global.x < kExtent.x && global.y < kExtent.y &&
         global.z < kExtent.z;
}

runtime::Real3 cell_center(runtime::Int3 global) noexcept {
  constexpr double h = 1.0 / 12.0;
  return {(static_cast<double>(global.x) + 0.5) * h,
          (static_cast<double>(global.y) + 0.5) * h,
          (static_cast<double>(global.z) + 0.5) * h};
}

void fill_fields(runtime::FieldStorage &storage,
                 const runtime::FieldAccessPlan &access,
                 const FieldIds &ids,
                 const runtime::StructuredDecomposition &decomposition,
                 const mesh::MeshTopology &topology, int ghost,
                 runtime::Real3 wall_intercept, const CaseInput &input) {
  auto pressure = storage.view<double>(ids.pressure);
  auto velocity = storage.view<double>(ids.velocity);
  auto gradient = storage.view<double>(ids.gradient);
  auto cell_viscosity = storage.view<double>(ids.cell_viscosity);
  const auto local = decomposition.local_extent();
  const auto box = decomposition.owned_box();
  for (int k = -ghost; k < local.z + ghost; ++k)
    for (int j = -ghost; j < local.y + ghost; ++j)
      for (int i = -ghost; i < local.x + ghost; ++i) {
        const runtime::Int3 global{box.begin.x + i, box.begin.y + j,
                                   box.begin.z + k};
        const auto point = cell_center(global);
        const auto displacement = subtract(point, wall_intercept);
        const bool valid = inside(global);
        pressure(i, j, k, 0) =
            valid ? input.pressure_offset +
                        dot(input.pressure_gradient, displacement)
                  : 0.0;
        const auto value = multiply(input.velocity_gradient, displacement);
        velocity(i, j, k, 0) = valid ? value.x : 0.0;
        velocity(i, j, k, 1) = valid ? value.y : 0.0;
        velocity(i, j, k, 2) = valid ? value.z : 0.0;
        for (std::size_t component = 0U; component < 9U; ++component)
          gradient(i, j, k, static_cast<int>(component)) =
              valid ? input.velocity_gradient[component] : 0.0;
        cell_viscosity(i, j, k, 0) = valid ? input.viscosity : 0.0;
      }
  auto face_velocity = storage.acquire_face_write<double>(
      access, kPhase, kActor, ids.face_velocity);
  auto face_viscosity = storage.acquire_face_write<double>(
      access, kPhase, kActor, ids.face_viscosity);
  auto mass_flux = storage.acquire_face_write<double>(
      access, kPhase, kActor, ids.mass_flux);
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    for (int component = 0; component < 3; ++component)
      face_velocity(face, component) = 0.0;
    face_viscosity(face, 0) = input.viscosity;
    mass_flux(face, 0) = 0.0;
  }
  auto residual = storage.view<double>(ids.residual);
  for (int k = 0; k < local.z; ++k)
    for (int j = 0; j < local.y; ++j)
      for (int i = 0; i < local.x; ++i)
        for (int component = 0; component < 3; ++component)
          residual(i, j, k, component) = 0.0;
}

void reduce_force(immersed::ForceComponents &value,
                  const runtime::MpiContext &mpi, const char *operation) {
  std::array<double, 9> packed{
      value.pressure_N.x, value.pressure_N.y, value.pressure_N.z,
      value.viscous_N.x,  value.viscous_N.y,  value.viscous_N.z,
      value.total_N.x,    value.total_N.y,    value.total_N.z};
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, packed.data(),
                          static_cast<int>(packed.size()), MPI_DOUBLE, MPI_SUM,
                          mpi.comm()),
            operation);
  value.pressure_N = {packed[0], packed[1], packed[2]};
  value.viscous_N = {packed[3], packed[4], packed[5]};
  value.total_N = {packed[6], packed[7], packed[8]};
}

struct CaseResult final {
  immersed::ForceComponents raw;
  immersed::ForceComponents budget;
  immersed::ForceComponents surface;
  immersed::ForceComponents full_point_sum;
  immersed::WallForceSample reduced_surface;
  flow::ForceAttemptReport report;
  bool row_update_uses_positive_raw{};
  bool trace_reduced_equals_point_sum{};
  bool affine_constant_sum_conserves{};
  bool linear_gradient_exact{};
};

CaseResult run_case(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost_plan,
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const immersed::WallForceIntegrator &integrator,
    runtime::FieldRegistry &registry, runtime::FieldStorage &storage,
    const runtime::FieldAccessPlan &access, const FieldIds &ids, int ghost,
    const Preflight &preflight, const CaseInput &input) {
  fill_fields(storage, access, ids, decomposition, topology, ghost,
              preflight.wall_intercept, input);
  const auto mass_flux = finite_volume::FaceMassFlux::acquire(
      registry, storage, access, kPhase, kActor, ids.mass_flux, topology);
  const auto &constant_storage =
      static_cast<const runtime::FieldStorage &>(storage);
  const auto pressure = constant_storage.view<double>(ids.pressure);
  const auto velocity = constant_storage.view<double>(ids.velocity);
  const auto gradient = constant_storage.view<double>(ids.gradient);
  const auto cell_viscosity =
      constant_storage.view<double>(ids.cell_viscosity);
  const auto face_velocity = constant_storage.acquire_face_read<double>(
      access, kPhase, kActor, ids.face_velocity);
  const auto face_viscosity = constant_storage.acquire_face_read<double>(
      access, kPhase, kActor, ids.face_viscosity);
  auto residual = storage.view<double>(ids.residual);

  bool local_gradient_ok = true;
  int local_gradient_count = 0;
  for (const auto &link : domain.links())
    if (link.id == preflight.link) {
      ++local_gradient_count;
      const auto &reconstruction = ghost_plan.reconstruction(link.id);
      for (std::size_t component = 0U; component < 3U; ++component) {
        const auto observed = reconstruction.gradient(
            preflight.wall_intercept, velocity, component);
        const runtime::Real3 expected{
            input.velocity_gradient[component * 3U],
            input.velocity_gradient[component * 3U + 1U],
            input.velocity_gradient[component * 3U + 2U]};
        local_gradient_ok = local_gradient_ok && near(observed, expected);
      }
    }
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_gradient_count, 1, MPI_INT,
                          MPI_SUM, mpi.comm()),
            "MPI_Allreduce selected gradient count");
  const bool gradient_ok =
      all_true(local_gradient_ok, mpi,
               "MPI_Allreduce selected gradient exactness") &&
      local_gradient_count == 1;

  adapter.accumulate_momentum(mass_flux, face_velocity, velocity, pressure,
                              gradient, face_viscosity, residual);
  const auto evaluations = finite_volume::test::ImmersedOperatorTestAccess::
      last_boundary_row_evaluations(adapter);
  std::array<double, 24> row_values{};
  int local_row_count = 0;
  bool local_affine_sum_ok = true;
  for (const auto &evaluation : evaluations) {
    if (evaluation.active_cell != preflight.active_cell)
      continue;
    ++local_row_count;
    for (std::size_t component = 0U; component < 3U; ++component) {
      row_values[component] = evaluation.residual_before_wall[component];
      row_values[3U + component] =
          evaluation.background_contribution.convective[component] +
          evaluation.background_contribution.pressure[component] +
          evaluation.background_contribution.viscous[component];
      row_values[6U + component] =
          evaluation.removed_background_contribution.convective[component] +
          evaluation.removed_background_contribution.pressure[component] +
          evaluation.removed_background_contribution.viscous[component];
      row_values[9U + component] =
          evaluation.wall_contribution.pressure[component];
      row_values[12U + component] =
          evaluation.wall_contribution.viscous[component];
      row_values[15U + component] = evaluation.residual_after_wall[component];
      row_values[18U + component] =
          evaluation.budget_reaction_delta.pressure[component];
      row_values[21U + component] =
          evaluation.budget_reaction_delta.viscous[component];
    }
    std::array<double, 3> transformed_sum{};
    std::array<double, 3> background_sum{};
    for (const auto &term : evaluation.affine_donor_terms)
      if (term.input_kind ==
          finite_volume::test::BoundaryAffineInputKind::pressure)
        transformed_sum[term.output_component] += term.coefficient;
    for (const auto &term : evaluation.background_affine_donor_terms)
      if (term.input_kind ==
          finite_volume::test::BoundaryAffineInputKind::pressure)
        background_sum[term.output_component] += term.coefficient;
    for (std::size_t component = 0U; component < 3U; ++component)
      local_affine_sum_ok =
          local_affine_sum_ok &&
          near(transformed_sum[component], background_sum[component]);
  }
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_row_count, 1, MPI_INT, MPI_SUM,
                          mpi.comm()),
            "MPI_Allreduce selected evaluation count");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, row_values.data(),
                          static_cast<int>(row_values.size()), MPI_DOUBLE,
                          MPI_SUM, mpi.comm()),
            "MPI_Allreduce selected row values");
  const bool affine_sum_ok =
      all_true(local_affine_sum_ok, mpi,
               "MPI_Allreduce affine constant conservation") &&
      local_row_count == 1;

  CaseResult result;
  result.raw = force({row_values[9], row_values[10], row_values[11]},
                     {row_values[12], row_values[13], row_values[14]});
  result.budget =
      force({row_values[18], row_values[19], row_values[20]},
            {row_values[21], row_values[22], row_values[23]});
  const runtime::Real3 before{row_values[0], row_values[1], row_values[2]};
  const runtime::Real3 background{row_values[3], row_values[4], row_values[5]};
  const runtime::Real3 removed{row_values[6], row_values[7], row_values[8]};
  const runtime::Real3 after{row_values[15], row_values[16], row_values[17]};
  const auto expected_after =
      add(subtract(add(before, background), removed), result.raw.total_N);
  result.row_update_uses_positive_raw =
      local_row_count == 1 && near(after, expected_after);

  const auto trace = immersed::detail::trace_wall_force_for_test(
      integrator, pressure, velocity, gradient, cell_viscosity);
  immersed::ForceComponents local_selected{};
  immersed::ForceComponents local_full{};
  for (const auto &point : trace.local_points) {
    auto &selected_target = local_selected;
    local_full.pressure_N = add(local_full.pressure_N, point.pressure_force_N);
    local_full.viscous_N = add(local_full.viscous_N, point.viscous_force_N);
    local_full.total_N = add(local_full.total_N, point.total_force_N);
    if (point.link != preflight.link)
      continue;
    selected_target.pressure_N =
        add(selected_target.pressure_N, point.pressure_force_N);
    selected_target.viscous_N =
        add(selected_target.viscous_N, point.viscous_force_N);
    selected_target.total_N = add(selected_target.total_N, point.total_force_N);
  }
  reduce_force(local_selected, mpi, "MPI_Allreduce selected traced force");
  reduce_force(local_full, mpi, "MPI_Allreduce full traced force");
  result.surface = local_selected;
  result.full_point_sum = local_full;
  result.reduced_surface = trace.reduced;
  result.trace_reduced_equals_point_sum =
      trace.reduced.lowest_failing_rank == -1 &&
      near_force(trace.reduced.surface_traction, local_full);
  result.report =
      flow::test::ImmersedFlowTestAccess::
          assemble_force_attempt_report_from_budget(result.budget,
                                                    result.surface);
  result.affine_constant_sum_conserves = affine_sum_ok;
  result.linear_gradient_exact = gradient_ok;
  return result;
}

struct AggregateProbe final {
  immersed::ForceComponents raw;
  immersed::ForceComponents budget;
  immersed::ForceComponents adapter_report;
  bool adapter_budget_equals_negative_raw{};
};

AggregateProbe run_aggregate_probe(
    const runtime::MpiContext &mpi,
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology,
    const finite_volume::ImmersedOperatorAdapter &adapter,
    runtime::FieldRegistry &registry, runtime::FieldStorage &storage,
    const runtime::FieldAccessPlan &access, const FieldIds &ids, int ghost,
    const Preflight &preflight) {
  const CaseInput input{kPressure, {0.31, -0.17, 0.23}, {}, 0.0};
  fill_fields(storage, access, ids, decomposition, topology, ghost,
              preflight.wall_intercept, input);
  const auto mass_flux = finite_volume::FaceMassFlux::acquire(
      registry, storage, access, kPhase, kActor, ids.mass_flux, topology);
  const auto &constant_storage =
      static_cast<const runtime::FieldStorage &>(storage);
  const auto pressure = constant_storage.view<double>(ids.pressure);
  const auto velocity = constant_storage.view<double>(ids.velocity);
  const auto gradient = constant_storage.view<double>(ids.gradient);
  const auto face_velocity = constant_storage.acquire_face_read<double>(
      access, kPhase, kActor, ids.face_velocity);
  const auto face_viscosity = constant_storage.acquire_face_read<double>(
      access, kPhase, kActor, ids.face_viscosity);
  auto residual = storage.view<double>(ids.residual);
  adapter.accumulate_momentum(mass_flux, face_velocity, velocity, pressure,
                              gradient, face_viscosity, residual);
  AggregateProbe result;
  for (const auto &evaluation : finite_volume::test::ImmersedOperatorTestAccess::
           last_boundary_row_evaluations(adapter)) {
    result.raw.pressure_N =
        add(result.raw.pressure_N,
            {evaluation.wall_contribution.pressure[0],
             evaluation.wall_contribution.pressure[1],
             evaluation.wall_contribution.pressure[2]});
    result.raw.viscous_N =
        add(result.raw.viscous_N,
            {evaluation.wall_contribution.viscous[0],
             evaluation.wall_contribution.viscous[1],
             evaluation.wall_contribution.viscous[2]});
    result.budget.pressure_N =
        add(result.budget.pressure_N,
            {evaluation.budget_reaction_delta.pressure[0],
             evaluation.budget_reaction_delta.pressure[1],
             evaluation.budget_reaction_delta.pressure[2]});
    result.budget.viscous_N =
        add(result.budget.viscous_N,
            {evaluation.budget_reaction_delta.viscous[0],
             evaluation.budget_reaction_delta.viscous[1],
             evaluation.budget_reaction_delta.viscous[2]});
  }
  result.raw.total_N = add(result.raw.pressure_N, result.raw.viscous_N);
  result.budget.total_N =
      add(result.budget.pressure_N, result.budget.viscous_N);
  const auto report = adapter.report().budget_reaction_N;
  result.adapter_report =
      force({report.pressure[0], report.pressure[1], report.pressure[2]},
            {report.viscous[0], report.viscous[1], report.viscous[2]});
  const bool local_match = near_force(result.adapter_report, result.budget) &&
                           near_force(result.budget, negate_force(result.raw));
  result.adapter_budget_equals_negative_raw =
      all_true(local_match, mpi, "MPI_Allreduce aggregate budget relation");
  reduce_force(result.raw, mpi, "MPI_Allreduce aggregate raw");
  reduce_force(result.budget, mpi, "MPI_Allreduce aggregate budget");
  reduce_force(result.adapter_report, mpi,
               "MPI_Allreduce aggregate adapter report");
  return result;
}

bool near_moment(const immersed::MomentComponents &actual,
                 const immersed::MomentComponents &expected) noexcept {
  return near(actual.pressure_N_m, expected.pressure_N_m) &&
         near(actual.viscous_N_m, expected.viscous_N_m) &&
         near(actual.total_N_m, expected.total_N_m);
}

bool failed_wall_sample_is_zero(const immersed::WallForceSample &sample,
                                int expected_rank) noexcept {
  return sample.lowest_failing_rank == expected_rank &&
         sample.quadrature_point_count == 0U &&
         near_force(sample.surface_traction, {}) &&
         near_moment(sample.moment_about_global_origin, {}) &&
         near(sample.area_vector_closure_m2, {});
}

struct AuthorityConsumerObservations final {
  int fallback_lowest_failing_rank{-1};
  int supplied_lowest_failing_rank{-1};
  int wrong_provider_lowest_failing_rank{-1};
  int triangle_zero_link77_point_count{};
  int triangle_zero_execution_owner{-1};
  int link77_authority_owner{-1};
  int link77_catalog_provider{-1};
  std::uint64_t supplied_catalog_count{};
  std::uint64_t expected_catalog_count{};
  bool catalog_complete{};
  bool triangle_zero_execution_consumes_link77{};
  bool link77_authority_remains_rank0{};
  bool fallback_success{};
  bool supplied_success{};
  bool force_matches{};
  bool moment_matches{};
  bool point_count_matches{};
  bool area_closure_matches{};
  bool exactly_once_point_count{};
  bool wrong_provider_fails_on_rank0{};

  bool valid() const noexcept {
    return catalog_complete && triangle_zero_execution_consumes_link77 &&
           link77_authority_remains_rank0 && fallback_success &&
           supplied_success && force_matches && moment_matches &&
           point_count_matches && area_closure_matches &&
           exactly_once_point_count && wrong_provider_fails_on_rank0;
  }
};

AuthorityConsumerObservations
run_authority_consumer(const runtime::MpiContext &mpi,
                       const runtime::StructuredDecomposition &decomposition,
                       const mesh::MeshTopology &topology,
                       const immersed::ImmersedSurface &surface,
                       const immersed::WallQuadraturePlan &wall_plan,
                       const immersed::WallForceIntegrator &integrator,
                       runtime::FieldStorage &storage,
                       const runtime::FieldAccessPlan &access,
                       const FieldIds &ids, int halo_width) {
  const CaseInput constant_pressure{kPressure, {}, {}, 0.0};
  fill_fields(storage, access, ids, decomposition, topology, halo_width, {},
              constant_pressure);
  const auto &constant_storage =
      static_cast<const runtime::FieldStorage &>(storage);
  const auto pressure = constant_storage.view<double>(ids.pressure);
  const auto velocity = constant_storage.view<double>(ids.velocity);
  const auto gradient = constant_storage.view<double>(ids.gradient);
  const auto viscosity = constant_storage.view<double>(ids.cell_viscosity);

  AuthorityConsumerObservations result;
  const auto fallback =
      integrator.integrate(pressure, velocity, gradient, viscosity);

  int local_catalog_root =
      wall_plan.local_points().empty() ? mpi.size() : mpi.rank();
  int catalog_root = mpi.size();
  check_mpi(MPI_Allreduce(&local_catalog_root, &catalog_root, 1, MPI_INT,
                          MPI_MIN, mpi.comm()),
            "MPI_Allreduce pressure authority catalog root");
  if (catalog_root < 0 || catalog_root >= mpi.size())
    throw runtime::Error("Task 11 pressure authority catalog is empty");
  int catalog_count = 0;
  const std::vector<immersed::detail::BoundaryAuthorityOwner> *root_catalog =
      nullptr;
  if (mpi.rank() == catalog_root) {
    root_catalog = &immersed::detail::boundary_authority_catalog(
        wall_plan.local_points().front().reconstruction);
    catalog_count = static_cast<int>(root_catalog->size());
  }
  check_mpi(MPI_Bcast(&catalog_count, 1, MPI_INT, catalog_root, mpi.comm()),
            "MPI_Bcast pressure authority catalog count");
  if (catalog_count <= 0)
    throw runtime::Error("Task 11 pressure authority catalog is invalid");
  std::vector<std::uint64_t> catalog_links(
      static_cast<std::size_t>(catalog_count));
  std::vector<int> catalog_owners(static_cast<std::size_t>(catalog_count));
  if (mpi.rank() == catalog_root) {
    for (int index = 0; index < catalog_count; ++index) {
      const auto &entry = (*root_catalog)[static_cast<std::size_t>(index)];
      catalog_links[static_cast<std::size_t>(index)] = entry.link;
      catalog_owners[static_cast<std::size_t>(index)] = entry.owner_rank;
    }
  }
  check_mpi(MPI_Bcast(catalog_links.data(), catalog_count, MPI_UINT64_T,
                      catalog_root, mpi.comm()),
            "MPI_Bcast pressure authority catalog links");
  check_mpi(MPI_Bcast(catalog_owners.data(), catalog_count, MPI_INT,
                      catalog_root, mpi.comm()),
            "MPI_Bcast pressure authority catalog owners");
  bool local_catalog_consistent = true;
  for (const auto &point : wall_plan.local_points()) {
    const auto &catalog =
        immersed::detail::boundary_authority_catalog(point.reconstruction);
    if (catalog.size() != catalog_links.size()) {
      local_catalog_consistent = false;
      continue;
    }
    for (std::size_t index = 0U; index < catalog.size(); ++index)
      local_catalog_consistent =
          local_catalog_consistent &&
          catalog[index].link == catalog_links[index] &&
          catalog[index].owner_rank == catalog_owners[index];
  }

  std::vector<immersed::detail::WallPressureNormalGradient> supplied;
  supplied.reserve(catalog_links.size());
  int local_link77_provider_count = 0;
  for (std::size_t index = 0U; index < catalog_links.size(); ++index) {
    if (catalog_owners[index] != mpi.rank())
      continue;
    supplied.push_back({catalog_links[index], 0.0});
    if (catalog_links[index] == kMixedAuthorityLink)
      ++local_link77_provider_count;
  }
  std::uint64_t local_supplied_count = supplied.size();
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_supplied_count, 1, MPI_UINT64_T,
                          MPI_SUM, mpi.comm()),
            "MPI_Allreduce supplied pressure authority count");
  result.supplied_catalog_count = local_supplied_count;
  result.expected_catalog_count = catalog_links.size();
  result.catalog_complete =
      all_true(local_catalog_consistent, mpi,
               "MPI_Allreduce pressure authority catalog consistency") &&
      result.expected_catalog_count > 0U &&
      result.supplied_catalog_count == result.expected_catalog_count;

  const auto supplied_sample =
      immersed::detail::integrate_with_wall_pressure_authority(
          integrator, pressure, velocity, gradient, viscosity, supplied);

  int local_link77_point_count = 0;
  int local_execution_owner_min = mpi.size();
  int local_execution_owner_max = -1;
  int local_authority_owner_min = mpi.size();
  int local_authority_owner_max = -1;
  for (const auto &point : wall_plan.local_points()) {
    if (point.triangle != 0U ||
        immersed::detail::boundary_authority_link(point.reconstruction) !=
            kMixedAuthorityLink)
      continue;
    ++local_link77_point_count;
    local_execution_owner_min =
        std::min(local_execution_owner_min, point.owner_rank);
    local_execution_owner_max =
        std::max(local_execution_owner_max, point.owner_rank);
    const int authority_owner =
        immersed::detail::boundary_authority_owner_rank(point.reconstruction);
    local_authority_owner_min =
        std::min(local_authority_owner_min, authority_owner);
    local_authority_owner_max =
        std::max(local_authority_owner_max, authority_owner);
  }
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_link77_point_count, 1, MPI_INT,
                          MPI_SUM, mpi.comm()),
            "MPI_Allreduce triangle zero link 77 point count");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_execution_owner_min, 1, MPI_INT,
                          MPI_MIN, mpi.comm()),
            "MPI_Allreduce triangle zero execution owner minimum");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_execution_owner_max, 1, MPI_INT,
                          MPI_MAX, mpi.comm()),
            "MPI_Allreduce triangle zero execution owner maximum");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_authority_owner_min, 1, MPI_INT,
                          MPI_MIN, mpi.comm()),
            "MPI_Allreduce link 77 authority owner minimum");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_authority_owner_max, 1, MPI_INT,
                          MPI_MAX, mpi.comm()),
            "MPI_Allreduce link 77 authority owner maximum");
  int local_provider_min =
      local_link77_provider_count == 0 ? mpi.size() : mpi.rank();
  int local_provider_max = local_link77_provider_count == 0 ? -1 : mpi.rank();
  int global_link77_provider_count = local_link77_provider_count;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &global_link77_provider_count, 1,
                          MPI_INT, MPI_SUM, mpi.comm()),
            "MPI_Allreduce link 77 provider count");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_provider_min, 1, MPI_INT,
                          MPI_MIN, mpi.comm()),
            "MPI_Allreduce link 77 provider minimum");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_provider_max, 1, MPI_INT,
                          MPI_MAX, mpi.comm()),
            "MPI_Allreduce link 77 provider maximum");

  result.triangle_zero_link77_point_count = local_link77_point_count;
  if (local_execution_owner_min == local_execution_owner_max)
    result.triangle_zero_execution_owner = local_execution_owner_min;
  if (local_authority_owner_min == local_authority_owner_max)
    result.link77_authority_owner = local_authority_owner_min;
  if (local_provider_min == local_provider_max)
    result.link77_catalog_provider = local_provider_min;
  result.triangle_zero_execution_consumes_link77 =
      result.triangle_zero_link77_point_count == 2 &&
      result.triangle_zero_execution_owner == 1;
  result.link77_authority_remains_rank0 = result.link77_authority_owner == 0 &&
                                          result.link77_catalog_provider == 0 &&
                                          global_link77_provider_count == 1;

  auto wrong_provider = supplied;
  if (mpi.rank() == 0)
    wrong_provider.erase(
        std::remove_if(wrong_provider.begin(), wrong_provider.end(),
                       [](const auto &datum) {
                         return datum.link == kMixedAuthorityLink;
                       }),
        wrong_provider.end());
  if (mpi.rank() == 1)
    wrong_provider.push_back({kMixedAuthorityLink, 0.0});
  const auto wrong_provider_sample =
      immersed::detail::integrate_with_wall_pressure_authority(
          integrator, pressure, velocity, gradient, viscosity, wrong_provider);

  result.fallback_lowest_failing_rank = fallback.lowest_failing_rank;
  result.supplied_lowest_failing_rank = supplied_sample.lowest_failing_rank;
  result.wrong_provider_lowest_failing_rank =
      wrong_provider_sample.lowest_failing_rank;
  result.fallback_success = fallback.lowest_failing_rank == -1;
  result.supplied_success = supplied_sample.lowest_failing_rank == -1;
  result.force_matches =
      near_force(supplied_sample.surface_traction, fallback.surface_traction);
  result.moment_matches =
      near_moment(supplied_sample.moment_about_global_origin,
                  fallback.moment_about_global_origin);
  result.point_count_matches =
      supplied_sample.quadrature_point_count == fallback.quadrature_point_count;
  result.area_closure_matches = near(supplied_sample.area_vector_closure_m2,
                                     fallback.area_vector_closure_m2);
  const auto expected_point_count =
      3U * static_cast<std::uint64_t>(surface.triangle_count());
  result.exactly_once_point_count =
      supplied_sample.quadrature_point_count == expected_point_count;
  result.wrong_provider_fails_on_rank0 =
      failed_wall_sample_is_zero(wrong_provider_sample, 0);

  if (mpi.rank() == 0) {
    std::cout << "TASK11_A4_AUTHORITY_CONSUMER triangle=0 link="
              << kMixedAuthorityLink << " execution_link_point_count="
              << result.triangle_zero_link77_point_count
              << " execution_owner=" << result.triangle_zero_execution_owner
              << " authority_owner=" << result.link77_authority_owner
              << " catalog_provider=" << result.link77_catalog_provider
              << " supplied_catalog_count=" << result.supplied_catalog_count
              << " expected_catalog_count=" << result.expected_catalog_count
              << " fallback_lowest_failing_rank="
              << result.fallback_lowest_failing_rank
              << " supplied_lowest_failing_rank="
              << result.supplied_lowest_failing_rank
              << " wrong_provider_lowest_failing_rank="
              << result.wrong_provider_lowest_failing_rank
              << " fallback_success=" << std::boolalpha
              << result.fallback_success
              << " supplied_success=" << result.supplied_success
              << " force_match=" << result.force_matches
              << " moment_match=" << result.moment_matches
              << " point_count_match=" << result.point_count_matches
              << " exactly_once=" << result.exactly_once_point_count
              << " area_closure_match=" << result.area_closure_matches
              << " wrong_provider_lowest_rank0="
              << result.wrong_provider_fails_on_rank0 << '\n';
    std::cout.flush();
  }
  return result;
}

immersed::ForceComponents pressure_only(
    const immersed::ForceComponents &value) noexcept {
  return force(value.pressure_N, {});
}

immersed::ForceComponents viscous_only(
    const immersed::ForceComponents &value) noexcept {
  return force({}, value.viscous_N);
}

void print_force(const char *name, const immersed::ForceComponents &value) {
  std::cout << name << ".pressure=";
  print_vector(std::cout, value.pressure_N);
  std::cout << ' ' << name << ".viscous=";
  print_vector(std::cout, value.viscous_N);
  std::cout << ' ' << name << ".total=";
  print_vector(std::cout, value.total_N);
  std::cout << '\n';
}

void run(const runtime::MpiContext &mpi, const std::string &mode) {
  const auto grid = process_grid(mpi.size());
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, kExtent, {false, false, false}, runtime::DecompositionOptions{grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}});
  auto boundaries = boundary::BoundaryRegistry::create(
      case_config(mpi.size(), grid), topology);
  const auto fixture = signed_force_fixture();
  print_pre_product_fixture_observations(fixture, mpi);
  const bool pre_product_fixture_valid = all_true(
      fixture.valid(), mpi, "MPI_Allreduce pre-product fixture observations");
  HUNDUN_CHECK(pre_product_fixture_valid);
  const auto &triangles = fixture.emitted;
  FixtureFile file(mpi, triangles);
  const auto surface =
      immersed::ImmersedSurface::load_collective(file.path(), 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto wall_plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  const int ghost_reach = static_cast<int>(ghost_plan.maximum_halo_reach());
  const int wall_reach = static_cast<int>(wall_plan.maximum_halo_reach());
  HUNDUN_CHECK(ghost_reach == kWallReach);
  HUNDUN_CHECK(wall_reach == kWallReach);
  const int halo_width = std::max(ghost_reach, wall_reach);
  auto halo = runtime::HaloExchange::create(
      decomposition,
      runtime::ExchangePlan::create(decomposition, decomposition.local_extent(),
                                    halo_width));
  HUNDUN_CHECK(halo.ghost_width() == halo_width);
  auto reconstruction = finite_volume::ImmersedReconstruction::create(
      topology, geometry, boundaries, domain, ghost_plan, decomposition, mpi,
      halo);
  const immersed::LocalFlowPatternTransform transform;
  auto adapter = finite_volume::ImmersedOperatorAdapter::create(
      topology, geometry, domain, ghost_plan, transform, reconstruction);
  const auto plan_fixture = plan_fixture_observations(
      mpi, decomposition, surface, wall_plan, triangles.size());
  print_plan_fixture_observations(plan_fixture, mpi);
  HUNDUN_CHECK(plan_fixture.valid());
  const auto preflight = select_preflight(mpi, topology, geometry, domain,
                                          adapter, wall_plan, triangles);
  print_preflight(preflight, mpi);
  if (!preflight.valid)
    throw runtime::Error(
        "Task 11 RED-S1 fixture preflight found no eligible global link");
  if (mode == "fixture_preflight")
    return;
  if (mode != "authority_consumer" && mode != "behavioral_red")
    throw runtime::Error("unknown Task 11 RED-S1 test mode");
  if (mode == "authority_consumer" && mpi.size() != 2)
    throw runtime::Error(
        "Task 11 authority consumer requires exactly two MPI ranks");

  const auto integrator = immersed::WallForceIntegrator::create(wall_plan, mpi);
  runtime::FieldRegistry registry;
  const auto ids = declare_fields(registry, halo_width);
  registry.freeze();
  runtime::FieldStorage storage(
      registry, runtime::FieldLayoutSet{decomposition.local_extent(),
                                        topology.local_face_count()});
  runtime::FieldAccessPlan access(registry);
  const std::array<runtime::FieldId, 8> fields{
      ids.pressure,       ids.velocity,       ids.gradient,
      ids.cell_viscosity, ids.face_velocity,  ids.face_viscosity,
      ids.mass_flux,      ids.residual};
  for (const auto field : fields)
    access.declare_access(kPhase, kActor, field, runtime::AccessMode::read_write);
  access.freeze();

  if (mode == "authority_consumer") {
    const auto authority =
        run_authority_consumer(mpi, decomposition, topology, surface, wall_plan,
                               integrator, storage, access, ids, halo_width);
    HUNDUN_CHECK(authority.valid());
    return;
  }

  const auto stress = analytic_stress(preflight.solid_to_fluid_normal);
  const CaseInput zero_input{};
  const CaseInput pressure_input{kPressure, {}, {}, 0.0};
  const CaseInput viscous_input{0.0, {}, stress.gradient, kViscosity};
  const CaseInput combined_input{kPressure, {}, stress.gradient, kViscosity};
  const auto zero = run_case(mpi, decomposition, topology, domain, ghost_plan,
                             adapter, integrator, registry, storage, access,
                             ids, halo_width, preflight, zero_input);
  const auto pressure = run_case(
      mpi, decomposition, topology, domain, ghost_plan, adapter, integrator,
      registry, storage, access, ids, halo_width, preflight, pressure_input);
  const auto viscous = run_case(
      mpi, decomposition, topology, domain, ghost_plan, adapter, integrator,
      registry, storage, access, ids, halo_width, preflight, viscous_input);
  const auto combined = run_case(
      mpi, decomposition, topology, domain, ghost_plan, adapter, integrator,
      registry, storage, access, ids, halo_width, preflight, combined_input);
  const auto aggregate =
      run_aggregate_probe(mpi, decomposition, topology, adapter, registry,
                          storage, access, ids, halo_width, preflight);

  const auto expected_operator_pressure =
      pressure_only(preflight.expected.operator_force);
  const auto expected_operator_viscous =
      viscous_only(preflight.expected.operator_force);
  const auto expected_surface_pressure =
      pressure_only(preflight.expected.surface_force);
  const auto expected_surface_viscous =
      viscous_only(preflight.expected.surface_force);
  const auto expected_consistency_pressure = subtract_force(
      expected_operator_pressure, expected_surface_pressure);
  const auto expected_consistency_viscous = subtract_force(
      expected_operator_viscous, expected_surface_viscous);
  const auto expected_consistency_combined = subtract_force(
      preflight.expected.operator_force, preflight.expected.surface_force);

  const auto traction =
      multiply(stress.tau, preflight.solid_to_fluid_normal);
  const auto normal_traction =
      multiply(dot(traction, preflight.solid_to_fluid_normal),
               preflight.solid_to_fluid_normal);
  const auto tangential_traction = subtract(traction, normal_traction);
  const double traction_floor =
      kArithmeticFactor * std::numeric_limits<double>::epsilon();

  const double perturb_scale =
      std::max(1.0, max_abs(preflight.expected.operator_force.total_N));
  auto operator_perturb_input = combined.budget;
  const runtime::Real3 operator_delta =
      multiply(perturb_scale, {0.125, -0.25, 0.375});
  operator_perturb_input.pressure_N =
      add(operator_perturb_input.pressure_N, operator_delta);
  operator_perturb_input.total_N =
      add(operator_perturb_input.pressure_N,
          operator_perturb_input.viscous_N);
  const auto operator_perturbed =
      flow::test::ImmersedFlowTestAccess::
          assemble_force_attempt_report_from_budget(operator_perturb_input,
                                                    combined.surface);
  auto surface_perturb_input = combined.surface;
  const runtime::Real3 surface_delta =
      multiply(perturb_scale, {-0.30, 0.20, 0.10});
  surface_perturb_input.viscous_N =
      add(surface_perturb_input.viscous_N, surface_delta);
  surface_perturb_input.total_N =
      add(surface_perturb_input.pressure_N,
          surface_perturb_input.viscous_N);
  const auto surface_perturbed =
      flow::test::ImmersedFlowTestAccess::
          assemble_force_attempt_report_from_budget(combined.budget,
                                                    surface_perturb_input);

  std::vector<std::pair<std::string, bool>> observations;
  const auto observe = [&](std::string name, bool value) {
    observations.emplace_back(std::move(name), value);
  };
  observe("fixture.children_finite_positive",
          fixture.children_finite_positive);
  observe("fixture.child_normal_orientation",
          fixture.child_normal_orientation);
  observe("fixture.no_strict_process_plane_crossing",
          fixture.no_strict_process_plane_crossing);
  observe("fixture.per_source_scalar_area_preserved",
          fixture.per_source_scalar_area_preserved);
  observe("fixture.per_source_oriented_area_preserved",
          fixture.per_source_oriented_area_preserved);
  observe("fixture.total_scalar_area_preserved",
          fixture.total_scalar_area_preserved);
  observe("fixture.total_oriented_area_preserved",
          fixture.total_oriented_area_preserved);
  observe("fixture.aabb_unchanged", fixture.aabb_unchanged);
  observe("fixture.bitwise_edge_incidence_closed",
          fixture.bitwise_edge_incidence_closed);
  observe("fixture.loaded_triangle_count_matches_emitted",
          plan_fixture.loaded_triangle_count_matches_emitted);
  observe("fixture.local_point_owner_matches_rank",
          plan_fixture.local_point_owner_matches_rank);
  observe("fixture.every_triangle_has_three_points",
          plan_fixture.every_triangle_has_three_points);
  observe("fixture.every_triangle_point_index_mask_exact",
          plan_fixture.every_triangle_point_index_mask_exact);
  observe("fixture.every_triangle_has_one_owner_rank",
          plan_fixture.every_triangle_has_one_owner_rank);
  observe("fixture.selected_one_link_row", preflight.selected_one_link_row);
  observe("fixture.oblique_normal", preflight.oblique_normal);
  observe("fixture.surface_points_present", preflight.surface_points_present);
  observe("fixture.planar_link_patch", preflight.planar_link_patch);
  observe("fixture.domain_normal_equals_triangle_normal",
          preflight.domain_normal_equals_triangle_normal);
  observe("fixture.wall_intercept_on_triangle_plane",
          preflight.wall_intercept_on_triangle_plane);
  observe("fixture.quadrature_geometry_equals_fixture",
          preflight.quadrature_geometry_equals_fixture);
  observe("fixture.lfp_full_affine_row_sum_conserves_constant",
          pressure.affine_constant_sum_conserves);
  observe("fixture.lfp_source_is_zero",
          near_force(zero.raw, {}) && near_force(zero.budget, {}));
  observe("fixture.linear_gradient_reconstruction_exact",
          viscous.linear_gradient_exact);
  observe("fixture.constant_stress_defects_cancel",
          near_force(pressure.raw, expected_operator_pressure) &&
              near_force(viscous.raw, expected_operator_viscous));
  observe("fixture.operator_measure_nonzero",
          preflight.operator_measure_nonzero);
  observe("fixture.surface_measure_nonzero",
          preflight.surface_measure_nonzero);
  observe("fixture.operator_surface_measures_distinct",
          preflight.operator_surface_measures_distinct);
  observe("trace.reduced_equals_point_sum",
          combined.trace_reduced_equals_point_sum);

  observe("pressure.raw_equals_analytic_operator",
          near_force(pressure.raw, expected_operator_pressure));
  observe("pressure.row_update_uses_positive_raw_wall_term",
          pressure.row_update_uses_positive_raw);
  observe("pressure.surface_equals_analytic_surface",
          near_force(pressure.surface, expected_surface_pressure));
  observe("pressure.budget_equals_negative_raw",
          near_force(pressure.report.budget_reaction, pressure.budget) &&
              near_force(pressure.budget, negate_force(pressure.raw)));
  observe("pressure.physical_report_equals_analytic_operator",
          near_force(pressure.report.operator_force,
                     expected_operator_pressure));
  observe("pressure.consistency_equals_operator_minus_surface",
          near_force(pressure.report.consistency,
                     expected_consistency_pressure));
  observe("pressure.isolated_budget_closure",
          near_force(add_force(pressure.budget, pressure.raw), {}) &&
              near(pressure.raw.viscous_N, {}));

  observe("viscous.normal_traction_nonzero",
          norm(normal_traction) > traction_floor);
  observe("viscous.tangential_traction_nonzero",
          norm(tangential_traction) > traction_floor);
  observe("viscous.raw_equals_analytic_operator",
          near_force(viscous.raw, expected_operator_viscous));
  observe("viscous.row_update_uses_positive_raw_wall_term",
          viscous.row_update_uses_positive_raw);
  observe("viscous.surface_equals_analytic_surface",
          near_force(viscous.surface, expected_surface_viscous));
  observe("viscous.budget_equals_negative_raw",
          near_force(viscous.report.budget_reaction, viscous.budget) &&
              near_force(viscous.budget, negate_force(viscous.raw)));
  observe("viscous.physical_report_equals_analytic_operator",
          near_force(viscous.report.operator_force,
                     expected_operator_viscous));
  observe("viscous.consistency_equals_operator_minus_surface",
          near_force(viscous.report.consistency,
                     expected_consistency_viscous));
  observe("viscous.isolated_budget_closure",
          near_force(add_force(viscous.budget, viscous.raw), {}) &&
              near(viscous.raw.pressure_N, {}));

  observe("combined.operator_pressure_viscous_nonparallel",
          norm(cross(preflight.expected.operator_force.pressure_N,
                     preflight.expected.operator_force.viscous_N)) >
              traction_floor);
  observe("combined.surface_pressure_viscous_nonparallel",
          norm(cross(preflight.expected.surface_force.pressure_N,
                     preflight.expected.surface_force.viscous_N)) >
              traction_floor);
  observe("combined.raw_equals_analytic_operator",
          near_force(combined.raw, preflight.expected.operator_force));
  observe("combined.surface_equals_analytic_surface",
          near_force(combined.surface, preflight.expected.surface_force));
  observe("combined.physical_report_equals_analytic_operator",
          near_force(combined.report.operator_force,
                     preflight.expected.operator_force));
  observe("combined.consistency_equals_operator_minus_surface",
          near_force(combined.report.consistency,
                     expected_consistency_combined));

  observe("aggregate.raw_wall_sum_nonzero",
          max_abs(aggregate.raw.total_N) > traction_floor &&
              max_abs(aggregate.budget.total_N) > traction_floor);
  observe("aggregate.adapter_budget_equals_negative_raw_sum",
          aggregate.adapter_budget_equals_negative_raw &&
              near_force(aggregate.adapter_report, aggregate.budget) &&
              near_force(aggregate.budget, negate_force(aggregate.raw)));
  observe("perturb.operator_only_isolated",
          near_force(operator_perturbed.surface_traction,
                     combined.report.surface_traction) &&
              !near_force(operator_perturbed.budget_reaction,
                          combined.report.budget_reaction) &&
              !near_force(operator_perturbed.operator_force,
                          combined.report.operator_force));
  observe("perturb.surface_only_isolated",
          near_force(surface_perturbed.budget_reaction,
                     combined.report.budget_reaction) &&
              near_force(surface_perturbed.operator_force,
                         combined.report.operator_force) &&
              !near_force(surface_perturbed.surface_traction,
                          combined.report.surface_traction));
  observe("perturb.assembly_depends_on_operator",
          !near_force(operator_perturbed.budget_reaction,
                      combined.report.budget_reaction) &&
              !near_force(operator_perturbed.operator_force,
                          combined.report.operator_force));
  observe("perturb.assembly_depends_on_surface",
          !near_force(surface_perturbed.surface_traction,
                      combined.report.surface_traction));
  observe("components.pressure_viscous_separate",
          near(combined.raw.pressure_N, pressure.raw.pressure_N) &&
              near(combined.raw.viscous_N, viscous.raw.viscous_N) &&
              near(combined.surface.pressure_N, pressure.surface.pressure_N) &&
              near(combined.surface.viscous_N, viscous.surface.viscous_N) &&
              near_force(combined.report.surface_traction, combined.surface) &&
              near_force(combined.report.operator_force, combined.raw) &&
              near_force(combined.report.budget_reaction, combined.budget));

  bool aggregate_ok = true;
  for (const auto &[name, value] : observations)
    aggregate_ok = aggregate_ok && value;
  if (mpi.rank() == 0) {
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    print_force("analytic.operator", preflight.expected.operator_force);
    print_force("analytic.surface", preflight.expected.surface_force);
    print_force("pressure.raw", pressure.raw);
    print_force("pressure.budget", pressure.budget);
    print_force("pressure.surface", pressure.surface);
    print_force("pressure.report.operator", pressure.report.operator_force);
    print_force("pressure.report.consistency", pressure.report.consistency);
    print_force("viscous.raw", viscous.raw);
    print_force("viscous.budget", viscous.budget);
    print_force("viscous.surface", viscous.surface);
    print_force("viscous.report.operator", viscous.report.operator_force);
    print_force("viscous.report.consistency", viscous.report.consistency);
    print_force("combined.raw", combined.raw);
    print_force("combined.budget", combined.budget);
    print_force("combined.surface", combined.surface);
    print_force("combined.report.operator", combined.report.operator_force);
    print_force("combined.report.consistency", combined.report.consistency);
    print_force("aggregate.raw", aggregate.raw);
    print_force("aggregate.budget", aggregate.budget);
    print_force("aggregate.adapter_report", aggregate.adapter_report);
    for (const auto &[name, value] : observations)
      std::cout << name << '=' << std::boolalpha << value << '\n';
    std::cout << "aggregate.all_named_observations=" << std::boolalpha
              << aggregate_ok << '\n';
  }
  HUNDUN_CHECK(aggregate_ok);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    const std::string mode = argc > 1 ? argv[1] : "behavioral_red";
    run(mpi, mode);
  });
}
