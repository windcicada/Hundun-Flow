// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/mesh/uniform_structured_mesh.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::mesh::AnalyticWarpedBoxMapping;
using hundun::mesh::FaceAxis;
using hundun::mesh::FaceSide;
using hundun::mesh::LocalCellId;
using hundun::mesh::LocalFaceId;
using hundun::mesh::LogicalFace;
using hundun::mesh::MappingJacobian;
using hundun::mesh::MappingKind;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::mesh::UniformStructuredMesh;
using hundun::runtime::Box3;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr Int3 kExtent{7, 5, 4};
constexpr Int3 kSeamExtent{1, 5, 4};
constexpr Int3 kMixedScaleExtent{1, 2, 2};
constexpr Int3 kClosureScaleExtent{4, 1, 1};
constexpr Int3 kNonorthScaleExtent{7, 5, 4};
constexpr Int3 kJacobianScaleExtent{1, 10, 10};
constexpr Real3 kOrigin{-1.25, 0.375, 2.5};
constexpr Real3 kLength{2.75, 1.625, 4.5};
constexpr Real3 kAmplitude{0.02, -0.015, 0.01};
constexpr std::array<bool, 3> kPeriodic{true, false, false};

static_assert(!std::is_copy_constructible_v<MeshGeometry>);
static_assert(std::is_nothrow_move_constructible_v<MeshGeometry>);

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(Box3 lhs, Box3 rhs) {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

Real3 subtract(Real3 lhs, Real3 rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 add(Real3 lhs, Real3 rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 multiply(Real3 value, double factor) {
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(Real3 lhs, Real3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Real3 cross(Real3 lhs, Real3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double norm(Real3 value) { return std::hypot(value.x, value.y, value.z); }

struct Long3 {
  long double x{};
  long double y{};
  long double z{};
};

Long3 subtract(Long3 lhs, Long3 rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Long3 add(Long3 lhs, Long3 rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Long3 multiply(Long3 value, long double factor) {
  return {value.x * factor, value.y * factor, value.z * factor};
}

long double dot(Long3 lhs, Long3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Long3 cross(Long3 lhs, Long3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y,
          lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

long double norm(Long3 value) {
  return std::hypot(value.x, value.y, value.z);
}

Real3 narrow(Long3 value) {
  const Real3 narrowed{static_cast<double>(value.x),
                       static_cast<double>(value.y),
                       static_cast<double>(value.z)};
  HUNDUN_CHECK(std::isfinite(narrowed.x));
  HUNDUN_CHECK(std::isfinite(narrowed.y));
  HUNDUN_CHECK(std::isfinite(narrowed.z));
  return narrowed;
}

std::uint64_t bits(double value) {
  std::uint64_t result = 0;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool bitwise_equal(Real3 lhs, Real3 rhs) {
  return bits(lhs.x) == bits(rhs.x) && bits(lhs.y) == bits(rhs.y) &&
         bits(lhs.z) == bits(rhs.z);
}

void check_near(double actual, double expected, double factor = 128.0) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  HUNDUN_CHECK(std::abs(actual - expected) <=
               factor * std::numeric_limits<double>::epsilon() * scale);
}

void check_near(Real3 actual, Real3 expected, double factor = 128.0) {
  check_near(actual.x, expected.x, factor);
  check_near(actual.y, expected.y, factor);
  check_near(actual.z, expected.z, factor);
}

void check_relative(double actual, double expected, double factor = 128.0) {
  HUNDUN_CHECK(std::isfinite(actual));
  HUNDUN_CHECK(std::isfinite(expected));
  HUNDUN_CHECK(expected != 0.0);
  HUNDUN_CHECK(std::abs(actual / expected - 1.0) <=
               factor * std::numeric_limits<double>::epsilon());
}

template <class Function>
void expect_error(Function&& function) {
  bool rejected = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    rejected = !std::string(error.what()).empty();
  }
  HUNDUN_CHECK(rejected);
}

Int3 process_grid_for(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {1, 2, 1};
    case 4:
      return {1, 2, 2};
    default:
      throw hundun::runtime::Error("unsupported mesh-geometry test size");
  }
}

Real3 oracle_map(Real3 logical) {
  return {
      kOrigin.x + kLength.x * (logical.x + kAmplitude.x *
                                               std::sin(2.0 * kPi * logical.x) *
                                               std::sin(kPi * logical.y) *
                                               std::sin(kPi * logical.z)),
      kOrigin.y +
          kLength.y * (logical.y + kAmplitude.y * std::sin(kPi * logical.x) *
                                       std::sin(2.0 * kPi * logical.y) *
                                       std::sin(kPi * logical.z)),
      kOrigin.z +
          kLength.z * (logical.z + kAmplitude.z * std::sin(kPi * logical.x) *
                                       std::sin(kPi * logical.y) *
                                       std::sin(2.0 * kPi * logical.z))};
}

hundun::mesh::MappingJacobian oracle_jacobian(Real3 q) {
  hundun::mesh::MappingJacobian result{};
  result.d_xi_m = {
      kLength.x * (1.0 + 2.0 * kPi * kAmplitude.x * std::cos(2.0 * kPi * q.x) *
                             std::sin(kPi * q.y) * std::sin(kPi * q.z)),
      kLength.y * kPi * kAmplitude.y * std::cos(kPi * q.x) *
          std::sin(2.0 * kPi * q.y) * std::sin(kPi * q.z),
      kLength.z * kPi * kAmplitude.z * std::cos(kPi * q.x) *
          std::sin(kPi * q.y) * std::sin(2.0 * kPi * q.z)};
  result.d_eta_m = {
      kLength.x * kPi * kAmplitude.x * std::sin(2.0 * kPi * q.x) *
          std::cos(kPi * q.y) * std::sin(kPi * q.z),
      kLength.y * (1.0 + 2.0 * kPi * kAmplitude.y * std::sin(kPi * q.x) *
                             std::cos(2.0 * kPi * q.y) * std::sin(kPi * q.z)),
      kLength.z * kPi * kAmplitude.z * std::sin(kPi * q.x) *
          std::cos(kPi * q.y) * std::sin(2.0 * kPi * q.z)};
  result.d_zeta_m = {
      kLength.x * kPi * kAmplitude.x * std::sin(2.0 * kPi * q.x) *
          std::sin(kPi * q.y) * std::cos(kPi * q.z),
      kLength.y * kPi * kAmplitude.y * std::sin(kPi * q.x) *
          std::sin(2.0 * kPi * q.y) * std::cos(kPi * q.z),
      kLength.z * (1.0 + 2.0 * kPi * kAmplitude.z * std::sin(kPi * q.x) *
                             std::sin(kPi * q.y) * std::cos(2.0 * kPi * q.z))};
  return result;
}

std::array<LogicalFace, 6> cell_faces(Int3 cell) {
  return {LogicalFace{FaceAxis::x, {cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::x, {cell.x + 1, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, {cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::y, {cell.x, cell.y + 1, cell.z}},
          LogicalFace{FaceAxis::z, {cell.x, cell.y, cell.z}},
          LogicalFace{FaceAxis::z, {cell.x, cell.y, cell.z + 1}}};
}

std::array<Int3, 4> oracle_face_vertices(LogicalFace face) {
  const Int3 c = face.coordinate;
  switch (face.axis) {
    case FaceAxis::x:
      return {Int3{c.x, c.y, c.z}, Int3{c.x, c.y + 1, c.z},
              Int3{c.x, c.y + 1, c.z + 1}, Int3{c.x, c.y, c.z + 1}};
    case FaceAxis::y:
      return {Int3{c.x, c.y, c.z}, Int3{c.x, c.y, c.z + 1},
              Int3{c.x + 1, c.y, c.z + 1}, Int3{c.x + 1, c.y, c.z}};
    case FaceAxis::z:
      return {Int3{c.x, c.y, c.z}, Int3{c.x + 1, c.y, c.z},
              Int3{c.x + 1, c.y + 1, c.z}, Int3{c.x, c.y + 1, c.z}};
  }
  throw hundun::runtime::Error("invalid oracle face axis");
}

Real3 oracle_vertex(Int3 vertex) {
  return oracle_map(
      {static_cast<double>(vertex.x) / static_cast<double>(kExtent.x),
       static_cast<double>(vertex.y) / static_cast<double>(kExtent.y),
       static_cast<double>(vertex.z) / static_cast<double>(kExtent.z)});
}

Real3 affine_vertex(Int3 vertex) {
  return {kOrigin.x + kLength.x * static_cast<double>(vertex.x) /
                          static_cast<double>(kExtent.x),
          kOrigin.y + kLength.y * static_cast<double>(vertex.y) /
                          static_cast<double>(kExtent.y),
          kOrigin.z + kLength.z * static_cast<double>(vertex.z) /
                          static_cast<double>(kExtent.z)};
}

bool oracle_minimum_face(LogicalFace face) {
  switch (face.axis) {
    case FaceAxis::x:
      return face.coordinate.x == 0;
    case FaceAxis::y:
      return face.coordinate.y == 0;
    case FaceAxis::z:
      return face.coordinate.z == 0;
  }
  throw hundun::runtime::Error("invalid oracle face axis");
}

std::pair<Real3, Real3> oracle_face_geometry(LogicalFace face) {
  const auto logical_vertices = oracle_face_vertices(face);
  std::array<Real3, 4> vertices{};
  Real3 center{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] = oracle_vertex(logical_vertices[index]);
    center.x += vertices[index].x;
    center.y += vertices[index].y;
    center.z += vertices[index].z;
  }
  center.x *= 0.25;
  center.y *= 0.25;
  center.z *= 0.25;
  const Real3 first = cross(subtract(vertices[1], vertices[0]),
                            subtract(vertices[2], vertices[0]));
  const Real3 second = cross(subtract(vertices[2], vertices[0]),
                             subtract(vertices[3], vertices[0]));
  Real3 area{0.5 * (first.x + second.x), 0.5 * (first.y + second.y),
             0.5 * (first.z + second.z)};
  if (oracle_minimum_face(face)) {
    area = {-area.x, -area.y, -area.z};
  }
  return {center, area};
}

struct OracleCellGeometry {
  Real3 center{};
  double volume{};
};

OracleCellGeometry oracle_warped_cell(Int3 cell) {
  const std::array<Int3, 8> coordinates{
      Int3{cell.x, cell.y, cell.z},
      Int3{cell.x + 1, cell.y, cell.z},
      Int3{cell.x + 1, cell.y + 1, cell.z},
      Int3{cell.x, cell.y + 1, cell.z},
      Int3{cell.x, cell.y, cell.z + 1},
      Int3{cell.x + 1, cell.y, cell.z + 1},
      Int3{cell.x + 1, cell.y + 1, cell.z + 1},
      Int3{cell.x, cell.y + 1, cell.z + 1}};
  std::array<Real3, 8> vertices{};
  Real3 reference{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] = oracle_vertex(coordinates[index]);
    reference = add(reference, vertices[index]);
  }
  reference = multiply(reference, 0.125);

  // These outward quads are derived directly from the frozen face orders.
  // Minimum faces are reversed while retaining the fixed v0-v2 diagonal.
  constexpr std::array<std::array<std::size_t, 4>, 6> quads{
      std::array<std::size_t, 4>{0, 4, 7, 3},
      std::array<std::size_t, 4>{1, 2, 6, 5},
      std::array<std::size_t, 4>{0, 1, 5, 4},
      std::array<std::size_t, 4>{3, 7, 6, 2},
      std::array<std::size_t, 4>{0, 3, 2, 1},
      std::array<std::size_t, 4>{4, 5, 6, 7}};

  double volume = 0.0;
  Real3 first_moment{};
  const auto add_triangle = [&](Real3 a, Real3 b, Real3 c) {
    const double signed_volume =
        dot(subtract(a, reference),
            cross(subtract(b, reference), subtract(c, reference))) /
        6.0;
    const Real3 tetra_center =
        multiply(add(add(add(reference, a), b), c), 0.25);
    volume += signed_volume;
    first_moment = add(first_moment, multiply(tetra_center, signed_volume));
  };
  for (const auto quad : quads) {
    add_triangle(vertices[quad[0]], vertices[quad[1]], vertices[quad[2]]);
    add_triangle(vertices[quad[0]], vertices[quad[2]], vertices[quad[3]]);
  }
  HUNDUN_CHECK(std::isfinite(volume));
  HUNDUN_CHECK(volume > 0.0);
  return {multiply(first_moment, 1.0 / volume), volume};
}

constexpr long double kLongPi =
    3.141592653589793238462643383279502884L;

Long3 nonorth_scale_map(Long3 logical) {
  constexpr Long3 length{1.0e96L, 1.0e109L, 1.0e96L};
  constexpr Long3 amplitude{0.0L, -0.015L, 0.0L};
  return {
      length.x *
          (logical.x + amplitude.x * std::sin(2.0L * kLongPi * logical.x) *
                           std::sin(kLongPi * logical.y) *
                           std::sin(kLongPi * logical.z)),
      length.y *
          (logical.y + amplitude.y * std::sin(kLongPi * logical.x) *
                           std::sin(2.0L * kLongPi * logical.y) *
                           std::sin(kLongPi * logical.z)),
      length.z *
          (logical.z + amplitude.z * std::sin(kLongPi * logical.x) *
                           std::sin(kLongPi * logical.y) *
                           std::sin(2.0L * kLongPi * logical.z))};
}

Long3 nonorth_scale_vertex(Int3 vertex) {
  return nonorth_scale_map(
      {static_cast<long double>(vertex.x) /
           static_cast<long double>(kNonorthScaleExtent.x),
       static_cast<long double>(vertex.y) /
           static_cast<long double>(kNonorthScaleExtent.y),
       static_cast<long double>(vertex.z) /
           static_cast<long double>(kNonorthScaleExtent.z)});
}

std::pair<Long3, Long3> nonorth_scale_face_geometry(LogicalFace face) {
  const auto coordinates = oracle_face_vertices(face);
  std::array<Long3, 4> vertices{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] = nonorth_scale_vertex(coordinates[index]);
  }
  const Long3 reference = vertices[0];
  const Long3 center =
      add(reference,
          multiply(add(add(subtract(vertices[1], reference),
                           subtract(vertices[2], reference)),
                       subtract(vertices[3], reference)),
                   0.25L));
  const Long3 first = cross(subtract(vertices[1], vertices[0]),
                            subtract(vertices[2], vertices[0]));
  const Long3 second = cross(subtract(vertices[2], vertices[0]),
                             subtract(vertices[3], vertices[0]));
  return {center, multiply(add(first, second), 0.5L)};
}

struct LongCellGeometry {
  Long3 center{};
  long double volume{};
};

LongCellGeometry nonorth_scale_cell_geometry(Int3 cell) {
  const std::array<Int3, 8> coordinates{
      Int3{cell.x, cell.y, cell.z},
      Int3{cell.x + 1, cell.y, cell.z},
      Int3{cell.x + 1, cell.y + 1, cell.z},
      Int3{cell.x, cell.y + 1, cell.z},
      Int3{cell.x, cell.y, cell.z + 1},
      Int3{cell.x + 1, cell.y, cell.z + 1},
      Int3{cell.x + 1, cell.y + 1, cell.z + 1},
      Int3{cell.x, cell.y + 1, cell.z + 1}};
  std::array<Long3, 8> vertices{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] = nonorth_scale_vertex(coordinates[index]);
  }
  const Long3 anchor = vertices[0];
  Long3 relative_sum{};
  for (const Long3 vertex : vertices) {
    relative_sum = add(relative_sum, subtract(vertex, anchor));
  }
  const Long3 reference = add(anchor, multiply(relative_sum, 0.125L));

  constexpr std::array<std::array<std::size_t, 4>, 6> quads{
      std::array<std::size_t, 4>{0, 4, 7, 3},
      std::array<std::size_t, 4>{1, 2, 6, 5},
      std::array<std::size_t, 4>{0, 1, 5, 4},
      std::array<std::size_t, 4>{3, 7, 6, 2},
      std::array<std::size_t, 4>{0, 3, 2, 1},
      std::array<std::size_t, 4>{4, 5, 6, 7}};
  long double volume = 0.0L;
  Long3 relative_first_moment{};
  const auto add_triangle = [&](Long3 a, Long3 b, Long3 c) {
    const long double signed_volume =
        dot(subtract(a, reference),
            cross(subtract(b, reference), subtract(c, reference))) /
        6.0L;
    const Long3 relative_tetra_center =
        multiply(add(add(subtract(a, reference), subtract(b, reference)),
                     subtract(c, reference)),
                 0.25L);
    volume += signed_volume;
    relative_first_moment =
        add(relative_first_moment,
            multiply(relative_tetra_center, signed_volume));
  };
  for (const auto quad : quads) {
    add_triangle(vertices[quad[0]], vertices[quad[1]], vertices[quad[2]]);
    add_triangle(vertices[quad[0]], vertices[quad[2]], vertices[quad[3]]);
  }
  HUNDUN_CHECK(std::isfinite(volume));
  HUNDUN_CHECK(volume > 0.0L);
  return {add(reference, multiply(relative_first_moment, 1.0L / volume)),
          volume};
}

double oracle_minimum_sampled_jacobian(Int3 cell) {
  double minimum = std::numeric_limits<double>::infinity();
  const auto sample = [&](double x, double y, double z) {
    const Real3 logical{
        (static_cast<double>(cell.x) + x) /
            static_cast<double>(kExtent.x),
        (static_cast<double>(cell.y) + y) /
            static_cast<double>(kExtent.y),
        (static_cast<double>(cell.z) + z) /
            static_cast<double>(kExtent.z)};
    const auto jacobian = oracle_jacobian(logical);
    minimum = std::min(
        minimum,
        dot(jacobian.d_xi_m, cross(jacobian.d_eta_m, jacobian.d_zeta_m)));
  };
  for (int z = 0; z <= 1; ++z) {
    for (int y = 0; y <= 1; ++y) {
      for (int x = 0; x <= 1; ++x) {
        sample(static_cast<double>(x), static_cast<double>(y),
               static_cast<double>(z));
      }
    }
  }
  sample(0.5, 0.5, 0.5);
  const double offset = 1.0 / (2.0 * std::sqrt(3.0));
  for (int z = -1; z <= 1; z += 2) {
    for (int y = -1; y <= 1; y += 2) {
      for (int x = -1; x <= 1; x += 2) {
        sample(0.5 + static_cast<double>(x) * offset,
               0.5 + static_cast<double>(y) * offset,
               0.5 + static_cast<double>(z) * offset);
      }
    }
  }
  HUNDUN_CHECK(std::isfinite(minimum));
  HUNDUN_CHECK(minimum > 0.0);
  return minimum;
}

void test_mapping_formulas() {
  const AnalyticWarpedBoxMapping mapping(kOrigin, kLength, kAmplitude);
  const Real3 point{0.137, 0.421, 0.733};
  check_near(mapping.map(point), oracle_map(point));
  const auto actual = mapping.jacobian(point);
  const auto expected = oracle_jacobian(point);
  check_near(actual.d_xi_m, expected.d_xi_m);
  check_near(actual.d_eta_m, expected.d_eta_m);
  check_near(actual.d_zeta_m, expected.d_zeta_m);
  check_near(actual.determinant_m3(),
             dot(expected.d_xi_m, cross(expected.d_eta_m, expected.d_zeta_m)));

  const UniformBoxMapping uniform(kOrigin, kLength);
  check_near(uniform.map(point),
             {kOrigin.x + kLength.x * point.x, kOrigin.y + kLength.y * point.y,
              kOrigin.z + kLength.z * point.z});
  check_near(uniform.jacobian(point).determinant_m3(),
             kLength.x * kLength.y * kLength.z);
}

void test_range_safe_mapping_jacobians() {
  constexpr Real3 origin{0.0, 0.0, 0.0};
  constexpr Real3 length{1.0e308, 1.0e-154, 1.0e-154};
  constexpr Real3 amplitude{0.005, 0.0, 0.0};
  constexpr Real3 logical{0.25, 0.0, 0.5};
  const AnalyticWarpedBoxMapping mapping(origin, length, amplitude);
  const Real3 point = mapping.map(logical);
  HUNDUN_CHECK(std::isfinite(point.x));
  HUNDUN_CHECK(std::isfinite(point.y));
  HUNDUN_CHECK(std::isfinite(point.z));
  const MappingJacobian jacobian = mapping.jacobian(logical);
  const double expected_off_diagonal =
      length.x * ((kPi * amplitude.x) *
                  std::sin(2.0 * kPi * logical.x) *
                  std::cos(kPi * logical.y) * std::sin(kPi * logical.z));
  HUNDUN_CHECK(std::isfinite(expected_off_diagonal));
  check_relative(jacobian.d_eta_m.x, expected_off_diagonal, 32.0);
  const double expected_determinant = static_cast<double>(
      static_cast<long double>(length.x) *
      static_cast<long double>(length.y) *
      static_cast<long double>(length.z));
  check_relative(jacobian.determinant_m3(), expected_determinant, 32.0);

  const MappingJacobian overflow_first{{1.0e-210, 0.0, 0.0},
                                       {0.0, 1.0e160, 0.0},
                                       {0.0, 0.0, 1.0e150}};
  check_relative(overflow_first.determinant_m3(), 1.0e100, 32.0);

  const MappingJacobian underflow_first{{1.0e210, 0.0, 0.0},
                                        {0.0, 1.0e-200, 0.0},
                                        {0.0, 0.0, 1.0e-200}};
  check_relative(underflow_first.determinant_m3(), 1.0e-190, 32.0);

  const MappingJacobian mixed_rows{{1.0e200, 0.0, 1.0e200},
                                   {0.0, 1.0e200, 0.0},
                                   {0.0, 0.0, 1.0e-300}};
  check_relative(mixed_rows.determinant_m3(), 1.0e100, 32.0);

  const MappingJacobian mixed_columns{{1.0e200, 0.0, 0.0},
                                      {0.0, 1.0e200, 0.0},
                                      {1.0e200, 0.0, 1.0e-300}};
  check_relative(mixed_columns.determinant_m3(), 1.0e100, 32.0);

  const MappingJacobian even_column_permutation{
      mixed_rows.d_eta_m, mixed_rows.d_zeta_m, mixed_rows.d_xi_m};
  check_relative(even_column_permutation.determinant_m3(), 1.0e100, 32.0);

  const MappingJacobian odd_column_permutation{
      mixed_rows.d_eta_m, mixed_rows.d_xi_m, mixed_rows.d_zeta_m};
  check_relative(odd_column_permutation.determinant_m3(), -1.0e100, 32.0);

  constexpr double exactly_integral = 4503599627370496.0;
  const MappingJacobian cancellation{
      {exactly_integral, exactly_integral + 1.0, 0.0},
      {exactly_integral - 1.0, exactly_integral, 0.0},
      {0.0, 0.0, 1.0}};
  HUNDUN_CHECK(cancellation.determinant_m3() == 1.0);

  const MappingJacobian exact_zero{{1.0, 2.0, 3.0},
                                   {1.0, 2.0, 3.0},
                                   {4.0, 5.0, 6.0}};
  HUNDUN_CHECK(exact_zero.determinant_m3() == 0.0);

  const MappingJacobian non_finite{
      {std::numeric_limits<double>::infinity(), 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}};
  HUNDUN_CHECK(std::isnan(non_finite.determinant_m3()));

  const MappingJacobian genuine_overflow{
      {std::numeric_limits<double>::max(), 0.0, 0.0},
      {0.0, 2.0, 0.0},
      {0.0, 0.0, 1.0}};
  HUNDUN_CHECK(std::isinf(genuine_overflow.determinant_m3()));

  const MappingJacobian genuine_underflow{
      {std::numeric_limits<double>::denorm_min(), 0.0, 0.0},
      {0.0, 0.5, 0.0},
      {0.0, 0.0, 1.0}};
  HUNDUN_CHECK(genuine_underflow.determinant_m3() == 0.0);

  const MappingJacobian rounded_subnormal{
      {std::numeric_limits<double>::denorm_min(), 0.0, 0.0},
      {0.0, 0.75, 0.0},
      {0.0, 0.0, 1.0}};
  HUNDUN_CHECK(rounded_subnormal.determinant_m3() ==
               std::numeric_limits<double>::denorm_min());
}

void test_mapping_rejections() {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  expect_error([&] { UniformBoxMapping({nan, 0.0, 0.0}, kLength); });
  expect_error([&] { UniformBoxMapping(kOrigin, {0.0, 1.0, 1.0}); });
  expect_error([&] {
    UniformBoxMapping({std::numeric_limits<double>::max(), 0.0, 0.0},
                      {std::numeric_limits<double>::max(), 1.0, 1.0});
  });
  expect_error([&] {
    AnalyticWarpedBoxMapping(kOrigin, kLength, {0.0200001, 0.0, 0.0});
  });
  expect_error([&] {
    AnalyticWarpedBoxMapping(kOrigin, kLength, {0.0, nan, 0.0});
  });
  const AnalyticWarpedBoxMapping mapping(kOrigin, kLength, kAmplitude);
  expect_error([&] { static_cast<void>(mapping.map({-0.01, 0.5, 0.5})); });
  expect_error([&] { static_cast<void>(mapping.map({0.5, infinity, 0.5})); });
  expect_error([&] { static_cast<void>(mapping.jacobian({0.5, 0.5, 1.01})); });
}

void test_uniform_adapter(const StructuredDecomposition& decomposition,
                          const MeshTopology& topology,
                          const MeshGeometry& geometry) {
  const UniformStructuredMesh legacy(kExtent, kOrigin, kLength, decomposition);
  HUNDUN_CHECK(geometry.mapping_kind() == MappingKind::uniform_box);
  HUNDUN_CHECK(same(geometry.global_extent(), kExtent));
  HUNDUN_CHECK(same(geometry.owned_global_box(), legacy.owned_global_box()));
  HUNDUN_CHECK(geometry.uniform_spacing_m().has_value());
  HUNDUN_CHECK(
      bitwise_equal(*geometry.uniform_spacing_m(), legacy.spacing_m()));
  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    const Int3 global = topology.global_cell(local);
    const Box3 box = topology.owned_global_box();
    const Int3 legacy_local{global.x - box.begin.x, global.y - box.begin.y,
                            global.z - box.begin.z};
    HUNDUN_CHECK(bitwise_equal(geometry.cell_center_m(local),
                               legacy.cell_center(legacy_local)));
    HUNDUN_CHECK(bits(geometry.cell_volume_m3(local)) ==
                 bits(legacy.cell_volume_m3()));
    HUNDUN_CHECK(geometry.minimum_jacobian_determinant_m3(local) > 0.0);
  }
  for (LocalFaceId local = 0; local < topology.local_face_count(); ++local) {
    HUNDUN_CHECK(geometry.face_skewness(local) == 0.0);
    HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(local) == 0.0);
  }

  HUNDUN_CHECK(!topology.patch(0U).local_faces().empty());
  const LocalFaceId x_min = topology.patch(0U).local_faces().front();
  HUNDUN_CHECK(topology.neighbour(x_min).has_value());
  HUNDUN_CHECK(topology.owner(x_min) != *topology.neighbour(x_min));
  HUNDUN_CHECK(geometry.face_area_vector_m2(x_min, FaceSide::owner).x < 0.0);
}

struct FaceRecord {
  std::uint64_t id{};
  std::array<std::uint64_t, 6> values{};
};

static_assert(std::is_trivially_copyable_v<FaceRecord>);

void test_duplicate_face_bits(const MpiContext& mpi,
                              const MeshTopology& topology,
                              const MeshGeometry& geometry) {
  std::vector<FaceRecord> local;
  local.reserve(topology.local_face_count());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const Real3 center = geometry.face_center_m(face);
    const Real3 area = geometry.face_area_vector_m2(face, FaceSide::owner);
    local.push_back({topology.global_face_id(face),
                     {bits(center.x), bits(center.y), bits(center.z),
                      bits(area.x), bits(area.y), bits(area.z)}});
  }
  HUNDUN_CHECK(local.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                   sizeof(FaceRecord));
  const int local_bytes = static_cast<int>(local.size() * sizeof(FaceRecord));
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  HUNDUN_CHECK(MPI_Allgather(&local_bytes, 1, MPI_INT, counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> offsets(static_cast<std::size_t>(mpi.size()));
  int total = 0;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    offsets[static_cast<std::size_t>(rank)] = total;
    HUNDUN_CHECK(counts[static_cast<std::size_t>(rank)] >= 0);
    HUNDUN_CHECK(total <= std::numeric_limits<int>::max() -
                              counts[static_cast<std::size_t>(rank)]);
    total += counts[static_cast<std::size_t>(rank)];
  }
  HUNDUN_CHECK(total % static_cast<int>(sizeof(FaceRecord)) == 0);
  std::vector<FaceRecord> gathered(static_cast<std::size_t>(total) /
                                   sizeof(FaceRecord));
  HUNDUN_CHECK(MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE,
                              gathered.data(), counts.data(), offsets.data(),
                              MPI_BYTE, mpi.comm()) == MPI_SUCCESS);
  std::sort(gathered.begin(), gathered.end(),
            [](const FaceRecord& lhs, const FaceRecord& rhs) {
              return lhs.id < rhs.id;
            });
  std::size_t duplicates = 0;
  for (std::size_t index = 1; index < gathered.size(); ++index) {
    if (gathered[index - 1].id == gathered[index].id) {
      HUNDUN_CHECK(gathered[index - 1].values == gathered[index].values);
      ++duplicates;
    }
  }
  if (mpi.size() > 1) {
    HUNDUN_CHECK(duplicates > 0);
  }
}

void test_warped_invariants(const MpiContext& mpi, const MeshTopology& topology,
                            const MeshGeometry& geometry) {
  HUNDUN_CHECK(geometry.mapping_kind() == MappingKind::analytic_warped_box);
  HUNDUN_CHECK(!geometry.uniform_spacing_m().has_value());

  int local_curved_vertex = 0;
  for (int k = 1; k < kExtent.z; ++k) {
    for (int j = 1; j < kExtent.y; ++j) {
      for (int i = 1; i < kExtent.x; ++i) {
        if (norm(subtract(geometry.vertex_position_m({i, j, k}),
                          affine_vertex({i, j, k}))) > 1.0e-12) {
          local_curved_vertex = 1;
        }
      }
    }
  }
  int global_curved_vertex = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_curved_vertex, &global_curved_vertex, 1,
                             MPI_INT, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_curved_vertex == 1);

  int local_curved_cell = 0;
  for (LocalCellId local = 0; local < topology.local_cell_count(); ++local) {
    const Int3 global = topology.global_cell(local);
    const auto expected = oracle_warped_cell(global);
    check_near(geometry.cell_center_m(local), expected.center, 512.0);
    check_near(geometry.cell_volume_m3(local), expected.volume, 512.0);
    check_near(geometry.minimum_jacobian_determinant_m3(local),
               oracle_minimum_sampled_jacobian(global), 256.0);

    const Real3 affine_center{
        kOrigin.x +
            (static_cast<double>(global.x) + 0.5) * kLength.x /
                static_cast<double>(kExtent.x),
        kOrigin.y +
            (static_cast<double>(global.y) + 0.5) * kLength.y /
                static_cast<double>(kExtent.y),
        kOrigin.z +
            (static_cast<double>(global.z) + 0.5) * kLength.z /
                static_cast<double>(kExtent.z)};
    if (norm(subtract(expected.center, affine_center)) > 1.0e-12) {
      local_curved_cell = 1;
    }
  }
  int global_curved_cell = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_curved_cell, &global_curved_cell, 1,
                             MPI_INT, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_curved_cell == 1);

  double local_volume = 0.0;
  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    HUNDUN_CHECK(geometry.cell_volume_m3(local) > 0.0);
    HUNDUN_CHECK(geometry.minimum_jacobian_determinant_m3(local) > 0.0);
    local_volume += geometry.cell_volume_m3(local);

    const Real3 closure = geometry.cell_closure_m2(local);
    double area_sum = 0.0;
    for (const auto logical : cell_faces(topology.global_cell(local))) {
      const auto face =
          topology.find_local_face(topology.global_face_id(logical));
      HUNDUN_CHECK(face.has_value());
      area_sum += geometry.face_area_m2(*face);
    }
    HUNDUN_CHECK(norm(closure) <=
                 256.0 * std::numeric_limits<double>::epsilon() * area_sum);
    const Real3 constant_flux{0.375, -0.625, 1.125};
    HUNDUN_CHECK(std::abs(dot(constant_flux, closure)) <=
                 256.0 * std::numeric_limits<double>::epsilon() * area_sum *
                     norm(constant_flux));
  }
  double global_volume = 0.0;
  HUNDUN_CHECK(MPI_Allreduce(&local_volume, &global_volume, 1, MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  const double box_volume = kLength.x * kLength.y * kLength.z;
  // Internal fixed-diagonal triangles cancel in the closed global polyhedron;
  // 1024 eps covers the remaining binary64 triangle and rank-sum roundoff.
  HUNDUN_CHECK(std::abs(global_volume - box_volume) <=
               1024.0 * std::numeric_limits<double>::epsilon() * box_volume);

  int local_curved_metric = 0;
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const auto expected = oracle_face_geometry(topology.logical_face(face));
    check_near(geometry.face_center_m(face), expected.first, 256.0);
    check_near(geometry.face_area_vector_m2(face, FaceSide::owner),
               expected.second, 256.0);
    HUNDUN_CHECK(std::isfinite(geometry.face_skewness(face)));
    HUNDUN_CHECK(geometry.face_skewness(face) >= 0.0);
    HUNDUN_CHECK(std::isfinite(geometry.face_non_orthogonality_degrees(face)));
    HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(face) >= 0.0);
    HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(face) < 90.0);
    if (geometry.face_skewness(face) > 1.0e-12 ||
        geometry.face_non_orthogonality_degrees(face) > 1.0e-10) {
      local_curved_metric = 1;
    }
    if (topology.neighbour(face).has_value()) {
      const Real3 owner = geometry.face_area_vector_m2(face, FaceSide::owner);
      const Real3 neighbour =
          geometry.face_area_vector_m2(face, FaceSide::neighbour);
      HUNDUN_CHECK(bits(neighbour.x) == bits(-owner.x));
      HUNDUN_CHECK(bits(neighbour.y) == bits(-owner.y));
      HUNDUN_CHECK(bits(neighbour.z) == bits(-owner.z));
    } else {
      HUNDUN_CHECK(geometry.face_skewness(face) == 0.0);
      expect_error([&] {
        static_cast<void>(
            geometry.face_area_vector_m2(face, FaceSide::neighbour));
      });
    }
  }
  int global_curved_metric = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_curved_metric, &global_curved_metric, 1,
                             MPI_INT, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_curved_metric == 1);
  test_duplicate_face_bits(mpi, topology, geometry);
}

void test_extent_one_periodic_seam(const MpiContext& context,
                                   std::optional<MeshGeometry>& detached) {
  const DecompositionOptions options{process_grid_for(context.size())};
  auto decomposition = StructuredDecomposition::create(
      context, kSeamExtent, kPeriodic, options);
  MeshTopology topology(decomposition);

  MeshGeometry uniform(topology, UniformBoxMapping(kOrigin, kLength));
  HUNDUN_CHECK(uniform.compatible(topology));
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(uniform.face_skewness(face) == 0.0);
    HUNDUN_CHECK(uniform.face_non_orthogonality_degrees(face) == 0.0);
  }

  detached.emplace(topology,
                   AnalyticWarpedBoxMapping(kOrigin, kLength, kAmplitude));
  HUNDUN_CHECK(detached->compatible(topology));
  HUNDUN_CHECK(!topology.patch(0U).local_faces().empty());
  const LocalFaceId x_min = topology.patch(0U).local_faces().front();
  HUNDUN_CHECK(topology.neighbour(x_min).has_value());
  HUNDUN_CHECK(topology.owner(x_min) == *topology.neighbour(x_min));
  const Real3 owner =
      detached->face_area_vector_m2(x_min, FaceSide::owner);
  const Real3 neighbour =
      detached->face_area_vector_m2(x_min, FaceSide::neighbour);
  HUNDUN_CHECK(owner.x < 0.0);
  HUNDUN_CHECK(bits(neighbour.x) == bits(-owner.x));
  HUNDUN_CHECK(bits(neighbour.y) == bits(-owner.y));
  HUNDUN_CHECK(bits(neighbour.z) == bits(-owner.z));
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(std::isfinite(detached->face_skewness(face)));
    HUNDUN_CHECK(detached->face_skewness(face) >= 0.0);
    HUNDUN_CHECK(
        std::isfinite(detached->face_non_orthogonality_degrees(face)));
    HUNDUN_CHECK(detached->face_non_orthogonality_degrees(face) >= 0.0);
    HUNDUN_CHECK(detached->face_non_orthogonality_degrees(face) < 90.0);
  }
  context.barrier();
}

void test_mixed_scale_geometry(const MpiContext& context) {
  constexpr Real3 origin{5.0e307, 0.0, 0.0};
  constexpr Real3 length{8.0e292, 1.0e-150, 1.0e-150};
  const DecompositionOptions options{process_grid_for(context.size())};
  auto decomposition = StructuredDecomposition::create(
      context, kMixedScaleExtent, {false, false, false}, options);
  MeshTopology topology(decomposition);
  const UniformStructuredMesh legacy(kMixedScaleExtent, origin, length,
                                     decomposition);
  MeshGeometry uniform(topology, UniformBoxMapping(origin, length));
  MeshGeometry analytic(
      topology, AnalyticWarpedBoxMapping(origin, length, {0.0, 0.0, 0.0}));

  HUNDUN_CHECK(bitwise_equal(*uniform.uniform_spacing_m(), legacy.spacing_m()));
  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    const Int3 global = topology.global_cell(local);
    const Box3 box = topology.owned_global_box();
    const Int3 legacy_local{global.x - box.begin.x, global.y - box.begin.y,
                            global.z - box.begin.z};
    HUNDUN_CHECK(bitwise_equal(uniform.cell_center_m(local),
                               legacy.cell_center(legacy_local)));
    HUNDUN_CHECK(bits(uniform.cell_volume_m3(local)) ==
                 bits(legacy.cell_volume_m3()));
    const Real3 analytic_center = analytic.cell_center_m(local);
    HUNDUN_CHECK(std::isfinite(analytic_center.x));
    HUNDUN_CHECK(std::isfinite(analytic_center.y));
    HUNDUN_CHECK(std::isfinite(analytic_center.z));
    HUNDUN_CHECK(std::isfinite(analytic.cell_volume_m3(local)));
    HUNDUN_CHECK(analytic.cell_volume_m3(local) > 0.0);
  }
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const Real3 center = analytic.face_center_m(face);
    HUNDUN_CHECK(std::isfinite(center.x));
    HUNDUN_CHECK(std::isfinite(center.y));
    HUNDUN_CHECK(std::isfinite(center.z));
    HUNDUN_CHECK(std::isfinite(analytic.face_area_m2(face)));
    HUNDUN_CHECK(analytic.face_area_m2(face) > 0.0);
    HUNDUN_CHECK(std::isfinite(analytic.face_skewness(face)));
    HUNDUN_CHECK(
        std::isfinite(analytic.face_non_orthogonality_degrees(face)));
  }
}

void check_scaled_closure(const MeshTopology& topology,
                          const MeshGeometry& geometry,
                          LocalCellId local) {
  const Real3 closure = geometry.cell_closure_m2(local);
  HUNDUN_CHECK(std::isfinite(closure.x));
  HUNDUN_CHECK(std::isfinite(closure.y));
  HUNDUN_CHECK(std::isfinite(closure.z));

  std::array<double, 6> areas{};
  double scale = 0.0;
  std::size_t index = 0;
  for (const auto logical : cell_faces(topology.global_cell(local))) {
    const auto face =
        topology.find_local_face(topology.global_face_id(logical));
    HUNDUN_CHECK(face.has_value());
    const double area = geometry.face_area_m2(*face);
    HUNDUN_CHECK(std::isfinite(area));
    HUNDUN_CHECK(area > 0.0);
    areas[index++] = area;
    scale = std::max(scale, area);
  }
  HUNDUN_CHECK(index == areas.size());
  HUNDUN_CHECK(std::isfinite(scale));
  HUNDUN_CHECK(scale > 0.0);

  double scaled_sum = 0.0;
  for (const double area : areas) {
    scaled_sum += area / scale;
  }
  const Real3 scaled_closure{closure.x / scale, closure.y / scale,
                             closure.z / scale};
  HUNDUN_CHECK(std::isfinite(scaled_sum));
  HUNDUN_CHECK(scaled_sum > 0.0);
  HUNDUN_CHECK(std::isfinite(norm(scaled_closure)));
  HUNDUN_CHECK(norm(scaled_closure) <=
               256.0 * std::numeric_limits<double>::epsilon() * scaled_sum);
}

void check_anisotropic_metrics(const MeshTopology& topology,
                               const MeshGeometry& geometry) {
  for (LocalCellId local = 0; local < topology.local_cell_count(); ++local) {
    const Real3 center = geometry.cell_center_m(local);
    HUNDUN_CHECK(std::isfinite(center.x));
    HUNDUN_CHECK(std::isfinite(center.y));
    HUNDUN_CHECK(std::isfinite(center.z));
    HUNDUN_CHECK(std::isfinite(geometry.cell_volume_m3(local)));
    HUNDUN_CHECK(geometry.cell_volume_m3(local) > 0.0);
    HUNDUN_CHECK(
        std::isfinite(geometry.minimum_jacobian_determinant_m3(local)));
    HUNDUN_CHECK(geometry.minimum_jacobian_determinant_m3(local) > 0.0);
  }
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const Real3 center = geometry.face_center_m(face);
    HUNDUN_CHECK(std::isfinite(center.x));
    HUNDUN_CHECK(std::isfinite(center.y));
    HUNDUN_CHECK(std::isfinite(center.z));
    HUNDUN_CHECK(std::isfinite(geometry.face_area_m2(face)));
    HUNDUN_CHECK(geometry.face_area_m2(face) > 0.0);
    HUNDUN_CHECK(std::isfinite(geometry.face_skewness(face)));
    HUNDUN_CHECK(geometry.face_skewness(face) >= 0.0);
    HUNDUN_CHECK(
        std::isfinite(geometry.face_non_orthogonality_degrees(face)));
    HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(face) >= 0.0);
    HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(face) < 90.0);
  }
  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    check_scaled_closure(topology, geometry, local);
  }
}

void test_nonorthogonality_scale_geometry(const MpiContext& context) {
  constexpr Real3 origin{0.0, 0.0, 0.0};
  constexpr Real3 length{1.0e96, 1.0e109, 1.0e96};
  constexpr Real3 amplitude{0.0, -0.015, 0.0};
  const DecompositionOptions options{process_grid_for(context.size())};
  auto decomposition = StructuredDecomposition::create(
      context, kNonorthScaleExtent, {false, false, false}, options);
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, AnalyticWarpedBoxMapping(origin, length, amplitude));
  check_anisotropic_metrics(topology, geometry);

  const LogicalFace target{FaceAxis::x, {1, 0, 0}};
  const auto local_face =
      topology.find_local_face(topology.global_face_id(target));
  int local_hits = 0;
  if (local_face.has_value()) {
    local_hits = 1;
    const LocalCellId owner = topology.owner(*local_face);
    HUNDUN_CHECK(topology.neighbour(*local_face).has_value());
    const LocalCellId neighbour = *topology.neighbour(*local_face);
    HUNDUN_CHECK(same(topology.global_cell(owner), Int3{0, 0, 0}));
    HUNDUN_CHECK(same(topology.global_cell(neighbour), Int3{1, 0, 0}));

    const auto expected_face = nonorth_scale_face_geometry(target);
    const auto expected_owner = nonorth_scale_cell_geometry({0, 0, 0});
    const auto expected_neighbour = nonorth_scale_cell_geometry({1, 0, 0});
    check_near(geometry.face_center_m(*local_face),
               narrow(expected_face.first), 1024.0);
    check_near(geometry.face_area_vector_m2(*local_face, FaceSide::owner),
               narrow(expected_face.second), 1024.0);
    check_near(geometry.cell_center_m(owner), narrow(expected_owner.center),
               4096.0);
    check_near(geometry.cell_center_m(neighbour),
               narrow(expected_neighbour.center), 4096.0);
    check_near(geometry.cell_volume_m3(owner),
               static_cast<double>(expected_owner.volume), 4096.0);

    const Long3 expected_displacement =
        subtract(expected_neighbour.center, expected_owner.center);
    const long double expected_area_norm = norm(expected_face.second);
    const long double expected_displacement_norm = norm(expected_displacement);
    const long double expected_projection =
        dot(expected_face.second, expected_displacement);
    HUNDUN_CHECK(std::isfinite(expected_area_norm));
    HUNDUN_CHECK(expected_area_norm > 0.0L);
    HUNDUN_CHECK(std::isfinite(expected_displacement_norm));
    HUNDUN_CHECK(expected_displacement_norm > 0.0L);
    HUNDUN_CHECK(std::isfinite(expected_projection));
    HUNDUN_CHECK(expected_projection > 0.0L);
    const long double expected_cosine = std::clamp(
        expected_projection /
            (expected_area_norm * expected_displacement_norm),
        -1.0L, 1.0L);
    HUNDUN_CHECK(std::isfinite(expected_cosine));
    const double expected_angle = static_cast<double>(
        std::acos(expected_cosine) * 180.0L / kLongPi);
    HUNDUN_CHECK(std::isfinite(expected_angle));
    HUNDUN_CHECK(expected_angle > 89.999999);
    HUNDUN_CHECK(expected_angle < 90.0);

    const Real3 displacement =
        subtract(geometry.cell_center_m(neighbour),
                 geometry.cell_center_m(owner));
    const double area = geometry.face_area_m2(*local_face);
    const double displacement_norm = norm(displacement);
    const double projection = dot(
        geometry.face_area_vector_m2(*local_face, FaceSide::owner),
        displacement);
    HUNDUN_CHECK(std::isfinite(area));
    HUNDUN_CHECK(area > 0.0);
    HUNDUN_CHECK(std::isfinite(displacement_norm));
    HUNDUN_CHECK(displacement_norm > 0.0);
    HUNDUN_CHECK(std::isfinite(projection));
    HUNDUN_CHECK(projection > 0.0);
    // The requested metric remains representable even though this ordinary
    // dimensional intermediate is outside binary64's range.
    HUNDUN_CHECK(!std::isfinite(area * displacement_norm));
    const double actual_angle =
        geometry.face_non_orthogonality_degrees(*local_face);
    HUNDUN_CHECK(actual_angle > 89.999999);
    HUNDUN_CHECK(actual_angle < 90.0);
    check_near(actual_angle, expected_angle, 1024.0);
  }

  int global_hits = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_hits, &global_hits, 1, MPI_INT, MPI_SUM,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_hits > 0);
}

void test_closure_scale_geometry(const MpiContext& context) {
  constexpr Real3 origin{0.0, 0.0, 0.0};
  constexpr Real3 length{1.0e-300, 1.0e154, 1.0e154};
  const DecompositionOptions options{Int3{context.size(), 1, 1}};
  auto decomposition = StructuredDecomposition::create(
      context, kClosureScaleExtent, {false, false, false}, options);
  MeshTopology topology(decomposition);
  const UniformStructuredMesh legacy(kClosureScaleExtent, origin, length,
                                     decomposition);
  MeshGeometry uniform(topology, UniformBoxMapping(origin, length));
  MeshGeometry analytic(
      topology, AnalyticWarpedBoxMapping(origin, length, {0.0, 0.0, 0.0}));

  HUNDUN_CHECK(std::isfinite(legacy.cell_volume_m3()));
  HUNDUN_CHECK(legacy.cell_volume_m3() > 0.0);
  HUNDUN_CHECK(bitwise_equal(*uniform.uniform_spacing_m(), legacy.spacing_m()));
  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    const Int3 global = topology.global_cell(local);
    const Box3 box = topology.owned_global_box();
    const Int3 legacy_local{global.x - box.begin.x, global.y - box.begin.y,
                            global.z - box.begin.z};
    HUNDUN_CHECK(bitwise_equal(uniform.cell_center_m(local),
                               legacy.cell_center(legacy_local)));
    HUNDUN_CHECK(bits(uniform.cell_volume_m3(local)) ==
                 bits(legacy.cell_volume_m3()));

    const auto faces = cell_faces(global);
    const auto x_min =
        topology.find_local_face(topology.global_face_id(faces[0]));
    const auto x_max =
        topology.find_local_face(topology.global_face_id(faces[1]));
    HUNDUN_CHECK(x_min.has_value());
    HUNDUN_CHECK(x_max.has_value());
    const double first = uniform.face_area_m2(*x_min);
    const double second = uniform.face_area_m2(*x_max);
    HUNDUN_CHECK(first > std::numeric_limits<double>::max() / 2.0);
    HUNDUN_CHECK(second > std::numeric_limits<double>::max() / 2.0);
    HUNDUN_CHECK(!std::isfinite(first + second));
    const double analytic_first = analytic.face_area_m2(*x_min);
    const double analytic_second = analytic.face_area_m2(*x_max);
    HUNDUN_CHECK(analytic_first >
                 std::numeric_limits<double>::max() / 2.0);
    HUNDUN_CHECK(analytic_second >
                 std::numeric_limits<double>::max() / 2.0);
    HUNDUN_CHECK(!std::isfinite(analytic_first + analytic_second));
  }

  check_anisotropic_metrics(topology, uniform);
  check_anisotropic_metrics(topology, analytic);
}

void test_jacobian_scale_geometry(const MpiContext& context) {
  constexpr Real3 origin{0.0, 0.0, 0.0};
  constexpr Real3 length{1.0e-210, 1.0e160, 1.0e150};
  const DecompositionOptions options{Int3{1, context.size(), 1}};
  auto decomposition = StructuredDecomposition::create(
      context, kJacobianScaleExtent, {false, false, false}, options);
  MeshTopology topology(decomposition);
  const UniformStructuredMesh legacy(kJacobianScaleExtent, origin, length,
                                     decomposition);
  MeshGeometry uniform(topology, UniformBoxMapping(origin, length));
  MeshGeometry analytic(
      topology, AnalyticWarpedBoxMapping(origin, length, {0.0, 0.0, 0.0}));

  HUNDUN_CHECK(bitwise_equal(*uniform.uniform_spacing_m(), legacy.spacing_m()));
  HUNDUN_CHECK(std::isfinite(legacy.cell_volume_m3()));
  HUNDUN_CHECK(legacy.cell_volume_m3() > 0.0);
  const double expected_jacobian = static_cast<double>(
      static_cast<long double>(length.x) *
      static_cast<long double>(length.y) *
      static_cast<long double>(length.z));
  check_relative(expected_jacobian, 1.0e100, 32.0);

  for (LocalCellId local = 0; local < topology.owned_cell_count(); ++local) {
    const Int3 global = topology.global_cell(local);
    const Box3 box = topology.owned_global_box();
    const Int3 legacy_local{global.x - box.begin.x, global.y - box.begin.y,
                            global.z - box.begin.z};
    HUNDUN_CHECK(bitwise_equal(uniform.cell_center_m(local),
                               legacy.cell_center(legacy_local)));
    HUNDUN_CHECK(bits(uniform.cell_volume_m3(local)) ==
                 bits(legacy.cell_volume_m3()));
    check_relative(uniform.minimum_jacobian_determinant_m3(local),
                   expected_jacobian, 32.0);
    check_relative(analytic.minimum_jacobian_determinant_m3(local),
                   expected_jacobian, 32.0);
    check_relative(analytic.cell_volume_m3(local), legacy.cell_volume_m3(),
                   512.0);
    check_scaled_closure(topology, uniform, local);
    check_scaled_closure(topology, analytic, local);
  }

  std::array<bool, 3> saw_axis{};
  const auto check_faces = [&](const MeshGeometry& geometry) {
    for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
      const auto axis = topology.logical_face(face).axis;
      double expected_area = 0.0;
      switch (axis) {
        case FaceAxis::x:
          expected_area = static_cast<double>(
              static_cast<long double>(legacy.spacing_m().y) *
              static_cast<long double>(legacy.spacing_m().z));
          saw_axis[0] = true;
          break;
        case FaceAxis::y:
          expected_area = static_cast<double>(
              static_cast<long double>(legacy.spacing_m().x) *
              static_cast<long double>(legacy.spacing_m().z));
          saw_axis[1] = true;
          break;
        case FaceAxis::z:
          expected_area = static_cast<double>(
              static_cast<long double>(legacy.spacing_m().x) *
              static_cast<long double>(legacy.spacing_m().y));
          saw_axis[2] = true;
          break;
      }
      check_relative(geometry.face_area_m2(face), expected_area, 128.0);
      const LocalCellId owner = topology.owner(face);
      const Real3 displacement = topology.neighbour(face).has_value()
                                    ? subtract(geometry.cell_center_m(
                                                   *topology.neighbour(face)),
                                               geometry.cell_center_m(owner))
                                    : subtract(geometry.face_center_m(face),
                                               geometry.cell_center_m(owner));
      const double projection =
          dot(geometry.face_area_vector_m2(face, FaceSide::owner),
              displacement);
      HUNDUN_CHECK(std::isfinite(projection));
      HUNDUN_CHECK(projection > 0.0);
      HUNDUN_CHECK(std::isfinite(geometry.face_skewness(face)));
      HUNDUN_CHECK(geometry.face_skewness(face) >= 0.0);
      HUNDUN_CHECK(
          std::isfinite(geometry.face_non_orthogonality_degrees(face)));
      HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(face) >= 0.0);
      HUNDUN_CHECK(geometry.face_non_orthogonality_degrees(face) < 90.0);
    }
  };
  check_faces(uniform);
  check_faces(analytic);
  HUNDUN_CHECK(saw_axis[0]);
  HUNDUN_CHECK(saw_axis[1]);
  HUNDUN_CHECK(saw_axis[2]);
}

void test_invalid_geometry_queries(const MeshTopology& topology,
                                   const MeshGeometry& geometry) {
  expect_error([&] {
    static_cast<void>(geometry.vertex_position_m({kExtent.x + 1, 0, 0}));
  });
  expect_error([&] {
    static_cast<void>(geometry.cell_center_m(topology.local_cell_count()));
  });
  expect_error([&] {
    static_cast<void>(geometry.face_area_m2(topology.local_face_count()));
  });
  expect_error([&] {
    static_cast<void>(
        geometry.face_area_vector_m2(0, static_cast<FaceSide>(200)));
  });
  expect_error([&] {
    static_cast<void>(geometry.cell_closure_m2(topology.local_cell_count()));
  });
  if (topology.ghost_cell_count() > 0) {
    expect_error([&] {
      static_cast<void>(geometry.cell_closure_m2(topology.owned_cell_count()));
    });
  }
}

void run_mpi_tests(const MpiContext& context,
                   std::optional<MeshGeometry>& detached) {
  test_mapping_formulas();
  test_range_safe_mapping_jacobians();
  test_mapping_rejections();

  const DecompositionOptions options{process_grid_for(context.size())};
  auto decomposition =
      StructuredDecomposition::create(context, kExtent, kPeriodic, options);
  MeshTopology topology(decomposition);
  MeshGeometry uniform(topology, UniformBoxMapping(kOrigin, kLength));
  test_uniform_adapter(decomposition, topology, uniform);
  HUNDUN_CHECK(uniform.compatible(topology));
  uniform.require_compatible(topology);

  MeshGeometry warped(topology,
                      AnalyticWarpedBoxMapping(kOrigin, kLength, kAmplitude));
  test_warped_invariants(context, topology, warped);
  test_invalid_geometry_queries(topology, warped);
  HUNDUN_CHECK(warped.compatible(topology));

  {
    auto incompatible_decomposition = StructuredDecomposition::create(
        context, kExtent, {false, false, false}, options);
    MeshTopology incompatible(incompatible_decomposition);
    HUNDUN_CHECK(!warped.compatible(incompatible));
    expect_error([&] { warped.require_compatible(incompatible); });
  }

  expect_error([&] {
    const double tiny = std::numeric_limits<double>::denorm_min();
    MeshGeometry invalid(topology,
                         UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, tiny, 1.0}));
    static_cast<void>(invalid);
  });

  test_extent_one_periodic_seam(context, detached);
  test_mixed_scale_geometry(context);
  test_closure_scale_geometry(context);
  test_nonorthogonality_scale_geometry(context);
  test_jacobian_scale_geometry(context);
  context.barrier();
}

void run_injected_geometry_sample_failure(const MpiContext& context) {
  if (context.rank() == 1) {
    throw std::runtime_error(
        "injected numerical geometry sample failure 2718");
  }
  context.barrier();
}

[[noreturn]] void abort_active_mpi_test(const MpiContext& context) noexcept {
  int initialized = 0;
  int finalized = 1;
  if (MPI_Initialized(&initialized) == MPI_SUCCESS && initialized != 0 &&
      MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0 &&
      context.comm() != MPI_COMM_NULL) {
    (void)MPI_Abort(context.comm(), EXIT_FAILURE);
  }
  std::abort();
}

}  // namespace

int main(int argc, char** argv) {
  const bool inject_sample_failure =
      argc == 2 && std::string(argv[1]) == "injected_sample_failure";
  std::optional<MeshGeometry> detached;
  int active_result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    MpiContext context = MpiContext::duplicate(MPI_COMM_WORLD);
    active_result = hundun::test::run([&] {
      if (inject_sample_failure) {
        run_injected_geometry_sample_failure(context);
      } else {
        run_mpi_tests(context, detached);
      }
    });
    if (active_result != EXIT_SUCCESS) {
      abort_active_mpi_test(context);
    }
  }

  const int finalized_result = hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);
    HUNDUN_CHECK(detached.has_value());
    HUNDUN_CHECK(detached->mapping_kind() == MappingKind::analytic_warped_box);
    HUNDUN_CHECK(detached->cell_volume_m3(0) > 0.0);
    HUNDUN_CHECK(std::isfinite(detached->face_area_m2(0)));
    check_near(detached->vertex_position_m({0, 0, 0}), kOrigin);
    detached.reset();
  });
  return active_result == EXIT_SUCCESS ? finalized_result : active_result;
}
