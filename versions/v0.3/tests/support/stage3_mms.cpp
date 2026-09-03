// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/stage3_mms.hpp"

#include "tests/support/flow_immersed_test_access.hpp"
#include "src/fvm_immersed_boundary_authority_detail.hpp"
#include "tests/support/fvm_immersed_operator_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "tests/support/stage3_dual3.hpp"
#include "tests/support/stage3_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace hundun::test::stage3 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr std::size_t kParallelWorkerBudget = 96U;
using DualVector = std::array<Dual3, 3>;

template <class T> using Vector = std::array<T, 3>;

template <class T> T dot(const Vector<double> &left, const Vector<T> &right) {
  T result{};
  for (std::size_t axis = 0U; axis < 3U; ++axis)
    result = result + left[axis] * right[axis];
  return result;
}

runtime::Real3 normalized(runtime::Real3 value) {
  const double magnitude =
      std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
  return {value.x / magnitude, value.y / magnitude, value.z / magnitude};
}

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

double magnitude(runtime::Real3 value) noexcept {
  return std::sqrt(dot(value, value));
}

std::array<runtime::Int3, 4> logical_face_vertices(mesh::LogicalFace face) {
  const auto c = face.coordinate;
  switch (face.axis) {
  case mesh::FaceAxis::x:
    return {runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x, c.y + 1, c.z},
            runtime::Int3{c.x, c.y + 1, c.z + 1},
            runtime::Int3{c.x, c.y, c.z + 1}};
  case mesh::FaceAxis::y:
    return {runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x, c.y, c.z + 1},
            runtime::Int3{c.x + 1, c.y, c.z + 1},
            runtime::Int3{c.x + 1, c.y, c.z}};
  case mesh::FaceAxis::z:
    return {runtime::Int3{c.x, c.y, c.z}, runtime::Int3{c.x + 1, c.y, c.z},
            runtime::Int3{c.x + 1, c.y + 1, c.z},
            runtime::Int3{c.x, c.y + 1, c.z}};
  }
  throw std::runtime_error("invalid immersed-flow manufactured logical face axis");
}

std::array<Vector<double>, 3> prism_axes() {
  constexpr double degree = kPi / 180.0;
  const double cx = std::cos(17.0 * degree);
  const double sx = std::sin(17.0 * degree);
  const double cy = std::cos(23.0 * degree);
  const double sy = std::sin(23.0 * degree);
  const double cz = std::cos(31.0 * degree);
  const double sz = std::sin(31.0 * degree);
  const std::array<std::array<double, 3>, 3> rotation{
      std::array<double, 3>{cz * cy, cz * sy * sx - sz * cx,
                            cz * sy * cx + sz * sx},
      std::array<double, 3>{sz * cy, sz * sy * sx + cz * cx,
                            sz * sy * cx - cz * sx},
      std::array<double, 3>{-sy, cy * sx, cy * cx}};
  return {Vector<double>{rotation[0][0], rotation[1][0], rotation[2][0]},
          Vector<double>{rotation[0][1], rotation[1][1], rotation[2][1]},
          Vector<double>{rotation[0][2], rotation[1][2], rotation[2][2]}};
}

std::array<runtime::Real3, 3> prism_real_axes() {
  const auto axes = prism_axes();
  return {runtime::Real3{axes[0][0], axes[0][1], axes[0][2]},
          runtime::Real3{axes[1][0], axes[1][1], axes[1][2]},
          runtime::Real3{axes[2][0], axes[2][1], axes[2][2]}};
}

runtime::Real3 orthogonal_unit(runtime::Real3 axis) {
  const runtime::Real3 seed = std::abs(axis.x) < 0.8
                                  ? runtime::Real3{1.0, 0.0, 0.0}
                                  : runtime::Real3{0.0, 1.0, 0.0};
  return normalized(cross(axis, seed));
}

StlFixtureTriangle triangle(runtime::Real3 a, runtime::Real3 b,
                            runtime::Real3 c, runtime::Real3 outward) {
  auto normal = cross(subtract(b, a), subtract(c, a));
  if (dot(normal, outward) < 0.0) {
    std::swap(b, c);
    normal = multiply(-1.0, normal);
  }
  return {normalized(normal), {a, b, c}};
}

double edge_length(runtime::Real3 a, runtime::Real3 b) {
  return magnitude(subtract(a, b));
}

double cylinder_surface_distance(const BodySpec &body, runtime::Real3 point) {
  const auto axis = normalized(body.axis);
  const auto displacement = subtract(point, body.centre_m);
  const double axial = dot(displacement, axis);
  const auto radial_vector = subtract(displacement, multiply(axial, axis));
  const double radial = magnitude(radial_vector);
  const double half = 0.5 * body.length_m;
  double distance = std::numeric_limits<double>::infinity();
  if (std::abs(axial) <= half + 1.0e-13)
    distance = std::abs(radial - body.radius_m);
  if (radial <= body.radius_m + 1.0e-13)
    distance = std::min(distance, std::abs(std::abs(axial) - half));
  return distance;
}

double chord_error(const BodySpec &body, runtime::Real3 point) {
  if (body.kind == BodyKind::sphere ||
      body.kind == BodyKind::inside_sphere_cavity)
    return std::abs(magnitude(subtract(point, body.centre_m)) - body.radius_m);
  if (body.kind == BodyKind::finite_cylinder)
    return cylinder_surface_distance(body, point);
  return 0.0;
}

ManufacturedSurface
finalize_surface(const BodySpec &body,
                 std::vector<StlFixtureTriangle> triangles) {
  ManufacturedSurface result;
  result.triangles = std::move(triangles);
  for (const auto &face : result.triangles) {
    for (std::size_t edge = 0U; edge < 3U; ++edge) {
      const auto a = face.vertices[edge];
      const auto b = face.vertices[(edge + 1U) % 3U];
      result.maximum_edge_m =
          std::max(result.maximum_edge_m, edge_length(a, b));
      result.maximum_chord_error_m =
          std::max(result.maximum_chord_error_m,
                   chord_error(body, multiply(0.5, add(a, b))));
    }
  }
  return result;
}

std::vector<std::pair<double, double>>
gauss_legendre(std::size_t order, double lower, double upper) {
  if (order == 0U)
    throw std::runtime_error("immersed-flow quadrature order must be positive");
  std::vector<std::pair<double, double>> result(order);
  const std::size_t roots = (order + 1U) / 2U;
  for (std::size_t root = 0U; root < roots; ++root) {
    double x = std::cos(kPi * (static_cast<double>(root) + 0.75) /
                        (static_cast<double>(order) + 0.5));
    double derivative = 0.0;
    for (int iteration = 0; iteration < 64; ++iteration) {
      double p0 = 1.0;
      double p1 = x;
      if (order == 1U)
        p1 = x;
      for (std::size_t degree = 2U; degree <= order; ++degree) {
        const double next =
            ((2.0 * static_cast<double>(degree) - 1.0) * x * p1 -
             (static_cast<double>(degree) - 1.0) * p0) /
            static_cast<double>(degree);
        p0 = p1;
        p1 = next;
      }
      const double polynomial = order == 0U ? p0 : p1;
      const double prior = order == 1U ? 1.0 : p0;
      derivative =
          static_cast<double>(order) * (x * polynomial - prior) / (x * x - 1.0);
      const double next = x - polynomial / derivative;
      if (std::abs(next - x) <= 8.0 * std::numeric_limits<double>::epsilon()) {
        x = next;
        break;
      }
      x = next;
    }
    const double weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
    const double midpoint = 0.5 * (lower + upper);
    const double half = 0.5 * (upper - lower);
    result[root] = {midpoint - half * x, half * weight};
    result[order - 1U - root] = {midpoint + half * x, half * weight};
  }
  return result;
}

runtime::Real3 viscous_traction(const MmsSample &sample,
                                runtime::Real3 normal) {
  const double divergence = sample.velocity_gradient_per_s[0][0] +
                            sample.velocity_gradient_per_s[1][1] +
                            sample.velocity_gradient_per_s[2][2];
  const std::array<double, 3> n{normal.x, normal.y, normal.z};
  std::array<double, 3> traction{};
  for (std::size_t row = 0U; row < 3U; ++row)
    for (std::size_t column = 0U; column < 3U; ++column) {
      const double strain = sample.velocity_gradient_per_s[row][column] +
                            sample.velocity_gradient_per_s[column][row] -
                            (row == column ? (2.0 / 3.0) * divergence : 0.0);
      traction[row] += kDynamicViscosityPaS * strain * n[column];
    }
  return {traction[0], traction[1], traction[2]};
}

runtime::Real3
viscous_traction_from_gradient(const std::array<double, 9> &gradient,
                               runtime::Real3 normal) {
  const double divergence = gradient[0] + gradient[4] + gradient[8];
  const std::array<double, 3> n{normal.x, normal.y, normal.z};
  std::array<double, 3> traction{};
  for (std::size_t row = 0U; row < 3U; ++row)
    for (std::size_t column = 0U; column < 3U; ++column) {
      const double strain = gradient[row * 3U + column] +
                            gradient[column * 3U + row] -
                            (row == column ? (2.0 / 3.0) * divergence : 0.0);
      traction[row] += kDynamicViscosityPaS * strain * n[column];
    }
  return {traction[0], traction[1], traction[2]};
}

void accumulate_force(AnalyticForceReference &result, const BodySpec &body,
                      runtime::Real3 point, runtime::Real3 normal,
                      double weight, double time_s) {
  const auto sample = evaluate_mms(body, point, time_s);
  const auto pressure_traction =
      multiply(-sample.mechanical_pressure_pa, normal);
  const auto viscous_traction_value = viscous_traction(sample, normal);
  const auto total_traction = add(pressure_traction, viscous_traction_value);
  const auto pressure = multiply(weight, pressure_traction);
  const auto viscous = multiply(weight, viscous_traction_value);
  result.force.pressure_N = add(result.force.pressure_N, pressure);
  result.force.viscous_N = add(result.force.viscous_N, viscous);
  result.force.total_N = add(result.force.total_N, add(pressure, viscous));
  result.viscous_absolute_component_traction_force_N =
      add(result.viscous_absolute_component_traction_force_N,
          absolute_viscous_traction_force_increment(viscous_traction_value,
                                                    weight));
  result.pressure_traction_rms_force_N +=
      weight * dot(pressure_traction, pressure_traction);
  result.viscous_traction_rms_force_N +=
      weight * dot(viscous_traction_value, viscous_traction_value);
  result.total_traction_rms_force_N +=
      weight * dot(total_traction, total_traction);
  result.surface_area_m2 += weight;
}

AnalyticForceReference finalize_force_reference(AnalyticForceReference result) {
  const auto absolute =
      result.viscous_absolute_component_traction_force_N;
  result.viscous_traction_l1_force_N =
      (absolute.x + absolute.y) + absolute.z;
  result.pressure_traction_rms_force_N =
      std::sqrt(result.surface_area_m2 * result.pressure_traction_rms_force_N);
  result.viscous_traction_rms_force_N =
      std::sqrt(result.surface_area_m2 * result.viscous_traction_rms_force_N);
  result.total_traction_rms_force_N =
      std::sqrt(result.surface_area_m2 * result.total_traction_rms_force_N);
  return result;
}

struct FactorGradient final {
  Dual3 factor;
  DualVector gradient;
};

FactorGradient factor_gradient(const BodySpec &body, const DualVector &point) {
  const DualVector displacement{point[0] - body.centre_m.x,
                                point[1] - body.centre_m.y,
                                point[2] - body.centre_m.z};
  if (body.kind == BodyKind::sphere ||
      body.kind == BodyKind::inside_sphere_cavity) {
    const Dual3 factor =
        displacement[0] * displacement[0] + displacement[1] * displacement[1] +
        displacement[2] * displacement[2] - body.radius_m * body.radius_m;
    return {
        factor,
        {2.0 * displacement[0], 2.0 * displacement[1], 2.0 * displacement[2]}};
  }
  if (body.kind == BodyKind::finite_cylinder) {
    const auto unit = normalized(body.axis);
    const Vector<double> axis{unit.x, unit.y, unit.z};
    const Dual3 axial = dot(axis, displacement);
    const Dual3 radius2 = displacement[0] * displacement[0] +
                          displacement[1] * displacement[1] +
                          displacement[2] * displacement[2] - axial * axial;
    const Dual3 radial = radius2 - body.radius_m * body.radius_m;
    const Dual3 cap = axial * axial - 0.25 * body.length_m * body.length_m;
    DualVector radial_gradient;
    DualVector cap_gradient;
    DualVector gradient;
    for (std::size_t direction = 0U; direction < 3U; ++direction) {
      radial_gradient[direction] =
          2.0 * (displacement[direction] - axial * axis[direction]);
      cap_gradient[direction] = 2.0 * axial * axis[direction];
      gradient[direction] =
          radial_gradient[direction] * cap + radial * cap_gradient[direction];
    }
    return {radial * cap, gradient};
  }

  const auto axes = prism_axes();
  const DualVector local{dot(axes[0], displacement), dot(axes[1], displacement),
                         dot(axes[2], displacement)};
  const std::array<double, 3> half{body.half_lengths_m.x, body.half_lengths_m.y,
                                   body.half_lengths_m.z};
  const DualVector factors{local[0] * local[0] - half[0] * half[0],
                           local[1] * local[1] - half[1] * half[1],
                           local[2] * local[2] - half[2] * half[2]};
  Dual3 factor = factors[0] * factors[1] * factors[2];
  DualVector gradient;
  for (std::size_t direction = 0U; direction < 3U; ++direction) {
    gradient[direction] =
        2.0 * local[0] * axes[0][direction] * factors[1] * factors[2] +
        2.0 * local[1] * axes[1][direction] * factors[0] * factors[2] +
        2.0 * local[2] * axes[2][direction] * factors[0] * factors[1];
    gradient[direction] =
        gradient[direction] * body.prism_mms_factor_multiplier;
  }
  factor = factor * body.prism_mms_factor_multiplier;
  return {factor, gradient};
}

DualVector psi_gradient(const BodySpec &body, const DualVector &point,
                        double amplitude) {
  const Dual3 x = point[0];
  const Dual3 y = point[1];
  const Dual3 z = point[2];
  const Dual3 bx = x * (1.0 - x);
  const Dual3 by = y * (1.0 - y);
  const Dual3 bz = z * (1.0 - z);
  const Dual3 base = bx * by * bz;
  const Dual3 envelope = base * base;
  const DualVector base_gradient{(1.0 - 2.0 * x) * by * bz,
                                 bx * (1.0 - 2.0 * y) * bz,
                                 bx * by * (1.0 - 2.0 * z)};
  DualVector envelope_gradient;
  for (std::size_t direction = 0U; direction < 3U; ++direction)
    envelope_gradient[direction] = 2.0 * base * base_gradient[direction];
  const Dual3 asymmetry = 1.0 + 0.17 * x + 0.11 * y - 0.07 * z;
  constexpr std::array<double, 3> asymmetry_gradient{0.17, 0.11, -0.07};
  const auto factor = factor_gradient(body, point);
  DualVector result;
  for (std::size_t direction = 0U; direction < 3U; ++direction) {
    result[direction] =
        amplitude * (envelope_gradient[direction] * asymmetry * factor.factor *
                         factor.factor +
                     envelope * asymmetry_gradient[direction] * factor.factor *
                         factor.factor +
                     2.0 * envelope * asymmetry * factor.factor *
                         factor.gradient[direction]);
  }
  return result;
}

DualVector velocity(const BodySpec &body, const DualVector &point,
                    double amplitude) {
  const auto gradient = psi_gradient(body, point, amplitude);
  return {gradient[1] + 0.4 * gradient[2], 0.7 * gradient[2] - gradient[0],
          -0.4 * gradient[0] - 0.7 * gradient[1]};
}

Dual3 pressure(const DualVector &point, double time_s) {
  const Dual3 bx = point[0] * (1.0 - point[0]);
  const Dual3 by = point[1] * (1.0 - point[1]);
  const Dual3 bz = point[2] * (1.0 - point[2]);
  const Dual3 boundary_flat_envelope =
      262144.0 * bx * bx * bx * by * by * by * bz * bz * bz;
  const Dual3 spatial =
      0.6 * cos(kPi * point[0]) + 0.4 * cos(kPi * point[1]) -
      0.3 * cos(kPi * point[2]) + 0.15 * cos(2.0 * kPi * point[0]) +
      0.10 * cos(2.0 * kPi * point[1]) - 0.08 * cos(2.0 * kPi * point[2]);
  return boundary_flat_envelope * spatial * std::cos(kPi * time_s);
}

runtime::Real3 as_real3(const DualVector &values) {
  return {values[0].value, values[1].value, values[2].value};
}

runtime::Real3 evaluate_mms_velocity(const BodySpec &body,
                                     runtime::Real3 point_m, double time_s) {
  const DualVector point{Dual3::variable(point_m.x, 0U),
                         Dual3::variable(point_m.y, 1U),
                         Dual3::variable(point_m.z, 2U)};
  const double amplitude = 1.0 + 0.1 * std::sin(2.0 * kPi * time_s);
  return multiply(amplitude, as_real3(velocity(body, point, 1.0)));
}

struct TriangleQuadratureNode final {
  std::array<double, 3> barycentric{};
  double normalized_weight{};
};

constexpr std::array<TriangleQuadratureNode, 3> kQuadraticTriangleRule{
    TriangleQuadratureNode{{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}, 1.0 / 3.0},
    TriangleQuadratureNode{{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}, 1.0 / 3.0},
    TriangleQuadratureNode{{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}, 1.0 / 3.0}};

constexpr std::array<TriangleQuadratureNode, 7> kDegreeFiveTriangleRule{
    TriangleQuadratureNode{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, 0.225},
    TriangleQuadratureNode{{0.059715871789770, 0.470142064105115,
                            0.470142064105115},
                           0.132394152788506},
    TriangleQuadratureNode{{0.470142064105115, 0.059715871789770,
                            0.470142064105115},
                           0.132394152788506},
    TriangleQuadratureNode{{0.470142064105115, 0.470142064105115,
                            0.059715871789770},
                           0.132394152788506},
    TriangleQuadratureNode{{0.797426985353087, 0.101286507323456,
                            0.101286507323456},
                           0.125939180544827},
    TriangleQuadratureNode{{0.101286507323456, 0.797426985353087,
                            0.101286507323456},
                           0.125939180544827},
    TriangleQuadratureNode{{0.101286507323456, 0.101286507323456,
                            0.797426985353087},
                           0.125939180544827}};

template <std::size_t Count, class Evaluator>
runtime::Real3 integrate_triangle_pressure_area(
    runtime::Real3 a, runtime::Real3 b, runtime::Real3 c,
    const std::array<TriangleQuadratureNode, Count> &rule,
    Evaluator &&evaluate) {
  const auto oriented_area =
      multiply(0.5, cross(subtract(b, a), subtract(c, a)));
  double average = 0.0;
  for (const auto &node : rule) {
    const auto point =
        add(multiply(node.barycentric[0], a),
            add(multiply(node.barycentric[1], b),
                multiply(node.barycentric[2], c)));
    average += node.normalized_weight * evaluate(point);
  }
  if (!std::isfinite(oriented_area.x) || !std::isfinite(oriented_area.y) ||
      !std::isfinite(oriented_area.z) || !std::isfinite(average))
    throw std::runtime_error(
        "manufactured pressure face quadrature is non-finite");
  return multiply(average, oriented_area);
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

} // namespace

BodySpec approved_body(BodyKind kind) {
  switch (kind) {
  case BodyKind::sphere:
    return {kind, {0.5, 0.5, 0.5}, 0.18, 0.0, {}, {}};
  case BodyKind::finite_cylinder:
    return {kind, {0.5, 0.5, 0.5}, 0.12, 0.36, normalized({1.0, 1.0, 0.5}), {}};
  case BodyKind::oblique_rectangular_prism:
    return {kind, {0.5, 0.5, 0.5}, 0.0, 0.0, {}, {0.14, 0.11, 0.09}};
  case BodyKind::inside_sphere_cavity:
    return {kind, {0.5, 0.5, 0.5}, 0.32, 0.0, {}, {}};
  }
  throw std::runtime_error("unknown immersed-flow manufactured body");
}

BodySpec translated_sphere() {
  auto result = approved_body(BodyKind::sphere);
  result.centre_m.x += 0.013;
  result.centre_m.y -= 0.009;
  result.centre_m.z += 0.007;
  return result;
}

BodySpec force_certified_oblique_prism() {
  auto result = approved_body(BodyKind::oblique_rectangular_prism);
  const auto half = result.half_lengths_m;
  result.prism_mms_factor_multiplier =
      3.0 * kReferenceLengthM * kReferenceLengthM /
      (half.x * half.y + half.y * half.z + half.z * half.x);
  return result;
}

bool pressure_face_quadrature_oracle_is_mutation_sensitive() noexcept {
  try {
    const runtime::Real3 a{0.0, 0.0, 0.0};
    const runtime::Real3 b{1.0, 0.0, 0.0};
    const runtime::Real3 c{0.0, 1.0, 0.0};
    const auto integrate = [&](const auto &rule, const auto &evaluate) {
      return integrate_triangle_pressure_area(a, b, c, rule, evaluate).z;
    };
    constexpr double tolerance =
        512.0 * std::numeric_limits<double>::epsilon();
    const std::array<double, 6> observed{
        integrate(kQuadraticTriangleRule,
                  [](runtime::Real3) { return 1.0; }),
        integrate(kQuadraticTriangleRule,
                  [](runtime::Real3 point) { return point.x; }),
        integrate(kQuadraticTriangleRule,
                  [](runtime::Real3 point) { return point.y; }),
        integrate(kQuadraticTriangleRule,
                  [](runtime::Real3 point) { return point.x * point.x; }),
        integrate(kQuadraticTriangleRule,
                  [](runtime::Real3 point) { return point.x * point.y; }),
        integrate(kQuadraticTriangleRule,
                  [](runtime::Real3 point) { return point.y * point.y; })};
    constexpr std::array<double, 6> expected{
        1.0 / 2.0, 1.0 / 6.0, 1.0 / 6.0,
        1.0 / 12.0, 1.0 / 24.0, 1.0 / 12.0};
    for (std::size_t index = 0U; index < observed.size(); ++index)
      if (std::abs(observed[index] - expected[index]) > tolerance)
        return false;

    auto mutated = kQuadraticTriangleRule;
    mutated[0].normalized_weight += 0.125;
    const double mutated_constant = integrate(
        mutated, [](runtime::Real3) { return 1.0; });
    return std::abs(mutated_constant - expected[0]) > tolerance;
  } catch (...) {
    return false;
  }
}

MmsSample evaluate_mms(const BodySpec &body, runtime::Real3 point_m,
                       double time_s) {
  const DualVector point{Dual3::variable(point_m.x, 0U),
                         Dual3::variable(point_m.y, 1U),
                         Dual3::variable(point_m.z, 2U)};
  const double omega = 2.0 * kPi;
  const double amplitude = 1.0 + 0.1 * std::sin(omega * time_s);
  const double amplitude_derivative = 0.1 * omega * std::cos(omega * time_s);
  const auto velocity_basis = velocity(body, point, 1.0);
  const auto pressure_value = pressure(point, time_s);

  MmsSample result;
  result.velocity_m_per_s = multiply(amplitude, as_real3(velocity_basis));
  result.mechanical_pressure_pa = pressure_value.value;
  result.mechanical_pressure_gradient_pa_per_m = {pressure_value.gradient[0],
                                                  pressure_value.gradient[1],
                                                  pressure_value.gradient[2]};
  result.mechanical_pressure_hessian_pa_per_m2 = pressure_value.hessian;
  std::array<double *, 3> source{&result.body_source_N_per_m3.x,
                                 &result.body_source_N_per_m3.y,
                                 &result.body_source_N_per_m3.z};
  for (std::size_t component = 0U; component < 3U; ++component) {
    double convection = 0.0;
    double laplacian = 0.0;
    for (std::size_t direction = 0U; direction < 3U; ++direction) {
      result.velocity_gradient_per_s[component][direction] =
          amplitude * velocity_basis[component].gradient[direction];
      convection += velocity_basis[direction].value *
                    velocity_basis[component].gradient[direction];
      laplacian += velocity_basis[component].hessian[direction][direction];
    }
    *source[component] =
        amplitude_derivative * velocity_basis[component].value +
        amplitude * amplitude * convection +
        pressure_value.gradient[component] -
        kDynamicViscosityPaS * amplitude * laplacian;
  }
  return result;
}

double implicit_factor(const BodySpec &body, runtime::Real3 point_m) {
  const DualVector point{Dual3(point_m.x), Dual3(point_m.y), Dual3(point_m.z)};
  return factor_gradient(body, point).factor.value;
}

ManufacturedSurface make_manufactured_surface(const BodySpec &body,
                                              double h_max_m) {
  if (!(h_max_m > 0.0) || !std::isfinite(h_max_m))
    throw std::runtime_error("immersed-flow manufactured mesh spacing is invalid");
  const double target = 0.45 * h_max_m;
  std::vector<StlFixtureTriangle> triangles;

  if (body.kind == BodyKind::sphere ||
      body.kind == BodyKind::inside_sphere_cavity) {
    constexpr double golden = 1.6180339887498948482045868343656381;
    std::array<runtime::Real3, 12> vertices{
        runtime::Real3{-1.0, golden, 0.0},  runtime::Real3{1.0, golden, 0.0},
        runtime::Real3{-1.0, -golden, 0.0}, runtime::Real3{1.0, -golden, 0.0},
        runtime::Real3{0.0, -1.0, golden},  runtime::Real3{0.0, 1.0, golden},
        runtime::Real3{0.0, -1.0, -golden}, runtime::Real3{0.0, 1.0, -golden},
        runtime::Real3{golden, 0.0, -1.0},  runtime::Real3{golden, 0.0, 1.0},
        runtime::Real3{-golden, 0.0, -1.0}, runtime::Real3{-golden, 0.0, 1.0}};
    for (auto &value : vertices)
      value = normalized(value);
    constexpr std::array<std::array<std::size_t, 3>, 20> faces{
        std::array<std::size_t, 3>{0, 11, 5},
        {0, 5, 1},
        {0, 1, 7},
        {0, 7, 10},
        {0, 10, 11},
        {1, 5, 9},
        {5, 11, 4},
        {11, 10, 2},
        {10, 7, 6},
        {7, 1, 8},
        {3, 9, 4},
        {3, 4, 2},
        {3, 2, 6},
        {3, 6, 8},
        {3, 8, 9},
        {4, 9, 5},
        {2, 4, 11},
        {6, 2, 10},
        {8, 6, 7},
        {9, 8, 1}};
    using UnitTriangle = std::array<runtime::Real3, 3>;
    std::vector<UnitTriangle> unit_faces;
    unit_faces.reserve(faces.size());
    for (const auto &face : faces)
      unit_faces.push_back(
          {vertices[face[0]], vertices[face[1]], vertices[face[2]]});
    const auto maximum_edge = [&] {
      double value = 0.0;
      for (const auto &face : unit_faces)
        for (std::size_t edge = 0U; edge < 3U; ++edge)
          value = std::max(value,
                           body.radius_m *
                               edge_length(face[edge], face[(edge + 1U) % 3U]));
      return value;
    };
    while (maximum_edge() > target) {
      std::vector<UnitTriangle> refined;
      refined.reserve(4U * unit_faces.size());
      for (const auto &face : unit_faces) {
        const auto ab = normalized(add(face[0], face[1]));
        const auto bc = normalized(add(face[1], face[2]));
        const auto ca = normalized(add(face[2], face[0]));
        refined.push_back({face[0], ab, ca});
        refined.push_back({ab, face[1], bc});
        refined.push_back({ca, bc, face[2]});
        refined.push_back({ab, bc, ca});
      }
      unit_faces = std::move(refined);
    }
    triangles.reserve(unit_faces.size());
    for (const auto &face : unit_faces) {
      const auto a = add(body.centre_m, multiply(body.radius_m, face[0]));
      const auto b = add(body.centre_m, multiply(body.radius_m, face[1]));
      const auto c = add(body.centre_m, multiply(body.radius_m, face[2]));
      triangles.push_back(
          triangle(a, b, c, add(add(face[0], face[1]), face[2])));
    }
    return finalize_surface(body, std::move(triangles));
  }

  if (body.kind == BodyKind::finite_cylinder) {
    const auto axis = normalized(body.axis);
    const auto first = orthogonal_unit(axis);
    const auto second = cross(axis, first);
    const std::size_t angular_count = std::max<std::size_t>(
        8U, static_cast<std::size_t>(std::ceil(std::sqrt(2.0) * 2.0 * kPi *
                                               body.radius_m / target)));
    const std::size_t axial_count = std::max<std::size_t>(
        1U, static_cast<std::size_t>(
                std::ceil(std::sqrt(2.0) * body.length_m / target)));
    const std::size_t radial_count = std::max<std::size_t>(
        1U, static_cast<std::size_t>(
                std::ceil(std::sqrt(2.0) * body.radius_m / target)));
    const auto radial = [&](std::size_t angle) {
      const double phi = 2.0 * kPi *
                         static_cast<double>(angle % angular_count) /
                         static_cast<double>(angular_count);
      return add(multiply(std::cos(phi), first),
                 multiply(std::sin(phi), second));
    };
    const auto side_point = [&](std::size_t axial_index, std::size_t angle) {
      const double axial = -0.5 * body.length_m +
                           body.length_m * static_cast<double>(axial_index) /
                               static_cast<double>(axial_count);
      return add(body.centre_m, add(multiply(axial, axis),
                                    multiply(body.radius_m, radial(angle))));
    };
    for (std::size_t axial_index = 0U; axial_index < axial_count; ++axial_index)
      for (std::size_t angle = 0U; angle < angular_count; ++angle) {
        const auto next = (angle + 1U) % angular_count;
        const auto p00 = side_point(axial_index, angle);
        const auto p01 = side_point(axial_index, next);
        const auto p10 = side_point(axial_index + 1U, angle);
        const auto p11 = side_point(axial_index + 1U, next);
        const auto outward = add(radial(angle), radial(next));
        triangles.push_back(triangle(p00, p10, p11, outward));
        triangles.push_back(triangle(p00, p11, p01, outward));
      }
    for (const double sign : {-1.0, 1.0}) {
      const auto cap_centre =
          add(body.centre_m, multiply(sign * 0.5 * body.length_m, axis));
      const auto outward = multiply(sign, axis);
      const auto ring_point = [&](std::size_t ring, std::size_t angle) {
        const double radius = body.radius_m * static_cast<double>(ring) /
                              static_cast<double>(radial_count);
        return add(cap_centre, multiply(radius, radial(angle)));
      };
      for (std::size_t angle = 0U; angle < angular_count; ++angle) {
        const auto next = (angle + 1U) % angular_count;
        triangles.push_back(triangle(cap_centre, ring_point(1U, angle),
                                     ring_point(1U, next), outward));
      }
      for (std::size_t ring = 1U; ring < radial_count; ++ring)
        for (std::size_t angle = 0U; angle < angular_count; ++angle) {
          const auto next = (angle + 1U) % angular_count;
          const auto p00 = ring_point(ring, angle);
          const auto p01 = ring_point(ring, next);
          const auto p10 = ring_point(ring + 1U, angle);
          const auto p11 = ring_point(ring + 1U, next);
          triangles.push_back(triangle(p00, p10, p11, outward));
          triangles.push_back(triangle(p00, p11, p01, outward));
        }
    }
    return finalize_surface(body, std::move(triangles));
  }

  const auto axes = prism_real_axes();
  const std::array<double, 3> half{body.half_lengths_m.x, body.half_lengths_m.y,
                                   body.half_lengths_m.z};
  for (std::size_t fixed = 0U; fixed < 3U; ++fixed) {
    const std::size_t first = (fixed + 1U) % 3U;
    const std::size_t second = (fixed + 2U) % 3U;
    const std::size_t first_count = std::max<std::size_t>(
        1U, static_cast<std::size_t>(
                std::ceil(std::sqrt(2.0) * 2.0 * half[first] / target)));
    const std::size_t second_count = std::max<std::size_t>(
        1U, static_cast<std::size_t>(
                std::ceil(std::sqrt(2.0) * 2.0 * half[second] / target)));
    for (const double sign : {-1.0, 1.0}) {
      const auto outward = multiply(sign, axes[fixed]);
      const auto point = [&](std::size_t first_index,
                             std::size_t second_index) {
        const double first_coordinate =
            -half[first] + 2.0 * half[first] *
                               static_cast<double>(first_index) /
                               static_cast<double>(first_count);
        const double second_coordinate =
            -half[second] + 2.0 * half[second] *
                                static_cast<double>(second_index) /
                                static_cast<double>(second_count);
        return add(body.centre_m,
                   add(multiply(sign * half[fixed], axes[fixed]),
                       add(multiply(first_coordinate, axes[first]),
                           multiply(second_coordinate, axes[second]))));
      };
      for (std::size_t i = 0U; i < first_count; ++i)
        for (std::size_t j = 0U; j < second_count; ++j) {
          const auto p00 = point(i, j);
          const auto p10 = point(i + 1U, j);
          const auto p01 = point(i, j + 1U);
          const auto p11 = point(i + 1U, j + 1U);
          triangles.push_back(triangle(p00, p10, p11, outward));
          triangles.push_back(triangle(p00, p11, p01, outward));
        }
    }
  }
  return finalize_surface(body, std::move(triangles));
}

AnalyticForceReference analytic_force_reference(const BodySpec &body,
                                                double time_s,
                                                std::size_t quadrature_order) {
  const auto unit = gauss_legendre(quadrature_order, -1.0, 1.0);
  const auto angle = gauss_legendre(quadrature_order, 0.0, 2.0 * kPi);
  AnalyticForceReference result;
  if (body.kind == BodyKind::sphere ||
      body.kind == BodyKind::inside_sphere_cavity) {
    for (const auto &[mu, mu_weight] : unit)
      for (const auto &[phi, phi_weight] : angle) {
        const double radial = std::sqrt(std::max(0.0, 1.0 - mu * mu));
        const runtime::Real3 geometric_normal{radial * std::cos(phi),
                                              radial * std::sin(phi), mu};
        const runtime::Real3 fluid_normal =
            body.kind == BodyKind::inside_sphere_cavity
                ? multiply(-1.0, geometric_normal)
                : geometric_normal;
        accumulate_force(
            result, body,
            add(body.centre_m, multiply(body.radius_m, geometric_normal)),
            fluid_normal,
            body.radius_m * body.radius_m * mu_weight * phi_weight, time_s);
      }
    return finalize_force_reference(result);
  }
  if (body.kind == BodyKind::finite_cylinder) {
    const auto axis = normalized(body.axis);
    const auto first = orthogonal_unit(axis);
    const auto second = cross(axis, first);
    const auto axial = gauss_legendre(quadrature_order, -0.5 * body.length_m,
                                      0.5 * body.length_m);
    for (const auto &[z, z_weight] : axial)
      for (const auto &[phi, phi_weight] : angle) {
        const auto normal = add(multiply(std::cos(phi), first),
                                multiply(std::sin(phi), second));
        const auto point =
            add(body.centre_m,
                add(multiply(z, axis), multiply(body.radius_m, normal)));
        accumulate_force(result, body, point, normal,
                         body.radius_m * z_weight * phi_weight, time_s);
      }
    const auto radial = gauss_legendre(quadrature_order, 0.0, body.radius_m);
    for (const double sign : {-1.0, 1.0})
      for (const auto &[radius, radius_weight] : radial)
        for (const auto &[phi, phi_weight] : angle) {
          const auto direction = add(multiply(std::cos(phi), first),
                                     multiply(std::sin(phi), second));
          const auto point =
              add(body.centre_m, add(multiply(sign * 0.5 * body.length_m, axis),
                                     multiply(radius, direction)));
          accumulate_force(result, body, point, multiply(sign, axis),
                           radius * radius_weight * phi_weight, time_s);
        }
    return finalize_force_reference(result);
  }

  const auto axes = prism_real_axes();
  const std::array<double, 3> half{body.half_lengths_m.x, body.half_lengths_m.y,
                                   body.half_lengths_m.z};
  for (std::size_t fixed = 0U; fixed < 3U; ++fixed) {
    const std::size_t first = (fixed + 1U) % 3U;
    const std::size_t second = (fixed + 2U) % 3U;
    const auto first_nodes =
        gauss_legendre(quadrature_order, -half[first], half[first]);
    const auto second_nodes =
        gauss_legendre(quadrature_order, -half[second], half[second]);
    for (const double sign : {-1.0, 1.0})
      for (const auto &[u, u_weight] : first_nodes)
        for (const auto &[v, v_weight] : second_nodes) {
          const auto point = add(
              body.centre_m,
              add(multiply(sign * half[fixed], axes[fixed]),
                  add(multiply(u, axes[first]), multiply(v, axes[second]))));
          accumulate_force(result, body, point, multiply(sign, axes[fixed]),
                           u_weight * v_weight, time_s);
        }
  }
  return finalize_force_reference(result);
}

namespace detail {

MmsFaceHistory evaluate_face_history_validated(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    mesh::LocalFaceId face, const BodySpec &body, double time_s) {
  const auto coordinates = logical_face_vertices(topology.logical_face(face));
  std::array<runtime::Real3, 4> vertices{};
  for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex)
    vertices[vertex] = geometry.vertex_position_m(coordinates[vertex]);

  const auto owner_area =
      geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
  const double owner_area_squared = dot(owner_area, owner_area);
  if (!(owner_area_squared > 0.0) || !std::isfinite(owner_area_squared))
    throw std::runtime_error(
        "immersed-flow face-history owner area is non-finite or degenerate");

  static const auto quadrature = gauss_legendre(5U, 0.0, 1.0);
  double positive_axis_flux = 0.0;
  runtime::Real3 integrated_positive_area{};
  const auto accumulate = [&](runtime::Real3 point,
                              runtime::Real3 area_weight) {
    const auto velocity = evaluate_mms_velocity(body, point, time_s);
    if (!finite(velocity))
      throw std::runtime_error(
          "immersed-flow face-history quadrature sample is non-finite");
    positive_axis_flux += dot(velocity, area_weight);
    integrated_positive_area = add(integrated_positive_area, area_weight);
  };

  if (geometry.mapping_kind() == mesh::MappingKind::uniform_box) {
    const auto edge_a = subtract(vertices[1], vertices[0]);
    const auto edge_b = subtract(vertices[3], vertices[0]);
    const auto positive_area = cross(edge_a, edge_b);
    for (const auto &[a, weight_a] : quadrature)
      for (const auto &[b, weight_b] : quadrature) {
        const auto point =
            add(vertices[0], add(multiply(a, edge_a), multiply(b, edge_b)));
        accumulate(point, multiply(weight_a * weight_b, positive_area));
      }
  } else {
    constexpr std::array<std::array<std::size_t, 3>, 2> triangles{
        std::array<std::size_t, 3>{0U, 1U, 2U},
        std::array<std::size_t, 3>{0U, 2U, 3U}};
    for (const auto &triangle : triangles) {
      const auto origin = vertices[triangle[0]];
      const auto edge_a = subtract(vertices[triangle[1]], origin);
      const auto edge_b = subtract(vertices[triangle[2]], origin);
      const auto positive_area = cross(edge_a, edge_b);
      for (const auto &[a, weight_a] : quadrature)
        for (const auto &[b, weight_b] : quadrature) {
          const double one_minus_a = 1.0 - a;
          const auto point =
              add(origin,
                  add(multiply(a, edge_a), multiply(one_minus_a * b, edge_b)));
          accumulate(point, multiply(weight_a * weight_b * one_minus_a,
                                     positive_area));
        }
    }
  }

  const double orientation = dot(integrated_positive_area, owner_area);
  if (!std::isfinite(positive_axis_flux) || !finite(integrated_positive_area) ||
      orientation == 0.0 || !std::isfinite(orientation))
    throw std::runtime_error("immersed-flow face-history quadrature is invalid");
  const double area_tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
                                std::max(1.0, std::sqrt(owner_area_squared));
  const auto oriented_integrated_area =
      orientation > 0.0 ? integrated_positive_area
                        : multiply(-1.0, integrated_positive_area);
  if (max_abs_difference(oriented_integrated_area, owner_area) > area_tolerance)
    throw std::runtime_error(
        "immersed-flow face-history quadrature area disagrees with geometry");

  const double owner_flux =
      orientation > 0.0 ? positive_axis_flux : -positive_axis_flux;
  const auto center_velocity =
      evaluate_mms_velocity(body, geometry.face_center_m(face), time_s);
  if (!finite(center_velocity))
    throw std::runtime_error(
        "immersed-flow face-history center sample is non-finite");
  const double correction =
      (owner_flux - dot(center_velocity, owner_area)) / owner_area_squared;
  const auto face_velocity =
      add(center_velocity, multiply(correction, owner_area));
  if (!finite(face_velocity) || !std::isfinite(owner_flux))
    throw std::runtime_error("immersed-flow face-history result is non-finite");
  return {face_velocity, owner_flux};
}

MmsCellAverage evaluate_cell_average_at_global(
    const mesh::MeshGeometry &geometry, runtime::Int3 global,
    const BodySpec &body, double time_s,
    std::optional<double> expected_volume_m3 = std::nullopt) {
  constexpr std::array<runtime::Int3, 8> offsets{
      runtime::Int3{0, 0, 0}, runtime::Int3{1, 0, 0}, runtime::Int3{1, 1, 0},
      runtime::Int3{0, 1, 0}, runtime::Int3{0, 0, 1}, runtime::Int3{1, 0, 1},
      runtime::Int3{1, 1, 1}, runtime::Int3{0, 1, 1}};
  std::array<runtime::Real3, 8> vertices{};
  runtime::Real3 reference{};
  for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex) {
    vertices[vertex] = geometry.vertex_position_m(
        {global.x + offsets[vertex].x, global.y + offsets[vertex].y,
         global.z + offsets[vertex].z});
    reference = add(reference, multiply(0.125, vertices[vertex]));
  }
  constexpr std::array<std::array<std::size_t, 3>, 12> triangles{{
      {{0, 4, 7}},
      {{0, 7, 3}},
      {{1, 2, 6}},
      {{1, 6, 5}},
      {{0, 1, 5}},
      {{0, 5, 4}},
      {{3, 7, 6}},
      {{3, 6, 2}},
      {{0, 3, 2}},
      {{0, 2, 1}},
      {{4, 5, 6}},
      {{4, 6, 7}},
  }};
  constexpr double root = 0.774596669241483377035853079956;
  constexpr std::array<double, 3> nodes{-root, 0.0, root};
  constexpr std::array<double, 3> weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  runtime::Real3 source_integral{};
  runtime::Real3 velocity_integral{};
  runtime::Real3 pressure_gradient_integral{};
  double pressure_integral = 0.0;
  double volume = 0.0;
  const auto accumulate_sample = [&](runtime::Real3 point, double weight) {
    const auto sample = evaluate_mms(body, point, time_s);
    source_integral =
        add(source_integral, multiply(weight, sample.body_source_N_per_m3));
    velocity_integral =
        add(velocity_integral, multiply(weight, sample.velocity_m_per_s));
    pressure_gradient_integral = add(
        pressure_gradient_integral,
        multiply(weight, sample.mechanical_pressure_gradient_pa_per_m));
    pressure_integral += weight * sample.mechanical_pressure_pa;
    volume += weight;
  };
  if (geometry.mapping_kind() == mesh::MappingKind::uniform_box) {
    for (std::size_t kz = 0U; kz < nodes.size(); ++kz)
      for (std::size_t jy = 0U; jy < nodes.size(); ++jy)
        for (std::size_t ix = 0U; ix < nodes.size(); ++ix) {
          runtime::Real3 point{};
          runtime::Real3 d_xi{};
          runtime::Real3 d_eta{};
          runtime::Real3 d_zeta{};
          for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex) {
            const double sx = offsets[vertex].x == 0 ? -1.0 : 1.0;
            const double sy = offsets[vertex].y == 0 ? -1.0 : 1.0;
            const double sz = offsets[vertex].z == 0 ? -1.0 : 1.0;
            const double nx = 1.0 + sx * nodes[ix];
            const double ny = 1.0 + sy * nodes[jy];
            const double nz = 1.0 + sz * nodes[kz];
            point =
                add(point, multiply(0.125 * nx * ny * nz, vertices[vertex]));
            d_xi = add(d_xi, multiply(0.125 * sx * ny * nz, vertices[vertex]));
            d_eta =
                add(d_eta, multiply(0.125 * nx * sy * nz, vertices[vertex]));
            d_zeta =
                add(d_zeta, multiply(0.125 * nx * ny * sz, vertices[vertex]));
          }
          const double jacobian = dot(d_xi, cross(d_eta, d_zeta));
          if (!(jacobian > 0.0) || !std::isfinite(jacobian))
            throw std::runtime_error(
                "immersed-flow source quadrature Jacobian is invalid");
          accumulate_sample(point,
                            weights[ix] * weights[jy] * weights[kz] * jacobian);
        }
  } else {
    for (const auto triangle : triangles) {
      const auto edge1 = subtract(vertices[triangle[0]], reference);
      const auto edge2 = subtract(vertices[triangle[1]], reference);
      const auto edge3 = subtract(vertices[triangle[2]], reference);
      const double determinant = dot(edge1, cross(edge2, edge3));
      if (!(determinant > 0.0) || !std::isfinite(determinant))
        throw std::runtime_error(
            "immersed-flow source quadrature tetrahedron is invalid");
      for (std::size_t ia = 0U; ia < nodes.size(); ++ia)
        for (std::size_t ib = 0U; ib < nodes.size(); ++ib)
          for (std::size_t ic = 0U; ic < nodes.size(); ++ic) {
            const double a = 0.5 * (1.0 + nodes[ia]);
            const double b = 0.5 * (1.0 + nodes[ib]);
            const double c = 0.5 * (1.0 + nodes[ic]);
            const double one_minus_a = 1.0 - a;
            const double one_minus_b = 1.0 - b;
            const auto point =
                add(reference,
                    add(multiply(a, edge1),
                        add(multiply(one_minus_a * b, edge2),
                            multiply(one_minus_a * one_minus_b * c, edge3))));
            const double weight = 0.125 * weights[ia] * weights[ib] *
                                  weights[ic] * determinant * one_minus_a *
                                  one_minus_a * one_minus_b;
            accumulate_sample(point, weight);
          }
    }
  }
  if (!(volume > 0.0) || !std::isfinite(volume))
    throw std::runtime_error("immersed-flow source quadrature volume is invalid");
  if (expected_volume_m3.has_value()) {
    const double geometry_volume = *expected_volume_m3;
    const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, geometry_volume);
    if (std::abs(volume - geometry_volume) > tolerance)
      throw std::runtime_error(
          "immersed-flow source quadrature volume disagrees with geometry");
  }
  return {multiply(1.0 / volume, velocity_integral), pressure_integral / volume,
          multiply(1.0 / volume, source_integral),
          multiply(1.0 / volume, pressure_gradient_integral)};
}

MmsCellAverage evaluate_cell_average_validated(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    mesh::LocalCellId cell, const BodySpec &body, double time_s) {
  geometry.require_compatible(topology);
  return evaluate_cell_average_at_global(
      geometry, topology.global_cell(cell), body, time_s,
      geometry.cell_volume_m3(cell));
}

double cell_volume_global(const mesh::MeshGeometry &geometry,
                          runtime::Int3 global) {
  constexpr std::array<runtime::Int3, 8> offsets{
      runtime::Int3{0, 0, 0}, runtime::Int3{1, 0, 0}, runtime::Int3{1, 1, 0},
      runtime::Int3{0, 1, 0}, runtime::Int3{0, 0, 1}, runtime::Int3{1, 0, 1},
      runtime::Int3{1, 1, 1}, runtime::Int3{0, 1, 1}};
  constexpr std::array<std::array<std::size_t, 3>, 12> triangles{{
      {{0, 4, 7}},
      {{0, 7, 3}},
      {{1, 2, 6}},
      {{1, 6, 5}},
      {{0, 1, 5}},
      {{0, 5, 4}},
      {{3, 7, 6}},
      {{3, 6, 2}},
      {{0, 3, 2}},
      {{0, 2, 1}},
      {{4, 5, 6}},
      {{4, 6, 7}},
  }};
  std::array<runtime::Real3, 8> vertices{};
  runtime::Real3 reference{};
  for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex) {
    vertices[vertex] = geometry.vertex_position_m(
        {global.x + offsets[vertex].x, global.y + offsets[vertex].y,
         global.z + offsets[vertex].z});
    reference = add(reference, multiply(0.125, vertices[vertex]));
  }
  double volume6 = 0.0;
  for (const auto &triangle : triangles) {
    const auto edge1 = subtract(vertices[triangle[0]], reference);
    const auto edge2 = subtract(vertices[triangle[1]], reference);
    const auto edge3 = subtract(vertices[triangle[2]], reference);
    volume6 += dot(edge1, cross(edge2, edge3));
  }
  const double volume = volume6 / 6.0;
  if (!(volume > 0.0) || !std::isfinite(volume))
    throw std::runtime_error("immersed-flow global cell volume is invalid");
  return volume;
}

} // namespace detail

MmsCellAverage evaluate_cell_average(const mesh::MeshTopology &topology,
                                     const mesh::MeshGeometry &geometry,
                                     mesh::LocalCellId cell,
                                     const BodySpec &body, double time_s) {
  geometry.require_compatible(topology);
  return detail::evaluate_cell_average_validated(topology, geometry, cell, body,
                                                 time_s);
}

MmsFaceHistory evaluate_face_history(const mesh::MeshTopology &topology,
                                     const mesh::MeshGeometry &geometry,
                                     mesh::LocalFaceId face,
                                     const BodySpec &body, double time_s) {
  geometry.require_compatible(topology);
  return detail::evaluate_face_history_validated(topology, geometry, face, body,
                                                 time_s);
}

runtime::Real3 source_cell_average(const mesh::MeshTopology &topology,
                                   const mesh::MeshGeometry &geometry,
                                   mesh::LocalCellId cell, const BodySpec &body,
                                   double time_s) {
  return evaluate_cell_average(topology, geometry, cell, body, time_s)
      .body_source_N_per_m3;
}

bool finite(const MmsSample &sample) noexcept {
  if (!finite(sample.velocity_m_per_s) ||
      !finite(sample.body_source_N_per_m3) ||
      !finite(sample.mechanical_pressure_gradient_pa_per_m) ||
      !std::isfinite(sample.mechanical_pressure_pa))
    return false;
  for (const auto &row : sample.velocity_gradient_per_s)
    for (const double value : row)
      if (!std::isfinite(value))
        return false;
  for (const auto &row : sample.mechanical_pressure_hessian_pa_per_m2)
    for (const double value : row)
      if (!std::isfinite(value))
        return false;
  return true;
}

double max_abs(runtime::Real3 value) noexcept {
  return std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
}

double max_abs_difference(runtime::Real3 left, runtime::Real3 right) noexcept {
  return max_abs(subtract(left, right));
}

runtime::Real3 absolute_viscous_traction_force_increment(
    runtime::Real3 traction, double positive_weight) noexcept {
  return {positive_weight * std::abs(traction.x),
          positive_weight * std::abs(traction.y),
          positive_weight * std::abs(traction.z)};
}

bool signed_viscous_force_reference_preflight_accepts(
    const AnalyticForceReference &coarse,
    const AnalyticForceReference &refined) noexcept {
  constexpr double epsilon = std::numeric_limits<double>::epsilon();
  constexpr double rms_discrimination = 256.0 * epsilon;
  constexpr double component_and_net_floor = 4096.0 * epsilon;
  const double cancellation_floor = 1024.0 * std::sqrt(epsilon);
  const auto positive_scale = [](const AnalyticForceReference &reference,
                                 double &scale) {
    if (!(reference.surface_area_m2 > 0.0) ||
        !std::isfinite(reference.surface_area_m2))
      return false;
    scale = kDynamicViscosityPaS * kReferenceVelocityMPerS *
            reference.surface_area_m2 / kReferenceLengthM;
    return scale > 0.0 && std::isfinite(scale);
  };
  double coarse_scale = 0.0;
  double refined_scale = 0.0;
  if (!positive_scale(coarse, coarse_scale) ||
      !positive_scale(refined, refined_scale))
    return false;
  const auto stable_vector = [](runtime::Real3 left, runtime::Real3 right,
                                double scale) {
    return max_abs_difference(left, right) <=
           1.0e-13 * std::max(1.0, scale);
  };
  const auto stable_scalar = [](double left, double right, double scale) {
    return std::abs(left - right) <= 1.0e-13 * std::max(1.0, scale);
  };
  if (!(refined.viscous_traction_rms_force_N >
        rms_discrimination * refined_scale) ||
      !stable_scalar(coarse.viscous_traction_rms_force_N,
                     refined.viscous_traction_rms_force_N, refined_scale))
    return false;

  const auto level_accepts = [&](const AnalyticForceReference &reference,
                                 double scale,
                                 double &traction_certificate) {
    const auto force = reference.force.viscous_N;
    if (!std::isfinite(force.x) || !std::isfinite(force.y) ||
        !std::isfinite(force.z))
      return false;
    const auto absolute =
        reference.viscous_absolute_component_traction_force_N;
    const double component_floor = component_and_net_floor * scale;
    if (!(absolute.x >= component_floor) ||
        !(absolute.y >= component_floor) ||
        !(absolute.z >= component_floor))
      return false;
    const double fixed_l1 = (absolute.x + absolute.y) + absolute.z;
    if (!(reference.viscous_traction_l1_force_N == fixed_l1) ||
        !std::isfinite(reference.viscous_traction_l1_force_N))
      return false;
    const double force_net = max_abs(force);
    if (!(force_net >= component_and_net_floor * scale))
      return false;
    traction_certificate =
        std::sqrt(3.0) * reference.viscous_traction_rms_force_N;
    return force_net >= cancellation_floor * traction_certificate;
  };
  double coarse_traction_certificate = 0.0;
  double refined_traction_certificate = 0.0;
  if (!level_accepts(coarse, coarse_scale, coarse_traction_certificate) ||
      !level_accepts(refined, refined_scale, refined_traction_certificate))
    return false;
  return stable_vector(coarse.force.viscous_N, refined.force.viscous_N,
                       refined_scale) &&
         stable_scalar(coarse_traction_certificate,
                       refined_traction_certificate, refined_scale);
}

namespace {

config::FlowCaseConfig
manufactured_flow_config(const ManufacturedCase &definition, int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-task11-manufactured";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = definition.process_grid;
  result.mesh.cells = {definition.cells, definition.cells, definition.cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = kReferenceDensityKgPerM3;
  result.physics.dynamic_viscosity_pa_s = kDynamicViscosityPaS;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    result.boundaries[patch].patch = names[patch];
    result.boundaries[patch].type = config::BoundaryType::no_slip_wall;
  }
  return result;
}

runtime::FieldDescriptor
manufactured_cell_descriptor(const char *name, std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor
manufactured_face_descriptor(const char *name, std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

double canonical_manufactured_zero(double value) {
  return std::abs(value) <= 1.0e-30 ? 0.0 : value;
}

MmsCellAverage rescale_average(const MmsCellAverage &reference,
                               double reference_time_s, double time_s) {
  const double reference_amplitude =
      1.0 + 0.1 * std::sin(2.0 * kPi * reference_time_s);
  const double amplitude = 1.0 + 0.1 * std::sin(2.0 * kPi * time_s);
  const double pressure_reference = std::cos(kPi * reference_time_s);
  if (!(std::abs(reference_amplitude) > 0.0) ||
      !(std::abs(pressure_reference) > 0.0))
    throw std::runtime_error("invalid manufactured history scale");
  const double velocity_scale = amplitude / reference_amplitude;
  const double pressure_scale = std::cos(kPi * time_s) / pressure_reference;
  return {multiply(velocity_scale, reference.velocity_m_per_s),
          pressure_scale * reference.mechanical_pressure_pa,
          {},
          multiply(pressure_scale,
                   reference.mechanical_pressure_gradient_pa_per_m)};
}

void check_mpi(int code, const char *operation) {
  if (code != MPI_SUCCESS)
    throw std::runtime_error(std::string(operation) + " failed");
}

void sum_in_place(const runtime::MpiContext &mpi, double *values,
                  std::size_t count) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error("manufactured reduction count exceeds MPI int");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured sum)");
}

void max_in_place(const runtime::MpiContext &mpi, double *values,
                  std::size_t count) {
  if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error("manufactured reduction count exceeds MPI int");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, MPI_MAX, mpi.comm()),
            "MPI_Allreduce(manufactured maximum)");
}

double vector_error(runtime::Real3 actual, runtime::Real3 expected,
                    double scale) {
  if (!(scale > 0.0) || !std::isfinite(scale))
    throw std::runtime_error("manufactured force scale is invalid");
  return max_abs_difference(actual, expected) / scale;
}

AnalyticForceReference verified_force_reference(const BodySpec &body,
                                                double time_s) {
  const auto coarse = analytic_force_reference(body, time_s, 48U);
  const auto refined = analytic_force_reference(body, time_s, 96U);
  const double pressure_scale =
      kReferenceDensityKgPerM3 * kReferenceVelocityMPerS *
      kReferenceVelocityMPerS * refined.surface_area_m2;
  const double viscous_scale = kDynamicViscosityPaS * kReferenceVelocityMPerS *
                               refined.surface_area_m2 / kReferenceLengthM;
  const double total_scale = std::max(pressure_scale, viscous_scale);
  constexpr double discrimination =
      256.0 * std::numeric_limits<double>::epsilon();
  if (!(refined.pressure_traction_rms_force_N >
        discrimination * pressure_scale) ||
      !(refined.viscous_traction_rms_force_N >
        discrimination * viscous_scale) ||
      !(refined.total_traction_rms_force_N > discrimination * total_scale))
    throw std::runtime_error(
        "manufactured force reference is non-discriminating");
  const auto stable_vector = [](runtime::Real3 left, runtime::Real3 right,
                                double scale) {
    return max_abs_difference(left, right) <= 1.0e-13 * std::max(1.0, scale);
  };
  const auto stable_scalar = [](double left, double right, double scale) {
    return std::abs(left - right) <= 1.0e-13 * std::max(1.0, scale);
  };
  if (!stable_vector(coarse.force.pressure_N, refined.force.pressure_N,
                     pressure_scale) ||
      !stable_vector(coarse.force.total_N, refined.force.total_N,
                     total_scale) ||
      !stable_scalar(coarse.pressure_traction_rms_force_N,
                     refined.pressure_traction_rms_force_N, pressure_scale) ||
      !stable_scalar(coarse.total_traction_rms_force_N,
                     refined.total_traction_rms_force_N, total_scale) ||
      !signed_viscous_force_reference_preflight_accepts(coarse, refined))
    throw std::runtime_error(
        "manufactured force reference failed quadrature stability");
  return refined;
}

runtime::Int3 logical_cell_from_global_id(mesh::GlobalCellId id,
                                          runtime::Int3 extent) {
  if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0)
    throw std::runtime_error(
        "manufactured pressure-extremum extent is invalid");
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto nz = static_cast<std::uint64_t>(extent.z);
  const auto plane = nx * ny;
  if (id >= plane * nz)
    throw std::runtime_error(
        "manufactured pressure-extremum global cell ID is invalid");
  return {static_cast<int>(id % nx), static_cast<int>((id / nx) % ny),
          static_cast<int>(id / plane)};
}

runtime::Int3 pressure_measure_neighbour_offset(std::size_t occurrence) {
  constexpr std::array<runtime::Int3, 6> offsets{
      runtime::Int3{-1, 0, 0}, runtime::Int3{1, 0, 0},
      runtime::Int3{0, -1, 0}, runtime::Int3{0, 1, 0},
      runtime::Int3{0, 0, -1}, runtime::Int3{0, 0, 1}};
  if (occurrence >= offsets.size())
    throw std::runtime_error(
        "manufactured pressure-measure occurrence is invalid");
  return offsets[occurrence];
}

PressureMeasureDiagnostics evaluate_pressure_measure_diagnostics(
    const runtime::MpiContext &mpi, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const immersed::WallQuadraturePlan &wall_plan,
    const immersed::GhostStencilPlan &ghost_plan,
    const immersed::LocalFlowPatternTransform &transform,
    const finite_volume::ImmersedOperatorAdapter &adapter,
    const runtime::FieldView<const double> &pressure,
    std::vector<flow::test::ImmersedFlowWallGradientSnapshot> wall_gradients,
    runtime::Real3 product_reaction_N,
    runtime::Real3 surface_pressure_force_N, const BodySpec &body,
    double time_s, runtime::Real3 analytic_pressure_force_N,
    PressureScalarDiagnostics &scalar_result) {
  const auto rows =
      finite_volume::test::ImmersedOperatorTestAccess::rows(adapter);
  std::map<immersed::ImmersedLinkId,
           finite_volume::test::ImmersedWallLinkSnapshot>
      link_geometry;
  for (const auto &row : rows)
    for (const auto &link : row.links) {
      const auto [position, inserted] = link_geometry.emplace(link.id, link);
      if (!inserted || position->second.id != link.id)
        throw std::runtime_error(
            "manufactured A22 link geometry is duplicated");
    }
  std::sort(wall_gradients.begin(), wall_gradients.end(),
            [](const auto &left, const auto &right) {
              return left.link < right.link;
            });
  if (std::adjacent_find(wall_gradients.begin(), wall_gradients.end(),
                         [](const auto &left, const auto &right) {
                           return left.link == right.link;
                         }) != wall_gradients.end())
    throw std::runtime_error(
        "manufactured pressure-measure wall gradient is duplicated");

  std::uint64_t local_max_link = 0U;
  std::uint64_t local_link_presence = 0U;
  for (const auto &row : rows)
    for (const auto &link : row.links) {
      local_max_link = std::max(local_max_link, link.id);
      ++local_link_presence;
    }
  for (const auto &point : wall_plan.local_points()) {
    local_max_link = std::max(
        local_max_link,
        immersed::detail::boundary_authority_link(point.reconstruction));
    ++local_link_presence;
  }
  std::uint64_t global_max_link = 0U;
  std::uint64_t global_link_presence = 0U;
  check_mpi(MPI_Allreduce(&local_max_link, &global_max_link, 1, MPI_UINT64_T,
                          MPI_MAX, mpi.comm()),
            "MPI_Allreduce(manufactured pressure-measure maximum link)");
  check_mpi(MPI_Allreduce(&local_link_presence, &global_link_presence, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured pressure-measure link presence)");
  if (global_link_presence == 0U ||
      global_max_link >= 6U * topology.global_cell_count() ||
      global_max_link >
          static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max() -
                                     1U))
    throw std::runtime_error(
        "manufactured pressure-measure link catalog is invalid");
  const std::size_t link_capacity =
      static_cast<std::size_t>(global_max_link + 1U);
  if (link_capacity >
      static_cast<std::size_t>(std::numeric_limits<int>::max()) / 5U)
    throw std::runtime_error(
        "manufactured pressure-measure catalog exceeds MPI int");

  constexpr std::size_t kSurfaceStride = 5U;
  std::vector<double> surface_by_link(link_capacity * kSurfaceStride, 0.0);
  std::vector<double> centroid_by_link(link_capacity * 3U, 0.0);
  std::vector<double> weight_by_link(link_capacity, 0.0);
  for (const auto &point : wall_plan.local_points()) {
    const auto link =
        immersed::detail::boundary_authority_link(point.reconstruction);
    const std::size_t base =
        static_cast<std::size_t>(link) * kSurfaceStride;
    const std::size_t centroid_base = static_cast<std::size_t>(link) * 3U;
    surface_by_link[base] -=
        point.weight_m2 * point.solid_to_fluid_normal.x;
    surface_by_link[base + 1U] -=
        point.weight_m2 * point.solid_to_fluid_normal.y;
    surface_by_link[base + 2U] -=
        point.weight_m2 * point.solid_to_fluid_normal.z;
    surface_by_link[base + 3U] += point.weight_m2;
    surface_by_link[base + 4U] += 1.0;
    centroid_by_link[centroid_base] += point.weight_m2 * point.position_m.x;
    centroid_by_link[centroid_base + 1U] +=
        point.weight_m2 * point.position_m.y;
    centroid_by_link[centroid_base + 2U] +=
        point.weight_m2 * point.position_m.z;
    weight_by_link[static_cast<std::size_t>(link)] += point.weight_m2;
  }
  sum_in_place(mpi, surface_by_link.data(), surface_by_link.size());
  sum_in_place(mpi, centroid_by_link.data(), centroid_by_link.size());
  sum_in_place(mpi, weight_by_link.data(), weight_by_link.size());

  std::vector<std::uint64_t> operator_presence(link_capacity, 0U);
  for (const auto &row : rows)
    for (const auto &link : row.links)
      ++operator_presence[static_cast<std::size_t>(link.id)];
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, operator_presence.data(),
                          static_cast<int>(operator_presence.size()),
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured pressure-measure operator links)");
  for (const auto count : operator_presence)
    if (count > 1U)
      throw std::runtime_error(
          "manufactured pressure-measure operator link is duplicated");

  std::array<double, 24> sums{};
  std::array<double, 24> scalar_sums{};
  std::vector<double> link_pressure_by_id(link_capacity, 0.0);
  std::uint64_t local_link_count = 0U;
  std::uint64_t local_supported_count = 0U;
  std::uint64_t local_missing_count = 0U;
  std::uint64_t local_multiple_count = 0U;
  const auto extent = topology.global_extent();
  for (const auto &row : rows) {
    if (row.links.empty())
      continue;
    const auto fluid_logical =
        logical_cell_from_global_id(row.active_cell, extent);
    const auto fluid_local = topology.find_local_cell(row.active_cell);
    if (!fluid_local.has_value())
      throw std::runtime_error(
          "manufactured pressure-measure fluid cell is unavailable");
    const auto fluid_center = geometry.cell_center_m(*fluid_local);
    std::array<runtime::Real3, 6> neighbour_centers{};
    for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
      const auto offset = pressure_measure_neighbour_offset(occurrence);
      const runtime::Int3 logical{fluid_logical.x + offset.x,
                                  fluid_logical.y + offset.y,
                                  fluid_logical.z + offset.z};
      const auto local =
          topology.find_local_cell(topology.global_cell_id(logical));
      if (!local.has_value())
        throw std::runtime_error(
            "manufactured pressure-measure neighbour is unavailable");
      neighbour_centers[occurrence] = geometry.cell_center_m(*local);
    }
    for (const auto &link : row.links) {
      const auto gradient = std::lower_bound(
          wall_gradients.begin(), wall_gradients.end(), link.id,
          [](const auto &candidate, std::uint64_t id) {
            return candidate.link < id;
          });
      if (gradient == wall_gradients.end() || gradient->link != link.id ||
          !std::isfinite(gradient->normal_gradient_pa_per_m))
        throw std::runtime_error(
            "manufactured pressure-measure wall gradient is unavailable");
      const auto solid_center = neighbour_centers[link.occurrence];
      const auto fluid_to_solid = subtract(solid_center, fluid_center);
      const auto fluid_to_face =
          subtract(link.pressure_quadrature_m, fluid_center);
      const double center_distance = magnitude(fluid_to_solid);
      const double owner_to_face = magnitude(fluid_to_face);
      if (!(center_distance > 0.0) || !(owner_to_face > 0.0) ||
          !(owner_to_face < center_distance))
        throw std::runtime_error(
            "manufactured pressure-measure background row is invalid");
      immersed::LocalCoefficientRow background{};
      background.neighbour[link.occurrence] =
          owner_to_face / center_distance;
      background.diagonal = 1.0 - background.neighbour[link.occurrence];
      const auto transformed = transform.transform_full(
          background, link.normal_scale, link.solid_to_fluid_normal);
      const auto &reconstruction = ghost_plan.reconstruction(link.id);
      const auto evaluate_link_pressure = [&] {
        const double wall_value = reconstruction.value(
            link.wall_intercept_m, pressure, 0U);
        const auto wall_gradient = reconstruction.gradient(
            link.wall_intercept_m, pressure, 0U);
        const auto remainder = [&](runtime::Real3 point_m) {
          return reconstruction.value(point_m, pressure, 0U) - wall_value -
                 dot(wall_gradient,
                     subtract(point_m, link.wall_intercept_m));
        };
        const auto evaluate_row =
            [&](const immersed::LocalCoefficientRow &coefficients) {
              double value =
                  coefficients.diagonal * remainder(fluid_center) +
                  coefficients.source;
              for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence)
                value += coefficients.neighbour[occurrence] *
                         remainder(neighbour_centers[occurrence]);
              return value;
            };
        return reconstruction.value(link.pressure_quadrature_m, pressure,
                                    0U) +
               evaluate_row(transformed) - evaluate_row(background);
      };
      const double numerical_face_pressure = reconstruction.value(
          link.pressure_quadrature_m, pressure, 0U);
      const double numerical_wall_pressure = reconstruction.value(
          link.wall_intercept_m, pressure, 0U);
      const double link_pressure = evaluate_link_pressure();
      if (!std::isfinite(link_pressure))
        throw std::runtime_error(
            "manufactured pressure-measure link pressure is non-finite");
      link_pressure_by_id[static_cast<std::size_t>(link.id)] = link_pressure;
      const runtime::Real3 background_area = link.area_from_fluid_m2;
      const runtime::Real3 projected_area =
          multiply(-link.signed_wall_measure_m2,
                   link.solid_to_fluid_normal);
      const std::size_t surface_base =
          static_cast<std::size_t>(link.id) * kSurfaceStride;
      const runtime::Real3 surface_area{surface_by_link[surface_base],
                                        surface_by_link[surface_base + 1U],
                                        surface_by_link[surface_base + 2U]};
      const auto background_reaction = multiply(-link_pressure, background_area);
      const auto projected_reaction = multiply(-link_pressure, projected_area);
      const auto surface_reaction = multiply(-link_pressure, surface_area);
      const auto exact_wall = evaluate_mms(body, link.wall_intercept_m, time_s);
      const double exact_face_pressure =
          evaluate_mms(body, link.pressure_quadrature_m, time_s)
              .mechanical_pressure_pa;
      const auto exact_remainder = [&](runtime::Real3 point_m) {
        const auto exact = evaluate_mms(body, point_m, time_s);
        return exact.mechanical_pressure_pa - exact_wall.mechanical_pressure_pa -
               dot(exact_wall.mechanical_pressure_gradient_pa_per_m,
                   subtract(point_m, link.wall_intercept_m));
      };
      const auto evaluate_exact_row =
          [&](const immersed::LocalCoefficientRow &coefficients) {
            double value =
                coefficients.diagonal * exact_remainder(fluid_center) +
                coefficients.source;
            for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence)
              value += coefficients.neighbour[occurrence] *
                       exact_remainder(neighbour_centers[occurrence]);
            return value;
          };
      const double exact_hybrid_pressure =
          exact_face_pressure + evaluate_exact_row(transformed) -
          evaluate_exact_row(background);
      const auto projected_defect_reaction =
          [&](double face_pressure, double corrected_pressure) {
            return add(multiply(-face_pressure, background_area),
                       multiply(-(corrected_pressure - face_pressure),
                                projected_area));
          };
      const std::array<runtime::Real3, 8> scalar_reactions{
          multiply(-numerical_face_pressure, background_area),
          multiply(-exact_face_pressure, background_area),
          multiply(-exact_wall.mechanical_pressure_pa, background_area),
          multiply(-exact_hybrid_pressure, background_area),
          projected_defect_reaction(numerical_face_pressure, link_pressure),
          projected_defect_reaction(exact_face_pressure,
                                    exact_hybrid_pressure),
          projected_defect_reaction(numerical_face_pressure,
                                    numerical_wall_pressure),
          projected_defect_reaction(exact_face_pressure,
                                    exact_wall.mechanical_pressure_pa)};
      for (std::size_t candidate = 0U; candidate < scalar_reactions.size();
           ++candidate) {
        scalar_sums[candidate * 3U] += scalar_reactions[candidate].x;
        scalar_sums[candidate * 3U + 1U] +=
            scalar_reactions[candidate].y;
        scalar_sums[candidate * 3U + 2U] +=
            scalar_reactions[candidate].z;
      }
      const std::array<runtime::Real3, 6> vectors{
          background_reaction, projected_reaction, surface_reaction,
          background_area, projected_area, surface_area};
      for (std::size_t vector = 0U; vector < vectors.size(); ++vector) {
        sums[vector * 3U] += vectors[vector].x;
        sums[vector * 3U + 1U] += vectors[vector].y;
        sums[vector * 3U + 2U] += vectors[vector].z;
      }
      sums[18] += magnitude(background_area);
      sums[19] += link.signed_wall_measure_m2;
      sums[20] += surface_by_link[surface_base + 3U];
      ++local_link_count;
      if (surface_by_link[surface_base + 3U] > 0.0)
        ++local_supported_count;
      else
        ++local_missing_count;
      if (surface_by_link[surface_base + 4U] > 1.0)
        ++local_multiple_count;
    }
  }
  sum_in_place(mpi, sums.data(), sums.size());
  sum_in_place(mpi, scalar_sums.data(), scalar_sums.size());
  sum_in_place(mpi, link_pressure_by_id.data(), link_pressure_by_id.size());
  std::uint64_t counts[4]{local_link_count, local_supported_count,
                          local_missing_count, local_multiple_count};
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, counts, 4, MPI_UINT64_T, MPI_SUM,
                          mpi.comm()),
            "MPI_Allreduce(manufactured pressure-measure counts)");
  const auto read_vector = [&](std::size_t offset) {
    return runtime::Real3{sums[offset], sums[offset + 1U], sums[offset + 2U]};
  };
  const auto background_reaction = read_vector(0U);
  const auto projected_reaction = read_vector(3U);
  const auto surface_reaction = read_vector(6U);
  runtime::Real3 authoritative_product_reaction{};
  runtime::Real3 exact_state_a22_reaction{};
  runtime::Real3 exact_state_a22_donor_reaction{};
  runtime::Real3 exact_state_a22_wall_reaction{};
  bool local_authority_discriminating = false;
  double local_authority_coefficient_mutation = 0.0;
  const auto owned_box = topology.owned_global_box();
  const auto pressure_at_global_cell = [&](mesh::GlobalCellId global_cell) {
    const auto logical = logical_cell_from_global_id(global_cell, extent);
    const runtime::Int3 index{logical.x - owned_box.begin.x,
                              logical.y - owned_box.begin.y,
                              logical.z - owned_box.begin.z};
    const auto interior = pressure.interior_extent();
    const int ghost = pressure.ghost_width();
    if (index.x < -ghost || index.x >= interior.x + ghost ||
        index.y < -ghost || index.y >= interior.y + ghost ||
        index.z < -ghost || index.z >= interior.z + ghost)
      throw std::runtime_error(
          "manufactured A22 pressure donor is outside the exchanged field");
    return pressure(index.x, index.y, index.z, 0);
  };
  std::map<mesh::GlobalCellId, double> exact_pressure_averages;
  const auto exact_pressure_average = [&](mesh::GlobalCellId global_cell) {
    const auto cached = exact_pressure_averages.find(global_cell);
    if (cached != exact_pressure_averages.end())
      return cached->second;
    const auto logical =
        logical_cell_from_global_id(global_cell, geometry.global_extent());
    const double value =
        detail::evaluate_cell_average_at_global(geometry, logical, body,
                                                time_s)
            .mechanical_pressure_pa;
    if (!std::isfinite(value))
      throw std::runtime_error(
          "manufactured exact A22 donor pressure is non-finite");
    exact_pressure_averages.emplace(global_cell, value);
    return value;
  };
  const auto interface_rows =
      finite_volume::test::ImmersedOperatorTestAccess::
          interface_pressure_force_rows(adapter);
  constexpr runtime::Real3 linear_gradient{0.7, -0.3, 0.5};
  const auto linear_center = [&](mesh::GlobalCellId cell_id) {
    const auto logical =
        logical_cell_from_global_id(cell_id, geometry.global_extent());
    constexpr std::array<runtime::Int3, 8> vertex_offsets{
        runtime::Int3{0, 0, 0}, runtime::Int3{1, 0, 0},
        runtime::Int3{1, 1, 0}, runtime::Int3{0, 1, 0},
        runtime::Int3{0, 0, 1}, runtime::Int3{1, 0, 1},
        runtime::Int3{1, 1, 1}, runtime::Int3{0, 1, 1}};
    runtime::Real3 center{};
    for (const auto &offset : vertex_offsets)
      center = add(
          center, geometry.vertex_position_m(
                      {logical.x + offset.x, logical.y + offset.y,
                       logical.z + offset.z}));
    return multiply(0.125, center);
  };
  double linear_row_residual_linf = 0.0;
  double linear_row_residual_l2 = 0.0;
  double mms_row_residual_linf = 0.0;
  double mms_row_residual_l2 = 0.0;
  std::uint64_t mms_row_count = 0U;
  double link_true_normal_mismatch_linf = 0.0;
  std::array<double, 3> mms_row_sum{};
  double mms_row_sum_abs = 0.0;
  std::array<double, 3> linear_row_sum{};
  std::array<double, 3> single_link_row_sum{};
  std::array<double, 3> multi_link_row_sum{};
  std::array<double, 3> single_link_donor_sum{};
  std::array<double, 3> single_link_wall_sum{};
  std::array<double, 3> single_link_linear_donor_sum{};
  std::array<double, 3> single_link_donor_constant_sum{};
  std::array<double, 3> single_link_donor_moment_x_sum{};
  std::array<double, 3> single_link_donor_moment_y_sum{};
  std::array<double, 3> single_link_donor_moment_z_sum{};
  std::array<double, 3> single_link_wall_coefficient_sum{};
  std::array<double, 3> single_link_donor_second_moment_sum{};
  std::array<double, 3> single_link_authority_wall_flux_sum{};
  std::array<double, 3> single_link_exact_face_flux_sum{};
  std::array<double, 3> multi_link_exact_face_flux_sum{};
  std::array<double, 3> single_link_authority_value_difference_l2{};
  std::array<double, 3> exact_centroid_quadrature_sum{};
  std::array<double, 3> authority_centroid_quadrature_sum{};
  std::array<double, 3> surface_vector_sum{};
  std::array<double, 3> ghost_measure_sum{};
  std::array<double, 3> ghost_authority_centroid_quadrature_sum{};
  std::uint64_t single_link_row_count = 0U;
  std::uint64_t multi_link_row_count = 0U;
  std::map<mesh::GlobalCellId, std::size_t> row_link_counts;
  for (const auto &row : rows)
    row_link_counts[row.active_cell] = row.links.size();
  for (const auto &row : interface_rows) {
    std::array<double, 3> residual{};
    std::array<double, 3> mms_residual{};
    for (const auto &term : row.a22_donor_terms) {
      if (term.output_component >= residual.size())
        throw std::runtime_error(
            "manufactured linear A22 donor component is invalid");
      residual[term.output_component] +=
          term.coefficient *
          dot(linear_gradient, linear_center(term.pressure_cell));
      mms_residual[term.output_component] +=
          term.coefficient * exact_pressure_average(term.pressure_cell);
    }
    for (const auto &term : row.a22_wall_terms) {
      if (term.output_component >= residual.size())
        throw std::runtime_error(
            "manufactured linear A22 wall component is invalid");
      const auto geometry_entry = link_geometry.find(term.link);
      if (geometry_entry == link_geometry.end())
        throw std::runtime_error(
            "manufactured linear A22 wall geometry is unavailable");
      residual[term.output_component] +=
          term.coefficient *
          dot(linear_gradient,
              geometry_entry->second.solid_to_fluid_normal);
      const auto exact_wall_term =
          evaluate_mms(body, geometry_entry->second.wall_intercept_m, time_s);
      mms_residual[term.output_component] +=
          term.coefficient *
          dot(exact_wall_term.mechanical_pressure_gradient_pa_per_m,
              geometry_entry->second.solid_to_fluid_normal);
    }
    const double row_norm =
        std::max({std::abs(residual[0]), std::abs(residual[1]),
                  std::abs(residual[2])});
    linear_row_residual_linf =
        std::max(linear_row_residual_linf, row_norm);
    linear_row_residual_l2 +=
        residual[0] * residual[0] + residual[1] * residual[1] +
        residual[2] * residual[2];
    const double mms_row_norm =
        std::max({std::abs(mms_residual[0]), std::abs(mms_residual[1]),
                  std::abs(mms_residual[2])});
    mms_row_residual_linf = std::max(mms_row_residual_linf, mms_row_norm);
    mms_row_residual_l2 +=
        mms_residual[0] * mms_residual[0] +
        mms_residual[1] * mms_residual[1] +
        mms_residual[2] * mms_residual[2];
    ++mms_row_count;
    mms_row_sum[0] += mms_residual[0];
    mms_row_sum[1] += mms_residual[1];
    mms_row_sum[2] += mms_residual[2];
    mms_row_sum_abs += std::abs(mms_residual[0]) + std::abs(mms_residual[1]) +
                       std::abs(mms_residual[2]);
    linear_row_sum[0] += residual[0];
    linear_row_sum[1] += residual[1];
    linear_row_sum[2] += residual[2];
    const auto found = row_link_counts.find(row.momentum_cell);
    const std::size_t link_count = found == row_link_counts.end() ? 0U
                                                                  : found->second;
    if (link_count <= 1U) {
      single_link_row_sum[0] += mms_residual[0];
      single_link_row_sum[1] += mms_residual[1];
      single_link_row_sum[2] += mms_residual[2];
      std::array<double, 3> donor_part{};
      std::array<double, 3> wall_part{};
      for (const auto &term : row.a22_donor_terms)
        if (term.output_component < 3U) {
          const auto center = linear_center(term.pressure_cell);
          donor_part[term.output_component] +=
              term.coefficient * exact_pressure_average(term.pressure_cell);
          single_link_donor_constant_sum[term.output_component] +=
              term.coefficient;
          single_link_donor_moment_x_sum[term.output_component] +=
              term.coefficient * center.x;
          single_link_donor_moment_y_sum[term.output_component] +=
              term.coefficient * center.y;
          single_link_donor_moment_z_sum[term.output_component] +=
              term.coefficient * center.z;
          single_link_donor_second_moment_sum[term.output_component] +=
              term.coefficient *
              (center.x * center.x + center.y * center.y +
               center.z * center.z);
        }
      for (const auto &term : row.a22_wall_terms) {
        if (term.output_component >= 3U)
          continue;
        const auto geometry_entry = link_geometry.find(term.link);
        if (geometry_entry == link_geometry.end())
          continue;
        const auto exact_wall = evaluate_mms(
            body, geometry_entry->second.wall_intercept_m, time_s);
        wall_part[term.output_component] +=
            term.coefficient *
            dot(exact_wall.mechanical_pressure_gradient_pa_per_m,
                geometry_entry->second.solid_to_fluid_normal);
        single_link_wall_coefficient_sum[term.output_component] +=
            term.coefficient;
      }
      single_link_donor_sum[0] += donor_part[0];
      single_link_donor_sum[1] += donor_part[1];
      single_link_donor_sum[2] += donor_part[2];
      std::array<double, 3> linear_donor_part{};
      for (const auto &term : row.a22_donor_terms)
        if (term.output_component < 3U)
          linear_donor_part[term.output_component] +=
              term.coefficient *
              dot(linear_gradient, linear_center(term.pressure_cell));
      single_link_linear_donor_sum[0] += linear_donor_part[0];
      single_link_linear_donor_sum[1] += linear_donor_part[1];
      single_link_linear_donor_sum[2] += linear_donor_part[2];
      single_link_wall_sum[0] += wall_part[0];
      single_link_wall_sum[1] += wall_part[1];
      single_link_wall_sum[2] += wall_part[2];
      if (!row.a22_wall_terms.empty()) {
        const auto &wall_term = row.a22_wall_terms.front();
        const auto geometry_entry = link_geometry.find(wall_term.link);
        if (geometry_entry != link_geometry.end()) {
          const auto &wall = geometry_entry->second;
          const auto authority =
              immersed::detail::QuadraticReconstructionWeights::
                  origin_normal_gradient_constrained_value_weights(
                      ghost_plan.reconstruction(wall_term.link),
                      wall.wall_intercept_m);
          const auto exact_wall =
              evaluate_mms(body, wall.wall_intercept_m, time_s);
          const double exact_g = dot(
              exact_wall.mechanical_pressure_gradient_pa_per_m,
              wall.solid_to_fluid_normal);
          double authority_value = authority.boundary_coefficient * exact_g;
          for (const auto &donor : authority.donors)
            authority_value +=
                donor.weight * exact_pressure_average(donor.global_cell);
          const double wall_measure =
              wall.signed_wall_measure_m2;
          single_link_authority_wall_flux_sum[0] +=
              wall_measure * wall.solid_to_fluid_normal.x * authority_value;
          single_link_authority_wall_flux_sum[1] +=
              wall_measure * wall.solid_to_fluid_normal.y * authority_value;
          single_link_authority_wall_flux_sum[2] +=
              wall_measure * wall.solid_to_fluid_normal.z * authority_value;
          const auto exact_face =
              evaluate_mms(body, wall.pressure_quadrature_m, time_s);
          single_link_exact_face_flux_sum[0] +=
              wall.area_from_fluid_m2.x *
              exact_face.mechanical_pressure_pa;
          single_link_exact_face_flux_sum[1] +=
              wall.area_from_fluid_m2.y *
              exact_face.mechanical_pressure_pa;
          single_link_exact_face_flux_sum[2] +=
              wall.area_from_fluid_m2.z *
              exact_face.mechanical_pressure_pa;
          single_link_authority_value_difference_l2[0] +=
              (authority_value - exact_face.mechanical_pressure_pa) *
              (authority_value - exact_face.mechanical_pressure_pa);
        }
      }
      ++single_link_row_count;
    } else {
      multi_link_row_sum[0] += mms_residual[0];
      multi_link_row_sum[1] += mms_residual[1];
      multi_link_row_sum[2] += mms_residual[2];
      std::vector<immersed::ImmersedLinkId> row_links;
      for (const auto &wall_term : row.a22_wall_terms)
        row_links.push_back(wall_term.link);
      std::sort(row_links.begin(), row_links.end());
      row_links.erase(std::unique(row_links.begin(), row_links.end()),
                      row_links.end());
      for (const auto link_id : row_links) {
        const auto link_geometry_entry = link_geometry.find(link_id);
        if (link_geometry_entry != link_geometry.end()) {
          const auto &wall = link_geometry_entry->second;
          const auto exact_face =
              evaluate_mms(body, wall.pressure_quadrature_m, time_s);
          multi_link_exact_face_flux_sum[0] +=
              wall.area_from_fluid_m2.x * exact_face.mechanical_pressure_pa;
          multi_link_exact_face_flux_sum[1] +=
              wall.area_from_fluid_m2.y * exact_face.mechanical_pressure_pa;
          multi_link_exact_face_flux_sum[2] +=
              wall.area_from_fluid_m2.z * exact_face.mechanical_pressure_pa;
        }
      }
      ++multi_link_row_count;
    }
  }
  const runtime::Real3 sphere_center{0.5, 0.5, 0.5};
  for (const auto &[link_id, wall] : link_geometry) {
    static_cast<void>(link_id);
    const auto offset = subtract(wall.wall_intercept_m, sphere_center);
    const auto true_normal =
        multiply(1.0 / magnitude(offset), offset);
    link_true_normal_mismatch_linf =
        std::max(link_true_normal_mismatch_linf,
                 1.0 - dot(wall.solid_to_fluid_normal, true_normal));
  }
  max_in_place(mpi, &link_true_normal_mismatch_linf, 1U);
  sum_in_place(mpi, mms_row_sum.data(), mms_row_sum.size());
  sum_in_place(mpi, linear_row_sum.data(), linear_row_sum.size());
  sum_in_place(mpi, &mms_row_sum_abs, 1U);
  sum_in_place(mpi, single_link_row_sum.data(), single_link_row_sum.size());
  sum_in_place(mpi, multi_link_row_sum.data(), multi_link_row_sum.size());
  sum_in_place(mpi, single_link_donor_sum.data(), single_link_donor_sum.size());
  sum_in_place(mpi, single_link_wall_sum.data(), single_link_wall_sum.size());
  sum_in_place(mpi, single_link_linear_donor_sum.data(),
               single_link_linear_donor_sum.size());
  sum_in_place(mpi, single_link_donor_constant_sum.data(),
               single_link_donor_constant_sum.size());
  sum_in_place(mpi, single_link_donor_moment_x_sum.data(),
               single_link_donor_moment_x_sum.size());
  sum_in_place(mpi, single_link_donor_moment_y_sum.data(),
               single_link_donor_moment_y_sum.size());
  sum_in_place(mpi, single_link_donor_moment_z_sum.data(),
               single_link_donor_moment_z_sum.size());
  sum_in_place(mpi, single_link_wall_coefficient_sum.data(),
               single_link_wall_coefficient_sum.size());
  sum_in_place(mpi, single_link_donor_second_moment_sum.data(),
               single_link_donor_second_moment_sum.size());
  sum_in_place(mpi, single_link_authority_wall_flux_sum.data(),
               single_link_authority_wall_flux_sum.size());
  sum_in_place(mpi, single_link_exact_face_flux_sum.data(),
               single_link_exact_face_flux_sum.size());
  sum_in_place(mpi, multi_link_exact_face_flux_sum.data(),
               multi_link_exact_face_flux_sum.size());
  sum_in_place(mpi, single_link_authority_value_difference_l2.data(),
               single_link_authority_value_difference_l2.size());
  std::uint64_t global_single_link_rows = single_link_row_count;
  std::uint64_t global_multi_link_rows = multi_link_row_count;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &global_single_link_rows, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured single-link row count)");
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &global_multi_link_rows, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured multi-link row count)");
  for (const auto &row : interface_rows) {
    std::array<double, 3> residual{};
    std::array<double, 3> exact_residual{};
    for (const auto &term : row.a22_donor_terms) {
      if (term.output_component >= residual.size() ||
          !std::isfinite(term.coefficient))
        throw std::runtime_error(
            "manufactured A22 pressure donor term is invalid");
      const double value = pressure_at_global_cell(term.pressure_cell);
      residual[term.output_component] += term.coefficient * value;
      exact_residual[term.output_component] +=
          term.coefficient * exact_pressure_average(term.pressure_cell);
      auto &donor_component =
          term.output_component == 0U   ? exact_state_a22_donor_reaction.x
          : term.output_component == 1U ? exact_state_a22_donor_reaction.y
                                        : exact_state_a22_donor_reaction.z;
      donor_component -=
          term.coefficient * exact_pressure_average(term.pressure_cell);
      if (term.coefficient != 0.0 && value != 0.0) {
        const double coefficient_delta =
            0.125 * std::max(1.0, std::abs(term.coefficient));
        local_authority_coefficient_mutation =
            std::max(local_authority_coefficient_mutation,
                     std::abs(coefficient_delta * value));
      }
      local_authority_discriminating =
          local_authority_discriminating || term.coefficient != 0.0;
    }
    for (const auto &term : row.a22_wall_terms) {
      if (term.output_component >= residual.size() ||
          !std::isfinite(term.coefficient))
        throw std::runtime_error(
            "manufactured A22 pressure wall term is invalid");
      const auto gradient = std::lower_bound(
          wall_gradients.begin(), wall_gradients.end(), term.link,
          [](const auto &candidate, std::uint64_t id) {
            return candidate.link < id;
          });
      if (gradient == wall_gradients.end() || gradient->link != term.link)
        throw std::runtime_error(
            "manufactured A22 pressure wall gradient is unavailable");
      residual[term.output_component] +=
          term.coefficient * gradient->normal_gradient_pa_per_m;
      const auto geometry_entry = link_geometry.find(term.link);
      if (geometry_entry == link_geometry.end())
        throw std::runtime_error(
            "manufactured exact A22 wall geometry is unavailable");
      const auto exact_wall = evaluate_mms(
          body, geometry_entry->second.wall_intercept_m, time_s);
      const double exact_normal_gradient = dot(
          exact_wall.mechanical_pressure_gradient_pa_per_m,
          geometry_entry->second.solid_to_fluid_normal);
      if (!std::isfinite(exact_normal_gradient))
        throw std::runtime_error(
            "manufactured exact A22 wall gradient is non-finite");
      exact_residual[term.output_component] +=
          term.coefficient * exact_normal_gradient;
      auto &wall_component =
          term.output_component == 0U ? exact_state_a22_wall_reaction.x
          : term.output_component == 1U
                ? exact_state_a22_wall_reaction.y
                : exact_state_a22_wall_reaction.z;
      wall_component -= term.coefficient * exact_normal_gradient;
      if (term.coefficient != 0.0 &&
          gradient->normal_gradient_pa_per_m != 0.0) {
        const double coefficient_delta =
            0.125 * std::max(1.0, std::abs(term.coefficient));
        local_authority_coefficient_mutation = std::max(
            local_authority_coefficient_mutation,
            std::abs(coefficient_delta *
                     gradient->normal_gradient_pa_per_m));
      }
      local_authority_discriminating =
          local_authority_discriminating || term.coefficient != 0.0;
    }
    authoritative_product_reaction.x -= residual[0];
    authoritative_product_reaction.y -= residual[1];
    authoritative_product_reaction.z -= residual[2];
    exact_state_a22_reaction.x -= exact_residual[0];
    exact_state_a22_reaction.y -= exact_residual[1];
    exact_state_a22_reaction.z -= exact_residual[2];
  }
  std::array<double, 12> authoritative_reaction{
      authoritative_product_reaction.x, authoritative_product_reaction.y,
      authoritative_product_reaction.z, exact_state_a22_reaction.x,
      exact_state_a22_reaction.y, exact_state_a22_reaction.z};
  authoritative_reaction[6] = exact_state_a22_donor_reaction.x;
  authoritative_reaction[7] = exact_state_a22_donor_reaction.y;
  authoritative_reaction[8] = exact_state_a22_donor_reaction.z;
  authoritative_reaction[9] = exact_state_a22_wall_reaction.x;
  authoritative_reaction[10] = exact_state_a22_wall_reaction.y;
  authoritative_reaction[11] = exact_state_a22_wall_reaction.z;
  sum_in_place(mpi, authoritative_reaction.data(),
               authoritative_reaction.size());
  authoritative_product_reaction = {authoritative_reaction[0],
                                    authoritative_reaction[1],
                                    authoritative_reaction[2]};
  exact_state_a22_reaction = {authoritative_reaction[3],
                              authoritative_reaction[4],
                              authoritative_reaction[5]};
  exact_state_a22_donor_reaction = {authoritative_reaction[6],
                                    authoritative_reaction[7],
                                    authoritative_reaction[8]};
  exact_state_a22_wall_reaction = {authoritative_reaction[9],
                                   authoritative_reaction[10],
                                   authoritative_reaction[11]};
  int authority_discriminating = local_authority_discriminating ? 1 : 0;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &authority_discriminating, 1, MPI_INT,
                          MPI_MAX, mpi.comm()),
            "MPI_Allreduce(manufactured A22 pressure authority evidence)");
  max_in_place(mpi, &local_authority_coefficient_mutation, 1U);
  if (authority_discriminating != 1 ||
      !(local_authority_coefficient_mutation > 0.0) ||
      !std::isfinite(local_authority_coefficient_mutation) ||
      !std::isfinite(authoritative_product_reaction.x) ||
      !std::isfinite(authoritative_product_reaction.y) ||
      !std::isfinite(authoritative_product_reaction.z) ||
      !std::isfinite(exact_state_a22_reaction.x) ||
      !std::isfinite(exact_state_a22_reaction.y) ||
      !std::isfinite(exact_state_a22_reaction.z))
    throw std::runtime_error(
        "manufactured A22 pressure authority evidence is invalid");
  const auto read_scalar_vector = [&](std::size_t offset) {
    return runtime::Real3{scalar_sums[offset], scalar_sums[offset + 1U],
                          scalar_sums[offset + 2U]};
  };
  scalar_result.numerical_face_reaction_N = read_scalar_vector(0U);
  scalar_result.exact_face_reaction_N = read_scalar_vector(3U);
  scalar_result.exact_wall_reaction_N = read_scalar_vector(6U);
  scalar_result.exact_hybrid_reaction_N = read_scalar_vector(9U);
  scalar_result.numerical_projected_defect_reaction_N =
      read_scalar_vector(12U);
  scalar_result.exact_projected_defect_reaction_N = read_scalar_vector(15U);
  scalar_result.numerical_wall_projected_reaction_N =
      read_scalar_vector(18U);
  scalar_result.exact_wall_projected_reaction_N = read_scalar_vector(21U);
  scalar_result.numerical_face_consistency_N =
      add(scalar_result.numerical_face_reaction_N, surface_pressure_force_N);
  scalar_result.exact_face_consistency_N =
      add(scalar_result.exact_face_reaction_N, analytic_pressure_force_N);
  scalar_result.exact_wall_consistency_N =
      add(scalar_result.exact_wall_reaction_N, analytic_pressure_force_N);
  scalar_result.exact_hybrid_consistency_N =
      add(scalar_result.exact_hybrid_reaction_N, analytic_pressure_force_N);
  scalar_result.numerical_projected_defect_consistency_N =
      add(scalar_result.numerical_projected_defect_reaction_N,
          surface_pressure_force_N);
  scalar_result.exact_projected_defect_consistency_N =
      add(scalar_result.exact_projected_defect_reaction_N,
          analytic_pressure_force_N);
  scalar_result.numerical_wall_projected_consistency_N =
      add(scalar_result.numerical_wall_projected_reaction_N,
          surface_pressure_force_N);
  scalar_result.exact_wall_projected_consistency_N =
      add(scalar_result.exact_wall_projected_reaction_N,
          analytic_pressure_force_N);
  runtime::Real3 exact_copy_surface_reaction{};
  bool sign_mutation_detected = false;
  bool measure_mutation_detected = false;
  bool association_mutation_detected = false;
  std::optional<std::size_t> association_anchor;
  for (std::size_t link = 0U; link < link_capacity; ++link) {
    if (operator_presence[link] == 0U)
      continue;
    const std::size_t surface_base = link * kSurfaceStride;
    const runtime::Real3 area{surface_by_link[surface_base],
                              surface_by_link[surface_base + 1U],
                              surface_by_link[surface_base + 2U]};
    const double link_pressure = link_pressure_by_id[link];
    exact_copy_surface_reaction =
        add(exact_copy_surface_reaction, multiply(-link_pressure, area));
    if (!sign_mutation_detected && link_pressure != 0.0 &&
        max_abs(area) > 0.0) {
      const auto mutated = add(surface_reaction,
                               multiply(2.0 * link_pressure, area));
      sign_mutation_detected =
          max_abs_difference(mutated, surface_reaction) > 0.0;
    }
    if (!measure_mutation_detected && link_pressure != 0.0) {
      const double delta =
          0.125 * std::max(surface_by_link[surface_base + 3U], 1.0e-12);
      const auto mutated =
          add(surface_reaction, {-link_pressure * delta, 0.0, 0.0});
      measure_mutation_detected =
          max_abs_difference(mutated, surface_reaction) > 0.0;
    }
    if (!association_anchor.has_value()) {
      association_anchor = link;
    } else if (!association_mutation_detected) {
      const std::size_t anchor = *association_anchor;
      const std::size_t anchor_base = anchor * kSurfaceStride;
      const runtime::Real3 anchor_area{
          surface_by_link[anchor_base], surface_by_link[anchor_base + 1U],
          surface_by_link[anchor_base + 2U]};
      const double anchor_pressure = link_pressure_by_id[anchor];
      const auto delta =
          add(multiply(anchor_pressure - link_pressure, anchor_area),
              multiply(link_pressure - anchor_pressure, area));
      association_mutation_detected = max_abs(delta) > 0.0;
    }
  }
  const double copy_scale =
      std::max({1.0, max_abs(surface_reaction),
                max_abs(exact_copy_surface_reaction)});
  if (max_abs_difference(exact_copy_surface_reaction, surface_reaction) >
          1024.0 * std::numeric_limits<double>::epsilon() * copy_scale ||
      !sign_mutation_detected || !measure_mutation_detected ||
      !association_mutation_detected)
    throw std::runtime_error(
        "manufactured pressure-measure oracle is not mutation-sensitive");
  const auto closure = [&](std::size_t vector_offset,
                           std::size_t scale_offset) {
    if (!(sums[scale_offset] > 0.0))
      throw std::runtime_error(
          "manufactured pressure-measure closure scale is invalid");
    return magnitude(read_vector(vector_offset)) / sums[scale_offset];
  };
  double local_row_wall_value_linf = 0.0;
  double local_row_wall_value_l2 = 0.0;
  double local_per_link_wall_value_linf = 0.0;
  double local_per_link_wall_value_l2 = 0.0;
  double local_per_link_origin_g_linf = 0.0;
  double local_per_link_origin_g_l2 = 0.0;
  double local_linear_wall_value_linf = 0.0;
  double local_linear_wall_value_l2 = 0.0;
  double local_linear_wall_value_hg_linf = 0.0;
  double local_plain_linear_wall_value_linf = 0.0;
  double local_linear_extrap_wall_linf = 0.0;
  double local_linear_extrap_wall_l2 = 0.0;
  double local_constrained_origin_linear_linf = 0.0;
  double local_constrained_origin_hg_linear_linf = 0.0;
  double local_constant_field_linf = 0.0;
  double local_abs_sum_w_t1 = 0.0;
  double local_abs_sum_w_t2 = 0.0;
  double local_abs_sum_w_n = 0.0;
  double local_authority_linear_linf = 0.0;
  double local_ghost_donor_wall_linf = 0.0;
  double local_ghost_donor_wall_l2 = 0.0;
  std::uint64_t local_ghost_donor_count = 0U;
  std::uint64_t local_linear_wall_value_count = 0U;
  double local_authority_wall_value_linf = 0.0;
  double local_authority_wall_value_l2 = 0.0;
  double local_authority_centroid_value_linf = 0.0;
  double local_authority_centroid_value_l2 = 0.0;
  double local_worst_link_obliqueness = 0.0;
  std::array<double, 17> correlation_sums{};
  std::uint64_t local_wall_value_count = 0U;
  for (const auto &row : rows) {
    const auto fluid_local = topology.find_local_cell(row.active_cell);
    if (!fluid_local.has_value())
      throw std::runtime_error(
          "manufactured wall-value fluid cell is absent");
    const auto fluid_center = geometry.cell_center_m(*fluid_local);
    const double h = std::cbrt(geometry.cell_volume_m3(*fluid_local));
    for (const auto &link : row.links) {
      const auto geometry_entry = link_geometry.find(link.id);
      if (geometry_entry == link_geometry.end())
        continue;
      const auto &wall = geometry_entry->second;
      const auto exact_wall =
          evaluate_mms(body, wall.wall_intercept_m, time_s);
      const double exact_p = exact_wall.mechanical_pressure_pa;
      const double exact_g = dot(
          exact_wall.mechanical_pressure_gradient_pa_per_m,
          wall.solid_to_fluid_normal);
      const auto row_weights =
          immersed::detail::QuadraticReconstructionWeights::
              origin_normal_gradient_constrained_value_weights(
                  finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                      row_reconstruction(adapter, link.id),
                  wall.wall_intercept_m);
      double row_value = row_weights.boundary_coefficient * exact_g;
      for (const auto &donor : row_weights.donors)
        row_value += donor.weight * exact_pressure_average(donor.global_cell);
      const auto link_weights =
          immersed::detail::QuadraticReconstructionWeights::
              origin_normal_gradient_constrained_value_weights(
                  ghost_plan.reconstruction(link.id), wall.wall_intercept_m);
      double link_value = link_weights.boundary_coefficient * exact_g;
      for (const auto &donor : link_weights.donors)
        link_value += donor.weight * exact_pressure_average(donor.global_cell);
      const double row_error = std::abs(row_value - exact_p);
      const double link_error = std::abs(link_value - exact_p);
      const auto exact_origin =
          evaluate_mms(body, fluid_center, time_s);
      const double origin_g = dot(
          exact_origin.mechanical_pressure_gradient_pa_per_m,
          wall.solid_to_fluid_normal);
      double origin_g_link_value =
          link_weights.boundary_coefficient * origin_g;
      for (const auto &donor : link_weights.donors)
        origin_g_link_value +=
            donor.weight * exact_pressure_average(donor.global_cell);
      const double origin_g_link_error =
          std::abs(origin_g_link_value - exact_p);
      const auto exact_origin_sample =
          evaluate_mms(body, fluid_center, time_s);
      const double linear_extrap_wall =
          exact_origin_sample.mechanical_pressure_pa - exact_g *
              magnitude(subtract(wall.wall_intercept_m, fluid_center));
      const double linear_extrap_error =
          std::abs(linear_extrap_wall - exact_p);
      local_linear_extrap_wall_linf =
          std::max(local_linear_extrap_wall_linf, linear_extrap_error);
      local_linear_extrap_wall_l2 += linear_extrap_error * linear_extrap_error;
      double linear_wall_value =
          row_weights.boundary_coefficient *
          dot(linear_gradient, wall.solid_to_fluid_normal);
      for (const auto &donor : row_weights.donors)
        linear_wall_value +=
            donor.weight * dot(linear_gradient, linear_center(donor.global_cell));
      double linear_wall_value_hg =
          row_weights.boundary_coefficient * h *
          dot(linear_gradient, wall.solid_to_fluid_normal);
      for (const auto &donor : row_weights.donors)
        linear_wall_value_hg +=
            donor.weight * dot(linear_gradient, linear_center(donor.global_cell));
      const auto plain_weights =
          immersed::detail::QuadraticReconstructionWeights::value_weights(
              finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                  row_reconstruction(adapter, link.id),
              wall.wall_intercept_m);
      double plain_linear_value = 0.0;
      for (const auto &donor : plain_weights)
        plain_linear_value +=
            donor.weight * dot(linear_gradient, linear_center(donor.global_cell));
      double constrained_linear_at_origin = 0.0;
      double constrained_linear_hg_at_origin = 0.0;
      {
        const auto origin_weights =
            immersed::detail::QuadraticReconstructionWeights::
                origin_normal_gradient_constrained_value_weights(
                    finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                        row_reconstruction(adapter, link.id),
                    fluid_center);
        constrained_linear_at_origin =
            origin_weights.boundary_coefficient *
            dot(linear_gradient, wall.solid_to_fluid_normal);
        for (const auto &donor : origin_weights.donors)
          constrained_linear_at_origin +=
              donor.weight * dot(linear_gradient, linear_center(donor.global_cell));
        constrained_linear_hg_at_origin =
            origin_weights.boundary_coefficient * h *
            dot(linear_gradient, wall.solid_to_fluid_normal);
        for (const auto &donor : origin_weights.donors)
          constrained_linear_hg_at_origin +=
              donor.weight * dot(linear_gradient, linear_center(donor.global_cell));
      }
      double constant_field_value = 0.0;
      {
        const auto origin_weights =
            immersed::detail::QuadraticReconstructionWeights::
                origin_normal_gradient_constrained_value_weights(
                    finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                        row_reconstruction(adapter, link.id),
                    wall.wall_intercept_m);
        constant_field_value = origin_weights.boundary_coefficient * 0.0;
        for (const auto &donor : origin_weights.donors)
          constant_field_value += donor.weight * 1.0;
      }
      const double linear_wall_error =
          std::abs(linear_wall_value - dot(linear_gradient, wall.wall_intercept_m));
      const double linear_wall_error_hg = std::abs(
          linear_wall_value_hg - dot(linear_gradient, wall.wall_intercept_m));
      const double plain_linear_error = std::abs(
          plain_linear_value - dot(linear_gradient, wall.wall_intercept_m));
      const double constrained_origin_error = std::abs(
          constrained_linear_at_origin - dot(linear_gradient, fluid_center));
      const double constrained_origin_hg_error = std::abs(
          constrained_linear_hg_at_origin - dot(linear_gradient, fluid_center));
      const double constant_field_error = std::abs(constant_field_value - 1.0);
      {
        const auto &row_reconstruction =
            finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                row_reconstruction(adapter, link.id);
        const auto &row_donors =
            immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
                row_reconstruction);
        const auto fluid_logical = logical_cell_from_global_id(
            row.active_cell, geometry.global_extent());
        const auto solid_offset =
            pressure_measure_neighbour_offset(link.occurrence);
        const auto solid_logical = runtime::Int3{
            fluid_logical.x + solid_offset.x, fluid_logical.y + solid_offset.y,
            fluid_logical.z + solid_offset.z};
        const auto solid_global = topology.global_cell_id(solid_logical);
        std::vector<mesh::GlobalCellId> ghost_donors(row_donors.begin(),
                                                     row_donors.end());
        ghost_donors.push_back(solid_global);
        std::sort(ghost_donors.begin(), ghost_donors.end());
        ghost_donors.erase(
            std::unique(ghost_donors.begin(), ghost_donors.end()),
            ghost_donors.end());
        std::vector<runtime::Int3> ghost_donor_cells;
        ghost_donor_cells.reserve(ghost_donors.size());
        for (const auto donor : ghost_donors)
          ghost_donor_cells.push_back(logical_cell_from_global_id(
              donor, geometry.global_extent()));
        const auto normal = normalized(wall.solid_to_fluid_normal);
        const std::array<double, 3> normal_magnitudes{
            std::abs(normal.x), std::abs(normal.y), std::abs(normal.z)};
        std::size_t selected = 0U;
        for (std::size_t axis = 1U; axis < normal_magnitudes.size(); ++axis)
          if (normal_magnitudes[axis] < normal_magnitudes[selected])
            selected = axis;
        const std::array<runtime::Real3, 3> axes{
            runtime::Real3{1.0, 0.0, 0.0}, runtime::Real3{0.0, 1.0, 0.0},
            runtime::Real3{0.0, 0.0, 1.0}};
        auto tangent1 = cross(axes[selected], normal);
        tangent1 = multiply(1.0 / magnitude(tangent1), tangent1);
        const auto tangent2 = cross(normal, tangent1);
        immersed::detail::BoundaryAuthorityCoverageScope ghost_coverage;
        const auto ghost_reconstruction =
            immersed::QuadraticReconstruction::create(
                wall.wall_intercept_m, normal, tangent1, tangent2, h,
                fluid_logical, ghost_donor_cells, topology, geometry);
        const auto ghost_weights =
            immersed::detail::QuadraticReconstructionWeights::
                origin_normal_gradient_constrained_value_weights(
                    ghost_reconstruction, wall.wall_intercept_m);
        const auto exact_fluid = evaluate_mms(body, fluid_center, time_s);
        const auto solid_center =
            geometry.cell_center_m(*topology.find_local_cell(solid_global));
        const double d_f =
            magnitude(subtract(wall.wall_intercept_m, fluid_center));
        const double d_g =
            magnitude(subtract(solid_center, wall.wall_intercept_m));
        const double ghost_value =
            exact_fluid.mechanical_pressure_pa + exact_g * (d_f + d_g);
        double ghost_wall_value =
            ghost_weights.boundary_coefficient * exact_g;
        for (const auto &donor : ghost_weights.donors) {
          const double donor_value =
              donor.global_cell == solid_global
                  ? ghost_value
                  : exact_pressure_average(donor.global_cell);
          ghost_wall_value += donor.weight * donor_value;
        }
        const double ghost_error = std::abs(ghost_wall_value - exact_p);
        local_ghost_donor_wall_linf =
            std::max(local_ghost_donor_wall_linf, ghost_error);
        local_ghost_donor_wall_l2 += ghost_error * ghost_error;
        ++local_ghost_donor_count;
      }
      double sum_w_t1 = 0.0;
      double sum_w_t2 = 0.0;
      double sum_w_n = 0.0;
      {
        const auto origin_weights =
            immersed::detail::QuadraticReconstructionWeights::
                origin_normal_gradient_constrained_value_weights(
                    finite_volume::detail::ImmersedBoundaryAuthorityAccess::
                        row_reconstruction(adapter, link.id),
                    fluid_center);
        for (const auto &donor : origin_weights.donors) {
          const auto center = linear_center(donor.global_cell);
          const auto offset = subtract(center, fluid_center);
          sum_w_t1 += donor.weight *
                      offset.y;
          sum_w_t2 += donor.weight *
                      offset.z;
          sum_w_n += donor.weight * offset.x;
        }
      }
      local_linear_wall_value_linf =
          std::max(local_linear_wall_value_linf, linear_wall_error);
      local_linear_wall_value_hg_linf =
          std::max(local_linear_wall_value_hg_linf, linear_wall_error_hg);
      local_plain_linear_wall_value_linf =
          std::max(local_plain_linear_wall_value_linf, plain_linear_error);
      local_constrained_origin_linear_linf =
          std::max(local_constrained_origin_linear_linf,
                   constrained_origin_error);
      local_constrained_origin_hg_linear_linf =
          std::max(local_constrained_origin_hg_linear_linf,
                   constrained_origin_hg_error);
      local_constant_field_linf =
          std::max(local_constant_field_linf, constant_field_error);
      local_abs_sum_w_t1 = std::max(local_abs_sum_w_t1, std::abs(sum_w_t1));
      local_abs_sum_w_t2 = std::max(local_abs_sum_w_t2, std::abs(sum_w_t2));
      local_abs_sum_w_n = std::max(local_abs_sum_w_n, std::abs(sum_w_n));
      local_linear_wall_value_l2 += linear_wall_error * linear_wall_error;
      ++local_linear_wall_value_count;
      local_row_wall_value_linf =
          std::max(local_row_wall_value_linf, row_error);
      local_row_wall_value_l2 += row_error * row_error;
      local_per_link_wall_value_linf =
          std::max(local_per_link_wall_value_linf, link_error);
      local_per_link_wall_value_l2 += link_error * link_error;
      local_per_link_origin_g_linf =
          std::max(local_per_link_origin_g_linf, origin_g_link_error);
      local_per_link_origin_g_l2 += origin_g_link_error * origin_g_link_error;
      const double obliqueness =
          1.0 - std::max({std::abs(wall.solid_to_fluid_normal.x),
                          std::abs(wall.solid_to_fluid_normal.y),
                          std::abs(wall.solid_to_fluid_normal.z)});
      const double delta_over_h =
          magnitude(subtract(wall.wall_intercept_m, fluid_center)) / h;
      const auto &normal = wall.solid_to_fluid_normal;
      const auto &hessian = exact_wall.mechanical_pressure_hessian_pa_per_m2;
      const auto component = [](runtime::Real3 value, std::size_t index) {
        return index == 0U ? value.x : index == 1U ? value.y : value.z;
      };
      double curvature = 0.0;
      for (std::size_t i = 0U; i < 3U; ++i)
        for (std::size_t j = 0U; j < 3U; ++j)
          curvature += component(normal, i) * hessian[i][j] *
                       component(normal, j);
      if (row_error >= local_row_wall_value_linf)
        local_worst_link_obliqueness =
            std::max(local_worst_link_obliqueness, obliqueness);
      correlation_sums[0] += link_error;
      correlation_sums[1] += delta_over_h;
      correlation_sums[2] += link_error * delta_over_h;
      correlation_sums[3] += delta_over_h * delta_over_h;
      correlation_sums[4] += link_error * link_error;
      correlation_sums[5] += obliqueness;
      correlation_sums[6] += link_error * obliqueness;
      correlation_sums[7] += obliqueness * obliqueness;
      correlation_sums[8] += curvature;
      correlation_sums[9] += link_error * curvature;
      correlation_sums[10] += curvature * curvature;
      if (link_error > 0.0 && delta_over_h > 0.0) {
        const double ln_delta = std::log(delta_over_h);
        const double ln_error = std::log(link_error);
        correlation_sums[11] += ln_delta;
        correlation_sums[12] += ln_error;
        correlation_sums[13] += ln_delta * ln_delta;
        correlation_sums[14] += ln_error * ln_error;
        correlation_sums[15] += ln_delta * ln_error;
      }
      ++local_wall_value_count;
    }
  }
  for (const auto &point : wall_plan.local_points()) {
    const auto authority =
        immersed::detail::boundary_authority_reconstruction(
            point.reconstruction);
    const auto weights =
        immersed::detail::QuadraticReconstructionWeights::
            origin_normal_gradient_constrained_value_weights(
                authority, point.position_m);
    const auto exact_wall = evaluate_mms(body, point.position_m, time_s);
    const double exact_p = exact_wall.mechanical_pressure_pa;
    const double exact_g =
        dot(exact_wall.mechanical_pressure_gradient_pa_per_m,
            point.solid_to_fluid_normal);
    double value = weights.boundary_coefficient * exact_g;
    for (const auto &donor : weights.donors)
      value += donor.weight * exact_pressure_average(donor.global_cell);
    const double error = std::abs(value - exact_p);
    local_authority_wall_value_linf =
        std::max(local_authority_wall_value_linf, error);
    local_authority_wall_value_l2 += error * error;
  }
  for (const auto &[link_id, wall] : link_geometry) {
    const std::size_t base = static_cast<std::size_t>(link_id) * 3U;
    const double patch_weight = weight_by_link[static_cast<std::size_t>(link_id)];
    if (!(patch_weight > 0.0))
      continue;
    const runtime::Real3 centroid{
        centroid_by_link[base] / patch_weight,
        centroid_by_link[base + 1U] / patch_weight,
        centroid_by_link[base + 2U] / patch_weight};
    const std::size_t surface_base =
        static_cast<std::size_t>(link_id) * kSurfaceStride;
    const runtime::Real3 surface_vector{
        surface_by_link[surface_base], surface_by_link[surface_base + 1U],
        surface_by_link[surface_base + 2U]};
    const auto exact_centroid_value =
        evaluate_mms(body, centroid, time_s).mechanical_pressure_pa;
    exact_centroid_quadrature_sum[0] +=
        surface_vector.x * exact_centroid_value;
    exact_centroid_quadrature_sum[1] +=
        surface_vector.y * exact_centroid_value;
    exact_centroid_quadrature_sum[2] +=
        surface_vector.z * exact_centroid_value;
    surface_vector_sum[0] += surface_vector.x;
    surface_vector_sum[1] += surface_vector.y;
    surface_vector_sum[2] += surface_vector.z;
    // Force-row wall value: the per-link authority reconstruction's
    // extrapolated value at the body-surface patch centroid (unconstrained;
    // the force row's wall-face value functional).
    const auto weights =
        immersed::detail::QuadraticReconstructionWeights::value_weights(
            ghost_plan.reconstruction(link_id), centroid);
    const auto exact_centroid = evaluate_mms(body, centroid, time_s);
    double value = 0.0;
    for (const auto &donor : weights)
      value += donor.weight * exact_pressure_average(donor.global_cell);
    const double error =
        std::abs(value - exact_centroid.mechanical_pressure_pa);
    authority_centroid_quadrature_sum[0] += surface_vector.x * value;
    authority_centroid_quadrature_sum[1] += surface_vector.y * value;
    authority_centroid_quadrature_sum[2] += surface_vector.z * value;
    ghost_measure_sum[0] += wall.surface_measure_m2.x;
    ghost_measure_sum[1] += wall.surface_measure_m2.y;
    ghost_measure_sum[2] += wall.surface_measure_m2.z;
    const auto ghost_weights =
        immersed::detail::QuadraticReconstructionWeights::value_weights(
            ghost_plan.reconstruction(link_id),
            wall.surface_patch_centroid_m);
    double ghost_value = 0.0;
    for (const auto &donor : ghost_weights)
      ghost_value += donor.weight * exact_pressure_average(donor.global_cell);
    ghost_authority_centroid_quadrature_sum[0] +=
        wall.surface_measure_m2.x * ghost_value;
    ghost_authority_centroid_quadrature_sum[1] +=
        wall.surface_measure_m2.y * ghost_value;
    ghost_authority_centroid_quadrature_sum[2] +=
        wall.surface_measure_m2.z * ghost_value;
    local_authority_centroid_value_linf =
        std::max(local_authority_centroid_value_linf, error);
    local_authority_centroid_value_l2 += error * error;
  }
  for (const auto &point : wall_plan.local_points()) {
    const auto authority =
        immersed::detail::boundary_authority_reconstruction(
            point.reconstruction);
    const auto weights =
        immersed::detail::QuadraticReconstructionWeights::
            origin_normal_gradient_constrained_value_weights(
                authority, point.position_m);
    double value = weights.boundary_coefficient *
                   dot(linear_gradient, point.solid_to_fluid_normal);
    for (const auto &donor : weights.donors)
      value +=
          donor.weight * dot(linear_gradient, linear_center(donor.global_cell));
    local_authority_linear_linf =
        std::max(local_authority_linear_linf,
                 std::abs(value - dot(linear_gradient, point.position_m)));
  }
  max_in_place(mpi, &local_row_wall_value_linf, 1U);
  sum_in_place(mpi, &local_row_wall_value_l2, 1U);
  max_in_place(mpi, &local_per_link_wall_value_linf, 1U);
  sum_in_place(mpi, &local_per_link_wall_value_l2, 1U);
  max_in_place(mpi, &local_per_link_origin_g_linf, 1U);
  sum_in_place(mpi, &local_per_link_origin_g_l2, 1U);
  max_in_place(mpi, &local_linear_wall_value_linf, 1U);
  max_in_place(mpi, &local_linear_wall_value_hg_linf, 1U);
  max_in_place(mpi, &local_plain_linear_wall_value_linf, 1U);
  max_in_place(mpi, &local_linear_extrap_wall_linf, 1U);
  sum_in_place(mpi, &local_linear_extrap_wall_l2, 1U);
  max_in_place(mpi, &local_constrained_origin_linear_linf, 1U);
  max_in_place(mpi, &local_constrained_origin_hg_linear_linf, 1U);
  max_in_place(mpi, &local_constant_field_linf, 1U);
  max_in_place(mpi, &local_abs_sum_w_t1, 1U);
  max_in_place(mpi, &local_abs_sum_w_t2, 1U);
  max_in_place(mpi, &local_abs_sum_w_n, 1U);
  max_in_place(mpi, &local_authority_linear_linf, 1U);
  max_in_place(mpi, &local_ghost_donor_wall_linf, 1U);
  sum_in_place(mpi, &local_ghost_donor_wall_l2, 1U);
  std::uint64_t global_ghost_donor_count = local_ghost_donor_count;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &global_ghost_donor_count, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured ghost-donor link count)");
  sum_in_place(mpi, &local_linear_wall_value_l2, 1U);
  max_in_place(mpi, &local_authority_wall_value_linf, 1U);
  sum_in_place(mpi, &local_authority_wall_value_l2, 1U);
  max_in_place(mpi, &local_authority_centroid_value_linf, 1U);
  sum_in_place(mpi, &local_authority_centroid_value_l2, 1U);
  sum_in_place(mpi, exact_centroid_quadrature_sum.data(),
               exact_centroid_quadrature_sum.size());
  sum_in_place(mpi, authority_centroid_quadrature_sum.data(),
               authority_centroid_quadrature_sum.size());
  sum_in_place(mpi, surface_vector_sum.data(), surface_vector_sum.size());
  sum_in_place(mpi, ghost_measure_sum.data(), ghost_measure_sum.size());
  sum_in_place(mpi, ghost_authority_centroid_quadrature_sum.data(),
               ghost_authority_centroid_quadrature_sum.size());
  max_in_place(mpi, &local_worst_link_obliqueness, 1U);
  sum_in_place(mpi, correlation_sums.data(), correlation_sums.size());
  std::uint64_t global_wall_value_count = local_wall_value_count;
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, &global_wall_value_count, 1,
                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(manufactured wall-value link count)");
  const auto pearson = [&](std::size_t x_offset, std::size_t xy_offset,
                           std::size_t x2_offset, std::size_t y2_offset,
                           std::uint64_t count) {
    if (count < 2U)
      return 0.0;
    const double n = static_cast<double>(count);
    const double sum_e = correlation_sums[0];
    const double sum_x = correlation_sums[x_offset];
    const double sum_xy = correlation_sums[xy_offset];
    const double sum_x2 = correlation_sums[x2_offset];
    const double sum_e2 = correlation_sums[y2_offset];
    const double numerator = n * sum_xy - sum_x * sum_e;
    const double denominator =
        std::sqrt(std::max(0.0, (n * sum_x2 - sum_x * sum_x) *
                                    (n * sum_e2 - sum_e * sum_e)));
    return denominator > 0.0 ? numerator / denominator : 0.0;
  };
  const double error_delta_pearson =
      pearson(1U, 2U, 3U, 4U, global_wall_value_count);
  const double error_obliqueness_pearson =
      pearson(5U, 6U, 7U, 4U, global_wall_value_count);
  const double error_curvature_pearson =
      pearson(8U, 9U, 10U, 4U, global_wall_value_count);
  double error_delta_loglog_slope = 0.0;
  if (correlation_sums[11] != 0.0) {
    const double n = static_cast<double>(global_wall_value_count);
    const double denominator =
        n * correlation_sums[13] - correlation_sums[11] * correlation_sums[11];
    if (denominator > 0.0)
      error_delta_loglog_slope =
          (n * correlation_sums[15] - correlation_sums[11] *
                                          correlation_sums[12]) /
          denominator;
  }
  PressureMeasureDiagnostics result;
  result.background_reaction_N = background_reaction;
  result.projected_reaction_N = projected_reaction;
  result.surface_partition_reaction_N = surface_reaction;
  result.exact_state_a22_reaction_N = exact_state_a22_reaction;
  result.exact_state_a22_donor_reaction_N = exact_state_a22_donor_reaction;
  result.exact_state_a22_wall_reaction_N = exact_state_a22_wall_reaction;
  result.background_consistency_N =
      add(background_reaction, surface_pressure_force_N);
  result.projected_consistency_N =
      add(projected_reaction, surface_pressure_force_N);
  result.surface_partition_consistency_N =
      add(surface_reaction, surface_pressure_force_N);
  result.exact_state_a22_consistency_N =
      add(exact_state_a22_reaction, analytic_pressure_force_N);
  result.numerical_state_a22_defect_N =
      subtract(product_reaction_N, exact_state_a22_reaction);
  result.row_wall_value_linf = local_row_wall_value_linf;
  result.row_wall_value_l2 = std::sqrt(local_row_wall_value_l2);
  result.per_link_wall_value_linf = local_per_link_wall_value_linf;
  result.per_link_wall_value_l2 = std::sqrt(local_per_link_wall_value_l2);
  result.per_link_origin_g_wall_value_linf = local_per_link_origin_g_linf;
  result.per_link_origin_g_wall_value_l2 =
      std::sqrt(local_per_link_origin_g_l2);
  result.linear_wall_value_linf = local_linear_wall_value_linf;
  result.linear_wall_value_hg_linf = local_linear_wall_value_hg_linf;
  result.plain_linear_wall_value_linf = local_plain_linear_wall_value_linf;
  result.linear_extrap_wall_linf = local_linear_extrap_wall_linf;
  result.linear_extrap_wall_l2 = std::sqrt(local_linear_extrap_wall_l2);
  result.constrained_origin_linear_linf =
      local_constrained_origin_linear_linf;
  result.constrained_origin_hg_linear_linf =
      local_constrained_origin_hg_linear_linf;
  result.constant_field_linf = local_constant_field_linf;
  result.abs_sum_w_t1_linf = local_abs_sum_w_t1;
  result.abs_sum_w_t2_linf = local_abs_sum_w_t2;
  result.abs_sum_w_n_linf = local_abs_sum_w_n;
  result.authority_linear_linf = local_authority_linear_linf;
  result.ghost_donor_wall_linf = local_ghost_donor_wall_linf;
  result.ghost_donor_wall_l2 = std::sqrt(local_ghost_donor_wall_l2);
  result.ghost_donor_link_count = global_ghost_donor_count;
  result.linear_row_residual_linf = linear_row_residual_linf;
  result.linear_row_residual_l2 = std::sqrt(linear_row_residual_l2);
  result.mms_row_residual_linf = mms_row_residual_linf;
  result.mms_row_residual_l2 = std::sqrt(mms_row_residual_l2);
  result.mms_row_count = mms_row_count;
  result.link_true_normal_mismatch_linf = link_true_normal_mismatch_linf;
  const double mms_row_sum_norm = magnitude(
      runtime::Real3{mms_row_sum[0], mms_row_sum[1], mms_row_sum[2]});
  result.mms_row_coherence_ratio =
      mms_row_sum_abs > 0.0 ? mms_row_sum_norm / mms_row_sum_abs : 0.0;
  result.linear_row_sum_N = {
      linear_row_sum[0], linear_row_sum[1], linear_row_sum[2]};
  result.single_link_row_sum_N = {
      single_link_row_sum[0], single_link_row_sum[1], single_link_row_sum[2]};
  result.multi_link_row_sum_N = {
      multi_link_row_sum[0], multi_link_row_sum[1], multi_link_row_sum[2]};
  result.single_link_donor_sum_N = {
      single_link_donor_sum[0], single_link_donor_sum[1],
      single_link_donor_sum[2]};
  result.single_link_wall_sum_N = {
      single_link_wall_sum[0], single_link_wall_sum[1], single_link_wall_sum[2]};
  result.single_link_linear_donor_sum_N = {
      single_link_linear_donor_sum[0], single_link_linear_donor_sum[1],
      single_link_linear_donor_sum[2]};
  result.single_link_donor_constant_sum_N = {
      single_link_donor_constant_sum[0], single_link_donor_constant_sum[1],
      single_link_donor_constant_sum[2]};
  result.single_link_donor_moment_x_sum_N = {
      single_link_donor_moment_x_sum[0], single_link_donor_moment_x_sum[1],
      single_link_donor_moment_x_sum[2]};
  result.single_link_donor_moment_y_sum_N = {
      single_link_donor_moment_y_sum[0], single_link_donor_moment_y_sum[1],
      single_link_donor_moment_y_sum[2]};
  result.single_link_donor_moment_z_sum_N = {
      single_link_donor_moment_z_sum[0], single_link_donor_moment_z_sum[1],
      single_link_donor_moment_z_sum[2]};
  result.single_link_wall_coefficient_sum_N = {
      single_link_wall_coefficient_sum[0],
      single_link_wall_coefficient_sum[1],
      single_link_wall_coefficient_sum[2]};
  result.single_link_donor_second_moment_sum_N = {
      single_link_donor_second_moment_sum[0],
      single_link_donor_second_moment_sum[1],
      single_link_donor_second_moment_sum[2]};
  result.single_link_authority_wall_flux_sum_N = {
      single_link_authority_wall_flux_sum[0],
      single_link_authority_wall_flux_sum[1],
      single_link_authority_wall_flux_sum[2]};
  result.single_link_exact_face_flux_sum_N = {
      single_link_exact_face_flux_sum[0], single_link_exact_face_flux_sum[1],
      single_link_exact_face_flux_sum[2]};
  result.multi_link_exact_face_flux_sum_N = {
      multi_link_exact_face_flux_sum[0], multi_link_exact_face_flux_sum[1],
      multi_link_exact_face_flux_sum[2]};
  result.single_link_authority_value_difference_l2_N = {
      std::sqrt(single_link_authority_value_difference_l2[0]),
      std::sqrt(single_link_authority_value_difference_l2[1]),
      std::sqrt(single_link_authority_value_difference_l2[2])};
  result.single_link_row_count = global_single_link_rows;
  result.multi_link_row_count = global_multi_link_rows;
  result.linear_wall_value_l2 = std::sqrt(local_linear_wall_value_l2);
  result.linear_wall_value_count = local_linear_wall_value_count;
  result.authority_wall_value_linf = local_authority_wall_value_linf;
  result.authority_wall_value_l2 = std::sqrt(local_authority_wall_value_l2);
  result.authority_centroid_value_linf =
      local_authority_centroid_value_linf;
  result.authority_centroid_value_l2 =
      std::sqrt(local_authority_centroid_value_l2);
  result.exact_centroid_quadrature_sum_N = {
      exact_centroid_quadrature_sum[0], exact_centroid_quadrature_sum[1],
      exact_centroid_quadrature_sum[2]};
  result.authority_centroid_quadrature_sum_N = {
      authority_centroid_quadrature_sum[0],
      authority_centroid_quadrature_sum[1],
      authority_centroid_quadrature_sum[2]};
  result.surface_vector_sum_N = {surface_vector_sum[0], surface_vector_sum[1],
                                 surface_vector_sum[2]};
  result.ghost_measure_sum_N = {ghost_measure_sum[0], ghost_measure_sum[1],
                                ghost_measure_sum[2]};
  result.ghost_authority_centroid_quadrature_sum_N = {
      ghost_authority_centroid_quadrature_sum[0],
      ghost_authority_centroid_quadrature_sum[1],
      ghost_authority_centroid_quadrature_sum[2]};
  result.error_delta_pearson = error_delta_pearson;
  result.error_obliqueness_pearson = error_obliqueness_pearson;
  result.error_curvature_pearson = error_curvature_pearson;
  result.error_delta_loglog_slope = error_delta_loglog_slope;
  result.worst_link_obliqueness = local_worst_link_obliqueness;
  result.wall_value_link_count = global_wall_value_count;
  result.background_area_closure = closure(9U, 18U);
  result.projected_area_closure = closure(12U, 19U);
  result.surface_partition_area_closure = closure(15U, 20U);
  result.link_count = counts[0];
  result.surface_supported_link_count = counts[1];
  result.missing_surface_link_count = counts[2];
  result.multiple_surface_point_link_count = counts[3];
  const double reproduction_scale = std::max(
      {1.0, max_abs(authoritative_product_reaction),
       max_abs(product_reaction_N)});
  result.product_reproduction_linf_N = max_abs_difference(
      authoritative_product_reaction, product_reaction_N);
  if (result.product_reproduction_linf_N >
      1024.0 * std::numeric_limits<double>::epsilon() * reproduction_scale)
    throw std::runtime_error(
        "manufactured A22 pressure row did not reproduce the product reaction");
  return result;
}

PressureErrorExtremum reduce_pressure_error_extremum(
    const runtime::MpiContext &mpi,
    const std::vector<PressureErrorExtremum> &local_candidates,
    runtime::Int3 global_extent) {
  const auto local = select_pressure_error_extremum(local_candidates);
  int local_valid = local_candidates.empty() || local.has_value() ? 1 : 0;
  int globally_valid = 0;
  check_mpi(MPI_Allreduce(&local_valid, &globally_valid, 1, MPI_INT, MPI_MIN,
                          mpi.comm()),
            "MPI_Allreduce(manufactured pressure-extremum validity)");
  if (globally_valid != 1)
    throw std::runtime_error(
        "manufactured pressure-extremum candidate is invalid");

  int local_present = local.has_value() ? 1 : 0;
  int present_ranks = 0;
  check_mpi(MPI_Allreduce(&local_present, &present_ranks, 1, MPI_INT, MPI_SUM,
                          mpi.comm()),
            "MPI_Allreduce(manufactured pressure-extremum presence)");
  if (present_ranks <= 0)
    throw std::runtime_error(
        "manufactured pressure-extremum active support is empty");

  const double local_absolute = local.has_value()
                                    ? local->absolute_error_pa
                                    : -std::numeric_limits<double>::infinity();
  double global_absolute = 0.0;
  check_mpi(MPI_Allreduce(&local_absolute, &global_absolute, 1, MPI_DOUBLE,
                          MPI_MAX, mpi.comm()),
            "MPI_Allreduce(manufactured pressure-extremum maximum)");
  if (!std::isfinite(global_absolute) || !(global_absolute >= 0.0))
    throw std::runtime_error(
        "manufactured pressure-extremum maximum is invalid");

  const mesh::GlobalCellId local_global_id =
      local.has_value() && local->absolute_error_pa == global_absolute
          ? local->global_cell_id
          : std::numeric_limits<mesh::GlobalCellId>::max();
  mesh::GlobalCellId global_id = 0U;
  check_mpi(MPI_Allreduce(&local_global_id, &global_id, 1, MPI_UINT64_T,
                          MPI_MIN, mpi.comm()),
            "MPI_Allreduce(manufactured pressure-extremum global ID)");

  const int local_owner = local.has_value() &&
                                  local->absolute_error_pa == global_absolute &&
                                  local->global_cell_id == global_id
                              ? 1
                              : 0;
  int owner_count = 0;
  check_mpi(MPI_Allreduce(&local_owner, &owner_count, 1, MPI_INT, MPI_SUM,
                          mpi.comm()),
            "MPI_Allreduce(manufactured pressure-extremum owner count)");
  if (owner_count != 1)
    throw std::runtime_error(
        "manufactured pressure-extremum ownership is not unique");
  const int local_owner_rank = local_owner == 1 ? mpi.rank() : mpi.size();
  int owner_rank = 0;
  check_mpi(MPI_Allreduce(&local_owner_rank, &owner_rank, 1, MPI_INT, MPI_MIN,
                          mpi.comm()),
            "MPI_Allreduce(manufactured pressure-extremum owner rank)");
  if (owner_rank < 0 || owner_rank >= mpi.size())
    throw std::runtime_error(
        "manufactured pressure-extremum owner rank is invalid");

  PressureErrorExtremum result;
  if (local_owner == 1)
    result = *local;
  std::array<double, 3> values{result.signed_error_pa, result.absolute_error_pa,
                               result.wall_distance_m};
  std::array<int, 3> logical{result.logical_cell.x, result.logical_cell.y,
                             result.logical_cell.z};
  check_mpi(MPI_Bcast(values.data(), static_cast<int>(values.size()),
                      MPI_DOUBLE, owner_rank, mpi.comm()),
            "MPI_Bcast(manufactured pressure-extremum values)");
  check_mpi(MPI_Bcast(&result.global_cell_id, 1, MPI_UINT64_T, owner_rank,
                      mpi.comm()),
            "MPI_Bcast(manufactured pressure-extremum global ID)");
  check_mpi(MPI_Bcast(logical.data(), static_cast<int>(logical.size()), MPI_INT,
                      owner_rank, mpi.comm()),
            "MPI_Bcast(manufactured pressure-extremum logical cell)");
  result.signed_error_pa = values[0];
  result.absolute_error_pa = values[1];
  result.wall_distance_m = values[2];
  result.logical_cell = {logical[0], logical[1], logical[2]};

  const auto validated = select_pressure_error_extremum(
      std::vector<PressureErrorExtremum>{result});
  const auto expected_logical =
      logical_cell_from_global_id(result.global_cell_id, global_extent);
  if (!validated.has_value() || result.absolute_error_pa != global_absolute ||
      result.global_cell_id != global_id ||
      result.logical_cell.x != expected_logical.x ||
      result.logical_cell.y != expected_logical.y ||
      result.logical_cell.z != expected_logical.z)
    throw std::runtime_error(
        "manufactured pressure-extremum broadcast is inconsistent");
  return result;
}

} // namespace

AnalyticForceReference verified_force_reference_preflight(const BodySpec &body,
                                                          double time_s) {
  return verified_force_reference(body, time_s);
}

std::optional<PressureErrorExtremum> select_pressure_error_extremum(
    const std::vector<PressureErrorExtremum> &candidates) noexcept {
  if (candidates.empty())
    return std::nullopt;
  std::optional<PressureErrorExtremum> selected;
  for (const auto &candidate : candidates) {
    if (!std::isfinite(candidate.signed_error_pa) ||
        !std::isfinite(candidate.absolute_error_pa) ||
        !(candidate.absolute_error_pa >= 0.0) ||
        candidate.absolute_error_pa != std::abs(candidate.signed_error_pa) ||
        !std::isfinite(candidate.wall_distance_m) ||
        !(candidate.wall_distance_m >= 0.0))
      return std::nullopt;
    if (!selected.has_value() ||
        candidate.absolute_error_pa > selected->absolute_error_pa ||
        (candidate.absolute_error_pa == selected->absolute_error_pa &&
         candidate.global_cell_id < selected->global_cell_id))
      selected = candidate;
  }
  return selected;
}

std::vector<PressureErrorExtremum> select_pressure_error_extrema(
    const std::vector<PressureErrorExtremum> &candidates,
    std::size_t limit) noexcept {
  if (limit == 0U || candidates.empty())
    return {};
  std::vector<PressureErrorExtremum> selected = candidates;
  for (const auto &candidate : selected) {
    if (!std::isfinite(candidate.signed_error_pa) ||
        !std::isfinite(candidate.absolute_error_pa) ||
        !(candidate.absolute_error_pa >= 0.0) ||
        candidate.absolute_error_pa != std::abs(candidate.signed_error_pa) ||
        !std::isfinite(candidate.wall_distance_m) ||
        !(candidate.wall_distance_m >= 0.0))
      return {};
  }
  std::sort(selected.begin(), selected.end(),
            [](const auto &left, const auto &right) {
              if (left.absolute_error_pa != right.absolute_error_pa)
                return left.absolute_error_pa > right.absolute_error_pa;
              return left.global_cell_id < right.global_cell_id;
            });
  if (selected.size() > limit)
    selected.resize(limit);
  return selected;
}

std::optional<std::size_t>
normalized_near_wall_pressure_band(double wall_distance_over_h) noexcept {
  if (!std::isfinite(wall_distance_over_h) || wall_distance_over_h < 0.0 ||
      wall_distance_over_h > 2.0)
    return std::nullopt;
  if (wall_distance_over_h <= 0.5)
    return 0U;
  if (wall_distance_over_h <= 1.0)
    return 1U;
  if (wall_distance_over_h <= 1.5)
    return 2U;
  return 3U;
}

bool manufactured_near_wall_band_contains(double wall_distance_m) noexcept {
  // The formal matrix starts at 12^3.  Freezing the original two-cell coarse
  // support gives a physical band that is identical across all refinements,
  // as required by the immersed-flow accuracy specification.
  constexpr double thickness_m = 2.0 * kReferenceLengthM / 12.0;
  return std::isfinite(wall_distance_m) && wall_distance_m >= 0.0 &&
         wall_distance_m <= thickness_m;
}

bool accumulate_pressure_error_moments(PressureErrorMoments &moments,
                                       double signed_error_pa,
                                       double volume_m3) noexcept {
  if (!std::isfinite(signed_error_pa) || !std::isfinite(volume_m3) ||
      !(volume_m3 > 0.0) ||
      moments.cell_count == std::numeric_limits<std::uint64_t>::max())
    return false;
  const double signed_increment = signed_error_pa * volume_m3;
  const double squared_increment =
      signed_error_pa * signed_error_pa * volume_m3;
  const double signed_sum =
      moments.signed_error_volume_pa_m3 + signed_increment;
  const double squared_sum =
      moments.squared_error_volume_pa2_m3 + squared_increment;
  const double volume_sum = moments.volume_m3 + volume_m3;
  if (!std::isfinite(signed_sum) || !std::isfinite(squared_sum) ||
      !std::isfinite(volume_sum))
    return false;
  moments.signed_error_volume_pa_m3 = signed_sum;
  moments.squared_error_volume_pa2_m3 = squared_sum;
  moments.volume_m3 = volume_sum;
  ++moments.cell_count;
  return true;
}

std::optional<PressureErrorStatistics>
summarize_pressure_error_moments(const PressureErrorMoments &moments) noexcept {
  if (!std::isfinite(moments.signed_error_volume_pa_m3) ||
      !std::isfinite(moments.squared_error_volume_pa2_m3) ||
      !std::isfinite(moments.volume_m3) || !(moments.volume_m3 > 0.0) ||
      moments.squared_error_volume_pa2_m3 < 0.0 || moments.cell_count == 0U)
    return std::nullopt;
  const double mean =
      moments.signed_error_volume_pa_m3 / moments.volume_m3;
  const double mean_square =
      moments.squared_error_volume_pa2_m3 / moments.volume_m3;
  if (!std::isfinite(mean) || !std::isfinite(mean_square) ||
      mean_square < 0.0)
    return std::nullopt;
  const double variance = std::max(0.0, mean_square - mean * mean);
  return PressureErrorStatistics{mean, std::sqrt(mean_square),
                                 std::sqrt(variance), moments.volume_m3,
                                 moments.cell_count};
}

ManufacturedRunResult
run_manufactured_case(const runtime::MpiContext &mpi,
                      const ManufacturedCase &definition) {
  if (definition.cells <= 0)
    throw std::runtime_error("manufactured cell count must be positive");
  if (definition.process_grid.x * definition.process_grid.y *
          definition.process_grid.z !=
      mpi.size())
    throw std::runtime_error(
        "manufactured process grid does not match communicator size");
  if (definition.final_force_failure_rank.has_value() &&
      (!definition.collect_force || *definition.final_force_failure_rank < 0 ||
       *definition.final_force_failure_rank >= mpi.size()))
    throw std::runtime_error(
        "manufactured final-force failure rank is invalid");

  const double surface_h =
      definition.surface_policy == ManufacturedSurfacePolicy::fixed_48
          ? 1.0 / 48.0
          : 1.0 / static_cast<double>(definition.cells);
  const auto surface_fixture =
      make_manufactured_surface(definition.body, surface_h);
  std::optional<test::Stage3TemporaryDirectory> root_directory;
  std::string surface_path_text;
  if (mpi.rank() == 0) {
    root_directory.emplace("task11-manufactured");
    const auto root_path = root_directory->path() / "body.stl";
    test::write_text(root_path,
                     test::ascii_stl(surface_fixture.triangles, "body"));
    surface_path_text = root_path.string();
  }
  std::uint64_t surface_path_size = surface_path_text.size();
  check_mpi(MPI_Bcast(&surface_path_size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(manufactured surface path size)");
  if (surface_path_size >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::runtime_error("manufactured surface path exceeds MPI int");
  surface_path_text.resize(static_cast<std::size_t>(surface_path_size));
  check_mpi(MPI_Bcast(surface_path_text.data(),
                      static_cast<int>(surface_path_text.size()), MPI_BYTE, 0,
                      mpi.comm()),
            "MPI_Bcast(manufactured surface path)");
  const std::filesystem::path surface_path(surface_path_text);

  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {definition.cells, definition.cells, definition.cells},
      {false, false, false},
      runtime::DecompositionOptions{definition.process_grid});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry =
      definition.mapping == ManufacturedMapping::uniform
          ? mesh::MeshGeometry(
                topology,
                mesh::UniformBoxMapping{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}})
          : mesh::MeshGeometry(
                topology,
                mesh::AnalyticWarpedBoxMapping{
                    {0.0, 0.0, 0.0},
                    {1.0, 1.0, 1.0},
                    {kWarpAmplitude[0], kWarpAmplitude[1], kWarpAmplitude[2]}});
  auto boundaries = boundary::BoundaryRegistry::create(
      manufactured_flow_config(definition, mpi.size()), topology);
  const auto surface =
      immersed::ImmersedSurface::load_collective(surface_path, 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain =
      immersed::ImmersedDomain::create(surface, query, definition.fluid_side,
                                       topology, geometry, boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto wall_plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  if (ghost_plan.maximum_halo_reach() > 4U)
    throw std::runtime_error("manufactured stencil exceeds four halo layers");

  double h_max = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    if (domain.region(cell) == immersed::CellRegion::fluid)
      h_max = std::max(h_max, std::cbrt(geometry.cell_volume_m3(cell)));
  max_in_place(mpi, &h_max, 1U);
  if (!(h_max > 0.0) || !std::isfinite(h_max))
    throw std::runtime_error("manufactured active h_max is invalid");
  const double dt =
      0.05 * h_max * h_max / (kReferenceVelocityMPerS * kReferenceLengthM);
  const auto force_reference =
      definition.collect_force
          ? std::optional<AnalyticForceReference>(
                verified_force_reference(definition.body, dt))
          : std::nullopt;

  geometry.require_compatible(topology);
  const auto make_cell_averages = [&](double time_s) {
    std::vector<MmsCellAverage> averages(topology.owned_cell_count());
    const auto hardware_threads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::size_t worker_count = std::min<std::size_t>(
        std::min<std::size_t>(
            hardware_threads,
            std::max<std::size_t>(1U,
                                  kParallelWorkerBudget /
                                      static_cast<std::size_t>(mpi.size()))),
        topology.owned_cell_count());
    std::vector<std::exception_ptr> failures(worker_count);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    std::atomic<bool> cancel{false};
    try {
      for (std::size_t worker = 0U; worker < worker_count; ++worker) {
        workers.emplace_back([&, worker] {
          try {
            const auto begin =
                topology.owned_cell_count() * worker / worker_count;
            const auto end =
                topology.owned_cell_count() * (worker + 1U) / worker_count;
            for (mesh::LocalCellId cell = begin;
                 cell < end && !cancel.load(std::memory_order_relaxed); ++cell)
              if (domain.region(cell) == immersed::CellRegion::fluid)
                averages[cell] = detail::evaluate_cell_average_validated(
                    topology, geometry, cell, definition.body, time_s);
          } catch (...) {
            failures[worker] = std::current_exception();
            cancel.store(true, std::memory_order_relaxed);
          }
        });
      }
    } catch (...) {
      cancel.store(true, std::memory_order_relaxed);
      for (auto &worker : workers)
        worker.join();
      throw;
    }
    for (auto &worker : workers)
      worker.join();
    for (const auto &failure : failures)
      if (failure)
        std::rethrow_exception(failure);
    return averages;
  };

  const auto next_averages = make_cell_averages(dt);
  auto previous_averages = next_averages;
  auto current_averages = next_averages;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    previous_averages[cell] = rescale_average(next_averages[cell], dt, -dt);
    current_averages[cell] = rescale_average(next_averages[cell], dt, 0.0);
  }

  const auto make_face_histories = [&](double time_s) {
    std::vector<MmsFaceHistory> histories(topology.local_face_count());
    const auto hardware_threads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::size_t worker_count = std::min<std::size_t>(
        std::min<std::size_t>(
            hardware_threads,
            std::max<std::size_t>(1U,
                                  kParallelWorkerBudget /
                                      static_cast<std::size_t>(mpi.size()))),
        topology.local_face_count());
    std::vector<std::exception_ptr> failures(worker_count);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    std::atomic<bool> cancel{false};
    try {
      for (std::size_t worker = 0U; worker < worker_count; ++worker) {
        workers.emplace_back([&, worker] {
          try {
            const auto begin =
                topology.local_face_count() * worker / worker_count;
            const auto end =
                topology.local_face_count() * (worker + 1U) / worker_count;
            for (mesh::LocalFaceId face = begin;
                 face < end && !cancel.load(std::memory_order_relaxed);
                 ++face) {
              const bool owner_active = domain.region(topology.owner(face)) ==
                                        immersed::CellRegion::fluid;
              const auto neighbour = topology.neighbour(face);
              const bool neighbour_active =
                  neighbour.has_value() &&
                  domain.region(*neighbour) == immersed::CellRegion::fluid;
              if (owner_active && (!neighbour.has_value() || neighbour_active))
                histories[face] = detail::evaluate_face_history_validated(
                    topology, geometry, face, definition.body, time_s);
            }
          } catch (...) {
            failures[worker] = std::current_exception();
            cancel.store(true, std::memory_order_relaxed);
          }
        });
      }
    } catch (...) {
      cancel.store(true, std::memory_order_relaxed);
      for (auto &worker : workers)
        worker.join();
      throw;
    }
    for (auto &worker : workers)
      worker.join();
    for (const auto &failure : failures)
      if (failure)
        std::rethrow_exception(failure);
    return histories;
  };
  const auto previous_face_histories = make_face_histories(-dt);
  const auto current_face_histories = make_face_histories(0.0);
  const auto next_face_histories =
      definition.collect_exact_momentum_residual
          ? make_face_histories(dt)
          : std::vector<MmsFaceHistory>{};

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density =
      registry.declare_field(manufactured_cell_descriptor("rho", 1U));
  fields.velocity =
      registry.declare_field(manufactured_cell_descriptor("velocity", 3U));
  fields.mechanical_pressure =
      registry.declare_field(manufactured_cell_descriptor("pi", 1U));
  fields.face_velocity =
      registry.declare_field(manufactured_face_descriptor("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();

  const auto make_layer =
      [&](const std::vector<MmsCellAverage> &averages,
          const std::vector<MmsFaceHistory> &face_histories) {
        flow::FlowLayerValues values;
        values.density.resize(topology.owned_cell_count(), 0.0);
        values.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
        values.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
        values.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
        values.face_mass_flux.resize(topology.local_face_count(), 0.0);
        for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
             ++cell) {
          if (domain.region(cell) != immersed::CellRegion::fluid)
            continue;
          const auto &exact = averages[cell];
          values.density[cell] = kReferenceDensityKgPerM3;
          values.velocity[cell * 3U] =
              canonical_manufactured_zero(exact.velocity_m_per_s.x);
          values.velocity[cell * 3U + 1U] =
              canonical_manufactured_zero(exact.velocity_m_per_s.y);
          values.velocity[cell * 3U + 2U] =
              canonical_manufactured_zero(exact.velocity_m_per_s.z);
          values.mechanical_pressure[cell] =
              canonical_manufactured_zero(exact.mechanical_pressure_pa);
        }
        for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
             ++face) {
          const bool owner_active = domain.region(topology.owner(face)) ==
                                    immersed::CellRegion::fluid;
          const auto neighbour = topology.neighbour(face);
          const bool neighbour_active =
              neighbour.has_value() &&
              domain.region(*neighbour) == immersed::CellRegion::fluid;
          if (!owner_active || (neighbour.has_value() && !neighbour_active))
            continue;
          const auto &history = face_histories[face];
          values.face_velocity[face * 3U] =
              canonical_manufactured_zero(history.velocity_m_per_s.x);
          values.face_velocity[face * 3U + 1U] =
              canonical_manufactured_zero(history.velocity_m_per_s.y);
          values.face_velocity[face * 3U + 2U] =
              canonical_manufactured_zero(history.velocity_m_per_s.z);
          values.face_mass_flux[face] = canonical_manufactured_zero(
              kReferenceDensityKgPerM3 *
              history.owner_normal_volume_flux_m3_per_s);
        }
        return values;
      };

  auto state = flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields, {1U, 0.0, dt, dt, flow::MomentumTimeOrder::bdf2});
  state.seed_accepted_layers(
      make_layer(previous_averages, previous_face_histories),
      make_layer(current_averages, current_face_histories));
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 4));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  immersed::LocalFlowPatternTransform transform;
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      definition.collect_force ? &wall_plan : nullptr, &transform, nullptr, mpi,
      execution, halo, momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_pc);
  ManufacturedRunResult result;
  result.surface_maximum_edge_m = surface_fixture.maximum_edge_m;
  result.surface_maximum_chord_error_m = surface_fixture.maximum_chord_error_m;
  result.classification_fingerprint = domain.classification_fingerprint();
  result.surface_coverage_fingerprint = domain.surface_coverage_fingerprint();
  result.ghost_plan_fingerprint = ghost_plan.fingerprint();
  result.wall_plan_fingerprint = wall_plan.fingerprint();
  result.operator_structure_fingerprint =
      flow::test::ImmersedFlowTestAccess::immersed_operator_structure_fingerprint(
          immersed_flow);

  AnalyticForceReference wall_plan_reference;
  for (const auto &point : wall_plan.local_points())
    accumulate_force(wall_plan_reference, definition.body, point.position_m,
                     point.solid_to_fluid_normal, point.weight_m2, dt);
  std::array<double, 9> wall_plan_force{wall_plan_reference.force.pressure_N.x,
                                        wall_plan_reference.force.pressure_N.y,
                                        wall_plan_reference.force.pressure_N.z,
                                        wall_plan_reference.force.viscous_N.x,
                                        wall_plan_reference.force.viscous_N.y,
                                        wall_plan_reference.force.viscous_N.z,
                                        wall_plan_reference.force.total_N.x,
                                        wall_plan_reference.force.total_N.y,
                                        wall_plan_reference.force.total_N.z};
  sum_in_place(mpi, wall_plan_force.data(), wall_plan_force.size());
  result.wall_plan_analytic_force = {
      {wall_plan_force[0], wall_plan_force[1], wall_plan_force[2]},
      {wall_plan_force[3], wall_plan_force[4], wall_plan_force[5]},
      {wall_plan_force[6], wall_plan_force[7], wall_plan_force[8]}};

  std::vector<double> source(domain.active_cells().owned_active_count() * 3U,
                             0.0);
  const auto &active_ids = domain.active_cells().ordered_global_ids();
  for (std::size_t row = 0U; row < domain.active_cells().owned_active_count();
       ++row) {
    const auto local = topology.find_local_cell(active_ids[row]);
    if (!local.has_value())
      throw std::runtime_error("manufactured active source cell is absent");
    const auto value = next_averages[*local].body_source_N_per_m3;
    source[row * 3U] = value.x;
    source[row * 3U + 1U] = value.y;
    source[row * 3U + 2U] = value.z;
  }
  flow::test::ImmersedFlowTestAccess::set_manufactured_body_source(immersed_flow,
                                                                 source);
  if (definition.collect_exact_momentum_residual) {
    if (mpi.size() != 1)
      throw std::runtime_error(
          "manufactured exact momentum residual fast diagnostic is one-rank");
    std::array<std::vector<flow::test::ImmersedFlowWallGradientSnapshot>, 3>
        exact_wall_gradients;
    const std::array<double, 3> times{-dt, 0.0, dt};
    for (std::size_t time_index = 0U; time_index < times.size(); ++time_index) {
      auto &values = exact_wall_gradients[time_index];
      for (const auto &link : domain.links()) {
        const auto fluid = topology.find_local_cell(link.fluid_cell);
        if (!fluid.has_value() ||
            topology.cell_ownership(*fluid) != mesh::EntityOwnership::owned)
          continue;
        const auto sample = evaluate_mms(definition.body,
                                         link.wall_intercept_m,
                                         times[time_index]);
        values.push_back(
            {link.id,
             dot(sample.mechanical_pressure_gradient_pa_per_m,
                 link.solid_to_fluid_normal)});
      }
      std::sort(values.begin(), values.end(),
                [](const auto &left, const auto &right) {
                  return left.link < right.link;
                });
    }
    const auto exact_trial = make_layer(next_averages, next_face_histories);
    const auto report =
        flow::test::ImmersedFlowTestAccess::exact_momentum_residual_terms(
            immersed_flow, state, exact_trial,
            flow::make_momentum_time_stencil(flow::MomentumTimeOrder::bdf2,
                                             dt, dt),
            kReferenceDensityKgPerM3, kDynamicViscosityPaS, source,
            exact_wall_gradients);
    if (!flow::test::ImmersedFlowTestAccess::
             exact_momentum_residual_report_is_self_consistent(report))
      throw std::runtime_error(
          "manufactured exact momentum residual report is inconsistent");
    auto mutation = report;
    mutation.convective_residual_N.front() += 1.0;
    if (flow::test::ImmersedFlowTestAccess::
            exact_momentum_residual_report_is_self_consistent(mutation))
      throw std::runtime_error(
          "manufactured exact momentum residual oracle missed a mutation");
    auto pressure_split_mutation = report;
    pressure_split_mutation.background_pressure_residual_N.front() += 1.0;
    if (flow::test::ImmersedFlowTestAccess::
            exact_momentum_residual_report_is_self_consistent(
                pressure_split_mutation))
      throw std::runtime_error(
          "manufactured exact pressure split oracle missed a mutation");

    auto reconstructed_face_pressure =
        report.background_pressure_residual_N;
    auto analytic_face_pressure = report.background_pressure_residual_N;
    auto unconstrained_center_pressure =
        report.background_pressure_residual_N;
    auto analytic_center_pressure = report.background_pressure_residual_N;
    auto exact_shared_center_pressure =
        report.background_pressure_residual_N;
    std::vector<double> reconstructed_face_correction(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> analytic_face_correction(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> unconstrained_center_correction(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> analytic_center_correction(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> exact_shared_center_correction(
        report.background_pressure_residual_N.size(), 0.0);
    std::map<std::pair<mesh::GlobalCellId, mesh::GlobalCellId>,
             mesh::LocalFaceId>
        faces_by_cell_pair;
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value())
        continue;
      const auto owner = topology.global_cell_id(topology.owner(face));
      const auto other = topology.global_cell_id(*neighbour);
      const auto key = std::minmax(owner, other);
      if (!faces_by_cell_pair.emplace(
               std::pair<mesh::GlobalCellId, mesh::GlobalCellId>{key.first,
                                                                 key.second},
               face)
               .second)
        throw std::runtime_error(
            "manufactured pressure face identity is duplicated");
    }
    const auto operator_rows =
        finite_volume::test::ImmersedOperatorTestAccess::rows(
            flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow));
    double face_area_closure_linf = 0.0;
    for (const auto &row : operator_rows) {
      if (row.links.empty())
        continue;
      const auto report_row = std::lower_bound(
          report.active_global_cell_ids.begin(),
          report.active_global_cell_ids.end(), row.active_cell);
      if (report_row == report.active_global_cell_ids.end() ||
          *report_row != row.active_cell)
        throw std::runtime_error(
            "manufactured pressure face row is unavailable");
      const std::size_t row_index = static_cast<std::size_t>(
          report_row - report.active_global_cell_ids.begin());
      for (const auto &link : row.links) {
        const auto domain_link = std::lower_bound(
            domain.links().begin(), domain.links().end(), link.id,
            [](const auto &candidate, immersed::ImmersedLinkId id) {
              return candidate.id < id;
            });
        if (domain_link == domain.links().end() || domain_link->id != link.id)
          throw std::runtime_error(
              "manufactured pressure face link is unavailable");
        const auto key = std::minmax(domain_link->fluid_cell,
                                     domain_link->solid_cell);
        const auto face = faces_by_cell_pair.find(
            {key.first, key.second});
        if (face == faces_by_cell_pair.end())
          throw std::runtime_error(
              "manufactured immersed pressure face is unavailable");
        const auto logical_vertices =
            logical_face_vertices(topology.logical_face(face->second));
        std::array<runtime::Real3, 4> vertices{};
        for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex)
          vertices[vertex] =
              geometry.vertex_position_m(logical_vertices[vertex]);
        const auto &reconstruction = ghost_plan.reconstruction(link.id);
        const auto gradient = std::lower_bound(
            exact_wall_gradients[2].begin(), exact_wall_gradients[2].end(),
            link.id, [](const auto &candidate, immersed::ImmersedLinkId id) {
              return candidate.link < id;
            });
        if (gradient == exact_wall_gradients[2].end() ||
            gradient->link != link.id)
          throw std::runtime_error(
              "manufactured pressure face gradient is unavailable");
        const auto reconstructed_value = [&](runtime::Real3 point_m) {
          const auto functional =
              immersed::detail::QuadraticReconstructionWeights::
                  origin_normal_gradient_constrained_value_weights(
                      reconstruction, point_m);
          double value = functional.boundary_coefficient *
                         gradient->normal_gradient_pa_per_m;
          for (const auto &donor : functional.donors) {
            const auto local = topology.find_local_cell(donor.global_cell);
            if (!local.has_value())
              throw std::runtime_error(
                  "manufactured pressure face donor is unavailable");
            value += donor.weight *
                     next_averages[*local].mechanical_pressure_pa;
          }
          if (!std::isfinite(value))
            throw std::runtime_error(
                "manufactured reconstructed face pressure is non-finite");
          return value;
        };
        const auto unconstrained_value = [&](runtime::Real3 point_m) {
          const auto weights =
              immersed::detail::QuadraticReconstructionWeights::value_weights(
                  reconstruction, point_m);
          double value = 0.0;
          for (const auto &donor : weights) {
            const auto local = topology.find_local_cell(donor.global_cell);
            if (!local.has_value())
              throw std::runtime_error(
                  "manufactured unconstrained pressure donor is unavailable");
            value += donor.weight *
                     next_averages[*local].mechanical_pressure_pa;
          }
          if (!std::isfinite(value))
            throw std::runtime_error(
                "manufactured unconstrained face pressure is non-finite");
          return value;
        };
        const auto analytic_value = [&](runtime::Real3 point_m) {
          return evaluate_mms(definition.body, point_m, dt)
              .mechanical_pressure_pa;
        };
        const auto integrate_face = [&](const auto &rule,
                                        const auto &evaluate) {
          return add(integrate_triangle_pressure_area(
                         vertices[0], vertices[1], vertices[2], rule,
                         evaluate),
                     integrate_triangle_pressure_area(
                         vertices[0], vertices[2], vertices[3], rule,
                         evaluate));
        };
        auto oriented_area = integrate_face(
            kQuadraticTriangleRule, [](runtime::Real3) { return 1.0; });
        auto reconstructed_integral =
            integrate_face(kQuadraticTriangleRule, reconstructed_value);
        auto analytic_integral =
            integrate_face(kDegreeFiveTriangleRule, analytic_value);
        if (dot(oriented_area, link.area_from_fluid_m2) < 0.0) {
          oriented_area = multiply(-1.0, oriented_area);
          reconstructed_integral = multiply(-1.0, reconstructed_integral);
          analytic_integral = multiply(-1.0, analytic_integral);
        }
        const double area_closure =
            magnitude(subtract(oriented_area, link.area_from_fluid_m2));
        const double area_scale =
            std::max(1.0, magnitude(link.area_from_fluid_m2));
        if (area_closure >
            256.0 * std::numeric_limits<double>::epsilon() * area_scale)
          throw std::runtime_error(
              "manufactured pressure face area does not match geometry");
        face_area_closure_linf =
            std::max(face_area_closure_linf, area_closure);
        const auto centered = multiply(
            reconstructed_value(link.pressure_quadrature_m),
            link.area_from_fluid_m2);
        const auto reconstructed_correction =
            subtract(reconstructed_integral, centered);
        const auto analytic_correction =
            subtract(analytic_integral, centered);
        const auto unconstrained_centered =
            multiply(unconstrained_value(link.pressure_quadrature_m),
                     link.area_from_fluid_m2);
        const auto analytic_centered =
            multiply(analytic_value(link.pressure_quadrature_m),
                     link.area_from_fluid_m2);
        const auto unconstrained_correction =
            subtract(unconstrained_centered, centered);
        const auto analytic_center_correction_value =
            subtract(analytic_centered, centered);
        const std::array<double, 3> reconstructed_components{
            reconstructed_correction.x, reconstructed_correction.y,
            reconstructed_correction.z};
        const std::array<double, 3> analytic_components{
            analytic_correction.x, analytic_correction.y,
            analytic_correction.z};
        const std::array<double, 3> unconstrained_components{
            unconstrained_correction.x, unconstrained_correction.y,
            unconstrained_correction.z};
        const std::array<double, 3> analytic_center_components{
            analytic_center_correction_value.x,
            analytic_center_correction_value.y,
            analytic_center_correction_value.z};
        for (std::size_t component = 0U; component < 3U; ++component) {
          const auto offset = row_index * 3U + component;
          reconstructed_face_correction[offset] +=
              reconstructed_components[component];
          analytic_face_correction[offset] += analytic_components[component];
          unconstrained_center_correction[offset] +=
              unconstrained_components[component];
          analytic_center_correction[offset] +=
              analytic_center_components[component];
          reconstructed_face_pressure[offset] +=
              reconstructed_components[component];
          analytic_face_pressure[offset] += analytic_components[component];
          unconstrained_center_pressure[offset] +=
              unconstrained_components[component];
          analytic_center_pressure[offset] +=
              analytic_center_components[component];
        }
      }
    }
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      const auto neighbour = topology.neighbour(face);
      if (!neighbour.has_value() ||
          domain.region(topology.owner(face)) != immersed::CellRegion::fluid ||
          domain.region(*neighbour) != immersed::CellRegion::fluid)
        continue;
      const auto owner = topology.owner(face);
      const auto owner_center = geometry.cell_center_m(owner);
      const auto face_center = geometry.face_center_m(face);
      const auto displacement = geometry.face_displacement_m(face);
      const double distance = magnitude(displacement);
      const double a = magnitude(subtract(face_center, owner_center));
      const double b = distance - a;
      if (!(distance > 0.0) || !(a > 0.0) || !(b > 0.0) ||
          !std::isfinite(distance))
        throw std::runtime_error(
            "manufactured shared pressure interpolation is invalid");
      const double linear_pressure =
          (b / distance) * next_averages[owner].mechanical_pressure_pa +
          (a / distance) * next_averages[*neighbour].mechanical_pressure_pa;
      const double exact_pressure =
          evaluate_mms(definition.body, face_center, dt)
              .mechanical_pressure_pa;
      const auto area =
          geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
      const auto correction =
          multiply(exact_pressure - linear_pressure, area);
      if (!finite(correction))
        throw std::runtime_error(
            "manufactured shared pressure correction is non-finite");
      const auto add_to_row = [&](mesh::LocalCellId cell,
                                  runtime::Real3 value) {
        const auto global = topology.global_cell_id(cell);
        const auto position = std::lower_bound(
            report.active_global_cell_ids.begin(),
            report.active_global_cell_ids.end(), global);
        if (position == report.active_global_cell_ids.end() ||
            *position != global)
          throw std::runtime_error(
              "manufactured shared pressure row is unavailable");
        const auto row_index = static_cast<std::size_t>(
            position - report.active_global_cell_ids.begin());
        const std::array<double, 3> components{value.x, value.y, value.z};
        for (std::size_t component = 0U; component < 3U; ++component) {
          const auto offset = row_index * 3U + component;
          exact_shared_center_correction[offset] += components[component];
          exact_shared_center_pressure[offset] += components[component];
        }
      };
      if (topology.cell_ownership(owner) == mesh::EntityOwnership::owned)
        add_to_row(owner, correction);
      const bool neighbour_owned =
          topology.cell_ownership(*neighbour) == mesh::EntityOwnership::owned &&
          !topology.periodic_pair(face).has_value();
      if (neighbour_owned)
        add_to_row(*neighbour, multiply(-1.0, correction));
    }
    runtime::Real3 shared_center_correction_sum{};
    double shared_center_correction_scale = 0.0;
    for (std::size_t row = 0U;
         row < report.active_global_cell_ids.size(); ++row) {
      shared_center_correction_sum =
          add(shared_center_correction_sum,
              {exact_shared_center_correction[row * 3U],
               exact_shared_center_correction[row * 3U + 1U],
               exact_shared_center_correction[row * 3U + 2U]});
      for (std::size_t component = 0U; component < 3U; ++component)
        shared_center_correction_scale +=
            std::abs(exact_shared_center_correction[row * 3U + component]);
    }
    if (max_abs(shared_center_correction_sum) >
        512.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, shared_center_correction_scale))
      throw std::runtime_error(
          "manufactured shared pressure correction is not conservative");
    const auto correction_closure = [](const std::vector<double> &candidate,
                                       const std::vector<double> &background,
                                       const std::vector<double> &correction) {
      if (candidate.size() != background.size() ||
          candidate.size() != correction.size())
        return std::numeric_limits<double>::infinity();
      double linf = 0.0;
      for (std::size_t index = 0U; index < candidate.size(); ++index)
        linf = std::max(
            linf,
            std::abs(candidate[index] -
                     (background[index] + correction[index])));
      return linf;
    };
    const double reconstructed_correction_closure = correction_closure(
        reconstructed_face_pressure, report.background_pressure_residual_N,
        reconstructed_face_correction);
    const double analytic_correction_closure = correction_closure(
        analytic_face_pressure, report.background_pressure_residual_N,
        analytic_face_correction);
    const double unconstrained_center_correction_closure = correction_closure(
        unconstrained_center_pressure,
        report.background_pressure_residual_N,
        unconstrained_center_correction);
    const double analytic_center_correction_closure = correction_closure(
        analytic_center_pressure, report.background_pressure_residual_N,
        analytic_center_correction);
    const double exact_shared_center_correction_closure = correction_closure(
        exact_shared_center_pressure, report.background_pressure_residual_N,
        exact_shared_center_correction);
    const double face_correction_closure_linf =
        std::max({reconstructed_correction_closure,
                  analytic_correction_closure,
                  unconstrained_center_correction_closure,
                  analytic_center_correction_closure,
                  exact_shared_center_correction_closure});
    if (!std::isfinite(face_correction_closure_linf) ||
        face_correction_closure_linf >
            2048.0 * std::numeric_limits<double>::epsilon())
      throw std::runtime_error(
          "manufactured pressure face correction does not close");
    auto correction_mutation = reconstructed_face_pressure;
    correction_mutation.front() += 1.0;
    if (correction_closure(correction_mutation,
                           report.background_pressure_residual_N,
                           reconstructed_face_correction) <=
        2048.0 * std::numeric_limits<double>::epsilon())
      throw std::runtime_error(
          "manufactured pressure face correction oracle missed a mutation");

    struct IncidentFace final {
      mesh::LocalFaceId face{};
      double owner_orientation{};
    };
    std::vector<std::vector<IncidentFace>> incident_faces(
        topology.local_cell_count());
    for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
         ++face) {
      incident_faces[topology.owner(face)].push_back({face, 1.0});
      const auto neighbour = topology.neighbour(face);
      if (neighbour.has_value())
        incident_faces[*neighbour].push_back({face, -1.0});
    }
    std::vector<double> full_exact_center_pressure(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> full_exact_integral_pressure(
        report.background_pressure_residual_N.size(), 0.0);
    for (std::size_t row = 0U;
         row < report.active_global_cell_ids.size(); ++row) {
      const auto local =
          topology.find_local_cell(report.active_global_cell_ids[row]);
      if (!local.has_value() || incident_faces[*local].size() != 6U)
        throw std::runtime_error(
            "manufactured complete pressure row does not have six faces");
      runtime::Real3 center_row{};
      runtime::Real3 integral_row{};
      for (const auto &incident : incident_faces[*local]) {
        const auto owner_area = geometry.face_area_vector_m2(
            incident.face, mesh::FaceSide::owner);
        const auto outward_area =
            multiply(incident.owner_orientation, owner_area);
        const auto face_center = geometry.face_center_m(incident.face);
        center_row =
            add(center_row,
                multiply(evaluate_mms(definition.body, face_center, dt)
                             .mechanical_pressure_pa,
                         outward_area));
        const auto logical_vertices =
            logical_face_vertices(topology.logical_face(incident.face));
        std::array<runtime::Real3, 4> vertices{};
        for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex)
          vertices[vertex] =
              geometry.vertex_position_m(logical_vertices[vertex]);
        const auto analytic_value = [&](runtime::Real3 point_m) {
          return evaluate_mms(definition.body, point_m, dt)
              .mechanical_pressure_pa;
        };
        auto face_integral =
            add(integrate_triangle_pressure_area(
                    vertices[0], vertices[1], vertices[2],
                    kDegreeFiveTriangleRule, analytic_value),
                integrate_triangle_pressure_area(
                    vertices[0], vertices[2], vertices[3],
                    kDegreeFiveTriangleRule, analytic_value));
        const auto raw_area =
            add(integrate_triangle_pressure_area(
                    vertices[0], vertices[1], vertices[2],
                    kQuadraticTriangleRule,
                    [](runtime::Real3) { return 1.0; }),
                integrate_triangle_pressure_area(
                    vertices[0], vertices[2], vertices[3],
                    kQuadraticTriangleRule,
                    [](runtime::Real3) { return 1.0; }));
        if (dot(raw_area, owner_area) < 0.0)
          face_integral = multiply(-1.0, face_integral);
        integral_row =
            add(integral_row,
                multiply(incident.owner_orientation, face_integral));
      }
      const std::array<double, 3> center_components{
          center_row.x, center_row.y, center_row.z};
      const std::array<double, 3> integral_components{
          integral_row.x, integral_row.y, integral_row.z};
      for (std::size_t component = 0U; component < 3U; ++component) {
        full_exact_center_pressure[row * 3U + component] =
            center_components[component];
        full_exact_integral_pressure[row * 3U + component] =
            integral_components[component];
      }
    }
    std::vector<double> coherent_single_link_pressure(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> coherent_single_link_integral_pressure(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<double> unconstrained_single_link_integral_pressure(
        report.background_pressure_residual_N.size(), 0.0);
    std::vector<std::uint8_t> wall_link_class(
        report.active_global_cell_ids.size(), 0U);
    for (const auto &operator_row : operator_rows) {
      if (operator_row.links.empty())
        continue;
      const auto report_row = std::lower_bound(
          report.active_global_cell_ids.begin(),
          report.active_global_cell_ids.end(), operator_row.active_cell);
      if (report_row == report.active_global_cell_ids.end() ||
          *report_row != operator_row.active_cell)
        throw std::runtime_error(
            "manufactured coherent pressure row is unavailable");
      const auto row = static_cast<std::size_t>(
          report_row - report.active_global_cell_ids.begin());
      wall_link_class[row] = operator_row.links.size() == 1U ? 1U : 2U;
      if (operator_row.links.size() != 1U)
        continue;
      const auto local = topology.find_local_cell(operator_row.active_cell);
      if (!local.has_value() || incident_faces[*local].size() != 6U)
        throw std::runtime_error(
            "manufactured coherent pressure row does not have six faces");
      const auto link = operator_row.links.front().id;
      const auto &reconstruction = ghost_plan.reconstruction(link);
      const auto gradient = std::lower_bound(
          exact_wall_gradients[2].begin(), exact_wall_gradients[2].end(), link,
          [](const auto &candidate, immersed::ImmersedLinkId id) {
            return candidate.link < id;
          });
      if (gradient == exact_wall_gradients[2].end() ||
          gradient->link != link)
        throw std::runtime_error(
            "manufactured coherent pressure gradient is unavailable");
      const auto coherent_value = [&](runtime::Real3 point_m) {
        const auto functional = immersed::detail::QuadraticReconstructionWeights::
            origin_normal_gradient_constrained_value_weights(reconstruction,
                                                             point_m);
        double value = functional.boundary_coefficient *
                       gradient->normal_gradient_pa_per_m;
        for (const auto &donor : functional.donors) {
          const auto donor_local =
              topology.find_local_cell(donor.global_cell);
          if (!donor_local.has_value())
            throw std::runtime_error(
                "manufactured coherent pressure donor is unavailable");
          value += donor.weight *
                   next_averages[*donor_local].mechanical_pressure_pa;
        }
        if (!std::isfinite(value))
          throw std::runtime_error(
              "manufactured coherent pressure value is non-finite");
        return value;
      };
      const auto unconstrained_value = [&](runtime::Real3 point_m) {
        const auto weights =
            immersed::detail::QuadraticReconstructionWeights::value_weights(
                reconstruction, point_m);
        double value = 0.0;
        for (const auto &donor : weights) {
          const auto donor_local =
              topology.find_local_cell(donor.global_cell);
          if (!donor_local.has_value())
            throw std::runtime_error(
                "manufactured unconstrained row donor is unavailable");
          value += donor.weight *
                   next_averages[*donor_local].mechanical_pressure_pa;
        }
        if (!std::isfinite(value))
          throw std::runtime_error(
              "manufactured unconstrained row value is non-finite");
        return value;
      };
      runtime::Real3 coherent_row{};
      runtime::Real3 coherent_integral_row{};
      runtime::Real3 unconstrained_integral_row{};
      for (const auto &incident : incident_faces[*local]) {
        const auto owner_area = geometry.face_area_vector_m2(
            incident.face, mesh::FaceSide::owner);
        const auto outward_area =
            multiply(incident.owner_orientation, owner_area);
        coherent_row =
            add(coherent_row,
                multiply(coherent_value(
                             geometry.face_center_m(incident.face)),
                         outward_area));
        const auto logical_vertices =
            logical_face_vertices(topology.logical_face(incident.face));
        std::array<runtime::Real3, 4> vertices{};
        for (std::size_t vertex = 0U; vertex < vertices.size(); ++vertex)
          vertices[vertex] =
              geometry.vertex_position_m(logical_vertices[vertex]);
        auto face_integral =
            add(integrate_triangle_pressure_area(
                    vertices[0], vertices[1], vertices[2],
                    kQuadraticTriangleRule, coherent_value),
                integrate_triangle_pressure_area(
                    vertices[0], vertices[2], vertices[3],
                    kQuadraticTriangleRule, coherent_value));
        auto unconstrained_face_integral =
            add(integrate_triangle_pressure_area(
                    vertices[0], vertices[1], vertices[2],
                    kQuadraticTriangleRule, unconstrained_value),
                integrate_triangle_pressure_area(
                    vertices[0], vertices[2], vertices[3],
                    kQuadraticTriangleRule, unconstrained_value));
        const auto raw_area =
            add(integrate_triangle_pressure_area(
                    vertices[0], vertices[1], vertices[2],
                    kQuadraticTriangleRule,
                    [](runtime::Real3) { return 1.0; }),
                integrate_triangle_pressure_area(
                    vertices[0], vertices[2], vertices[3],
                    kQuadraticTriangleRule,
                    [](runtime::Real3) { return 1.0; }));
        if (dot(raw_area, owner_area) < 0.0) {
          face_integral = multiply(-1.0, face_integral);
          unconstrained_face_integral =
              multiply(-1.0, unconstrained_face_integral);
        }
        coherent_integral_row =
            add(coherent_integral_row,
                multiply(incident.owner_orientation, face_integral));
        unconstrained_integral_row =
            add(unconstrained_integral_row,
                multiply(incident.owner_orientation,
                         unconstrained_face_integral));
      }
      const std::array<double, 3> components{
          coherent_row.x, coherent_row.y, coherent_row.z};
      const std::array<double, 3> integral_components{
          coherent_integral_row.x, coherent_integral_row.y,
          coherent_integral_row.z};
      const std::array<double, 3> unconstrained_integral_components{
          unconstrained_integral_row.x, unconstrained_integral_row.y,
          unconstrained_integral_row.z};
      for (std::size_t component = 0U; component < 3U; ++component) {
        coherent_single_link_pressure[row * 3U + component] =
            components[component];
        coherent_single_link_integral_pressure[row * 3U + component] =
            integral_components[component];
        unconstrained_single_link_integral_pressure[row * 3U + component] =
            unconstrained_integral_components[component];
      }
    }
    const auto vector_difference_linf = [](const std::vector<double> &left,
                                           const std::vector<double> &right) {
      if (left.size() != right.size())
        return std::numeric_limits<double>::infinity();
      double result = 0.0;
      for (std::size_t index = 0U; index < left.size(); ++index)
        result = std::max(result, std::abs(left[index] - right[index]));
      return result;
    };
    if (vector_difference_linf(full_exact_integral_pressure,
                               full_exact_integral_pressure) != 0.0)
      throw std::runtime_error(
          "manufactured complete pressure exact-copy oracle failed");
    auto full_row_mutation = full_exact_integral_pressure;
    full_row_mutation.front() += 1.0;
    if (!(vector_difference_linf(full_row_mutation,
                                 full_exact_integral_pressure) > 0.0))
      throw std::runtime_error(
          "manufactured complete pressure oracle missed a mutation");
    if (vector_difference_linf(coherent_single_link_pressure,
                               coherent_single_link_pressure) != 0.0)
      throw std::runtime_error(
          "manufactured coherent pressure exact-copy oracle failed");
    const auto coherent_row = std::find(wall_link_class.begin(),
                                        wall_link_class.end(), 1U);
    if (coherent_row == wall_link_class.end())
      throw std::runtime_error(
          "manufactured coherent pressure row class is empty");
    auto coherent_mutation = coherent_single_link_pressure;
    coherent_mutation[static_cast<std::size_t>(
                          coherent_row - wall_link_class.begin()) *
                      3U] += 1.0;
    if (!(vector_difference_linf(coherent_mutation,
                                 coherent_single_link_pressure) > 0.0))
      throw std::runtime_error(
          "manufactured coherent pressure oracle missed a value mutation");
    if (vector_difference_linf(coherent_single_link_integral_pressure,
                               coherent_single_link_integral_pressure) != 0.0)
      throw std::runtime_error(
          "manufactured coherent pressure integral exact-copy oracle failed");
    auto integral_mutation = coherent_single_link_integral_pressure;
    integral_mutation[static_cast<std::size_t>(
                          coherent_row - wall_link_class.begin()) *
                      3U] += 1.0;
    if (!(vector_difference_linf(
              integral_mutation, coherent_single_link_integral_pressure) >
          0.0))
      throw std::runtime_error(
          "manufactured coherent pressure integral oracle missed a mutation");
    if (vector_difference_linf(unconstrained_single_link_integral_pressure,
                               unconstrained_single_link_integral_pressure) !=
        0.0)
      throw std::runtime_error(
          "manufactured unconstrained pressure exact-copy oracle failed");
    auto unconstrained_mutation = unconstrained_single_link_integral_pressure;
    unconstrained_mutation[static_cast<std::size_t>(
                               coherent_row - wall_link_class.begin()) *
                           3U] += 1.0;
    if (!(vector_difference_linf(
              unconstrained_mutation,
              unconstrained_single_link_integral_pressure) > 0.0))
      throw std::runtime_error(
          "manufactured unconstrained pressure oracle missed a mutation");
    const auto classification_is_consistent =
        [&](const std::vector<std::uint8_t> &candidate) {
          if (candidate.size() != report.active_global_cell_ids.size())
            return false;
          for (const auto &operator_row : operator_rows) {
            const auto report_row = std::lower_bound(
                report.active_global_cell_ids.begin(),
                report.active_global_cell_ids.end(), operator_row.active_cell);
            if (report_row == report.active_global_cell_ids.end() ||
                *report_row != operator_row.active_cell)
              return false;
            const auto row = static_cast<std::size_t>(
                report_row - report.active_global_cell_ids.begin());
            const std::uint8_t expected = operator_row.links.empty()
                                              ? 0U
                                              : operator_row.links.size() == 1U
                                                    ? 1U
                                                    : 2U;
            if (candidate[row] != expected)
              return false;
          }
          return true;
        };
    if (!classification_is_consistent(wall_link_class))
      throw std::runtime_error(
          "manufactured coherent pressure row classification is invalid");
    auto class_mutation = wall_link_class;
    class_mutation[static_cast<std::size_t>(coherent_row -
                                            wall_link_class.begin())] = 2U;
    if (classification_is_consistent(class_mutation))
      throw std::runtime_error(
          "manufactured coherent pressure oracle missed a class mutation");

    std::array<const std::vector<double> *, 7> terms{
        &report.time_residual_N,
        &report.convective_residual_N,
        &report.viscous_remainder_residual_N,
        &report.pressure_residual_N,
        &report.implicit_viscous_reference_residual_N,
        &report.source_residual_N,
        &report.total_residual_N};
    std::array<double, 7> interface_squares{};
    std::array<double, 7> bulk_squares{};
    double interface_pressure_balance_square = 0.0;
    double bulk_pressure_balance_square = 0.0;
    double interface_background_pressure_balance_square = 0.0;
    double bulk_background_pressure_balance_square = 0.0;
    double interface_reconstructed_face_pressure_balance_square = 0.0;
    double bulk_reconstructed_face_pressure_balance_square = 0.0;
    double interface_analytic_face_pressure_balance_square = 0.0;
    double bulk_analytic_face_pressure_balance_square = 0.0;
    double interface_reconstructed_to_analytic_face_square = 0.0;
    double interface_unconstrained_center_pressure_balance_square = 0.0;
    double bulk_unconstrained_center_pressure_balance_square = 0.0;
    double interface_analytic_center_pressure_balance_square = 0.0;
    double bulk_analytic_center_pressure_balance_square = 0.0;
    double interface_constrained_to_unconstrained_center_square = 0.0;
    double interface_exact_shared_center_pressure_balance_square = 0.0;
    double bulk_exact_shared_center_pressure_balance_square = 0.0;
    double interface_shared_center_correction_square = 0.0;
    double interface_full_exact_center_pressure_balance_square = 0.0;
    double bulk_full_exact_center_pressure_balance_square = 0.0;
    double interface_full_exact_integral_pressure_balance_square = 0.0;
    double bulk_full_exact_integral_pressure_balance_square = 0.0;
    double interface_background_to_full_exact_integral_square = 0.0;
    double single_link_current_background_pressure_balance_square = 0.0;
    double single_link_full_exact_center_pressure_balance_square = 0.0;
    double single_link_coherent_pressure_balance_square = 0.0;
    double single_link_coherent_integral_pressure_balance_square = 0.0;
    double single_link_unconstrained_integral_pressure_balance_square = 0.0;
    double single_link_constraint_integral_difference_square = 0.0;
    double single_link_volume = 0.0;
    double interface_pressure_wall_defect_square = 0.0;
    double bulk_pressure_wall_defect_square = 0.0;
    std::array<double, 4> interface_pressure_balance_linf{};
    double interface_volume = 0.0;
    double bulk_volume = 0.0;
    for (std::size_t row = 0U;
         row < report.active_global_cell_ids.size(); ++row) {
      const double volume = report.cell_volume_m3[row];
      auto &squares = report.immersed_interface_row[row] != 0U
                          ? interface_squares
                          : bulk_squares;
      auto &volume_sum = report.immersed_interface_row[row] != 0U
                             ? interface_volume
                             : bulk_volume;
      auto &pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_pressure_balance_square
              : bulk_pressure_balance_square;
      auto &background_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_background_pressure_balance_square
              : bulk_background_pressure_balance_square;
      auto &pressure_wall_defect_square =
          report.immersed_interface_row[row] != 0U
              ? interface_pressure_wall_defect_square
              : bulk_pressure_wall_defect_square;
      auto &reconstructed_face_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_reconstructed_face_pressure_balance_square
              : bulk_reconstructed_face_pressure_balance_square;
      auto &analytic_face_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_analytic_face_pressure_balance_square
              : bulk_analytic_face_pressure_balance_square;
      auto &unconstrained_center_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_unconstrained_center_pressure_balance_square
              : bulk_unconstrained_center_pressure_balance_square;
      auto &analytic_center_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_analytic_center_pressure_balance_square
              : bulk_analytic_center_pressure_balance_square;
      auto &exact_shared_center_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_exact_shared_center_pressure_balance_square
              : bulk_exact_shared_center_pressure_balance_square;
      auto &full_exact_center_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_full_exact_center_pressure_balance_square
              : bulk_full_exact_center_pressure_balance_square;
      auto &full_exact_integral_pressure_balance_square =
          report.immersed_interface_row[row] != 0U
              ? interface_full_exact_integral_pressure_balance_square
              : bulk_full_exact_integral_pressure_balance_square;
      const auto local =
          topology.find_local_cell(report.active_global_cell_ids[row]);
      if (!local.has_value())
        throw std::runtime_error(
            "manufactured exact pressure balance cell is absent");
      const auto &pressure_gradient =
          next_averages[*local].mechanical_pressure_gradient_pa_per_m;
      const std::array<double, 3> pressure_gradient_components{
          pressure_gradient.x, pressure_gradient.y, pressure_gradient.z};
      volume_sum += volume;
      if (wall_link_class[row] == 1U) {
        single_link_volume += volume;
        ++result.exact_momentum_residual.single_link_row_count;
      } else if (wall_link_class[row] == 2U) {
        ++result.exact_momentum_residual.multi_link_row_count;
      }
      if (report.immersed_interface_row[row] != 0U)
        ++result.exact_momentum_residual.interface_row_count;
      else
        ++result.exact_momentum_residual.bulk_row_count;
      for (std::size_t term_index = 0U; term_index < terms.size();
           ++term_index)
        for (std::size_t component_index = 0U; component_index < 3U;
             ++component_index) {
          const double value =
              (*terms[term_index])[row * 3U + component_index];
          squares[term_index] += value * value / volume;
        }
      for (std::size_t component_index = 0U; component_index < 3U;
           ++component_index) {
        const double defect =
            report.pressure_residual_N[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double background_defect =
            report.background_pressure_residual_N[row * 3U +
                                                  component_index] -
            volume * pressure_gradient_components[component_index];
        const double wall_defect =
            report.pressure_wall_defect_residual_N[row * 3U +
                                                   component_index];
        const double reconstructed_face_defect =
            reconstructed_face_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double analytic_face_defect =
            analytic_face_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double unconstrained_center_defect =
            unconstrained_center_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double analytic_center_defect =
            analytic_center_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double exact_shared_center_defect =
            exact_shared_center_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double full_exact_center_defect =
            full_exact_center_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        const double full_exact_integral_defect =
            full_exact_integral_pressure[row * 3U + component_index] -
            volume * pressure_gradient_components[component_index];
        if (report.immersed_interface_row[row] != 0U) {
          interface_pressure_balance_linf[0] = std::max(
              interface_pressure_balance_linf[0], std::abs(defect) / volume);
          interface_pressure_balance_linf[1] =
              std::max(interface_pressure_balance_linf[1],
                       std::abs(background_defect) / volume);
          interface_pressure_balance_linf[2] =
              std::max(interface_pressure_balance_linf[2],
                       std::abs(full_exact_center_defect) / volume);
          interface_pressure_balance_linf[3] =
              std::max(interface_pressure_balance_linf[3],
                       std::abs(wall_defect) / volume);
        }
        pressure_balance_square += defect * defect / volume;
        background_pressure_balance_square +=
            background_defect * background_defect / volume;
        pressure_wall_defect_square += wall_defect * wall_defect / volume;
        reconstructed_face_pressure_balance_square +=
            reconstructed_face_defect * reconstructed_face_defect / volume;
        analytic_face_pressure_balance_square +=
            analytic_face_defect * analytic_face_defect / volume;
        unconstrained_center_pressure_balance_square +=
            unconstrained_center_defect * unconstrained_center_defect / volume;
        analytic_center_pressure_balance_square +=
            analytic_center_defect * analytic_center_defect / volume;
        exact_shared_center_pressure_balance_square +=
            exact_shared_center_defect * exact_shared_center_defect / volume;
        full_exact_center_pressure_balance_square +=
            full_exact_center_defect * full_exact_center_defect / volume;
        full_exact_integral_pressure_balance_square +=
            full_exact_integral_defect * full_exact_integral_defect / volume;
        if (wall_link_class[row] == 1U) {
          const double coherent_defect =
              coherent_single_link_pressure[row * 3U + component_index] -
              volume * pressure_gradient_components[component_index];
          const double coherent_integral_defect =
              coherent_single_link_integral_pressure
                  [row * 3U + component_index] -
              volume * pressure_gradient_components[component_index];
          const double unconstrained_integral_defect =
              unconstrained_single_link_integral_pressure
                  [row * 3U + component_index] -
              volume * pressure_gradient_components[component_index];
          const double constraint_integral_difference =
              coherent_single_link_integral_pressure
                  [row * 3U + component_index] -
              unconstrained_single_link_integral_pressure
                  [row * 3U + component_index];
          single_link_current_background_pressure_balance_square +=
              background_defect * background_defect / volume;
          single_link_full_exact_center_pressure_balance_square +=
              full_exact_center_defect * full_exact_center_defect / volume;
          single_link_coherent_pressure_balance_square +=
              coherent_defect * coherent_defect / volume;
          single_link_coherent_integral_pressure_balance_square +=
              coherent_integral_defect * coherent_integral_defect / volume;
          single_link_unconstrained_integral_pressure_balance_square +=
              unconstrained_integral_defect * unconstrained_integral_defect /
              volume;
          single_link_constraint_integral_difference_square +=
              constraint_integral_difference * constraint_integral_difference /
              volume;
        }
        if (report.immersed_interface_row[row] != 0U) {
          const double reconstruction_difference =
              reconstructed_face_pressure[row * 3U + component_index] -
              analytic_face_pressure[row * 3U + component_index];
          interface_reconstructed_to_analytic_face_square +=
              reconstruction_difference * reconstruction_difference / volume;
          const double constraint_difference =
              report.background_pressure_residual_N
                  [row * 3U + component_index] -
              unconstrained_center_pressure[row * 3U + component_index];
          interface_constrained_to_unconstrained_center_square +=
              constraint_difference * constraint_difference / volume;
          const double shared_correction =
              exact_shared_center_correction[row * 3U + component_index];
          interface_shared_center_correction_square +=
              shared_correction * shared_correction / volume;
          const double background_difference =
              report.background_pressure_residual_N
                  [row * 3U + component_index] -
              full_exact_integral_pressure[row * 3U + component_index];
          interface_background_to_full_exact_integral_square +=
              background_difference * background_difference / volume;
        }
      }
    }
    if (!(interface_volume > 0.0) || !(bulk_volume > 0.0))
      throw std::runtime_error(
          "manufactured exact momentum residual row classes are empty");
    if (!(single_link_volume > 0.0) ||
        result.exact_momentum_residual.multi_link_row_count == 0U)
      throw std::runtime_error(
          "manufactured coherent pressure link classes are empty");
    max_in_place(mpi, interface_pressure_balance_linf.data(),
                 interface_pressure_balance_linf.size());
    for (std::size_t term_index = 0U; term_index < terms.size();
         ++term_index) {
      result.exact_momentum_residual.interface_rms_N_per_m3[term_index] =
          std::sqrt(interface_squares[term_index] /
                    (3.0 * interface_volume));
      result.exact_momentum_residual.bulk_rms_N_per_m3[term_index] =
          std::sqrt(bulk_squares[term_index] / (3.0 * bulk_volume));
    }
    result.exact_momentum_residual.interface_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual.bulk_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_pressure_balance_square / (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_background_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_background_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_background_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_background_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_reconstructed_face_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_reconstructed_face_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_reconstructed_face_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_reconstructed_face_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_analytic_face_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_analytic_face_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_analytic_face_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_analytic_face_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_reconstructed_to_analytic_face_rms_N_per_m3 =
        std::sqrt(interface_reconstructed_to_analytic_face_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .interface_unconstrained_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_unconstrained_center_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_unconstrained_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_unconstrained_center_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_analytic_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_analytic_center_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_analytic_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_analytic_center_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_constrained_to_unconstrained_center_rms_N_per_m3 =
        std::sqrt(interface_constrained_to_unconstrained_center_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .interface_exact_shared_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_exact_shared_center_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_exact_shared_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_exact_shared_center_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_shared_center_correction_rms_N_per_m3 =
        std::sqrt(interface_shared_center_correction_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .shared_center_correction_conservation_linf_N =
        max_abs(shared_center_correction_sum);
    result.exact_momentum_residual
        .interface_full_exact_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_full_exact_center_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_full_exact_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_full_exact_center_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_full_exact_integral_pressure_balance_rms_N_per_m3 =
        std::sqrt(interface_full_exact_integral_pressure_balance_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .bulk_full_exact_integral_pressure_balance_rms_N_per_m3 =
        std::sqrt(bulk_full_exact_integral_pressure_balance_square /
                  (3.0 * bulk_volume));
    result.exact_momentum_residual
        .interface_background_to_full_exact_integral_rms_N_per_m3 =
        std::sqrt(interface_background_to_full_exact_integral_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual
        .single_link_current_background_pressure_balance_rms_N_per_m3 =
        std::sqrt(single_link_current_background_pressure_balance_square /
                  (3.0 * single_link_volume));
    result.exact_momentum_residual
        .single_link_full_exact_center_pressure_balance_rms_N_per_m3 =
        std::sqrt(single_link_full_exact_center_pressure_balance_square /
                  (3.0 * single_link_volume));
    result.exact_momentum_residual
        .single_link_coherent_pressure_balance_rms_N_per_m3 =
        std::sqrt(single_link_coherent_pressure_balance_square /
                  (3.0 * single_link_volume));
    result.exact_momentum_residual
        .single_link_coherent_integral_pressure_balance_rms_N_per_m3 =
        std::sqrt(single_link_coherent_integral_pressure_balance_square /
                  (3.0 * single_link_volume));
    result.exact_momentum_residual
        .single_link_unconstrained_integral_pressure_balance_rms_N_per_m3 =
        std::sqrt(single_link_unconstrained_integral_pressure_balance_square /
                  (3.0 * single_link_volume));
    result.exact_momentum_residual
        .single_link_constraint_integral_difference_rms_N_per_m3 =
        std::sqrt(single_link_constraint_integral_difference_square /
                  (3.0 * single_link_volume));
    result.exact_momentum_residual
        .interface_pressure_wall_defect_rms_N_per_m3 =
        std::sqrt(interface_pressure_wall_defect_square /
                  (3.0 * interface_volume));
    result.exact_momentum_residual.bulk_pressure_wall_defect_rms_N_per_m3 =
        std::sqrt(bulk_pressure_wall_defect_square / (3.0 * bulk_volume));
    result.exact_momentum_residual.interface_pressure_balance_linf_N_per_m3 =
        interface_pressure_balance_linf;
    result.exact_momentum_residual.face_area_closure_linf_m2 =
        face_area_closure_linf;
    result.exact_momentum_residual.face_correction_closure_linf_N =
        face_correction_closure_linf;
    result.exact_momentum_residual.pointwise_split_closure_linf_N =
        report.pointwise_split_closure_linf_N;
    result.exact_momentum_residual.available = true;
  }
  const auto rollback_reference =
      definition.final_force_failure_rank.has_value()
          ? std::optional<Stage3StateSnapshot>(snapshot_stage3_state(state))
          : std::nullopt;
  if (definition.final_force_failure_rank.has_value())
    flow::test::ImmersedFlowTestAccess::set_failure_stage(
        flow::test::ImmersedFlowAttemptFailureStage::final_force_reconstruction,
        *definition.final_force_failure_rank);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    kReferenceDensityKgPerM3,
                                    kDynamicViscosityPaS,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  const auto stencil =
      flow::make_momentum_time_stencil(flow::MomentumTimeOrder::bdf2, dt, dt);
  linear::SolveControl momentum_control;
  momentum_control.atol = kManufacturedSolveAtol;
  momentum_control.rtol = kManufacturedSolveRtol;
  momentum_control.max_iterations = kManufacturedSolveMaxIterations;
  momentum_control.residual_recompute_interval =
      kManufacturedResidualRecomputeInterval;
  auto pressure_control = momentum_control;
  const auto attempt = [&] {
    try {
      auto value = immersed_flow.attempt(state, physics, stencil, momentum_control,
                                  pressure_control);
      flow::test::ImmersedFlowTestAccess::clear_failure_stage();
      flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(immersed_flow);
      return value;
    } catch (...) {
      flow::test::ImmersedFlowTestAccess::clear_failure_stage();
      flow::test::ImmersedFlowTestAccess::clear_manufactured_body_source(immersed_flow);
      throw;
    }
  }();
  const auto &base = std::get<flow::StepAttemptReport>(attempt.base);
  if (definition.final_force_failure_rank.has_value()) {
    result.committed = false;
    result.lowest_failing_rank = base.lowest_failing_rank;
    result.rollback_bitwise_equal =
        rollback_reference.has_value() &&
        stage3_state_bitwise_equal(*rollback_reference,
                                   snapshot_stage3_state(state));
    if (base.disposition != flow::StepAttemptDisposition::recoverable_failure ||
        base.reason != flow::StepFailureReason::final_conservation_defect ||
        base.pressure_corrector_count != 2U || attempt.force.has_value() ||
        base.lowest_failing_rank != *definition.final_force_failure_rank ||
        !result.rollback_bitwise_equal)
      throw std::runtime_error(
          "manufactured final-force failure contract is inconsistent");
    return result;
  }
  if (base.disposition != flow::StepAttemptDisposition::committed ||
      base.reason != flow::StepFailureReason::none ||
      base.pressure_corrector_count != 2U) {
    std::ostringstream message;
    message.precision(17);
    message << "manufactured immersed-flow attempt did not commit: disposition="
            << static_cast<int>(base.disposition)
            << " reason=" << static_cast<int>(base.reason)
            << " rank=" << base.lowest_failing_rank
            << " correctors=" << base.pressure_corrector_count
            << " continuity=" << base.final_continuity_normalized_l2
            << " pressure=" << base.final_pressure_residual_l2
            << " momentum=" << base.final_momentum_normalized_l2[0] << ','
            << base.final_momentum_normalized_l2[1] << ','
            << base.final_momentum_normalized_l2[2]
            << " p_solve0_reason="
            << static_cast<int>(base.pressure[0].reason)
            << " p_solve0_iterations=" << base.pressure[0].iterations
            << " p_solve0_final=" << base.pressure[0].final_residual
            << " p_solve1_reason="
            << static_cast<int>(base.pressure[1].reason)
            << " p_solve1_iterations=" << base.pressure[1].iterations
            << " p_solve1_final=" << base.pressure[1].final_residual
            << " p_solve1_init=" << base.pressure[1].initial_residual
            << " p_solve1_lowest=" << base.pressure[1].lowest_failing_rank
            << " p_solve1_matvec=" << base.pressure[1].matvec_count;
    throw std::runtime_error(message.str());
  }
  if (attempt.force.has_value() != definition.collect_force)
    throw std::runtime_error("manufactured force provenance is inconsistent");
  const auto committed = state.snapshot(flow::FlowLayer::committed);
  std::array<double, 8> sums{};
  double local_exact_velocity_linf = 0.0;
  double local_exact_pressure_linf = 0.0;
  double local_velocity_linf = 0.0;
  double local_pressure_linf = 0.0;
  double local_pressure_error_integral = 0.0;
  double local_active_volume = 0.0;
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const double volume = geometry.cell_volume_m3(cell);
    local_pressure_error_integral +=
        (committed.mechanical_pressure[cell] -
         next_averages[cell].mechanical_pressure_pa) *
        volume;
    local_active_volume += volume;
  }
  std::array<double, 2> gauge_values{local_pressure_error_integral,
                                     local_active_volume};
  sum_in_place(mpi, gauge_values.data(), gauge_values.size());
  if (!(gauge_values[1] > 0.0))
    throw std::runtime_error("manufactured active volume is empty");
  const double pressure_error_mean = gauge_values[0] / gauge_values[1];

  double local_numerical_pressure_integral = 0.0;
  NearWallPressureDiagnostics near_wall_pressure;
  std::map<mesh::GlobalCellId, std::uint64_t> incident_wall_link_counts;
  for (const auto &link : domain.links())
    ++incident_wall_link_counts[link.fluid_cell];
  std::vector<PressureErrorExtremum> local_pressure_extrema;
  local_pressure_extrema.reserve(domain.active_cells().owned_active_count());
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto &exact = next_averages[cell];
    const double volume = geometry.cell_volume_m3(cell);
    const double wall_distance = std::sqrt(
        query.closest_point(geometry.cell_center_m(cell)).squared_distance_m2);
    const double pressure_error = committed.mechanical_pressure[cell] -
                                  exact.mechanical_pressure_pa -
                                  pressure_error_mean;
    const auto global_id = topology.global_cell_id(cell);
    local_pressure_extrema.push_back(
        {pressure_error, std::abs(pressure_error), global_id,
         logical_cell_from_global_id(global_id, topology.global_extent()),
         wall_distance});
    const auto distance_band =
        normalized_near_wall_pressure_band(wall_distance / h_max);
    if (distance_band.has_value()) {
      const bool incident =
          incident_wall_link_counts.find(global_id) !=
          incident_wall_link_counts.end();
      auto &class_total = incident ? near_wall_pressure.incident_total
                                   : near_wall_pressure.nonincident_total;
      auto &class_band =
          incident ? near_wall_pressure.incident_distance_bands[*distance_band]
                   : near_wall_pressure
                         .nonincident_distance_bands[*distance_band];
      if (!accumulate_pressure_error_moments(near_wall_pressure.total,
                                             pressure_error, volume) ||
          !accumulate_pressure_error_moments(
              near_wall_pressure.distance_bands[*distance_band],
              pressure_error, volume) ||
          !accumulate_pressure_error_moments(class_total, pressure_error,
                                             volume) ||
          !accumulate_pressure_error_moments(class_band, pressure_error,
                                             volume))
        throw std::runtime_error(
            "manufactured near-wall pressure diagnostics are invalid");
    }
    sums[3] += pressure_error * pressure_error * volume;
    if (manufactured_near_wall_band_contains(wall_distance)) {
      sums[4] += pressure_error * pressure_error * volume;
      sums[5] += volume;
    }
    local_pressure_linf =
        std::max(local_pressure_linf, std::abs(pressure_error));
    local_exact_pressure_linf =
        std::max(local_exact_pressure_linf,
                 std::abs(exact.mechanical_pressure_pa + pressure_error_mean));
    local_numerical_pressure_integral +=
        committed.mechanical_pressure[cell] * volume;
    const std::array<double, 3> expected{exact.velocity_m_per_s.x,
                                         exact.velocity_m_per_s.y,
                                         exact.velocity_m_per_s.z};
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double numerical = committed.velocity[cell * 3U + component];
      const double error = numerical - expected[component];
      sums[0] += error * error * volume;
      sums[1] += volume;
      if (manufactured_near_wall_band_contains(wall_distance)) {
        sums[2] += error * error * volume;
        sums[6] += volume;
      }
      local_velocity_linf = std::max(local_velocity_linf, std::abs(error));
      local_exact_velocity_linf =
          std::max(local_exact_velocity_linf, std::abs(expected[component]));
    }
  }
  const auto pressure_error_extremum = reduce_pressure_error_extremum(
      mpi, local_pressure_extrema, topology.global_extent());
  sums[7] = local_numerical_pressure_integral;
  sum_in_place(mpi, sums.data(), sums.size());
  {
    std::array<PressureErrorMoments *, 15> moment_sets{
        &near_wall_pressure.distance_bands[0],
        &near_wall_pressure.distance_bands[1],
        &near_wall_pressure.distance_bands[2],
        &near_wall_pressure.distance_bands[3],
        &near_wall_pressure.incident_distance_bands[0],
        &near_wall_pressure.incident_distance_bands[1],
        &near_wall_pressure.incident_distance_bands[2],
        &near_wall_pressure.incident_distance_bands[3],
        &near_wall_pressure.nonincident_distance_bands[0],
        &near_wall_pressure.nonincident_distance_bands[1],
        &near_wall_pressure.nonincident_distance_bands[2],
        &near_wall_pressure.nonincident_distance_bands[3],
        &near_wall_pressure.total,
        &near_wall_pressure.incident_total,
        &near_wall_pressure.nonincident_total};
    std::array<double, 45> real_values{};
    std::array<std::uint64_t, 15> counts{};
    for (std::size_t index = 0U; index < moment_sets.size(); ++index) {
      real_values[index * 3U] =
          moment_sets[index]->signed_error_volume_pa_m3;
      real_values[index * 3U + 1U] =
          moment_sets[index]->squared_error_volume_pa2_m3;
      real_values[index * 3U + 2U] = moment_sets[index]->volume_m3;
      counts[index] = moment_sets[index]->cell_count;
    }
    sum_in_place(mpi, real_values.data(), real_values.size());
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, counts.data(),
                            static_cast<int>(counts.size()), MPI_UINT64_T,
                            MPI_SUM, mpi.comm()),
              "MPI_Allreduce(manufactured near-wall pressure counts)");
    for (std::size_t index = 0U; index < moment_sets.size(); ++index) {
      moment_sets[index]->signed_error_volume_pa_m3 =
          real_values[index * 3U];
      moment_sets[index]->squared_error_volume_pa2_m3 =
          real_values[index * 3U + 1U];
      moment_sets[index]->volume_m3 = real_values[index * 3U + 2U];
      moment_sets[index]->cell_count = counts[index];
    }
    near_wall_pressure.available = true;
  }
  std::array<double, 4> maxima{local_velocity_linf, local_pressure_linf,
                               local_exact_velocity_linf,
                               local_exact_pressure_linf};
  max_in_place(mpi, maxima.data(), maxima.size());
  if (!(sums[1] > 0.0) || !(sums[5] > 0.0) || !(sums[6] > 0.0))
    throw std::runtime_error("manufactured error support is empty");

  double penetration_square = 0.0;
  double penetration_weight = 0.0;
  double penetration_linf = 0.0;
  const auto velocity_view =
      state.layer(flow::FlowLayer::committed).view<double>(fields.velocity);
  for (const auto &point : wall_plan.local_points()) {
    runtime::Real3 wall_velocity{};
    wall_velocity.x = immersed::detail::value_with_origin_constraint(
        point.reconstruction, point.position_m, velocity_view, 0U, 0.0);
    wall_velocity.y = immersed::detail::value_with_origin_constraint(
        point.reconstruction, point.position_m, velocity_view, 1U, 0.0);
    wall_velocity.z = immersed::detail::value_with_origin_constraint(
        point.reconstruction, point.position_m, velocity_view, 2U, 0.0);
    const double normal_velocity =
        dot(wall_velocity, point.solid_to_fluid_normal);
    penetration_square += point.weight_m2 * normal_velocity * normal_velocity;
    penetration_weight += point.weight_m2;
    penetration_linf = std::max(penetration_linf, std::abs(normal_velocity));
  }
  std::array<double, 2> penetration_sums{penetration_square,
                                         penetration_weight};
  sum_in_place(mpi, penetration_sums.data(), penetration_sums.size());
  max_in_place(mpi, &penetration_linf, 1U);
  if (!(penetration_sums[1] > 0.0))
    throw std::runtime_error("manufactured wall quadrature is empty");

  result.errors.velocity_l2 = std::sqrt(sums[0] / sums[1]);
  result.errors.velocity_linf = maxima[0];
  result.errors.near_wall_velocity_l2 = std::sqrt(sums[2] / sums[6]);
  result.errors.pressure_l2 = std::sqrt(sums[3] / gauge_values[1]);
  result.errors.pressure_linf = maxima[1];
  result.pressure_error_extremum = pressure_error_extremum;
  if (result.errors.pressure_linf !=
      result.pressure_error_extremum.absolute_error_pa)
    throw std::runtime_error(
        "manufactured pressure-extremum maximum is inconsistent");
  result.errors.near_wall_pressure_l2 = std::sqrt(sums[4] / sums[5]);
  result.near_wall_pressure = near_wall_pressure;
  const auto near_wall_statistics =
      summarize_pressure_error_moments(result.near_wall_pressure.total);
  if (!near_wall_statistics.has_value())
    throw std::runtime_error(
        "manufactured near-wall pressure diagnostic closure failed");
  result.errors.penetration_l2 =
      std::sqrt(penetration_sums[0] / penetration_sums[1]);
  result.errors.penetration_linf = penetration_linf;
  result.exact_velocity_linf = maxima[2];
  result.exact_pressure_linf = maxima[3];
  result.active_pressure_mean = sums[7] / gauge_values[1];

  {
    const auto extremum_local =
        topology.find_local_cell(pressure_error_extremum.global_cell_id);
    const bool owns_extremum =
        extremum_local.has_value() &&
        topology.cell_ownership(*extremum_local) ==
            mesh::EntityOwnership::owned;
    int owner_count = owns_extremum ? 1 : 0;
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, &owner_count, 1, MPI_INT, MPI_SUM,
                            mpi.comm()),
              "MPI_Allreduce(pressure-extremum authority owner count)");
    if (owner_count != 1)
      throw std::runtime_error(
          "manufactured pressure-extremum authority owner is not unique");

    std::array<std::uint64_t, 4> integer_values{};
    std::array<double, 4> real_values{};
    int owner_rank = 0;
    int direct_fluid_cell = 0;
    if (owns_extremum) {
      const auto cell_center = geometry.cell_center_m(*extremum_local);
      const immersed::ImmersedLink *nearest = nullptr;
      double nearest_distance_squared =
          std::numeric_limits<double>::infinity();
      for (const auto &link : domain.links()) {
        const auto offset = subtract(cell_center, link.wall_intercept_m);
        const double distance_squared = dot(offset, offset);
        if (!std::isfinite(distance_squared))
          throw std::runtime_error(
              "manufactured pressure-extremum link distance is non-finite");
        if (nearest == nullptr ||
            std::tie(distance_squared, link.id) <
                std::tie(nearest_distance_squared, nearest->id)) {
          nearest = &link;
          nearest_distance_squared = distance_squared;
        }
      }
      if (nearest == nullptr)
        throw std::runtime_error(
            "manufactured pressure-extremum has no immersed authority link");
      const auto fluid_local = topology.find_local_cell(nearest->fluid_cell);
      if (!fluid_local.has_value())
        throw std::runtime_error(
            "manufactured pressure-extremum authority fluid cell is absent");
      const double h = std::cbrt(geometry.cell_volume_m3(*fluid_local));
      if (!(h > 0.0) || !std::isfinite(h))
        throw std::runtime_error(
            "manufactured pressure-extremum authority scale is invalid");
      const auto &reconstruction = ghost_plan.reconstruction(nearest->id);
      const auto pressure_functional =
          immersed::detail::QuadraticReconstructionWeights::
              origin_normal_gradient_constrained_value_weights(
                  reconstruction, nearest->wall_intercept_m);
      const auto gradient_functional =
          immersed::detail::QuadraticReconstructionWeights::
              origin_normal_gradient_constrained_directional_gradient_weights(
                  reconstruction, nearest->wall_intercept_m,
                  nearest->solid_to_fluid_normal);
      const auto donor_l1 = [](const auto &donors) {
        double value = 0.0;
        for (const auto &donor : donors)
          value += std::abs(donor.weight);
        return value;
      };
      const auto &donor_ids =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              reconstruction);
      std::uint64_t donor_fingerprint = UINT64_C(14695981039346656037);
      for (const auto donor : donor_ids) {
        donor_fingerprint ^= donor;
        donor_fingerprint *= UINT64_C(1099511628211);
      }
      donor_fingerprint ^= static_cast<std::uint64_t>(donor_ids.size());
      donor_fingerprint *= UINT64_C(1099511628211);
      if (donor_fingerprint == 0U)
        donor_fingerprint = 1U;
      const auto row_link_count = static_cast<std::uint64_t>(std::count_if(
          domain.links().begin(), domain.links().end(), [&](const auto &link) {
            return link.fluid_cell == pressure_error_extremum.global_cell_id;
          }));
      if (row_link_count == 0U)
        throw std::runtime_error(
            "manufactured pressure-extremum row has no immersed link");
      integer_values = {nearest->id,
                        static_cast<std::uint64_t>(donor_ids.size()),
                        row_link_count,
                        donor_fingerprint};
      real_values = {
          std::sqrt(nearest_distance_squared),
          reconstruction.quality().condition_estimate,
          donor_l1(pressure_functional.donors) +
              std::abs(pressure_functional.boundary_coefficient) / h,
          h * donor_l1(gradient_functional.donors) +
              std::abs(gradient_functional.boundary_coefficient)};
      owner_rank = mpi.rank();
      direct_fluid_cell =
          nearest->fluid_cell == pressure_error_extremum.global_cell_id ? 1
                                                                         : 0;
    }
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, integer_values.data(),
                            static_cast<int>(integer_values.size()),
                            MPI_UINT64_T, MPI_SUM, mpi.comm()),
              "MPI_Allreduce(pressure-extremum authority integers)");
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, real_values.data(),
                            static_cast<int>(real_values.size()), MPI_DOUBLE,
                            MPI_SUM, mpi.comm()),
              "MPI_Allreduce(pressure-extremum authority values)");
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, &owner_rank, 1, MPI_INT, MPI_SUM,
                            mpi.comm()),
              "MPI_Allreduce(pressure-extremum authority rank)");
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, &direct_fluid_cell, 1, MPI_INT,
                            MPI_SUM, mpi.comm()),
              "MPI_Allreduce(pressure-extremum authority association)");
    if (owner_rank < 0 || owner_rank >= mpi.size() ||
        integer_values[1] == 0U || integer_values[2] == 0U ||
        integer_values[3] == 0U ||
        !std::all_of(real_values.begin(), real_values.end(),
                     [](double value) {
                       return std::isfinite(value) && value >= 0.0;
                     }))
      throw std::runtime_error(
          "manufactured pressure-extremum authority is invalid");
    result.pressure_extremum_authority = {
        integer_values[0], owner_rank, direct_fluid_cell == 1, real_values[0],
        real_values[1], integer_values[1], integer_values[2],
        integer_values[3], real_values[2], real_values[3]};
  }

  if (definition.collect_pressure_extrema_diagnostics) {
    if (mpi.size() != 1)
      throw std::runtime_error(
          "pressure-extrema diagnostics require one MPI rank");
    const auto top_extrema =
        select_pressure_error_extrema(local_pressure_extrema, 8U);
    result.pressure_error_extrema_top_k.reserve(top_extrema.size());
    const auto donor_l1 = [](const auto &donors) {
      double value = 0.0;
      for (const auto &donor : donors)
        value += std::abs(donor.weight);
      return value;
    };
    for (const auto &candidate : top_extrema) {
      PressureErrorExtremumDetail detail;
      detail.error = candidate;
      const auto local = topology.find_local_cell(candidate.global_cell_id);
      if (!local.has_value())
        throw std::runtime_error(
            "pressure-extrema diagnostic cell is not locally owned");
      const double local_h = std::cbrt(geometry.cell_volume_m3(*local));
      if (!(local_h > 0.0) || !std::isfinite(local_h))
        throw std::runtime_error(
            "pressure-extrema diagnostic cell scale is invalid");
      detail.wall_distance_over_h = candidate.wall_distance_m / local_h;
      detail.incident_wall_link_count = static_cast<std::uint64_t>(
          std::count_if(domain.links().begin(), domain.links().end(),
                        [&](const auto &link) {
                          return link.fluid_cell == candidate.global_cell_id;
                        }));

      const auto cell_center = geometry.cell_center_m(*local);
      const immersed::ImmersedLink *nearest = nullptr;
      double nearest_distance_squared =
          std::numeric_limits<double>::infinity();
      for (const auto &link : domain.links()) {
        const auto offset = subtract(cell_center, link.wall_intercept_m);
        const double distance_squared = dot(offset, offset);
        if (!std::isfinite(distance_squared))
          throw std::runtime_error(
              "pressure-extrema diagnostic link distance is non-finite");
        if (nearest == nullptr ||
            std::tie(distance_squared, link.id) <
                std::tie(nearest_distance_squared, nearest->id)) {
          nearest = &link;
          nearest_distance_squared = distance_squared;
        }
      }
      if (nearest == nullptr) {
        result.pressure_error_extrema_top_k.push_back(std::move(detail));
        continue;
      }
      const auto fluid_local = topology.find_local_cell(nearest->fluid_cell);
      if (!fluid_local.has_value())
        throw std::runtime_error(
            "pressure-extrema diagnostic authority cell is absent");
      const double authority_h =
          std::cbrt(geometry.cell_volume_m3(*fluid_local));
      const auto &reconstruction = ghost_plan.reconstruction(nearest->id);
      const auto pressure_functional =
          immersed::detail::QuadraticReconstructionWeights::
              origin_normal_gradient_constrained_value_weights(
                  reconstruction, nearest->wall_intercept_m);
      const auto gradient_functional =
          immersed::detail::QuadraticReconstructionWeights::
              origin_normal_gradient_constrained_directional_gradient_weights(
                  reconstruction, nearest->wall_intercept_m,
                  nearest->solid_to_fluid_normal);
      const auto &donor_ids =
          immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
              reconstruction);
      std::uint64_t donor_fingerprint = UINT64_C(14695981039346656037);
      for (const auto donor : donor_ids) {
        donor_fingerprint ^= donor;
        donor_fingerprint *= UINT64_C(1099511628211);
      }
      donor_fingerprint ^= static_cast<std::uint64_t>(donor_ids.size());
      donor_fingerprint *= UINT64_C(1099511628211);
      if (donor_fingerprint == 0U)
        donor_fingerprint = 1U;
      detail.authority_available = true;
      detail.authority = {
          nearest->id,
          mpi.rank(),
          nearest->fluid_cell == candidate.global_cell_id,
          std::sqrt(nearest_distance_squared),
          reconstruction.quality().condition_estimate,
          static_cast<std::uint64_t>(donor_ids.size()),
          detail.incident_wall_link_count,
          donor_fingerprint,
          donor_l1(pressure_functional.donors) +
              std::abs(pressure_functional.boundary_coefficient) /
                  authority_h,
          authority_h * donor_l1(gradient_functional.donors) +
              std::abs(gradient_functional.boundary_coefficient)};
      result.pressure_error_extrema_top_k.push_back(std::move(detail));
    }
  }

  std::uint64_t local_counts[3]{
      static_cast<std::uint64_t>(domain.active_cells().owned_active_count()),
      static_cast<std::uint64_t>(domain.links().size()),
      static_cast<std::uint64_t>(wall_plan.local_points().size())};
  std::uint64_t global_counts[3]{};
  check_mpi(MPI_Allreduce(local_counts, global_counts, 3, MPI_UINT64_T, MPI_SUM,
                          mpi.comm()),
            "MPI_Allreduce(manufactured counts)");
  result.active_cell_count = global_counts[0];
  result.immersed_link_count = global_counts[1];
  result.wall_point_count = global_counts[2];

  if (attempt.force) {
    if (!force_reference)
      throw std::runtime_error(
          "manufactured force reference was not frozen before the solve");
    result.operator_force = attempt.force->operator_force;
    result.budget_reaction = attempt.force->budget_reaction;
    result.surface_traction = attempt.force->surface_traction;
    result.consistency = attempt.force->consistency;
    result.pressure_measure = evaluate_pressure_measure_diagnostics(
        mpi, topology, geometry, wall_plan, ghost_plan, transform,
        flow::test::ImmersedFlowTestAccess::immersed_operator(immersed_flow),
        state.layer(flow::FlowLayer::committed).view<double>(
            fields.mechanical_pressure),
        flow::test::ImmersedFlowTestAccess::final_pressure_wall_gradients(immersed_flow),
        attempt.force->budget_reaction.pressure_N,
        attempt.force->surface_traction.pressure_N, definition.body, dt,
        force_reference->force.pressure_N, result.pressure_scalar);
    const auto &reference = *force_reference;
    const double pressure_scale =
        kReferenceDensityKgPerM3 * kReferenceVelocityMPerS *
        kReferenceVelocityMPerS * reference.surface_area_m2;
    const double viscous_scale = kDynamicViscosityPaS *
                                 kReferenceVelocityMPerS *
                                 reference.surface_area_m2 / kReferenceLengthM;
    const double total_scale = std::max(pressure_scale, viscous_scale);
    result.analytic_force = reference.force;
    result.pressure_force_scale_N = pressure_scale;

    double local_viscous_traction_error_square = 0.0;
    double local_exact_cell_traction_error_square = 0.0;
    double local_numerical_cell_traction_error_square = 0.0;
    double local_traction_split_closure_linf = 0.0;
    double local_functional_amplification_square = 0.0;
    double local_functional_weight = 0.0;
    double local_functional_amplification_max = 0.0;
    double local_functional_condition_max = 0.0;
    runtime::Real3 local_signed_traction_error{};
    runtime::Real3 local_absolute_component_traction_error{};
    double local_absolute_traction_error = 0.0;
    std::uint64_t local_functional_point_count = 0U;
    runtime::Real3 local_plan_full_gradient_force{};
    runtime::Real3 local_link_projected_force{};
    runtime::Real3 local_link_full_gradient_force{};
    const auto committed_velocity_view =
        state.layer(flow::FlowLayer::committed).view<double>(fields.velocity);
    std::map<mesh::GlobalCellId, MmsCellAverage> remote_exact_donor_averages;
    const auto exact_donor_average = [&](mesh::GlobalCellId global_cell)
        -> const MmsCellAverage & {
      const auto local = topology.find_local_cell(global_cell);
      if (local.has_value() &&
          domain.region(*local) != immersed::CellRegion::fluid)
        throw std::runtime_error(
            "manufactured exact wall donor is unavailable");
      if (local.has_value() &&
          topology.cell_ownership(*local) == mesh::EntityOwnership::owned) {
        if (*local >= next_averages.size())
          throw std::runtime_error(
              "manufactured exact owned wall donor is out of range");
        return next_averages[*local];
      }
      const auto found = remote_exact_donor_averages.find(global_cell);
      if (found != remote_exact_donor_averages.end())
        return found->second;
      const auto logical =
          logical_cell_from_global_id(global_cell, geometry.global_extent());
      const auto inserted = remote_exact_donor_averages.emplace(
          global_cell,
          detail::evaluate_cell_average_at_global(geometry, logical,
                                                  definition.body, dt));
      return inserted.first->second;
    };
    struct TractionErrorPoint final {
      immersed::ImmersedLinkId link{};
      std::vector<mesh::GlobalCellId> donor_global_ids;
      runtime::Real3 error{};
      double weight_m2{};
    };
    std::map<immersed::TriangleId, std::vector<TractionErrorPoint>>
        traction_error_by_triangle;
    std::map<immersed::ImmersedLinkId, mesh::GlobalCellId>
        global_link_fluid_cells;
    {
      const auto &local_links = domain.links();
      std::vector<std::uint64_t> local_pairs;
      local_pairs.reserve(local_links.size() * 2U);
      for (const auto &link : local_links) {
        local_pairs.push_back(link.id);
        local_pairs.push_back(link.fluid_cell);
      }
      std::vector<int> counts(static_cast<std::size_t>(mpi.size()), 0);
      std::vector<int> displs(static_cast<std::size_t>(mpi.size()), 0);
      const int local_count = static_cast<int>(local_pairs.size());
      check_mpi(MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1,
                              MPI_INT, mpi.comm()),
                "MPI_Allgather(manufactured wall link counts)");
      for (std::size_t index = 1U; index < counts.size(); ++index)
        displs[index] = displs[index - 1U] + counts[index - 1U];
      std::vector<std::uint64_t> all_pairs(
          static_cast<std::size_t>(displs.back() + counts.back()), 0U);
      check_mpi(MPI_Allgatherv(local_pairs.data(), local_count, MPI_UINT64_T,
                               all_pairs.data(), counts.data(), displs.data(),
                               MPI_UINT64_T, mpi.comm()),
                "MPI_Allgatherv(manufactured wall link pairs)");
      for (std::size_t index = 0U; index + 1U < all_pairs.size();
           index += 2U)
        global_link_fluid_cells.emplace(all_pairs[index], all_pairs[index + 1U]);
    }
    for (std::size_t point_index = 0U;
         point_index < wall_plan.local_points().size(); ++point_index) {
      const auto &point = wall_plan.local_points()[point_index];
      const auto &reconstruction = point.reconstruction;
      const auto exact = evaluate_mms(definition.body, point.position_m, dt);
      std::array<double, 9> numerical_gradient{};
      std::array<double, 9> exact_cell_gradient{};
      std::array<double, 9> plan_full_gradient{};
      std::array<double, 9> link_projected_gradient{};
      std::array<double, 9> link_full_gradient{};
      const auto &link_reconstruction = ghost_plan.reconstruction(
          immersed::detail::boundary_authority_link(point.reconstruction));
      const auto exact_gradient_weights =
          immersed::detail::QuadraticReconstructionWeights::
              origin_constrained_directional_gradient_weights(
                  point.reconstruction, point.position_m,
                  point.solid_to_fluid_normal);
      for (std::size_t component = 0U; component < 3U; ++component) {
        const auto raw = immersed::detail::gradient_with_origin_constraint(
            reconstruction, point.position_m, committed_velocity_view,
            component, 0.0);
        const auto link_raw = immersed::detail::gradient_with_origin_constraint(
            link_reconstruction, point.position_m, committed_velocity_view,
            component, 0.0);
        const double numerical_normal = dot(raw, point.solid_to_fluid_normal);
        const double link_normal = dot(link_raw, point.solid_to_fluid_normal);
        double exact_cell_normal = 0.0;
        for (const auto &donor : exact_gradient_weights) {
          const auto &average = exact_donor_average(donor.global_cell);
          const double value = component == 0U   ? average.velocity_m_per_s.x
                               : component == 1U ? average.velocity_m_per_s.y
                                                 : average.velocity_m_per_s.z;
          exact_cell_normal += donor.weight * value;
        }
        if (!std::isfinite(exact_cell_normal))
          throw std::runtime_error(
              "manufactured exact wall gradient is non-finite");
        for (std::size_t axis = 0U; axis < 3U; ++axis) {
          const double normal_component =
              axis == 0U   ? point.solid_to_fluid_normal.x
              : axis == 1U ? point.solid_to_fluid_normal.y
                           : point.solid_to_fluid_normal.z;
          numerical_gradient[component * 3U + axis] =
              numerical_normal * normal_component;
          exact_cell_gradient[component * 3U + axis] =
              exact_cell_normal * normal_component;
          plan_full_gradient[component * 3U + axis] = axis == 0U   ? raw.x
                                                      : axis == 1U ? raw.y
                                                                   : raw.z;
          link_projected_gradient[component * 3U + axis] =
              link_normal * normal_component;
          link_full_gradient[component * 3U + axis] = axis == 0U   ? link_raw.x
                                                      : axis == 1U ? link_raw.y
                                                                   : link_raw.z;
        }
      }
      const auto numerical_traction = viscous_traction_from_gradient(
          numerical_gradient, point.solid_to_fluid_normal);
      const auto exact_cell_traction = viscous_traction_from_gradient(
          exact_cell_gradient, point.solid_to_fluid_normal);
      local_plan_full_gradient_force =
          add(local_plan_full_gradient_force,
              multiply(point.weight_m2,
                       viscous_traction_from_gradient(
                           plan_full_gradient, point.solid_to_fluid_normal)));
      local_link_projected_force =
          add(local_link_projected_force,
              multiply(point.weight_m2, viscous_traction_from_gradient(
                                            link_projected_gradient,
                                            point.solid_to_fluid_normal)));
      local_link_full_gradient_force =
          add(local_link_full_gradient_force,
              multiply(point.weight_m2,
                       viscous_traction_from_gradient(
                           link_full_gradient, point.solid_to_fluid_normal)));
      const auto true_surface_traction =
          viscous_traction(exact, point.solid_to_fluid_normal);
      const auto traction_error =
          subtract(numerical_traction, true_surface_traction);
      const auto exact_cell_traction_error =
          subtract(exact_cell_traction, true_surface_traction);
      const auto numerical_cell_traction_error =
          subtract(numerical_traction, exact_cell_traction);
      const auto split_traction_error =
          add(numerical_cell_traction_error, exact_cell_traction_error);
      const double split_scale =
          std::max({std::numeric_limits<double>::min(),
                    max_abs(traction_error), max_abs(split_traction_error),
                    max_abs(numerical_traction),
                    max_abs(true_surface_traction)});
      const double split_defect =
          max_abs_difference(traction_error, split_traction_error);
      if (!std::isfinite(split_defect) ||
          split_defect > 2048.0 * std::numeric_limits<double>::epsilon() *
                             split_scale)
        throw std::runtime_error(
            "manufactured wall traction error split is inconsistent");
      local_traction_split_closure_linf =
          std::max(local_traction_split_closure_linf, split_defect);
      local_exact_cell_traction_error_square +=
          point.weight_m2 *
          dot(exact_cell_traction_error, exact_cell_traction_error);
      local_numerical_cell_traction_error_square +=
          point.weight_m2 *
          dot(numerical_cell_traction_error, numerical_cell_traction_error);
      const auto authority_link =
          immersed::detail::boundary_authority_link(point.reconstruction);
      const auto link_fluid = global_link_fluid_cells.find(authority_link);
      if (link_fluid == global_link_fluid_cells.end())
        throw std::runtime_error(
            "manufactured wall functional authority link is absent");
      const auto fluid_logical = logical_cell_from_global_id(
          link_fluid->second, geometry.global_extent());
      const double h =
          std::cbrt(detail::cell_volume_global(geometry, fluid_logical));
      const auto normal_gradient_weights =
          immersed::detail::QuadraticReconstructionWeights::
              origin_constrained_directional_gradient_weights(
                  point.reconstruction, point.position_m,
                  point.solid_to_fluid_normal);
      double normal_gradient_amplification = 0.0;
      for (const auto &donor : normal_gradient_weights)
        normal_gradient_amplification += std::abs(donor.weight);
      normal_gradient_amplification *= h;
      const double condition =
          point.reconstruction.quality().condition_estimate;
      if (!(h > 0.0) || !std::isfinite(h) ||
          !std::isfinite(normal_gradient_amplification) ||
          !(normal_gradient_amplification >= 0.0) ||
          !std::isfinite(condition) || !(condition >= 0.0))
        throw std::runtime_error(
            "manufactured wall functional conditioning is invalid");
      local_functional_amplification_square +=
          point.weight_m2 * normal_gradient_amplification *
          normal_gradient_amplification;
      local_functional_weight += point.weight_m2;
      local_functional_amplification_max = std::max(
          local_functional_amplification_max, normal_gradient_amplification);
      local_functional_condition_max =
          std::max(local_functional_condition_max, condition);
      local_signed_traction_error =
          add(local_signed_traction_error,
              multiply(point.weight_m2, traction_error));
      local_absolute_component_traction_error =
          add(local_absolute_component_traction_error,
              {point.weight_m2 * std::abs(traction_error.x),
               point.weight_m2 * std::abs(traction_error.y),
               point.weight_m2 * std::abs(traction_error.z)});
      local_absolute_traction_error +=
          point.weight_m2 * std::sqrt(dot(traction_error, traction_error));
      ++local_functional_point_count;
      local_viscous_traction_error_square +=
          point.weight_m2 * dot(traction_error, traction_error);
      traction_error_by_triangle[point.triangle].push_back(
          {immersed::detail::boundary_authority_link(point.reconstruction),
           immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
               point.reconstruction),
           traction_error, point.weight_m2});
    }
    constexpr std::size_t kSameLink = 0U;
    constexpr std::size_t kCrossLinkSameDonor = 1U;
    constexpr std::size_t kCrossDonor = 2U;
    std::array<double, 3> local_jump_square{};
    std::array<double, 3> local_jump_weight{};
    std::array<double, 3> local_jump_maximum{};
    std::array<std::uint64_t, 3> local_jump_count{};
    const double traction_scale = viscous_scale / reference.surface_area_m2;
    if (!(traction_scale > 0.0) || !std::isfinite(traction_scale))
      throw std::runtime_error(
          "manufactured traction-jump scale is invalid");
    for (const auto &[triangle, points] : traction_error_by_triangle) {
      static_cast<void>(triangle);
      if (points.size() != 3U)
        throw std::runtime_error(
            "manufactured triangle traction support is incomplete");
      for (std::size_t first = 0U; first < points.size(); ++first) {
        for (std::size_t second = first + 1U; second < points.size(); ++second) {
          const auto difference = subtract(points[first].error,
                                           points[second].error);
          const double normalized_jump =
              std::sqrt(dot(difference, difference)) / traction_scale;
          const double pair_weight =
              0.5 * (points[first].weight_m2 + points[second].weight_m2);
          if (!std::isfinite(normalized_jump) || !(pair_weight > 0.0) ||
              !std::isfinite(pair_weight))
            throw std::runtime_error(
                "manufactured traction-jump sample is invalid");
          const std::size_t category =
              points[first].link == points[second].link
                  ? kSameLink
                  : points[first].donor_global_ids ==
                            points[second].donor_global_ids
                        ? kCrossLinkSameDonor
                        : kCrossDonor;
          local_jump_square[category] +=
              pair_weight * normalized_jump * normalized_jump;
          local_jump_weight[category] += pair_weight;
          local_jump_maximum[category] =
              std::max(local_jump_maximum[category], normalized_jump);
          ++local_jump_count[category];
        }
      }
    }
    std::array<double, 3> traction_error_squares{
        local_viscous_traction_error_square,
        local_exact_cell_traction_error_square,
        local_numerical_cell_traction_error_square};
    sum_in_place(mpi, traction_error_squares.data(),
                 traction_error_squares.size());
    max_in_place(mpi, &local_traction_split_closure_linf, 1U);
    std::array<double, 9> conditioning_sums{
        local_functional_amplification_square,
        local_functional_weight,
        local_signed_traction_error.x,
        local_signed_traction_error.y,
        local_signed_traction_error.z,
        local_absolute_component_traction_error.x,
        local_absolute_component_traction_error.y,
        local_absolute_component_traction_error.z,
        local_absolute_traction_error};
    sum_in_place(mpi, conditioning_sums.data(), conditioning_sums.size());
    std::array<double, 2> conditioning_maxima{
        local_functional_amplification_max,
        local_functional_condition_max};
    max_in_place(mpi, conditioning_maxima.data(), conditioning_maxima.size());
    check_mpi(MPI_Allreduce(MPI_IN_PLACE, &local_functional_point_count, 1,
                            MPI_UINT64_T, MPI_SUM, mpi.comm()),
              "MPI_Allreduce(manufactured wall functional count)");
    if (!(conditioning_sums[1] > 0.0) ||
        !(conditioning_sums[0] >= 0.0) ||
        local_functional_point_count != result.wall_point_count ||
        !std::all_of(conditioning_sums.begin(), conditioning_sums.end(),
                     [](double value) { return std::isfinite(value); }) ||
        !std::all_of(conditioning_maxima.begin(), conditioning_maxima.end(),
                     [](double value) {
                       return std::isfinite(value) && value >= 0.0;
                     }))
      throw std::runtime_error(
          "manufactured wall functional reduction is invalid");
    const runtime::Real3 signed_error{conditioning_sums[2],
                                      conditioning_sums[3],
                                      conditioning_sums[4]};
    const auto independent_signed_error = subtract(
        result.surface_traction.viscous_N,
        result.wall_plan_analytic_force.viscous_N);
    const double signed_scale =
        std::max({1.0, max_abs(signed_error),
                  max_abs(independent_signed_error)});
    if (max_abs_difference(signed_error, independent_signed_error) >
        2048.0 * std::numeric_limits<double>::epsilon() * signed_scale)
      throw std::runtime_error(
          "manufactured wall functional signed error is inconsistent");
    result.wall_functional_conditioning = {
        std::sqrt(conditioning_sums[0] / conditioning_sums[1]),
        conditioning_maxima[0],
        conditioning_maxima[1],
        signed_error,
        {conditioning_sums[5], conditioning_sums[6], conditioning_sums[7]},
        conditioning_sums[8],
        conditioning_sums[8] > 0.0
            ? magnitude(signed_error) / conditioning_sums[8]
            : 0.0,
        local_functional_point_count};
    const auto normalized_traction_l2 = [&](double error_square) {
      if (!std::isfinite(error_square) || !(error_square >= 0.0))
        throw std::runtime_error(
            "manufactured wall traction decomposition is invalid");
      return std::sqrt(reference.surface_area_m2 * error_square) /
             viscous_scale;
    };
    result.wall_traction_error_decomposition = {
        normalized_traction_l2(traction_error_squares[1]),
        normalized_traction_l2(traction_error_squares[2]),
        normalized_traction_l2(traction_error_squares[0]),
        local_traction_split_closure_linf,
        local_functional_point_count};
    sum_in_place(mpi, local_jump_square.data(), local_jump_square.size());
    sum_in_place(mpi, local_jump_weight.data(), local_jump_weight.size());
    max_in_place(mpi, local_jump_maximum.data(), local_jump_maximum.size());
    std::array<std::uint64_t, 3> global_jump_count{};
    check_mpi(MPI_Allreduce(local_jump_count.data(), global_jump_count.data(),
                            static_cast<int>(local_jump_count.size()),
                            MPI_UINT64_T, MPI_SUM, mpi.comm()),
              "MPI_Allreduce(manufactured traction-jump counts)");
    const auto rms = [&](std::size_t category) {
      if (local_jump_weight[category] == 0.0) {
        if (global_jump_count[category] != 0U ||
            local_jump_square[category] != 0.0 ||
            local_jump_maximum[category] != 0.0)
          throw std::runtime_error(
              "manufactured traction-jump empty category is inconsistent");
        return 0.0;
      }
      if (global_jump_count[category] == 0U ||
          !(local_jump_square[category] >= 0.0))
        throw std::runtime_error(
            "manufactured traction-jump category is inconsistent");
      return std::sqrt(local_jump_square[category] /
                       local_jump_weight[category]);
    };
    result.wall_traction_jump = {
        rms(kSameLink),
        local_jump_maximum[kSameLink],
        rms(kCrossLinkSameDonor),
        local_jump_maximum[kCrossLinkSameDonor],
        rms(kCrossDonor),
        local_jump_maximum[kCrossDonor],
        global_jump_count[kSameLink],
        global_jump_count[kCrossLinkSameDonor],
        global_jump_count[kCrossDonor]};
    std::array<double, 9> diagnostic_forces{
        local_plan_full_gradient_force.x, local_plan_full_gradient_force.y,
        local_plan_full_gradient_force.z, local_link_projected_force.x,
        local_link_projected_force.y,     local_link_projected_force.z,
        local_link_full_gradient_force.x, local_link_full_gradient_force.y,
        local_link_full_gradient_force.z};
    sum_in_place(mpi, diagnostic_forces.data(), diagnostic_forces.size());
    result.wall_plan_full_gradient_viscous_force_N = {
        diagnostic_forces[0], diagnostic_forces[1], diagnostic_forces[2]};
    result.wall_link_projected_viscous_force_N = {
        diagnostic_forces[3], diagnostic_forces[4], diagnostic_forces[5]};
    result.wall_link_full_gradient_viscous_force_N = {
        diagnostic_forces[6], diagnostic_forces[7], diagnostic_forces[8]};
    result.errors.viscous_traction_l2 =
        normalized_traction_l2(traction_error_squares[0]);
    result.errors.pressure_force =
        vector_error(result.surface_traction.pressure_N,
                     reference.force.pressure_N, pressure_scale);
    result.errors.viscous_force =
        vector_error(result.surface_traction.viscous_N,
                     reference.force.viscous_N, viscous_scale);
    result.errors.total_force = vector_error(
        result.surface_traction.total_N, reference.force.total_N, total_scale);
    result.errors.pressure_consistency =
        max_abs(result.consistency.pressure_N) / pressure_scale;
    result.errors.viscous_consistency =
        max_abs(result.consistency.viscous_N) / viscous_scale;
    result.errors.total_consistency =
        max_abs(result.consistency.total_N) / total_scale;
  }

  const std::size_t global_cells = static_cast<std::size_t>(definition.cells) *
                                   static_cast<std::size_t>(definition.cells) *
                                   static_cast<std::size_t>(definition.cells);
  result.global_velocity.assign(global_cells * 3U, 0.0);
  result.global_pressure.assign(global_cells, 0.0);
  std::vector<double> active_mask(global_cells, 0.0);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto global_id =
        static_cast<std::size_t>(topology.global_cell_id(cell));
    active_mask[global_id] = 1.0;
    result.global_pressure[global_id] =
        committed.mechanical_pressure[cell] - pressure_error_mean;
    for (std::size_t component = 0U; component < 3U; ++component)
      result.global_velocity[global_id * 3U + component] =
          committed.velocity[cell * 3U + component];
  }
  sum_in_place(mpi, active_mask.data(), active_mask.size());
  sum_in_place(mpi, result.global_velocity.data(),
               result.global_velocity.size());
  sum_in_place(mpi, result.global_pressure.data(),
               result.global_pressure.size());
  result.global_active_cell_ids.reserve(
      static_cast<std::size_t>(result.active_cell_count));
  for (std::size_t global_id = 0U; global_id < global_cells; ++global_id) {
    if (active_mask[global_id] == 1.0)
      result.global_active_cell_ids.push_back(
          static_cast<mesh::GlobalCellId>(global_id));
    else if (active_mask[global_id] != 0.0)
      throw std::runtime_error("manufactured global active mask is not unique");
  }
  if (result.global_active_cell_ids.size() != result.active_cell_count)
    throw std::runtime_error(
        "manufactured global active-cell count is inconsistent");
  return result;
}

} // namespace hundun::test::stage3
