// SPDX-License-Identifier: Apache-2.0

#include "../support/ibm_force_fixture.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr VelocityGradient kGradient{{
    0.2, 0.3, -0.1, -0.4, -0.05, 0.2, 0.1, -0.15, -0.15}};

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

double pressure(Real3 point, double reference = 0.0) {
  return reference + 2.0 + 0.7 * point.x - 0.4 * point.y + 0.2 * point.z +
         0.1 * (point.x * point.x + point.y * point.y + point.z * point.z);
}

Real3 add(Real3 left, Real3 right) {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Real3 scale(double factor, Real3 value) {
  return {factor * value.x, factor * value.y, factor * value.z};
}

Real3 cross(Real3 left, Real3 right) {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

SurfaceForce independent_oracle(Span<const SurfaceTractionPoint> points,
                                Real3 origin) {
  SurfaceForce force;
  for (std::size_t point_index = 0U; point_index < points.size;
       ++point_index) {
    const SurfaceTractionPoint& point = points.data[point_index];
    const Real3 pressure_force =
        scale(-point.absolute_pressure, point.solid_to_fluid_normal);
    const auto& g = point.velocity_gradient.value;
    const double divergence = g[0U] + g[4U] + g[8U];
    double stress[3][3]{};
    for (std::size_t i = 0U; i < 3U; ++i) {
      for (std::size_t j = 0U; j < 3U; ++j) {
        stress[i][j] = point.effective_viscosity *
                       (g[i * 3U + j] + g[j * 3U + i] -
                        (i == j ? (2.0 / 3.0) * divergence : 0.0));
      }
    }
    const Real3 n = point.solid_to_fluid_normal;
    const Real3 viscous{
        stress[0][0] * n.x + stress[0][1] * n.y + stress[0][2] * n.z,
        stress[1][0] * n.x + stress[1][1] * n.y + stress[1][2] * n.z,
        stress[2][0] * n.x + stress[2][1] * n.y + stress[2][2] * n.z};
    const Real3 total = add(pressure_force, viscous);
    const Real3 arm{point.position.x - origin.x,
                    point.position.y - origin.y,
                    point.position.z - origin.z};
    force.pressure = add(force.pressure,
                         scale(point.weight, pressure_force));
    force.viscous = add(force.viscous, scale(point.weight, viscous));
    force.moment = add(force.moment,
                       scale(point.weight, cross(arm, total)));
  }
  force.total = add(force.pressure, force.viscous);
  return force;
}

bool close(Real3 left, Real3 right, double tolerance = 2.0e-12) {
  return std::abs(left.x - right.x) <= tolerance *
             std::max({1.0, std::abs(left.x), std::abs(right.x)}) &&
         std::abs(left.y - right.y) <= tolerance *
             std::max({1.0, std::abs(left.y), std::abs(right.y)}) &&
         std::abs(left.z - right.z) <= tolerance *
             std::max({1.0, std::abs(left.z), std::abs(right.z)});
}

bool compare(const SurfaceForce& product, const SurfaceForce& oracle,
             std::string_view description, double tolerance = 2.0e-12) {
  const bool matched = close(product.pressure, oracle.pressure, tolerance) &&
                       close(product.viscous, oracle.viscous, tolerance) &&
                       close(product.total, oracle.total, tolerance) &&
                       close(product.moment, oracle.moment, tolerance);
  if (!matched) {
    std::cerr << "product p=" << product.pressure.x << ','
              << product.pressure.y << ',' << product.pressure.z
              << " v=" << product.viscous.x << ',' << product.viscous.y
              << ',' << product.viscous.z << " m=" << product.moment.x
              << ',' << product.moment.y << ',' << product.moment.z
              << " oracle p=" << oracle.pressure.x << ','
              << oracle.pressure.y << ',' << oracle.pressure.z
              << " v=" << oracle.viscous.x << ',' << oracle.viscous.y
              << ',' << oracle.viscous.z << " m=" << oracle.moment.x << ','
              << oracle.moment.y << ',' << oracle.moment.z << '\n';
  }
  return expect(matched, description);
}

SurfaceTractionPoint analytic_point(Real3 position, Real3 normal,
                                    double weight, double reference = 0.0) {
  return {position, normal, weight, pressure(position, reference), kGradient,
          1.9e-5};
}

bool test_plane_sphere_cylinder_point_oracles() {
  bool passed = true;
  const std::array plane{
      analytic_point({-0.5, -0.5, 0.0}, {0.0, 0.0, 1.0}, 0.25),
      analytic_point({0.5, -0.5, 0.0}, {0.0, 0.0, 1.0}, 0.25),
      analytic_point({-0.5, 0.5, 0.0}, {0.0, 0.0, 1.0}, 0.25),
      analytic_point({0.5, 0.5, 0.0}, {0.0, 0.0, 1.0}, 0.25)};
  SurfaceForce product;
  passed &= expect(static_cast<bool>(integrate_surface_traction(
                       {plane.data(), plane.size()}, {}, product)),
                   "product plane-point traction integrates");
  passed &= compare(product,
                    independent_oracle({plane.data(), plane.size()}, {}),
                    "plane pressure/viscous/total/moment match independent oracle");

  const double sphere_weight = 4.0 * kPi / 6.0;
  const std::array sphere{
      analytic_point({1, 0, 0}, {1, 0, 0}, sphere_weight),
      analytic_point({-1, 0, 0}, {-1, 0, 0}, sphere_weight),
      analytic_point({0, 1, 0}, {0, 1, 0}, sphere_weight),
      analytic_point({0, -1, 0}, {0, -1, 0}, sphere_weight),
      analytic_point({0, 0, 1}, {0, 0, 1}, sphere_weight),
      analytic_point({0, 0, -1}, {0, 0, -1}, sphere_weight)};
  passed &= expect(static_cast<bool>(integrate_surface_traction(
                       {sphere.data(), sphere.size()}, {0.1, -0.2, 0.3},
                       product)),
                   "product sphere-point traction integrates");
  passed &= compare(product,
                    independent_oracle({sphere.data(), sphere.size()},
                                       {0.1, -0.2, 0.3}),
                    "sphere pressure/viscous/total/moment match independent oracle");

  std::vector<SurfaceTractionPoint> cylinder;
  constexpr std::size_t count = 64U;
  cylinder.reserve(count);
  for (std::size_t index = 0U; index < count; ++index) {
    const double angle = 2.0 * kPi * static_cast<double>(index) /
                         static_cast<double>(count);
    const Real3 normal{std::cos(angle), std::sin(angle), 0.0};
    cylinder.push_back(analytic_point(normal, normal, 4.0 * kPi / count));
  }
  passed &= expect(static_cast<bool>(integrate_surface_traction(
                       {cylinder.data(), cylinder.size()}, {}, product)),
                   "product cylinder-point traction integrates");
  passed &= compare(product,
                    independent_oracle({cylinder.data(), cylinder.size()}, {}),
                    "cylinder pressure/viscous/total/moment match independent oracle");
  return passed;
}

bool test_product_quadratic_reconstruction_against_oracle() {
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(),
                       "production cube quadrature fixture compiles");
  if (!passed) return false;
  fixture.fill_analytic();
  const FinalSurfaceState state =
      fixture.state(fixture.committed_flux.flux());
  SurfaceForce product;
  passed &= expect(static_cast<bool>(evaluate_surface_force(
                       MPI_COMM_SELF, fixture.quadrature, state, product)),
                   "product final-state force evaluates");
  std::vector<SurfaceTractionPoint> analytic;
  const Span<const SurfaceQuadraturePoint> points =
      fixture.quadrature.local_points();
  analytic.reserve(points.size);
  for (std::size_t index = 0U; index < points.size; ++index) {
    analytic.push_back(analytic_point(
        points.data[index].position,
        points.data[index].solid_to_fluid_normal, points.data[index].weight,
        state.pressure_reference));
  }
  const SurfaceForce oracle =
      independent_oracle({analytic.data(), analytic.size()}, {});
  passed &= compare(product, oracle,
                    "quadratic cell-average reconstruction matches independent final-state oracle",
                    1.0e-10);
  return passed;
}

double cubic_pressure(Real3 point) {
  return point.x * point.x * point.x +
         0.5 * point.y * point.y * point.y -
         0.25 * point.z * point.z * point.z;
}

bool cubic_force_error(std::int32_t cells, double& error) {
  IbmForceFixture fixture;
  if (!fixture.initialize(MPI_COMM_SELF, cells)) return false;
  fixture.fill_analytic(0.0);
  const std::uint8_t ghosts = fixture.pressure.view.ghosts.x;
  const double width = fixture.geometry.x().uniform_width();
  const double correction = width * width / 4.0;
  for (std::int32_t z = -ghosts; z < fixture.patch.cells.z + ghosts; ++z) {
    for (std::int32_t y = -ghosts; y < fixture.patch.cells.y + ghosts; ++y) {
      for (std::int32_t x = -ghosts; x < fixture.patch.cells.x + ghosts; ++x) {
        const Int3 global{fixture.patch.begin.x + x,
                          fixture.patch.begin.y + y,
                          fixture.patch.begin.z + z};
        const double px = fixture.extrapolated_centre(fixture.geometry.x(),
                                                      global.x);
        const double py = fixture.extrapolated_centre(fixture.geometry.y(),
                                                      global.y);
        const double pz = fixture.extrapolated_centre(fixture.geometry.z(),
                                                      global.z);
        fixture.pressure.view.unchecked({x, y, z}, 0U) =
            px * px * px + px * correction +
            0.5 * (py * py * py + py * correction) -
            0.25 * (pz * pz * pz + pz * correction);
      }
    }
  }
  const FinalSurfaceState state = fixture.state(fixture.committed_flux.flux(),
                                                0.0);
  SurfaceForce product;
  if (!evaluate_surface_force(MPI_COMM_SELF, fixture.quadrature, state,
                              product)) {
    return false;
  }
  std::vector<SurfaceTractionPoint> analytic;
  const auto points = fixture.quadrature.local_points();
  analytic.reserve(points.size);
  for (std::size_t index = 0U; index < points.size; ++index) {
    SurfaceTractionPoint point = analytic_point(
        points.data[index].position,
        points.data[index].solid_to_fluid_normal, points.data[index].weight,
        0.0);
    point.absolute_pressure = cubic_pressure(point.position);
    analytic.push_back(point);
  }
  const SurfaceForce oracle =
      independent_oracle({analytic.data(), analytic.size()}, {});
  const double dx = product.pressure.x - oracle.pressure.x;
  const double dy = product.pressure.y - oracle.pressure.y;
  const double dz = product.pressure.z - oracle.pressure.z;
  error = std::sqrt(dx * dx + dy * dy + dz * dz);
  return std::isfinite(error) && error > 0.0;
}

bool test_force_reconstruction_order() {
  const std::array<std::int32_t, 3U> cells{12, 18, 24};
  std::array<double, cells.size()> errors{};
  bool passed = true;
  for (std::size_t level = 0U; level < cells.size(); ++level) {
    passed &= expect(cubic_force_error(cells[level], errors[level]),
                     "cubic pressure force level evaluates");
  }
  if (!passed) return false;
  const double order_0 =
      std::log(errors[0U] / errors[1U]) /
      std::log(static_cast<double>(cells[1U]) / cells[0U]);
  const double order_1 =
      std::log(errors[1U] / errors[2U]) /
      std::log(static_cast<double>(cells[2U]) / cells[1U]);
  passed &= expect(order_0 >= 1.8 && order_1 >= 1.8,
                   "final pressure-force reconstruction converges at order >= 1.8");
  if (!passed) {
    std::cerr << "force errors=" << errors[0U] << ',' << errors[1U] << ','
              << errors[2U] << " orders=" << order_0 << ',' << order_1
              << '\n';
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  bool passed = test_plane_sphere_cylinder_point_oracles();
  passed &= test_product_quadratic_reconstruction_against_oracle();
  passed &= test_force_reconstruction_order();
  MPI_Finalize();
  return passed ? 0 : 1;
}
