// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/turbulence_fixture.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <string_view>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

WallFunctionSample sample() {
  WallFunctionSample value;
  value.surface = WallSurfaceKind::external;
  value.solid_to_fluid_normal = {1.0, 0.0, 0.0};
  value.fluid_velocity = {3.0, 12.0, -2.0};
  value.wall_velocity = {1.0, 2.0, -2.0};
  value.wall_distance = 0.005;
  value.density = 1.2;
  value.molecular_viscosity = 1.8e-5;
  value.heat_capacity = 1005.0;
  value.molecular_conductivity = 0.026;
  value.fluid_temperature = 300.0;
  value.wall_temperature = 330.0;
  value.molecular_mass_diffusivity = 2.0e-5;
  value.fluid_scalar = 0.1;
  value.wall_scalar = 0.2;
  return value;
}

bool test_external_and_immersed_wall_fluxes() {
  TurbulenceFixture fixture;
  bool passed = expect(fixture.initialize(TurbulencePlanSpec{}),
                       "wall-function turbulence fixture compiles");
  WallFunctionResult external;
  WallFunctionSample input = sample();
  passed &= expect(fixture.plan.evaluate_wall_function(input, external) &&
                       external.friction_velocity > 0.0 &&
                       external.y_plus > 0.0 &&
                       external.wall_kinematic_viscosity >= 0.0,
                   "external smooth wall function is finite and nonnegative");
  passed &= expect(std::abs(external.shear_on_fluid.x) < 1.0e-14 &&
                       external.shear_on_fluid.y < 0.0 &&
                       std::abs(external.shear_on_fluid.z) < 1.0e-14,
                   "moving-wall traction uses only relative tangential velocity");
  passed &= expect(external.heat_flux_into_fluid > 0.0 &&
                       external.scalar_flux_into_fluid > 0.0,
                   "hot/high-scalar wall fluxes point consistently into fluid");

  input.surface = WallSurfaceKind::immersed;
  WallFunctionResult immersed;
  passed &= expect(fixture.plan.evaluate_wall_function(input, immersed) &&
                       std::abs(immersed.friction_velocity -
                                external.friction_velocity) < 1.0e-14 &&
                       std::abs(immersed.heat_flux_into_fluid -
                                external.heat_flux_into_fluid) < 1.0e-12,
                   "external and IBM wall faces share one wall-law kernel");

  const WallFunctionResult marker = immersed;
  input.roughness_height = 1.0e-6;
  passed &= expect(fixture.plan.evaluate_wall_function(input, immersed).code ==
                           StatusCode::invalid_plan &&
                       immersed.friction_velocity == marker.friction_velocity,
                   "unsupported rough wall rejects atomically");
  return passed;
}

bool test_smooth_transition_and_resolved_rejection() {
  TurbulenceFixture fixture;
  bool passed = expect(fixture.initialize(TurbulencePlanSpec{}),
                       "transition fixture compiles");
  WallFunctionSample input = sample();
  input.wall_distance = 0.01;
  input.wall_velocity = {};
  input.fluid_velocity = {};
  const double nu = input.molecular_viscosity / input.density;
  constexpr double intersection = 11.5301072732;
  const double switch_speed = intersection * intersection * nu /
                              input.wall_distance;
  input.fluid_velocity.y = 0.99 * switch_speed;
  WallFunctionResult below;
  input.wall_temperature = input.fluid_temperature;
  input.wall_scalar = input.fluid_scalar;
  passed &= expect(fixture.plan.evaluate_wall_function(input, below) &&
                       below.heat_flux_into_fluid == 0.0 &&
                       below.scalar_flux_into_fluid == 0.0,
                   "viscous-side wall state preserves zero heat/scalar jump");
  input.fluid_velocity.y = 1.01 * switch_speed;
  WallFunctionResult above;
  passed &= expect(fixture.plan.evaluate_wall_function(input, above) &&
                       above.friction_velocity > below.friction_velocity &&
                       (above.friction_velocity - below.friction_velocity) /
                               below.friction_velocity <
                           0.05,
                   "viscous/log-law transition is monotone and smooth");

  TurbulencePlanSpec wale;
  wale.kind = TurbulenceKind::wale;
  TurbulenceFixture resolved;
  WallFunctionResult rejected;
  passed &= expect(resolved.initialize(wale) &&
                       resolved.plan.evaluate_wall_function(input, rejected)
                               .code == StatusCode::invalid_plan,
                   "resolved WALE plan cannot enter wall-function branch");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_external_and_immersed_wall_fluxes();
  passed &= test_smooth_transition_and_resolved_rejection();
  MPI_Finalize();
  return passed ? 0 : 1;
}
