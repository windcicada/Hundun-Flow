// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_boundary.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected, double relative = 2.0e-13,
           double absolute = 2.0e-13) {
  return std::abs(actual - expected) <=
         absolute + relative * std::max(std::abs(actual), std::abs(expected));
}

NscbcPrimitive primitive(double normal_velocity) {
  return NscbcPrimitive{1.2, 101325.0, 300.0, normal_velocity, 3.0, -2.0,
                        340.0, 1.4};
}

NscbcTarget target() {
  NscbcTarget value;
  value.pressure = 101000.0;
  value.temperature = 290.0;
  value.normal_velocity = -34.0;
  value.tangent_velocity_1 = 1.0;
  value.tangent_velocity_2 = -1.0;
  value.pressure_relaxation = 2.0;
  value.velocity_relaxation = 0.25;
  value.temperature_relaxation = 0.5;
  value.mach_limit = 0.95;
  value.reversal_tolerance = 1.0e-8;
  return value;
}

NscbcNormalGradient gradient() {
  return NscbcNormalGradient{0.002, 5.0, 0.01, 0.02, -0.03};
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec value;
  value.flow_kind = BoundaryKind::no_slip_wall;
  value.thermal_kind = BoundaryKind::adiabatic_wall;
  return value;
}

bool compile_nscbc_plan(BoundaryPlan& plan) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = Real3{0.0, 0.0, 0.0};
  mesh.upper = Real3{4.0, 4.0, 4.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = Int3{4, 4, 4};
  mesh.minimum_spacing = Real3{1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 64U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 20U;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch)) {
    return false;
  }
  ValidatedModel model;
  model.fingerprint = 901U;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face = wall();
  }
  BoundaryFaceSpec& inlet = model.boundaries[0U];
  inlet.flow_kind = BoundaryKind::nscbc_inlet;
  inlet.thermal_kind = BoundaryKind::none;
  inlet.velocity = Real3{34.0, 1.0, -1.0};
  inlet.pressure = 101000.0;
  inlet.temperature = 290.0;
  inlet.relaxation = 0.25;
  inlet.mach_limit = 0.9;
  BoundaryFaceSpec& outlet = model.boundaries[1U];
  outlet.flow_kind = BoundaryKind::nscbc_outlet;
  outlet.thermal_kind = BoundaryKind::none;
  outlet.pressure = 101000.0;
  outlet.temperature = 290.0;
  outlet.backflow_velocity = Real3{-2.0, 0.5, -0.25};
  outlet.backflow_temperature = 285.0;
  outlet.relaxation = 0.25;
  outlet.mach_limit = 0.9;
  outlet.allow_backflow = true;
  FieldRegistry registry;
  SchemePlan schemes;
  TimeSchemePlan time;
  return static_cast<bool>(BoundaryCompiler::compile(
      MPI_COMM_SELF, model, geometry, patch, registry, plan, schemes, time));
}

bool test_outlet_analytic() {
  const NscbcPrimitive state = primitive(68.0);
  const NscbcNormalGradient normal = gradient();
  const NscbcTarget condition = target();
  NscbcWaves waves;
  NscbcRates rates;
  bool passed = expect(static_cast<bool>(evaluate_nscbc(
                           BoundaryKind::nscbc_outlet, state, normal,
                           condition, waves, rates)),
                       "subsonic outlet evaluates");
  const double outgoing = (68.0 + 340.0) * (5.0 + 1.2 * 340.0 * 0.01);
  const double incoming = 2.0 * (101325.0 - 101000.0);
  const double entropy =
      68.0 * (340.0 * 340.0 * 0.002 - 5.0);
  const double sum = 0.5 * (incoming + outgoing);
  passed &= expect(close(waves.acoustic_incoming, incoming) &&
                       close(waves.acoustic_outgoing, outgoing) &&
                       close(waves.entropy, entropy) &&
                       close(waves.tangent_1, 68.0 * 0.02) &&
                       close(waves.tangent_2, 68.0 * -0.03),
                   "outlet replaces only the incoming acoustic wave");
  passed &= expect(close(rates.pressure, -sum) &&
                       close(rates.density,
                             -(entropy + sum) / (340.0 * 340.0)) &&
                       close(rates.normal_velocity,
                             (incoming - outgoing) /
                                 (2.0 * 1.2 * 340.0)) &&
                       close(rates.tangent_velocity_1, -68.0 * 0.02) &&
                       close(rates.tangent_velocity_2, 68.0 * 0.03) &&
                       !rates.used_backflow,
                   "outlet primitive rates follow the analytic inverse");
  const double expected_temperature =
      300.0 / (1.4 * 101325.0) * (entropy - 0.4 * sum);
  passed &= expect(close(rates.temperature, expected_temperature),
                   "outlet pressure/temperature closure uses absolute pressure");
  return passed;
}

bool test_outlet_target_temperature_is_unused_for_normal_outflow() {
  NscbcTarget condition = target();
  condition.temperature = 0.0;
  NscbcWaves waves;
  NscbcRates rates;
  return expect(static_cast<bool>(evaluate_nscbc(
                    BoundaryKind::nscbc_outlet, primitive(10.0), gradient(),
                    condition, waves, rates)) &&
                    !rates.used_backflow,
                "normal NSCBC outflow does not require an unused target temperature");
}

bool test_inlet_analytic() {
  const NscbcPrimitive state = primitive(-34.0);
  const NscbcNormalGradient normal = gradient();
  const NscbcTarget condition = target();
  NscbcWaves waves;
  NscbcRates rates;
  bool passed = expect(static_cast<bool>(evaluate_nscbc(
                           BoundaryKind::nscbc_inlet, state, normal,
                           condition, waves, rates)),
                       "subsonic inlet evaluates");
  const double outgoing = (-34.0 + 340.0) * (5.0 + 1.2 * 340.0 * 0.01);
  const double requested_normal = -0.25 * (-34.0 - -34.0);
  const double incoming = outgoing + 2.0 * 1.2 * 340.0 * requested_normal;
  const double sum = 0.5 * (incoming + outgoing);
  const double requested_temperature = -0.5 * (300.0 - 290.0);
  const double entropy = 1.4 * 101325.0 / 300.0 * requested_temperature +
                         0.4 * sum;
  passed &= expect(close(waves.acoustic_outgoing, outgoing) &&
                       close(waves.acoustic_incoming, incoming),
                   "inlet retains the outgoing acoustic wave");
  passed &= expect(close(waves.entropy, entropy) &&
                       close(waves.tangent_1, 0.25 * (3.0 - 1.0)) &&
                       close(waves.tangent_2, 0.25 * (-2.0 - -1.0)),
                   "inlet supplies entropy and tangential incoming waves");
  passed &= expect(close(rates.normal_velocity, requested_normal) &&
                       close(rates.temperature, requested_temperature) &&
                       close(rates.tangent_velocity_1,
                             -0.25 * (3.0 - 1.0)) &&
                       close(rates.tangent_velocity_2,
                             -0.25 * (-2.0 - -1.0)),
                   "inlet rates reproduce requested target relaxation");
  return passed;
}

bool test_normal_sign_and_backflow() {
  NscbcTarget condition = target();
  NscbcWaves retained_waves{1.0, 2.0, 3.0, 4.0, 5.0};
  NscbcRates retained_rates{6.0, 7.0, 8.0, 9.0, 10.0, 11.0, false};
  const Status reversed_inlet = evaluate_nscbc(
      BoundaryKind::nscbc_inlet, primitive(1.0), gradient(), condition,
      retained_waves, retained_rates);
  bool passed = expect(reversed_inlet.code == StatusCode::invalid_plan,
                       "positive outward velocity is rejected at inlet");
  const Status reversed_outlet = evaluate_nscbc(
      BoundaryKind::nscbc_outlet, primitive(-1.0), gradient(), condition,
      retained_waves, retained_rates);
  passed &= expect(reversed_outlet.code == StatusCode::invalid_plan,
                   "negative outward velocity is rejected at outlet by default");
  condition.allow_backflow = true;
  condition.backflow_temperature = 285.0;
  condition.backflow_normal_velocity = -2.0;
  condition.backflow_tangent_velocity_1 = 0.5;
  condition.backflow_tangent_velocity_2 = -0.25;
  NscbcWaves waves;
  NscbcRates rates;
  passed &= expect(static_cast<bool>(evaluate_nscbc(
                       BoundaryKind::nscbc_outlet, primitive(-1.0), gradient(),
                       condition, waves, rates)) &&
                       rates.used_backflow &&
                       close(rates.normal_velocity, -0.25 * (-1.0 - -2.0)) &&
                       close(rates.temperature, -0.5 * (300.0 - 285.0)),
                   "declared outlet backflow switches to inlet target closure");
  return passed;
}

bool test_zero_mach_and_local_sound_speed() {
  NscbcTarget condition = target();
  condition.normal_velocity = 0.0;
  NscbcPrimitive zero = primitive(0.0);
  NscbcWaves inlet_waves;
  NscbcRates inlet_rates;
  bool passed = expect(static_cast<bool>(evaluate_nscbc(
                           BoundaryKind::nscbc_inlet, zero, gradient(),
                           condition, inlet_waves, inlet_rates)),
                       "zero-Mach inlet has a finite declared-kind limit");
  passed &= expect(std::isfinite(inlet_rates.density) &&
                       std::isfinite(inlet_rates.pressure) &&
                       std::isfinite(inlet_rates.temperature) &&
                       std::isfinite(inlet_rates.normal_velocity),
                   "zero-Mach evaluation never divides by velocity");

  NscbcWaves outlet_waves;
  NscbcRates outlet_rates;
  passed &= expect(static_cast<bool>(evaluate_nscbc(
                       BoundaryKind::nscbc_outlet, zero, gradient(), condition,
                       outlet_waves, outlet_rates)),
                   "zero-Mach outlet has a finite declared-kind limit");
  NscbcPrimitive slower_sound = zero;
  slower_sound.sound_speed = 300.0;
  NscbcWaves changed;
  NscbcRates changed_rates;
  passed &= expect(static_cast<bool>(evaluate_nscbc(
                       BoundaryKind::nscbc_outlet, slower_sound, gradient(),
                       condition, changed, changed_rates)) &&
                       !close(changed.acoustic_outgoing,
                              outlet_waves.acoustic_outgoing, 1.0e-15,
                              1.0e-15),
                   "characteristics use the supplied local sound speed");
  return passed;
}

bool test_rejections_and_atomic_outputs() {
  const NscbcWaves original_waves{1.0, 2.0, 3.0, 4.0, 5.0};
  const NscbcRates original_rates{6.0, 7.0, 8.0, 9.0, 10.0, 11.0, true};
  NscbcWaves waves = original_waves;
  NscbcRates rates = original_rates;
  NscbcTarget condition = target();
  NscbcPrimitive sonic = primitive(340.0 * condition.mach_limit);
  Status result = evaluate_nscbc(BoundaryKind::nscbc_outlet, sonic, gradient(),
                                 condition, waves, rates);
  bool passed = expect(result.code == StatusCode::invalid_plan,
                       "sonic/supersonic-capable state is rejected");
  passed &= expect(waves.acoustic_incoming == original_waves.acoustic_incoming &&
                       rates.density == original_rates.density &&
                       rates.used_backflow == original_rates.used_backflow,
                   "rejection publishes no partial output");
  NscbcPrimitive invalid = primitive(10.0);
  invalid.pressure = std::numeric_limits<double>::quiet_NaN();
  result = evaluate_nscbc(BoundaryKind::nscbc_outlet, invalid, gradient(),
                          condition, waves, rates);
  passed &= expect(result.code == StatusCode::numerical_failure,
                   "non-finite local thermo state is rejected");
  result = evaluate_nscbc(BoundaryKind::pressure_outlet, primitive(10.0),
                          gradient(), condition, waves, rates);
  passed &= expect(result.code == StatusCode::invalid_plan,
                   "non-NSCBC boundary kind is rejected");
  return passed;
}

bool test_compiled_plan_integration() {
  BoundaryPlan plan;
  bool passed = expect(compile_nscbc_plan(plan),
                       "NSCBC boundary descriptors compile");
  if (!passed) {
    return false;
  }
  const BoundaryFacePlan* inlet = nullptr;
  const BoundaryFacePlan* outlet = nullptr;
  passed &= expect(plan.face(CartesianFace::x_min, inlet) &&
                       plan.face(CartesianFace::x_max, outlet) &&
                       inlet != nullptr && outlet != nullptr,
                   "compiled NSCBC face descriptors are addressable");
  if (inlet == nullptr || outlet == nullptr) {
    return false;
  }
  passed &= expect(inlet->flow_kernel == BoundaryFlowKernel::characteristic &&
                       outlet->flow_kernel ==
                           BoundaryFlowKernel::characteristic,
                   "compiled NSCBC faces retain characteristic dispatch");
  for (std::size_t index = 0U; index < plan.spans().size; ++index) {
    const BoundaryIndexSpan& span = plan.spans().data[index];
    if ((span.face == CartesianFace::x_min ||
         span.face == CartesianFace::x_max) &&
        (span.stage == BoundaryStage::momentum ||
         span.stage == BoundaryStage::pressure)) {
      passed &= expect(span.relation != BoundaryRelation::dirichlet,
                       "NSCBC stencil ghosts do not double-impose U/pi");
    }
  }
  NscbcTarget compiled_target;
  passed &= expect(static_cast<bool>(
                       plan.nscbc_target(CartesianFace::x_max,
                                         compiled_target)) &&
                       compiled_target.allow_backflow &&
                       close(compiled_target.pressure, 101000.0) &&
                       close(compiled_target.backflow_temperature, 285.0) &&
                       close(compiled_target.backflow_normal_velocity, -2.0) &&
                       close(compiled_target.backflow_tangent_velocity_1,
                             0.5) &&
                       close(compiled_target.backflow_tangent_velocity_2,
                             -0.25),
                   "compiled NSCBC target retains backflow authority");
  NscbcWaves waves;
  NscbcRates rates;
  passed &= expect(static_cast<bool>(plan.evaluate_nscbc(
                       CartesianFace::x_max, primitive(10.0), gradient(),
                       waves, rates)) &&
                       std::isfinite(rates.pressure),
                   "compiled descriptor directly dispatches characteristic kernel");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_outlet_analytic();
  passed &= test_outlet_target_temperature_is_unused_for_normal_outflow();
  passed &= test_inlet_analytic();
  passed &= test_normal_sign_and_backflow();
  passed &= test_zero_mach_and_local_sound_speed();
  passed &= test_rejections_and_atomic_outputs();
  passed &= test_compiled_plan_integration();
  if (passed) {
    std::cout << "v0.4 NSCBC analytic tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
