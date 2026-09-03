// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_boundary.hpp"

#include <cmath>
#include <cstdint>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kNscbcKind = 401U;
constexpr std::uint32_t kNscbcState = 402U;
constexpr std::uint32_t kNscbcTarget = 403U;
constexpr std::uint32_t kNscbcRegime = 404U;
constexpr std::uint32_t kNscbcResult = 405U;
constexpr std::uint32_t kNscbcPlan = 406U;

constexpr std::size_t face_index(CartesianFace face) noexcept {
  return static_cast<std::size_t>(face);
}

bool valid_parameter(const BoundaryPlan& plan,
                     std::uint32_t parameter) noexcept {
  const std::size_t count = plan.parameter_count();
  return parameter != kInvalidBoundaryParameter && parameter < count &&
         plan.velocity_x().size == count &&
         plan.velocity_y().size == count &&
         plan.velocity_z().size == count &&
         plan.backflow_velocity_x().size == count &&
         plan.backflow_velocity_y().size == count &&
         plan.backflow_velocity_z().size == count &&
         plan.pressure_targets().size == count &&
         plan.temperature_targets().size == count &&
         plan.backflow_temperature_targets().size == count &&
         plan.relaxation_rates().size == count &&
         plan.mach_limits().size == count &&
         plan.allow_backflow().size == count;
}

void cartesian_components(CartesianFace face, Real3 vector,
                          double& normal, double& tangent_1,
                          double& tangent_2) noexcept {
  const std::size_t index = face_index(face);
  const double sign = (index & 1U) == 0U ? -1.0 : 1.0;
  if (index < 2U) {
    normal = sign * vector.x;
    tangent_1 = vector.y;
    tangent_2 = vector.z;
  } else if (index < 4U) {
    normal = sign * vector.y;
    tangent_1 = vector.x;
    tangent_2 = vector.z;
  } else {
    normal = sign * vector.z;
    tangent_1 = vector.x;
    tangent_2 = vector.y;
  }
}

bool finite(double value) noexcept { return std::isfinite(value); }

bool finite_primitive(const NscbcPrimitive& value) noexcept {
  return finite(value.density) && finite(value.pressure) &&
         finite(value.temperature) && finite(value.normal_velocity) &&
         finite(value.tangent_velocity_1) &&
         finite(value.tangent_velocity_2) && finite(value.sound_speed) &&
         finite(value.gamma);
}

bool finite_gradient(const NscbcNormalGradient& value) noexcept {
  return finite(value.density) && finite(value.pressure) &&
         finite(value.normal_velocity) &&
         finite(value.tangent_velocity_1) &&
         finite(value.tangent_velocity_2);
}

bool finite_target(const NscbcTarget& value) noexcept {
  return finite(value.pressure) && finite(value.temperature) &&
         finite(value.normal_velocity) &&
         finite(value.tangent_velocity_1) &&
         finite(value.tangent_velocity_2) &&
         finite(value.pressure_relaxation) &&
         finite(value.velocity_relaxation) &&
         finite(value.temperature_relaxation) && finite(value.mach_limit) &&
         finite(value.reversal_tolerance) &&
         finite(value.backflow_temperature) &&
         finite(value.backflow_normal_velocity) &&
         finite(value.backflow_tangent_velocity_1) &&
         finite(value.backflow_tangent_velocity_2);
}

bool finite_waves(const NscbcWaves& value) noexcept {
  return finite(value.acoustic_incoming) && finite(value.entropy) &&
         finite(value.tangent_1) && finite(value.tangent_2) &&
         finite(value.acoustic_outgoing);
}

bool finite_rates(const NscbcRates& value) noexcept {
  return finite(value.density) && finite(value.pressure) &&
         finite(value.temperature) && finite(value.normal_velocity) &&
         finite(value.tangent_velocity_1) &&
         finite(value.tangent_velocity_2);
}

void primitive_rates(const NscbcPrimitive& primitive,
                     const NscbcWaves& waves, NscbcRates& rates) noexcept {
  const double acoustic_sum =
      0.5 * (waves.acoustic_incoming + waves.acoustic_outgoing);
  const double sound_speed_squared =
      primitive.sound_speed * primitive.sound_speed;
  rates.pressure = -acoustic_sum;
  rates.density = -(waves.entropy + acoustic_sum) / sound_speed_squared;
  rates.normal_velocity =
      (waves.acoustic_incoming - waves.acoustic_outgoing) /
      (2.0 * primitive.density * primitive.sound_speed);
  rates.tangent_velocity_1 = -waves.tangent_1;
  rates.tangent_velocity_2 = -waves.tangent_2;
  rates.temperature =
      primitive.temperature /
      (primitive.gamma * primitive.pressure) *
      (waves.entropy - (primitive.gamma - 1.0) * acoustic_sum);
}

void set_inlet_waves(const NscbcPrimitive& primitive,
                     const NscbcTarget& target, NscbcWaves& waves) noexcept {
  const double requested_normal_rate =
      -target.velocity_relaxation *
      (primitive.normal_velocity - target.normal_velocity);
  waves.acoustic_incoming =
      waves.acoustic_outgoing +
      2.0 * primitive.density * primitive.sound_speed *
          requested_normal_rate;
  const double acoustic_sum =
      0.5 * (waves.acoustic_incoming + waves.acoustic_outgoing);
  const double requested_temperature_rate =
      -target.temperature_relaxation *
      (primitive.temperature - target.temperature);
  waves.entropy =
      primitive.gamma * primitive.pressure / primitive.temperature *
          requested_temperature_rate +
      (primitive.gamma - 1.0) * acoustic_sum;
  waves.tangent_1 =
      target.velocity_relaxation *
      (primitive.tangent_velocity_1 - target.tangent_velocity_1);
  waves.tangent_2 =
      target.velocity_relaxation *
      (primitive.tangent_velocity_2 - target.tangent_velocity_2);
}

}  // namespace

Status BoundaryPlan::pressure_perturbation_target(CartesianFace selected,
                                                  double p_ref,
                                                  double& out) const noexcept {
  if (static_cast<std::uint8_t>(selected) >
          static_cast<std::uint8_t>(CartesianFace::z_max) ||
      !std::isfinite(p_ref)) {
    return {StatusCode::invalid_plan, kNscbcPlan};
  }
  const BoundaryFacePlan* selected_face = nullptr;
  const Status face_status = face(selected, selected_face);
  if (!face_status || selected_face == nullptr) {
    return {StatusCode::invalid_plan, kNscbcPlan};
  }
  const bool has_direct_pressure_target =
      selected_face->flow_kind == BoundaryKind::pressure_outlet ||
      selected_face->flow_kind == BoundaryKind::static_state_inlet;
  if (!valid_parameter(*this, selected_face->flow_parameter) ||
      !has_direct_pressure_target) {
    return {StatusCode::invalid_plan, kNscbcPlan};
  }
  const double target = pressure_targets().data[selected_face->flow_parameter] -
                        p_ref;
  if (!std::isfinite(target)) {
    return {StatusCode::numerical_failure, kNscbcResult};
  }
  out = target;
  return {};
}

Status BoundaryPlan::nscbc_target(CartesianFace selected,
                                  NscbcTarget& out) const noexcept {
  if (static_cast<std::uint8_t>(selected) >
      static_cast<std::uint8_t>(CartesianFace::z_max)) {
    return {StatusCode::invalid_plan, kNscbcPlan};
  }
  const BoundaryFacePlan* selected_face = nullptr;
  const Status face_status = face(selected, selected_face);
  if (!face_status || selected_face == nullptr ||
      selected_face->flow_kernel != BoundaryFlowKernel::characteristic ||
      !valid_parameter(*this, selected_face->flow_parameter)) {
    return {StatusCode::invalid_plan, kNscbcPlan};
  }
  const std::uint32_t parameter = selected_face->flow_parameter;
  Real3 velocity{velocity_x().data[parameter], velocity_y().data[parameter],
                 velocity_z().data[parameter]};
  Real3 backflow{backflow_velocity_x().data[parameter],
                 backflow_velocity_y().data[parameter],
                 backflow_velocity_z().data[parameter]};
  NscbcTarget candidate;
  candidate.pressure = pressure_targets().data[parameter];
  candidate.temperature = temperature_targets().data[parameter];
  cartesian_components(selected, velocity, candidate.normal_velocity,
                       candidate.tangent_velocity_1,
                       candidate.tangent_velocity_2);
  candidate.pressure_relaxation = relaxation_rates().data[parameter];
  candidate.velocity_relaxation = relaxation_rates().data[parameter];
  candidate.temperature_relaxation = relaxation_rates().data[parameter];
  candidate.mach_limit = mach_limits().data[parameter];
  candidate.allow_backflow = allow_backflow().data[parameter] != 0U;
  candidate.backflow_temperature =
      backflow_temperature_targets().data[parameter];
  cartesian_components(selected, backflow,
                       candidate.backflow_normal_velocity,
                       candidate.backflow_tangent_velocity_1,
                       candidate.backflow_tangent_velocity_2);
  out = candidate;
  return {};
}

Status BoundaryPlan::evaluate_nscbc(
    CartesianFace selected, const NscbcPrimitive& primitive,
    const NscbcNormalGradient& gradient, NscbcWaves& waves,
    NscbcRates& rates) const noexcept {
  NscbcTarget target;
  const Status status = nscbc_target(selected, target);
  if (!status) {
    return status;
  }
  const BoundaryFacePlan* selected_face = nullptr;
  const Status face_status = face(selected, selected_face);
  if (!face_status || selected_face == nullptr) {
    return {StatusCode::invalid_plan, kNscbcPlan};
  }
  return hundun::v04::evaluate_nscbc(selected_face->flow_kind, primitive,
                                     gradient, target, waves, rates);
}

Status evaluate_nscbc(BoundaryKind kind, const NscbcPrimitive& primitive,
                      const NscbcNormalGradient& gradient,
                      const NscbcTarget& target, NscbcWaves& waves,
                      NscbcRates& rates) noexcept {
  if (kind != BoundaryKind::nscbc_inlet &&
      kind != BoundaryKind::nscbc_outlet) {
    return {StatusCode::invalid_plan, kNscbcKind};
  }
  if (!finite_primitive(primitive) || !finite_gradient(gradient) ||
      primitive.density <= 0.0 || primitive.pressure <= 0.0 ||
      primitive.temperature <= 0.0 || primitive.sound_speed <= 0.0 ||
      primitive.gamma <= 1.0) {
    return {StatusCode::numerical_failure, kNscbcState};
  }
  if (!finite_target(target) || target.pressure <= 0.0 ||
      (kind == BoundaryKind::nscbc_inlet && target.temperature <= 0.0) ||
      target.pressure_relaxation < 0.0 ||
      target.velocity_relaxation < 0.0 ||
      target.temperature_relaxation < 0.0 || target.mach_limit <= 0.0 ||
      target.mach_limit >= 1.0 || target.reversal_tolerance < 0.0 ||
      (target.allow_backflow && target.backflow_temperature <= 0.0)) {
    return {StatusCode::invalid_plan, kNscbcTarget};
  }

  const double mach =
      std::abs(primitive.normal_velocity) / primitive.sound_speed;
  if (!finite(mach) || mach >= target.mach_limit) {
    return {StatusCode::invalid_plan, kNscbcRegime};
  }
  const bool inlet_reversed =
      kind == BoundaryKind::nscbc_inlet &&
      primitive.normal_velocity > target.reversal_tolerance;
  const bool outlet_backflow =
      kind == BoundaryKind::nscbc_outlet &&
      primitive.normal_velocity < -target.reversal_tolerance;
  if (inlet_reversed || (outlet_backflow && !target.allow_backflow)) {
    return {StatusCode::invalid_plan, kNscbcRegime};
  }

  const double density_sound = primitive.density * primitive.sound_speed;
  NscbcWaves candidate_waves;
  candidate_waves.acoustic_incoming =
      (primitive.normal_velocity - primitive.sound_speed) *
      (gradient.pressure - density_sound * gradient.normal_velocity);
  candidate_waves.entropy =
      primitive.normal_velocity *
      (primitive.sound_speed * primitive.sound_speed * gradient.density -
       gradient.pressure);
  candidate_waves.tangent_1 =
      primitive.normal_velocity * gradient.tangent_velocity_1;
  candidate_waves.tangent_2 =
      primitive.normal_velocity * gradient.tangent_velocity_2;
  candidate_waves.acoustic_outgoing =
      (primitive.normal_velocity + primitive.sound_speed) *
      (gradient.pressure + density_sound * gradient.normal_velocity);

  NscbcTarget effective = target;
  const bool inlet_closure =
      kind == BoundaryKind::nscbc_inlet || outlet_backflow;
  if (outlet_backflow) {
    effective.temperature = target.backflow_temperature;
    effective.normal_velocity = target.backflow_normal_velocity;
    effective.tangent_velocity_1 = target.backflow_tangent_velocity_1;
    effective.tangent_velocity_2 = target.backflow_tangent_velocity_2;
  }
  if (inlet_closure) {
    set_inlet_waves(primitive, effective, candidate_waves);
  } else {
    candidate_waves.acoustic_incoming =
        target.pressure_relaxation * (primitive.pressure - target.pressure);
  }

  NscbcRates candidate_rates;
  primitive_rates(primitive, candidate_waves, candidate_rates);
  candidate_rates.used_backflow = outlet_backflow;
  if (!finite_waves(candidate_waves) || !finite_rates(candidate_rates)) {
    return {StatusCode::numerical_failure, kNscbcResult};
  }
  waves = candidate_waves;
  rates = candidate_rates;
  return {};
}

}  // namespace hundun::v04
