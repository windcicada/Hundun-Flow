// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <cmath>

namespace hundun::v04 {
namespace {

constexpr double kKappa = 0.41;
constexpr double kLogE = 9.8;
constexpr double kViscousLogIntersection = 11.5301072732;
constexpr double kVanDriestA = 26.0;

double dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double magnitude(Real3 value) noexcept {
  return std::sqrt(dot(value, value));
}

double friction_velocity(double tangential_speed, double distance,
                         double kinematic_viscosity) noexcept {
  if (tangential_speed == 0.0) {
    return 0.0;
  }
  const double viscous =
      std::sqrt(tangential_speed * kinematic_viscosity / distance);
  if (distance * viscous / kinematic_viscosity <=
      kViscousLogIntersection) {
    return viscous;
  }
  double lower = kViscousLogIntersection * kinematic_viscosity / distance;
  double upper = std::max(tangential_speed, 2.0 * lower);
  const auto residual = [&](double value) noexcept {
    const double y_plus = distance * value / kinematic_viscosity;
    return value * std::log(kLogE * y_plus) / kKappa -
           tangential_speed;
  };
  while (residual(upper) < 0.0 && std::isfinite(upper) &&
         upper < 1.0e150) {
    upper *= 2.0;
  }
  for (std::uint8_t iteration = 0U; iteration < 64U; ++iteration) {
    const double midpoint = 0.5 * (lower + upper);
    if (residual(midpoint) < 0.0) {
      lower = midpoint;
    } else {
      upper = midpoint;
    }
  }
  return 0.5 * (lower + upper);
}

}  // namespace

Status TurbulencePlan::evaluate_wall_function(
    const WallFunctionSample& sample,
    WallFunctionResult& result) const noexcept {
  if (fingerprint_ == 0U ||
      wall_ != WallTreatmentKind::equilibrium_wall_function ||
      (sample.surface != WallSurfaceKind::external &&
       sample.surface != WallSurfaceKind::immersed) ||
      !finite(sample.solid_to_fluid_normal) ||
      !finite(sample.fluid_velocity) || !finite(sample.wall_velocity) ||
      !std::isfinite(sample.wall_distance) ||
      !(sample.wall_distance > 0.0) || !std::isfinite(sample.density) ||
      !(sample.density > 0.0) ||
      !std::isfinite(sample.molecular_viscosity) ||
      !(sample.molecular_viscosity > 0.0) ||
      !std::isfinite(sample.heat_capacity) ||
      !(sample.heat_capacity > 0.0) ||
      !std::isfinite(sample.molecular_conductivity) ||
      sample.molecular_conductivity < 0.0 ||
      !std::isfinite(sample.fluid_temperature) ||
      !std::isfinite(sample.wall_temperature) ||
      !std::isfinite(sample.molecular_mass_diffusivity) ||
      sample.molecular_mass_diffusivity < 0.0 ||
      !std::isfinite(sample.fluid_scalar) ||
      !std::isfinite(sample.wall_scalar) ||
      !std::isfinite(sample.roughness_height) ||
      sample.roughness_height != 0.0 ||
      std::abs(dot(sample.solid_to_fluid_normal,
                   sample.solid_to_fluid_normal) -
               1.0) >
          1.0e-10) {
    return {StatusCode::invalid_plan, 2501U};
  }
  const Real3 relative{sample.fluid_velocity.x - sample.wall_velocity.x,
                       sample.fluid_velocity.y - sample.wall_velocity.y,
                       sample.fluid_velocity.z - sample.wall_velocity.z};
  const double normal_speed = dot(relative, sample.solid_to_fluid_normal);
  const Real3 tangential{
      relative.x - normal_speed * sample.solid_to_fluid_normal.x,
      relative.y - normal_speed * sample.solid_to_fluid_normal.y,
      relative.z - normal_speed * sample.solid_to_fluid_normal.z};
  const double tangential_speed = magnitude(tangential);
  const double molecular_kinematic =
      sample.molecular_viscosity / sample.density;
  const double u_tau = friction_velocity(tangential_speed,
                                         sample.wall_distance,
                                         molecular_kinematic);
  const double y_plus = sample.wall_distance * u_tau / molecular_kinematic;
  const double damping = 1.0 - std::exp(-y_plus / kVanDriestA);
  const double wall_kinematic =
      kKappa * sample.wall_distance * u_tau * damping * damping;
  WallFunctionResult candidate;
  if (tangential_speed > 0.0) {
    const double traction = -sample.density * u_tau * u_tau /
                            tangential_speed;
    candidate.shear_on_fluid = {traction * tangential.x,
                                traction * tangential.y,
                                traction * tangential.z};
  }
  candidate.friction_velocity = u_tau;
  candidate.y_plus = y_plus;
  candidate.wall_kinematic_viscosity = wall_kinematic;
  const double thermal_conductivity =
      sample.molecular_conductivity +
      sample.density * sample.heat_capacity * wall_kinematic /
          turbulent_prandtl_;
  const double scalar_mass_conductivity =
      sample.density *
      (sample.molecular_mass_diffusivity +
       wall_kinematic / turbulent_schmidt_);
  candidate.heat_flux_into_fluid =
      thermal_conductivity *
      (sample.wall_temperature - sample.fluid_temperature) /
      sample.wall_distance;
  candidate.scalar_flux_into_fluid =
      scalar_mass_conductivity *
      (sample.wall_scalar - sample.fluid_scalar) / sample.wall_distance;
  if (!std::isfinite(u_tau) || u_tau < 0.0 || !std::isfinite(y_plus) ||
      y_plus < 0.0 || !std::isfinite(wall_kinematic) ||
      wall_kinematic < 0.0 || !finite(candidate.shear_on_fluid) ||
      !std::isfinite(candidate.heat_flux_into_fluid) ||
      !std::isfinite(candidate.scalar_flux_into_fluid)) {
    return {StatusCode::numerical_failure, 2502U};
  }
  result = candidate;
  return {};
}

}  // namespace hundun::v04
