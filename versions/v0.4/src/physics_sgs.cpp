// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>

namespace hundun::v04 {
namespace {
struct ProductScale { double mantissa; int exponent; };

ProductScale product_scale(std::initializer_list<double> factors) noexcept {
  ProductScale result{1.0, 0};
  for (double value : factors) {
    if (value == 0.0) return {0.0, 0};
    int exponent{};
    result.mantissa *= std::frexp(value, &exponent);
    result.exponent += exponent;
  }
  return result;
}

double two_thirds_power(ProductScale value) noexcept {
  // Integer exponent splitting avoids forming the possibly unrepresentable
  // nu_t * Delta * |S|^2 product before taking its two-thirds power.
  int quotient = value.exponent / 3;
  int remainder = value.exponent % 3;
  if (remainder < 0) { --quotient; remainder += 3; }
  const double root = std::cbrt(std::ldexp(value.mantissa, remainder));
  return std::ldexp(root * root, 2 * quotient);
}
}  // namespace

Status sgs_state_from_viscosity(const VelocityGradient& gradient,
                               double filter_width_m, double density_kg_m3,
                               double molecular_viscosity_pa_s,
                               double sgs_kinematic_viscosity_m2_s,
                               SgsState& out) noexcept {
  if (!std::isfinite(filter_width_m) || !(filter_width_m > 0.0) ||
      !std::isfinite(density_kg_m3) || !(density_kg_m3 > 0.0) ||
      !std::isfinite(molecular_viscosity_pa_s) || !(molecular_viscosity_pa_s > 0.0) ||
      !std::isfinite(sgs_kinematic_viscosity_m2_s) || sgs_kinematic_viscosity_m2_s < 0.0)
    return {StatusCode::invalid_plan, 2401U};
  double scale = 0.0;
  for (double value : gradient.value) {
    if (!std::isfinite(value)) return {StatusCode::numerical_failure, 2404U};
    scale = std::max(scale, std::abs(value));
  }
  const double mu = molecular_viscosity_pa_s + density_kg_m3 * sgs_kinematic_viscosity_m2_s;
  if (!std::isfinite(mu)) return {StatusCode::numerical_failure, 2404U};
  SgsState candidate{sgs_kinematic_viscosity_m2_s, 0.0, 0.0, 0.0};
  if (scale == 0.0) { out = candidate; return {}; }
  double contraction = 0.0;
  for (unsigned i = 0; i < 3; ++i) {
    const double diagonal = gradient.value[3*i+i] / scale;
    contraction += 2.0 * diagonal * diagonal;
    for (unsigned j = i+1; j < 3; ++j) {
      const double strain = .5*(gradient.value[3*i+j]/scale) +
                            .5*(gradient.value[3*j+i]/scale);
      contraction += 4.0 * strain * strain;
    }
  }
  candidate.kinetic_energy_m2_s2 = two_thirds_power(product_scale(
      {sgs_kinematic_viscosity_m2_s, filter_width_m, contraction, scale, scale}));
  const auto dissipation = product_scale({mu, contraction, scale, scale});
  candidate.dissipation_w_m3 = std::ldexp(dissipation.mantissa, dissipation.exponent);
  int density_exponent{};
  const double density_mantissa = std::frexp(density_kg_m3, &density_exponent);
  candidate.specific_dissipation_m2_s3 = std::ldexp(
      dissipation.mantissa / density_mantissa, dissipation.exponent - density_exponent);
  if (!std::isfinite(candidate.kinetic_energy_m2_s2) ||
      !std::isfinite(candidate.dissipation_w_m3) ||
      !std::isfinite(candidate.specific_dissipation_m2_s3))
    return {StatusCode::numerical_failure, 2404U};
  out = candidate;
  return {};
}
}  // namespace hundun::v04
