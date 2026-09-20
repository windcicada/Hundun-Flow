// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <cmath>

namespace hundun::v04 {

Status smagorinsky_kinematic_viscosity(const VelocityGradient& gradient,
                                      double filter_width, double coefficient,
                                      double& out) noexcept {
  if (!std::isfinite(filter_width) || !(filter_width > 0.0) ||
      !std::isfinite(coefficient) || coefficient < 0.0)
    return {StatusCode::invalid_plan, 2401U};
  double maximum_gradient = 0.0;
  for (double value : gradient.value) {
    if (!std::isfinite(value))
      return {StatusCode::numerical_failure, 2404U};
    maximum_gradient = std::max(maximum_gradient, std::abs(value));
  }
  const double gradient_scale = maximum_gradient > 0.0 &&
      (maximum_gradient < 1e-30 || maximum_gradient > 1e30) ? maximum_gradient : 1.0;
  VelocityGradient scaled = gradient;
  if (gradient_scale != 1.0)
    for (double& value : scaled.value) value /= gradient_scale;
  // REFERENCE gamma_smagorinsky uses the full symmetric strain, including
  // dilatation, and Delta = cbrt(cell volume). Gradients are velocity-major.
  double contraction = 0.0;
  for (unsigned i = 0; i < 3; ++i) {
    const double diagonal = scaled.value[3*i+i];
    contraction += diagonal * diagonal;
    for (unsigned j = i+1; j < 3; ++j) {
      const double strain = .5*scaled.value[3*i+j] + .5*scaled.value[3*j+i];
      contraction += 2.0 * strain * strain;
    }
  }
  const double length = coefficient * filter_width;
  const double candidate = (length * length * std::sqrt(2.0 * contraction)) * gradient_scale;
  if (!std::isfinite(contraction) || !std::isfinite(candidate) || candidate < 0.0)
    return {StatusCode::numerical_failure, 2404U};
  out = candidate;
  return {};
}

}  // namespace hundun::v04
