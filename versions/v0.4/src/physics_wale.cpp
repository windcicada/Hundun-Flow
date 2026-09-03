// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <cmath>

namespace hundun::v04 {

Status wale_kinematic_viscosity(const VelocityGradient& gradient,
                                double filter_width, double coefficient,
                                double& out) noexcept {
  if (!std::isfinite(filter_width) || !(filter_width > 0.0) ||
      !std::isfinite(coefficient) || coefficient < 0.0) {
    return {StatusCode::invalid_plan, 2301U};
  }
  double g[3][3]{};
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      g[i][j] = gradient.value[i * 3U + j];
      if (!std::isfinite(g[i][j])) {
        return {StatusCode::numerical_failure, 2302U};
      }
    }
  }
  double strain_squared = 0.0;
  double squared_gradient[3][3]{};
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      const double strain = 0.5 * (g[i][j] + g[j][i]);
      strain_squared += strain * strain;
      for (std::size_t k = 0U; k < 3U; ++k) {
        squared_gradient[i][j] += g[i][k] * g[k][j];
      }
    }
  }
  const double trace = (squared_gradient[0][0] + squared_gradient[1][1] +
                        squared_gradient[2][2]) /
                       3.0;
  double traceless_squared = 0.0;
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      double value =
          0.5 * (squared_gradient[i][j] + squared_gradient[j][i]);
      if (i == j) {
        value -= trace;
      }
      traceless_squared += value * value;
    }
  }
  if (!std::isfinite(strain_squared) || strain_squared < 0.0 ||
      !std::isfinite(traceless_squared) || traceless_squared < 0.0) {
    return {StatusCode::numerical_failure, 2302U};
  }
  if (coefficient == 0.0 || traceless_squared == 0.0) {
    out = 0.0;
    return {};
  }
  const double numerator = std::pow(traceless_squared, 1.5);
  const double denominator = std::pow(strain_squared, 2.5) +
                             std::pow(traceless_squared, 1.25);
  const double length = coefficient * filter_width;
  const double candidate = length * length * numerator / denominator;
  if (!std::isfinite(numerator) || !std::isfinite(denominator) ||
      !(denominator > 0.0) || !std::isfinite(candidate) || candidate < 0.0) {
    return {StatusCode::numerical_failure, 2302U};
  }
  out = candidate;
  return {};
}

}  // namespace hundun::v04
