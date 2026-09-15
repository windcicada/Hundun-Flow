// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <array>

namespace hundun::v04::detail {

// Keep the low-temperature formation reference and continue h across Tmid.
// Both native PH inversion and the chemistry provider use this operation.
inline double nasa_high_enthalpy_constant(double temperature,
    const std::array<double, 7>& low,
    const std::array<double, 7>& high) noexcept {
  const double t2 = temperature * temperature;
  const double t3 = t2 * temperature;
  const double t4 = t3 * temperature;
  const auto polynomial = [&](const std::array<double, 7>& a) noexcept {
    return a[0] * temperature + a[1] * t2 * 0.5 + a[2] * t3 / 3.0 +
           a[3] * t4 * 0.25 + a[4] * t4 * temperature * 0.2;
  };
  return (polynomial(low) + low[5]) - polynomial(high);
}

} // namespace hundun::v04::detail
