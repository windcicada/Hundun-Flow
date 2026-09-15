// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <cmath>

namespace hundun::v04::detail {
// Rachner, DLR-Mitteilung 98-01, p.48, equations 3.10.3/3.10.4.
// Saturated-liquid approximation; callers admit a subcritical material range.
// The scientific coefficient identity was accepted for 624CF on 2026-09-15.
inline double kerosene_liquid_density(double temperature) noexcept {
  return 1.037096e3 - 7.233865e-1 * temperature - 9.255437e3 / (733.0 - temperature);
}
// COAST's gaseous kerosene relation, shared by bulk flow and the film query.
// Callers validate the temperature domain and the positive returned viscosity.
inline double kerosene_vapor_viscosity(double temperature) noexcept {
  return temperature > 400.0
      ? -2.872095e-5 - 1.444046e-8 * temperature + 2.104857e-6 * std::sqrt(temperature)
      : 1.051939e-6 + 1.325433e-8 * temperature + 1.901751e-27 * std::pow(temperature, 8.0);
}
}  // namespace hundun::v04::detail
