// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace hundun::v04::chemistry::detail {
// A stiff integrator can cross an exhausted species' zero by roundoff.
// Bound the total correction, retain the existing conservation audits after
// this operation, and preserve the caller's entire candidate on rejection.
inline bool bound_chemistry_roundoff(std::vector<double> &fractions,
                                     double absolute_tolerance) noexcept {
  if (!(absolute_tolerance > 0) || !std::isfinite(absolute_tolerance))
    return false;
  // The solver controls local error; an interpolated interval endpoint may
  // exceed that tolerance by a few units. The total repair is also capped by
  // FP64 roundoff at the unit composition scale, even for loose solver controls.
  const double budget = 8 * std::min(absolute_tolerance,
      8 * std::numeric_limits<double>::epsilon());
  double correction = 0;
  for (double value : fractions) {
    if (!std::isfinite(value) || value < -budget || value > 1 + budget)
      return false;
    correction += std::abs(value - std::clamp(value, 0., 1.));
    if (correction > budget)
      return false;
  }
  for (double &value : fractions)
    value = std::clamp(value, 0., 1.);
  return true;
}
} // namespace hundun::v04::chemistry::detail
