// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace hundun::v04::detail {
// Close a complete species tuple after a linear solve. The dependent species
// absorbs only FP64 roundoff; larger changes reject the whole candidate.
// Enthalpy remains a separate coordinate, and conservation is audited later.
inline bool close_statistical_composition(Span<double> fractions,
                                         std::size_t dependent) noexcept {
  if (!fractions.data || fractions.size < 2 || fractions.size > UINT8_MAX ||
      dependent >= fractions.size) return false;
  constexpr long double budget = 64 * std::numeric_limits<double>::epsilon();
  std::array<double, UINT8_MAX> candidate{};
  long double sum = 0, raw_sum = 0, correction = 0;
  std::size_t largest = dependent == 0 ? 1 : 0;
  for (std::size_t s = 0; s < fractions.size; ++s) {
    if (s == dependent) continue;
    const double q = fractions.data[s];
    if (!std::isfinite(q)) return false;
    raw_sum += q;
    candidate[s] = std::clamp(q, 0., 1.);
    correction += std::abs(static_cast<long double>(candidate[s]) - q);
    if (2 * correction > budget) return false;
    sum += candidate[s];
    if (candidate[s] > candidate[largest]) largest = s;
  }
  if (sum > 1) {
    const long double excess = sum - 1;
    if (2 * (correction + excess) > budget) return false;
    const double before = candidate[largest];
    candidate[largest] = static_cast<double>(static_cast<long double>(before) - excess);
    // Round toward the simplex interior when nearest rounding leaves a tiny
    // excess. A long-double sum also matches accepted-history validation.
    sum = 0;
    for (std::size_t s = 0; s < fractions.size; ++s)
      if (s != dependent) sum += candidate[s];
    if (sum > 1) {
      candidate[largest] = std::nextafter(candidate[largest], 0.);
      sum = 0;
      for (std::size_t s = 0; s < fractions.size; ++s)
        if (s != dependent) sum += candidate[s];
    }
    correction += static_cast<long double>(before) - candidate[largest];
    if (sum > 1 || candidate[largest] < 0 || 2 * correction > budget) return false;
  }
  candidate[dependent] = static_cast<double>(1 - sum);
  if (correction + std::abs(static_cast<long double>(candidate[dependent]) -
                            (1 - raw_sum)) > budget) return false;
  for (std::size_t s = 0; s < fractions.size; ++s) fractions.data[s] = candidate[s];
  return true;
}

// A component solve can temporarily leave the composition simplex while the
// other components still hold their previous iterate. Accepted histories stay
// physical; the enclosing tuple transaction audits the complete trial before
// thermodynamic evaluation or chemistry.
class StatisticalSpecies {
 public:
  static Status assemble(const SpeciesEquationPlan&, std::size_t,
      const EquationStateView&, const EquationMaterialView&,
      Span<const EquationContributionView>, const EquationAssemblyContext&,
      EquationSystemView, EquationAssemblyCertificate&) noexcept;
};
// Statistical total-enthalpy transport shares the scalar matrix. Pressure
// work, ignition and other physical sources enter via compiled contributions;
// the mean total-energy equation retains its own coupled closure.
class StatisticalEnthalpy {
 public:
  static Status assemble(const EnthalpyEquationPlan&, const EquationStateView&,
      ConstFieldView diffusivity, Span<const EquationContributionView>,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&) noexcept;
};
} // namespace hundun::v04::detail
