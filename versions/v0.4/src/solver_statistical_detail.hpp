// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <cstring>

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

// A basis lives for one frozen field sweep. Only geometry/material face
// conductances are shared; source diagonals, boundary rows, RHS and residuals
// are independently assembled for every component. The enclosing scheduler
// owns these output faces and resets the basis before any unrelated write.
class StatisticalFaceBasis {
 public:
  using Key=std::array<std::uint64_t,25>;
  static Key identify(const EquationAssemblyContext& c,ConstFieldView gamma,
                      const EquationSystemView& out) noexcept {
    std::uint64_t dt{};std::memcpy(&dt,&c.dt,sizeof(dt));
    return {c.geometry,c.boundary,c.thermo,c.transport,c.time,c.face_flux,
      c.face_flux_authority,c.face_flux_storage,c.face_flux_revision_domain,dt,
      gamma.revision,gamma.storage_identity,gamma.revision_domain,
      static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(gamma.base)),
      static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(c.immersed_interface)),
      static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(c.mixture_transport)),
      static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(out.x_coefficient.base)),
      static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(out.y_coefficient.base)),
      static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(out.z_coefficient.base)),
      static_cast<std::uint64_t>(c.box.begin.x),static_cast<std::uint64_t>(c.box.begin.y),
      static_cast<std::uint64_t>(c.box.begin.z),static_cast<std::uint64_t>(c.scope),
      c.mixture_transport ? c.mixture_transport->linearization : 0,
      c.mixture_transport ? c.mixture_transport->face_flux : 0};
  }
  bool take(const Key& key) noexcept {
    const bool reuse=ready_ && key==key_;
    ready_=false;
    if(reuse)++reuses;
    return reuse;
  }
  void publish(const Key& key) noexcept {key_=key;ready_=true;}
  unsigned long long reuses{};
 private:
  Key key_{};
  bool ready_{};
};

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
