// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"

namespace hundun::v04::detail {

// The contribution registry groups views by conserved field. Select a group
// without allocating; assembly validates its stage, units, identity and count.
inline bool select_species_sources(Span<const EquationContributionView> all,
    FieldId field, Span<const EquationContributionView>& selected) noexcept {
  selected = {};
  if (all.size && !all.data) return false;
  std::size_t begin = all.size, count = 0;
  bool ended = false;
  for (std::size_t i = 0; i < all.size; ++i) {
    if (all.data[i].conserved_quantity == field) {
      if (ended) return false;
      if (begin == all.size) begin = i;
      ++count;
    } else if (count) ended = true;
  }
  if (count) selected = {all.data + begin, count};
  return true;
}

// Controls work on rejected coupling iterates, never time-step acceptance.
// A stalled composition update requests the complete inner solve next time.
class SpeciesCouplingForcing {
 public:
  void begin(unsigned sweep) noexcept {
    if (sweep == 1U || sweep != sweep_ + 1U) {
      sampled_ = false;
      full_ = false;
      previous_ = 0.0;
    }
    sweep_ = sweep;
  }
  void observe(double residual) noexcept {
    full_ = !std::isfinite(residual) || residual < 0.0 ||
            (sampled_ && residual >= 0.9 * previous_);
    previous_ = residual;
    sampled_ = true;
  }
  bool can_update_composition(double continuity, double energy,
                              double final_continuity,
                              double final_energy) const noexcept {
    if (full_ || !std::isfinite(continuity) || continuity < 0.0 ||
        !std::isfinite(energy) || energy < 0.0 ||
        !std::isfinite(final_continuity) || final_continuity <= 0.0 ||
        !std::isfinite(final_energy) || final_energy <= 0.0 ||
        (sampled_ && previous_ <= final_energy))
      return false;
    const double c = sampled_ ? 0.01 * previous_ : std::sqrt(final_continuity);
    const double e = sampled_ ? 0.01 * previous_ : std::sqrt(final_energy);
    return c > final_continuity && e > final_energy && continuity <= c &&
           energy <= e;
  }

 private:
  unsigned sweep_{};
  double previous_{};
  bool sampled_{}, full_{};
};

// A single FP64 subtraction preserves a direction just beyond half an ULP.
// An intermediate extended-precision subtraction can round that direction
// to the midpoint before conversion back to the stored FP64 composition.
inline double species_search_update(double value, double correction) noexcept {
  return value-correction;
}

// A private search may exhaust the represented value's half-ULP response.
// This predicate concerns an uncommitted guess, separately from final audit.
inline bool species_search_quantized(double value, double proposal,
    double diagonal, double residual) noexcept {
  if (!std::isfinite(value) || value < 0 || value > 1 || proposal != value ||
      !std::isfinite(diagonal) || diagonal <= 0 ||
      !std::isfinite(residual) || residual == 0) return false;
  const double adjacent=std::nextafter(value, residual > 0
      ? -std::numeric_limits<double>::infinity()
      : std::numeric_limits<double>::infinity());
  if (adjacent < 0 || adjacent > 1) return false;
  const long double spacing=std::abs(static_cast<long double>(adjacent)-value);
  return std::abs(static_cast<long double>(residual)) <= .5L*diagonal*spacing;
}

// Positive local convective response used by the private nonlinear search.
// The physical residual and public species assembly are independent of it.
inline double species_coupling_search_diagonal(const CartesianKernelPlan& kernels,
    ConvectionScheme scheme,ConstFaceFluxView flux,Int3 cell,double diagonal) noexcept {
  const std::array<ConstFaceFieldView,3U> faces{flux.x,flux.y,flux.z};
  for(unsigned a=0U;a<3U;++a) {
    Int3 upper=cell; (a==0U ? upper.x : a==1U ? upper.y : upper.z)++;
    double high_response=1.0,low_response=1.0;
    if(scheme!=ConvectionScheme::central2) {
      const auto axis=static_cast<CartesianAxis>(a);
      const int normal=a==0U ? cell.x : a==1U ? cell.y : cell.z;
      const double centre=centre_coordinate(kernels,axis,normal);
      const double lower_distance=centre-face_coordinate(kernels,axis,normal);
      const double upper_distance=face_coordinate(kernels,axis,normal+1)-centre;
      // MC can reconstruct 2*q on a donor face, not just q. Its exact
      // one-sided donor response is bounded by 1+limiter*d_out/d_opposite.
      // This includes branch ties at q=0. The centred branch is no larger;
      // the neighbour reconstruction's response to this cell is <=limiter.
      high_response=1.0+kernels.limiter()*upper_distance/lower_distance;
      low_response=1.0+kernels.limiter()*lower_distance/upper_distance;
      if(scheme==ConvectionScheme::limited_central2) {
        high_response=0.5*(high_response+kernels.limiter());
        low_response=0.5*(low_response+kernels.limiter());
      }
    }
    diagonal+=(high_response*std::max(faces[a].unchecked(upper),0.0)+
        low_response*std::max(-faces[a].unchecked(cell),0.0))/cell_volume(kernels,cell);
  }
  return diagonal;
}

// Private initial-guess assembly on the lagged, uncertified predictor flux.
// It cannot certify final conservation. The public species assembler still
// accepts only final_conservative with a matching final-flux certificate.
Status assemble_species_guess(const SpeciesEquationPlan& plan,
    std::size_t species, const EquationStateView& state,
    const EquationMaterialView& material, const EquationAssemblyContext& context,
    EquationSystemView system, EquationAssemblyCertificate& certificate) noexcept;

// Private nonlinear solve rows, divided by cell volume BEFORE subtracting
// their terms, to avoid integral underflow for representable trace species.
// Accepts only the same full-domain provisional-guess or certified-final
// scopes. Does not expose a certificate with a different unit interpretation.
// Diagonal/rhs/residual are rate densities; face diffusion coefficients retain
// their ordinary integrated transmissibility. Public assembly stays integral.
Status assemble_species_coupling_rows(const SpeciesEquationPlan& plan,
    std::size_t species, const EquationStateView& state,
    const EquationMaterialView& material, const EquationAssemblyContext& context,
    EquationSystemView system,
    Span<const EquationContributionView> sources = {}) noexcept;

// Re-evaluate the current species equation using system.diagonal prepared by
// assemble_species_coupling_rows. The private caller holds density, material,
// flux, time, geometry and boundary fixed until this solve ends. Only the
// species iterate changes. Retains face coefficients and all physical checks.
Status assemble_species_coupling_residual(const SpeciesEquationPlan& plan,
    std::size_t species, const EquationStateView& state,
    const EquationMaterialView& material, const EquationAssemblyContext& context,
    EquationSystemView system,
    Span<const EquationContributionView> sources = {}) noexcept;

} // namespace hundun::v04::detail
