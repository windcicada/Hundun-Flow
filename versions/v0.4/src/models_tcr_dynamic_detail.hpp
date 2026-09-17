// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <array>
#include <cstdint>

namespace hundun::v04::tcr::detail {
// The native dynamic scalar model is distinct from signed-root diagnostics.
enum class SpeciesControlState : std::uint8_t { inactive, direct, projected, similar };
struct SpeciesControl {
  bool available{};
  SpeciesControlState state{SpeciesControlState::inactive};
  double effective{1}, selected{1}, rate_ratio{};
  bool neighbor_completion() const noexcept {
    return available && state == SpeciesControlState::similar;
  }
};
SpeciesControl cdphyso_species_control(double eta, double pdf_rate,
                                       double psr_rate, double weak_rate) noexcept;

// Test-filter moments use one common density/volume normalization. The
// caller owns coordinates, gradients, molecular-weight conversion and halos.
struct DynamicFilterMoments {
  double density{}, delta_squared{}, gradient_squared{};
  double density_delta_squared_gradient_squared{};
  double density_scalar_squared{}, scalar{};
};
struct DynamicFilterDonor {
  double density{}, volume{}, scalar{};
  std::array<double, 3> gradient{};
};
// All coordinates are scalar-specific. For the mixture-fraction channel they
// are dimensionless; species molar/mass conversion is an explicit caller duty.
bool dynamic_filter_moments(const DynamicFilterDonor *donors, unsigned count,
    DynamicFilterMoments &moments) noexcept;
struct DynamicFilterProducts {
  bool available{};
  double m_squared{}, l_times_m{};
};
DynamicFilterProducts dynamic_filter_products(const DynamicFilterMoments &) noexcept;
double dynamic_cd_from_products(double m_squared, double l_times_m) noexcept;
enum class MixingGroup : std::uint8_t { product, reactant, radical };
double dynamic_group_cd(const std::array<double, 3> &groups,
                         MixingGroup group) noexcept;
} // namespace hundun::v04::tcr::detail
