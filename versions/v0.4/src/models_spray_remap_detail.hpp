// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_spray_properties_detail.hpp"

namespace hundun::v04::spray::detail {
struct LiquidRemapReport {
  portable::Status status{portable::Status::invalid_input};
  bool available{};
  double temperature_k{}, droplet_diameter_m{}, density_kg_per_m3{};
  double specific_enthalpy_j_per_kg{}, enthalpy_residual_j_per_kg{};
  double geometry_relative_mass_error{};
};
// Cold-path thermodynamic remap. Mass and absolute specific enthalpy are
// explicit source invariants. The caller preserves weight, velocity, identity
// and history; mesh ownership is a separate operation. Failure is canonical.
[[nodiscard]] LiquidRemapReport remap_liquid_mass_enthalpy(
    const LiquidAsset&, double droplet_mass_kg,
    double source_specific_enthalpy_j_per_kg) noexcept;
} // namespace hundun::v04::spray::detail
