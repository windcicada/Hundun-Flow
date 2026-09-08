// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_spray.hpp"

namespace hundun::v04::spray::detail {

// Signed, new-minus-old, extensive changes. Multiplicity and the normalized
// deposition weight have already been applied exactly once. These are actual
// endpoint changes, not instantaneous rates or a force-only approximation.
struct ParcelExchangeBudget {
  double mass_delta_kg{};
  Vector3 momentum_delta_kg_m_per_s{};
  double thermochemical_enthalpy_delta_j{};
  double kinetic_energy_delta_j{};
  // Independently integrated vapor absolute enthalpy minus heat to liquid.
  // Its sum with actual liquid ΔH must satisfy the explicit tolerance below.
  double thermal_exchange_to_gas_j{};
};

struct GasCellExchangeInput {
  double gas_mass_kg{};
  Vector3 gas_momentum_kg_m_per_s{};
  ParcelExchangeBudget parcel{};
  double thermal_absolute_tolerance_j{};
  double thermal_relative_tolerance{};
};

enum class GasCellExchangeStatus : std::uint8_t {
  success,
  invalid_input,
  exhausted_gas_mass,
  inconsistent_thermal_budget,
  non_finite_output
};

struct GasCellExchangeCandidate {
  GasCellExchangeStatus status{GasCellExchangeStatus::invalid_input};
  std::string_view model_id{"closed_isobaric_exchange_v1"};
  bool available{};
  double gas_mass_delta_kg{};
  Vector3 gas_momentum_delta_kg_m_per_s{};
  double gas_mass_candidate_kg{};
  Vector3 gas_momentum_candidate_kg_m_per_s{};
  double gas_kinetic_energy_delta_j{};
  double gas_thermochemical_enthalpy_delta_j{};
  double thermal_budget_residual_j{};
  double total_energy_residual_j{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == GasCellExchangeStatus::success && available;
  }
};

// Closed interphase exchange at fixed thermodynamic pressure: the pair's
// extensive H_tc + K is conserved. H_tc includes formation enthalpy, not K.
// Gas ΔH = -parcel ΔH - parcel ΔK - gas ΔK; using only -parcel ΔH would
// discard the kinetic energy converted by drag and vapor/gas mixing.
//
// No external force, wall/outlet work, radiation, or unsteady pressure-work
// contribution may be included in the parcel budget. A driver must combine
// those separately in its governing-equation convention. All budgets assigned
// to this cell must be accumulated before this nonlinear kinetic correction.
// Inconsistent flux/state thermal budgets are rejected, not corrected silently.
// The default zero tolerances require an exact thermal budget; an adapter must
// explicitly supply any allowed integration-error budget and handle rejection.
// Species mapping, element constraints, pressure work, EOS/positivity checks,
// source admission and state publication remain accepted-head adapter work.
[[nodiscard]] GasCellExchangeCandidate make_gas_cell_exchange_candidate(
    const GasCellExchangeInput& input) noexcept;

}  // namespace hundun::v04::spray::detail
