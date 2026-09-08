// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_spray.hpp"

#include <array>
#include <cstdint>
#include <string_view>

namespace hundun::v04::spray::detail {

inline constexpr std::uint32_t kMaximumBreakupChildParcels = 16U;

enum class BreakupChildStatus : std::uint8_t {
  success,
  invalid_input,
  breakup_not_requested,
  inconsistent_parent_geometry,
  child_id_collision,
  non_finite_output,
  conservation_failure
};

enum class BreakupEnergyDisposition : std::uint8_t {
  balanced,
  supplied_deformation_energy_deficit,
  supplied_deformation_energy_surplus
};

struct SuppliedDiameterSplitParameters {
  std::uint32_t child_parcel_count{};
  std::uint32_t maximum_child_parcels{};
  double supplied_uniform_child_diameter_m{};
  double liquid_density_kg_per_m3{};
  double surface_tension_n_per_m{};
  double liquid_absolute_thermochemical_enthalpy_j_per_kg{};
  double supplied_deformation_energy_per_parent_droplet_j{};
};

struct BreakupChildCandidate {
  bool available{};
  std::uint32_t child_count{};
  std::array<SprayParcelState, kMaximumBreakupChildParcels> children{};
};

struct BreakupConservationReport {
  double parent_total_mass_kg{};
  double children_total_mass_kg{};
  double mass_residual_kg{};

  Vector3 parent_total_momentum_kg_m_per_s{};
  Vector3 children_total_momentum_kg_m_per_s{};
  Vector3 momentum_residual_kg_m_per_s{};

  double parent_bulk_kinetic_energy_j{};
  double children_bulk_kinetic_energy_j{};
  double bulk_kinetic_energy_residual_j{};

  double parent_thermochemical_enthalpy_j{};
  double children_thermochemical_enthalpy_j{};
  double thermochemical_enthalpy_residual_j{};

  double parent_temperature_k{};
  double minimum_child_temperature_k{};
  double maximum_child_temperature_k{};

  // Multiplicity is physical droplet count represented by computational
  // parcels.  It increases as d_parent^3/d_child^3; equality with the parent
  // multiplicity is neither expected nor claimed.
  double parent_multiplicity{};
  double expected_children_total_multiplicity{};
  double children_total_multiplicity{};
  double multiplicity_residual{};

  // These are deliberately separate from bulk kinetic and thermochemical
  // enthalpy.  The supplied-diameter split does not silently convert a deficit
  // or surplus into child velocity, heat, or a gas source.
  double parent_spherical_surface_energy_j{};
  double children_spherical_surface_energy_j{};
  double surface_energy_increase_j{};
  double supplied_deformation_energy_j{};
  double unassigned_deformation_energy_j{};
  BreakupEnergyDisposition energy_disposition{
      BreakupEnergyDisposition::balanced};
};

struct BreakupChildInput {
  SprayParcelState parent{};
  TabBreakupReport tab_trigger{};
  std::uint64_t accepted_step{};
  std::uint64_t breakup_ordinal{};
  SuppliedDiameterSplitParameters parameters{};
};

struct BreakupChildReport {
  BreakupChildStatus status{BreakupChildStatus::invalid_input};
  std::string_view model_id{
      "conservative_supplied_uniform_diameter_split_v1"};
  std::string_view trigger_model_id{};
  std::string_view rng_model_id{};
  ParcelId parent_id{};
  std::uint64_t accepted_step{};
  std::uint64_t breakup_ordinal{};
  BreakupChildCandidate candidate{};
  BreakupConservationReport conservation{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == BreakupChildStatus::success && candidate.available;
  }
};

[[nodiscard]] BreakupChildReport generate_supplied_diameter_children(
    const BreakupChildInput& input) noexcept;

}  // namespace hundun::v04::spray::detail
