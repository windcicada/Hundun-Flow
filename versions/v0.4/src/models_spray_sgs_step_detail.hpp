// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "models_spray_breakup_detail.hpp"
#include "models_spray_migration_detail.hpp"

namespace hundun::v04::spray::detail {
// Birth initialization is explicit. Imported histories enter through their
// versioned restoration path, with their physical units resolved by import.
[[nodiscard]] PersistentSgsBreakupState initialize_sgs_lineage(
    ParcelId id, std::uint64_t seed, std::uint64_t birth_step) noexcept;

struct SgsStepEnvironment {
  Vector3 gas_velocity_m_per_s{};
  double gas_density_kg_per_m3{};
  double gas_dynamic_viscosity_pa_s{};
  double dissipation_m2_per_s3{};
  double liquid_density_kg_per_m3{};
  double surface_tension_n_per_m{};
  double liquid_absolute_enthalpy_j_per_kg{};
  double deterministic_time_coefficient{0.57735026918962576451};
  double stochastic_coefficient{2.0};
  bool stochastic_enabled{true};
};
struct SgsParcelStepInput {
  // Parcel after the retained trajectory/evaporation step; history still at
  // its accepted step origin. Invoke once per fluid step, after substep trials.
  ParcelMigrationValue retained{};
  portable::Revision revision{};
  std::uint64_t seed{};
  double duration_s{};
  SgsStepEnvironment environment{};
};
struct SgsParcelStepReport {
  bool available{};
  SgsBreakupStatus status{SgsBreakupStatus::invalid_input};
  std::string_view model_id{"sgs_post_evaporation_step_v1"};
  portable::Revision revision{};
  bool parent_removed{};
  std::uint32_t candidate_count{};
  std::array<ParcelMigrationValue,2U> candidates{};
  SgsBreakupReport evolution{};
  BreakupConservationReport breakup_budget{};
};
[[nodiscard]] SgsParcelStepReport advance_sgs_parcel_step(
    const SgsParcelStepInput&) noexcept;
} // namespace hundun::v04::spray::detail
