// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#include "models_spray_sgs_step_detail.hpp"
#include <cmath>

namespace hundun::v04::spray::detail {
namespace {
ParcelRandomAddress sgs_address(ParcelId id, std::uint64_t seed,
    std::uint64_t step, std::uint64_t ordinal, ParcelRandomPurpose purpose) noexcept {
  return {seed ^ id.high, step, id.low, ordinal, purpose};
}
} // namespace

PersistentSgsBreakupState initialize_sgs_lineage(
    ParcelId id, std::uint64_t seed, std::uint64_t birth_step) noexcept {
  if (id == ParcelId{}) return {};
  const double uniform = parcel_uniform_01(
      sgs_address(id, seed, birth_step, 0U, ParcelRandomPurpose::breakup_lifetime), 0U);
  // Inverse CDF for min(Poisson(1),7). The last bin owns the complete upper
  // tail. A single counter draw gives partition- and retry-stable lifetimes.
  double probability = std::exp(-1.0), cumulative = probability;
  std::uint8_t multiplier = 0U;
  while (uniform >= cumulative && multiplier < 7U) {
    ++multiplier;
    probability /= multiplier;
    cumulative += probability;
  }
  return {1U,{0.0,0.0,0.0,0.0,multiplier}};
}

SgsParcelStepReport advance_sgs_parcel_step(const SgsParcelStepInput& input) noexcept {
  const auto fail = [&](SgsBreakupStatus status) {
    SgsParcelStepReport out;
    out.status = status;
    out.revision = input.revision;
    return out;
  };
  const auto& retained = input.retained;
  const auto& parcel = retained.parcel;
  const auto& env = input.environment;
  if (input.revision.algorithm_version != 1U || !input.revision.input_revision ||
      input.revision.accepted_step == UINT64_MAX ||
      retained.sgs.version != 1U || !valid_sgs_state(retained.sgs) ||
      validate_parcel_state(parcel) != ParcelStateStatus::success ||
      !std::isfinite(retained.tab_deformation) ||
      !std::isfinite(retained.tab_deformation_rate_per_s) ||
      !std::isfinite(env.liquid_absolute_enthalpy_j_per_kg))
    return fail(SgsBreakupStatus::invalid_input);
  double relative_speed = 0.0;
  for (std::size_t axis=0;axis<3U;++axis) {
    if (!std::isfinite(env.gas_velocity_m_per_s[axis]))
      return fail(SgsBreakupStatus::invalid_input);
    relative_speed = std::hypot(relative_speed,
        env.gas_velocity_m_per_s[axis]-parcel.velocity_m_per_s[axis]);
  }
  SgsBreakupInput kernel;
  kernel.history = retained.sgs.history;
  kernel.droplet_diameter_m = parcel.droplet_diameter_m;
  kernel.relative_speed_m_per_s = relative_speed;
  kernel.gas_density_kg_per_m3 = env.gas_density_kg_per_m3;
  kernel.liquid_density_kg_per_m3 = env.liquid_density_kg_per_m3;
  kernel.surface_tension_n_per_m = env.surface_tension_n_per_m;
  kernel.gas_dynamic_viscosity_pa_s = env.gas_dynamic_viscosity_pa_s;
  kernel.dissipation_m2_per_s3 = env.dissipation_m2_per_s3;
  kernel.duration_s = input.duration_s;
  kernel.deterministic_time_coefficient = env.deterministic_time_coefficient;
  kernel.stochastic_coefficient = env.stochastic_coefficient;
  kernel.stochastic_enabled = env.stochastic_enabled;
  kernel.daughter_uniform_01 = parcel_uniform_01(sgs_address(parcel.id,input.seed,
      input.revision.accepted_step,retained.breakup_ordinal,
      ParcelRandomPurpose::breakup_quantile),0U);
  const auto evolved = evaluate_sgs_breakup(kernel);
  if (!evolved.succeeded()) return fail(evolved.status);
  SgsParcelStepReport out;
  out.revision = input.revision;
  out.evolution = evolved;
  if (evolved.breakup_requested) {
    SgsChildInput split_input;
    split_input.parent = parcel;
    split_input.trigger = evolved;
    split_input.accepted_step = input.revision.accepted_step;
    split_input.breakup_ordinal = retained.breakup_ordinal;
    split_input.liquid_density_kg_per_m3 = env.liquid_density_kg_per_m3;
    split_input.surface_tension_n_per_m = env.surface_tension_n_per_m;
    split_input.liquid_absolute_thermochemical_enthalpy_j_per_kg =
        env.liquid_absolute_enthalpy_j_per_kg;
    const auto split = generate_sgs_children(split_input);
    if (!split.succeeded()) return fail(SgsBreakupStatus::numerical_failure);
    out.parent_removed = true;
    out.candidate_count = 2U;
    out.breakup_budget = split.conservation;
    for (std::size_t i=0;i<2U;++i) {
      auto& child = out.candidates[i];
      child.parcel = split.candidate.children[i];
      child.sgs = initialize_sgs_lineage(child.parcel.id,input.seed,
                                        input.revision.accepted_step+1U);
    }
  } else {
    out.candidate_count = 1U;
    out.candidates[0] = retained;
    out.candidates[0].sgs.history = evolved.candidate;
  }
  out.status = SgsBreakupStatus::success;
  out.available = true;
  return out;
}
} // namespace hundun::v04::spray::detail
