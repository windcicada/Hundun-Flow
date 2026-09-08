// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_combustion.hpp"

#include <new>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <utility>

namespace hundun::v04::combustion {

MixingTimeReport evaluate_mixing_time(const MixingTimeInput& input) noexcept {
  MixingTimeReport report;
  if (!std::isfinite(input.filter_width_m) ||
      !(input.filter_width_m > 0.0) ||
      !std::isfinite(input.molecular_diffusivity_m2_per_s) ||
      input.molecular_diffusivity_m2_per_s < 0.0 ||
      !std::isfinite(input.turbulent_kinematic_viscosity_m2_per_s) ||
      input.turbulent_kinematic_viscosity_m2_per_s < 0.0 ||
      !std::isfinite(input.turbulent_schmidt) ||
      !(input.turbulent_schmidt > 0.0) || !std::isfinite(input.c_z) ||
      input.c_z < 0.0) {
    return report;
  }

  report.filter_width_m = input.filter_width_m;
  report.molecular_diffusivity_m2_per_s =
      input.molecular_diffusivity_m2_per_s;
  report.turbulent_kinematic_viscosity_m2_per_s =
      input.turbulent_kinematic_viscosity_m2_per_s;
  report.turbulent_schmidt = input.turbulent_schmidt;
  report.c_z = input.c_z;
  report.turbulent_diffusivity_m2_per_s =
      input.turbulent_kinematic_viscosity_m2_per_s /
      input.turbulent_schmidt;
  const double total_diffusivity =
      input.molecular_diffusivity_m2_per_s +
      report.turbulent_diffusivity_m2_per_s;
  if (!std::isfinite(report.turbulent_diffusivity_m2_per_s) ||
      !std::isfinite(total_diffusivity) || !(total_diffusivity > 0.0)) {
    report.status = MixingTimeStatus::timescale_unavailable;
    report.turbulent_diffusivity_m2_per_s = 0.0;
    return report;
  }
  report.tau_mix_s =
      input.c_z * input.filter_width_m * input.filter_width_m /
      (2.0 * total_diffusivity);
  if (!std::isfinite(report.tau_mix_s) || report.tau_mix_s < 0.0) {
    report.status = MixingTimeStatus::timescale_unavailable;
    report.tau_mix_s = 0.0;
    return report;
  }
  report.status = MixingTimeStatus::success;
  return report;
}

namespace {

constexpr double kSimplexTolerance = 1.0e-12;
constexpr double kConservationAbsoluteTolerance = 1.0e-12;
constexpr double kConservationRelativeTolerance = 1.0e-10;
constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

void hash_byte(std::uint64_t& hash, std::uint8_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
}

void hash_u64(std::uint64_t& hash, std::uint64_t value) noexcept {
  for (std::uint32_t shift = 0U; shift < 64U; shift += 8U) {
    hash_byte(hash, static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }
}

void hash_double(std::uint64_t& hash, double value) noexcept {
  std::uint64_t bits{};
  static_assert(sizeof(bits) == sizeof(value), "double/u64 size mismatch");
  std::memcpy(&bits, &value, sizeof(bits));
  hash_u64(hash, bits);
}

bool valid_identity(const ChemistryIdentity& identity) noexcept {
  if (identity.element_count == 0U || identity.species.empty() ||
      identity.fingerprint == 0U ||
      identity.fingerprint != chemistry_identity_fingerprint(identity)) {
    return false;
  }
  for (const ChemicalSpeciesIdentity& species : identity.species) {
    if (!std::isfinite(species.molecular_weight_kg_per_kmol) ||
        !(species.molecular_weight_kg_per_kmol > 0.0) ||
        species.element_counts.size() != identity.element_count ||
        !std::isfinite(species.formation_enthalpy_j_per_kg)) {
      return false;
    }
  }
  return true;
}

bool valid_state(const ThermochemicalState& state,
                 const ChemistryIdentity& identity) noexcept {
  if (!std::isfinite(state.pressure_pa) || !(state.pressure_pa > 0.0) ||
      !std::isfinite(state.density_kg_per_m3) ||
      !(state.density_kg_per_m3 > 0.0) ||
      !std::isfinite(state.total_thermochemical_enthalpy_j_per_kg) ||
      state.mass_fractions.size() != identity.species.size()) {
    return false;
  }
  double sum = 0.0;
  for (const double value : state.mass_fractions) {
    if (!std::isfinite(value) || value < 0.0) {
      return false;
    }
    sum += value;
  }
  return std::isfinite(sum) &&
         std::abs(sum - 1.0) <=
             kSimplexTolerance *
                 static_cast<double>(state.mass_fractions.size());
}

bool conservative_values(const std::vector<double>& values,
                         const ChemistryIdentity& identity) noexcept {
  if (values.size() != identity.species.size()) {
    return false;
  }
  double mass_residual = 0.0;
  double mass_scale = 0.0;
  for (const double value : values) {
    if (!std::isfinite(value)) {
      return false;
    }
    mass_residual += value;
    mass_scale += std::abs(value);
  }
  if (std::abs(mass_residual) >
      kConservationAbsoluteTolerance +
          kConservationRelativeTolerance * mass_scale) {
    return false;
  }
  for (std::size_t element = 0U; element < identity.element_count;
       ++element) {
    double residual = 0.0;
    double scale = 0.0;
    for (std::size_t species = 0U; species < identity.species.size();
         ++species) {
      const double contribution =
          values[species] *
          static_cast<double>(identity.species[species]
                                  .element_counts[element]) /
          identity.species[species].molecular_weight_kg_per_kmol;
      residual += contribution;
      scale += std::abs(contribution);
    }
    if (std::abs(residual) >
        kConservationAbsoluteTolerance +
            kConservationRelativeTolerance * scale) {
      return false;
    }
  }
  return true;
}

bool same_elements(const ThermochemicalState& initial,
                   const ThermochemicalState& final,
                   const ChemistryIdentity& identity) noexcept {
  for (std::size_t element = 0U; element < identity.element_count;
       ++element) {
    double residual = 0.0;
    double scale = 0.0;
    for (std::size_t species = 0U; species < identity.species.size();
         ++species) {
      const double coefficient =
          static_cast<double>(identity.species[species]
                                  .element_counts[element]) /
          identity.species[species].molecular_weight_kg_per_kmol;
      const double before = coefficient * initial.mass_fractions[species];
      const double after = coefficient * final.mass_fractions[species];
      residual += after - before;
      scale += std::abs(after) + std::abs(before);
    }
    if (std::abs(residual) >
        kConservationAbsoluteTolerance +
            kConservationRelativeTolerance * scale) {
      return false;
    }
  }
  return true;
}

bool heat_release_consistent(const ChemistryAdvanceReport& report,
                             const ChemistryIdentity& identity) noexcept {
  double expected = 0.0;
  double scale = 0.0;
  for (std::size_t species = 0U; species < identity.species.size();
       ++species) {
    const double contribution =
        identity.species[species].formation_enthalpy_j_per_kg *
        report.integrated_species_density_delta_kg_per_m3[species];
    expected -= contribution;
    scale += std::abs(contribution);
  }
  return std::isfinite(report.integrated_heat_release_j_per_m3) &&
         std::abs(report.integrated_heat_release_j_per_m3 - expected) <=
             kConservationAbsoluteTolerance +
                 kConservationRelativeTolerance * scale;
}

bool valid_advance_report(const ChemistryAdvanceRequest& request,
                          const ChemistryAdvanceReport& report,
                          const ChemistryIdentity& identity,
                          CombustionClosureStatus& failure) noexcept {
  if (!report.succeeded()) {
    failure = CombustionClosureStatus::backend_failure;
    return false;
  }
  if (!std::isfinite(report.completed_duration_s) ||
      report.completed_duration_s != request.duration_s ||
      (request.duration_s > 0.0 && report.internal_step_count == 0U) ||
      !valid_state(report.final_state, identity) ||
      report.final_state.pressure_pa != request.state.pressure_pa ||
      report.final_state.total_thermochemical_enthalpy_j_per_kg !=
          request.state.total_thermochemical_enthalpy_j_per_kg ||
      report.integrated_species_density_delta_kg_per_m3.size() !=
          identity.species.size()) {
    failure = CombustionClosureStatus::backend_contract_failure;
    return false;
  }
  if (!conservative_values(
          report.integrated_species_density_delta_kg_per_m3, identity) ||
      !same_elements(request.state, report.final_state, identity) ||
      !heat_release_consistent(report, identity)) {
    failure = CombustionClosureStatus::conservation_failure;
    return false;
  }
  return true;
}

struct MeanAdvance {
  CombustionClosureStatus status{CombustionClosureStatus::backend_failure};
  ThermochemicalState final_state;
  std::vector<double> species_delta;
  double heat_release{};
  std::uint32_t call_count{};
  std::uint64_t internal_step_count{};
};

MeanAdvance advance_mean_interval(const ThermochemicalState& initial,
                                  double start_time_s, double duration_s,
                                  const ChemistryIdentity& identity,
                                  ChemistryRepresentationAdapter& chemistry) {
  MeanAdvance result;
  result.species_delta.assign(identity.species.size(), 0.0);
  ThermochemicalState state = initial;
  const double half_duration = 0.5 * duration_s;
  for (std::uint32_t half = 0U; half < 2U; ++half) {
    const double half_start = start_time_s +
                              static_cast<double>(half) * half_duration;
    if (!std::isfinite(half_start)) {
      result.status = CombustionClosureStatus::invalid_input;
      return result;
    }
    const ChemistryAdvanceRequest request{state, half_start, half_duration};
    ChemistryAdvanceReport report = chemistry.integrate(request);
    ++result.call_count;
    CombustionClosureStatus failure{};
    if (!valid_advance_report(request, report, identity, failure)) {
      result.status = failure;
      return result;
    }
    for (std::size_t species = 0U; species < result.species_delta.size();
         ++species) {
      result.species_delta[species] +=
          report.integrated_species_density_delta_kg_per_m3[species];
    }
    result.heat_release += report.integrated_heat_release_j_per_m3;
    result.internal_step_count += report.internal_step_count;
    state = std::move(report.final_state);
  }
  if (!conservative_values(result.species_delta, identity) ||
      !std::isfinite(result.heat_release)) {
    result.status = CombustionClosureStatus::conservation_failure;
    return result;
  }
  result.final_state = std::move(state);
  result.status = CombustionClosureStatus::success;
  return result;
}

ChemicalTimescaleReport evaluate_chemical_timescale(
    const ChemicalRateQuery& query, const ChemistryIdentity& identity,
    ChemicalRateQueryProvider& provider) noexcept {
  ChemicalTimescaleReport result;
  result.model = query.model;
  const ChemicalRateQueryReport provider_report = provider.query(query);
  if (!provider_report.succeeded()) {
    result.status =
        provider_report.status == ChemicalRateQueryStatus::unavailable
            ? ChemicalTimescaleStatus::query_unavailable
            : ChemicalTimescaleStatus::query_failure;
    return result;
  }
  if (provider_report.model != query.model ||
      provider_report.composition_fingerprint !=
          query.composition_fingerprint ||
      provider_report.net_species_mass_rates_kg_per_m3_s.size() !=
          identity.species.size() ||
      !conservative_values(
          provider_report.net_species_mass_rates_kg_per_m3_s, identity)) {
    result.status = ChemicalTimescaleStatus::query_contract_failure;
    return result;
  }
  for (std::size_t species = 0U; species < identity.species.size();
       ++species) {
    const double rate =
        provider_report.net_species_mass_rates_kg_per_m3_s[species];
    if (!std::isfinite(rate)) {
      result.status = ChemicalTimescaleStatus::query_contract_failure;
      return result;
    }
    if (rate < 0.0) {
      result.consuming_species_density_kg_per_m3 +=
          query.state.density_kg_per_m3 * query.state.mass_fractions[species];
      result.consumption_rate_l1_kg_per_m3_s -= rate;
    }
  }
  if (!std::isfinite(result.consuming_species_density_kg_per_m3) ||
      result.consuming_species_density_kg_per_m3 < 0.0 ||
      !std::isfinite(result.consumption_rate_l1_kg_per_m3_s) ||
      !(result.consumption_rate_l1_kg_per_m3_s > 0.0)) {
    result.status = ChemicalTimescaleStatus::timescale_unavailable;
    result.consuming_species_density_kg_per_m3 = 0.0;
    result.consumption_rate_l1_kg_per_m3_s = 0.0;
    return result;
  }
  result.tau_chem_s = result.consuming_species_density_kg_per_m3 /
                      result.consumption_rate_l1_kg_per_m3_s;
  if (!std::isfinite(result.tau_chem_s) || result.tau_chem_s < 0.0) {
    result.status = ChemicalTimescaleStatus::timescale_unavailable;
    result.tau_chem_s = 0.0;
    return result;
  }
  result.status = ChemicalTimescaleStatus::success;
  return result;
}

std::uint64_t closure_identity(const CombustionClosureConfig& config,
                               const ChemistryIdentity& chemistry) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash_u64(hash, UINT64_C(0x48554e44554e4331));
  hash_u64(hash, static_cast<std::uint64_t>(config.tci_closure));
  hash_u64(hash,
           static_cast<std::uint64_t>(config.chemistry_representation));
  hash_u64(hash,
           static_cast<std::uint64_t>(config.chemical_timescale_model));
  hash_u64(hash, chemistry.fingerprint);
  return hash == 0U ? 1U : hash;
}

std::string_view evaluator_model_id(TciClosureKind kind) noexcept {
  switch (kind) {
  case TciClosureKind::finite_rate_mean:
    return "finite_rate_mean_shadow";
  case TciClosureKind::pasr_algebraic_v1:
    return "pasr_algebraic_v1";
  case TciClosureKind::esf_tpdf:
    return "esf_tpdf";
  }
  return "unknown_tci_closure";
}

std::string_view representation_model_id(
    ChemistryRepresentationKind kind) noexcept {
  switch (kind) {
  case ChemistryRepresentationKind::direct_cantera:
    return "direct_cantera";
  case ChemistryRepresentationKind::fgm_table:
    return "fgm_table";
  }
  return "unknown_chemistry_representation";
}

std::string_view timescale_model_id(ChemicalTimescaleModel model) noexcept {
  switch (model) {
  case ChemicalTimescaleModel::reactant_depletion_l1_v1:
    return "reactant_depletion_l1_v1";
  }
  return "unknown_chemical_timescale";
}

} // namespace

std::uint64_t
chemistry_identity_fingerprint(const ChemistryIdentity& identity) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash_u64(hash, static_cast<std::uint64_t>(identity.element_count));
  hash_u64(hash, static_cast<std::uint64_t>(identity.species.size()));
  for (const ChemicalSpeciesIdentity& species : identity.species) {
    hash_double(hash, species.molecular_weight_kg_per_kmol);
    hash_u64(hash, static_cast<std::uint64_t>(species.element_counts.size()));
    for (const std::uint16_t count : species.element_counts) {
      hash_u64(hash, count);
    }
    hash_double(hash, species.formation_enthalpy_j_per_kg);
  }
  return hash == 0U ? 1U : hash;
}

PasrReactingFractionReport evaluate_pasr_reacting_fraction(
    double tau_mix_s, double tau_chem_s) noexcept {
  PasrReactingFractionReport report;
  if (!std::isfinite(tau_mix_s) || tau_mix_s < 0.0 ||
      !std::isfinite(tau_chem_s) || tau_chem_s < 0.0) {
    return report;
  }
  report.tau_mix_s = tau_mix_s;
  report.tau_chem_s = tau_chem_s;
  const double denominator = tau_chem_s + tau_mix_s;
  if (!std::isfinite(denominator) || !(denominator > 0.0)) {
    report.status = PasrReactingFractionStatus::timescale_unavailable;
    return report;
  }
  report.kappa_raw = tau_chem_s / denominator;
  report.kappa = report.kappa_raw;
  if (!std::isfinite(report.kappa) || report.kappa < 0.0 ||
      report.kappa > 1.0) {
    report.kappa_raw = 0.0;
    report.kappa = 0.0;
    report.status = PasrReactingFractionStatus::timescale_unavailable;
    return report;
  }
  report.status = PasrReactingFractionStatus::success;
  return report;
}

CombustionClosureReport evaluate_combustion_closure(
    const CombustionClosureRequest& request,
    const CombustionClosureConfig& config,
    ChemistryRepresentationAdapter& chemistry,
    ChemicalRateQueryProvider* rate_query) noexcept try {
  CombustionClosureReport report;
  report.configured_tci = config.tci_closure;
  report.configured_chemistry = config.chemistry_representation;
  report.evaluator_model_id = evaluator_model_id(config.tci_closure);
  report.chemistry_representation_id =
      representation_model_id(config.chemistry_representation);
  report.chemical_timescale_model_id =
      timescale_model_id(config.chemical_timescale_model);
  report.source_state_revision = request.source_state_revision;
  report.accepted_step = request.accepted_step;

  const ChemistryIdentity& identity = chemistry.identity();
  if (!valid_identity(identity) || request.composition_fingerprint == 0U ||
      !std::isfinite(request.start_time_s) || request.start_time_s < 0.0 ||
      !std::isfinite(request.duration_s) || request.duration_s < 0.0 ||
      !std::isfinite(request.start_time_s + request.duration_s) ||
      !valid_state(request.mean_state, identity)) {
    return report;
  }
  report.candidate.species_density_delta_kg_per_m3.assign(
      identity.species.size(), 0.0);
  if (request.composition_fingerprint != identity.fingerprint) {
    report.status = CombustionClosureStatus::composition_mismatch;
    return report;
  }
  if (chemistry.representation() != config.chemistry_representation) {
    report.status = CombustionClosureStatus::representation_mismatch;
    return report;
  }
  report.model_identity = closure_identity(config, identity);

  if (config.tci_closure == TciClosureKind::finite_rate_mean) {
    const MeanAdvance advance = advance_mean_interval(
        request.mean_state, request.start_time_s, request.duration_s,
        identity, chemistry);
    report.chemistry_call_count = advance.call_count;
    report.chemistry_internal_step_count = advance.internal_step_count;
    report.finite_rate_mean_shadow = true;
    if (advance.status != CombustionClosureStatus::success) {
      report.status = advance.status;
      return report;
    }
    report.candidate.species_density_delta_kg_per_m3 =
        advance.species_delta;
    report.candidate.integrated_heat_release_j_per_m3 =
        advance.heat_release;
    report.candidate.available = true;
    report.status = CombustionClosureStatus::success;
    return report;
  }

  if (config.tci_closure == TciClosureKind::esf_tpdf) {
    const std::size_t field_count = request.stochastic_fields.size();
    if (field_count < 2U || field_count % 2U != 0U ||
        field_count > std::numeric_limits<std::uint32_t>::max()) {
      return report;
    }
    for (const ThermochemicalState& state : request.stochastic_fields) {
      if (!valid_state(state, identity)) {
        return report;
      }
    }
    report.stochastic_field_count =
        static_cast<std::uint32_t>(field_count);
    std::vector<double> mean_delta(identity.species.size(), 0.0);
    double mean_heat_release = 0.0;
    const double weight = 1.0 / static_cast<double>(field_count);
    for (const ThermochemicalState& state : request.stochastic_fields) {
      const MeanAdvance advance = advance_mean_interval(
          state, request.start_time_s, request.duration_s, identity,
          chemistry);
      report.chemistry_call_count += advance.call_count;
      report.chemistry_internal_step_count += advance.internal_step_count;
      if (advance.status != CombustionClosureStatus::success) {
        report.status = advance.status;
        return report;
      }
      for (std::size_t species = 0U; species < mean_delta.size();
           ++species) {
        mean_delta[species] += weight * advance.species_delta[species];
      }
      mean_heat_release += weight * advance.heat_release;
    }
    if (!conservative_values(mean_delta, identity) ||
        !std::isfinite(mean_heat_release)) {
      report.status = CombustionClosureStatus::conservation_failure;
      return report;
    }
    report.candidate.species_density_delta_kg_per_m3 =
        std::move(mean_delta);
    report.candidate.integrated_heat_release_j_per_m3 = mean_heat_release;
    report.candidate.available = true;
    report.status = CombustionClosureStatus::success;
    return report;
  }

  if (config.tci_closure != TciClosureKind::pasr_algebraic_v1) {
    report.status = CombustionClosureStatus::unsupported_portable_model;
    return report;
  }

  const MixingTimeReport mixing = evaluate_mixing_time(config.mixing_time);
  report.timescales.mixing_model = mixing.model;
  report.timescales.chemical_model = config.chemical_timescale_model;
  report.timescales.mixing_status = mixing.status;
  report.timescales.filter_width_m = mixing.filter_width_m;
  report.timescales.molecular_diffusivity_m2_per_s =
      mixing.molecular_diffusivity_m2_per_s;
  report.timescales.turbulent_kinematic_viscosity_m2_per_s =
      mixing.turbulent_kinematic_viscosity_m2_per_s;
  report.timescales.turbulent_schmidt = mixing.turbulent_schmidt;
  report.timescales.turbulent_diffusivity_m2_per_s =
      mixing.turbulent_diffusivity_m2_per_s;
  report.timescales.c_z = mixing.c_z;
  report.timescales.tau_mix_s = mixing.tau_mix_s;
  if (!mixing.succeeded()) {
    report.status = CombustionClosureStatus::mixing_timescale_unavailable;
    return report;
  }
  if (rate_query == nullptr) {
    report.status = CombustionClosureStatus::chemical_timescale_unavailable;
    return report;
  }
  const ChemicalRateQuery query{config.chemical_timescale_model,
                                request.composition_fingerprint,
                                request.mean_state};
  const ChemicalTimescaleReport chemical =
      evaluate_chemical_timescale(query, identity, *rate_query);
  report.timescales.chemical_status = chemical.status;
  report.timescales.tau_chem_s = chemical.tau_chem_s;
  if (!chemical.succeeded()) {
    report.status = CombustionClosureStatus::chemical_timescale_unavailable;
    return report;
  }
  const PasrReactingFractionReport fraction =
      evaluate_pasr_reacting_fraction(mixing.tau_mix_s,
                                      chemical.tau_chem_s);
  if (fraction.status != PasrReactingFractionStatus::success) {
    report.status = CombustionClosureStatus::chemical_timescale_unavailable;
    return report;
  }
  report.timescales.kappa_raw = fraction.kappa_raw;
  report.timescales.kappa = fraction.kappa;
  report.timescales.clamped = fraction.clamped;

  const MeanAdvance advance = advance_mean_interval(
      request.mean_state, request.start_time_s, request.duration_s, identity,
      chemistry);
  report.chemistry_call_count = advance.call_count;
  report.chemistry_internal_step_count = advance.internal_step_count;
  if (advance.status != CombustionClosureStatus::success) {
    report.status = advance.status;
    return report;
  }
  if (report.timescales.kappa == 0.0) {
    report.candidate.species_density_delta_kg_per_m3.assign(
        identity.species.size(), 0.0);
    report.candidate.integrated_heat_release_j_per_m3 = 0.0;
  } else if (report.timescales.kappa == 1.0) {
    report.candidate.species_density_delta_kg_per_m3 =
        advance.species_delta;
    report.candidate.integrated_heat_release_j_per_m3 =
        advance.heat_release;
  } else {
    for (std::size_t species = 0U; species < advance.species_delta.size();
         ++species) {
      report.candidate.species_density_delta_kg_per_m3[species] =
          report.timescales.kappa * advance.species_delta[species];
    }
    report.candidate.integrated_heat_release_j_per_m3 =
        report.timescales.kappa * advance.heat_release;
  }
  report.candidate.integrated_thermochemical_enthalpy_delta_j_per_m3 = 0.0;
  if (!conservative_values(
          report.candidate.species_density_delta_kg_per_m3, identity) ||
      !std::isfinite(report.candidate.integrated_heat_release_j_per_m3)) {
    report.candidate.species_density_delta_kg_per_m3.assign(
        identity.species.size(), 0.0);
    report.candidate.integrated_heat_release_j_per_m3 = 0.0;
    report.status = CombustionClosureStatus::conservation_failure;
    return report;
  }
  report.candidate.available = true;
  report.status = CombustionClosureStatus::success;
  return report;
} catch (const std::bad_alloc&) {
  // The portable value interface owns vectors. Allocation failure is an
  // explicit status, not noexcept termination or a partially formed source.
  CombustionClosureReport report;
  report.status = CombustionClosureStatus::workspace_failure;
  report.configured_tci = config.tci_closure;
  report.configured_chemistry = config.chemistry_representation;
  report.evaluator_model_id = evaluator_model_id(config.tci_closure);
  report.chemistry_representation_id =
      representation_model_id(config.chemistry_representation);
  report.chemical_timescale_model_id =
      timescale_model_id(config.chemical_timescale_model);
  report.model_identity = closure_identity(config, chemistry.identity());
  report.source_state_revision = request.source_state_revision;
  report.accepted_step = request.accepted_step;
  return report;
}

EsfCommonSourceReport apply_esf_common_source(
    const ThermochemicalState& mean,
    const std::vector<ThermochemicalState>& fields,
    const ChemistryIdentity& identity, const EsfCommonSource& source) noexcept try {
  auto failure = [](EsfCommonSourceStatus status) {
    EsfCommonSourceReport report;
    report.status = status;
    return report;
  };
  if (!valid_identity(identity) || !valid_state(mean, identity) ||
      (fields.size() != 2U && fields.size() != 4U) ||
      !std::isfinite(source.density_delta_kg_per_m3) ||
      !std::isfinite(source.thermochemical_enthalpy_density_delta_j_per_m3) ||
      source.species_density_delta_kg_per_m3.size() != identity.species.size())
    return failure(EsfCommonSourceStatus::invalid_input);
  if (source.composition_fingerprint != identity.fingerprint)
    return failure(EsfCommonSourceStatus::composition_mismatch);

  constexpr long double rounding = 256.0L * std::numeric_limits<double>::epsilon();
  const long double count = static_cast<long double>(fields.size());
  long double source_mass = 0.0L, source_mass_scale = 0.0L;
  for (double value : source.species_density_delta_kg_per_m3) {
    if (!std::isfinite(value)) return failure(EsfCommonSourceStatus::invalid_input);
    source_mass += value;
    source_mass_scale += std::abs(static_cast<long double>(value));
  }
  if (std::abs(source_mass - source.density_delta_kg_per_m3) > rounding *
      (source_mass_scale + std::abs(static_cast<long double>(source.density_delta_kg_per_m3))))
    return failure(EsfCommonSourceStatus::source_mass_mismatch);
  for (const auto& field : fields)
    if (!valid_state(field, identity) || field.pressure_pa != mean.pressure_pa ||
        field.density_kg_per_m3 != mean.density_kg_per_m3)
      return failure(EsfCommonSourceStatus::invalid_input);
  for (std::size_t species = 0; species < identity.species.size(); ++species) {
    long double average = 0.0L;
    for (const auto& field : fields) average += field.mass_fractions[species] / count;
    if (std::abs(average - mean.mass_fractions[species]) > rounding)
      return failure(EsfCommonSourceStatus::mean_mismatch);
  }
  long double average_h = 0.0L, h_scale = 0.0L;
  for (const auto& field : fields) {
    const long double h = field.total_thermochemical_enthalpy_j_per_kg;
    average_h += h / count;
    h_scale += std::abs(h) / count;
  }
  if (std::abs(average_h - mean.total_thermochemical_enthalpy_j_per_kg) >
      rounding * (1.0L + h_scale +
        std::abs(static_cast<long double>(mean.total_thermochemical_enthalpy_j_per_kg))))
    return failure(EsfCommonSourceStatus::mean_mismatch);

  const long double rho = mean.density_kg_per_m3;
  const long double new_rho = rho + source.density_delta_kg_per_m3;
  if (new_rho <= 0.0L || !std::isfinite(static_cast<double>(new_rho)))
    return failure(EsfCommonSourceStatus::inadmissible_state);
  auto advance = [&](const ThermochemicalState& initial, ThermochemicalState& candidate) {
    candidate = initial;
    candidate.density_kg_per_m3 = static_cast<double>(new_rho);
    candidate.total_thermochemical_enthalpy_j_per_kg = static_cast<double>(
        (rho * initial.total_thermochemical_enthalpy_j_per_kg +
         source.thermochemical_enthalpy_density_delta_j_per_m3) / new_rho);
    for (std::size_t species = 0; species < identity.species.size(); ++species) {
      const long double inventory = rho * initial.mass_fractions[species] +
          source.species_density_delta_kg_per_m3[species];
      const double value = static_cast<double>(inventory / new_rho);
      if (!std::isfinite(value) || value < 0.0 || value > 1.0) return false;
      candidate.mass_fractions[species] = value;
    }
    return valid_state(candidate, identity);
  };
  EsfCommonSourceReport report;
  report.fields.resize(fields.size());
  if (!advance(mean, report.mean_state))
    return failure(EsfCommonSourceStatus::inadmissible_state);
  for (std::size_t i = 0; i < fields.size(); ++i)
    if (!advance(fields[i], report.fields[i]))
      return failure(EsfCommonSourceStatus::inadmissible_state);
  report.status = EsfCommonSourceStatus::success;
  report.available = true;
  return report;
} catch (const std::bad_alloc&) {
  EsfCommonSourceReport report;
  report.status = EsfCommonSourceStatus::workspace_failure;
  return report;
}

} // namespace hundun::v04::combustion
