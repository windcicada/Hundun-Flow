// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_combustion.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04::combustion;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool test_scalar_dissipation_mixing_time() {
  const MixingTimeInput input{
      0.02,   // filter_width_m
      1.0e-5, // molecular_diffusivity_m2_per_s
      3.0e-5, // turbulent_kinematic_viscosity_m2_per_s
      0.75,   // turbulent_schmidt
      0.5};   // c_z

  const MixingTimeReport report = evaluate_mixing_time(input);
  return expect(report.status == MixingTimeStatus::success,
                "valid scalar-dissipation inputs succeed") &&
         expect(report.model ==
                    MixingTimeModel::les_scalar_dissipation_v1,
                "mixing-time model identity is explicit") &&
         expect(report.turbulent_diffusivity_m2_per_s == 4.0e-5,
                "D_t is nu_t/Sc_t") &&
         expect(std::abs(report.tau_mix_s - 2.0) < 1.0e-14,
                "tau_mix=C_Z*Delta^2/(2*(D+D_t))") &&
         expect(report.c_z == input.c_z,
                "C_Z is preserved in the report");
}

bool test_pasr_fraction_exact_limits() {
  const PasrReactingFractionReport off =
      evaluate_pasr_reacting_fraction(2.0, 0.0);
  const PasrReactingFractionReport finite_rate =
      evaluate_pasr_reacting_fraction(0.0, 3.0);
  const PasrReactingFractionReport interior =
      evaluate_pasr_reacting_fraction(2.0, 3.0);
  return expect(off.status == PasrReactingFractionStatus::success &&
                    off.kappa == 0.0,
                "tau_chem=0 is the exact kappa=0 limit") &&
         expect(finite_rate.status ==
                        PasrReactingFractionStatus::success &&
                    finite_rate.kappa == 1.0,
                "tau_mix=0 is the exact finite-rate kappa=1 limit") &&
         expect(interior.status ==
                        PasrReactingFractionStatus::success &&
                    interior.kappa == 0.6 && interior.kappa >= 0.0 &&
                    interior.kappa <= 1.0,
                "valid PaSR fraction is bounded without a hidden fallback");
}

ChemistryIdentity chemistry_identity() {
  ChemistryIdentity identity;
  identity.element_count = 2U;
  identity.species = {
      {2.0, {1U, 1U}, 100.0},
      {1.0, {1U, 0U}, 0.0},
      {1.0, {0U, 1U}, 0.0}};
  identity.fingerprint = chemistry_identity_fingerprint(identity);
  return identity;
}

ChemistryIdentity chemistry_identity_with_isomers() {
  ChemistryIdentity identity;
  identity.element_count = 2U;
  identity.species = {
      {2.0, {1U, 1U}, 100.0},
      {1.0, {1U, 0U}, 0.0},
      {1.0, {0U, 1U}, 0.0},
      {1.0, {1U, 0U}, 0.0},
      {1.0, {1U, 0U}, 0.0}};
  identity.fingerprint = chemistry_identity_fingerprint(identity);
  return identity;
}

enum class FakeChemistryMode {
  normal,
  fail_second_call,
  inconsistent_heat_release
};

class FakeChemistry final : public ChemistryRepresentationAdapter {
public:
  explicit FakeChemistry(
      ChemistryIdentity identity,
      FakeChemistryMode mode = FakeChemistryMode::normal)
      : identity_(std::move(identity)), mode_(mode) {}

  const ChemistryIdentity& identity() const noexcept override {
    return identity_;
  }

  ChemistryRepresentationKind representation() const noexcept override {
    return ChemistryRepresentationKind::direct_cantera;
  }

  ChemistryAdvanceReport
  integrate(const ChemistryAdvanceRequest& request) noexcept override {
    ++call_count;
    ChemistryAdvanceReport report;
    if (mode_ == FakeChemistryMode::fail_second_call && call_count == 2U) {
      report.status = ChemistryAdvanceStatus::integration_failure;
      report.final_state = request.state;
      report.integrated_species_density_delta_kg_per_m3.assign(
          identity_.species.size(), 99.0);
      report.integrated_heat_release_j_per_m3 = 99.0;
      return report;
    }
    report.status = ChemistryAdvanceStatus::success;
    report.final_state = request.state;
    report.final_state.mass_fractions[0U] -= 0.1;
    report.final_state.mass_fractions[1U] += 0.05;
    report.final_state.mass_fractions[2U] += 0.05;
    report.integrated_species_density_delta_kg_per_m3.assign(
        identity_.species.size(), 0.0);
    report.integrated_species_density_delta_kg_per_m3[0U] = -0.1;
    report.integrated_species_density_delta_kg_per_m3[1U] = 0.05;
    report.integrated_species_density_delta_kg_per_m3[2U] = 0.05;
    report.integrated_heat_release_j_per_m3 =
        mode_ == FakeChemistryMode::inconsistent_heat_release ? 9.0 : 10.0;
    report.completed_duration_s = request.duration_s;
    report.internal_step_count = 3U;
    return report;
  }

  std::size_t call_count{};

private:
  ChemistryIdentity identity_;
  FakeChemistryMode mode_;
};

class FakeRateQuery final : public ChemicalRateQueryProvider {
public:
  explicit FakeRateQuery(std::uint64_t fingerprint)
      : fingerprint_(fingerprint) {}

  ChemicalRateQueryReport
  query(const ChemicalRateQuery& request) noexcept override {
    ++call_count;
    ChemicalRateQueryReport report;
    report.status = ChemicalRateQueryStatus::success;
    report.model = request.model;
    report.composition_fingerprint = fingerprint_;
    report.net_species_mass_rates_kg_per_m3_s = {-0.2, 0.1, 0.1};
    return report;
  }

  std::size_t call_count{};

private:
  std::uint64_t fingerprint_{};
};

class ZeroInventoryRateQuery final : public ChemicalRateQueryProvider {
public:
  explicit ZeroInventoryRateQuery(std::uint64_t fingerprint)
      : fingerprint_(fingerprint) {}

  ChemicalRateQueryReport
  query(const ChemicalRateQuery& request) noexcept override {
    ++call_count;
    return {ChemicalRateQueryStatus::success,
            request.model,
            fingerprint_,
            {0.0, 0.0, 0.0, -0.1, 0.1}};
  }

  std::size_t call_count{};

private:
  std::uint64_t fingerprint_{};
};

class UnavailableRateQuery final : public ChemicalRateQueryProvider {
public:
  ChemicalRateQueryReport
  query(const ChemicalRateQuery& request) noexcept override {
    ++call_count;
    return {ChemicalRateQueryStatus::unavailable,
            request.model,
            request.composition_fingerprint,
            {}};
  }

  std::size_t call_count{};
};

bool near(double left, double right, double tolerance = 1.0e-14) {
  return std::abs(left - right) <= tolerance;
}

bool test_pasr_scales_one_mean_advance_consistently() {
  const ChemistryIdentity identity = chemistry_identity();
  FakeChemistry backend(identity);
  FakeRateQuery rate_query(identity.fingerprint);

  CombustionClosureRequest request;
  request.composition_fingerprint = identity.fingerprint;
  request.source_state_revision = 41U;
  request.accepted_step = 9U;
  request.start_time_s = 1.25;
  request.duration_s = 1.0;
  request.mean_state =
      {101325.0, 1.0, 420000.0, {0.6, 0.2, 0.2}};

  CombustionClosureConfig config;
  config.tci_closure = TciClosureKind::pasr_algebraic_v1;
  config.chemistry_representation =
      ChemistryRepresentationKind::direct_cantera;
  config.chemical_timescale_model =
      ChemicalTimescaleModel::reactant_depletion_l1_v1;
  config.mixing_time = {0.02, 1.0e-5, 3.0e-5, 0.75, 0.5};

  const CombustionClosureReport report = evaluate_combustion_closure(
      request, config, backend, &rate_query);

  bool passed = true;
  passed &= expect(report.status == CombustionClosureStatus::success,
                   "PaSR closure succeeds through the portable seam");
  passed &= expect(report.configured_tci ==
                       TciClosureKind::pasr_algebraic_v1 &&
                       report.configured_chemistry ==
                           ChemistryRepresentationKind::direct_cantera,
                   "TCI and chemistry representation remain orthogonal");
  passed &= expect(report.evaluator_model_id == "pasr_algebraic_v1" &&
                       report.chemistry_representation_id ==
                           "direct_cantera" &&
                       report.chemical_timescale_model_id ==
                           "reactant_depletion_l1_v1",
                   "PaSR report carries stable readable model identities");
  passed &= expect(backend.call_count == 2U &&
                       report.chemistry_call_count == 2U,
                   "one mean interval uses two chemistry half advances");
  passed &= expect(rate_query.call_count == 1U,
                   "PaSR makes one named chemical-rate query");
  passed &= expect(near(report.timescales.tau_mix_s, 2.0) &&
                       near(report.timescales.tau_chem_s, 3.0) &&
                       near(report.timescales.kappa, 0.6),
                   "HUNDUN reports explicit tau_mix, tau_chem and kappa");
  passed &= expect(report.candidate.available &&
                       report.candidate.species_density_delta_kg_per_m3.size() ==
                           3U,
                   "successful closure returns one candidate value");
  passed &= expect(
      near(report.candidate.species_density_delta_kg_per_m3[0U], -0.12) &&
          near(report.candidate.species_density_delta_kg_per_m3[1U], 0.06) &&
          near(report.candidate.species_density_delta_kg_per_m3[2U], 0.06),
      "one kappa scales every reacting-species increment");
  passed &= expect(
      near(report.candidate.integrated_heat_release_j_per_m3, 12.0) &&
          report.candidate.integrated_thermochemical_enthalpy_delta_j_per_m3 ==
              0.0,
      "the same kappa scales heat release while h_tc remains conserved");
  passed &= expect(!report.source_published &&
                       report.source_state_revision == 41U &&
                       report.accepted_step == 9U,
                   "portable closure reports provenance and never publishes state");
  return passed;
}

bool canonical_zero_candidate(const CombustionCandidate& candidate,
                              std::size_t species_count) {
  if (candidate.available ||
      candidate.species_density_delta_kg_per_m3.size() != species_count ||
      candidate.integrated_thermochemical_enthalpy_delta_j_per_m3 != 0.0 ||
      candidate.integrated_heat_release_j_per_m3 != 0.0) {
    return false;
  }
  for (const double value :
       candidate.species_density_delta_kg_per_m3) {
    if (value != 0.0 || std::signbit(value)) {
      return false;
    }
  }
  return true;
}

CombustionClosureRequest closure_request(const ChemistryIdentity& identity);

CombustionClosureConfig pasr_config() {
  CombustionClosureConfig config;
  config.tci_closure = TciClosureKind::pasr_algebraic_v1;
  config.chemistry_representation =
      ChemistryRepresentationKind::direct_cantera;
  config.chemical_timescale_model =
      ChemicalTimescaleModel::reactant_depletion_l1_v1;
  config.mixing_time = {0.02, 1.0e-5, 3.0e-5, 0.75, 0.5};
  return config;
}

bool test_pasr_candidate_exact_limits() {
  bool passed = true;
  const ChemistryIdentity finite_identity = chemistry_identity();
  CombustionClosureRequest finite_request = closure_request(finite_identity);
  CombustionClosureConfig finite_config = pasr_config();
  finite_config.mixing_time.c_z = 0.0;
  FakeChemistry finite_backend(finite_identity);
  FakeRateQuery finite_rates(finite_identity.fingerprint);
  const CombustionClosureReport finite_rate = evaluate_combustion_closure(
      finite_request, finite_config, finite_backend, &finite_rates);
  passed &= expect(finite_rate.succeeded() &&
                       finite_rate.timescales.kappa == 1.0 &&
                       finite_rate.candidate
                               .species_density_delta_kg_per_m3[0U] ==
                           -0.2 &&
                       finite_rate.candidate
                               .species_density_delta_kg_per_m3[1U] ==
                           0.1 &&
                       finite_rate.candidate
                               .species_density_delta_kg_per_m3[2U] ==
                           0.1 &&
                       finite_rate.candidate
                               .integrated_heat_release_j_per_m3 ==
                           20.0,
                   "kappa=1 is exactly the finite-rate candidate");

  const ChemistryIdentity off_identity = chemistry_identity_with_isomers();
  CombustionClosureRequest off_request;
  off_request.composition_fingerprint = off_identity.fingerprint;
  off_request.source_state_revision = 18U;
  off_request.accepted_step = 6U;
  off_request.start_time_s = 0.5;
  off_request.duration_s = 1.0;
  off_request.mean_state =
      {101325.0, 1.0, 420000.0, {0.5, 0.15, 0.25, 0.0, 0.1}};
  FakeChemistry off_backend(off_identity);
  ZeroInventoryRateQuery off_rates(off_identity.fingerprint);
  const CombustionClosureReport off = evaluate_combustion_closure(
      off_request, pasr_config(), off_backend, &off_rates);
  passed &= expect(off.succeeded() && off.timescales.kappa == 0.0 &&
                       off_backend.call_count == 2U &&
                       off.candidate.available,
                   "kappa=0 still performs one mean chemistry advance");
  for (const double value : off.candidate.species_density_delta_kg_per_m3) {
    passed &= expect(value == 0.0 && !std::signbit(value),
                     "kappa=0 gives canonical zero for every species");
  }
  passed &= expect(off.candidate.integrated_heat_release_j_per_m3 == 0.0 &&
                       !std::signbit(
                           off.candidate.integrated_heat_release_j_per_m3),
                   "kappa=0 gives canonical zero heat release");
  return passed;
}

bool test_failures_never_expose_partial_candidate() {
  const ChemistryIdentity identity = chemistry_identity();
  const CombustionClosureRequest request = closure_request(identity);
  bool passed = true;

  CombustionClosureConfig no_diffusivity = pasr_config();
  no_diffusivity.mixing_time.molecular_diffusivity_m2_per_s = 0.0;
  no_diffusivity.mixing_time.turbulent_kinematic_viscosity_m2_per_s = 0.0;
  FakeChemistry mixing_backend(identity);
  FakeRateQuery unused_rates(identity.fingerprint);
  const CombustionClosureReport mixing_failure = evaluate_combustion_closure(
      request, no_diffusivity, mixing_backend, &unused_rates);
  passed &= expect(
      mixing_failure.status ==
              CombustionClosureStatus::mixing_timescale_unavailable &&
          mixing_failure.timescales.mixing_status ==
              MixingTimeStatus::timescale_unavailable &&
          mixing_backend.call_count == 0U && unused_rates.call_count == 0U &&
          canonical_zero_candidate(mixing_failure.candidate, 3U),
      "unavailable mixing time is explicit and fails before chemistry");

  FakeChemistry query_backend(identity);
  UnavailableRateQuery unavailable_rates;
  const CombustionClosureReport query_failure = evaluate_combustion_closure(
      request, pasr_config(), query_backend, &unavailable_rates);
  passed &= expect(
      query_failure.status ==
              CombustionClosureStatus::chemical_timescale_unavailable &&
          query_failure.timescales.chemical_status ==
              ChemicalTimescaleStatus::query_unavailable &&
          query_backend.call_count == 0U && unavailable_rates.call_count == 1U &&
          canonical_zero_candidate(query_failure.candidate, 3U),
      "unavailable chemical time is explicit and has no model fallback");

  FakeChemistry failing_backend(identity,
                                FakeChemistryMode::fail_second_call);
  FakeRateQuery valid_rates(identity.fingerprint);
  const CombustionClosureReport backend_failure =
      evaluate_combustion_closure(request, pasr_config(), failing_backend,
                                  &valid_rates);
  passed &= expect(
      backend_failure.status == CombustionClosureStatus::backend_failure &&
          failing_backend.call_count == 2U &&
          canonical_zero_candidate(backend_failure.candidate, 3U) &&
          !backend_failure.source_published,
      "second-half backend failure discards the first-half candidate");

  FakeChemistry inconsistent_energy(
      identity, FakeChemistryMode::inconsistent_heat_release);
  FakeRateQuery energy_rates(identity.fingerprint);
  const CombustionClosureReport energy_failure = evaluate_combustion_closure(
      request, pasr_config(), inconsistent_energy, &energy_rates);
  passed &= expect(
      energy_failure.status == CombustionClosureStatus::conservation_failure &&
          inconsistent_energy.call_count == 1U &&
          canonical_zero_candidate(energy_failure.candidate, 3U),
      "species-inconsistent heat release fails atomically");
  return passed;
}

CombustionClosureRequest closure_request(const ChemistryIdentity& identity) {
  CombustionClosureRequest request;
  request.composition_fingerprint = identity.fingerprint;
  request.source_state_revision = 17U;
  request.accepted_step = 5U;
  request.start_time_s = 0.5;
  request.duration_s = 1.0;
  request.mean_state =
      {101325.0, 1.0, 420000.0, {0.6, 0.2, 0.2}};
  return request;
}

bool test_portable_call_topology_and_model_identity() {
  const ChemistryIdentity identity = chemistry_identity();
  CombustionClosureConfig mean_config;
  mean_config.tci_closure = TciClosureKind::finite_rate_mean;
  mean_config.chemistry_representation =
      ChemistryRepresentationKind::direct_cantera;

  FakeChemistry first_mean_backend(identity);
  const CombustionClosureReport first_mean = evaluate_combustion_closure(
      closure_request(identity), mean_config, first_mean_backend, nullptr);
  FakeChemistry second_mean_backend(identity);
  const CombustionClosureReport second_mean = evaluate_combustion_closure(
      closure_request(identity), mean_config, second_mean_backend, nullptr);

  bool passed = true;
  passed &= expect(first_mean.succeeded() &&
                       first_mean.finite_rate_mean_shadow &&
                       first_mean.evaluator_model_id ==
                           "finite_rate_mean_shadow" &&
                       first_mean.chemistry_call_count == 2U &&
                       first_mean_backend.call_count == 2U,
                   "finite_rate_mean_shadow uses two half calls");
  passed &= expect(first_mean.model_identity != 0U &&
                       first_mean.model_identity == second_mean.model_identity,
                   "equal mean-model reports have deterministic identity");
  passed &= expect(
      near(first_mean.candidate.species_density_delta_kg_per_m3[0U], -0.2) &&
          near(first_mean.candidate.species_density_delta_kg_per_m3[1U], 0.1) &&
          near(first_mean.candidate.species_density_delta_kg_per_m3[2U], 0.1) &&
          near(first_mean.candidate.integrated_heat_release_j_per_m3, 20.0) &&
          !first_mean.source_published,
      "mean shadow returns the unscaled conservative candidate only");

  CombustionClosureRequest esf_request = closure_request(identity);
  esf_request.stochastic_fields.assign(4U, esf_request.mean_state);
  CombustionClosureConfig esf_config = mean_config;
  esf_config.tci_closure = TciClosureKind::esf_tpdf;
  FakeChemistry esf_backend(identity);
  const CombustionClosureReport esf = evaluate_combustion_closure(
      esf_request, esf_config, esf_backend, nullptr);

  passed &= expect(esf.succeeded() && esf.stochastic_field_count == 4U &&
                       esf.chemistry_call_count == 8U &&
                       esf_backend.call_count == 8U,
                   "N=4 TPDF chemistry topology is exactly 2N");
  passed &= expect(esf.model_identity != first_mean.model_identity &&
                       esf.evaluator_model_id == "esf_tpdf" &&
                       !esf.finite_rate_mean_shadow &&
                       !esf.source_published,
                   "TPDF identity is distinct and remains candidate-only");
  passed &= expect(
      near(esf.candidate.species_density_delta_kg_per_m3[0U], -0.2) &&
          near(esf.candidate.species_density_delta_kg_per_m3[1U], 0.1) &&
          near(esf.candidate.species_density_delta_kg_per_m3[2U], 0.1) &&
          near(esf.candidate.integrated_heat_release_j_per_m3, 20.0),
      "TPDF chemistry-only ensemble uses deterministic equal weights");
  return passed;
}

} // namespace

int main() {
  bool passed = test_scalar_dissipation_mixing_time();
  passed &= test_pasr_fraction_exact_limits();
  passed &= test_pasr_scales_one_mean_advance_consistently();
  passed &= test_pasr_candidate_exact_limits();
  passed &= test_failures_never_expose_partial_candidate();
  passed &= test_portable_call_topology_and_model_identity();
  return passed ? 0 : 1;
}
