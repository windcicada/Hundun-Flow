// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include <cstdint>
#include <cstddef>
#include <string_view>
#include <vector>

namespace hundun::v04::combustion {

enum class MixingTimeModel : std::uint8_t {
  les_scalar_dissipation_v1
};

enum class MixingTimeStatus : std::uint8_t {
  success,
  invalid_input,
  timescale_unavailable
};

struct MixingTimeInput {
  double filter_width_m{};
  double molecular_diffusivity_m2_per_s{};
  double turbulent_kinematic_viscosity_m2_per_s{};
  double turbulent_schmidt{};
  double c_z{};
};

struct MixingTimeReport {
  MixingTimeModel model{MixingTimeModel::les_scalar_dissipation_v1};
  MixingTimeStatus status{MixingTimeStatus::invalid_input};
  double filter_width_m{};
  double molecular_diffusivity_m2_per_s{};
  double turbulent_kinematic_viscosity_m2_per_s{};
  double turbulent_schmidt{};
  double turbulent_diffusivity_m2_per_s{};
  double c_z{};
  double tau_mix_s{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == MixingTimeStatus::success;
  }
};

[[nodiscard]] MixingTimeReport
evaluate_mixing_time(const MixingTimeInput& input) noexcept;

enum class TciClosureKind : std::uint8_t {
  finite_rate_mean,
  pasr_algebraic_v1,
  esf_tpdf
};

enum class ChemistryRepresentationKind : std::uint8_t {
  direct_cantera,
  fgm_table
};

enum class ChemicalTimescaleModel : std::uint8_t {
  // tau = sum(rho*Y_i for omega_i < 0) / sum(-omega_i for omega_i < 0).
  // The query adapter supplies omega at the exact requested state; HUNDUN
  // owns this reduction, its validity rules, and the versioned identity.
  reactant_depletion_l1_v1
};

struct ChemicalSpeciesIdentity {
  double molecular_weight_kg_per_kmol{};
  std::vector<std::uint16_t> element_counts;
  double formation_enthalpy_j_per_kg{};
};

struct ChemistryIdentity {
  std::size_t element_count{};
  std::vector<ChemicalSpeciesIdentity> species;
  std::uint64_t fingerprint{};
};

[[nodiscard]] std::uint64_t
chemistry_identity_fingerprint(const ChemistryIdentity& identity) noexcept;

struct ThermochemicalState {
  double pressure_pa{};
  double density_kg_per_m3{};
  double total_thermochemical_enthalpy_j_per_kg{};
  std::vector<double> mass_fractions;
};

struct ChemistryAdvanceRequest {
  ThermochemicalState state;
  double start_time_s{};
  double duration_s{};
};

enum class ChemistryAdvanceStatus : std::uint8_t {
  success,
  invalid_input,
  integration_failure,
  workspace_failure,
  non_finite_output
};

struct ChemistryAdvanceReport {
  ChemistryAdvanceStatus status{ChemistryAdvanceStatus::invalid_input};
  ThermochemicalState final_state;
  std::vector<double> integrated_species_density_delta_kg_per_m3;
  double integrated_heat_release_j_per_m3{};
  double completed_duration_s{};
  std::uint32_t internal_step_count{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == ChemistryAdvanceStatus::success;
  }
};

class ChemistryRepresentationAdapter {
public:
  virtual ~ChemistryRepresentationAdapter() = default;
  virtual const ChemistryIdentity& identity() const noexcept = 0;
  virtual ChemistryRepresentationKind representation() const noexcept = 0;
  virtual ChemistryAdvanceReport
  integrate(const ChemistryAdvanceRequest& request) noexcept = 0;
};

struct ChemicalRateQuery {
  ChemicalTimescaleModel model{
      ChemicalTimescaleModel::reactant_depletion_l1_v1};
  std::uint64_t composition_fingerprint{};
  ThermochemicalState state;
};

enum class ChemicalRateQueryStatus : std::uint8_t {
  success,
  unavailable,
  invalid_input,
  provider_failure
};

struct ChemicalRateQueryReport {
  ChemicalRateQueryStatus status{ChemicalRateQueryStatus::unavailable};
  ChemicalTimescaleModel model{
      ChemicalTimescaleModel::reactant_depletion_l1_v1};
  std::uint64_t composition_fingerprint{};
  std::vector<double> net_species_mass_rates_kg_per_m3_s;

  [[nodiscard]] bool succeeded() const noexcept {
    return status == ChemicalRateQueryStatus::success;
  }
};

class ChemicalRateQueryProvider {
public:
  virtual ~ChemicalRateQueryProvider() = default;
  virtual ChemicalRateQueryReport
  query(const ChemicalRateQuery& request) noexcept = 0;
};

enum class ChemicalTimescaleStatus : std::uint8_t {
  success,
  query_unavailable,
  query_failure,
  query_contract_failure,
  timescale_unavailable
};

struct ChemicalTimescaleReport {
  ChemicalTimescaleModel model{
      ChemicalTimescaleModel::reactant_depletion_l1_v1};
  ChemicalTimescaleStatus status{
      ChemicalTimescaleStatus::query_unavailable};
  double consuming_species_density_kg_per_m3{};
  double consumption_rate_l1_kg_per_m3_s{};
  double tau_chem_s{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == ChemicalTimescaleStatus::success;
  }
};

struct PasrTimescaleReport {
  MixingTimeModel mixing_model{
      MixingTimeModel::les_scalar_dissipation_v1};
  ChemicalTimescaleModel chemical_model{
      ChemicalTimescaleModel::reactant_depletion_l1_v1};
  MixingTimeStatus mixing_status{MixingTimeStatus::invalid_input};
  ChemicalTimescaleStatus chemical_status{
      ChemicalTimescaleStatus::query_unavailable};
  double filter_width_m{};
  double molecular_diffusivity_m2_per_s{};
  double turbulent_kinematic_viscosity_m2_per_s{};
  double turbulent_schmidt{};
  double turbulent_diffusivity_m2_per_s{};
  double c_z{};
  double tau_mix_s{};
  double tau_chem_s{};
  double kappa_raw{};
  double kappa{};
  bool clamped{};
};

enum class PasrReactingFractionStatus : std::uint8_t {
  success,
  invalid_input,
  timescale_unavailable
};

struct PasrReactingFractionReport {
  PasrReactingFractionStatus status{
      PasrReactingFractionStatus::invalid_input};
  double tau_mix_s{};
  double tau_chem_s{};
  double kappa_raw{};
  double kappa{};
  bool clamped{};
};

[[nodiscard]] PasrReactingFractionReport
evaluate_pasr_reacting_fraction(double tau_mix_s,
                                double tau_chem_s) noexcept;

struct CombustionCandidate {
  bool available{};
  std::vector<double> species_density_delta_kg_per_m3;
  double integrated_thermochemical_enthalpy_delta_j_per_m3{};
  double integrated_heat_release_j_per_m3{};
};

enum class CombustionClosureStatus : std::uint8_t {
  success,
  invalid_input,
  composition_mismatch,
  representation_mismatch,
  mixing_timescale_unavailable,
  chemical_timescale_unavailable,
  backend_failure,
  backend_contract_failure,
  conservation_failure,
  unsupported_portable_model,
  workspace_failure
};

struct CombustionClosureRequest {
  std::uint64_t composition_fingerprint{};
  std::uint64_t source_state_revision{};
  std::uint64_t accepted_step{};
  double start_time_s{};
  double duration_s{};
  ThermochemicalState mean_state;
  std::vector<ThermochemicalState> stochastic_fields;
};

struct CombustionClosureConfig {
  TciClosureKind tci_closure{TciClosureKind::finite_rate_mean};
  ChemistryRepresentationKind chemistry_representation{
      ChemistryRepresentationKind::direct_cantera};
  ChemicalTimescaleModel chemical_timescale_model{
      ChemicalTimescaleModel::reactant_depletion_l1_v1};
  MixingTimeInput mixing_time{};
};

struct CombustionClosureReport {
  CombustionClosureStatus status{CombustionClosureStatus::invalid_input};
  TciClosureKind configured_tci{TciClosureKind::finite_rate_mean};
  ChemistryRepresentationKind configured_chemistry{
      ChemistryRepresentationKind::direct_cantera};
  std::string_view evaluator_model_id{"finite_rate_mean_shadow"};
  std::string_view chemistry_representation_id{"direct_cantera"};
  std::string_view chemical_timescale_model_id{
      "reactant_depletion_l1_v1"};
  std::uint64_t model_identity{};
  std::uint64_t source_state_revision{};
  std::uint64_t accepted_step{};
  std::uint32_t stochastic_field_count{};
  std::uint32_t chemistry_call_count{};
  std::uint64_t chemistry_internal_step_count{};
  PasrTimescaleReport timescales{};
  CombustionCandidate candidate;
  bool finite_rate_mean_shadow{};
  bool source_published{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == CombustionClosureStatus::success &&
           candidate.available && !source_published;
  }
};

// Vector-backed portable value evaluator, not an allocation-free field loop.
// Closure-owned allocation failure returns workspace_failure and an empty
// unavailable candidate. noexcept providers must report their own failures.
[[nodiscard]] CombustionClosureReport evaluate_combustion_closure(
    const CombustionClosureRequest& request,
    const CombustionClosureConfig& config,
    ChemistryRepresentationAdapter& chemistry,
    ChemicalRateQueryProvider* rate_query) noexcept;

// One physical, interval-integrated source per unit volume. The source is
// added identically to every field; it is neither multiplied nor divided by N.
struct EsfCommonSource {
  std::uint64_t composition_fingerprint{};
  double density_delta_kg_per_m3{};
  std::vector<double> species_density_delta_kg_per_m3;
  double thermochemical_enthalpy_density_delta_j_per_m3{};
};

enum class EsfCommonSourceStatus : std::uint8_t {
  success,
  invalid_input,
  composition_mismatch,
  mean_mismatch,
  source_mass_mismatch,
  inadmissible_state,
  workspace_failure
};

struct EsfCommonSourceReport {
  EsfCommonSourceStatus status{EsfCommonSourceStatus::invalid_input};
  std::string_view model_id{"esf_common_conservative_source_v1"};
  bool available{};
  ThermochemicalState mean_state;
  std::vector<ThermochemicalState> fields;
  [[nodiscard]] bool succeeded() const noexcept {
    return status == EsfCommonSourceStatus::success && available;
  }
};

// Equal-weight N=2/4 value map at a common pressure and density. Validates the
// input mean against its ensemble, conservatively adds rho/rhoY/rhoh, and
// rejects an inadmissible field without clipping or partial publication.
// There is no parcel loop, chemistry call, field storage or MPI side effect.
// EOS, molecular/element identity mapping and the product common-source
// transaction remain responsibilities of the accepted-head adapter.
[[nodiscard]] EsfCommonSourceReport apply_esf_common_source(
    const ThermochemicalState& mean,
    const std::vector<ThermochemicalState>& fields,
    const ChemistryIdentity& identity, const EsfCommonSource& source) noexcept;

} // namespace hundun::v04::combustion
