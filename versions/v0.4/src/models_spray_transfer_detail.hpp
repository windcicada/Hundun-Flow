// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_spray.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace hundun::v04::spray::detail {

// Polynomial correlations use theta = T - reference_temperature_k and are
// evaluated by Horner's rule.  A constant correlation consumes only c[0].
enum class TemperatureCorrelationKind : std::uint8_t {
  constant,
  polynomial_cubic,
  // Fixed Rachner saturated-liquid coefficients; density only. c[] stays zero.
  kerosene_density_v1,
  // Fixed phase correlations, used only for cp and latent heat respectively.
  // Their four user coefficients stay zero.
  kerosene_cp_v1,
  kerosene_latent_v1
};

struct TemperatureCorrelation {
  TemperatureCorrelationKind kind{TemperatureCorrelationKind::constant};
  double reference_temperature_k{};
  std::array<double, 4U> c{};
};

enum class SaturationPressureCorrelationKind : std::uint8_t {
  antoine_kelvin,
  clausius_clapeyron,
  kerosene_v1
};

struct SaturationPressureCorrelation {
  SaturationPressureCorrelationKind kind{
      SaturationPressureCorrelationKind::antoine_kelvin};

  // Antoine: p_sat = pressure_scale_pa * 10^(a-b/(T+c)).
  double antoine_a{};
  double antoine_b_k{};
  double antoine_c_k{};
  double pressure_scale_pa{};

  // Clausius--Clapeyron:
  // p_sat = p_ref exp[L M/R (1/T_ref - 1/T)], R=8314.462618... J/kmol/K.
  double reference_pressure_pa{};
  double reference_temperature_k{};
  double latent_heat_j_per_kg{};
  double molecular_weight_kg_per_kmol{};
};

struct LiquidPropertyPack {
  std::uint64_t material_fingerprint{};
  double minimum_temperature_k{};
  double maximum_temperature_k{};
  TemperatureCorrelation density_kg_per_m3{};
  TemperatureCorrelation cp_j_per_kg_k{};
  TemperatureCorrelation latent_heat_j_per_kg{};
  TemperatureCorrelation surface_tension_n_per_m{};
  TemperatureCorrelation viscosity_pa_s{};
  SaturationPressureCorrelation saturation_pressure{};
};

class LiquidPropertyService final : public LiquidPropertyProvider {
 public:
  constexpr LiquidPropertyService(const LiquidPropertyPack* packs,
                                  std::size_t count) noexcept
      : packs_(packs), count_(count) {}

  [[nodiscard]] LiquidPropertyReport
  evaluate(const LiquidPropertyQuery& query) const noexcept override;

 private:
  const LiquidPropertyPack* packs_{};
  std::size_t count_{};
};

enum class FilmSampleStatus : std::uint8_t {
  success,
  invalid_input,
  non_finite_output
};

// The transport values are supplied by the future Stage-4 service adapter at
// the reference (film) state.  This pure slice owns only the one-third state
// algebra and validates that adapter payload.
struct OneThirdFilmInput {
  Vector3 far_gas_velocity_m_per_s{};
  double droplet_surface_temperature_k{};
  double far_gas_temperature_k{};
  double surface_vapor_mass_fraction{};
  double far_gas_vapor_mass_fraction{};
  double thermodynamic_pressure_pa{};
  double reference_density_kg_per_m3{};
  double reference_dynamic_viscosity_pa_s{};
  double reference_thermal_conductivity_w_per_m_k{};
  double reference_vapor_diffusivity_m2_per_s{};
  double reference_cp_j_per_kg_k{};
};

struct OneThirdFilmSample {
  FilmSampleStatus status{FilmSampleStatus::invalid_input};
  std::string_view model_id{"one_third_film_v1"};
  double surface_weight{};
  double far_gas_weight{};
  Vector3 gas_velocity_m_per_s{};
  double surface_temperature_k{};
  double far_gas_temperature_k{};
  double reference_temperature_k{};
  double surface_vapor_mass_fraction{};
  double far_gas_vapor_mass_fraction{};
  double reference_vapor_mass_fraction{};
  double thermodynamic_pressure_pa{};
  double density_kg_per_m3{};
  double dynamic_viscosity_pa_s{};
  double thermal_conductivity_w_per_m_k{};
  double vapor_diffusivity_m2_per_s{};
  double cp_j_per_kg_k{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == FilmSampleStatus::success;
  }
};

[[nodiscard]] OneThirdFilmSample
sample_one_third_film(const OneThirdFilmInput& input) noexcept;

struct TransferExchangeCandidate {
  bool available{};

  // Actual parcel changes, new minus old, after multiplicity is applied once.
  double parcel_liquid_mass_delta_kg{};
  Vector3 parcel_momentum_delta_kg_m_per_s{};
  double parcel_thermochemical_enthalpy_delta_j{};
  double parcel_kinetic_energy_delta_j{};

  // Resolved physical categories.  These remain separate because the product
  // enthalpy and pressure-work admission algebra is owned by the driver seam.
  double vapor_mass_to_gas_kg{};
  Vector3 vapor_momentum_to_gas_kg_m_per_s{};
  double vapor_absolute_thermochemical_enthalpy_to_gas_j{};
  Vector3 drag_momentum_to_parcel_kg_m_per_s{};
  double convective_heat_to_parcel_j{};
  // Integral of F_drag,on-parcel dot u_parcel.  It is a mechanical-work
  // category, not part of thermal_exchange_to_gas_j.
  double drag_work_to_parcel_j{};

  // Resolved mass/momentum transfers.  No gas endpoint or field is touched by
  // this module; the later deposition kernel must combine these with its
  // actual gas state and pressure-work convention.
  double gas_mass_delta_kg{};
  Vector3 gas_momentum_delta_kg_m_per_s{};
  double mass_closure_residual_kg{};
  Vector3 momentum_closure_residual_kg_m_per_s{};

  // Phase/thermal exchange only: vapor absolute thermochemical enthalpy minus
  // heat delivered to the parcel.  This is not a complete product rho*h-p
  // source and contains neither kinetic energy nor a pressure-work repair.
  double thermal_exchange_to_gas_j{};
  // parcel_delta_H_tc + thermal_exchange_to_gas.  A nonzero value is reported
  // to the later admission/refinement layer. A-S reports its independently
  // integrated heat flux; THICK_EX defines its source from the same endpoint
  // increment as the liquid update and its explicit temperature limit.
  double thermal_exchange_state_residual_j{};

  // Momentum needed to force an exact endpoint-state balance beyond the
  // resolved vapor/drag quadrature.  It is observable and never silently
  // deposited.
  Vector3 momentum_quadrature_residual_kg_m_per_s{};
};

enum class AbramzonSirignanoStatus : std::uint8_t {
  success,
  invalid_input,
  inconsistent_droplet_geometry,
  phase_equilibrium_unavailable,
  condensation_unsupported,
  iteration_not_converged,
  non_finite_output
};

struct AbramzonSirignanoInput {
  SprayParcelState parcel{};
  // Must be sampled for this exact parcel surface temperature. Properties at
  // an older predictor temperature cannot be reused as a valid film sample.
  OneThirdFilmSample film{};
  LiquidProperties liquid_properties{};
  double vapor_absolute_thermochemical_enthalpy_j_per_kg{};
  double duration_s{};
  std::uint32_t maximum_bt_iterations{64U};
  double bt_relative_tolerance{1.0e-10};
};

struct AbramzonSirignanoReport {
  AbramzonSirignanoStatus status{AbramzonSirignanoStatus::invalid_input};
  std::string_view model_id{"abramzon_sirignano_single_component_v1"};
  double reynolds_number{};
  double prandtl_number{};
  double schmidt_number{};
  double spalding_mass_number{};
  double spalding_heat_number{};
  double base_nusselt_number{};
  double base_sherwood_number{};
  double modified_nusselt_number{};
  double modified_sherwood_number{};
  std::uint32_t bt_iterations{};
  bool bt_converged{};
  bool saturated{};
  bool heating_only{};
  double evaporation_mass_rate_one_droplet_kg_per_s{};
  double convective_heat_rate_one_droplet_w{};
  double diameter_squared_rate_m2_per_s{};
  double temperature_rate_k_per_s{};
  bool has_terminal_event{};
  bool complete_evaporation{};
  double event_time_s{};
  double advanced_duration_s{};
  SprayParcelState candidate_parcel{};
  TransferExchangeCandidate exchange{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == AbramzonSirignanoStatus::success && exchange.available;
  }
};

[[nodiscard]] AbramzonSirignanoReport evaluate_abramzon_sirignano(
    const AbramzonSirignanoInput& input) noexcept;

enum class ThickExchangeStatus : std::uint8_t {
  success,
  invalid_input,
  inconsistent_droplet_geometry,
  non_finite_output
};

// Frozen-property THICK_EX evaporation interval. The gas-property adapter
// supplies the local-pressure boiling temperature separately from drop T.
// Drag, trajectory and subsequent density/diameter refresh belong to the
// parcel interval driver; all exchange entries here include multiplicity once.
struct ThickExchangeInput {
  SprayParcelState parcel{};
  LiquidProperties liquid_properties{};
  double gas_temperature_k{};
  double gas_cp_j_per_kg_k{};
  double gas_dynamic_viscosity_pa_s{};
  double prandtl_number{};
  double boiling_temperature_k{};
  double vapor_absolute_thermochemical_enthalpy_j_per_kg{};
  double duration_s{};
};

struct ThickExchangeReport {
  ThickExchangeStatus status{ThickExchangeStatus::invalid_input};
  std::string_view model_id{"thick_exchange_layer_v1"};
  double spalding_heat_number{};
  double layer_conductivity_w_per_m_k{};
  double nusselt_number{};
  double diameter_squared_loss_rate_m2_per_s{};
  double convective_heat_transfer_w_per_k{};
  // Requested-interval averages, matching the reference dmdt/dhdt outputs.
  double mean_evaporation_rate_one_droplet_kg_per_s{};
  double mean_liquid_enthalpy_rate_one_droplet_w{};
  bool temperature_limited{};
  double temperature_limit_energy_one_droplet_j{};
  bool complete_evaporation{};
  double event_time_s{};
  double advanced_duration_s{};
  SprayParcelState candidate_parcel{};
  TransferExchangeCandidate exchange{};
  [[nodiscard]] bool succeeded() const noexcept {
    return status == ThickExchangeStatus::success && exchange.available;
  }
};

[[nodiscard]] ThickExchangeReport evaluate_thick_exchange(
    const ThickExchangeInput& input) noexcept;

struct GasTransferEnvironment {
  Vector3 gas_velocity_m_per_s{};
  double gas_temperature_k{};
  double gas_vapor_mass_fraction{};
  double thermodynamic_pressure_pa{};
  double film_density_kg_per_m3{};
  double film_dynamic_viscosity_pa_s{};
  double film_thermal_conductivity_w_per_m_k{};
  double film_vapor_diffusivity_m2_per_s{};
  double film_cp_j_per_kg_k{};
  double vapor_molecular_weight_kg_per_kmol{};
  double carrier_molecular_weight_kg_per_kmol{};
  double vapor_absolute_thermochemical_enthalpy_j_per_kg{};
};

enum class FixedExchangeStatus : std::uint8_t {
  success,
  invalid_input,
  liquid_property_failure,
  film_sample_failure,
  transfer_failure,
  substep_limit_exceeded,
  non_finite_output
};

struct FixedExchangeInput {
  SprayParcelState committed_parcel{};
  const LiquidPropertyService* liquid_properties{};
  GasTransferEnvironment predictor_environment{};
  GasTransferEnvironment corrector_environment{};
  double duration_s{};
  std::uint32_t internal_substeps{1U};
  std::uint32_t maximum_internal_substeps{1U};
  std::uint32_t maximum_bt_iterations{64U};
  double bt_relative_tolerance{1.0e-10};
};

struct FixedExchangeReport {
  FixedExchangeStatus status{FixedExchangeStatus::invalid_input};
  std::string_view model_id{"fixed_predictor_corrector_exchange_v1"};
  std::uint32_t predictor_passes{};
  std::uint32_t corrector_passes{};
  std::uint32_t predictor_evaluations{};
  std::uint32_t corrector_evaluations{};
  std::uint32_t completed_substeps{};
  bool has_terminal_event{};
  double event_time_s{};
  double advanced_duration_s{};
  SprayParcelState candidate_parcel{};
  TransferExchangeCandidate exchange{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == FixedExchangeStatus::success && exchange.available;
  }
};

[[nodiscard]] FixedExchangeReport integrate_fixed_exchange(
    const FixedExchangeInput& input) noexcept;

}  // namespace hundun::v04::spray::detail
