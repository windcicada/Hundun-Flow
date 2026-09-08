// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include <array>
#include <cstdint>
#include <string_view>

namespace hundun::v04::spray {

using Vector3 = std::array<double, 3U>;

struct ParcelId {
  std::uint64_t high{};
  std::uint64_t low{};

  friend constexpr bool operator==(ParcelId left, ParcelId right) noexcept {
    return left.high == right.high && left.low == right.low;
  }
  friend constexpr bool operator!=(ParcelId left, ParcelId right) noexcept {
    return !(left == right);
  }
  friend constexpr bool operator<(ParcelId left, ParcelId right) noexcept {
    return left.high < right.high ||
           (left.high == right.high && left.low < right.low);
  }
};

enum class ParcelRandomPurpose : std::uint64_t {
  stable_id = 0U,
  injection_direction = 1U,
  breakup_child = 2U
};

struct ParcelRandomAddress {
  std::uint64_t seed{};
  std::uint64_t accepted_step{};
  std::uint64_t stream_identity{};
  std::uint64_t ordinal{};
  ParcelRandomPurpose purpose{ParcelRandomPurpose::stable_id};
};

[[nodiscard]] std::string_view parcel_rng_model_id() noexcept;
[[nodiscard]] std::uint64_t
parcel_random_u64(const ParcelRandomAddress& address,
                  std::uint64_t lane) noexcept;
[[nodiscard]] double parcel_uniform_01(const ParcelRandomAddress& address,
                                       std::uint64_t lane) noexcept;
[[nodiscard]] ParcelId
make_stable_parcel_id(const ParcelRandomAddress& address) noexcept;

struct SprayParcelState {
  ParcelId id{};
  Vector3 position_m{};
  Vector3 velocity_m_per_s{};
  double droplet_mass_kg{};
  double droplet_diameter_m{};
  double multiplicity{};
  double temperature_k{};
  std::uint64_t liquid_material_fingerprint{};
  std::uint64_t owner_global_cell{};
  double age_s{};
};

enum class ParcelStateStatus : std::uint8_t { success, invalid_input };

[[nodiscard]] ParcelStateStatus
validate_parcel_state(const SprayParcelState& parcel) noexcept;

struct LiquidProperties {
  double density_kg_per_m3{};
  double cp_j_per_kg_k{};
  double latent_heat_j_per_kg{};
  double saturation_pressure_pa{};
  double surface_tension_n_per_m{};
  double viscosity_pa_s{};

  friend constexpr bool operator==(const LiquidProperties& left,
                                   const LiquidProperties& right) noexcept {
    return left.density_kg_per_m3 == right.density_kg_per_m3 &&
           left.cp_j_per_kg_k == right.cp_j_per_kg_k &&
           left.latent_heat_j_per_kg == right.latent_heat_j_per_kg &&
           left.saturation_pressure_pa == right.saturation_pressure_pa &&
           left.surface_tension_n_per_m == right.surface_tension_n_per_m &&
           left.viscosity_pa_s == right.viscosity_pa_s;
  }
  friend constexpr bool operator!=(const LiquidProperties& left,
                                   const LiquidProperties& right) noexcept {
    return !(left == right);
  }
};

struct LiquidPropertyQuery {
  std::uint64_t material_fingerprint{};
  double temperature_k{};
};

enum class LiquidPropertyStatus : std::uint8_t {
  success,
  invalid_input,
  unknown_material,
  temperature_out_of_range,
  correlation_domain_error,
  non_finite_output,
  provider_contract_failure
};

struct LiquidPropertyReport {
  LiquidPropertyStatus status{LiquidPropertyStatus::invalid_input};
  std::uint64_t material_fingerprint{};
  double evaluated_temperature_k{};
  LiquidProperties properties{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == LiquidPropertyStatus::success;
  }
};

class LiquidPropertyProvider {
 public:
  virtual ~LiquidPropertyProvider() = default;
  virtual LiquidPropertyReport
  evaluate(const LiquidPropertyQuery& query) const noexcept = 0;
};

[[nodiscard]] LiquidPropertyStatus validate_liquid_property_report(
    const LiquidPropertyQuery& query,
    const LiquidPropertyReport& report) noexcept;

struct ParcelGasExchangeCandidate {
  bool available{};
  double parcel_mass_delta_kg{};
  Vector3 parcel_momentum_delta_kg_m_per_s{};
  double parcel_energy_delta_j{};
  double gas_mass_delta_kg{};
  Vector3 gas_momentum_delta_kg_m_per_s{};
  double gas_energy_delta_j{};
};

enum class ExchangeConservationStatus : std::uint8_t {
  success,
  invalid_input,
  conservation_failure
};

struct ExchangeConservationReport {
  ExchangeConservationStatus status{
      ExchangeConservationStatus::invalid_input};
  double mass_residual_kg{};
  Vector3 momentum_residual_kg_m_per_s{};
  double energy_residual_j{};
};

[[nodiscard]] ParcelGasExchangeCandidate make_conservative_exchange(
    double parcel_mass_delta_kg,
    Vector3 parcel_momentum_delta_kg_m_per_s,
    double parcel_energy_delta_j) noexcept;

[[nodiscard]] ExchangeConservationReport evaluate_exchange_conservation(
    const ParcelGasExchangeCandidate& candidate) noexcept;

struct SchillerNaumannDragInput {
  Vector3 gas_velocity_m_per_s{};
  Vector3 parcel_velocity_m_per_s{};
  double gas_density_kg_per_m3{};
  double gas_dynamic_viscosity_pa_s{};
  double droplet_diameter_m{};
  double droplet_mass_kg{};
  double multiplicity{};
  double duration_s{};
};

enum class DragKernelStatus : std::uint8_t {
  success,
  invalid_input,
  non_finite_output
};

struct SchillerNaumannDragReport {
  DragKernelStatus status{DragKernelStatus::invalid_input};
  std::string_view model_id{"schiller_naumann_v1"};
  double slip_speed_m_per_s{};
  double reynolds_number{};
  double drag_coefficient{};
  bool used_stokes_limit{};
  bool used_constant_drag_branch{};
  Vector3 force_on_one_droplet_n{};
  Vector3 acceleration_m_per_s2{};
  ParcelGasExchangeCandidate exchange{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == DragKernelStatus::success && exchange.available;
  }
};

[[nodiscard]] SchillerNaumannDragReport
evaluate_schiller_naumann_drag(
    const SchillerNaumannDragInput& input) noexcept;

struct RanzMarshallInput {
  double reynolds_number{};
  double prandtl_number{};
  double schmidt_number{};
  double gas_thermal_conductivity_w_per_m_k{};
  double vapor_diffusivity_m2_per_s{};
  double droplet_diameter_m{};
};

enum class HeatMassTransferStatus : std::uint8_t {
  success,
  invalid_input,
  non_finite_output
};

struct RanzMarshallReport {
  HeatMassTransferStatus status{HeatMassTransferStatus::invalid_input};
  std::string_view model_id{"ranz_marshall_v1"};
  double nusselt_number{};
  double sherwood_number{};
  double heat_transfer_coefficient_w_per_m2_k{};
  double mass_transfer_coefficient_m_per_s{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == HeatMassTransferStatus::success;
  }
};

[[nodiscard]] RanzMarshallReport
evaluate_ranz_marshall(const RanzMarshallInput& input) noexcept;

struct SingleComponentEvaporationInput {
  Vector3 parcel_velocity_m_per_s{};
  double droplet_mass_kg{};
  double droplet_diameter_m{};
  double multiplicity{};
  double droplet_temperature_k{};
  double gas_temperature_k{};
  double gas_pressure_pa{};
  double gas_density_kg_per_m3{};
  double gas_vapor_mass_fraction{};
  double vapor_molecular_weight_kg_per_kmol{};
  double carrier_molecular_weight_kg_per_kmol{};
  double reynolds_number{};
  double prandtl_number{};
  double schmidt_number{};
  double gas_thermal_conductivity_w_per_m_k{};
  double vapor_diffusivity_m2_per_s{};
  double duration_s{};
  LiquidProperties liquid_properties{};
};

enum class EvaporationKernelStatus : std::uint8_t {
  success,
  invalid_input,
  inconsistent_droplet_geometry,
  phase_equilibrium_unavailable,
  condensation_unsupported,
  heat_mass_transfer_failure,
  non_finite_output
};

struct SingleComponentEvaporationReport {
  EvaporationKernelStatus status{EvaporationKernelStatus::invalid_input};
  std::string_view model_id{
      "single_component_spalding_ranz_marshall_v1"};
  RanzMarshallReport transfer{};
  double surface_vapor_mole_fraction{};
  double surface_vapor_mass_fraction{};
  double spalding_mass_number{};
  double evaporation_mass_rate_kg_per_s{};
  double convective_heat_rate_w{};
  double latent_heat_rate_w{};
  double droplet_temperature_rate_k_per_s{};
  double diameter_squared_rate_m2_per_s{};
  bool evaporation_active{};
  bool has_terminal_event{};
  bool complete_evaporation{};
  double event_time_s{};
  double advanced_duration_s{};
  double final_droplet_mass_kg{};
  double final_droplet_diameter_m{};
  double final_droplet_temperature_k{};
  ParcelGasExchangeCandidate exchange{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == EvaporationKernelStatus::success &&
           transfer.succeeded() && exchange.available;
  }
};

[[nodiscard]] SingleComponentEvaporationReport
evaluate_single_component_evaporation(
    const SingleComponentEvaporationInput& input) noexcept;

struct TabCoefficients {
  double forcing_coefficient{};
  double breakup_coefficient{};
  double stiffness_coefficient{};
  double damping_coefficient{};
};

struct TabBreakupInput {
  double initial_deformation{};
  double initial_deformation_rate_per_s{};
  double relative_speed_m_per_s{};
  double gas_density_kg_per_m3{};
  double liquid_density_kg_per_m3{};
  double liquid_viscosity_pa_s{};
  double surface_tension_n_per_m{};
  double droplet_radius_m{};
  double breakup_threshold{};
  double duration_s{};
  TabCoefficients coefficients{};
};

enum class TabKernelStatus : std::uint8_t {
  success,
  invalid_input,
  event_search_unavailable,
  non_finite_output
};

struct TabDeformationCandidate {
  bool available{};
  double deformation{};
  double deformation_rate_per_s{};
};

struct TabBreakupReport {
  TabKernelStatus status{TabKernelStatus::invalid_input};
  std::string_view model_id{"tab_linear_oscillator_v1"};
  double forcing_per_s2{};
  double damping_per_s{};
  double stiffness_per_s2{};
  bool threshold_crossed{};
  bool breakup_requested{};
  double event_time_s{};
  double advanced_duration_s{};
  TabDeformationCandidate candidate{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == TabKernelStatus::success && candidate.available;
  }
};

[[nodiscard]] TabBreakupReport
evaluate_tab_breakup(const TabBreakupInput& input) noexcept;

}  // namespace hundun::v04::spray
