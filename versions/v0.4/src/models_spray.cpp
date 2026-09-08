// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_spray.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::spray {
namespace {

constexpr double kAbsoluteTolerance = 1.0e-12;
constexpr double kRelativeTolerance = 1.0e-10;

bool finite_vector(const Vector3& value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
}

bool within_tolerance(double residual, double left, double right) noexcept {
  return std::abs(residual) <=
         kAbsoluteTolerance +
             kRelativeTolerance * (std::abs(left) + std::abs(right));
}

struct TabLinearState {
  double deformation{};
  double rate_per_s{};
};

TabLinearState evaluate_tab_linear_state(double initial_deformation,
                                         double initial_rate_per_s,
                                         double forcing_per_s2,
                                         double damping_per_s,
                                         double stiffness_per_s2,
                                         double time_s) noexcept {
  if (time_s == 0.0) {
    return {initial_deformation, initial_rate_per_s};
  }
  if (stiffness_per_s2 == 0.0) {
    if (damping_per_s == 0.0) {
      return {initial_deformation + initial_rate_per_s * time_s +
                  0.5 * forcing_per_s2 * time_s * time_s,
              initial_rate_per_s + forcing_per_s2 * time_s};
    }
    const double terminal_rate = forcing_per_s2 / damping_per_s;
    const double transient_rate = initial_rate_per_s - terminal_rate;
    const double decay = std::exp(-damping_per_s * time_s);
    return {initial_deformation + terminal_rate * time_s +
                transient_rate * (1.0 - decay) / damping_per_s,
            terminal_rate + transient_rate * decay};
  }

  const double equilibrium = forcing_per_s2 / stiffness_per_s2;
  const double initial_offset = initial_deformation - equilibrium;
  const double half_damping = 0.5 * damping_per_s;
  const double discriminant =
      stiffness_per_s2 - half_damping * half_damping;
  const double regime_scale =
      std::max({1.0, stiffness_per_s2,
                half_damping * half_damping});
  const double regime_tolerance =
      32.0 * std::numeric_limits<double>::epsilon() * regime_scale;
  const double decay = std::exp(-half_damping * time_s);
  if (discriminant > regime_tolerance) {
    const double frequency = std::sqrt(discriminant);
    const double sine_coefficient =
        (initial_rate_per_s + half_damping * initial_offset) /
        frequency;
    const double cosine = std::cos(frequency * time_s);
    const double sine = std::sin(frequency * time_s);
    const double oscillation =
        initial_offset * cosine + sine_coefficient * sine;
    const double oscillation_rate =
        -initial_offset * frequency * sine +
        sine_coefficient * frequency * cosine;
    return {equilibrium + decay * oscillation,
            decay * (oscillation_rate - half_damping * oscillation)};
  }
  if (discriminant < -regime_tolerance) {
    const double growth = std::sqrt(-discriminant);
    const double sinh_coefficient =
        (initial_rate_per_s + half_damping * initial_offset) /
        growth;
    const double hyperbolic_cosine = std::cosh(growth * time_s);
    const double hyperbolic_sine = std::sinh(growth * time_s);
    const double transient =
        initial_offset * hyperbolic_cosine +
        sinh_coefficient * hyperbolic_sine;
    const double transient_rate =
        initial_offset * growth * hyperbolic_sine +
        sinh_coefficient * growth * hyperbolic_cosine;
    return {equilibrium + decay * transient,
            decay * (transient_rate - half_damping * transient)};
  }

  const double linear_coefficient =
      initial_rate_per_s + half_damping * initial_offset;
  const double transient =
      initial_offset + linear_coefficient * time_s;
  return {equilibrium + decay * transient,
          decay * (linear_coefficient - half_damping * transient)};
}

bool finite_tab_state(const TabLinearState& state) noexcept {
  return std::isfinite(state.deformation) &&
         std::isfinite(state.rate_per_s);
}

bool opposite_sign(double left, double right) noexcept {
  return (left < 0.0 && right > 0.0) ||
         (left > 0.0 && right < 0.0);
}

std::uint64_t splitmix64(std::uint64_t value) noexcept {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

}  // namespace

std::string_view parcel_rng_model_id() noexcept {
  return "parcel_counter_splitmix64_v1";
}

std::uint64_t parcel_random_u64(const ParcelRandomAddress& address,
                                std::uint64_t lane) noexcept {
  std::uint64_t counter = address.seed;
  counter += address.accepted_step * UINT64_C(0xd2b74407b1ce6e93);
  counter += address.stream_identity * UINT64_C(0xca5a826395121157);
  counter += address.ordinal * UINT64_C(0x9e3779b185ebca87);
  counter += static_cast<std::uint64_t>(address.purpose) *
             UINT64_C(0x8cb92ba72f3d8dd7);
  counter += lane * UINT64_C(0xdb4f0b9175ae2165);
  return splitmix64(counter);
}

double parcel_uniform_01(const ParcelRandomAddress& address,
                         std::uint64_t lane) noexcept {
  constexpr double inverse_two_to_53 = 1.0 / 9007199254740992.0;
  return static_cast<double>(parcel_random_u64(address, lane) >> 11U) *
         inverse_two_to_53;
}

ParcelId make_stable_parcel_id(const ParcelRandomAddress& address) noexcept {
  ParcelRandomAddress id_address = address;
  id_address.purpose = ParcelRandomPurpose::stable_id;
  ParcelId result{parcel_random_u64(id_address, 0U),
                  parcel_random_u64(id_address, 1U)};
  if (result.high == 0U && result.low == 0U) {
    result.low = 1U;
  }
  return result;
}

ParcelStateStatus
validate_parcel_state(const SprayParcelState& parcel) noexcept {
  if ((parcel.id.high == 0U && parcel.id.low == 0U) ||
      !finite_vector(parcel.position_m) ||
      !finite_vector(parcel.velocity_m_per_s) ||
      !std::isfinite(parcel.droplet_mass_kg) ||
      !(parcel.droplet_mass_kg > 0.0) ||
      !std::isfinite(parcel.droplet_diameter_m) ||
      !(parcel.droplet_diameter_m > 0.0) ||
      !std::isfinite(parcel.multiplicity) || !(parcel.multiplicity > 0.0) ||
      !std::isfinite(parcel.temperature_k) ||
      !(parcel.temperature_k > 0.0) ||
      parcel.liquid_material_fingerprint == 0U ||
      !std::isfinite(parcel.age_s) || parcel.age_s < 0.0) {
    return ParcelStateStatus::invalid_input;
  }
  return ParcelStateStatus::success;
}

LiquidPropertyStatus validate_liquid_property_report(
    const LiquidPropertyQuery& query,
    const LiquidPropertyReport& report) noexcept {
  if (query.material_fingerprint == 0U ||
      !std::isfinite(query.temperature_k) || !(query.temperature_k > 0.0)) {
    return LiquidPropertyStatus::invalid_input;
  }
  const LiquidProperties& properties = report.properties;
  if (!report.succeeded()) {
    const double values[]{properties.density_kg_per_m3,
                          properties.cp_j_per_kg_k,
                          properties.latent_heat_j_per_kg,
                          properties.saturation_pressure_pa,
                          properties.surface_tension_n_per_m,
                          properties.viscosity_pa_s,
                          report.evaluated_temperature_k};
    bool canonical = report.material_fingerprint == 0U;
    for (const double value : values) {
      canonical &= value == 0.0 && !std::signbit(value);
    }
    return canonical ? report.status
                     : LiquidPropertyStatus::provider_contract_failure;
  }
  if (report.material_fingerprint != query.material_fingerprint ||
      report.evaluated_temperature_k != query.temperature_k ||
      !std::isfinite(properties.density_kg_per_m3) ||
      !(properties.density_kg_per_m3 > 0.0) ||
      !std::isfinite(properties.cp_j_per_kg_k) ||
      !(properties.cp_j_per_kg_k > 0.0) ||
      !std::isfinite(properties.latent_heat_j_per_kg) ||
      properties.latent_heat_j_per_kg < 0.0 ||
      !std::isfinite(properties.saturation_pressure_pa) ||
      properties.saturation_pressure_pa < 0.0 ||
      !std::isfinite(properties.surface_tension_n_per_m) ||
      properties.surface_tension_n_per_m < 0.0 ||
      !std::isfinite(properties.viscosity_pa_s) ||
      properties.viscosity_pa_s < 0.0) {
    return LiquidPropertyStatus::provider_contract_failure;
  }
  return LiquidPropertyStatus::success;
}

ParcelGasExchangeCandidate make_conservative_exchange(
    double parcel_mass_delta_kg,
    Vector3 parcel_momentum_delta_kg_m_per_s,
    double parcel_energy_delta_j) noexcept {
  ParcelGasExchangeCandidate result;
  if (!std::isfinite(parcel_mass_delta_kg) ||
      !finite_vector(parcel_momentum_delta_kg_m_per_s) ||
      !std::isfinite(parcel_energy_delta_j)) {
    return result;
  }
  result.available = true;
  result.parcel_mass_delta_kg = parcel_mass_delta_kg;
  result.parcel_momentum_delta_kg_m_per_s =
      parcel_momentum_delta_kg_m_per_s;
  result.parcel_energy_delta_j = parcel_energy_delta_j;
  result.gas_mass_delta_kg = -parcel_mass_delta_kg;
  for (std::size_t component = 0U; component < 3U; ++component) {
    result.gas_momentum_delta_kg_m_per_s[component] =
        -parcel_momentum_delta_kg_m_per_s[component];
  }
  result.gas_energy_delta_j = -parcel_energy_delta_j;
  return result;
}

ExchangeConservationReport evaluate_exchange_conservation(
    const ParcelGasExchangeCandidate& candidate) noexcept {
  ExchangeConservationReport report;
  if (!candidate.available ||
      !std::isfinite(candidate.parcel_mass_delta_kg) ||
      !finite_vector(candidate.parcel_momentum_delta_kg_m_per_s) ||
      !std::isfinite(candidate.parcel_energy_delta_j) ||
      !std::isfinite(candidate.gas_mass_delta_kg) ||
      !finite_vector(candidate.gas_momentum_delta_kg_m_per_s) ||
      !std::isfinite(candidate.gas_energy_delta_j)) {
    return report;
  }
  report.mass_residual_kg =
      candidate.parcel_mass_delta_kg + candidate.gas_mass_delta_kg;
  report.energy_residual_j =
      candidate.parcel_energy_delta_j + candidate.gas_energy_delta_j;
  bool conservative = within_tolerance(
      report.mass_residual_kg, candidate.parcel_mass_delta_kg,
      candidate.gas_mass_delta_kg);
  conservative &= within_tolerance(
      report.energy_residual_j, candidate.parcel_energy_delta_j,
      candidate.gas_energy_delta_j);
  for (std::size_t component = 0U; component < 3U; ++component) {
    report.momentum_residual_kg_m_per_s[component] =
        candidate.parcel_momentum_delta_kg_m_per_s[component] +
        candidate.gas_momentum_delta_kg_m_per_s[component];
    conservative &= within_tolerance(
        report.momentum_residual_kg_m_per_s[component],
        candidate.parcel_momentum_delta_kg_m_per_s[component],
        candidate.gas_momentum_delta_kg_m_per_s[component]);
  }
  report.status = conservative ? ExchangeConservationStatus::success
                               : ExchangeConservationStatus::conservation_failure;
  return report;
}

SchillerNaumannDragReport evaluate_schiller_naumann_drag(
    const SchillerNaumannDragInput& input) noexcept {
  SchillerNaumannDragReport report;
  if (!finite_vector(input.gas_velocity_m_per_s) ||
      !finite_vector(input.parcel_velocity_m_per_s) ||
      !std::isfinite(input.gas_density_kg_per_m3) ||
      !(input.gas_density_kg_per_m3 > 0.0) ||
      !std::isfinite(input.gas_dynamic_viscosity_pa_s) ||
      !(input.gas_dynamic_viscosity_pa_s > 0.0) ||
      !std::isfinite(input.droplet_diameter_m) ||
      !(input.droplet_diameter_m > 0.0) ||
      !std::isfinite(input.droplet_mass_kg) ||
      !(input.droplet_mass_kg > 0.0) ||
      !std::isfinite(input.multiplicity) || !(input.multiplicity > 0.0) ||
      !std::isfinite(input.duration_s) || input.duration_s < 0.0) {
    return report;
  }

  Vector3 slip{};
  for (std::size_t component = 0U; component < 3U; ++component) {
    slip[component] = input.gas_velocity_m_per_s[component] -
                      input.parcel_velocity_m_per_s[component];
  }
  report.slip_speed_m_per_s =
      std::hypot(slip[0U], slip[1U], slip[2U]);
  if (!std::isfinite(report.slip_speed_m_per_s)) {
    report = {};
    report.status = DragKernelStatus::non_finite_output;
    return report;
  }
  if (report.slip_speed_m_per_s == 0.0) {
    report.used_stokes_limit = true;
    report.exchange = make_conservative_exchange(0.0, {}, 0.0);
    report.status = DragKernelStatus::success;
    return report;
  }

  report.reynolds_number =
      input.gas_density_kg_per_m3 * report.slip_speed_m_per_s *
      input.droplet_diameter_m / input.gas_dynamic_viscosity_pa_s;
  if (!std::isfinite(report.reynolds_number) ||
      !(report.reynolds_number > 0.0)) {
    report = {};
    report.status = DragKernelStatus::non_finite_output;
    return report;
  }

  constexpr double stokes_reynolds_limit = 1.0e-8;
  constexpr double pi = 3.141592653589793238462643383279502884;
  if (report.reynolds_number <= stokes_reynolds_limit) {
    report.used_stokes_limit = true;
    report.drag_coefficient = 24.0 / report.reynolds_number;
    const double stokes_scale = 3.0 * pi *
                                input.gas_dynamic_viscosity_pa_s *
                                input.droplet_diameter_m;
    for (std::size_t component = 0U; component < 3U; ++component) {
      report.force_on_one_droplet_n[component] =
          stokes_scale * slip[component];
    }
  } else {
    if (report.reynolds_number <= 1000.0) {
      report.drag_coefficient =
          24.0 / report.reynolds_number *
          (1.0 + 0.15 * std::pow(report.reynolds_number, 0.687));
    } else {
      report.drag_coefficient = 0.44;
      report.used_constant_drag_branch = true;
    }
    const double projected_area_m2 =
        0.25 * pi * input.droplet_diameter_m *
        input.droplet_diameter_m;
    const double force_scale =
        0.5 * report.drag_coefficient * input.gas_density_kg_per_m3 *
        projected_area_m2 * report.slip_speed_m_per_s;
    for (std::size_t component = 0U; component < 3U; ++component) {
      report.force_on_one_droplet_n[component] =
          force_scale * slip[component];
    }
  }

  Vector3 parcel_momentum_delta{};
  double parcel_energy_delta = 0.0;
  for (std::size_t component = 0U; component < 3U; ++component) {
    report.acceleration_m_per_s2[component] =
        report.force_on_one_droplet_n[component] /
        input.droplet_mass_kg;
    parcel_momentum_delta[component] =
        report.force_on_one_droplet_n[component] * input.duration_s *
        input.multiplicity;
    parcel_energy_delta +=
        report.force_on_one_droplet_n[component] *
        input.parcel_velocity_m_per_s[component] * input.duration_s *
        input.multiplicity;
  }
  if (!std::isfinite(report.drag_coefficient) ||
      !finite_vector(report.force_on_one_droplet_n) ||
      !finite_vector(report.acceleration_m_per_s2) ||
      !finite_vector(parcel_momentum_delta) ||
      !std::isfinite(parcel_energy_delta)) {
    report = {};
    report.status = DragKernelStatus::non_finite_output;
    return report;
  }
  report.exchange =
      make_conservative_exchange(0.0, parcel_momentum_delta,
                                 parcel_energy_delta);
  if (!report.exchange.available) {
    report = {};
    report.status = DragKernelStatus::non_finite_output;
    return report;
  }
  report.status = DragKernelStatus::success;
  return report;
}

RanzMarshallReport
evaluate_ranz_marshall(const RanzMarshallInput& input) noexcept {
  RanzMarshallReport report;
  if (!std::isfinite(input.reynolds_number) ||
      input.reynolds_number < 0.0 ||
      !std::isfinite(input.prandtl_number) ||
      !(input.prandtl_number > 0.0) ||
      !std::isfinite(input.schmidt_number) ||
      !(input.schmidt_number > 0.0) ||
      !std::isfinite(input.gas_thermal_conductivity_w_per_m_k) ||
      !(input.gas_thermal_conductivity_w_per_m_k > 0.0) ||
      !std::isfinite(input.vapor_diffusivity_m2_per_s) ||
      !(input.vapor_diffusivity_m2_per_s > 0.0) ||
      !std::isfinite(input.droplet_diameter_m) ||
      !(input.droplet_diameter_m > 0.0)) {
    return report;
  }
  const double root_reynolds = std::sqrt(input.reynolds_number);
  report.nusselt_number =
      2.0 + 0.6 * root_reynolds * std::cbrt(input.prandtl_number);
  report.sherwood_number =
      2.0 + 0.6 * root_reynolds * std::cbrt(input.schmidt_number);
  report.heat_transfer_coefficient_w_per_m2_k =
      report.nusselt_number *
      input.gas_thermal_conductivity_w_per_m_k /
      input.droplet_diameter_m;
  report.mass_transfer_coefficient_m_per_s =
      report.sherwood_number * input.vapor_diffusivity_m2_per_s /
      input.droplet_diameter_m;
  if (!std::isfinite(report.nusselt_number) ||
      !std::isfinite(report.sherwood_number) ||
      !std::isfinite(report.heat_transfer_coefficient_w_per_m2_k) ||
      !std::isfinite(report.mass_transfer_coefficient_m_per_s)) {
    report = {};
    report.status = HeatMassTransferStatus::non_finite_output;
    return report;
  }
  report.status = HeatMassTransferStatus::success;
  return report;
}

SingleComponentEvaporationReport evaluate_single_component_evaporation(
    const SingleComponentEvaporationInput& input) noexcept {
  SingleComponentEvaporationReport report;
  const LiquidProperties& liquid = input.liquid_properties;
  if (!finite_vector(input.parcel_velocity_m_per_s) ||
      !std::isfinite(input.droplet_mass_kg) ||
      !(input.droplet_mass_kg > 0.0) ||
      !std::isfinite(input.droplet_diameter_m) ||
      !(input.droplet_diameter_m > 0.0) ||
      !std::isfinite(input.multiplicity) || !(input.multiplicity > 0.0) ||
      !std::isfinite(input.droplet_temperature_k) ||
      !(input.droplet_temperature_k > 0.0) ||
      !std::isfinite(input.gas_temperature_k) ||
      !(input.gas_temperature_k > 0.0) ||
      !std::isfinite(input.gas_pressure_pa) ||
      !(input.gas_pressure_pa > 0.0) ||
      !std::isfinite(input.gas_density_kg_per_m3) ||
      !(input.gas_density_kg_per_m3 > 0.0) ||
      !std::isfinite(input.gas_vapor_mass_fraction) ||
      input.gas_vapor_mass_fraction < 0.0 ||
      !(input.gas_vapor_mass_fraction < 1.0) ||
      !std::isfinite(input.vapor_molecular_weight_kg_per_kmol) ||
      !(input.vapor_molecular_weight_kg_per_kmol > 0.0) ||
      !std::isfinite(input.carrier_molecular_weight_kg_per_kmol) ||
      !(input.carrier_molecular_weight_kg_per_kmol > 0.0) ||
      !std::isfinite(input.duration_s) || input.duration_s < 0.0 ||
      !std::isfinite(liquid.density_kg_per_m3) ||
      !(liquid.density_kg_per_m3 > 0.0) ||
      !std::isfinite(liquid.cp_j_per_kg_k) ||
      !(liquid.cp_j_per_kg_k > 0.0) ||
      !std::isfinite(liquid.latent_heat_j_per_kg) ||
      liquid.latent_heat_j_per_kg < 0.0 ||
      !std::isfinite(liquid.saturation_pressure_pa) ||
      liquid.saturation_pressure_pa < 0.0 ||
      !std::isfinite(liquid.surface_tension_n_per_m) ||
      liquid.surface_tension_n_per_m < 0.0 ||
      !std::isfinite(liquid.viscosity_pa_s) ||
      liquid.viscosity_pa_s < 0.0) {
    return report;
  }

  constexpr double pi = 3.141592653589793238462643383279502884;
  const double diameter_squared =
      input.droplet_diameter_m * input.droplet_diameter_m;
  const double geometric_mass =
      liquid.density_kg_per_m3 * pi / 6.0 * diameter_squared *
      input.droplet_diameter_m;
  const double geometry_scale =
      std::max(std::abs(input.droplet_mass_kg), std::abs(geometric_mass));
  if (!std::isfinite(geometric_mass) ||
      std::abs(input.droplet_mass_kg - geometric_mass) >
          1.0e-10 * geometry_scale) {
    report.status =
        EvaporationKernelStatus::inconsistent_droplet_geometry;
    return report;
  }
  if (!(liquid.saturation_pressure_pa < input.gas_pressure_pa)) {
    report.status =
        EvaporationKernelStatus::phase_equilibrium_unavailable;
    return report;
  }

  report.transfer = evaluate_ranz_marshall(
      {input.reynolds_number, input.prandtl_number,
       input.schmidt_number,
       input.gas_thermal_conductivity_w_per_m_k,
       input.vapor_diffusivity_m2_per_s,
       input.droplet_diameter_m});
  if (!report.transfer.succeeded()) {
    report = {};
    report.status =
        EvaporationKernelStatus::heat_mass_transfer_failure;
    return report;
  }

  report.surface_vapor_mole_fraction =
      liquid.saturation_pressure_pa / input.gas_pressure_pa;
  const double surface_mass_denominator =
      report.surface_vapor_mole_fraction *
          input.vapor_molecular_weight_kg_per_kmol +
      (1.0 - report.surface_vapor_mole_fraction) *
          input.carrier_molecular_weight_kg_per_kmol;
  report.surface_vapor_mass_fraction =
      report.surface_vapor_mole_fraction *
      input.vapor_molecular_weight_kg_per_kmol /
      surface_mass_denominator;
  if (!std::isfinite(report.surface_vapor_mass_fraction) ||
      report.surface_vapor_mass_fraction < 0.0 ||
      !(report.surface_vapor_mass_fraction < 1.0)) {
    report = {};
    report.status =
        EvaporationKernelStatus::phase_equilibrium_unavailable;
    return report;
  }
  if (input.gas_vapor_mass_fraction >
      report.surface_vapor_mass_fraction) {
    report = {};
    report.status = EvaporationKernelStatus::condensation_unsupported;
    return report;
  }

  report.spalding_mass_number =
      (report.surface_vapor_mass_fraction -
       input.gas_vapor_mass_fraction) /
      (1.0 - report.surface_vapor_mass_fraction);
  const double log_mass_transfer =
      std::log1p(report.spalding_mass_number);
  report.evaporation_active = log_mass_transfer > 0.0;
  report.evaporation_mass_rate_kg_per_s =
      pi * input.droplet_diameter_m * input.gas_density_kg_per_m3 *
      input.vapor_diffusivity_m2_per_s *
      report.transfer.sherwood_number * log_mass_transfer;
  report.convective_heat_rate_w =
      pi * input.droplet_diameter_m *
      input.gas_thermal_conductivity_w_per_m_k *
      report.transfer.nusselt_number *
      (input.gas_temperature_k - input.droplet_temperature_k);
  report.latent_heat_rate_w =
      report.evaporation_mass_rate_kg_per_s *
      liquid.latent_heat_j_per_kg;
  report.droplet_temperature_rate_k_per_s =
      (report.convective_heat_rate_w - report.latent_heat_rate_w) /
      (input.droplet_mass_kg * liquid.cp_j_per_kg_k);
  if (report.evaporation_active) {
    const double diameter_squared_loss_rate =
        4.0 * input.gas_density_kg_per_m3 *
        input.vapor_diffusivity_m2_per_s *
        report.transfer.sherwood_number * log_mass_transfer /
        liquid.density_kg_per_m3;
    report.diameter_squared_rate_m2_per_s =
        -diameter_squared_loss_rate;
    report.event_time_s =
        diameter_squared / diameter_squared_loss_rate;
    report.has_terminal_event = true;
  }
  if (!std::isfinite(report.spalding_mass_number) ||
      report.spalding_mass_number < 0.0 ||
      !std::isfinite(report.evaporation_mass_rate_kg_per_s) ||
      report.evaporation_mass_rate_kg_per_s < 0.0 ||
      !std::isfinite(report.convective_heat_rate_w) ||
      !std::isfinite(report.latent_heat_rate_w) ||
      !std::isfinite(report.droplet_temperature_rate_k_per_s) ||
      !std::isfinite(report.diameter_squared_rate_m2_per_s) ||
      (report.has_terminal_event &&
       (!std::isfinite(report.event_time_s) ||
        !(report.event_time_s > 0.0)))) {
    report = {};
    report.status = EvaporationKernelStatus::non_finite_output;
    return report;
  }

  report.advanced_duration_s = input.duration_s;
  if (report.has_terminal_event &&
      input.duration_s >= report.event_time_s) {
    report.advanced_duration_s = report.event_time_s;
    report.complete_evaporation = true;
  }
  if (input.duration_s == 0.0) {
    report.final_droplet_mass_kg = input.droplet_mass_kg;
    report.final_droplet_diameter_m = input.droplet_diameter_m;
    report.final_droplet_temperature_k = input.droplet_temperature_k;
  } else {
    report.final_droplet_temperature_k =
        input.droplet_temperature_k +
        report.droplet_temperature_rate_k_per_s *
            report.advanced_duration_s;
    if (report.complete_evaporation) {
      report.final_droplet_mass_kg = 0.0;
      report.final_droplet_diameter_m = 0.0;
    } else {
      const double final_diameter_squared =
          std::fma(report.diameter_squared_rate_m2_per_s,
                   report.advanced_duration_s, diameter_squared);
      report.final_droplet_diameter_m =
          std::sqrt(final_diameter_squared);
      report.final_droplet_mass_kg =
          liquid.density_kg_per_m3 * pi / 6.0 *
          final_diameter_squared *
          report.final_droplet_diameter_m;
    }
  }
  if (!std::isfinite(report.final_droplet_mass_kg) ||
      report.final_droplet_mass_kg < 0.0 ||
      !std::isfinite(report.final_droplet_diameter_m) ||
      report.final_droplet_diameter_m < 0.0 ||
      !std::isfinite(report.final_droplet_temperature_k) ||
      !(report.final_droplet_temperature_k > 0.0)) {
    report = {};
    report.status = EvaporationKernelStatus::non_finite_output;
    return report;
  }

  const double parcel_mass_delta =
      (report.final_droplet_mass_kg - input.droplet_mass_kg) *
      input.multiplicity;
  Vector3 parcel_momentum_delta{};
  for (std::size_t component = 0U; component < 3U; ++component) {
    parcel_momentum_delta[component] =
        parcel_mass_delta * input.parcel_velocity_m_per_s[component];
  }
  const double initial_sensible_energy =
      input.droplet_mass_kg * liquid.cp_j_per_kg_k *
      input.droplet_temperature_k;
  const double final_sensible_energy =
      report.final_droplet_mass_kg * liquid.cp_j_per_kg_k *
      report.final_droplet_temperature_k;
  const double parcel_energy_delta =
      (final_sensible_energy - initial_sensible_energy) *
      input.multiplicity;
  report.exchange =
      make_conservative_exchange(parcel_mass_delta, parcel_momentum_delta,
                                 parcel_energy_delta);
  if (!report.exchange.available) {
    report = {};
    report.status = EvaporationKernelStatus::non_finite_output;
    return report;
  }
  report.status = EvaporationKernelStatus::success;
  return report;
}

TabBreakupReport evaluate_tab_breakup(const TabBreakupInput& input) noexcept {
  TabBreakupReport report;
  const TabCoefficients& coefficients = input.coefficients;
  if (!std::isfinite(input.initial_deformation) ||
      !std::isfinite(input.initial_deformation_rate_per_s) ||
      !std::isfinite(input.relative_speed_m_per_s) ||
      input.relative_speed_m_per_s < 0.0 ||
      !std::isfinite(input.gas_density_kg_per_m3) ||
      !(input.gas_density_kg_per_m3 > 0.0) ||
      !std::isfinite(input.liquid_density_kg_per_m3) ||
      !(input.liquid_density_kg_per_m3 > 0.0) ||
      !std::isfinite(input.liquid_viscosity_pa_s) ||
      input.liquid_viscosity_pa_s < 0.0 ||
      !std::isfinite(input.surface_tension_n_per_m) ||
      input.surface_tension_n_per_m < 0.0 ||
      !std::isfinite(input.droplet_radius_m) ||
      !(input.droplet_radius_m > 0.0) ||
      !std::isfinite(input.breakup_threshold) ||
      !(input.breakup_threshold > 0.0) ||
      !std::isfinite(input.duration_s) || input.duration_s < 0.0 ||
      !std::isfinite(coefficients.forcing_coefficient) ||
      coefficients.forcing_coefficient < 0.0 ||
      !std::isfinite(coefficients.breakup_coefficient) ||
      !(coefficients.breakup_coefficient > 0.0) ||
      !std::isfinite(coefficients.stiffness_coefficient) ||
      coefficients.stiffness_coefficient < 0.0 ||
      !std::isfinite(coefficients.damping_coefficient) ||
      coefficients.damping_coefficient < 0.0) {
    return report;
  }

  const double radius_squared =
      input.droplet_radius_m * input.droplet_radius_m;
  const double radius_cubed = radius_squared * input.droplet_radius_m;
  report.forcing_per_s2 =
      coefficients.forcing_coefficient /
      coefficients.breakup_coefficient *
      input.gas_density_kg_per_m3 / input.liquid_density_kg_per_m3 *
      input.relative_speed_m_per_s * input.relative_speed_m_per_s /
      radius_squared;
  report.damping_per_s =
      coefficients.damping_coefficient * input.liquid_viscosity_pa_s /
      (input.liquid_density_kg_per_m3 * radius_squared);
  report.stiffness_per_s2 =
      coefficients.stiffness_coefficient *
      input.surface_tension_n_per_m /
      (input.liquid_density_kg_per_m3 * radius_cubed);
  if (!std::isfinite(report.forcing_per_s2) ||
      report.forcing_per_s2 < 0.0 ||
      !std::isfinite(report.damping_per_s) ||
      report.damping_per_s < 0.0 ||
      !std::isfinite(report.stiffness_per_s2) ||
      report.stiffness_per_s2 < 0.0) {
    report = {};
    report.status = TabKernelStatus::non_finite_output;
    return report;
  }

  const auto state_at = [&](double time_s) noexcept {
    return evaluate_tab_linear_state(
        input.initial_deformation,
        input.initial_deformation_rate_per_s, report.forcing_per_s2,
        report.damping_per_s, report.stiffness_per_s2, time_s);
  };
  const TabLinearState initial{input.initial_deformation,
                               input.initial_deformation_rate_per_s};
  if (std::abs(initial.deformation) >= input.breakup_threshold) {
    report.breakup_requested = true;
    report.advanced_duration_s = 0.0;
    report.event_time_s = 0.0;
    report.candidate =
        {true, initial.deformation, initial.rate_per_s};
    report.status = TabKernelStatus::success;
    return report;
  }
  if (input.duration_s == 0.0) {
    report.candidate =
        {true, initial.deformation, initial.rate_per_s};
    report.status = TabKernelStatus::success;
    return report;
  }

  TabLinearState left_state = initial;
  double left_time = 0.0;
  bool non_finite_state = false;
  const auto process_monotone_boundary = [&](double right_time) noexcept {
    const TabLinearState right_state = state_at(right_time);
    if (!finite_tab_state(right_state)) {
      non_finite_state = true;
      return false;
    }
    if (std::abs(right_state.deformation) < input.breakup_threshold) {
      left_time = right_time;
      left_state = right_state;
      return false;
    }

    const double target = right_state.deformation >= 0.0
                              ? input.breakup_threshold
                              : -input.breakup_threshold;
    double lower = left_time;
    double upper = right_time;
    for (std::uint32_t iteration = 0U; iteration < 80U; ++iteration) {
      const double middle = lower + 0.5 * (upper - lower);
      const TabLinearState middle_state = state_at(middle);
      if (!finite_tab_state(middle_state)) {
        non_finite_state = true;
        return false;
      }
      const bool beyond = target > 0.0
                              ? middle_state.deformation >= target
                              : middle_state.deformation <= target;
      if (beyond) {
        upper = middle;
      } else {
        lower = middle;
      }
    }
    const TabLinearState event_state = state_at(upper);
    if (!finite_tab_state(event_state)) {
      non_finite_state = true;
      return false;
    }
    report.threshold_crossed = true;
    report.breakup_requested = true;
    report.event_time_s = upper;
    report.advanced_duration_s = upper;
    report.candidate = {true, target, event_state.rate_per_s};
    return true;
  };

  const double half_damping = 0.5 * report.damping_per_s;
  const double discriminant =
      report.stiffness_per_s2 - half_damping * half_damping;
  const double regime_scale =
      std::max({1.0, report.stiffness_per_s2,
                half_damping * half_damping});
  const double regime_tolerance =
      32.0 * std::numeric_limits<double>::epsilon() * regime_scale;
  if (discriminant > regime_tolerance) {
    constexpr double pi =
        3.141592653589793238462643383279502884;
    constexpr std::uint32_t maximum_extrema = 65536U;
    const double frequency = std::sqrt(discriminant);
    const double equilibrium =
        report.forcing_per_s2 / report.stiffness_per_s2;
    const double initial_offset =
        input.initial_deformation - equilibrium;
    const double sine_coefficient =
        (input.initial_deformation_rate_per_s +
         half_damping * initial_offset) /
        frequency;
    const double derivative_cosine =
        input.initial_deformation_rate_per_s;
    const double derivative_sine =
        -initial_offset * frequency -
        half_damping * sine_coefficient;
    if (derivative_cosine != 0.0 || derivative_sine != 0.0) {
      const double phase =
          std::atan2(derivative_sine, derivative_cosine);
      const double first_phase = phase + 0.5 * pi;
      double index = std::ceil(-first_phase / pi);
      double extremum_time =
          (first_phase + index * pi) / frequency;
      while (extremum_time <= 0.0) {
        index += 1.0;
        extremum_time =
            (first_phase + index * pi) / frequency;
      }
      std::uint32_t extrema_count = 0U;
      while (extremum_time < input.duration_s) {
        if (extrema_count == maximum_extrema) {
          report = {};
          report.status = TabKernelStatus::event_search_unavailable;
          return report;
        }
        if (process_monotone_boundary(extremum_time)) {
          report.status = TabKernelStatus::success;
          return report;
        }
        if (non_finite_state) {
          report = {};
          report.status = TabKernelStatus::non_finite_output;
          return report;
        }
        ++extrema_count;
        index += 1.0;
        const double next_extremum =
            (first_phase + index * pi) / frequency;
        if (!(next_extremum > extremum_time)) {
          report = {};
          report.status = TabKernelStatus::event_search_unavailable;
          return report;
        }
        extremum_time = next_extremum;
      }
    }
  } else {
    const TabLinearState final_state = state_at(input.duration_s);
    if (!finite_tab_state(final_state)) {
      report = {};
      report.status = TabKernelStatus::non_finite_output;
      return report;
    }
    if (opposite_sign(initial.rate_per_s, final_state.rate_per_s)) {
      double lower = 0.0;
      double upper = input.duration_s;
      double lower_rate = initial.rate_per_s;
      for (std::uint32_t iteration = 0U; iteration < 80U; ++iteration) {
        const double middle = lower + 0.5 * (upper - lower);
        const TabLinearState middle_state = state_at(middle);
        if (!finite_tab_state(middle_state)) {
          report = {};
          report.status = TabKernelStatus::non_finite_output;
          return report;
        }
        if (middle_state.rate_per_s == 0.0) {
          lower = middle;
          upper = middle;
          break;
        }
        if (opposite_sign(lower_rate, middle_state.rate_per_s)) {
          upper = middle;
        } else {
          lower = middle;
          lower_rate = middle_state.rate_per_s;
        }
      }
      const double extremum_time = lower + 0.5 * (upper - lower);
      if (process_monotone_boundary(extremum_time)) {
        report.status = TabKernelStatus::success;
        return report;
      }
      if (non_finite_state) {
        report = {};
        report.status = TabKernelStatus::non_finite_output;
        return report;
      }
    }
  }

  if (process_monotone_boundary(input.duration_s)) {
    report.status = TabKernelStatus::success;
    return report;
  }
  if (non_finite_state) {
    report = {};
    report.status = TabKernelStatus::non_finite_output;
    return report;
  }
  report.advanced_duration_s = input.duration_s;
  report.candidate =
      {true, left_state.deformation, left_state.rate_per_s};
  report.status = TabKernelStatus::success;
  return report;
}

}  // namespace hundun::v04::spray
