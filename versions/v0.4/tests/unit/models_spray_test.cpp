// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_spray.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04::spray;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool test_parcel_gas_exchange_conservation() {
  const ParcelGasExchangeCandidate exchange = make_conservative_exchange(
      -0.25, {1.0, -2.0, 0.5}, 3.5);
  const ExchangeConservationReport report =
      evaluate_exchange_conservation(exchange);
  bool passed = true;
  passed &= expect(exchange.available &&
                       exchange.gas_mass_delta_kg == 0.25 &&
                       exchange.gas_momentum_delta_kg_m_per_s[0U] == -1.0 &&
                       exchange.gas_momentum_delta_kg_m_per_s[1U] == 2.0 &&
                       exchange.gas_momentum_delta_kg_m_per_s[2U] == -0.5 &&
                       exchange.gas_energy_delta_j == -3.5,
                   "gas delta is the exact opposite parcel delta");
  passed &= expect(report.status == ExchangeConservationStatus::success &&
                       report.mass_residual_kg == 0.0 &&
                       report.momentum_residual_kg_m_per_s ==
                           Vector3{0.0, 0.0, 0.0} &&
                       report.energy_residual_j == 0.0,
                   "closed parcel+gas totals have exact zero residuals");

  ParcelGasExchangeCandidate tampered = exchange;
  tampered.gas_mass_delta_kg += 1.0e-3;
  passed &= expect(evaluate_exchange_conservation(tampered).status ==
                       ExchangeConservationStatus::conservation_failure,
                   "a one-sided exchange mutation is rejected");
  return passed;
}

bool test_stable_ids_and_pure_rng_domain() {
  const ParcelRandomAddress zero_address{
      0U, 0U, 0U, 0U, ParcelRandomPurpose::stable_id};
  const std::uint64_t golden = parcel_random_u64(zero_address, 0U);
  const ParcelId first = make_stable_parcel_id(zero_address);
  const ParcelId repeated = make_stable_parcel_id(zero_address);

  ParcelRandomAddress next_step = zero_address;
  next_step.accepted_step = 1U;
  ParcelRandomAddress injection = zero_address;
  injection.purpose = ParcelRandomPurpose::injection_direction;

  bool passed = true;
  passed &= expect(parcel_rng_model_id() ==
                       "parcel_counter_splitmix64_v1" &&
                       golden == UINT64_C(0xe220a8397b1dcdaf),
                   "parcel RNG identity and published zero vector are frozen");
  passed &= expect(first == repeated && first.high == golden &&
                       (first.high != 0U || first.low != 0U),
                   "stable ID is deterministic and never the null ID");
  passed &= expect(parcel_random_u64(zero_address, 1U) != golden &&
                       parcel_random_u64(next_step, 0U) != golden &&
                       parcel_random_u64(injection, 0U) != golden,
                   "lane, accepted step and purpose are independent domains");
  const double uniform = parcel_uniform_01(zero_address, 2U);
  passed &= expect(uniform >= 0.0 && uniform < 1.0 &&
                       uniform == parcel_uniform_01(zero_address, 2U),
                   "uniform variate is stateless and half-open");

  SprayParcelState parcel;
  parcel.id = first;
  parcel.position_m = {1.0, 2.0, 3.0};
  parcel.velocity_m_per_s = {4.0, 5.0, 6.0};
  parcel.droplet_mass_kg = 2.0e-9;
  parcel.droplet_diameter_m = 1.0e-4;
  parcel.multiplicity = 50.0;
  parcel.temperature_k = 350.0;
  parcel.liquid_material_fingerprint = 91U;
  parcel.owner_global_cell = 17U;
  parcel.age_s = 0.01;
  passed &= expect(validate_parcel_state(parcel) == ParcelStateStatus::success,
                   "positive finite parcel value type is valid");
  parcel.multiplicity = 0.0;
  passed &= expect(validate_parcel_state(parcel) ==
                       ParcelStateStatus::invalid_input,
                   "zero multiplicity is rejected without state mutation");
  return passed;
}

class FakeLiquidProperties final : public LiquidPropertyProvider {
 public:
  LiquidPropertyReport evaluate(const LiquidPropertyQuery& query) const
      noexcept override {
    if (query.material_fingerprint != 91U) {
      return {LiquidPropertyStatus::unknown_material, 0U, 0.0, {}};
    }
    return {LiquidPropertyStatus::success,
            query.material_fingerprint,
            query.temperature_k,
            {750.0, 2200.0, 2.5e5, 2.0e4, 0.025, 8.0e-4}};
  }
};

bool test_liquid_property_interface_contract() {
  const FakeLiquidProperties provider;
  const LiquidPropertyQuery query{91U, 350.0};
  const LiquidPropertyReport success = provider.evaluate(query);
  bool passed = true;
  passed &= expect(validate_liquid_property_report(query, success) ==
                       LiquidPropertyStatus::success &&
                       success.properties.density_kg_per_m3 == 750.0 &&
                       success.properties.cp_j_per_kg_k == 2200.0 &&
                       success.properties.latent_heat_j_per_kg == 2.5e5,
                   "fingerprint query returns finite SI liquid properties");

  const LiquidPropertyQuery unknown_query{92U, 350.0};
  const LiquidPropertyReport unknown = provider.evaluate(unknown_query);
  passed &= expect(
      validate_liquid_property_report(unknown_query, unknown) ==
              LiquidPropertyStatus::unknown_material &&
          unknown.material_fingerprint == 0U &&
          unknown.evaluated_temperature_k == 0.0 &&
          unknown.properties == LiquidProperties{},
      "unavailable material has explicit status and canonical zero payload");

  LiquidPropertyReport malformed = success;
  malformed.properties.surface_tension_n_per_m = -1.0;
  passed &= expect(validate_liquid_property_report(query, malformed) ==
                       LiquidPropertyStatus::provider_contract_failure,
                   "invalid successful property payload is rejected");
  return passed;
}

bool test_schiller_naumann_drag_limits_and_exchange() {
  SchillerNaumannDragInput zero;
  zero.gas_velocity_m_per_s = {4.0, -2.0, 1.0};
  zero.parcel_velocity_m_per_s = zero.gas_velocity_m_per_s;
  zero.gas_density_kg_per_m3 = 1.2;
  zero.gas_dynamic_viscosity_pa_s = 1.8e-5;
  zero.droplet_diameter_m = 1.0e-4;
  zero.droplet_mass_kg = 4.0e-10;
  zero.multiplicity = 25.0;
  zero.duration_s = 2.0e-4;
  const SchillerNaumannDragReport zero_report =
      evaluate_schiller_naumann_drag(zero);

  bool passed = true;
  passed &= expect(
      zero_report.succeeded() &&
          zero_report.model_id == "schiller_naumann_v1" &&
          zero_report.reynolds_number == 0.0 &&
          zero_report.force_on_one_droplet_n == Vector3{0.0, 0.0, 0.0} &&
          zero_report.acceleration_m_per_s2 == Vector3{0.0, 0.0, 0.0} &&
          zero_report.exchange.available &&
          zero_report.exchange.parcel_momentum_delta_kg_m_per_s ==
              Vector3{0.0, 0.0, 0.0} &&
          zero_report.exchange.parcel_energy_delta_j == 0.0,
      "zero slip has an exact finite zero drag limit");

  SchillerNaumannDragInput stokes = zero;
  stokes.gas_velocity_m_per_s = {1.0e-10, 0.0, 0.0};
  stokes.parcel_velocity_m_per_s = {0.0, 0.0, 0.0};
  stokes.gas_density_kg_per_m3 = 1.0;
  stokes.gas_dynamic_viscosity_pa_s = 1.0;
  stokes.droplet_diameter_m = 1.0;
  stokes.droplet_mass_kg = 2.0;
  stokes.multiplicity = 3.0;
  stokes.duration_s = 0.5;
  const SchillerNaumannDragReport stokes_report =
      evaluate_schiller_naumann_drag(stokes);
  constexpr double pi = 3.141592653589793238462643383279502884;
  const double expected_stokes_force = 3.0 * pi * 1.0e-10;
  passed &= expect(
      stokes_report.succeeded() && stokes_report.used_stokes_limit &&
          std::abs(stokes_report.force_on_one_droplet_n[0U] -
                   expected_stokes_force) < 1.0e-24 &&
          std::abs(stokes_report.acceleration_m_per_s2[0U] -
                   expected_stokes_force / 2.0) < 1.0e-24 &&
          std::abs(
              stokes_report.exchange.parcel_momentum_delta_kg_m_per_s[0U] -
              expected_stokes_force * 0.5 * 3.0) < 1.0e-24,
      "low Reynolds drag is the Stokes force and multiplicity is applied once");

  SchillerNaumannDragInput finite = zero;
  finite.gas_velocity_m_per_s = {3.0, 4.0, 0.0};
  finite.parcel_velocity_m_per_s = {0.0, 0.0, 0.0};
  finite.gas_density_kg_per_m3 = 1.2;
  finite.gas_dynamic_viscosity_pa_s = 1.8e-5;
  finite.droplet_diameter_m = 1.0e-3;
  finite.droplet_mass_kg = 8.0e-7;
  finite.multiplicity = 7.0;
  finite.duration_s = 1.0e-3;
  const SchillerNaumannDragReport finite_report =
      evaluate_schiller_naumann_drag(finite);
  const double expected_reynolds = 1.2 * 5.0 * 1.0e-3 / 1.8e-5;
  const double expected_cd =
      24.0 / expected_reynolds *
      (1.0 + 0.15 * std::pow(expected_reynolds, 0.687));
  const double expected_force_magnitude =
      0.5 * expected_cd * 1.2 * (pi * 1.0e-6 / 4.0) * 25.0;
  passed &= expect(
      finite_report.succeeded() && !finite_report.used_stokes_limit &&
          !finite_report.used_constant_drag_branch &&
          std::abs(finite_report.reynolds_number - expected_reynolds) <
              1.0e-12 &&
          std::abs(finite_report.drag_coefficient - expected_cd) < 1.0e-14 &&
          std::abs(finite_report.force_on_one_droplet_n[0U] -
                   expected_force_magnitude * 0.6) < 1.0e-14 &&
          std::abs(finite_report.force_on_one_droplet_n[1U] -
                   expected_force_magnitude * 0.8) < 1.0e-14,
      "finite Reynolds Schiller--Naumann reference value and direction match");
  passed &= expect(
      evaluate_exchange_conservation(finite_report.exchange).status ==
          ExchangeConservationStatus::success,
      "drag candidate is a closed parcel+gas momentum and energy exchange");

  SchillerNaumannDragInput shifted_rotated = finite;
  shifted_rotated.gas_velocity_m_per_s = {15.0, -10.0, 8.0};
  shifted_rotated.parcel_velocity_m_per_s = {11.0, -7.0, 8.0};
  const SchillerNaumannDragReport shifted_report =
      evaluate_schiller_naumann_drag(shifted_rotated);
  passed &= expect(
      shifted_report.succeeded() &&
          shifted_report.reynolds_number == finite_report.reynolds_number &&
          std::abs(shifted_report.force_on_one_droplet_n[0U] -
                   finite_report.force_on_one_droplet_n[1U]) < 1.0e-14 &&
          std::abs(shifted_report.force_on_one_droplet_n[1U] +
                   finite_report.force_on_one_droplet_n[0U]) < 1.0e-14 &&
          shifted_report.force_on_one_droplet_n[2U] == 0.0,
      "drag force is Galilean invariant and rotates with slip");

  SchillerNaumannDragInput high_re = finite;
  high_re.gas_velocity_m_per_s = {100.0, 0.0, 0.0};
  const SchillerNaumannDragReport high_re_report =
      evaluate_schiller_naumann_drag(high_re);
  passed &= expect(high_re_report.succeeded() &&
                       high_re_report.used_constant_drag_branch &&
                       high_re_report.drag_coefficient == 0.44,
                   "Reynolds number above 1000 uses the declared constant branch");

  SchillerNaumannDragInput invalid = finite;
  invalid.gas_dynamic_viscosity_pa_s =
      std::numeric_limits<double>::quiet_NaN();
  const SchillerNaumannDragReport invalid_report =
      evaluate_schiller_naumann_drag(invalid);
  passed &= expect(!invalid_report.succeeded() &&
                       !invalid_report.exchange.available &&
                       invalid_report.force_on_one_droplet_n ==
                           Vector3{0.0, 0.0, 0.0},
                   "invalid drag input returns no partial exchange");
  return passed;
}

bool test_ranz_marshall_heat_mass_transfer() {
  RanzMarshallInput stagnant;
  stagnant.reynolds_number = 0.0;
  stagnant.prandtl_number = 0.7;
  stagnant.schmidt_number = 1.1;
  stagnant.gas_thermal_conductivity_w_per_m_k = 0.08;
  stagnant.vapor_diffusivity_m2_per_s = 2.0e-5;
  stagnant.droplet_diameter_m = 1.0e-3;
  const RanzMarshallReport stagnant_report =
      evaluate_ranz_marshall(stagnant);
  bool passed = true;
  passed &= expect(
      stagnant_report.succeeded() &&
          stagnant_report.model_id == "ranz_marshall_v1" &&
          stagnant_report.nusselt_number == 2.0 &&
          stagnant_report.sherwood_number == 2.0 &&
          stagnant_report.heat_transfer_coefficient_w_per_m2_k == 160.0 &&
          stagnant_report.mass_transfer_coefficient_m_per_s == 0.04,
      "Ranz--Marshall has exact stagnant-sphere limits and SI coefficients");

  RanzMarshallInput independent = stagnant;
  independent.reynolds_number = 16.0;
  independent.prandtl_number = 8.0;
  independent.schmidt_number = 27.0;
  independent.gas_thermal_conductivity_w_per_m_k = 0.1;
  independent.vapor_diffusivity_m2_per_s = 3.0e-5;
  independent.droplet_diameter_m = 0.02;
  const RanzMarshallReport independent_report =
      evaluate_ranz_marshall(independent);
  passed &= expect(
      independent_report.succeeded() &&
          std::abs(independent_report.nusselt_number - 6.8) < 1.0e-14 &&
          std::abs(independent_report.sherwood_number - 9.2) < 1.0e-14 &&
          std::abs(
              independent_report.heat_transfer_coefficient_w_per_m2_k -
              34.0) < 1.0e-13 &&
          std::abs(independent_report.mass_transfer_coefficient_m_per_s -
                   0.0138) < 1.0e-15,
      "Prandtl and Schmidt channels use independent cube-root factors");

  RanzMarshallInput invalid = independent;
  invalid.schmidt_number = 0.0;
  const RanzMarshallReport invalid_report =
      evaluate_ranz_marshall(invalid);
  passed &= expect(!invalid_report.succeeded() &&
                       invalid_report.nusselt_number == 0.0 &&
                       invalid_report.sherwood_number == 0.0 &&
                       invalid_report.mass_transfer_coefficient_m_per_s ==
                           0.0,
                   "invalid transfer input has explicit status and zero payload");
  return passed;
}

SingleComponentEvaporationInput make_evaporation_input() {
  constexpr double pi = 3.141592653589793238462643383279502884;
  SingleComponentEvaporationInput input;
  input.parcel_velocity_m_per_s = {2.0, -1.0, 0.5};
  input.droplet_diameter_m = 1.0e-3;
  input.liquid_properties =
      {800.0, 2000.0, 2.0e5, 0.0, 0.02, 1.0e-3};
  input.droplet_mass_kg =
      input.liquid_properties.density_kg_per_m3 * pi / 6.0 *
      input.droplet_diameter_m * input.droplet_diameter_m *
      input.droplet_diameter_m;
  input.multiplicity = 10.0;
  input.droplet_temperature_k = 350.0;
  input.gas_temperature_k = 450.0;
  input.gas_pressure_pa = 1.0e5;
  input.gas_density_kg_per_m3 = 1.0;
  input.gas_vapor_mass_fraction = 0.0;
  input.vapor_molecular_weight_kg_per_kmol = 100.0;
  input.carrier_molecular_weight_kg_per_kmol = 29.0;
  input.reynolds_number = 0.0;
  input.prandtl_number = 1.0;
  input.schmidt_number = 1.0;
  input.gas_thermal_conductivity_w_per_m_k = 0.1;
  input.vapor_diffusivity_m2_per_s = 1.0e-5;
  input.duration_s = 1.0e-3;
  return input;
}

bool test_single_component_evaporation_limits_and_conservation() {
  constexpr double pi = 3.141592653589793238462643383279502884;
  const SingleComponentEvaporationInput heating = make_evaporation_input();
  const SingleComponentEvaporationReport heating_report =
      evaluate_single_component_evaporation(heating);
  const double expected_heat_rate =
      pi * heating.droplet_diameter_m *
      heating.gas_thermal_conductivity_w_per_m_k * 2.0 *
      (heating.gas_temperature_k - heating.droplet_temperature_k);
  const double expected_temperature =
      heating.droplet_temperature_k +
      expected_heat_rate /
          (heating.droplet_mass_kg *
           heating.liquid_properties.cp_j_per_kg_k) *
          heating.duration_s;
  bool passed = true;
  passed &= expect(
      heating_report.succeeded() &&
          heating_report.model_id ==
              "single_component_spalding_ranz_marshall_v1" &&
          !heating_report.evaporation_active &&
          !heating_report.has_terminal_event &&
          heating_report.evaporation_mass_rate_kg_per_s == 0.0 &&
          heating_report.diameter_squared_rate_m2_per_s == 0.0 &&
          heating_report.final_droplet_mass_kg ==
              heating.droplet_mass_kg &&
          heating_report.final_droplet_diameter_m ==
              heating.droplet_diameter_m &&
          std::abs(heating_report.final_droplet_temperature_k -
                   expected_temperature) < 1.0e-12,
      "zero saturation pressure degenerates exactly to pure droplet heating");
  passed &= expect(
      evaluate_exchange_conservation(heating_report.exchange).status ==
              ExchangeConservationStatus::success &&
          heating_report.exchange.parcel_mass_delta_kg == 0.0 &&
          std::abs(heating_report.exchange.parcel_energy_delta_j -
                   expected_heat_rate * heating.duration_s *
                       heating.multiplicity) < 1.0e-12,
      "heating candidate is an exactly paired parcel+gas energy exchange");

  SingleComponentEvaporationInput evaporating = make_evaporation_input();
  evaporating.droplet_temperature_k = 500.0;
  evaporating.gas_temperature_k = 500.0;
  evaporating.liquid_properties.latent_heat_j_per_kg = 0.0;
  evaporating.liquid_properties.saturation_pressure_pa = 2.0e4;
  const double surface_mole_fraction = 0.2;
  const double surface_mass_fraction =
      surface_mole_fraction *
      evaporating.vapor_molecular_weight_kg_per_kmol /
      (surface_mole_fraction *
           evaporating.vapor_molecular_weight_kg_per_kmol +
       (1.0 - surface_mole_fraction) *
           evaporating.carrier_molecular_weight_kg_per_kmol);
  const double spalding =
      surface_mass_fraction / (1.0 - surface_mass_fraction);
  const double expected_d2_loss_rate =
      4.0 * evaporating.gas_density_kg_per_m3 *
      evaporating.vapor_diffusivity_m2_per_s * 2.0 *
      std::log1p(spalding) /
      evaporating.liquid_properties.density_kg_per_m3;
  const double expected_event_time =
      evaporating.droplet_diameter_m *
      evaporating.droplet_diameter_m / expected_d2_loss_rate;
  evaporating.duration_s = 0.25 * expected_event_time;
  const SingleComponentEvaporationReport d2_report =
      evaluate_single_component_evaporation(evaporating);
  const double expected_final_d2 =
      evaporating.droplet_diameter_m *
          evaporating.droplet_diameter_m -
      expected_d2_loss_rate * evaporating.duration_s;
  passed &= expect(
      d2_report.succeeded() && d2_report.evaporation_active &&
          d2_report.has_terminal_event &&
          !d2_report.complete_evaporation &&
          std::abs(d2_report.surface_vapor_mass_fraction -
                   surface_mass_fraction) < 1.0e-15 &&
          std::abs(d2_report.spalding_mass_number - spalding) < 1.0e-15 &&
          std::abs(d2_report.diameter_squared_rate_m2_per_s +
                   expected_d2_loss_rate) < 1.0e-20 &&
          std::abs(d2_report.event_time_s - expected_event_time) <
              1.0e-12 &&
          std::abs(d2_report.final_droplet_diameter_m *
                       d2_report.final_droplet_diameter_m -
                   expected_final_d2) < 1.0e-20,
      "single-component evaporation follows the analytic frozen d-squared law");

  evaporating.duration_s = 2.0 * expected_event_time;
  const SingleComponentEvaporationReport complete_report =
      evaluate_single_component_evaporation(evaporating);
  passed &= expect(
      complete_report.succeeded() && complete_report.complete_evaporation &&
          complete_report.advanced_duration_s ==
              complete_report.event_time_s &&
          complete_report.final_droplet_mass_kg == 0.0 &&
          complete_report.final_droplet_diameter_m == 0.0 &&
          complete_report.exchange.parcel_mass_delta_kg ==
              -evaporating.droplet_mass_kg * evaporating.multiplicity &&
          complete_report.exchange.parcel_momentum_delta_kg_m_per_s[0U] ==
              complete_report.exchange.parcel_mass_delta_kg *
                  evaporating.parcel_velocity_m_per_s[0U] &&
          evaluate_exchange_conservation(complete_report.exchange).status ==
              ExchangeConservationStatus::success,
      "terminal event reaches exact zero with closed mass, momentum and energy");

  SingleComponentEvaporationInput zero_duration = evaporating;
  zero_duration.duration_s = 0.0;
  const SingleComponentEvaporationReport zero_duration_report =
      evaluate_single_component_evaporation(zero_duration);
  passed &= expect(
      zero_duration_report.succeeded() &&
          zero_duration_report.advanced_duration_s == 0.0 &&
          zero_duration_report.final_droplet_mass_kg ==
              zero_duration.droplet_mass_kg &&
          zero_duration_report.final_droplet_diameter_m ==
              zero_duration.droplet_diameter_m &&
          zero_duration_report.final_droplet_temperature_k ==
              zero_duration.droplet_temperature_k &&
          zero_duration_report.exchange.parcel_mass_delta_kg == 0.0,
      "zero duration preserves the input droplet exactly");

  SingleComponentEvaporationInput unavailable = evaporating;
  unavailable.liquid_properties.saturation_pressure_pa =
      unavailable.gas_pressure_pa;
  const SingleComponentEvaporationReport unavailable_report =
      evaluate_single_component_evaporation(unavailable);
  passed &= expect(
      unavailable_report.status ==
              EvaporationKernelStatus::phase_equilibrium_unavailable &&
          !unavailable_report.exchange.available &&
          unavailable_report.final_droplet_mass_kg == 0.0,
      "saturated pure-vapor surface is explicit failure, never a clipped result");
  return passed;
}

TabBreakupInput make_tab_input() {
  TabBreakupInput input;
  input.initial_deformation = 0.2;
  input.initial_deformation_rate_per_s = 0.3;
  input.relative_speed_m_per_s = 0.0;
  input.gas_density_kg_per_m3 = 1.0;
  input.liquid_density_kg_per_m3 = 1.0;
  input.liquid_viscosity_pa_s = 0.0;
  input.surface_tension_n_per_m = 4.0;
  input.droplet_radius_m = 1.0;
  input.breakup_threshold = 10.0;
  input.duration_s = 0.4;
  input.coefficients = {1.0, 1.0, 1.0, 1.0};
  return input;
}

bool test_tab_breakup_analytic_limits_and_event() {
  const TabBreakupInput harmonic = make_tab_input();
  const TabBreakupReport harmonic_report = evaluate_tab_breakup(harmonic);
  const double omega = 2.0;
  const double expected_y =
      harmonic.initial_deformation *
          std::cos(omega * harmonic.duration_s) +
      harmonic.initial_deformation_rate_per_s / omega *
          std::sin(omega * harmonic.duration_s);
  const double expected_rate =
      -harmonic.initial_deformation * omega *
          std::sin(omega * harmonic.duration_s) +
      harmonic.initial_deformation_rate_per_s *
          std::cos(omega * harmonic.duration_s);
  bool passed = true;
  passed &= expect(
      harmonic_report.succeeded() &&
          harmonic_report.model_id == "tab_linear_oscillator_v1" &&
          harmonic_report.forcing_per_s2 == 0.0 &&
          harmonic_report.damping_per_s == 0.0 &&
          harmonic_report.stiffness_per_s2 == 4.0 &&
          std::abs(harmonic_report.candidate.deformation - expected_y) <
              1.0e-14 &&
          std::abs(harmonic_report.candidate.deformation_rate_per_s -
                   expected_rate) < 1.0e-14 &&
          !harmonic_report.breakup_requested,
      "unforced undamped TAB state is the analytic harmonic oscillator");

  TabBreakupInput damped = make_tab_input();
  damped.initial_deformation = 1.0;
  damped.initial_deformation_rate_per_s = 0.0;
  damped.liquid_viscosity_pa_s = 1.0;
  damped.duration_s = 0.3;
  const TabBreakupReport damped_report = evaluate_tab_breakup(damped);
  const double decay = 0.5;
  const double damped_omega = std::sqrt(4.0 - decay * decay);
  const double expected_damped_y =
      std::exp(-decay * damped.duration_s) *
      (std::cos(damped_omega * damped.duration_s) +
       decay / damped_omega *
           std::sin(damped_omega * damped.duration_s));
  passed &= expect(
      damped_report.succeeded() && damped_report.damping_per_s == 1.0 &&
          std::abs(damped_report.candidate.deformation -
                   expected_damped_y) < 1.0e-14,
      "TAB damping has the analytic underdamped decay and restoring sign");

  TabBreakupInput equilibrium = make_tab_input();
  equilibrium.initial_deformation = 0.5;
  equilibrium.initial_deformation_rate_per_s = 0.0;
  equilibrium.relative_speed_m_per_s = 1.0;
  equilibrium.gas_density_kg_per_m3 = 2.0;
  equilibrium.duration_s = 3.0;
  const TabBreakupReport equilibrium_report =
      evaluate_tab_breakup(equilibrium);
  passed &= expect(
      equilibrium_report.succeeded() &&
          equilibrium_report.forcing_per_s2 == 2.0 &&
          equilibrium_report.stiffness_per_s2 == 4.0 &&
          equilibrium_report.candidate.deformation == 0.5 &&
          equilibrium_report.candidate.deformation_rate_per_s == 0.0,
      "forcing/stiffness equilibrium remains exactly stationary");

  TabBreakupInput interior_event = make_tab_input();
  interior_event.initial_deformation = 0.0;
  interior_event.initial_deformation_rate_per_s = 1.0;
  interior_event.surface_tension_n_per_m = 1.0;
  interior_event.breakup_threshold = 0.9;
  constexpr double pi = 3.141592653589793238462643383279502884;
  interior_event.duration_s = pi;
  const TabBreakupReport event_report =
      evaluate_tab_breakup(interior_event);
  const double expected_event_time = std::asin(0.9);
  passed &= expect(
      event_report.succeeded() && event_report.breakup_requested &&
          event_report.threshold_crossed &&
          std::abs(event_report.event_time_s - expected_event_time) <
              1.0e-12 &&
          event_report.advanced_duration_s == event_report.event_time_s &&
          std::abs(event_report.candidate.deformation - 0.9) < 1.0e-14,
      "TAB finds the first interior threshold crossing even when endpoint is safe");

  TabBreakupInput full = make_tab_input();
  full.duration_s = 0.7;
  const TabBreakupReport full_report = evaluate_tab_breakup(full);
  TabBreakupInput first = full;
  first.duration_s = 0.3;
  const TabBreakupReport first_report = evaluate_tab_breakup(first);
  TabBreakupInput second = full;
  second.initial_deformation = first_report.candidate.deformation;
  second.initial_deformation_rate_per_s =
      first_report.candidate.deformation_rate_per_s;
  second.duration_s = 0.4;
  const TabBreakupReport second_report = evaluate_tab_breakup(second);
  passed &= expect(
      full_report.succeeded() && first_report.succeeded() &&
          second_report.succeeded() &&
          std::abs(full_report.candidate.deformation -
                   second_report.candidate.deformation) < 1.0e-14 &&
          std::abs(full_report.candidate.deformation_rate_per_s -
                   second_report.candidate.deformation_rate_per_s) <
              1.0e-14,
      "analytic TAB advance is consistent under timestep subdivision");

  TabBreakupInput zero_duration = make_tab_input();
  zero_duration.duration_s = 0.0;
  const TabBreakupReport zero_report = evaluate_tab_breakup(zero_duration);
  passed &= expect(
      zero_report.succeeded() && zero_report.advanced_duration_s == 0.0 &&
          zero_report.candidate.deformation ==
              zero_duration.initial_deformation &&
          zero_report.candidate.deformation_rate_per_s ==
              zero_duration.initial_deformation_rate_per_s,
      "zero-duration TAB candidate is the exact input state");

  TabBreakupInput invalid = make_tab_input();
  invalid.droplet_radius_m =
      std::numeric_limits<double>::infinity();
  const TabBreakupReport invalid_report = evaluate_tab_breakup(invalid);
  passed &= expect(!invalid_report.succeeded() &&
                       !invalid_report.candidate.available &&
                       invalid_report.forcing_per_s2 == 0.0,
                   "invalid TAB input returns no trial state");
  return passed;
}

}  // namespace

int main() {
  bool passed = test_parcel_gas_exchange_conservation();
  passed &= test_stable_ids_and_pure_rng_domain();
  passed &= test_liquid_property_interface_contract();
  passed &= test_schiller_naumann_drag_limits_and_exchange();
  passed &= test_ranz_marshall_heat_mass_transfer();
  passed &= test_single_component_evaporation_limits_and_conservation();
  passed &= test_tab_breakup_analytic_limits_and_event();
  return passed ? 0 : 1;
}
