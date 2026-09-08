// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_transfer_detail.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

constexpr double kPi = 3.141592653589793238462643383279502884;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool near(double left, double right, double relative = 1.0e-10,
          double absolute = 1.0e-14) {
  return std::abs(left - right) <=
         absolute + relative * (std::abs(left) + std::abs(right));
}

double dot(const Vector3& left, const Vector3& right) {
  return left[0U] * right[0U] + left[1U] * right[1U] +
         left[2U] * right[2U];
}

bool canonical_failure(const LiquidPropertyReport& report) {
  const double values[]{report.evaluated_temperature_k,
                        report.properties.density_kg_per_m3,
                        report.properties.cp_j_per_kg_k,
                        report.properties.latent_heat_j_per_kg,
                        report.properties.saturation_pressure_pa,
                        report.properties.surface_tension_n_per_m,
                        report.properties.viscosity_pa_s};
  bool canonical = report.material_fingerprint == 0U;
  for (double value : values) canonical &= value == 0.0 && !std::signbit(value);
  return canonical;
}

bool canonical_parcel(const SprayParcelState& state) {
  return state.id == ParcelId{} && state.position_m == Vector3{} &&
         state.velocity_m_per_s == Vector3{} &&
         state.droplet_mass_kg == 0.0 &&
         !std::signbit(state.droplet_mass_kg) &&
         state.droplet_diameter_m == 0.0 &&
         !std::signbit(state.droplet_diameter_m) &&
         state.multiplicity == 0.0 && !std::signbit(state.multiplicity) &&
         state.temperature_k == 0.0 && !std::signbit(state.temperature_k) &&
         state.liquid_material_fingerprint == 0U &&
         state.owner_global_cell == 0U && state.age_s == 0.0 &&
         !std::signbit(state.age_s);
}

TemperatureCorrelation constant(double value) {
  TemperatureCorrelation correlation;
  correlation.kind = TemperatureCorrelationKind::constant;
  correlation.c[0U] = value;
  return correlation;
}

LiquidPropertyPack base_pack(std::uint64_t fingerprint) {
  LiquidPropertyPack pack;
  pack.material_fingerprint = fingerprint;
  pack.minimum_temperature_k = 250.0;
  pack.maximum_temperature_k = 650.0;
  pack.density_kg_per_m3 = constant(750.0);
  pack.cp_j_per_kg_k = constant(2200.0);
  pack.latent_heat_j_per_kg = constant(2.5e5);
  pack.surface_tension_n_per_m = constant(0.025);
  pack.viscosity_pa_s = constant(8.0e-4);
  pack.saturation_pressure.kind =
      SaturationPressureCorrelationKind::antoine_kelvin;
  pack.saturation_pressure.antoine_a = 4.0;
  pack.saturation_pressure.antoine_b_k = 0.0;
  pack.saturation_pressure.antoine_c_k = 0.0;
  pack.saturation_pressure.pressure_scale_pa = 1.0;
  return pack;
}

SprayParcelState parcel(double multiplicity = 4.0) {
  SprayParcelState state;
  state.id = {1U, 2U};
  state.velocity_m_per_s = {1.0, -0.5, 0.25};
  state.droplet_diameter_m = 1.0e-4;
  state.droplet_mass_kg =
      750.0 * kPi / 6.0 * std::pow(state.droplet_diameter_m, 3.0);
  state.multiplicity = multiplicity;
  state.temperature_k = 350.0;
  state.liquid_material_fingerprint = 91U;
  return state;
}

OneThirdFilmInput film_input(double surface_y = 0.2,
                             double far_y = 0.02) {
  OneThirdFilmInput input;
  input.far_gas_velocity_m_per_s = {4.0, -0.5, 0.25};
  input.droplet_surface_temperature_k = 350.0;
  input.far_gas_temperature_k = 950.0;
  input.surface_vapor_mass_fraction = surface_y;
  input.far_gas_vapor_mass_fraction = far_y;
  input.thermodynamic_pressure_pa = 101325.0;
  input.reference_density_kg_per_m3 = 0.8;
  input.reference_dynamic_viscosity_pa_s = 2.0e-5;
  input.reference_thermal_conductivity_w_per_m_k = 0.05;
  input.reference_vapor_diffusivity_m2_per_s = 2.0e-5;
  input.reference_cp_j_per_kg_k = 1100.0;
  return input;
}

bool test_liquid_property_dispatch_and_correlations() {
  LiquidPropertyPack packs[2U]{base_pack(91U), base_pack(92U)};
  packs[0U].density_kg_per_m3.kind =
      TemperatureCorrelationKind::polynomial_cubic;
  packs[0U].density_kg_per_m3.reference_temperature_k = 300.0;
  packs[0U].density_kg_per_m3.c = {800.0, -0.5, 1.0e-3, 0.0};
  packs[1U].saturation_pressure.kind =
      SaturationPressureCorrelationKind::clausius_clapeyron;
  packs[1U].saturation_pressure.reference_pressure_pa = 101325.0;
  packs[1U].saturation_pressure.reference_temperature_k = 400.0;
  packs[1U].saturation_pressure.latent_heat_j_per_kg = 2.5e5;
  packs[1U].saturation_pressure.molecular_weight_kg_per_kmol = 100.0;
  const LiquidPropertyService service(packs, 2U);

  const LiquidPropertyReport polynomial = service.evaluate({91U, 350.0});
  const LiquidPropertyReport cc_reference = service.evaluate({92U, 400.0});
  const LiquidPropertyReport boundary = service.evaluate({91U, 250.0});
  const LiquidPropertyReport unknown = service.evaluate({93U, 350.0});
  const LiquidPropertyReport outside = service.evaluate({91U, 249.999});

  bool passed = true;
  passed &= expect(polynomial.succeeded() &&
                       near(polynomial.properties.density_kg_per_m3, 777.5) &&
                       polynomial.properties.saturation_pressure_pa == 1.0e4,
                   "fingerprint dispatch evaluates polynomial and Kelvin Antoine law");
  passed &= expect(cc_reference.succeeded() &&
                       near(cc_reference.properties.saturation_pressure_pa,
                            101325.0) &&
                       boundary.succeeded(),
                   "Clausius reference and inclusive validity boundary are exact");
  passed &= expect(unknown.status == LiquidPropertyStatus::unknown_material &&
                       canonical_failure(unknown) &&
                       outside.status ==
                           LiquidPropertyStatus::temperature_out_of_range &&
                       canonical_failure(outside),
                   "unknown and out-of-range queries have explicit canonical failures");

  LiquidPropertyPack invalid = base_pack(77U);
  invalid.cp_j_per_kg_k = constant(-1.0);
  const LiquidPropertyService invalid_service(&invalid, 1U);
  const LiquidPropertyReport invalid_report =
      invalid_service.evaluate({77U, 350.0});
  passed &= expect(invalid_report.status ==
                           LiquidPropertyStatus::correlation_domain_error &&
                       canonical_failure(invalid_report),
                   "negative finite property is rejected without partial payload");
  return passed;
}

bool test_one_third_film_rule() {
  const OneThirdFilmSample sample =
      sample_one_third_film(film_input(0.2, 0.02));
  bool passed = true;
  passed &= expect(sample.succeeded() &&
                       sample.model_id == "one_third_film_v1" &&
                       near(sample.surface_weight, 2.0 / 3.0, 0.0) &&
                       near(sample.far_gas_weight, 1.0 / 3.0, 0.0) &&
                       near(sample.reference_temperature_k, 550.0) &&
                       near(sample.reference_vapor_mass_fraction, 0.14) &&
                       sample.thermodynamic_pressure_pa == 101325.0,
                   "one-third film uses surface 2/3 and far gas 1/3");
  OneThirdFilmInput invalid = film_input();
  invalid.reference_cp_j_per_kg_k =
      std::numeric_limits<double>::quiet_NaN();
  const OneThirdFilmSample failed = sample_one_third_film(invalid);
  passed &= expect(!failed.succeeded() && failed.surface_weight == 0.0 &&
                       !std::signbit(failed.surface_weight) &&
                       failed.gas_velocity_m_per_s == Vector3{},
                   "film service failure returns no partial sample");
  return passed;
}

AbramzonSirignanoInput as_input(double surface_y, double far_y,
                                double duration_s) {
  AbramzonSirignanoInput input;
  input.parcel = parcel();
  input.film = sample_one_third_film(film_input(surface_y, far_y));
  input.liquid_properties =
      {750.0, 2200.0, 2.5e5, 2.0e4, 0.025, 8.0e-4};
  input.vapor_absolute_thermochemical_enthalpy_j_per_kg = 1.1e6;
  input.duration_s = duration_s;
  return input;
}

bool mass_momentum_consistent(const TransferExchangeCandidate& exchange) {
  bool closed = exchange.available &&
                exchange.parcel_liquid_mass_delta_kg +
                        exchange.gas_mass_delta_kg ==
                    0.0;
  for (std::size_t i = 0; i < 3U; ++i) {
    closed &= near(exchange.parcel_momentum_delta_kg_m_per_s[i] +
                       exchange.gas_momentum_delta_kg_m_per_s[i],
                   0.0, 0.0, 1.0e-18);
  }
  return closed;
}

bool test_as_limits_and_terminal_event() {
  AbramzonSirignanoInput zero_re = as_input(0.2, 0.02, 1.0e-5);
  zero_re.film.gas_velocity_m_per_s = zero_re.parcel.velocity_m_per_s;
  const AbramzonSirignanoReport evaporating =
      evaluate_abramzon_sirignano(zero_re);
  AbramzonSirignanoInput zero_duration_input = zero_re;
  zero_duration_input.duration_s = 0.0;
  const AbramzonSirignanoReport zero_duration =
      evaluate_abramzon_sirignano(zero_duration_input);
  const AbramzonSirignanoReport saturated =
      evaluate_abramzon_sirignano(as_input(0.02, 0.02, 1.0e-5));
  const AbramzonSirignanoInput finite_re_input =
      as_input(0.4, 0.02, 1.0e-5);
  const AbramzonSirignanoReport finite_re =
      evaluate_abramzon_sirignano(finite_re_input);
  const AbramzonSirignanoReport terminal =
      evaluate_abramzon_sirignano(as_input(0.8, 0.0, 100.0));

  bool passed = true;
  passed &= expect(evaporating.succeeded() && evaporating.bt_converged &&
                       evaporating.reynolds_number == 0.0 &&
                       evaporating.base_nusselt_number == 2.0 &&
                       evaporating.base_sherwood_number == 2.0 &&
                       evaporating.modified_nusselt_number == 2.0 &&
                       evaporating.modified_sherwood_number == 2.0 &&
                       evaporating.evaporation_mass_rate_one_droplet_kg_per_s >
                           0.0 &&
                       evaporating.diameter_squared_rate_m2_per_s < 0.0 &&
                       near(std::pow(
                                evaporating.candidate_parcel
                                    .droplet_diameter_m,
                                2.0),
                            std::pow(zero_re.parcel.droplet_diameter_m, 2.0) +
                                evaporating.diameter_squared_rate_m2_per_s *
                                    evaporating.advanced_duration_s) &&
                       mass_momentum_consistent(evaporating.exchange) &&
                       near(evaporating.exchange.thermal_exchange_to_gas_j,
                            evaporating.exchange
                                    .vapor_absolute_thermochemical_enthalpy_to_gas_j -
                                evaporating.exchange
                                    .convective_heat_to_parcel_j) &&
                       near(evaporating.exchange
                                .thermal_exchange_state_residual_j,
                            evaporating.exchange
                                    .parcel_thermochemical_enthalpy_delta_j +
                                evaporating.exchange
                                    .thermal_exchange_to_gas_j),
                   "Re=0 and d-squared A-S limits preserve explicit transfer budgets");
  passed &= expect(zero_duration.succeeded() &&
                       zero_duration.candidate_parcel.droplet_mass_kg ==
                           zero_duration_input.parcel.droplet_mass_kg &&
                       zero_duration.candidate_parcel.droplet_diameter_m ==
                           zero_duration_input.parcel.droplet_diameter_m &&
                       zero_duration.exchange.gas_mass_delta_kg == 0.0 &&
                       zero_duration.exchange.thermal_exchange_to_gas_j ==
                           0.0 &&
                       zero_duration.exchange.parcel_kinetic_energy_delta_j ==
                           0.0,
                   "zero duration degenerates exactly without roundoff source");
  passed &= expect(saturated.succeeded() && saturated.saturated &&
                       saturated.heating_only &&
                       saturated.evaporation_mass_rate_one_droplet_kg_per_s ==
                           0.0 &&
                       saturated.exchange.gas_mass_delta_kg == 0.0 &&
                       saturated.convective_heat_rate_one_droplet_w > 0.0,
                   "saturated gas is an explicit heating-only limit");
  const double finite_phi =
      finite_re_input.film.cp_j_per_kg_k *
      finite_re_input.film.density_kg_per_m3 *
      finite_re_input.film.vapor_diffusivity_m2_per_s *
      finite_re.modified_sherwood_number /
      (finite_re_input.film.thermal_conductivity_w_per_m_k *
       finite_re.modified_nusselt_number);
  const double fixed_point_bt = std::expm1(
      finite_phi * std::log1p(finite_re.spalding_mass_number));
  passed &= expect(
      finite_re.succeeded() && finite_re.bt_iterations > 0U &&
          near(finite_re.spalding_heat_number, fixed_point_bt, 2.0e-10) &&
          finite_re.base_nusselt_number > 2.0 &&
          finite_re.modified_nusselt_number >= 2.0 &&
          finite_re.modified_nusselt_number <
              finite_re.base_nusselt_number &&
          finite_re.modified_sherwood_number >= 2.0 &&
          finite_re.modified_sherwood_number <
              finite_re.base_sherwood_number,
      "finite-Re Stefan correction converges to the coupled B_T fixed point");
  passed &= expect(terminal.succeeded() && terminal.has_terminal_event &&
                       terminal.complete_evaporation &&
                       terminal.event_time_s > 0.0 &&
                       terminal.event_time_s < terminal.advanced_duration_s +
                                                   1.0e-15 &&
                       terminal.candidate_parcel.droplet_mass_kg == 0.0 &&
                       terminal.candidate_parcel.droplet_diameter_m == 0.0 &&
                       mass_momentum_consistent(terminal.exchange) &&
                       terminal.exchange.thermal_exchange_to_gas_j ==
                           terminal.exchange
                                   .vapor_absolute_thermochemical_enthalpy_to_gas_j -
                               terminal.exchange.convective_heat_to_parcel_j &&
                       terminal.exchange.thermal_exchange_state_residual_j !=
                           0.0,
                   "terminal event lands exactly at zero without negative clipping");

  const double speed_squared = dot(zero_re.parcel.velocity_m_per_s,
                                   zero_re.parcel.velocity_m_per_s);
  const double expected_kinetic_delta =
      0.5 * evaporating.exchange.parcel_liquid_mass_delta_kg *
      speed_squared;
  passed &= expect(
      near(evaporating.exchange.parcel_kinetic_energy_delta_j,
           expected_kinetic_delta) &&
          evaporating.exchange.parcel_kinetic_energy_delta_j < 0.0,
      "actual parcel kinetic-energy change is explicit and not hidden in h_tc");

  const AbramzonSirignanoReport condensation =
      evaluate_abramzon_sirignano(as_input(0.01, 0.02, 1.0e-5));
  passed &= expect(condensation.status ==
                       AbramzonSirignanoStatus::condensation_unsupported &&
                       !condensation.exchange.available &&
                       condensation.exchange.thermal_exchange_to_gas_j ==
                           0.0 &&
                       condensation.exchange.parcel_kinetic_energy_delta_j ==
                           0.0 &&
                       canonical_parcel(condensation.candidate_parcel),
                   "supersaturated condensation is not silently converted to no-transfer");

  AbramzonSirignanoInput no_iteration =
      as_input(0.6, 0.0, 1.0e-5);
  no_iteration.maximum_bt_iterations = 1U;
  no_iteration.bt_relative_tolerance = 1.0e-30;
  const AbramzonSirignanoReport iteration_failure =
      evaluate_abramzon_sirignano(no_iteration);
  passed &= expect(iteration_failure.status ==
                       AbramzonSirignanoStatus::iteration_not_converged &&
                       !iteration_failure.exchange.available &&
                       iteration_failure.exchange
                               .thermal_exchange_to_gas_j ==
                           0.0 &&
                       canonical_parcel(iteration_failure.candidate_parcel),
                   "unconverged B_T iteration exposes no partial candidate");

  SingleComponentEvaporationInput oracle;
  oracle.parcel_velocity_m_per_s = zero_re.parcel.velocity_m_per_s;
  oracle.droplet_mass_kg = zero_re.parcel.droplet_mass_kg;
  oracle.droplet_diameter_m = zero_re.parcel.droplet_diameter_m;
  oracle.multiplicity = zero_re.parcel.multiplicity;
  oracle.droplet_temperature_k = zero_re.parcel.temperature_k;
  oracle.gas_temperature_k = zero_re.film.far_gas_temperature_k;
  oracle.gas_pressure_pa = zero_re.film.thermodynamic_pressure_pa;
  oracle.gas_density_kg_per_m3 = zero_re.film.density_kg_per_m3;
  oracle.gas_vapor_mass_fraction =
      zero_re.film.far_gas_vapor_mass_fraction;
  oracle.vapor_molecular_weight_kg_per_kmol = 100.0;
  oracle.carrier_molecular_weight_kg_per_kmol = 29.0;
  oracle.reynolds_number = 0.0;
  oracle.prandtl_number = evaporating.prandtl_number;
  oracle.schmidt_number = evaporating.schmidt_number;
  oracle.gas_thermal_conductivity_w_per_m_k =
      zero_re.film.thermal_conductivity_w_per_m_k;
  oracle.vapor_diffusivity_m2_per_s =
      zero_re.film.vapor_diffusivity_m2_per_s;
  oracle.duration_s = zero_re.duration_s;
  oracle.liquid_properties = zero_re.liquid_properties;
  const SingleComponentEvaporationReport simplified_oracle =
      evaluate_single_component_evaporation(oracle);
  passed &= expect(simplified_oracle.succeeded() &&
                       simplified_oracle.diameter_squared_rate_m2_per_s < 0.0 &&
                       evaporating.diameter_squared_rate_m2_per_s < 0.0,
                   "legacy simplified kernel is consulted only as a monotone d2 oracle");
  return passed;
}

GasTransferEnvironment environment() {
  GasTransferEnvironment gas;
  gas.gas_velocity_m_per_s = {4.0, -0.5, 0.25};
  gas.gas_temperature_k = 950.0;
  gas.gas_vapor_mass_fraction = 0.01;
  gas.thermodynamic_pressure_pa = 101325.0;
  gas.film_density_kg_per_m3 = 0.8;
  gas.film_dynamic_viscosity_pa_s = 2.0e-5;
  gas.film_thermal_conductivity_w_per_m_k = 0.05;
  gas.film_vapor_diffusivity_m2_per_s = 2.0e-5;
  gas.film_cp_j_per_kg_k = 1100.0;
  gas.vapor_molecular_weight_kg_per_kmol = 100.0;
  gas.carrier_molecular_weight_kg_per_kmol = 29.0;
  gas.vapor_absolute_thermochemical_enthalpy_j_per_kg = 1.1e6;
  return gas;
}

FixedExchangeReport integrated(const LiquidPropertyService& service,
                               std::uint32_t substeps,
                               double multiplicity = 4.0) {
  FixedExchangeInput input;
  input.committed_parcel = parcel(multiplicity);
  input.liquid_properties = &service;
  input.predictor_environment = environment();
  input.corrector_environment = environment();
  input.corrector_environment.gas_temperature_k = 930.0;
  input.duration_s = 2.0e-5;
  input.internal_substeps = substeps;
  input.maximum_internal_substeps = 64U;
  return integrate_fixed_exchange(input);
}

bool test_fixed_two_pass_exchange_and_refinement() {
  const LiquidPropertyPack pack = base_pack(91U);
  const LiquidPropertyService service(&pack, 1U);
  const FixedExchangeReport coarse = integrated(service, 4U);
  const FixedExchangeReport fine = integrated(service, 16U);
  const FixedExchangeReport one = integrated(service, 16U, 1.0);
  const FixedExchangeReport five = integrated(service, 16U, 5.0);

  bool passed = true;
  passed &= expect(coarse.succeeded() && fine.succeeded() &&
                       coarse.predictor_passes == 1U &&
                       coarse.corrector_passes == 1U &&
                       coarse.predictor_evaluations == 4U &&
                       coarse.corrector_evaluations == 4U &&
                       fine.predictor_evaluations == 16U &&
                       fine.corrector_evaluations == 16U &&
                       mass_momentum_consistent(coarse.exchange) &&
                       mass_momentum_consistent(fine.exchange) &&
                       coarse.exchange.thermal_exchange_to_gas_j ==
                           coarse.exchange
                                   .vapor_absolute_thermochemical_enthalpy_to_gas_j -
                               coarse.exchange.convective_heat_to_parcel_j &&
                       near(coarse.exchange
                                .thermal_exchange_state_residual_j,
                            coarse.exchange
                                    .parcel_thermochemical_enthalpy_delta_j +
                                coarse.exchange
                                    .thermal_exchange_to_gas_j),
                   "bounded integration has exactly one predictor and one corrector pass");
  passed &= expect(
      near(coarse.candidate_parcel.droplet_mass_kg,
           fine.candidate_parcel.droplet_mass_kg, 2.0e-3) &&
          near(coarse.candidate_parcel.temperature_k,
               fine.candidate_parcel.temperature_k, 2.0e-3),
      "internal-substep refinement approaches the same finite candidate");
  passed &= expect(one.succeeded() && five.succeeded() &&
                       near(five.exchange.vapor_mass_to_gas_kg,
                            5.0 * one.exchange.vapor_mass_to_gas_kg) &&
                       near(five.exchange.vapor_momentum_to_gas_kg_m_per_s[0U],
                            5.0 * one.exchange.vapor_momentum_to_gas_kg_m_per_s[0U]) &&
                       near(five.exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j,
                            5.0 * one.exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j),
                   "multiplicity scales vapor mass, momentum and absolute enthalpy once");
  const double initial_kinetic =
      0.5 * parcel().droplet_mass_kg *
      dot(parcel().velocity_m_per_s, parcel().velocity_m_per_s) * 4.0;
  const double final_kinetic =
      0.5 * coarse.candidate_parcel.droplet_mass_kg *
      dot(coarse.candidate_parcel.velocity_m_per_s,
          coarse.candidate_parcel.velocity_m_per_s) *
      4.0;
  passed &= expect(
      near(coarse.exchange.parcel_kinetic_energy_delta_j,
           final_kinetic - initial_kinetic) &&
          std::isfinite(coarse.exchange.drag_work_to_parcel_j),
      "integrated report separates actual parcel delta-K from drag work");

  FixedExchangeInput limit;
  limit.committed_parcel = parcel();
  limit.liquid_properties = &service;
  limit.predictor_environment = environment();
  limit.corrector_environment = environment();
  limit.duration_s = 1.0e-5;
  limit.internal_substeps = 5U;
  limit.maximum_internal_substeps = 4U;
  const FixedExchangeReport rejected = integrate_fixed_exchange(limit);
  passed &= expect(rejected.status ==
                       FixedExchangeStatus::substep_limit_exceeded &&
                       !rejected.exchange.available &&
                       rejected.exchange.thermal_exchange_to_gas_j == 0.0 &&
                       canonical_parcel(rejected.candidate_parcel),
                   "substep overflow and failure publish no candidate");

  FixedExchangeInput missing = limit;
  missing.internal_substeps = 1U;
  missing.maximum_internal_substeps = 1U;
  missing.committed_parcel.liquid_material_fingerprint = 999U;
  const FixedExchangeReport missing_report = integrate_fixed_exchange(missing);
  passed &= expect(missing_report.status ==
                       FixedExchangeStatus::liquid_property_failure &&
                       !missing_report.exchange.available &&
                       missing_report.exchange.thermal_exchange_to_gas_j ==
                           0.0 &&
                       canonical_parcel(missing_report.candidate_parcel),
                   "property failure cannot expose partial parcel or gas state");
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  auto stale_film = as_input(0.2, 0.02, 1.0e-5);
  stale_film.film.surface_temperature_k += 5.0;
  const auto stale = evaluate_abramzon_sirignano(stale_film);
  passed &= expect(!stale.succeeded() && !stale.exchange.available,
                   "film sample must belong to the exact parcel surface temperature");
  passed &= test_liquid_property_dispatch_and_correlations();
  passed &= test_one_third_film_rule();
  passed &= test_as_limits_and_terminal_event();
  passed &= test_fixed_two_pass_exchange_and_refinement();
  if (!passed) return 1;
  std::cout << "models_spray_transfer_test: PASS\n";
  return 0;
}
