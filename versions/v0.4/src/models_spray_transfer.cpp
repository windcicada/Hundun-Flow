// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_transfer_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::spray::detail {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kUniversalGasConstantJPerKmolK = 8314.46261815324;
constexpr double kGeometryRelativeTolerance = 1.0e-10;

bool finite_vector(const Vector3& value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
}

double dot(const Vector3& left, const Vector3& right) noexcept {
  return left[0U] * right[0U] + left[1U] * right[1U] +
         left[2U] * right[2U];
}

LiquidPropertyReport liquid_failure(LiquidPropertyStatus status) noexcept {
  LiquidPropertyReport report;
  report.status = status;
  return report;
}

bool evaluate_temperature_correlation(const TemperatureCorrelation& law,
                                      double temperature_k,
                                      double& value) noexcept {
  value = 0.0;
  if (!std::isfinite(temperature_k) || !(temperature_k > 0.0)) return false;
  switch (law.kind) {
    case TemperatureCorrelationKind::constant:
      if (!std::isfinite(law.c[0U])) return false;
      value = law.c[0U];
      return true;
    case TemperatureCorrelationKind::polynomial_cubic: {
      if (!std::isfinite(law.reference_temperature_k) ||
          !(law.reference_temperature_k > 0.0)) {
        return false;
      }
      for (double coefficient : law.c) {
        if (!std::isfinite(coefficient)) return false;
      }
      const double theta = temperature_k - law.reference_temperature_k;
      value = ((law.c[3U] * theta + law.c[2U]) * theta + law.c[1U]) *
                  theta +
              law.c[0U];
      return true;
    }
  }
  return false;
}

bool evaluate_saturation_pressure(const SaturationPressureCorrelation& law,
                                  double temperature_k,
                                  double& pressure_pa) noexcept {
  pressure_pa = 0.0;
  if (!std::isfinite(temperature_k) || !(temperature_k > 0.0)) return false;
  switch (law.kind) {
    case SaturationPressureCorrelationKind::antoine_kelvin: {
      if (!std::isfinite(law.antoine_a) ||
          !std::isfinite(law.antoine_b_k) ||
          !std::isfinite(law.antoine_c_k) ||
          !std::isfinite(law.pressure_scale_pa) ||
          !(law.pressure_scale_pa > 0.0)) {
        return false;
      }
      const double denominator = temperature_k + law.antoine_c_k;
      if (!std::isfinite(denominator) || denominator == 0.0) return false;
      const double exponent =
          law.antoine_a - law.antoine_b_k / denominator;
      pressure_pa = law.pressure_scale_pa * std::pow(10.0, exponent);
      return true;
    }
    case SaturationPressureCorrelationKind::clausius_clapeyron: {
      if (!std::isfinite(law.reference_pressure_pa) ||
          !(law.reference_pressure_pa > 0.0) ||
          !std::isfinite(law.reference_temperature_k) ||
          !(law.reference_temperature_k > 0.0) ||
          !std::isfinite(law.latent_heat_j_per_kg) ||
          !(law.latent_heat_j_per_kg > 0.0) ||
          !std::isfinite(law.molecular_weight_kg_per_kmol) ||
          !(law.molecular_weight_kg_per_kmol > 0.0)) {
        return false;
      }
      const double exponent =
          law.latent_heat_j_per_kg * law.molecular_weight_kg_per_kmol /
          kUniversalGasConstantJPerKmolK *
          (1.0 / law.reference_temperature_k - 1.0 / temperature_k);
      pressure_pa = law.reference_pressure_pa * std::exp(exponent);
      return true;
    }
  }
  return false;
}

OneThirdFilmSample film_failure(FilmSampleStatus status) noexcept {
  OneThirdFilmSample report;
  report.status = status;
  return report;
}

AbramzonSirignanoReport as_failure(
    AbramzonSirignanoStatus status) noexcept {
  AbramzonSirignanoReport report;
  report.status = status;
  return report;
}

FixedExchangeReport fixed_failure(FixedExchangeStatus status) noexcept {
  FixedExchangeReport report;
  report.status = status;
  return report;
}

bool valid_liquid_properties(const LiquidProperties& liquid) noexcept {
  return std::isfinite(liquid.density_kg_per_m3) &&
         liquid.density_kg_per_m3 > 0.0 &&
         std::isfinite(liquid.cp_j_per_kg_k) &&
         liquid.cp_j_per_kg_k > 0.0 &&
         std::isfinite(liquid.latent_heat_j_per_kg) &&
         liquid.latent_heat_j_per_kg >= 0.0 &&
         std::isfinite(liquid.saturation_pressure_pa) &&
         liquid.saturation_pressure_pa >= 0.0 &&
         std::isfinite(liquid.surface_tension_n_per_m) &&
         liquid.surface_tension_n_per_m >= 0.0 &&
         std::isfinite(liquid.viscosity_pa_s) &&
         liquid.viscosity_pa_s >= 0.0;
}

bool geometry_consistent(const SprayParcelState& parcel,
                         double liquid_density_kg_per_m3) noexcept {
  const double geometric_mass = liquid_density_kg_per_m3 * kPi / 6.0 *
                                parcel.droplet_diameter_m *
                                parcel.droplet_diameter_m *
                                parcel.droplet_diameter_m;
  return std::isfinite(geometric_mass) && geometric_mass > 0.0 &&
         std::abs(geometric_mass - parcel.droplet_mass_kg) <=
             kGeometryRelativeTolerance *
                 std::max(geometric_mass, parcel.droplet_mass_kg);
}

double stefan_factor(double spalding_number) noexcept {
  if (spalding_number == 0.0) return 1.0;
  return std::pow(1.0 + spalding_number, 0.7) *
         std::log1p(spalding_number) / spalding_number;
}

bool finite_exchange(const TransferExchangeCandidate& exchange) noexcept {
  return std::isfinite(exchange.parcel_liquid_mass_delta_kg) &&
         finite_vector(exchange.parcel_momentum_delta_kg_m_per_s) &&
         std::isfinite(exchange.parcel_thermochemical_enthalpy_delta_j) &&
         std::isfinite(exchange.parcel_kinetic_energy_delta_j) &&
         std::isfinite(exchange.vapor_mass_to_gas_kg) &&
         finite_vector(exchange.vapor_momentum_to_gas_kg_m_per_s) &&
         std::isfinite(
             exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j) &&
         finite_vector(exchange.drag_momentum_to_parcel_kg_m_per_s) &&
         std::isfinite(exchange.convective_heat_to_parcel_j) &&
         std::isfinite(exchange.drag_work_to_parcel_j) &&
         std::isfinite(exchange.gas_mass_delta_kg) &&
         finite_vector(exchange.gas_momentum_delta_kg_m_per_s) &&
         std::isfinite(exchange.mass_closure_residual_kg) &&
         finite_vector(exchange.momentum_closure_residual_kg_m_per_s) &&
         std::isfinite(exchange.thermal_exchange_to_gas_j) &&
         std::isfinite(exchange.thermal_exchange_state_residual_j) &&
         finite_vector(exchange.momentum_quadrature_residual_kg_m_per_s);
}

TransferExchangeCandidate make_evaporation_exchange(
    const SprayParcelState& initial, const SprayParcelState& final,
    const LiquidProperties& liquid, double vapor_absolute_h_j_per_kg,
    double convective_heat_one_droplet_j) noexcept {
  TransferExchangeCandidate exchange;
  const double multiplicity = initial.multiplicity;
  exchange.parcel_liquid_mass_delta_kg =
      (final.droplet_mass_kg - initial.droplet_mass_kg) * multiplicity;
  exchange.vapor_mass_to_gas_kg =
      -exchange.parcel_liquid_mass_delta_kg;
  for (std::size_t component = 0U; component < 3U; ++component) {
    exchange.parcel_momentum_delta_kg_m_per_s[component] =
        exchange.parcel_liquid_mass_delta_kg *
        initial.velocity_m_per_s[component];
    exchange.vapor_momentum_to_gas_kg_m_per_s[component] =
        exchange.vapor_mass_to_gas_kg * initial.velocity_m_per_s[component];
  }

  const double initial_liquid_h =
      vapor_absolute_h_j_per_kg - liquid.latent_heat_j_per_kg;
  const double initial_h =
      initial.droplet_mass_kg * initial_liquid_h * multiplicity;
  double final_h = 0.0;
  if (final.droplet_mass_kg > 0.0) {
    const double final_liquid_h =
        initial_liquid_h + liquid.cp_j_per_kg_k *
                               (final.temperature_k - initial.temperature_k);
    final_h = final.droplet_mass_kg * final_liquid_h * multiplicity;
  }
  exchange.parcel_thermochemical_enthalpy_delta_j = final_h - initial_h;
  exchange.parcel_kinetic_energy_delta_j =
      0.5 * multiplicity *
      (final.droplet_mass_kg * dot(final.velocity_m_per_s,
                                   final.velocity_m_per_s) -
       initial.droplet_mass_kg * dot(initial.velocity_m_per_s,
                                     initial.velocity_m_per_s));
  exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j =
      exchange.vapor_mass_to_gas_kg * vapor_absolute_h_j_per_kg;
  exchange.convective_heat_to_parcel_j =
      convective_heat_one_droplet_j * multiplicity;
  exchange.thermal_exchange_to_gas_j =
      exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j -
      exchange.convective_heat_to_parcel_j;

  exchange.gas_mass_delta_kg =
      -exchange.parcel_liquid_mass_delta_kg;
  for (std::size_t component = 0U; component < 3U; ++component) {
    exchange.gas_momentum_delta_kg_m_per_s[component] =
        exchange.vapor_momentum_to_gas_kg_m_per_s[component] -
        exchange.drag_momentum_to_parcel_kg_m_per_s[component];
    exchange.momentum_closure_residual_kg_m_per_s[component] =
        exchange.parcel_momentum_delta_kg_m_per_s[component] +
        exchange.gas_momentum_delta_kg_m_per_s[component];
    exchange.momentum_quadrature_residual_kg_m_per_s[component] =
        -exchange.parcel_momentum_delta_kg_m_per_s[component] -
        exchange.gas_momentum_delta_kg_m_per_s[component];
  }
  exchange.mass_closure_residual_kg =
      exchange.parcel_liquid_mass_delta_kg + exchange.gas_mass_delta_kg;
  exchange.thermal_exchange_state_residual_j =
      exchange.parcel_thermochemical_enthalpy_delta_j +
      exchange.thermal_exchange_to_gas_j;
  exchange.available = finite_exchange(exchange);
  if (!exchange.available) return {};
  return exchange;
}

bool valid_environment(const GasTransferEnvironment& gas) noexcept {
  return finite_vector(gas.gas_velocity_m_per_s) &&
         std::isfinite(gas.gas_temperature_k) &&
         gas.gas_temperature_k > 0.0 &&
         std::isfinite(gas.gas_vapor_mass_fraction) &&
         gas.gas_vapor_mass_fraction >= 0.0 &&
         gas.gas_vapor_mass_fraction < 1.0 &&
         std::isfinite(gas.thermodynamic_pressure_pa) &&
         gas.thermodynamic_pressure_pa > 0.0 &&
         std::isfinite(gas.film_density_kg_per_m3) &&
         gas.film_density_kg_per_m3 > 0.0 &&
         std::isfinite(gas.film_dynamic_viscosity_pa_s) &&
         gas.film_dynamic_viscosity_pa_s > 0.0 &&
         std::isfinite(gas.film_thermal_conductivity_w_per_m_k) &&
         gas.film_thermal_conductivity_w_per_m_k > 0.0 &&
         std::isfinite(gas.film_vapor_diffusivity_m2_per_s) &&
         gas.film_vapor_diffusivity_m2_per_s > 0.0 &&
         std::isfinite(gas.film_cp_j_per_kg_k) &&
         gas.film_cp_j_per_kg_k > 0.0 &&
         std::isfinite(gas.vapor_molecular_weight_kg_per_kmol) &&
         gas.vapor_molecular_weight_kg_per_kmol > 0.0 &&
         std::isfinite(gas.carrier_molecular_weight_kg_per_kmol) &&
         gas.carrier_molecular_weight_kg_per_kmol > 0.0 &&
         std::isfinite(
             gas.vapor_absolute_thermochemical_enthalpy_j_per_kg);
}

bool surface_vapor_fraction(double saturation_pressure_pa,
                            const GasTransferEnvironment& gas,
                            double& mass_fraction) noexcept {
  mass_fraction = 0.0;
  if (!std::isfinite(saturation_pressure_pa) ||
      saturation_pressure_pa < 0.0 ||
      !(saturation_pressure_pa < gas.thermodynamic_pressure_pa)) {
    return false;
  }
  const double mole_fraction =
      saturation_pressure_pa / gas.thermodynamic_pressure_pa;
  const double denominator =
      mole_fraction * gas.vapor_molecular_weight_kg_per_kmol +
      (1.0 - mole_fraction) * gas.carrier_molecular_weight_kg_per_kmol;
  mass_fraction =
      mole_fraction * gas.vapor_molecular_weight_kg_per_kmol / denominator;
  return std::isfinite(mass_fraction) && mass_fraction >= 0.0 &&
         mass_fraction < 1.0;
}

struct LocalRates {
  AbramzonSirignanoReport evaporation{};
  SchillerNaumannDragReport drag{};
  LiquidProperties liquid{};
  double vapor_absolute_h_j_per_kg{};
};

FixedExchangeStatus evaluate_local_rates(
    const SprayParcelState& state, const LiquidPropertyService& service,
    const GasTransferEnvironment& gas, std::uint32_t maximum_bt_iterations,
    double bt_relative_tolerance, LocalRates& rates) noexcept {
  rates = {};
  const LiquidPropertyReport properties = service.evaluate(
      {state.liquid_material_fingerprint, state.temperature_k});
  if (!properties.succeeded() ||
      validate_liquid_property_report(
          {state.liquid_material_fingerprint, state.temperature_k},
          properties) != LiquidPropertyStatus::success) {
    return FixedExchangeStatus::liquid_property_failure;
  }

  double surface_y = 0.0;
  if (!surface_vapor_fraction(properties.properties.saturation_pressure_pa,
                              gas, surface_y)) {
    return FixedExchangeStatus::transfer_failure;
  }
  OneThirdFilmInput film_input;
  film_input.far_gas_velocity_m_per_s = gas.gas_velocity_m_per_s;
  film_input.droplet_surface_temperature_k = state.temperature_k;
  film_input.far_gas_temperature_k = gas.gas_temperature_k;
  film_input.surface_vapor_mass_fraction = surface_y;
  film_input.far_gas_vapor_mass_fraction = gas.gas_vapor_mass_fraction;
  film_input.thermodynamic_pressure_pa = gas.thermodynamic_pressure_pa;
  film_input.reference_density_kg_per_m3 = gas.film_density_kg_per_m3;
  film_input.reference_dynamic_viscosity_pa_s =
      gas.film_dynamic_viscosity_pa_s;
  film_input.reference_thermal_conductivity_w_per_m_k =
      gas.film_thermal_conductivity_w_per_m_k;
  film_input.reference_vapor_diffusivity_m2_per_s =
      gas.film_vapor_diffusivity_m2_per_s;
  film_input.reference_cp_j_per_kg_k = gas.film_cp_j_per_kg_k;
  const OneThirdFilmSample film = sample_one_third_film(film_input);
  if (!film.succeeded()) return FixedExchangeStatus::film_sample_failure;

  AbramzonSirignanoInput evaporation_input;
  evaporation_input.parcel = state;
  evaporation_input.film = film;
  evaporation_input.liquid_properties = properties.properties;
  evaporation_input.vapor_absolute_thermochemical_enthalpy_j_per_kg =
      gas.vapor_absolute_thermochemical_enthalpy_j_per_kg;
  evaporation_input.duration_s = 0.0;
  evaporation_input.maximum_bt_iterations = maximum_bt_iterations;
  evaporation_input.bt_relative_tolerance = bt_relative_tolerance;
  rates.evaporation = evaluate_abramzon_sirignano(evaporation_input);
  if (!rates.evaporation.succeeded()) {
    return FixedExchangeStatus::transfer_failure;
  }

  SchillerNaumannDragInput drag_input;
  drag_input.gas_velocity_m_per_s = gas.gas_velocity_m_per_s;
  drag_input.parcel_velocity_m_per_s = state.velocity_m_per_s;
  drag_input.gas_density_kg_per_m3 = film.density_kg_per_m3;
  drag_input.gas_dynamic_viscosity_pa_s =
      film.dynamic_viscosity_pa_s;
  drag_input.droplet_diameter_m = state.droplet_diameter_m;
  drag_input.droplet_mass_kg = state.droplet_mass_kg;
  drag_input.multiplicity = state.multiplicity;
  drag_input.duration_s = 0.0;
  rates.drag = evaluate_schiller_naumann_drag(drag_input);
  if (!rates.drag.succeeded()) return FixedExchangeStatus::transfer_failure;
  rates.liquid = properties.properties;
  rates.vapor_absolute_h_j_per_kg =
      gas.vapor_absolute_thermochemical_enthalpy_j_per_kg;
  return FixedExchangeStatus::success;
}

bool make_euler_state(const SprayParcelState& initial,
                      const LocalRates& rates, double duration_s,
                      SprayParcelState& candidate) noexcept {
  candidate = initial;
  candidate.droplet_mass_kg -=
      rates.evaporation.evaporation_mass_rate_one_droplet_kg_per_s *
      duration_s;
  if (!std::isfinite(candidate.droplet_mass_kg) ||
      !(candidate.droplet_mass_kg > 0.0)) {
    return false;
  }
  candidate.droplet_diameter_m = std::cbrt(
      6.0 * candidate.droplet_mass_kg /
      (kPi * rates.liquid.density_kg_per_m3));
  candidate.temperature_k +=
      rates.evaporation.temperature_rate_k_per_s * duration_s;
  for (std::size_t component = 0U; component < 3U; ++component) {
    candidate.velocity_m_per_s[component] +=
        rates.drag.acceleration_m_per_s2[component] * duration_s;
  }
  candidate.age_s += duration_s;
  return std::isfinite(candidate.droplet_mass_kg) &&
         candidate.droplet_mass_kg > 0.0 &&
         std::isfinite(candidate.droplet_diameter_m) &&
         candidate.droplet_diameter_m > 0.0 &&
         std::isfinite(candidate.temperature_k) &&
         candidate.temperature_k > 0.0 &&
         finite_vector(candidate.velocity_m_per_s) &&
         std::isfinite(candidate.age_s);
}

}  // namespace

LiquidPropertyReport LiquidPropertyService::evaluate(
    const LiquidPropertyQuery& query) const noexcept {
  if (query.material_fingerprint == 0U ||
      !std::isfinite(query.temperature_k) ||
      !(query.temperature_k > 0.0) || (count_ > 0U && packs_ == nullptr)) {
    return liquid_failure(LiquidPropertyStatus::invalid_input);
  }

  const LiquidPropertyPack* selected = nullptr;
  for (std::size_t index = 0U; index < count_; ++index) {
    if (packs_[index].material_fingerprint == query.material_fingerprint) {
      if (selected != nullptr) {
        return liquid_failure(LiquidPropertyStatus::provider_contract_failure);
      }
      selected = &packs_[index];
    }
  }
  if (selected == nullptr) {
    return liquid_failure(LiquidPropertyStatus::unknown_material);
  }
  if (!std::isfinite(selected->minimum_temperature_k) ||
      !std::isfinite(selected->maximum_temperature_k) ||
      !(selected->minimum_temperature_k > 0.0) ||
      selected->maximum_temperature_k < selected->minimum_temperature_k) {
    return liquid_failure(LiquidPropertyStatus::correlation_domain_error);
  }
  if (query.temperature_k < selected->minimum_temperature_k ||
      query.temperature_k > selected->maximum_temperature_k) {
    return liquid_failure(LiquidPropertyStatus::temperature_out_of_range);
  }

  LiquidProperties properties;
  bool evaluated = evaluate_temperature_correlation(
      selected->density_kg_per_m3, query.temperature_k,
      properties.density_kg_per_m3);
  evaluated &= evaluate_temperature_correlation(
      selected->cp_j_per_kg_k, query.temperature_k,
      properties.cp_j_per_kg_k);
  evaluated &= evaluate_temperature_correlation(
      selected->latent_heat_j_per_kg, query.temperature_k,
      properties.latent_heat_j_per_kg);
  evaluated &= evaluate_temperature_correlation(
      selected->surface_tension_n_per_m, query.temperature_k,
      properties.surface_tension_n_per_m);
  evaluated &= evaluate_temperature_correlation(
      selected->viscosity_pa_s, query.temperature_k,
      properties.viscosity_pa_s);
  evaluated &= evaluate_saturation_pressure(
      selected->saturation_pressure, query.temperature_k,
      properties.saturation_pressure_pa);
  if (!evaluated) {
    return liquid_failure(LiquidPropertyStatus::correlation_domain_error);
  }
  const double values[]{properties.density_kg_per_m3,
                        properties.cp_j_per_kg_k,
                        properties.latent_heat_j_per_kg,
                        properties.saturation_pressure_pa,
                        properties.surface_tension_n_per_m,
                        properties.viscosity_pa_s};
  for (double value : values) {
    if (!std::isfinite(value)) {
      return liquid_failure(LiquidPropertyStatus::non_finite_output);
    }
  }
  if (!valid_liquid_properties(properties)) {
    return liquid_failure(LiquidPropertyStatus::correlation_domain_error);
  }
  return {LiquidPropertyStatus::success, query.material_fingerprint,
          query.temperature_k, properties};
}

OneThirdFilmSample sample_one_third_film(
    const OneThirdFilmInput& input) noexcept {
  if (!finite_vector(input.far_gas_velocity_m_per_s) ||
      !std::isfinite(input.droplet_surface_temperature_k) ||
      !(input.droplet_surface_temperature_k > 0.0) ||
      !std::isfinite(input.far_gas_temperature_k) ||
      !(input.far_gas_temperature_k > 0.0) ||
      !std::isfinite(input.surface_vapor_mass_fraction) ||
      input.surface_vapor_mass_fraction < 0.0 ||
      !(input.surface_vapor_mass_fraction < 1.0) ||
      !std::isfinite(input.far_gas_vapor_mass_fraction) ||
      input.far_gas_vapor_mass_fraction < 0.0 ||
      !(input.far_gas_vapor_mass_fraction < 1.0) ||
      !std::isfinite(input.thermodynamic_pressure_pa) ||
      !(input.thermodynamic_pressure_pa > 0.0) ||
      !std::isfinite(input.reference_density_kg_per_m3) ||
      !(input.reference_density_kg_per_m3 > 0.0) ||
      !std::isfinite(input.reference_dynamic_viscosity_pa_s) ||
      !(input.reference_dynamic_viscosity_pa_s > 0.0) ||
      !std::isfinite(input.reference_thermal_conductivity_w_per_m_k) ||
      !(input.reference_thermal_conductivity_w_per_m_k > 0.0) ||
      !std::isfinite(input.reference_vapor_diffusivity_m2_per_s) ||
      !(input.reference_vapor_diffusivity_m2_per_s > 0.0) ||
      !std::isfinite(input.reference_cp_j_per_kg_k) ||
      !(input.reference_cp_j_per_kg_k > 0.0)) {
    return film_failure(FilmSampleStatus::invalid_input);
  }
  OneThirdFilmSample sample;
  sample.surface_weight = 2.0 / 3.0;
  sample.far_gas_weight = 1.0 / 3.0;
  sample.gas_velocity_m_per_s = input.far_gas_velocity_m_per_s;
  sample.surface_temperature_k = input.droplet_surface_temperature_k;
  sample.far_gas_temperature_k = input.far_gas_temperature_k;
  sample.reference_temperature_k =
      sample.surface_weight * input.droplet_surface_temperature_k +
      sample.far_gas_weight * input.far_gas_temperature_k;
  sample.surface_vapor_mass_fraction = input.surface_vapor_mass_fraction;
  sample.far_gas_vapor_mass_fraction = input.far_gas_vapor_mass_fraction;
  sample.reference_vapor_mass_fraction =
      sample.surface_weight * input.surface_vapor_mass_fraction +
      sample.far_gas_weight * input.far_gas_vapor_mass_fraction;
  sample.thermodynamic_pressure_pa = input.thermodynamic_pressure_pa;
  sample.density_kg_per_m3 = input.reference_density_kg_per_m3;
  sample.dynamic_viscosity_pa_s =
      input.reference_dynamic_viscosity_pa_s;
  sample.thermal_conductivity_w_per_m_k =
      input.reference_thermal_conductivity_w_per_m_k;
  sample.vapor_diffusivity_m2_per_s =
      input.reference_vapor_diffusivity_m2_per_s;
  sample.cp_j_per_kg_k = input.reference_cp_j_per_kg_k;
  if (!std::isfinite(sample.reference_temperature_k) ||
      !(sample.reference_temperature_k > 0.0) ||
      !std::isfinite(sample.reference_vapor_mass_fraction) ||
      sample.reference_vapor_mass_fraction < 0.0 ||
      !(sample.reference_vapor_mass_fraction < 1.0)) {
    return film_failure(FilmSampleStatus::non_finite_output);
  }
  sample.status = FilmSampleStatus::success;
  return sample;
}

AbramzonSirignanoReport evaluate_abramzon_sirignano(
    const AbramzonSirignanoInput& input) noexcept {
  if (validate_parcel_state(input.parcel) != ParcelStateStatus::success ||
      !input.film.succeeded() ||
      input.film.surface_temperature_k != input.parcel.temperature_k ||
      !valid_liquid_properties(input.liquid_properties) ||
      !std::isfinite(
          input.vapor_absolute_thermochemical_enthalpy_j_per_kg) ||
      !std::isfinite(input.duration_s) || input.duration_s < 0.0 ||
      input.maximum_bt_iterations == 0U ||
      !std::isfinite(input.bt_relative_tolerance) ||
      !(input.bt_relative_tolerance > 0.0)) {
    return as_failure(AbramzonSirignanoStatus::invalid_input);
  }
  if (!geometry_consistent(input.parcel,
                           input.liquid_properties.density_kg_per_m3)) {
    return as_failure(
        AbramzonSirignanoStatus::inconsistent_droplet_geometry);
  }
  if (!(input.film.surface_vapor_mass_fraction < 1.0) ||
      input.film.surface_vapor_mass_fraction < 0.0 ||
      !(input.film.far_gas_vapor_mass_fraction < 1.0) ||
      input.film.far_gas_vapor_mass_fraction < 0.0) {
    return as_failure(
        AbramzonSirignanoStatus::phase_equilibrium_unavailable);
  }
  if (input.film.surface_vapor_mass_fraction <
      input.film.far_gas_vapor_mass_fraction) {
    return as_failure(AbramzonSirignanoStatus::condensation_unsupported);
  }

  AbramzonSirignanoReport report;
  Vector3 slip{};
  for (std::size_t component = 0U; component < 3U; ++component) {
    slip[component] = input.film.gas_velocity_m_per_s[component] -
                      input.parcel.velocity_m_per_s[component];
  }
  const double slip_speed = std::hypot(slip[0U], slip[1U], slip[2U]);
  report.reynolds_number =
      input.film.density_kg_per_m3 * slip_speed *
      input.parcel.droplet_diameter_m /
      input.film.dynamic_viscosity_pa_s;
  report.prandtl_number =
      input.film.dynamic_viscosity_pa_s * input.film.cp_j_per_kg_k /
      input.film.thermal_conductivity_w_per_m_k;
  report.schmidt_number =
      input.film.dynamic_viscosity_pa_s /
      (input.film.density_kg_per_m3 *
       input.film.vapor_diffusivity_m2_per_s);
  if (!std::isfinite(report.reynolds_number) ||
      report.reynolds_number < 0.0 ||
      !std::isfinite(report.prandtl_number) ||
      !(report.prandtl_number > 0.0) ||
      !std::isfinite(report.schmidt_number) ||
      !(report.schmidt_number > 0.0)) {
    return as_failure(AbramzonSirignanoStatus::non_finite_output);
  }

  const double reynolds_factor =
      std::max(1.0,
               std::pow(std::min(400.0, report.reynolds_number), 0.077));
  report.base_nusselt_number =
      1.0 + std::cbrt(1.0 +
                       report.reynolds_number * report.prandtl_number) *
                reynolds_factor;
  report.base_sherwood_number =
      1.0 + std::cbrt(1.0 +
                       report.reynolds_number * report.schmidt_number) *
                reynolds_factor;
  report.spalding_mass_number =
      (input.film.surface_vapor_mass_fraction -
       input.film.far_gas_vapor_mass_fraction) /
      (1.0 - input.film.surface_vapor_mass_fraction);
  if (!std::isfinite(report.base_nusselt_number) ||
      !std::isfinite(report.base_sherwood_number) ||
      !std::isfinite(report.spalding_mass_number) ||
      report.spalding_mass_number < 0.0) {
    return as_failure(AbramzonSirignanoStatus::non_finite_output);
  }

  report.saturated = report.spalding_mass_number == 0.0;
  report.heating_only = report.saturated;
  if (report.saturated) {
    report.modified_nusselt_number = report.base_nusselt_number;
    report.modified_sherwood_number = report.base_sherwood_number;
    report.bt_converged = true;
  } else {
    const double mass_factor = stefan_factor(report.spalding_mass_number);
    report.modified_sherwood_number =
        2.0 + (report.base_sherwood_number - 2.0) / mass_factor;
    double heat_number = report.spalding_mass_number;
    bool converged = false;
    for (std::uint32_t iteration = 1U;
         iteration <= input.maximum_bt_iterations; ++iteration) {
      const double heat_factor = stefan_factor(heat_number);
      const double modified_nusselt =
          2.0 + (report.base_nusselt_number - 2.0) / heat_factor;
      const double phi =
          input.film.cp_j_per_kg_k *
          input.film.density_kg_per_m3 *
          input.film.vapor_diffusivity_m2_per_s *
          report.modified_sherwood_number /
          (input.film.thermal_conductivity_w_per_m_k *
           modified_nusselt);
      const double next_heat_number = std::expm1(
          phi * std::log1p(report.spalding_mass_number));
      if (!std::isfinite(heat_factor) || !(heat_factor > 0.0) ||
          !std::isfinite(modified_nusselt) ||
          !(modified_nusselt > 0.0) || !std::isfinite(phi) ||
          phi < 0.0 || !std::isfinite(next_heat_number) ||
          next_heat_number < 0.0) {
        return as_failure(AbramzonSirignanoStatus::non_finite_output);
      }
      report.bt_iterations = iteration;
      if (std::abs(next_heat_number - heat_number) <=
          input.bt_relative_tolerance *
              std::max(1.0, std::abs(next_heat_number))) {
        heat_number = next_heat_number;
        report.modified_nusselt_number =
            2.0 + (report.base_nusselt_number - 2.0) /
                      stefan_factor(heat_number);
        converged = true;
        break;
      }
      heat_number = next_heat_number;
    }
    if (!converged) {
      return as_failure(
          AbramzonSirignanoStatus::iteration_not_converged);
    }
    report.bt_converged = true;
    report.spalding_heat_number = heat_number;
  }

  const double heat_spalding_ratio =
      report.spalding_heat_number == 0.0
          ? 1.0
          : std::log1p(report.spalding_heat_number) /
                report.spalding_heat_number;
  report.evaporation_mass_rate_one_droplet_kg_per_s =
      kPi * input.film.density_kg_per_m3 *
      input.film.vapor_diffusivity_m2_per_s *
      input.parcel.droplet_diameter_m * report.modified_sherwood_number *
      std::log1p(report.spalding_mass_number);
  report.convective_heat_rate_one_droplet_w =
      kPi * input.film.thermal_conductivity_w_per_m_k *
      input.parcel.droplet_diameter_m *
      (input.film.far_gas_temperature_k -
       input.parcel.temperature_k) *
      report.modified_nusselt_number * heat_spalding_ratio;
  report.diameter_squared_rate_m2_per_s =
      -4.0 * report.evaporation_mass_rate_one_droplet_kg_per_s /
      (kPi * input.liquid_properties.density_kg_per_m3 *
       input.parcel.droplet_diameter_m);
  report.temperature_rate_k_per_s =
      (report.convective_heat_rate_one_droplet_w -
       report.evaporation_mass_rate_one_droplet_kg_per_s *
           input.liquid_properties.latent_heat_j_per_kg) /
      (input.parcel.droplet_mass_kg *
       input.liquid_properties.cp_j_per_kg_k);
  if (!std::isfinite(report.modified_nusselt_number) ||
      !(report.modified_nusselt_number > 0.0) ||
      !std::isfinite(report.modified_sherwood_number) ||
      !(report.modified_sherwood_number > 0.0) ||
      !std::isfinite(report.evaporation_mass_rate_one_droplet_kg_per_s) ||
      report.evaporation_mass_rate_one_droplet_kg_per_s < 0.0 ||
      !std::isfinite(report.convective_heat_rate_one_droplet_w) ||
      !std::isfinite(report.diameter_squared_rate_m2_per_s) ||
      report.diameter_squared_rate_m2_per_s > 0.0 ||
      !std::isfinite(report.temperature_rate_k_per_s)) {
    return as_failure(AbramzonSirignanoStatus::non_finite_output);
  }

  if (input.duration_s == 0.0) {
    report.candidate_parcel = input.parcel;
    report.exchange = make_evaporation_exchange(
        input.parcel, input.parcel, input.liquid_properties,
        input.vapor_absolute_thermochemical_enthalpy_j_per_kg, 0.0);
    if (!report.exchange.available) {
      return as_failure(AbramzonSirignanoStatus::non_finite_output);
    }
    report.status = AbramzonSirignanoStatus::success;
    return report;
  }

  report.advanced_duration_s = input.duration_s;
  const double initial_d2 = input.parcel.droplet_diameter_m *
                            input.parcel.droplet_diameter_m;
  double final_d2 =
      initial_d2 + report.diameter_squared_rate_m2_per_s * input.duration_s;
  if (report.diameter_squared_rate_m2_per_s < 0.0) {
    const double event_time =
        -initial_d2 / report.diameter_squared_rate_m2_per_s;
    if (event_time <= input.duration_s) {
      report.has_terminal_event = true;
      report.complete_evaporation = true;
      report.event_time_s = event_time;
      report.advanced_duration_s = event_time;
      final_d2 = 0.0;
    }
  }

  report.candidate_parcel = input.parcel;
  if (report.complete_evaporation) {
    report.candidate_parcel.droplet_mass_kg = 0.0;
    report.candidate_parcel.droplet_diameter_m = 0.0;
    report.candidate_parcel.temperature_k =
        input.parcel.temperature_k +
        report.temperature_rate_k_per_s * report.advanced_duration_s;
  } else {
    if (!std::isfinite(final_d2) || !(final_d2 > 0.0)) {
      return as_failure(AbramzonSirignanoStatus::non_finite_output);
    }
    if (report.evaporation_mass_rate_one_droplet_kg_per_s == 0.0) {
      // Preserve the saturated/no-mass-transfer limit bit for bit.  Rebuilding
      // the same sphere from d^2 would otherwise manufacture a roundoff-sized
      // vapor source.
      report.candidate_parcel.droplet_diameter_m =
          input.parcel.droplet_diameter_m;
      report.candidate_parcel.droplet_mass_kg =
          input.parcel.droplet_mass_kg;
    } else {
      report.candidate_parcel.droplet_diameter_m = std::sqrt(final_d2);
      report.candidate_parcel.droplet_mass_kg =
          input.liquid_properties.density_kg_per_m3 * kPi / 6.0 * final_d2 *
          report.candidate_parcel.droplet_diameter_m;
    }
    const double evaporated_mass =
        input.parcel.droplet_mass_kg -
        report.candidate_parcel.droplet_mass_kg;
    report.candidate_parcel.temperature_k =
        input.parcel.temperature_k +
        (report.convective_heat_rate_one_droplet_w *
             report.advanced_duration_s -
         evaporated_mass * input.liquid_properties.latent_heat_j_per_kg) /
            (report.candidate_parcel.droplet_mass_kg *
             input.liquid_properties.cp_j_per_kg_k);
  }
  report.candidate_parcel.age_s += report.advanced_duration_s;
  if (!std::isfinite(report.candidate_parcel.temperature_k) ||
      !(report.candidate_parcel.temperature_k > 0.0) ||
      !std::isfinite(report.candidate_parcel.age_s)) {
    return as_failure(AbramzonSirignanoStatus::non_finite_output);
  }
  report.exchange = make_evaporation_exchange(
      input.parcel, report.candidate_parcel, input.liquid_properties,
      input.vapor_absolute_thermochemical_enthalpy_j_per_kg,
      report.convective_heat_rate_one_droplet_w *
          report.advanced_duration_s);
  if (!report.exchange.available) {
    return as_failure(AbramzonSirignanoStatus::non_finite_output);
  }
  report.status = AbramzonSirignanoStatus::success;
  return report;
}

FixedExchangeReport integrate_fixed_exchange(
    const FixedExchangeInput& input) noexcept {
  if (validate_parcel_state(input.committed_parcel) !=
          ParcelStateStatus::success ||
      input.liquid_properties == nullptr ||
      !valid_environment(input.predictor_environment) ||
      !valid_environment(input.corrector_environment) ||
      !std::isfinite(input.duration_s) || input.duration_s < 0.0 ||
      input.internal_substeps == 0U ||
      input.maximum_internal_substeps == 0U ||
      input.maximum_bt_iterations == 0U ||
      !std::isfinite(input.bt_relative_tolerance) ||
      !(input.bt_relative_tolerance > 0.0)) {
    return fixed_failure(FixedExchangeStatus::invalid_input);
  }
  if (input.internal_substeps > input.maximum_internal_substeps) {
    return fixed_failure(FixedExchangeStatus::substep_limit_exceeded);
  }

  const LiquidPropertyReport initial_properties =
      input.liquid_properties->evaluate(
          {input.committed_parcel.liquid_material_fingerprint,
           input.committed_parcel.temperature_k});
  if (!initial_properties.succeeded()) {
    return fixed_failure(FixedExchangeStatus::liquid_property_failure);
  }

  if (input.duration_s == 0.0) {
    FixedExchangeReport report;
    report.candidate_parcel = input.committed_parcel;
    report.exchange.available = true;
    report.status = FixedExchangeStatus::success;
    return report;
  }

  FixedExchangeReport report;
  report.predictor_passes = 1U;
  report.corrector_passes = 1U;
  SprayParcelState current = input.committed_parcel;
  const double substep_duration =
      input.duration_s / static_cast<double>(input.internal_substeps);
  double vapor_absolute_h_integral = 0.0;
  Vector3 vapor_momentum_integral{};
  Vector3 drag_impulse_integral{};
  double heat_integral = 0.0;
  double drag_work_integral = 0.0;
  double elapsed = 0.0;

  for (std::uint32_t substep = 0U; substep < input.internal_substeps;
       ++substep) {
    LocalRates predictor;
    const FixedExchangeStatus predictor_status = evaluate_local_rates(
        current, *input.liquid_properties, input.predictor_environment,
        input.maximum_bt_iterations, input.bt_relative_tolerance, predictor);
    if (predictor_status != FixedExchangeStatus::success) {
      return fixed_failure(predictor_status);
    }
    ++report.predictor_evaluations;

    double advance = substep_duration;
    if (predictor.evaporation.evaporation_mass_rate_one_droplet_kg_per_s >
        0.0) {
      const double predictor_event =
          current.droplet_mass_kg /
          predictor.evaporation.evaporation_mass_rate_one_droplet_kg_per_s;
      if (predictor_event <= advance) {
        advance = 0.5 * predictor_event;
      }
    }
    SprayParcelState predicted;
    if (!make_euler_state(current, predictor, advance, predicted)) {
      return fixed_failure(FixedExchangeStatus::non_finite_output);
    }
    const LiquidPropertyReport predicted_properties =
        input.liquid_properties->evaluate(
            {predicted.liquid_material_fingerprint,
             predicted.temperature_k});
    if (!predicted_properties.succeeded()) {
      return fixed_failure(FixedExchangeStatus::liquid_property_failure);
    }
    predicted.droplet_diameter_m = std::cbrt(
        6.0 * predicted.droplet_mass_kg /
        (kPi * predicted_properties.properties.density_kg_per_m3));
    if (!std::isfinite(predicted.droplet_diameter_m) ||
        !(predicted.droplet_diameter_m > 0.0)) {
      return fixed_failure(FixedExchangeStatus::non_finite_output);
    }

    LocalRates corrector;
    const FixedExchangeStatus corrector_status = evaluate_local_rates(
        predicted, *input.liquid_properties, input.corrector_environment,
        input.maximum_bt_iterations, input.bt_relative_tolerance, corrector);
    if (corrector_status != FixedExchangeStatus::success) {
      return fixed_failure(corrector_status);
    }
    ++report.corrector_evaluations;

    const double average_mass_rate =
        0.5 *
        (predictor.evaporation.evaporation_mass_rate_one_droplet_kg_per_s +
         corrector.evaporation.evaporation_mass_rate_one_droplet_kg_per_s);
    const double average_temperature_rate =
        0.5 * (predictor.evaporation.temperature_rate_k_per_s +
               corrector.evaporation.temperature_rate_k_per_s);
    double actual_advance = substep_duration;
    bool terminal = false;
    if (average_mass_rate > 0.0) {
      const double event_time = current.droplet_mass_kg / average_mass_rate;
      if (event_time <= actual_advance) {
        actual_advance = event_time;
        terminal = true;
      }
    }

    SprayParcelState next = current;
    if (terminal) {
      next.droplet_diameter_m = 0.0;
      next.droplet_mass_kg = 0.0;
    } else {
      next.droplet_mass_kg =
          current.droplet_mass_kg - average_mass_rate * actual_advance;
      if (!std::isfinite(next.droplet_mass_kg) ||
          !(next.droplet_mass_kg > 0.0)) {
        return fixed_failure(FixedExchangeStatus::non_finite_output);
      }
    }
    next.temperature_k += average_temperature_rate * actual_advance;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double average_acceleration =
          0.5 * (predictor.drag.acceleration_m_per_s2[component] +
                 corrector.drag.acceleration_m_per_s2[component]);
      next.velocity_m_per_s[component] +=
          average_acceleration * actual_advance;
    }
    next.age_s += actual_advance;
    if (!std::isfinite(next.temperature_k) || !(next.temperature_k > 0.0) ||
        !finite_vector(next.velocity_m_per_s) ||
        !std::isfinite(next.age_s)) {
      return fixed_failure(FixedExchangeStatus::non_finite_output);
    }
    if (!terminal) {
      const LiquidPropertyReport next_properties =
          input.liquid_properties->evaluate(
              {next.liquid_material_fingerprint, next.temperature_k});
      if (!next_properties.succeeded()) {
        return fixed_failure(FixedExchangeStatus::liquid_property_failure);
      }
      next.droplet_diameter_m = std::cbrt(
          6.0 * next.droplet_mass_kg /
          (kPi * next_properties.properties.density_kg_per_m3));
      if (!std::isfinite(next.droplet_diameter_m) ||
          !(next.droplet_diameter_m > 0.0)) {
        return fixed_failure(FixedExchangeStatus::non_finite_output);
      }
    }

    const double evaporated_mass_from_state =
        current.droplet_mass_kg - next.droplet_mass_kg;
    const double resolved_mass = evaporated_mass_from_state;
    const double average_vapor_h =
        0.5 * (predictor.vapor_absolute_h_j_per_kg +
               corrector.vapor_absolute_h_j_per_kg);
    vapor_absolute_h_integral += resolved_mass * average_vapor_h;
    for (std::size_t component = 0U; component < 3U; ++component) {
      const double emission_velocity =
          0.5 * (current.velocity_m_per_s[component] +
                 predicted.velocity_m_per_s[component]);
      vapor_momentum_integral[component] +=
          resolved_mass * emission_velocity;
      const double average_force =
          0.5 * (predictor.drag.force_on_one_droplet_n[component] +
                 corrector.drag.force_on_one_droplet_n[component]);
      drag_impulse_integral[component] += average_force * actual_advance;
    }
    heat_integral +=
        0.5 *
        (predictor.evaporation.convective_heat_rate_one_droplet_w +
         corrector.evaporation.convective_heat_rate_one_droplet_w) *
        actual_advance;
    drag_work_integral +=
        0.5 *
        (dot(predictor.drag.force_on_one_droplet_n,
             current.velocity_m_per_s) +
         dot(corrector.drag.force_on_one_droplet_n,
             predicted.velocity_m_per_s)) *
        actual_advance;

    elapsed += actual_advance;
    current = next;
    ++report.completed_substeps;
    if (terminal) {
      report.has_terminal_event = true;
      report.event_time_s = elapsed;
      break;
    }
  }

  report.advanced_duration_s = elapsed;
  report.candidate_parcel = current;
  TransferExchangeCandidate exchange;
  const double multiplicity = input.committed_parcel.multiplicity;
  exchange.parcel_liquid_mass_delta_kg =
      (current.droplet_mass_kg -
       input.committed_parcel.droplet_mass_kg) *
      multiplicity;
  exchange.vapor_mass_to_gas_kg =
      -exchange.parcel_liquid_mass_delta_kg;
  for (std::size_t component = 0U; component < 3U; ++component) {
    exchange.parcel_momentum_delta_kg_m_per_s[component] =
        (current.droplet_mass_kg * current.velocity_m_per_s[component] -
         input.committed_parcel.droplet_mass_kg *
             input.committed_parcel.velocity_m_per_s[component]) *
        multiplicity;
    exchange.vapor_momentum_to_gas_kg_m_per_s[component] =
        vapor_momentum_integral[component] * multiplicity;
    exchange.drag_momentum_to_parcel_kg_m_per_s[component] =
        drag_impulse_integral[component] * multiplicity;
  }

  const double initial_liquid_h =
      input.predictor_environment
          .vapor_absolute_thermochemical_enthalpy_j_per_kg -
      initial_properties.properties.latent_heat_j_per_kg;
  const double initial_total_h =
      input.committed_parcel.droplet_mass_kg * initial_liquid_h *
      multiplicity;
  double final_total_h = 0.0;
  if (current.droplet_mass_kg > 0.0) {
    const LiquidPropertyReport final_properties =
        input.liquid_properties->evaluate(
            {current.liquid_material_fingerprint, current.temperature_k});
    if (!final_properties.succeeded()) {
      return fixed_failure(FixedExchangeStatus::liquid_property_failure);
    }
    const double average_cp =
        0.5 * (initial_properties.properties.cp_j_per_kg_k +
               final_properties.properties.cp_j_per_kg_k);
    const double final_liquid_h =
        initial_liquid_h +
        average_cp *
            (current.temperature_k -
             input.committed_parcel.temperature_k);
    final_total_h =
        current.droplet_mass_kg * final_liquid_h * multiplicity;
  }
  exchange.parcel_thermochemical_enthalpy_delta_j =
      final_total_h - initial_total_h;
  exchange.parcel_kinetic_energy_delta_j =
      0.5 * multiplicity *
      (current.droplet_mass_kg *
           dot(current.velocity_m_per_s, current.velocity_m_per_s) -
       input.committed_parcel.droplet_mass_kg *
           dot(input.committed_parcel.velocity_m_per_s,
               input.committed_parcel.velocity_m_per_s));
  exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j =
      vapor_absolute_h_integral * multiplicity;
  exchange.convective_heat_to_parcel_j = heat_integral * multiplicity;
  exchange.drag_work_to_parcel_j = drag_work_integral * multiplicity;
  exchange.thermal_exchange_to_gas_j =
      exchange.vapor_absolute_thermochemical_enthalpy_to_gas_j -
      exchange.convective_heat_to_parcel_j;
  exchange.gas_mass_delta_kg =
      -exchange.parcel_liquid_mass_delta_kg;
  for (std::size_t component = 0U; component < 3U; ++component) {
    exchange.gas_momentum_delta_kg_m_per_s[component] =
        exchange.vapor_momentum_to_gas_kg_m_per_s[component] -
        exchange.drag_momentum_to_parcel_kg_m_per_s[component];
    exchange.momentum_closure_residual_kg_m_per_s[component] =
        exchange.parcel_momentum_delta_kg_m_per_s[component] +
        exchange.gas_momentum_delta_kg_m_per_s[component];
    exchange.momentum_quadrature_residual_kg_m_per_s[component] =
        -exchange.parcel_momentum_delta_kg_m_per_s[component] -
        exchange.gas_momentum_delta_kg_m_per_s[component];
  }
  exchange.mass_closure_residual_kg =
      exchange.parcel_liquid_mass_delta_kg + exchange.gas_mass_delta_kg;
  exchange.thermal_exchange_state_residual_j =
      exchange.parcel_thermochemical_enthalpy_delta_j +
      exchange.thermal_exchange_to_gas_j;
  exchange.available = finite_exchange(exchange);
  if (!exchange.available) {
    return fixed_failure(FixedExchangeStatus::non_finite_output);
  }
  report.exchange = exchange;
  report.status = FixedExchangeStatus::success;
  return report;
}

}  // namespace hundun::v04::spray::detail
