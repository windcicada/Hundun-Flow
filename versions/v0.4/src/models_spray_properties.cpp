// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_transfer_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::spray::detail {
namespace {

constexpr double kUniversalGasConstantJPerKmolK = 8314.46261815324;
bool finite_vector(const Vector3& value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
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

}  // namespace hundun::v04::spray::detail
