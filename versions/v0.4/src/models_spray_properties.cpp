// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include "models_spray_properties_detail.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace hundun::v04::spray::detail {
namespace {

constexpr double kUniversalGasConstantJPerKmolK = 8314.46261815324;
bool finite_vector(const Vector3 &value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
}
LiquidPropertyReport liquid_failure(LiquidPropertyStatus status) noexcept {
  LiquidPropertyReport report;
  report.status = status;
  return report;
}

bool evaluate_temperature_correlation(const TemperatureCorrelation &law,
                                      double temperature_k,
                                      double &value) noexcept {
  value = 0.0;
  if (!std::isfinite(temperature_k) || !(temperature_k > 0.0))
    return false;
  switch (law.kind) {
  case TemperatureCorrelationKind::constant:
    if (!std::isfinite(law.c[0U]))
      return false;
    value = law.c[0U];
    return true;
  case TemperatureCorrelationKind::polynomial_cubic: {
    if (!std::isfinite(law.reference_temperature_k) ||
        !(law.reference_temperature_k > 0.0)) {
      return false;
    }
    for (double coefficient : law.c) {
      if (!std::isfinite(coefficient))
        return false;
    }
    const double theta = temperature_k - law.reference_temperature_k;
    value = ((law.c[3U] * theta + law.c[2U]) * theta + law.c[1U]) * theta +
            law.c[0U];
    return true;
  }
  }
  return false;
}

bool evaluate_saturation_pressure(const SaturationPressureCorrelation &law,
                                  double temperature_k,
                                  double &pressure_pa) noexcept {
  pressure_pa = 0.0;
  if (!std::isfinite(temperature_k) || !(temperature_k > 0.0))
    return false;
  switch (law.kind) {
  case SaturationPressureCorrelationKind::antoine_kelvin: {
    if (!std::isfinite(law.antoine_a) || !std::isfinite(law.antoine_b_k) ||
        !std::isfinite(law.antoine_c_k) ||
        !std::isfinite(law.pressure_scale_pa) ||
        !(law.pressure_scale_pa > 0.0)) {
      return false;
    }
    const double denominator = temperature_k + law.antoine_c_k;
    if (!std::isfinite(denominator) || denominator == 0.0)
      return false;
    const double exponent = law.antoine_a - law.antoine_b_k / denominator;
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
bool valid_liquid_properties(const LiquidProperties &liquid) noexcept {
  return std::isfinite(liquid.density_kg_per_m3) &&
         liquid.density_kg_per_m3 > 0.0 &&
         std::isfinite(liquid.cp_j_per_kg_k) && liquid.cp_j_per_kg_k > 0.0 &&
         std::isfinite(liquid.latent_heat_j_per_kg) &&
         liquid.latent_heat_j_per_kg >= 0.0 &&
         std::isfinite(liquid.saturation_pressure_pa) &&
         liquid.saturation_pressure_pa >= 0.0 &&
         std::isfinite(liquid.surface_tension_n_per_m) &&
         liquid.surface_tension_n_per_m >= 0.0 &&
         std::isfinite(liquid.viscosity_pa_s) && liquid.viscosity_pa_s >= 0.0;
}

} // namespace

LiquidAssetReport load_liquid_asset(const std::filesystem::path &path,
                                    std::uint64_t expected,
                                    const portable::GasIdentity &gas) noexcept {
  LiquidAssetReport out;
  try {
    std::ifstream file(path, std::ios::binary);
    if (!file)
      return out;
    std::string content;
    char ch;
    while (file.get(ch)) {
      if (content.size() == 65536)
        return out;
      content.push_back(ch);
    }
    if (!file.eof())
      return out;
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char c : content) {
      hash ^= c;
      hash *= UINT64_C(1099511628211);
    }
    if (!expected || hash != expected) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    std::istringstream input(content);
    std::string token;
    const auto key = [&](const char *name) {
      if (!(input >> token) || token != name)
        throw std::invalid_argument("liquid asset field");
    };
    LiquidAsset candidate;
    key("HUNDUN_LIQUID_ASSET_V1");
    key("units");
    key("SI");
    key("source");
    input >> candidate.source;
    key("gas_sha");
    input >> token;
    if (token != gas.mechanism_sha256) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    key("gas_phase");
    input >> token;
    if (token != gas.phase) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    key("gas_species");
    std::size_t n;
    input >> n;
    if (!input || n == 0 || n != gas.species_names.size() || n > 65536 ||
        gas.molecular_weights_kg_per_kmol.size() != n ||
        gas.element_names.empty() ||
        gas.element_counts.size() != n * gas.element_names.size() ||
        !gas.composition_fingerprint || !gas.closure_fingerprint) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    for (std::size_t i = 0; i < n; ++i) {
      input >> token;
      if (token != gas.species_names[i]) {
        out.status = portable::Status::identity_mismatch;
        return out;
      }
    }
    key("enthalpy_reference");
    input >> token;
    if (token != gas.enthalpy_reference) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    key("vapor");
    input >> candidate.vapor_species_name;
    auto found = std::find(gas.species_names.begin(), gas.species_names.end(),
                           candidate.vapor_species_name);
    if (found == gas.species_names.end()) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    candidate.vapor_species_index =
        static_cast<std::size_t>(found - gas.species_names.begin());
    key("molecular_weight");
    input >> candidate.vapor_molecular_weight_kg_per_kmol;
    if (candidate.vapor_molecular_weight_kg_per_kmol !=
        gas.molecular_weights_kg_per_kmol[candidate.vapor_species_index]) {
      out.status = portable::Status::identity_mismatch;
      return out;
    }
    auto &pack = candidate.pack;
    key("temperature_range");
    input >> pack.minimum_temperature_k >> pack.maximum_temperature_k;
    key("liquid_reference");
    input >> candidate.reference_temperature_k >>
        candidate.reference_liquid_enthalpy_j_per_kg;
    const auto correlation = [&](const char *name,
                                 TemperatureCorrelation &law) {
      key(name);
      input >> token;
      if (token == "constant")
        law.kind = TemperatureCorrelationKind::constant;
      else if (token == "cubic")
        law.kind = TemperatureCorrelationKind::polynomial_cubic;
      else
        throw std::invalid_argument("unknown liquid correlation");
      input >> law.reference_temperature_k;
      for (double &c : law.c)
        input >> c;
      if (!input)
        throw std::invalid_argument("invalid liquid coefficient");
      for (double c : law.c)
        if (!std::isfinite(c))
          throw std::invalid_argument("nonfinite coefficient");
      if (law.kind == TemperatureCorrelationKind::constant &&
          (law.c[1] != 0 || law.c[2] != 0 || law.c[3] != 0))
        throw std::invalid_argument("unused nonzero constant coefficient");
    };
    correlation("density", pack.density_kg_per_m3);
    correlation("cp", pack.cp_j_per_kg_k);
    correlation("latent", pack.latent_heat_j_per_kg);
    correlation("surface_tension", pack.surface_tension_n_per_m);
    correlation("viscosity", pack.viscosity_pa_s);
    key("saturation");
    input >> token;
    auto &sat = pack.saturation_pressure;
    if (token == "antoine") {
      sat.kind = SaturationPressureCorrelationKind::antoine_kelvin;
      input >> sat.antoine_a >> sat.antoine_b_k >> sat.antoine_c_k >>
          sat.pressure_scale_pa;
      if (pack.minimum_temperature_k + sat.antoine_c_k <= 0 &&
          pack.maximum_temperature_k + sat.antoine_c_k >= 0)
        return out;
    } else if (token == "clausius") {
      sat.kind = SaturationPressureCorrelationKind::clausius_clapeyron;
      input >> sat.reference_pressure_pa >> sat.reference_temperature_k >>
          sat.latent_heat_j_per_kg >> sat.molecular_weight_kg_per_kmol;
      if (sat.molecular_weight_kg_per_kmol !=
          candidate.vapor_molecular_weight_kg_per_kmol)
        return out;
    } else
      return out;
    key("end");
    if (input >> token)
      return out;
    candidate.content_fingerprint = hash;
    pack.material_fingerprint = hash;
    if (candidate.source.empty() ||
        !std::isfinite(pack.minimum_temperature_k) ||
        pack.minimum_temperature_k <= 0 ||
        !std::isfinite(pack.maximum_temperature_k) ||
        pack.maximum_temperature_k < pack.minimum_temperature_k)
      return out;
    // Check all extrema of every cubic, not only sampled temperatures.
    LiquidPropertyService service(&pack, 1);
    const auto valid_t = [&](double t) {
      return service.evaluate({hash, t}).succeeded() &&
             evaluate_liquid_enthalpy(candidate, t).available;
    };
    if (!valid_t(pack.minimum_temperature_k) ||
        !valid_t(pack.maximum_temperature_k))
      return out;
    const TemperatureCorrelation *laws[]{
        &pack.density_kg_per_m3, &pack.cp_j_per_kg_k,
        &pack.latent_heat_j_per_kg, &pack.surface_tension_n_per_m,
        &pack.viscosity_pa_s};
    for (const auto *law : laws) {
      if (law->kind != TemperatureCorrelationKind::polynomial_cubic)
        continue;
      const double a = 3 * law->c[3], b = 2 * law->c[2], c = law->c[1];
      double roots[2]{};
      std::size_t count = 0;
      if (a == 0) {
        if (b != 0)
          roots[count++] = -c / b;
      } else {
        const double disc = b * b - 4 * a * c;
        if (disc >= 0) {
          roots[count++] = (-b - std::sqrt(disc)) / (2 * a);
          roots[count++] = (-b + std::sqrt(disc)) / (2 * a);
        }
      }
      for (std::size_t i = 0; i < count; ++i) {
        const double t = law->reference_temperature_k + roots[i];
        if (t >= pack.minimum_temperature_k &&
            t <= pack.maximum_temperature_k && !valid_t(t))
          return out;
      }
    }
    candidate.gas_identity = gas;
    out.asset = std::move(candidate);
    out.status = portable::Status::success;
    out.available = true;
    return out;
  } catch (const std::bad_alloc &) {
    out.status = portable::Status::capacity_exceeded;
  } catch (...) {
    out.status = portable::Status::invalid_input;
  }
  return out;
}

LiquidEnthalpyReport evaluate_liquid_enthalpy(const LiquidAsset &asset,
                                              double t) noexcept {
  LiquidEnthalpyReport out;
  const auto &law = asset.pack.cp_j_per_kg_k;
  if (!std::isfinite(t) || t < asset.pack.minimum_temperature_k ||
      t > asset.pack.maximum_temperature_k ||
      !std::isfinite(asset.reference_temperature_k) ||
      asset.reference_temperature_k < asset.pack.minimum_temperature_k ||
      asset.reference_temperature_k > asset.pack.maximum_temperature_k ||
      !std::isfinite(asset.reference_liquid_enthalpy_j_per_kg))
    return out;
  double cp;
  if (!evaluate_temperature_correlation(law, t, cp) || cp <= 0)
    return out;
  double integral = 0;
  if (law.kind == TemperatureCorrelationKind::constant)
    integral = law.c[0] * (t - asset.reference_temperature_k);
  else {
    const double a = t - law.reference_temperature_k,
                 b = asset.reference_temperature_k -
                     law.reference_temperature_k;
    for (std::size_t i = 0; i < 4; ++i)
      integral += law.c[i] / static_cast<double>(i + 1) *
                  (std::pow(a, i + 1) - std::pow(b, i + 1));
  }
  const double h = asset.reference_liquid_enthalpy_j_per_kg + integral;
  if (!std::isfinite(h))
    return out;
  out.status = portable::Status::success;
  out.available = true;
  out.liquid_enthalpy_j_per_kg = h;
  out.cp_j_per_kg_k = cp;
  return out;
}

FilmQueryWorkspace::FilmQueryWorkspace(std::size_t n)
    : film_y_(n), surface_y_(n), pure_y_(n), d_(n), h_(n), w_(n) {
  if (!n)
    throw std::invalid_argument("film workspace needs species");
}
const double *FilmQueryWorkspace::film_mass_fractions() const noexcept {
  return available_ ? film_y_.data() : nullptr;
}
FilmQueryReport
FilmQueryWorkspace::query(const LiquidAsset &asset,
                          portable::GasQueryProvider &provider,
                          const FilmQueryInput &input) noexcept {
  available_ = false;
  FilmQueryReport out;
  const auto &id = provider.gas_identity();
  const auto &expected = asset.gas_identity;
  const auto n = id.species_names.size(), v = asset.vapor_species_index;
  if (input.far_gas.revision != input.expected_revision) {
    out.status = portable::Status::stale_revision;
    return out;
  }
  if (id.mechanism_sha256 != expected.mechanism_sha256 ||
      id.phase != expected.phase ||
      id.species_names != expected.species_names ||
      id.element_names != expected.element_names ||
      id.element_counts != expected.element_counts ||
      id.molecular_weights_kg_per_kmol !=
          expected.molecular_weights_kg_per_kmol ||
      id.enthalpy_reference != expected.enthalpy_reference ||
      id.composition_fingerprint != expected.composition_fingerprint ||
      id.closure_fingerprint != expected.closure_fingerprint ||
      input.far_gas.composition_fingerprint != id.composition_fingerprint) {
    out.status = portable::Status::identity_mismatch;
    return out;
  }
  if (n != film_y_.size()) {
    out.status = portable::Status::capacity_exceeded;
    return out;
  }
  if (v >= n || input.far_gas.species_count != n ||
      !input.far_gas.mass_fractions ||
      input.far_gas.coordinates !=
          portable::GasStateCoordinates::pressure_enthalpy ||
      !finite_vector(input.velocity_m_per_s))
    return out;
  const auto liquid_h =
      evaluate_liquid_enthalpy(asset, input.surface_temperature_k);
  LiquidPropertyService liquid_service(&asset.pack, 1);
  const auto liquid = liquid_service.evaluate(
      {asset.pack.material_fingerprint, input.surface_temperature_k});
  if (!liquid_h.available || !liquid.succeeded())
    return out;
  portable::GasQueryOutput response{{}, d_.data(), h_.data(), w_.data(), n};
  const auto call = [&](const portable::GasQuery &q) {
    const auto status = provider.query_gas(q, response);
    if (status != portable::Status::success)
      return status;
    if (response.sample.revision != input.expected_revision)
      return portable::Status::stale_revision;
    if (response.sample.composition_fingerprint != id.composition_fingerprint)
      return portable::Status::identity_mismatch;
    return portable::Status::success;
  };
  out.status = call(input.far_gas);
  if (out.status != portable::Status::success)
    return out;
  const auto far = response.sample;
  const double far_v = input.far_gas.mass_fractions[v];
  if (!std::isfinite(far_v) || far_v < 0 || far_v >= 1) {
    out.status = portable::Status::invalid_input;
    return out;
  }
  double carrier_inverse_mw = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double y = input.far_gas.mass_fractions[i],
                 mw = id.molecular_weights_kg_per_kmol[i];
    if (!std::isfinite(y) || y < 0 || y > 1 || !std::isfinite(mw) || mw <= 0) {
      out.status = portable::Status::invalid_input;
      return out;
    }
    if (i != v)
      carrier_inverse_mw += y / mw;
  }
  const double mole_v =
      liquid.properties.saturation_pressure_pa / far.pressure_pa;
  if (!std::isfinite(mole_v) || mole_v < 0 || mole_v >= 1 ||
      carrier_inverse_mw <= 0) {
    out.status = portable::Status::unavailable;
    return out; // boiling/critical conditions outside scope
  }
  const double carrier_mw = (1 - far_v) / carrier_inverse_mw;
  const double surface_v = mole_v * asset.vapor_molecular_weight_kg_per_kmol /
                           (mole_v * asset.vapor_molecular_weight_kg_per_kmol +
                            (1 - mole_v) * carrier_mw);
  if (surface_v < far_v) {
    out.status = portable::Status::unavailable;
    return out;
  } // condensation not this model
  std::fill(pure_y_.begin(), pure_y_.end(), 0);
  pure_y_[v] = 1;
  auto query = input.far_gas;
  query.coordinates = portable::GasStateCoordinates::pressure_temperature;
  query.temperature_k = input.surface_temperature_k;
  query.mass_fractions = pure_y_.data();
  out.status = call(query);
  if (out.status != portable::Status::success)
    return out;
  const double vapor_h = h_[v], vapor_cp = response.sample.cp_j_per_kg_k;
  const double residual = vapor_h - liquid_h.liquid_enthalpy_j_per_kg -
                          liquid.properties.latent_heat_j_per_kg;
  const auto &latent = asset.pack.latent_heat_j_per_kg;
  const double theta =
      input.surface_temperature_k - latent.reference_temperature_k;
  const double dlatent = latent.kind == TemperatureCorrelationKind::constant
                             ? 0
                             : latent.c[1] + 2 * latent.c[2] * theta +
                                   3 * latent.c[3] * theta * theta;
  if (!std::isfinite(residual) ||
      std::abs(residual) >
          1e-6 + 1e-9 * std::max(std::abs(vapor_h),
                                 std::abs(liquid_h.liquid_enthalpy_j_per_kg)) ||
      !std::isfinite(vapor_cp) ||
      std::abs(vapor_cp - liquid_h.cp_j_per_kg_k - dlatent) >
          1e-6 + 1e-9 * std::abs(vapor_cp)) {
    out.status = portable::Status::conservation_failure;
    return out;
  }
  for (std::size_t i = 0; i < n; ++i) {
    surface_y_[i] = i == v ? surface_v
                           : input.far_gas.mass_fractions[i] * (1 - surface_v) /
                                 (1 - far_v);
    film_y_[i] = (2 * surface_y_[i] + input.far_gas.mass_fractions[i]) / 3;
  }
  query.mass_fractions = film_y_.data();
  query.temperature_k =
      (2 * input.surface_temperature_k + far.temperature_k) / 3;
  out.status = call(query);
  if (out.status != portable::Status::success)
    return out;
  const OneThirdFilmInput values{input.velocity_m_per_s,
                                 input.surface_temperature_k,
                                 far.temperature_k,
                                 surface_v,
                                 far_v,
                                 far.pressure_pa,
                                 response.sample.density_kg_per_m3,
                                 response.sample.viscosity_pa_s,
                                 response.sample.conductivity_w_per_m_k,
                                 d_[v],
                                 response.sample.cp_j_per_kg_k};
  out.film = sample_one_third_film(values);
  if (!out.film.succeeded()) {
    out.film = {};
    out.status = portable::Status::provider_failure;
    return out;
  }
  out.liquid_enthalpy_j_per_kg = liquid_h.liquid_enthalpy_j_per_kg;
  out.vapor_enthalpy_j_per_kg = vapor_h;
  out.latent_heat_consistency_residual_j_per_kg = residual;
  out.far_gas_density_kg_per_m3 = far.density_kg_per_m3;
  out.available = true;
  available_ = true;
  out.status = portable::Status::success;
  return out;
}

LiquidPropertyReport LiquidPropertyService::evaluate(
    const LiquidPropertyQuery &query) const noexcept {
  if (query.material_fingerprint == 0U || !std::isfinite(query.temperature_k) ||
      !(query.temperature_k > 0.0) || (count_ > 0U && packs_ == nullptr)) {
    return liquid_failure(LiquidPropertyStatus::invalid_input);
  }

  const LiquidPropertyPack *selected = nullptr;
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
      selected->cp_j_per_kg_k, query.temperature_k, properties.cp_j_per_kg_k);
  evaluated &= evaluate_temperature_correlation(
      selected->latent_heat_j_per_kg, query.temperature_k,
      properties.latent_heat_j_per_kg);
  evaluated &= evaluate_temperature_correlation(
      selected->surface_tension_n_per_m, query.temperature_k,
      properties.surface_tension_n_per_m);
  evaluated &= evaluate_temperature_correlation(
      selected->viscosity_pa_s, query.temperature_k, properties.viscosity_pa_s);
  evaluated &= evaluate_saturation_pressure(selected->saturation_pressure,
                                            query.temperature_k,
                                            properties.saturation_pressure_pa);
  if (!evaluated) {
    return liquid_failure(LiquidPropertyStatus::correlation_domain_error);
  }
  const double values[]{
      properties.density_kg_per_m3,       properties.cp_j_per_kg_k,
      properties.latent_heat_j_per_kg,    properties.saturation_pressure_pa,
      properties.surface_tension_n_per_m, properties.viscosity_pa_s};
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

OneThirdFilmSample
sample_one_third_film(const OneThirdFilmInput &input) noexcept {
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
  sample.dynamic_viscosity_pa_s = input.reference_dynamic_viscosity_pa_s;
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

} // namespace hundun::v04::spray::detail
