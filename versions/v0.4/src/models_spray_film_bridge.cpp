// SPDX-License-Identifier: Apache-2.0
#include "models_spray_film_bridge_detail.hpp"
#include <cmath>
namespace hundun::v04::spray::detail {
FilmEnvironmentBridge::FilmEnvironmentBridge(
    const LiquidAsset &asset, portable::GasQueryProvider &gas,
    const ParcelGasStateProvider &sampler, portable::Revision revision,
    EvaporationModel evaporation, const FilmTransportProvider *transport)
    : asset_(asset), gas_(gas), sampler_(sampler), revision_(revision),
      liquid_(&asset.pack, 1), film_(asset.gas_identity.species_names.size()),
      thick_film_(asset.gas_identity.species_names.size()),
      evaporation_(evaporation), transport_(transport),
      y_(asset.gas_identity.species_names.size()) {}
bool FilmEnvironmentBridge::compatible() const noexcept {
  if (evaporation_ == EvaporationModel::abramzon_sirignano) return true;
  return evaporation_ == EvaporationModel::thick_exchange && transport_ &&
      asset_.pack.density_kg_per_m3.kind == TemperatureCorrelationKind::kerosene_density_v1 &&
      asset_.pack.cp_j_per_kg_k.kind == TemperatureCorrelationKind::kerosene_cp_v1 &&
      asset_.pack.latent_heat_j_per_kg.kind == TemperatureCorrelationKind::kerosene_latent_v1 &&
      asset_.pack.saturation_pressure.kind == SaturationPressureCorrelationKind::kerosene_v1;
}
FilmEnvironmentReport
FilmEnvironmentBridge::query(const SprayParcelState &parcel, double elapsed,
                             ParcelPass pass,
                             portable::Revision revision) const noexcept {
  FilmEnvironmentReport out;
  if (revision != revision_) {
    out.status = portable::Status::stale_revision;
    return out;
  }
  if (!compatible() || !std::isfinite(elapsed) || elapsed < 0 ||
      parcel.liquid_material_fingerprint != asset_.pack.material_fingerprint ||
      (pass != ParcelPass::predictor && pass != ParcelPass::corrector))
    return out;
  auto far =
      sampler_.sample(parcel, elapsed, pass, revision, y_.data(), y_.size());
  if (far.status != portable::Status::success) {
    out.status = far.status;
    return out;
  }
  if (far.revision != revision) {
    out.status = portable::Status::stale_revision;
    return out;
  }
  if (far.species_count != y_.size()) {
    out.status = portable::Status::identity_mismatch;
    return out;
  }
  FilmQueryInput request;
  request.expected_revision = revision;
  request.surface_temperature_k = parcel.temperature_k;
  request.velocity_m_per_s = far.velocity_m_per_s;
  request.far_gas = {revision,
                     far.composition_fingerprint,
                     portable::GasStateCoordinates::pressure_enthalpy,
                     far.pressure_pa,
                     far.enthalpy_j_per_kg,
                     0,
                     y_.data(),
                     y_.size()};
  if (evaporation_ == EvaporationModel::thick_exchange) {
    KeroseneFilmInput input;
    input.far_gas = request.far_gas;
    input.expected_revision = input.transport_revision = revision;
    input.surface_temperature_k = parcel.temperature_k;
    input.transport = transport_;
    const auto f = thick_film_.query(asset_, gas_, input);
    const auto h = evaluate_liquid_enthalpy(asset_, parcel.temperature_k);
    if (!f.available || !h.available) {
      out.status = !f.available ? f.status : h.status;
      return out;
    }
    for (double u : far.velocity_m_per_s)
      if (!std::isfinite(u)) return out;
    auto &e = out.environment;
    e.revision = revision;
    e.gas.gas_velocity_m_per_s = far.velocity_m_per_s;
    e.gas.gas_temperature_k = f.far_temperature_k;
    e.gas.gas_vapor_mass_fraction = y_[asset_.vapor_species_index];
    e.gas.thermodynamic_pressure_pa = far.pressure_pa;
    e.gas.film_dynamic_viscosity_pa_s = f.gas_dynamic_viscosity_pa_s;
    e.gas.film_cp_j_per_kg_k = f.gas_cp_j_per_kg_k;
    e.gas.vapor_absolute_thermochemical_enthalpy_j_per_kg = f.vapor_absolute_enthalpy_j_per_kg;
    e.liquid = &liquid_;
    e.liquid_absolute_enthalpy_j_per_kg = h.liquid_enthalpy_j_per_kg;
    e.far_gas_density_kg_per_m3 = f.far_density_kg_per_m3;
    e.far_gas_dynamic_viscosity_pa_s = f.far_dynamic_viscosity_pa_s;
    e.boiling_temperature_k = f.liquid.boiling_temperature_k;
    e.vapor_prandtl_number = f.vapor_prandtl_number;
    e.available = true;
    out.status = portable::Status::success;
    return out;
  }
  const auto sampled = film_.query(asset_, gas_, request);
  if (!sampled.available) {
    out.status = sampled.status;
    return out;
  }
  const auto v = asset_.vapor_species_index;
  double inverse_carrier_mw = 0;
  for (std::size_t i = 0; i < y_.size(); ++i)
    if (i != v)
      inverse_carrier_mw +=
          y_[i] / asset_.gas_identity.molecular_weights_kg_per_kmol[i];
  if (!std::isfinite(inverse_carrier_mw) || inverse_carrier_mw <= 0)
    return out;
  const auto &f = sampled.film;
  auto &environment = out.environment;
  environment.revision = revision;
  environment.gas = {f.gas_velocity_m_per_s,
                     f.far_gas_temperature_k,
                     f.far_gas_vapor_mass_fraction,
                     f.thermodynamic_pressure_pa,
                     f.density_kg_per_m3,
                     f.dynamic_viscosity_pa_s,
                     f.thermal_conductivity_w_per_m_k,
                     f.vapor_diffusivity_m2_per_s,
                     f.cp_j_per_kg_k,
                     asset_.vapor_molecular_weight_kg_per_kmol,
                     (1 - y_[v]) / inverse_carrier_mw,
                     sampled.vapor_enthalpy_j_per_kg};
  environment.liquid = &liquid_;
  environment.liquid_absolute_enthalpy_j_per_kg =
      sampled.liquid_enthalpy_j_per_kg;
  environment.far_gas_density_kg_per_m3 = sampled.far_gas_density_kg_per_m3;
  environment.available = true;
  out.status = portable::Status::success;
  return out;
}
ParcelTransferEnvironment
FilmEnvironmentBridge::sample(const SprayParcelState &parcel, double elapsed,
                              ParcelPass pass,
                              portable::Revision revision) const noexcept {
  return query(parcel, elapsed, pass, revision).environment;
}
} // namespace hundun::v04::spray::detail
