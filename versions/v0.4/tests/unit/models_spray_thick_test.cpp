// SPDX-License-Identifier: Apache-2.0
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "models_spray_transfer_detail.hpp"
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
ThickExchangeInput input() {
  ThickExchangeInput q;
  q.parcel.id = {1, 1};
  q.parcel.liquid_material_fingerprint = 1;
  q.parcel.temperature_k = 350;
  q.parcel.droplet_diameter_m = 30e-6;
  q.parcel.droplet_mass_kg =
      800 * pi / 6 * std::pow(q.parcel.droplet_diameter_m, 3);
  q.parcel.multiplicity = 7;
  q.parcel.velocity_m_per_s = {2, -3, 4};
  q.liquid_properties = {800, 2100, 3e5, 1000, .02, .001};
  q.gas_temperature_k = 900;
  q.gas_cp_j_per_kg_k = 1150;
  q.gas_dynamic_viscosity_pa_s = 3e-5;
  q.prandtl_number = .7;
  q.boiling_temperature_k = 520;
  q.vapor_absolute_thermochemical_enthalpy_j_per_kg = 5e5;
  q.duration_s = 4e-10;
  return q;
}
bool unit() {
  auto q = input();
  const auto initial = q.parcel;
  auto hot = evaluate_thick_exchange(q);
  if (!hot.succeeded() ||
      hot.candidate_parcel.droplet_mass_kg >= initial.droplet_mass_kg ||
      q.parcel.droplet_mass_kg != initial.droplet_mass_kg)
    return false;
  const auto& e = hot.exchange;
  if (e.gas_mass_delta_kg != -e.parcel_liquid_mass_delta_kg ||
      e.thermal_exchange_to_gas_j != -e.parcel_thermochemical_enthalpy_delta_j)
    return false;
  for (unsigned d = 0; d < 3; ++d)
    if (e.gas_momentum_delta_kg_m_per_s[d] !=
        -e.parcel_momentum_delta_kg_m_per_s[d])
      return false;
  auto one = q;
  one.parcel.multiplicity = 1;
  const auto single = evaluate_thick_exchange(one);
  if (!single.succeeded() ||
      single.candidate_parcel.droplet_mass_kg !=
          hot.candidate_parcel.droplet_mass_kg ||
      std::abs(e.gas_mass_delta_kg / 7 - single.exchange.gas_mass_delta_kg) >
          1e-25)
    return false;
  q.gas_temperature_k = q.parcel.temperature_k;
  const auto equilibrium = evaluate_thick_exchange(q);
  if (!equilibrium.succeeded() ||
      equilibrium.candidate_parcel.droplet_mass_kg !=
          q.parcel.droplet_mass_kg ||
      equilibrium.candidate_parcel.temperature_k != q.parcel.temperature_k ||
      equilibrium.exchange.gas_mass_delta_kg != 0)
    return false;
  q.gas_temperature_k = 250;
  const auto cold = evaluate_thick_exchange(q);
  if (!cold.succeeded() ||
      cold.candidate_parcel.droplet_mass_kg != q.parcel.droplet_mass_kg ||
      !(cold.candidate_parcel.temperature_k < q.parcel.temperature_k) ||
      cold.exchange.thermal_exchange_to_gas_j <= 0)
    return false;
  const double capacity =
      q.parcel.droplet_mass_kg * q.liquid_properties.cp_j_per_kg_k;
  const double conductance_dt =
      cold.convective_heat_transfer_w_per_k * q.duration_s;
  const double expected_cold = (capacity * q.parcel.temperature_k +
                                conductance_dt * q.gas_temperature_k) /
                               (capacity + conductance_dt);
  if (std::abs(cold.candidate_parcel.temperature_k - expected_cold) > 1e-12)
    return false;
  q = input();
  q.boiling_temperature_k = q.parcel.temperature_k;
  q.duration_s = 1e-11;
  const auto limited = evaluate_thick_exchange(q);
  if (!limited.succeeded() || !limited.temperature_limited ||
      limited.candidate_parcel.temperature_k !=
          .999 * q.boiling_temperature_k ||
      !(limited.temperature_limit_energy_one_droplet_j > 0) ||
      limited.exchange.thermal_exchange_state_residual_j != 0)
    return false;
  q = input();
  q.duration_s = 2 * q.parcel.droplet_diameter_m * q.parcel.droplet_diameter_m /
                 hot.diameter_squared_loss_rate_m2_per_s;
  const auto vanished = evaluate_thick_exchange(q);
  if (!vanished.succeeded() || !vanished.complete_evaporation ||
      vanished.candidate_parcel.droplet_mass_kg != 0 ||
      vanished.exchange.gas_mass_delta_kg !=
          q.parcel.droplet_mass_kg * q.parcel.multiplicity ||
      !(vanished.event_time_s < q.duration_s))
    return false;
  q = input();
  q.duration_s = 0;
  const auto zero = evaluate_thick_exchange(q);
  if (!zero.succeeded() ||
      zero.candidate_parcel.droplet_mass_kg != q.parcel.droplet_mass_kg ||
      zero.exchange.thermal_exchange_to_gas_j != 0)
    return false;
  q = input();
  q.parcel.droplet_mass_kg *= 1.01;
  if (evaluate_thick_exchange(q).status !=
      ThickExchangeStatus::inconsistent_droplet_geometry)
    return false;
  q = input();
  q.gas_temperature_k = std::numeric_limits<double>::quiet_NaN();
  const auto bad = evaluate_thick_exchange(q);
  return !bad.succeeded() && !bad.exchange.available &&
         bad.candidate_parcel.droplet_mass_kg == 0;
}
}  // namespace
int main(int argc, char**) {
  if (argc == 1) return unit() ? 0 : 1;
  auto q = input();
  unsigned count = 0;
  std::cout << std::scientific << std::setprecision(17);
  while (std::cin >> q.parcel.droplet_diameter_m >> q.parcel.temperature_k >>
         q.gas_temperature_k >> q.gas_cp_j_per_kg_k >>
         q.liquid_properties.cp_j_per_kg_k >>
         q.liquid_properties.latent_heat_j_per_kg >>
         q.liquid_properties.density_kg_per_m3 >>
         q.gas_dynamic_viscosity_pa_s >> q.prandtl_number >>
         q.boiling_temperature_k >>
         q.vapor_absolute_thermochemical_enthalpy_j_per_kg >> q.duration_s) {
    q.parcel.droplet_mass_kg = q.liquid_properties.density_kg_per_m3 * pi / 6 *
                               std::pow(q.parcel.droplet_diameter_m, 3);
    const auto r = evaluate_thick_exchange(q);
    if (!r.succeeded()) return 2;
    const auto& p = r.candidate_parcel;
    const auto& e = r.exchange;
    std::cout << ++count << ' ' << p.droplet_diameter_m << ' '
              << p.droplet_mass_kg << ' ' << p.temperature_k << ' '
              << r.spalding_heat_number << ' '
              << r.diameter_squared_loss_rate_m2_per_s << ' '
              << -r.mean_evaporation_rate_one_droplet_kg_per_s << ' '
              << r.mean_liquid_enthalpy_rate_one_droplet_w << ' '
              << e.mass_closure_residual_kg << ' '
              << e.thermal_exchange_state_residual_j << ' '
              << r.temperature_limit_energy_one_droplet_j << '\n';
  }
  return count ? 0 : 3;
}
