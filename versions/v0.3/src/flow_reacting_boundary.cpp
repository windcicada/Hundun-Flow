// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_reacting_boundary_detail.hpp"

#include "chem_thermodynamics_service_detail.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hundun::flow::detail {
namespace {

double enthalpy_at_temperature(
    const chemistry::ThermodynamicsService &service, double p0_pa,
    double temperature_k, const std::vector<double> &mass_fractions) {
  if (!std::isfinite(p0_pa) || p0_pa <= 0.0 ||
      !std::isfinite(temperature_k) || temperature_k <= 0.0)
    throw std::invalid_argument("reacting boundary thermo input is invalid");
  double h = temperature_k * 1000.0;
  for (std::size_t iteration = 0; iteration < 40U; ++iteration) {
    const auto value = service.evaluate({p0_pa, h, mass_fractions});
    const double error = value.temperature_k - temperature_k;
    if (std::abs(error) <= 1.0e-11 * std::max(1.0, temperature_k))
      return h;
    h -= error * value.cp_j_per_kg_k;
  }
  throw std::invalid_argument("reacting boundary enthalpy inversion failed");
}

void require_composition(const chemistry::ThermodynamicsService &service,
                         const std::vector<double> &mass_fractions) {
  if (mass_fractions.size() != service.composition().species.size())
    throw std::invalid_argument("reacting boundary composition shape mismatch");
  double sum{};
  for (double value : mass_fractions) {
    if (!std::isfinite(value) || value < 0.0)
      throw std::invalid_argument("reacting boundary composition is invalid");
    sum += value;
  }
  if (std::abs(sum - 1.0) > 1.0e-12)
    throw std::invalid_argument("reacting boundary composition is not normalized");
}

} // namespace

ReactingWallBoundaryState resolve_reacting_wall_boundary(
    const mesh::BoundaryPatch &patch, const ReactingWallBoundarySpec &spec,
    double p0_pa, const std::vector<double> &mass_fractions,
    const chemistry::ThermodynamicsService &service) {
  require_composition(service, mass_fractions);
  if (spec.catalytic)
    throw std::invalid_argument("catalytic reacting walls are unsupported");
  ReactingWallBoundaryState result;
  result.patch_id = patch.stable_id();
  result.species_flux_kg_per_m2_s.assign(mass_fractions.size(), 0.0);
  if (spec.thermal_policy == ReactingWallThermalPolicy::adiabatic) {
    if (spec.wall_temperature_k)
      throw std::invalid_argument("adiabatic wall cannot specify temperature");
    result.heat_flux_w_per_m2 = 0.0;
    return result;
  }
  if (!spec.wall_temperature_k)
    throw std::invalid_argument("isothermal wall requires temperature");
  result.wall_h_tc_j_per_kg = enthalpy_at_temperature(
      service, p0_pa, *spec.wall_temperature_k, mass_fractions);
  return result;
}

ReactingInletState resolve_reacting_inlet(
    const mesh::BoundaryPatch &patch, const ReactingInletSpec &spec,
    const chemistry::ThermodynamicsService &service) {
  require_composition(service, spec.mass_fractions);
  if (!std::isfinite(spec.p0_pa) || spec.p0_pa <= 0.0 ||
      (!spec.temperature_k && !spec.h_tc_j_per_kg))
    throw std::invalid_argument("reacting inlet state is incomplete");
  const double h = spec.h_tc_j_per_kg.value_or(enthalpy_at_temperature(
      service, spec.p0_pa, *spec.temperature_k, spec.mass_fractions));
  const auto properties =
      service.evaluate({spec.p0_pa, h, spec.mass_fractions});
  if (spec.temperature_k &&
      std::abs(properties.temperature_k - *spec.temperature_k) >
          1.0e-10 * std::max(1.0, *spec.temperature_k))
    throw std::invalid_argument("reacting inlet redundant state is inconsistent");
  return {patch.stable_id(), {spec.p0_pa, h, spec.mass_fractions},
          properties.temperature_k};
}

chemistry::ThermochemicalPoint reacting_pressure_outlet_state(
    const mesh::BoundaryPatch &patch, double p0_pa,
    double mechanical_pressure_pa, double h_tc_j_per_kg,
    const std::vector<double> &mass_fractions,
    const chemistry::ThermodynamicsService &service) {
  static_cast<void>(patch.stable_id());
  require_composition(service, mass_fractions);
  if (!std::isfinite(mechanical_pressure_pa))
    throw std::invalid_argument("reacting outlet mechanical pressure is invalid");
  chemistry::ThermochemicalPoint point{p0_pa, h_tc_j_per_kg, mass_fractions};
  static_cast<void>(service.evaluate(point));
  return point;
}

} // namespace hundun::flow::detail
