// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/chem_composition.hpp"
#include "hundun/mesh_topology.hpp"

#include <optional>
#include <vector>

namespace hundun::chemistry {
class ThermodynamicsService;
}

namespace hundun::flow::detail {

enum class ReactingWallThermalPolicy { adiabatic, isothermal };

struct ReactingWallBoundarySpec final {
  bool catalytic{};
  ReactingWallThermalPolicy thermal_policy{
      ReactingWallThermalPolicy::adiabatic};
  std::optional<double> wall_temperature_k;
};

struct ReactingWallBoundaryState final {
  std::uint32_t patch_id{};
  std::vector<double> species_flux_kg_per_m2_s;
  double heat_flux_w_per_m2{};
  std::optional<double> wall_h_tc_j_per_kg;
};

ReactingWallBoundaryState resolve_reacting_wall_boundary(
    const mesh::BoundaryPatch &, const ReactingWallBoundarySpec &,
    double p0_pa, const std::vector<double> &mass_fractions,
    const chemistry::ThermodynamicsService &);

struct ReactingInletSpec final {
  double p0_pa{};
  std::vector<double> mass_fractions;
  std::optional<double> temperature_k;
  std::optional<double> h_tc_j_per_kg;
};

struct ReactingInletState final {
  std::uint32_t patch_id{};
  chemistry::ThermochemicalPoint thermochemical;
  double temperature_k{};
};

ReactingInletState resolve_reacting_inlet(
    const mesh::BoundaryPatch &, const ReactingInletSpec &,
    const chemistry::ThermodynamicsService &);

chemistry::ThermochemicalPoint reacting_pressure_outlet_state(
    const mesh::BoundaryPatch &, double p0_pa, double mechanical_pressure_pa,
    double h_tc_j_per_kg, const std::vector<double> &mass_fractions,
    const chemistry::ThermodynamicsService &);

} // namespace hundun::flow::detail
