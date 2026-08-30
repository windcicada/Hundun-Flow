// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_product.hpp"

#include <array>
#include <cstdint>

namespace hundun::v04::test {

inline ValidatedModel product_model(Int3 cells = {17, 11, 7}) {
  ValidatedModel model;
  model.mesh.kind = GeometryKind::uniform;
  model.mesh.lower = {0.0, 0.0, 0.0};
  model.mesh.upper = {2.0, 1.0, 0.5};
  model.mesh.has_exact_cells = true;
  model.mesh.exact_cells = cells;
  model.mesh.minimum_spacing = {
      2.0 / static_cast<double>(cells.x),
      1.0 / static_cast<double>(cells.y),
      0.5 / static_cast<double>(cells.z)};
  model.mesh.max_growth_ratio = 1.0;
  model.mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(cells.x) *
      static_cast<std::uint64_t>(cells.y) *
      static_cast<std::uint64_t>(cells.z);
  model.mesh.limits.max_memory_bytes_per_rank = UINT64_C(1073741824);
  model.pressure_reference = PressureReferenceKind::closed_mass;
  model.turbulence = TurbulenceKind::vreman_wall_function;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
    face.relaxation = 1.0;
  }
  model.schemes.momentum = ConvectionScheme::limited_central2;
  model.schemes.enthalpy = ConvectionScheme::limited_central2;
  model.schemes.species = ConvectionScheme::tvd2;
  model.schemes.passive_scalar = ConvectionScheme::tvd2;
  model.schemes.diffusion = DiffusionScheme::central2;
  model.schemes.limiter = 1.0;
  model.time.control = TimeControlKind::adaptive_flow;
  model.time.scheme = TimeScheme::variable_bdf2;
  model.time.initial_dt = 1.0e-3;
  model.time.minimum_dt = 1.0e-8;
  model.time.maximum_dt = 0.1;
  model.time.convective_cfl = 0.8;
  model.time.viscous_cfl = 0.5;
  model.time.thermal_cfl = 0.5;
  model.time.species_cfl = 0.5;
  model.time.acoustic_cfl = 0.8;
  model.time.maximum_growth = 1.2;
  model.time.retry_factor = 0.5;
  model.time.maximum_retries = 6U;
  model.time.minimum_bdf_ratio = 0.25;
  model.time.maximum_bdf_ratio = 4.0;
  model.thermophysics.data_file = "analytic.d";
  model.thermophysics.minimum_temperature = 200.0;
  model.thermophysics.maximum_temperature = 2000.0;
  model.thermophysics.temperature_relative_tolerance = 1.0e-12;
  model.thermophysics.maximum_temperature_iterations = 64U;
  model.thermophysics.closed_mass_relative_tolerance = 1.0e-12;
  model.thermophysics.maximum_closed_mass_iterations = 32U;
  model.thermophysics.maximum_closed_mass_relative_step = 0.2;
  SpeciesThermophysicalSpec air;
  air.stable_name = "air";
  air.molecular_weight = 28.96546;
  air.temperature_switch = 1000.0;
  air.nasa7_low[0U] = 3.5;
  air.nasa7_high[0U] = 3.5;
  air.viscosity_reference = 1.8e-5;
  air.conductivity = 0.026;
  model.thermophysics.species.push_back(air);
  model.fingerprint = UINT64_C(0x180000001);
  return model;
}

inline constexpr std::array<ProductFreezePhase, 10U> kFreezeOrder{{
    ProductFreezePhase::geometry_and_decomposition,
    ProductFreezePhase::capability_registration,
    ProductFreezePhase::logical_analysis,
    ProductFreezePhase::schema_and_allocation,
    ProductFreezePhase::plan_instantiation,
    ProductFreezePhase::numeric_capacity,
    ProductFreezePhase::communication_binding,
    ProductFreezePhase::view_and_graph_binding,
    ProductFreezePhase::validation,
    ProductFreezePhase::sealed}};

}  // namespace hundun::v04::test
