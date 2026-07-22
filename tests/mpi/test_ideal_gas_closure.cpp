// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/ideal_gas_closure.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

namespace {

hundun::runtime::Int3 grid(int ranks) { return {ranks, 1, 1}; }

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::ideal_gas;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.physics.cp_J_per_kg_K = 1000.0;
  config.physics.gas_constant_J_per_kg_K = 287.05;
  config.physics.thermodynamic_pressure_pa = 101325.0;
  constexpr std::array names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t patch = 0; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::config::FlowCaseConfig open_case() {
  auto config = periodic_case();
  config.boundaries[0].type = hundun::config::BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s = {1.0, 0.0, 0.0};
  config.boundaries[0].thermal_authority =
      hundun::config::InletThermalAuthority::temperature;
  config.boundaries[0].temperature_K = 300.0;
  config.boundaries[0].enthalpy_J_per_kg = 300000.0;
  config.boundaries[0].density_kg_per_m3 = 101325.0 / (287.05 * 300.0);
  config.boundaries[1].type = hundun::config::BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = 0.0;
  for (std::size_t patch = 2U; patch < config.boundaries.size(); ++patch)
    config.boundaries[patch].type = hundun::config::BoundaryType::symmetry;
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name,
          unit,
          "task21",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name, const char *unit,
                                      std::uint32_t components) {
  return {name,
          unit,
          "task21",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

void run(const hundun::runtime::MpiContext &mpi) {
  constexpr hundun::runtime::Int3 extent{8, 4, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
  fields.velocity = registry.declare_field(cell("velocity", "m/s", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell("pi", "Pa", 1U, false));
  fields.face_velocity =
      registry.declare_field(face("face_velocity", "m/s", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3", 1U, true));
  fields.transported_cell_fields = {rho_h};
  registry.freeze();

  auto state = hundun::flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields,
      {0U, 0.0, 1.0e-3, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  constexpr double cp = 1000.0;
  constexpr double gas_constant = 287.05;
  constexpr double temperature = 300.0;
  constexpr double pressure = 101325.0;
  const double density = pressure / (gas_constant * temperature);
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), density);
  initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields = {std::vector<double>(
      topology.owned_cell_count(), density * cp * temperature)};
  state.seed_accepted_layers(initial, initial);

  const auto before = mpi.fp64_reduction_counters();
  auto closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, boundaries, mpi, registry, fields, state,
      {rho_h, cp, gas_constant, pressure});
  const auto after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(after.collective_calls - before.collective_calls == 2U);
  HUNDUN_CHECK(after.reduced_scalars - before.reduced_scalars == 14U);
  HUNDUN_CHECK(after.logical_payload_bytes - before.logical_payload_bytes ==
               112U);
  const auto closure_state = closure.state();
  HUNDUN_CHECK(closure_state.mode ==
               hundun::flow::IdealGasPressureMode::closed_dynamic);
  HUNDUN_CHECK(closure_state.revision == 0U);
  HUNDUN_CHECK(closure_state.target_mass_kg.has_value());
  HUNDUN_CHECK_NEAR(closure_state.thermodynamic_pressure_pa, pressure,
                    1.0e-12 * pressure);
  HUNDUN_CHECK_NEAR(*closure_state.target_mass_kg, density, 5.0e-12 * density);

  auto moved = std::move(closure);
  bool moved_from_rejected = false;
  try {
    static_cast<void>(closure.state());
  } catch (const hundun::runtime::Error &error) {
    moved_from_rejected = std::string_view(error.what()) ==
                          "ideal-gas closure has been moved from";
  }
  HUNDUN_CHECK(moved_from_rejected);
  HUNDUN_CHECK(moved.state().revision == 0U);

  auto open_decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {false, false, false},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology open_topology(open_decomposition);
  hundun::mesh::MeshGeometry open_geometry(
      open_topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 0.25, 0.25}));
  auto open_boundaries =
      hundun::boundary::BoundaryRegistry::create(open_case(), open_topology);
  auto open_state = hundun::flow::FlowState::create(
      registry,
      {open_decomposition.local_extent(), open_topology.local_face_count()},
      fields,
      {0U, 0.0, 6.25e-3, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  auto open_initial = initial;
  open_initial.density.assign(open_topology.owned_cell_count(), density);
  open_initial.velocity.assign(open_topology.owned_cell_count() * 3U, 0.0);
  for (std::size_t cell_index = 0;
       cell_index < open_topology.owned_cell_count(); ++cell_index)
    open_initial.velocity[cell_index * 3U] = 1.0;
  open_initial.mechanical_pressure.assign(open_topology.owned_cell_count(),
                                          0.0);
  open_initial.face_velocity.assign(open_topology.local_face_count() * 3U, 0.0);
  open_initial.face_mass_flux.assign(open_topology.local_face_count(), 0.0);
  open_initial.transported_cell_fields = {std::vector<double>(
      open_topology.owned_cell_count(), density * cp * temperature)};
  open_state.seed_accepted_layers(open_initial, open_initial);
  const auto open_before = mpi.fp64_reduction_counters();
  auto open_closure = hundun::flow::IdealGasClosure::create(
      open_topology, open_geometry, open_boundaries, mpi, registry, fields,
      open_state, {rho_h, cp, gas_constant, pressure});
  const auto open_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(open_after.collective_calls - open_before.collective_calls ==
               2U);
  HUNDUN_CHECK(open_after.reduced_scalars - open_before.reduced_scalars == 14U);
  const auto open_closure_state = open_closure.state();
  HUNDUN_CHECK(open_closure_state.mode ==
               hundun::flow::IdealGasPressureMode::open_fixed);
  HUNDUN_CHECK(!open_closure_state.target_mass_kg.has_value());
  HUNDUN_CHECK(open_closure_state.thermodynamic_pressure_pa == pressure);
  HUNDUN_CHECK(open_closure_state.revision == 0U);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    run(mpi);
  });
}
