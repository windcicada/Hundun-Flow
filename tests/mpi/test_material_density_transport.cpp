// SPDX-License-Identifier: Apache-2.0

#include "flow/src/material_density_transport_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/flow_state.hpp"
#include "hundun/flow/material_density_transport.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace {

hundun::runtime::Int3 grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {4, 1, 1};
  throw hundun::runtime::Error("unsupported task19 rank count");
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::material;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back({"alpha", 0.0});
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::config::FlowCaseConfig open_case() {
  auto config = periodic_case();
  config.boundaries[0].type = hundun::config::BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s = {1.0, 0.0, 0.0};
  config.boundaries[0].thermal_authority =
      hundun::config::InletThermalAuthority::enthalpy;
  config.boundaries[0].enthalpy_J_per_kg = 1.0;
  config.boundaries[0].density_kg_per_m3 = 1.0;
  config.boundaries[0].scalar_values = {{"alpha", 0.2}};
  config.boundaries[1].type = hundun::config::BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = 0.0;
  for (std::size_t patch = 2U; patch < config.boundaries.size(); ++patch)
    config.boundaries[patch].type = hundun::config::BoundaryType::symmetry;
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit) {
  return {name,
          unit,
          "material-density-transport",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          1U,
          2,
          true,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor cell_components(const char *name,
                                                 std::uint32_t components) {
  auto result = cell(name, "1");
  result.components = components;
  result.conservative = false;
  return result;
}

hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name,
          "1",
          "material-density-transport",
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
  fields.density = registry.declare_field(cell("rho", "kg/m3"));
  fields.velocity = registry.declare_field(cell_components("velocity", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("pi", "Pa"));
  fields.face_velocity = registry.declare_field(face("u_face", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3"));
  const auto rho_phi = registry.declare_field(cell("rho_phi_alpha", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_phi};
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), 1.0);
  initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 19.0);
  initial.transported_cell_fields.resize(2U);
  initial.transported_cell_fields[0].assign(topology.owned_cell_count(), 3.0);
  initial.transported_cell_fields[1].assign(topology.owned_cell_count(), 0.4);
  state.seed_accepted_layers(initial, initial);

  hundun::runtime::FieldStorage supplied_flux(
      registry,
      hundun::runtime::FieldLayoutSet{local, topology.local_face_count()});
  constexpr hundun::runtime::PhaseId phase = 1901U;
  constexpr hundun::runtime::ActorId actor = 1901U;
  hundun::runtime::FieldAccessPlan flux_access(registry);
  flux_access.declare_access(phase, actor, fields.face_mass_flux,
                             hundun::runtime::AccessMode::read_write);
  flux_access.freeze();
  auto flux_values = supplied_flux.acquire_face_write<double>(
      flux_access, phase, actor, fields.face_mass_flux);
  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index) {
    flux_values(face_index, 0) = 0.0;
  }

  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::flow::MaterialDensityTransportSpec spec;
  spec.enthalpy_density = rho_h;
  spec.enthalpy_diffusivity_kg_per_m_s = 0.0;
  spec.scalar_densities = {rho_phi};
  spec.scalar_diffusivities_kg_per_m_s = {0.0};
  auto transport = hundun::flow::MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, halo,
      fields, spec);
  const auto expect_create_rejected =
      [&](hundun::flow::FlowFieldIds candidate_fields,
          hundun::flow::MaterialDensityTransportSpec candidate_spec) {
        bool rejected = false;
        try {
          static_cast<void>(hundun::flow::MaterialDensityTransport::create(
              registry, decomposition, topology, geometry, boundaries, mpi,
              halo, std::move(candidate_fields), std::move(candidate_spec)));
        } catch (const hundun::runtime::Error &) {
          rejected = true;
        }
        HUNDUN_CHECK(rejected);
      };
  auto invalid_spec = spec;
  invalid_spec.enthalpy_diffusivity_kg_per_m_s =
      std::numeric_limits<double>::quiet_NaN();
  expect_create_rejected(fields, invalid_spec);
  invalid_spec = spec;
  invalid_spec.scalar_diffusivities_kg_per_m_s.clear();
  expect_create_rejected(fields, invalid_spec);
  invalid_spec = spec;
  invalid_spec.scalar_diffusivities_kg_per_m_s.front() = -1.0;
  expect_create_rejected(fields, invalid_spec);
  auto invalid_order = fields;
  std::swap(invalid_order.transported_cell_fields[0],
            invalid_order.transported_cell_fields[1]);
  expect_create_rejected(invalid_order, spec);
  const auto stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);

  state.begin_attempt();
  auto final_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto report = transport.finalize_trial(state, final_flux, stencil);
  HUNDUN_CHECK(report.disposition ==
               hundun::flow::MaterialTransportDisposition::finalized);
  HUNDUN_CHECK(report.reason ==
               hundun::flow::MaterialTransportFailureReason::none);
  HUNDUN_CHECK(report.lowest_failing_rank == -1);
  HUNDUN_CHECK(report.flux_provenance ==
               hundun::flow::MaterialFluxProvenance::final_corrected);
  HUNDUN_CHECK(report.shared_face_mass_flux_field == fields.face_mass_flux);
  HUNDUN_CHECK(report.minimum_density_kg_per_m3 > 0.0);
  const auto staged = state.snapshot(hundun::flow::FlowLayer::trial);
  for (double value : staged.density)
    HUNDUN_CHECK(value == 1.0);
  for (double value : staged.transported_cell_fields[0])
    HUNDUN_CHECK(value == 3.0);
  for (double value : staged.transported_cell_fields[1])
    HUNDUN_CHECK(value == 0.4);
  for (double value : staged.face_mass_flux)
    HUNDUN_CHECK(value == 0.0);
  state.rollback_attempt();

  auto structural_over_numerical = initial;
  if (mpi.rank() == 0)
    structural_over_numerical.density.front() = -1.0;
  state.seed_accepted_layers(structural_over_numerical,
                             structural_over_numerical);
  const auto structural_before = state.snapshot(hundun::flow::FlowLayer::trial);
  state.begin_attempt();
  auto provisional_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::provisional);
  const auto rejected =
      transport.finalize_trial(state, provisional_flux, stencil);
  HUNDUN_CHECK(
      rejected.disposition ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(rejected.reason ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  HUNDUN_CHECK(rejected.lowest_failing_rank == 0);
  state.rollback_attempt();
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      structural_before, state.snapshot(hundun::flow::FlowLayer::trial)));

  auto invalid = initial;
  if (mpi.rank() == 0)
    invalid.density.front() = -1.0;
  state.seed_accepted_layers(invalid, invalid);
  const auto history_before = state.snapshot(hundun::flow::FlowLayer::history);
  const auto committed_before =
      state.snapshot(hundun::flow::FlowLayer::committed);
  const auto trial_before = state.snapshot(hundun::flow::FlowLayer::trial);
  const auto metadata_before = state.metadata();
  state.begin_attempt();
  auto positive_failure_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto positive_failure =
      transport.finalize_trial(state, positive_failure_flux, stencil);
  HUNDUN_CHECK(positive_failure.disposition ==
               hundun::flow::MaterialTransportDisposition::recoverable_failure);
  HUNDUN_CHECK(
      positive_failure.reason ==
      hundun::flow::MaterialTransportFailureReason::non_positive_density);
  HUNDUN_CHECK(positive_failure.lowest_failing_rank == 0);
  state.rollback_attempt();
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      history_before, state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      committed_before, state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      trial_before, state.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      metadata_before, state.metadata()));

  state.seed_accepted_layers(initial, initial);
  const auto finite_trial_before =
      state.snapshot(hundun::flow::FlowLayer::trial);
  state.begin_attempt();
  auto nonfinite_write = supplied_flux.acquire_face_write<double>(
      flux_access, phase, actor, fields.face_mass_flux);
  if (mpi.rank() == 0)
    nonfinite_write(0U, 0) = std::numeric_limits<double>::infinity();
  if (mpi.rank() == 1)
    nonfinite_write(0U, 0) = 1.0e6;
  auto nonfinite_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto nonfinite =
      transport.finalize_trial(state, nonfinite_flux, stencil);
  HUNDUN_CHECK(nonfinite.disposition ==
               hundun::flow::MaterialTransportDisposition::recoverable_failure);
  HUNDUN_CHECK(nonfinite.reason ==
               hundun::flow::MaterialTransportFailureReason::non_finite_state);
  HUNDUN_CHECK(nonfinite.lowest_failing_rank == 0);
  state.rollback_attempt();
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      finite_trial_before, state.snapshot(hundun::flow::FlowLayer::trial)));

  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    nonfinite_write(face_index, 0) = 0.0;
  using TestAccess = hundun::flow::test::MaterialDensityTransportTestAccess;
  const auto check_post_stage_gate = [&](auto arm, auto expected_reason) {
    TestAccess::reset();
    state.seed_accepted_layers(initial, initial);
    const auto history = state.snapshot(hundun::flow::FlowLayer::history);
    const auto committed = state.snapshot(hundun::flow::FlowLayer::committed);
    const auto metadata = state.metadata();
    state.begin_attempt();
    arm();
    auto gate_flux = hundun::flow::MaterialFaceMassFlux::acquire(
        registry, supplied_flux, flux_access, phase, actor,
        fields.face_mass_flux, topology,
        hundun::flow::MaterialFluxProvenance::final_corrected);
    const auto gate = transport.finalize_trial(state, gate_flux, stencil);
    HUNDUN_CHECK(
        gate.disposition ==
        hundun::flow::MaterialTransportDisposition::recoverable_failure);
    HUNDUN_CHECK(gate.reason == expected_reason);
    HUNDUN_CHECK(gate.lowest_failing_rank == 0);
    const auto dirty = state.snapshot(hundun::flow::FlowLayer::trial);
    HUNDUN_CHECK(dirty.face_mass_flux.front() == 0.0);
    HUNDUN_CHECK(
        !hundun::test::flow_layer_values_bitwise_equal(dirty, initial));
    state.rollback_attempt();
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        history, state.snapshot(hundun::flow::FlowLayer::history)));
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        committed, state.snapshot(hundun::flow::FlowLayer::committed)));
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        committed, state.snapshot(hundun::flow::FlowLayer::trial)));
    HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
        metadata, state.metadata()));
  };
  check_post_stage_gate(
      [&] { TestAccess::set_density_residual(2.0e-10); },
      hundun::flow::MaterialTransportFailureReason::final_density_residual);
  check_post_stage_gate(
      [&] { TestAccess::set_transport_residual(0U, 2.0e-9); },
      hundun::flow::MaterialTransportFailureReason::final_transport_residual);
  check_post_stage_gate(
      [&] { TestAccess::set_mass_conservation_defect(1.0e-10); },
      hundun::flow::MaterialTransportFailureReason::final_conservation_defect);
  check_post_stage_gate(
      [&] { TestAccess::set_transport_conservation_defect(0U, 1.0e-10); },
      hundun::flow::MaterialTransportFailureReason::final_conservation_defect);
  if (mpi.size() > 1) {
    check_post_stage_gate(
        [&] {
          TestAccess::set_density_residual(2.0e-10, 1);
          TestAccess::set_transport_residual(0U, 2.0e-9, 0);
        },
        hundun::flow::MaterialTransportFailureReason::final_transport_residual);
  }

  state.seed_accepted_layers(initial, initial);
  state.begin_attempt();
  hundun::runtime::FieldStorage stale_storage(
      registry, {local, topology.local_face_count()});
  auto stale_write = stale_storage.acquire_face_write<double>(
      flux_access, phase, actor, fields.face_mass_flux);
  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    stale_write(face_index, 0) = 0.0;
  auto stale_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, stale_storage, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  stale_storage.begin_rebuild();
  const auto stale_report =
      transport.finalize_trial(state, stale_flux, stencil);
  HUNDUN_CHECK(
      stale_report.disposition ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(stale_report.reason ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  HUNDUN_CHECK(stale_report.lowest_failing_rank == 0);
  state.rollback_attempt();

  auto reordered_fields = fields;
  std::swap(reordered_fields.transported_cell_fields[0],
            reordered_fields.transported_cell_fields[1]);
  auto reordered_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, reordered_fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  reordered_state.seed_accepted_layers(initial, initial);
  reordered_state.begin_attempt();
  auto reordered_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto reordered_report =
      transport.finalize_trial(reordered_state, reordered_flux, stencil);
  HUNDUN_CHECK(
      reordered_report.disposition ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(reordered_report.reason ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  reordered_state.rollback_attempt();

  hundun::runtime::FieldRegistry alternate_registry;
  const auto alternate_density =
      alternate_registry.declare_field(cell("rho", "kg/m3"));
  static_cast<void>(alternate_density);
  static_cast<void>(
      alternate_registry.declare_field(cell_components("velocity", 3U)));
  static_cast<void>(alternate_registry.declare_field(cell("pi", "Pa")));
  static_cast<void>(alternate_registry.declare_field(face("u_face", 3U)));
  const auto alternate_flux_id =
      hundun::finite_volume::declare_face_mass_flux(alternate_registry);
  static_cast<void>(alternate_registry.declare_field(cell("rho_h", "J/m3")));
  static_cast<void>(
      alternate_registry.declare_field(cell("rho_phi_alpha", "kg/m3")));
  alternate_registry.freeze();
  hundun::runtime::FieldStorage alternate_storage(
      alternate_registry, {local, topology.local_face_count()});
  hundun::runtime::FieldAccessPlan alternate_access(alternate_registry);
  alternate_access.declare_access(phase, actor, alternate_flux_id,
                                  hundun::runtime::AccessMode::read_write);
  alternate_access.freeze();
  auto alternate_write = alternate_storage.acquire_face_write<double>(
      alternate_access, phase, actor, alternate_flux_id);
  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    alternate_write(face_index, 0) = 0.0;
  state.seed_accepted_layers(initial, initial);
  state.begin_attempt();
  auto wrong_registry_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      alternate_registry, alternate_storage, alternate_access, phase, actor,
      alternate_flux_id, topology,
      hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto wrong_registry_report =
      transport.finalize_trial(state, wrong_registry_flux, stencil);
  HUNDUN_CHECK(
      wrong_registry_report.disposition ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(wrong_registry_report.reason ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  state.rollback_attempt();

  state.seed_accepted_layers(initial, initial);
  state.begin_attempt();
  const auto before_finalization_wrap =
      state.snapshot(hundun::flow::FlowLayer::trial);
  TestAccess::force_finalization_identity_wrap(transport);
  bool wrap_rejected = false;
  try {
    auto wrap_flux = hundun::flow::MaterialFaceMassFlux::acquire(
        registry, supplied_flux, flux_access, phase, actor,
        fields.face_mass_flux, topology,
        hundun::flow::MaterialFluxProvenance::final_corrected);
    static_cast<void>(transport.finalize_trial(state, wrap_flux, stencil));
  } catch (const hundun::runtime::Error &) {
    wrap_rejected = true;
  }
  HUNDUN_CHECK(wrap_rejected);
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      before_finalization_wrap,
      state.snapshot(hundun::flow::FlowLayer::trial)));
  state.rollback_attempt();

  const auto before_attempt_wrap =
      state.snapshot(hundun::flow::FlowLayer::trial);
  TestAccess::force_attempt_identity_wrap(state);
  wrap_rejected = false;
  try {
    state.begin_attempt();
  } catch (const hundun::runtime::Error &) {
    wrap_rejected = true;
  }
  HUNDUN_CHECK(wrap_rejected);
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      before_attempt_wrap, state.snapshot(hundun::flow::FlowLayer::trial)));
}

void run_open_boundary_case(const hundun::runtime::MpiContext &mpi,
                            double enthalpy_diffusivity,
                            double scalar_diffusivity) {
  constexpr hundun::runtime::Int3 extent{8, 4, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {false, false, false},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(open_case(), topology);
  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", "kg/m3"));
  fields.velocity = registry.declare_field(cell_components("velocity", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("pi", "Pa"));
  fields.face_velocity = registry.declare_field(face("u_face", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3"));
  const auto rho_phi = registry.declare_field(cell("rho_phi", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_phi};
  registry.freeze();
  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  constexpr double dt = 0.01;
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, dt, dt, hundun::flow::MomentumTimeOrder::bdf2});
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues current;
  history.density.assign(topology.owned_cell_count(), 1.0);
  current.density = history.density;
  history.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  current.velocity = history.velocity;
  history.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  current.mechanical_pressure = history.mechanical_pressure;
  history.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  current.face_velocity = history.face_velocity;
  history.face_mass_flux.assign(topology.local_face_count(), 17.0);
  current.face_mass_flux = history.face_mass_flux;
  history.transported_cell_fields = {
      std::vector<double>(topology.owned_cell_count(), 2.0),
      std::vector<double>(topology.owned_cell_count(), 0.3)};
  current.transported_cell_fields = {
      std::vector<double>(topology.owned_cell_count(), 3.0),
      std::vector<double>(topology.owned_cell_count(), 0.4)};
  state.seed_accepted_layers(history, current);
  hundun::runtime::FieldStorage flux_storage(
      registry, {local, topology.local_face_count()});
  constexpr hundun::runtime::PhaseId phase = 1912U;
  constexpr hundun::runtime::ActorId actor = 1912U;
  hundun::runtime::FieldAccessPlan access(registry);
  access.declare_access(phase, actor, fields.face_mass_flux,
                        hundun::runtime::AccessMode::read_write);
  access.freeze();
  auto writer = flux_storage.acquire_face_write<double>(access, phase, actor,
                                                        fields.face_mass_flux);
  for (hundun::mesh::LocalFaceId face_id = 0;
       face_id < topology.local_face_count(); ++face_id) {
    writer(face_id, 0) =
        geometry.face_area_vector_m2(face_id, hundun::mesh::FaceSide::owner).x;
  }
  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::flow::MaterialDensityTransportSpec spec;
  spec.enthalpy_density = rho_h;
  spec.enthalpy_diffusivity_kg_per_m_s = enthalpy_diffusivity;
  spec.scalar_densities = {rho_phi};
  spec.scalar_diffusivities_kg_per_m_s = {scalar_diffusivity};
  auto transport = hundun::flow::MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, halo,
      fields, spec);
  state.begin_attempt();
  auto final_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::bdf2, dt, dt);
  const auto report = transport.finalize_trial(state, final_flux, stencil);
  HUNDUN_CHECK(report.disposition ==
               hundun::flow::MaterialTransportDisposition::finalized);
  HUNDUN_CHECK(report.mass_relative_conservation_defect <= 5.0e-11);
  for (double value : report.transport_relative_conservation_defect)
    HUNDUN_CHECK(value <= 5.0e-11);
  const auto trial = state.snapshot(hundun::flow::FlowLayer::trial);
  double local_h{};
  double local_phi{};
  for (hundun::mesh::LocalCellId cell_id = 0;
       cell_id < topology.owned_cell_count(); ++cell_id) {
    const double volume = geometry.cell_volume_m3(cell_id);
    local_h += volume * trial.transported_cell_fields[0][cell_id];
    local_phi += volume * trial.transported_cell_fields[1][cell_id];
  }
  double totals[2]{local_h, local_phi};
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  constexpr double convective_h = 3.3133333333333333333;
  constexpr double convective_phi = 0.4313333333333333333;
  constexpr double diffusive_h = 3.2973333333333333333;
  constexpr double diffusive_phi = 0.4306933333333333333;
  if (enthalpy_diffusivity == 0.0 && scalar_diffusivity == 0.0) {
    HUNDUN_CHECK(std::abs(totals[0] - convective_h) <= 2.0e-13);
    HUNDUN_CHECK(std::abs(totals[1] - convective_phi) <= 2.0e-13);
  } else {
    HUNDUN_CHECK(std::abs(totals[0] - diffusive_h) <= 2.0e-13);
    HUNDUN_CHECK(std::abs(totals[1] - diffusive_phi) <= 2.0e-13);
  }
  state.rollback_attempt();
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    run(mpi);
    run_open_boundary_case(mpi, 0.0, 0.0);
    run_open_boundary_case(mpi, 0.05, 0.02);
  });
}
