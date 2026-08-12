// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_material_density_transport_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_state.hpp"
#include "hundun/flow_material_density_transport.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_mpi_operation_error.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
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
  config.scalars.push_back({"beta", 0.0});
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
  config.scalars.pop_back();
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

struct ExactFlowState final {
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues committed;
  hundun::flow::FlowLayerValues trial;
  hundun::flow::AcceptedStepMetadata metadata;
};

ExactFlowState capture_exact_state(const hundun::flow::FlowState &state) {
  return {state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata()};
}

void require_exact_state(const ExactFlowState &expected,
                         const hundun::flow::FlowState &state) {
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.history, state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.committed, state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.trial, state.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      expected.metadata, state.metadata()));
}

bool provenance_matches(
    const hundun::flow::MaterialDensityTransportReport &report,
    double expected_value, hundun::mesh::GlobalCellId expected_global_cell,
    int expected_rank) {
  return report.minimum_density_available() &&
         hundun::test::fp64_bits(report.minimum_density_kg_per_m3()) ==
             hundun::test::fp64_bits(expected_value) &&
         report.minimum_density_global_cell() == expected_global_cell &&
         report.minimum_density_rank() == expected_rank;
}

void run(const hundun::runtime::MpiContext &mpi) {
  using TestAccess = hundun::flow::test::MaterialDensityTransportTestAccess;
  using ScaleInput = hundun::flow::test::MaterialConservationScaleInput;
  const double epsilon = std::numeric_limits<double>::epsilon();
  HUNDUN_CHECK(TestAccess::conservation_denominator(ScaleInput{
                   3.0, 4.0, 5.0, -0.5, 2.0, 30.0, 40.0, 50.0}) == 4.0);
  HUNDUN_CHECK(TestAccess::conservation_denominator(ScaleInput{
                   2.0, 2.0, -6.0, -0.25, 0.0, 2.0, 2.0, 6.0}) == 2.0);
  HUNDUN_CHECK(TestAccess::conservation_denominator(ScaleInput{
                   63.0 * epsilon, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}) == 1.0);
  HUNDUN_CHECK(TestAccess::conservation_denominator(ScaleInput{
                   65.0 * epsilon, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}) ==
               65.0 * epsilon);
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
  const auto rho_phi_beta =
      registry.declare_field(cell("rho_phi_beta", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_phi, rho_phi_beta};
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
  initial.transported_cell_fields.resize(3U);
  initial.transported_cell_fields[0].assign(topology.owned_cell_count(), 3.0);
  initial.transported_cell_fields[1].assign(topology.owned_cell_count(), 0.4);
  initial.transported_cell_fields[2].assign(topology.owned_cell_count(), 0.7);
  state.seed_accepted_layers(initial, initial);

  {
    auto rollback_history = initial;
    rollback_history.density.assign(topology.owned_cell_count(), 1.25);
    rollback_history.transported_cell_fields[0].assign(
        topology.owned_cell_count(), 3.75);
    rollback_history.transported_cell_fields[1].assign(
        topology.owned_cell_count(), 0.5);
    rollback_history.transported_cell_fields[2].assign(
        topology.owned_cell_count(), 0.875);
    auto rollback_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    rollback_state.seed_accepted_layers(rollback_history, initial);
    hundun::runtime::FieldAccessPlan rollback_access(registry);
    constexpr hundun::runtime::PhaseId rollback_phase = 1919U;
    constexpr hundun::runtime::ActorId rollback_actor = 1919U;
    rollback_access.declare_access(rollback_phase, rollback_actor,
                                   fields.density,
                                   hundun::runtime::AccessMode::read_write);
    rollback_access.freeze();
    rollback_state.begin_attempt();
    auto first_trial = rollback_state.trial_layer().acquire_write<double>(
        rollback_access, rollback_phase, rollback_actor, fields.density);
    first_trial(0, 0, 0, 0) = 2.0;
    rollback_state.commit_attempt(
        {1U, 0.01, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    const auto arbitrary_inactive_trial =
        rollback_state.snapshot(hundun::flow::FlowLayer::trial);
    auto pre_attempt_view =
        rollback_state.layer(hundun::flow::FlowLayer::trial)
            .acquire_read<double>(rollback_access, rollback_phase,
                                  rollback_actor, fields.density);
    rollback_state.begin_attempt();
    bool stale = false;
    try {
      static_cast<void>(pre_attempt_view(0, 0, 0, 0));
    } catch (const hundun::runtime::Error &) {
      stale = true;
    }
    HUNDUN_CHECK(stale);
    auto active_view = rollback_state.trial_layer().acquire_write<double>(
        rollback_access, rollback_phase, rollback_actor, fields.density);
    active_view(0, 0, 0, 0) = 7.0;
    rollback_state.rollback_attempt();
    stale = false;
    try {
      static_cast<void>(active_view(0, 0, 0, 0));
    } catch (const hundun::runtime::Error &) {
      stale = true;
    }
    HUNDUN_CHECK(stale);
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        arbitrary_inactive_trial,
        rollback_state.snapshot(hundun::flow::FlowLayer::trial)));
  }

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
  spec.scalar_densities = {rho_phi, rho_phi_beta};
  spec.scalar_diffusivities_kg_per_m_s = {0.0, 0.0};
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

  const auto success_before = capture_exact_state(state);
  state.begin_attempt();
  auto final_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  auto moved_final_flux = std::move(final_flux);
  const auto moved_from_report =
      transport.finalize_trial(state, final_flux, stencil);
  HUNDUN_CHECK(
      moved_from_report.disposition() ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(moved_from_report.reason() ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  const auto report =
      transport.finalize_trial(state, moved_final_flux, stencil);
  HUNDUN_CHECK(report.disposition() ==
               hundun::flow::MaterialTransportDisposition::finalized);
  HUNDUN_CHECK(report.reason() ==
               hundun::flow::MaterialTransportFailureReason::none);
  HUNDUN_CHECK(report.lowest_failing_rank() == -1);
  HUNDUN_CHECK(report.flux_provenance() ==
               hundun::flow::MaterialFluxProvenance::final_corrected);
  HUNDUN_CHECK(report.shared_face_mass_flux_field() == fields.face_mass_flux);
  HUNDUN_CHECK(report.minimum_density_kg_per_m3() > 0.0);
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      success_before.history,
      state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      success_before.committed,
      state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      success_before.metadata, state.metadata()));
  const auto staged = state.snapshot(hundun::flow::FlowLayer::trial);
  for (double value : staged.density)
    HUNDUN_CHECK(value == 1.0);
  for (double value : staged.transported_cell_fields[0])
    HUNDUN_CHECK(value == 3.0);
  for (double value : staged.transported_cell_fields[1])
    HUNDUN_CHECK(value == 0.4);
  for (double value : staged.transported_cell_fields[2])
    HUNDUN_CHECK(value == 0.7);
  for (double value : staged.face_mass_flux)
    HUNDUN_CHECK(value == 0.0);
  state.rollback_attempt();
  require_exact_state(success_before, state);

  auto provenance_state = initial;
  provenance_state.density.assign(topology.owned_cell_count(), 2.0);
  constexpr double tied_minimum = 0.5;
  if (mpi.size() == 1) {
    provenance_state.density.front() = tied_minimum;
    provenance_state.density.back() = tied_minimum;
  } else if (mpi.rank() == 0) {
    provenance_state.density.back() = tied_minimum;
  } else if (mpi.rank() == 1) {
    provenance_state.density.front() = tied_minimum;
  }
  state.seed_accepted_layers(provenance_state, provenance_state);
  const auto provenance_before = capture_exact_state(state);
  state.begin_attempt();
  auto provenance_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto provenance_report =
      transport.finalize_trial(state, provenance_flux, stencil);
  HUNDUN_CHECK(provenance_report.disposition() ==
               hundun::flow::MaterialTransportDisposition::finalized);
  HUNDUN_CHECK(provenance_report.reason() ==
               hundun::flow::MaterialTransportFailureReason::none);
  const auto expected_minimum_id = static_cast<hundun::mesh::GlobalCellId>(
      mpi.size() == 1 ? 0 : extent.x / mpi.size());
  const int expected_minimum_rank = mpi.size() == 1 ? 0 : 1;
  HUNDUN_CHECK(provenance_matches(provenance_report, tied_minimum,
                                  expected_minimum_id, expected_minimum_rank));
  const auto rank_first_id = static_cast<hundun::mesh::GlobalCellId>(
      ((extent.z - 1) * extent.y + (extent.y - 1)) * extent.x +
      (extent.x / mpi.size() - 1));
  HUNDUN_CHECK(
      !provenance_matches(provenance_report, tied_minimum, rank_first_id, 0));
  HUNDUN_CHECK(provenance_report.density_residual_available());
  HUNDUN_CHECK(provenance_report.density_normalized_l2() <= 1.0e-10);
  HUNDUN_CHECK(provenance_report.mass_conservation_available());
  HUNDUN_CHECK(provenance_report.mass_relative_conservation_defect() <=
               5.0e-11);
  HUNDUN_CHECK(
      std::all_of(provenance_report.transport_residual_availability().begin(),
                  provenance_report.transport_residual_availability().end(),
                  [](std::uint8_t available) { return available != 0U; }));
  HUNDUN_CHECK(std::all_of(
      provenance_report.transport_conservation_availability().begin(),
      provenance_report.transport_conservation_availability().end(),
      [](std::uint8_t available) { return available != 0U; }));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      provenance_before.history,
      state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      provenance_before.committed,
      state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      provenance_before.metadata, state.metadata()));
  state.rollback_attempt();
  require_exact_state(provenance_before, state);

  auto structural_over_numerical = initial;
  if (mpi.rank() == 0)
    structural_over_numerical.density.front() = -1.0;
  state.seed_accepted_layers(structural_over_numerical,
                             structural_over_numerical);
  const auto structural_before = capture_exact_state(state);
  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    flux_values(face_index, 0) = 5.0;
  state.begin_attempt();
  auto provisional_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::provisional);
  const auto rejected =
      transport.finalize_trial(state, provisional_flux, stencil);
  HUNDUN_CHECK(
      rejected.disposition() ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(rejected.reason() ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  HUNDUN_CHECK(rejected.lowest_failing_rank() == 0);
  state.rollback_attempt();
  require_exact_state(structural_before, state);

  TestAccess::reset();
  TestAccess::set_density_residual(2.0e-10);
  state.seed_accepted_layers(initial, initial);
  const auto predictor_before = capture_exact_state(state);
  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    flux_values(face_index, 0) = 7.0;
  state.begin_attempt();
  auto early_rejected_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::predictor);
  const auto early_rejected =
      transport.finalize_trial(state, early_rejected_flux, stencil);
  HUNDUN_CHECK(early_rejected.reason() ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  state.rollback_attempt();
  require_exact_state(predictor_before, state);
  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    flux_values(face_index, 0) = 0.0;
  const auto final_after_provenance_before = capture_exact_state(state);
  state.begin_attempt();
  auto clean_after_early_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto clean_after_early =
      transport.finalize_trial(state, clean_after_early_flux, stencil);
  HUNDUN_CHECK(clean_after_early.disposition() ==
               hundun::flow::MaterialTransportDisposition::finalized);
  for (double value :
       state.snapshot(hundun::flow::FlowLayer::trial).face_mass_flux)
    HUNDUN_CHECK(hundun::test::fp64_bits(value) ==
                 hundun::test::fp64_bits(0.0));
  state.rollback_attempt();
  require_exact_state(final_after_provenance_before, state);

  auto invalid = initial;
  if (mpi.rank() == 0)
    invalid.density.front() = -1.0;
  state.seed_accepted_layers(invalid, invalid);
  const auto positive_before = capture_exact_state(state);
  state.begin_attempt();
  auto positive_failure_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto positive_failure =
      transport.finalize_trial(state, positive_failure_flux, stencil);
  HUNDUN_CHECK(positive_failure.disposition() ==
               hundun::flow::MaterialTransportDisposition::recoverable_failure);
  HUNDUN_CHECK(
      positive_failure.reason() ==
      hundun::flow::MaterialTransportFailureReason::non_positive_density);
  HUNDUN_CHECK(positive_failure.lowest_failing_rank() == 0);
  state.rollback_attempt();
  require_exact_state(positive_before, state);

  if (mpi.size() > 1) {
    const auto check_phase_two_arbitration =
        [&](int density_rank, int transported_rank,
            hundun::flow::MaterialTransportFailureReason expected_reason,
            int expected_rank) {
          auto mixed = initial;
          if (mpi.rank() == density_rank)
            mixed.density.front() = -1.0;
          state.seed_accepted_layers(mixed, mixed);
          auto before = capture_exact_state(state);
          state.begin_attempt();
          if (mpi.rank() == transported_rank) {
            TestAccess::set_accepted_transport_value(
                state, true, 2U, 0U, std::numeric_limits<double>::infinity());
            TestAccess::set_accepted_transport_value(
                state, false, 2U, 0U, std::numeric_limits<double>::infinity());
          }
          before.history = state.snapshot(hundun::flow::FlowLayer::history);
          before.committed = state.snapshot(hundun::flow::FlowLayer::committed);
          auto mixed_flux = hundun::flow::MaterialFaceMassFlux::acquire(
              registry, supplied_flux, flux_access, phase, actor,
              fields.face_mass_flux, topology,
              hundun::flow::MaterialFluxProvenance::final_corrected);
          const auto mixed_report =
              transport.finalize_trial(state, mixed_flux, stencil);
          HUNDUN_CHECK(mixed_report.reason() == expected_reason);
          HUNDUN_CHECK(mixed_report.lowest_failing_rank() == expected_rank);
          state.rollback_attempt();
          require_exact_state(before, state);
        };
    check_phase_two_arbitration(
        0, 1,
        hundun::flow::MaterialTransportFailureReason::non_positive_density, 0);
    check_phase_two_arbitration(
        1, 0, hundun::flow::MaterialTransportFailureReason::non_finite_state,
        0);
    if (mpi.size() == 4)
      check_phase_two_arbitration(
          3, 2, hundun::flow::MaterialTransportFailureReason::non_finite_state,
          2);
  }

  state.seed_accepted_layers(initial, initial);
  const auto finite_before = capture_exact_state(state);
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
  HUNDUN_CHECK(nonfinite.disposition() ==
               hundun::flow::MaterialTransportDisposition::recoverable_failure);
  HUNDUN_CHECK(nonfinite.reason() ==
               hundun::flow::MaterialTransportFailureReason::non_finite_state);
  HUNDUN_CHECK(nonfinite.lowest_failing_rank() == 0);
  state.rollback_attempt();
  require_exact_state(finite_before, state);

  for (std::size_t face_index = 0; face_index < topology.local_face_count();
       ++face_index)
    nonfinite_write(face_index, 0) = 0.0;
  auto finite_overflow_state = initial;
  if (mpi.rank() == 0)
    finite_overflow_state.transported_cell_fields[2].front() =
        std::numeric_limits<double>::max();
  state.seed_accepted_layers(finite_overflow_state, finite_overflow_state);
  const auto overflow_before = capture_exact_state(state);
  state.begin_attempt();
  auto finite_overflow_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto finite_overflow_report =
      transport.finalize_trial(state, finite_overflow_flux, stencil);
  HUNDUN_CHECK(finite_overflow_report.disposition() ==
               hundun::flow::MaterialTransportDisposition::recoverable_failure);
  HUNDUN_CHECK(finite_overflow_report.reason() ==
               hundun::flow::MaterialTransportFailureReason::non_finite_state);
  HUNDUN_CHECK(finite_overflow_report.lowest_failing_rank() == 0);
  state.rollback_attempt();
  require_exact_state(overflow_before, state);

  const auto check_post_stage_gate = [&](auto arm, auto expected_reason,
                                         int expected_rank = 0) {
    TestAccess::reset();
    state.seed_accepted_layers(initial, initial);
    const auto before = capture_exact_state(state);
    state.begin_attempt();
    arm();
    auto gate_flux = hundun::flow::MaterialFaceMassFlux::acquire(
        registry, supplied_flux, flux_access, phase, actor,
        fields.face_mass_flux, topology,
        hundun::flow::MaterialFluxProvenance::final_corrected);
    const auto gate = transport.finalize_trial(state, gate_flux, stencil);
    HUNDUN_CHECK(
        gate.disposition() ==
        hundun::flow::MaterialTransportDisposition::recoverable_failure);
    HUNDUN_CHECK(gate.reason() == expected_reason);
    HUNDUN_CHECK(gate.lowest_failing_rank() == expected_rank);
    const auto dirty = state.snapshot(hundun::flow::FlowLayer::trial);
    HUNDUN_CHECK(dirty.face_mass_flux.front() == 0.0);
    HUNDUN_CHECK(
        !hundun::test::flow_layer_values_bitwise_equal(dirty, initial));
    state.rollback_attempt();
    require_exact_state(before, state);
  };
  check_post_stage_gate(
      [&] { TestAccess::set_density_residual(2.0e-10); },
      hundun::flow::MaterialTransportFailureReason::final_density_residual);
  check_post_stage_gate(
      [&] { TestAccess::set_transport_residual(2U, 2.0e-9); },
      hundun::flow::MaterialTransportFailureReason::final_transport_residual);
  check_post_stage_gate(
      [&] { TestAccess::set_mass_conservation_defect(1.0e-10); },
      hundun::flow::MaterialTransportFailureReason::final_conservation_defect);
  check_post_stage_gate(
      [&] { TestAccess::set_transport_conservation_defect(2U, 1.0e-10); },
      hundun::flow::MaterialTransportFailureReason::final_conservation_defect);
  if (mpi.size() > 1) {
    check_post_stage_gate(
        [&] {
          TestAccess::set_density_residual(2.0e-10, 1);
          TestAccess::set_transport_residual(2U, 2.0e-9, 0);
        },
        hundun::flow::MaterialTransportFailureReason::final_transport_residual);
    check_post_stage_gate(
        [&] {
          TestAccess::set_density_residual(2.0e-10, 0);
          TestAccess::set_transport_conservation_defect(2U, 1.0e-10, 1);
        },
        hundun::flow::MaterialTransportFailureReason::final_density_residual);
    if (mpi.size() == 4)
      check_post_stage_gate(
          [&] {
            TestAccess::set_density_residual(2.0e-10, 3);
            TestAccess::set_transport_residual(2U, 2.0e-9, 2);
          },
          hundun::flow::MaterialTransportFailureReason::
              final_transport_residual,
          2);
  }

  state.seed_accepted_layers(initial, initial);
  const auto stale_before = capture_exact_state(state);
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
      stale_report.disposition() ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(stale_report.reason() ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  HUNDUN_CHECK(stale_report.lowest_failing_rank() == 0);
  state.rollback_attempt();
  require_exact_state(stale_before, state);

  auto reordered_fields = fields;
  std::swap(reordered_fields.transported_cell_fields[0],
            reordered_fields.transported_cell_fields[1]);
  auto reordered_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, reordered_fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  reordered_state.seed_accepted_layers(initial, initial);
  const auto reordered_before = capture_exact_state(reordered_state);
  reordered_state.begin_attempt();
  auto reordered_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      registry, supplied_flux, flux_access, phase, actor, fields.face_mass_flux,
      topology, hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto reordered_report =
      transport.finalize_trial(reordered_state, reordered_flux, stencil);
  HUNDUN_CHECK(
      reordered_report.disposition() ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(reordered_report.reason() ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  reordered_state.rollback_attempt();
  require_exact_state(reordered_before, reordered_state);

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
  static_cast<void>(
      alternate_registry.declare_field(cell("rho_phi_beta", "kg/m3")));
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
  const auto wrong_registry_before = capture_exact_state(state);
  state.begin_attempt();
  auto wrong_registry_flux = hundun::flow::MaterialFaceMassFlux::acquire(
      alternate_registry, alternate_storage, alternate_access, phase, actor,
      alternate_flux_id, topology,
      hundun::flow::MaterialFluxProvenance::final_corrected);
  const auto wrong_registry_report =
      transport.finalize_trial(state, wrong_registry_flux, stencil);
  HUNDUN_CHECK(
      wrong_registry_report.disposition() ==
      hundun::flow::MaterialTransportDisposition::non_retryable_failure);
  HUNDUN_CHECK(wrong_registry_report.reason() ==
               hundun::flow::MaterialTransportFailureReason::invalid_input);
  state.rollback_attempt();
  require_exact_state(wrong_registry_before, state);

  state.seed_accepted_layers(initial, initial);
  const auto finalization_wrap_before = capture_exact_state(state);
  state.begin_attempt();
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
  state.rollback_attempt();
  require_exact_state(finalization_wrap_before, state);

  const auto attempt_wrap_before = capture_exact_state(state);
  TestAccess::force_attempt_identity_wrap(state);
  wrap_rejected = false;
  try {
    state.begin_attempt();
  } catch (const hundun::runtime::Error &) {
    wrap_rejected = true;
  }
  HUNDUN_CHECK(wrap_rejected);
  require_exact_state(attempt_wrap_before, state);
}

void run_nontrivial_be_case(const hundun::runtime::MpiContext &mpi) {
  using namespace hundun;
  constexpr runtime::Int3 extent{8, 4, 4};
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      runtime::DecompositionOptions{grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(), topology);
  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("be_rho", "kg/m3"));
  fields.velocity = registry.declare_field(cell_components("be_velocity", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("be_pi", "Pa"));
  fields.face_velocity = registry.declare_field(face("be_u_face", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("be_rho_h", "J/m3"));
  const auto rho_phi_0 = registry.declare_field(cell("be_rho_phi_0", "kg/m3"));
  const auto rho_phi_1 = registry.declare_field(cell("be_rho_phi_1", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_phi_0, rho_phi_1};
  const auto oracle_q = registry.declare_field(cell("be_oracle_q", "1"));
  const auto oracle_residual =
      registry.declare_field(cell("be_oracle_residual", "1"));
  const auto oracle_face = registry.declare_field(face("be_oracle_face", 1U));
  registry.freeze();
  const auto box = decomposition.owned_box();
  const runtime::Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                            box.end.z - box.begin.z};
  constexpr double dt = 1.0e-3;
  auto state = flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, dt, 0.0, flow::MomentumTimeOrder::backward_euler});
  flow::FlowLayerValues initial;
  initial.transported_cell_fields.resize(3U);
  for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
       ++cell_id) {
    const auto center = geometry.cell_center_m(cell_id);
    const double rho = 1.0 + 0.1 * center.x + 0.03 * center.y;
    const std::array<double, 3> intensive{
        2.0 + 0.4 * center.x - 0.1 * center.y,
        0.2 + 0.3 * center.x + 0.05 * center.z,
        0.7 - 0.2 * center.y + 0.04 * center.z};
    initial.density.push_back(rho);
    initial.velocity.insert(initial.velocity.end(), {0.0, 0.0, 0.0});
    initial.mechanical_pressure.push_back(0.0);
    for (std::size_t field = 0; field < intensive.size(); ++field)
      initial.transported_cell_fields[field].push_back(rho * intensive[field]);
  }
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 91.0);
  state.seed_accepted_layers(initial, initial);

  runtime::FieldStorage workspace(registry,
                                  {local, topology.local_face_count()});
  constexpr runtime::PhaseId flux_phase = 1920U;
  constexpr runtime::ActorId flux_actor = 1920U;
  runtime::FieldAccessPlan flux_access(registry);
  flux_access.declare_access(flux_phase, flux_actor, fields.face_mass_flux,
                             runtime::AccessMode::read_write);
  flux_access.freeze();
  auto flux_write = workspace.acquire_face_write<double>(
      flux_access, flux_phase, flux_actor, fields.face_mass_flux);
  for (mesh::LocalFaceId face_id = 0; face_id < topology.local_face_count();
       ++face_id) {
    const auto area =
        geometry.face_area_vector_m2(face_id, mesh::FaceSide::owner);
    flux_write(face_id, 0) = 0.05 * area.x + 0.02 * area.y;
  }

  runtime::FieldAccessPlan oracle_access(registry);
  constexpr runtime::PhaseId oracle_phase = 1921U;
  constexpr runtime::ActorId oracle_actor = 1921U;
  for (const auto field : {oracle_q, oracle_residual, oracle_face})
    oracle_access.declare_access(oracle_phase, oracle_actor, field,
                                 runtime::AccessMode::read_write);
  oracle_access.freeze();
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(decomposition, local, 2));
  auto operators =
      finite_volume::CellCenteredFvmOperators::create(topology, geometry);
  auto raw_flux = finite_volume::FaceMassFlux::acquire(
      registry, workspace, flux_access, flux_phase, flux_actor,
      fields.face_mass_flux, topology);
  const auto local_index = [&](mesh::LocalCellId cell_id) {
    const auto global = topology.global_cell(cell_id);
    const auto owned = topology.owned_global_box();
    return runtime::Int3{global.x - owned.begin.x, global.y - owned.begin.y,
                         global.z - owned.begin.z};
  };
  const auto residual_for = [&](std::size_t field,
                                finite_volume::FiniteVolumeQuantity quantity,
                                bool divide_by_density) {
    auto q = workspace.acquire_write<double>(oracle_access, oracle_phase,
                                             oracle_actor, oracle_q);
    for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
         ++cell_id) {
      const auto index = local_index(cell_id);
      const double conserved = initial.transported_cell_fields[field][cell_id];
      q(index.x, index.y, index.z, 0) =
          divide_by_density ? conserved / initial.density[cell_id] : conserved;
    }
    halo.exchange(workspace, oracle_q);
    const auto q_read = std::as_const(workspace).acquire_read<double>(
        oracle_access, oracle_phase, oracle_actor, oracle_q);
    auto face_values = workspace.acquire_face_write<double>(
        oracle_access, oracle_phase, oracle_actor, oracle_face);
    operators.reconstruct_transport_faces(quantity, boundaries, raw_flux,
                                          q_read, face_values);
    const auto face_read = std::as_const(workspace).acquire_face_read<double>(
        oracle_access, oracle_phase, oracle_actor, oracle_face);
    auto residual = workspace.acquire_write<double>(
        oracle_access, oracle_phase, oracle_actor, oracle_residual);
    for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
         ++cell_id) {
      const auto index = local_index(cell_id);
      residual(index.x, index.y, index.z, 0) = 0.0;
    }
    operators.accumulate_convective_residual(raw_flux, face_read, residual);
    std::vector<double> result;
    result.reserve(topology.owned_cell_count());
    for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
         ++cell_id) {
      const auto index = local_index(cell_id);
      result.push_back(residual(index.x, index.y, index.z, 0));
    }
    return result;
  };

  std::array<std::vector<double>, 3> expected_residuals;
  std::array<std::vector<double>, 3> wrong_residuals;
  for (std::size_t field = 0; field < 3U; ++field) {
    const auto quantity =
        field == 0U ? finite_volume::FiniteVolumeQuantity::enthalpy()
                    : finite_volume::FiniteVolumeQuantity::scalar(field - 1U);
    expected_residuals[field] = residual_for(field, quantity, true);
    wrong_residuals[field] = residual_for(field, quantity, false);
  }
  auto mass_residual = workspace.acquire_write<double>(
      oracle_access, oracle_phase, oracle_actor, oracle_residual);
  for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
       ++cell_id) {
    const auto index = local_index(cell_id);
    mass_residual(index.x, index.y, index.z, 0) = 0.0;
  }
  operators.accumulate_mass_residual(raw_flux, mass_residual);

  flow::MaterialDensityTransportSpec spec;
  spec.enthalpy_density = rho_h;
  spec.scalar_densities = {rho_phi_0, rho_phi_1};
  spec.scalar_diffusivities_kg_per_m_s = {0.0, 0.0};
  auto transport = flow::MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, halo,
      fields, spec);
  state.begin_attempt();
  auto final_flux = flow::MaterialFaceMassFlux::acquire(
      registry, workspace, flux_access, flux_phase, flux_actor,
      fields.face_mass_flux, topology,
      flow::MaterialFluxProvenance::final_corrected);
  const auto report = transport.finalize_trial(
      state, final_flux,
      flow::make_momentum_time_stencil(flow::MomentumTimeOrder::backward_euler,
                                       dt, 0.0));
  HUNDUN_CHECK(report.disposition() ==
               flow::MaterialTransportDisposition::finalized);
  const auto trial = state.snapshot(flow::FlowLayer::trial);
  double local_wrong_difference{};
  for (mesh::LocalCellId cell_id = 0; cell_id < topology.owned_cell_count();
       ++cell_id) {
    const auto index = local_index(cell_id);
    const double volume = geometry.cell_volume_m3(cell_id);
    const double expected_density =
        initial.density[cell_id] -
        dt * mass_residual(index.x, index.y, index.z, 0) / volume;
    HUNDUN_CHECK(std::abs(trial.density[cell_id] - expected_density) <=
                 64.0 * std::numeric_limits<double>::epsilon() *
                     std::max(1.0, std::abs(expected_density)));
    for (std::size_t field = 0; field < 3U; ++field) {
      const double expected = initial.transported_cell_fields[field][cell_id] -
                              dt * expected_residuals[field][cell_id] / volume;
      const double wrong = initial.transported_cell_fields[field][cell_id] -
                           dt * wrong_residuals[field][cell_id] / volume;
      HUNDUN_CHECK(
          std::abs(trial.transported_cell_fields[field][cell_id] - expected) <=
          128.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, std::abs(expected)));
      local_wrong_difference =
          std::max(local_wrong_difference, std::abs(expected - wrong));
    }
  }
  double global_wrong_difference{};
  runtime::check_mpi_result(
      MPI_Allreduce(&local_wrong_difference, &global_wrong_difference, 1,
                    MPI_DOUBLE, MPI_MAX, mpi.comm()),
      "MPI_Allreduce(nontrivial BE wrong reconstruction oracle)");
  HUNDUN_CHECK(global_wrong_difference > 1.0e-8);
  state.rollback_attempt();
}

void run_descriptor_rejections(const hundun::runtime::MpiContext &mpi) {
  using namespace hundun;
  constexpr runtime::Int3 extent{8, 4, 4};
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      runtime::DecompositionOptions{grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(), topology);
  const auto box = decomposition.owned_box();
  const runtime::Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                            box.end.z - box.begin.z};
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(decomposition, local, 2));
  const auto expect_rejected = [&](runtime::FieldDescriptor density_descriptor,
                                   runtime::FieldDescriptor enthalpy_descriptor,
                                   runtime::FieldDescriptor scalar_descriptor,
                                   bool duplicate, bool wrong_flux) {
    runtime::FieldRegistry registry;
    flow::FlowFieldIds fields;
    fields.density = registry.declare_field(std::move(density_descriptor));
    fields.velocity =
        registry.declare_field(cell_components("descriptor_velocity", 3U));
    fields.mechanical_pressure =
        registry.declare_field(cell("descriptor_pi", "Pa"));
    fields.face_velocity =
        registry.declare_field(face("descriptor_u_face", 3U));
    fields.face_mass_flux =
        wrong_flux ? registry.declare_field(face("descriptor_bad_flux", 1U))
                   : finite_volume::declare_face_mass_flux(registry);
    const auto enthalpy =
        registry.declare_field(std::move(enthalpy_descriptor));
    const auto scalar0 = registry.declare_field(std::move(scalar_descriptor));
    const auto scalar1 =
        registry.declare_field(cell("descriptor_scalar_1", "kg/m3"));
    flow::MaterialDensityTransportSpec spec;
    spec.enthalpy_density = enthalpy;
    if (duplicate) {
      spec.scalar_densities = {enthalpy, scalar0};
      fields.transported_cell_fields = {enthalpy, enthalpy, scalar0};
    } else {
      spec.scalar_densities = {scalar0, scalar1};
      fields.transported_cell_fields = {enthalpy, scalar0, scalar1};
    }
    spec.scalar_diffusivities_kg_per_m_s = {0.0, 0.0};
    registry.freeze();
    bool rejected = false;
    try {
      static_cast<void>(flow::MaterialDensityTransport::create(
          registry, decomposition, topology, geometry, boundaries, mpi, halo,
          fields, spec));
    } catch (const runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  };
  const auto valid_density = cell("descriptor_rho", "kg/m3");
  const auto valid_enthalpy = cell("descriptor_rho_h", "J/m3");
  const auto valid_scalar = cell("descriptor_scalar_0", "kg/m3");
  {
    auto bad = valid_enthalpy;
    bad.unit = "J/m^3";
    expect_rejected(valid_density, bad, valid_scalar, false, false);
  }
  {
    auto bad = valid_enthalpy;
    bad.conservative = false;
    expect_rejected(valid_density, bad, valid_scalar, false, false);
  }
  {
    auto bad = valid_enthalpy;
    bad.ghost_width = 1;
    expect_rejected(valid_density, bad, valid_scalar, false, false);
  }
  {
    auto bad = valid_enthalpy;
    bad.components = 2U;
    expect_rejected(valid_density, bad, valid_scalar, false, false);
  }
  {
    auto bad = valid_enthalpy;
    bad.scalar_type = runtime::ScalarType::int32;
    expect_rejected(valid_density, bad, valid_scalar, false, false);
  }
  {
    auto bad = valid_enthalpy;
    bad.space = runtime::FunctionSpace::face_value;
    expect_rejected(valid_density, bad, valid_scalar, false, false);
  }
  {
    auto bad = valid_scalar;
    bad.unit = "1";
    expect_rejected(valid_density, valid_enthalpy, bad, false, false);
  }
  {
    auto bad = valid_density;
    bad.restart = runtime::RestartPolicy::transient;
    expect_rejected(bad, valid_enthalpy, valid_scalar, false, false);
  }
  expect_rejected(valid_density, valid_enthalpy, valid_scalar, true, false);
  expect_rejected(valid_density, valid_enthalpy, valid_scalar, false, true);
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
  const auto oracle_q = registry.declare_field(cell("boundary_oracle_q", "1"));
  const auto oracle_gradient =
      registry.declare_field(cell_components("boundary_oracle_gradient", 3U));
  const auto oracle_face =
      registry.declare_field(face("boundary_oracle_face", 1U));
  const auto oracle_gamma =
      registry.declare_field(face("boundary_oracle_gamma", 1U));
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
  for (const auto field :
       {oracle_q, oracle_gradient, oracle_face, oracle_gamma})
    access.declare_access(phase, actor, field,
                          hundun::runtime::AccessMode::read_write);
  access.freeze();
  auto writer = flux_storage.acquire_face_write<double>(access, phase, actor,
                                                        fields.face_mass_flux);
  for (hundun::mesh::LocalFaceId face_id = 0;
       face_id < topology.local_face_count(); ++face_id) {
    writer(face_id, 0) =
        geometry.face_area_vector_m2(face_id, hundun::mesh::FaceSide::owner).x;
  }
  double signed_mass_boundary[2]{};
  for (hundun::mesh::LocalFaceId face_id = 0;
       face_id < topology.local_face_count(); ++face_id) {
    const auto patch = topology.patch_id(face_id);
    if (patch && *patch < 2U &&
        topology.cell_ownership(topology.owner(face_id)) ==
            hundun::mesh::EntityOwnership::owned)
      signed_mass_boundary[*patch] += writer(face_id, 0);
  }
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, signed_mass_boundary, 2, MPI_DOUBLE,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(std::abs(signed_mass_boundary[0] + 1.0) <= 2.0e-14);
  HUNDUN_CHECK(std::abs(signed_mass_boundary[1] - 1.0) <= 2.0e-14);
  HUNDUN_CHECK(std::abs(signed_mass_boundary[0] + signed_mass_boundary[1]) <=
               2.0e-14);
  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  auto raw_flux = hundun::finite_volume::FaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology);
  auto operators = hundun::finite_volume::CellCenteredFvmOperators::create(
      topology, geometry);
  const auto local_index = [&](hundun::mesh::LocalCellId cell_id) {
    const auto global = topology.global_cell(cell_id);
    const auto owned = topology.owned_global_box();
    return hundun::runtime::Int3{global.x - owned.begin.x,
                                 global.y - owned.begin.y,
                                 global.z - owned.begin.z};
  };
  struct BoundaryTerms final {
    double convective{};
    double diffusive{};
  };
  const auto boundary_terms = [&](std::size_t field,
                                  const std::vector<double> &density,
                                  const std::vector<double> &conserved,
                                  double diffusivity) {
    auto q = flux_storage.acquire_write<double>(access, phase, actor, oracle_q);
    for (hundun::mesh::LocalCellId cell_id = 0;
         cell_id < topology.owned_cell_count(); ++cell_id) {
      const auto index = local_index(cell_id);
      q(index.x, index.y, index.z, 0) = conserved[cell_id] / density[cell_id];
    }
    halo.exchange(flux_storage, oracle_q);
    const auto q_read =
        std::as_const(flux_storage)
            .acquire_read<double>(access, phase, actor, oracle_q);
    auto gradient = flux_storage.acquire_write<double>(access, phase, actor,
                                                       oracle_gradient);
    const auto quantity =
        field == 0U ? hundun::finite_volume::FiniteVolumeQuantity::enthalpy()
                    : hundun::finite_volume::FiniteVolumeQuantity::scalar(0U);
    operators.compute_gradient(
        hundun::finite_volume::GradientScheme::green_gauss, quantity,
        boundaries, q_read, gradient);
    halo.exchange(flux_storage, oracle_gradient);
    const auto gradient_read =
        std::as_const(flux_storage)
            .acquire_read<double>(access, phase, actor, oracle_gradient);
    auto face_values = flux_storage.acquire_face_write<double>(
        access, phase, actor, oracle_face);
    operators.reconstruct_transport_faces(quantity, boundaries, raw_flux,
                                          q_read, face_values);
    auto gamma = flux_storage.acquire_face_write<double>(access, phase, actor,
                                                         oracle_gamma);
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id)
      gamma(face_id, 0) = diffusivity;
    const auto face_read =
        std::as_const(flux_storage)
            .acquire_face_read<double>(access, phase, actor, oracle_face);
    const auto gamma_read =
        std::as_const(flux_storage)
            .acquire_face_read<double>(access, phase, actor, oracle_gamma);
    std::vector<hundun::finite_volume::PhysicalBoundaryTransportContribution>
        contributions;
    operators.physical_boundary_transport_contributions(
        quantity, boundaries, raw_flux, face_read, q_read, gradient_read,
        gamma_read, contributions);
    BoundaryTerms local_terms;
    for (const auto &entry : contributions) {
      local_terms.convective += entry.convective;
      local_terms.diffusive += entry.diffusive;
    }
    double global_terms[2]{local_terms.convective, local_terms.diffusive};
    hundun::runtime::check_mpi_result(
        MPI_Allreduce(MPI_IN_PLACE, global_terms, 2, MPI_DOUBLE, MPI_SUM,
                      mpi.comm()),
        "MPI_Allreduce(material boundary term oracle)");
    return BoundaryTerms{global_terms[0], global_terms[1]};
  };
  const std::array<BoundaryTerms, 2> current_terms{
      boundary_terms(0U, current.density, current.transported_cell_fields[0],
                     enthalpy_diffusivity),
      boundary_terms(1U, current.density, current.transported_cell_fields[1],
                     scalar_diffusivity)};
  const std::array<BoundaryTerms, 2> history_terms{
      boundary_terms(0U, history.density, history.transported_cell_fields[0],
                     enthalpy_diffusivity),
      boundary_terms(1U, history.density, history.transported_cell_fields[1],
                     scalar_diffusivity)};
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
  HUNDUN_CHECK(report.disposition() ==
               hundun::flow::MaterialTransportDisposition::finalized);
  HUNDUN_CHECK(report.mass_relative_conservation_defect() <= 5.0e-11);
  for (double value : report.transport_relative_conservation_defect())
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
  constexpr std::array<double, 2> current_integrals{3.0, 0.4};
  constexpr std::array<double, 2> history_integrals{2.0, 0.3};
  for (std::size_t field = 0; field < 2U; ++field) {
    const double current_boundary =
        current_terms[field].convective + current_terms[field].diffusive;
    const double history_boundary =
        history_terms[field].convective + history_terms[field].diffusive;
    const double effective_boundary = 2.0 * current_boundary - history_boundary;
    const double expected =
        current_integrals[field] -
        (1.0 / 3.0) * (history_integrals[field] - current_integrals[field]) -
        (dt / 1.5) * effective_boundary;
    HUNDUN_CHECK(std::abs(totals[field] - expected) <= 2.0e-13);
    HUNDUN_CHECK(std::abs(current_terms[field].convective) > 1.0e-12);
    HUNDUN_CHECK(std::abs(history_terms[field].convective) > 1.0e-12);
    const double omit_history =
        current_integrals[field] -
        (1.0 / 3.0) * (history_integrals[field] - current_integrals[field]) -
        (dt / 1.5) * (2.0 * current_boundary);
    const double omit_factor_two =
        current_integrals[field] -
        (1.0 / 3.0) * (history_integrals[field] - current_integrals[field]) -
        (dt / 1.5) * (current_boundary - history_boundary);
    const double reversed_sign =
        current_integrals[field] -
        (1.0 / 3.0) * (history_integrals[field] - current_integrals[field]) +
        (dt / 1.5) * effective_boundary;
    HUNDUN_CHECK(std::abs(totals[field] - omit_history) > 1.0e-6);
    HUNDUN_CHECK(std::abs(totals[field] - omit_factor_two) > 1.0e-6);
    HUNDUN_CHECK(std::abs(totals[field] - reversed_sign) > 1.0e-6);
    if (enthalpy_diffusivity != 0.0 && scalar_diffusivity != 0.0) {
      HUNDUN_CHECK(std::abs(current_terms[field].diffusive) > 1.0e-12);
      HUNDUN_CHECK(std::abs(history_terms[field].diffusive) > 1.0e-12);
      const double omit_diffusion =
          current_integrals[field] -
          (1.0 / 3.0) * (history_integrals[field] - current_integrals[field]) -
          (dt / 1.5) * (2.0 * current_terms[field].convective -
                        history_terms[field].convective);
      HUNDUN_CHECK(std::abs(totals[field] - omit_diffusion) > 1.0e-6);
    }
  }
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
    run_nontrivial_be_case(mpi);
    run_descriptor_rejections(mpi);
    run_open_boundary_case(mpi, 0.0, 0.0);
    run_open_boundary_case(mpi, 0.05, 0.02);
  });
}
