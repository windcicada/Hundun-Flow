// SPDX-License-Identifier: Apache-2.0

#include "flow/src/material_density_piso_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/material_density_piso.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
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

namespace {

hundun::runtime::Int3 grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task 20 rank count");
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t patch = 0; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name,
          unit,
          "task20",
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
          "task20",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

struct ExactState final {
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues committed;
  hundun::flow::FlowLayerValues trial;
  hundun::flow::AcceptedStepMetadata metadata;
};

ExactState capture(const hundun::flow::FlowState &state) {
  return {state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata()};
}

void check_equal(const ExactState &expected,
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

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    constexpr hundun::runtime::Int3 extent{8, 4, 4};
    auto decomposition = hundun::runtime::StructuredDecomposition::create(
        mpi, extent, {true, true, true},
        hundun::runtime::DecompositionOptions{grid(mpi.size())});
    hundun::mesh::MeshTopology topology(decomposition);
    hundun::mesh::MeshGeometry geometry(
        topology,
        hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    auto boundaries = hundun::boundary::BoundaryRegistry::create(
        periodic_case(), topology);

    hundun::runtime::FieldRegistry registry;
    hundun::flow::FlowFieldIds fields;
    fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
    fields.velocity =
        registry.declare_field(cell("velocity", "m/s", 3U, false));
    fields.mechanical_pressure =
        registry.declare_field(cell("pi", "Pa", 1U, false));
    fields.face_velocity =
        registry.declare_field(face("face_velocity", "m/s", 3U));
    fields.face_mass_flux =
        hundun::finite_volume::declare_face_mass_flux(registry);
    const auto rho_h =
        registry.declare_field(cell("rho_h", "J/m3", 1U, true));
    fields.transported_cell_fields = {rho_h};
    registry.freeze();

    const auto box = decomposition.owned_box();
    const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                      box.end.y - box.begin.y,
                                      box.end.z - box.begin.z};
    const auto make_state = [&](double density) {
      auto state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      hundun::flow::FlowLayerValues initial;
      initial.density.assign(topology.owned_cell_count(), density);
      initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
      initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
      initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
      initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
      initial.transported_cell_fields = {
          std::vector<double>(topology.owned_cell_count(), 3.0)};
      state.seed_accepted_layers(initial, initial);
      return state;
    };

    auto halo = hundun::runtime::HaloExchange::create(
        decomposition,
        hundun::runtime::ExchangePlan::create(decomposition, local, 2));
    hundun::execution::CpuReferenceContext execution;
    hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
    hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
    hundun::linear::JacobiPreconditioner mx(execution);
    hundun::linear::JacobiPreconditioner my(execution);
    hundun::linear::JacobiPreconditioner mz(execution);
    hundun::linear::JacobiPreconditioner pressure_preconditioner(execution);
    hundun::flow::MaterialDensityTransportSpec specification;
    specification.enthalpy_density = rho_h;
    specification.enthalpy_diffusivity_kg_per_m_s = 0.0;
    auto flow = hundun::flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, specification);
    const auto be = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);

    auto state = make_state(1.0);
    const auto report = flow.attempt(state, 0.0, be, {}, {});
    HUNDUN_CHECK(report.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(report.flow().reason == hundun::flow::StepFailureReason::none);
    HUNDUN_CHECK(report.flow().pressure_corrector_count == 2U);
    HUNDUN_CHECK(report.material_report_available());
    HUNDUN_CHECK(report.material_report().disposition() ==
                 hundun::flow::MaterialTransportDisposition::finalized);
    HUNDUN_CHECK(report.flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(report.material_field_count() == 1U);
    HUNDUN_CHECK(report.final_continuity_residual_available());
    HUNDUN_CHECK(report.final_pressure_residual_available());
    HUNDUN_CHECK(state.metadata().step == 1U);

    auto source = flow.diagnostic_source(state, report);
    HUNDUN_CHECK(source.fingerprint_field_count() == 5U);
    HUNDUN_CHECK(source.fingerprint_field_id(0) == "face_mass_flux");
    HUNDUN_CHECK(source.fingerprint_field_id(4) == "velocity");
    HUNDUN_CHECK(source.committed_step() == 1U);
    HUNDUN_CHECK(source.field_item_count(3) == topology.owned_cell_count());
    HUNDUN_CHECK(source.field_value(3, 0U, 0U) == 1.0);

    auto moved_flow = std::move(flow);
    bool pre_move_source_stale = false;
    try {
      static_cast<void>(source.report());
    } catch (const hundun::runtime::Error &error) {
      pre_move_source_stale = std::string_view(error.what()) ==
                              "material flow diagnostic source is stale";
    }
    HUNDUN_CHECK(pre_move_source_stale);
    auto fresh_source = moved_flow.diagnostic_source(state, report);
    HUNDUN_CHECK(fresh_source.committed_step() == 1U);
    bool moved_from_rejected = false;
    try {
      static_cast<void>(flow.diagnostic_source(state, report));
    } catch (const hundun::runtime::Error &error) {
      moved_from_rejected = std::string_view(error.what()) ==
                            "material flow object has been moved from";
    }
    HUNDUN_CHECK(moved_from_rejected);

    auto invalid_state = make_state(1.0);
    const auto invalid_before = capture(invalid_state);
    const auto invalid =
        moved_flow.attempt(invalid_state, -1.0, be, {}, {});
    HUNDUN_CHECK(invalid.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(invalid.flow().reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    check_equal(invalid_before, invalid_state);
    bool stale = false;
    try {
      static_cast<void>(fresh_source.report());
    } catch (const hundun::runtime::Error &) {
      stale = true;
    }
    HUNDUN_CHECK(stale);

    auto negative_state = make_state(-1.0);
    const auto negative_before = capture(negative_state);
    const auto negative =
        moved_flow.attempt(negative_state, 0.0, be, {}, {});
    HUNDUN_CHECK(negative.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(negative.flow().reason ==
                 hundun::flow::StepFailureReason::transport_failure);
    HUNDUN_CHECK(negative.material_failure_reason() ==
                 hundun::flow::MaterialTransportFailureReason::
                     non_positive_density);
    HUNDUN_CHECK(negative.flow().suggested_dt_s == 0.005);
    check_equal(negative_before, negative_state);

    const auto bdf2 = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::bdf2, 0.01, 0.01);
    const auto second = moved_flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(second.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(state.metadata().step == 2U);

    const auto exhaustion_before = capture(state);
    hundun::flow::test::MaterialDensityPisoTestAccess::
        force_flow_attempt_identity(moved_flow,
                                    std::numeric_limits<std::uint64_t>::max());
    bool exhausted = false;
    try {
      static_cast<void>(moved_flow.attempt(state, 0.0, bdf2, {}, {}));
    } catch (const hundun::runtime::Error &error) {
      exhausted = std::string_view(error.what()) ==
                  "material flow attempt identity would wrap";
    }
    HUNDUN_CHECK(exhausted);
    check_equal(exhaustion_before, state);
  });
}
