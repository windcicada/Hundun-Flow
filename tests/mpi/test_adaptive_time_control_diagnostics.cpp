// SPDX-License-Identifier: Apache-2.0

#include "diagnostics/src/time_control_diagnostics_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/diagnostics/time_control_diagnostics.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <string>

namespace {

struct RecordingSink final : hundun::diagnostics::DiagnosticSink {
  void submit(const hundun::diagnostics::DiagnosticRecord &value) override {
    record = value;
    ++calls;
  }
  hundun::diagnostics::DiagnosticRecord record;
  int calls{};
};

hundun::runtime::Int3 grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task22 rank count");
}

hundun::runtime::FieldDescriptor cell(const char *name,
                                      std::uint32_t components) {
  return {name, "1", "task22", hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, components, 2, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name, "1", "task22", hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64, components, 0, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back({"alpha", 0.0});
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t i = 0; i < names.size(); ++i) {
    config.boundaries[i].patch = names[i];
    config.boundaries[i].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

struct HistoricalSource final {
  hundun::flow::TimeControlDiagnosticSource source;
  std::uint64_t step{};
  double time_s{};
};

HistoricalSource
make_historical_source(const hundun::runtime::MpiContext &mpi) {
  constexpr hundun::runtime::Int3 extent{8, 6, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", 1U));
  fields.velocity = registry.declare_field(cell("u", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("pi", 1U));
  fields.face_velocity = registry.declare_field(face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields.push_back(
      registry.declare_field(cell("alpha", 1U)));
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  const auto cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(cells, 1.0);
  initial.velocity.assign(cells * 3U, 0.0);
  initial.mechanical_pressure.assign(cells, 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields = {std::vector<double>(cells, 1.0)};
  state.seed_accepted_layers(initial, initial);

  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution), my(execution),
      mz(execution), pressure_preconditioner(execution);
  auto facade = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner,
      {{fields.transported_cell_fields.front(),
        hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  const hundun::config::FlowTimeConfig time{
      hundun::config::TimeMode::adaptive, 2, 0.01, 0.00125, 0.02,
      0.5, 0.25, 1.25, 0.5, 8};
  auto controller = hundun::flow::Bdf2RetryController::create(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state);
  const auto first = controller.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(first.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  auto restored = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state, controller.state());
  const auto second = restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(second.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(state.metadata().step == 2U);
  const auto metadata = state.metadata();
  return {restored.diagnostic_source(state, second), metadata.step,
          metadata.time_s};
}

void run(const hundun::runtime::MpiContext &mpi) {
  using DiagnosticAccess =
      hundun::diagnostics::test::TimeControlDiagnosticsTestAccess;
  using Fault = hundun::diagnostics::test::TimeControlDiagnosticFault;
  DiagnosticAccess::reset();

  // The source is owning: all controller, report, state, mesh, and solver
  // owners used to create it have already been destroyed on return.
  auto historical = make_historical_source(mpi);
  auto &source = historical.source;
  const hundun::diagnostics::DiagnosticFrame frame{
      mpi.rank(), historical.step, historical.time_s,
      "time-control.advance-result"};
  const auto counters_before = mpi.fp64_reduction_counters();
  for (const auto level :
       {hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticLevel::invariants,
        hundun::diagnostics::DiagnosticLevel::counters,
        hundun::diagnostics::DiagnosticLevel::bounded_state_sample}) {
    RecordingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        level, hundun::diagnostics::DiagnosticScope::local, frame, {},
        level == hundun::diagnostics::DiagnosticLevel::bounded_state_sample
            ? 256U
            : 0U};
    hundun::diagnostics::collect_diagnostics(source, request, sink);
    HUNDUN_CHECK(sink.calls == 1);
    HUNDUN_CHECK(sink.record.identities.size() == 4U);
    if (level == hundun::diagnostics::DiagnosticLevel::summary)
      HUNDUN_CHECK(sink.record.metrics.size() == 5U);
    if (level == hundun::diagnostics::DiagnosticLevel::invariants)
      HUNDUN_CHECK(sink.record.invariants.size() == 7U);
    if (level == hundun::diagnostics::DiagnosticLevel::counters)
      HUNDUN_CHECK(sink.record.counters.size() == 5U);
    if (level ==
        hundun::diagnostics::DiagnosticLevel::bounded_state_sample) {
      const std::size_t expected = mpi.rank() == 0 ? 9U : 0U;
      HUNDUN_CHECK(sink.record.samples.size() == expected);
      HUNDUN_CHECK(sink.record.eligible_sample_count == expected);
    }
  }
  const auto counters_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(counters_after.collective_calls ==
               counters_before.collective_calls);
  HUNDUN_CHECK(counters_after.reduced_scalars ==
               counters_before.reduced_scalars);
  HUNDUN_CHECK(counters_after.logical_payload_bytes ==
               counters_before.logical_payload_bytes);

  const hundun::diagnostics::DiagnosticRequest collective_request{
      hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
      hundun::diagnostics::DiagnosticScope::collective, frame, {}, 256U};
  RecordingSink collective_sink;
  hundun::diagnostics::collect_diagnostics(source, mpi, collective_request,
                                           collective_sink);
  HUNDUN_CHECK(collective_sink.calls == 1);
  HUNDUN_CHECK(collective_sink.record.samples.size() == 9U);
  HUNDUN_CHECK(collective_sink.record.eligible_sample_count == 9U);

  struct FaultCase final {
    Fault fault;
    const char *code;
  };
  constexpr std::array<FaultCase, 7> faults{{
      {Fault::phase1_layout, "time-control.diagnostics.local-layout"},
      {Fault::phase2_request,
       "time-control.diagnostics.request-preparation"},
      {Fault::phase3_provider,
       "time-control.diagnostics.provider-agreement"},
      {Fault::phase4_payload,
       "time-control.diagnostics.sample-preparation"},
      {Fault::phase4_wire, "time-control.diagnostics.sample-wire"},
      {Fault::phase5_record, "time-control.diagnostics.record"},
      {Fault::phase6_sink, "diagnostics.sink.submit"},
  }};
  const int injection_rank = mpi.size() == 1 ? 0 : 1;
  for (const auto &fault : faults) {
    DiagnosticAccess::reset();
    DiagnosticAccess::set_fault(fault.fault, injection_rank);
    RecordingSink sink;
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, mpi, collective_request,
                                               sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(error.code() == fault.code);
      HUNDUN_CHECK(error.lowest_failing_rank() == injection_rank);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(DiagnosticAccess::submission_count() ==
                 (fault.fault == Fault::phase6_sink ? 1U : 0U));
  }

  if (mpi.size() > 1) {
    const int mismatch_rank = 1;
    for (const std::size_t offset : {0U, 8U, 24U, 48U}) {
      DiagnosticAccess::reset();
      DiagnosticAccess::set_request_projection_mutation(offset, mismatch_rank);
      RecordingSink sink;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            source, mpi, collective_request, sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = true;
        HUNDUN_CHECK(error.code() ==
                     "time-control.diagnostics.request-agreement");
        HUNDUN_CHECK(error.lowest_failing_rank() == mismatch_rank);
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(sink.calls == 0);
    }
    for (const std::size_t offset : {0U, 16U, 32U, 48U, 64U}) {
      DiagnosticAccess::reset();
      DiagnosticAccess::set_provider_projection_mutation(offset,
                                                         mismatch_rank);
      RecordingSink sink;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            source, mpi, collective_request, sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = true;
        HUNDUN_CHECK(error.code() ==
                     "time-control.diagnostics.provider-agreement");
        HUNDUN_CHECK(error.lowest_failing_rank() == mismatch_rank);
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(sink.calls == 0);
    }
  }

  for (std::size_t ordinal = 1U; ordinal <= 9U; ++ordinal) {
    DiagnosticAccess::reset();
    DiagnosticAccess::set_raw_fault(ordinal);
    RecordingSink sink;
    bool typed = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, mpi, collective_request,
                                               sink);
    } catch (const hundun::runtime::MpiOperationError &) {
      typed = true;
    }
    HUNDUN_CHECK(typed);
    HUNDUN_CHECK(sink.calls == 0);
  }
  DiagnosticAccess::reset();
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run(mpi); });
}
