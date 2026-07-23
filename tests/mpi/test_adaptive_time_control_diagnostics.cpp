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
#include <cmath>
#include <cstring>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct RecordingSink final : hundun::diagnostics::DiagnosticSink {
  void submit(const hundun::diagnostics::DiagnosticRecord &value) override {
    record = value;
    ++calls;
  }
  hundun::diagnostics::DiagnosticRecord record;
  int calls{};
};

struct ThrowingSink final : hundun::diagnostics::DiagnosticSink {
  void submit(const hundun::diagnostics::DiagnosticRecord &) override {
    ++calls;
    throw std::runtime_error("injected Task22 sink rejection");
  }
  int calls{};
};

bool exact_fp64(const hundun::diagnostics::DiagnosticFp64 &actual,
                double expected) noexcept {
  std::uint64_t expected_bits{};
  std::memcpy(&expected_bits, &expected, sizeof(expected_bits));
  return actual.status ==
             hundun::diagnostics::DiagnosticValueStatus::finite &&
         actual.bits == expected_bits;
}

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
  const auto first_metadata = state.metadata();
  auto first_source = controller.diagnostic_source(state, first);
  const auto second = controller.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(second.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(state.metadata().step == 2U);
  HUNDUN_CHECK(first_metadata.step == 1U);
  return {std::move(first_source), first_metadata.step, first_metadata.time_s};
}

void run_fast(const hundun::runtime::MpiContext &mpi) {
  using DiagnosticAccess =
      hundun::diagnostics::test::TimeControlDiagnosticsTestAccess;
  using Fault = hundun::diagnostics::test::TimeControlDiagnosticFault;
  DiagnosticAccess::reset();

  // The source is owning: all controller, report, state, mesh, and solver
  // owners used to create it have already been destroyed on return.
  auto historical = make_historical_source(mpi);
  auto &source = historical.source;
  const auto descriptor =
      hundun::diagnostics::describe_time_control_diagnostics();
  HUNDUN_CHECK(descriptor.schema_version == 1U);
  HUNDUN_CHECK(descriptor.module_kind ==
               hundun::diagnostics::DiagnosticModuleKind::time_control);
  HUNDUN_CHECK(descriptor.module_id ==
               "hundun.flow.bdf2-retry-controller");
  HUNDUN_CHECK(descriptor.instance_id == "primary");
  HUNDUN_CHECK(descriptor.capabilities == 31U);
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
    HUNDUN_CHECK(sink.record.module_id ==
                 "hundun.flow.bdf2-retry-controller");
    HUNDUN_CHECK(sink.record.instance_id == "primary");
    HUNDUN_CHECK(sink.record.phase == "time-control.advance-result");
    HUNDUN_CHECK(sink.record.step == historical.step);
    HUNDUN_CHECK(sink.record.failure.code == "none");
    HUNDUN_CHECK(sink.record.status ==
                 hundun::diagnostics::DiagnosticStatus::ok);
    HUNDUN_CHECK(sink.record.failure.classification ==
                 hundun::diagnostics::DiagnosticFailureClass::none);
    HUNDUN_CHECK(sink.record.failure.lowest_failing_rank == -1);
    HUNDUN_CHECK(sink.record.identities[0].subject_id ==
                 "flow-state.accepted.cells");
    HUNDUN_CHECK(sink.record.identities[1].subject_id ==
                 "flow-state.accepted.faces");
    HUNDUN_CHECK(sink.record.identities[2].subject_id ==
                 "time-control.controller");
    HUNDUN_CHECK(sink.record.identities[3].subject_id ==
                 "time-control.report");
    HUNDUN_CHECK(
        sink.record.identities[0].layout_fingerprint.has_value());
    HUNDUN_CHECK(
        sink.record.identities[1].layout_fingerprint.has_value());
    HUNDUN_CHECK(
        sink.record.identities[0].layout_fingerprint->find(
            "cell.f64.owned.") == 0U);
    HUNDUN_CHECK(
        sink.record.identities[1].layout_fingerprint->find(
            "face.f64.owned.") == 0U);
    HUNDUN_CHECK(sink.record.identities[0].revision.has_value());
    HUNDUN_CHECK(sink.record.identities[1].revision ==
                 sink.record.identities[0].revision);
    HUNDUN_CHECK(sink.record.identities[2].revision == historical.step);
    HUNDUN_CHECK(sink.record.identities[3].revision == 1U);
    for (const auto &identity : sink.record.identities) {
      HUNDUN_CHECK(!identity.generation.has_value());
      HUNDUN_CHECK(!identity.allocation_identity.has_value());
    }
    if (level == hundun::diagnostics::DiagnosticLevel::summary) {
      HUNDUN_CHECK(sink.record.metrics.size() == 5U);
      constexpr std::array<const char *, 5> ids{
          "time-control.accepted-dt", "time-control.attempted-dt",
          "time-control.convective-number", "time-control.diffusion-number",
          "time-control.next-dt"};
      constexpr std::array<const char *, 5> units{"s", "s", "1", "1", "s"};
      constexpr std::array<double, 5> values{0.01, 0.01, 0.0, 0.0, 0.0125};
      for (std::size_t index = 0; index < ids.size(); ++index) {
        HUNDUN_CHECK(sink.record.metrics[index].id == ids[index]);
        HUNDUN_CHECK(
            sink.record.metrics[index].kind ==
            hundun::diagnostics::DiagnosticMetricKind::state_summary);
        HUNDUN_CHECK(sink.record.metrics[index].unit == units[index]);
        HUNDUN_CHECK(exact_fp64(sink.record.metrics[index].value,
                               values[index]));
      }
    }
    if (level == hundun::diagnostics::DiagnosticLevel::invariants) {
      HUNDUN_CHECK(sink.record.invariants.size() == 7U);
      constexpr std::array<const char *, 7> ids{
          "time-control.adaptive-limit-or-minimum",
          "time-control.attempt-count-bounded",
          "time-control.controller-state-valid",
          "time-control.next-dt-at-least-minimum",
          "time-control.next-dt-at-most-maximum",
          "time-control.next-dt-positive",
          "time-control.order-history-consistent"};
      constexpr std::array<const char *, 7> units{
          "1", "count", "1", "s", "s", "s", "1"};
      constexpr std::array<double, 7> observed{
          1.0, 1.0, 1.0, 0.0125, 0.0125, 0.0125, 1.0};
      constexpr std::array<double, 7> limits{
          1.0, 9.0, 1.0, 0.00125, 0.02, 0.0, 1.0};
      constexpr std::array<hundun::diagnostics::InvariantRelation, 7>
          relations{
              hundun::diagnostics::InvariantRelation::equal,
              hundun::diagnostics::InvariantRelation::less_equal,
              hundun::diagnostics::InvariantRelation::equal,
              hundun::diagnostics::InvariantRelation::greater_equal,
              hundun::diagnostics::InvariantRelation::less_equal,
              hundun::diagnostics::InvariantRelation::positive,
              hundun::diagnostics::InvariantRelation::equal};
      for (std::size_t index = 0; index < ids.size(); ++index) {
        HUNDUN_CHECK(sink.record.invariants[index].id == ids[index]);
        HUNDUN_CHECK(sink.record.invariants[index].unit == units[index]);
        HUNDUN_CHECK(exact_fp64(sink.record.invariants[index].observed,
                               observed[index]));
        if (relations[index] ==
            hundun::diagnostics::InvariantRelation::positive)
          HUNDUN_CHECK(
              sink.record.invariants[index].limit.status ==
              hundun::diagnostics::DiagnosticValueStatus::unavailable);
        else
          HUNDUN_CHECK(exact_fp64(sink.record.invariants[index].limit,
                                 limits[index]));
        HUNDUN_CHECK(sink.record.invariants[index].relation ==
                     relations[index]);
        HUNDUN_CHECK(sink.record.invariants[index].passed);
      }
    }
    if (level == hundun::diagnostics::DiagnosticLevel::counters) {
      HUNDUN_CHECK(sink.record.counters.size() == 5U);
      constexpr std::array<const char *, 5> ids{
          "time-control.accepted-step", "time-control.attempt-count",
          "time-control.controller-revision", "time-control.retry-count",
          "time-control.stability-reductions"};
      constexpr std::array<std::uint64_t, 5> values{1U, 1U, 1U, 0U, 1U};
      for (std::size_t index = 0; index < ids.size(); ++index) {
        HUNDUN_CHECK(sink.record.counters[index].id == ids[index]);
        HUNDUN_CHECK(sink.record.counters[index].unit == "count");
        HUNDUN_CHECK(sink.record.counters[index].value == values[index]);
      }
    }
    if (level ==
        hundun::diagnostics::DiagnosticLevel::bounded_state_sample) {
      const std::size_t expected = mpi.rank() == 0 ? 9U : 0U;
      HUNDUN_CHECK(sink.record.samples.size() == expected);
      HUNDUN_CHECK(sink.record.eligible_sample_count == expected);
      if (mpi.rank() == 0) {
        constexpr std::array<const char *, 9> fields{
            "time-control.accepted-step", "time-control.accepted-step",
            "time-control.history-ready", "time-control.last-accepted-dt",
            "time-control.last-order", "time-control.last-work-gate",
            "time-control.next-dt", "time-control.revision",
            "time-control.revision"};
        constexpr std::array<std::uint32_t, 9> components{
            0U, 1U, 0U, 0U, 0U, 0U, 0U, 0U, 1U};
        constexpr std::array<const char *, 9> units{
            "count", "count", "1", "s", "count",
            "1",     "s",     "count", "count"};
        constexpr std::array<double, 9> values{
            1.0, 0.0, 1.0, 0.01, 1.0, 1.0, 0.0125, 1.0, 0.0};
        hundun::diagnostics::DiagnosticFingerprintAccumulator
            expected_fingerprint;
        for (std::size_t index = 0; index < fields.size(); ++index) {
          HUNDUN_CHECK(sink.record.samples[index].field_id == fields[index]);
          HUNDUN_CHECK(sink.record.samples[index].global_id == 0U);
          HUNDUN_CHECK(sink.record.samples[index].component ==
                       components[index]);
          HUNDUN_CHECK(sink.record.samples[index].unit == units[index]);
          std::uint64_t expected_bits{};
          std::memcpy(&expected_bits, &values[index],
                      sizeof(expected_bits));
          HUNDUN_CHECK(sink.record.samples[index].value.bits ==
                       expected_bits);
          expected_fingerprint.add(
              fields[index], 0U, components[index],
              hundun::diagnostics::describe_fp64(values[index]));
        }
        HUNDUN_CHECK(sink.record.state_fingerprint.hex ==
                     expected_fingerprint.finish().hex);
      }
    }
    const auto json = hundun::diagnostics::to_canonical_json(sink.record);
    HUNDUN_CHECK(json ==
                 hundun::diagnostics::to_canonical_json(sink.record));
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

  if (mpi.rank() == 0) {
    std::string previous_fingerprint;
    for (const auto value :
         {std::uint64_t{1} << 53U, (std::uint64_t{1} << 53U) + 1U,
          std::numeric_limits<std::uint64_t>::max()}) {
      DiagnosticAccess::set_state_counters(source, value, value);
      RecordingSink limb_sink;
      const hundun::diagnostics::DiagnosticRequest limb_request{
          hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
          hundun::diagnostics::DiagnosticScope::local, frame, {}, 256U};
      hundun::diagnostics::collect_diagnostics(source, limb_request,
                                               limb_sink);
      HUNDUN_CHECK(limb_sink.record.samples.size() == 9U);
      const auto low =
          static_cast<std::uint64_t>(
              limb_sink.record.samples[0].value.bits == 0U
                  ? 0.0
                  : [&] {
                      double decoded{};
                      const auto bits = limb_sink.record.samples[0].value.bits;
                      std::memcpy(&decoded, &bits, sizeof(decoded));
                      return decoded;
                    }());
      double high_value{};
      const auto high_bits = limb_sink.record.samples[1].value.bits;
      std::memcpy(&high_value, &high_bits, sizeof(high_value));
      const auto reconstructed =
          low | (static_cast<std::uint64_t>(high_value) << 32U);
      HUNDUN_CHECK(reconstructed == value);
      HUNDUN_CHECK(previous_fingerprint.empty() ||
                   previous_fingerprint !=
                       limb_sink.record.state_fingerprint.hex);
      previous_fingerprint = limb_sink.record.state_fingerprint.hex;
      HUNDUN_CHECK(
          !hundun::diagnostics::to_canonical_json(limb_sink.record).empty());
    }
    DiagnosticAccess::set_state_counters(source, historical.step,
                                         historical.step);
  }

  struct FaultCase final {
    Fault fault;
    hundun::diagnostics::DiagnosticFailureClass classification;
    const char *code;
  };
  constexpr std::array<FaultCase, 13> faults{{
      {Fault::phase1_layout,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.local-layout"},
      {Fault::phase2_request,
       hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.diagnostics.request-preparation"},
      {Fault::projection_root_size,
       hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.diagnostics.request-preparation"},
      {Fault::projection_payload,
       hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.diagnostics.request-preparation"},
      {Fault::phase3_provider,
       hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.diagnostics.provider-agreement"},
      {Fault::phase4_payload,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.sample-preparation"},
      {Fault::fingerprint_aggregation,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.aggregation"},
      {Fault::eligible_aggregation,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.aggregation"},
      {Fault::wire_root_size,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.sample-preparation"},
      {Fault::wire_payload,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.sample-preparation"},
      {Fault::phase4_wire,
       hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.diagnostics.sample-wire"},
      {Fault::phase5_record,
       hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.diagnostics.record"},
      {Fault::phase6_sink,
       hundun::diagnostics::DiagnosticFailureClass::sink_failure,
       "diagnostics.sink.submit"},
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
      HUNDUN_CHECK(error.classification() == fault.classification);
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
        HUNDUN_CHECK(
            error.classification() ==
            hundun::diagnostics::DiagnosticFailureClass::invalid_request);
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
        HUNDUN_CHECK(
            error.classification() ==
            hundun::diagnostics::DiagnosticFailureClass::invalid_input);
        HUNDUN_CHECK(error.code() ==
                     "time-control.diagnostics.provider-agreement");
        HUNDUN_CHECK(error.lowest_failing_rank() == mismatch_rank);
      }
      HUNDUN_CHECK(rejected);
      HUNDUN_CHECK(sink.calls == 0);
    }
    for (const auto mutation :
         {hundun::diagnostics::test::TimeControlProviderMutation::state_seal,
          hundun::diagnostics::test::TimeControlProviderMutation::summary,
          hundun::diagnostics::test::TimeControlProviderMutation::config,
          hundun::diagnostics::test::TimeControlProviderMutation::model,
          hundun::diagnostics::test::TimeControlProviderMutation::frame,
          hundun::diagnostics::test::TimeControlProviderMutation::identities,
          hundun::diagnostics::test::TimeControlProviderMutation::
              global_cell_authority,
          hundun::diagnostics::test::TimeControlProviderMutation::
              global_cell_layout,
          hundun::diagnostics::test::TimeControlProviderMutation::
              global_face_authority,
          hundun::diagnostics::test::TimeControlProviderMutation::
              global_face_layout}) {
      DiagnosticAccess::reset();
      DiagnosticAccess::set_provider_mutation(mutation, mismatch_rank);
      RecordingSink sink;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            source, mpi, collective_request, sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = true;
        HUNDUN_CHECK(
            error.classification() ==
            hundun::diagnostics::DiagnosticFailureClass::invalid_input);
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
    DiagnosticAccess::set_raw_fault(ordinal, injection_rank);
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

std::uint64_t stable_text_hash(std::string_view text) noexcept {
  std::uint64_t hash = 1469598103934665603ULL;
  for (const auto byte : text) {
    hash ^= static_cast<unsigned char>(byte);
    hash *= 1099511628211ULL;
  }
  return hash;
}

void run_acceptance_only(const hundun::runtime::MpiContext &mpi) {
  using DiagnosticAccess =
      hundun::diagnostics::test::TimeControlDiagnosticsTestAccess;
  using WireMutation =
      hundun::diagnostics::test::TimeControlWireMutation;
  DiagnosticAccess::reset();
  auto historical = make_historical_source(mpi);
  auto &source = historical.source;
  const hundun::diagnostics::DiagnosticFrame frame{
      mpi.rank(), historical.step, historical.time_s,
      "time-control.advance-result"};

  const auto counters_before = mpi.fp64_reduction_counters();
  const hundun::diagnostics::DiagnosticRequest collective_request{
      hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
      hundun::diagnostics::DiagnosticScope::collective, frame, {}, 256U};
  RecordingSink collective_sink;
  hundun::diagnostics::collect_diagnostics(source, mpi, collective_request,
                                           collective_sink);
  const auto counters_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(counters_after.collective_calls ==
               counters_before.collective_calls);
  HUNDUN_CHECK(counters_after.reduced_scalars ==
               counters_before.reduced_scalars);
  HUNDUN_CHECK(counters_after.logical_payload_bytes ==
               counters_before.logical_payload_bytes);

  auto rank_neutral_record = collective_sink.record;
  rank_neutral_record.rank = 0;
  const auto canonical =
      hundun::diagnostics::to_canonical_json(rank_neutral_record);
  const auto local_hash = stable_text_hash(canonical);
  HUNDUN_CHECK(local_hash == 406274725881973949ULL);
  std::uint64_t minimum_hash{};
  std::uint64_t maximum_hash{};
  HUNDUN_CHECK(MPI_Allreduce(&local_hash, &minimum_hash, 1, MPI_UINT64_T,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&local_hash, &maximum_hash, 1, MPI_UINT64_T,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum_hash == maximum_hash);

  struct FailureExpectation final {
    hundun::diagnostics::DiagnosticFailureClass classification;
    const char *code;
  };
  const int failure_rank = mpi.size() == 1 ? 0 : 1;
  constexpr std::array<FailureExpectation, 8> preflight_failures{{
      {hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.preflight.config"},
      {hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.preflight.identity"},
      {hundun::diagnostics::DiagnosticFailureClass::layout,
       "time-control.preflight.layout"},
      {hundun::diagnostics::DiagnosticFailureClass::capability,
       "time-control.preflight.capability"},
      {hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.preflight.state"},
      {hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.preflight.transport-authority"},
      {hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.preflight.preparation"},
      {hundun::diagnostics::DiagnosticFailureClass::invalid_input,
       "time-control.preflight.report"},
  }};
  for (std::size_t index = 0; index < preflight_failures.size(); ++index) {
    DiagnosticAccess::set_outcome(
        source, hundun::flow::TimeAdvanceDisposition::non_retryable_failure,
        hundun::flow::StepFailureReason::invalid_input, failure_rank,
        static_cast<std::uint8_t>(index + 1U), 0U);
    const auto &expected = preflight_failures[index];
    for (const auto scope :
         {hundun::diagnostics::DiagnosticScope::local,
          hundun::diagnostics::DiagnosticScope::collective}) {
      RecordingSink sink;
      const hundun::diagnostics::DiagnosticRequest request{
          hundun::diagnostics::DiagnosticLevel::summary, scope, frame, {}, 0U};
      if (scope == hundun::diagnostics::DiagnosticScope::local)
        hundun::diagnostics::collect_diagnostics(source, request, sink);
      else
        hundun::diagnostics::collect_diagnostics(source, mpi, request, sink);
      HUNDUN_CHECK(sink.calls == 1);
      HUNDUN_CHECK(sink.record.status ==
                   hundun::diagnostics::DiagnosticStatus::failed);
      HUNDUN_CHECK(sink.record.failure.classification ==
                   expected.classification);
      HUNDUN_CHECK(sink.record.failure.code == expected.code);
      HUNDUN_CHECK(
          sink.record.failure.lowest_failing_rank ==
          (scope == hundun::diagnostics::DiagnosticScope::collective
               ? failure_rank
               : -1));
    }
  }

  struct ReasonExpectation final {
    hundun::flow::StepFailureReason reason;
    hundun::diagnostics::DiagnosticFailureClass classification;
    const char *code;
  };
  constexpr std::array<ReasonExpectation, 12> reason_failures{{
      {hundun::flow::StepFailureReason::momentum_linear_solve,
       hundun::diagnostics::DiagnosticFailureClass::non_convergence,
       "time-control.momentum-linear-solve"},
      {hundun::flow::StepFailureReason::pressure_linear_solve,
       hundun::diagnostics::DiagnosticFailureClass::non_convergence,
       "time-control.pressure-linear-solve"},
      {hundun::flow::StepFailureReason::non_finite_trial,
       hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
       "time-control.non-finite-trial"},
      {hundun::flow::StepFailureReason::boundary_backflow,
       hundun::diagnostics::DiagnosticFailureClass::boundary,
       "time-control.boundary-backflow"},
      {hundun::flow::StepFailureReason::transport_failure,
       hundun::diagnostics::DiagnosticFailureClass::numerical_breakdown,
       "time-control.transport-failure"},
      {hundun::flow::StepFailureReason::final_momentum_residual,
       hundun::diagnostics::DiagnosticFailureClass::non_convergence,
       "time-control.final-momentum-residual"},
      {hundun::flow::StepFailureReason::final_transport_residual,
       hundun::diagnostics::DiagnosticFailureClass::non_convergence,
       "time-control.final-transport-residual"},
      {hundun::flow::StepFailureReason::final_conservation_defect,
       hundun::diagnostics::DiagnosticFailureClass::conservation,
       "time-control.final-conservation-defect"},
      {hundun::flow::StepFailureReason::final_continuity_residual,
       hundun::diagnostics::DiagnosticFailureClass::non_convergence,
       "time-control.final-continuity-residual"},
      {hundun::flow::StepFailureReason::final_pressure_residual,
       hundun::diagnostics::DiagnosticFailureClass::non_convergence,
       "time-control.final-pressure-residual"},
      {hundun::flow::StepFailureReason::density_closure_failure,
       hundun::diagnostics::DiagnosticFailureClass::numerical_breakdown,
       "time-control.density-closure-failure"},
      {hundun::flow::StepFailureReason::collective_operation,
       hundun::diagnostics::DiagnosticFailureClass::collective_operation,
       "time-control.collective-operation"},
  }};
  for (const auto &expected : reason_failures) {
    DiagnosticAccess::set_outcome(
        source, hundun::flow::TimeAdvanceDisposition::retry_limit_reached,
        expected.reason, failure_rank, 0U, 1U);
    RecordingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::collective, frame, {}, 0U};
    hundun::diagnostics::collect_diagnostics(source, mpi, request, sink);
    HUNDUN_CHECK(sink.calls == 1);
    HUNDUN_CHECK(sink.record.status ==
                 hundun::diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(sink.record.failure.classification ==
                 expected.classification);
    HUNDUN_CHECK(sink.record.failure.code == expected.code);
    HUNDUN_CHECK(sink.record.failure.lowest_failing_rank == failure_rank);
  }
  DiagnosticAccess::set_outcome(
      source, hundun::flow::TimeAdvanceDisposition::committed,
      hundun::flow::StepFailureReason::none, -1, 0U, 1U);

  if (mpi.rank() == 0) {
    RecordingSink budget_sink;
    const hundun::diagnostics::DiagnosticRequest budget_request{
        hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
        hundun::diagnostics::DiagnosticScope::local, frame,
        {"time-control.next-dt"}, 1U};
    hundun::diagnostics::collect_diagnostics(source, budget_request,
                                             budget_sink);
    HUNDUN_CHECK(budget_sink.record.eligible_sample_count == 1U);
    HUNDUN_CHECK(budget_sink.record.samples.size() == 1U);
    HUNDUN_CHECK(budget_sink.record.samples.front().field_id ==
                 "time-control.next-dt");
    HUNDUN_CHECK(!budget_sink.record.samples_truncated);

    RecordingSink truncated_sink;
    auto truncated_request = budget_request;
    truncated_request.selected_fields = {"time-control.accepted-step"};
    hundun::diagnostics::collect_diagnostics(source, truncated_request,
                                             truncated_sink);
    HUNDUN_CHECK(truncated_sink.record.eligible_sample_count == 2U);
    HUNDUN_CHECK(truncated_sink.record.samples.size() == 1U);
    HUNDUN_CHECK(truncated_sink.record.samples_truncated);
  }

  for (const auto &selected :
       {std::vector<std::string_view>{"time-control.next-dt",
                                      "time-control.next-dt"},
        std::vector<std::string_view>{"time-control.revision",
                                      "time-control.accepted-step"},
        std::vector<std::string_view>{"time-control.unknown"}}) {
    RecordingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
        hundun::diagnostics::DiagnosticScope::local, frame, selected, 1U};
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_request);
      HUNDUN_CHECK(error.code() ==
                   "time-control.diagnostics.selected-field");
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);
  }

  {
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
        hundun::diagnostics::DiagnosticScope::collective, frame,
        {"time-control.next-dt", "time-control.next-dt"}, 1U};
    RecordingSink sink;
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, mpi, request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_request);
      HUNDUN_CHECK(error.code() ==
                   "time-control.diagnostics.selected-field");
      HUNDUN_CHECK(error.lowest_failing_rank() == 0);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);

  }

  {
    RecordingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        static_cast<hundun::diagnostics::DiagnosticLevel>(255U),
        hundun::diagnostics::DiagnosticScope::local, frame, {}, 0U};
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::capability);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.capability");
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);
  }
  {
    RecordingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::collective, frame, {}, 0U};
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::capability);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.capability");
      HUNDUN_CHECK(error.lowest_failing_rank() == -1);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);
  }

  for (const auto mutation :
       {hundun::diagnostics::test::TimeControlLocalMutation::local_box,
        hundun::diagnostics::test::TimeControlLocalMutation::
            canonical_owned_faces,
        hundun::diagnostics::test::TimeControlLocalMutation::local_faces,
        hundun::diagnostics::test::TimeControlLocalMutation::local_cell_layout,
        hundun::diagnostics::test::TimeControlLocalMutation::
            global_cell_layout,
        hundun::diagnostics::test::TimeControlLocalMutation::local_face_layout,
        hundun::diagnostics::test::TimeControlLocalMutation::
            global_face_layout}) {
    auto malformed = make_historical_source(mpi);
    DiagnosticAccess::set_local_mutation(malformed.source, mutation);
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local,
        {mpi.rank(), malformed.step, malformed.time_s,
         "time-control.advance-result"},
        {}, 0U};
    const auto before = mpi.fp64_reduction_counters();
    RecordingSink sink;
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(malformed.source, request,
                                               sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(error.classification() ==
                   hundun::diagnostics::DiagnosticFailureClass::layout);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.local-layout");
      HUNDUN_CHECK(error.lowest_failing_rank() == -1);
    }
    const auto after = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);
    HUNDUN_CHECK(after.collective_calls == before.collective_calls);
    HUNDUN_CHECK(after.reduced_scalars == before.reduced_scalars);
  }

  for (int field = 0; field < 4; ++field) {
    auto bad_frame = frame;
    if (field == 0)
      ++bad_frame.rank;
    else if (field == 1)
      ++bad_frame.step;
    else if (field == 2)
      bad_frame.time_s =
          std::nextafter(bad_frame.time_s,
                         std::numeric_limits<double>::infinity());
    else
      bad_frame.phase = "time-control.wrong-phase";
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local, bad_frame, {}, 0U};
    RecordingSink sink;
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_request);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.frame");
      HUNDUN_CHECK(error.lowest_failing_rank() == -1);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);

    const hundun::diagnostics::DiagnosticRequest collective_bad_request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::collective, bad_frame, {}, 0U};
    RecordingSink collective_bad_sink;
    rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(
          source, mpi, collective_bad_request, collective_bad_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_request);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.frame");
      HUNDUN_CHECK(error.lowest_failing_rank() == 0);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(collective_bad_sink.calls == 0);
  }

  {
    ThrowingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local, frame, {}, 0U};
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::sink_failure);
      HUNDUN_CHECK(error.code() == "diagnostics.sink.submit");
      HUNDUN_CHECK(error.lowest_failing_rank() == -1);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 1);
  }

  {
    auto moved_fixture = make_historical_source(mpi);
    auto moved = std::move(moved_fixture.source);
    RecordingSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local,
        {mpi.rank(), moved_fixture.step, moved_fixture.time_s,
         "time-control.advance-result"},
        {}, 0U};
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(moved_fixture.source, request,
                                               sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_input);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.stale-source");
      HUNDUN_CHECK(error.lowest_failing_rank() == -1);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);
    RecordingSink retained_sink;
    hundun::diagnostics::collect_diagnostics(moved, request, retained_sink);
    HUNDUN_CHECK(retained_sink.calls == 1);
  }

  const int injection_rank = mpi.size() == 1 ? 0 : 1;
  constexpr std::array<WireMutation, 12> mutations{
      WireMutation::short_payload, WireMutation::trailing_payload,
      WireMutation::count,         WireMutation::field,
      WireMutation::global_id,     WireMutation::component,
      WireMutation::unit,          WireMutation::value,
      WireMutation::duplicate,      WireMutation::missing_limb,
      WireMutation::swapped_limbs,  WireMutation::tuple_order};
  for (const auto mutation : mutations) {
    DiagnosticAccess::reset();
    DiagnosticAccess::set_wire_mutation(mutation, injection_rank);
    RecordingSink sink;
    bool rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, mpi,
                                               collective_request, sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = true;
      HUNDUN_CHECK(error.classification() ==
                   hundun::diagnostics::DiagnosticFailureClass::layout);
      HUNDUN_CHECK(error.code() == "time-control.diagnostics.sample-wire");
      HUNDUN_CHECK(error.lowest_failing_rank() == injection_rank);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.calls == 0);
  }
  DiagnosticAccess::reset();
}

void run_acceptance(const hundun::runtime::MpiContext &mpi) {
  run_fast(mpi);
  run_acceptance_only(mpi);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  if (argc != 2)
    return 2;
  const std::string mode(argv[1]);
  if (mode == "fast")
    return hundun::test::run([&] { run_fast(mpi); });
  if (mode == "acceptance")
    return hundun::test::run([&] { run_acceptance(mpi); });
  return 2;
}
