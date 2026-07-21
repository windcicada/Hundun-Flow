// SPDX-License-Identifier: Apache-2.0

#include "diagnostics/src/material_density_transport_diagnostics_test_access.hpp"
#include "flow/src/material_density_transport_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/diagnostics/material_density_transport_diagnostics.hpp"
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

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

hundun::runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {4, 1, 1};
  throw hundun::runtime::Error("unsupported diagnostic rank count");
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
  for (std::size_t i = 0; i < names.size(); ++i) {
    config.boundaries[i].patch = names[i];
    config.boundaries[i].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components = 1U,
                                      bool conservative = true) {
  return {name,
          unit,
          "material-density-diagnostics",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name,
          "1",
          "material-density-diagnostics",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

class CaptureSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &record) override {
    ++calls;
    value = record;
    if (throw_classified)
      throw hundun::diagnostics::DiagnosticCollectionError(
          hundun::diagnostics::DiagnosticFailureClass::conservation,
          "sink.chosen-class", 42, "sink classified fixture");
    if (throw_after_copy)
      throw std::runtime_error("sink fixture");
  }
  std::size_t calls{};
  bool throw_after_copy{};
  bool throw_classified{};
  hundun::diagnostics::DiagnosticRecord value;
};

const hundun::diagnostics::DiagnosticMetric &
find_metric(const hundun::diagnostics::DiagnosticRecord &record,
            std::string_view id) {
  const auto found =
      std::find_if(record.metrics.begin(), record.metrics.end(),
                   [&](const auto &value) { return value.id == id; });
  HUNDUN_CHECK(found != record.metrics.end());
  return *found;
}

const hundun::diagnostics::DiagnosticInvariant &
find_invariant(const hundun::diagnostics::DiagnosticRecord &record,
               std::string_view id) {
  const auto found =
      std::find_if(record.invariants.begin(), record.invariants.end(),
                   [&](const auto &value) { return value.id == id; });
  HUNDUN_CHECK(found != record.invariants.end());
  return *found;
}

const hundun::diagnostics::DiagnosticCounter &
find_counter(const hundun::diagnostics::DiagnosticRecord &record,
             std::string_view id) {
  const auto found =
      std::find_if(record.counters.begin(), record.counters.end(),
                   [&](const auto &value) { return value.id == id; });
  HUNDUN_CHECK(found != record.counters.end());
  return *found;
}

void run(const hundun::runtime::MpiContext &mpi) {
  using namespace hundun;
  constexpr runtime::Int3 extent{8, 4, 4};
  auto decomposition = runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  mesh::MeshTopology topology(decomposition);
  mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      boundary::BoundaryRegistry::create(periodic_case(), topology);

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", "kg/m3"));
  fields.velocity = registry.declare_field(cell("velocity", "1", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(cell("pi", "Pa", 1U, false));
  fields.face_velocity = registry.declare_field(face("u_face", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  const auto rho_h = registry.declare_field(cell("rho_h", "J/m3"));
  const auto rho_phi = registry.declare_field(cell("rho_phi", "kg/m3"));
  const auto rho_phi_beta =
      registry.declare_field(cell("rho_phi_beta", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_phi, rho_phi_beta};
  registry.freeze();
  const auto box = decomposition.owned_box();
  const runtime::Int3 local{box.end.x - box.begin.x, box.end.y - box.begin.y,
                            box.end.z - box.begin.z};
  auto state = flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
  flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), 1.0);
  initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 99.0);
  initial.transported_cell_fields = {
      std::vector<double>(topology.owned_cell_count(), 3.0),
      std::vector<double>(topology.owned_cell_count(), 0.4),
      std::vector<double>(topology.owned_cell_count(), 0.7)};
  state.seed_accepted_layers(initial, initial);

  runtime::FieldStorage flux_storage(
      registry, runtime::FieldLayoutSet{local, topology.local_face_count()});
  constexpr runtime::PhaseId phase = 1910U;
  constexpr runtime::ActorId actor = 1910U;
  runtime::FieldAccessPlan access(registry);
  access.declare_access(phase, actor, fields.face_mass_flux,
                        runtime::AccessMode::read_write);
  access.freeze();
  auto flux_write = flux_storage.acquire_face_write<double>(
      access, phase, actor, fields.face_mass_flux);
  for (std::size_t face_id = 0; face_id < topology.local_face_count();
       ++face_id)
    flux_write(face_id, 0) = 0.0;
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(decomposition, local, 2));
  flow::MaterialDensityTransportSpec spec;
  spec.enthalpy_density = rho_h;
  spec.scalar_densities = {rho_phi, rho_phi_beta};
  spec.scalar_diffusivities_kg_per_m_s = {0.0, 0.0};
  auto transport = flow::MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, halo,
      fields, spec);
  const auto stencil = flow::make_momentum_time_stencil(
      flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
  state.begin_attempt();
  auto flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::final_corrected);
  const auto report = transport.finalize_trial(state, flux, stencil);
  HUNDUN_CHECK(report.disposition() ==
               flow::MaterialTransportDisposition::finalized);
  auto source = [&]() {
    const auto transient_copy = report;
    return transport.diagnostic_source(state, transient_copy);
  }();
  const auto descriptor = diagnostics::describe_diagnostics(source);
  diagnostics::validate(descriptor);
  const auto ids = diagnostics::diagnostic_fingerprint_field_ids(source);
  HUNDUN_CHECK(ids.size() == 4U);
  HUNDUN_CHECK(ids[0] == "rho");
  HUNDUN_CHECK(ids[1] == "rho_h");
  HUNDUN_CHECK(ids[2] == "rho_phi.s00000000000000000000");
  HUNDUN_CHECK(ids[3] == "rho_phi.s00000000000000000001");
  const auto history_state_before = state.snapshot(flow::FlowLayer::history);
  const auto committed_state_before =
      state.snapshot(flow::FlowLayer::committed);
  const auto state_before = state.snapshot(flow::FlowLayer::trial);
  const auto metadata_before = state.metadata();
  const auto counters_before = mpi.fp64_reduction_counters();
  std::string fingerprint;
  for (const auto level :
       {diagnostics::DiagnosticLevel::summary,
        diagnostics::DiagnosticLevel::invariants,
        diagnostics::DiagnosticLevel::counters,
        diagnostics::DiagnosticLevel::bounded_state_sample}) {
    diagnostics::DiagnosticRequest request;
    request.level = level;
    request.scope = diagnostics::DiagnosticScope::local;
    request.frame = {mpi.rank(), 0U, -0.0, "material.finalized-trial"};
    if (level == diagnostics::DiagnosticLevel::bounded_state_sample)
      request.sample_budget = 2U;
    CaptureSink sink;
    diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::reset();
    diagnostics::collect_diagnostics(source, request, sink);
    const auto work = diagnostics::test::
        MaterialDensityTransportDiagnosticsTestAccess::work_counts();
    HUNDUN_CHECK(sink.calls == 1U);
    HUNDUN_CHECK(sink.value.status == diagnostics::DiagnosticStatus::ok);
    HUNDUN_CHECK(sink.value.metrics.empty() ==
                 (level != diagnostics::DiagnosticLevel::summary));
    HUNDUN_CHECK(sink.value.invariants.empty() ==
                 (level != diagnostics::DiagnosticLevel::invariants));
    HUNDUN_CHECK(sink.value.counters.empty() ==
                 (level != diagnostics::DiagnosticLevel::counters));
    if (fingerprint.empty())
      fingerprint = sink.value.state_fingerprint.hex;
    HUNDUN_CHECK(sink.value.state_fingerprint.hex == fingerprint);
    HUNDUN_CHECK(work.fingerprint_items ==
                 source.fingerprint_field_count() * source.owned_cell_count());
    HUNDUN_CHECK(work.raw_collectives == 0U);
    if (level == diagnostics::DiagnosticLevel::counters) {
      HUNDUN_CHECK(work.summary_accumulations == 0U);
      HUNDUN_CHECK(work.sample_candidates == 0U);
      HUNDUN_CHECK(work.volume_reads == 0U);
    } else if (level == diagnostics::DiagnosticLevel::bounded_state_sample) {
      HUNDUN_CHECK(work.summary_accumulations == 0U);
      HUNDUN_CHECK(work.sample_candidates > 0U);
      HUNDUN_CHECK(work.volume_reads == 0U);
    } else if (level == diagnostics::DiagnosticLevel::summary) {
      HUNDUN_CHECK(work.summary_accumulations > 0U);
      HUNDUN_CHECK(work.sample_candidates == 0U);
      HUNDUN_CHECK(work.volume_reads > 0U);
    } else {
      HUNDUN_CHECK(work.summary_accumulations > 0U);
      HUNDUN_CHECK(work.sample_candidates == 0U);
      HUNDUN_CHECK(work.volume_reads == 0U);
    }
    diagnostics::validate(sink.value, descriptor, request);
  }
  HUNDUN_CHECK(
      runtime::Fp64ReductionCounters{counters_before}.collective_calls ==
      mpi.fp64_reduction_counters().collective_calls);
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      state_before, state.snapshot(flow::FlowLayer::trial)));
  HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(metadata_before,
                                                          state.metadata()));

  diagnostics::DiagnosticRequest collective;
  collective.level = diagnostics::DiagnosticLevel::bounded_state_sample;
  collective.scope = diagnostics::DiagnosticScope::collective;
  collective.frame = {mpi.rank(), 0U, -0.0, "material.finalized-trial"};
  collective.selected_fields = {"rho_h"};
  collective.sample_budget = 2U;
  CaptureSink collective_sink;
  diagnostics::collect_diagnostics(source, mpi, collective, collective_sink);
  HUNDUN_CHECK(collective_sink.calls == 1U);
  HUNDUN_CHECK(collective_sink.value.samples.size() == 2U);
  for (const auto &sample : collective_sink.value.samples)
    HUNDUN_CHECK(sample.field_id == "rho_h");
  HUNDUN_CHECK(mpi.fp64_reduction_counters().collective_calls ==
               counters_before.collective_calls);

  std::string collective_fingerprint;
  std::string repeated_collective_json;
  for (const std::size_t budget : {1U, 256U}) {
    for (const bool narrowed : {false, true}) {
      auto sampled = collective;
      sampled.sample_budget = budget;
      sampled.selected_fields = narrowed
                                    ? std::vector<std::string_view>{"rho_h"}
                                    : std::vector<std::string_view>{};
      CaptureSink sampled_sink;
      diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::reset();
      diagnostics::collect_diagnostics(source, mpi, sampled, sampled_sink);
      const auto sampled_work = diagnostics::test::
          MaterialDensityTransportDiagnosticsTestAccess::work_counts();
      const std::uint64_t eligible =
          static_cast<std::uint64_t>(extent.x * extent.y * extent.z) *
          (narrowed ? 1U : ids.size());
      HUNDUN_CHECK(sampled_sink.value.eligible_sample_count == eligible);
      HUNDUN_CHECK(sampled_sink.value.samples.size() ==
                   std::min<std::uint64_t>(budget, eligible));
      HUNDUN_CHECK(sampled_sink.value.samples_truncated == (eligible > budget));
      HUNDUN_CHECK(sampled_work.summary_accumulations == 0U);
      HUNDUN_CHECK(sampled_work.volume_reads == 0U);
      HUNDUN_CHECK(sampled_work.sample_candidates > 0U);
      if (collective_fingerprint.empty())
        collective_fingerprint = sampled_sink.value.state_fingerprint.hex;
      HUNDUN_CHECK(sampled_sink.value.state_fingerprint.hex ==
                   collective_fingerprint);
      if (budget == 1U && narrowed) {
        repeated_collective_json =
            diagnostics::to_canonical_json(sampled_sink.value);
        CaptureSink repeated_sink;
        diagnostics::collect_diagnostics(source, mpi, sampled, repeated_sink);
        HUNDUN_CHECK(diagnostics::to_canonical_json(repeated_sink.value) ==
                     repeated_collective_json);
      }
    }
  }
  diagnostics::DiagnosticFingerprintAccumulator expected_fingerprint;
  const std::array<double, 4> expected_values{1.0, 3.0, 0.4, 0.7};
  for (std::size_t field = 0; field < ids.size(); ++field) {
    for (std::uint64_t global = 0U;
         global < static_cast<std::uint64_t>(extent.x * extent.y * extent.z);
         ++global)
      expected_fingerprint.add(
          ids[field], global, 0U,
          diagnostics::describe_fp64(expected_values[field]));
  }
  HUNDUN_CHECK(collective_fingerprint == expected_fingerprint.finish().hex);

  runtime::FieldAccessPlan mutation_access(registry);
  constexpr runtime::PhaseId mutation_phase = 1918U;
  constexpr runtime::ActorId mutation_actor = 1918U;
  mutation_access.declare_access(mutation_phase, mutation_actor, rho_phi_beta,
                                 runtime::AccessMode::read_write);
  mutation_access.freeze();
  auto beta_write = state.trial_layer().acquire_write<double>(
      mutation_access, mutation_phase, mutation_actor, rho_phi_beta);
  const double authoritative_before = beta_write(0, 0, 0, 0);
  if (mpi.rank() == 0)
    beta_write(0, 0, 0, 0) = std::nextafter(authoritative_before, 1.0);
  auto narrow_prefix = collective;
  narrow_prefix.sample_budget = 1U;
  narrow_prefix.selected_fields = {"rho_h"};
  CaptureSink changed_fingerprint_sink;
  diagnostics::collect_diagnostics(source, mpi, narrow_prefix,
                                   changed_fingerprint_sink);
  HUNDUN_CHECK(changed_fingerprint_sink.value.state_fingerprint.hex !=
               collective_fingerprint);
  if (mpi.rank() == 0)
    beta_write(0, 0, 0, 0) = authoritative_before;
  const double ghost_before = beta_write(-1, 0, 0, 0);
  if (mpi.rank() == 0)
    beta_write(-1, 0, 0, 0) = std::nextafter(ghost_before, 1.0);
  auto external_flux_write = flux_storage.acquire_face_write<double>(
      access, phase, actor, fields.face_mass_flux);
  const double external_flux_before = external_flux_write(0U, 0);
  external_flux_write(0U, 0) = 77.0;
  CaptureSink excluded_state_sink;
  diagnostics::collect_diagnostics(source, mpi, narrow_prefix,
                                   excluded_state_sink);
  HUNDUN_CHECK(excluded_state_sink.value.state_fingerprint.hex ==
               collective_fingerprint);
  if (mpi.rank() == 0)
    beta_write(-1, 0, 0, 0) = ghost_before;
  external_flux_write(0U, 0) = external_flux_before;

  for (const auto level : {diagnostics::DiagnosticLevel::summary,
                           diagnostics::DiagnosticLevel::invariants,
                           diagnostics::DiagnosticLevel::counters}) {
    auto level_request = collective;
    level_request.level = level;
    level_request.selected_fields.clear();
    level_request.sample_budget = 0U;
    CaptureSink level_sink;
    diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::reset();
    diagnostics::collect_diagnostics(source, mpi, level_request, level_sink);
    const auto work = diagnostics::test::
        MaterialDensityTransportDiagnosticsTestAccess::work_counts();
    HUNDUN_CHECK(work.sample_candidates == 0U);
    if (level == diagnostics::DiagnosticLevel::counters) {
      HUNDUN_CHECK(work.summary_accumulations == 0U);
      HUNDUN_CHECK(work.volume_reads == 0U);
    } else if (level == diagnostics::DiagnosticLevel::summary) {
      HUNDUN_CHECK(work.summary_accumulations > 0U);
      HUNDUN_CHECK(work.volume_reads > 0U);
    } else {
      HUNDUN_CHECK(work.summary_accumulations > 0U);
      HUNDUN_CHECK(work.volume_reads == 0U);
    }
  }

  if (mpi.size() > 1) {
    for (int dimension = 0; dimension < 7; ++dimension) {
      auto disagreement = collective;
      if (mpi.rank() == 1) {
        switch (dimension) {
        case 0:
          disagreement.frame.rank = 0;
          break;
        case 1:
          disagreement.level = diagnostics::DiagnosticLevel::summary;
          disagreement.selected_fields.clear();
          disagreement.sample_budget = 0U;
          break;
        case 2:
          ++disagreement.frame.step;
          break;
        case 3:
          disagreement.frame.time_s = 0.0;
          break;
        case 4:
          disagreement.frame.phase = "material.other";
          break;
        case 5:
          disagreement.selected_fields = {"rho"};
          break;
        default:
          disagreement.sample_budget = 1U;
          break;
        }
      }
      CaptureSink disagreement_sink;
      bool mismatch_rejected = false;
      try {
        diagnostics::collect_diagnostics(source, mpi, disagreement,
                                         disagreement_sink);
      } catch (const diagnostics::DiagnosticCollectionError &error) {
        mismatch_rejected = true;
        HUNDUN_CHECK(error.classification() ==
                     diagnostics::DiagnosticFailureClass::invalid_request);
        HUNDUN_CHECK(error.lowest_failing_rank() == 1);
      }
      HUNDUN_CHECK(mismatch_rejected);
      HUNDUN_CHECK(disagreement_sink.calls == 0U);
    }
  }

  diagnostics::DiagnosticRequest sink_request;
  sink_request.level = diagnostics::DiagnosticLevel::summary;
  sink_request.scope = diagnostics::DiagnosticScope::collective;
  sink_request.frame = {mpi.rank(), 0U, -0.0, "material.finalized-trial"};
  CaptureSink failing_sink;
  const int sink_failure_rank = mpi.size() > 1 ? 1 : 0;
  failing_sink.throw_after_copy = mpi.rank() == sink_failure_rank;
  bool sink_rejected = false;
  try {
    diagnostics::collect_diagnostics(source, mpi, sink_request, failing_sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    sink_rejected = true;
    HUNDUN_CHECK(error.classification() ==
                 diagnostics::DiagnosticFailureClass::sink_failure);
    HUNDUN_CHECK(error.lowest_failing_rank() == sink_failure_rank);
  }
  HUNDUN_CHECK(sink_rejected);
  HUNDUN_CHECK(failing_sink.calls == 1U);

  diagnostics::DiagnosticRequest local_sink_request = sink_request;
  local_sink_request.scope = diagnostics::DiagnosticScope::local;
  for (const bool classified : {false, true}) {
    CaptureSink local_failing_sink;
    local_failing_sink.throw_after_copy = !classified;
    local_failing_sink.throw_classified = classified;
    bool local_sink_rejected = false;
    try {
      diagnostics::collect_diagnostics(source, local_sink_request,
                                       local_failing_sink);
    } catch (const diagnostics::DiagnosticCollectionError &error) {
      local_sink_rejected = true;
      HUNDUN_CHECK(error.classification() ==
                   diagnostics::DiagnosticFailureClass::sink_failure);
      HUNDUN_CHECK(error.code() == "diagnostics.sink.submit");
      HUNDUN_CHECK(error.lowest_failing_rank() == -1);
    }
    HUNDUN_CHECK(local_sink_rejected);
    HUNDUN_CHECK(local_failing_sink.calls == 1U);
  }

  CaptureSink classified_collective_sink;
  classified_collective_sink.throw_classified = mpi.rank() == sink_failure_rank;
  sink_rejected = false;
  try {
    diagnostics::collect_diagnostics(source, mpi, sink_request,
                                     classified_collective_sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    sink_rejected = true;
    HUNDUN_CHECK(error.classification() ==
                 diagnostics::DiagnosticFailureClass::sink_failure);
    HUNDUN_CHECK(error.code() == "diagnostics.sink.submit");
    HUNDUN_CHECK(error.lowest_failing_rank() == sink_failure_rank);
  }
  HUNDUN_CHECK(sink_rejected);

  auto require_injected_failure = [&](auto arm,
                                      diagnostics::DiagnosticFailureClass klass,
                                      std::string_view code,
                                      int expected_rank) {
    diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::reset();
    arm();
    CaptureSink injected_sink;
    bool injected_rejected = false;
    try {
      diagnostics::collect_diagnostics(source, mpi, collective, injected_sink);
    } catch (const diagnostics::DiagnosticCollectionError &error) {
      injected_rejected = true;
      HUNDUN_CHECK(error.classification() == klass);
      HUNDUN_CHECK(error.code() == code);
      HUNDUN_CHECK(error.lowest_failing_rank() == expected_rank);
    }
    HUNDUN_CHECK(injected_rejected);
    HUNDUN_CHECK(injected_sink.calls == 0U);
  };
  const auto first_global_id = source.global_cell_id(0U);
  require_injected_failure(
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            override_global_id(1U, first_global_id, 0);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.duplicate-owned-sample", 0);
  const auto outside_global_id = source.global_cell_id(9U);
  require_injected_failure(
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            override_global_id(10U, outside_global_id, 0);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.duplicate-owned-sample", 0);
  const int injection_rank = mpi.size() > 1 ? 1 : 0;
  require_injected_failure(
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_provider_failure(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.injected-provider", injection_rank);
  require_injected_failure(
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_record_failure(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.injected-record", injection_rank);
  require_injected_failure(
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_request_size_overflow(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::invalid_request,
      "material.diagnostics.request-size", injection_rank);
  require_injected_failure(
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_sample_wire_overflow(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.sample-wire-size", 0);

  for (const auto scope : {diagnostics::DiagnosticScope::local,
                           diagnostics::DiagnosticScope::collective}) {
    auto wrong_phase = sink_request;
    wrong_phase.scope = scope;
    wrong_phase.frame.phase = "material.other";
    CaptureSink phase_sink;
    bool phase_rejected = false;
    try {
      if (scope == diagnostics::DiagnosticScope::local)
        diagnostics::collect_diagnostics(source, wrong_phase, phase_sink);
      else
        diagnostics::collect_diagnostics(source, mpi, wrong_phase, phase_sink);
    } catch (const diagnostics::DiagnosticCollectionError &error) {
      phase_rejected = true;
      HUNDUN_CHECK(error.classification() ==
                   diagnostics::DiagnosticFailureClass::invalid_request);
      HUNDUN_CHECK(error.code() == "material.diagnostics.phase");
    }
    HUNDUN_CHECK(phase_rejected);
    HUNDUN_CHECK(phase_sink.calls == 0U);
  }
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      history_state_before, state.snapshot(flow::FlowLayer::history)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      committed_state_before, state.snapshot(flow::FlowLayer::committed)));
  HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
      state_before, state.snapshot(flow::FlowLayer::trial)));
  HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(metadata_before,
                                                          state.metadata()));
  HUNDUN_CHECK(mpi.fp64_reduction_counters().collective_calls ==
               counters_before.collective_calls);

  auto bad = collective;
  bad.selected_fields = {"unknown"};
  CaptureSink no_submit;
  bool rejected = false;
  try {
    diagnostics::collect_diagnostics(source, mpi, bad, no_submit);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    rejected = true;
    HUNDUN_CHECK(error.classification() ==
                 diagnostics::DiagnosticFailureClass::invalid_request);
    HUNDUN_CHECK(error.code() == "material.diagnostics.unknown-field");
    HUNDUN_CHECK(error.lowest_failing_rank() == 0);
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(no_submit.calls == 0U);

  for (const auto corruption :
       {flow::test::MaterialReportCorruption::scalar,
        flow::test::MaterialReportCorruption::vector_size,
        flow::test::MaterialReportCorruption::vector_element,
        flow::test::MaterialReportCorruption::availability,
        flow::test::MaterialReportCorruption::seal}) {
    auto mismatched_report = report;
    flow::test::MaterialDensityTransportTestAccess::corrupt_report(
        mismatched_report, corruption);
    rejected = false;
    try {
      static_cast<void>(transport.diagnostic_source(state, mismatched_report));
    } catch (const runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }

  auto later_flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::final_corrected);
  const auto later_report =
      transport.finalize_trial(state, later_flux, stencil);
  HUNDUN_CHECK(later_report.disposition() ==
               flow::MaterialTransportDisposition::finalized);
  auto later_source = transport.diagnostic_source(state, later_report);
  CaptureSink superseded_sink;
  diagnostics::DiagnosticRequest stale;
  stale.level = diagnostics::DiagnosticLevel::summary;
  stale.scope = diagnostics::DiagnosticScope::local;
  stale.frame = {mpi.rank(), 0U, 0.0, "material.finalized-trial"};
  rejected = false;
  try {
    diagnostics::collect_diagnostics(source, stale, superseded_sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    rejected = true;
    HUNDUN_CHECK(error.classification() ==
                 diagnostics::DiagnosticFailureClass::invalid_input);
    HUNDUN_CHECK(error.code() == "material.diagnostics.provider");
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(superseded_sink.calls == 0U);

  auto other_state = flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
  other_state.seed_accepted_layers(initial, initial);
  other_state.begin_attempt();
  rejected = false;
  try {
    static_cast<void>(transport.diagnostic_source(other_state, later_report));
  } catch (const runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  other_state.rollback_attempt();

  state.rollback_attempt();
  CaptureSink stale_sink;
  rejected = false;
  try {
    diagnostics::collect_diagnostics(later_source, stale, stale_sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    rejected = true;
    HUNDUN_CHECK(error.classification() ==
                 diagnostics::DiagnosticFailureClass::invalid_input);
    HUNDUN_CHECK(error.code() == "material.diagnostics.provider");
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(stale_sink.calls == 0U);

  state.begin_attempt();
  auto commit_flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::final_corrected);
  const auto commit_report =
      transport.finalize_trial(state, commit_flux, stencil);
  HUNDUN_CHECK(commit_report.disposition() ==
               flow::MaterialTransportDisposition::finalized);
  auto commit_source = transport.diagnostic_source(state, commit_report);
  state.commit_attempt(
      {1U, 0.01, 0.01, 0.0, flow::MomentumTimeOrder::backward_euler});
  CaptureSink committed_sink;
  rejected = false;
  try {
    diagnostics::collect_diagnostics(commit_source, stale, committed_sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    rejected = true;
    HUNDUN_CHECK(error.classification() ==
                 diagnostics::DiagnosticFailureClass::invalid_input);
    HUNDUN_CHECK(error.code() == "material.diagnostics.provider");
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(committed_sink.calls == 0U);

  const auto check_failed_record_levels =
      [&](const flow::MaterialDensityTransportReport &failed_report,
          diagnostics::DiagnosticFailureClass expected_class,
          bool assessments_available) {
        auto failed_source = transport.diagnostic_source(state, failed_report);
        for (const auto scope : {diagnostics::DiagnosticScope::local,
                                 diagnostics::DiagnosticScope::collective}) {
          for (const auto level :
               {diagnostics::DiagnosticLevel::summary,
                diagnostics::DiagnosticLevel::invariants,
                diagnostics::DiagnosticLevel::counters,
                diagnostics::DiagnosticLevel::bounded_state_sample}) {
            diagnostics::DiagnosticRequest request;
            request.level = level;
            request.scope = scope;
            request.frame = {mpi.rank(), 0U, 0.0, "material.finalized-trial"};
            if (level == diagnostics::DiagnosticLevel::bounded_state_sample)
              request.sample_budget = 1U;
            CaptureSink failed_sink;
            if (scope == diagnostics::DiagnosticScope::local)
              diagnostics::collect_diagnostics(failed_source, request,
                                               failed_sink);
            else
              diagnostics::collect_diagnostics(failed_source, mpi, request,
                                               failed_sink);
            HUNDUN_CHECK(failed_sink.value.status ==
                         diagnostics::DiagnosticStatus::failed);
            HUNDUN_CHECK(failed_sink.value.failure.classification ==
                         expected_class);
            if (level == diagnostics::DiagnosticLevel::summary) {
              const auto &density_residual =
                  find_metric(failed_sink.value, "density.residual");
              HUNDUN_CHECK((density_residual.value.status !=
                            diagnostics::DiagnosticValueStatus::unavailable) ==
                           assessments_available);
              const auto &last_transport =
                  find_metric(failed_sink.value,
                              "transport.s00000000000000000002.residual");
              HUNDUN_CHECK((last_transport.value.status !=
                            diagnostics::DiagnosticValueStatus::unavailable) ==
                           assessments_available);
            } else if (level == diagnostics::DiagnosticLevel::invariants) {
              const auto &availability = find_invariant(
                  failed_sink.value, "density.residual.available");
              HUNDUN_CHECK(availability.passed == assessments_available);
              const auto &last_availability = find_invariant(
                  failed_sink.value,
                  "transport.s00000000000000000002.residual.available");
              HUNDUN_CHECK(last_availability.passed == assessments_available);
            } else if (level == diagnostics::DiagnosticLevel::counters) {
              HUNDUN_CHECK(
                  find_counter(failed_sink.value, "transport.fields").value ==
                  3U);
            }
          }
        }
      };

  state.seed_accepted_layers(initial, initial);
  state.begin_attempt();
  auto predictor_failure_flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::predictor);
  const auto predictor_failure_report =
      transport.finalize_trial(state, predictor_failure_flux, stencil);
  check_failed_record_levels(predictor_failure_report,
                             diagnostics::DiagnosticFailureClass::invalid_input,
                             false);
  state.rollback_attempt();

  auto invalid_initial = initial;
  if (mpi.rank() == 0)
    invalid_initial.density.front() = -1.0;
  state.seed_accepted_layers(invalid_initial, invalid_initial);
  state.begin_attempt();
  auto state_failure_flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::final_corrected);
  const auto state_failure_report =
      transport.finalize_trial(state, state_failure_flux, stencil);
  check_failed_record_levels(
      state_failure_report,
      diagnostics::DiagnosticFailureClass::non_positive_state, false);
  state.rollback_attempt();

  state.seed_accepted_layers(initial, initial);
  state.begin_attempt();
  flow::test::MaterialDensityTransportTestAccess::set_density_residual(2.0e-10);
  auto final_gate_flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::final_corrected);
  const auto final_gate_report =
      transport.finalize_trial(state, final_gate_flux, stencil);
  check_failed_record_levels(
      final_gate_report, diagnostics::DiagnosticFailureClass::non_convergence,
      true);
  state.rollback_attempt();
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    run(mpi);
  });
}
