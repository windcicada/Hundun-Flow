// SPDX-License-Identifier: Apache-2.0

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

#include <array>
#include <stdexcept>
#include <string>
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
    if (throw_after_copy)
      throw std::runtime_error("sink fixture");
  }
  std::size_t calls{};
  bool throw_after_copy{};
  hundun::diagnostics::DiagnosticRecord value;
};

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
  fields.transported_cell_fields = {rho_h, rho_phi};
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
      std::vector<double>(topology.owned_cell_count(), 0.4)};
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
  spec.scalar_densities = {rho_phi};
  spec.scalar_diffusivities_kg_per_m_s = {0.0};
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
  HUNDUN_CHECK(report.disposition ==
               flow::MaterialTransportDisposition::finalized);
  auto source = transport.diagnostic_source(state, report);
  const auto descriptor = diagnostics::describe_diagnostics(source);
  diagnostics::validate(descriptor);
  const auto ids = diagnostics::diagnostic_fingerprint_field_ids(source);
  HUNDUN_CHECK(ids.size() == 3U);
  HUNDUN_CHECK(ids[0] == "rho");
  HUNDUN_CHECK(ids[1] == "rho_h");
  HUNDUN_CHECK(ids[2] == "rho_phi.s00000000000000000000");
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
    diagnostics::collect_diagnostics(source, request, sink);
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

  if (mpi.size() > 1) {
    auto disagreement = collective;
    if (mpi.rank() == 1)
      disagreement.frame.time_s = 0.0;
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

  auto mismatched_report = report;
  ++mismatched_report.finalization_identity;
  rejected = false;
  try {
    static_cast<void>(transport.diagnostic_source(state, mismatched_report));
  } catch (const runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto later_flux = flow::MaterialFaceMassFlux::acquire(
      registry, flux_storage, access, phase, actor, fields.face_mass_flux,
      topology, flow::MaterialFluxProvenance::final_corrected);
  const auto later_report =
      transport.finalize_trial(state, later_flux, stencil);
  HUNDUN_CHECK(later_report.disposition ==
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
  HUNDUN_CHECK(commit_report.disposition ==
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
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    run(mpi);
  });
}
