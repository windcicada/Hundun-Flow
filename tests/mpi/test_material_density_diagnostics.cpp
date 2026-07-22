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
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
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

struct ExpectedSampleSet final {
  std::uint64_t eligible_count{};
  bool truncated{};
  std::vector<hundun::diagnostics::DiagnosticSample> samples;
};

ExpectedSampleSet
expected_samples(const std::vector<std::string_view> &field_ids,
                 const std::vector<std::uint64_t> &global_ids, bool narrowed,
                 std::size_t budget) {
  using namespace hundun::diagnostics;
  constexpr std::array<double, 4> values{1.0, 3.0, 0.4, 0.7};
  ExpectedSampleSet result;
  for (std::size_t field = 0; field < field_ids.size(); ++field) {
    if (narrowed && field_ids[field] != "rho_h")
      continue;
    const std::string_view unit =
        field == 0U ? "kg/m3" : (field == 1U ? "J/m3" : "kg/m3");
    for (const auto global : global_ids) {
      result.samples.push_back({std::string(field_ids[field]), global, 0U,
                                std::string(unit),
                                describe_fp64(values[field])});
    }
  }
  std::sort(result.samples.begin(), result.samples.end(),
            [](const auto &left, const auto &right) {
              return std::tie(left.field_id, left.global_id, left.component) <
                     std::tie(right.field_id, right.global_id, right.component);
            });
  result.eligible_count = static_cast<std::uint64_t>(result.samples.size());
  if (result.samples.size() > budget)
    result.samples.resize(budget);
  result.truncated = result.eligible_count > result.samples.size();
  return result;
}

bool exact_samples_equal(
    const std::vector<hundun::diagnostics::DiagnosticSample> &left,
    const std::vector<hundun::diagnostics::DiagnosticSample> &right) {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (left[index].field_id != right[index].field_id ||
        left[index].global_id != right[index].global_id ||
        left[index].component != right[index].component ||
        left[index].unit != right[index].unit ||
        left[index].value.status != right[index].value.status ||
        left[index].value.bits != right[index].value.bits)
      return false;
  }
  return true;
}

void require_expected_samples(
    const hundun::diagnostics::DiagnosticRecord &record,
    const ExpectedSampleSet &expected) {
  HUNDUN_CHECK(record.eligible_sample_count == expected.eligible_count);
  HUNDUN_CHECK(record.samples_truncated == expected.truncated);
  HUNDUN_CHECK(exact_samples_equal(record.samples, expected.samples));
}

std::string owned_layout(const hundun::runtime::Box3 &box) {
  return "cell.f64.c1.g2plus.owned." + std::to_string(box.begin.x) + "." +
         std::to_string(box.begin.y) + "." + std::to_string(box.begin.z) + "." +
         std::to_string(box.end.x) + "." + std::to_string(box.end.y) + "." +
         std::to_string(box.end.z);
}

std::string global_layout(hundun::runtime::Int3 extent) {
  return "cell.f64.c1.g2plus.global." + std::to_string(extent.x) + "." +
         std::to_string(extent.y) + "." + std::to_string(extent.z);
}

void require_identities(
    const hundun::diagnostics::DiagnosticRecord &record,
    const std::vector<std::string_view> &field_ids,
    std::string_view expected_layout,
    const hundun::flow::MaterialDensityTransportReport &report) {
  const std::size_t expected_count = field_ids.size() + 3U;
  HUNDUN_CHECK(record.identities.size() == expected_count);
  for (std::size_t field = 0; field < field_ids.size(); ++field) {
    const auto &identity = record.identities[field];
    HUNDUN_CHECK(identity.subject_id ==
                 std::string("field.") + std::string(field_ids[field]));
    HUNDUN_CHECK(identity.layout_fingerprint.has_value());
    HUNDUN_CHECK(*identity.layout_fingerprint == expected_layout);
    HUNDUN_CHECK(!identity.revision.has_value());
    HUNDUN_CHECK(!identity.generation.has_value());
    HUNDUN_CHECK(!identity.allocation_identity.has_value());
  }
  const auto &attempt = record.identities[field_ids.size()];
  HUNDUN_CHECK(attempt.subject_id == "flow_state.attempt");
  HUNDUN_CHECK(!attempt.layout_fingerprint.has_value());
  HUNDUN_CHECK(attempt.revision == report.attempt_identity());
  HUNDUN_CHECK(!attempt.generation.has_value());
  HUNDUN_CHECK(!attempt.allocation_identity.has_value());
  const auto &layout = record.identities[field_ids.size() + 1U];
  HUNDUN_CHECK(layout.subject_id == "layout.cells");
  HUNDUN_CHECK(layout.layout_fingerprint.has_value());
  HUNDUN_CHECK(*layout.layout_fingerprint == expected_layout);
  HUNDUN_CHECK(!layout.revision.has_value());
  HUNDUN_CHECK(!layout.generation.has_value());
  HUNDUN_CHECK(!layout.allocation_identity.has_value());
  const auto &finalization = record.identities[field_ids.size() + 2U];
  HUNDUN_CHECK(finalization.subject_id == "material.finalization");
  HUNDUN_CHECK(!finalization.layout_fingerprint.has_value());
  HUNDUN_CHECK(finalization.revision == report.finalization_identity());
  HUNDUN_CHECK(!finalization.generation.has_value());
  HUNDUN_CHECK(!finalization.allocation_identity.has_value());
}

std::vector<std::uint64_t>
owned_global_ids_for_rank(hundun::runtime::Int3 extent, int ranks, int rank) {
  const int quotient = extent.x / ranks;
  const int remainder = extent.x % ranks;
  const int begin = rank * quotient + std::min(rank, remainder);
  const int width = quotient + (rank < remainder ? 1 : 0);
  std::vector<std::uint64_t> result;
  result.reserve(static_cast<std::size_t>(width * extent.y * extent.z));
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = begin; i < begin + width; ++i) {
        result.push_back(
            static_cast<std::uint64_t>((k * extent.y + j) * extent.x + i));
      }
    }
  }
  return result;
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
  const std::string expected_owned_layout = owned_layout(box);
  const std::string expected_global_layout = global_layout(extent);
  HUNDUN_CHECK(source.owned_cell_layout_fingerprint() == expected_owned_layout);
  HUNDUN_CHECK(source.global_cell_layout_fingerprint() ==
               expected_global_layout);
  std::vector<std::uint64_t> local_global_ids;
  local_global_ids.reserve(topology.owned_cell_count());
  for (std::size_t cell_id = 0; cell_id < topology.owned_cell_count();
       ++cell_id)
    local_global_ids.push_back(topology.global_cell_id(cell_id));
  std::vector<std::uint64_t> all_global_ids;
  all_global_ids.reserve(
      static_cast<std::size_t>(extent.x * extent.y * extent.z));
  for (std::uint64_t global_id = 0U;
       global_id < static_cast<std::uint64_t>(extent.x * extent.y * extent.z);
       ++global_id)
    all_global_ids.push_back(global_id);
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
    require_identities(sink.value, ids, expected_owned_layout, report);
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

  for (const std::size_t budget : {1U, 256U}) {
    for (const bool narrowed : {false, true}) {
      diagnostics::DiagnosticRequest request;
      request.level = diagnostics::DiagnosticLevel::bounded_state_sample;
      request.scope = diagnostics::DiagnosticScope::local;
      request.frame = {mpi.rank(), 0U, -0.0, "material.finalized-trial"};
      request.sample_budget = budget;
      if (narrowed)
        request.selected_fields = {"rho_h"};
      CaptureSink sink;
      diagnostics::collect_diagnostics(source, request, sink);
      const auto expected =
          expected_samples(ids, local_global_ids, narrowed, budget);
      require_expected_samples(sink.value, expected);
      require_identities(sink.value, ids, expected_owned_layout, report);
      if (budget == 1U && !narrowed) {
        auto tuple_mutant = expected.samples;
        tuple_mutant.front().global_id += 1U;
        HUNDUN_CHECK(!exact_samples_equal(expected.samples, tuple_mutant));

        std::vector<diagnostics::DiagnosticSample> per_field_budget_mutant;
        const auto all_local_samples = expected_samples(
            ids, local_global_ids, false, ids.size() * local_global_ids.size());
        for (const auto field : ids) {
          const auto first = std::find_if(
              all_local_samples.samples.begin(),
              all_local_samples.samples.end(),
              [&](const auto &sample) { return sample.field_id == field; });
          HUNDUN_CHECK(first != all_local_samples.samples.end());
          per_field_budget_mutant.push_back(*first);
        }
        std::sort(
            per_field_budget_mutant.begin(), per_field_budget_mutant.end(),
            [](const auto &left, const auto &right) {
              return std::tie(left.field_id, left.global_id, left.component) <
                     std::tie(right.field_id, right.global_id, right.component);
            });
        HUNDUN_CHECK(
            !exact_samples_equal(expected.samples, per_field_budget_mutant));
      }
    }
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
  require_identities(collective_sink.value, ids, expected_global_layout,
                     report);
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
      const auto expected =
          expected_samples(ids, all_global_ids, narrowed, budget);
      require_expected_samples(sampled_sink.value, expected);
      require_identities(sampled_sink.value, ids, expected_global_layout,
                         report);
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
      if (budget == 1U && !narrowed && mpi.size() > 1) {
        std::vector<diagnostics::DiagnosticSample> per_rank_budget_mutant;
        for (int rank = 0; rank < mpi.size(); ++rank) {
          const auto rank_expected = expected_samples(
              ids, owned_global_ids_for_rank(extent, mpi.size(), rank), false,
              budget);
          per_rank_budget_mutant.insert(per_rank_budget_mutant.end(),
                                        rank_expected.samples.begin(),
                                        rank_expected.samples.end());
        }
        std::sort(
            per_rank_budget_mutant.begin(), per_rank_budget_mutant.end(),
            [](const auto &left, const auto &right) {
              return std::tie(left.field_id, left.global_id, left.component) <
                     std::tie(right.field_id, right.global_id, right.component);
            });
        HUNDUN_CHECK(
            !exact_samples_equal(expected.samples, per_rank_budget_mutant));
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
    require_identities(level_sink.value, ids, expected_global_layout, report);
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

  runtime::FieldAccessPlan checked_read_access(registry);
  constexpr runtime::PhaseId checked_read_phase = 1919U;
  constexpr runtime::ActorId checked_read_actor = 1919U;
  checked_read_access.declare_access(checked_read_phase, checked_read_actor,
                                     fields.density, runtime::AccessMode::read);
  checked_read_access.freeze();
  auto checked_density = state.trial_layer().acquire_read<double>(
      checked_read_access, checked_read_phase, checked_read_actor,
      fields.density);
  const double checked_density_before = checked_density(0, 0, 0, 0);

  auto require_injected_failure =
      [&](const diagnostics::DiagnosticRequest &request, auto arm,
          diagnostics::DiagnosticFailureClass klass, std::string_view code,
          int expected_rank,
          std::optional<diagnostics::test::MaterialDiagnosticRawCollectivePoint>
              expected_last = std::nullopt,
          std::optional<std::uint64_t> expected_collectives = std::nullopt) {
        const auto history_before = state.snapshot(flow::FlowLayer::history);
        const auto committed_before =
            state.snapshot(flow::FlowLayer::committed);
        const auto trial_before = state.snapshot(flow::FlowLayer::trial);
        const auto metadata_before_each = state.metadata();
        const auto counters_before_each = mpi.fp64_reduction_counters();
        const auto attempt_before = report.attempt_identity();
        const auto finalization_before = report.finalization_identity();
        const auto source_attempt_before = source.report().attempt_identity();
        const auto source_finalization_before =
            source.report().finalization_identity();
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            reset();
        arm();
        CaptureSink injected_sink;
        bool injected_rejected = false;
        try {
          diagnostics::collect_diagnostics(source, mpi, request, injected_sink);
        } catch (const diagnostics::DiagnosticCollectionError &error) {
          injected_rejected = true;
          HUNDUN_CHECK(error.classification() == klass);
          HUNDUN_CHECK(error.code() == code);
          HUNDUN_CHECK(error.lowest_failing_rank() == expected_rank);
        }
        HUNDUN_CHECK(injected_rejected);
        HUNDUN_CHECK(injected_sink.calls == 0U);
        const auto work = diagnostics::test::
            MaterialDensityTransportDiagnosticsTestAccess::work_counts();
        if (expected_last)
          HUNDUN_CHECK(work.last_raw_collective == *expected_last);
        if (expected_collectives)
          HUNDUN_CHECK(work.raw_collectives == *expected_collectives);
        HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
            history_before, state.snapshot(flow::FlowLayer::history)));
        HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
            committed_before, state.snapshot(flow::FlowLayer::committed)));
        HUNDUN_CHECK(test::flow_layer_values_bitwise_equal(
            trial_before, state.snapshot(flow::FlowLayer::trial)));
        HUNDUN_CHECK(test::accepted_step_metadata_bitwise_equal(
            metadata_before_each, state.metadata()));
        HUNDUN_CHECK(mpi.fp64_reduction_counters().collective_calls ==
                     counters_before_each.collective_calls);
        HUNDUN_CHECK(report.attempt_identity() == attempt_before);
        HUNDUN_CHECK(report.finalization_identity() == finalization_before);
        HUNDUN_CHECK(source.report().attempt_identity() ==
                     source_attempt_before);
        HUNDUN_CHECK(source.report().finalization_identity() ==
                     source_finalization_before);
        HUNDUN_CHECK(checked_density(0, 0, 0, 0) == checked_density_before);
      };
  const auto first_global_id = source.global_cell_id(0U);
  require_injected_failure(
      collective,
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            override_global_id(1U, first_global_id, 0);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.duplicate-owned-sample", 0);
  const auto outside_global_id = source.global_cell_id(9U);
  require_injected_failure(
      collective,
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            override_global_id(10U, outside_global_id, 0);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.duplicate-owned-sample", 0);
  const int injection_rank = mpi.size() > 1 ? 1 : 0;
  require_injected_failure(
      collective,
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_provider_failure(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.injected-provider", injection_rank);
  require_injected_failure(
      collective,
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_record_failure(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.injected-record", injection_rank);
  require_injected_failure(
      collective,
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_request_size_overflow(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::invalid_request,
      "material.diagnostics.request-size", injection_rank);
  require_injected_failure(
      collective,
      [&] {
        diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
            inject_sample_wire_overflow(injection_rank);
      },
      diagnostics::DiagnosticFailureClass::layout,
      "material.diagnostics.sample-wire-size", injection_rank,
      diagnostics::test::MaterialDiagnosticRawCollectivePoint::
          sample_size_exchange,
      21U);

  if (mpi.size() > 1) {
    std::vector<std::uint64_t> cumulative_counts(
        static_cast<std::size_t>(mpi.size()), 0U);
    cumulative_counts[0] =
        static_cast<std::uint64_t>(std::numeric_limits<int>::max());
    cumulative_counts[1] = 1U;
    require_injected_failure(
        collective,
        [&] {
          diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
              override_reported_sample_wire_bytes(cumulative_counts.data(),
                                                  cumulative_counts.size());
        },
        diagnostics::DiagnosticFailureClass::layout,
        "material.diagnostics.sample-wire-size", 1,
        diagnostics::test::MaterialDiagnosticRawCollectivePoint::
            sample_size_exchange,
        21U);
  }

  using AllocationPoint = diagnostics::test::MaterialDiagnosticAllocationPoint;
  using RawPoint = diagnostics::test::MaterialDiagnosticRawCollectivePoint;
  struct AllocationFailureCase final {
    AllocationPoint point;
    RawPoint raw_point;
    std::uint64_t raw_collectives;
    bool summary_request;
    bool required_distributed;
  };
  const std::array allocation_failures{
      AllocationFailureCase{AllocationPoint::summary_gather,
                            RawPoint::preparation_summary_gather, 9U, true,
                            true},
      AllocationFailureCase{AllocationPoint::transport_totals,
                            RawPoint::preparation_transport_totals, 11U, true,
                            true},
      AllocationFailureCase{AllocationPoint::owned_id_local,
                            RawPoint::preparation_owned_id_local, 12U, false,
                            false},
      AllocationFailureCase{AllocationPoint::ownership_counts,
                            RawPoint::preparation_ownership_counts, 14U, false,
                            true},
      AllocationFailureCase{AllocationPoint::ownership_gather,
                            RawPoint::preparation_ownership_gather, 16U, false,
                            true},
      AllocationFailureCase{AllocationPoint::eligible_counts,
                            RawPoint::preparation_eligible_counts, 18U, false,
                            false},
      AllocationFailureCase{
          AllocationPoint::local_sample_wire_and_size_counts,
          RawPoint::preparation_local_sample_wire_and_size_counts, 20U, false,
          true},
      AllocationFailureCase{AllocationPoint::sample_exchange_buffers,
                            RawPoint::preparation_sample_exchange_buffers, 22U,
                            false, true},
      AllocationFailureCase{AllocationPoint::decoded_and_retained_samples,
                            RawPoint::preparation_decoded_and_retained_samples,
                            24U, false, true}};
  for (const auto &failure : allocation_failures) {
    if (mpi.size() > 1 && !failure.required_distributed)
      continue;
    const auto &request = failure.summary_request ? sink_request : collective;
    require_injected_failure(
        request,
        [&] {
          diagnostics::test::MaterialDensityTransportDiagnosticsTestAccess::
              inject_allocation_failure(failure.point, injection_rank);
        },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "material.diagnostics.aggregation-preparation", injection_rank,
        failure.raw_point, failure.raw_collectives);
  }

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
