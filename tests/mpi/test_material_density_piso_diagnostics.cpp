// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/diagnostics/material_density_piso_diagnostics.hpp"
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
#include "tests/support/test_main.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

hundun::runtime::Int3 grid(int ranks) {
  return ranks == 1 ? hundun::runtime::Int3{1, 1, 1}
                    : ranks == 2 ? hundun::runtime::Int3{2, 1, 1}
                                 : hundun::runtime::Int3{2, 2, 1};
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
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name, unit, "task20-diagnostics",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, components, 2, conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name, const char *unit,
                                      std::uint32_t components) {
  return {name, unit, "task20-diagnostics",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64, components, 0, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

class Sink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &record) override {
    ++calls;
    if (fail)
      throw std::runtime_error("injected sink failure");
    records.push_back(record);
  }
  bool fail{};
  std::size_t calls{};
  std::vector<hundun::diagnostics::DiagnosticRecord> records;
};

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
    const auto make_state = [&](double rho) {
      auto state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      hundun::flow::FlowLayerValues initial;
      initial.density.assign(topology.owned_cell_count(), rho);
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
    hundun::linear::JacobiPreconditioner mx(execution), my(execution),
        mz(execution), pressure_preconditioner(execution);
    hundun::flow::MaterialDensityTransportSpec specification;
    specification.enthalpy_density = rho_h;
    auto flow = hundun::flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, specification);
    const auto stencil = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
    auto state = make_state(1.0);
    const auto report = flow.attempt(state, 0.0, stencil, {}, {});
    HUNDUN_CHECK(report.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    auto source = flow.diagnostic_source(state, report);
    const auto descriptor = hundun::diagnostics::describe_diagnostics(source);
    HUNDUN_CHECK(descriptor.module_kind ==
                 hundun::diagnostics::DiagnosticModuleKind::piso);
    HUNDUN_CHECK(descriptor.module_id ==
                 "hundun.flow.fixed_step_material_density");
    constexpr std::array<std::string_view, 5> expected_fields{
        "face_mass_flux", "face_velocity", "pi", "rho", "velocity"};
    const auto actual_fields =
        hundun::diagnostics::diagnostic_fingerprint_field_ids(source);
    HUNDUN_CHECK(std::equal(actual_fields.begin(), actual_fields.end(),
                            expected_fields.begin(), expected_fields.end()));
    const auto source_step = source.committed_step();
    const double source_time = source.committed_time_s();

    const auto request = [&](hundun::diagnostics::DiagnosticLevel level,
                             hundun::diagnostics::DiagnosticScope scope) {
      return hundun::diagnostics::DiagnosticRequest{
          level,
          scope,
          {mpi.rank(), source_step, source_time,
           "material-density.attempt-result"},
          {},
          level == hundun::diagnostics::DiagnosticLevel::bounded_state_sample
              ? 7U
              : 0U};
    };
    std::array<std::string, 4> local_json;
    for (std::size_t level = 0; level < 4U; ++level) {
      Sink first;
      auto local_request = request(
          static_cast<hundun::diagnostics::DiagnosticLevel>(level),
          hundun::diagnostics::DiagnosticScope::local);
      hundun::diagnostics::collect_diagnostics(source, local_request, first);
      HUNDUN_CHECK(first.calls == 1U && first.records.size() == 1U);
      HUNDUN_CHECK(first.records[0].state_fingerprint.algorithm ==
                   hundun::diagnostics::kStateFingerprintAlgorithmV1);
      local_json[level] =
          hundun::diagnostics::to_canonical_json(first.records[0]);
      Sink repeated;
      hundun::diagnostics::collect_diagnostics(source, local_request, repeated);
      HUNDUN_CHECK(local_json[level] ==
                   hundun::diagnostics::to_canonical_json(repeated.records[0]));
      if (level == 0U)
        HUNDUN_CHECK(first.records[0].metrics.size() == 20U);
      if (level == 1U)
        HUNDUN_CHECK(!first.records[0].invariants.empty());
      if (level == 2U)
        HUNDUN_CHECK(first.records[0].counters.size() == 12U);
      if (level == 3U) {
        HUNDUN_CHECK(first.records[0].samples.size() <= 7U);
        HUNDUN_CHECK(first.records[0].sample_budget == 7U);
      }

      Sink collective;
      auto collective_request = request(
          static_cast<hundun::diagnostics::DiagnosticLevel>(level),
          hundun::diagnostics::DiagnosticScope::collective);
      hundun::diagnostics::collect_diagnostics(source, mpi, collective_request,
                                               collective);
      HUNDUN_CHECK(collective.calls == 1U && collective.records.size() == 1U);
      HUNDUN_CHECK(collective.records[0].scope ==
                   hundun::diagnostics::DiagnosticScope::collective);
      if (level == 3U)
        HUNDUN_CHECK(collective.records[0].samples.size() <= 7U);
      Sink collective_repeated;
      hundun::diagnostics::collect_diagnostics(
          source, mpi, collective_request, collective_repeated);
      HUNDUN_CHECK(
          hundun::diagnostics::to_canonical_json(collective.records[0]) ==
          hundun::diagnostics::to_canonical_json(
              collective_repeated.records[0]));
    }

    {
      Sink sink;
      auto invalid = request(hundun::diagnostics::DiagnosticLevel::summary,
                             hundun::diagnostics::DiagnosticScope::local);
      ++invalid.frame.step;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(source, invalid, sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = error.classification() ==
                       hundun::diagnostics::DiagnosticFailureClass::
                           invalid_request &&
                   error.code() == "flow.diagnostics.frame";
      }
      HUNDUN_CHECK(rejected && sink.calls == 0U);
    }
    {
      Sink sink;
      auto unknown = request(hundun::diagnostics::DiagnosticLevel::summary,
                             hundun::diagnostics::DiagnosticScope::local);
      unknown.selected_fields = {"not-a-field"};
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(source, unknown, sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = error.code() == "flow.diagnostics.unknown-field";
      }
      HUNDUN_CHECK(rejected && sink.calls == 0U);
    }
    {
      Sink sink;
      sink.fail = true;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            source,
            request(hundun::diagnostics::DiagnosticLevel::summary,
                    hundun::diagnostics::DiagnosticScope::local),
            sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = error.classification() ==
                       hundun::diagnostics::DiagnosticFailureClass::
                           sink_failure &&
                   error.code() == "diagnostics.sink.submit";
      }
      HUNDUN_CHECK(rejected && sink.calls == 1U);
    }
    {
      Sink sink;
      auto invalid = request(hundun::diagnostics::DiagnosticLevel::summary,
                             hundun::diagnostics::DiagnosticScope::collective);
      const int failing_rank = mpi.size() - 1;
      if (mpi.rank() == failing_rank)
        ++invalid.frame.step;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(source, mpi, invalid, sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = error.classification() ==
                       hundun::diagnostics::DiagnosticFailureClass::
                           invalid_request &&
                   error.code() == "flow.diagnostics.frame" &&
                   error.lowest_failing_rank() == failing_rank;
      }
      HUNDUN_CHECK(rejected && sink.calls == 0U);
    }
    {
      Sink sink;
      sink.fail = mpi.rank() == mpi.size() - 1;
      bool rejected = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            source, mpi,
            request(hundun::diagnostics::DiagnosticLevel::summary,
                    hundun::diagnostics::DiagnosticScope::collective),
            sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = error.classification() ==
                       hundun::diagnostics::DiagnosticFailureClass::
                           sink_failure &&
                   error.code() == "diagnostics.sink.submit" &&
                   error.lowest_failing_rank() == mpi.size() - 1;
      }
      HUNDUN_CHECK(rejected && sink.calls == 1U);
    }

    auto failed_state = make_state(-1.0);
    const auto failed_report =
        flow.attempt(failed_state, 0.0, stencil, {}, {});
    HUNDUN_CHECK(failed_report.flow().reason ==
                 hundun::flow::StepFailureReason::transport_failure);
    auto failed_source = flow.diagnostic_source(failed_state, failed_report);
    Sink failed_sink;
    hundun::diagnostics::DiagnosticRequest failed_request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local,
        {mpi.rank(), failed_source.committed_step(),
         failed_source.committed_time_s(), "material-density.attempt-result"},
        {}, 0U};
    hundun::diagnostics::collect_diagnostics(failed_source, failed_request,
                                             failed_sink);
    HUNDUN_CHECK(failed_sink.records[0].status ==
                 hundun::diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(failed_sink.records[0].failure.classification ==
                 hundun::diagnostics::DiagnosticFailureClass::
                     non_positive_state);
    HUNDUN_CHECK(failed_sink.records[0].failure.code ==
                 "flow.transport-non-positive-density");
    Sink failed_collective;
    failed_request.scope = hundun::diagnostics::DiagnosticScope::collective;
    hundun::diagnostics::collect_diagnostics(
        failed_source, mpi, failed_request, failed_collective);
    HUNDUN_CHECK(failed_collective.records[0].failure.classification ==
                 hundun::diagnostics::DiagnosticFailureClass::
                     non_positive_state);
    HUNDUN_CHECK(failed_collective.records[0].failure.lowest_failing_rank == 0);

    bool stale = false;
    try {
      Sink stale_sink;
      hundun::diagnostics::collect_diagnostics(
          source,
          request(hundun::diagnostics::DiagnosticLevel::summary,
                  hundun::diagnostics::DiagnosticScope::local),
          stale_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      stale = error.code() == "flow.diagnostics.stale-source";
    }
    HUNDUN_CHECK(stale);
  });
}
