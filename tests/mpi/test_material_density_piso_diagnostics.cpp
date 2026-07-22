// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/diagnostics/material_density_piso_diagnostics.hpp"
#include "diagnostics/src/material_density_piso_diagnostics_test_access.hpp"
#include "flow/src/material_density_piso_test_access.hpp"
#include "flow/src/material_density_transport_test_access.hpp"
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
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
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
  config.scalars.push_back({"scalar_0", 0.0});
  config.scalars.push_back({"scalar_1", 0.0});
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

class DiagnosticFaultReset final {
public:
  ~DiagnosticFaultReset() noexcept {
    hundun::diagnostics::test::MaterialDensityPisoDiagnosticTestAccess::reset();
  }
};

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

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::uint64_t fnv1a(std::string_view text) noexcept {
  std::uint64_t result = UINT64_C(14695981039346656037);
  for (const char character : text) {
    const auto value = static_cast<unsigned char>(character);
    result ^= value;
    result *= UINT64_C(1099511628211);
  }
  return result;
}

bool ends_with(std::string_view value, std::string_view suffix) noexcept {
  return value.size() >= suffix.size() &&
         value.substr(value.size() - suffix.size()) == suffix;
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    using DiagnosticAccess = hundun::diagnostics::test::
        MaterialDensityPisoDiagnosticTestAccess;
    using SampleWireItem =
        hundun::diagnostics::test::MaterialDensityPisoSampleWireItem;
    using ValueStatus = hundun::diagnostics::DiagnosticValueStatus;
    const std::vector<SampleWireItem> wire_items{
        {3U, UINT64_C(0x0102030405060708), 0U,
         UINT64_C(0x3ff0000000000000), ValueStatus::finite},
        {1U, UINT64_C(9), 2U, UINT64_C(0x7ff8000000000001),
         ValueStatus::quiet_nan}};
    const std::vector<unsigned char> expected_wire{
        0x01U, 0x00U, 0x00U, 0x00U,
        0x02U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U,
        0x03U, 0x00U, 0x00U, 0x00U,
        0x08U, 0x07U, 0x06U, 0x05U, 0x04U, 0x03U, 0x02U, 0x01U,
        0x00U, 0x00U, 0x00U, 0x00U,
        0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0xf0U, 0x3fU,
        0x00U,
        0x01U, 0x00U, 0x00U, 0x00U,
        0x09U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U,
        0x02U, 0x00U, 0x00U, 0x00U,
        0x01U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0xf8U, 0x7fU,
        0x03U};
    const auto encoded_wire = DiagnosticAccess::encode_sample_wire(wire_items);
    HUNDUN_CHECK(encoded_wire == expected_wire);
    const auto decoded_wire = DiagnosticAccess::decode_sample_wire(encoded_wire);
    HUNDUN_CHECK(decoded_wire.size() == wire_items.size());
    for (std::size_t item = 0; item < wire_items.size(); ++item) {
      HUNDUN_CHECK(decoded_wire[item].field == wire_items[item].field);
      HUNDUN_CHECK(decoded_wire[item].global_id == wire_items[item].global_id);
      HUNDUN_CHECK(decoded_wire[item].component == wire_items[item].component);
      HUNDUN_CHECK(decoded_wire[item].value_bits == wire_items[item].value_bits);
      HUNDUN_CHECK(decoded_wire[item].status == wire_items[item].status);
    }
    const auto empty_wire = DiagnosticAccess::encode_sample_wire({});
    const std::vector<unsigned char> expected_empty_wire{
        0x01U, 0x00U, 0x00U, 0x00U,
        0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U, 0x00U};
    HUNDUN_CHECK(empty_wire == expected_empty_wire);
    HUNDUN_CHECK(DiagnosticAccess::decode_sample_wire(empty_wire).empty());

    const auto expect_wire_rejected = [](std::vector<unsigned char> bytes) {
      bool rejected = false;
      try {
        static_cast<void>(DiagnosticAccess::decode_sample_wire(bytes));
      } catch (const std::exception &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
    };
    auto truncated_wire = expected_wire;
    truncated_wire.pop_back();
    expect_wire_rejected(std::move(truncated_wire));
    auto trailing_wire = expected_wire;
    trailing_wire.push_back(0U);
    expect_wire_rejected(std::move(trailing_wire));
    auto invalid_schema_wire = expected_wire;
    invalid_schema_wire[0] = 2U;
    expect_wire_rejected(std::move(invalid_schema_wire));
    auto impossible_count_wire = expected_empty_wire;
    std::fill(impossible_count_wire.begin() + 4, impossible_count_wire.end(),
              0xffU);
    expect_wire_rejected(std::move(impossible_count_wire));
    auto invalid_field_wire = expected_wire;
    invalid_field_wire[12] = 5U;
    expect_wire_rejected(std::move(invalid_field_wire));
    auto invalid_component_wire = expected_wire;
    invalid_component_wire[24] = 1U;
    expect_wire_rejected(std::move(invalid_component_wire));
    auto invalid_status_wire = expected_wire;
    invalid_status_wire[36] =
        static_cast<unsigned char>(ValueStatus::unavailable);
    expect_wire_rejected(std::move(invalid_status_wire));
    auto inconsistent_status_wire = expected_wire;
    inconsistent_status_wire[35] = 0x7fU;
    expect_wire_rejected(std::move(inconsistent_status_wire));

    auto changed_items = wire_items;
    changed_items[0].value_bits = UINT64_C(0x4000000000000000);
    const auto changed_wire = DiagnosticAccess::encode_sample_wire(changed_items);
    HUNDUN_CHECK(changed_wire != encoded_wire);
    const auto changed_decoded =
        DiagnosticAccess::decode_sample_wire(changed_wire);
    HUNDUN_CHECK(changed_decoded.size() == changed_items.size());
    HUNDUN_CHECK(changed_decoded[0].value_bits == changed_items[0].value_bits);

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
    const auto rho_scalar_0 =
        registry.declare_field(cell("rho_scalar_0", "kg/m3", 1U, true));
    const auto rho_scalar_1 =
        registry.declare_field(cell("rho_scalar_1", "kg/m3", 1U, true));
    fields.transported_cell_fields = {rho_h, rho_scalar_0, rho_scalar_1};
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
          std::vector<double>(topology.owned_cell_count(), 3.0),
          std::vector<double>(topology.owned_cell_count(), 0.25),
          std::vector<double>(topology.owned_cell_count(), 0.75)};
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
    specification.scalar_densities = {rho_scalar_0, rho_scalar_1};
    specification.scalar_diffusivities_kg_per_m_s = {0.0, 0.0};
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
    const auto diagnostic_state = capture(state);

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
    std::optional<hundun::diagnostics::DiagnosticRecord> summary_record;
    std::optional<hundun::diagnostics::DiagnosticRecord> invariant_record;
    std::optional<hundun::diagnostics::DiagnosticRecord> counter_record;
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
      check_equal(diagnostic_state, state);
      HUNDUN_CHECK(hundun::flow::test::MaterialDensityPisoTestAccess::
                       report_authenticated(report));
      if (level == 0U)
      {
        constexpr std::array<std::string_view, 24> ids{
            "attempt.dt", "continuity.residual", "density.maximum",
            "density.minimum", "mass.relative-defect", "mass.total",
            "material.density.residual", "material.mass.relative-defect",
            "momentum.x.relative-defect", "momentum.x.residual",
            "momentum.x.total", "momentum.y.relative-defect",
            "momentum.y.residual", "momentum.y.total",
            "momentum.z.relative-defect", "momentum.z.residual",
            "momentum.z.total", "pressure.residual",
            "transport.s00000000000000000000.relative-defect",
            "transport.s00000000000000000000.residual",
            "transport.s00000000000000000001.relative-defect",
            "transport.s00000000000000000001.residual",
            "transport.s00000000000000000002.relative-defect",
            "transport.s00000000000000000002.residual"};
        HUNDUN_CHECK(first.records[0].metrics.size() == ids.size());
        for (std::size_t index = 0; index < ids.size(); ++index) {
          const auto &metric = first.records[0].metrics[index];
          HUNDUN_CHECK(metric.id == ids[index]);
          HUNDUN_CHECK(metric.value.status ==
                       hundun::diagnostics::DiagnosticValueStatus::finite);
          const bool state_summary = metric.id == "attempt.dt" ||
                                     metric.id == "density.maximum" ||
                                     metric.id == "density.minimum";
          const bool residual = metric.id.find("residual") != std::string::npos;
          HUNDUN_CHECK(metric.kind ==
                       (state_summary
                            ? hundun::diagnostics::DiagnosticMetricKind::
                                  state_summary
                            : residual
                                  ? hundun::diagnostics::DiagnosticMetricKind::
                                        residual
                                  : hundun::diagnostics::DiagnosticMetricKind::
                                        conservation));
          const std::string_view expected_unit =
              metric.id == "attempt.dt"
                  ? "s"
                  : (metric.id == "density.maximum" ||
                     metric.id == "density.minimum")
                        ? "kg/m3"
                        : metric.id == "mass.total"
                              ? "kg"
                              : metric.id.find("momentum.") == 0U &&
                                        ends_with(metric.id, ".total")
                                    ? "kg*m/s"
                                    : "1";
          HUNDUN_CHECK(metric.unit == expected_unit);
          const double expected_value =
              metric.id == "attempt.dt"
                  ? 0.01
                  : (metric.id == "density.maximum" ||
                     metric.id == "density.minimum" ||
                     metric.id == "mass.total")
                        ? 1.0
                        : 0.0;
          HUNDUN_CHECK(
              metric.value.bits ==
              bits(metric.id == "mass.total"
                       ? expected_value / static_cast<double>(mpi.size())
                       : expected_value));
        }
        summary_record = first.records[0];
      }
      if (level == 1U) {
        constexpr std::array<std::string_view, 38> ids{
            "continuity.residual", "continuity.residual.available",
            "density.positive", "flux.final-corrected",
            "mass.relative-defect", "mass.relative-defect.available",
            "material.density.residual",
            "material.density.residual.available", "material.finalized",
            "material.mass.relative-defect",
            "material.mass.relative-defect.available",
            "momentum.x.relative-defect",
            "momentum.x.relative-defect.available", "momentum.x.residual",
            "momentum.x.residual.available", "momentum.y.relative-defect",
            "momentum.y.relative-defect.available", "momentum.y.residual",
            "momentum.y.residual.available", "momentum.z.relative-defect",
            "momentum.z.relative-defect.available", "momentum.z.residual",
            "momentum.z.residual.available", "pressure.correctors",
            "pressure.residual", "pressure.residual.available",
            "transport.s00000000000000000000.relative-defect",
            "transport.s00000000000000000000.relative-defect.available",
            "transport.s00000000000000000000.residual",
            "transport.s00000000000000000000.residual.available",
            "transport.s00000000000000000001.relative-defect",
            "transport.s00000000000000000001.relative-defect.available",
            "transport.s00000000000000000001.residual",
            "transport.s00000000000000000001.residual.available",
            "transport.s00000000000000000002.relative-defect",
            "transport.s00000000000000000002.relative-defect.available",
            "transport.s00000000000000000002.residual",
            "transport.s00000000000000000002.residual.available"};
        HUNDUN_CHECK(first.records[0].invariants.size() == ids.size());
        for (std::size_t index = 0; index < ids.size(); ++index) {
          const auto &invariant = first.records[0].invariants[index];
          HUNDUN_CHECK(invariant.id == ids[index]);
          HUNDUN_CHECK(invariant.observed.status ==
                       hundun::diagnostics::DiagnosticValueStatus::finite);
          HUNDUN_CHECK(invariant.limit.status ==
                       (invariant.id == "density.positive"
                            ? hundun::diagnostics::DiagnosticValueStatus::
                                  unavailable
                            : hundun::diagnostics::DiagnosticValueStatus::
                                  finite));
          HUNDUN_CHECK(invariant.passed);
          const bool equality = ends_with(invariant.id, ".available") ||
                                invariant.id == "flux.final-corrected" ||
                                invariant.id == "material.finalized" ||
                                invariant.id == "pressure.correctors";
          HUNDUN_CHECK(invariant.relation ==
                       (invariant.id == "density.positive"
                            ? hundun::diagnostics::InvariantRelation::positive
                            : equality
                                  ? hundun::diagnostics::InvariantRelation::equal
                                  : hundun::diagnostics::InvariantRelation::
                                        less_equal));
          HUNDUN_CHECK(invariant.unit ==
                       (invariant.id == "density.positive"
                            ? "kg/m3"
                            : invariant.id == "pressure.correctors" ? "count"
                                                                    : "1"));
        }
        invariant_record = first.records[0];
      }
      if (level == 2U) {
        constexpr std::array<std::string_view, 12> ids{
            "attempt.identity", "material.fields",
            "material.finalization.identity", "pressure.correctors",
            "solver.momentum.iterations", "solver.momentum.matvec",
            "solver.momentum.preconditioner", "solver.momentum.reductions",
            "solver.pressure.iterations", "solver.pressure.matvec",
            "solver.pressure.preconditioner", "solver.pressure.reductions"};
        HUNDUN_CHECK(first.records[0].counters.size() == ids.size());
        std::uint64_t momentum_iterations{};
        std::uint64_t momentum_matvec{};
        std::uint64_t momentum_preconditioner{};
        std::uint64_t momentum_reductions{};
        for (const auto &solve : report.flow().momentum.components) {
          momentum_iterations += solve.iterations;
          momentum_matvec += solve.matvec_count;
          momentum_preconditioner += solve.preconditioner_apply_count;
          momentum_reductions += solve.global_reduction_count;
        }
        std::uint64_t pressure_iterations{};
        std::uint64_t pressure_matvec{};
        std::uint64_t pressure_preconditioner_count{};
        std::uint64_t pressure_reductions{};
        for (const auto &solve : report.flow().pressure) {
          pressure_iterations += solve.iterations;
          pressure_matvec += solve.matvec_count;
          pressure_preconditioner_count += solve.preconditioner_apply_count;
          pressure_reductions += solve.global_reduction_count;
        }
        const std::array expected_values{
            report.attempt_identity(), UINT64_C(3),
            report.material_report().finalization_identity(), UINT64_C(2),
            momentum_iterations, momentum_matvec, momentum_preconditioner,
            momentum_reductions, pressure_iterations, pressure_matvec,
            pressure_preconditioner_count, pressure_reductions};
        for (std::size_t index = 0; index < ids.size(); ++index) {
          HUNDUN_CHECK(first.records[0].counters[index].id == ids[index]);
          HUNDUN_CHECK(first.records[0].counters[index].unit == "count");
          HUNDUN_CHECK(first.records[0].counters[index].value ==
                       expected_values[index]);
        }
        counter_record = first.records[0];
      }
      if (level == 3U) {
        HUNDUN_CHECK(first.records[0].samples.size() == 7U);
        HUNDUN_CHECK(first.records[0].sample_budget == 7U);
        HUNDUN_CHECK(std::is_sorted(
            first.records[0].samples.begin(), first.records[0].samples.end(),
            [](const auto &left, const auto &right) {
              return std::tie(left.field_id, left.global_id, left.component) <
                     std::tie(right.field_id, right.global_id,
                             right.component);
            }));
        for (std::size_t index = 0; index < 7U; ++index) {
          const auto &sample = first.records[0].samples[index];
          HUNDUN_CHECK(sample.field_id == "face_mass_flux");
          HUNDUN_CHECK(sample.global_id == source.field_global_id(0U, index));
          HUNDUN_CHECK(sample.component == 0U);
          HUNDUN_CHECK(sample.unit == "kg/s");
          HUNDUN_CHECK(sample.value.status ==
                       hundun::diagnostics::DiagnosticValueStatus::finite);
          HUNDUN_CHECK(sample.value.bits ==
                       bits(source.field_value(0U, index, 0U)));
        }
      }

      constexpr std::array<std::string_view, 9> identity_ids{
          "field.face_mass_flux", "field.face_velocity", "field.pi",
          "field.rho", "field.velocity", "flow.attempt", "layout.cells",
          "layout.faces", "material.finalization"};
      HUNDUN_CHECK(first.records[0].identities.size() == identity_ids.size());
      for (std::size_t index = 0; index < identity_ids.size(); ++index)
        HUNDUN_CHECK(first.records[0].identities[index].subject_id ==
                     identity_ids[index]);
      for (std::size_t index = 0; index < 5U; ++index) {
        HUNDUN_CHECK(first.records[0].identities[index].layout_fingerprint);
        HUNDUN_CHECK(!first.records[0].identities[index].revision);
      }
      HUNDUN_CHECK(first.records[0].identities[5].revision ==
                   report.attempt_identity());
      HUNDUN_CHECK(first.records[0].identities[8].revision ==
                   report.material_report().finalization_identity());
      for (const auto &identity : first.records[0].identities) {
        HUNDUN_CHECK(!identity.generation);
        HUNDUN_CHECK(!identity.allocation_identity);
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
      check_equal(diagnostic_state, state);
    }

    const auto expect_record_mutation_rejected = [&](const auto &original,
                                                      auto mutate) {
      auto changed = original;
      mutate(changed);
      bool rejected = false;
      try {
        hundun::diagnostics::validate(changed);
      } catch (const std::exception &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
    };
    HUNDUN_CHECK(summary_record && invariant_record && counter_record);
    expect_record_mutation_rejected(*summary_record, [](auto &record) {
      std::swap(record.metrics[0], record.metrics[1]);
    });
    expect_record_mutation_rejected(*summary_record, [](auto &record) {
      record.metrics[1] = record.metrics[0];
    });
    expect_record_mutation_rejected(*invariant_record, [](auto &record) {
      std::swap(record.invariants[0], record.invariants[1]);
    });
    expect_record_mutation_rejected(*invariant_record, [](auto &record) {
      record.invariants[1] = record.invariants[0];
    });
    expect_record_mutation_rejected(*counter_record, [](auto &record) {
      std::swap(record.counters[0], record.counters[1]);
    });
    expect_record_mutation_rejected(*counter_record, [](auto &record) {
      record.counters[1] = record.counters[0];
    });

    const auto check_budget = [&](hundun::diagnostics::DiagnosticScope scope,
                                  std::size_t budget,
                                  std::uint64_t eligible) {
      auto selected_request = request(
          hundun::diagnostics::DiagnosticLevel::bounded_state_sample, scope);
      selected_request.selected_fields = {"rho"};
      selected_request.sample_budget = budget;
      Sink selected_sink;
      if (scope == hundun::diagnostics::DiagnosticScope::local)
        hundun::diagnostics::collect_diagnostics(source, selected_request,
                                                 selected_sink);
      else
        hundun::diagnostics::collect_diagnostics(
            source, mpi, selected_request, selected_sink);
      const auto &selected_record = selected_sink.records[0];
      HUNDUN_CHECK(selected_record.eligible_sample_count == eligible);
      HUNDUN_CHECK(selected_record.samples.size() ==
                   std::min<std::uint64_t>(budget, eligible));
      HUNDUN_CHECK(selected_record.samples_truncated ==
                   (eligible > selected_record.samples.size()));
      for (std::size_t index = 0; index < selected_record.samples.size();
           ++index) {
        HUNDUN_CHECK(selected_record.samples[index].field_id == "rho");
        HUNDUN_CHECK(selected_record.samples[index].component == 0U);
        const std::uint64_t expected_global =
            scope == hundun::diagnostics::DiagnosticScope::local
                ? source.field_global_id(3U, index)
                : static_cast<std::uint64_t>(index);
        HUNDUN_CHECK(selected_record.samples[index].global_id ==
                     expected_global);
      }
      check_equal(diagnostic_state, state);
    };
    const std::uint64_t local_eligible = source.owned_cell_count();
    for (const auto budget :
         {std::size_t{1},
          static_cast<std::size_t>(local_eligible / 2U),
          static_cast<std::size_t>(local_eligible),
          static_cast<std::size_t>(local_eligible + 1U)})
      check_budget(hundun::diagnostics::DiagnosticScope::local, budget,
                   local_eligible);
    constexpr std::uint64_t global_eligible = 8U * 4U * 4U;
    for (const auto budget : {std::size_t{1},
                              std::size_t{63}, std::size_t{128},
                              std::size_t{129}})
      check_budget(hundun::diagnostics::DiagnosticScope::collective, budget,
                   global_eligible);
    for (const auto scope : {hundun::diagnostics::DiagnosticScope::local,
                             hundun::diagnostics::DiagnosticScope::collective}) {
      auto zero_budget = request(
          hundun::diagnostics::DiagnosticLevel::bounded_state_sample, scope);
      zero_budget.selected_fields = {"rho"};
      zero_budget.sample_budget = 0U;
      Sink zero_sink;
      bool rejected = false;
      try {
        if (scope == hundun::diagnostics::DiagnosticScope::local)
          hundun::diagnostics::collect_diagnostics(source, zero_budget,
                                                   zero_sink);
        else
          hundun::diagnostics::collect_diagnostics(source, mpi, zero_budget,
                                                   zero_sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        rejected = error.classification() ==
                       hundun::diagnostics::DiagnosticFailureClass::
                           invalid_request &&
                   error.code() == "flow.diagnostics.frame";
      }
      HUNDUN_CHECK(rejected && zero_sink.calls == 0U);
    }

    const std::array invalid_selections{
        std::vector<std::string_view>{"velocity", "rho"},
        std::vector<std::string_view>{"rho", "rho"}};
    for (const auto &invalid_selection : invalid_selections) {
      Sink sink;
      auto invalid = request(
          hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
          hundun::diagnostics::DiagnosticScope::local);
      invalid.selected_fields = invalid_selection;
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
      check_equal(diagnostic_state, state);
    }

    {
      hundun::linear::JacobiPreconditioner fx(execution), fy(execution),
          fz(execution), fp(execution);
      auto fixture_flow =
          hundun::flow::FixedStepMaterialDensityFlow::create(
              decomposition, topology, geometry, boundaries, mpi, execution,
              halo, momentum_solver, {&fx, &fy, &fz}, pressure_solver, fp,
              registry, fields, specification);
      auto pre_state = make_state(1.0);
      const auto pre_before = capture(pre_state);
      const auto pre_report =
          fixture_flow.attempt(pre_state, -1.0, stencil, {}, {});
      HUNDUN_CHECK(!pre_report.material_report_available());
      HUNDUN_CHECK(pre_report.flow().reason ==
                   hundun::flow::StepFailureReason::invalid_input);
      check_equal(pre_before, pre_state);
      auto pre_source = fixture_flow.diagnostic_source(pre_state, pre_report);
      Sink pre_sink;
      hundun::diagnostics::DiagnosticRequest pre_request{
          hundun::diagnostics::DiagnosticLevel::summary,
          hundun::diagnostics::DiagnosticScope::local,
          {mpi.rank(), pre_source.committed_step(),
           pre_source.committed_time_s(), "material-density.attempt-result"},
          {}, 0U};
      hundun::diagnostics::collect_diagnostics(pre_source, pre_request,
                                               pre_sink);
      const auto pre_json =
          hundun::diagnostics::to_canonical_json(pre_sink.records[0]);
      HUNDUN_CHECK(pre_sink.records[0].failure.code == "flow.invalid-input");

      auto post_state = make_state(1.0);
      const auto post_before = capture(post_state);
      using Access = hundun::flow::test::MaterialDensityPisoTestAccess;
      Access::reset_terminal_fault();
      Access::set_terminal_fault(
          hundun::flow::test::MaterialTerminalPointForTest::
              final_pressure_entry,
          hundun::flow::test::MaterialTerminalModeForTest::returned_reliable);
      const auto post_report =
          fixture_flow.attempt(post_state, 0.0, stencil, {}, {});
      Access::reset_terminal_fault();
      HUNDUN_CHECK(post_report.material_report_available());
      HUNDUN_CHECK(post_report.flow().reason ==
                   hundun::flow::StepFailureReason::collective_operation);
      HUNDUN_CHECK(post_report.flow().lowest_failing_rank == 0);
      check_equal(post_before, post_state);
      auto post_source =
          fixture_flow.diagnostic_source(post_state, post_report);
      Sink post_sink;
      hundun::diagnostics::DiagnosticRequest post_request{
          hundun::diagnostics::DiagnosticLevel::summary,
          hundun::diagnostics::DiagnosticScope::local,
          {mpi.rank(), post_source.committed_step(),
           post_source.committed_time_s(), "material-density.attempt-result"},
          {}, 0U};
      hundun::diagnostics::collect_diagnostics(post_source, post_request,
                                               post_sink);
      const auto post_json =
          hundun::diagnostics::to_canonical_json(post_sink.records[0]);
      HUNDUN_CHECK(post_sink.records[0].failure.code ==
                   "flow.collective-operation");
      const std::size_t fixture_index =
          mpi.size() == 1 ? 0U : mpi.size() == 2 ? 1U : 2U;
      if (mpi.rank() == 0) {
        constexpr std::array<std::size_t, 3> success_sizes{5076U, 5076U,
                                                           5077U};
        constexpr std::array<std::uint64_t, 3> success_hashes{
            UINT64_C(10999181260519817470), UINT64_C(2976880533104871064),
            UINT64_C(6173095363178543331)};
        constexpr std::array<std::size_t, 3> pre_sizes{4510U, 4510U, 4511U};
        constexpr std::array<std::uint64_t, 3> pre_hashes{
            UINT64_C(8594279570245361814), UINT64_C(13666731215613088441),
            UINT64_C(129156281885931287)};
        constexpr std::array<std::size_t, 3> post_sizes{4871U, 4871U, 4872U};
        constexpr std::array<std::uint64_t, 3> post_hashes{
            UINT64_C(10519832359024824119), UINT64_C(5165974121610780774),
            UINT64_C(1326800360729308900)};
        HUNDUN_CHECK(local_json[0].size() == success_sizes[fixture_index]);
        HUNDUN_CHECK(fnv1a(local_json[0]) == success_hashes[fixture_index]);
        HUNDUN_CHECK(pre_json.size() == pre_sizes[fixture_index]);
        HUNDUN_CHECK(fnv1a(pre_json) == pre_hashes[fixture_index]);
        HUNDUN_CHECK(post_json.size() == post_sizes[fixture_index]);
        HUNDUN_CHECK(fnv1a(post_json) == post_hashes[fixture_index]);
      }

      auto nonfinite_record = *summary_record;
      nonfinite_record.metrics[0].value = hundun::diagnostics::describe_fp64(
          std::numeric_limits<double>::quiet_NaN());
      HUNDUN_CHECK(nonfinite_record.metrics[0].value.status ==
                   hundun::diagnostics::DiagnosticValueStatus::quiet_nan);
      HUNDUN_CHECK(nonfinite_record.metrics[0].value.bits ==
                   bits(std::numeric_limits<double>::quiet_NaN()));
      hundun::diagnostics::validate(nonfinite_record);
      const auto nonfinite_json =
          hundun::diagnostics::to_canonical_json(nonfinite_record);
      HUNDUN_CHECK(nonfinite_json.find("quiet_nan") != std::string::npos);
    }

    for (std::size_t mutation = 0; mutation < 4U; ++mutation) {
      Sink sink;
      auto invalid = request(hundun::diagnostics::DiagnosticLevel::summary,
                             hundun::diagnostics::DiagnosticScope::local);
      if (mutation == 0U)
        ++invalid.frame.step;
      else if (mutation == 1U)
        ++invalid.frame.rank;
      else if (mutation == 2U)
        invalid.frame.time_s = std::nextafter(invalid.frame.time_s, 1.0);
      else
        invalid.frame.phase = "material-density.wrong-phase";
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
      check_equal(diagnostic_state, state);
    }
    {
      using Fault = hundun::diagnostics::test::
          MaterialDensityPisoDiagnosticFault;
      using Access = hundun::diagnostics::test::
          MaterialDensityPisoDiagnosticTestAccess;
      const auto capture_local_json =
          [&](hundun::diagnostics::DiagnosticLevel level) {
        Sink local_sink;
        hundun::diagnostics::collect_diagnostics(
            source,
            request(level, hundun::diagnostics::DiagnosticScope::local),
            local_sink);
        HUNDUN_CHECK(local_sink.calls == 1U);
        return hundun::diagnostics::to_canonical_json(local_sink.records[0]);
      };
      const auto summary_before =
          capture_local_json(hundun::diagnostics::DiagnosticLevel::summary);
      const auto counters_before =
          capture_local_json(hundun::diagnostics::DiagnosticLevel::counters);
      const auto attempt_before = source.report().attempt_identity();
      const auto step_before = source.committed_step();
      const auto time_before = bits(source.committed_time_s());
      const auto expect_fault = [&](int failing_rank, Fault fault,
                                    hundun::diagnostics::DiagnosticLevel level,
                                    std::string_view code) {
        Access::reset();
        DiagnosticFaultReset reset_on_exit;
        if (mpi.rank() == failing_rank)
          Access::set_fault(fault);
        Sink sink;
        bool rejected = false;
        try {
          hundun::diagnostics::collect_diagnostics(
              source, mpi,
              request(level,
                      hundun::diagnostics::DiagnosticScope::collective),
              sink);
        } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
          rejected = error.classification() ==
                         hundun::diagnostics::DiagnosticFailureClass::layout &&
                     error.code() == code &&
                     error.lowest_failing_rank() == failing_rank;
        }
        Access::reset();
        HUNDUN_CHECK(rejected && sink.calls == 0U);
        check_equal(diagnostic_state, state);
        HUNDUN_CHECK(source.report().attempt_identity() == attempt_before);
        HUNDUN_CHECK(source.committed_step() == step_before);
        HUNDUN_CHECK(bits(source.committed_time_s()) == time_before);
        HUNDUN_CHECK(capture_local_json(
                         hundun::diagnostics::DiagnosticLevel::summary) ==
                     summary_before);
        HUNDUN_CHECK(capture_local_json(
                         hundun::diagnostics::DiagnosticLevel::counters) ==
                     counters_before);
      };
      std::array failing_ranks{0, mpi.size() - 1};
      for (std::size_t rank_index = 0;
           rank_index < (mpi.size() == 1 ? 1U : 2U); ++rank_index) {
        const int failing_rank = failing_ranks[rank_index];
        // Rank zero is the byte-agreement reference.  A rank-zero-only
        // perturbation is therefore observed first on rank one, whereas a
        // non-root perturbation is attributed to that rank.  Exercise the
        // seam only where its injected rank and the reported mismatch rank
        // have the same semantics; the other collective faults below cover
        // both the lowest and highest ranks.
        if (mpi.size() > 1 && failing_rank != 0)
          expect_fault(failing_rank, Fault::provider_agreement,
                       hundun::diagnostics::DiagnosticLevel::summary,
                       "flow.diagnostics.provider-agreement");
        expect_fault(failing_rank, Fault::cell_exact_cover,
                     hundun::diagnostics::DiagnosticLevel::summary,
                     "flow.diagnostics.exact-cover");
        expect_fault(failing_rank, Fault::face_exact_cover,
                     hundun::diagnostics::DiagnosticLevel::summary,
                     "flow.diagnostics.exact-cover");
        expect_fault(
            failing_rank, Fault::sample_send_preparation,
            hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
            "flow.diagnostics.sample-preparation");
        expect_fault(
            failing_rank, Fault::sample_receive_preparation,
            hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
            "flow.diagnostics.sample-preparation");
        expect_fault(
            failing_rank, Fault::sample_wire_malformed,
            hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
            "flow.diagnostics.sample-wire");
        expect_fault(failing_rank, Fault::record_validation,
                     hundun::diagnostics::DiagnosticLevel::summary,
                     "flow.diagnostics.record");
      }
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
    std::array collective_failing_ranks{0, mpi.size() - 1};
    for (std::size_t rank_index = 0;
         rank_index < (mpi.size() == 1 ? 1U : 2U); ++rank_index) {
      const int failing_rank = collective_failing_ranks[rank_index];
      {
        Sink sink;
        auto invalid = request(
            hundun::diagnostics::DiagnosticLevel::summary,
            hundun::diagnostics::DiagnosticScope::collective);
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
        sink.fail = mpi.rank() == failing_rank;
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
                     error.lowest_failing_rank() == failing_rank;
        }
        HUNDUN_CHECK(rejected && sink.calls == 1U);
      }
      check_equal(diagnostic_state, state);
    }

    {
      using TransportAccess =
          hundun::flow::test::MaterialDensityTransportTestAccess;
      hundun::linear::JacobiPreconditioner lx(execution), ly(execution),
          lz(execution), lp(execution);
      auto late_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&lx, &ly, &lz}, pressure_solver, lp, registry,
          fields, specification);
      auto late_state = make_state(1.0);
      const auto late_before = capture(late_state);
      TransportAccess::reset();
      TransportAccess::set_transport_residual(0U, 2.0e-9);
      std::optional<hundun::flow::MaterialDensityStepAttemptReport> late_report;
      try {
        late_report.emplace(
            late_flow.attempt(late_state, 0.0, stencil, {}, {}));
      } catch (...) {
        TransportAccess::reset();
        throw;
      }
      TransportAccess::reset();
      HUNDUN_CHECK(late_report->flow().reason ==
                   hundun::flow::StepFailureReason::final_transport_residual);
      HUNDUN_CHECK(late_report->flow().pressure_corrector_count == 2U);
      HUNDUN_CHECK(late_report->material_report_available());
      HUNDUN_CHECK(late_report->material_report().flux_provenance() ==
                   hundun::flow::MaterialFluxProvenance::final_corrected);
      HUNDUN_CHECK(bits(late_report->flow().final_transport_normalized_l2[0]) ==
                   bits(late_report->material_report()
                            .transport_normalized_l2()[0]));
      HUNDUN_CHECK(
          hundun::flow::test::MaterialDensityPisoTestAccess::
              report_authenticated(*late_report));
      check_equal(late_before, late_state);

      auto late_source =
          late_flow.diagnostic_source(late_state, *late_report);
      hundun::diagnostics::DiagnosticRequest late_request{
          hundun::diagnostics::DiagnosticLevel::summary,
          hundun::diagnostics::DiagnosticScope::local,
          {mpi.rank(), late_source.committed_step(),
           late_source.committed_time_s(), "material-density.attempt-result"},
          {}, 0U};
      Sink late_first;
      Sink late_second;
      hundun::diagnostics::collect_diagnostics(late_source, late_request,
                                               late_first);
      hundun::diagnostics::collect_diagnostics(late_source, late_request,
                                               late_second);
      HUNDUN_CHECK(late_first.records[0].status ==
                   hundun::diagnostics::DiagnosticStatus::failed);
      HUNDUN_CHECK(late_first.records[0].failure.classification ==
                   hundun::diagnostics::DiagnosticFailureClass::
                       non_convergence);
      HUNDUN_CHECK(late_first.records[0].failure.code ==
                   "flow.final-transport-residual");
      const auto residual = std::find_if(
          late_first.records[0].metrics.begin(),
          late_first.records[0].metrics.end(), [](const auto &metric) {
            return metric.id ==
                   "transport.s00000000000000000000.residual";
          });
      HUNDUN_CHECK(residual != late_first.records[0].metrics.end());
      HUNDUN_CHECK(residual->value.status ==
                   hundun::diagnostics::DiagnosticValueStatus::finite);
      HUNDUN_CHECK(residual->value.bits == bits(2.0e-9));
      HUNDUN_CHECK(hundun::diagnostics::to_canonical_json(
                       late_first.records[0]) ==
                   hundun::diagnostics::to_canonical_json(
                       late_second.records[0]));
      auto late_counter_request = late_request;
      late_counter_request.level =
          hundun::diagnostics::DiagnosticLevel::counters;
      Sink late_counters_first;
      Sink late_counters_second;
      hundun::diagnostics::collect_diagnostics(
          late_source, late_counter_request, late_counters_first);
      hundun::diagnostics::collect_diagnostics(
          late_source, late_counter_request, late_counters_second);
      HUNDUN_CHECK(hundun::diagnostics::to_canonical_json(
                       late_counters_first.records[0]) ==
                   hundun::diagnostics::to_canonical_json(
                       late_counters_second.records[0]));
      HUNDUN_CHECK(
          hundun::flow::test::MaterialDensityPisoTestAccess::
              report_authenticated(late_source.report()));
      check_equal(late_before, late_state);
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
