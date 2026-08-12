// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_module.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/lin_preconditioners.hpp"

#include <algorithm>
#include <limits>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

class CaptureSink final : public hundun::diagnostics::DiagnosticSink {
 public:
  void submit(const hundun::diagnostics::DiagnosticRecord& record) override {
    records.push_back(record);
  }
  std::vector<hundun::diagnostics::DiagnosticRecord> records;
};

class ThrowingSink final : public hundun::diagnostics::DiagnosticSink {
 public:
  void submit(const hundun::diagnostics::DiagnosticRecord&) override {
    throw std::runtime_error("injected sink failure");
  }
};

class QueryOnlyOperator final : public hundun::linear::LinearOperator {
 public:
  explicit QueryOnlyOperator(hundun::execution::ExecutionContext& context)
      : QueryOnlyOperator(context,
                          hundun::linear::VectorLayout(2U, {0U, 1U}),
                          hundun::linear::VectorLayout(2U, {0U, 1U})) {}

  QueryOnlyOperator(hundun::execution::ExecutionContext& context,
                    hundun::linear::VectorLayout domain,
                    hundun::linear::VectorLayout range)
      : context_(&context), domain_(std::move(domain)),
        range_(std::move(range)) {}

  hundun::linear::VectorLayout domain_layout() const override {
    return domain_;
  }
  hundun::linear::VectorLayout range_layout() const override {
    return range_;
  }
  const hundun::execution::ExecutionContext& context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return 7U; }
  hundun::execution::ExecutionEvent apply(
      hundun::execution::VectorView<const double>,
      hundun::execution::VectorView<double>) const override {
    ++apply_calls;
    return hundun::execution::ExecutionEvent::completed();
  }
  bool has_diagonal() const override { return true; }
  hundun::execution::ExecutionEvent diagonal(
      hundun::execution::VectorView<double>) const override {
    ++diagonal_calls;
    return hundun::execution::ExecutionEvent::completed();
  }

  mutable std::size_t apply_calls{};
  mutable std::size_t diagonal_calls{};

 private:
  hundun::execution::ExecutionContext* context_;
  hundun::linear::VectorLayout domain_;
  hundun::linear::VectorLayout range_;
};

hundun::diagnostics::DiagnosticRequest request(
    hundun::diagnostics::DiagnosticLevel level =
        hundun::diagnostics::DiagnosticLevel::summary) {
  return {level,
          hundun::diagnostics::DiagnosticScope::local,
          {0, 3U, 0.25, "accepted-step"},
          {},
          0U};
}
template <class Source>
hundun::diagnostics::DiagnosticRecord collect_one(
    const Source& source,
    hundun::diagnostics::DiagnosticLevel level =
        hundun::diagnostics::DiagnosticLevel::summary) {
  CaptureSink sink;
  hundun::diagnostics::collect_diagnostics(source, request(level), sink);
  HUNDUN_CHECK(sink.records.size() == 1U);
  hundun::diagnostics::validate(sink.records.front());
  const auto first =
      hundun::diagnostics::to_canonical_json(sink.records.front());
  CaptureSink repeated;
  hundun::diagnostics::collect_diagnostics(source, request(level), repeated);
  HUNDUN_CHECK(repeated.records.size() == 1U);
  HUNDUN_CHECK(
      hundun::diagnostics::to_canonical_json(repeated.records.front()) ==
      first);
  return sink.records.front();
}

template <class Source>
hundun::diagnostics::DiagnosticRecord collect_valid_empty(
    const Source& source) {
  try {
    return collect_one(source);
  } catch (const hundun::diagnostics::DiagnosticCollectionError& error) {
    throw std::runtime_error(
        "valid empty diagnostic collection rejected with " +
        std::string(error.code()));
  }
}

template <class Source>
void check_fingerprint_ids(const Source &source,
                           std::initializer_list<std::string_view> expected) {
  const auto actual =
      hundun::diagnostics::diagnostic_fingerprint_field_ids(source);
  HUNDUN_CHECK(std::is_sorted(actual.begin(), actual.end()));
  HUNDUN_CHECK(std::adjacent_find(actual.begin(), actual.end()) ==
               actual.end());
  HUNDUN_CHECK(std::vector<std::string_view>(expected) == actual);
}

} // namespace

int main() {
  return hundun::test::run([] {
    hundun::execution::CpuReferenceContext execution;
    std::string empty_collection_failures;
    const auto record_empty_failure = [&](std::string_view family,
                                          const auto& source) {
      try {
        static_cast<void>(collect_one(source));
      } catch (const hundun::diagnostics::DiagnosticCollectionError& error) {
        if (!empty_collection_failures.empty())
          empty_collection_failures += "; ";
        empty_collection_failures += std::string(family) + "=" +
                                     std::string(error.code());
      }
    };
    {
      hundun::runtime::FieldRegistry empty_registry;
      empty_registry.freeze();
      const hundun::diagnostics::FieldLayoutDiagnosticSource source{
          &empty_registry, {{0, 0, 0}, 0U}, {}};
      record_empty_failure("field-layout", source);
    }
    {
      hundun::linear::GhostedVector source(
          execution, hundun::linear::VectorLayout(0U, {}));
      record_empty_failure("ghosted-vector", source);
    }
    {
      QueryOnlyOperator source(
          execution, hundun::linear::VectorLayout(0U, {}),
          hundun::linear::VectorLayout(0U, {}));
      record_empty_failure(
          "linear-operator",
          static_cast<const hundun::linear::LinearOperator&>(source));
    }
    {
      hundun::flow::StepAttemptReport report;
      report.disposition =
          hundun::flow::StepAttemptDisposition::non_retryable_failure;
      report.reason = hundun::flow::StepFailureReason::invalid_input;
      const hundun::diagnostics::ConstantDensityPisoDiagnosticSource source{
          &report};
      record_empty_failure("constant-density-piso", source);
    }
    if (!empty_collection_failures.empty())
      throw std::runtime_error("valid empty diagnostic collections rejected: " +
                               empty_collection_failures);

    check_fingerprint_ids(execution, {"backend_identity", "backend_name",
                                      "backend_name.length", "capabilities",
                                      "ordered", "space"});
    const auto context_record = collect_one(execution);
    HUNDUN_CHECK(context_record.module_id == "hundun.execution.context");
    bool rejected = false;
    try {
      static_cast<void>(collect_one(
          execution, hundun::diagnostics::DiagnosticLevel::invariants));
    } catch (const hundun::diagnostics::DiagnosticCollectionError& error) {
      rejected =
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::capability;
    }
    HUNDUN_CHECK(rejected);

    hundun::execution::Buffer buffer(execution, 4U * sizeof(double));
    check_fingerprint_ids(buffer, {"allocation_identity", "backend_identity",
                                   "byte_size", "epoch", "space"});
    const auto allocation = buffer.allocation_identity();
    const auto epoch = buffer.epoch();
    const auto buffer_record = collect_one(buffer);
    HUNDUN_CHECK(buffer_record.module_id == "hundun.execution.buffer");
    HUNDUN_CHECK(buffer.allocation_identity() == allocation);
    HUNDUN_CHECK(buffer.epoch() == epoch);

    const auto& const_buffer = buffer;
    const auto view = const_buffer.view(1U, 2U);
    check_fingerprint_ids(view,
                          {"allocation_identity", "backend_identity",
                           "element_count", "epoch", "offset_bytes",
                           "scalar_format", "space", "stride", "writable"});
    const auto view_record = collect_one(view);
    HUNDUN_CHECK(view_record.module_id == "hundun.execution.vector_view");

    QueryOnlyOperator linear_operator(execution);
    check_fingerprint_ids(
        static_cast<const hundun::linear::LinearOperator &>(linear_operator),
        {"context.backend_identity", "context.backend_name",
         "context.backend_name.length", "context.space", "diagonal_available",
         "domain.ghost", "domain.global_id", "domain.global_id.count",
         "domain.local", "domain.owned", "range.ghost", "range.global_id",
         "range.global_id.count", "range.local", "range.owned", "revision"});
    const auto operator_record = collect_one(
        static_cast<const hundun::linear::LinearOperator &>(linear_operator));
    HUNDUN_CHECK(operator_record.module_id == "hundun.linear.operator");
    HUNDUN_CHECK(linear_operator.apply_calls == 0U);
    HUNDUN_CHECK(linear_operator.diagonal_calls == 0U);

    hundun::linear::GhostedVector empty_vector(
        execution, hundun::linear::VectorLayout(0U, {}));
    check_fingerprint_ids(
        empty_vector,
        {"allocation_identity", "backend_identity", "epoch", "ghost_count",
         "layout.global_id.count", "local_count", "owned_count", "space"});
    const auto empty_vector_record = collect_valid_empty(empty_vector);
    hundun::linear::GhostedVector one_value_vector(
        execution, hundun::linear::VectorLayout(1U, {17U}));
    check_fingerprint_ids(
        one_value_vector,
        {"allocation_identity", "backend_identity", "epoch", "ghost_count",
         "layout.global_id", "layout.global_id.count", "local_count",
         "owned_count", "space"});
    HUNDUN_CHECK(collect_one(one_value_vector).state_fingerprint.hex !=
                 empty_vector_record.state_fingerprint.hex);

    QueryOnlyOperator empty_operator(
        execution, hundun::linear::VectorLayout(0U, {}),
        hundun::linear::VectorLayout(0U, {}));
    check_fingerprint_ids(
        static_cast<const hundun::linear::LinearOperator&>(empty_operator),
        {"context.backend_identity", "context.backend_name",
         "context.backend_name.length", "context.space", "diagonal_available",
         "domain.ghost", "domain.global_id.count", "domain.local",
         "domain.owned", "range.ghost", "range.global_id.count",
         "range.local", "range.owned", "revision"});
    const auto empty_operator_record = collect_valid_empty(
        static_cast<const hundun::linear::LinearOperator&>(empty_operator));
    QueryOnlyOperator empty_domain_operator(
        execution, hundun::linear::VectorLayout(0U, {}),
        hundun::linear::VectorLayout(1U, {29U}));
    check_fingerprint_ids(
        static_cast<const hundun::linear::LinearOperator&>(
            empty_domain_operator),
        {"context.backend_identity", "context.backend_name",
         "context.backend_name.length", "context.space", "diagonal_available",
         "domain.ghost", "domain.global_id.count", "domain.local",
         "domain.owned", "range.ghost", "range.global_id",
         "range.global_id.count", "range.local", "range.owned", "revision"});
    const auto empty_domain_record = collect_valid_empty(
        static_cast<const hundun::linear::LinearOperator&>(
            empty_domain_operator));
    HUNDUN_CHECK(empty_domain_record.state_fingerprint.hex !=
                 empty_operator_record.state_fingerprint.hex);
    QueryOnlyOperator empty_range_operator(
        execution, hundun::linear::VectorLayout(1U, {31U}),
        hundun::linear::VectorLayout(0U, {}));
    check_fingerprint_ids(
        static_cast<const hundun::linear::LinearOperator&>(
            empty_range_operator),
        {"context.backend_identity", "context.backend_name",
         "context.backend_name.length", "context.space", "diagonal_available",
         "domain.ghost", "domain.global_id", "domain.global_id.count",
         "domain.local", "domain.owned", "range.ghost",
         "range.global_id.count", "range.local", "range.owned", "revision"});
    static_cast<void>(collect_valid_empty(
        static_cast<const hundun::linear::LinearOperator&>(
            empty_range_operator)));

    hundun::linear::SolveReport solve_report;
    solve_report.reason =
        hundun::linear::SolveTerminationReason::converged;
    solve_report.iterations = 4U;
    solve_report.matvec_count = 5U;
    solve_report.preconditioner_apply_count = 4U;
    solve_report.global_reduction_count = 6U;
    const hundun::diagnostics::LinearSolveDiagnosticSource solve_source{
        "momentum-x", &solve_report};
    check_fingerprint_ids(
        solve_source,
        {"final_residual", "global_reduction_count", "initial_residual",
         "instance_id", "instance_id.length", "iterations",
         "lowest_failing_rank", "matvec_count", "preconditioner_apply_count",
         "recursive_residual", "termination"});
    const auto solve_summary = collect_one(solve_source);
    const auto solve_counters =
        collect_one(solve_source,
                    hundun::diagnostics::DiagnosticLevel::counters);
    HUNDUN_CHECK(solve_summary.metrics.size() == 3U);
    HUNDUN_CHECK(solve_counters.counters.size() == 3U);
    HUNDUN_CHECK(solve_summary.state_fingerprint.hex ==
                 solve_counters.state_fingerprint.hex);
    auto mutated_solve = solve_report;
    ++mutated_solve.global_reduction_count;
    const hundun::diagnostics::LinearSolveDiagnosticSource mutated_source{
        "momentum-x", &mutated_solve};
    HUNDUN_CHECK(collect_one(mutated_source).state_fingerprint.hex !=
                 solve_summary.state_fingerprint.hex);

    const hundun::diagnostics::FlowDriverDiagnosticSource driver_source{
        hundun::config::DensityModel::material,
        3U,
        0.25,
        2U,
        hundun::flow::TimeAdvanceDisposition::committed,
        hundun::flow::StepFailureReason::none,
        -1};
    check_fingerprint_ids(driver_source,
                          {"attempt_count", "density_model", "disposition",
                           "lowest_failing_rank", "reason", "step", "time"});
    const auto driver_summary = collect_one(driver_source);
    const auto driver_counters =
        collect_one(driver_source,
                    hundun::diagnostics::DiagnosticLevel::counters);
    HUNDUN_CHECK(driver_summary.status ==
                 hundun::diagnostics::DiagnosticStatus::ok);
    HUNDUN_CHECK(driver_summary.state_fingerprint.hex ==
                 driver_counters.state_fingerprint.hex);

    hundun::runtime::FieldRegistry registry;
    const auto density = registry.declare_field(
        {"rho", "kg/m3", "flow", hundun::runtime::FunctionSpace::cell_average,
         hundun::runtime::ScalarType::float64, 1U, 2, true,
         hundun::runtime::RestartPolicy::persistent,
         hundun::runtime::OutputPolicy::always});
    const auto flux = registry.declare_field(
        {"face_mass_flux", "kg/s", "flow",
         hundun::runtime::FunctionSpace::face_value,
         hundun::runtime::ScalarType::float64, 1U, 0, true,
         hundun::runtime::RestartPolicy::persistent,
         hundun::runtime::OutputPolicy::selected});
    registry.freeze();
    const hundun::diagnostics::FieldLayoutDiagnosticSource field_source{
        &registry, {{2, 2, 2}, 36U},
        {{"density", density}, {"face-mass-flux", flux}}};
    check_fingerprint_ids(
        field_source,
        {"field.components", "field.conservative", "field.count",
         "field.ghost_width", "field.id", "field.name", "field.name.length",
         "field.output", "field.owner", "field.owner.length", "field.restart",
         "field.scalar_type", "field.space", "field.unit", "field.unit.length",
         "layout.cell_extent", "layout.face_count", "role.count", "role.field",
         "role.name", "role.name.length"});
    const auto field_record = collect_one(field_source);
    HUNDUN_CHECK(field_record.module_id == "hundun.runtime.field_layout");
    HUNDUN_CHECK(field_record.identities.size() == 4U);

    hundun::runtime::FieldRegistry empty_registry;
    empty_registry.freeze();
    const hundun::diagnostics::FieldLayoutDiagnosticSource empty_field_source{
        &empty_registry, {{0, 0, 0}, 0U}, {}};
    check_fingerprint_ids(
        empty_field_source,
        {"field.count", "layout.cell_extent", "layout.face_count",
         "role.count"});
    const auto empty_field_record = collect_valid_empty(empty_field_source);
    HUNDUN_CHECK(empty_field_record.identities.empty());
    const hundun::diagnostics::FieldLayoutDiagnosticSource empty_roles_source{
        &registry, {{2, 2, 2}, 36U}, {}};
    check_fingerprint_ids(
        empty_roles_source,
        {"field.components", "field.conservative", "field.count",
         "field.ghost_width", "field.id", "field.name", "field.name.length",
         "field.output", "field.owner", "field.owner.length", "field.restart",
         "field.scalar_type", "field.space", "field.unit", "field.unit.length",
         "layout.cell_extent", "layout.face_count", "role.count"});
    static_cast<void>(collect_valid_empty(empty_roles_source));
    HUNDUN_CHECK(empty_field_record.state_fingerprint.hex !=
                 field_record.state_fingerprint.hex);

    const hundun::diagnostics::SharedFluxDiagnosticSource shared_flux{flux, 36U,
                                                                      true};
    check_fingerprint_ids(shared_flux,
                          {"face_count", "field_id", "final_flux"});
    const auto shared_record = collect_one(shared_flux);
    HUNDUN_CHECK(shared_record.module_id ==
                 "hundun.finite_volume.shared_flux");

    hundun::flow::StepAttemptReport piso_report;
    piso_report.disposition =
        hundun::flow::StepAttemptDisposition::committed;
    piso_report.reason = hundun::flow::StepFailureReason::none;
    piso_report.pressure_corrector_count = 2U;
    piso_report.attempted_dt_s = 0.125;
    piso_report.suggested_dt_s = 0.125;
    for (auto& solve : piso_report.momentum.components)
      solve.reason =
          hundun::linear::SolveTerminationReason::converged;
    for (auto& solve : piso_report.pressure)
      solve.reason =
          hundun::linear::SolveTerminationReason::converged;
    piso_report.final_transport_normalized_l2 = {0.0};
    piso_report.final_transport_relative_conservation_defect = {0.0};
    const hundun::diagnostics::ConstantDensityPisoDiagnosticSource piso_source{
        &piso_report};
    check_fingerprint_ids(piso_source, {"attempted_dt",
                                        "backflow.present",
                                        "continuity",
                                        "disposition",
                                        "lowest_failing_rank",
                                        "mass_conservation",
                                        "momentum_conservation",
                                        "momentum_residual",
                                        "pressure",
                                        "pressure_corrector_count",
                                        "reason",
                                        "solve.final_residual",
                                        "solve.global_reduction_count",
                                        "solve.initial_residual",
                                        "solve.iterations",
                                        "solve.lowest_failing_rank",
                                        "solve.matvec_count",
                                        "solve.preconditioner_apply_count",
                                        "solve.prefix",
                                        "solve.prefix.length",
                                        "solve.reason",
                                        "solve.recursive_residual",
                                        "suggested_dt",
                                        "transport.count",
                                        "transport_conservation",
                                        "transport_residual"});
    const auto piso_summary = collect_one(piso_source);
    const auto piso_counters =
        collect_one(piso_source,
                    hundun::diagnostics::DiagnosticLevel::counters);
    HUNDUN_CHECK(piso_summary.state_fingerprint.hex ==
                 piso_counters.state_fingerprint.hex);
    auto mutated_piso = piso_report;
    ++mutated_piso.pressure[1].iterations;
    const hundun::diagnostics::ConstantDensityPisoDiagnosticSource
        mutated_piso_source{&mutated_piso};
    HUNDUN_CHECK(collect_one(mutated_piso_source).state_fingerprint.hex !=
                 piso_summary.state_fingerprint.hex);

    hundun::flow::StepAttemptReport early_failed_piso;
    early_failed_piso.disposition =
        hundun::flow::StepAttemptDisposition::non_retryable_failure;
    early_failed_piso.reason =
        hundun::flow::StepFailureReason::invalid_input;
    early_failed_piso.attempted_dt_s = 0.125;
    early_failed_piso.final_continuity_normalized_l2 =
        std::numeric_limits<double>::infinity();
    early_failed_piso.final_pressure_residual_l2 =
        std::numeric_limits<double>::infinity();
    early_failed_piso.final_momentum_normalized_l2.fill(
        std::numeric_limits<double>::infinity());
    early_failed_piso.final_mass_relative_conservation_defect =
        std::numeric_limits<double>::infinity();
    early_failed_piso.final_momentum_relative_conservation_defect.fill(
        std::numeric_limits<double>::infinity());
    const hundun::diagnostics::ConstantDensityPisoDiagnosticSource
        early_failed_source{&early_failed_piso};
    check_fingerprint_ids(early_failed_source,
                          {"attempted_dt",
                           "backflow.present",
                           "continuity",
                           "disposition",
                           "lowest_failing_rank",
                           "mass_conservation",
                           "momentum_conservation",
                           "momentum_residual",
                           "pressure",
                           "pressure_corrector_count",
                           "reason",
                           "solve.final_residual",
                           "solve.global_reduction_count",
                           "solve.initial_residual",
                           "solve.iterations",
                           "solve.lowest_failing_rank",
                           "solve.matvec_count",
                           "solve.preconditioner_apply_count",
                           "solve.prefix",
                           "solve.prefix.length",
                           "solve.reason",
                           "solve.recursive_residual",
                           "suggested_dt",
                           "transport.count"});
    const auto early_failed_record = collect_valid_empty(early_failed_source);
    HUNDUN_CHECK(early_failed_record.status ==
                 hundun::diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(early_failed_record.failure.code == "flow.invalid-input");
    HUNDUN_CHECK(early_failed_record.state_fingerprint.hex !=
                 piso_summary.state_fingerprint.hex);
    auto zero_transport_committed = piso_report;
    zero_transport_committed.final_transport_normalized_l2.clear();
    zero_transport_committed.final_transport_relative_conservation_defect
        .clear();
    const hundun::diagnostics::ConstantDensityPisoDiagnosticSource
        zero_transport_source{&zero_transport_committed};
    check_fingerprint_ids(zero_transport_source,
                          {"attempted_dt",
                           "backflow.present",
                           "continuity",
                           "disposition",
                           "lowest_failing_rank",
                           "mass_conservation",
                           "momentum_conservation",
                           "momentum_residual",
                           "pressure",
                           "pressure_corrector_count",
                           "reason",
                           "solve.final_residual",
                           "solve.global_reduction_count",
                           "solve.initial_residual",
                           "solve.iterations",
                           "solve.lowest_failing_rank",
                           "solve.matvec_count",
                           "solve.preconditioner_apply_count",
                           "solve.prefix",
                           "solve.prefix.length",
                           "solve.reason",
                           "solve.recursive_residual",
                           "suggested_dt",
                           "transport.count"});
    static_cast<void>(collect_valid_empty(zero_transport_source));
    auto malformed_piso = piso_report;
    malformed_piso.final_transport_relative_conservation_defect.clear();
    rejected = false;
    try {
      static_cast<void>(collect_one(
          hundun::diagnostics::ConstantDensityPisoDiagnosticSource{
              &malformed_piso}));
    } catch (const hundun::diagnostics::DiagnosticCollectionError& error) {
      rejected =
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_input;
    }
    HUNDUN_CHECK(rejected);

    hundun::runtime::FieldRegistry unfrozen;
    const hundun::diagnostics::FieldLayoutDiagnosticSource invalid{
        &unfrozen, {{1, 1, 1}, 1U}, {}};
    rejected = false;
    try {
      static_cast<void>(collect_one(invalid));
    } catch (const hundun::diagnostics::DiagnosticCollectionError& error) {
      rejected =
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::invalid_input;
    }
    HUNDUN_CHECK(rejected);

    ThrowingSink throwing;
    rejected = false;
    try {
      hundun::diagnostics::collect_diagnostics(
          driver_source, request(), throwing);
    } catch (const hundun::diagnostics::DiagnosticCollectionError& error) {
      rejected =
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::sink_failure;
    }
    HUNDUN_CHECK(rejected);
  });
}
