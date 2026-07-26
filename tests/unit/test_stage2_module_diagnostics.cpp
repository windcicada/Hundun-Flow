// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/stage2_module_diagnostics.hpp"
#include "tests/support/test_main.hpp"

#include "hundun/linear/preconditioners.hpp"

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
      : context_(&context), layout_(2U, {0U, 1U}) {}

  hundun::linear::VectorLayout domain_layout() const override {
    return layout_;
  }
  hundun::linear::VectorLayout range_layout() const override {
    return layout_;
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
  hundun::linear::VectorLayout layout_;
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

}  // namespace

int main() {
  return hundun::test::run([] {
    hundun::execution::CpuReferenceContext execution;
    const auto context_record = collect_one(execution);
    HUNDUN_CHECK(context_record.module_id == "hundun.execution.context");

    hundun::execution::Buffer buffer(execution, 4U * sizeof(double));
    const auto allocation = buffer.allocation_identity();
    const auto epoch = buffer.epoch();
    const auto buffer_record = collect_one(buffer);
    HUNDUN_CHECK(buffer_record.module_id == "hundun.execution.buffer");
    HUNDUN_CHECK(buffer.allocation_identity() == allocation);
    HUNDUN_CHECK(buffer.epoch() == epoch);

    const auto& const_buffer = buffer;
    const auto view = const_buffer.view(1U, 2U);
    const auto view_record = collect_one(view);
    HUNDUN_CHECK(view_record.module_id == "hundun.execution.vector_view");

    QueryOnlyOperator linear_operator(execution);
    const auto operator_record =
        collect_one(static_cast<const hundun::linear::LinearOperator&>(
            linear_operator));
    HUNDUN_CHECK(operator_record.module_id == "hundun.linear.operator");
    HUNDUN_CHECK(linear_operator.apply_calls == 0U);
    HUNDUN_CHECK(linear_operator.diagonal_calls == 0U);

    hundun::linear::SolveReport solve_report;
    solve_report.reason =
        hundun::linear::SolveTerminationReason::converged;
    solve_report.iterations = 4U;
    solve_report.matvec_count = 5U;
    solve_report.preconditioner_apply_count = 4U;
    solve_report.global_reduction_count = 6U;
    const hundun::diagnostics::LinearSolveDiagnosticSource solve_source{
        "momentum-x", &solve_report};
    HUNDUN_CHECK(collect_one(solve_source).metrics.size() == 2U);
    HUNDUN_CHECK(
        collect_one(solve_source,
                    hundun::diagnostics::DiagnosticLevel::counters)
            .counters.size() == 3U);

    const hundun::diagnostics::FlowDriverDiagnosticSource driver_source{
        hundun::config::DensityModel::material,
        3U,
        0.25,
        2U,
        hundun::flow::TimeAdvanceDisposition::committed,
        hundun::flow::StepFailureReason::none,
        -1};
    HUNDUN_CHECK(collect_one(driver_source).status ==
                 hundun::diagnostics::DiagnosticStatus::ok);

    hundun::runtime::FieldRegistry unfrozen;
    const hundun::diagnostics::FieldLayoutDiagnosticSource invalid{
        &unfrozen, {{1, 1, 1}, 1U}};
    bool rejected = false;
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
