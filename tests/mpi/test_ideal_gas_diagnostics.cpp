// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
#include "diagnostics/src/ideal_gas_closure_diagnostics_test_access.hpp"

namespace {
void run_task21_diagnostic_matrix(
    const hundun::flow::IdealGasClosureDiagnosticSource &,
    const hundun::runtime::MpiContext &, const hundun::flow::FlowState &,
    const hundun::flow::IdealGasStepAttemptReport &);
}

#define HUNDUN_TASK21_IDEAL_GAS_FIXTURE_ONLY 1
#define HUNDUN_TASK21_DIAGNOSTIC_MATRIX(source, mpi, state, report)           \
  run_task21_diagnostic_matrix(source, mpi, state, report)
#include "tests/mpi/test_ideal_gas_piso.cpp"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace {

namespace diagnostics = hundun::diagnostics;
using Access = diagnostics::test::IdealGasClosureDiagnosticTestAccess;
using Fault = diagnostics::test::IdealGasClosureDiagnosticFault;

class RecordingSink final : public diagnostics::DiagnosticSink {
public:
  void submit(const diagnostics::DiagnosticRecord &record) override {
    ++calls;
    if (throw_on_submit)
      throw std::runtime_error("injected ideal-gas diagnostic sink failure");
    records.push_back(record);
  }

  bool throw_on_submit{};
  std::size_t calls{};
  std::vector<diagnostics::DiagnosticRecord> records;
};

class FaultReset final {
public:
  ~FaultReset() noexcept { Access::reset(); }
};

struct ExactState final {
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues committed;
  hundun::flow::FlowLayerValues trial;
  hundun::flow::AcceptedStepMetadata metadata;
  hundun::flow::IdealGasClosureState closure;
};

ExactState capture(
    const hundun::flow::FlowState &state,
    const hundun::flow::IdealGasClosureDiagnosticSource &source) {
  return {state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          source.closure_state()};
}

void require_pure(
    const ExactState &expected, const hundun::flow::FlowState &state,
    const hundun::flow::IdealGasClosureDiagnosticSource &source) {
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.history, state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.committed, state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.trial, state.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      expected.metadata, state.metadata()));
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      expected.closure, source.closure_state()));
}

diagnostics::DiagnosticRequest request_for(
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    diagnostics::DiagnosticLevel level, diagnostics::DiagnosticScope scope,
    std::size_t budget = 7U,
    std::vector<std::string_view> selected_fields = {}) {
  return {level,
          scope,
          {source.relative_rank(), source.committed_step(),
           source.committed_time_s(), "ideal-gas-closure.attempt-result"},
          std::move(selected_fields),
          level == diagnostics::DiagnosticLevel::bounded_state_sample ? budget
                                                                      : 0U};
}

void collect(const hundun::flow::IdealGasClosureDiagnosticSource &source,
             const hundun::runtime::MpiContext &mpi,
             const diagnostics::DiagnosticRequest &request,
             diagnostics::DiagnosticSink &sink) {
  if (request.scope == diagnostics::DiagnosticScope::local)
    diagnostics::collect_diagnostics(source, request, sink);
  else
    diagnostics::collect_diagnostics(source, mpi, request, sink);
}

template <class Items>
void require_ids(const Items &items,
                 const std::vector<std::string_view> &expected) {
  HUNDUN_CHECK(items.size() == expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index)
    HUNDUN_CHECK(items[index].id == expected[index]);
}

void require_payload(const diagnostics::DiagnosticRecord &record) {
  static const std::vector<std::string_view> metric_ids{
      "closure.actual-mass",
      "closure.candidate-pressure",
      "closure.committed-mass",
      "closure.committed-pressure",
      "closure.density-maximum",
      "closure.density-minimum",
      "closure.enthalpy-maximum",
      "closure.enthalpy-minimum",
      "closure.enthalpy-temperature-error",
      "closure.eos-error",
      "closure.rho-h-remap-conservation",
      "closure.rho-h-remap-l2",
      "closure.rho-remap-conservation",
      "closure.rho-remap-l2",
      "closure.target-mass",
      "closure.temperature-maximum",
      "closure.temperature-minimum"};
  static const std::vector<std::string_view> invariant_ids{
      "closure.committed-eos", "closure.final-evidence-available",
      "closure.pressure-positive", "density.positive", "enthalpy.positive",
      "temperature.positive"};
  static const std::vector<std::string_view> counter_ids{
      "closure.collectives", "closure.evaluations", "closure.revision"};
  static const std::array<std::string_view, 6> identity_ids{
      "closure.attempt", "closure.state", "field.p0", "field.rho",
      "field.rho_h", "layout.cells"};

  HUNDUN_CHECK(record.identities.size() == identity_ids.size());
  for (std::size_t index = 0; index < identity_ids.size(); ++index)
    HUNDUN_CHECK(record.identities[index].subject_id == identity_ids[index]);
  if (record.level == diagnostics::DiagnosticLevel::summary) {
    require_ids(record.metrics, metric_ids);
    HUNDUN_CHECK(record.invariants.empty());
    HUNDUN_CHECK(record.counters.empty());
    HUNDUN_CHECK(record.samples.empty());
  } else if (record.level == diagnostics::DiagnosticLevel::invariants) {
    require_ids(record.invariants, invariant_ids);
    HUNDUN_CHECK(std::all_of(record.invariants.begin(),
                             record.invariants.end(),
                             [](const auto &item) { return item.passed; }));
    HUNDUN_CHECK(record.metrics.empty());
    HUNDUN_CHECK(record.counters.empty());
    HUNDUN_CHECK(record.samples.empty());
  } else if (record.level == diagnostics::DiagnosticLevel::counters) {
    require_ids(record.counters, counter_ids);
    HUNDUN_CHECK(record.metrics.empty());
    HUNDUN_CHECK(record.invariants.empty());
    HUNDUN_CHECK(record.samples.empty());
  } else {
    HUNDUN_CHECK(record.metrics.empty());
    HUNDUN_CHECK(record.invariants.empty());
    HUNDUN_CHECK(record.counters.empty());
    HUNDUN_CHECK(record.samples.size() <= record.sample_budget);
    HUNDUN_CHECK(record.eligible_sample_count >= record.samples.size());
    HUNDUN_CHECK(std::is_sorted(
        record.samples.begin(), record.samples.end(),
        [](const auto &left, const auto &right) {
          return std::tie(left.field_id, left.global_id, left.component) <
                 std::tie(right.field_id, right.global_id, right.component);
        }));
  }
}

void require_work(
    diagnostics::DiagnosticLevel level, std::size_t budget,
    const hundun::flow::IdealGasClosureDiagnosticSource &source) {
  const auto work = Access::work();
  const auto cells = source.owned_cell_count();
  HUNDUN_CHECK(work.observations == 1U);
  HUNDUN_CHECK(work.fingerprint_items ==
               2U * cells + (source.relative_rank() == 0 ? 1U : 0U));
  HUNDUN_CHECK(work.summary_items ==
               (level == diagnostics::DiagnosticLevel::summary ||
                        level == diagnostics::DiagnosticLevel::invariants
                    ? cells
                    : 0U));
  HUNDUN_CHECK(work.sample_items ==
               (level == diagnostics::DiagnosticLevel::bounded_state_sample
                    ? std::min(budget, cells)
                    : 0U));
}

template <class Operation>
void require_error(Operation &&operation,
                   diagnostics::DiagnosticFailureClass classification,
                   std::string_view code, int rank, RecordingSink &sink,
                   std::size_t expected_calls = 0U,
                   std::size_t expected_records = 0U) {
  bool caught = false;
  try {
    operation();
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    caught = true;
    HUNDUN_CHECK(error.classification() == classification);
    HUNDUN_CHECK(error.code() == code);
    HUNDUN_CHECK(error.lowest_failing_rank() == rank);
  }
  HUNDUN_CHECK(caught);
  HUNDUN_CHECK(sink.calls == expected_calls);
  HUNDUN_CHECK(sink.records.size() == expected_records);
}

void run_task21_diagnostic_matrix(
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    const hundun::runtime::MpiContext &mpi,
    const hundun::flow::FlowState &state,
    const hundun::flow::IdealGasStepAttemptReport &report) {
  FaultReset reset;
  const auto before = capture(state, source);
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   report_authenticated(report));

  for (const auto scope : {diagnostics::DiagnosticScope::local,
                           diagnostics::DiagnosticScope::collective}) {
    std::optional<std::string> fingerprint;
    for (std::uint8_t encoded = 0U; encoded < 4U; ++encoded) {
      const auto level = static_cast<diagnostics::DiagnosticLevel>(encoded);
      const auto request = request_for(source, level, scope);
      Access::reset();
      RecordingSink first;
      collect(source, mpi, request, first);
      HUNDUN_CHECK(first.calls == 1U && first.records.size() == 1U);
      const auto &record = first.records.front();
      HUNDUN_CHECK(record.schema_version ==
                   diagnostics::kDiagnosticRecordSchemaV1);
      HUNDUN_CHECK(record.module_kind ==
                   diagnostics::DiagnosticModuleKind::density_closure);
      HUNDUN_CHECK(record.module_id == "flow.ideal-gas-closure");
      HUNDUN_CHECK(record.instance_id == "primary");
      HUNDUN_CHECK(record.level == level && record.scope == scope);
      HUNDUN_CHECK(record.rank == mpi.rank());
      HUNDUN_CHECK(record.step == source.committed_step());
      HUNDUN_CHECK(record.phase == "ideal-gas-closure.attempt-result");
      HUNDUN_CHECK(record.status == diagnostics::DiagnosticStatus::ok);
      HUNDUN_CHECK(record.failure.classification ==
                   diagnostics::DiagnosticFailureClass::none);
      HUNDUN_CHECK(record.failure.code == "none");
      HUNDUN_CHECK(record.failure.lowest_failing_rank == -1);
      HUNDUN_CHECK(record.state_fingerprint.algorithm ==
                   diagnostics::kStateFingerprintAlgorithmV1);
      HUNDUN_CHECK(!record.state_fingerprint.hex.empty());
      if (!fingerprint)
        fingerprint = record.state_fingerprint.hex;
      HUNDUN_CHECK(record.state_fingerprint.hex == *fingerprint);
      require_payload(record);
      diagnostics::validate(record, diagnostics::describe_diagnostics(source),
                            request);
      const auto canonical = diagnostics::to_canonical_json(record);
      HUNDUN_CHECK(!canonical.empty());
      require_work(level, 7U, source);
      require_pure(before, state, source);

      Access::reset();
      RecordingSink repeated;
      collect(source, mpi, request, repeated);
      HUNDUN_CHECK(repeated.calls == 1U && repeated.records.size() == 1U);
      HUNDUN_CHECK(diagnostics::to_canonical_json(repeated.records.front()) ==
                   canonical);
      require_work(level, 7U, source);
      require_pure(before, state, source);
    }

    for (const std::size_t budget : {1U, 7U}) {
      auto sampled = request_for(
          source, diagnostics::DiagnosticLevel::bounded_state_sample, scope,
          budget, {"rho", "temperature"});
      Access::reset();
      RecordingSink sink;
      collect(source, mpi, sampled, sink);
      HUNDUN_CHECK(sink.calls == 1U && sink.records.size() == 1U);
      HUNDUN_CHECK(sink.records.front().sample_budget == budget);
      HUNDUN_CHECK(sink.records.front().samples.size() <= budget);
      for (const auto &sample : sink.records.front().samples)
        HUNDUN_CHECK(sample.field_id == "rho" ||
                     sample.field_id == "temperature");
      require_work(diagnostics::DiagnosticLevel::bounded_state_sample, budget,
                   source);
      require_pure(before, state, source);
    }

    for (const auto &[fault, status, classification, code] :
         {std::tuple{Fault::status_warning, diagnostics::DiagnosticStatus::warning,
                     diagnostics::DiagnosticFailureClass::none,
                     std::string_view("none")},
          std::tuple{Fault::status_failed, diagnostics::DiagnosticStatus::failed,
                     diagnostics::DiagnosticFailureClass::non_convergence,
                     std::string_view("closure.injected-failure")},
          std::tuple{Fault::status_unavailable,
                     diagnostics::DiagnosticStatus::unavailable,
                     diagnostics::DiagnosticFailureClass::unavailable,
                     std::string_view("closure.injected-unavailable")}}) {
      Access::set_fault(fault);
      RecordingSink sink;
      const auto request =
          request_for(source, diagnostics::DiagnosticLevel::summary, scope);
      collect(source, mpi, request, sink);
      HUNDUN_CHECK(sink.calls == 1U && sink.records.size() == 1U);
      HUNDUN_CHECK(sink.records.front().status == status);
      HUNDUN_CHECK(sink.records.front().failure.classification ==
                   classification);
      HUNDUN_CHECK(sink.records.front().failure.code == code);
      HUNDUN_CHECK(sink.records.front().failure.lowest_failing_rank ==
                   (scope == diagnostics::DiagnosticScope::collective &&
                            status != diagnostics::DiagnosticStatus::warning
                        ? 0
                        : -1));
      diagnostics::validate(sink.records.front(),
                            diagnostics::describe_diagnostics(source), request);
      HUNDUN_CHECK(!diagnostics::to_canonical_json(sink.records.front()).empty());
      require_pure(before, state, source);
      Access::reset();
    }
  }

  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    ++request.frame.step;
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", -1, sink);
  }
  {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::local, 0U);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", -1, sink);
  }
  for (auto fields : {std::vector<std::string_view>{"unknown"},
                      std::vector<std::string_view>{"temperature", "rho"}}) {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::local, 7U, std::move(fields));
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.selected-field", -1, sink);
  }
  {
    Access::set_fault(Fault::record_validation, mpi.rank());
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "closure.diagnostics.record", -1, sink);
    Access::reset();
  }
  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    RecordingSink sink;
    sink.throw_on_submit = true;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::sink_failure,
        "diagnostics.sink.submit", -1, sink, 1U);
  }

  const int target = mpi.size() > 1 ? 1 : 0;
  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    if (mpi.rank() == target)
      ++request.frame.step;
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", target, sink);
  }
  {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 0U);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", 0, sink);
  }
  {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 7U);
    if (mpi.rank() == target)
      request.selected_fields = {"unknown"};
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.selected-field", target, sink);
  }
  if (mpi.size() > 1) {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 7U);
    if (mpi.rank() == target)
      request.selected_fields = {"rho"};
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.request-agreement", target, sink);
  }
  {
    Access::set_fault(Fault::provider_agreement, target);
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "closure.diagnostics.provider-agreement", target, sink);
    Access::reset();
  }
  {
    Access::set_fault(Fault::ownership_layout, target);
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::layout,
        "closure.diagnostics.ownership", target, sink);
    Access::reset();
  }
  {
    Access::set_fault(Fault::sample_wire, target);
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 7U);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::layout,
        "closure.diagnostics.sample-wire", target, sink);
    Access::reset();
  }
  {
    Access::set_fault(Fault::record_validation, target);
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "closure.diagnostics.record", target, sink);
    Access::reset();
  }
  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    sink.throw_on_submit = mpi.rank() == target;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::sink_failure,
        "diagnostics.sink.submit", target, sink, 1U,
        mpi.rank() == target ? 0U : 1U);
  }
  require_pure(before, state, source);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const auto descriptor =
        diagnostics::describe_ideal_gas_closure_diagnostics();
    HUNDUN_CHECK(descriptor.module_kind ==
                 diagnostics::DiagnosticModuleKind::density_closure);
    HUNDUN_CHECK(descriptor.module_id ==
                 std::string_view("flow.ideal-gas-closure"));
    HUNDUN_CHECK(descriptor.instance_id == std::string_view("primary"));
    const auto positive = Access::positive_invariant(
        "density.positive", "kg/m3", 1.0);
    HUNDUN_CHECK(positive.relation ==
                 diagnostics::InvariantRelation::positive);
    HUNDUN_CHECK(positive.limit.status ==
                 diagnostics::DiagnosticValueStatus::unavailable);
    HUNDUN_CHECK(positive.limit.bits == 0U);
    HUNDUN_CHECK(positive.passed);
    run(mpi, false, false);
  });
}
