// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_reacting.hpp"

#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <vector>

namespace {

class Sink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &value) override {
    records.push_back(value);
  }
  std::vector<hundun::diagnostics::DiagnosticRecord> records;
};

hundun::diagnostics::ReactingDiagnosticsSnapshot accepted() {
  return {1U,
          77U,
          std::string(64U, 'a'),
          1.0e-12,
          2.0e-13,
          {0.0, 0.0},
          {0.0, 0.0, 0.0},
          0.0,
          2U,
          0U,
          2U,
          4U,
          1U,
          true,
          hundun::diagnostics::DiagnosticFailureClass::none,
          "none"};
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    auto snapshot = accepted();
    hundun::diagnostics::DiagnosticProviderRegistry registry;
    hundun::diagnostics::DiagnosticRequest request;
    request.scope = hundun::diagnostics::DiagnosticScope::local;
    request.level = hundun::diagnostics::DiagnosticLevel::counters;
    request.frame = {mpi.rank(), 8U, 0.08, "accepted"};
    Sink sink;
    HUNDUN_CHECK(!registry.collect_local(
        hundun::diagnostics::DiagnosticModuleKind::reacting, request, sink));
    registry.register_provider(
        hundun::diagnostics::make_reacting_diagnostic_provider(snapshot));
    const auto counters_before = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(registry.collect_local(
        hundun::diagnostics::DiagnosticModuleKind::reacting, request, sink));
    const auto counters_after = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(counters_before.collective_calls ==
                 counters_after.collective_calls);
    HUNDUN_CHECK(sink.records.size() == 1U);
    const auto &record = sink.records.front();
    HUNDUN_CHECK(record.module_kind ==
                 hundun::diagnostics::DiagnosticModuleKind::reacting);
    HUNDUN_CHECK(record.counters.size() == 5U);
    HUNDUN_CHECK(record.counters[0].id == "chemistry.calls");
    HUNDUN_CHECK(record.counters[0].value == 2U);
    HUNDUN_CHECK(record.counters[2].id == "piso.calls");
    HUNDUN_CHECK(record.counters[2].value == 2U);
    HUNDUN_CHECK(record.identities[0].subject_id == "composition");
    HUNDUN_CHECK(record.identities[1].subject_id == "mechanism");

    request.level = hundun::diagnostics::DiagnosticLevel::summary;
    Sink summary;
    HUNDUN_CHECK(registry.collect_local(
        hundun::diagnostics::DiagnosticModuleKind::reacting, request,
        summary));
    HUNDUN_CHECK(summary.records[0].metrics[0].id == "element-budget.0");
    HUNDUN_CHECK(summary.records[0].metrics[4].id == "final-residual");
    HUNDUN_CHECK(summary.records[0].metrics[5].id == "species-budget.0");

    auto failed = accepted();
    failed.revision = 2U;
    failed.accepted_step = false;
    failed.chemistry_calls = 1U;
    failed.pressure_corrector_calls = 0U;
    failed.chemistry_failures = 1U;
    failed.failure_class =
        hundun::diagnostics::DiagnosticFailureClass::non_convergence;
    failed.failure_code = "chemistry.integration";
    hundun::diagnostics::DiagnosticProviderRegistry failed_registry;
    failed_registry.register_provider(
        hundun::diagnostics::make_reacting_diagnostic_provider(failed));
    Sink failed_sink;
    HUNDUN_CHECK(failed_registry.collect_local(
        hundun::diagnostics::DiagnosticModuleKind::reacting, request,
        failed_sink));
    HUNDUN_CHECK(failed_sink.records[0].failure.classification ==
                 hundun::diagnostics::DiagnosticFailureClass::non_convergence);
    HUNDUN_CHECK(failed_sink.records[0].failure.code ==
                 "chemistry.integration");
  });
}
