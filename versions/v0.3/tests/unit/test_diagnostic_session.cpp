// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_session.hpp"
#include "hundun/diag_structured.hpp"
#include "tests/support/test_main.hpp"

#include <filesystem>
#include <string>

namespace {

hundun::diagnostics::DiagnosticRecord record() {
  using namespace hundun::diagnostics;
  DiagnosticRecord result;
  result.module_kind = DiagnosticModuleKind::flow_driver;
  result.module_id = "hundun.application.flow_driver";
  result.instance_id = "flow-driver";
  result.rank = 0;
  result.step = 2;
  result.time_s = describe_fp64(0.25);
  result.phase = "accepted-step";
  result.metrics.push_back(
      {"attempt-count", DiagnosticMetricKind::state_summary, "count",
       describe_fp64(1.0)});
  DiagnosticFingerprintAccumulator fingerprint;
  fingerprint.add("attempt-count", 0, 0, describe_fp64(1.0));
  result.state_fingerprint = fingerprint.finish();
  return result;
}

}  // namespace

int main() {
  return hundun::test::run([] {
    using namespace hundun::diagnostics;
    const auto directory = std::filesystem::temp_directory_path() /
                           "hundun-task24-session-unit";
    std::filesystem::remove_all(directory);
    DiagnosticSession session(directory, 2, 0);
    HUNDUN_CHECK(!session.due(1));
    HUNDUN_CHECK(session.due(2));
    HUNDUN_CHECK(!std::filesystem::exists(directory));

    const DiagnosticDescriptor descriptor{
        kDiagnosticRecordSchemaV1, DiagnosticModuleKind::flow_driver,
        "hundun.application.flow_driver", "flow-driver",
        static_cast<DiagnosticCapabilityFlags>(
            DiagnosticCapability::summary)};
    const DiagnosticRequest request{
        DiagnosticLevel::summary, DiagnosticScope::local,
        {0, 2, 0.25, "accepted-step"}, {}, 0};
    DiagnosticBatch batch;
    DiagnosticBatchSink sink(descriptor, request, batch);
    sink.submit(record());
    HUNDUN_CHECK(batch.size() == 1U);
    HUNDUN_CHECK(batch.canonical_json_lines().find(
                     "\"module_id\":\"hundun.application.flow_driver\"") !=
                 std::string::npos);

    auto mismatched = record();
    mismatched.step = 3;
    bool rejected = false;
    try {
      DiagnosticBatchSink bad_sink(descriptor, request, batch);
      bad_sink.submit(mismatched);
    } catch (const DiagnosticCollectionError& error) {
      rejected =
          error.classification() == DiagnosticFailureClass::invalid_input;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(batch.size() == 1U);
    std::filesystem::remove_all(directory);
  });
}
