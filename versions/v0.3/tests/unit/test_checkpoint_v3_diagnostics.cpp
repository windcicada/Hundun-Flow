// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_checkpoint_v3.hpp"
#include "tests/support/test_main.hpp"

#include <vector>

namespace {

class Sink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &record) override {
    records.push_back(record);
  }
  std::vector<hundun::diagnostics::DiagnosticRecord> records;
};

} // namespace

int main() {
  return hundun::test::run([] {
    const hundun::flow::CheckpointV3Report report;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary,
        hundun::diagnostics::DiagnosticScope::local,
        {0, 0U, 0.0, "preflight"}, {}, 0U};
    Sink sink;
    hundun::diagnostics::collect_diagnostics(report, request, sink);
    HUNDUN_CHECK(sink.records.size() == 1U);
    HUNDUN_CHECK(sink.records[0].status ==
                 hundun::diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(sink.records[0].failure.classification ==
                 hundun::diagnostics::DiagnosticFailureClass::invalid_input);
    HUNDUN_CHECK(sink.records[0].metrics.size() == 4U);
    HUNDUN_CHECK(!sink.records[0].state_fingerprint.hex.empty());
  });
}
