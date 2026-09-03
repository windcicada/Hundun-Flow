// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_stage3_evidence_manifest.hpp"
#include "tests/support/test_main.hpp"

#include <stdexcept>

namespace {

hundun::diagnostics::Stage3EvidenceRecord record() {
  return {"row", "0123456789012345678901234567890123456789", "tree", "clean",
          "debug", "/build", "cache", "binary", 7U, "clang", "libc++",
          "openmpi", "test --arg", "OMP_NUM_THREADS=1", "0-3", 1,
          "start", "end", 0, 1.25, 1024U, "log", {"artifact"}};
}

template <class Function> void rejected(Function function) {
  bool threw = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  HUNDUN_CHECK(threw);
}

} // namespace

int main() {
  return hundun::test::run([] {
    auto value = record();
    const auto json =
        hundun::diagnostics::stage3_evidence_manifest_json({value});
    HUNDUN_CHECK(json.find("\"exit_status\":0") != std::string::npos);
    auto missing = value;
    missing.log_sha256.clear();
    rejected([&] { hundun::diagnostics::validate_stage3_evidence_manifest({missing}); });
    missing = value;
    missing.binary_sha256.clear();
    rejected([&] { hundun::diagnostics::validate_stage3_evidence_manifest({missing}); });
    missing = value;
    missing.exit_status = -1;
    rejected([&] { hundun::diagnostics::validate_stage3_evidence_manifest({missing}); });
    rejected([&] { hundun::diagnostics::validate_stage3_evidence_manifest({value, value}); });
  });
}
