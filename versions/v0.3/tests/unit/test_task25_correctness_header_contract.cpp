// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_performance_correctness.hpp"

#include <string>
#include <string_view>
#include <type_traits>

static_assert(std::is_same_v<
              decltype(&hundun::diagnostics::parse_performance_correctness),
              hundun::diagnostics::PerformanceCorrectnessRecord (*)(
                  std::string_view)>);
static_assert(std::is_same_v<
              decltype(&hundun::diagnostics::serialize_performance_correctness),
              std::string (*)(
                  const hundun::diagnostics::PerformanceCorrectnessRecord&)>);

int task25_correctness_header_contract() {
  hundun::diagnostics::PerformanceWorkRecord work;
  return work.phase == 'W' ? 0 : 1;
}
