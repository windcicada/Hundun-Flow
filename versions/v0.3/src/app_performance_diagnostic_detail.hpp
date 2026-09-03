// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstdint>
#include <utility>

namespace hundun::application::detail {

struct FlowDiagnosticBranchResult final {
  bool submitted{};
  std::uint64_t logical_bytes{};
};

template <class DueAction>
FlowDiagnosticBranchResult execute_flow_diagnostic_branch(
    bool enabled, std::uint64_t write_interval, std::uint64_t step,
    DueAction&& due_action) {
  if (!enabled || step == 0U || write_interval == 0U ||
      step % write_interval != 0U) {
    return {};
  }
  return {true, std::forward<DueAction>(due_action)()};
}

}  // namespace hundun::application::detail
