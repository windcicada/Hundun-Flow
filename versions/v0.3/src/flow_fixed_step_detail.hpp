// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/flow_state.hpp"

namespace hundun::flow::detail {

struct FlowStateSolverAccess final {
  static const runtime::FieldRegistry &registry(const FlowState &state) {
    return state.solver_registry();
  }
  static const runtime::FieldAccessPlan &access(const FlowState &state) {
    return state.solver_access_plan();
  }
  static runtime::FieldStorage &layer(FlowState &state, FlowLayer selected) {
    return state.solver_layer(selected);
  }
  static std::uint64_t instance_identity(const FlowState &) noexcept;
};

} // namespace hundun::flow::detail
