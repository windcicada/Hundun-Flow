// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/flow_state.hpp"

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
};

} // namespace hundun::flow::detail
