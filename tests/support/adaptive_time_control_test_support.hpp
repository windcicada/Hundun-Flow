// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flow/src/material_density_piso_test_access.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "runtime/src/field_epoch_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"

#include <cstdint>

namespace hundun::test {

inline AdaptiveFlowStateSnapshot capture_adaptive_flow_state(
    const flow::FlowState &state, const flow::TimeControlState &controller,
    flow::IdealGasClosureState closure = {},
    std::uint64_t cache_identity = 0U, std::size_t cache_capacity = 0U,
    std::uint64_t cache_revision = 0U,
    std::array<std::uint64_t, 4> operation_counters = {}) {
  AdaptiveFlowStateSnapshot result;
  result.history = state.snapshot(flow::FlowLayer::history);
  result.committed = state.snapshot(flow::FlowLayer::committed);
  result.trial = state.snapshot(flow::FlowLayer::trial);
  result.metadata = state.metadata();
  result.controller = controller;
  result.closure = closure;
  result.history_allocation_identity =
      flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
          state, flow::FlowLayer::history);
  result.committed_allocation_identity =
      flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
          state, flow::FlowLayer::committed);
  result.trial_allocation_identity =
      flow::test::MaterialDensityPisoTestAccess::state_allocation_identity(
          state, flow::FlowLayer::trial);
  const auto &history = state.layer(flow::FlowLayer::history);
  const auto &committed = state.layer(flow::FlowLayer::committed);
  const auto &trial = state.layer(flow::FlowLayer::trial);
  result.history_generation =
      runtime::detail::FieldEpochTestAccess::generation(history);
  result.committed_generation =
      runtime::detail::FieldEpochTestAccess::generation(committed);
  result.trial_generation =
      runtime::detail::FieldEpochTestAccess::generation(trial);
  result.attempt_identity =
      flow::test::MaterialDensityPisoTestAccess::state_attempt_identity(state);
  result.diagnostic_identity =
      flow::test::MaterialDensityPisoTestAccess::state_diagnostic_identity(
          state);
  result.accepted_cache_identity = cache_identity;
  result.accepted_cache_capacity = cache_capacity;
  result.accepted_cache_revision = cache_revision;
  result.operation_counters = operation_counters;
  return result;
}

inline bool adaptive_flow_state_failed_attempt_preserved(
    const AdaptiveFlowStateSnapshot &before,
    const AdaptiveFlowStateSnapshot &after,
    std::uint64_t attempted_work_count,
    std::array<std::uint64_t, 4> expected_operation_delta) noexcept {
  auto normalized = after;
  normalized.trial_generation = before.trial_generation;
  normalized.attempt_identity = before.attempt_identity;
  normalized.diagnostic_identity = before.diagnostic_identity;
  normalized.operation_counters = before.operation_counters;
  std::array<std::uint64_t, 4> observed_operation_delta{};
  for (std::size_t index = 0; index < observed_operation_delta.size();
       ++index) {
    if (after.operation_counters[index] < before.operation_counters[index])
      return false;
    observed_operation_delta[index] =
        after.operation_counters[index] - before.operation_counters[index];
  }
  return adaptive_flow_state_bitwise_equal(before, normalized) &&
         after.trial_generation ==
             before.trial_generation + attempted_work_count &&
         after.attempt_identity ==
             before.attempt_identity + attempted_work_count &&
         after.diagnostic_identity ==
             before.diagnostic_identity + attempted_work_count &&
         observed_operation_delta == expected_operation_delta;
}

} // namespace hundun::test
