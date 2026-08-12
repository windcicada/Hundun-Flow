// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_adaptive_time_control.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace hundun::diagnostics::detail {

void time_control_diagnostic_reset_raw() noexcept;
void time_control_diagnostic_set_fault_raw(std::uint8_t, int) noexcept;
void time_control_diagnostic_set_raw_fault_raw(std::size_t, int) noexcept;
void time_control_diagnostic_set_request_mutation_raw(std::size_t,
                                                      int) noexcept;
void time_control_diagnostic_set_provider_mutation_offset_raw(std::size_t,
                                                              int) noexcept;
void time_control_diagnostic_set_provider_mutation_raw(std::uint8_t,
                                                       int) noexcept;
void time_control_diagnostic_set_wire_mutation_raw(std::uint8_t,
                                                   int) noexcept;
std::size_t time_control_diagnostic_phase_count_raw() noexcept;
std::size_t time_control_diagnostic_raw_count_raw() noexcept;
std::array<std::size_t, 3>
time_control_diagnostic_raw_fault_observation_raw() noexcept;
int time_control_diagnostic_raw_fault_rank_raw() noexcept;
std::size_t time_control_diagnostic_submission_count_raw() noexcept;
void time_control_diagnostic_set_state_counters_raw(
    flow::TimeControlDiagnosticSource &, std::uint64_t,
    std::uint64_t) noexcept;
void time_control_diagnostic_set_outcome_raw(
    flow::TimeControlDiagnosticSource &, flow::TimeAdvanceDisposition,
    flow::StepFailureReason, int, std::uint8_t, std::size_t) noexcept;
void time_control_diagnostic_set_local_mutation_raw(
    flow::TimeControlDiagnosticSource &, std::uint8_t);

} // namespace hundun::diagnostics::detail

namespace hundun::diagnostics::test {

struct TimeControlRawFaultObservation final {
  std::size_t raw_operations{};
  std::size_t local_origins{};
  std::size_t fault_ordinal{};
  int requested_rank{-1};
};

enum class TimeControlDiagnosticFault : std::uint8_t {
  none,
  phase1_layout,
  phase2_request,
  projection_root_size,
  projection_payload,
  phase3_provider,
  phase4_payload,
  fingerprint_aggregation,
  eligible_aggregation,
  wire_root_size,
  wire_payload,
  phase4_wire,
  phase5_record,
  phase6_sink
};

enum class TimeControlWireMutation : std::uint8_t {
  none,
  short_payload,
  trailing_payload,
  count,
  field,
  global_id,
  component,
  unit,
  value,
  duplicate,
  missing_limb,
  swapped_limbs,
  tuple_order
};

enum class TimeControlProviderMutation : std::uint8_t {
  none,
  state_seal,
  summary,
  config,
  model,
  frame,
  identities,
  global_cell_authority,
  global_cell_layout,
  global_face_authority,
  global_face_layout
};

enum class TimeControlLocalMutation : std::uint8_t {
  local_box,
  canonical_owned_faces,
  local_faces,
  local_cell_layout,
  global_cell_layout,
  local_face_layout,
  global_face_layout
};

class TimeControlDiagnosticsTestAccess final {
public:
  static void reset() noexcept { detail::time_control_diagnostic_reset_raw(); }
  static void set_fault(TimeControlDiagnosticFault fault, int rank) noexcept {
    detail::time_control_diagnostic_set_fault_raw(
        static_cast<std::uint8_t>(fault), rank);
  }
  static void set_raw_fault(std::size_t ordinal, int rank) noexcept {
    detail::time_control_diagnostic_set_raw_fault_raw(ordinal, rank);
  }
  static void set_request_projection_mutation(std::size_t offset,
                                              int rank) noexcept {
    detail::time_control_diagnostic_set_request_mutation_raw(offset, rank);
  }
  static void set_provider_projection_mutation(std::size_t offset,
                                               int rank) noexcept {
    detail::time_control_diagnostic_set_provider_mutation_offset_raw(offset,
                                                                    rank);
  }
  static void set_provider_mutation(TimeControlProviderMutation,
                                    int rank) noexcept;
  static void set_wire_mutation(TimeControlWireMutation mutation,
                                int rank) noexcept {
    detail::time_control_diagnostic_set_wire_mutation_raw(
        static_cast<std::uint8_t>(mutation), rank);
  }
  static void set_local_mutation(flow::TimeControlDiagnosticSource &,
                                 TimeControlLocalMutation);
  static std::size_t phase_count() noexcept {
    return detail::time_control_diagnostic_phase_count_raw();
  }
  static std::size_t raw_operation_count() noexcept {
    return detail::time_control_diagnostic_raw_count_raw();
  }
  static TimeControlRawFaultObservation raw_fault_observation() noexcept {
    const auto raw =
        detail::time_control_diagnostic_raw_fault_observation_raw();
    return {raw[0], raw[1], raw[2],
            detail::time_control_diagnostic_raw_fault_rank_raw()};
  }
  static std::size_t submission_count() noexcept {
    return detail::time_control_diagnostic_submission_count_raw();
  }
  static void set_state_counters(flow::TimeControlDiagnosticSource &source,
                                 std::uint64_t accepted_step,
                                 std::uint64_t revision) noexcept {
    detail::time_control_diagnostic_set_state_counters_raw(
        source, accepted_step, revision);
  }
  static void set_outcome(flow::TimeControlDiagnosticSource &source,
                          flow::TimeAdvanceDisposition disposition,
                          flow::StepFailureReason reason,
                          int lowest_failing_rank,
                          std::uint8_t preflight_category,
                          std::size_t attempt_count) noexcept {
    detail::time_control_diagnostic_set_outcome_raw(
        source, disposition, reason, lowest_failing_rank, preflight_category,
        attempt_count);
  }
};

inline void TimeControlDiagnosticsTestAccess::set_provider_mutation(
    TimeControlProviderMutation mutation, int rank) noexcept {
  detail::time_control_diagnostic_set_provider_mutation_raw(
      static_cast<std::uint8_t>(mutation), rank);
}

inline void TimeControlDiagnosticsTestAccess::set_local_mutation(
    flow::TimeControlDiagnosticSource &source,
    TimeControlLocalMutation mutation) {
  detail::time_control_diagnostic_set_local_mutation_raw(
      source, static_cast<std::uint8_t>(mutation));
}

} // namespace hundun::diagnostics::test
