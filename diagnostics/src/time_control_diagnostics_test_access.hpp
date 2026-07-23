// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "time-control diagnostics test access is tests-only"
#endif

#include <cstddef>
#include <cstdint>

namespace hundun::flow {
class TimeControlDiagnosticSource;
enum class StepFailureReason : std::uint8_t;
enum class TimeAdvanceDisposition : std::uint8_t;
}

namespace hundun::diagnostics::test {

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
  static void reset() noexcept;
  static void set_fault(TimeControlDiagnosticFault, int rank) noexcept;
  static void set_raw_fault(std::size_t ordinal, int rank) noexcept;
  static void set_request_projection_mutation(std::size_t offset,
                                              int rank) noexcept;
  static void set_provider_projection_mutation(std::size_t offset,
                                               int rank) noexcept;
  static void set_provider_mutation(TimeControlProviderMutation,
                                    int rank) noexcept;
  static void set_wire_mutation(TimeControlWireMutation, int rank) noexcept;
  static void set_local_mutation(flow::TimeControlDiagnosticSource &,
                                 TimeControlLocalMutation);
  static std::size_t phase_count() noexcept;
  static std::size_t raw_operation_count() noexcept;
  static std::size_t submission_count() noexcept;
  static void set_state_counters(flow::TimeControlDiagnosticSource &,
                                 std::uint64_t accepted_step,
                                 std::uint64_t revision) noexcept;
  static void set_outcome(flow::TimeControlDiagnosticSource &,
                          flow::TimeAdvanceDisposition,
                          flow::StepFailureReason, int lowest_failing_rank,
                          std::uint8_t preflight_category,
                          std::size_t attempt_count) noexcept;
};

} // namespace hundun::diagnostics::test
