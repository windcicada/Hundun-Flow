// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "time-control diagnostics test access is tests-only"
#endif

#include <cstddef>
#include <cstdint>

namespace hundun::diagnostics::test {

enum class TimeControlDiagnosticFault : std::uint8_t {
  none,
  phase1_layout,
  phase2_request,
  phase3_provider,
  phase4_payload,
  phase4_wire,
  phase5_record,
  phase6_sink
};

class TimeControlDiagnosticsTestAccess final {
public:
  static void reset() noexcept;
  static void set_fault(TimeControlDiagnosticFault, int rank) noexcept;
  static void set_raw_fault(std::size_t ordinal) noexcept;
  static void set_request_projection_mutation(std::size_t offset,
                                              int rank) noexcept;
  static void set_provider_projection_mutation(std::size_t offset,
                                               int rank) noexcept;
  static std::size_t phase_count() noexcept;
  static std::size_t raw_operation_count() noexcept;
  static std::size_t submission_count() noexcept;
};

} // namespace hundun::diagnostics::test
