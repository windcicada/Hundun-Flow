// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "checkpoint-v2 diagnostic test access is unavailable in tests-off builds"
#endif

#include <cstdint>

namespace hundun::diagnostics::test {

enum class CheckpointV2DiagnosticFault : std::uint8_t {
  none,
  raw_mpi,
  request_preparation,
  local_record_preparation,
  post_gather_preparation,
  final_record_preparation
};

struct CheckpointV2DiagnosticWork final {
  std::uint64_t collective_calls{};
};

class CheckpointV2DiagnosticsTestAccess final {
public:
  static void set_fault(CheckpointV2DiagnosticFault, int rank = -1) noexcept;
  static void reset() noexcept;
  static CheckpointV2DiagnosticWork work() noexcept;
};
}
