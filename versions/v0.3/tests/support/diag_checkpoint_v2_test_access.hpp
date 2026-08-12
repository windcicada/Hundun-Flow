// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>

namespace hundun::diagnostics::detail {

void checkpoint_v2_diagnostics_set_fault_for_test(std::uint8_t,
                                                  int) noexcept;
void checkpoint_v2_diagnostics_reset_for_test() noexcept;
std::uint64_t checkpoint_v2_diagnostics_collective_calls_for_test() noexcept;

} // namespace hundun::diagnostics::detail

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
  static void set_fault(CheckpointV2DiagnosticFault fault,
                        int rank = -1) noexcept {
    detail::checkpoint_v2_diagnostics_set_fault_for_test(
        static_cast<std::uint8_t>(fault), rank);
  }
  static void reset() noexcept {
    detail::checkpoint_v2_diagnostics_reset_for_test();
  }
  static CheckpointV2DiagnosticWork work() noexcept {
    return {detail::checkpoint_v2_diagnostics_collective_calls_for_test()};
  }
};

} // namespace hundun::diagnostics::test
