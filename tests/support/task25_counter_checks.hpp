// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/execution/execution.hpp"
#include "hundun/runtime/halo_performance_counters.hpp"

#include <cstring>
#include <cstdint>

namespace hundun::test {

inline std::uint64_t task25_fp64_bits(double value) noexcept {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline bool task25_counters_equal(
    const runtime::HaloPerformanceCounters& left,
    const runtime::HaloPerformanceCounters& right) noexcept {
  return left.completed_exchanges == right.completed_exchanges &&
         left.begin_calls == right.begin_calls &&
         left.wait_calls == right.wait_calls &&
         left.send_payload_bytes == right.send_payload_bytes &&
         left.receive_payload_bytes == right.receive_payload_bytes &&
         left.pack_bytes == right.pack_bytes &&
         left.unpack_bytes == right.unpack_bytes &&
         left.send_messages == right.send_messages &&
         left.receive_messages == right.receive_messages &&
         task25_fp64_bits(left.completed_wait_seconds) ==
             task25_fp64_bits(right.completed_wait_seconds);
}

inline bool task25_counters_equal(
    const execution::AllocationCounters& left,
    const execution::AllocationCounters& right) noexcept {
  return left.allocation_events == right.allocation_events &&
         left.allocated_bytes == right.allocated_bytes &&
         left.deallocation_events == right.deallocation_events &&
         left.deallocated_bytes == right.deallocated_bytes &&
         left.live_bytes == right.live_bytes &&
         left.peak_live_bytes == right.peak_live_bytes;
}

}  // namespace hundun::test
