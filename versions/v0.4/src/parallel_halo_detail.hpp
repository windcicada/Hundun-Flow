// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace hundun::v04::detail {

// Stable detail values are intentionally kept out of the public Status API.
// Focused tests use them to distinguish state-machine phases without adding a
// branch or callback to the production HaloEngine interface.
inline constexpr std::uint32_t halo_detail_state = 301U;
inline constexpr std::uint32_t halo_detail_input = 302U;
inline constexpr std::uint32_t halo_detail_topology = 303U;
inline constexpr std::uint32_t halo_detail_field = 304U;
inline constexpr std::uint32_t halo_detail_overflow = 305U;
inline constexpr std::uint32_t halo_detail_collective = 306U;
inline constexpr std::uint32_t halo_detail_request = 307U;
inline constexpr std::uint32_t halo_detail_pack_failure = 308U;
inline constexpr std::uint32_t halo_detail_start_failure = 309U;
inline constexpr std::uint32_t halo_detail_completion_failure = 310U;
inline constexpr std::uint32_t halo_detail_unpack_failure = 311U;
inline constexpr std::uint32_t halo_detail_reserve_allocation = 312U;

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)

enum class HaloFailurePoint : std::uint8_t {
  none,
  pack,
  start,
  completion,
  unpack,
  reserve_before_contract,
  reserve_before_alltoall
};

// Set on every participating process.  The selected phase fails only on
// failing_rank; HaloEngine then exercises its normal collective consensus.
void set_halo_failure_for_test(HaloFailurePoint point,
                               int failing_rank) noexcept;
void clear_halo_failure_for_test() noexcept;
void set_halo_maximum_chunk_doubles_for_test(
    std::uint64_t maximum) noexcept;
void clear_halo_maximum_chunk_doubles_for_test() noexcept;

#endif

}  // namespace hundun::v04::detail
