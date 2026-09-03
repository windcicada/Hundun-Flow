// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_halo_performance_counters.hpp"

#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>

namespace hundun::runtime::detail {

struct HaloTestOptions {
  std::size_t chunk_limit{static_cast<std::size_t>(INT_MAX)};
  std::size_t waitall_limit{static_cast<std::size_t>(INT_MAX)};
  int inject_post_error_rank{-1};
  int inject_wait_error_rank{-1};
  bool observe{};
  int inject_cleanup_wait_error_rank{-1};
  int inject_plan_width_first_collective_error_rank{-1};
  int inject_wire_first_collective_error_rank{-1};
  bool use_initial_performance_counters{};
  HaloPerformanceCounters initial_performance_counters{};
  int inject_unpack_failure_rank{-1};
};

struct HaloTestSnapshot {
  bool communicator_is_distinct_congruent{};
  bool communicator_uses_errors_return{};
  bool all_receives_preceded_sends{true};
  bool chunk_offsets_ordered{true};
  std::size_t receive_posts{};
  std::size_t send_posts{};
  std::size_t send_buffer_capacity{};
  std::size_t receive_buffer_capacity{};
  std::size_t request_capacity{};
  std::size_t chunk_metadata_capacity{};
  std::size_t wait_batch_capacity{};
  std::size_t destructor_drains{};
  std::size_t cleanup_wait_errors_injected{};
  std::size_t post_errors_injected{};
  std::size_t wait_errors_injected{};
  std::size_t cancel_calls{};
  std::size_t cancellation_status_checks{};
  std::size_t cancelled_requests{};
  std::size_t completed_requests{};
  std::size_t context_generation{};
  std::size_t context_replacements{};
  bool last_context_replacement_distinct_congruent{};
  std::size_t plan_width_first_collective_errors_injected{};
  std::size_t wire_first_collective_errors_injected{};
  std::size_t plan_width_second_collective_entries{};
  std::size_t wire_second_collective_entries{};
  std::size_t pack_row_copy_events{};
  std::size_t unpack_row_copy_events{};
  std::size_t first_send_sequence{};
};

void set_halo_options_raw(
    std::size_t chunk_limit, std::size_t waitall_limit,
    int inject_post_error_rank, int inject_wait_error_rank, bool observe,
    int inject_cleanup_wait_error_rank,
    int inject_plan_width_first_collective_error_rank,
    int inject_wire_first_collective_error_rank,
    bool use_initial_performance_counters,
    HaloPerformanceCounters initial_performance_counters,
    int inject_unpack_failure_rank);
void reset_halo_observation_raw() noexcept;
std::array<std::uint64_t, 29> halo_snapshot_raw() noexcept;

inline void set_halo_test_options(HaloTestOptions options) {
  set_halo_options_raw(
      options.chunk_limit, options.waitall_limit, options.inject_post_error_rank,
      options.inject_wait_error_rank, options.observe,
      options.inject_cleanup_wait_error_rank,
      options.inject_plan_width_first_collective_error_rank,
      options.inject_wire_first_collective_error_rank,
      options.use_initial_performance_counters,
      options.initial_performance_counters,
      options.inject_unpack_failure_rank);
}

inline void reset_halo_test_observation() noexcept {
  reset_halo_observation_raw();
}

inline HaloTestSnapshot halo_test_snapshot() noexcept {
  const auto values = halo_snapshot_raw();
  return HaloTestSnapshot{
      values[0] != 0U,  values[1] != 0U,  values[2] != 0U,
      values[3] != 0U,  static_cast<std::size_t>(values[4]),
      static_cast<std::size_t>(values[5]),
      static_cast<std::size_t>(values[6]),
      static_cast<std::size_t>(values[7]),
      static_cast<std::size_t>(values[8]),
      static_cast<std::size_t>(values[9]),
      static_cast<std::size_t>(values[10]),
      static_cast<std::size_t>(values[11]),
      static_cast<std::size_t>(values[12]),
      static_cast<std::size_t>(values[13]),
      static_cast<std::size_t>(values[14]),
      static_cast<std::size_t>(values[15]),
      static_cast<std::size_t>(values[16]),
      static_cast<std::size_t>(values[17]),
      static_cast<std::size_t>(values[18]),
      static_cast<std::size_t>(values[19]),
      static_cast<std::size_t>(values[20]), values[21] != 0U,
      static_cast<std::size_t>(values[22]),
      static_cast<std::size_t>(values[23]),
      static_cast<std::size_t>(values[24]),
      static_cast<std::size_t>(values[25]),
      static_cast<std::size_t>(values[26]),
      static_cast<std::size_t>(values[27]),
      static_cast<std::size_t>(values[28])};
}

} // namespace hundun::runtime::detail
