// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <climits>
#include <cstddef>

namespace hundun::runtime::detail {

struct HaloTestOptions {
  std::size_t chunk_limit{static_cast<std::size_t>(INT_MAX)};
  std::size_t waitall_limit{static_cast<std::size_t>(INT_MAX)};
  int inject_post_error_rank{-1};
  int inject_wait_error_rank{-1};
  bool observe{};
  int inject_cleanup_wait_error_rank{-1};
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
};

void set_halo_test_options(HaloTestOptions options);
void reset_halo_test_observation() noexcept;
HaloTestSnapshot halo_test_snapshot() noexcept;

}  // namespace hundun::runtime::detail
