// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace hundun::runtime::detail {

struct CountRange {
  std::size_t offset{};
  int count{};
};

struct HaloRowLayout {
  int ghost_width{};
  std::uint32_t components{};
  std::size_t scalar_bytes{};
  std::size_t x_stride{};
  std::size_t y_stride{};
  std::size_t z_stride{};
  std::size_t field_bytes{};
};

enum class PostEvent { receive, send };

struct PostEventState {
  std::size_t expected_receive_posts{};
  std::size_t receive_posts{};
  std::size_t send_posts{};
  std::size_t sequence{};
  std::size_t first_send_sequence{static_cast<std::size_t>(-1)};
  bool send_seen{};
  bool valid{true};
};

enum class HaloFailureCategory : std::int64_t {
  none = 0,
  post = 1,
  completion = 2,
  performance_counter = 3
};

enum class HaloMpiOperation : std::int64_t {
  none = 0,
  irecv = 1,
  isend = 2,
  waitall = 3,
  wait = 4
};

struct HaloFailureRecord {
  std::int64_t category{};
  std::int64_t failing_rank{-1};
  std::int64_t operation{};
  std::int64_t mpi_result{};
  std::int64_t region_index{-1};
  std::int64_t chunk_offset{-1};
  std::int64_t chunk_count{};
  std::int64_t tag{-1};
};

void observe_post_event(PostEventState& state, PostEvent event) noexcept;
bool post_event_sequence_valid(const PostEventState& state) noexcept;

int halo_offset_code(Int3 offset);
Int3 reverse_offset(Int3 offset);
int halo_receive_tag(Int3 offset);

int effective_halo_tag_upper_bound(bool attribute_present,
                                   const int* upper_bound);

std::size_t checked_region_payload_bytes(Box3 box, std::uint32_t components,
                                         std::size_t scalar_bytes);

std::size_t pack_halo_region_rows(
    const std::byte* field_data, const HaloRowLayout& layout, Box3 box,
    std::byte* wire_data, std::size_t wire_bytes);

std::size_t unpack_halo_region_rows(
    std::byte* field_data, const HaloRowLayout& layout, Box3 box,
    const std::byte* wire_data, std::size_t wire_bytes);

std::vector<CountRange> split_count_ranges(std::size_t total,
                                           std::size_t maximum_count,
                                           std::string_view quantity);

void split_count_ranges_into(std::size_t total, std::size_t maximum_count,
                             std::string_view quantity,
                             std::vector<CountRange>& ranges);

enum class ActiveDestructionAction {
  clear_without_mpi,
  drain_requests,
  terminate_process
};

ActiveDestructionAction active_destruction_action(
    bool active, bool mpi_live, bool requests_proven_complete) noexcept;

enum class CompletionFailureAction { discard_and_throw, terminate_process };

enum class FailureRecoveryAction {
  return_idle,
  replace_context,
  terminate_process
};

struct CompletionOutcome {
  bool mpi_error_seen{};
  bool requests_proven_null{};
};

bool completion_succeeded(CompletionOutcome outcome) noexcept;

CompletionFailureAction completion_failure_action(
    CompletionOutcome outcome) noexcept;

FailureRecoveryAction failure_recovery_action(
    bool cancellation_calls_succeeded, bool requests_proven_null,
    bool channel_quiescent, bool context_replaced) noexcept;

}  // namespace hundun::runtime::detail
