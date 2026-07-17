// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/types.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace hundun::runtime::detail {

struct CountRange {
  std::size_t offset{};
  int count{};
};

int halo_offset_code(Int3 offset);
Int3 reverse_offset(Int3 offset);
int halo_receive_tag(Int3 offset);

void validate_halo_tag_upper_bound(bool attribute_present, int upper_bound);

std::size_t checked_region_payload_bytes(Box3 box, std::uint32_t components,
                                         std::size_t scalar_bytes);

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

CompletionFailureAction completion_failure_action(
    bool requests_proven_complete) noexcept;

}  // namespace hundun::runtime::detail
