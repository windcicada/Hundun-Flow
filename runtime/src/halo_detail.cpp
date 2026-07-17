// SPDX-License-Identifier: Apache-2.0

#include "halo_detail.hpp"

#include "hundun/runtime/error.hpp"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace hundun::runtime::detail {
namespace {

bool legal_offset(Int3 offset) noexcept {
  const bool in_range = offset.x >= -1 && offset.x <= 1 &&
                        offset.y >= -1 && offset.y <= 1 &&
                        offset.z >= -1 && offset.z <= 1;
  const bool nonzero = offset.x != 0 || offset.y != 0 || offset.z != 0;
  return in_range && nonzero;
}

std::size_t checked_multiply(std::size_t left, std::size_t right,
                             std::string_view quantity) {
  constexpr std::size_t limit = std::numeric_limits<std::size_t>::max();
  if (right != 0U && left > limit / right) {
    throw Error("halo " + std::string(quantity) + " overflow");
  }
  return left * right;
}

std::size_t axis_span(int begin, int end) {
  if (end < begin) {
    throw Error("halo region box has a negative axis span");
  }
  const auto span = static_cast<std::int64_t>(end) -
                    static_cast<std::int64_t>(begin);
  return static_cast<std::size_t>(span);
}

}  // namespace

int halo_offset_code(Int3 offset) {
  if (!legal_offset(offset)) {
    throw Error(
        "halo offset must identify a nonzero direction in the 3x3x3 stencil");
  }
  return (offset.z + 1) * 9 + (offset.y + 1) * 3 + (offset.x + 1);
}

Int3 reverse_offset(Int3 offset) {
  static_cast<void>(halo_offset_code(offset));
  return Int3{-offset.x, -offset.y, -offset.z};
}

int halo_receive_tag(Int3 offset) {
  return halo_offset_code(reverse_offset(offset));
}

int effective_halo_tag_upper_bound(bool attribute_present,
                                   const int* upper_bound) {
  constexpr int mpi3_guaranteed_minimum = 32767;
  if (!attribute_present) {
    return mpi3_guaranteed_minimum;
  }
  if (upper_bound == nullptr) {
    throw Error("MPI_TAG_UB attribute pointer is null");
  }
  if (*upper_bound < 26) {
    throw Error("MPI tag upper bound is below the halo direction range");
  }
  return *upper_bound;
}

std::size_t checked_region_payload_bytes(Box3 box, std::uint32_t components,
                                         std::size_t scalar_bytes) {
  if (components == 0U) {
    throw Error("halo payload requires a positive component count");
  }
  if (scalar_bytes == 0U) {
    throw Error("halo payload requires a positive scalar byte size");
  }
  const std::size_t nx = axis_span(box.begin.x, box.end.x);
  const std::size_t ny = axis_span(box.begin.y, box.end.y);
  const std::size_t nz = axis_span(box.begin.z, box.end.z);
  std::size_t count = checked_multiply(nx, ny, "payload size");
  count = checked_multiply(count, nz, "payload size");
  count = checked_multiply(count, static_cast<std::size_t>(components),
                           "payload size");
  return checked_multiply(count, scalar_bytes, "payload byte size");
}

std::vector<CountRange> split_count_ranges(std::size_t total,
                                           std::size_t maximum_count,
                                           std::string_view quantity) {
  std::vector<CountRange> ranges;
  split_count_ranges_into(total, maximum_count, quantity, ranges);
  return ranges;
}

void split_count_ranges_into(std::size_t total, std::size_t maximum_count,
                             std::string_view quantity,
                             std::vector<CountRange>& ranges) {
  if (maximum_count == 0U) {
    throw Error("halo " + std::string(quantity) +
                " limit must be positive");
  }
  if (maximum_count > static_cast<std::size_t>(INT_MAX)) {
    throw Error("halo " + std::string(quantity) +
                " limit exceeds the MPI int count range");
  }
  if (total == 0U) {
    ranges.clear();
    return;
  }

  const std::size_t range_count = 1U + (total - 1U) / maximum_count;
  if (range_count > ranges.max_size()) {
    throw Error("halo " + std::string(quantity) +
                " range metadata exceeds local capacity");
  }
  if (range_count > ranges.capacity()) {
    ranges.reserve(range_count);
  }
  ranges.clear();
  std::size_t offset = 0U;
  for (std::size_t index = 0U; index < range_count; ++index) {
    const std::size_t remaining = total - offset;
    const std::size_t count =
        remaining < maximum_count ? remaining : maximum_count;
    ranges.push_back(
        CountRange{offset, static_cast<int>(count)});
    offset += count;
  }
}

ActiveDestructionAction active_destruction_action(
    bool active, bool mpi_live, bool requests_proven_complete) noexcept {
  if (!active || requests_proven_complete) {
    return ActiveDestructionAction::clear_without_mpi;
  }
  return mpi_live ? ActiveDestructionAction::drain_requests
                  : ActiveDestructionAction::terminate_process;
}

CompletionFailureAction completion_failure_action(
    bool requests_proven_complete) noexcept {
  return requests_proven_complete
             ? CompletionFailureAction::discard_and_throw
             : CompletionFailureAction::terminate_process;
}

}  // namespace hundun::runtime::detail
