// SPDX-License-Identifier: Apache-2.0

#include "halo_detail.hpp"

#include "hundun/runtime/error.hpp"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

std::size_t checked_add(std::size_t left, std::size_t right,
                        std::string_view quantity) {
  constexpr std::size_t limit = std::numeric_limits<std::size_t>::max();
  if (left > limit - right) {
    throw Error("halo " + std::string(quantity) + " overflow");
  }
  return left + right;
}

std::size_t padded_coordinate(int coordinate, int ghost_width) {
  if (ghost_width < 0) {
    throw Error("halo row layout ghost width must not be negative");
  }
  const auto padded = static_cast<std::int64_t>(coordinate) +
                      static_cast<std::int64_t>(ghost_width);
  if (padded < 0) {
    throw Error("halo row coordinate precedes allocated field storage");
  }
  return static_cast<std::size_t>(padded);
}

std::size_t row_field_byte_offset(const HaloRowLayout& layout, int i, int j,
                                  int k, std::size_t row_bytes) {
  const std::size_t ii = padded_coordinate(i, layout.ghost_width);
  const std::size_t jj = padded_coordinate(j, layout.ghost_width);
  const std::size_t kk = padded_coordinate(k, layout.ghost_width);
  std::size_t linear = checked_multiply(kk, layout.z_stride,
                                        "row field offset");
  linear = checked_add(
      linear, checked_multiply(jj, layout.y_stride, "row field offset"),
      "row field offset");
  linear = checked_add(
      linear, checked_multiply(ii, layout.x_stride, "row field offset"),
      "row field offset");
  const std::size_t byte_offset =
      checked_multiply(linear, layout.scalar_bytes, "row byte offset");
  if (byte_offset > layout.field_bytes ||
      row_bytes > layout.field_bytes - byte_offset) {
    throw Error("halo row copy exceeds allocated field storage");
  }
  return byte_offset;
}

template <class CopyRow>
std::size_t for_each_halo_row(const HaloRowLayout& layout, Box3 box,
                              std::size_t wire_bytes, CopyRow&& copy_row) {
  if (layout.components == 0U || layout.scalar_bytes == 0U ||
      layout.x_stride != static_cast<std::size_t>(layout.components)) {
    throw Error("halo row layout is not contiguous in x/component order");
  }
  const std::size_t nx = axis_span(box.begin.x, box.end.x);
  const std::size_t ny = axis_span(box.begin.y, box.end.y);
  const std::size_t nz = axis_span(box.begin.z, box.end.z);
  std::size_t row_bytes = checked_multiply(
      nx, static_cast<std::size_t>(layout.components), "row scalar count");
  row_bytes =
      checked_multiply(row_bytes, layout.scalar_bytes, "row byte size");
  const std::size_t row_count =
      checked_multiply(ny, nz, "row copy count");
  const std::size_t required_wire_bytes =
      checked_multiply(row_count, row_bytes, "row wire byte size");
  if (wire_bytes != required_wire_bytes) {
    throw Error("halo row wire buffer size does not match its region");
  }
  if (row_bytes == 0U || row_count == 0U) {
    return 0U;
  }

  std::size_t cursor = 0U;
  std::size_t events = 0U;
  for (int k = box.begin.z; k < box.end.z; ++k) {
    for (int j = box.begin.y; j < box.end.y; ++j) {
      const std::size_t field_offset =
          row_field_byte_offset(layout, box.begin.x, j, k, row_bytes);
      copy_row(field_offset, cursor, row_bytes);
      cursor = checked_add(cursor, row_bytes, "row wire cursor");
      ++events;
    }
  }
  if (cursor != wire_bytes || events != row_count) {
    throw Error("halo row traversal did not cover its complete region");
  }
  return events;
}

}  // namespace

int halo_offset_code(Int3 offset) {
  if (!legal_offset(offset)) {
    throw Error(
        "halo offset must identify a nonzero direction in the 3x3x3 stencil");
  }
  return (offset.z + 1) * 9 + (offset.y + 1) * 3 + (offset.x + 1);
}

void observe_post_event(PostEventState& state, PostEvent event) noexcept {
  if (event == PostEvent::receive) {
    if (state.send_seen) {
      state.valid = false;
    }
    if (state.receive_posts == std::numeric_limits<std::size_t>::max()) {
      state.valid = false;
    } else {
      ++state.receive_posts;
    }
  } else {
    if (!state.send_seen) {
      state.first_send_sequence = state.sequence;
      state.send_seen = true;
      if (state.receive_posts != state.expected_receive_posts) {
        state.valid = false;
      }
    }
    if (state.send_posts == std::numeric_limits<std::size_t>::max()) {
      state.valid = false;
    } else {
      ++state.send_posts;
    }
  }
  if (state.sequence == std::numeric_limits<std::size_t>::max()) {
    state.valid = false;
  } else {
    ++state.sequence;
  }
}

bool post_event_sequence_valid(const PostEventState& state) noexcept {
  return state.valid &&
         state.receive_posts == state.expected_receive_posts;
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
    throw Error("MPI_TAG_UB=" + std::to_string(*upper_bound) +
                ", required >=26 for halo direction tags");
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

std::size_t pack_halo_region_rows(
    const std::byte* field_data, const HaloRowLayout& layout, Box3 box,
    std::byte* wire_data, std::size_t wire_bytes) {
  if (wire_bytes != 0U && (field_data == nullptr || wire_data == nullptr)) {
    throw Error("halo row pack requires non-null buffers");
  }
  return for_each_halo_row(
      layout, box, wire_bytes,
      [field_data, wire_data](std::size_t field_offset,
                              std::size_t wire_offset,
                              std::size_t row_bytes) {
        std::memcpy(wire_data + wire_offset, field_data + field_offset,
                    row_bytes);
      });
}

std::size_t unpack_halo_region_rows(
    std::byte* field_data, const HaloRowLayout& layout, Box3 box,
    const std::byte* wire_data, std::size_t wire_bytes) {
  if (wire_bytes != 0U && (field_data == nullptr || wire_data == nullptr)) {
    throw Error("halo row unpack requires non-null buffers");
  }
  return for_each_halo_row(
      layout, box, wire_bytes,
      [field_data, wire_data](std::size_t field_offset,
                              std::size_t wire_offset,
                              std::size_t row_bytes) {
        std::memcpy(field_data + field_offset, wire_data + wire_offset,
                    row_bytes);
      });
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

bool completion_succeeded(CompletionOutcome outcome) noexcept {
  return !outcome.mpi_error_seen && outcome.requests_proven_null;
}

CompletionFailureAction completion_failure_action(
    CompletionOutcome outcome) noexcept {
  return outcome.requests_proven_null
             ? CompletionFailureAction::discard_and_throw
             : CompletionFailureAction::terminate_process;
}

FailureRecoveryAction failure_recovery_action(
    bool cancellation_calls_succeeded, bool requests_proven_null,
    bool channel_quiescent, bool context_replaced) noexcept {
  if (!cancellation_calls_succeeded || !requests_proven_null) {
    return FailureRecoveryAction::terminate_process;
  }
  if (channel_quiescent || context_replaced) {
    return FailureRecoveryAction::return_idle;
  }
  return FailureRecoveryAction::replace_context;
}

}  // namespace hundun::runtime::detail
