// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/types.hpp"
#include "runtime/src/halo_detail.hpp"
#include "tests/support/test_main.hpp"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using hundun::runtime::Box3;
using hundun::runtime::Error;
using hundun::runtime::Int3;
using hundun::runtime::detail::ActiveDestructionAction;
using hundun::runtime::detail::CompletionOutcome;
using hundun::runtime::detail::CompletionFailureAction;
using hundun::runtime::detail::FailureRecoveryAction;
using hundun::runtime::detail::HaloRowLayout;
using hundun::runtime::detail::HaloFailureRecord;
using hundun::runtime::detail::PostEvent;
using hundun::runtime::detail::PostEventState;
using hundun::runtime::detail::active_destruction_action;
using hundun::runtime::detail::checked_region_payload_bytes;
using hundun::runtime::detail::completion_failure_action;
using hundun::runtime::detail::completion_succeeded;
using hundun::runtime::detail::halo_offset_code;
using hundun::runtime::detail::halo_receive_tag;
using hundun::runtime::detail::pack_halo_region_rows;
using hundun::runtime::detail::observe_post_event;
using hundun::runtime::detail::post_event_sequence_valid;
using hundun::runtime::detail::reverse_offset;
using hundun::runtime::detail::split_count_ranges;
using hundun::runtime::detail::effective_halo_tag_upper_bound;
using hundun::runtime::detail::failure_recovery_action;
using hundun::runtime::detail::unpack_halo_region_rows;

template <class Function>
void expect_error(Function&& function, const std::string& text) {
  bool threw = false;
  try {
    function();
  } catch (const Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()).find(text) != std::string::npos);
  }
  HUNDUN_CHECK(threw);
}

void test_direction_codes_and_tags() {
  bool seen[27]{};
  int count = 0;
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        const Int3 offset{x, y, z};
        if (x == 0 && y == 0 && z == 0) {
          expect_error([&] { static_cast<void>(halo_offset_code(offset)); },
                       "nonzero");
          continue;
        }
        const int code = halo_offset_code(offset);
        HUNDUN_CHECK(code >= 0);
        HUNDUN_CHECK(code <= 26);
        HUNDUN_CHECK(code != 13);
        HUNDUN_CHECK(!seen[code]);
        seen[code] = true;
        ++count;

        const Int3 reversed = reverse_offset(offset);
        HUNDUN_CHECK(reversed.x == -x);
        HUNDUN_CHECK(reversed.y == -y);
        HUNDUN_CHECK(reversed.z == -z);
        HUNDUN_CHECK(halo_receive_tag(offset) ==
                     halo_offset_code(reversed));
        HUNDUN_CHECK(reverse_offset(reversed).x == x);
        HUNDUN_CHECK(reverse_offset(reversed).y == y);
        HUNDUN_CHECK(reverse_offset(reversed).z == z);
      }
    }
  }
  HUNDUN_CHECK(count == 26);
  for (int code = 0; code <= 26; ++code) {
    HUNDUN_CHECK(seen[code] == (code != 13));
  }

  expect_error(
      [] { static_cast<void>(halo_offset_code(Int3{-2, 0, 0})); },
      "offset");
  expect_error(
      [] { static_cast<void>(halo_offset_code(Int3{0, 2, 0})); },
      "offset");
  expect_error(
      [] { static_cast<void>(reverse_offset(Int3{0, 0, 0})); },
      "nonzero");
}

void test_tag_upper_bound() {
  HUNDUN_CHECK(effective_halo_tag_upper_bound(false, nullptr) == 32767);
  expect_error(
      [] {
        static_cast<void>(effective_halo_tag_upper_bound(true, nullptr));
      },
      "null");

  const int too_small = 25;
  const int minimum_required = 26;
  const int maximum = INT_MAX;
  bool detailed_bound_error = false;
  try {
    static_cast<void>(effective_halo_tag_upper_bound(true, &too_small));
  } catch (const Error& error) {
    const std::string message = error.what();
    detailed_bound_error =
        message.find("MPI_TAG_UB=25") != std::string::npos &&
        message.find("required >=26") != std::string::npos;
  }
  HUNDUN_CHECK(detailed_bound_error);
  HUNDUN_CHECK(effective_halo_tag_upper_bound(
                   true, &minimum_required) == minimum_required);
  HUNDUN_CHECK(effective_halo_tag_upper_bound(true, &maximum) == maximum);
}

void check_range(const hundun::runtime::detail::CountRange& range,
                 std::size_t offset, int count) {
  HUNDUN_CHECK(range.offset == offset);
  HUNDUN_CHECK(range.count == count);
}

void test_count_chunking_and_wait_batches() {
  const auto empty = split_count_ranges(0U, INT_MAX, "byte chunks");
  HUNDUN_CHECK(empty.empty());

  const auto one = split_count_ranges(
      static_cast<std::size_t>(INT_MAX), INT_MAX, "byte chunks");
  HUNDUN_CHECK(one.size() == 1U);
  check_range(one[0], 0U, INT_MAX);

  const auto two = split_count_ranges(
      static_cast<std::size_t>(INT_MAX) + 1U, INT_MAX, "byte chunks");
  HUNDUN_CHECK(two.size() == 2U);
  check_range(two[0], 0U, INT_MAX);
  check_range(two[1], static_cast<std::size_t>(INT_MAX), 1);

  const auto small = split_count_ranges(17U, 7U, "byte chunks");
  HUNDUN_CHECK(small.size() == 3U);
  check_range(small[0], 0U, 7);
  check_range(small[1], 7U, 7);
  check_range(small[2], 14U, 3);

  const auto wait_batches = split_count_ranges(
      static_cast<std::size_t>(INT_MAX) + 1U, INT_MAX,
      "MPI_Waitall batches");
  HUNDUN_CHECK(wait_batches.size() == 2U);
  check_range(wait_batches[0], 0U, INT_MAX);
  check_range(wait_batches[1], static_cast<std::size_t>(INT_MAX), 1);

  expect_error(
      [] { static_cast<void>(split_count_ranges(1U, 0U, "ranges")); },
      "positive");
  expect_error(
      [] {
        static_cast<void>(split_count_ranges(
            1U, static_cast<std::size_t>(INT_MAX) + 1U, "ranges"));
      },
      "MPI int");

  // This input would require an impossible local metadata allocation. It must
  // be rejected by checked arithmetic/capacity guards, not wrapped.
  expect_error(
      [] {
        static_cast<void>(split_count_ranges(
            std::numeric_limits<std::size_t>::max(), 1U, "ranges"));
      },
      "range metadata");
}

void test_payload_arithmetic() {
  HUNDUN_CHECK(checked_region_payload_bytes(
                   Box3{Int3{0, 0, 0}, Int3{2, 3, 4}}, 3U,
                   sizeof(double)) ==
               2U * 3U * 4U * 3U * sizeof(double));
  HUNDUN_CHECK(checked_region_payload_bytes(
                   Box3{Int3{5, 0, 0}, Int3{5, 7, 9}}, 1U, 1U) == 0U);

  expect_error(
      [] {
        static_cast<void>(checked_region_payload_bytes(
            Box3{Int3{2, 0, 0}, Int3{1, 1, 1}}, 1U, 1U));
      },
      "box");
  expect_error(
      [] {
        static_cast<void>(checked_region_payload_bytes(
            Box3{Int3{0, 0, 0}, Int3{1, 1, 1}}, 0U, 1U));
      },
      "component");
  expect_error(
      [] {
        static_cast<void>(checked_region_payload_bytes(
            Box3{Int3{0, 0, 0}, Int3{1, 1, 1}}, 1U, 0U));
      },
      "scalar");
  expect_error(
      [] {
        static_cast<void>(checked_region_payload_bytes(
            Box3{Int3{0, 0, 0},
                 Int3{std::numeric_limits<int>::max(),
                      std::numeric_limits<int>::max(),
                      std::numeric_limits<int>::max()}},
            std::numeric_limits<std::uint32_t>::max(), sizeof(double)));
      },
      "overflow");
}

void test_cleanup_policies() {
  HUNDUN_CHECK(active_destruction_action(false, false, false) ==
               ActiveDestructionAction::clear_without_mpi);
  HUNDUN_CHECK(active_destruction_action(false, true, false) ==
               ActiveDestructionAction::clear_without_mpi);
  HUNDUN_CHECK(active_destruction_action(true, false, true) ==
               ActiveDestructionAction::clear_without_mpi);
  HUNDUN_CHECK(active_destruction_action(true, true, true) ==
               ActiveDestructionAction::clear_without_mpi);
  HUNDUN_CHECK(active_destruction_action(true, true, false) ==
               ActiveDestructionAction::drain_requests);
  HUNDUN_CHECK(active_destruction_action(true, false, false) ==
               ActiveDestructionAction::terminate_process);

  const CompletionOutcome successful_completion{false, true};
  HUNDUN_CHECK(completion_succeeded(successful_completion));

  const CompletionOutcome failed_but_proven_complete{true, true};
  HUNDUN_CHECK(!completion_succeeded(failed_but_proven_complete));
  HUNDUN_CHECK(completion_failure_action(failed_but_proven_complete) ==
               CompletionFailureAction::discard_and_throw);

  const CompletionOutcome failed_and_unproven{true, false};
  HUNDUN_CHECK(!completion_succeeded(failed_and_unproven));
  HUNDUN_CHECK(completion_failure_action(failed_and_unproven) ==
               CompletionFailureAction::terminate_process);
}

void test_failure_recovery_requires_a_quiescent_or_replaced_context() {
  HUNDUN_CHECK(failure_recovery_action(true, true, false, false) ==
               FailureRecoveryAction::replace_context);
  HUNDUN_CHECK(failure_recovery_action(true, true, true, false) ==
               FailureRecoveryAction::return_idle);
  HUNDUN_CHECK(failure_recovery_action(true, true, false, true) ==
               FailureRecoveryAction::return_idle);
  HUNDUN_CHECK(failure_recovery_action(true, false, false, false) ==
               FailureRecoveryAction::terminate_process);
  HUNDUN_CHECK(failure_recovery_action(false, true, false, false) ==
               FailureRecoveryAction::terminate_process);
}

std::size_t test_field_byte_offset(const HaloRowLayout& layout, int i, int j,
                                   int k, std::uint32_t component) {
  const auto ii = static_cast<std::size_t>(i + layout.ghost_width);
  const auto jj = static_cast<std::size_t>(j + layout.ghost_width);
  const auto kk = static_cast<std::size_t>(k + layout.ghost_width);
  return (kk * layout.z_stride + jj * layout.y_stride +
          ii * layout.x_stride + static_cast<std::size_t>(component)) *
         layout.scalar_bytes;
}

void test_row_wise_face_edge_corner_copy() {
  constexpr Int3 extent{3, 2, 2};
  constexpr int ghost = 1;
  constexpr std::uint32_t components = 3U;
  constexpr std::size_t scalar_bytes = 2U;
  constexpr std::size_t padded_x = 5U;
  constexpr std::size_t padded_y = 4U;
  constexpr std::size_t padded_z = 4U;
  const HaloRowLayout layout{
      ghost,
      components,
      scalar_bytes,
      static_cast<std::size_t>(components),
      padded_x * static_cast<std::size_t>(components),
      padded_y * padded_x * static_cast<std::size_t>(components),
      padded_z * padded_y * padded_x *
          static_cast<std::size_t>(components) * scalar_bytes};
  std::vector<std::byte> source(layout.field_bytes);
  for (std::size_t index = 0U; index < source.size(); ++index) {
    source[index] = static_cast<std::byte>((index * 37U + 11U) % 251U);
  }

  const std::vector<Box3> boxes{
      Box3{Int3{0, 0, 0}, Int3{1, extent.y, extent.z}},
      Box3{Int3{0, 0, 0}, Int3{extent.x, 1, 1}},
      Box3{Int3{extent.x - 1, extent.y - 1, extent.z - 1}, extent}};
  for (const Box3 box : boxes) {
    const std::size_t payload = checked_region_payload_bytes(
        box, components, scalar_bytes);
    std::vector<std::byte> wire(payload);
    const std::size_t pack_events = pack_halo_region_rows(
        source.data(), layout, box, wire.data(), wire.size());
    const std::size_t expected_rows =
        static_cast<std::size_t>(box.end.y - box.begin.y) *
        static_cast<std::size_t>(box.end.z - box.begin.z);
    HUNDUN_CHECK(pack_events == expected_rows);

    std::size_t cursor = 0U;
    for (int k = box.begin.z; k < box.end.z; ++k) {
      for (int j = box.begin.y; j < box.end.y; ++j) {
        for (int i = box.begin.x; i < box.end.x; ++i) {
          for (std::uint32_t component = 0U; component < components;
               ++component) {
            const std::size_t offset =
                test_field_byte_offset(layout, i, j, k, component);
            for (std::size_t byte = 0U; byte < scalar_bytes; ++byte) {
              HUNDUN_CHECK(wire[cursor] == source[offset + byte]);
              wire[cursor] ^= static_cast<std::byte>(0x5aU);
              ++cursor;
            }
          }
        }
      }
    }
    HUNDUN_CHECK(cursor == wire.size());

    std::vector<std::byte> destination(
        layout.field_bytes, static_cast<std::byte>(0xeeU));
    const std::size_t unpack_events = unpack_halo_region_rows(
        destination.data(), layout, box, wire.data(), wire.size());
    HUNDUN_CHECK(unpack_events == expected_rows);
    cursor = 0U;
    for (int k = -ghost; k < extent.z + ghost; ++k) {
      for (int j = -ghost; j < extent.y + ghost; ++j) {
        for (int i = -ghost; i < extent.x + ghost; ++i) {
          const bool in_box = i >= box.begin.x && i < box.end.x &&
                              j >= box.begin.y && j < box.end.y &&
                              k >= box.begin.z && k < box.end.z;
          for (std::uint32_t component = 0U; component < components;
               ++component) {
            const std::size_t offset =
                test_field_byte_offset(layout, i, j, k, component);
            for (std::size_t byte = 0U; byte < scalar_bytes; ++byte) {
              const std::byte expected =
                  in_box ? wire[cursor++] : static_cast<std::byte>(0xeeU);
              HUNDUN_CHECK(destination[offset + byte] == expected);
            }
          }
        }
      }
    }
    HUNDUN_CHECK(cursor == wire.size());
  }
}

PostEventState observe_events(std::initializer_list<PostEvent> events) {
  PostEventState state;
  state.expected_receive_posts = 2U;
  for (const PostEvent event : events) {
    observe_post_event(state, event);
  }
  return state;
}

void test_receive_before_send_event_observer() {
  const PostEventState rrss =
      observe_events({PostEvent::receive, PostEvent::receive,
                      PostEvent::send, PostEvent::send});
  HUNDUN_CHECK(post_event_sequence_valid(rrss));
  HUNDUN_CHECK(rrss.first_send_sequence == 2U);

  const PostEventState rsrs =
      observe_events({PostEvent::receive, PostEvent::send,
                      PostEvent::receive, PostEvent::send});
  HUNDUN_CHECK(!post_event_sequence_valid(rsrs));

  const PostEventState ssrr =
      observe_events({PostEvent::send, PostEvent::send,
                      PostEvent::receive, PostEvent::receive});
  HUNDUN_CHECK(!post_event_sequence_valid(ssrr));
}

}  // namespace

int main() {
  static_assert(std::is_trivially_copyable_v<
                hundun::runtime::detail::CountRange>);
  static_assert(std::is_trivially_copyable_v<HaloFailureRecord>);
  static_assert(sizeof(HaloFailureRecord) == 8U * sizeof(std::int64_t));
  return hundun::test::run([] {
    test_direction_codes_and_tags();
    test_tag_upper_bound();
    test_count_chunking_and_wait_batches();
    test_payload_arithmetic();
    test_cleanup_policies();
    test_failure_recovery_requires_a_quiescent_or_replaced_context();
    test_row_wise_face_edge_corner_copy();
    test_receive_before_send_event_observer();
  });
}
