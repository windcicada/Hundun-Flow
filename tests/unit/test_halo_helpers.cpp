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

namespace {

using hundun::runtime::Box3;
using hundun::runtime::Error;
using hundun::runtime::Int3;
using hundun::runtime::detail::ActiveDestructionAction;
using hundun::runtime::detail::CompletionFailureAction;
using hundun::runtime::detail::active_destruction_action;
using hundun::runtime::detail::checked_region_payload_bytes;
using hundun::runtime::detail::completion_failure_action;
using hundun::runtime::detail::halo_offset_code;
using hundun::runtime::detail::halo_receive_tag;
using hundun::runtime::detail::reverse_offset;
using hundun::runtime::detail::split_count_ranges;
using hundun::runtime::detail::effective_halo_tag_upper_bound;

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
  expect_error(
      [&] {
        static_cast<void>(effective_halo_tag_upper_bound(true, &too_small));
      },
      "upper bound");
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

  HUNDUN_CHECK(completion_failure_action(true) ==
               CompletionFailureAction::discard_and_throw);
  HUNDUN_CHECK(completion_failure_action(false) ==
               CompletionFailureAction::terminate_process);
}

}  // namespace

int main() {
  static_assert(std::is_trivially_copyable_v<
                hundun::runtime::detail::CountRange>);
  return hundun::test::run([] {
    test_direction_codes_and_tags();
    test_tag_upper_bound();
    test_count_chunking_and_wait_batches();
    test_payload_arithmetic();
    test_cleanup_policies();
  });
}
