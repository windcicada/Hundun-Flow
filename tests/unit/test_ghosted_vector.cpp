// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/linear/ghosted_vector_halo.hpp"
#include "hundun/runtime/error.hpp"
#include "linear/src/ghosted_vector_halo_detail.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::execution::CpuReferenceContext;
using hundun::execution::ExecutionSpace;
using hundun::linear::BufferHaloPath;
using hundun::linear::GhostedVector;
using hundun::linear::VectorLayout;
using hundun::linear::detail::DevicePathCapabilities;
using hundun::linear::detail::DeviceTraceStep;
using hundun::linear::detail::NonblockingPostIssueAction;
using hundun::linear::detail::NonblockingPostIssueOrigin;
using hundun::linear::detail::OrderedRequest;
using hundun::runtime::Error;

template <class Function>
std::string expect_error(Function&& function) {
  try {
    function();
  } catch (const Error& error) {
    return error.what();
  }
  throw std::runtime_error("expected hundun::runtime::Error");
}

template <class Function>
void expect_error_containing(Function&& function, const std::string& text) {
  const std::string message = expect_error(std::forward<Function>(function));
  HUNDUN_CHECK(message.find(text) != std::string::npos);
}

std::vector<double> values(hundun::execution::VectorView<const double> view) {
  std::vector<double> result;
  result.reserve(view.size());
  for (std::size_t index = 0; index < view.size(); ++index) {
    result.push_back(view[index]);
  }
  return result;
}

void test_layout_and_vector() {
  VectorLayout empty;
  HUNDUN_CHECK(empty.owned_count() == 0U);
  HUNDUN_CHECK(empty.ghost_count() == 0U);
  HUNDUN_CHECK(empty.local_count() == 0U);
  HUNDUN_CHECK(empty.global_ids().empty());
  HUNDUN_CHECK(empty == VectorLayout{});

  const VectorLayout layout(3U, {11U, 14U, 19U, 22U, 31U});
  HUNDUN_CHECK(layout.owned_count() == 3U);
  HUNDUN_CHECK(layout.ghost_count() == 2U);
  HUNDUN_CHECK(layout.local_count() == 5U);
  HUNDUN_CHECK(layout.global_ids() ==
               std::vector<hundun::mesh::GlobalCellId>(
                   {11U, 14U, 19U, 22U, 31U}));
  VectorLayout cheap_copy = layout;
  HUNDUN_CHECK(cheap_copy == layout);
  HUNDUN_CHECK(cheap_copy.global_ids().data() == layout.global_ids().data());
  HUNDUN_CHECK(VectorLayout(3U, {11U, 14U, 19U, 22U, 32U}) != layout);
  HUNDUN_CHECK(VectorLayout(2U, {11U, 14U, 19U, 22U, 31U}) != layout);
  expect_error_containing([&] { VectorLayout invalid(2U, {1U}); }, "owned");
  expect_error_containing([&] { VectorLayout invalid(1U, {1U, 1U}); },
                          "duplicate");

  CpuReferenceContext context;
  GhostedVector vector(context, layout);
  HUNDUN_CHECK(vector.layout() == layout);
  HUNDUN_CHECK(vector.owned_count() == 3U);
  HUNDUN_CHECK(vector.ghost_count() == 2U);
  HUNDUN_CHECK(vector.local_count() == 5U);
  HUNDUN_CHECK(vector.allocation_identity() != 0U);
  HUNDUN_CHECK(vector.epoch() != 0U);
  HUNDUN_CHECK(vector.backend_identity() == context.backend_identity());
  HUNDUN_CHECK(vector.space() == ExecutionSpace::host);
  HUNDUN_CHECK(vector.local_view().allocation_identity() ==
               vector.allocation_identity());
  HUNDUN_CHECK(vector.owned_view().offset_bytes() == 0U);
  HUNDUN_CHECK(vector.ghost_view().offset_bytes() == 3U * sizeof(double));
  HUNDUN_CHECK(vector.owned_view().size() == 3U);
  HUNDUN_CHECK(vector.ghost_view().size() == 2U);
  auto local = vector.local_view();
  for (std::size_t index = 0; index < local.size(); ++index) {
    local[index] = 100.0 + static_cast<double>(index);
  }
  HUNDUN_CHECK(values(vector.owned_view()) ==
               std::vector<double>({100.0, 101.0, 102.0}));
  HUNDUN_CHECK(values(vector.ghost_view()) ==
               std::vector<double>({103.0, 104.0}));
  const GhostedVector& constant = vector;
  static_assert(std::is_same_v<decltype(constant.local_view()),
                               hundun::execution::VectorView<const double>>);
  HUNDUN_CHECK(!constant.local_view().writable());

  auto surviving = vector.local_view();
  const auto identity = vector.allocation_identity();
  GhostedVector moved(std::move(vector));
  HUNDUN_CHECK(moved.allocation_identity() == identity);
  HUNDUN_CHECK(values(surviving) ==
               std::vector<double>({100.0, 101.0, 102.0, 103.0, 104.0}));
  HUNDUN_CHECK(vector.local_count() == 0U);
  HUNDUN_CHECK(vector.layout() == VectorLayout{});
  expect_error_containing([&] { static_cast<void>(vector.local_view()); },
                          "moved-from");

  GhostedVector destination(context, VectorLayout(1U, {90U}));
  auto stale_destination = destination.local_view();
  auto source_view = moved.local_view();
  destination = std::move(moved);
  expect_error_containing(
      [&] { static_cast<void>(stale_destination.data()); }, "live");
  HUNDUN_CHECK(values(source_view) ==
               std::vector<double>({100.0, 101.0, 102.0, 103.0, 104.0}));
  HUNDUN_CHECK(destination.layout() == layout);
  HUNDUN_CHECK(moved.layout() == VectorLayout{});

  GhostedVector zero(context, VectorLayout{});
  HUNDUN_CHECK(zero.local_view().data() == nullptr);
  HUNDUN_CHECK(zero.owned_view().data() == nullptr);
  HUNDUN_CHECK(zero.ghost_view().data() == nullptr);

  HUNDUN_CHECK(hundun::linear::detail::checked_vector_bytes(0U) == 0U);
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::linear::detail::checked_vector_bytes(
            std::numeric_limits<std::size_t>::max()));
      },
      "overflow");
}

DevicePathCapabilities complete_device_capabilities() {
  DevicePathCapabilities capabilities{};
  capabilities.space = ExecutionSpace::device;
  capabilities.context_device_access = true;
  capabilities.buffer_device_lifetime = true;
  capabilities.device_mpi = true;
  capabilities.device_to_host_transfer = true;
  capabilities.host_to_device_transfer = true;
  capabilities.host_staging_allocation = true;
  capabilities.host_mpi = true;
  capabilities.event_lifetime = true;
  return capabilities;
}

void test_path_selection_and_trace() {
  DevicePathCapabilities host{};
  host.space = ExecutionSpace::host;
  host.context_host_access = true;
  host.buffer_host_access = true;
  host.host_mpi = true;
  HUNDUN_CHECK(hundun::linear::detail::select_buffer_halo_path(host) ==
               BufferHaloPath::host_direct);
  host.host_mpi = false;
  expect_error_containing(
      [&] { static_cast<void>(
                hundun::linear::detail::select_buffer_halo_path(host)); },
      "path");

  auto device = complete_device_capabilities();
  HUNDUN_CHECK(hundun::linear::detail::select_buffer_halo_path(device) ==
               BufferHaloPath::device_direct);
  device.device_mpi = false;
  HUNDUN_CHECK(hundun::linear::detail::select_buffer_halo_path(device) ==
               BufferHaloPath::device_host_staged);
  device.host_to_device_transfer = false;
  expect_error_containing(
      [&] { static_cast<void>(
                hundun::linear::detail::select_buffer_halo_path(device)); },
      "path");

  device = complete_device_capabilities();
  device.device_mpi = false;
  std::vector<DeviceTraceStep> trace;
  hundun::linear::detail::run_device_path_script(device, std::nullopt, trace);
  HUNDUN_CHECK(trace == std::vector<DeviceTraceStep>(
                            {DeviceTraceStep::device_to_host_submit,
                             DeviceTraceStep::pre_communication_event_wait,
                             DeviceTraceStep::mpi_begin,
                             DeviceTraceStep::mpi_wait,
                             DeviceTraceStep::host_to_device_submit,
                             DeviceTraceStep::post_communication_event_wait}));
  for (const DeviceTraceStep failure : trace) {
    std::vector<DeviceTraceStep> failed_trace;
    expect_error_containing(
        [&] {
          hundun::linear::detail::run_device_path_script(
              device, failure, failed_trace);
        },
        "injected");
    const auto found = std::find(trace.begin(), trace.end(), failure);
    const auto expected_count =
        static_cast<std::size_t>(found - trace.begin()) + 1U;
    HUNDUN_CHECK(failed_trace.size() == expected_count);
    HUNDUN_CHECK(failed_trace.back() == failure);
  }
}

void test_pure_plan_helpers() {
  HUNDUN_CHECK(
      hundun::linear::detail::injection_selection_collective_count(-1) ==
      0U);
  HUNDUN_CHECK(
      hundun::linear::detail::injection_selection_collective_count(0) ==
      1U);
  HUNDUN_CHECK(
      hundun::linear::detail::injection_selection_collective_count(
          hundun::linear::detail::kInjectAllEligibleRanks) == 1U);

  HUNDUN_CHECK(hundun::linear::detail::nonblocking_post_issue_action(
                   NonblockingPostIssueOrigin::synthetic_before_call) ==
               NonblockingPostIssueAction::recover_known_prefix);
  HUNDUN_CHECK(hundun::linear::detail::nonblocking_post_issue_action(
                   NonblockingPostIssueOrigin::mpi_call_error) ==
               NonblockingPostIssueAction::terminate_process);

  const auto chunks = hundun::linear::detail::split_count_ranges(8U, 3U);
  HUNDUN_CHECK(chunks.size() == 3U);
  HUNDUN_CHECK(chunks[0].offset == 0U && chunks[0].count == 3);
  HUNDUN_CHECK(chunks[1].offset == 3U && chunks[1].count == 3);
  HUNDUN_CHECK(chunks[2].offset == 6U && chunks[2].count == 2);
  HUNDUN_CHECK(hundun::linear::detail::split_count_ranges(0U, 1U).empty());
  expect_error_containing(
      [&] { static_cast<void>(
                hundun::linear::detail::split_count_ranges(1U, 0U)); },
      "positive");
  expect_error_containing(
      [&] { static_cast<void>(hundun::linear::detail::split_count_ranges(
                1U, static_cast<std::size_t>(INT_MAX) + 1U)); },
      "INT_MAX");

  const int upper = 4096;
  HUNDUN_CHECK(hundun::linear::detail::effective_tag_upper_bound(
                   true, &upper, 7) == upper);
  HUNDUN_CHECK(hundun::linear::detail::effective_tag_upper_bound(
                   false, nullptr, 7) >= 7);
  const int too_small = 6;
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::linear::detail::effective_tag_upper_bound(
            true, &too_small, 7));
      },
      "tag");

  HUNDUN_CHECK(hundun::linear::detail::owner_coordinate(10, 3, 0) == 0);
  HUNDUN_CHECK(hundun::linear::detail::owner_coordinate(10, 3, 3) == 0);
  HUNDUN_CHECK(hundun::linear::detail::owner_coordinate(10, 3, 4) == 1);
  HUNDUN_CHECK(hundun::linear::detail::owner_coordinate(10, 3, 9) == 2);
  expect_error_containing(
      [&] { static_cast<void>(
                hundun::linear::detail::owner_coordinate(10, 3, 10)); },
      "coordinate");

  const auto ordered = hundun::linear::detail::order_requests(
      {{2, 8U, 4U}, {1, 9U, 7U}, {2, 3U, 2U}});
  HUNDUN_CHECK(ordered == std::vector<OrderedRequest>(
                              {{1, 9U, 7U}, {2, 3U, 2U}, {2, 8U, 4U}}));
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::linear::detail::order_requests(
            {{2, 8U, 4U}, {2, 8U, 5U}}));
      },
      "duplicate");
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_layout_and_vector();
    test_path_selection_and_trace();
    test_pure_plan_helpers();
  });
}
