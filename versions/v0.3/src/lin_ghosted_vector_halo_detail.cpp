// SPDX-License-Identifier: Apache-2.0

#include "lin_ghosted_vector_halo_detail.hpp"

#include "hundun/rt_error.hpp"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>

namespace hundun::linear::detail {
namespace {

thread_local VectorHaloTestOptions test_options;
thread_local VectorHaloTestSnapshot creation_snapshot;
thread_local VectorHaloTestSnapshot* active_snapshot;

[[noreturn]] void throw_injected_trace_failure() {
  throw runtime::Error("injected device-path test-double failure");
}

void reset_snapshot_preserving_storage(
    VectorHaloTestSnapshot& snapshot) noexcept {
  auto send_wire_identities = std::move(snapshot.send_wire_identities);
  auto receive_wire_identities = std::move(snapshot.receive_wire_identities);
  auto post_events = std::move(snapshot.post_events);
  const bool storage_prepared = snapshot.observation_storage_prepared;
  snapshot = VectorHaloTestSnapshot{};
  snapshot.observation_storage_prepared = storage_prepared;
  snapshot.send_wire_identities = std::move(send_wire_identities);
  snapshot.receive_wire_identities = std::move(receive_wire_identities);
  snapshot.post_events = std::move(post_events);
}

}  // namespace

std::size_t checked_vector_bytes(std::size_t values) {
  if (values > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw runtime::Error("GhostedVector byte-size multiplication overflow");
  }
  return values * sizeof(double);
}

std::vector<CountRange> split_count_ranges(std::size_t total,
                                           std::size_t maximum_count) {
  if (maximum_count == 0U) {
    throw runtime::Error("message chunk limit must be positive");
  }
  if (maximum_count > static_cast<std::size_t>(INT_MAX)) {
    throw runtime::Error("message chunk limit exceeds MPI INT_MAX");
  }
  std::vector<CountRange> result;
  if (total == 0U) {
    return result;
  }
  const std::size_t chunks = total / maximum_count +
                             (total % maximum_count == 0U ? 0U : 1U);
  try {
    result.reserve(chunks);
    for (std::size_t offset = 0U; offset < total;) {
      const std::size_t count = std::min(maximum_count, total - offset);
      result.push_back(CountRange{offset, static_cast<int>(count)});
      offset += count;
    }
  } catch (const std::bad_alloc&) {
    throw runtime::Error("message chunk metadata allocation failed");
  } catch (const std::length_error&) {
    throw runtime::Error("message chunk metadata size is unsupported");
  }
  return result;
}

int effective_tag_upper_bound(bool attribute_present, const int* upper_bound,
                              int required_tag) {
  constexpr int mpi_minimum_tag_upper_bound = 32767;
  if (required_tag < 0) {
    throw runtime::Error("required MPI tag must be nonnegative");
  }
  if (!attribute_present) {
    if (required_tag > mpi_minimum_tag_upper_bound) {
      throw runtime::Error("required MPI tag exceeds the MPI-3 fallback");
    }
    return mpi_minimum_tag_upper_bound;
  }
  if (upper_bound == nullptr) {
    throw runtime::Error("MPI_TAG_UB attribute is present but null");
  }
  if (*upper_bound < required_tag) {
    throw runtime::Error("MPI tag upper bound is too small");
  }
  return *upper_bound;
}

int owner_coordinate(int cells, int processes, int global_coordinate) {
  if (cells <= 0 || processes <= 0 || processes > cells) {
    throw runtime::Error("owner-coordinate partition shape is invalid");
  }
  if (global_coordinate < 0 || global_coordinate >= cells) {
    throw runtime::Error("owner-coordinate global coordinate is invalid");
  }
  const int quotient = cells / processes;
  const int remainder = cells % processes;
  const int larger_span = (quotient + 1) * remainder;
  if (global_coordinate < larger_span) {
    return global_coordinate / (quotient + 1);
  }
  return remainder + (global_coordinate - larger_span) / quotient;
}

bool operator==(const OrderedRequest& left,
                const OrderedRequest& right) noexcept {
  return left.peer == right.peer && left.global_id == right.global_id &&
         left.local_index == right.local_index;
}

std::vector<OrderedRequest> order_requests(
    std::vector<OrderedRequest> requests) {
  std::sort(requests.begin(), requests.end(),
            [](const OrderedRequest& left, const OrderedRequest& right) {
              if (left.peer != right.peer) {
                return left.peer < right.peer;
              }
              if (left.global_id != right.global_id) {
                return left.global_id < right.global_id;
              }
              return left.local_index < right.local_index;
            });
  for (std::size_t index = 1U; index < requests.size(); ++index) {
    if (requests[index - 1U].peer == requests[index].peer &&
        requests[index - 1U].global_id == requests[index].global_id) {
      throw runtime::Error("compact Halo request contains a duplicate ID");
    }
  }
  return requests;
}

BufferHaloPath select_buffer_halo_path(
    const DevicePathCapabilities& capabilities) {
  if (capabilities.space == execution::ExecutionSpace::host) {
    if (capabilities.context_host_access && capabilities.buffer_host_access &&
        capabilities.host_mpi) {
      return BufferHaloPath::host_direct;
    }
    throw runtime::Error("no complete host Buffer Halo path is available");
  }
  if (capabilities.context_device_access &&
      capabilities.buffer_device_lifetime && capabilities.device_mpi) {
    return BufferHaloPath::device_direct;
  }
  if (capabilities.context_device_access &&
      capabilities.buffer_device_lifetime &&
      capabilities.device_to_host_transfer &&
      capabilities.host_to_device_transfer &&
      capabilities.host_staging_allocation && capabilities.host_mpi &&
      capabilities.event_lifetime) {
    return BufferHaloPath::device_host_staged;
  }
  throw runtime::Error("no complete device Buffer Halo path is available");
}

void run_device_path_script(
    const DevicePathCapabilities& capabilities,
    std::optional<DeviceTraceStep> injected_failure,
    std::vector<DeviceTraceStep>& trace) {
  trace.clear();
  const BufferHaloPath path = select_buffer_halo_path(capabilities);
  const auto record = [&](DeviceTraceStep step) {
    trace.push_back(step);
    if (injected_failure.has_value() && *injected_failure == step) {
      throw_injected_trace_failure();
    }
  };
  if (path == BufferHaloPath::device_host_staged) {
    record(DeviceTraceStep::device_to_host_submit);
    record(DeviceTraceStep::pre_communication_event_wait);
  }
  record(DeviceTraceStep::mpi_begin);
  record(DeviceTraceStep::mpi_wait);
  if (path == BufferHaloPath::device_host_staged) {
    record(DeviceTraceStep::host_to_device_submit);
    record(DeviceTraceStep::post_communication_event_wait);
  }
}

NonblockingPostIssueAction nonblocking_post_issue_action(
    NonblockingPostIssueOrigin origin) noexcept {
  return origin == NonblockingPostIssueOrigin::synthetic_before_call
             ? NonblockingPostIssueAction::recover_known_prefix
             : NonblockingPostIssueAction::terminate_process;
}

std::size_t injection_selection_collective_count(
    int configured_rank) noexcept {
  return configured_rank == -1 ? 0U : 1U;
}

void set_vector_halo_test_options(VectorHaloTestOptions options) noexcept {
  test_options = options;
}

VectorHaloTestOptions current_vector_halo_test_options() noexcept {
  return test_options;
}

void reset_vector_halo_test_observation() noexcept {
  reset_snapshot_preserving_storage(
      active_snapshot == nullptr ? creation_snapshot : *active_snapshot);
}

void begin_vector_halo_creation_observation() noexcept {
  active_snapshot = &creation_snapshot;
  creation_snapshot = VectorHaloTestSnapshot{};
}

void prepare_vector_halo_creation_observation(std::size_t peer_count,
                                              std::size_t request_count) {
  creation_snapshot.send_wire_identities.resize(peer_count);
  creation_snapshot.receive_wire_identities.resize(peer_count);
  creation_snapshot.post_events.resize(request_count);
  creation_snapshot.observation_storage_prepared = true;
}

VectorHaloTestSnapshot take_vector_halo_creation_observation() noexcept {
  active_snapshot = nullptr;
  return std::move(creation_snapshot);
}

void activate_vector_halo_test_observation(
    VectorHaloTestSnapshot* snapshot) noexcept {
  active_snapshot = snapshot;
}

void deactivate_vector_halo_test_observation(
    const VectorHaloTestSnapshot* snapshot) noexcept {
  if (active_snapshot == snapshot) {
    active_snapshot = nullptr;
  }
}

bool vector_halo_test_observation_storage_prepared() noexcept {
  return active_snapshot != nullptr &&
         active_snapshot->observation_storage_prepared;
}

VectorHaloTestSnapshot vector_halo_test_snapshot() {
  return active_snapshot == nullptr ? creation_snapshot : *active_snapshot;
}

VectorHaloTestSnapshot& mutable_vector_halo_test_snapshot() noexcept {
  return active_snapshot == nullptr ? creation_snapshot : *active_snapshot;
}

}  // namespace hundun::linear::detail
