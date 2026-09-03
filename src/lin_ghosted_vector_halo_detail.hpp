// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/exec_execution.hpp"
#include "hundun/lin_ghosted_vector_halo.hpp"
#include "hundun/mesh_topology.hpp"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace hundun::linear::detail {

struct CountRange {
  std::size_t offset{};
  int count{};
};

std::size_t checked_vector_bytes(std::size_t values);
std::vector<CountRange> split_count_ranges(std::size_t total,
                                           std::size_t maximum_count);
int effective_tag_upper_bound(bool attribute_present, const int* upper_bound,
                              int required_tag);
int owner_coordinate(int cells, int processes, int global_coordinate);

struct OrderedRequest {
  int peer{};
  mesh::GlobalCellId global_id{};
  std::size_t local_index{};
};

bool operator==(const OrderedRequest& left,
                const OrderedRequest& right) noexcept;
std::vector<OrderedRequest> order_requests(
    std::vector<OrderedRequest> requests);

struct DevicePathCapabilities {
  execution::ExecutionSpace space{execution::ExecutionSpace::host};
  bool context_host_access{};
  bool buffer_host_access{};
  bool host_mpi{};
  bool context_device_access{};
  bool buffer_device_lifetime{};
  bool device_mpi{};
  bool device_to_host_transfer{};
  bool host_to_device_transfer{};
  bool host_staging_allocation{};
  bool event_lifetime{};
};

BufferHaloPath select_buffer_halo_path(
    const DevicePathCapabilities& capabilities);

enum class DeviceTraceStep {
  device_to_host_submit,
  pre_communication_event_wait,
  mpi_begin,
  mpi_wait,
  host_to_device_submit,
  post_communication_event_wait
};

void run_device_path_script(
    const DevicePathCapabilities& capabilities,
    std::optional<DeviceTraceStep> injected_failure,
    std::vector<DeviceTraceStep>& trace);

enum class NonblockingPostIssueOrigin {
  synthetic_before_call,
  mpi_call_error
};

enum class NonblockingPostIssueAction {
  recover_known_prefix,
  terminate_process
};

NonblockingPostIssueAction nonblocking_post_issue_action(
    NonblockingPostIssueOrigin origin) noexcept;

inline constexpr int kInjectAllEligibleRanks = -2;

std::size_t injection_selection_collective_count(
    int configured_rank) noexcept;

struct VectorHaloTestOptions {
  std::size_t chunk_limit{static_cast<std::size_t>(INT_MAX)};
  bool observe{};
  int inject_request_id_mismatch_rank{-1};
  int inject_metadata_post_failure_rank{-1};
  int inject_metadata_completion_failure_rank{-1};
  int inject_post_failure_rank{-1};
  int inject_completion_failure_rank{-1};
#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
  int inject_target_copy_failure_rank{-1};
#endif
};

struct WirePostEvent {
  bool receive{};
  int peer{};
  std::size_t offset{};
  int count{};
  int tag{};
};

struct FailureDiagnosticSnapshot {
  bool valid{};
  int category{};
  int rank{-1};
  int operation{};
  int result{};
  int peer{-1};
  std::size_t value_offset{};
  int value_count{};
  int tag{-1};
};

struct VectorHaloTestSnapshot {
  bool observation_storage_prepared{};
  std::size_t injection_selection_collectives{};
  bool receives_preceded_sends{true};
  bool chunk_offsets_ordered{true};
  std::size_t receive_posts{};
  std::size_t send_posts{};
  std::size_t request_capacity{};
  std::size_t context_replacements{};
  std::size_t metadata_context_replacements{};
  std::size_t metadata_post_calls{};
  std::size_t metadata_wait_calls{};
  std::size_t metadata_posts_before_failure{};
  std::size_t metadata_completion_prefix{};
  std::size_t metadata_non_null_before_cleanup{};
  std::size_t metadata_non_null_after_cleanup{};
  std::size_t runtime_posts_before_failure{};
  std::size_t runtime_completion_prefix{};
  std::size_t runtime_non_null_before_cleanup{};
  std::size_t runtime_non_null_after_cleanup{};
  FailureDiagnosticSnapshot failure;
  std::vector<execution::AllocationIdentity> send_wire_identities;
  std::vector<execution::AllocationIdentity> receive_wire_identities;
  std::vector<WirePostEvent> post_events;
};

void set_vector_halo_test_options(VectorHaloTestOptions options) noexcept;
VectorHaloTestOptions current_vector_halo_test_options() noexcept;
void reset_vector_halo_test_observation() noexcept;
void begin_vector_halo_creation_observation() noexcept;
void prepare_vector_halo_creation_observation(std::size_t peer_count,
                                              std::size_t request_count);
VectorHaloTestSnapshot take_vector_halo_creation_observation() noexcept;
void activate_vector_halo_test_observation(
    VectorHaloTestSnapshot* snapshot) noexcept;
void deactivate_vector_halo_test_observation(
    const VectorHaloTestSnapshot* snapshot) noexcept;
bool vector_halo_test_observation_storage_prepared() noexcept;
VectorHaloTestSnapshot vector_halo_test_snapshot();
VectorHaloTestSnapshot& mutable_vector_halo_test_snapshot() noexcept;

}  // namespace hundun::linear::detail
