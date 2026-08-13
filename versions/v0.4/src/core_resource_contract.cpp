// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kAllocationOverflow = 401U;
constexpr std::uint32_t kHaloMessageOverflow = 402U;
constexpr std::uint32_t kHaloByteOverflow = 403U;
constexpr std::uint32_t kNumericRefillOverflow = 404U;
constexpr std::uint32_t kHierarchyRebuildOverflow = 405U;
constexpr std::uint32_t kCachePublishOverflow = 406U;
constexpr std::uint32_t kLinearIterationOverflow = 407U;
constexpr std::uint32_t kStageWallOverflow = 408U;

constexpr std::uint32_t kWorkspaceLimit = 411U;
constexpr std::uint32_t kAllocationLimit = 412U;
constexpr std::uint32_t kHaloMessageLimit = 413U;
constexpr std::uint32_t kHaloByteLimit = 414U;
constexpr std::uint32_t kNumericRefillLimit = 415U;
constexpr std::uint32_t kHierarchyRebuildLimit = 416U;
constexpr std::uint32_t kCachePublishLimit = 417U;
constexpr std::uint32_t kLinearIterationLimit = 418U;
constexpr std::uint32_t kStageWallLimit = 419U;

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& result) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    return false;
  }
  result = left + right;
  return true;
}

}  // namespace

Status add_resource_counters(ResourceCounters& counters,
                             ResourceCounters increment) noexcept {
  ResourceCounters candidate = counters;
  candidate.peak_workspace_bytes =
      std::max(candidate.peak_workspace_bytes, increment.peak_workspace_bytes);

  if (!checked_add(candidate.allocations, increment.allocations,
                   candidate.allocations)) {
    return {StatusCode::invalid_plan, kAllocationOverflow};
  }
  if (!checked_add(candidate.merged_halo_messages,
                   increment.merged_halo_messages,
                   candidate.merged_halo_messages)) {
    return {StatusCode::invalid_plan, kHaloMessageOverflow};
  }
  if (!checked_add(candidate.merged_halo_bytes, increment.merged_halo_bytes,
                   candidate.merged_halo_bytes)) {
    return {StatusCode::invalid_plan, kHaloByteOverflow};
  }
  if (!checked_add(candidate.numeric_refills, increment.numeric_refills,
                   candidate.numeric_refills)) {
    return {StatusCode::invalid_plan, kNumericRefillOverflow};
  }
  if (!checked_add(candidate.hierarchy_rebuilds,
                   increment.hierarchy_rebuilds,
                   candidate.hierarchy_rebuilds)) {
    return {StatusCode::invalid_plan, kHierarchyRebuildOverflow};
  }
  if (!checked_add(candidate.cache_publishes, increment.cache_publishes,
                   candidate.cache_publishes)) {
    return {StatusCode::invalid_plan, kCachePublishOverflow};
  }
  if (!checked_add(candidate.linear_iterations, increment.linear_iterations,
                   candidate.linear_iterations)) {
    return {StatusCode::invalid_plan, kLinearIterationOverflow};
  }
  if (!checked_add(candidate.stage_wall_nanoseconds,
                   increment.stage_wall_nanoseconds,
                   candidate.stage_wall_nanoseconds)) {
    return {StatusCode::invalid_plan, kStageWallOverflow};
  }

  counters = candidate;
  return {};
}

Status validate_resource_counters(const ResourceContract& contract,
                                  const ResourceCounters& counters) noexcept {
  if (counters.peak_workspace_bytes > contract.max_live_workspace_bytes) {
    return {StatusCode::invalid_plan, kWorkspaceLimit};
  }
  if (counters.allocations > contract.allocation_allowance) {
    return {StatusCode::invalid_plan, kAllocationLimit};
  }
  if (counters.merged_halo_messages > contract.merged_halo_messages) {
    return {StatusCode::invalid_plan, kHaloMessageLimit};
  }
  if (counters.merged_halo_bytes > contract.merged_halo_bytes) {
    return {StatusCode::invalid_plan, kHaloByteLimit};
  }
  if (counters.numeric_refills > contract.numeric_refills) {
    return {StatusCode::invalid_plan, kNumericRefillLimit};
  }
  if (counters.hierarchy_rebuilds > contract.hierarchy_rebuilds) {
    return {StatusCode::invalid_plan, kHierarchyRebuildLimit};
  }
  if (counters.cache_publishes > contract.cache_publishes) {
    return {StatusCode::invalid_plan, kCachePublishLimit};
  }
  if (counters.linear_iterations > contract.linear_iterations) {
    return {StatusCode::invalid_plan, kLinearIterationLimit};
  }
  if (counters.stage_wall_nanoseconds > contract.stage_wall_nanoseconds) {
    return {StatusCode::invalid_plan, kStageWallLimit};
  }
  return {};
}

}  // namespace hundun::v04
