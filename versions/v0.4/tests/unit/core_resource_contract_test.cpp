// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_execution.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using hundun::v04::ResourceContract;
using hundun::v04::ResourceCounters;
using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::add_resource_counters;
using hundun::v04::validate_resource_counters;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool equal(ResourceCounters left, ResourceCounters right) noexcept {
  return left.peak_workspace_bytes == right.peak_workspace_bytes &&
         left.allocations == right.allocations &&
         left.merged_halo_messages == right.merged_halo_messages &&
         left.merged_halo_bytes == right.merged_halo_bytes &&
         left.numeric_refills == right.numeric_refills &&
         left.hierarchy_rebuilds == right.hierarchy_rebuilds &&
         left.cache_publishes == right.cache_publishes &&
         left.linear_iterations == right.linear_iterations &&
         left.stage_wall_nanoseconds == right.stage_wall_nanoseconds;
}

bool test_exact_accumulation_and_peak() {
  ResourceCounters counters{10U, 1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U};
  const ResourceCounters increment{12U, 7U, 8U, 9U, 10U, 11U, 12U, 13U,
                                   14U};

  bool passed = true;
  passed &= expect(static_cast<bool>(add_resource_counters(counters, increment)),
                   "valid resource counters accumulate");
  passed &= expect(equal(counters,
                         ResourceCounters{12U, 8U, 10U, 12U, 14U, 16U, 18U,
                                          20U, 22U}),
                   "resource sums are exact and workspace uses max");

  passed &= expect(
      static_cast<bool>(add_resource_counters(
          counters,
          ResourceCounters{9U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U})),
      "a smaller workspace observation is accepted");
  passed &= expect(counters.peak_workspace_bytes == 12U,
                   "workspace peak never decreases");

  ResourceCounters maximum_peak{};
  passed &= expect(
      static_cast<bool>(add_resource_counters(
          maximum_peak,
          ResourceCounters{std::numeric_limits<std::uint64_t>::max(), 0U, 0U,
                           0U, 0U, 0U, 0U, 0U, 0U})),
      "UINT64_MAX is a valid peak observation because peaks are not summed");
  passed &= expect(maximum_peak.peak_workspace_bytes ==
                       std::numeric_limits<std::uint64_t>::max(),
                   "maximum workspace peak is preserved exactly");
  return passed;
}

bool test_contract_boundaries() {
  const ResourceContract contract{10U, 11U, 12U, 13U, 14U, 15U, 16U, 17U,
                                  18U};
  const ResourceCounters exact{10U, 11U, 12U, 13U, 14U, 15U, 16U, 17U, 18U};

  bool passed = true;
  passed &= expect(static_cast<bool>(validate_resource_counters(contract, exact)),
                   "every resource may equal its contract ceiling");
  passed &= expect(
      static_cast<bool>(validate_resource_counters(ResourceContract{},
                                                   ResourceCounters{})),
      "the default zero-allocation contract accepts zero use");

  ResourceCounters one_allocation{};
  one_allocation.allocations = 1U;
  passed &= expect(
      validate_resource_counters(ResourceContract{}, one_allocation).code ==
          StatusCode::invalid_plan,
      "allocation allowance is enforced by the contract, including default zero");
  return passed;
}

bool test_each_contract_limit() {
  const ResourceContract contract{20U, 20U, 20U, 20U, 20U, 20U, 20U, 20U,
                                  20U};
  std::array<ResourceCounters, 9U> excess{};
  excess[0].peak_workspace_bytes = 21U;
  excess[1].allocations = 21U;
  excess[2].merged_halo_messages = 21U;
  excess[3].merged_halo_bytes = 21U;
  excess[4].numeric_refills = 21U;
  excess[5].hierarchy_rebuilds = 21U;
  excess[6].cache_publishes = 21U;
  excess[7].linear_iterations = 21U;
  excess[8].stage_wall_nanoseconds = 21U;

  bool passed = true;
  for (std::size_t index = 0U; index < excess.size(); ++index) {
    const Status status = validate_resource_counters(contract, excess[index]);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "each resource counter has an independent ceiling");
    passed &= expect(status.detail != 0U,
                     "a resource ceiling failure identifies its field");
  }
  return passed;
}

bool test_overflow_is_atomic() {
  constexpr std::uint64_t maximum =
      std::numeric_limits<std::uint64_t>::max();
  std::array<ResourceCounters, 8U> bases{};
  std::array<ResourceCounters, 8U> increments{};

  bases[0].allocations = maximum;
  increments[0].allocations = 1U;
  bases[1].merged_halo_messages = maximum;
  increments[1].merged_halo_messages = 1U;
  bases[2].merged_halo_bytes = maximum;
  increments[2].merged_halo_bytes = 1U;
  bases[3].numeric_refills = maximum;
  increments[3].numeric_refills = 1U;
  bases[4].hierarchy_rebuilds = maximum;
  increments[4].hierarchy_rebuilds = 1U;
  bases[5].cache_publishes = maximum;
  increments[5].cache_publishes = 1U;
  bases[6].linear_iterations = maximum;
  increments[6].linear_iterations = 1U;
  bases[7].stage_wall_nanoseconds = maximum;
  increments[7].stage_wall_nanoseconds = 1U;

  bool passed = true;
  for (std::size_t index = 0U; index < bases.size(); ++index) {
    // A changed peak and an earlier valid addition expose partial updates if the
    // implementation mutates the destination before discovering overflow.
    bases[index].peak_workspace_bytes = 7U;
    increments[index].peak_workspace_bytes = 9U;
    if (index != 0U) {
      bases[index].allocations = 3U;
      increments[index].allocations = 4U;
    }
    const ResourceCounters before = bases[index];
    const Status status = add_resource_counters(bases[index], increments[index]);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "UINT64 counter overflow is rejected");
    passed &= expect(status.detail != 0U,
                     "counter overflow carries a diagnostic detail");
    passed &= expect(equal(bases[index], before),
                     "overflow leaves the complete counter set unchanged");
  }
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  passed &= test_exact_accumulation_and_peak();
  passed &= test_contract_boundaries();
  passed &= test_each_contract_limit();
  passed &= test_overflow_is_atomic();
  return passed ? 0 : 1;
}
