// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/rt_halo_performance_counters.hpp"

#include <type_traits>

static_assert(std::is_trivially_copyable_v<
              hundun::runtime::HaloPerformanceCounters>);
static_assert(std::is_same_v<
              decltype(hundun::runtime::HaloPerformanceCounters::
                           completed_wait_seconds),
              double>);

int task25_correctness_header_contract();

int main() {
  hundun::runtime::HaloPerformanceCounters counters;
  return counters.completed_exchanges == 0U
             ? task25_correctness_header_contract()
             : 1;
}
