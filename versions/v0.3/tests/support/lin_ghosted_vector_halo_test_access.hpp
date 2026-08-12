// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
#error "GhostedVectorHalo test access is available only in tests-on builds"
#endif

#include "hundun/rt_halo_performance_counters.hpp"

namespace hundun::linear::test {

void set_next_halo_performance_counters(
    runtime::HaloPerformanceCounters counters) noexcept;

class GhostedVectorHaloTestAccess final {
 public:
  static void set_initial_performance_counters(
      runtime::HaloPerformanceCounters counters) noexcept {
    set_next_halo_performance_counters(counters);
  }
};

}  // namespace hundun::linear::test
