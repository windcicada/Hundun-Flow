// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace hundun::application::detail {

template <class Clock, class Operation>
double measure_performance_region(Clock&& clock, Operation&& operation) {
  const double start = clock();
  operation();
  return clock() - start;
}

}  // namespace hundun::application::detail
