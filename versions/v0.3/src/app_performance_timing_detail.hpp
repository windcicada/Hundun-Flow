// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

namespace hundun::application::detail {

template <class Clock, class Operation>
double measure_performance_region(Clock&& clock, Operation&& operation) {
  const double start = clock();
  operation();
  return clock() - start;
}

}  // namespace hundun::application::detail
