// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstddef>

namespace hundun::flow::detail {

void momentum_predictor_reset_collective_selection_calls_raw() noexcept;
std::size_t momentum_predictor_collective_selection_calls_raw() noexcept;

} // namespace hundun::flow::detail

namespace hundun::flow::test {

// Private test seam. This header is not part of the installed/public ABI.
class MomentumPredictorTestAccess final {
public:
  static void reset_collective_selection_calls() noexcept {
    detail::momentum_predictor_reset_collective_selection_calls_raw();
  }
  static std::size_t collective_selection_calls() noexcept {
    return detail::momentum_predictor_collective_selection_calls_raw();
  }
};

} // namespace hundun::flow::test
