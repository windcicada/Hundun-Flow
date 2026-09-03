// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/flow_adaptive_time_control.hpp"
#include "tests/support/test_main.hpp"

#include <type_traits>

int main() {
  return hundun::test::run([] {
    using namespace hundun::flow;
    static_assert(!std::is_copy_constructible_v<TimeAdvanceReport>);
    static_assert(std::is_nothrow_move_constructible_v<TimeAdvanceReport>);
    static_assert(!std::is_move_assignable_v<TimeAdvanceReport>);
    static_assert(
        !std::is_copy_constructible_v<TimeControlDiagnosticSource>);
    static_assert(
        std::is_nothrow_move_constructible_v<TimeControlDiagnosticSource>);
    static_assert(!std::is_copy_constructible_v<Bdf2RetryController>);
    static_assert(std::is_nothrow_move_constructible_v<Bdf2RetryController>);
    static_assert(std::is_same_v<
                  decltype(std::declval<const Bdf2RetryController &>().state()),
                  TimeControlState>);
  });
}
