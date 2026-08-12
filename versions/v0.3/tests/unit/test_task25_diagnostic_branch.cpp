// SPDX-License-Identifier: Apache-2.0

#include "src/app_performance_diagnostic_detail.hpp"
#include "hundun/exec_execution.hpp"
#include "tests/support/test_main.hpp"

#include <cstdint>

int main() {
  return hundun::test::run([] {
    std::uint64_t traversals = 0U;
    const auto before = hundun::execution::allocation_counters();
    const auto disabled = hundun::application::detail::
        execute_flow_diagnostic_branch(
            false, 2U, 2U, [&] {
              ++traversals;
              return std::uint64_t{17U};
            });
    HUNDUN_CHECK(!disabled.submitted);
    HUNDUN_CHECK(disabled.logical_bytes == 0U);
    HUNDUN_CHECK(traversals == 0U);
    HUNDUN_CHECK(hundun::execution::allocation_counters().allocation_events ==
                 before.allocation_events);

    const auto not_due = hundun::application::detail::
        execute_flow_diagnostic_branch(
            true, 2U, 1U, [&] {
              ++traversals;
              return std::uint64_t{17U};
            });
    HUNDUN_CHECK(!not_due.submitted);
    HUNDUN_CHECK(not_due.logical_bytes == 0U);
    HUNDUN_CHECK(traversals == 0U);

    const auto due = hundun::application::detail::
        execute_flow_diagnostic_branch(
            true, 2U, 2U, [&] {
              ++traversals;
              return std::uint64_t{17U};
            });
    HUNDUN_CHECK(due.submitted);
    HUNDUN_CHECK(due.logical_bytes == 17U);
    HUNDUN_CHECK(traversals == 1U);
  });
}
