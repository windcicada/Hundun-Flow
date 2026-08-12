// SPDX-License-Identifier: Apache-2.0

#include "flow_reacting_pressure_constraint_detail.hpp"

#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    using namespace hundun::flow::detail;
    ReactingPressureState state(100000.0, 101000.0);
    state.begin_attempt();
    const auto counters_before = mpi.fp64_reduction_counters();
    const double ranks = static_cast<double>(mpi.size());
    const auto closed = attempt_reacting_pressure_constraint(
        state, mpi,
        {ReactingPressureDomain::closed, 0.1, 2.0 / ranks,
         20.0 / ranks, 24.0 / ranks, 0.0});
    const auto counters_after = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(counters_after.collective_calls ==
                 counters_before.collective_calls + 1U);
    HUNDUN_CHECK(closed.global_reduction_count == 1U);
    HUNDUN_CHECK(closed.predictor_p0_pa == 101010.0);
    HUNDUN_CHECK(closed.corrected_p0_pa == 101011.0);
    HUNDUN_CHECK(closed.midpoint_dp0_dt_pa_per_s == 110.0);
    HUNDUN_CHECK(state.trial_p0_pa() == 101011.0);
    HUNDUN_CHECK(state.commit({true, -1, {}}));
    HUNDUN_CHECK(state.history_p0_pa() == 101000.0);
    HUNDUN_CHECK(state.committed_p0_pa() == 101011.0);

    state.begin_attempt();
    const auto partial = attempt_reacting_pressure_constraint(
        state, mpi,
        {ReactingPressureDomain::partially_closed, 0.1, 2.0 / ranks,
         20.0 / ranks, 24.0 / ranks, 10.0 / ranks});
    HUNDUN_CHECK(partial.corrected_p0_pa == 101017.0);
    HUNDUN_CHECK(!state.commit({false, 0, "injected"}));
    HUNDUN_CHECK(state.committed_p0_pa() == 101011.0);

    state.begin_attempt();
    const auto open = attempt_reacting_pressure_constraint(
        state, mpi, {ReactingPressureDomain::open, 0.1, 0.0, 0.0, 0.0, 0.0});
    HUNDUN_CHECK(open.global_reduction_count == 0U);
    HUNDUN_CHECK(open.corrected_p0_pa == 101011.0);
    state.rollback();
  });
}
