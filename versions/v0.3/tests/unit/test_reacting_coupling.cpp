// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/flow_reacting_coupling_detail.hpp"

#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <vector>

namespace {

using hundun::flow::OpenReactingStepOperators;
using hundun::flow::ReactingOperatorUpdate;
using hundun::flow::detail::ReactingAttemptState;
using hundun::flow::detail::ReactingLayerState;

ReactingLayerState layer(double first, double second) {
  return {{first, second}, {10.0}, 101325.0};
}

void production_step_contract(const hundun::runtime::MpiContext &mpi) {
  ReactingAttemptState state(1U, 2U, layer(0.8, 0.2), layer(0.75, 0.25));
  const auto history = state.history();
  std::vector<std::string> calls;
  OpenReactingStepOperators operators;
  operators.chemistry = [&](std::uint32_t call, double time, double duration,
                            const ReactingLayerState &working) {
    calls.push_back("C" + std::to_string(call));
    HUNDUN_CHECK(duration == 0.05);
    HUNDUN_CHECK(time == (call == 1U ? 2.0 : 2.05));
    if (call == 2U)
      HUNDUN_CHECK(working.rho_y_kg_per_m3 ==
                   std::vector<double>({0.73, 0.27}));
    return ReactingOperatorUpdate{true, {}, {-0.01, 0.01}, {}};
  };
  operators.scalar_transport =
      [&](double time, double duration, const ReactingLayerState &working) {
        calls.push_back("T");
        HUNDUN_CHECK(time == 2.05);
        HUNDUN_CHECK(duration == 0.1);
        HUNDUN_CHECK(working.rho_y_kg_per_m3 ==
                     std::vector<double>({0.74, 0.26}));
        return ReactingOperatorUpdate{true, {}, {-0.01, 0.01}, {2.0}};
      };
  operators.pressure_corrector =
      [&](std::uint32_t call, double, std::uint64_t epoch,
          const ReactingLayerState &working, const std::vector<double> &species,
          const std::vector<double> &enthalpy, std::string &) {
        calls.push_back("P" + std::to_string(call));
        HUNDUN_CHECK(epoch == call);
        if (call == 2U) {
          HUNDUN_CHECK(working.rho_y_kg_per_m3 ==
                       std::vector<double>({0.72, 0.28}));
          HUNDUN_CHECK(species == std::vector<double>({-0.03, 0.03}));
          HUNDUN_CHECK(enthalpy == std::vector<double>({2.0}));
        }
        return true;
      };
  operators.collective_validation =
      [&](bool ok, const std::string &message) {
        return hundun::runtime::collective_status(mpi, ok, message);
      };

  const auto report = hundun::flow::attempt_open_reacting_step(
      state, 2.0, 0.1, operators);
  HUNDUN_CHECK(report.accepted);
  HUNDUN_CHECK(hundun::flow::validate_reacting_step_report(
      hundun::flow::ReactingCouplingSchedule::second_order(), report));
  HUNDUN_CHECK(calls ==
               std::vector<std::string>({"C1", "T", "P1", "C2", "P2"}));
  HUNDUN_CHECK(state.history() == layer(0.75, 0.25));
  HUNDUN_CHECK(state.history() != history);
  HUNDUN_CHECK(state.committed().rho_y_kg_per_m3 ==
               std::vector<double>({0.72, 0.28}));
  HUNDUN_CHECK(state.committed().rho_h_tc_j_per_m3 ==
               std::vector<double>({12.0}));
  HUNDUN_CHECK(report.p0_before_pa == report.p0_after_pa);
}

void failure_rollback_contract(const hundun::runtime::MpiContext &mpi) {
  for (const std::string failing : {"C1", "T", "P1", "C2", "P2"}) {
    ReactingAttemptState state(1U, 2U, layer(0.8, 0.2), layer(0.75, 0.25));
    const auto history = state.history();
    const auto committed = state.committed();
    OpenReactingStepOperators operators;
    operators.chemistry = [&](std::uint32_t call, double, double,
                              const ReactingLayerState &) {
      const std::string stage = "C" + std::to_string(call);
      const bool injected_rank =
          mpi.size() == 1 || mpi.rank() == 1;
      return ReactingOperatorUpdate{stage != failing || !injected_rank, stage,
                                    {-0.01, 0.01}, {}};
    };
    operators.scalar_transport = [&](double, double,
                                     const ReactingLayerState &) {
      const bool injected_rank =
          mpi.size() == 1 || mpi.rank() == 1;
      return ReactingOperatorUpdate{failing != "T" || !injected_rank, "T",
                                    {}, {1.0}};
    };
    operators.pressure_corrector =
        [&](std::uint32_t call, double, std::uint64_t,
            const ReactingLayerState &, const std::vector<double> &,
            const std::vector<double> &, std::string &message) {
          const std::string stage = "P" + std::to_string(call);
          message = stage;
          const bool injected_rank =
              mpi.size() == 1 || mpi.rank() == 1;
          return stage != failing || !injected_rank;
        };
    operators.collective_validation =
        [&](bool ok, const std::string &message) {
          return hundun::runtime::collective_status(mpi, ok, message);
        };
    const auto report = hundun::flow::attempt_open_reacting_step(
        state, 0.0, 0.1, operators);
    HUNDUN_CHECK(!report.accepted);
    HUNDUN_CHECK(state.history() == history);
    HUNDUN_CHECK(state.committed() == committed);
    HUNDUN_CHECK(!state.attempt_active());
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    using hundun::flow::ReactingCouplingEvent;
    using hundun::flow::ReactingCouplingEventKind;
    const auto schedule = hundun::flow::ReactingCouplingSchedule::second_order();
    HUNDUN_CHECK(schedule.chemistry_call_count == 2U);
    HUNDUN_CHECK(schedule.scalar_transport_count == 1U);
    HUNDUN_CHECK(schedule.pressure_corrector_count == 2U);
    HUNDUN_CHECK(schedule.chemistry_duration_fractions[0] == 0.5);
    HUNDUN_CHECK(schedule.chemistry_duration_fractions[1] == 0.5);
    HUNDUN_CHECK(schedule.piso2_consumes_post_chemistry2);
    HUNDUN_CHECK(schedule.requires_predictor_to_final_delta_flux);

    const std::vector<ReactingCouplingEvent> trace{
        {ReactingCouplingEventKind::chemistry, 2.0, 0.05, 1U},
        {ReactingCouplingEventKind::scalar_transport, 2.05, 0.1, 1U},
        {ReactingCouplingEventKind::piso, 2.05, 0.1, 1U},
        {ReactingCouplingEventKind::chemistry, 2.05, 0.05, 2U},
        {ReactingCouplingEventKind::piso, 2.1, 0.1, 2U}};
    HUNDUN_CHECK(hundun::flow::validate_reacting_coupling_trace(
        schedule, trace, 2.0, 0.1));

    auto wrong = trace;
    wrong[3].duration_s = 0.1;
    HUNDUN_CHECK(!hundun::flow::validate_reacting_coupling_trace(
        schedule, wrong, 2.0, 0.1));
    wrong = trace;
    wrong[4].state_epoch = 1U;
    HUNDUN_CHECK(!hundun::flow::validate_reacting_coupling_trace(
        schedule, wrong, 2.0, 0.1));
    wrong = trace;
    wrong.erase(wrong.begin() + 3);
    HUNDUN_CHECK(!hundun::flow::validate_reacting_coupling_trace(
        schedule, wrong, 2.0, 0.1));
    production_step_contract(mpi);
    failure_rollback_contract(mpi);
  });
}
