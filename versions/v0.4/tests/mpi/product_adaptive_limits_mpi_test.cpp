// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"

namespace {
using namespace hundun::v04;

bool uniform_flow(double speed, double external_scale, bool fixed) {
  ValidatedModel model = test::product_model({8, 8, 8});
  model.turbulence = TurbulenceKind::none;
  model.time.control =
      fixed ? TimeControlKind::fixed : TimeControlKind::adaptive_flow;
  model.time.initial_dt = model.time.maximum_dt = 0.1;
  model.time.maximum_retries = 1U;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 101325.0;
  initial.temperature = 300.0;
  initial.velocity = {10.0 * speed, -4.0 * speed, 2.0 * speed};
  if (status) status = driver.initialize(initial);
  const double rate = speed * (10.0 / 0.25 + 4.0 / 0.125 + 2.0 / 0.0625);
  const double flow_dt =
      rate == 0.0 ? model.time.maximum_dt : model.time.convective_cfl / rate;
  const double expected =
      std::min(flow_dt, std::min(model.time.maximum_dt,
                                 model.time.convective_cfl * external_scale));
  DriverStepReport step;
  LocalTimeLimits limits{external_scale, 1.0, 1.0, 1.0, 1.0};
  if (status) status = driver.constrain_convective_time_limit(limits);
  if (status) status = driver.advance(limits, step);
  bool passed =
      fixed ? !status && step.proposal.dt == model.time.initial_dt
            : status && step.accepted && step.attempts == 1U &&
                  std::abs(step.proposal.dt - expected) <= 1.0e-13 * expected &&
                  step.piso.committed_convective_cfl_out_max <=
                      model.time.convective_cfl * (1.0 + 1.0e-12);
  if (!passed)
    std::cerr << "FAIL accepted-flow dt speed=" << speed << " fixed=" << fixed
              << " status=" << static_cast<unsigned>(status.code) << '/'
              << status.detail << " attempts=" << step.attempts
              << " dt=" << step.proposal.dt << " expected=" << expected << '\n';
  return passed;
}
bool minimum_dt_rejection_preserves_state_and_reports_failure() {
  ValidatedModel model = test::product_model({8, 8, 8});
  model.turbulence = TurbulenceKind::none;
  model.time.control = TimeControlKind::adaptive_flow;
  model.time.initial_dt = model.time.maximum_dt = 0.01;
  model.time.minimum_dt = 0.001;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  if (status) status = driver.initialize({});
  DriverStepReport report;
  if (status) status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, report);
  int okay = status && report.accepted && report.accepted_step == 1U ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!okay) return false;
  const double committed_time = report.accepted_time;
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  // Only one rank constrains the proposal below minimum_dt. Reuse the last
  // report, as a caller naturally does in a time loop.
  const double limit = rank == ranks - 1 ? 1.0e-6 : 1.0;
  status = driver.advance({limit, 1.0, 1.0, 1.0, 1.0}, report);
  bool passed = status.code == StatusCode::rejected_step &&
                status.detail == 454U && report.failure.code == status.code &&
                report.failure.detail == status.detail && !report.accepted &&
                report.attempts == 0U && report.failed_stage == 0U &&
                report.accepted_step == 1U &&
                report.accepted_time == committed_time;
  passed &= report.completion.outcome.detail == status.detail &&
            report.completion.stop_reason.detail == 454U &&
            report.completion.first_failure.attempt == 0U &&
            report.completion.last_failure.attempt == 0U;
  passed &=
      report.initial_time_proposal.evaluated &&
      std::abs(report.initial_time_proposal.proposed_dt - 8.0e-7) < 1.0e-20 &&
      report.initial_time_proposal.minimum_dt == 0.001 &&
      report.proposal.generation == 0U;
  if (!passed)
    std::cerr << "FAIL pre-solve time rejection: status="
              << static_cast<int>(status.code) << '/' << status.detail
              << " report_failure=" << static_cast<int>(report.failure.code)
              << '/' << report.failure.detail << " accepted=" << report.accepted
              << " attempts=" << report.attempts << '\n';
  CommittedOutputSnapshot snapshot;
  const Status snapshot_status = driver.committed_output_snapshot(snapshot);
  passed &=
      snapshot_status && snapshot.step == 1U && snapshot.time == committed_time;
  // An invalid local input has no trustworthy globally evaluated proposal.
  // It must clear the previous rejection diagnostic instead of reusing it.
  status = driver.advance({rank == ranks - 1 ? -1.0 : 1.0, 1.0, 1.0, 1.0, 1.0},
                          report);
  passed &= !status && !report.initial_time_proposal.evaluated &&
            !report.accepted && report.attempts == 0U &&
            report.failure.code == status.code &&
            report.failure.detail == status.detail;
  DriverStepReport recovered;
  status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, recovered);
  passed &= status && recovered.accepted && recovered.accepted_step == 2U &&
            recovered.attempts == 1U &&
            recovered.accepted_time > committed_time;
  passed &=
      recovered.initial_time_proposal.evaluated &&
      recovered.initial_time_proposal.proposed_dt == recovered.proposal.dt;
  return passed;
}
bool uninitialized_driver_resets_report() {
  ProductDriver driver;
  DriverStepReport report;
  report.accepted = true;
  report.accepted_step = 123U;
  report.attempts = 5U;
  report.failed_stage = 50U;
  const Status status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, report);
  const bool passed = !status && !report.accepted &&
                      report.accepted_step == 0U && report.attempts == 0U &&
                      report.failed_stage == 0U &&
                      report.failure.code == status.code &&
                      report.failure.detail == status.detail;
  if (!passed)
    std::cerr << "FAIL uninitialized driver retained previous report\n";
  return passed;
}
}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed =
      uninitialized_driver_resets_report() && uniform_flow(1.0, 1.0, false) &&
      uniform_flow(0.5, 1.0, false) && uniform_flow(0.0, 1.0, false) &&
      uniform_flow(1.0, 0.002, false) && uniform_flow(1.0, 1.0, true) &&
      minimum_dt_rejection_preserves_state_and_reports_failure();
  int result = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return result != 0 ? 0 : 1;
}
