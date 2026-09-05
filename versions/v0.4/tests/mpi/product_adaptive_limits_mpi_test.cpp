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
}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed =
      uniform_flow(1.0, 1.0, false) && uniform_flow(0.5, 1.0, false) &&
      uniform_flow(0.0, 1.0, false) && uniform_flow(1.0, 0.002, false) &&
      uniform_flow(1.0, 1.0, true);
  int result = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return result != 0 ? 0 : 1;
}
