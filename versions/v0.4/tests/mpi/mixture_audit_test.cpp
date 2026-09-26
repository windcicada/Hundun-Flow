// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "hundun/v04_app.hpp"
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace hundun::v04;

// Exercise the final conservative audits in a hot mixture at a small step.
// Pre-correction temperature guesses are not the final equation residuals.
bool run(const char* fixture, int rank) {
  constexpr double dt = 1e-9;
  constexpr double temperature = 1200;
  ValidatedModel model;
  auto status = CaseCompiler::load_and_compile(MPI_COMM_WORLD, fixture, model);
  model.solver.cold_stopping = ColdStoppingSpec{1e-6, 1e-8, 1e-12, 1e-12};
  model.time.control = TimeControlKind::fixed;
  model.time.initial_dt = model.time.minimum_dt = model.time.maximum_dt = dt;
  CompiledCasePlan plan;
  if (status) status = ProductCompiler::compile(MPI_COMM_WORLD, model, fixture, plan);
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  std::vector<double> composition(model.transported_scalars.size());
  for (std::size_t i = 0; i < composition.size(); ++i) {
    const auto& name = model.transported_scalars[i].stable_name;
    if (name == "H2") composition[i] = 0.028;
    if (name == "O2") composition[i] = 0.225;
    if (name == "H" || name == "OH") composition[i] = 1e-5;
  }
  DriverInitialState initial;
  initial.pressure_reference = 101325;
  initial.temperature = temperature;
  initial.velocity = {2.0, 0.0, 0.0};
  initial.transported_scalars = {composition.data(), composition.size()};
  if (status) status = driver.initialize(initial);
  DriverStepReport report;

  if (status) status = driver.advance({dt, dt, dt, dt, dt}, report);
  const auto& cold = report.piso.cold;
  bool okay = status && report.accepted && cold.active && cold.stopping &&
              cold.solid_velocity_max == 0;
  const double limits[]{1e-8, 1e-12, 1e-12};
  for (unsigned i = 0; i < 3; ++i) {
    const double residual = cold.reference_residual[i];
    okay = okay && std::isfinite(residual) && residual >= 0 && residual < limits[i];
  }
  if (rank == 0)
    std::cout << "mixture_audit status=" << unsigned(status.code) << '/'
              << status.detail << " accepted=" << report.accepted
              << " outer=" << cold.outer_iterations
              << " M=" << cold.reference_residual[0]
              << " E=" << cold.reference_residual[1]
              << " Y=" << cold.reference_residual[2] << '\n';
  return okay;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank{};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int okay = argc == 2 && run(argv[1], rank);
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return okay ? 0 : 1;
}
