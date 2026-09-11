// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include <cmath>
#include <iostream>
#include <mpi.h>

using namespace hundun::v04;

bool run(double speed, double composition = 0.3, bool inlet_reference = false) {
  auto model = test::product_model({12, 8, 8});
  model.time.scheme = TimeScheme::cn_be;
  model.time.initial_dt = 9.7088612375381536e-8;
  if (inlet_reference) model.solver.cold_stopping = ColdStoppingSpec{0.001};
  model.legacy_time_fingerprint = model.fingerprint + 1;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  auto &inlet = model.boundaries[0];
  inlet.flow_kind = BoundaryKind::mass_flow_inlet;
  inlet.direction = {1, 0, 0};
  inlet.mass_flow_rate = 0.01;
  inlet.temperature = 300;
  inlet.scalars = {{"air", ScalarBoundaryKind::dirichlet, 0.3}};
  auto &outlet = model.boundaries[1];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = 100000;
  outlet.allow_backflow = true;
  outlet.backflow_temperature = 310;
  outlet.backflow_velocity = {0, 0.01, 0};
  outlet.scalars = {{"air", ScalarBoundaryKind::zero_gradient}};
  outlet.scalars[0].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[0].backflow_value = 0.1;
  auto &gas = model.thermophysics.species.front();
  gas.transport_law = TransportLaw::coast_perry;
  gas.viscosity_reference = gas.conductivity = 0;
  gas.prandtl = 0.70;
  gas.critical_temperature = 126.2;
  gas.critical_pressure = 33.5;
  auto dependent = gas;
  dependent.stable_name = "balance";
  model.thermophysics.species.push_back(dependent);
  model.transported_scalars = {
      {"air", TransportedScalarRole::species, 0.70, 0.70}};
  CompiledCasePlan plan;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 100100;
  initial.temperature = 300;
  initial.velocity = {speed, 0, 0};
  initial.transported_scalars = {&composition, 1};
  if (status)
    status = driver.initialize(initial);
  DriverStepReport report;
  for (int i = 0; i < 3 && status; ++i)
    status = driver.advance(
        {model.time.initial_dt / model.time.convective_cfl, 1, 1, 1, 1},
        report);
  RestartSnapshot snapshot;
  if (status)
    status = driver.committed_restart_snapshot(snapshot);
  bool passed = status && report.accepted && snapshot.step == 3;
  if (inlet_reference)
    passed &= report.piso.cold.species_reference_scales.size() == 2 &&
        report.piso.cold.species_reference_scales[0] > 1.0 &&
        report.piso.cold.reference_residual[2] < model.solver.cold_stopping->species;
  if (status && snapshot.patch.begin.x + snapshot.patch.cells.x == 12) {
    ConstFieldView u, p, h, y;
    for (std::size_t i = 0; i < snapshot.fields.size; ++i) {
      const auto &f = snapshot.fields.data[i];
      if (f.role == RestartFieldRole::velocity)
        u = f.values;
      if (f.role == RestartFieldRole::pressure_perturbation)
        p = f.values;
      if (f.role == RestartFieldRole::enthalpy)
        h = f.values;
      if (f.role == RestartFieldRole::transported_scalar ||
          f.role == RestartFieldRole::independent_species)
        y = f.values;
    }
    passed &= u.base && p.base && h.base && y.base;
    if (passed) {
      const Int3 c{snapshot.patch.cells.x - 1, 2, 2}, g{c.x + 1, c.y, c.z};
      const auto face = [&](ConstFieldView f, unsigned component = 0) {
        return 0.5 * (f.unchecked(c, component) + f.unchecked(g, component));
      };
      const double un = face(u);
      passed &= std::abs(snapshot.pressure_reference + face(p) - 100000) < 1e-8;
      passed &= std::abs(un - u.unchecked(c, 0)) < 1e-12 && un * speed > 0;
      passed &= un > speed + 1e-5;
      if (speed < 0) {
        passed &= std::abs(face(y) - 0.1) < 1e-12;
        passed &= std::abs(face(u, 1) - 0.01) < 1e-12;
        const double cp = 3.5 * 8314.46261815324 /
                          model.thermophysics.species.front().molecular_weight;
        passed &= std::abs(face(h) / cp - 310) < 1e-8;
        ConstFieldView old_u;
        for (std::size_t i = 0; i < snapshot.previous_fields.size; ++i)
          if (snapshot.previous_fields.data[i].role ==
              RestartFieldRole::velocity)
            old_u = snapshot.previous_fields.data[i].values;
        const double density =
            100000 * model.thermophysics.species.front().molecular_weight /
            (8314.46261815324 * 310);
        const double expected_flux =
            density * 0.5 * (u.unchecked(c, 0) + old_u.unchecked(c, 0)) /
            (8.0 * 16.0);
        const double actual_flux = snapshot.final_mass_flux.x.unchecked(g);
        passed &= std::abs(actual_flux - expected_flux) <
                  1e-6 * std::abs(expected_flux);
        std::cerr << "cold_outlet EOS_flux=" << actual_flux
                  << " expected=" << expected_flux << '\n';
      } else {
        passed &= std::abs(face(y) - y.unchecked(c, 0)) < 1e-12;
        passed &= std::abs(face(h) - h.unchecked(c, 0)) < 1e-8;
      }
      std::cerr << "cold_outlet speed=" << speed
                << " p=" << snapshot.pressure_reference + face(p) << " U=" << un
                << " Y=" << face(y) << " passed=" << passed << '\n';
    }
  }
  if (!status)
    std::cerr << "cold_outlet status=" << unsigned(status.code) << '/'
              << status.detail << '\n';
  int all = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return all;
}

bool run_air() {
  auto model = test::product_model({12, 8, 8});
  model.time.scheme = TimeScheme::cn_be;
  model.time.initial_dt = 1.3808912271980336e-5;
  model.legacy_time_fingerprint = model.fingerprint + 1;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  auto &inlet = model.boundaries[0];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.velocity = {0.1, 0, 0};
  inlet.direction = {1, 0, 0};
  inlet.temperature = 300;
  auto &outlet = model.boundaries[1];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = 100000;
  outlet.backflow_temperature = 300;
  CompiledCasePlan plan;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 100000;
  initial.temperature = 300;
  initial.velocity = {0.1, 0, 0};
  if (status) status = driver.initialize(initial);
  DriverStepReport report;
  for (int i = 0; i < 3 && status; ++i)
    status = driver.advance({model.time.initial_dt / model.time.convective_cfl,
                             1, 1, 1, 1}, report);
  bool passed = status && report.accepted && report.accepted_step == 3 &&
      report.piso.cold.active && report.piso.cold.species_solve_calls == 0 &&
      report.piso.cold.species_iterations == 0 &&
      report.piso.eos_residual < 1e-12 && report.piso.continuity_residual < 1e-10;
  if (!passed)
    std::cerr << "cold_air status=" << unsigned(status.code) << '/'
              << status.detail << " step=" << report.accepted_step << '\n';
  int all = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return all;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  const bool out = run(0.02), back = run(-0.02), air = run_air();
  const bool zero_species = run(0.02, 0.0, true);
  MPI_Finalize();
  return out && back && air && zero_species ? 0 : 1;
}
