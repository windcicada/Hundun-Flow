// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string_view>
using namespace hundun::v04;

bool run(int axis, ImmersedFluidSide side, bool mass_inlet, bool reverse,
         bool immersed, int rank, double initial_temperature = 300) {
  auto model = test::product_model({16, 16, 16});
  model.mesh.lower = {-2, -2, -2};
  model.mesh.upper = {2, 2, 2};
  model.mesh.minimum_spacing = {.25, .25, .25};
  const auto component = [axis](auto& value) -> auto& {
    return axis == 0 ? value.x : axis == 1 ? value.y : value.z;
  };
  component(model.mesh.lower) = -.5;
  component(model.mesh.upper) = .5;
  component(model.mesh.exact_cells) = 8;
  component(model.mesh.minimum_spacing) = .125;
  if (immersed)
    model.immersed_boundary = ImmersedBoundarySpec{"cube_binary.stl", side};
  model.fingerprint += 32U * axis + 8U * unsigned(side) + 4U * mass_inlet +
                       2U * reverse + immersed;
  model.legacy_time_fingerprint = model.fingerprint + 100U;
  model.time.scheme = TimeScheme::cn_be;
  model.time.initial_dt = model.time.maximum_dt = 1e-3;
  model.time.maximum_growth = 1.0;
  model.time.convective_cfl = .30;
  model.time.convective_cfl_margin = .05;
  model.solver.coupling = CouplingKind::outer_corrected;
  model.solver.cold_stopping = ColdStoppingSpec{1.0, 1e-8, 1e-11, 1e-8};
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (auto& boundary : model.boundaries)
    boundary.flow_kind = BoundaryKind::symmetry;
  const int inlet_end = reverse ? 1 : 0;
  const double direction = reverse ? -1. : 1.;
  auto& inlet = model.boundaries[2 * axis + inlet_end];
  inlet.flow_kind = mass_inlet ? BoundaryKind::mass_flow_inlet : BoundaryKind::velocity_inlet;
  inlet.temperature = 300;
  if (mass_inlet) {
    inlet.mass_flow_rate = .4;
    component(inlet.direction) = direction;
  } else {
    component(inlet.velocity) = .1 * direction;
  }
  auto& outlet = model.boundaries[2 * axis + 1 - inlet_end];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.pressure = 100000;
  outlet.allow_backflow = true;
  outlet.backflow_temperature = 300;
  CompiledCasePlan plan;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model,
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests/data", plan);
  const char* phase = "compile";
  ProductDriver driver;
  if (status) { phase = "create"; status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver); }
  DriverInitialState initial;
  initial.pressure_reference = 100000;
  initial.temperature = initial_temperature;
  if (status) { phase = "initialize"; status = driver.initialize(initial); }
  DriverStepReport report;
  for (int step = 0; status && step < 2; ++step) {
    phase = "advance";
    status = driver.advance({1e-3 / .30, 1, 1, 1, 1}, report);
  }
  RestartSnapshot restart;
  CommittedOutputSnapshot output;
  if (status) { phase = "snapshot"; status = driver.committed_restart_snapshot(restart); }
  if (status) status = driver.committed_output_snapshot(output);
  const RestartFieldView* velocity = nullptr;
  const RestartFieldView* pressure = nullptr;
  for (std::size_t i = 0; status && i < restart.fields.size; ++i) {
    const auto& field = restart.fields.data[i];
    if (field.role == RestartFieldRole::velocity) velocity = &field;
    if (field.role == RestartFieldRole::pressure_perturbation) pressure = &field;
  }
  bool valid = status && report.piso.cold.midpoint_enthalpy && velocity && pressure && output.committed &&
      restart.final_mass_flux.certificate.valid();
  std::array<double, 6U> sums{}; // inlet/outlet flow, opening areas, fluid/solid cells
  std::array<double, 2U> errors{}; // inactive-face flux and inlet EOS flux error
  double solid_speed = 0;
  if (valid) {
    const auto n = output.patch.cells;
    const auto flux = axis == 0 ? restart.final_mass_flux.x :
                      axis == 1 ? restart.final_mass_flux.y : restart.final_mass_flux.z;
    constexpr double area = .25 * .25;
    constexpr double gas_constant = 8314.46261815324 / 28.96546;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x) {
          const Int3 c{x, y, z};
          const auto index = std::size_t(x) + std::size_t(n.x) * (y + std::size_t(n.y) * z);
          const bool fluid = !immersed || output.cell_activity.data[index] == 1U;
          sums[fluid ? 4 : 5] += 1;
          for (std::size_t a = 0; a < 3; ++a) {
            const double u = velocity->values.unchecked(c, a);
            valid &= std::isfinite(u);
            if (!fluid) solid_speed = std::max(solid_speed, std::abs(u));
          }
          const int global = component(output.patch.begin) + component(c);
          for (int end = 0; end < 2; ++end) {
            if (global != (end == 0 ? 0 : component(model.mesh.exact_cells) - 1)) continue;
            auto face = c;
            if (end) ++component(face);
            const double actual = flux.unchecked(face);
            valid &= std::isfinite(actual);
            const int slot = end == inlet_end ? 0 : 1;
            sums[slot] += direction * actual;
            if (fluid) {
              sums[2 + slot] += area;
              if (slot == 0 && !mass_inlet) {
                const double p = restart.pressure_reference + pressure->values.unchecked(c, 0);
                const double expected = p / (gas_constant * inlet.temperature) * .1 * direction * area;
                errors[1] = std::max(errors[1], std::abs(actual - expected));
              }
            } else errors[0] = std::max(errors[0], std::abs(actual));
          }
        }
  }
  MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, errors.data(), errors.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &solid_speed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  const double expected_area = !immersed ? 16.0 :
      side == ImmersedFluidSide::inside ? 4.0 : 12.0;
  valid &= report.accepted && report.accepted_step == 2U && report.conservation.valid &&
      sums[2] == expected_area && sums[3] == expected_area &&
      sums[4] == expected_area * 128 && sums[4] + sums[5] == 2048 &&
      sums[0] > 0 && sums[1] > 0 && errors[0] == 0 && solid_speed == 0 &&
      (mass_inlet ? std::abs(sums[0] - .4) < 1e-12 : errors[1] < 1e-12) &&
      report.piso.eos_residual < 1e-10 && report.piso.continuity_residual < 1e-10 &&
      std::abs(report.conservation.mass_balance_defect) < 1e-6 * sums[0] &&
      std::abs(report.conservation.total_energy_balance_defect) < 1e-6 * sums[0] * 3.5 * (8314.46261815324 / 28.96546) * 300;
  int passed = valid ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &passed, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << std::setprecision(17) << "ibm-open axis=" << axis
              << " immersed=" << immersed << " inside=" << (side == ImmersedFluidSide::inside)
              << " mass=" << mass_inlet << " reverse=" << reverse << " phase=" << phase
              << " status=" << unsigned(status.code) << '/' << status.detail
              << " area=" << sums[2] << '/' << sums[3] << " flow=" << sums[0] << '/' << sums[1]
              << " inactive/EOS-flux-error=" << errors[0] << '/' << errors[1]
              << " solid-speed=" << solid_speed << " mass/energy=" << report.conservation.mass_balance_defect
              << '/' << report.conservation.total_energy_balance_defect
              << " EOS/continuity=" << report.piso.eos_residual << '/' << report.piso.continuity_residual
              << " pass=" << passed << '\n';
  return passed != 0;
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool passed = true;
  if (argc == 2 && std::string_view(argv[1]) == "--heat") {
    for (int axis = 0; axis < 3; ++axis)
      for (auto side : {ImmersedFluidSide::inside, ImmersedFluidSide::outside})
        passed &= run(axis, side, true, false, true, rank, 299.9);
    MPI_Finalize();
    return passed ? 0 : 1;
  }
  for (int axis = 0; axis < 3; ++axis)
    for (auto side : {ImmersedFluidSide::inside, ImmersedFluidSide::outside})
      for (bool mass : {false, true})
        for (bool reverse : {false, true})
          passed &= run(axis, side, mass, reverse, true, rank);
  for (bool reverse : {false, true})
    passed &= run(0, ImmersedFluidSide::inside, false, reverse, false, rank);
  MPI_Finalize();
  return passed ? 0 : 1;
}
