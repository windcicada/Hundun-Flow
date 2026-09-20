// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>

using namespace hundun::v04;

bool run(unsigned face, double flux, bool perry, bool backward_euler = false,
         bool immersed = false) {
  auto model = test::product_model(immersed ? Int3{8, 16, 16} : Int3{8, 8, 8});
  if (immersed) {
    // The [-1,1]^3 solid intersects x_min. Its blocked area is exactly 4
    // inside the 4-by-4 physical face; the remaining fluid area is 12.
    model.mesh.lower = {0, -2, -2};
    model.mesh.upper = {2, 2, 2};
    model.mesh.exact_cells = {8, 16, 16};
    model.mesh.minimum_spacing = {.25, .25, .25};
    model.immersed_boundary = ImmersedBoundarySpec{
        "cube_binary.stl", ImmersedFluidSide::outside};
  }
  model.time.scheme = backward_euler ? TimeScheme::backward_euler : TimeScheme::cn_be;
  model.solver.coupling = backward_euler ? CouplingKind::piso : CouplingKind::outer_corrected;
  model.time.initial_dt = 1e-3;
  model.solver.terminal.continuity = 1e-13;
  model.solver.terminal.eos = 1e-13;
  model.solver.terminal.gauge = 1e-13;
  model.solver.pressure.relative_tolerance = 1e-14;
  model.solver.pressure.absolute_tolerance = 1e-14;
  if (!backward_euler)
    model.solver.cold_stopping = ColdStoppingSpec{1.0, 1e-4, 1e-11, 1e-4};
  model.legacy_time_fingerprint = model.fingerprint + 1;
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (auto& boundary : model.boundaries) {
    boundary.flow_kind = BoundaryKind::pressure_outlet;
    boundary.pressure = 100000;
    boundary.backflow_temperature = 300;
    boundary.allow_backflow = true;
  }
  model.boundaries[face].flow_kind = BoundaryKind::no_slip_wall;
  model.boundaries[face].thermal_kind = BoundaryKind::heat_flux_wall;
  model.boundaries[face].heat_flux = flux;
  if (perry) {
    auto& gas = model.thermophysics.species.front();
    gas.transport_law = TransportLaw::perry;
    gas.viscosity_reference = gas.conductivity = 0;
    gas.prandtl = 0.70;
    gas.critical_temperature = 126.2;
    gas.critical_pressure = 33.5;
  }
  CompiledCasePlan plan;
  auto status = ProductCompiler::compile(MPI_COMM_WORLD, model,
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data", plan);
  const bool schedule_matches = status && plan.summary().coupling == model.solver.coupling &&
      plan.summary().time_scheme == model.time.scheme;
  ProductDriver driver;
  if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 100000;
  initial.temperature = 300;
  if (status) status = driver.initialize(initial);
  DriverStepReport report;
  bool quiescent_scale = backward_euler;
  for (int step = 0; step < 2 && status; ++step) {
    status = driver.advance({model.time.initial_dt / model.time.convective_cfl,
                            1, 1, 1, 1}, report);
    if (step == 0 && status && !backward_euler)
      quiescent_scale = report.piso.cold.momentum_reference_from_corrector &&
          report.piso.cold.momentum_reference_scale > 1e-12;
  }
  const double area = immersed ? 12.0 : (face < 2 ? 0.5 : (face < 4 ? 1.0 : 2.0));
  RestartSnapshot snapshot;
  if (status) status = driver.committed_restart_snapshot(snapshot);
  CommittedOutputSnapshot output;
  if (status && immersed) status = driver.committed_output_snapshot(output);
  double enthalpy_change{};
  double solid_enthalpy_drift{};
  int solid_cells{};
  bool found_enthalpy = false;
  if (status) {
    const double initial_h = 3.5 * 8314.46261815324 / 28.96546 * 300;
    for (std::size_t i = 0; i < snapshot.fields.size; ++i) {
      const auto& field = snapshot.fields.data[i];
      if (field.role != RestartFieldRole::enthalpy) continue;
      found_enthalpy = true;
      const auto cells = field.values.interior;
      for (int z = 0; z < cells.z; ++z)
        for (int y = 0; y < cells.y; ++y)
          for (int x = 0; x < cells.x; ++x) {
            const double dh = field.values.unchecked({x, y, z}, 0) - initial_h;
            const auto flat = static_cast<std::size_t>(x) +
                static_cast<std::size_t>(cells.x) * (y + cells.y * z);
            if (immersed && output.cell_activity.data[flat] == 0) {
              ++solid_cells;
              solid_enthalpy_drift = std::max(solid_enthalpy_drift, std::abs(dh));
            } else {
              enthalpy_change += dh;
            }
          }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &enthalpy_change, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &solid_enthalpy_drift, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &solid_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  bool passed = status && report.accepted && report.accepted_step == 2 &&
      found_enthalpy && flux * enthalpy_change < 0 && quiescent_scale && schedule_matches &&
      (!immersed || (solid_cells == 256 && solid_enthalpy_drift == 0)) &&
      report.conservation.valid &&
      // The 300 K inlet reservoir also conducts heat during cooling. The
      // exact prescribed wall contribution is verified by the independent
      // enthalpy operator test; this bound includes the open-face exchange.
      std::abs(report.conservation.conductive_heat_input + flux * area) <
          1e-4 * std::abs(flux * area) &&
      std::abs(report.conservation.total_energy_balance_defect) <= 1e-6 * std::abs(flux * area) &&
      report.piso.eos_residual < 1e-12 && report.piso.continuity_residual < 1e-10;
  std::cerr << "heat_wall face=" << face << " q=" << flux << " perry=" << perry
            << " BE=" << backward_euler
            << " IBM=" << immersed << " solid_cells=" << solid_cells
            << " solid_h_drift=" << solid_enthalpy_drift
            << " status=" << unsigned(status.code) << '/' << status.detail
            << " heat=" << report.conservation.conductive_heat_input
            << " dh_sum=" << enthalpy_change
            << " energy_defect=" << report.conservation.total_energy_balance_defect
            << " expected=" << -flux * area << " passed=" << passed << '\n';
  if (!status && backward_euler) {
    const auto& g = report.pressure_energy_globalization;
    std::cerr << std::setprecision(17) << "heat_globalization valid=" << g.valid
              << " corrector=" << unsigned(g.corrector) << " samples=" << unsigned(g.sample_count)
              << " baseline_C=" << g.baseline.global_normalized_continuity
              << " baseline_E=" << g.baseline.global_normalized_energy
              << " weight=" << g.baseline.energy_merit_weight
              << " dp=" << g.maximum_absolute_pressure_correction
              << " dh=" << g.maximum_absolute_enthalpy_correction << '\n';
    for (unsigned i = 0; i < g.sample_count; ++i) {
      const auto& c = g.candidates[i];
      std::cerr << "heat_candidate alpha=" << c.alpha
                << " status=" << unsigned(g.candidate_evaluation_status[i].code)
                << '/' << g.candidate_evaluation_status[i].detail
                << " C=" << c.global_normalized_continuity
                << " E=" << c.global_normalized_energy
                << " admissible=" << c.thermodynamically_admissible
                << " finite=" << c.state_and_flux_finite << '\n';
    }
  }
  int all = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return all;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  for (bool backward_euler : {false, true})
    for (unsigned face = 0; face < 6; ++face)
      for (bool perry : {false, true})
        passed = run(face, face % 2 ? -10.0 : 10.0, perry, backward_euler) && passed;
  for (bool backward_euler : {false, true})
    for (double flux : {-10.0, 10.0})
      passed = run(0, flux, false, backward_euler, true) && passed;
  MPI_Finalize();
  return passed ? 0 : 1;
}
