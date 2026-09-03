// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_case.hpp"

#include "app_case_detail.hpp"

#include <mpi.h>

#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

namespace {

namespace fs = std::filesystem;
using hundun::v04::CaseCompiler;
using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::ValidatedModel;

constexpr std::string_view kThermophysics = R"data(
HUNDUN_THERMOPHYSICS_V1
temperature_bounds 200 3000
temperature_inversion 1e-12 32
closed_mass_newton 1e-12 24 0.25
species_count 2
species O2
molecular_weight 31.9988
temperature_switch 1000
nasa7_low 3.5 0 0 0 0 0 0
nasa7_high 3.5 0 0 0 0 0 0
transport_sutherland 1.919e-5 273.15 139 0.72
end_species
species N2
molecular_weight 28.0134
temperature_switch 1000
nasa7_low 3.5 0 0 0 0 0 0
nasa7_high 3.5 0 0 0 0 0 0
transport_sutherland 1.663e-5 273.15 107 0.72
end_species
end
)data";

int g_open_count = 0;
int g_observed_rank = -1;

void observe_open(int rank) {
  ++g_open_count;
  g_observed_rank = rank;
}

void reset_observer() {
  g_open_count = 0;
  g_observed_rank = -1;
  hundun::v04::detail::set_file_open_observer_for_test(observe_open);
}

bool expect(bool condition, int rank, std::string_view message) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << message << '\n';
  }
  return condition;
}

void write_file(const fs::path& path, std::string_view content) {
  std::ofstream output(path, std::ios::binary);
  output.write(content.data(), static_cast<std::streamsize>(content.size()));
}

bool same_u64(std::uint64_t value, MPI_Comm communicator) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool same_status(Status status, MPI_Comm communicator) {
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  return same_u64(packed, communicator);
}

bool collective_result(bool local, MPI_Comm communicator) {
  const int value = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }

  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  fs::path root;
  if (rank == 0) {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    root = fs::temp_directory_path() /
           ("hundun-v04-case-broadcast-" + std::to_string(::getpid()) +
            "-" + std::to_string(stamp));
    std::error_code error;
    fs::remove_all(root, error);
    fs::create_directories(root);
    write_file(root / "profile.d", "0 1\n1 2\n");
    write_file(root / "thermophysics.d", kThermophysics);
    write_file(root / "shape.stl", "solid shape\nendsolid shape\n");
    write_file(root / "case.json", R"json({
      "schema_version": 1,
      "units": "SI",
      "mesh": {
        "kind": "tensor_stretched",
        "domain": {"lower": [0, 0, 0], "upper": [2, 1, 1]},
        "exact_cells": [16, 10, 8],
        "base_spacing": [0.25, 0.2, 0.2],
        "minimum_spacing": [0.05, 0.05, 0.05],
        "max_growth_ratio": 1.2,
        "focus_regions": [
          {"lower": [-1, 0.2, 0.2], "upper": [0.5, 0.8, 0.8],
           "target_spacing": [0.1, 0.1, 0.1]},
          {"lower": [-2, 0.2, 0.2], "upper": [0.5, 0.8, 0.8],
           "target_spacing": [0.1, 0.1, 0.1]}
        ],
        "limits": {
          "max_global_cells": 100000,
          "max_memory_bytes_per_rank": 134217728
        },
        "data_files": ["profile.d"],
        "immersed_boundary": {"stl_file": "shape.stl", "fluid_side": "inside"}
      },
      "flow": {
        "model": "single_phase_low_mach_compressible",
        "pressure_reference": "boundary_absolute",
        "reacting": false
      },
      "solver": {
        "coupling": "SIMPLE",
        "pressure_correctors": 2,
        "pressure_linear": {
          "absolute_tolerance": 1e-8,
          "relative_tolerance": 2e-7,
          "maximum_iterations": 333,
          "true_residual_interval": 7,
          "krylov_restart": 16
        },
        "terminal_tolerances": {
          "eos": 3e-6,
          "continuity": 4e-6,
          "closed_mass": 5e-6,
          "gauge": 6e-6
        }
      },
      "turbulence": {"model": "vreman_wall_function"},
      "thermophysics": {"data_file": "thermophysics.d"},
      "transported_scalars": [
        {"stable_name":"O2","role":"species","molecular_schmidt":0.7,"turbulent_schmidt":0.9},
        {"stable_name":"mixture_fraction","role":"passive_scalar","molecular_schmidt":0.8,"turbulent_schmidt":0.6}
      ],
      "boundaries": {
        "x_min": {"flow_kind":"velocity_inlet","thermal_kind":"none","velocity":[1,0,0],"direction":[1,0,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[{"stable_name":"O2","kind":"dirichlet","value":0.21,"backflow_kind":"zero_gradient","backflow_value":0},{"stable_name":"mixture_fraction","kind":"dirichlet","value":1,"backflow_kind":"zero_gradient","backflow_value":0}]},
        "x_max": {"flow_kind":"pressure_outlet","thermal_kind":"none","velocity":[0,0,0],"direction":[1,0,0],"backflow_velocity":[-1,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":0.8,"mach_limit":0.9,"allow_backflow":true,"scalars":[{"stable_name":"O2","kind":"convective","value":0,"backflow_kind":"dirichlet","backflow_value":0.21},{"stable_name":"mixture_fraction","kind":"convective","value":0,"backflow_kind":"dirichlet","backflow_value":0}]},
        "y_min": {"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,1,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[{"stable_name":"O2","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0},{"stable_name":"mixture_fraction","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0}]},
        "y_max": {"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,1,0],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[{"stable_name":"O2","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0},{"stable_name":"mixture_fraction","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0}]},
        "z_min": {"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,0,1],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[{"stable_name":"O2","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0},{"stable_name":"mixture_fraction","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0}]},
        "z_max": {"flow_kind":"symmetry","thermal_kind":"none","velocity":[0,0,0],"direction":[0,0,1],"backflow_velocity":[0,0,0],"mass_flow_rate":0,"pressure":101325,"temperature":300,"total_pressure":101325,"total_temperature":300,"backflow_temperature":300,"heat_flux":0,"relaxation":1,"mach_limit":0.95,"allow_backflow":false,"scalars":[{"stable_name":"O2","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0},{"stable_name":"mixture_fraction","kind":"zero_gradient","value":0,"backflow_kind":"zero_gradient","backflow_value":0}]}
      },
      "schemes": {"momentum":"limited_central2","enthalpy":"central2","species":"tvd2","passive_scalar":"tvd2","diffusion":"central2","limiter":0.75},
      "time": {"control":"adaptive_flow","scheme":"variable_bdf2","initial_dt":0.001,"minimum_dt":1e-8,"maximum_dt":0.1,"convective_cfl":0.8,"viscous_cfl":0.5,"thermal_cfl":0.5,"species_cfl":0.4,"acoustic_cfl":0.9,"maximum_growth":1.2,"retry_factor":0.5,"maximum_retries":6,"minimum_bdf_ratio":0.25,"maximum_bdf_ratio":4.0}
    })json");
  }
  MPI_Barrier(MPI_COMM_WORLD);

  const fs::path nonroot_poison =
      "/nonroot-must-not-open-or-canonicalize-this-case";
  reset_observer();
  ValidatedModel first;
  const Status first_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison, first);

  bool passed = true;
  passed &= expect(static_cast<bool>(first_status), rank,
                   "valid root case compiles collectively");
  passed &= expect(same_status(first_status, MPI_COMM_WORLD), rank,
                   "success status is identical");
  passed &= expect(same_u64(first.fingerprint, MPI_COMM_WORLD), rank,
                   "fingerprint is identical");
  const std::uint64_t model_signature =
      (static_cast<std::uint64_t>(first.mesh.kind) << 24U) |
      (static_cast<std::uint64_t>(first.turbulence) << 16U) |
      (static_cast<std::uint64_t>(first.time.control) << 8U) |
      (static_cast<std::uint64_t>(first.pressure_reference) << 48U) |
      static_cast<std::uint64_t>(first.data_files.size()) |
      (first.immersed_boundary.has_value() ? (1ULL << 40U) : 0U);
  passed &= expect(same_u64(model_signature, MPI_COMM_WORLD), rank,
                   "typed model is identical");
  passed &= expect(first.mesh.kind == hundun::v04::GeometryKind::tensor_stretched &&
                       first.solver.coupling ==
                           hundun::v04::CouplingKind::simple &&
                       first.mesh.lower.x == 0.0 && first.mesh.upper.x == 2.0 &&
                       first.mesh.has_exact_cells &&
                       first.mesh.exact_cells.x == 16 &&
                       first.mesh.exact_cells.y == 10 &&
                       first.mesh.exact_cells.z == 8 &&
                       first.mesh.has_base_spacing &&
                       first.mesh.base_spacing.x == 0.25 &&
                       first.mesh.minimum_spacing.x == 0.05 &&
                       first.mesh.max_growth_ratio == 1.2 &&
                       first.mesh.limits.max_global_cells == 100000 &&
                       first.mesh.limits.max_memory_bytes_per_rank == 134217728,
                   rank, "wire publishes the complete typed mesh");
  passed &= expect(first.mesh.focus_regions.size() == 1U &&
                       first.mesh.focus_regions[0].lower.x == 0.0 &&
                       first.mesh.focus_regions[0].upper.x == 0.5 &&
                       first.mesh.focus_regions[0].target_spacing.x == 0.1,
                   rank, "wire publishes canonical clipped focus regions");
  passed &= expect(
      first.pressure_reference ==
              hundun::v04::PressureReferenceKind::boundary_absolute &&
          first.transported_scalars.size() == 2U &&
          first.transported_scalars[0].stable_name == "O2" &&
          first.transported_scalars[0].role ==
              hundun::v04::TransportedScalarRole::species &&
          first.transported_scalars[0].molecular_schmidt == 0.7 &&
          first.transported_scalars[0].turbulent_schmidt == 0.9 &&
          first.transported_scalars[1].stable_name == "mixture_fraction" &&
          first.transported_scalars[1].role ==
              hundun::v04::TransportedScalarRole::passive_scalar &&
          first.transported_scalars[1].molecular_schmidt == 0.8 &&
          first.transported_scalars[1].turbulent_schmidt == 0.6 &&
          first.boundaries[0].flow_kind ==
              hundun::v04::BoundaryKind::velocity_inlet &&
          first.boundaries[0].scalars.size() == 2U &&
          first.boundaries[0].scalars[0].stable_name == "O2" &&
          first.boundaries[0].scalars[1].stable_name == "mixture_fraction" &&
          first.boundaries[1].scalars[0].backflow_kind ==
              hundun::v04::ScalarBoundaryKind::dirichlet &&
          first.boundaries[1].scalars[0].backflow_value == 0.21 &&
          first.boundaries[1].allow_backflow &&
          first.boundaries[1].relaxation == 0.8,
      rank, "wire publishes scalar catalog pressure reference and boundaries");
  passed &= expect(
      first.solver.pressure.absolute_tolerance == 1.0e-8 &&
          first.solver.pressure.relative_tolerance == 2.0e-7 &&
          first.solver.pressure.maximum_iterations == 333U &&
          first.solver.pressure.true_residual_interval == 7U &&
          first.solver.pressure.krylov_restart == 16U &&
          first.solver.terminal.eos == 3.0e-6 &&
          first.solver.terminal.continuity == 4.0e-6 &&
          first.solver.terminal.closed_mass == 5.0e-6 &&
          first.solver.terminal.gauge == 6.0e-6 &&
          first.schemes.momentum ==
              hundun::v04::ConvectionScheme::limited_central2 &&
          first.schemes.enthalpy == hundun::v04::ConvectionScheme::central2 &&
          first.schemes.limiter == 0.75 &&
          first.time.scheme == hundun::v04::TimeScheme::variable_bdf2 &&
          first.time.maximum_retries == 6U &&
          first.time.minimum_bdf_ratio == 0.25,
      rank, "wire publishes typed schemes and time control");
  passed &= expect(first.data_files.size() == 1U &&
                       first.data_files[0] == "profile.d" &&
                       first.thermophysics.data_file ==
                           fs::path("thermophysics.d") &&
                       first.thermophysics.species.size() == 2U &&
                       first.immersed_boundary.has_value() &&
                       first.immersed_boundary->stl_file ==
                           fs::path("shape.stl") &&
                       first.immersed_boundary->fluid_side ==
                           hundun::v04::ImmersedFluidSide::inside,
                   rank, "relative references survive the wire");
  passed &= expect(rank == 0 ? g_open_count == 4 : g_open_count == 0, rank,
                   "only rank zero opens JSON/data/STL");
  passed &= expect(rank == 0 ? g_observed_rank == 0 : g_observed_rank == -1,
                   rank, "observer records only rank zero");

  const std::uint64_t first_fingerprint = first.fingerprint;
  if (rank == 0) {
    write_file(root / "profile.d", "0 1\n1 3\n");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  reset_observer();
  ValidatedModel mutated;
  const Status mutated_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison, mutated);
  passed &= expect(static_cast<bool>(mutated_status), rank,
                   "mutated referenced data recompiles");
  passed &= expect(same_u64(mutated.fingerprint, MPI_COMM_WORLD), rank,
                   "mutated fingerprint is identical");
  passed &= expect(mutated.fingerprint != first_fingerprint, rank,
                   "referenced bytes affect fingerprint");
  passed &= expect(rank == 0 ? g_open_count == 4 : g_open_count == 0, rank,
                   "recompile remains root-only I/O");

  if (rank == 0) {
    write_file(root / "profile.d", "0 1\n1 2\n");
    std::string changed_thermophysics{kThermophysics};
    const std::string original_viscosity{"1.919e-5"};
    const std::size_t viscosity =
        changed_thermophysics.find(original_viscosity);
    if (viscosity != std::string::npos) {
      changed_thermophysics.replace(viscosity, original_viscosity.size(),
                                    "2.019e-5");
    }
    write_file(root / "thermophysics.d", changed_thermophysics);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  reset_observer();
  ValidatedModel thermophysics_mutated;
  const Status thermophysics_mutated_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison,
      thermophysics_mutated);
  passed &= expect(static_cast<bool>(thermophysics_mutated_status), rank,
                   "mutated thermophysics recompiles");
  passed &= expect(
      same_u64(thermophysics_mutated.fingerprint, MPI_COMM_WORLD), rank,
      "thermophysics-mutated fingerprint is identical");
  passed &= expect(thermophysics_mutated.fingerprint != first_fingerprint &&
                       thermophysics_mutated.thermophysics.species.size() ==
                           2U &&
                       thermophysics_mutated.thermophysics.species[0]
                               .viscosity_reference ==
                           2.019e-5,
                   rank,
                   "thermophysics bytes and typed values affect fingerprint");
  passed &= expect(rank == 0 ? g_open_count == 4 : g_open_count == 0, rank,
                   "thermophysics is opened exactly once on rank zero");

  if (rank == 0) {
    write_file(root / "thermophysics.d", kThermophysics);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    write_file(root / "case.json", "{ invalid json");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  reset_observer();
  ValidatedModel rejected;
  rejected.fingerprint = 777U;
  const Status rejected_status = CaseCompiler::load_and_compile(
      MPI_COMM_WORLD, rank == 0 ? root : nonroot_poison, rejected);
  passed &= expect(rejected_status.code == StatusCode::invalid_case, rank,
                   "invalid root JSON fails collectively");
  passed &= expect(same_status(rejected_status, MPI_COMM_WORLD), rank,
                   "failure category and detail are identical");
  passed &= expect(rejected.fingerprint == 777U, rank,
                   "failed compile does not publish a partial model");
  passed &= expect(
      hundun::v04::detail::last_lowest_failing_rank_for_test() == 0, rank,
      "lowest failing rank is broadcast");
  passed &= expect(rank == 0 ? g_open_count == 1 : g_open_count == 0, rank,
                   "failure parsing remains root-only I/O");

  passed = collective_result(passed, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    std::error_code error;
    fs::remove_all(root, error);
  }
  hundun::v04::detail::set_file_open_observer_for_test(nullptr);

  if (rank == 0 && !passed) {
    std::cerr << "v04 app case broadcast failed for " << size << " ranks\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
