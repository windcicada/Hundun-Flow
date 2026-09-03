// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr double kMolecularWeightA = 24.0;
constexpr double kMolecularWeightB = 32.0;
constexpr double kCpA = 1050.0;
constexpr double kCpB = 1150.0;
constexpr std::size_t kLocalCells = 7U;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool near_relative(double actual, double expected, double tolerance) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             tolerance * std::max(1.0, std::abs(expected));
}

bool identical_double(MPI_Comm communicator, double value) {
  double minimum = 0.0;
  double maximum = 0.0;
  return MPI_Allreduce(&value, &minimum, 1, MPI_DOUBLE, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_DOUBLE, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool identical_status(MPI_Comm communicator, Status status) {
  const std::array<std::uint64_t, 2U> local{
      static_cast<std::uint64_t>(status.code), status.detail};
  std::array<std::uint64_t, 2U> minimum{};
  std::array<std::uint64_t, 2U> maximum{};
  return MPI_Allreduce(local.data(), minimum.data(),
                       static_cast<int>(local.size()), MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(local.data(), maximum.data(),
                       static_cast<int>(local.size()), MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool expected_collective_status(MPI_Comm communicator, Status status,
                                StatusCode expected, int rank,
                                std::string_view description) {
  const bool same = identical_status(communicator, status);
  return expect(status.code == expected && same, rank, description);
}

SpeciesThermophysicalSpec species(std::string_view name, double mw,
                                  double cp) {
  SpeciesThermophysicalSpec value;
  value.stable_name.assign(name.data(), name.size());
  value.molecular_weight = mw;
  value.temperature_switch = 1000.0;
  const double gas_constant = kUniversalGasConstant / mw;
  value.nasa7_low[0U] = cp / gas_constant;
  value.nasa7_high[0U] = cp / gas_constant;
  value.viscosity_reference = 1.0e-5;
  value.conductivity = 0.02;
  return value;
}

ThermophysicalSpec spec() {
  ThermophysicalSpec value;
  value.data_file = "thermo.d";
  value.minimum_temperature = 200.0;
  value.maximum_temperature = 2500.0;
  value.temperature_relative_tolerance = 1.0e-12;
  value.maximum_temperature_iterations = 64U;
  value.closed_mass_relative_tolerance = 1.0e-12;
  value.maximum_closed_mass_iterations = 32U;
  value.maximum_closed_mass_relative_step = 0.2;
  value.species.push_back(species("A", kMolecularWeightA, kCpA));
  value.species.push_back(species("B", kMolecularWeightB, kCpB));
  return value;
}

struct Fixture {
  std::vector<double> temperature;
  std::vector<double> pressure_perturbation;
  std::vector<double> enthalpy;
  std::vector<double> independent_mass_fraction;
  std::vector<double> volume;
  std::vector<std::uint8_t> active;

  ClosedMassCellView view() const noexcept {
    return {
        Span<const double>{pressure_perturbation.data(),
                           pressure_perturbation.size()},
        Span<const double>{enthalpy.data(), enthalpy.size()},
        Span<const double>{independent_mass_fraction.data(),
                           independent_mass_fraction.size()},
        Span<const double>{volume.data(), volume.size()},
        Span<const std::uint8_t>{active.data(), active.size()}};
  }
};

Fixture make_fixture(int rank, int size) {
  Fixture fixture;
  fixture.temperature.resize(kLocalCells);
  fixture.pressure_perturbation.resize(kLocalCells);
  fixture.enthalpy.resize(kLocalCells);
  fixture.independent_mass_fraction.resize(kLocalCells);
  fixture.volume.resize(kLocalCells);
  fixture.active.assign(kLocalCells, 1U);
  for (std::size_t cell = 0U; cell < kLocalCells; ++cell) {
    const double cell_value = static_cast<double>(cell);
    const double rank_value = static_cast<double>(rank);
    const double temperature =
        405.0 + 23.0 * rank_value + 31.0 * cell_value;
    const double y_a = 0.08 + 0.013 * rank_value + 0.027 * cell_value;
    const double cp_mix = y_a * kCpA + (1.0 - y_a) * kCpB;
    fixture.temperature[cell] = temperature;
    fixture.pressure_perturbation[cell] =
        -1730.0 + 117.0 * rank_value + 281.0 * cell_value;
    // Independent constant-cp NASA7 oracle: a6=0, hence h=cp_mix*T.
    fixture.enthalpy[cell] = cp_mix * temperature;
    fixture.independent_mass_fraction[cell] = y_a;
    fixture.volume[cell] =
        0.19 + 0.021 * rank_value + 0.037 * cell_value;
  }
  if (rank == size - 1) {
    fixture.active.back() = 0U;
  }
  return fixture;
}

// This oracle deliberately does not invoke ThermodynamicsPlan. It evaluates
// the known binary ideal-gas mixture directly from T, Y, W, p_ref+pi, and V.
bool oracle_global_mass(MPI_Comm communicator, const Fixture& fixture,
                        double pressure_reference, double& mass) {
  long double local_mass = 0.0L;
  long double correction = 0.0L;
  for (std::size_t cell = 0U; cell < fixture.temperature.size(); ++cell) {
    if (fixture.active[cell] == 0U) {
      continue;
    }
    const long double y_a = fixture.independent_mass_fraction[cell];
    const long double gas_constant =
        static_cast<long double>(kUniversalGasConstant) *
        (y_a / static_cast<long double>(kMolecularWeightA) +
         (1.0L - y_a) / static_cast<long double>(kMolecularWeightB));
    const long double pressure =
        static_cast<long double>(pressure_reference) +
        static_cast<long double>(fixture.pressure_perturbation[cell]);
    const long double contribution =
        pressure * static_cast<long double>(fixture.volume[cell]) /
        (gas_constant * static_cast<long double>(fixture.temperature[cell]));
    const long double corrected = contribution - correction;
    const long double next = local_mass + corrected;
    correction = (next - local_mass) - corrected;
    local_mass = next;
  }
  long double global_mass = 0.0L;
  if (MPI_Allreduce(&local_mass, &global_mass, 1, MPI_LONG_DOUBLE, MPI_SUM,
                    communicator) != MPI_SUCCESS ||
      !std::isfinite(global_mass)) {
    return false;
  }
  mass = static_cast<double>(global_mass);
  return std::isfinite(mass);
}

ClosedMassResult sentinel(double seed) {
  ClosedMassResult value;
  value.pressure_reference = seed;
  value.mass = seed + 1.0;
  value.residual = seed + 2.0;
  value.iterations = static_cast<std::uint32_t>(seed) + 3U;
  value.lowest_failing_rank = static_cast<int>(seed) + 4;
  return value;
}

bool same_result(const ClosedMassResult& left,
                 const ClosedMassResult& right) {
  return left.pressure_reference == right.pressure_reference &&
         left.mass == right.mass && left.residual == right.residual &&
         left.iterations == right.iterations &&
         left.lowest_failing_rank == right.lowest_failing_rank;
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

  const ThermophysicalSpec thermo_spec = spec();
  std::array<TransportedScalarSpec, 1U> catalog{
      TransportedScalarSpec{"A", TransportedScalarRole::species}};
  ThermodynamicsPlan thermodynamics;
  ClosedMassPlan closed;
  ClosedMassPlan open;
  bool passed = expect(static_cast<bool>(ThermodynamicsPlan::compile(
                           thermo_spec,
                           Span<const TransportedScalarSpec>{catalog.data(),
                                                            catalog.size()},
                           thermodynamics)),
                       rank, "closed-mass thermo plan compiles");
  passed &= expect(static_cast<bool>(ClosedMassPlan::compile(
                       PressureReferenceKind::closed_mass, thermo_spec,
                       closed)) &&
                       static_cast<bool>(ClosedMassPlan::compile(
                           PressureReferenceKind::boundary_absolute,
                           thermo_spec, open)),
                   rank, "open and closed pressure policies compile");

  const Fixture fixture = make_fixture(rank, size);
  const ClosedMassCellView cells = fixture.view();
  constexpr double target_pressure = 146000.0;
  constexpr double initial_pressure = 80000.0;
  double target_mass = 0.0;
  passed &= expect(oracle_global_mass(MPI_COMM_WORLD, fixture,
                                     target_pressure, target_mass) &&
                       target_mass > 0.0,
                   rank, "independent analytic EOS oracle produces target mass");

  ClosedMassResult result;
  const Status solve = closed.solve(MPI_COMM_WORLD, thermodynamics, cells,
                                    target_mass, initial_pressure, result);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, solve, StatusCode::ok, rank,
      "collective bounded Newton converges on nonuniform fields");
  passed &= expect(near_relative(result.pressure_reference, target_pressure,
                                 5.0e-12) &&
                       near_relative(result.mass, target_mass, 5.0e-12) &&
                       std::abs(result.residual) <=
                           thermo_spec.closed_mass_relative_tolerance *
                               target_mass &&
                       result.iterations == 5U &&
                       result.lowest_failing_rank == -1,
                   rank,
                   "nonuniform T/Y/pi/V reaches 5e-12 mass and the 20-percent "
                   "bound requires four updates plus the convergence check");

  CartesianMeshSpec direct_mesh;
  direct_mesh.kind = GeometryKind::uniform;
  direct_mesh.lower = {0.0, 0.0, 0.0};
  direct_mesh.upper = {static_cast<double>(kLocalCells * size), 1.0, 1.0};
  direct_mesh.has_exact_cells = true;
  direct_mesh.exact_cells = {
      static_cast<std::int32_t>(kLocalCells * size), 1, 1};
  direct_mesh.minimum_spacing = {1.0, 1.0, 1.0};
  direct_mesh.max_growth_ratio = 1.0;
  direct_mesh.limits.max_global_cells = kLocalCells * size;
  direct_mesh.limits.max_memory_bytes_per_rank = 1U << 24U;
  CartesianGeometryPlan direct_geometry;
  MeshPatch direct_patch;
  passed &= expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, direct_mesh, {}, direct_geometry, direct_patch)) &&
          direct_patch.cells.x == static_cast<std::int32_t>(kLocalCells) &&
          direct_patch.cells.y == 1 && direct_patch.cells.z == 1,
      rank, "direct-view closed-mass geometry decomposes without packing");
  const auto direct_view = [&](const std::vector<double>& values,
                               FieldId field, RevisionToken revision,
                               StorageIdentity storage) {
    return ConstFieldView{values.data(), direct_patch.cells, {0, 0, 0}, 1U,
                          static_cast<std::size_t>(direct_patch.cells.x),
                          static_cast<std::size_t>(direct_patch.cells.x),
                          values.size(), 0U, field, revision, storage, 991U};
  };
  const ConstFieldView direct_pi =
      direct_view(fixture.pressure_perturbation, 20U, 901U, 1901U);
  const ConstFieldView direct_h =
      direct_view(fixture.enthalpy, 21U, 902U, 1902U);
  const ConstFieldView direct_y =
      direct_view(fixture.independent_mass_fraction, 22U, 903U, 1903U);
  const std::array<ConstFieldView, 1U> direct_species{direct_y};
  ClosedMassFieldView direct_cells{
      direct_pi,
      direct_h,
      {direct_species.data(), direct_species.size()},
      &direct_geometry,
      direct_patch,
      {fixture.active.data(), fixture.active.size()},
      904U};
  Fixture unit_volume = fixture;
  std::fill(unit_volume.volume.begin(), unit_volume.volume.end(), 1.0);
  double direct_target_mass = 0.0;
  passed &= expect(oracle_global_mass(MPI_COMM_WORLD, unit_volume,
                                     target_pressure, direct_target_mass),
                   rank, "direct-view target uses geometry-owned volume");
  ClosedMassResult direct_result;
  const Status direct_status = closed.solve_fields(
      MPI_COMM_WORLD, thermodynamics, direct_cells, direct_target_mass,
      initial_pressure, direct_result);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, direct_status, StatusCode::ok, rank,
      "strided field-view closed-mass solve succeeds collectively");
  passed &= expect(
      near_relative(direct_result.pressure_reference, target_pressure,
                    5.0e-12) &&
          near_relative(direct_result.mass, direct_target_mass, 5.0e-12),
      rank, "field-view solve matches the independent EOS oracle");
  passed &= expect(identical_double(MPI_COMM_WORLD,
                                    result.pressure_reference),
                   rank, "all ranks publish one bitwise pressure result");

  // With an affine EOS an unbounded Newton step would converge on iteration
  // two. Requiring failure with only two iterations proves the 20% first-step
  // bound is active (80000 -> at most 96000), rather than a direct jump.
  ThermophysicalSpec two_iteration_spec = thermo_spec;
  two_iteration_spec.maximum_closed_mass_iterations = 2U;
  ClosedMassPlan two_iteration_plan;
  passed &= expect(static_cast<bool>(ClosedMassPlan::compile(
                       PressureReferenceKind::closed_mass,
                       two_iteration_spec, two_iteration_plan)),
                   rank, "two-iteration bounded plan compiles");
  ClosedMassResult unchanged = sentinel(31.0);
  const ClosedMassResult expected_unchanged = unchanged;
  const Status bounded_status = two_iteration_plan.solve(
      MPI_COMM_WORLD, thermodynamics, cells, target_mass, initial_pressure,
      unchanged);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, bounded_status, StatusCode::numerical_failure, rank,
      "two iterations cannot bypass the configured relative-step bound");
  passed &= expect(same_result(unchanged, expected_unchanged), rank,
                   "bounded non-convergence leaves every output member atomic");

  unchanged = sentinel(41.0);
  const ClosedMassResult expected_open = unchanged;
  const Status open_status = open.solve(MPI_COMM_WORLD, thermodynamics, cells,
                                        target_mass, initial_pressure,
                                        unchanged);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, open_status, StatusCode::invalid_plan, rank,
      "open pressure authority cannot invoke mass correction");
  passed &= expect(same_result(unchanged, expected_open), rank,
                   "invalid pressure authority preserves the output atomically");

  const int mutation_rank = std::min(1, size - 1);

  Fixture invalid_composition = fixture;
  if (rank == mutation_rank) {
    invalid_composition.independent_mass_fraction[2U] = 1.125;
  }
  unchanged = sentinel(51.0);
  const ClosedMassResult expected_composition = unchanged;
  const Status composition_status = closed.solve(
      MPI_COMM_WORLD, thermodynamics, invalid_composition.view(), target_mass,
      initial_pressure, unchanged);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, composition_status, StatusCode::numerical_failure, rank,
      "one-rank invalid composition fails collectively");
  passed &= expect(same_result(unchanged, expected_composition), rank,
                   "composition failure leaves every output member atomic");

  Fixture invalid_active = fixture;
  if (rank == mutation_rank) {
    invalid_active.active[1U] = 2U;
  }
  unchanged = sentinel(61.0);
  const ClosedMassResult expected_active = unchanged;
  const Status active_status = closed.solve(
      MPI_COMM_WORLD, thermodynamics, invalid_active.view(), target_mass,
      initial_pressure, unchanged);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, active_status, StatusCode::invalid_plan, rank,
      "one-rank invalid active marker fails collectively");
  passed &= expect(same_result(unchanged, expected_active), rank,
                   "active-marker failure leaves every output member atomic");

  ClosedMassCellView invalid_shape = cells;
  if (rank == mutation_rank) {
    invalid_shape.volume =
        Span<const double>{fixture.volume.data(), fixture.volume.size() - 1U};
  }
  unchanged = sentinel(71.0);
  const ClosedMassResult expected_shape = unchanged;
  const Status shape_status = closed.solve(
      MPI_COMM_WORLD, thermodynamics, invalid_shape, target_mass,
      initial_pressure, unchanged);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, shape_status, StatusCode::invalid_plan, rank,
      "one-rank shape mutation fails collectively");
  passed &= expect(same_result(unchanged, expected_shape), rank,
                   "shape failure leaves every output member atomic");

  // The initial p_abs remains positive, but the requested final pressure makes
  // one active cell negative. This exercises the post-Newton collective gate.
  Fixture invalid_final_pressure = fixture;
  if (rank == mutation_rank) {
    invalid_final_pressure.pressure_perturbation[0U] = -60000.0;
    invalid_final_pressure.volume[0U] = 1.0e-4;
    invalid_final_pressure.active[0U] = 1U;
  }
  constexpr double final_invalid_target_pressure = 50000.0;
  constexpr double final_invalid_initial_pressure = 100000.0;
  double final_invalid_target_mass = 0.0;
  passed &= expect(oracle_global_mass(MPI_COMM_WORLD, invalid_final_pressure,
                                     final_invalid_target_pressure,
                                     final_invalid_target_mass) &&
                       final_invalid_target_mass > 0.0,
                   rank, "final-p_abs mutation has a positive global target mass");
  unchanged = sentinel(81.0);
  const ClosedMassResult expected_final_pressure = unchanged;
  const Status final_pressure_status = closed.solve(
      MPI_COMM_WORLD, thermodynamics, invalid_final_pressure.view(),
      final_invalid_target_mass, final_invalid_initial_pressure, unchanged);
  passed &= expected_collective_status(
      MPI_COMM_WORLD, final_pressure_status, StatusCode::numerical_failure,
      rank, "one-rank nonpositive final p_abs fails collectively");
  passed &= expect(same_result(unchanged, expected_final_pressure), rank,
                   "final-p_abs failure leaves every output member atomic");

  if (size > 1) {
    const double mismatched_target =
        target_mass + (rank == mutation_rank ? 1.0 : 0.0);
    unchanged = sentinel(91.0);
    const ClosedMassResult expected_target = unchanged;
    const Status target_status = closed.solve(
        MPI_COMM_WORLD, thermodynamics, cells, mismatched_target,
        initial_pressure, unchanged);
    passed &= expected_collective_status(
        MPI_COMM_WORLD, target_status, StatusCode::invalid_plan, rank,
        "rank-mismatched target is rejected without deadlock");
    passed &= expect(same_result(unchanged, expected_target), rank,
                     "target identity failure preserves output atomically");

    const double mismatched_initial =
        initial_pressure + (rank == mutation_rank ? 1.0 : 0.0);
    unchanged = sentinel(101.0);
    const ClosedMassResult expected_initial = unchanged;
    const Status initial_status = closed.solve(
        MPI_COMM_WORLD, thermodynamics, cells, target_mass,
        mismatched_initial, unchanged);
    passed &= expected_collective_status(
        MPI_COMM_WORLD, initial_status, StatusCode::invalid_plan, rank,
        "rank-mismatched initial pressure is rejected without deadlock");
    passed &= expect(same_result(unchanged, expected_initial), rank,
                     "initial-pressure identity failure preserves output");

    ThermophysicalSpec local_controls = thermo_spec;
    if (rank == mutation_rank) {
      local_controls.maximum_closed_mass_relative_step = 0.25;
    }
    ClosedMassPlan mismatched_plan;
    passed &= expect(static_cast<bool>(ClosedMassPlan::compile(
                         PressureReferenceKind::closed_mass, local_controls,
                         mismatched_plan)),
                     rank, "rank-local control mutation compiles");
    unchanged = sentinel(111.0);
    const ClosedMassResult expected_plan = unchanged;
    const Status plan_status = mismatched_plan.solve(
        MPI_COMM_WORLD, thermodynamics, cells, target_mass, initial_pressure,
        unchanged);
    passed &= expected_collective_status(
        MPI_COMM_WORLD, plan_status, StatusCode::invalid_plan, rank,
        "rank-mismatched ClosedMass fingerprint is rejected without deadlock");
    passed &= expect(same_result(unchanged, expected_plan), rank,
                     "control identity failure preserves output atomically");

    ThermophysicalSpec local_thermo_spec = thermo_spec;
    if (rank == mutation_rank) {
      constexpr double mutated_molecular_weight_b = 31.5;
      local_thermo_spec.species[1U].molecular_weight =
          mutated_molecular_weight_b;
      local_thermo_spec.species[1U].nasa7_low[0U] =
          kCpB / (kUniversalGasConstant / mutated_molecular_weight_b);
      local_thermo_spec.species[1U].nasa7_high[0U] =
          local_thermo_spec.species[1U].nasa7_low[0U];
    }
    ThermodynamicsPlan mismatched_thermodynamics;
    passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                         local_thermo_spec,
                         Span<const TransportedScalarSpec>{catalog.data(),
                                                          catalog.size()},
                         mismatched_thermodynamics)),
                     rank, "rank-local thermodynamic mutation compiles");
    unchanged = sentinel(121.0);
    const ClosedMassResult expected_thermodynamics = unchanged;
    const Status thermodynamics_status = closed.solve(
        MPI_COMM_WORLD, mismatched_thermodynamics, cells, target_mass,
        initial_pressure, unchanged);
    passed &= expected_collective_status(
        MPI_COMM_WORLD, thermodynamics_status, StatusCode::invalid_plan, rank,
        "rank-mismatched thermodynamics fingerprint is rejected");
    passed &= expect(same_result(unchanged, expected_thermodynamics), rank,
                     "thermodynamics identity failure preserves output");
  }

  int local = passed ? 1 : 0;
  int all = 0;
  MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0 && all != 0) {
    std::cout << "v0.4 closed-mass tests passed\n";
  }
  MPI_Finalize();
  return all != 0 ? 0 : 1;
}
