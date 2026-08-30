// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"

#include "../support/product_fixture.hpp"
#include "core_product_freeze_detail.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr Int3 kImmersedCells{16, 16, 16};

std::uint64_t wire_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool same_bits(double left, double right) noexcept {
  return wire_bits(left) == wire_bits(right);
}

bool collective(bool local) {
  const int value = local ? 1 : 0;
  int global = 0;
  return MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ==
             MPI_SUCCESS &&
         global != 0;
}

bool any_rank(bool local) {
  const int value = local ? 1 : 0;
  int global = 0;
  return MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ==
             MPI_SUCCESS &&
         global != 0;
}

bool expect(bool condition, int rank, const char *description) {
  if (!condition)
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  return condition;
}

void print_fresh_diagnostic(int rank) {
  detail::FreshInitializationDiagnostic diagnostic;
  if (!detail::fresh_initialization_diagnostic_for_test(diagnostic)) {
    std::cerr << "rank " << rank << " Fresh diagnostic=unavailable\n";
    return;
  }
  std::cerr << "rank " << rank << " Fresh diagnostic status="
            << static_cast<unsigned>(diagnostic.terminal_status.code) << '/'
            << diagnostic.terminal_status.detail
            << " attempted/audited/commit=" << diagnostic.projection_attempted
            << '/' << diagnostic.audited << '/' << diagnostic.committed
            << " lineages=" << diagnostic.prepared_lineage << '/'
            << diagnostic.solved_lineage << '/' << diagnostic.candidate_lineage
            << " solve=" << static_cast<unsigned>(diagnostic.solve.status.code)
            << '/' << diagnostic.solve.status.detail << '/'
            << static_cast<unsigned>(diagnostic.solve.termination) << '/'
            << diagnostic.solve.iterations
            << " continuity=" << diagnostic.initial_continuity_maximum << '/'
            << diagnostic.final_continuity_maximum << '/'
            << diagnostic.continuity_limit << " envelope/cut/cut-neg/solid="
            << diagnostic.maximum_face_envelope << '/'
            << diagnostic.cut_face_nonzero_count << '/'
            << diagnostic.cut_face_negative_zero_count << '/'
            << diagnostic.solid_nonpositive_zero_component_count
            << " layers=" << diagnostic.velocity_layers_bitwise_equal << '\n';
}

bool valid_ibm_success_diagnostic(PlanFingerprint driver_plan, int rank) {
  detail::FreshInitializationDiagnostic diagnostic;
  const bool available =
      detail::fresh_initialization_diagnostic_for_test(diagnostic);
  const bool linear_success =
      available && diagnostic.solve.status &&
      (diagnostic.solve.termination == LinearTermination::converged ||
       diagnostic.solve.termination == LinearTermination::zero_rhs);
  const bool local =
      available && diagnostic.driver_plan == driver_plan &&
      diagnostic.generation != 0U && diagnostic.immersed &&
      diagnostic.projection_attempted && !diagnostic.no_ibm_bypassed &&
      diagnostic.audited && diagnostic.committed &&
      diagnostic.terminal_status && diagnostic.red_plan != 0U &&
      diagnostic.prepared_lineage != 0U && diagnostic.solved_lineage != 0U &&
      diagnostic.candidate_lineage != 0U && linear_success &&
      std::isfinite(diagnostic.initial_continuity_maximum) &&
      diagnostic.initial_continuity_maximum > 0.0 &&
      std::isfinite(diagnostic.final_continuity_maximum) &&
      diagnostic.final_continuity_maximum >= 0.0 &&
      std::isfinite(diagnostic.continuity_limit) &&
      diagnostic.continuity_limit > 0.0 &&
      diagnostic.final_continuity_maximum <= diagnostic.continuity_limit &&
      std::isfinite(diagnostic.maximum_face_envelope) &&
      diagnostic.maximum_face_envelope > 0.0 &&
      diagnostic.cut_face_nonzero_count == 0U &&
      diagnostic.cut_face_negative_zero_count == 0U &&
      diagnostic.solid_nonpositive_zero_component_count == 0U &&
      diagnostic.velocity_layers_bitwise_equal &&
      diagnostic.derived_velocity_dependents_rebuilt &&
      diagnostic.derived_velocity_lineage != 0U &&
      diagnostic.maximum_h_by_a_velocity_difference == 0.0 &&
      diagnostic.velocity_gradient_nonfinite_count == 0U &&
      diagnostic.effective_viscosity_nonpositive_count == 0U &&
      diagnostic.solid_velocity_gradient_nonfinite_count == 0U &&
      diagnostic.velocity_dependent_rate_layers_bitwise_equal;
  const bool passed = collective(local);
  if (!local)
    print_fresh_diagnostic(rank);
  return expect(passed, rank,
                "distributed IBM Fresh diagnostic certifies continuity, "
                "cut flux, three U layers, and RED lineage");
}

bool valid_no_ibm_bypass_diagnostic(int rank) {
  detail::FreshInitializationDiagnostic diagnostic;
  const bool available =
      detail::fresh_initialization_diagnostic_for_test(diagnostic);
  const bool local =
      available && !diagnostic.immersed && !diagnostic.projection_attempted &&
      diagnostic.no_ibm_bypassed && diagnostic.audited &&
      diagnostic.committed && diagnostic.terminal_status &&
      diagnostic.red_plan != 0U && diagnostic.prepared_lineage != 0U &&
      diagnostic.solved_lineage != 0U && diagnostic.candidate_lineage != 0U &&
      diagnostic.solve.status &&
      diagnostic.solve.termination == LinearTermination::zero_rhs &&
      diagnostic.initial_continuity_maximum == 0.0 &&
      diagnostic.final_continuity_maximum == 0.0 &&
      diagnostic.cut_face_nonzero_count == 0U &&
      diagnostic.solid_nonpositive_zero_component_count == 0U &&
      diagnostic.velocity_layers_bitwise_equal &&
      !diagnostic.derived_velocity_dependents_rebuilt &&
      diagnostic.derived_velocity_lineage == 0U;
  const bool passed = collective(local);
  if (!local)
    print_fresh_diagnostic(rank);
  return expect(passed, rank,
                "distributed non-IBM Fresh diagnostic certifies a literal "
                "bitwise bypass");
}

ValidatedModel immersed_model() {
  ValidatedModel model = hundun::v04::test::product_model(kImmersedCells);
  model.mesh.lower = {-2.0, -2.0, -2.0};
  model.mesh.upper = {2.0, 2.0, 2.0};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  model.time.initial_dt = 2.5e-4;
  model.immersed_boundary =
      ImmersedBoundarySpec{"cylinder_ascii.stl", ImmersedFluidSide::outside};
  return model;
}

std::filesystem::path data_root() {
  return std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";
}

bool cube_solid(const ValidatedModel &model, const MeshPatch &patch,
                Int3 local) noexcept {
  const auto centre = [](double lower, double upper, std::int32_t count,
                         std::int32_t global) noexcept {
    return lower + (static_cast<double>(global) + 0.5) * (upper - lower) /
                       static_cast<double>(count);
  };
  const double x = centre(model.mesh.lower.x, model.mesh.upper.x,
                          model.mesh.exact_cells.x, patch.begin.x + local.x);
  const double y = centre(model.mesh.lower.y, model.mesh.upper.y,
                          model.mesh.exact_cells.y, patch.begin.y + local.y);
  const double z = centre(model.mesh.lower.z, model.mesh.upper.z,
                          model.mesh.exact_cells.z, patch.begin.z + local.z);
  return std::abs(x) < 1.0 && std::abs(y) < 1.0 && std::abs(z) < 1.0;
}

const SnapshotFieldView *find_field(const CommittedOutputSnapshot &snapshot,
                                    std::string_view name) noexcept {
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index)
    if (snapshot.fields.data[index].stable_name == name)
      return &snapshot.fields.data[index];
  return nullptr;
}

bool copy_output_field(const SnapshotFieldView &source,
                       std::vector<double> &target) {
  const ConstFieldView values = source.values;
  if (values.base == nullptr || values.components == 0U)
    return false;
  target.clear();
  target.reserve(static_cast<std::size_t>(values.interior.x) *
                 values.interior.y * values.interior.z * values.components);
  for (std::int32_t z = 0; z < values.interior.z; ++z)
    for (std::int32_t y = 0; y < values.interior.y; ++y)
      for (std::int32_t x = 0; x < values.interior.x; ++x)
        for (std::uint8_t component = 0U; component < values.components;
             ++component)
          target.push_back(values.unchecked({x, y, z}, component));
  return true;
}

bool make_restart_image(const ValidatedModel &model,
                        const CommittedOutputSnapshot &fresh,
                        const RestartExpected &expected, RestartImage &image) {
  const SnapshotFieldView *enthalpy = find_field(fresh, "h");
  const SnapshotFieldView *pressure = find_field(fresh, "pi");
  if (enthalpy == nullptr || pressure == nullptr)
    return false;
  std::vector<double> enthalpy_values;
  std::vector<double> pressure_values;
  if (!copy_output_field(*enthalpy, enthalpy_values) ||
      !copy_output_field(*pressure, pressure_values))
    return false;

  const MeshPatch patch = expected.target_patch;
  const std::size_t cell_count =
      static_cast<std::size_t>(patch.cells.x) * patch.cells.y * patch.cells.z;
  RestartImage candidate;
  candidate.global_cells = expected.global_cells;
  candidate.patch = patch;
  candidate.plan = expected.plan;
  candidate.schema = expected.schema;
  candidate.geometry = expected.geometry;
  candidate.time = model.time.initial_dt;
  candidate.dt = model.time.initial_dt;
  candidate.pressure_reference = 98000.0;
  candidate.step = 1U;
  candidate.backward_euler_recovery = true;
  candidate.fields.reserve(expected.fields.size);
  for (std::size_t index = 0U; index < expected.fields.size; ++index) {
    const RestartExpectedField descriptor = expected.fields.data[index];
    RestartImageField field;
    field.role = descriptor.role;
    field.field = descriptor.field;
    field.components = descriptor.components;
    field.values.resize(cell_count * descriptor.components);
    switch (descriptor.role) {
    case RestartFieldRole::velocity:
      if (descriptor.components != 3U)
        return false;
      for (std::size_t cell = 0U; cell < cell_count; ++cell) {
        field.values[cell * 3U] = 0.375;
        field.values[cell * 3U + 1U] = -0.25;
        field.values[cell * 3U + 2U] = 0.125;
      }
      break;
    case RestartFieldRole::pressure_perturbation:
      if (descriptor.components != 1U || pressure_values.size() != cell_count)
        return false;
      field.values = pressure_values;
      break;
    case RestartFieldRole::enthalpy:
      if (descriptor.components != 1U || enthalpy_values.size() != cell_count)
        return false;
      field.values = enthalpy_values;
      break;
    default:
      return false;
    }
    candidate.fields.push_back(std::move(field));
  }
  const Int3 cells = patch.cells;
  const std::array<std::size_t, 3U> face_counts{
      static_cast<std::size_t>(cells.x + 1) * cells.y * cells.z,
      static_cast<std::size_t>(cells.x) * (cells.y + 1) * cells.z,
      static_cast<std::size_t>(cells.x) * cells.y * (cells.z + 1)};
  for (std::size_t axis = 0U; axis < face_counts.size(); ++axis) {
    candidate.final_mass_flux[axis].resize(face_counts[axis]);
    for (std::size_t index = 0U; index < face_counts[axis]; ++index)
      candidate.final_mass_flux[axis][index] =
          ((index + axis) & 1U) == 0U ? 0.0 : -0.0;
  }
  image = std::move(candidate);
  return true;
}

bool distributed_restart_skips_fresh(const ValidatedModel &model,
                                     const CommittedOutputSnapshot &fresh,
                                     int rank) {
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, data_root(), plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  RestartExpected expected;
  if (status)
    status = driver.restart_expected(expected);
  RestartImage image;
  if (status && !make_restart_image(model, fresh, expected, image))
    status = {StatusCode::invalid_plan, 1U};
  detail::clear_fresh_initialization_diagnostic_for_test();
  if (status)
    status = driver.initialize_restart(image);
  bool passed = collective(static_cast<bool>(status));
  passed &= expect(passed, rank,
                   "distributed IBM Restart accepts the checkpoint image");
  if (!passed)
    return false;

  CommittedOutputSnapshot output;
  status = driver.committed_output_snapshot(output);
  const SnapshotFieldView *velocity =
      status ? find_field(output, "U") : nullptr;
  const RestartImageField *expected_velocity = nullptr;
  for (const RestartImageField &field : image.fields)
    if (field.role == RestartFieldRole::velocity)
      expected_velocity = &field;
  bool local_exact =
      status && velocity != nullptr && expected_velocity != nullptr;
  if (velocity != nullptr && expected_velocity != nullptr) {
    std::size_t cell = 0U;
    for (std::int32_t z = 0; z < velocity->values.interior.z; ++z)
      for (std::int32_t y = 0; y < velocity->values.interior.y; ++y)
        for (std::int32_t x = 0; x < velocity->values.interior.x; ++x, ++cell)
          for (std::uint8_t component = 0U; component < 3U; ++component)
            local_exact &=
                same_bits(velocity->values.unchecked({x, y, z}, component),
                          expected_velocity->values[cell * 3U + component]);
  }

  RestartSnapshot restart;
  if (status)
    status = driver.committed_restart_snapshot(restart);
  local_exact &= static_cast<bool>(status);
  if (status) {
    const std::array<ConstFaceFieldView, 3U> faces{restart.final_mass_flux.x,
                                                   restart.final_mass_flux.y,
                                                   restart.final_mass_flux.z};
    for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
      std::size_t face = 0U;
      for (std::int32_t z = 0; z < faces[axis].extents.z; ++z)
        for (std::int32_t y = 0; y < faces[axis].extents.y; ++y)
          for (std::int32_t x = 0; x < faces[axis].extents.x; ++x, ++face)
            local_exact &= same_bits(faces[axis].unchecked({x, y, z}),
                                     image.final_mass_flux[axis][face]);
    }
  }
  detail::FreshInitializationDiagnostic unexpected;
  local_exact &= !detail::fresh_initialization_diagnostic_for_test(unexpected);
  detail::ColdVelocityDependentsDiagnostic derived;
  local_exact &=
      detail::cold_velocity_dependents_diagnostic_for_test(derived) &&
      derived.restart && derived.rebuilt &&
      derived.driver_plan == expected.plan && derived.lineage != 0U &&
      derived.maximum_h_by_a_velocity_difference == 0.0 &&
      derived.velocity_gradient_nonfinite_count == 0U &&
      derived.effective_viscosity_nonpositive_count == 0U &&
      derived.solid_velocity_gradient_nonfinite_count == 0U &&
      derived.rate_layers_bitwise_equal;
  passed = collective(local_exact);
  return expect(passed, rank,
                "distributed Restart skips Fresh and restores U/final-flux "
                "bytes exactly");
}

bool distributed_ibm_fresh(int rank) {
  const ValidatedModel model = immersed_model();
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, data_root(), plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.25, -0.125};
  if (status)
    status = driver.initialize(initial);
  if (!status)
    std::cerr << "rank " << rank << " IBM Fresh initialize status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
  if (!status)
    print_fresh_diagnostic(rank);
  bool passed = collective(static_cast<bool>(status));
  passed &=
      expect(passed, rank, "distributed IBM Fresh initialization succeeds");
  if (!passed)
    return false;

  CommittedOutputSnapshot snapshot;
  status = driver.committed_output_snapshot(snapshot);
  passed = collective(static_cast<bool>(status));
  passed &= expect(passed, rank,
                   "distributed IBM Fresh step-zero snapshot is available");
  if (!passed)
    return false;
  const SnapshotFieldView *velocity = find_field(snapshot, "U");
  const SnapshotFieldView *pressure = find_field(snapshot, "pi");
  const bool views_valid =
      velocity != nullptr && velocity->values.components == 3U &&
      pressure != nullptr && pressure->values.components == 1U;
  passed = collective(views_valid);
  passed &= expect(passed, rank, "distributed IBM Fresh snapshot has U and pi");
  if (!passed)
    return false;

  const std::array<double, 3U> base{initial.velocity.x, initial.velocity.y,
                                    initial.velocity.z};
  std::uint64_t local_solid = 0U;
  bool local_wall_zero = true;
  bool local_fluid_changed = false;
  bool local_pressure_untouched = true;
  for (std::int32_t z = 0; z < snapshot.patch.cells.z; ++z)
    for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y)
      for (std::int32_t x = 0; x < snapshot.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const bool solid = cube_solid(model, snapshot.patch, cell);
        local_solid += solid ? 1U : 0U;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double value = velocity->values.unchecked(cell, component);
          if (solid)
            local_wall_zero &= same_bits(value, 0.0);
          else
            local_fluid_changed |= !same_bits(value, base[component]);
        }
        local_pressure_untouched &=
            same_bits(pressure->values.unchecked(cell, 0U), 0.0);
      }
  std::uint64_t global_solid = 0U;
  const bool reduced =
      MPI_Allreduce(&local_solid, &global_solid, 1, MPI_UINT64_T, MPI_SUM,
                    MPI_COMM_WORLD) == MPI_SUCCESS;
  const bool global_wall_zero = collective(local_wall_zero);
  const bool global_fluid_changed = any_rank(local_fluid_changed);
  const bool global_pressure_untouched = collective(local_pressure_untouched);
  passed &= expect(reduced && global_solid != 0U, rank,
                   "distributed IBM fixture contains solid cells");
  passed &= expect(global_wall_zero, rank,
                   "distributed IBM Fresh writes exact stationary-wall +0 "
                   "on every rank");
  passed &= expect(global_fluid_changed, rank,
                   "distributed IBM Fresh projects the fluid plug across "
                   "the rank decomposition");
  passed &= expect(global_pressure_untouched, rank,
                   "distributed IBM Fresh leaves pressure bitwise untouched");
  passed &= valid_ibm_success_diagnostic(snapshot.plan, rank);
  passed &= distributed_restart_skips_fresh(model, snapshot, rank);
  return collective(passed);
}

bool distributed_no_ibm_bypass(int rank) {
  const ValidatedModel model = hundun::v04::test::product_model({12, 8, 6});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {-0.0, 0.375, -0.625};
  if (status)
    status = driver.initialize(initial);
  CommittedOutputSnapshot snapshot;
  if (status)
    status = driver.committed_output_snapshot(snapshot);
  bool passed = collective(static_cast<bool>(status));
  passed &=
      expect(passed, rank, "distributed non-IBM Fresh initialization succeeds");
  if (!passed)
    return false;
  const SnapshotFieldView *velocity = find_field(snapshot, "U");
  passed = collective(velocity != nullptr && velocity->values.components == 3U);
  if (!passed)
    return false;
  const std::array<double, 3U> expected{initial.velocity.x, initial.velocity.y,
                                        initial.velocity.z};
  bool local_exact = true;
  for (std::int32_t z = 0; z < snapshot.patch.cells.z; ++z)
    for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y)
      for (std::int32_t x = 0; x < snapshot.patch.cells.x; ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component)
          local_exact &=
              same_bits(velocity->values.unchecked({x, y, z}, component),
                        expected[component]);
  passed = collective(local_exact);
  passed &=
      expect(passed, rank, "distributed non-IBM Fresh is a bitwise U bypass");
  passed &= valid_no_ibm_bypass_diagnostic(rank);
  return passed;
}

bool distributed_projection_failure_retry(int rank) {
  const ValidatedModel model = immersed_model();
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_WORLD, model, data_root(), plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  if (!collective(static_cast<bool>(status)))
    return false;

  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.25, -0.125};
  detail::arm_fresh_initialization_candidate_poison_once_for_test(
      0, detail::FreshInitializationPoisonKind::flux);
  status = driver.initialize(initial);
  CommittedOutputSnapshot unpublished;
  const Status output_status = driver.committed_output_snapshot(unpublished);
  detail::FreshInitializationDiagnostic rejected;
  const bool rejected_local =
      !status && !driver.initialized() && !output_status &&
      driver.pressure_reference() == 0.0 &&
      driver.closed_mass_target() == 0.0 &&
      detail::fresh_initialization_diagnostic_for_test(rejected) &&
      rejected.immersed && rejected.projection_attempted && rejected.audited &&
      !rejected.committed && !rejected.terminal_status &&
      rejected.prepared_lineage != 0U && rejected.solved_lineage != 0U &&
      rejected.candidate_lineage != 0U &&
      rejected.velocity_layers_bitwise_equal;
  bool passed = collective(rejected_local);
  if (!rejected_local)
    print_fresh_diagnostic(rank);
  passed &= expect(
      passed, rank,
      "rank-0 candidate poison is rejected collectively before Fresh publish");
  if (!passed)
    return false;

  status = driver.initialize(initial);
  detail::FreshInitializationDiagnostic recovered;
  const bool recovered_local =
      status && driver.initialized() &&
      detail::fresh_initialization_diagnostic_for_test(recovered) &&
      recovered.committed && recovered.terminal_status &&
      recovered.generation > rejected.generation &&
      recovered.derived_velocity_dependents_rebuilt &&
      recovered.derived_velocity_lineage != 0U &&
      recovered.velocity_dependent_rate_layers_bitwise_equal;
  const bool recovered_collective = collective(recovered_local);
  passed &= expect(recovered_collective, rank,
                   "collective poison is consumed and the same Fresh input "
                   "retries successfully");
  return passed;
}

bool distributed_failure_retry(int rank) {
  const ValidatedModel model = hundun::v04::test::product_model({12, 8, 6});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
  if (!collective(static_cast<bool>(status)))
    return false;

  DriverInitialState invalid;
  invalid.pressure_reference = 98000.0;
  invalid.temperature = std::numeric_limits<double>::quiet_NaN();
  status = driver.initialize(invalid);
  CommittedOutputSnapshot output;
  const Status output_status = driver.committed_output_snapshot(output);
  bool passed =
      collective(!status && status.code == StatusCode::invalid_plan &&
                 !driver.initialized() && driver.pressure_reference() == 0.0 &&
                 driver.closed_mass_target() == 0.0 && !output_status);
  passed &= expect(passed, rank,
                   "distributed failed Fresh publishes no committed state");

  DriverInitialState valid;
  valid.pressure_reference = 98000.0;
  valid.temperature = 315.0;
  valid.velocity = {0.125, -0.25, 0.5};
  status = driver.initialize(valid);
  const bool retried =
      collective(static_cast<bool>(status) && driver.initialized());
  passed &= expect(retried, rank,
                   "distributed failed Fresh can be retried collectively");
  return passed;
}

bool run(int rank) {
  bool passed = distributed_ibm_fresh(rank);
  passed &= distributed_projection_failure_retry(rank);
  passed &= distributed_no_ibm_bypass(rank);
  passed &= distributed_failure_retry(rank);
  return collective(passed);
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const bool passed = run(rank);
  MPI_Finalize();
  return passed ? 0 : 1;
}
