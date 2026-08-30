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
constexpr Real3 kRestartVelocity{0.375, -0.25, 0.125};

std::uint64_t wire_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool same_bits(double left, double right) noexcept {
  return wire_bits(left) == wire_bits(right);
}

bool expect(bool condition, const char *description) {
  if (!condition)
    std::cerr << "FAIL: " << description << '\n';
  return condition;
}

void print_fresh_diagnostic(const char *prefix) {
  detail::FreshInitializationDiagnostic diagnostic;
  if (!detail::fresh_initialization_diagnostic_for_test(diagnostic)) {
    std::cerr << prefix << " diagnostic=unavailable\n";
    return;
  }
  std::cerr << prefix << " diagnostic status="
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

bool valid_ibm_success_diagnostic(PlanFingerprint driver_plan) {
  detail::FreshInitializationDiagnostic diagnostic;
  bool passed =
      expect(detail::fresh_initialization_diagnostic_for_test(diagnostic),
             "IBM Fresh publishes its test-only transaction diagnostic");
  if (!passed)
    return false;
  const bool linear_success =
      diagnostic.solve.status &&
      (diagnostic.solve.termination == LinearTermination::converged ||
       diagnostic.solve.termination == LinearTermination::zero_rhs);
  passed &= expect(
      diagnostic.driver_plan == driver_plan && diagnostic.generation != 0U &&
          diagnostic.immersed && diagnostic.projection_attempted &&
          !diagnostic.no_ibm_bypassed && diagnostic.audited &&
          diagnostic.committed && diagnostic.terminal_status &&
          diagnostic.red_plan != 0U && diagnostic.prepared_lineage != 0U &&
          diagnostic.solved_lineage != 0U &&
          diagnostic.candidate_lineage != 0U && linear_success,
      "IBM Fresh diagnostic carries complete RED/solve/transaction lineage");
  passed &= expect(
      std::isfinite(diagnostic.initial_continuity_maximum) &&
          diagnostic.initial_continuity_maximum > 0.0 &&
          std::isfinite(diagnostic.final_continuity_maximum) &&
          diagnostic.final_continuity_maximum >= 0.0 &&
          std::isfinite(diagnostic.continuity_limit) &&
          diagnostic.continuity_limit > 0.0 &&
          diagnostic.final_continuity_maximum <= diagnostic.continuity_limit &&
          std::isfinite(diagnostic.maximum_face_envelope) &&
          diagnostic.maximum_face_envelope > 0.0,
      "IBM Fresh fluid continuity satisfies the certified normalized gate");
  passed &= expect(diagnostic.cut_face_nonzero_count == 0U &&
                       diagnostic.cut_face_negative_zero_count == 0U,
                   "IBM Fresh cut-face final flux is exact positive zero");
  passed &= expect(
      diagnostic.solid_nonpositive_zero_component_count == 0U &&
          diagnostic.velocity_layers_bitwise_equal &&
          diagnostic.derived_velocity_dependents_rebuilt &&
          diagnostic.derived_velocity_lineage != 0U &&
          diagnostic.maximum_h_by_a_velocity_difference == 0.0 &&
          diagnostic.velocity_gradient_nonfinite_count == 0U &&
          diagnostic.effective_viscosity_nonpositive_count == 0U &&
          diagnostic.solid_velocity_gradient_nonfinite_count == 0U &&
          diagnostic.velocity_dependent_rate_layers_bitwise_equal,
      "IBM Fresh accepted/previous/trial U are bitwise equal and solid-wall "
      "exact +0 with rebuilt velocity dependents");
  if (!passed)
    print_fresh_diagnostic("IBM Fresh success");
  return passed;
}

bool valid_no_ibm_bypass_diagnostic() {
  detail::FreshInitializationDiagnostic diagnostic;
  bool passed =
      expect(detail::fresh_initialization_diagnostic_for_test(diagnostic),
             "non-IBM Fresh publishes its bypass diagnostic");
  if (!passed)
    return false;
  passed &= expect(
      !diagnostic.immersed && !diagnostic.projection_attempted &&
          diagnostic.no_ibm_bypassed && diagnostic.audited &&
          diagnostic.committed && diagnostic.terminal_status &&
          diagnostic.red_plan != 0U && diagnostic.prepared_lineage != 0U &&
          diagnostic.solved_lineage != 0U &&
          diagnostic.candidate_lineage != 0U && diagnostic.solve.status &&
          diagnostic.solve.termination == LinearTermination::zero_rhs &&
          diagnostic.initial_continuity_maximum == 0.0 &&
          diagnostic.final_continuity_maximum == 0.0 &&
          diagnostic.cut_face_nonzero_count == 0U &&
          diagnostic.solid_nonpositive_zero_component_count == 0U &&
          diagnostic.velocity_layers_bitwise_equal &&
          !diagnostic.derived_velocity_dependents_rebuilt &&
          diagnostic.derived_velocity_lineage == 0U,
      "non-IBM Fresh records a literal bitwise bypass with no projection");
  if (!passed)
    print_fresh_diagnostic("non-IBM Fresh success");
  return passed;
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
  // cylinder_ascii.stl is the closed [-1,1]^3 regression cube.  Cell centres
  // never lie on its surface for this 0.25 grid, so the analytic predicate is
  // identical to the EBTopology region authority.
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

bool fresh_ibm_projects_and_constrains(ProductDriver &driver,
                                       const ValidatedModel &model,
                                       CommittedOutputSnapshot &snapshot) {
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.25, -0.125};
  Status status = driver.initialize(initial);
  if (!status)
    std::cerr << "IBM Fresh initialize status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
  if (!status)
    print_fresh_diagnostic("IBM Fresh");
  bool passed =
      expect(static_cast<bool>(status), "IBM Fresh initialization succeeds");
  if (!status)
    return false;
  status = driver.committed_output_snapshot(snapshot);
  passed &= expect(static_cast<bool>(status) && snapshot.committed &&
                       snapshot.step == 0U,
                   "IBM Fresh publishes a committed step-zero snapshot");
  if (!status)
    return false;

  const SnapshotFieldView *velocity = find_field(snapshot, "U");
  const SnapshotFieldView *pressure = find_field(snapshot, "pi");
  passed &= expect(velocity != nullptr && velocity->values.components == 3U,
                   "IBM Fresh snapshot exposes accepted velocity");
  passed &= expect(pressure != nullptr && pressure->values.components == 1U,
                   "IBM Fresh snapshot exposes pressure perturbation");
  if (velocity == nullptr || pressure == nullptr)
    return false;

  std::uint64_t solid_cells = 0U;
  bool stationary_wall_zero = true;
  bool fluid_velocity_changed = false;
  bool pressure_untouched = true;
  const std::array<double, 3U> base_velocity{
      initial.velocity.x, initial.velocity.y, initial.velocity.z};
  for (std::int32_t z = 0; z < snapshot.patch.cells.z; ++z)
    for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y)
      for (std::int32_t x = 0; x < snapshot.patch.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const bool solid = cube_solid(model, snapshot.patch, cell);
        solid_cells += solid ? 1U : 0U;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double value = velocity->values.unchecked(cell, component);
          if (solid)
            stationary_wall_zero &= same_bits(value, 0.0);
          else
            fluid_velocity_changed |=
                !same_bits(value, base_velocity[component]);
        }
        pressure_untouched &=
            same_bits(pressure->values.unchecked(cell, 0U), 0.0);
      }
  passed &= expect(solid_cells != 0U,
                   "IBM Fresh fixture contains solid control volumes");
  passed &= expect(stationary_wall_zero,
                   "IBM Fresh writes exact stationary-wall +0 in every "
                   "solid accepted velocity component");
  passed &= expect(fluid_velocity_changed,
                   "IBM Fresh compatibility projection changes the fluid "
                   "plug velocity adjacent to the solid");
  passed &= expect(pressure_untouched,
                   "IBM Fresh kinematic projection does not write pressure");
  passed &= valid_ibm_success_diagnostic(snapshot.plan);
  return passed;
}

bool no_ibm_is_bitwise_bypassed() {
  const ValidatedModel model = hundun::v04::test::product_model({8, 7, 6});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {-0.0, 0.375, -0.625};
  if (status)
    status = driver.initialize(initial);
  CommittedOutputSnapshot snapshot;
  if (status)
    status = driver.committed_output_snapshot(snapshot);
  bool passed = expect(static_cast<bool>(status),
                       "non-IBM Fresh initialization succeeds");
  if (!status)
    return false;
  const SnapshotFieldView *velocity = find_field(snapshot, "U");
  passed &= expect(velocity != nullptr,
                   "non-IBM Fresh snapshot exposes accepted velocity");
  if (velocity == nullptr)
    return false;
  const std::array<double, 3U> expected{initial.velocity.x, initial.velocity.y,
                                        initial.velocity.z};
  for (std::int32_t z = 0; z < snapshot.patch.cells.z; ++z)
    for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y)
      for (std::int32_t x = 0; x < snapshot.patch.cells.x; ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component)
          passed &= same_bits(velocity->values.unchecked({x, y, z}, component),
                              expected[component]);
  passed &=
      expect(passed, "non-IBM Fresh is a bitwise accepted-velocity bypass");
  passed &= valid_no_ibm_bypass_diagnostic();
  return passed;
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
  candidate.controller_state = 0U;
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
      for (std::int32_t z = 0; z < patch.cells.z; ++z)
        for (std::int32_t y = 0; y < patch.cells.y; ++y)
          for (std::int32_t x = 0; x < patch.cells.x; ++x) {
            const std::size_t cell =
                static_cast<std::size_t>(x) +
                static_cast<std::size_t>(patch.cells.x) *
                    (static_cast<std::size_t>(y) +
                     static_cast<std::size_t>(patch.cells.y) * z);
            // A non-wall sentinel in solid cells makes any accidental
            // Fresh projection on Restart observable.  Restart owns the
            // checkpoint bytes and must restore them literally.
            const bool solid = cube_solid(model, patch, {x, y, z});
            field.values[cell * 3U] = solid ? kRestartVelocity.x : 0.0625;
            field.values[cell * 3U + 1U] =
                solid ? kRestartVelocity.y : -0.03125;
            field.values[cell * 3U + 2U] =
                solid ? kRestartVelocity.z : 0.015625;
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

bool restart_matches_image(ProductDriver &driver, const RestartImage &image) {
  RestartSnapshot restored;
  Status status = driver.committed_restart_snapshot(restored);
  bool passed = expect(static_cast<bool>(status),
                       "Restart publishes a committed restart snapshot");
  if (!status)
    return false;
  passed &= restored.fields.size == image.fields.size();
  for (std::size_t field_index = 0U;
       field_index < restored.fields.size && passed; ++field_index) {
    const RestartFieldView source = restored.fields.data[field_index];
    const RestartImageField &expected = image.fields[field_index];
    passed &= source.role == expected.role &&
              source.values.components == expected.components;
    std::size_t cell = 0U;
    for (std::int32_t z = 0; z < source.values.interior.z; ++z)
      for (std::int32_t y = 0; y < source.values.interior.y; ++y)
        for (std::int32_t x = 0; x < source.values.interior.x; ++x, ++cell)
          for (std::uint8_t component = 0U;
               component < source.values.components; ++component)
            passed &= same_bits(
                source.values.unchecked({x, y, z}, component),
                expected.values[cell * expected.components + component]);
  }
  const std::array<ConstFaceFieldView, 3U> faces{restored.final_mass_flux.x,
                                                 restored.final_mass_flux.y,
                                                 restored.final_mass_flux.z};
  for (std::size_t axis = 0U; axis < faces.size() && passed; ++axis) {
    std::size_t face = 0U;
    for (std::int32_t z = 0; z < faces[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < faces[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < faces[axis].extents.x; ++x, ++face)
          passed &= same_bits(faces[axis].unchecked({x, y, z}),
                              image.final_mass_flux[axis][face]);
  }
  return expect(passed,
                "Restart skips Fresh projection and restores velocity and "
                "signed-zero final flux bytes exactly");
}

bool restart_skips_fresh_projection(const ValidatedModel &model,
                                    const CommittedOutputSnapshot &fresh) {
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_SELF, model, data_root(), plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  RestartExpected expected;
  if (status)
    status = driver.restart_expected(expected);
  RestartImage image;
  if (status && !make_restart_image(model, fresh, expected, image))
    status = {StatusCode::invalid_plan, 1U};
  detail::clear_fresh_initialization_diagnostic_for_test();
  if (status)
    status = driver.initialize_restart(image);
  bool passed =
      expect(static_cast<bool>(status),
             "valid IBM Restart initializes without Fresh projection");
  if (!status)
    return false;
  detail::FreshInitializationDiagnostic unexpected;
  passed &=
      expect(!detail::fresh_initialization_diagnostic_for_test(unexpected),
             "Restart does not enter or publish the Fresh projection path");
  detail::ColdVelocityDependentsDiagnostic derived;
  passed &= expect(
      detail::cold_velocity_dependents_diagnostic_for_test(derived) &&
          derived.restart && derived.rebuilt &&
          derived.driver_plan == expected.plan && derived.lineage != 0U &&
          derived.maximum_h_by_a_velocity_difference == 0.0 &&
          derived.velocity_gradient_nonfinite_count == 0U &&
          derived.effective_viscosity_nonpositive_count == 0U &&
          derived.solid_velocity_gradient_nonfinite_count == 0U &&
          derived.rate_layers_bitwise_equal,
      "Restart rebuilds velocity dependents without rewriting checkpoint "
      "U/phi");
  return restart_matches_image(driver, image) && passed;
}

bool projection_failure_is_atomic_and_retryable() {
  const ValidatedModel model = immersed_model();
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_SELF, model, data_root(), plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  if (!status)
    return false;

  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.25, -0.125};
  detail::arm_fresh_initialization_candidate_poison_once_for_test(
      0, detail::FreshInitializationPoisonKind::velocity);
  status = driver.initialize(initial);
  CommittedOutputSnapshot unpublished;
  const Status output_status = driver.committed_output_snapshot(unpublished);
  detail::FreshInitializationDiagnostic rejected;
  bool passed = expect(
      !status && !driver.initialized() && !output_status &&
          driver.pressure_reference() == 0.0 &&
          driver.closed_mass_target() == 0.0 &&
          detail::fresh_initialization_diagnostic_for_test(rejected) &&
          rejected.immersed && rejected.projection_attempted &&
          rejected.audited && !rejected.committed &&
          !rejected.terminal_status && rejected.prepared_lineage != 0U &&
          rejected.solved_lineage != 0U && rejected.candidate_lineage != 0U &&
          rejected.velocity_layers_bitwise_equal,
      "audit-after candidate poison rejects before any Fresh publication");
  if (!passed)
    print_fresh_diagnostic("IBM Fresh poisoned");

  status = driver.initialize(initial);
  detail::FreshInitializationDiagnostic recovered;
  passed &= expect(
      status && driver.initialized() &&
          detail::fresh_initialization_diagnostic_for_test(recovered) &&
          recovered.committed && recovered.terminal_status &&
          recovered.generation > rejected.generation &&
          recovered.derived_velocity_dependents_rebuilt &&
          recovered.derived_velocity_lineage != 0U &&
          recovered.velocity_dependent_rate_layers_bitwise_equal,
      "poison authority is consumed and the same valid Fresh input retries "
      "successfully");
  return passed;
}

bool failed_initialize_is_retryable() {
  const ValidatedModel model = hundun::v04::test::product_model({8, 7, 6});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  if (!status)
    return false;

  DriverInitialState invalid;
  invalid.pressure_reference = 98000.0;
  invalid.temperature = 315.0;
  invalid.velocity.x = std::numeric_limits<double>::quiet_NaN();
  status = driver.initialize(invalid);
  CommittedOutputSnapshot unpublished;
  const Status output_status = driver.committed_output_snapshot(unpublished);
  bool passed =
      expect(!status && status.code == StatusCode::invalid_plan &&
                 !driver.initialized() && driver.pressure_reference() == 0.0 &&
                 driver.closed_mass_target() == 0.0 && !output_status,
             "failed Fresh initialization publishes no committed driver state");

  DriverInitialState valid;
  valid.pressure_reference = 98000.0;
  valid.temperature = 315.0;
  valid.velocity = {0.125, -0.25, 0.5};
  status = driver.initialize(valid);
  passed &= expect(static_cast<bool>(status) && driver.initialized(),
                   "failed Fresh initialization can be retried successfully");
  return passed;
}

bool run() {
  const ValidatedModel model = immersed_model();
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_SELF, model, data_root(), plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  bool passed = expect(static_cast<bool>(status),
                       "IBM Fresh ProductDriver fixture compiles");
  if (!status)
    return false;

  CommittedOutputSnapshot fresh;
  passed &= fresh_ibm_projects_and_constrains(driver, model, fresh);
  passed &= no_ibm_is_bitwise_bypassed();
  passed &= restart_skips_fresh_projection(model, fresh);
  passed &= projection_failure_is_atomic_and_retryable();
  passed &= failed_initialize_is_retryable();
  return passed;
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  const bool passed = run();
  MPI_Finalize();
  return passed ? 0 : 1;
}
