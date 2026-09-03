// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"

#include "../support/product_fixture.hpp"
#include "core_product_freeze_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr Int3 kGlobalCells{8, 8, 8};
constexpr double kCheckpointTime = 0.006;
constexpr double kFullDt = 0.009;
constexpr double kHalfDt = 0.0045;
constexpr double kReferencePressure = 101325.0;
constexpr double kMolecularWeight = 28.96546;
constexpr double kGasConstant = kUniversalGasConstant / kMolecularWeight;
constexpr double kMinimumEnthalpy = 3.5 * kGasConstant * 200.0;
constexpr double kBaseEnthalpy = kMinimumEnthalpy + 1.0;
constexpr double kLowEnthalpy = kMinimumEnthalpy + 1.0e-10;
constexpr double kHighEnthalpy = kMinimumEnthalpy + 2.0;
constexpr double kStreamwiseVelocity = 3.0e-3;
constexpr double kConductivity = 1.625e-5;

std::uint64_t wire_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool collective(bool local) noexcept {
  const int value = local ? 1 : 0;
  int global = 0;
  return MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ==
             MPI_SUCCESS &&
         global != 0;
}

bool same_u64(std::uint64_t value) noexcept {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

bool expect(bool condition, int rank, const char *description) {
  if (!condition)
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  return condition;
}

std::size_t cell_offset(Int3 cells, Int3 local) noexcept {
  return static_cast<std::size_t>(local.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(local.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(local.z));
}

bool owns_cell(const MeshPatch &patch, Int3 global, Int3 &local) noexcept {
  local = {global.x - patch.begin.x, global.y - patch.begin.y,
           global.z - patch.begin.z};
  return local.x >= 0 && local.y >= 0 && local.z >= 0 &&
         local.x < patch.cells.x && local.y < patch.cells.y &&
         local.z < patch.cells.z;
}

ValidatedModel retry_model(double minimum_dt, double maximum_dt,
                           std::uint32_t maximum_retries,
                           PlanFingerprint fingerprint) {
  ValidatedModel model = test::product_model(kGlobalCells);
  model.fingerprint = fingerprint;
  model.solver.pressure.krylov_restart = 64U;
  model.time.initial_dt = maximum_dt;
  model.time.minimum_dt = minimum_dt;
  model.time.maximum_dt = maximum_dt;
  model.time.maximum_growth = 1.0;
  model.time.retry_factor = 0.5;
  model.time.maximum_retries = maximum_retries;
  model.thermophysics.species[0U].conductivity = kConductivity;
  return model;
}

RestartImage make_restart_image(const RestartExpected &expected,
                                double checkpoint_dt = kFullDt) {
  RestartImage image;
  image.global_cells = expected.global_cells;
  image.patch = expected.target_patch;
  image.plan = expected.plan;
  image.schema = expected.schema;
  image.geometry = expected.geometry;
  image.time = kCheckpointTime;
  image.dt = checkpoint_dt;
  image.pressure_reference = kReferencePressure;
  image.step = 1U;
  image.controller_state = 1U;
  image.backward_euler_recovery = true;

  const std::size_t cells = static_cast<std::size_t>(image.patch.cells.x) *
                            static_cast<std::size_t>(image.patch.cells.y) *
                            static_cast<std::size_t>(image.patch.cells.z);
  image.fields.reserve(expected.fields.size);
  for (std::size_t index = 0U; index < expected.fields.size; ++index) {
    const RestartExpectedField descriptor = expected.fields.data[index];
    RestartImageField field;
    field.role = descriptor.role;
    field.field = descriptor.field;
    field.components = descriptor.components;
    field.values.assign(cells * descriptor.components, 0.0);
    if (descriptor.role == RestartFieldRole::enthalpy)
      std::fill(field.values.begin(), field.values.end(), kBaseEnthalpy);
    if (descriptor.role == RestartFieldRole::velocity)
      for (std::size_t cell = 0U; cell < cells; ++cell)
        field.values[cell * 3U] = kStreamwiseVelocity;
    image.fields.push_back(std::move(field));
  }

  // These perturbations are global physical cells, not rank-local offsets.
  // The fixture is therefore identical under the 1/2/4-rank decompositions.
  for (RestartImageField &field : image.fields) {
    if (field.role != RestartFieldRole::enthalpy)
      continue;
    Int3 local{};
    if (owns_cell(image.patch, {0, 0, 0}, local))
      field.values[cell_offset(image.patch.cells, local)] = kLowEnthalpy;
    if (owns_cell(image.patch, {1, 0, 0}, local))
      field.values[cell_offset(image.patch.cells, local)] = kHighEnthalpy;
  }

  const Int3 local = image.patch.cells;
  image.final_mass_flux[0U].assign(
      static_cast<std::size_t>(local.x + 1) * local.y * local.z, 0.0);
  image.final_mass_flux[1U].assign(
      static_cast<std::size_t>(local.x) * (local.y + 1) * local.z, 0.0);
  image.final_mass_flux[2U].assign(
      static_cast<std::size_t>(local.x) * local.y * (local.z + 1), 0.0);
  // The periodic mass flux is exactly divergence-free and consistent on
  // every decomposition.  Thus the rejection below is a pressure--enthalpy
  // positivity/terminal issue, not an incompatible-flux cold start.
  constexpr double density = kReferencePressure / (kGasConstant * 200.0);
  constexpr double x_face_area = (1.0 / 8.0) * (0.5 / 8.0);
  std::fill(image.final_mass_flux[0U].begin(), image.final_mass_flux[0U].end(),
            density * kStreamwiseVelocity * x_face_area);
  return image;
}

struct DriverHarness {
  ValidatedModel model{};
  ProductDriver driver{};
  Status status{};
};

DriverHarness make_driver(ValidatedModel model,
                          double checkpoint_dt = kFullDt) {
  DriverHarness harness;
  harness.model = std::move(model);
  CompiledCasePlan plan;
  harness.status =
      ProductCompiler::compile(MPI_COMM_WORLD, harness.model, {}, plan);
  if (harness.status)
    harness.status =
        ProductDriver::create(MPI_COMM_WORLD, std::move(plan), harness.driver);
  RestartExpected expected;
  if (harness.status)
    harness.status = harness.driver.restart_expected(expected);
  if (harness.status) {
    RestartImage image = make_restart_image(expected, checkpoint_dt);
    harness.status = harness.driver.initialize_restart(image);
  }
  return harness;
}

struct ExactCommittedBits {
  std::uint64_t time{};
  std::uint64_t dt{};
  std::uint64_t pressure_reference{};
  std::uint64_t step{};
  std::uint64_t controller_state{};
  std::vector<std::uint64_t> field_metadata;
  std::vector<std::uint64_t> field_values;
  std::array<std::vector<std::uint64_t>, 3U> flux_values;
  RevisionToken flux_revision{};
  FaceFluxCertificate flux_certificate{};
};

bool capture_exact_committed(ProductDriver &driver, ExactCommittedBits &out) {
  RestartSnapshot snapshot;
  const Status status = driver.committed_restart_snapshot(snapshot);
  if (!status || !snapshot.final_mass_flux.certificate.valid() ||
      !snapshot.final_mass_flux.certificate.matches(snapshot.final_mass_flux))
    return false;
  ExactCommittedBits candidate;
  candidate.time = wire_bits(snapshot.time);
  candidate.dt = wire_bits(snapshot.dt);
  candidate.pressure_reference = wire_bits(snapshot.pressure_reference);
  candidate.step = snapshot.step;
  candidate.controller_state = snapshot.controller_state;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView &field = snapshot.fields.data[index];
    candidate.field_metadata.push_back(static_cast<std::uint64_t>(field.role));
    candidate.field_metadata.push_back(field.values.field);
    candidate.field_metadata.push_back(field.values.revision);
    candidate.field_metadata.push_back(field.values.storage_identity);
    candidate.field_metadata.push_back(field.values.revision_domain);
    for (std::uint8_t component = 0U; component < field.values.components;
         ++component)
      for (std::int32_t z = 0; z < field.values.interior.z; ++z)
        for (std::int32_t y = 0; y < field.values.interior.y; ++y)
          for (std::int32_t x = 0; x < field.values.interior.x; ++x)
            candidate.field_values.push_back(
                wire_bits(field.values.unchecked({x, y, z}, component)));
  }
  const std::array<ConstFaceFieldView, 3U> flux{{snapshot.final_mass_flux.x,
                                                 snapshot.final_mass_flux.y,
                                                 snapshot.final_mass_flux.z}};
  for (std::size_t axis = 0U; axis < flux.size(); ++axis)
    for (std::int32_t z = 0; z < flux[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < flux[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < flux[axis].extents.x; ++x)
          candidate.flux_values[axis].push_back(
              wire_bits(flux[axis].unchecked({x, y, z})));
  candidate.flux_revision = snapshot.final_mass_flux.revision;
  candidate.flux_certificate = snapshot.final_mass_flux.certificate;
  out = std::move(candidate);
  return true;
}

bool exact_committed_equal(const ExactCommittedBits &left,
                           const ExactCommittedBits &right) noexcept {
  return left.time == right.time && left.dt == right.dt &&
         left.pressure_reference == right.pressure_reference &&
         left.step == right.step &&
         left.controller_state == right.controller_state &&
         left.field_metadata == right.field_metadata &&
         left.field_values == right.field_values &&
         left.flux_values == right.flux_values &&
         left.flux_revision == right.flux_revision &&
         left.flux_certificate == right.flux_certificate;
}

struct PhysicalCommittedBits {
  std::uint64_t time{};
  std::uint64_t dt{};
  std::uint64_t pressure_reference{};
  std::uint64_t step{};
  std::vector<std::string> field_names;
  std::vector<std::uint64_t> field_components;
  std::vector<std::uint64_t> field_values;
  std::array<std::vector<std::uint64_t>, 3U> flux_values;
};

bool capture_physical_committed(ProductDriver &driver,
                                PhysicalCommittedBits &out) {
  CommittedOutputSnapshot output;
  if (!driver.committed_output_snapshot(output) || !output.committed)
    return false;
  RestartSnapshot snapshot;
  if (!driver.committed_restart_snapshot(snapshot) ||
      !snapshot.final_mass_flux.certificate.valid() ||
      !snapshot.final_mass_flux.certificate.matches(snapshot.final_mass_flux) ||
      wire_bits(output.time) != wire_bits(snapshot.time) ||
      output.step != snapshot.step)
    return false;
  PhysicalCommittedBits candidate;
  candidate.time = wire_bits(output.time);
  candidate.dt = wire_bits(snapshot.dt);
  candidate.pressure_reference = wire_bits(snapshot.pressure_reference);
  candidate.step = output.step;
  for (std::size_t index = 0U; index < output.fields.size; ++index) {
    const SnapshotFieldView &field = output.fields.data[index];
    candidate.field_names.emplace_back(field.stable_name);
    candidate.field_components.push_back(field.values.components);
    for (std::uint8_t component = 0U; component < field.values.components;
         ++component)
      for (std::int32_t z = 0; z < field.values.interior.z; ++z)
        for (std::int32_t y = 0; y < field.values.interior.y; ++y)
          for (std::int32_t x = 0; x < field.values.interior.x; ++x)
            candidate.field_values.push_back(
                wire_bits(field.values.unchecked({x, y, z}, component)));
  }
  const std::array<ConstFaceFieldView, 3U> flux{{snapshot.final_mass_flux.x,
                                                 snapshot.final_mass_flux.y,
                                                 snapshot.final_mass_flux.z}};
  for (std::size_t axis = 0U; axis < flux.size(); ++axis)
    for (std::int32_t z = 0; z < flux[axis].extents.z; ++z)
      for (std::int32_t y = 0; y < flux[axis].extents.y; ++y)
        for (std::int32_t x = 0; x < flux[axis].extents.x; ++x)
          candidate.flux_values[axis].push_back(
              wire_bits(flux[axis].unchecked({x, y, z})));
  out = std::move(candidate);
  return true;
}

bool physical_committed_equal(const PhysicalCommittedBits &left,
                              const PhysicalCommittedBits &right) noexcept {
  return left.time == right.time && left.dt == right.dt &&
         left.pressure_reference == right.pressure_reference &&
         left.step == right.step && left.field_names == right.field_names &&
         left.field_components == right.field_components &&
         left.field_values == right.field_values &&
         left.flux_values == right.flux_values;
}

bool finite_positive_terminal_state(ProductDriver &driver,
                                    const ValidatedModel &model,
                                    const DriverStepReport &report) {
  RestartSnapshot snapshot;
  if (!driver.committed_restart_snapshot(snapshot))
    return false;
  const ConstFieldView *pressure = nullptr;
  const ConstFieldView *enthalpy = nullptr;
  const ConstFieldView *velocity = nullptr;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView &field = snapshot.fields.data[index];
    if (field.role == RestartFieldRole::pressure_perturbation)
      pressure = &field.values;
    else if (field.role == RestartFieldRole::enthalpy)
      enthalpy = &field.values;
    else if (field.role == RestartFieldRole::velocity)
      velocity = &field.values;
  }
  if (pressure == nullptr || enthalpy == nullptr || velocity == nullptr)
    return false;
  // For this one-species NASA7 fixture h=3.5*R*T exactly.  Reconstructing the
  // EOS from the public p/h snapshot independently certifies positive T/rho;
  // the ProductDriver report supplies the production EOS residual gate.
  bool valid = true;
  for (std::int32_t z = 0; z < pressure->interior.z; ++z)
    for (std::int32_t y = 0; y < pressure->interior.y; ++y)
      for (std::int32_t x = 0; x < pressure->interior.x; ++x) {
        const Int3 cell{x, y, z};
        const double absolute_pressure =
            snapshot.pressure_reference + pressure->unchecked(cell, 0U);
        const double h = enthalpy->unchecked(cell, 0U);
        const double temperature = h / (3.5 * kGasConstant);
        const double rho = absolute_pressure / (kGasConstant * temperature);
        valid &= std::isfinite(absolute_pressure) && absolute_pressure > 0.0 &&
                 std::isfinite(h) && h > 0.0 && std::isfinite(temperature) &&
                 temperature >= model.thermophysics.minimum_temperature &&
                 temperature <= model.thermophysics.maximum_temperature &&
                 std::isfinite(rho) && rho > 0.0;
        for (std::uint8_t component = 0U; component < 3U; ++component)
          valid &= std::isfinite(velocity->unchecked(cell, component));
      }
  const std::array<ConstFaceFieldView, 3U> flux{{snapshot.final_mass_flux.x,
                                                 snapshot.final_mass_flux.y,
                                                 snapshot.final_mass_flux.z}};
  for (const ConstFaceFieldView &axis : flux)
    for (std::int32_t z = 0; z < axis.extents.z; ++z)
      for (std::int32_t y = 0; y < axis.extents.y; ++y)
        for (std::int32_t x = 0; x < axis.extents.x; ++x)
          valid &= std::isfinite(axis.unchecked({x, y, z}));
  valid &=
      std::isfinite(report.piso.eos_residual) &&
      report.piso.eos_residual <= model.solver.terminal.eos &&
      std::isfinite(report.piso.continuity_residual) &&
      report.piso.continuity_residual <= model.solver.terminal.continuity &&
      std::isfinite(report.piso.energy_residual) &&
      report.piso.energy_residual <= model.solver.terminal.continuity &&
      std::isfinite(report.piso.closed_mass_residual) &&
      report.piso.closed_mass_residual <= model.solver.terminal.closed_mass &&
      std::isfinite(report.piso.gauge_residual) &&
      report.piso.gauge_residual <= model.solver.terminal.gauge;
  return valid;
}

bool refinement_prefix_certificate(
    const DriverStepReport &report, std::uint8_t expected_calls,
    PressureEnergyRefinementTermination expected_termination) {
  const PisoAttemptReport &piso = report.piso;
  bool valid =
      PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(piso) &&
      piso.pressure_energy_refinement_solve_calls == expected_calls &&
      piso.pressure_energy_refinement_termination == expected_termination;
  RevisionToken target_generation = 0U;
  PlanFingerprint prior_collective_lineage = 0U;
  RevisionToken prior_pressure_state = 0U;
  RevisionToken prior_numeric = 0U;
  LinearIdentity static_identity{};
  for (std::uint8_t index = 0U; index < expected_calls; ++index) {
    const PisoPressureEnergyRefinementSolveReport &refinement =
        piso.pressure_energy_refinement[index];
    const bool accepted_solve =
        refinement.solve.status &&
        (refinement.solve.termination == LinearTermination::converged ||
         refinement.solve.termination == LinearTermination::zero_rhs) &&
        std::isfinite(refinement.solve.initial_true_residual) &&
        std::isfinite(refinement.solve.final_true_residual) &&
        refinement.solve.final_true_residual <=
            refinement.solve.initial_true_residual;
    const bool no_dead_recycle_capture =
        refinement.solve.recycle_cycle_corrections == 0U &&
        refinement.solve.recycle_capture_vector_passes == 0U &&
        refinement.solve.recycle_capture_cycle_attempts == 0U &&
        refinement.solve.recycle_capture_reduction_calls == 0U &&
        refinement.solve.recycle_capture_blocking_operations == 0U;
    valid &= refinement.valid() && accepted_solve &&
             no_dead_recycle_capture &&
             refinement.ordinal == static_cast<std::uint8_t>(index + 1U) &&
             same_u64(refinement.target_generation) &&
             same_u64(refinement.collective_lineage) &&
             refinement.linear_identity.numeric == refinement.pressure_state;
    for (std::uint8_t prior = 0U; prior < index; ++prior)
      valid &= refinement.collective_lineage !=
               piso.pressure_energy_refinement[prior].collective_lineage;
    if (index == 0U) {
      target_generation = refinement.target_generation;
      static_identity = refinement.linear_identity;
    } else {
      // A refinement stays at the same physical target but must be assembled
      // from the newly accepted provisional state and a fresh numeric
      // linearization on every round.
      valid &=
          refinement.target_generation == target_generation &&
          refinement.collective_lineage != prior_collective_lineage &&
          refinement.pressure_state != prior_pressure_state &&
          refinement.linear_identity.numeric != prior_numeric &&
          refinement.linear_identity.fingerprint !=
              piso.pressure_energy_refinement[index - 1U]
                  .linear_identity.fingerprint &&
          refinement.linear_identity.symbolic == static_identity.symbolic &&
          refinement.linear_identity.hierarchy == static_identity.hierarchy &&
          refinement.linear_identity.workspace == static_identity.workspace;
    }
    prior_collective_lineage = refinement.collective_lineage;
    prior_pressure_state = refinement.pressure_state;
    prior_numeric = refinement.linear_identity.numeric;
  }
  for (std::size_t index = expected_calls;
       index < piso.pressure_energy_refinement.size(); ++index)
    valid &= !piso.pressure_energy_refinement[index].valid();
  return valid;
}

bool refinement_trajectory_certificate(const DriverStepReport &report,
                                       std::uint8_t expected_calls) {
  const PressureEnergyGlobalizationAttemptReport &globalization =
      report.pressure_energy_globalization;
  const std::uint8_t expected_trajectory =
      static_cast<std::uint8_t>(2U + expected_calls);
  bool valid = globalization.valid &&
               globalization.trajectory_count == expected_trajectory;
  RevisionToken target_time = 0U;
  const auto mix = [](std::uint64_t hash, std::uint64_t value) {
    hash ^= value;
    return hash * UINT64_C(1099511628211);
  };
  const auto mix_sample = [&](std::uint64_t hash,
                              const PressureEnergyGlobalizationSample &sample) {
    hash = mix(hash, wire_bits(sample.alpha));
    hash = mix(hash, wire_bits(sample.global_normalized_continuity));
    hash = mix(hash, wire_bits(sample.global_normalized_energy));
    hash = mix(hash, sample.thermodynamically_admissible);
    hash = mix(hash, sample.state_and_flux_finite);
    hash = mix(hash, sample.corrector);
    hash = mix(hash, sample.target_time);
    hash = mix(hash, sample.correction_direction);
    hash = mix(hash, sample.state_provenance);
    return mix(hash, sample.mass_flux_provenance);
  };
  std::uint64_t trajectory_hash = UINT64_C(1469598103934665603);
  for (std::uint8_t index = 0U; index < expected_trajectory; ++index) {
    const PressureEnergyGlobalizationIterationReport &iteration =
        globalization.trajectory[index];
    const std::uint8_t expected_corrector = index == 0U ? 1U : 2U;
    const std::uint8_t expected_refinement =
        index < 2U ? 0U : static_cast<std::uint8_t>(index - 1U);
    const double baseline_merit =
        std::hypot(iteration.baseline.global_normalized_continuity,
                   iteration.baseline.global_normalized_energy);
    const double selected_merit =
        std::hypot(iteration.selected.global_normalized_continuity,
                   iteration.selected.global_normalized_energy);
    valid &= iteration.valid && iteration.corrector == expected_corrector &&
             iteration.refinement_iteration == expected_refinement &&
             iteration.baseline.alpha == 0.0 &&
             iteration.selected.alpha > 0.0 &&
             iteration.selected.alpha <= 1.0 &&
             iteration.baseline.thermodynamically_admissible &&
             iteration.baseline.state_and_flux_finite &&
             iteration.selected.thermodynamically_admissible &&
             iteration.selected.state_and_flux_finite &&
             std::isfinite(baseline_merit) && std::isfinite(selected_merit) &&
             selected_merit < baseline_merit &&
             iteration.baseline.correction_direction != 0U &&
             iteration.baseline.state_provenance != 0U &&
             iteration.baseline.mass_flux_provenance != 0U &&
             iteration.selected.correction_direction ==
                 iteration.baseline.correction_direction &&
             iteration.selected.state_provenance !=
                 iteration.baseline.state_provenance &&
             iteration.selected.mass_flux_provenance !=
                 iteration.baseline.mass_flux_provenance;
    if (index == 0U)
      target_time = iteration.baseline.target_time;
    valid &= iteration.baseline.corrector == expected_corrector &&
             iteration.selected.corrector == expected_corrector &&
             iteration.baseline.target_time == target_time &&
             iteration.selected.target_time == target_time;
    if (index >= 2U) {
      const PressureEnergyGlobalizationIterationReport &prior =
          globalization.trajectory[index - 1U];
      valid &= iteration.baseline.correction_direction !=
                   prior.baseline.correction_direction &&
               iteration.baseline.state_provenance !=
                   prior.baseline.state_provenance &&
               iteration.baseline.mass_flux_provenance !=
                   prior.baseline.mass_flux_provenance &&
               wire_bits(iteration.baseline.global_normalized_continuity) ==
                   wire_bits(prior.selected.global_normalized_continuity) &&
               wire_bits(iteration.baseline.global_normalized_energy) ==
                   wire_bits(prior.selected.global_normalized_energy);
    }
    trajectory_hash = mix(trajectory_hash, iteration.valid);
    trajectory_hash = mix(trajectory_hash, iteration.corrector);
    trajectory_hash = mix(trajectory_hash, iteration.refinement_iteration);
    trajectory_hash =
        mix(trajectory_hash,
            wire_bits(iteration.maximum_absolute_pressure_correction));
    trajectory_hash =
        mix(trajectory_hash,
            wire_bits(iteration.maximum_absolute_enthalpy_correction));
    trajectory_hash = mix_sample(trajectory_hash, iteration.baseline);
    trajectory_hash = mix_sample(trajectory_hash, iteration.selected);
  }
  for (std::size_t index = expected_trajectory;
       index < globalization.trajectory.size(); ++index)
    valid &= !globalization.trajectory[index].valid;
  return valid && same_u64(trajectory_hash);
}

std::uint64_t diagnostic_hash(
    const detail::PressureEnergyCandidateGlobalizationDiagnostic &value) {
  const auto mix = [](std::uint64_t hash, std::uint64_t item) {
    hash ^= item;
    return hash * UINT64_C(1099511628211);
  };
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, value.valid);
  hash = mix(hash, value.production_candidate_loop);
  hash = mix(hash, value.selection_valid);
  hash = mix(hash, value.replay_valid);
  hash = mix(hash, value.committed);
  hash = mix(hash, value.attempted_step);
  hash = mix(hash, value.generation);
  hash = mix(hash, value.attempt);
  hash = mix(hash, value.corrector);
  hash = mix(hash, value.sample_count);
  hash = mix(hash, value.first_admissible_sample);
  hash = mix(hash, wire_bits(value.maximum_absolute_pressure_correction));
  hash = mix(hash, wire_bits(value.maximum_absolute_enthalpy_correction));
  const std::size_t count = std::min<std::size_t>(
      value.sample_count,
      detail::kPressureEnergyCandidateGlobalizationSampleCapacity);
  for (std::size_t index = 0U; index < count; ++index) {
    const auto &sample = value.samples[index];
    hash = mix(hash, wire_bits(sample.alpha));
    hash = mix(hash, sample.admissible);
    hash = mix(hash, sample.state_and_flux_finite);
    hash = mix(hash, sample.first_failing_global_cell);
    hash = mix(hash, static_cast<std::uint8_t>(sample.first_failure_reason));
    hash = mix(hash, wire_bits(sample.minimum_absolute_pressure));
    hash = mix(hash, wire_bits(sample.minimum_temperature));
    hash = mix(hash, wire_bits(sample.minimum_density));
    hash = mix(hash, wire_bits(sample.normalized_continuity));
    hash = mix(hash, wire_bits(sample.normalized_energy));
  }
  return hash;
}

struct ProbeResult {
  bool passed{};
  Status status{};
  DriverStepReport report{};
  detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic{};
};

ProbeResult run_full_probe(int rank) {
  DriverHarness probe =
      make_driver(retry_model(kFullDt, kFullDt, 1U, UINT64_C(0x18000c101)));
  ExactCommittedBits before;
  const bool captured_before =
      probe.status && capture_exact_committed(probe.driver, before);
  if (probe.status)
    detail::arm_pressure_energy_candidate_globalization_once_for_test();
  DriverStepReport report;
  if (probe.status)
    probe.status = probe.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, report);
  detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic;
  const bool observed =
      detail::pressure_energy_candidate_globalization_diagnostic_for_test(
          diagnostic);
  detail::clear_pressure_energy_candidate_globalization_for_test();
  ExactCommittedBits after;
  const bool captured_after = capture_exact_committed(probe.driver, after);

  const std::size_t count = std::min<std::size_t>(
      diagnostic.sample_count,
      detail::kPressureEnergyCandidateGlobalizationSampleCapacity);
  bool alpha_sequence = count == diagnostic.sample_count;
  bool smaller_admissible = false;
  for (std::size_t index = 0U; index < count; ++index) {
    alpha_sequence &= diagnostic.samples[index].alpha ==
                      std::ldexp(1.0, -static_cast<int>(index));
    if (index != 0U)
      smaller_admissible |= diagnostic.samples[index].admissible;
  }
  const bool first_admissible_valid =
      diagnostic.first_admissible_sample > 0U &&
      diagnostic.first_admissible_sample < count;
  bool first_admissible_is_first = first_admissible_valid;
  for (std::size_t index = 0U;
       index < diagnostic.first_admissible_sample && index < count; ++index)
    first_admissible_is_first &= !diagnostic.samples[index].admissible;
  if (first_admissible_valid)
    first_admissible_is_first &=
        diagnostic.samples[diagnostic.first_admissible_sample].admissible &&
        diagnostic.samples[diagnostic.first_admissible_sample]
            .state_and_flux_finite;
  const bool full_invalid = count != 0U && !diagnostic.samples[0U].admissible &&
                            !diagnostic.samples[0U].state_and_flux_finite &&
                            diagnostic.samples[0U].first_failure_reason ==
                                detail::PressureEnergyCandidateFailureReason::
                                    production_candidate_evaluation;
  const bool rejected = !probe.status && !report.accepted &&
                        probe.status.code == StatusCode::rejected_step &&
                        probe.status.detail == 10210U &&
                        report.failed_stage == 54U && report.attempts == 1U &&
                        report.piso.final_flux_revision == 0U &&
                        !report.piso.continuity_witness.valid;
  // A fatal proposal retires its one-shot controller generation while every
  // committed physical value and revision rolls back exactly.  Normalize the
  // expected ticket transition before comparing the committed snapshot.
  const bool fatal_controller_recovery =
      captured_before && captured_after &&
      before.controller_state != std::numeric_limits<std::uint64_t>::max() &&
      after.controller_state == before.controller_state + 1U;
  if (captured_before && captured_after)
    before.controller_state = after.controller_state;
  const bool rollback = captured_before && captured_after &&
                        fatal_controller_recovery &&
                        exact_committed_equal(before, after);
  const bool diagnostic_identical = same_u64(diagnostic_hash(diagnostic));
  const bool refinement_exhausted =
      refinement_prefix_certificate(
          report, static_cast<std::uint8_t>(kPressureEnergyRefinementCapacity),
          PressureEnergyRefinementTermination::iteration_capacity_exhausted) &&
      refinement_trajectory_certificate(
          report, static_cast<std::uint8_t>(kPressureEnergyRefinementCapacity));
  const bool local =
      rejected && rollback && observed && diagnostic.valid &&
      diagnostic.production_candidate_loop && diagnostic.corrector == 2U &&
      diagnostic.selection_valid && diagnostic.replay_valid &&
      !diagnostic.committed &&
      diagnostic.sample_count > diagnostic.first_admissible_sample &&
      diagnostic.sample_count <=
          detail::kPressureEnergyCandidateGlobalizationSampleCapacity &&
      alpha_sequence && full_invalid && smaller_admissible &&
      first_admissible_valid && first_admissible_is_first &&
      std::isfinite(diagnostic.maximum_absolute_pressure_correction) &&
      diagnostic.maximum_absolute_pressure_correction > 0.0 &&
      std::isfinite(diagnostic.maximum_absolute_enthalpy_correction) &&
      diagnostic.maximum_absolute_enthalpy_correction > 0.0 &&
      diagnostic.selected_alpha ==
          diagnostic.samples[diagnostic.sample_count - 1U].alpha &&
      std::isfinite(diagnostic.selected_normalized_continuity) &&
      diagnostic.selected_normalized_continuity <=
          probe.model.solver.terminal.continuity &&
      std::isfinite(diagnostic.selected_normalized_energy) &&
      diagnostic.selected_normalized_energy >
          probe.model.solver.terminal.continuity &&
      diagnostic_identical && refinement_exhausted;
  if (rank == 0) {
    std::cout << std::setprecision(17)
              << "full-probe U/k/h=" << kStreamwiseVelocity << '/'
              << kConductivity << '/' << kLowEnthalpy << ':' << kBaseEnthalpy
              << ':' << kHighEnthalpy
              << " status=" << static_cast<unsigned>(probe.status.code) << '/'
              << probe.status.detail << " stage=" << report.failed_stage
              << " attempts/accepted=" << report.attempts << '/'
              << report.accepted << " diag-corrector/select/commit="
              << static_cast<unsigned>(diagnostic.corrector) << '/'
              << diagnostic.selection_valid << '/' << diagnostic.committed
              << " first-admissible="
              << static_cast<unsigned>(diagnostic.first_admissible_sample)
              << " dp/dh=" << diagnostic.maximum_absolute_pressure_correction
              << '/' << diagnostic.maximum_absolute_enthalpy_correction
              << " terminal=" << report.piso.eos_residual << '/'
              << report.piso.continuity_residual << '/'
              << report.piso.energy_residual << '/'
              << report.piso.closed_mass_residual << '/'
              << report.piso.gauge_residual << " refinement="
              << static_cast<unsigned>(
                     report.piso.pressure_energy_refinement_solve_calls)
              << '/'
              << static_cast<unsigned>(
                     report.piso.pressure_energy_refinement_termination)
              << " rollback=" << rollback << '\n';
    if (!local)
      for (std::size_t index = 0U; index < count; ++index) {
        const auto &sample = diagnostic.samples[index];
        std::cout << "  alpha[" << index << "]=" << sample.alpha
                  << " admissible/finite=" << sample.admissible << '/'
                  << sample.state_and_flux_finite << " reason="
                  << static_cast<unsigned>(sample.first_failure_reason)
                  << " min=" << sample.minimum_absolute_pressure << '/'
                  << sample.minimum_temperature << '/' << sample.minimum_density
                  << " RC/RE=" << sample.normalized_continuity << '/'
                  << sample.normalized_energy << '\n';
      }
  }
  return {collective(local), probe.status, report, diagnostic};
}

bool run_refinement_certificate(int rank) {
  const bool replay_roundoff_contract =
      detail::alpha_zero_energy_replay_equivalent(false, true, 128.0 *
              std::numeric_limits<double>::epsilon()) &&
      !detail::alpha_zero_energy_replay_equivalent(false, true, 1.0e-7) &&
      !detail::alpha_zero_energy_replay_equivalent(
          false, false, std::numeric_limits<double>::denorm_min()) &&
      detail::alpha_zero_energy_replay_equivalent(true, false, 1.0);
  // The first six refreshed solves still form a strict descent sequence but
  // do not cross this gate.  The twelfth does.  This is the focused RED for
  // the production failure in which a useful same-target direction was
  // truncated by the original six-entry hot-resource contract.
  constexpr double kRefinementFixtureContinuityGate = 1.0460408e-10;
  ValidatedModel model =
      retry_model(kFullDt, kFullDt, 1U, UINT64_C(0x18000c401));
  // The ordinary C2 replay is just outside this fixed component gate.  The
  // first eleven refreshed solves remain outside and the twelfth crosses it,
  // giving a deterministic accepted product step with bounded extra headroom.
  model.solver.terminal.continuity = kRefinementFixtureContinuityGate;
  DriverHarness refined = make_driver(std::move(model));
  DriverStepReport report;
  if (refined.status)
    refined.status = refined.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, report);

  constexpr std::uint8_t expected_calls = 12U;
  constexpr std::size_t expected_terminal_index = expected_calls + 1U;
  const bool expected_prefix_available =
      report.piso.pressure_energy_refinement_solve_calls == expected_calls &&
      report.pressure_energy_globalization.trajectory_count ==
          2U + expected_calls;
  const bool prefix =
      expected_prefix_available &&
      refinement_prefix_certificate(
          report, expected_calls,
          PressureEnergyRefinementTermination::component_residuals_converged);
  const bool trajectory = expected_prefix_available &&
                          refinement_trajectory_certificate(report,
                                                            expected_calls);
  const auto &path = report.pressure_energy_globalization.trajectory;
  const double component_tolerance = refined.model.solver.terminal.continuity;
  const bool first_crosses_component_gate =
      expected_prefix_available &&
      path[expected_terminal_index - 1U]
              .selected.global_normalized_energy > component_tolerance &&
      path[expected_terminal_index].selected.global_normalized_energy <=
          component_tolerance &&
      path[expected_terminal_index].selected.global_normalized_continuity <=
          component_tolerance;
  const bool terminal_matches_selected_trajectory =
      expected_prefix_available &&
      wire_bits(path[expected_terminal_index]
                    .selected.global_normalized_continuity) ==
          wire_bits(report.piso.continuity_residual) &&
      wire_bits(path[expected_terminal_index]
                    .selected.global_normalized_energy) ==
          wire_bits(report.piso.energy_residual);
  const bool terminal =
      refined.status && report.accepted && report.attempts == 1U &&
      report.effective_bdf.order == 1U && report.failed_stage == 0U &&
      report.piso.final_flux_revision != 0U &&
      !report.piso.continuity_witness.valid &&
      finite_positive_terminal_state(refined.driver, refined.model, report);

  // Preserve a sub-capacity success certificate as a separate falsifier for
  // the early-exit contract.  Raising the bound must not turn refinement into
  // an unconditional twelve-solve loop.
  constexpr double kEarlyExitContinuityGate = 1.0460433e-10;
  ValidatedModel early_model =
      retry_model(kFullDt, kFullDt, 1U, UINT64_C(0x18000c501));
  early_model.solver.terminal.continuity = kEarlyExitContinuityGate;
  DriverHarness early = make_driver(std::move(early_model));
  DriverStepReport early_report;
  if (early.status)
    early.status =
        early.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, early_report);
  constexpr std::uint8_t kEarlyCalls = 2U;
  constexpr std::size_t kEarlyTerminalIndex = kEarlyCalls + 1U;
  const auto& early_path =
      early_report.pressure_energy_globalization.trajectory;
  const bool early_prefix_available =
      early_report.piso.pressure_energy_refinement_solve_calls ==
          kEarlyCalls &&
      early_report.pressure_energy_globalization.trajectory_count ==
          2U + kEarlyCalls;
  const bool early_exit =
      early_prefix_available && early.status && early_report.accepted &&
      early_report.attempts == 1U && early_report.failed_stage == 0U &&
      refinement_prefix_certificate(
          early_report, kEarlyCalls,
          PressureEnergyRefinementTermination::component_residuals_converged) &&
      refinement_trajectory_certificate(early_report, kEarlyCalls) &&
      early_path[kEarlyTerminalIndex - 1U]
              .selected.global_normalized_energy >
          kEarlyExitContinuityGate &&
      early_path[kEarlyTerminalIndex].selected.global_normalized_energy <=
          kEarlyExitContinuityGate &&
      early_path[kEarlyTerminalIndex]
              .selected.global_normalized_continuity <=
          kEarlyExitContinuityGate &&
      finite_positive_terminal_state(early.driver, early.model, early_report);
  const bool local =
      replay_roundoff_contract && prefix && trajectory &&
      first_crosses_component_gate &&
      terminal_matches_selected_trajectory && terminal && early_exit;
  if (rank == 0)
    std::cout << std::setprecision(17) << "refinement status="
              << static_cast<unsigned>(refined.status.code) << '/'
              << refined.status.detail
              << " accepted/attempts=" << report.accepted << '/'
              << report.attempts << " calls/term="
              << static_cast<unsigned>(
                     report.piso.pressure_energy_refinement_solve_calls)
              << '/'
              << static_cast<unsigned>(
                     report.piso.pressure_energy_refinement_termination)
              << " trajectory="
              << static_cast<unsigned>(
                     report.pressure_energy_globalization.trajectory_count)
              << " terminal=" << report.piso.eos_residual << '/'
              << report.piso.continuity_residual << '/'
              << report.piso.energy_residual << '/'
              << report.piso.closed_mass_residual << '/'
              << report.piso.gauge_residual << " early-exit="
              << static_cast<unsigned>(
                     early_report.piso.pressure_energy_refinement_solve_calls)
              << '/' << early_exit << '\n';
  return expect(collective(local), rank,
                "real ProductDriver C2 refinement refreshes state and "
                "linearization, then passes every independent terminal gate");
}

bool run_retry_certificate(int rank) {
  const ProbeResult probe = run_full_probe(rank);

  DriverHarness retry =
      make_driver(retry_model(kHalfDt, kFullDt, 2U, UINT64_C(0x18000c201)));
  DriverStepReport retry_first;
  if (retry.status)
    retry.status = retry.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, retry_first);

  DriverHarness control = make_driver(
      retry_model(kHalfDt, kHalfDt, 1U, UINT64_C(0x18000c301)), kHalfDt);
  DriverStepReport control_first;
  if (control.status)
    control.status =
        control.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, control_first);

  PhysicalCommittedBits retry_at_first;
  PhysicalCommittedBits control_at_first;
  const bool captured_first =
      capture_physical_committed(retry.driver, retry_at_first) &&
      capture_physical_committed(control.driver, control_at_first);
  const bool same_first =
      captured_first &&
      physical_committed_equal(retry_at_first, control_at_first);
  const bool first_terminal =
      retry.status && control.status &&
      finite_positive_terminal_state(retry.driver, retry.model, retry_first) &&
      finite_positive_terminal_state(control.driver, control.model,
                                     control_first);

  // A/B now expose the same accepted_n state and final phi.  Advancing both
  // through the next BDF2 target makes any leaked rejected-attempt previous
  // layer, rate history, or warm authority observable in the public result.
  DriverStepReport retry_second;
  if (retry.status)
    retry.status =
        retry.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, retry_second);
  DriverStepReport control_second;
  if (control.status)
    control.status =
        control.driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, control_second);

  PhysicalCommittedBits retry_at_second;
  PhysicalCommittedBits control_at_second;
  const bool captured_second =
      capture_physical_committed(retry.driver, retry_at_second) &&
      capture_physical_committed(control.driver, control_at_second);
  const bool same_second =
      captured_second &&
      physical_committed_equal(retry_at_second, control_at_second);

  const bool first_retry_semantics =
      retry_first.accepted && retry_first.attempts == 2U &&
      retry_first.proposal.origin == StepOrigin::retry &&
      wire_bits(retry_first.proposal.dt) == wire_bits(kHalfDt) &&
      retry_first.effective_bdf.order == 1U && !retry_first.failure &&
      retry_first.failure.code == StatusCode::rejected_step &&
      retry_first.failure.detail == 10210U && retry_first.failed_stage == 54U;
  const bool first_control_semantics =
      control_first.accepted && control_first.attempts == 1U &&
      control_first.proposal.origin == StepOrigin::restart &&
      wire_bits(control_first.proposal.dt) == wire_bits(kHalfDt) &&
      control_first.effective_bdf.order == 1U;
  const bool second_semantics =
      retry_second.accepted && control_second.accepted &&
      retry_second.attempts == 1U && control_second.attempts == 1U &&
      wire_bits(retry_second.proposal.dt) == wire_bits(kHalfDt) &&
      wire_bits(control_second.proposal.dt) == wire_bits(kHalfDt) &&
      retry_second.effective_bdf.order == 2U &&
      control_second.effective_bdf.order == 2U;
  const bool terminal =
      retry.status && control.status &&
      finite_positive_terminal_state(retry.driver, retry.model, retry_second) &&
      finite_positive_terminal_state(control.driver, control.model,
                                     control_second);
  const bool local = probe.passed && first_retry_semantics &&
                     first_control_semantics && same_first && first_terminal &&
                     second_semantics && same_second && terminal;
  if (rank == 0)
    std::cout << std::setprecision(17) << "retry/control first status="
              << static_cast<unsigned>(retry.status.code) << '/'
              << retry.status.detail << ','
              << static_cast<unsigned>(control.status.code) << '/'
              << control.status.detail << " attempts=" << retry_first.attempts
              << '/' << control_first.attempts << " origin="
              << static_cast<unsigned>(retry_first.proposal.origin) << '/'
              << static_cast<unsigned>(control_first.proposal.origin)
              << " dt=" << retry_first.proposal.dt << '/'
              << control_first.proposal.dt << " same=" << same_first
              << " terminal1=" << retry_first.piso.eos_residual << '/'
              << retry_first.piso.continuity_residual << '/'
              << retry_first.piso.energy_residual << '/'
              << retry_first.piso.closed_mass_residual << '/'
              << retry_first.piso.gauge_residual
              << " second same=" << same_second << " bdf="
              << static_cast<unsigned>(retry_second.effective_bdf.order) << '/'
              << static_cast<unsigned>(control_second.effective_bdf.order)
              << " terminal2=" << retry_second.piso.eos_residual << '/'
              << retry_second.piso.continuity_residual << '/'
              << retry_second.piso.energy_residual << '/'
              << retry_second.piso.closed_mass_residual << '/'
              << retry_second.piso.gauge_residual << '\n';
  return expect(collective(local), rank,
                "full-dt rejection rolls back exactly, half-dt retry matches "
                "a direct BE control, and the following target remains "
                "bitwise identical and terminal-positive");
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool passed = run_refinement_certificate(rank);
  passed &= run_retry_certificate(rank);
  MPI_Finalize();
  return passed ? 0 : 1;
}
