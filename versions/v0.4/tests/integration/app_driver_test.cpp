// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"

#include "core_product_freeze_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include "../support/product_fixture.hpp"

#include <mpi.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

bool shared_cell_convective_cfl_kernel_has_one_formula() {
  detail::CellConvectiveCflResult result;
  const std::array<double, 6U> fluxes{{-2.0, 3.0, -4.0, 5.0, -6.0,
                                       7.0}};
  const std::array<std::uint8_t, 6U> active{{1U, 1U, 1U, 1U, 1U, 1U}};
  const detail::CellConvectiveCflStatus status =
      detail::evaluate_cell_convective_cfl(2.0, 4.0, fluxes, active, 0.5,
                                           result);
  detail::CellConvectiveCflResult sealed_result;
  auto sealed = active;
  sealed[0U] = 0U;
  const detail::CellConvectiveCflStatus sealed_status =
      detail::evaluate_cell_convective_cfl(2.0, 4.0, fluxes, sealed, 0.5,
                                           sealed_result);
  detail::CellConvectiveCflResult nonphysical_result;
  const detail::CellConvectiveCflStatus nonphysical_status =
      detail::evaluate_cell_convective_cfl(0.0, 4.0, fluxes, active, 0.5,
                                           nonphysical_result);
  return status == detail::CellConvectiveCflStatus::success &&
         result.density_volume == 8.0 &&
         result.outgoing_mass_flow == 27.0 &&
         result.absolute_mass_flow == 27.0 && result.out == 1.6875 &&
         result.absolute == 0.84375 &&
         sealed_status ==
             detail::CellConvectiveCflStatus::inactive_nonzero_flux &&
         nonphysical_status ==
             detail::CellConvectiveCflStatus::nonphysical_state;
}

bool copy_restart_image(const RestartSnapshot& snapshot,
                        const RestartExpected& expected,
                        RestartImage& image) {
  if (snapshot.fields.size != expected.fields.size ||
      snapshot.patch.cells.x != expected.target_patch.cells.x ||
      snapshot.patch.cells.y != expected.target_patch.cells.y ||
      snapshot.patch.cells.z != expected.target_patch.cells.z)
    return false;
  RestartImage candidate;
  candidate.global_cells = snapshot.global_cells;
  candidate.patch = snapshot.patch;
  candidate.plan = snapshot.plan;
  candidate.schema = snapshot.schema;
  candidate.geometry = snapshot.geometry;
  candidate.time = snapshot.time;
  candidate.dt = snapshot.dt;
  candidate.pressure_reference = snapshot.pressure_reference;
  candidate.step = snapshot.step;
  candidate.controller_state = snapshot.controller_state;
  candidate.backward_euler_recovery = true;
  candidate.fields.reserve(snapshot.fields.size);
  const Int3 cells = snapshot.patch.cells;
  for (std::size_t field_index = 0U; field_index < snapshot.fields.size;
       ++field_index) {
    const RestartFieldView source = snapshot.fields.data[field_index];
    const RestartExpectedField descriptor = expected.fields.data[field_index];
    if (source.role != descriptor.role ||
        source.values.components != descriptor.components)
      return false;
    RestartImageField target;
    target.role = descriptor.role;
    target.field = descriptor.field;
    target.components = descriptor.components;
    target.values.resize(static_cast<std::size_t>(cells.x) * cells.y *
                         cells.z * descriptor.components);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const std::size_t cell = static_cast<std::size_t>(x) +
                                   static_cast<std::size_t>(cells.x) *
                                       (static_cast<std::size_t>(y) +
                                        static_cast<std::size_t>(cells.y) * z);
          for (std::uint8_t component = 0U;
               component < descriptor.components; ++component)
            target.values[cell * descriptor.components + component] =
                source.values.unchecked({x, y, z}, component);
        }
    candidate.fields.push_back(std::move(target));
  }
  const std::array<ConstFaceFieldView, 3U> source_faces{
      snapshot.final_mass_flux.x, snapshot.final_mass_flux.y,
      snapshot.final_mass_flux.z};
  for (std::size_t axis = 0U; axis < source_faces.size(); ++axis) {
    const ConstFaceFieldView source = source_faces[axis];
    std::vector<double>& target = candidate.final_mass_flux[axis];
    target.reserve(static_cast<std::size_t>(source.extents.x) *
                   source.extents.y * source.extents.z);
    for (std::int32_t z = 0; z < source.extents.z; ++z)
      for (std::int32_t y = 0; y < source.extents.y; ++y)
        for (std::int32_t x = 0; x < source.extents.x; ++x)
          target.push_back(source.unchecked({x, y, z}));
  }
  image = std::move(candidate);
  return true;
}

bool terminal_physical_certificate(ProductDriver& driver,
                                   const ValidatedModel& model,
                                   const DriverStepReport& report,
                                   RestartSnapshot* published = nullptr) {
  RestartSnapshot snapshot;
  Status status = driver.committed_restart_snapshot(snapshot);
  const auto accepted_solve = [](const LinearSolveResult& solve) {
    return solve.status &&
           (solve.termination == LinearTermination::converged ||
            solve.termination == LinearTermination::zero_rhs) &&
           std::isfinite(solve.initial_true_residual) &&
           std::isfinite(solve.final_true_residual) &&
           solve.final_true_residual <= solve.initial_true_residual;
  };
  if (!status || !report.accepted || report.piso.pressure_solve_calls != 2U ||
      !accepted_solve(report.piso.pressure[0U]) ||
      !accepted_solve(report.piso.pressure[1U]) ||
      report.piso.final_flux_revision == 0U ||
      !report.momentum_predictor_limiter.advective_cfl.valid() ||
      report.momentum_predictor_limiter.advective_cfl.face_flux ==
          report.piso.final_flux_revision ||
      report.momentum_predictor_limiter.advective_cfl.limit !=
          model.time.convective_cfl ||
      report.momentum_predictor_limiter.advective_cfl.out_max >
          report.momentum_predictor_limiter.advective_cfl.limit *
              (1.0 + 64.0 * std::numeric_limits<double>::epsilon()) ||
      !std::isfinite(report.piso.committed_convective_cfl_out_max) ||
      report.piso.committed_convective_cfl_out_max < 0.0 ||
      !std::isfinite(report.piso.committed_convective_cfl_abs_max) ||
      report.piso.committed_convective_cfl_abs_max < 0.0 ||
      !std::isfinite(report.piso.committed_convective_cfl_limit) ||
      report.piso.committed_convective_cfl_limit !=
          model.time.convective_cfl ||
      report.piso.committed_convective_cfl_out_max >
          report.piso.committed_convective_cfl_limit *
              (1.0 + 64.0 * std::numeric_limits<double>::epsilon()) ||
      !snapshot.final_mass_flux.certificate.valid() ||
      snapshot.final_mass_flux.revision != report.piso.final_flux_revision ||
      snapshot.step != report.accepted_step ||
      snapshot.time != report.accepted_time ||
      !std::isfinite(report.proposal.time) ||
      !std::isfinite(report.proposal.dt) || !(report.proposal.dt > 0.0) ||
      report.accepted_time != report.proposal.time + report.proposal.dt ||
      report.effective_bdf.order != report.proposal.bdf.order)
    return false;

  const auto accepted_residual = [](double value, double tolerance) {
    return std::isfinite(value) && value >= 0.0 && value <= tolerance;
  };
  if (!accepted_residual(report.piso.eos_residual,
                         model.solver.terminal.eos) ||
      !accepted_residual(report.piso.continuity_residual,
                         model.solver.terminal.continuity) ||
      !accepted_residual(report.piso.energy_residual,
                         model.solver.terminal.continuity) ||
      !accepted_residual(report.piso.closed_mass_residual,
                         model.solver.terminal.closed_mass) ||
      !accepted_residual(report.piso.gauge_residual,
                         model.solver.terminal.gauge))
    return false;

  const ConstFieldView* pressure = nullptr;
  const ConstFieldView* enthalpy = nullptr;
  const ConstFieldView* velocity = nullptr;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView& field = snapshot.fields.data[index];
    if (field.role == RestartFieldRole::pressure_perturbation)
      pressure = &field.values;
    else if (field.role == RestartFieldRole::enthalpy)
      enthalpy = &field.values;
    else if (field.role == RestartFieldRole::velocity)
      velocity = &field.values;
  }
  if (pressure == nullptr || enthalpy == nullptr || velocity == nullptr ||
      pressure->components != 1U || enthalpy->components != 1U ||
      velocity->components != 3U)
    return false;

  ThermodynamicsPlan thermodynamics;
  status = ThermodynamicsPlan::compile(
      model.thermophysics,
      {model.transported_scalars.data(), model.transported_scalars.size()},
      thermodynamics);
  if (!status) return false;
  for (std::int32_t z = 0; z < pressure->interior.z; ++z)
    for (std::int32_t y = 0; y < pressure->interior.y; ++y)
      for (std::int32_t x = 0; x < pressure->interior.x; ++x) {
        const Int3 cell{x, y, z};
        const double pressure_absolute =
            snapshot.pressure_reference + pressure->unchecked(cell, 0U);
        const Real3 cell_velocity{velocity->unchecked(cell, 0U),
                                  velocity->unchecked(cell, 1U),
                                  velocity->unchecked(cell, 2U)};
        ThermoState state;
        status = thermodynamics.evaluate(
            pressure_absolute, enthalpy->unchecked(cell, 0U), {},
            cell_velocity, state);
        if (!status || !std::isfinite(pressure_absolute) ||
            !(pressure_absolute > 0.0) || !std::isfinite(state.rho) ||
            !(state.rho > 0.0) || !std::isfinite(state.temperature) ||
            state.temperature < model.thermophysics.minimum_temperature ||
            state.temperature > model.thermophysics.maximum_temperature ||
            !std::isfinite(cell_velocity.x) ||
            !std::isfinite(cell_velocity.y) ||
            !std::isfinite(cell_velocity.z))
          return false;
      }

  const std::array<ConstFaceFieldView, 3U> flux{{snapshot.final_mass_flux.x,
                                                 snapshot.final_mass_flux.y,
                                                 snapshot.final_mass_flux.z}};
  for (const ConstFaceFieldView& axis : flux)
    for (std::int32_t z = 0; z < axis.extents.z; ++z)
      for (std::int32_t y = 0; y < axis.extents.y; ++y)
        for (std::int32_t x = 0; x < axis.extents.x; ++x)
          if (!std::isfinite(axis.unchecked({x, y, z}))) return false;

  if (published != nullptr) *published = snapshot;
  return true;
}

bool open_x_boundary_flux_certificate(const RestartSnapshot& snapshot,
                                      double mass_flow_target = 0.0) {
  if (snapshot.patch.begin.x != 0 || snapshot.patch.begin.y != 0 ||
      snapshot.patch.begin.z != 0 ||
      snapshot.patch.cells.x != snapshot.global_cells.x ||
      snapshot.patch.cells.y != snapshot.global_cells.y ||
      snapshot.patch.cells.z != snapshot.global_cells.z)
    return false;
  double inlet_mass_flow = 0.0;
  for (std::int32_t z = 0; z < snapshot.patch.cells.z; ++z)
    for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y) {
      const double inlet =
          snapshot.final_mass_flux.x.unchecked({0, y, z});
      const double outlet = snapshot.final_mass_flux.x.unchecked(
          {snapshot.patch.cells.x, y, z});
      if (!(inlet > 0.0) || outlet < 0.0) return false;
      inlet_mass_flow += inlet;
    }
  for (std::int32_t z = 0; z < snapshot.patch.cells.z; ++z)
    for (std::int32_t x = 0; x < snapshot.patch.cells.x; ++x)
      for (const std::int32_t y : {0, snapshot.patch.cells.y}) {
        const double value = snapshot.final_mass_flux.y.unchecked({x, y, z});
        if (value != 0.0 || std::signbit(value)) return false;
      }
  for (std::int32_t y = 0; y < snapshot.patch.cells.y; ++y)
    for (std::int32_t x = 0; x < snapshot.patch.cells.x; ++x)
      for (const std::int32_t z : {0, snapshot.patch.cells.z}) {
        const double value = snapshot.final_mass_flux.z.unchecked({x, y, z});
        if (value != 0.0 || std::signbit(value)) return false;
      }
  if (mass_flow_target == 0.0) return true;
  const double tolerance =
      128.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::abs(mass_flow_target));
  return std::abs(inlet_mass_flow - mass_flow_target) <= tolerance;
}

LocalTimeLimits limits_for_exact_dt(const ValidatedModel& model,
                                    double dt) {
  const double smallest_target =
      std::min({model.time.convective_cfl, model.time.viscous_cfl,
                model.time.thermal_cfl, model.time.species_cfl});
  const double scale = dt / smallest_target;
  return {scale, scale, scale, scale, scale};
}

struct DriverFieldImage {
  std::string stable_name;
  ConstFieldView view{};
  RevisionToken accepted_revision{};
  std::vector<double> values;
};

struct DriverStateImage {
  const CartesianGeometryPlan* geometry{};
  MeshPatch patch{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  double time{};
  std::uint64_t step{};
  std::vector<DriverFieldImage> fields;
  RevisionToken flux_revision{};
  FaceFluxCertificate flux_certificate{};
  std::array<std::vector<double>, 3U> flux_values;
  double pressure_reference{};
  double closed_mass{};
};

bool capture_driver_state(ProductDriver& driver, DriverStateImage& out) {
  CommittedOutputSnapshot snapshot;
  ConstFaceFluxView flux;
  Status status = driver.committed_output_snapshot(snapshot);
  if (status) status = driver.committed_final_mass_flux_for_test(flux);
  if (!status || !snapshot.committed || !flux.certificate.valid())
    return false;
  DriverStateImage candidate;
  candidate.geometry = snapshot.geometry;
  candidate.patch = snapshot.patch;
  candidate.plan = snapshot.plan;
  candidate.schema = snapshot.schema;
  candidate.time = snapshot.time;
  candidate.step = snapshot.step;
  candidate.fields.reserve(snapshot.fields.size);
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const SnapshotFieldView& source = snapshot.fields.data[index];
    DriverFieldImage field;
    field.stable_name = source.stable_name;
    field.view = source.values;
    field.accepted_revision = source.accepted_revision;
    const ConstFieldView view = source.values;
    field.values.reserve(static_cast<std::size_t>(view.interior.x) *
                         view.interior.y * view.interior.z * view.components);
    for (std::uint8_t component = 0U; component < view.components; ++component)
      for (std::int32_t z = 0; z < view.interior.z; ++z)
        for (std::int32_t y = 0; y < view.interior.y; ++y)
          for (std::int32_t x = 0; x < view.interior.x; ++x)
            field.values.push_back(view.unchecked({x, y, z}, component));
    candidate.fields.push_back(std::move(field));
  }
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
    const ConstFaceFieldView face = faces[axis];
    std::vector<double>& values = candidate.flux_values[axis];
    values.reserve(static_cast<std::size_t>(face.extents.x) * face.extents.y *
                   face.extents.z);
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          values.push_back(face.unchecked({x, y, z}));
  }
  candidate.flux_revision = flux.revision;
  candidate.flux_certificate = flux.certificate;
  candidate.pressure_reference = driver.pressure_reference();
  candidate.closed_mass = driver.closed_mass_target();
  out = std::move(candidate);
  return true;
}

bool exact_double_vector(const std::vector<double>& left,
                         const std::vector<double>& right) {
  return left.size() == right.size() &&
         (left.empty() ||
          std::memcmp(left.data(), right.data(),
                      left.size() * sizeof(double)) == 0);
}

bool same_driver_fields(const DriverStateImage& left,
                        const DriverStateImage& right) {
  if (left.fields.size() != right.fields.size()) return false;
  for (std::size_t index = 0U; index < left.fields.size(); ++index) {
    const DriverFieldImage& a = left.fields[index];
    const DriverFieldImage& b = right.fields[index];
    if (a.stable_name != b.stable_name ||
        a.view.field != b.view.field ||
        a.view.components != b.view.components ||
        a.view.interior.x != b.view.interior.x ||
        a.view.interior.y != b.view.interior.y ||
        a.view.interior.z != b.view.interior.z ||
        !exact_double_vector(a.values, b.values))
      return false;
  }
  for (std::size_t axis = 0U; axis < left.flux_values.size(); ++axis)
    if (!exact_double_vector(left.flux_values[axis], right.flux_values[axis]))
      return false;
  return true;
}

bool same_patch(MeshPatch left, MeshPatch right) {
  const auto same = [](Int3 a, Int3 b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  };
  return same(left.begin, right.begin) && same(left.cells, right.cells) &&
         same(left.process_grid, right.process_grid) &&
         same(left.process_coord, right.process_coord);
}

bool same_physical_state(const DriverStateImage& left,
                         const DriverStateImage& right) {
  return left.time == right.time && left.step == right.step &&
         left.pressure_reference == right.pressure_reference &&
         left.closed_mass == right.closed_mass &&
         same_patch(left.patch, right.patch) && same_driver_fields(left, right);
}

bool same_exact_state(const DriverStateImage& left,
                      const DriverStateImage& right) {
  if (!same_physical_state(left, right) || left.geometry != right.geometry ||
      left.plan != right.plan || left.schema != right.schema ||
      left.flux_revision != right.flux_revision ||
      !(left.flux_certificate == right.flux_certificate))
    return false;
  for (std::size_t index = 0U; index < left.fields.size(); ++index) {
    const DriverFieldImage& a = left.fields[index];
    const DriverFieldImage& b = right.fields[index];
    const ConstFieldView& x = a.view;
    const ConstFieldView& y = b.view;
    if (a.accepted_revision != b.accepted_revision || x.base != y.base ||
        x.ghosts.x != y.ghosts.x || x.ghosts.y != y.ghosts.y ||
        x.ghosts.z != y.ghosts.z || x.stride_y != y.stride_y ||
        x.stride_z != y.stride_z ||
        x.component_stride != y.component_stride ||
        x.replica != y.replica || x.revision != y.revision ||
        x.storage_identity != y.storage_identity ||
        x.revision_domain != y.revision_domain)
      return false;
  }
  return true;
}

bool run() {
  const fs::path root = fs::temp_directory_path() /
                        ("hundun-v04-driver-plan-" +
                         std::to_string(::getpid()));
  std::error_code error;
  fs::remove_all(root, error);
  if (!ApplicationService::initialize_case_directory(root)) return false;
  CaseValidationReport first;
  CaseValidationReport second;
  bool passed = static_cast<bool>(
      ApplicationService::validate(MPI_COMM_SELF, root, first));
  passed &= static_cast<bool>(
      ApplicationService::validate(MPI_COMM_SELF, root, second));
  passed &= first.case_model == second.case_model &&
            first.product == second.product && first.summary.sealed &&
            first.summary.graph_stage_count == 15U &&
            first.summary.pressure_correctors == 2U &&
            !first.summary.exact_numeric_certified &&
            !first.summary.preconditioner_setup_certified;
  fs::remove_all(root, error);
  return passed;
}

bool run_initialized_product() {
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(
      MPI_COMM_SELF, test::product_model({8, 7, 6}), {}, plan);
  if (!status) return false;
  ProductDriver driver;
  status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  if (!status || plan.fingerprint() != 0U || driver.initialized()) return false;

  DriverInitialState invalid;
  invalid.temperature = std::numeric_limits<double>::quiet_NaN();
  status = driver.initialize(invalid);
  if (status.code != StatusCode::invalid_plan || driver.initialized())
    return false;

  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {0.0, 0.0, 0.0};
  initial.start_time = 2.0;
  status = driver.initialize(initial);
  if (!status || !driver.initialized() ||
      driver.pressure_reference() != initial.pressure_reference ||
      !(driver.closed_mass_target() > 0.0)) {
    return false;
  }

  CommittedOutputSnapshot snapshot;
  status = driver.committed_output_snapshot(snapshot);
  if (!status || !snapshot.committed || snapshot.step != 0U ||
      snapshot.time != initial.start_time || snapshot.plan == 0U ||
      snapshot.schema == 0U || snapshot.fields.size != 3U) {
    return false;
  }
  bool saw_velocity = false;
  bool saw_pressure = false;
  bool saw_enthalpy = false;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const SnapshotFieldView& field = snapshot.fields.data[index];
    if (field.values.base == nullptr || field.accepted_revision == 0U)
      return false;
    if (field.stable_name == "U") {
      saw_velocity = field.values.components == 3U &&
                     field.values.unchecked({0, 0, 0}, 0U) ==
                         initial.velocity.x &&
                     field.values.unchecked({0, 0, 0}, 1U) ==
                         initial.velocity.y &&
                     field.values.unchecked({0, 0, 0}, 2U) ==
                         initial.velocity.z;
    } else if (field.stable_name == "pi") {
      saw_pressure = field.values.components == 1U &&
                     field.values.unchecked({0, 0, 0}, 0U) == 0.0;
    } else if (field.stable_name == "h") {
      saw_enthalpy = field.values.components == 1U &&
                     std::isfinite(field.values.unchecked({0, 0, 0}, 0U));
    }
  }
  if (!saw_velocity || !saw_pressure || !saw_enthalpy) return false;

  RestartSnapshot restart;
  status = driver.committed_restart_snapshot(restart);
  if (status.code != StatusCode::invalid_plan) return false;
  status = driver.initialize(initial);
  if (status.code != StatusCode::invalid_plan || !driver.initialized())
    return false;
  CommittedOutputSnapshot after_rejected_reinitialize;
  status = driver.committed_output_snapshot(after_rejected_reinitialize);
  if (!(status && after_rejected_reinitialize.plan == snapshot.plan &&
         after_rejected_reinitialize.schema == snapshot.schema &&
         after_rejected_reinitialize.time == snapshot.time &&
         after_rejected_reinitialize.step == snapshot.step))
    return false;

  DriverStepReport step;
  const LocalTimeLimits limits{1.0, 1.0, 1.0, 1.0, 1.0};
  status = driver.advance(limits, step);
  if (!status) {
    std::cerr << "advance status=" << static_cast<unsigned>(status.code)
              << " detail=" << status.detail
              << " stage=" << step.failed_stage
              << " attempts=" << step.attempts
              << " eos=" << step.piso.eos_residual
              << " continuity=" << step.piso.continuity_residual
              << " mass=" << step.piso.closed_mass_residual
              << " gauge=" << step.piso.gauge_residual << '\n';
    return false;
  }
  if (!step.accepted || step.accepted_step != 1U ||
      step.accepted_time != initial.start_time + 1.0e-3 ||
      step.piso.pressure_solve_calls != 2U ||
      !step.piso.pressure[0U].status || !step.piso.pressure[1U].status) {
    std::cerr << "accepted=" << step.accepted
              << " accepted_step=" << step.accepted_step
              << " accepted_time=" << step.accepted_time
              << " attempts=" << step.attempts
              << " prior_failure="
              << static_cast<unsigned>(step.failure.code)
              << "/" << step.failure.detail
              << " prior_stage=" << step.failed_stage
              << " solve_calls="
              << static_cast<unsigned>(step.piso.pressure_solve_calls)
              << " solve0="
              << static_cast<unsigned>(step.piso.pressure[0U].status.code)
              << " solve1="
              << static_cast<unsigned>(step.piso.pressure[1U].status.code)
              << '\n';
    return false;
  }
  CommittedOutputSnapshot advanced;
  status = driver.committed_output_snapshot(advanced);
  if (!status || advanced.step != 1U ||
      advanced.time != step.accepted_time) {
    std::cerr << "output status=" << static_cast<unsigned>(status.code)
              << " step=" << advanced.step << " time=" << advanced.time
              << '\n';
    return false;
  }
  status = driver.committed_restart_snapshot(restart);
  const bool restart_ok =
      status && restart.step == 1U && restart.time == step.accepted_time &&
      restart.fields.size == 3U && restart.final_mass_flux.revision != 0U;
  if (!restart_ok)
    std::cerr << "restart status=" << static_cast<unsigned>(status.code)
              << " detail=" << status.detail << " step=" << restart.step
              << " time=" << restart.time
              << " fields=" << restart.fields.size
              << " flux=" << restart.final_mass_flux.revision << '\n';
  return restart_ok;
}

bool run_initialized_species_product() {
  ValidatedModel model = test::product_model({8, 7, 6});
  model.transported_scalars.push_back(
      TransportedScalarSpec{"tracer", TransportedScalarRole::passive_scalar});
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  if (!status) {
    std::cerr << "species compile=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << '\n';
    return false;
  }
  ProductDriver driver;
  status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  if (!status) {
    std::cerr << "species create=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << '\n';
    return false;
  }
  const std::array<double, 1U> composition{0.25};
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.transported_scalars = {composition.data(), composition.size()};
  status = driver.initialize(initial);
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  if (!status) {
    std::cerr << "species advance=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << " stage=" << step.failed_stage
              << " attempts=" << step.attempts << '\n';
    return false;
  }
  const bool first_physical =
      terminal_physical_certificate(driver, model, step);
  if (!first_physical) return false;
  DriverStepReport second;
  status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);
  if (!status) {
    std::cerr << "species bdf2 advance="
              << static_cast<unsigned>(status.code) << "/" << status.detail
              << " stage=" << second.failed_stage
              << " attempts=" << second.attempts << '\n';
    return false;
  }
  RestartSnapshot restart;
  const bool second_physical =
      terminal_physical_certificate(driver, model, second, &restart);
  const bool passed = second_physical && second.accepted_step == 2U &&
                      restart.step == 2U && restart.fields.size == 4U;
  if (!passed)
    std::cerr << "species final=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << " accepted=" << step.accepted
              << " attempts=" << step.attempts
              << " fields=" << restart.fields.size << '\n';
  return passed;
}

bool reject_invalid_thermodynamic_restart_atomically() {
  const ValidatedModel model = test::product_model({8, 7, 6});
  CompiledCasePlan source_plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {},
                                            source_plan);
  ProductDriver source;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(source_plan),
                                   source);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = source.initialize(initial);
  DriverStepReport step;
  if (status)
    status = source.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  RestartSnapshot snapshot;
  if (status) status = source.committed_restart_snapshot(snapshot);

  CompiledCasePlan target_plan;
  if (status)
    status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, target_plan);
  ProductDriver target;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(target_plan),
                                   target);
  RestartExpected expected;
  if (status) status = target.restart_expected(expected);
  RestartImage valid;
  if (status && !copy_restart_image(snapshot, expected, valid))
    status = {StatusCode::invalid_plan, 0U};
  if (!status) {
    std::cerr << "restart atomic fixture="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
    return false;
  }

  RestartImage invalid = valid;
  bool mutated = false;
  for (RestartImageField& field : invalid.fields) {
    if (field.role == RestartFieldRole::enthalpy && !field.values.empty()) {
      field.values[0U] = -std::numeric_limits<double>::max();
      mutated = true;
      break;
    }
  }
  if (!mutated) return false;
  status = target.initialize_restart(invalid);
  if (status || status.code != StatusCode::numerical_failure ||
      target.initialized() || target.pressure_reference() != 0.0 ||
      target.closed_mass_target() != 0.0) {
    std::cerr << "invalid restart published state=" << target.initialized()
              << " status=" << static_cast<unsigned>(status.code) << '/'
              << status.detail << '\n';
    return false;
  }

  status = target.initialize_restart(valid);
  CommittedOutputSnapshot restored;
  if (status) status = target.committed_output_snapshot(restored);
  if (!status)
    std::cerr << "valid restart after rejection="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
  return status && target.initialized() && restored.committed &&
         restored.step == snapshot.step && restored.time == snapshot.time;
}

bool restart_preserves_direct_flux_lineage() {
  ValidatedModel model = test::product_model({8, 7, 6});
  // A restart carries its own thermodynamic state.  Keep 500 K admissible but
  // exclude DriverInitialState's 300 K default so any accidental delegation
  // through the fresh-start initializer is observable at product level.
  model.thermophysics.minimum_temperature = 400.0;

  CompiledCasePlan source_plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {},
                                            source_plan);
  ProductDriver source;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(source_plan),
                                   source);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 500.0;
  initial.velocity = {0.2, -0.1, 0.0};
  if (status) status = source.initialize(initial);
  DriverStepReport source_step;
  if (status)
    status = source.advance({1.0, 1.0, 1.0, 1.0, 1.0}, source_step);
  RestartSnapshot source_snapshot;
  if (status) status = source.committed_restart_snapshot(source_snapshot);

  CompiledCasePlan target_plan;
  if (status)
    status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, target_plan);
  ProductDriver target;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(target_plan),
                                   target);
  RestartExpected expected;
  if (status) status = target.restart_expected(expected);
  RestartImage image;
  if (status && !copy_restart_image(source_snapshot, expected, image))
    status = {StatusCode::invalid_plan, 0U};
  if (status && image.final_mass_flux[2U].size() < 2U)
    status = {StatusCode::invalid_plan, 0U};
  if (status) {
    // Both values are physically zero, but their object representations are
    // different.  Together with the nonzero x/y fluxes this makes the exact
    // restore assertion sensitive to byte-level copying, not double ==.
    image.final_mass_flux[2U][0U] = -0.0;
    image.final_mass_flux[2U][1U] = +0.0;
  }
  if (status) status = target.initialize_restart(image);
  RestartSnapshot restored;
  if (status) status = target.committed_restart_snapshot(restored);
  if (!status) {
    std::cerr << "direct restart status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " initialized=" << target.initialized() << '\n';
    return false;
  }

  const std::array<ConstFaceFieldView, 3U> restored_faces{
      restored.final_mass_flux.x, restored.final_mass_flux.y,
      restored.final_mass_flux.z};
  bool exact_flux = restored.final_mass_flux.revision == 2U &&
                    std::signbit(image.final_mass_flux[2U][0U]) &&
                    !std::signbit(image.final_mass_flux[2U][1U]);
  bool saw_nonzero_flux = false;
  for (std::size_t axis = 0U; axis < restored_faces.size() && exact_flux;
       ++axis) {
    const ConstFaceFieldView face = restored_faces[axis];
    std::size_t offset = 0U;
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x, ++offset) {
          const double restored_value = face.unchecked({x, y, z});
          saw_nonzero_flux = saw_nonzero_flux || restored_value != 0.0;
          exact_flux = std::memcmp(&restored_value,
                                   &image.final_mass_flux[axis][offset],
                                   sizeof(double)) == 0;
        }
  }
  if (!exact_flux || !saw_nonzero_flux ||
      !std::signbit(restored.final_mass_flux.z.unchecked({0, 0, 0})) ||
      std::signbit(restored.final_mass_flux.z.unchecked({1, 0, 0}))) {
    std::cerr << "direct restart flux revision="
              << restored.final_mass_flux.revision << '\n';
    return false;
  }

  DriverStepReport restarted_step;
  status = target.advance({1.0, 1.0, 1.0, 1.0, 1.0}, restarted_step);
  RestartSnapshot advanced;
  if (status) status = target.committed_restart_snapshot(advanced);
  if (!status || !restarted_step.accepted ||
      restarted_step.proposal.origin != StepOrigin::restart ||
      restarted_step.proposal.bdf.order != 1U ||
      advanced.final_mass_flux.revision != 3U) {
    std::cerr << "post-restart flux status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted=" << restarted_step.accepted
              << " origin/order="
              << static_cast<unsigned>(restarted_step.proposal.origin) << '/'
              << static_cast<unsigned>(restarted_step.proposal.bdf.order)
              << " revision=" << advanced.final_mass_flux.revision << '\n';
    return false;
  }
  return true;
}

bool exact_restart_preserves_bdf2_history() {
  const fs::path restart_root =
      fs::temp_directory_path() /
      ("hundun-v04-exact-restart-" + std::to_string(::getpid()));
  std::error_code cleanup_error;
  fs::remove_all(restart_root, cleanup_error);

  ValidatedModel model = test::product_model({8, 7, 6});
  CompiledCasePlan source_plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_SELF, model, {}, source_plan);
  ProductDriver source;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(source_plan),
                                   source);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 500.0;
  initial.velocity = {0.2, -0.1, 0.0};
  if (status) status = source.initialize(initial);
  DriverStepReport first;
  DriverStepReport second;
  if (status)
    status = source.advance({1.0, 1.0, 1.0, 1.0, 1.0}, first);
  if (status)
    status = source.advance({1.0, 1.0, 1.0, 1.0, 1.0}, second);
  RestartSnapshot checkpoint;
  if (status) status = source.committed_restart_snapshot(checkpoint);
  if (status)
    status = RestartWriter::write(MPI_COMM_SELF, restart_root, checkpoint,
                                  {1U});

  CompiledCasePlan resumed_plan;
  if (status)
    status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, resumed_plan);
  ProductDriver resumed;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(resumed_plan),
                                   resumed);
  RestartExpected expected;
  if (status) status = resumed.restart_expected(expected);
  RestartImage image;
  if (status)
    status = RestartReader::load(MPI_COMM_SELF, restart_root, expected,
                                 image);
  const bool exact_image =
      status && !image.backward_euler_recovery &&
      image.previous_fields.size() == image.fields.size() &&
      !image.accepted_rate_fields.empty() &&
      image.previous_rate_fields.size() ==
          image.accepted_rate_fields.size() &&
      image.previous_mass_flux_revision != 0U &&
      image.previous_mass_flux_revision < image.final_mass_flux_revision;
  bool failure_atomic = false;
  if (status && exact_image) {
    RestartImage invalid = image;
    invalid.previous_rate_fields[0U].values[0U] =
        std::numeric_limits<double>::quiet_NaN();
    const Status rejected = resumed.initialize_restart(invalid);
    failure_atomic = !rejected &&
                     rejected.code == StatusCode::numerical_failure &&
                     !resumed.initialized() &&
                     resumed.pressure_reference() == 0.0 &&
                     resumed.closed_mass_target() == 0.0;
    if (!failure_atomic) status = rejected;
  }
  if (status) status = resumed.initialize_restart(image);
  RestartSnapshot restored;
  if (status) status = resumed.committed_restart_snapshot(restored);

  const auto same_field_payload = [](Span<const RestartFieldView> left,
                                     Span<const RestartFieldView> right) {
    if (left.size != right.size) return false;
    for (std::size_t field = 0U; field < left.size; ++field) {
      const RestartFieldView a = left.data[field];
      const RestartFieldView b = right.data[field];
      if (a.role != b.role || a.values.field != b.values.field ||
          a.values.components != b.values.components ||
          a.values.interior.x != b.values.interior.x ||
          a.values.interior.y != b.values.interior.y ||
          a.values.interior.z != b.values.interior.z)
        return false;
      for (std::uint8_t component = 0U;
           component < a.values.components; ++component)
        for (std::int32_t z = 0; z < a.values.interior.z; ++z)
          for (std::int32_t y = 0; y < a.values.interior.y; ++y)
            for (std::int32_t x = 0; x < a.values.interior.x; ++x) {
              const double av = a.values.unchecked({x, y, z}, component);
              const double bv = b.values.unchecked({x, y, z}, component);
              if (std::memcmp(&av, &bv, sizeof(double)) != 0) return false;
            }
    }
    return true;
  };
  const auto same_flux_payload = [](ConstFaceFluxView left,
                                    ConstFaceFluxView right) {
    if (left.revision != right.revision) return false;
    const std::array<ConstFaceFieldView, 3U> a{left.x, left.y, left.z};
    const std::array<ConstFaceFieldView, 3U> b{right.x, right.y, right.z};
    for (std::size_t axis = 0U; axis < a.size(); ++axis) {
      if (a[axis].extents.x != b[axis].extents.x ||
          a[axis].extents.y != b[axis].extents.y ||
          a[axis].extents.z != b[axis].extents.z)
        return false;
      for (std::int32_t z = 0; z < a[axis].extents.z; ++z)
        for (std::int32_t y = 0; y < a[axis].extents.y; ++y)
          for (std::int32_t x = 0; x < a[axis].extents.x; ++x) {
            const double av = a[axis].unchecked({x, y, z});
            const double bv = b[axis].unchecked({x, y, z});
            if (std::memcmp(&av, &bv, sizeof(double)) != 0) return false;
          }
    }
    return true;
  };
  const bool exact_restore =
      status && exact_image && restored.controller_state ==
                                   checkpoint.controller_state &&
      restored.previous_pressure_reference ==
          checkpoint.previous_pressure_reference &&
      restored.closed_mass_target == checkpoint.closed_mass_target &&
      same_field_payload(restored.fields, checkpoint.fields) &&
      same_field_payload(restored.previous_fields,
                         checkpoint.previous_fields) &&
      same_field_payload(restored.accepted_rate_fields,
                         checkpoint.accepted_rate_fields) &&
      same_field_payload(restored.previous_rate_fields,
                         checkpoint.previous_rate_fields) &&
      same_flux_payload(restored.final_mass_flux,
                        checkpoint.final_mass_flux) &&
      same_flux_payload(restored.previous_mass_flux,
                        checkpoint.previous_mass_flux);

  DriverStepReport continuous_step;
  DriverStepReport resumed_step;
  if (status)
    status = source.advance({1.0, 1.0, 1.0, 1.0, 1.0}, continuous_step);
  if (status)
    status = resumed.advance({1.0, 1.0, 1.0, 1.0, 1.0}, resumed_step);
  RestartSnapshot continuous_endpoint;
  RestartSnapshot resumed_endpoint;
  if (status)
    status = source.committed_restart_snapshot(continuous_endpoint);
  if (status)
    status = resumed.committed_restart_snapshot(resumed_endpoint);

  double maximum_relative_difference = 0.0;
  const auto compare_field_payload = [&](Span<const RestartFieldView> left,
                                         Span<const RestartFieldView> right) {
    if (left.size != right.size) return false;
    for (std::size_t field = 0U; field < left.size; ++field) {
      const ConstFieldView a = left.data[field].values;
      const ConstFieldView b = right.data[field].values;
      if (left.data[field].role != right.data[field].role ||
          a.components != b.components) return false;
      for (std::uint8_t component = 0U; component < a.components; ++component)
        for (std::int32_t z = 0; z < a.interior.z; ++z)
          for (std::int32_t y = 0; y < a.interior.y; ++y)
            for (std::int32_t x = 0; x < a.interior.x; ++x) {
              const double av = a.unchecked({x, y, z}, component);
              const double bv = b.unchecked({x, y, z}, component);
              maximum_relative_difference = std::max(
                  maximum_relative_difference,
                  std::abs(av - bv) /
                      std::max({1.0, std::abs(av), std::abs(bv)}));
            }
    }
    return true;
  };
  const bool endpoint_layout =
      status && compare_field_payload(continuous_endpoint.fields,
                                      resumed_endpoint.fields);
  const bool passed =
      endpoint_layout && failure_atomic && exact_restore && first.accepted &&
      second.accepted && second.effective_bdf.order == 2U &&
      continuous_step.accepted && resumed_step.accepted &&
      continuous_step.effective_bdf.order == 2U &&
      resumed_step.proposal.origin == StepOrigin::restart &&
      resumed_step.proposal.bdf.order == 2U &&
      resumed_step.effective_bdf.order == 2U &&
      !resumed_step.temporal_method_fallback &&
      resumed_step.proposal.generation ==
          continuous_step.proposal.generation &&
      continuous_endpoint.step == resumed_endpoint.step &&
      continuous_endpoint.time == resumed_endpoint.time &&
      continuous_endpoint.controller_state ==
          resumed_endpoint.controller_state &&
      maximum_relative_difference <= 1.0e-9;
  if (!passed) {
    std::cerr << std::setprecision(17)
              << "exact restart status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " image/atomic/restore=" << exact_image << '/'
              << failure_atomic << '/' << exact_restore
              << " orders=" << static_cast<unsigned>(second.effective_bdf.order)
              << '/'
              << static_cast<unsigned>(continuous_step.effective_bdf.order)
              << '/' << static_cast<unsigned>(resumed_step.proposal.bdf.order)
              << '/' << static_cast<unsigned>(resumed_step.effective_bdf.order)
              << " fallback=" << resumed_step.temporal_method_fallback
              << " max_relative=" << maximum_relative_difference << '\n';
  }
  fs::remove_all(restart_root, cleanup_error);
  return passed;
}

bool run_open_boundary_product() {
  ValidatedModel model = test::product_model({8, 7, 6});
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
  model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  model.boundaries[1U].allow_backflow = false;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.0, 0.0};
  if (status) status = driver.initialize(initial);
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  if (!status) {
    std::cerr << "open advance=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << " stage=" << step.failed_stage
              << " attempts=" << step.attempts << '\n';
    return false;
  }
  RestartSnapshot snapshot;
  const bool passed =
      step.proposal.bdf.order == 1U && step.effective_bdf.order == 1U &&
      terminal_physical_certificate(driver, model, step, &snapshot) &&
      open_x_boundary_flux_certificate(snapshot);
  if (!passed)
    std::cerr << "open semantic certificate attempts/origin/order="
              << step.attempts << '/'
              << static_cast<unsigned>(step.proposal.origin) << '/'
              << static_cast<unsigned>(step.proposal.bdf.order)
              << " residuals=" << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '\n';
  return passed;
}

bool fixed_advective_cfl_gate_rolls_back_product() {
  ValidatedModel model = test::product_model({8, 7, 6});
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
  model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  model.boundaries[1U].allow_backflow = false;
  model.time.control = TimeControlKind::fixed;
  model.time.convective_cfl = 1.0e-12;

  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  if (!status) {
    std::cerr << "fixed CFL fixture rejected before ProductDriver="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
    return false;
  }
  ProductDriver driver;
  status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.0, 0.0};
  if (status) status = driver.initialize(initial);
  if (!status) return false;

  const auto field_bytes = [](ConstFieldView field) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(field.interior.x) *
                   field.interior.y * field.interior.z * field.components);
    for (std::uint8_t component = 0U; component < field.components;
         ++component)
      for (std::int32_t z = 0; z < field.interior.z; ++z)
        for (std::int32_t y = 0; y < field.interior.y; ++y)
          for (std::int32_t x = 0; x < field.interior.x; ++x)
            values.push_back(field.unchecked({x, y, z}, component));
    return values;
  };
  const auto face_bytes = [](ConstFaceFieldView face) {
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(face.extents.x) *
                   face.extents.y * face.extents.z);
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          values.push_back(face.unchecked({x, y, z}));
    return values;
  };
  const auto exact_values = [](const std::vector<double>& left,
                               const std::vector<double>& right) {
    return left.size() == right.size() &&
           (left.empty() ||
            std::memcmp(left.data(), right.data(),
                        left.size() * sizeof(double)) == 0);
  };

  CommittedOutputSnapshot before;
  status = driver.committed_output_snapshot(before);
  ConstFaceFluxView before_flux;
  if (status) status = driver.committed_final_mass_flux_for_test(before_flux);
  if (!status || !before.committed || !before_flux.certificate.valid())
    return false;
  std::vector<SnapshotFieldView> before_fields(
      before.fields.data, before.fields.data + before.fields.size);
  std::vector<std::vector<double>> before_field_values;
  before_field_values.reserve(before_fields.size());
  for (const SnapshotFieldView& field : before_fields)
    before_field_values.push_back(field_bytes(field.values));
  const std::array<std::vector<double>, 3U> before_flux_values{{
      face_bytes(before_flux.x), face_bytes(before_flux.y),
      face_bytes(before_flux.z)}};
  const double pressure_reference_before = driver.pressure_reference();
  const double closed_mass_before = driver.closed_mass_target();

  DriverStepReport step;
  status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  const MomentumAdvectiveCflCertificate& cfl =
      step.momentum_predictor_limiter.advective_cfl;
  const bool advective_audit_completed =
      cfl.valid() && cfl.out_max > cfl.limit &&
      cfl.failure_witness.valid && cfl.failure_witness.rank == 0 &&
      cfl.failure_witness.out == cfl.out_max &&
      cfl.face_flux != 0U && step.piso.final_flux_revision == 0U;
  const bool fatal_gate =
      !status && status.code == StatusCode::invalid_case &&
      status.detail == 10211U && !step.accepted && step.attempts == 1U &&
      step.proposal.origin == StepOrigin::fresh_start &&
      step.proposal.attempt == 0U &&
      step.failure.code == StatusCode::invalid_case &&
      step.failure.detail == 10211U && step.failed_stage == 30U;

  CommittedOutputSnapshot after;
  const Status output_status = driver.committed_output_snapshot(after);
  ConstFaceFluxView after_flux;
  const Status flux_status =
      driver.committed_final_mass_flux_for_test(after_flux);
  bool exact_state = output_status && after.committed &&
                     after.time == before.time && after.step == before.step &&
                     after.plan == before.plan && after.schema == before.schema &&
                     after.fields.size == before_fields.size();
  for (std::size_t index = 0U;
       index < before_fields.size() && exact_state; ++index) {
    const SnapshotFieldView& expected = before_fields[index];
    const SnapshotFieldView& actual = after.fields.data[index];
    const ConstFieldView& left = expected.values;
    const ConstFieldView& right = actual.values;
    exact_state = actual.stable_name == expected.stable_name &&
                  actual.accepted_revision == expected.accepted_revision &&
                  right.base == left.base && right.interior.x == left.interior.x &&
                  right.interior.y == left.interior.y &&
                  right.interior.z == left.interior.z &&
                  right.ghosts.x == left.ghosts.x &&
                  right.ghosts.y == left.ghosts.y &&
                  right.ghosts.z == left.ghosts.z &&
                  right.components == left.components &&
                  right.stride_y == left.stride_y &&
                  right.stride_z == left.stride_z &&
                  right.component_stride == left.component_stride &&
                  right.replica == left.replica && right.field == left.field &&
                  right.revision == left.revision &&
                  right.storage_identity == left.storage_identity &&
                  right.revision_domain == left.revision_domain &&
                  exact_values(before_field_values[index],
                               field_bytes(right));
  }
  const bool exact_flux =
      flux_status && after_flux.revision == before_flux.revision &&
      after_flux.certificate == before_flux.certificate &&
      exact_values(before_flux_values[0U], face_bytes(after_flux.x)) &&
      exact_values(before_flux_values[1U], face_bytes(after_flux.y)) &&
      exact_values(before_flux_values[2U], face_bytes(after_flux.z));
  const bool exact_scalars =
      driver.pressure_reference() == pressure_reference_before &&
      driver.closed_mass_target() == closed_mass_before;
  const bool passed = advective_audit_completed && fatal_gate && exact_state &&
                      exact_flux && exact_scalars;
  if (!passed)
    std::cerr << "fixed advective CFL gate status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " failure=" << static_cast<unsigned>(step.failure.code)
              << '/' << step.failure.detail << " stage/attempts="
              << step.failed_stage << '/' << step.attempts << " audit/rev="
              << advective_audit_completed << '/' << cfl.face_flux
              << " final=" << step.piso.final_flux_revision << " CFL="
              << cfl.out_max << '/' << cfl.absolute_max << '/' << cfl.limit
              << " rollback=" << exact_state << '/' << exact_flux << '/'
              << exact_scalars << '\n';
  return passed;
}

bool convective_cfl_policy_classifies_terminal_failure() {
  constexpr double limit = 0.8;
  constexpr double slack =
      64.0 * std::numeric_limits<double>::epsilon();
  const double edge = limit * (1.0 + slack);
  const double over = std::nextafter(
      edge, std::numeric_limits<double>::infinity());
  const Status accepted = detail::convective_cfl_acceptance_status_for_test(
      TimeControlKind::fixed, edge, limit);
  const Status fixed = detail::convective_cfl_acceptance_status_for_test(
      TimeControlKind::fixed, over, limit);
  const Status adaptive_flow =
      detail::convective_cfl_acceptance_status_for_test(
          TimeControlKind::adaptive_flow, over, limit);
  const Status adaptive_acoustic =
      detail::convective_cfl_acceptance_status_for_test(
          TimeControlKind::adaptive_acoustic, over, limit);
  const Status invalid = detail::convective_cfl_acceptance_status_for_test(
      TimeControlKind::adaptive_flow,
      std::numeric_limits<double>::infinity(), limit);
  return accepted && fixed.code == StatusCode::invalid_case &&
         fixed.detail == 10211U &&
         adaptive_flow.code == StatusCode::rejected_step &&
         adaptive_flow.detail == 10211U &&
         adaptive_acoustic.code == StatusCode::rejected_step &&
         adaptive_acoustic.detail == 10211U &&
         invalid.code == StatusCode::invalid_plan && invalid.detail == 10211U;
}

ValidatedModel mass_flow_boundary_model() {
  ValidatedModel model = test::product_model({8, 7, 6});
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  model.boundaries[0U].flow_kind = BoundaryKind::mass_flow_inlet;
  model.boundaries[0U].direction = {1.0, 0.0, 0.0};
  model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  model.boundaries[0U].mass_flow_rate = 1.0;
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  return model;
}

bool terminal_committed_cfl_gate_is_transactional() {
  constexpr double full_dt = 1.0e-3;
  constexpr double half_dt = 5.0e-4;
  constexpr double cfl_limit = 3.0e-3;
  const auto configured_model = [](TimeControlKind control, double initial_dt,
                                   double minimum_dt, double maximum_dt) {
    ValidatedModel model = mass_flow_boundary_model();
    model.time.control = control;
    model.time.initial_dt = initial_dt;
    model.time.minimum_dt = minimum_dt;
    model.time.maximum_dt = maximum_dt;
    model.time.convective_cfl = cfl_limit;
    model.time.maximum_growth = 1.0;
    model.time.retry_factor = 0.5;
    model.time.maximum_retries = 1U;
    return model;
  };
  const auto initialize = [](const ValidatedModel& model,
                             ProductDriver& driver) {
    CompiledCasePlan plan;
    Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
    if (status)
      status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
    DriverInitialState initial;
    initial.pressure_reference = 98000.0;
    initial.temperature = 315.0;
    if (status) status = driver.initialize(initial);
    return status;
  };
  const auto terminal_components_pass = [](const ValidatedModel& model,
                                           const DriverStepReport& report) {
    const auto accepted = [](double value, double limit) {
      return std::isfinite(value) && value >= 0.0 && value <= limit;
    };
    return report.piso.final_flux_revision != 0U &&
           accepted(report.piso.eos_residual, model.solver.terminal.eos) &&
           accepted(report.piso.continuity_residual,
                    model.solver.terminal.continuity) &&
           accepted(report.piso.energy_residual,
                    model.solver.terminal.continuity) &&
           accepted(report.piso.closed_mass_residual,
                    model.solver.terminal.closed_mass) &&
           accepted(report.piso.gauge_residual,
                    model.solver.terminal.gauge);
  };

  const ValidatedModel fixed_model = configured_model(
      TimeControlKind::fixed, full_dt, half_dt, full_dt);
  ProductDriver fixed_driver;
  Status status = initialize(fixed_model, fixed_driver);
  DriverStateImage fixed_before;
  if (!status || !capture_driver_state(fixed_driver, fixed_before))
    return false;
  DriverStepReport fixed_report;
  status = fixed_driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, fixed_report);
  DriverStateImage fixed_after;
  const bool fixed_snapshot = capture_driver_state(fixed_driver, fixed_after);
  const MomentumAdvectiveCflCertificate& fixed_advective =
      fixed_report.momentum_predictor_limiter.advective_cfl;
  const CommittedConvectiveCflCertificate& fixed_committed =
      fixed_report.piso.committed_convective_cfl;
  const bool fixed_committed_authority =
      fixed_committed.valid() &&
      fixed_committed.final_flux == fixed_report.piso.final_flux_revision &&
      fixed_committed.density != 0U &&
      fixed_committed.density_view_identity.valid() &&
      fixed_committed.face_flux_view_identity.valid() &&
      fixed_committed.dt == fixed_report.proposal.dt &&
      fixed_committed.out_max ==
          fixed_report.piso.committed_convective_cfl_out_max &&
      fixed_committed.absolute_max ==
          fixed_report.piso.committed_convective_cfl_abs_max &&
      fixed_committed.limit ==
          fixed_report.piso.committed_convective_cfl_limit &&
      fixed_committed.out_winner.valid &&
      fixed_committed.absolute_winner.valid &&
      fixed_committed.out_winner.rank == 0 &&
      fixed_committed.absolute_winner.rank == 0 &&
      fixed_committed.out_winner.out == fixed_committed.out_max &&
      fixed_committed.absolute_winner.absolute ==
          fixed_committed.absolute_max &&
      fixed_committed.out_winner.out ==
          (fixed_committed.dt /
           fixed_committed.out_winner.density_volume) *
              fixed_committed.out_winner.outgoing_mass_flow &&
      fixed_committed.absolute_winner.absolute ==
          (fixed_committed.dt /
           fixed_committed.absolute_winner.density_volume) *
              (0.5 *
               fixed_committed.absolute_winner.absolute_mass_flow) &&
      fixed_committed.failure_witness.valid &&
      fixed_committed.failure_witness.out_winner.out ==
          fixed_committed.out_max &&
      fixed_committed.failure_witness.absolute_winner.absolute ==
          fixed_committed.absolute_max;
  CommittedConvectiveCflCertificate wrong_density = fixed_committed;
  ++wrong_density.density;
  CommittedConvectiveCflCertificate wrong_flux = fixed_committed;
  ++wrong_flux.final_flux;
  CommittedConvectiveCflCertificate wrong_dt = fixed_committed;
  wrong_dt.dt *= 0.5;
  CommittedConvectiveCflCertificate wrong_out_winner = fixed_committed;
  wrong_out_winner.out_winner.out = std::nextafter(
      wrong_out_winner.out_winner.out,
      std::numeric_limits<double>::infinity());
  CommittedConvectiveCflCertificate wrong_abs_winner = fixed_committed;
  wrong_abs_winner.absolute_winner.absolute_mass_flow = std::nextafter(
      wrong_abs_winner.absolute_winner.absolute_mass_flow,
      std::numeric_limits<double>::infinity());
  CommittedConvectiveCflCertificate missing_failure = fixed_committed;
  missing_failure.failure_witness = {};
  const bool committed_mutations_rejected =
      !wrong_density.valid() && !wrong_flux.valid() && !wrong_dt.valid() &&
      !wrong_out_winner.valid() && !wrong_abs_winner.valid() &&
      !missing_failure.valid();
  const bool fixed_passed =
      !status && status.code == StatusCode::invalid_case &&
      status.detail == 10211U && !fixed_report.accepted &&
      fixed_report.attempts == 1U && fixed_report.failed_stage == 60U &&
      fixed_report.failure.code == StatusCode::invalid_case &&
      fixed_report.failure.detail == 10211U && fixed_advective.valid() &&
      fixed_advective.out_max <= fixed_advective.limit &&
      terminal_components_pass(fixed_model, fixed_report) &&
      fixed_report.piso.committed_convective_cfl_out_max >
          fixed_report.piso.committed_convective_cfl_limit &&
      fixed_committed_authority && committed_mutations_rejected &&
      fixed_snapshot && same_exact_state(fixed_before, fixed_after);
  if (!fixed_passed) {
    std::cerr << std::setprecision(17)
              << "terminal fixed CFL status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " failure/stage/attempts="
              << static_cast<unsigned>(fixed_report.failure.code) << '/'
              << fixed_report.failure.detail << '/' << fixed_report.failed_stage
              << '/' << fixed_report.attempts << " provisional/terminal="
              << fixed_advective.out_max << '/' << fixed_advective.limit << ' '
              << fixed_report.piso.committed_convective_cfl_out_max << '/'
              << fixed_report.piso.committed_convective_cfl_limit
              << " rollback=" << fixed_snapshot << '/'
              << (fixed_snapshot && same_exact_state(fixed_before, fixed_after))
              << " authority/mutations=" << fixed_committed_authority << '/'
              << committed_mutations_rejected
              << " cert=" << fixed_committed.valid() << ':'
              << fixed_committed.density << ':' << fixed_committed.final_flux
              << ':' << fixed_committed.dt << ':'
              << fixed_committed.out_winner.formula_valid(fixed_committed.dt)
              << ':'
              << fixed_committed.absolute_winner.formula_valid(
                     fixed_committed.dt)
              << ':' << fixed_committed.failure_witness.valid
              << '\n';
    return false;
  }

  const ValidatedModel retry_model = configured_model(
      TimeControlKind::adaptive_flow, full_dt, half_dt, full_dt);
  ProductDriver retry_driver;
  status = initialize(retry_model, retry_driver);
  DriverStepReport retry_report;
  if (status)
    status = retry_driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, retry_report);
  DriverStateImage retry_state;
  const bool retry_snapshot = status &&
                              capture_driver_state(retry_driver, retry_state);
  const MomentumAdvectiveCflCertificate& retry_advective =
      retry_report.momentum_predictor_limiter.advective_cfl;
  const CommittedConvectiveCflCertificate& retry_committed =
      retry_report.piso.committed_convective_cfl;
  const bool retry_passed =
      retry_snapshot && retry_report.accepted && retry_report.attempts == 2U &&
      retry_report.proposal.origin == StepOrigin::retry &&
      retry_report.proposal.attempt == 1U &&
      retry_report.proposal.dt == half_dt &&
      retry_report.failure.code == StatusCode::rejected_step &&
      retry_report.failure.detail == 10211U &&
      retry_report.failed_stage == 60U && retry_advective.valid() &&
      retry_advective.out_max <= retry_advective.limit &&
      retry_report.piso.committed_convective_cfl_out_max <=
          retry_report.piso.committed_convective_cfl_limit &&
      retry_committed.valid() &&
      retry_committed.final_flux == retry_report.piso.final_flux_revision &&
      retry_committed.dt == retry_report.proposal.dt &&
      retry_committed.out_winner.valid &&
      retry_committed.absolute_winner.valid &&
      !retry_committed.failure_witness.valid &&
      terminal_physical_certificate(retry_driver, retry_model, retry_report);
  if (!retry_passed) {
    std::cerr << std::setprecision(17)
              << "terminal adaptive CFL status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/attempts/origin=" << retry_report.accepted << '/'
              << retry_report.attempts << '/'
              << static_cast<unsigned>(retry_report.proposal.origin)
              << " failure/stage="
              << static_cast<unsigned>(retry_report.failure.code) << '/'
              << retry_report.failure.detail << '/' << retry_report.failed_stage
              << " dt=" << retry_report.proposal.dt << " CFL="
              << retry_advective.out_max << '/' << retry_advective.limit << ' '
              << retry_report.piso.committed_convective_cfl_out_max << '/'
              << retry_report.piso.committed_convective_cfl_limit << '\n';
    return false;
  }

  const ValidatedModel half_model = configured_model(
      TimeControlKind::adaptive_flow, half_dt, half_dt, half_dt);
  ProductDriver half_driver;
  status = initialize(half_model, half_driver);
  DriverStepReport half_report;
  if (status)
    status = half_driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, half_report);
  DriverStateImage half_state;
  const bool half_snapshot = status &&
                             capture_driver_state(half_driver, half_state);
  const bool half_passed =
      half_snapshot && half_report.accepted && half_report.attempts == 1U &&
      half_report.proposal.origin == StepOrigin::fresh_start &&
      half_report.proposal.dt == half_dt &&
      terminal_physical_certificate(half_driver, half_model, half_report) &&
      same_physical_state(retry_state, half_state);
  if (!half_passed)
    std::cerr << std::setprecision(17)
              << "terminal direct-half status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/attempts/dt=" << half_report.accepted << '/'
              << half_report.attempts << '/' << half_report.proposal.dt
              << " retry/direct state="
              << (half_snapshot && same_physical_state(retry_state, half_state))
              << '\n';
  return half_passed;
}

bool run_mass_flow_boundary_product() {
  ValidatedModel model = mass_flow_boundary_model();
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport step;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  if (!status) {
    std::cerr << "mass-flow status=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << " stage=" << step.failed_stage
              << " attempts=" << step.attempts << '\n';
    return false;
  }
  RestartSnapshot snapshot;
  const bool passed =
      step.proposal.bdf.order == 1U && step.effective_bdf.order == 1U &&
      terminal_physical_certificate(driver, model, step, &snapshot) &&
      open_x_boundary_flux_certificate(
          snapshot, model.boundaries[0U].mass_flow_rate);
  if (!passed)
    std::cerr << "mass-flow semantic certificate attempts/origin/order="
              << step.attempts << '/'
              << static_cast<unsigned>(step.proposal.origin) << '/'
              << static_cast<unsigned>(step.proposal.bdf.order)
              << " residuals=" << step.piso.eos_residual << '/'
              << step.piso.continuity_residual << '/'
              << step.piso.energy_residual << '\n';
  return passed;
}

bool retry_consumes_warm_seed_and_restart_starts_cold() {
  ValidatedModel base_model = test::product_model({17, 11, 7});
  base_model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : base_model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  base_model.boundaries[0U].flow_kind = BoundaryKind::mass_flow_inlet;
  base_model.boundaries[0U].direction = {1.0, 0.0, 0.0};
  base_model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  base_model.boundaries[0U].mass_flow_rate = 0.25;
  base_model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;

  // The strict linear-work policy makes the initial proposal recover through
  // one or more BE retries.  The exact count and rejected residual are not a
  // contract; only rollback, retry origin, final acceptance, and physics are.
  ValidatedModel retry_model = base_model;
  retry_model.solver.pressure.absolute_tolerance = 1.0e-16;
  retry_model.solver.pressure.relative_tolerance = 1.0e-16;
  retry_model.solver.pressure.maximum_iterations = 24U;
  retry_model.time.maximum_growth = 2.0;

  CompiledCasePlan retry_plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, retry_model, {},
                                            retry_plan);
  ProductDriver retry_driver;
  if (status) {
    status = ProductDriver::create(MPI_COMM_SELF, std::move(retry_plan),
                                   retry_driver);
  }
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = retry_driver.initialize(initial);
  DriverStepReport retried;
  if (status)
    status = retry_driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, retried);
  RestartSnapshot retry_restart;
  const bool retry_physical =
      status && terminal_physical_certificate(
                    retry_driver, retry_model, retried, &retry_restart);
  const bool retry_certificate =
      retry_physical && retried.attempts > 1U &&
      retried.proposal.origin == StepOrigin::retry &&
      retried.proposal.attempt > 0U && retried.proposal.bdf.order == 1U &&
      retried.effective_bdf.order == 1U &&
      open_x_boundary_flux_certificate(
          retry_restart, retry_model.boundaries[0U].mass_flow_rate);
  if (!retry_certificate) {
    std::cerr << "warm retry fixture first="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/attempts=" << retried.accepted << '/'
              << retried.attempts << " origin/attempt/order="
              << static_cast<unsigned>(retried.proposal.origin) << '/'
              << retried.proposal.attempt << '/'
              << static_cast<unsigned>(retried.proposal.bdf.order)
              << " residuals=" << retried.piso.eos_residual << '/'
              << retried.piso.continuity_residual << '/'
              << retried.piso.energy_residual << '\n';
    return false;
  }

  // Warm-start provenance is an independent lifecycle contract.  Use a stable
  // Product budget so an accepted-origin BDF2 step is observed directly;
  // coupling this observation to the intentionally strict retry fixture would
  // only test its work limit a second time.
  ValidatedModel authority_model = test::product_model({17, 11, 7});
  authority_model.time.initial_dt = 1.0e-5;
  authority_model.time.maximum_growth = 1.0;
  CompiledCasePlan live_plan;
  status = ProductCompiler::compile(MPI_COMM_SELF, authority_model, {},
                                    live_plan);
  ProductDriver live;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(live_plan), live);
  if (status) status = live.initialize(initial);
  DriverStepReport branch_step;
  if (status)
    status = live.advance({1.0, 1.0, 1.0, 1.0, 1.0}, branch_step);
  RestartSnapshot branch_restart;
  const bool branch_physical =
      status && terminal_physical_certificate(
                    live, authority_model, branch_step, &branch_restart) &&
      branch_step.proposal.bdf.order == 1U &&
      branch_step.effective_bdf.order == 1U;
  if (!branch_physical) {
    std::cerr << "warm authority branch="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/attempts/origin/order="
              << branch_step.accepted << '/' << branch_step.attempts << '/'
              << static_cast<unsigned>(branch_step.proposal.origin) << '/'
              << static_cast<unsigned>(branch_step.proposal.bdf.order)
              << " residuals=" << branch_step.piso.eos_residual << '/'
              << branch_step.piso.continuity_residual << '/'
              << branch_step.piso.energy_residual << '\n';
    return false;
  }

  CompiledCasePlan cold_plan;
  status = ProductCompiler::compile(MPI_COMM_SELF, authority_model, {},
                                    cold_plan);
  ProductDriver cold;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(cold_plan), cold);
  RestartExpected expected;
  if (status) status = cold.restart_expected(expected);
  RestartImage image;
  if (status && !copy_restart_image(branch_restart, expected, image))
    status = {StatusCode::invalid_plan, 0U};
  if (status) status = cold.initialize_restart(image);
  if (!status) {
    std::cerr << "warm retry restart setup="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
    return false;
  }

  const double common_target_dt = branch_step.proposal.dt;
  const LocalTimeLimits same_target_limits =
      limits_for_exact_dt(authority_model, common_target_dt);
  DriverStepReport live_step;
  status = live.advance(same_target_limits, live_step);
  if (!status) {
    std::cerr << "warm authority live=" << static_cast<unsigned>(status.code)
              << '/' << status.detail << " stage=" << live_step.failed_stage
              << " attempts=" << live_step.attempts << '\n';
    return false;
  }
  detail::PressureCorrectionWarmStartDiagnostic live_warm_start;
  const bool observed_live_warm_start =
      detail::pressure_correction_warm_start_diagnostic_for_test(
          live_warm_start);

  DriverStepReport cold_step;
  status = cold.advance(same_target_limits, cold_step);
  if (!status) {
    std::cerr << "warm retry cold=" << static_cast<unsigned>(status.code)
              << '/' << status.detail << " stage=" << cold_step.failed_stage
              << " attempts=" << cold_step.attempts << '\n';
    return false;
  }
  detail::PressureCorrectionWarmStartDiagnostic cold_warm_start;
  const bool observed_cold_warm_start =
      detail::pressure_correction_warm_start_diagnostic_for_test(
          cold_warm_start);

  RestartSnapshot live_restart;
  RestartSnapshot cold_restart;
  const bool live_physical = terminal_physical_certificate(
      live, authority_model, live_step, &live_restart);
  const bool cold_physical = terminal_physical_certificate(
      cold, authority_model, cold_step, &cold_restart);

  const bool passed =
      live_physical && live_step.proposal.origin == StepOrigin::accepted &&
      live_step.proposal.attempt == 0U &&
      live_step.proposal.bdf.order == 2U &&
      live_step.effective_bdf.order == 2U && observed_live_warm_start &&
      live_warm_start.valid &&
      live_warm_start.origin == StepOrigin::accepted &&
      live_warm_start.attempt == 0U &&
      live_warm_start.authority_available && live_warm_start.used &&
      cold_physical &&
      cold_step.proposal.origin == StepOrigin::restart &&
      cold_step.proposal.attempt == 0U &&
      cold_step.proposal.bdf.order == 1U &&
      cold_step.effective_bdf.order == 1U && observed_cold_warm_start &&
      cold_warm_start.valid &&
      cold_warm_start.origin == StepOrigin::restart &&
      cold_warm_start.attempt == 0U &&
      !cold_warm_start.authority_available && !cold_warm_start.used &&
      live_step.proposal.time == cold_step.proposal.time &&
      live_step.proposal.dt == cold_step.proposal.dt &&
      live_step.accepted_step == cold_step.accepted_step &&
      live_step.accepted_time == cold_step.accepted_time;
  if (!passed) {
    std::cerr << std::setprecision(17)
              << "warm retry comparison status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " attempts=" << live_step.attempts << '/'
              << cold_step.attempts << " failure="
              << static_cast<unsigned>(live_step.failure.code) << '/'
              << live_step.failure.detail << '/' << live_step.failed_stage
              << " origin=" << static_cast<unsigned>(live_step.proposal.origin)
              << '/' << static_cast<unsigned>(cold_step.proposal.origin)
              << " order="
              << static_cast<unsigned>(live_step.proposal.bdf.order) << '/'
              << static_cast<unsigned>(cold_step.proposal.bdf.order)
              << " dt/time=" << live_step.proposal.dt << '/'
              << cold_step.proposal.dt << ' ' << live_step.accepted_time << '/'
              << cold_step.accepted_time << " warm="
              << observed_live_warm_start << '/' << live_warm_start.valid
              << '/' << live_warm_start.authority_available << '/'
              << live_warm_start.used << " cold="
              << observed_cold_warm_start << '/' << cold_warm_start.valid
              << '/' << cold_warm_start.authority_available << '/'
              << cold_warm_start.used << " physical=" << live_physical << '/'
              << cold_physical << " residuals="
              << live_step.piso.eos_residual << '/'
              << live_step.piso.continuity_residual << '/'
              << live_step.piso.energy_residual << ' '
              << cold_step.piso.eos_residual << '/'
              << cold_step.piso.continuity_residual << '/'
              << cold_step.piso.energy_residual << '\n';
  }
  return passed;
}

bool run_immersed_final_force_product() {
  ValidatedModel model = test::product_model({16, 16, 16});
  model.mesh.lower = {-2.0, -2.0, -2.0};
  model.mesh.upper = {2.0, 2.0, 2.0};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  // Keep the product integration focused on lifecycle/authority at a stable
  // coarse-mesh target step; retry count is an adaptive-controller detail.
  model.time.initial_dt = 2.5e-4;
  model.immersed_boundary = ImmersedBoundarySpec{
      "cylinder_ascii.stl", ImmersedFluidSide::outside};

  CompiledCasePlan plan;
  const fs::path data_root =
      fs::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, data_root,
                                            plan);
  if (!status) {
    std::cerr << "immersed compile=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << '\n';
    return false;
  }
  ProductDriver driver;
  status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = driver.initialize(initial);
  if (!status) {
    std::cerr << "immersed initialize="
              << static_cast<unsigned>(status.code) << "/" << status.detail
              << '\n';
    return false;
  }

  SurfaceForce force;
  FinalForceCertificate certificate;
  status = driver.committed_surface_force(force, certificate);
  if (status || status.code != StatusCode::invalid_plan ||
      certificate.valid()) {
    std::cerr << "immersed pre-step force unexpectedly available\n";
    return false;
  }

  DriverStepReport step;
  status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  if (!status) {
    std::cerr << "immersed advance=" << static_cast<unsigned>(status.code)
              << "/" << status.detail << " stage=" << step.failed_stage
              << " attempts=" << step.attempts << '\n';
    return false;
  }
  status = driver.committed_surface_force(force, certificate);
  if (!status) {
    std::cerr << "immersed committed force="
              << static_cast<unsigned>(status.code) << "/" << status.detail
              << '\n';
  }
  RestartSnapshot snapshot;
  const bool terminal_physical =
      status && terminal_physical_certificate(driver, model, step, &snapshot);
  bool passed =
      terminal_physical && certificate.valid() &&
      std::isfinite(force.pressure.x) && std::isfinite(force.pressure.y) &&
      std::isfinite(force.pressure.z) && std::isfinite(force.viscous.x) &&
      std::isfinite(force.viscous.y) && std::isfinite(force.viscous.z);
  if (!passed)
    std::cerr << "immersed final accepted=" << step.accepted
              << " attempts=" << step.attempts
              << " prior_failure="
              << static_cast<unsigned>(step.failure.code) << '/'
              << step.failure.detail << " prior_stage=" << step.failed_stage
              << " certificate=" << certificate.valid()
              << " pressure=" << force.pressure.x << ',' << force.pressure.y
              << ',' << force.pressure.z << " viscous=" << force.viscous.x
              << ',' << force.viscous.y << ',' << force.viscous.z << '\n';
  if (!passed) return false;

  CompiledCasePlan target_plan;
  if (status)
    status = ProductCompiler::compile(MPI_COMM_SELF, model, data_root,
                                      target_plan);
  ProductDriver target;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(target_plan),
                                   target);
  RestartExpected expected;
  if (status) status = target.restart_expected(expected);
  RestartImage valid;
  if (status && !copy_restart_image(snapshot, expected, valid))
    status = {StatusCode::invalid_plan, 0U};
  if (!status) return false;
  RestartImage invalid = valid;
  for (std::vector<double>& face : invalid.final_mass_flux)
    std::fill(face.begin(), face.end(), 1.0);
  status = target.initialize_restart(invalid);
  passed &= !status && status.code == StatusCode::numerical_failure &&
            !target.initialized() && target.pressure_reference() == 0.0 &&
            target.closed_mass_target() == 0.0;
  if (!passed) {
    std::cerr << "immersed invalid Restart flux published state status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " initialized=" << target.initialized() << '\n';
    return false;
  }
  status = target.initialize_restart(valid);
  passed &= static_cast<bool>(status) && target.initialized();
  return passed;
}

bool reject_unresolved_characteristic_boundary() {
  ValidatedModel model = test::product_model({8, 7, 6});
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
    face.pressure = 98000.0;
    face.temperature = 315.0;
  }
  model.boundaries[0U].flow_kind = BoundaryKind::static_state_inlet;
  model.boundaries[0U].direction = {1.0, 0.0, 0.0};
  model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
  model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  if (!status) {
    const bool rejected = status.code == StatusCode::invalid_plan &&
                          status.detail != 0U && plan.fingerprint() == 0U;
    if (!rejected)
      std::cerr << "unresolved characteristic compile="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << " plan=" << plan.fingerprint() << '\n';
    return rejected;
  }
  ProductDriver driver;
  status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  if (!status) {
    const bool rejected = status.code == StatusCode::invalid_plan &&
                          status.detail != 0U && !driver.initialized();
    if (!rejected)
      std::cerr << "unresolved characteristic create="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << " initialized=" << driver.initialized() << '\n';
    return rejected;
  }
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  initial.velocity = {1.0, 0.0, 0.0};
  status = driver.initialize(initial);
  if (!status) {
    const bool rejected = status.code == StatusCode::invalid_plan &&
                          status.detail != 0U && !driver.initialized() &&
                          driver.pressure_reference() == 0.0 &&
                          driver.closed_mass_target() == 0.0;
    if (!rejected)
      std::cerr << "unresolved characteristic initialize="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << " initialized=" << driver.initialized() << '\n';
    return rejected;
  }
  CommittedOutputSnapshot before;
  status = driver.committed_output_snapshot(before);
  if (!status) return false;
  DriverStepReport step;
  status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
  if (status || status.code != StatusCode::invalid_plan || step.accepted) {
    std::cerr << "unresolved characteristic status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " accepted/stage/attempts=" << step.accepted << '/'
              << step.failed_stage << '/' << step.attempts << '\n';
    return false;
  }
  CommittedOutputSnapshot after;
  const Status snapshot = driver.committed_output_snapshot(after);
  const bool rolled_back = snapshot && before.step == 0U && after.step == 0U &&
                           before.time == after.time;
  if (!rolled_back)
    std::cerr << "unresolved characteristic rollback="
              << static_cast<unsigned>(snapshot.code) << '/'
              << snapshot.detail << " step=" << before.step << '/'
              << after.step << " time=" << before.time << '/' << after.time
              << '\n';
  return rolled_back;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  bool passed = true;
  const auto check = [&](const char* name, bool result) {
    if (!result) std::cerr << "driver fixture failed: " << name << '\n';
    passed &= result;
  };
  check("shared-cell-CFL", shared_cell_convective_cfl_kernel_has_one_formula());
  check("cold-plan", run());
  check("initialized-product", run_initialized_product());
  check("initialized-species", run_initialized_species_product());
  check("invalid-restart-rollback",
        reject_invalid_thermodynamic_restart_atomically());
  check("restart-flux-lineage", restart_preserves_direct_flux_lineage());
  check("exact-restart-BDF2", exact_restart_preserves_bdf2_history());
  check("open-boundary", run_open_boundary_product());
  check("fixed-provisional-CFL",
        fixed_advective_cfl_gate_rolls_back_product());
  check("terminal-CFL-policy",
        convective_cfl_policy_classifies_terminal_failure());
  check("terminal-CFL-transaction",
        terminal_committed_cfl_gate_is_transactional());
  check("mass-flow-boundary", run_mass_flow_boundary_product());
  check("immersed-force", run_immersed_final_force_product());
  check("warm-seed-retry-restart",
        retry_consumes_warm_seed_and_restart_starts_cold());
  check("unresolved-characteristic",
        reject_unresolved_characteristic_boundary());
  if (!passed) std::cerr << "driver cold-plan failure\n";
  MPI_Finalize();
  return passed ? 0 : 1;
}
