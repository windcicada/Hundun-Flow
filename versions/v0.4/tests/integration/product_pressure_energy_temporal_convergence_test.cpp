// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"

#include "../support/product_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr Int3 kCells{8, 4, 4};
constexpr double kPressureReference = 101325.0;
constexpr double kBaseTemperature = 320.0;
constexpr double kPressureAmplitude = 50.0;
constexpr double kStreamwiseVelocity = 0.50;
constexpr double kCrossflowAmplitude = 0.20;
constexpr double kCoarseDt = 1.25e-4;
constexpr double kEvolutionTime = 1.0e-3;
constexpr double kPressureLinearTolerance = 1.0e-13;
constexpr std::array<double, 3U> kRefinement{{1.0, 2.0, 4.0}};
constexpr double kMolecularWeight = 28.96546;
constexpr double kGasConstant = kUniversalGasConstant / kMolecularWeight;
constexpr double kHeatCapacity = 3.5 * kGasConstant;
constexpr double kDomainX = 2.0;
constexpr double kDomainY = 1.0;
constexpr double kDomainZ = 0.5;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool close_coefficient(double actual, double expected) noexcept {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             32.0 * std::numeric_limits<double>::epsilon() * scale;
}

bool backward_euler_coefficients(BdfCoefficients bdf, double dt) noexcept {
  return bdf.order == 1U && close_coefficient(bdf.a0, 1.0 / dt) &&
         close_coefficient(bdf.a1, -1.0 / dt) &&
         close_coefficient(bdf.a2, 0.0);
}

bool constant_bdf2_coefficients(BdfCoefficients bdf, double dt) noexcept {
  return bdf.order == 2U && close_coefficient(bdf.a0, 1.5 / dt) &&
         close_coefficient(bdf.a1, -2.0 / dt) &&
         close_coefficient(bdf.a2, 0.5 / dt);
}

std::size_t cell_offset(Int3 cells, Int3 cell,
                        std::uint8_t components = 1U,
                        std::uint8_t component = 0U) noexcept {
  return (static_cast<std::size_t>(cell.x) +
          static_cast<std::size_t>(cells.x) *
              (static_cast<std::size_t>(cell.y) +
               static_cast<std::size_t>(cells.y) *
                   static_cast<std::size_t>(cell.z))) *
             components +
         component;
}

std::size_t face_offset(Int3 extents, Int3 face) noexcept {
  return static_cast<std::size_t>(face.x) +
         static_cast<std::size_t>(extents.x) *
             (static_cast<std::size_t>(face.y) +
              static_cast<std::size_t>(extents.y) *
                  static_cast<std::size_t>(face.z));
}

ValidatedModel temporal_model() {
  ValidatedModel model = test::product_model(kCells);
  model.fingerprint = UINT64_C(0x18000c201);
  model.turbulence = TurbulenceKind::none;
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.time.initial_dt = kCoarseDt;
  model.time.minimum_dt = kCoarseDt / 16.0;
  model.time.maximum_dt = kCoarseDt;
  model.time.convective_cfl = 1.0;
  model.time.viscous_cfl = 1.0;
  model.time.thermal_cfl = 1.0;
  model.time.species_cfl = 1.0;
  model.time.acoustic_cfl = 1.0;
  model.time.maximum_growth = 4.0;
  model.time.minimum_bdf_ratio = 0.2;
  model.time.maximum_bdf_ratio = 5.0;
  model.time.maximum_retries = 1U;
  model.solver.pressure.absolute_tolerance = kPressureLinearTolerance;
  model.solver.pressure.relative_tolerance = kPressureLinearTolerance;
  model.solver.pressure.maximum_iterations = 800U;
  model.solver.pressure.true_residual_interval = 4U;
  model.solver.pressure.krylov_restart = 64U;
  model.solver.terminal.eos = 1.0e-12;
  model.solver.terminal.continuity = 1.0e-12;
  model.solver.terminal.closed_mass = 1.0e-12;
  model.solver.terminal.gauge = 1.0e-12;
  model.thermophysics.species[0U].viscosity_reference = 1.0e-5;
  model.thermophysics.species[0U].conductivity = 1.0e-2;
  // Match the internal closed-mass Newton stopping criterion to the terminal
  // mass gate, avoiding a stricter double-roundoff stagnation criterion while
  // leaving the accepted-state physical gate unchanged.
  model.thermophysics.closed_mass_relative_tolerance =
      model.solver.terminal.closed_mass;
  return model;
}

RestartImage compatible_wave(const RestartExpected& expected,
                             double restart_dt) {
  RestartImage image;
  image.global_cells = expected.global_cells;
  image.patch = expected.target_patch;
  image.plan = expected.plan;
  image.schema = expected.schema;
  image.geometry = expected.geometry;
  image.time = 0.0;
  image.dt = restart_dt;
  image.pressure_reference = kPressureReference;
  image.step = 4U;
  image.controller_state = 1U;
  image.backward_euler_recovery = true;

  const Int3 local = image.patch.cells;
  const std::size_t local_cell_count =
      static_cast<std::size_t>(local.x) * local.y * local.z;
  image.fields.reserve(expected.fields.size);
  for (std::size_t field_index = 0U; field_index < expected.fields.size;
       ++field_index) {
    const RestartExpectedField descriptor = expected.fields.data[field_index];
    RestartImageField field;
    field.role = descriptor.role;
    field.field = descriptor.field;
    field.components = descriptor.components;
    field.values.assign(local_cell_count * descriptor.components, 0.0);
    image.fields.push_back(std::move(field));
  }

  const double rho = kPressureReference / (kGasConstant * kBaseTemperature);
  constexpr double gamma = 1.4;
  const double sound_speed =
      std::sqrt(gamma * kGasConstant * kBaseTemperature);
  for (RestartImageField& field : image.fields) {
    for (std::int32_t z = 0; z < local.z; ++z) {
      for (std::int32_t y = 0; y < local.y; ++y) {
        for (std::int32_t x = 0; x < local.x; ++x) {
          const double global_x =
              (static_cast<double>(image.patch.begin.x + x) + 0.5) *
              kDomainX / static_cast<double>(kCells.x);
          const double phase = 2.0 * std::acos(-1.0) * global_x / kDomainX;
          const double pressure_perturbation =
              kPressureAmplitude * std::sin(phase);
          const double absolute_pressure =
              kPressureReference + pressure_perturbation;
          const double pressure_ratio =
              absolute_pressure / kPressureReference;
          const double temperature =
              kBaseTemperature *
              std::pow(pressure_ratio, (gamma - 1.0) / gamma);
          const Int3 cell{x, y, z};
          const std::size_t offset =
              cell_offset(local, cell, field.components);
          if (field.role == RestartFieldRole::pressure_perturbation)
            field.values[offset] = pressure_perturbation;
          else if (field.role == RestartFieldRole::pressure_absolute)
            field.values[offset] = absolute_pressure;
          else if (field.role == RestartFieldRole::enthalpy)
            field.values[offset] = kHeatCapacity * temperature;
          else if (field.role == RestartFieldRole::velocity) {
            field.values[offset] =
                kStreamwiseVelocity +
                pressure_perturbation / (rho * sound_speed);
            field.values[offset + 1U] =
                kCrossflowAmplitude * std::sin(phase);
            field.values[offset + 2U] = 0.0;
          }
        }
      }
    }
  }

  const Int3 x_faces{local.x + 1, local.y, local.z};
  const Int3 y_faces{local.x, local.y + 1, local.z};
  const Int3 z_faces{local.x, local.y, local.z + 1};
  image.final_mass_flux[0U].assign(
      static_cast<std::size_t>(x_faces.x) * x_faces.y * x_faces.z,
      0.0);
  image.final_mass_flux[1U].assign(
      static_cast<std::size_t>(y_faces.x) * y_faces.y * y_faces.z, 0.0);
  image.final_mass_flux[2U].assign(
      static_cast<std::size_t>(z_faces.x) * z_faces.y * z_faces.z, 0.0);
  const double x_area = (kDomainY / static_cast<double>(kCells.y)) *
                        (kDomainZ / static_cast<double>(kCells.z));
  const auto streamwise_momentum = [&](std::int32_t global_cell_x) {
    global_cell_x %= kCells.x;
    if (global_cell_x < 0) global_cell_x += kCells.x;
    const double global_x =
        (static_cast<double>(global_cell_x) + 0.5) * kDomainX /
        static_cast<double>(kCells.x);
    const double phase = 2.0 * std::acos(-1.0) * global_x / kDomainX;
    const double pressure_perturbation = kPressureAmplitude * std::sin(phase);
    const double absolute_pressure =
        kPressureReference + pressure_perturbation;
    const double temperature =
        kBaseTemperature *
        std::pow(absolute_pressure / kPressureReference,
                 (gamma - 1.0) / gamma);
    const double local_density =
        absolute_pressure / (kGasConstant * temperature);
    const double velocity =
        kStreamwiseVelocity + pressure_perturbation / (rho * sound_speed);
    return local_density * velocity;
  };
  for (std::int32_t z = 0; z < x_faces.z; ++z) {
    for (std::int32_t y = 0; y < x_faces.y; ++y) {
      for (std::int32_t x = 0; x < x_faces.x; ++x) {
        const std::int32_t right = image.patch.begin.x + x;
        const std::int32_t left = right - 1;
        image.final_mass_flux[0U][face_offset(x_faces, {x, y, z})] =
            0.5 * (streamwise_momentum(left) +
                   streamwise_momentum(right)) *
            x_area;
      }
    }
  }
  const double y_area = (kDomainX / static_cast<double>(kCells.x)) *
                        (kDomainZ / static_cast<double>(kCells.z));
  for (std::int32_t z = 0; z < y_faces.z; ++z) {
    for (std::int32_t y = 0; y < y_faces.y; ++y) {
      for (std::int32_t x = 0; x < y_faces.x; ++x) {
        const double global_x =
            (static_cast<double>(image.patch.begin.x + x) + 0.5) *
            kDomainX / static_cast<double>(kCells.x);
        const double phase = 2.0 * std::acos(-1.0) * global_x / kDomainX;
        const double absolute_pressure =
            kPressureReference + kPressureAmplitude * std::sin(phase);
        const double temperature =
            kBaseTemperature *
            std::pow(absolute_pressure / kPressureReference,
                     (gamma - 1.0) / gamma);
        const double local_density =
            absolute_pressure / (kGasConstant * temperature);
        image.final_mass_flux[1U][face_offset(y_faces, {x, y, z})] =
            local_density * kCrossflowAmplitude * std::sin(phase) * y_area;
      }
    }
  }
  return image;
}

struct Endpoint {
  double time{};
  double dt{};
  double pressure_reference{};
  std::uint64_t step{};
  Int3 cells{};
  std::vector<double> pressure_absolute;
  std::vector<double> enthalpy;
  std::vector<double> density;
  std::vector<double> temperature;
  std::vector<double> velocity;
  std::array<std::vector<double>, 3U> flux;
};

bool capture(ProductDriver& driver, Endpoint& out,
             RevisionToken expected_flux_revision = 0U) {
  RestartSnapshot snapshot;
  CommittedOutputSnapshot output;
  if (!driver.committed_restart_snapshot(snapshot) ||
      !driver.committed_output_snapshot(output) || !output.committed ||
      output.time != snapshot.time || output.step != snapshot.step ||
      !snapshot.final_mass_flux.certificate.valid() ||
      !snapshot.final_mass_flux.certificate.matches(snapshot.final_mass_flux) ||
      (expected_flux_revision != 0U &&
       snapshot.final_mass_flux.revision != expected_flux_revision))
    return false;

  const ConstFieldView* pressure = nullptr;
  const ConstFieldView* enthalpy = nullptr;
  const ConstFieldView* velocity = nullptr;
  bool pressure_is_absolute = false;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView& field = snapshot.fields.data[index];
    if (field.role == RestartFieldRole::pressure_perturbation)
      pressure = &field.values;
    else if (field.role == RestartFieldRole::pressure_absolute) {
      pressure = &field.values;
      pressure_is_absolute = true;
    } else if (field.role == RestartFieldRole::enthalpy)
      enthalpy = &field.values;
    else if (field.role == RestartFieldRole::velocity)
      velocity = &field.values;
  }
  if (pressure == nullptr || enthalpy == nullptr || velocity == nullptr)
    return false;

  Endpoint candidate;
  candidate.time = snapshot.time;
  candidate.dt = snapshot.dt;
  candidate.pressure_reference = snapshot.pressure_reference;
  candidate.step = snapshot.step;
  candidate.cells = snapshot.patch.cells;
  const std::size_t count = static_cast<std::size_t>(candidate.cells.x) *
                            candidate.cells.y * candidate.cells.z;
  candidate.pressure_absolute.reserve(count);
  candidate.enthalpy.reserve(count);
  candidate.density.reserve(count);
  candidate.temperature.reserve(count);
  candidate.velocity.reserve(3U * count);
  for (std::int32_t z = 0; z < candidate.cells.z; ++z) {
    for (std::int32_t y = 0; y < candidate.cells.y; ++y) {
      for (std::int32_t x = 0; x < candidate.cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double stored_pressure = pressure->unchecked(cell, 0U);
        const double p = pressure_is_absolute
                             ? stored_pressure
                             : snapshot.pressure_reference + stored_pressure;
        const double h = enthalpy->unchecked(cell, 0U);
        const double temperature_value = h / kHeatCapacity;
        const double rho = p / (kGasConstant * temperature_value);
        candidate.pressure_absolute.push_back(p);
        candidate.enthalpy.push_back(h);
        candidate.temperature.push_back(temperature_value);
        candidate.density.push_back(rho);
        for (std::uint8_t component = 0U; component < 3U; ++component)
          candidate.velocity.push_back(
              velocity->unchecked(cell, component));
      }
    }
  }
  const std::array<ConstFaceFieldView, 3U> faces{
      snapshot.final_mass_flux.x, snapshot.final_mass_flux.y,
      snapshot.final_mass_flux.z};
  for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
    const ConstFaceFieldView face = faces[axis];
    candidate.flux[axis].reserve(static_cast<std::size_t>(face.extents.x) *
                                 face.extents.y * face.extents.z);
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          candidate.flux[axis].push_back(face.unchecked({x, y, z}));
  }
  out = std::move(candidate);
  return true;
}

bool finite_positive(const Endpoint& endpoint) {
  const auto all_finite = [](const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
  };
  if (!all_finite(endpoint.pressure_absolute) ||
      !all_finite(endpoint.enthalpy) || !all_finite(endpoint.density) ||
      !all_finite(endpoint.temperature) || !all_finite(endpoint.velocity))
    return false;
  for (const std::vector<double>& face : endpoint.flux)
    if (!all_finite(face)) return false;
  for (std::size_t cell = 0U; cell < endpoint.density.size(); ++cell)
    if (!(endpoint.pressure_absolute[cell] > 0.0) ||
        !(endpoint.enthalpy[cell] > 0.0) ||
        !(endpoint.density[cell] > 0.0) ||
        endpoint.temperature[cell] < 200.0 ||
        endpoint.temperature[cell] > 2000.0)
      return false;
  return true;
}

double independent_continuity(const Endpoint& current,
                              const Endpoint& accepted,
                              const Endpoint& previous,
                              BdfCoefficients bdf) {
  if (current.cells.x != accepted.cells.x ||
      current.cells.y != accepted.cells.y ||
      current.cells.z != accepted.cells.z ||
      current.cells.x != previous.cells.x ||
      current.cells.y != previous.cells.y ||
      current.cells.z != previous.cells.z)
    return std::numeric_limits<double>::infinity();
  const Int3 cells = current.cells;
  const Int3 x_faces{cells.x + 1, cells.y, cells.z};
  const Int3 y_faces{cells.x, cells.y + 1, cells.z};
  const Int3 z_faces{cells.x, cells.y, cells.z + 1};
  const double volume = kDomainX * kDomainY * kDomainZ /
                        static_cast<double>(cells.x * cells.y * cells.z);
  double maximum = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::size_t offset = cell_offset(cells, cell);
        const double rho = current.density[offset];
        const double rho_n = accepted.density[offset];
        const double rho_nm1 = previous.density[offset];
        const double fxm = current.flux[0U][face_offset(x_faces, cell)];
        const double fxp =
            current.flux[0U][face_offset(x_faces, {x + 1, y, z})];
        const double fym = current.flux[1U][face_offset(y_faces, cell)];
        const double fyp =
            current.flux[1U][face_offset(y_faces, {x, y + 1, z})];
        const double fzm = current.flux[2U][face_offset(z_faces, cell)];
        const double fzp =
            current.flux[2U][face_offset(z_faces, {x, y, z + 1})];
        const double residual =
            volume * (bdf.a0 * rho + bdf.a1 * rho_n + bdf.a2 * rho_nm1) +
            (fxp - fxm) + (fyp - fym) + (fzp - fzm);
        const double scale =
            std::abs(volume * bdf.a0 * rho) +
            std::abs(volume * bdf.a1 * rho_n) +
            std::abs(volume * bdf.a2 * rho_nm1) + std::abs(fxm) +
            std::abs(fxp) + std::abs(fym) + std::abs(fyp) + std::abs(fzm) +
            std::abs(fzp);
        maximum = std::max(
            maximum,
            std::abs(residual) /
                std::max(scale, std::numeric_limits<double>::min()));
      }
    }
  }
  return maximum;
}

struct PhysicalGateSummary {
  double maximum_eos{};
  double maximum_continuity{};
  double maximum_energy{};
  double maximum_closed_mass{};
  double maximum_gauge{};
  double maximum_independent_continuity{};
  double minimum_pressure{std::numeric_limits<double>::infinity()};
  double minimum_enthalpy{std::numeric_limits<double>::infinity()};
  double minimum_density{std::numeric_limits<double>::infinity()};
  double minimum_temperature{std::numeric_limits<double>::infinity()};
  bool nonstationary_linear_solve{};
};

bool physical_step(const ValidatedModel& model, const DriverStepReport& report,
                   const Endpoint& current, const Endpoint& accepted,
                   const Endpoint& previous, PhysicalGateSummary& summary) {
  const bool report_valid =
      report.accepted && report.attempts == 1U &&
      report.failure.code == StatusCode::ok &&
      !report.temporal_method_fallback &&
      report.piso.pressure_solve_calls == 2U &&
      report.piso.pressure[0U].status && report.piso.pressure[1U].status &&
      std::isfinite(report.piso.eos_residual) &&
      report.piso.eos_residual <= model.solver.terminal.eos &&
      std::isfinite(report.piso.continuity_residual) &&
      report.piso.continuity_residual <= model.solver.terminal.continuity &&
      std::isfinite(report.piso.energy_residual) &&
      report.piso.energy_residual <= model.solver.terminal.continuity &&
      std::isfinite(report.piso.closed_mass_residual) &&
      report.piso.closed_mass_residual <= model.solver.terminal.closed_mass &&
      std::isfinite(report.piso.gauge_residual) &&
      report.piso.gauge_residual <= model.solver.terminal.gauge &&
      current.time == report.accepted_time &&
      current.step == report.accepted_step;
  summary.maximum_eos =
      std::max(summary.maximum_eos, report.piso.eos_residual);
  summary.maximum_continuity =
      std::max(summary.maximum_continuity, report.piso.continuity_residual);
  summary.maximum_energy =
      std::max(summary.maximum_energy, report.piso.energy_residual);
  summary.maximum_closed_mass = std::max(summary.maximum_closed_mass,
                                         report.piso.closed_mass_residual);
  summary.maximum_gauge =
      std::max(summary.maximum_gauge, report.piso.gauge_residual);
  summary.nonstationary_linear_solve =
      summary.nonstationary_linear_solve ||
      report.piso.pressure[0U].initial_true_residual > 0.0 ||
      report.piso.pressure[1U].initial_true_residual > 0.0;
  const double continuity = independent_continuity(
      current, accepted, previous, report.effective_bdf);
  summary.maximum_independent_continuity =
      std::max(summary.maximum_independent_continuity, continuity);
  for (const double value : current.pressure_absolute)
    summary.minimum_pressure = std::min(summary.minimum_pressure, value);
  for (const double value : current.enthalpy)
    summary.minimum_enthalpy = std::min(summary.minimum_enthalpy, value);
  for (const double value : current.density)
    summary.minimum_density = std::min(summary.minimum_density, value);
  for (const double value : current.temperature)
    summary.minimum_temperature = std::min(summary.minimum_temperature, value);
  if (!report_valid || !finite_positive(current) ||
      !std::isfinite(continuity) || continuity > 2.0e-9) {
    std::cerr << std::setprecision(17)
              << "step status=" << static_cast<unsigned>(report.failure.code)
              << '/' << report.failure.detail << " order="
              << static_cast<unsigned>(report.effective_bdf.order)
              << " accepted=" << report.accepted
              << " attempts=" << report.attempts
              << " fallback=" << report.temporal_method_fallback
              << " solve-status="
              << static_cast<unsigned>(report.piso.pressure[0U].status.code)
              << '/'
              << static_cast<unsigned>(report.piso.pressure[1U].status.code)
              << " time/step=" << current.time << '/' << current.step << ':'
              << report.accepted_time << '/' << report.accepted_step
              << " finite=" << finite_positive(current) << " solves="
              << static_cast<unsigned>(report.piso.pressure_solve_calls)
              << " residuals=" << report.piso.eos_residual << ','
              << report.piso.continuity_residual << ','
              << report.piso.energy_residual << ','
              << report.piso.closed_mass_residual << ','
              << report.piso.gauge_residual
              << " independent-continuity=" << continuity << '\n';
  }
  return report_valid && finite_positive(current) &&
         std::isfinite(continuity) && continuity <= 2.0e-9;
}

struct StepAudit {
  double time{};
  BdfCoefficients bdf{};
  StepOrigin origin{StepOrigin::fresh_start};
  std::uint32_t attempts{};
  bool temporal_method_fallback{};
  std::uint8_t pressure_solve_calls{};
  std::array<LinearTermination, 2U> linear_termination{};
  std::array<std::uint32_t, 2U> linear_iterations{};
  std::array<double, 2U> initial_true_residual{};
  std::array<double, 2U> final_true_residual{};
  bool thermophysical_limited{};
  double thermophysical_theta{1.0};
  std::uint32_t thermophysical_low_order_substeps{};
  bool momentum_limited{};
  double momentum_theta{1.0};
  std::uint32_t momentum_limiter_activations{};
  double eos_residual{};
  double continuity_residual{};
  double energy_residual{};
  double closed_mass_residual{};
  double gauge_residual{};
};

StepAudit audit_step(const DriverStepReport& report) {
  StepAudit audit;
  audit.time = report.accepted_time;
  audit.bdf = report.effective_bdf;
  audit.origin = report.proposal.origin;
  audit.attempts = report.attempts;
  audit.temporal_method_fallback = report.temporal_method_fallback;
  audit.pressure_solve_calls = report.piso.pressure_solve_calls;
  for (std::size_t corrector = 0U; corrector < 2U; ++corrector) {
    audit.linear_termination[corrector] =
        report.piso.pressure[corrector].termination;
    audit.linear_iterations[corrector] =
        report.piso.pressure[corrector].iterations;
    audit.initial_true_residual[corrector] =
        report.piso.pressure[corrector].initial_true_residual;
    audit.final_true_residual[corrector] =
        report.piso.pressure[corrector].final_true_residual;
  }
  audit.thermophysical_limited = report.thermophysical_predictor.limited;
  audit.thermophysical_theta = report.thermophysical_predictor.theta;
  audit.thermophysical_low_order_substeps =
      report.thermophysical_predictor.low_order_substeps;
  audit.momentum_limited = report.momentum_predictor_limiter.limited;
  audit.momentum_theta = report.momentum_predictor_limiter.theta;
  audit.momentum_limiter_activations =
      report.momentum_predictor_limiter.activations;
  audit.eos_residual = report.piso.eos_residual;
  audit.continuity_residual = report.piso.continuity_residual;
  audit.energy_residual = report.piso.energy_residual;
  audit.closed_mass_residual = report.piso.closed_mass_residual;
  audit.gauge_residual = report.piso.gauge_residual;
  return audit;
}

struct Trajectory {
  Endpoint common_start;
  Endpoint terminal;
  PhysicalGateSummary physical{};
  std::vector<StepAudit> common_steps;
  std::vector<StepAudit> production_steps;
};

bool run_level(const ValidatedModel& model, double dt, Trajectory& out) {
  CompiledCasePlan plan;
  Status status = ProductCompiler::compile(MPI_COMM_SELF, model, {}, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  RestartExpected expected;
  if (status) status = driver.restart_expected(expected);
  RestartImage image;
  if (status) image = compatible_wave(expected, dt);
  if (status) status = driver.initialize_restart(image);
  Endpoint initial;
  if (status && !capture(driver, initial))
    status = {StatusCode::invalid_plan, 0U};
  if (!status) {
    std::cerr << "level setup=" << static_cast<unsigned>(status.code) << '/'
              << status.detail << '\n';
    return false;
  }
  out.common_start = initial;

  const std::uint64_t step_count = static_cast<std::uint64_t>(
      std::llround(kEvolutionTime / dt));
  if (step_count < 8U ||
      std::abs(static_cast<double>(step_count) * dt - kEvolutionTime) >
          64.0 * std::numeric_limits<double>::epsilon() * kEvolutionTime)
    return false;

  Endpoint previous = initial;
  Endpoint accepted = initial;
  DriverStepReport startup;
  const LocalTimeLimits limits{dt, dt, dt, dt, dt};
  status = driver.advance(limits, startup);
  Endpoint current;
  if (status && !capture(driver, current, startup.piso.final_flux_revision))
    status = {StatusCode::invalid_plan, 0U};
  bool passed = static_cast<bool>(status) &&
                startup.proposal.origin == StepOrigin::restart &&
                backward_euler_coefficients(startup.effective_bdf, dt) &&
                startup.proposal.dt == dt &&
                physical_step(model, startup, current, accepted, previous,
                              out.physical);
  if (!passed) {
    std::cerr << "startup status=" << static_cast<unsigned>(status.code) << '/'
              << status.detail << " origin="
              << static_cast<unsigned>(startup.proposal.origin) << " dt="
              << std::setprecision(17) << startup.proposal.dt
              << " attempts/fallback/order=" << startup.attempts << '/'
              << startup.temporal_method_fallback << '/'
              << static_cast<unsigned>(startup.effective_bdf.order)
              << " failure/stage="
              << static_cast<unsigned>(startup.failure.code) << '/'
              << startup.failure.detail << '/' << startup.failed_stage
              << " residuals=" << startup.piso.eos_residual << ','
              << startup.piso.continuity_residual << ','
              << startup.piso.energy_residual << ','
              << startup.piso.closed_mass_residual << ','
              << startup.piso.gauge_residual << '\n';
    return false;
  }
  out.common_steps.push_back(audit_step(startup));
  previous = accepted;
  accepted = current;

  // Each refinement starts from the same exact Restart state, takes exactly
  // one BE recovery step at its own dt, then uses constant-step BDF2.  The BE
  // state has O(dt^2) error and therefore supplies a second-order-consistent
  // starting value without a large variable-ratio jump.
  for (std::uint64_t step_index = 1U; step_index < step_count; ++step_index) {
    DriverStepReport report;
    status = driver.advance(limits, report);
    current = {};
    if (status && !capture(driver, current, report.piso.final_flux_revision))
      status = {StatusCode::invalid_plan, 0U};
    passed = static_cast<bool>(status) && report.accepted &&
             report.proposal.origin == StepOrigin::accepted &&
             constant_bdf2_coefficients(report.effective_bdf, dt) &&
             report.proposal.dt == dt &&
             physical_step(model, report, current, accepted, previous,
                           out.physical);
    if (!passed) {
      std::cerr << "level dt=" << std::setprecision(17) << dt
                << " index=" << step_index << " status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << " accepted=" << report.accepted << " origin/order="
                << static_cast<unsigned>(report.proposal.origin) << '/'
                << static_cast<unsigned>(report.effective_bdf.order)
                << " proposed-dt=" << report.proposal.dt
                << " attempts=" << report.attempts << " prior="
                << static_cast<unsigned>(report.failure.code) << '/'
                << report.failure.detail << " prior-stage="
                << report.failed_stage << " residuals="
                << report.piso.eos_residual << ','
                << report.piso.continuity_residual << ','
                << report.piso.energy_residual << ','
                << report.piso.closed_mass_residual << ','
                << report.piso.gauge_residual << '\n';
      return false;
    }
    out.production_steps.push_back(audit_step(report));
    previous = accepted;
    accepted = std::move(current);
  }
  out.terminal = std::move(accepted);
  return out.production_steps.size() == step_count - 1U &&
         std::abs(out.terminal.time - kEvolutionTime) <=
             32.0 * std::numeric_limits<double>::epsilon() &&
         out.terminal.step == initial.step + step_count;
}

double rms_difference(const std::vector<double>& left,
                      const std::vector<double>& right) {
  if (left.size() != right.size() || left.empty())
    return std::numeric_limits<double>::infinity();
  long double sum = 0.0L;
  for (std::size_t index = 0U; index < left.size(); ++index) {
    const long double difference =
        static_cast<long double>(left[index]) - right[index];
    sum += difference * difference;
  }
  return std::sqrt(static_cast<double>(sum / left.size()));
}

double flux_rms_difference(const Endpoint& left, const Endpoint& right) {
  long double sum = 0.0L;
  std::size_t count = 0U;
  for (std::size_t axis = 0U; axis < left.flux.size(); ++axis) {
    if (left.flux[axis].size() != right.flux[axis].size())
      return std::numeric_limits<double>::infinity();
    for (std::size_t index = 0U; index < left.flux[axis].size(); ++index) {
      const long double difference =
          static_cast<long double>(left.flux[axis][index]) -
          right.flux[axis][index];
      sum += difference * difference;
      ++count;
    }
  }
  return count == 0U
             ? std::numeric_limits<double>::infinity()
             : std::sqrt(static_cast<double>(sum / count));
}

enum class FluxRegion : std::uint8_t { all, periodic_duplicate, interior };

double flux_axis_rms_difference(const Endpoint& left, const Endpoint& right,
                                std::size_t axis, FluxRegion region) {
  if (axis >= left.flux.size() || left.cells.x != right.cells.x ||
      left.cells.y != right.cells.y || left.cells.z != right.cells.z ||
      left.flux[axis].size() != right.flux[axis].size())
    return std::numeric_limits<double>::infinity();
  Int3 extents = left.cells;
  if (axis == 0U)
    ++extents.x;
  else if (axis == 1U)
    ++extents.y;
  else
    ++extents.z;
  long double sum = 0.0L;
  std::size_t count = 0U;
  for (std::int32_t z = 0; z < extents.z; ++z) {
    for (std::int32_t y = 0; y < extents.y; ++y) {
      for (std::int32_t x = 0; x < extents.x; ++x) {
        const std::int32_t normal = axis == 0U ? x : (axis == 1U ? y : z);
        const std::int32_t normal_extent =
            axis == 0U ? extents.x : (axis == 1U ? extents.y : extents.z);
        const bool duplicate = normal == 0 || normal + 1 == normal_extent;
        if ((region == FluxRegion::periodic_duplicate && !duplicate) ||
            (region == FluxRegion::interior && duplicate))
          continue;
        const std::size_t offset = face_offset(extents, {x, y, z});
        const long double difference =
            static_cast<long double>(left.flux[axis][offset]) -
            right.flux[axis][offset];
        sum += difference * difference;
        ++count;
      }
    }
  }
  return count == 0U
             ? std::numeric_limits<double>::infinity()
             : std::sqrt(static_cast<double>(sum / count));
}

void dump_flux_axis_diagnostics(const Endpoint& coarse, const Endpoint& fine,
                                const Endpoint& quarter) {
  constexpr std::array<std::string_view, 3U> names{{"x", "y", "z"}};
  for (std::size_t axis = 0U; axis < names.size(); ++axis) {
    const double all_coarse_fine = flux_axis_rms_difference(
        coarse, fine, axis, FluxRegion::all);
    const double all_fine_quarter = flux_axis_rms_difference(
        fine, quarter, axis, FluxRegion::all);
    const double duplicate_coarse_fine = flux_axis_rms_difference(
        coarse, fine, axis, FluxRegion::periodic_duplicate);
    const double duplicate_fine_quarter = flux_axis_rms_difference(
        fine, quarter, axis, FluxRegion::periodic_duplicate);
    const double interior_coarse_fine = flux_axis_rms_difference(
        coarse, fine, axis, FluxRegion::interior);
    const double interior_fine_quarter = flux_axis_rms_difference(
        fine, quarter, axis, FluxRegion::interior);
    const auto order = [](double coarse_fine, double fine_quarter) {
      return coarse_fine > 0.0 && fine_quarter > 0.0
                 ? std::log(coarse_fine / fine_quarter) / std::log(2.0)
                 : -std::numeric_limits<double>::infinity();
    };
    std::cerr << std::setprecision(12) << "phi-" << names[axis]
              << " all(E1,E2,q)=" << all_coarse_fine << ','
              << all_fine_quarter << ','
              << order(all_coarse_fine, all_fine_quarter)
              << " duplicate(E1,E2,q)=" << duplicate_coarse_fine << ','
              << duplicate_fine_quarter << ','
              << order(duplicate_coarse_fine, duplicate_fine_quarter)
              << " interior(E1,E2,q)=" << interior_coarse_fine << ','
              << interior_fine_quarter << ','
              << order(interior_coarse_fine, interior_fine_quarter) << '\n';
  }
}

using FaceVectors = std::array<std::vector<double>, 3U>;

FaceVectors reconstruct_velocity_flux(const Endpoint& endpoint) {
  FaceVectors reconstructed;
  const Int3 cells = endpoint.cells;
  const auto wrap = [](std::int32_t index, std::int32_t extent) {
    index %= extent;
    return index < 0 ? index + extent : index;
  };
  const std::array<double, 3U> areas{{
      (kDomainY / cells.y) * (kDomainZ / cells.z),
      (kDomainX / cells.x) * (kDomainZ / cells.z),
      (kDomainX / cells.x) * (kDomainY / cells.y),
  }};
  for (std::size_t axis = 0U; axis < reconstructed.size(); ++axis) {
    Int3 extents = cells;
    if (axis == 0U)
      ++extents.x;
    else if (axis == 1U)
      ++extents.y;
    else
      ++extents.z;
    reconstructed[axis].resize(static_cast<std::size_t>(extents.x) *
                               extents.y * extents.z);
    for (std::int32_t z = 0; z < extents.z; ++z) {
      for (std::int32_t y = 0; y < extents.y; ++y) {
        for (std::int32_t x = 0; x < extents.x; ++x) {
          Int3 right{wrap(x, cells.x), wrap(y, cells.y), wrap(z, cells.z)};
          Int3 left = right;
          if (axis == 0U)
            left.x = wrap(right.x - 1, cells.x);
          else if (axis == 1U)
            left.y = wrap(right.y - 1, cells.y);
          else
            left.z = wrap(right.z - 1, cells.z);
          const std::size_t left_cell = cell_offset(cells, left);
          const std::size_t right_cell = cell_offset(cells, right);
          const double left_momentum =
              endpoint.density[left_cell] *
              endpoint.velocity[3U * left_cell + axis];
          const double right_momentum =
              endpoint.density[right_cell] *
              endpoint.velocity[3U * right_cell + axis];
          reconstructed[axis][face_offset(extents, {x, y, z})] =
              0.5 * (left_momentum + right_momentum) * areas[axis];
        }
      }
    }
  }
  return reconstructed;
}

FaceVectors subtract_flux(const FaceVectors& left, const FaceVectors& right) {
  FaceVectors difference;
  for (std::size_t axis = 0U; axis < difference.size(); ++axis) {
    if (left[axis].size() != right[axis].size()) continue;
    difference[axis].resize(left[axis].size());
    for (std::size_t face = 0U; face < left[axis].size(); ++face)
      difference[axis][face] = left[axis][face] - right[axis][face];
  }
  return difference;
}

double vector_rms_difference(const std::vector<double>& left,
                             const std::vector<double>& right) {
  return rms_difference(left, right);
}

double face_vectors_rms_difference(const FaceVectors& left,
                                   const FaceVectors& right) {
  long double sum = 0.0L;
  std::size_t count = 0U;
  for (std::size_t axis = 0U; axis < left.size(); ++axis) {
    if (left[axis].size() != right[axis].size())
      return std::numeric_limits<double>::infinity();
    for (std::size_t face = 0U; face < left[axis].size(); ++face) {
      const long double difference =
          static_cast<long double>(left[axis][face]) - right[axis][face];
      sum += difference * difference;
      ++count;
    }
  }
  return count == 0U
             ? std::numeric_limits<double>::infinity()
             : std::sqrt(static_cast<double>(sum / count));
}

double vector_mean(const std::vector<double>& values) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  long double sum = 0.0L;
  for (const double value : values) sum += value;
  return static_cast<double>(sum / values.size());
}

double demeaned_rms_difference(const std::vector<double>& left,
                               const std::vector<double>& right) {
  if (left.size() != right.size() || left.empty())
    return std::numeric_limits<double>::infinity();
  const double left_mean = vector_mean(left);
  const double right_mean = vector_mean(right);
  long double sum = 0.0L;
  for (std::size_t face = 0U; face < left.size(); ++face) {
    const long double difference =
        (static_cast<long double>(left[face]) - left_mean) -
        (static_cast<long double>(right[face]) - right_mean);
    sum += difference * difference;
  }
  return std::sqrt(static_cast<double>(sum / left.size()));
}

struct DivergenceDiagnostic {
  double rms{};
  double maximum{};
  double maximum_normalized{};
};

DivergenceDiagnostic flux_divergence_diagnostic(const FaceVectors& flux,
                                                 Int3 cells) {
  const Int3 x_faces{cells.x + 1, cells.y, cells.z};
  const Int3 y_faces{cells.x, cells.y + 1, cells.z};
  const Int3 z_faces{cells.x, cells.y, cells.z + 1};
  long double sum = 0.0L;
  std::size_t count = 0U;
  DivergenceDiagnostic diagnostic;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const double fxm = flux[0U][face_offset(x_faces, {x, y, z})];
        const double fxp = flux[0U][face_offset(x_faces, {x + 1, y, z})];
        const double fym = flux[1U][face_offset(y_faces, {x, y, z})];
        const double fyp = flux[1U][face_offset(y_faces, {x, y + 1, z})];
        const double fzm = flux[2U][face_offset(z_faces, {x, y, z})];
        const double fzp = flux[2U][face_offset(z_faces, {x, y, z + 1})];
        const double divergence =
            (fxp - fxm) + (fyp - fym) + (fzp - fzm);
        const double scale = std::abs(fxp) + std::abs(fxm) + std::abs(fyp) +
                             std::abs(fym) + std::abs(fzp) + std::abs(fzm);
        sum += static_cast<long double>(divergence) * divergence;
        ++count;
        diagnostic.maximum =
            std::max(diagnostic.maximum, std::abs(divergence));
        diagnostic.maximum_normalized = std::max(
            diagnostic.maximum_normalized,
            std::abs(divergence) /
                std::max(scale, std::numeric_limits<double>::min()));
      }
    }
  }
  diagnostic.rms =
      count == 0U ? std::numeric_limits<double>::infinity()
                  : std::sqrt(static_cast<double>(sum / count));
  return diagnostic;
}

void dump_final_flux_decomposition(const Endpoint& coarse,
                                   const Endpoint& fine,
                                   const Endpoint& quarter) {
  const std::array<FaceVectors, 3U> advective{{
      reconstruct_velocity_flux(coarse), reconstruct_velocity_flux(fine),
      reconstruct_velocity_flux(quarter),
  }};
  const std::array<FaceVectors, 3U> final{{coarse.flux, fine.flux,
                                          quarter.flux}};
  std::array<FaceVectors, 3U> delta{};
  for (std::size_t level = 0U; level < delta.size(); ++level)
    delta[level] = subtract_flux(final[level], advective[level]);
  const auto order = [](double coarse_fine, double fine_quarter) {
    return coarse_fine > 0.0 && fine_quarter > 0.0
               ? std::log(coarse_fine / fine_quarter) / std::log(2.0)
               : -std::numeric_limits<double>::infinity();
  };
  const double advective_e1 =
      face_vectors_rms_difference(advective[0U], advective[1U]);
  const double advective_e2 =
      face_vectors_rms_difference(advective[1U], advective[2U]);
  std::cerr << std::setprecision(12) << "phi_U combined(E1,E2,q)="
            << advective_e1 << ',' << advective_e2 << ','
            << order(advective_e1, advective_e2) << '\n';
  constexpr std::array<std::string_view, 3U> names{{"x", "y", "z"}};
  for (std::size_t axis = 0U; axis < names.size(); ++axis) {
    const double phi_u_e1 = vector_rms_difference(
        advective[0U][axis], advective[1U][axis]);
    const double phi_u_e2 = vector_rms_difference(
        advective[1U][axis], advective[2U][axis]);
    const std::array<double, 3U> means{{
        vector_mean(delta[0U][axis]), vector_mean(delta[1U][axis]),
        vector_mean(delta[2U][axis]),
    }};
    const double mean_e1 = std::abs(means[0U] - means[1U]);
    const double mean_e2 = std::abs(means[1U] - means[2U]);
    const double demeaned_e1 = demeaned_rms_difference(
        delta[0U][axis], delta[1U][axis]);
    const double demeaned_e2 = demeaned_rms_difference(
        delta[1U][axis], delta[2U][axis]);
    std::cerr << "phi_U-" << names[axis] << "(E1,E2,q)=" << phi_u_e1
              << ',' << phi_u_e2 << ',' << order(phi_u_e1, phi_u_e2)
              << " delta-mean(levels/E1,E2,q)=" << means[0U] << ','
              << means[1U] << ',' << means[2U] << '/' << mean_e1 << ','
              << mean_e2 << ',' << order(mean_e1, mean_e2)
              << " delta-demeaned(E1,E2,q)=" << demeaned_e1 << ','
              << demeaned_e2 << ',' << order(demeaned_e1, demeaned_e2)
              << '\n';
  }
  constexpr std::array<std::string_view, 3U> levels{{"coarse", "fine",
                                                      "quarter"}};
  for (std::size_t level = 0U; level < delta.size(); ++level) {
    const DivergenceDiagnostic divergence =
        flux_divergence_diagnostic(delta[level], coarse.cells);
    std::cerr << "delta-div-" << levels[level] << "(rms,max,normalized-max)="
              << divergence.rms << ',' << divergence.maximum << ','
              << divergence.maximum_normalized << '\n';
  }
}

struct ConvergenceMetric {
  std::string_view name;
  double coarse_fine{};
  double fine_quarter{};
  double order{};
};

ConvergenceMetric metric(std::string_view name, double coarse_fine,
                         double fine_quarter) {
  return {name, coarse_fine, fine_quarter,
          coarse_fine > 0.0 && fine_quarter > 0.0
              ? std::log(coarse_fine / fine_quarter) / std::log(2.0)
              : -std::numeric_limits<double>::infinity()};
}

bool advance_ode(const std::vector<StepAudit>& steps, double decay_rate,
                 double& previous, double& accepted,
                 bool require_bdf2 = false) {
  for (const StepAudit& step : steps) {
    const BdfCoefficients bdf = step.bdf;
    if ((require_bdf2 && bdf.order != 2U) ||
        (!require_bdf2 && bdf.order != 1U && bdf.order != 2U) ||
        !std::isfinite(bdf.a0) || !std::isfinite(bdf.a1) ||
        !std::isfinite(bdf.a2) || !(bdf.a0 + decay_rate > 0.0))
      return false;
    const double target =
        -(bdf.a1 * accepted + bdf.a2 * previous) /
        (bdf.a0 + decay_rate);
    if (!std::isfinite(target)) return false;
    previous = accepted;
    accepted = target;
  }
  return true;
}

double mirror_history_ode_terminal(const Trajectory& trajectory,
                                   double decay_rate) {
  if (trajectory.common_steps.empty() ||
      trajectory.production_steps.empty() || !(decay_rate > 0.0))
    return std::numeric_limits<double>::infinity();
  // Exactly mirror the Restart path: y(0) is duplicated into n and n-1,
  // then the recorded level-dt BE recovery and production BDF2 coefficients
  // are replayed.  Any startup/history pollution in the Product fixture is
  // therefore also present in this scalar control.
  double previous = 1.0;
  double accepted = 1.0;
  if (!advance_ode(trajectory.common_steps, decay_rate, previous, accepted) ||
      !advance_ode(trajectory.production_steps, decay_rate, previous,
                   accepted, true))
    return std::numeric_limits<double>::infinity();
  return accepted;
}

double exact_branch_history_ode_terminal(const Trajectory& trajectory,
                                         double decay_rate) {
  if (trajectory.common_steps.empty() || trajectory.production_steps.empty() ||
      !(decay_rate > 0.0))
    return std::numeric_limits<double>::infinity();
  // Replace the branch operands by the analytic solution at the exact two
  // history times, then replay only Product's production BDF coefficients.
  // This removes all Restart/BE state error without changing the constant-dt
  // time-coefficient path under test.
  const double accepted_time = trajectory.common_steps.back().time;
  const double previous_time = trajectory.common_start.time;
  double previous = std::exp(-decay_rate * previous_time);
  double accepted = std::exp(-decay_rate * accepted_time);
  if (!advance_ode(trajectory.production_steps, decay_rate, previous,
                   accepted, true))
    return std::numeric_limits<double>::infinity();
  return accepted;
}

bool ode_self_convergence(std::string_view label,
                          const std::array<double, 3U>& terminals) {
  const double coarse_fine = std::abs(terminals[0U] - terminals[1U]);
  const double fine_quarter = std::abs(terminals[1U] - terminals[2U]);
  const double order =
      coarse_fine > 0.0 && fine_quarter > 0.0
          ? std::log(coarse_fine / fine_quarter) / std::log(2.0)
          : -std::numeric_limits<double>::infinity();
  const bool passed =
      std::all_of(terminals.begin(), terminals.end(), [](double value) {
        return std::isfinite(value);
      }) &&
      std::isfinite(order) && order >= 1.8;
  std::cerr << std::setprecision(12) << label << " ODE terminals="
            << terminals[0U] << ',' << terminals[1U] << ','
            << terminals[2U] << " differences=" << coarse_fine << ','
            << fine_quarter << " self-order=" << order << '\n';
  return passed;
}

bool bdf_history_ode_oracle(const std::array<Trajectory, 3U>& trajectories) {
  constexpr double decay_rate = 100.0;
  std::array<double, 3U> mirror{};
  std::array<double, 3U> analytic_branch{};
  for (std::size_t level = 0U; level < trajectories.size(); ++level) {
    mirror[level] =
        mirror_history_ode_terminal(trajectories[level], decay_rate);
    analytic_branch[level] =
        exact_branch_history_ode_terminal(trajectories[level], decay_rate);
  }
  return ode_self_convergence("mirror-startup", mirror) &&
         ode_self_convergence("analytic-branch-history", analytic_branch);
}

void dump_step_diagnostics(const std::array<Trajectory, 3U>& trajectories) {
  for (std::size_t level = 0U; level < trajectories.size(); ++level) {
    const double requested_dt = kCoarseDt / kRefinement[level];
    const auto dump = [&](std::string_view phase,
                          const std::vector<StepAudit>& steps) {
      for (std::size_t index = 0U; index < steps.size(); ++index) {
        const StepAudit& step = steps[index];
        std::cerr << std::setprecision(17) << "level=" << level
                  << " requested-dt=" << requested_dt << " phase=" << phase
                  << " step=" << index << " time=" << step.time
                  << " origin/order=" << static_cast<unsigned>(step.origin)
                  << '/' << static_cast<unsigned>(step.bdf.order) << " bdf="
                  << step.bdf.a0 << ',' << step.bdf.a1 << ',' << step.bdf.a2
                  << " attempts/fallback/solves=" << step.attempts << '/'
                  << step.temporal_method_fallback << '/'
                  << static_cast<unsigned>(step.pressure_solve_calls)
                  << " C1(term,iter,r0,rf)="
                  << static_cast<unsigned>(step.linear_termination[0U]) << ','
                  << step.linear_iterations[0U] << ','
                  << step.initial_true_residual[0U] << ','
                  << step.final_true_residual[0U]
                  << " C2(term,iter,r0,rf)="
                  << static_cast<unsigned>(step.linear_termination[1U]) << ','
                  << step.linear_iterations[1U] << ','
                  << step.initial_true_residual[1U] << ','
                  << step.final_true_residual[1U] << " terminal="
                  << step.eos_residual << ',' << step.continuity_residual
                  << ',' << step.energy_residual << ','
                  << step.closed_mass_residual << ',' << step.gauge_residual
                  << " limiter(thermo,theta,substeps/momentum,theta,count)="
                  << step.thermophysical_limited << ','
                  << step.thermophysical_theta << ','
                  << step.thermophysical_low_order_substeps << '/'
                  << step.momentum_limited << ',' << step.momentum_theta << ','
                  << step.momentum_limiter_activations
                  << '\n';
      }
    };
    dump("startup", trajectories[level].common_steps);
    dump("BDF2", trajectories[level].production_steps);
  }
}

void dump_gate_summaries(const std::array<Trajectory, 3U>& trajectories) {
  for (std::size_t level = 0U; level < trajectories.size(); ++level) {
    const Trajectory& trajectory = trajectories[level];
    std::array<std::uint32_t, 2U> maximum_iterations{};
    std::array<double, 2U> maximum_final_true_residual{};
    std::array<std::size_t, 2U> absolute_stops{};
    std::array<std::size_t, 2U> relative_stops{};
    std::size_t thermophysical_limited_steps = 0U;
    std::size_t momentum_limited_steps = 0U;
    double minimum_thermophysical_theta = 1.0;
    double minimum_momentum_theta = 1.0;
    const auto accumulate = [&](const std::vector<StepAudit>& steps) {
      for (const StepAudit& step : steps)
        for (std::size_t corrector = 0U; corrector < 2U; ++corrector) {
          maximum_iterations[corrector] =
              std::max(maximum_iterations[corrector],
                       step.linear_iterations[corrector]);
          maximum_final_true_residual[corrector] =
              std::max(maximum_final_true_residual[corrector],
                       step.final_true_residual[corrector]);
          absolute_stops[corrector] +=
              step.final_true_residual[corrector] <=
              kPressureLinearTolerance;
          relative_stops[corrector] +=
              step.final_true_residual[corrector] <=
              kPressureLinearTolerance *
                  step.initial_true_residual[corrector];
        }
      for (const StepAudit& step : steps) {
        thermophysical_limited_steps += step.thermophysical_limited;
        momentum_limited_steps += step.momentum_limited;
        minimum_thermophysical_theta =
            std::min(minimum_thermophysical_theta,
                     step.thermophysical_theta);
        minimum_momentum_theta =
            std::min(minimum_momentum_theta, step.momentum_theta);
      }
    };
    accumulate(trajectory.common_steps);
    accumulate(trajectory.production_steps);
    const PhysicalGateSummary& gate = trajectory.physical;
    std::cerr << std::setprecision(12) << "level=" << level
              << " dt/N=" << kCoarseDt / kRefinement[level] << '/'
              << trajectory.common_steps.size() +
                     trajectory.production_steps.size()
              << " BE/BDF2=" << trajectory.common_steps.size() << '/'
              << trajectory.production_steps.size()
              << " max-terminal(eos,C,E,mass,gauge)=" << gate.maximum_eos
              << ',' << gate.maximum_continuity << ',' << gate.maximum_energy
              << ',' << gate.maximum_closed_mass << ',' << gate.maximum_gauge
              << " max-independent-C="
              << gate.maximum_independent_continuity
              << " max-C1/C2(iter,rf)=" << maximum_iterations[0U] << ','
              << maximum_final_true_residual[0U] << '/'
              << maximum_iterations[1U] << ','
              << maximum_final_true_residual[1U]
              << " C1/C2-stop(abs,rel)=" << absolute_stops[0U] << ','
              << relative_stops[0U] << '/' << absolute_stops[1U] << ','
              << relative_stops[1U]
              << " limited(thermo,min-theta/momentum,min-theta)="
              << thermophysical_limited_steps << ','
              << minimum_thermophysical_theta << '/'
              << momentum_limited_steps << ',' << minimum_momentum_theta
              << " min(p,h,rho,T)=" << gate.minimum_pressure << ','
              << gate.minimum_enthalpy << ',' << gate.minimum_density << ','
              << gate.minimum_temperature << '\n';
  }
}

bool coherent_face_local_limiter(
    const std::array<Trajectory, 3U>& trajectories) {
  for (const Trajectory& trajectory : trajectories) {
    const auto consistent = [](const std::vector<StepAudit>& steps) {
      return std::all_of(steps.begin(), steps.end(), [](const StepAudit& step) {
        // The thermophysical IDP path must remain wholly inactive.  The
        // face-local momentum limiter may become inactive as dt is refined,
        // but it must never collapse the complete anti-diffusive correction.
        // The independent U/phi convergence gates below prove that this
        // smooth-extremum transition preserves the intended temporal order.
        return !step.thermophysical_limited &&
               step.thermophysical_theta == 1.0 &&
               step.thermophysical_low_order_substeps == 0U &&
               std::isfinite(step.momentum_theta) &&
               step.momentum_theta > 0.0 && step.momentum_theta <= 1.0 &&
               step.momentum_limited ==
                   (step.momentum_limiter_activations != 0U) &&
               (step.momentum_limited ? step.momentum_theta < 1.0
                                      : step.momentum_theta == 1.0);
      });
    };
    if (!consistent(trajectory.common_steps) ||
        !consistent(trajectory.production_steps))
      return false;
  }
  return true;
}

bool test_product_pressure_energy_temporal_convergence() {
  const ValidatedModel model = temporal_model();
  std::array<Trajectory, 3U> trajectories{};
  bool passed = true;
  for (std::size_t level = 0U; level < trajectories.size(); ++level)
    passed &= expect(run_level(model, kCoarseDt / kRefinement[level],
                               trajectories[level]),
                     "real ProductDriver trajectory advances");
  if (!passed) return false;
  dump_gate_summaries(trajectories);
  passed &= expect(
      coherent_face_local_limiter(trajectories),
      "face-local limiter reports remain coherent without global collapse");

  // The physical Restart fields are bitwise identical.  Only the Restart dt
  // metadata differs, so every level performs its own one-step BE recovery
  // followed by an entirely constant-step BDF2 trajectory.
  const auto exact_vector = [](const std::vector<double>& left,
                               const std::vector<double>& right) {
    return left == right;
  };
  bool common_start = true;
  for (std::size_t level = 1U; level < trajectories.size(); ++level) {
    const Endpoint& reference = trajectories[0U].common_start;
    const Endpoint& candidate = trajectories[level].common_start;
    common_start &= reference.time == candidate.time &&
                    reference.step == candidate.step &&
                    exact_vector(reference.pressure_absolute,
                                 candidate.pressure_absolute) &&
                    exact_vector(reference.enthalpy, candidate.enthalpy) &&
                    exact_vector(reference.density, candidate.density) &&
                    exact_vector(reference.temperature,
                                 candidate.temperature) &&
                    exact_vector(reference.velocity, candidate.velocity);
    for (std::size_t axis = 0U; axis < reference.flux.size(); ++axis)
      common_start &= exact_vector(reference.flux[axis],
                                   candidate.flux[axis]);
  }
  passed &= expect(common_start,
                   "physical Restart state is common to all levels");

  const Endpoint& coarse = trajectories[0U].terminal;
  const Endpoint& fine = trajectories[1U].terminal;
  const Endpoint& quarter = trajectories[2U].terminal;
  const std::array<ConvergenceMetric, 6U> metrics{{
      metric("p_abs", rms_difference(coarse.pressure_absolute,
                                      fine.pressure_absolute),
             rms_difference(fine.pressure_absolute,
                            quarter.pressure_absolute)),
      metric("h", rms_difference(coarse.enthalpy, fine.enthalpy),
             rms_difference(fine.enthalpy, quarter.enthalpy)),
      metric("rho", rms_difference(coarse.density, fine.density),
             rms_difference(fine.density, quarter.density)),
      metric("T", rms_difference(coarse.temperature, fine.temperature),
             rms_difference(fine.temperature, quarter.temperature)),
      metric("U", rms_difference(coarse.velocity, fine.velocity),
             rms_difference(fine.velocity, quarter.velocity)),
      metric("final_phi", flux_rms_difference(coarse, fine),
             flux_rms_difference(fine, quarter)),
  }};

  bool endpoint_converged = true;
  for (const ConvergenceMetric& value : metrics) {
    const bool converged = std::isfinite(value.coarse_fine) &&
                           std::isfinite(value.fine_quarter) &&
                           value.coarse_fine > 64.0 *
                                                   std::numeric_limits<double>::
                                                       epsilon() &&
                           value.fine_quarter > 64.0 *
                                                    std::numeric_limits<double>::
                                                        epsilon() &&
                           std::isfinite(value.order) && value.order >= 1.8;
    endpoint_converged &= converged;
    passed &= expect(converged, "endpoint field converges at order >= 1.8");
    std::cerr << std::setprecision(12) << value.name
              << " E(dt,dt/2)=" << value.coarse_fine
              << " E(dt/2,dt/4)=" << value.fine_quarter
              << " observed-order=" << value.order << '\n';
  }
  dump_flux_axis_diagnostics(coarse, fine, quarter);
  dump_final_flux_decomposition(coarse, fine, quarter);
  const bool ode_converged = bdf_history_ode_oracle(trajectories);
  passed &= expect(ode_converged,
                   "mirror-startup and analytic-history ODE controls with "
                   "ProductDriver BDF sequences are second order");
  if (!endpoint_converged || !ode_converged)
    dump_step_diagnostics(trajectories);

  const Endpoint& common = trajectories[0U].common_start;
  const double pressure_change =
      rms_difference(common.pressure_absolute, quarter.pressure_absolute);
  const double enthalpy_change =
      rms_difference(common.enthalpy, quarter.enthalpy);
  const double velocity_change =
      rms_difference(common.velocity, quarter.velocity);
  const double flux_change = flux_rms_difference(common, quarter);
  double maximum_energy = 0.0;
  bool nonstationary_linear = false;
  for (const Trajectory& trajectory : trajectories) {
    maximum_energy =
        std::max(maximum_energy, trajectory.physical.maximum_energy);
    nonstationary_linear =
        nonstationary_linear ||
        trajectory.physical.nonstationary_linear_solve;
  }
  const bool coupled_witness =
      pressure_change > 1.0e-5 && enthalpy_change > 1.0e-5 &&
      velocity_change > 1.0e-8 && flux_change > 1.0e-12 &&
      maximum_energy > std::numeric_limits<double>::denorm_min() &&
      nonstationary_linear;
  passed &= expect(
      coupled_witness,
      "nonstationary p/h/flux response and nonzero energy audit exclude a "
      "pressure-only or stationary bypass");
  std::cerr << "coupling witness=" << pressure_change << ','
            << enthalpy_change << ',' << velocity_change << ',' << flux_change
            << " max-energy=" << maximum_energy
            << " nonstationary-linear=" << nonstationary_linear << '\n';
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed = test_product_pressure_energy_temporal_convergence();
  MPI_Finalize();
  return passed ? 0 : 1;
}
