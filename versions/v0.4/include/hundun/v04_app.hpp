// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_product.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <cstdint>

namespace hundun::v04 {

struct CaseValidationReport {
  PlanFingerprint case_model{};
  PlanFingerprint product{};
  PlanSummary summary{};
};

struct ApplicationRunOptions {
  std::filesystem::path case_root;
  std::filesystem::path run_directory;
  std::filesystem::path source_root;
  std::filesystem::path restart_directory;
  std::uint64_t steps{};
  std::uint64_t output_interval{1U};
  std::uint64_t restart_interval{1U};
  LocalTimeLimits time_limits{1.0, 1.0, 1.0, 1.0, 1.0};
};

inline constexpr std::size_t kNumericalFailureMassFractionCapacity = 64U;

enum class NumericalFailureField : std::uint8_t {
  none,
  density,
  enthalpy,
  temperature,
  species,
  passive_scalar,
  pressure
};

enum class NumericalFailureContributor : std::uint8_t {
  none,
  density_predictor,
  conservative_history,
  accepted_advection,
  accepted_nonadvective,
  previous_advection,
  previous_nonadvective,
  thermodynamic_inversion,
  derived_state_lifecycle
};

struct NumericalFailureContext {
  bool valid{};
  NumericalFailureField field{NumericalFailureField::none};
  NumericalFailureContributor first_bad_contributor{
      NumericalFailureContributor::none};
  FieldId field_id{};
  Status failure{};
  StageId stage{};
  std::uint64_t global_cell{};
  Int3 global_index{};
  std::int32_t rank{-1};
  std::uint64_t attempted_step{};
  std::uint64_t generation{};
  double time{};
  double target_time{};
  double dt_before{};
  BdfCoefficients method_before{};
  StepOrigin origin_before{StepOrigin::fresh_start};
  bool retry_proposed{};
  double dt_after{};
  BdfCoefficients method_after{};
  StepOrigin origin_after{StepOrigin::fresh_start};

  double failed_value{};
  double allowed_minimum{};
  double allowed_maximum{};
  double pressure_absolute{};
  double temperature_before{};
  double temperature_estimate{};
  double density_before{};
  double density_previous{};
  double density_predicted{};
  double enthalpy_before{};
  double enthalpy_previous{};
  double cp_before{};
  std::array<double, kNumericalFailureMassFractionCapacity> mass_fractions{};
  std::size_t mass_fraction_count{};
  bool mass_fractions_truncated{};
  RevisionToken enthalpy_accepted_revision{};
  RevisionToken enthalpy_previous_revision{};
  RevisionToken nonadvective_accepted_revision{};
  RevisionToken nonadvective_previous_revision{};
  RevisionToken mass_flux_accepted_revision{};
  RevisionToken mass_flux_previous_revision{};
  RevisionToken temperature_accepted_revision{};
  RevisionToken conductivity_revision{};
  RevisionToken velocity_gradient_revision{};
  RevisionToken effective_viscosity_revision{};
  bool physical_boundary_stencil{};
  bool mpi_boundary_stencil{};

  double mass_divergence_accepted{};
  double mass_divergence_previous{};
  double advection_accepted{};
  double advection_previous{};
  double nonadvective_accepted{};
  double nonadvective_previous{};
  double conservative_history_value{};
  double accepted_advection_delta{};
  double accepted_nonadvective_delta{};
  double previous_advection_delta{};
  double previous_nonadvective_delta{};
  double reconstructed_value{};

  double diffusion_accepted{};
  double pressure_work_accepted{};
  double viscous_dissipation_accepted{};
  double explicit_source_accepted{};
  double implicit_sink_accepted{};
  bool rate_breakdown_complete{};
  bool immersed_interface_cell{};

  bool face_envelope_checked{};
  bool face_envelope_valid{};
  double maximum_face_envelope_violation{};
  double selected_face_value{};
  double selected_donor_minimum{};
  double selected_donor_maximum{};
  double counterfactual_be_density{};
  double counterfactual_be_enthalpy{};
  bool counterfactual_be_admissible{};
};

// Failure-visible record of the exact frozen-momentum continuity--energy
// globalization.  Samples already contain global residual norms; copying
// them into the step report adds no diagnostic collective or mutable alias.
struct PressureEnergyGlobalizationIterationReport {
  bool valid{};
  std::uint8_t corrector{};
  std::uint8_t refinement_iteration{};
  double maximum_absolute_pressure_correction{};
  double maximum_absolute_enthalpy_correction{};
  PressureEnergyGlobalizationSample baseline{};
  PressureEnergyGlobalizationSample selected{};
};

inline constexpr std::size_t kPressureEnergyGlobalizationTrajectoryCapacity =
    2U + kPressureEnergyRefinementCapacity;

struct PressureEnergyGlobalizationAttemptReport {
  bool valid{};
  std::uint8_t corrector{};
  std::uint8_t sample_count{};
  std::uint8_t trajectory_count{};
  double maximum_absolute_pressure_correction{};
  double maximum_absolute_enthalpy_correction{};
  PressureEnergyGlobalizationSample baseline{};
  std::array<PressureEnergyGlobalizationSample,
             kPressureEnergyGlobalizationCandidateCount>
      candidates{};
  std::array<PressureEnergyGlobalizationIterationReport,
             kPressureEnergyGlobalizationTrajectoryCapacity>
      trajectory{};
};

struct ApplicationRunReport {
  PlanFingerprint case_model{};
  PlanFingerprint product{};
  std::uint64_t accepted_steps{};
  double final_time{};
  StageId failed_stage{};
  std::uint32_t attempts{};
  Status failure{};
  BdfCoefficients requested_bdf{};
  BdfCoefficients effective_bdf{};
  std::uint8_t thermophysical_predictor_calls{};
  bool temporal_method_fallback{};
  PisoAttemptReport piso{};
  MomentumPredictorSolveReport momentum_predictor_solve{};
  NumericalFailureContext numerical_failure{};
  ThermophysicalPredictorDiagnostics thermophysical_predictor{};
  PressureEnergyGlobalizationAttemptReport pressure_energy_globalization{};
  std::uint64_t predictor_limiter_activations{};
  std::uint64_t predictor_low_order_transport_passes{};
  std::uint64_t predictor_low_order_halo_exchanges{};
  double minimum_predictor_theta{1.0};
};

class ApplicationService {
 public:
  static Status validate(MPI_Comm communicator,
                         const std::filesystem::path& case_root,
                         CaseValidationReport& out) noexcept;
  static Status initialize_case_directory(
      const std::filesystem::path& output_directory) noexcept;
  static Status validate_run_directories(
      const std::filesystem::path& case_root,
      const std::filesystem::path& run_directory,
      const std::filesystem::path& source_root) noexcept;
  static Status run(MPI_Comm communicator,
                    const ApplicationRunOptions& options,
                    ApplicationRunReport& report) noexcept;
};

struct DriverInitialState {
  double pressure_reference{101325.0};
  double temperature{300.0};
  Real3 velocity{};
  Span<const double> transported_scalars{};
  double start_time{};
};

struct DriverResourceReport {
  std::uint64_t structured_exchanges{};
  std::uint64_t structured_messages{};
  std::uint64_t structured_bytes{};
  std::uint64_t ibm_exchanges{};
  std::uint64_t ibm_messages{};
  std::uint64_t ibm_bytes{};
  std::uint64_t reduction_collectives{};
  std::uint64_t predictor_blocking_collectives{};
  std::uint64_t reduction_nanoseconds{};
  std::uint64_t reduction_logical_bytes{};
  std::uint64_t reduction_tree_messages{};
  std::uint64_t mg_blocking_collectives{};
  std::uint64_t mg_collective_logical_bytes{};
  std::uint64_t linear_iterations{};
  std::uint64_t exact_numeric_refills{};
  std::uint64_t hierarchy_rebuilds{};
  std::uint64_t preconditioner_applications{};
};

inline constexpr std::size_t kDriverTimedStageCapacity = 10U;

struct DriverStageTiming {
  StageId stage{};
  std::uint64_t nanoseconds{};
};

struct DriverStepReport {
  StepTime proposal{};
  BdfCoefficients effective_bdf{};
  PisoAttemptReport piso{};
  std::uint64_t accepted_step{};
  double accepted_time{};
  std::uint32_t attempts{};
  StageId failed_stage{};
  Status failure{};
  NumericalFailureContext numerical_failure{};
  ThermophysicalPredictorDiagnostics thermophysical_predictor{};
  PressureEnergyGlobalizationAttemptReport pressure_energy_globalization{};
  std::uint8_t thermophysical_predictor_calls{};
  bool temporal_method_fallback{};
  MomentumPredictorLimiterReport momentum_predictor_limiter{};
  MomentumPredictorSolveReport momentum_predictor_solve{};
  DriverResourceReport resources{};
  std::array<DriverStageTiming, kDriverTimedStageCapacity> stage_timings{};
  std::size_t stage_timing_count{};
  bool accepted{};
};

class ProductDriver {
 public:
  ProductDriver() noexcept = default;
  ~ProductDriver() noexcept;
  ProductDriver(const ProductDriver&) = delete;
  ProductDriver& operator=(const ProductDriver&) = delete;
  ProductDriver(ProductDriver&&) noexcept;
  ProductDriver& operator=(ProductDriver&&) noexcept;

  static Status create(MPI_Comm communicator, CompiledCasePlan&& plan,
                       ProductDriver& out) noexcept;
  Status restart_expected(RestartExpected& out) noexcept;
  Status initialize(const DriverInitialState& initial) noexcept;
  Status initialize_restart(const RestartImage& image) noexcept;
  Status advance(LocalTimeLimits limits, DriverStepReport& report) noexcept;
  Status committed_output_snapshot(CommittedOutputSnapshot& out) noexcept;
  Status committed_restart_snapshot(RestartSnapshot& out) noexcept;
  Status committed_surface_force(SurfaceForce& force,
                                 FinalForceCertificate& certificate) const
      noexcept;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  // Read-only test seam for the committed Fresh final-flux transaction.
  // Production Restart remains unavailable until at least one step commits.
  Status committed_final_mass_flux_for_test(ConstFaceFluxView& out) const
      noexcept;
#endif
  bool initialized() const noexcept;
  double pressure_reference() const noexcept;
  double closed_mass_target() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

}  // namespace hundun::v04
