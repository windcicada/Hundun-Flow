// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

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
  RestartStorageCompatibility restart_storage_compatibility{
      RestartStorageCompatibility::strict};
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
struct PressureEnergyCandidateWorkReport {
  std::uint32_t baseline_evaluations{};
  std::uint32_t extrapolation_evaluations{};
  std::uint32_t ladder_evaluations{};
  std::uint32_t incomplete_evaluations{};
  std::uint32_t rejected_extrapolations{};
  // Local inclusive evaluation time, not a rank maximum or complete step time.
  std::uint64_t local_evaluation_nanoseconds{};
  // Disjoint inclusive subphases: state copy/revision, velocity/correction
  // halo, thermophysics, state halo/boundary/derived transport, flux,
  // certificate, residual assembly/audit, final equivalence/hash. Not pure
  // kernel costs.
  std::array<std::uint64_t, 8U> local_phase_nanoseconds{};
};

struct PressureEnergyExtrapolationReport {
  bool attempted{};
  bool complete{};
  bool selection_attempted{};
  bool selected{};
  double alpha{};
  Status evaluation_status{};
  Status selection_status{};
  PressureEnergyGlobalizationSample sample{};
};

struct PressureEnergyGlobalizationIterationReport {
  bool valid{};
  std::uint8_t corrector{};
  std::uint8_t refinement_iteration{};
  double maximum_absolute_pressure_correction{};
  double maximum_absolute_enthalpy_correction{};
  PressureEnergyGlobalizationSample baseline{};
  PressureEnergyGlobalizationSample selected{};
  bool jacobian_scope_valid{};
  PressureEnergyJacobianScope jacobian_scope{
      PressureEnergyJacobianScope::generic_algebraic_quasi_newton};
  PressureEnergyCandidateWorkReport work{};
  PressureEnergyExtrapolationReport extrapolation{};
};

inline constexpr std::size_t kPressureEnergyGlobalizationTrajectoryCapacity =
    2U + kPressureEnergyRefinementCapacity;

struct PressureEnergyGlobalizationAttemptReport {
  bool valid{};
  std::uint8_t corrector{};
  // Occupied legacy candidate slots in the latest nonstationary loop:
  // an accepted extrapolation, or the final backtracking ladder. Not total
  // work.
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
  // Includes unsuccessful evaluations and all loops in this numerical attempt.
  // A subsequent time-step retry starts a new attempt report.
  PressureEnergyCandidateWorkReport work{};
  PressureEnergyCandidateWorkReport last_loop_work{};
  PressureEnergyExtrapolationReport last_extrapolation{};
  // Disjoint local preparation, Krylov solve, and close/enthalpy recovery.
  // Includes all C1/C2/refinement solves in this numerical attempt.
  std::array<std::uint64_t, 3U> local_solve_nanoseconds{};
  std::array<Status, kPressureEnergyGlobalizationCandidateCount>
      candidate_evaluation_status{};
};

enum class ApplicationFailurePhase : std::uint8_t {
  none,
  advance,
  visit,
  screen,
  monitor,
  restart,
  resources,
  evidence,
  input,
  case_compile,
  product_compile,
  driver_create,
  initialize,
  restart_load,
  runtime_identity,
  time_control
};

struct StepAttemptFailure {
  Status failure{};
  StageId stage{};
  std::uint32_t attempt{};  // One-based; zero means no failed attempt occurred.
  double dt{};
};

struct StepCompletionReport {
  Status outcome{};  // Exactly the status returned by advance.
  // The controller/commit decision, before preserving a legacy attempt error
  // in the return value. For example minimum_dt (454) or retry capacity (455).
  // Both outcome and stop_reason are OK after an accepted retry.
  Status stop_reason{};
  StepAttemptFailure first_failure{};
  StepAttemptFailure last_failure{};
};

struct ApplicationRunReport {
  PlanFingerprint case_model{};
  PlanFingerprint product{};
  std::uint64_t accepted_steps{};
  double final_time{};
  StageId failed_stage{};
  std::uint32_t attempts{};
  Status failure{};
  ApplicationFailurePhase failure_phase{ApplicationFailurePhase::none};
  StepCompletionReport step_completion{};
  // First proposal in this advance; subsequent retry dt remains in the
  // numerical attempt diagnostics. Present even for a minimum-dt rejection.
  TimeProposalDiagnostic initial_time_proposal{};
  BdfCoefficients requested_bdf{};
  BdfCoefficients effective_bdf{};
  std::uint8_t thermophysical_predictor_calls{};
  bool temporal_method_fallback{};
  PisoAttemptReport piso{};
  MomentumPredictorLimiterReport momentum_predictor_limiter{};
  MomentumPredictorSolveReport momentum_predictor_solve{};
  NumericalFailureContext numerical_failure{};
  ThermophysicalPredictorDiagnostics thermophysical_predictor{};
  PressureEnergyGlobalizationAttemptReport pressure_energy_globalization{};
  std::uint64_t predictor_limiter_activations{};
  std::uint64_t predictor_low_order_transport_passes{};
  std::uint64_t predictor_low_order_halo_exchanges{};
  double minimum_predictor_theta{1.0};
  double maximum_advective_convective_cfl_out{};
  double maximum_advective_convective_cfl_abs{};
  double advective_convective_cfl_limit{};
  double maximum_committed_convective_cfl_out{};
  double maximum_committed_convective_cfl_abs{};
  double committed_convective_cfl_limit{};
  bool restart_storage_migrated{};
  PlanFingerprint restart_source_plan{};
  PlanFingerprint restart_source_schema{};
  // Local complete ApplicationService::run interval, including setup, output,
  // evidence, and cleanup. Not the legacy advance-only evidence step timer.
  std::uint64_t local_run_nanoseconds{};
  // Disjoint totals: setup, time control, advance, bookkeeping, Visit,
  // screen/monitor, restart, resource sampling, evidence encoding/writing.
  std::array<std::uint64_t, 9U> local_phase_nanoseconds{};
  IoFailureContext io_failure{};
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
  std::uint64_t structured_control_collectives{};
  std::uint64_t ibm_control_collectives{};
};

inline constexpr std::size_t kDriverTimedStageCapacity = 10U;

struct DriverStageTiming {
  StageId stage{};
  std::uint64_t nanoseconds{};
};

struct DriverStepReport {
  StepCompletionReport completion{};
  TimeProposalDiagnostic initial_time_proposal{};
  StepTime proposal{};
  BdfCoefficients effective_bdf{};
  PisoAttemptReport piso{};
  std::uint64_t accepted_step{};
  double accepted_time{};
  std::uint32_t attempts{};
  StageId failed_stage{};
  Status failure{};  // Legacy last-attempt error; use completion.outcome for
                     // the result.
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
  Status restart_expected(RestartExpected& out,
                          RestartStorageCompatibility compatibility =
                              RestartStorageCompatibility::strict) noexcept;
  Status initialize(const DriverInitialState& initial) noexcept;
  Status initialize_restart(const RestartImage& image,
                            RestartStorageCompatibility compatibility =
                                RestartStorageCompatibility::strict) noexcept;
  // Add the accepted-flow convective bound without changing the caller's
  // explicit diffusion/acoustic bounds. No-op for a fixed-dt product.
  Status constrain_convective_time_limit(LocalTimeLimits& limits) noexcept;
  // Programmatic callers retain authority over explicit physical time scales.
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
