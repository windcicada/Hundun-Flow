// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_product.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <cstdint>
#include <optional>
#include <iosfwd>

namespace hundun::v04 {

struct CaseValidationReport {
  PlanFingerprint case_model{};
  PlanFingerprint product{};
  PlanSummary summary{};
};

struct DriverInitialState {
  double pressure_reference{101325.0};
  double temperature{300.0};
  Real3 velocity{};
  // Species are mass fractions in [0,1] (including an admissible balance
  // species); passive scalars are finite signed quantities, not fractions.
  Span<const double> transported_scalars{};
  double start_time{};
};

struct ApplicationRunOptions {
  // Paths must refer to the same logical shared inputs/run/checkpoint across
  // ranks; local mount spellings are not compared as raw bytes. At the cold
  // run entry, steps, intervals, restart/initial-state presence and recovery
  // policies must agree. time_limits are rank-local physical candidates.
  std::filesystem::path case_root;
  std::filesystem::path run_directory;
  std::filesystem::path source_root;
  std::filesystem::path restart_directory;
  std::uint64_t steps{};
  std::uint64_t output_interval{1U};
  std::uint64_t restart_interval{1U};
  // Caller-supplied operator time scales. The application tightens convection
  // from the accepted state; remaining entries are NOT automatic physical
  // diffusion/acoustic estimates. Fixed-dt control does not use these entries.
  LocalTimeLimits time_limits{1.0, 1.0, 1.0, 1.0, 1.0};
  RestartStorageCompatibility restart_storage_compatibility{
      RestartStorageCompatibility::strict};
  // Explicit consumer policy. Unknown history remains rejected by default;
  // rebuilding never relaxes source integrity or physical compatibility.
  RestartHistoryPolicy restart_history_policy{
      RestartHistoryPolicy::require_compatible};
  // Synchronous borrowed scalar values, in the frozen catalog order. Fresh
  // application evidence starts at t=0. Mutually exclusive with a restart.
  // If absent, consistent boundary hints/defaults are used, never last-wins.
  std::optional<DriverInitialState> initial_state{};
  // Compile-generated identity of the final CLI target. Empty means the
  // library core manifest for callers embedding ApplicationService.
  std::string_view target_build_manifest{};
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

enum class NumericalCellRegion : std::uint8_t {
  unknown, fluid, solid_placeholder
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
  NumericalCellRegion runtime_region{NumericalCellRegion::unknown};
  ThermoInversionDiagnostic inversion{};
};

// The caller owns the existing MPI failure boundary and prints on root only.
// A failed log is secondary: this function never changes the numerical record.
Status write_numerical_failure(std::ostream& stream,
                               const NumericalFailureContext& failure) noexcept;

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

enum class PressureEnergySolveKind : std::uint8_t {
  pressure_continuity, diagonal_schur, spatial_schur
};

struct PressureEnergySolveObservation {
  std::uint8_t corrector{};
  std::uint8_t refinement{};
  PressureEnergySolveKind kind{PressureEnergySolveKind::pressure_continuity};
  bool invoked{};
  std::array<std::uint64_t, 3U> local_nanoseconds{};
  std::uint64_t mg_refill_nanoseconds{}, mg_copy_nanoseconds{};
  std::uint64_t structured_wait_nanoseconds{}, structured_control_nanoseconds{};
};

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
  std::uint8_t solve_observation_count{};
  std::array<PressureEnergySolveObservation,
             kPressureEnergyGlobalizationTrajectoryCapacity> solve_observations{};
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

Status write_step_completion_failure(std::ostream& stream,
                                     const StepCompletionReport& completion) noexcept;

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
  RestartWriteReport restart_output{};
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
  std::uint64_t structured_wait_nanoseconds{}, structured_control_nanoseconds{};
};

inline constexpr std::size_t kDriverTimedStageCapacity = 10U;

struct DriverStageTiming {
  StageId stage{};
  std::uint64_t nanoseconds{};
};

struct PressureEnergyPerformanceTotals {
  // Sums across ALL numerical attempts in this advance, including rejected
  // attempts. Last-attempt trajectories above remain a separate diagnostic.
  PressureEnergyCandidateWorkReport candidate{};
  std::array<std::uint64_t, 3U> solve_nanoseconds{};
  std::array<std::uint64_t, 3U> calls_by_kind{};
  std::array<std::array<std::uint64_t, 3U>, 3U> nanoseconds_by_kind{};
  // A apply, M apply, main Arnoldi dot, reduction, update (local partial costs).
  std::array<std::uint64_t, 5U> krylov_nanoseconds{};
  // Fixed diagnostic storage, never grows RSS. Overflow is explicit and makes
  // the detailed observation unusable for attribution, not a solver failure.
  struct Loop {
    std::uint32_t attempt{};
    std::uint32_t scalar_coupling_sweep{1U};
    double dt{};
    Status attempt_status{};
    PressureEnergySolveObservation solve{};
    LinearSolveResult linear{};
    PressureEnergyGlobalizationIterationReport globalization{};
  };
  std::array<Loop, 64U> loops{};
  std::uint32_t loop_count{}, dropped_loops{};
  // Disjoint final momentum assembly, terminal metrics, boundary ledger costs.
  std::array<std::uint64_t, 3U> final_audit_nanoseconds{};
};

// Same-state, same-final-flux observations, not additional acceptance
// tolerances. Momentum values are cell-integrated forces (N), energy values
// are powers (W). The total-equation defect is R_E + U.R_m - |U|^2 R_C/2;
// it is not the physical boundary total-energy balance (discrete product-rule
// defects must be measured separately).
struct DriverTerminalEquationReport {
  bool valid{};
  RevisionToken final_flux{};
  std::array<double, 3U> momentum_linf{};
  std::array<double, 3U> momentum_l1{};
  double continuity_signed{};
  double energy_signed{};
  double energy_l1{};
  double total_equation_defect{};
  double mass{};
  double internal_energy{};
  double kinetic_energy{};
  // Normalization: |R_m|/(a0*rho*V*U_rms), U_rms=sqrt(2*K/M).
  // A temporal inertial-force reference, NOT a convergence gate/backward error.
  bool momentum_normalization_valid{};
  double momentum_reference_velocity{};
  std::array<ReductionMaximumLocation, 3U> momentum_worst{};
  std::array<double, 3U> momentum_normalized_linf{};
  // 0: active rows touching an IBM control link; 1: all other active rows.
  std::array<std::uint64_t, 2U> momentum_region_cells{};
  std::array<std::array<double, 3U>, 2U> momentum_region_normalized_linf{};
  std::array<std::array<double, 3U>, 2U> momentum_region_normalized_rms{};
};

// Physical external-boundary balance for the current fixed, source-free
// product: stationary adiabatic IBM walls perform no physical heat/work.
// All rates are positive outwards except heat/stress input (positive inwards).
// The cumulative defects compare inventory with a BDF-integrated boundary
// ledger. A restart starts a new explicitly labelled observation epoch.
struct DriverConservationReport {
  bool valid{};
  std::uint64_t epoch_start_step{};
  double mass_outflow{};                 // kg/s
  double enthalpy_outflow{};             // W
  double kinetic_energy_outflow{};       // W
  double conductive_heat_input{};        // W, external physical faces
  double viscous_work_input{};           // W, external physical faces
  double mass_bdf_rate{};                // kg/s
  double total_energy_bdf_rate{};        // W
  double mass_balance_defect{};          // kg/s
  double total_energy_balance_defect{};  // W
  double cumulative_mass_defect{};       // kg since epoch_start_step
  double cumulative_energy_defect{};     // J since epoch_start_step
};

struct DriverScalarTransportReport {
  bool active{};
  std::uint64_t owned_payload_bytes{}; // Separate from halo/MPI and RSS.
  std::uint32_t coupling_sweeps{};
  std::uint32_t remap_iterations{};
  std::uint64_t remap_nanoseconds{};
  double final_species_residual{};
  double final_remap_residual{};
  double mass_pairing_residual{};
};

struct DriverCellTraceWindow {
  std::array<Int3, 2U> cells{};
  std::size_t count{}; // zero disables tracing
  std::uint64_t first_step{};
  std::uint64_t last_step{};
};

struct DriverCellTraceSample {
  Int3 global_index{};
  std::uint64_t step{};
  std::uint32_t attempt{}, composition_sweep{};
  StageId stage{};
  std::uint8_t active{};
  double rho{}, h{}, temperature{}, pressure{}, rate{};
  double rho_accepted{}, rho_previous{}, h_accepted{}, h_previous{};
};

struct DriverCellTrace {
  std::array<DriverCellTraceSample, 96U> samples{};
  std::size_t count{};
  std::size_t dropped{};
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
  PressureEnergyPerformanceTotals pressure_energy_performance{};
  DriverTerminalEquationReport terminal_equations{};
  DriverConservationReport conservation{};
  DriverScalarTransportReport scalar_transport{};
  DriverCellTrace cell_trace{};
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
                              RestartStorageCompatibility::strict,
                          RestartHistoryPolicy history_policy =
                              RestartHistoryPolicy::require_compatible) noexcept;
  Status initialize(const DriverInitialState& initial) noexcept;
  // Collective cold configuration. Same global cells/window on every rank;
  // samples remain rank-local until the caller's ordinary reporting boundary.
  Status set_cell_trace_window(const DriverCellTraceWindow& window) noexcept;
  Status initialize_restart(const RestartImage& image,
                            RestartStorageCompatibility compatibility =
                                RestartStorageCompatibility::strict,
                            RestartHistoryPolicy history_policy =
                                RestartHistoryPolicy::require_compatible) noexcept;
  // Add the accepted-flow convective bound without changing the caller's
  // explicit diffusion/acoustic bounds. No-op for a fixed-dt product.
  Status constrain_convective_time_limit(LocalTimeLimits& limits) noexcept;
  // Programmatic callers retain authority over explicit physical time scales.
  Status advance(LocalTimeLimits limits, DriverStepReport& report) noexcept;
  // Synchronous borrowed views. Consume before any advance/initialization,
  // another snapshot of the same kind, or destruction of the storage owner.
  // This conservative lifetime also applies after a rejected advance. Writers
  // validate view metadata, not an owner/epoch lease; stale use is unsupported.
  // Move construction transfers the owner without moving its storage. Move
  // assignment invalidates snapshots of the destination's former storage.
  // Asynchronous use requires a separately owned, budgeted copy of ALL data.
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
