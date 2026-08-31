// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_execution.hpp"
#include "hundun/v04_mesh.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

namespace hundun::v04 {

class HaloEngine;
class ReductionEngine;
class PisoPressureSolveEpoch;
struct LinearSolveResult;
struct ResourceCounters;
namespace detail {
struct NativeMgTestAccess;
}

inline constexpr std::size_t kLinearRecycleMaximumDirections = 4U;

enum class ReductionMode : std::uint8_t {
  reproducible_tree,
  mpi_allreduce
};
enum class LinearBackend : std::uint8_t {
  native_cartesian,
  hypre_struct
};
enum class LinearLocation : std::uint8_t {
  cell,
  face_x,
  face_y,
  face_z
};
enum class HierarchyUpdate : std::uint8_t {
  retain,
  refresh_numeric,
  rebuild
};
enum class LinearOperatorClass : std::uint8_t { spd, nonsymmetric };
// Operator validation and arithmetic failures are rank-local by default.
// A callback that already performed its own collective agreement may attach
// matching failure provenance so Krylov does not repeat that agreement and
// erase the callback-owned lowest failing rank.
enum class LinearOperatorStatusScope : std::uint8_t {
  rank_local,
  collective
};
enum class LinearPreconditionerClass : std::uint8_t {
  fixed_spd,
  fixed_general,
  flexible
};
// The default is deliberately conservative: a preconditioner callback may
// return a rank-local status, so the Krylov driver must perform its outer
// status consensus.  `collective` is a stronger callback contract: the
// callback owns that consensus and returns the same status on every rank.
enum class LinearPreconditionerStatusScope : std::uint8_t {
  rank_local,
  collective
};
// A per-call checked callback validates its input/output contract on every
// application.  A prepared batch callback validates a solver-owned, immutable
// set of vector slots once before the solve and then returns a collective
// status from each hot application.
enum class LinearPreconditionerApplyLifecycle : std::uint8_t {
  per_call_checked,
  prepared_batch
};
enum class LinearTermination : std::uint8_t {
  converged,
  zero_rhs,
  maximum_iterations,
  breakdown,
  non_finite,
  operator_failure,
  preconditioner_failure,
  collective_failure,
  invalid_plan,
  convergence_audit_failure
};
enum class CoarseningKind : std::uint8_t {
  full_xyz,
  semi_xy,
  semi_xz,
  semi_yz,
  x_only,
  y_only,
  z_only
};
enum class LineRelaxationDirection : std::uint8_t { none, x, y, z };
enum class MgBoundaryKind : std::uint8_t {
  dirichlet,
  neumann,
  periodic
};
enum class MgNullSpace : std::uint8_t { none, constant };
enum class MgPointSmootherKind : std::uint8_t {
  red_black,
  chebyshev_jacobi
};
enum class MgCycleKind : std::uint8_t { v_cycle, f_cycle };
enum class MgOperatorClass : std::uint8_t {
  general,
  symmetric_diagonally_dominant_m_matrix
};

struct MgHierarchyPolicy {
  double anisotropy_threshold{4.0};
  double coefficient_change_rebuild_ratio{0.25};
  std::uint8_t pre_sweeps{2U};
  std::uint8_t post_sweeps{2U};
  std::uint8_t maximum_levels{32U};
  std::uint16_t coarse_sweeps{24U};
  std::uint8_t minimum_coarse_extent{3U};
  std::uint16_t line_relaxation_maximum_extent{4096U};
  MgPointSmootherKind point_smoother{MgPointSmootherKind::red_black};
  MgCycleKind cycle{MgCycleKind::v_cycle};
  double chebyshev_lower_spectrum_fraction{0.3};
};

struct MgBoundarySet {
  MgBoundaryKind x_min{MgBoundaryKind::dirichlet};
  MgBoundaryKind x_max{MgBoundaryKind::dirichlet};
  MgBoundaryKind y_min{MgBoundaryKind::dirichlet};
  MgBoundaryKind y_max{MgBoundaryKind::dirichlet};
  MgBoundaryKind z_min{MgBoundaryKind::dirichlet};
  MgBoundaryKind z_max{MgBoundaryKind::dirichlet};
};

struct MgCoefficientIdentity {
  RevisionToken revision{};
  PlanFingerprint fingerprint{};
  double maximum_relative_change{};
};

struct MgPlanCounters {
  std::uint64_t symbolic_builds{};
  std::uint64_t numeric_refreshes{};
  std::uint64_t hierarchy_rebuilds{};
  std::uint64_t applications{};
  std::uint64_t blocking_collectives{};
  std::uint64_t collective_logical_bytes{};
};

struct MgLevelView {
  Int3 global_shape{};
  Int3 local_shape{};
  CoarseningKind coarsening{CoarseningKind::full_xyz};
  std::uint8_t line_axis_mask{};
};

struct MgCoefficientViews {
  ConstFieldView diagonal{};
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
};

// Optional immutable activity map for an embedded/internal boundary.  Values
// are exactly zero for inactive cells/faces and one for active cells/faces;
// empty spans select the ordinary full Cartesian domain.  The MG hierarchy
// uses the map only as an approximate preconditioner operator.  The Krylov
// LinearOperator remains the exact discretization authority.
struct MgDomainActivityView {
  Span<const std::uint8_t> cells{};
  Span<const std::uint8_t> x_faces{};
  Span<const std::uint8_t> y_faces{};
  Span<const std::uint8_t> z_faces{};
  PlanFingerprint local_fingerprint{};
  PlanFingerprint collective_fingerprint{};
};

inline constexpr std::size_t kMgMaximumLevels = 32U;

enum class MgWorkspaceSlot : std::uint8_t {
  solution,
  rhs,
  residual,
  temporary
};

struct MgWorkspaceLevelRequirements {
  Int3 global_shape{};
  MeshPatch patch{};
  CoarseningKind coarsening{CoarseningKind::full_xyz};
  std::uint8_t coarsen_axis_mask{};
  std::uint8_t line_axis_mask{};
  std::size_t offset_doubles{};
  std::size_t slot_stride_doubles{};
  std::size_t stride_y{};
  std::size_t stride_z{};
};

struct MgWorkspaceRequirements {
  std::array<MgWorkspaceLevelRequirements, kMgMaximumLevels> levels{};
  std::size_t level_count{};
  Int3 arena_shape{};
  std::uint8_t ghost_width{1U};
  std::uint8_t slots_per_level{4U};
  std::size_t total_doubles{};
  RevisionToken execution_revision{};
  PlanFingerprint collective_fingerprint{};
  PlanFingerprint fingerprint{};
};

struct MgWorkspaceFieldIds {
  FieldId vector_bundle{};
};

struct SymbolicSpec {
  std::uint32_t operator_kind{};
  LinearLocation location{LinearLocation::cell};
  RevisionToken topology{};
  RevisionToken boundary_layout{};
  RevisionToken partition{};
  RevisionToken eb_interface{};
  RevisionToken stencil{};
  LinearBackend backend{LinearBackend::native_cartesian};
};

struct CoefficientRevisions {
  RevisionToken diagonal{};
  RevisionToken off_diagonal{};
  RevisionToken time{};
  RevisionToken material_transport{};
  RevisionToken numeric_boundary{};
  RevisionToken constraint{};
};

struct HierarchyPolicyIdentity {
  RevisionToken coarsening{};
  RevisionToken transfer_smoother{};
  RevisionToken policy_epoch{};
  LinearBackend backend{LinearBackend::native_cartesian};
};

struct LinearLifecycleCounters {
  std::uint64_t symbolic_builds{};
  std::uint64_t numeric_refills{};
  std::uint64_t hierarchy_refreshes{};
  std::uint64_t hierarchy_rebuilds{};
  std::uint64_t workspace_bindings{};
  std::uint64_t workspace_replacements{};
};

struct LinearWorkspaceRequirements {
  LinearAlgorithm algorithm{LinearAlgorithm::pcg};
  Int3 maximum_shape{};
  std::uint8_t ghost_width{};
  std::uint32_t maximum_restart{};
  ReductionMode reduction_mode{ReductionMode::mpi_allreduce};
  RevisionToken execution_revision{};
  std::uint8_t vector_slots{};
  std::size_t scalar_doubles{};
  std::size_t reduction_capacity{};
  PlanFingerprint fingerprint{};
};

struct LinearWorkspaceFieldIds {
  FieldId vector_bundle{};
  FieldId scalar_buffer{};
};

struct LinearIdentity {
  PlanFingerprint symbolic{};
  PlanFingerprint numeric{};
  PlanFingerprint hierarchy{};
  PlanFingerprint workspace{};
  PlanFingerprint fingerprint{};
};

struct NativeCartesianMgSpec {
  MPI_Comm communicator{MPI_COMM_NULL};
  const CartesianGeometryPlan* geometry{};
  MeshPatch patch{};
  MgBoundarySet boundaries{};
  MgNullSpace null_space{MgNullSpace::none};
  MgOperatorClass operator_class{MgOperatorClass::general};
  MgHierarchyPolicy policy{};
  MgCorrectionScaling correction_scaling{
      MgCorrectionScaling::residual_minimizing};
  LinearIdentity identity{};
  MgCoefficientIdentity coefficients{};
  MgDomainActivityView activity{};
};

class NumericState;
class HierarchyState;
class SolverWorkspace;
class MgWorkspace;

struct LinearPreconditionerBatchDescriptor {
  const SolverWorkspace* workspace{};
  Int3 shape{};
  std::uint8_t input_slot_begin{};
  std::uint8_t output_slot_begin{};
  std::uint8_t slot_count{};
  std::uint32_t maximum_applications{};
};

// The ticket is intentionally opaque to callers.  It is a per-solve value
// produced by prepare_batch() and accepted only by apply_prepared().
class LinearPreconditionerBatchTicket {
 public:
  LinearPreconditionerBatchTicket() noexcept = default;

 private:
  friend class LinearPreconditioner;
  friend class NativeCartesianMgPlan;
  std::uintptr_t owner{};
  const SolverWorkspace* workspace{};
  Int3 shape{};
  std::uint8_t input_slot_begin{};
  std::uint8_t output_slot_begin{};
  std::uint8_t slot_count{};
  std::uint32_t maximum_applications{};
  PlanFingerprint preconditioner_fingerprint{};
  RevisionToken preconditioner_generation{};
  std::uint64_t preconditioner_application_base{};
};

struct LinearOperatorCertificate {
  LinearIdentity identity{};
  PlanFingerprint collective_fingerprint{};
  Int3 local_shape{};
  LinearOperatorClass operator_class{LinearOperatorClass::nonsymmetric};
};

// This metadata is valid only for the exact Status returned by the most
// recent failed apply and remains stable until the next apply.  A default or
// non-matching record is conservatively treated as rank-local.
struct LinearOperatorFailureProvenance {
  Status status{};
  LinearOperatorStatusScope status_scope{
      LinearOperatorStatusScope::rank_local};
  int lowest_failing_rank{-1};
};

class LinearOperator {
 public:
  virtual LinearOperatorCertificate certificate() const noexcept = 0;
  virtual Status apply(FieldView x, FieldView y) const noexcept = 0;
  virtual LinearOperatorFailureProvenance failure_provenance() const
      noexcept {
    return {};
  }
  virtual ~LinearOperator() = default;
};

struct LinearPreconditionerCertificate {
  LinearIdentity identity{};
  PlanFingerprint collective_fingerprint{};
  LinearPreconditionerClass preconditioner_class{
      LinearPreconditionerClass::fixed_general};
  LinearPreconditionerStatusScope status_scope{
      LinearPreconditionerStatusScope::rank_local};
  LinearPreconditionerApplyLifecycle apply_lifecycle{
      LinearPreconditionerApplyLifecycle::per_call_checked};
};

class LinearPreconditioner {
 public:
  virtual LinearPreconditionerCertificate certificate() const noexcept = 0;
  virtual Status apply(ConstFieldView input, FieldView output,
                       std::uint32_t iteration) noexcept = 0;
  virtual Status prepare_batch(
      const LinearPreconditionerBatchDescriptor&,
      LinearPreconditionerBatchTicket&) noexcept {
    return {StatusCode::invalid_plan, 603U};
  }
  virtual Status apply_prepared(
      ConstFieldView, FieldView, std::uint32_t,
      const LinearPreconditionerBatchTicket&) noexcept {
    return {StatusCode::invalid_plan, 603U};
  }
  virtual ~LinearPreconditioner() = default;

 protected:
  static void issue_batch_ticket(
      LinearPreconditionerBatchTicket& ticket,
      const LinearPreconditioner* owner,
      const LinearPreconditionerBatchDescriptor& descriptor,
      PlanFingerprint preconditioner_fingerprint = 0U,
      RevisionToken preconditioner_generation = 0U) noexcept {
    ticket.owner = reinterpret_cast<std::uintptr_t>(owner);
    ticket.workspace = descriptor.workspace;
    ticket.shape = descriptor.shape;
    ticket.input_slot_begin = descriptor.input_slot_begin;
    ticket.output_slot_begin = descriptor.output_slot_begin;
    ticket.slot_count = descriptor.slot_count;
    ticket.maximum_applications = descriptor.maximum_applications;
    ticket.preconditioner_fingerprint = preconditioner_fingerprint;
    ticket.preconditioner_generation = preconditioner_generation;
  }
};

struct LinearSolveControl {
  double absolute_tolerance{};
  double relative_tolerance{};
  std::uint32_t maximum_iterations{};
  std::uint32_t true_residual_interval{};
  std::uint32_t restart{};
};

struct LinearConvergenceAuditCertificate {
  PlanFingerprint collective_fingerprint{};
};

struct LinearConvergenceFailureProvenance {
  bool valid{};
  std::uint64_t global_cell{};
  Int3 global_index{-1, -1, -1};
  int owner_rank{-1};
  std::uint8_t pressure_activity{};
  double predecessor_application_scale{1.0};
  double predecessor_maximum_depletion{};
  double density{};
  double psi{};
  double raw_correction{};
  double depletion{};
  double bdf_storage{};
  double flux_divergence{};
  double reconstructed_rhs{};
  double sealed_rhs{};
  double rhs_absolute_mismatch{};
  double rhs_relative_mismatch{};
};

struct LinearConvergenceAuditResult {
  bool accepted{};
  double metric{};
  double limit{};
  double application_scale{1.0};
  double unscaled_metric{};
  double maximum_depletion{};
  double operator_parity_error{};
  // A rejected supplemental criterion normally asks the Krylov method to
  // continue.  Set this only when further iterations of the unchanged linear
  // system cannot alter the application-level decision.
  bool terminal_rejection{};
  LinearConvergenceFailureProvenance failure_provenance{};
};

// Optional native convergence obligation evaluated only after the canonical
// FP64 true-residual criterion has passed.  The callback owns any reductions
// needed to return one collective metric/decision on every rank.  It may make
// a solver continue or terminate an application-level rejection, but it can
// never replace or relax the true-residual gate.
class LinearConvergenceAudit {
 public:
  virtual LinearConvergenceAuditCertificate certificate() const noexcept = 0;
  virtual Status evaluate(ConstFieldView solution,
                          ConstFieldView true_residual,
                          ReductionEngine& reductions,
                          LinearConvergenceAuditResult& result) noexcept = 0;
  virtual ~LinearConvergenceAudit() = default;
};

struct LinearSolveInvocation {
  ConstFieldView rhs{};
  FieldView solution{};
  LinearIdentity expected_identity{};
  LinearSolveControl control{};
  LinearConvergenceAudit* convergence_audit{};
};

struct LinearSolveResult {
  Status status{};
  LinearTermination termination{LinearTermination::invalid_plan};
  std::uint32_t iterations{};
  double initial_true_residual{};
  double final_true_residual{};
  double recursive_residual{};
  std::uint64_t reduction_calls{};
  std::uint64_t operator_applies{};
  std::uint64_t preconditioner_applies{};
  std::uint64_t norm_breakdown_restarts{};
  std::uint64_t convergence_audits{};
  std::uint64_t convergence_rejections{};
  double final_convergence_metric{};
  double convergence_limit{};
  double convergence_application_scale{1.0};
  double convergence_unscaled_metric{};
  double convergence_maximum_depletion{};
  double convergence_operator_parity_error{};
  LinearConvergenceFailureProvenance convergence_failure_provenance{};
  int lowest_failing_rank{-1};
  std::uint64_t recycle_offered_directions{};
  std::uint64_t recycle_retained_directions{};
  std::uint64_t recycle_operator_applies{};
  std::uint64_t recycle_reduction_calls{};
  bool recycle_projection_attempted{};
  bool recycle_projection_accepted{};
  double recycle_projected_true_residual{};
  std::uint64_t recycle_cycle_corrections{};
  std::uint64_t recycle_capture_vector_passes{};
  std::uint64_t recycle_capture_cycle_attempts{};
  std::uint64_t recycle_capture_reduction_calls{};
  std::uint64_t recycle_capture_blocking_operations{};
};

struct LinearReductionCounters {
  std::uint64_t calls{};
  std::uint64_t scalars{};
  std::uint64_t logical_bytes{};
  std::uint64_t tree_messages{};
  std::uint64_t blocking_operations{};
  std::uint64_t wall_nanoseconds{};
};

inline constexpr std::size_t kReductionMaximumLocationPayload = 5U;

// Deterministic max-location reduction.  Equal values select the smallest
// global location and then the smallest rank, so the result is independent
// of the MPI partition.  Payload travels with the selected record.
struct ReductionMaximumLocation {
  bool valid{};
  double value{};
  std::uint64_t global_location{};
  std::int32_t rank{-1};
  std::array<double, kReductionMaximumLocationPayload> payload{};
};

class ReductionEngine {
 public:
  ReductionEngine() noexcept = default;
  ~ReductionEngine() noexcept;
  ReductionEngine(const ReductionEngine&) = delete;
  ReductionEngine& operator=(const ReductionEngine&) = delete;
  ReductionEngine(ReductionEngine&&) noexcept;
  ReductionEngine& operator=(ReductionEngine&&) noexcept;

  static Status compile(MPI_Comm communicator, ReductionMode mode,
                        std::size_t maximum_scalars,
                        ReductionEngine& out) noexcept;
  Status checked_sum(Span<const double> local, Span<double> global,
                     Status local_status = {}) noexcept;
  Status checked_max(Span<const double> local, Span<double> global,
                     Status local_status = {}) noexcept;
  Status checked_max_locations(
      Span<const ReductionMaximumLocation> local,
      Span<ReductionMaximumLocation> global,
      Status local_status = {}) noexcept;
  Status consensus(Status local_status) noexcept;
  Status consensus_contract(PlanFingerprint local_fingerprint) noexcept;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  Status arm_checked_sum_fault_for_test(std::uint64_t ordinal,
                                        int rank) noexcept;
  void clear_checked_sum_fault_for_test() noexcept;
#endif
  Status validate_communicator(MPI_Comm communicator) const noexcept;
  int lowest_failing_rank() const noexcept;
  LinearReductionCounters counters() const noexcept;
  std::uintptr_t send_storage_address() const noexcept;
  std::uintptr_t receive_storage_address() const noexcept;
  std::uintptr_t instance_identity() const noexcept {
    return reinterpret_cast<std::uintptr_t>(implementation_);
  }
  std::size_t capacity() const noexcept;
  ReductionMode mode() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

class SymbolicPlan {
 public:
  static Status compile(SymbolicSpec spec, SymbolicPlan& out,
                        LinearLifecycleCounters* counters = nullptr) noexcept;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  RevisionToken generation() const noexcept { return generation_; }
  bool valid() const noexcept { return fingerprint_ != 0U; }

 private:
  friend class NumericState;
  friend class HierarchyState;
  friend LinearIdentity compose_linear_identity(
      const SymbolicPlan&, const NumericState&, const HierarchyState&,
      const SolverWorkspace&) noexcept;
  SymbolicSpec spec_{};
  PlanFingerprint fingerprint_{};
  RevisionToken generation_{};
};

class NumericState {
 public:
  Status refill(const SymbolicPlan& symbolic,
                CoefficientRevisions revisions,
                PlanFingerprint content_fingerprint,
                LinearLifecycleCounters* counters = nullptr) noexcept;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  RevisionToken generation() const noexcept { return generation_; }
  bool valid_for(const SymbolicPlan& symbolic) const noexcept;

 private:
  friend class HierarchyState;
  friend LinearIdentity compose_linear_identity(
      const SymbolicPlan&, const NumericState&, const HierarchyState&,
      const SolverWorkspace&) noexcept;
  CoefficientRevisions revisions_{};
  PlanFingerprint symbolic_fingerprint_{};
  RevisionToken symbolic_generation_{};
  PlanFingerprint content_fingerprint_{};
  PlanFingerprint fingerprint_{};
  RevisionToken generation_{};
};

class HierarchyState {
 public:
  Status update(const SymbolicPlan& symbolic, const NumericState& numeric,
                HierarchyPolicyIdentity policy, HierarchyUpdate decision,
                PlanFingerprint content_fingerprint,
                LinearLifecycleCounters* counters = nullptr) noexcept;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  RevisionToken generation() const noexcept { return generation_; }
  bool valid_for(const SymbolicPlan& symbolic,
                 const NumericState& numeric) const noexcept;

 private:
  friend LinearIdentity compose_linear_identity(
      const SymbolicPlan&, const NumericState&, const HierarchyState&,
      const SolverWorkspace&) noexcept;
  HierarchyPolicyIdentity policy_{};
  PlanFingerprint symbolic_fingerprint_{};
  RevisionToken symbolic_generation_{};
  PlanFingerprint numeric_fingerprint_{};
  RevisionToken numeric_generation_{};
  PlanFingerprint content_fingerprint_{};
  PlanFingerprint fingerprint_{};
  RevisionToken generation_{};
};

class SolverWorkspace {
 public:
  static Status bind(const LinearWorkspaceRequirements& requirements,
                     FieldView vector_bundle, FieldView scalar_buffer,
                     SolverWorkspace& out,
                     LinearLifecycleCounters* counters = nullptr) noexcept;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  RevisionToken binding_identity() const noexcept { return binding_identity_; }
  std::uintptr_t vector_storage_address() const noexcept {
    return reinterpret_cast<std::uintptr_t>(vector_bundle_.base);
  }
  std::uintptr_t scalar_storage_address() const noexcept {
    return reinterpret_cast<std::uintptr_t>(scalar_buffer_.base);
  }
  const LinearWorkspaceRequirements& requirements() const noexcept {
    return requirements_;
  }
  bool valid_for(const LinearWorkspaceRequirements& requirements) const
      noexcept;
  FieldView vector(std::uint8_t slot, Int3 active_shape) const noexcept;
  Status revise_vector(std::uint8_t slot) noexcept;
  bool overlaps_storage(ConstFieldView view) const noexcept;
  bool overlaps_storage(FieldView view) const noexcept;
  Span<double> scalars(std::size_t offset, std::size_t count) const noexcept;

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  Status recycle_begin_capture_for_test(
      Int3 shape, PlanFingerprint source_identity) noexcept;
  Status recycle_begin_projection_for_test(
      Int3 shape, PlanFingerprint current_identity) noexcept;
  Status recycle_capture_cycle_start_for_test(
      ConstFieldView solution, ReductionEngine& reductions,
      Status prerequisite = {}) noexcept;
  Status recycle_capture_cycle_publish_for_test(
      ConstFieldView solution, ReductionEngine& reductions) noexcept;
  void recycle_clear_for_test() noexcept;
  std::size_t recycle_correction_count_for_test() const noexcept;
  std::uint64_t recycle_capture_vector_passes_for_test() const noexcept;
  std::uint64_t recycle_capture_cycle_attempts_for_test() const noexcept;
  std::uint64_t recycle_capture_reduction_calls_for_test() const noexcept;
  std::uint64_t recycle_capture_blocking_operations_for_test() const noexcept;
  std::uint64_t recycle_capture_cycle_corrections_for_test() const noexcept;
  std::uint8_t recycle_snapshot_slot_for_test() const noexcept;
  std::uint8_t recycle_correction_physical_slot_for_test(
      std::size_t index) const noexcept;
  ConstFieldView recycle_correction_for_test(
      std::size_t index, Int3 shape) const noexcept;
#endif

 private:
  friend class PisoPressureSolveEpoch;
  friend LinearSolveResult solve_fgmres(
      const LinearOperator&, LinearPreconditioner&,
      const LinearSolveInvocation&, SolverWorkspace&, ReductionEngine&,
      ResourceCounters*) noexcept;
  friend LinearIdentity compose_linear_identity(
      const SymbolicPlan&, const NumericState&, const HierarchyState&,
      const SolverWorkspace&) noexcept;

  struct RecycleState {
    bool capture_active{};
    bool cycle_active{};
    bool projection_pending{};
    Int3 shape{};
    PlanFingerprint source_identity{};
    PlanFingerprint current_identity{};
    std::uint8_t snapshot_slot{};
    std::uint8_t correction_count{};
    std::uint8_t oldest_correction{};
    std::array<std::uint8_t, kLinearRecycleMaximumDirections>
        correction_order{};
    std::uint64_t capture_vector_passes{};
    std::uint64_t capture_cycle_attempts{};
    std::uint64_t capture_reduction_calls{};
    std::uint64_t capture_blocking_operations{};
    std::uint64_t published_cycle_corrections{};
  };

  Status recycle_begin_capture(Int3 shape,
                               PlanFingerprint source_identity) noexcept;
  Status recycle_begin_projection(Int3 shape,
                                  PlanFingerprint current_identity) noexcept;
  void recycle_clear() noexcept;
  bool recycle_capture_active() const noexcept {
    return recycle_.capture_active;
  }
  Status recycle_capture_cycle_start(ConstFieldView solution,
                                     ReductionEngine& reductions,
                                     Status prerequisite = {}) noexcept;
  Status recycle_capture_cycle_publish(ConstFieldView solution,
                                       ReductionEngine& reductions) noexcept;
  std::size_t recycle_correction_count() const noexcept {
    return recycle_.correction_count;
  }
  std::uint64_t recycle_capture_vector_passes() const noexcept {
    return recycle_.capture_vector_passes;
  }
  std::uint64_t recycle_capture_reduction_calls() const noexcept {
    return recycle_.capture_reduction_calls;
  }
  std::uint64_t recycle_capture_cycle_attempts() const noexcept {
    return recycle_.capture_cycle_attempts;
  }
  std::uint64_t recycle_capture_blocking_operations() const noexcept {
    return recycle_.capture_blocking_operations;
  }
  std::uint64_t recycle_capture_cycle_corrections() const noexcept {
    return recycle_.published_cycle_corrections;
  }
  ConstFieldView recycle_correction(std::size_t index,
                                    Int3 shape) const noexcept;
  FieldView recycle_snapshot(Int3 shape) const noexcept;
  FieldView recycle_correction_storage(std::size_t index,
                                       Int3 shape) const noexcept;
  std::uint8_t recycle_pool_slot(std::size_t index) const noexcept;
  std::uint8_t recycle_snapshot_slot() const noexcept;
  std::uint8_t recycle_correction_logical_slot(
      std::size_t index) const noexcept;
  bool recycle_projection_pending() const noexcept {
    return recycle_.projection_pending;
  }
  void recycle_skip_projection(LinearSolveResult* result = nullptr) noexcept {
    if (result != nullptr && recycle_.projection_pending) {
      result->recycle_offered_directions = recycle_.correction_count;
      result->recycle_retained_directions = 0U;
      result->recycle_operator_applies = 0U;
      result->recycle_reduction_calls = 0U;
      result->recycle_projection_attempted = false;
      result->recycle_projection_accepted = false;
      result->recycle_projected_true_residual = 0.0;
    }
    recycle_.projection_pending = false;
  }
  void recycle_set_projection_result(LinearSolveResult& result,
                                     std::uint64_t offered,
                                     std::uint64_t retained,
                                     std::uint64_t operator_applies,
                                     std::uint64_t reduction_calls,
                                     bool attempted, bool accepted,
                                     double projected_residual) noexcept;

  LinearWorkspaceRequirements requirements_{};
  FieldView vector_bundle_{};
  FieldView scalar_buffer_{};
  RevisionToken vector_revisions_[std::numeric_limits<std::uint8_t>::max()]{};
  RevisionToken next_vector_revision_{};
  PlanFingerprint fingerprint_{};
  RevisionToken binding_identity_{};
  RecycleState recycle_{};
};

class MgWorkspace {
 public:
  static Status bind(const MgWorkspaceRequirements& requirements,
                     FieldView vector_bundle, MgWorkspace& out) noexcept;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  PlanFingerprint collective_fingerprint() const noexcept {
    return requirements_.collective_fingerprint;
  }
  RevisionToken binding_identity() const noexcept { return binding_identity_; }
  std::uintptr_t storage_address() const noexcept {
    return reinterpret_cast<std::uintptr_t>(raw_base_);
  }
  std::size_t storage_doubles() const noexcept {
    return requirements_.total_doubles;
  }
  std::size_t level_count() const noexcept {
    return requirements_.level_count;
  }
  const MgWorkspaceRequirements& requirements() const noexcept {
    return requirements_;
  }
  const MgWorkspaceLevelRequirements* level_requirements(
      std::size_t level) const noexcept;
  bool valid_for(const MgWorkspaceRequirements& requirements) const noexcept;
  FieldView level(std::size_t level, MgWorkspaceSlot slot) const noexcept;
  Status revise_level(std::size_t level, MgWorkspaceSlot slot) noexcept;
  bool overlaps_storage(ConstFieldView view) const noexcept;
  bool overlaps_storage(FieldView view) const noexcept;

 private:
  MgWorkspaceRequirements requirements_{};
  FieldView vector_bundle_{};
  double* raw_base_{};
  std::array<RevisionToken, kMgMaximumLevels * 4U> revisions_{};
  RevisionToken next_revision_{};
  PlanFingerprint fingerprint_{};
  RevisionToken binding_identity_{};
};

struct MgRuntimeServices {
  HaloEngine* finest_halo{};
  ReductionEngine* reductions{};
  MgWorkspace* workspace{};
  Span<HaloEngine* const> level_halos{};
};

class NativeCartesianMgPlan final : public LinearPreconditioner {
 public:
  NativeCartesianMgPlan() noexcept = default;
  ~NativeCartesianMgPlan() noexcept;
  NativeCartesianMgPlan(const NativeCartesianMgPlan&) = delete;
  NativeCartesianMgPlan& operator=(const NativeCartesianMgPlan&) = delete;
  NativeCartesianMgPlan(NativeCartesianMgPlan&&) noexcept;
  NativeCartesianMgPlan& operator=(NativeCartesianMgPlan&&) noexcept;
  static Status compile(const NativeCartesianMgSpec& spec,
                        MgRuntimeServices services,
                        MgCoefficientViews coefficients,
                        NativeCartesianMgPlan& out,
                        MgPlanCounters* counters = nullptr) noexcept;
  Status update_coefficients(LinearIdentity next_identity,
                             MgCoefficientIdentity identity,
                             MgCoefficientViews coefficients,
                             MgPlanCounters* counters = nullptr) noexcept;
  LinearPreconditionerCertificate certificate() const noexcept override;
  Status apply(ConstFieldView residual, FieldView correction,
               std::uint32_t iteration) noexcept override;
  Status prepare_batch(const LinearPreconditionerBatchDescriptor& descriptor,
                       LinearPreconditionerBatchTicket& ticket) noexcept override;
  Status apply_prepared(ConstFieldView residual, FieldView correction,
                        std::uint32_t iteration,
                        const LinearPreconditionerBatchTicket& ticket) noexcept override;
  std::size_t level_count() const noexcept;
  Status level(std::size_t index, MgLevelView& out) const noexcept;
  CoarseningKind finest_coarsening() const noexcept;
  std::uint8_t line_axis_mask() const noexcept;
  PlanFingerprint symbolic_fingerprint() const noexcept;
  PlanFingerprint numeric_fingerprint() const noexcept;
  PlanFingerprint hierarchy_fingerprint() const noexcept;
  RevisionToken generation() const noexcept;
  std::uintptr_t hierarchy_storage_address() const noexcept;
  std::uintptr_t workspace_storage_address() const noexcept;
  std::size_t workspace_doubles() const noexcept;
  double last_cycle_initial_residual() const noexcept;
  double last_cycle_final_residual() const noexcept;
  MgPlanCounters counters() const noexcept;
  int lowest_failing_rank() const noexcept;

 private:
  friend struct detail::NativeMgTestAccess;
  struct Impl;
  Status apply_impl(ConstFieldView residual, FieldView correction,
                    bool prepared) noexcept;
  void release() noexcept;
  Impl* implementation_{};
};

class HypreStructAdapter final : public LinearPreconditioner {
 public:
  HypreStructAdapter() noexcept = default;
  ~HypreStructAdapter() noexcept;
  HypreStructAdapter(const HypreStructAdapter&) = delete;
  HypreStructAdapter& operator=(const HypreStructAdapter&) = delete;
  HypreStructAdapter(HypreStructAdapter&&) noexcept;
  HypreStructAdapter& operator=(HypreStructAdapter&&) noexcept;
  static bool available() noexcept;
  static Status compile(const NativeCartesianMgSpec& spec,
                        MgRuntimeServices services,
                        HypreStructAdapter& out) noexcept;
  Status update_coefficients(LinearIdentity next_identity,
                             MgCoefficientIdentity identity,
                             MgCoefficientViews coefficients) noexcept;
  LinearPreconditionerCertificate certificate() const noexcept override;
  Status apply(ConstFieldView residual, FieldView correction,
               std::uint32_t iteration) noexcept override;
  PlanFingerprint numeric_fingerprint() const noexcept;
  PlanFingerprint fingerprint() const noexcept;
  std::uintptr_t native_handle_storage_address() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

Status make_linear_workspace_requirements(
    LinearAlgorithm algorithm, Int3 maximum_shape, std::uint8_t ghost_width,
    std::uint32_t maximum_restart, ReductionMode reduction_mode,
    RevisionToken execution_revision,
    LinearWorkspaceRequirements& out) noexcept;

Status register_linear_workspace(
    FieldRegistry& registry, std::string_view stable_prefix,
    const LinearWorkspaceRequirements& requirements,
    LinearWorkspaceFieldIds& out) noexcept;

Status make_mg_workspace_requirements(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    MeshPatch finest_patch, const MgHierarchyPolicy& policy,
    RevisionToken execution_revision,
    MgWorkspaceRequirements& out) noexcept;

Status register_mg_workspace(
    FieldRegistry& registry, std::string_view stable_prefix,
    const MgWorkspaceRequirements& requirements,
    MgWorkspaceFieldIds& out) noexcept;

LinearIdentity compose_linear_identity(
    const SymbolicPlan& symbolic, const NumericState& numeric,
    const HierarchyState& hierarchy,
    const SolverWorkspace& workspace) noexcept;

LinearSolveResult solve_pcg(
    const LinearOperator& linear_operator,
    LinearPreconditioner& preconditioner,
    const LinearSolveInvocation& invocation, SolverWorkspace& workspace,
    ReductionEngine& reductions,
    ResourceCounters* resources = nullptr) noexcept;
LinearSolveResult solve_fgmres(
    const LinearOperator& linear_operator,
    LinearPreconditioner& preconditioner,
    const LinearSolveInvocation& invocation, SolverWorkspace& workspace,
    ReductionEngine& reductions,
    ResourceCounters* resources = nullptr) noexcept;
LinearSolveResult solve_bicgstab(
    const LinearOperator& linear_operator,
    LinearPreconditioner& preconditioner,
    const LinearSolveInvocation& invocation, SolverWorkspace& workspace,
    ReductionEngine& reductions,
    ResourceCounters* resources = nullptr) noexcept;

}  // namespace hundun::v04
