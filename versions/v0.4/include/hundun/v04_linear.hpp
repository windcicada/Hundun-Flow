// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_execution.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

namespace hundun::v04 {

enum class LinearAlgorithm : std::uint8_t { pcg, fgmres, bicgstab };
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
enum class LinearPreconditionerClass : std::uint8_t {
  fixed_spd,
  fixed_general,
  flexible
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
  invalid_plan
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

class NumericState;
class HierarchyState;
class SolverWorkspace;

struct LinearOperatorCertificate {
  LinearIdentity identity{};
  PlanFingerprint collective_fingerprint{};
  Int3 local_shape{};
  LinearOperatorClass operator_class{LinearOperatorClass::nonsymmetric};
};

class LinearOperator {
 public:
  virtual LinearOperatorCertificate certificate() const noexcept = 0;
  virtual Status apply(FieldView x, FieldView y) const noexcept = 0;
  virtual ~LinearOperator() = default;
};

struct LinearPreconditionerCertificate {
  LinearIdentity identity{};
  PlanFingerprint collective_fingerprint{};
  LinearPreconditionerClass preconditioner_class{
      LinearPreconditionerClass::fixed_general};
};

class LinearPreconditioner {
 public:
  virtual LinearPreconditionerCertificate certificate() const noexcept = 0;
  virtual Status apply(ConstFieldView input, FieldView output,
                       std::uint32_t iteration) noexcept = 0;
  virtual ~LinearPreconditioner() = default;
};

struct LinearSolveControl {
  double absolute_tolerance{};
  double relative_tolerance{};
  std::uint32_t maximum_iterations{};
  std::uint32_t true_residual_interval{};
  std::uint32_t restart{};
};

struct LinearSolveInvocation {
  ConstFieldView rhs{};
  FieldView solution{};
  LinearIdentity expected_identity{};
  LinearSolveControl control{};
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
  int lowest_failing_rank{-1};
};

struct LinearReductionCounters {
  std::uint64_t calls{};
  std::uint64_t scalars{};
  std::uint64_t logical_bytes{};
  std::uint64_t tree_messages{};
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
  Status consensus(Status local_status) noexcept;
  Status consensus_contract(PlanFingerprint local_fingerprint) noexcept;
  int lowest_failing_rank() const noexcept;
  LinearReductionCounters counters() const noexcept;
  std::uintptr_t send_storage_address() const noexcept;
  std::uintptr_t receive_storage_address() const noexcept;
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

 private:
  friend LinearIdentity compose_linear_identity(
      const SymbolicPlan&, const NumericState&, const HierarchyState&,
      const SolverWorkspace&) noexcept;
  LinearWorkspaceRequirements requirements_{};
  FieldView vector_bundle_{};
  FieldView scalar_buffer_{};
  RevisionToken vector_revisions_[std::numeric_limits<std::uint8_t>::max()]{};
  RevisionToken next_vector_revision_{};
  PlanFingerprint fingerprint_{};
  RevisionToken binding_identity_{};
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
