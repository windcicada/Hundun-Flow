// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_flow.hpp"

#include <cstddef>
#include <cstdint>

namespace hundun::v04 {

// The production route keeps the exact anchored Cartesian operator and uses
// the existing scalable Native MG lifecycle.  PCG cannot legally consume the
// current Native MG certificate (fixed_general/prepared_batch), so the PCG
// route is deliberately an independent fixed-SPD Jacobi oracle for small
// verification problems; it is not a disguised production fallback.
enum class FreshStartProjectionLinearRoute : std::uint8_t {
  native_mg_fgmres,
  jacobi_pcg_oracle,
};

struct FreshStartKinematicProjectionSpec {
  MPI_Comm communicator{MPI_COMM_NULL};
  const CartesianGeometryPlan *geometry{};
  const CartesianKernelPlan *kernels{};
  const BoundaryPlan *boundary{};
  MeshPatch patch{};
  MgDomainActivityView activity{};
  FreshStartProjectionLinearRoute route{
      FreshStartProjectionLinearRoute::native_mg_fgmres};
  LinearSolveControl solve{};
  MgHierarchyPolicy mg_policy{};
  double compatibility_absolute_tolerance{1.0e-13};
  double compatibility_relative_tolerance{1.0e-12};
  double continuity_absolute_tolerance{1.0e-12};
  double continuity_relative_tolerance{1.0e-10};
};

// All services are borrowed and must outlive the plan.  operator_halo is
// reserved for the Krylov vector field.  correction_halo is a distinct halo
// for chi after the solve and before candidate reconstruction.  Native-MG
// services are mandatory only on the production route.
struct FreshStartKinematicProjectionServices {
  HaloEngine *operator_halo{};
  StageId operator_halo_stage{};
  FieldId krylov_field{};
  HaloEngine *correction_halo{};
  StageId correction_halo_stage{};
  SolverWorkspace *solver_workspace{};
  ReductionEngine *reductions{};
  MgRuntimeServices mg{};
};

// physical_mobility is the unanchored a_f used by the physical U/phi update.
// solver_coefficient is a separate anchored graph: every gauge-anchor row is
// identity and every off-diagonal incident on an anchor is removed.  Keeping
// both sets of faces is essential; a gauge choice must never become an IBM
// wall or alter the physical flux Jacobian.
struct FreshStartKinematicProjectionWorkspace {
  FieldView chi{};
  FieldView rhs{};
  FieldView diagonal{};
  FaceFieldView x_physical_mobility{};
  FaceFieldView y_physical_mobility{};
  FaceFieldView z_physical_mobility{};
  FaceFieldView x_solver_coefficient{};
  FaceFieldView y_solver_coefficient{};
  FaceFieldView z_solver_coefficient{};
  FieldView candidate_velocity{};
  FaceFluxView candidate_mass_flux{};
};

struct FreshStartKinematicProjectionInput {
  ConstFieldView density{};
  // Fresh starts have three velocity replicas but only accepted_n supplies
  // the kinematic base.  A successful joint commit publishes the projected
  // velocity into all three, preventing an incompatible BDF history from
  // being manufactured immediately after initialization.
  ConstFieldView velocity{};
  ConstFieldView velocity_accepted_n_minus_one{};
  ConstFieldView velocity_trial{};
  ConstFaceFluxView mass_flux{};
  RevisionToken state{};
};

// Closed-form RED/Jacobian/Schur authority for the cold projection:
//
//   R(U,phi) = D(phi),
//   J_chi    = D[-a_f jump_B(.)],
//   S_anchor = anchored(-J_chi),
//   U*       = U - G_IBM chi,
//   phi*     = phi - a_f jump_B(chi).
//
// `chi_writes_pressure == false` is an executable interface property: no
// method in this module accepts a pressure destination.
struct FreshStartKinematicProjectionRedCertificate {
  PlanFingerprint plan{};
  PlanFingerprint boundary_semantics{};
  PlanFingerprint activity_graph{};
  PlanFingerprint component_graph{};
  PlanFingerprint physical_flux_jacobian{};
  PlanFingerprint anchored_schur{};
  std::uint64_t active_cells{};
  // Compile-time uint64 payload used to reconcile rank-local components.
  // This is the maximum per-rank peak; the following field is the global sum
  // of those peaks.  Both contain only neighbour-boundary records and
  // distributed-component owner routes.  Interior and purely local
  // components never enter a collective payload.
  std::uint64_t component_collective_payload_u64{};
  std::uint64_t component_collective_global_payload_u64{};
  std::uint32_t component_count{};
  std::uint32_t anchored_component_count{};
  bool physical_mobility_separate{};
  bool per_component_compatibility{};
  bool global_minimum_gid_anchors{};
  bool mg_null_space_none{};
  bool ibm_gradient_from_face_activity{};
  bool cut_face_zero_mobility{};
  bool inactive_velocity_stationary_wall_zero{};
  bool component_collective_volume_independent{};
  bool runtime_workspace_preallocated{};
  bool no_immersed_bitwise_bypass{};
  bool three_layer_joint_commit{};
  bool chi_writes_pressure{};

  bool valid() const noexcept {
    return plan != 0U && boundary_semantics != 0U && activity_graph != 0U &&
           component_graph != 0U && physical_flux_jacobian != 0U &&
           anchored_schur != 0U && active_cells != 0U &&
           component_count != 0U &&
           anchored_component_count <= component_count &&
           physical_mobility_separate && per_component_compatibility &&
           global_minimum_gid_anchors && mg_null_space_none &&
           ibm_gradient_from_face_activity && cut_face_zero_mobility &&
           inactive_velocity_stationary_wall_zero &&
           component_collective_volume_independent &&
           runtime_workspace_preallocated && no_immersed_bitwise_bypass &&
           three_layer_joint_commit && !chi_writes_pressure;
  }
};

class FreshStartKinematicProjectionPreparedCertificate {
public:
  bool valid() const noexcept;
  std::uint32_t component_count() const noexcept { return component_count_; }
  std::uint32_t anchor_count() const noexcept { return anchor_count_; }
  double maximum_compatibility_defect() const noexcept {
    return maximum_compatibility_defect_;
  }
  PlanFingerprint lineage() const noexcept { return lineage_; }

private:
  friend class FreshStartKinematicProjectionPlan;
  const void *issuer_{};
  PlanFingerprint lineage_{};
  PlanFingerprint numeric_{};
  RevisionToken state_{};
  std::uint32_t component_count_{};
  std::uint32_t anchor_count_{};
  double maximum_compatibility_defect_{};
};

class FreshStartKinematicProjectionSolvedCertificate {
public:
  bool valid() const noexcept;
  const LinearSolveResult &result() const noexcept { return result_; }
  PlanFingerprint lineage() const noexcept { return lineage_; }

private:
  friend class FreshStartKinematicProjectionPlan;
  const void *issuer_{};
  PlanFingerprint prepared_lineage_{};
  PlanFingerprint lineage_{};
  RevisionToken chi_revision_{};
  LinearSolveResult result_{};
};

class FreshStartKinematicProjectionCandidateCertificate {
public:
  bool valid() const noexcept;
  double initial_continuity_maximum() const noexcept {
    return initial_continuity_maximum_;
  }
  double final_continuity_maximum() const noexcept {
    return final_continuity_maximum_;
  }
  double final_continuity_l2() const noexcept { return final_continuity_l2_; }
  PlanFingerprint lineage() const noexcept { return lineage_; }

private:
  friend class FreshStartKinematicProjectionPlan;
  const void *issuer_{};
  PlanFingerprint solved_lineage_{};
  PlanFingerprint candidate_state_{};
  PlanFingerprint candidate_flux_{};
  PlanFingerprint lineage_{};
  double initial_continuity_maximum_{};
  double final_continuity_maximum_{};
  double final_continuity_l2_{};
  bool fixed_flux_bitwise_unchanged_{};
  bool joint_candidate_{};
};

class FreshStartKinematicProjectionPlan {
public:
  FreshStartKinematicProjectionPlan() noexcept = default;
  ~FreshStartKinematicProjectionPlan() noexcept;
  FreshStartKinematicProjectionPlan(const FreshStartKinematicProjectionPlan &) =
      delete;
  FreshStartKinematicProjectionPlan &
  operator=(const FreshStartKinematicProjectionPlan &) = delete;
  FreshStartKinematicProjectionPlan(
      FreshStartKinematicProjectionPlan &&) noexcept;
  FreshStartKinematicProjectionPlan &
  operator=(FreshStartKinematicProjectionPlan &&) noexcept;

  static Status compile(const FreshStartKinematicProjectionSpec &spec,
                        FreshStartKinematicProjectionServices services,
                        FreshStartKinematicProjectionWorkspace workspace,
                        FreshStartKinematicProjectionPlan &out) noexcept;

  const FreshStartKinematicProjectionRedCertificate &red() const noexcept;

  Status prepare(
      const FreshStartKinematicProjectionInput &input,
      FreshStartKinematicProjectionPreparedCertificate &certificate) noexcept;

  Status solve(const FreshStartKinematicProjectionPreparedCertificate &prepared,
               FreshStartKinematicProjectionSolvedCertificate &certificate,
               ResourceCounters *resources = nullptr) noexcept;

  Status audit(
      const FreshStartKinematicProjectionSolvedCertificate &solved,
      FreshStartKinematicProjectionCandidateCertificate &certificate) noexcept;

  // Revalidates the complete candidate transaction collectively.  All
  // fallible work ends before the first store; the final region performs only
  // deterministic copies of U and phi.  Pressure is intentionally absent.
  Status
  commit(const FreshStartKinematicProjectionCandidateCertificate &candidate,
         Span<FieldView> velocity_layers, FaceFluxView mass_flux) noexcept;

private:
  struct Impl;
  void release() noexcept;
  Impl *implementation_{};
};

} // namespace hundun::v04
