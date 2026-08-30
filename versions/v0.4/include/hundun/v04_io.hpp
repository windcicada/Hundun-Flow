// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_execution.hpp"
#include "hundun/v04_linear.hpp"
#include "hundun/v04_mesh.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::v04 {

enum class RuntimeServiceKind : std::uint8_t {
  restart,
  visit,
  screen,
  monitor,
  evidence
};

struct SnapshotFieldSpec {
  FieldId field{};
  std::uint8_t components{};
};

struct RuntimeServiceCapacity {
  RuntimeServiceKind kind{RuntimeServiceKind::restart};
  StageId stage{};
  std::size_t maximum_snapshot_bytes_per_rank{};
  std::size_t maximum_staging_bytes_per_rank{};
  std::uint32_t maximum_collectives{};
};

class IoServicePlan {
 public:
  static Status compile(Span<const SnapshotFieldSpec> snapshot_fields,
                        Span<const RuntimeServiceCapacity> services,
                        std::size_t local_cells,
                        IoServicePlan& out) noexcept;

  Span<const SnapshotFieldSpec> snapshot_fields() const noexcept {
    return {snapshot_fields_.data(), snapshot_field_count_};
  }
  Span<const RuntimeServiceCapacity> services() const noexcept {
    return {services_.data(), service_count_};
  }
  std::size_t maximum_staging_bytes() const noexcept {
    return maximum_staging_bytes_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  static constexpr std::size_t kMaximumSnapshotFields = 64U;
  static constexpr std::size_t kMaximumServices = 5U;
  std::array<SnapshotFieldSpec, kMaximumSnapshotFields> snapshot_fields_{};
  std::array<RuntimeServiceCapacity, kMaximumServices> services_{};
  std::size_t snapshot_field_count_{};
  std::size_t service_count_{};
  std::size_t maximum_staging_bytes_{};
  PlanFingerprint fingerprint_{};
};

enum class RestartFieldRole : std::uint8_t {
  velocity,
  pressure_perturbation,
  pressure_absolute,
  enthalpy,
  independent_species,
  transported_scalar
};

struct RestartFieldView {
  RestartFieldRole role{RestartFieldRole::velocity};
  ConstFieldView values{};
};

struct RestartSnapshot {
  Int3 global_cells{};
  MeshPatch patch{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  PlanFingerprint geometry{};
  double time{};
  double dt{};
  double pressure_reference{};
  std::uint64_t step{};
  std::uint64_t controller_state{};
  Span<const RestartFieldView> fields{};
  ConstFaceFluxView final_mass_flux{};
};

struct RestartWriteOptions {
  std::uint32_t keep_last{1U};
};

struct RestartExpectedField {
  RestartFieldRole role{RestartFieldRole::velocity};
  FieldId field{};
  std::uint8_t components{};
};

struct RestartExpected {
  Int3 global_cells{};
  MeshPatch target_patch{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  PlanFingerprint geometry{};
  Span<const RestartExpectedField> fields{};
};

struct RestartImageField {
  RestartFieldRole role{RestartFieldRole::velocity};
  FieldId field{};
  std::uint8_t components{};
  std::vector<double> values;
};

struct RestartImage {
  Int3 global_cells{};
  MeshPatch patch{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  PlanFingerprint geometry{};
  double time{};
  double dt{};
  double pressure_reference{};
  std::uint64_t step{};
  std::uint64_t controller_state{};
  std::vector<RestartImageField> fields;
  std::array<std::vector<double>, 3U> final_mass_flux;
  bool backward_euler_recovery{true};

  void clear() noexcept;
};

class RestartWriter {
 public:
  static Status write(MPI_Comm communicator,
                      const std::filesystem::path& restart_directory,
                      const RestartSnapshot& snapshot,
                      RestartWriteOptions options = {}) noexcept;
};

class RestartReader {
 public:
  static Status load(MPI_Comm communicator,
                     const std::filesystem::path& restart_directory,
                     const RestartExpected& expected,
                     RestartImage& out) noexcept;
};

struct SnapshotFieldView {
  std::string_view stable_name;
  ConstFieldView values{};
  RevisionToken accepted_revision{};
};

struct CommittedOutputSnapshot {
  const CartesianGeometryPlan* geometry{};
  MeshPatch patch{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  double time{};
  std::uint64_t step{};
  Span<const SnapshotFieldView> fields{};
  bool committed{};
};

struct StageTimingRecord {
  StageId stage{};
  std::uint64_t minimum_nanoseconds{};
  std::uint64_t mean_nanoseconds{};
  std::uint64_t maximum_nanoseconds{};
};

// One pressure-solve contract applies to both correctors in a committed
// product attempt.  Invalid is the default so a runtime evidence row cannot
// infer linear-audit obligations from termination alone.
enum class RuntimePressureSolveContract : std::uint8_t {
  invalid,
  pressure_continuity,
  continuity_energy_coupled,
};

// Runtime evidence deliberately projects only the rank-invariant part of a
// same-target pressure--enthalpy refinement solve.  The rank-local pressure
// state and linear-system identities remain internal solver provenance.
inline constexpr std::size_t kRuntimePressureEnergyRefinementCapacity = 6U;

enum class RuntimePressureEnergyRefinementTermination : std::uint8_t {
  none,
  component_residuals_converged,
  iteration_capacity_exhausted,
  rejected_candidate,
};

struct RuntimePressureEnergyRefinementSolve {
  LinearSolveResult solve{};
  RevisionToken target_generation{};
  PlanFingerprint collective_lineage{};
  std::uint8_t ordinal{};
};

struct RuntimeTerminalPhysicalAudit {
  bool present{};
  RevisionToken final_flux_revision{};
  double eos_residual{};
  double eos_tolerance{};
  double continuity_residual{};
  double continuity_tolerance{};
  double energy_residual{};
  double energy_tolerance{};
  double closed_mass_residual{};
  double closed_mass_tolerance{};
  double gauge_residual{};
  double gauge_tolerance{};
};

struct RuntimeEvidenceRecord {
  PlanFingerprint build{};
  PlanFingerprint binary{};
  PlanFingerprint case_model{};
  PlanFingerprint stl{};
  PlanFingerprint product{};
  PlanFingerprint cpu_plan{};
  std::uint64_t step{};
  double time{};
  std::uint8_t requested_bdf_order{};
  std::uint8_t bdf_order{};
  std::uint8_t thermophysical_predictor_calls{};
  bool temporal_method_fallback{};
  std::uint64_t launcher_nanoseconds{};
  std::uint64_t maximum_rank_step_nanoseconds{};
  std::uint64_t maximum_rank_rss_bytes{};
  std::uint64_t maximum_node_rss_bytes{};
  std::uint64_t structured_messages{};
  std::uint64_t structured_bytes{};
  std::uint64_t ibm_messages{};
  std::uint64_t ibm_bytes{};
  std::uint64_t blocking_collectives{};
  std::uint64_t nonblocking_collectives{};
  std::uint64_t reduction_nanoseconds{};
  std::uint64_t linear_iterations{};
  std::uint64_t exact_numeric_refills{};
  std::uint64_t coarse_numeric_refills{};
  std::uint64_t preconditioner_setups{};
  std::uint64_t preconditioner_reuses{};
  std::uint64_t heap_allocations{};
  std::array<LinearSolveResult, 2U> pressure{};
  std::uint8_t pressure_solve_calls{};
  RuntimePressureSolveContract pressure_solve_contract{
      RuntimePressureSolveContract::invalid};
  std::array<RuntimePressureEnergyRefinementSolve,
             kRuntimePressureEnergyRefinementCapacity>
      pressure_energy_refinement{};
  std::uint8_t pressure_energy_refinement_solve_calls{};
  RuntimePressureEnergyRefinementTermination
      pressure_energy_refinement_termination{
          RuntimePressureEnergyRefinementTermination::none};
  RuntimeTerminalPhysicalAudit terminal_physical_audit{};
  std::array<LinearSolveResult, 3U> momentum_predictor{};
  std::uint8_t momentum_predictor_solve_calls{};
  LinearSolveResult predictor_enthalpy_endpoint{};
  double predictor_enthalpy_endpoint_alpha{1.0};
  double predictor_bdf_endpoint_alpha{1.0};
  double predictor_source_endpoint_alpha{1.0};
  std::uint8_t predictor_enthalpy_solve_calls{};
  double momentum_predictor_theta{1.0};
  std::uint32_t momentum_predictor_activations{};
  bool momentum_predictor_limited{};
  double predictor_theta{1.0};
  double predictor_mass_flux_scale{1.0};
  double predictor_low_margin{};
  double predictor_high_margin{};
  std::uint64_t predictor_low_order_transport_passes{};
  std::uint32_t predictor_low_order_substeps{};
  std::uint32_t predictor_low_order_halo_exchanges{};
  std::uint32_t predictor_blocking_collectives{};
  std::int32_t predictor_limiting_cell_x{};
  std::int32_t predictor_limiting_cell_y{};
  std::int32_t predictor_limiting_cell_z{};
  std::int32_t predictor_limiting_rank{-1};
  std::uint8_t predictor_constraint{};
  std::uint8_t predictor_low_state{};
  bool predictor_limited{};
  Span<const StageTimingRecord> stages{};
  bool startup{};
  bool retry{};
  bool restart_recovery{};
  bool statistics_eligible{};
};

class VisitWriter {
 public:
  static Status write(MPI_Comm communicator,
                      const std::filesystem::path& visit_directory,
                      const IoServicePlan& services,
                      const CommittedOutputSnapshot& snapshot) noexcept;
};

class ScreenWriter {
 public:
  static Status append(MPI_Comm communicator,
                       const std::filesystem::path& screen_file,
                       const IoServicePlan& services,
                       const CommittedOutputSnapshot& snapshot,
                       std::string_view summary) noexcept;
};

class MonitorWriter {
 public:
  static Status append(MPI_Comm communicator,
                       const std::filesystem::path& monitor_file,
                       const IoServicePlan& services,
                       const CommittedOutputSnapshot& snapshot,
                       std::string_view json_payload) noexcept;
};

class EvidenceWriter {
 public:
  static Status append(MPI_Comm communicator,
                       const std::filesystem::path& evidence_file,
                       const IoServicePlan& services,
                       const RuntimeEvidenceRecord& record) noexcept;
};

}  // namespace hundun::v04
