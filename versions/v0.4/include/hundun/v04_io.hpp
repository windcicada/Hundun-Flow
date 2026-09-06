// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

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

inline constexpr std::size_t kRuntimeSha256HexCharacters = 64U;
inline constexpr std::size_t kRuntimeGitObjectHexCharacters = 40U;
using RuntimeSha256Digest =
    std::array<char, kRuntimeSha256HexCharacters + 1U>;

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
  transported_scalar,
  enthalpy_nonadvective_rate,
  scalar_nonadvective_rate
};

struct RestartFieldView {
  RestartFieldRole role{RestartFieldRole::velocity};
  ConstFieldView values{};
};

// Borrowed synchronous snapshot: neither the metadata spans nor field/flux
// storage are owned here. See ProductDriver::committed_*_snapshot for lifetime.
// Copying this structure does NOT freeze the accepted state for asynchronous I/O.
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
  // Version-two Restart images carry the complete accepted BDF history.
  // The original fields/flux above remain the t_n authority so legacy
  // aggregate initializers and version-one readers retain their meaning.
  Span<const RestartFieldView> previous_fields{};
  Span<const RestartFieldView> accepted_rate_fields{};
  Span<const RestartFieldView> previous_rate_fields{};
  ConstFaceFluxView previous_mass_flux{};
  double previous_pressure_reference{};
  double closed_mass_target{};
  // Zero writes legacy unsigned V2 history. A nonzero semantic signature
  // writes V3 (same physical history plus this integrity-protected identity).
  PlanFingerprint method_history_signature{};
};

enum class IoFailureOperation : std::uint8_t {
  none, open, write, sync, close, create_directory,
  read, stat, rename, remove
};

struct IoFailureContext {
  bool valid{};
  IoFailureOperation operation{IoFailureOperation::none};
  int system_error{};
  int rank{-1};
  bool path_truncated{};
  std::array<char, 512U> path{};
};

enum class RestartPublicationState : std::uint8_t {
  not_switched,
  visible_not_durable,
  durable
};

struct RestartWriteReport {
  RestartPublicationState publication{RestartPublicationState::not_switched};
  IoFailureContext failure{};
  // Cleanup runs only after durable publication. Failure is a warning, not
  // evidence that current is unchanged; callers may continue using the new generation.
  Status cleanup_status{};
  IoFailureContext cleanup_failure{};
  std::size_t rank_payload_bytes{};
  // Peak bulk byte buffers, gather/record arrays; excludes borrowed fields,
  // small path strings, allocator bookkeeping, MPI and the solver resident set.
  std::size_t peak_bulk_staging_bytes{};
};

struct RestartWriteOptions {
  std::uint32_t keep_last{1U};
  // Rank-local consumer choice; does not alter the collective sequence.
  RestartWriteReport* report{};
  // Zero preserves the legacy caller's unbudgeted contract. Production callers
  // pass the sealed restart service capacity. Checked before bulk allocation.
  std::size_t maximum_bulk_staging_bytes{};
};

struct RestartReadReport {
  IoFailureContext failure{};
  std::uint32_t integrity_blocks{};
  std::uint32_t restoration_blocks{};
  std::uint64_t rank_file_bytes_read{};
};

struct RestartExpectedField {
  RestartFieldRole role{RestartFieldRole::velocity};
  FieldId field{};
  std::uint8_t components{};
};

// Opt-in compatibility for the pre-0b0eff raw MG bundle's outer ghost width.
// This never relaxes physical, geometry, field, or payload-integrity checks.
enum class RestartStorageCompatibility : std::uint8_t {
  strict,
  mg_bundle_ghost_v1
};

// A caller policy, not a statement about what the source file contains.
enum class RestartHistoryPolicy : std::uint8_t {
  require_compatible,
  rebuild_method_history
};

enum class RestartHistoryCompatibility : std::uint8_t {
  missing, unknown, compatible, incompatible
};

struct RestartExpected {
  Int3 global_cells{};
  MeshPatch target_patch{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  PlanFingerprint geometry{};
  Span<const RestartExpectedField> fields{};
  Span<const RestartExpectedField> rate_fields{};
  PlanFingerprint compatible_storage_plan{};
  PlanFingerprint compatible_storage_schema{};
  PlanFingerprint method_history_signature{};
  // Read-only, known historical method identity. Populated only for explicit
  // method recovery; geometry, physical configuration and storage still match.
  PlanFingerprint compatible_method_plan{};
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
  RuntimeSha256Digest source_manifest_sha256{};
  std::vector<RestartImageField> fields;
  std::array<std::vector<double>, 3U> final_mass_flux;
  // True only for a legacy version-one image whose absent t_{n-1} state is
  // synthesized from t_n.  Callers must insert one backward-Euler recovery
  // step in that case. False denotes complete V2/V3 source history, not proof
  // of same-method compatibility; use history_compatibility and a separate
  // RestartHistoryPolicy. Callers must not repurpose this file-format fact.
  bool backward_euler_recovery{true};
  std::vector<RestartImageField> previous_fields;
  std::vector<RestartImageField> accepted_rate_fields;
  std::vector<RestartImageField> previous_rate_fields;
  std::array<std::vector<double>, 3U> previous_mass_flux;
  double previous_pressure_reference{};
  double closed_mass_target{};
  RevisionToken final_mass_flux_revision{};
  RevisionToken previous_mass_flux_revision{};
  // Source plan/schema above remain unchanged when this is true.
  bool storage_layout_migrated{};
  // File facts remain unchanged when a caller requests method recovery.
  std::uint32_t source_format_version{1U};
  PlanFingerprint method_history_signature{};

  RestartHistoryCompatibility history_compatibility(
      PlanFingerprint expected_signature) const noexcept {
    if (backward_euler_recovery) return RestartHistoryCompatibility::missing;
    if (method_history_signature == 0U || expected_signature == 0U)
      return RestartHistoryCompatibility::unknown;
    return method_history_signature == expected_signature
        ? RestartHistoryCompatibility::compatible
        : RestartHistoryCompatibility::incompatible;
  }

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
                     RestartImage& out,
                     RestartReadReport* report = nullptr) noexcept;
};

struct SnapshotFieldView {
  std::string_view stable_name;
  ConstFieldView values{};
  RevisionToken accepted_revision{};
};

// Borrowed metadata and fields, not an owning image or a checked epoch lease.
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

// One terminal pressure--energy acceptance contract applies to a committed
// product attempt. Internal linear directions may differ by coupling; invalid
// is the default so evidence cannot infer obligations from termination alone.
enum class RuntimePressureSolveContract : std::uint8_t {
  invalid,
  pressure_continuity,
  continuity_energy_coupled,
};

enum class RuntimeCouplingKind : std::uint8_t {
  invalid,
  piso,
  simple,
};

// Runtime evidence deliberately projects only the rank-invariant part of a
// same-target pressure--enthalpy refinement solve.  The rank-local pressure
// state and linear-system identities remain internal solver provenance.
inline constexpr std::size_t kRuntimePressureEnergyRefinementCapacity = 12U;

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

struct RuntimeConvectiveCflWinner {
  bool valid{};
  Int3 global_cell{};
  std::int32_t rank{-1};
  double out{};
  double absolute{};
  double density_volume{};
  double outgoing_mass_flow{};
  double absolute_mass_flow{};
};

struct RuntimeCommittedConvectiveCflAudit {
  bool valid{};
  RevisionToken density_revision{};
  RevisionToken final_flux_revision{};
  FieldId density_field{};
  // These communicator-wide fingerprints bind every rank-local logical view
  // identity (storage/revision domain and layout) without publishing a
  // rank-local token in an evidence row that must be byte-identical on all
  // ranks.
  PlanFingerprint density_view_collective{};
  PlanFingerprint final_flux_view_collective{};
  PlanFingerprint activity_collective{};
  double dt{};
  double out_max{};
  double abs_max{};
  double limit{};
  RuntimeConvectiveCflWinner out_winner{};
  RuntimeConvectiveCflWinner abs_winner{};
};

struct RuntimeAdvectiveCflAudit {
  bool present{};
  PlanFingerprint plan{};
  PlanFingerprint time_revision_collective{};
  PlanFingerprint density_view_collective{};
  PlanFingerprint face_flux_view_collective{};
  PlanFingerprint activity_collective{};
  double dt{};
  double out_max{};
  double abs_max{};
  double limit{};
};

// Runtime candidate identity contains only immutable content digests.  Paths,
// mtimes and configure timestamps are deliberately excluded.  `head` and
// `tree` are the exact Git object ids embedded by CMake; `executable` is
// computed from the bytes of the image that is actually running.  `identity`
// binds the four preceding values together with the runtime-evidence schema
// under the versioned canonical payload used by both the writer and the
// external validator.  This prevents a frozen row from being relabelled as a
// newer evidence policy without changing its immutable candidate identity.
struct RuntimeCandidateIdentity {
  std::array<char, kRuntimeGitObjectHexCharacters + 1U> head{};
  std::array<char, kRuntimeGitObjectHexCharacters + 1U> tree{};
  std::array<char, kRuntimeSha256HexCharacters + 1U> build_manifest{};
  std::array<char, kRuntimeSha256HexCharacters + 1U> executable{};
  std::array<char, kRuntimeSha256HexCharacters + 1U> identity{};
};

enum class RuntimeRunStartKind : std::uint8_t { fresh, restart };

// The first evidence row is anchored independently of that row's proposed
// time.  Restart runs bind the prior step/time to the exact integrity-checked
// manifest loaded by RestartReader; fresh runs bind the canonical zero state.
struct RuntimeRunStartAnchor {
  RuntimeRunStartKind kind{RuntimeRunStartKind::fresh};
  std::uint64_t previous_step{};
  double previous_time{};
  RuntimeSha256Digest restart_manifest_sha256{};
  // Optional for historical evidence; new restart producers fill all fields.
  std::uint32_t source_format_version{};
  PlanFingerprint source_history_signature{};
  PlanFingerprint target_history_signature{};
  RestartHistoryPolicy history_policy{RestartHistoryPolicy::require_compatible};
};

struct RuntimeEvidenceRecord {
  PlanFingerprint build{};
  PlanFingerprint binary{};
  RuntimeCandidateIdentity candidate_identity{};
  PlanFingerprint case_model{};
  PlanFingerprint stl{};
  PlanFingerprint product{};
  PlanFingerprint cpu_plan{};
  RuntimeRunStartAnchor run_start{};
  std::uint64_t step{};
  double previous_committed_time{};
  double time{};
  std::uint8_t requested_bdf_order{};
  std::uint8_t bdf_order{};
  RuntimeCouplingKind coupling{RuntimeCouplingKind::invalid};
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
  RuntimeCommittedConvectiveCflAudit committed_convective_cfl{};
  std::array<LinearSolveResult, 3U> momentum_predictor{};
  std::uint8_t momentum_predictor_solve_calls{};
  std::uint8_t momentum_predictor_passes{};
  LinearSolveResult predictor_enthalpy_endpoint{};
  double predictor_enthalpy_endpoint_alpha{1.0};
  double predictor_bdf_endpoint_alpha{1.0};
  double predictor_source_endpoint_alpha{1.0};
  std::uint8_t predictor_enthalpy_solve_calls{};
  double momentum_predictor_theta{1.0};
  std::uint32_t momentum_predictor_activations{};
  bool momentum_correction_metrics_applicable{};
  double momentum_minimum_face_alpha{};
  std::uint32_t momentum_active_correction_faces{};
  double momentum_limited_face_fraction{};
  bool momentum_predictor_limited{};
  RuntimeAdvectiveCflAudit momentum_advective_cfl{};
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
                      const CommittedOutputSnapshot& snapshot,
                      IoFailureContext* failure = nullptr) noexcept;
};

class ScreenWriter {
 public:
  static Status append(MPI_Comm communicator,
                       const std::filesystem::path& screen_file,
                       const IoServicePlan& services,
                       const CommittedOutputSnapshot& snapshot,
                       std::string_view summary,
                       IoFailureContext* failure = nullptr) noexcept;
};

class MonitorWriter {
 public:
  static Status append(MPI_Comm communicator,
                       const std::filesystem::path& monitor_file,
                       const IoServicePlan& services,
                       const CommittedOutputSnapshot& snapshot,
                       std::string_view json_payload,
                       IoFailureContext* failure = nullptr) noexcept;
};

class EvidenceWriter {
 public:
  static Status append(MPI_Comm communicator,
                       const std::filesystem::path& evidence_file,
                       const IoServicePlan& services,
                       const RuntimeEvidenceRecord& record,
                       IoFailureContext* failure = nullptr) noexcept;
};

}  // namespace hundun::v04
