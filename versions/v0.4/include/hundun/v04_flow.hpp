// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_linear.hpp"
#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04 {

class IbmEquationInterfacePlan;
class EBTopology;
class RemoteDonorExchangePlan;
class PressureEnergyCandidateBoundaryFinalizer;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
class FinalBoundaryFluxCertificateTestAccess;
#endif

enum class EquationAssemblyScope : std::uint8_t {
  momentum_predictor,
  target_coupled,
  final_conservative
};

enum class MomentumStressForm : std::uint8_t {
  full_newtonian,
  laplacian_only
};

struct ScalarEquationSpec {
  FieldId field{};
  TransportedScalarRole role{TransportedScalarRole::passive_scalar};
  double molecular_schmidt{1.0};
  double turbulent_schmidt{1.0};
};

struct EquationPlanSpec {
  FieldId density{};
  FieldId velocity{};
  FieldId pressure_perturbation{};
  FieldId enthalpy{};
  FieldId temperature{};
  FieldId effective_viscosity{};
  FieldId velocity_gradient{};
  FieldId pressure_compressibility{};
  PressureReferenceKind pressure_reference{
      PressureReferenceKind::boundary_absolute};
  Span<const ScalarEquationSpec> scalars{};
  std::size_t maximum_cells_per_rank{};
  StageId closed_mass_service_stage{};
};

struct EquationCompileDiagnostics {
  int lowest_failing_rank{-1};
};

struct PrimitiveHistory {
  ConstFieldView trial{};
  ConstFieldView accepted{};
  ConstFieldView previous{};
};

struct EquationStateView {
  PrimitiveHistory density{};
  PrimitiveHistory velocity{};
  PrimitiveHistory pressure_perturbation{};
  PrimitiveHistory enthalpy{};
  PrimitiveHistory temperature{};
  double pressure_reference{};
  Span<const PrimitiveHistory> independent_species{};
  Span<const PrimitiveHistory> passive_scalars{};
  double accepted_pressure_reference{};
  double previous_pressure_reference{};
};

struct EquationMaterialView {
  ConstFieldView molecular_viscosity{};
  ConstFieldView turbulent_viscosity{};
  ConstFieldView effective_viscosity{};
  ConstFieldView thermal_conductivity{};
  ConstFieldView heat_capacity{};
  // Frozen local dT/dh=1/cp linearization: lambda_over_cp=lambda/cp.
  ConstFieldView enthalpy_diffusivity{};
  // Species first, then passive scalars.  Values are rho*D [kg/(m s)].
  Span<const ConstFieldView> scalar_mass_diffusivity{};
  ConstFieldView pressure_compressibility{};
};

enum class PressureGauge : std::uint8_t {
  absolute_boundary_dirichlet,
  compressibility_weighted_zero_mean
};

struct EquationAssemblyContext;
struct EquationAssemblyCertificate;
struct ThermophysicalPredictorCertificate;
struct PressureReferenceCertificate;
struct ClosedGaugeCorrectionPrepareInput;
struct ClosedGaugeCorrectionCertificate;

class PressureReferencePlan {
 public:
  PressureReferenceKind kind() const noexcept { return kind_; }
  PressureGauge gauge() const noexcept { return gauge_; }
  StageId service_stage() const noexcept { return service_stage_; }
  FieldId compressibility_field() const noexcept {
    return compressibility_field_;
  }
  FieldId pressure_perturbation_field() const noexcept {
    return pressure_perturbation_field_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Int3 cells() const noexcept { return cells_; }
  Status solve_closed_mass(
      MPI_Comm communicator, const ClosedMassPlan& closed_mass,
      const ThermodynamicsPlan& thermodynamics,
      const ThermophysicalPredictorCertificate& predictor, StageId stage,
      const ClosedMassCellView& cells, double target_mass,
      double current_pressure_reference, ClosedMassResult& result,
      PressureReferenceCertificate& certificate) const noexcept;
  Status solve_closed_mass_fields(
      MPI_Comm communicator, const ClosedMassPlan& closed_mass,
      const ThermodynamicsPlan& thermodynamics,
      const ThermophysicalPredictorCertificate& predictor, StageId stage,
      const ClosedMassFieldView& cells, double target_mass,
      double current_pressure_reference, ClosedMassResult& result,
      PressureReferenceCertificate& certificate) const noexcept;
  Status certify_closed_mass_density_fields(
      MPI_Comm communicator, const ClosedMassPlan& closed_mass,
      const ThermodynamicsPlan& thermodynamics,
      const ThermophysicalPredictorCertificate& predictor, StageId stage,
      const ClosedMassDensityFieldView& cells, double target_mass,
      double current_pressure_reference, ClosedMassResult& result,
      PressureReferenceCertificate& certificate) const noexcept;
  Status normalize_closed_gauge(FieldView pressure_perturbation,
                                ConstFieldView drho_dp_h_y,
                                ReductionEngine& reductions,
                                double& pressure_reference) const noexcept;
  Status prepare_closed_gauge_correction(
      const ClosedGaugeCorrectionPrepareInput& input,
      ReductionEngine& reductions,
      ClosedGaugeCorrectionCertificate& certificate) const noexcept;
  // Allocation-free, rank-local revalidation for the write transaction.  A
  // caller publishes a false result through its own existing collective
  // preflight; this seam deliberately performs no communication.
  bool matches_closed_gauge_correction(
      const ClosedGaugeCorrectionPrepareInput& input,
      const ClosedGaugeCorrectionCertificate& certificate) const noexcept;
  Status certify_open_boundary(
      double pressure_reference, RevisionToken boundary_value_revision,
      const ThermophysicalPredictorCertificate& predictor,
      PressureReferenceCertificate& certificate) const noexcept;

 private:
  friend class EquationPlanSet;
  friend Status assemble_pressure_storage(
      const PressureReferencePlan&, ConstFieldView,
      const PressureReferenceCertificate&, const EquationAssemblyContext&,
      FieldView,
      EquationAssemblyCertificate&) noexcept;
  Int3 cells_{};
  const CartesianKernelPlan* kernels_{};
  PressureReferenceKind kind_{PressureReferenceKind::boundary_absolute};
  PressureGauge gauge_{PressureGauge::absolute_boundary_dirichlet};
  FieldId pressure_perturbation_field_{};
  FieldId compressibility_field_{};
  StageId service_stage_{};
  RevisionToken geometry_revision_{};
  PlanFingerprint thermodynamics_fingerprint_{};
  PlanFingerprint predictor_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

struct EquationContributionView {
  // Per-volume RHS source and non-negative implicit sink coefficient.
  ConstFieldView explicit_source_density{};
  ConstFieldView implicit_sink_density{};
  bool has_implicit_sink{};
  FieldId conserved_quantity{};
  UnitDimension units{};
  StageId stage{};
  FieldId explicit_source_field{};
  FieldId implicit_sink_field{};
};

struct EquationAssemblyContext {
  double dt{};
  BdfCoefficients bdf{};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken boundary{};
  RevisionToken thermo{};
  RevisionToken transport{};
  RevisionToken face_flux{};
  PlanFingerprint face_flux_authority{};
  StorageIdentity face_flux_storage{};
  RevisionDomainIdentity face_flux_revision_domain{};
  StageId contribution_stage{};
  EquationAssemblyScope scope{EquationAssemblyScope::momentum_predictor};
  ConstFaceFluxView mass_flux{};
  bool provisional_mass_flux{};
  KernelBox box{};
  KernelCounters* counters{};
  const IbmEquationInterfacePlan* immersed_interface{};
  const TurbulencePlan* wall_treatment{};
};

struct EquationSystemView {
  FieldView diagonal{};
  FieldView rhs{};
  FieldView residual{};
  FaceFieldView x_coefficient{};
  FaceFieldView y_coefficient{};
  FaceFieldView z_coefficient{};
};

struct MomentumPredictorLimiterReport {
  double theta{1.0};
  std::uint32_t activations{};
  bool limited{};
};

struct MomentumPredictorSolveReport {
  std::array<LinearSolveResult, 3U> components{};
  std::uint8_t solve_calls{};
};

struct EquationAssemblyCertificate {
  PlanFingerprint plan{};
  EquationAssemblyScope scope{EquationAssemblyScope::momentum_predictor};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken face_flux{};
  RevisionToken state{};

  bool valid() const noexcept {
    return plan != 0U && time != 0U && geometry != 0U && face_flux != 0U &&
           state != 0U;
  }
};

struct PredictorRateHistory {
  // Positive conservative RHS density [quantity/(m^3 s)].  An empty view
  // means an exactly zero registered non-advective rate.
  ConstFieldView accepted{};
  ConstFieldView previous{};
};

struct ThermophysicalGhostAuthority {
  std::uintptr_t exchange_plan{};
  FieldId field{};
  RevisionToken state{};
  StorageIdentity storage{};
  RevisionDomainIdentity revision_domain{};
  RevisionToken geometry{};
  RevisionToken boundary{};
  std::uint8_t reach{};

  bool valid() const noexcept {
    return exchange_plan != 0U && state != 0U && storage != 0U &&
           revision_domain != 0U && geometry != 0U && boundary != 0U &&
           reach != 0U;
  }
};

struct ThermophysicalGhostHistory {
  ThermophysicalGhostAuthority accepted{};
  ThermophysicalGhostAuthority previous{};
};

struct ThermophysicalPredictorInput {
  double dt{};
  BdfCoefficients bdf{};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken boundary{};
  RevisionToken transport{};
  ConstFieldView density_accepted{};
  ConstFieldView density_previous{};
  ConstFieldView enthalpy_accepted{};
  ConstFieldView enthalpy_previous{};
  Span<const ConstFieldView> species_accepted{};
  Span<const ConstFieldView> species_previous{};
  Span<const ConstFieldView> passive_scalars_accepted{};
  Span<const ConstFieldView> passive_scalars_previous{};
  ConstFaceFluxView mass_flux_accepted{};
  ConstFaceFluxView mass_flux_previous{};
  PredictorRateHistory enthalpy_nonadvective_rhs{};
  Span<const PredictorRateHistory> species_nonadvective_rhs{};
  Span<const PredictorRateHistory> passive_scalar_nonadvective_rhs{};
  ThermophysicalGhostHistory enthalpy_ghosts{};
  Span<const ThermophysicalGhostHistory> species_ghosts{};
  Span<const ThermophysicalGhostHistory> passive_scalar_ghosts{};
  KernelCounters* counters{};
};

struct ThermophysicalPredictorOutput {
  FieldView enthalpy{};
  Span<FieldView> independent_species{};
  // Attempt-local rho* is formed once and reused by h and every Y component.
  FieldView density_workspace{};
  // Reused for each transported quantity; both are attempt-local workspace.
  FieldView accepted_advection_workspace{};
  FieldView previous_advection_workspace{};
  Span<FieldView> passive_scalars{};
  // Cold-allocated cell bundle used only when the high BDF/EX candidate is
  // inadmissible.  These views never alias the published high output.
  FieldView low_order_density_workspace{};
  FieldView low_order_enthalpy_workspace{};
  Span<FieldView> low_order_independent_species{};
  Span<FieldView> low_order_passive_scalars{};
  // Attempt-local mass-flux authority paired with density_workspace.  The
  // caller owns the preallocated face replica; the predictor only fills it
  // and binds its storage/revision into the returned certificate.
  FaceFluxView paired_mass_flux{};
};

struct ConservativeEnthalpyEndpointInput {
  BdfCoefficients bdf{};
  RevisionToken time{};
  ConstFieldView predicted_density{};
  ConstFieldView accepted_enthalpy{};
  ConstFieldView high_enthalpy{};
  ConstFaceFluxView accepted_mass_flux{};
  FieldView diagonal_workspace{};
  FieldView rhs_workspace{};
  FieldView endpoint{};
  ResourceCounters* resources{};
};

struct ConservativeEnthalpyEndpointServices {
  MPI_Comm communicator{MPI_COMM_NULL};
  const CartesianKernelPlan* kernels{};
  const BoundaryPlan* boundary{};
  MeshPatch patch{};
  MgDomainActivityView activity{};
  HaloEngine* halo{};
  ReductionEngine* reductions{};
  StageId halo_stage{};
  LinearWorkspaceRequirements workspace_requirements{};
  FieldView krylov_vectors{};
  FieldView krylov_scalars{};
  FieldId enthalpy{};
};

class ConservativeEnthalpyEndpoint {
 public:
  ConservativeEnthalpyEndpoint() noexcept = default;
  ConservativeEnthalpyEndpoint(const ConservativeEnthalpyEndpoint&) = delete;
  ConservativeEnthalpyEndpoint& operator=(
      const ConservativeEnthalpyEndpoint&) = delete;
  ConservativeEnthalpyEndpoint(ConservativeEnthalpyEndpoint&&) = delete;
  ConservativeEnthalpyEndpoint& operator=(
      ConservativeEnthalpyEndpoint&&) = delete;

  static Status bind(const ConservativeEnthalpyEndpointServices& services,
                     ConservativeEnthalpyEndpoint& out) noexcept;
  Status solve(const ConservativeEnthalpyEndpointInput& input,
               LinearSolveResult& result) noexcept;
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  MPI_Comm communicator_{MPI_COMM_NULL};
  const CartesianKernelPlan* kernels_{};
  const BoundaryPlan* boundary_{};
  MeshPatch patch_{};
  MgDomainActivityView activity_{};
  HaloEngine* halo_{};
  ReductionEngine* reductions_{};
  StageId halo_stage_{};
  SolverWorkspace workspace_{};
  FieldId enthalpy_{};
  std::array<BoundaryRelation, 6U> boundary_relations_{};
  std::array<bool, 6U> physical_boundaries_{};
  PlanFingerprint collective_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

struct ThermophysicalPredictorSlowPath {
  // Dedicated one-scalar donor-factor exchange.  It is used only when the
  // accepted-flux BDF endpoint requires local outgoing-flux limiting.
  HaloEngine* halo{};
  StageId halo_stage{};
  // h, then independent species, then passive scalars.  The views alias the
  // corresponding output fields and are preallocated by the product runtime.
  Span<FieldView> exchange_fields{};
  BoundaryResolvedValues boundary_values{};
  ConservativeEnthalpyEndpoint* enthalpy_endpoint{};
  ResourceCounters* resources{};
  const IbmEquationInterfacePlan* immersed_interface{};
};

enum class ThermophysicalAdmissibilityConstraint : std::uint8_t {
  none,
  density,
  independent_species,
  dependent_species,
  enthalpy_lower,
  enthalpy_upper
};

enum class ThermophysicalLowStateKind : std::uint8_t {
  none = 0U,
  bdf_accepted_rate = 1U,
  scaled_euler = 2U,
  implicit_upwind = 3U,
  implicit_upwind_source_limited = 4U,
  bdf_accepted_rate_homotopy = 5U,
  bdf_local_donor_flux = 6U,
  bdf_local_donor_flux_source_limited = 7U
};

enum class ThermophysicalPredictorFailureReason : std::uint8_t {
  none,
  high_input_preflight,
  high_density_bdf,
  high_quantity_transport,
  high_quantity_bdf,
  low_density_endpoint,
  low_cfl,
  substep_arithmetic,
  transport_pass_overflow,
  low_divergence,
  low_conserved_update,
  low_tuple_admissibility,
  implicit_enthalpy_endpoint,
  limiter_theta,
  limiter_metadata,
  blended_density,
  final_tuple_admissibility,
  low_bdf_source_base_admissibility
};

enum class ThermophysicalPredictorFailureField : std::uint8_t {
  none,
  density,
  enthalpy,
  mass_flux,
  independent_species,
  passive_scalar,
  nonadvective_rhs
};

struct ThermophysicalPredictorFailure {
  // This record is caller-owned fixed storage.  The predictor only publishes
  // it into caller diagnostics on a numerical failure; successful calls use
  // the existing diagnostics reset path.
  bool valid{};
  ThermophysicalPredictorFailureReason reason{
      ThermophysicalPredictorFailureReason::none};
  ThermophysicalPredictorFailureField field{
      ThermophysicalPredictorFailureField::none};
  ThermophysicalAdmissibilityConstraint constraint{
      ThermophysicalAdmissibilityConstraint::none};
  std::uint32_t field_index{UINT32_MAX};
  std::int32_t rank{-1};
  Int3 global_index{-1, -1, -1};
  std::uint32_t substep{};
  bool has_cell{};
  bool has_substep{};
  bool has_margins{};
  double low_margin{};
  double high_margin{};

  // Presence bits make a finite zero distinguishable from an unavailable
  // operand.  The payload is intentionally fixed-size and scalar-only so it
  // can be sent through MPI without transmitting a padded C++ object.
  std::uint32_t scalar_mask{};
  static constexpr std::uint32_t scalar_dt_substep = 1U << 0U;
  static constexpr std::uint32_t scalar_maximum_cfl = 1U << 1U;
  static constexpr std::uint32_t scalar_local_cfl = 1U << 2U;
  static constexpr std::uint32_t scalar_density_current = 1U << 3U;
  static constexpr std::uint32_t scalar_density_next = 1U << 4U;
  static constexpr std::uint32_t scalar_quantity_current = 1U << 5U;
  static constexpr std::uint32_t scalar_quantity_previous = 1U << 6U;
  static constexpr std::uint32_t scalar_divergence = 1U << 7U;
  static constexpr std::uint32_t scalar_nonadvective_rhs = 1U << 8U;
  static constexpr std::uint32_t scalar_conserved_next = 1U << 9U;
  static constexpr std::uint32_t scalar_observed_value = 1U << 10U;
  static constexpr std::uint32_t scalar_allowed_lower = 1U << 11U;
  static constexpr std::uint32_t scalar_allowed_upper = 1U << 12U;

  double dt_substep{};
  double maximum_cfl{};
  double local_cfl{};
  double density_current{};
  double density_next{};
  double quantity_current{};
  double quantity_previous{};
  double divergence{};
  double nonadvective_rhs{};
  double conserved_next{};
  double observed_value{};
  double allowed_lower{};
  double allowed_upper{};
};

struct ThermophysicalPredictorDiagnostics {
  double theta{1.0};
  double mass_flux_scale{1.0};
  double low_margin{};
  double high_margin{};
  std::uint32_t low_order_substeps{};
  std::uint32_t low_order_halo_exchanges{};
  std::uint32_t blocking_collectives{};
  std::uint64_t low_order_transport_passes{};
  // Global mesh index of the selected limiting cell.
  Int3 limiting_cell{};
  std::int32_t limiting_rank{-1};
  ThermophysicalAdmissibilityConstraint constraint{
      ThermophysicalAdmissibilityConstraint::none};
  ThermophysicalLowStateKind low_state{ThermophysicalLowStateKind::none};
  LinearSolveResult enthalpy_endpoint{};
  double enthalpy_endpoint_alpha{1.0};
  double bdf_endpoint_alpha{1.0};
  double source_endpoint_alpha{1.0};
  std::uint8_t enthalpy_solve_calls{};
  bool limited{};
  ThermophysicalPredictorFailure failure{};
};

struct ThermophysicalPredictorCertificate {
  PlanFingerprint plan{};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken accepted_face_flux{};
  RevisionToken previous_face_flux{};
  PlanFingerprint committed_face_flux_authority{};
  StorageIdentity committed_face_flux_storage{};
  RevisionDomainIdentity committed_face_flux_revision_domain{};
  RevisionToken predicted_density{};
  StorageIdentity predicted_density_storage{};
  RevisionDomainIdentity predicted_density_revision_domain{};
  RevisionToken paired_face_flux{};
  StorageIdentity paired_face_flux_storage{};
  RevisionDomainIdentity paired_face_flux_revision_domain{};
  RevisionToken state{};
  std::uint8_t order{};

  bool valid() const noexcept {
    return plan != 0U && time != 0U && geometry != 0U &&
           accepted_face_flux != 0U && predicted_density != 0U &&
           committed_face_flux_authority != 0U &&
           committed_face_flux_storage != 0U &&
           committed_face_flux_revision_domain != 0U &&
           predicted_density_storage != 0U &&
           predicted_density_revision_domain != 0U &&
           paired_face_flux != 0U && paired_face_flux_storage != 0U &&
           paired_face_flux_revision_domain != 0U && state != 0U &&
           ((order == 1U && previous_face_flux == 0U) ||
            (order == 2U && previous_face_flux != 0U));
  }
};

struct ThermophysicalRateInput {
  EquationStateView state{};
  EquationMaterialView material{};
  ConstFieldView velocity_gradient{};
  BdfCoefficients bdf{};
  StageId contribution_stage{};
  Span<const EquationContributionView> contributions{};
  const IbmEquationInterfacePlan* immersed_interface{};
};

struct ThermophysicalRateOutput {
  FieldView enthalpy_nonadvective_rhs{};
  Span<FieldView> species_nonadvective_rhs{};
  Span<FieldView> passive_scalar_nonadvective_rhs{};
  FieldView diffusion_scratch{};
  FieldView scalar_diffusivity_workspace{};
};

struct ThermophysicalRateCertificate {
  PlanFingerprint equations{};
  RevisionToken state{};
  RevisionToken rates{};

  bool valid() const noexcept {
    return equations != 0U && state != 0U && rates != 0U;
  }
};

struct PressureReferenceCertificate {
  PlanFingerprint plan{};
  PlanFingerprint predictor{};
  PlanFingerprint thermodynamics{};
  PlanFingerprint closure{};
  RevisionToken time{};
  RevisionToken pressure_reference{};
  PressureReferenceKind kind{PressureReferenceKind::boundary_absolute};

  bool valid() const noexcept {
    return plan != 0U && predictor != 0U && thermodynamics != 0U &&
           closure != 0U && time != 0U && pressure_reference != 0U;
  }
};

class ThermophysicalPredictorPlan {
 public:
  ThermophysicalPredictorPlan() noexcept = default;
  ThermophysicalPredictorPlan(const ThermophysicalPredictorPlan&) = delete;
  ThermophysicalPredictorPlan& operator=(const ThermophysicalPredictorPlan&) =
      delete;
  ThermophysicalPredictorPlan(ThermophysicalPredictorPlan&&) = delete;
  ThermophysicalPredictorPlan& operator=(ThermophysicalPredictorPlan&&) =
      delete;
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  std::size_t species_count() const noexcept { return species_.size(); }
  std::size_t passive_scalar_count() const noexcept {
    return passive_scalars_.size();
  }
  ConvectionScheme enthalpy_convection() const noexcept {
    return enthalpy_convection_;
  }
  std::uint8_t enthalpy_reach() const noexcept { return enthalpy_reach_; }
  std::uint8_t species_reach() const noexcept { return species_reach_; }
  std::uint8_t passive_scalar_reach() const noexcept {
    return passive_scalar_reach_;
  }
  Status predict(MPI_Comm communicator,
                 Status prerequisite,
                 const ThermophysicalPredictorInput& input,
                 ThermophysicalPredictorOutput output,
                 const ThermophysicalPredictorSlowPath& slow_path,
                 ThermophysicalPredictorDiagnostics& diagnostics,
                 ThermophysicalPredictorCertificate& certificate) const
      noexcept;

 private:
  friend class EquationPlanSet;
  Status predict_high_local(
      const ThermophysicalPredictorInput& input,
      ThermophysicalPredictorOutput output,
      ThermophysicalPredictorCertificate& certificate,
      ThermophysicalPredictorFailure& failure) const noexcept;
  const CartesianKernelPlan* kernels_{};
  Int3 cells_{};
  Int3 patch_begin_{};
  FieldId density_{};
  FieldId enthalpy_{};
  std::vector<FieldId> species_;
  std::vector<FieldId> passive_scalars_;
  std::vector<double> species_enthalpy_minimum_;
  std::vector<double> species_enthalpy_maximum_;
  double dependent_enthalpy_minimum_{};
  double dependent_enthalpy_maximum_{};
  ConvectionScheme enthalpy_convection_{ConvectionScheme::limited_central2};
  ConvectionScheme species_convection_{ConvectionScheme::tvd2};
  ConvectionScheme passive_scalar_convection_{ConvectionScheme::tvd2};
  std::uint8_t enthalpy_reach_{2U};
  std::uint8_t species_reach_{2U};
  std::uint8_t passive_scalar_reach_{2U};
  const BoundaryPlan* boundary_{};
  RevisionToken geometry_revision_{};
  RevisionToken boundary_revision_{};
  PlanFingerprint transport_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

class EquationPlanSet;

class ContinuityEquationPlan {
 public:
  ContinuityEquationPlan() noexcept = default;
  ContinuityEquationPlan(const ContinuityEquationPlan&) = delete;
  ContinuityEquationPlan& operator=(const ContinuityEquationPlan&) = delete;
  ContinuityEquationPlan(ContinuityEquationPlan&&) = delete;
  ContinuityEquationPlan& operator=(ContinuityEquationPlan&&) = delete;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Int3 cells() const noexcept { return cells_; }
  PressureReferenceKind pressure_reference() const noexcept {
    return pressure_reference_;
  }

 private:
  friend class EquationPlanSet;
  friend Status assemble_continuity(
      const ContinuityEquationPlan&, const EquationStateView&,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&) noexcept;
  friend Status assemble_continuity_impl(
      const ContinuityEquationPlan&, const EquationStateView&,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&, bool) noexcept;
  const CartesianKernelPlan* kernels_{};
  Int3 cells_{};
  FieldId density_{};
  FieldId velocity_{};
  PressureReferenceKind pressure_reference_{
      PressureReferenceKind::boundary_absolute};
  RevisionToken geometry_revision_{};
  PlanFingerprint fingerprint_{};
};

class MomentumEquationPlan {
 public:
  MomentumEquationPlan() noexcept = default;
  MomentumEquationPlan(const MomentumEquationPlan&) = delete;
  MomentumEquationPlan& operator=(const MomentumEquationPlan&) = delete;
  MomentumEquationPlan(MomentumEquationPlan&&) = delete;
  MomentumEquationPlan& operator=(MomentumEquationPlan&&) = delete;

  MomentumStressForm stress_form() const noexcept { return stress_form_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Int3 cells() const noexcept { return cells_; }

 private:
  friend class EquationPlanSet;
  friend Status assemble_momentum(
      const MomentumEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, const EquationAssemblyContext&,
      EquationSystemView, EquationAssemblyCertificate&) noexcept;
  friend Status assemble_momentum_impl(
      const MomentumEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, const EquationAssemblyContext&,
      EquationSystemView, FieldView, EquationAssemblyCertificate&,
      bool) noexcept;
  const CartesianKernelPlan* kernels_{};
  Int3 cells_{};
  FieldId density_{};
  FieldId velocity_{};
  FieldId pressure_{};
  FieldId viscosity_{};
  FieldId velocity_gradient_{};
  ConvectionScheme convection_{ConvectionScheme::limited_central2};
  std::uint8_t convection_reach_{2U};
  MomentumStressForm stress_form_{MomentumStressForm::full_newtonian};
  std::vector<CompiledContribution> contributions_;
  RevisionToken geometry_revision_{};
  RevisionToken boundary_revision_{};
  PlanFingerprint transport_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

class EnthalpyEquationPlan {
 public:
  EnthalpyEquationPlan() noexcept = default;
  EnthalpyEquationPlan(const EnthalpyEquationPlan&) = delete;
  EnthalpyEquationPlan& operator=(const EnthalpyEquationPlan&) = delete;
  EnthalpyEquationPlan(EnthalpyEquationPlan&&) = delete;
  EnthalpyEquationPlan& operator=(EnthalpyEquationPlan&&) = delete;

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Int3 cells() const noexcept { return cells_; }

 private:
  friend class EquationPlanSet;
  friend Status assemble_enthalpy(
      const EnthalpyEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, const EquationAssemblyContext&,
      EquationSystemView, EquationAssemblyCertificate&) noexcept;
  friend Status assemble_enthalpy_impl(
      const EnthalpyEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, const EquationAssemblyContext&,
      EquationSystemView, EquationAssemblyCertificate&, bool) noexcept;
  friend Status evaluate_thermophysical_rates(
      const EquationPlanSet&, const ThermophysicalRateInput&,
      ThermophysicalRateOutput, ThermophysicalRateCertificate&) noexcept;
  const CartesianKernelPlan* kernels_{};
  Int3 cells_{};
  FieldId density_{};
  FieldId velocity_{};
  FieldId pressure_{};
  FieldId enthalpy_{};
  FieldId temperature_{};
  FieldId velocity_gradient_{};
  ConvectionScheme convection_{ConvectionScheme::limited_central2};
  std::uint8_t convection_reach_{2U};
  std::vector<CompiledContribution> contributions_;
  RevisionToken geometry_revision_{};
  RevisionToken boundary_revision_{};
  PlanFingerprint thermodynamics_fingerprint_{};
  PlanFingerprint transport_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

class SpeciesEquationPlan {
 public:
  SpeciesEquationPlan() noexcept = default;
  SpeciesEquationPlan(const SpeciesEquationPlan&) = delete;
  SpeciesEquationPlan& operator=(const SpeciesEquationPlan&) = delete;
  SpeciesEquationPlan(SpeciesEquationPlan&&) = delete;
  SpeciesEquationPlan& operator=(SpeciesEquationPlan&&) = delete;

  std::size_t size() const noexcept { return specs_.size(); }
  const ScalarEquationSpec* spec(std::size_t index) const noexcept {
    return index < specs_.size() ? &specs_[index] : nullptr;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Int3 cells() const noexcept { return cells_; }

 private:
  friend class EquationPlanSet;
  friend Status assemble_species(
      const SpeciesEquationPlan&, std::size_t, const EquationStateView&,
      const EquationMaterialView&, Span<const EquationContributionView>,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&) noexcept;
  friend Status assemble_species_impl(
      const SpeciesEquationPlan&, std::size_t, const EquationStateView&,
      const EquationMaterialView&, Span<const EquationContributionView>,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&, bool) noexcept;
  friend Status evaluate_thermophysical_rates(
      const EquationPlanSet&, const ThermophysicalRateInput&,
      ThermophysicalRateOutput, ThermophysicalRateCertificate&) noexcept;
  const CartesianKernelPlan* kernels_{};
  Int3 cells_{};
  FieldId density_{};
  ConvectionScheme convection_{ConvectionScheme::tvd2};
  std::uint8_t convection_reach_{2U};
  std::vector<ScalarEquationSpec> specs_;
  std::vector<std::size_t> contribution_counts_;
  std::vector<std::size_t> contribution_begins_;
  std::vector<CompiledContribution> contributions_;
  RevisionToken geometry_revision_{};
  RevisionToken boundary_revision_{};
  PlanFingerprint thermodynamics_fingerprint_{};
  PlanFingerprint transport_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

class ScalarEquationPlan {
 public:
  ScalarEquationPlan() noexcept = default;
  ScalarEquationPlan(const ScalarEquationPlan&) = delete;
  ScalarEquationPlan& operator=(const ScalarEquationPlan&) = delete;
  ScalarEquationPlan(ScalarEquationPlan&&) = delete;
  ScalarEquationPlan& operator=(ScalarEquationPlan&&) = delete;

  std::size_t size() const noexcept { return specs_.size(); }
  const ScalarEquationSpec* spec(std::size_t index) const noexcept {
    return index < specs_.size() ? &specs_[index] : nullptr;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Int3 cells() const noexcept { return cells_; }

 private:
  friend class EquationPlanSet;
  friend Status assemble_scalar(
      const ScalarEquationPlan&, std::size_t, const EquationStateView&,
      const EquationMaterialView&, Span<const EquationContributionView>,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&) noexcept;
  friend Status assemble_scalar_impl(
      const ScalarEquationPlan&, std::size_t, const EquationStateView&,
      const EquationMaterialView&, Span<const EquationContributionView>,
      const EquationAssemblyContext&, EquationSystemView,
      EquationAssemblyCertificate&, bool) noexcept;
  friend Status evaluate_thermophysical_rates(
      const EquationPlanSet&, const ThermophysicalRateInput&,
      ThermophysicalRateOutput, ThermophysicalRateCertificate&) noexcept;
  const CartesianKernelPlan* kernels_{};
  Int3 cells_{};
  FieldId density_{};
  ConvectionScheme convection_{ConvectionScheme::tvd2};
  std::uint8_t convection_reach_{2U};
  std::vector<ScalarEquationSpec> specs_;
  std::vector<std::size_t> contribution_counts_;
  std::vector<std::size_t> contribution_begins_;
  std::vector<CompiledContribution> contributions_;
  RevisionToken geometry_revision_{};
  RevisionToken boundary_revision_{};
  PlanFingerprint thermodynamics_fingerprint_{};
  PlanFingerprint transport_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

class AssemblyEpoch {
 public:
  static constexpr std::size_t maximum_tiles = 128U;

  Status begin(const EquationAssemblyContext& context,
               EquationSystemView system) noexcept;
  Status finalize(EquationAssemblyCertificate& certificate) noexcept;
  bool active() const noexcept { return active_; }
  std::size_t completed_tiles() const noexcept { return tile_count_; }

 private:
  friend Status assemble_tile(
      AssemblyEpoch&, const ContinuityEquationPlan&,
      const EquationStateView&, KernelBox) noexcept;
  friend Status assemble_tile(
      AssemblyEpoch&, const MomentumEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, KernelBox) noexcept;
  friend Status assemble_tile(
      AssemblyEpoch&, const MomentumEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, KernelBox, FieldView) noexcept;
  friend Status assemble_tile(
      AssemblyEpoch&, const EnthalpyEquationPlan&, const EquationStateView&,
      const EquationMaterialView&, ConstFieldView,
      Span<const EquationContributionView>, KernelBox) noexcept;
  friend Status assemble_tile(
      AssemblyEpoch&, const SpeciesEquationPlan&, std::size_t,
      const EquationStateView&, const EquationMaterialView&,
      Span<const EquationContributionView>, KernelBox) noexcept;
  friend Status assemble_tile(
      AssemblyEpoch&, const ScalarEquationPlan&, std::size_t,
      const EquationStateView&, const EquationMaterialView&,
      Span<const EquationContributionView>, KernelBox) noexcept;

  Status record(KernelBox box, Int3 cells,
                const EquationAssemblyCertificate& certificate) noexcept;
  Status fail(Status status) noexcept;
  struct TileRecord {
    std::int32_t begin_x;
    std::int32_t begin_y;
    std::int32_t begin_z;
    std::int32_t cells_x;
    std::int32_t cells_y;
    std::int32_t cells_z;
  };
  EquationAssemblyContext context_{};
  EquationSystemView system_{};
  // Only [0,tile_count_) is live; avoiding value-initialization keeps the
  // common one-tile path from clearing 3 KiB of unused schedule storage.
  std::array<TileRecord, maximum_tiles> tiles_;
  EquationAssemblyCertificate candidate_{};
  Int3 cells_{};
  std::uint64_t covered_cells_{};
  std::size_t tile_count_{};
  Status failure_{};
  bool active_{};
};

Status assemble_tile(AssemblyEpoch& epoch,
                     const ContinuityEquationPlan& plan,
                     const EquationStateView& state, KernelBox box) noexcept;
Status assemble_tile(
    AssemblyEpoch& epoch, const MomentumEquationPlan& plan,
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept;
Status assemble_tile(
    AssemblyEpoch& epoch, const EnthalpyEquationPlan& plan,
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept;
Status assemble_tile(
    AssemblyEpoch& epoch, const SpeciesEquationPlan& plan,
    std::size_t species, const EquationStateView& state,
    const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept;
Status assemble_tile(
    AssemblyEpoch& epoch, const ScalarEquationPlan& plan,
    std::size_t scalar, const EquationStateView& state,
    const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept;

class EquationPlanSet {
 public:
  EquationPlanSet() noexcept = default;
  EquationPlanSet(const EquationPlanSet&) = delete;
  EquationPlanSet& operator=(const EquationPlanSet&) = delete;
  EquationPlanSet(EquationPlanSet&& other) noexcept;
  EquationPlanSet& operator=(EquationPlanSet&& other) noexcept;

  static Status compile(
      MPI_Comm communicator, const SchemePlan& schemes,
      const CartesianGeometryPlan& geometry, const MeshPatch& patch,
      const BoundaryPlan& boundary, const ContributionRegistry& contributions,
      const ThermodynamicsPlan& thermodynamics, const TransportPlan& transport,
      const EquationPlanSpec& spec, EquationPlanSet& out,
      EquationCompileDiagnostics* diagnostics = nullptr) noexcept;

  const CartesianKernelPlan& kernels() const noexcept { return kernels_; }
  const ContinuityEquationPlan& continuity() const noexcept {
    return continuity_;
  }
  const MomentumEquationPlan& momentum() const noexcept { return momentum_; }
  const EnthalpyEquationPlan& enthalpy() const noexcept { return enthalpy_; }
  const SpeciesEquationPlan& species() const noexcept { return species_; }
  const ScalarEquationPlan& scalars() const noexcept { return scalars_; }
  const PressureReferencePlan& pressure_reference() const noexcept {
    return pressure_reference_;
  }
  const ThermophysicalPredictorPlan& thermophysical_predictor() const noexcept {
    return thermophysical_predictor_;
  }
  Int3 global_cells() const noexcept { return global_cells_; }
  std::size_t local_cells() const noexcept { return local_cells_; }
  PlanFingerprint semantic_fingerprint() const noexcept {
    return semantic_fingerprint_;
  }
  PlanFingerprint thermodynamics_fingerprint() const noexcept {
    return thermodynamics_fingerprint_;
  }
  PlanFingerprint transport_fingerprint() const noexcept {
    return transport_fingerprint_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  void reset() noexcept;
  void move_from(EquationPlanSet&& other) noexcept;
  void rebind() noexcept;

  CartesianKernelPlan kernels_;
  ContinuityEquationPlan continuity_;
  MomentumEquationPlan momentum_;
  EnthalpyEquationPlan enthalpy_;
  SpeciesEquationPlan species_;
  ScalarEquationPlan scalars_;
  PressureReferencePlan pressure_reference_;
  ThermophysicalPredictorPlan thermophysical_predictor_;
  Int3 global_cells_{};
  Int3 patch_begin_{};
  std::size_t local_cells_{};
  PlanFingerprint semantic_fingerprint_{};
  PlanFingerprint thermodynamics_fingerprint_{};
  PlanFingerprint transport_fingerprint_{};
  PlanFingerprint fingerprint_{};
};

// Thermodynamic part of the target-time continuity/energy Jacobian.  The
// density derivatives come from ThermodynamicsPlan; this helper only forms
// the derivatives of q=rho*h-p_abs used by the coupled energy residual.
struct PressureEnergyThermoJacobian {
  double density{};
  double drho_dp_hY{};
  double drho_dh_pY{};
  double dq_dp_hY{};
  double dq_dh_pY{};
};

Status form_pressure_energy_thermo_jacobian(
    double pressure_absolute, double enthalpy, const ThermoState& state,
    PressureEnergyThermoJacobian& jacobian) noexcept;

struct PressureEnergyTemporalPoint {
  BdfCoefficients bdf{};
  double cell_volume{};
  double pressure_absolute{};
  double density{};
  double enthalpy{};
  double accepted_pressure_absolute{};
  double accepted_density{};
  double accepted_enthalpy{};
  double previous_pressure_absolute{};
  double previous_density{};
  double previous_enthalpy{};
  PressureEnergyThermoJacobian target_thermo{};
};

struct PressureEnergyTemporalLinearization {
  double target_q{};
  double accepted_q{};
  double previous_q{};
  double continuity_residual{};
  double energy_residual{};
  double continuity_pressure{};
  double continuity_enthalpy{};
  double energy_pressure{};
  double energy_enthalpy{};
};

Status linearize_pressure_energy_temporal(
    const PressureEnergyTemporalPoint& point,
    PressureEnergyTemporalLinearization& linearization) noexcept;

inline constexpr std::size_t kPressureEnergyGlobalizationCandidateCount = 25U;
inline constexpr double kPressureEnergyGlobalizationArmijoCoefficient = 1.0e-4;

// A pure, already-globally-reduced sample.  The policy below neither evaluates
// candidate fields nor performs MPI communication; its caller is responsible
// for synchronously producing p/h/rho/T/U and final mass-flux provenance.
// Merit is the Euclidean norm of the two normalized residuals.  This permits
// a coupled Newton direction to exchange residual between blocks while still
// requiring joint descent.  A candidate must satisfy
// merit(alpha) <= (1-c*alpha)*merit(0), as well as strict decrease.
struct PressureEnergyGlobalizationSample {
  double alpha{};
  double global_normalized_continuity{};
  double global_normalized_energy{};
  bool thermodynamically_admissible{};
  bool state_and_flux_finite{};
  std::uint8_t corrector{};
  RevisionToken target_time{};
  PlanFingerprint correction_direction{};
  PlanFingerprint state_provenance{};
  PlanFingerprint mass_flux_provenance{};
};

enum class PressureEnergyGlobalizationScope : std::uint8_t {
  invalid,
  frozen_momentum_continuity_energy_globalization,
};

// This certificate signs only globalization of one frozen-momentum
// continuity-energy correction.  It is deliberately not a full-Newton
// certificate and does not authorize publishing candidate fields by itself.
struct PressureEnergyGlobalizationSelectionCertificate {
  PressureEnergyGlobalizationScope scope{
      PressureEnergyGlobalizationScope::invalid};
  double alpha{};
  double baseline_normalized_continuity{};
  double baseline_normalized_energy{};
  double candidate_normalized_continuity{};
  double candidate_normalized_energy{};
  double baseline_merit{};
  double candidate_merit{};
  double armijo_upper_bound{};
  RevisionToken target_time{};
  PlanFingerprint correction_direction{};
  PlanFingerprint baseline_state_provenance{};
  PlanFingerprint baseline_mass_flux_provenance{};
  PlanFingerprint candidate_state_provenance{};
  PlanFingerprint candidate_mass_flux_provenance{};
  PlanFingerprint selection_provenance{};
  std::uint8_t corrector{};
  std::uint8_t selected_halvings{};
  bool thermodynamically_admissible{};
  bool state_and_flux_finite{};
  bool strict_merit_decrease{};
  bool armijo_sufficient_decrease{};
  bool full_nonlinear_newton{};

  bool valid() const noexcept;
};

// Deterministically considers an already-evaluated contiguous prefix of
// alpha=1,1/2,...,2^-24 (between one and 25 samples).  Product evaluates in
// that order and stops at the first acceptable candidate.  Malformed or
// foreign sample sets are invalid plans; an empty prefix or a well-formed
// evaluated prefix with no acceptable candidate returns rejected_step and no
// certificate.
Status select_pressure_energy_globalization(
    const PressureEnergyGlobalizationSample& baseline,
    Span<const PressureEnergyGlobalizationSample> candidates,
    PressureEnergyGlobalizationSelectionCertificate& certificate) noexcept;

enum class PressureEnergySchurSignClass : std::uint8_t { general };

// Immutable 0/1 cell activity certificate used by the coupled block.  Empty
// selects the ordinary Cartesian domain.  A non-empty map uses one for fluid
// rows and zero for IBM solid rows.
struct PressureEnergyCellActivity {
  Span<const std::uint8_t> cells{};
  PlanFingerprint local_fingerprint{};
  PlanFingerprint collective_fingerprint{};
};

// Read-only authority for decomposing one raw pressure correction as
//
//   dp_abs = (dp_raw - shift) + (p_ref,next - p_ref,old).
//
// The candidate compressibility is evaluated at the same target layer as the
// raw correction.  Only active fluid cells contribute to the weighted gauge,
// while absolute pressure is checked in every local cell (including IBM
// inactive rows) before a transaction can be certified.
struct ClosedGaugeCorrectionPrepareInput {
  PressureReferenceCertificate predecessor{};
  double pressure_reference{};
  std::uint8_t corrector{};
  RevisionToken time{};
  RevisionToken geometry{};
  // Exact rank-local producer lineages.  Their presence enters the collective
  // contract, while their values enter the rank-local transaction.
  RevisionToken pressure_correction_authority{};
  PlanFingerprint target_thermodynamic_closure{};
  ConstFieldView pressure_perturbation{};
  ConstFieldView raw_pressure_correction{};
  ConstFieldView candidate_pressure_compressibility{};
  PressureEnergyCellActivity activity{};
};

// The collective token contains only decomposition-independent semantics.
// The rank-local token additionally binds concrete pointers, layouts,
// storage/revision domains, exact revisions, cell values, and the local IBM
// activity map.  output_pressure_reference.pressure_reference is the latter:
// pressure-reference consumers are rank-local field transactions, while the
// collective token proves that all ranks prepared one common gauge shift.
struct ClosedGaugeCorrectionCertificate {
  double shift{};
  double next_pressure_reference{};
  double local_moment{};
  double local_weight{};
  double global_moment{};
  double global_weight{};
  double local_post_shift_moment{};
  double local_post_shift_absolute_moment{};
  double global_post_shift_moment{};
  double global_post_shift_absolute_moment{};
  double post_shift_gauge_residual{};
  double post_shift_gauge_tolerance{};
  PressureReferenceCertificate output_pressure_reference{};
  RevisionToken predecessor_pressure_reference{};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken pressure_correction_authority{};
  PlanFingerprint target_thermodynamic_closure{};
  PlanFingerprint activity_local_fingerprint{};
  PlanFingerprint activity_collective_fingerprint{};
  PlanFingerprint collective_transaction{};
  RevisionToken rank_local_transaction{};
  std::size_t local_active_cells{};
  std::uint8_t corrector{};

  bool valid() const noexcept {
    return std::isfinite(shift) &&
           std::isfinite(next_pressure_reference) &&
           next_pressure_reference > 0.0 && std::isfinite(local_moment) &&
           std::isfinite(local_weight) && local_weight >= 0.0 &&
           std::isfinite(global_moment) && std::isfinite(global_weight) &&
           global_weight > 0.0 &&
           std::isfinite(local_post_shift_moment) &&
           std::isfinite(local_post_shift_absolute_moment) &&
           local_post_shift_absolute_moment >= 0.0 &&
           std::isfinite(global_post_shift_moment) &&
           std::isfinite(global_post_shift_absolute_moment) &&
           global_post_shift_absolute_moment >= 0.0 &&
           std::isfinite(post_shift_gauge_residual) &&
           post_shift_gauge_residual >= 0.0 &&
           std::isfinite(post_shift_gauge_tolerance) &&
           post_shift_gauge_tolerance > 0.0 &&
           post_shift_gauge_residual <= post_shift_gauge_tolerance &&
           output_pressure_reference.valid() &&
           predecessor_pressure_reference != 0U && time != 0U &&
           geometry != 0U && pressure_correction_authority != 0U &&
           target_thermodynamic_closure != 0U &&
           collective_transaction != 0U &&
           rank_local_transaction != 0U &&
           output_pressure_reference.pressure_reference ==
               rank_local_transaction &&
           output_pressure_reference.time == time &&
           ((activity_local_fingerprint == 0U) ==
            (activity_collective_fingerprint == 0U)) &&
           (corrector == 1U || corrector == 2U);
  }
};

// Exact active-cell/face map consumed by continuity and both spatial energy
// blocks. Empty spans select the ordinary full Cartesian domain.
struct PressureContinuityActivityView {
  Span<const std::uint8_t> cells{};
  Span<const std::uint8_t> x_faces{};
  Span<const std::uint8_t> y_faces{};
  Span<const std::uint8_t> z_faces{};
  PlanFingerprint local_fingerprint{};
  PlanFingerprint collective_fingerprint{};
};

// IBM-owned, immutable authority for Cartesian physical control faces whose
// owner cell is inactive. It exposes neither EBTopology nor raw activity to
// the physical-boundary engine. The local fingerprint binds the exact local
// topology and boundary planes; the collective fingerprint binds the ordered
// rank-local activities to one communicator-wide topology contract.
class IbmPhysicalBoundaryFluxAuthority {
 public:
  IbmPhysicalBoundaryFluxAuthority() noexcept = default;
  ~IbmPhysicalBoundaryFluxAuthority() noexcept;
  IbmPhysicalBoundaryFluxAuthority(
      const IbmPhysicalBoundaryFluxAuthority&) = delete;
  IbmPhysicalBoundaryFluxAuthority& operator=(
      const IbmPhysicalBoundaryFluxAuthority&) = delete;
  IbmPhysicalBoundaryFluxAuthority(
      IbmPhysicalBoundaryFluxAuthority&&) noexcept;
  IbmPhysicalBoundaryFluxAuthority& operator=(
      IbmPhysicalBoundaryFluxAuthority&&) noexcept;

  static Status compile(
      MPI_Comm communicator, const CartesianGeometryPlan& geometry,
      MeshPatch patch, const EBTopology& topology,
      const IbmEquationInterfacePlan& immersed_interface,
      IbmPhysicalBoundaryFluxAuthority& out) noexcept;
  Status physical_face_active(CartesianAxis axis, Int3 face,
                              bool& active) const noexcept;
  Status zero_inactive_physical_boundary_flux(FaceFluxView flux) const
      noexcept;
  Status validate_inactive_physical_boundary_flux(
      ConstFaceFluxView flux, double absolute_tolerance = 0.0) const noexcept;
  bool matches(MPI_Comm communicator,
               const CartesianGeometryPlan* geometry, MeshPatch patch,
               const IbmEquationInterfacePlan* immersed_interface) const
      noexcept;
  bool ready() const noexcept;
  PlanFingerprint local_fingerprint() const noexcept;
  PlanFingerprint collective_fingerprint() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

struct PressureEnergyDiagonalBinding {
  ConstFieldView diagonal{};
  PressureEnergyCellActivity activity{};
  LinearIdentity identity{};
  // Applied on inactive rows without reading diagonal.  Production Ep uses
  // zero and Eh uses one, forming the closed solid block dp=dh=0.
  double inactive_diagonal{};
};

struct PressureEnergyDiagonalCertificate {
  LinearOperatorCertificate linear{};
  RevisionToken diagonal_revision{};
  double inactive_diagonal{};
  PlanFingerprint activity_local_fingerprint{};
  PlanFingerprint activity_collective_fingerprint{};
  std::size_t active_cells{};
  std::size_t inactive_cells{};
  bool cell_local{};

  bool valid() const noexcept {
    return linear.identity.symbolic != 0U && linear.identity.numeric != 0U &&
           linear.identity.hierarchy != 0U &&
           linear.identity.workspace != 0U &&
           linear.identity.fingerprint != 0U &&
           linear.collective_fingerprint != 0U &&
           linear.local_shape.x > 0 && linear.local_shape.y > 0 &&
           linear.local_shape.z > 0 &&
           linear.operator_class == LinearOperatorClass::nonsymmetric &&
           diagonal_revision != 0U &&
           active_cells + inactive_cells ==
               static_cast<std::size_t>(linear.local_shape.x) *
                   static_cast<std::size_t>(linear.local_shape.y) *
                   static_cast<std::size_t>(linear.local_shape.z) &&
           ((activity_local_fingerprint == 0U) ==
            (activity_collective_fingerprint == 0U)) &&
           std::isfinite(inactive_diagonal) && cell_local;
  }
};

// Exact cell-local action used for the frozen-spatial Ep/Eh blocks.  The
// binding is immutable: a new diagonal revision or activity certificate must
// produce a newly bound operator.
class PressureEnergyDiagonalOperator final : public LinearOperator {
 public:
  PressureEnergyDiagonalOperator() noexcept = default;
  ~PressureEnergyDiagonalOperator() noexcept override = default;
  PressureEnergyDiagonalOperator(const PressureEnergyDiagonalOperator&) =
      delete;
  PressureEnergyDiagonalOperator& operator=(
      const PressureEnergyDiagonalOperator&) = delete;
  PressureEnergyDiagonalOperator(PressureEnergyDiagonalOperator&&) noexcept =
      default;
  PressureEnergyDiagonalOperator& operator=(
      PressureEnergyDiagonalOperator&&) noexcept = default;

  static Status bind(const PressureEnergyDiagonalBinding& binding,
                     PressureEnergyDiagonalOperator& out,
                     PressureEnergyDiagonalCertificate& certificate) noexcept;

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_.linear;
  }
  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return failure_;
  }
  const PressureEnergyDiagonalCertificate& diagonal_certificate()
      const noexcept {
    return certificate_;
  }
  Status apply(FieldView input, FieldView output) const noexcept override;

 private:
  ConstFieldView diagonal_{};
  PressureEnergyCellActivity activity_{};
  double inactive_diagonal_{};
  PressureEnergyDiagonalCertificate certificate_{};
  mutable LinearOperatorFailureProvenance failure_{};
};

// Immutable same-target authorities for the frozen-spatial E_h action.  The
// semantic fingerprints name the producers, while revisions name this
// concrete target layer.  Keeping both prevents a numerically current field
// from being substituted from a foreign thermodynamics/transport plan.
struct PressureEnergyEnthalpyAuthority {
  BdfCoefficients bdf{};
  RevisionToken target_time{};
  RevisionToken geometry{};
  RevisionToken numeric_boundary{};
  RevisionToken thermodynamics{};
  RevisionToken transport{};
  PlanFingerprint equation_semantics{};
  PlanFingerprint thermodynamics_semantics{};
  PlanFingerprint transport_semantics{};
};

struct PressureEnergyEnthalpyServices {
  MPI_Comm communicator{MPI_COMM_NULL};
  HaloEngine* halo{};
  StageId halo_stage{};
  FieldId enthalpy_variation_field{};
  FieldId temperature_variation_field{};
};

struct PressureEnergyEnthalpyWorkspace {
  FieldView delta_temperature{};
  FrozenConvectionFaceOutput directional_enthalpy{};
};

struct PressureEnergyEnthalpyBinding {
  const CartesianGeometryPlan* geometry{};
  const CartesianKernelPlan* kernels{};
  const BoundaryPlan* boundary{};
  MeshPatch patch{};
  ConvectionScheme convection{ConvectionScheme::limited_central2};
  PressureEnergyEnthalpyServices services{};
  PressureEnergyEnthalpyAuthority authority{};
  ConstFieldView assembled_diagonal{};
  ConstFieldView target_enthalpy{};
  ConstFieldView density_enthalpy_derivative{};
  ConstFieldView heat_capacity{};
  ConstFieldView thermal_conductivity{};
  ConstFieldView enthalpy_diffusivity{};
  ConstFaceFluxView target_flux{};
  FrozenConvectionContext convection_context{};
  FrozenConvectionFaceField frozen_face_enthalpy{};
  PressureEnergyEnthalpyWorkspace workspace{};
  PressureContinuityActivityView activity{};
  LinearIdentity identity{};
  FrozenConvectionLinearizationPolicy linearization_policy{
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope};
};

struct PressureEnergyEnthalpyCertificate {
  LinearOperatorCertificate linear{};
  RevisionToken binding_revision{};
  PressureEnergyEnthalpyAuthority authority{};
  RevisionToken assembled_diagonal{};
  RevisionToken target_enthalpy{};
  RevisionToken density_enthalpy_derivative{};
  RevisionToken heat_capacity{};
  RevisionToken thermal_conductivity{};
  RevisionToken enthalpy_diffusivity{};
  RevisionToken target_flux{};
  RevisionToken frozen_face_enthalpy{};
  PlanFingerprint frozen_reconstruction{};
  PlanFingerprint frozen_local_binding{};
  PlanFingerprint directional_reconstruction{};
  PlanFingerprint directional_branch_authority{};
  PlanFingerprint geometry_fingerprint{};
  PlanFingerprint kernels{};
  PlanFingerprint boundary_semantics{};
  PlanFingerprint boundary_rank_local_layout{};
  StageId halo_stage{};
  FieldId enthalpy_variation_field{};
  FieldId temperature_variation_field{};
  RevisionToken delta_temperature_revision{};
  StorageIdentity delta_temperature_storage{};
  RevisionDomainIdentity delta_temperature_revision_domain{};
  StorageIdentity directional_enthalpy_storage{};
  RevisionDomainIdentity directional_enthalpy_revision_domain{};
  std::uintptr_t halo_instance{};
  PlanFingerprint activity_local_fingerprint{};
  PlanFingerprint activity_collective_fingerprint{};
  std::size_t active_cells{};
  std::size_t inactive_cells{};
  std::uint64_t generalized_face_count{};
  FrozenConvectionLinearizationPolicy linearization_policy{
      FrozenConvectionLinearizationPolicy::classical_active_branch};
  bool exact_cartesian_spatial_response{};
  bool exact_temperature_space_conduction{};
  bool ibm_spatial_derivative{};
  bool inactive_rows_identity{};
  bool inactive_interfaces_zero{};
  bool allocation_free_apply{};

  bool valid() const noexcept;
};

// Complete frozen-spatial target-layer energy response
//
//   E_h dh = [A_h - diag(lambda/cp) + a0 V h rho_h] dh
//            + sum_f(outward phi_f * dh_f)
//            - V div(lambda grad(dh/cp)).
//
// A_h is the diagonal published by assemble_enthalpy.  The old scalar
// lambda/cp proxy is removed before the exact temperature-space response is
// added.  Apply exchanges dh(reach=2) and deltaT(reach=1) together and does
// not allocate.
class PressureEnergyEnthalpyOperator final : public LinearOperator {
 public:
  PressureEnergyEnthalpyOperator() noexcept = default;
  ~PressureEnergyEnthalpyOperator() noexcept override = default;
  PressureEnergyEnthalpyOperator(const PressureEnergyEnthalpyOperator&) =
      delete;
  PressureEnergyEnthalpyOperator& operator=(
      const PressureEnergyEnthalpyOperator&) = delete;
  PressureEnergyEnthalpyOperator(PressureEnergyEnthalpyOperator&&) noexcept =
      default;
  PressureEnergyEnthalpyOperator& operator=(
      PressureEnergyEnthalpyOperator&&) noexcept = default;

  static Status bind(const PressureEnergyEnthalpyBinding& binding,
                     PressureEnergyEnthalpyOperator& out,
                     PressureEnergyEnthalpyCertificate& certificate) noexcept;

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_.linear;
  }
  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return failure_;
  }
  const PressureEnergyEnthalpyCertificate& enthalpy_certificate()
      const noexcept {
    return certificate_;
  }
  Status apply(FieldView input, FieldView output) const noexcept override;

 private:
  const CartesianGeometryPlan* geometry_{};
  const CartesianKernelPlan* kernels_{};
  const BoundaryPlan* boundary_{};
  MeshPatch patch_{};
  ConvectionScheme convection_{ConvectionScheme::limited_central2};
  PressureEnergyEnthalpyServices services_{};
  ConstFieldView assembled_diagonal_{};
  ConstFieldView target_enthalpy_{};
  ConstFieldView density_enthalpy_derivative_{};
  ConstFieldView heat_capacity_{};
  ConstFieldView thermal_conductivity_{};
  ConstFieldView enthalpy_diffusivity_{};
  ConstFaceFluxView target_flux_{};
  FrozenConvectionContext convection_context_{};
  FrozenConvectionFaceField frozen_face_enthalpy_{};
  PressureEnergyEnthalpyWorkspace workspace_{};
  PressureContinuityActivityView activity_{};
  PressureEnergyEnthalpyCertificate certificate_{};
  mutable LinearOperatorFailureProvenance failure_{};
};

struct PressureEnergySchurWorkspace {
  FieldView continuity_response{};
  FieldView eliminated_enthalpy{};
  FieldView energy_response{};
};

struct PressureEnergyPressureFluxCertificate;
class PressureEnergyPressureFluxOperator;

enum class PressureEnergySchurBlockScope : std::uint8_t {
  invalid,
  exact_cartesian_frozen_spatial,
  ibm_cartesian_spatial_quasi_newton,
  ibm_double_diagonal_quasi_newton,
  generic_algebraic_quasi_newton,
};

// Typed provenance for the two energy blocks.  Production authorities are
// derived only from complete Cartesian E_p/E_h certificates, the masked
// Cartesian IBM spatial pair, or the two legacy IBM fallback diagonals
// (inactive E_p=0, inactive E_h=1).  No raw boolean can promote a flux-only
// Cartesian block to an exact Jacobian.
class PressureEnergySchurBlockAuthority {
 public:
  PressureEnergySchurBlockAuthority() noexcept = default;

  static Status exact_cartesian(
      const PressureEnergyPressureFluxOperator& energy_pressure,
      const PressureEnergyEnthalpyOperator& energy_enthalpy,
      PressureEnergySchurBlockAuthority& out) noexcept;
  static Status ibm_cartesian_spatial_quasi_newton(
      const PressureEnergyPressureFluxOperator& energy_pressure,
      const PressureEnergyEnthalpyOperator& energy_enthalpy,
      PressureEnergySchurBlockAuthority& out) noexcept;
  static Status ibm_double_diagonal(
      const PressureEnergyDiagonalOperator& energy_pressure,
      const PressureEnergyDiagonalOperator& energy_enthalpy,
      PressureEnergySchurBlockAuthority& out) noexcept;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  static Status generic_algebraic_quasi_newton_for_test(
      LinearOperatorCertificate energy_pressure,
      LinearOperatorCertificate energy_enthalpy,
      PressureEnergyCellActivity activity,
      PressureEnergySchurBlockAuthority& out) noexcept;
#endif

  bool valid() const noexcept;
  PressureEnergySchurBlockScope scope() const noexcept { return scope_; }

 private:
  friend class PressureEnergySchurOperator;

  bool matches_operators(const LinearOperator* energy_pressure,
                         const LinearOperator* energy_enthalpy) const noexcept;

  PressureEnergySchurBlockScope scope_{
      PressureEnergySchurBlockScope::invalid};
  const LinearOperator* energy_pressure_operator_{};
  const LinearOperator* energy_enthalpy_operator_{};
  LinearOperatorCertificate energy_pressure_{};
  LinearOperatorCertificate energy_enthalpy_{};
  PlanFingerprint activity_local_fingerprint_{};
  PlanFingerprint activity_collective_fingerprint_{};
  PlanFingerprint collective_fingerprint_{};
  RevisionToken rank_local_revision_{};
};

// C_h is a certified cell-local diagonal.  The other three block actions are
// exact matrix-free operators over the same target state:
//
//   [ C_p  C_h ] [dp] = -[R_C]
//   [ E_p  E_h ] [dh]   [R_E]
//
// The bound operator applies S_p=E_p-E_h inv(C_h) C_p.  Native pressure MG,
// when retained, is only an approximate preconditioner for this operator.
struct PressureEnergySchurBinding {
  const LinearOperator* continuity_pressure{};
  const LinearOperator* energy_pressure{};
  const LinearOperator* energy_enthalpy{};
  ConstFieldView continuity_enthalpy_diagonal{};
  ConstFieldView continuity_enthalpy_row_scale{};
  PressureEnergySchurWorkspace workspace{};
  PressureEnergyCellActivity activity{};
  double scaled_pivot_floor{};
  PressureEnergySchurBlockAuthority block_authority{};
};

enum class PressureEnergyJacobianScope : std::uint8_t {
  exact_cartesian_frozen_spatial,
  ibm_cartesian_spatial_quasi_newton,
  ibm_double_diagonal_quasi_newton,
  generic_algebraic_quasi_newton,
};

struct PressureEnergyJacobianCertificate {
  LinearOperatorCertificate schur{};
  RevisionToken block_jacobian{};
  RevisionToken block_scope_authority{};
  RevisionToken continuity_enthalpy{};
  double minimum_scaled_abs_c_h{};
  double maximum_scaled_abs_c_h{};
  std::size_t active_cells{};
  std::size_t inactive_cells{};
  PlanFingerprint activity_local_fingerprint{};
  PlanFingerprint activity_collective_fingerprint{};
  PressureEnergySchurSignClass sign_class{
      PressureEnergySchurSignClass::general};
  PressureEnergyJacobianScope jacobian_scope{
      PressureEnergyJacobianScope::generic_algebraic_quasi_newton};
  // "Exact" is deliberately scoped to the assembled quasi-Newton block.  A
  // terminal nonlinear residual audit is still required.
  bool exact_algebraic_schur{};
  bool full_nonlinear_jacobian{};
  bool exact_block_equivalent{};
  bool cell_local_continuity_enthalpy{};
  bool native_mg_preconditioner_only{};

  bool valid() const noexcept {
    return schur.identity.symbolic != 0U && schur.identity.numeric != 0U &&
           schur.identity.hierarchy != 0U &&
           schur.identity.workspace != 0U &&
           schur.identity.fingerprint != 0U &&
           schur.collective_fingerprint != 0U &&
           schur.local_shape.x > 0 && schur.local_shape.y > 0 &&
           schur.local_shape.z > 0 &&
           schur.operator_class == LinearOperatorClass::nonsymmetric &&
           block_jacobian != 0U && block_scope_authority != 0U &&
           continuity_enthalpy != 0U &&
           active_cells + inactive_cells ==
               static_cast<std::size_t>(schur.local_shape.x) *
                   static_cast<std::size_t>(schur.local_shape.y) *
                   static_cast<std::size_t>(schur.local_shape.z) &&
           ((activity_local_fingerprint == 0U) ==
            (activity_collective_fingerprint == 0U)) &&
           ((active_cells == 0U && minimum_scaled_abs_c_h == 0.0 &&
             maximum_scaled_abs_c_h == 0.0) ||
            (active_cells > 0U && minimum_scaled_abs_c_h > 0.0 &&
             maximum_scaled_abs_c_h >= minimum_scaled_abs_c_h)) &&
           (jacobian_scope ==
                PressureEnergyJacobianScope::
                    exact_cartesian_frozen_spatial ||
            jacobian_scope == PressureEnergyJacobianScope::
                                  ibm_cartesian_spatial_quasi_newton ||
            jacobian_scope == PressureEnergyJacobianScope::
                                  ibm_double_diagonal_quasi_newton ||
            jacobian_scope == PressureEnergyJacobianScope::
                                  generic_algebraic_quasi_newton) &&
           exact_algebraic_schur && !full_nonlinear_jacobian &&
           exact_block_equivalent && cell_local_continuity_enthalpy &&
           native_mg_preconditioner_only;
  }
};

class PressureEnergySchurOperator final : public LinearOperator {
 public:
  PressureEnergySchurOperator() noexcept = default;
  ~PressureEnergySchurOperator() noexcept override = default;
  PressureEnergySchurOperator(const PressureEnergySchurOperator&) = delete;
  PressureEnergySchurOperator& operator=(
      const PressureEnergySchurOperator&) = delete;
  PressureEnergySchurOperator(PressureEnergySchurOperator&&) noexcept =
      default;
  PressureEnergySchurOperator& operator=(
      PressureEnergySchurOperator&&) noexcept = default;

  static Status bind(const PressureEnergySchurBinding& binding,
                     PressureEnergySchurOperator& out,
                     PressureEnergyJacobianCertificate& certificate) noexcept;

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_.schur;
  }
  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return failure_;
  }
  const PressureEnergyJacobianCertificate& jacobian_certificate()
      const noexcept {
    return certificate_;
  }

  Status apply(FieldView pressure, FieldView output) const noexcept override;
  Status form_pressure_rhs(ConstFieldView continuity_residual,
                           ConstFieldView energy_residual,
                           FieldView output) const noexcept;
  Status recover_enthalpy(ConstFieldView continuity_residual,
                          FieldView pressure_correction,
                          FieldView enthalpy_correction) const noexcept;
  // PressureCorrectionSystemView stores rhs=-R_C.  These variants consume it
  // directly, so production does not have to overwrite or reconstruct R_C.
  Status form_pressure_rhs_from_continuity_system_rhs(
      ConstFieldView continuity_system_rhs,
      ConstFieldView energy_residual, FieldView output) const noexcept;
  Status recover_enthalpy_from_continuity_system_rhs(
      ConstFieldView continuity_system_rhs, FieldView pressure_correction,
      FieldView enthalpy_correction) const noexcept;

 private:
  const LinearOperator* continuity_pressure_{};
  const LinearOperator* energy_pressure_{};
  const LinearOperator* energy_enthalpy_{};
  LinearOperatorCertificate continuity_pressure_certificate_{};
  LinearOperatorCertificate energy_pressure_certificate_{};
  LinearOperatorCertificate energy_enthalpy_certificate_{};
  ConstFieldView continuity_enthalpy_diagonal_{};
  ConstFieldView continuity_enthalpy_row_scale_{};
  PressureEnergySchurWorkspace workspace_{};
  PressureEnergyCellActivity activity_{};
  PressureEnergyJacobianCertificate certificate_{};
  mutable LinearOperatorFailureProvenance failure_{};
};

struct PisoPlanSpec {
  std::uint8_t pressure_correctors{2U};
  StageId pressure_stage{};
  RevisionSlotId final_flux_slot{};
  LinearAlgorithm pressure_algorithm{LinearAlgorithm::fgmres};
  MgCorrectionScaling mg_correction_scaling{
      MgCorrectionScaling::residual_minimizing};
  LinearSolveControl pressure_solve{};
  double eos_tolerance{};
  double continuity_tolerance{};
  // Zero preserves the legacy four-gate terminal audit. A positive value
  // enables the same-target normalized energy-residual gate.
  double energy_tolerance{};
  double closed_mass_tolerance{};
  double gauge_tolerance{};
};

struct PisoCompileDiagnostics {
  int lowest_failing_rank{-1};
};

class PisoPlan {
 public:
  PisoPlan() noexcept = default;
  PisoPlan(const PisoPlan&) = delete;
  PisoPlan& operator=(const PisoPlan&) = delete;
  PisoPlan(PisoPlan&&) noexcept = default;
  PisoPlan& operator=(PisoPlan&&) noexcept = default;

  static Status compile(MPI_Comm communicator,
                        const EquationPlanSet& equations,
                        const PisoPlanSpec& spec, PisoPlan& out,
                        PisoCompileDiagnostics* diagnostics = nullptr) noexcept;
  std::uint8_t pressure_correctors() const noexcept { return 2U; }
  LinearAlgorithm pressure_algorithm() const noexcept {
    return pressure_algorithm_;
  }
  MgCorrectionScaling pressure_correction_scaling() const noexcept {
    return mg_correction_scaling_;
  }
  StageId pressure_stage() const noexcept { return pressure_stage_; }
  RevisionSlotId final_flux_slot() const noexcept { return final_flux_slot_; }
  const LinearSolveControl& pressure_solve() const noexcept {
    return pressure_solve_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class PressureVelocityCoupler;
  Int3 cells_{};
  PlanFingerprint equations_fingerprint_{};
  PlanFingerprint predictor_fingerprint_{};
  PlanFingerprint pressure_reference_fingerprint_{};
  StageId pressure_stage_{};
  RevisionSlotId final_flux_slot_{};
  LinearAlgorithm pressure_algorithm_{LinearAlgorithm::fgmres};
  MgCorrectionScaling mg_correction_scaling_{
      MgCorrectionScaling::residual_minimizing};
  LinearSolveControl pressure_solve_{};
  double eos_tolerance_{};
  double continuity_tolerance_{};
  double energy_tolerance_{};
  double closed_mass_tolerance_{};
  double gauge_tolerance_{};
  PlanFingerprint fingerprint_{};
};

struct PisoContinuityWitness {
  bool valid{};
  std::uint64_t global_cell{};
  Int3 global_index{};
  std::int32_t rank{-1};
  double normalized_residual{};
  double raw_balance{};
  double unsteady{};
  double flux_divergence{};
  double scale{};
  double density{};
  double accepted_density{};
  double previous_density{};
  std::array<double, 6U> face_fluxes{};
  bool global_domain_edge{};
};

// Extra same-target pressure--enthalpy solves are nonlinear refinements of
// PISO corrector two, not additional PISO correctors.  The capacity is a
// compile-time hot-resource contract and is deliberately not case input.
inline constexpr std::size_t kPressureEnergyRefinementCapacity = 6U;

enum class PressureEnergyRefinementTermination : std::uint8_t {
  none,
  component_residuals_converged,
  iteration_capacity_exhausted,
  rejected_candidate,
};

struct PisoPressureEnergyRefinementSolveReport {
  LinearSolveResult solve{};
  // Same physical target generation as the ordinary C2 solve.
  RevisionToken target_generation{};
  // Rank-invariant lineage of the exact candidate state/flux that generated
  // this refinement system.
  PlanFingerprint collective_lineage{};
  // Rank-local pressure-system and linear-lifecycle identities.  These are
  // diagnostic capabilities and must never be compared across ranks.
  RevisionToken pressure_state{};
  LinearIdentity linear_identity{};
  std::uint8_t ordinal{};

  bool valid() const noexcept {
    return target_generation != 0U && collective_lineage != 0U &&
           pressure_state != 0U && linear_identity.symbolic != 0U &&
           linear_identity.numeric != 0U && linear_identity.hierarchy != 0U &&
           linear_identity.workspace != 0U &&
           linear_identity.fingerprint != 0U && ordinal != 0U &&
           ordinal <= kPressureEnergyRefinementCapacity;
  }
};

struct PisoAttemptReport {
  std::array<LinearSolveResult, 2U> pressure{};
  std::array<PisoPressureEnergyRefinementSolveReport,
             kPressureEnergyRefinementCapacity>
      pressure_energy_refinement{};
  double eos_residual{};
  double continuity_residual{};
  double energy_residual{};
  double closed_mass_residual{};
  double gauge_residual{};
  RevisionToken final_flux_revision{};
  std::uint8_t pressure_solve_calls{};
  std::uint8_t pressure_energy_refinement_solve_calls{};
  PressureEnergyRefinementTermination pressure_energy_refinement_termination{
      PressureEnergyRefinementTermination::none};
  PisoContinuityWitness continuity_witness{};
};

struct PisoCommittedFaceFluxHistory {
  // Certified committed mass flux at t_n. Corrector one always requires it.
  ConstFaceFluxView accepted{};
  // Certified committed mass flux at t_{n-1}. It is present only for BDF2.
  ConstFaceFluxView previous{};
};

struct PisoIntermediateInput {
  EquationAssemblyCertificate momentum{};
  ThermophysicalPredictorCertificate predictor{};
  PressureReferenceCertificate pressure_reference{};
  FieldView density{};
  ConstFieldView trial_velocity{};
  // Current attempt-local face mass flux. Corrector one consumes the
  // momentum-predictor revision; corrector two must name the distinct
  // revision published by the first pressure correction.
  ConstFaceFluxView trial_flux{};
  EquationSystemView momentum_system{};
  BdfCoefficients bdf{};
  RevisionToken numeric_boundary{};
  BoundaryResolvedValues boundary_values{};
  RevisionToken prior_corrector{};
  std::uint8_t corrector{};
  const IbmEquationInterfacePlan* immersed_interface{};
  // Corrector one consumes this attempt-bound temporal mass-flux reference
  // in place.  It must name the coupler's replica-0 storage.  Corrector two
  // leaves the view empty and consumes the C1-corrected trial_flux.
  ConstFaceFluxView temporal_reference{};
  // Corrector one forms its internal Rhie-Chow temporal correction from the
  // normalized BDF history of these committed fluxes. The attempt-local
  // trial_flux remains authoritative only for fixed physical boundaries.
  // Corrector two leaves both views empty.
  PisoCommittedFaceFluxHistory committed_face_history{};
  BoundaryThermophysicalGhostUse thermophysical_boundary{};
};

struct PisoCouplerWorkspace {
  // Cell-volume-scaled reciprocal of the integrated momentum diagonal,
  // V/A_U.  This makes r_au*grad(p) a velocity and gives the pressure-face
  // coefficient mass-flow-per-pressure units on non-unit cells.
  FieldView r_au{};
  FieldView h_by_a{};
  FieldView pressure_gradient{};
  FaceFieldView x_pressure_coefficient{};
  FaceFieldView y_pressure_coefficient{};
  FaceFieldView z_pressure_coefficient{};
  FaceFluxView phi_h_by_a{};
};

struct PisoCouplerServices {
  MPI_Comm communicator{MPI_COMM_NULL};
  const CartesianGeometryPlan* geometry{};
  MeshPatch patch{};
  const BoundaryPlan* boundary{};
  // Retained cold dependency used by the exact-candidate issuer to replay
  // periodic EOS closure rather than trusting caller-provided rho/T.
  const ThermodynamicsPlan* thermodynamics{};
  HaloEngine* halo{};
  StageId halo_stage{};
  FieldId density_field{};
  HaloEngine* correction_halo{};
  StageId correction_halo_stage{};
  FieldId correction_field{};
  PressureContinuityActivityView continuity_activity{};
  RemoteDonorExchangePlan* pressure_correction_donors{};
  StageId pressure_correction_donor_stage{};
  // Same-target candidates must never alias the live pressure-correction
  // donor plan: RemoteDonorExchangePlan seals FieldId and StageId.  The
  // explicit reach/fingerprint make the injected candidate gather part of
  // the coupler's frozen mechanical contract.
  RemoteDonorExchangePlan* candidate_pressure_correction_donors{};
  StageId candidate_pressure_correction_donor_stage{};
  FieldId candidate_pressure_correction_field{};
  std::uint8_t candidate_pressure_correction_donor_reach{};
  PlanFingerprint candidate_pressure_correction_donor_fingerprint{};
  // IBM identity is a cold mechanical dependency. Candidate boundary
  // finalization is also cold-bound, so refresh() may not discover it after
  // allocation and resource seals have been issued.
  const IbmEquationInterfacePlan* immersed_interface{};
};

struct PisoIntermediateCertificate {
  PlanFingerprint plan{};
  RevisionToken r_au{};
  RevisionToken h_by_a{};
  RevisionToken pressure_face_coefficient{};
  RevisionToken phi_h_by_a{};
  RevisionToken trial_face_flux{};
  RevisionToken temporal_face_flux{};
  RevisionToken committed_face_history{};
  RevisionToken dependency{};
  std::uint8_t corrector{};
  PlanFingerprint thermophysical_boundary_semantics{};
  RevisionToken thermophysical_boundary_target{};
  PlanFingerprint thermophysical_boundary_rank_local_binding{};
  PlanFingerprint thermophysical_boundary_collective_lineage{};
  PlanFingerprint thermophysical_boundary_rank_local_lineage{};
  // Zero names ordinary C1/C2.  A positive ordinal is issued only after a
  // typed C2 provisional-refinement predecessor has been consumed.
  std::uint8_t pressure_energy_refinement{};
  PlanFingerprint pressure_energy_refinement_collective_lineage{};
  // Rank-local capability lineage; never compare this value across ranks.
  PlanFingerprint pressure_energy_refinement_lineage{};

  bool valid() const noexcept {
    return plan != 0U && r_au != 0U && h_by_a != 0U &&
           pressure_face_coefficient != 0U && phi_h_by_a != 0U &&
           trial_face_flux != 0U &&
           ((corrector == 1U && temporal_face_flux != 0U &&
             committed_face_history != 0U) ||
            (corrector == 2U && temporal_face_flux == 0U &&
             committed_face_history == 0U)) &&
           dependency != 0U &&
           (corrector == 1U || corrector == 2U) &&
           thermophysical_boundary_semantics != 0U &&
           thermophysical_boundary_target != 0U &&
           thermophysical_boundary_rank_local_binding != 0U &&
           thermophysical_boundary_collective_lineage != 0U &&
           thermophysical_boundary_rank_local_lineage != 0U &&
           ((pressure_energy_refinement == 0U &&
             pressure_energy_refinement_collective_lineage == 0U &&
             pressure_energy_refinement_lineage == 0U) ||
            (corrector == 2U && pressure_energy_refinement != 0U &&
             pressure_energy_refinement <=
                 kPressureEnergyRefinementCapacity &&
             pressure_energy_refinement_collective_lineage != 0U &&
             pressure_energy_refinement_lineage != 0U));
  }
};

struct PressureCorrectionInput {
  PisoIntermediateCertificate intermediate{};
  PressureReferenceCertificate pressure_reference{};
  ConstFieldView density_trial{};
  ConstFieldView density_accepted{};
  ConstFieldView density_previous{};
  ConstFieldView drho_dp_h_y{};
  BdfCoefficients bdf{};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken numeric_boundary{};
};

struct PressureCorrectionSystemView {
  FieldView diagonal{};
  FieldView rhs{};
};

struct PisoTrialStateView {
  FieldView velocity{};
  FieldView pressure_perturbation{};
  FieldView density{};
  ConstFieldView drho_dp_h_y{};
};

// Rank-local identity of one immutable field revision.  Exact-EOS closure
// authorities carry these identities so a same-numbered revision from a
// foreign storage replica cannot be substituted after thermodynamic
// evaluation.
struct PisoFieldRevisionIdentity {
  const double* base{};
  std::size_t replica{};
  FieldId field{};
  RevisionToken revision{};
  StorageIdentity storage{};
  RevisionDomainIdentity revision_domain{};

  bool valid() const noexcept {
    return base != nullptr && revision != 0U && storage != 0U &&
           revision_domain != 0U;
  }
};

inline PisoFieldRevisionIdentity make_piso_field_revision_identity(
    ConstFieldView field) noexcept {
  return {field.base, field.replica, field.field, field.revision,
          field.storage_identity, field.revision_domain};
}

// Caller-produced identity for an exact thermodynamic closure evaluated at
// p^{k+1}=p^k+dp and h^{k+1}=h^k+dh.  The PISO coupler cannot evaluate an EOS;
// it therefore verifies this complete provenance tuple before atomically
// publishing the candidate.  A terminal EOS oracle remains the independent
// numerical acceptance gate.
struct PisoExactEosClosureIdentity {
  PlanFingerprint thermodynamics{};
  RevisionToken pressure_reference{};
  PlanFingerprint composition{};
  PisoFieldRevisionIdentity pressure_state{};
  PisoFieldRevisionIdentity pressure_correction{};
  PisoFieldRevisionIdentity enthalpy_state{};
  PisoFieldRevisionIdentity enthalpy_correction{};
  PisoFieldRevisionIdentity candidate_enthalpy{};
  PisoFieldRevisionIdentity candidate_density{};
  PisoFieldRevisionIdentity candidate_temperature{};
  PlanFingerprint closure{};

  bool valid() const noexcept {
    return thermodynamics != 0U && pressure_reference != 0U &&
           composition != 0U && pressure_state.valid() &&
           pressure_correction.valid() && enthalpy_state.valid() &&
           enthalpy_correction.valid() && candidate_enthalpy.valid() &&
           candidate_density.valid() && candidate_temperature.valid() &&
           closure != 0U;
  }
};

struct PisoCoupledStateView {
  FieldView velocity{};
  FieldView pressure_perturbation{};
  FieldView enthalpy{};
  FieldView density{};
  FieldView temperature{};
};

struct PisoExactThermodynamicCandidateView {
  ConstFieldView enthalpy{};
  ConstFieldView density{};
  ConstFieldView temperature{};
  PisoExactEosClosureIdentity closure{};
  // Closed-mass exact corrections must carry the target-layer EOS
  // compressibility and the prepare-only gauge transaction that was derived
  // from the same raw pressure correction.  Open-boundary corrections leave
  // both members exactly empty.
  ConstFieldView pressure_compressibility{};
  ClosedGaugeCorrectionCertificate closed_gauge{};
  // Semantic target-layer composition aliases in the BoundaryPlan species
  // order.  Frozen candidate APIs carry their raw scratch composition in the
  // outer PisoFrozenMomentumExactCandidateInput and prove numeric equality
  // here; FieldId order is therefore canonical without making scratch FieldId
  // part of the physical provenance.  Callers must keep the span alive
  // through the synchronous certify/commit transaction.
  Span<const ConstFieldView> independent_species{};
};

enum class PisoStateClosure : std::uint8_t {
  pressure_affine,
  exact_eos,
};

struct PisoStateCorrectionCertificate {
  PlanFingerprint plan{};
  RevisionToken pressure_system{};
  RevisionToken correction{};
  RevisionToken enthalpy_correction{};
  RevisionToken velocity{};
  RevisionToken pressure{};
  RevisionToken enthalpy{};
  RevisionToken density{};
  RevisionToken temperature{};
  RevisionToken face_flux{};
  PlanFingerprint exact_eos_closure{};
  RevisionToken state{};
  std::uint8_t corrector{};
  PisoStateClosure closure{PisoStateClosure::pressure_affine};
  PlanFingerprint thermophysical_boundary_semantics{};
  RevisionToken thermophysical_boundary_target{};
  PlanFingerprint thermophysical_boundary_rank_local_binding{};
  PlanFingerprint thermophysical_boundary_collective_lineage{};
  PlanFingerprint thermophysical_boundary_rank_local_lineage{};
  PressureReferenceCertificate input_pressure_reference{};
  PressureReferenceCertificate output_pressure_reference{};
  PlanFingerprint closed_gauge_collective_transaction{};
  RevisionToken closed_gauge_rank_local_transaction{};
  // Exact-EOS candidates bind both the field-independent numeric composition
  // used by EOS/physical-boundary mass fluxes and the semantic target-layer
  // complete-view identity (FieldId order, base/storage/revision/domain and
  // strides).  The latter is deliberately rank-local and must not enter the
  // canonical physical lineage.  Raw candidate scratch never appears here.
  // Pressure-affine corrections leave all three members empty.
  PlanFingerprint composition_numeric_provenance{};
  PlanFingerprint composition_rank_local_binding{};
  std::size_t independent_species_count{};

  bool valid() const noexcept {
    const bool closure_valid =
        closure == PisoStateClosure::pressure_affine
            ? enthalpy_correction == 0U && enthalpy == 0U &&
                  temperature == 0U && exact_eos_closure == 0U
            : closure == PisoStateClosure::exact_eos &&
                  enthalpy_correction != 0U && enthalpy != 0U &&
                  temperature != 0U && exact_eos_closure != 0U;
    const bool reference_transition_valid =
        input_pressure_reference.valid() &&
        output_pressure_reference.valid() &&
        input_pressure_reference.plan == output_pressure_reference.plan &&
        input_pressure_reference.predictor ==
            output_pressure_reference.predictor &&
        input_pressure_reference.thermodynamics ==
            output_pressure_reference.thermodynamics &&
        input_pressure_reference.closure ==
            output_pressure_reference.closure &&
        input_pressure_reference.time == output_pressure_reference.time &&
        input_pressure_reference.kind == output_pressure_reference.kind &&
        (input_pressure_reference.kind == PressureReferenceKind::closed_mass
             ? closure == PisoStateClosure::exact_eos &&
                   closed_gauge_collective_transaction != 0U &&
                   closed_gauge_rank_local_transaction != 0U &&
                   output_pressure_reference.pressure_reference ==
                       closed_gauge_rank_local_transaction
             : input_pressure_reference.kind ==
                       PressureReferenceKind::boundary_absolute &&
                   input_pressure_reference.pressure_reference ==
                       output_pressure_reference.pressure_reference &&
                   closed_gauge_collective_transaction == 0U &&
                   closed_gauge_rank_local_transaction == 0U);
    const bool composition_valid =
        closure == PisoStateClosure::pressure_affine
            ? composition_numeric_provenance == 0U &&
                  composition_rank_local_binding == 0U &&
                  independent_species_count == 0U
            : composition_numeric_provenance != 0U &&
                  composition_rank_local_binding != 0U;
    return plan != 0U && pressure_system != 0U && correction != 0U &&
           velocity != 0U && pressure != 0U && density != 0U &&
           face_flux != 0U && state != 0U && closure_valid &&
           reference_transition_valid && composition_valid &&
           (corrector == 1U || corrector == 2U) &&
           thermophysical_boundary_semantics != 0U &&
           thermophysical_boundary_target != 0U &&
           thermophysical_boundary_rank_local_binding != 0U &&
           thermophysical_boundary_collective_lineage != 0U &&
           thermophysical_boundary_rank_local_lineage != 0U;
  }
};

// Single-use authority for continuing corrector two at the same target time.
// Only an exact, selected C2 candidate copied atomically into live trial
// p/h/rho/T/U and the provisional final-physical mass flux may issue it.
class PisoPressureEnergyRefinementStateCertificate {
 public:
  PisoPressureEnergyRefinementStateCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && predecessor_.valid() &&
           predecessor_.corrector == 2U && target_time_ != 0U &&
           iteration_ != 0U && iteration_ <= kPressureEnergyRefinementCapacity &&
           velocity_.valid() && pressure_.valid() && enthalpy_.valid() &&
           density_.valid() && temperature_.valid() && flux_.x.base != nullptr &&
           flux_.y.base != nullptr && flux_.z.base != nullptr &&
           flux_.revision != 0U && flux_.x.storage_identity != 0U &&
           flux_.x.revision_domain != 0U && state_provenance_ != 0U &&
           mass_flux_provenance_ != 0U && rank_local_state_numeric_ != 0U &&
           rank_local_flux_numeric_ != 0U && collective_lineage_ != 0U &&
           lineage_ != 0U;
  }
  std::uint8_t iteration() const noexcept { return iteration_; }
  PlanFingerprint collective_lineage() const noexcept {
    return collective_lineage_;
  }
  // Rank-local single-use capability lineage.
  PlanFingerprint lineage() const noexcept { return lineage_; }

 private:
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PisoStateCorrectionCertificate predecessor_{};
  PisoFieldRevisionIdentity velocity_{};
  PisoFieldRevisionIdentity pressure_{};
  PisoFieldRevisionIdentity enthalpy_{};
  PisoFieldRevisionIdentity density_{};
  PisoFieldRevisionIdentity temperature_{};
  ConstFaceFluxView flux_{};
  RevisionToken target_time_{};
  PlanFingerprint state_provenance_{};
  PlanFingerprint mass_flux_provenance_{};
  PlanFingerprint rank_local_state_numeric_{};
  PlanFingerprint rank_local_flux_numeric_{};
  PlanFingerprint collective_lineage_{};
  PlanFingerprint lineage_{};
  std::uint8_t iteration_{};
};

struct PisoTerminalAuditInput {
  PisoStateCorrectionCertificate correction{};
  PressureReferenceCertificate pressure_reference{};
  ConstFieldView density{};
  ConstFieldView eos_density{};
  ConstFieldView density_accepted{};
  ConstFieldView density_previous{};
  ConstFieldView pressure_perturbation{};
  ConstFieldView drho_dp_h_y{};
  BdfCoefficients bdf{};
  // Rank-local normalized energy residual. audit_pending_final takes the
  // collective maximum before applying an enabled energy_tolerance.
  double energy_residual{};
  double closed_mass_target{};
  double boundary_closure_residual{};
  std::uint64_t boundary_closure_samples{};
  // Empty means every local cell is active; immersed products pass the
  // topology's 0/1 fluid mask so solid control volumes do not enter norms.
  Span<const std::uint8_t> active{};
  BoundaryThermophysicalGhostUse thermophysical_boundary{};
};

struct PisoTerminalCertificate {
  PlanFingerprint plan{};
  RevisionToken correction_state{};
  RevisionToken final_flux{};
  RevisionToken audit_state{};
  PlanFingerprint thermophysical_boundary_semantics{};
  RevisionToken thermophysical_boundary_target{};
  PlanFingerprint thermophysical_boundary_rank_local_binding{};
  PlanFingerprint thermophysical_boundary_collective_lineage{};
  PlanFingerprint thermophysical_boundary_rank_local_lineage{};
  PressureReferenceCertificate pressure_reference{};
  PlanFingerprint closed_gauge_collective_transaction{};
  PlanFingerprint composition_numeric_provenance{};
  // Rank-local identity of the exact species storage audited at stage 60.
  // Numeric provenance alone cannot distinguish an equal-valued foreign
  // buffer with forged FieldId/revision metadata.
  PlanFingerprint composition_rank_local_binding{};
  std::size_t independent_species_count{};

  bool valid() const noexcept {
    return plan != 0U && correction_state != 0U && final_flux != 0U &&
           audit_state != 0U && thermophysical_boundary_semantics != 0U &&
           thermophysical_boundary_target != 0U &&
           thermophysical_boundary_rank_local_binding != 0U &&
           thermophysical_boundary_collective_lineage != 0U &&
           thermophysical_boundary_rank_local_lineage != 0U &&
           composition_numeric_provenance != 0U &&
           composition_rank_local_binding != 0U &&
           pressure_reference.valid() &&
           (pressure_reference.kind == PressureReferenceKind::closed_mass
                ? closed_gauge_collective_transaction != 0U
                : pressure_reference.kind ==
                          PressureReferenceKind::boundary_absolute &&
                      closed_gauge_collective_transaction == 0U);
  }
};

struct PressureCorrectionCertificate {
  PlanFingerprint plan{};
  RevisionToken intermediate{};
  RevisionToken time{};
  RevisionToken geometry{};
  RevisionToken numeric_boundary{};
  RevisionToken state{};
  std::uint8_t corrector{};
  PlanFingerprint thermophysical_boundary_semantics{};
  RevisionToken thermophysical_boundary_target{};
  PlanFingerprint thermophysical_boundary_rank_local_binding{};
  PlanFingerprint thermophysical_boundary_collective_lineage{};
  PlanFingerprint thermophysical_boundary_rank_local_lineage{};
  std::uint8_t pressure_energy_refinement{};
  PlanFingerprint pressure_energy_refinement_collective_lineage{};
  // Rank-local capability lineage; never compare this value across ranks.
  PlanFingerprint pressure_energy_refinement_lineage{};

  bool valid() const noexcept {
    return plan != 0U && intermediate != 0U && time != 0U &&
           geometry != 0U && numeric_boundary != 0U && state != 0U &&
           (corrector == 1U || corrector == 2U) &&
           thermophysical_boundary_semantics != 0U &&
           thermophysical_boundary_target != 0U &&
           thermophysical_boundary_rank_local_binding != 0U &&
           thermophysical_boundary_collective_lineage != 0U &&
           thermophysical_boundary_rank_local_lineage != 0U &&
           ((pressure_energy_refinement == 0U &&
             pressure_energy_refinement_collective_lineage == 0U &&
             pressure_energy_refinement_lineage == 0U) ||
            (corrector == 2U && pressure_energy_refinement != 0U &&
             pressure_energy_refinement <=
                 kPressureEnergyRefinementCapacity &&
             pressure_energy_refinement_collective_lineage != 0U &&
             pressure_energy_refinement_lineage != 0U));
  }
};

enum class PressureCorrectionFaceKind : std::uint8_t {
  invalid,
  exchanged,
  homogeneous_neumann,
  homogeneous_dirichlet,
  periodic
};

struct PressureCorrectionFaceRule {
  PressureCorrectionFaceKind kind{PressureCorrectionFaceKind::invalid};
  bool physical{};
  bool high{};
};

// The collective semantic identity is decomposition independent.  The
// rank-local layout identity additionally binds the concrete patch and its
// locally owned physical faces; neither token can stand in for the other.
struct PressureCorrectionBoundaryCertificate {
  PlanFingerprint semantic{};
  PlanFingerprint rank_local_layout{};
  RevisionToken geometry{};
  // Exact BoundaryPlan source revision used at compile time.  This is not
  // the attempt-local PressureCorrectionCertificate::numeric_boundary token.
  RevisionToken source_revision{};
  Int3 global_cells{};
  Int3 patch_begin{};
  Int3 local_cells{};
  PlanFingerprint geometry_fingerprint{};

  bool valid() const noexcept {
    return semantic != 0U && rank_local_layout != 0U && geometry != 0U &&
           source_revision != 0U && global_cells.x > 0 &&
           global_cells.y > 0 && global_cells.z > 0 && patch_begin.x >= 0 &&
           patch_begin.y >= 0 && patch_begin.z >= 0 && local_cells.x > 0 &&
           local_cells.y > 0 && local_cells.z > 0 &&
           patch_begin.x <= global_cells.x - local_cells.x &&
           patch_begin.y <= global_cells.y - local_cells.y &&
           patch_begin.z <= global_cells.z - local_cells.z &&
           geometry_fingerprint != 0U;
  }
};

// Immutable authority for every Cartesian pressure-correction boundary
// decision.  Compile performs no allocation.  Hot consumers bind once and
// reuse the same homogeneous Dirichlet/Neumann/periodic semantics for ghost
// values, MG, operator neighbours, and corrected face mass flux.
class PressureCorrectionBoundaryPlan {
 public:
  PressureCorrectionBoundaryPlan() noexcept = default;
  PressureCorrectionBoundaryPlan(const PressureCorrectionBoundaryPlan&) =
      delete;
  PressureCorrectionBoundaryPlan& operator=(
      const PressureCorrectionBoundaryPlan&) = delete;
  PressureCorrectionBoundaryPlan(PressureCorrectionBoundaryPlan&&) noexcept =
      default;
  PressureCorrectionBoundaryPlan& operator=(
      PressureCorrectionBoundaryPlan&&) noexcept = default;

  static Status compile(const CartesianGeometryPlan& geometry,
                        MeshPatch patch, const BoundaryPlan& boundary,
                        PressureCorrectionBoundaryPlan& out) noexcept;

  const PressureCorrectionBoundaryCertificate& certificate() const noexcept {
    return certificate_;
  }
  bool current() const noexcept;
  PressureCorrectionFaceKind kind(CartesianFace face) const noexcept;
  Status face_rule(CartesianAxis axis, Int3 face,
                   PressureCorrectionFaceRule& out) const noexcept;
  Status mg_boundaries(MgBoundarySet& out) const noexcept;
  Status fill_ghosts(FieldView correction) const noexcept;
  double neighbor_value(ConstFieldView correction, Int3 cell,
                        CartesianAxis axis, int direction) const noexcept;
  double jump(ConstFieldView correction, CartesianAxis axis,
              Int3 face) const noexcept;
  double mass_flux_response(ConstFieldView correction, CartesianAxis axis,
                            Int3 face, double pressure_face_coefficient,
                            double correction_scale = 1.0) const noexcept;

 private:
  friend class PressureVelocityCoupler;
  friend class PressureLinearOperator;
  friend class PressureEnergyPressureFluxOperator;

  Status face_rule_unchecked(CartesianAxis axis, Int3 face,
                             PressureCorrectionFaceRule& out) const noexcept;
  double neighbor_value_unchecked(ConstFieldView correction, Int3 cell,
                                  CartesianAxis axis,
                                  int direction) const noexcept;
  double jump_unchecked(ConstFieldView correction, CartesianAxis axis,
                        Int3 face) const noexcept;
  double mass_flux_response_unchecked(
      ConstFieldView correction, CartesianAxis axis, Int3 face,
      double pressure_face_coefficient,
      double correction_scale = 1.0) const noexcept;

  const CartesianGeometryPlan* geometry_{};
  const BoundaryPlan* boundary_{};
  MeshPatch patch_{};
  std::array<PressureCorrectionFaceKind, 6U> kinds_{};
  std::array<bool, 6U> local_physical_{};
  PlanFingerprint geometry_fingerprint_{};
  PlanFingerprint boundary_semantic_{};
  PlanFingerprint boundary_local_layout_{};
  PressureCorrectionBoundaryCertificate certificate_{};
};

struct PressureOperatorServices {
  MPI_Comm communicator{MPI_COMM_NULL};
  HaloEngine* halo{};
  StageId halo_stage{};
  FieldId solution_field{};
  std::uint8_t halo_width{1U};
};

struct PressureOperatorRevision {
  PressureCorrectionCertificate pressure{};
  LinearIdentity identity{};
  PlanFingerprint collective_fingerprint{};
};

class PressureLinearOperator final : public LinearOperator {
 public:
  PressureLinearOperator() noexcept = default;
  ~PressureLinearOperator() noexcept override;
  PressureLinearOperator(const PressureLinearOperator&) = delete;
  PressureLinearOperator& operator=(const PressureLinearOperator&) = delete;
  PressureLinearOperator(PressureLinearOperator&&) noexcept;
  PressureLinearOperator& operator=(PressureLinearOperator&&) noexcept;
  Status refresh(PressureOperatorRevision revision) noexcept;
  LinearOperatorCertificate certificate() const noexcept override;
  Status apply(FieldView x, FieldView y) const noexcept override;
  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override;
  PlanFingerprint fingerprint() const noexcept;
  std::uintptr_t coefficient_storage_address() const noexcept;

 private:
  friend class PressureVelocityCoupler;
  static Status bind_internal(
      const CartesianGeometryPlan& geometry, MeshPatch patch,
      const BoundaryPlan& boundary, PlanFingerprint coupler,
      ConstFaceFieldView x_coefficient,
      ConstFaceFieldView y_coefficient,
      ConstFaceFieldView z_coefficient,
      PressureOperatorServices services,
      PressureCorrectionSystemView system,
      PressureLinearOperator& out) noexcept;
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

// Single authority for the pressure-correction face jump used by both PISO
// flux publication and the pressure-to-energy Jacobian.  At a physical
// pressure-Dirichlet boundary it applies the homogeneous correction value;
// at a physical pressure-Neumann boundary the jump is exactly zero.  Periodic
// and inter-rank faces consume the already exchanged correction ghost.
double pressure_correction_jump(ConstFieldView correction,
                                const CartesianGeometryPlan& geometry,
                                MeshPatch patch, const BoundaryPlan& boundary,
                                CartesianAxis axis, Int3 face) noexcept;

// Linear mass-flux response d(phi_f)=-a_f jump(dp). correction_scale is used
// by the legacy line-search path and is one for a Jacobian action.
double pressure_correction_mass_flux_response(
    ConstFieldView correction, const CartesianGeometryPlan& geometry,
    MeshPatch patch, const BoundaryPlan& boundary, CartesianAxis axis,
    Int3 face, double pressure_face_coefficient,
    double correction_scale = 1.0) noexcept;

// h_f is not an arbitrary interpolation.  Its producer must freeze the exact
// face value used by the target-time enthalpy convection residual
// D(phi_f*h_f), including the convection scheme, flux direction, limiter
// active set, boundary closure, and IBM face activity.  reconstruction is the
// producer's collective semantic fingerprint for that complete choice.
// revision names the exact numeric target on this rank and may therefore be
// rank-local; it is never a collective contract token.
struct PressureEnergyFrozenFaceEnthalpy {
  ConstFaceFieldView x{};
  ConstFaceFieldView y{};
  ConstFaceFieldView z{};
  RevisionToken revision{};
  PlanFingerprint reconstruction{};
  // Rank-local binding of the three concrete face workspaces.  The producer
  // fills this with pressure_energy_frozen_face_enthalpy_local_binding(); the
  // operator recomputes it instead of trusting a transferable revision token.
  PlanFingerprint local_binding{};
};

PlanFingerprint pressure_energy_frozen_face_enthalpy_local_binding(
    const PressureEnergyFrozenFaceEnthalpy& frozen) noexcept;

class PressureVelocityCoupler;

// Coupler-issued authority for the Cartesian pressure-work derivative
//
//   -V U^k . G(dp) + V [rAU G(dp)] . G(pi^k).
//
// Callers can only transport this capability.  The issuer privately binds it
// to the exact current intermediate/pressure systems, the validated target
// pi/U views, and the coupler-owned HbyA/rAU workspaces.  An empty capability
// deliberately selects the legacy flux-only quasi-Newton action.
class PisoCartesianPressureWorkLinearization {
 public:
  PisoCartesianPressureWorkLinearization() noexcept = default;
  bool valid() const noexcept;

 private:
  friend class PressureVelocityCoupler;
  friend class PressureEnergyPressureFluxOperator;

  const PressureVelocityCoupler* issuer_{};
  const CartesianKernelPlan* kernels_{};
  const CartesianGeometryPlan* geometry_{};
  MeshPatch patch_{};
  PisoIntermediateCertificate intermediate_{};
  PressureCorrectionCertificate pressure_{};
  ConstFieldView target_pressure_perturbation_{};
  ConstFieldView target_velocity_{};
  ConstFieldView h_by_a_{};
  ConstFieldView r_au_{};
  RevisionToken authority_{};
  PlanFingerprint numeric_fingerprint_{};
};

enum class PressureEnergyPressureWorkScope : std::uint8_t {
  flux_only_quasi_newton,
  exact_cartesian,
};

struct PressureEnergyPressureFluxBinding {
  const CartesianGeometryPlan* geometry{};
  const BoundaryPlan* boundary{};
  MeshPatch patch{};
  PressureOperatorServices services{};
  PisoIntermediateCertificate intermediate{};
  PressureCorrectionCertificate pressure{};
  ConstFieldView temporal_diagonal{};
  ConstFaceFieldView x_pressure_coefficient{};
  ConstFaceFieldView y_pressure_coefficient{};
  ConstFaceFieldView z_pressure_coefficient{};
  // The exact target flux whose sign and limiter active set produced h_f.
  // Intermediate PISO flux deliberately carries no final-flux certificate.
  ConstFaceFluxView target_flux{};
  PressureEnergyFrozenFaceEnthalpy frozen_face_enthalpy{};
  PressureContinuityActivityView activity{};
  LinearIdentity identity{};
  // Appended after the legacy aggregate prefix so flux-only producers remain
  // source compatible while exact Cartesian callers opt in explicitly.
  PisoCartesianPressureWorkLinearization pressure_work{};
};

struct PressureEnergyPressureFluxCertificate {
  LinearOperatorCertificate linear{};
  RevisionToken binding_revision{};
  RevisionToken pressure_system{};
  RevisionToken intermediate_dependency{};
  RevisionToken pressure_face_coefficients{};
  RevisionToken target_flux_identity{};
  RevisionToken target_flux_revision{};
  RevisionToken temporal_diagonal{};
  RevisionToken frozen_face_enthalpy{};
  PlanFingerprint frozen_face_enthalpy_local_binding{};
  RevisionToken geometry{};
  RevisionToken numeric_boundary{};
  PlanFingerprint pressure_boundary_semantic{};
  PlanFingerprint pressure_boundary_rank_local_layout{};
  PlanFingerprint frozen_reconstruction{};
  PlanFingerprint activity_local_fingerprint{};
  PlanFingerprint activity_collective_fingerprint{};
  std::size_t active_cells{};
  std::size_t inactive_cells{};
  PressureEnergyPressureWorkScope pressure_work_scope{
      PressureEnergyPressureWorkScope::flux_only_quasi_newton};
  RevisionToken pressure_work_linearization{};
  bool conservative_face_response{};
  bool exact_piso_flux_jump{};
  bool frozen_target_enthalpy{};
  bool full_cartesian_pressure_work{};
  bool flux_only_quasi_newton{};
  bool allocation_free_apply{};
  PlanFingerprint thermophysical_boundary_semantics{};
  RevisionToken thermophysical_boundary_target{};
  PlanFingerprint thermophysical_boundary_rank_local_binding{};
  PlanFingerprint thermophysical_boundary_collective_lineage{};
  PlanFingerprint thermophysical_boundary_rank_local_lineage{};

  bool valid() const noexcept {
    return linear.identity.symbolic != 0U && linear.identity.numeric != 0U &&
           linear.identity.hierarchy != 0U && linear.identity.workspace != 0U &&
           linear.identity.fingerprint != 0U &&
           linear.collective_fingerprint != 0U && linear.local_shape.x > 0 &&
           linear.local_shape.y > 0 && linear.local_shape.z > 0 &&
           linear.operator_class == LinearOperatorClass::nonsymmetric &&
           binding_revision != 0U && pressure_system != 0U &&
           intermediate_dependency != 0U && pressure_face_coefficients != 0U &&
           target_flux_identity != 0U && target_flux_revision != 0U &&
           temporal_diagonal != 0U && frozen_face_enthalpy != 0U &&
           frozen_face_enthalpy_local_binding != 0U && geometry != 0U &&
           numeric_boundary != 0U && pressure_boundary_semantic != 0U &&
           pressure_boundary_rank_local_layout != 0U &&
           frozen_reconstruction != 0U &&
           active_cells + inactive_cells ==
               static_cast<std::size_t>(linear.local_shape.x) *
                   static_cast<std::size_t>(linear.local_shape.y) *
                   static_cast<std::size_t>(linear.local_shape.z) &&
           ((activity_local_fingerprint == 0U) ==
            (activity_collective_fingerprint == 0U)) &&
           ((pressure_work_scope ==
                     PressureEnergyPressureWorkScope::exact_cartesian &&
             pressure_work_linearization != 0U &&
             full_cartesian_pressure_work && !flux_only_quasi_newton &&
             activity_local_fingerprint == 0U &&
             activity_collective_fingerprint == 0U) ||
            (pressure_work_scope == PressureEnergyPressureWorkScope::
                                        flux_only_quasi_newton &&
             pressure_work_linearization == 0U &&
             !full_cartesian_pressure_work && flux_only_quasi_newton)) &&
           conservative_face_response && exact_piso_flux_jump &&
           frozen_target_enthalpy && allocation_free_apply &&
           thermophysical_boundary_semantics != 0U &&
           thermophysical_boundary_target != 0U &&
           thermophysical_boundary_rank_local_binding != 0U &&
           thermophysical_boundary_collective_lineage != 0U &&
           thermophysical_boundary_rank_local_lineage != 0U;
  }
};

// Certified same-target E_p action
//
//   E_p dp = diag(E_p^time) dp
//            + D(h_f [-a_f jump(dp)]).
//
// The persistent HaloEngine is reserved at bind time; apply performs no
// allocation.  Cell and face activity are the exact IBM continuity masks, so
// a solid row and every inactive interface contribute identically zero.
class PressureEnergyPressureFluxOperator final : public LinearOperator {
 public:
  PressureEnergyPressureFluxOperator() noexcept = default;
  ~PressureEnergyPressureFluxOperator() noexcept override = default;
  PressureEnergyPressureFluxOperator(
      const PressureEnergyPressureFluxOperator&) = delete;
  PressureEnergyPressureFluxOperator& operator=(
      const PressureEnergyPressureFluxOperator&) = delete;
  PressureEnergyPressureFluxOperator(
      PressureEnergyPressureFluxOperator&&) noexcept = default;
  PressureEnergyPressureFluxOperator& operator=(
      PressureEnergyPressureFluxOperator&&) noexcept = default;

  static Status bind(
      const PressureEnergyPressureFluxBinding& binding,
      PressureEnergyPressureFluxOperator& out,
      PressureEnergyPressureFluxCertificate& certificate) noexcept;

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_.linear;
  }
  LinearOperatorFailureProvenance failure_provenance() const noexcept override {
    return failure_;
  }
  const PressureEnergyPressureFluxCertificate& pressure_flux_certificate()
      const noexcept {
    return certificate_;
  }
  Status apply(FieldView input, FieldView output) const noexcept override;

 private:
  friend class PressureEnergySchurBlockAuthority;

  bool exact_pressure_work_current() const noexcept;

  PressureCorrectionBoundaryPlan pressure_boundary_{};
  MeshPatch patch_{};
  PressureOperatorServices services_{};
  ConstFieldView temporal_diagonal_{};
  ConstFaceFieldView x_pressure_coefficient_{};
  ConstFaceFieldView y_pressure_coefficient_{};
  ConstFaceFieldView z_pressure_coefficient_{};
  ConstFaceFluxView target_flux_{};
  PressureEnergyFrozenFaceEnthalpy frozen_face_enthalpy_{};
  PressureContinuityActivityView activity_{};
  PisoCartesianPressureWorkLinearization pressure_work_{};
  PressureEnergyPressureFluxCertificate certificate_{};
  mutable LinearOperatorFailureProvenance failure_{};
};

// Replayable same-target candidate evaluation keeps the converged momentum
// predictor frozen while refreshing the pressure-dependent velocity and mass
// flux.  The first production scope is deliberately Cartesian: immersed
// interfaces and physical fixed-flux faces must use their own compatible
// candidate closure before they can acquire this authority.
enum class PisoFrozenMomentumStageScope : std::uint8_t {
  invalid,
  cartesian_periodic,
  cartesian_open_boundary_ibm,
};

class PisoFrozenMomentumStageAuthority {
 public:
  PisoFrozenMomentumStageAuthority() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr &&
           (scope_ == PisoFrozenMomentumStageScope::cartesian_periodic ||
            scope_ ==
                PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm) &&
           intermediate_.valid() && pressure_.valid() &&
           pressure_.intermediate == intermediate_.dependency &&
           (corrector_ == 1U || corrector_ == 2U) &&
           corrector_ == intermediate_.corrector &&
           corrector_ == pressure_.corrector && baseline_ != 0U &&
           canonical_lineage_ != 0U;
  }
  PisoFrozenMomentumStageScope scope() const noexcept { return scope_; }
  std::uint8_t corrector() const noexcept { return corrector_; }
  PlanFingerprint canonical_lineage() const noexcept {
    return canonical_lineage_;
  }
 private:
  friend class PressureEnergyCandidateBoundaryFinalizer;
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PisoIntermediateCertificate intermediate_{};
  PressureCorrectionCertificate pressure_{};
  RevisionToken baseline_{};
  PlanFingerprint canonical_lineage_{};
  PisoFrozenMomentumStageScope scope_{PisoFrozenMomentumStageScope::invalid};
  std::uint8_t corrector_{};
};

// The scaled correction is intentionally a separate stage.  It writes owned
// cells only; callers must complete the registered correction halo before the
// velocity stage may consume its ghosts.
class PisoFrozenMomentumPressureStageCertificate {
 public:
  PisoFrozenMomentumPressureStageCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && stage_lineage_ != 0U &&
           pressure_direction_.valid() && scaled_pressure_correction_.valid() &&
           halo_instance_ != 0U &&
           std::isfinite(alpha_) && alpha_ >= 0.0 && alpha_ <= 1.0 &&
           correction_direction_ != 0U && canonical_lineage_ != 0U &&
           scratch_binding_ != 0U &&
           (corrector_ == 1U || corrector_ == 2U);
  }
  double alpha() const noexcept { return alpha_; }
  std::uint8_t corrector() const noexcept { return corrector_; }
  PlanFingerprint canonical_lineage() const noexcept {
    return canonical_lineage_;
  }
  PlanFingerprint scratch_binding() const noexcept { return scratch_binding_; }
  PlanFingerprint correction_direction() const noexcept {
    return correction_direction_;
  }
  RevisionToken scaled_pressure_revision() const noexcept {
    return scaled_pressure_correction_.revision;
  }

 private:
  friend class PressureEnergyCandidateBoundaryFinalizer;
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PlanFingerprint stage_lineage_{};
  PisoFieldRevisionIdentity pressure_direction_{};
  PisoFieldRevisionIdentity scaled_pressure_correction_{};
  ConstFieldView pressure_direction_view_{};
  ConstFieldView scaled_pressure_correction_view_{};
  std::uintptr_t halo_instance_{};
  double alpha_{};
  PlanFingerprint correction_direction_{};
  PlanFingerprint canonical_lineage_{};
  PlanFingerprint scratch_binding_{};
  std::uint8_t corrector_{};
};

class PisoFrozenMomentumVelocityStageCertificate {
 public:
  PisoFrozenMomentumVelocityStageCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && stage_lineage_ != 0U &&
           pressure_stage_lineage_ != 0U &&
           scaled_pressure_correction_.valid() && candidate_velocity_.valid() &&
           correction_halo_ != nullptr && halo_instance_ != 0U &&
           halo_ghost_revision_ != 0U &&
           candidate_velocity_numeric_ != 0U &&
           std::isfinite(alpha_) &&
           alpha_ >= 0.0 && alpha_ <= 1.0 && canonical_lineage_ != 0U &&
           scratch_binding_ != 0U && (corrector_ == 1U || corrector_ == 2U);
  }
  double alpha() const noexcept { return alpha_; }
  std::uint8_t corrector() const noexcept { return corrector_; }
  PlanFingerprint canonical_lineage() const noexcept {
    return canonical_lineage_;
  }
  PlanFingerprint scratch_binding() const noexcept { return scratch_binding_; }

 private:
  friend class PressureEnergyCandidateBoundaryFinalizer;
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PlanFingerprint stage_lineage_{};
  PlanFingerprint pressure_stage_lineage_{};
  PisoFieldRevisionIdentity scaled_pressure_correction_{};
  PisoFieldRevisionIdentity candidate_velocity_{};
  ConstFieldView scaled_pressure_correction_view_{};
  ConstFieldView candidate_velocity_view_{};
  const HaloEngine* correction_halo_{};
  std::uintptr_t halo_instance_{};
  RevisionToken halo_ghost_revision_{};
  PlanFingerprint candidate_donor_fingerprint_{};
  StageId candidate_donor_stage_{};
  FieldId candidate_donor_field_{};
  std::uint8_t candidate_donor_reach_{};
  PlanFingerprint ibm_geometry_fingerprint_{};
  PlanFingerprint candidate_velocity_numeric_{};
  double alpha_{};
  PlanFingerprint canonical_lineage_{};
  PlanFingerprint scratch_binding_{};
  std::uint8_t corrector_{};
};

class PisoFrozenMomentumFluxStageCertificate {
 public:
  PisoFrozenMomentumFluxStageCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && stage_lineage_ != 0U &&
           velocity_lineage_ != 0U && candidate_density_.valid() &&
           face_flux_revision_ != 0U && face_flux_storage_ != 0U &&
           face_flux_revision_domain_ != 0U && std::isfinite(alpha_) &&
           candidate_density_numeric_ != 0U &&
           mechanical_flux_numeric_ != 0U &&
           alpha_ >= 0.0 && alpha_ <= 1.0 && canonical_lineage_ != 0U &&
           scratch_binding_ != 0U &&
           (scope_ == PisoFrozenMomentumStageScope::cartesian_periodic ||
            (scope_ == PisoFrozenMomentumStageScope::
                           cartesian_open_boundary_ibm &&
             nonphysical_flux_provenance_ != 0U &&
             pressure_outlet_provisional_provenance_ != 0U &&
             ibm_interface_provenance_ != 0U)) &&
           (corrector_ == 1U || corrector_ == 2U);
  }
  double alpha() const noexcept { return alpha_; }
  std::uint8_t corrector() const noexcept { return corrector_; }
  PlanFingerprint canonical_lineage() const noexcept {
    return canonical_lineage_;
  }
  PlanFingerprint scratch_binding() const noexcept { return scratch_binding_; }
  PlanFingerprint nonphysical_flux_provenance() const noexcept {
    return nonphysical_flux_provenance_;
  }
  PlanFingerprint pressure_outlet_provisional_provenance() const noexcept {
    return pressure_outlet_provisional_provenance_;
  }
  PisoFrozenMomentumStageScope scope() const noexcept { return scope_; }

 private:
  friend class PressureEnergyCandidateBoundaryFinalizer;
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PlanFingerprint stage_lineage_{};
  PlanFingerprint velocity_lineage_{};
  PisoFieldRevisionIdentity candidate_density_{};
  ConstFieldView candidate_density_view_{};
  ConstFaceFluxView candidate_flux_view_{};
  RevisionToken face_flux_revision_{};
  StorageIdentity face_flux_storage_{};
  RevisionDomainIdentity face_flux_revision_domain_{};
  double alpha_{};
  PlanFingerprint nonphysical_flux_provenance_{};
  PlanFingerprint pressure_outlet_provisional_provenance_{};
  PlanFingerprint ibm_interface_provenance_{};
  PlanFingerprint ibm_geometry_fingerprint_{};
  PlanFingerprint candidate_density_numeric_{};
  PlanFingerprint mechanical_flux_numeric_{};
  PlanFingerprint canonical_lineage_{};
  PlanFingerprint scratch_binding_{};
  PisoFrozenMomentumStageScope scope_{
      PisoFrozenMomentumStageScope::invalid};
  std::uint8_t corrector_{};
};

// Opaque publication boundary between the frozen mechanical mapping and the
// candidate boundary finalizer.  Only the finalizer may issue production
// instances; the pressure/velocity coupler consumes and replays them but
// cannot self-sign a physical boundary closure.
class FinalBoundaryFluxCertificate {
 public:
  FinalBoundaryFluxCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && coupler_ != nullptr &&
           scope_ ==
               PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm &&
           stage_lineage_ != 0U && velocity_lineage_ != 0U &&
           mechanical_flux_lineage_ != 0U &&
           nonphysical_flux_provenance_ != 0U &&
           pressure_outlet_provisional_provenance_ != 0U &&
           boundary_semantic_ != 0U && boundary_layout_ != 0U &&
           pressure_reference_.valid() &&
           pressure_reference_.kind ==
               PressureReferenceKind::boundary_absolute &&
           std::isfinite(absolute_pressure_reference_) &&
           absolute_pressure_reference_ > 0.0 && target_time_ != 0U &&
           (corrector_ == 1U || corrector_ == 2U) &&
           std::isfinite(alpha_) && alpha_ >= 0.0 && alpha_ <= 1.0 &&
           candidate_pressure_.valid() && candidate_enthalpy_.valid() &&
           candidate_density_.valid() && candidate_temperature_.valid() &&
           candidate_velocity_.valid() && composition_identity_ != 0U &&
           composition_numeric_provenance_ != 0U &&
           independent_species_count_ == independent_species_views_.size &&
           (independent_species_count_ == 0U ||
            independent_species_views_.data != nullptr) &&
           candidate_state_provenance_ != 0U &&
           candidate_state_binding_ != 0U && state_halo_instance_ != 0U &&
           state_halo_lineage_ != 0U &&
           thermophysical_boundary_semantics_ != 0U &&
           thermophysical_boundary_target_ != 0U &&
           thermophysical_boundary_collective_lineage_ != 0U &&
           thermophysical_boundary_rank_local_lineage_ != 0U &&
           thermophysical_boundary_rank_local_binding_ != 0U &&
           physical_state_halo_lineage_ != 0U &&
           face_closure_lineage_ != 0U && inlet_flux_provenance_ != 0U &&
           outlet_backflow_provenance_ != 0U &&
           ibm_donor_lineage_ != 0U && ibm_geometry_lineage_ != 0U &&
           ibm_zero_interface_lineage_ != 0U &&
           std::isfinite(inlet_mass_target_) &&
           std::isfinite(inlet_mass_achieved_) &&
           final_flux_view_.x.base != nullptr &&
           final_flux_view_.y.base != nullptr &&
           final_flux_view_.z.base != nullptr && final_flux_revision_ != 0U &&
           final_flux_storage_ != 0U && final_flux_revision_domain_ != 0U &&
           final_flux_provenance_ != 0U && canonical_lineage_ != 0U &&
           scratch_binding_ != 0U;
  }
  PlanFingerprint canonical_lineage() const noexcept {
    return canonical_lineage_;
  }
  PlanFingerprint scratch_binding() const noexcept {
    return scratch_binding_;
  }
  PlanFingerprint ibm_donor_lineage() const noexcept {
    return ibm_donor_lineage_;
  }
  PlanFingerprint ibm_geometry_lineage() const noexcept {
    return ibm_geometry_lineage_;
  }
  PlanFingerprint ibm_zero_interface_lineage() const noexcept {
    return ibm_zero_interface_lineage_;
  }
  std::uint32_t outlet_fixed_point_iterations() const noexcept {
    return outlet_fixed_point_iterations_;
  }
  double inlet_mass_target() const noexcept { return inlet_mass_target_; }
  double inlet_mass_achieved() const noexcept {
    return inlet_mass_achieved_;
  }
  PlanFingerprint final_flux_provenance() const noexcept {
    return final_flux_provenance_;
  }

 private:
  friend class PressureEnergyCandidateBoundaryFinalizer;
  friend class PressureVelocityCoupler;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  friend class FinalBoundaryFluxCertificateTestAccess;
#endif
  const PressureEnergyCandidateBoundaryFinalizer* issuer_{};
  const PressureVelocityCoupler* coupler_{};
  PisoFrozenMomentumStageScope scope_{
      PisoFrozenMomentumStageScope::invalid};
  PlanFingerprint stage_lineage_{};
  PlanFingerprint velocity_lineage_{};
  PlanFingerprint mechanical_flux_lineage_{};
  PlanFingerprint nonphysical_flux_provenance_{};
  PlanFingerprint pressure_outlet_provisional_provenance_{};
  PlanFingerprint boundary_semantic_{};
  PlanFingerprint boundary_layout_{};
  PressureReferenceCertificate pressure_reference_{};
  double absolute_pressure_reference_{};
  RevisionToken target_time_{};
  double alpha_{};
  PisoFieldRevisionIdentity candidate_pressure_{};
  PisoFieldRevisionIdentity candidate_enthalpy_{};
  PisoFieldRevisionIdentity candidate_density_{};
  PisoFieldRevisionIdentity candidate_temperature_{};
  PisoFieldRevisionIdentity candidate_velocity_{};
  PlanFingerprint composition_identity_{};
  PlanFingerprint composition_numeric_provenance_{};
  Span<const ConstFieldView> independent_species_views_{};
  std::size_t independent_species_count_{};
  PlanFingerprint candidate_state_provenance_{};
  PlanFingerprint candidate_state_binding_{};
  std::uintptr_t state_halo_instance_{};
  PlanFingerprint state_halo_lineage_{};
  PlanFingerprint thermophysical_boundary_semantics_{};
  RevisionToken thermophysical_boundary_target_{};
  PlanFingerprint thermophysical_boundary_collective_lineage_{};
  PlanFingerprint thermophysical_boundary_rank_local_lineage_{};
  PlanFingerprint thermophysical_boundary_rank_local_binding_{};
  PlanFingerprint physical_state_halo_lineage_{};
  PlanFingerprint face_closure_lineage_{};
  PlanFingerprint inlet_flux_provenance_{};
  PlanFingerprint outlet_backflow_provenance_{};
  PlanFingerprint ibm_donor_lineage_{};
  PlanFingerprint ibm_geometry_lineage_{};
  PlanFingerprint ibm_zero_interface_lineage_{};
  double inlet_mass_target_{};
  double inlet_mass_achieved_{};
  std::uint32_t outlet_fixed_point_iterations_{};
  ConstFaceFluxView final_flux_view_{};
  RevisionToken final_flux_revision_{};
  StorageIdentity final_flux_storage_{};
  RevisionDomainIdentity final_flux_revision_domain_{};
  PlanFingerprint final_flux_provenance_{};
  PlanFingerprint canonical_lineage_{};
  PlanFingerprint scratch_binding_{};
  std::uint8_t corrector_{};
};

struct PressureEnergyCandidateBoundaryFinalizeInput {
  PisoFrozenMomentumStageAuthority authority{};
  PisoFrozenMomentumPressureStageCertificate pressure_stage{};
  PisoFrozenMomentumVelocityStageCertificate velocity_stage{};
  PisoFrozenMomentumFluxStageCertificate flux_stage{};
  PressureReferenceCertificate pressure_reference{};
  double absolute_pressure_reference{};
  ConstFieldView pressure_perturbation{};
  ConstFieldView enthalpy{};
  ConstFieldView density{};
  ConstFieldView temperature{};
  ConstFieldView velocity{};
  Span<const ConstFieldView> independent_species{};
  PlanFingerprint composition_identity{};
  BoundaryThermophysicalGhostUse thermophysical_boundary{};
  const HaloEngine* state_halo{};
  ConstFaceFluxView mechanical_flux{};
  FaceFluxView final_flux{};
};

// Certificate-free Fresh-start adapter for the same physical boundary-flux
// engine used by the hot pressure--enthalpy candidate finalizer.  It may
// close caller-owned scratch, but it cannot issue FinalBoundaryFluxCertificate
// publication authority.  The pressure and independent-species views carry
// the semantic BoundaryPlan FieldIds used to resolve configured T/Y tuples.
struct FreshPhysicalBoundaryFluxClosureInput {
  double absolute_pressure_reference{};
  ConstFieldView pressure_perturbation{};
  ConstFieldView velocity{};
  Span<const ConstFieldView> independent_species{};
  ConstFaceFluxView mechanical_flux{};
  FaceFluxView final_flux{};
};

struct PressureEnergyCandidateBoundaryFinalizerBinding {
  MPI_Comm communicator{MPI_COMM_NULL};
  const CartesianGeometryPlan* geometry{};
  MeshPatch patch{};
  const BoundaryPlan* boundary{};
  const CartesianKernelPlan* kernels{};
  const ThermodynamicsPlan* thermodynamics{};
  const TransportPlan* transport{};
  const PressureVelocityCoupler* coupler{};
  const IbmEquationInterfacePlan* immersed_interface{};
  const IbmPhysicalBoundaryFluxAuthority* immersed_physical_boundary_flux{};
  RemoteDonorExchangePlan* candidate_pressure_correction_donors{};
  StageId candidate_pressure_correction_donor_stage{};
  FieldId candidate_pressure_correction_field{};
  std::uint8_t candidate_pressure_correction_donor_reach{};
};

// Deep same-target physical-face module.  It has its own issuer identity:
// the coupler certifies only the frozen mechanical mapping and can consume,
// but can never manufacture, the final physical-boundary transaction.
class PressureEnergyCandidateBoundaryFinalizer {
 public:
  PressureEnergyCandidateBoundaryFinalizer() noexcept = default;
  ~PressureEnergyCandidateBoundaryFinalizer() noexcept;
  PressureEnergyCandidateBoundaryFinalizer(
      const PressureEnergyCandidateBoundaryFinalizer&) = delete;
  PressureEnergyCandidateBoundaryFinalizer& operator=(
      const PressureEnergyCandidateBoundaryFinalizer&) = delete;
  PressureEnergyCandidateBoundaryFinalizer(
      PressureEnergyCandidateBoundaryFinalizer&&) noexcept;
  PressureEnergyCandidateBoundaryFinalizer& operator=(
      PressureEnergyCandidateBoundaryFinalizer&&) noexcept;

  static Status bind(
      const PressureEnergyCandidateBoundaryFinalizerBinding& binding,
      PressureEnergyCandidateBoundaryFinalizer& out) noexcept;
  Status finalize(
      const PressureEnergyCandidateBoundaryFinalizeInput& input,
      ReductionEngine& reductions,
      FinalBoundaryFluxCertificate& certificate) noexcept;
  Status close_fresh_physical_flux(
      const FreshPhysicalBoundaryFluxClosureInput& input,
      ReductionEngine& reductions) noexcept;
  bool ready() const noexcept;
  PlanFingerprint fingerprint() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
class FinalBoundaryFluxCertificateTestAccess {
 public:
  static FinalBoundaryFluxCertificate make_foreign(
      const PressureVelocityCoupler* claimed_coupler,
      ConstFaceFluxView final_flux) noexcept;
};
#endif

struct PisoFrozenMomentumExactCandidateInput {
  ConstFieldView raw_enthalpy_direction{};
  ConstFieldView scaled_pressure_correction{};
  ConstFieldView scaled_enthalpy_correction{};
  PisoCoupledStateView base_state{};
  ConstFieldView candidate_pressure{};
  PisoExactThermodynamicCandidateView thermodynamic{};
  ConstFieldView candidate_velocity{};
  ConstFaceFluxView candidate_flux{};
  FinalBoundaryFluxCertificate final_boundary_flux{};
  Span<const ConstFieldView> independent_species{};
};

// Complete non-consuming binding of one staged pressure/enthalpy/EOS state.
// It is suitable for globalization sampling, but deliberately grants no
// publication authority without a matching selection certificate.
class PisoFrozenMomentumExactCandidateCertificate {
 public:
  PisoFrozenMomentumExactCandidateCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && stage_lineage_ != 0U &&
           pressure_lineage_ != 0U && velocity_lineage_ != 0U &&
           flux_lineage_ != 0U && pressure_stage_.valid() &&
           velocity_stage_.valid() && flux_stage_.valid() &&
           ((flux_stage_.scope() ==
                 PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm &&
             final_boundary_flux_.valid()) ||
            (flux_stage_.scope() ==
                 PisoFrozenMomentumStageScope::cartesian_periodic &&
             !final_boundary_flux_.valid())) &&
           exact_lineage_ != 0U &&
           scratch_binding_ != 0U &&
           correction_direction_ != 0U &&
           baseline_state_provenance_ != 0U &&
           baseline_mass_flux_provenance_ != 0U &&
           candidate_state_provenance_ != 0U &&
           candidate_mass_flux_provenance_ != 0U &&
           base_state_provenance_ != 0U &&
           composition_numeric_provenance_ != 0U &&
           composition_binding_ != 0U &&
           independent_species_count_ == independent_species_views_.size &&
           (independent_species_count_ == 0U ||
            independent_species_views_.data != nullptr) &&
           target_time_ != 0U &&
           std::isfinite(alpha_) && alpha_ >= 0.0 && alpha_ <= 1.0 &&
           (corrector_ == 1U || corrector_ == 2U);
  }
  double alpha() const noexcept { return alpha_; }
  std::uint8_t corrector() const noexcept { return corrector_; }
  RevisionToken target_time() const noexcept { return target_time_; }
  PlanFingerprint correction_direction() const noexcept {
    return correction_direction_;
  }
  PlanFingerprint baseline_state_provenance() const noexcept {
    return baseline_state_provenance_;
  }
  PlanFingerprint baseline_mass_flux_provenance() const noexcept {
    return baseline_mass_flux_provenance_;
  }
  PlanFingerprint candidate_state_provenance() const noexcept {
    return candidate_state_provenance_;
  }
  PlanFingerprint candidate_mass_flux_provenance() const noexcept {
    return candidate_mass_flux_provenance_;
  }
  PlanFingerprint base_state_provenance() const noexcept {
    return base_state_provenance_;
  }
  PlanFingerprint canonical_lineage() const noexcept { return exact_lineage_; }
  PlanFingerprint scratch_binding() const noexcept { return scratch_binding_; }
  bool final_boundary_flux_certified() const noexcept {
    return final_boundary_flux_.valid();
  }
  PlanFingerprint final_boundary_flux_provenance() const noexcept {
    return final_boundary_flux_.valid()
               ? final_boundary_flux_.final_flux_provenance()
               : PlanFingerprint{};
  }

 private:
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PlanFingerprint stage_lineage_{};
  PlanFingerprint pressure_lineage_{};
  PlanFingerprint velocity_lineage_{};
  PlanFingerprint flux_lineage_{};
  PisoFrozenMomentumPressureStageCertificate pressure_stage_{};
  PisoFrozenMomentumVelocityStageCertificate velocity_stage_{};
  PisoFrozenMomentumFluxStageCertificate flux_stage_{};
  FinalBoundaryFluxCertificate final_boundary_flux_{};
  PisoFieldRevisionIdentity base_velocity_{};
  PisoFieldRevisionIdentity base_pressure_{};
  PisoFieldRevisionIdentity base_enthalpy_{};
  PisoFieldRevisionIdentity base_density_{};
  PisoFieldRevisionIdentity base_temperature_{};
  PisoFieldRevisionIdentity raw_pressure_direction_{};
  PisoFieldRevisionIdentity scaled_pressure_correction_{};
  PisoFieldRevisionIdentity raw_enthalpy_direction_{};
  PisoFieldRevisionIdentity enthalpy_correction_{};
  PisoFieldRevisionIdentity candidate_pressure_{};
  PisoFieldRevisionIdentity candidate_velocity_{};
  ConstFieldView base_velocity_view_{};
  ConstFieldView base_pressure_view_{};
  ConstFieldView base_enthalpy_view_{};
  ConstFieldView base_density_view_{};
  ConstFieldView base_temperature_view_{};
  ConstFieldView raw_pressure_direction_view_{};
  ConstFieldView scaled_pressure_correction_view_{};
  ConstFieldView raw_enthalpy_direction_view_{};
  ConstFieldView enthalpy_correction_view_{};
  ConstFieldView candidate_pressure_view_{};
  ConstFieldView candidate_velocity_view_{};
  PisoExactThermodynamicCandidateView thermodynamic_candidate_{};
  ConstFaceFluxView candidate_flux_view_{};
  Span<const ConstFieldView> independent_species_views_{};
  double alpha_{};
  RevisionToken target_time_{};
  PlanFingerprint correction_direction_{};
  PlanFingerprint baseline_state_provenance_{};
  PlanFingerprint baseline_mass_flux_provenance_{};
  PlanFingerprint candidate_state_provenance_{};
  PlanFingerprint candidate_mass_flux_provenance_{};
  PlanFingerprint base_state_provenance_{};
  PlanFingerprint composition_numeric_provenance_{};
  PlanFingerprint composition_binding_{};
  std::size_t independent_species_count_{};
  PlanFingerprint exact_lineage_{};
  PlanFingerprint scratch_binding_{};
  std::uint8_t corrector_{};
};

// Separate publication authority for a converged alpha=0 baseline.  It is
// issued only after the exact baseline state is replayed, both scaled
// corrections are globally zero, and the compiled continuity/energy terminal
// gates are met.  It cannot be substituted for an Armijo selection.
class PressureEnergyStationaryCertificate {
 public:
  PressureEnergyStationaryCertificate() noexcept = default;
  bool valid() const noexcept {
    return issuer_ != nullptr && exact_lineage_ != 0U && target_time_ != 0U &&
           correction_direction_ != 0U && state_provenance_ != 0U &&
           mass_flux_provenance_ != 0U && stationary_lineage_ != 0U &&
           (corrector_ == 1U || corrector_ == 2U) &&
           maximum_pressure_correction_ == 0.0 &&
           maximum_enthalpy_correction_ == 0.0 &&
           std::isfinite(normalized_continuity_) &&
           normalized_continuity_ >= 0.0 &&
           std::isfinite(normalized_energy_) && normalized_energy_ >= 0.0 &&
           std::isfinite(continuity_limit_) && continuity_limit_ > 0.0 &&
           std::isfinite(energy_limit_) && energy_limit_ > 0.0 &&
           normalized_continuity_ <= continuity_limit_ &&
           normalized_energy_ <= energy_limit_;
  }
  std::uint8_t corrector() const noexcept { return corrector_; }
  RevisionToken target_time() const noexcept { return target_time_; }
  PlanFingerprint canonical_lineage() const noexcept {
    return stationary_lineage_;
  }

 private:
  friend class PressureVelocityCoupler;
  const PressureVelocityCoupler* issuer_{};
  PlanFingerprint exact_lineage_{};
  RevisionToken target_time_{};
  PlanFingerprint correction_direction_{};
  PlanFingerprint state_provenance_{};
  PlanFingerprint mass_flux_provenance_{};
  double maximum_pressure_correction_{};
  double maximum_enthalpy_correction_{};
  double normalized_continuity_{};
  double normalized_energy_{};
  double continuity_limit_{};
  double energy_limit_{};
  PlanFingerprint stationary_lineage_{};
  std::uint8_t corrector_{};
};

class PressureVelocityCoupler {
 public:
  PressureVelocityCoupler() noexcept = default;
  ~PressureVelocityCoupler() noexcept;
  PressureVelocityCoupler(const PressureVelocityCoupler&) = delete;
  PressureVelocityCoupler& operator=(const PressureVelocityCoupler&) = delete;
  PressureVelocityCoupler(PressureVelocityCoupler&&) noexcept;
  PressureVelocityCoupler& operator=(PressureVelocityCoupler&&) noexcept;
  static Status bind(const PisoPlan& plan,
                     const EquationPlanSet& equations,
                     PisoCouplerServices services,
                     PisoCouplerWorkspace workspace,
                     PressureVelocityCoupler& out) noexcept;
  Status refresh(const PisoIntermediateInput& input,
                 PisoIntermediateCertificate& certificate,
                 Status prerequisite = {}) noexcept;
  Status refresh_pressure_energy_refinement(
      const PisoPressureEnergyRefinementStateCertificate& refinement,
      const PisoIntermediateInput& input, PisoCoupledStateView state,
      PisoIntermediateCertificate& certificate,
      Status prerequisite = {}) noexcept;
  Status inspect_intermediate_flux(
      const PisoIntermediateCertificate& intermediate,
      ConstFaceFluxView& flux) const noexcept;
  Status inspect_cartesian_pressure_work_linearization(
      const PisoIntermediateCertificate& intermediate,
      const PressureCorrectionCertificate& pressure,
      ConstFieldView target_pressure_perturbation,
      ConstFieldView target_velocity,
      PisoCartesianPressureWorkLinearization& linearization) const noexcept;
  Status assemble_pressure_system(
      const PressureCorrectionInput& input,
      PressureCorrectionSystemView system,
      PressureCorrectionCertificate& certificate) noexcept;
  Status make_frozen_momentum_stage_authority(
      const PisoIntermediateCertificate& intermediate,
      const PressureCorrectionCertificate& pressure,
      PisoFrozenMomentumStageAuthority& authority) const noexcept;
  Status form_frozen_momentum_scaled_pressure(
      const PisoFrozenMomentumStageAuthority& authority,
      ConstFieldView pressure_direction,
      const HaloEngine& candidate_correction_halo,
      double alpha, FieldView& scaled_pressure_correction,
      PisoFrozenMomentumPressureStageCertificate& certificate) const noexcept;
  Status stage_frozen_momentum_velocity(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
      const HaloEngine& candidate_correction_halo,
      ConstFieldView scaled_pressure_correction, FieldView candidate_velocity,
      PisoFrozenMomentumVelocityStageCertificate& certificate) const noexcept;
  Status stage_frozen_momentum_velocity(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
      const HaloEngine& candidate_correction_halo,
      FieldView& scaled_pressure_correction, FieldView candidate_velocity,
      PisoFrozenMomentumVelocityStageCertificate& certificate) const noexcept;
  Status stage_frozen_momentum_flux(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumVelocityStageCertificate& velocity,
      ConstFieldView candidate_density, FaceFluxView candidate_flux,
      PisoFrozenMomentumFluxStageCertificate& certificate) const noexcept;
  Status certify_frozen_momentum_exact_baseline(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
      const PisoFrozenMomentumVelocityStageCertificate& velocity_stage,
      const PisoFrozenMomentumFluxStageCertificate& flux_stage,
      PisoFrozenMomentumExactCandidateInput input,
      ReductionEngine& reductions,
      PisoFrozenMomentumExactCandidateCertificate& certificate) const noexcept;
  Status certify_frozen_momentum_exact_candidate(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& baseline,
      const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
      const PisoFrozenMomentumVelocityStageCertificate& velocity_stage,
      const PisoFrozenMomentumFluxStageCertificate& flux_stage,
      PisoFrozenMomentumExactCandidateInput input, ReductionEngine& reductions,
      PisoFrozenMomentumExactCandidateCertificate& certificate) const noexcept;
  Status certify_frozen_momentum_stationary(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& baseline,
      double global_normalized_continuity,
      double global_normalized_energy, ReductionEngine& reductions,
      PressureEnergyStationaryCertificate& certificate) const noexcept;
  Status commit_frozen_momentum_coupled_trial_state(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& candidate,
      const PressureEnergyGlobalizationSelectionCertificate& selection,
      PisoCoupledStateView state, FaceFluxView trial_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status commit_frozen_momentum_coupled_pending_state(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& candidate,
      const PressureEnergyGlobalizationSelectionCertificate& selection,
      PisoCoupledStateView state, PendingFaceFluxView& pending_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status commit_frozen_momentum_coupled_refinement_trial_state(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& candidate,
      const PressureEnergyGlobalizationSelectionCertificate& selection,
      PisoCoupledStateView state, FaceFluxView trial_flux,
      std::uint8_t refinement_iteration, ReductionEngine& reductions,
      PisoStateCorrectionCertificate& correction,
      PisoPressureEnergyRefinementStateCertificate& refinement) noexcept;
  Status commit_frozen_momentum_stationary_trial_state(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& baseline,
      const PressureEnergyStationaryCertificate& stationary,
      PisoCoupledStateView state, FaceFluxView trial_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status commit_frozen_momentum_stationary_pending_state(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& baseline,
      const PressureEnergyStationaryCertificate& stationary,
      PisoCoupledStateView state, PendingFaceFluxView& pending_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status bind_pressure_operator(
      PressureOperatorServices services,
      PressureCorrectionSystemView system,
      PressureLinearOperator& out) const noexcept;
  Status make_native_pressure_mg_spec(
      MPI_Comm communicator, LinearIdentity identity,
      MgCoefficientIdentity coefficients,
      NativeCartesianMgSpec& out) const noexcept;
  Status compile_native_pressure_mg(
      const PressureCorrectionCertificate& pressure,
      const NativeCartesianMgSpec& spec,
      MgRuntimeServices services,
      PressureCorrectionSystemView system,
      NativeCartesianMgPlan& out,
      MgPlanCounters* counters = nullptr) const noexcept;
  Status refresh_pressure_linear_lifecycle(
      const PressureCorrectionCertificate& pressure,
      LinearIdentity identity,
      MgCoefficientIdentity coefficients,
      PressureCorrectionSystemView system,
      PressureLinearOperator& linear_operator,
      NativeCartesianMgPlan& preconditioner,
      MgPlanCounters* counters = nullptr,
      bool collective_fail_close = false) const noexcept;
  Status correct_trial_state(
      const PressureCorrectionCertificate& pressure,
      FieldView correction,
      PisoTrialStateView state,
      FaceFluxView trial_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status correct_coupled_trial_state(
      const PressureCorrectionCertificate& pressure,
      FieldView correction,
      ConstFieldView enthalpy_correction,
      PisoCoupledStateView state,
      PisoExactThermodynamicCandidateView candidate,
      FaceFluxView trial_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status correct_pending_state(
      const PressureCorrectionCertificate& pressure,
      FieldView correction,
      PisoTrialStateView state,
      PendingFaceFluxView& pending_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status correct_coupled_pending_state(
      const PressureCorrectionCertificate& pressure,
      FieldView correction,
      ConstFieldView enthalpy_correction,
      PisoCoupledStateView state,
      PisoExactThermodynamicCandidateView candidate,
      PendingFaceFluxView& pending_flux,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  // Read the exact C2-corrected pending mass flux before terminal audit so
  // callers can assemble final_conservative residuals. This grants no final
  // publication authority and returns no PisoTerminalCertificate.
  Status inspect_corrected_pending(
      const PisoStateCorrectionCertificate& correction,
      const PendingFaceFluxView& pending_flux,
      ConstFaceFluxView& view) const noexcept;
  Status audit_pending_final(
      const PisoTerminalAuditInput& input,
      const PendingFaceFluxView& pending_flux,
      ReductionEngine& reductions,
      PisoAttemptReport& report,
      PisoTerminalCertificate& certificate,
      Status prerequisite = {}) noexcept;
  Status inspect_pending_final(
      const PisoTerminalCertificate& terminal,
      const PendingFaceFluxView& pending_flux,
      ConstFaceFluxView& view) const noexcept;
  Status publish_pending_final(
      const PisoTerminalCertificate& terminal,
      Span<const RevisionDependency> dependencies,
      Span<const ConstFieldView> independent_species,
      ReductionEngine& reductions,
      FinalFaceFluxWriter& writer,
      PendingFaceFluxView& pending_flux) noexcept;
  PlanFingerprint fingerprint() const noexcept;
  std::uintptr_t workspace_storage_address() const noexcept;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  void set_frozen_stationary_tolerances_for_test(
      double continuity_tolerance, double energy_tolerance) noexcept;
#endif

 private:
  friend class PressureEnergyCandidateBoundaryFinalizer;
  friend class PressureEnergyPressureFluxOperator;
  enum class StateCorrectionContract : std::uint8_t {
    pressure_unsealed,
    pressure_sealed,
  };
  friend class PisoPlan;
  friend class PisoPressureSolveEpoch;
  friend class PressureContinuityConvergenceAudit;
  struct Impl;
  LinearConvergenceAuditCertificate pressure_convergence_certificate(
      const PressureCorrectionCertificate& pressure) const noexcept;
  Status audit_pressure_convergence(
      const PressureCorrectionCertificate& pressure,
      ConstFieldView correction,
      ConstFieldView true_residual,
      ReductionEngine& reductions,
      double& selected_alpha,
      LinearConvergenceAuditResult& result) noexcept;
  Status capture_pressure_failure_provenance(
      const PressureCorrectionCertificate& pressure,
      ConstFieldView correction,
      ConstFieldView sealed_rhs,
      ReductionEngine& reductions,
      LinearConvergenceAuditResult& result) noexcept;
  Status stage_frozen_momentum_velocity_impl(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
      const HaloEngine& candidate_correction_halo,
      ConstFieldView scaled_pressure_correction,
      FieldView* mutable_scaled_pressure_correction,
      FieldView candidate_velocity,
      PisoFrozenMomentumVelocityStageCertificate& certificate) const noexcept;
  Status correct_state_impl(
      const PressureCorrectionCertificate& pressure,
      FieldView correction,
      PisoTrialStateView state,
      FaceFluxView output_flux,
      std::uint8_t corrector,
      StateCorrectionContract contract,
      double sealed_alpha,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status correct_exact_coupled_state_impl(
      const PressureCorrectionCertificate& pressure,
      FieldView pressure_correction,
      ConstFieldView enthalpy_correction,
      PisoCoupledStateView state,
      PisoExactThermodynamicCandidateView candidate,
      FaceFluxView output_flux,
      std::uint8_t corrector,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status commit_frozen_momentum_exact_state_impl(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate& candidate,
      const PressureEnergyGlobalizationSelectionCertificate* selection,
      const PressureEnergyStationaryCertificate* stationary,
      PisoCoupledStateView state,
      FaceFluxView output_flux, std::uint8_t corrector,
      ReductionEngine& reductions,
      PisoStateCorrectionCertificate& certificate) noexcept;
  Status certify_frozen_momentum_exact_candidate_impl(
      const PisoFrozenMomentumStageAuthority& authority,
      const PisoFrozenMomentumExactCandidateCertificate* baseline,
      const PisoFrozenMomentumPressureStageCertificate& pressure_stage,
      const PisoFrozenMomentumVelocityStageCertificate& velocity_stage,
      const PisoFrozenMomentumFluxStageCertificate& flux_stage,
      PisoFrozenMomentumExactCandidateInput input, ReductionEngine& reductions,
      PisoFrozenMomentumExactCandidateCertificate& certificate) const noexcept;
  Status refresh_impl(
      const PisoPressureEnergyRefinementStateCertificate* refinement,
      const PisoIntermediateInput& input, const PisoCoupledStateView* state,
      PisoIntermediateCertificate& certificate,
      Status prerequisite) noexcept;
  bool matches_cartesian_pressure_work_linearization(
      const PisoCartesianPressureWorkLinearization& linearization,
      bool verify_numeric_fingerprint) const noexcept;
  bool matches_candidate_boundary_finalizer_binding(
      const CartesianGeometryPlan* geometry, MeshPatch patch,
      const BoundaryPlan* boundary, const CartesianKernelPlan* kernels,
      const ThermodynamicsPlan* thermodynamics,
      const TransportPlan* transport,
      const IbmEquationInterfacePlan* immersed_interface,
      const RemoteDonorExchangePlan* candidate_donors,
      StageId candidate_donor_stage, FieldId candidate_field,
      std::uint8_t candidate_reach) const noexcept;
  void release() noexcept;
  Impl* implementation_{};
};

enum class PisoPressureSolveContract : std::uint8_t {
  pressure_continuity,
  continuity_energy_coupled,
};

class PisoPressureSolveEpoch {
 public:
  PisoPressureSolveEpoch() noexcept = default;
  ~PisoPressureSolveEpoch() noexcept;
  PisoPressureSolveEpoch(const PisoPressureSolveEpoch&) = delete;
  PisoPressureSolveEpoch& operator=(const PisoPressureSolveEpoch&) = delete;
  PisoPressureSolveEpoch(PisoPressureSolveEpoch&&) = delete;
  PisoPressureSolveEpoch& operator=(PisoPressureSolveEpoch&&) = delete;
  Status begin(const PisoPlan& plan) noexcept;
  Status prepare_linear_lifecycle(
      const PisoPlan& plan, std::uint8_t corrector,
      const PressureCorrectionCertificate& pressure,
      LinearIdentity identity,
      MgCoefficientIdentity coefficients,
      PressureVelocityCoupler& coupler,
      PressureLinearOperator& lifecycle_operator,
      NativeCartesianMgPlan& preconditioner,
      PressureCorrectionSystemView system,
      FieldView correction,
      SolverWorkspace& workspace,
      MgPlanCounters* mg_counters = nullptr) noexcept;
  Status prepare_pressure_energy_refinement_lifecycle(
      const PisoPlan& plan, std::uint8_t refinement_iteration,
      const PressureCorrectionCertificate& pressure,
      LinearIdentity identity, MgCoefficientIdentity coefficients,
      PressureVelocityCoupler& coupler,
      PressureLinearOperator& lifecycle_operator,
      NativeCartesianMgPlan& preconditioner,
      PressureCorrectionSystemView system, FieldView correction,
      SolverWorkspace& workspace,
      MgPlanCounters* mg_counters = nullptr) noexcept;
  Status solve_prepared(
      LinearOperator& exact_operator,
      PisoPressureSolveContract contract,
      ReductionEngine& reductions,
      ResourceCounters* resources = nullptr) noexcept;
  Status solve(
      const PisoPlan& plan, std::uint8_t corrector,
      const PressureCorrectionCertificate& pressure,
      LinearIdentity identity,
      MgCoefficientIdentity coefficients,
      PressureVelocityCoupler& coupler,
      PressureLinearOperator& linear_operator,
      NativeCartesianMgPlan& preconditioner,
      PressureCorrectionSystemView system,
      FieldView correction,
      SolverWorkspace& workspace,
      ReductionEngine& reductions,
      ResourceCounters* resources = nullptr,
      MgPlanCounters* mg_counters = nullptr) noexcept;
  Status solve(
      const PisoPlan& plan, std::uint8_t corrector,
      const PressureCorrectionCertificate& pressure,
      LinearIdentity identity,
      MgCoefficientIdentity coefficients,
      PressureVelocityCoupler& coupler,
      PressureLinearOperator& lifecycle_operator,
      LinearOperator& exact_operator,
      NativeCartesianMgPlan& preconditioner,
      PressureCorrectionSystemView system,
      FieldView correction,
      SolverWorkspace& workspace,
      ReductionEngine& reductions,
      ResourceCounters* resources = nullptr,
      MgPlanCounters* mg_counters = nullptr) noexcept;
  Status finalize(PisoAttemptReport& report) noexcept;
  Status observe(PisoAttemptReport& report) const noexcept;
  std::uint8_t solve_calls() const noexcept { return solve_calls_; }
  std::uint8_t refinement_solve_calls() const noexcept {
    return refinement_solve_calls_;
  }
  bool active() const noexcept { return active_; }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  static bool validate_pressure_energy_refinement_report_for_test(
      const PisoAttemptReport& report) noexcept;
#endif

 private:
  struct PreparedLinearLifecycle {
    const PisoPlan* plan{};
    PressureCorrectionCertificate pressure{};
    LinearIdentity identity{};
    PressureVelocityCoupler* coupler{};
    PressureLinearOperator* lifecycle_operator{};
    NativeCartesianMgPlan* preconditioner{};
    PressureCorrectionSystemView system{};
    FieldView correction{};
    SolverWorkspace* workspace{};
    LinearOperatorCertificate lifecycle_certificate{};
    LinearPreconditionerCertificate preconditioner_certificate{};
    std::uint8_t corrector{};
    std::uint8_t refinement_iteration{};
    bool pressure_energy_refinement{};
    bool valid{};
  };
  void discard_workspace() noexcept;
  std::array<LinearSolveResult, 2U> results_{};
  std::array<PisoPressureEnergyRefinementSolveReport,
             kPressureEnergyRefinementCapacity>
      refinement_results_{};
  PreparedLinearLifecycle prepared_{};
  PressureVelocityCoupler* epoch_coupler_{};
  SolverWorkspace* epoch_workspace_{};
  SolverWorkspace* workspace_{};
  MPI_Comm epoch_communicator_{MPI_COMM_NULL};
  int epoch_rank_{-1};
  int epoch_size_{};
  PlanFingerprint plan_{};
  RevisionToken refinement_target_generation_{};
  std::uint8_t solve_calls_{};
  std::uint8_t refinement_solve_calls_{};
  bool active_{};
  bool failed_{};
};

Status assemble_continuity(
    const ContinuityEquationPlan& plan, const EquationStateView& state,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept;

Status assemble_momentum(
    const MomentumEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept;

Status assemble_momentum_predictor(
    const MomentumEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    FieldView low_order_rhs_delta,
    EquationAssemblyCertificate& certificate) noexcept;

Status limit_momentum_predictor_correction(
    ConstFieldView velocity, ConstFaceFluxView mass_flux,
    MgDomainActivityView activity, EquationSystemView system,
    ConstFieldView low_order_rhs_delta, ReductionEngine& reductions,
    MomentumPredictorLimiterReport& report) noexcept;

Status solve_momentum_predictor(
    MPI_Comm communicator, const MomentumEquationPlan& plan,
    const BoundaryPlan& boundary,
    MeshPatch patch, const EquationAssemblyCertificate& assembly,
    ConstFaceFluxView mass_flux, MgDomainActivityView activity,
    EquationSystemView system, FieldView velocity, HaloEngine& krylov_halo,
    SolverWorkspace& workspace, ReductionEngine& reductions,
    ResourceCounters* resources,
    MomentumPredictorSolveReport& report) noexcept;

Status assemble_enthalpy(
    const EnthalpyEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept;

Status assemble_species(
    const SpeciesEquationPlan& plan, std::size_t species,
    const EquationStateView& state, const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept;

Status assemble_scalar(
    const ScalarEquationPlan& plan, std::size_t scalar,
    const EquationStateView& state, const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept;

Status assemble_pressure_storage(
    const PressureReferencePlan& plan, ConstFieldView drho_dp_hY,
    const PressureReferenceCertificate& pressure_reference,
    const EquationAssemblyContext& context, FieldView cell_integral_diagonal,
    EquationAssemblyCertificate& certificate) noexcept;

Status evaluate_thermophysical_rates(
    const EquationPlanSet& plans, const ThermophysicalRateInput& input,
    ThermophysicalRateOutput output,
    ThermophysicalRateCertificate& certificate) noexcept;

Status form_scalar_mass_diffusivity(
    const ScalarEquationSpec& spec, ConstFieldView molecular_viscosity,
    ConstFieldView turbulent_viscosity, KernelBox box,
    FieldView mass_diffusivity) noexcept;

struct PressureWorkPoint {
  double pressure{};
  double accepted_pressure{};
  double previous_pressure{};
  Real3 pressure_gradient{};
  Real3 velocity{};
  double velocity_divergence{};
};

Status evaluate_pressure_material_derivative(
    BdfCoefficients bdf, const PressureWorkPoint& point,
    double& material_derivative) noexcept;

enum class EnthalpyTermKind : std::uint8_t {
  unsteady,
  advection,
  pressure_work,
  conduction,
  viscous_dissipation,
  source
};

struct EnthalpyTermPoint {
  EnthalpyTermKind kind{EnthalpyTermKind::unsteady};
  double value{};
  double auxiliary0{};
  double auxiliary1{};
  double auxiliary2{};
};

Status assemble_isolated_enthalpy_term(const EnthalpyTermPoint& term,
                                       double cell_volume,
                                       double& residual) noexcept;
Status newtonian_viscous_dissipation(const VelocityGradient& gradient,
                                     double effective_viscosity,
                                     double& dissipation) noexcept;

struct SpeciesClosure {
  double dependent_mass_fraction{};
  double dependent_diffusive_flux{};
};

Status close_independent_species(Span<const double> independent,
                                 Span<const double> diffusive_fluxes,
                                 SpeciesClosure& closure) noexcept;

Status resolve_heat_flux_normal_gradient(double outward_heat_flux,
                                         double conductivity,
                                         double& normal_gradient) noexcept;
Status resolve_scalar_flux_normal_gradient(double outward_scalar_flux,
                                           double mass_diffusivity,
                                           double& normal_gradient) noexcept;

}  // namespace hundun::v04
