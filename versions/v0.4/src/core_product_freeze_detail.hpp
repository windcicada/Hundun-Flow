// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_product.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04::detail {

enum ProductCapability : std::uint32_t {
  product_fields = 1U << 0U,
  product_structured_ghosts = 1U << 1U,
  product_ibm_donors = 1U << 2U,
  product_cache_dependencies = 1U << 3U,
  product_exact_numeric = 1U << 4U,
  product_coarse_numeric = 1U << 5U,
  product_preconditioner = 1U << 6U,
  product_workspace = 1U << 7U,
  product_service_snapshot = 1U << 8U,
  product_collective_epochs = 1U << 9U
};

inline constexpr std::uint32_t kProductRequiredCapabilities =
    product_fields | product_structured_ghosts | product_cache_dependencies |
    product_exact_numeric | product_coarse_numeric | product_preconditioner |
    product_workspace | product_service_snapshot | product_collective_epochs;

// A closed-mass gauge change may alter the final rounded bit of an otherwise
// identical energy-residual replay.  This is a numerical reproducibility
// allowance, not a physical residual tolerance: an O(EOS-tolerance) mismatch
// must never be relabelled as gauge equivalence.
inline constexpr double kAlphaZeroEnergyReplayRoundoffTolerance =
    256.0 * std::numeric_limits<double>::epsilon();

inline bool alpha_zero_energy_replay_equivalent(
    bool bitwise_equal, bool closed_mass_gauge_shifted,
    double relative_error) noexcept {
  return bitwise_equal ||
         (closed_mass_gauge_shifted && std::isfinite(relative_error) &&
          relative_error <= kAlphaZeroEnergyReplayRoundoffTolerance);
}

inline Status validate_product_capabilities(std::uint32_t capabilities,
                                            bool immersed) noexcept {
  const std::uint32_t required =
      kProductRequiredCapabilities | (immersed ? product_ibm_donors : 0U);
  return (capabilities & required) == required
             ? Status{}
             : Status{StatusCode::invalid_plan, 10203U};
}

struct ProductPressureEnergyTemporalState {
  double density{};
  double enthalpy{};
  double absolute_pressure{};
};

// Backward-error scale for the temporal part of
//   d(rho h)/dt - dp/dt.
//
// The BDF identity cancels six independently rounded operands at order two.
// Scaling by the already-cancelled target a0 term alone makes a legitimate
// small-dt history cancellation look like an O(1) energy defect.  Retain the
// magnitude of every rho*h and pressure operand instead.  BE deliberately
// ignores the previous state so an unbound n-1 layer cannot affect recovery.
inline bool product_pressure_energy_temporal_operand_scale(
    BdfCoefficients bdf, double volume,
    ProductPressureEnergyTemporalState target,
    ProductPressureEnergyTemporalState accepted,
    ProductPressureEnergyTemporalState previous,
    double& scale) noexcept {
  scale = 0.0;
  const bool bdf2 = bdf.order == 2U;
  const auto valid_state = [](ProductPressureEnergyTemporalState state) {
    return std::isfinite(state.density) && state.density > 0.0 &&
           std::isfinite(state.enthalpy) &&
           std::isfinite(state.absolute_pressure) &&
           state.absolute_pressure > 0.0;
  };
  if ((bdf.order != 1U && !bdf2) || !std::isfinite(volume) ||
      !(volume > 0.0) || !std::isfinite(bdf.a0) || !(bdf.a0 > 0.0) ||
      !std::isfinite(bdf.a1) || !std::isfinite(bdf.a2) ||
      !valid_state(target) || !valid_state(accepted) ||
      (bdf2 && !valid_state(previous)))
    return false;

  const auto operand = [](double coefficient,
                          ProductPressureEnergyTemporalState state,
                          double& magnitude) {
    const double rho_h = state.density * state.enthalpy;
    magnitude = std::abs(coefficient) *
                (std::abs(rho_h) + std::abs(state.absolute_pressure));
    return std::isfinite(rho_h) && std::isfinite(magnitude);
  };
  double target_magnitude = 0.0;
  double accepted_magnitude = 0.0;
  double previous_magnitude = 0.0;
  if (!operand(bdf.a0, target, target_magnitude) ||
      !operand(bdf.a1, accepted, accepted_magnitude) ||
      (bdf2 && !operand(bdf.a2, previous, previous_magnitude)))
    return false;
  const double temporal_magnitude =
      target_magnitude + accepted_magnitude + previous_magnitude;
  const double candidate =
      volume * std::max(std::abs(bdf.a0), temporal_magnitude);
  if (!std::isfinite(temporal_magnitude) ||
      !std::isfinite(candidate) || !(candidate > 0.0))
    return false;
  scale = candidate;
  return true;
}

struct ProductPressureInexactForcingState {
  double normalized_continuity{};
  double normalized_energy{};
  double previous_merit{};
  double terminal_tolerance{};
  bool residual_available{};
  bool previous_merit_available{};
};

inline bool product_pressure_coupled_merit(
    double normalized_continuity, double normalized_energy,
    double& merit) noexcept {
  merit = 0.0;
  if (!std::isfinite(normalized_continuity) ||
      normalized_continuity < 0.0 ||
      !std::isfinite(normalized_energy) || normalized_energy < 0.0)
    return false;
  const double candidate =
      std::hypot(normalized_continuity, normalized_energy);
  if (!std::isfinite(candidate)) return false;
  merit = candidate;
  return true;
}

inline bool product_use_simple_diagonal_schur(
    CouplingKind coupling, std::uint8_t corrector,
    std::uint8_t refinement_iteration, bool immersed) noexcept {
  // Reserve half of the existing refinement budget for the spatial response
  // if the cheap diagonal approximation has not reached the physical gates.
  // This is still the existing IBM quasi-Newton operator, not a claim of an
  // exact IBM Jacobian, and does not extend the iteration/acceptance limits.
  return coupling == CouplingKind::simple && corrector == 2U && immersed &&
         refinement_iteration <= kPressureEnergyRefinementCapacity / 2U;
}

// Select only the inner coupled linear tolerance.  The case control remains
// the lower bound and is restored near the terminal gate, on stagnation, or
// whenever the nonlinear evidence is invalid.
inline LinearSolveControl product_pressure_inexact_forcing_control(
    LinearSolveControl base,
    ProductPressureInexactForcingState state) noexcept {
  const bool valid_base = std::isfinite(base.absolute_tolerance) &&
                          std::isfinite(base.relative_tolerance) &&
                          base.absolute_tolerance > 0.0 &&
                          base.relative_tolerance > 0.0 &&
                          base.relative_tolerance < 1.0 &&
                          base.maximum_iterations > 0U &&
                          base.true_residual_interval > 0U;
  if (!valid_base ||
      base.relative_tolerance >=
          kPressureInexactForcingRelativeToleranceCeiling ||
      !std::isfinite(state.terminal_tolerance) ||
      !(state.terminal_tolerance > 0.0) ||
      base.relative_tolerance < state.terminal_tolerance)
    return base;

  if (!state.residual_available) {
    base.relative_tolerance =
        kPressureInexactForcingRelativeToleranceCeiling;
    return base;
  }
  double merit = 0.0;
  if (!product_pressure_coupled_merit(state.normalized_continuity,
                                      state.normalized_energy, merit))
    return base;
  const double terminal_component =
      std::max(state.normalized_continuity, state.normalized_energy);
  const double terminal_band = 10.0 * state.terminal_tolerance;
  if (!std::isfinite(terminal_band) || terminal_component <= terminal_band)
    return base;

  if (state.previous_merit_available) {
    if (!std::isfinite(state.previous_merit) ||
        !(state.previous_merit > 0.0) ||
        merit >= 0.9 * state.previous_merit)
      return base;
  }
  const double requested = 0.1 * merit;
  if (!std::isfinite(requested)) return base;
  base.relative_tolerance = std::min(
      kPressureInexactForcingRelativeToleranceCeiling,
      std::max(base.relative_tolerance, requested));
  return base;
}

// Scalar Aitken relaxation for a contracting fixed-point merit.  The observed
// contraction already includes the previous relaxation, so carry that factor
// forward instead of lagging one accelerated update behind.  The exact
// candidate layer remains the acceptance authority; invalid or stagnating
// evidence keeps the legacy full step, and extrapolation is capped at two.
inline double product_pressure_aitken_initial_alpha(
    double previous_merit, double current_merit,
    double previous_alpha) noexcept {
  if (!std::isfinite(previous_merit) || !(previous_merit > 0.0) ||
      !std::isfinite(current_merit) || !(current_merit > 0.0) ||
      !std::isfinite(previous_alpha) || previous_alpha < 1.0 ||
      previous_alpha > kPressureEnergyAitkenMaximumAlpha ||
      current_merit >= 0.9 * previous_merit)
    return 1.0;
  const double contraction = current_merit / previous_merit;
  if (!std::isfinite(contraction) || !(contraction > 0.0) ||
      !(contraction < 0.9))
    return 1.0;
  const double alpha = previous_alpha / (1.0 - contraction);
  if (!std::isfinite(alpha) || !(alpha > 1.0)) return 1.0;
  return std::min(kPressureEnergyAitkenMaximumAlpha, alpha);
}

#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
inline constexpr std::size_t
    kPressureEnergyCandidateGlobalizationSampleCapacity = 25U;
inline constexpr std::size_t kPressureEnergyCandidateFieldCount = 5U;
inline constexpr std::size_t kPressureEnergyCandidateScratchFieldCount = 9U;
inline constexpr std::size_t kPressureEnergyCandidateFluxReplicaCount = 4U;

struct PressureEnergyCandidateFieldStorageDiagnostic {
  FieldId candidate_field{};
  FieldId trial_field{};
  std::uintptr_t candidate_base{};
  std::uintptr_t candidate_begin{};
  std::uintptr_t candidate_end{};
  std::uintptr_t trial_base{};
  std::uintptr_t trial_begin{};
  std::uintptr_t trial_end{};
  RevisionToken candidate_revision{};
  RevisionToken trial_revision{};
  StorageIdentity candidate_storage{};
  StorageIdentity trial_storage{};
  RevisionDomainIdentity candidate_revision_domain{};
  RevisionDomainIdentity trial_revision_domain{};
  std::uint8_t components{};
  std::uint8_t ghost_width{};
};

struct PressureEnergyCandidateFluxStorageDiagnostic {
  std::array<std::uintptr_t, 3U> bases{};
  std::uintptr_t replica_begin{};
  std::uintptr_t replica_end{};
  StorageIdentity storage{};
  RevisionDomainIdentity revision_domain{};
  RevisionToken revision{};
};

struct PressureEnergyCandidateStorageDiagnostic {
  bool valid{};
  PlanFingerprint plan{};
  PlanFingerprint lineage_fingerprint{};
  std::array<PressureEnergyCandidateFieldStorageDiagnostic,
             kPressureEnergyCandidateFieldCount>
      fields{};
  std::array<PressureEnergyCandidateFieldStorageDiagnostic,
             kPressureEnergyCandidateScratchFieldCount>
      scratch_fields{};
  std::uintptr_t coupled_state_halo{};
  std::uintptr_t candidate_state_halo{};
  std::uintptr_t coupled_thermal_halo{};
  std::uintptr_t candidate_thermal_halo{};
  std::uintptr_t candidate_finalizer_state_halo{};
  std::uintptr_t correction_halo{};
  std::uintptr_t candidate_correction_halo{};
  HaloPlanStats coupled_state_halo_plan{};
  HaloPlanStats candidate_state_halo_plan{};
  HaloPlanStats coupled_thermal_halo_plan{};
  HaloPlanStats candidate_thermal_halo_plan{};
  HaloPlanStats candidate_finalizer_state_halo_plan{};
  HaloPlanStats correction_halo_plan{};
  HaloPlanStats candidate_correction_halo_plan{};
  PlanFingerprint candidate_pressure_donor_plan{};
  PlanFingerprint candidate_velocity_donor_plan{};
  PlanFingerprint candidate_rate_donor_plan{};
  std::uint8_t local_ibm_reconstruction_reach{};
  std::uint8_t candidate_pressure_donor_reach{};
  RemoteDonorExchangeStats candidate_pressure_donor_stats{};
  RemoteDonorExchangeStats candidate_velocity_donor_stats{};
  RemoteDonorExchangeStats candidate_rate_donor_stats{};
  PlanFingerprint execution_graph{};
  StageResourceSpec corrector_one_resource_contract{};
  StageResourceSpec corrector_two_resource_contract{};
  std::array<PressureEnergyCandidateFluxStorageDiagnostic,
             kPressureEnergyCandidateFluxReplicaCount>
      flux_replicas{};
  FaceFluxStorageCounters workspace_flux_capacity{};
  FaceFluxStorageCounters final_flux_capacity{};
  ExecutionCounters field_storage_capacity{};
  std::size_t arena_doubles{};
};

enum class PressureEnergyCandidateFailureReason : std::uint8_t {
  none,
  nonfinite_pressure,
  nonfinite_enthalpy,
  nonpositive_absolute_pressure,
  inactive_nonzero_correction,
  inactive_density,
  inactive_temperature,
  thermodynamic_evaluation,
  density,
  temperature,
  heat_capacity,
  pressure_compressibility,
  enthalpy_compressibility,
  transport_evaluation,
  viscosity,
  conductivity,
  production_candidate_evaluation
};

struct PressureEnergyCandidateAlphaDiagnostic {
  double alpha{};
  bool admissible{};
  std::uint64_t first_failing_global_cell{
      std::numeric_limits<std::uint64_t>::max()};
  PressureEnergyCandidateFailureReason first_failure_reason{
      PressureEnergyCandidateFailureReason::none};
  double minimum_absolute_pressure{std::numeric_limits<double>::infinity()};
  double minimum_temperature{std::numeric_limits<double>::infinity()};
  double minimum_density{std::numeric_limits<double>::infinity()};
  double normalized_continuity{std::numeric_limits<double>::infinity()};
  double normalized_energy{std::numeric_limits<double>::infinity()};
  bool state_and_flux_finite{};
};

struct PressureEnergyCandidateGlobalizationDiagnostic {
  bool valid{};
  bool production_candidate_loop{};
  bool baseline_commit{};
  bool selection_valid{};
  bool replay_valid{};
  bool committed{};
  std::uint64_t attempted_step{};
  std::uint64_t generation{};
  std::uint32_t attempt{};
  std::uint8_t corrector{};
  std::uint8_t sample_count{};
  std::uint8_t first_admissible_sample{
      std::numeric_limits<std::uint8_t>::max()};
  double maximum_absolute_pressure_correction{};
  double maximum_absolute_enthalpy_correction{};
  double alpha_zero_maximum_velocity_difference{};
  double alpha_zero_maximum_density_difference{};
  double alpha_zero_maximum_temperature_difference{};
  double alpha_zero_maximum_pressure_difference{};
  double alpha_zero_maximum_pressure_relative_error{};
  double alpha_zero_maximum_density_relative_error{};
  double alpha_zero_maximum_temperature_relative_error{};
  std::uint64_t alpha_zero_maximum_density_ulp_difference{};
  std::uint64_t alpha_zero_maximum_temperature_ulp_difference{};
  double alpha_zero_material_oracle_error{};
  double alpha_zero_gradient_oracle_error{};
  double alpha_zero_energy_residual_oracle_error{};
  double alpha_zero_effective_viscosity_oracle_error{};
  PlanFingerprint alpha_zero_oracle_numeric_lineage{};
  bool alpha_zero_velocity_bitwise_equal{};
  bool alpha_zero_density_bitwise_equal{};
  bool alpha_zero_temperature_bitwise_equal{};
  double baseline_normalized_continuity{};
  double baseline_normalized_energy{};
  double selected_normalized_continuity{};
  double selected_normalized_energy{};
  double selected_alpha{};
  double corrector_one_baseline_normalized_continuity{
      std::numeric_limits<double>::infinity()};
  double corrector_one_baseline_normalized_energy{
      std::numeric_limits<double>::infinity()};
  double corrector_one_selected_normalized_continuity{
      std::numeric_limits<double>::infinity()};
  double corrector_one_selected_normalized_energy{
      std::numeric_limits<double>::infinity()};
  double corrector_one_selected_alpha{};
  double corrector_one_linear_predicted_normalized_continuity{
      std::numeric_limits<double>::infinity()};
  double corrector_one_linear_predicted_normalized_energy{
      std::numeric_limits<double>::infinity()};
  double corrector_two_linear_predicted_normalized_continuity{
      std::numeric_limits<double>::infinity()};
  double corrector_two_linear_predicted_normalized_energy{
      std::numeric_limits<double>::infinity()};
  PlanFingerprint correction_direction{};
  PlanFingerprint selected_state_provenance{};
  PlanFingerprint selected_mass_flux_provenance{};
  bool final_boundary_flux_certified{};
  PlanFingerprint final_boundary_canonical_lineage{};
  PlanFingerprint final_physical_flux_provenance{};
  double final_inlet_mass_target{};
  double final_inlet_mass_achieved{};
  std::uint64_t local_ibm_interface_links{};
  std::uint64_t local_ibm_interface_nonzero_count{};
  std::uint64_t local_ibm_interface_negative_zero_count{};
  std::uint64_t candidate_runtime_halo_messages{};
  std::uint64_t candidate_runtime_halo_bytes{};
  std::uint64_t candidate_sealed_halo_messages{};
  std::uint64_t candidate_sealed_halo_bytes{};
  std::array<PressureEnergyCandidateAlphaDiagnostic,
             kPressureEnergyCandidateGlobalizationSampleCapacity>
      samples{};
};

enum class PressureEnergyCandidatePoisonKind : std::uint8_t {
  density_without_revision,
};

struct ProductFinalFluxHistoryDiagnostic {
  bool valid{};
  RevisionToken accepted_revision{};
  RevisionToken previous_revision{};
  FaceFluxCertificate accepted_certificate{};
  FaceFluxCertificate previous_certificate{};
  PlanFingerprint accepted_payload{};
  PlanFingerprint previous_payload{};
};

struct PressureCorrectionWarmStartDiagnostic {
  bool valid{};
  StepOrigin origin{StepOrigin::fresh_start};
  std::uint32_t attempt{};
  bool authority_available{};
  bool used{};
};

enum class FreshInitializationPoisonKind : std::uint8_t {
  velocity,
  flux,
};

enum class RestartRestoreAllocationPoint : std::uint8_t {
  thermophysical_staging,
  scalar_seed,
};

// Test-only observation of the production Fresh initialization transaction.
// It intentionally publishes numeric summaries and typed lineages rather than
// exposing mutable state storage; no production decision consumes this data.
struct FreshInitializationDiagnostic {
  bool valid{};
  PlanFingerprint driver_plan{};
  RevisionToken generation{};
  bool immersed{};
  bool projection_attempted{};
  bool no_ibm_bypassed{};
  bool audited{};
  bool committed{};
  Status terminal_status{};
  LinearSolveResult solve{};
  PlanFingerprint red_plan{};
  PlanFingerprint prepared_lineage{};
  PlanFingerprint solved_lineage{};
  PlanFingerprint candidate_lineage{};
  double initial_continuity_maximum{};
  double final_continuity_maximum{};
  double maximum_face_envelope{};
  double continuity_limit{};
  std::uint64_t cut_face_nonzero_count{};
  std::uint64_t cut_face_negative_zero_count{};
  std::uint64_t solid_nonpositive_zero_component_count{};
  bool velocity_layers_bitwise_equal{};
  bool derived_velocity_dependents_rebuilt{};
  PlanFingerprint derived_velocity_lineage{};
  double maximum_h_by_a_velocity_difference{};
  std::uint64_t velocity_gradient_nonfinite_count{};
  std::uint64_t effective_viscosity_nonpositive_count{};
  std::uint64_t solid_velocity_gradient_nonfinite_count{};
  bool velocity_dependent_rate_layers_bitwise_equal{};
};

struct ColdVelocityDependentsDiagnostic {
  bool valid{};
  bool restart{};
  bool rebuilt{};
  PlanFingerprint driver_plan{};
  PlanFingerprint lineage{};
  double maximum_h_by_a_velocity_difference{};
  std::uint64_t velocity_gradient_nonfinite_count{};
  std::uint64_t effective_viscosity_nonpositive_count{};
  std::uint64_t solid_velocity_gradient_nonfinite_count{};
  bool rate_layers_bitwise_equal{};
};

// Arm on every rank in the communicator before ProductDriver::advance.  The
// next successfully recovered pressure--enthalpy correction is inspected
// synchronously and the arm is consumed before any diagnostic collective.
// The observer never scales, clips, or writes the live candidate state.
void arm_pressure_energy_candidate_globalization_once_for_test() noexcept;
void clear_pressure_energy_candidate_globalization_for_test() noexcept;
bool pressure_energy_candidate_globalization_diagnostic_for_test(
    PressureEnergyCandidateGlobalizationDiagnostic& out) noexcept;
bool pressure_energy_candidate_storage_diagnostic_for_test(
    PressureEnergyCandidateStorageDiagnostic& out) noexcept;
// Inject one rank-local PISO cold-bind contract failure. ProductCompiler must
// globalize it before any other rank enters the candidate-finalizer collective.
void arm_product_piso_bind_failure_once_for_test(int failing_rank) noexcept;
void clear_product_piso_bind_failure_for_test() noexcept;
void arm_pressure_energy_candidate_poison_once_for_test(
    int poison_rank, PressureEnergyCandidatePoisonKind kind,
    std::uint8_t corrector = 0U) noexcept;
void clear_pressure_energy_candidate_poison_for_test() noexcept;
bool product_final_flux_history_diagnostic_for_test(
    ProductFinalFluxHistoryDiagnostic& out) noexcept;
bool pressure_correction_warm_start_diagnostic_for_test(
    PressureCorrectionWarmStartDiagnostic& out) noexcept;
// Suppress only the next otherwise-eligible C1 warm seed.  This permits a
// same-state cold comparison without changing the time proposal or physics.
void suppress_pressure_correction_warm_start_once_for_test() noexcept;
void arm_fresh_initialization_candidate_poison_once_for_test(
    int poison_rank, FreshInitializationPoisonKind kind) noexcept;
void clear_fresh_initialization_diagnostic_for_test() noexcept;
bool fresh_initialization_diagnostic_for_test(
    FreshInitializationDiagnostic& out) noexcept;
bool cold_velocity_dependents_diagnostic_for_test(
    ColdVelocityDependentsDiagnostic& out) noexcept;
void arm_restart_restore_allocation_failure_once_for_test(
    int failing_rank, RestartRestoreAllocationPoint point) noexcept;
void clear_restart_restore_allocation_failure_for_test() noexcept;
int restart_restore_allocation_lowest_failing_rank_for_test() noexcept;
Status convective_cfl_acceptance_status_for_test(
    TimeControlKind control, double outward_max, double limit) noexcept;

inline Status validate_product_capabilities_for_test(
    std::uint32_t capabilities, bool immersed) noexcept {
  return validate_product_capabilities(capabilities, immersed);
}
#endif

inline bool product_cell_count(Int3 cells, std::size_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) return false;
  const std::size_t x = static_cast<std::size_t>(cells.x);
  const std::size_t y = static_cast<std::size_t>(cells.y);
  const std::size_t z = static_cast<std::size_t>(cells.z);
  if (x > std::numeric_limits<std::size_t>::max() / y ||
      x * y > std::numeric_limits<std::size_t>::max() / z) {
    return false;
  }
  out = x * y * z;
  return true;
}

inline bool product_checked_add(std::size_t left, std::size_t right,
                                std::size_t& out) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left) return false;
  out = left + right;
  return true;
}

inline bool product_checked_multiply(std::size_t left, std::size_t right,
                                     std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  out = left * right;
  return true;
}

inline bool product_field_bytes(std::size_t cells, std::size_t components,
                                std::size_t& out) noexcept {
  std::size_t doubles = 0U;
  return product_checked_multiply(cells, components, doubles) &&
         product_checked_multiply(doubles, sizeof(double), out);
}

inline bool product_face_doubles(Int3 cells, std::size_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0 ||
      cells.x == std::numeric_limits<std::int32_t>::max() ||
      cells.y == std::numeric_limits<std::int32_t>::max() ||
      cells.z == std::numeric_limits<std::int32_t>::max())
    return false;
  std::size_t x_faces = 0U;
  std::size_t y_faces = 0U;
  std::size_t z_faces = 0U;
  std::size_t plane = 0U;
  const std::size_t x = static_cast<std::size_t>(cells.x);
  const std::size_t y = static_cast<std::size_t>(cells.y);
  const std::size_t z = static_cast<std::size_t>(cells.z);
  if (!product_checked_multiply(x + 1U, y, plane) ||
      !product_checked_multiply(plane, z, x_faces) ||
      !product_checked_multiply(x, y + 1U, plane) ||
      !product_checked_multiply(plane, z, y_faces) ||
      !product_checked_multiply(x, y, plane) ||
      !product_checked_multiply(plane, z + 1U, z_faces) ||
      !product_checked_add(x_faces, y_faces, plane) ||
      !product_checked_add(plane, z_faces, out))
    return false;
  return true;
}

inline bool product_halo_bytes(Int3 cells, std::size_t components,
                               std::uint8_t width,
                               std::size_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0 || components == 0U ||
      width == 0U)
    return false;
  const std::size_t x = static_cast<std::size_t>(cells.x);
  const std::size_t y = static_cast<std::size_t>(cells.y);
  const std::size_t z = static_cast<std::size_t>(cells.z);
  std::size_t yz = 0U;
  std::size_t xz = 0U;
  std::size_t xy = 0U;
  std::size_t faces = 0U;
  std::size_t values = 0U;
  if (!product_checked_multiply(y, z, yz) ||
      !product_checked_multiply(x, z, xz) ||
      !product_checked_multiply(x, y, xy) ||
      !product_checked_add(yz, xz, faces) ||
      !product_checked_add(faces, xy, faces) ||
      !product_checked_multiply(faces, 2U, faces) ||
      !product_checked_multiply(faces, width, values) ||
      !product_checked_multiply(values, components, values) ||
      !product_checked_multiply(values, sizeof(double), out))
    return false;
  return true;
}

inline std::uint64_t product_mix(std::uint64_t hash,
                                 std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
  return hash;
}

}  // namespace hundun::v04::detail
