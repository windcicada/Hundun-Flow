// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/flow_immersed.hpp"
#include "flow_immersed_access_detail.hpp"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace hundun::finite_volume {
class ImmersedOperatorAdapter;
}

namespace hundun::flow::test {

enum class ImmersedFlowAttemptFailureStage : std::uint8_t {
  none,
  after_corrector_1,
  after_corrector_2,
  final_momentum_residual,
  after_final_transport,
  final_wall_penetration,
  final_continuity_residual,
  final_pressure_residual,
  final_force_reconstruction,
  before_commit
};

enum class ImmersedFlowWallInputFailure : std::uint8_t {
  none,
  non_positive_density,
  non_finite_density,
  non_positive_coefficient,
  non_finite_coefficient,
  stale_preconditioner_revision
};

using ImmersedFlowWallMeasureSnapshot = detail::ImmersedFlowWallMeasure;
using ImmersedFlowWallGradientSnapshot = detail::ImmersedFlowWallGradient;
using ImmersedFlowExactPredictorWorkspaceSnapshot =
    detail::ImmersedFlowPredictorWorkspaceData;
using ImmersedFlowOperatorDiagonalContractReport =
    detail::ImmersedFlowOperatorDiagonalData;
using ImmersedFlowActiveExchangeReport =
    detail::ImmersedFlowActiveExchangeData;
using ImmersedFlowPressureAuthorityStateSnapshot =
    detail::ImmersedFlowPressureAuthorityData;
using ImmersedFlowPhysicalOutletPredictorSnapshot =
    detail::ImmersedFlowOutletPredictorValue;
using ImmersedFlowWallPredictorSnapshot = detail::ImmersedFlowWallPredictorValue;
using ImmersedFlowFinalMomentumPressureResidualReport =
    detail::ImmersedFlowFinalPressureResidualData;
using ImmersedFlowPressureFluxIdentityReport =
    detail::ImmersedFlowPressureFluxIdentityData;
using ImmersedFlowCellPressureCorrectionRecord =
    detail::ImmersedFlowPressureCorrectionRecord;
using ImmersedFlowCellPressureCorrectionReport =
    detail::ImmersedFlowPressureCorrectionData;
using ImmersedFlowSpatialEnergyReport = detail::ImmersedFlowSpatialEnergyData;
using ImmersedFlowExactMomentumResidualReport =
    detail::ImmersedFlowExactMomentumResidualData;

class ImmersedFlowTestAccess final {
public:
  static void set_linear_solve_failure(ImmersedLinearSolvePhase,
                                       std::uint32_t component) noexcept;
  static void clear_linear_solve_failure() noexcept;
  static void set_failure_stage(ImmersedFlowAttemptFailureStage,
                                int failing_rank = -1) noexcept;
  static void clear_failure_stage() noexcept;
  static void set_wall_input_failure(ImmersedFlowWallInputFailure,
                                     int failing_rank) noexcept;
  static void clear_wall_input_failure() noexcept;
  static void set_momentum_time_diagonal_scale(double) noexcept;
  static void clear_momentum_time_diagonal_scale() noexcept;
  static void
  set_manufactured_body_source(FixedStepImmersedFlow &,
                               std::vector<double> owned_active_source);
  static void clear_manufactured_body_source(FixedStepImmersedFlow &) noexcept;
  static std::uint64_t pressure_revision(const FixedStepImmersedFlow &) noexcept;
  static bool
  has_active_pressure_reference(const FixedStepImmersedFlow &) noexcept;
  static std::uint64_t
  pressure_apply_schedule(const FixedStepImmersedFlow &) noexcept;
  static std::uint32_t
  last_corrector_count(const FixedStepImmersedFlow &) noexcept;
  static std::uint64_t
  wale_evaluation_count(const FixedStepImmersedFlow &) noexcept;
  static les::WaleCoefficientIdentity
  wale_coefficient_identity(const FixedStepImmersedFlow &) noexcept;
  static std::uint64_t
  wall_effective_viscosity_fingerprint(
      const FixedStepImmersedFlow &) noexcept;
  static double last_wall_pressure_gradient_application_norm(
      const FixedStepImmersedFlow &) noexcept;
  static std::vector<ImmersedFlowWallMeasureSnapshot>
  wall_effective_measures(const FixedStepImmersedFlow &);
  static std::uint64_t
  immersed_operator_structure_fingerprint(const FixedStepImmersedFlow &);
  static const finite_volume::ImmersedOperatorAdapter &
  immersed_operator(const FixedStepImmersedFlow &);
  static std::vector<ImmersedFlowWallGradientSnapshot>
  final_pressure_wall_gradients(const FixedStepImmersedFlow &);
  static std::array<double, 3>
  active_pressure_algebra_probe(const FixedStepImmersedFlow &);
  static std::array<double, 3>
  active_pressure_current_algebra_probe(const FixedStepImmersedFlow &);
  static ImmersedFlowActiveExchangeReport
  active_exchange_performance(const FixedStepImmersedFlow &) noexcept;
  static std::array<double, 9> active_momentum_reference_probe(
      const FixedStepImmersedFlow &, double rho_ref_kg_per_m3,
      double dynamic_viscosity_pa_s, const MomentumTimeStencil &);
  static std::vector<double> active_pressure_operator_values(
      const FixedStepImmersedFlow &,
      std::vector<double> owned_active_pressure_pa,
      std::vector<ImmersedFlowWallGradientSnapshot> wall_gradients);
  static ImmersedFlowExactPredictorWorkspaceSnapshot
  exact_predictor_workspace(const FixedStepImmersedFlow &) noexcept;
  static ImmersedFlowOperatorDiagonalContractReport
  exact_predictor_diagonal_contract(const FixedStepImmersedFlow &);
  static ImmersedFlowPressureFluxIdentityReport exact_predictor_schur_identity(
      FixedStepImmersedFlow &, const FlowState &, double rho_ref_kg_per_m3,
      const MomentumTimeStencil &,
      std::vector<double> owned_active_pressure_pa,
      std::vector<ImmersedFlowWallGradientSnapshot> wall_gradients);
  static ImmersedFlowCellPressureCorrectionReport
  cell_pressure_correction_authority(const FixedStepImmersedFlow &);
  static ImmersedFlowSpatialEnergyReport
  spatial_energy_terms(FixedStepImmersedFlow &, const FlowState &, FlowLayer,
                       double rho_ref_kg_per_m3, double dynamic_viscosity_pa_s);
  static ImmersedFlowExactMomentumResidualReport exact_momentum_residual_terms(
      FixedStepImmersedFlow &, const FlowState &,
      const FlowLayerValues &exact_trial,
      const MomentumTimeStencil &, double rho_ref_kg_per_m3,
      double dynamic_viscosity_pa_s,
      const std::vector<double> &owned_active_source_N_per_m3,
      std::array<std::vector<ImmersedFlowWallGradientSnapshot>, 3>
          wall_pressure_gradients);
  static bool exact_momentum_residual_report_is_self_consistent(
      const ImmersedFlowExactMomentumResidualReport &) noexcept;
  static std::vector<ImmersedFlowPhysicalOutletPredictorSnapshot>
  physical_outlet_predictor_values(FixedStepImmersedFlow &, FlowState &,
                                   double rho_ref_kg_per_m3,
                                   const MomentumTimeStencil &,
                                   double momentum_diagonal);
  static std::vector<ImmersedFlowWallPredictorSnapshot>
  wall_predictor_values(FixedStepImmersedFlow &, FlowState &,
                        double rho_ref_kg_per_m3,
                        const MomentumTimeStencil &, double momentum_diagonal);
  static ImmersedFlowFinalMomentumPressureResidualReport
  final_momentum_pressure_residual_routes(FixedStepImmersedFlow &, FlowState &);
  static double inter_step_pressure_authority_difference_l2(
      const FixedStepImmersedFlow &) noexcept;
  static double bdf2_history_pressure_authority_difference_l2(
      const FixedStepImmersedFlow &) noexcept;
  static std::array<std::uint64_t, 2>
  pressure_authority_fingerprints(const FixedStepImmersedFlow &) noexcept;
  static ImmersedFlowPressureAuthorityStateSnapshot
  pressure_authority_state(const FixedStepImmersedFlow &);
  static ForceAttemptReport assemble_force_attempt_report_from_budget(
      const immersed::ForceComponents &budget_reaction,
      const immersed::ForceComponents &surface_traction);
};

inline void ImmersedFlowTestAccess::set_linear_solve_failure(
    ImmersedLinearSolvePhase phase, std::uint32_t component) noexcept {
  detail::ImmersedFlowAccess::set_linear_solve_failure(phase, component);
}

inline void ImmersedFlowTestAccess::clear_linear_solve_failure() noexcept {
  detail::ImmersedFlowAccess::clear_linear_solve_failure();
}

inline void ImmersedFlowTestAccess::set_failure_stage(
    ImmersedFlowAttemptFailureStage stage, int failing_rank) noexcept {
  detail::ImmersedFlowAccess::set_failure_stage(
      static_cast<std::uint8_t>(stage), failing_rank);
}

inline void ImmersedFlowTestAccess::clear_failure_stage() noexcept {
  detail::ImmersedFlowAccess::clear_failure_stage();
}

inline void ImmersedFlowTestAccess::set_wall_input_failure(
    ImmersedFlowWallInputFailure failure, int failing_rank) noexcept {
  detail::ImmersedFlowAccess::set_wall_input_failure(
      static_cast<std::uint8_t>(failure), failing_rank);
}

inline void ImmersedFlowTestAccess::clear_wall_input_failure() noexcept {
  detail::ImmersedFlowAccess::clear_wall_input_failure();
}

inline void
ImmersedFlowTestAccess::set_momentum_time_diagonal_scale(double scale) noexcept {
  detail::ImmersedFlowAccess::set_momentum_time_diagonal_scale(scale);
}

inline void
ImmersedFlowTestAccess::clear_momentum_time_diagonal_scale() noexcept {
  detail::ImmersedFlowAccess::clear_momentum_time_diagonal_scale();
}

inline void ImmersedFlowTestAccess::set_manufactured_body_source(
    FixedStepImmersedFlow &flow, std::vector<double> source) {
  detail::ImmersedFlowAccess::set_manufactured_body_source(flow,
                                                         std::move(source));
}

inline void ImmersedFlowTestAccess::clear_manufactured_body_source(
    FixedStepImmersedFlow &flow) noexcept {
  detail::ImmersedFlowAccess::clear_manufactured_body_source(flow);
}

inline std::uint64_t ImmersedFlowTestAccess::pressure_revision(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::pressure_revision(flow);
}

inline bool ImmersedFlowTestAccess::has_active_pressure_reference(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::has_active_pressure_reference(flow);
}

inline std::uint64_t ImmersedFlowTestAccess::pressure_apply_schedule(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::pressure_apply_schedule(flow);
}

inline std::uint32_t ImmersedFlowTestAccess::last_corrector_count(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::last_corrector_count(flow);
}

inline std::uint64_t ImmersedFlowTestAccess::wale_evaluation_count(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::wale_evaluation_count(flow);
}

inline les::WaleCoefficientIdentity
ImmersedFlowTestAccess::wale_coefficient_identity(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::wale_coefficient_identity(flow);
}

inline std::uint64_t
ImmersedFlowTestAccess::wall_effective_viscosity_fingerprint(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::wall_effective_viscosity_fingerprint(flow);
}

inline double
ImmersedFlowTestAccess::last_wall_pressure_gradient_application_norm(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::
      last_wall_pressure_gradient_application_norm(flow);
}

inline std::vector<ImmersedFlowWallMeasureSnapshot>
ImmersedFlowTestAccess::wall_effective_measures(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::wall_effective_measures(flow);
}

inline std::uint64_t
ImmersedFlowTestAccess::immersed_operator_structure_fingerprint(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::immersed_operator_structure_fingerprint(
      flow);
}

inline const finite_volume::ImmersedOperatorAdapter &
ImmersedFlowTestAccess::immersed_operator(const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::immersed_operator(flow);
}

inline std::vector<ImmersedFlowWallGradientSnapshot>
ImmersedFlowTestAccess::final_pressure_wall_gradients(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::final_pressure_wall_gradients(flow);
}

inline std::array<double, 3>
ImmersedFlowTestAccess::active_pressure_algebra_probe(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::active_pressure_algebra_probe(flow);
}

inline std::array<double, 3>
ImmersedFlowTestAccess::active_pressure_current_algebra_probe(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::active_pressure_current_algebra_probe(
      flow);
}

inline ImmersedFlowActiveExchangeReport
ImmersedFlowTestAccess::active_exchange_performance(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::active_exchange_performance(flow);
}

inline std::array<double, 9>
ImmersedFlowTestAccess::active_momentum_reference_probe(
    const FixedStepImmersedFlow &flow, double rho_ref_kg_per_m3,
    double dynamic_viscosity_pa_s, const MomentumTimeStencil &stencil) {
  return detail::ImmersedFlowAccess::active_momentum_reference_probe(
      flow, rho_ref_kg_per_m3, dynamic_viscosity_pa_s, stencil);
}

inline std::vector<double> ImmersedFlowTestAccess::active_pressure_operator_values(
    const FixedStepImmersedFlow &flow,
    std::vector<double> owned_active_pressure_pa,
    std::vector<ImmersedFlowWallGradientSnapshot> wall_gradients) {
  return detail::ImmersedFlowAccess::active_pressure_operator_values(
      flow, std::move(owned_active_pressure_pa), std::move(wall_gradients));
}

inline ImmersedFlowExactPredictorWorkspaceSnapshot
ImmersedFlowTestAccess::exact_predictor_workspace(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::exact_predictor_workspace(flow);
}

inline ImmersedFlowOperatorDiagonalContractReport
ImmersedFlowTestAccess::exact_predictor_diagonal_contract(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::exact_predictor_diagonal_contract(flow);
}

inline ImmersedFlowPressureFluxIdentityReport
ImmersedFlowTestAccess::exact_predictor_schur_identity(
    FixedStepImmersedFlow &flow, const FlowState &state,
    double rho_ref_kg_per_m3, const MomentumTimeStencil &stencil,
    std::vector<double> owned_active_pressure_pa,
    std::vector<ImmersedFlowWallGradientSnapshot> wall_gradients) {
  return detail::ImmersedFlowAccess::exact_predictor_schur_identity(
      flow, state, rho_ref_kg_per_m3, stencil,
      std::move(owned_active_pressure_pa), std::move(wall_gradients));
}

inline ImmersedFlowCellPressureCorrectionReport
ImmersedFlowTestAccess::cell_pressure_correction_authority(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::cell_pressure_correction_authority(flow);
}

inline ImmersedFlowSpatialEnergyReport ImmersedFlowTestAccess::spatial_energy_terms(
    FixedStepImmersedFlow &flow, const FlowState &state, FlowLayer layer,
    double rho_ref_kg_per_m3, double dynamic_viscosity_pa_s) {
  return detail::ImmersedFlowAccess::spatial_energy_terms(
      flow, state, layer, rho_ref_kg_per_m3, dynamic_viscosity_pa_s);
}

inline ImmersedFlowExactMomentumResidualReport
ImmersedFlowTestAccess::exact_momentum_residual_terms(
    FixedStepImmersedFlow &flow, const FlowState &state,
    const FlowLayerValues &exact_trial, const MomentumTimeStencil &stencil,
    double rho_ref_kg_per_m3, double dynamic_viscosity_pa_s,
    const std::vector<double> &owned_active_source_N_per_m3,
    std::array<std::vector<ImmersedFlowWallGradientSnapshot>, 3>
        wall_pressure_gradients) {
  return detail::ImmersedFlowAccess::exact_momentum_residual_terms(
      flow, state, exact_trial, stencil, rho_ref_kg_per_m3,
      dynamic_viscosity_pa_s, owned_active_source_N_per_m3,
      std::move(wall_pressure_gradients));
}

inline bool ImmersedFlowTestAccess::
    exact_momentum_residual_report_is_self_consistent(
        const ImmersedFlowExactMomentumResidualReport &report) noexcept {
  return detail::ImmersedFlowAccess::
      exact_momentum_residual_report_is_self_consistent(report);
}

inline std::vector<ImmersedFlowPhysicalOutletPredictorSnapshot>
ImmersedFlowTestAccess::physical_outlet_predictor_values(
    FixedStepImmersedFlow &flow, FlowState &state, double rho_ref_kg_per_m3,
    const MomentumTimeStencil &stencil, double momentum_diagonal) {
  return detail::ImmersedFlowAccess::physical_outlet_predictor_values(
      flow, state, rho_ref_kg_per_m3, stencil, momentum_diagonal);
}

inline std::vector<ImmersedFlowWallPredictorSnapshot>
ImmersedFlowTestAccess::wall_predictor_values(
    FixedStepImmersedFlow &flow, FlowState &state, double rho_ref_kg_per_m3,
    const MomentumTimeStencil &stencil, double momentum_diagonal) {
  return detail::ImmersedFlowAccess::wall_predictor_values(
      flow, state, rho_ref_kg_per_m3, stencil, momentum_diagonal);
}

inline ImmersedFlowFinalMomentumPressureResidualReport
ImmersedFlowTestAccess::final_momentum_pressure_residual_routes(
    FixedStepImmersedFlow &flow, FlowState &state) {
  return detail::ImmersedFlowAccess::final_momentum_pressure_residual_routes(
      flow, state);
}

inline double
ImmersedFlowTestAccess::inter_step_pressure_authority_difference_l2(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::
      inter_step_pressure_authority_difference_l2(flow);
}

inline double
ImmersedFlowTestAccess::bdf2_history_pressure_authority_difference_l2(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::
      bdf2_history_pressure_authority_difference_l2(flow);
}

inline std::array<std::uint64_t, 2>
ImmersedFlowTestAccess::pressure_authority_fingerprints(
    const FixedStepImmersedFlow &flow) noexcept {
  return detail::ImmersedFlowAccess::pressure_authority_fingerprints(flow);
}

inline ImmersedFlowPressureAuthorityStateSnapshot
ImmersedFlowTestAccess::pressure_authority_state(
    const FixedStepImmersedFlow &flow) {
  return detail::ImmersedFlowAccess::pressure_authority_state(flow);
}

inline ForceAttemptReport
ImmersedFlowTestAccess::assemble_force_attempt_report_from_budget(
    const immersed::ForceComponents &budget_reaction,
    const immersed::ForceComponents &surface_traction) {
  return detail::ImmersedFlowAccess::assemble_force_attempt_report_from_budget(
      budget_reaction, surface_traction);
}

} // namespace hundun::flow::test
