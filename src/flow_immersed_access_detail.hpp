// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_immersed.hpp"
#include "hundun/fvm_immersed_operator.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace hundun::flow::detail {

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS

struct ImmersedFlowWallMeasure final {
  std::uint64_t link{};
  double effective_measure_m2{};
};

struct ImmersedFlowWallGradient final {
  std::uint64_t link{};
  double normal_gradient_pa_per_m{};
};

struct ImmersedFlowPredictorWorkspaceData final {
  std::uint64_t probe_state_creation_count{};
  std::uint64_t response_capacity_growth_count{};
  std::uint64_t cached_response_consumption_count{};
  bool cached_input_mutation_rejected{};
  std::uint64_t last_response_allocation_event_count{};
  std::uint64_t accepted_inner_response_count{};
  std::uint64_t accepted_inner_solve_count{};
  std::uint64_t accepted_inner_iteration_count{};
  std::uint64_t accepted_inner_matvec_count{};
  std::uint64_t accepted_inner_preconditioner_apply_count{};
  std::uint64_t accepted_inner_global_reduction_count{};
};

struct ImmersedFlowOperatorDiagonalData final {
  bool advertised{};
  double declared_value{};
  double applied_basis_value{};
  double relative_difference{};
};

struct ImmersedFlowPressureAuthorityData final {
  std::vector<ImmersedFlowWallGradient> history;
  std::vector<ImmersedFlowWallGradient> committed;
  std::vector<ImmersedFlowWallGradient> pending;
  bool history_available{};
  bool committed_available{};
};

struct ImmersedFlowOutletPredictorValue final {
  std::uint64_t global_face_id{};
  runtime::Real3 velocity_m_per_s{};
  double mass_flux_kg_per_s{};
};

struct ImmersedFlowWallPredictorValue final {
  std::uint64_t link{};
  double predictor_mass_flux_kg_per_s{};
};

struct ImmersedFlowFinalPressureResidualData final {
  std::vector<std::uint64_t> active_global_cell_ids;
  std::vector<double> used_residual;
  std::vector<double> direct_residual;
  double used_l2{};
  double direct_l2{};
  double difference_l2{};
  double used_linf{};
  double direct_linf{};
  double difference_linf{};
  std::uint64_t maximum_difference_global_cell_id{};
  std::uint32_t maximum_difference_component{};
  double maximum_difference_used_value{};
  double maximum_difference_direct_value{};
};

struct ImmersedFlowPressureFluxIdentityData final {
  std::vector<std::uint64_t> active_global_cell_ids;
  std::vector<double> operator_residual_per_volume;
  std::vector<double> routed_flux_divergence_per_volume;
  double operator_l2{};
  double routed_l2{};
  double difference_l2{};
  double operator_linf{};
  double routed_linf{};
  double difference_linf{};
  std::uint64_t maximum_difference_global_cell_id{};
  double maximum_difference_operator_value{};
  double maximum_difference_routed_value{};
};

struct ImmersedFlowPressureCorrectionRecord final {
  std::vector<std::uint64_t> active_global_cell_ids;
  std::vector<double> exact_velocity_change;
  std::vector<double> momentum_operator_velocity_change;
  std::vector<double> lfp_pressure_residual_change;
  std::vector<double> pressure_before_pa;
  std::vector<double> pressure_correction_pa;
  std::vector<double> pressure_after_pa;
  std::vector<double> pressure_rhs_per_volume;
  std::vector<double> compact_pressure_action_per_volume;
  std::vector<double> interface_pressure_action_per_volume;
  std::vector<double> hybrid_pressure_action_per_volume;
  std::vector<double> affine_wall_source_per_volume;
  std::vector<double> compact_predictor_defect_per_volume;
  std::vector<double> hybrid_predictor_defect_per_volume;
  double momentum_l2{};
  double pressure_l2{};
  double closure_l2{};
  double momentum_linf{};
  double pressure_linf{};
  double closure_linf{};
  std::uint64_t maximum_closure_global_cell_id{};
  std::uint32_t maximum_closure_component{};
  double maximum_closure_momentum_value{};
  double maximum_closure_pressure_value{};
  double wall_predictor_mass_flux_l2_kg_per_s{};
  double wall_predictor_mass_flux_linf_kg_per_s{};
  double wall_correction_gradient_l2_pa_per_m{};
  double wall_correction_gradient_linf_pa_per_m{};
};

struct ImmersedFlowPressureCorrectionData final {
  std::vector<ImmersedFlowPressureCorrectionRecord> correctors;
  double inter_corrector_authority_difference_l2{};
};

struct ImmersedFlowSpatialEnergyData final {
  std::vector<std::uint64_t> active_global_cell_ids;
  std::vector<double> cell_volume_m3;
  std::vector<double> velocity_m_per_s;
  std::vector<double> pressure_pa;
  std::vector<double> mass_divergence_kg_per_s;
  std::vector<double> total_residual_N;
  std::vector<double> convective_residual_N;
  std::vector<double> pressure_residual_N;
  std::vector<double> viscous_residual_N;
  std::vector<double> implicit_viscous_reference_residual_N;
  double kinetic_energy_J{};
  double total_power_W{};
  double convective_power_W{};
  double centered_convective_power_W{};
  double reconstruction_power_W{};
  double pressure_power_W{};
  double pressure_continuity_power_W{};
  double pressure_adjoint_defect_W{};
  double viscous_power_W{};
  double implicit_viscous_reference_power_W{};
  double residual_closure_l2_N{};
  double residual_closure_linf_N{};
  double mass_divergence_l2_kg_per_s{};
  double mass_divergence_linf_kg_per_s{};
  double stationary_wall_flux_linf_kg_per_s{};
  double physical_boundary_flux_linf_kg_per_s{};
  std::uint64_t maximum_closure_global_cell_id{};
  std::uint32_t maximum_closure_component{};
  double maximum_closure_total_value_N{};
  double maximum_closure_parts_value_N{};
};

struct ImmersedFlowExactMomentumResidualData final {
  std::vector<std::uint64_t> active_global_cell_ids;
  std::vector<std::uint8_t> immersed_interface_row;
  std::vector<double> cell_volume_m3;
  std::vector<double> time_residual_N;
  std::vector<double> convective_residual_N;
  std::vector<double> viscous_remainder_residual_N;
  std::vector<double> pressure_residual_N;
  std::vector<double> background_pressure_residual_N;
  std::vector<double> pressure_wall_defect_residual_N;
  std::vector<double> implicit_viscous_reference_residual_N;
  std::vector<double> source_residual_N;
  std::vector<double> total_residual_N;
  double pointwise_split_closure_linf_N{};
};

struct ImmersedFlowAccess final {
  static void set_failure_stage(std::uint8_t, int) noexcept;
  static void clear_failure_stage() noexcept;
  static void set_wall_input_failure(std::uint8_t, int) noexcept;
  static void clear_wall_input_failure() noexcept;
  static void set_momentum_time_diagonal_scale(double) noexcept;
  static void clear_momentum_time_diagonal_scale() noexcept;
  static void set_manufactured_body_source(FixedStepImmersedFlow &,
                                           std::vector<double>);
  static void clear_manufactured_body_source(FixedStepImmersedFlow &) noexcept;
  static std::uint64_t pressure_revision(const FixedStepImmersedFlow &) noexcept;
  static bool
  has_active_pressure_reference(const FixedStepImmersedFlow &) noexcept;
  static std::uint64_t
  pressure_apply_schedule(const FixedStepImmersedFlow &) noexcept;
  static std::uint32_t
  last_corrector_count(const FixedStepImmersedFlow &) noexcept;
  static double last_wall_pressure_gradient_application_norm(
      const FixedStepImmersedFlow &) noexcept;
  static std::vector<ImmersedFlowWallMeasure>
  wall_effective_measures(const FixedStepImmersedFlow &);
  static std::uint64_t
  immersed_operator_structure_fingerprint(const FixedStepImmersedFlow &);
  static const finite_volume::ImmersedOperatorAdapter &
  immersed_operator(const FixedStepImmersedFlow &);
  static std::vector<ImmersedFlowWallGradient>
  final_pressure_wall_gradients(const FixedStepImmersedFlow &);
  static std::array<double, 3>
  active_pressure_algebra_probe(const FixedStepImmersedFlow &);
  static std::array<double, 9> active_momentum_reference_probe(
      const FixedStepImmersedFlow &, double, double, const MomentumTimeStencil &);
  static std::vector<double> active_pressure_operator_values(
      const FixedStepImmersedFlow &, std::vector<double>,
      std::vector<ImmersedFlowWallGradient>);
  static ImmersedFlowPredictorWorkspaceData
  exact_predictor_workspace(const FixedStepImmersedFlow &) noexcept;
  static ImmersedFlowOperatorDiagonalData
  exact_predictor_diagonal_contract(const FixedStepImmersedFlow &);
  static ImmersedFlowPressureFluxIdentityData exact_predictor_schur_identity(
      FixedStepImmersedFlow &, const FlowState &, double,
      const MomentumTimeStencil &, std::vector<double>,
      std::vector<ImmersedFlowWallGradient>);
  static ImmersedFlowPressureCorrectionData
  cell_pressure_correction_authority(const FixedStepImmersedFlow &);
  static ImmersedFlowSpatialEnergyData
  spatial_energy_terms(FixedStepImmersedFlow &, const FlowState &, FlowLayer,
                       double, double);
  static ImmersedFlowExactMomentumResidualData exact_momentum_residual_terms(
      FixedStepImmersedFlow &, const FlowState &, const FlowLayerValues &,
      const MomentumTimeStencil &, double, double, const std::vector<double> &,
      std::array<std::vector<ImmersedFlowWallGradient>, 3>);
  static bool exact_momentum_residual_report_is_self_consistent(
      const ImmersedFlowExactMomentumResidualData &) noexcept;
  static std::vector<ImmersedFlowOutletPredictorValue>
  physical_outlet_predictor_values(FixedStepImmersedFlow &, FlowState &, double,
                                   const MomentumTimeStencil &, double);
  static std::vector<ImmersedFlowWallPredictorValue>
  wall_predictor_values(FixedStepImmersedFlow &, FlowState &, double,
                        const MomentumTimeStencil &, double);
  static ImmersedFlowFinalPressureResidualData
  final_momentum_pressure_residual_routes(FixedStepImmersedFlow &, FlowState &);
  static double inter_step_pressure_authority_difference_l2(
      const FixedStepImmersedFlow &) noexcept;
  static double bdf2_history_pressure_authority_difference_l2(
      const FixedStepImmersedFlow &) noexcept;
  static std::array<std::uint64_t, 2>
  pressure_authority_fingerprints(const FixedStepImmersedFlow &) noexcept;
  static ImmersedFlowPressureAuthorityData
  pressure_authority_state(const FixedStepImmersedFlow &);
  static ForceAttemptReport assemble_force_attempt_report_from_budget(
      const immersed::ForceComponents &, const immersed::ForceComponents &);
};

#endif

} // namespace hundun::flow::detail
