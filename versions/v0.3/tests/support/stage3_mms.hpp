// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/ib_wall_force.hpp"
#include "hundun/mesh_geometry.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/stage3_test_contracts.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace hundun::runtime {
class MpiContext;
}

namespace hundun::test::stage3 {

struct MmsSample final {
  runtime::Real3 velocity_m_per_s{};
  std::array<std::array<double, 3>, 3> velocity_gradient_per_s{};
  double mechanical_pressure_pa{};
  runtime::Real3 mechanical_pressure_gradient_pa_per_m{};
  std::array<std::array<double, 3>, 3> mechanical_pressure_hessian_pa_per_m2{};
  runtime::Real3 body_source_N_per_m3{};
};

struct MmsCellAverage final {
  runtime::Real3 velocity_m_per_s{};
  double mechanical_pressure_pa{};
  runtime::Real3 body_source_N_per_m3{};
  runtime::Real3 mechanical_pressure_gradient_pa_per_m{};
};

struct MmsFaceHistory final {
  runtime::Real3 velocity_m_per_s{};
  double owner_normal_volume_flux_m3_per_s{};
};

struct ManufacturedSurface final {
  std::vector<StlFixtureTriangle> triangles;
  double maximum_edge_m{};
  double maximum_chord_error_m{};
};

struct AnalyticForceComponents final {
  runtime::Real3 pressure_N{};
  runtime::Real3 viscous_N{};
  runtime::Real3 total_N{};
};

struct AnalyticForceReference final {
  AnalyticForceComponents force;
  runtime::Real3 viscous_absolute_component_traction_force_N{};
  double viscous_traction_l1_force_N{};
  double surface_area_m2{};
  double pressure_traction_rms_force_N{};
  double viscous_traction_rms_force_N{};
  double total_traction_rms_force_N{};
};

enum class ManufacturedMapping : std::uint8_t { uniform, warped };

enum class ManufacturedSurfacePolicy : std::uint8_t {
  per_level,
  fixed_48,
};

struct ManufacturedCase final {
  BodySpec body;
  int cells{};
  ManufacturedMapping mapping{ManufacturedMapping::uniform};
  config::ImmersedFluidSide fluid_side{config::ImmersedFluidSide::outside};
  runtime::Int3 process_grid{1, 1, 1};
  bool collect_force{true};
  bool collect_exact_momentum_residual{};
  std::optional<int> final_force_failure_rank;
  ManufacturedSurfacePolicy surface_policy{
      ManufacturedSurfacePolicy::per_level};
  bool collect_pressure_extrema_diagnostics{};
};

struct ExactMomentumResidualDiagnostics final {
  // Term order: time, convection, viscous remainder, pressure,
  // implicit viscous reference, source, total.
  std::array<double, 7> interface_rms_N_per_m3{};
  std::array<double, 7> bulk_rms_N_per_m3{};
  double interface_pressure_balance_rms_N_per_m3{};
  double bulk_pressure_balance_rms_N_per_m3{};
  double interface_background_pressure_balance_rms_N_per_m3{};
  double bulk_background_pressure_balance_rms_N_per_m3{};
  double interface_reconstructed_face_pressure_balance_rms_N_per_m3{};
  double bulk_reconstructed_face_pressure_balance_rms_N_per_m3{};
  double interface_analytic_face_pressure_balance_rms_N_per_m3{};
  double bulk_analytic_face_pressure_balance_rms_N_per_m3{};
  double interface_reconstructed_to_analytic_face_rms_N_per_m3{};
  double interface_unconstrained_center_pressure_balance_rms_N_per_m3{};
  double bulk_unconstrained_center_pressure_balance_rms_N_per_m3{};
  double interface_analytic_center_pressure_balance_rms_N_per_m3{};
  double bulk_analytic_center_pressure_balance_rms_N_per_m3{};
  double interface_constrained_to_unconstrained_center_rms_N_per_m3{};
  double interface_exact_shared_center_pressure_balance_rms_N_per_m3{};
  double bulk_exact_shared_center_pressure_balance_rms_N_per_m3{};
  double interface_shared_center_correction_rms_N_per_m3{};
  double shared_center_correction_conservation_linf_N{};
  double interface_full_exact_center_pressure_balance_rms_N_per_m3{};
  double bulk_full_exact_center_pressure_balance_rms_N_per_m3{};
  double interface_full_exact_integral_pressure_balance_rms_N_per_m3{};
  double bulk_full_exact_integral_pressure_balance_rms_N_per_m3{};
  double interface_background_to_full_exact_integral_rms_N_per_m3{};
  double single_link_current_background_pressure_balance_rms_N_per_m3{};
  double single_link_full_exact_center_pressure_balance_rms_N_per_m3{};
  double single_link_coherent_pressure_balance_rms_N_per_m3{};
  double single_link_coherent_integral_pressure_balance_rms_N_per_m3{};
  double single_link_unconstrained_integral_pressure_balance_rms_N_per_m3{};
  double single_link_constraint_integral_difference_rms_N_per_m3{};
  std::uint64_t single_link_row_count{};
  std::uint64_t multi_link_row_count{};
  double interface_pressure_wall_defect_rms_N_per_m3{};
  double bulk_pressure_wall_defect_rms_N_per_m3{};
  // Order: assembled pressure, background pressure, full exact face-centre
  // pressure, and immersed-wall correction.
  std::array<double, 4> interface_pressure_balance_linf_N_per_m3{};
  double face_area_closure_linf_m2{};
  double face_correction_closure_linf_N{};
  double pointwise_split_closure_linf_N{};
  std::uint64_t interface_row_count{};
  std::uint64_t bulk_row_count{};
  bool available{};
};

struct ManufacturedErrors final {
  double velocity_l2{};
  double velocity_linf{};
  double near_wall_velocity_l2{};
  double pressure_l2{};
  double pressure_linf{};
  double near_wall_pressure_l2{};
  double penetration_l2{};
  double penetration_linf{};
  double pressure_force{};
  double viscous_force{};
  double viscous_traction_l2{};
  double total_force{};
  double pressure_consistency{};
  double viscous_consistency{};
  double total_consistency{};
};

struct PressureErrorExtremum final {
  double signed_error_pa{};
  double absolute_error_pa{};
  mesh::GlobalCellId global_cell_id{};
  runtime::Int3 logical_cell{};
  double wall_distance_m{};
};

struct PressureErrorMoments final {
  double signed_error_volume_pa_m3{};
  double squared_error_volume_pa2_m3{};
  double volume_m3{};
  std::uint64_t cell_count{};
};

struct PressureErrorStatistics final {
  double mean_error_pa{};
  double rms_error_pa{};
  double centered_rms_error_pa{};
  double volume_m3{};
  std::uint64_t cell_count{};
};

struct NearWallPressureDiagnostics final {
  std::array<PressureErrorMoments, 4> distance_bands;
  std::array<PressureErrorMoments, 4> incident_distance_bands;
  std::array<PressureErrorMoments, 4> nonincident_distance_bands;
  PressureErrorMoments total;
  PressureErrorMoments incident_total;
  PressureErrorMoments nonincident_total;
  bool available{};
};

struct WallTractionJumpDiagnostics final {
  double same_link_rms{};
  double same_link_max{};
  double cross_link_same_donor_rms{};
  double cross_link_same_donor_max{};
  double cross_donor_rms{};
  double cross_donor_max{};
  std::uint64_t same_link_pairs{};
  std::uint64_t cross_link_same_donor_pairs{};
  std::uint64_t cross_donor_pairs{};
};

struct WallFunctionalConditioningDiagnostics final {
  double normal_gradient_amplification_rms{};
  double normal_gradient_amplification_max{};
  double condition_estimate_max{};
  runtime::Real3 signed_traction_error_N{};
  runtime::Real3 absolute_component_traction_error_N{};
  double absolute_traction_error_N{};
  double cancellation_ratio{};
  std::uint64_t point_count{};
};

struct WallTractionErrorDecomposition final {
  double exact_cell_reconstruction_l2{};
  double numerical_cell_contamination_l2{};
  double total_error_l2{};
  double pointwise_split_closure_linf{};
  std::uint64_t point_count{};
};

struct PressureExtremumAuthorityDiagnostics final {
  immersed::ImmersedLinkId nearest_link{};
  int owner_rank{-1};
  bool direct_fluid_cell{};
  double cell_to_wall_distance_m{};
  double condition_estimate{};
  std::uint64_t donor_count{};
  std::uint64_t row_link_count{};
  std::uint64_t donor_fingerprint{};
  double pressure_value_amplification{};
  double normal_gradient_amplification{};
};

struct PressureErrorExtremumDetail final {
  PressureErrorExtremum error;
  double wall_distance_over_h{};
  std::uint64_t incident_wall_link_count{};
  bool authority_available{};
  PressureExtremumAuthorityDiagnostics authority;
};

struct PressureMeasureDiagnostics final {
  runtime::Real3 background_reaction_N{};
  runtime::Real3 projected_reaction_N{};
  runtime::Real3 surface_partition_reaction_N{};
  runtime::Real3 exact_state_a22_reaction_N{};
  runtime::Real3 exact_state_a22_donor_reaction_N{};
  runtime::Real3 exact_state_a22_wall_reaction_N{};
  runtime::Real3 background_consistency_N{};
  runtime::Real3 projected_consistency_N{};
  runtime::Real3 surface_partition_consistency_N{};
  runtime::Real3 exact_state_a22_consistency_N{};
  runtime::Real3 numerical_state_a22_defect_N{};
  double product_reproduction_linf_N{};
  double background_area_closure{};
  double projected_area_closure{};
  double surface_partition_area_closure{};
  double row_wall_value_linf{};
  double row_wall_value_l2{};
  double per_link_wall_value_linf{};
  double per_link_wall_value_l2{};
  double per_link_origin_g_wall_value_linf{};
  double per_link_origin_g_wall_value_l2{};
  double linear_wall_value_linf{};
  double linear_wall_value_l2{};
  double linear_wall_value_hg_linf{};
  double plain_linear_wall_value_linf{};
  double linear_extrap_wall_linf{};
  double linear_extrap_wall_l2{};
  double constrained_origin_linear_linf{};
  double constrained_origin_hg_linear_linf{};
  double constant_field_linf{};
  double abs_sum_w_t1_linf{};
  double abs_sum_w_t2_linf{};
  double abs_sum_w_n_linf{};
  double authority_linear_linf{};
  double ghost_donor_wall_linf{};
  double ghost_donor_wall_l2{};
  std::uint64_t ghost_donor_link_count{};
  double linear_row_residual_linf{};
  double linear_row_residual_l2{};
  double mms_row_residual_linf{};
  double mms_row_residual_l2{};
  std::uint64_t mms_row_count{};
  double link_true_normal_mismatch_linf{};
  double mms_row_coherence_ratio{};
  runtime::Real3 linear_row_sum_N{};
  runtime::Real3 single_link_row_sum_N{};
  runtime::Real3 multi_link_row_sum_N{};
  runtime::Real3 single_link_donor_sum_N{};
  runtime::Real3 single_link_wall_sum_N{};
  runtime::Real3 single_link_linear_donor_sum_N{};
  runtime::Real3 single_link_donor_constant_sum_N{};
  runtime::Real3 single_link_donor_moment_x_sum_N{};
  runtime::Real3 single_link_donor_moment_y_sum_N{};
  runtime::Real3 single_link_donor_moment_z_sum_N{};
  runtime::Real3 single_link_wall_coefficient_sum_N{};
  runtime::Real3 single_link_donor_second_moment_sum_N{};
  runtime::Real3 single_link_authority_wall_flux_sum_N{};
  runtime::Real3 single_link_exact_face_flux_sum_N{};
  runtime::Real3 multi_link_exact_face_flux_sum_N{};
  runtime::Real3 single_link_authority_value_difference_l2_N{};
  std::uint64_t single_link_row_count{};
  std::uint64_t multi_link_row_count{};
  std::uint64_t linear_wall_value_count{};
  double authority_wall_value_linf{};
  double authority_wall_value_l2{};
  double authority_centroid_value_linf{};
  double authority_centroid_value_l2{};
  runtime::Real3 exact_centroid_quadrature_sum_N{};
  runtime::Real3 authority_centroid_quadrature_sum_N{};
  runtime::Real3 surface_vector_sum_N{};
  runtime::Real3 ghost_measure_sum_N{};
  runtime::Real3 ghost_authority_centroid_quadrature_sum_N{};
  double error_delta_pearson{};
  double error_obliqueness_pearson{};
  double error_curvature_pearson{};
  double error_delta_loglog_slope{};
  double worst_link_obliqueness{};
  std::uint64_t wall_value_link_count{};
  std::uint64_t link_count{};
  std::uint64_t surface_supported_link_count{};
  std::uint64_t missing_surface_link_count{};
  std::uint64_t multiple_surface_point_link_count{};
};

struct PressureScalarDiagnostics final {
  runtime::Real3 numerical_face_reaction_N{};
  runtime::Real3 exact_face_reaction_N{};
  runtime::Real3 exact_wall_reaction_N{};
  runtime::Real3 exact_hybrid_reaction_N{};
  runtime::Real3 numerical_projected_defect_reaction_N{};
  runtime::Real3 exact_projected_defect_reaction_N{};
  runtime::Real3 numerical_wall_projected_reaction_N{};
  runtime::Real3 exact_wall_projected_reaction_N{};
  runtime::Real3 numerical_face_consistency_N{};
  runtime::Real3 exact_face_consistency_N{};
  runtime::Real3 exact_wall_consistency_N{};
  runtime::Real3 exact_hybrid_consistency_N{};
  runtime::Real3 numerical_projected_defect_consistency_N{};
  runtime::Real3 exact_projected_defect_consistency_N{};
  runtime::Real3 numerical_wall_projected_consistency_N{};
  runtime::Real3 exact_wall_projected_consistency_N{};
};

struct ManufacturedRunResult final {
  ManufacturedErrors errors;
  NearWallPressureDiagnostics near_wall_pressure;
  PressureErrorExtremum pressure_error_extremum;
  std::vector<PressureErrorExtremumDetail> pressure_error_extrema_top_k;
  WallTractionJumpDiagnostics wall_traction_jump;
  WallFunctionalConditioningDiagnostics wall_functional_conditioning;
  WallTractionErrorDecomposition wall_traction_error_decomposition;
  PressureExtremumAuthorityDiagnostics pressure_extremum_authority;
  PressureMeasureDiagnostics pressure_measure;
  PressureScalarDiagnostics pressure_scalar;
  ExactMomentumResidualDiagnostics exact_momentum_residual;
  immersed::ForceComponents operator_force;
  immersed::ForceComponents budget_reaction;
  immersed::ForceComponents surface_traction;
  immersed::ForceComponents consistency;
  AnalyticForceComponents analytic_force;
  AnalyticForceComponents wall_plan_analytic_force;
  runtime::Real3 wall_plan_full_gradient_viscous_force_N{};
  runtime::Real3 wall_link_projected_viscous_force_N{};
  runtime::Real3 wall_link_full_gradient_viscous_force_N{};
  double pressure_force_scale_N{};
  double surface_maximum_edge_m{};
  double surface_maximum_chord_error_m{};
  double exact_velocity_linf{};
  double exact_pressure_linf{};
  double active_pressure_mean{};
  std::uint64_t active_cell_count{};
  std::uint64_t immersed_link_count{};
  std::uint64_t wall_point_count{};
  std::uint64_t classification_fingerprint{};
  std::uint64_t surface_coverage_fingerprint{};
  std::uint64_t ghost_plan_fingerprint{};
  std::uint64_t wall_plan_fingerprint{};
  std::uint64_t operator_structure_fingerprint{};
  int lowest_failing_rank{-1};
  bool committed{true};
  bool rollback_bitwise_equal{true};
  std::vector<mesh::GlobalCellId> global_active_cell_ids;
  std::vector<double> global_velocity;
  std::vector<double> global_pressure;
};

MmsSample evaluate_mms(const BodySpec &, runtime::Real3 point_m, double time_s);
bool pressure_face_quadrature_oracle_is_mutation_sensitive() noexcept;
double implicit_factor(const BodySpec &, runtime::Real3 point_m);
ManufacturedSurface make_manufactured_surface(const BodySpec &, double h_max_m);
AnalyticForceReference analytic_force_reference(const BodySpec &, double time_s,
                                                std::size_t quadrature_order);
runtime::Real3 absolute_viscous_traction_force_increment(
    runtime::Real3 traction, double positive_weight) noexcept;
bool signed_viscous_force_reference_preflight_accepts(
    const AnalyticForceReference &coarse,
    const AnalyticForceReference &refined) noexcept;
AnalyticForceReference verified_force_reference_preflight(const BodySpec &,
                                                          double time_s);
runtime::Real3 source_cell_average(const mesh::MeshTopology &,
                                   const mesh::MeshGeometry &,
                                   mesh::LocalCellId, const BodySpec &,
                                   double time_s);
MmsCellAverage evaluate_cell_average(const mesh::MeshTopology &,
                                     const mesh::MeshGeometry &,
                                     mesh::LocalCellId, const BodySpec &,
                                     double time_s);
MmsFaceHistory evaluate_face_history(const mesh::MeshTopology &,
                                     const mesh::MeshGeometry &,
                                     mesh::LocalFaceId, const BodySpec &,
                                     double time_s);
namespace detail {
// Precondition: geometry compatibility was established by the enclosing
// test-support operation.
MmsCellAverage evaluate_cell_average_validated(const mesh::MeshTopology &,
                                               const mesh::MeshGeometry &,
                                               mesh::LocalCellId,
                                               const BodySpec &, double time_s);
MmsFaceHistory evaluate_face_history_validated(const mesh::MeshTopology &,
                                               const mesh::MeshGeometry &,
                                               mesh::LocalFaceId,
                                               const BodySpec &, double time_s);
} // namespace detail
bool finite(const MmsSample &) noexcept;
double max_abs(runtime::Real3) noexcept;
double max_abs_difference(runtime::Real3, runtime::Real3) noexcept;
std::optional<PressureErrorExtremum> select_pressure_error_extremum(
    const std::vector<PressureErrorExtremum> &candidates) noexcept;
std::vector<PressureErrorExtremum> select_pressure_error_extrema(
    const std::vector<PressureErrorExtremum> &candidates,
    std::size_t limit) noexcept;
std::optional<std::size_t>
normalized_near_wall_pressure_band(double wall_distance_over_h) noexcept;
bool manufactured_near_wall_band_contains(double wall_distance_m) noexcept;
bool accumulate_pressure_error_moments(PressureErrorMoments &moments,
                                       double signed_error_pa,
                                       double volume_m3) noexcept;
std::optional<PressureErrorStatistics>
summarize_pressure_error_moments(const PressureErrorMoments &moments) noexcept;
ManufacturedRunResult run_manufactured_case(const runtime::MpiContext &,
                                            const ManufacturedCase &);

} // namespace hundun::test::stage3
