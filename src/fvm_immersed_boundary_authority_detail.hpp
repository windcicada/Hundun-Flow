// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/fvm_immersed_operator.hpp"
#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_quadratic_reconstruction.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::finite_volume::detail {

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS

struct ImmersedBoundaryWallLinkRecord final {
  immersed::ImmersedLinkId id{};
  std::size_t occurrence{};
  double normal_scale{};
  runtime::Real3 solid_to_fluid_normal{};
  std::uint64_t replacement_fingerprint{};
  runtime::Real3 wall_intercept_m{};
  runtime::Real3 area_from_fluid_m2{};
  double signed_wall_measure_m2{};
  runtime::Real3 pressure_quadrature_m{};
  runtime::Real3 surface_measure_m2{};
  runtime::Real3 surface_patch_centroid_m{};
};

struct ImmersedBoundaryPhysicalTermRecord final {
  std::uint64_t stable_term_id{};
  immersed::ImmersedLinkId link{};
  std::uint8_t kind{};
  std::uint32_t algebraic_occurrence{};
  std::uint32_t output_component{};
  double coefficient{};
  std::uint64_t evaluation_group_id{};
  std::vector<std::uint64_t> source_term_ids;
};

struct ImmersedBoundaryRowRecord final {
  mesh::GlobalCellId active_cell{};
  std::vector<ImmersedBoundaryWallLinkRecord> links;
  std::uint64_t row_replacement_fingerprint{};
  std::uint64_t replacement_group_count{};
  std::vector<ImmersedBoundaryPhysicalTermRecord> covered_physical_terms;
};

struct ImmersedBoundaryReplacementTermRecord final {
  std::uint64_t stable_term_id{};
  immersed::ImmersedLinkId link{};
  std::uint8_t kind{};
  std::uint32_t occurrence{};
  std::uint32_t component{};
  double value{};
  std::uint64_t evaluation_group_id{};
};

struct ImmersedBoundaryResidualPartsRecord final {
  std::array<double, 3> convective{};
  std::array<double, 3> pressure{};
  std::array<double, 3> viscous{};
};

struct ImmersedBoundaryAffineDonorTermRecord final {
  std::uint64_t stable_term_id{};
  std::uint8_t input_kind{};
  mesh::GlobalCellId donor_global_cell{};
  std::uint32_t input_component{};
  std::uint32_t output_component{};
  double coefficient{};
  std::uint32_t contributing_link_count{};
  std::vector<std::uint64_t> source_term_ids;
};

struct ImmersedBoundaryAffineWallGradientTermRecord final {
  std::uint64_t stable_term_id{};
  immersed::ImmersedLinkId link{};
  std::uint32_t output_component{};
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

struct ImmersedBoundaryRowEvaluationRecord final {
  mesh::GlobalCellId active_cell{};
  std::uint64_t row_fingerprint{};
  std::vector<ImmersedBoundaryReplacementTermRecord> replacement_terms;
  std::array<double, 3> residual_before_wall{};
  ImmersedBoundaryResidualPartsRecord background_contribution;
  ImmersedBoundaryResidualPartsRecord removed_background_contribution;
  ImmersedBoundaryResidualPartsRecord wall_contribution;
  ImmersedBoundaryResidualPartsRecord budget_reaction_delta;
  std::array<double, 3> residual_after_wall{};
  std::uint64_t evaluated_group_count{};
  std::uint64_t simultaneous_substitution_count{};
  std::uint64_t affine_plan_fingerprint{};
  std::vector<ImmersedBoundaryAffineDonorTermRecord>
      background_affine_donor_terms;
  std::vector<ImmersedBoundaryAffineWallGradientTermRecord>
      background_affine_wall_gradient_terms;
  std::vector<ImmersedBoundaryAffineDonorTermRecord> affine_donor_terms;
  std::vector<ImmersedBoundaryAffineWallGradientTermRecord>
      affine_wall_gradient_terms;
  std::uint64_t canonical_affine_row_evaluation_count{};
  std::uint64_t link_local_runtime_evaluation_count{};
  std::uint64_t immutable_input_snapshot_count{};
  std::uint64_t background_functional_evaluation_count{};
  std::uint64_t background_removal_count{};
  std::uint64_t final_row_write_count{};
};

struct ImmersedInterfacePressureDonorTermRecord final {
  mesh::GlobalCellId pressure_cell{};
  std::uint32_t output_component{};
  double coefficient{};
};

struct ImmersedInterfacePressureWallTermRecord final {
  immersed::ImmersedLinkId link{};
  std::uint32_t output_component{};
  double coefficient{};
};

struct ImmersedInterfacePressureRowRecord final {
  mesh::GlobalCellId momentum_cell{};
  std::uint64_t authority_fingerprint{};
  std::vector<ImmersedInterfacePressureDonorTermRecord>
      background_donor_terms;
  std::vector<ImmersedInterfacePressureDonorTermRecord> a22_donor_terms;
  std::vector<ImmersedInterfacePressureDonorTermRecord>
      legacy_unconstrained_lfp_donor_terms;
  std::vector<ImmersedInterfacePressureDonorTermRecord> difference_donor_terms;
  std::vector<ImmersedInterfacePressureWallTermRecord> background_wall_terms;
  std::vector<ImmersedInterfacePressureWallTermRecord> a22_wall_terms;
  std::vector<ImmersedInterfacePressureWallTermRecord>
      legacy_unconstrained_lfp_wall_terms;
  std::vector<ImmersedInterfacePressureWallTermRecord> difference_wall_terms;
};

#endif

class ImmersedBoundaryAuthorityAccess final {
public:
  static const immersed::QuadraticReconstruction &
  row_reconstruction(const ImmersedOperatorAdapter &,
                     immersed::ImmersedLinkId);
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  static std::vector<ImmersedBoundaryRowRecord>
  rows(const ImmersedOperatorAdapter &);
  static std::uint64_t
  last_wall_functional_evaluation_count(const ImmersedOperatorAdapter &);
  static std::uint64_t last_boundary_authority_lookup_probe_count(
      const ImmersedOperatorAdapter &);
  static std::vector<ImmersedBoundaryRowEvaluationRecord>
  last_boundary_row_evaluations(const ImmersedOperatorAdapter &);
  static std::vector<ImmersedInterfacePressureRowRecord>
  interface_pressure_rows(const ImmersedOperatorAdapter &);
  static std::vector<ImmersedInterfacePressureRowRecord>
  interface_pressure_force_rows(const ImmersedOperatorAdapter &);
#endif
};

struct ImmersedWallNormalGradient final {
  std::uint64_t link{};
  double value{};
};

class ForceAuthorityEvaluationScope final {
public:
  ForceAuthorityEvaluationScope();
  ~ForceAuthorityEvaluationScope();

  ForceAuthorityEvaluationScope(const ForceAuthorityEvaluationScope &) =
      delete;
  ForceAuthorityEvaluationScope &
  operator=(const ForceAuthorityEvaluationScope &) = delete;
};

void compute_pressure_gradient_with_wall_normal_constraints(
    const ImmersedReconstruction &, GradientScheme,
    const ReconstructionFieldBinding &cell_values,
    const ReconstructionFieldBinding &cell_gradients,
    const std::vector<ImmersedWallNormalGradient> &);

void accumulate_momentum_with_wall_normal_constraints(
    const ImmersedOperatorAdapter &, const FaceMassFlux &,
    const runtime::FaceFieldView<const double> &face_velocity,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &pressure,
    const std::vector<ImmersedWallNormalGradient> &,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FaceFieldView<const double> &dynamic_viscosity_by_face,
    const runtime::FieldView<double> &residual);

} // namespace hundun::finite_volume::detail
