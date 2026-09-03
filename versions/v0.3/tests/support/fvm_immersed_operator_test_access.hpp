// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#ifndef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
#error "immersed operator test access is unavailable when tests are disabled"
#endif

#include "hundun/ib_surface.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_types.hpp"
#include "fvm_immersed_boundary_authority_detail.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hundun::finite_volume {
class ImmersedOperatorAdapter;
}

namespace hundun::finite_volume::test {

struct ImmersedWallLinkSnapshot final {
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

enum class ImmersedPhysicalTermKind : std::uint8_t {
  neighbour,
  diagonal,
  source,
  convective_direct,
  pressure_direct,
  viscous_orthogonal,
  viscous_deferred_gradient
};

struct ImmersedPhysicalTermSnapshot final {
  std::uint64_t stable_term_id{};
  immersed::ImmersedLinkId link{};
  ImmersedPhysicalTermKind kind{ImmersedPhysicalTermKind::neighbour};
  std::uint32_t algebraic_occurrence{};
  std::uint32_t output_component{};
  double coefficient{};
  std::uint64_t evaluation_group_id{};
  std::vector<std::uint64_t> source_term_ids;
};

struct ImmersedOperatorRowSnapshot final {
  mesh::GlobalCellId active_cell{};
  std::vector<ImmersedWallLinkSnapshot> links;
  std::uint64_t row_replacement_fingerprint{};
  std::uint64_t replacement_group_count{};
  std::vector<ImmersedPhysicalTermSnapshot> covered_physical_terms;
};

enum class BoundaryReplacementTermKind : std::uint8_t {
  pressure_face,
  pressure_diagonal_defect,
  pressure_neighbour_defect,
  viscous_wall,
  viscous_diagonal_defect,
  viscous_neighbour_defect
};

struct BoundaryReplacementTermSnapshot final {
  std::uint64_t stable_term_id{};
  immersed::ImmersedLinkId link{};
  BoundaryReplacementTermKind kind{BoundaryReplacementTermKind::pressure_face};
  std::uint32_t occurrence{};
  std::uint32_t component{};
  double value{};
  std::uint64_t evaluation_group_id{};
};

struct BoundaryResidualPartsSnapshot final {
  std::array<double, 3> convective{};
  std::array<double, 3> pressure{};
  std::array<double, 3> viscous{};
};

enum class BoundaryAffineInputKind : std::uint8_t { pressure, velocity };

struct BoundaryAffineDonorTermSnapshot final {
  std::uint64_t stable_term_id{};
  BoundaryAffineInputKind input_kind{BoundaryAffineInputKind::pressure};
  mesh::GlobalCellId donor_global_cell{};
  std::uint32_t input_component{};
  std::uint32_t output_component{};
  double coefficient{};
  std::uint32_t contributing_link_count{};
  std::vector<std::uint64_t> source_term_ids;
};

struct BoundaryAffineWallGradientTermSnapshot final {
  std::uint64_t stable_term_id{};
  immersed::ImmersedLinkId link{};
  std::uint32_t output_component{};
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

struct BoundaryRowEvaluationSnapshot final {
  mesh::GlobalCellId active_cell{};
  std::uint64_t row_fingerprint{};
  std::vector<BoundaryReplacementTermSnapshot> replacement_terms;
  std::array<double, 3> residual_before_wall{};
  BoundaryResidualPartsSnapshot background_contribution;
  BoundaryResidualPartsSnapshot removed_background_contribution;
  BoundaryResidualPartsSnapshot wall_contribution;
  BoundaryResidualPartsSnapshot budget_reaction_delta;
  std::array<double, 3> residual_after_wall{};
  std::uint64_t evaluated_group_count{};
  std::uint64_t simultaneous_substitution_count{};
  std::uint64_t affine_plan_fingerprint{};
  std::vector<BoundaryAffineDonorTermSnapshot> background_affine_donor_terms;
  std::vector<BoundaryAffineWallGradientTermSnapshot>
      background_affine_wall_gradient_terms;
  std::vector<BoundaryAffineDonorTermSnapshot> affine_donor_terms;
  std::vector<BoundaryAffineWallGradientTermSnapshot>
      affine_wall_gradient_terms;
  std::uint64_t canonical_affine_row_evaluation_count{};
  std::uint64_t link_local_runtime_evaluation_count{};
  std::uint64_t immutable_input_snapshot_count{};
  std::uint64_t background_functional_evaluation_count{};
  std::uint64_t background_removal_count{};
  std::uint64_t final_row_write_count{};
};

struct PressureCouplingDonorTermSnapshot final {
  mesh::GlobalCellId pressure_cell{};
  std::uint32_t output_component{};
  double coefficient{};
};

struct PressureCouplingWallTermSnapshot final {
  immersed::ImmersedLinkId link{};
  std::uint32_t output_component{};
  double coefficient{};
};

struct InterfacePressureRowSnapshot final {
  mesh::GlobalCellId momentum_cell{};
  std::uint64_t authority_fingerprint{};
  std::vector<PressureCouplingDonorTermSnapshot> background_donor_terms;
  std::vector<PressureCouplingDonorTermSnapshot> a22_donor_terms;
  std::vector<PressureCouplingDonorTermSnapshot>
      legacy_unconstrained_lfp_donor_terms;
  std::vector<PressureCouplingDonorTermSnapshot> difference_donor_terms;
  std::vector<PressureCouplingWallTermSnapshot> background_wall_terms;
  std::vector<PressureCouplingWallTermSnapshot> a22_wall_terms;
  std::vector<PressureCouplingWallTermSnapshot>
      legacy_unconstrained_lfp_wall_terms;
  std::vector<PressureCouplingWallTermSnapshot> difference_wall_terms;
};

class ImmersedOperatorTestAccess final {
public:
  static std::vector<ImmersedOperatorRowSnapshot>
  rows(const ImmersedOperatorAdapter &adapter) {
    std::vector<ImmersedOperatorRowSnapshot> result;
    const auto records = detail::ImmersedBoundaryAuthorityAccess::rows(adapter);
    result.reserve(records.size());
    for (const auto &record : records) {
      ImmersedOperatorRowSnapshot snapshot{};
      snapshot.active_cell = record.active_cell;
      snapshot.row_replacement_fingerprint =
          record.row_replacement_fingerprint;
      snapshot.replacement_group_count = record.replacement_group_count;
      snapshot.links.reserve(record.links.size());
      for (const auto &link : record.links)
        snapshot.links.push_back(
            {link.id, link.occurrence, link.normal_scale,
             link.solid_to_fluid_normal, link.replacement_fingerprint,
             link.wall_intercept_m, link.area_from_fluid_m2,
             link.signed_wall_measure_m2, link.pressure_quadrature_m,
             link.surface_measure_m2, link.surface_patch_centroid_m});
      snapshot.covered_physical_terms.reserve(
          record.covered_physical_terms.size());
      for (const auto &term : record.covered_physical_terms)
        snapshot.covered_physical_terms.push_back(
            {term.stable_term_id, term.link, physical_term_kind(term.kind),
             term.algebraic_occurrence, term.output_component,
             term.coefficient, term.evaluation_group_id,
             term.source_term_ids});
      result.push_back(std::move(snapshot));
    }
    return result;
  }
  static std::uint64_t
  last_wall_functional_evaluation_count(
      const ImmersedOperatorAdapter &adapter) {
    return detail::ImmersedBoundaryAuthorityAccess::
        last_wall_functional_evaluation_count(adapter);
  }
  static std::uint64_t
  last_boundary_authority_lookup_probe_count(
      const ImmersedOperatorAdapter &adapter) {
    return detail::ImmersedBoundaryAuthorityAccess::
        last_boundary_authority_lookup_probe_count(adapter);
  }
  static std::vector<BoundaryRowEvaluationSnapshot>
  last_boundary_row_evaluations(const ImmersedOperatorAdapter &adapter) {
    std::vector<BoundaryRowEvaluationSnapshot> result;
    const auto records = detail::ImmersedBoundaryAuthorityAccess::
        last_boundary_row_evaluations(adapter);
    result.reserve(records.size());
    for (const auto &record : records) {
      BoundaryRowEvaluationSnapshot snapshot{};
      snapshot.active_cell = record.active_cell;
      snapshot.row_fingerprint = record.row_fingerprint;
      snapshot.residual_before_wall = record.residual_before_wall;
      snapshot.background_contribution =
          residual_parts(record.background_contribution);
      snapshot.removed_background_contribution =
          residual_parts(record.removed_background_contribution);
      snapshot.wall_contribution = residual_parts(record.wall_contribution);
      snapshot.budget_reaction_delta =
          residual_parts(record.budget_reaction_delta);
      snapshot.residual_after_wall = record.residual_after_wall;
      snapshot.evaluated_group_count = record.evaluated_group_count;
      snapshot.simultaneous_substitution_count =
          record.simultaneous_substitution_count;
      snapshot.affine_plan_fingerprint = record.affine_plan_fingerprint;
      snapshot.canonical_affine_row_evaluation_count =
          record.canonical_affine_row_evaluation_count;
      snapshot.link_local_runtime_evaluation_count =
          record.link_local_runtime_evaluation_count;
      snapshot.immutable_input_snapshot_count =
          record.immutable_input_snapshot_count;
      snapshot.background_functional_evaluation_count =
          record.background_functional_evaluation_count;
      snapshot.background_removal_count = record.background_removal_count;
      snapshot.final_row_write_count = record.final_row_write_count;
      snapshot.replacement_terms.reserve(record.replacement_terms.size());
      for (const auto &term : record.replacement_terms)
        snapshot.replacement_terms.push_back(
            {term.stable_term_id, term.link, replacement_term_kind(term.kind),
             term.occurrence, term.component, term.value,
             term.evaluation_group_id});
      append_affine_donor_terms(record.background_affine_donor_terms,
                                snapshot.background_affine_donor_terms);
      append_affine_wall_terms(
          record.background_affine_wall_gradient_terms,
          snapshot.background_affine_wall_gradient_terms);
      append_affine_donor_terms(record.affine_donor_terms,
                                snapshot.affine_donor_terms);
      append_affine_wall_terms(record.affine_wall_gradient_terms,
                               snapshot.affine_wall_gradient_terms);
      result.push_back(std::move(snapshot));
    }
    return result;
  }
  static std::vector<InterfacePressureRowSnapshot>
  interface_pressure_rows(const ImmersedOperatorAdapter &adapter) {
    return interface_rows(
        detail::ImmersedBoundaryAuthorityAccess::interface_pressure_rows(
            adapter));
  }
  static std::vector<InterfacePressureRowSnapshot>
  interface_pressure_force_rows(const ImmersedOperatorAdapter &adapter) {
    return interface_rows(detail::ImmersedBoundaryAuthorityAccess::
                              interface_pressure_force_rows(adapter));
  }

private:
  static ImmersedPhysicalTermKind physical_term_kind(std::uint8_t kind) {
    switch (kind) {
    case 0U:
      return ImmersedPhysicalTermKind::convective_direct;
    case 1U:
      return ImmersedPhysicalTermKind::pressure_direct;
    case 2U:
      return ImmersedPhysicalTermKind::viscous_orthogonal;
    case 3U:
      return ImmersedPhysicalTermKind::viscous_deferred_gradient;
    }
    throw std::logic_error("immersed operator physical record is invalid");
  }

  static BoundaryReplacementTermKind replacement_term_kind(
      std::uint8_t kind) {
    switch (kind) {
    case 0U:
      return BoundaryReplacementTermKind::pressure_face;
    case 1U:
      return BoundaryReplacementTermKind::pressure_diagonal_defect;
    case 2U:
      return BoundaryReplacementTermKind::pressure_neighbour_defect;
    case 3U:
      return BoundaryReplacementTermKind::viscous_wall;
    case 4U:
      return BoundaryReplacementTermKind::viscous_diagonal_defect;
    case 5U:
      return BoundaryReplacementTermKind::viscous_neighbour_defect;
    }
    throw std::logic_error("immersed operator replacement record is invalid");
  }

  static BoundaryAffineInputKind affine_input_kind(std::uint8_t kind) {
    switch (kind) {
    case 0U:
      return BoundaryAffineInputKind::pressure;
    case 1U:
      return BoundaryAffineInputKind::velocity;
    }
    throw std::logic_error("immersed operator affine record is invalid");
  }

  static BoundaryResidualPartsSnapshot residual_parts(
      const detail::ImmersedBoundaryResidualPartsRecord &record) {
    return {record.convective, record.pressure, record.viscous};
  }

  static void append_affine_donor_terms(
      const std::vector<detail::ImmersedBoundaryAffineDonorTermRecord> &source,
      std::vector<BoundaryAffineDonorTermSnapshot> &target) {
    target.reserve(source.size());
    for (const auto &term : source)
      target.push_back(
          {term.stable_term_id, affine_input_kind(term.input_kind),
           term.donor_global_cell, term.input_component, term.output_component,
           term.coefficient, term.contributing_link_count,
           term.source_term_ids});
  }

  static void append_affine_wall_terms(
      const std::vector<detail::ImmersedBoundaryAffineWallGradientTermRecord>
          &source,
      std::vector<BoundaryAffineWallGradientTermSnapshot> &target) {
    target.reserve(source.size());
    for (const auto &term : source)
      target.push_back({term.stable_term_id, term.link, term.output_component,
                        term.coefficient, term.source_term_ids});
  }

  static std::vector<InterfacePressureRowSnapshot> interface_rows(
      const std::vector<detail::ImmersedInterfacePressureRowRecord> &records) {
    std::vector<InterfacePressureRowSnapshot> result;
    result.reserve(records.size());
    for (const auto &record : records) {
      InterfacePressureRowSnapshot snapshot{};
      snapshot.momentum_cell = record.momentum_cell;
      snapshot.authority_fingerprint = record.authority_fingerprint;
      append_pressure_donor_terms(record.background_donor_terms,
                                  snapshot.background_donor_terms);
      append_pressure_donor_terms(record.a22_donor_terms,
                                  snapshot.a22_donor_terms);
      append_pressure_donor_terms(
          record.legacy_unconstrained_lfp_donor_terms,
          snapshot.legacy_unconstrained_lfp_donor_terms);
      append_pressure_donor_terms(record.difference_donor_terms,
                                  snapshot.difference_donor_terms);
      append_pressure_wall_terms(record.background_wall_terms,
                                 snapshot.background_wall_terms);
      append_pressure_wall_terms(record.a22_wall_terms,
                                 snapshot.a22_wall_terms);
      append_pressure_wall_terms(
          record.legacy_unconstrained_lfp_wall_terms,
          snapshot.legacy_unconstrained_lfp_wall_terms);
      append_pressure_wall_terms(record.difference_wall_terms,
                                 snapshot.difference_wall_terms);
      result.push_back(std::move(snapshot));
    }
    return result;
  }

  static void append_pressure_donor_terms(
      const std::vector<detail::ImmersedInterfacePressureDonorTermRecord>
          &source,
      std::vector<PressureCouplingDonorTermSnapshot> &target) {
    target.reserve(source.size());
    for (const auto &term : source)
      target.push_back(
          {term.pressure_cell, term.output_component, term.coefficient});
  }

  static void append_pressure_wall_terms(
      const std::vector<detail::ImmersedInterfacePressureWallTermRecord>
          &source,
      std::vector<PressureCouplingWallTermSnapshot> &target) {
    target.reserve(source.size());
    for (const auto &term : source)
      target.push_back({term.link, term.output_component, term.coefficient});
  }
};

} // namespace hundun::finite_volume::test
