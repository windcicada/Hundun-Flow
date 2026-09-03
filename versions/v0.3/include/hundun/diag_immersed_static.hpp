// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_types.hpp"

#include <cstdint>
#include <string_view>
#include <vector>

namespace hundun::immersed {
class GhostStencilPlan;
class ImmersedDomain;
class ImmersedSurface;
class LocalFlowPatternTransform;
class SurfaceQuery;
class WallQuadraturePlan;
} // namespace hundun::immersed

namespace hundun::flow {
class ImmersedFlowDiagnosticSource;
}

namespace hundun::diagnostics {

// Owns only value evidence captured from immutable Stage 3 construction
// authorities. It does not borrow any numerical object.
struct ImmersedStaticDiagnosticSummary final {
  std::uint64_t vertex_count{};
  std::uint64_t triangle_count{};
  std::uint64_t connected_component_count{};
  runtime::Real3 bounding_box_min_m{};
  runtime::Real3 bounding_box_max_m{};
  double surface_area_m2{};
  double closed_volume_m3{};
  runtime::Real3 area_vector_closure_m2{};
  std::int64_t orientation{};
  std::uint64_t surface_fingerprint{};
  std::uint64_t query_fingerprint{};
  std::uint64_t classification_fingerprint{};
  std::uint64_t surface_coverage_fingerprint{};

  std::uint64_t immersed_link_count{};
  std::uint64_t donor_reference_count{};
  std::uint64_t minimum_reconstruction_rank{};
  std::uint64_t maximum_reconstruction_rank{};
  double maximum_reconstruction_condition{};
  std::uint64_t maximum_halo_reach{};
  std::uint64_t wall_quadrature_point_count{};
  std::uint64_t covered_triangle_count{};
  std::uint64_t ghost_plan_fingerprint{};
  std::uint64_t wall_plan_fingerprint{};

  bool local_flow_pattern_snapshot_available{};
  std::uint64_t local_flow_pattern_algorithm_fingerprint{};
  std::uint64_t local_flow_pattern_row_fingerprint{};
  std::uint64_t replacement_group_count{};
  std::uint64_t algebraic_occurrence_count{};
  double replacement_coefficient_l2{};
  std::uint64_t limiting_case_status{};
};

ImmersedStaticDiagnosticSummary summarize_immersed_static(
    const immersed::ImmersedSurface &, const immersed::SurfaceQuery &,
    const immersed::ImmersedDomain &, const immersed::GhostStencilPlan &,
    const immersed::WallQuadraturePlan &,
    const immersed::LocalFlowPatternTransform &);
ImmersedStaticDiagnosticSummary with_local_flow_pattern_snapshot(
    ImmersedStaticDiagnosticSummary,
    const flow::ImmersedFlowDiagnosticSource &);

std::vector<DiagnosticModuleKind>
stage3_added_provider_inventory(flow::CheckpointV3Presence);

DiagnosticDescriptor describe_diagnostics(
    const ImmersedStaticDiagnosticSummary &, DiagnosticModuleKind);
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const ImmersedStaticDiagnosticSummary &, DiagnosticModuleKind);
void collect_diagnostics(const ImmersedStaticDiagnosticSummary &,
                         DiagnosticModuleKind, const DiagnosticRequest &,
                         DiagnosticSink &);
void collect_diagnostics(const ImmersedStaticDiagnosticSummary &,
                         DiagnosticModuleKind, const runtime::MpiContext &,
                         const DiagnosticRequest &, DiagnosticSink &);

} // namespace hundun::diagnostics
