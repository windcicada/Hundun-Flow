// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_immersed_static.hpp"

#include "hundun/diag_immersed_module.hpp"

#include "hundun/ib_domain.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/rt_mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>

namespace hundun::diagnostics {
namespace {

constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);
constexpr std::string_view kInstanceId = "primary";

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
    hash ^= (value >> (byte * 8U)) & UINT64_C(0xff);
    hash *= prime;
  }
  return hash;
}

runtime::Real3 add(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

runtime::Real3 scale(double factor, runtime::Real3 value) noexcept {
  return {factor * value.x, factor * value.y, factor * value.z};
}

runtime::Real3 cross(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool static_kind(DiagnosticModuleKind kind) noexcept {
  return kind == DiagnosticModuleKind::immersed_surface ||
         kind == DiagnosticModuleKind::ghost_stencil ||
         kind == DiagnosticModuleKind::local_flow_pattern;
}

std::string_view module_id(DiagnosticModuleKind kind) {
  switch (kind) {
  case DiagnosticModuleKind::immersed_surface:
    return "hundun.immersed.surface";
  case DiagnosticModuleKind::ghost_stencil:
    return "hundun.immersed.ghost-stencil";
  case DiagnosticModuleKind::local_flow_pattern:
    return "hundun.immersed.local-flow-pattern";
  default:
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "stage3.immersed-static.diagnostics.kind", -1,
        "Diagnostic kind is not an immersed static provider");
  }
}

std::string_view phase_name(DiagnosticModuleKind kind) {
  switch (kind) {
  case DiagnosticModuleKind::immersed_surface:
    return "immersed-static.surface";
  case DiagnosticModuleKind::ghost_stencil:
    return "immersed-static.ghost-stencil";
  case DiagnosticModuleKind::local_flow_pattern:
    return "immersed-static.local-flow-pattern";
  default:
    return {};
  }
}

void require_summary(const ImmersedStaticDiagnosticSummary &value) {
  const bool empty_links = value.immersed_link_count == 0U;
  const bool valid_link_evidence =
      empty_links
          ? value.donor_reference_count == 0U &&
                value.minimum_reconstruction_rank == 0U &&
                value.maximum_reconstruction_rank == 0U &&
                bits(value.maximum_reconstruction_condition) == bits(0.0)
          : value.donor_reference_count > 0U &&
                value.minimum_reconstruction_rank > 0U &&
                value.maximum_reconstruction_rank >=
                    value.minimum_reconstruction_rank &&
                value.maximum_reconstruction_condition > 0.0 &&
                std::isfinite(value.maximum_reconstruction_condition);
  if (value.vertex_count == 0U || value.triangle_count == 0U ||
      value.connected_component_count == 0U ||
      !finite(value.bounding_box_min_m) ||
      !finite(value.bounding_box_max_m) ||
      !(value.surface_area_m2 > 0.0) ||
      !(value.closed_volume_m3 > 0.0) ||
      !std::isfinite(value.surface_area_m2) ||
      !std::isfinite(value.closed_volume_m3) ||
      !finite(value.area_vector_closure_m2) || value.orientation == 0 ||
      value.surface_fingerprint == 0U || value.query_fingerprint == 0U ||
      value.classification_fingerprint == 0U ||
      value.surface_coverage_fingerprint == 0U || !valid_link_evidence ||
      value.maximum_halo_reach == 0U ||
      value.covered_triangle_count > value.triangle_count ||
      value.ghost_plan_fingerprint == 0U ||
      value.wall_plan_fingerprint == 0U ||
      value.local_flow_pattern_algorithm_fingerprint == 0U) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input,
        "stage3.immersed-static.summary.input", -1,
        "Stage 3 immersed static summary is invalid");
  }
}

void require_local_flow_pattern_snapshot(
    const ImmersedStaticDiagnosticSummary &value) {
  if (!value.local_flow_pattern_snapshot_available ||
      value.local_flow_pattern_row_fingerprint == 0U ||
      value.replacement_coefficient_l2 < 0.0 ||
      !std::isfinite(value.replacement_coefficient_l2) ||
      value.limiting_case_status == 0U) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::unavailable,
        "stage3.local-flow-pattern.diagnostics.unavailable", -1,
        "Accepted local-flow-pattern snapshot is unavailable");
  }
}

void require_request(const ImmersedStaticDiagnosticSummary &source,
                     DiagnosticModuleKind kind,
                     const DiagnosticRequest &request,
                     DiagnosticScope scope) {
  try {
    require_summary(source);
    if (kind == DiagnosticModuleKind::local_flow_pattern)
      require_local_flow_pattern_snapshot(source);
    if (!static_kind(kind))
      static_cast<void>(module_id(kind));
    validate(request, describe_diagnostics(source, kind));
    if (request.scope != scope || request.frame.phase != phase_name(kind) ||
        request.frame.step != 0U || bits(request.frame.time_s) != bits(0.0) ||
        !request.selected_fields.empty() || request.sample_budget != 0U ||
        (request.level != DiagnosticLevel::summary &&
         request.level != DiagnosticLevel::counters)) {
      throw std::runtime_error("immersed static request mismatch");
    }
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "stage3.immersed-static.diagnostics.request", -1,
        "Stage 3 immersed static diagnostic request is invalid");
  }
}

void add_u64(DiagnosticFingerprintAccumulator &accumulator,
             std::string_view field, std::uint64_t value) {
  accumulator.add(field, 0U, 0U,
                  describe_fp64(static_cast<double>(value >> 32U)));
  accumulator.add(field, 0U, 1U,
                  describe_fp64(static_cast<double>(
                      static_cast<std::uint32_t>(value))));
}

void add_real3(DiagnosticFingerprintAccumulator &accumulator,
               std::string_view field, runtime::Real3 value) {
  accumulator.add(field, 0U, 0U, describe_fp64(value.x));
  accumulator.add(field, 0U, 1U, describe_fp64(value.y));
  accumulator.add(field, 0U, 2U, describe_fp64(value.z));
}

DiagnosticStateFingerprint fingerprint(
    const ImmersedStaticDiagnosticSummary &source, DiagnosticModuleKind kind) {
  DiagnosticFingerprintAccumulator accumulator;
  if (kind == DiagnosticModuleKind::immersed_surface) {
    add_u64(accumulator, "vertices", source.vertex_count);
    add_u64(accumulator, "triangles", source.triangle_count);
    add_u64(accumulator, "components", source.connected_component_count);
    add_real3(accumulator, "bbox", source.bounding_box_min_m);
    add_real3(accumulator, "bbox", source.bounding_box_max_m);
    accumulator.add("area", 0U, 0U, describe_fp64(source.surface_area_m2));
    accumulator.add("closed-volume", 0U, 0U,
                    describe_fp64(source.closed_volume_m3));
    add_real3(accumulator, "area-vector-closure",
              source.area_vector_closure_m2);
    add_u64(accumulator, "orientation",
            source.orientation > 0 ? 1U : 0U);
    add_u64(accumulator, "surface-fingerprint",
            source.surface_fingerprint);
  } else if (kind == DiagnosticModuleKind::ghost_stencil) {
    add_u64(accumulator, "links", source.immersed_link_count);
    add_u64(accumulator, "donors", source.donor_reference_count);
    add_u64(accumulator, "qr-rank", source.minimum_reconstruction_rank);
    add_u64(accumulator, "qr-rank", source.maximum_reconstruction_rank);
    accumulator.add("condition", 0U, 0U,
                    describe_fp64(source.maximum_reconstruction_condition));
    add_u64(accumulator, "halo-reach", source.maximum_halo_reach);
    add_u64(accumulator, "wall-points",
            source.wall_quadrature_point_count);
    add_u64(accumulator, "triangle-coverage",
            source.covered_triangle_count);
    add_u64(accumulator, "ghost-plan-fingerprint",
            source.ghost_plan_fingerprint);
    add_u64(accumulator, "wall-plan-fingerprint",
            source.wall_plan_fingerprint);
  } else {
    add_u64(accumulator, "algorithm-fingerprint",
            source.local_flow_pattern_algorithm_fingerprint);
    add_u64(accumulator, "row-fingerprint",
            source.local_flow_pattern_row_fingerprint);
    add_u64(accumulator, "replacement-groups",
            source.replacement_group_count);
    add_u64(accumulator, "occurrences",
            source.algebraic_occurrence_count);
    accumulator.add("coefficient-norm", 0U, 0U,
                    describe_fp64(source.replacement_coefficient_l2));
    add_u64(accumulator, "limiting-case", source.limiting_case_status);
  }
  return accumulator.finish();
}

void add_xyz_metrics(std::vector<DiagnosticMetric> &metrics,
                     std::string_view prefix, std::string_view unit,
                     runtime::Real3 value) {
  metrics.push_back({std::string(prefix) + ".x",
                     DiagnosticMetricKind::state_summary, std::string(unit),
                     describe_fp64(value.x)});
  metrics.push_back({std::string(prefix) + ".y",
                     DiagnosticMetricKind::state_summary, std::string(unit),
                     describe_fp64(value.y)});
  metrics.push_back({std::string(prefix) + ".z",
                     DiagnosticMetricKind::state_summary, std::string(unit),
                     describe_fp64(value.z)});
}

DiagnosticRecord build_record(const ImmersedStaticDiagnosticSummary &source,
                              DiagnosticModuleKind kind,
                              const DiagnosticRequest &request) {
  DiagnosticRecord record;
  record.schema_version = kDiagnosticRecordSchemaV1;
  record.module_kind = kind;
  record.module_id = std::string(module_id(kind));
  record.instance_id = std::string(kInstanceId);
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = 0U;
  record.time_s = describe_fp64(0.0);
  record.phase = std::string(phase_name(kind));
  record.status = DiagnosticStatus::ok;
  record.failure = {};
  record.state_fingerprint = fingerprint(source, kind);

  if (kind == DiagnosticModuleKind::immersed_surface) {
    record.identities.push_back(
        {"immersed-domain-classification", std::nullopt,
         source.classification_fingerprint, std::nullopt, std::nullopt});
    record.identities.push_back(
        {"immersed-domain-surface-coverage", std::nullopt,
         source.surface_coverage_fingerprint, std::nullopt, std::nullopt});
    record.identities.push_back(
        {"immersed-surface", std::nullopt, source.surface_fingerprint,
         std::nullopt, std::nullopt});
    record.identities.push_back(
        {"immersed-surface-query", std::nullopt, source.query_fingerprint,
         std::nullopt, std::nullopt});
    if (request.level == DiagnosticLevel::summary) {
      add_xyz_metrics(record.metrics, "area-vector-closure", "m2",
                      source.area_vector_closure_m2);
      add_xyz_metrics(record.metrics, "bounding-box.maximum", "m",
                      source.bounding_box_max_m);
      add_xyz_metrics(record.metrics, "bounding-box.minimum", "m",
                      source.bounding_box_min_m);
      record.metrics.push_back(
          {"closed-volume", DiagnosticMetricKind::state_summary, "m3",
           describe_fp64(source.closed_volume_m3)});
      record.metrics.push_back(
          {"surface-area", DiagnosticMetricKind::state_summary, "m2",
           describe_fp64(source.surface_area_m2)});
    } else {
      record.counters = {
          {"connected-components", "count", source.connected_component_count},
          {"orientation", "count",
           source.orientation > 0 ? 1U : 0U},
          {"surface-fingerprint", "count", source.surface_fingerprint},
          {"triangles", "count", source.triangle_count},
          {"vertices", "count", source.vertex_count},
      };
    }
  } else if (kind == DiagnosticModuleKind::ghost_stencil) {
    record.identities = {
        {"ghost-stencil", std::nullopt, source.ghost_plan_fingerprint,
         std::nullopt, std::nullopt},
        {"immersed-domain-classification", std::nullopt,
         source.classification_fingerprint, std::nullopt, std::nullopt},
        {"immersed-domain-surface-coverage", std::nullopt,
         source.surface_coverage_fingerprint, std::nullopt, std::nullopt},
        {"wall-quadrature", std::nullopt, source.wall_plan_fingerprint,
         std::nullopt, std::nullopt},
    };
    if (request.level == DiagnosticLevel::summary) {
      record.metrics.push_back(
          {"reconstruction.maximum-condition",
           DiagnosticMetricKind::state_summary, "1",
           describe_fp64(source.maximum_reconstruction_condition)});
    } else {
      record.counters = {
          {"donor-references", "count", source.donor_reference_count},
          {"ghost-plan-fingerprint", "count", source.ghost_plan_fingerprint},
          {"immersed-links", "count", source.immersed_link_count},
          {"maximum-halo-reach", "cell", source.maximum_halo_reach},
          {"reconstruction.maximum-rank", "count",
           source.maximum_reconstruction_rank},
          {"reconstruction.minimum-rank", "count",
           source.minimum_reconstruction_rank},
          {"triangle-coverage", "count", source.covered_triangle_count},
          {"wall-plan-fingerprint", "count", source.wall_plan_fingerprint},
          {"wall-quadrature-points", "count",
           source.wall_quadrature_point_count},
      };
    }
  } else {
    record.identities.push_back(
        {"local-flow-pattern", std::nullopt,
         source.local_flow_pattern_algorithm_fingerprint, std::nullopt,
         std::nullopt});
    if (request.level == DiagnosticLevel::summary) {
      record.metrics.push_back(
          {"replacement-coefficient-l2",
           DiagnosticMetricKind::state_summary, "1",
           describe_fp64(source.replacement_coefficient_l2)});
    } else {
      record.counters = {
          {"algebraic-occurrences", "count",
           source.algebraic_occurrence_count},
          {"algorithm-fingerprint", "count",
           source.local_flow_pattern_algorithm_fingerprint},
          {"limiting-case-status", "count", source.limiting_case_status},
          {"replacement-groups", "count", source.replacement_group_count},
          {"row-fingerprint", "count",
           source.local_flow_pattern_row_fingerprint},
      };
    }
  }
  return record;
}

void submit(DiagnosticSink &sink, const DiagnosticRecord &record,
            int failing_rank = -1) {
  try {
    sink.submit(record);
  } catch (...) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", failing_rank,
                                    "Diagnostic sink rejected a record");
  }
}

std::uint64_t collective_fingerprint(const runtime::MpiContext &mpi,
                                     std::uint64_t local,
                                     std::string_view operation) {
  std::array<std::uint64_t, 2> input{local, local};
  std::array<std::uint64_t, 2> reduced{};
  runtime::check_mpi_result(
      MPI_Allreduce(input.data(), reduced.data(), 1, MPI_UINT64_T, MPI_BXOR,
                    mpi.comm()),
      std::string(operation) + " xor");
  runtime::check_mpi_result(
      MPI_Allreduce(input.data() + 1, reduced.data() + 1, 1, MPI_UINT64_T,
                    MPI_SUM, mpi.comm()),
      std::string(operation) + " sum");
  std::uint64_t result = UINT64_C(1469598103934665603);
  result = mix(result, reduced[0]);
  result = mix(result, reduced[1]);
  result = mix(result, static_cast<std::uint64_t>(mpi.size()));
  return result == 0U ? 1U : result;
}

ImmersedStaticDiagnosticSummary collective_summary(
    const ImmersedStaticDiagnosticSummary &local, DiagnosticModuleKind kind,
    const runtime::MpiContext &mpi) {
  auto result = local;
  if (kind == DiagnosticModuleKind::immersed_surface) {
    const std::array<std::uint64_t, 17> identity{
        local.vertex_count,
        local.triangle_count,
        local.connected_component_count,
        bits(local.bounding_box_min_m.x),
        bits(local.bounding_box_min_m.y),
        bits(local.bounding_box_min_m.z),
        bits(local.bounding_box_max_m.x),
        bits(local.bounding_box_max_m.y),
        bits(local.bounding_box_max_m.z),
        bits(local.surface_area_m2),
        bits(local.closed_volume_m3),
        bits(local.area_vector_closure_m2.x),
        bits(local.area_vector_closure_m2.y),
        bits(local.area_vector_closure_m2.z),
        local.orientation > 0 ? 1U : 0U,
        local.surface_fingerprint,
        local.query_fingerprint,
    };
    std::array<std::uint64_t, 17> minimum{};
    std::array<std::uint64_t, 17> maximum{};
    runtime::check_mpi_result(
        MPI_Allreduce(identity.data(), minimum.data(),
                      static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MIN,
                      mpi.comm()),
        "MPI_Allreduce(Stage 3 surface identity minimum)");
    runtime::check_mpi_result(
        MPI_Allreduce(identity.data(), maximum.data(),
                      static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MAX,
                      mpi.comm()),
        "MPI_Allreduce(Stage 3 surface identity maximum)");
    if (minimum != maximum)
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::layout,
          "stage3.immersed-static.surface.agreement", 0,
          "Collective immersed-surface summaries disagree");
    result.classification_fingerprint = collective_fingerprint(
        mpi, local.classification_fingerprint,
        "MPI_Allreduce(Stage 3 classification fingerprint)");
    result.surface_coverage_fingerprint = collective_fingerprint(
        mpi, local.surface_coverage_fingerprint,
        "MPI_Allreduce(Stage 3 surface-coverage fingerprint)");
    return result;
  }

  if (kind == DiagnosticModuleKind::ghost_stencil) {
    const std::array<std::uint64_t, 4> local_counts{
        local.immersed_link_count, local.donor_reference_count,
        local.wall_quadrature_point_count, local.covered_triangle_count};
    std::array<std::uint64_t, 4> global_counts{};
    runtime::check_mpi_result(
        MPI_Allreduce(local_counts.data(), global_counts.data(),
                      static_cast<int>(local_counts.size()), MPI_UINT64_T,
                      MPI_SUM, mpi.comm()),
        "MPI_Allreduce(Stage 3 ghost counts)");
    result.immersed_link_count = global_counts[0];
    result.donor_reference_count = global_counts[1];
    result.wall_quadrature_point_count = global_counts[2];
    result.covered_triangle_count = global_counts[3];
    const std::uint64_t local_minimum_rank =
        local.immersed_link_count == 0U
            ? std::numeric_limits<std::uint64_t>::max()
            : local.minimum_reconstruction_rank;
    runtime::check_mpi_result(
        MPI_Allreduce(&local_minimum_rank,
                      &result.minimum_reconstruction_rank, 1, MPI_UINT64_T,
                      MPI_MIN, mpi.comm()),
        "MPI_Allreduce(Stage 3 ghost minimum rank)");
    if (result.immersed_link_count == 0U)
      result.minimum_reconstruction_rank = 0U;
    const std::array<std::uint64_t, 2> local_maximum{
        local.maximum_reconstruction_rank, local.maximum_halo_reach};
    std::array<std::uint64_t, 2> global_maximum{};
    runtime::check_mpi_result(
        MPI_Allreduce(local_maximum.data(), global_maximum.data(),
                      static_cast<int>(local_maximum.size()), MPI_UINT64_T,
                      MPI_MAX, mpi.comm()),
        "MPI_Allreduce(Stage 3 ghost maxima)");
    result.maximum_reconstruction_rank = global_maximum[0];
    result.maximum_halo_reach = global_maximum[1];
    runtime::check_mpi_result(
        MPI_Allreduce(&local.maximum_reconstruction_condition,
                      &result.maximum_reconstruction_condition, 1, MPI_DOUBLE,
                      MPI_MAX, mpi.comm()),
        "MPI_Allreduce(Stage 3 ghost condition)");
    result.ghost_plan_fingerprint = collective_fingerprint(
        mpi, local.ghost_plan_fingerprint,
        "MPI_Allreduce(Stage 3 ghost fingerprint)");
    result.wall_plan_fingerprint = collective_fingerprint(
        mpi, local.wall_plan_fingerprint,
        "MPI_Allreduce(Stage 3 wall-plan fingerprint)");
    result.classification_fingerprint = collective_fingerprint(
        mpi, local.classification_fingerprint,
        "MPI_Allreduce(Stage 3 classification fingerprint)");
    result.surface_coverage_fingerprint = collective_fingerprint(
        mpi, local.surface_coverage_fingerprint,
        "MPI_Allreduce(Stage 3 surface-coverage fingerprint)");
    return result;
  }

  const std::array<std::uint64_t, 2> algorithm{
      local.local_flow_pattern_algorithm_fingerprint,
      local.local_flow_pattern_algorithm_fingerprint};
  std::array<std::uint64_t, 2> algorithm_bounds{};
  runtime::check_mpi_result(
      MPI_Allreduce(algorithm.data(), algorithm_bounds.data(), 1,
                    MPI_UINT64_T, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Stage 3 LFP algorithm minimum)");
  runtime::check_mpi_result(
      MPI_Allreduce(algorithm.data() + 1, algorithm_bounds.data() + 1, 1,
                    MPI_UINT64_T, MPI_MAX, mpi.comm()),
      "MPI_Allreduce(Stage 3 LFP algorithm maximum)");
  if (algorithm_bounds[0] != algorithm_bounds[1])
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "stage3.local-flow-pattern.diagnostics.agreement", 0,
        "Collective local-flow-pattern algorithms disagree");
  const std::array<std::uint64_t, 2> local_counts{
      local.replacement_group_count, local.algebraic_occurrence_count};
  std::array<std::uint64_t, 2> global_counts{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_counts.data(), global_counts.data(),
                    static_cast<int>(local_counts.size()), MPI_UINT64_T,
                    MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Stage 3 LFP counts)");
  result.replacement_group_count = global_counts[0];
  result.algebraic_occurrence_count = global_counts[1];
  const double local_square = local.replacement_coefficient_l2 *
                              local.replacement_coefficient_l2;
  double global_square{};
  runtime::check_mpi_result(
      MPI_Allreduce(&local_square, &global_square, 1, MPI_DOUBLE, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(Stage 3 LFP coefficient norm)");
  result.replacement_coefficient_l2 = std::sqrt(global_square);
  runtime::check_mpi_result(
      MPI_Allreduce(&local.limiting_case_status,
                    &result.limiting_case_status, 1, MPI_UINT64_T, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(Stage 3 LFP limiting status)");
  result.local_flow_pattern_row_fingerprint = collective_fingerprint(
      mpi, local.local_flow_pattern_row_fingerprint,
      "MPI_Allreduce(Stage 3 LFP row fingerprint)");
  return result;
}

} // namespace

ImmersedStaticDiagnosticSummary summarize_immersed_static(
    const immersed::ImmersedSurface &surface,
    const immersed::SurfaceQuery &query,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost,
    const immersed::WallQuadraturePlan &wall,
    const immersed::LocalFlowPatternTransform &transform) {
  ImmersedStaticDiagnosticSummary result;
  result.vertex_count = surface.vertex_count();
  result.triangle_count = surface.triangle_count();
  result.bounding_box_min_m = surface.bounding_box_min_m();
  result.bounding_box_max_m = surface.bounding_box_max_m();
  result.closed_volume_m3 = surface.closed_volume_m3();
  result.surface_fingerprint = surface.fingerprint();
  result.query_fingerprint = query.fingerprint();
  result.classification_fingerprint = domain.classification_fingerprint();
  result.surface_coverage_fingerprint =
      domain.surface_coverage_fingerprint();
  result.ghost_plan_fingerprint = ghost.fingerprint();
  result.wall_plan_fingerprint = wall.fingerprint();
  result.local_flow_pattern_algorithm_fingerprint =
      transform.algorithm_fingerprint();

  using VertexKey = std::array<std::uint64_t, 3>;
  std::map<VertexKey, std::size_t> first_triangle_by_vertex;
  std::vector<std::size_t> parent(result.triangle_count);
  for (std::size_t triangle = 0U; triangle < parent.size(); ++triangle)
    parent[triangle] = triangle;
  const auto find = [&](std::size_t value) {
    while (parent[value] != value) {
      parent[value] = parent[parent[value]];
      value = parent[value];
    }
    return value;
  };
  const auto unite = [&](std::size_t left, std::size_t right) {
    left = find(left);
    right = find(right);
    if (left != right)
      parent[right] = left;
  };
  double signed_volume = 0.0;
  for (std::size_t index = 0U; index < result.triangle_count; ++index) {
    const auto &triangle = surface.triangle(index);
    result.surface_area_m2 += triangle.area_m2;
    result.area_vector_closure_m2 =
        add(result.area_vector_closure_m2,
            scale(triangle.area_m2, triangle.geometric_outward_normal));
    signed_volume += dot(triangle.vertices_m[0],
                         cross(triangle.vertices_m[1],
                               triangle.vertices_m[2])) /
                     6.0;
    for (const auto vertex : triangle.vertices_m) {
      const VertexKey key{bits(vertex.x), bits(vertex.y), bits(vertex.z)};
      const auto [position, inserted] =
          first_triangle_by_vertex.emplace(key, index);
      if (!inserted)
        unite(index, position->second);
    }
  }
  std::set<std::size_t> components;
  for (std::size_t index = 0U; index < parent.size(); ++index)
    components.insert(find(index));
  result.connected_component_count = components.size();
  result.orientation = signed_volume > 0.0 ? 1 : signed_volume < 0.0 ? -1 : 0;

  result.immersed_link_count = domain.links().size();
  result.minimum_reconstruction_rank =
      std::numeric_limits<std::uint64_t>::max();
  std::set<immersed::TriangleId> covered_triangles;
  for (const auto &link : domain.links()) {
    const auto add_donors = [&](const auto &constraint) {
      result.donor_reference_count += constraint.donors.size();
    };
    for (std::size_t component = 0U; component < 3U; ++component)
      add_donors(ghost.velocity_constraint(link.id, component));
    add_donors(ghost.zero_normal_constraint(link.id));
    add_donors(ghost.density_extrapolation(link.id));
    const auto &quality = ghost.reconstruction(link.id).quality();
    result.minimum_reconstruction_rank =
        std::min(result.minimum_reconstruction_rank,
                 static_cast<std::uint64_t>(quality.rank));
    result.maximum_reconstruction_rank =
        std::max(result.maximum_reconstruction_rank,
                 static_cast<std::uint64_t>(quality.rank));
    result.maximum_reconstruction_condition =
        std::max(result.maximum_reconstruction_condition,
                 quality.condition_estimate);
    result.maximum_halo_reach =
        std::max(result.maximum_halo_reach,
                 static_cast<std::uint64_t>(quality.halo_reach));
    covered_triangles.insert(link.triangle);
  }
  result.maximum_halo_reach =
      std::max({result.maximum_halo_reach,
                static_cast<std::uint64_t>(ghost.maximum_halo_reach()),
                static_cast<std::uint64_t>(wall.maximum_halo_reach())});
  result.wall_quadrature_point_count = wall.local_points().size();
  for (const auto &point : wall.local_points())
    covered_triangles.insert(point.triangle);
  result.covered_triangle_count = covered_triangles.size();
  if (result.immersed_link_count == 0U)
    result.minimum_reconstruction_rank = 0U;
  require_summary(result);
  return result;
}

ImmersedStaticDiagnosticSummary with_local_flow_pattern_snapshot(
    ImmersedStaticDiagnosticSummary result,
    const flow::ImmersedFlowDiagnosticSource &source) {
  require_summary(result);
  if (!source.local_flow_pattern_available() || source.snapshot_seal() == 0U ||
      result.local_flow_pattern_algorithm_fingerprint !=
          source.local_flow_pattern_algorithm_fingerprint()) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::unavailable,
        "stage3.local-flow-pattern.diagnostics.unavailable", -1,
        "Accepted local-flow-pattern snapshot is unavailable");
  }
  result.local_flow_pattern_snapshot_available = true;
  result.local_flow_pattern_row_fingerprint =
      source.local_flow_pattern_row_fingerprint();
  result.replacement_group_count =
      source.local_flow_pattern_replacement_group_count();
  result.algebraic_occurrence_count =
      source.local_flow_pattern_algebraic_occurrence_count();
  result.replacement_coefficient_l2 =
      source.local_flow_pattern_replacement_coefficient_l2();
  result.limiting_case_status =
      source.local_flow_pattern_limiting_case_status();
  require_local_flow_pattern_snapshot(result);
  return result;
}

DiagnosticDescriptor describe_diagnostics(
    const ImmersedStaticDiagnosticSummary &, DiagnosticModuleKind kind) {
  return {kDiagnosticRecordSchemaV1, kind, module_id(kind), kInstanceId,
          kCapabilities};
}

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const ImmersedStaticDiagnosticSummary &, DiagnosticModuleKind kind) {
  switch (kind) {
  case DiagnosticModuleKind::immersed_surface:
    return {"area", "area-vector-closure", "bbox", "closed-volume",
            "components", "orientation", "surface-fingerprint", "triangles",
            "vertices"};
  case DiagnosticModuleKind::ghost_stencil:
    return {"condition", "donors", "ghost-plan-fingerprint", "halo-reach",
            "links", "qr-rank", "triangle-coverage",
            "wall-plan-fingerprint", "wall-points"};
  case DiagnosticModuleKind::local_flow_pattern:
    return {"algorithm-fingerprint", "coefficient-norm", "limiting-case",
            "occurrences", "replacement-groups", "row-fingerprint"};
  default:
    static_cast<void>(module_id(kind));
  }
  return {};
}

void collect_diagnostics(const ImmersedStaticDiagnosticSummary &source,
                         DiagnosticModuleKind kind,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  require_request(source, kind, request, DiagnosticScope::local);
  const auto record = build_record(source, kind, request);
  validate(record, describe_diagnostics(source, kind), request);
  submit(sink, record);
}

void collect_diagnostics(const ImmersedStaticDiagnosticSummary &source,
                         DiagnosticModuleKind kind,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  bool ready = true;
  try {
    require_request(source, kind, request, DiagnosticScope::collective);
  } catch (...) {
    ready = false;
  }
  const int candidate = ready ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Stage 3 static diagnostic preflight)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "stage3.immersed-static.diagnostics.preflight", lowest,
        "Collective Stage 3 static diagnostic preflight failed");

  std::optional<DiagnosticRecord> record;
  bool record_failed = false;
  try {
    const auto aggregate = collective_summary(source, kind, mpi);
    record.emplace(build_record(aggregate, kind, request));
    validate(*record, describe_diagnostics(aggregate, kind), request);
  } catch (...) {
    record_failed = true;
  }
  const int record_candidate = record_failed ? mpi.rank() : mpi.size();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&record_candidate, &lowest, 1, MPI_INT, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(Stage 3 static diagnostic record)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "stage3.immersed-static.diagnostics.record", lowest,
        "Collective Stage 3 static diagnostic record failed");
  bool sink_failed = false;
  try {
    sink.submit(*record);
  } catch (...) {
    sink_failed = true;
  }
  const int sink_candidate = sink_failed ? mpi.rank() : mpi.size();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&sink_candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Stage 3 static diagnostic sink)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", lowest,
                                    "Collective diagnostic sink failed");
}

} // namespace hundun::diagnostics
