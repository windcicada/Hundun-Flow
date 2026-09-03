// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_immersed_static.hpp"
#include "hundun/diag_checkpoint_v3.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace hundun;

class RecordingSink final : public diagnostics::DiagnosticSink {
public:
  void submit(const diagnostics::DiagnosticRecord &record) override {
    records.push_back(record);
  }
  std::vector<diagnostics::DiagnosticRecord> records;
};

diagnostics::ImmersedStaticDiagnosticSummary summary() {
  diagnostics::ImmersedStaticDiagnosticSummary result;
  result.vertex_count = 8U;
  result.triangle_count = 12U;
  result.connected_component_count = 1U;
  result.bounding_box_min_m = {-0.5, -0.5, -0.5};
  result.bounding_box_max_m = {0.5, 0.5, 0.5};
  result.surface_area_m2 = 6.0;
  result.closed_volume_m3 = 1.0;
  result.area_vector_closure_m2 = {0.0, 0.0, 0.0};
  result.orientation = 1;
  result.surface_fingerprint = UINT64_C(0x101);
  result.query_fingerprint = UINT64_C(0x102);
  result.classification_fingerprint = UINT64_C(0x103);
  result.surface_coverage_fingerprint = UINT64_C(0x104);
  result.immersed_link_count = 48U;
  result.donor_reference_count = 384U;
  result.minimum_reconstruction_rank = 6U;
  result.maximum_reconstruction_rank = 6U;
  result.maximum_reconstruction_condition = 17.0;
  result.maximum_halo_reach = 4U;
  result.wall_quadrature_point_count = 36U;
  result.covered_triangle_count = 12U;
  result.ghost_plan_fingerprint = UINT64_C(0x105);
  result.wall_plan_fingerprint = UINT64_C(0x106);
  result.local_flow_pattern_snapshot_available = true;
  result.local_flow_pattern_algorithm_fingerprint = UINT64_C(0x107);
  result.local_flow_pattern_row_fingerprint = UINT64_C(0x108);
  result.replacement_group_count = 24U;
  result.algebraic_occurrence_count = 72U;
  result.replacement_coefficient_l2 = 3.25;
  result.limiting_case_status = 1U;
  return result;
}

std::vector<diagnostics::DiagnosticModuleKind>
expected_inventory(std::uint8_t profile) {
  using K = diagnostics::DiagnosticModuleKind;
  const bool ibm = profile == 1U || profile == 3U || profile == 4U ||
                   profile == 6U || profile == 7U || profile == 9U;
  const bool wale = profile == 2U || profile == 3U || profile == 5U ||
                    profile == 6U || profile == 8U || profile == 9U;
  std::vector<K> result;
  if (ibm) {
    result = {K::immersed_surface, K::ghost_stencil,
              K::local_flow_pattern, K::wall_force};
  }
  if (wale)
    result.push_back(K::les);
  return result;
}

void check_inventory() {
  for (std::uint8_t value = 1U; value <= 9U; ++value) {
    const auto actual = diagnostics::stage3_added_provider_inventory(
        static_cast<flow::CheckpointV3Presence>(value));
    HUNDUN_CHECK(actual == expected_inventory(value));
  }
  HUNDUN_CHECK(diagnostics::stage3_added_provider_inventory(
                   static_cast<flow::CheckpointV3Presence>(0U))
                   .empty());
  HUNDUN_CHECK(diagnostics::stage3_added_provider_inventory(
                   static_cast<flow::CheckpointV3Presence>(10U))
                   .empty());
  const flow::CheckpointV3Report default_report;
  HUNDUN_CHECK(diagnostics::stage3_added_provider_inventory(default_report) ==
               expected_inventory(1U));
}

diagnostics::DiagnosticRequest request(diagnostics::DiagnosticLevel level,
                                       std::string_view phase) {
  diagnostics::DiagnosticRequest result;
  result.level = level;
  result.scope = diagnostics::DiagnosticScope::local;
  result.frame = {0, 0U, 0.0, phase};
  return result;
}

std::vector<std::string> ids(const diagnostics::DiagnosticRecord &record) {
  std::vector<std::string> result;
  for (const auto &value : record.metrics)
    result.push_back(value.id + ":" + value.unit);
  for (const auto &value : record.counters)
    result.push_back(value.id + ":" + value.unit);
  return result;
}

void check_static_provider_contract() {
  using K = diagnostics::DiagnosticModuleKind;
  const auto source = summary();
  struct Expected final {
    K kind;
    std::string_view module;
    std::string_view phase;
    std::vector<std::string_view> fingerprint_fields;
    std::vector<std::string> summary_fields;
    std::vector<std::string> counter_fields;
  };
  const std::array<Expected, 3> expected{{
      {K::immersed_surface,
       "hundun.immersed.surface",
       "immersed-static.surface",
       {"area", "area-vector-closure", "bbox", "closed-volume",
        "components", "orientation", "surface-fingerprint", "triangles",
        "vertices"},
       {"area-vector-closure.x:m2", "area-vector-closure.y:m2",
        "area-vector-closure.z:m2", "bounding-box.maximum.x:m",
        "bounding-box.maximum.y:m", "bounding-box.maximum.z:m",
        "bounding-box.minimum.x:m", "bounding-box.minimum.y:m",
        "bounding-box.minimum.z:m", "closed-volume:m3", "surface-area:m2"},
       {"connected-components:count", "orientation:count",
        "surface-fingerprint:count", "triangles:count", "vertices:count"}},
      {K::ghost_stencil,
       "hundun.immersed.ghost-stencil",
       "immersed-static.ghost-stencil",
       {"condition", "donors", "ghost-plan-fingerprint", "halo-reach",
        "links", "qr-rank", "triangle-coverage",
        "wall-plan-fingerprint", "wall-points"},
       {"reconstruction.maximum-condition:1"},
       {"donor-references:count", "ghost-plan-fingerprint:count",
        "immersed-links:count", "maximum-halo-reach:cell",
        "reconstruction.maximum-rank:count",
        "reconstruction.minimum-rank:count",
        "triangle-coverage:count", "wall-plan-fingerprint:count",
        "wall-quadrature-points:count"}},
      {K::local_flow_pattern,
       "hundun.immersed.local-flow-pattern",
       "immersed-static.local-flow-pattern",
       {"algorithm-fingerprint", "coefficient-norm", "limiting-case",
        "occurrences", "replacement-groups", "row-fingerprint"},
       {"replacement-coefficient-l2:1"},
       {"algebraic-occurrences:count", "algorithm-fingerprint:count",
        "limiting-case-status:count", "replacement-groups:count",
        "row-fingerprint:count"}},
  }};

  for (const auto &entry : expected) {
    const auto descriptor = diagnostics::describe_diagnostics(source, entry.kind);
    HUNDUN_CHECK(descriptor.module_kind == entry.kind);
    HUNDUN_CHECK(descriptor.module_id == entry.module);
    HUNDUN_CHECK(descriptor.instance_id == "primary");
    diagnostics::validate(descriptor);
    HUNDUN_CHECK(diagnostics::diagnostic_fingerprint_field_ids(source,
                                                               entry.kind) ==
                 entry.fingerprint_fields);

    RecordingSink first;
    diagnostics::collect_diagnostics(
        source, entry.kind,
        request(diagnostics::DiagnosticLevel::summary, entry.phase), first);
    HUNDUN_CHECK(first.records.size() == 1U);
    HUNDUN_CHECK(ids(first.records.front()) == entry.summary_fields);
    RecordingSink second;
    diagnostics::collect_diagnostics(
        source, entry.kind,
        request(diagnostics::DiagnosticLevel::summary, entry.phase), second);
    HUNDUN_CHECK(diagnostics::to_canonical_json(first.records.front()) ==
                 diagnostics::to_canonical_json(second.records.front()));

    RecordingSink counters;
    diagnostics::collect_diagnostics(
        source, entry.kind,
        request(diagnostics::DiagnosticLevel::counters, entry.phase),
        counters);
    HUNDUN_CHECK(counters.records.size() == 1U);
    HUNDUN_CHECK(ids(counters.records.front()) == entry.counter_fields);
  }

  bool rejected = false;
  try {
    static_cast<void>(diagnostics::describe_diagnostics(
        source, diagnostics::DiagnosticModuleKind::les));
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    rejected = error.classification() ==
                   diagnostics::DiagnosticFailureClass::capability &&
               error.code() == "stage3.immersed-static.diagnostics.kind";
  }
  HUNDUN_CHECK(rejected);

  auto empty_local = source;
  empty_local.immersed_link_count = 0U;
  empty_local.donor_reference_count = 0U;
  empty_local.minimum_reconstruction_rank = 0U;
  empty_local.maximum_reconstruction_rank = 0U;
  empty_local.maximum_reconstruction_condition = 0.0;
  empty_local.wall_quadrature_point_count = 0U;
  empty_local.covered_triangle_count = 0U;
  empty_local.replacement_group_count = 0U;
  empty_local.algebraic_occurrence_count = 0U;
  empty_local.replacement_coefficient_l2 = 0.0;
  empty_local.limiting_case_status = 1U;
  for (const auto kind : {K::ghost_stencil, K::local_flow_pattern}) {
    const auto phase = kind == K::ghost_stencil
                           ? "immersed-static.ghost-stencil"
                           : "immersed-static.local-flow-pattern";
    RecordingSink sink;
    diagnostics::collect_diagnostics(
        empty_local, kind,
        request(diagnostics::DiagnosticLevel::summary, phase), sink);
    HUNDUN_CHECK(sink.records.size() == 1U);
  }
}

} // namespace

int main() {
  return hundun::test::run([] {
    check_inventory();
    check_static_provider_contract();
  });
}
