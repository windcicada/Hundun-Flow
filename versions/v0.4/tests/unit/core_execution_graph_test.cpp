// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_execution.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using hundun::v04::ExecutionGraphCompiler;
using hundun::v04::FieldAccessSpec;
using hundun::v04::FrozenExecutionGraph;
using hundun::v04::FrozenStage;
using hundun::v04::GraphEdge;
using hundun::v04::GraphFieldSpec;
using hundun::v04::GraphNode;
using hundun::v04::GraphNodeKind;
using hundun::v04::ResourceContract;
using hundun::v04::Span;
using hundun::v04::StageKind;
using hundun::v04::StageSpec;
using hundun::v04::StateVisibility;
using hundun::v04::Status;
using hundun::v04::StatusCode;

constexpr FieldAccessSpec accepted(std::uint16_t field) noexcept {
  return {field, StateVisibility::accepted};
}

constexpr FieldAccessSpec trial(std::uint16_t field) noexcept {
  return {field, StateVisibility::trial};
}

constexpr FieldAccessSpec committed(std::uint16_t field) noexcept {
  return {field, StateVisibility::committed_snapshot};
}

template <class T, std::size_t Size>
Span<const T> span(const std::array<T, Size>& values) noexcept {
  return {values.data(), values.size()};
}

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool expect_ok(Status status, std::string_view description) {
  return expect(static_cast<bool>(status), description);
}

bool expect_invalid(Status status, std::string_view description) {
  return expect(status.code == StatusCode::invalid_plan, description);
}

bool graph_is_pristine(const FrozenExecutionGraph& graph) noexcept {
  const ResourceContract& resources = graph.resources();
  return graph.fingerprint() == 0U && graph.stages().size == 0U &&
         graph.nodes().size == 0U && graph.edges().size == 0U &&
         resources.max_live_workspace_bytes == 0U &&
         resources.allocation_allowance == 0U &&
         resources.merged_halo_messages == 0U &&
         resources.merged_halo_bytes == 0U &&
         resources.numeric_refills == 0U &&
         resources.hierarchy_rebuilds == 0U &&
         resources.cache_publishes == 0U &&
         resources.linear_iterations == 0U &&
         resources.stage_wall_nanoseconds == 0U;
}

bool same_stage(const FrozenStage& left, const FrozenStage& right) noexcept {
  return left.id == right.id && left.kind == right.kind &&
         left.registration_ordinal == right.registration_ordinal &&
         left.node_begin == right.node_begin &&
         left.node_count == right.node_count &&
         left.read_begin == right.read_begin &&
         left.read_count == right.read_count &&
         left.write_begin == right.write_begin &&
         left.write_count == right.write_count &&
         left.ghost_begin == right.ghost_begin &&
         left.ghost_count == right.ghost_count &&
         left.invalidation_begin == right.invalidation_begin &&
         left.invalidation_count == right.invalidation_count &&
         left.workspace_offset == right.workspace_offset &&
         left.workspace_bytes == right.workspace_bytes &&
         left.workspace_alignment == right.workspace_alignment &&
         left.workspace_live_through == right.workspace_live_through &&
         left.resources.merged_halo_messages ==
             right.resources.merged_halo_messages &&
         left.resources.merged_halo_bytes ==
             right.resources.merged_halo_bytes &&
         left.resources.numeric_refills == right.resources.numeric_refills &&
         left.resources.hierarchy_rebuilds ==
             right.resources.hierarchy_rebuilds &&
         left.resources.cache_publishes ==
             right.resources.cache_publishes &&
         left.resources.linear_iterations ==
             right.resources.linear_iterations &&
         left.resources.stage_wall_nanoseconds ==
             right.resources.stage_wall_nanoseconds &&
         left.collective_consensus == right.collective_consensus;
}

bool same_node(GraphNode left, GraphNode right) noexcept {
  return left.stage == right.stage && left.kind == right.kind &&
         left.ordinal == right.ordinal;
}

bool same_resources(const ResourceContract& left,
                    const ResourceContract& right) noexcept {
  return left.max_live_workspace_bytes == right.max_live_workspace_bytes &&
         left.allocation_allowance == right.allocation_allowance &&
         left.merged_halo_messages == right.merged_halo_messages &&
         left.merged_halo_bytes == right.merged_halo_bytes &&
         left.numeric_refills == right.numeric_refills &&
         left.hierarchy_rebuilds == right.hierarchy_rebuilds &&
         left.cache_publishes == right.cache_publishes &&
         left.linear_iterations == right.linear_iterations &&
         left.stage_wall_nanoseconds == right.stage_wall_nanoseconds;
}

bool same_accesses(Span<const FieldAccessSpec> left,
                   Span<const FieldAccessSpec> right) noexcept {
  if (left.size != right.size) {
    return false;
  }
  for (std::size_t index = 0U; index < left.size; ++index) {
    if (!(left.data[index] == right.data[index])) {
      return false;
    }
  }
  return true;
}

bool same_graph(const FrozenExecutionGraph& left,
                const FrozenExecutionGraph& right) noexcept {
  if (left.fingerprint() != right.fingerprint() ||
      !same_resources(left.resources(), right.resources()) ||
      left.stages().size != right.stages().size ||
      left.nodes().size != right.nodes().size ||
      left.edges().size != right.edges().size) {
    return false;
  }
  for (std::size_t index = 0U; index < left.stages().size; ++index) {
    const FrozenStage& lhs = left.stages().data[index];
    const FrozenStage& rhs = right.stages().data[index];
    if (!same_stage(lhs, rhs) ||
        !same_accesses(left.reads(lhs.id), right.reads(rhs.id)) ||
        !same_accesses(left.writes(lhs.id), right.writes(rhs.id)) ||
        !same_accesses(left.ghosts(lhs.id), right.ghosts(rhs.id)) ||
        !same_accesses(left.invalidations(lhs.id),
                       right.invalidations(rhs.id))) {
      return false;
    }
    const Span<const std::uint8_t> lhs_widths = left.ghost_widths(lhs.id);
    const Span<const std::uint8_t> rhs_widths = right.ghost_widths(rhs.id);
    if (lhs_widths.size != rhs_widths.size) {
      return false;
    }
    for (std::size_t width = 0U; width < lhs_widths.size; ++width) {
      if (lhs_widths.data[width] != rhs_widths.data[width]) {
        return false;
      }
    }
  }
  for (std::size_t index = 0U; index < left.nodes().size; ++index) {
    if (!same_node(left.nodes().data[index], right.nodes().data[index])) {
      return false;
    }
  }
  for (std::size_t index = 0U; index < left.edges().size; ++index) {
    if (!(left.edges().data[index] == right.edges().data[index])) {
      return false;
    }
  }
  return true;
}

bool build_reference_graph(FrozenExecutionGraph& out) {
  const std::array<GraphFieldSpec, 3U> fields{{
      {accepted(1U), 2U, true},
      {trial(2U), 0U, false},
      {committed(3U), 0U, true},
  }};
  ExecutionGraphCompiler compiler;
  if (!compiler.configure(span(fields))) {
    return false;
  }

  const std::array<FieldAccessSpec, 1U> numerical_reads{{accepted(1U)}};
  const std::array<FieldAccessSpec, 1U> numerical_writes{{trial(2U)}};
  const std::array<FieldAccessSpec, 1U> numerical_ghosts{{accepted(1U)}};
  const std::array<std::uint8_t, 1U> numerical_widths{{2U}};
  StageSpec numerical;
  numerical.id = 10U;
  numerical.reads = span(numerical_reads);
  numerical.writes = span(numerical_writes);
  numerical.ghosts = span(numerical_ghosts);
  numerical.ghost_widths = span(numerical_widths);
  numerical.workspace_bytes = 128U;
  numerical.workspace_alignment = 64U;
  numerical.resources = {3U, 1024U, 2U, 1U, 4U};
  numerical.collective_consensus = true;
  if (!compiler.register_stage(numerical)) {
    return false;
  }

  const std::array<FieldAccessSpec, 1U> service_reads{{committed(3U)}};
  StageSpec service;
  service.id = 20U;
  service.reads = span(service_reads);
  service.kind = StageKind::service;
  service.collective_consensus = true;
  if (!compiler.register_stage(service)) {
    return false;
  }

  StageSpec commit;
  commit.id = 30U;
  commit.kind = StageKind::commit;
  return compiler.register_stage(commit) && compiler.freeze_for_test(out);
}

bool test_happy_graph_has_exact_order_and_edges() {
  FrozenExecutionGraph graph;
  bool passed = true;
  passed &= expect(build_reference_graph(graph),
                   "reference synthetic graph compiles");
  if (!passed) {
    return false;
  }

  constexpr std::array<GraphNode, 8U> expected_nodes{{
      {10U, GraphNodeKind::halo_begin, 0U},
      {10U, GraphNodeKind::compute_interior, 0U},
      {10U, GraphNodeKind::halo_finish, 0U},
      {10U, GraphNodeKind::compute_boundary, 0U},
      {10U, GraphNodeKind::collective_consensus, 0U},
      {20U, GraphNodeKind::service, 1U},
      {20U, GraphNodeKind::collective_consensus, 1U},
      {30U, GraphNodeKind::commit, 2U},
  }};
  constexpr std::array<GraphEdge, 7U> expected_edges{{
      {0U, 1U}, {1U, 2U}, {2U, 3U},
      {3U, 4U}, {4U, 5U}, {5U, 6U}, {6U, 7U},
  }};
  passed &= expect(graph.stages().size == 3U,
                   "declared stage count is preserved");
  passed &= expect(graph.nodes().size == expected_nodes.size(),
                   "halo, consensus, service, and commit node count is exact");
  passed &= expect(graph.edges().size == expected_edges.size(),
                   "declared serial dependency edge count is exact");
  if (graph.nodes().size == expected_nodes.size()) {
    for (std::size_t index = 0U; index < expected_nodes.size(); ++index) {
      passed &= expect(same_node(graph.nodes().data[index],
                                 expected_nodes[index]),
                       "node order and stage ordinal are deterministic");
    }
  }
  if (graph.edges().size == expected_edges.size()) {
    for (std::size_t index = 0U; index < expected_edges.size(); ++index) {
      passed &= expect(graph.edges().data[index] == expected_edges[index],
                       "dependency edge endpoints are exact");
    }
  }
  const FrozenStage* numerical = graph.stage(10U);
  passed &= expect(numerical != nullptr, "stage lookup finds numerical stage");
  if (numerical != nullptr) {
    passed &= expect(graph.reads(numerical->id).size == 1U &&
                         graph.reads(numerical->id).data[0] == accepted(1U),
                     "frozen numerical reads are exact");
    passed &= expect(graph.writes(numerical->id).size == 1U &&
                         graph.writes(numerical->id).data[0] == trial(2U),
                     "frozen numerical writes are exact");
    passed &= expect(graph.ghost_widths(numerical->id).size == 1U &&
                         graph.ghost_widths(numerical->id).data[0] == 2U,
                     "frozen ghost width is exact");
  }
  const ResourceContract& resources = graph.resources();
  passed &= expect(resources.max_live_workspace_bytes == 128U &&
                       resources.allocation_allowance == 0U &&
                       resources.merged_halo_messages == 3U &&
                       resources.merged_halo_bytes == 1024U &&
                       resources.numeric_refills == 2U &&
                       resources.hierarchy_rebuilds == 1U &&
                       resources.cache_publishes == 4U,
                   "frozen graph resource contract is exact");
  passed &= expect(graph.fingerprint() != 0U,
                   "successful graph has a nonzero fingerprint");
  return passed;
}

bool test_reads_and_ghosts_require_production() {
  bool passed = true;
  {
    const std::array<GraphFieldSpec, 1U> fields{{
        {trial(1U), 1U, false},
    }};
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "read-before-produce fixture configures");
    const std::array<FieldAccessSpec, 1U> reads{{trial(1U)}};
    StageSpec stage;
    stage.id = 1U;
    stage.reads = span(reads);
    passed &= expect_ok(compiler.register_stage(stage),
                        "availability is checked at freeze, not registration");
    FrozenExecutionGraph out;
    passed &= expect_invalid(compiler.freeze_for_test(out),
                             "read-before-produce is rejected");
    passed &= expect(graph_is_pristine(out),
                     "read-before-produce failure leaves output pristine");
  }
  {
    const std::array<GraphFieldSpec, 1U> fields{{
        {trial(1U), 1U, false},
    }};
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "missing ghost production fixture configures");
    const std::array<FieldAccessSpec, 1U> reads{{trial(1U)}};
    const std::array<FieldAccessSpec, 1U> ghosts{{trial(1U)}};
    const std::array<std::uint8_t, 1U> widths{{1U}};
    StageSpec stage;
    stage.id = 1U;
    stage.reads = span(reads);
    stage.ghosts = span(ghosts);
    stage.ghost_widths = span(widths);
    passed &= expect_ok(compiler.register_stage(stage),
                        "ghost stage shape registers before liveness check");
    FrozenExecutionGraph out;
    passed &= expect_invalid(compiler.freeze_for_test(out),
                             "ghost read without a producer is rejected");
  }
  return passed;
}

bool test_writer_version_rules() {
  bool passed = true;
  {
    const std::array<GraphFieldSpec, 1U> fields{{
        {trial(1U), 0U, false},
    }};
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "two-writer fixture configures");
    const std::array<FieldAccessSpec, 1U> writes{{trial(1U)}};
    StageSpec first;
    first.id = 1U;
    first.writes = span(writes);
    first.collective_consensus = true;
    StageSpec second = first;
    second.id = 2U;
    passed &= expect_ok(compiler.register_stage(first),
                        "first writer registers");
    passed &= expect_ok(compiler.register_stage(second),
                        "second writer registers for whole-graph validation");
    FrozenExecutionGraph out;
    passed &= expect_invalid(
        compiler.freeze_for_test(out),
        "second writer without read and invalidation is rejected");
  }
  {
    const std::array<GraphFieldSpec, 1U> fields{{
        {trial(1U), 0U, false},
    }};
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "explicit version-chain fixture configures");
    const std::array<FieldAccessSpec, 1U> access{{trial(1U)}};
    StageSpec first;
    first.id = 1U;
    first.writes = span(access);
    first.collective_consensus = true;
    StageSpec second;
    second.id = 2U;
    second.reads = span(access);
    second.writes = span(access);
    second.invalidates = span(access);
    second.collective_consensus = true;
    passed &= expect_ok(compiler.register_stage(first),
                        "initial version writer registers");
    passed &= expect_ok(compiler.register_stage(second),
                        "read-invalidate-write replacement registers");
    FrozenExecutionGraph out;
    passed &= expect_ok(compiler.freeze_for_test(out),
                        "explicit read-invalidate-write version chain compiles");
    const FrozenStage* replacement = out.stage(2U);
    passed &= expect(replacement != nullptr &&
                         out.invalidations(replacement->id).size == 1U &&
                         out.invalidations(replacement->id).data[0] == trial(1U),
                     "replacement invalidation is retained in frozen graph");
  }
  return passed;
}

bool test_ghost_contract_validation() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {accepted(1U), 1U, true},
  }};
  const std::array<FieldAccessSpec, 1U> access{{accepted(1U)}};
  bool passed = true;
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "ghost-not-read fixture configures");
    const std::array<std::uint8_t, 1U> widths{{1U}};
    StageSpec stage;
    stage.id = 1U;
    stage.ghosts = span(access);
    stage.ghost_widths = span(widths);
    passed &= expect_invalid(compiler.register_stage(stage),
                             "ghost access must also be declared as a read");
  }
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "ghost-width mismatch fixture configures");
    StageSpec stage;
    stage.id = 1U;
    stage.reads = span(access);
    stage.ghosts = span(access);
    passed &= expect_invalid(compiler.register_stage(stage),
                             "ghost and width span sizes must match");
  }
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "ghost-capacity fixture configures");
    const std::array<std::uint8_t, 1U> widths{{2U}};
    StageSpec stage;
    stage.id = 1U;
    stage.reads = span(access);
    stage.ghosts = span(access);
    stage.ghost_widths = span(widths);
    passed &= expect_invalid(compiler.register_stage(stage),
                             "ghost width above field capacity is rejected");
  }
  return passed;
}

bool test_invalidation_contract_validation() {
  bool passed = true;
  {
    const std::array<GraphFieldSpec, 1U> fields{{
        {accepted(1U), 0U, true},
    }};
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "undeclared-invalidation fixture configures");
    const std::array<FieldAccessSpec, 1U> undeclared{{accepted(2U)}};
    StageSpec stage;
    stage.id = 1U;
    stage.invalidates = span(undeclared);
    passed &= expect_invalid(compiler.register_stage(stage),
                             "invalidation of an undeclared field is rejected");
  }
  {
    const std::array<GraphFieldSpec, 1U> fields{{
        {accepted(1U), 0U, true},
    }};
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "invalid-invalidation fixture configures");
    const std::array<FieldAccessSpec, 1U> access{{accepted(1U)}};
    StageSpec stage;
    stage.id = 1U;
    stage.invalidates = span(access);
    passed &= expect_ok(compiler.register_stage(stage),
                        "declared invalidation defers version validation");
    FrozenExecutionGraph out;
    passed &= expect_invalid(
        compiler.freeze_for_test(out),
        "invalidation without matching read and replacement write is rejected");
  }
  return passed;
}

bool test_service_reads_only_committed_snapshot() {
  const std::array<GraphFieldSpec, 2U> fields{{
      {accepted(1U), 0U, true},
      {committed(1U), 0U, true},
  }};
  bool passed = true;
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "trial-service fixture configures");
    const std::array<FieldAccessSpec, 1U> reads{{accepted(1U)}};
    StageSpec service;
    service.id = 1U;
    service.reads = span(reads);
    service.kind = StageKind::service;
    service.collective_consensus = true;
    passed &= expect_invalid(
        compiler.register_stage(service),
        "service cannot read accepted/trial transaction state");
  }
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "committed-service fixture configures");
    const std::array<FieldAccessSpec, 1U> reads{{committed(1U)}};
    StageSpec service;
    service.id = 1U;
    service.reads = span(reads);
    service.kind = StageKind::service;
    service.resources.stage_wall_nanoseconds = 17U;
    service.collective_consensus = true;
    passed &= expect_ok(compiler.register_stage(service),
                        "service may read committed snapshot");
    FrozenExecutionGraph out;
    passed &= expect_ok(compiler.freeze_for_test(out),
                        "committed-snapshot service graph freezes");
    passed &= expect(
        out.nodes().size == 2U &&
            out.nodes().data[0].kind == GraphNodeKind::service &&
            out.nodes().data[1].kind == GraphNodeKind::collective_consensus &&
            out.resources().stage_wall_nanoseconds == 17U,
        "service compiles to a service and consensus node with timing budget");
  }
  return passed;
}

bool test_service_resource_isolation() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {committed(1U), 0U, true},
  }};
  const std::array<FieldAccessSpec, 1U> reads{{committed(1U)}};
  bool passed = true;
  for (std::size_t resource = 0U; resource < 6U; ++resource) {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "service-resource fixture configures");
    StageSpec service;
    service.id = 1U;
    service.kind = StageKind::service;
    service.reads = span(reads);
    service.collective_consensus = true;
    switch (resource) {
      case 0U:
        service.resources.merged_halo_messages = 1U;
        break;
      case 1U:
        service.resources.merged_halo_bytes = 1U;
        break;
      case 2U:
        service.resources.numeric_refills = 1U;
        break;
      case 3U:
        service.resources.hierarchy_rebuilds = 1U;
        break;
      case 4U:
        service.resources.cache_publishes = 1U;
        break;
      case 5U:
        service.resources.linear_iterations = 1U;
        break;
      default:
        break;
    }
    passed &= expect_invalid(
        compiler.register_stage(service),
        "committed-snapshot service cannot claim numerical mutation resources");
  }
  return passed;
}

bool test_zero_length_nonnull_spans_are_empty() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {committed(1U), 0U, true},
  }};
  const FieldAccessSpec unused = committed(1U);
  const Span<const FieldAccessSpec> empty{&unused, 0U};
  ExecutionGraphCompiler compiler;
  bool passed = expect_ok(compiler.configure(span(fields)),
                          "nonnull-empty-span fixture configures");
  StageSpec stage;
  stage.id = 1U;
  stage.reads = empty;
  stage.writes = empty;
  stage.ghosts = empty;
  stage.ghost_widths = {reinterpret_cast<const std::uint8_t*>(&unused), 0U};
  stage.invalidates = empty;
  stage.collective_consensus = true;
  passed &= expect_ok(compiler.register_stage(stage),
                      "zero-length spans are empty regardless of pointer value");
  FrozenExecutionGraph out;
  passed &= expect_ok(compiler.freeze_for_test(out),
                      "nonnull zero-length spans freeze as empty partitions");
  return passed;
}

bool test_commit_requires_prior_collective_consensus() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {trial(1U), 0U, false},
  }};
  ExecutionGraphCompiler compiler;
  bool passed = true;
  passed &= expect_ok(compiler.configure(span(fields)),
                      "commit-consensus fixture configures");
  const std::array<FieldAccessSpec, 1U> writes{{trial(1U)}};
  StageSpec writer;
  writer.id = 1U;
  writer.writes = span(writes);
  StageSpec commit;
  commit.id = 2U;
  commit.kind = StageKind::commit;
  passed &= expect_ok(compiler.register_stage(writer), "writer registers");
  passed &= expect_ok(compiler.register_stage(commit), "commit registers");
  FrozenExecutionGraph out;
  passed &= expect_invalid(compiler.freeze_for_test(out),
                           "commit before collective consensus is rejected");
  passed &= expect(graph_is_pristine(out),
                   "failed commit compilation does not modify output");
  passed &= expect(!compiler.frozen(),
                   "failed commit compilation does not freeze compiler");
  return passed;
}

bool test_commit_requires_consensus_after_last_work() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {committed(1U), 0U, true},
  }};
  const std::array<FieldAccessSpec, 1U> committed_reads{{committed(1U)}};
  bool passed = true;
  for (const StageKind trailing_kind :
       {StageKind::compute, StageKind::service}) {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "post-consensus-work fixture configures");

    StageSpec agreed;
    agreed.id = 1U;
    agreed.collective_consensus = true;
    StageSpec trailing;
    trailing.id = 2U;
    trailing.kind = trailing_kind;
    if (trailing_kind == StageKind::service) {
      trailing.reads = span(committed_reads);
    }
    StageSpec commit_stage;
    commit_stage.id = 3U;
    commit_stage.kind = StageKind::commit;

    passed &= expect_ok(compiler.register_stage(agreed),
                        "initially agreed stage registers");
    passed &= expect_ok(compiler.register_stage(trailing),
                        "unconsensed trailing work registers cold");
    passed &= expect_ok(compiler.register_stage(commit_stage),
                        "commit registers cold");
    FrozenExecutionGraph out;
    passed &= expect_invalid(
        compiler.freeze_for_test(out),
        "commit rejects compute or service work after the last consensus");
    passed &= expect(graph_is_pristine(out),
                     "post-consensus-work rejection is atomic");
  }
  return passed;
}

bool test_commit_is_a_pure_transition() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {committed(1U), 0U, true},
  }};
  const std::array<FieldAccessSpec, 1U> reads{{committed(1U)}};
  bool passed = true;
  for (std::size_t mutation = 0U; mutation < 13U; ++mutation) {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "pure-commit fixture configures");
    StageSpec commit_stage;
    commit_stage.id = 1U;
    commit_stage.kind = StageKind::commit;
    switch (mutation) {
      case 0U:
        commit_stage.reads = span(reads);
        break;
      case 1U:
        commit_stage.resources.merged_halo_messages = 1U;
        break;
      case 2U:
        commit_stage.resources.merged_halo_bytes = 1U;
        break;
      case 3U:
        commit_stage.resources.numeric_refills = 1U;
        break;
      case 4U:
        commit_stage.resources.hierarchy_rebuilds = 1U;
        break;
      case 5U:
        commit_stage.resources.cache_publishes = 1U;
        break;
      case 6U:
        commit_stage.resources.linear_iterations = 1U;
        break;
      case 7U:
        commit_stage.resources.stage_wall_nanoseconds = 1U;
        break;
      case 8U:
        commit_stage.workspace_bytes = 64U;
        break;
      case 9U:
        commit_stage.workspace_live_through = 2U;
        break;
      case 10U:
        commit_stage.has_fixed_workspace_offset = true;
        break;
      case 11U:
        commit_stage.fixed_workspace_offset = 64U;
        break;
      case 12U:
        commit_stage.workspace_alignment = 128U;
        break;
      default:
        break;
    }
    passed &= expect_invalid(
        compiler.register_stage(commit_stage),
        "commit cannot read fields or claim numerical resources");
  }
  return passed;
}

bool test_registration_spans_are_deep_copied() {
  std::array<GraphFieldSpec, 2U> fields{{
      {accepted(1U), 1U, true},
      {trial(2U), 0U, false},
  }};
  ExecutionGraphCompiler compiler;
  bool passed = true;
  passed &= expect_ok(compiler.configure(span(fields)),
                      "deep-copy fixture configures");
  fields[0] = {accepted(99U), 0U, false};
  fields[1] = {trial(98U), 0U, true};

  std::array<FieldAccessSpec, 1U> reads{{accepted(1U)}};
  std::array<FieldAccessSpec, 1U> writes{{trial(2U)}};
  std::array<FieldAccessSpec, 1U> ghosts{{accepted(1U)}};
  std::array<std::uint8_t, 1U> widths{{1U}};
  StageSpec stage;
  stage.id = 1U;
  stage.reads = span(reads);
  stage.writes = span(writes);
  stage.ghosts = span(ghosts);
  stage.ghost_widths = span(widths);
  stage.collective_consensus = true;
  passed &= expect_ok(compiler.register_stage(stage),
                      "stage registers against copied field declarations");
  reads[0] = accepted(99U);
  writes[0] = trial(98U);
  ghosts[0] = accepted(99U);
  widths[0] = 7U;

  FrozenExecutionGraph out;
  passed &= expect_ok(compiler.freeze_for_test(out),
                      "mutating caller spans cannot change registered stage");
  const FrozenStage* frozen = out.stage(1U);
  passed &= expect(frozen != nullptr, "deep-copied stage remains addressable");
  if (frozen != nullptr) {
    passed &= expect(out.reads(frozen->id).data[0] == accepted(1U) &&
                         out.writes(frozen->id).data[0] == trial(2U) &&
                         out.ghosts(frozen->id).data[0] == accepted(1U) &&
                         out.ghost_widths(frozen->id).data[0] == 1U,
                     "field, stage, and ghost-width spans are owned copies");
  }
  return passed;
}

bool test_failed_freeze_is_atomic_and_recoverable() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {trial(1U), 0U, false},
  }};
  ExecutionGraphCompiler compiler;
  bool passed = true;
  passed &= expect_ok(compiler.configure(span(fields)),
                      "recoverable-freeze fixture configures");
  const std::array<FieldAccessSpec, 1U> writes{{trial(1U)}};
  StageSpec writer;
  writer.id = 1U;
  writer.writes = span(writes);
  passed &= expect_ok(compiler.register_stage(writer),
                      "unconsensed writer registers");
  FrozenExecutionGraph out;
  passed &= expect_invalid(compiler.freeze_for_test(out),
                           "trailing unconsensed mutation fails freeze");
  passed &= expect(graph_is_pristine(out),
                   "failed freeze leaves every output property unchanged");
  passed &= expect(!compiler.frozen(),
                   "failed freeze leaves compiler open for completion");

  StageSpec consensus;
  consensus.id = 2U;
  consensus.collective_consensus = true;
  passed &= expect_ok(compiler.register_stage(consensus),
                      "a missing consensus stage can be registered after failure");
  passed &= expect_ok(compiler.freeze_for_test(out),
                      "completed registration freezes after prior failure");
  passed &= expect(compiler.frozen(), "successful freeze seals compiler");
  StageSpec late;
  late.id = 3U;
  passed &= expect_invalid(compiler.register_stage(late),
                           "registration after successful freeze is rejected");
  return passed;
}

bool test_registration_replay_is_identical() {
  FrozenExecutionGraph first;
  FrozenExecutionGraph second;
  bool passed = true;
  passed &= expect(build_reference_graph(first),
                   "first deterministic replay compiles");
  passed &= expect(build_reference_graph(second),
                   "second deterministic replay compiles");
  if (passed) {
    passed &= expect(same_graph(first, second),
                     "same registration replay yields identical graph and fingerprint");
  }
  return passed;
}

bool test_fingerprint_covers_compiled_properties() {
  constexpr std::array<GraphFieldSpec, 2U> fields{{
      {accepted(1U), 2U, true},
      {trial(2U), 0U, false},
  }};
  constexpr std::array<FieldAccessSpec, 1U> reads{{accepted(1U)}};
  constexpr std::array<FieldAccessSpec, 1U> writes{{trial(2U)}};
  constexpr std::array<FieldAccessSpec, 1U> ghosts{{accepted(1U)}};
  const auto compile = [&](std::size_t alignment, std::uint8_t ghost_width,
                           std::uint64_t messages,
                           bool writer_consensus,
                           FrozenExecutionGraph& out) {
    ExecutionGraphCompiler compiler;
    if (!compiler.configure(span(fields))) {
      return false;
    }
    const std::array<std::uint8_t, 1U> widths{{ghost_width}};
    StageSpec writer;
    writer.id = 10U;
    writer.reads = span(reads);
    writer.writes = span(writes);
    writer.ghosts = span(ghosts);
    writer.ghost_widths = span(widths);
    writer.workspace_bytes = 128U;
    writer.workspace_alignment = alignment;
    writer.resources.merged_halo_messages = messages;
    writer.collective_consensus = writer_consensus;
    if (!compiler.register_stage(writer)) {
      return false;
    }
    StageSpec consensus;
    consensus.id = 20U;
    consensus.kind = StageKind::service;
    consensus.collective_consensus = true;
    if (!compiler.register_stage(consensus)) {
      return false;
    }
    return static_cast<bool>(compiler.freeze_for_test(out));
  };

  FrozenExecutionGraph baseline;
  FrozenExecutionGraph alignment;
  FrozenExecutionGraph width;
  FrozenExecutionGraph resources;
  FrozenExecutionGraph consensus;
  bool passed = expect(compile(64U, 1U, 3U, true, baseline),
                       "fingerprint baseline compiles");
  passed &= expect(compile(128U, 1U, 3U, true, alignment),
                   "alignment fingerprint mutation compiles");
  passed &= expect(compile(64U, 2U, 3U, true, width),
                   "ghost-width fingerprint mutation compiles");
  passed &= expect(compile(64U, 1U, 4U, true, resources),
                   "resource fingerprint mutation compiles");
  passed &= expect(compile(64U, 1U, 3U, false, consensus),
                   "consensus-placement fingerprint mutation compiles");
  if (passed) {
    passed &= expect(baseline.fingerprint() != alignment.fingerprint() &&
                         baseline.fingerprint() != width.fingerprint() &&
                         baseline.fingerprint() != resources.fingerprint() &&
                         baseline.fingerprint() != consensus.fingerprint(),
                     "every compiled property changes the graph fingerprint");
  }
  return passed;
}

bool test_resource_and_workspace_overflow_is_atomic() {
  constexpr std::array<GraphFieldSpec, 1U> fields{{
      {committed(1U), 0U, true},
  }};
  bool passed = true;
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "resource-overflow fixture configures");
    StageSpec first;
    first.id = 1U;
    first.resources.merged_halo_bytes =
        std::numeric_limits<std::uint64_t>::max();
    StageSpec second;
    second.id = 2U;
    second.resources.merged_halo_bytes = 1U;
    passed &= expect_ok(compiler.register_stage(first),
                        "maximum resource stage registers");
    passed &= expect_ok(compiler.register_stage(second),
                        "overflowing resource stage registers");
    FrozenExecutionGraph out;
    passed &= expect_invalid(compiler.freeze_for_test(out),
                             "resource sum overflow is rejected");
    passed &= expect(graph_is_pristine(out),
                     "resource overflow leaves graph output atomic");
  }
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "workspace-overflow fixture configures");
    StageSpec stage;
    stage.id = 1U;
    stage.workspace_bytes = 64U;
    stage.workspace_alignment = 64U;
    stage.has_fixed_workspace_offset = true;
    stage.fixed_workspace_offset =
        std::numeric_limits<std::size_t>::max() - 63U;
    passed &= expect_ok(compiler.register_stage(stage),
                        "overflowing workspace metadata registers cold");
    FrozenExecutionGraph out;
    passed &= expect_invalid(compiler.freeze_for_test(out),
                             "workspace end overflow is rejected");
    passed &= expect(graph_is_pristine(out),
                     "workspace overflow leaves graph output atomic");
  }
  return passed;
}

bool test_fixed_workspace_overlap_and_reuse() {
  const std::array<GraphFieldSpec, 1U> fields{{
      {committed(1U), 0U, true},
  }};
  bool passed = true;
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "overlapping-workspace fixture configures");
    StageSpec first;
    first.id = 10U;
    first.workspace_bytes = 64U;
    first.workspace_live_through = 20U;
    first.fixed_workspace_offset = 0U;
    first.has_fixed_workspace_offset = true;
    StageSpec second;
    second.id = 20U;
    second.workspace_bytes = 64U;
    second.fixed_workspace_offset = 0U;
    second.has_fixed_workspace_offset = true;
    second.collective_consensus = true;
    passed &= expect_ok(compiler.register_stage(first),
                        "first fixed workspace registers");
    passed &= expect_ok(compiler.register_stage(second),
                        "second fixed workspace registers");
    FrozenExecutionGraph out;
    passed &= expect_invalid(compiler.freeze_for_test(out),
                             "overlapping live fixed workspaces are rejected");
    passed &= expect(graph_is_pristine(out),
                     "workspace-overlap failure leaves output pristine");
  }
  {
    ExecutionGraphCompiler compiler;
    passed &= expect_ok(compiler.configure(span(fields)),
                        "reused-workspace fixture configures");
    StageSpec first;
    first.id = 10U;
    first.workspace_bytes = 64U;
    first.fixed_workspace_offset = 0U;
    first.has_fixed_workspace_offset = true;
    StageSpec second = first;
    second.id = 20U;
    second.collective_consensus = true;
    passed &= expect_ok(compiler.register_stage(first),
                        "first non-overlapping workspace registers");
    passed &= expect_ok(compiler.register_stage(second),
                        "second non-overlapping workspace registers");
    FrozenExecutionGraph out;
    passed &= expect_ok(compiler.freeze_for_test(out),
                        "non-overlapping stages may reuse fixed workspace");
    passed &= expect(out.stages().size == 2U &&
                         out.stages().data[0].workspace_offset == 0U &&
                         out.stages().data[1].workspace_offset == 0U &&
                         out.resources().max_live_workspace_bytes == 64U,
                     "fixed workspace reuse preserves a 64-byte high-water mark");
  }
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  passed &= test_happy_graph_has_exact_order_and_edges();
  passed &= test_reads_and_ghosts_require_production();
  passed &= test_writer_version_rules();
  passed &= test_ghost_contract_validation();
  passed &= test_invalidation_contract_validation();
  passed &= test_service_reads_only_committed_snapshot();
  passed &= test_service_resource_isolation();
  passed &= test_zero_length_nonnull_spans_are_empty();
  passed &= test_commit_requires_prior_collective_consensus();
  passed &= test_commit_requires_consensus_after_last_work();
  passed &= test_commit_is_a_pure_transition();
  passed &= test_registration_spans_are_deep_copied();
  passed &= test_failed_freeze_is_atomic_and_recoverable();
  passed &= test_registration_replay_is_identical();
  passed &= test_fingerprint_covers_compiled_properties();
  passed &= test_resource_and_workspace_overflow_is_atomic();
  passed &= test_fixed_workspace_overlap_and_reuse();
  return passed ? 0 : 1;
}
