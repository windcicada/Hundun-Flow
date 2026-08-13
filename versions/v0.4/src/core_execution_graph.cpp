// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include "core_arena_detail.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kGraphConfigure = 941U;
constexpr std::uint32_t kGraphStage = 942U;
constexpr std::uint32_t kGraphAccess = 943U;
constexpr std::uint32_t kGraphGhost = 944U;
constexpr std::uint32_t kGraphInvalidation = 945U;
constexpr std::uint32_t kGraphConsensus = 946U;
constexpr std::uint32_t kGraphWorkspace = 947U;
constexpr std::uint32_t kGraphOverflow = 948U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

static_assert(sizeof(double) == sizeof(std::uint64_t));

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

template <class Enum>
bool valid_enum(Enum value, Enum maximum) noexcept {
  using Underlying = std::underlying_type_t<Enum>;
  return static_cast<Underlying>(value) <= static_cast<Underlying>(maximum);
}

bool checked_u64_add(std::uint64_t& destination,
                     std::uint64_t increment) noexcept {
  if (increment > std::numeric_limits<std::uint64_t>::max() - destination) {
    return false;
  }
  destination += increment;
  return true;
}

bool valid_alignment(std::size_t alignment) noexcept {
  return alignment != 0U && (alignment & (alignment - 1U)) == 0U;
}

bool checked_align(std::size_t value, std::size_t alignment,
                   std::size_t& out) noexcept {
  if (!valid_alignment(alignment)) {
    return false;
  }
  return detail::checked_align(value, alignment, out);
}

bool span_shape_valid(const void* data, std::size_t size) noexcept {
  return size == 0U || data != nullptr;
}

bool same_access(FieldAccessSpec left, FieldAccessSpec right) noexcept {
  return left == right;
}

bool contains(Span<const FieldAccessSpec> accesses,
              FieldAccessSpec selected) noexcept {
  for (std::size_t index = 0U; index < accesses.size; ++index) {
    if (same_access(accesses.data[index], selected)) {
      return true;
    }
  }
  return false;
}

bool append_edge(std::vector<GraphEdge>& edges, std::uint32_t from,
                 std::uint32_t to) {
  if (from == to) {
    return true;
  }
  if (!edges.empty() && edges.back().from == from && edges.back().to == to) {
    return true;
  }
  edges.push_back({from, to});
  return true;
}

struct FieldState {
  FieldAccessSpec access{};
  bool available{};
  bool initialized{};
  bool invalidated{};
};

FieldState* find_state(std::vector<FieldState>& states,
                       FieldAccessSpec access) noexcept {
  for (FieldState& state : states) {
    if (same_access(state.access, access)) {
      return &state;
    }
  }
  return nullptr;
}

struct WorkspaceInterval {
  std::size_t begin{};
  std::size_t end{};
  std::size_t offset{};
  std::size_t bytes{};
  std::size_t byte_end{};
};

bool stages_overlap(const WorkspaceInterval& left,
                    const WorkspaceInterval& right) noexcept {
  return left.begin <= right.end && right.begin <= left.end;
}

bool byte_ranges_overlap(const WorkspaceInterval& left,
                         const WorkspaceInterval& right) noexcept {
  if (left.bytes == 0U || right.bytes == 0U) {
    return false;
  }
  return left.offset < right.byte_end && right.offset < left.byte_end;
}

template <class Stage>
Status allocate_workspace(
    const Stage& stage, std::size_t registration,
    std::size_t live_through,
    std::vector<WorkspaceInterval>& active, std::size_t& offset,
    std::size_t& high_water) noexcept {
  offset = 0U;
  if (stage.workspace_bytes == 0U) {
    return {};
  }
  if (live_through < registration ||
      !valid_alignment(stage.workspace_alignment)) {
    return {StatusCode::invalid_plan, kGraphWorkspace};
  }
  WorkspaceInterval candidate{registration, live_through, 0U,
                              stage.workspace_bytes, 0U};
  if (stage.has_fixed_workspace_offset) {
    if (stage.fixed_workspace_offset % stage.workspace_alignment != 0U) {
      return {StatusCode::invalid_plan, kGraphWorkspace};
    }
    candidate.offset = stage.fixed_workspace_offset;
    if (!detail::checked_add(candidate.offset, candidate.bytes,
                             candidate.byte_end)) {
      return {StatusCode::invalid_plan, kGraphOverflow};
    }
    for (const WorkspaceInterval& prior : active) {
      if (stages_overlap(candidate, prior) &&
          byte_ranges_overlap(candidate, prior)) {
        return {StatusCode::invalid_plan, kGraphWorkspace};
      }
    }
  } else {
    std::size_t search = 0U;
    while (true) {
      if (!checked_align(search, stage.workspace_alignment,
                         candidate.offset)) {
        return {StatusCode::invalid_plan, kGraphOverflow};
      }
      if (!detail::checked_add(candidate.offset, candidate.bytes,
                               candidate.byte_end)) {
        return {StatusCode::invalid_plan, kGraphOverflow};
      }
      bool collided = false;
      std::size_t next = candidate.offset;
      for (const WorkspaceInterval& prior : active) {
        if (!stages_overlap(candidate, prior) ||
            !byte_ranges_overlap(candidate, prior)) {
          continue;
        }
        next = prior.byte_end;
        collided = true;
        break;
      }
      if (!collided) {
        break;
      }
      search = next;
    }
  }
  high_water = std::max(high_water, candidate.byte_end);
  offset = candidate.offset;
  active.push_back(candidate);
  return {};
}

void hash_access(std::uint64_t& hash, FieldAccessSpec access) noexcept {
  hash = hash_mix(hash, access.field);
  hash = hash_mix(hash, static_cast<std::uint8_t>(access.visibility));
}

}  // namespace

const FrozenStage* FrozenExecutionGraph::stage(StageId id) const noexcept {
  for (const FrozenStage& selected : stages_) {
    if (selected.id == id) {
      return &selected;
    }
  }
  return nullptr;
}

Span<const FieldAccessSpec> FrozenExecutionGraph::reads(
    StageId stage_id) const noexcept {
  const FrozenStage* const selected = stage(stage_id);
  if (selected == nullptr) {
    return {};
  }
  const std::size_t begin = selected->read_begin;
  const std::size_t count = selected->read_count;
  return begin > reads_.size() || count > reads_.size() - begin || count == 0U
             ? Span<const FieldAccessSpec>{}
             : Span<const FieldAccessSpec>{reads_.data() + begin, count};
}

Span<const FieldAccessSpec> FrozenExecutionGraph::writes(
    StageId stage_id) const noexcept {
  const FrozenStage* const selected = stage(stage_id);
  if (selected == nullptr) {
    return {};
  }
  const std::size_t begin = selected->write_begin;
  const std::size_t count = selected->write_count;
  return begin > writes_.size() || count > writes_.size() - begin ||
                 count == 0U
             ? Span<const FieldAccessSpec>{}
             : Span<const FieldAccessSpec>{writes_.data() + begin, count};
}

Span<const FieldAccessSpec> FrozenExecutionGraph::ghosts(
    StageId stage_id) const noexcept {
  const FrozenStage* const selected = stage(stage_id);
  if (selected == nullptr) {
    return {};
  }
  const std::size_t begin = selected->ghost_begin;
  const std::size_t count = selected->ghost_count;
  return begin > ghosts_.size() || count > ghosts_.size() - begin ||
                 count == 0U
             ? Span<const FieldAccessSpec>{}
             : Span<const FieldAccessSpec>{ghosts_.data() + begin, count};
}

Span<const std::uint8_t> FrozenExecutionGraph::ghost_widths(
    StageId stage_id) const noexcept {
  const FrozenStage* const selected = stage(stage_id);
  if (selected == nullptr) {
    return {};
  }
  const std::size_t begin = selected->ghost_begin;
  const std::size_t count = selected->ghost_count;
  return begin > ghost_widths_.size() ||
                 count > ghost_widths_.size() - begin || count == 0U
             ? Span<const std::uint8_t>{}
             : Span<const std::uint8_t>{ghost_widths_.data() + begin, count};
}

Span<const FieldAccessSpec> FrozenExecutionGraph::invalidations(
    StageId stage_id) const noexcept {
  const FrozenStage* const selected = stage(stage_id);
  if (selected == nullptr) {
    return {};
  }
  const std::size_t begin = selected->invalidation_begin;
  const std::size_t count = selected->invalidation_count;
  return begin > invalidations_.size() ||
                 count > invalidations_.size() - begin || count == 0U
             ? Span<const FieldAccessSpec>{}
             : Span<const FieldAccessSpec>{invalidations_.data() + begin,
                                           count};
}

Status ExecutionGraphCompiler::configure(
    Span<const GraphFieldSpec> fields) noexcept {
  if (configured_ || frozen_ || !span_shape_valid(fields.data, fields.size) ||
      fields.size == 0U) {
    return {StatusCode::invalid_plan, kGraphConfigure};
  }
  try {
    std::vector<GraphFieldSpec> candidate(fields.data,
                                          fields.data + fields.size);
    for (std::size_t index = 0U; index < candidate.size(); ++index) {
      const GraphFieldSpec& field = candidate[index];
      if (!valid_enum(field.access.visibility, StateVisibility::workspace)) {
        return {StatusCode::invalid_plan, kGraphConfigure};
      }
      for (std::size_t prior = 0U; prior < index; ++prior) {
        if (same_access(candidate[prior].access, field.access)) {
          return {StatusCode::invalid_plan, kGraphConfigure};
        }
      }
    }
    fields_.swap(candidate);
    configured_ = true;
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kGraphConfigure};
  }
}

Status ExecutionGraphCompiler::register_stage(
    const StageSpec& stage) noexcept {
  if (!configured_ || frozen_ || stage.id == 0U ||
      !valid_enum(stage.kind, StageKind::commit) ||
      !span_shape_valid(stage.reads.data, stage.reads.size) ||
      !span_shape_valid(stage.writes.data, stage.writes.size) ||
      !span_shape_valid(stage.ghosts.data, stage.ghosts.size) ||
      !span_shape_valid(stage.ghost_widths.data,
                        stage.ghost_widths.size) ||
      !span_shape_valid(stage.invalidates.data, stage.invalidates.size) ||
      stage.ghosts.size != stage.ghost_widths.size ||
      (stage.workspace_bytes != 0U &&
       !valid_alignment(stage.workspace_alignment))) {
    return {StatusCode::invalid_plan, kGraphStage};
  }
  for (const OwnedStageSpec& prior : stages_) {
    if (prior.id == stage.id) {
      return {StatusCode::invalid_plan, kGraphStage};
    }
  }
  const auto valid_access_span = [this](Span<const FieldAccessSpec> accesses) {
    for (std::size_t index = 0U; index < accesses.size; ++index) {
      if (!valid_enum(accesses.data[index].visibility,
                      StateVisibility::workspace) ||
          std::none_of(fields_.begin(), fields_.end(),
                       [&](const GraphFieldSpec& field) {
                         return same_access(field.access,
                                            accesses.data[index]);
                       })) {
        return false;
      }
      for (std::size_t prior = 0U; prior < index; ++prior) {
        if (same_access(accesses.data[prior], accesses.data[index])) {
          return false;
        }
      }
    }
    return true;
  };
  if (!valid_access_span(stage.reads) || !valid_access_span(stage.writes) ||
      !valid_access_span(stage.ghosts) ||
      !valid_access_span(stage.invalidates)) {
    return {StatusCode::invalid_plan, kGraphAccess};
  }
  for (std::size_t index = 0U; index < stage.ghosts.size; ++index) {
    const FieldAccessSpec ghost = stage.ghosts.data[index];
    const auto declared = std::find_if(
        fields_.begin(), fields_.end(), [&](const GraphFieldSpec& field) {
          return same_access(field.access, ghost);
        });
    if (!contains(stage.reads, ghost) || stage.ghost_widths.data[index] == 0U ||
        declared == fields_.end() ||
        stage.ghost_widths.data[index] > declared->ghost_capacity) {
      return {StatusCode::invalid_plan, kGraphGhost};
    }
  }
  if (stage.kind == StageKind::service) {
    for (std::size_t index = 0U; index < stage.reads.size; ++index) {
      if (stage.reads.data[index].visibility !=
          StateVisibility::committed_snapshot) {
        return {StatusCode::invalid_plan, kGraphAccess};
      }
    }
    if (stage.writes.size != 0U || stage.invalidates.size != 0U ||
        stage.ghosts.size != 0U ||
        stage.resources.merged_halo_messages != 0U ||
        stage.resources.merged_halo_bytes != 0U ||
        stage.resources.numeric_refills != 0U ||
        stage.resources.hierarchy_rebuilds != 0U ||
        stage.resources.cache_publishes != 0U ||
        stage.resources.linear_iterations != 0U) {
      return {StatusCode::invalid_plan, kGraphAccess};
    }
  }
  if (stage.kind == StageKind::commit &&
      (stage.reads.size != 0U || stage.writes.size != 0U ||
       stage.invalidates.size != 0U || stage.ghosts.size != 0U ||
       stage.resources.merged_halo_messages != 0U ||
       stage.resources.merged_halo_bytes != 0U ||
       stage.resources.numeric_refills != 0U ||
       stage.resources.hierarchy_rebuilds != 0U ||
       stage.resources.cache_publishes != 0U ||
       stage.resources.linear_iterations != 0U ||
       stage.resources.stage_wall_nanoseconds != 0U ||
       stage.workspace_bytes != 0U ||
       stage.workspace_alignment != 64U ||
       stage.workspace_live_through != 0U ||
       stage.fixed_workspace_offset != 0U ||
       stage.has_fixed_workspace_offset ||
       stage.collective_consensus)) {
    return {StatusCode::invalid_plan, kGraphStage};
  }
  try {
    OwnedStageSpec candidate;
    candidate.id = stage.id;
    if (stage.reads.size != 0U) {
      candidate.reads.assign(stage.reads.data,
                             stage.reads.data + stage.reads.size);
    }
    if (stage.writes.size != 0U) {
      candidate.writes.assign(stage.writes.data,
                              stage.writes.data + stage.writes.size);
    }
    if (stage.ghosts.size != 0U) {
      candidate.ghosts.assign(stage.ghosts.data,
                              stage.ghosts.data + stage.ghosts.size);
      candidate.ghost_widths.assign(stage.ghost_widths.data,
                                    stage.ghost_widths.data +
                                        stage.ghost_widths.size);
    }
    if (stage.invalidates.size != 0U) {
      candidate.invalidates.assign(stage.invalidates.data,
                                   stage.invalidates.data +
                                       stage.invalidates.size);
    }
    candidate.workspace_bytes = stage.workspace_bytes;
    candidate.workspace_alignment = stage.workspace_alignment;
    candidate.workspace_live_through = stage.workspace_live_through;
    candidate.fixed_workspace_offset = stage.fixed_workspace_offset;
    candidate.resources = stage.resources;
    candidate.kind = stage.kind;
    candidate.has_fixed_workspace_offset = stage.has_fixed_workspace_offset;
    candidate.collective_consensus = stage.collective_consensus;
    stages_.push_back(std::move(candidate));
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kGraphStage};
  }
}

Status ExecutionGraphCompiler::freeze_for_test(
    FrozenExecutionGraph& out) noexcept {
  if (!configured_ || frozen_ || stages_.empty() ||
      out.fingerprint_ != 0U || !out.stages_.empty() ||
      !out.nodes_.empty() || !out.edges_.empty()) {
    return {StatusCode::invalid_plan, kGraphConfigure};
  }
  try {
    FrozenExecutionGraph candidate;
    std::vector<FieldState> states;
    states.reserve(fields_.size());
    for (const GraphFieldSpec& field : fields_) {
      states.push_back({field.access, field.initially_available,
                        field.initially_available, false});
    }
    std::vector<WorkspaceInterval> workspaces;
    std::size_t workspace_high_water = 0U;
    std::uint32_t last_node = std::numeric_limits<std::uint32_t>::max();
    bool unconsensed_work = false;
    bool has_consensus = false;
    bool seen_commit = false;
    ResourceContract resources{};

    for (std::size_t ordinal = 0U; ordinal < stages_.size(); ++ordinal) {
      const OwnedStageSpec& registered = stages_[ordinal];
      if (ordinal > std::numeric_limits<std::uint32_t>::max() ||
          (seen_commit && registered.kind != StageKind::service) ||
          (registered.kind == StageKind::commit &&
           (seen_commit || unconsensed_work || !has_consensus))) {
        return {StatusCode::invalid_plan, kGraphConsensus};
      }
      const auto as_span = [](const std::vector<FieldAccessSpec>& values) {
        return Span<const FieldAccessSpec>{values.data(), values.size()};
      };
      const Span<const FieldAccessSpec> reads = as_span(registered.reads);
      const Span<const FieldAccessSpec> writes = as_span(registered.writes);
      const Span<const FieldAccessSpec> invalidates =
          as_span(registered.invalidates);
      for (const FieldAccessSpec read : registered.reads) {
        const FieldState* state = find_state(states, read);
        if (state == nullptr || !state->available || state->invalidated) {
          return {StatusCode::invalid_plan, kGraphAccess};
        }
      }
      for (const FieldAccessSpec invalidated : registered.invalidates) {
        FieldState* state = find_state(states, invalidated);
        if (state == nullptr || !state->available || state->invalidated ||
            !contains(reads, invalidated) || !contains(writes, invalidated)) {
          return {StatusCode::invalid_plan, kGraphInvalidation};
        }
      }
      for (const FieldAccessSpec write : registered.writes) {
        FieldState* state = find_state(states, write);
        if (state == nullptr) {
          return {StatusCode::invalid_plan, kGraphAccess};
        }
        if (state->initialized &&
            (!contains(reads, write) || !contains(invalidates, write))) {
          return {StatusCode::invalid_plan, kGraphAccess};
        }
      }
      std::size_t live_through_ordinal = ordinal;
      if (registered.workspace_live_through != 0U) {
        const auto lifetime_end = std::find_if(
            stages_.begin() + static_cast<std::ptrdiff_t>(ordinal),
            stages_.end(), [&](const OwnedStageSpec& selected) {
              return selected.id == registered.workspace_live_through;
            });
        if (lifetime_end == stages_.end()) {
          return {StatusCode::invalid_plan, kGraphWorkspace};
        }
        live_through_ordinal = static_cast<std::size_t>(
            std::distance(stages_.begin(), lifetime_end));
      }
      std::size_t workspace_offset = 0U;
      const Status allocated = allocate_workspace(
          registered, ordinal, live_through_ordinal, workspaces,
          workspace_offset,
          workspace_high_water);
      if (!allocated) {
        return allocated;
      }
      const std::uint32_t node_begin =
          static_cast<std::uint32_t>(candidate.nodes_.size());
      const auto add_node = [&](GraphNodeKind kind) {
        if (candidate.nodes_.size() >=
            static_cast<std::size_t>(
                std::numeric_limits<std::uint32_t>::max())) {
          return false;
        }
        const std::uint32_t index =
            static_cast<std::uint32_t>(candidate.nodes_.size());
        candidate.nodes_.push_back(
            {registered.id, kind, static_cast<std::uint32_t>(ordinal)});
        // Registration order is the numerical order. Keeping this explicit
        // chain prevents algebraic stages from being reordered while still
        // allowing nonblocking halo begin/interior/finish overlap.
        if (last_node != std::numeric_limits<std::uint32_t>::max()) {
          append_edge(candidate.edges_, last_node, index);
        }
        last_node = index;
        return true;
      };
      if (!registered.ghosts.empty()) {
        if (!add_node(GraphNodeKind::halo_begin) ||
            !add_node(GraphNodeKind::compute_interior) ||
            !add_node(GraphNodeKind::halo_finish) ||
            !add_node(GraphNodeKind::compute_boundary)) {
          return {StatusCode::invalid_plan, kGraphOverflow};
        }
      } else {
        const GraphNodeKind kind =
            registered.kind == StageKind::service
                ? GraphNodeKind::service
                : (registered.kind == StageKind::commit
                       ? GraphNodeKind::commit
                       : GraphNodeKind::compute);
        if (!add_node(kind)) {
          return {StatusCode::invalid_plan, kGraphOverflow};
        }
      }
      for (const FieldAccessSpec invalidated : registered.invalidates) {
        FieldState* state = find_state(states, invalidated);
        state->invalidated = true;
      }
      for (const FieldAccessSpec write : registered.writes) {
        FieldState* state = find_state(states, write);
        state->available = true;
        state->initialized = true;
        state->invalidated = false;
      }
      const bool may_fail_or_mutate =
          registered.kind != StageKind::commit ||
          !registered.writes.empty() ||
          !registered.invalidates.empty() ||
          registered.resources.numeric_refills != 0U ||
          registered.resources.hierarchy_rebuilds != 0U ||
          registered.resources.cache_publishes != 0U;
      if (may_fail_or_mutate) {
        unconsensed_work = true;
      }
      if (registered.collective_consensus) {
        if (!add_node(GraphNodeKind::collective_consensus)) {
          return {StatusCode::invalid_plan, kGraphOverflow};
        }
        unconsensed_work = false;
        has_consensus = true;
      }
      seen_commit = seen_commit || registered.kind == StageKind::commit;

      const auto fit_u16 = [](std::size_t value) {
        return value <= std::numeric_limits<std::uint16_t>::max();
      };
      const auto fits_appended_u32 = [](std::size_t current,
                                        std::size_t appended) {
        constexpr std::size_t maximum =
            static_cast<std::size_t>(
                std::numeric_limits<std::uint32_t>::max());
        return current <= maximum && appended <= maximum - current;
      };
      if (!fits_appended_u32(candidate.reads_.size(),
                             registered.reads.size()) ||
          !fits_appended_u32(candidate.writes_.size(),
                             registered.writes.size()) ||
          !fits_appended_u32(candidate.ghosts_.size(),
                             registered.ghosts.size()) ||
          !fits_appended_u32(candidate.invalidations_.size(),
                             registered.invalidates.size()) ||
          !fit_u16(registered.reads.size()) ||
          !fit_u16(registered.writes.size()) ||
          !fit_u16(registered.ghosts.size()) ||
          !fit_u16(registered.invalidates.size()) ||
          candidate.nodes_.size() - node_begin >
              std::numeric_limits<std::uint16_t>::max()) {
        return {StatusCode::invalid_plan, kGraphOverflow};
      }
      FrozenStage frozen;
      frozen.id = registered.id;
      frozen.kind = registered.kind;
      frozen.registration_ordinal = static_cast<std::uint32_t>(ordinal);
      frozen.node_begin = node_begin;
      frozen.node_count = static_cast<std::uint16_t>(
          candidate.nodes_.size() - node_begin);
      frozen.read_begin = static_cast<std::uint32_t>(candidate.reads_.size());
      frozen.read_count = static_cast<std::uint16_t>(registered.reads.size());
      frozen.write_begin =
          static_cast<std::uint32_t>(candidate.writes_.size());
      frozen.write_count =
          static_cast<std::uint16_t>(registered.writes.size());
      frozen.ghost_begin =
          static_cast<std::uint32_t>(candidate.ghosts_.size());
      frozen.ghost_count =
          static_cast<std::uint16_t>(registered.ghosts.size());
      frozen.invalidation_begin =
          static_cast<std::uint32_t>(candidate.invalidations_.size());
      frozen.invalidation_count =
          static_cast<std::uint16_t>(registered.invalidates.size());
      frozen.workspace_offset = workspace_offset;
      frozen.workspace_bytes = registered.workspace_bytes;
      frozen.workspace_alignment = registered.workspace_alignment;
      frozen.workspace_live_through = registered.workspace_live_through;
      frozen.resources = registered.resources;
      frozen.collective_consensus = registered.collective_consensus;
      candidate.reads_.insert(candidate.reads_.end(), registered.reads.begin(),
                              registered.reads.end());
      candidate.writes_.insert(candidate.writes_.end(),
                               registered.writes.begin(),
                               registered.writes.end());
      candidate.ghosts_.insert(candidate.ghosts_.end(),
                               registered.ghosts.begin(),
                               registered.ghosts.end());
      candidate.ghost_widths_.insert(candidate.ghost_widths_.end(),
                                     registered.ghost_widths.begin(),
                                     registered.ghost_widths.end());
      candidate.invalidations_.insert(candidate.invalidations_.end(),
                                      registered.invalidates.begin(),
                                      registered.invalidates.end());
      candidate.stages_.push_back(frozen);
      if (!checked_u64_add(resources.merged_halo_messages,
                           registered.resources.merged_halo_messages) ||
          !checked_u64_add(resources.merged_halo_bytes,
                           registered.resources.merged_halo_bytes) ||
          !checked_u64_add(resources.numeric_refills,
                           registered.resources.numeric_refills) ||
          !checked_u64_add(resources.hierarchy_rebuilds,
                           registered.resources.hierarchy_rebuilds) ||
          !checked_u64_add(resources.cache_publishes,
                           registered.resources.cache_publishes) ||
          !checked_u64_add(resources.linear_iterations,
                           registered.resources.linear_iterations) ||
          !checked_u64_add(resources.stage_wall_nanoseconds,
                           registered.resources.stage_wall_nanoseconds)) {
        return {StatusCode::invalid_plan, kGraphOverflow};
      }
    }
    if (unconsensed_work) {
      return {StatusCode::invalid_plan, kGraphConsensus};
    }
    resources.max_live_workspace_bytes = workspace_high_water;
    resources.allocation_allowance = 0U;
    candidate.resources_ = resources;

    std::uint64_t hash = kFnvOffset;
    hash = hash_mix(hash, UINT64_C(0x6669656c6473));
    hash = hash_mix(hash, fields_.size());
    for (const GraphFieldSpec& field : fields_) {
      hash_access(hash, field.access);
      hash = hash_mix(hash, field.ghost_capacity);
      hash = hash_mix(hash, field.initially_available ? 1U : 0U);
    }
    hash = hash_mix(hash, UINT64_C(0x737461676573));
    hash = hash_mix(hash, candidate.stages_.size());
    for (const FrozenStage& stage : candidate.stages_) {
      hash = hash_mix(hash, stage.id);
      hash = hash_mix(hash, static_cast<std::uint8_t>(stage.kind));
      hash = hash_mix(hash, stage.registration_ordinal);
      hash = hash_mix(hash, stage.node_begin);
      hash = hash_mix(hash, stage.node_count);
      hash = hash_mix(hash, stage.workspace_offset);
      hash = hash_mix(hash, stage.workspace_bytes);
      hash = hash_mix(hash, stage.workspace_alignment);
      hash = hash_mix(hash, stage.workspace_live_through);
      hash = hash_mix(hash, stage.collective_consensus ? 1U : 0U);
      hash = hash_mix(hash, stage.resources.merged_halo_messages);
      hash = hash_mix(hash, stage.resources.merged_halo_bytes);
      hash = hash_mix(hash, stage.resources.numeric_refills);
      hash = hash_mix(hash, stage.resources.hierarchy_rebuilds);
      hash = hash_mix(hash, stage.resources.cache_publishes);
      hash = hash_mix(hash, stage.resources.linear_iterations);
      hash = hash_mix(hash, stage.resources.stage_wall_nanoseconds);
      hash = hash_mix(hash, UINT64_C(0x7265616473));
      hash = hash_mix(hash, stage.read_count);
      const Span<const FieldAccessSpec> stage_reads = candidate.reads(stage.id);
      for (std::size_t index = 0U; index < stage_reads.size; ++index) {
        hash_access(hash, stage_reads.data[index]);
      }
      const Span<const FieldAccessSpec> stage_writes =
          candidate.writes(stage.id);
      hash = hash_mix(hash, UINT64_C(0x777269746573));
      hash = hash_mix(hash, stage.write_count);
      for (std::size_t index = 0U; index < stage_writes.size; ++index) {
        hash_access(hash, stage_writes.data[index]);
      }
      const Span<const FieldAccessSpec> ghosts = candidate.ghosts(stage.id);
      const Span<const std::uint8_t> widths =
          candidate.ghost_widths(stage.id);
      hash = hash_mix(hash, UINT64_C(0x67686f737473));
      hash = hash_mix(hash, stage.ghost_count);
      for (std::size_t ghost = 0U; ghost < ghosts.size; ++ghost) {
        hash_access(hash, ghosts.data[ghost]);
        hash = hash_mix(hash, widths.data[ghost]);
      }
      const Span<const FieldAccessSpec> stage_invalidations =
          candidate.invalidations(stage.id);
      hash = hash_mix(hash, UINT64_C(0x696e76616c6964));
      hash = hash_mix(hash, stage.invalidation_count);
      for (std::size_t index = 0U; index < stage_invalidations.size;
           ++index) {
        hash_access(hash, stage_invalidations.data[index]);
      }
    }
    hash = hash_mix(hash, UINT64_C(0x6e6f646573));
    hash = hash_mix(hash, candidate.nodes_.size());
    for (const GraphNode& node : candidate.nodes_) {
      hash = hash_mix(hash, node.stage);
      hash = hash_mix(hash, static_cast<std::uint8_t>(node.kind));
      hash = hash_mix(hash, node.ordinal);
    }
    hash = hash_mix(hash, UINT64_C(0x6564676573));
    hash = hash_mix(hash, candidate.edges_.size());
    for (const GraphEdge& edge : candidate.edges_) {
      hash = hash_mix(hash, edge.from);
      hash = hash_mix(hash, edge.to);
    }
    hash = hash_mix(hash, UINT64_C(0x7265736f75726365));
    hash = hash_mix(hash, resources.max_live_workspace_bytes);
    hash = hash_mix(hash, resources.allocation_allowance);
    hash = hash_mix(hash, resources.merged_halo_messages);
    hash = hash_mix(hash, resources.merged_halo_bytes);
    hash = hash_mix(hash, resources.numeric_refills);
    hash = hash_mix(hash, resources.hierarchy_rebuilds);
    hash = hash_mix(hash, resources.cache_publishes);
    hash = hash_mix(hash, resources.linear_iterations);
    hash = hash_mix(hash, resources.stage_wall_nanoseconds);
    candidate.fingerprint_ = hash == 0U ? 1U : hash;
    out = std::move(candidate);
    frozen_ = true;
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kGraphStage};
  }
}

}  // namespace hundun::v04
