// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "field_view_interval_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kMgWorkspacePlan = 7301U;
constexpr std::uint32_t kMgWorkspaceCollective = 7302U;
constexpr std::uint32_t kMgWorkspaceRegistry = 7303U;
constexpr std::uint8_t kAxisX = 1U;
constexpr std::uint8_t kAxisY = 2U;
constexpr std::uint8_t kAxisZ = 4U;
constexpr std::size_t kAlignmentDoubles = 8U;
constexpr std::uint8_t kSlotsPerLevel = 4U;
constexpr std::uint8_t kGhostWidth = 1U;

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(MeshPatch left, MeshPatch right) noexcept {
  return same(left.begin, right.begin) && same(left.cells, right.cells) &&
         same(left.process_grid, right.process_grid) &&
         same(left.process_coord, right.process_coord);
}

bool positive(Int3 value) noexcept {
  return value.x > 0 && value.y > 0 && value.z > 0;
}

bool valid_cycle(MgCycleKind cycle) noexcept {
  return cycle == MgCycleKind::v_cycle || cycle == MgCycleKind::f_cycle;
}

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& out) noexcept {
  if (left > std::numeric_limits<std::size_t>::max() - right) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

bool checked_align(std::size_t value, std::size_t& out) noexcept {
  std::size_t candidate = 0U;
  if (!checked_add(value, kAlignmentDoubles - 1U, candidate)) {
    return false;
  }
  out = candidate & ~(kAlignmentDoubles - 1U);
  return true;
}

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
  return hash;
}

PlanFingerprint finish(std::uint64_t hash) noexcept {
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

void hash_int3(std::uint64_t& hash, Int3 value) noexcept {
  hash = mix(hash, static_cast<std::uint32_t>(value.x));
  hash = mix(hash, static_cast<std::uint32_t>(value.y));
  hash = mix(hash, static_cast<std::uint32_t>(value.z));
}

struct AxisSpacing {
  double minimum{std::numeric_limits<double>::infinity()};
  double maximum{};
  double representative{};
};

AxisSpacing spacing(const AxisMetrics& axis) noexcept {
  const Span<const double> widths = axis.widths();
  if (widths.data == nullptr || widths.size == 0U) {
    return {};
  }
  AxisSpacing result;
  for (std::size_t index = 0U; index < widths.size; ++index) {
    result.minimum = std::min(result.minimum, widths.data[index]);
    result.maximum = std::max(result.maximum, widths.data[index]);
    result.representative += widths.data[index];
  }
  result.representative /= static_cast<double>(widths.size);
  return result;
}

std::uint8_t strategy_mask(const CartesianGeometryPlan& geometry,
                           const MgHierarchyPolicy& policy,
                           std::uint8_t& line_mask) noexcept {
  const AxisSpacing axes[3]{spacing(geometry.x()), spacing(geometry.y()),
                            spacing(geometry.z())};
  const double largest_representative =
      std::max({axes[0].representative, axes[1].representative,
                axes[2].representative});
  line_mask = 0U;
  const std::uint8_t masks[3]{kAxisX, kAxisY, kAxisZ};
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const bool internally_strong =
        axes[axis].minimum > 0.0 &&
        axes[axis].maximum / axes[axis].minimum >=
            policy.anisotropy_threshold;
    const bool cross_axis_strong =
        axes[axis].representative > 0.0 &&
        largest_representative / axes[axis].representative >=
            policy.anisotropy_threshold;
    if (internally_strong || cross_axis_strong) {
      line_mask |= masks[axis];
    }
  }
  const std::uint8_t all = kAxisX | kAxisY | kAxisZ;
  const std::uint8_t candidate = static_cast<std::uint8_t>(all & ~line_mask);
  return candidate == 0U ? all : candidate;
}

CoarseningKind kind(std::uint8_t mask) noexcept {
  switch (mask) {
    case kAxisX | kAxisY | kAxisZ:
      return CoarseningKind::full_xyz;
    case kAxisX | kAxisY:
      return CoarseningKind::semi_xy;
    case kAxisX | kAxisZ:
      return CoarseningKind::semi_xz;
    case kAxisY | kAxisZ:
      return CoarseningKind::semi_yz;
    case kAxisX:
      return CoarseningKind::x_only;
    case kAxisY:
      return CoarseningKind::y_only;
    case kAxisZ:
      return CoarseningKind::z_only;
    default:
      return CoarseningKind::full_xyz;
  }
}

std::int32_t component(Int3 value, std::uint8_t axis) noexcept {
  return axis == kAxisX ? value.x : (axis == kAxisY ? value.y : value.z);
}

void set_component(Int3& value, std::uint8_t axis,
                   std::int32_t selected) noexcept {
  if (axis == kAxisX) {
    value.x = selected;
  } else if (axis == kAxisY) {
    value.y = selected;
  } else {
    value.z = selected;
  }
}

bool coarsen_axis(Int3 global, MeshPatch patch, std::uint8_t axis,
                  std::uint8_t minimum_extent, Int3& next_global,
                  MeshPatch& next_patch) noexcept {
  const std::int32_t extent = component(global, axis);
  const std::int32_t begin = component(patch.begin, axis);
  const std::int32_t end = begin + component(patch.cells, axis);
  if (extent <= static_cast<std::int32_t>(minimum_extent)) {
    return false;
  }
  const std::int32_t coarse_begin = (begin + 1) / 2;
  const std::int32_t coarse_end = (end + 1) / 2;
  if (coarse_end <= coarse_begin) {
    return false;
  }
  set_component(next_global, axis, (extent + 1) / 2);
  set_component(next_patch.begin, axis, coarse_begin);
  set_component(next_patch.cells, axis, coarse_end - coarse_begin);
  return true;
}

bool append_storage(MgWorkspaceLevelRequirements& level,
                    std::size_t& cursor) noexcept {
  const auto ghost = static_cast<std::size_t>(kGhostWidth);
  std::size_t padded_x = 0U;
  std::size_t padded_y = 0U;
  std::size_t padded_z = 0U;
  if (!checked_add(static_cast<std::size_t>(level.patch.cells.x), 2U * ghost,
                   padded_x) ||
      !checked_add(static_cast<std::size_t>(level.patch.cells.y), 2U * ghost,
                   padded_y) ||
      !checked_add(static_cast<std::size_t>(level.patch.cells.z), 2U * ghost,
                   padded_z) ||
      !checked_align(padded_x, level.stride_y) ||
      !checked_multiply(level.stride_y, padded_y, level.stride_z) ||
      !checked_multiply(level.stride_z, padded_z,
                        level.slot_stride_doubles) ||
      !checked_align(cursor, level.offset_doubles)) {
    return false;
  }
  std::size_t level_span = 0U;
  return checked_multiply(level.slot_stride_doubles, kSlotsPerLevel,
                          level_span) &&
         checked_add(level.offset_doubles, level_span, cursor);
}

bool same_level(const MgWorkspaceLevelRequirements& left,
                const MgWorkspaceLevelRequirements& right) noexcept {
  return same(left.global_shape, right.global_shape) &&
         same(left.patch, right.patch) && left.coarsening == right.coarsening &&
         left.coarsen_axis_mask == right.coarsen_axis_mask &&
         left.line_axis_mask == right.line_axis_mask &&
         left.offset_doubles == right.offset_doubles &&
         left.slot_stride_doubles == right.slot_stride_doubles &&
         left.stride_y == right.stride_y && left.stride_z == right.stride_z;
}

bool valid_requirements(const MgWorkspaceRequirements& requirements) noexcept {
  if (requirements.level_count < 2U ||
      requirements.level_count > kMgMaximumLevels ||
      requirements.ghost_width != kGhostWidth ||
      requirements.slots_per_level != kSlotsPerLevel ||
      requirements.total_doubles == 0U ||
      requirements.total_doubles >
          static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) ||
      requirements.arena_shape.x !=
          static_cast<std::int32_t>(requirements.total_doubles) ||
      requirements.arena_shape.y != 1 || requirements.arena_shape.z != 1 ||
      requirements.execution_revision == 0U ||
      requirements.collective_fingerprint == 0U ||
      requirements.fingerprint == 0U) {
    return false;
  }
  std::size_t expected = 0U;
  for (std::size_t index = 0U; index < requirements.level_count; ++index) {
    const MgWorkspaceLevelRequirements& level = requirements.levels[index];
    MgWorkspaceLevelRequirements checked = level;
    if (!positive(level.global_shape) || !positive(level.patch.cells) ||
        !append_storage(checked, expected) || !same_level(level, checked)) {
      return false;
    }
  }
  return expected == requirements.total_doubles;
}

bool same_requirements(const MgWorkspaceRequirements& left,
                       const MgWorkspaceRequirements& right) noexcept {
  if (left.level_count != right.level_count ||
      !same(left.arena_shape, right.arena_shape) ||
      left.ghost_width != right.ghost_width ||
      left.slots_per_level != right.slots_per_level ||
      left.total_doubles != right.total_doubles ||
      left.execution_revision != right.execution_revision ||
      left.collective_fingerprint != right.collective_fingerprint ||
      left.fingerprint != right.fingerprint) {
    return false;
  }
  for (std::size_t index = 0U; index < left.level_count; ++index) {
    if (!same_level(left.levels[index], right.levels[index])) {
      return false;
    }
  }
  return true;
}

bool sufficient_bundle(const MgWorkspaceRequirements& requirements,
                       FieldView bundle) noexcept {
  return bundle.base != nullptr && bundle.components == 1U &&
         bundle.replica == 0U && bundle.revision != 0U &&
         bundle.storage_identity != 0U && bundle.revision_domain != 0U &&
         bundle.ghosts.x >= 0 && bundle.ghosts.y >= 0 &&
         bundle.ghosts.z >= 0 &&
         bundle.interior.x >= requirements.arena_shape.x &&
         bundle.interior.y >= requirements.arena_shape.y &&
         bundle.interior.z >= requirements.arena_shape.z &&
         bundle.stride_y >= static_cast<std::size_t>(bundle.interior.x) &&
         bundle.stride_z >=
             bundle.stride_y * static_cast<std::size_t>(bundle.interior.y) &&
         bundle.component_stride >= requirements.total_doubles;
}

PlanFingerprint binding_fingerprint(const MgWorkspaceRequirements& requirements,
                                    FieldView bundle) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, requirements.fingerprint);
  hash = mix(hash, bundle.field);
  hash = mix(hash, bundle.storage_identity);
  hash = mix(hash, bundle.revision_domain);
  return finish(hash);
}

std::size_t slot_index(std::size_t level, MgWorkspaceSlot slot) noexcept {
  return level * kSlotsPerLevel + static_cast<std::size_t>(slot);
}

}  // namespace

Status make_mg_workspace_requirements(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    MeshPatch finest_patch, const MgHierarchyPolicy& policy,
    RevisionToken execution_revision,
    MgWorkspaceRequirements& out) noexcept {
  const bool local_valid =
      communicator != MPI_COMM_NULL && positive(geometry.global_cells()) &&
      positive(finest_patch.cells) && execution_revision != 0U &&
      std::isfinite(policy.anisotropy_threshold) &&
      policy.anisotropy_threshold > 1.0 && policy.maximum_levels >= 2U &&
      policy.maximum_levels <= kMgMaximumLevels &&
      policy.minimum_coarse_extent >= 2U && valid_cycle(policy.cycle);
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  int local_ok = local_valid ? 1 : 0;
  int global_ok = 0;
  if (MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kMgWorkspaceCollective};
  }
  if (global_ok == 0) {
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  std::uint64_t semantic = UINT64_C(1469598103934665603);
  semantic = mix(semantic, geometry.fingerprint());
  semantic = mix(semantic, execution_revision);
  semantic = mix(semantic, policy.maximum_levels);
  semantic = mix(semantic, policy.minimum_coarse_extent);
  semantic = mix(semantic, static_cast<std::uint64_t>(policy.cycle));
  std::uint64_t threshold_bits = 0U;
  static_assert(sizeof(threshold_bits) == sizeof(policy.anisotropy_threshold));
  std::memcpy(&threshold_bits, &policy.anisotropy_threshold,
              sizeof(threshold_bits));
  semantic = mix(semantic, threshold_bits);
  std::uint64_t minimum_semantic = 0U;
  std::uint64_t maximum_semantic = 0U;
  const int minimum_status = MPI_Allreduce(
      &semantic, &minimum_semantic, 1, MPI_UINT64_T, MPI_MIN, communicator);
  const int maximum_status = MPI_Allreduce(
      &semantic, &maximum_semantic, 1, MPI_UINT64_T, MPI_MAX, communicator);
  if (minimum_status != MPI_SUCCESS || maximum_status != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kMgWorkspaceCollective};
  }
  if (minimum_semantic != maximum_semantic) {
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  MgWorkspaceRequirements candidate;
  candidate.ghost_width = kGhostWidth;
  candidate.slots_per_level = kSlotsPerLevel;
  candidate.execution_revision = execution_revision;
  std::uint8_t line_mask = 0U;
  const std::uint8_t preferred_mask = strategy_mask(geometry, policy,
                                                    line_mask);
  Int3 global = geometry.global_cells();
  MeshPatch patch = finest_patch;
  std::size_t cursor = 0U;
  for (std::size_t level_index = 0U;
       level_index < policy.maximum_levels; ++level_index) {
    MgWorkspaceLevelRequirements level;
    level.global_shape = global;
    level.patch = patch;
    level.line_axis_mask = line_mask;
    Int3 next_global = global;
    MeshPatch next_patch = patch;
    std::uint8_t local_mask = 0U;
    for (const std::uint8_t axis : {kAxisX, kAxisY, kAxisZ}) {
      if ((preferred_mask & axis) != 0U &&
          coarsen_axis(global, patch, axis, policy.minimum_coarse_extent,
                       next_global, next_patch)) {
        local_mask |= axis;
      }
    }
    int local_masks[3]{(local_mask & kAxisX) != 0U ? 1 : 0,
                       (local_mask & kAxisY) != 0U ? 1 : 0,
                       (local_mask & kAxisZ) != 0U ? 1 : 0};
    int global_masks[3]{};
    if (MPI_Allreduce(local_masks, global_masks, 3, MPI_INT, MPI_MIN,
                      communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kMgWorkspaceCollective};
    }
    std::uint8_t actual_mask = 0U;
    actual_mask |= global_masks[0] != 0 ? kAxisX : 0U;
    actual_mask |= global_masks[1] != 0 ? kAxisY : 0U;
    actual_mask |= global_masks[2] != 0 ? kAxisZ : 0U;
    if (actual_mask != local_mask) {
      next_global = global;
      next_patch = patch;
      for (const std::uint8_t axis : {kAxisX, kAxisY, kAxisZ}) {
        if ((actual_mask & axis) != 0U) {
          static_cast<void>(coarsen_axis(
              global, patch, axis, policy.minimum_coarse_extent,
              next_global, next_patch));
        }
      }
    }
    level.coarsening = kind(actual_mask == 0U ? preferred_mask : actual_mask);
    level.coarsen_axis_mask = actual_mask;
    if (!append_storage(level, cursor)) {
      return {StatusCode::invalid_plan, kMgWorkspacePlan};
    }
    candidate.levels[candidate.level_count++] = level;
    if (actual_mask == 0U || candidate.level_count == policy.maximum_levels) {
      break;
    }
    global = next_global;
    patch = next_patch;
  }
  if (candidate.level_count < 2U || cursor == 0U ||
      cursor >
          static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  candidate.total_doubles = cursor;
  candidate.arena_shape = {static_cast<std::int32_t>(cursor), 1, 1};
  std::uint64_t collective = semantic;
  collective = mix(collective, candidate.level_count);
  std::uint64_t local = collective;
  for (std::size_t index = 0U; index < candidate.level_count; ++index) {
    const MgWorkspaceLevelRequirements& level = candidate.levels[index];
    hash_int3(collective, level.global_shape);
    collective = mix(collective, level.coarsen_axis_mask);
    collective = mix(collective, level.line_axis_mask);
    hash_int3(local, level.global_shape);
    hash_int3(local, level.patch.begin);
    hash_int3(local, level.patch.cells);
    local = mix(local, level.coarsen_axis_mask);
    local = mix(local, level.line_axis_mask);
    local = mix(local, level.offset_doubles);
    local = mix(local, level.slot_stride_doubles);
  }
  candidate.collective_fingerprint = finish(collective);
  candidate.fingerprint = finish(local);
  out = candidate;
  return {};
}

Status register_mg_workspace(
    FieldRegistry& registry, std::string_view stable_prefix,
    const MgWorkspaceRequirements& requirements,
    MgWorkspaceFieldIds& out) noexcept {
  if (stable_prefix.empty() || !valid_requirements(requirements)) {
    return {StatusCode::invalid_plan, kMgWorkspaceRegistry};
  }
  try {
    FieldRegistry candidate = registry;
    std::string name(stable_prefix);
    name += ".mg_vectors";
    FieldId field = 0U;
    const Status declared = candidate.declare_field(
        name, 1U, 0U, field);
    if (!declared) {
      return declared;
    }
    registry = std::move(candidate);
    out.vector_bundle = field;
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kMgWorkspaceRegistry};
  }
}

Status MgWorkspace::bind(const MgWorkspaceRequirements& requirements,
                         FieldView vector_bundle,
                         MgWorkspace& out) noexcept {
  if (!valid_requirements(requirements) ||
      !sufficient_bundle(requirements, vector_bundle)) {
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  if (out.fingerprint_ != 0U) {
    if (same_requirements(out.requirements_, requirements) &&
        out.vector_bundle_.base == vector_bundle.base &&
        out.vector_bundle_.field == vector_bundle.field &&
        out.vector_bundle_.storage_identity == vector_bundle.storage_identity &&
        out.vector_bundle_.revision_domain == vector_bundle.revision_domain) {
      return {};
    }
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  MgWorkspace candidate;
  candidate.requirements_ = requirements;
  candidate.vector_bundle_ = vector_bundle;
  candidate.raw_base_ = vector_bundle.base;
  const std::size_t revision_count =
      requirements.level_count * requirements.slots_per_level;
  for (std::size_t index = 0U; index < revision_count; ++index) {
    candidate.revisions_[index] = static_cast<RevisionToken>(index) + 1U;
  }
  candidate.next_revision_ = static_cast<RevisionToken>(revision_count) + 1U;
  candidate.fingerprint_ = binding_fingerprint(requirements, vector_bundle);
  candidate.binding_identity_ = 1U;
  out = candidate;
  return {};
}

const MgWorkspaceLevelRequirements* MgWorkspace::level_requirements(
    std::size_t level_index) const noexcept {
  return fingerprint_ != 0U && level_index < requirements_.level_count
             ? &requirements_.levels[level_index]
             : nullptr;
}

bool MgWorkspace::valid_for(
    const MgWorkspaceRequirements& requirements) const noexcept {
  return fingerprint_ != 0U && same_requirements(requirements_, requirements);
}

FieldView MgWorkspace::level(std::size_t level_index,
                             MgWorkspaceSlot slot) const noexcept {
  const std::size_t selected_slot = static_cast<std::size_t>(slot);
  if (fingerprint_ == 0U || level_index >= requirements_.level_count ||
      selected_slot >= requirements_.slots_per_level) {
    return {};
  }
  const MgWorkspaceLevelRequirements& level = requirements_.levels[level_index];
  const std::size_t ghost = requirements_.ghost_width;
  FieldView result = vector_bundle_;
  result.base = raw_base_ + level.offset_doubles +
                selected_slot * level.slot_stride_doubles + ghost +
                ghost * level.stride_y + ghost * level.stride_z;
  result.interior = level.patch.cells;
  result.ghosts = {requirements_.ghost_width, requirements_.ghost_width,
                   requirements_.ghost_width};
  result.components = 1U;
  result.stride_y = level.stride_y;
  result.stride_z = level.stride_z;
  result.component_stride = level.slot_stride_doubles;
  result.revision = revisions_[slot_index(level_index, slot)];
  return result;
}

Status MgWorkspace::revise_level(std::size_t level_index,
                                 MgWorkspaceSlot slot) noexcept {
  const std::size_t selected = slot_index(level_index, slot);
  if (fingerprint_ == 0U || level_index >= requirements_.level_count ||
      static_cast<std::size_t>(slot) >= requirements_.slots_per_level ||
      selected >= revisions_.size() || revisions_[selected] == 0U ||
      next_revision_ == 0U ||
      next_revision_ == std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kMgWorkspacePlan};
  }
  revisions_[selected] = next_revision_++;
  return {};
}

bool MgWorkspace::overlaps_storage(ConstFieldView view) const noexcept {
  if (fingerprint_ == 0U || view.base == nullptr) {
    return false;
  }
  return detail::field_view_overlaps_storage(
      view, raw_base_, requirements_.total_doubles);
}

bool MgWorkspace::overlaps_storage(FieldView view) const noexcept {
  return overlaps_storage(as_const(view));
}

}  // namespace hundun::v04
