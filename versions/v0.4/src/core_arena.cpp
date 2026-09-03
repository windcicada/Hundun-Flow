// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_execution.hpp"

#include "core_arena_detail.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kArenaSchema = 1U;
constexpr std::uint32_t kArenaRequest = 2U;
constexpr std::uint32_t kArenaExtent = 3U;
constexpr std::uint32_t kArenaOverflow = 4U;
constexpr std::uint32_t kArenaOwner = 5U;
constexpr std::uint32_t kArenaView = 6U;
constexpr std::uint32_t kArenaRevision = 7U;

bool checked_ghost_extent(std::int32_t interior, std::int32_t ghost,
                          std::size_t& out) noexcept {
  if (interior <= 0 || ghost < 0) {
    return false;
  }
  std::size_t two_ghosts = 0U;
  return detail::checked_multiply(static_cast<std::size_t>(ghost), 2U,
                                  two_ghosts) &&
         detail::checked_add(static_cast<std::size_t>(interior), two_ghosts,
                             out);
}

std::size_t role_slot(StateRole role) noexcept {
  switch (role) {
    case StateRole::accepted_n:
      return 0U;
    case StateRole::accepted_n_minus_one:
      return 1U;
    case StateRole::trial:
      return 2U;
  }
  return 3U;
}

template <class T>
Status make_view(T* owner_base, const ArenaFieldLayout& layout,
                 std::size_t replica, RevisionToken revision,
                 StorageIdentity storage_identity,
                 RevisionDomainIdentity revision_domain,
                 BasicFieldView<T>& out) noexcept {
  if (owner_base == nullptr || replica >= layout.replicas || revision == 0U ||
      storage_identity == 0U || revision_domain == 0U) {
    return {StatusCode::invalid_plan, kArenaView};
  }
  std::size_t replica_offset = 0U;
  if (!detail::checked_multiply(replica, layout.replica_stride_doubles,
                                replica_offset)) {
    return {StatusCode::invalid_plan, kArenaOverflow};
  }
  BasicFieldView<T> candidate;
  candidate.base = owner_base + layout.offset_doubles + replica_offset;
  candidate.interior = layout.interior;
  candidate.ghosts = layout.ghosts;
  candidate.components = layout.components;
  candidate.stride_y = layout.stride_y;
  candidate.stride_z = layout.stride_z;
  candidate.component_stride = layout.component_stride;
  candidate.replica = replica;
  candidate.field = layout.id;
  candidate.revision = revision;
  candidate.storage_identity = storage_identity;
  candidate.revision_domain = revision_domain;
  out = candidate;
  return {};
}

}  // namespace

namespace detail {

std::uint64_t issue_identity() noexcept {
  static std::atomic<std::uint64_t> next{1U};
  std::uint64_t candidate = next.load(std::memory_order_relaxed);
  while (candidate != std::numeric_limits<std::uint64_t>::max()) {
    if (next.compare_exchange_weak(candidate, candidate + 1U,
                                   std::memory_order_relaxed,
                                   std::memory_order_relaxed)) {
      return candidate;
    }
  }
  return 0U;
}

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& out) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
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

bool checked_align(std::size_t value, std::size_t alignment,
                   std::size_t& out) noexcept {
  if (alignment == 0U || (alignment & (alignment - 1U)) != 0U) {
    return false;
  }
  std::size_t candidate = 0U;
  if (!checked_add(value, alignment - 1U, candidate)) {
    return false;
  }
  out = candidate & ~(alignment - 1U);
  return true;
}

}  // namespace detail

Status ArenaLayout::compile(const FieldSchema& schema,
                            Span<const ArenaFieldRequest> requests,
                            ArenaLayout& out) {
  if (schema.empty() || requests.data == nullptr ||
      requests.size != schema.size()) {
    return {StatusCode::invalid_plan, kArenaSchema};
  }
  try {
    ArenaLayout candidate;
    candidate.fields_.resize(schema.size());
    std::vector<std::uint8_t> seen(schema.size(), 0U);
    std::vector<ArenaFieldRequest> by_field(requests.data,
                                            requests.data + requests.size);
    std::sort(by_field.begin(), by_field.end(),
              [](const ArenaFieldRequest& left,
                 const ArenaFieldRequest& right) {
                return left.id < right.id;
              });
    for (const ArenaFieldRequest& request : by_field) {
      const std::size_t id = static_cast<std::size_t>(request.id);
      if (id >= schema.size() || schema[id].id != request.id ||
          seen[id] != 0U ||
          request.lifetime > FieldLifetime::step_scratch) {
        return {StatusCode::invalid_plan, kArenaRequest};
      }
      seen[id] = 1U;
    }
    std::vector<ArenaFieldRequest> by_owner = by_field;
    std::sort(by_owner.begin(), by_owner.end(),
              [](const ArenaFieldRequest& left,
                 const ArenaFieldRequest& right) {
                return left.placement.owner < right.placement.owner ||
                       (left.placement.owner == right.placement.owner &&
                        left.id < right.id);
              });
    candidate.owners_.reserve(by_owner.size());

    std::size_t request_index = 0U;
    while (request_index < by_owner.size()) {
      const ArenaOwnerId owner_id = by_owner[request_index].placement.owner;
      const std::size_t owner_index = candidate.owners_.size();
      std::size_t owner_cursor = 0U;
      while (request_index < by_owner.size() &&
             by_owner[request_index].placement.owner == owner_id) {
        const ArenaFieldRequest& request = by_owner[request_index];
        const FieldDescriptor& descriptor =
            schema[static_cast<std::size_t>(request.id)];
        std::size_t x_ghosted = 0U;
        std::size_t y_ghosted = 0U;
        std::size_t z_ghosted = 0U;
        const auto ghost = static_cast<std::int32_t>(descriptor.ghost_width);
        if (!checked_ghost_extent(request.interior.x, ghost, x_ghosted) ||
            !checked_ghost_extent(request.interior.y, ghost, y_ghosted) ||
            !checked_ghost_extent(request.interior.z, ghost, z_ghosted)) {
          return {StatusCode::invalid_plan, kArenaExtent};
        }

        std::size_t leading_x = 0U;
        if (!detail::checked_align(static_cast<std::size_t>(ghost),
                                   detail::kDoublesPerCacheLine, leading_x)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        std::size_t x_without_leading = 0U;
        if (!detail::checked_add(static_cast<std::size_t>(request.interior.x),
                                 static_cast<std::size_t>(ghost),
                                 x_without_leading)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        std::size_t row_width = 0U;
        if (!detail::checked_add(leading_x, x_without_leading, row_width)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        std::size_t stride_y = 0U;
        if (!detail::checked_align(row_width, detail::kDoublesPerCacheLine,
                                   stride_y)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        std::size_t stride_z = 0U;
        std::size_t component_stride = 0U;
        std::size_t raw_component_span = 0U;
        std::size_t replica_stride = 0U;
        if (!detail::checked_multiply(stride_y, y_ghosted, stride_z) ||
            !detail::checked_multiply(stride_z, z_ghosted,
                                      component_stride) ||
            !detail::checked_multiply(
                component_stride,
                static_cast<std::size_t>(descriptor.components),
                raw_component_span) ||
            !detail::checked_align(raw_component_span,
                                   detail::kDoublesPerCacheLine,
                                   replica_stride)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        const std::size_t replicas =
            request.lifetime == FieldLifetime::state_layer ? 3U : 1U;
        std::size_t span = 0U;
        if (!detail::checked_multiply(replica_stride, replicas, span)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }

        std::size_t raw_start = 0U;
        std::size_t ghost_y_offset = 0U;
        std::size_t ghost_z_offset = 0U;
        if (!detail::checked_align(owner_cursor,
                                   detail::kDoublesPerCacheLine, raw_start) ||
            !detail::checked_multiply(static_cast<std::size_t>(ghost), stride_y,
                                      ghost_y_offset) ||
            !detail::checked_multiply(static_cast<std::size_t>(ghost), stride_z,
                                      ghost_z_offset)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        std::size_t interior_offset = raw_start;
        if (!detail::checked_add(interior_offset, leading_x,
                                 interior_offset) ||
            !detail::checked_add(interior_offset, ghost_y_offset,
                                 interior_offset) ||
            !detail::checked_add(interior_offset, ghost_z_offset,
                                 interior_offset)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        std::size_t field_end = 0U;
        if (!detail::checked_add(raw_start, span, field_end) ||
            field_end > static_cast<std::size_t>(
                            std::numeric_limits<std::ptrdiff_t>::max()) /
                            sizeof(double)) {
          return {StatusCode::invalid_plan, kArenaOverflow};
        }
        owner_cursor = field_end;
        candidate.fields_[static_cast<std::size_t>(request.id)] =
            ArenaFieldLayout{request.id,
                             request.interior,
                             Int3{ghost, ghost, ghost},
                             descriptor.components,
                             stride_y,
                             stride_z,
                             component_stride,
                             replicas,
                             replica_stride,
                             owner_id,
                             owner_index,
                             raw_start,
                             interior_offset,
                             span,
                             request.lifetime};
        static_cast<void>(x_ghosted);
        ++request_index;
      }
      std::size_t owner_span = 0U;
      if (!detail::checked_align(owner_cursor, detail::kDoublesPerCacheLine,
                                 owner_span) ||
          owner_span == 0U) {
        return {StatusCode::invalid_plan, kArenaOwner};
      }
      std::size_t global_offset = 0U;
      if (!detail::checked_align(candidate.total_doubles_,
                                 detail::kDoublesPerCacheLine,
                                 global_offset)) {
        return {StatusCode::invalid_plan, kArenaOverflow};
      }
      candidate.owners_.push_back(
          ArenaOwnerLayout{owner_id, global_offset, owner_span});
      if (!detail::checked_add(global_offset, owner_span,
                               candidate.total_doubles_)) {
        return {StatusCode::invalid_plan, kArenaOverflow};
      }
    }
    std::size_t total_bytes = 0U;
    if (!detail::checked_multiply(candidate.total_doubles_, sizeof(double),
                                  total_bytes) ||
        total_bytes >
            static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())) {
      return {StatusCode::invalid_plan, kArenaOverflow};
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kArenaRequest};
  }
}

const ArenaFieldLayout* ArenaLayout::field(FieldId id) const noexcept {
  const std::size_t index = static_cast<std::size_t>(id);
  return index < fields_.size() ? &fields_[index] : nullptr;
}

const ArenaOwnerLayout* ArenaLayout::owner(std::size_t index) const noexcept {
  return index < owners_.size() ? &owners_[index] : nullptr;
}

FieldStorage::~FieldStorage() noexcept { release(); }

FieldStorage::FieldStorage(FieldStorage&& other) noexcept
    : fields_(std::move(other.fields_)),
      owners_(std::move(other.owners_)),
      counters_(other.counters_),
      identity_(other.identity_) {
  other.counters_ = {};
  other.identity_ = 0U;
}

FieldStorage& FieldStorage::operator=(FieldStorage&& other) noexcept {
  if (this != &other) {
    release();
    fields_ = std::move(other.fields_);
    owners_ = std::move(other.owners_);
    counters_ = other.counters_;
    identity_ = other.identity_;
    other.counters_ = {};
    other.identity_ = 0U;
  }
  return *this;
}

void FieldStorage::release() noexcept {
  for (OwnerAllocation& owner : owners_) {
    ::operator delete(owner.data,
                      std::align_val_t{detail::kCacheLineBytes});
    owner.data = nullptr;
    owner.doubles = 0U;
  }
  owners_.clear();
  fields_.clear();
  counters_ = {};
  identity_ = 0U;
}

Status FieldStorage::allocate(const ArenaLayout& layout, FieldStorage& out) {
  if (layout.field_count() == 0U || layout.owner_count() == 0U) {
    return {StatusCode::invalid_plan, kArenaSchema};
  }
  try {
    FieldStorage candidate;
    candidate.identity_ = detail::issue_identity();
    if (candidate.identity_ == 0U) {
      return {StatusCode::invalid_plan, kArenaRevision};
    }
    candidate.fields_ = layout.fields_;
    candidate.owners_.reserve(layout.owners_.size());
    for (const ArenaOwnerLayout& owner_layout : layout.owners_) {
      std::size_t bytes = 0U;
      if (!detail::checked_multiply(owner_layout.span_doubles, sizeof(double),
                                    bytes) || bytes == 0U) {
        return {StatusCode::invalid_plan, kArenaOverflow};
      }
      double* data = static_cast<double*>(::operator new(
          bytes, std::align_val_t{detail::kCacheLineBytes}));
      std::uninitialized_fill_n(data, owner_layout.span_doubles, 0.0);
      try {
        candidate.owners_.push_back(
            OwnerAllocation{owner_layout.owner, data,
                            owner_layout.span_doubles});
      } catch (...) {
        ::operator delete(data,
                          std::align_val_t{detail::kCacheLineBytes});
        throw;
      }
      ++candidate.counters_.aligned_payload_allocations;
      candidate.counters_.aligned_payload_bytes +=
          static_cast<std::uint64_t>(bytes);
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kArenaOwner};
  }
}

const ArenaFieldLayout* FieldStorage::field_layout(FieldId field) const noexcept {
  const std::size_t index = static_cast<std::size_t>(field);
  return index < fields_.size() ? &fields_[index] : nullptr;
}

Status FieldStorage::view(FieldId field, std::size_t replica,
                          const RevisionSet& revisions,
                          FieldView& out) noexcept {
  return view(field, replica, revisions.token(field), revisions.identity_, out);
}

Status FieldStorage::view(FieldId field, std::size_t replica,
                          const RevisionSet& revisions,
                          ConstFieldView& out) const noexcept {
  return view(field, replica, revisions.token(field), revisions.identity_, out);
}

Status FieldStorage::view(FieldId field, std::size_t replica,
                          RevisionToken revision,
                          RevisionDomainIdentity revision_domain,
                          FieldView& out) noexcept {
  const ArenaFieldLayout* layout = field_layout(field);
  if (layout == nullptr || layout->owner_index >= owners_.size() ||
      owners_[layout->owner_index].owner != layout->owner) {
    return {StatusCode::invalid_plan, kArenaView};
  }
  return make_view(owners_[layout->owner_index].data, *layout, replica,
                   revision, identity_, revision_domain, out);
}

Status FieldStorage::view(FieldId field, std::size_t replica,
                          RevisionToken revision,
                          RevisionDomainIdentity revision_domain,
                          ConstFieldView& out) const noexcept {
  const ArenaFieldLayout* layout = field_layout(field);
  if (layout == nullptr || layout->owner_index >= owners_.size() ||
      owners_[layout->owner_index].owner != layout->owner) {
    return {StatusCode::invalid_plan, kArenaView};
  }
  return make_view(static_cast<const double*>(
                       owners_[layout->owner_index].data),
                   *layout, replica, revision, identity_, revision_domain,
                   out);
}

template <class T>
T* FieldStorage::checked_ptr_impl(
    const BasicFieldView<T>& view, Int3 index, std::uint8_t component,
    RevisionToken expected_revision,
    RevisionDomainIdentity expected_domain) const noexcept {
  const ArenaFieldLayout* layout = field_layout(view.field);
  if (layout == nullptr || layout->owner_index >= owners_.size() ||
      owners_[layout->owner_index].owner != layout->owner ||
      owners_[layout->owner_index].data == nullptr ||
      view.replica >= layout->replicas || expected_revision == 0U ||
      expected_domain == 0U || view.revision != expected_revision ||
      view.storage_identity != identity_ ||
      view.revision_domain != expected_domain ||
      view.interior.x != layout->interior.x ||
      view.interior.y != layout->interior.y ||
      view.interior.z != layout->interior.z ||
      view.ghosts.x != layout->ghosts.x ||
      view.ghosts.y != layout->ghosts.y ||
      view.ghosts.z != layout->ghosts.z ||
      view.components != layout->components || view.stride_y != layout->stride_y ||
      view.stride_z != layout->stride_z ||
      view.component_stride != layout->component_stride) {
    return nullptr;
  }
  std::size_t replica_offset = 0U;
  if (!detail::checked_multiply(view.replica,
                                layout->replica_stride_doubles,
                                replica_offset)) {
    return nullptr;
  }
  T* const expected_base =
      owners_[layout->owner_index].data + layout->offset_doubles +
      replica_offset;
  const auto x = static_cast<std::int64_t>(index.x);
  const auto y = static_cast<std::int64_t>(index.y);
  const auto z = static_cast<std::int64_t>(index.z);
  const auto gx = static_cast<std::int64_t>(layout->ghosts.x);
  const auto gy = static_cast<std::int64_t>(layout->ghosts.y);
  const auto gz = static_cast<std::int64_t>(layout->ghosts.z);
  if (view.base != expected_base || component >= layout->components ||
      x < -gx || y < -gy || z < -gz ||
      x >= static_cast<std::int64_t>(layout->interior.x) + gx ||
      y >= static_cast<std::int64_t>(layout->interior.y) + gy ||
      z >= static_cast<std::int64_t>(layout->interior.z) + gz) {
    return nullptr;
  }
  const auto offset = static_cast<std::ptrdiff_t>(index.x) +
                      static_cast<std::ptrdiff_t>(index.y) *
                          static_cast<std::ptrdiff_t>(layout->stride_y) +
                      static_cast<std::ptrdiff_t>(index.z) *
                          static_cast<std::ptrdiff_t>(layout->stride_z) +
                      static_cast<std::ptrdiff_t>(component) *
                          static_cast<std::ptrdiff_t>(layout->component_stride);
  return expected_base + offset;
}

double* FieldStorage::checked_ptr(const FieldView& view, Int3 index,
                                  std::uint8_t component,
                                  const RevisionSet& revisions) noexcept {
  return checked_ptr(view, index, component, revisions.token(view.field),
                     revisions.identity_);
}

const double* FieldStorage::checked_ptr(
    const ConstFieldView& view, Int3 index, std::uint8_t component,
    const RevisionSet& revisions) const noexcept {
  return checked_ptr(view, index, component, revisions.token(view.field),
                     revisions.identity_);
}

double* FieldStorage::checked_ptr(const FieldView& view, Int3 index,
                                  std::uint8_t component,
                                  RevisionToken expected_revision,
                                  RevisionDomainIdentity expected_domain) noexcept {
  return checked_ptr_impl(view, index, component, expected_revision,
                          expected_domain);
}

const double* FieldStorage::checked_ptr(
    const ConstFieldView& view, Int3 index, std::uint8_t component,
    RevisionToken expected_revision,
    RevisionDomainIdentity expected_domain) const noexcept {
  return checked_ptr_impl(view, index, component, expected_revision,
                          expected_domain);
}

StateLayers::StateLayers(StateLayers&& other) noexcept
    : storage_(std::move(other.storage_)),
      revisions_(std::move(other.revisions_)),
      state_revisions_{other.state_revisions_[0U],
                       other.state_revisions_[1U],
                       other.state_revisions_[2U]},
      state_fields_(std::move(other.state_fields_)),
      field_count_(other.field_count_),
      role_handles_{other.role_handles_[0U], other.role_handles_[1U],
                    other.role_handles_[2U]},
      next_revision_(other.next_revision_) {
  other.state_revisions_[0U] = 0U;
  other.state_revisions_[1U] = 0U;
  other.state_revisions_[2U] = 0U;
  other.field_count_ = 0U;
  other.role_handles_[0U] = 0U;
  other.role_handles_[1U] = 1U;
  other.role_handles_[2U] = 2U;
  other.next_revision_ = 1U;
}

StateLayers& StateLayers::operator=(StateLayers&& other) noexcept {
  if (this != &other) {
    storage_ = std::move(other.storage_);
    revisions_ = std::move(other.revisions_);
    state_fields_ = std::move(other.state_fields_);
    for (std::size_t role = 0U; role < kRoleCount; ++role) {
      state_revisions_[role] = other.state_revisions_[role];
      role_handles_[role] = other.role_handles_[role];
      other.state_revisions_[role] = 0U;
    }
    field_count_ = other.field_count_;
    next_revision_ = other.next_revision_;
    other.field_count_ = 0U;
    other.role_handles_[0U] = 0U;
    other.role_handles_[1U] = 1U;
    other.role_handles_[2U] = 2U;
    other.next_revision_ = 1U;
  }
  return *this;
}

Status StateLayers::allocate(const ArenaLayout& layout, StateLayers& out) {
  try {
    StateLayers candidate;
    const Status allocated = FieldStorage::allocate(layout, candidate.storage_);
    if (!allocated) {
      return allocated;
    }
    candidate.field_count_ = layout.field_count();
    std::size_t revision_count = 0U;
    if (!detail::checked_multiply(candidate.field_count_, kRoleCount,
                                  revision_count)) {
      return {StatusCode::invalid_plan, kArenaOverflow};
    }
    candidate.revisions_.assign(revision_count, RevisionToken{0U});
    candidate.state_fields_.assign(candidate.field_count_, 0U);
    std::size_t state_field_count = 0U;
    for (const ArenaFieldLayout& field : layout.fields_) {
      const std::size_t index = static_cast<std::size_t>(field.id);
      if (field.lifetime == FieldLifetime::state_layer) {
        if (field.replicas != kRoleCount) {
          return {StatusCode::invalid_plan, kArenaSchema};
        }
        candidate.state_fields_[index] = 1U;
        ++state_field_count;
        for (std::size_t replica = 0U; replica < kRoleCount; ++replica) {
          const RevisionToken issued = candidate.next_revision_++;
          candidate.revisions_[replica * candidate.field_count_ + index] =
              issued;
          candidate.state_revisions_[replica] = issued;
        }
      } else {
        if (field.replicas != 1U ||
            candidate.next_revision_ ==
                std::numeric_limits<RevisionToken>::max()) {
          return {StatusCode::invalid_plan, kArenaSchema};
        }
        candidate.revisions_[index] = candidate.next_revision_++;
      }
    }
    if (state_field_count == 0U) {
      return {StatusCode::invalid_plan, kArenaSchema};
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kArenaSchema};
  }
}

std::size_t StateLayers::role_index(StateRole role) const noexcept {
  return role_slot(role);
}

std::size_t StateLayers::handle(StateRole role) const noexcept {
  const std::size_t slot = role_index(role);
  return slot < kRoleCount ? role_handles_[slot] : kRoleCount;
}

RevisionToken StateLayers::revision(StateRole role, FieldId field) const noexcept {
  const std::size_t id = static_cast<std::size_t>(field);
  const std::size_t replica = handle(role);
  if (id >= field_count_ || id >= state_fields_.size() ||
      state_fields_[id] == 0U || replica >= kRoleCount) {
    return RevisionToken{0U};
  }
  return revisions_[replica * field_count_ + id];
}

RevisionToken StateLayers::state_revision(StateRole role) const noexcept {
  const std::size_t replica = handle(role);
  return replica < kRoleCount ? state_revisions_[replica]
                              : RevisionToken{0U};
}

Status StateLayers::view(StateRole role, FieldId field, FieldView& out) noexcept {
  const RevisionToken token = revision(role, field);
  return token == 0U
             ? Status{StatusCode::invalid_plan, kArenaView}
             : storage_.view(field, handle(role), token, storage_.identity_,
                             out);
}

Status StateLayers::view(StateRole role, FieldId field,
                         ConstFieldView& out) const noexcept {
  const RevisionToken token = revision(role, field);
  return token == 0U
             ? Status{StatusCode::invalid_plan, kArenaView}
             : storage_.view(field, handle(role), token, storage_.identity_,
                             out);
}

RevisionToken StateLayers::runtime_revision(FieldLifetime lifetime,
                                            FieldId field) const noexcept {
  const ArenaFieldLayout* const layout = storage_.field_layout(field);
  const std::size_t id = static_cast<std::size_t>(field);
  if (lifetime == FieldLifetime::state_layer || layout == nullptr ||
      layout->lifetime != lifetime || layout->replicas != 1U ||
      id >= field_count_ || id >= revisions_.size()) {
    return RevisionToken{0U};
  }
  return revisions_[id];
}

Status StateLayers::revise_runtime(FieldLifetime lifetime,
                                   FieldId field) noexcept {
  const ArenaFieldLayout* const layout = storage_.field_layout(field);
  const std::size_t id = static_cast<std::size_t>(field);
  if (lifetime == FieldLifetime::state_layer || layout == nullptr ||
      layout->lifetime != lifetime || layout->replicas != 1U ||
      id >= field_count_ || id >= revisions_.size() ||
      id >= state_fields_.size() || state_fields_[id] != 0U ||
      next_revision_ == 0U ||
      next_revision_ == std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kArenaRevision};
  }
  revisions_[id] = next_revision_++;
  return {};
}

Status StateLayers::runtime_view(FieldLifetime lifetime, FieldId field,
                                 FieldView& out) noexcept {
  const RevisionToken token = runtime_revision(lifetime, field);
  return token == 0U
             ? Status{StatusCode::invalid_plan, kArenaView}
             : storage_.view(field, 0U, token, storage_.identity_, out);
}

Status StateLayers::runtime_view(FieldLifetime lifetime, FieldId field,
                                 ConstFieldView& out) const noexcept {
  const RevisionToken token = runtime_revision(lifetime, field);
  return token == 0U
             ? Status{StatusCode::invalid_plan, kArenaView}
             : storage_.view(field, 0U, token, storage_.identity_, out);
}

double* StateLayers::checked_ptr(StateRole role, const FieldView& view,
                                 Int3 index,
                                 std::uint8_t component) noexcept {
  return storage_.checked_ptr(view, index, component,
                              revision(role, view.field),
                              storage_.identity_);
}

const double* StateLayers::checked_ptr(StateRole role,
                                       const ConstFieldView& view, Int3 index,
                                       std::uint8_t component) const noexcept {
  return storage_.checked_ptr(view, index, component,
                              revision(role, view.field),
                              storage_.identity_);
}

ConstFieldView as_const(FieldView view) noexcept {
  return ConstFieldView{view.base,
                        view.interior,
                        view.ghosts,
                        view.components,
                        view.stride_y,
                        view.stride_z,
                        view.component_stride,
                        view.replica,
                        view.field,
                        view.revision,
                        view.storage_identity,
                        view.revision_domain};
}

Status StateLayers::issue_revision(std::size_t replica, FieldId field) noexcept {
  const std::size_t id = static_cast<std::size_t>(field);
  if (replica >= kRoleCount || id >= field_count_ ||
      state_fields_[id] == 0U || next_revision_ == 0U ||
      next_revision_ == std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kArenaRevision};
  }
  revisions_[replica * field_count_ + id] = next_revision_;
  state_revisions_[replica] = next_revision_;
  ++next_revision_;
  return {};
}

void StateLayers::rotate_commit() noexcept {
  const std::size_t old_previous = role_handles_[1U];
  role_handles_[1U] = role_handles_[0U];
  role_handles_[0U] = role_handles_[2U];
  role_handles_[2U] = old_previous;
}

}  // namespace hundun::v04
