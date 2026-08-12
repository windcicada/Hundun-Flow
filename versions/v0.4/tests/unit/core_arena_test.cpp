// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"
#include "hundun/v04_field.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

using hundun::v04::ArenaFieldLayout;
using hundun::v04::ArenaFieldRequest;
using hundun::v04::ArenaLayout;
using hundun::v04::ArenaOwnerLayout;
using hundun::v04::FieldId;
using hundun::v04::FieldLifetime;
using hundun::v04::FieldPlacement;
using hundun::v04::FieldRegistry;
using hundun::v04::FieldSchema;
using hundun::v04::FieldStorage;
using hundun::v04::FieldView;
using hundun::v04::Int3;
using hundun::v04::RevisionSet;
using hundun::v04::Span;
using hundun::v04::Status;
using hundun::v04::StatusCode;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool same(const Int3& left, const Int3& right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(const ArenaFieldLayout& left, const ArenaFieldLayout& right) {
  return left.id == right.id && same(left.interior, right.interior) &&
         same(left.ghosts, right.ghosts) &&
         left.components == right.components &&
         left.stride_y == right.stride_y &&
         left.stride_z == right.stride_z &&
         left.component_stride == right.component_stride &&
         left.replicas == right.replicas &&
         left.replica_stride_doubles == right.replica_stride_doubles &&
         left.owner == right.owner && left.owner_index == right.owner_index &&
         left.raw_offset_doubles == right.raw_offset_doubles &&
         left.offset_doubles == right.offset_doubles &&
         left.span_doubles == right.span_doubles &&
         left.lifetime == right.lifetime;
}

bool same(const ArenaOwnerLayout& left, const ArenaOwnerLayout& right) {
  return left.owner == right.owner &&
         left.offset_doubles == right.offset_doubles &&
         left.span_doubles == right.span_doubles;
}

bool make_schema(FieldSchema& schema, FieldId& velocity, FieldId& pressure,
                 FieldId& workspace) {
  FieldRegistry registry;
  bool passed = true;
  passed &= expect(
      static_cast<bool>(registry.declare_field("velocity", 3U, 2U, velocity)),
      "velocity field declaration succeeds");
  passed &= expect(
      static_cast<bool>(registry.declare_field("pressure", 1U, 1U, pressure)),
      "pressure field declaration succeeds");
  passed &= expect(static_cast<bool>(
                       registry.declare_field("operator_workspace", 2U, 0U,
                                              workspace)),
                   "workspace field declaration succeeds");
  passed &= expect(static_cast<bool>(registry.freeze_for_test(schema)),
                   "synthetic schema freezes");
  return passed;
}

bool test_rejects_invalid_and_overflow_plans(const FieldSchema& schema,
                                             FieldId velocity,
                                             FieldId pressure,
                                             FieldId workspace) {
  bool passed = true;
  ArenaLayout ignored;

  const std::array duplicate{
      ArenaFieldRequest{velocity, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{velocity, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{workspace, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::persistent_workspace},
  };
  const Status duplicate_status = ArenaLayout::compile(
      schema, Span<const ArenaFieldRequest>{duplicate.data(), duplicate.size()},
      ignored);
  passed &= expect(duplicate_status.code == StatusCode::invalid_plan,
                   "duplicate field requests are rejected");

  const std::array missing{
      ArenaFieldRequest{velocity, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{pressure, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
  };
  const Status missing_status = ArenaLayout::compile(
      schema, Span<const ArenaFieldRequest>{missing.data(), missing.size()},
      ignored);
  passed &= expect(missing_status.code == StatusCode::invalid_plan,
                   "a request set missing a schema field is rejected");

  const std::array nonpositive{
      ArenaFieldRequest{velocity, Int3{13, 0, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{pressure, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{workspace, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::persistent_workspace},
  };
  const Status nonpositive_status = ArenaLayout::compile(
      schema,
      Span<const ArenaFieldRequest>{nonpositive.data(), nonpositive.size()},
      ignored);
  passed &= expect(nonpositive_status.code == StatusCode::invalid_plan,
                   "non-positive interior extents are rejected");

  constexpr std::int32_t largest = std::numeric_limits<std::int32_t>::max();
  const std::array overflow{
      ArenaFieldRequest{velocity, Int3{largest, largest, largest},
                        FieldPlacement{0U}, FieldLifetime::state_layer},
      ArenaFieldRequest{pressure, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{workspace, Int3{13, 5, 3}, FieldPlacement{0U},
                        FieldLifetime::persistent_workspace},
  };
  const Status overflow_status = ArenaLayout::compile(
      schema, Span<const ArenaFieldRequest>{overflow.data(), overflow.size()},
      ignored);
  passed &= expect(overflow_status.code == StatusCode::invalid_plan,
                   "overflowing padded SoA sizes are rejected before allocation");
  return passed;
}

bool test_layout_and_storage(const FieldSchema& schema, FieldId velocity,
                             FieldId pressure, FieldId workspace) {
  const std::array requests{
      ArenaFieldRequest{velocity, Int3{13, 5, 3}, FieldPlacement{3U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{pressure, Int3{13, 5, 3}, FieldPlacement{1U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{workspace, Int3{9, 4, 2}, FieldPlacement{3U},
                        FieldLifetime::persistent_workspace},
  };

  ArenaLayout layout;
  bool passed = expect(
      static_cast<bool>(ArenaLayout::compile(
          schema, Span<const ArenaFieldRequest>{requests.data(), requests.size()},
          layout)),
      "valid arena plan compiles");
  passed &= expect(layout.field_count() == requests.size(),
                   "layout records every field exactly once");
  passed &= expect(layout.owner_count() == 2U,
                   "layout records one aggregate range per active owner");
  passed &= expect(layout.total_doubles() > 0U,
                   "layout publishes a checked aggregate size");

  ArenaLayout repeated;
  passed &= expect(
      static_cast<bool>(ArenaLayout::compile(
          schema, Span<const ArenaFieldRequest>{requests.data(), requests.size()},
          repeated)),
      "the same plan compiles repeatedly");
  passed &= expect(repeated.field_count() == layout.field_count() &&
                       repeated.owner_count() == layout.owner_count() &&
                       repeated.total_doubles() == layout.total_doubles(),
                   "aggregate layout is deterministic");

  for (std::size_t index = 0; index < layout.field_count(); ++index) {
    const ArenaFieldLayout* field = layout.field(requests[index].id);
    const ArenaFieldLayout* again = repeated.field(requests[index].id);
    passed &= expect(field != nullptr && again != nullptr,
                     "compiled field is addressable by stable FieldId");
    if (field != nullptr && again != nullptr) {
      passed &= expect(same(*field, *again),
                       "per-field offsets and lifetimes are deterministic");
      passed &= expect(field->stride_y % 8U == 0U,
                       "every y row is padded to eight doubles");
      passed &= expect(field->stride_z % field->stride_y == 0U,
                       "z stride is an integral padded-row count");
      passed &= expect(field->component_stride % field->stride_y == 0U,
                       "each SoA component occupies complete padded rows");
      passed &= expect(field->offset_doubles % 8U == 0U,
                       "each field interior base is 64-byte aligned");
      passed &= expect(field->raw_offset_doubles % 8U == 0U &&
                           field->raw_offset_doubles <= field->offset_doubles,
                       "raw and interior field offsets are distinct and aligned");
      passed &= expect(field->span_doubles >= field->component_stride,
                       "field span contains its ghosted component storage");
      const std::size_t expected_replicas =
          field->lifetime == FieldLifetime::state_layer ? 3U : 1U;
      passed &= expect(field->replicas == expected_replicas,
                       "state layers have three replicas and workspaces one");
      passed &= expect(field->replica_stride_doubles % 8U == 0U &&
                           field->span_doubles ==
                               field->replica_stride_doubles * field->replicas,
                       "replica bases are aligned and total span is exact");
    }
  }

  for (std::size_t index = 0; index < layout.owner_count(); ++index) {
    const ArenaOwnerLayout* owner = layout.owner(index);
    const ArenaOwnerLayout* again = repeated.owner(index);
    passed &= expect(owner != nullptr && again != nullptr,
                     "active owner aggregate is addressable");
    if (owner != nullptr && again != nullptr) {
      passed &= expect(same(*owner, *again),
                       "owner aggregate offsets are deterministic");
      passed &= expect(owner->span_doubles > 0U,
                       "active owner aggregate is non-empty");
      if (index > 0U) {
        const ArenaOwnerLayout* previous = layout.owner(index - 1U);
        passed &= expect(previous != nullptr && previous->owner < owner->owner,
                         "active owner aggregates have stable owner ordering");
        passed &= expect(previous != nullptr &&
                             previous->offset_doubles +
                                     previous->span_doubles <=
                                 owner->offset_doubles,
                         "owner aggregate ranges never overlap");
      }
    }
  }

  const ArenaFieldLayout* u_layout = layout.field(velocity);
  const ArenaFieldLayout* p_layout = layout.field(pressure);
  const ArenaFieldLayout* w_layout = layout.field(workspace);
  passed &= expect(u_layout != nullptr && u_layout->owner == 3U &&
                       u_layout->lifetime == FieldLifetime::state_layer,
                   "velocity placement and lifetime survive planning");
  passed &= expect(p_layout != nullptr && p_layout->owner == 1U &&
                       p_layout->lifetime == FieldLifetime::state_layer,
                   "pressure placement and lifetime survive planning");
  passed &= expect(w_layout != nullptr && w_layout->owner == 3U &&
                       w_layout->lifetime ==
                           FieldLifetime::persistent_workspace,
                   "workspace placement and lifetime survive planning");

  FieldStorage storage;
  passed &= expect(static_cast<bool>(FieldStorage::allocate(layout, storage)),
                   "field storage allocates from a complete plan");
  passed &= expect(storage.counters().aligned_payload_allocations ==
                       layout.owner_count(),
                   "storage performs exactly one payload allocation per owner");
  passed &= expect(storage.counters().aligned_payload_bytes ==
                       layout.total_doubles() * sizeof(double),
                   "payload-byte counter matches the planned owner aggregates");
  passed &= expect(storage.counters().whole_field_copies == 0U,
                   "cold storage construction copies no whole field");

  RevisionSet revisions;
  passed &= expect(static_cast<bool>(RevisionSet::create(schema, revisions)),
                   "revision slots are created before views");
  FieldView u;
  FieldView p;
  passed &= expect(static_cast<bool>(
                       storage.view(velocity, 0U, revisions, u)),
                   "velocity borrowed view is created");
  passed &= expect(static_cast<bool>(
                       storage.view(pressure, 0U, revisions, p)),
                   "pressure borrowed view is created");
  static_assert(std::is_trivially_copyable_v<FieldView>,
                "borrowed views must not own storage or revision containers");
  static_assert(std::is_trivially_copyable_v<hundun::v04::ConstFieldView>,
                "const borrowed views must remain non-owning");

  passed &= expect(reinterpret_cast<std::uintptr_t>(u.base) % 64U == 0U &&
                       reinterpret_cast<std::uintptr_t>(p.base) % 64U == 0U,
                   "each borrowed field interior base is 64-byte aligned");
  passed &= expect(&u.unchecked(Int3{1, 0, 0}, 0U) -
                           &u.unchecked(Int3{0, 0, 0}, 0U) ==
                       1,
                   "x is unit stride");
  passed &= expect(static_cast<std::size_t>(
                       &u.unchecked(Int3{0, 1, 0}, 0U) -
                       &u.unchecked(Int3{0, 0, 0}, 0U)) == u.stride_y,
                   "y address arithmetic uses the padded row stride");
  passed &= expect(static_cast<std::size_t>(
                       &u.unchecked(Int3{0, 0, 1}, 0U) -
                       &u.unchecked(Int3{0, 0, 0}, 0U)) == u.stride_z,
                   "z address arithmetic uses the planned plane stride");
  passed &= expect(static_cast<std::size_t>(
                       &u.unchecked(Int3{0, 0, 0}, 1U) -
                       &u.unchecked(Int3{0, 0, 0}, 0U)) ==
                       u.component_stride,
                   "SoA components occupy separate contiguous spans");

  double* low_ghost = storage.checked_ptr(u, Int3{-2, -2, -2}, 0U, revisions);
  double* high_ghost =
      storage.checked_ptr(u, Int3{14, 6, 4}, 2U, revisions);
  passed &= expect(low_ghost != nullptr && high_ghost != nullptr,
                   "interior and all ghost layers share the owner allocation");
  passed &= expect(
      storage.checked_ptr(u, Int3{-3, 0, 0}, 0U, revisions) == nullptr &&
          storage.checked_ptr(u, Int3{15, 0, 0}, 0U, revisions) == nullptr &&
          storage.checked_ptr(u, Int3{0, 0, 0}, 3U, revisions) == nullptr,
                   "checked access rejects invalid logical/component indices");
  const auto const_u = hundun::v04::as_const(u);
  passed &= expect(const_u.base == u.base &&
                       storage.checked_ptr(const_u, Int3{0, 0, 0}, 0U,
                                           revisions) != nullptr,
                   "mutable views convert to read-only borrowed views");

  u.unchecked(Int3{0, 0, 0}, 0U) = 9.0;
  const FieldView stale = u;
  revisions.revise(velocity);
  passed &= expect(storage.checked_ptr(stale, Int3{0, 0, 0}, 0U,
                                       revisions) == nullptr,
                   "checked access rejects a stale captured revision");
  FieldView refreshed;
  passed &= expect(static_cast<bool>(
                       storage.view(velocity, 0U, revisions, refreshed)) &&
                       storage.checked_ptr(refreshed, Int3{0, 0, 0}, 0U,
                                           revisions) != nullptr,
                   "a view rebound to the live revision is valid");

  FieldStorage other_storage;
  passed &= expect(static_cast<bool>(FieldStorage::allocate(layout,
                                                             other_storage)),
                   "a second storage instance allocates independently");
  passed &= expect(other_storage.checked_ptr(refreshed, Int3{0, 0, 0}, 0U,
                                             revisions) == nullptr,
                   "checked access rejects a view from another storage");
  RevisionSet other_revisions;
  passed &= expect(static_cast<bool>(RevisionSet::create(schema,
                                                          other_revisions)) &&
                       storage.checked_ptr(refreshed, Int3{0, 0, 0}, 0U,
                                           other_revisions) == nullptr,
                   "checked access rejects another revision domain");

  FieldView forged = refreshed;
  forged.base = nullptr;
  passed &= expect(storage.checked_ptr(forged, Int3{0, 0, 0}, 0U,
                                       revisions) == nullptr,
                   "checked access rejects a forged null base");
  forged = refreshed;
  ++forged.stride_y;
  passed &= expect(storage.checked_ptr(forged, Int3{0, 0, 0}, 0U,
                                       revisions) == nullptr,
                   "checked access rejects forged layout metadata");

  FieldStorage moved_target;
  passed &= expect(static_cast<bool>(FieldStorage::allocate(layout,
                                                             moved_target)),
                   "move-assignment target storage allocates");
  FieldView old_target;
  passed &= expect(static_cast<bool>(moved_target.view(
                       velocity, 0U, revisions, old_target)),
                   "move-assignment target view is captured");
  moved_target = std::move(other_storage);
  passed &= expect(moved_target.checked_ptr(old_target, Int3{0, 0, 0}, 0U,
                                            revisions) == nullptr,
                   "move assignment rejects a view of released target storage");
  return passed;
}

}  // namespace

int main() {
  FieldSchema schema;
  FieldId velocity{};
  FieldId pressure{};
  FieldId workspace{};
  bool passed = make_schema(schema, velocity, pressure, workspace);
  passed &= test_rejects_invalid_and_overflow_plans(
      schema, velocity, pressure, workspace);
  passed &= test_layout_and_storage(schema, velocity, pressure, workspace);
  return passed ? 0 : 1;
}
