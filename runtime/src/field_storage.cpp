// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/field_storage.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "hundun/runtime/error.hpp"
#include "field_epoch_test_access.hpp"

namespace hundun::runtime {
namespace detail {

class FieldEpochControl final {
 public:
  std::atomic<std::uint64_t> generation{1U};
  std::atomic<bool> alive{true};
};

std::uint64_t field_epoch_generation(
    const std::shared_ptr<FieldEpochControl> &epoch) noexcept {
  if (!epoch) {
    return 0U;
  }
  return epoch->generation.load(std::memory_order_acquire);
}

void validate_field_epoch(
    const std::shared_ptr<const FieldEpochControl> &epoch,
    std::uint64_t captured_generation) {
  if (!epoch || !epoch->alive.load(std::memory_order_acquire)) {
    throw Error("field view owner is no longer alive");
  }
  if (epoch->generation.load(std::memory_order_acquire) !=
      captured_generation) {
    throw Error("field view generation is stale");
  }
}

}  // namespace detail
namespace {

void invalidate_epoch(
    const std::shared_ptr<detail::FieldEpochControl> &epoch) noexcept {
  if (epoch) {
    epoch->alive.store(false, std::memory_order_release);
  }
}

void advance_epoch(const std::shared_ptr<detail::FieldEpochControl> &epoch) {
  if (!epoch || !epoch->alive.load(std::memory_order_acquire)) {
    throw Error("field storage has no active generation");
  }
  const auto generation =
      epoch->generation.load(std::memory_order_acquire);
  if (generation == std::numeric_limits<std::uint64_t>::max()) {
    throw Error("field storage generation would wrap");
  }
  epoch->generation.store(generation + 1U, std::memory_order_release);
}

std::size_t checked_add(std::size_t left, std::size_t right,
                        const char *quantity) {
  constexpr auto limit = std::numeric_limits<std::size_t>::max();
  if (left > limit - right) {
    throw Error(std::string("field storage ") + quantity + " overflow");
  }
  return left + right;
}

std::size_t checked_multiply(std::size_t left, std::size_t right,
                             const char *quantity) {
  constexpr auto limit = std::numeric_limits<std::size_t>::max();
  if (right != 0U && left > limit / right) {
    throw Error(std::string("field storage ") + quantity + " overflow");
  }
  return left * right;
}

std::size_t scalar_size(ScalarType scalar_type) {
  switch (scalar_type) {
    case ScalarType::float64:
      return sizeof(double);
    case ScalarType::int32:
      return sizeof(std::int32_t);
    case ScalarType::uint8:
      return sizeof(std::uint8_t);
  }
  throw Error("field storage received an unrecognized scalar type");
}

std::size_t scalar_alignment(ScalarType scalar_type) {
  switch (scalar_type) {
    case ScalarType::float64:
      return alignof(double);
    case ScalarType::int32:
      return alignof(std::int32_t);
    case ScalarType::uint8:
      return alignof(std::uint8_t);
  }
  throw Error("field storage received an unrecognized scalar type");
}

template <class T>
void begin_scalar_lifetimes(std::byte *storage, std::size_t element_count) {
  static_assert(std::is_trivially_destructible_v<T>,
                "field scalars must not require destructor dispatch");
  for (std::size_t index = 0; index < element_count; ++index) {
    ::new (static_cast<void *>(storage + index * sizeof(T))) T{};
  }
}

void begin_scalar_lifetimes(ScalarType scalar_type, std::byte *storage,
                            std::size_t element_count) {
  // The byte vector provides raw storage. std::align establishes a suitably
  // aligned base, then placement new starts every supported scalar lifetime
  // under C++17. All supported scalars are trivially destructible.
  switch (scalar_type) {
    case ScalarType::float64:
      begin_scalar_lifetimes<double>(storage, element_count);
      return;
    case ScalarType::int32:
      begin_scalar_lifetimes<std::int32_t>(storage, element_count);
      return;
    case ScalarType::uint8:
      begin_scalar_lifetimes<std::uint8_t>(storage, element_count);
      return;
  }
  throw Error("field storage received an unrecognized scalar type");
}

struct AllocationPlan {
  FunctionSpace space{};
  ScalarType scalar_type{};
  int ghost_width{};
  std::uint32_t components{};
  std::size_t face_count{};
  std::size_t x_stride{};
  std::size_t y_stride{};
  std::size_t z_stride{};
  std::size_t element_count{};
  std::size_t byte_count{};
  std::size_t alignment{};
  std::size_t storage_byte_count{};
};

void validate_component_index_range(const FieldDescriptor &descriptor_value) {
  constexpr auto max_indexable_components =
      static_cast<std::uintmax_t>(std::numeric_limits<int>::max()) + 1U;
  if (static_cast<std::uintmax_t>(descriptor_value.components) >
      max_indexable_components) {
    throw Error("field storage component count exceeds the view index range");
  }
}

AllocationPlan make_cell_plan(const FieldDescriptor &descriptor_value,
                              Int3 interior_extent) {
  validate_component_index_range(descriptor_value);

  constexpr auto max_indexable_axis_span =
      static_cast<std::uintmax_t>(std::numeric_limits<int>::max()) + 1U;
  const auto coordinate_ghost_width =
      static_cast<std::uintmax_t>(descriptor_value.ghost_width);
  if (static_cast<std::uintmax_t>(interior_extent.x) + coordinate_ghost_width >
          max_indexable_axis_span ||
      static_cast<std::uintmax_t>(interior_extent.y) + coordinate_ghost_width >
          max_indexable_axis_span ||
      static_cast<std::uintmax_t>(interior_extent.z) + coordinate_ghost_width >
          max_indexable_axis_span) {
    throw Error(
        "field storage upper ghost coordinate exceeds the view index range");
  }

  const auto ghost_width =
      static_cast<std::size_t>(descriptor_value.ghost_width);
  const auto double_ghost =
      checked_multiply(ghost_width, 2U, "padded-axis width");
  const auto padded_x = checked_add(static_cast<std::size_t>(interior_extent.x),
                                    double_ghost, "padded x extent");
  const auto padded_y = checked_add(static_cast<std::size_t>(interior_extent.y),
                                    double_ghost, "padded y extent");
  const auto padded_z = checked_add(static_cast<std::size_t>(interior_extent.z),
                                    double_ghost, "padded z extent");
  const auto components = static_cast<std::size_t>(descriptor_value.components);

  AllocationPlan plan;
  plan.space = descriptor_value.space;
  plan.scalar_type = descriptor_value.scalar_type;
  plan.ghost_width = descriptor_value.ghost_width;
  plan.components = descriptor_value.components;
  plan.x_stride = components;
  plan.y_stride = checked_multiply(padded_x, plan.x_stride, "y stride");
  plan.z_stride = checked_multiply(padded_y, plan.y_stride, "z stride");
  plan.element_count =
      checked_multiply(padded_z, plan.z_stride, "element count");
  plan.byte_count = checked_multiply(
      plan.element_count, scalar_size(plan.scalar_type), "byte count");
  plan.alignment = scalar_alignment(plan.scalar_type);
  plan.storage_byte_count =
      checked_add(plan.byte_count, plan.alignment - 1U, "aligned byte count");
  return plan;
}

AllocationPlan make_face_plan(const FieldDescriptor &descriptor_value,
                              std::size_t face_count) {
  if (face_count == 0U) {
    throw Error("face field storage requires a positive face count");
  }
  if (descriptor_value.ghost_width != 0) {
    throw Error("face_value fields require zero ghost width");
  }
  validate_component_index_range(descriptor_value);

  AllocationPlan plan;
  plan.space = descriptor_value.space;
  plan.scalar_type = descriptor_value.scalar_type;
  plan.ghost_width = descriptor_value.ghost_width;
  plan.components = descriptor_value.components;
  plan.face_count = face_count;
  plan.x_stride = static_cast<std::size_t>(descriptor_value.components);
  plan.element_count = checked_multiply(
      face_count, plan.x_stride, "face element count");
  plan.byte_count = checked_multiply(
      plan.element_count, scalar_size(plan.scalar_type), "face byte count");
  plan.alignment = scalar_alignment(plan.scalar_type);
  plan.storage_byte_count = checked_add(
      plan.byte_count, plan.alignment - 1U, "aligned face byte count");
  return plan;
}

}  // namespace

FieldStorage::FieldStorage(const FieldRegistry &registry, Int3 interior_extent)
    : FieldStorage(registry, FieldLayoutSet{interior_extent, 0U},
                   ConstructionMode::legacy_cell_only) {}

FieldStorage::FieldStorage(const FieldRegistry &registry,
                           FieldLayoutSet layout_set)
    : FieldStorage(registry, layout_set, ConstructionMode::cell_and_face) {}

FieldStorage::FieldStorage(const FieldRegistry &registry,
                           FieldLayoutSet layout_set,
                           ConstructionMode construction_mode)
    : interior_extent_(layout_set.cell_interior_extent),
      layout_set_(layout_set) {
  if (!registry.frozen()) {
    throw Error("field storage requires a frozen registry");
  }
  if (interior_extent_.x <= 0 || interior_extent_.y <= 0 ||
      interior_extent_.z <= 0) {
    throw Error("field storage interior extents must be positive");
  }

  std::vector<AllocationPlan> plans;
  plans.reserve(registry.size());
  for (std::size_t index = 0; index < registry.size(); ++index) {
    const auto &descriptor_value =
        registry.descriptor(static_cast<FieldId>(index));
    if (construction_mode == ConstructionMode::legacy_cell_only) {
      if (descriptor_value.space != FunctionSpace::cell_average) {
        throw Error(
            "Stage 1 field storage supports cell_average fields only");
      }
      plans.push_back(make_cell_plan(descriptor_value, interior_extent_));
      continue;
    }

    switch (descriptor_value.space) {
      case FunctionSpace::cell_average:
        plans.push_back(make_cell_plan(descriptor_value, interior_extent_));
        break;
      case FunctionSpace::face_value:
        plans.push_back(
            make_face_plan(descriptor_value, layout_set_.face_count));
        break;
      case FunctionSpace::vertex_value:
      case FunctionSpace::element_dof:
      case FunctionSpace::quadrature_point:
      case FunctionSpace::particle:
        throw Error("field storage does not support this function space");
    }
  }

  auto entries = std::make_unique<std::vector<Entry>>();
  entries->reserve(plans.size());
  for (const auto &plan : plans) {
    Entry field;
    field.space = plan.space;
    field.scalar_type = plan.scalar_type;
    field.ghost_width = plan.ghost_width;
    field.components = plan.components;
    field.face_count = plan.face_count;
    field.x_stride = plan.x_stride;
    field.y_stride = plan.y_stride;
    field.z_stride = plan.z_stride;
    if (plan.storage_byte_count > field.bytes.max_size()) {
      throw Error("field storage byte count exceeds vector capacity");
    }
    field.bytes.resize(plan.storage_byte_count);

    void *candidate = static_cast<void *>(field.bytes.data());
    std::size_t available = field.bytes.size();
    void *aligned =
        std::align(plan.alignment, plan.byte_count, candidate, available);
    if (aligned == nullptr) {
      throw Error("field storage could not align its scalar buffer");
    }
    auto *aligned_bytes = static_cast<std::byte *>(aligned);
    field.data_offset =
        static_cast<std::size_t>(aligned_bytes - field.bytes.data());
    begin_scalar_lifetimes(plan.scalar_type, aligned_bytes, plan.element_count);
    entries->push_back(std::move(field));
  }
  epoch_ = std::make_shared<detail::FieldEpochControl>();
  entries_ = std::move(entries);
}

FieldStorage::~FieldStorage() noexcept { invalidate_epoch(epoch_); }

FieldStorage::FieldStorage(FieldStorage &&other) noexcept
    : interior_extent_(other.interior_extent_),
      layout_set_(other.layout_set_),
      epoch_(std::move(other.epoch_)),
      entries_(std::move(other.entries_)) {
  other.interior_extent_ = Int3{};
  other.layout_set_ = FieldLayoutSet{};
}

FieldStorage &FieldStorage::operator=(FieldStorage &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  invalidate_epoch(epoch_);
  entries_.reset();
  epoch_.reset();

  interior_extent_ = other.interior_extent_;
  layout_set_ = other.layout_set_;
  epoch_ = std::move(other.epoch_);
  entries_ = std::move(other.entries_);
  other.interior_extent_ = Int3{};
  other.layout_set_ = FieldLayoutSet{};
  return *this;
}

Int3 FieldStorage::interior_extent() const noexcept { return interior_extent_; }

FieldLayoutSet FieldStorage::layout_set() const noexcept { return layout_set_; }

void FieldStorage::begin_rebuild() { advance_epoch(epoch_); }

void FieldStorage::begin_repartition() { advance_epoch(epoch_); }

void FieldStorage::begin_restart_v2_read_transaction() {
  advance_epoch(epoch_);
}

FieldStorage::Entry &FieldStorage::entry(FieldId id) {
  const auto index = static_cast<std::size_t>(id);
  if (!entries_ || index >= entries_->size()) {
    throw Error("field storage ID is out of bounds");
  }
  return (*entries_)[index];
}

const FieldStorage::Entry &FieldStorage::entry(FieldId id) const {
  const auto index = static_cast<std::size_t>(id);
  if (!entries_ || index >= entries_->size()) {
    throw Error("field storage ID is out of bounds");
  }
  return (*entries_)[index];
}

std::size_t FieldStorage::field_count() const noexcept {
  return entries_ ? entries_->size() : 0U;
}

std::uint64_t detail::FieldEpochTestAccess::generation(
    const FieldStorage &storage) noexcept {
  return detail::field_epoch_generation(storage.epoch_);
}

void detail::FieldEpochTestAccess::force_generation(
    FieldStorage &storage, std::uint64_t generation) noexcept {
  if (storage.epoch_) {
    storage.epoch_->generation.store(generation, std::memory_order_release);
  }
}

}  // namespace hundun::runtime
