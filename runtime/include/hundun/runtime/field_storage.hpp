// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_view.hpp"
#include "hundun/runtime/types.hpp"

namespace hundun::flow {
class IdealGasClosure;
namespace detail {
struct FlowStateCheckpointAccess;
}
}

namespace hundun::runtime {

class HaloExchange;

namespace detail {

class FieldEpochControl;
struct FieldEpochTestAccess;

template <class T>
struct FieldScalarType;

template <>
struct FieldScalarType<double> {
  static constexpr ScalarType value = ScalarType::float64;
};

template <>
struct FieldScalarType<std::int32_t> {
  static constexpr ScalarType value = ScalarType::int32;
};

template <>
struct FieldScalarType<std::uint8_t> {
  static constexpr ScalarType value = ScalarType::uint8;
};

}  // namespace detail

struct FieldLayoutSet final {
  Int3 cell_interior_extent;
  std::size_t face_count{};
};

class FieldStorage final {
 public:
  FieldStorage(const FieldRegistry &, Int3 interior_extent);
  FieldStorage(const FieldRegistry &, FieldLayoutSet layout_set);
  ~FieldStorage() noexcept;
  FieldStorage(const FieldStorage &) = delete;
  FieldStorage &operator=(const FieldStorage &) = delete;
  FieldStorage(FieldStorage &&other) noexcept;
  FieldStorage &operator=(FieldStorage &&other) noexcept;

  Int3 interior_extent() const noexcept;
  FieldLayoutSet layout_set() const noexcept;

  void begin_rebuild();
  void begin_repartition();
  void begin_restart_v2_read_transaction();
  static void begin_restart_v2_read_transactions(
      FieldStorage *const *storages, std::size_t count);

  template <class T>
  FieldView<T> view(FieldId id);

  template <class T>
  FieldView<const T> view(FieldId id) const;

  template <class T>
  FieldView<const T> acquire_read(const FieldAccessPlan &plan, PhaseId phase,
                                  ActorId actor, FieldId id) const;

  template <class T>
  FieldView<T> acquire_write(const FieldAccessPlan &plan, PhaseId phase,
                             ActorId actor, FieldId id);

  template <class T>
  FaceFieldView<const T> acquire_face_read(const FieldAccessPlan &plan,
                                           PhaseId phase, ActorId actor,
                                           FieldId id) const;

  template <class T>
  FaceFieldView<T> acquire_face_write(const FieldAccessPlan &plan,
                                      PhaseId phase, ActorId actor,
                                      FieldId id);

 private:
  friend class HaloExchange;
  friend struct detail::FieldEpochTestAccess;
  friend class ::hundun::flow::IdealGasClosure;
  friend struct ::hundun::flow::detail::FlowStateCheckpointAccess;

  struct Entry {
    FunctionSpace space{};
    ScalarType scalar_type{};
    int ghost_width{};
    std::uint32_t components{};
    std::size_t face_count{};
    std::size_t x_stride{};
    std::size_t y_stride{};
    std::size_t z_stride{};
    std::size_t data_offset{};
    std::vector<std::byte> bytes;
  };

  Entry &entry(FieldId id);
  const Entry &entry(FieldId id) const;
  std::size_t field_count() const noexcept;
  static bool restart_v2_read_transactions_ready(
      FieldStorage *const *storages, std::size_t count) noexcept;
  void publish_validated_cell_interior_double(
      FieldId id, const double *candidate,
      std::size_t candidate_count) noexcept;

  enum class ConstructionMode { legacy_cell_only, cell_and_face };
  FieldStorage(const FieldRegistry &, FieldLayoutSet,
               ConstructionMode construction_mode);

  Int3 interior_extent_{};
  FieldLayoutSet layout_set_{};
  // Declaration order keeps the control state alive while entries_ is torn
  // down; checked views may retain it after the storage itself is gone.
  std::shared_ptr<detail::FieldEpochControl> epoch_;
  std::unique_ptr<std::vector<Entry>> entries_;
};

template <class T>
FieldView<T> FieldStorage::view(FieldId id) {
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  Entry &field = entry(id);
  if (field.space != FunctionSpace::cell_average) {
    throw Error("field acquisition requires a cell_average field");
  }
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  auto *data = field.bytes.data() + field.data_offset;
  return FieldView<T>(data, epoch_, detail::field_epoch_generation(epoch_),
                      interior_extent_, field.ghost_width, field.components,
                      field.x_stride, field.y_stride, field.z_stride);
}

template <class T>
FieldView<const T> FieldStorage::view(FieldId id) const {
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  const Entry &field = entry(id);
  if (field.space != FunctionSpace::cell_average) {
    throw Error("field acquisition requires a cell_average field");
  }
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  const auto *data = field.bytes.data() + field.data_offset;
  return FieldView<const T>(
      data, epoch_, detail::field_epoch_generation(epoch_), interior_extent_,
      field.ghost_width, field.components, field.x_stride, field.y_stride,
      field.z_stride);
}

template <class T>
FieldView<const T> FieldStorage::acquire_read(const FieldAccessPlan &plan,
                                              PhaseId phase, ActorId actor,
                                              FieldId id) const {
  plan.require_read(phase, actor, id, field_count());
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  const Entry &field = entry(id);
  if (field.space != FunctionSpace::cell_average) {
    throw Error("field acquisition requires a cell_average field");
  }
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  const auto *data = field.bytes.data() + field.data_offset;
  return FieldView<const T>(
      data, epoch_, detail::field_epoch_generation(epoch_), interior_extent_,
      field.ghost_width, field.components, field.x_stride, field.y_stride,
      field.z_stride);
}

template <class T>
FieldView<T> FieldStorage::acquire_write(const FieldAccessPlan &plan,
                                         PhaseId phase, ActorId actor,
                                         FieldId id) {
  plan.require_write(phase, actor, id, field_count());
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  Entry &field = entry(id);
  if (field.space != FunctionSpace::cell_average) {
    throw Error("field acquisition requires a cell_average field");
  }
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  auto *data = field.bytes.data() + field.data_offset;
  return FieldView<T>(data, epoch_, detail::field_epoch_generation(epoch_),
                      interior_extent_, field.ghost_width, field.components,
                      field.x_stride, field.y_stride, field.z_stride);
}

template <class T>
FaceFieldView<const T> FieldStorage::acquire_face_read(
    const FieldAccessPlan &plan, PhaseId phase, ActorId actor,
    FieldId id) const {
  plan.require_read(phase, actor, id, field_count());
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  const Entry &field = entry(id);
  if (field.space != FunctionSpace::face_value) {
    throw Error("field acquisition requires a face_value field");
  }
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  const auto *data = field.bytes.data() + field.data_offset;
  return FaceFieldView<const T>(
      data, epoch_, detail::field_epoch_generation(epoch_), field.face_count,
      field.components);
}

template <class T>
FaceFieldView<T> FieldStorage::acquire_face_write(
    const FieldAccessPlan &plan, PhaseId phase, ActorId actor, FieldId id) {
  plan.require_write(phase, actor, id, field_count());
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  Entry &field = entry(id);
  if (field.space != FunctionSpace::face_value) {
    throw Error("field acquisition requires a face_value field");
  }
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  auto *data = field.bytes.data() + field.data_offset;
  return FaceFieldView<T>(data, epoch_, detail::field_epoch_generation(epoch_),
                          field.face_count, field.components);
}

}  // namespace hundun::runtime
