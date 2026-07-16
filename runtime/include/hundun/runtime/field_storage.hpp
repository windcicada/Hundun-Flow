// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_view.hpp"
#include "hundun/runtime/types.hpp"

namespace hundun::runtime {

namespace detail {

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

class FieldStorage final {
 public:
  FieldStorage(const FieldRegistry &, Int3 interior_extent);
  Int3 interior_extent() const noexcept;

  template <class T>
  FieldView<T> view(FieldId id);

  template <class T>
  FieldView<const T> view(FieldId id) const;

 private:
  struct Entry {
    ScalarType scalar_type{};
    int ghost_width{};
    std::uint32_t components{};
    std::size_t x_stride{};
    std::size_t y_stride{};
    std::size_t z_stride{};
    std::size_t data_offset{};
    std::vector<std::byte> bytes;
  };

  Entry &entry(FieldId id);
  const Entry &entry(FieldId id) const;

  Int3 interior_extent_{};
  std::unique_ptr<std::vector<Entry>> entries_;
};

template <class T>
FieldView<T> FieldStorage::view(FieldId id) {
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  Entry &field = entry(id);
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  auto *data = field.bytes.data() + field.data_offset;
  return FieldView<T>(data, interior_extent_, field.ghost_width,
                      field.components, field.x_stride, field.y_stride,
                      field.z_stride);
}

template <class T>
FieldView<const T> FieldStorage::view(FieldId id) const {
  constexpr ScalarType requested_type = detail::FieldScalarType<T>::value;
  const Entry &field = entry(id);
  if (field.scalar_type != requested_type) {
    throw Error("field view scalar type does not match its descriptor");
  }
  const auto *data = field.bytes.data() + field.data_offset;
  return FieldView<const T>(data, interior_extent_, field.ghost_width,
                            field.components, field.x_stride, field.y_stride,
                            field.z_stride);
}

}  // namespace hundun::runtime
