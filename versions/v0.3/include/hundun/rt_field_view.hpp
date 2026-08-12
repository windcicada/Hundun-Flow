// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <utility>

#include "hundun/rt_error.hpp"
#include "hundun/rt_types.hpp"

namespace hundun::runtime {

class FieldStorage;
template <class T>
class FieldView;
template <class T>
class FaceFieldView;

namespace detail {

class FieldEpochControl;
class KernelFieldViewAccess;

std::uint64_t field_epoch_generation(
    const std::shared_ptr<FieldEpochControl> &epoch) noexcept;
void validate_field_epoch(
    const std::shared_ptr<const FieldEpochControl> &epoch,
    std::uint64_t captured_generation);

}  // namespace detail

template <class T>
class FieldView final {
 public:
  Int3 interior_extent() const noexcept;
  int ghost_width() const noexcept;
  std::uint32_t components() const noexcept;
  T &operator()(int i, int j, int k, int component) const;

 private:
  friend class FieldStorage;
  friend class detail::KernelFieldViewAccess;

  using Byte =
      std::conditional_t<std::is_const_v<T>, const std::byte, std::byte>;
  using Value = std::remove_const_t<T>;

  FieldView(Byte *data,
            std::shared_ptr<const detail::FieldEpochControl> epoch,
            std::uint64_t generation, Int3 interior_extent, int ghost_width,
            std::uint32_t components, std::size_t x_stride,
            std::size_t y_stride, std::size_t z_stride) noexcept;

  Byte *data_{};
  std::shared_ptr<const detail::FieldEpochControl> epoch_;
  std::uint64_t generation_{};
  Int3 interior_extent_{};
  int ghost_width_{};
  std::uint32_t components_{};
  std::size_t x_stride_{};
  std::size_t y_stride_{};
  std::size_t z_stride_{};
};

template <class T>
class FaceFieldView final {
 public:
  std::size_t face_count() const noexcept;
  std::uint32_t components() const noexcept;
  T &operator()(std::size_t face, int component) const;

 private:
  friend class FieldStorage;
  friend class detail::KernelFieldViewAccess;

  using Byte =
      std::conditional_t<std::is_const_v<T>, const std::byte, std::byte>;
  using Value = std::remove_const_t<T>;

  FaceFieldView(Byte *data,
                std::shared_ptr<const detail::FieldEpochControl> epoch,
                std::uint64_t generation, std::size_t face_count,
                std::uint32_t components) noexcept;

  Byte *data_{};
  std::shared_ptr<const detail::FieldEpochControl> epoch_;
  std::uint64_t generation_{};
  std::size_t face_count_{};
  std::uint32_t components_{};
};

template <class T>
FieldView<T>::FieldView(
    Byte *data, std::shared_ptr<const detail::FieldEpochControl> epoch,
    std::uint64_t generation, Int3 interior_extent, int ghost_width,
    std::uint32_t components, std::size_t x_stride, std::size_t y_stride,
    std::size_t z_stride) noexcept
    : data_(data),
      epoch_(std::move(epoch)),
      generation_(generation),
      interior_extent_(interior_extent),
      ghost_width_(ghost_width),
      components_(components),
      x_stride_(x_stride),
      y_stride_(y_stride),
      z_stride_(z_stride) {}

template <class T>
Int3 FieldView<T>::interior_extent() const noexcept {
  return interior_extent_;
}

template <class T>
int FieldView<T>::ghost_width() const noexcept {
  return ghost_width_;
}

template <class T>
std::uint32_t FieldView<T>::components() const noexcept {
  return components_;
}

template <class T>
T &FieldView<T>::operator()(int i, int j, int k, int component) const {
  detail::validate_field_epoch(epoch_, generation_);

  const auto coordinate_in_bounds = [this](int coordinate, int extent) {
    const auto wide_coordinate = static_cast<std::int64_t>(coordinate);
    const auto lower = -static_cast<std::int64_t>(ghost_width_);
    const auto upper = static_cast<std::int64_t>(extent) + ghost_width_;
    return wide_coordinate >= lower && wide_coordinate < upper;
  };
  if (!coordinate_in_bounds(i, interior_extent_.x) ||
      !coordinate_in_bounds(j, interior_extent_.y) ||
      !coordinate_in_bounds(k, interior_extent_.z)) {
    throw Error("field view spatial index is outside the ghosted extent");
  }
  if (component < 0 || static_cast<std::uint64_t>(component) >= components_) {
    throw Error("field view component index is out of bounds");
  }

  const auto ii =
      static_cast<std::size_t>(static_cast<std::int64_t>(i) + ghost_width_);
  const auto jj =
      static_cast<std::size_t>(static_cast<std::int64_t>(j) + ghost_width_);
  const auto kk =
      static_cast<std::size_t>(static_cast<std::int64_t>(k) + ghost_width_);
  const auto component_index = static_cast<std::size_t>(component);
  const auto linear =
      kk * z_stride_ + jj * y_stride_ + ii * x_stride_ + component_index;
  auto *value =
      std::launder(reinterpret_cast<T *>(data_ + linear * sizeof(Value)));
  return *value;
}

template <class T>
FaceFieldView<T>::FaceFieldView(
    Byte *data, std::shared_ptr<const detail::FieldEpochControl> epoch,
    std::uint64_t generation, std::size_t face_count,
    std::uint32_t components) noexcept
    : data_(data),
      epoch_(std::move(epoch)),
      generation_(generation),
      face_count_(face_count),
      components_(components) {}

template <class T>
std::size_t FaceFieldView<T>::face_count() const noexcept {
  return face_count_;
}

template <class T>
std::uint32_t FaceFieldView<T>::components() const noexcept {
  return components_;
}

template <class T>
T &FaceFieldView<T>::operator()(std::size_t face, int component) const {
  detail::validate_field_epoch(epoch_, generation_);
  if (face >= face_count_) {
    throw Error("face field view face index is out of bounds");
  }
  if (component < 0 || static_cast<std::uint64_t>(component) >= components_) {
    throw Error("face field view component index is out of bounds");
  }

  const auto linear =
      face * static_cast<std::size_t>(components_) +
      static_cast<std::size_t>(component);
  auto *value =
      std::launder(reinterpret_cast<T *>(data_ + linear * sizeof(Value)));
  return *value;
}

}  // namespace hundun::runtime
