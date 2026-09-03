// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <new>
#include <type_traits>
#include <utility>

#include "hundun/rt_field_view.hpp"
#include "hundun/rt_types.hpp"

namespace hundun::runtime {

// A non-owning, unchecked borrow for one lexically scoped numerical callback.
// Every coordinate and component must be valid; callers must not retain a
// copy after the callback returns.
template <class T>
class KernelCellView final {
 public:
  Int3 interior_extent() const noexcept { return interior_extent_; }
  int ghost_width() const noexcept { return ghost_width_; }
  std::uint32_t components() const noexcept { return components_; }

  T &operator()(int i, int j, int k, int component) const noexcept {
    const auto ii =
        static_cast<std::size_t>(static_cast<std::int64_t>(i) + ghost_width_);
    const auto jj =
        static_cast<std::size_t>(static_cast<std::int64_t>(j) + ghost_width_);
    const auto kk =
        static_cast<std::size_t>(static_cast<std::int64_t>(k) + ghost_width_);
    const auto linear =
        kk * z_stride_ + jj * y_stride_ + ii * x_stride_ +
        static_cast<std::size_t>(component);
    auto *value =
        std::launder(reinterpret_cast<T *>(data_ + linear * sizeof(Value)));
    return *value;
  }

 private:
  friend class detail::KernelFieldViewAccess;

  using Byte =
      std::conditional_t<std::is_const_v<T>, const std::byte, std::byte>;
  using Value = std::remove_const_t<T>;

  KernelCellView(Byte *data, Int3 interior_extent, int ghost_width,
                 std::uint32_t components, std::size_t x_stride,
                 std::size_t y_stride, std::size_t z_stride) noexcept
      : data_(data),
        interior_extent_(interior_extent),
        ghost_width_(ghost_width),
        components_(components),
        x_stride_(x_stride),
        y_stride_(y_stride),
        z_stride_(z_stride) {}

  Byte *data_;
  Int3 interior_extent_;
  int ghost_width_;
  std::uint32_t components_;
  std::size_t x_stride_;
  std::size_t y_stride_;
  std::size_t z_stride_;
};

// Face-layout analogue of KernelCellView, with the same unchecked lifetime
// and valid-index contract.
template <class T>
class KernelFaceView final {
 public:
  std::size_t face_count() const noexcept { return face_count_; }
  std::uint32_t components() const noexcept { return components_; }

  T &operator()(std::size_t face, int component) const noexcept {
    const auto linear =
        face * static_cast<std::size_t>(components_) +
        static_cast<std::size_t>(component);
    auto *value =
        std::launder(reinterpret_cast<T *>(data_ + linear * sizeof(Value)));
    return *value;
  }

 private:
  friend class detail::KernelFieldViewAccess;

  using Byte =
      std::conditional_t<std::is_const_v<T>, const std::byte, std::byte>;
  using Value = std::remove_const_t<T>;

  KernelFaceView(Byte *data, std::size_t face_count,
                 std::uint32_t components) noexcept
      : data_(data), face_count_(face_count), components_(components) {}

  Byte *data_;
  std::size_t face_count_;
  std::uint32_t components_;
};

template <class T, class Callback>
void with_kernel_cell_view(const FieldView<T> &checked,
                           Callback &&callback);

template <class T, class Callback>
void with_kernel_face_view(const FaceFieldView<T> &checked,
                           Callback &&callback);

namespace detail {

// This private bridge is the only type trusted by both checked and kernel
// views. Keeping its factories private prevents a second public construction
// path while allowing field_view.hpp to remain independent of the unchecked
// API.
class KernelFieldViewAccess final {
 private:
  KernelFieldViewAccess() = delete;

  template <class T>
  static KernelCellView<T> cell(const FieldView<T> &checked) {
    validate_field_epoch(checked.epoch_, checked.generation_);
    return KernelCellView<T>(
        checked.data_, checked.interior_extent_, checked.ghost_width_,
        checked.components_, checked.x_stride_, checked.y_stride_,
        checked.z_stride_);
  }

  template <class T>
  static KernelFaceView<T> face(const FaceFieldView<T> &checked) {
    validate_field_epoch(checked.epoch_, checked.generation_);
    return KernelFaceView<T>(checked.data_, checked.face_count_,
                             checked.components_);
  }

  template <class T, class Callback>
  friend void ::hundun::runtime::with_kernel_cell_view(
      const FieldView<T> &checked, Callback &&callback);

  template <class T, class Callback>
  friend void ::hundun::runtime::with_kernel_face_view(
      const FaceFieldView<T> &checked, Callback &&callback);
};

}  // namespace detail

// Each helper validates the checked borrow once, then invokes the callback
// once with a by-value kernel borrow. Callback results are deliberately
// discarded and exceptions propagate unchanged.
template <class T, class Callback>
void with_kernel_cell_view(const FieldView<T> &checked,
                           Callback &&callback) {
  static_cast<void>(std::invoke(
      std::forward<Callback>(callback),
      detail::KernelFieldViewAccess::cell(checked)));
}

template <class T, class Callback>
void with_kernel_face_view(const FaceFieldView<T> &checked,
                           Callback &&callback) {
  static_cast<void>(std::invoke(
      std::forward<Callback>(callback),
      detail::KernelFieldViewAccess::face(checked)));
}

}  // namespace hundun::runtime
