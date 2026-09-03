// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_execution.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04::detail {

struct FieldStorageInterval {
  std::uintptr_t begin{};
  std::uintptr_t end{};
};

inline bool interval_checked_add(std::size_t left, std::size_t right,
                                 std::size_t& out) noexcept {
  if (left > std::numeric_limits<std::size_t>::max() - right) {
    return false;
  }
  out = left + right;
  return true;
}

inline bool interval_checked_multiply(std::size_t left, std::size_t right,
                                      std::size_t& out) noexcept {
  if (left != 0U &&
      right > std::numeric_limits<std::size_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

template <class T>
bool field_storage_interval(BasicFieldView<T> view,
                            FieldStorageInterval& out) noexcept {
  if (view.base == nullptr || view.interior.x <= 0 ||
      view.interior.y <= 0 || view.interior.z <= 0 || view.ghosts.x < 0 ||
      view.ghosts.y < 0 || view.ghosts.z < 0 || view.components == 0U ||
      view.stride_y == 0U || view.stride_z == 0U ||
      view.component_stride == 0U) {
    return false;
  }

  const std::size_t gx = static_cast<std::size_t>(view.ghosts.x);
  const std::size_t gy = static_cast<std::size_t>(view.ghosts.y);
  const std::size_t gz = static_cast<std::size_t>(view.ghosts.z);
  const std::size_t nx = static_cast<std::size_t>(view.interior.x);
  const std::size_t ny = static_cast<std::size_t>(view.interior.y);
  const std::size_t nz = static_cast<std::size_t>(view.interior.z);

  std::size_t padded_x = 0U;
  std::size_t padded_y = 0U;
  std::size_t padded_z = 0U;
  std::size_t twice_ghost = 0U;
  if (!interval_checked_multiply(gx, 2U, twice_ghost) ||
      !interval_checked_add(nx, twice_ghost, padded_x) ||
      !interval_checked_multiply(gy, 2U, twice_ghost) ||
      !interval_checked_add(ny, twice_ghost, padded_y) ||
      !interval_checked_multiply(gz, 2U, twice_ghost) ||
      !interval_checked_add(nz, twice_ghost, padded_z)) {
    return false;
  }

  std::size_t minimum_stride_z = 0U;
  std::size_t minimum_component_stride = 0U;
  std::size_t storage_span = 0U;
  if (view.stride_y < padded_x ||
      !interval_checked_multiply(view.stride_y, padded_y,
                                 minimum_stride_z) ||
      view.stride_z < minimum_stride_z ||
      !interval_checked_multiply(view.stride_z, padded_z,
                                 minimum_component_stride) ||
      view.component_stride < minimum_component_stride ||
      !interval_checked_multiply(
          view.component_stride, static_cast<std::size_t>(view.components),
          storage_span)) {
    return false;
  }

  std::size_t prefix_y = 0U;
  std::size_t prefix_z = 0U;
  std::size_t prefix = 0U;
  if (!interval_checked_multiply(gy, view.stride_y, prefix_y) ||
      !interval_checked_multiply(gz, view.stride_z, prefix_z) ||
      !interval_checked_add(gx, prefix_y, prefix) ||
      !interval_checked_add(prefix, prefix_z, prefix)) {
    return false;
  }

  std::size_t last_x = 0U;
  std::size_t last_y_index = 0U;
  std::size_t last_z_index = 0U;
  std::size_t last_y = 0U;
  std::size_t last_z = 0U;
  std::size_t last_component = 0U;
  std::size_t maximum_offset = 0U;
  std::size_t end_offset = 0U;
  if (!interval_checked_add(nx, gx, last_x) ||
      !interval_checked_add(ny, gy, last_y_index) ||
      !interval_checked_add(nz, gz, last_z_index) || last_x == 0U ||
      last_y_index == 0U || last_z_index == 0U ||
      !interval_checked_multiply(last_y_index - 1U, view.stride_y, last_y) ||
      !interval_checked_multiply(last_z_index - 1U, view.stride_z, last_z) ||
      !interval_checked_multiply(static_cast<std::size_t>(view.components - 1U),
                                 view.component_stride, last_component) ||
      !interval_checked_add(last_x - 1U, last_y, maximum_offset) ||
      !interval_checked_add(maximum_offset, last_z, maximum_offset) ||
      !interval_checked_add(maximum_offset, last_component, maximum_offset) ||
      !interval_checked_add(maximum_offset, 1U, end_offset)) {
    return false;
  }

  const auto ptrdiff_max =
      static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max());
  if (view.stride_y > ptrdiff_max || view.stride_z > ptrdiff_max ||
      view.component_stride > ptrdiff_max || prefix > ptrdiff_max ||
      maximum_offset > ptrdiff_max || storage_span == 0U) {
    return false;
  }

  std::size_t prefix_bytes = 0U;
  std::size_t end_bytes = 0U;
  if (!interval_checked_multiply(prefix, sizeof(double), prefix_bytes) ||
      !interval_checked_multiply(end_offset, sizeof(double), end_bytes)) {
    return false;
  }
  const auto base = reinterpret_cast<std::uintptr_t>(view.base);
  if (prefix_bytes > base ||
      end_bytes > std::numeric_limits<std::uintptr_t>::max() - base) {
    return false;
  }
  out = {base - prefix_bytes, base + end_bytes};
  return out.begin < out.end;
}

template <class T>
bool face_storage_interval(BasicFaceFieldView<T> view,
                           FieldStorageInterval& out) noexcept {
  if (view.base == nullptr || view.extents.x <= 0 || view.extents.y <= 0 ||
      view.extents.z <= 0 || view.stride_y == 0U || view.stride_z == 0U) {
    return false;
  }

  const std::size_t nx = static_cast<std::size_t>(view.extents.x);
  const std::size_t ny = static_cast<std::size_t>(view.extents.y);
  const std::size_t nz = static_cast<std::size_t>(view.extents.z);
  std::size_t minimum_stride_z = 0U;
  std::size_t storage_span = 0U;
  if (view.stride_y < nx ||
      !interval_checked_multiply(view.stride_y, ny, minimum_stride_z) ||
      view.stride_z < minimum_stride_z ||
      !interval_checked_multiply(view.stride_z, nz, storage_span)) {
    return false;
  }

  std::size_t last_y = 0U;
  std::size_t last_z = 0U;
  std::size_t maximum_offset = nx - 1U;
  std::size_t end_offset = 0U;
  if (!interval_checked_multiply(ny - 1U, view.stride_y, last_y) ||
      !interval_checked_multiply(nz - 1U, view.stride_z, last_z) ||
      !interval_checked_add(maximum_offset, last_y, maximum_offset) ||
      !interval_checked_add(maximum_offset, last_z, maximum_offset) ||
      !interval_checked_add(maximum_offset, 1U, end_offset)) {
    return false;
  }

  const auto ptrdiff_max =
      static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max());
  if (view.stride_y > ptrdiff_max || view.stride_z > ptrdiff_max ||
      maximum_offset > ptrdiff_max || storage_span == 0U) {
    return false;
  }

  std::size_t end_bytes = 0U;
  if (!interval_checked_multiply(end_offset, sizeof(double), end_bytes)) {
    return false;
  }
  const auto base = reinterpret_cast<std::uintptr_t>(view.base);
  if (end_bytes > std::numeric_limits<std::uintptr_t>::max() - base) {
    return false;
  }
  out = {base, base + end_bytes};
  return out.begin < out.end;
}

inline bool storage_intervals_overlap(FieldStorageInterval left,
                                      FieldStorageInterval right) noexcept {
  return left.begin < right.end && right.begin < left.end;
}

template <class Left, class Right>
bool field_views_overlap(BasicFieldView<Left> left,
                         BasicFieldView<Right> right) noexcept {
  FieldStorageInterval left_interval{};
  FieldStorageInterval right_interval{};
  if (!field_storage_interval(left, left_interval) ||
      !field_storage_interval(right, right_interval)) {
    return true;
  }
  return storage_intervals_overlap(left_interval, right_interval);
}

template <class CellT, class FaceT>
bool cell_face_views_overlap(BasicFieldView<CellT> cell,
                             BasicFaceFieldView<FaceT> face) noexcept {
  FieldStorageInterval cell_interval{};
  FieldStorageInterval face_interval{};
  if (!field_storage_interval(cell, cell_interval) ||
      !face_storage_interval(face, face_interval)) {
    return true;
  }
  return storage_intervals_overlap(cell_interval, face_interval);
}

template <class Left, class Right>
bool face_views_overlap(BasicFaceFieldView<Left> left,
                        BasicFaceFieldView<Right> right) noexcept {
  FieldStorageInterval left_interval{};
  FieldStorageInterval right_interval{};
  if (!face_storage_interval(left, left_interval) ||
      !face_storage_interval(right, right_interval)) {
    return true;
  }
  return storage_intervals_overlap(left_interval, right_interval);
}

template <class T>
bool field_view_overlaps_storage(BasicFieldView<T> view,
                                 const double* storage,
                                 std::size_t storage_doubles) noexcept {
  FieldStorageInterval interval{};
  if (!field_storage_interval(view, interval) || storage == nullptr ||
      storage_doubles == 0U ||
      storage_doubles >
          std::numeric_limits<std::uintptr_t>::max() / sizeof(double)) {
    return true;
  }
  const auto begin = reinterpret_cast<std::uintptr_t>(storage);
  const auto bytes = static_cast<std::uintptr_t>(storage_doubles) *
                     static_cast<std::uintptr_t>(sizeof(double));
  if (bytes > std::numeric_limits<std::uintptr_t>::max() - begin) {
    return true;
  }
  const auto end = begin + bytes;
  return storage_intervals_overlap(interval, {begin, end});
}

}  // namespace hundun::v04::detail
