// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_ibm.hpp"

#include "field_view_interval_detail.hpp"

#include <cmath>

namespace hundun::v04::detail {

inline bool same_force_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

inline bool valid_force_field(ConstFieldView field, Int3 cells,
                              std::uint8_t components,
                              std::uint8_t ghosts) noexcept {
  FieldStorageInterval interval{};
  return field.base != nullptr && same_force_shape(field.interior, cells) &&
         field.components == components && field.ghosts.x >= ghosts &&
         field.ghosts.y >= ghosts && field.ghosts.z >= ghosts &&
         field.revision != 0U && field.storage_identity != 0U &&
         field.revision_domain != 0U &&
         field_storage_interval(field, interval);
}

inline bool finite_force_vector(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

inline Real3 force_cross(Real3 left, Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

inline Real3 force_add(Real3 left, Real3 right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

inline Real3 force_scale(double scale, Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

}  // namespace hundun::v04::detail
