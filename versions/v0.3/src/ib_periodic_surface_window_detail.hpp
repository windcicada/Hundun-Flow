// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_surface.hpp"
#include "hundun/rt_types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>

namespace hundun::immersed::detail {

struct PeriodicCellImage final {
  runtime::Int3 canonical{};
  runtime::Int3 image{};
  runtime::Real3 shift_m{};
};

class PeriodicCellMapper final {
public:
  PeriodicCellMapper(runtime::Int3 extent, std::array<bool, 3> periodic,
                     runtime::Real3 length_m) noexcept
      : extent_(extent), periodic_(periodic), length_m_(length_m) {}

  std::optional<PeriodicCellImage>
  image(runtime::Int3 image_cell) const noexcept {
    PeriodicCellImage result;
    result.image = image_cell;
    result.canonical = image_cell;
    for (std::size_t axis = 0U; axis < periodic_.size(); ++axis) {
      const int count = int_coordinate(extent_, axis);
      const int value = int_coordinate(image_cell, axis);
      if (count <= 0) {
        return std::nullopt;
      }
      if (!periodic_[axis]) {
        if (value < 0 || value >= count) {
          return std::nullopt;
        }
        continue;
      }
      int period = value / count;
      int canonical = value % count;
      if (canonical < 0) {
        canonical += count;
        --period;
      }
      set_int_coordinate(result.canonical, axis, canonical);
      set_real_coordinate(result.shift_m, axis,
                          static_cast<double>(period) *
                              real_coordinate(length_m_, axis));
    }
    return result;
  }

  PeriodicCellImage nearest_image(runtime::Int3 canonical,
                                  runtime::Int3 anchor) const noexcept {
    runtime::Int3 image_cell = canonical;
    for (std::size_t axis = 0U; axis < periodic_.size(); ++axis) {
      if (!periodic_[axis]) {
        continue;
      }
      const int count = int_coordinate(extent_, axis);
      int value = int_coordinate(image_cell, axis);
      const int anchor_value = int_coordinate(anchor, axis);
      if (value - anchor_value > count / 2) {
        value -= count;
      } else if (value - anchor_value < -count / 2) {
        value += count;
      }
      set_int_coordinate(image_cell, axis, value);
    }
    const auto result = image(image_cell);
    return result.has_value() ? *result : PeriodicCellImage{};
  }

  PeriodicCellImage nearest_image_to_box(runtime::Int3 canonical,
                                         runtime::Box3 box) const noexcept {
    runtime::Int3 image_cell = canonical;
    for (std::size_t axis = 0U; axis < periodic_.size(); ++axis) {
      if (!periodic_[axis]) {
        continue;
      }
      const int count = int_coordinate(extent_, axis);
      const int value = int_coordinate(canonical, axis);
      const int begin = int_coordinate(box.begin, axis);
      const int end = int_coordinate(box.end, axis);
      const auto distance = [begin, end](int candidate) noexcept {
        if (candidate < begin)
          return begin - candidate;
        if (candidate >= end)
          return candidate - (end - 1);
        return 0;
      };
      int selected = value;
      int selected_distance = distance(value);
      for (const int candidate : {value - count, value + count}) {
        const int candidate_distance = distance(candidate);
        if (candidate_distance < selected_distance) {
          selected = candidate;
          selected_distance = candidate_distance;
        }
      }
      set_int_coordinate(image_cell, axis, selected);
    }
    const auto result = image(image_cell);
    return result.has_value() ? *result : PeriodicCellImage{};
  }

private:
  static int int_coordinate(runtime::Int3 value, std::size_t axis) noexcept {
    return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
  }

  static void set_int_coordinate(runtime::Int3 &value, std::size_t axis,
                                 int coordinate) noexcept {
    if (axis == 0U)
      value.x = coordinate;
    else if (axis == 1U)
      value.y = coordinate;
    else
      value.z = coordinate;
  }

  static double real_coordinate(runtime::Real3 value,
                                std::size_t axis) noexcept {
    return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
  }

  static void set_real_coordinate(runtime::Real3 &value, std::size_t axis,
                                  double coordinate) noexcept {
    if (axis == 0U)
      value.x = coordinate;
    else if (axis == 1U)
      value.y = coordinate;
    else
      value.z = coordinate;
  }

  runtime::Int3 extent_{};
  std::array<bool, 3> periodic_{};
  runtime::Real3 length_m_{};
};

class PeriodicSurfaceWindow final {
public:
  PeriodicSurfaceWindow(runtime::Real3 minimum, runtime::Real3 maximum,
                        std::array<bool, 3> periodic) noexcept
      : minimum_(minimum), maximum_(maximum), periodic_(periodic) {
    const double scale = std::max(
        {1.0, maximum_.x - minimum_.x, maximum_.y - minimum_.y,
         maximum_.z - minimum_.z});
    plane_tolerance_ =
        16.0 * static_cast<double>(std::numeric_limits<float>::epsilon()) *
        scale;
  }

  bool bounding_box_is_admissible(runtime::Real3 surface_minimum,
                                  runtime::Real3 surface_maximum,
                                  double maximum_cell_scale) const noexcept {
    if (!(maximum_cell_scale > 0.0) ||
        !std::isfinite(maximum_cell_scale)) {
      return false;
    }
    for (std::size_t axis = 0U; axis < periodic_.size(); ++axis) {
      const double domain_minimum = coordinate(minimum_, axis);
      const double domain_maximum = coordinate(maximum_, axis);
      const double surface_lower = coordinate(surface_minimum, axis);
      const double surface_upper = coordinate(surface_maximum, axis);
      const bool separated =
          surface_lower - domain_minimum >= 2.0 * maximum_cell_scale &&
          domain_maximum - surface_upper >= 2.0 * maximum_cell_scale;
      const bool spans_period =
          periodic_[axis] &&
          surface_lower < domain_minimum - plane_tolerance_ &&
          surface_upper > domain_maximum + plane_tolerance_;
      if (!separated && !spans_period) {
        return false;
      }
    }
    return true;
  }

  bool triangle_is_split_at_periodic_planes(
      const SurfaceTriangle &triangle) const noexcept {
    for (std::size_t axis = 0U; axis < periodic_.size(); ++axis) {
      if (!periodic_[axis]) {
        continue;
      }
      double lower = coordinate(triangle.vertices_m.front(), axis);
      double upper = lower;
      for (const runtime::Real3 vertex : triangle.vertices_m) {
        lower = std::min(lower, coordinate(vertex, axis));
        upper = std::max(upper, coordinate(vertex, axis));
      }
      for (const double plane :
           {coordinate(minimum_, axis), coordinate(maximum_, axis)}) {
        if (lower < plane - plane_tolerance_ &&
            upper > plane + plane_tolerance_) {
          return false;
        }
      }
    }
    return true;
  }

  bool triangle_is_active(const SurfaceTriangle &triangle) const noexcept {
    runtime::Real3 centroid{};
    for (const runtime::Real3 vertex : triangle.vertices_m) {
      centroid.x += vertex.x / 3.0;
      centroid.y += vertex.y / 3.0;
      centroid.z += vertex.z / 3.0;
    }
    for (std::size_t axis = 0U; axis < periodic_.size(); ++axis) {
      if (!periodic_[axis]) {
        continue;
      }
      const double value = coordinate(centroid, axis);
      if (value < coordinate(minimum_, axis) ||
          value >= coordinate(maximum_, axis)) {
        return false;
      }
    }
    return true;
  }

private:
  static double coordinate(runtime::Real3 value,
                           std::size_t axis) noexcept {
    return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
  }

  runtime::Real3 minimum_{};
  runtime::Real3 maximum_{};
  std::array<bool, 3> periodic_{};
  double plane_tolerance_{};
};

} // namespace hundun::immersed::detail
