// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_surface.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::immersed::detail {

struct IndexedTriangle final {
  TriangleId id{};
  std::array<std::uint64_t, 3> vertices{};
};

struct SurfaceStorage final {
  std::vector<runtime::Real3> vertices;
  std::vector<IndexedTriangle> indexed_triangles;
  std::vector<SurfaceTriangle> triangles;
  runtime::Real3 bounding_box_min{};
  runtime::Real3 bounding_box_max{};
  double reference_length_m{};
  double weld_tolerance_m{};
  double minimum_triangle_area_m2{};
  double intersection_coincidence_m{};
  double closed_volume_m3{};
  std::uint64_t fingerprint{};
};

runtime::Real3 add(runtime::Real3 a, runtime::Real3 b) noexcept;
runtime::Real3 subtract(runtime::Real3 a, runtime::Real3 b) noexcept;
runtime::Real3 multiply(runtime::Real3 value, double factor) noexcept;
double dot(runtime::Real3 a, runtime::Real3 b) noexcept;
runtime::Real3 cross(runtime::Real3 a, runtime::Real3 b) noexcept;
double squared_norm(runtime::Real3 value) noexcept;
double norm(runtime::Real3 value) noexcept;
bool finite(runtime::Real3 value) noexcept;

std::shared_ptr<const SurfaceStorage>
read_and_normalize_stl(const std::vector<std::uint8_t> &bytes,
                       double length_scale_to_m);

std::vector<std::uint8_t>
encode_normalized_surface(const SurfaceStorage &surface);

std::shared_ptr<const SurfaceStorage>
decode_normalized_surface(const std::vector<std::uint8_t> &bytes);

} // namespace hundun::immersed::detail
