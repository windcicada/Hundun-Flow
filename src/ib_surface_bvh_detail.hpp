// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_surface_query.hpp"

#include "ib_stl_reader_detail.hpp"

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::immersed::detail {

struct AxisAlignedBox final {
  runtime::Real3 minimum{};
  runtime::Real3 maximum{};
};

struct BvhNode final {
  AxisAlignedBox bounds{};
  std::uint32_t left{};
  std::uint32_t right{};
  std::uint64_t begin{};
  std::uint64_t end{};
  bool leaf{};
};

struct SurfaceQueryStorage final {
  std::shared_ptr<const SurfaceStorage> surface;
  std::vector<AxisAlignedBox> triangle_bounds;
  std::vector<TriangleId> triangle_order;
  std::vector<BvhNode> nodes;
  std::uint64_t fingerprint{};
  mutable std::atomic<std::uint64_t> closest_calls{};
  mutable std::atomic<std::uint64_t> segment_calls{};
};

std::shared_ptr<const SurfaceQueryStorage>
build_surface_query(std::shared_ptr<const SurfaceStorage> surface);

#if defined(HUNDUN_IMMERSED_ENABLE_TEST_ACCESS)
std::shared_ptr<const SurfaceQueryStorage>
build_surface_query_with_initial_order(
    std::shared_ptr<const SurfaceStorage> surface,
    std::vector<TriangleId> initial_order);
#endif

std::vector<TriangleId> bounded_candidates(const SurfaceQueryStorage &query,
                                           runtime::Real3 minimum_m,
                                           runtime::Real3 maximum_m);

bool require_consistent_parity(const std::array<bool, 3> &votes);

} // namespace hundun::immersed::detail
