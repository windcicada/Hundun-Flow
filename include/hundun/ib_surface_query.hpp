// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/ib_surface.hpp"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace hundun::config {
enum class ImmersedFluidSide : std::uint8_t;
}

namespace hundun::immersed {

struct ClosestPointResult final {
  TriangleId triangle{};
  runtime::Real3 point_m{};
  runtime::Real3 geometric_outward_normal{};
  double squared_distance_m2{};
};

struct SegmentIntersection final {
  TriangleId triangle{};
  runtime::Real3 point_m{};
  double segment_fraction{};
};

enum class CellRegion : std::uint8_t { fluid, solid };

namespace detail {
struct SurfaceQueryStorage;
}

class SurfaceQuery final {
public:
  static SurfaceQuery create(const ImmersedSurface &surface);

  ClosestPointResult closest_point(runtime::Real3 point_m) const;
  std::vector<SegmentIntersection>
  segment_intersections(runtime::Real3 a_m, runtime::Real3 b_m) const;
  CellRegion classify(runtime::Real3 point_m,
                      config::ImmersedFluidSide fluid_side) const;
  std::uint64_t fingerprint() const noexcept;

private:
  friend class test::ImmersedTestAccess;

  explicit SurfaceQuery(
      std::shared_ptr<const detail::SurfaceQueryStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::SurfaceQueryStorage> storage_;
};

} // namespace hundun::immersed
