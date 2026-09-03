// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_surface.hpp"
#include "hundun/ib_surface_query.hpp"

#include "hundun/cfg_resolved_case_v3.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <type_traits>
#include <utility>
#include <vector>

static_assert(std::is_same_v<hundun::immersed::TriangleId, std::uint64_t>);
static_assert(std::is_same_v<hundun::immersed::ImmersedLinkId, std::uint64_t>);
static_assert(
    std::is_same_v<decltype(hundun::immersed::ImmersedSurface::load_collective(
                       std::declval<const std::filesystem::path &>(), 1.0,
                       std::declval<const hundun::runtime::MpiContext &>(), 0)),
                   hundun::immersed::ImmersedSurface>);
static_assert(std::is_same_v<
              decltype(std::declval<const hundun::immersed::ImmersedSurface &>()
                           .triangle(hundun::immersed::TriangleId{})),
              const hundun::immersed::SurfaceTriangle &>);
static_assert(std::is_same_v<
              decltype(std::declval<const hundun::immersed::ImmersedSurface &>()
                           .vertex_count()),
              std::size_t>);
static_assert(std::is_same_v<
              decltype(std::declval<const hundun::immersed::ImmersedSurface &>()
                           .triangle_count()),
              std::size_t>);
static_assert(noexcept(
    std::declval<const hundun::immersed::ImmersedSurface &>().vertex_count()));
static_assert(noexcept(std::declval<const hundun::immersed::ImmersedSurface &>()
                           .triangle_count()));
static_assert(noexcept(
    std::declval<const hundun::immersed::ImmersedSurface &>().fingerprint()));
static_assert(
    std::is_same_v<decltype(hundun::immersed::SurfaceTriangle::vertices_m),
                   std::array<hundun::runtime::Real3, 3>>);
static_assert(std::is_same_v<
              decltype(hundun::immersed::SurfaceQuery::create(
                  std::declval<const hundun::immersed::ImmersedSurface &>())),
              hundun::immersed::SurfaceQuery>);
static_assert(std::is_same_v<
              decltype(std::declval<const hundun::immersed::SurfaceQuery &>()
                           .closest_point({})),
              hundun::immersed::ClosestPointResult>);
static_assert(std::is_same_v<
              decltype(std::declval<const hundun::immersed::SurfaceQuery &>()
                           .segment_intersections({}, {})),
              std::vector<hundun::immersed::SegmentIntersection>>);
static_assert(
    std::is_same_v<
        decltype(std::declval<const hundun::immersed::SurfaceQuery &>()
                     .classify({}, hundun::config::ImmersedFluidSide::outside)),
        hundun::immersed::CellRegion>);
static_assert(sizeof(hundun::immersed::CellRegion) == sizeof(std::uint8_t));
static_assert(noexcept(
    std::declval<const hundun::immersed::SurfaceQuery &>().fingerprint()));
static_assert(std::is_final_v<hundun::immersed::ImmersedSurface>);
static_assert(std::is_final_v<hundun::immersed::SurfaceQuery>);

int main() { return 0; }
