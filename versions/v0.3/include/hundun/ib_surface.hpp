// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <utility>

namespace hundun::immersed {

using TriangleId = std::uint64_t;
using ImmersedLinkId = std::uint64_t;

struct SurfaceTriangle final {
  TriangleId id{};
  std::array<runtime::Real3, 3> vertices_m{};
  runtime::Real3 geometric_outward_normal{};
  double area_m2{};
};

namespace detail {
struct SurfaceStorage;
}

namespace test {
class ImmersedTestAccess;
}

class ImmersedSurface final {
public:
  static ImmersedSurface load_collective(const std::filesystem::path &path,
                                         double length_scale_to_m,
                                         const runtime::MpiContext &mpi,
                                         int root);

  std::size_t vertex_count() const noexcept;
  std::size_t triangle_count() const noexcept;
  const SurfaceTriangle &triangle(TriangleId id) const;
  runtime::Real3 bounding_box_min_m() const noexcept;
  runtime::Real3 bounding_box_max_m() const noexcept;
  double closed_volume_m3() const noexcept;
  std::uint64_t fingerprint() const noexcept;

private:
  friend class SurfaceQuery;
  friend class test::ImmersedTestAccess;

  explicit ImmersedSurface(
      std::shared_ptr<const detail::SurfaceStorage> storage)
      : storage_(std::move(storage)) {}

  static ImmersedSurface load_collective_impl(const std::filesystem::path &path,
                                              double length_scale_to_m,
                                              const runtime::MpiContext &mpi,
                                              int root,
                                              std::size_t maximum_chunk_bytes);

  std::shared_ptr<const detail::SurfaceStorage> storage_;
};

} // namespace hundun::immersed
