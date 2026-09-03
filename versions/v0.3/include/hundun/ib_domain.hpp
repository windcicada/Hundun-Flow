// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_surface_query.hpp"
#include "hundun/mesh_topology.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace hundun::boundary {
class BoundaryRegistry;
}

namespace hundun::mesh {
class MeshGeometry;
}

namespace hundun::immersed {

struct ImmersedLink final {
  ImmersedLinkId id{};
  mesh::GlobalCellId fluid_cell{};
  mesh::GlobalCellId solid_cell{};
  TriangleId triangle{};
  runtime::Real3 wall_intercept_m{};
  runtime::Real3 solid_to_fluid_normal{};
  double fluid_to_wall_fraction{};
};

namespace detail {
struct ActiveCellLayoutStorage;
struct ActiveBoundaryLayoutStorage;
struct ImmersedDomainStorage;
} // namespace detail

class ImmersedDomain;

class ActiveCellLayout final {
public:
  std::size_t owned_active_count() const noexcept;
  std::size_t local_active_count() const noexcept;
  bool active(mesh::LocalCellId local_cell) const;
  std::optional<std::size_t> active_index(mesh::LocalCellId local_cell) const;
  const std::vector<mesh::GlobalCellId> &ordered_global_ids() const noexcept;
  std::uint64_t fingerprint() const noexcept;

private:
  friend class ImmersedDomain;

  explicit ActiveCellLayout(
      std::shared_ptr<const detail::ActiveCellLayoutStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::ActiveCellLayoutStorage> storage_;
};

class ActiveBoundaryLayout final {
public:
  const std::vector<mesh::GlobalFaceId> &
  patch_faces(std::uint32_t stable_patch_id) const;
  bool open_domain() const noexcept;
  bool has_pressure_reference() const noexcept;
  std::uint64_t fingerprint() const noexcept;

private:
  friend class ImmersedDomain;

  explicit ActiveBoundaryLayout(
      std::shared_ptr<const detail::ActiveBoundaryLayoutStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::ActiveBoundaryLayoutStorage> storage_;
};

class ImmersedDomain final {
public:
  static ImmersedDomain create(const ImmersedSurface &surface,
                               const SurfaceQuery &query,
                               config::ImmersedFluidSide fluid_side,
                               const mesh::MeshTopology &topology,
                               const mesh::MeshGeometry &geometry,
                               const boundary::BoundaryRegistry &boundaries,
                               const runtime::MpiContext &mpi);

  CellRegion region(mesh::LocalCellId local_cell) const;
  const std::vector<ImmersedLink> &links() const noexcept;
  const ActiveCellLayout &active_cells() const noexcept;
  const ActiveBoundaryLayout &active_boundaries() const noexcept;
  std::uint64_t classification_fingerprint() const noexcept;
  std::uint64_t surface_coverage_fingerprint() const noexcept;
  diagnostics::Stage3PerformanceCounters
  performance_counters() const noexcept;

private:
  explicit ImmersedDomain(
      std::shared_ptr<const detail::ImmersedDomainStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::ImmersedDomainStorage> storage_;
};

} // namespace hundun::immersed
