// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_domain.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::config::ImmersedFluidSide;
using hundun::immersed::ActiveBoundaryLayout;
using hundun::immersed::ActiveCellLayout;
using hundun::immersed::CellRegion;
using hundun::immersed::ImmersedDomain;
using hundun::immersed::ImmersedLink;
using hundun::immersed::ImmersedLinkId;
using hundun::immersed::ImmersedSurface;
using hundun::immersed::SurfaceQuery;
using hundun::immersed::TriangleId;
using hundun::mesh::GlobalCellId;
using hundun::mesh::GlobalFaceId;
using hundun::mesh::LocalCellId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;

using CreateDomain = ImmersedDomain (*)(const ImmersedSurface &,
                                        const SurfaceQuery &, ImmersedFluidSide,
                                        const MeshTopology &,
                                        const MeshGeometry &,
                                        const BoundaryRegistry &,
                                        const MpiContext &);
using RegionQuery = CellRegion (ImmersedDomain::*)(LocalCellId) const;

static_assert(std::is_same_v<ImmersedLinkId, std::uint64_t>);
static_assert(std::is_final_v<ImmersedLink>);
static_assert(std::is_same_v<decltype(ImmersedLink::id), ImmersedLinkId>);
static_assert(std::is_same_v<decltype(ImmersedLink::fluid_cell), GlobalCellId>);
static_assert(std::is_same_v<decltype(ImmersedLink::solid_cell), GlobalCellId>);
static_assert(std::is_same_v<decltype(ImmersedLink::triangle), TriangleId>);
static_assert(std::is_same_v<decltype(ImmersedLink::wall_intercept_m), Real3>);
static_assert(
    std::is_same_v<decltype(ImmersedLink::solid_to_fluid_normal), Real3>);
static_assert(
    std::is_same_v<decltype(ImmersedLink::fluid_to_wall_fraction), double>);
static_assert(std::is_same_v<CreateDomain, decltype(&ImmersedDomain::create)>);
static_assert(std::is_same_v<RegionQuery, decltype(&ImmersedDomain::region)>);

static_assert(std::is_same_v<decltype(std::declval<const ActiveCellLayout &>()
                                          .owned_active_count()),
                             std::size_t>);
static_assert(std::is_same_v<decltype(std::declval<const ActiveCellLayout &>()
                                          .local_active_count()),
                             std::size_t>);
static_assert(std::is_same_v<decltype(std::declval<const ActiveCellLayout &>()
                                          .active(LocalCellId{})),
                             bool>);
static_assert(std::is_same_v<decltype(std::declval<const ActiveCellLayout &>()
                                          .active_index(LocalCellId{})),
                             std::optional<std::size_t>>);
static_assert(std::is_same_v<decltype(std::declval<const ActiveCellLayout &>()
                                          .ordered_global_ids()),
                             const std::vector<GlobalCellId> &>);
static_assert(std::is_same_v<
              decltype(std::declval<const ActiveCellLayout &>().fingerprint()),
              std::uint64_t>);
static_assert(
    noexcept(std::declval<const ActiveCellLayout &>().owned_active_count()));
static_assert(
    noexcept(std::declval<const ActiveCellLayout &>().local_active_count()));
static_assert(
    noexcept(std::declval<const ActiveCellLayout &>().ordered_global_ids()));
static_assert(noexcept(std::declval<const ActiveCellLayout &>().fingerprint()));

static_assert(
    std::is_same_v<decltype(std::declval<const ActiveBoundaryLayout &>()
                                .patch_faces(std::uint32_t{})),
                   const std::vector<GlobalFaceId> &>);
static_assert(
    std::is_same_v<
        decltype(std::declval<const ActiveBoundaryLayout &>().open_domain()),
        bool>);
static_assert(
    std::is_same_v<decltype(std::declval<const ActiveBoundaryLayout &>()
                                .has_pressure_reference()),
                   bool>);
static_assert(
    std::is_same_v<
        decltype(std::declval<const ActiveBoundaryLayout &>().fingerprint()),
        std::uint64_t>);
static_assert(
    noexcept(std::declval<const ActiveBoundaryLayout &>().open_domain()));
static_assert(noexcept(
    std::declval<const ActiveBoundaryLayout &>().has_pressure_reference()));
static_assert(
    noexcept(std::declval<const ActiveBoundaryLayout &>().fingerprint()));

static_assert(
    std::is_same_v<decltype(std::declval<const ImmersedDomain &>().links()),
                   const std::vector<ImmersedLink> &>);
static_assert(std::is_same_v<decltype(std::declval<const ImmersedDomain &>()
                                          .region(LocalCellId{})),
                             CellRegion>);
static_assert(std::is_same_v<
              decltype(std::declval<const ImmersedDomain &>().active_cells()),
              const ActiveCellLayout &>);
static_assert(std::is_same_v<decltype(std::declval<const ImmersedDomain &>()
                                          .active_boundaries()),
                             const ActiveBoundaryLayout &>);
static_assert(std::is_same_v<decltype(std::declval<const ImmersedDomain &>()
                                          .classification_fingerprint()),
                             std::uint64_t>);
static_assert(std::is_same_v<decltype(std::declval<const ImmersedDomain &>()
                                          .surface_coverage_fingerprint()),
                             std::uint64_t>);
static_assert(noexcept(std::declval<const ImmersedDomain &>().links()));
static_assert(noexcept(std::declval<const ImmersedDomain &>().active_cells()));
static_assert(
    noexcept(std::declval<const ImmersedDomain &>().active_boundaries()));
static_assert(noexcept(
    std::declval<const ImmersedDomain &>().classification_fingerprint()));
static_assert(noexcept(
    std::declval<const ImmersedDomain &>().surface_coverage_fingerprint()));
static_assert(std::is_same_v<
              decltype(std::declval<const ImmersedDomain &>()
                           .performance_counters()),
              hundun::diagnostics::Stage3PerformanceCounters>);
static_assert(noexcept(
    std::declval<const ImmersedDomain &>().performance_counters()));

static_assert(std::is_final_v<ActiveCellLayout>);
static_assert(std::is_final_v<ActiveBoundaryLayout>);
static_assert(std::is_final_v<ImmersedDomain>);
static_assert(!std::is_default_constructible_v<ActiveCellLayout>);
static_assert(!std::is_default_constructible_v<ActiveBoundaryLayout>);
static_assert(!std::is_default_constructible_v<ImmersedDomain>);
static_assert(std::is_copy_constructible_v<ActiveCellLayout>);
static_assert(std::is_copy_assignable_v<ActiveCellLayout>);
static_assert(std::is_nothrow_move_constructible_v<ActiveCellLayout>);
static_assert(std::is_nothrow_move_assignable_v<ActiveCellLayout>);
static_assert(std::is_copy_constructible_v<ActiveBoundaryLayout>);
static_assert(std::is_copy_assignable_v<ActiveBoundaryLayout>);
static_assert(std::is_nothrow_move_constructible_v<ActiveBoundaryLayout>);
static_assert(std::is_nothrow_move_assignable_v<ActiveBoundaryLayout>);
static_assert(std::is_copy_constructible_v<ImmersedDomain>);
static_assert(std::is_copy_assignable_v<ImmersedDomain>);
static_assert(std::is_nothrow_move_constructible_v<ImmersedDomain>);
static_assert(std::is_nothrow_move_assignable_v<ImmersedDomain>);

} // namespace

int main() { return 0; }
