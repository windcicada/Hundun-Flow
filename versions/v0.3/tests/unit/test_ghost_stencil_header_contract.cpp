// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_ghost_stencil_plan.hpp"

#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::AffineGhostConstraint;
using hundun::immersed::FluidExtrapolation;
using hundun::immersed::GhostStencilPlan;
using hundun::immersed::ImmersedDomain;
using hundun::immersed::ImmersedLinkId;
using hundun::immersed::ImmersedSurface;
using hundun::immersed::QuadraticReconstruction;
using hundun::immersed::SurfaceQuery;
using hundun::immersed::TriangleId;
using hundun::immersed::WallQuadraturePlan;
using hundun::immersed::WallQuadraturePoint;
using hundun::immersed::WeightedDonor;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

using GhostCreate = GhostStencilPlan (*)(
    const ImmersedSurface &, const SurfaceQuery &, const ImmersedDomain &,
    const MeshTopology &, const MeshGeometry &, const StructuredDecomposition &,
    const MpiContext &);
using WallCreate = WallQuadraturePlan (*)(
    const ImmersedSurface &, const SurfaceQuery &, const ImmersedDomain &,
    const MeshTopology &, const MeshGeometry &, const MpiContext &);

static_assert(std::is_final_v<AffineGhostConstraint>);
static_assert(
    std::is_same_v<decltype(AffineGhostConstraint::link), ImmersedLinkId>);
static_assert(std::is_same_v<decltype(AffineGhostConstraint::donors),
                             std::vector<WeightedDonor>>);
static_assert(
    std::is_same_v<decltype(AffineGhostConstraint::wall_value_weight), double>);
static_assert(std::is_same_v<
              decltype(AffineGhostConstraint::wall_normal_gradient_weight_m),
              double>);
static_assert(std::is_final_v<FluidExtrapolation>);
static_assert(
    std::is_same_v<decltype(FluidExtrapolation::link), ImmersedLinkId>);
static_assert(std::is_same_v<decltype(FluidExtrapolation::donors),
                             std::vector<WeightedDonor>>);
static_assert(std::is_same_v<GhostCreate, decltype(&GhostStencilPlan::create)>);
static_assert(
    std::is_same_v<decltype(std::declval<const GhostStencilPlan &>()
                                .velocity_constraint(ImmersedLinkId{}, 0U)),
                   const AffineGhostConstraint &>);
static_assert(
    std::is_same_v<decltype(std::declval<const GhostStencilPlan &>()
                                .zero_normal_constraint(ImmersedLinkId{})),
                   const AffineGhostConstraint &>);
static_assert(
    std::is_same_v<decltype(std::declval<const GhostStencilPlan &>()
                                .density_extrapolation(ImmersedLinkId{})),
                   const FluidExtrapolation &>);
static_assert(std::is_same_v<decltype(std::declval<const GhostStencilPlan &>()
                                          .reconstruction(ImmersedLinkId{})),
                             const QuadraticReconstruction &>);
static_assert(
    noexcept(std::declval<const GhostStencilPlan &>().maximum_halo_reach()));
static_assert(noexcept(std::declval<const GhostStencilPlan &>().fingerprint()));
static_assert(std::is_same_v<
              decltype(std::declval<const GhostStencilPlan &>()
                           .performance_counters()),
              hundun::diagnostics::Stage3PerformanceCounters>);
static_assert(std::is_copy_constructible_v<GhostStencilPlan>);
static_assert(!std::is_default_constructible_v<GhostStencilPlan>);

static_assert(std::is_final_v<WallQuadraturePoint>);
static_assert(
    std::is_same_v<decltype(WallQuadraturePoint::triangle), TriangleId>);
static_assert(
    std::is_same_v<decltype(WallQuadraturePoint::point_index), std::uint32_t>);
static_assert(std::is_same_v<decltype(WallQuadraturePoint::position_m), Real3>);
static_assert(std::is_same_v<
              decltype(WallQuadraturePoint::solid_to_fluid_normal), Real3>);
static_assert(std::is_same_v<decltype(WallQuadraturePoint::weight_m2), double>);
static_assert(std::is_same_v<decltype(WallQuadraturePoint::owner_rank), int>);
static_assert(std::is_same_v<decltype(WallQuadraturePoint::reconstruction),
                             QuadraticReconstruction>);
static_assert(
    std::is_same_v<WallCreate, decltype(&WallQuadraturePlan::create)>);
static_assert(std::is_same_v<decltype(std::declval<const WallQuadraturePlan &>()
                                          .local_points()),
                             const std::vector<WallQuadraturePoint> &>);
static_assert(
    noexcept(std::declval<const WallQuadraturePlan &>().local_points()));
static_assert(
    std::is_same_v<decltype(std::declval<const WallQuadraturePlan &>()
                                .maximum_halo_reach()),
                   std::uint32_t>);
static_assert(noexcept(std::declval<const WallQuadraturePlan &>()
                           .maximum_halo_reach()));
static_assert(
    noexcept(std::declval<const WallQuadraturePlan &>().fingerprint()));
static_assert(std::is_same_v<
              decltype(std::declval<const WallQuadraturePlan &>()
                           .performance_counters()),
              hundun::diagnostics::Stage3PerformanceCounters>);
static_assert(std::is_copy_constructible_v<WallQuadraturePlan>);
static_assert(!std::is_default_constructible_v<WallQuadraturePlan>);

} // namespace

int main() { return 0; }
