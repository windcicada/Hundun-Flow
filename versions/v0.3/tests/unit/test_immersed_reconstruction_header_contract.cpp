// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_immersed_reconstruction.hpp"

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"

#include <cstdint>
#include <type_traits>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::finite_volume::FaceMassFlux;
using hundun::finite_volume::FiniteVolumeQuantity;
using hundun::finite_volume::GradientScheme;
using hundun::finite_volume::ImmersedReconstruction;
using hundun::finite_volume::ReconstructionFieldBinding;
using hundun::immersed::GhostStencilPlan;
using hundun::immersed::ImmersedDomain;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::ActorId;
using hundun::runtime::FieldAccessPlan;
using hundun::runtime::FieldId;
using hundun::runtime::FieldStorage;
using hundun::runtime::HaloExchange;
using hundun::runtime::MpiContext;
using hundun::runtime::PhaseId;
using hundun::runtime::StructuredDecomposition;

using Create = ImmersedReconstruction (*)(
    const MeshTopology &, const MeshGeometry &, const BoundaryRegistry &,
    const ImmersedDomain &, const GhostStencilPlan &,
    const StructuredDecomposition &, const MpiContext &, HaloExchange &);
using ComputeGradient = void (ImmersedReconstruction::*)(
    GradientScheme, FiniteVolumeQuantity, const ReconstructionFieldBinding &,
    const ReconstructionFieldBinding &) const;
template <class T, class = void>
struct HasConstrainedPressureGradient : std::false_type {};
template <class T>
struct HasConstrainedPressureGradient<
    T,
    std::void_t<
        decltype(&T::compute_pressure_gradient_with_wall_normal_constraints)>>
    : std::true_type {};
using ReconstructTransport =
    void (ImmersedReconstruction::*)(FiniteVolumeQuantity, const FaceMassFlux &,
                                     const ReconstructionFieldBinding &,
                                     const ReconstructionFieldBinding &) const;
using ReconstructMomentum = void (ImmersedReconstruction::*)(
    const FaceMassFlux &, const ReconstructionFieldBinding &,
    const ReconstructionFieldBinding &) const;

static_assert(std::is_final_v<ReconstructionFieldBinding>);
static_assert(std::is_same_v<decltype(ReconstructionFieldBinding::storage),
                             FieldStorage &>);
static_assert(std::is_same_v<decltype(ReconstructionFieldBinding::access_plan),
                             const FieldAccessPlan &>);
static_assert(
    std::is_same_v<decltype(ReconstructionFieldBinding::phase), PhaseId>);
static_assert(
    std::is_same_v<decltype(ReconstructionFieldBinding::actor), ActorId>);
static_assert(
    std::is_same_v<decltype(ReconstructionFieldBinding::field), FieldId>);
static_assert(std::is_final_v<ImmersedReconstruction>);
static_assert(!std::is_copy_constructible_v<ImmersedReconstruction>);
static_assert(!std::is_copy_assignable_v<ImmersedReconstruction>);
static_assert(std::is_nothrow_move_constructible_v<ImmersedReconstruction>);
static_assert(!std::is_move_assignable_v<ImmersedReconstruction>);
static_assert(std::is_nothrow_destructible_v<ImmersedReconstruction>);
static_assert(
    std::is_same_v<decltype(&ImmersedReconstruction::create), Create>);
static_assert(
    std::is_same_v<decltype(&ImmersedReconstruction::compute_gradient),
                   ComputeGradient>);
static_assert(!HasConstrainedPressureGradient<ImmersedReconstruction>::value);
static_assert(std::is_same_v<
              decltype(&ImmersedReconstruction::reconstruct_transport_faces),
              ReconstructTransport>);
static_assert(std::is_same_v<
              decltype(&ImmersedReconstruction::reconstruct_momentum_faces),
              ReconstructMomentum>);
static_assert(noexcept(
    std::declval<const ImmersedReconstruction &>().dependency_fingerprint()));
static_assert(noexcept(FiniteVolumeQuantity::pressure()));

} // namespace

int main() {
  const auto pressure = FiniteVolumeQuantity::pressure();
  return pressure.kind == hundun::finite_volume::FiniteVolumeQuantityKind::
                              pressure &&
                 pressure.scalar_index == 0U
             ? 0
             : 1;
}
