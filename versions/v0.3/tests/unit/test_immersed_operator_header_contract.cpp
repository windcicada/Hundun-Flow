// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_immersed_operator.hpp"

#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"

#include <type_traits>

namespace {

using hundun::finite_volume::FaceMassFlux;
using hundun::finite_volume::FiniteVolumeQuantity;
using hundun::finite_volume::ImmersedOperatorAdapter;
using hundun::finite_volume::ImmersedOperatorReport;
using hundun::finite_volume::ImmersedReconstruction;
using hundun::immersed::GhostStencilPlan;
using hundun::immersed::ImmersedDomain;
using hundun::immersed::LocalFlowPatternTransform;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::runtime::FaceFieldView;
using hundun::runtime::FieldView;

using Create = ImmersedOperatorAdapter (*)(const MeshTopology &,
                                           const MeshGeometry &,
                                           const ImmersedDomain &,
                                           const GhostStencilPlan &,
                                           const LocalFlowPatternTransform &,
                                           const ImmersedReconstruction &);
using MomentumUnconstrained = void (ImmersedOperatorAdapter::*)(
    const FaceMassFlux &, const FaceFieldView<const double> &,
    const FieldView<const double> &, const FieldView<const double> &,
    const FieldView<const double> &, const FaceFieldView<const double> &,
    const FieldView<double> &) const;
using Transport = void (ImmersedOperatorAdapter::*)(
    FiniteVolumeQuantity, const FaceMassFlux &,
    const FaceFieldView<const double> &, const FieldView<const double> &,
    const FieldView<const double> &, const FaceFieldView<const double> &,
    const FieldView<double> &) const;

static_assert(std::is_final_v<hundun::finite_volume::ImmersedResidualParts>);
static_assert(std::is_final_v<ImmersedOperatorReport>);
static_assert(std::is_final_v<ImmersedOperatorAdapter>);
static_assert(!std::is_copy_constructible_v<ImmersedOperatorAdapter>);
static_assert(!std::is_copy_assignable_v<ImmersedOperatorAdapter>);
static_assert(std::is_nothrow_move_constructible_v<ImmersedOperatorAdapter>);
static_assert(!std::is_move_assignable_v<ImmersedOperatorAdapter>);
static_assert(std::is_nothrow_destructible_v<ImmersedOperatorAdapter>);
static_assert(
    std::is_same_v<decltype(&ImmersedOperatorAdapter::create), Create>);
static_assert(
    std::is_same_v<decltype(&ImmersedOperatorAdapter::accumulate_momentum),
                   MomentumUnconstrained>);
static_assert(
    std::is_same_v<decltype(&ImmersedOperatorAdapter::accumulate_transport),
                   Transport>);
static_assert(
    std::is_same_v<
        decltype(std::declval<const ImmersedOperatorAdapter &>().report()),
        ImmersedOperatorReport>);

} // namespace

int main() {
  const auto create = &ImmersedOperatorAdapter::create;
  return create == nullptr ? 1 : 0;
}
