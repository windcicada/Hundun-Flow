// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/flow_immersed.hpp"

#include <type_traits>

namespace {

using hundun::flow::FixedStepImmersedFlow;
using hundun::flow::ImmersedFlowDensitySetup;
using hundun::flow::ImmersedFlowPhysics;
using hundun::flow::ImmersedFlowStepAttemptReport;
using hundun::flow::ImmersedLinearSolvePhase;

static_assert(std::is_final_v<ImmersedFlowDensitySetup>);
static_assert(std::is_final_v<ImmersedFlowPhysics>);
static_assert(std::is_final_v<ImmersedFlowStepAttemptReport>);
static_assert(std::is_final_v<FixedStepImmersedFlow>);
static_assert(!std::is_copy_constructible_v<FixedStepImmersedFlow>);
static_assert(!std::is_copy_assignable_v<FixedStepImmersedFlow>);
static_assert(std::is_nothrow_move_constructible_v<FixedStepImmersedFlow>);
static_assert(!std::is_move_assignable_v<FixedStepImmersedFlow>);
static_assert(std::is_nothrow_destructible_v<FixedStepImmersedFlow>);
static_assert(
    static_cast<unsigned>(ImmersedLinearSolvePhase::predictor_momentum) == 0U);
static_assert(static_cast<unsigned>(
                  ImmersedLinearSolvePhase::pressure_affine_momentum) == 1U);
static_assert(static_cast<unsigned>(
                  ImmersedLinearSolvePhase::pressure_homogeneous_momentum) ==
              2U);
static_assert(static_cast<unsigned>(ImmersedLinearSolvePhase::pressure_outer) ==
              3U);
static_assert(static_cast<unsigned>(
                  ImmersedLinearSolvePhase::pressure_independent_residual) ==
              4U);
static_assert(static_cast<unsigned>(
                  ImmersedLinearSolvePhase::pressure_compact_preconditioner) ==
              5U);
static_assert(std::is_same_v<
              decltype(std::declval<const FixedStepImmersedFlow &>()
                           .performance_counters()),
              hundun::diagnostics::Stage3PerformanceCounters>);

using Create = FixedStepImmersedFlow (*)(
    const hundun::runtime::StructuredDecomposition &,
    const hundun::mesh::MeshTopology &, const hundun::mesh::MeshGeometry &,
    const hundun::boundary::BoundaryRegistry &,
    const hundun::immersed::ImmersedDomain *,
    const hundun::immersed::GhostStencilPlan *,
    const hundun::immersed::WallQuadraturePlan *,
    const hundun::immersed::LocalFlowPatternTransform *,
    const hundun::les::WaleModel *, const hundun::runtime::MpiContext &,
    hundun::execution::ExecutionContext &, hundun::runtime::HaloExchange &,
    const hundun::linear::LinearSolver &,
    std::array<hundun::linear::Preconditioner *, 3>,
    const hundun::linear::LinearSolver &, hundun::linear::Preconditioner &);

using CreateWithDensity = FixedStepImmersedFlow (*)(
    const hundun::runtime::StructuredDecomposition &,
    const hundun::mesh::MeshTopology &, const hundun::mesh::MeshGeometry &,
    const hundun::boundary::BoundaryRegistry &,
    const hundun::immersed::ImmersedDomain *,
    const hundun::immersed::GhostStencilPlan *,
    const hundun::immersed::WallQuadraturePlan *,
    const hundun::immersed::LocalFlowPatternTransform *,
    const hundun::les::WaleModel *, ImmersedFlowDensitySetup,
    const hundun::runtime::MpiContext &,
    hundun::execution::ExecutionContext &, hundun::runtime::HaloExchange &,
    const hundun::linear::LinearSolver &,
    std::array<hundun::linear::Preconditioner *, 3>,
    const hundun::linear::LinearSolver &, hundun::linear::Preconditioner &);

static_assert(std::is_same_v<
              decltype(static_cast<Create>(&FixedStepImmersedFlow::create)),
              Create>);
static_assert(
    std::is_same_v<decltype(static_cast<CreateWithDensity>(
                       &FixedStepImmersedFlow::create)),
                   CreateWithDensity>);

} // namespace

int main() {
  const auto create = static_cast<Create>(&FixedStepImmersedFlow::create);
  const auto create_with_density =
      static_cast<CreateWithDensity>(&FixedStepImmersedFlow::create);
  return create == nullptr || create_with_density == nullptr ? 1 : 0;
}
