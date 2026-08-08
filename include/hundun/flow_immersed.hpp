// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/ib_wall_force.hpp"

#include <array>
#include <memory>
#include <optional>

namespace hundun::boundary {
class BoundaryRegistry;
}

namespace hundun::immersed {
class GhostStencilPlan;
class ImmersedDomain;
class LocalFlowPatternTransform;
class WallQuadraturePlan;
} // namespace hundun::immersed

namespace hundun::les {
class WaleModel;
}

namespace hundun::flow {

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace detail {
struct ImmersedFlowAccess;
}
#endif

struct ForceAttemptReport final {
  immersed::ForceComponents operator_force;
  immersed::ForceComponents budget_reaction;
  immersed::ForceComponents surface_traction;
  immersed::ForceComponents consistency;
};

struct ImmersedFlowStepAttemptReport final {
  DensityStepAttemptReport base;
  std::optional<ForceAttemptReport> force;
};

struct ImmersedFlowPhysics final {
  config::DensityModel density_model{config::DensityModel::constant};
  double rho_ref_kg_per_m3{};
  double dynamic_viscosity_pa_s{};
  std::optional<double> cp_J_per_kg_K;
  std::optional<double> gas_constant_J_per_kg_K;
  std::optional<double> thermodynamic_pressure_pa;
};

class FixedStepImmersedFlow final {
public:
  static FixedStepImmersedFlow
  create(const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
         const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
         const immersed::ImmersedDomain *, const immersed::GhostStencilPlan *,
         const immersed::WallQuadraturePlan *,
         const immersed::LocalFlowPatternTransform *, const les::WaleModel *,
         const runtime::MpiContext &, execution::ExecutionContext &,
         runtime::HaloExchange &, const linear::LinearSolver &,
         std::array<linear::Preconditioner *, 3>, const linear::LinearSolver &,
         linear::Preconditioner &);

  ~FixedStepImmersedFlow() noexcept;
  FixedStepImmersedFlow(FixedStepImmersedFlow &&) noexcept;
  FixedStepImmersedFlow &operator=(FixedStepImmersedFlow &&) = delete;
  FixedStepImmersedFlow(const FixedStepImmersedFlow &) = delete;
  FixedStepImmersedFlow &operator=(const FixedStepImmersedFlow &) = delete;

  ImmersedFlowStepAttemptReport attempt(FlowState &, const ImmersedFlowPhysics &,
                                  const MomentumTimeStencil &,
                                  const linear::SolveControl &,
                                  const linear::SolveControl &) const;

private:
  struct Impl;
  explicit FixedStepImmersedFlow(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend struct detail::ImmersedFlowAccess;
#endif
};

} // namespace hundun::flow
