// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case.hpp"
#include "hundun/diag_stage3_performance.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/ib_wall_force.hpp"
#include "hundun/les_wale.hpp"

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

class ImmersedFlowDiagnosticSource;

namespace detail {
struct ImmersedFlowCheckpointAccess;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
struct ImmersedFlowAccess;
#endif
}

struct ForceAttemptReport final {
  immersed::ForceComponents operator_force;
  immersed::ForceComponents budget_reaction;
  immersed::ForceComponents surface_traction;
  immersed::ForceComponents consistency;
};

struct ImmersedFlowStepAttemptReport final {
  DensityStepAttemptReport base;
  std::optional<ForceAttemptReport> force;
  // Present only after the same-attempt WALE authority commits successfully.
  std::optional<les::WaleSummary> wale;
};

struct ImmersedFlowPhysics final {
  config::DensityModel density_model{config::DensityModel::constant};
  double rho_ref_kg_per_m3{};
  double dynamic_viscosity_pa_s{};
  std::optional<double> cp_J_per_kg_K;
  std::optional<double> gas_constant_J_per_kg_K;
  std::optional<double> thermodynamic_pressure_pa;
};

struct ImmersedFlowDensitySetup final {
  config::DensityModel model{config::DensityModel::constant};
  const runtime::FieldRegistry *registry{};
  FlowFieldIds fields{};
  std::optional<MaterialDensityTransportSpec> material_transport;
  IdealGasClosure *ideal_gas_closure{};
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

  static FixedStepImmersedFlow
  create(const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
         const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
         const immersed::ImmersedDomain *, const immersed::GhostStencilPlan *,
         const immersed::WallQuadraturePlan *,
         const immersed::LocalFlowPatternTransform *, const les::WaleModel *,
         ImmersedFlowDensitySetup, const runtime::MpiContext &,
         execution::ExecutionContext &, runtime::HaloExchange &,
         const linear::LinearSolver &,
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

  // The returned source borrows this flow and is valid only until the next
  // attempt or checkpoint restore. No diagnostic sampling occurs otherwise.
  ImmersedFlowDiagnosticSource
  diagnostic_source(const FlowState &,
                    const ImmersedFlowStepAttemptReport &) const;
  diagnostics::Stage3PerformanceCounters
  performance_counters() const noexcept;

private:
  struct Impl;
  explicit FixedStepImmersedFlow(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;

  friend struct detail::ImmersedFlowCheckpointAccess;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend struct detail::ImmersedFlowAccess;
#endif
};

} // namespace hundun::flow
