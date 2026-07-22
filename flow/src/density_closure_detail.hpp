// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/material_density_piso.hpp"
#include "hundun/runtime/field_view.hpp"

#include <cstdint>

namespace hundun::flow {

class IdealGasClosure;
class IdealGasClosureDiagnosticSource;

namespace detail {

enum class DensityClosureStage : std::uint8_t { predictor, provisional, final };

struct DensityClosureEvaluation final {
  bool accepted{};
  bool recoverable{};
  int lowest_failing_rank{-1};
};

struct DensityClosureHooks final {
  void *object{};
  runtime::FieldId enthalpy_density{};
  double enthalpy_rate_J_per_kg_s{};
  void (*begin)(void *, const FlowState &, std::uint64_t){};
  DensityClosureEvaluation (*evaluate)(void *, FlowState &,
                                       DensityClosureStage){};
  void (*prepare)(void *){};
  void (*publish)(void *) noexcept {};
  void (*rollback)(void *) noexcept {};

  explicit operator bool() const noexcept { return object != nullptr; }
};

struct DensityClosureBridge final {
  static FixedStepMaterialDensityFlow create_open_capable(
      const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
      const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
      const runtime::MpiContext &, execution::ExecutionContext &,
      runtime::HaloExchange &, const linear::LinearSolver &,
      std::array<linear::Preconditioner *, 3>, const linear::LinearSolver &,
      linear::Preconditioner &, const runtime::FieldRegistry &, FlowFieldIds,
      MaterialDensityTransportSpec);

  static MaterialDensityStepAttemptReport
  attempt(const FixedStepMaterialDensityFlow &, FlowState &, double,
          const MomentumTimeStencil &, const linear::SolveControl &,
          const linear::SolveControl &, const DensityClosureHooks &);

  static bool
  report_authenticated(const MaterialDensityStepAttemptReport &) noexcept;
  static std::uint64_t
  report_seal(const MaterialDensityStepAttemptReport &) noexcept;
  static MaterialDensityStepAttemptReport make_report();
  static bool post_eos_evidence_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
};

struct DensityClosureAdapter final {
  static DensityClosureHooks bind(IdealGasClosure &, runtime::FieldId,
                                  double enthalpy_rate_J_per_kg_s);
  static void begin(void *, const FlowState &, std::uint64_t);
  static DensityClosureEvaluation evaluate(void *, FlowState &,
                                           DensityClosureStage);
  static void prepare(void *);
  static void publish(void *) noexcept;
  static void rollback(void *) noexcept;
  static double gas_constant_J_per_kg_K(const IdealGasClosure &) noexcept;
  static double cp_J_per_kg_K(const IdealGasClosure &) noexcept;
};

struct DensityClosureReadSession final {
  runtime::FieldView<const double> density;
  runtime::FieldView<const double> enthalpy_density;
  runtime::Box3 owned_box{};
  runtime::Int3 global_extent{};
};

struct DensityClosureDiagnosticAccess final {
  static DensityClosureReadSession
  acquire_committed(const IdealGasClosureDiagnosticSource &);
  static double
  gas_constant_J_per_kg_K(const IdealGasClosureDiagnosticSource &);
  static double cp_J_per_kg_K(const IdealGasClosureDiagnosticSource &);
  static const runtime::MpiContext &
  mpi(const IdealGasClosureDiagnosticSource &);
};

} // namespace detail
} // namespace hundun::flow
