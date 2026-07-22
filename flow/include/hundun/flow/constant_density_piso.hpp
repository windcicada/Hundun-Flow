// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/flow_state.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace hundun::flow {

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class ConstantDensityPisoTestAccess;
class MaterialDensityPisoTestAccess;
}
#endif

enum class StepAttemptDisposition : std::uint8_t {
  committed,
  recoverable_failure,
  non_retryable_failure
};

enum class StepFailureReason : std::uint8_t {
  none,
  invalid_input,
  momentum_linear_solve,
  pressure_linear_solve,
  non_finite_trial,
  boundary_backflow,
  transport_failure,
  final_momentum_residual,
  final_transport_residual,
  final_conservation_defect,
  final_continuity_residual,
  final_pressure_residual,
  collective_operation,
  density_closure_failure
};

enum class PressureCorrectionDisposition : std::uint8_t {
  accepted,
  recoverable_failure,
  non_retryable_failure
};

struct PressureCorrectionReport final {
  // Successful corrections report {accepted, none, -1, accepted=true}.
  // Expected numerical solve/residual failures are returned as recoverable
  // pressure_linear_solve results.  Invalid input, layout, capability, or
  // construction failures are non-retryable invalid_input results; MPI
  // operation failures are non-retryable collective_operation results.
  // lowest_failing_rank is -1 only when no reliable rank evidence exists.
  PressureCorrectionDisposition disposition{
      PressureCorrectionDisposition::non_retryable_failure};
  StepFailureReason reason{StepFailureReason::invalid_input};
  int lowest_failing_rank{-1};
  linear::SolveReport solve;
  double independent_residual_l2{};
  double rhs_l2{};
  bool accepted{};
};

struct StepAttemptReport final {
  StepAttemptDisposition disposition{
      StepAttemptDisposition::non_retryable_failure};
  StepFailureReason reason{StepFailureReason::invalid_input};
  int lowest_failing_rank{-1};
  std::uint32_t pressure_corrector_count{};
  double attempted_dt_s{};
  double suggested_dt_s{};
  MomentumPredictorReport momentum{};
  std::array<linear::SolveReport, 2> pressure{};
  double final_continuity_normalized_l2{};
  double final_pressure_residual_l2{};
  std::array<double, 3> final_momentum_normalized_l2{};
  std::vector<double> final_transport_normalized_l2;
  double final_mass_relative_conservation_defect{};
  std::array<double, 3> final_momentum_relative_conservation_defect{};
  std::vector<double> final_transport_relative_conservation_defect;
  std::optional<boundary::OutletBackflowEvidence> final_backflow_evidence;
};

struct ConstantDensityTransportSpec final {
  runtime::FieldId field{};
  finite_volume::FiniteVolumeQuantity quantity{};
  double diffusivity_kg_per_m_s{};
};

class PisoCoupler final {
public:
  // All collaborators are borrowed and must outlive the returned coupler.
  // Calls on one coupler object must not overlap.
  static PisoCoupler
  create(const runtime::StructuredDecomposition &decomposition,
         const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
         const boundary::BoundaryRegistry &boundaries,
         const runtime::MpiContext &mpi,
         execution::ExecutionContext &execution_context,
         runtime::HaloExchange &cell_halo,
         const linear::LinearSolver &pressure_solver,
         linear::Preconditioner &pressure_preconditioner);

  ~PisoCoupler() noexcept;
  PisoCoupler(PisoCoupler &&) noexcept;
  PisoCoupler &operator=(PisoCoupler &&) = delete;
  PisoCoupler(const PisoCoupler &) = delete;
  PisoCoupler &operator=(const PisoCoupler &) = delete;

  PressureCorrectionReport
  correct(FlowState &state, double rho_ref,
          const runtime::FieldView<const double> &actual_momentum_diagonal,
          const linear::SolveControl &control) const;

  PressureCorrectionReport correct_material_density(
      FlowState &state, const MomentumTimeStencil &stencil,
      const runtime::FieldView<const double> &actual_momentum_diagonal,
      const linear::SolveControl &control) const;

private:
  struct MaterialPressureAssessment final {
    PressureCorrectionDisposition disposition{
        PressureCorrectionDisposition::non_retryable_failure};
    StepFailureReason reason{StepFailureReason::invalid_input};
    int lowest_failing_rank{-1};
    double independent_residual_l2{};
    double rhs_l2{};
    double normalized_residual{};
    bool residual_available{};
    bool accepted{};
  };
  struct Impl;
  explicit PisoCoupler(std::unique_ptr<Impl>) noexcept;
  PressureCorrectionReport correct_throwing(
      FlowState &state, double rho_ref,
      const runtime::FieldView<const double> &actual_momentum_diagonal,
      const linear::SolveControl &control,
      PressureCorrectionReport &result) const;
  PressureCorrectionReport correct_common_throwing(
      FlowState &state, double rho_ref,
      const MomentumTimeStencil *material_stencil,
      const runtime::FieldView<const double> &actual_momentum_diagonal,
      const linear::SolveControl &control,
      PressureCorrectionReport &result) const;
  MaterialPressureAssessment assess_final_material_density_pressure(
      FlowState &state, const MomentumTimeStencil &stencil,
      const linear::SolveControl &control) const;
  void prepare_material_density_assessment();
  std::unique_ptr<Impl> impl_;

  friend class FixedStepMaterialDensityFlow;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::ConstantDensityPisoTestAccess;
  friend class test::MaterialDensityPisoTestAccess;
#endif
};

class FixedStepConstantDensityFlow final {
public:
  // All collaborators are borrowed and must outlive the returned flow.
  // Calls on one flow object must not overlap.
  static FixedStepConstantDensityFlow
  create(const runtime::StructuredDecomposition &decomposition,
         const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
         const boundary::BoundaryRegistry &boundaries,
         const runtime::MpiContext &mpi,
         execution::ExecutionContext &execution_context,
         runtime::HaloExchange &cell_halo,
         const linear::LinearSolver &momentum_solver,
         std::array<linear::Preconditioner *, 3> momentum_preconditioners,
         const linear::LinearSolver &pressure_solver,
         linear::Preconditioner &pressure_preconditioner,
         std::vector<ConstantDensityTransportSpec> transported_fields = {});

  ~FixedStepConstantDensityFlow() noexcept;
  FixedStepConstantDensityFlow(FixedStepConstantDensityFlow &&) noexcept;
  FixedStepConstantDensityFlow &
  operator=(FixedStepConstantDensityFlow &&) = delete;
  FixedStepConstantDensityFlow(const FixedStepConstantDensityFlow &) = delete;
  FixedStepConstantDensityFlow &
  operator=(const FixedStepConstantDensityFlow &) = delete;

  StepAttemptReport attempt(FlowState &state, double rho_ref, double mu,
                            const MomentumTimeStencil &stencil,
                            const linear::SolveControl &momentum_control,
                            const linear::SolveControl &pressure_control) const;

private:
  struct Impl;
  explicit FixedStepConstantDensityFlow(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::ConstantDensityPisoTestAccess;
#endif
};

} // namespace hundun::flow
