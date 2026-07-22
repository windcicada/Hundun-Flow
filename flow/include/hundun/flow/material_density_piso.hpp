// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/flow/material_density_transport.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string_view>

namespace hundun::flow {

namespace detail {
struct DensityClosureBridge;
struct DensityClosureHooks;
} // namespace detail

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class MaterialDensityPisoTestAccess;
}
#endif

class MaterialDensityStepAttemptReport final {
public:
  MaterialDensityStepAttemptReport(const MaterialDensityStepAttemptReport &) =
      default;
  MaterialDensityStepAttemptReport(
      MaterialDensityStepAttemptReport &&) noexcept = default;
  MaterialDensityStepAttemptReport &
  operator=(const MaterialDensityStepAttemptReport &) = default;
  MaterialDensityStepAttemptReport &
  operator=(MaterialDensityStepAttemptReport &&) noexcept = default;

  const StepAttemptReport &flow() const noexcept;
  bool material_report_available() const noexcept;
  const MaterialDensityTransportReport &material_report() const;
  MaterialTransportFailureReason material_failure_reason() const noexcept;
  std::uint64_t material_field_count() const noexcept;
  runtime::FieldId shared_face_mass_flux_field() const noexcept;
  MaterialFluxProvenance flux_provenance() const noexcept;
  std::uint64_t attempt_identity() const noexcept;
  bool final_continuity_residual_available() const noexcept;
  bool final_pressure_residual_available() const noexcept;
  double final_pressure_normalized_residual() const noexcept;
  const std::array<std::uint8_t, 3> &
  final_momentum_residual_availability() const noexcept;
  bool mass_conservation_available() const noexcept;
  const std::array<std::uint8_t, 3> &
  momentum_conservation_availability() const noexcept;

private:
  MaterialDensityStepAttemptReport() = default;
  std::uint64_t compute_seal() const noexcept;
  bool semantic_valid() const noexcept;
  void seal() noexcept;
  bool authenticated() const noexcept;

  StepAttemptReport flow_{};
  std::optional<MaterialDensityTransportReport> material_report_;
  MaterialTransportFailureReason material_failure_reason_{
      MaterialTransportFailureReason::none};
  std::uint64_t material_field_count_{};
  runtime::FieldId shared_face_mass_flux_field_{};
  MaterialFluxProvenance flux_provenance_{MaterialFluxProvenance::predictor};
  std::uint64_t attempt_identity_{};
  std::uint64_t material_attempt_identity_{};
  std::uint64_t material_finalization_identity_{};
  bool final_continuity_residual_available_{};
  bool final_pressure_residual_available_{};
  double final_pressure_normalized_residual_{};
  std::array<std::uint8_t, 3> final_momentum_residual_available_{};
  bool mass_conservation_available_{};
  std::array<std::uint8_t, 3> momentum_conservation_available_{};
  bool closure_origin_{};
  bool post_closure_evidence_available_{};
  std::optional<MaterialDensityTransportReport> post_closure_report_;
  std::uint64_t seal_{};

  friend class FixedStepMaterialDensityFlow;
  friend struct detail::DensityClosureBridge;
  friend class MaterialDensityFlowDiagnosticSource;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::MaterialDensityPisoTestAccess;
#endif
};

enum class MaterialDensityDiagnosticEntity : std::uint8_t { cell, face };

class MaterialDensityFlowDiagnosticSource final {
public:
  ~MaterialDensityFlowDiagnosticSource() noexcept;
  MaterialDensityFlowDiagnosticSource(
      MaterialDensityFlowDiagnosticSource &&) noexcept;
  MaterialDensityFlowDiagnosticSource &
  operator=(MaterialDensityFlowDiagnosticSource &&) = delete;
  MaterialDensityFlowDiagnosticSource(
      const MaterialDensityFlowDiagnosticSource &) = delete;
  MaterialDensityFlowDiagnosticSource &
  operator=(const MaterialDensityFlowDiagnosticSource &) = delete;

  std::size_t fingerprint_field_count() const;
  std::string_view fingerprint_field_id(std::size_t field) const;
  std::string_view field_unit(std::size_t field) const;
  MaterialDensityDiagnosticEntity field_entity(std::size_t field) const;
  std::size_t field_component_count(std::size_t field) const;
  std::size_t field_item_count(std::size_t field) const;
  std::uint64_t field_global_id(std::size_t field, std::size_t item) const;
  double field_value(std::size_t field, std::size_t item,
                     std::size_t component) const;
  std::string_view owned_cell_layout_fingerprint() const;
  std::string_view global_cell_layout_fingerprint() const;
  std::string_view owned_face_layout_fingerprint() const;
  std::string_view global_face_layout_fingerprint() const;
  int relative_rank() const;
  std::uint64_t committed_step() const;
  double committed_time_s() const;
  runtime::Int3 global_cell_extent() const;
  runtime::Box3 owned_global_box() const;
  std::size_t owned_cell_count() const;
  std::size_t canonical_owned_face_count() const;
  double cell_volume_m3(std::size_t owned_cell) const;
  const MaterialDensityStepAttemptReport &report() const;

private:
  struct Impl;
  explicit MaterialDensityFlowDiagnosticSource(std::unique_ptr<Impl>) noexcept;
  void validate() const;
  std::unique_ptr<Impl> impl_;
  friend class FixedStepMaterialDensityFlow;
};

class FixedStepMaterialDensityFlow final {
public:
  static FixedStepMaterialDensityFlow
  create(const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
         const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
         const runtime::MpiContext &, execution::ExecutionContext &,
         runtime::HaloExchange &, const linear::LinearSolver &momentum_solver,
         std::array<linear::Preconditioner *, 3> momentum_preconditioners,
         const linear::LinearSolver &pressure_solver,
         linear::Preconditioner &pressure_preconditioner,
         const runtime::FieldRegistry &, FlowFieldIds,
         MaterialDensityTransportSpec);

  ~FixedStepMaterialDensityFlow() noexcept;
  FixedStepMaterialDensityFlow(FixedStepMaterialDensityFlow &&) noexcept;
  FixedStepMaterialDensityFlow &
  operator=(FixedStepMaterialDensityFlow &&) = delete;
  FixedStepMaterialDensityFlow(const FixedStepMaterialDensityFlow &) = delete;
  FixedStepMaterialDensityFlow &
  operator=(const FixedStepMaterialDensityFlow &) = delete;

  MaterialDensityStepAttemptReport
  attempt(FlowState &, double mu, const MomentumTimeStencil &,
          const linear::SolveControl &momentum_control,
          const linear::SolveControl &pressure_control) const;

  MaterialDensityFlowDiagnosticSource
  diagnostic_source(const FlowState &,
                    const MaterialDensityStepAttemptReport &) const;

private:
  static FixedStepMaterialDensityFlow create_open_capable(
      const runtime::StructuredDecomposition &, const mesh::MeshTopology &,
      const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
      const runtime::MpiContext &, execution::ExecutionContext &,
      runtime::HaloExchange &, const linear::LinearSolver &,
      std::array<linear::Preconditioner *, 3>, const linear::LinearSolver &,
      linear::Preconditioner &, const runtime::FieldRegistry &, FlowFieldIds,
      MaterialDensityTransportSpec);
  MaterialDensityStepAttemptReport
  attempt_common(FlowState &, double mu, const MomentumTimeStencil &,
                 const linear::SolveControl &momentum_control,
                 const linear::SolveControl &pressure_control,
                 const detail::DensityClosureHooks *) const;
  struct Impl;
  explicit FixedStepMaterialDensityFlow(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class MaterialDensityFlowDiagnosticSource;
  friend struct detail::DensityClosureBridge;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::MaterialDensityPisoTestAccess;
#endif
};

} // namespace hundun::flow
