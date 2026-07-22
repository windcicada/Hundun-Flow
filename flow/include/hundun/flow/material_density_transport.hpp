// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/momentum_predictor.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/field_descriptor.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace hundun::boundary {
class BoundaryRegistry;
}
namespace hundun::mesh {
class MeshGeometry;
}
namespace hundun::runtime {
class FieldAccessPlan;
class FieldRegistry;
class FieldStorage;
class HaloExchange;
class MpiContext;
class StructuredDecomposition;
} // namespace hundun::runtime

namespace hundun::flow {

class FlowState;
struct FlowFieldIds;
class FixedStepMaterialDensityFlow;
class MaterialDensityStepAttemptReport;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace test {
class MaterialDensityTransportTestAccess;
}
#endif

enum class MaterialFluxProvenance : std::uint8_t {
  predictor,
  provisional,
  final_corrected
};

class MaterialFaceMassFlux final {
public:
  static MaterialFaceMassFlux
  acquire(const runtime::FieldRegistry &registry,
          const runtime::FieldStorage &storage,
          const runtime::FieldAccessPlan &access_plan, runtime::PhaseId phase,
          runtime::ActorId actor, runtime::FieldId field,
          const mesh::MeshTopology &topology,
          MaterialFluxProvenance provenance);

  ~MaterialFaceMassFlux() noexcept;
  MaterialFaceMassFlux(MaterialFaceMassFlux &&) noexcept;
  MaterialFaceMassFlux &operator=(MaterialFaceMassFlux &&) = delete;
  MaterialFaceMassFlux(const MaterialFaceMassFlux &) = delete;
  MaterialFaceMassFlux &operator=(const MaterialFaceMassFlux &) = delete;

  runtime::FieldId field_id() const noexcept;
  std::size_t face_count() const noexcept;
  MaterialFluxProvenance provenance() const noexcept;

private:
  struct Impl;
  explicit MaterialFaceMassFlux(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class MaterialDensityTransport;
  friend class FixedStepMaterialDensityFlow;
};

struct MaterialDensityTransportSpec final {
  runtime::FieldId enthalpy_density{};
  double enthalpy_diffusivity_kg_per_m_s{};
  std::vector<runtime::FieldId> scalar_densities;
  std::vector<double> scalar_diffusivities_kg_per_m_s;
};

enum class MaterialTransportDisposition : std::uint8_t {
  finalized,
  recoverable_failure,
  non_retryable_failure
};
enum class MaterialTransportFailureReason : std::uint8_t {
  none,
  invalid_input,
  non_finite_state,
  non_positive_density,
  final_density_residual,
  final_transport_residual,
  final_conservation_defect,
  collective_operation
};

class MaterialDensityTransportReport final {
public:
  MaterialDensityTransportReport(const MaterialDensityTransportReport &) =
      default;
  MaterialDensityTransportReport(MaterialDensityTransportReport &&) noexcept =
      default;
  MaterialDensityTransportReport &
  operator=(const MaterialDensityTransportReport &) = default;
  MaterialDensityTransportReport &
  operator=(MaterialDensityTransportReport &&) noexcept = default;
  ~MaterialDensityTransportReport() noexcept = default;

  MaterialTransportDisposition disposition() const noexcept;
  MaterialTransportFailureReason reason() const noexcept;
  int lowest_failing_rank() const noexcept;
  const MomentumTimeStencil &stencil() const noexcept;
  MaterialFluxProvenance flux_provenance() const noexcept;
  std::uint64_t attempt_identity() const noexcept;
  std::uint64_t finalization_identity() const noexcept;
  runtime::FieldId shared_face_mass_flux_field() const noexcept;
  bool density_residual_available() const noexcept;
  double density_normalized_l2() const noexcept;
  const std::vector<std::uint8_t> &
  transport_residual_availability() const noexcept;
  const std::vector<double> &transport_normalized_l2() const noexcept;
  bool mass_conservation_available() const noexcept;
  double mass_relative_conservation_defect() const noexcept;
  const std::vector<std::uint8_t> &
  transport_conservation_availability() const noexcept;
  const std::vector<double> &
  transport_relative_conservation_defect() const noexcept;
  bool minimum_density_available() const noexcept;
  double minimum_density_kg_per_m3() const noexcept;
  mesh::GlobalCellId minimum_density_global_cell() const noexcept;
  int minimum_density_rank() const noexcept;

private:
  MaterialDensityTransportReport() = default;
  std::uint64_t compute_seal() const noexcept;
  void seal() noexcept;
  bool authenticated() const noexcept;

  MaterialTransportDisposition disposition_{
      MaterialTransportDisposition::non_retryable_failure};
  MaterialTransportFailureReason reason_{
      MaterialTransportFailureReason::invalid_input};
  int lowest_failing_rank_{-1};
  MomentumTimeStencil stencil_{};
  MaterialFluxProvenance flux_provenance_{MaterialFluxProvenance::predictor};
  std::uint64_t attempt_identity_{};
  std::uint64_t finalization_identity_{};
  runtime::FieldId shared_face_mass_flux_field_{};
  bool density_residual_available_{};
  double density_normalized_l2_{};
  std::vector<std::uint8_t> transport_residual_available_;
  std::vector<double> transport_normalized_l2_;
  bool mass_conservation_available_{};
  double mass_relative_conservation_defect_{};
  std::vector<std::uint8_t> transport_conservation_available_;
  std::vector<double> transport_relative_conservation_defect_;
  bool minimum_density_available_{};
  double minimum_density_kg_per_m3_{};
  mesh::GlobalCellId minimum_density_global_cell_{};
  int minimum_density_rank_{-1};
  std::uint64_t seal_{};

  friend class MaterialDensityTransport;
  friend class MaterialDensityDiagnosticSource;
  friend class MaterialDensityStepAttemptReport;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::MaterialDensityTransportTestAccess;
#endif
};

class MaterialDensityDiagnosticSource;

class MaterialDensityTransport final {
public:
  static MaterialDensityTransport
  create(const runtime::FieldRegistry &registry,
         const runtime::StructuredDecomposition &decomposition,
         const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
         const boundary::BoundaryRegistry &boundaries,
         const runtime::MpiContext &mpi, runtime::HaloExchange &halo,
         FlowFieldIds expected_fields,
         MaterialDensityTransportSpec specification);

  ~MaterialDensityTransport() noexcept;
  MaterialDensityTransport(MaterialDensityTransport &&) noexcept;
  MaterialDensityTransport &operator=(MaterialDensityTransport &&) = delete;
  MaterialDensityTransport(const MaterialDensityTransport &) = delete;
  MaterialDensityTransport &
  operator=(const MaterialDensityTransport &) = delete;

  MaterialDensityTransportReport
  finalize_trial(FlowState &state, const MaterialFaceMassFlux &final_mass_flux,
                 const MomentumTimeStencil &stencil) const;

  MaterialDensityDiagnosticSource
  diagnostic_source(const FlowState &state,
                    const MaterialDensityTransportReport &report) const;

private:
  struct StagingResult final {
    MaterialTransportDisposition disposition{
        MaterialTransportDisposition::non_retryable_failure};
    MaterialTransportFailureReason reason{
        MaterialTransportFailureReason::invalid_input};
    int lowest_failing_rank{-1};
  };
  StagingResult stage_trial(FlowState &, const MaterialFaceMassFlux &,
                            const MomentumTimeStencil &) const;
  void prepare_task20_attempt() const;
  struct Impl;
  explicit MaterialDensityTransport(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class MaterialDensityDiagnosticSource;
  friend class FixedStepMaterialDensityFlow;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  friend class test::MaterialDensityTransportTestAccess;
#endif
};

class MaterialDensityDiagnosticSource final {
public:
  ~MaterialDensityDiagnosticSource() noexcept;
  MaterialDensityDiagnosticSource(MaterialDensityDiagnosticSource &&) noexcept;
  MaterialDensityDiagnosticSource &
  operator=(MaterialDensityDiagnosticSource &&) = delete;
  MaterialDensityDiagnosticSource(const MaterialDensityDiagnosticSource &) =
      delete;
  MaterialDensityDiagnosticSource &
  operator=(const MaterialDensityDiagnosticSource &) = delete;

  std::size_t fingerprint_field_count() const;
  std::string_view fingerprint_field_id(std::size_t index) const;
  std::string_view field_unit(std::size_t index) const;
  std::string_view owned_cell_layout_fingerprint() const;
  std::string_view global_cell_layout_fingerprint() const;
  runtime::Int3 global_cell_extent() const;
  runtime::Box3 owned_global_box() const;
  std::size_t owned_cell_count() const;
  mesh::GlobalCellId global_cell_id(std::size_t local_cell) const;
  double cell_volume_m3(std::size_t local_cell) const;
  double field_value(std::size_t field, std::size_t local_cell) const;
  const MaterialDensityTransportReport &report() const;

private:
  struct Impl;
  explicit MaterialDensityDiagnosticSource(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class MaterialDensityTransport;
};

} // namespace hundun::flow
