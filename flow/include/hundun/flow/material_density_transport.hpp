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

struct MaterialDensityTransportReport final {
  MaterialTransportDisposition disposition{
      MaterialTransportDisposition::non_retryable_failure};
  MaterialTransportFailureReason reason{
      MaterialTransportFailureReason::invalid_input};
  int lowest_failing_rank{-1};
  MomentumTimeStencil stencil{};
  MaterialFluxProvenance flux_provenance{MaterialFluxProvenance::predictor};
  std::uint64_t attempt_identity{};
  std::uint64_t finalization_identity{};
  runtime::FieldId shared_face_mass_flux_field{};
  double density_normalized_l2{};
  std::vector<double> transport_normalized_l2;
  double mass_relative_conservation_defect{};
  std::vector<double> transport_relative_conservation_defect;
  double minimum_density_kg_per_m3{};
  mesh::GlobalCellId minimum_density_global_cell{};
  int minimum_density_rank{-1};
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
  struct Impl;
  explicit MaterialDensityTransport(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend class MaterialDensityDiagnosticSource;
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
