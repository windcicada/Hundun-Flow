// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_field_view.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::flow {
class FixedStepMaterialDensityFlow;
class MaterialFaceMassFlux;
class PisoCoupler;
} // namespace hundun::flow
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
namespace hundun::finite_volume::test {
class PreparedFaceMassFluxForTest;
}
#endif

namespace hundun::finite_volume {

class ImmersedReconstruction;
class ImmersedOperatorAdapter;

enum class GradientScheme : std::uint8_t {
  green_gauss,
  weighted_least_squares
};

enum class FiniteVolumeQuantityKind : std::uint8_t {
  density,
  enthalpy,
  scalar,
  velocity,
  pressure
};

struct FiniteVolumeQuantity final {
  FiniteVolumeQuantityKind kind;
  std::size_t scalar_index;

  static FiniteVolumeQuantity density() noexcept;
  static FiniteVolumeQuantity enthalpy() noexcept;
  static FiniteVolumeQuantity scalar(std::size_t index) noexcept;
  static FiniteVolumeQuantity velocity() noexcept;
  static FiniteVolumeQuantity pressure() noexcept;
};

double monotonized_central(double left, double right) noexcept;

runtime::FieldDescriptor face_mass_flux_descriptor();
runtime::FieldId declare_face_mass_flux(runtime::FieldRegistry &registry);
void require_face_mass_flux_field(const runtime::FieldRegistry &registry,
                                  runtime::FieldId field);

struct PhysicalBoundaryMomentumContribution final {
  mesh::GlobalFaceId global_face_id{};
  std::array<double, 3> convective{};
  std::array<double, 3> viscous{};
};

struct PhysicalBoundaryPressureContribution final {
  mesh::GlobalFaceId global_face_id{};
  std::array<double, 3> pressure{};
};

struct PhysicalBoundaryTransportContribution final {
  mesh::GlobalFaceId global_face_id{};
  double convective{};
  double diffusive{};
};

// Physical-boundary contribution buffers contain each unique, owner-owned,
// nonperiodic physical face in stable topology order.  Area orientation is
// owner-outward.  Convective values use equation units (momentum: N;
// transport: quantity flux), pressure uses N, and viscous/diffusive values
// are equation contributions: the negatives of the physical traction or
// diffusive flux under the residual convention used by this class.

class FaceMassFlux final {
public:
  static FaceMassFlux acquire(const runtime::FieldRegistry &registry,
                              const runtime::FieldStorage &storage,
                              const runtime::FieldAccessPlan &access_plan,
                              runtime::PhaseId phase, runtime::ActorId actor,
                              runtime::FieldId field,
                              const mesh::MeshTopology &topology);

  ~FaceMassFlux() noexcept;
  FaceMassFlux(FaceMassFlux &&) noexcept;
  FaceMassFlux &operator=(FaceMassFlux &&) = delete;
  FaceMassFlux(const FaceMassFlux &) = delete;
  FaceMassFlux &operator=(const FaceMassFlux &) = delete;

  runtime::FieldId field_id() const noexcept;
  std::size_t face_count() const noexcept;

private:
  struct State;
  struct OwnedState;
  struct PreparedState;
  using PreparedStatePtr =
      std::unique_ptr<PreparedState, void (*)(PreparedState *)>;
  static PreparedStatePtr prepare(const mesh::MeshTopology &topology);
  static void destroy_prepared(PreparedState *) noexcept;
  static FaceMassFlux
  bind_prepared(PreparedState &, const runtime::FieldRegistry &,
                const runtime::FieldStorage &, const runtime::FieldAccessPlan &,
                runtime::PhaseId, runtime::ActorId, runtime::FieldId,
                const mesh::MeshTopology &);
  explicit FaceMassFlux(std::unique_ptr<OwnedState>) noexcept;
  FaceMassFlux(State *, PreparedState *) noexcept;
  void
  require_immersed_reconstruction_compatible(const mesh::MeshTopology &) const;
  double value_for_immersed_reconstruction(mesh::LocalFaceId local_face) const;
  void require_immersed_operator_compatible(
      const mesh::MeshTopology &) const;
  double value_for_immersed_operator(mesh::LocalFaceId local_face) const;
  runtime::FieldId field_{};
  std::size_t face_count_{};
  std::unique_ptr<OwnedState> owned_state_;
  State *state_{};
  PreparedState *prepared_{};
  friend class CellCenteredFvmOperators;
  friend class ImmersedOperatorAdapter;
  friend class ImmersedReconstruction;
  friend class flow::FixedStepMaterialDensityFlow;
  friend class flow::MaterialFaceMassFlux;
  friend class flow::PisoCoupler;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  friend class test::PreparedFaceMassFluxForTest;
#endif
};

class CellCenteredFvmOperators final {
public:
  static CellCenteredFvmOperators create(const mesh::MeshTopology &topology,
                                         const mesh::MeshGeometry &geometry);

  ~CellCenteredFvmOperators() noexcept;
  CellCenteredFvmOperators(CellCenteredFvmOperators &&) noexcept;
  CellCenteredFvmOperators &operator=(CellCenteredFvmOperators &&) = delete;
  CellCenteredFvmOperators(const CellCenteredFvmOperators &) = delete;
  CellCenteredFvmOperators &
  operator=(const CellCenteredFvmOperators &) = delete;

  void compute_gradient(GradientScheme scheme, FiniteVolumeQuantity quantity,
                        const boundary::BoundaryRegistry &boundaries,
                        const runtime::FieldView<const double> &cell_values,
                        const runtime::FieldView<double> &cell_gradients) const;

  void reconstruct_transport_faces(
      FiniteVolumeQuantity quantity,
      const boundary::BoundaryRegistry &boundaries,
      const FaceMassFlux &mass_flux,
      const runtime::FieldView<const double> &cell_values,
      const runtime::FaceFieldView<double> &face_values) const;

  void reconstruct_momentum_faces(
      const boundary::BoundaryRegistry &boundaries,
      const FaceMassFlux &mass_flux,
      const runtime::FieldView<const double> &velocity,
      const runtime::FaceFieldView<double> &face_velocity) const;

  void assemble_provisional_mass_flux(
      const boundary::BoundaryRegistry &boundaries,
      const runtime::FieldView<const double> &density,
      const runtime::FieldView<const double> &velocity,
      const runtime::FieldRegistry &registry, runtime::FieldStorage &storage,
      const runtime::FieldAccessPlan &access_plan, runtime::PhaseId phase,
      runtime::ActorId actor, runtime::FieldId mass_flux_field) const;

  void accumulate_mass_residual(
      const FaceMassFlux &mass_flux,
      const runtime::FieldView<double> &raw_residual) const;

  void accumulate_convective_residual(
      const FaceMassFlux &mass_flux,
      const runtime::FaceFieldView<const double> &transported_face_values,
      const runtime::FieldView<double> &raw_residual) const;

  void accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity quantity,
      const boundary::BoundaryRegistry &boundaries,
      const runtime::FieldView<const double> &cell_values,
      const runtime::FieldView<const double> &cell_gradients,
      const runtime::FaceFieldView<const double> &gamma_by_face,
      const runtime::FieldView<double> &raw_residual) const;

  void accumulate_viscous_residual(
      const boundary::BoundaryRegistry &boundaries,
      const runtime::FieldView<const double> &velocity,
      const runtime::FieldView<const double> &velocity_gradients,
      double dynamic_viscosity_pa_s,
      const runtime::FieldView<double> &raw_momentum_residual) const;

  // `output` is cleared and refilled; callers may retain its capacity for
  // repeated, non-overlapping calls on this non-concurrent operator object.
  void physical_boundary_momentum_contributions(
      const boundary::BoundaryRegistry &boundaries,
      const FaceMassFlux &mass_flux,
      const runtime::FaceFieldView<const double> &face_velocity,
      const runtime::FieldView<const double> &velocity,
      const runtime::FieldView<const double> &velocity_gradients,
      double dynamic_viscosity_pa_s,
      std::vector<PhysicalBoundaryMomentumContribution> &output) const;

  void physical_boundary_pressure_contributions(
      const boundary::BoundaryRegistry &boundaries,
      const runtime::FieldView<const double> &pressure,
      std::vector<PhysicalBoundaryPressureContribution> &output) const;

  void physical_boundary_transport_contributions(
      FiniteVolumeQuantity quantity,
      const boundary::BoundaryRegistry &boundaries,
      const FaceMassFlux &mass_flux,
      const runtime::FaceFieldView<const double> &face_values,
      const runtime::FieldView<const double> &cell_values,
      const runtime::FieldView<const double> &cell_gradients,
      const runtime::FaceFieldView<const double> &gamma_by_face,
      std::vector<PhysicalBoundaryTransportContribution> &output) const;

private:
  struct Impl;
  explicit CellCenteredFvmOperators(std::unique_ptr<Impl>) noexcept;
  Impl &require_impl() const;
  std::unique_ptr<Impl> impl_;
};

} // namespace hundun::finite_volume
