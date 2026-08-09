// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flow_density_closure_detail.hpp"

#include "hundun/flow_immersed.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_quadratic_reconstruction.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace hundun::flow::detail {

struct ImmersedWallDensityValue final {
  immersed::ImmersedLinkId link{};
  double rho_wall_kg_per_m3{};
  double normal_derivative_kg_per_m4{};
};

struct ImmersedDensityAuthorityFailure final {
  MaterialTransportFailureReason reason{
      MaterialTransportFailureReason::non_finite_state};
};

inline ImmersedWallDensityValue reconstruct_immersed_wall_density(
    const immersed::ImmersedLink &link,
    const immersed::QuadraticReconstruction &row,
    const runtime::FieldView<const double> &density) {
  const double rho_wall = row.value(link.wall_intercept_m, density, 0U);
  const auto gradient = row.gradient(link.wall_intercept_m, density, 0U);
  const double normal_derivative =
      gradient.x * link.solid_to_fluid_normal.x +
      gradient.y * link.solid_to_fluid_normal.y +
      gradient.z * link.solid_to_fluid_normal.z;
  if (!std::isfinite(rho_wall) || !std::isfinite(normal_derivative))
    throw ImmersedDensityAuthorityFailure{
        MaterialTransportFailureReason::non_finite_state};
  if (!(rho_wall > 0.0))
    throw ImmersedDensityAuthorityFailure{
        MaterialTransportFailureReason::non_positive_density};
  return {link.id, rho_wall, normal_derivative};
}

struct ImmersedDensityAttemptAuthority final {
  std::vector<double> owned_active_density;
  std::vector<double> face_density;
  std::vector<ImmersedWallDensityValue> wall_density;
  std::uint64_t fingerprint{};
};

struct ImmersedConstantDensityAdapter final {};

class ImmersedMaterialDensityAdapter final {
public:
  static ImmersedMaterialDensityAdapter create(
      const runtime::FieldRegistry &registry,
      const runtime::StructuredDecomposition &decomposition,
      const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
      const boundary::BoundaryRegistry &boundaries,
      const immersed::ImmersedDomain &domain, const runtime::MpiContext &mpi,
      runtime::HaloExchange &halo, FlowFieldIds fields,
      MaterialDensityTransportSpec specification) {
    const std::uint64_t field_count =
        1U + static_cast<std::uint64_t>(specification.scalar_densities.size());
    auto transport = DensityClosureBridge::create_immersed_transport(
        registry, decomposition, topology, geometry, boundaries, domain, mpi,
        halo, fields, std::move(specification));
    return ImmersedMaterialDensityAdapter(
        registry, topology, std::move(fields), field_count,
        std::move(transport));
  }

  ImmersedMaterialDensityAdapter(ImmersedMaterialDensityAdapter &&) noexcept =
      default;
  ImmersedMaterialDensityAdapter &
  operator=(ImmersedMaterialDensityAdapter &&) = delete;
  ImmersedMaterialDensityAdapter(const ImmersedMaterialDensityAdapter &) =
      delete;
  ImmersedMaterialDensityAdapter &
  operator=(const ImmersedMaterialDensityAdapter &) = delete;

  void prepare_attempt() {
    if (attempt_identity_ == std::numeric_limits<std::uint64_t>::max())
      throw runtime::Error("immersed material attempt identity would wrap");
    ++attempt_identity_;
    DensityClosureBridge::prepare_material_transport(transport_);
    prepared_ = true;
    final_report_.reset();
    post_closure_report_.reset();
    failure_reason_ = MaterialTransportFailureReason::none;
  }

  MaterialDensityStageResult stage_predictor_transport(
      FlowState &state, const MomentumTimeStencil &stencil) {
    return stage(state, stencil, MaterialFluxProvenance::predictor);
  }

  MaterialDensityStageResult stage_after_corrector_one(
      FlowState &state, const MomentumTimeStencil &stencil) {
    return stage(state, stencil, MaterialFluxProvenance::provisional);
  }

  MaterialDensityStageResult finalize_from_corrector_two_flux(
      FlowState &state, const MomentumTimeStencil &stencil) {
    if (!prepared_)
      throw runtime::Error("immersed material adapter is not prepared");
    auto flux = MaterialFaceMassFlux::acquire(
        *registry_, state.trial_layer(), access_, kPhase, kActor,
        fields_.face_mass_flux, *topology_,
        MaterialFluxProvenance::final_corrected);
    auto report = DensityClosureBridge::finalize_material_transport(
        transport_, state, flux, stencil);
    DensityClosureBridge::rebind_material_report_attempt_identity(
        report, attempt_identity_);
    failure_reason_ = report.reason();
    final_report_.emplace(std::move(report));
    return {failure_reason_, final_report_->lowest_failing_rank()};
  }

  MaterialDensityStageResult assess_after_closure(
      FlowState &state, const MomentumTimeStencil &stencil,
      double enthalpy_rate_J_per_kg_s) {
    if (!prepared_ || !final_report_)
      throw runtime::Error(
          "immersed material post-closure adapter is not prepared");
    auto flux = MaterialFaceMassFlux::acquire(
        *registry_, state.trial_layer(), access_, kPhase, kActor,
        fields_.face_mass_flux, *topology_,
        MaterialFluxProvenance::final_corrected);
    auto report = DensityClosureBridge::assess_material_after_closure(
        transport_, state, flux, stencil, enthalpy_rate_J_per_kg_s);
    DensityClosureBridge::rebind_material_report_attempt_identity(
        report, attempt_identity_);
    failure_reason_ = report.reason();
    post_closure_report_.emplace(std::move(report));
    return {failure_reason_, post_closure_report_->lowest_failing_rank()};
  }

  void prepare_commit() const {
    if (!prepared_ || !final_report_.has_value() ||
        final_report_->disposition() !=
            MaterialTransportDisposition::finalized ||
        final_report_->reason() != MaterialTransportFailureReason::none ||
        final_report_->attempt_identity() != attempt_identity_ ||
        final_report_->flux_provenance() !=
            MaterialFluxProvenance::final_corrected)
      throw runtime::Error("immersed material adapter is not final");
  }

  void publish() noexcept { prepared_ = false; }

  void rollback() noexcept {
    prepared_ = false;
    final_report_.reset();
    post_closure_report_.reset();
  }

  std::optional<MaterialDensityTransportReport> take_final_report() {
    return final_report_;
  }

  std::optional<MaterialDensityTransportReport> take_post_closure_report() {
    return post_closure_report_;
  }

  MaterialTransportFailureReason failure_reason() const noexcept {
    return failure_reason_;
  }
  std::uint64_t material_field_count() const noexcept { return field_count_; }
  runtime::FieldId shared_face_mass_flux_field() const noexcept {
    return fields_.face_mass_flux;
  }
  std::uint64_t attempt_identity() const noexcept { return attempt_identity_; }

private:
  static constexpr runtime::PhaseId kPhase = 0x53334431U;
  static constexpr runtime::ActorId kActor = 0x53334431U;

  ImmersedMaterialDensityAdapter(
      const runtime::FieldRegistry &registry,
      const mesh::MeshTopology &topology, FlowFieldIds fields,
      std::uint64_t field_count, MaterialDensityTransport transport)
      : registry_(&registry), topology_(&topology), fields_(std::move(fields)),
        field_count_(field_count), access_(registry),
        transport_(std::move(transport)) {
    access_.declare_access(kPhase, kActor, fields_.face_mass_flux,
                           runtime::AccessMode::read);
    access_.freeze();
  }

  MaterialDensityStageResult stage(FlowState &state,
                                   const MomentumTimeStencil &stencil,
                                   MaterialFluxProvenance provenance) {
    if (!prepared_)
      throw runtime::Error("immersed material adapter is not prepared");
    auto flux = MaterialFaceMassFlux::acquire(
        *registry_, state.trial_layer(), access_, kPhase, kActor,
        fields_.face_mass_flux, *topology_, provenance);
    const auto result = DensityClosureBridge::stage_material_transport(
        transport_, state, flux, stencil);
    failure_reason_ = result.reason;
    return result;
  }

  const runtime::FieldRegistry *registry_{};
  const mesh::MeshTopology *topology_{};
  FlowFieldIds fields_;
  std::uint64_t field_count_{};
  runtime::FieldAccessPlan access_;
  MaterialDensityTransport transport_;
  std::optional<MaterialDensityTransportReport> final_report_;
  std::optional<MaterialDensityTransportReport> post_closure_report_;
  MaterialTransportFailureReason failure_reason_{
      MaterialTransportFailureReason::none};
  std::uint64_t attempt_identity_{};
  bool prepared_{};
};

class ImmersedIdealGasDensityAdapter final {
public:
  static ImmersedIdealGasDensityAdapter create(
      const runtime::FieldRegistry &registry,
      const runtime::StructuredDecomposition &decomposition,
      const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
      const boundary::BoundaryRegistry &boundaries,
      const immersed::ImmersedDomain &domain, const runtime::MpiContext &mpi,
      runtime::HaloExchange &halo, FlowFieldIds fields,
      MaterialDensityTransportSpec specification,
      IdealGasClosure &closure) {
    if (!DensityClosureAdapter::matches_immersed(
            closure, topology, geometry, boundaries, domain, mpi, registry,
            fields))
      throw runtime::Error(
          "immersed ideal-gas closure collaborator does not match");
    const auto enthalpy_density = specification.enthalpy_density;
    auto material = ImmersedMaterialDensityAdapter::create(
        registry, decomposition, topology, geometry, boundaries, domain, mpi,
        halo, fields, std::move(specification));
    auto hooks = DensityClosureAdapter::bind(closure, enthalpy_density, 0.0);
    return ImmersedIdealGasDensityAdapter(std::move(material),
                                          std::move(hooks), std::move(fields),
                                          closure);
  }

  ImmersedIdealGasDensityAdapter(ImmersedIdealGasDensityAdapter &&) noexcept =
      default;
  ImmersedIdealGasDensityAdapter &
  operator=(ImmersedIdealGasDensityAdapter &&) = delete;
  ImmersedIdealGasDensityAdapter(const ImmersedIdealGasDensityAdapter &) =
      delete;
  ImmersedIdealGasDensityAdapter &
  operator=(const ImmersedIdealGasDensityAdapter &) = delete;

  void prepare_attempt() {
    material_.prepare_attempt();
    closure_report_.reset();
    closure_active_ = false;
  }

  void begin_closure(const FlowState &state) {
    hooks_.begin(hooks_.object, state, material_.attempt_identity());
    closure_active_ = true;
  }

  MaterialDensityStageResult stage_predictor_transport(
      FlowState &state, const MomentumTimeStencil &stencil) {
    return material_.stage_predictor_transport(state, stencil);
  }

  MaterialDensityStageResult stage_after_corrector_one(
      FlowState &state, const MomentumTimeStencil &stencil) {
    return material_.stage_after_corrector_one(state, stencil);
  }

  MaterialDensityStageResult finalize_from_corrector_two_flux(
      FlowState &state, const MomentumTimeStencil &stencil) {
    return material_.finalize_from_corrector_two_flux(state, stencil);
  }

  DensityClosureEvaluation evaluate(FlowState &state,
                                    DensityClosureStage stage) {
    const auto result = hooks_.evaluate(hooks_.object, state, stage);
    closure_report_ = DensityClosureAdapter::latest_report(*closure_);
    return result;
  }

  void after_halo(DensityClosureStage stage) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (hooks_.after_halo != nullptr)
      hooks_.after_halo(hooks_.object, stage, density_field(),
                        hooks_.enthalpy_density);
#else
    static_cast<void>(stage);
#endif
  }

  MaterialDensityStageResult assess_after_closure(
      FlowState &state, const MomentumTimeStencil &stencil) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (hooks_.before_post_assessment != nullptr)
      hooks_.before_post_assessment(hooks_.object, state);
#endif
    return material_.assess_after_closure(
        state, stencil, hooks_.enthalpy_rate_J_per_kg_s);
  }

  void before_outlet(FlowState &state) {
    hooks_.before_outlet(hooks_.object, state);
  }

  int before_prepare(FlowState &state, AcceptedStepMetadata accepted) {
    return hooks_.before_prepare(hooks_.object, state, accepted);
  }

  int prepare_commit() {
    const auto post = material_.take_post_closure_report();
    if (!post || post->disposition() !=
                     MaterialTransportDisposition::finalized ||
        post->reason() != MaterialTransportFailureReason::none ||
        post->attempt_identity() != material_.attempt_identity() ||
        post->flux_provenance() !=
            MaterialFluxProvenance::final_corrected)
      throw runtime::Error(
          "immersed ideal-gas post-closure adapter is not final");
    material_.prepare_commit();
    return hooks_.prepare(hooks_.object);
  }

  void publish() noexcept {
    material_.publish();
    hooks_.publish(hooks_.object);
    closure_active_ = false;
  }

  void rollback() noexcept {
    material_.rollback();
    if (closure_active_)
      hooks_.rollback(hooks_.object);
    closure_active_ = false;
  }

  std::optional<MaterialDensityTransportReport> take_final_report() {
    return material_.take_final_report();
  }

  std::optional<MaterialDensityTransportReport> take_post_closure_report() {
    return material_.take_post_closure_report();
  }

  std::optional<IdealGasClosureReport> take_closure_report() const {
    return closure_report_;
  }

  std::uint64_t material_field_count() const noexcept {
    return material_.material_field_count();
  }
  runtime::FieldId shared_face_mass_flux_field() const noexcept {
    return material_.shared_face_mass_flux_field();
  }
  runtime::FieldId density_field() const noexcept {
    return fields_.density;
  }
  std::uint64_t attempt_identity() const noexcept {
    return material_.attempt_identity();
  }
  double cp_J_per_kg_K() const noexcept {
    return DensityClosureAdapter::cp_J_per_kg_K(*closure_);
  }
  double gas_constant_J_per_kg_K() const noexcept {
    return DensityClosureAdapter::gas_constant_J_per_kg_K(*closure_);
  }
  double configured_pressure_pa() const noexcept {
    return DensityClosureAdapter::configured_pressure_pa(*closure_);
  }
  IdealGasClosureState closure_state() const { return closure_->state(); }

private:
  ImmersedIdealGasDensityAdapter(ImmersedMaterialDensityAdapter material,
                                 DensityClosureHooks hooks, FlowFieldIds fields,
                                 IdealGasClosure &closure)
      : material_(std::move(material)), hooks_(std::move(hooks)),
        fields_(std::move(fields)), closure_(&closure) {}

  ImmersedMaterialDensityAdapter material_;
  DensityClosureHooks hooks_;
  FlowFieldIds fields_;
  IdealGasClosure *closure_{};
  std::optional<IdealGasClosureReport> closure_report_;
  bool closure_active_{};
};

using ImmersedDensityAdapter =
    std::variant<ImmersedConstantDensityAdapter,
                 ImmersedMaterialDensityAdapter,
                 ImmersedIdealGasDensityAdapter>;

} // namespace hundun::flow::detail
