// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kEnthalpyAssembly = 1431U;
constexpr std::uint32_t kEnthalpyNumerical = 1432U;
constexpr std::uint32_t kEnthalpyPoint = 1433U;
constexpr std::uint32_t kEnthalpyEndpointPlan = 1434U;
constexpr std::uint32_t kEnthalpyEndpointOperator = 1435U;
constexpr std::uint32_t kEnthalpyEndpointSolve = 1436U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

bool finite_bdf(BdfCoefficients bdf) noexcept {
  return detail::valid_bdf_coefficients(bdf);
}

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

RevisionToken state_revision(
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions) noexcept {
  std::uint64_t hash = kFnvOffset;
  const PrimitiveHistory histories[]{state.density, state.velocity,
                                     state.pressure_perturbation,
                                     state.enthalpy, state.temperature};
  for (const PrimitiveHistory& history : histories) {
    const ConstFieldView views[]{history.trial, history.accepted,
                                 history.previous};
    for (ConstFieldView view : views) {
      hash = hash_mix(hash, view.revision);
      hash = hash_mix(hash, view.storage_identity);
      hash = hash_mix(hash, view.revision_domain);
    }
  }
  std::uint64_t pressure_bits = 0U;
  std::memcpy(&pressure_bits, &state.pressure_reference,
              sizeof(pressure_bits));
  hash = hash_mix(hash, pressure_bits);
  std::memcpy(&pressure_bits, &state.accepted_pressure_reference,
              sizeof(pressure_bits));
  hash = hash_mix(hash, pressure_bits);
  std::memcpy(&pressure_bits, &state.previous_pressure_reference,
              sizeof(pressure_bits));
  hash = hash_mix(hash, pressure_bits);
  const ConstFieldView material_fields[]{material.thermal_conductivity,
                                         material.effective_viscosity,
                                         material.enthalpy_diffusivity,
                                         velocity_gradient};
  for (ConstFieldView field : material_fields) {
    hash = hash_mix(hash, field.revision);
    hash = hash_mix(hash, field.storage_identity);
    hash = hash_mix(hash, field.revision_domain);
  }
  for (std::size_t index = 0U; index < contributions.size; ++index) {
    const ConstFieldView explicit_source =
        contributions.data[index].explicit_source_density;
    hash = hash_mix(hash, explicit_source.revision);
    hash = hash_mix(hash, explicit_source.storage_identity);
    hash = hash_mix(hash, explicit_source.revision_domain);
    if (contributions.data[index].has_implicit_sink) {
      const ConstFieldView implicit_sink =
          contributions.data[index].implicit_sink_density;
      hash = hash_mix(hash, implicit_sink.revision);
      hash = hash_mix(hash, implicit_sink.storage_identity);
      hash = hash_mix(hash, implicit_sink.revision_domain);
    }
  }
  return hash == 0U ? RevisionToken{1U} : hash;
}

bool compatible_history(PrimitiveHistory history, Int3 cells, FieldId field,
                        std::uint8_t components,
                        std::uint8_t trial_ghosts) noexcept {
  return history.trial.field == field && history.accepted.field == field &&
         history.previous.field == field &&
         detail::valid_cell_view(history.trial, cells, 0U, components,
                                 trial_ghosts) &&
         detail::valid_cell_view(history.accepted, cells, 0U, components,
                                 0U) &&
         detail::valid_cell_view(history.previous, cells, 0U, components,
                                 0U);
}

KernelBox resolved_box(KernelBox requested, Int3 cells) noexcept {
  if (requested.begin.x == 0 && requested.begin.y == 0 &&
      requested.begin.z == 0 && requested.cells.x == 0 &&
      requested.cells.y == 0 && requested.cells.z == 0) {
    return {{0, 0, 0}, cells};
  }
  return requested;
}

bool valid_flux_context(const CartesianKernelPlan& kernels,
                        const EquationAssemblyContext& context,
                        KernelBox box) noexcept {
  if (!detail::bdf_matches_time_step(context.bdf, context.dt) ||
      context.time == 0U ||
      context.geometry == 0U ||
      context.face_flux == 0U ||
      context.face_flux != context.mass_flux.revision ||
      !detail::valid_kernel_box(box, kernels.cells()) ||
      !detail::valid_flux_view(context.mass_flux, kernels.cells(),
                               context.face_flux)) {
    return false;
  }
  if (context.scope == EquationAssemblyScope::final_conservative) {
    return !context.provisional_mass_flux &&
           context.face_flux_authority != 0U &&
           context.face_flux_storage != 0U &&
           context.face_flux_revision_domain != 0U &&
           context.face_flux_authority ==
               context.mass_flux.certificate.authority() &&
           context.face_flux_storage ==
               context.mass_flux.certificate.storage() &&
           context.face_flux_revision_domain ==
               context.mass_flux.certificate.revision_domain() &&
           context.mass_flux.certificate.matches(context.mass_flux);
  }
  if (context.scope == EquationAssemblyScope::target_coupled) {
    return !context.provisional_mass_flux &&
           context.face_flux_authority == 0U &&
           context.face_flux_storage == 0U &&
           context.face_flux_revision_domain == 0U &&
           !context.mass_flux.certificate.valid();
  }
  return context.scope == EquationAssemblyScope::momentum_predictor &&
         context.provisional_mass_flux &&
         !context.mass_flux.certificate.valid();
}

bool valid_linear_system(EquationSystemView system, Int3 cells,
                         bool& linear) noexcept {
  linear = system.diagonal.base != nullptr || system.rhs.base != nullptr ||
           system.x_coefficient.base != nullptr ||
           system.y_coefficient.base != nullptr ||
           system.z_coefficient.base != nullptr;
  if (!detail::valid_cell_view(system.residual, cells, 0U, 1U)) {
    return false;
  }
  if (!linear) {
    return true;
  }
  return detail::valid_cell_view(system.diagonal, cells, 0U, 1U) &&
         detail::valid_cell_view(system.rhs, cells, 0U, 1U) &&
         detail::valid_equation_faces(system, cells) &&
         detail::equation_face_views_disjoint(system) &&
         !detail::field_views_overlap(as_const(system.diagonal),
                                      as_const(system.rhs)) &&
         !detail::field_views_overlap(as_const(system.diagonal),
                                      as_const(system.residual)) &&
         !detail::field_views_overlap(as_const(system.rhs),
                                      as_const(system.residual));
}

double pressure_gradient_component(const CartesianKernelPlan& kernels,
                                   ConstFieldView pressure, Int3 cell,
                                   std::size_t axis) noexcept {
  const std::int32_t normal =
      axis == 0U ? cell.x : (axis == 1U ? cell.y : cell.z);
  const detail::DerivativeWeights weights =
      kernels.geometry_kind() == GeometryKind::uniform
          ? detail::metric_derivative_weights<true>(kernels, axis, normal)
          : detail::metric_derivative_weights<false>(kernels, axis, normal);
  Int3 minus = cell;
  Int3 plus = cell;
  if (axis == 0U) {
    --minus.x;
    ++plus.x;
  } else if (axis == 1U) {
    --minus.y;
    ++plus.y;
  } else {
    --minus.z;
    ++plus.z;
  }
  return weights.minus * pressure.unchecked(minus, 0U) +
         weights.centre * pressure.unchecked(cell, 0U) +
         weights.plus * pressure.unchecked(plus, 0U);
}

bool valid_contributions(Span<const EquationContributionView> contributions,
                         Span<const CompiledContribution> descriptors,
                         StageId stage, Int3 cells) noexcept {
  const std::size_t expected = descriptors.size;
  if (contributions.size != expected ||
      (contributions.size != 0U && contributions.data == nullptr)) {
    return false;
  }
  for (std::size_t i = 0U; i < contributions.size; ++i) {
    const EquationContributionView view = contributions.data[i];
    if (!detail::valid_cell_view(view.explicit_source_density, cells, 0U, 1U,
                                 0U) ||
        view.stage != stage || view.stage != descriptors.data[i].stage ||
        view.conserved_quantity != descriptors.data[i].conserved_quantity ||
        !(view.units == descriptors.data[i].units) ||
        view.explicit_source_field != descriptors.data[i].explicit_source ||
        view.explicit_source_density.field !=
            descriptors.data[i].explicit_source ||
        view.has_implicit_sink !=
            descriptors.data[i].supplies_implicit_diagonal ||
        (view.has_implicit_sink &&
         (!detail::valid_cell_view(view.implicit_sink_density, cells, 0U, 1U,
                                   0U) ||
          view.implicit_sink_field != descriptors.data[i].implicit_diagonal ||
          view.implicit_sink_density.field !=
              descriptors.data[i].implicit_diagonal))) {
      return false;
    }
  }
  return true;
}

bool output_aliases(ConstFieldView input, EquationSystemView system) noexcept {
  return detail::output_aliases_input(input, system, true);
}

std::size_t endpoint_cell_offset(Int3 cells, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(cell.z));
}

bool valid_endpoint_activity(MgDomainActivityView activity,
                             Int3 cells) noexcept {
  const Int3 x_faces{cells.x + 1, cells.y, cells.z};
  const Int3 y_faces{cells.x, cells.y + 1, cells.z};
  const Int3 z_faces{cells.x, cells.y, cells.z + 1};
  const auto count = [](Int3 shape) noexcept {
    return static_cast<std::size_t>(shape.x) *
           static_cast<std::size_t>(shape.y) *
           static_cast<std::size_t>(shape.z);
  };
  const bool empty = activity.cells.size == 0U &&
                     activity.x_faces.size == 0U &&
                     activity.y_faces.size == 0U &&
                     activity.z_faces.size == 0U;
  if (empty)
    return activity.local_fingerprint == 0U &&
           activity.collective_fingerprint == 0U;
  return activity.cells.data != nullptr && activity.x_faces.data != nullptr &&
         activity.y_faces.data != nullptr && activity.z_faces.data != nullptr &&
         activity.cells.size == count(cells) &&
         activity.x_faces.size == count(x_faces) &&
         activity.y_faces.size == count(y_faces) &&
         activity.z_faces.size == count(z_faces) &&
         activity.local_fingerprint != 0U &&
         activity.collective_fingerprint != 0U;
}

bool endpoint_fluid_cell(MgDomainActivityView activity, Int3 cells,
                         Int3 cell) noexcept {
  return activity.cells.size == 0U ||
         activity.cells.data[endpoint_cell_offset(cells, cell)] ==
             static_cast<std::uint8_t>(RegionFlag::fluid);
}

bool endpoint_active_face(MgDomainActivityView activity, Int3 cells,
                          CartesianAxis axis, Int3 face) noexcept {
  if (activity.cells.size == 0U) return true;
  Int3 extents = cells;
  Span<const std::uint8_t> faces;
  if (axis == CartesianAxis::x) {
    ++extents.x;
    faces = activity.x_faces;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
    faces = activity.y_faces;
  } else {
    ++extents.z;
    faces = activity.z_faces;
  }
  return faces.data[endpoint_cell_offset(extents, face)] != 0U;
}

std::size_t endpoint_face_index(CartesianFace face) noexcept {
  return static_cast<std::size_t>(face);
}

CartesianFace endpoint_cartesian_face(CartesianAxis axis,
                                      bool high) noexcept {
  return axis == CartesianAxis::x
             ? (high ? CartesianFace::x_max : CartesianFace::x_min)
             : (axis == CartesianAxis::y
                    ? (high ? CartesianFace::y_max : CartesianFace::y_min)
                    : (high ? CartesianFace::z_max
                            : CartesianFace::z_min));
}

std::int32_t endpoint_axis_extent(Int3 cells,
                                  CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? cells.x
             : (axis == CartesianAxis::y ? cells.y : cells.z);
}

std::int32_t endpoint_axis_coordinate(Int3 cell,
                                      CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? cell.x
             : (axis == CartesianAxis::y ? cell.y : cell.z);
}

Int3 endpoint_shifted(Int3 cell, CartesianAxis axis,
                      int direction) noexcept {
  if (axis == CartesianAxis::x)
    cell.x += direction;
  else if (axis == CartesianAxis::y)
    cell.y += direction;
  else
    cell.z += direction;
  return cell;
}

ConstFaceFieldView endpoint_face_flux(ConstFaceFluxView flux,
                                      CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

double endpoint_homogeneous_neighbor(
    ConstFieldView field, Int3 cells,
    const std::array<BoundaryRelation, 6U>& relations,
    const std::array<bool, 6U>& physical, Int3 cell,
    CartesianAxis axis, int direction) noexcept {
  const Int3 neighbor = endpoint_shifted(cell, axis, direction);
  const std::int32_t coordinate =
      endpoint_axis_coordinate(neighbor, axis);
  if (coordinate >= 0 && coordinate < endpoint_axis_extent(cells, axis))
    return field.unchecked(neighbor, 0U);
  const CartesianFace selected =
      endpoint_cartesian_face(axis, direction > 0);
  const std::size_t index = endpoint_face_index(selected);
  if (!physical[index]) return field.unchecked(neighbor, 0U);
  const double centre = field.unchecked(cell, 0U);
  switch (relations[index]) {
    case BoundaryRelation::dirichlet:
    case BoundaryRelation::convective:
      return -centre;
    case BoundaryRelation::zero_gradient:
    case BoundaryRelation::normal_gradient:
    case BoundaryRelation::reflect_normal:
      return centre;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double endpoint_incoming_flux(ConstFaceFluxView flux, Int3 cell,
                              CartesianAxis axis, int direction) noexcept {
  const Int3 face = direction < 0 ? cell : endpoint_shifted(cell, axis, 1);
  const double mass = endpoint_face_flux(flux, axis).unchecked(face);
  return direction < 0 ? std::max(mass, 0.0) : std::max(-mass, 0.0);
}

double endpoint_outgoing_flux(ConstFaceFluxView flux, Int3 cell,
                              CartesianAxis axis, int direction) noexcept {
  const Int3 face = direction < 0 ? cell : endpoint_shifted(cell, axis, 1);
  const double mass = endpoint_face_flux(flux, axis).unchecked(face);
  return direction < 0 ? std::max(-mass, 0.0) : std::max(mass, 0.0);
}

double endpoint_neighbor_sum(
    ConstFieldView field, ConstFaceFluxView flux,
    MgDomainActivityView activity,
    const std::array<BoundaryRelation, 6U>& relations,
    const std::array<bool, 6U>& physical, const CartesianKernelPlan& kernels,
    Int3 cell) noexcept {
  const Int3 cells = kernels.cells();
  const double inverse_volume = 1.0 / detail::cell_volume(kernels, cell);
  double sum = 0.0;
  for (CartesianAxis axis : {CartesianAxis::x, CartesianAxis::y,
                             CartesianAxis::z}) {
    for (int direction : {-1, 1}) {
      const Int3 face = direction < 0
                            ? cell
                            : endpoint_shifted(cell, axis, 1);
      if (!endpoint_active_face(activity, cells, axis, face)) continue;
      sum += endpoint_incoming_flux(flux, cell, axis, direction) *
             inverse_volume *
             endpoint_homogeneous_neighbor(field, cells, relations,
                                            physical, cell, axis, direction);
    }
  }
  return sum;
}

Status resolve_endpoint_boundaries(
    const BoundaryPlan& boundary, FieldId enthalpy,
    std::array<BoundaryRelation, 6U>& relations,
    std::array<bool, 6U>& physical) noexcept {
  std::array<bool, 6U> found{};
  for (std::size_t index = 0U; index < physical.size(); ++index) {
    const BoundaryFacePlan* face = nullptr;
    const Status status =
        boundary.face(static_cast<CartesianFace>(index), face);
    if (!status || face == nullptr)
      return status ? Status{StatusCode::invalid_plan,
                             kEnthalpyEndpointPlan}
                    : status;
    physical[index] = face->local_owner && !face->periodic;
  }
  const Span<const BoundaryIndexSpan> spans = boundary.spans();
  for (std::size_t index = 0U; index < spans.size; ++index) {
    const BoundaryIndexSpan span = spans.data[index];
    if (span.stage != BoundaryStage::enthalpy || span.field != enthalpy ||
        span.component_begin != 0U || span.component_count != 1U)
      continue;
    const std::size_t selected = endpoint_face_index(span.face);
    if (selected >= relations.size() || found[selected])
      return {StatusCode::invalid_plan, kEnthalpyEndpointPlan};
    relations[selected] = span.relation;
    found[selected] = true;
  }
  for (std::size_t index = 0U; index < physical.size(); ++index) {
    if (physical[index] && !found[index])
      return {StatusCode::invalid_plan, kEnthalpyEndpointPlan};
  }
  return {};
}

class EnthalpyEndpointOperator final : public LinearOperator {
 public:
  EnthalpyEndpointOperator(
      LinearIdentity identity, PlanFingerprint collective_fingerprint,
      ConstFieldView diagonal, ConstFaceFluxView flux,
      MgDomainActivityView activity,
      const std::array<BoundaryRelation, 6U>& relations,
      const std::array<bool, 6U>& physical,
      const CartesianKernelPlan& kernels, FieldId workspace_field,
      HaloEngine& halo, StageId halo_stage) noexcept
      : diagonal_(diagonal),
        flux_(flux),
        activity_(activity),
        relations_(relations),
        physical_(physical),
        kernels_(&kernels),
        workspace_field_(workspace_field),
        halo_(&halo),
        halo_stage_(halo_stage) {
    certificate_.identity = identity;
    certificate_.collective_fingerprint = collective_fingerprint;
    certificate_.local_shape = diagonal.interior;
    certificate_.operator_class = LinearOperatorClass::nonsymmetric;
  }

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_;
  }

  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return failure_;
  }

  Status apply(FieldView x, FieldView y) const noexcept override {
    failure_ = {};
    const Int3 cells = diagonal_.interior;
    if (kernels_ == nullptr || halo_ == nullptr ||
        x.field != workspace_field_ || x.components != 1U ||
        y.components != 1U || x.base == nullptr || y.base == nullptr ||
        x.interior.x != cells.x || x.interior.y != cells.y ||
        x.interior.z != cells.z || y.interior.x != cells.x ||
        y.interior.y != cells.y || y.interior.z != cells.z ||
        x.ghosts.x < 1 || x.ghosts.y < 1 || x.ghosts.z < 1 ||
        detail::field_views_overlap(as_const(x), as_const(y))) {
      return {StatusCode::invalid_plan, kEnthalpyEndpointOperator};
    }
    std::array<FieldView, 1U> fields{x};
    HaloTicket ticket;
    Status status =
        halo_->begin(halo_stage_, {fields.data(), fields.size()}, ticket);
    if (status)
      status = halo_->finish(ticket, {fields.data(), fields.size()});
    if (!status) {
      failure_ = {status, LinearOperatorStatusScope::collective,
                  halo_->lowest_failing_rank()};
      return status;
    }
    x = fields[0U];
    Status arithmetic;
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y_index = 0; y_index < cells.y; ++y_index) {
        for (std::int32_t x_index = 0; x_index < cells.x; ++x_index) {
          const Int3 cell{x_index, y_index, z};
          const double diagonal = diagonal_.unchecked(cell, 0U);
          const double neighbor = endpoint_fluid_cell(activity_, cells, cell)
                                      ? endpoint_neighbor_sum(
                                            as_const(x), flux_, activity_,
                                            relations_, physical_, *kernels_,
                                            cell)
                                      : 0.0;
          const double value = diagonal * x.unchecked(cell, 0U) - neighbor;
          y.unchecked(cell, 0U) = value;
          if (!std::isfinite(diagonal) || !(diagonal > 0.0) ||
              !std::isfinite(neighbor) || !std::isfinite(value))
            arithmetic = {StatusCode::numerical_failure,
                          kEnthalpyEndpointOperator};
        }
      }
    }
    if (!arithmetic)
      failure_ = {arithmetic, LinearOperatorStatusScope::rank_local, -1};
    return arithmetic;
  }

 private:
  ConstFieldView diagonal_{};
  ConstFaceFluxView flux_{};
  MgDomainActivityView activity_{};
  std::array<BoundaryRelation, 6U> relations_{};
  std::array<bool, 6U> physical_{};
  const CartesianKernelPlan* kernels_{};
  FieldId workspace_field_{};
  HaloEngine* halo_{};
  StageId halo_stage_{};
  LinearOperatorCertificate certificate_{};
  mutable LinearOperatorFailureProvenance failure_{};
};

class EnthalpyEndpointJacobi final : public LinearPreconditioner {
 public:
  EnthalpyEndpointJacobi(LinearIdentity identity,
                         PlanFingerprint collective_fingerprint,
                         ConstFieldView diagonal) noexcept
      : diagonal_(diagonal) {
    certificate_.identity = identity;
    certificate_.collective_fingerprint = collective_fingerprint;
    certificate_.preconditioner_class =
        LinearPreconditionerClass::fixed_general;
    certificate_.status_scope = LinearPreconditionerStatusScope::rank_local;
    certificate_.apply_lifecycle =
        LinearPreconditionerApplyLifecycle::per_call_checked;
  }

  LinearPreconditionerCertificate certificate() const noexcept override {
    return certificate_;
  }

  Status apply(ConstFieldView input, FieldView output,
               std::uint32_t) noexcept override {
    const Int3 cells = diagonal_.interior;
    if (input.base == nullptr || output.base == nullptr ||
        input.components != 1U || output.components != 1U ||
        input.interior.x != cells.x || input.interior.y != cells.y ||
        input.interior.z != cells.z || output.interior.x != cells.x ||
        output.interior.y != cells.y || output.interior.z != cells.z ||
        detail::field_views_overlap(input, as_const(output)))
      return {StatusCode::invalid_plan, kEnthalpyEndpointOperator};
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const double diagonal = diagonal_.unchecked(cell, 0U);
          const double value = input.unchecked(cell, 0U) / diagonal;
          if (!std::isfinite(diagonal) || !(diagonal > 0.0) ||
              !std::isfinite(value))
            return {StatusCode::numerical_failure,
                    kEnthalpyEndpointOperator};
          output.unchecked(cell, 0U) = value;
        }
      }
    }
    return {};
  }

 private:
  ConstFieldView diagonal_{};
  LinearPreconditionerCertificate certificate_{};
};

}  // namespace

Status ConservativeEnthalpyEndpoint::bind(
    const ConservativeEnthalpyEndpointServices& services,
    ConservativeEnthalpyEndpoint& out) noexcept {
  if (out.fingerprint_ != 0U || services.communicator == MPI_COMM_NULL ||
      services.kernels == nullptr || services.boundary == nullptr ||
      services.halo == nullptr || services.reductions == nullptr ||
      services.halo_stage == 0U ||
      services.patch.cells.x <= 0 || services.patch.cells.y <= 0 ||
      services.patch.cells.z <= 0 ||
      services.kernels->cells().x != services.patch.cells.x ||
      services.kernels->cells().y != services.patch.cells.y ||
      services.kernels->cells().z != services.patch.cells.z ||
      services.boundary->local_cells().x != services.patch.cells.x ||
      services.boundary->local_cells().y != services.patch.cells.y ||
      services.boundary->local_cells().z != services.patch.cells.z ||
      services.workspace_requirements.algorithm != LinearAlgorithm::fgmres ||
      services.workspace_requirements.ghost_width < 1U ||
      services.workspace_requirements.maximum_restart == 0U ||
      !valid_endpoint_activity(services.activity, services.patch.cells)) {
    return {StatusCode::invalid_plan, kEnthalpyEndpointPlan};
  }
  Status status = SolverWorkspace::bind(
      services.workspace_requirements, services.krylov_vectors,
      services.krylov_scalars, out.workspace_);
  if (!status) return status;
  const FieldView probe = out.workspace_.vector(0U, services.patch.cells);
  const std::array<HaloFieldSpec, 1U> halo_fields{{
      {probe.field, 1U, 1U},
  }};
  if (probe.base == nullptr ||
      !services.halo->validate_contract(
          services.communicator, services.patch,
          {halo_fields.data(), halo_fields.size()},
          services.boundary->halo_topology())) {
    return {StatusCode::invalid_plan, kEnthalpyEndpointPlan};
  }
  status = resolve_endpoint_boundaries(
      *services.boundary, services.enthalpy, out.boundary_relations_,
      out.physical_boundaries_);
  if (!status) return status;
  status = services.reductions->validate_communicator(services.communicator);
  if (!status) return status;

  out.communicator_ = services.communicator;
  out.kernels_ = services.kernels;
  out.boundary_ = services.boundary;
  out.patch_ = services.patch;
  out.activity_ = services.activity;
  out.halo_ = services.halo;
  out.reductions_ = services.reductions;
  out.halo_stage_ = services.halo_stage;
  out.enthalpy_ = services.enthalpy;
  std::uint64_t collective =
      hash_mix(kFnvOffset, UINT64_C(0x524f55544542));
  collective =
      hash_mix(collective, services.boundary->semantic_fingerprint());
  collective =
      hash_mix(collective, services.activity.collective_fingerprint);
  collective = hash_mix(collective, services.enthalpy);
  out.collective_fingerprint_ = collective == 0U ? 1U : collective;
  std::uint64_t fingerprint =
      hash_mix(out.collective_fingerprint_, services.kernels->fingerprint());
  fingerprint = hash_mix(fingerprint, out.workspace_.fingerprint());
  out.fingerprint_ = fingerprint == 0U ? 1U : fingerprint;
  return {};
}

Status ConservativeEnthalpyEndpoint::solve(
    const ConservativeEnthalpyEndpointInput& input,
    LinearSolveResult& result) noexcept {
  result = {};
  const Int3 cells = patch_.cells;
  const auto valid_scalar = [&](ConstFieldView view,
                                std::uint8_t ghosts) noexcept {
    return detail::valid_cell_view(view, cells, 0U, 1U, ghosts);
  };
  Status local;
  if (fingerprint_ == 0U || communicator_ == MPI_COMM_NULL ||
      kernels_ == nullptr || boundary_ == nullptr || halo_ == nullptr ||
      reductions_ == nullptr || !detail::valid_bdf_coefficients(input.bdf) ||
      input.time == 0U || !valid_scalar(input.predicted_density, 0U) ||
      input.accepted_enthalpy.field != enthalpy_ ||
      !valid_scalar(input.accepted_enthalpy, 1U) ||
      input.high_enthalpy.field != enthalpy_ ||
      !valid_scalar(input.high_enthalpy, 0U) ||
      input.diagonal_workspace.components != 1U ||
      !detail::valid_cell_view(as_const(input.diagonal_workspace), cells, 0U,
                               1U, 0U) ||
      input.rhs_workspace.components != 1U ||
      !detail::valid_cell_view(as_const(input.rhs_workspace), cells, 0U, 1U,
                               0U) ||
      input.endpoint.field != enthalpy_ ||
      !detail::valid_cell_view(as_const(input.endpoint), cells, 0U, 1U, 1U) ||
      !detail::valid_flux_view(input.accepted_mass_flux, cells,
                               input.accepted_mass_flux.revision) ||
      !input.accepted_mass_flux.certificate.valid() ||
      !input.accepted_mass_flux.certificate.matches(
          input.accepted_mass_flux) ||
      detail::field_views_overlap(as_const(input.diagonal_workspace),
                                  as_const(input.rhs_workspace)) ||
      detail::field_views_overlap(as_const(input.endpoint),
                                  input.accepted_enthalpy) ||
      detail::field_views_overlap(as_const(input.endpoint),
                                  input.high_enthalpy) ||
      workspace_.overlaps_storage(input.diagonal_workspace) ||
      workspace_.overlaps_storage(input.rhs_workspace) ||
      workspace_.overlaps_storage(input.endpoint)) {
    local = {StatusCode::invalid_plan, kEnthalpyEndpointPlan};
  }
  Status status = reductions_->consensus(local);
  if (!status) return status;

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double accepted = input.accepted_enthalpy.unchecked(cell, 0U);
        input.endpoint.unchecked(cell, 0U) = accepted;
        if (!endpoint_fluid_cell(activity_, cells, cell)) {
          input.diagonal_workspace.unchecked(cell, 0U) = 1.0;
          input.rhs_workspace.unchecked(cell, 0U) = accepted;
          if (!std::isfinite(accepted))
            local = {StatusCode::numerical_failure,
                     kEnthalpyEndpointOperator};
          continue;
        }
        const double density =
            input.predicted_density.unchecked(cell, 0U);
        const double high = input.high_enthalpy.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*kernels_, cell);
        double outgoing = 0.0;
        double incoming = 0.0;
        for (CartesianAxis axis : {CartesianAxis::x, CartesianAxis::y,
                                   CartesianAxis::z}) {
          for (int direction : {-1, 1}) {
            const Int3 face = direction < 0
                                  ? cell
                                  : endpoint_shifted(cell, axis, 1);
            if (!endpoint_active_face(activity_, cells, axis, face)) continue;
            outgoing += endpoint_outgoing_flux(
                input.accepted_mass_flux, cell, axis, direction);
            incoming += endpoint_incoming_flux(
                input.accepted_mass_flux, cell, axis, direction);
          }
        }
        const double inverse_volume = 1.0 / volume;
        const double base_diagonal =
            input.bdf.a0 * density + outgoing * inverse_volume;
        const double diagonal =
            std::max(base_diagonal, incoming * inverse_volume);
        const double neighbor = endpoint_neighbor_sum(
            input.accepted_enthalpy, input.accepted_mass_flux, activity_,
            boundary_relations_, physical_boundaries_, *kernels_, cell);
        const double high_residual =
            input.bdf.a0 * density * (accepted - high);
        const double rhs = diagonal * accepted - neighbor - high_residual;
        input.diagonal_workspace.unchecked(cell, 0U) = diagonal;
        input.rhs_workspace.unchecked(cell, 0U) = rhs;
        if (!std::isfinite(density) || !(density > 0.0) ||
            !std::isfinite(accepted) || !std::isfinite(high) ||
            !std::isfinite(volume) || !(volume > 0.0) ||
            !std::isfinite(outgoing) || outgoing < 0.0 ||
            !std::isfinite(incoming) || incoming < 0.0 ||
            !std::isfinite(base_diagonal) || !(base_diagonal > 0.0) ||
            !std::isfinite(diagonal) || !(diagonal > 0.0) ||
            !std::isfinite(neighbor) || !std::isfinite(high_residual) ||
            !std::isfinite(rhs)) {
          local = {StatusCode::numerical_failure,
                   kEnthalpyEndpointOperator};
        }
      }
    }
  }
  status = reductions_->consensus(local);
  if (!status) return status;

  const std::uint64_t symbolic = collective_fingerprint_;
  std::uint64_t numeric = hash_mix(symbolic, input.time);
  numeric = hash_mix(numeric, input.bdf.order);
  numeric = hash_mix(numeric, double_bits(input.bdf.a0));
  numeric = hash_mix(numeric, double_bits(input.bdf.a1));
  numeric = hash_mix(numeric, double_bits(input.bdf.a2));
  numeric = hash_mix(numeric, input.predicted_density.revision);
  numeric = hash_mix(numeric, input.accepted_enthalpy.revision);
  numeric = hash_mix(numeric, input.high_enthalpy.revision);
  numeric = hash_mix(numeric, input.accepted_mass_flux.revision);
  numeric = hash_mix(numeric, input.diagonal_workspace.revision);
  numeric = hash_mix(numeric, activity_.local_fingerprint);
  numeric = numeric == 0U ? 1U : numeric;
  std::uint64_t hierarchy = hash_mix(symbolic, UINT64_C(1));
  hierarchy = hierarchy == 0U ? 1U : hierarchy;
  LinearIdentity identity{symbolic, numeric, hierarchy,
                          workspace_.fingerprint(), 0U};
  std::uint64_t identity_fingerprint =
      hash_mix(kFnvOffset, identity.symbolic);
  identity_fingerprint = hash_mix(identity_fingerprint, identity.numeric);
  identity_fingerprint = hash_mix(identity_fingerprint, identity.hierarchy);
  identity_fingerprint = hash_mix(identity_fingerprint, identity.workspace);
  identity.fingerprint =
      identity_fingerprint == 0U ? 1U : identity_fingerprint;
  const FieldView workspace_probe = workspace_.vector(0U, cells);
  EnthalpyEndpointOperator linear_operator(
      identity, symbolic, as_const(input.diagonal_workspace),
      input.accepted_mass_flux, activity_, boundary_relations_,
      physical_boundaries_, *kernels_, workspace_probe.field, *halo_,
      halo_stage_);
  EnthalpyEndpointJacobi preconditioner(
      identity, symbolic, as_const(input.diagonal_workspace));
  const std::uint32_t restart = std::min<std::uint32_t>(
      12U, workspace_.requirements().maximum_restart);
  const LinearSolveControl control{1.0e-10, 1.0e-6, 128U, 4U, restart};
  const LinearSolveInvocation invocation{
      as_const(input.rhs_workspace), input.endpoint, identity, control,
      nullptr};
  result = solve_fgmres(linear_operator, preconditioner, invocation,
                        workspace_, *reductions_, input.resources);
  const bool accepted =
      result.status &&
      (result.termination == LinearTermination::converged ||
       result.termination == LinearTermination::zero_rhs);
  if (!accepted) {
    if (result.status.code == StatusCode::ok ||
        result.status.code == StatusCode::rejected_step) {
      return {StatusCode::numerical_failure, kEnthalpyEndpointSolve};
    }
    return result.status;
  }
  return {};
}

Status evaluate_pressure_material_derivative(
    BdfCoefficients bdf, const PressureWorkPoint& point,
    double& material_derivative) noexcept {
  if (!finite_bdf(bdf) || !std::isfinite(point.pressure) ||
      !std::isfinite(point.accepted_pressure) ||
      !std::isfinite(point.previous_pressure) ||
      !std::isfinite(point.pressure_gradient.x) ||
      !std::isfinite(point.pressure_gradient.y) ||
      !std::isfinite(point.pressure_gradient.z) ||
      !std::isfinite(point.velocity.x) || !std::isfinite(point.velocity.y) ||
      !std::isfinite(point.velocity.z)) {
    return {StatusCode::numerical_failure, kEnthalpyPoint};
  }
  const double candidate =
      bdf.a0 * point.pressure + bdf.a1 * point.accepted_pressure +
      bdf.a2 * point.previous_pressure +
      point.velocity.x * point.pressure_gradient.x +
      point.velocity.y * point.pressure_gradient.y +
      point.velocity.z * point.pressure_gradient.z;
  if (!std::isfinite(candidate)) {
    return {StatusCode::numerical_failure, kEnthalpyPoint};
  }
  material_derivative = candidate;
  return {};
}

Status assemble_isolated_enthalpy_term(const EnthalpyTermPoint& term,
                                       double cell_volume,
                                       double& residual) noexcept {
  if (!std::isfinite(term.value) || !std::isfinite(cell_volume) ||
      cell_volume <= 0.0 ||
      static_cast<std::uint8_t>(term.kind) >
          static_cast<std::uint8_t>(EnthalpyTermKind::source)) {
    return {StatusCode::numerical_failure, kEnthalpyPoint};
  }
  double sign = 1.0;
  if (term.kind == EnthalpyTermKind::pressure_work ||
      term.kind == EnthalpyTermKind::conduction ||
      term.kind == EnthalpyTermKind::viscous_dissipation ||
      term.kind == EnthalpyTermKind::source) {
    sign = -1.0;
  }
  const double candidate = sign * term.value * cell_volume;
  if (!std::isfinite(candidate)) {
    return {StatusCode::numerical_failure, kEnthalpyPoint};
  }
  residual = candidate;
  return {};
}

Status newtonian_viscous_dissipation(const VelocityGradient& gradient,
                                     double effective_viscosity,
                                     double& dissipation) noexcept {
  if (!std::isfinite(effective_viscosity) || effective_viscosity < 0.0) {
    return {StatusCode::numerical_failure, kEnthalpyPoint};
  }
  for (double value : gradient.value) {
    if (!std::isfinite(value)) {
      return {StatusCode::numerical_failure, kEnthalpyPoint};
    }
  }
  const double trace =
      gradient.value[0U] + gradient.value[4U] + gradient.value[8U];
  long double value = 0.0L;
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      const double gij = gradient.value[3U * i + j];
      const double gji = gradient.value[3U * j + i];
      const double tau = effective_viscosity *
                         (gij + gji -
                          (i == j ? (2.0 / 3.0) * trace : 0.0));
      value += static_cast<long double>(tau) * gij;
    }
  }
  const double candidate = static_cast<double>(value);
  if (!std::isfinite(candidate)) {
    return {StatusCode::numerical_failure, kEnthalpyPoint};
  }
  dissipation = candidate;
  return {};
}

namespace {

struct EnthalpyCellTerms {
  double pressure_work{};
  double viscous_dissipation{};
  double source{};
  double sink{};
  double diffusion_diagonal{};
  double diagonal{};
};

struct EnthalpyCellSystem {
  double residual{};
  double diagonal{};
  double rhs{};
};

void form_enthalpy_linear_terms(
    const CartesianKernelPlan& kernels,
    const EquationMaterialView& material,
    const EquationAssemblyContext& context, Int3 cell, double rho,
    EnthalpyCellTerms& terms) noexcept {
  const double volume = detail::cell_volume(kernels, cell);
  terms.diffusion_diagonal =
      detail::diffusion_diagonal(kernels, material.enthalpy_diffusivity, cell);
  terms.diagonal = (context.bdf.a0 * rho + terms.sink) * volume +
                   terms.diffusion_diagonal;
}

Status evaluate_enthalpy_cell_terms(
    const CartesianKernelPlan& kernels, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, Int3 cell,
    bool validate_individual_sinks, EnthalpyCellTerms& terms) noexcept {
  const double rho = state.density.trial.unchecked(cell, 0U);
  const double rho_n = state.density.accepted.unchecked(cell, 0U);
  const double rho_nm1 = state.density.previous.unchecked(cell, 0U);
  const double h = state.enthalpy.trial.unchecked(cell, 0U);
  const double hn = state.enthalpy.accepted.unchecked(cell, 0U);
  const double hnm1 = state.enthalpy.previous.unchecked(cell, 0U);
  const double temperature = state.temperature.trial.unchecked(cell, 0U);
  const double conductivity =
      material.thermal_conductivity.unchecked(cell, 0U);
  const double enthalpy_diffusivity =
      material.enthalpy_diffusivity.unchecked(cell, 0U);
  const double viscosity = material.effective_viscosity.unchecked(cell, 0U);
  if (!std::isfinite(rho) || !std::isfinite(rho_n) ||
      !std::isfinite(rho_nm1) || rho <= 0.0 || rho_n <= 0.0 ||
      rho_nm1 <= 0.0 || !std::isfinite(h) || !std::isfinite(hn) ||
      !std::isfinite(hnm1) || !std::isfinite(temperature) ||
      !std::isfinite(conductivity) || conductivity <= 0.0 ||
      !std::isfinite(enthalpy_diffusivity) || enthalpy_diffusivity <= 0.0 ||
      !std::isfinite(viscosity) || viscosity < 0.0) {
    return {StatusCode::numerical_failure, kEnthalpyNumerical};
  }

  VelocityGradient local_gradient;
  for (std::uint8_t component = 0U; component < 9U; ++component) {
    local_gradient.value[component] =
        velocity_gradient.unchecked(cell, component);
  }
  const PressureWorkPoint work{
      state.pressure_reference +
          state.pressure_perturbation.trial.unchecked(cell, 0U),
      state.accepted_pressure_reference +
          state.pressure_perturbation.accepted.unchecked(cell, 0U),
      state.previous_pressure_reference +
          state.pressure_perturbation.previous.unchecked(cell, 0U),
      {pressure_gradient_component(kernels,
                                   state.pressure_perturbation.trial, cell,
                                   0U),
       pressure_gradient_component(kernels,
                                   state.pressure_perturbation.trial, cell,
                                   1U),
       pressure_gradient_component(kernels,
                                   state.pressure_perturbation.trial, cell,
                                   2U)},
      {state.velocity.trial.unchecked(cell, 0U),
       state.velocity.trial.unchecked(cell, 1U),
       state.velocity.trial.unchecked(cell, 2U)},
      0.0};
  EnthalpyCellTerms candidate;
  if (!evaluate_pressure_material_derivative(
          context.bdf, work, candidate.pressure_work) ||
      !newtonian_viscous_dissipation(
          local_gradient, viscosity, candidate.viscous_dissipation)) {
    return {StatusCode::numerical_failure, kEnthalpyNumerical};
  }
  for (std::size_t index = 0U; index < contributions.size; ++index) {
    candidate.source += contributions.data[index]
                            .explicit_source_density.unchecked(cell, 0U);
    if (contributions.data[index].has_implicit_sink) {
      const double value = contributions.data[index]
                               .implicit_sink_density.unchecked(cell, 0U);
      if (validate_individual_sinks &&
          (!std::isfinite(value) || value < 0.0)) {
        return {StatusCode::numerical_failure, kEnthalpyNumerical};
      }
      candidate.sink += value;
    }
  }
  form_enthalpy_linear_terms(kernels, material, context, cell, rho,
                             candidate);
  if (!std::isfinite(candidate.source) ||
      !std::isfinite(candidate.sink) || candidate.sink < 0.0 ||
      !std::isfinite(candidate.diagonal) || candidate.diagonal <= 0.0 ||
      !std::isfinite(candidate.pressure_work) ||
      !std::isfinite(candidate.viscous_dissipation)) {
    return {StatusCode::numerical_failure, kEnthalpyNumerical};
  }
  terms = candidate;
  return {};
}

Status combine_enthalpy_cell_system(
    const CartesianKernelPlan& kernels, const EquationStateView& state,
    const EquationAssemblyContext& context, Int3 cell, double convection,
    double diffusion, const EnthalpyCellTerms& terms,
    EnthalpyCellSystem& system) noexcept {
  const double rho = state.density.trial.unchecked(cell, 0U);
  const double h = state.enthalpy.trial.unchecked(cell, 0U);
  const double unsteady =
      context.bdf.a0 * rho * h +
      context.bdf.a1 * state.density.accepted.unchecked(cell, 0U) *
          state.enthalpy.accepted.unchecked(cell, 0U) +
      context.bdf.a2 * state.density.previous.unchecked(cell, 0U) *
          state.enthalpy.previous.unchecked(cell, 0U);
  const double volume = detail::cell_volume(kernels, cell);
  EnthalpyCellSystem candidate;
  candidate.residual =
      (unsteady + convection - terms.pressure_work - diffusion -
       terms.viscous_dissipation - terms.source + terms.sink * h) *
      volume;
  candidate.diagonal = terms.diagonal;
  candidate.rhs = candidate.diagonal * h - candidate.residual;
  if (!std::isfinite(candidate.residual) ||
      !std::isfinite(candidate.diagonal) || candidate.diagonal <= 0.0 ||
      !std::isfinite(terms.diffusion_diagonal) ||
      terms.diffusion_diagonal <= 0.0 || !std::isfinite(candidate.rhs)) {
    return {StatusCode::numerical_failure, kEnthalpyNumerical};
  }
  system = candidate;
  return {};
}

}  // namespace

Status assemble_enthalpy_impl(
    const EnthalpyEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate, bool allow_partial) noexcept {
  Span<const CompiledContribution> selected_descriptors{};
  if (!detail::select_contribution_stage(
          {plan.contributions_.data(), plan.contributions_.size()},
          context.contribution_stage, selected_descriptors)) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }
  if (plan.kernels_ == nullptr || plan.fingerprint_ == 0U ||
      context.geometry != plan.geometry_revision_ ||
      context.boundary != plan.boundary_revision_ ||
      context.thermo != plan.thermodynamics_fingerprint_ ||
      context.transport != plan.transport_fingerprint_ ||
      (context.scope != EquationAssemblyScope::final_conservative &&
       context.scope != EquationAssemblyScope::target_coupled) ||
      context.provisional_mass_flux ||
      !compatible_history(state.density, plan.cells_, plan.density_, 1U, 0U) ||
      !compatible_history(state.velocity, plan.cells_, plan.velocity_, 3U, 1U) ||
      !compatible_history(state.pressure_perturbation, plan.cells_,
                          plan.pressure_, 1U, 1U) ||
      !compatible_history(state.enthalpy, plan.cells_, plan.enthalpy_, 1U,
                          plan.convection_reach_) ||
      !compatible_history(state.temperature, plan.cells_, plan.temperature_,
                          1U, 1U) ||
      !detail::valid_cell_view(material.thermal_conductivity, plan.cells_, 0U,
                               1U, 1U) ||
      !detail::valid_cell_view(material.effective_viscosity, plan.cells_, 0U,
                               1U, 0U) ||
      !detail::valid_cell_view(material.enthalpy_diffusivity, plan.cells_, 0U,
                               1U, 1U) ||
      velocity_gradient.field != plan.velocity_gradient_ ||
      !detail::valid_cell_view(velocity_gradient, plan.cells_, 0U, 9U, 0U) ||
      !valid_contributions(contributions, selected_descriptors,
          context.contribution_stage, plan.cells_)) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }
  bool linear = false;
  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (!valid_flux_context(*plan.kernels_, context, box) ||
      (!allow_partial && !detail::full_equation_box(box, plan.cells_)) ||
      !valid_linear_system(system, plan.cells_, linear) || !linear ||
      !detail::finite_face_flux(context.mass_flux, box) ||
      !detail::finite_face_neighbour_slabs(
          state.enthalpy.trial, box, 0U, 1U, plan.convection_reach_) ||
      !detail::finite_face_neighbour_slabs(
          state.temperature.trial, box, 0U, 1U) ||
      !detail::finite_face_neighbour_slabs(
          state.pressure_perturbation.trial, box, 0U, 1U) ||
      detail::output_aliases_flux(system, true, context.mass_flux)) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }
  const PrimitiveHistory histories[]{
      state.density, state.velocity, state.pressure_perturbation,
      state.enthalpy, state.temperature};
  for (const PrimitiveHistory& history : histories) {
    if (output_aliases(history.trial, system) ||
        output_aliases(history.accepted, system) ||
        output_aliases(history.previous, system)) {
      return {StatusCode::invalid_plan, kEnthalpyAssembly};
    }
  }
  const ConstFieldView material_fields[]{material.thermal_conductivity,
                                         material.effective_viscosity,
                                         material.enthalpy_diffusivity,
                                         velocity_gradient};
  for (ConstFieldView field : material_fields) {
    if (output_aliases(field, system)) {
      return {StatusCode::invalid_plan, kEnthalpyAssembly};
    }
  }
  for (std::size_t index = 0U; index < contributions.size; ++index) {
    if (output_aliases(
            contributions.data[index].explicit_source_density, system) ||
        (contributions.data[index].has_implicit_sink &&
         output_aliases(contributions.data[index].implicit_sink_density,
                        system))) {
      return {StatusCode::invalid_plan, kEnthalpyAssembly};
    }
  }

  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        EnthalpyCellTerms terms;
        const Status cell_status = evaluate_enthalpy_cell_terms(
            *plan.kernels_, state, material, velocity_gradient,
            contributions, context, {x, y, z}, false, terms);
        if (!cell_status) return cell_status;
      }
    }
  }

  const std::array<ConstFieldView, 1U> reads{state.enthalpy.trial};
  const std::array<FieldView, 1U> writes{system.residual};
  const KernelInvocation convection{{reads.data(), reads.size()},
                                    {writes.data(), writes.size()}, box,
                                    0U, 0U, 1U, context.face_flux,
                                    context.counters};
  Status status;
  if (context.scope == EquationAssemblyScope::final_conservative) {
    status = cartesian_convection(*plan.kernels_, plan.convection_,
                                  context.mass_flux, convection);
  } else if (context.scope == EquationAssemblyScope::target_coupled) {
    status = cartesian_target_convection(*plan.kernels_, plan.convection_,
                                         context.mass_flux, convection);
  } else {
    status = cartesian_provisional_convection(
        *plan.kernels_, plan.convection_, context.mass_flux, convection);
  }
  if (!status) {
    return status;
  }

  // Temperature-space conduction is evaluated exactly.  The Task 15 driver
  // may later freeze dT/dh into an implicit operator; Task 14 never replaces
  // this term by a generic h diffusion.
  const std::array<ConstFieldView, 1U> thermal_reads{state.temperature.trial};
  const KernelInvocation conduction{{thermal_reads.data(),
                                     thermal_reads.size()},
                                    {writes.data(), writes.size()}, box,
                                    0U, 0U, 1U, 0U, context.counters};

  // Preserve convection in rhs when linear outputs exist; otherwise no
  // caller-provided second scratch exists, so exact fused enthalpy assembly
  // requires the linear pair.
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        system.rhs.unchecked(cell, 0U) =
            system.residual.unchecked(cell, 0U);
      }
    }
  }
  status = cartesian_diffusion(*plan.kernels_, material.thermal_conductivity,
                               conduction);
  if (!status) {
    return status;
  }

  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        EnthalpyCellTerms terms;
        Status cell_status = evaluate_enthalpy_cell_terms(
            *plan.kernels_, state, material, velocity_gradient,
            contributions, context, cell, true, terms);
        EnthalpyCellSystem cell_system;
        if (cell_status)
          cell_status = combine_enthalpy_cell_system(
              *plan.kernels_, state, context, cell,
              system.rhs.unchecked(cell, 0U),
              system.residual.unchecked(cell, 0U), terms, cell_system);
        if (!cell_status) return cell_status;
        system.residual.unchecked(cell, 0U) = cell_system.residual;
        system.diagonal.unchecked(cell, 0U) = cell_system.diagonal;
        system.rhs.unchecked(cell, 0U) = cell_system.rhs;
      }
    }
  }
  status = detail::fill_equation_face_coefficients(
      *plan.kernels_, material.enthalpy_diffusivity, box, system,
      kEnthalpyNumerical);
  if (!status) {
    return status;
  }
  certificate = {plan.fingerprint_, context.scope, context.time,
                 context.geometry, context.face_flux,
                 state_revision(state, material, velocity_gradient,
                                contributions),
                 context.dt};
  return {};
}

Status assemble_tile(
    AssemblyEpoch& epoch, const EnthalpyEquationPlan& plan,
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept {
  if (!epoch.active_) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }
  EquationAssemblyContext context = epoch.context_;
  context.box = box;
  EquationAssemblyCertificate candidate;
  const Status status = assemble_enthalpy_impl(
      plan, state, material, velocity_gradient, contributions, context,
      epoch.system_, candidate, true);
  if (!status) {
    return epoch.fail(status);
  }
  return epoch.record(box, plan.cells(), candidate);
}

Status assemble_enthalpy(
    const EnthalpyEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept {
  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (!detail::full_equation_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }
  EquationAssemblyContext epoch_context = context;
  epoch_context.box = {};
  AssemblyEpoch epoch;
  Status status = epoch.begin(epoch_context, system);
  if (status) {
    status = assemble_tile(epoch, plan, state, material, velocity_gradient,
                           contributions, box);
  }
  if (!status) {
    return status;
  }
  return epoch.finalize(certificate);
}

Status assemble_target_coupled_enthalpy_residual(
    const EnthalpyEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    const EquationAssemblyContext& context, FieldView residual,
    TargetCoupledEnthalpyResidualWorkspace workspace,
    EquationAssemblyCertificate& certificate) noexcept {
  Span<const CompiledContribution> selected_descriptors{};
  if (!detail::select_contribution_stage(
          {plan.contributions_.data(), plan.contributions_.size()},
          context.contribution_stage, selected_descriptors) ||
      selected_descriptors.size != 0U) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }

  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (plan.kernels_ == nullptr || plan.fingerprint_ == 0U ||
      context.geometry != plan.geometry_revision_ ||
      context.boundary != plan.boundary_revision_ ||
      context.thermo != plan.thermodynamics_fingerprint_ ||
      context.transport != plan.transport_fingerprint_ ||
      context.scope != EquationAssemblyScope::target_coupled ||
      context.provisional_mass_flux ||
      !detail::full_equation_box(box, plan.cells_) ||
      !compatible_history(state.density, plan.cells_, plan.density_, 1U, 0U) ||
      !compatible_history(state.velocity, plan.cells_, plan.velocity_, 3U,
                          1U) ||
      !compatible_history(state.pressure_perturbation, plan.cells_,
                          plan.pressure_, 1U, 1U) ||
      !compatible_history(state.enthalpy, plan.cells_, plan.enthalpy_, 1U,
                          plan.convection_reach_) ||
      !compatible_history(state.temperature, plan.cells_, plan.temperature_,
                          1U, 1U) ||
      !detail::valid_cell_view(material.thermal_conductivity, plan.cells_, 0U,
                               1U, 1U) ||
      !detail::valid_cell_view(material.effective_viscosity, plan.cells_, 0U,
                               1U, 0U) ||
      !detail::valid_cell_view(material.enthalpy_diffusivity, plan.cells_, 0U,
                               1U, 1U) ||
      velocity_gradient.field != plan.velocity_gradient_ ||
      !detail::valid_cell_view(velocity_gradient, plan.cells_, 0U, 9U, 0U) ||
      !valid_flux_context(*plan.kernels_, context, box) ||
      !detail::finite_face_flux(context.mass_flux, box) ||
      !detail::finite_face_neighbour_slabs(
          state.enthalpy.trial, box, 0U, 1U, plan.convection_reach_) ||
      !detail::finite_face_neighbour_slabs(state.temperature.trial, box, 0U,
                                           1U) ||
      !detail::finite_face_neighbour_slabs(
          state.pressure_perturbation.trial, box, 0U, 1U) ||
      !detail::valid_cell_view(residual, plan.cells_, 0U, 1U) ||
      !detail::valid_cell_view(workspace.pressure_work, plan.cells_, 0U, 1U) ||
      !detail::valid_cell_view(workspace.viscous_dissipation, plan.cells_, 0U,
                               1U) ||
      !detail::valid_cell_view(workspace.diffusion, plan.cells_, 0U, 1U)) {
    return {StatusCode::invalid_plan, kEnthalpyAssembly};
  }

  const std::array<FieldView, 4U> outputs{
      residual, workspace.pressure_work, workspace.viscous_dissipation,
      workspace.diffusion};
  for (std::size_t left = 0U; left < outputs.size(); ++left) {
    for (std::size_t right = left + 1U; right < outputs.size(); ++right) {
      if (detail::field_views_overlap(as_const(outputs[left]),
                                      as_const(outputs[right]))) {
        return {StatusCode::invalid_plan, kEnthalpyAssembly};
      }
    }
    if (detail::cell_face_views_overlap(outputs[left], context.mass_flux.x) ||
        detail::cell_face_views_overlap(outputs[left], context.mass_flux.y) ||
        detail::cell_face_views_overlap(outputs[left], context.mass_flux.z)) {
      return {StatusCode::invalid_plan, kEnthalpyAssembly};
    }
  }

  const auto aliases_output = [&](ConstFieldView input) noexcept {
    for (FieldView output : outputs) {
      if (detail::field_views_overlap(input, as_const(output))) return true;
    }
    return false;
  };
  const PrimitiveHistory histories[]{
      state.density, state.velocity, state.pressure_perturbation,
      state.enthalpy, state.temperature};
  for (const PrimitiveHistory& history : histories) {
    if (aliases_output(history.trial) || aliases_output(history.accepted) ||
        aliases_output(history.previous)) {
      return {StatusCode::invalid_plan, kEnthalpyAssembly};
    }
  }
  const ConstFieldView material_fields[]{material.thermal_conductivity,
                                         material.effective_viscosity,
                                         material.enthalpy_diffusivity,
                                         velocity_gradient};
  for (ConstFieldView field : material_fields) {
    if (aliases_output(field)) {
      return {StatusCode::invalid_plan, kEnthalpyAssembly};
    }
  }

  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        EnthalpyCellTerms terms;
        const Status cell_status = evaluate_enthalpy_cell_terms(
            *plan.kernels_, state, material, velocity_gradient, {}, context,
            cell, false, terms);
        if (!cell_status) return cell_status;
        workspace.pressure_work.unchecked(cell, 0U) = terms.pressure_work;
        workspace.viscous_dissipation.unchecked(cell, 0U) =
            terms.viscous_dissipation;
      }
    }
  }

  const std::array<ConstFieldView, 1U> enthalpy_reads{state.enthalpy.trial};
  const std::array<FieldView, 1U> convection_writes{residual};
  const KernelInvocation convection{
      {enthalpy_reads.data(), enthalpy_reads.size()},
      {convection_writes.data(), convection_writes.size()}, box,
      0U, 0U, 1U, context.face_flux, context.counters};
  Status status = cartesian_target_convection(
      *plan.kernels_, plan.convection_, context.mass_flux, convection);
  if (!status) return status;

  const std::array<ConstFieldView, 1U> thermal_reads{state.temperature.trial};
  const std::array<FieldView, 1U> diffusion_writes{workspace.diffusion};
  const KernelInvocation conduction{
      {thermal_reads.data(), thermal_reads.size()},
      {diffusion_writes.data(), diffusion_writes.size()}, box,
      0U, 0U, 1U, 0U, context.counters};
  status = cartesian_diffusion(*plan.kernels_, material.thermal_conductivity,
                               conduction);
  if (!status) return status;

  // residual already owns the convection term.  Complete it in place from
  // the preflighted expensive terms and the exact temperature diffusion.
  // It is attempt-local and may therefore be partial on failure; the typed
  // certificate remains unpublished until every cell succeeds.
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho = state.density.trial.unchecked(cell, 0U);
        EnthalpyCellTerms terms;
        terms.pressure_work = workspace.pressure_work.unchecked(cell, 0U);
        terms.viscous_dissipation =
            workspace.viscous_dissipation.unchecked(cell, 0U);
        form_enthalpy_linear_terms(*plan.kernels_, material, context, cell,
                                   rho, terms);
        EnthalpyCellSystem cell_system;
        const Status cell_status = combine_enthalpy_cell_system(
            *plan.kernels_, state, context, cell,
            residual.unchecked(cell, 0U),
            workspace.diffusion.unchecked(cell, 0U), terms, cell_system);
        if (!cell_status) return cell_status;
        residual.unchecked(cell, 0U) = cell_system.residual;
      }
    }
  }
  certificate = {plan.fingerprint_, context.scope, context.time,
                 context.geometry, context.face_flux,
                 state_revision(state, material, velocity_gradient, {}),
                 context.dt};
  return {};
}

}  // namespace hundun::v04
