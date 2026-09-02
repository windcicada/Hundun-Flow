// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kMomentumAssembly = 1411U;
constexpr std::uint32_t kMomentumNumerical = 1412U;
constexpr std::uint32_t kMomentumSolve = 1413U;
constexpr std::uint32_t kMomentumOperator = 1414U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
bool momentum_partition_alpha_mutation_armed = false;
double momentum_partition_alpha_mutation = 1.0;
#endif

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

Status collective_status(MPI_Comm communicator, Status local) noexcept {
  if (communicator == MPI_COMM_NULL)
    return {StatusCode::invalid_plan, kMomentumAssembly};
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || rank < 0 ||
      rank >= size || size <= 0) {
    return {StatusCode::mpi_failure, kMomentumAssembly};
  }
  const int candidate = local ? size : rank;
  int lowest = size;
  if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kMomentumAssembly};
  }
  if (lowest == size) return {};
  std::uint64_t encoded = 0U;
  if (rank == lowest) {
    encoded = (static_cast<std::uint64_t>(local.code) << 32U) |
              static_cast<std::uint64_t>(local.detail);
  }
  if (MPI_Bcast(&encoded, 1, MPI_UINT64_T, lowest, communicator) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kMomentumAssembly};
  }
  return {static_cast<StatusCode>(encoded >> 32U),
          static_cast<std::uint32_t>(encoded)};
}

RevisionToken state_revision(
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView gradient,
    Span<const EquationContributionView> contributions) noexcept {
  std::uint64_t hash = kFnvOffset;
  const PrimitiveHistory histories[]{state.density, state.velocity,
                                     state.pressure_perturbation};
  for (const PrimitiveHistory& history : histories) {
    const ConstFieldView views[]{history.trial, history.accepted,
                                 history.previous};
    for (ConstFieldView view : views) {
      hash = hash_mix(hash, view.revision);
      hash = hash_mix(hash, view.storage_identity);
      hash = hash_mix(hash, view.revision_domain);
    }
  }
  hash = hash_mix(hash, material.effective_viscosity.revision);
  hash = hash_mix(hash, material.effective_viscosity.storage_identity);
  hash = hash_mix(hash, material.effective_viscosity.revision_domain);
  hash = hash_mix(hash, gradient.revision);
  hash = hash_mix(hash, gradient.storage_identity);
  hash = hash_mix(hash, gradient.revision_domain);
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

bool valid_history(PrimitiveHistory history, Int3 cells, FieldId field,
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

bool valid_context(const CartesianKernelPlan& kernels,
                   const EquationAssemblyContext& context,
                   KernelBox box) noexcept {
  if (!detail::bdf_matches_time_step(context.bdf, context.dt) ||
      context.time == 0U || context.geometry == 0U ||
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
               context.mass_flux.x.storage_identity &&
           context.face_flux_revision_domain ==
               context.mass_flux.x.revision_domain &&
           context.mass_flux.certificate.matches(context.mass_flux);
  }
  return context.scope == EquationAssemblyScope::momentum_predictor &&
         context.provisional_mass_flux &&
         !context.mass_flux.certificate.valid();
}

bool valid_system(EquationSystemView system, Int3 cells,
                  bool& linear) noexcept {
  linear = system.diagonal.base != nullptr || system.rhs.base != nullptr ||
           !detail::equation_faces_empty(system);
  if (!detail::valid_cell_view(system.residual, cells, 0U, 3U)) {
    return false;
  }
  if (!linear) {
    return true;
  }
  return detail::valid_cell_view(system.diagonal, cells, 0U, 3U) &&
         detail::valid_cell_view(system.rhs, cells, 0U, 3U) &&
         detail::valid_equation_faces(system, cells) &&
         detail::equation_face_views_disjoint(system) &&
         !detail::field_views_overlap(as_const(system.diagonal),
                                      as_const(system.rhs)) &&
         !detail::field_views_overlap(as_const(system.diagonal),
                                      as_const(system.residual)) &&
         !detail::field_views_overlap(as_const(system.rhs),
                                      as_const(system.residual));
}

bool output_aliases(ConstFieldView input, EquationSystemView system,
                    bool linear) noexcept {
  return detail::output_aliases_input(input, system, linear);
}

bool valid_contribution(EquationContributionView view,
                        const CompiledContribution& descriptor,
                        StageId stage, Int3 cells,
                        EquationSystemView system, bool linear) noexcept {
  if (view.stage != stage || descriptor.stage != stage ||
      view.conserved_quantity != descriptor.conserved_quantity ||
      !(view.units == descriptor.units) ||
      view.explicit_source_field != descriptor.explicit_source ||
      view.explicit_source_density.field != descriptor.explicit_source ||
      view.has_implicit_sink != descriptor.supplies_implicit_diagonal ||
      !detail::valid_cell_view(view.explicit_source_density, cells, 0U, 3U,
                               0U) ||
      output_aliases(view.explicit_source_density, system, linear)) {
    return false;
  }
  return !view.has_implicit_sink ||
         (view.implicit_sink_field == descriptor.implicit_diagonal &&
          view.implicit_sink_density.field == descriptor.implicit_diagonal &&
          detail::valid_cell_view(view.implicit_sink_density, cells, 0U, 3U,
                                  0U) &&
          !output_aliases(view.implicit_sink_density, system, linear));
}

double derivative_component(const CartesianKernelPlan& kernels,
                            ConstFieldView field, Int3 cell,
                            std::uint8_t component,
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
  return weights.minus * field.unchecked(minus, component) +
         weights.centre * field.unchecked(cell, component) +
         weights.plus * field.unchecked(plus, component);
}

Int3 axis_offset(Int3 value, std::size_t axis, std::int32_t offset) noexcept {
  if (axis == 0U) {
    value.x += offset;
  } else if (axis == 1U) {
    value.y += offset;
  } else {
    value.z += offset;
  }
  return value;
}

double face_cross_traction(const CartesianKernelPlan& kernels,
                           ConstFieldView gradient,
                           ConstFieldView viscosity, std::size_t axis,
                           Int3 face,
                           std::uint8_t momentum_component) noexcept {
  const Int3 left = axis_offset(face, axis, -1);
  const double mu_face =
      kernels.geometry_kind() == GeometryKind::uniform
          ? detail::metric_interpolate_face<true>(
                kernels, axis,
                axis == 0U ? face.x : (axis == 1U ? face.y : face.z),
                viscosity.unchecked(left, 0U),
                viscosity.unchecked(face, 0U))
          : detail::metric_interpolate_face<false>(
                kernels, axis,
                axis == 0U ? face.x : (axis == 1U ? face.y : face.z),
                viscosity.unchecked(left, 0U),
                viscosity.unchecked(face, 0U));
  const auto interpolate_gradient = [&](std::uint8_t component) noexcept {
    return kernels.geometry_kind() == GeometryKind::uniform
               ? detail::metric_interpolate_face<true>(
                     kernels, axis,
                     axis == 0U ? face.x : (axis == 1U ? face.y : face.z),
                     gradient.unchecked(left, component),
                     gradient.unchecked(face, component))
               : detail::metric_interpolate_face<false>(
                     kernels, axis,
                     axis == 0U ? face.x : (axis == 1U ? face.y : face.z),
                     gradient.unchecked(left, component),
                     gradient.unchecked(face, component));
  };
  const double divergence = interpolate_gradient(0U) +
                            interpolate_gradient(4U) +
                            interpolate_gradient(8U);
  // G(component, derivative) is stored at 3*component+derivative.
  const double transpose = interpolate_gradient(
      static_cast<std::uint8_t>(3U * axis + momentum_component));
  return mu_face *
         (transpose - (axis == momentum_component
                           ? (2.0 / 3.0) * divergence
                           : 0.0));
}

double cross_stress_cell_integral(const CartesianKernelPlan& kernels,
                                  ConstFieldView gradient,
                                  ConstFieldView viscosity, Int3 cell,
                                  std::uint8_t component) noexcept {
  double integral = 0.0;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    const Int3 plus = axis_offset(cell, axis, 1);
    const double plus_traction =
        face_cross_traction(kernels, gradient, viscosity, axis, plus,
                            component);
    const double minus_traction =
        face_cross_traction(kernels, gradient, viscosity, axis, cell,
                            component);
    integral += plus_traction *
                    detail::face_area(kernels,
                                      static_cast<CartesianAxis>(axis), plus) -
                minus_traction *
                    detail::face_area(kernels,
                                      static_cast<CartesianAxis>(axis), cell);
  }
  return integral;
}

double outgoing_mass_coefficient(ConstFaceFluxView flux, Int3 cell) noexcept {
  return std::max(flux.x.unchecked({cell.x + 1, cell.y, cell.z}), 0.0) +
         std::max(-flux.x.unchecked(cell), 0.0) +
         std::max(flux.y.unchecked({cell.x, cell.y + 1, cell.z}), 0.0) +
         std::max(-flux.y.unchecked(cell), 0.0) +
         std::max(flux.z.unchecked({cell.x, cell.y, cell.z + 1}), 0.0) +
         std::max(-flux.z.unchecked(cell), 0.0);
}

double first_order_convection_integral(ConstFaceFluxView flux,
                                       ConstFieldView velocity, Int3 cell,
                                       std::uint8_t component) noexcept {
  const Int3 xp{cell.x + 1, cell.y, cell.z};
  const Int3 yp{cell.x, cell.y + 1, cell.z};
  const Int3 zp{cell.x, cell.y, cell.z + 1};
  const Int3 xm{cell.x - 1, cell.y, cell.z};
  const Int3 ym{cell.x, cell.y - 1, cell.z};
  const Int3 zm{cell.x, cell.y, cell.z - 1};
  const double mxm = flux.x.unchecked(cell);
  const double mxp = flux.x.unchecked(xp);
  const double mym = flux.y.unchecked(cell);
  const double myp = flux.y.unchecked(yp);
  const double mzm = flux.z.unchecked(cell);
  const double mzp = flux.z.unchecked(zp);
  const double centre = velocity.unchecked(cell, component);
  return mxp * (mxp >= 0.0 ? centre
                           : velocity.unchecked(xp, component)) -
         mxm * (mxm >= 0.0 ? velocity.unchecked(xm, component) : centre) +
         myp * (myp >= 0.0 ? centre
                           : velocity.unchecked(yp, component)) -
         mym * (mym >= 0.0 ? velocity.unchecked(ym, component) : centre) +
         mzp * (mzp >= 0.0 ? centre
                           : velocity.unchecked(zp, component)) -
         mzm * (mzm >= 0.0 ? velocity.unchecked(zm, component) : centre);
}

std::size_t cell_offset(Int3 cells, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(cell.z));
}

bool fluid_cell(Span<const std::uint8_t> active, Int3 cells,
                Int3 cell) noexcept {
  return active.size == 0U ||
         active.data[cell_offset(cells, cell)] ==
             static_cast<std::uint8_t>(RegionFlag::fluid);
}

bool valid_activity(MgDomainActivityView activity, Int3 cells) noexcept {
  const Int3 x_faces{cells.x + 1, cells.y, cells.z};
  const Int3 y_faces{cells.x, cells.y + 1, cells.z};
  const Int3 z_faces{cells.x, cells.y, cells.z + 1};
  const auto count = [](Int3 shape) noexcept {
    return static_cast<std::size_t>(shape.x) * shape.y * shape.z;
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

bool active_face(MgDomainActivityView activity, Int3 cells,
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
  return faces.data[cell_offset(extents, face)] != 0U;
}

std::size_t face_index(CartesianFace face) noexcept {
  return static_cast<std::size_t>(face);
}

CartesianFace cartesian_face(CartesianAxis axis, bool high) noexcept {
  return axis == CartesianAxis::x
             ? (high ? CartesianFace::x_max : CartesianFace::x_min)
             : (axis == CartesianAxis::y
                    ? (high ? CartesianFace::y_max : CartesianFace::y_min)
                    : (high ? CartesianFace::z_max
                            : CartesianFace::z_min));
}

std::int32_t axis_extent(Int3 cells, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? cells.x
             : (axis == CartesianAxis::y ? cells.y : cells.z);
}

std::int32_t axis_coordinate(Int3 cell, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? cell.x
             : (axis == CartesianAxis::y ? cell.y : cell.z);
}

Int3 shifted(Int3 cell, CartesianAxis axis, int direction) noexcept {
  if (axis == CartesianAxis::x)
    cell.x += direction;
  else if (axis == CartesianAxis::y)
    cell.y += direction;
  else
    cell.z += direction;
  return cell;
}

struct MomentumBoundaryRelations {
  std::array<BoundaryRelation, 6U> relation{};
  std::array<bool, 6U> physical{};
};

Status resolve_momentum_boundary_relations(
    const BoundaryPlan& boundary, FieldId velocity,
    MomentumBoundaryRelations& result) noexcept {
  MomentumBoundaryRelations candidate;
  std::array<bool, 6U> found{};
  for (std::size_t index = 0U; index < candidate.relation.size(); ++index) {
    const auto face = static_cast<CartesianFace>(index);
    const BoundaryFacePlan* plan = nullptr;
    Status status = boundary.face(face, plan);
    if (!status || plan == nullptr) return status ? Status{StatusCode::invalid_plan,
                                                           kMomentumSolve}
                                                 : status;
    candidate.physical[index] = plan->local_owner && !plan->periodic;
  }
  const Span<const BoundaryIndexSpan> spans = boundary.spans();
  for (std::size_t index = 0U; index < spans.size; ++index) {
    const BoundaryIndexSpan span = spans.data[index];
    if (span.stage != BoundaryStage::momentum || span.field != velocity ||
        span.component_begin != 0U || span.component_count != 3U)
      continue;
    const std::size_t selected = face_index(span.face);
    if (selected >= candidate.relation.size() || found[selected])
      return {StatusCode::invalid_plan, kMomentumSolve};
    candidate.relation[selected] = span.relation;
    found[selected] = true;
  }
  for (std::size_t index = 0U; index < candidate.relation.size(); ++index) {
    if (candidate.physical[index] && !found[index])
      return {StatusCode::invalid_plan, kMomentumSolve};
  }
  result = candidate;
  return {};
}

double homogeneous_neighbor(ConstFieldView field, Int3 cells,
                            const MomentumBoundaryRelations& boundary,
                            Int3 cell, CartesianAxis axis, int direction,
                            std::uint8_t velocity_component,
                            std::uint8_t field_component) noexcept {
  const Int3 neighbor = shifted(cell, axis, direction);
  const std::int32_t coordinate = axis_coordinate(neighbor, axis);
  if (coordinate >= 0 && coordinate < axis_extent(cells, axis))
    return field.unchecked(neighbor, field_component);
  const CartesianFace selected = cartesian_face(axis, direction > 0);
  const std::size_t selected_index = face_index(selected);
  if (!boundary.physical[selected_index])
    return field.unchecked(neighbor, field_component);
  const double centre = field.unchecked(cell, field_component);
  switch (boundary.relation[selected_index]) {
    case BoundaryRelation::dirichlet:
    case BoundaryRelation::convective:
      return -centre;
    case BoundaryRelation::zero_gradient:
    case BoundaryRelation::normal_gradient:
      return centre;
    case BoundaryRelation::reflect_normal:
      return velocity_component == static_cast<std::uint8_t>(axis) ? -centre
                                                                   : centre;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

ConstFaceFieldView face_coefficient(EquationSystemView system,
                                    CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? as_const(system.x_coefficient)
             : (axis == CartesianAxis::y
                    ? as_const(system.y_coefficient)
                    : as_const(system.z_coefficient));
}

ConstFaceFieldView face_flux(ConstFaceFluxView flux,
                             CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

double neighbor_coefficient(EquationSystemView system,
                            ConstFaceFluxView flux, Int3 cell,
                            CartesianAxis axis, int direction) noexcept {
  const Int3 face = direction < 0 ? cell : shifted(cell, axis, 1);
  const double diffusion = face_coefficient(system, axis).unchecked(face);
  const double mass = face_flux(flux, axis).unchecked(face);
  const double incoming = direction < 0 ? std::max(mass, 0.0)
                                        : std::max(-mass, 0.0);
  return diffusion + incoming;
}

double momentum_neighbor_sum(
    ConstFieldView field, std::uint8_t field_component,
    std::uint8_t velocity_component, EquationSystemView system,
    ConstFaceFluxView flux, MgDomainActivityView activity,
    const MomentumBoundaryRelations& boundary, Int3 cell) noexcept {
  const Int3 cells = system.diagonal.interior;
  double sum = 0.0;
  const CartesianAxis axes[]{CartesianAxis::x, CartesianAxis::y,
                             CartesianAxis::z};
  for (CartesianAxis axis : axes) {
    for (int direction : {-1, 1}) {
      const Int3 face = direction < 0 ? cell : shifted(cell, axis, 1);
      if (!active_face(activity, cells, axis, face)) continue;
      sum += neighbor_coefficient(system, flux, cell, axis, direction) *
             homogeneous_neighbor(field, cells, boundary, cell, axis,
                                  direction, velocity_component,
                                  field_component);
    }
  }
  return sum;
}

FieldView scalar_component(FieldView field, std::uint8_t component) noexcept {
  field.base += static_cast<std::size_t>(component) * field.component_stride;
  field.components = 1U;
  return field;
}

ConstFieldView scalar_component(ConstFieldView field,
                                std::uint8_t component) noexcept {
  field.base += static_cast<std::size_t>(component) * field.component_stride;
  field.components = 1U;
  return field;
}

class MomentumLinearOperator final : public LinearOperator {
 public:
  MomentumLinearOperator(LinearIdentity identity,
                         PlanFingerprint collective_fingerprint,
                         EquationSystemView system, ConstFaceFluxView flux,
                         MgDomainActivityView activity,
                         MomentumBoundaryRelations boundary,
                         std::uint8_t component, FieldId solution_field,
                         HaloEngine& halo) noexcept
      : system_(system),
        flux_(flux),
        activity_(activity),
        boundary_(boundary),
        component_(component),
        solution_field_(solution_field),
        halo_(&halo) {
    certificate_.identity = identity;
    certificate_.collective_fingerprint = collective_fingerprint;
    certificate_.local_shape = system.diagonal.interior;
    certificate_.operator_class = LinearOperatorClass::nonsymmetric;
  }

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_;
  }

  LinearOperatorFailureProvenance failure_provenance() const noexcept override {
    return failure_;
  }

  Status apply(FieldView x, FieldView y) const noexcept override {
    failure_ = {};
    const Int3 cells = system_.diagonal.interior;
    if (halo_ == nullptr || x.field != solution_field_ || x.components != 1U ||
        y.components != 1U || x.base == nullptr || y.base == nullptr ||
        x.interior.x != cells.x || x.interior.y != cells.y ||
        x.interior.z != cells.z || y.interior.x != cells.x ||
        y.interior.y != cells.y || y.interior.z != cells.z ||
        x.ghosts.x < 1 || x.ghosts.y < 1 || x.ghosts.z < 1 ||
        detail::field_views_overlap(as_const(x), as_const(y))) {
      return {StatusCode::invalid_plan, kMomentumOperator};
    }
    std::array<FieldView, 1U> fields{x};
    HaloTicket ticket;
    Status status = halo_->begin(31U, {fields.data(), fields.size()}, ticket);
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
          const double diagonal =
              system_.diagonal.unchecked(cell, component_);
          double value = diagonal * x.unchecked(cell, 0U);
          if (fluid_cell(activity_.cells, cells, cell)) {
            value -= momentum_neighbor_sum(
                as_const(x), 0U, component_, system_, flux_, activity_,
                boundary_, cell);
          }
          y.unchecked(cell, 0U) = value;
          if (!std::isfinite(diagonal) || !(diagonal > 0.0) ||
              !std::isfinite(value))
            arithmetic = {StatusCode::numerical_failure, kMomentumOperator};
        }
      }
    }
    if (!arithmetic)
      failure_ = {arithmetic, LinearOperatorStatusScope::rank_local, -1};
    return arithmetic;
  }

 private:
  EquationSystemView system_{};
  ConstFaceFluxView flux_{};
  MgDomainActivityView activity_{};
  MomentumBoundaryRelations boundary_{};
  std::uint8_t component_{};
  FieldId solution_field_{};
  HaloEngine* halo_{};
  LinearOperatorCertificate certificate_{};
  mutable LinearOperatorFailureProvenance failure_{};
};

class MomentumJacobiPreconditioner final : public LinearPreconditioner {
 public:
  MomentumJacobiPreconditioner(LinearIdentity identity,
                               PlanFingerprint collective_fingerprint,
                               ConstFieldView diagonal,
                               std::uint8_t component) noexcept
      : diagonal_(diagonal), component_(component) {
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
      return {StatusCode::invalid_plan, kMomentumOperator};
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const double diagonal = diagonal_.unchecked(cell, component_);
          const double value = input.unchecked(cell, 0U) / diagonal;
          if (!std::isfinite(diagonal) || !(diagonal > 0.0) ||
              !std::isfinite(value))
            return {StatusCode::numerical_failure, kMomentumOperator};
          output.unchecked(cell, 0U) = value;
        }
      }
    }
    return {};
  }

 private:
  ConstFieldView diagonal_{};
  std::uint8_t component_{};
  LinearPreconditionerCertificate certificate_{};
};

}  // namespace

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

void arm_momentum_partition_alpha_mutation_for_test(
    bool armed, double alpha) noexcept {
  momentum_partition_alpha_mutation_armed = armed;
  momentum_partition_alpha_mutation = alpha;
}

void clear_momentum_partition_alpha_mutation_for_test() noexcept {
  momentum_partition_alpha_mutation_armed = false;
  momentum_partition_alpha_mutation = 1.0;
}

}  // namespace detail
#endif

Status assemble_momentum_impl(
    const MomentumEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    FieldView low_order_rhs_delta,
    EquationAssemblyCertificate& certificate, bool allow_partial) noexcept {
  Span<const CompiledContribution> selected_descriptors{};
  if (!detail::select_contribution_stage(
          {plan.contributions_.data(), plan.contributions_.size()},
          context.contribution_stage, selected_descriptors)) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  if (plan.kernels_ == nullptr || plan.fingerprint_ == 0U ||
      context.geometry != plan.geometry_revision_ ||
      context.boundary != plan.boundary_revision_ ||
      context.transport != plan.transport_fingerprint_ ||
      !valid_history(state.density, plan.cells_, plan.density_, 1U, 0U) ||
      !valid_history(state.velocity, plan.cells_, plan.velocity_, 3U,
                     plan.convection_reach_) ||
      !valid_history(state.pressure_perturbation, plan.cells_, plan.pressure_,
                     1U, 1U) ||
      material.effective_viscosity.field != plan.viscosity_ ||
      !detail::valid_cell_view(material.effective_viscosity, plan.cells_, 0U,
                               1U, 1U) ||
      velocity_gradient.field != plan.velocity_gradient_ ||
      !detail::valid_cell_view(velocity_gradient, plan.cells_, 0U, 9U, 1U) ||
      contributions.size != selected_descriptors.size ||
      (contributions.size != 0U && contributions.data == nullptr)) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  bool linear = false;
  const KernelBox box = resolved_box(context.box, plan.cells_);
  const bool produce_low_order_delta = low_order_rhs_delta.base != nullptr;
  if (!valid_system(system, plan.cells_, linear) ||
      !valid_context(*plan.kernels_, context, box) ||
      (!allow_partial && !detail::full_equation_box(box, plan.cells_)) ||
      (produce_low_order_delta &&
       (context.scope != EquationAssemblyScope::momentum_predictor ||
        !linear ||
        !detail::valid_cell_view(low_order_rhs_delta, plan.cells_, 0U, 3U) ||
        detail::output_aliases_input(as_const(low_order_rhs_delta), system,
                                     linear) ||
        detail::cell_face_views_overlap(low_order_rhs_delta,
                                        context.mass_flux.x) ||
        detail::cell_face_views_overlap(low_order_rhs_delta,
                                        context.mass_flux.y) ||
        detail::cell_face_views_overlap(low_order_rhs_delta,
                                        context.mass_flux.z))) ||
      !detail::finite_face_flux(context.mass_flux, box) ||
      !detail::finite_face_neighbour_slabs(
          state.velocity.trial, box, 0U, 3U, plan.convection_reach_) ||
      !detail::finite_face_neighbour_slabs(
          state.pressure_perturbation.trial, box, 0U, 1U) ||
      detail::output_aliases_flux(system, linear, context.mass_flux) ||
      output_aliases(material.effective_viscosity, system, linear) ||
      output_aliases(velocity_gradient, system, linear)) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  const PrimitiveHistory histories[]{state.density, state.velocity,
                                     state.pressure_perturbation};
  for (const PrimitiveHistory& history : histories) {
    if (output_aliases(history.trial, system, linear) ||
        output_aliases(history.accepted, system, linear) ||
        output_aliases(history.previous, system, linear) ||
        (produce_low_order_delta &&
         (detail::field_views_overlap(history.trial,
                                      as_const(low_order_rhs_delta)) ||
          detail::field_views_overlap(history.accepted,
                                      as_const(low_order_rhs_delta)) ||
          detail::field_views_overlap(history.previous,
                                      as_const(low_order_rhs_delta))))) {
      return {StatusCode::invalid_plan, kMomentumAssembly};
    }
  }
  if (produce_low_order_delta &&
      (detail::field_views_overlap(material.effective_viscosity,
                                   as_const(low_order_rhs_delta)) ||
       detail::field_views_overlap(velocity_gradient,
                                   as_const(low_order_rhs_delta)))) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  for (std::size_t index = 0U; index < contributions.size; ++index) {
    if (index >= selected_descriptors.size ||
        !valid_contribution(contributions.data[index],
                            selected_descriptors.data[index],
                            context.contribution_stage, plan.cells_, system,
                            linear) ||
        (produce_low_order_delta &&
         (detail::field_views_overlap(
              contributions.data[index].explicit_source_density,
              as_const(low_order_rhs_delta)) ||
          (contributions.data[index].has_implicit_sink &&
           detail::field_views_overlap(
               contributions.data[index].implicit_sink_density,
               as_const(low_order_rhs_delta)))))) {
      return {StatusCode::invalid_plan, kMomentumAssembly};
    }
  }

  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  // Full fail-before-write preflight. Only the six face-neighbour slabs are
  // consumed; corner and edge ghost values are deliberately not required.
  const auto valid_face_cell = [&](Int3 cell) noexcept {
    const double mu = material.effective_viscosity.unchecked(cell, 0U);
    if (!std::isfinite(mu) || mu <= 0.0) {
      return false;
    }
    for (std::uint8_t component = 0U; component < 9U; ++component) {
      if (!std::isfinite(velocity_gradient.unchecked(cell, component))) {
        return false;
      }
    }
    return true;
  };
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        if (!valid_face_cell(cell) ||
            (x == box.begin.x && !valid_face_cell({x - 1, y, z})) ||
            (x + 1 == end.x && !valid_face_cell({x + 1, y, z})) ||
            (y == box.begin.y && !valid_face_cell({x, y - 1, z})) ||
            (y + 1 == end.y && !valid_face_cell({x, y + 1, z})) ||
            (z == box.begin.z && !valid_face_cell({x, y, z - 1})) ||
            (z + 1 == end.z && !valid_face_cell({x, y, z + 1}))) {
          return {StatusCode::numerical_failure, kMomentumNumerical};
        }
      }
    }
  }
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho = state.density.trial.unchecked(cell, 0U);
        const double rho_n = state.density.accepted.unchecked(cell, 0U);
        const double rho_nm1 = state.density.previous.unchecked(cell, 0U);
        const double diffusion = detail::diffusion_diagonal(
            *plan.kernels_, material.effective_viscosity, cell);
        if (!std::isfinite(rho) || !std::isfinite(rho_n) ||
            !std::isfinite(rho_nm1) || rho <= 0.0 || rho_n <= 0.0 ||
            rho_nm1 <= 0.0 || !std::isfinite(diffusion) ||
            diffusion <= 0.0) {
          return {StatusCode::numerical_failure, kMomentumNumerical};
        }
        const double pressure_probe = derivative_component(
            *plan.kernels_, state.pressure_perturbation.trial, cell, 0U,
            0U);
        if (!std::isfinite(pressure_probe) ||
            !std::isfinite(derivative_component(
                *plan.kernels_, state.pressure_perturbation.trial, cell, 0U,
                1U)) ||
            !std::isfinite(derivative_component(
                *plan.kernels_, state.pressure_perturbation.trial, cell, 0U,
                2U))) {
          return {StatusCode::numerical_failure, kMomentumNumerical};
        }
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          if (!std::isfinite(state.velocity.trial.unchecked(cell, component)) ||
              !std::isfinite(
                  state.velocity.accepted.unchecked(cell, component)) ||
              !std::isfinite(
                  state.velocity.previous.unchecked(cell, component))) {
            return {StatusCode::numerical_failure, kMomentumNumerical};
          }
          double source = 0.0;
          double sink = 0.0;
          for (std::size_t index = 0U; index < contributions.size; ++index) {
            source += contributions.data[index]
                          .explicit_source_density.unchecked(cell, component);
            if (contributions.data[index].has_implicit_sink) {
              sink += contributions.data[index]
                          .implicit_sink_density.unchecked(cell, component);
            }
          }
          if (!std::isfinite(source) || !std::isfinite(sink) || sink < 0.0) {
            return {StatusCode::numerical_failure, kMomentumNumerical};
          }
          const double velocity =
              state.velocity.trial.unchecked(cell, component);
          const double unsteady =
              context.bdf.a0 * rho * velocity +
              context.bdf.a1 * rho_n *
                  state.velocity.accepted.unchecked(cell, component) +
              context.bdf.a2 * rho_nm1 *
                  state.velocity.previous.unchecked(cell, component);
          const double pressure = derivative_component(
              *plan.kernels_, state.pressure_perturbation.trial, cell, 0U,
              component);
          const double laplacian = detail::negative_diffusion_operator(
              *plan.kernels_, material.effective_viscosity,
              state.velocity.trial, cell, component);
          const double cross_integral = cross_stress_cell_integral(
              *plan.kernels_, velocity_gradient,
              material.effective_viscosity, cell, component);
          const double volume = detail::cell_volume(*plan.kernels_, cell);
          const double residual =
              (unsteady + pressure - source + sink * velocity) * volume +
              laplacian - cross_integral;
          const double diagonal =
              (context.bdf.a0 * rho + sink) * volume + diffusion;
          if (!std::isfinite(residual) || !std::isfinite(diagonal) ||
              diagonal <= 0.0) {
            return {StatusCode::numerical_failure, kMomentumNumerical};
          }
        }
      }
    }
  }

  const std::array<ConstFieldView, 1U> reads{state.velocity.trial};
  const std::array<FieldView, 1U> writes{system.residual};
  const KernelInvocation convection{{reads.data(), reads.size()},
                                    {writes.data(), writes.size()}, box,
                                    0U, 0U, 3U, context.face_flux,
                                    context.counters};
  Status status =
      context.scope == EquationAssemblyScope::final_conservative
          ? cartesian_convection(*plan.kernels_, plan.convection_,
                                 context.mass_flux, convection)
          : cartesian_provisional_convection(
                *plan.kernels_, plan.convection_, context.mass_flux,
                convection);
  if (!status) {
    return status;
  }

  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho = state.density.trial.unchecked(cell, 0U);
        const double rho_n = state.density.accepted.unchecked(cell, 0U);
        const double rho_nm1 = state.density.previous.unchecked(cell, 0U);
        const double volume = detail::cell_volume(*plan.kernels_, cell);
        const double diffusion_diagonal = detail::diffusion_diagonal(
            *plan.kernels_, material.effective_viscosity, cell);
        const double convective_diagonal =
            linear && context.scope == EquationAssemblyScope::momentum_predictor
                ? outgoing_mass_coefficient(context.mass_flux, cell)
                : 0.0;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double velocity =
              state.velocity.trial.unchecked(cell, component);
          const double unsteady =
              context.bdf.a0 * rho * velocity +
              context.bdf.a1 * rho_n *
                  state.velocity.accepted.unchecked(cell, component) +
              context.bdf.a2 * rho_nm1 *
                  state.velocity.previous.unchecked(cell, component);
          const double pressure = derivative_component(
              *plan.kernels_, state.pressure_perturbation.trial, cell, 0U,
              component);
          const double laplacian = detail::negative_diffusion_operator(
              *plan.kernels_, material.effective_viscosity,
              state.velocity.trial, cell, component);
          const double cross_integral = cross_stress_cell_integral(
              *plan.kernels_, velocity_gradient,
              material.effective_viscosity, cell, component);
          double source = 0.0;
          double sink = 0.0;
          for (std::size_t index = 0U; index < contributions.size; ++index) {
            const EquationContributionView entry = contributions.data[index];
            source +=
                entry.explicit_source_density.unchecked(cell, component);
            if (entry.has_implicit_sink) {
              sink +=
                  entry.implicit_sink_density.unchecked(cell, component);
            }
          }
          const double residual =
              (unsteady + system.residual.unchecked(cell, component) +
               pressure - source + sink * velocity) *
                  volume +
              laplacian - cross_integral;
          const double diagonal =
              (context.bdf.a0 * rho + sink) * volume +
              diffusion_diagonal + convective_diagonal;
          const double rhs = diagonal * velocity - residual;
          const double low_order_delta =
              produce_low_order_delta
                  ? system.residual.unchecked(cell, component) * volume -
                        first_order_convection_integral(
                            context.mass_flux, state.velocity.trial, cell,
                            component)
                  : 0.0;
          if (!std::isfinite(residual) || !std::isfinite(diagonal) ||
              diagonal <= 0.0 || !std::isfinite(rhs) ||
              !std::isfinite(low_order_delta)) {
            return {StatusCode::numerical_failure, kMomentumNumerical};
          }
          system.residual.unchecked(cell, component) = residual;
          if (linear) {
            system.diagonal.unchecked(cell, component) = diagonal;
            system.rhs.unchecked(cell, component) = rhs;
            if (produce_low_order_delta) {
              low_order_rhs_delta.unchecked(cell, component) =
                  low_order_delta;
            }
          }
        }
      }
    }
  }
  if (linear) {
    status = detail::fill_equation_face_coefficients(
        *plan.kernels_, material.effective_viscosity, box, system,
        kMomentumNumerical);
    if (!status) {
      return status;
    }
  }
  RevisionToken assembled_state =
      state_revision(state, material, velocity_gradient, contributions);
  if (context.immersed_interface != nullptr) {
    if (!linear) return {StatusCode::invalid_plan, kMomentumAssembly};
    status = context.immersed_interface->constrain_momentum(
        state.velocity.trial, velocity_gradient,
        state.pressure_perturbation.trial,
        state.density.trial, material.molecular_viscosity,
        material.effective_viscosity, context.wall_treatment,
        {system.diagonal, system.rhs, system.residual});
    if (!status) return status;
    assembled_state = context.immersed_interface->constrain_certificate(
        assembled_state, state.velocity.trial.revision,
        material.effective_viscosity.revision);
    if (assembled_state == 0U)
      return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  certificate = {plan.fingerprint_, context.scope, context.time,
                 context.geometry, context.face_flux, assembled_state,
                 context.dt};
  return {};
}

Status assemble_tile(
    AssemblyEpoch& epoch, const MomentumEquationPlan& plan,
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    KernelBox box, FieldView low_order_rhs_delta) noexcept {
  if (!epoch.active_) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  EquationAssemblyContext context = epoch.context_;
  context.box = box;
  EquationAssemblyCertificate candidate;
  const Status status = assemble_momentum_impl(
      plan, state, material, velocity_gradient, contributions, context,
      epoch.system_, low_order_rhs_delta, candidate, true);
  if (!status) {
    return epoch.fail(status);
  }
  return epoch.record(box, plan.cells(), candidate);
}

Status assemble_tile(
    AssemblyEpoch& epoch, const MomentumEquationPlan& plan,
    const EquationStateView& state, const EquationMaterialView& material,
    ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept {
  return assemble_tile(epoch, plan, state, material, velocity_gradient,
                       contributions, box, {});
}

Status assemble_momentum(
    const MomentumEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept {
  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (!detail::full_equation_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
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

Status assemble_momentum_predictor(
    const MomentumEquationPlan& plan, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView velocity_gradient,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    FieldView low_order_rhs_delta,
    EquationAssemblyCertificate& certificate) noexcept {
  const KernelBox box = resolved_box(context.box, plan.cells());
  if (context.scope != EquationAssemblyScope::momentum_predictor ||
      !detail::full_equation_box(box, plan.cells())) {
    return {StatusCode::invalid_plan, kMomentumAssembly};
  }
  EquationAssemblyContext epoch_context = context;
  epoch_context.box = {};
  AssemblyEpoch epoch;
  Status status = epoch.begin(epoch_context, system);
  if (status) {
    status = assemble_tile(epoch, plan, state, material, velocity_gradient,
                           contributions, box, low_order_rhs_delta);
  }
  if (!status) {
    return status;
  }
  return epoch.finalize(certificate);
}

Status limit_momentum_predictor_correction(
    MPI_Comm communicator, const MomentumEquationPlan& plan,
    const BoundaryPlan& boundary, MeshPatch patch,
    const EquationAssemblyCertificate& assembly, ConstFieldView velocity,
    ConstFieldView density, double step_dt, double convective_cfl_limit,
    ConstFaceFluxView mass_flux, MgDomainActivityView activity,
    EquationSystemView system, MomentumPredictorLimiterWorkspace workspace,
    HaloEngine& limiter_halo,
    ReductionEngine& reductions,
    MomentumPredictorLimiterReport& report) noexcept {
  report = {};
  FieldView cell_workspace = workspace.cell_ratios;
  const std::array<FaceFluxView, 3U> face_workspace =
      workspace.high_order_faces;
  FaceFluxView common_alpha = workspace.common_face_alpha;
  const Int3 cells = system.rhs.interior;
  const std::size_t cell_count =
      cells.x > 0 && cells.y > 0 && cells.z > 0
          ? static_cast<std::size_t>(cells.x) *
                static_cast<std::size_t>(cells.y) *
                static_cast<std::size_t>(cells.z)
          : 0U;
  const std::array<HaloFieldSpec, 1U> halo_fields{{
      {cell_workspace.field, 1U, 3U}}};
  bool linear = false;
  const auto overlaps_flux = [](ConstFieldView cell,
                                ConstFaceFluxView flux) noexcept {
    return detail::cell_face_views_overlap(cell, flux.x) ||
           detail::cell_face_views_overlap(cell, flux.y) ||
           detail::cell_face_views_overlap(cell, flux.z);
  };
  const std::array<ConstFaceFluxView, 3U> limiter_faces{
      as_const(face_workspace[0U]), as_const(face_workspace[1U]),
      as_const(face_workspace[2U])};
  const ConstFaceFluxView alpha_faces = as_const(common_alpha);
  const ConstFaceFluxView equation_faces{
      as_const(system.x_coefficient), as_const(system.y_coefficient),
      as_const(system.z_coefficient), 0U, {}};
  const auto fluxes_overlap = [](ConstFaceFluxView left,
                                 ConstFaceFluxView right) noexcept {
    const ConstFaceFieldView left_faces[]{left.x, left.y, left.z};
    const ConstFaceFieldView right_faces[]{right.x, right.y, right.z};
    for (ConstFaceFieldView left_face : left_faces)
      for (ConstFaceFieldView right_face : right_faces)
        if (detail::face_views_overlap(left_face, right_face)) return true;
    return false;
  };
  const auto overlaps_limiter_faces = [&](ConstFieldView cell) noexcept {
    for (const ConstFaceFluxView faces : limiter_faces)
      if (overlaps_flux(cell, faces)) return true;
    return false;
  };
  const auto limiter_faces_overlap = [&](ConstFaceFluxView faces) noexcept {
    for (const ConstFaceFluxView cached : limiter_faces)
      if (fluxes_overlap(cached, faces)) return true;
    return false;
  };
  bool valid_limiter_faces = true;
  for (std::size_t component = 0U; component < limiter_faces.size();
       ++component) {
    valid_limiter_faces =
        valid_limiter_faces &&
        detail::valid_flux_view(limiter_faces[component], cells,
                                face_workspace[component].revision);
    for (std::size_t prior = 0U; prior < component; ++prior) {
      valid_limiter_faces =
          valid_limiter_faces &&
          !fluxes_overlap(limiter_faces[component], limiter_faces[prior]);
    }
  }
  Status preflight = reductions.validate_communicator(communicator);
  if (preflight && (communicator == MPI_COMM_NULL || cell_count == 0U ||
      plan.kernels_ == nullptr || plan.kernels_->fingerprint() == 0U ||
      plan.cells().x != cells.x || plan.cells().y != cells.y ||
      plan.cells().z != cells.z || boundary.local_cells().x != cells.x ||
      boundary.local_cells().y != cells.y ||
      boundary.local_cells().z != cells.z || patch.cells.x != cells.x ||
      patch.cells.y != cells.y || patch.cells.z != cells.z ||
      !assembly.valid() || assembly.plan != plan.fingerprint() ||
      assembly.scope != EquationAssemblyScope::momentum_predictor ||
      assembly.geometry != plan.geometry_revision_ ||
      boundary.revision() != plan.boundary_revision_ ||
      assembly.dt != step_dt ||
      !valid_system(system, cells, linear) || !linear ||
      velocity.field != plan.velocity_ ||
      !detail::valid_cell_view(velocity, cells, 0U, 3U, 1U) ||
      density.field != plan.density_ ||
      !detail::valid_cell_view(density, cells, 0U, 1U, 0U) ||
      !std::isfinite(step_dt) || !(step_dt > 0.0) ||
      !std::isfinite(convective_cfl_limit) ||
      !(convective_cfl_limit > 0.0) ||
      !detail::valid_cell_view(as_const(cell_workspace), cells, 0U, 3U, 1U) ||
      !detail::valid_flux_view(mass_flux, cells, assembly.face_flux) ||
      !valid_limiter_faces ||
      !detail::valid_flux_view(as_const(common_alpha), cells,
                               common_alpha.revision) ||
      !detail::finite_face_flux(mass_flux, {{0, 0, 0}, cells}) ||
      !valid_activity(activity, cells) ||
      detail::field_views_overlap(velocity, as_const(system.diagonal)) ||
      detail::field_views_overlap(velocity, as_const(system.rhs)) ||
      detail::field_views_overlap(velocity, as_const(system.residual)) ||
      detail::field_views_overlap(velocity, as_const(cell_workspace)) ||
      detail::field_views_overlap(density, as_const(system.diagonal)) ||
      detail::field_views_overlap(density, as_const(system.rhs)) ||
      detail::field_views_overlap(density, as_const(system.residual)) ||
      detail::field_views_overlap(density, as_const(cell_workspace)) ||
      detail::field_views_overlap(as_const(system.diagonal),
                                  as_const(cell_workspace)) ||
      detail::field_views_overlap(as_const(system.rhs),
                                  as_const(cell_workspace)) ||
      detail::field_views_overlap(as_const(system.residual),
                                  as_const(cell_workspace)) ||
      overlaps_limiter_faces(velocity) ||
      overlaps_limiter_faces(as_const(system.diagonal)) ||
      overlaps_limiter_faces(as_const(system.rhs)) ||
      overlaps_limiter_faces(as_const(system.residual)) ||
      overlaps_limiter_faces(as_const(cell_workspace)) ||
      overlaps_limiter_faces(density) ||
      overlaps_flux(velocity, alpha_faces) ||
      overlaps_flux(as_const(system.diagonal), alpha_faces) ||
      overlaps_flux(as_const(system.rhs), alpha_faces) ||
      overlaps_flux(as_const(system.residual), alpha_faces) ||
      overlaps_flux(as_const(cell_workspace), alpha_faces) ||
      overlaps_flux(density, alpha_faces) ||
      overlaps_flux(velocity, mass_flux) ||
      overlaps_flux(as_const(system.diagonal), mass_flux) ||
      overlaps_flux(as_const(system.rhs), mass_flux) ||
      overlaps_flux(as_const(system.residual), mass_flux) ||
      overlaps_flux(as_const(cell_workspace), mass_flux) ||
      overlaps_flux(density, mass_flux) ||
      limiter_faces_overlap(mass_flux) ||
      fluxes_overlap(alpha_faces, mass_flux) ||
      limiter_faces_overlap(alpha_faces) ||
      limiter_faces_overlap(equation_faces) ||
      fluxes_overlap(equation_faces, alpha_faces) ||
      fluxes_overlap(equation_faces, mass_flux) ||
      !limiter_halo.ready() ||
      !limiter_halo.validate_contract(
          communicator, patch, {halo_fields.data(), halo_fields.size()},
          boundary.halo_topology()) ||
      reductions.capacity() < 4U))
    preflight = {StatusCode::invalid_plan, kMomentumAssembly};
  preflight = collective_status(communicator, preflight);
  if (!preflight) return preflight;
  const CartesianKernelPlan& kernels = *plan.kernels_;
  const double certified_dt = assembly.dt;

  MomentumBoundaryRelations relations;
  Status status = resolve_momentum_boundary_relations(
      boundary, velocity.field, relations);
  status = collective_status(communicator, status);
  if (!status) return status;

  const auto fluid = [&](Int3 cell) noexcept {
    return fluid_cell(activity.cells, cells, cell);
  };
  const auto offset = [](Int3 cell, std::size_t axis,
                         std::int32_t amount) noexcept {
    if (axis == 0U) {
      cell.x += amount;
    } else if (axis == 1U) {
      cell.y += amount;
    } else {
      cell.z += amount;
    }
    return cell;
  };
  const auto face_flux = [&](std::size_t axis, Int3 face) noexcept {
    return axis == 0U
               ? mass_flux.x.unchecked(face)
               : (axis == 1U ? mass_flux.y.unchecked(face)
                             : mass_flux.z.unchecked(face));
  };

  const auto physical_face = [&](std::size_t axis, Int3 face) noexcept {
    const std::int32_t coordinate =
        axis == 0U ? face.x : (axis == 1U ? face.y : face.z);
    const std::int32_t extent =
        axis == 0U ? cells.x : (axis == 1U ? cells.y : cells.z);
    if (coordinate != 0 && coordinate != extent) return false;
    const CartesianFace selected = cartesian_face(
        static_cast<CartesianAxis>(axis), coordinate == extent);
    return relations.physical[face_index(selected)];
  };
  const auto frozen_face = [](const FrozenConvectionFaceField& frozen,
                              std::size_t axis) noexcept {
    return axis == 0U ? frozen.x : (axis == 1U ? frozen.y : frozen.z);
  };
  const auto alpha_face = [&](std::size_t axis) noexcept {
    return axis == 0U ? common_alpha.x
                      : (axis == 1U ? common_alpha.y : common_alpha.z);
  };
  const auto inside = [&](Int3 cell) noexcept {
    return cell.x >= 0 && cell.y >= 0 && cell.z >= 0 &&
           cell.x < cells.x && cell.y < cells.y && cell.z < cells.z;
  };
  const auto anti_diffusive_flux = [&](const FrozenConvectionFaceField& frozen,
                                       std::size_t axis, Int3 face,
                                       std::uint8_t component) noexcept {
    const auto selected_axis = static_cast<CartesianAxis>(axis);
    if (!active_face(activity, cells, selected_axis, face))
      return 0.0;
    const double mass = face_flux(axis, face);
    const Int3 donor = mass >= 0.0 ? offset(face, axis, -1) : face;
    // At physical inflow/backflow the frozen face value is the typed
    // boundary authority.  Replacing it by an upwind ghost would alter the
    // prescribed boundary state, so there is no anti-diffusive correction
    // to limit.  Physical outflow has an interior donor and is budgeted
    // against its one adjacent fluid cell below.
    if (physical_face(axis, face) && !inside(donor)) return 0.0;
    const double high = frozen_face(frozen, axis).unchecked(face);
    const double low = velocity.unchecked(donor, component);
    return mass * (high - low);
  };
  const auto neighbour_sum = [&](Int3 cell, Status& local) noexcept {
    double sum = 0.0;
    for (std::size_t axis = 0U; axis < 3U && local; ++axis) {
      const auto selected_axis = static_cast<CartesianAxis>(axis);
      for (int direction : {-1, 1}) {
        const Int3 face = direction < 0 ? cell : offset(cell, axis, 1);
        if (!active_face(activity, cells, selected_axis, face)) continue;
        const double coefficient = neighbor_coefficient(
            system, mass_flux, cell, selected_axis, direction);
        if (!std::isfinite(coefficient) || coefficient < 0.0 ||
            !std::isfinite(sum + coefficient)) {
          local = {StatusCode::numerical_failure, kMomentumNumerical};
          break;
        }
        sum += coefficient;
      }
    }
    return sum;
  };

  double local_cfl_max[2U]{};
  Int3 local_cfl_cell{};
  double local_cfl_rho_volume = 0.0;
  double local_cfl_outgoing = 0.0;
  double local_cfl_absolute = 0.0;
  bool local_cfl_cell_valid = false;
  Status cfl_status;
  for (std::int32_t z = 0; z < cells.z && cfl_status; ++z) {
    for (std::int32_t y = 0; y < cells.y && cfl_status; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (!fluid(cell)) continue;
        const Int3 xp{x + 1, y, z};
        const Int3 yp{x, y + 1, z};
        const Int3 zp{x, y, z + 1};
        const double fxm = mass_flux.x.unchecked(cell);
        const double fxp = mass_flux.x.unchecked(xp);
        const double fym = mass_flux.y.unchecked(cell);
        const double fyp = mass_flux.y.unchecked(yp);
        const double fzm = mass_flux.z.unchecked(cell);
        const double fzp = mass_flux.z.unchecked(zp);
        const double rho = density.unchecked(cell, 0U);
        const double volume = detail::cell_volume(kernels, cell);
        const std::array<double, 6U> cell_flux{
            fxm, fxp, fym, fyp, fzm, fzp};
        const std::array<std::uint8_t, 6U> cell_active{{
            static_cast<std::uint8_t>(
                active_face(activity, cells, CartesianAxis::x, cell)),
            static_cast<std::uint8_t>(
                active_face(activity, cells, CartesianAxis::x, xp)),
            static_cast<std::uint8_t>(
                active_face(activity, cells, CartesianAxis::y, cell)),
            static_cast<std::uint8_t>(
                active_face(activity, cells, CartesianAxis::y, yp)),
            static_cast<std::uint8_t>(
                active_face(activity, cells, CartesianAxis::z, cell)),
            static_cast<std::uint8_t>(
                active_face(activity, cells, CartesianAxis::z, zp)),
        }};
        detail::CellConvectiveCflResult cell_cfl;
        const detail::CellConvectiveCflStatus cell_cfl_status =
            detail::evaluate_cell_convective_cfl(
                rho, volume, cell_flux, cell_active, certified_dt, cell_cfl);
        if (cell_cfl_status != detail::CellConvectiveCflStatus::success) {
          const bool authority_failure =
              cell_cfl_status == detail::CellConvectiveCflStatus::
                                     inactive_nonzero_flux ||
              cell_cfl_status ==
                  detail::CellConvectiveCflStatus::invalid_activity;
          cfl_status = {authority_failure ? StatusCode::invalid_plan
                                          : StatusCode::numerical_failure,
                        authority_failure ? kMomentumAssembly
                                          : kMomentumNumerical};
          break;
        }
        if (!local_cfl_cell_valid || cell_cfl.out > local_cfl_max[0U]) {
          local_cfl_cell_valid = true;
          local_cfl_cell = cell;
          local_cfl_rho_volume = cell_cfl.density_volume;
          local_cfl_outgoing = cell_cfl.outgoing_mass_flow;
          local_cfl_absolute = cell_cfl.absolute_mass_flow;
        }
        local_cfl_max[0U] = std::max(local_cfl_max[0U], cell_cfl.out);
        local_cfl_max[1U] =
            std::max(local_cfl_max[1U], cell_cfl.absolute);
      }
    }
  }
  double global_cfl_max[2U]{};
  status = reductions.checked_max({local_cfl_max, 2U},
                                  {global_cfl_max, 2U}, cfl_status);
  if (!status) return status;
  MomentumAdvectiveCflCertificate advective_cfl;
  advective_cfl.plan = plan.fingerprint();
  advective_cfl.time = assembly.time;
  advective_cfl.density = density.revision;
  advective_cfl.face_flux = mass_flux.revision;
  advective_cfl.density_storage = density.storage_identity;
  advective_cfl.density_revision_domain = density.revision_domain;
  advective_cfl.face_flux_storage = mass_flux.x.storage_identity;
  advective_cfl.face_flux_revision_domain = mass_flux.x.revision_domain;
  advective_cfl.density_view_identity = {density};
  advective_cfl.face_flux_view_identity = {mass_flux};
  advective_cfl.activity_collective = activity.collective_fingerprint;
  advective_cfl.dt = certified_dt;
  advective_cfl.out_max = global_cfl_max[0U];
  advective_cfl.absolute_max = global_cfl_max[1U];
  advective_cfl.limit = convective_cfl_limit;
  constexpr double kCflRoundoffSlack =
      64.0 * std::numeric_limits<double>::epsilon();
  if (global_cfl_max[0U] >
      convective_cfl_limit * (1.0 + kCflRoundoffSlack)) {
    struct MaxLocation {
      double value;
      int rank;
    } local_location{local_cfl_cell_valid ? local_cfl_max[0U] : -1.0, 0},
        global_location{};
    if (MPI_Comm_rank(communicator, &local_location.rank) != MPI_SUCCESS ||
        MPI_Allreduce(&local_location, &global_location, 1, MPI_DOUBLE_INT,
                      MPI_MAXLOC, communicator) != MPI_SUCCESS ||
        global_location.rank < 0) {
      return {StatusCode::mpi_failure, kMomentumAssembly};
    }
    ConvectiveCflFailureWitness witness;
    if (local_location.rank == global_location.rank &&
        local_cfl_cell_valid) {
      witness.valid = true;
      witness.global_cell = {patch.begin.x + local_cfl_cell.x,
                             patch.begin.y + local_cfl_cell.y,
                             patch.begin.z + local_cfl_cell.z};
      witness.rank = local_location.rank;
      witness.out = local_cfl_max[0U];
      witness.absolute = certified_dt * local_cfl_absolute /
                         (2.0 * local_cfl_rho_volume);
      witness.density_volume = local_cfl_rho_volume;
      witness.outgoing_mass_flow = local_cfl_outgoing;
      witness.absolute_mass_flow = local_cfl_absolute;
    }
    if (MPI_Bcast(&witness, static_cast<int>(sizeof(witness)), MPI_BYTE,
                  global_location.rank, communicator) != MPI_SUCCESS ||
        !witness.valid) {
      return {StatusCode::mpi_failure, kMomentumAssembly};
    }
    advective_cfl.failure_witness = witness;
  }
  if (!advective_cfl.valid())
    return {StatusCode::numerical_failure, kMomentumNumerical};
  report.advective_cfl = advective_cfl;

  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    Int3 extents = cells;
    if (axis == 0U)
      ++extents.x;
    else if (axis == 1U)
      ++extents.y;
    else
      ++extents.z;
    FaceFieldView alpha = alpha_face(axis);
    for (std::int32_t z = 0; z < extents.z; ++z)
      for (std::int32_t y = 0; y < extents.y; ++y)
        for (std::int32_t x = 0; x < extents.x; ++x)
          alpha.unchecked({x, y, z}) = 1.0;
  }

  Status local_status;
  constexpr StageId kLimiterHaloStage = 33U;
  constexpr std::uint8_t kNeighbourSumCacheComponent = 2U;
  const RevisionToken base_workspace_revision = cell_workspace.revision;
  std::array<FrozenConvectionFaceField, 3U> frozen_components{};
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    if (local_status) {
      const FrozenConvectionContext context{plan.fingerprint(),
                                            assembly.state};
      const FrozenConvectionFaceOutput output{
          face_workspace[component].x, face_workspace[component].y,
          face_workspace[component].z};
      local_status = freeze_cartesian_target_convection_faces(
          kernels, plan.convection_, mass_flux, velocity, component,
          context, output, frozen_components[component]);
    }
    const FrozenConvectionFaceField& frozen = frozen_components[component];

    if (local_status) {
      for (std::int32_t z = 0; z < cells.z && local_status; ++z) {
        for (std::int32_t y = 0; y < cells.y && local_status; ++y) {
          for (std::int32_t x = 0; x < cells.x && local_status; ++x) {
            const Int3 cell{x, y, z};
            cell_workspace.unchecked(cell, 0U) = 1.0;
            cell_workspace.unchecked(cell, 1U) = 1.0;
            if (component == 0U)
              cell_workspace.unchecked(cell,
                                       kNeighbourSumCacheComponent) = 1.0;
            if (!fluid(cell)) continue;

            const double a = system.diagonal.unchecked(cell, component);
            const double high_rhs = system.rhs.unchecked(cell, component);
            const double centre = velocity.unchecked(cell, component);
            const double row_neighbour_sum =
                component == 0U
                    ? neighbour_sum(cell, local_status)
                    : cell_workspace.unchecked(
                          cell, kNeighbourSumCacheComponent);
            if (!local_status || !std::isfinite(row_neighbour_sum) ||
                row_neighbour_sum < 0.0) {
              if (local_status)
                local_status = {StatusCode::numerical_failure,
                                kMomentumNumerical};
              break;
            }
            // Positive/negative ratios occupy components 0/1.  Component 2
            // is otherwise unused by the limiter, and its interior values
            // survive each ratio halo, so cache the common row sum there.
            if (component == 0U)
              cell_workspace.unchecked(cell,
                                       kNeighbourSumCacheComponent) =
                  row_neighbour_sum;
            const double majorant = std::max(a, row_neighbour_sum);
            const double majorant_rhs = high_rhs + (majorant - a) * centre;
            double delta = 0.0;
            double positive = 0.0;
            double negative = 0.0;
            for (std::size_t axis = 0U; axis < 3U; ++axis) {
              const Int3 plus = offset(cell, axis, 1);
              const double minus_flux = anti_diffusive_flux(
                  frozen, axis, cell, component);
              const double plus_flux = anti_diffusive_flux(
                  frozen, axis, plus, component);
              delta += plus_flux - minus_flux;
              const double contributions[2U]{minus_flux, -plus_flux};
              for (double value : contributions) {
                if (value > 0.0)
                  positive += value;
                else
                  negative += value;
              }
            }
            const double high = majorant_rhs / majorant;
            const double low = (majorant_rhs + delta) / majorant;
            if (!std::isfinite(a) || !(a > 0.0) ||
                !std::isfinite(high_rhs) || !std::isfinite(centre) ||
                !std::isfinite(majorant_rhs) || !std::isfinite(delta) ||
                !std::isfinite(positive) || !std::isfinite(negative) ||
                !std::isfinite(high) || !std::isfinite(low)) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }

            double lower = std::min(low, centre);
            double upper = std::max(low, centre);
            for (std::size_t axis = 0U; axis < 3U; ++axis) {
              const Int3 minus = offset(cell, axis, -1);
              const Int3 plus = offset(cell, axis, 1);
              const auto include = [&](Int3 neighbour, Int3 face) noexcept {
                if (!active_face(activity, cells,
                                 static_cast<CartesianAxis>(axis), face) ||
                    physical_face(axis, face))
                  return true;
                const double value =
                    velocity.unchecked(neighbour, component);
                if (!std::isfinite(value)) return false;
                lower = std::min(lower, value);
                upper = std::max(upper, value);
                return true;
              };
              if (!include(minus, cell) || !include(plus, plus)) {
                local_status = {StatusCode::numerical_failure,
                                kMomentumNumerical};
                break;
              }
            }
            if (!local_status) break;

            // A strict old-time neighbour envelope clips a translating
            // smooth extremum by O(dx^2), which makes the limiter branch
            // change under temporal refinement and destroys BDF2 order.
            // Extend only the extremum-facing side by the exact quadratic
            // three-point curvature allowance.  This is an O(dx^2) spatial
            // truncation allowance (coefficient 1/2 follows from the centred
            // second difference), not a time-step- or Re3900-tuned floor.
            for (std::size_t axis = 0U; axis < 3U; ++axis) {
              const Int3 minus = offset(cell, axis, -1);
              const Int3 plus = offset(cell, axis, 1);
              if (!active_face(activity, cells,
                               static_cast<CartesianAxis>(axis), cell) ||
                  !active_face(activity, cells,
                               static_cast<CartesianAxis>(axis), plus) ||
                  physical_face(axis, cell) || physical_face(axis, plus))
                continue;
              const double minus_value =
                  velocity.unchecked(minus, component);
              const double plus_value = velocity.unchecked(plus, component);
              const double scale = std::max(
                  {1.0, std::abs(centre), std::abs(minus_value),
                   std::abs(plus_value)});
              const double roundoff =
                  64.0 * std::numeric_limits<double>::epsilon() * scale;
              const double curvature =
                  minus_value - 2.0 * centre + plus_value;
              if (!std::isfinite(curvature)) {
                local_status = {StatusCode::numerical_failure,
                                kMomentumNumerical};
                break;
              }
              const double allowance = 0.5 * std::abs(curvature);
              if (curvature < 0.0 &&
                  centre + roundoff >=
                      std::max(minus_value, plus_value))
                upper += allowance;
              else if (curvature > 0.0 &&
                       centre - roundoff <=
                           std::min(minus_value, plus_value))
                lower -= allowance;
            }
            if (!local_status || !std::isfinite(lower) ||
                !std::isfinite(upper)) {
              if (local_status)
                local_status = {StatusCode::numerical_failure,
                                kMomentumNumerical};
              break;
            }

            const double positive_budget =
                std::max(0.0, majorant * (upper - low));
            const double negative_budget =
                std::min(0.0, majorant * (lower - low));
            double positive_ratio =
                positive > 0.0 ? positive_budget / positive : 1.0;
            double negative_ratio =
                negative < 0.0 ? negative_budget / negative : 1.0;
            positive_ratio =
                std::max(0.0, std::min(1.0, positive_ratio));
            negative_ratio =
                std::max(0.0, std::min(1.0, negative_ratio));
            if (!std::isfinite(positive_ratio) ||
                !std::isfinite(negative_ratio)) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }
            cell_workspace.unchecked(cell, 0U) = positive_ratio;
            cell_workspace.unchecked(cell, 1U) = negative_ratio;
          }
        }
      }
    }

    FieldView ratio_view = cell_workspace;
    ratio_view.components = 3U;
    ratio_view.revision = hash_mix(kFnvOffset, assembly.plan);
    ratio_view.revision = hash_mix(ratio_view.revision, assembly.time);
    ratio_view.revision = hash_mix(ratio_view.revision, assembly.geometry);
    ratio_view.revision = hash_mix(ratio_view.revision, assembly.face_flux);
    ratio_view.revision = hash_mix(ratio_view.revision, assembly.state);
    ratio_view.revision =
        hash_mix(ratio_view.revision, double_bits(assembly.dt));
    ratio_view.revision =
        hash_mix(ratio_view.revision, system.diagonal.revision);
    ratio_view.revision =
        hash_mix(ratio_view.revision, base_workspace_revision);
    ratio_view.revision = hash_mix(ratio_view.revision, component + 1U);
    if (ratio_view.revision == 0U) ratio_view.revision = 1U;
    std::array<FieldView, 1U> halo_views{ratio_view};
    HaloTicket ticket;
    Status halo_status = limiter_halo.begin(
        kLimiterHaloStage + component,
        {halo_views.data(), halo_views.size()}, local_status, ticket);
    if (halo_status)
      halo_status = limiter_halo.finish(
          ticket, {halo_views.data(), halo_views.size()});
    if (!halo_status) return halo_status;
    if (halo_views[0U].revision != ratio_view.revision ||
        limiter_halo.ghost_revision(cell_workspace.field) !=
            ratio_view.revision)
      return {StatusCode::invalid_plan, kMomentumAssembly};

    for (std::size_t axis = 0U; axis < 3U && local_status; ++axis) {
      Int3 extents = cells;
      if (axis == 0U)
        ++extents.x;
      else if (axis == 1U)
        ++extents.y;
      else
        ++extents.z;
      for (std::int32_t z = 0; z < extents.z && local_status; ++z) {
        for (std::int32_t y = 0; y < extents.y && local_status; ++y) {
          for (std::int32_t x = 0; x < extents.x && local_status; ++x) {
            const Int3 face{x, y, z};
            const double correction = anti_diffusive_flux(
                frozen, axis, face, component);
            if (!std::isfinite(correction)) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }
            if (correction == 0.0) continue;
            const Int3 left = offset(face, axis, -1);
            const Int3 right = face;
            const double left_contribution = -correction;
            const double right_contribution = correction;
            const auto ratio = [&](Int3 cell, double contribution) noexcept {
              return cell_workspace.unchecked(
                  cell, contribution > 0.0 ? 0U : 1U);
            };
            double alpha = 1.0;
            if (physical_face(axis, face)) {
              if (inside(left) && fluid(left))
                alpha = ratio(left, left_contribution);
              else if (inside(right) && fluid(right))
                alpha = ratio(right, right_contribution);
            } else {
              alpha = std::min(ratio(left, left_contribution),
                               ratio(right, right_contribution));
            }
            if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }
            FaceFieldView common = alpha_face(axis);
            common.unchecked(face) =
                std::min(common.unchecked(face), alpha);
          }
        }
      }
    }
  }

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  // The positive-side copy is deliberately perturbed only in tests.  The
  // following owner publication must replace it with the negative-side
  // authority, proving that shared-face equality is not an accidental
  // consequence of equal cell-ratio arithmetic on the two ranks.
  if (momentum_partition_alpha_mutation_armed && cells.x > 0 &&
      !physical_face(0U, {0, 0, 0})) {
    common_alpha.x.unchecked({0, 0, 0}) =
        momentum_partition_alpha_mutation;
  }
#endif

  // Publish one final alpha authority for every partition and periodic face.
  // The negative-side cell owns its plus face and stores that coefficient in
  // the matching workspace component.  Halo exchange sends the owner value
  // to the positive-side rank, whose local low-face copy consumes it.  This
  // makes the subsequent equal-and-opposite RHS update independent of each
  // rank's redundant local alpha calculation.
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        cell_workspace.unchecked(cell, 0U) =
            common_alpha.x.unchecked({x + 1, y, z});
        cell_workspace.unchecked(cell, 1U) =
            common_alpha.y.unchecked({x, y + 1, z});
        cell_workspace.unchecked(cell, 2U) =
            common_alpha.z.unchecked({x, y, z + 1});
      }
    }
  }
  FieldView published_alpha = cell_workspace;
  published_alpha.components = 3U;
  published_alpha.revision = hash_mix(kFnvOffset, assembly.plan);
  published_alpha.revision =
      hash_mix(published_alpha.revision, assembly.time);
  published_alpha.revision =
      hash_mix(published_alpha.revision, assembly.geometry);
  published_alpha.revision =
      hash_mix(published_alpha.revision, assembly.face_flux);
  published_alpha.revision =
      hash_mix(published_alpha.revision, assembly.state);
  published_alpha.revision =
      hash_mix(published_alpha.revision, double_bits(assembly.dt));
  published_alpha.revision =
      hash_mix(published_alpha.revision, common_alpha.revision);
  published_alpha.revision = hash_mix(published_alpha.revision,
                                      base_workspace_revision);
  if (published_alpha.revision == 0U) published_alpha.revision = 1U;
  std::array<FieldView, 1U> published_views{{published_alpha}};
  HaloTicket published_ticket;
  Status published_status = limiter_halo.begin(
      kLimiterHaloStage + 3U,
      {published_views.data(), published_views.size()}, local_status,
      published_ticket);
  if (published_status)
    published_status = limiter_halo.finish(
        published_ticket,
        {published_views.data(), published_views.size()});
  if (!published_status) return published_status;
  if (published_views[0U].revision != published_alpha.revision ||
      limiter_halo.ghost_revision(cell_workspace.field) !=
          published_alpha.revision)
    return {StatusCode::invalid_plan, kMomentumAssembly};
  if (local_status) {
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      if (physical_face(axis, {0, 0, 0})) continue;
      Int3 extents = cells;
      if (axis == 0U)
        extents.x = 1;
      else if (axis == 1U)
        extents.y = 1;
      else
        extents.z = 1;
      for (std::int32_t z = 0; z < extents.z; ++z) {
        for (std::int32_t y = 0; y < extents.y; ++y) {
          for (std::int32_t x = 0; x < extents.x; ++x) {
            Int3 face{x, y, z};
            if (axis == 0U)
              face.x = 0;
            else if (axis == 1U)
              face.y = 0;
            else
              face.z = 0;
            const Int3 owner = offset(face, axis, -1);
            const double alpha =
                published_alpha.unchecked(owner,
                                          static_cast<std::uint8_t>(axis));
            if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }
            alpha_face(axis).unchecked(face) = alpha;
          }
          if (!local_status) break;
        }
        if (!local_status) break;
      }
      if (!local_status) break;
    }
  }

  // Reuse the haloed cell workspace as one unique active-face marker per
  // axis.  The marker coordinate is the face's negative-side cell; a
  // physical low face therefore uses the one-cell ghost rather than
  // colliding with the first interior plus-face marker.
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component)
          cell_workspace.unchecked({x, y, z}, component) = 0.0;

  // Apply the single synchronized face coefficient to the full momentum
  // correction vector.  Starting from the high-order equation and adding
  // only the rejected fraction preserves its direction and retains the
  // byte-identical high-order fast path when alpha=1 and no majorization is
  // required.
  std::array<double, 4U> local_metrics{};
  double local_minimum_alpha = 1.0;
  // Diffusion and incoming mass flux are shared by all three momentum rows,
  // so compute their row sum once per cell and majorize every component.
  if (local_status) {
    for (std::int32_t z = 0; z < cells.z && local_status; ++z) {
      for (std::int32_t y = 0; y < cells.y && local_status; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (!fluid(cell)) continue;
          const double row_neighbour_sum =
              neighbour_sum(cell, local_status);
          if (!local_status) break;
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const double a =
                system.diagonal.unchecked(cell, component);
            const double majorant = std::max(a, row_neighbour_sum);
            const double adjusted_rhs =
                system.rhs.unchecked(cell, component) +
                (majorant - a) * velocity.unchecked(cell, component);
            if (!std::isfinite(majorant) || !(majorant > 0.0) ||
                !std::isfinite(adjusted_rhs)) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }
            if (majorant != a) {
              system.rhs.unchecked(cell, component) = adjusted_rhs;
              system.diagonal.unchecked(cell, component) = majorant;
            }
          }
          if (!local_status) break;
        }
      }
    }
  }
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    const FrozenConvectionFaceField& frozen = frozen_components[component];
    for (std::size_t axis = 0U; axis < 3U && local_status; ++axis) {
      Int3 extents = cells;
      if (axis == 0U)
        ++extents.x;
      else if (axis == 1U)
        ++extents.y;
      else
        ++extents.z;
      FaceFieldView common = alpha_face(axis);
      for (std::int32_t z = 0; z < extents.z && local_status; ++z) {
        for (std::int32_t y = 0; y < extents.y && local_status; ++y) {
          for (std::int32_t x = 0; x < extents.x; ++x) {
            const Int3 face{x, y, z};
            const double correction = anti_diffusive_flux(
                frozen, axis, face, component);
            const double alpha = common.unchecked(face);
            if (!std::isfinite(correction) || !std::isfinite(alpha) ||
                alpha < 0.0 || alpha > 1.0) {
              local_status = {StatusCode::numerical_failure,
                              kMomentumNumerical};
              break;
            }
            if (correction == 0.0) continue;
            const Int3 left = offset(face, axis, -1);
            const Int3 right = face;
            if (alpha < 1.0) {
              const double rejected = (1.0 - alpha) * correction;
              if (inside(left) && fluid(left))
                system.rhs.unchecked(left, component) += rejected;
              if (inside(right) && fluid(right))
                system.rhs.unchecked(right, component) -= rejected;
            }

            // The negative-side owner uniquely reports this physical face.
            // At an MPI face the positive-side rank sees the owner as a
            // ghost; periodic duplicates are owned by the global high-face
            // representation.
            const Int3 owner =
                physical_face(axis, face) && !inside(left) ? right : left;
            if (inside(owner) && fluid(owner)) {
              const double magnitude = std::abs(correction);
              local_metrics[0U] += magnitude;
              local_metrics[1U] += alpha * magnitude;
              // Count a vector face only when its first active component is
              // encountered; retained L1 remains component-wise above.
              double& active_marker = cell_workspace.unchecked(
                  left, static_cast<std::uint8_t>(axis));
              if (active_marker == 0.0) {
                active_marker = 1.0;
                local_metrics[3U] += 1.0;
                local_minimum_alpha =
                    std::min(local_minimum_alpha, alpha);
                if (alpha < 1.0) local_metrics[2U] += 1.0;
              }
            }
          }
        }
      }
    }
  }

  std::array<double, 4U> global_metrics{};
  status = reductions.checked_sum(
      {local_metrics.data(), local_metrics.size()},
      {global_metrics.data(), global_metrics.size()}, local_status);
  if (!status) return status;
  double global_negative_minimum_alpha = -1.0;
  const double local_negative_minimum_alpha = -local_minimum_alpha;
  status = reductions.checked_max(
      {&local_negative_minimum_alpha, 1U},
      {&global_negative_minimum_alpha, 1U}, {});
  if (!status) return status;
  const double total = global_metrics[0U];
  const double retained = global_metrics[1U];
  const double limited_count = global_metrics[2U];
  const double active_count = global_metrics[3U];
  const double rounded_limited_count = std::nearbyint(limited_count);
  const double rounded_active_count = std::nearbyint(active_count);
  if (!std::isfinite(total) || !std::isfinite(retained) ||
      !std::isfinite(limited_count) || !std::isfinite(active_count) ||
      !std::isfinite(global_negative_minimum_alpha) || total < 0.0 ||
      retained < 0.0 || retained > total || limited_count < 0.0 ||
      active_count < limited_count ||
      limited_count != rounded_limited_count ||
      active_count != rounded_active_count ||
      rounded_limited_count > static_cast<double>(UINT32_MAX) ||
      rounded_active_count > static_cast<double>(UINT32_MAX) ||
      global_negative_minimum_alpha < -1.0 ||
      global_negative_minimum_alpha > 0.0) {
    return {StatusCode::numerical_failure, kMomentumNumerical};
  }
  const std::uint32_t activations =
      static_cast<std::uint32_t>(rounded_limited_count);
  const std::uint32_t active_faces =
      static_cast<std::uint32_t>(rounded_active_count);
  double theta = total > 0.0 ? retained / total : 1.0;
  theta = std::max(0.0, std::min(1.0, theta));
  if (activations == 0U) {
    theta = 1.0;
  } else if (theta == 1.0) {
    theta = std::nextafter(1.0, 0.0);
  }
  report.theta = theta;
  report.activations = activations;
  report.limited = activations != 0U;
  report.correction_metrics_applicable = active_faces != 0U;
  report.minimum_face_alpha =
      active_faces != 0U ? -global_negative_minimum_alpha : 0.0;
  report.active_correction_faces = active_faces;
  report.limited_face_fraction =
      active_faces != 0U
          ? static_cast<double>(activations) /
                static_cast<double>(active_faces)
          : 0.0;
  if (activations != 0U && total > 0.0 && retained == 0.0)
    return {StatusCode::rejected_step, kMomentumNumerical};
  return {};
}

Status solve_momentum_predictor(
    MPI_Comm communicator, const MomentumEquationPlan& plan,
    const BoundaryPlan& boundary,
    MeshPatch patch, const EquationAssemblyCertificate& assembly,
    ConstFaceFluxView mass_flux, MgDomainActivityView activity,
    EquationSystemView system, FieldView velocity, HaloEngine& krylov_halo,
    SolverWorkspace& workspace, ReductionEngine& reductions,
    ResourceCounters* resources,
    MomentumPredictorSolveReport& report) noexcept {
  report = {};
  const Int3 cells = plan.cells();
  const std::size_t count =
      cells.x > 0 && cells.y > 0 && cells.z > 0
          ? static_cast<std::size_t>(cells.x) *
                static_cast<std::size_t>(cells.y) *
                static_cast<std::size_t>(cells.z)
          : 0U;
  const LinearWorkspaceRequirements& requirements = workspace.requirements();
  const FieldView workspace_probe = workspace.vector(0U, cells);
  const std::array<HaloFieldSpec, 1U> halo_fields{{
      {workspace_probe.field, 1U, 1U}}};
  bool linear = false;
  Status preflight = reductions.validate_communicator(communicator);
  if (preflight &&
      (communicator == MPI_COMM_NULL || count == 0U || !assembly.valid() ||
      assembly.plan != plan.fingerprint() ||
      assembly.scope != EquationAssemblyScope::momentum_predictor ||
      assembly.geometry != plan.geometry_revision_ ||
      boundary.revision() != plan.boundary_revision_ ||
      patch.cells.x != cells.x || patch.cells.y != cells.y ||
      patch.cells.z != cells.z || boundary.local_cells().x != cells.x ||
      boundary.local_cells().y != cells.y ||
      boundary.local_cells().z != cells.z ||
      !valid_system(system, cells, linear) || !linear ||
      velocity.field != plan.velocity_ ||
      !detail::valid_cell_view(as_const(velocity), cells, 0U, 3U, 1U) ||
      !detail::valid_flux_view(mass_flux, cells, assembly.face_flux) ||
      !valid_activity(activity, cells) ||
      requirements.algorithm != LinearAlgorithm::fgmres ||
      requirements.maximum_restart == 0U || workspace_probe.base == nullptr ||
      !krylov_halo.ready() ||
      !krylov_halo.validate_contract(communicator, patch,
                                     {halo_fields.data(), halo_fields.size()},
                                     boundary.halo_topology())))
    preflight = {StatusCode::invalid_plan, kMomentumSolve};
  preflight = collective_status(communicator, preflight);
  if (!preflight) return preflight;
  MomentumBoundaryRelations relations;
  Status status = resolve_momentum_boundary_relations(
      boundary, velocity.field, relations);
  status = collective_status(communicator, status);
  if (!status) return status;

  // The limiter publishes H(U_old).  Remove only the homogeneous neighbour
  // action to recover the constant RHS of the actual implicit-upwind system;
  // affine physical-boundary values and the limited high-order correction
  // remain on the right-hand side.  Validate the complete local transform,
  // then reach rank consensus before mutating either U or the assembled RHS.
  Status rhs_status;
  for (std::int32_t z = 0; z < cells.z && rhs_status; ++z) {
    for (std::int32_t y = 0; y < cells.y && rhs_status; ++y) {
      for (std::int32_t x = 0; x < cells.x && rhs_status; ++x) {
        const Int3 cell{x, y, z};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double neighbor =
              fluid_cell(activity.cells, cells, cell)
                  ? momentum_neighbor_sum(as_const(velocity), component,
                                          component, system, mass_flux,
                                          activity, relations, cell)
                  : 0.0;
          const double rhs = system.rhs.unchecked(cell, component) - neighbor;
          if (!std::isfinite(neighbor) || !std::isfinite(rhs)) {
            rhs_status = {StatusCode::numerical_failure, kMomentumSolve};
            break;
          }
        }
      }
    }
  }
  rhs_status = collective_status(communicator, rhs_status);
  if (!rhs_status) return rhs_status;

  // The immersed solid is outside the momentum solve.  Mask the caller's
  // startup/restart value only after the collective validation above, so an
  // inactive identity row cannot seed the Krylov norm and a failing rank
  // cannot publish a partial state.
  if (activity.cells.size != 0U) {
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (fluid_cell(activity.cells, cells, cell)) continue;
          for (std::uint8_t component = 0U; component < 3U; ++component)
            velocity.unchecked(cell, component) = 0.0;
        }
      }
    }
  }

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double neighbor =
              fluid_cell(activity.cells, cells, cell)
                  ? momentum_neighbor_sum(as_const(velocity), component,
                                          component, system, mass_flux,
                                          activity, relations, cell)
                  : 0.0;
          const double rhs = system.rhs.unchecked(cell, component) - neighbor;
          system.rhs.unchecked(cell, component) = rhs;
        }
      }
    }
  }

  const std::uint32_t restart =
      std::min<std::uint32_t>(12U, requirements.maximum_restart);
  const LinearSolveControl control{1.0e-10, 1.0e-4, 64U, 4U, restart};
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    // LinearIdentity is a collective contract.  Boundary ownership and the
    // assembly state certificate intentionally contain rank-local storage
    // facts, so neither may enter this fingerprint.  Their rank-invariant
    // semantic authority is already sealed into the equation plan, while
    // time/geometry/face-flux revisions identify this numeric refill.
    std::uint64_t symbolic = hash_mix(kFnvOffset, plan.fingerprint());
    symbolic = hash_mix(symbolic, component + 1U);
    symbolic = hash_mix(symbolic, activity.collective_fingerprint);
    symbolic = symbolic == 0U ? 1U : symbolic;
    std::uint64_t numeric = hash_mix(symbolic, assembly.time);
    numeric = hash_mix(numeric, assembly.geometry);
    numeric = hash_mix(numeric, assembly.face_flux);
    numeric = hash_mix(numeric, assembly.state);
    numeric = hash_mix(numeric, double_bits(assembly.dt));
    numeric = hash_mix(numeric, system.diagonal.revision);
    numeric = hash_mix(numeric, activity.local_fingerprint);
    numeric = numeric == 0U ? 1U : numeric;
    std::uint64_t hierarchy = hash_mix(symbolic, UINT64_C(1));
    hierarchy = hierarchy == 0U ? 1U : hierarchy;
    LinearIdentity identity{symbolic, numeric, hierarchy,
                            workspace.fingerprint(), 0U};
    std::uint64_t fingerprint = hash_mix(kFnvOffset, identity.symbolic);
    fingerprint = hash_mix(fingerprint, identity.numeric);
    fingerprint = hash_mix(fingerprint, identity.hierarchy);
    fingerprint = hash_mix(fingerprint, identity.workspace);
    identity.fingerprint = fingerprint == 0U ? 1U : fingerprint;
    MomentumLinearOperator linear_operator(
        identity, symbolic, system, mass_flux, activity, relations,
        component, workspace_probe.field, krylov_halo);
    MomentumJacobiPreconditioner preconditioner(
        identity, symbolic, as_const(system.diagonal), component);
    const LinearSolveInvocation invocation{
        scalar_component(as_const(system.rhs), component),
        scalar_component(velocity, component), identity, control, nullptr};
    const LinearSolveResult result = solve_fgmres(
        linear_operator, preconditioner, invocation, workspace, reductions,
        resources);
    report.components[component] = result;
    report.solve_calls = static_cast<std::uint8_t>(component + 1U);
    const bool accepted =
        result.status &&
        (result.termination == LinearTermination::converged ||
         result.termination == LinearTermination::zero_rhs);
    if (!accepted)
      return result.status.code == StatusCode::ok
                 ? Status{StatusCode::rejected_step, kMomentumSolve}
                 : result.status;
  }
  return {};
}

}  // namespace hundun::v04
