// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kScalarPlan = 1451U;
constexpr std::uint32_t kScalarAssembly = 1452U;
constexpr std::uint32_t kScalarNumerical = 1453U;
constexpr std::uint32_t kScalarComposition = 1454U;
constexpr std::uint32_t kScalarBoundary = 1455U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

bool same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool valid_scalar_role(TransportedScalarRole role) noexcept {
  return role == TransportedScalarRole::species ||
         role == TransportedScalarRole::passive_scalar;
}

KernelBox resolved_box(KernelBox requested, Int3 cells) noexcept {
  if (requested.begin.x == 0 && requested.begin.y == 0 &&
      requested.begin.z == 0 && requested.cells.x == 0 &&
      requested.cells.y == 0 && requested.cells.z == 0) {
    return {{0, 0, 0}, cells};
  }
  return requested;
}

bool compatible_history(PrimitiveHistory history, Int3 cells, FieldId field,
                        std::uint8_t trial_ghosts) noexcept {
  return history.trial.field == field && history.accepted.field == field &&
         history.previous.field == field &&
         detail::valid_cell_view(history.trial, cells, 0U, 1U,
                                 trial_ghosts) &&
         detail::valid_cell_view(history.accepted, cells, 0U, 1U, 0U) &&
         detail::valid_cell_view(history.previous, cells, 0U, 1U, 0U);
}

bool valid_scalar_context(const CartesianKernelPlan& kernels,
                          const EquationAssemblyContext& context,
                          KernelBox& box) noexcept {
  box = resolved_box(context.box, kernels.cells());
  if (!detail::bdf_matches_time_step(context.bdf, context.dt) ||
      context.geometry == 0U ||
      context.time == 0U ||
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
  return context.scope == EquationAssemblyScope::momentum_predictor &&
         context.provisional_mass_flux &&
         !context.mass_flux.certificate.valid();
}

bool valid_system(EquationSystemView system, Int3 cells) noexcept {
  const bool all_faces_empty = system.x_coefficient.base == nullptr &&
                               system.y_coefficient.base == nullptr &&
                               system.z_coefficient.base == nullptr;
  const bool all_faces_valid =
      detail::valid_equation_faces(system, cells) &&
      detail::equation_face_views_disjoint(system);
  return detail::valid_cell_view(system.diagonal, cells, 0U, 1U) &&
         detail::valid_cell_view(system.rhs, cells, 0U, 1U) &&
         detail::valid_cell_view(system.residual, cells, 0U, 1U) &&
         (all_faces_empty || all_faces_valid) &&
         !detail::field_views_overlap(as_const(system.diagonal),
                                      as_const(system.rhs)) &&
         !detail::field_views_overlap(as_const(system.diagonal),
                                      as_const(system.residual)) &&
         !detail::field_views_overlap(as_const(system.rhs),
                                      as_const(system.residual));
}

bool valid_material(ConstFieldView diffusivity, Int3 cells) noexcept {
  return detail::valid_cell_view(diffusivity, cells, 0U, 1U, 1U);
}

bool valid_contribution(EquationContributionView contribution, Int3 cells,
                        ConstFieldView scalar, ConstFieldView density,
                        EquationSystemView system) noexcept {
  if (!detail::valid_cell_view(contribution.explicit_source_density, cells,
                               0U, 1U, 0U) ||
      detail::output_aliases_input(contribution.explicit_source_density,
                                   system, true)) {
    return false;
  }
  if (detail::field_views_overlap(scalar, as_const(system.diagonal)) ||
      detail::field_views_overlap(scalar, as_const(system.rhs)) ||
      detail::field_views_overlap(scalar, as_const(system.residual)) ||
      detail::field_views_overlap(density, as_const(system.diagonal)) ||
      detail::field_views_overlap(density, as_const(system.rhs)) ||
      detail::field_views_overlap(density, as_const(system.residual))) {
    return false;
  }
  return !contribution.has_implicit_sink ||
         (detail::valid_cell_view(contribution.implicit_sink_density, cells,
                                  0U, 1U, 0U) &&
          !detail::output_aliases_input(contribution.implicit_sink_density,
                                        system, true));
}

bool matches_contribution(EquationContributionView view,
                          const CompiledContribution& descriptor,
                          StageId stage) noexcept {
  return view.stage == stage && descriptor.stage == stage &&
         view.conserved_quantity == descriptor.conserved_quantity &&
         view.units == descriptor.units &&
         view.explicit_source_field == descriptor.explicit_source &&
         view.explicit_source_density.field == descriptor.explicit_source &&
         view.has_implicit_sink == descriptor.supplies_implicit_diagonal &&
         (!view.has_implicit_sink ||
          (view.implicit_sink_field == descriptor.implicit_diagonal &&
           view.implicit_sink_density.field == descriptor.implicit_diagonal));
}

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

RevisionToken scalar_state_revision(const EquationStateView& state,
                                    PrimitiveHistory scalar,
                                    ConstFieldView diffusivity,
                                    Span<const EquationContributionView>
                                        contributions) noexcept {
  std::uint64_t hash = kFnvOffset;
  const std::array histories{state.density, scalar};
  for (const PrimitiveHistory& history : histories) {
    const ConstFieldView views[]{history.trial, history.accepted,
                                 history.previous};
    for (ConstFieldView view : views) {
      hash = hash_mix(hash, view.revision);
      hash = hash_mix(hash, view.storage_identity);
      hash = hash_mix(hash, view.revision_domain);
    }
  }
  hash = hash_mix(hash, diffusivity.revision);
  hash = hash_mix(hash, diffusivity.storage_identity);
  hash = hash_mix(hash, diffusivity.revision_domain);
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

double positive_transmissibility(const CartesianKernelPlan& kernels,
                                 ConstFieldView diffusivity,
                                 CartesianAxis axis, Int3 face) noexcept {
  const std::int32_t normal = axis == CartesianAxis::x
                                  ? face.x
                                  : (axis == CartesianAxis::y ? face.y
                                                              : face.z);
  Int3 left = face;
  if (axis == CartesianAxis::x) {
    --left.x;
  } else if (axis == CartesianAxis::y) {
    --left.y;
  } else {
    --left.z;
  }
  const double gamma_left = diffusivity.unchecked(left, 0U);
  const double gamma_right = diffusivity.unchecked(face, 0U);
  const double face_coordinate = detail::face_coordinate(kernels, axis, normal);
  const double left_distance =
      face_coordinate - detail::centre_coordinate(kernels, axis, normal - 1);
  const double right_distance =
      detail::centre_coordinate(kernels, axis, normal) - face_coordinate;
  if (!finite_positive(gamma_left) || !finite_positive(gamma_right) ||
      !finite_positive(left_distance) || !finite_positive(right_distance)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return detail::face_area(kernels, axis, face) /
         (left_distance / gamma_left + right_distance / gamma_right);
}

template <CartesianAxis Axis>
Status fill_face_coefficients(const CartesianKernelPlan& kernels,
                              ConstFieldView diffusivity, KernelBox box,
                              FaceFieldView output) noexcept {
  const Int3 cell_end{box.begin.x + box.cells.x,
                      box.begin.y + box.cells.y,
                      box.begin.z + box.cells.z};
  Int3 face_begin = box.begin;
  Int3 face_end = cell_end;
  if constexpr (Axis == CartesianAxis::x) {
    face_begin.x = box.begin.x == 0 ? 0 : box.begin.x + 1;
    ++face_end.x;
  } else if constexpr (Axis == CartesianAxis::y) {
    face_begin.y = box.begin.y == 0 ? 0 : box.begin.y + 1;
    ++face_end.y;
  } else {
    face_begin.z = box.begin.z == 0 ? 0 : box.begin.z + 1;
    ++face_end.z;
  }
  for (std::int32_t z = face_begin.z; z < face_end.z; ++z) {
    for (std::int32_t y = face_begin.y; y < face_end.y; ++y) {
      for (std::int32_t x = face_begin.x; x < face_end.x; ++x) {
        const Int3 face{x, y, z};
        const double coefficient =
            positive_transmissibility(kernels, diffusivity, Axis, face);
        if (!std::isfinite(coefficient)) {
          return {StatusCode::numerical_failure, kScalarNumerical};
        }
        output.unchecked(face) = coefficient;
      }
    }
  }
  return {};
}

Status assemble_transport(
    const CartesianKernelPlan& kernels, Int3 cells, PlanFingerprint fingerprint,
    FieldId density_field, ConvectionScheme convection,
    const ScalarEquationSpec& spec,
    Span<const CompiledContribution> all_descriptors,
    PrimitiveHistory scalar, const EquationStateView& state,
    ConstFieldView diffusivity,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate, bool allow_partial) noexcept {
  Span<const CompiledContribution> descriptors{};
  if (!detail::select_contribution_stage(
          all_descriptors, context.contribution_stage, descriptors)) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  const std::uint8_t required_ghosts =
      convection == ConvectionScheme::central2 ? 1U : 2U;
  KernelBox box{};
  if (fingerprint == 0U || kernels.fingerprint() == 0U ||
      !valid_scalar_role(spec.role) || !finite_positive(spec.molecular_schmidt) ||
      !finite_positive(spec.turbulent_schmidt) ||
      contributions.size != descriptors.size ||
      (contributions.size != 0U && contributions.data == nullptr) ||
      !compatible_history(state.density, cells, density_field, 0U) ||
      !compatible_history(scalar, cells, spec.field, required_ghosts) ||
      !valid_scalar_context(kernels, context, box) ||
      (!allow_partial && !detail::full_equation_box(box, cells)) ||
      !detail::finite_face_flux(context.mass_flux, box) ||
      !detail::finite_face_neighbour_slabs(scalar.trial, box, 0U, 1U,
                                           required_ghosts) ||
      !valid_material(diffusivity, cells) || !valid_system(system, cells) ||
      detail::output_aliases_input(diffusivity, system, true) ||
      detail::output_aliases_flux(system, true, context.mass_flux)) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  const PrimitiveHistory histories[]{state.density, scalar};
  for (const PrimitiveHistory& history : histories) {
    if (detail::output_aliases_input(history.trial, system, true) ||
        detail::output_aliases_input(history.accepted, system, true) ||
        detail::output_aliases_input(history.previous, system, true)) {
      return {StatusCode::invalid_plan, kScalarAssembly};
    }
  }
  for (std::size_t index = 0U; index < contributions.size; ++index) {
    if (!valid_contribution(contributions.data[index], cells, scalar.trial,
                            state.density.trial, system) ||
        !matches_contribution(contributions.data[index],
                              descriptors.data[index],
                              context.contribution_stage)) {
      return {StatusCode::invalid_plan, kScalarAssembly};
    }
  }

  // The per-cell pass below validates every quantity needed by the kernels
  // and assembly before any caller-owned output is modified.
  const Int3 validation_end{
      box.begin.x + box.cells.x, box.begin.y + box.cells.y,
      box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < validation_end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < validation_end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < validation_end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho_trial = state.density.trial.unchecked(cell, 0U);
        const double rho_accepted =
            state.density.accepted.unchecked(cell, 0U);
        const double rho_previous =
            state.density.previous.unchecked(cell, 0U);
        const double q_trial = scalar.trial.unchecked(cell, 0U);
        const double q_accepted = scalar.accepted.unchecked(cell, 0U);
        const double q_previous = scalar.previous.unchecked(cell, 0U);
        double explicit_source = 0.0;
        double implicit_sink = 0.0;
        for (std::size_t index = 0U; index < contributions.size; ++index) {
          const EquationContributionView contribution =
              contributions.data[index];
          explicit_source +=
              contribution.explicit_source_density.unchecked(cell, 0U);
          if (contribution.has_implicit_sink) {
            implicit_sink +=
                contribution.implicit_sink_density.unchecked(cell, 0U);
          }
        }
        const double diffusion_diagonal =
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::x, cell) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::x, {x + 1, y, z}) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::y, cell) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::y, {x, y + 1, z}) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::z, cell) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::z, {x, y, z + 1});
        if (!finite_positive(rho_trial) || !finite_positive(rho_accepted) ||
            !finite_positive(rho_previous) || !std::isfinite(q_trial) ||
            !std::isfinite(q_accepted) || !std::isfinite(q_previous) ||
            !std::isfinite(explicit_source) || !std::isfinite(implicit_sink) ||
            implicit_sink < 0.0 || !finite_positive(diffusion_diagonal)) {
          return {StatusCode::numerical_failure, kScalarNumerical};
        }
        const double volume = detail::cell_volume(kernels, cell);
        const double unsteady =
            context.bdf.a0 * rho_trial * q_trial +
            context.bdf.a1 * rho_accepted * q_accepted +
            context.bdf.a2 * rho_previous * q_previous;
        const double diagonal =
            (context.bdf.a0 * rho_trial + implicit_sink) * volume +
            diffusion_diagonal;
        const double non_diffusive_without_convection =
            (unsteady - explicit_source + implicit_sink * q_trial) * volume;
        if (!finite_positive(diagonal) ||
            !std::isfinite(non_diffusive_without_convection)) {
          return {StatusCode::numerical_failure, kScalarNumerical};
        }
      }
    }
  }

  // Both kernels write the same scratch view.  Convection is consumed before
  // diffusion overwrites it, keeping the hot path allocation-free.
  const std::array<ConstFieldView, 1U> reads{scalar.trial};
  const std::array<FieldView, 1U> writes{system.residual};
  const KernelInvocation invocation{{reads.data(), reads.size()},
                                    {writes.data(), writes.size()}, box,
                                    0U, 0U, 1U, context.face_flux,
                                    context.counters};
  Status evaluated =
      context.scope == EquationAssemblyScope::final_conservative
          ? cartesian_convection(kernels, convection, context.mass_flux,
                                 invocation)
          : cartesian_provisional_convection(
                kernels, convection, context.mass_flux, invocation);
  if (!evaluated) {
    return evaluated;
  }

  // Assemble all non-diffusive pieces while the convection density is live.
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho_trial = state.density.trial.unchecked(cell, 0U);
        const double rho_accepted =
            state.density.accepted.unchecked(cell, 0U);
        const double rho_previous =
            state.density.previous.unchecked(cell, 0U);
        const double q_trial = scalar.trial.unchecked(cell, 0U);
        const double q_accepted = scalar.accepted.unchecked(cell, 0U);
        const double q_previous = scalar.previous.unchecked(cell, 0U);
        double explicit_source = 0.0;
        double implicit_sink = 0.0;
        for (std::size_t index = 0U; index < contributions.size; ++index) {
          const EquationContributionView contribution =
              contributions.data[index];
          explicit_source +=
              contribution.explicit_source_density.unchecked(cell, 0U);
          if (contribution.has_implicit_sink) {
            implicit_sink +=
                contribution.implicit_sink_density.unchecked(cell, 0U);
          }
        }
        const double unsteady =
            context.bdf.a0 * rho_trial * q_trial +
            context.bdf.a1 * rho_accepted * q_accepted +
            context.bdf.a2 * rho_previous * q_previous;
        const double volume = detail::cell_volume(kernels, cell);
        const double diffusion_diagonal =
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::x, cell) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::x, {x + 1, y, z}) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::y, cell) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::y, {x, y + 1, z}) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::z, cell) +
            positive_transmissibility(kernels, diffusivity,
                                      CartesianAxis::z, {x, y, z + 1});
        const double diagonal =
            (context.bdf.a0 * rho_trial + implicit_sink) * volume +
            diffusion_diagonal;
        const double non_diffusive =
            (unsteady + system.residual.unchecked(cell, 0U) -
             explicit_source + implicit_sink * q_trial) *
            volume;
        if (!finite_positive(rho_trial) || !std::isfinite(rho_accepted) ||
            !std::isfinite(rho_previous) || !std::isfinite(q_trial) ||
            !std::isfinite(q_accepted) || !std::isfinite(q_previous) ||
            !std::isfinite(explicit_source) || !std::isfinite(implicit_sink) ||
            implicit_sink < 0.0 || !finite_positive(diffusion_diagonal) ||
            !finite_positive(diagonal) ||
            !std::isfinite(non_diffusive)) {
          return {StatusCode::numerical_failure, kScalarNumerical};
        }
        system.diagonal.unchecked(cell, 0U) = diagonal;
        system.rhs.unchecked(cell, 0U) = non_diffusive;
      }
    }
  }

  KernelInvocation diffusion_call = invocation;
  diffusion_call.required_face_flux_revision = 0U;
  evaluated = cartesian_diffusion(kernels, diffusivity, diffusion_call);
  if (!evaluated) {
    return evaluated;
  }
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double volume = detail::cell_volume(kernels, cell);
        const double residual =
            system.rhs.unchecked(cell, 0U) -
            system.residual.unchecked(cell, 0U) * volume;
        const double rhs =
            system.diagonal.unchecked(cell, 0U) *
                scalar.trial.unchecked(cell, 0U) -
            residual;
        if (!std::isfinite(residual) || !std::isfinite(rhs)) {
          return {StatusCode::numerical_failure, kScalarNumerical};
        }
        system.residual.unchecked(cell, 0U) = residual;
        system.rhs.unchecked(cell, 0U) = rhs;
      }
    }
  }

  if (system.x_coefficient.base != nullptr) {
    evaluated = fill_face_coefficients<CartesianAxis::x>(
        kernels, diffusivity, box, system.x_coefficient);
    if (evaluated) {
      evaluated = fill_face_coefficients<CartesianAxis::y>(
          kernels, diffusivity, box, system.y_coefficient);
    }
    if (evaluated) {
      evaluated = fill_face_coefficients<CartesianAxis::z>(
          kernels, diffusivity, box, system.z_coefficient);
    }
  }
  if (!evaluated) {
    return evaluated;
  }

  certificate = {fingerprint,
                 context.scope,
                 context.time,
                 context.geometry,
                 context.face_flux,
                 scalar_state_revision(state, scalar, diffusivity,
                                       contributions),
                 context.dt};
  return {};
}

}  // namespace

Status close_independent_species(Span<const double> independent,
                                 Span<const double> diffusive_fluxes,
                                 SpeciesClosure& closure) noexcept {
  if ((independent.size != 0U && independent.data == nullptr) ||
      (diffusive_fluxes.size != 0U &&
       (diffusive_fluxes.data == nullptr ||
        diffusive_fluxes.size != independent.size))) {
    return {StatusCode::invalid_plan, kScalarComposition};
  }
  long double independent_sum = 0.0L;
  long double independent_correction = 0.0L;
  long double flux_sum = 0.0L;
  long double flux_correction = 0.0L;
  for (std::size_t index = 0U; index < independent.size; ++index) {
    const double value = independent.data[index];
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
      return {StatusCode::numerical_failure, kScalarComposition};
    }
    const long double corrected =
        static_cast<long double>(value) - independent_correction;
    const long double next = independent_sum + corrected;
    independent_correction = (next - independent_sum) - corrected;
    independent_sum = next;
    if (diffusive_fluxes.size != 0U) {
      const double flux = diffusive_fluxes.data[index];
      if (!std::isfinite(flux)) {
        return {StatusCode::numerical_failure, kScalarComposition};
      }
      const long double corrected_flux =
          static_cast<long double>(flux) - flux_correction;
      const long double next_flux = flux_sum + corrected_flux;
      flux_correction = (next_flux - flux_sum) - corrected_flux;
      flux_sum = next_flux;
    }
  }
  const long double dependent = 1.0L - independent_sum;
  if (dependent < 0.0L || dependent > 1.0L ||
      !std::isfinite(static_cast<double>(dependent)) ||
      !std::isfinite(static_cast<double>(flux_sum))) {
    return {StatusCode::numerical_failure, kScalarComposition};
  }
  const SpeciesClosure candidate{static_cast<double>(dependent),
                                 -static_cast<double>(flux_sum)};
  if (!std::isfinite(candidate.dependent_mass_fraction) ||
      !std::isfinite(candidate.dependent_diffusive_flux)) {
    return {StatusCode::numerical_failure, kScalarComposition};
  }
  closure = candidate;
  return {};
}

Status form_scalar_mass_diffusivity(
    const ScalarEquationSpec& spec, ConstFieldView molecular_viscosity,
    ConstFieldView turbulent_viscosity, KernelBox box,
    FieldView mass_diffusivity) noexcept {
  const Int3 cells = molecular_viscosity.interior;
  if (!valid_scalar_role(spec.role) || !finite_positive(spec.molecular_schmidt) ||
      !finite_positive(spec.turbulent_schmidt) ||
      !same_cells(cells, turbulent_viscosity.interior) ||
      !same_cells(cells, mass_diffusivity.interior) ||
      !detail::valid_cell_view(molecular_viscosity, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(turbulent_viscosity, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(mass_diffusivity, cells, 0U, 1U) ||
      !detail::valid_kernel_box(box, cells) ||
      detail::field_views_overlap(molecular_viscosity,
                                  as_const(mass_diffusivity)) ||
      detail::field_views_overlap(turbulent_viscosity,
                                  as_const(mass_diffusivity))) {
    return {StatusCode::invalid_plan, kScalarPlan};
  }
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  // Validate the entire box before publishing any output.
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        const double molecular = molecular_viscosity.unchecked(cell, 0U);
        const double turbulent = turbulent_viscosity.unchecked(cell, 0U);
        const double value = molecular / spec.molecular_schmidt +
                             turbulent / spec.turbulent_schmidt;
        if (!std::isfinite(molecular) || molecular < 0.0 ||
            !std::isfinite(turbulent) || turbulent < 0.0 ||
            !finite_positive(value)) {
          return {StatusCode::numerical_failure, kScalarNumerical};
        }
      }
    }
  }
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        mass_diffusivity.unchecked(cell, 0U) =
            molecular_viscosity.unchecked(cell, 0U) /
                spec.molecular_schmidt +
            turbulent_viscosity.unchecked(cell, 0U) /
                spec.turbulent_schmidt;
      }
    }
  }
  return {};
}

Status resolve_heat_flux_normal_gradient(double outward_heat_flux,
                                         double conductivity,
                                         double& normal_gradient) noexcept {
  if (!std::isfinite(outward_heat_flux) || !finite_positive(conductivity)) {
    return {StatusCode::numerical_failure, kScalarBoundary};
  }
  const double candidate = -outward_heat_flux / conductivity;
  if (!std::isfinite(candidate)) {
    return {StatusCode::numerical_failure, kScalarBoundary};
  }
  normal_gradient = candidate;
  return {};
}

Status resolve_scalar_flux_normal_gradient(double outward_scalar_flux,
                                           double mass_diffusivity,
                                           double& normal_gradient) noexcept {
  if (!std::isfinite(outward_scalar_flux) ||
      !finite_positive(mass_diffusivity)) {
    return {StatusCode::numerical_failure, kScalarBoundary};
  }
  const double candidate = -outward_scalar_flux / mass_diffusivity;
  if (!std::isfinite(candidate)) {
    return {StatusCode::numerical_failure, kScalarBoundary};
  }
  normal_gradient = candidate;
  return {};
}

Status assemble_species_impl(
    const SpeciesEquationPlan& plan, std::size_t species,
    const EquationStateView& state, const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate, bool allow_partial) noexcept {
  if (plan.kernels_ == nullptr ||
      context.geometry != plan.geometry_revision_ ||
      context.boundary != plan.boundary_revision_ ||
      context.thermo != plan.thermodynamics_fingerprint_ ||
      context.transport != plan.transport_fingerprint_ ||
      context.scope != EquationAssemblyScope::final_conservative ||
      context.provisional_mass_flux ||
      state.independent_species.size != plan.specs_.size() ||
      species >= plan.specs_.size() ||
      species >= plan.contribution_counts_.size() ||
      species >= state.independent_species.size ||
      (state.independent_species.size != 0U &&
       state.independent_species.data == nullptr) ||
      material.scalar_mass_diffusivity.data == nullptr ||
      species >= material.scalar_mass_diffusivity.size ||
      plan.specs_[species].role != TransportedScalarRole::species) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (!detail::valid_kernel_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  for (std::size_t independent = 0U;
       independent < state.independent_species.size; ++independent) {
    if (independent >= plan.specs_.size() ||
        !compatible_history(state.independent_species.data[independent],
                            plan.cells_, plan.specs_[independent].field,
                            plan.convection_ == ConvectionScheme::central2
                                ? 1U
                                : 2U)) {
      return {StatusCode::invalid_plan, kScalarAssembly};
    }
  }
  // Enforce the N-1 composition contract before any scalar output is
  // modified.  This deliberately rejects instead of clipping or normalizing.
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        const Int3 cell{x, y, z};
        long double trial_sum = 0.0L;
        long double accepted_sum = 0.0L;
        long double previous_sum = 0.0L;
        for (std::size_t independent = 0U;
             independent < state.independent_species.size; ++independent) {
          const PrimitiveHistory history =
              state.independent_species.data[independent];
          const double trial = history.trial.unchecked(cell, 0U);
          const double accepted = history.accepted.unchecked(cell, 0U);
          const double previous = history.previous.unchecked(cell, 0U);
          if (!std::isfinite(trial) || !std::isfinite(accepted) ||
              !std::isfinite(previous) || trial < 0.0 || trial > 1.0 ||
              accepted < 0.0 || accepted > 1.0 || previous < 0.0 ||
              previous > 1.0) {
            return {StatusCode::numerical_failure, kScalarComposition};
          }
          trial_sum += trial;
          accepted_sum += accepted;
          previous_sum += previous;
        }
        if (trial_sum > 1.0L || accepted_sum > 1.0L ||
            previous_sum > 1.0L) {
          return {StatusCode::numerical_failure, kScalarComposition};
        }
      }
    }
  }
  return assemble_transport(
      *plan.kernels_, plan.cells_, plan.fingerprint_, plan.density_,
      plan.convection_, plan.specs_[species],
      plan.contribution_counts_[species] == 0U
          ? Span<const CompiledContribution>{}
          : Span<const CompiledContribution>{
                plan.contributions_.data() +
                    plan.contribution_begins_[species],
                plan.contribution_counts_[species]},
      state.independent_species.data[species],
      state, material.scalar_mass_diffusivity.data[species], contributions,
      context, system, certificate, allow_partial);
}

Status assemble_scalar_impl(
    const ScalarEquationPlan& plan, std::size_t scalar,
    const EquationStateView& state, const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate, bool allow_partial) noexcept {
  const std::size_t diffusivity = state.independent_species.size + scalar;
  if (plan.kernels_ == nullptr ||
      context.geometry != plan.geometry_revision_ ||
      context.boundary != plan.boundary_revision_ ||
      context.thermo != plan.thermodynamics_fingerprint_ ||
      context.transport != plan.transport_fingerprint_ ||
      context.scope != EquationAssemblyScope::final_conservative ||
      context.provisional_mass_flux ||
      state.passive_scalars.size != plan.specs_.size() ||
      scalar >= plan.specs_.size() ||
      scalar >= plan.contribution_counts_.size() ||
      scalar >= state.passive_scalars.size ||
      (state.passive_scalars.size != 0U &&
       state.passive_scalars.data == nullptr) ||
      material.scalar_mass_diffusivity.data == nullptr ||
      diffusivity >= material.scalar_mass_diffusivity.size ||
      plan.specs_[scalar].role != TransportedScalarRole::passive_scalar) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  return assemble_transport(
      *plan.kernels_, plan.cells_, plan.fingerprint_, plan.density_,
      plan.convection_, plan.specs_[scalar],
      plan.contribution_counts_[scalar] == 0U
          ? Span<const CompiledContribution>{}
          : Span<const CompiledContribution>{
                plan.contributions_.data() +
                    plan.contribution_begins_[scalar],
                plan.contribution_counts_[scalar]},
      state.passive_scalars.data[scalar], state,
      material.scalar_mass_diffusivity.data[diffusivity], contributions,
      context, system, certificate, allow_partial);
}

Status assemble_tile(
    AssemblyEpoch& epoch, const SpeciesEquationPlan& plan,
    std::size_t species, const EquationStateView& state,
    const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept {
  if (!epoch.active_) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  EquationAssemblyContext context = epoch.context_;
  context.box = box;
  EquationAssemblyCertificate candidate;
  const Status status = assemble_species_impl(
      plan, species, state, material, contributions, context, epoch.system_,
      candidate, true);
  if (!status) {
    return epoch.fail(status);
  }
  return epoch.record(box, plan.cells(), candidate);
}

Status assemble_tile(
    AssemblyEpoch& epoch, const ScalarEquationPlan& plan,
    std::size_t scalar, const EquationStateView& state,
    const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    KernelBox box) noexcept {
  if (!epoch.active_) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  EquationAssemblyContext context = epoch.context_;
  context.box = box;
  EquationAssemblyCertificate candidate;
  const Status status = assemble_scalar_impl(
      plan, scalar, state, material, contributions, context, epoch.system_,
      candidate, true);
  if (!status) {
    return epoch.fail(status);
  }
  return epoch.record(box, plan.cells(), candidate);
}

Status assemble_species(
    const SpeciesEquationPlan& plan, std::size_t species,
    const EquationStateView& state, const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept {
  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (!detail::full_equation_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  EquationAssemblyContext epoch_context = context;
  epoch_context.box = {};
  AssemblyEpoch epoch;
  Status status = epoch.begin(epoch_context, system);
  if (status) {
    status = assemble_tile(epoch, plan, species, state, material,
                           contributions, box);
  }
  if (!status) {
    return status;
  }
  return epoch.finalize(certificate);
}

Status assemble_scalar(
    const ScalarEquationPlan& plan, std::size_t scalar,
    const EquationStateView& state, const EquationMaterialView& material,
    Span<const EquationContributionView> contributions,
    const EquationAssemblyContext& context, EquationSystemView system,
    EquationAssemblyCertificate& certificate) noexcept {
  const KernelBox box = resolved_box(context.box, plan.cells_);
  if (!detail::full_equation_box(box, plan.cells_)) {
    return {StatusCode::invalid_plan, kScalarAssembly};
  }
  EquationAssemblyContext epoch_context = context;
  epoch_context.box = {};
  AssemblyEpoch epoch;
  Status status = epoch.begin(epoch_context, system);
  if (status) {
    status = assemble_tile(epoch, plan, scalar, state, material,
                           contributions, box);
  }
  if (!status) {
    return status;
  }
  return epoch.finalize(certificate);
}

}  // namespace hundun::v04
