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

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kRatePlan = 1501U;
constexpr std::uint32_t kRateNumerical = 1502U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool descriptor_matches(const EquationContributionView& view,
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

double pressure_gradient(const CartesianKernelPlan& kernels,
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

bool aliases_output(ConstFieldView input,
                    ThermophysicalRateOutput output) noexcept {
  if (input.base == nullptr) return false;
  const ConstFieldView candidates[]{
      as_const(output.enthalpy_nonadvective_rhs),
      as_const(output.diffusion_scratch),
      as_const(output.scalar_diffusivity_workspace)};
  for (ConstFieldView candidate : candidates) {
    if (detail::field_views_overlap(input, candidate)) return true;
  }
  for (std::size_t index = 0U;
       index < output.species_nonadvective_rhs.size; ++index) {
    if (detail::field_views_overlap(
            input, as_const(output.species_nonadvective_rhs.data[index]))) {
      return true;
    }
  }
  for (std::size_t index = 0U;
       index < output.passive_scalar_nonadvective_rhs.size; ++index) {
    if (detail::field_views_overlap(
            input,
            as_const(output.passive_scalar_nonadvective_rhs.data[index]))) {
      return true;
    }
  }
  return false;
}

}  // namespace

Status evaluate_thermophysical_rates(
    const EquationPlanSet& plans, const ThermophysicalRateInput& input,
    ThermophysicalRateOutput output,
    ThermophysicalRateCertificate& certificate) noexcept {
  const EnthalpyEquationPlan& enthalpy_plan = plans.enthalpy();
  const SpeciesEquationPlan& species_plan = plans.species();
  const ScalarEquationPlan& scalar_plan = plans.scalars();
  const CartesianKernelPlan& kernels = plans.kernels();
  const Int3 cells = enthalpy_plan.cells_;
  const std::size_t scalar_count = species_plan.specs_.size() +
                                   scalar_plan.specs_.size();
  if (plans.fingerprint() == 0U || kernels.fingerprint() == 0U ||
      !same_cells(cells, species_plan.cells_) ||
      !same_cells(cells, scalar_plan.cells_) ||
      input.contribution_stage == 0U ||
      !detail::valid_bdf_coefficients(input.bdf) ||
      output.species_nonadvective_rhs.size != species_plan.specs_.size() ||
      output.passive_scalar_nonadvective_rhs.size !=
          scalar_plan.specs_.size() ||
      (scalar_count != 0U &&
       ((species_plan.specs_.size() != 0U &&
         output.species_nonadvective_rhs.data == nullptr) ||
        (scalar_plan.specs_.size() != 0U &&
         output.passive_scalar_nonadvective_rhs.data == nullptr))) ||
      !detail::valid_cell_view(output.enthalpy_nonadvective_rhs, cells, 0U,
                               1U) ||
      !detail::valid_cell_view(output.diffusion_scratch, cells, 0U, 1U) ||
      !detail::valid_cell_view(as_const(output.scalar_diffusivity_workspace),
                               cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(input.velocity_gradient, cells, 0U, 9U, 0U) ||
      !detail::valid_cell_view(input.material.molecular_viscosity, cells, 0U,
                               1U, 1U) ||
      !detail::valid_cell_view(input.material.effective_viscosity, cells, 0U,
                               1U, 1U) ||
      !detail::valid_cell_view(input.material.thermal_conductivity, cells, 0U,
                               1U, 1U)) {
    return {StatusCode::invalid_plan, kRatePlan};
  }

  const PrimitiveHistory histories[]{
      input.state.density, input.state.velocity,
      input.state.pressure_perturbation, input.state.enthalpy,
      input.state.temperature};
  const std::uint8_t components[]{1U, 3U, 1U, 1U, 1U};
  const FieldId fields[]{enthalpy_plan.density_, enthalpy_plan.velocity_,
                         enthalpy_plan.pressure_, enthalpy_plan.enthalpy_,
                         enthalpy_plan.temperature_};
  for (std::size_t index = 0U; index < 5U; ++index) {
    const PrimitiveHistory history = histories[index];
    if (history.trial.field != fields[index] ||
        history.accepted.field != fields[index] ||
        history.previous.field != fields[index] ||
        !detail::valid_cell_view(history.trial, cells, 0U,
                                 components[index], 1U) ||
        !detail::valid_cell_view(history.accepted, cells, 0U,
                                 components[index], 0U) ||
        !detail::valid_cell_view(history.previous, cells, 0U,
                                 components[index], 0U) ||
        aliases_output(history.trial, output) ||
        aliases_output(history.accepted, output) ||
        aliases_output(history.previous, output)) {
      return {StatusCode::invalid_plan, kRatePlan};
    }
  }

  const FieldView outputs[]{output.enthalpy_nonadvective_rhs,
                            output.diffusion_scratch,
                            output.scalar_diffusivity_workspace};
  for (std::size_t left = 0U; left < 3U; ++left) {
    for (std::size_t right = left + 1U; right < 3U; ++right) {
      if (detail::field_views_overlap(as_const(outputs[left]),
                                      as_const(outputs[right]))) {
        return {StatusCode::invalid_plan, kRatePlan};
      }
    }
  }
  for (std::size_t index = 0U; index < species_plan.specs_.size(); ++index) {
    const FieldView rate = output.species_nonadvective_rhs.data[index];
    if (!detail::valid_cell_view(rate, cells, 0U, 1U)) {
      return {StatusCode::invalid_plan, kRatePlan};
    }
    for (FieldView primary : outputs) {
      if (detail::field_views_overlap(as_const(rate), as_const(primary)))
        return {StatusCode::invalid_plan, kRatePlan};
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (detail::field_views_overlap(
              as_const(rate),
              as_const(output.species_nonadvective_rhs.data[prior])))
        return {StatusCode::invalid_plan, kRatePlan};
    }
  }
  for (std::size_t index = 0U; index < scalar_plan.specs_.size(); ++index) {
    const FieldView rate = output.passive_scalar_nonadvective_rhs.data[index];
    if (!detail::valid_cell_view(rate, cells, 0U, 1U)) {
      return {StatusCode::invalid_plan, kRatePlan};
    }
    for (FieldView primary : outputs) {
      if (detail::field_views_overlap(as_const(rate), as_const(primary)))
        return {StatusCode::invalid_plan, kRatePlan};
    }
    for (std::size_t species = 0U;
         species < output.species_nonadvective_rhs.size; ++species) {
      if (detail::field_views_overlap(
              as_const(rate),
              as_const(output.species_nonadvective_rhs.data[species])))
        return {StatusCode::invalid_plan, kRatePlan};
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (detail::field_views_overlap(
              as_const(rate),
              as_const(output.passive_scalar_nonadvective_rhs.data[prior])))
        return {StatusCode::invalid_plan, kRatePlan};
    }
  }

  if (input.state.independent_species.size != species_plan.specs_.size() ||
      input.state.passive_scalars.size != scalar_plan.specs_.size() ||
      (species_plan.specs_.size() != 0U &&
       input.state.independent_species.data == nullptr) ||
      (scalar_plan.specs_.size() != 0U &&
       input.state.passive_scalars.data == nullptr)) {
    return {StatusCode::invalid_plan, kRatePlan};
  }
  const auto valid_scalar_history = [&](PrimitiveHistory history,
                                        FieldId field) {
    return history.trial.field == field && history.accepted.field == field &&
           history.previous.field == field &&
           detail::valid_cell_view(history.trial, cells, 0U, 1U, 1U) &&
           detail::valid_cell_view(history.accepted, cells, 0U, 1U, 0U) &&
           detail::valid_cell_view(history.previous, cells, 0U, 1U, 0U) &&
           !aliases_output(history.trial, output) &&
           !aliases_output(history.accepted, output) &&
           !aliases_output(history.previous, output);
  };
  for (std::size_t index = 0U; index < species_plan.specs_.size(); ++index) {
    if (!valid_scalar_history(input.state.independent_species.data[index],
                              species_plan.specs_[index].field))
      return {StatusCode::invalid_plan, kRatePlan};
  }
  for (std::size_t index = 0U; index < scalar_plan.specs_.size(); ++index) {
    if (!valid_scalar_history(input.state.passive_scalars.data[index],
                              scalar_plan.specs_[index].field))
      return {StatusCode::invalid_plan, kRatePlan};
  }

  const auto descriptor_count = [&](const CompiledContribution& descriptor) {
    std::size_t matches = 0U;
    for (std::size_t index = 0U; index < input.contributions.size; ++index) {
      if (descriptor_matches(input.contributions.data[index], descriptor,
                             input.contribution_stage)) {
        ++matches;
      }
    }
    return matches;
  };
  std::size_t expected_contributions = 0U;
  const auto validate_descriptors = [&](Span<const CompiledContribution> all) {
    Span<const CompiledContribution> selected;
    if (!detail::select_contribution_stage(
            all, input.contribution_stage, selected)) {
      return false;
    }
    expected_contributions += selected.size;
    for (std::size_t index = 0U; index < selected.size; ++index) {
      if (descriptor_count(selected.data[index]) != 1U) return false;
    }
    return true;
  };
  if (!validate_descriptors(
          {enthalpy_plan.contributions_.data(),
           enthalpy_plan.contributions_.size()}) ||
      !validate_descriptors(
          {species_plan.contributions_.data(),
           species_plan.contributions_.size()}) ||
      !validate_descriptors({scalar_plan.contributions_.data(),
                             scalar_plan.contributions_.size()}) ||
      expected_contributions != input.contributions.size) {
    return {StatusCode::invalid_plan, kRatePlan};
  }
  for (std::size_t index = 0U; index < input.contributions.size; ++index) {
    const EquationContributionView view = input.contributions.data[index];
    if (!detail::valid_cell_view(view.explicit_source_density, cells, 0U,
                                 1U, 0U) ||
        (view.has_implicit_sink &&
         !detail::valid_cell_view(view.implicit_sink_density, cells, 0U,
                                  1U, 0U)) ||
        aliases_output(view.explicit_source_density, output) ||
        (view.has_implicit_sink &&
         aliases_output(view.implicit_sink_density, output))) {
      return {StatusCode::invalid_plan, kRatePlan};
    }
  }
  const ConstFieldView material_inputs[]{
      input.material.molecular_viscosity,
      input.material.effective_viscosity,
      input.material.thermal_conductivity,
      input.velocity_gradient};
  for (ConstFieldView material : material_inputs) {
    if (aliases_output(material, output))
      return {StatusCode::invalid_plan, kRatePlan};
  }

  const KernelBox box{{0, 0, 0}, cells};
  const std::array<ConstFieldView, 1U> thermal_read{
      input.state.temperature.trial};
  const std::array<FieldView, 1U> diffusion_write{output.diffusion_scratch};
  const KernelInvocation conduction{{thermal_read.data(), thermal_read.size()},
                                    {diffusion_write.data(),
                                     diffusion_write.size()},
                                    box, 0U, 0U, 1U, 0U, nullptr};
  Status status = cartesian_diffusion(
      kernels, input.material.thermal_conductivity, conduction);
  if (!status) return status;
  if (input.immersed_interface != nullptr) {
    status = input.immersed_interface->correct_positive_bounded_zero_normal_diffusion(
        input.state.temperature.trial,
        input.material.thermal_conductivity, output.diffusion_scratch);
    if (!status) return status;
  }

  std::uint64_t state_hash = kFnvOffset;
  const auto mix_view = [&](ConstFieldView view) {
    state_hash = mix(state_hash, view.field);
    state_hash = mix(state_hash, view.revision);
    state_hash = mix(state_hash, view.storage_identity);
    state_hash = mix(state_hash, view.revision_domain);
  };
  for (PrimitiveHistory history : histories) {
    mix_view(history.trial);
    mix_view(history.accepted);
    mix_view(history.previous);
  }
  for (std::size_t index = 0U;
       index < input.state.independent_species.size; ++index) {
    const PrimitiveHistory history =
        input.state.independent_species.data[index];
    mix_view(history.trial);
    mix_view(history.accepted);
    mix_view(history.previous);
  }
  for (std::size_t index = 0U; index < input.state.passive_scalars.size;
       ++index) {
    const PrimitiveHistory history = input.state.passive_scalars.data[index];
    mix_view(history.trial);
    mix_view(history.accepted);
    mix_view(history.previous);
  }
  mix_view(input.velocity_gradient);
  mix_view(input.material.molecular_viscosity);
  mix_view(input.material.effective_viscosity);
  mix_view(input.material.thermal_conductivity);
  state_hash = mix(state_hash, bits(input.bdf.a0));
  state_hash = mix(state_hash, bits(input.bdf.a1));
  state_hash = mix(state_hash, bits(input.bdf.a2));
  state_hash = mix(state_hash, input.bdf.order);
  state_hash = mix(state_hash, bits(input.state.pressure_reference));
  state_hash = mix(state_hash, bits(input.state.accepted_pressure_reference));
  state_hash = mix(state_hash, bits(input.state.previous_pressure_reference));
  state_hash = mix(state_hash, input.contribution_stage);
  if (input.immersed_interface != nullptr)
    state_hash = mix(state_hash, input.immersed_interface->fingerprint());
  for (std::size_t index = 0U; index < input.contributions.size; ++index) {
    const EquationContributionView contribution =
        input.contributions.data[index];
    mix_view(contribution.explicit_source_density);
    if (contribution.has_implicit_sink)
      mix_view(contribution.implicit_sink_density);
  }

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        VelocityGradient gradient;
        for (std::uint8_t component = 0U; component < 9U; ++component) {
          gradient.value[component] =
              input.velocity_gradient.unchecked(cell, component);
        }
        // This is a persisted *spatial* rate used by the next predictor.
        // The target-time BDF(p) term is owned exclusively by the coupled
        // BDF(rho*h-p) residual assembled in solver_enthalpy.cpp.  Storing it
        // here as well would EX2-extrapolate a BDF history and recreate the
        // pressure-history feedback loop.
        const std::array<double, 3U> pressure_gradient_value{
            pressure_gradient(kernels,
                              input.state.pressure_perturbation.trial, cell,
                              0U),
            pressure_gradient(kernels,
                              input.state.pressure_perturbation.trial, cell,
                              1U),
            pressure_gradient(kernels,
                              input.state.pressure_perturbation.trial, cell,
                              2U)};
        const std::array<double, 3U> velocity{
            input.state.velocity.trial.unchecked(cell, 0U),
            input.state.velocity.trial.unchecked(cell, 1U),
            input.state.velocity.trial.unchecked(cell, 2U)};
        const double pressure_work =
            velocity[0U] * pressure_gradient_value[0U] +
            velocity[1U] * pressure_gradient_value[1U] +
            velocity[2U] * pressure_gradient_value[2U];
        double dissipation = 0.0;
        if (!std::isfinite(pressure_work) ||
            !newtonian_viscous_dissipation(
                gradient,
                input.material.effective_viscosity.unchecked(cell, 0U),
                dissipation)) {
          return {StatusCode::numerical_failure, kRateNumerical};
        }
        double source = 0.0;
        double sink = 0.0;
        for (std::size_t index = 0U; index < input.contributions.size;
             ++index) {
          const EquationContributionView contribution =
              input.contributions.data[index];
          if (contribution.conserved_quantity != enthalpy_plan.enthalpy_)
            continue;
          source += contribution.explicit_source_density.unchecked(cell, 0U);
          if (contribution.has_implicit_sink) {
            sink += contribution.implicit_sink_density.unchecked(cell, 0U);
          }
        }
        const double rate = pressure_work +
                            output.diffusion_scratch.unchecked(cell, 0U) +
                            dissipation + source -
                            sink * input.state.enthalpy.trial.unchecked(cell,
                                                                        0U);
        if (!std::isfinite(rate))
          return {StatusCode::numerical_failure, kRateNumerical};
        output.enthalpy_nonadvective_rhs.unchecked(cell, 0U) = rate;
      }
    }
  }

  if (input.immersed_interface != nullptr) {
    status = input.immersed_interface->correct_pressure_work(
        input.state.pressure_perturbation.trial, input.state.velocity.trial,
        output.enthalpy_nonadvective_rhs);
    if (!status) return status;
  }

  const auto evaluate_scalar = [&](const ScalarEquationSpec& spec,
                                   ConstFieldView scalar,
                                   FieldView rate) noexcept -> Status {
    const auto fill_coefficient = [&](Int3 cell) {
      const double molecular =
          input.material.molecular_viscosity.unchecked(cell, 0U);
      const double effective =
          input.material.effective_viscosity.unchecked(cell, 0U);
      double turbulent = effective - molecular;
      if (turbulent < 0.0 &&
          turbulent > -1.0e-12 * std::max(1.0, molecular))
        turbulent = 0.0;
      const double coefficient = molecular / spec.molecular_schmidt +
                                 turbulent / spec.turbulent_schmidt;
      if (!std::isfinite(molecular) || molecular < 0.0 ||
          !std::isfinite(effective) || effective < 0.0 || turbulent < 0.0 ||
          !std::isfinite(coefficient) || coefficient <= 0.0) {
        return false;
      }
      output.scalar_diffusivity_workspace.unchecked(cell, 0U) = coefficient;
      return true;
    };
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          if (!fill_coefficient({x, y, z}))
            return {StatusCode::numerical_failure, kRateNumerical};
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        if (!fill_coefficient({-1, y, z}) ||
            !fill_coefficient({cells.x, y, z}))
          return {StatusCode::numerical_failure, kRateNumerical};
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (!fill_coefficient({x, -1, z}) ||
            !fill_coefficient({x, cells.y, z}))
          return {StatusCode::numerical_failure, kRateNumerical};
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        if (!fill_coefficient({x, y, -1}) ||
            !fill_coefficient({x, y, cells.z}))
          return {StatusCode::numerical_failure, kRateNumerical};
    const std::array<ConstFieldView, 1U> reads{scalar};
    const KernelInvocation diffusion{{reads.data(), reads.size()},
                                     {diffusion_write.data(),
                                      diffusion_write.size()},
                                     box, 0U, 0U, 1U, 0U, nullptr};
    Status evaluated = cartesian_diffusion(
        kernels, as_const(output.scalar_diffusivity_workspace), diffusion);
    if (!evaluated) return evaluated;
    if (input.immersed_interface != nullptr) {
      evaluated = input.immersed_interface->correct_zero_normal_diffusion(
          scalar, as_const(output.scalar_diffusivity_workspace),
          output.diffusion_scratch);
      if (!evaluated) return evaluated;
    }
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          double source = 0.0;
          double sink = 0.0;
          for (std::size_t index = 0U; index < input.contributions.size;
               ++index) {
            const EquationContributionView contribution =
                input.contributions.data[index];
            if (contribution.conserved_quantity != spec.field) continue;
            source +=
                contribution.explicit_source_density.unchecked(cell, 0U);
            if (contribution.has_implicit_sink) {
              sink +=
                  contribution.implicit_sink_density.unchecked(cell, 0U);
            }
          }
          const double value = output.diffusion_scratch.unchecked(cell, 0U) +
                               source - sink * scalar.unchecked(cell, 0U);
          if (!std::isfinite(value))
            return {StatusCode::numerical_failure, kRateNumerical};
          rate.unchecked(cell, 0U) = value;
        }
      }
    }
    return {};
  };

  for (std::size_t index = 0U; index < species_plan.specs_.size(); ++index) {
    status = evaluate_scalar(
        species_plan.specs_[index],
        input.state.independent_species.data[index].trial,
        output.species_nonadvective_rhs.data[index]);
    if (!status) return status;
  }
  for (std::size_t index = 0U; index < scalar_plan.specs_.size(); ++index) {
    status = evaluate_scalar(
        scalar_plan.specs_[index],
        input.state.passive_scalars.data[index].trial,
        output.passive_scalar_nonadvective_rhs.data[index]);
    if (!status) return status;
  }

  std::uint64_t rates = mix(state_hash, plans.fingerprint());
  rates = mix(rates, output.enthalpy_nonadvective_rhs.revision);
  rates = mix(rates, output.diffusion_scratch.revision);
  rates = mix(rates, output.scalar_diffusivity_workspace.revision);
  for (std::size_t index = 0U;
       index < output.species_nonadvective_rhs.size; ++index) {
    rates = mix(rates, output.species_nonadvective_rhs.data[index].revision);
  }
  for (std::size_t index = 0U;
       index < output.passive_scalar_nonadvective_rhs.size; ++index) {
    rates =
        mix(rates, output.passive_scalar_nonadvective_rhs.data[index].revision);
  }
  certificate = {plans.fingerprint(),
                 state_hash == 0U ? 1U : state_hash,
                 rates == 0U ? 1U : rates};
  return {};
}

}  // namespace hundun::v04
