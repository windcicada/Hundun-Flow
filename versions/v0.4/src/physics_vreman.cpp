// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_physics.hpp"

#include "field_view_interval_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTurbulencePlan = 2401U;
constexpr std::uint32_t kTurbulenceCollective = 2402U;
constexpr std::uint32_t kTurbulenceView = 2403U;
constexpr std::uint32_t kTurbulenceNumerical = 2404U;
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

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool valid_scalar(ConstFieldView field, Int3 cells) noexcept {
  detail::FieldStorageInterval interval{};
  return field.base != nullptr && same_shape(field.interior, cells) &&
         field.components == 1U && field.revision != 0U &&
         field.storage_identity != 0U && field.revision_domain != 0U &&
         detail::field_storage_interval(field, interval);
}

bool valid_scalar(FieldView field, Int3 cells) noexcept {
  return valid_scalar(as_const(field), cells);
}

bool valid_gradient(ConstFieldView field, Int3 cells) noexcept {
  detail::FieldStorageInterval interval{};
  return field.base != nullptr && same_shape(field.interior, cells) &&
         field.components == 9U && field.revision != 0U &&
         field.storage_identity != 0U && field.revision_domain != 0U &&
         detail::field_storage_interval(field, interval);
}

TurbulenceCandidateFieldBinding candidate_binding(
    ConstFieldView field) noexcept {
  return {field.base,
          field.interior,
          field.ghosts,
          field.components,
          field.stride_y,
          field.stride_z,
          field.component_stride,
          field.replica,
          field.field,
          field.revision,
          field.storage_identity,
          field.revision_domain};
}

std::uint64_t mix_candidate_binding(
    std::uint64_t hash,
    const TurbulenceCandidateFieldBinding& binding) noexcept {
  hash = mix(hash, static_cast<std::uint64_t>(
                       reinterpret_cast<std::uintptr_t>(binding.base)));
  hash = mix(hash, static_cast<std::uint64_t>(binding.interior.x));
  hash = mix(hash, static_cast<std::uint64_t>(binding.interior.y));
  hash = mix(hash, static_cast<std::uint64_t>(binding.interior.z));
  hash = mix(hash, static_cast<std::uint64_t>(binding.ghosts.x));
  hash = mix(hash, static_cast<std::uint64_t>(binding.ghosts.y));
  hash = mix(hash, static_cast<std::uint64_t>(binding.ghosts.z));
  hash = mix(hash, binding.components);
  hash = mix(hash, binding.stride_y);
  hash = mix(hash, binding.stride_z);
  hash = mix(hash, binding.component_stride);
  hash = mix(hash, binding.replica);
  hash = mix(hash, binding.field);
  hash = mix(hash, binding.revision);
  hash = mix(hash, binding.storage);
  hash = mix(hash, binding.revision_domain);
  return hash;
}

Status collective_status(MPI_Comm communicator, Status local) noexcept {
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(local.code) << 32U) | local.detail;
  std::uint64_t global = 0U;
  if (MPI_Allreduce(&packed, &global, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kTurbulenceCollective};
  }
  return {static_cast<StatusCode>(global >> 32U),
          static_cast<std::uint32_t>(global)};
}

bool valid_patch(const CartesianGeometryPlan& geometry,
                 const MeshPatch& patch) noexcept {
  const Int3 global = geometry.global_cells();
  return geometry.fingerprint() != 0U &&
         geometry.topology_revision() != 0U && patch.begin.x >= 0 &&
         patch.begin.y >= 0 && patch.begin.z >= 0 && patch.cells.x > 0 &&
         patch.cells.y > 0 && patch.cells.z > 0 &&
         patch.begin.x <= global.x - patch.cells.x &&
         patch.begin.y <= global.y - patch.cells.y &&
         patch.begin.z <= global.z - patch.cells.z;
}

bool contribution_output_conflict(const ContributionRegistry& registry,
                                  FieldId output) noexcept {
  const Span<const CompiledContribution> contributions =
      registry.contributions();
  return std::any_of(
      contributions.data, contributions.data + contributions.size,
      [output](const CompiledContribution& contribution) {
        return contribution.explicit_source == output ||
               (contribution.supplies_implicit_diagonal &&
                contribution.implicit_diagonal == output);
      });
}

}  // namespace

Status vreman_kinematic_viscosity(const VelocityGradient& gradient,
                                  Real3 filter_widths, double coefficient,
                                  double& out) noexcept {
  if (!std::isfinite(filter_widths.x) || !(filter_widths.x > 0.0) ||
      !std::isfinite(filter_widths.y) || !(filter_widths.y > 0.0) ||
      !std::isfinite(filter_widths.z) || !(filter_widths.z > 0.0) ||
      !std::isfinite(coefficient) || coefficient < 0.0) {
    return {StatusCode::invalid_plan, kTurbulencePlan};
  }
  double alpha[3][3]{};
  double alpha_squared = 0.0;
  for (std::size_t derivative = 0U; derivative < 3U; ++derivative) {
    for (std::size_t velocity = 0U; velocity < 3U; ++velocity) {
      // Stored layout is velocity-major: d U_velocity / d x_derivative.
      const double value = gradient.value[velocity * 3U + derivative];
      if (!std::isfinite(value)) {
        return {StatusCode::numerical_failure, kTurbulenceNumerical};
      }
      alpha[derivative][velocity] = value;
      alpha_squared += value * value;
    }
  }
  if (!std::isfinite(alpha_squared)) {
    return {StatusCode::numerical_failure, kTurbulenceNumerical};
  }
  if (coefficient == 0.0 || alpha_squared == 0.0) {
    out = 0.0;
    return {};
  }
  const std::array<double, 3U> delta_squared{
      filter_widths.x * filter_widths.x,
      filter_widths.y * filter_widths.y,
      filter_widths.z * filter_widths.z};
  double beta[3][3]{};
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      for (std::size_t derivative = 0U; derivative < 3U; ++derivative) {
        beta[i][j] += delta_squared[derivative] *
                      alpha[derivative][i] * alpha[derivative][j];
      }
    }
  }
  double invariant = beta[0][0] * beta[1][1] - beta[0][1] * beta[0][1] +
                     beta[0][0] * beta[2][2] - beta[0][2] * beta[0][2] +
                     beta[1][1] * beta[2][2] - beta[1][2] * beta[1][2];
  const double scale = std::max(
      {1.0, std::abs(beta[0][0] * beta[1][1]),
       std::abs(beta[0][0] * beta[2][2]),
       std::abs(beta[1][1] * beta[2][2])});
  if (invariant < 0.0 &&
      invariant >= -128.0 * std::numeric_limits<double>::epsilon() * scale) {
    invariant = 0.0;
  }
  const double candidate =
      coefficient * std::sqrt(invariant / alpha_squared);
  if (!std::isfinite(invariant) || invariant < 0.0 ||
      !std::isfinite(candidate) || candidate < 0.0) {
    return {StatusCode::numerical_failure, kTurbulenceNumerical};
  }
  out = candidate;
  return {};
}

Status TurbulencePlan::compile(MPI_Comm communicator,
                               const TurbulencePlanSpec& spec,
                               const CartesianGeometryPlan& geometry,
                               const MeshPatch& patch,
                               FieldId effective_viscosity_output,
                               StageId update_stage,
                               ContributionRegistry& contributions,
                               TurbulencePlan& out) noexcept {
  SubgridKind subgrid = SubgridKind::none;
  WallTreatmentKind wall = WallTreatmentKind::resolved;
  double coefficient = 0.0;
  if (spec.kind == TurbulenceKind::wale) {
    subgrid = SubgridKind::wale;
    coefficient = spec.wale_coefficient;
  } else if (spec.kind == TurbulenceKind::vreman_wall_function) {
    subgrid = SubgridKind::vreman;
    wall = WallTreatmentKind::equilibrium_wall_function;
    coefficient = spec.vreman_coefficient;
  }
  const bool supported = spec.kind == TurbulenceKind::none ||
                         spec.kind == TurbulenceKind::wale ||
                         spec.kind == TurbulenceKind::vreman_wall_function;
  const bool declared = std::binary_search(contributions.declared_fields_.begin(),
                                           contributions.declared_fields_.end(),
                                           effective_viscosity_output);
  Status local =
      communicator == MPI_COMM_NULL || !supported ||
              !valid_patch(geometry, patch) || update_stage == 0U ||
              !std::isfinite(spec.wale_coefficient) ||
              spec.wale_coefficient < 0.0 ||
              !std::isfinite(spec.vreman_coefficient) ||
              spec.vreman_coefficient < 0.0 ||
              !std::isfinite(spec.turbulent_prandtl) ||
              !(spec.turbulent_prandtl > 0.0) ||
              !std::isfinite(spec.turbulent_schmidt) ||
              !(spec.turbulent_schmidt > 0.0) || contributions.frozen_ ||
              contributions.declared_fields_.empty() || !declared ||
              contributions.effective_viscosity_claimed_ ||
              contribution_output_conflict(contributions,
                                           effective_viscosity_output)
          ? Status{StatusCode::invalid_plan, kTurbulencePlan}
          : Status{};
  local = collective_status(communicator, local);
  if (!local) {
    return local;
  }

  std::uint64_t semantic = kFnvOffset;
  semantic = mix(semantic, geometry.fingerprint());
  semantic = mix(semantic, static_cast<std::uint64_t>(spec.kind));
  semantic = mix(semantic, bits(spec.wale_coefficient));
  semantic = mix(semantic, bits(spec.vreman_coefficient));
  semantic = mix(semantic, bits(spec.turbulent_prandtl));
  semantic = mix(semantic, bits(spec.turbulent_schmidt));
  semantic = mix(semantic, effective_viscosity_output);
  semantic = mix(semantic, update_stage);
  semantic = semantic == 0U ? 1U : semantic;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  if (MPI_Allreduce(&semantic, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(&semantic, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kTurbulenceCollective};
  }
  if (minimum != maximum) {
    return {StatusCode::invalid_plan, kTurbulenceCollective};
  }

  TurbulencePlan candidate;
  const std::size_t cell_count = static_cast<std::size_t>(patch.cells.x) *
                                 patch.cells.y * patch.cells.z;
  try {
    candidate.filter_metrics_.resize(cell_count);
    candidate.pending_effective_viscosity_.resize(cell_count);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kTurbulencePlan};
  } catch (...) {
    local = {StatusCode::invalid_plan, kTurbulencePlan};
  }
  local = collective_status(communicator, local);
  if (!local) {
    return local;
  }
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x, ++flat) {
        const std::size_t gx = static_cast<std::size_t>(patch.begin.x + x);
        const std::size_t gy = static_cast<std::size_t>(patch.begin.y + y);
        const std::size_t gz = static_cast<std::size_t>(patch.begin.z + z);
        const double dx = geometry.x().widths().data[gx];
        const double dy = geometry.y().widths().data[gy];
        const double dz = geometry.z().widths().data[gz];
        candidate.filter_metrics_[flat] =
            {dx, dy, dz, std::cbrt(dx * dy * dz)};
      }
    }
  }
  local = candidate.authority_.claim(effective_viscosity_output, update_stage,
                                     contributions);
  local = collective_status(communicator, local);
  if (!local) {
    return local;
  }
  candidate.cells_ = patch.cells;
  candidate.kind_ = spec.kind;
  candidate.subgrid_ = subgrid;
  candidate.wall_ = wall;
  candidate.coefficient_ = coefficient;
  candidate.turbulent_prandtl_ = spec.turbulent_prandtl;
  candidate.turbulent_schmidt_ = spec.turbulent_schmidt;
  candidate.fingerprint_ = semantic;
  out = std::move(candidate);
  return {};
}

Status TurbulencePlan::update(const TurbulenceUpdateInput& input,
                              FieldView effective_viscosity,
                              TurbulenceCertificate& certificate) noexcept {
  const std::size_t cell_count = filter_metrics_.size();
  const bool gradient_span = input.velocity_gradient.data != nullptr ||
                             input.velocity_gradient.size != 0U;
  const bool gradient_field = input.velocity_gradient_field.base != nullptr;
  if (fingerprint_ == 0U || !authority_.claimed() ||
      !valid_scalar(input.density, cells_) ||
      !valid_scalar(input.molecular_viscosity, cells_) ||
      !valid_scalar(effective_viscosity, cells_) ||
      gradient_span == gradient_field ||
      (gradient_span && (input.velocity_gradient.data == nullptr ||
                         input.velocity_gradient.size != cell_count)) ||
      (gradient_field &&
       (!valid_gradient(input.velocity_gradient_field, cells_) ||
        input.velocity_gradient_field.revision != input.gradient_revision)) ||
      input.gradient_revision == 0U ||
      detail::field_views_overlap(input.density,
                                  as_const(effective_viscosity)) ||
      detail::field_views_overlap(input.molecular_viscosity,
                                  as_const(effective_viscosity))) {
    return {StatusCode::invalid_plan, kTurbulenceView};
  }
  std::uint64_t state = kFnvOffset;
  state = mix(state, fingerprint_);
  const ConstFieldView views[]{input.density, input.molecular_viscosity,
                               as_const(effective_viscosity)};
  for (ConstFieldView view : views) {
    state = mix(state, view.field);
    state = mix(state, view.revision);
    state = mix(state, view.storage_identity);
    state = mix(state, view.revision_domain);
  }
  state = mix(state, input.gradient_revision);
  if (gradient_field) {
    state = mix(state, input.velocity_gradient_field.field);
    state = mix(state, input.velocity_gradient_field.storage_identity);
    state = mix(state, input.velocity_gradient_field.revision_domain);
  }
  state = state == 0U ? 1U : state;
  if (active_.valid() && active_.state == state) {
    certificate = active_;
    return {};
  }

  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells_.z; ++z) {
    for (std::int32_t y = 0; y < cells_.y; ++y) {
      for (std::int32_t x = 0; x < cells_.x; ++x, ++flat) {
        const Int3 cell{x, y, z};
        const double rho = input.density.unchecked(cell, 0U);
        const double molecular =
            input.molecular_viscosity.unchecked(cell, 0U);
        if (!std::isfinite(rho) || !(rho > 0.0) ||
            !std::isfinite(molecular) || !(molecular > 0.0)) {
          return {StatusCode::numerical_failure, kTurbulenceNumerical};
        }
        double kinematic = 0.0;
        VelocityGradient field_gradient;
        const VelocityGradient* gradient = nullptr;
        if (gradient_field) {
          for (std::uint8_t component = 0U; component < 9U; ++component) {
            field_gradient.value[component] =
                input.velocity_gradient_field.unchecked(cell, component);
          }
          gradient = &field_gradient;
        } else {
          gradient = &input.velocity_gradient.data[flat];
        }
        Status evaluated;
        if (subgrid_ == SubgridKind::wale) {
          evaluated = wale_kinematic_viscosity(
              *gradient,
              filter_metrics_[flat].isotropic, coefficient_, kinematic);
        } else if (subgrid_ == SubgridKind::vreman) {
          evaluated = vreman_kinematic_viscosity(
              *gradient,
              {filter_metrics_[flat].x, filter_metrics_[flat].y,
               filter_metrics_[flat].z},
              coefficient_, kinematic);
        }
        const double value = molecular + rho * kinematic;
        if (!evaluated || !std::isfinite(value) || value < molecular) {
          return evaluated ? Status{StatusCode::numerical_failure,
                                    kTurbulenceNumerical}
                           : evaluated;
        }
        pending_effective_viscosity_[flat] = value;
      }
    }
  }
  flat = 0U;
  for (std::int32_t z = 0; z < cells_.z; ++z) {
    for (std::int32_t y = 0; y < cells_.y; ++y) {
      for (std::int32_t x = 0; x < cells_.x; ++x, ++flat) {
        effective_viscosity.unchecked({x, y, z}, 0U) =
            pending_effective_viscosity_[flat];
      }
    }
  }
  TurbulenceCertificate next{fingerprint_,
                             input.density.revision,
                             input.molecular_viscosity.revision,
                             input.gradient_revision,
                             effective_viscosity.revision,
                             state};
  active_ = next;
  certificate = next;
  ++update_count_;
  return {};
}

Status TurbulencePlan::evaluate_candidate_effective_viscosity(
    const TurbulenceCandidateInput& input, FieldView effective_viscosity,
    TurbulenceCandidateCertificate& certificate) const noexcept {
  certificate = {};
  if (fingerprint_ == 0U || !authority_.claimed() ||
      !valid_scalar(input.density, cells_) ||
      !valid_scalar(input.molecular_viscosity, cells_) ||
      !valid_gradient(input.velocity_gradient, cells_) ||
      !valid_scalar(effective_viscosity, cells_) ||
      input.gradient_revision == 0U ||
      input.velocity_gradient.revision != input.gradient_revision ||
      effective_viscosity.field == authority_.output() ||
      detail::field_views_overlap(input.density,
                                  as_const(effective_viscosity)) ||
      detail::field_views_overlap(input.molecular_viscosity,
                                  as_const(effective_viscosity)) ||
      detail::field_views_overlap(input.velocity_gradient,
                                  as_const(effective_viscosity))) {
    return {StatusCode::invalid_plan, kTurbulenceView};
  }

  const auto evaluate_cell = [&](Int3 cell, std::size_t flat,
                                 double& value) noexcept -> Status {
    const double rho = input.density.unchecked(cell, 0U);
    const double molecular =
        input.molecular_viscosity.unchecked(cell, 0U);
    if (!std::isfinite(rho) || !(rho > 0.0) ||
        !std::isfinite(molecular) || !(molecular > 0.0)) {
      return {StatusCode::numerical_failure, kTurbulenceNumerical};
    }
    VelocityGradient gradient;
    for (std::uint8_t component = 0U; component < 9U; ++component) {
      gradient.value[component] =
          input.velocity_gradient.unchecked(cell, component);
    }
    double kinematic = 0.0;
    Status evaluated;
    if (subgrid_ == SubgridKind::wale) {
      evaluated = wale_kinematic_viscosity(
          gradient, filter_metrics_[flat].isotropic, coefficient_, kinematic);
    } else if (subgrid_ == SubgridKind::vreman) {
      evaluated = vreman_kinematic_viscosity(
          gradient,
          {filter_metrics_[flat].x, filter_metrics_[flat].y,
           filter_metrics_[flat].z},
          coefficient_, kinematic);
    }
    value = molecular + rho * kinematic;
    if (!evaluated || !std::isfinite(value) || value < molecular) {
      return evaluated
                 ? Status{StatusCode::numerical_failure,
                          kTurbulenceNumerical}
                 : evaluated;
    }
    return {};
  };

  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells_.z; ++z) {
    for (std::int32_t y = 0; y < cells_.y; ++y) {
      for (std::int32_t x = 0; x < cells_.x; ++x, ++flat) {
        double value = 0.0;
        const Status status = evaluate_cell({x, y, z}, flat, value);
        if (!status) {
          return status;
        }
      }
    }
  }

  const TurbulenceCandidateFieldBinding density =
      candidate_binding(input.density);
  const TurbulenceCandidateFieldBinding molecular =
      candidate_binding(input.molecular_viscosity);
  const TurbulenceCandidateFieldBinding gradient =
      candidate_binding(input.velocity_gradient);
  const TurbulenceCandidateFieldBinding output =
      candidate_binding(as_const(effective_viscosity));
  std::uint64_t state = mix(kFnvOffset, fingerprint_);
  state = mix_candidate_binding(state, density);
  state = mix_candidate_binding(state, molecular);
  state = mix_candidate_binding(state, gradient);
  state = mix_candidate_binding(state, output);
  state = mix(state, input.gradient_revision);
  state = mix(state, authority_.output());
  state = state == 0U ? 1U : state;
  const TurbulenceCandidateCertificate next{
      fingerprint_, density, molecular, gradient, output, authority_.output(),
      state, true, true, true};
  if (!next.valid() ||
      !next.matches(fingerprint_, input, as_const(effective_viscosity))) {
    return {StatusCode::invalid_plan, kTurbulenceView};
  }

  flat = 0U;
  for (std::int32_t z = 0; z < cells_.z; ++z) {
    for (std::int32_t y = 0; y < cells_.y; ++y) {
      for (std::int32_t x = 0; x < cells_.x; ++x, ++flat) {
        double value = 0.0;
        const Status status = evaluate_cell({x, y, z}, flat, value);
        if (!status) {
          return status;
        }
        effective_viscosity.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
  certificate = next;
  return {};
}

}  // namespace hundun::v04
