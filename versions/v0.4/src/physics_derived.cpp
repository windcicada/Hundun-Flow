// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kDerivedInput = 2201U;
constexpr std::uint32_t kDerivedShape = 2202U;
constexpr std::uint32_t kDerivedRevision = 2203U;
constexpr std::uint32_t kDerivedTransaction = 2204U;
constexpr std::uint32_t kDerivedNumerical = 2205U;

bool same_tuple(DerivedRevisionTuple left,
                DerivedRevisionTuple right) noexcept {
  return left.velocity == right.velocity && left.geometry == right.geometry &&
         left.boundary == right.boundary &&
         left.turbulence == right.turbulence &&
         left.velocity_storage == right.velocity_storage &&
         left.velocity_revision_domain == right.velocity_revision_domain &&
         left.patch_identity == right.patch_identity &&
         left.geometry_plan_identity == right.geometry_plan_identity &&
         left.boundary_plan_identity == right.boundary_plan_identity &&
         left.turbulence_plan_identity == right.turbulence_plan_identity;
}

bool valid_tuple(DerivedRevisionTuple value) noexcept {
  return value.velocity != 0U && value.geometry != 0U &&
         value.boundary != 0U && value.turbulence != 0U &&
         value.velocity_storage != 0U &&
         value.velocity_revision_domain != 0U &&
         value.patch_identity != 0U && value.geometry_plan_identity != 0U &&
         value.boundary_plan_identity != 0U &&
         value.turbulence_plan_identity != 0U;
}

PlanFingerprint patch_identity(const MeshPatch& patch,
                               PlanFingerprint geometry) noexcept {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  auto mix = [&](std::uint64_t value) noexcept {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
  };
  mix(geometry);
  mix(static_cast<std::uint32_t>(patch.begin.x));
  mix(static_cast<std::uint32_t>(patch.begin.y));
  mix(static_cast<std::uint32_t>(patch.begin.z));
  mix(static_cast<std::uint32_t>(patch.cells.x));
  mix(static_cast<std::uint32_t>(patch.cells.y));
  mix(static_cast<std::uint32_t>(patch.cells.z));
  return hash == 0U ? 1U : hash;
}

bool valid_field_view(ConstFieldView view, FieldId expected,
                      std::size_t capacity) noexcept {
  if (view.base == nullptr || view.field != expected || view.components != 3U ||
      view.interior.x <= 0 || view.interior.y <= 0 || view.interior.z <= 0 ||
      view.ghosts.x < 1 || view.ghosts.y < 1 || view.ghosts.z < 1 ||
      view.stride_y == 0U || view.stride_z == 0U ||
      view.component_stride == 0U || view.revision == 0U ||
      view.storage_identity == 0U || view.revision_domain == 0U) {
    return false;
  }
  const auto x = static_cast<std::size_t>(view.interior.x);
  const auto y = static_cast<std::size_t>(view.interior.y);
  const auto z = static_cast<std::size_t>(view.interior.z);
  return x <= std::numeric_limits<std::size_t>::max() / y &&
         x * y <= std::numeric_limits<std::size_t>::max() / z &&
         x * y * z == capacity;
}

bool valid_patch(const MeshPatch& patch, ConstFieldView view,
                 const CartesianGeometryPlan& geometry) noexcept {
  const Int3 global = geometry.global_cells();
  return patch.cells.x == view.interior.x && patch.cells.y == view.interior.y &&
         patch.cells.z == view.interior.z && patch.begin.x >= 0 &&
         patch.begin.y >= 0 && patch.begin.z >= 0 && patch.cells.x > 0 &&
         patch.cells.y > 0 && patch.cells.z > 0 &&
         patch.begin.x <= global.x - patch.cells.x &&
         patch.begin.y <= global.y - patch.cells.y &&
         patch.begin.z <= global.z - patch.cells.z;
}

struct DerivativeWeights {
  double minus{};
  double centre{};
  double plus{};
};

DerivativeWeights derivative_weights(const AxisMetrics& axis,
                                     std::int32_t global_index) noexcept {
  if (axis.uniform()) {
    const double half_inverse_width = 0.5 * axis.uniform_inverse_width();
    return {-half_inverse_width, 0.0, half_inverse_width};
  }

  const Span<const double> centres = axis.centres();
  const Span<const double> widths = axis.widths();
  const std::size_t index = static_cast<std::size_t>(global_index);
  const double x_centre = centres.data[index];
  const double x_minus =
      global_index == 0 ? x_centre - widths.data[0U]
                        : centres.data[index - 1U];
  const double x_plus =
      index + 1U == centres.size
          ? x_centre + widths.data[centres.size - 1U]
          : centres.data[index + 1U];
  const double left_distance = x_centre - x_minus;
  const double right_distance = x_plus - x_centre;
  const double total_distance = left_distance + right_distance;
  return {-right_distance / (left_distance * total_distance),
          (right_distance - left_distance) /
              (left_distance * right_distance),
          left_distance / (right_distance * total_distance)};
}

double apply_derivative(ConstFieldView velocity, Int3 minus, Int3 centre,
                        Int3 plus, std::uint8_t component,
                        DerivativeWeights weights) noexcept {
  return weights.minus * velocity.unchecked(minus, component) +
         weights.centre * velocity.unchecked(centre, component) +
         weights.plus * velocity.unchecked(plus, component);
}

}  // namespace

DerivedRevisionTuple make_derived_revision_tuple(
    ConstFieldView velocity, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, RevisionToken boundary,
    PlanFingerprint boundary_plan_identity, RevisionToken turbulence,
    PlanFingerprint turbulence_plan_identity) noexcept {
  return {velocity.revision,
          geometry.topology_revision(),
          boundary,
          turbulence,
          velocity.storage_identity,
          velocity.revision_domain,
          patch_identity(patch, geometry.fingerprint()),
          geometry.fingerprint(),
          boundary_plan_identity,
          turbulence_plan_identity};
}

Status DerivedFieldPlan::compile(
    FieldId velocity_field, Span<const FieldId> declared_fields,
    RevisionSlotId gradient_cache_slot, RevisionSourceId geometry_source,
    RevisionSourceId boundary_source, RevisionSourceId turbulence_source,
    std::size_t cell_capacity, DerivedFieldPlan& out) noexcept {
  if (declared_fields.data == nullptr || declared_fields.size == 0U ||
      std::find(declared_fields.data,
                declared_fields.data + declared_fields.size,
                velocity_field) == declared_fields.data +
                                       declared_fields.size ||
      geometry_source == 0U || boundary_source == 0U ||
      turbulence_source == 0U ||
      geometry_source == boundary_source ||
      geometry_source == turbulence_source ||
      boundary_source == turbulence_source ||
      cell_capacity == 0U) {
    return {StatusCode::invalid_plan, kDerivedInput};
  }
  for (std::size_t first = 0U; first < declared_fields.size; ++first) {
    const RevisionSourceId field_source =
        AttemptTransaction::field_revision_source(declared_fields.data[first]);
    if (geometry_source == field_source || boundary_source == field_source ||
        turbulence_source == field_source) {
      return {StatusCode::invalid_plan, kDerivedInput};
    }
    for (std::size_t second = first + 1U; second < declared_fields.size;
         ++second) {
      if (declared_fields.data[first] == declared_fields.data[second]) {
        return {StatusCode::invalid_plan, kDerivedInput};
      }
    }
  }
  try {
    DerivedFieldPlan candidate;
    candidate.active_gradient_.resize(cell_capacity);
    candidate.pending_gradient_.resize(cell_capacity);
    candidate.velocity_field_ = velocity_field;
    candidate.cache_slot_ = gradient_cache_slot;
    candidate.geometry_source_ = geometry_source;
    candidate.boundary_source_ = boundary_source;
    candidate.turbulence_source_ = turbulence_source;
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kDerivedInput};
  } catch (...) {
    return {StatusCode::invalid_plan, kDerivedInput};
  }
}

Status DerivedFieldPlan::prepare_velocity_gradient(
    ConstFieldView velocity, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, DerivedRevisionTuple revisions,
    AttemptTransaction& transaction,
    Span<const VelocityGradient>& out) noexcept {
  out = {};
  if (!valid_tuple(revisions) || revisions.velocity != velocity.revision ||
      revisions.geometry != geometry.topology_revision() ||
      revisions.velocity_storage != velocity.storage_identity ||
      revisions.velocity_revision_domain != velocity.revision_domain ||
      revisions.patch_identity != patch_identity(patch, geometry.fingerprint()) ||
      revisions.geometry_plan_identity != geometry.fingerprint() ||
      !valid_field_view(velocity, velocity_field_, pending_gradient_.size()) ||
      !valid_patch(patch, velocity, geometry) || !transaction.active() ||
      transaction.attempt_identity() == 0U) {
    return {StatusCode::invalid_plan, kDerivedShape};
  }
  if (pending_valid_) {
    if (!same_tuple(revisions, pending_revisions_) ||
        transaction.attempt_identity() != pending_attempt_identity_) {
      return {StatusCode::invalid_plan, kDerivedTransaction};
    }
    out = {pending_gradient_.data(), pending_gradient_.size()};
    return {};
  }
  if (active_valid_ && same_tuple(revisions, active_revisions_)) {
    out = {active_gradient_.data(), active_gradient_.size()};
    return {};
  }

  const std::array dependencies{
      RevisionDependency{
          AttemptTransaction::field_revision_source(velocity_field_),
          revisions.velocity},
      RevisionDependency{geometry_source_, revisions.geometry},
      RevisionDependency{boundary_source_, revisions.boundary},
      RevisionDependency{turbulence_source_, revisions.turbulence}};
  if (next_cache_revision_ == 0U ||
      next_cache_revision_ == std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kDerivedRevision};
  }

  const AxisMetrics& x_axis = geometry.x();
  const AxisMetrics& y_axis = geometry.y();
  const AxisMetrics& z_axis = geometry.z();
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < velocity.interior.z; ++z) {
    const DerivativeWeights z_weights =
        derivative_weights(z_axis, patch.begin.z + z);
    for (std::int32_t y = 0; y < velocity.interior.y; ++y) {
      const DerivativeWeights y_weights =
          derivative_weights(y_axis, patch.begin.y + y);
      for (std::int32_t x = 0; x < velocity.interior.x; ++x, ++flat) {
        const DerivativeWeights x_weights =
            derivative_weights(x_axis, patch.begin.x + x);
        VelocityGradient gradient;
        const Int3 local{x, y, z};
        const Int3 x_minus{x - 1, y, z};
        const Int3 x_plus{x + 1, y, z};
        const Int3 y_minus{x, y - 1, z};
        const Int3 y_plus{x, y + 1, z};
        const Int3 z_minus{x, y, z - 1};
        const Int3 z_plus{x, y, z + 1};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const std::size_t offset =
              static_cast<std::size_t>(component) * 3U;
          gradient.value[offset] = apply_derivative(
              velocity, x_minus, local, x_plus, component, x_weights);
          gradient.value[offset + 1U] = apply_derivative(
              velocity, y_minus, local, y_plus, component, y_weights);
          gradient.value[offset + 2U] = apply_derivative(
              velocity, z_minus, local, z_plus, component, z_weights);
          if (!std::isfinite(gradient.value[offset]) ||
              !std::isfinite(gradient.value[offset + 1U]) ||
              !std::isfinite(gradient.value[offset + 2U])) {
            return {StatusCode::numerical_failure, kDerivedNumerical};
          }
        }
        pending_gradient_[flat] = gradient;
      }
    }
  }
  Status bound = transaction.bind_dependency(dependencies[1]);
  if (bound) {
    bound = transaction.bind_dependency(dependencies[2]);
  }
  if (bound) {
    bound = transaction.bind_dependency(dependencies[3]);
  }
  if (!bound) {
    return bound;
  }
  const RevisionToken candidate_revision = next_cache_revision_;
  const Status published = transaction.publish_pending_cache(
      cache_slot_,
      Span<const RevisionDependency>{dependencies.data(), dependencies.size()},
      PendingCacheStamp{candidate_revision});
  if (!published) {
    return published;
  }
  pending_revisions_ = revisions;
  pending_cache_revision_ = candidate_revision;
  pending_attempt_identity_ = transaction.attempt_identity();
  pending_valid_ = true;
  ++next_cache_revision_;
  ++gradient_compute_count_;
  out = {pending_gradient_.data(), pending_gradient_.size()};
  return {};
}

Status DerivedFieldPlan::finalize_velocity_gradient(
    const AttemptTransaction& transaction) noexcept {
  if (!pending_valid_) {
    return {StatusCode::invalid_plan, kDerivedTransaction};
  }
  if (!transaction.finished() || transaction.active() ||
      transaction.attempt_identity() != pending_attempt_identity_) {
    return {StatusCode::invalid_plan, kDerivedTransaction};
  }
  if (transaction.committed()) {
    if (!transaction.pending_cache_valid(cache_slot_) ||
        transaction.pending_cache(cache_slot_) != pending_cache_revision_) {
      pending_valid_ = false;
      pending_cache_revision_ = 0U;
      pending_attempt_identity_ = 0U;
      return {StatusCode::invalid_plan, kDerivedTransaction};
    }
    active_gradient_.swap(pending_gradient_);
    active_revisions_ = pending_revisions_;
    cache_revision_ = pending_cache_revision_;
    active_valid_ = true;
  }
  pending_valid_ = false;
  pending_cache_revision_ = 0U;
  pending_attempt_identity_ = 0U;
  pending_revisions_ = {};
  return {};
}

Status DerivedFieldPlan::velocity_gradient(
    DerivedRevisionTuple revisions,
    Span<const VelocityGradient>& out) const noexcept {
  out = {};
  if (!active_valid_ || !valid_tuple(revisions) ||
      !same_tuple(revisions, active_revisions_) || cache_revision_ == 0U) {
    return {StatusCode::invalid_plan, kDerivedRevision};
  }
  out = {active_gradient_.data(), active_gradient_.size()};
  return {};
}

}  // namespace hundun::v04
