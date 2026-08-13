// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kApplyStage = 301U;
constexpr std::uint32_t kApplyField = 302U;
constexpr std::uint32_t kApplyView = 303U;
constexpr std::uint32_t kApplyParameter = 304U;

constexpr std::size_t face_index(CartesianFace face) noexcept {
  return static_cast<std::size_t>(face);
}

FieldView* find_field(Span<FieldView> fields, FieldId id) noexcept {
  for (std::size_t index = 0U; index < fields.size; ++index) {
    if (fields.data[index].field == id) {
      return fields.data + index;
    }
  }
  return nullptr;
}

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& out) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

bool checked_ghosted_extent(std::int32_t interior, std::int32_t ghosts,
                            std::size_t& out) noexcept {
  if (interior <= 0 || ghosts < 0) {
    return false;
  }
  std::size_t twice_ghosts = 0U;
  return checked_multiply(static_cast<std::size_t>(ghosts), 2U,
                          twice_ghosts) &&
         checked_add(static_cast<std::size_t>(interior), twice_ghosts, out);
}

bool checked_affine_offset(const FieldView& view,
                           std::size_t extent_x,
                           std::size_t extent_y,
                           std::size_t extent_z) noexcept {
  const std::size_t ptrdiff_max = static_cast<std::size_t>(
      std::numeric_limits<std::ptrdiff_t>::max());
  if (view.stride_y > ptrdiff_max || view.stride_z > ptrdiff_max ||
      view.component_stride > ptrdiff_max || extent_x == 0U ||
      extent_y == 0U || extent_z == 0U || view.components == 0U) {
    return false;
  }
  std::size_t offset_y = 0U;
  std::size_t offset_z = 0U;
  std::size_t offset_component = 0U;
  std::size_t maximum_offset = extent_x - 1U;
  return checked_multiply(extent_y - 1U, view.stride_y, offset_y) &&
         checked_add(maximum_offset, offset_y, maximum_offset) &&
         checked_multiply(extent_z - 1U, view.stride_z, offset_z) &&
         checked_add(maximum_offset, offset_z, maximum_offset) &&
         checked_multiply(static_cast<std::size_t>(view.components - 1U),
                          view.component_stride, offset_component) &&
         checked_add(maximum_offset, offset_component, maximum_offset) &&
         maximum_offset <= ptrdiff_max;
}

bool valid_view(const FieldView& view,
                const BoundaryIndexSpan& span,
                Int3 expected_cells) noexcept {
  if (view.base == nullptr || view.interior.x <= 0 || view.interior.y <= 0 ||
      view.interior.z <= 0 || view.ghosts.x < 0 || view.ghosts.y < 0 ||
      view.ghosts.z < 0 || span.ghost_layers == 0U ||
      span.component_count == 0U ||
      static_cast<unsigned int>(span.component_begin) +
              static_cast<unsigned int>(span.component_count) >
          static_cast<unsigned int>(view.components)) {
    return false;
  }
  if (view.interior.x != expected_cells.x ||
      view.interior.y != expected_cells.y ||
      view.interior.z != expected_cells.z) {
    return false;
  }

  const auto layers = static_cast<std::int32_t>(span.ghost_layers);
  const std::size_t face = face_index(span.face);
  const std::int32_t normal_interior =
      face < 2U ? view.interior.x
                : (face < 4U ? view.interior.y : view.interior.z);
  const std::int32_t normal_ghosts =
      face < 2U ? view.ghosts.x
                : (face < 4U ? view.ghosts.y : view.ghosts.z);
  if (normal_interior < layers || normal_ghosts < layers) {
    return false;
  }

  const std::uint32_t expected_inner =
      face < 2U
          ? static_cast<std::uint32_t>(view.interior.y)
          : static_cast<std::uint32_t>(view.interior.x);
  const std::uint32_t expected_outer =
      face < 4U
          ? static_cast<std::uint32_t>(view.interior.z)
          : static_cast<std::uint32_t>(view.interior.y);
  if (span.tangent_inner_count != expected_inner ||
      span.tangent_outer_count != expected_outer) {
    return false;
  }

  std::size_t extent_x = 0U;
  std::size_t extent_y = 0U;
  std::size_t extent_z = 0U;
  std::size_t minimum_stride_z = 0U;
  std::size_t minimum_component_stride = 0U;
  return checked_ghosted_extent(view.interior.x, view.ghosts.x, extent_x) &&
         checked_ghosted_extent(view.interior.y, view.ghosts.y, extent_y) &&
         checked_ghosted_extent(view.interior.z, view.ghosts.z, extent_z) &&
         view.stride_y >= extent_x &&
         checked_multiply(view.stride_y, extent_y, minimum_stride_z) &&
         view.stride_z >= minimum_stride_z &&
         checked_multiply(view.stride_z, extent_z,
                          minimum_component_stride) &&
         view.component_stride >= minimum_component_stride &&
         checked_affine_offset(view, extent_x, extent_y, extent_z);
}

bool valid_parameter_storage(const BoundaryPlan& plan) noexcept {
  const std::size_t count = plan.parameter_count();
  return plan.velocity_x().size == count &&
         plan.velocity_y().size == count &&
         plan.velocity_z().size == count &&
         plan.direction_x().size == count &&
         plan.direction_y().size == count &&
         plan.direction_z().size == count &&
         plan.backflow_velocity_x().size == count &&
         plan.backflow_velocity_y().size == count &&
         plan.backflow_velocity_z().size == count &&
         plan.scalar_targets().size == count &&
         plan.pressure_targets().size == count &&
         plan.temperature_targets().size == count &&
         plan.total_pressure_targets().size == count &&
         plan.total_temperature_targets().size == count &&
         plan.backflow_temperature_targets().size == count &&
         plan.heat_flux_targets().size == count &&
         plan.mass_flow_targets().size == count &&
         plan.relaxation_rates().size == count &&
         plan.mach_limits().size == count &&
         plan.allow_backflow().size == count &&
         plan.parameter_roles().size == count &&
         plan.scalar_kinds().size == count &&
         plan.scalar_roles().size == count &&
         plan.normal_distance_1().size == count &&
         plan.normal_distance_2().size == count;
}

bool finite(double value) noexcept { return std::isfinite(value); }

bool finite(Real3 value) noexcept {
  return finite(value.x) && finite(value.y) && finite(value.z);
}

template <class T, class Finite>
bool valid_resolved_slice(const BoundaryIndexSpan& span,
                          std::size_t plan_count, Span<const T> values,
                          Finite is_finite) noexcept {
  std::size_t face_cells = 0U;
  std::size_t end = 0U;
  if (span.resolved_begin == kInvalidBoundaryParameter ||
      !checked_multiply(static_cast<std::size_t>(span.tangent_inner_count),
                        static_cast<std::size_t>(span.tangent_outer_count),
                        face_cells) ||
      face_cells == 0U ||
      face_cells > std::numeric_limits<std::uint32_t>::max() ||
      span.resolved_stride != face_cells ||
      !checked_add(static_cast<std::size_t>(span.resolved_begin), face_cells,
                   end) ||
      end > plan_count || end > values.size || values.data == nullptr) {
    return false;
  }
  for (std::size_t index = span.resolved_begin; index < end; ++index) {
    if (!is_finite(values.data[index])) {
      return false;
    }
  }
  return true;
}

bool valid_relation_source(const BoundaryIndexSpan& span) noexcept {
  switch (span.relation) {
    case BoundaryRelation::dirichlet:
      return (span.component_count == 1U &&
              (span.value_source == BoundaryValueSource::compiled_scalar ||
               span.value_source == BoundaryValueSource::resolved_scalar)) ||
             (span.component_count == 3U &&
              (span.value_source == BoundaryValueSource::compiled_vector ||
               span.value_source == BoundaryValueSource::resolved_vector));
    case BoundaryRelation::zero_gradient:
    case BoundaryRelation::reflect_normal:
      return span.value_source == BoundaryValueSource::none;
    case BoundaryRelation::normal_gradient:
      return span.component_count == 1U &&
             span.value_source ==
                 BoundaryValueSource::resolved_normal_gradient;
    case BoundaryRelation::convective:
      return span.component_count == 1U &&
             span.value_source == BoundaryValueSource::resolved_scalar;
  }
  return false;
}

bool valid_stage_plan(BoundaryStage stage, const BoundaryPlan& plan,
                      Span<FieldView> fields,
                      BoundaryResolvedValues values) noexcept {
  if (static_cast<std::uint8_t>(stage) >
      static_cast<std::uint8_t>(BoundaryStage::diagnostic)) {
    return false;
  }
  const Span<const BoundaryIndexSpan> spans = plan.spans();
  const Span<const BoundaryKernelBatch> batches = plan.batches();
  std::size_t prior_end = 0U;
  for (std::size_t batch_index = 0U; batch_index < batches.size;
       ++batch_index) {
    const BoundaryKernelBatch& batch = batches.data[batch_index];
    const std::size_t begin = batch.span_begin;
    const std::size_t count = batch.span_count;
    std::size_t end = 0U;
    if (static_cast<std::uint8_t>(batch.stage) >
            static_cast<std::uint8_t>(BoundaryStage::diagnostic) ||
        static_cast<std::uint8_t>(batch.relation) >
            static_cast<std::uint8_t>(BoundaryRelation::convective)) {
      return false;
    }
    if (batch.stage != stage || count == 0U) {
      continue;
    }
    if (begin < prior_end || !checked_add(begin, count, end) ||
        begin > spans.size || end > spans.size) {
      return false;
    }
    prior_end = end;
    for (std::size_t index = begin; index < end; ++index) {
      const BoundaryIndexSpan& span = spans.data[index];
      if (span.stage != batch.stage || span.relation != batch.relation ||
          static_cast<std::uint8_t>(span.face) >
              static_cast<std::uint8_t>(CartesianFace::z_max) ||
          span.parameter >= plan.parameter_count() ||
          !valid_relation_source(span)) {
        return false;
      }
      switch (span.value_source) {
        case BoundaryValueSource::none:
        case BoundaryValueSource::compiled_scalar:
        case BoundaryValueSource::compiled_vector:
          if (span.resolved_begin != kInvalidBoundaryParameter ||
              span.resolved_stride != 0U) {
            return false;
          }
          if (span.value_source == BoundaryValueSource::compiled_scalar &&
              !finite(plan.scalar_targets().data[span.parameter])) {
            return false;
          }
          if (span.value_source == BoundaryValueSource::compiled_vector &&
              (!finite(plan.velocity_x().data[span.parameter]) ||
               !finite(plan.velocity_y().data[span.parameter]) ||
               !finite(plan.velocity_z().data[span.parameter]))) {
            return false;
          }
          break;
        case BoundaryValueSource::resolved_scalar:
          if (!valid_resolved_slice(
                  span, plan.resolved_scalar_count(), values.scalar,
                  [](double value) noexcept { return finite(value); })) {
            return false;
          }
          break;
        case BoundaryValueSource::resolved_vector:
          if (!valid_resolved_slice(
                  span, plan.resolved_vector_count(), values.vector,
                  [](Real3 value) noexcept { return finite(value); })) {
            return false;
          }
          break;
        case BoundaryValueSource::resolved_normal_gradient:
          if (!valid_resolved_slice(
                  span, plan.resolved_normal_gradient_count(),
                  values.normal_gradient,
                  [](double value) noexcept { return finite(value); })) {
            return false;
          }
          break;
        default:
          return false;
      }
      FieldView* view = find_field(fields, span.field);
      if (view == nullptr || !valid_view(*view, span, plan.local_cells())) {
        return false;
      }
    }
  }
  return true;
}

double resolved_scalar(const BoundaryPlan& plan,
                       BoundaryResolvedValues values,
                       const BoundaryIndexSpan& span,
                       std::size_t face_cell, bool& ok) noexcept {
  const std::uint32_t parameter = span.parameter;
  if (parameter >= plan.parameter_count()) {
    ok = false;
    return 0.0;
  }
  if (span.value_source == BoundaryValueSource::resolved_scalar) {
    std::size_t index = 0U;
    if (!checked_add(static_cast<std::size_t>(span.resolved_begin), face_cell,
                     index) ||
        face_cell >= span.resolved_stride || index >= values.scalar.size) {
      ok = false;
      return 0.0;
    }
    return values.scalar.data[index];
  }
  if (span.value_source == BoundaryValueSource::compiled_scalar) {
    return plan.scalar_targets().data[parameter];
  }
  return 0.0;
}

Real3 resolved_vector(const BoundaryPlan& plan,
                      BoundaryResolvedValues values,
                      const BoundaryIndexSpan& span,
                      std::size_t face_cell, bool& ok) noexcept {
  const std::uint32_t parameter = span.parameter;
  if (parameter >= plan.parameter_count()) {
    ok = false;
    return {};
  }
  if (span.value_source == BoundaryValueSource::resolved_vector) {
    std::size_t index = 0U;
    if (!checked_add(static_cast<std::size_t>(span.resolved_begin), face_cell,
                     index) ||
        face_cell >= span.resolved_stride || index >= values.vector.size) {
      ok = false;
      return {};
    }
    return values.vector.data[index];
  }
  if (span.value_source != BoundaryValueSource::compiled_vector) {
    return {};
  }
  return Real3{plan.velocity_x().data[parameter],
               plan.velocity_y().data[parameter],
               plan.velocity_z().data[parameter]};
}

double resolved_gradient(const BoundaryPlan& plan,
                         BoundaryResolvedValues values,
                         const BoundaryIndexSpan& span,
                         std::size_t face_cell, bool& ok) noexcept {
  const std::uint32_t parameter = span.parameter;
  if (parameter >= plan.parameter_count()) {
    ok = false;
    return 0.0;
  }
  std::size_t index = 0U;
  if (span.value_source != BoundaryValueSource::resolved_normal_gradient ||
      !checked_add(static_cast<std::size_t>(span.resolved_begin), face_cell,
                   index) ||
      face_cell >= span.resolved_stride ||
      index >= values.normal_gradient.size) {
    ok = false;
    return 0.0;
  }
  return values.normal_gradient.data[index];
}

double normal_distance(const BoundaryPlan& plan, std::uint32_t parameter,
                       std::int32_t layer, bool& ok) noexcept {
  const Span<const double> distances =
      layer == 1 ? plan.normal_distance_1() : plan.normal_distance_2();
  if (parameter >= distances.size) {
    ok = false;
    return 0.0;
  }
  return distances.data[parameter];
}

double vector_component(Real3 value, std::uint8_t component) noexcept {
  return component == 0U ? value.x : (component == 1U ? value.y : value.z);
}

Int3 interior_index(CartesianFace face, const FieldView& view,
                    std::int32_t layer, std::int32_t inner,
                    std::int32_t outer) noexcept {
  switch (face) {
    case CartesianFace::x_min:
      return Int3{layer - 1, inner, outer};
    case CartesianFace::x_max:
      return Int3{view.interior.x - layer, inner, outer};
    case CartesianFace::y_min:
      return Int3{inner, layer - 1, outer};
    case CartesianFace::y_max:
      return Int3{inner, view.interior.y - layer, outer};
    case CartesianFace::z_min:
      return Int3{inner, outer, layer - 1};
    case CartesianFace::z_max:
      return Int3{inner, outer, view.interior.z - layer};
  }
  return {};
}

Int3 ghost_index(CartesianFace face, const FieldView& view,
                 std::int32_t layer, std::int32_t inner,
                 std::int32_t outer) noexcept {
  switch (face) {
    case CartesianFace::x_min:
      return Int3{-layer, inner, outer};
    case CartesianFace::x_max:
      return Int3{view.interior.x - 1 + layer, inner, outer};
    case CartesianFace::y_min:
      return Int3{inner, -layer, outer};
    case CartesianFace::y_max:
      return Int3{inner, view.interior.y - 1 + layer, outer};
    case CartesianFace::z_min:
      return Int3{inner, outer, -layer};
    case CartesianFace::z_max:
      return Int3{inner, outer, view.interior.z - 1 + layer};
  }
  return {};
}

template <BoundaryRelation Relation>
Status apply_span(const BoundaryPlan& plan, const BoundaryIndexSpan& span,
                  FieldView& view, BoundaryResolvedValues values) noexcept {
  const auto inner_count =
      static_cast<std::int32_t>(span.tangent_inner_count);
  const auto outer_count =
      static_cast<std::int32_t>(span.tangent_outer_count);
  for (std::int32_t outer = 0; outer < outer_count; ++outer) {
    for (std::int32_t inner = 0; inner < inner_count; ++inner) {
      const std::size_t face_cell =
          static_cast<std::size_t>(outer) * span.tangent_inner_count +
          static_cast<std::size_t>(inner);
      bool parameter_ok = true;
      double scalar = 0.0;
      Real3 vector{};
      double gradient = 0.0;
      if constexpr (Relation == BoundaryRelation::dirichlet ||
                    Relation == BoundaryRelation::convective) {
        if (span.component_count == 3U) {
          vector =
              resolved_vector(plan, values, span, face_cell, parameter_ok);
        } else {
          scalar =
              resolved_scalar(plan, values, span, face_cell, parameter_ok);
        }
      } else if constexpr (Relation == BoundaryRelation::normal_gradient) {
        gradient =
            resolved_gradient(plan, values, span, face_cell, parameter_ok);
      }
      if (!parameter_ok) {
        return {StatusCode::invalid_plan, kApplyParameter};
      }
      for (std::int32_t layer = 1;
           layer <= static_cast<std::int32_t>(span.ghost_layers); ++layer) {
        const Int3 source =
            interior_index(span.face, view, layer, inner, outer);
        const Int3 destination =
            ghost_index(span.face, view, layer, inner, outer);
        for (std::uint8_t offset = 0U; offset < span.component_count;
             ++offset) {
          const std::uint8_t component =
              static_cast<std::uint8_t>(span.component_begin + offset);
          const double interior = view.unchecked(source, component);
          double result = interior;
          if constexpr (Relation == BoundaryRelation::dirichlet) {
            const double face_value = span.component_count == 3U
                                          ? vector_component(vector, component)
                                          : scalar;
            result = 2.0 * face_value - interior;
          } else if constexpr (Relation == BoundaryRelation::normal_gradient) {
            const double distance =
                normal_distance(plan, span.parameter, layer, parameter_ok);
            if (!parameter_ok) {
              return {StatusCode::invalid_plan, kApplyParameter};
            }
            result = interior + gradient * distance;
          } else if constexpr (Relation == BoundaryRelation::reflect_normal) {
            const std::uint8_t normal_component =
                face_index(span.face) < 2U
                    ? 0U
                    : (face_index(span.face) < 4U ? 1U : 2U);
            result = component == normal_component ? -interior : interior;
          } else if constexpr (Relation == BoundaryRelation::convective) {
            // The transport solver supplies its already-discretized
            // field-space convective target; ghost application remains a
            // branch-free Dirichlet mirror over the compiled batch.
            result = 2.0 * scalar - interior;
          }
          view.unchecked(destination, component) = result;
        }
      }
    }
  }
  return {};
}

template <BoundaryRelation Relation>
Status apply_batch(const BoundaryPlan& plan, const BoundaryKernelBatch& batch,
                   Span<FieldView> fields,
                   BoundaryResolvedValues values) noexcept {
  const Span<const BoundaryIndexSpan> spans = plan.spans();
  const std::size_t begin = batch.span_begin;
  const std::size_t end = begin + batch.span_count;
  if (begin > spans.size || end > spans.size) {
    return {StatusCode::invalid_plan, kApplyStage};
  }
  for (std::size_t index = begin; index < end; ++index) {
    const BoundaryIndexSpan& span = spans.data[index];
    FieldView* view = find_field(fields, span.field);
    if (view == nullptr) {
      return {StatusCode::invalid_plan, kApplyField};
    }
    if (!valid_view(*view, span, plan.local_cells())) {
      return {StatusCode::invalid_plan, kApplyView};
    }
    const Status status = apply_span<Relation>(plan, span, *view, values);
    if (!status) {
      return status;
    }
  }
  return {};
}

}  // namespace

Status BoundaryPlan::face(CartesianFace selected,
                          const BoundaryFacePlan*& out) const noexcept {
  const std::size_t index = face_index(selected);
  if (index >= faces_.size()) {
    out = nullptr;
    return {StatusCode::invalid_plan, kApplyParameter};
  }
  out = &faces_[index];
  return {};
}

Status apply_boundary_ghosts(BoundaryStage stage, const BoundaryPlan& plan,
                             Span<FieldView> fields,
                             BoundaryResolvedValues values) noexcept {
  if ((fields.size != 0U && fields.data == nullptr) ||
      !valid_parameter_storage(plan)) {
    return {StatusCode::invalid_plan, kApplyView};
  }
  // This complete stage preflight is deliberately separate from mutation.
  // A later invalid span must never leave an earlier field half-updated.
  if (!valid_stage_plan(stage, plan, fields, values)) {
    return {StatusCode::invalid_plan, kApplyView};
  }
  const Span<const BoundaryKernelBatch> batches = plan.batches();
  for (std::size_t index = 0U; index < batches.size; ++index) {
    const BoundaryKernelBatch& batch = batches.data[index];
    if (batch.stage != stage) {
      continue;
    }
    Status status;
    switch (batch.relation) {
      case BoundaryRelation::dirichlet:
        status = apply_batch<BoundaryRelation::dirichlet>(
            plan, batch, fields, values);
        break;
      case BoundaryRelation::zero_gradient:
        status = apply_batch<BoundaryRelation::zero_gradient>(
            plan, batch, fields, values);
        break;
      case BoundaryRelation::normal_gradient:
        status = apply_batch<BoundaryRelation::normal_gradient>(
            plan, batch, fields, values);
        break;
      case BoundaryRelation::reflect_normal:
        status = apply_batch<BoundaryRelation::reflect_normal>(
            plan, batch, fields, values);
        break;
      case BoundaryRelation::convective:
        status = apply_batch<BoundaryRelation::convective>(
            plan, batch, fields, values);
        break;
    }
    if (!status) {
      return status;
    }
  }
  return {};
}

}  // namespace hundun::v04
