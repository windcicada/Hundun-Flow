// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_mesh.hpp"

#include "core_arena_detail.hpp"
#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kKernelPlan = 901U;
constexpr std::uint32_t kKernelInvocation = 902U;
constexpr std::uint32_t kKernelNumerical = 903U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

template <class T>
bool component_byte_interval(BasicFieldView<T> view,
                             std::uint8_t component_begin,
                             std::uint8_t component_count,
                             std::uintptr_t& begin,
                             std::uintptr_t& end) noexcept {
  const unsigned component_end = static_cast<unsigned>(component_begin) +
                                 static_cast<unsigned>(component_count);
  detail::FieldStorageInterval storage_interval{};
  if (!detail::field_storage_interval(view, storage_interval) ||
      component_count == 0U ||
      component_end > view.components || view.stride_y == 0U ||
      view.stride_z == 0U || view.component_stride == 0U) {
    return false;
  }

  std::size_t prefix_y = 0U;
  std::size_t prefix_z = 0U;
  std::size_t prefix = 0U;
  std::size_t first_component = 0U;
  if (!detail::checked_multiply(static_cast<std::size_t>(view.ghosts.y),
                                view.stride_y, prefix_y) ||
      !detail::checked_multiply(static_cast<std::size_t>(view.ghosts.z),
                                view.stride_z, prefix_z) ||
      !detail::checked_add(static_cast<std::size_t>(view.ghosts.x), prefix_y,
                           prefix) ||
      !detail::checked_add(prefix, prefix_z, prefix) ||
      !detail::checked_multiply(static_cast<std::size_t>(component_begin),
                                view.component_stride, first_component)) {
    return false;
  }

  std::size_t last_x = 0U;
  std::size_t last_y_index = 0U;
  std::size_t last_z_index = 0U;
  std::size_t last_y = 0U;
  std::size_t last_z = 0U;
  std::size_t last_component = 0U;
  std::size_t maximum_offset = 0U;
  std::size_t end_offset = 0U;
  std::size_t end_bytes = 0U;
  if (!detail::checked_add(static_cast<std::size_t>(view.interior.x),
                           static_cast<std::size_t>(view.ghosts.x), last_x) ||
      !detail::checked_add(static_cast<std::size_t>(view.interior.y),
                           static_cast<std::size_t>(view.ghosts.y),
                           last_y_index) ||
      !detail::checked_add(static_cast<std::size_t>(view.interior.z),
                           static_cast<std::size_t>(view.ghosts.z),
                           last_z_index) ||
      last_x == 0U || last_y_index == 0U || last_z_index == 0U ||
      !detail::checked_multiply(last_y_index - 1U, view.stride_y, last_y) ||
      !detail::checked_multiply(last_z_index - 1U, view.stride_z, last_z) ||
      !detail::checked_multiply(
          static_cast<std::size_t>(component_end - 1U),
          view.component_stride, last_component) ||
      !detail::checked_add(last_x - 1U, last_y, maximum_offset) ||
      !detail::checked_add(maximum_offset, last_z, maximum_offset) ||
      !detail::checked_add(maximum_offset, last_component, maximum_offset) ||
      !detail::checked_add(maximum_offset, 1U, end_offset) ||
      !detail::checked_multiply(end_offset, sizeof(double), end_bytes)) {
    return false;
  }

  const auto base = reinterpret_cast<std::uintptr_t>(view.base);
  std::size_t begin_delta = 0U;
  if (first_component >= prefix) {
    std::size_t begin_bytes = 0U;
    begin_delta = first_component - prefix;
    if (!detail::checked_multiply(begin_delta, sizeof(double), begin_bytes) ||
        begin_bytes > std::numeric_limits<std::uintptr_t>::max() - base) {
      return false;
    }
    begin = base + begin_bytes;
  } else {
    std::size_t begin_bytes = 0U;
    begin_delta = prefix - first_component;
    if (!detail::checked_multiply(begin_delta, sizeof(double), begin_bytes) ||
        begin_bytes > base) {
      return false;
    }
    begin = base - begin_bytes;
  }
  if (end_bytes > std::numeric_limits<std::uintptr_t>::max() - base) {
    return false;
  }
  end = base + end_bytes;
  return begin < end;
}

}  // namespace

namespace detail {

bool valid_cells(Int3 cells) noexcept {
  return cells.x > 0 && cells.y > 0 && cells.z > 0;
}

bool same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool valid_kernel_box(KernelBox box, Int3 cells) noexcept {
  if (!valid_cells(box.cells) || box.begin.x < 0 || box.begin.y < 0 ||
      box.begin.z < 0 || !valid_cells(cells)) {
    return false;
  }
  return static_cast<std::int64_t>(box.begin.x) + box.cells.x <= cells.x &&
         static_cast<std::int64_t>(box.begin.y) + box.cells.y <= cells.y &&
         static_cast<std::int64_t>(box.begin.z) + box.cells.z <= cells.z;
}

std::uint64_t box_cell_count(KernelBox box) noexcept {
  if (!valid_cells(box.cells)) {
    return 0U;
  }
  const auto x = static_cast<std::uint64_t>(box.cells.x);
  const auto y = static_cast<std::uint64_t>(box.cells.y);
  const auto z = static_cast<std::uint64_t>(box.cells.z);
  if (x > std::numeric_limits<std::uint64_t>::max() / y ||
      x * y > std::numeric_limits<std::uint64_t>::max() / z) {
    return 0U;
  }
  return x * y * z;
}

bool valid_cell_view(ConstFieldView view, Int3 cells,
                     std::uint8_t component_begin,
                     std::uint8_t component_count,
                     std::uint8_t required_ghost_width) noexcept {
  const unsigned end = static_cast<unsigned>(component_begin) +
                       static_cast<unsigned>(component_count);
  FieldStorageInterval interval{};
  return field_storage_interval(view, interval) &&
         same_cells(view.interior, cells) &&
         view.components != 0U && component_count != 0U &&
         end <= view.components && view.stride_y != 0U &&
         view.stride_z != 0U && view.component_stride != 0U &&
         view.revision != 0U && view.storage_identity != 0U &&
         view.revision_domain != 0U &&
         view.ghosts.x >= required_ghost_width &&
         view.ghosts.y >= required_ghost_width &&
         view.ghosts.z >= required_ghost_width;
}

bool valid_cell_view(FieldView view, Int3 cells,
                     std::uint8_t component_begin,
                     std::uint8_t component_count) noexcept {
  return valid_cell_view(as_const(view), cells, component_begin,
                         component_count, 0U);
}

bool valid_face(ConstFaceFieldView view, CartesianAxis expected,
                Int3 extents) noexcept {
  FieldStorageInterval interval{};
  return face_storage_interval(view, interval) && view.axis == expected &&
         same_cells(view.extents, extents) && view.storage_identity != 0U &&
         view.revision_domain != 0U;
}

bool valid_flux_view(ConstFaceFluxView flux, Int3 cells,
                     RevisionToken required_revision) noexcept {
  if (!valid_cells(cells) ||
      cells.x == std::numeric_limits<std::int32_t>::max() ||
      cells.y == std::numeric_limits<std::int32_t>::max() ||
      cells.z == std::numeric_limits<std::int32_t>::max()) {
    return false;
  }
  return flux.revision != 0U &&
         (required_revision == 0U || flux.revision == required_revision) &&
         valid_face(flux.x, CartesianAxis::x,
                    {cells.x + 1, cells.y, cells.z}) &&
         valid_face(flux.y, CartesianAxis::y,
                    {cells.x, cells.y + 1, cells.z}) &&
         valid_face(flux.z, CartesianAxis::z,
                    {cells.x, cells.y, cells.z + 1}) &&
         !face_views_overlap(flux.x, flux.y) &&
         !face_views_overlap(flux.x, flux.z) &&
         !face_views_overlap(flux.y, flux.z) &&
         flux.x.storage_identity == flux.y.storage_identity &&
         flux.x.storage_identity == flux.z.storage_identity &&
         flux.x.revision_domain == flux.y.revision_domain &&
         flux.x.revision_domain == flux.z.revision_domain;
}

bool valid_flux_view(FaceFluxView flux, Int3 cells) noexcept {
  return valid_flux_view(as_const(flux), cells, flux.revision);
}

bool component_ranges_overlap(ConstFieldView input,
                              std::uint8_t input_component_begin,
                              std::uint8_t input_component_count,
                              FieldView output,
                              std::uint8_t output_component_begin,
                              std::uint8_t output_component_count) noexcept {
  std::uintptr_t input_begin = 0U;
  std::uintptr_t input_end = 0U;
  std::uintptr_t output_begin = 0U;
  std::uintptr_t output_end = 0U;
  if (!component_byte_interval(input, input_component_begin,
                               input_component_count, input_begin,
                               input_end) ||
      !component_byte_interval(output, output_component_begin,
                               output_component_count, output_begin,
                               output_end)) {
    // Invalid or overflowing metadata is not safe to use as a no-alias proof.
    return true;
  }
  return input_begin < output_end && output_begin < input_end;
}

double cell_width(const CartesianKernelPlan& plan, CartesianAxis axis,
                  std::int32_t local_index) noexcept {
  const std::size_t selected = static_cast<std::size_t>(axis);
  return plan.geometry_kind() == GeometryKind::uniform
             ? metric_width<true>(plan, selected, local_index)
             : metric_width<false>(plan, selected, local_index);
}

double inverse_cell_width(const CartesianKernelPlan& plan, CartesianAxis axis,
                          std::int32_t local_index) noexcept {
  const std::size_t selected = static_cast<std::size_t>(axis);
  return plan.geometry_kind() == GeometryKind::uniform
             ? metric_inverse_width<true>(plan, selected, local_index)
             : metric_inverse_width<false>(plan, selected, local_index);
}

double face_area(const CartesianKernelPlan& plan, CartesianAxis axis,
                 Int3 local_face) noexcept {
  if (axis == CartesianAxis::x) {
    return cell_width(plan, CartesianAxis::y, local_face.y) *
           cell_width(plan, CartesianAxis::z, local_face.z);
  }
  if (axis == CartesianAxis::y) {
    return cell_width(plan, CartesianAxis::x, local_face.x) *
           cell_width(plan, CartesianAxis::z, local_face.z);
  }
  return cell_width(plan, CartesianAxis::x, local_face.x) *
         cell_width(plan, CartesianAxis::y, local_face.y);
}

double cell_volume(const CartesianKernelPlan& plan, Int3 local) noexcept {
  return cell_width(plan, CartesianAxis::x, local.x) *
         cell_width(plan, CartesianAxis::y, local.y) *
         cell_width(plan, CartesianAxis::z, local.z);
}

double centre_coordinate(const CartesianKernelPlan& plan, CartesianAxis axis,
                         std::int32_t local_index) noexcept {
  const std::size_t selected = static_cast<std::size_t>(axis);
  return plan.geometry_kind() == GeometryKind::uniform
             ? metric_centre<true>(plan, selected, local_index)
             : metric_centre<false>(plan, selected, local_index);
}

double face_coordinate(const CartesianKernelPlan& plan, CartesianAxis axis,
                       std::int32_t local_face) noexcept {
  const std::size_t selected = static_cast<std::size_t>(axis);
  return plan.geometry_kind() == GeometryKind::uniform
             ? metric_face<true>(plan, selected, local_face)
             : metric_face<false>(plan, selected, local_face);
}

DerivativeWeights derivative_weights(const CartesianKernelPlan& plan,
                                     CartesianAxis axis,
                                     std::int32_t local_index) noexcept {
  const std::size_t selected = static_cast<std::size_t>(axis);
  return plan.geometry_kind() == GeometryKind::uniform
             ? metric_derivative_weights<true>(plan, selected, local_index)
             : metric_derivative_weights<false>(plan, selected, local_index);
}

double interpolate_face(const CartesianKernelPlan& plan, CartesianAxis axis,
                        std::int32_t local_face, double left,
                        double right) noexcept {
  const std::size_t selected = static_cast<std::size_t>(axis);
  return plan.geometry_kind() == GeometryKind::uniform
             ? metric_interpolate_face<true>(plan, selected, local_face, left,
                                             right)
             : metric_interpolate_face<false>(plan, selected, local_face,
                                              left, right);
}

ConstFaceFieldView select(ConstFaceFluxView flux,
                          CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

FaceFieldView select(FaceFluxView flux, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

bool checked_counter_add(std::uint64_t& destination,
                         std::uint64_t increment) noexcept {
  if (increment > std::numeric_limits<std::uint64_t>::max() - destination) {
    return false;
  }
  destination += increment;
  return true;
}

Status add_kernel_counters(KernelCounters* counters, KernelBox box,
                           std::uint64_t faces_per_cell,
                           std::uint64_t read_doubles_per_cell,
                           std::uint64_t written_doubles_per_cell) noexcept {
  const std::uint64_t cells = box_cell_count(box);
  if (cells == 0U ||
      faces_per_cell > std::numeric_limits<std::uint64_t>::max() / cells ||
      read_doubles_per_cell >
          std::numeric_limits<std::uint64_t>::max() / cells ||
      written_doubles_per_cell >
          std::numeric_limits<std::uint64_t>::max() / cells) {
    return {StatusCode::invalid_plan, kKernelInvocation};
  }
  return add_kernel_counter_totals(counters, cells,
                                   cells * faces_per_cell,
                                   cells * read_doubles_per_cell,
                                   cells * written_doubles_per_cell);
}

Status add_kernel_counter_totals(KernelCounters* counters,
                                 std::uint64_t cells,
                                 std::uint64_t faces,
                                 std::uint64_t reads,
                                 std::uint64_t writes) noexcept {
  if (counters == nullptr) {
    return {};
  }
  if (cells == 0U) {
    return {StatusCode::invalid_plan, kKernelInvocation};
  }
  if (reads > std::numeric_limits<std::uint64_t>::max() / sizeof(double) ||
      writes > std::numeric_limits<std::uint64_t>::max() / sizeof(double)) {
    return {StatusCode::invalid_plan, kKernelInvocation};
  }
  KernelCounters candidate = *counters;
  if (!checked_counter_add(candidate.invocations, 1U) ||
      !checked_counter_add(candidate.cells, cells) ||
      !checked_counter_add(candidate.faces, faces) ||
      !checked_counter_add(candidate.logical_bytes_read,
                           reads * sizeof(double)) ||
      !checked_counter_add(candidate.logical_bytes_written,
                           writes * sizeof(double))) {
    return {StatusCode::invalid_plan, kKernelInvocation};
  }
  *counters = candidate;
  return {};
}

}  // namespace detail

namespace {

template <bool Uniform>
Status gradient_kernel(const CartesianKernelPlan& plan,
                       const KernelInvocation& invocation,
                       ConstFieldView input, FieldView output) noexcept {
  KernelCounters prepared{};
  if (invocation.counters != nullptr) {
    prepared = *invocation.counters;
    const Status counted = detail::add_kernel_counters(
        &prepared, invocation.box, 0U,
        static_cast<std::uint64_t>(7U) * invocation.component_count,
        static_cast<std::uint64_t>(3U) * invocation.component_count);
    if (!counted) {
      return counted;
    }
  }
  const Int3 end{invocation.box.begin.x + invocation.box.cells.x,
                 invocation.box.begin.y + invocation.box.cells.y,
                 invocation.box.begin.z + invocation.box.cells.z};
  for (std::uint8_t component_id = 0U;
       component_id < invocation.component_count; ++component_id) {
    const std::uint8_t read_component = static_cast<std::uint8_t>(
        invocation.read_component_begin + component_id);
    const std::uint8_t write_component = static_cast<std::uint8_t>(
        invocation.write_component_begin + 3U * component_id);
    for (std::int32_t z = invocation.box.begin.z; z < end.z; ++z) {
      const detail::DerivativeWeights wz =
          detail::metric_derivative_weights<Uniform>(plan, 2U, z);
      for (std::int32_t y = invocation.box.begin.y; y < end.y; ++y) {
        const detail::DerivativeWeights wy =
            detail::metric_derivative_weights<Uniform>(plan, 1U, y);
        for (std::int32_t x = invocation.box.begin.x; x < end.x; ++x) {
          const detail::DerivativeWeights wx =
              detail::metric_derivative_weights<Uniform>(plan, 0U, x);
          const Int3 centre{x, y, z};
          const double q_centre = input.unchecked(centre, read_component);
          const double gx =
              wx.minus * input.unchecked({x - 1, y, z}, read_component) +
              wx.centre * q_centre +
              wx.plus * input.unchecked({x + 1, y, z}, read_component);
          const double gy =
              wy.minus * input.unchecked({x, y - 1, z}, read_component) +
              wy.centre * q_centre +
              wy.plus * input.unchecked({x, y + 1, z}, read_component);
          const double gz =
              wz.minus * input.unchecked({x, y, z - 1}, read_component) +
              wz.centre * q_centre +
              wz.plus * input.unchecked({x, y, z + 1}, read_component);
          if (!std::isfinite(gx) || !std::isfinite(gy) ||
              !std::isfinite(gz)) {
            return {StatusCode::numerical_failure, kKernelNumerical};
          }
          output.unchecked(centre, write_component) = gx;
          output.unchecked(
              centre, static_cast<std::uint8_t>(write_component + 1U)) = gy;
          output.unchecked(
              centre, static_cast<std::uint8_t>(write_component + 2U)) = gz;
        }
      }
    }
  }
  if (invocation.counters != nullptr) {
    *invocation.counters = prepared;
  }
  return {};
}

}  // namespace

CartesianKernelPlan::CartesianKernelPlan(CartesianKernelPlan&& other) noexcept {
  move_from(std::move(other));
}

CartesianKernelPlan& CartesianKernelPlan::operator=(
    CartesianKernelPlan&& other) noexcept {
  if (this != &other) {
    reset();
    move_from(std::move(other));
  }
  return *this;
}

void CartesianKernelPlan::reset() noexcept {
  patch_begin_ = {};
  cells_ = {};
  process_grid_ = {};
  process_coord_ = {};
  boundary_identity_ = 0U;
  scheme_identity_ = 0U;
  fingerprint_ = 0U;
  limiter_ = 1.0;
  reach_ = 0U;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    metric_faces_[axis].clear();
    metric_centres_[axis].clear();
    metric_widths_[axis].clear();
    metric_inverse_widths_[axis].clear();
    metrics_[axis] = {};
  }
}

void CartesianKernelPlan::rebind_metrics() noexcept {
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    metrics_[axis].faces = metric_faces_[axis].data();
    metrics_[axis].centres = metric_centres_[axis].data();
    metrics_[axis].widths = metric_widths_[axis].data();
    metrics_[axis].inverse_widths = metric_inverse_widths_[axis].data();
  }
}

void CartesianKernelPlan::move_from(CartesianKernelPlan&& other) noexcept {
  patch_begin_ = other.patch_begin_;
  cells_ = other.cells_;
  process_grid_ = other.process_grid_;
  process_coord_ = other.process_coord_;
  boundary_identity_ = other.boundary_identity_;
  scheme_identity_ = other.scheme_identity_;
  fingerprint_ = other.fingerprint_;
  geometry_kind_ = other.geometry_kind_;
  limiter_ = other.limiter_;
  reach_ = other.reach_;
  for (std::size_t axis = 0U; axis < 3U; ++axis) {
    metric_faces_[axis] = std::move(other.metric_faces_[axis]);
    metric_centres_[axis] = std::move(other.metric_centres_[axis]);
    metric_widths_[axis] = std::move(other.metric_widths_[axis]);
    metric_inverse_widths_[axis] =
        std::move(other.metric_inverse_widths_[axis]);
    metrics_[axis] = other.metrics_[axis];
  }
  rebind_metrics();
  other.reset();
}

Status CartesianKernelPlan::compile(const SchemePlan& schemes,
                                    const CartesianGeometryPlan& geometry,
                                    const MeshPatch& patch,
                                    const BoundaryPlan& boundary,
                                    CartesianKernelPlan& out) noexcept {
  const Int3 global = geometry.global_cells();
  const std::uint8_t reach = schemes.required_ghost_width();
  const std::int64_t end_x = static_cast<std::int64_t>(patch.begin.x) +
                             static_cast<std::int64_t>(patch.cells.x);
  const std::int64_t end_y = static_cast<std::int64_t>(patch.begin.y) +
                             static_cast<std::int64_t>(patch.cells.y);
  const std::int64_t end_z = static_cast<std::int64_t>(patch.begin.z) +
                             static_cast<std::int64_t>(patch.cells.z);
  if (geometry.fingerprint() == 0U || geometry.topology_revision() == 0U ||
      boundary.semantic_fingerprint() == 0U ||
      boundary.local_layout_fingerprint() == 0U ||
      schemes.fingerprint() == 0U || !detail::valid_cells(patch.cells) ||
      patch.begin.x < 0 || patch.begin.y < 0 || patch.begin.z < 0 ||
      end_x > global.x || end_y > global.y || end_z > global.z ||
      !detail::same_cells(boundary.local_cells(), patch.cells) || reach == 0U ||
      boundary.required_ghost_width() < reach) {
    return {StatusCode::invalid_plan, kKernelPlan};
  }
  std::uint64_t fingerprint = kFnvOffset;
  fingerprint = hash_mix(fingerprint, geometry.fingerprint());
  fingerprint = hash_mix(fingerprint, geometry.topology_revision());
  fingerprint = hash_mix(fingerprint, schemes.fingerprint());
  fingerprint = hash_mix(fingerprint, boundary.semantic_fingerprint());
  fingerprint = hash_mix(fingerprint, boundary.local_layout_fingerprint());
  fingerprint = hash_mix(fingerprint, static_cast<std::uint64_t>(patch.begin.x));
  fingerprint = hash_mix(fingerprint, static_cast<std::uint64_t>(patch.begin.y));
  fingerprint = hash_mix(fingerprint, static_cast<std::uint64_t>(patch.begin.z));
  fingerprint = hash_mix(fingerprint, static_cast<std::uint64_t>(patch.cells.x));
  fingerprint = hash_mix(fingerprint, static_cast<std::uint64_t>(patch.cells.y));
  fingerprint = hash_mix(fingerprint, static_cast<std::uint64_t>(patch.cells.z));
  fingerprint = hash_mix(fingerprint, reach);
  try {
    CartesianKernelPlan candidate;
    candidate.patch_begin_ = patch.begin;
    candidate.cells_ = patch.cells;
    candidate.process_grid_ = patch.process_grid;
    candidate.process_coord_ = patch.process_coord;
    candidate.boundary_identity_ = boundary.local_layout_fingerprint();
    candidate.scheme_identity_ = schemes.fingerprint();
    candidate.fingerprint_ = fingerprint == 0U ? 1U : fingerprint;
    candidate.geometry_kind_ = geometry.kind();
    candidate.limiter_ = schemes.limiter();
    candidate.reach_ = reach;
    const AxisMetrics* const source[3]{&geometry.x(), &geometry.y(),
                                       &geometry.z()};
    const std::int32_t begins[3]{patch.begin.x, patch.begin.y,
                                 patch.begin.z};
    const std::int32_t local_cells[3]{patch.cells.x, patch.cells.y,
                                      patch.cells.z};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const AxisMetrics& metric = *source[axis];
      const Span<const double> faces = metric.faces();
      const Span<const double> centres = metric.centres();
      const Span<const double> widths = metric.widths();
      const Span<const double> inverse_widths = metric.inverse_widths();
      if (faces.data == nullptr || centres.data == nullptr ||
          widths.data == nullptr || inverse_widths.data == nullptr ||
          centres.size == 0U || faces.size != centres.size + 1U ||
          widths.size != centres.size ||
          inverse_widths.size != centres.size ||
          static_cast<std::size_t>(begins[axis]) >= centres.size) {
        return {StatusCode::invalid_plan, kKernelPlan};
      }
      const std::int64_t requested_begin =
          static_cast<std::int64_t>(begins[axis]) - reach;
      const std::int64_t requested_end =
          static_cast<std::int64_t>(begins[axis]) + local_cells[axis] + reach;
      const std::size_t slice_begin = static_cast<std::size_t>(
          std::max<std::int64_t>(0, requested_begin));
      const std::size_t slice_end = static_cast<std::size_t>(
          std::min<std::int64_t>(static_cast<std::int64_t>(centres.size),
                                 requested_end));
      if (slice_begin >= slice_end ||
          static_cast<std::size_t>(begins[axis]) < slice_begin) {
        return {StatusCode::invalid_plan, kKernelPlan};
      }
      candidate.metric_faces_[axis].assign(faces.data + slice_begin,
                                            faces.data + slice_end + 1U);
      candidate.metric_centres_[axis].assign(centres.data + slice_begin,
                                              centres.data + slice_end);
      candidate.metric_widths_[axis].assign(widths.data + slice_begin,
                                             widths.data + slice_end);
      candidate.metric_inverse_widths_[axis].assign(
          inverse_widths.data + slice_begin,
          inverse_widths.data + slice_end);
      candidate.metrics_[axis] = {
          nullptr,
          nullptr,
          nullptr,
          nullptr,
          slice_end - slice_begin,
          static_cast<std::int32_t>(
              static_cast<std::size_t>(begins[axis]) - slice_begin),
          faces.data[static_cast<std::size_t>(begins[axis])],
          centres.data[static_cast<std::size_t>(begins[axis])],
          metric.uniform_width(),
          metric.uniform_inverse_width()};
    }
    candidate.rebind_metrics();
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kKernelPlan};
  }
}

Status cartesian_gradient(const CartesianKernelPlan& plan,
                          const KernelInvocation& invocation) noexcept {
  if (plan.fingerprint_ == 0U ||
      invocation.reads.data == nullptr || invocation.reads.size != 1U ||
      invocation.writes.data == nullptr || invocation.writes.size != 1U ||
      !detail::valid_kernel_box(invocation.box, plan.cells_) ||
      invocation.component_count == 0U ||
      static_cast<unsigned>(invocation.write_component_begin) +
              3U * invocation.component_count >
          invocation.writes.data[0U].components ||
      !detail::valid_cell_view(invocation.reads.data[0U], plan.cells_,
                               invocation.read_component_begin,
                               invocation.component_count, 1U) ||
      !detail::valid_cell_view(invocation.writes.data[0U], plan.cells_,
                               invocation.write_component_begin,
                               static_cast<std::uint8_t>(
                                   3U * invocation.component_count)) ||
      detail::component_ranges_overlap(
          invocation.reads.data[0U], invocation.read_component_begin,
          invocation.component_count, invocation.writes.data[0U],
          invocation.write_component_begin,
          static_cast<std::uint8_t>(3U * invocation.component_count))) {
    return {StatusCode::invalid_plan, kKernelInvocation};
  }
  const ConstFieldView input = invocation.reads.data[0U];
  const FieldView output = invocation.writes.data[0U];
  return plan.geometry_kind_ == GeometryKind::uniform
             ? gradient_kernel<true>(plan, invocation, input, output)
             : gradient_kernel<false>(plan, invocation, input, output);
}

}  // namespace hundun::v04
