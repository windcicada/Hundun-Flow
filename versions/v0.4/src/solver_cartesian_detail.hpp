// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_execution.hpp"
#include "hundun/v04_mesh.hpp"

#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

struct DerivativeWeights {
  double minus{};
  double centre{};
  double plus{};
};

bool valid_cells(Int3 cells) noexcept;
bool same_cells(Int3 left, Int3 right) noexcept;
bool valid_kernel_box(KernelBox box, Int3 cells) noexcept;
std::uint64_t box_cell_count(KernelBox box) noexcept;
bool valid_cell_view(ConstFieldView view, Int3 cells,
                     std::uint8_t component_begin,
                     std::uint8_t component_count,
                     std::uint8_t required_ghost_width) noexcept;
bool valid_cell_view(FieldView view, Int3 cells,
                     std::uint8_t component_begin,
                     std::uint8_t component_count) noexcept;
bool valid_flux_view(ConstFaceFluxView flux, Int3 cells,
                     RevisionToken required_revision) noexcept;
bool valid_flux_view(FaceFluxView flux, Int3 cells) noexcept;

bool component_ranges_overlap(ConstFieldView input,
                              std::uint8_t input_component_begin,
                              std::uint8_t input_component_count,
                              FieldView output,
                              std::uint8_t output_component_begin,
                              std::uint8_t output_component_count) noexcept;

double cell_width(const CartesianKernelPlan& plan, CartesianAxis axis,
                  std::int32_t local_index) noexcept;
double inverse_cell_width(const CartesianKernelPlan& plan, CartesianAxis axis,
                          std::int32_t local_index) noexcept;
double face_area(const CartesianKernelPlan& plan, CartesianAxis axis,
                 Int3 local_face) noexcept;
double cell_volume(const CartesianKernelPlan& plan, Int3 local) noexcept;
double centre_coordinate(const CartesianKernelPlan& plan, CartesianAxis axis,
                         std::int32_t local_index) noexcept;
double face_coordinate(const CartesianKernelPlan& plan, CartesianAxis axis,
                       std::int32_t local_face) noexcept;
DerivativeWeights derivative_weights(const CartesianKernelPlan& plan,
                                     CartesianAxis axis,
                                     std::int32_t local_index) noexcept;
double interpolate_face(const CartesianKernelPlan& plan, CartesianAxis axis,
                        std::int32_t local_face, double left,
                        double right) noexcept;

ConstFaceFieldView select(ConstFaceFluxView flux,
                          CartesianAxis axis) noexcept;
FaceFieldView select(FaceFluxView flux, CartesianAxis axis) noexcept;

bool checked_counter_add(std::uint64_t& destination,
                         std::uint64_t increment) noexcept;
Status add_kernel_counters(KernelCounters* counters, KernelBox box,
                           std::uint64_t faces_per_cell,
                           std::uint64_t read_doubles_per_cell,
                           std::uint64_t written_doubles_per_cell) noexcept;
Status add_kernel_counter_totals(KernelCounters* counters,
                                 std::uint64_t cells,
                                 std::uint64_t faces,
                                 std::uint64_t read_doubles,
                                 std::uint64_t written_doubles) noexcept;

template <bool Uniform>
inline double metric_width(const CartesianKernelPlan& plan,
                           std::size_t axis,
                           std::int32_t local_index) noexcept {
  const CartesianMetricPacket& metric = plan.metric(axis);
  if constexpr (Uniform) {
    return metric.uniform_width;
  }
  return metric.widths[static_cast<std::size_t>(metric.global_begin +
                                                local_index)];
}

template <bool Uniform>
inline double metric_inverse_width(const CartesianKernelPlan& plan,
                                   std::size_t axis,
                                   std::int32_t local_index) noexcept {
  const CartesianMetricPacket& metric = plan.metric(axis);
  if constexpr (Uniform) {
    return metric.uniform_inverse_width;
  }
  return metric.inverse_widths[static_cast<std::size_t>(metric.global_begin +
                                                        local_index)];
}

template <bool Uniform>
inline double metric_centre(const CartesianKernelPlan& plan,
                            std::size_t axis,
                            std::int32_t local_index) noexcept {
  const CartesianMetricPacket& metric = plan.metric(axis);
  if constexpr (Uniform) {
    return metric.local_centre_origin +
           static_cast<double>(local_index) * metric.uniform_width;
  }
  const std::int32_t global = metric.global_begin + local_index;
  if (global < 0) {
    return metric.centres[0U] + static_cast<double>(global) * metric.widths[0U];
  }
  if (static_cast<std::size_t>(global) >= metric.cells) {
    const auto last = static_cast<std::int32_t>(metric.cells - 1U);
    return metric.centres[metric.cells - 1U] +
           static_cast<double>(global - last) * metric.widths[metric.cells - 1U];
  }
  return metric.centres[static_cast<std::size_t>(global)];
}

template <bool Uniform>
inline double metric_face(const CartesianKernelPlan& plan, std::size_t axis,
                          std::int32_t local_face) noexcept {
  const CartesianMetricPacket& metric = plan.metric(axis);
  if constexpr (Uniform) {
    return metric.local_face_origin +
           static_cast<double>(local_face) * metric.uniform_width;
  }
  const std::int32_t global = metric.global_begin + local_face;
  if (global < 0) {
    return metric.faces[0U] + static_cast<double>(global) * metric.widths[0U];
  }
  if (static_cast<std::size_t>(global) > metric.cells) {
    const auto last = static_cast<std::int32_t>(metric.cells);
    return metric.faces[metric.cells] +
           static_cast<double>(global - last) * metric.widths[metric.cells - 1U];
  }
  return metric.faces[static_cast<std::size_t>(global)];
}

template <bool Uniform>
inline DerivativeWeights metric_derivative_weights(
    const CartesianKernelPlan& plan, std::size_t axis,
    std::int32_t local_index) noexcept {
  if constexpr (Uniform) {
    const double half_inverse =
        0.5 * plan.metric(axis).uniform_inverse_width;
    return {-half_inverse, 0.0, half_inverse};
  }
  const double centre = metric_centre<false>(plan, axis, local_index);
  const double minus = metric_centre<false>(plan, axis, local_index - 1);
  const double plus = metric_centre<false>(plan, axis, local_index + 1);
  const double left = centre - minus;
  const double right = plus - centre;
  const double total = left + right;
  return {-right / (left * total), (right - left) / (left * right),
          left / (right * total)};
}

template <bool Uniform>
inline double metric_interpolate_face(const CartesianKernelPlan& plan,
                                      std::size_t axis,
                                      std::int32_t local_face, double left,
                                      double right) noexcept {
  if constexpr (Uniform) {
    return 0.5 * (left + right);
  }
  const double face = metric_face<false>(plan, axis, local_face);
  const double left_centre =
      metric_centre<false>(plan, axis, local_face - 1);
  const double right_centre = metric_centre<false>(plan, axis, local_face);
  const double left_distance = face - left_centre;
  const double right_distance = right_centre - face;
  return (right_distance * left + left_distance * right) /
         (left_distance + right_distance);
}

template <bool Uniform>
inline double metric_inverse_volume(const CartesianKernelPlan& plan,
                                    Int3 local) noexcept {
  return metric_inverse_width<Uniform>(plan, 0U, local.x) *
         metric_inverse_width<Uniform>(plan, 1U, local.y) *
         metric_inverse_width<Uniform>(plan, 2U, local.z);
}

template <bool Uniform>
inline double metric_face_area(const CartesianKernelPlan& plan,
                               std::size_t axis, Int3 face) noexcept {
  if (axis == 0U) {
    return metric_width<Uniform>(plan, 1U, face.y) *
           metric_width<Uniform>(plan, 2U, face.z);
  }
  if (axis == 1U) {
    return metric_width<Uniform>(plan, 0U, face.x) *
           metric_width<Uniform>(plan, 2U, face.z);
  }
  return metric_width<Uniform>(plan, 0U, face.x) *
         metric_width<Uniform>(plan, 1U, face.y);
}

}  // namespace hundun::v04::detail
