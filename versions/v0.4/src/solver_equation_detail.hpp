// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_flow.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04::detail {

inline bool select_contribution_stage(
    Span<const CompiledContribution> descriptors, StageId stage,
    Span<const CompiledContribution>& selected) noexcept {
  selected = {};
  if (stage == 0U ||
      (descriptors.size != 0U && descriptors.data == nullptr)) {
    return false;
  }
  std::size_t begin = descriptors.size;
  std::size_t count = 0U;
  bool left_selected_stage = false;
  for (std::size_t index = 0U; index < descriptors.size; ++index) {
    if (descriptors.data[index].stage == stage) {
      if (left_selected_stage) {
        return false;
      }
      if (begin == descriptors.size) {
        begin = index;
      }
      ++count;
    } else if (begin != descriptors.size) {
      left_selected_stage = true;
    }
  }
  if (count != 0U) {
    selected = {descriptors.data + begin, count};
  }
  return true;
}

inline bool equation_same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

inline bool full_equation_box(KernelBox box, Int3 cells) noexcept {
  return box.begin.x == 0 && box.begin.y == 0 && box.begin.z == 0 &&
         equation_same_cells(box.cells, cells);
}

inline bool valid_bdf_coefficients(BdfCoefficients bdf) noexcept {
  if (!std::isfinite(bdf.a0) || !std::isfinite(bdf.a1) ||
      !std::isfinite(bdf.a2) || bdf.a0 <= 0.0 || bdf.order < 1U ||
      bdf.order > 2U) {
    return false;
  }
  const double scale = std::max(
      {1.0, std::abs(bdf.a0), std::abs(bdf.a1), std::abs(bdf.a2)});
  if (std::abs(bdf.a0 + bdf.a1 + bdf.a2) >
      64.0 * std::numeric_limits<double>::epsilon() * scale) {
    return false;
  }
  return bdf.order == 1U ? bdf.a1 < 0.0 && bdf.a2 == 0.0
                         : bdf.a1 < 0.0 && bdf.a2 > 0.0;
}

inline bool finite_field_box(ConstFieldView field, KernelBox box,
                             std::uint8_t component_begin,
                             std::uint8_t component_count) noexcept {
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        for (std::uint8_t component = 0U; component < component_count;
             ++component) {
          if (!std::isfinite(field.unchecked(
                  {x, y, z},
                  static_cast<std::uint8_t>(component_begin + component)))) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

inline bool finite_face_neighbour_slabs(ConstFieldView field, KernelBox box,
                                        std::uint8_t component_begin,
                                        std::uint8_t component_count,
                                        std::uint8_t reach = 1U) noexcept {
  if (reach == 0U) {
    return finite_field_box(field, box, component_begin, component_count);
  }
  if (!finite_field_box(field, box, component_begin, component_count)) {
    return false;
  }
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  const auto finite_cell = [&](Int3 cell) noexcept {
    for (std::uint8_t component = 0U; component < component_count;
         ++component) {
      if (!std::isfinite(field.unchecked(
              cell, static_cast<std::uint8_t>(component_begin + component)))) {
        return false;
      }
    }
    return true;
  };
  for (std::uint8_t layer = 1U; layer <= reach; ++layer) {
    const auto offset = static_cast<std::int32_t>(layer);
    for (std::int32_t z = box.begin.z; z < end.z; ++z) {
      for (std::int32_t y = box.begin.y; y < end.y; ++y) {
        if (!finite_cell({box.begin.x - offset, y, z}) ||
            !finite_cell({end.x - 1 + offset, y, z})) {
          return false;
        }
      }
    }
    for (std::int32_t z = box.begin.z; z < end.z; ++z) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        if (!finite_cell({x, box.begin.y - offset, z}) ||
            !finite_cell({x, end.y - 1 + offset, z})) {
          return false;
        }
      }
    }
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        if (!finite_cell({x, y, box.begin.z - offset}) ||
            !finite_cell({x, y, end.z - 1 + offset})) {
          return false;
        }
      }
    }
  }
  return true;
}

inline bool finite_face_flux(ConstFaceFluxView flux, KernelBox box) noexcept {
  const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                 box.begin.z + box.cells.z};
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x <= end.x; ++x) {
        if (!std::isfinite(flux.x.unchecked({x, y, z}))) {
          return false;
        }
      }
    }
  }
  for (std::int32_t z = box.begin.z; z < end.z; ++z) {
    for (std::int32_t y = box.begin.y; y <= end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        if (!std::isfinite(flux.y.unchecked({x, y, z}))) {
          return false;
        }
      }
    }
  }
  for (std::int32_t z = box.begin.z; z <= end.z; ++z) {
    for (std::int32_t y = box.begin.y; y < end.y; ++y) {
      for (std::int32_t x = box.begin.x; x < end.x; ++x) {
        if (!std::isfinite(flux.z.unchecked({x, y, z}))) {
          return false;
        }
      }
    }
  }
  return true;
}

inline bool valid_equation_face_view(FaceFieldView view, CartesianAxis axis,
                                     Int3 cells) noexcept {
  Int3 expected = cells;
  if (axis == CartesianAxis::x) {
    ++expected.x;
  } else if (axis == CartesianAxis::y) {
    ++expected.y;
  } else {
    ++expected.z;
  }
  if (view.base == nullptr || view.axis != axis ||
      !equation_same_cells(view.extents, expected) || view.stride_y == 0U ||
      view.stride_z == 0U || view.storage_identity == 0U ||
      view.revision_domain == 0U) {
    return false;
  }
  const std::size_t nx = static_cast<std::size_t>(expected.x);
  const std::size_t ny = static_cast<std::size_t>(expected.y);
  return view.stride_y >= nx &&
         (ny == 0U || view.stride_y <=
                            std::numeric_limits<std::size_t>::max() / ny) &&
         view.stride_z >= view.stride_y * ny;
}

inline bool equation_faces_empty(EquationSystemView system) noexcept {
  return system.x_coefficient.base == nullptr &&
         system.y_coefficient.base == nullptr &&
         system.z_coefficient.base == nullptr;
}

inline bool valid_equation_faces(EquationSystemView system,
                                 Int3 cells) noexcept {
  return valid_equation_face_view(system.x_coefficient, CartesianAxis::x,
                                  cells) &&
         valid_equation_face_view(system.y_coefficient, CartesianAxis::y,
                                  cells) &&
         valid_equation_face_view(system.z_coefficient, CartesianAxis::z,
                                  cells);
}

inline bool equation_face_storage_interval(FaceFieldView view,
                                           std::uintptr_t& begin,
                                           std::uintptr_t& end) noexcept {
  if (view.base == nullptr || view.extents.x <= 0 || view.extents.y <= 0 ||
      view.extents.z <= 0 || view.stride_y == 0U || view.stride_z == 0U) {
    return false;
  }
  const std::size_t last =
      static_cast<std::size_t>(view.extents.x - 1) +
      static_cast<std::size_t>(view.extents.y - 1) * view.stride_y +
      static_cast<std::size_t>(view.extents.z - 1) * view.stride_z;
  if (last == std::numeric_limits<std::size_t>::max() ||
      last + 1U > std::numeric_limits<std::uintptr_t>::max() /
                       sizeof(double)) {
    return false;
  }
  begin = reinterpret_cast<std::uintptr_t>(view.base);
  const auto bytes = static_cast<std::uintptr_t>(last + 1U) * sizeof(double);
  if (bytes > std::numeric_limits<std::uintptr_t>::max() - begin) {
    return false;
  }
  end = begin + bytes;
  return begin < end;
}

inline bool equation_face_views_overlap(FaceFieldView left,
                                        FaceFieldView right) noexcept {
  std::uintptr_t left_begin = 0U;
  std::uintptr_t left_end = 0U;
  std::uintptr_t right_begin = 0U;
  std::uintptr_t right_end = 0U;
  return !equation_face_storage_interval(left, left_begin, left_end) ||
         !equation_face_storage_interval(right, right_begin, right_end) ||
         (left_begin < right_end && right_begin < left_end);
}

inline bool equation_face_views_disjoint(EquationSystemView system) noexcept {
  return !equation_face_views_overlap(system.x_coefficient,
                                      system.y_coefficient) &&
         !equation_face_views_overlap(system.x_coefficient,
                                      system.z_coefficient) &&
         !equation_face_views_overlap(system.y_coefficient,
                                      system.z_coefficient);
}

inline bool cell_output_aliases(ConstFieldView input,
                                EquationSystemView system,
                                bool linear) noexcept {
  return field_views_overlap(input, as_const(system.residual)) ||
         (linear &&
          (field_views_overlap(input, as_const(system.diagonal)) ||
           field_views_overlap(input, as_const(system.rhs))));
}

inline bool cell_input_aliases_faces(ConstFieldView input,
                                     EquationSystemView system,
                                     bool linear) noexcept {
  return linear &&
         (cell_face_views_overlap(input, system.x_coefficient) ||
          cell_face_views_overlap(input, system.y_coefficient) ||
          cell_face_views_overlap(input, system.z_coefficient));
}

inline bool output_aliases_input(ConstFieldView input,
                                 EquationSystemView system,
                                 bool linear) noexcept {
  return cell_output_aliases(input, system, linear) ||
         cell_input_aliases_faces(input, system, linear);
}

inline bool output_aliases_flux(EquationSystemView system, bool linear,
                                ConstFaceFluxView flux) noexcept {
  const ConstFaceFieldView faces[]{flux.x, flux.y, flux.z};
  for (ConstFaceFieldView face : faces) {
    if (cell_face_views_overlap(system.residual, face) ||
        (linear &&
         (cell_face_views_overlap(system.diagonal, face) ||
          cell_face_views_overlap(system.rhs, face) ||
          face_views_overlap(system.x_coefficient, face) ||
          face_views_overlap(system.y_coefficient, face) ||
          face_views_overlap(system.z_coefficient, face)))) {
      return true;
    }
  }
  return false;
}

inline double positive_transmissibility(const CartesianKernelPlan& kernels,
                                        ConstFieldView coefficient,
                                        CartesianAxis axis,
                                        Int3 face) noexcept {
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
  const double left_value = coefficient.unchecked(left, 0U);
  const double right_value = coefficient.unchecked(face, 0U);
  const double face_location = face_coordinate(kernels, axis, normal);
  const double left_distance =
      face_location - centre_coordinate(kernels, axis, normal - 1);
  const double right_distance =
      centre_coordinate(kernels, axis, normal) - face_location;
  if (!std::isfinite(left_value) || !std::isfinite(right_value) ||
      left_value <= 0.0 || right_value <= 0.0 ||
      !std::isfinite(left_distance) || !std::isfinite(right_distance) ||
      left_distance <= 0.0 || right_distance <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return face_area(kernels, axis, face) /
         (left_distance / left_value + right_distance / right_value);
}

inline double diffusion_diagonal(const CartesianKernelPlan& kernels,
                                 ConstFieldView coefficient,
                                 Int3 cell) noexcept {
  return positive_transmissibility(kernels, coefficient, CartesianAxis::x,
                                   cell) +
         positive_transmissibility(kernels, coefficient, CartesianAxis::x,
                                   {cell.x + 1, cell.y, cell.z}) +
         positive_transmissibility(kernels, coefficient, CartesianAxis::y,
                                   cell) +
         positive_transmissibility(kernels, coefficient, CartesianAxis::y,
                                   {cell.x, cell.y + 1, cell.z}) +
         positive_transmissibility(kernels, coefficient, CartesianAxis::z,
                                   cell) +
         positive_transmissibility(kernels, coefficient, CartesianAxis::z,
                                   {cell.x, cell.y, cell.z + 1});
}

inline double negative_diffusion_operator(
    const CartesianKernelPlan& kernels, ConstFieldView coefficient,
    ConstFieldView field, Int3 cell, std::uint8_t component) noexcept {
  const double centre = field.unchecked(cell, component);
  const double txm = positive_transmissibility(
      kernels, coefficient, CartesianAxis::x, cell);
  const double txp = positive_transmissibility(
      kernels, coefficient, CartesianAxis::x,
      {cell.x + 1, cell.y, cell.z});
  const double tym = positive_transmissibility(
      kernels, coefficient, CartesianAxis::y, cell);
  const double typ = positive_transmissibility(
      kernels, coefficient, CartesianAxis::y,
      {cell.x, cell.y + 1, cell.z});
  const double tzm = positive_transmissibility(
      kernels, coefficient, CartesianAxis::z, cell);
  const double tzp = positive_transmissibility(
      kernels, coefficient, CartesianAxis::z,
      {cell.x, cell.y, cell.z + 1});
  return txm * (centre - field.unchecked({cell.x - 1, cell.y, cell.z},
                                        component)) +
         txp * (centre - field.unchecked({cell.x + 1, cell.y, cell.z},
                                        component)) +
         tym * (centre - field.unchecked({cell.x, cell.y - 1, cell.z},
                                        component)) +
         typ * (centre - field.unchecked({cell.x, cell.y + 1, cell.z},
                                        component)) +
         tzm * (centre - field.unchecked({cell.x, cell.y, cell.z - 1},
                                        component)) +
         tzp * (centre - field.unchecked({cell.x, cell.y, cell.z + 1},
                                        component));
}

template <CartesianAxis Axis>
Status fill_equation_face_coefficients(const CartesianKernelPlan& kernels,
                                       ConstFieldView coefficient,
                                       KernelBox box, FaceFieldView output,
                                       std::uint32_t detail_code) noexcept {
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
        const double value =
            positive_transmissibility(kernels, coefficient, Axis, face);
        if (!std::isfinite(value) || value <= 0.0) {
          return {StatusCode::numerical_failure, detail_code};
        }
        output.unchecked(face) = value;
      }
    }
  }
  return {};
}

inline Status fill_equation_face_coefficients(
    const CartesianKernelPlan& kernels, ConstFieldView coefficient,
    KernelBox box, EquationSystemView system,
    std::uint32_t detail_code) noexcept {
  Status status = fill_equation_face_coefficients<CartesianAxis::x>(
      kernels, coefficient, box, system.x_coefficient, detail_code);
  if (status) {
    status = fill_equation_face_coefficients<CartesianAxis::y>(
        kernels, coefficient, box, system.y_coefficient, detail_code);
  }
  if (status) {
    status = fill_equation_face_coefficients<CartesianAxis::z>(
        kernels, coefficient, box, system.z_coefficient, detail_code);
  }
  return status;
}

}  // namespace hundun::v04::detail
