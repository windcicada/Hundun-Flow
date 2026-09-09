// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_flow.hpp"
#include "solver_equation_detail.hpp"

namespace hundun::v04::detail {

inline Int3 viscous_offset(Int3 cell, std::size_t axis, int delta) noexcept {
  if (axis == 0U) cell.x += delta;
  else if (axis == 1U) cell.y += delta;
  else cell.z += delta;
  return cell;
}

// The same transpose/deviatoric traction for momentum, kinetic work and
// compatible enthalpy heating. Preserve the momentum kernel's operation order.
inline double viscous_cross_traction(const CartesianKernelPlan& kernels,
                                     ConstFieldView gradient,
                                     ConstFieldView viscosity, std::size_t axis,
                                     Int3 face, std::uint8_t component) noexcept {
  const Int3 left = viscous_offset(face, axis, -1);
  const auto normal = axis == 0U ? face.x : (axis == 1U ? face.y : face.z);
  const double mu_face = kernels.geometry_kind() == GeometryKind::uniform
      ? metric_interpolate_material_face<true>(kernels, axis, normal,
            viscosity.unchecked(left, 0U), viscosity.unchecked(face, 0U))
      : metric_interpolate_material_face<false>(kernels, axis, normal,
            viscosity.unchecked(left, 0U), viscosity.unchecked(face, 0U));
  const auto interpolate_gradient = [&](std::uint8_t c) noexcept {
    return kernels.geometry_kind() == GeometryKind::uniform
        ? metric_interpolate_face<true>(kernels, axis, normal,
              gradient.unchecked(left, c), gradient.unchecked(face, c))
        : metric_interpolate_face<false>(kernels, axis, normal,
              gradient.unchecked(left, c), gradient.unchecked(face, c));
  };
  const double divergence = interpolate_gradient(0U) +
                            interpolate_gradient(4U) + interpolate_gradient(8U);
  const double transpose = interpolate_gradient(static_cast<std::uint8_t>(3U * axis + component));
  return mu_face * (transpose - (axis == component ? (2.0 / 3.0) * divergence : 0.0));
}

inline double viscous_face_traction_area(const CartesianKernelPlan& kernels,
                                         ConstFieldView velocity,
                                         ConstFieldView gradient,
                                         ConstFieldView viscosity,
                                         CartesianAxis axis, Int3 face,
                                         std::uint8_t component) noexcept {
  const Int3 left = viscous_offset(face, static_cast<std::size_t>(axis), -1);
  const double ml = viscosity.unchecked(left, 0U), mr = viscosity.unchecked(face, 0U);
  // The physical zero-viscosity limit is needed by inviscid/term-isolation
  // fixtures. Production Newtonian materials remain strictly positive.
  const double weight = ml >= 0.0 && mr >= 0.0 && (ml == 0.0 || mr == 0.0)
      ? 0.0 : positive_transmissibility(kernels, viscosity, axis, face);
  return weight * (velocity.unchecked(face, component) - velocity.unchecked(left, component)) +
         viscous_cross_traction(kernels, gradient, viscosity,
                                static_cast<std::size_t>(axis), face, component) *
             face_area(kernels, axis, face);
}

// Discrete product rule: Phi_V = sum_faces (U_face-U_cell).(tau.n A).
// Hence U_cell.(div tau)_V + Phi_V is exactly the face work divergence.
// This includes shear modes invisible to a central cell-gradient square.
inline Status cartesian_viscous_heating(const CartesianKernelPlan& kernels,
                                        ConstFieldView velocity,
                                        ConstFieldView gradient,
                                        ConstFieldView viscosity, Int3 cell,
                                        double& out, bool total_energy_work = false) noexcept {
  double integral = 0.0;
  for (std::size_t a = 0U; a < 3U; ++a) {
    const auto axis = static_cast<CartesianAxis>(a);
    for (int side : {0, 1}) {
      const Int3 face = viscous_offset(cell, a, side);
      const Int3 left = viscous_offset(face, a, -1);
      const auto normal = a == 0U ? face.x : (a == 1U ? face.y : face.z);
      for (std::uint8_t c = 0U; c < 3U; ++c) {
        const double face_velocity = interpolate_face(kernels, axis, normal,
            velocity.unchecked(left, c), velocity.unchecked(face, c));
        const double traction = viscous_face_traction_area(
            kernels, velocity, gradient, viscosity, axis, face, c);
        integral += (side == 1 ? 1.0 : -1.0) *
                    (face_velocity - (total_energy_work ? 0.0 :
                                         velocity.unchecked(cell, c))) * traction;
      }
    }
  }
  const double value = integral / cell_volume(kernels, cell);
  if (!std::isfinite(value)) return {StatusCode::numerical_failure, 1432U};
  // Do not clip a local value: that would destroy the discrete product rule.
  out = value;
  return {};
}

}  // namespace hundun::v04::detail
