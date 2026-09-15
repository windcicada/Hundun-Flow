// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "solver_equation_detail.hpp"

namespace hundun::v04::detail {

inline double mixture_face_extra(const MixtureTransportFaces* mixture,
                                 unsigned axis, Int3 face) noexcept {
  if (mixture == nullptr) return 0.0;
  const std::array<ConstFaceFieldView,3> values{mixture->x,mixture->y,mixture->z};
  return values[axis].unchecked(face);
}

inline double mixture_diffusion_diagonal(const MixtureTransportFaces* mixture,
                                         Int3 cell) noexcept {
  if (mixture == nullptr) return 0.0;
  double diagonal{};
  for (unsigned a=0; a<3; ++a) {
    auto upper=cell;
    (a==0 ? upper.x : a==1 ? upper.y : upper.z)++;
    diagonal += mixture_face_extra(mixture,a,cell)+mixture_face_extra(mixture,a,upper);
  }
  return diagonal;
}

inline bool valid_mixture_transport(const MixtureTransportFaces* mixture,
    Int3 cells, RevisionToken flux, EquationSystemView system) noexcept {
  if (mixture == nullptr) return true;
  if (mixture->linearization == 0U || mixture->face_flux != flux) return false;
  const std::array<ConstFaceFieldView,3> faces{mixture->x,mixture->y,mixture->z};
  const std::array<FaceFieldView,3> output{
      system.x_coefficient,system.y_coefficient,system.z_coefficient};
  for (unsigned a=0; a<3; ++a) {
    if (!valid_equation_face_view(faces[a],static_cast<CartesianAxis>(a),cells) ||
        cell_face_views_overlap(as_const(system.diagonal),faces[a]) ||
        cell_face_views_overlap(as_const(system.rhs),faces[a]) ||
        cell_face_views_overlap(as_const(system.residual),faces[a])) return false;
    for (const auto& field : output)
      if (field.base && face_views_overlap(field,faces[a])) return false;
  }
  return true;
}

inline void add_mixture_face_coefficients(const MixtureTransportFaces* mixture,
    KernelBox box, EquationSystemView system) noexcept {
  if (mixture == nullptr || system.x_coefficient.base == nullptr) return;
  const std::array<FaceFieldView,3> coefficients{
      system.x_coefficient,system.y_coefficient,system.z_coefficient};
  for (unsigned a=0; a<3; ++a) {
    Int3 end{box.begin.x+box.cells.x,box.begin.y+box.cells.y,box.begin.z+box.cells.z};
    (a==0 ? end.x : a==1 ? end.y : end.z)++;
    for (int z=box.begin.z; z<end.z; ++z)
      for (int y=box.begin.y; y<end.y; ++y)
        for (int x=box.begin.x; x<end.x; ++x) {
          const Int3 face{x,y,z};
          coefficients[a].unchecked(face) += mixture_face_extra(mixture,a,face);
        }
  }
}

} // namespace hundun::v04::detail
