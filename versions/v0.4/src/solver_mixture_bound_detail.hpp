// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "solver_equation_detail.hpp"
#include "hundun/v04_ibm.hpp"

namespace hundun::v04::detail {
// Repair a frozen high-order candidate by raising the COMMON face conductance
// around inadmissible cells to its upwind value. Both face owners exchange the
// cell mask before this call. Every species and h then solve again from the
// same accepted history; no cell value or conserved inventory is clipped.
inline Status constrain_mixture_bounds(const CartesianKernelPlan& kernels,
    ConstFieldView marked, ConstFieldView gamma, ConstFaceFluxView flux,
    FaceFluxView extra, const IbmEquationInterfacePlan* immersed) noexcept {
  const auto cells=kernels.cells();
  if(!valid_cell_view(marked,cells,0,1,1) || !valid_cell_view(gamma,cells,0,1,1))
    return {StatusCode::invalid_plan,17876};
  const std::array<ConstFaceFieldView,3> mass{flux.x,flux.y,flux.z};
  const std::array<FaceFieldView,3> faces{extra.x,extra.y,extra.z};
  // Inputs and the writable coefficient storage belong to the immediately
  // preceding prepare_cartesian_mixture_transport call.
  for(unsigned a=0;a<3;++a) {
    const auto axis=static_cast<CartesianAxis>(a);
    for(int z=0;z<faces[a].extents.z;++z)
      for(int y=0;y<faces[a].extents.y;++y)
        for(int x=0;x<faces[a].extents.x;++x) {
          const Int3 face{x,y,z};auto lower=face;
          (a==0 ? lower.x : a==1 ? lower.y : lower.z)--;
          if(marked.unchecked(lower,0)==0 && marked.unchecked(face,0)==0)continue;
          const auto normal=a==0 ? x : a==1 ? y : z;
          const double phi=mass[a].unchecked(face);double imposed{};
          if(phi==0 || kernels.physical_inlet_material(a,normal) ||
              (immersed && immersed->prescribed_face_flux(axis,face,imposed)))continue;
          const double weight=interpolate_face(kernels,axis,normal,1.,0.);
          const double physical=positive_transmissibility(kernels,gamma,axis,face);
          const double upwind=(phi>0 ? 1-weight : weight)*std::abs(phi);
          const double value=std::max(faces[a].unchecked(face),upwind-physical);
          if(!std::isfinite(value) || value<0)return {StatusCode::numerical_failure,17877};
          faces[a].unchecked(face)=value;
        }
  }
  return {};
}
} // namespace hundun::v04::detail
