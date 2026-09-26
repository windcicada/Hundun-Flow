// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "solver_equation_detail.hpp"
#include <vector>

namespace hundun::v04::detail {
// One component's immutable material during a bounded scalar solve. Refilled
// at every outer correction; it never retains solution, source or limiter data.
class ScalarDiffusionFaces {
 public:
  Status allocate(Int3 cells) {
    if(!valid_cells(cells))return {StatusCode::invalid_plan,1452};
    cells_=cells;
    const auto sy=std::size_t(cells.x)+2U, sz=sy*(std::size_t(cells.y)+2U);
    material_.resize(sz*(std::size_t(cells.z)+2U));
    gamma_.base=material_.data()+1U+sy+sz;
    gamma_.interior=cells;gamma_.ghosts={1,1,1};gamma_.components=1;
    gamma_.stride_y=sy;gamma_.stride_z=sz;gamma_.component_stride=material_.size();
    gamma_.storage_identity=reinterpret_cast<StorageIdentity>(material_.data());
    gamma_.revision_domain=gamma_.storage_identity;
    auto status=FaceFluxStorage::allocate_workspace(cells,1U,storage_);
    if(status)status=storage_.workspace_view(0U,1U,faces_);
    return status;
  }
  Status prepare(const CartesianKernelPlan& kernels,ConstFieldView input) noexcept {
    fingerprint_=0;
    if(!valid_cell_view(input,cells_,0,1,1) || !same_cells(kernels.cells(),cells_))
      return {StatusCode::invalid_plan,1452};
    gamma_.field=input.field;gamma_.revision=input.revision;
    for(int z=-1;z<=cells_.z;++z)for(int y=-1;y<=cells_.y;++y)for(int x=-1;x<=cells_.x;++x) {
      if((x<0 || x>=cells_.x)+(y<0 || y>=cells_.y)+(z<0 || z>=cells_.z)>1)continue;
      gamma_.unchecked({x,y,z},0)=input.unchecked({x,y,z},0);
    }
    auto status=fill_equation_face_coefficients<CartesianAxis::x>(kernels,as_const(gamma_),{{0,0,0},cells_},faces_.x,1453);
    if(status)status=fill_equation_face_coefficients<CartesianAxis::y>(kernels,as_const(gamma_),{{0,0,0},cells_},faces_.y,1453);
    if(status)status=fill_equation_face_coefficients<CartesianAxis::z>(kernels,as_const(gamma_),{{0,0,0},cells_},faces_.z,1453);
    if(status)fingerprint_=kernels.fingerprint();
    return status;
  }
  bool matches(const CartesianKernelPlan& kernels,ConstFieldView gamma) const noexcept {
    return fingerprint_ && fingerprint_==kernels.fingerprint() && gamma.base==gamma_.base &&
        gamma.field==gamma_.field && gamma.revision==gamma_.revision &&
        gamma.storage_identity==gamma_.storage_identity && gamma.revision_domain==gamma_.revision_domain &&
        gamma.components==1 && same_cells(gamma.interior,cells_) && same_cells(gamma.ghosts,gamma_.ghosts) &&
        gamma.stride_y==gamma_.stride_y && gamma.stride_z==gamma_.stride_z &&
        gamma.component_stride==gamma_.component_stride;
  }
  ConstFieldView material() const noexcept {return as_const(gamma_);}
  double face(unsigned axis,Int3 cell) const noexcept {
    return (axis==0 ? faces_.x : axis==1 ? faces_.y : faces_.z).unchecked(cell);
  }
  double diagonal(Int3 cell,std::uint8_t blocked=0) const noexcept {
    double value=0;
    for(unsigned d=0;d<6;++d) {
      if(blocked&(1U<<d))continue;
      auto at=cell;if(d%2)++(d/2==0 ? at.x : d/2==1 ? at.y : at.z);
      value+=face(d/2,at);
    }
    return value;
  }
  void copy_faces(EquationSystemView output) const noexcept {
    const std::array<FaceFieldView,3> out{output.x_coefficient,output.y_coefficient,output.z_coefficient};
    for(unsigned a=0;a<3;++a) {
      auto end=cells_;++(a==0 ? end.x : a==1 ? end.y : end.z);
      for(int z=0;z<end.z;++z)for(int y=0;y<end.y;++y)for(int x=0;x<end.x;++x)
        out[a].unchecked({x,y,z})=face(a,{x,y,z});
    }
  }
  std::uint64_t owned_payload_bytes() const noexcept {
    return material_.capacity()*sizeof(double)+storage_.counters().aligned_payload_bytes;
  }
 private:
  Int3 cells_{};
  std::vector<double> material_;
  FieldView gamma_{};
  FaceFluxStorage storage_;
  FaceFluxView faces_{};
  PlanFingerprint fingerprint_{};
};

Status cached_mixture_transport(const CartesianKernelPlan&,const MixtureTransportFaces&,
    const ScalarDiffusionFaces&,ConstFaceFluxView,const KernelInvocation&) noexcept;
} // namespace hundun::v04::detail
