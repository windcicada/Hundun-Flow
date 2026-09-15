// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "solver_ibm_scalar_transport_detail.hpp"
#include "solver_scalar_boundary_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"
#include <array>
#include <cmath>

namespace hundun::v04::detail {
struct StatisticalEnergyBalanceView {
  std::array<ConstFaceFieldView,3> correction{};
  ConstFieldView source{};
  RevisionToken generation{};
};

// Reusable face workspace with attempt-scoped contents. Freeze the actual pre-chemistry field fluxes,
// then repeat pressure corrections from that same frozen state. Only the
// enclosing accepted transaction persists physical fields and balance history.
class StatisticalEnergyLedger {
 public:
  Status reserve(Int3 cells) noexcept {
    discard_attempt();
    if (owned_payload_bytes())
      return cells.x == cells_.x && cells.y == cells_.y && cells.z == cells_.z
          ? Status{} : invalid();
    auto status = FaceFluxStorage::allocate_workspace(cells,5,storage_);
    if (status) cells_ = cells;
    return status;
  }
  std::size_t owned_payload_bytes() const noexcept {
    return static_cast<std::size_t>(storage_.counters().aligned_payload_bytes);
  }
  void discard_attempt() noexcept {
    bank_ = {}; enthalpy_ = nullptr; generation_ = 0;
    fields_ = frozen_fields_ = corrected_fields_ = 0; ready_ = false;
  }
  Status initialize(Int3 cells,std::size_t fields,RevisionToken generation,
      const EnthalpyEquationPlan& enthalpy) noexcept {
    discard_attempt();
    if(!fields || !generation)return invalid();
    auto status=reserve(cells);
    for(unsigned b=0;b<5 && status;++b)status=storage_.workspace_view(b,generation,bank_[b]);
    if(!status)return status;
    cells_=cells;fields_=fields;generation_=generation;enthalpy_=&enthalpy;
    frozen_fields_=corrected_fields_=0;ready_=false;
    visit([&](unsigned a,Int3 f){for(unsigned b=0;b<5;++b)face(b,a).unchecked(f)=0;});
    return {};
  }
  Status freeze(const CartesianKernelPlan& kernels,ConstFieldView h,
      ConstFieldView gamma,ConstFaceFluxView flux,const MixtureTransportFaces* mixture,
      const IbmEquationInterfacePlan* ibm) noexcept {
    if(!generation_ || frozen_fields_>=fields_)return invalid();
    const MixtureTransportFaces centered{as_const(bank_[3].x),as_const(bank_[3].y),as_const(bank_[3].z),flux.revision,generation_};
    auto status=export_flux(kernels,h,gamma,flux,mixture ? *mixture : centered,ibm);
    if(!status)return status;
    const auto& selected=mixture ? *mixture : centered;
    const std::array<ConstFaceFieldView,3> extra{selected.x,selected.y,selected.z};
    visit([&](unsigned a,Int3 f){
      face(0,a).unchecked(f)+=face(2,a).unchecked(f)/fields_;
      face(4,a).unchecked(f)+=extra[a].unchecked(f)/fields_;
    });
    ++frozen_fields_;return {};
  }
  // A common pressure-work increment inherits the mean frozen field operator
  // and the same upwind pressure-flux increment. A limiter evaluated on total
  // h would act on delta-h after anchor subtraction, creating a discontinuous
  // flux when grad(h) vanishes while grad(delta-h) stays finite.
  Status correction_mixture(const CartesianKernelPlan& kernels,
      const BoundaryPlan& boundary,ConstFaceFluxView old_flux,
      ConstFaceFluxView flux,FaceFluxView output,RevisionToken linearization,
      MixtureTransportFaces& out) const noexcept {
    if(!generation_ || frozen_fields_!=fields_ || !linearization ||
        !valid_flux_view(old_flux,cells_,old_flux.revision) ||
        !valid_flux_view(flux,cells_,flux.revision) || !valid_flux_view(output,cells_))return invalid();
    const std::array<ConstFaceFieldView,3> saved{as_const(bank_[4].x),as_const(bank_[4].y),as_const(bank_[4].z)};
    for(unsigned a=0;a<3;++a)for(unsigned b=0;b<3;++b) {
      const auto target=select(output,static_cast<CartesianAxis>(a));
      const auto axis=static_cast<CartesianAxis>(b);
      if(face_views_overlap(target,select(old_flux,axis)) || face_views_overlap(target,select(flux,axis)) ||
          face_views_overlap(target,saved[b]) || (b<a && face_views_overlap(target,select(output,axis))))return invalid();
    }
    Status status;
    const auto value=[&](unsigned a,Int3 f) {
      const auto axis=static_cast<CartesianAxis>(a);
      const double dm=select(flux,axis).unchecked(f)-select(old_flux,axis).unchecked(f);
      const int normal=a==0 ? f.x : a==1 ? f.y : f.z;
      const int extent=a==0 ? cells_.x : a==1 ? cells_.y : cells_.z;
      double added=0;
      bool face_donor=false;
      if((normal==0 && dm>0) || (normal==extent && dm<0)) {
        const BoundaryFacePlan* rule{};
        const auto found=boundary.face(static_cast<CartesianFace>(2*a+(normal==extent)),rule);
        face_donor=found && rule && rule->local_owner && !rule->periodic;
      }
      if(dm!=0 && !face_donor && !kernels.physical_inlet_material(a,normal)) {
        const double w=interpolate_face(kernels,axis,normal,1.,0.);
        added=std::abs(dm)*(dm>0 ? 1.-w : w);
      }
      return saved[a].unchecked(f)+added;
    };
    visit([&](unsigned a,Int3 f){const double v=value(a,f);if(!std::isfinite(v) || v<0)status=invalid();});
    if(!status)return status;
    visit([&](unsigned a,Int3 f){select(output,static_cast<CartesianAxis>(a)).unchecked(f)=value(a,f);});
    out={as_const(output.x),as_const(output.y),as_const(output.z),flux.revision,linearization};
    return {};
  }
  Status start_correction() noexcept {
    if(frozen_fields_!=fields_ || !generation_)return invalid();
    corrected_fields_=0;ready_=false;
    visit([&](unsigned a,Int3 f){face(1,a).unchecked(f)=face(0,a).unchecked(f);});
    return {};
  }
  Status add_pressure(const BoundaryPlan& boundary,ConstFieldView h,
      ConstFieldView velocity,ConstFaceFluxView old_flux,ConstFaceFluxView flux) noexcept {
    if(!generation_ || corrected_fields_>=fields_ || !valid_cell_view(h,cells_,0,1,1) ||
        !valid_flux_view(old_flux,cells_,old_flux.revision) || !valid_flux_view(flux,cells_,flux.revision))return invalid();
    Status status;
    visit([&](unsigned a,Int3 f){
      if(!status)return;
      const auto axis=static_cast<CartesianAxis>(a);
      const double dm=select(flux,axis).unchecked(f)-select(old_flux,axis).unchecked(f);
      if(dm==0)return;
      Int3 donor=f;if(dm>0)--(a==0 ? donor.x : a==1 ? donor.y : donor.z);
      const double value=face(1,a).unchecked(f)+dm*scalar_upwind_donor(boundary,h,donor,velocity)/fields_;
      if(!std::isfinite(value))status={StatusCode::numerical_failure,17858};
      else face(1,a).unchecked(f)=value;
    });
    if(status)++corrected_fields_;
    return status;
  }
  Status prepare(const CartesianKernelPlan& kernels,ConstFieldView anchor,
      ConstFieldView gamma,ConstFaceFluxView flux,const MixtureTransportFaces& mixture,
      const IbmEquationInterfacePlan* ibm) noexcept {
    ready_=false;
    if(!generation_ || corrected_fields_!=fields_)return invalid();
    auto status=export_flux(kernels,anchor,gamma,flux,mixture,ibm);
    if(!status)return status;
    visit([&](unsigned a,Int3 f){face(3,a).unchecked(f)=face(1,a).unchecked(f)-face(2,a).unchecked(f);});
    ready_=true;return {};
  }
  Status add_residual(FieldView residual,Span<const std::uint8_t> active) const noexcept {
    if(!ready_ || !valid_cell_view(as_const(residual),cells_,0,1,0) ||
        (active.size && (!active.data || active.size!=std::size_t(cells_.x)*cells_.y*cells_.z)))return invalid();
    std::size_t i{};
    for(int z=0;z<cells_.z;++z)for(int y=0;y<cells_.y;++y)for(int x=0;x<cells_.x;++x,++i) {
      if(active.size && !active.data[i])continue;
      const Int3 c{x,y,z};long double increment{};
      for(unsigned a=0;a<3;++a) {
        Int3 upper=c;++(a==0 ? upper.x : a==1 ? upper.y : upper.z);
        const auto value=select(as_const(bank_[3]),static_cast<CartesianAxis>(a));
        increment+=static_cast<long double>(value.unchecked(upper))-value.unchecked(c);
      }
      const double value=residual.unchecked(c,0)+static_cast<double>(increment);
      if(!std::isfinite(value))return {StatusCode::numerical_failure,17858};
      residual.unchecked(c,0)=value;
    }
    return {};
  }
  StatisticalEnergyBalanceView balance(ConstFieldView source) const noexcept {
    if(!ready_)return {};
    return {{as_const(bank_[3].x),as_const(bank_[3].y),as_const(bank_[3].z)},source,generation_};
  }
 private:
  static Status invalid() noexcept {return {StatusCode::invalid_plan,17859};}
  FaceFieldView face(unsigned b,unsigned a) noexcept {return select(bank_[b],static_cast<CartesianAxis>(a));}
  template<class F> void visit(F&& f) const noexcept {
    for(unsigned a=0;a<3;++a) {
      auto end=cells_;++(a==0 ? end.x : a==1 ? end.y : end.z);
      for(int z=0;z<end.z;++z)for(int y=0;y<end.y;++y)for(int x=0;x<end.x;++x)f(a,{x,y,z});
    }
  }
  Status export_flux(const CartesianKernelPlan& kernels,ConstFieldView h,ConstFieldView gamma,
      ConstFaceFluxView flux,const MixtureTransportFaces& mixture,const IbmEquationInterfacePlan* ibm) noexcept {
    std::array<FaceFieldView,3> output{bank_[2].x,bank_[2].y,bank_[2].z};
    auto status=form_cartesian_mixture_transport_flux(kernels,mixture,gamma,flux,h,output);
    if(status && ibm)status=IbmScalarTransport::constrain_flux(*ibm,{IbmScalarTransport::Quantity::enthalpy,0},flux,output);
    if(status && enthalpy_->has_prescribed_heat_flux()) {
      const auto active=ibm ? ibm->cell_activity() : Span<const std::uint8_t>{};
      for(unsigned f=0;f<6 && status;++f) {
        double outward{};if(!enthalpy_->prescribed_heat_flux(static_cast<CartesianFace>(f),outward))continue;
        const unsigned a=f/2;const bool high=f%2;
        const int n=a==0 ? cells_.x : a==1 ? cells_.y : cells_.z;
        auto end=cells_;(a==0 ? end.x : a==1 ? end.y : end.z)=1;
        for(int z=0;z<end.z;++z)for(int y=0;y<end.y;++y)for(int x=0;x<end.x;++x) {
          Int3 face{x,y,z},owner=face;
          (a==0 ? face.x : a==1 ? face.y : face.z)=high ? n : 0;
          (a==0 ? owner.x : a==1 ? owner.y : owner.z)=high ? n-1 : 0;
          const auto i=std::size_t(owner.x)+std::size_t(cells_.x)*(owner.y+std::size_t(cells_.y)*owner.z);
          if(active.size && !active.data[i])continue;
          Int3 lower=face;--(a==0 ? lower.x : a==1 ? lower.y : lower.z);
          const auto axis=static_cast<CartesianAxis>(a);
          const double correction=positive_transmissibility(kernels,gamma,axis,face)*
              (h.unchecked(face,0)-h.unchecked(lower,0))+(high ? 1. : -1.)*outward*face_area(kernels,axis,face);
          const double value=output[a].unchecked(face)+correction;
          if(!std::isfinite(value))status={StatusCode::numerical_failure,17858};
          else output[a].unchecked(face)=value;
        }
      }
    }
    return status;
  }
  const EnthalpyEquationPlan* enthalpy_{};
  FaceFluxStorage storage_;
  std::array<FaceFluxView,5> bank_{};
  Int3 cells_{};
  std::size_t fields_{},frozen_fields_{},corrected_fields_{};
  RevisionToken generation_{};
  bool ready_{};
};
} // namespace hundun::v04::detail
