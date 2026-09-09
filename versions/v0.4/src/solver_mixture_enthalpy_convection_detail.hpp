// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_conservative_energy_detail.hpp"

namespace hundun::v04::detail {

// A nonlinear limiter does not commute with summing species enthalpies.
// At each face use a common temperature anchor T_f and reconstruct
//   h_f = R_h[h - sum_s (h_s(T_f)-h_dep(T_f))*Y_s]
//         + sum_s (h_s(T_f)-h_dep(T_f))*R_s[Y_s].
// Species reference shifts then ride the ACTUAL species convective flux.
// At constant T the reconstructed thermal field is constant for any number
// of independently varying species. No physical table or cell state changes.
class MixtureEnthalpyConvection {
 public:
  static bool active(const EnthalpyEquationPlan& plan) noexcept {
    return plan.conservative_total_energy_ && !plan.species_specs_.empty();
  }
  static bool valid_temperature(const EnthalpyEquationPlan& plan,
                                ConstFieldView temperature) noexcept {
    return temperature.field==plan.temperature_ &&
        valid_cell_view(temperature,plan.cells_,0U,1U,1U);
  }
  static Status validate(const EnthalpyEquationPlan& plan,
                         const EquationStateView& state) noexcept {
    if(!active(plan)) return {};
    const auto cells=plan.cells_;
    const auto reach=static_cast<std::uint8_t>(std::max(
        plan.convection_reach_,static_cast<std::uint8_t>(plan.species_convection_==ConvectionScheme::central2 ? 1U : 2U)));
    if(plan.kernels_==nullptr || plan.fingerprint_==0U ||
        state.independent_species.size!=plan.species_specs_.size() ||
        state.independent_species.data==nullptr ||
        !valid_cell_view(state.enthalpy.trial,cells,0U,1U,plan.convection_reach_) ||
        !valid_cell_view(state.temperature.trial,cells,0U,1U,1U))
      return {StatusCode::invalid_plan,1431U};
    for(std::size_t s=0U;s<plan.species_specs_.size();++s)
      if(state.independent_species.data[s].trial.field!=plan.species_specs_[s].field ||
          !valid_cell_view(state.independent_species.data[s].trial,cells,0U,1U,reach))
        return {StatusCode::invalid_plan,1431U};
    return {};
  }

  // Mass-rate times the difference from ordinary scalar-h reconstruction.
  // The correction is a conservative face-flux difference, not a heat source
  // or a post-solve global energy adjustment. Inputs have passed validate().
  static Status face_delta(const EnthalpyEquationPlan& plan,
      const EquationStateView& state,CartesianAxis axis,Int3 face,
      double mass_rate,double& out) noexcept {
    return face_delta_impl(plan,state.enthalpy.trial,state.temperature.trial,
        [&](std::size_t s) { return state.independent_species.data[s].trial; },
        axis,face,mass_rate,out);
  }

  static Status add_predictor_correction(const EnthalpyEquationPlan& plan,
      ConstFieldView enthalpy,ConstFieldView temperature,
      Span<const ConstFieldView> species,ConstFaceFluxView flux,
      const IbmEquationInterfacePlan* immersed,KernelBox box,FieldView rate) noexcept {
    if(!active(plan)) return {};
    if(species.size!=plan.species_specs_.size() || species.data==nullptr)
      return {StatusCode::invalid_plan,1431U};
    return add_impl(plan,enthalpy,temperature,
        [&](std::size_t s) { return species.data[s]; },flux,immersed,box,rate);
  }

  static Status add_correction(const EnthalpyEquationPlan& plan,
      const EquationStateView& state,ConstFaceFluxView flux,
      const IbmEquationInterfacePlan* immersed,KernelBox box,FieldView rate) noexcept {
    if(!active(plan)) return {};
    const auto status=validate(plan,state);
    if(!status) return status;
    return add_impl(plan,state.enthalpy.trial,state.temperature.trial,
        [&](std::size_t s) { return state.independent_species.data[s].trial; },
        flux,immersed,box,rate);
  }

 private:
  template<class SpeciesView>
  static Status face_delta_impl(const EnthalpyEquationPlan& plan,
      ConstFieldView enthalpy,ConstFieldView temperature_field,SpeciesView species,
      CartesianAxis axis,Int3 face,double mass_rate,double& out) noexcept {
    out=0.0;
    if(mass_rate==0.0 || !active(plan)) return {};
    const auto& kernels=*plan.kernels_;
    const auto offset=[&](Int3 c,int n) {
      (axis==CartesianAxis::x ? c.x : axis==CartesianAxis::y ? c.y : c.z)+=n;
      return c;
    };
    const int normal=axis==CartesianAxis::x ? face.x : axis==CartesianAxis::y ? face.y : face.z;
    const double temperature=interpolate_face(kernels,axis,normal,
        temperature_field.unchecked(offset(face,-1),0U),temperature_field.unchecked(face,0U));
    const int first=plan.convection_==ConvectionScheme::central2 ? 1 : 0;
    const int last=plan.convection_==ConvectionScheme::central2 ? 3 : 4;
    std::array<long double,4U> thermal{};
    for(int i=first;i<last;++i) thermal[i]=enthalpy.unchecked(offset(face,i-2),0U);
    long double species_enthalpy=0.0L;
    for(std::size_t s=0U;s<plan.species_specs_.size();++s) {
      double difference=0.0;
      auto status=plan.thermodynamics_.independent_species_enthalpy_difference(s,temperature,difference);
      if(!status) return status;
      if(difference==0.0) continue;
      const auto y=species(s);
      double y_face=0.0;
      status=reconstruct_cartesian_convection_face(kernels,plan.species_convection_,y,0U,axis,face,mass_rate,y_face);
      if(!status) return status;
      species_enthalpy+=static_cast<long double>(difference)*y_face;
      for(int i=first;i<last;++i)
        thermal[i]-=static_cast<long double>(difference)*y.unchecked(offset(face,i-2),0U);
    }
    std::array<double,4U> samples{};
    for(int i=first;i<last;++i) samples[i]=static_cast<double>(thermal[i]);
    const double thermal_face=sampled_convection_face(kernels,plan.convection_,samples,axis,face,mass_rate);
    double ordinary=0.0;
    const auto status=reconstruct_cartesian_convection_face(kernels,plan.convection_,enthalpy,0U,axis,face,mass_rate,ordinary);
    if(!status) return status;
    const double value=static_cast<double>(static_cast<long double>(mass_rate)*
        (static_cast<long double>(thermal_face)+species_enthalpy-ordinary));
    if(!std::isfinite(value)) return {StatusCode::numerical_failure,1432U};
    out=value;
    return {};
  }

  // Add after the ordinary h convection and its prescribed IBM inlet fix.
  // Remove this difference on IBM cut faces so their already authoritative
  // prescribed h_in (or zero wall flux) is retained unchanged.
  template<class SpeciesView>
  static Status add_impl(const EnthalpyEquationPlan& plan,
      ConstFieldView enthalpy,ConstFieldView temperature,SpeciesView species,ConstFaceFluxView flux,
      const IbmEquationInterfacePlan* immersed,KernelBox box,FieldView rate) noexcept {
    Status status;
    if(!valid_cell_view(enthalpy,plan.cells_,0U,1U,plan.convection_reach_) ||
        !valid_cell_view(temperature,plan.cells_,0U,1U,1U) ||
        !valid_kernel_box(box,plan.cells_) || !valid_cell_view(rate,plan.cells_,0U,1U) ||
        !valid_flux_view(flux,plan.cells_,flux.revision) ||
        (immersed!=nullptr && (immersed->fingerprint()==0U || immersed->topology_==nullptr)) ||
        field_views_overlap(enthalpy,as_const(rate)) ||
        field_views_overlap(temperature,as_const(rate)))
      return {StatusCode::invalid_plan,1431U};
    const auto reach=static_cast<std::uint8_t>(std::max(
        plan.convection_reach_,static_cast<std::uint8_t>(plan.species_convection_==ConvectionScheme::central2 ? 1U : 2U)));
    for(std::size_t s=0U;s<plan.species_specs_.size();++s)
      if(!valid_cell_view(species(s),plan.cells_,0U,1U,reach) ||
          species(s).field!=plan.species_specs_[s].field ||
          field_views_overlap(species(s),as_const(rate)))
        return {StatusCode::invalid_plan,1431U};
    const std::array<ConstFaceFieldView,3U> faces{flux.x,flux.y,flux.z};
    const Int3 end{box.begin.x+box.cells.x,box.begin.y+box.cells.y,box.begin.z+box.cells.z};
    const auto activity=immersed==nullptr ? Span<const std::uint8_t>{} : immersed->cell_activity();
    for(int z=box.begin.z;z<end.z;++z) for(int y=box.begin.y;y<end.y;++y) for(int x=box.begin.x;x<end.x;++x) {
      const Int3 c{x,y,z};
      const auto index=(static_cast<std::size_t>(z)*plan.cells_.y+y)*plan.cells_.x+x;
      if(activity.size!=0U && activity.data[index]==0U) continue;
      long double difference=0.0L;
      for(unsigned a=0U;a<3U;++a) for(unsigned side=0U;side<2U;++side) {
        Int3 face=c; (a==0U ? face.x : a==1U ? face.y : face.z)+=side;
        double delta=0.0;
        status=face_delta_impl(plan,enthalpy,temperature,species,static_cast<CartesianAxis>(a),face,faces[a].unchecked(face),delta);
        if(!status) return status;
        difference+=(side==0U ? -1.0L : 1.0L)*delta;
      }
      const double value=rate.unchecked(c,0U)+static_cast<double>(difference/cell_volume(*plan.kernels_,c));
      if(!std::isfinite(value)) return {StatusCode::numerical_failure,1432U};
      rate.unchecked(c,0U)=value;
    }
    if(immersed!=nullptr) {
      const auto links=immersed->topology_->links();
      for(std::size_t i=0U;i<links.size;++i) {
        const auto c=links.data[i].fluid_local_index, s=links.data[i].solid_local_index;
        if(c.x<box.begin.x || c.x>=end.x || c.y<box.begin.y || c.y>=end.y || c.z<box.begin.z || c.z>=end.z) continue;
        const unsigned a=c.x!=s.x ? 0U : c.y!=s.y ? 1U : 2U;
        const bool positive=a==0U ? s.x>c.x : a==1U ? s.y>c.y : s.z>c.z;
        const auto face=positive ? s : c;
        double delta=0.0;
        status=face_delta_impl(plan,enthalpy,temperature,species,static_cast<CartesianAxis>(a),face,faces[a].unchecked(face),delta);
        if(!status) return status;
        rate.unchecked(c,0U)-=(positive ? delta : -delta)/cell_volume(*plan.kernels_,c);
      }
    }
    return {};
  }
};
} // namespace hundun::v04::detail
