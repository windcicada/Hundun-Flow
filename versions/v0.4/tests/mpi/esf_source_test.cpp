// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "../support/candidate_boundary_fixture.hpp"
#include "../support/product_fixture.hpp"
#include "../../src/solver_equation_detail.hpp"
#include "../../src/core_esf_detail.hpp"
#include "../../src/solver_esf_flux_detail.hpp"
#include <iomanip>
#include <iostream>

using namespace hundun::v04;
using namespace hundun::v04::test;

bool run(unsigned nf) {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;spec.cells_per_axis=9;
  if(!fixture.initialize(MPI_COMM_WORLD,spec))return false;
  const auto cells=fixture.patch.cells;
  constexpr double dt=.02,pressure=1e5,mu=1.8e-5;
  const double R=kUniversalGasConstant/28,cp=1000,rho=pressure/(R*300);
  const double base_h=cp*(300-298.15)+30000;
  chemistry::detail::AnalyticIsomerBackend gas(2,false,cp);
  auto model=product_model(fixture.geometry.global_cells());
  model.schemes.species=model.schemes.enthalpy=ConvectionScheme::central2;
  const auto air=model.thermophysics.species.front();
  model.thermophysics.species.assign(2,air);
  model.thermophysics.species[0].stable_name="A";model.thermophysics.species[1].stable_name="B";
  for(unsigned s=0;s<2;++s) {
    auto& property=model.thermophysics.species[s];property.molecular_weight=28;
    property.nasa7_low[0]=property.nasa7_high[0]=cp/R;
    property.nasa7_low[5]=property.nasa7_high[5]=((s==0 ? 100000. : 0.)-cp*298.15)/R;
  }
  model.transported_scalars={{"A",TransportedScalarRole::species,.7,.7}};
  model.reaction.mode=ReactionMode::esf_tpdf;
  model.reaction.mechanism_sha256=gas.gas_identity().mechanism_sha256;
  model.reaction.phase=gas.gas_identity().phase;
  model.reaction.esf=EsfSpec{};model.reaction.esf->fields=nf;
  detail::ProductReactionSources reaction;
  auto status=reaction.configure(model,{&gas,&gas,&gas.closure_identity()},{});
  detail::ProductEsf esf;
  if(status)status=esf.configure(model,reaction,fixture.patch,fixture.geometry.global_cells());
  ThermodynamicsPlan thermo;
  if(status)status=ThermodynamicsPlan::compile(model.thermophysics,
      {model.transported_scalars.data(),model.transported_scalars.size()},thermo);
  if(!status){std::cerr<<"configure "<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  const auto owned=[&](FieldId id,Int3 c,std::uint8_t n,std::uint8_t g,RevisionToken r) {
    return make_field(id,c,n,g,r,10000+2*id+r);
  };
  auto density=owned(800,cells,1,2,7),pi=owned(801,cells,1,2,7);
  auto mean=owned(802,cells,1,2,7),heat=owned(803,cells,1,2,7);
  auto cache=owned(804,cells,4,2,7),aux=owned(805,cells,3,0,9),source=owned(806,cells,1,0,9);
  fill(density,rho);fill(pi,0);fill(mean,.3);fill(heat,base_h);fill(cache,0);
  for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x) {
    cache.view.unchecked({x,y,z},0)=mu/.7;cache.view.unchecked({x,y,z},2)=mu;
  }
  std::vector<OwnedField> accepted(nf),trial(nf);
  std::vector<ConstFieldView> accepted_views(nf);
  std::vector<FieldView> trial_views(nf);
  for(unsigned f=0;f<nf;++f) {
    accepted[f]=owned(810+f,cells,3,2,7);trial[f]=owned(810+f,cells,3,2,9);
    const double shift=.01*(2*int(f)-int(nf)+1);
    for(int z=-2;z<cells.z+2;++z)for(int y=-2;y<cells.y+2;++y)for(int x=-2;x<cells.x+2;++x) {
      accepted[f].view.unchecked({x,y,z},0)=.3+shift;
      accepted[f].view.unchecked({x,y,z},1)=.7-shift;
      accepted[f].view.unchecked({x,y,z},2)=base_h+(100*cp+100000)*shift;
    }
    accepted_views[f]=as_const(accepted[f].view);trial_views[f]=trial[f].view;
  }
  FaceFluxStorage storage;FaceFluxView flux;
  status=FaceFluxStorage::allocate_workspace(cells,1,storage);
  if(status)status=storage.workspace_view(0,81,flux);
  if(!status)return false;
  for(auto face:{flux.x,flux.y,flux.z})
    for(int z=0;z<face.extents.z;++z)for(int y=0;y<face.extents.y;++y)for(int x=0;x<face.extents.x;++x)
      face.unchecked({x,y,z})=0;
  const auto global=fixture.geometry.global_cells();
  const auto mass_flux=[&](int gx) {
    return .04*rho*std::sin(2*std::acos(-1.)*gx/global.x)/(global.y*global.z);
  };
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<=cells.x;++x)
    flux.x.unchecked({x,y,z})=mass_flux(x+fixture.patch.begin.x);
  ConstFieldView mean_view=as_const(mean.view);
  const auto prepare=[&](RevisionToken generation) {
    auto s=esf.prepare_transport(MPI_COMM_WORLD,fixture.kernels,fixture.equations.enthalpy(),
        reaction,thermo,{accepted_views.data(),nf},{trial_views.data(),nf},as_const(cache.view),
        as_const(density.view),as_const(pi.view),{&mean_view,1},as_const(heat.view),
        pressure,as_const(flux),0,dt,4,generation,{&source.view,1});
    if(!s)return s;
    // Uniform per-field coordinates have zero advective derivative and
    // stochastic gradient. Solve their actual implicit IEM scalar rows.
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
      const Int3 c{x,y,z};
      esf::detail::IemSource mixing;
      if(esf::detail::iem_source(fixture.cell_volume(c),mu,0,esf.mixing_cd(),1,.3,mixing)!=portable::Status::success)
        return Status{StatusCode::numerical_failure,1};
      const double factor=dt*mixing.implicit_sink_density/rho;
      const double targets[3]{.3,.7,base_h};
      for(unsigned j=0;j<3;++j) {
        aux.view.unchecked(c,j)=targets[j];
        for(unsigned f=0;f<nf;++f)
          trial_views[f].unchecked(c,j)=(accepted_views[f].unchecked(c,j)+factor*targets[j])/(1+factor);
      }
    }
    return esf.finish_implicit_transport(reaction,thermo,{trial_views.data(),nf},as_const(pi.view),pressure,0,dt,4,generation);
  };
  auto corrected=owned(830,cells,3,1,9),next=owned(831,cells,3,0,9);
  HaloEngine correction_halo;
  const HaloFieldSpec correction_spec{830,1,3};
  if(!correction_halo.reserve(MPI_COMM_WORLD,fixture.patch,{&correction_spec,1},fixture.boundary.halo_topology()))return false;
  const auto original_bytes=esf.owned_bytes();
  double error=0,gap=0;bool passed=true;
  std::vector<double> first;
  for(unsigned attempt=0;attempt<2;++attempt) {
    status=prepare(9+attempt);
    if(status && attempt==1) {
      const auto saved_source=source.storage;
      auto stale=as_const(density.view);++stale.revision;
      const auto rejected=esf.react(reaction,thermo,{trial_views.data(),nf},aux.view,
          stale,as_const(pi.view),pressure,0,dt,4,9+attempt,{&source.view,1});
      passed &= rejected.code==StatusCode::invalid_plan && source.storage==saved_source;
      status=prepare(9+attempt);
    }
    if(status)status=esf.react(reaction,thermo,{trial_views.data(),nf},aux.view,
        as_const(density.view),as_const(pi.view),pressure,0,dt,4,9+attempt,{&source.view,1});
    if(!status){std::cerr<<"reaction "<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
    std::size_t i=0;
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
      const Int3 c{x,y,z};const int gx=x+fixture.patch.begin.x;
      const double basis=rho-dt*(mass_flux(gx+1)-mass_flux(gx))/fixture.cell_volume(c);
      const double delta=.3*std::expm1(-2*dt),expected=basis*delta/dt;
      error=std::max(error,std::abs(source.view.unchecked(c,0)-expected)/std::max(1.,std::abs(expected)));
      gap=std::max(gap,std::abs(source.view.unchecked(c,0)-rho*delta/dt));
      passed &= std::abs(aux.view.unchecked(c,0)-(.3+delta))<2e-14;
      if(attempt==0)first.push_back(source.view.unchecked(c,0));
      else passed &= first[i]==source.view.unchecked(c,0);
    }
    for(unsigned f=0;f<=nf;++f) {
      ConstFieldView snapshot;
      status=esf.correction_seed(f,830,{4,9+attempt,1},snapshot);
      if(!status)return false;
      const auto source_tuple=f==nf ? aux.view : trial_views[f];
      for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
        for(unsigned c=0;c<3;++c)passed &= snapshot.unchecked({x,y,z},c)==source_tuple.unchecked({x,y,z},c);
      const double original=source_tuple.unchecked({0,0,0},0);
      source_tuple.unchecked({0,0,0},0)=-19;
      passed &= snapshot.unchecked({0,0,0},0)==original;
      source_tuple.unchecked({0,0,0},0)=original;
      auto unchanged=snapshot;
      passed &= !esf.correction_seed(f,830,{4,8+attempt,1},unchanged) && unchanged.base==snapshot.base;
      const auto close=[&](FieldView v) {
        for(unsigned face=0;face<6;++face) {
          if(!fixture.local_face_owner(static_cast<CartesianFace>(face)))continue;
          const unsigned axis=face/2;
          const int n=axis==0 ? cells.x : axis==1 ? cells.y : cells.z;
          for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
            Int3 c{x,y,z},g=c;
            int& normal=axis==0 ? g.x : axis==1 ? g.y : g.z;
            if(normal!=(face%2 ? n-1 : 0))continue;
            normal=face%2 ? n : -1;
            for(unsigned s=0;s<3;++s)v.unchecked(g,s)=v.unchecked(c,s);
          }
        }
        return Status{};
      };
      detail::StatisticalFluxReport flux_report;
      status=detail::correct_statistical_flux(fixture.kernels,as_const(density.view),snapshot,
          as_const(flux),as_const(flux),dt,{},f==nf,fixture.boundary,
          as_const(fixture.velocity.view),corrected.view,next.view,correction_halo,178,
          fixture.reductions,close,flux_report);
      if(!status)return false;
      passed &= flux_report.iterations==0;
      for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
        for(unsigned c=0;c<3;++c)passed &= corrected.view.unchecked({x,y,z},c)==snapshot.unchecked({x,y,z},c);
    }
    passed &= esf.owned_bytes()==original_bytes;
    const auto duplicate=esf.react(reaction,thermo,{trial_views.data(),nf},aux.view,
        as_const(density.view),as_const(pi.view),pressure,0,dt,4,9+attempt,{&source.view,1});
    passed &= duplicate.code==StatusCode::invalid_plan;
    ConstFieldView stale_seed;
    passed &= !esf.correction_seed(0,830,{4,9+attempt,1},stale_seed);
  }
  // A nonpositive transported carrier fails in preparation, before the
  // chemistry boundary. A fresh proposal reconstructs the frozen basis.
  const double saved_face=flux.x.unchecked({1,0,0});
  flux.x.unchecked({1,0,0})=1e6;
  const auto invalid_carrier=prepare(11);
  passed &= invalid_carrier.code==StatusCode::numerical_failure;
  ConstFieldView stale_seed;
  passed &= !esf.correction_seed(0,830,{4,10,1},stale_seed);
  flux.x.unchecked({1,0,0})=saved_face;
  status=prepare(12);
  if(status)status=esf.react(reaction,thermo,{trial_views.data(),nf},aux.view,
      as_const(density.view),as_const(pi.view),pressure,0,dt,4,12,{&source.view,1});
  passed &= bool(status);
  std::size_t check=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++check)
    passed &= source.view.unchecked({x,y,z},0)==first[check];
  ConstFieldView current_seed;
  passed &= bool(esf.correction_seed(0,830,{4,12,1},current_seed));
  passed &= !esf.correction_seed(nf+1,830,{4,12,1},current_seed);
  esf.discard();
  passed &= !esf.correction_seed(0,830,{4,12,1},current_seed);
  status=prepare(13);
  if(status)status=esf.react(reaction,thermo,{trial_views.data(),nf},aux.view,
      as_const(density.view),as_const(pi.view),pressure,0,dt,4,13,{&source.view,1});
  passed &= bool(status) && bool(esf.correction_seed(nf,830,{4,13,1},current_seed));
  // Joint physical-mean/auxiliary-pressure closure on the actual reactor
  // output. The signed field0 anchor has a distinct enthalpy and density.
  const Int3 probe{0,0,0};
  const auto saved_auxiliary=aux.storage;
  double anchor_h{};
  for(unsigned f=0;f<nf;++f)anchor_h+=trial_views[f].unchecked(probe,2)/nf;
  aux.view.unchecked(probe,0)=-.0003;
  aux.view.unchecked(probe,1)=1.0003;
  aux.view.unchecked(probe,2)=anchor_h+12000;
  std::vector<std::vector<double>> immutable_fields;
  for(const auto& f:trial)immutable_fields.push_back(f.storage);
  const auto signed_auxiliary=aux.storage;
  std::array<double,3> physical_mean{};
  detail::ProductEsf::DualState state;
  const auto evaluate=[&](double p,double dh,detail::ProductEsf::DualState& output) {
    return esf.evaluate_dual_state_cell(probe,{trial_views.data(),nf},as_const(aux.view),
        p,dh,{},reaction,thermo,{4,13,1},{physical_mean.data(),physical_mean.size()},output);
  };
  passed &= !esf.evaluate_dual_state_cell(probe,{trial_views.data(),nf},as_const(aux.view),
      pressure,0,{},reaction,thermo,{4,13,1},{trial_views[0].base+1,3},state);
  for(unsigned f=0;f<nf;++f)passed &= trial[f].storage==immutable_fields[f];
  constexpr double weight=1.0003,increment=1250;
  status=evaluate(pressure,increment,state);
  if(!status)return false;
  const double auxiliary_T=298.15+(anchor_h+12000+increment)/(weight*cp);
  const double auxiliary_rho=pressure/(weight*R*auxiliary_T);
  const double physical_T=298.15+(anchor_h+increment-100000*physical_mean[0])/cp;
  passed &= std::abs(state.pressure.density_kg_per_m3-auxiliary_rho)<1e-12 &&
      std::abs(state.physical.temperature-physical_T)<1e-10 &&
      std::abs(state.physical.rho-pressure/(R*physical_T))<1e-12 &&
      std::abs(state.physical_enthalpy-anchor_h-increment)<1e-9 &&
      std::abs(state.auxiliary_enthalpy-state.physical_enthalpy-12000)<1e-9 &&
      std::abs(state.physical.rho-state.pressure.density_kg_per_m3)>.01;
  double inverse_density{};
  for(unsigned f=0;f<nf;++f) {
    const double temperature=298.15+(trial_views[f].unchecked(probe,2)+increment
        -100000*trial_views[f].unchecked(probe,0))/cp;
    inverse_density+=R*temperature/(pressure*nf);
  }
  passed &= std::abs(state.statistical_density-1/inverse_density)<1e-12;
  const auto finite_difference=[&](bool enthalpy,double spacing) {
    detail::ProductEsf::DualState plus,minus;
    if(!evaluate(pressure+(enthalpy?0:spacing),increment+(enthalpy?spacing:0),plus) ||
       !evaluate(pressure-(enthalpy?0:spacing),increment-(enthalpy?spacing:0),minus))
      return std::array<double,2>{NAN,NAN};
    const double ep=plus.pressure.density_kg_per_m3*plus.physical_enthalpy-(pressure+(enthalpy?0:spacing));
    const double em=minus.pressure.density_kg_per_m3*minus.physical_enthalpy-(pressure-(enthalpy?0:spacing));
    return std::array<double,2>{(plus.pressure.density_kg_per_m3-minus.pressure.density_kg_per_m3)/(2*spacing),(ep-em)/(2*spacing)};
  };
  PressureEnergyThermoJacobian jacobian;
  passed &= bool(form_pressure_energy_density_jacobian(pressure,state.physical_enthalpy,
      {state.pressure.density_kg_per_m3,state.pressure.density_pressure_derivative,
       state.pressure.density_enthalpy_derivative},jacobian));
  double tangent_error{};
  for(bool enthalpy:{false,true}) {
    const auto full=finite_difference(enthalpy,10),half=finite_difference(enthalpy,5);
    const double expected[2]{enthalpy ? state.pressure.density_enthalpy_derivative : state.pressure.density_pressure_derivative,
        enthalpy ? jacobian.dq_dh_pY : jacobian.dq_dp_hY};
    for(unsigned c=0;c<2;++c) {
      if(!std::isfinite(full[c]) || !std::isfinite(half[c]))tangent_error=INFINITY;
      else tangent_error=std::max(tangent_error,
          std::abs((4*half[c]-full[c])/3-expected[c])/std::abs(expected[c]));
    }
  }
  passed &= std::isfinite(tangent_error) && tangent_error<1e-6;
  // Solve the common public pressure/energy block with a separately known
  // endpoint. Every Newton evaluation uses the same reactor anchor.
  detail::ProductEsf::DualState endpoint;
  if(!evaluate(120000,15000,endpoint))return false;
  const double target_rho=endpoint.pressure.density_kg_per_m3;
  const double target_energy=target_rho*endpoint.physical_enthalpy-120000;
  double solved_p=pressure,solved_dh=0;
  for(unsigned iteration=0;iteration<8;++iteration) {
    if(!evaluate(solved_p,solved_dh,state))return false;
    if(!form_pressure_energy_density_jacobian(solved_p,state.physical_enthalpy,
        {state.pressure.density_kg_per_m3,state.pressure.density_pressure_derivative,
         state.pressure.density_enthalpy_derivative},jacobian))return false;
    PressureEnergyTemporalPoint point;
    point.bdf={1/dt,-1/dt,0,1};point.cell_volume=fixture.cell_volume(probe);
    point.pressure_absolute=solved_p;point.enthalpy=state.physical_enthalpy;
    point.density=state.pressure.density_kg_per_m3;
    point.accepted_pressure_absolute=120000;point.accepted_enthalpy=endpoint.physical_enthalpy;
    point.accepted_density=target_rho;point.target_thermo=jacobian;
    PressureEnergyTemporalLinearization block;
    if(!linearize_pressure_energy_temporal(point,block))return false;
    const double a=block.continuity_pressure,b=block.continuity_enthalpy,
        c=block.energy_pressure,d=block.energy_enthalpy,det=a*d-b*c;
    solved_p+=(-block.continuity_residual*d+b*block.energy_residual)/det;
    solved_dh+=(c*block.continuity_residual-a*block.energy_residual)/det;
  }
  if(!evaluate(solved_p,solved_dh,state))return false;
  const double mass_error=std::abs(state.pressure.density_kg_per_m3-target_rho)/target_rho;
  const double energy_error=std::abs(state.pressure.density_kg_per_m3*state.physical_enthalpy-solved_p-target_energy)/std::max(1.,std::abs(target_energy));
  passed &= mass_error<1e-12 && energy_error<1e-12 && std::abs(solved_p-120000)<1e-7 && std::abs(solved_dh-15000)<1e-7;
  const auto mean_before=physical_mean;
  const auto state_before=state;
  const double late_value=trial_views.back().unchecked(probe,2);
  trial_views.back().unchecked(probe,2)=NAN;
  passed &= !evaluate(pressure,0,state) && physical_mean==mean_before &&
      state.physical_enthalpy==state_before.physical_enthalpy &&
      state.pressure.density_kg_per_m3==state_before.pressure.density_kg_per_m3;
  trial_views.back().unchecked(probe,2)=late_value;
  for(unsigned f=0;f<nf;++f)passed &= immutable_fields[f]==trial[f].storage;
  passed &= aux.storage==signed_auxiliary;
  auto published=owned(840,cells,3,0,13);
  fill(published,-81);
  detail::ProductEsf::DualState publication;
  const auto before_publication=published.storage;
  const auto publish=[&](double dh) {
    return esf.publish_dual_state_cell(probe,{trial_views.data(),nf},aux.view,published.view,
        pressure,dh,{},reaction,thermo,{4,13,1},publication);
  };
  passed &= !publish(-1e20) && published.storage==before_publication && aux.storage==signed_auxiliary;
  for(unsigned f=0;f<nf;++f)passed &= trial[f].storage==immutable_fields[f];
  if(!evaluate(pressure,2500,state) || !publish(2500))return false;
  for(unsigned c=0;c<3;++c)passed &= published.view.unchecked(probe,c)==physical_mean[c];
  for(unsigned f=0;f<nf;++f) {
    for(unsigned c=0;c<2;++c)passed &= trial_views[f].unchecked(probe,c)==immutable_fields[f][
        std::size_t(trial_views[f].base-trial[f].storage.data())+c*trial_views[f].component_stride];
    const auto h_index=std::size_t(trial_views[f].base-trial[f].storage.data())+2*trial_views[f].component_stride;
    passed &= trial_views[f].unchecked(probe,2)==immutable_fields[f][h_index]+2500;
    trial[f].storage=immutable_fields[f];
  }
  passed &= aux.view.unchecked(probe,0)==-.0003 && aux.view.unchecked(probe,1)==1.0003 &&
      publication.physical_enthalpy==state.physical_enthalpy &&
      publication.pressure.density_kg_per_m3==state.pressure.density_kg_per_m3 &&
      aux.view.unchecked(probe,2)==signed_auxiliary[2*aux.view.component_stride]+2500;
  const auto native_publish=[&]() {
    FieldRegistry registry;FieldSchema schema;
    std::vector<FieldId> ids(nf+2);
    for(unsigned f=0;f<nf+2;++f)
      if(!registry.declare_field("dual"+std::to_string(f),3,0,ids[f]))return false;
    if(!registry.freeze(schema))return false;
    std::vector<ArenaFieldRequest> requests;
    for(auto id:ids)requests.push_back({id,cells,{0},FieldLifetime::state_layer});
    ArenaLayout layout;StateLayers layers;
    if(!ArenaLayout::compile(schema,{requests.data(),requests.size()},layout) || !StateLayers::allocate(layout,layers))return false;
    std::vector<FieldView> views(nf+2);
    for(unsigned f=0;f<nf+2;++f)if(!layers.view(StateRole::trial,ids[f],views[f]))return false;
    for(unsigned f=0;f<nf;++f)for(unsigned c=0;c<3;++c)
      views[f].unchecked(probe,c)=trial_views[f].unchecked(probe,c);
    for(unsigned c=0;c<3;++c)views[nf].unchecked(probe,c)=signed_auxiliary[c*aux.view.component_stride];
    const auto result=esf.publish_dual_state_cell(probe,{views.data(),nf},views[nf],views[nf+1],
        pressure,2500,{},reaction,thermo,{4,13,1},publication);
    if(!result)std::cerr<<"native_dual shared_arena status="<<unsigned(result.code)<<'/'<<result.detail<<'\n';
    bool valid=bool(result) && views[0].storage_identity==views[nf+1].storage_identity;
    if(result)for(unsigned c=0;c<3;++c)valid &= views[nf+1].unchecked(probe,c)==physical_mean[c];
    auto aliased_mean=views[0];aliased_mean.field=views[nf+1].field;
    valid &= !esf.publish_dual_state_cell(probe,{views.data(),nf},views[nf],aliased_mean,
        pressure,0,{},reaction,thermo,{4,13,1},publication);
    return valid;
  };
  passed &= native_publish();
  aux.storage=saved_auxiliary;
  int diagnostic_rank;MPI_Comm_rank(MPI_COMM_WORLD,&diagnostic_rank);
  if(diagnostic_rank==0)std::cout<<std::setprecision(17)<<"esf_dual fields="<<nf
      <<" tangent_error="<<tangent_error<<" mass_error="<<mass_error<<" energy_error="<<energy_error
      <<" pressure="<<solved_p<<" dh="<<solved_dh<<'\n';
  esf.commit();
  passed &= !esf.correction_seed(nf,830,{4,13,1},current_seed);
  passed &= bool(evaluate(pressure,0,state));

  double local[2]{error,gap},result[2]{};
  MPI_Allreduce(local,result,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0)std::cout<<std::setprecision(17)<<"esf_source fields="<<nf<<" paired_error="<<result[0]<<" old_weight_gap="<<result[1]<<'\n';
  return passed && result[0]<1e-12 && result[1]>1e-5;
}
int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  int local=run(2) && run(4),global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Finalize();return global ? 0 : 1;
}
