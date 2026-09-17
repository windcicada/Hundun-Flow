// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_app.hpp"
#include "../../src/core_product_freeze_detail.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <vector>
using namespace hundun::v04;
namespace {
bool all(bool okay) {
  int value=okay;
  MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  return value!=0;
}
std::vector<double> values(const RestartSnapshot& state) {
  std::vector<double> result{state.time,state.dt,state.pressure_reference,
      state.previous_pressure_reference,double(state.step),double(state.controller_state)};
  for(auto fields:{state.fields,state.previous_fields,state.accepted_rate_fields,state.previous_rate_fields})
    for(std::size_t i=0;i<fields.size;++i) {
      const auto& f=fields.data[i];const auto v=f.values;
      for(int z=0;z<v.interior.z;++z)for(int y=0;y<v.interior.y;++y)for(int x=0;x<v.interior.x;++x)
        for(unsigned c=0;c<v.components;++c)result.push_back(v.unchecked({x,y,z},c));
    }
  for(auto flux:{state.final_mass_flux,state.previous_mass_flux})
    for(auto v:{flux.x,flux.y,flux.z})
      for(int z=0;z<v.extents.z;++z)for(int y=0;y<v.extents.y;++y)for(int x=0;x<v.extents.x;++x)
        result.push_back(v.unchecked({x,y,z}));
  return result;
}
bool run(const std::filesystem::path& assets,const std::filesystem::path& restart,int rank,bool backflow=false) {
  ValidatedModel model;
  auto status=CaseCompiler::load_and_compile(MPI_COMM_WORLD,assets,model);
  if(!all(bool(status)))return false;
  model.time.scheme=TimeScheme::cn_be;model.solver.coupling=CouplingKind::outer_corrected;
  model.pressure_reference=PressureReferenceKind::boundary_absolute;
  model.legacy_time_fingerprint=model.fingerprint;
  model.fingerprint^=0x434e455346484541ULL;
  model.boundaries[2].flow_kind=BoundaryKind::symmetry;
  model.boundaries[2].scalars={{"A",ScalarBoundaryKind::zero_gradient}};
  auto& outlet=model.boundaries[3];outlet.flow_kind=BoundaryKind::pressure_outlet;
  outlet.pressure=101325;outlet.allow_backflow=true;outlet.backflow_temperature=300;
  outlet.scalars={{"A",ScalarBoundaryKind::zero_gradient}};
  outlet.scalars[0].backflow_kind=ScalarBoundaryKind::dirichlet;outlet.scalars[0].backflow_value=backflow ? 0. : .25;
  const auto create=[&](ProductDriver& driver) {
    CompiledCasePlan plan;auto s=ProductCompiler::compile(MPI_COMM_WORLD,model,assets,plan);
    if(s && (!plan.summary().conservative_total_energy || plan.summary().midpoint_enthalpy))s={StatusCode::invalid_plan,19130};
    if(s)s=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
    return s;
  };
  ProductDriver driver;status=create(driver);
  DriverInitialState initial;initial.pressure_reference=101325;initial.temperature=300;
  if(backflow)initial.velocity.y=-.01;
  double fraction=.25;initial.transported_scalars={&fraction,1};
  if(status)status=driver.initialize(initial);
  if(!all(bool(status))) {
    if(rank==0)std::cerr<<"CN heat setup="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
    return false;
  }
  double max_mass{},max_energy{},max_restart{};
  auto advance=[&](ProductDriver& d) {
    DriverStepReport report;auto s=d.advance({1,1,1,1,1},report);
    const auto& c=report.conservation;
    const double mass=std::abs(c.mass_balance_defect)/std::max({1.,std::abs(c.mass_bdf_rate),std::abs(c.mass_outflow)});
    const double energy=std::abs(c.total_energy_balance_defect)/std::max({1.,std::abs(c.total_energy_bdf_rate),std::abs(c.enthalpy_outflow),std::abs(c.kinetic_energy_outflow)});
    max_mass=std::max(max_mass,mass);max_energy=std::max(max_energy,energy);
    bool okay=s && report.accepted && c.valid && mass<1e-6 && energy<1e-6;
    if(rank==0)std::cout<<std::setprecision(17)<<"CN_heat backflow="<<backflow<<" step="<<report.accepted_step
      <<" status="<<unsigned(s.code)<<'/'<<s.detail<<" stage="<<report.failed_stage
      <<" accepted="<<report.accepted<<" mass="<<mass<<" energy="<<energy<<" passed="<<okay<<'\n';
    return all(okay);
  };
  if(!advance(driver))return false;
  // Zero-fuel backflow exercises mirrors after reaction changes the owner.
  // The full runtime mass/energy/EOS acceptance above owns this regression.
  if(backflow)return true;
  RestartSnapshot state;status=driver.committed_restart_snapshot(state);
  if(!all(bool(status)))return false;
  double tmin=1e100,tmax=-1e100,ymin=1e100,ymax=-1e100;
  ConstFieldView h{},a{},auxiliary{};
  for(std::size_t i=0;i<state.fields.size;++i) {
    const auto& f=state.fields.data[i];
    if(f.role==RestartFieldRole::stochastic_auxiliary)auxiliary=f.values;
    if(f.role==RestartFieldRole::enthalpy)h=f.values;
    if(f.role==RestartFieldRole::independent_species)a=f.values;
  }
  if(!all(h.base && a.base && auxiliary.base && auxiliary.components==3 &&
      auxiliary.ghosts.x==0 && auxiliary.ghosts.y==0 && auxiliary.ghosts.z==0))return false;
  for(int z=0;z<h.interior.z;++z)for(int y=0;y<h.interior.y;++y)for(int x=0;x<h.interior.x;++x) {
    const double ya=a.unchecked({x,y,z},0),temperature=298.15+(h.unchecked({x,y,z},0)-1e5*ya)/1000.;
    tmin=std::min(tmin,temperature);tmax=std::max(tmax,temperature);
    ymin=std::min(ymin,ya);ymax=std::max(ymax,ya);
  }
  MPI_Allreduce(MPI_IN_PLACE,&tmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ymin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ymax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(rank==0)std::cout<<"CN_heat physical Tmin="<<tmin<<" Tmax="<<tmax<<" Ymin="<<ymin<<" Ymax="<<ymax<<'\n';
  if(!all(tmin>300 && tmax<301 && ymin>.249 && ymax<.25))return false;
  const auto saved=values(state);
  status=RestartWriter::write(MPI_COMM_WORLD,restart,state);
  ProductDriver restored;if(status)status=create(restored);
  RestartExpected expected;if(status)status=restored.restart_expected(expected);
  RestartImage image;if(status)status=RestartReader::load(MPI_COMM_WORLD,restart,expected,image);
  if(status)status=restored.initialize_restart(image);
  if(status)status=restored.committed_restart_snapshot(state);
  if(!all(bool(status) && saved==values(state)))return false;
  // A signed field0 tuple represents a distinct pressure-density state.
  // Restore it with the physical ensemble unchanged and inspect the public
  // density/temperature fields, independently of a subsequent time step.
  {
    auto split=image;
    for(auto* level:{&split.fields,&split.previous_fields})
      for(auto& field:*level)if(field.role==RestartFieldRole::stochastic_auxiliary)
        for(std::size_t i=0;i<field.values.size();i+=3) {
          field.values[i]=-.01;field.values[i+1]=1.01;
          field.values[i+2]=(field.values[i+2]+1000)*1.01;
        }
    ProductDriver probe;status=create(probe);
    if(status)status=probe.initialize_restart(split);
    ThermodynamicsPlan thermo;
    if(status)status=ThermodynamicsPlan::compile(model.thermophysics,
        {model.transported_scalars.data(),model.transported_scalars.size()},thermo);
    ConstFieldView density{},temperature{};
    if(status)status=probe.committed_thermo_for_test(StateRole::accepted_n,density,temperature);
    if(!all(bool(status))) {
      if(rank==0)std::cerr<<"split_restore setup="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
      return false;
    }
    double error{},gap{};
    std::size_t i{};
    const auto cells=density.interior;
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
      double zero=0,physical_y=image.fields[3].values[i];
      ThermoState selected,physical;
      auto q=thermo.evaluate(image.pressure_reference+image.fields[1].values[i],
          split.fields.back().values[3*i+2]/1.01,{&zero,1},{},selected,300);
      if(q)q=thermo.evaluate(image.pressure_reference+image.fields[1].values[i],
          image.fields[2].values[i],{&physical_y,1},{},physical,300);
      if(!q) {status=q;break;}
      const double expected_density=selected.rho/1.01;
      error=std::max({error,std::abs(density.unchecked({x,y,z},0)/expected_density-1),
          std::abs(temperature.unchecked({x,y,z},0)/physical.temperature-1)});
      gap=std::max(gap,std::abs(expected_density/physical.rho-1));
    }
    MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&gap,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(rank==0)std::cout<<"CN_heat split_restore error="<<error<<" density_gap="<<gap<<'\n';
    if(!all(bool(status) && error<1e-11 && gap>.01))return false;
  }
  for(unsigned step=0;step<2;++step) {
    detail::suppress_pressure_correction_warm_start_once_for_test();
    if(!advance(driver))return false;
    status=driver.committed_restart_snapshot(state);
    if(!all(bool(status)))return false;
    const auto uninterrupted=values(state);
    detail::suppress_pressure_correction_warm_start_once_for_test();
    if(!advance(restored))return false;
    status=restored.committed_restart_snapshot(state);
    if(!all(bool(status)))return false;
    const auto resumed=values(state);
    if(!all(uninterrupted.size()==resumed.size()))return false;
    for(std::size_t i=0;i<resumed.size();++i)
      max_restart=std::max(max_restart,std::abs(uninterrupted[i]-resumed[i])/
          std::max({1.,std::abs(uninterrupted[i]),std::abs(resumed[i])}));
  }
  MPI_Allreduce(MPI_IN_PLACE,&max_restart,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(rank==0)std::cout<<"CN_heat final restart_error="<<max_restart<<" mass="<<max_mass<<" energy="<<max_energy<<'\n';
  return all(max_restart<1e-10);
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  int id=int(getpid());MPI_Bcast(&id,1,MPI_INT,0,MPI_COMM_WORLD);
  const auto restart=std::filesystem::temp_directory_path()/("hf-cn-"+std::to_string(id));
  bool okay=argc==2 && run(argv[1],restart,rank);
  if(okay)okay=run(argv[1],restart,rank,true);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)std::filesystem::remove_all(restart);
  MPI_Finalize();return okay?0:1;
}
