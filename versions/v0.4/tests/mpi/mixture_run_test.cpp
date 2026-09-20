// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using namespace hundun::v04;

namespace {
std::array<double,2> composition(int x, bool variable_density=false) {
  if (variable_density) {
    const std::array<std::array<double,2>,4> profile{{
        {0,.7},{.1,.66},{.3,.51},{.4,.47}}};
    return profile[std::clamp(x-1,0,3)];
  }
  return x<4 ? std::array<double,2>{0.0,.6} :
         x==4 ? std::array<double,2>{.2,.7} : std::array<double,2>{.5,.5};
}

bool run(double formation, std::vector<double>& fractions, bool variable_density=false) {
  const Int3 cells{12,8,8};
  constexpr double dt=1e-3, pressure=100000.0, temperature=300.0, speed=1.0;
  auto model=test::product_model(cells);
  model.turbulence=TurbulenceKind::none;
  model.time.scheme=TimeScheme::cn_be;
  model.time.control=TimeControlKind::fixed;
  model.time.initial_dt=model.time.minimum_dt=model.time.maximum_dt=dt;
  model.time.maximum_retries=1;
  // Resolve outer-iteration error below the unchanged 1e-10 temperature
  // comparison when testing formation-enthalpy reference invariance.
  model.solver.cold_stopping=ColdStoppingSpec{dt,1e-12,1e-14,1e-14};
  model.legacy_time_fingerprint=model.fingerprint+1;
  model.pressure_reference=PressureReferenceKind::boundary_absolute;
  auto& inlet=model.boundaries[0];
  inlet.flow_kind=BoundaryKind::velocity_inlet;
  inlet.velocity={speed,0,0}; inlet.direction={1,0,0}; inlet.temperature=temperature;
  inlet.scalars={{"air",ScalarBoundaryKind::dirichlet,0.0},
                 {"B",ScalarBoundaryKind::dirichlet,variable_density ? .7 : .6}};
  auto& outlet=model.boundaries[1];
  outlet.flow_kind=BoundaryKind::pressure_outlet;
  outlet.pressure=pressure; outlet.backflow_temperature=temperature;
  outlet.scalars={{"air",ScalarBoundaryKind::zero_gradient},
                  {"B",ScalarBoundaryKind::zero_gradient}};
  auto gas=model.thermophysics.species.front();
  gas.transport_law=TransportLaw::perry;
  gas.viscosity_reference=gas.conductivity=0;
  gas.prandtl=.7; gas.critical_temperature=126.2; gas.critical_pressure=33.5;
  const double R=kUniversalGasConstant/gas.molecular_weight;
  model.thermophysics.species={gas,gas,gas};
  model.thermophysics.species[0].nasa7_low[5]=
      model.thermophysics.species[0].nasa7_high[5]=formation/R;
  model.thermophysics.species[1].stable_name="B";
  model.thermophysics.species[2].stable_name="balance";
  if (variable_density) {
    model.thermophysics.species[1].molecular_weight=16.0;
    model.thermophysics.species[2].molecular_weight=44.0;
  }
  const auto mixture_cp=[&](const std::array<double,2>& Y) {
    return 3.5*kUniversalGasConstant*(Y[0]/gas.molecular_weight+
        Y[1]/model.thermophysics.species[1].molecular_weight+
        (1-Y[0]-Y[1])/model.thermophysics.species[2].molecular_weight);
  };
  model.transported_scalars={{"air",TransportedScalarRole::species,.7,.7},
                             {"B",TransportedScalarRole::species,.7,.7}};
  CompiledCasePlan plan;
  ThermodynamicsPlan thermo;
  auto status=ThermodynamicsPlan::compile(model.thermophysics,
      {model.transported_scalars.data(),model.transported_scalars.size()},thermo);
  for (int x : {0,4,5}) {
    const auto Y=composition(x,variable_density);
    ThermoState state;
    if (status) status=thermo.evaluate(pressure,mixture_cp(Y)*temperature+formation*Y[0],
        {Y.data(),Y.size()},{speed,0,0},state,temperature);
    if (!status) std::cerr << "mixture input x=" << x << " status="
        << unsigned(status.code) << '/' << status.detail << '\n';
  }
  if (status) status=ProductCompiler::compile(MPI_COMM_WORLD,model,{},plan);
  ProductDriver driver;
  if (status) status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;
  if (status) status=driver.restart_expected(expected);
  RestartImage start;
  start.global_cells=expected.global_cells; start.patch=expected.target_patch;
  start.plan=expected.plan; start.schema=expected.schema; start.geometry=expected.geometry;
  start.time=0.0; start.dt=dt; start.step=4; start.controller_state=1;
  start.pressure_reference=pressure; start.backward_euler_recovery=true;
  const auto n=start.patch.cells;
  std::size_t species{};
  for (std::size_t f=0; f<expected.fields.size && status; ++f) {
    const auto& field=expected.fields.data[f];
    RestartImageField out;
    out.role=field.role; out.field=field.field; out.components=field.components;
    out.values.resize(std::size_t(n.x)*n.y*n.z*field.components);
    for (int z=0; z<n.z; ++z) for (int y=0; y<n.y; ++y) for (int x=0; x<n.x; ++x) {
      const auto i=std::size_t((z*n.y+y)*n.x+x)*field.components;
      const auto Y=composition(x+start.patch.begin.x,variable_density);
      if (field.role==RestartFieldRole::velocity) out.values[i]=speed;
      else if (field.role==RestartFieldRole::pressure_absolute) out.values[i]=pressure;
      else if (field.role==RestartFieldRole::enthalpy)
        out.values[i]=mixture_cp(Y)*temperature+formation*Y[0];
      else if (field.role==RestartFieldRole::independent_species) out.values[i]=Y[species];
    }
    if (field.role==RestartFieldRole::independent_species) ++species;
    start.fields.push_back(std::move(out));
  }
  for (unsigned a=0; a<3; ++a) {
    auto extent=n; (a==0 ? extent.x : a==1 ? extent.y : extent.z)++;
    start.final_mass_flux[a].assign(std::size_t(extent.x)*extent.y*extent.z,
        a==0 ? pressure/(R*temperature)*speed/(8.0*16.0) : 0.0);
    if (a==0 && variable_density)
      for (int z=0;z<extent.z;++z) for (int y=0;y<extent.y;++y)
        for (int x=0;x<extent.x;++x) {
          const auto left=composition(x+start.patch.begin.x-1,true);
          const auto right=composition(x+start.patch.begin.x,true);
          const double face_density=.5*pressure/temperature*3.5*
              (1/mixture_cp(left)+1/mixture_cp(right));
          start.final_mass_flux[a][std::size_t((z*extent.y+y)*extent.x+x)]=
              face_density*speed/(8.0*16.0);
        }
  }
  if (status) status=driver.initialize_restart(start);
  if (!status) std::cerr << "mixture initialization=" << unsigned(status.code)
      << '/' << status.detail << " formation=" << formation << '\n';
  DriverStepReport report;
  for (unsigned step=0; step<3 && status; ++step)
    status=driver.advance({dt/model.time.convective_cfl,1,1,1,1},report);
  RestartSnapshot snapshot;
  if (status) status=driver.committed_restart_snapshot(snapshot);
  std::array<ConstFieldView,2> Y;
  ConstFieldView h,p;
  species=0;
  if (status) for (std::size_t f=0; f<snapshot.fields.size; ++f) {
    const auto& field=snapshot.fields.data[f];
    if (field.role==RestartFieldRole::independent_species && species<2)
      Y[species++]=field.values;
    if (field.role==RestartFieldRole::enthalpy) h=field.values;
    if (field.role==RestartFieldRole::pressure_perturbation) p=field.values;
  }
  bool passed=status && report.accepted && species==2 && h.base && p.base &&
      std::abs(snapshot.time-3*dt)<1e-15;
  double error[3]{};
  if (passed) for (int z=0; z<n.z; ++z) for (int y=0; y<n.y; ++y) for (int x=0; x<n.x; ++x) {
    const Int3 c{x,y,z};
    const double a=Y[0].unchecked(c,0), b=Y[1].unchecked(c,0), dependent=1-a-b;
    passed &= a>=0 && b>=0 && dependent>=-1e-13;
    fractions.push_back(a); fractions.push_back(b);
    const double local_temperature=(h.unchecked(c,0)-formation*a)/mixture_cp({a,b});
    error[0]=std::max(error[0],std::abs(local_temperature-temperature));
    error[1]=std::max(error[1],std::abs(snapshot.pressure_reference+p.unchecked(c,0)-pressure)/pressure);
    if (variable_density) {
      passed &= std::isfinite(local_temperature) && local_temperature>0;
      fractions.push_back(local_temperature/temperature);
      fractions.push_back((snapshot.pressure_reference+p.unchecked(c,0))/pressure);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,error,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if (!variable_density) passed &= error[0]<1e-8 && error[1]<1e-11;
  int all=passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  int rank{}; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank==0) std::cout << "mixture_run formation=" << formation
      << " variable_density=" << variable_density
      << " status=" << unsigned(status.code) << '/' << status.detail
      << " temperature_error=" << error[0] << " pressure_relative=" << error[1]
      << " passed=" << all << '\n';
  return all;
}
} // namespace

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  std::vector<double> ordinary;
  bool passed=run(0.0,ordinary);
  double error=0;
  for (double formation : {-1.2e7,1.2e7}) {
    std::vector<double> shifted;
    passed=run(formation,shifted) && passed;
    if (ordinary.size()!=shifted.size()) passed=false;
    else for (std::size_t i=0; i<ordinary.size(); ++i)
      error=std::max(error,std::abs(ordinary[i]-shifted[i]));
  }
  MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  passed &= error<1e-11;
  int rank{}; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (rank==0) std::cout << "mixture_reference composition_error=" << error << '\n';
  std::vector<double> variable_reference;
  passed=run(0,variable_reference,true) && passed;
  error=0;
  std::array<double,4> component_error{};
  for (double formation : {-1.2e7,1.2e7}) {
    std::vector<double> shifted;
    passed=run(formation,shifted,true) && passed;
    if (shifted.size()!=variable_reference.size()) passed=false;
    else for (std::size_t i=0;i<shifted.size();++i) {
      error=std::max(error,std::abs(shifted[i]-variable_reference[i]));
      component_error[i%4]=std::max(component_error[i%4],
          std::abs(shifted[i]-variable_reference[i]));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,component_error.data(),4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  passed &= error<1e-10;
  if (rank==0) std::cout << "mixture_variable_reference normalized_state_error=" << error
      << " Y0="<<component_error[0]<<" Y1="<<component_error[1]
      << " temperature="<<component_error[2]<<" pressure="<<component_error[3]<<'\n';
  MPI_Finalize();
  return passed ? 0 : 1;
}
