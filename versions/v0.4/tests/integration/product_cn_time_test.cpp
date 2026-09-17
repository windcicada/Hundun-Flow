// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_app.hpp"
#include "../support/product_fixture.hpp"
#include <mpi.h>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>

using namespace hundun::v04;
namespace {
constexpr Int3 cells{8,4,4};
constexpr double amplitude=1e-3, temperature=300, pressure=101325;
constexpr double viscosity=1, duration=.02;
constexpr double gas_constant=kUniversalGasConstant/28.96546;
constexpr double cp=3.5*gas_constant;
double wave(int x) {return std::sin(2*std::acos(-1.)*(x+.5)/cells.x);}

bool run(double dt, bool be, bool thermal, double& error) {
  const auto scheme=be ? TimeScheme::backward_euler : TimeScheme::cn_be;
  const auto coupling=be ? CouplingKind::piso : CouplingKind::outer_corrected;
  auto model=test::product_model(cells);
  model.fingerprint=UINT64_C(0x434e534845415231);
  model.legacy_time_fingerprint=model.fingerprint+1;
  model.turbulence=TurbulenceKind::none;
  model.pressure_reference=PressureReferenceKind::boundary_absolute;
  model.time.scheme=scheme;
  model.time.control=TimeControlKind::fixed;
  model.time.initial_dt=model.time.minimum_dt=model.time.maximum_dt=dt;
  model.time.maximum_retries=1;
  model.solver.coupling=coupling;
  model.solver.cold_stopping=ColdStoppingSpec{1.,1e-7,1e-10,1e-10};
  if(be) model.solver.cold_stopping.reset();
  model.solver.pressure.relative_tolerance=1e-13;
  model.solver.pressure.absolute_tolerance=1e-13;
  model.schemes.momentum=ConvectionScheme::central2;
  model.schemes.enthalpy=ConvectionScheme::central2;
  model.thermophysics.species[0].viscosity_reference=viscosity;
  const double conductivity=pressure/(gas_constant*temperature)*cp;
  if(thermal) {
    model.thermophysics.fixed_pressure_pa=pressure;
    model.thermophysics.species[0].conductivity=conductivity;
    model.thermophysics.species[0].viscosity_reference=1.8e-5;
    model.solver.cold_stopping=ColdStoppingSpec{1.,1e-10,1e-13,1e-10};
    model.solver.pressure.absolute_tolerance=1e-17;
    model.solver.pressure.relative_tolerance=1e-14;
  }
  // The exact shear has zero outlet-normal velocity. Keep the transverse
  // extrapolation condition continuous through FP64 normal-velocity noise.
  for(unsigned face: {4U,5U}) {
    auto& b=model.boundaries[face];
    b.flow_kind=BoundaryKind::pressure_outlet;
    b.pressure=pressure; b.backflow_temperature=temperature; b.allow_backflow=false;
  }
  CompiledCasePlan plan;
  auto status=ProductCompiler::compile(MPI_COMM_SELF,model,{},plan);
  if(status && (plan.summary().time_scheme!=scheme ||
                plan.summary().coupling!=coupling)) return false;
  if(!status) std::cerr<<"CN shear compile failure\n";
  ProductDriver driver;
  if(status) status=ProductDriver::create(MPI_COMM_SELF,std::move(plan),driver);
  if(!status) std::cerr<<"CN shear driver setup failure\n";
  RestartExpected expected;
  if(status) status=driver.restart_expected(expected);
  RestartImage seed;
  if(status) {
    seed.global_cells=expected.global_cells; seed.patch=expected.target_patch;
    seed.plan=expected.plan; seed.schema=expected.schema; seed.geometry=expected.geometry;
    seed.dt=dt; seed.step=1; seed.controller_state=1;
    seed.pressure_reference=pressure; seed.backward_euler_recovery=true;
    for(std::size_t f=0;f<expected.fields.size;++f) {
      const auto d=expected.fields.data[f]; RestartImageField out;
      out.role=d.role; out.field=d.field; out.components=d.components;
      out.values.resize(std::size_t(cells.x)*cells.y*cells.z*d.components);
      for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
        const auto i=(std::size_t(z)*cells.y+y)*cells.x+x;
        if(d.role==RestartFieldRole::enthalpy) out.values[i]=cp*(temperature+(thermal ? amplitude*wave(x) : 0.));
        if(d.role==RestartFieldRole::pressure_absolute) out.values[i]=pressure;
        if(d.role==RestartFieldRole::velocity && !thermal) out.values[3*i+1]=amplitude*wave(x);
      }
      seed.fields.push_back(std::move(out));
    }
    const double rho=pressure/(gas_constant*temperature);
    seed.final_mass_flux[0].assign(std::size_t(cells.x+1)*cells.y*cells.z,0);
    seed.final_mass_flux[1].resize(std::size_t(cells.x)*(cells.y+1)*cells.z);
    seed.final_mass_flux[2].assign(std::size_t(cells.x)*cells.y*(cells.z+1),0);
    const double area=(2./cells.x)*(.5/cells.z);
    for(int z=0;z<cells.z;++z) for(int y=0;y<=cells.y;++y) for(int x=0;x<cells.x;++x)
      seed.final_mass_flux[1][(std::size_t(z)*(cells.y+1)+y)*cells.x+x]=thermal ? 0. : rho*amplitude*wave(x)*area;
    status=driver.initialize_restart(seed);
    if(!status) std::cerr<<"CN shear restart failure\n";
  }
  const int steps=std::lround(duration/dt);
  DriverStepReport report;
  for(int i=0;status && i<steps;++i) {
    status=driver.advance({dt,dt,dt,dt,dt},report);
    if(status && (!report.accepted || report.proposal.dt!=dt || report.effective_bdf.order!=1 ||
                  !report.conservation.valid || !report.terminal_equations.valid)) return false;
    if(status) {
      const double energy=report.terminal_equations.internal_energy+report.terminal_equations.kinetic_energy;
      if(!(energy>0) || std::abs(report.conservation.total_energy_balance_defect)*dt>1e-10*energy ||
         std::abs(report.conservation.cumulative_energy_defect)>1e-10*energy) return false;
    }
  }
  RestartSnapshot snapshot;
  if(status) status=driver.committed_restart_snapshot(snapshot);
  if(!status) {
    std::cerr<<(thermal ? "CN heat dt=" : "CN shear dt=")<<dt<<" status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
    return false;
  }
  // The periodic central Laplacian has this exact Fourier eigenvalue.
  // This semi-discrete reference separates time error from spatial error.
  // O(amplitude^2) thermodynamic feedback is resolved by the energy equation.
  const double dx=2./cells.x, rho=pressure/(gas_constant*temperature);
  const double lambda=4*(thermal ? conductivity/(rho*cp) : viscosity/rho)*std::pow(std::sin(std::acos(-1.)/cells.x),2)/(dx*dx);
  const double exact_amplitude=amplitude*std::exp(-lambda*duration);
  // COAST cmod halves each spatial row and adds its accepted-state action
  // before step adds rho/dt. For this eigenmode that gives the CN rational
  // factor below; the BE control uses the full implicit spatial row.
  const double factor=be ? 1/(1+lambda*dt) : (1-.5*lambda*dt)/(1+.5*lambda*dt);
  const double discrete_amplitude=amplitude*std::pow(factor,steps);
  bool found=false; long double square{}, projection{}, norm{}; double method_error{};
  for(std::size_t f=0;f<snapshot.fields.size;++f) {
    const auto& field=snapshot.fields.data[f];
    if(field.role!=(thermal ? RestartFieldRole::enthalpy : RestartFieldRole::velocity)) continue;
    found=true;
    for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
      const double value=thermal ? field.values.unchecked({x,y,z},0)/cp-temperature
                                 : field.values.unchecked({x,y,z},1);
      projection+=value*wave(x);norm+=wave(x)*wave(x);
      const double residual=value-exact_amplitude*wave(x);
      square+=residual*residual;
      method_error=std::max(method_error,std::abs(value-discrete_amplitude*wave(x)));
    }
  }
  error=std::sqrt(double(square)/(cells.x*cells.y*cells.z));
  if(thermal) {
    // Project the fundamental mode: quadratic expansion/advection feeds the
    // second harmonic, while this amplitude isolates the linear time scheme.
    const double resolved=static_cast<double>(projection/norm);
    error=std::abs(resolved-exact_amplitude);
    method_error=std::abs(resolved-discrete_amplitude);
  }
  std::cout<<std::setprecision(17)<<(be?"BE":"CN")<<(thermal ? " heat dt=" : " shear dt=")<<dt<<" steps="<<steps<<" lambda="<<lambda
           <<" error="<<error<<" discrete_method_error="<<method_error<<" energy_defect="<<report.conservation.cumulative_energy_defect<<'\n';
  return found && std::isfinite(error) && error>1e-13 &&
      method_error<=(thermal ? 3e-7 : 1e-11)*amplitude && std::abs(snapshot.time-duration)<1e-14;
}
}
int main(int argc,char** argv) {
  if(MPI_Init(&argc,&argv)!=MPI_SUCCESS) return 2;
  const bool be=argc==2 && std::string_view(argv[1])=="--be";
  const bool thermal=argc==2 && std::string_view(argv[1])=="--heat";
  bool passed=true; std::array<double,3> errors{};
  for(unsigned i=0;i<3;++i) passed=run(.005/std::pow(2.,i),be,thermal,errors[i]) && passed;
  if(passed) {
    const double first=std::log2(errors[0]/errors[1]),second=std::log2(errors[1]/errors[2]);
    std::cout<<(be?"BE":"CN")<<(thermal ? " enthalpy time order=" : " momentum time order=")<<first<<','<<second<<'\n';
    passed=be ? first>=.9 && second>=.9 && first<1.2 && second<1.2
              : first>=1.8 && second>=1.8;
  }
  MPI_Finalize(); return passed?0:1;
}
