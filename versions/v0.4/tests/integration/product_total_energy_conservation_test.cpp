// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
// Independent public-snapshot total-energy oracle for unforced periodic shear.
#include "hundun/v04_app.hpp"
#include "../support/product_fixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
using namespace hundun::v04;
namespace {
constexpr Int3 cells{8,4,4};
constexpr double pressure=101325.0;
constexpr double r_a=kUniversalGasConstant/28.96546;
int rank=0;
double width(int axis,int) {
  return axis==0 ? 2.0/cells.x : axis==1 ? 1.0/cells.y : 0.5/cells.z;
}
ValidatedModel model(bool,double dt) {
  auto m=test::product_model(cells);
  m.turbulence=TurbulenceKind::none;
  m.time.control=TimeControlKind::fixed;
  m.time.initial_dt=m.time.minimum_dt=m.time.maximum_dt=dt;
  m.time.maximum_growth=1.0; m.time.maximum_retries=1U;
  m.solver.pressure.absolute_tolerance=m.solver.pressure.relative_tolerance=1e-13;
  m.solver.pressure.maximum_iterations=800U;
  m.solver.pressure.true_residual_interval=4U;
  m.solver.pressure.krylov_restart=64U;
  m.solver.terminal.eos=m.solver.terminal.continuity=1e-12;
  m.solver.terminal.closed_mass=m.solver.terminal.gauge=1e-12;
  m.thermophysics.species[0].viscosity_reference=1e-5;
  m.thermophysics.species[0].conductivity=1e-2;
  m.transported_scalars={{"constant",TransportedScalarRole::passive_scalar,1.0,1.0},
                        {"wave",TransportedScalarRole::passive_scalar,1.0,1.0}};
  return m;
}
RestartImage image(const RestartExpected& expected,bool,bool,double dt) {
  RestartImage r;
  r.global_cells=expected.global_cells; r.patch=expected.target_patch;
  r.plan=expected.plan; r.schema=expected.schema; r.geometry=expected.geometry;
  r.time=0.0; r.dt=dt; r.step=4U; r.controller_state=1U;
  r.pressure_reference=pressure; r.backward_euler_recovery=true;
  const auto n=r.patch.cells;
  for(std::size_t i=0;i<expected.fields.size;++i) {
    const auto& d=expected.fields.data[i];
    RestartImageField f; f.role=d.role; f.field=d.field; f.components=d.components;
    f.values.resize(static_cast<std::size_t>(n.x)*n.y*n.z*d.components);
    r.fields.push_back(std::move(f));
  }
  for(int a=0;a<3;++a) {
    auto ext=n; (a==0 ? ext.x : a==1 ? ext.y : ext.z)++;
    r.final_mass_flux[a].resize(static_cast<std::size_t>(ext.x)*ext.y*ext.z);
  }
  return r;
}
} // namespace

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const int outcome=[&]() {
    const double dt=argc>1 ? std::stod(argv[1]) : 0.001;
    auto m=model(false,dt);
    m.solver.coupling=CouplingKind::piso;
    m.schemes.momentum=ConvectionScheme::tvd2;
    m.schemes.enthalpy=ConvectionScheme::tvd2;
    CompiledCasePlan plan;
    Status status=ProductCompiler::compile(MPI_COMM_WORLD,m,{},plan);
    ProductDriver driver;
    if(status) status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
    RestartExpected expected;
    if(status) status=driver.restart_expected(expected);
    auto start=image(expected,false,true,dt);
    const auto n=start.patch.cells;
    const double temp=295.0,cp=3.5*r_a,rho=pressure/(r_a*temp),u=2.0;
    const auto v=[&](int x){return 10.0*std::sin(2.0*std::acos(-1.0)*(x+0.5)/cells.x);};
    for(auto& field:start.fields)
      for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x) {
        const auto i=static_cast<std::size_t>((z*n.y+y)*n.x+x)*field.components;
        if(field.role==RestartFieldRole::velocity) {
          field.values[i]=u;field.values[i+1]=v(x+start.patch.begin.x);field.values[i+2]=0.0;
        } else if(field.role==RestartFieldRole::enthalpy) field.values[i]=cp*temp;
        else if(field.role==RestartFieldRole::pressure_perturbation) field.values[i]=0.0;
        else if(field.role==RestartFieldRole::pressure_absolute) field.values[i]=pressure;
        else if(field.role==RestartFieldRole::transported_scalar) field.values[i]=0.2;
      }
    for(int axis=0;axis<3;++axis) {
      auto ext=n;(axis==0 ? ext.x : axis==1 ? ext.y : ext.z)++;
      std::size_t i=0;
      for(int z=0;z<ext.z;++z) for(int y=0;y<ext.y;++y) for(int x=0;x<ext.x;++x) {
        const int gx=x+start.patch.begin.x,gy=y+start.patch.begin.y,gz=z+start.patch.begin.z;
        start.final_mass_flux[axis][i++]=axis==0 ? rho*u*width(1,gy)*width(2,gz) :
          axis==1 ? rho*v(gx)*width(0,gx)*width(2,gz) : 0.0;
      }
    }
    if(status) status=driver.initialize_restart(start);
    const auto inventory=[&](long double& total)->Status {
      RestartSnapshot snapshot;
      auto s=driver.committed_restart_snapshot(snapshot);
      ConstFieldView h{},p{},vel{};
      if(s) for(std::size_t f=0;f<snapshot.fields.size;++f) {
        const auto& field=snapshot.fields.data[f];
        if(field.role==RestartFieldRole::enthalpy) h=field.values;
        if(field.role==RestartFieldRole::pressure_perturbation) p=field.values;
        if(field.role==RestartFieldRole::velocity) vel=field.values;
      }
      long double local=0.0;
      if(s) for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x) {
        const Int3 cell{x,y,z};
        const long double enthalpy=h.unchecked(cell,0U);
        const long double absolute=snapshot.pressure_reference+p.unchecked(cell,0U);
        const long double density=absolute/(r_a*(enthalpy/cp));
        long double kinetic=0.0;
        for(int c=0;c<3;++c) {const long double speed=vel.unchecked(cell,c);kinetic+=0.5L*speed*speed;}
        local+=width(0,x+start.patch.begin.x)*width(1,y+start.patch.begin.y)*width(2,z+start.patch.begin.z)*
          (density*(enthalpy+kinetic)-absolute);
      }
      MPI_Allreduce(&local,&total,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      return s;
    };
    long double initial{},before{},previous{},after{};
    if(status) status=inventory(initial);
    before=previous=initial;
    bool passed=static_cast<bool>(status);
    // One explicit BE history recovery followed by genuine BDF2 steps.
    // Both the BDF residual and the unchanged closed-system inventory are
    // independent oracles; neither uses the driver's reported defect.
    for (int step=0;step<3 && passed;++step) {
      DriverStepReport report;
      status=driver.advance({dt,dt,dt,dt,dt},report);
      if(status) status=inventory(after);
      const double defect=static_cast<double>(step==0 ? (after-before)/dt :
          (1.5L*after-2.0L*before+0.5L*previous)/dt);
      if(rank==0) std::cout<<std::setprecision(17)<<"PERIODIC_ENERGY step="<<step
        <<" dt="<<dt<<" status="<<unsigned(status.code)<<'/'<<status.detail
        <<" accepted="<<report.accepted<<" independent_defect_W="<<defect
        <<" report_defect_W="<<report.conservation.total_energy_balance_defect
        <<" inventory_change_J="<<(after-initial)<<'\n';
      passed=status && report.accepted && std::abs(defect)<1e-4 &&
          std::abs(after-initial)<1e-6L;
      previous=before; before=after;
    }
    return passed ? 0 : 1;
  }();
  MPI_Finalize();
  return outcome;
}
