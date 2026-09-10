// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

// A large, nonreacting formation reference exposes slow composition feedback
// through the real p/h/EOS solve. The original Picard map exhausts 16 sweeps.
// This is a convergence/inventory regression, NOT a contact-preservation claim.
#include "hundun/v04_app.hpp"
#include "../support/product_fixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
using namespace hundun::v04;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank{};
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int result = [&]() {
    const bool wave=argc==2 && std::string(argv[1])=="--wave";
    const Int3 cells=wave ? Int3{16,4,4} : Int3{8,4,4};
    const double dt=wave ? 0.0001 : 0.001;
    constexpr double pressure=101325.0, temperature=295.0;
    constexpr double formation=-1.2e7, speed=0.5;
    auto model=test::product_model(cells);
    model.turbulence=TurbulenceKind::none;
    model.solver.coupling=CouplingKind::piso;
    model.schemes.momentum=ConvectionScheme::central2;
    model.schemes.enthalpy=model.schemes.species=ConvectionScheme::tvd2;
    model.time.control=TimeControlKind::fixed;
    model.time.initial_dt=model.time.minimum_dt=model.time.maximum_dt=dt;
    model.time.maximum_growth=1.0;
    model.time.maximum_retries=1U;
    model.solver.pressure.absolute_tolerance=model.solver.pressure.relative_tolerance=1e-13;
    model.solver.pressure.maximum_iterations=800U;
    model.solver.pressure.true_residual_interval=4U;
    model.solver.pressure.krylov_restart=64U;
    model.solver.terminal.eos=model.solver.terminal.continuity=1e-12;
    model.solver.terminal.closed_mass=model.solver.terminal.gauge=1e-12;
    auto& a=model.thermophysics.species.front();
    a.stable_name="A"; a.viscosity_reference=1e-5; a.conductivity=1e-2;
    const double gas=kUniversalGasConstant/a.molecular_weight;
    auto b=a; b.stable_name="B";
    a.nasa7_low[5]=a.nasa7_high[5]=formation/gas;
    model.thermophysics.species.push_back(b);
    model.transported_scalars={{"constant",TransportedScalarRole::passive_scalar,1.0,1.0},
                              {"A",TransportedScalarRole::species,1.0,1.0}};
    ThermodynamicsPlan thermo;
    Status status=ThermodynamicsPlan::compile(model.thermophysics,
        {model.transported_scalars.data(),model.transported_scalars.size()},thermo);
    CompiledCasePlan plan;
    if(status) status=ProductCompiler::compile(MPI_COMM_WORLD,model,{},plan);
    ProductDriver driver;
    if(status) status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
    RestartExpected expected;
    if(status) status=driver.restart_expected(expected);
    RestartImage start;
    start.global_cells=expected.global_cells; start.patch=expected.target_patch;
    start.plan=expected.plan; start.schema=expected.schema; start.geometry=expected.geometry;
    start.time=0.0; start.dt=dt; start.step=4U; start.controller_state=1U;
    start.pressure_reference=pressure; start.backward_euler_recovery=true;
    const auto n=start.patch.cells;
    const double rho=pressure/(gas*temperature);
    for(std::size_t f=0;f<expected.fields.size && status;++f) {
      const auto& field=expected.fields.data[f];
      RestartImageField out;
      out.role=field.role; out.field=field.field; out.components=field.components;
      out.values.resize(static_cast<std::size_t>(n.x)*n.y*n.z*field.components);
      for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x) {
        const auto i=static_cast<std::size_t>((z*n.y+y)*n.x+x)*field.components;
        const double q=0.2+0.05*std::sin(2.0*std::acos(-1.0)*
            (x+start.patch.begin.x+0.5)/cells.x);
        const double pi=wave ? pressure*0.03*std::cos(2.0*std::acos(-1.0)*
            (x+start.patch.begin.x+0.5)/cells.x) : 0.0;
        if(field.role==RestartFieldRole::velocity) out.values[i]=speed;
        else if(field.role==RestartFieldRole::pressure_absolute) out.values[i]=pressure+pi;
        else if(field.role==RestartFieldRole::pressure_perturbation) out.values[i]=pi;
        else if(field.role==RestartFieldRole::enthalpy) {
          double cp{}, r{};
          status=thermo.mixture_enthalpy(temperature,{&q,1U},out.values[i],cp,r);
        } else if(field.role==RestartFieldRole::independent_species) out.values[i]=q;
        else if(field.role==RestartFieldRole::transported_scalar) out.values[i]=0.2;
      }
      start.fields.push_back(std::move(out));
    }
    for(int axis=0;axis<3;++axis) {
      auto ext=n; (axis==0 ? ext.x : axis==1 ? ext.y : ext.z)++;
      start.final_mass_flux[axis].assign(static_cast<std::size_t>(ext.x)*ext.y*ext.z,
          axis==0 ? rho*speed*(1.0/cells.y)*(0.5/cells.z) : 0.0);
    }
    if(status) status=driver.initialize_restart(start);
    const auto inventory=[&](long double& mass,long double& energy) {
      RestartSnapshot snapshot;
      auto s=driver.committed_restart_snapshot(snapshot);
      ConstFieldView h{},p{},u{},q{};
      if(s) for(std::size_t f=0;f<snapshot.fields.size;++f) {
        const auto& field=snapshot.fields.data[f];
        if(field.role==RestartFieldRole::enthalpy) h=field.values;
        if(field.role==RestartFieldRole::pressure_perturbation) p=field.values;
        if(field.role==RestartFieldRole::velocity) u=field.values;
        if(field.role==RestartFieldRole::independent_species) q=field.values;
      }
      long double local[2]{},global[2]{};
      if(s) for(int z=0;z<n.z;++z) for(int y=0;y<n.y;++y) for(int x=0;x<n.x;++x) {
        const Int3 c{x,y,z};
        const double fraction=q.unchecked(c,0U);
        ThermoState material;
        const auto evaluated=thermo.evaluate(snapshot.pressure_reference+p.unchecked(c,0U),
            h.unchecked(c,0U),{&fraction,1U},{},material,temperature);
        if(!evaluated) {s=evaluated; continue;}
        long double kinetic=0.0;
        for(int a=0;a<3;++a) {const long double v=u.unchecked(c,a); kinetic+=0.5L*v*v;}
        const long double volume=(2.0/cells.x)*(1.0/cells.y)*(0.5/cells.z);
        local[0]+=volume*material.rho*fraction;
        local[1]+=volume*(material.rho*(h.unchecked(c,0U)+kinetic)-
            (snapshot.pressure_reference+p.unchecked(c,0U)));
      }
      MPI_Allreduce(local,global,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      mass=global[0]; energy=global[1];
      const int local_ok=s ? 1 : 0;
      int all_ok{};
      MPI_Allreduce(&local_ok,&all_ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
      return all_ok ? Status{} : Status{StatusCode::numerical_failure,1U};
    };
    long double before_mass{},before_energy{},after_mass{},after_energy{};
    if(status) status=inventory(before_mass,before_energy);
    DriverStepReport report;
    if(status) status=driver.advance({dt,dt,dt,dt,dt},report);
    if(status) status=inventory(after_mass,after_energy);
    const long double mass_error=std::abs(after_mass-before_mass)/before_mass;
    const long double energy_error=std::abs(after_energy-before_energy)/std::abs(before_energy);
    // A 3% pressure wave retains C2 refinement after composition recoupling.
    // The 1% wave now converges with only the initial C1 and C2 directions.
    // C1's merit belongs to another momentum predictor and cannot seed it.
    bool first_refinement_observed=false, c2_history_local=true;
    const auto& path=report.pressure_energy_globalization;
    for(std::size_t i=0U;i<path.trajectory_count;++i) {
      const auto& item=path.trajectory[i];
      if(item.corrector==2U && item.refinement_iteration==1U) {
        first_refinement_observed=true;
        c2_history_local &= !item.extrapolation.attempted;
        if(rank==0) std::cout<<"C2_FIRST_REFINEMENT extrapolation="
          <<item.extrapolation.attempted<<" alpha="<<item.selected.alpha<<'\n';
      }
    }
    if(rank==0) std::cout<<std::setprecision(17)<<"FORMATION_COUPLING status="
      <<unsigned(status.code)<<'/'<<status.detail<<" accepted="<<report.accepted
      <<" trajectory="<<unsigned(path.trajectory_count)
      <<" sweeps="<<report.scalar_transport.coupling_sweeps
      <<" residual="<<report.scalar_transport.final_species_residual
      <<" mass_relative_error="<<mass_error<<" energy_relative_error="<<energy_error<<'\n';
    const bool passed=status && report.accepted && mass_error<1e-12L && energy_error<1e-12L &&
        (!wave || (first_refinement_observed && c2_history_local));
    int local_pass=passed ? 1 : 0, all_pass{};
    MPI_Allreduce(&local_pass,&all_pass,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    return all_pass ? 0 : 1;
  }();
  MPI_Finalize();
  return result;
}
