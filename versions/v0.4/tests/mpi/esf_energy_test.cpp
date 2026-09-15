// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
// Isothermal affine-mixture transport exercises the complete CN/BE driver.
// A species formation-reference change must preserve every field temperature,
// the physical composition, the independent auxiliary state and acceptance.
// The inlet/outlet carry positive normal flow to isolate the energy operator
// from the separate zero-normal-velocity backflow switching regression.
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string_view>

using namespace hundun::v04;
namespace {
double formation=0;
bool pressure_flux=false,trace_species=false;
double species_b(double a,int y) {return trace_species ? (y>=2 ? 1e-24 : 0.) : .7-a;}
double normal_speed(double y) {return .1+(pressure_flux ? .02*std::sin(std::acos(-1.)*y/4.) : 0.);}
constexpr int nx=12;
constexpr double dt=.02, speed=.2, mu=.03, pressure=1e5, temperature=300;
constexpr double mw=28.96546, cp=3.5*kUniversalGasConstant/mw;
constexpr double rho=pressure*mw/(kUniversalGasConstant*temperature);
constexpr double gamma=mu/.7;
double cosine(int x) {return std::cos(2*std::acos(-1.)*(x+.5)/nx);}
double offset(unsigned f,unsigned fields) {return (2*int(f)-int(fields)+1)*.02;}
double mean_a(int x) {return .3+.03*cosine(x);}
double field_a(int x,unsigned f,unsigned fields) {return mean_a(x)+offset(f,fields)*cosine(x);}

class Gas final : public portable::GasQueryProvider, public portable::GasAdvanceProvider {
 public:
  portable::GasIdentity identity;
  combustion::ChemistryIdentity closure;
  unsigned calls{};
  double tuple_error{};
  Gas() {
    identity.mechanism_sha256=std::string(64,'b');identity.phase="inert-mixture";
    identity.species_names={"A","B","C"};identity.element_names={"X"};
    identity.element_counts={1,1,1};identity.molecular_weights_kg_per_kmol={mw,mw,mw};
    identity.enthalpy_reference="nasa7-total-enthalpy";identity.composition_fingerprint=4101;
    closure.element_count=1;closure.species={{mw,{1},0},{mw,{1},0},{mw,{1},0}};
    closure.species[0].formation_enthalpy_j_per_kg=formation;
    closure.fingerprint=combustion::chemistry_identity_fingerprint(closure);
    identity.closure_fingerprint=closure.fingerprint;
  }
  const portable::GasIdentity& gas_identity() const noexcept override {return identity;}
  portable::Status query_gas(const portable::GasQuery& q,portable::GasQueryOutput& out) noexcept override {
    if(q.species_count!=3 || out.capacity!=3)return portable::Status::invalid_input;
    const double t=q.coordinates==portable::GasStateCoordinates::pressure_enthalpy ? (q.enthalpy_j_per_kg-formation*q.mass_fractions[0])/cp : q.temperature_k;
    out.sample={q.revision,q.composition_fingerprint,q.pressure_pa,t,
        q.pressure_pa*mw/(kUniversalGasConstant*t),cp*t+formation*q.mass_fractions[0],cp,mu,gamma*cp};
    for(unsigned s=0;s<3;++s) {
      out.diffusivities_m2_per_s[s]=gamma/rho;out.species_enthalpies_j_per_kg[s]=cp*t+(s==0 ? formation : 0);
      out.net_mass_rates_kg_per_m3_s[s]=0;
    }
    return portable::Status::success;
  }
  portable::Status advance_gas(const portable::GasAdvanceQuery& q,portable::GasAdvanceOutput& out) noexcept override {
    ++calls;
    tuple_error=std::max(tuple_error,std::abs(
        (q.state.enthalpy_j_per_kg-formation*q.state.mass_fractions[0])/cp-temperature));
    double d[3],h[3],w[3];portable::GasQueryOutput sample{{},d,h,w,3};
    auto status=query_gas(q.state,sample);if(status!=portable::Status::success)return status;
    for(unsigned s=0;s<3;++s) {
      out.final_mass_fractions[s]=q.state.mass_fractions[s];
      out.integrated_species_density_delta_kg_per_m3[s]=0;
    }
    out.final_sample=sample.sample;out.completed_duration_s=q.duration_s;out.internal_step_count=1;
    out.integrated_heat_release_j_per_m3=0;return portable::Status::success;
  }
};

bool run(unsigned fields,std::vector<double>& composition,std::vector<double>& thermal) {
  Gas gas;
  auto model=test::product_model({nx,4,4});
  {
    model.time.scheme=TimeScheme::cn_be;model.solver.coupling=CouplingKind::outer_corrected;
    model.pressure_reference=PressureReferenceKind::boundary_absolute;
    model.legacy_time_fingerprint=model.fingerprint+1;
  }
  model.turbulence=TurbulenceKind::none;
  {
    auto& outlet=model.boundaries[3];outlet.flow_kind=BoundaryKind::pressure_outlet;
    outlet.pressure=pressure;outlet.allow_backflow=true;outlet.backflow_temperature=temperature;
    outlet.scalars={{"A",ScalarBoundaryKind::zero_gradient},{"B",ScalarBoundaryKind::zero_gradient}};
    for(unsigned s=0;s<2;++s) {outlet.scalars[s].backflow_kind=ScalarBoundaryKind::dirichlet;outlet.scalars[s].backflow_value=s==0 ? .3 : trace_species ? 0. : .4;}
  }
  {
    auto& inlet=model.boundaries[2];inlet.flow_kind=BoundaryKind::velocity_inlet;
    inlet.velocity={speed,.1,0};inlet.temperature=temperature;
    inlet.scalars={{"A",ScalarBoundaryKind::dirichlet,.3},{"B",ScalarBoundaryKind::dirichlet,trace_species ? 0. : .4}};
    model.boundaries[3].backflow_velocity={speed,.1,0};
  }
  model.schemes.species=model.schemes.enthalpy=ConvectionScheme::central2;
  model.time.control=TimeControlKind::fixed;
  model.time.initial_dt=model.time.minimum_dt=model.time.maximum_dt=dt;
  model.time.maximum_retries=1;
  const auto property=model.thermophysics.species.front();
  model.thermophysics.species.assign(3,property);
  for(unsigned s=0;s<3;++s) {
    auto& p=model.thermophysics.species[s];p.stable_name=std::string(1,char('A'+s));
    p.viscosity_reference=mu;p.conductivity=gamma*cp;
    p.nasa7_low[5]=p.nasa7_high[5]=s==0 ? formation*mw/kUniversalGasConstant : 0;
  }
  model.transported_scalars={{"A",TransportedScalarRole::species,.7,.7},{"B",TransportedScalarRole::species,.7,.7}};
  model.reaction.mode=ReactionMode::esf_tpdf;model.reaction.mechanism_sha256=gas.identity.mechanism_sha256;
  model.reaction.phase=gas.identity.phase;model.reaction.esf=EsfSpec{};model.reaction.esf->fields=fields;
  CompiledCasePlan plan;auto status=ProductCompiler::compile(MPI_COMM_WORLD,model,{},plan,{&gas,&gas,&gas.closure});
  ProductDriver driver;if(status)status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;if(status)status=driver.restart_expected(expected);
  if(!status) {
    int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0)std::cerr<<"native_ESF setup status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
    return false;
  }
  RestartImage start;start.global_cells=expected.global_cells;start.patch=expected.target_patch;
  start.plan=expected.plan;start.schema=expected.schema;start.geometry=expected.geometry;
  start.dt=dt;start.step=4;start.controller_state=1;start.pressure_reference=pressure;start.backward_euler_recovery=true;
  const auto n=start.patch.cells;
  unsigned independent=0,field_index=0;
  for(std::size_t f=0;f<expected.fields.size;++f) {
    const auto& spec=expected.fields.data[f];RestartImageField field;
    field.role=spec.role;field.field=spec.field;field.components=spec.components;
    field.values.resize(std::size_t(n.x)*n.y*n.z*spec.components);
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      const auto i=std::size_t((z*n.y+y)*n.x+x)*spec.components;const int gx=x+start.patch.begin.x;
      if(spec.role==RestartFieldRole::velocity){field.values[i]=speed;field.values[i+1]=normal_speed(y+start.patch.begin.y+.5);}
      if(spec.role==RestartFieldRole::pressure_absolute)field.values[i]=pressure;
      if(spec.role==RestartFieldRole::enthalpy)field.values[i]=formation*mean_a(gx)+cp*temperature;
      if(spec.role==RestartFieldRole::independent_species)field.values[i]=independent==0 ? mean_a(gx) : species_b(mean_a(gx),y+start.patch.begin.y);
      if(spec.role==RestartFieldRole::stochastic_field) {
        const double a=field_a(gx,field_index,fields);
        field.values[i]=a;field.values[i+1]=species_b(a,y+start.patch.begin.y);field.values[i+2]=trace_species ? 1-a-field.values[i+1] : .3;field.values[i+3]=formation*a+cp*temperature;
      }
      if(spec.role==RestartFieldRole::stochastic_auxiliary) {
        const double a=mean_a(gx);
        field.values[i]=a;field.values[i+1]=species_b(a,y+start.patch.begin.y);field.values[i+2]=trace_species ? 1-a-field.values[i+1] : .3;
        field.values[i+3]=formation*a+cp*temperature;
      }
      if(spec.role==RestartFieldRole::stochastic_transport) {field.values[i]=gamma;field.values[i+2]=mu;}
    }
    if(spec.role==RestartFieldRole::independent_species)++independent;
    if(spec.role==RestartFieldRole::stochastic_field)++field_index;
    start.fields.push_back(std::move(field));
  }
  for(unsigned a=0;a<3;++a) {
    auto extent=n;(a==0 ? extent.x : a==1 ? extent.y : extent.z)++;
    start.final_mass_flux[a].assign(std::size_t(extent.x)*extent.y*extent.z,a==0 ? rho*speed/32 : 0);
    if(a==1)for(int z=0;z<extent.z;++z)for(int y=0;y<extent.y;++y)for(int x=0;x<extent.x;++x)
      start.final_mass_flux[a][x+extent.x*(y+extent.y*z)]=rho*normal_speed(y+start.patch.begin.y)/(4*nx);
  }
  status=driver.initialize_restart(start);DriverStepReport report;

  if(status)status=driver.advance({1,1,1,1,1},report);
  double final_temperature_error=0;
  if(status) {
    RestartSnapshot accepted;status=driver.committed_restart_snapshot(accepted);
    ConstFieldView h{},a{};unsigned random_fields{},auxiliaries{};
    for(std::size_t k=0;k<accepted.fields.size;++k) {
      const auto& f=accepted.fields.data[k];
      if(f.role==RestartFieldRole::enthalpy)h=f.values;
      random_fields+=f.role==RestartFieldRole::stochastic_field;
      auxiliaries+=f.role==RestartFieldRole::stochastic_auxiliary;
      if(f.role==RestartFieldRole::independent_species && !a.base)a=f.values;
      if(f.role==RestartFieldRole::independent_species || f.role==RestartFieldRole::stochastic_field)
        for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x)
          for(unsigned q=0;q<(f.role==RestartFieldRole::stochastic_field ? 3U : 1U);++q)
            composition.push_back(f.values.unchecked({x,y,z},q));
      if(f.role==RestartFieldRole::stochastic_field || f.role==RestartFieldRole::stochastic_auxiliary)
        for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
          Int3 c{x,y,z};
          thermal.push_back((f.values.unchecked(c,3)-formation*f.values.unchecked(c,0))/cp);
          final_temperature_error=std::max(final_temperature_error,std::abs((f.values.unchecked(c,3)-formation*f.values.unchecked(c,0))/cp-temperature));
        }
    }
    if(!h.base || !a.base || random_fields!=fields || auxiliaries!=1)status={StatusCode::invalid_plan,19125};
    if(h.base && a.base)for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      Int3 c{x,y,z};thermal.push_back((h.unchecked(c,0)-formation*a.unchecked(c,0))/cp);final_temperature_error=std::max(final_temperature_error,std::abs((h.unchecked(c,0)-formation*a.unchecked(c,0))/cp-temperature));
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,&final_temperature_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&gas.tuple_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  int okay=status && report.accepted && gas.calls==fields*unsigned(n.x*n.y*n.z) &&
      gas.tuple_error<3e-8 && std::isfinite(final_temperature_error) &&
      (pressure_flux ? final_temperature_error>1e-4 && final_temperature_error<.01 : final_temperature_error<3e-8);
  MPI_Allreduce(MPI_IN_PLACE,&okay,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0)std::cout<<std::setprecision(17)<<"ESF_energy fields="<<fields
      <<" formation="<<formation<<" status="<<unsigned(status.code)<<'/'<<status.detail
      <<" accepted="<<report.accepted<<" transported_T_error="<<gas.tuple_error
      <<" final_T_error="<<final_temperature_error<<" passed="<<okay<<'\n';
  return okay;
}
} // namespace
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  trace_species=argc==2 && std::string_view(argv[1])=="--trace";
  pressure_flux=trace_species || (argc==2 && std::string_view(argv[1])=="--pressure-flux");
  bool passed=true;
  for(unsigned fields:{2U,4U}) {
    formation=0;std::vector<double> baseline,baseT;
    passed=run(fields,baseline,baseT)&&passed;
    for(double shift:{-4000000.,4000000.}) {
      formation=shift;std::vector<double> candidate,testT;
      passed=run(fields,candidate,testT)&&passed;
      double error=baseline.size()==candidate.size() ? 0 : 1;
      if(error==0)for(std::size_t i=0;i<baseline.size();++i)
        error=std::max(error,std::abs(baseline[i]-candidate[i]));
      MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      double thermal_error=baseT.size()==testT.size() && !baseT.empty() ? 0 : 1;
      if(thermal_error==0)for(std::size_t i=0;i<baseT.size();++i)thermal_error=std::max(thermal_error,std::abs(baseT[i]-testT[i]));
      MPI_Allreduce(MPI_IN_PLACE,&thermal_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      passed=error<1e-11 && thermal_error<3e-8 && passed;
      int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if(rank==0)std::cout<<"ESF_energy reference_shift_composition_error="<<error<<" thermal_shift_K="<<thermal_error<<'\n';
    }
  }
  MPI_Finalize();return passed ? 0 : 1;
}
