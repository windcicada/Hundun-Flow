// SPDX-License-Identifier: Apache-2.0
#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include "../../src/models_esf_detail.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

using namespace hundun::v04;
namespace {
constexpr int nx=12;
constexpr double dt=.02, speed=.2, mu=.03, pressure=1e5, temperature=300;
constexpr double mw=28.96546, cp=3.5*kUniversalGasConstant/mw;
constexpr double rho=pressure*mw/(kUniversalGasConstant*temperature);
constexpr double dx=2./nx, gamma=mu/.7, volume=1./(nx*4*4);
using Row=std::array<double,nx>;
using Matrix=std::array<Row,nx>;

// Independent dense elimination of the periodic one-dimensional FV matrix.
// All y/z derivatives vanish in the manufactured three-dimensional field.
Row solve(Matrix a,Row b) {
  for (int k=0;k<nx;++k) {
    int pivot=k;
    for (int i=k+1;i<nx;++i) if(std::abs(a[i][k])>std::abs(a[pivot][k])) pivot=i;
    std::swap(a[k],a[pivot]);std::swap(b[k],b[pivot]);
    for (int i=k+1;i<nx;++i) {
      const double factor=a[i][k]/a[k][k];
      for (int j=k;j<nx;++j) a[i][j]-=factor*a[k][j];
      b[i]-=factor*b[k];
    }
  }
  Row x{};
  for(int i=nx-1;i>=0;--i) {
    double value=b[i];for(int j=i+1;j<nx;++j)value-=a[i][j]*x[j];
    x[i]=value/a[i][i];
  }
  return x;
}
Matrix matrix(double mixing) {
  Matrix a{};
  const double diffusion=dt*gamma/(rho*dx*dx), convection=dt*speed/(2*dx);
  for(int i=0;i<nx;++i) {
    a[i][i]=1+2*diffusion+mixing;
    a[i][(i+nx-1)%nx]=-diffusion-convection;
    a[i][(i+1)%nx]=-diffusion+convection;
  }
  return a;
}
double cosine(int x) {return std::cos(2*std::acos(-1.)*(x+.5)/nx);}
double offset(unsigned f,unsigned fields) {return (2*int(f)-int(fields)+1)*.02;}
double mean_a(int x) {return .3+.03*cosine(x);}
double field_a(int x,unsigned f,unsigned fields) {return mean_a(x)+offset(f,fields)*cosine(x);}

class Gas final : public portable::GasQueryProvider, public portable::GasAdvanceProvider {
 public:
  portable::GasIdentity identity;
  combustion::ChemistryIdentity closure;
  std::vector<Row> implicit, legacy, expected_h;
  bool pressure_case{}, pressure_history{}, reacting{};
  Row expected_mean{};
  static constexpr double pressure_amplitude=1e-5;
  unsigned fields{},local_nx{},begin_x{},calls{};
  double implicit_error{},legacy_error{},matrix_error{},tuple_error{};
  double split_gap{},mean_gap{},history_gap{};
  explicit Gas(unsigned n):fields(n) {
    identity.mechanism_sha256=std::string(64,'b');identity.phase="inert-mixture";
    identity.species_names={"A","B","C"};identity.element_names={"X"};
    identity.element_counts={1,1,1};identity.molecular_weights_kg_per_kmol={mw,mw,mw};
    identity.enthalpy_reference="nasa7-total-enthalpy";identity.composition_fingerprint=4101;
    closure.element_count=1;closure.species={{mw,{1},0},{mw,{1},0},{mw,{1},0}};
    closure.fingerprint=combustion::chemistry_identity_fingerprint(closure);
    identity.closure_fingerprint=closure.fingerprint;
  }
  const portable::GasIdentity& gas_identity() const noexcept override {return identity;}
  portable::Status query_gas(const portable::GasQuery& q,portable::GasQueryOutput& out) noexcept override {
    if(q.species_count!=3 || out.capacity!=3)return portable::Status::invalid_input;
    const double t=q.coordinates==portable::GasStateCoordinates::pressure_enthalpy ? q.enthalpy_j_per_kg/cp : q.temperature_k;
    out.sample={q.revision,q.composition_fingerprint,q.pressure_pa,t,
        q.pressure_pa*mw/(kUniversalGasConstant*t),cp*t,cp,mu,gamma*cp};
    for(unsigned s=0;s<3;++s) {
      out.diffusivities_m2_per_s[s]=gamma/rho;out.species_enthalpies_j_per_kg[s]=cp*t;
      out.net_mass_rates_kg_per_m3_s[s]=0;
    }
    return portable::Status::success;
  }
  portable::Status advance_gas(const portable::GasAdvanceQuery& q,portable::GasAdvanceOutput& out) noexcept override {
    const unsigned f=calls%fields,x=(calls/fields)%local_nx+begin_x;
    ++calls;
    implicit_error=std::max(implicit_error,std::abs(q.state.mass_fractions[0]-implicit[f][x]));
    legacy_error=std::max(legacy_error,std::abs(q.state.mass_fractions[0]-legacy[f][x]));
    tuple_error=std::max({tuple_error,std::abs(q.state.mass_fractions[0]+q.state.mass_fractions[1]-.7),
        std::abs(q.state.mass_fractions[2]-.3),std::abs(q.state.enthalpy_j_per_kg/(cp*temperature)-(pressure_case ? expected_h[f][x] : 1))});
    double d[3],h[3],w[3];portable::GasQueryOutput sample{{},d,h,w,3};
    auto status=query_gas(q.state,sample);if(status!=portable::Status::success)return status;
    for(unsigned s=0;s<3;++s) {const double change=reacting ? q.state.mass_fractions[0]*std::expm1(-2*q.duration_s) : 0;
      const double delta=s==0 ? change : s==1 ? -change : 0;
      out.final_mass_fractions[s]=q.state.mass_fractions[s]+delta;
      out.integrated_species_density_delta_kg_per_m3[s]=sample.sample.density_kg_per_m3*delta;}
    out.final_sample=sample.sample;out.completed_duration_s=q.duration_s;out.internal_step_count=1;
    out.integrated_heat_release_j_per_m3=0;return portable::Status::success;
  }
  bool reference(std::uint64_t seed,std::uint64_t step) {
    const auto w=esf::detail::balanced_wiener(fields,dt,{seed,step,1,0,0,1});
    if(w.status!=portable::Status::success)return false;
    const double beta=mu/std::pow(volume,2./3.),mix=dt*beta/rho;
    Row old_mean{};for(int i=0;i<nx;++i)old_mean[i]=reacting ? .3 : mean_a(i);
    const Row transported_mean=solve(matrix(0),old_mean);
    implicit.resize(fields);legacy.resize(fields);
    std::vector<Row> moved(fields),rhs(fields),unmixed(fields),noise(fields);
    for(unsigned f=0;f<fields;++f) {
      Row old{};for(int i=0;i<nx;++i)old[i]=reacting ? .3+offset(f,fields) : field_a(i,f,fields);
      const double lo=*std::min_element(old.begin(),old.end()),hi=*std::max_element(old.begin(),old.end());
      for(int i=0;i<nx;++i) {
        const double gradient=(old[(i+1)%nx]-old[(i+nx-1)%nx])/(2*dx);
        double delta=std::sqrt(2*gamma/rho)*w.increments[f][0]*gradient;
        // A and B have opposite gradients and complementary bounds; C and
        // h are constant. Their common tuple limiter reduces to this bound.
        if(delta>0)delta=std::min(delta,hi-old[i]);
        if(delta<0)delta=std::max(delta,lo-old[i]);
        noise[f][i]=delta;
        rhs[f][i]=old[i]+delta+mix*transported_mean[i];
        moved[f][i]=old[i]+dt*(-speed*gradient+gamma/rho*
            (old[(i+1)%nx]-2*old[i]+old[(i+nx-1)%nx])/(dx*dx))+delta;
      }
      implicit[f]=solve(matrix(mix),rhs[f]);
      Row transport_rhs{};
      for(int i=0;i<nx;++i)transport_rhs[i]=old[i]+noise[f][i];
      unmixed[f]=solve(matrix(0),transport_rhs);
      const auto a=matrix(mix);
      for(int i=0;i<nx;++i) {
        double r=-rhs[f][i];for(int j=0;j<nx;++j)r+=a[i][j]*implicit[f][j];
        matrix_error=std::max(matrix_error,std::abs(r));
      }
    }
    for(int i=0;i<nx;++i) {
      double m=0;for(unsigned f=0;f<fields;++f)m+=moved[f][i]/fields;
      for(unsigned f=0;f<fields;++f)legacy[f][i]=m+std::exp(-mix)*(moved[f][i]-m);
    }
    // Sensitivity controls: the oracle must distinguish three plausible
    // wiring defects, including advancing accepted history at each jstep.
    Row noisy_mean{};
    for(int i=0;i<nx;++i)for(unsigned f=0;f<fields;++f)noisy_mean[i]+=unmixed[f][i]/fields;
    const Row second_mean=solve(matrix(0),transported_mean);
    for(unsigned f=0;f<fields;++f) {
      Row wrong_mean_rhs{},second_rhs{};
      for(int i=0;i<nx;++i) {
        const double split=transported_mean[i]+std::exp(-mix)*(unmixed[f][i]-transported_mean[i]);
        split_gap=std::max(split_gap,std::abs(split-implicit[f][i]));
        wrong_mean_rhs[i]=(reacting ? .3+offset(f,fields) : field_a(i,f,fields))+noise[f][i]+mix*noisy_mean[i];
        second_rhs[i]=implicit[f][i]+noise[f][i]+mix*second_mean[i];
      }
      const Row wrong_mean=solve(matrix(mix),wrong_mean_rhs),second=solve(matrix(mix),second_rhs);
      for(int i=0;i<nx;++i) {
        mean_gap=std::max(mean_gap,std::abs(wrong_mean[i]-implicit[f][i]));
        history_gap=std::max(history_gap,std::abs(second[i]-implicit[f][i]));
      }
    }
    if(reacting) {
      Row rhs=old_mean;
      for(int x=0;x<nx;++x)for(unsigned f=0;f<fields;++f)
        rhs[x]+=std::expm1(-2*dt)*implicit[f][x]/fields;
      expected_mean=solve(matrix(0),rhs);
    }
    if(pressure_case) {
      Row old{},work{},rhs{};
      for(int x=0;x<nx;++x) {
        old[x]=1+pressure_amplitude*cosine(x);
        // dp/dt is zero for recovered history; the periodic advective work
        // uses an independently written central pressure difference.
        work[x]=dt*speed*pressure*pressure_amplitude*
            (cosine((x+1)%nx)-cosine((x+nx-1)%nx))/(2*dx*rho*cp*temperature);
        if(pressure_history)work[x]+=dt*.1/(.01*rho*cp*temperature);
        rhs[x]=old[x]+work[x];
      }
      const auto mean_h=solve(matrix(0),rhs);
      expected_h.resize(fields);
      const double lo=*std::min_element(old.begin(),old.end()),hi=*std::max_element(old.begin(),old.end());
      for(unsigned f=0;f<fields;++f) {
        for(int x=0;x<nx;++x) {
          double noise=std::sqrt(2*gamma/rho)*w.increments[f][0]*
              (old[(x+1)%nx]-old[(x+nx-1)%nx])/(2*dx);
          noise=std::max(lo-old[x],std::min(hi-old[x],noise));
          rhs[x]=old[x]+noise+mix*mean_h[x]+work[x];
        }
        expected_h[f]=solve(matrix(mix),rhs);
      }
    }
    return matrix_error<1e-13 && (reacting || (split_gap>1e-8 && mean_gap>1e-8 && history_gap>1e-8));
  }
};

bool run(unsigned fields,bool require_implicit,bool pressure_case=false,bool pressure_history=false,bool cn=false,bool reacting=false) {
  Gas gas(fields);gas.pressure_case=pressure_case;gas.pressure_history=pressure_history;gas.reacting=reacting;
  auto model=test::product_model({nx,4,4});
  if(cn) {
    model.time.scheme=TimeScheme::cn_be;model.solver.coupling=CouplingKind::outer_corrected;
    model.pressure_reference=PressureReferenceKind::boundary_absolute;
    model.legacy_time_fingerprint=model.fingerprint+1;
  }
  model.turbulence=TurbulenceKind::none;
  if(reacting) {
    model.boundaries[2].flow_kind=BoundaryKind::symmetry;
    model.boundaries[2].scalars={{"A",ScalarBoundaryKind::zero_gradient},{"B",ScalarBoundaryKind::zero_gradient}};
    auto& outlet=model.boundaries[3];outlet.flow_kind=BoundaryKind::pressure_outlet;
    outlet.pressure=pressure;outlet.allow_backflow=true;outlet.backflow_temperature=temperature;
    outlet.scalars={{"A",ScalarBoundaryKind::zero_gradient},{"B",ScalarBoundaryKind::zero_gradient}};
    for(unsigned s=0;s<2;++s) {outlet.scalars[s].backflow_kind=ScalarBoundaryKind::dirichlet;outlet.scalars[s].backflow_value=s==0 ? .3 : .4;}
  }
  model.schemes.species=model.schemes.enthalpy=ConvectionScheme::central2;
  model.time.control=TimeControlKind::fixed;
  model.time.initial_dt=model.time.minimum_dt=model.time.maximum_dt=dt;
  model.time.maximum_retries=1;
  if(pressure_history) {model.time.minimum_dt=.01;model.time.maximum_growth=2;}
  const auto property=model.thermophysics.species.front();
  model.thermophysics.species.assign(3,property);
  for(unsigned s=0;s<3;++s) {
    auto& p=model.thermophysics.species[s];p.stable_name=std::string(1,char('A'+s));
    p.viscosity_reference=mu;p.conductivity=gamma*cp;
  }
  model.transported_scalars={{"A",TransportedScalarRole::species,.7,.7},{"B",TransportedScalarRole::species,.7,.7}};
  model.reaction.mode=ReactionMode::esf_tpdf;model.reaction.mechanism_sha256=gas.identity.mechanism_sha256;
  model.reaction.phase=gas.identity.phase;model.reaction.esf=EsfSpec{};model.reaction.esf->fields=fields;
  CompiledCasePlan plan;auto status=ProductCompiler::compile(MPI_COMM_WORLD,model,{},plan,{&gas,&gas,&gas.closure});
  std::vector<FieldId> transport_work;
  if(status) {
    // Production registration owns one source/sink pair per independent
    // species and h. Source ownership and units survive scalar view aliasing.
    unsigned sources=0,means=0,heat_sources=0;std::vector<FieldId> outputs;
    for(const auto& field:*plan.field_schema()) {
      if(field.stable_name=="esf_mean" || field.stable_name=="esf_old" || field.stable_name=="esf_iter") {
        means+=field.components==4 && field.ghost_width==2;
        transport_work.push_back(field.id);
      }
    }
    const auto contributions=plan.contribution_plan()->contributions();
    for(std::size_t i=0;i<contributions.size;++i) {
      const auto& c=contributions.data[i];
      if(c.stage!=178)continue;
      bool heat=false;
      for(const auto& field:*plan.field_schema())
        if(field.id==c.conserved_quantity)heat=field.stable_name=="h";
      heat_sources+=heat;
      const std::array<std::int8_t,7> units=heat ? std::array<std::int8_t,7>{1,-1,-3,0,0,0,0} :
          std::array<std::int8_t,7>{1,-3,-1,0,0,0,0};
      if(!c.supplies_implicit_diagonal || c.units.si_exponents!=units ||
          c.capability!=ContributionCapability::reacting || !c.source_identity)
        status={StatusCode::invalid_plan,19120};
      outputs.push_back(c.explicit_source);outputs.push_back(c.implicit_diagonal);++sources;
    }
    std::sort(outputs.begin(),outputs.end());
    if(sources!=3 || means!=3 || heat_sources!=1 ||
        std::adjacent_find(outputs.begin(),outputs.end())!=outputs.end())
      status={StatusCode::invalid_plan,19121};
    transport_work.insert(transport_work.end(),outputs.begin(),outputs.end());
    for(auto id:transport_work) {
      const auto* allocation=plan.arena_layout()->field(id);
      if(!allocation || allocation->lifetime!=FieldLifetime::persistent_workspace || allocation->replicas!=1)
        status={StatusCode::invalid_plan,19122};
    }
  }
  ProductDriver driver;if(status)status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;if(status)status=driver.restart_expected(expected);
  if(!status) {
    int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0)std::cerr<<"native_ESF setup status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
    return false;
  }
  for(std::size_t f=0;f<expected.fields.size;++f)
    if(std::find(transport_work.begin(),transport_work.end(),expected.fields.data[f].field)!=transport_work.end())
      return false;
  RestartImage start;start.global_cells=expected.global_cells;start.patch=expected.target_patch;
  start.plan=expected.plan;start.schema=expected.schema;start.geometry=expected.geometry;
  start.dt=dt;start.step=4;start.controller_state=1;start.pressure_reference=pressure;start.backward_euler_recovery=true;
  const auto n=start.patch.cells;gas.local_nx=n.x;gas.begin_x=start.patch.begin.x;
  if(!gas.reference(model.reaction.esf->seed,start.step))return false;
  unsigned independent=0,field_index=0;
  for(std::size_t f=0;f<expected.fields.size;++f) {
    const auto& spec=expected.fields.data[f];RestartImageField field;
    field.role=spec.role;field.field=spec.field;field.components=spec.components;
    field.values.resize(std::size_t(n.x)*n.y*n.z*spec.components);
    for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
      const auto i=std::size_t((z*n.y+y)*n.x+x)*spec.components;const int gx=x+start.patch.begin.x;
      if(spec.role==RestartFieldRole::velocity)field.values[i]=speed;
      if(spec.role==RestartFieldRole::pressure_absolute)field.values[i]=pressure*
          (1+(pressure_case ? Gas::pressure_amplitude*cosine(gx) : 0));
      if(spec.role==RestartFieldRole::pressure_perturbation && pressure_case)
        field.values[i]=pressure*Gas::pressure_amplitude*cosine(gx);
      if(spec.role==RestartFieldRole::enthalpy)field.values[i]=cp*temperature*
          (1+(pressure_case ? Gas::pressure_amplitude*cosine(gx) : 0));
      if(spec.role==RestartFieldRole::independent_species)field.values[i]=independent==0 ? (reacting ? .3 : mean_a(gx)) : .7-(reacting ? .3 : mean_a(gx));
      if(spec.role==RestartFieldRole::stochastic_field) {
        const double a=reacting ? .3+offset(field_index,fields) : field_a(gx,field_index,fields);
        field.values[i]=a;field.values[i+1]=.7-a;field.values[i+2]=.3;field.values[i+3]=cp*temperature*(1+(pressure_case ? Gas::pressure_amplitude*cosine(gx) : 0));
      }
      if(spec.role==RestartFieldRole::stochastic_auxiliary) {
        const double a=reacting ? .3 : mean_a(gx);
        field.values[i]=a;field.values[i+1]=.7-a;field.values[i+2]=.3;
        field.values[i+3]=cp*temperature*(1+(pressure_case ? Gas::pressure_amplitude*cosine(gx) : 0));
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
  }
  if(pressure_history) {
    start.source_format_version=3;start.backward_euler_recovery=false;
    start.method_history_signature=expected.method_history_signature;
    start.time=.04;start.dt=.01;
    start.previous_pressure_reference=pressure-.1;
    start.closed_mass_target=rho;
    start.previous_mass_flux_revision=1;start.final_mass_flux_revision=2;
    start.previous_fields=start.fields;start.previous_mass_flux=start.final_mass_flux;
    for(std::size_t f=0;f<expected.rate_fields.size;++f) {
      const auto& spec=expected.rate_fields.data[f];RestartImageField rate;
      rate.role=spec.role;rate.field=spec.field;rate.components=spec.components;
      rate.values.assign(std::size_t(n.x)*n.y*n.z*spec.components,0.);
      start.accepted_rate_fields.push_back(rate);start.previous_rate_fields.push_back(rate);
    }
  }
  status=driver.initialize_restart(start);DriverStepReport report;

  if(status)status=driver.advance({1,1,1,1,1},report);
  double mean_error=0,stored_rate_error=0;
  if(status && reacting) {
    RestartSnapshot accepted;status=driver.committed_restart_snapshot(accepted);
    bool checked=false;
    if(status)for(std::size_t i=0;i<accepted.fields.size;++i) {
      const auto& field=accepted.fields.data[i];
      if(field.role!=RestartFieldRole::independent_species || checked)continue;
      checked=true;
      for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x)
        mean_error=std::max(mean_error,std::abs(field.values.unchecked({x,y,z},0)-gas.expected_mean[x+start.patch.begin.x]));
    }
    if(!checked)mean_error=1;
    // Uniform inert-isomer transport has zero endpoint EX2 rates. The
    // integrated chemistry source belongs only to this step's mean BE solve.
    unsigned rate_fields=0;
    if(status)for(std::size_t i=0;i<accepted.accepted_rate_fields.size;++i) {
      const auto& field=accepted.accepted_rate_fields.data[i];
      if(field.role!=RestartFieldRole::scalar_nonadvective_rate)continue;
      ++rate_fields;
      for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x)
        stored_rate_error=std::max(stored_rate_error,std::abs(field.values.unchecked({x,y,z},0)));
    }
    if(rate_fields!=2)stored_rate_error=1;
    MPI_Allreduce(MPI_IN_PLACE,&mean_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&stored_rate_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  }
  double errors[]{gas.implicit_error,gas.legacy_error,gas.matrix_error,gas.tuple_error};
  MPI_Allreduce(MPI_IN_PLACE,errors,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  const unsigned batch=fields*unsigned(n.x*n.y*n.z);
  // Pressure-dependent outer candidates may repeat chemistry from the same
  // accepted state. Every complete evaluation is checked by the callback.
  const unsigned evaluations=gas.calls/batch;
  unsigned minimum=evaluations,maximum=evaluations;
  MPI_Allreduce(MPI_IN_PLACE,&minimum,1,MPI_UNSIGNED,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&maximum,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);
  int okay=status && report.accepted && gas.calls%batch==0 && minimum==maximum &&
      ((pressure_case || reacting) ? evaluations>=1 : evaluations==1) && mean_error<2e-12 && stored_rate_error<2e-12 &&
      errors[require_implicit ? 0 : 1]<2e-12 && errors[3]<2e-12;
  MPI_Allreduce(MPI_IN_PLACE,&okay,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0)std::cout<<std::setprecision(17)<<"native_ESF fields="<<fields<<" status="<<unsigned(status.code)<<'/'<<status.detail
      <<" accepted="<<report.accepted<<" attempts="<<report.attempts<<" dt="<<report.proposal.dt<<" calls_local="<<gas.calls
      <<" implicit_error="<<errors[0]<<" legacy_error="<<errors[1]<<" matrix_residual="<<errors[2]
      <<" tuple_error="<<errors[3]<<" mean_error="<<mean_error<<" stored_rate_error="<<stored_rate_error
      <<" split_gap="<<gas.split_gap<<" mean_gap="<<gas.mean_gap
      <<" history_gap="<<gas.history_gap<<" passed="<<okay<<'\n';
  return okay;
}
} // namespace
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  const bool implicit=!(argc==2 && std::string_view(argv[1])=="--legacy");
  const bool reacting=argc==2 && std::string_view(argv[1])=="--reacting-cn";
  const bool cn=reacting || (argc==2 && std::string_view(argv[1])=="--pressure-cn");
  const bool pressure_history=(cn && !reacting) || (argc==2 && std::string_view(argv[1])=="--pressure-history");
  const bool pressure_case=pressure_history || (argc==2 && std::string_view(argv[1])=="--pressure");
  bool passed=run(2,implicit,pressure_case,pressure_history,cn,reacting);
  passed=run(4,implicit,pressure_case,pressure_history,cn,reacting)&&passed;
  MPI_Finalize();return passed ? 0 : 1;
}
