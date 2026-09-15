// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "models_spray_properties_detail.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string_view>
using namespace hundun::v04;
using namespace spray::detail;
namespace {
constexpr double gas_r=8314.46261815324;
class Transport final : public FilmTransportProvider {
 public:
  bool stale{}, fail{};
  mutable unsigned calls{};
  mutable double temperature{};
  FilmTransportReport query(const portable::GasQuery& q) const noexcept override {
    ++calls; temperature=q.temperature_k;
    FilmTransportReport out;
    if (fail || q.coordinates!=portable::GasStateCoordinates::pressure_temperature ||
        q.composition_fingerprint!=42 || q.species_count!=3) return out;
    out.status=portable::Status::success;out.revision=q.revision;
    if(stale) ++out.revision.input_revision;
    out.dynamic_viscosity_pa_s=2e-5;
    return out;
  }
};
class Gas final : public portable::GasQueryProvider {
 public:
  portable::GasIdentity id;
  unsigned calls{},fail_at{},stale_at{},nan_at{};
  std::array<double,4> queried_temperature{},queried_vapor{};
  Gas() {
    id.mechanism_sha256=std::string(64,'c');id.phase="kerosene-operator";
    id.species_names={"O2","N2","C12H23"};id.molecular_weights_kg_per_kmol={32.,28.,167.};
    id.enthalpy_reference="nasa7-total-enthalpy";id.composition_fingerprint=42;
  }
  const portable::GasIdentity& gas_identity() const noexcept override {return id;}
  double a(unsigned i) const {return gas_r*(i==2 ? 30. : 3.5)/id.molecular_weights_kg_per_kmol[i];}
  double b(unsigned i) const {return gas_r*(i==2 ? .002 : i==1 ? .0002 : .0001)/id.molecular_weights_kg_per_kmol[i];}
  double cp(unsigned i,double t) const {return a(i)+b(i)*t;}
  double h(unsigned i,double t) const {return a(i)*t+.5*b(i)*t*t-(i==2 ? 1e6 : 0.);}
  double mixture_cp(const double* y,double t) const {return y[0]*cp(0,t)+y[1]*cp(1,t)+y[2]*cp(2,t);}
  double mixture_h(const double* y,double t) const {return y[0]*h(0,t)+y[1]*h(1,t)+y[2]*h(2,t);}
  portable::Status query_gas(const portable::GasQuery& q,portable::GasQueryOutput& out) noexcept override {
    const unsigned index=calls++%4;
    if(calls==fail_at)return portable::Status::provider_failure;
    double ca=0,cb=0;
    for(unsigned i=0;i<3;++i) {ca+=q.mass_fractions[i]*a(i);cb+=q.mass_fractions[i]*b(i);}
    const double sensible=q.enthalpy_j_per_kg+1e6*q.mass_fractions[2];
    const double t=q.coordinates==portable::GasStateCoordinates::pressure_enthalpy
        ? 2*sensible/(ca+std::sqrt(ca*ca+2*cb*sensible)) : q.temperature_k;
    const double c=mixture_cp(q.mass_fractions,t);
    queried_temperature[index]=t;queried_vapor[index]=q.mass_fractions[2];
    out.sample={q.revision,q.composition_fingerprint,q.pressure_pa,t,1.,mixture_h(q.mass_fractions,t),c,9e-6,.025};
    if(calls==stale_at)++out.sample.revision.input_revision;
    if(calls==nan_at)out.sample.cp_j_per_kg_k=std::numeric_limits<double>::quiet_NaN();
    for(unsigned i=0;i<3;++i) {out.diffusivities_m2_per_s[i]=1e-5;out.species_enthalpies_j_per_kg[i]=h(i,t);out.net_mass_rates_kg_per_m3_s[i]=0;}
    return portable::Status::success;
  }
};
}
int main(int argc,char** argv) {
  Gas gas; LiquidAsset asset;asset.gas_identity=gas.id;
  asset.vapor_species_index=2;asset.vapor_species_name="C12H23";
  asset.vapor_molecular_weight_kg_per_kmol=167.;
  asset.pack.minimum_temperature_k=200.;asset.pack.maximum_temperature_k=680.;
  double y[]{.23,.76,.01};
  KeroseneFilmInput input;
  input.expected_revision=input.transport_revision={32004,3,1};
  input.far_gas={input.expected_revision,42,portable::GasStateCoordinates::pressure_enthalpy,
      101325.,gas.mixture_h(y,900.),0.,y,3};
  input.surface_temperature_k=350.;input.far_dynamic_viscosity_pa_s=2e-5;
  KeroseneFilmWorkspace workspace(3);
  if(argc==2 && std::string_view(argv[1])=="--sample") {
    double t,p,unused;std::cout<<std::setprecision(17);
    while(std::cin>>t>>p>>unused) {
      input.surface_temperature_k=t;input.far_gas.pressure_pa=p;gas.calls=0;
      const auto r=workspace.query(asset,gas,input);if(!r.available)return 2;
      std::cout<<r.film_temperature_k<<' '<<r.gas_cp_j_per_kg_k<<' '
          <<r.gas_dynamic_viscosity_pa_s<<' '<<r.vapor_cp_j_per_kg_k<<' '
          <<r.vapor_prandtl_number<<' '<<r.vapor_absolute_enthalpy_j_per_kg<<' '
          <<r.surface_vapor_mass_fraction<<'\n';
    }
    return 0;
  }
  auto r=workspace.query(asset,gas,input);
  if(!r.available){std::cerr<<"THICK_EX gas-property query unavailable\n";return 1;}
  bool passed=r.revision==input.expected_revision && gas.calls==4 &&
      std::abs(r.film_temperature_k-(350.+550./3.))<1e-12 &&
      gas.queried_temperature[1]==298.15 && gas.queried_temperature[2]==350. &&
      gas.queried_temperature[3]==r.film_temperature_k && gas.queried_vapor[1]==1. && gas.queried_vapor[2]==1. && gas.queried_vapor[3]==.01 &&
      std::abs(r.gas_cp_j_per_kg_k-gas.mixture_cp(y,r.film_temperature_k))<1e-12 &&
      r.vapor_absolute_enthalpy_j_per_kg==gas.h(2,350.) && workspace.owned_payload_bytes()==96;
  const auto reference=r;
  input.far_dynamic_viscosity_pa_s=3e-5;gas.calls=0;
  r=workspace.query(asset,gas,input);
  passed &= r.available && r.gas_dynamic_viscosity_pa_s>reference.gas_dynamic_viscosity_pa_s &&
      r.vapor_prandtl_number==reference.vapor_prandtl_number;
  input.far_dynamic_viscosity_pa_s=2e-5;
  for(unsigned at=1;at<=4;++at) {
    for(unsigned kind=0;kind<3;++kind) {
      gas.calls=0;gas.fail_at=kind==0 ? at : 0;gas.stale_at=kind==1 ? at : 0;gas.nan_at=kind==2 ? at : 0;
      r=workspace.query(asset,gas,input);
      passed &= !r.available && r.revision.input_revision==0 && r.gas_cp_j_per_kg_k==0;
    }
  }
  gas.calls=gas.fail_at=gas.stale_at=gas.nan_at=0;
  input.transport_revision.input_revision++;
  passed &= workspace.query(asset,gas,input).status==portable::Status::stale_revision && gas.calls==0;
  input.transport_revision=input.expected_revision;
  asset.vapor_species_index=0;
  passed &= workspace.query(asset,gas,input).status==portable::Status::identity_mismatch;
  asset.vapor_species_index=2;
  KeroseneFilmWorkspace short_workspace(2);
  passed &= short_workspace.query(asset,gas,input).status==portable::Status::capacity_exceeded;
  y[2]=.02;passed &= !workspace.query(asset,gas,input).available;y[2]=.01;
  gas.calls=0;r=workspace.query(asset,gas,input);
  passed &= r.available && r.gas_dynamic_viscosity_pa_s==reference.gas_dynamic_viscosity_pa_s && gas.calls==4;
  Transport transport;
  input.transport=&transport;input.far_dynamic_viscosity_pa_s=0.;gas.calls=0;
  r=workspace.query(asset,gas,input);
  passed &= r.available && gas.calls==4 && transport.calls==1 &&
      std::abs(transport.temperature-900.)<1e-12 &&
      r.gas_dynamic_viscosity_pa_s==reference.gas_dynamic_viscosity_pa_s &&
      r.far_dynamic_viscosity_pa_s==2e-5 && r.far_density_kg_per_m3==1.;
  transport.stale=true;
  passed &= workspace.query(asset,gas,input).status==portable::Status::stale_revision;
  transport.stale=false;transport.fail=true;
  passed &= !workspace.query(asset,gas,input).available;
  transport.fail=false;input.far_dynamic_viscosity_pa_s=2e-5;
  passed &= !workspace.query(asset,gas,input).available;
  std::cout<<"thick_gas_queries_identity_and_failure passed="<<passed<<'\n';
  return passed ? 0 : 1;
}
