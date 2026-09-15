// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "hundun/v04_cantera.hpp"
#include "../support/chemistry_test_support.hpp"
#include <array>
#include <numeric>
using namespace hundun::v04;
int main(int argc,char** argv) {
 return hundun::test::run([&] {
  const bool reject=argc==4 && std::string(argv[1])=="--reject";
  HUNDUN_CHECK(argc==3 || reject);
  chemistry::CanteraBackendConfig config;
  config.mechanism={argv[reject ? 2 : 1],argv[reject ? 3 : 2],"kerosene"};
  config.species_names={"H2","H2O","CO","CO2","O2","N2","C12H23"};
  config.chemistry={1e-10,1e-16,20000};config.continuous_enthalpy=true;
  if(reject) {
    bool rejected=false;
    try { chemistry::CanteraBackendRuntime runtime(config); }
    catch(const std::exception& error) {
      rejected=true;std::cout<<error.what()<<'\n';
    }
    HUNDUN_CHECK(rejected);
    return;
  }
  auto runtime=std::make_shared<chemistry::CanteraBackendRuntime>(config);
  HUNDUN_CHECK(runtime->reaction_model()=="kerosene_4_step_v1/frozen_material");
  chemistry::CanteraWorkspacePool pool(runtime,2);
  auto backend=chemistry::make_cantera_backend(config,pool);
  auto other=chemistry::make_cantera_backend(config,pool);
  HUNDUN_CHECK(pool.workspaces_are_distinct());
  std::array<double,7> ys{.001,.05,.005,.01,.2,.714,.02},d{},h{},rates{};
  portable::GasQuery q;q.revision.algorithm_version=1;
  q.composition_fingerprint=backend->composition().fingerprint;
  q.coordinates=portable::GasStateCoordinates::pressure_temperature;
  q.pressure_pa=100000;q.temperature_k=1200;q.mass_fractions=ys.data();q.species_count=7;
  portable::GasQueryOutput sample{{},d.data(),h.data(),rates.data(),7};
  HUNDUN_CHECK(backend->query_gas(q,sample)==portable::Status::success);
  const double magnitude=std::accumulate(rates.begin(),rates.end(),0.,[](double s,double v){return s+std::abs(v);});
  std::cout<<"kerosene gas_query source_magnitude="<<magnitude<<'\n';
  HUNDUN_CHECK(magnitude>1.);
  HUNDUN_CHECK(std::abs(std::accumulate(rates.begin(),rates.end(),0.))<1e-12*magnitude);
  const double initial_h=sample.sample.enthalpy_j_per_kg;
  q.coordinates=portable::GasStateCoordinates::pressure_enthalpy;q.enthalpy_j_per_kg=initial_h;
  portable::GasAdvanceQuery request{q,1000.,1e-5};
  std::array<double,7> final{},delta{};
  portable::GasAdvanceOutput result{{},final.data(),delta.data(),7};
  HUNDUN_CHECK(backend->advance_gas(request,result)==portable::Status::success);
  HUNDUN_CHECK(result.completed_duration_s==request.duration_s && result.internal_step_count>0);
  HUNDUN_CHECK(final[6]<ys[6] && result.final_sample.temperature_k>q.temperature_k);
  HUNDUN_CHECK(result.final_sample.enthalpy_j_per_kg==initial_h);
  HUNDUN_CHECK(std::abs(std::accumulate(final.begin(),final.end(),0.)-1)<2e-12);
  HUNDUN_CHECK(std::abs(std::accumulate(delta.begin(),delta.end(),0.))<1e-12);
  for(double value:final)HUNDUN_CHECK(value>=0. && value<=1.);
  const auto saved=final;
  auto changed=request;changed.state.temperature_k=1800.;changed.state.enthalpy_j_per_kg+=2e5;
  HUNDUN_CHECK(other->advance_gas(changed,result)==portable::Status::success);
  HUNDUN_CHECK(backend->advance_gas(request,result)==portable::Status::success);
  for(unsigned i=0;i<7;++i)HUNDUN_CHECK(std::abs(final[i]-saved[i])<1e-12);
  auto restricted=config;restricted.chemistry.maximum_internal_steps=1;
  chemistry::CanteraWorkspacePool restricted_pool(runtime,1);
  auto limited=chemistry::make_cantera_backend(restricted,restricted_pool);
  final.fill(71.);delta.fill(73.);
  auto too_long=request;too_long.duration_s=1e-3;
  HUNDUN_CHECK(limited->advance_gas(too_long,result)!=portable::Status::success);
  HUNDUN_CHECK(result.completed_duration_s==0. && result.final_sample.temperature_k==0.);
  for(unsigned i=0;i<7;++i)HUNDUN_CHECK(final[i]==71. && delta[i]==73.);
  HUNDUN_CHECK(backend->advance_gas(request,result)==portable::Status::success);
  for(unsigned i=0;i<7;++i)HUNDUN_CHECK(std::abs(final[i]-saved[i])<1e-12);
  chemistry::ChemistryIntervalRequest legacy;
  legacy.state.p0_pa=q.pressure_pa;legacy.state.h_tc_j_per_kg=initial_h;
  legacy.state.mass_fractions.assign(ys.begin(),ys.end());legacy.start_time_s=1000.;legacy.duration_s=1e-5;
  const auto old=backend->integrate(legacy);HUNDUN_CHECK(old.succeeded());
  for(unsigned i=0;i<7;++i)HUNDUN_CHECK(std::abs(old.final_state.mass_fractions[i]-saved[i])<1e-10);
 });
}
