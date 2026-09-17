// SPDX-License-Identifier: Apache-2.0
// Check output closures independently from accepted native Restart tuples.
#include "hundun/v04_app.hpp"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace hundun::v04;
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  struct End { ~End(){MPI_Finalize();} } end;
  if(argc!=3)return 2;
  int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const auto agree=[](Status s) {
    int ok=bool(s);MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    return bool(ok);
  };
  ValidatedModel model;
  auto status=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[1],model);
  CompiledCasePlan plan;
  if(status)status=ProductCompiler::compile(MPI_COMM_WORLD,model,argv[1],plan);
  ProductDriver driver;
  if(status)status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;RestartImage image;
  if(status)status=driver.restart_expected(expected);
  if(status)status=RestartReader::load(MPI_COMM_WORLD,argv[2],expected,image);
  if(status)status=driver.initialize_restart(image);
  ThermodynamicsPlan thermo;
  if(status)status=ThermodynamicsPlan::compile(model.thermophysics,
      {model.transported_scalars.data(),model.transported_scalars.size()},thermo);
  if(!agree(status))return 3;
  CommittedOutputSnapshot snapshot;
  status=driver.committed_sgs_output_snapshot(snapshot);
  if(!agree(status)) {std::cerr<<unsigned(status.code)<<'/'<<status.detail<<'\n';return 4;}
  const auto field=[&](std::string_view name) {
    for(std::size_t i=0;i<snapshot.fields.size;++i) {
      const auto& f=snapshot.fields.data[i];
      if(f.stable_name==name)return f.values;
    }
    throw std::runtime_error("missing output field");
  };
  const auto ns=model.thermophysics.species.size();
  std::vector<std::size_t> mapping;
  for(const auto& scalar:model.transported_scalars)if(scalar.role==TransportedScalarRole::species)
    for(std::size_t j=0;j<ns;++j)if(model.thermophysics.species[j].stable_name==scalar.stable_name)mapping.push_back(j);
  std::vector<const RestartImageField*> pdf;
  const RestartImageField* auxiliary=nullptr;
  for(const auto& f:image.fields) {
    if(f.role==RestartFieldRole::stochastic_field)pdf.push_back(&f);
    if(f.role==RestartFieldRole::stochastic_auxiliary)auxiliary=&f;
  }
  if(!auxiliary || pdf.empty())return 5;
  double error{},gap{};
  const auto check=[&](double a,double b) {
    error=std::max(error,std::isfinite(a)&&std::isfinite(b)
        ? std::abs(a-b)/std::max({1.,std::abs(a),std::abs(b)}) : INFINITY);
  };
  const auto pi=field("pi"),rho=field("rho"),temperature=field("T");
  const auto mean_rho=field("rho_mean_eos"),stat_rho=field("rho_pdf_mean"),rho0=field("rho_field0");
  const auto t0=field("T_field0"),h0=field("h_field0"),y0=field("Y_field0");
  std::vector<double> y(mapping.size()),mean(ns+1U);
  const auto n=image.patch.cells;
  std::size_t cell{};
  for(int z=0;z<n.z;++z)for(int yy=0;yy<n.y;++yy)for(int x=0;x<n.x;++x,++cell) {
    const Int3 c{x,yy,z};
    if(snapshot.cell_activity.size && snapshot.cell_activity.data[cell]==0U)continue;
    const auto p=thermo.eos_pressure(image.pressure_reference+pi.unchecked(c,0));
    const auto evaluate=[&](const double* row,double weight) {
      for(std::size_t s=0;s<y.size();++s)y[s]=std::max(0.,row[mapping[s]])/weight;
      ThermoState state;
      auto result=thermo.evaluate(p,row[ns]/weight,{y.data(),y.size()},{},state);
      if(!result)error=INFINITY;
      return state;
    };
    std::fill(mean.begin(),mean.end(),0.);
    long double specific_volume{};
    for(const auto* f:pdf) {
      const auto* row=f->values.data()+cell*(ns+1);
      specific_volume+=1.L/evaluate(row,1.).rho/pdf.size();
      for(std::size_t s=0;s<=ns;++s)mean[s]+=row[s]/pdf.size();
    }
    const auto physical=evaluate(mean.data(),1.);
    const auto* row=auxiliary->values.data()+cell*(ns+1);
    long double weight{};for(std::size_t s=0;s<ns;++s)weight+=std::max(0.,row[s]);
    const auto auxiliary_state=evaluate(row,static_cast<double>(weight));
    check(mean_rho.unchecked(c,0),physical.rho);
    check(temperature.unchecked(c,0),physical.temperature);
    check(stat_rho.unchecked(c,0),static_cast<double>(1.L/specific_volume));
    check(rho0.unchecked(c,0),auxiliary_state.rho/static_cast<double>(weight));
    check(rho0.unchecked(c,0),rho.unchecked(c,0));
    check(t0.unchecked(c,0),auxiliary_state.temperature);
    check(h0.unchecked(c,0),row[ns]);
    for(std::size_t s=0;s<ns;++s)check(y0.unchecked(c,s),row[s]);
    gap=std::max(gap,std::abs(mean_rho.unchecked(c,0)-stat_rho.unchecked(c,0)));
  }
  MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&gap,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(!rank)std::cout<<"ESF output closure error="<<error<<" distinct_density_gap="<<gap<<'\n';
  return error<1e-11 && gap>1e-12 ? 0 : 6;
}
