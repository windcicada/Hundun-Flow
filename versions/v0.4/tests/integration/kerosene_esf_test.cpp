// SPDX-License-Identifier: Apache-2.0
// Small native 624CF gas fixture: full histories and absolute face-flux audit.
#include "hundun/v04_app.hpp"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
using namespace hundun::v04;
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);struct End{~End(){MPI_Finalize();}} end;int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const auto all=[&](bool value){int ok=value;MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);return bool(ok);};
  if(argc!=4)return 2;
  ValidatedModel model;auto s=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[1],model);
  CompiledCasePlan plan;if(s)s=ProductCompiler::compile(MPI_COMM_WORLD,model,argv[1],plan);
  ProductDriver driver;if(s)s=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
  RestartExpected expected;if(s)s=driver.restart_expected(expected);
  RestartImage a,b;if(s)s=RestartReader::load(MPI_COMM_WORLD,argv[2],expected,a);
  if(s)s=RestartReader::load(MPI_COMM_WORLD,argv[3],expected,b);
  if(!all(bool(s))) {if(!rank)std::cout<<"load "<<unsigned(s.code)<<'/'<<s.detail<<std::endl;return 3;}
  const bool metadata=a.step==b.step && a.time==b.time && a.dt==b.dt && a.method_history_signature==b.method_history_signature && a.plan==b.plan && a.schema==b.schema && a.geometry==b.geometry && a.source_format_version==b.source_format_version && !a.backward_euler_recovery && !b.backward_euler_recovery && a.controller_state==b.controller_state && a.pressure_reference==b.pressure_reference && a.previous_pressure_reference==b.previous_pressure_reference && a.closed_mass_target==b.closed_mass_target && a.cell_records==b.cell_records;
  if(!all(metadata))return 4;
  bool passed=true;double worst{};
  const auto role_count=[&](RestartFieldRole role) {
    return std::count_if(a.fields.begin(),a.fields.end(),[&](const auto& f){return f.role==role;});
  };
  if(!all(role_count(RestartFieldRole::stochastic_field)==2 &&
          role_count(RestartFieldRole::stochastic_auxiliary)==1))return 8;
  double maximum_flux_difference{};unsigned long long total{};
  const auto fields=[&](const char* level,const std::vector<RestartImageField>& x,const std::vector<RestartImageField>& y) {
    if(!all(x.size()==y.size()))return false;
    for(std::size_t f=0;f<x.size();++f) {
      if(!all(x[f].field==y[f].field && x[f].role==y[f].role && x[f].components==y[f].components && x[f].values.size()==y[f].values.size()))return false;
      for(unsigned c=0;c<x[f].components;++c) {
        double v[2]{};unsigned long long differences{};
        for(std::size_t i=c;i<x[f].values.size();i+=x[f].components) {
          const double aa=x[f].values[i],bb=y[f].values[i];v[0]=std::max(v[0],std::abs(aa-bb));v[1]=std::max({v[1],std::abs(aa),std::abs(bb)});differences+=aa!=bb;
          if(!std::isfinite(aa)||!std::isfinite(bb))v[0]=INFINITY;
          ++total;
        }
        MPI_Allreduce(MPI_IN_PLACE,v,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&differences,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
        const double relative=v[0]/std::max(1.,v[1]);worst=std::max(worst,relative);passed &= relative<1e-11;
        if(!rank)std::printf("compare level=%s field=%u role=%u component=%u abs=%.17g scale=%.17g relative=%.17g unequal=%llu\n",level,x[f].field,unsigned(x[f].role),c,v[0],std::max(1.,v[1]),relative,differences);
      }
    }
    return true;
  };
  if(!fields("current",a.fields,b.fields) || !fields("previous",a.previous_fields,b.previous_fields) || !fields("rate",a.accepted_rate_fields,b.accepted_rate_fields) || !fields("previous_rate",a.previous_rate_fields,b.previous_rate_fields))return 5;
  for(unsigned level=0;level<2;++level)for(unsigned axis=0;axis<3;++axis) {
    const auto& x=level ? a.previous_mass_flux[axis] : a.final_mass_flux[axis];const auto& y=level ? b.previous_mass_flux[axis] : b.final_mass_flux[axis];
    if(!all(x.size()==y.size()))return 6;
    double v[2]{};
    for(std::size_t i=0;i<x.size();++i){v[0]=std::max(v[0],std::abs(x[i]-y[i]));v[1]=std::max({v[1],std::abs(x[i]),std::abs(y[i])});if(!std::isfinite(x[i]) || !std::isfinite(y[i]))v[0]=INFINITY;++total;}
    MPI_Allreduce(MPI_IN_PLACE,v,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    const double relative=v[0]/std::max(1e-30,v[1]);maximum_flux_difference=std::max(maximum_flux_difference,v[0]);passed &= v[0]<1e-11;
    if(!rank)std::printf("compare_flux level=%u axis=%u abs=%.17g scale=%.17g relative=%.17g\n",level,axis,v[0],v[1],relative);
  }
  MPI_Allreduce(MPI_IN_PLACE,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if(!rank)std::printf("kerosene_esf_compare step=%llu dt=%.17g metadata=exact field_worst=%.17g field_limit=1e-11 flux_abs_max_kg_s=%.17g flux_abs_limit_kg_s=1e-11 values=%llu passed=%d\n",static_cast<unsigned long long>(a.step),a.dt,worst,maximum_flux_difference,total,int(passed));
  return passed ? 0 : 7;
}
