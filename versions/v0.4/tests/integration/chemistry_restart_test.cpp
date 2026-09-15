// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_app.hpp"
#include <mpi.h>
#include <iostream>
using namespace hundun::v04;
int main(int argc,char** argv) {
 MPI_Init(&argc,&argv);
 struct Lifetime {~Lifetime(){MPI_Finalize();}} lifetime;
 if(argc!=4)return 2;
 ValidatedModel source,target;
 auto status=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[1],source);
 if(status)status=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[2],target);
 CompiledCasePlan plan;
 if(status)status=ProductCompiler::compile_chemistry_restart(MPI_COMM_WORLD,source,argv[1],target,argv[2],plan);
 ProductDriver driver;if(status)status=ProductDriver::create(MPI_COMM_WORLD,std::move(plan),driver);
 RestartExpected expected;
 if(status)status=driver.restart_expected(expected,RestartStorageCompatibility::strict,RestartHistoryPolicy::refine_chemistry);
 RestartImage image;
 if(status)status=RestartReader::load(MPI_COMM_WORLD,argv[3],expected,image);
 if(status)status=driver.initialize_restart(image,RestartStorageCompatibility::strict,RestartHistoryPolicy::refine_chemistry);
 RestartSnapshot snapshot;if(status)status=driver.committed_restart_snapshot(snapshot);
 int ok=bool(status);MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
 if(!ok){std::cerr<<"restore "<<unsigned(status.code)<<'/'<<status.detail<<'\n';return 3;}
 ok=snapshot.plan==expected.plan && snapshot.schema==image.schema && snapshot.geometry==image.geometry &&
    snapshot.time==image.time && snapshot.dt==image.dt && snapshot.step==image.step &&
    snapshot.controller_state==image.controller_state && snapshot.pressure_reference==image.pressure_reference &&
    snapshot.previous_pressure_reference==image.previous_pressure_reference && snapshot.closed_mass_target==image.closed_mass_target &&
    snapshot.method_history_signature==image.method_history_signature &&
    snapshot.final_mass_flux.revision==image.final_mass_flux_revision &&
    snapshot.previous_mass_flux.revision==image.previous_mass_flux_revision;
 unsigned long long values=0;
 const auto fields=[&](Span<const RestartFieldView> a,const std::vector<RestartImageField>& b) {
  if(a.size!=b.size())return false;
  for(std::size_t i=0;i<a.size;++i) {
   const auto v=a.data[i].values;
   if(a.data[i].role!=b[i].role || v.field!=b[i].field || v.components!=b[i].components)return false;
   std::size_t k=0;
   for(int z=0;z<v.interior.z;++z)for(int y=0;y<v.interior.y;++y)for(int x=0;x<v.interior.x;++x)
    for(unsigned c=0;c<v.components;++c) {
     if(k>=b[i].values.size() || v.unchecked({x,y,z},c)!=b[i].values[k++])return false;
     ++values;
    }
   if(k!=b[i].values.size())return false;
  }
  return true;
 };
 ok &= fields(snapshot.fields,image.fields) && fields(snapshot.previous_fields,image.previous_fields) &&
       fields(snapshot.accepted_rate_fields,image.accepted_rate_fields) && fields(snapshot.previous_rate_fields,image.previous_rate_fields);
 for(unsigned level=0;level<2;++level) {
  const auto flux=level ? snapshot.previous_mass_flux : snapshot.final_mass_flux;
  unsigned axis=0;
  for(const auto face:{flux.x,flux.y,flux.z}) {
   const auto& raw=level ? image.previous_mass_flux[axis++] : image.final_mass_flux[axis++];std::size_t k=0;
   for(int z=0;z<face.extents.z;++z)for(int y=0;y<face.extents.y;++y)for(int x=0;x<face.extents.x;++x) {
    ok &= k<raw.size() && face.unchecked({x,y,z})==raw[k++];++values;
   }
   ok &= k==raw.size();
  }
 }
 const auto records=snapshot.cell_records;
 ok &= records.identity==image.cell_record_identity && records.record_bytes==image.cell_record_bytes &&
       records.values.size==image.cell_records.size() && records.variable_cell_bytes.size==image.cell_record_lengths.size();
 if(ok) {
  for(std::size_t i=0;i<records.values.size;++i)ok &= records.values.data[i]==image.cell_records[i];
  for(std::size_t i=0;i<records.variable_cell_bytes.size;++i)
   ok &= records.variable_cell_bytes.data[i]==image.cell_record_lengths[i];
 }
 MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
 MPI_Allreduce(MPI_IN_PLACE,&values,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
 int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 if(!rank)std::cout<<"chemistry_refinement exact="<<ok<<" values="<<values<<'\n';
 return ok ? 0 : 4;
}
