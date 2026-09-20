// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_status.hpp"
#include "core_chem_batch_detail.hpp"
#undef NDEBUG
#include <cassert>
#include <cmath>
using namespace hundun::v04;
struct Provider : portable::GasAdvanceProvider {
  portable::GasIdentity identity; bool fail{}; unsigned calls{};
  Provider(){identity.composition_fingerprint=7;}
  const portable::GasIdentity& gas_identity()const noexcept override{return identity;}
  portable::Status advance_gas(const portable::GasAdvanceQuery& q,portable::GasAdvanceOutput& o)noexcept override {
    ++calls;if(fail)return portable::Status::provider_failure;
    const double change=q.duration_s*q.state.mass_fractions[0];
    o.final_mass_fractions[0]=q.state.mass_fractions[0]-change;
    o.final_mass_fractions[1]=q.state.mass_fractions[1]+change;
    o.integrated_species_density_delta_kg_per_m3[0]=-change;
    o.integrated_species_density_delta_kg_per_m3[1]=change;
    o.final_sample={q.state.revision,7,q.state.pressure_pa,1000.,1.,q.state.enthalpy_j_per_kg,1000.,1e-5,.1};
    o.completed_duration_s=q.duration_s;o.internal_step_count=q.duration_s>0?1:0;
    return portable::Status::success;
  }
};
int main(int argc,char** argv){
  MPI_Init(&argc,&argv);int rank{},size{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  Provider provider;detail::ChemistryBatch batch(provider);const portable::Revision revision{5,6,1};
  // All hot work starts on rank zero. Others still participate in both exchanges.
  batch.reset(2,revision,1.,.1);
  const int count=rank==0?17:0;
  for(int i=0;i<count;++i){double y[]{.1+i*.01,.9-i*.01};portable::GasSample sample{revision,7,1e5,1000.,1.,100.+i+1e-8,1000.,1e-5,.1};batch.add(1e5,100.+i,y,i%3!=0,&sample);}
  assert(batch.execute(MPI_COMM_WORLD));
  for(int i=0;i<count;++i){double y[]{.1+i*.01,.9-i*.01},end[2],delta[2];
    portable::GasAdvanceQuery q{{revision,7,portable::GasStateCoordinates::pressure_enthalpy,1e5,100.+i,0,y,2},1.,.1};
    portable::GasAdvanceOutput out{{},end,delta,2};assert(batch.advance_gas(q,out)==portable::Status::success);
    assert(std::abs(end[0]-y[0]*(i%3!=0?.9:1.))<1e-15);assert(out.completed_duration_s==.1);assert(out.final_sample.enthalpy_j_per_kg==q.state.enthalpy_j_per_kg);
  }
  if(size>1 && rank==1)assert(provider.calls>0);
  batch.reset(2,revision,1.,.1);assert(batch.execute(MPI_COMM_WORLD));
  batch.reset(2,revision,1.,.1);double y[]{.5,.5};if(rank==0)batch.add(1e5,100.,y,true);
  provider.fail=rank==0;assert(!batch.execute(MPI_COMM_WORLD));
  MPI_Finalize();
}
