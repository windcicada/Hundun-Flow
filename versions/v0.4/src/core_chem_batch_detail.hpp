// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_portable.hpp"
#include "hundun/v04_status.hpp"
#include <mpi.h>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <limits>
#include <vector>
namespace hundun::v04::detail {
// Owner order is (local cell, stochastic field). MPI transports only immutable
// PH/Y tasks and results; field state and source transactions stay at the owner.
class ChemistryBatch final : public portable::GasAdvanceProvider {
public:
  explicit ChemistryBatch(portable::GasAdvanceProvider& provider,std::size_t budget=SIZE_MAX):provider_(provider),budget_(budget) {}
  const portable::GasIdentity& gas_identity() const noexcept override {return provider_.gas_identity();}
  void reset(std::size_t species, portable::Revision revision,double time,double dt) {
    ns_=species;revision_=revision;time_=time;dt_=dt;cursor_=0;
    input_.clear();hot_.clear();output_.clear();samples_.clear();
  }
  void add(double p,double h,const double* y,bool hot,const portable::GasSample* sample=nullptr) {
    if((hot_.size()+1)*(ns_+2)*sizeof(double)>budget_/4)throw std::bad_alloc();
    input_.push_back(p);input_.push_back(h);
    input_.insert(input_.end(),y,y+ns_);hot_.push_back(hot);samples_.push_back(sample ? *sample : portable::GasSample{});
  }
  Status execute(MPI_Comm comm) noexcept {
    int rank{},size{};MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const auto begin=MPI_Wtime();
    int bad=0;std::size_t selected=0;unsigned long long steps=0;
    const auto agree=[&]() {int all=0;if(MPI_Allreduce(&bad,&all,1,MPI_INT,MPI_MAX,comm)!=MPI_SUCCESS)return false;bad=all;return !bad;};
    const std::size_t iw=ns_+2,ow=2*ns_+10;
    try {
      selected=std::count(hot_.begin(),hot_.end(),true);
      if(input_.size()>INT_MAX || hot_.size()>std::size_t(INT_MAX)/ow)bad=1;
      if((input_.size()+hot_.size()*ow+selected*(iw+ow))*sizeof(double)+samples_.size()*sizeof(portable::GasSample)>budget_)bad=1;
      if(!bad)output_.resize(hot_.size()*ow);
      counts.assign(size,0);received.resize(size);displacements.resize(size);offsets.resize(size);position.resize(size);totals.resize(size);
    } catch(...) {bad=1;}
    if(!agree())return failure();
    unsigned long long local=selected;
    if(MPI_Allgather(&local,1,MPI_UNSIGNED_LONG_LONG,totals.data(),1,MPI_UNSIGNED_LONG_LONG,comm)!=MPI_SUCCESS)return failure();
    unsigned long long prefix=0,total=0;
    for(int r=0;r<size;++r){if(r<rank)prefix+=totals[r];total+=totals[r];}
    std::size_t task=0;
    for(std::size_t i=0;i<hot_.size();++i)if(hot_[i])++counts[(prefix+task++)%size];
    try {
      int offset=0;
      for(int r=0;r<size;++r){displacements[r]=offset;offset+=counts[r];}
      position=displacements;indices.resize(selected);send.resize(selected*iw);
      task=0;
      for(std::size_t i=0;i<hot_.size();++i)if(hot_[i]) {
        const auto slot=position[(prefix+task++)%size]++;indices[slot]=i;
        std::copy_n(input_.data()+i*iw,iw,send.data()+slot*iw);
      }
    } catch(...) {bad=1;}
    if(!agree())return failure();
    if(MPI_Alltoall(counts.data(),1,MPI_INT,received.data(),1,MPI_INT,comm)!=MPI_SUCCESS)return failure();
    int number=0;
    try {
      for(int r=0;r<size;++r){offsets[r]=number;if(received[r]>INT_MAX-number){bad=1;break;}number+=received[r];}
      if(std::size_t(number)>std::size_t(INT_MAX)/ow)bad=1;
      if((input_.size()+output_.size()+send.size()+selected*ow+std::size_t(number)*(iw+ow))*sizeof(double)+samples_.size()*sizeof(portable::GasSample)>budget_)bad=1;
      if(!bad){recv.resize(std::size_t(number)*iw);result.resize(std::size_t(number)*ow);returned.resize(selected*ow);}
    } catch(...) {bad=1;}
    if(!agree())return failure();
    const auto scale=[&](int width) {for(int r=0;r<size;++r){counts[r]*=width;received[r]*=width;displacements[r]*=width;offsets[r]*=width;}};
    scale(iw);
    if(MPI_Alltoallv(send.data(),counts.data(),displacements.data(),MPI_DOUBLE,recv.data(),received.data(),offsets.data(),MPI_DOUBLE,comm)!=MPI_SUCCESS)return failure();
    for(int r=0;r<size;++r){counts[r]/=iw;received[r]/=iw;displacements[r]/=iw;offsets[r]/=iw;}
    const auto integration_begin=MPI_Wtime();
    for(int i=0;i<number;++i)if(!evaluate(recv.data()+i*iw,dt_,result.data()+i*ow,steps))bad=1;
    const double integration=MPI_Wtime()-integration_begin;
    // Non-reacting cells retain their exact composition and enthalpy. A zero
    // interval still obtains a thermodynamically closed sample for the ledger.
    for(std::size_t i=0;i<hot_.size();++i)if(!hot_[i]) {
      if(samples_[i].density_kg_per_m3>0) {
        double* data=output_.data()+i*ow;const auto& s=samples_[i];
        std::copy_n(input_.data()+i*iw+2,ns_,data);std::fill_n(data+ns_,ns_,0.);
        data+=2*ns_;data[0]=s.pressure_pa;data[1]=s.temperature_k;data[2]=s.density_kg_per_m3;
        data[3]=input_[i*iw+1];data[4]=s.cp_j_per_kg_k;data[5]=s.viscosity_pa_s;data[6]=s.conductivity_w_per_m_k;
        data[7]=data[8]=data[9]=0.;
      } else if(!evaluate(input_.data()+i*iw,0.,output_.data()+i*ow,steps))bad=1;
    }
    if(!agree())return failure();
    scale(ow);
    if(MPI_Alltoallv(result.data(),received.data(),offsets.data(),MPI_DOUBLE,returned.data(),counts.data(),displacements.data(),MPI_DOUBLE,comm)!=MPI_SUCCESS)return failure();
    for(std::size_t i=0;i<selected;++i)std::copy_n(returned.data()+i*ow,ow,output_.data()+indices[i]*ow);
    double local_times[2]{integration,MPI_Wtime()-begin},maximum[2]{},sum[2]{};
    MPI_Reduce(local_times,maximum,2,MPI_DOUBLE,MPI_MAX,0,comm);
    MPI_Reduce(local_times,sum,2,MPI_DOUBLE,MPI_SUM,0,comm);
    unsigned long long global_steps=0;MPI_Reduce(&steps,&global_steps,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,comm);
    if(rank==0)std::fprintf(stdout,"chemistry_batch tasks=%llu internal_steps=%llu integration_max=%.9g integration_mean=%.9g total_max=%.9g ranks=%d\n",total,global_steps,maximum[0],sum[0]/size,maximum[1],size);
    return {};
  }
  portable::Status advance_gas(const portable::GasAdvanceQuery& q,portable::GasAdvanceOutput& out) noexcept override {
    if(!q.state.mass_fractions || q.state.coordinates!=portable::GasStateCoordinates::pressure_enthalpy || q.state.composition_fingerprint!=gas_identity().composition_fingerprint || cursor_>=hot_.size() || q.state.revision!=revision_ || q.start_time_s!=time_ || q.duration_s!=dt_ || q.state.species_count!=ns_ || out.capacity<ns_ || !out.final_mass_fractions || !out.integrated_species_density_delta_kg_per_m3)return portable::Status::invalid_input;
    const double* input=input_.data()+cursor_*(ns_+2);
    if(q.state.pressure_pa!=input[0] || q.state.enthalpy_j_per_kg!=input[1] || !std::equal(input+2,input+2+ns_,q.state.mass_fractions))return portable::Status::stale_revision;
    const double* data=output_.data()+cursor_++*(2*ns_+10);
    std::copy_n(data,ns_,out.final_mass_fractions);std::copy_n(data+ns_,ns_,out.integrated_species_density_delta_kg_per_m3);
    data+=2*ns_;
    out.final_sample={revision_,gas_identity().composition_fingerprint,data[0],data[1],data[2],data[3],data[4],data[5],data[6]};
    out.internal_step_count=static_cast<std::uint32_t>(data[7]);out.integrated_heat_release_j_per_m3=data[8];out.completed_duration_s=dt_;
    return portable::Status::success;
  }
  std::size_t owned_bytes() const noexcept {
    return sizeof(*this)+(input_.capacity()+output_.capacity()+send.capacity()+recv.capacity()+result.capacity()+returned.capacity())*sizeof(double)
      +samples_.capacity()*sizeof(portable::GasSample)+(hot_.capacity()+7)/8
      +(counts.capacity()+displacements.capacity()+received.capacity()+offsets.capacity()+position.capacity())*sizeof(int)
      +indices.capacity()*sizeof(std::size_t)+totals.capacity()*sizeof(unsigned long long);
  }
private:
    std::vector<int> counts,displacements,received,offsets,position;
    std::vector<double> send,recv,result,returned;
    std::vector<std::size_t> indices;
    std::vector<unsigned long long> totals;
  bool evaluate(const double* input,double duration,double* result,unsigned long long& steps) noexcept {
    portable::GasAdvanceQuery request{{revision_,gas_identity().composition_fingerprint,portable::GasStateCoordinates::pressure_enthalpy,input[0],input[1],0,input+2,ns_},time_,duration};
    portable::GasAdvanceOutput out{{},result,result+ns_,ns_};
    if(provider_.advance_gas(request,out)!=portable::Status::success)return false;
    const auto& s=out.final_sample;double* tail=result+2*ns_;
    tail[0]=s.pressure_pa;tail[1]=s.temperature_k;tail[2]=s.density_kg_per_m3;tail[3]=s.enthalpy_j_per_kg;tail[4]=s.cp_j_per_kg_k;tail[5]=s.viscosity_pa_s;tail[6]=s.conductivity_w_per_m_k;tail[7]=out.internal_step_count;tail[8]=out.integrated_heat_release_j_per_m3;tail[9]=duration;
    steps+=out.internal_step_count;return true;
  }
  static Status failure() noexcept {return {StatusCode::numerical_failure,10242};}
  portable::GasAdvanceProvider& provider_;
  std::size_t budget_{},ns_{},cursor_{};portable::Revision revision_{};double time_{},dt_{};
  std::vector<portable::GasSample> samples_;
  std::vector<double> input_,output_;std::vector<bool> hot_;
};
}
