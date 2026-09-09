// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

// Deliberately not linked to Hundun: a broken unsigned minimum can swallow
// the UINT64_MAX success sentinel and invalidate rank-local failure handling.
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <limits>

int main(int argc,char** argv) {
  if(MPI_Init(&argc,&argv)!=MPI_SUCCESS) return 2;
  int rank{},size{};
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  const std::array<std::uint64_t,5U> samples{
      0U,1U,UINT64_C(1)<<63U,std::numeric_limits<std::uint64_t>::max()-1U,
      std::numeric_limits<std::uint64_t>::max()};
  bool passed=true;
  const auto check=[&](const char* operation,std::uint64_t actual,
                       std::uint64_t expected,int status) {
    const bool correct=status==MPI_SUCCESS && actual==expected;
    if(!correct) std::fprintf(stderr,
        "MPI_UNSIGNED_CONTRACT rank=%d operation=%s actual=%llu expected=%llu status=%d\n",
        rank,operation,static_cast<unsigned long long>(actual),
        static_cast<unsigned long long>(expected),status);
    passed &= correct;
  };
  for(const auto a : samples) for(const auto b : samples) {
    auto minimum=b,maximum=b;
    const int min_status=MPI_Reduce_local(&a,&minimum,1,MPI_UINT64_T,MPI_MIN);
    const int max_status=MPI_Reduce_local(&a,&maximum,1,MPI_UINT64_T,MPI_MAX);
    check("local-min",minimum,std::min(a,b),min_status);
    check("local-max",maximum,std::max(a,b),max_status);
  }
  for(int failing=0;failing<size;++failing) {
    const auto input=rank==failing ? static_cast<std::uint64_t>(rank) : samples.back();
    std::uint64_t output{};
    const int status=MPI_Allreduce(&input,&output,1,MPI_UINT64_T,MPI_MIN,MPI_COMM_WORLD);
    check("failure-sentinel-min",output,static_cast<std::uint64_t>(failing),status);
  }
  for(unsigned rotation=0U;rotation<samples.size();++rotation) {
    std::array<std::uint64_t,5U> input{},minimum{},maximum{};
    for(unsigned i=0U;i<samples.size();++i)
      input[i]=samples[(i+rotation+static_cast<unsigned>(rank))%samples.size()];
    MPI_Request requests[2U]{MPI_REQUEST_NULL,MPI_REQUEST_NULL};
    const int min_start=MPI_Iallreduce(input.data(),minimum.data(),5,MPI_UINT64_T,
        MPI_MIN,MPI_COMM_WORLD,&requests[0U]);
    const int max_start=MPI_Iallreduce(input.data(),maximum.data(),5,MPI_UINT64_T,
        MPI_MAX,MPI_COMM_WORLD,&requests[1U]);
    const int finish=MPI_Waitall(2,requests,MPI_STATUSES_IGNORE);
    for(unsigned i=0U;i<samples.size();++i) {
      auto lo=input[i],hi=input[i];
      for(int peer=0;peer<size;++peer) {
        const auto value=samples[(i+rotation+static_cast<unsigned>(peer))%samples.size()];
        lo=std::min(lo,value); hi=std::max(hi,value);
      }
      check("iallreduce-min",minimum[i],lo,min_start==MPI_SUCCESS ? finish : min_start);
      check("iallreduce-max",maximum[i],hi,max_start==MPI_SUCCESS ? finish : max_start);
    }
  }
  const int local=passed ? 1 : 0;
  int global{};
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(rank==0) std::printf("MPI_UNSIGNED_CONTRACT ranks=%d passed=%d\n",size,global);
  MPI_Finalize(); return global ? 0 : 1;
}
