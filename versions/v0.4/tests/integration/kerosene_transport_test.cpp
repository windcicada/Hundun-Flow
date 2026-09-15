// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_app.hpp"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <array>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>
using namespace hundun::v04;
int main(int argc,char** argv) {
 MPI_Init(&argc,&argv);struct End{~End(){MPI_Finalize();}} end;
 if(argc!=3)return 2;
 std::ifstream input(argv[2]);if(!input)return 2;
 ValidatedModel model;auto s=CaseCompiler::load_and_compile(MPI_COMM_WORLD,argv[1],model);
 ThermophysicalSpec spec;ThermodynamicsPlan thermo;TransportPlan transport;
 if(s)s=ThermophysicalCompiler::load_and_compile(MPI_COMM_WORLD,model,spec,thermo,transport);
 if(!s){std::cerr<<"compile "<<unsigned(s.code)<<'/'<<s.detail<<'\n';return 3;}
 unsigned count{};input>>count;if(count!=110)return 2;double worst{};std::cout<<std::setprecision(17);
 while(count--) {
  double t,y[7],reference{};input>>t;for(auto& x:y)input>>x;input>>reference;
  if(!input || !std::isfinite(reference) || reference<=0)return 2;
  const std::array<std::string,7> names{"H2","H2O","CO","CO2","O2","N2","C12H23"};
  std::vector<double> q;
  for(const auto& scalar:model.transported_scalars) {
   auto found=std::find(names.begin(),names.end(),scalar.stable_name);
   if(found==names.end())return 5;
   q.push_back(y[std::size_t(found-names.begin())]);
  }
  MolecularTransportState v;s=transport.evaluate(t,{q.data(),q.size()},v);
  if(!s){std::cerr<<"query "<<unsigned(s.code)<<'/'<<s.detail<<'\n';return 4;}
  const double error=std::abs(v.viscosity-reference)/std::max(v.viscosity,reference);
  if(!std::isfinite(error))return 6;
  worst=std::max(worst,error);
 }
 int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 if(!rank)std::cout<<"kerosene transport cases=110 max_relative="<<worst<<" limit=1e-11\n";
 return worst<1e-11 ? 0 : 7;
}
