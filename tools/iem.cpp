// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "models_esf_detail.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
using namespace hundun::v04;
int main(int argc,char** argv) {
  if(argc!=3)return 1;
  const double cz=std::stod(argv[1]),schmidt=std::stod(argv[2]);
  if(!std::isfinite(cz) || cz<=0 || !std::isfinite(schmidt) || schmidt<=0)return 1;
  double volume,mu,mut,cd,rho,dt,mean,old;
  combustion::ChemistryIdentity identity{1,{{1,{1},0},{1,{1},0}},42};
  esf::detail::Workspace workspace(2);
  std::cout<<std::setprecision(17);
  while(std::cin>>volume>>mu>>mut>>cd>>rho>>dt>>mean>>old) {
    // Uniform composition; enthalpy is the scalar with an arbitrary datum.
    double values[]{.5,.5,old,.5,.5,2*mean-old};
    esf::detail::Request request;
    request.accepted={{0,7,1},42,2,2,values};
    request.expected_revision=request.accepted.revision;
    request.identity=&identity;
    request.dt_s=dt;
    const double delta2=std::pow(volume,2./3.);
    // Unit-Lewis, constant-Pr=Sc material isolates the product cache scale.
    request.mixing_time_s=cz*delta2/(2*(mu+mut)/(schmidt*rho));
    request.tcr_control=std::pow(cd/2,3);
    auto current=workspace.advance(request);
    if(current.status!=portable::Status::success)return 2;
    const double first=current.candidate.values[2];
    const double factor=current.relaxation_factor;
    // Match the COAST IEM rate; retain the native exponential update.
    request.mixing_time_s=cd==0 ? 1 : rho*delta2/(cd*(mu+mut));
    request.tcr_control=cd==0 ? 0 : 1;
    auto matched=workspace.advance(request);
    if(matched.status!=portable::Status::success)return 3;
    std::cout<<first<<' '<<factor<<' '<<matched.candidate.values[2]<<' '
             <<matched.relaxation_factor<<' '<<matched.max_mean_residual<<'\n';
  }
  return std::cin.eof() ? 0 : 1;
}
