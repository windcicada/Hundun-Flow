// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "models_esf_detail.hpp"
#include "../support/chemistry_test_support.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
using namespace hundun::v04;
int main(int argc,char** argv) {
 return hundun::test::run([&] {
  HUNDUN_CHECK(argc==2 || argc==3);
  const double expected_attenuation=argc==3?std::stod(argv[2]):1.;
  std::ifstream file(argv[1]);HUNDUN_CHECK(file.good());
  std::array<std::array<double,8>,2> source{};
  std::array<double,2> factors{};
  for(unsigned sample=0;sample<2;++sample) {
   esf::detail::StochasticSourceRequest q;
   std::array<double,8> values{},lower{},upper{};
   std::array<double,24> gradients{};
   HUNDUN_CHECK(file>>q.density>>q.molecular_viscosity>>q.turbulent_viscosity
       >>q.molecular_schmidt>>q.turbulent_schmidt>>q.dt);
   for(auto& v:q.wiener)HUNDUN_CHECK(file>>v);
   for(auto& v:values)HUNDUN_CHECK(file>>v);
   for(auto& v:gradients)HUNDUN_CHECK(file>>v);
   for(auto& v:lower)HUNDUN_CHECK(file>>v);
   for(auto& v:upper)HUNDUN_CHECK(file>>v);
   q.values=values.data();q.gradients=gradients.data();q.lower=lower.data();q.upper=upper.data();q.components=8;
   q.canonicalize_roundoff=true;
   const auto report=esf::detail::stochastic_source(q,source[sample].data());
   HUNDUN_CHECK(report.status==portable::Status::success);
   factors[sample]=report.attenuation;
   for(unsigned i=0;i<8;++i) {
    const double endpoint=values[i]+q.dt*source[sample][i]/q.density;
    if(endpoint<lower[i] || endpoint>upper[i])
      std::cout<<std::setprecision(17)<<"stochastic_endpoint sample="<<sample
          <<" component="<<i<<" value="<<endpoint<<" lower="<<lower[i]<<" upper="<<upper[i]<<'\n';
    HUNDUN_CHECK(endpoint>=lower[i] && endpoint<=upper[i]);
   }
   if(expected_attenuation==0.) {
    const double original=values[0];
    values[0]=upper[0]-(original-lower[0]);
    for(unsigned d=0;d<3;++d)gradients[d]=-gradients[d];
    std::array<double,8> mirrored{};
    const auto upper_bound=esf::detail::stochastic_source(q,mirrored.data());
    HUNDUN_CHECK(upper_bound.status==portable::Status::success);
    HUNDUN_CHECK(upper_bound.attenuation==0.);
    HUNDUN_CHECK((mirrored==std::array<double,8>{}));
    values[0]=lower[0];gradients[0]=1.;gradients[1]=gradients[2]=0.;
    const auto inward=esf::detail::stochastic_source(q,mirrored.data());
    HUNDUN_CHECK(inward.status==portable::Status::success);
    HUNDUN_CHECK(inward.attenuation==1. && mirrored[0]>0.);
   }
   // A resolvable outward increment still enforces the shared bound exactly.
   values[0]=0.;gradients[0]=1.;gradients[1]=gradients[2]=0.;q.wiener={-1,0,0};
   std::array<double,8> limited{};
   const auto bounded=esf::detail::stochastic_source(q,limited.data());
   HUNDUN_CHECK(bounded.status==portable::Status::success && bounded.attenuation==0.);
   HUNDUN_CHECK((limited==std::array<double,8>{}));
  }
  double worst{};
  for(unsigned i=0;i<8;++i)worst=std::max(worst,std::abs(source[0][i]-source[1][i])/
      std::max({1.,std::abs(source[0][i]),std::abs(source[1][i])}));
  std::cout<<std::setprecision(17)<<"captured_stochastic_pair source_relative="<<worst
      <<" attenuation="<<factors[0]<<','<<factors[1]<<" limit=1e-9\n";
  HUNDUN_CHECK(worst<1e-9);
  HUNDUN_CHECK(factors[0]==expected_attenuation && factors[1]==expected_attenuation);
 });
}
