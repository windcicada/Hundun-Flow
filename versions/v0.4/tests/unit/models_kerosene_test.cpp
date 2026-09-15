// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "models_kerosene_detail.hpp"
#include "../support/chemistry_test_support.hpp"
#include <fstream>
#include <limits>
using Kernel=hundun::v04::chemistry::detail::KeroseneFourStep;
int main(int argc,char**argv) {
 return hundun::test::run([&] {
  HUNDUN_CHECK(argc==2);std::ifstream input(argv[1]);HUNDUN_CHECK(input.good());
  unsigned count=0;double worst=0,t,rho,gascon,mw;
  Kernel previous;Kernel::Tuple previous_trial{},previous_rates{};double previous_rho=0;
  while(input>>t>>rho>>gascon>>mw) {
   Kernel::Tuple accepted{},g{},trial{},expected{},rates{};
   for(auto* tuple:{&accepted,&g,&trial,&expected})for(auto&v:*tuple)HUNDUN_CHECK(input>>v);
   Kernel kernel;rates.fill(71.);const auto untouched=rates;
   HUNDUN_CHECK(!kernel.evaluate(rho,trial,rates));HUNDUN_CHECK(rates==untouched);
   HUNDUN_CHECK(kernel.prepare(t,gascon,mw,accepted,g));
   HUNDUN_CHECK(kernel.evaluate(rho,trial,rates));
   // Interleaved cell lanes keep independent prepared coefficients.
   if(count) {
    Kernel::Tuple earlier{};
    HUNDUN_CHECK(previous.evaluate(previous_rho,previous_trial,earlier));
    HUNDUN_CHECK(earlier==previous_rates);
   }
   previous=kernel;previous_trial=trial;previous_rates=rates;previous_rho=rho;
   for(unsigned i=0;i<7;++i)worst=std::max(worst,std::abs(rates[i]-expected[i])/
       std::max({1.,std::abs(rates[i]),std::abs(expected[i])}));
   const auto valid_rates=rates;
   auto bad=trial;bad[4]=std::numeric_limits<double>::quiet_NaN();
   HUNDUN_CHECK(!kernel.evaluate(rho,bad,rates));HUNDUN_CHECK(rates==valid_rates);
   HUNDUN_CHECK(!kernel.prepare(0,gascon,mw,accepted,g));
   HUNDUN_CHECK(kernel.evaluate(rho,trial,rates));HUNDUN_CHECK(rates==valid_rates);
   // Molar element ledgers close also for signed internal Newton trials.
   const double atoms[4][7]={{0,1,1,2,2,0,0},{2,2,0,0,0,0,23},
                            {0,0,1,1,0,0,12},{0,0,0,0,0,2,0}};
   for(const auto& element:atoms) {
    long double sum=0,scale=1;
    for(unsigned i=0;i<7;++i){const long double term=static_cast<long double>(rates[i])*element[i];sum+=term;scale+=std::abs(term);}
    HUNDUN_CHECK(std::abs(sum)/scale<1e-11);
   }
   ++count;
  }
  HUNDUN_CHECK(input.eof());HUNDUN_CHECK(count==192);HUNDUN_CHECK(worst<1e-11);
  std::cout<<"kerosene_four_step states="<<count<<" maximum_relative="<<worst<<'\n';
 });
}
