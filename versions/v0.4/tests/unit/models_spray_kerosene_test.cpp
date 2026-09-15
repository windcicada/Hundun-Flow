// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "models_spray_properties_detail.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
using namespace hundun::v04::spray::detail;
int main(int argc,char** argv) {
  if(argc==2 && std::string_view(argv[1])=="--sample") {
    double t,p,reference;
    std::cout<<std::setprecision(17);
    while(std::cin>>t>>p>>reference) {
      const auto r=evaluate_kerosene_phase(t,p,reference);
      if(!r.available)return 2;
      std::cout<<t<<' '<<r.saturation_pressure_pa<<' '<<r.boiling_temperature_k
          <<' '<<r.liquid_cp_j_per_kg_k<<' '<<r.latent_heat_j_per_kg
          <<' '<<r.sensible_enthalpy_increment_j_per_kg<<'\n';
    }
    return 0;
  }
  bool passed=true;
  for(double t:{280.,350.,477.95,500.,650.}) {
    const auto r=evaluate_kerosene_phase(t,101325.,298.15);
    if(!r.available) {std::cerr<<"phase query unavailable at "<<t<<'\n';return 1;}
    const auto low=evaluate_kerosene_phase(t,50000.,298.15);
    const auto high=evaluate_kerosene_phase(t,200000.,298.15);
    const auto back=evaluate_kerosene_phase(298.15,101325.,t);
    const auto zero=evaluate_kerosene_phase(t,101325.,t);
    passed &= low.available && high.available && back.available && zero.available;
    passed &= low.boiling_temperature_k<r.boiling_temperature_k &&
              r.boiling_temperature_k<high.boiling_temperature_k &&
              r.saturation_pressure_pa<=.999*101325. &&
              std::abs(r.sensible_enthalpy_increment_j_per_kg+back.sensible_enthalpy_increment_j_per_kg)<1e-8 &&
              zero.sensible_enthalpy_increment_j_per_kg==0;
    const double h=1e-3;
    const auto a=evaluate_kerosene_phase(t+h,101325.,298.15);
    const auto b=evaluate_kerosene_phase(t-h,101325.,298.15);
    passed &= std::abs((a.sensible_enthalpy_increment_j_per_kg-b.sensible_enthalpy_increment_j_per_kg)/(2*h)-r.liquid_cp_j_per_kg_k)<1e-5;
  }
  for(double t:{0.,43.,684.26,std::numeric_limits<double>::quiet_NaN()})
    passed &= !evaluate_kerosene_phase(t,101325.,298.15).available;
  for(double p:{0.,-1.,1e30,std::numeric_limits<double>::infinity()})
    passed &= !evaluate_kerosene_phase(350.,p,298.15).available;
  passed &= !evaluate_kerosene_phase(350.,101325.,684.26).available;
  const double near_critical=std::nextafter(684.26,0.);
  const auto warm=evaluate_kerosene_phase(near_critical,101325.,350.);
  const auto cool=evaluate_kerosene_phase(350.,101325.,near_critical);
  passed &= warm.available && cool.available &&
      warm.sensible_enthalpy_increment_j_per_kg == -cool.sensible_enthalpy_increment_j_per_kg;
  std::cout<<"kerosene_phase relations_and_caloric_integral passed="<<passed<<'\n';
  return passed ? 0 : 1;
}
