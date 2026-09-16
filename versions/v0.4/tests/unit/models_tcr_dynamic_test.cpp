// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dynamic_detail.hpp"
#include "core_tcr_dynamic_history_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04;
using namespace hundun::v04::tcr::detail;
int main() {
  bool ok=true;
  const auto near=[](double a,double b){return std::abs(a-b)<2e-14*std::max(1.,std::abs(b));};
  auto c=cdphyso_species_control(.3,0,0,1e-12);
  ok &= c.available && c.state==SpeciesControlState::inactive && c.effective==1;
  c=cdphyso_species_control(.3,1,1,1e-12);
  ok &= c.available && c.state==SpeciesControlState::direct && near(c.effective,3./7.);
  c=cdphyso_species_control(.3,-1,1,1e-12);
  ok &= c.available && near(c.effective,3./7.);
  c=cdphyso_species_control(.3,2,1,1e-12);
  ok &= c.available && c.state==SpeciesControlState::projected && near(c.effective,2/(2+std::sqrt(2.)));
  c=cdphyso_species_control(.3,1,0,1e-12);
  ok &= c.neighbor_completion() && c.effective==0;
  c=cdphyso_species_control(.3,1e-15,1,1e-12);
  ok &= c.available && near(c.selected/1e-15,.3) && c.effective==1e-4;
  ok &= !cdphyso_species_control(.3,NAN,1,1e-12).available;
  const auto m=dynamic_filter_products({2,3,4,20,22,3});
  ok &= m.available && near(m.m_squared,1) && near(m.l_times_m,1);
  ok &= dynamic_cd_from_products(0,0)==2 && dynamic_cd_from_products(1,4)==1 &&
        dynamic_cd_from_products(8,1)==8 && dynamic_cd_from_products(100,1)==16;
  ok &= dynamic_group_cd({4,9,4},MixingGroup::product)==4 &&
        dynamic_group_cd({4,9,4},MixingGroup::reactant)==6;
  hundun::v04::detail::DynamicTcrHistory h,r;
  h.configure(77,2,7,8);r.configure(77,2,7,8);
  const auto initial=h.snapshot();
  const std::vector<std::uint8_t> bytes(initial.values.data,initial.values.data+initial.values.size);
  ok &= bool(h.begin(0));
  h.candidate(0)[3]=.4;
  ok &= bool(h.seal());h.discard();
  ok &= std::equal(bytes.begin(),bytes.end(),h.snapshot().values.data);
  ok &= bool(h.begin(0));h.candidate(0)[3]=.4;
  ok &= bool(h.seal());h.commit();
  const auto image=h.snapshot();
  ok &= bool(r.restore(image.values,1));r.commit();
  ok &= std::equal(image.values.data,image.values.data+image.values.size,r.snapshot().values.data);
  ok &= !r.begin(0) && !r.restore(image.values,2);
  ok &= bool(r.begin(1));r.candidate(0)[3]=NAN;
  ok &= !r.seal() && r.prepared_snapshot().values.size==0;
  if(!ok)std::cerr<<"dynamic TCR algebra/filter/history failure\n";
  return ok?0:1;
}
