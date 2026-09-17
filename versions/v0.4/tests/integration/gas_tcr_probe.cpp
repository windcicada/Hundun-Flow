// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dyn711_detail.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace hundun::v04::tcr::detail;
int main() {
  unsigned reset,step{};
  double dt,eta,epsilon,velocity,atime{},mean{};
  std::array<double,6> pdf,psr,composition,times;
  std::array<Dyn711RateState,6> states;
  Dyn711Clock clock;
  std::cout<<std::setprecision(17);
  while(std::cin>>reset>>dt>>eta>>epsilon>>velocity) {
    for(auto *row:{&pdf,&psr,&composition})for(double &v:*row)
      if(!(std::cin>>v))return 1;
    if(reset){step=0;atime=mean=0;states={};clock={};}
    // Uniform fixture: these are the source's complete statistics outputs
    // for time scales. Production spatial statistics require their own audit.
    mean=(atime*mean+velocity*dt)/(atime+dt);
    const double flow=std::sqrt((1/(std::abs(mean)+1e-30))*
        std::sqrt((1e-4+1e-30)/(epsilon+1e-30)));
    for(unsigned s=0;s<6;++s) {
      const double raw=std::abs(psr[s])<1e-30 || composition[s]<1e-30 ? -1 :
          composition[s]/std::abs(psr[s]);
      times[s]=raw<0 ? std::min(1.,std::max(1e-12,raw)) :
          psr[s]>0 ? .5*composition[s]/psr[s] : raw;
      Dyn711RateState next;
      if(!dyn711_advance_rate(states[s],clock,dt,pdf[s],psr[s],eta,times[s],flow,1e-30,next))return 2;
      states[s]=next;
    }
    clock=dyn711_tick(clock).next;
    ++step;atime+=dt;
    std::cout<<step<<' '<<clock.rate_intervals<<' '<<clock.rates_initialized<<' '<<clock.cphi_count;
    for(const auto &s:states)std::cout<<' '<<s.selected;
    for(const auto &s:states)std::cout<<' '<<s.pdf_sum;
    for(const auto &s:states)std::cout<<' '<<s.psr_sum;
    for(double t:times)std::cout<<' '<<t;
    std::cout<<' '<<flow<<' '<<2.<<'\n';
  }
}
