// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dyn711_detail.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>
using namespace hundun::v04::tcr::detail;
namespace {
int mix_probe(bool cf) {
  constexpr std::array<std::array<int,3>,8> offsets{{
      {{0,0,0}},{{-1,0,0}},{{0,-1,0}},{{0,0,-1}},
      {{0,-1,-1}},{{-1,0,-1}},{{-1,-1,0}},{{-1,-1,-1}}}};
  const auto index=[](int x,int y,int z){return x+4*y+16*z;};
  std::array<double,64> rho,volume,scalar;
  std::cout<<std::setprecision(17);
  while(std::cin>>rho[0]) {
    for(unsigned i=1;i<64;++i)if(!(std::cin>>rho[i]))return 1;
    for(auto *a:{&volume,&scalar})for(double &v:*a)if(!(std::cin>>v))return 1;
    std::array<DynamicFilterDonor,64> grid;
    std::array<DynamicFilterProducts,64> products;
    std::array<double,64> cphi,ratio;
    cphi.fill(2.);
    for(int z=0;z<4;++z)for(int y=0;y<4;++y)for(int x=0;x<4;++x) {
      const auto i=index(x,y,z);
      auto &d=grid[i];d.density=rho[i];d.volume=volume[i];d.scalar=scalar[i];
      const std::array<int,3> p{x,y,z};
      for(unsigned axis=0;axis<3;++axis) {
        auto a=p,b=p;a[axis]=std::max(0,p[axis]-1);b[axis]=std::min(3,p[axis]+1);
        d.gradient[axis]=(scalar[index(b[0],b[1],b[2])]-scalar[index(a[0],a[1],a[2])])/
            (b[axis]-a[axis]);
      }
    }
    for(int z=1;z<3;++z)for(int y=1;y<3;++y)for(int x=1;x<3;++x) {
      std::array<DynamicFilterDonor,8> donors;
      for(unsigned j=0;j<8;++j) {
        const auto &o=offsets[j];donors[j]=grid[index(x+o[0],y+o[1],z+o[2])];
      }
      DynamicFilterMoments moments;
      if(!dynamic_filter_moments(donors.data(),8,moments))return 2;
      products[index(x,y,z)]=cf ? dynamic_filter_products(moments) : dyn711_filter_products(moments);
      if(!products[index(x,y,z)].available)return 3;
    }
    for(int z=1;z<3;++z)for(int y=1;y<3;++y)for(int x=1;x<3;++x) {
      double m2{},lm{};
      for(const auto &o:offsets) {
        const auto &d=products[index(x+(cf&&x==1 ? -o[0] : o[0]),
            y+(cf&&y==1 ? -o[1] : o[1]),z+(cf&&z==1 ? -o[2] : o[2]))];
        m2+=d.m_squared;lm+=d.l_times_m;
      }
      const auto i=index(x,y,z);
      ratio[i]=dyn711_filter_ratio(m2,lm);
      cphi[i]=cf ? dynamic_cd_from_products(m2,lm) : dyn711_select_cphi({ratio[i],-1.,-1.});
    }
    for(int z=1;z<3;++z)for(int y=1;y<3;++y)for(int x=1;x<3;++x) {
      std::array<double,7> neighbors;
      for(unsigned j=1;j<8;++j) {
        const auto &o=offsets[j];neighbors[j-1]=cphi[index(x+o[0],y+o[1],z+o[2])];
      }
      const auto i=index(x,y,z);
      if(cf)std::cout<<cphi[i]<<' '<<cphi[i]<<' '<<cphi[i]<<'\n';
      else std::cout<<ratio[i]<<' '<<dyn711_smooth_cphi(cphi[i],neighbors)<<'\n';
    }
  }
  return 0;
}
}
int main(int argc,char **argv) {
  if(argc==2 && std::string_view(argv[1])=="mix")return mix_probe(false);
  if(argc==2 && std::string_view(argv[1])=="cf_mix")return mix_probe(true);
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
