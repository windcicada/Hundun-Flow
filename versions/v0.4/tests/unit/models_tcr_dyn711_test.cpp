// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dyn711_detail.hpp"
#include "core_tcr_dyn711_history_detail.hpp"
#include <cmath>
#include <iostream>
#include <limits>
using namespace hundun::v04::tcr::detail;
int main() {
  bool ok = true;
  const auto near=[](double a,double b) {
    return std::isfinite(a) && std::abs(a-b) < 3e-14*std::max(1.,std::abs(b));
  };
  const auto root=[](double eta,double ratio,double tc=1.,double tf=1.) {
    return dyn711_species_control(eta,ratio,1,tc,tf,1e-30);
  };
  auto c=root(.3,1);
  ok &= c.available && !c.upper_branch && near(c.selected,3./7.);
  c=root(.3,-2);
  ok &= c.available && c.ratio==1 && near(c.selected,3./7.);
  c=root(.3,10);
  ok &= c.available && near(c.ratio,25./21.) && near(c.selected,5./7.);
  c=root(.3,1,100);
  ok &= c.available && !c.upper_branch;
  c=root(.3,1,std::nextafter(100.,101.));
  ok &= c.available && c.upper_branch && near(c.selected,1) && near(c.effective,1);
  c=root(.8,1,101);
  ok &= c.available && c.upper_branch && near(c.selected,4) && c.effective==1;
  c=root(1,2,101);
  ok &= c.available && c.linear_endpoint && !c.upper_branch && c.selected==2 && c.effective==1;
  c=root(0,1);
  ok &= c.available && c.selected==0 && c.effective==0;
  c=root(0,1,101);
  ok &= c.available && c.upper_branch && c.selected==1;
  c=dyn711_species_control(.3,10,0,1,1,1e-30);
  ok &= c.available && c.ratio==1;
  c=dyn711_species_control(.3,10,1e-30,1,1,1e-30);
  ok &= c.available && c.ratio==1;
  for(double eta:{1e-20,1e-8,.1,.5,.9,1-1e-8})
    for(double ratio:{1e-20,.5,1.,2.,1e100})
      for(double tc:{1.,101.}) {
        c=root(eta,ratio,tc);
        const long double k=c.selected;
        const long double defect=(1-static_cast<long double>(eta))*k*k-k+eta*c.ratio;
        ok &= c.available && std::abs(defect)<1e-14L*std::max(1.L,k*k);
      }
  ok &= !root(NAN,1).available && !root(-.01,1).available && !root(1.01,1).available &&
      !root(.3,NAN).available && !root(.3,1,1,0).available && !root(.3,1,-1,1).available;

  Dyn711Clock clock;
  Dyn711RateState accepted;
  for(unsigned call=1;call<=27;++call) {
    const auto tick=dyn711_tick(clock);
    ok &= tick.available && tick.evaluate_rates==(call%9==0) &&
        tick.update_cphi==(call==1 || call%7==2);
    Dyn711RateState trial;
    // The ninth call's very large rates must remain outside the eight-rate sum.
    const double rate=call%9==0 ? 1e10 : 2;
    ok &= dyn711_advance_rate(accepted,clock,.001,rate,1,.3,1,1,1e-30,trial);
    Dyn711RateState retry;
    ok &= dyn711_advance_rate(accepted,clock,.0005,rate,1,.3,1,1,1e-30,retry);
    if(call%9) {
      ok &= near(trial.pdf_sum-accepted.pdf_sum,200) &&
          near(retry.pdf_sum-accepted.pdf_sum,100);
      if(call<9)ok &= trial.selected==.2;
    } else {
      ok &= trial.pdf_sum==0 && trial.psr_sum==0 && near(trial.selected,5./7.);
    }
    accepted=trial;clock=tick.next;
  }
  Dyn711RateState unchanged{2,1,.7,.7}, trial=unchanged;
  ok &= !dyn711_advance_rate(accepted,clock,-1,1,1,.3,1,1,1e-30,trial) &&
      trial.pdf_sum==unchanged.pdf_sum && trial.selected==unchanged.selected;
  ok &= !dyn711_tick({9,6,false}).available && !dyn711_tick({0,5,false}).available;
  ok &= !dyn711_advance_rate(accepted,clock,1e308,1e308,1,.3,1,1,1e-30,trial) &&
      trial.pdf_sum==unchanged.pdf_sum;

  const auto a=dyn711_filter_products({1,1,0,2,4,0});
  const auto b=dyn711_filter_products({1,1,0,10,1,0});
  ok &= a.available && b.available && near(dyn711_filter_ratio(
      a.m_squared+b.m_squared,a.l_times_m+b.l_times_m),104./18.);
  ok &= dyn711_filter_ratio(0,0)==-1 && std::isnan(dyn711_filter_ratio(-1,1));
  ok &= dyn711_select_cphi({3,5,7})==7 && dyn711_select_cphi({3,5,12})==5 &&
      dyn711_select_cphi({3,1.5,12})==3 && dyn711_select_cphi({-1,-1,-1})==2 &&
      dyn711_select_cphi({20,22,24})==16 && dyn711_select_cphi({1,1,1})==1;
  ok &= near(dyn711_smooth_cphi(16,{1,1.5,12,16,2,4,8}),46./5.);
  ok &= std::isnan(dyn711_select_cphi({2,NAN,2})) &&
      std::isnan(dyn711_smooth_cphi(2,{2,2,2,2,2,2,NAN}));
  std::array<DynamicFilterDonor,8> donors;
  donors.fill({2,8,.5,{1,2,3}});
  DynamicFilterMoments moments;
  ok &= dynamic_filter_moments(donors.data(),8,moments) && near(moments.delta_squared,4) &&
      near(moments.density,2e6) && near(moments.gradient_squared,14) && near(moments.scalar,.5) &&
      near(moments.density_delta_squared_gradient_squared,112e6);
  const auto good=moments;
  donors[7].volume=-1;
  ok &= !dynamic_filter_moments(donors.data(),8,moments) && moments.density==good.density &&
      !dynamic_filter_moments(nullptr,8,moments) && !dynamic_filter_moments(donors.data(),0,moments);

  using hundun::v04::detail::Dyn711History;
  Dyn711History history, restored;
  history.configure(73,2,6,2);
  restored.configure(73,2,6,2);
  const auto image=[](const Dyn711History &h) {
    const auto b=h.snapshot().values;
    return std::vector<std::uint8_t>(b.data,b.data+b.size);
  };
  const auto stage=[&](Dyn711History &h,std::uint64_t step,double dt) {
    ok &= bool(h.begin(step));
    for(unsigned cell=0;cell<2;++cell) {
      if(cell==1) { ok &= bool(h.stage_inactive(cell)); continue; }
      for(unsigned q=0;q<6;++q)
        ok &= bool(h.stage_rate(cell,q,dt,2,1,.3,1,1,1e-30));
      if(dyn711_tick(h.clock()).update_cphi)
        ok &= bool(h.stage_cphi(cell,3+step%5));
    }
    ok &= bool(h.seal());
  };
  for(unsigned step=0;step<20;++step) {
    const auto before=image(history);
    stage(history,step,.001);
    history.discard();
    ok &= image(history)==before && history.prepared_snapshot().values.size==0;
    // A retry recomputes from the accepted window with its own dt.
    stage(history,step,.0005);
    history.commit();
    ok &= bool(restored.restore(history.snapshot().values,step+1));
    restored.commit();
    ok &= image(history)==image(restored) && restored.statistics_calls()==step+1;
    ok &= history.accepted(1,0).pdf_sum==0 && history.cphi(1)==2;
    if(step==7)ok &= near(history.accepted(0,0).pdf_sum,800);
    if(step==8)ok &= history.accepted(0,0).pdf_sum==0 &&
        near(history.accepted(0,0).selected,5./7.);
  }
  const auto before=image(history);
  ok &= bool(history.begin(20)) && !history.seal() &&
      history.prepared_snapshot().values.size==0;
  history.discard();
  ok &= image(history)==before;
  // A generic redistribution moves whole global-cell records, including the
  // partial window. Verify a single-cell receiving owner and corrupt clocks.
  Dyn711History one;
  one.configure(73,1,6,2);
  const auto record=history.snapshot();
  ok &= one.snapshot().identity==record.identity && bool(one.restore(
      {record.values.data,record.record_bytes},20));
  one.commit();
  ok &= one.accepted(0,0).pdf_sum==history.accepted(0,0).pdf_sum &&
      one.cphi(0)==history.cphi(0);
  auto corrupt=before;
  corrupt[record.record_bytes+24]^=1;
  ok &= !history.restore({corrupt.data(),corrupt.size()},20) && image(history)==before &&
      history.prepared_snapshot().values.size==0;
  corrupt=before;
  corrupt[40+32]=3; // Mutually exclusive upper-root and linear-endpoint flags.
  ok &= !history.restore({corrupt.data(),corrupt.size()},20) && image(history)==before;
  if(!ok)std::cerr<<"dyn711 root/window/filter contract failure\n";
  return ok?0:1;
}
