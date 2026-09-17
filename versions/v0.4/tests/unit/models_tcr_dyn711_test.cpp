// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dyn711_detail.hpp"
#include "core_tcr_dyn711_history_detail.hpp"
#include "core_tcr_history_detail.hpp"
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
      !root(.3,NAN).available && !root(.3,1,1,-1).available && !root(.3,1,-1,1).available;
  // Local dimensional times, including their laminar limits, feed the actual
  // upper/lower root switch. Scaling all volumetric material quantities leaves
  // the specific SGS state and both times invariant.
  const auto times=dyn711_flow_times(2,8,4,.02);
  ok &= times.available && near(times.integral,1) && near(times.kolmogorov,.05) &&
      near(times.flow,std::sqrt(.05));
  for(double scale:{.01,2.,1e6}) {
    const auto scaled=dyn711_flow_times(2,8*scale,4*scale,.02*scale);
    ok &= scaled.available && near(scaled.integral,times.integral) &&
        near(scaled.kolmogorov,times.kolmogorov) && near(scaled.flow,times.flow);
  }
  const auto still=dyn711_flow_times(0,0,1,1e-5);
  const auto shear=dyn711_flow_times(0,1,1,1e-5);
  ok &= still.available && std::isinf(still.flow) && shear.available && shear.flow==0;
  c=root(.3,1,1,still.flow);
  ok &= c.available && !c.upper_branch && near(c.selected,3./7.);
  c=root(.3,1,1,shear.flow);
  ok &= c.available && c.upper_branch && c.selected==1;
  c=root(.3,1,0,0);
  ok &= c.available && !c.upper_branch;
  ok &= !dyn711_flow_times(-1,1,1,1).available &&
      !dyn711_flow_times(1,-1,1,1).available &&
      !dyn711_flow_times(1,1,0,1).available &&
      !dyn711_flow_times(1,1,1,0).available &&
      !dyn711_flow_times(1,INFINITY,1,1).available &&
      !dyn711_flow_times(NAN,1,1,1).available && !root(.3,1,1,NAN).available;
  const auto wide=dyn711_flow_times(1e200,1e100,1e200,1e100);
  ok &= wide.available && near(wide.integral/1e300,1) &&
      near(wide.flow/1e150,1);

  // Same step-eight checkpoint, two/four ranks: roundoff-sized transverse
  // gradients straddle the raw clock's ninth-call branch threshold while
  // both momentum coefficients equal the molecular viscosity in FP64.
  constexpr double rho=.29187575647980235, mu=.000042722622661398233;
  constexpr double eps=.000011682891580905665, chemical=.00055889988844045094;
  const double energies[]{6.0303835783520853e-16,6.610474907196156e-16};
  const double viscosities[]{4.3322643819689027e-22,4.972176509199209e-22};
  for(unsigned partition=0;partition<2;++partition) {
    const auto raw=dyn711_flow_times(energies[partition],eps,rho,mu);
    const auto resolved=dyn711_resolved_flow_times(energies[partition],eps,rho,mu,
        viscosities[partition]);
    const auto select=[&](double flow) {
      return dyn711_species_control(.27057497707038558,-5.6204569200725123,
          -5.6204403098872344,chemical,flow,1e-30);
    };
    const auto before=select(raw.flow),after=select(resolved.flow);
    ok &= raw.available && before.available && before.upper_branch==(partition==0) &&
        mu+rho*viscosities[partition]==mu && resolved.available && resolved.flow==0 &&
        after.available && after.upper_branch && near(after.selected,.99999825730176661);
  }
  const auto resolved=dyn711_resolved_flow_times(2,8,4,.02,.001);
  const auto quiet=dyn711_resolved_flow_times(0,0,1,1e-5,0);
  ok &= resolved.available && resolved.flow==times.flow &&
      quiet.available && std::isinf(quiet.flow) &&
      !dyn711_resolved_flow_times(-1,1,1,1,0).available &&
      !dyn711_resolved_flow_times(1,1,1,1,-1).available &&
      !dyn711_resolved_flow_times(1,1,1,1,NAN).available &&
      !dyn711_resolved_flow_times(1,1,1,1,INFINITY).available &&
      !dyn711_resolved_flow_times(1,1,1e300,1,1e300).available;

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
  const auto stage=[&](Dyn711History &h,std::uint64_t step,double dt, bool seal_history=true) {
    ok &= bool(h.begin(step));
    for(unsigned cell=0;cell<2;++cell) {
      if(cell==1) { ok &= bool(h.stage_inactive(cell)); continue; }
      for(unsigned q=0;q<6;++q)
        ok &= bool(h.stage_rate(cell,q,dt,2,1,.3,1,1,1e-30));
      if(dyn711_tick(h.clock()).update_cphi)
        ok &= bool(h.stage_cphi(cell,3+step%5));
    }
    if(seal_history)ok &= bool(h.seal());
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
  // The native owner uses model-specific record identity and delegates the
  // full transaction, including fixed V4 and enclosing V5 byte contracts.
  using hundun::v04::detail::ProductTcrHistory;
  ProductTcrHistory owner, receiver, other;
  owner.configure_dyn711(73,2,6,2);
  receiver.configure_dyn711(73,2,6,2);
  other.configure_dynamic(73,2,6,2);
  ok &= owner.enabled() && owner.dyn711() && !owner.dynamic() &&
      owner.snapshot().identity != other.snapshot().identity &&
      owner.owned_bytes() >= owner.dyn711()->owned_bytes();
  hundun::v04::RestartImage native;
  native.source_format_version=4;native.step=20;native.backward_euler_recovery=false;
  native.cell_record_identity=record.identity;native.cell_record_bytes=record.record_bytes;
  native.cell_records=before;
  ok &= bool(owner.stage_restore(native)) && !other.stage_restore(native);
  const auto pending=owner.prepared_snapshot();
  ok &= pending.identity==record.identity && pending.record_bytes==record.record_bytes &&
      std::equal(before.begin(),before.end(),pending.values.data);
  owner.discard();
  ok &= owner.prepared_snapshot().values.size==0 && owner.dyn711()->statistics_calls()==0;
  ok &= bool(owner.stage_restore(native));owner.commit();
  ok &= owner.dyn711()->statistics_calls()==20 && image(*owner.dyn711())==before;
  // V5 combined owners have already validated outer identity and redistribute
  // complete per-cell records through this same typed restore operation.
  ok &= bool(receiver.stage_restore_records(owner.snapshot().values,20));receiver.commit();
  ok &= image(*receiver.dyn711())==before;
  native.source_format_version=6;
  native.backward_euler_recovery=true;
  ok &= bool(receiver.stage_restore(native));receiver.commit();
  ok &= image(*receiver.dyn711())==before && receiver.dyn711()->statistics_calls()==20;
  native.backward_euler_recovery=false;
  ok &= !receiver.stage_restore(native) && image(*receiver.dyn711())==before;
  native.backward_euler_recovery=true;
  native.step=21;
  ok &= !receiver.stage_restore(native) && image(*receiver.dyn711())==before;
  native.step=20;
  stage(*owner.dyn711(),20,.0005,false);
  ok &= bool(owner.seal());
  owner.commit();
  ok &= owner.dyn711()->statistics_calls()==21 && owner.prepared_snapshot().values.size==0;
  native.cell_records[24]^=1;
  ok &= !receiver.stage_restore(native) && image(*receiver.dyn711())==before &&
      receiver.prepared_snapshot().values.size==0;
  if(!ok)std::cerr<<"dyn711 root/window/filter contract failure\n";
  return ok?0:1;
}
