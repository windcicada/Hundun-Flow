// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#include "models_spray_sgs_step_detail.hpp"
#include "core_spray_history_detail.hpp"
#include <cmath>
#include <iostream>
#include <limits>
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

namespace {
bool expect(bool condition, const char* text) {
  if (!condition) std::cerr << "FAIL: " << text << '\n';
  return condition;
}
bool same_history(const PersistentSgsBreakupState& a,const PersistentSgsBreakupState& b) {
  std::uint64_t left[kSgsHistoryLanes]{},right[kSgsHistoryLanes]{};
  encode_sgs_state(a,left);encode_sgs_state(b,right);
  for (std::size_t i=0;i<kSgsHistoryLanes;++i) if (left[i]!=right[i]) return false;
  return true;
}
bool same(const ParcelMigrationValue& a,const ParcelMigrationValue& b) {
  const auto& x=a.parcel;const auto& y=b.parcel;
  return x.id==y.id && x.position_m==y.position_m && x.velocity_m_per_s==y.velocity_m_per_s &&
      x.droplet_mass_kg==y.droplet_mass_kg && x.droplet_diameter_m==y.droplet_diameter_m &&
      x.multiplicity==y.multiplicity && x.temperature_k==y.temperature_k && x.age_s==y.age_s &&
      x.owner_global_cell==y.owner_global_cell && x.liquid_material_fingerprint==y.liquid_material_fingerprint &&
      a.tab_deformation==b.tab_deformation && a.tab_deformation_rate_per_s==b.tab_deformation_rate_per_s &&
      a.breakup_ordinal==b.breakup_ordinal && same_history(a.sgs,b.sgs);
}
SgsParcelStepInput request() {
  SgsParcelStepInput in;
  in.revision={12U,3U,1U};in.seed=123456U;in.duration_s=.01;
  auto& p=in.retained.parcel;
  p.id={17U,41U};p.position_m={.1,.2,.3};p.velocity_m_per_s={0.,0.,0.};
  p.droplet_diameter_m=.001;p.droplet_mass_kg=750.*3.14159265358979323846/6.*1.e-9;
  p.multiplicity=37.;p.temperature_k=400.;p.age_s=.04;p.owner_global_cell=9;
  p.liquid_material_fingerprint=77;
  in.retained.sgs={1U,{2.e7,.03,0.,0.,7U}};
  auto& e=in.environment;
  e.gas_velocity_m_per_s={20.,0.,0.};e.gas_density_kg_per_m3=1.;
  e.gas_dynamic_viscosity_pa_s=2.e-5;e.dissipation_m2_per_s3=1.e8;
  e.liquid_density_kg_per_m3=750.;e.surface_tension_n_per_m=.025;
  e.liquid_absolute_enthalpy_j_per_kg=-2.e5;
  return in;
}
bool restart_and_rollback() {
  using namespace hundun::v04;
  using History=hundun::v04::detail::ProductSprayHistory;
  History history,restored;
  const Int3 cells{128,1,1};
  const MeshPatch patch{{0,0,0},cells,{1,1,1},{0,0,0}};
  if (!history.configure(1234U,patch,cells,4U,77U,{},{},1U<<20) ||
      !restored.configure(1234U,patch,cells,4U,77U,{},{},1U<<20)) return false;
  auto in=request();in.revision.accepted_step=0U;
  const auto split=advance_sgs_parcel_step(in);
  if (!split.available || !split.parent_removed ||
      !history.stage_next({split.candidates.data(),split.candidate_count},0U,{}) ||
      !history.preflight_commit()) return false;
  history.commit();
  const auto snapshot=history.snapshot();
  RestartImage image;
  image.source_format_version=5U;image.step=1U;image.backward_euler_recovery=false;
  image.cell_record_identity=snapshot.identity;image.cell_record_bytes=0U;
  image.cell_records.assign(snapshot.values.data,snapshot.values.data+snapshot.values.size);
  image.cell_record_lengths.assign(snapshot.variable_cell_bytes.data,
      snapshot.variable_cell_bytes.data+snapshot.variable_cell_bytes.size);
  if (!restored.stage_restore(image) || !restored.preflight_commit()) return false;
  restored.commit();
  const auto continuous=history.accepted_parcels();
  const auto resumed=restored.accepted_parcels();
  if (continuous.size!=2U || resumed.size!=2U) return false;
  std::array<ParcelMigrationValue,4U> next{};
  std::size_t count=0U;
  for (std::size_t i=0;i<2U;++i) {
    if (!same(continuous.data[i],resumed.data[i])) return false;
    in.revision.accepted_step=1U;in.retained=continuous.data[i];
    const auto a=advance_sgs_parcel_step(in);
    in.retained=resumed.data[i];const auto b=advance_sgs_parcel_step(in);
    if (!a.available || !b.available || a.candidate_count!=b.candidate_count) return false;
    for (std::size_t j=0;j<a.candidate_count;++j) {
      if (!same(a.candidates[j],b.candidates[j])) return false;
      next[count++]=b.candidates[j];
    }
  }
  if (!restored.stage_next({next.data(),count},1U,{}) || !restored.preflight_commit()) return false;
  restored.discard();
  const auto rolled_back=restored.snapshot();
  if (rolled_back.values.size!=image.cell_records.size() ||
      !std::equal(image.cell_records.begin(),image.cell_records.end(),rolled_back.values.data)) return false;
  if (!restored.stage_next({next.data(),count},1U,{}) || !restored.preflight_commit()) return false;
  restored.commit();
  return restored.accepted_parcels().size==count;
}
} // namespace
int main() {
  bool ok=true;
  constexpr unsigned draws=1048576U;
  std::array<unsigned,8U> counts{};
  for (unsigned i=1;i<=draws;++i) {
    const auto state=initialize_sgs_lineage({123U,i},456U,7U);
    if (!valid_sgs_state(state) || state.version!=1U) return 2;
    ++counts[state.history.poisson_multiplier];
  }
  double probability=std::exp(-1.),cdf=0.,worst_z=0.;
  for (unsigned k=0;k<8U;++k) {
    if (k==7U) probability=1.-cdf;
    const double mean=draws*probability;
    const double sigma=std::sqrt(draws*probability*(1.-probability));
    worst_z=std::max(worst_z,std::abs(counts[k]-mean)/sigma);
    ok &= expect(std::abs(counts[k]-mean)<=6.*sigma+1.,"capped Poisson(1) population");
    cdf+=probability;probability/=k+1U;
  }
  ok &= expect(initialize_sgs_lineage({},1,1).version==0,"birth requires a lineage identity");
  auto in=request();const auto before=in.retained;
  const auto split=advance_sgs_parcel_step(in);
  const auto retry=advance_sgs_parcel_step(in);
  ok &= expect(split.available && split.parent_removed && split.candidate_count==2U &&
      same(in.retained,before) && same(split.candidates[0],retry.candidates[0]) &&
      same(split.candidates[1],retry.candidates[1]),"one step yields repeatable children and preserves accepted state");
  for (const auto& child:split.candidates) {
    const auto& h=child.sgs.history;
    ok &= expect(child.sgs.version==1U && h.mean_dissipation_m2_per_s3==0. && h.dissipation_age_s==0. &&
        h.mean_rate_per_s==0. && h.rate_age_s==0. && child.breakup_ordinal==0U &&
        child.parcel.age_s==in.retained.parcel.age_s && h.poisson_multiplier<=7U,
        "daughter birth resets model exposure clocks while retaining parcel age");
  }
  auto moved=in;moved.retained.parcel.owner_global_cell=100U;
  const auto moved_result=advance_sgs_parcel_step(moved);
  ok &= expect(moved_result.available && moved_result.candidates[0].parcel.id==split.candidates[0].parcel.id &&
      same_history(moved_result.candidates[0].sgs,split.candidates[0].sgs),"ownership leaves random lineage state unchanged");
  in.retained.sgs.history.poisson_multiplier=0U;
  const auto retained=advance_sgs_parcel_step(in);
  ok &= expect(retained.available && !retained.parent_removed && retained.candidate_count==1U &&
      std::abs(retained.candidates[0].sgs.history.mean_dissipation_m2_per_s3-4.e7)<1.e-8 &&
      retained.candidates[0].sgs.history.dissipation_age_s==.04 &&
      retained.candidates[0].sgs.history.poisson_multiplier==0U &&
      retained.candidates[0].parcel.age_s==.04,"accepted history advances one interval without resampling its lifetime");
  auto child_step=request();
  child_step.retained=split.candidates[0];++child_step.revision.accepted_step;
  const auto later=advance_sgs_parcel_step(child_step);
  ok &= expect(later.available && later.evolution.candidate.dissipation_age_s==.01,
               "daughter evolution begins in the next fluid step");
  for (unsigned defect=0;defect<8U;++defect) {
    auto bad=request();
    if (defect==0U) bad.retained.sgs={};
    if (defect==1U) bad.revision.accepted_step=UINT64_MAX;
    if (defect==2U) bad.revision.input_revision=0;
    if (defect==3U) bad.duration_s=0.;
    if (defect==4U) bad.environment.dissipation_m2_per_s3=-1.;
    if (defect==5U) bad.environment.gas_velocity_m_per_s[0]=std::numeric_limits<double>::quiet_NaN();
    if (defect==6U) bad.retained.sgs.history.poisson_multiplier=8U;
    if (defect==7U) bad.retained.parcel.droplet_mass_kg*=2.;
    const auto rejected=advance_sgs_parcel_step(bad);
    ok &= expect(!rejected.available && rejected.candidate_count==0U && !rejected.parent_removed &&
        same(rejected.candidates[0],{}) && same(rejected.candidates[1],{}),"failed candidates expose an empty canonical result");
  }
  ok &= expect(restart_and_rollback(),"native restart and rollback preserve subsequent SGS candidates");
  std::cout << "sgs_step population=" << draws << " poisson_max_z=" << worst_z << " counts=";
  for (auto count:counts) std::cout << count << ',';
  std::cout << " retry=pass clocks=pass daughters=pass restart=pass rollback=pass invalid=8 passed=" << ok << '\n';
  return ok?0:1;
}
