// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn
#include "models_esf_detail.hpp"
#include "hundun/v04_flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>

namespace { std::size_t allocations{}; }
void* operator new(std::size_t n) {
  ++allocations;
  if (auto* p=std::malloc(n ? n : 1)) return p;
  throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }

using namespace hundun::v04;
namespace {
bool near(double a,double b) {
  return std::isfinite(a) && std::abs(a-b)<=2e-14*std::max(1.,std::abs(b));
}
bool successful(esf::detail::DualStateReport r) {
  return r.status==portable::Status::success;
}
}
int main() {
  namespace esf=esf::detail;
  using portable::Status;
  const portable::Revision revision{32000,7,1};
  std::array<double,12> fields{{.2,.8,100,.8,.2,300,.2,.8,100,.8,.2,300}};
  const auto saved=fields;
  std::array<double,3> auxiliary{{.1,.9,175}},mean{},variance{};
  const auto auxiliary_saved=auxiliary;
  double rho[]{1,4,1,4};
  esf::DualStateRequest q{{revision,42,4,2,fields.data()},
      {revision,42,1,2,auxiliary.data()},revision,rho,2.5};
  const esf::DualStateOutput out{mean.data(),variance.data(),3};
  for (auto count:{2U,4U}) {
    q.stochastic.fields=count;
    const auto before=allocations;
    const auto r=esf::dual_state_moments(q,out);
    if (!successful(r) || allocations!=before || r.revision!=revision ||
        r.composition_fingerprint!=42 || r.model_identity!=0x4553464455410001ULL ||
        !near(mean[0],.5) || !near(mean[1],.5) || mean[2]!=200 ||
        !near(variance[0],.09) || !near(variance[1],.09) || variance[2]!=10000 ||
        !near(r.statistical_density_kg_per_m3,1.6) || r.pressure_density_kg_per_m3!=2.5 ||
        !near(r.max_species_mean_gap,.4) || r.enthalpy_mean_gap_j_per_kg!=25 ||
        fields!=saved || auxiliary!=auxiliary_saved) return 1;
    // A density-weighted composition would produce .68, while Favre
    // ensemble moments use the existing equal field weights.
    if (near(mean[0],(.2+4*.8)/5)) return 2;
  }
  const auto valid=q;
  const auto rejected=[&](Status expected,esf::DualStateOutput destination) {
    mean.fill(71);variance.fill(83);
    const auto r=esf::dual_state_moments(q,destination);
    return r.status==expected && mean==std::array<double,3>{{71,71,71}} &&
        variance==std::array<double,3>{{83,83,83}};
  };
  ++q.auxiliary.revision.input_revision;
  if (!rejected(Status::stale_revision,out)) return 3;
  q=valid;++q.stochastic.composition_fingerprint;
  if (!rejected(Status::identity_mismatch,out)) return 4;
  q=valid;
  if (!rejected(Status::capacity_exceeded,{mean.data(),variance.data(),2})) return 5;
  if (!rejected(Status::invalid_input,{fields.data(),variance.data(),3}) || fields!=saved) return 6;
  if (!rejected(Status::invalid_input,{mean.data(),auxiliary.data(),3}) || auxiliary!=auxiliary_saved) return 7;
  if (!rejected(Status::invalid_input,{mean.data(),mean.data()+1,3})) return 8;
  rho[3]=0;
  if (!rejected(Status::invalid_input,out)) return 9;
  rho[3]=4;fields[10]=-1e-5;
  if (!rejected(Status::invalid_input,out)) return 10;
  fields=saved;auxiliary[0]=-.0003146342157678451;
  if (!rejected(Status::invalid_input,out)) return 11;
  auxiliary=auxiliary_saved;fields[11]=std::numeric_limits<double>::max();
  if (!rejected(Status::invalid_input,out)) return 12;
  fields=saved;
  // A borrowed input stays valid after a read-only summary. A subsequent
  // mutation by its owner invalidates that input in both state positions.
  esf::Workspace work(2);
  double target[]{.5,.5,200};
  const auto centered=work.recenter(q.stochastic,target);
  if (centered.status!=Status::success) return 13;
  q.stochastic=centered.candidate;
  if (!successful(esf::dual_state_moments(q,out)) || !work.valid(centered.candidate)) return 14;
  work.recenter(q.stochastic,target);
  if (!rejected(Status::stale_revision,out)) return 15;
  q=valid;q.auxiliary=centered.candidate;q.auxiliary.fields=1;
  if (!rejected(Status::stale_revision,out)) return 16;
  // Scale-safe harmonic reduction also covers subnormal positive density.
  q=valid;
  const double tiny=std::numeric_limits<double>::denorm_min();
  for (auto& r:rho) r=tiny;
  auto r=esf::dual_state_moments(q,out);
  if (!successful(r) || r.statistical_density_kg_per_m3!=tiny) return 17;

  // GTMC captured cell (71,118,94), transported four-field state at the
  // original first proposal. These are full JL4 fractions and total h.
  const double gtmc[]{
    .00088776958640251992,.045964734973690566,.09337334823944235,.023503917013272743,
    .10624205320045242,.000048751866576052752,.72997942512016334,-223156.57924449528,
    .000029621985456657821,.027367642871466091,.090893525705309314,.036238929467601981,
    .11956972616088404,.00016327942730240237,.72573727438197955,-248596.70771350342,
    .0025731932567788066,.033545342864482235,.089111484350296269,.035154548271952767,
    .11443695817781384,.00041830345628565154,.72476016962239043,-254914.66433910499,
    .0000058985966336856907,.022438717600629238,.11450207028777074,.018890222502643696,
    .11730612435394319,.000082232803981034831,.72677473385439839,-241796.96075271681};
  const double field0[]{.00017642873757605626,.031198793619182291,.096918548221474166,
      .029171552185479808,.11531664209700802,.00017388886713953325,.72704414627214009,
      -240739.43221865833};
  const double field_rho[]{.17483343024895043,.16424746983496011,
      .16910760077904316,.15940007657139027};
  double gmean[8],gvariance[8];
  q={{revision,99,4,7,gtmc},{revision,99,1,7,field0},revision,field_rho,.16580727378433224};
  r=esf::dual_state_moments(q,{gmean,gvariance,8});
  long double volume=0;
  for (auto d:field_rho) volume+=1.L/d;
  if (!successful(r) || !near(gmean[0],.0008741208563179175) ||
      !near(r.statistical_density_kg_per_m3,static_cast<double>(4/volume)) ||
      !(gvariance[0]>1e-6) || near(gmean[0],field0[0])) return 18;
  // Auxiliary signed coordinates follow COAST's positive-weight EOS.
  // For ideal mixtures this PH normalization yields the identical T and
  // pressure density, including the extensive weight in the gas law.
  double signed_aux[]{-.1,1.1,200000},normalized[]{71,83};
  const esf::View raw{revision,42,1,2,signed_aux};
  const auto before=allocations;
  auto coordinates=esf::auxiliary_eos_coordinates(raw,revision,normalized,2);
  if (coordinates.status!=Status::success || allocations!=before ||
      normalized[0]!=0 || normalized[1]!=1 ||
      !near(coordinates.positive_weight_sum,1.1) ||
      !near(coordinates.added_positive_weight,.1) || signed_aux[0]!=-.1) return 19;
  const double cp=1200,hf=-50000,pressure=100000,gas_constant=300;
  const double temperature=300+(coordinates.query_enthalpy_j_per_kg-hf)/cp;
  const double direct_temperature=300+(signed_aux[2]/1.1-hf)/cp;
  const double query_density=pressure/(gas_constant*temperature);
  if (!near(temperature,direct_temperature) ||
      !near(query_density/coordinates.positive_weight_sum,
            pressure/(gas_constant*direct_temperature*1.1))) return 20;
  q=valid;q.auxiliary=raw;q.auxiliary_density_kg_per_m3=query_density/1.1;
  r=esf::dual_state_moments(q,out);
  if (!successful(r) || !near(r.max_species_mean_gap,.6) || mean[0]!=.5 ||
      signed_aux[0]!=-.1) return 21;
  normalized[0]=71;normalized[1]=83;
  auto stale=revision;++stale.input_revision;
  if (esf::auxiliary_eos_coordinates(raw,stale,normalized,2).status!=Status::stale_revision ||
      normalized[0]!=71 || normalized[1]!=83) return 22;
  if (esf::auxiliary_eos_coordinates(raw,revision,signed_aux,2).status!=Status::invalid_input ||
      signed_aux[0]!=-.1) return 23;
  if (esf::auxiliary_eos_coordinates(raw,revision,normalized,1).status!=Status::capacity_exceeded ||
      normalized[0]!=71 || normalized[1]!=83) return 24;
  signed_aux[1]=1.2;
  if (esf::auxiliary_eos_coordinates(raw,revision,normalized,2).status!=Status::invalid_input ||
      normalized[0]!=71 || normalized[1]!=83) return 25;
  portable::GasSample normalized_sample{revision,42,pressure,temperature,query_density,
      coordinates.query_enthalpy_j_per_kg,cp,1e-5,.1};
  esf::AuxiliaryPressureState tangent;
  if(esf::auxiliary_pressure_state(coordinates,normalized_sample,revision,42,tangent)!=Status::success)
    return 26;
  // Independent extensive caloric EOS: h0=M*(hf+cp*(T-Tref)).
  const double weight=1.1;
  const auto raw_density=[&](double p,double h0) {
    const double extensive_cp=weight*cp;
    const double extensive_offset=weight*(hf-cp*300);
    return p*extensive_cp/(weight*gas_constant*(h0-extensive_offset));
  };
  const double dh=10,dp=1;
  const double finite_p=(raw_density(pressure+dp,200000)-raw_density(pressure-dp,200000))/(2*dp);
  const double finite_h=(raw_density(pressure,200000+dh)-raw_density(pressure,200000-dh))/(2*dh);
  if(std::abs(tangent.density_pressure_derivative/finite_p-1)>2e-11 ||
     std::abs(tangent.density_enthalpy_derivative/finite_h-1)>1e-8 ||
     !near(pressure*tangent.density_pressure_derivative,tangent.density_kg_per_m3) ||
     !near(tangent.density_kg_per_m3,raw_density(pressure,200000))) return 27;
  const auto saved_tangent=tangent;
  normalized_sample.revision=stale;
  if(esf::auxiliary_pressure_state(coordinates,normalized_sample,revision,42,tangent)!=Status::stale_revision ||
     tangent.density_enthalpy_derivative!=saved_tangent.density_enthalpy_derivative) return 28;
  normalized_sample.revision=revision;
  if(esf::auxiliary_pressure_state(coordinates,normalized_sample,revision,43,tangent)!=Status::identity_mismatch)
    return 29;
  normalized_sample.cp_j_per_kg_k=0;
  if(esf::auxiliary_pressure_state(coordinates,normalized_sample,revision,42,tangent)!=Status::invalid_input ||
     tangent.density_kg_per_m3!=saved_tangent.density_kg_per_m3) return 30;
  normalized_sample.cp_j_per_kg_k=cp;
  normalized_sample.enthalpy_j_per_kg+=1;
  if(esf::auxiliary_pressure_state(coordinates,normalized_sample,revision,42,tangent)!=Status::invalid_input ||
     tangent.density_pressure_derivative!=saved_tangent.density_pressure_derivative) return 31;
  // A common pressure-work correction may translate physical h and raw h0
  // together while preserving their distinct base values. Form the public
  // coupled Jacobian without representing physical properties as auxiliary.
  constexpr double physical_h=175000;
  const PressureThermoState density{saved_tangent.density_kg_per_m3,
      saved_tangent.density_pressure_derivative,saved_tangent.density_enthalpy_derivative};
  PressureEnergyThermoJacobian jacobian;
  if(!form_pressure_energy_density_jacobian(pressure,physical_h,density,jacobian))return 32;
  const auto energy=[&](double p,double physical) {
    return raw_density(p,physical+25000)*physical-p;
  };
  const double energy_p=(energy(pressure+dp,physical_h)-energy(pressure-dp,physical_h))/(2*dp);
  const double energy_h=(energy(pressure,physical_h+dh)-energy(pressure,physical_h-dh))/(2*dh);
  if(std::abs(jacobian.dq_dp_hY-energy_p)>2e-11 ||
     std::abs(jacobian.dq_dh_pY-energy_h)>1e-8)return 33;
  const auto saved_jacobian=jacobian;
  auto mixed=density;mixed.drho_dp_hY*=1.01;
  if(form_pressure_energy_density_jacobian(pressure,physical_h,mixed,jacobian) ||
     jacobian.density!=saved_jacobian.density || jacobian.dq_dh_pY!=saved_jacobian.dq_dh_pY)
    return 34;
  normalized_sample.enthalpy_j_per_kg=coordinates.query_enthalpy_j_per_kg;
  auto old_coordinates=coordinates;old_coordinates.revision=stale;
  if(esf::auxiliary_pressure_state(old_coordinates,normalized_sample,revision,42,tangent)!=Status::stale_revision ||
     tangent.density_kg_per_m3!=saved_tangent.density_kg_per_m3)return 35;
  // Closed-cell mass and energy correction with two distinct enthalpy
  // states. The public temporal block, rather than a density-only pressure
  // update, must preserve rho*h-p while retaining the h0-h offset.
  double coupled_p=100000,coupled_h=175000;
  for(unsigned iteration=0;iteration<8;++iteration) {
    double auxiliary_row[]{-.1,1.1,coupled_h+25000},query_y[2];
    const auto c=esf::auxiliary_eos_coordinates(
        {revision,42,1,2,auxiliary_row},revision,query_y,2);
    const double t=300+(c.query_enthalpy_j_per_kg-hf)/cp;
    const portable::GasSample sample{revision,42,coupled_p,t,coupled_p/(gas_constant*t),
        c.query_enthalpy_j_per_kg,cp,1e-5,.1};
    if(esf::auxiliary_pressure_state(c,sample,revision,42,tangent)!=Status::success)return 36;
    if(!form_pressure_energy_density_jacobian(coupled_p,coupled_h,
        {tangent.density_kg_per_m3,tangent.density_pressure_derivative,tangent.density_enthalpy_derivative},jacobian))return 37;
    PressureEnergyTemporalPoint point;
    point.bdf={1000,-1000,0,1};point.cell_volume=.125;
    point.pressure_absolute=coupled_p;point.enthalpy=coupled_h;point.density=jacobian.density;
    point.accepted_pressure_absolute=100000;point.accepted_enthalpy=170000;point.accepted_density=.7;
    point.target_thermo=jacobian;
    PressureEnergyTemporalLinearization block;
    if(!linearize_pressure_energy_temporal(point,block))return 38;
    const double a=block.continuity_pressure,b=block.continuity_enthalpy;
    const double c_p=block.energy_pressure,d=block.energy_enthalpy;
    const double determinant=a*d-b*c_p;
    if(!std::isfinite(determinant) || determinant==0)return 39;
    coupled_p+=(-block.continuity_residual*d+b*block.energy_residual)/determinant;
    coupled_h+=(c_p*block.continuity_residual-a*block.energy_residual)/determinant;
  }
  if(std::abs(coupled_p-117400)>1e-7 || std::abs(coupled_h-194857.14285714287)>1e-7 ||
     std::abs(raw_density(coupled_p,coupled_h+25000)-.7)>1e-12 ||
     std::abs(.7*coupled_h-coupled_p-19000)>1e-7)return 40;
  std::cout << "ESF dual-state moments PASS: 2/4 fields, GTMC, bounds, identity, lifetime, zero allocation\n";
}
