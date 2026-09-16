// SPDX-License-Identifier: Apache-2.0
#include "models_spray_breakup_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace hundun::v04::spray;

namespace {
using namespace hundun::v04::spray::detail;
constexpr double kPi = 3.141592653589793238462643383279502884;

bool close(long double a, long double b, long double tolerance = 2.e-12L) {
  return std::isfinite(a) && std::isfinite(b) &&
      std::abs(a-b) <= tolerance*std::max({1.e-300L,std::abs(a),std::abs(b)});
}

bool same(const SprayParcelState& a, const SprayParcelState& b) {
  return a.id == b.id && a.position_m == b.position_m &&
      a.velocity_m_per_s == b.velocity_m_per_s && a.droplet_mass_kg == b.droplet_mass_kg &&
      a.droplet_diameter_m == b.droplet_diameter_m && a.multiplicity == b.multiplicity &&
      a.temperature_k == b.temperature_k && a.age_s == b.age_s &&
      a.owner_global_cell == b.owner_global_cell &&
      a.liquid_material_fingerprint == b.liquid_material_fingerprint;
}

SgsChildInput child_input(const SgsBreakupInput& kernel, const SgsBreakupReport& trigger) {
  SgsChildInput input;
  input.parent.id = {1234U,5678U};
  input.parent.position_m = {.1,.2,.3};
  input.parent.velocity_m_per_s = {20.,-3.,4.};
  input.parent.droplet_diameter_m = kernel.droplet_diameter_m;
  input.parent.droplet_mass_kg = kernel.liquid_density_kg_per_m3*kPi/6.*
      std::pow(kernel.droplet_diameter_m,3.);
  input.parent.multiplicity = 37.25;
  input.parent.temperature_k = 420.;
  input.parent.liquid_material_fingerprint = 91U;
  input.parent.owner_global_cell = 17U;
  input.parent.age_s = .025;
  input.trigger = trigger;
  input.accepted_step = 12U;
  input.breakup_ordinal = 9U;
  input.liquid_density_kg_per_m3 = kernel.liquid_density_kg_per_m3;
  input.surface_tension_n_per_m = kernel.surface_tension_n_per_m;
  input.liquid_absolute_thermochemical_enthalpy_j_per_kg = -2.e5;
  return input;
}

bool check_children(const SgsChildInput& input) {
  const auto before = input.parent;
  const auto result = generate_sgs_children(input);
  const auto retry = generate_sgs_children(input);
  if (!result.succeeded() || !retry.succeeded() || result.candidate.child_count != 2U ||
      result.model_id != "conservative_sgs_binary_split_v1" ||
      result.trigger_model_id != input.trigger.model_id ||
      result.parent_id != input.parent.id || result.accepted_step != input.accepted_step ||
      result.breakup_ordinal != input.breakup_ordinal || !same(before,input.parent)) return false;
  long double mass{},surface{},kinetic{},enthalpy{};
  std::array<long double,3> momentum{};
  auto alternate = input;
  ++alternate.breakup_ordinal;
  const auto next = generate_sgs_children(alternate);
  alternate = input;
  ++alternate.accepted_step;
  const auto next_step = generate_sgs_children(alternate);
  if (!next.succeeded() || !next_step.succeeded()) return false;
  for (unsigned i=0;i<2U;++i) {
    const auto& child = result.candidate.children[i];
    if (!same(child,retry.candidate.children[i]) || child.id == input.parent.id ||
        child.id == next.candidate.children[i].id || child.id == next_step.candidate.children[i].id ||
        child.id == result.candidate.children[1U-i].id ||
        child.multiplicity != input.parent.multiplicity || child.age_s != input.parent.age_s ||
        child.position_m != input.parent.position_m || child.velocity_m_per_s != input.parent.velocity_m_per_s ||
        child.temperature_k != input.parent.temperature_k ||
        child.owner_global_cell != input.parent.owner_global_cell ||
        child.liquid_material_fingerprint != input.parent.liquid_material_fingerprint ||
        !(child.droplet_mass_kg > 0.) || !(child.droplet_diameter_m > 0.) ||
        !close(child.droplet_diameter_m/input.parent.droplet_diameter_m,
               input.trigger.daughter_diameter_ratios[i])) return false;
    const long double d=child.droplet_diameter_m;
    if (!close(child.droplet_mass_kg,input.liquid_density_kg_per_m3*kPi/6.L*d*d*d)) return false;
    const long double represented=static_cast<long double>(child.droplet_mass_kg)*child.multiplicity;
    mass += represented;
    surface += child.multiplicity*input.surface_tension_n_per_m*kPi*d*d;
    enthalpy += represented*input.liquid_absolute_thermochemical_enthalpy_j_per_kg;
    for (unsigned c=0;c<3U;++c) {
      momentum[c] += represented*child.velocity_m_per_s[c];
      kinetic += .5L*represented*child.velocity_m_per_s[c]*child.velocity_m_per_s[c];
    }
  }
  const long double old_mass=static_cast<long double>(input.parent.droplet_mass_kg)*input.parent.multiplicity;
  if (!close(mass,old_mass) || !close(enthalpy,old_mass*input.liquid_absolute_thermochemical_enthalpy_j_per_kg) ||
      !close(kinetic,old_mass*.5L*(400.L+9.L+16.L))) return false;
  for (unsigned c=0;c<3U;++c)
    if (!close(momentum[c],old_mass*input.parent.velocity_m_per_s[c])) return false;
  const auto& budget=result.conservation;
  if (!close(budget.children_spherical_surface_energy_j,surface) ||
      !close(budget.children_total_mass_kg,mass) ||
      budget.children_total_multiplicity != 2.*input.parent.multiplicity ||
      budget.supplied_deformation_energy_j != 0. ||
      budget.unassigned_deformation_energy_j != -budget.surface_energy_increase_j) return false;
  alternate=input;
  alternate.parent.owner_global_cell=999U;
  const auto moved=generate_sgs_children(alternate);
  if (!moved.succeeded()) return false;
  for (unsigned i=0;i<2U;++i)
    if (moved.candidate.children[i].id != result.candidate.children[i].id) return false;
  return true;
}

bool invalid_children(const SgsChildInput& valid) {
  for (unsigned i=0;i<13U;++i) {
    auto input=valid;
    auto expected=BreakupChildStatus::invalid_input;
    if (i==0U) input.trigger.daughter_diameter_ratios={.5,.5};
    if (i==1U) input.trigger.daughter_diameter_ratios[0]=0.;
    if (i==2U) input.trigger.daughter_diameter_ratios[0]=std::numeric_limits<double>::quiet_NaN();
    if (i==3U) input.trigger.model_id="unknown";
    if (i==4U) input.trigger.candidate.poisson_multiplier=0U;
    if (i==5U) input.trigger.candidate.rate_age_s=0.;
    if (i==6U) input.trigger.candidate.mean_dissipation_m2_per_s3=-1.;
    if (i==7U) input.liquid_density_kg_per_m3=0.;
    if (i==8U) input.liquid_absolute_thermochemical_enthalpy_j_per_kg=std::numeric_limits<double>::infinity();
    if (i==9U) {input.parent.droplet_mass_kg*=1.1;expected=BreakupChildStatus::inconsistent_parent_geometry;}
    if (i==10U) {input.trigger.breakup_requested=false;expected=BreakupChildStatus::breakup_not_requested;}
    if (i==11U) {input.parent.multiplicity=std::numeric_limits<double>::max();expected=BreakupChildStatus::non_finite_output;}
    if (i==12U) input.trigger.status=SgsBreakupStatus::numerical_failure;
    const auto result=generate_sgs_children(input);
    if (result.status != expected || result.candidate.available || result.candidate.child_count != 0U ||
        result.model_id != "conservative_sgs_binary_split_v1" || result.conservation.parent_total_mass_kg != 0.) {
      std::cerr << "invalid child input " << i << " status " << unsigned(result.status) << '\n';
      return false;
    }
    for (const auto& child:result.candidate.children) if (!same(child,{})) return false;
  }
  return true;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  std::ifstream data(argv[1]);
  unsigned count{};
  data >> count;
  if (!data || count != 324U) return 2;
  double worst{};
  unsigned events{};
  for (unsigned i=0;i<count;++i) {
    SgsBreakupInput input;
    unsigned poisson{},stochastic{};
    auto& h=input.history;
    data >> h.mean_dissipation_m2_per_s3 >> h.dissipation_age_s >> h.mean_rate_per_s
         >> h.rate_age_s >> poisson >> input.droplet_diameter_m >> input.relative_speed_m_per_s
         >> input.gas_density_kg_per_m3 >> input.liquid_density_kg_per_m3
         >> input.surface_tension_n_per_m >> input.gas_dynamic_viscosity_pa_s
         >> input.dissipation_m2_per_s3 >> input.duration_s >> input.deterministic_time_coefficient
         >> input.stochastic_coefficient >> stochastic >> input.daughter_uniform_01;
    h.poisson_multiplier=static_cast<std::uint8_t>(poisson);
    input.stochastic_enabled=stochastic != 0U;
    std::array<double,13> expected{};
    for (auto& value:expected) data >> value;
    if (!data) return 2;
    const auto result=evaluate_sgs_breakup(input);
    if (!result.succeeded()) {std::cerr << "request " << i << " status " << unsigned(result.status) << '\n';return 3;}
    const auto& next=result.candidate;
    const std::array<double,13> actual{next.mean_dissipation_m2_per_s3,next.dissipation_age_s,
        next.mean_rate_per_s,next.rate_age_s,result.weber_number,result.deterministic_rate_per_s,
        result.stochastic_rate_per_s,result.kolmogorov_length_m,result.breakup_requested ? 2.0 : 0.0,
        result.daughter_diameter_ratios[0],result.daughter_diameter_ratios[1],
        result.breakup_requested ? 1.0 : 0.0,result.breakup_requested ? 1.0 : 0.0};
    for (std::size_t c=0;c<actual.size();++c) {
      const double scale=std::max({1.e-300,std::abs(actual[c]),std::abs(expected[c])});
      const double error=std::abs(actual[c]-expected[c])/scale;
      worst=std::max(worst,error);
      if (!std::isfinite(actual[c]) || error > 1.e-11) {
        std::cerr << "request " << i << " component " << c << " relative " << error
                  << " actual " << actual[c] << " expected " << expected[c] << '\n';return 4;
      }
    }
    if (result.breakup_requested) {
      ++events;
      if (!check_children(child_input(input,result))) {std::cerr << "children request " << i << '\n';return 11;}
      const auto d=result.daughter_diameter_ratios;
      if (std::abs(d[0]*d[0]*d[0]+d[1]*d[1]*d[1]-1.) > 2.e-15) return 5;
    }
  }
  if (events == 0U || events == count) return 6;
  data >> std::ws;
  if (!data.eof()) return 2;
  SgsBreakupInput input;
  input.droplet_diameter_m=1.e-3;input.relative_speed_m_per_s=20.;
  input.gas_density_kg_per_m3=1.;input.liquid_density_kg_per_m3=750.;
  input.surface_tension_n_per_m=.025;input.gas_dynamic_viscosity_pa_s=2.e-5;
  input.dissipation_m2_per_s3=1.e8;input.duration_s=.01;input.history.poisson_multiplier=7U;
  for (double u:{0.,1.}) {
    input.daughter_uniform_01=u;
    const auto result=evaluate_sgs_breakup(input);
    if (!result.succeeded() || !result.breakup_requested ||
        !(result.daughter_diameter_ratios[0] > 0. && result.daughter_diameter_ratios[1] > 0.) ||
        !check_children(child_input(input,result))) return 7;
  }
  input.daughter_uniform_01=.5;
  if (!invalid_children(child_input(input,evaluate_sgs_breakup(input)))) return 12;
  auto boundary=input;
  boundary.droplet_diameter_m=evaluate_sgs_breakup(input).critical_diameter_m;
  boundary.history.mean_rate_per_s=17.;boundary.history.rate_age_s=.02;
  auto result=evaluate_sgs_breakup(boundary);
  if (!result.succeeded() || result.rate_active || result.breakup_requested ||
      result.candidate.mean_rate_per_s != 17. || result.candidate.rate_age_s != .02) return 8;
  auto zero=input;
  zero.dissipation_m2_per_s3=0.;zero.relative_speed_m_per_s=0.;
  result=evaluate_sgs_breakup(zero);
  if (!result.succeeded() || result.kolmogorov_length_available || result.breakup_requested ||
      result.candidate.mean_dissipation_m2_per_s3 != 1.e-8) return 9;
  for (unsigned i=0;i<7U;++i) {
    auto bad=input;
    if (i==0U) bad.gas_density_kg_per_m3=0.;
    if (i==1U) bad.dissipation_m2_per_s3=-1.;
    if (i==2U) bad.daughter_uniform_01=1.1;
    if (i==3U) bad.deterministic_time_coefficient=0.;
    if (i==4U) bad.history.poisson_multiplier=8U;
    if (i==5U) bad.duration_s=0.;
    if (i==6U) bad.history.mean_rate_per_s=std::numeric_limits<double>::quiet_NaN();
    if (evaluate_sgs_breakup(bad).status != SgsBreakupStatus::invalid_input) return 10;
  }
  std::cout << "sgs_breakup requests=" << count << " events=" << events << " max_relative=" << worst
            << " limit=1e-11 endpoints=pass clocks=pass zero_dissipation=pass invalid=7 binary_children=90 child_invalid=13 passed=1\n";
}
