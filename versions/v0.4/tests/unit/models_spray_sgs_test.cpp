// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_spray.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace hundun::v04::spray;

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
        !(result.daughter_diameter_ratios[0] > 0. && result.daughter_diameter_ratios[1] > 0.)) return 7;
  }
  input.daughter_uniform_01=.5;
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
            << " limit=1e-11 endpoints=pass clocks=pass zero_dissipation=pass invalid=7 passed=1\n";
}
