// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#include "models_esf_detail.hpp"
#include <cmath>
#include <iostream>
#include <limits>
using namespace hundun::v04;
bool near(double actual, double expected, double absolute = 1e-14,
          double relative = 1e-13) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             absolute + relative * std::abs(expected);
}
int main() {
  // Molecular diffusion drives noise with zero SGS viscosity. The limiting
  // species scales the complete Y/Y/h tuple by one common factor.
  std::array<double,3> noise_values{.2,.8,100},noise_lower{0,0,0},noise_upper{1,1,200};
  std::array<double,9> noise_gradient{2,0,0,-2,0,0,10,0,0};
  std::array<double,3> noise_output{99,99,99};
  esf::detail::StochasticSourceRequest noise_request{
      1.,.35,0.,.7,.5,1.,{1,0,0},noise_values.data(),noise_gradient.data(),
      noise_lower.data(),noise_upper.data(),3};
  auto noise_report=esf::detail::stochastic_source(noise_request,noise_output.data());
  if(noise_report.status!=portable::Status::success || !near(noise_report.attenuation,.4) ||
      !near(noise_output[0],.8) || !near(noise_output[1],-.8) ||
      !near(noise_output[2],4) || !near(noise_output[0]+noise_output[1],0)) return 24;
  const auto noise_saved=noise_output;
  noise_gradient[8]=NAN;
  if(esf::detail::stochastic_source(noise_request,noise_output.data()).status!=portable::Status::invalid_input ||
      noise_output!=noise_saved) return 25;
  noise_gradient[8]=0;
  if(esf::detail::stochastic_source(noise_request,noise_values.data()).status!=portable::Status::invalid_input ||
      noise_values[0]!=.2) return 26;
  noise_values[0]=0;noise_request.wiener[0]=-1;
  noise_report=esf::detail::stochastic_source(noise_request,noise_output.data());
  if(noise_report.status!=portable::Status::success || noise_report.attenuation!=0 ||
      noise_output!=std::array<double,3>{0,0,0}) return 27;
  esf::detail::IemSource source;
  if (esf::detail::iem_source(1e-9,1.8e-5,1e-4,2,8,.25,source)!=portable::Status::success ||
      !near(source.implicit_sink_density,236) || !near(source.explicit_source_density,59))
    return 20;
  const auto sentinel=source;
  const double largest=std::numeric_limits<double>::max();
  for (const auto& invalid:std::array<std::array<double,6>,8>{{
      {{0,1.8e-5,1e-4,2,8,.25}}, {{1e-9,-1,1e-4,2,8,.25}},
      {{1e-9,1.8e-5,-1,2,8,.25}}, {{1e-9,1.8e-5,1e-4,-2,8,.25}},
      {{1e-9,1.8e-5,1e-4,2,-8,.25}}, {{1e-9,1.8e-5,1e-4,2,8,NAN}},
      {{1e-9,largest,largest,2,8,.25}}, {{1e-9,1.8e-5,1e-4,2,8,largest}}}}) {
    if (esf::detail::iem_source(invalid[0],invalid[1],invalid[2],invalid[3],
          invalid[4],invalid[5],source)!=portable::Status::invalid_input ||
        source.explicit_source_density!=sentinel.explicit_source_density ||
        source.implicit_sink_density!=sentinel.implicit_sink_density) return 21;
  }
  if (esf::detail::iem_source(1e-9,1.8e-5,1e-4,2,8,-1e6,source)!=portable::Status::success ||
      !near(source.explicit_source_density,-236e6)) return 22;
  if (esf::detail::iem_source(1e-9,1.8e-5,1e-4,0,8,.25,source)!=portable::Status::success ||
      source.explicit_source_density!=0 || source.implicit_sink_density!=0) return 23;
  // Random123 Philox4x32-10 published all-zero known-answer vector.
  auto w = esf::detail::philox_words({0, 0, 0, 0}, {0, 0});
  if (w != std::array<std::uint32_t, 4>{0x6627e8d5, 0xe169c58d, 0xbc57ac4c,
                                        0x9b00dbd8})
    return 1;
  double f{};
  if (esf::detail::iem_factor(2, 1, 1, f) != portable::Status::success ||
      !near(f, 0.36787944117144233, 1e-15, 1e-14))
    return 2;
  combustion::ChemistryIdentity id{1, {{1, {1}, 0}, {1, {1}, 0}}, 42};
  double values[]{0.2, 0.8, 100, 0.8, 0.2, 300};
  esf::detail::Request q;
  q.accepted = {{0, 7, 1}, 42, 2, 2, values};
  q.expected_revision = q.accepted.revision;
  q.identity = &id;
  q.dt_s = 2;
  q.mixing_time_s = 1;
  esf::detail::Workspace work(2);
  double shifted[]{0,1,100,.4,.6,200,0,1,100,.4,.6,200},target[]{.1,.9,180};
  auto centered=work.recenter({{0,7,1},42,2,2,shifted},target);
  if(centered.status!=portable::Status::success ||
      !near(centered.relaxation_factor,.5) ||
      !near(centered.candidate.values[0],0) ||
      !near(centered.candidate.values[3],.2) ||
      !near(centered.means[0],.1) || !near(centered.means[2],180) ||
      !near(centered.candidate.values[2],155) || !near(centered.candidate.values[5],205)) return 28;
  // An interior target preserves the complete fluctuation amplitude.
  target[0]=.3;target[1]=.7;
  centered=work.recenter({{0,7,1},42,2,2,shifted},target);
  if(centered.status!=portable::Status::success || centered.relaxation_factor!=1 ||
      !near(centered.candidate.values[0],.1) || !near(centered.candidate.values[3],.5)) return 29;
  target[0]=0;target[1]=1;
  centered=work.recenter({{0,7,1},42,2,2,shifted},target);
  if(centered.status!=portable::Status::success || centered.relaxation_factor!=0 ||
      centered.candidate.values[0]!=0 || centered.candidate.values[3]!=0 ||
      centered.variances[0]!=0 || centered.means[2]!=180) return 30;
  target[0]=.1;target[1]=.9;
  centered=work.recenter({{0,7,1},42,4,2,shifted},target);
  if(centered.status!=portable::Status::success || !near(centered.relaxation_factor,.5) ||
      !near(centered.variances[0],.01)) return 31;
  // Borrowing both state and mean from the same workspace remains valid.
  centered=work.recenter(centered.candidate,centered.means);
  if(centered.status!=portable::Status::success || centered.relaxation_factor!=1 ||
      !near(centered.variances[0],.01) || !near(centered.means[0],.1)) return 32;
  const auto centered_view=centered.candidate;
  target[0]=-1;
  centered=work.recenter({{0,7,1},42,4,2,shifted},target);
  if(centered.status!=portable::Status::invalid_input || centered.candidate.values ||
      work.valid(centered_view) || shifted[0]!=0) return 33;
  auto r = work.advance(q);
  if (r.status != portable::Status::success ||
      !near(r.candidate.values[0], 0.3896361676485673) ||
      !near(r.variances[0], 0.01218017549129514) || r.means[2] != 200)
    return 3;
  auto old = r.candidate;
  q.accepted.revision.input_revision = 8;
  auto bad = work.advance(q);
  if (bad.status != portable::Status::stale_revision || bad.candidate.values ||
      work.valid(old))
    return 4;
  q.accepted.revision = q.expected_revision;
  for (std::size_t n : {2U, 4U, 6U, 8U, 16U}) {
    std::vector<std::array<double,3>> noise(n), retry(n);
    auto a = esf::detail::balanced_wiener(n, 0.25, {71,0,2,0,0,1},noise.data(),n);
    auto b = esf::detail::balanced_wiener(n, 0.25, {71,0,2,0,0,1},retry.data(),n);
    if (a != portable::Status::success || b != a || noise != retry)
      return 5;
    for (std::size_t p = 0; p < n / 2; ++p)
      for (unsigned d = 0; d < 3; ++d)
        if (noise[2 * p][d] + noise[2 * p + 1][d] != 0 ||
            std::abs(noise[2 * p][d]) != 0.5)
          return 6;
  }
  std::array<std::array<double,3>,16> buffer{};
  for (auto& row : buffer) row.fill(71);
  const auto original=buffer;
  for (std::size_t n : {0U,1U,3U,65U}) {
    if (esf::detail::balanced_wiener(n,.25,{},buffer.data(),buffer.size())!=
          portable::Status::invalid_input || buffer!=original) return 28;
  }
  if (esf::detail::balanced_wiener(16,.25,{},buffer.data(),8)!=
        portable::Status::capacity_exceeded || buffer!=original) return 29;
  if (esf::detail::balanced_wiener(16,NAN,{},buffer.data(),16)!=
        portable::Status::invalid_input || buffer!=original) return 30;
  double raw[]{2, -1}, y[]{0.25, 0.75}, corrected[]{99, 99};
  if (esf::detail::correct_species_flux(y, raw, 2, corrected) !=
          portable::Status::success ||
      corrected[0] != 1.75 || corrected[1] != -1.75)
    return 7;
  raw[1] = NAN;
  if (esf::detail::correct_species_flux(y, raw, 2, corrected) !=
          portable::Status::invalid_input ||
      corrected[0] != 1.75)
    return 8;
  // Zero turbulent diffusivity exactly ignores stochastic gradients;
  // deterministic molecular/source rate remains. No hidden clipping or
  // normalization.
  double rates[]{0.05, -0.05, 10, -0.05, 0.05, -10}, gradients[18];
  for (auto &x : gradients)
    x = NAN;
  q.dt_s = 1;
  q.tcr_control = 0;
  q.deterministic_rates = rates;
  q.gradients = gradients;
  r = work.advance(q);
  if (r.status != portable::Status::success || r.candidate.values[0] != 0.25 ||
      r.candidate.values[2] != 110 || r.candidate.values[5] != 290)
    return 9;
  q.turbulent_diffusivity_m2_s = 0.01;
  if (work.advance(q).candidate.values)
    return 10;
  q.turbulent_diffusivity_m2_s = 0;
  q.dt_s = 100;
  if (work.advance(q).status != portable::Status::conservation_failure ||
      values[0] != 0.2)
    return 11;
  // Four fields, same mean, independently tabulated exp(-2) variance ratio.
  double four[]{0.2, 0.8, 100, 0.8, 0.2, 300, 0.2, 0.8, 100, 0.8, 0.2, 300};
  q.accepted.values = four;
  q.accepted.fields = 4;
  q.dt_s = 2;
  q.tcr_control = 1;
  q.deterministic_rates = nullptr;
  q.gradients = nullptr;
  r = work.advance(q);
  if (r.status != portable::Status::success ||
      !near(r.variances[0], 0.01218017549129514) ||
      r.max_element_residual > 1e-14)
    return 12;
  esf::detail::Workspace small(1);
  if (small.advance(q).status != portable::Status::capacity_exceeded)
    return 13;
  // Uniform initial fields with complementary species gradients. Balanced
  // Wiener gives +/-0.05 exactly in each pair, hence variance 0.0025.
  double uniform[12]{.5, .5, 100, .5, .5, 100, .5, .5, 100, .5, .5, 100};
  double stochastic_gradient[36]{};
  for (std::size_t field = 0; field < 4; ++field) {
    stochastic_gradient[field * 9] = .1;
    stochastic_gradient[field * 9 + 3] = -.1;
  }
  q.accepted.values = uniform;
  q.gradients = stochastic_gradient;
  q.dt_s = 1;
  q.tcr_control = 0;
  q.turbulent_diffusivity_m2_s = .125;
  r = work.advance(q);
  if (r.status != portable::Status::success || !near(r.means[0], .5) ||
      !near(r.variances[0], .0025))
    return 14;
  std::cout << "ESF focused PASS\n";
}
