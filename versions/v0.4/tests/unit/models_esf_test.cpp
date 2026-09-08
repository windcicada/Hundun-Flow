// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#include "models_esf_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04;
bool near(double actual, double expected, double absolute = 1e-14,
          double relative = 1e-13) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             absolute + relative * std::abs(expected);
}
int main() {
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
  for (std::size_t n : {2U, 4U}) {
    auto noise = esf::detail::balanced_wiener(n, 0.25, {71, 0, 2, 0, 0, 1});
    auto retry = esf::detail::balanced_wiener(n, 0.25, {71, 0, 2, 0, 0, 1});
    if (noise.status != portable::Status::success ||
        noise.increments != retry.increments)
      return 5;
    for (std::size_t p = 0; p < n / 2; ++p)
      for (unsigned d = 0; d < 3; ++d)
        if (noise.increments[2 * p][d] + noise.increments[2 * p + 1][d] != 0 ||
            std::abs(noise.increments[2 * p][d]) != 0.5)
          return 6;
  }
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
