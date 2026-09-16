// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_spray.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace hundun::v04::spray {
namespace {

bool nonnegative(double value) noexcept {
  return std::isfinite(value) && value >= 0.0;
}

bool quantile(double lower, double upper, double lambda, double uniform,
              double& result) noexcept {
  constexpr std::size_t count = 20U;
  std::array<double, count> diameter{}, density{}, cumulative{};
  const double spacing = (upper-lower)/static_cast<double>(count-1U);
  if (!(spacing > 0.0)) return false;
  for (std::size_t i=0U;i<count;++i)
    diameter[i] = lower+static_cast<double>(i)*spacing;
  for (std::size_t i=1U;i+1U<count;++i) {
    const double d = diameter[i];
    density[i] = d*d*(std::pow(d,2.0/3.0)-lambda)*
        (std::pow(1.0-d*d*d,2.0/9.0)-lambda);
    if (!std::isfinite(density[i])) return false;
    density[i] = std::max(0.0,density[i]);
  }
  double integral{};
  for (std::size_t i=1U;i<count;++i)
    integral += 0.5*(density[i-1U]+density[i])*spacing;
  if (!(integral > 0.0) || !std::isfinite(integral)) return false;
  for (auto& value:density) value /= integral;
  for (std::size_t i=1U;i+1U<count;++i)
    cumulative[i] = cumulative[i-1U]+0.5*(density[i-1U]+density[i])*spacing;
  cumulative.back() = 1.0;
  if (uniform == 0.0) {result=lower;return true;}
  if (uniform == 1.0) {result=upper;return true;}
  for (std::size_t i=0U;i+1U<count;++i) {
    if (uniform > cumulative[i] && uniform <= cumulative[i+1U]) {
      result = diameter[i]+(uniform-cumulative[i])*(diameter[i+1U]-diameter[i])/
          (cumulative[i+1U]-cumulative[i]);
      return std::isfinite(result) && result > 0.0 && result < 1.0;
    }
  }
  return false;
}

}  // namespace

SgsBreakupReport evaluate_sgs_breakup(const SgsBreakupInput& input) noexcept {
  SgsBreakupReport report;
  const auto& old = input.history;
  const double positive[]{input.droplet_diameter_m,input.gas_density_kg_per_m3,
      input.liquid_density_kg_per_m3,input.surface_tension_n_per_m,
      input.duration_s,input.deterministic_time_coefficient};
  const double nonnegative_values[]{input.relative_speed_m_per_s,input.gas_dynamic_viscosity_pa_s,
      input.dissipation_m2_per_s3,input.stochastic_coefficient,old.mean_dissipation_m2_per_s3,
      old.dissipation_age_s,old.mean_rate_per_s,old.rate_age_s,input.daughter_uniform_01};
  for (double value:positive) if (!std::isfinite(value) || !(value > 0.0)) return report;
  for (double value:nonnegative_values) if (!nonnegative(value)) return report;
  if (old.poisson_multiplier > 7U || input.daughter_uniform_01 > 1.0) return report;
  report.status = SgsBreakupStatus::numerical_failure;
  auto& history = report.candidate;
  history = old;
  const double d = input.droplet_diameter_m,dt = input.duration_s;
  history.dissipation_age_s += dt;
  history.mean_dissipation_m2_per_s3 =
      (old.mean_dissipation_m2_per_s3*old.dissipation_age_s+input.dissipation_m2_per_s3*dt)/
      history.dissipation_age_s;
  if (history.mean_dissipation_m2_per_s3 == 0.0)
    history.mean_dissipation_m2_per_s3 = 1.0e-8;
  const double epsilon = history.mean_dissipation_m2_per_s3;
  constexpr double beta = 8.2;
  const double capillary = 12.0*input.surface_tension_n_per_m/(input.liquid_density_kg_per_m3*beta);
  report.weber_number = input.gas_density_kg_per_m3*input.relative_speed_m_per_s*
      input.relative_speed_m_per_s*d/input.surface_tension_n_per_m;
  report.critical_diameter_m = std::pow(capillary,3.0/5.0)*std::pow(epsilon,-2.0/5.0);
  report.distribution_lambda = std::pow(report.critical_diameter_m/d,5.0/3.0);
  report.minimum_diameter_ratio = std::pow(capillary,1.5)*std::pow(d,-2.5)/epsilon;
  if (report.minimum_diameter_ratio < 1.0)
    report.maximum_diameter_ratio = std::cbrt(1.0-std::pow(report.minimum_diameter_ratio,3.0));
  if (input.dissipation_m2_per_s3 > 0.0) {
    report.kolmogorov_length_m = std::pow(std::pow(input.gas_dynamic_viscosity_pa_s/
        input.gas_density_kg_per_m3,3.0)/input.dissipation_m2_per_s3,0.25);
    report.kolmogorov_length_available = true;
  }
  report.rate_active = d > report.critical_diameter_m;
  if (report.rate_active) {
    if (input.relative_speed_m_per_s > 0.0) {
      const double breakup_time = input.deterministic_time_coefficient*
          std::sqrt(input.liquid_density_kg_per_m3/input.gas_density_kg_per_m3)*
          ((d/2.0)/input.relative_speed_m_per_s);
      report.deterministic_rate_per_s = 1.0/breakup_time;
    }
    if (input.stochastic_enabled) {
      const double driving = beta*std::pow(epsilon*d,2.0/3.0)-
          12.0*input.surface_tension_n_per_m/(input.liquid_density_kg_per_m3*d);
      if (driving > 0.0)
        report.stochastic_rate_per_s = input.stochastic_coefficient*std::sqrt(driving)/d;
    }
    history.rate_age_s += dt;
    history.mean_rate_per_s = (old.mean_rate_per_s*old.rate_age_s+
        (report.deterministic_rate_per_s+report.stochastic_rate_per_s)*dt)/history.rate_age_s;
    report.breakup_requested = old.poisson_multiplier > 0U &&
        report.maximum_diameter_ratio > report.minimum_diameter_ratio &&
        history.mean_rate_per_s > 0.0 && history.rate_age_s >=
          1.0/(static_cast<double>(old.poisson_multiplier)*history.mean_rate_per_s);
    if (report.breakup_requested) {
      double ratio{};
      if (!quantile(report.minimum_diameter_ratio,report.maximum_diameter_ratio,
          report.distribution_lambda,input.daughter_uniform_01,ratio)) return report;
      report.daughter_diameter_ratios = {ratio,std::cbrt(1.0-ratio*ratio*ratio)};
      // The upper endpoint can round to one while its complementary daughter
      // remains representable. Preserve the known lower endpoint explicitly.
      if (input.daughter_uniform_01 == 1.0)
        report.daughter_diameter_ratios[1] = report.minimum_diameter_ratio;
    }
  }
  const double values[]{history.dissipation_age_s,history.mean_dissipation_m2_per_s3,
      history.rate_age_s,history.mean_rate_per_s,report.weber_number,report.critical_diameter_m,
      report.distribution_lambda,report.minimum_diameter_ratio,report.maximum_diameter_ratio,
      report.deterministic_rate_per_s,report.stochastic_rate_per_s,report.kolmogorov_length_m};
  for (double value:values) if (!nonnegative(value)) return report;
  if (report.breakup_requested && !(report.daughter_diameter_ratios[0] > 0.0 &&
      report.daughter_diameter_ratios[1] > 0.0)) return report;
  report.status = SgsBreakupStatus::success;
  return report;
}

}  // namespace hundun::v04::spray
