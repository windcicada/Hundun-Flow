// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dyn711_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::tcr::detail {
Dyn711Tick dyn711_tick(Dyn711Clock c) noexcept {
  if (c.rate_intervals > 8 || c.cphi_count > 12 ||
      (c.cphi_count > 0 && c.cphi_count < 6)) return {};
  Dyn711Tick t;
  t.available = true;
  t.evaluate_rates = c.rate_intervals == 8;
  t.update_cphi = c.cphi_count <= 6;
  t.next = {t.evaluate_rates ? 0u : c.rate_intervals + 1,
            t.update_cphi ? c.cphi_count + 6 : c.cphi_count - 1,
            c.rates_initialized || t.evaluate_rates};
  return t;
}

Dyn711Control dyn711_species_control(double eta, double pdf, double psr,
    double chemical_time, double flow_time, double weak) noexcept {
  Dyn711Control out;
  if (!std::isfinite(eta) || eta < 0 || eta > 1 ||
      !std::isfinite(pdf) || !std::isfinite(psr) ||
      !std::isfinite(chemical_time) || chemical_time < 0 ||
      !std::isfinite(flow_time) || flow_time <= 0 ||
      !std::isfinite(weak) || weak <= 0) return out;
  const long double e = eta, one_minus_e = 1 - e;
  long double ratio = std::abs(psr) <= weak ? 1 :
      static_cast<long double>(pdf) / psr;
  if (ratio <= 0) ratio = 1;
  if (eta > 0 && eta < 1)
    ratio = std::min(ratio, 1 / (4 * e * one_minus_e));
  const bool upper = ratio >= 1 &&
      static_cast<long double>(chemical_time) > 100.L * flow_time;
  long double selected;
  if (eta == 1) {
    // The quadratic becomes linear: kappa = ratio. The diverging upper root
    // is outside that equation's finite solution set.
    selected = ratio;
    out.linear_endpoint = true;
  } else {
    const auto root = std::sqrt(std::max(0.L, 1 - 4 * e * one_minus_e * ratio));
    selected = upper ? (1 + root) / (2 * one_minus_e) :
        (2 * e * ratio) / (1 + root);
    out.upper_branch = upper;
  }
  out.ratio = static_cast<double>(ratio);
  out.selected = static_cast<double>(selected);
  out.effective = std::min(1., out.selected);
  out.available = std::isfinite(out.ratio) && std::isfinite(out.selected);
  return out;
}

bool dyn711_advance_rate(const Dyn711RateState &accepted, Dyn711Clock clock,
    double dt, double pdf, double psr, double eta, double chemical_time,
    double flow_time, double weak, Dyn711RateState &candidate) noexcept {
  const auto tick = dyn711_tick(clock);
  if (!tick.available || !std::isfinite(dt) || dt <= 0 ||
      !std::isfinite(pdf) || !std::isfinite(psr) ||
      !std::isfinite(accepted.pdf_sum) || !std::isfinite(accepted.psr_sum) ||
      !std::isfinite(accepted.selected) || accepted.selected < 0 ||
      !std::isfinite(accepted.effective) ||
      accepted.effective != std::min(1., accepted.selected)) return false;
  auto next = accepted;
  if (tick.evaluate_rates) {
    const auto c = dyn711_species_control(eta, accepted.pdf_sum, accepted.psr_sum,
        chemical_time, flow_time, weak);
    if (!c.available) return false;
    next = {0, 0, c.selected, c.effective, c.upper_branch, c.linear_endpoint};
  } else {
    next.pdf_sum += dt * 1e5 * pdf;
    next.psr_sum += dt * 1e5 * psr;
    if (!clock.rates_initialized) {
      next.selected = next.effective = .2;
      next.upper_branch = next.linear_endpoint = false;
    }
    if (!std::isfinite(next.pdf_sum) || !std::isfinite(next.psr_sum)) return false;
  }
  candidate = next;
  return true;
}

bool dyn711_filter_moments(const Dyn711FilterDonor *donors, unsigned count,
    DynamicFilterMoments &moments) noexcept {
  if (!donors || count == 0 || count > 8) return false;
  DynamicFilterMoments m;
  double mass{};
  std::array<double,3> gradient{};
  for (unsigned i=0;i<count;++i) {
    const auto &d=donors[i];
    if (!std::isfinite(d.density) || d.density<=0 ||
        !std::isfinite(d.volume) || d.volume<=0 || !std::isfinite(d.scalar)) return false;
    const double w=d.density*d.volume;
    mass+=w;
    m.density+=d.density*w;
    m.delta_squared+=d.density*std::pow(d.volume,5./3.);
    m.scalar+=w*d.scalar;
    m.density_scalar_squared+=d.density*w*d.scalar*d.scalar;
    double g2{};
    for (unsigned axis=0;axis<3;++axis) {
      if (!std::isfinite(d.gradient[axis])) return false;
      gradient[axis]+=w*d.gradient[axis];
      g2+=d.gradient[axis]*d.gradient[axis];
    }
    m.density_delta_squared_gradient_squared+=
        d.density*d.density*std::pow(d.volume,5./3.)*g2;
  }
  if (!std::isfinite(mass) || mass<=0) return false;
  m.density=1e6*m.density/mass;
  m.delta_squared/=mass;
  m.scalar/=mass;
  m.density_scalar_squared=1e6*m.density_scalar_squared/mass;
  m.density_delta_squared_gradient_squared=1e6*m.density_delta_squared_gradient_squared/mass;
  for (double g:gradient) m.gradient_squared+=(g/mass)*(g/mass);
  for (double v:{m.density,m.delta_squared,m.scalar,m.density_scalar_squared,
      m.density_delta_squared_gradient_squared,m.gradient_squared})
    if (!std::isfinite(v)) return false;
  moments=m;
  return true;
}
DynamicFilterProducts dyn711_filter_products(const DynamicFilterMoments &v) noexcept {
  for (double x : {v.density, v.delta_squared, v.gradient_squared,
      v.density_delta_squared_gradient_squared, v.density_scalar_squared, v.scalar})
    if (!std::isfinite(x)) return {};
  if (v.density <= 0 || v.delta_squared <= 0 || v.gradient_squared < 0 ||
      v.density_delta_squared_gradient_squared < 0 || v.density_scalar_squared < 0)
    return {};
  const double l = std::abs(v.density_scalar_squared - v.density * v.scalar * v.scalar);
  const double m = std::abs(v.density * v.delta_squared * v.gradient_squared -
      v.density_delta_squared_gradient_squared);
  const double m2 = m * m, lm = l * m;
  return {std::isfinite(m2) && std::isfinite(lm), m2, lm};
}
double dyn711_filter_ratio(double m2, double lm) noexcept {
  if (!std::isfinite(m2) || !std::isfinite(lm) || m2 < 0 || lm < 0)
    return std::numeric_limits<double>::quiet_NaN();
  return m2 == 0 || lm == 0 ? -1 : m2 / lm;
}
double dyn711_select_cphi(const std::array<double, 3> &ratios) noexcept {
  for (double c : ratios)
    if (!std::isfinite(c)) return std::numeric_limits<double>::quiet_NaN();
  double c;
  if (ratios[2] > 1.5 && ratios[2] < 12) c = ratios[2];
  else if (ratios[1] > 1.5 && ratios[1] < 12) c = ratios[1];
  else if (ratios[0] > 1.5 && ratios[0] < 12) c = ratios[0];
  else {
    c = std::abs(ratios[2]-4) <= std::abs(ratios[1]-4) ? ratios[2] : ratios[1];
    if (std::abs(ratios[0]-4) < std::abs(c-4)) c = ratios[0];
  }
  return std::clamp(c <= 0 ? 2. : c, 1., 16.);
}
double dyn711_smooth_cphi(double centre, const std::array<double, 7> &neighbors) noexcept {
  if (!std::isfinite(centre) || centre < 1 || centre > 16)
    return std::numeric_limits<double>::quiet_NaN();
  double sum = 2 * centre;
  unsigned count = 2;
  for (double c : neighbors) {
    if (!std::isfinite(c) || c < 1 || c > 16)
      return std::numeric_limits<double>::quiet_NaN();
    if (c > 1.5 && c < 12) { sum += c; ++count; }
  }
  return sum / count;
}
} // namespace hundun::v04::tcr::detail
