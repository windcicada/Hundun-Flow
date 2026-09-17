// SPDX-License-Identifier: Apache-2.0
#include "models_tcr_dynamic_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::tcr::detail {
SpeciesControl cdphyso_species_control(double eta, double pdf, double psr,
                                       double weak) noexcept {
  SpeciesControl out;
  if (!std::isfinite(eta) || !std::isfinite(pdf) || !std::isfinite(psr) ||
      !std::isfinite(weak) || weak <= 0) return out;
  pdf = std::abs(pdf);
  psr = std::abs(psr);
  if (pdf <= weak && psr <= weak) {
    out.available = true;
    return out;
  }
  out.rate_ratio = pdf / std::max(psr, weak);
  if (!std::isfinite(out.rate_ratio)) return out;
  if (psr <= weak) {
    out.state = SpeciesControlState::similar;
    out.effective = out.selected = 0;
  } else {
    eta = std::clamp(eta, 1e-6, 1 - 1e-6);
    const double discriminant = 1 - (4 * eta * (1 - eta)) * out.rate_ratio;
    if (discriminant >= 0) {
      // Rationalized small root retains trace-reaction precision.
      out.selected = (2 * eta * out.rate_ratio) / (1 + std::sqrt(discriminant));
      out.state = SpeciesControlState::direct;
    } else {
      const double inverse = 1 / out.rate_ratio;
      const double projected = inverse / (2 * (1 + std::sqrt(std::max(0., 1 - inverse))));
      eta = std::clamp(projected, 1e-6, 1 - 1e-6);
      out.selected = 1 / (2 * (1 - eta));
      out.state = SpeciesControlState::projected;
    }
    out.effective = std::clamp(out.selected, 1e-4, 1.);
  }
  out.available = std::isfinite(out.selected);
  return out;
}
bool dynamic_filter_moments(const DynamicFilterDonor *donors, unsigned count,
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
DynamicFilterProducts dynamic_filter_products(const DynamicFilterMoments &v) noexcept {
  for (double x : {v.density, v.delta_squared, v.gradient_squared,
                   v.density_delta_squared_gradient_squared,
                   v.density_scalar_squared, v.scalar})
    if (!std::isfinite(x)) return {};
  if (v.density <= 0 || v.delta_squared <= 0 || v.gradient_squared < 0 ||
      v.density_delta_squared_gradient_squared < 0 || v.density_scalar_squared < 0)
    return {};
  const double l = std::abs(v.density_scalar_squared - v.density * v.scalar * v.scalar);
  const double m = std::abs(v.density * v.delta_squared * v.gradient_squared -
                             v.density_delta_squared_gradient_squared);
  if (!std::isfinite(l) || !std::isfinite(m)) return {};
  // The frozen reference bounds each donor's products before its second
  // spatial filter; retain that weighting as part of this model identity.
  const double scale = std::max({1., l, m});
  const double ls = l / scale, ms = m / scale;
  return {true, ms * ms, ls * ms};
}
double dynamic_cd_from_products(double m2, double lm) noexcept {
  if (!std::isfinite(m2) || !std::isfinite(lm) || m2 < 0 || lm < 0)
    return std::numeric_limits<double>::quiet_NaN();
  if (lm <= 0 || m2 <= 0) return 2.;
  const double denominator = std::max(lm, 1e-30);
  const double ratio = m2 >= 64 * denominator ? 64 : m2 / denominator;
  if (ratio < .1) return 2.;
  return std::clamp(ratio, 1., 16.);
}
double dynamic_group_cd(const std::array<double, 3> &groups,
                         MixingGroup group) noexcept {
  for (double c : groups)
    if (!std::isfinite(c) || c < 1 || c > 16)
      return std::numeric_limits<double>::quiet_NaN();
  const double c = group == MixingGroup::radical ? groups[2] :
      std::sqrt(groups[group == MixingGroup::product ? 0 : 1] * groups[2]);
  return std::clamp(c, 1., 10.);
}
} // namespace hundun::v04::tcr::detail
