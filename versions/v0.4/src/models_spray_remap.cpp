// SPDX-License-Identifier: Apache-2.0
#include "models_spray_remap_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::spray::detail {
namespace {
// Positive cp over the whole bracket makes the inversion single-valued.
// Cubic extrema are checked analytically, including an interior minimum.
bool monotone(const LiquidAsset& a) noexcept {
  const double lo = a.pack.minimum_temperature_k;
  const double hi = a.pack.maximum_temperature_k;
  if (!std::isfinite(lo) || !std::isfinite(hi) || !(lo > 0) || !(hi > lo))
    return false;
  const auto valid = [&](double t) {
    return evaluate_liquid_enthalpy(a, t).available;
  };
  if (!valid(lo) || !valid(hi)) return false;
  const auto& c = a.pack.cp_j_per_kg_k;
  if (c.kind == TemperatureCorrelationKind::constant) return true;
  // This fixed correlation has positive derivative throughout its admitted
  // subcritical domain: 7.4680836 - 0.0084247974*T + 14948.931/(684.26-T)^2.
  if (c.kind == TemperatureCorrelationKind::kerosene_cp_v1) return true;
  if (c.kind != TemperatureCorrelationKind::polynomial_cubic ||
      !std::isfinite(c.reference_temperature_k)) return false;
  for (double x : c.c) if (!std::isfinite(x)) return false;
  const long double aa = 3.L*c.c[3], bb = 2.L*c.c[2], cc = c.c[1];
  const auto root = [&](long double x) {
    const long double t = x+c.reference_temperature_k;
    return !(t > lo && t < hi) || valid(static_cast<double>(t));
  };
  if (aa == 0) return bb == 0 || root(-cc/bb);
  const long double disc = bb*bb-4*aa*cc;
  if (disc < 0) return true;
  const long double q = -.5L*(bb+std::copysign(std::sqrt(disc), bb));
  if (q == 0) return root(-bb/(2*aa));
  return root(q/aa) && root(cc/q);
}
} // namespace
LiquidRemapReport remap_liquid_mass_enthalpy(
    const LiquidAsset& a, double mass, double h) noexcept {
  LiquidRemapReport out;
  if (!std::isfinite(mass) || !(mass > 0) || !std::isfinite(h) || !monotone(a))
    return out;
  double lo = a.pack.minimum_temperature_k, hi = a.pack.maximum_temperature_k;
  auto lower = evaluate_liquid_enthalpy(a, lo);
  auto upper = evaluate_liquid_enthalpy(a, hi);
  if (h < lower.liquid_enthalpy_j_per_kg || h > upper.liquid_enthalpy_j_per_kg)
    return out;
  for (unsigned i = 0; i < 80; ++i) {
    const double mid = lo+.5*(hi-lo);
    if (mid == lo || mid == hi) break;
    const auto sample = evaluate_liquid_enthalpy(a, mid);
    if (!sample.available) return out;
    if (sample.liquid_enthalpy_j_per_kg < h) { lo = mid; lower = sample; }
    else { hi = mid; upper = sample; }
  }
  const bool use_lower = std::abs(lower.liquid_enthalpy_j_per_kg-h) <=
                         std::abs(upper.liquid_enthalpy_j_per_kg-h);
  const double t = use_lower ? lo : hi;
  const double final_h = (use_lower ? lower : upper).liquid_enthalpy_j_per_kg;
  const auto property = LiquidPropertyService(&a.pack, 1).evaluate(
      {a.pack.material_fingerprint, t});
  if (!property.succeeded()) return out;
  const double rho = property.properties.density_kg_per_m3;
  constexpr double pi = 3.141592653589793238462643383279502884;
  const double d = std::cbrt(mass/rho)*std::cbrt(6./pi);
  const double geometric_mass = (pi/6.)*rho*d*d*d;
  const double error = (geometric_mass-mass)/mass;
  const double eps = std::numeric_limits<double>::epsilon();
  if (!std::isfinite(d) || !(d > 0) || !std::isfinite(error) ||
      std::abs(error) > 64*eps ||
      std::abs(final_h-h) > 64*eps*std::max({1.,std::abs(h),std::abs(final_h)}))
    return out;
  out.status = portable::Status::success;
  out.available = true;
  out.temperature_k = t;
  out.droplet_diameter_m = d;
  out.density_kg_per_m3 = rho;
  out.specific_enthalpy_j_per_kg = final_h;
  out.enthalpy_residual_j_per_kg = final_h-h;
  out.geometry_relative_mass_error = error;
  return out;
}
} // namespace hundun::v04::spray::detail
