// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include <algorithm>
#include <array>
#include <cmath>

namespace hundun::v04::chemistry::detail {
// Fixed semantic order: H2, H2O, CO, CO2, O2, N2, C12H23. The asset adapter
// maps its species before calling this kernel. Mole amounts and their rates
// have units kmol/kg and kmol/(kg s); Gibbs inputs are dimensionless g/(RT).
// Preparing on an admitted state and evaluating on an internal trial are
// separate operations, matching the reference's coefficient lifetime.
class KeroseneFourStep final {
public:
  using Tuple = std::array<double, 7>;
  bool prepare(double temperature_k, double gas_constant_j_per_kmol_k,
               double fuel_molecular_weight, const Tuple &accepted_moles,
               const Tuple &gibbs_rt) noexcept {
    if (!(temperature_k > 0) || !std::isfinite(temperature_k) ||
        !(gas_constant_j_per_kmol_k > 0) ||
        !std::isfinite(gas_constant_j_per_kmol_k) ||
        !(fuel_molecular_weight > 0) || !std::isfinite(fuel_molecular_weight))
      return false;
    double sum = 0;
    for (unsigned i = 0; i < 7; ++i) {
      if (!(accepted_moles[i] >= 0) || !std::isfinite(accepted_moles[i]) ||
          !std::isfinite(gibbs_rt[i])) return false;
      sum += accepted_moles[i];
    }
    KeroseneFourStep next;
    next.threshold_ = 1e-5 * sum;
    const double fuel_equivalent =
        (accepted_moles[2] + accepted_moles[3] + 12 * accepted_moles[6]) *
        fuel_molecular_weight / 12;
    const double correction = std::clamp(1 - 20 * (fuel_equivalent - .07), 0., 1.);
    const double rt = gas_constant_j_per_kmol_k / 4184 * temperature_k;
    next.k1_ = correction * (1e15 * std::exp(-40000 / rt) / temperature_k);
    next.k2_ = 2.75e9 * std::exp(-20000 / rt);
    next.k3_ = 3e8 * std::exp(-30000 / rt);
    next.k4_ = correction * (8.5e9 * std::exp(-30000 / rt));
    next.eq1_ = std::exp(gibbs_rt[0] + .5 * gibbs_rt[4] - gibbs_rt[1] +
        .5 * std::log(gas_constant_j_per_kmol_k * temperature_k / 1.013e5));
    const double eq2 = std::exp(gibbs_rt[2] + gibbs_rt[1] - gibbs_rt[3] - gibbs_rt[0]);
    next.kb2_ = next.k2_ / eq2;
    if (!(next.threshold_ > 0) || !std::isfinite(next.threshold_) ||
        !(next.eq1_ > 0) || !std::isfinite(next.eq1_) ||
        !(eq2 > 0) || !std::isfinite(eq2)) return false;
    for (const double value : {next.k1_, next.k2_, next.k3_, next.k4_, next.kb2_})
      if (!std::isfinite(value)) return false;
    next.prepared_ = true;
    *this = next;
    return true;
  }

  bool evaluate(double density_kg_per_m3, const Tuple &trial_moles,
                Tuple &rates_kmol_per_kg_s) const noexcept {
    if (!prepared_ || !(density_kg_per_m3 > 0) ||
        !std::isfinite(density_kg_per_m3)) return false;
    for (double value : trial_moles) if (!std::isfinite(value)) return false;
    const auto &f = trial_moles;
    constexpr double power = .9;
    const double oxygen = std::max(f[4], 0.);
    const double rho_forward = std::pow(density_kg_per_m3, power - .5);
    const double rho_backward = std::pow(density_kg_per_m3, power - 1);
    const double forward = k1_ * rho_forward * std::pow(oxygen, power + 1);
    const double backward = k1_ / eq1_ * rho_backward * std::pow(oxygen, power + .5);
    double r1, rb1;
    if (f[1] < threshold_ && f[0] < threshold_) {
      const double denominator = std::pow(threshold_, 1.5);
      r1 = forward * f[0] / denominator;
      rb1 = backward * f[1] / denominator;
    } else if (f[0] < threshold_) {
      r1 = forward * f[0] / (std::sqrt(threshold_) * f[1]);
      rb1 = backward / std::sqrt(threshold_);
    } else if (f[1] < threshold_) {
      r1 = forward * std::sqrt(f[0]) / threshold_;
      rb1 = backward * f[1] / (std::sqrt(f[0]) * threshold_);
    } else {
      r1 = forward * std::sqrt(f[0]) / f[1];
      rb1 = backward / std::sqrt(f[0]);
    }
    const double r2 = k2_ * density_kg_per_m3 * f[2] * f[1];
    const double rb2 = kb2_ * density_kg_per_m3 * f[3] * f[0];
    const double r3 = k3_ * density_kg_per_m3 * f[6] * f[1];
    const double r4 = f[6] < threshold_
        ? k4_ * rho_forward * f[6] * std::pow(oxygen, power) / std::sqrt(threshold_)
        : k4_ * rho_forward * std::sqrt(f[6]) * std::pow(oxygen, power);
    const Tuple next{rb1 - r1 + r2 - rb2 + 23.5 * r3 + 11.5 * r4,
                     r1 - rb1 + rb2 - r2 - 12 * r3,
                     12 * (r3 + r4) + rb2 - r2, r2 - rb2,
                     .5 * (rb1 - r1) - 6 * r4, 0., -(r3 + r4)};
    for (double value : next) if (!std::isfinite(value)) return false;
    rates_kmol_per_kg_s = next;
    return true;
  }

private:
  double threshold_{}, k1_{}, k2_{}, k3_{}, k4_{}, kb2_{}, eq1_{};
  bool prepared_{};
};
} // namespace hundun::v04::chemistry::detail
