// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_tcr_dynamic_detail.hpp"

namespace hundun::v04::tcr::detail {
// dyn711 statistics calls have their own clocks. Eight accumulation calls
// precede a separate evaluation call; Cphi updates on calls 1, 2, 9, 16, ... .
// These values belong to accepted model history, not an outer-iteration clock.
struct Dyn711Clock {
  unsigned rate_intervals{}, cphi_count{};
  bool rates_initialized{};
};
struct Dyn711Tick {
  bool available{}, evaluate_rates{}, update_cphi{};
  Dyn711Clock next{};
};
Dyn711Tick dyn711_tick(Dyn711Clock accepted) noexcept;

struct Dyn711Control {
  bool available{}, upper_branch{}, linear_endpoint{};
  double ratio{}, selected{}, effective{};
};
Dyn711Control dyn711_species_control(double eta, double pdf_sum, double psr_sum,
    double chemical_time, double flow_time, double weak_sum) noexcept;

struct Dyn711RateState {
  double pdf_sum{}, psr_sum{}, selected{.2}, effective{.2};
  bool upper_branch{}, linear_endpoint{};
};
// Candidate-only operation: every retry starts from the same accepted window.
// The source's 1e5 scale is retained so weak_sum has its original units.
// On an evaluation call, that call's rates are deliberately outside the window.
bool dyn711_advance_rate(const Dyn711RateState &accepted, Dyn711Clock clock,
    double dt, double pdf_rate, double psr_rate, double eta,
    double chemical_time, double flow_time, double weak_sum,
    Dyn711RateState &candidate) noexcept;

// Unnormalized products preserve the relative weights of the second spatial
// filter. This is a separate identity from the 624CF bounded-donor products.
struct Dyn711FilterDonor {
  double density{}, volume{}, scalar{};
  std::array<double, 3> gradient{};
};
// All coordinates are scalar-specific. For the mixture-fraction channel they
// are dimensionless; species molar/mass conversion is an explicit caller duty.
bool dyn711_filter_moments(const Dyn711FilterDonor *donors, unsigned count,
    DynamicFilterMoments &moments) noexcept;
DynamicFilterProducts dyn711_filter_products(const DynamicFilterMoments &) noexcept;
double dyn711_filter_ratio(double m_squared, double l_times_m) noexcept;
// Ordered input: mixture fraction, fuel, OH. The source prefers an admissible
// OH coefficient, then fuel, then mixture fraction, else the nearest to four.
double dyn711_select_cphi(const std::array<double, 3> &) noexcept;
double dyn711_smooth_cphi(double centre, const std::array<double, 7> &neighbors) noexcept;
} // namespace hundun::v04::tcr::detail
