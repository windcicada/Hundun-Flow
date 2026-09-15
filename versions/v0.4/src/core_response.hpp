// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/v04_status.hpp"
#include "hundun/v04_types.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <vector>

namespace hundun::v04::detail {

// An interval reaction increment is common to the coupled equation updates.
// Each proposed reuse is checked against a freshly integrated endpoint. This
// budget concerns outer coupling, separately from the integrator's tolerances.
class IntervalReactionResponse {
public:
  struct Report {
    bool reused{};
    double error_ratio{};
  };

  Status prepare(std::size_t cells, std::size_t species, std::size_t dependent,
                 double relative, double absolute) noexcept {
    if (!cells || species < 2 || dependent >= species ||
        !std::isfinite(relative) || relative <= 0 || relative >= 1 ||
        !std::isfinite(absolute) || absolute < 0 ||
        species > (SIZE_MAX - 1) / sizeof(double) ||
        cells > SIZE_MAX / (species * sizeof(double) + 1) ||
        cells > valid_.max_size() || cells * species > delta_.max_size())
      return {StatusCode::invalid_plan, 10300};
    try {
      delta_.resize(cells * species);
      valid_.assign(cells, 0);
    } catch (const std::bad_alloc&) {
      return {StatusCode::allocation_failure, 10300};
    }
    species_ = species;
    dependent_ = dependent;
    relative_ = relative;
    absolute_ = absolute;
    return {};
  }

  void reset() noexcept { std::fill(valid_.begin(), valid_.end(), 0); }
  std::size_t bytes() const noexcept {
    return sizeof(double) * delta_.capacity() + valid_.capacity();
  }

  Status select(std::size_t cell, Span<const double> initial,
                Span<double> endpoint, Report& report,
                double interval_scale = 1.0) noexcept {
    report = {};
    if (cell >= valid_.size() || initial.size != species_ ||
        endpoint.size != species_ || !initial.data || !endpoint.data ||
        !std::isfinite(interval_scale) || interval_scale < 0 || interval_scale > 1)
      return {StatusCode::invalid_plan, 10300};
    for (std::size_t j = 0; j < species_; ++j)
      if (!std::isfinite(initial.data[j]) || initial.data[j] < 0 ||
          initial.data[j] > 1 || !std::isfinite(endpoint.data[j]) ||
          endpoint.data[j] < 0 || endpoint.data[j] > 1)
        return {StatusCode::numerical_failure, 10301};
    double* delta = delta_.data() + cell * species_;
    if (valid_[cell]) {
      bool reusable = true;
      double maximum{};
      long double independent_sum{};
      for (std::size_t j = 0; j < species_; ++j) {
        const double retained = initial.data[j] + delta[j];
        if (!std::isfinite(retained) || retained < 0 || retained > 1) {
          reusable = false;
          break;
        }
        if (j != dependent_) independent_sum += retained;
        const double error = std::abs(retained - endpoint.data[j]);
        const double limit = interval_scale * (absolute_ + relative_ *
            std::max({std::abs(initial.data[j]), std::abs(endpoint.data[j]),
                      std::abs(retained)}));
        if (error > limit) {
          reusable = false;
          break;
        }
        maximum = std::max(maximum, error == 0 ? 0.0 : error / limit);
      }
      if (reusable && independent_sum <= 1.0L) {
        for (std::size_t j = 0; j < species_; ++j)
          endpoint.data[j] = initial.data[j] + delta[j];
        report = {true, maximum};
        return {};
      }
    }
    for (std::size_t j = 0; j < species_; ++j)
      delta[j] = endpoint.data[j] - initial.data[j];
    valid_[cell] = 1;
    return {};
  }

private:
  std::size_t species_{}, dependent_{};
  double relative_{}, absolute_{};
  std::vector<double> delta_;
  std::vector<std::uint8_t> valid_;
};
} // namespace hundun::v04::detail
