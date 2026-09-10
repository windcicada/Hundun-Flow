// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/v04_flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace hundun::v04::detail {

// Four previous fixed-point directions improve only an uncommitted species
// guess. The caller still solves and certifies the original final equations.
class SpeciesCouplingHistory {
public:
  static constexpr std::size_t reduction_capacity = 21U;
  Status reserve(std::size_t count) {
    if (count > d_.max_size() / 5U)
      return {StatusCode::invalid_plan, 10214U};
    stride_ = count;
    d_.resize(5U * count);
    g_.resize(5U * count);
    guess_.resize(count);
    return {};
  }
  void begin(unsigned sweep) noexcept {
    if (sweep == 1U || sweep != sweep_ + 1U)
      reset();
    sweep_ = sweep;
  }
  bool valid() const noexcept { return valid_; }
  double guess(std::size_t index) const noexcept { return guess_[index]; }
  std::uint64_t owned_payload_bytes() const noexcept {
    return sizeof(double) * (d_.capacity() + g_.capacity() + guess_.capacity());
  }

  // Gram and rhs layout: 4x4 Gram, four rhs values, current residual norm^2.
  // Column normalization makes rank decisions independent of direction size.
  static bool coefficients(const std::array<double, 21U> &gram, unsigned m,
                           std::array<double, 4U> &gamma) noexcept {
    gamma = {};
    if (m == 0U || m > 4U || !std::isfinite(gram[20]) || gram[20] <= 0.0)
      return false;
    double matrix[4][5]{}, column[4]{};
    const double residual = std::sqrt(gram[20]);
    for (unsigned a = 0U; a < m; ++a) {
      if (!std::isfinite(gram[4U * a + a]) || gram[4U * a + a] <= 0.0)
        return false;
      column[a] = std::sqrt(gram[4U * a + a]);
    }
    for (unsigned a = 0U; a < m; ++a) {
      for (unsigned b = 0U; b < m; ++b) {
        matrix[a][b] = gram[4U * a + b] / (column[a] * column[b]);
        if (!std::isfinite(matrix[a][b]))
          return false;
      }
      matrix[a][a] += 1e-14;
      matrix[a][m] = gram[16U + a] / (column[a] * residual);
      if (!std::isfinite(matrix[a][m]))
        return false;
    }
    for (unsigned a = 0U; a < m; ++a) {
      unsigned pivot = a;
      for (unsigned b = a + 1U; b < m; ++b)
        if (std::abs(matrix[b][a]) > std::abs(matrix[pivot][a]))
          pivot = b;
      for (unsigned b = a; b <= m; ++b)
        std::swap(matrix[a][b], matrix[pivot][b]);
      if (std::abs(matrix[a][a]) < 1e-13)
        return false;
      for (unsigned b = a + 1U; b < m; ++b) {
        const double factor = matrix[b][a] / matrix[a][a];
        for (unsigned c = a; c <= m; ++c)
          matrix[b][c] -= factor * matrix[a][c];
      }
    }
    for (int a = static_cast<int>(m) - 1; a >= 0; --a) {
      double value = matrix[a][m];
      for (unsigned b = a + 1U; b < m; ++b)
        value -= matrix[a][b] * gamma[b];
      gamma[a] = value / matrix[a][a];
    }
    double length = 0.0;
    for (unsigned a = 0U; a < m; ++a) {
      gamma[a] *= residual / column[a];
      length += std::abs(gamma[a]);
    }
    if (!std::isfinite(length) || length > 4.0)
      return false;
    double estimate = gram[20];
    for (unsigned a = 0U; a < m; ++a) {
      estimate -= 2.0 * gamma[a] * gram[16U + a];
      for (unsigned b = 0U; b < m; ++b)
        estimate += gamma[a] * gamma[b] * gram[4U * a + b];
    }
    return std::isfinite(estimate) && estimate >= -1e-12 * gram[20] &&
           estimate < 0.99 * gram[20];
  }

  Status update(Int3 cells, Span<const FieldView> input,
                Span<const FieldView> solved, Span<const std::size_t> species,
                Span<const std::uint8_t> activity, double residual,
                ReductionEngine &reductions, bool &accelerated) noexcept {
    accelerated = false;
    if (stride_ == 0U)
      return {};
    valid_ = false;
    if (!std::isfinite(residual) || residual < 0.0)
      return {StatusCode::rejected_step, 10214U};
    if (count_ && residual >= 0.9 * previous_)
      reset();
    previous_ = residual;
    const std::size_t cells_count =
        static_cast<std::size_t>(cells.x) * cells.y * cells.z;
    const unsigned current = count_ % 5U, m = std::min(count_, 4U);
    auto *d = d_.data() + current * stride_;
    auto *g = g_.data() + current * stride_;
    for (std::size_t s = 0U; s < species.size; ++s) {
      const auto slot = species.data[s];
      std::size_t i = 0U;
      for (int z = 0; z < cells.z; ++z)
        for (int y = 0; y < cells.y; ++y)
          for (int x = 0; x < cells.x; ++x, ++i) {
            const auto j = slot * cells_count + i;
            g[j] = solved.data[slot].unchecked({x, y, z}, 0U);
            d[j] = g[j] - input.data[slot].unchecked({x, y, z}, 0U);
            guess_[j] = g[j];
          }
    }
    std::array<double, 21U> local{}, global{};
    for (std::size_t s = 0U; s < species.size; ++s)
      for (std::size_t i = 0U; i < cells_count; ++i) {
        if (activity.size && activity.data[i] == 0U)
          continue;
        const auto j = species.data[s] * cells_count + i;
        double diff[4]{};
        for (unsigned a = 0U; a < m; ++a)
          diff[a] = d[j] - d_[((count_ - a - 1U) % 5U) * stride_ + j];
        for (unsigned a = 0U; a < m; ++a) {
          for (unsigned b = 0U; b < m; ++b)
            local[4U * a + b] += diff[a] * diff[b];
          local[16U + a] += diff[a] * d[j];
        }
        local[20] += d[j] * d[j];
      }
    auto status = reductions.checked_sum({local.data(), local.size()},
                                         {global.data(), global.size()});
    if (!status) {
      reset();
      return status;
    }
    std::array<double, 4U> gamma{};
    accelerated = coefficients(global, m, gamma);
    if (accelerated)
      for (std::size_t i = 0U; i < cells_count; ++i) {
        if (activity.size && activity.data[i] == 0U)
          continue;
        long double sum = 0.0L;
        bool admissible = true;
        for (std::size_t s = 0U; s < species.size; ++s) {
          const auto j = species.data[s] * cells_count + i;
          double value = g[j];
          for (unsigned a = 0U; a < m; ++a)
            value -=
                gamma[a] * (g[j] - g_[((count_ - a - 1U) % 5U) * stride_ + j]);
          guess_[j] = value;
          sum += value;
          admissible &= std::isfinite(value) && value >= 0.0 && value <= 1.0;
        }
        admissible &= sum <= 1.0L;
        if (!admissible)
          for (std::size_t s = 0U; s < species.size; ++s) {
            const auto j = species.data[s] * cells_count + i;
            guess_[j] = g[j];
          }
      }
    valid_ = true;
    ++count_;
    return {};
  }

private:
  void reset() noexcept {
    count_ = 0U;
    valid_ = false;
    previous_ = 0.0;
  }
  std::vector<double> d_, g_, guess_;
  std::size_t stride_{};
  unsigned sweep_{}, count_{};
  double previous_{};
  bool valid_{};
};
} // namespace hundun::v04::detail
