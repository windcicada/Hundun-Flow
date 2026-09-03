// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04::detail::ibm_qr {

inline constexpr std::size_t kMaximumRows = 32U;
inline constexpr std::size_t kMaximumColumns = 10U;

struct Factorization {
  std::array<double, kMaximumRows * kMaximumColumns> qr{};
  std::array<double, kMaximumColumns> tau{};
  std::array<std::uint8_t, kMaximumColumns> permutation{};
  std::size_t rows{};
  std::size_t columns{};
  std::uint8_t rank{};
  double condition_estimate{};
};

inline double& entry(Factorization& factor, std::size_t row,
                     std::size_t column) noexcept {
  return factor.qr[row * kMaximumColumns + column];
}

inline double entry(const Factorization& factor, std::size_t row,
                    std::size_t column) noexcept {
  return factor.qr[row * kMaximumColumns + column];
}

inline double stable_norm(const Factorization& factor, std::size_t begin,
                          std::size_t column) noexcept {
  double scale = 0.0;
  double sum_squares = 1.0;
  for (std::size_t row = begin; row < factor.rows; ++row) {
    const double magnitude = std::abs(entry(factor, row, column));
    if (magnitude == 0.0) {
      continue;
    }
    if (scale < magnitude) {
      const double ratio = scale / magnitude;
      sum_squares = 1.0 + sum_squares * ratio * ratio;
      scale = magnitude;
    } else {
      const double ratio = magnitude / scale;
      sum_squares += ratio * ratio;
    }
  }
  return scale == 0.0 ? 0.0 : scale * std::sqrt(sum_squares);
}

inline bool factorize(const double* matrix, std::size_t rows,
                      std::size_t columns, double condition_limit,
                      Factorization& out) noexcept {
  if (matrix == nullptr || rows < columns || rows > kMaximumRows ||
      columns == 0U || columns > kMaximumColumns ||
      !std::isfinite(condition_limit) || condition_limit < 1.0) {
    return false;
  }

  Factorization candidate;
  candidate.rows = rows;
  candidate.columns = columns;
  for (std::size_t row = 0U; row < rows; ++row) {
    for (std::size_t column = 0U; column < columns; ++column) {
      const double value = matrix[row * columns + column];
      if (!std::isfinite(value)) {
        return false;
      }
      entry(candidate, row, column) = value;
    }
  }
  for (std::size_t column = 0U; column < columns; ++column) {
    candidate.permutation[column] = static_cast<std::uint8_t>(column);
  }

  constexpr double kTieMultiplier = 64.0;
  const double epsilon = std::numeric_limits<double>::epsilon();
  for (std::size_t pivot = 0U; pivot < columns; ++pivot) {
    std::array<double, kMaximumColumns> norms{};
    double maximum_norm = 0.0;
    for (std::size_t column = pivot; column < columns; ++column) {
      norms[column] = stable_norm(candidate, pivot, column);
      maximum_norm = std::max(maximum_norm, norms[column]);
    }
    const double tie_tolerance =
        kTieMultiplier * epsilon * std::max(1.0, maximum_norm);
    std::size_t selected = pivot;
    for (std::size_t column = pivot + 1U; column < columns; ++column) {
      if (norms[column] > norms[selected] + tie_tolerance ||
          (std::abs(norms[column] - norms[selected]) <= tie_tolerance &&
           candidate.permutation[column] <
               candidate.permutation[selected])) {
        selected = column;
      }
    }
    if (selected != pivot) {
      for (std::size_t row = 0U; row < rows; ++row) {
        std::swap(entry(candidate, row, pivot),
                  entry(candidate, row, selected));
      }
      std::swap(candidate.permutation[pivot],
                candidate.permutation[selected]);
    }

    const double norm = stable_norm(candidate, pivot, pivot);
    const double alpha = entry(candidate, pivot, pivot);
    if (!std::isfinite(norm) || norm == 0.0) {
      candidate.tau[pivot] = 0.0;
      entry(candidate, pivot, pivot) = 0.0;
      continue;
    }
    const double beta = std::copysign(norm, alpha == 0.0 ? -1.0 : -alpha);
    const double denominator = alpha - beta;
    if (!std::isfinite(beta) || denominator == 0.0) {
      return false;
    }
    const double tau = (beta - alpha) / beta;
    candidate.tau[pivot] = tau;
    entry(candidate, pivot, pivot) = beta;
    for (std::size_t row = pivot + 1U; row < rows; ++row) {
      entry(candidate, row, pivot) /= denominator;
    }

    for (std::size_t column = pivot + 1U; column < columns; ++column) {
      double projection = entry(candidate, pivot, column);
      for (std::size_t row = pivot + 1U; row < rows; ++row) {
        projection +=
            entry(candidate, row, pivot) * entry(candidate, row, column);
      }
      projection *= tau;
      entry(candidate, pivot, column) -= projection;
      for (std::size_t row = pivot + 1U; row < rows; ++row) {
        entry(candidate, row, column) -=
            entry(candidate, row, pivot) * projection;
      }
    }
  }

  double maximum_diagonal = 0.0;
  for (std::size_t diagonal = 0U; diagonal < columns; ++diagonal) {
    maximum_diagonal =
        std::max(maximum_diagonal,
                 std::abs(entry(candidate, diagonal, diagonal)));
  }
  const double rank_tolerance =
      kTieMultiplier * epsilon * std::max(1.0, maximum_diagonal);
  double minimum_accepted_diagonal = std::numeric_limits<double>::infinity();
  for (std::size_t diagonal = 0U; diagonal < columns; ++diagonal) {
    const double magnitude =
        std::abs(entry(candidate, diagonal, diagonal));
    if (std::isfinite(magnitude) && magnitude > rank_tolerance) {
      ++candidate.rank;
      minimum_accepted_diagonal =
          std::min(minimum_accepted_diagonal, magnitude);
    }
  }
  if (candidate.rank != columns || maximum_diagonal == 0.0 ||
      !std::isfinite(minimum_accepted_diagonal)) {
    candidate.condition_estimate =
        std::numeric_limits<double>::infinity();
    out = candidate;
    return false;
  }
  candidate.condition_estimate =
      maximum_diagonal / minimum_accepted_diagonal;
  if (!std::isfinite(candidate.condition_estimate) ||
      candidate.condition_estimate > condition_limit) {
    out = candidate;
    return false;
  }
  out = candidate;
  return true;
}

inline bool solve_transpose_minimum_norm(
    const Factorization& factor, const double* functional,
    double* weights) noexcept {
  if (functional == nullptr || weights == nullptr || factor.rows == 0U ||
      factor.columns == 0U || factor.rank != factor.columns ||
      factor.rows > kMaximumRows || factor.columns > kMaximumColumns) {
    return false;
  }
  std::array<double, kMaximumRows> transformed{};
  for (std::size_t pivot = 0U; pivot < factor.columns; ++pivot) {
    const std::size_t original = factor.permutation[pivot];
    double value = functional[original];
    if (!std::isfinite(value)) {
      return false;
    }
    for (std::size_t previous = 0U; previous < pivot; ++previous) {
      value -= entry(factor, previous, pivot) * transformed[previous];
    }
    const double diagonal = entry(factor, pivot, pivot);
    if (!std::isfinite(value) || !std::isfinite(diagonal) ||
        diagonal == 0.0) {
      return false;
    }
    transformed[pivot] = value / diagonal;
  }

  for (std::size_t reverse = factor.columns; reverse-- > 0U;) {
    double projection = transformed[reverse];
    for (std::size_t row = reverse + 1U; row < factor.rows; ++row) {
      projection += entry(factor, row, reverse) * transformed[row];
    }
    projection *= factor.tau[reverse];
    transformed[reverse] -= projection;
    for (std::size_t row = reverse + 1U; row < factor.rows; ++row) {
      transformed[row] -= entry(factor, row, reverse) * projection;
    }
  }

  for (std::size_t row = 0U; row < factor.rows; ++row) {
    if (!std::isfinite(transformed[row])) {
      return false;
    }
    weights[row] = transformed[row];
  }
  return true;
}

}  // namespace hundun::v04::detail::ibm_qr
