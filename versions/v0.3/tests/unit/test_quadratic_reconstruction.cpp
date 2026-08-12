// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_error.hpp"
#include "src/ib_deterministic_qr_detail.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::detail::DeterministicQr;
using hundun::immersed::detail::factorize_design_matrix;
using hundun::runtime::Error;

template <class Function> std::string expect_error(Function &&function) {
  try {
    function();
  } catch (const Error &error) {
    return error.what();
  }
  throw std::runtime_error("expected hundun::runtime::Error");
}

std::vector<double> full_rank_matrix(std::size_t rows) {
  HUNDUN_CHECK(rows >= 10U);
  std::vector<double> matrix(rows * 10U, 0.0);
  for (std::size_t diagonal = 0; diagonal < 10U; ++diagonal) {
    matrix[diagonal * 10U + diagonal] = 1.0;
  }
  for (std::size_t row = 10U; row < rows; ++row) {
    for (std::size_t column = 0; column < 10U; ++column) {
      matrix[row * 10U + column] =
          static_cast<double>((row + 1U) * (column + 3U) % 11U) / 32.0;
    }
  }
  return matrix;
}

std::vector<double> full_rank_matrix(std::size_t rows, std::size_t columns) {
  HUNDUN_CHECK(columns > 0U && columns <= 10U && rows >= columns);
  std::vector<double> matrix(rows * columns, 0.0);
  for (std::size_t diagonal = 0; diagonal < columns; ++diagonal)
    matrix[diagonal * columns + diagonal] = 1.0;
  for (std::size_t row = columns; row < rows; ++row)
    for (std::size_t column = 0; column < columns; ++column)
      matrix[row * columns + column] =
          static_cast<double>((row + 1U) * (column + 3U) % 11U) / 32.0;
  return matrix;
}

void test_rank_condition_and_functional_weights() {
  const auto matrix = full_rank_matrix(14U);
  const DeterministicQr qr = factorize_design_matrix(matrix, 14U, 10U);
  HUNDUN_CHECK(qr.rank == 10U);
  HUNDUN_CHECK(std::isfinite(qr.condition_estimate));
  HUNDUN_CHECK(qr.condition_estimate > 0.0);
  HUNDUN_CHECK(qr.pivots.size() == 10U);
  HUNDUN_CHECK(qr.pivot_fingerprint != 0U);

  for (std::size_t selected_functional = 0; selected_functional < 10U;
       ++selected_functional) {
    std::vector<double> functional(10U, 0.0);
    functional[selected_functional] = 1.0;
    const auto weights = qr.functional_weights(functional);
    HUNDUN_CHECK(weights.size() == 14U);
    double weight_l1 = 0.0;
    for (const double weight : weights) {
      weight_l1 += std::abs(weight);
    }
    const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, weight_l1);
    for (std::size_t column = 0; column < 10U; ++column) {
      double observed = 0.0;
      for (std::size_t row = 0; row < 14U; ++row) {
        observed += matrix[row * 10U + column] * weights[row];
      }
      HUNDUN_CHECK_NEAR(observed, functional[column], tolerance);
    }
  }
}

void test_rank_and_input_rejection() {
  auto rank_nine = full_rank_matrix(14U);
  for (std::size_t row = 0; row < 14U; ++row) {
    rank_nine[row * 10U + 9U] = rank_nine[row * 10U + 8U];
  }
  HUNDUN_CHECK(expect_error([&] {
                 factorize_design_matrix(rank_nine, 14U, 10U);
               }).find("rank") != std::string::npos);
  HUNDUN_CHECK(expect_error([&] {
                 factorize_design_matrix(full_rank_matrix(13U), 14U, 10U);
               }).find("size") != std::string::npos);
  auto nonfinite = full_rank_matrix(14U);
  nonfinite[3] = std::numeric_limits<double>::quiet_NaN();
  HUNDUN_CHECK(expect_error([&] {
                 factorize_design_matrix(nonfinite, 14U, 10U);
               }).find("finite") != std::string::npos);
}

void test_condition_threshold() {
  std::vector<double> below(100U, 0.0);
  std::vector<double> above(100U, 0.0);
  for (std::size_t diagonal = 0; diagonal < 10U; ++diagonal) {
    below[diagonal * 10U + diagonal] = 1.0;
    above[diagonal * 10U + diagonal] = 1.0;
  }
  below[9U * 10U + 9U] = 1.0 / (1.0e8 * (1.0 - 1.0e-6));
  above[9U * 10U + 9U] = 1.0 / (1.0e8 * (1.0 + 1.0e-6));
  const auto accepted = factorize_design_matrix(below, 10U, 10U);
  HUNDUN_CHECK(accepted.condition_estimate < 1.0e8);
  HUNDUN_CHECK(expect_error([&] {
                 factorize_design_matrix(above, 10U, 10U);
               }).find("condition") != std::string::npos);
}

void test_pivot_tie_and_row_order_identity() {
  const auto matrix = full_rank_matrix(10U);
  const auto first = factorize_design_matrix(matrix, 10U, 10U);
  HUNDUN_CHECK(first.pivots == std::vector<std::uint32_t>(
                                   {0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U}));

  std::vector<double> reversed(matrix.size());
  for (std::size_t row = 0; row < 10U; ++row) {
    for (std::size_t column = 0; column < 10U; ++column) {
      reversed[(9U - row) * 10U + column] = matrix[row * 10U + column];
    }
  }
  const auto second = factorize_design_matrix(reversed, 10U, 10U);
  HUNDUN_CHECK(first.pivots == second.pivots);
  HUNDUN_CHECK(first.pivot_fingerprint == second.pivot_fingerprint);

  auto inside_tie = matrix;
  inside_tie[11U] = 1.0 + 8.0 * std::numeric_limits<double>::epsilon();
  HUNDUN_CHECK(factorize_design_matrix(inside_tie, 10U, 10U).pivots.front() ==
               0U);
  auto outside_tie = matrix;
  outside_tie[11U] = 1.0 + 128.0 * std::numeric_limits<double>::epsilon();
  HUNDUN_CHECK(factorize_design_matrix(outside_tie, 10U, 10U).pivots.front() ==
               1U);
}

void test_constrained_quadratic_free_subspace_dimensions() {
  for (std::size_t columns = 4U; columns <= 9U; ++columns) {
    const std::size_t rows = columns + 5U;
    const auto matrix = full_rank_matrix(rows, columns);
    const auto qr = factorize_design_matrix(matrix, rows, columns);
    HUNDUN_CHECK(qr.rank == columns);
    for (std::size_t selected = 0U; selected < columns; ++selected) {
      std::vector<double> functional(columns, 0.0);
      functional[selected] = 1.0;
      const auto weights = qr.functional_weights(functional);
      HUNDUN_CHECK(weights.size() == rows);
      for (std::size_t column = 0U; column < columns; ++column) {
        double observed = 0.0;
        for (std::size_t row = 0U; row < rows; ++row)
          observed += matrix[row * columns + column] * weights[row];
        HUNDUN_CHECK_NEAR(
            observed, functional[column],
            256.0 * std::numeric_limits<double>::epsilon());
      }
    }
  }
}

void run_tests() {
  test_rank_condition_and_functional_weights();
  test_rank_and_input_rejection();
  test_condition_threshold();
  test_pivot_tie_and_row_order_identity();
  test_constrained_quadratic_free_subspace_dimensions();
}

} // namespace

int main() { return hundun::test::run(run_tests); }
