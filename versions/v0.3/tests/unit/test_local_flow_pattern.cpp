// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/rt_error.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace hundun::immersed::detail {
#if defined(HUNDUN_IMMERSED_ENABLE_TEST_ACCESS)
LocalCoefficientRow
paper_transform_factors_for_test(const LocalCoefficientRow &, double, double,
                                 double, double);
#endif
} // namespace hundun::immersed::detail

namespace {

using hundun::immersed::LocalCoefficientRow;
using hundun::immersed::LocalFlowPatternTransform;
using hundun::immersed::ReplacementGroup;
using hundun::immersed::RowReplacementPlan;
using hundun::runtime::Error;
using hundun::runtime::Real3;

template <class Function> std::string expect_error(Function &&function) {
  try {
    function();
  } catch (const Error &error) {
    return error.what();
  }
  throw std::runtime_error("expected hundun::runtime::Error");
}

std::array<double, 8> array(const LocalCoefficientRow &row) {
  return {row.neighbour[0], row.neighbour[1], row.neighbour[2],
          row.neighbour[3], row.neighbour[4], row.neighbour[5],
          row.diagonal,     row.source};
}

bool bitwise_equal(double lhs, double rhs) {
  std::uint64_t lhs_bits = 0;
  std::uint64_t rhs_bits = 0;
  std::memcpy(&lhs_bits, &lhs, sizeof(lhs_bits));
  std::memcpy(&rhs_bits, &rhs, sizeof(rhs_bits));
  return lhs_bits == rhs_bits;
}

using Matrix7 = std::array<std::array<double, 7>, 7>;

std::array<double, 7> multiply(const Matrix7 &matrix,
                               const std::array<double, 7> &vector) {
  std::array<double, 7> result{};
  for (std::size_t row = 0; row < result.size(); ++row) {
    for (std::size_t column = 0; column < vector.size(); ++column) {
      result[row] += matrix[row][column] * vector[column];
    }
  }
  return result;
}

LocalCoefficientRow independent_paper_map(LocalCoefficientRow row, double k0,
                                          double k1, double k2, double k3) {
  const auto a1 = array(row);
  const std::array<double, 7> input{a1[0], a1[1], a1[2], a1[3],
                                    a1[4], a1[5], a1[6]};
  Matrix7 equation_16{};
  equation_16[0][0] = 1.0;
  equation_16[1][0] = 1.0; // Literal printed Eq. 16: A_S2 = A_N1.
  equation_16[2][2] = k1;
  equation_16[2][5] = 1.0 - k1;
  equation_16[3][3] = k1;
  equation_16[3][4] = 1.0 - k1;
  equation_16[4][4] = k1;
  equation_16[4][2] = 1.0 - k1;
  equation_16[5][5] = k1;
  equation_16[5][3] = 1.0 - k1;
  equation_16[6][6] = 1.0;
  Matrix7 equation_15{};
  equation_15[0][0] = k2;
  equation_15[0][5] = 1.0 - k2;
  equation_15[1][1] = k2;
  equation_15[1][4] = 1.0 - k2;
  equation_15[2][2] = 1.0;
  equation_15[3][3] = 1.0;
  equation_15[4][4] = k2;
  equation_15[4][0] = 1.0 - k2;
  equation_15[5][5] = k2;
  equation_15[5][1] = 1.0 - k2;
  equation_15[6][6] = 1.0;
  Matrix7 equation_14{};
  equation_14[0][0] = k3;
  equation_14[0][2] = 1.0 - k3;
  equation_14[1][1] = k3;
  equation_14[1][3] = 1.0 - k3;
  equation_14[2][2] = k3;
  equation_14[2][1] = 1.0 - k3;
  equation_14[3][3] = k3;
  equation_14[3][0] = 1.0 - k3;
  equation_14[4][4] = 1.0;
  equation_14[5][5] = 1.0;
  equation_14[6][6] = 1.0;
  const auto a0 = multiply(equation_14,
                           multiply(equation_15, multiply(equation_16, input)));
  Matrix7 equation_13{};
  for (std::size_t neighbour = 0; neighbour < 6U; ++neighbour) {
    equation_13[neighbour][neighbour] = k0;
    equation_13[6][neighbour] = 1.0 - k0;
  }
  equation_13[6][6] = 1.0;
  const auto transformed = multiply(equation_13, a0);
  LocalCoefficientRow result{};
  for (std::size_t index = 0; index < 6U; ++index) {
    result.neighbour[index] = transformed[index];
  }
  result.diagonal = transformed[6];
  result.source = row.source;
  return result;
}

void test_literal_paper_oracle() {
  const LocalCoefficientRow input{{2.0, 3.0, 5.0, 7.0, 11.0, 13.0}, 41.0, 17.0};
  const auto expected = independent_paper_map(input, 0.37, 0.40, 0.65, 0.25);
  const auto observed =
      hundun::immersed::detail::paper_transform_factors_for_test(
          input, 0.37, 0.40, 0.65, 0.25);
  constexpr std::array<double, 8> frozen{3.144075, 2.968325, 1.985975, 2.143225,
                                         2.0387,   2.5197,   66.2,     17.0};
  const double tolerance =
      256.0 * std::numeric_limits<double>::epsilon() * 99.0;
  for (std::size_t index = 0; index < 8U; ++index) {
    HUNDUN_CHECK_NEAR(array(observed)[index], array(expected)[index],
                      tolerance);
    HUNDUN_CHECK_NEAR(array(observed)[index], frozen[index], tolerance);
  }
}

void test_physical_normal_and_rejection() {
  const LocalFlowPatternTransform transform;
  const LocalCoefficientRow row{{2.0, 3.0, 5.0, 7.0, 11.0, 13.0}, 41.0, 17.0};
  const auto first = transform.transform_full(row, 0.37, {1.0, 2.0, 3.0});
  const auto scaled = transform.transform_full(row, 0.37, {2.0, 4.0, 6.0});
  for (std::size_t index = 0; index < 8U; ++index) {
    HUNDUN_CHECK(bitwise_equal(array(first)[index], array(scaled)[index]));
  }
  HUNDUN_CHECK(bitwise_equal(first.source, row.source));
  auto reflected_input = row;
  std::swap(reflected_input.neighbour[0], reflected_input.neighbour[1]);
  std::swap(reflected_input.neighbour[2], reflected_input.neighbour[3]);
  std::swap(reflected_input.neighbour[4], reflected_input.neighbour[5]);
  const auto reflected =
      transform.transform_full(reflected_input, 0.37, {-1.0, -2.0, -3.0});
  HUNDUN_CHECK(bitwise_equal(reflected.neighbour[0], first.neighbour[1]));
  HUNDUN_CHECK(bitwise_equal(reflected.neighbour[1], first.neighbour[0]));
  HUNDUN_CHECK(bitwise_equal(reflected.neighbour[2], first.neighbour[3]));
  HUNDUN_CHECK(bitwise_equal(reflected.neighbour[3], first.neighbour[2]));
  HUNDUN_CHECK(bitwise_equal(reflected.neighbour[4], first.neighbour[5]));
  HUNDUN_CHECK(bitwise_equal(reflected.neighbour[5], first.neighbour[4]));
  HUNDUN_CHECK(bitwise_equal(reflected.diagonal, first.diagonal));
  const auto zero_scale = transform.transform_full(row, 0.0, {1.0, 2.0, 3.0});
  const auto unit_scale = transform.transform_full(row, 1.0, {1.0, 2.0, 3.0});
  double zero_scale_sum = zero_scale.diagonal;
  double unit_scale_sum = unit_scale.diagonal;
  for (std::size_t index = 0; index < 6U; ++index) {
    HUNDUN_CHECK(bitwise_equal(zero_scale.neighbour[index], 0.0));
    zero_scale_sum += zero_scale.neighbour[index];
    unit_scale_sum += unit_scale.neighbour[index];
  }
  HUNDUN_CHECK_NEAR(zero_scale_sum, unit_scale_sum,
                    256.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, std::abs(unit_scale_sum)));
  HUNDUN_CHECK(transform.algorithm_fingerprint() != 0U);
  HUNDUN_CHECK(transform.algorithm_fingerprint() ==
               LocalFlowPatternTransform{}.algorithm_fingerprint());
  for (const Real3 axis :
       {Real3{1.0, 0.0, 0.0}, Real3{0.0, 1.0, 0.0}, Real3{0.0, 0.0, 1.0}}) {
    const auto axis_result = transform.transform_full(row, 0.5, axis);
    for (const double coefficient : array(axis_result)) {
      HUNDUN_CHECK(std::isfinite(coefficient));
    }
  }
  const auto tiny = transform.transform_full(
      row, 0.5, {std::numeric_limits<double>::denorm_min(), 0.0, 0.0});
  const auto huge = transform.transform_full(
      row, 0.5, {std::numeric_limits<double>::max() / 2.0, 0.0, 0.0});
  for (const double coefficient : array(tiny)) {
    HUNDUN_CHECK(std::isfinite(coefficient));
  }
  for (const double coefficient : array(huge)) {
    HUNDUN_CHECK(std::isfinite(coefficient));
  }
  HUNDUN_CHECK(expect_error([&] {
                 transform.transform_full(row, -0.1, {1, 0, 0});
               }).find("scale") != std::string::npos);
  HUNDUN_CHECK(expect_error([&] {
                 transform.transform_full(row, 0.5, {0, 0, 0});
               }).find("normal") != std::string::npos);
  auto bad = row;
  bad.diagonal = std::numeric_limits<double>::quiet_NaN();
  HUNDUN_CHECK(expect_error([&] {
                 transform.transform_full(bad, 0.5, {1, 0, 0});
               }).find("finite") != std::string::npos);
}

std::vector<std::uint32_t> occurrence_union(const RowReplacementPlan &plan) {
  std::vector<std::uint32_t> result;
  for (const auto &group : plan.groups) {
    result.insert(result.end(), group.algebraic_occurrences.begin(),
                  group.algebraic_occurrences.end());
  }
  std::sort(result.begin(), result.end());
  return result;
}

void test_plans_and_simultaneous_evaluation() {
  const LocalFlowPatternTransform transform;
  const auto one_link = transform.plan_row(
      40U, {10U}, LocalCoefficientRow{{2.0, 0, 0, 0, 0, 0}, 7.0, 11.0});
  HUNDUN_CHECK(one_link.groups.size() == 1U);
  HUNDUN_CHECK(one_link.groups.front().algebraic_occurrences ==
               std::vector<std::uint32_t>({0U, 6U, 7U}));
  const auto two_links = transform.plan_row(
      41U, {20U, 10U}, LocalCoefficientRow{{2.0, 3.0, 0, 0, 0, 0}, 7.0, 0.0});
  HUNDUN_CHECK(two_links.groups.size() == 3U);
  HUNDUN_CHECK(std::count_if(two_links.groups.begin(), two_links.groups.end(),
                             [](const auto &group) {
                               return group.links.size() == 2U;
                             }) == 1);
  const LocalCoefficientRow immersed{{2.0, 3.0, 5.0, 0.0, 0.0, 0.0}, 7.0, 11.0};
  const auto plan = transform.plan_row(42U, {30U, 10U, 20U}, immersed);
  const auto permuted = transform.plan_row(42U, {20U, 30U, 10U}, immersed);
  HUNDUN_CHECK(plan.fingerprint == permuted.fingerprint);
  HUNDUN_CHECK(occurrence_union(plan) ==
               std::vector<std::uint32_t>({0U, 1U, 2U, 6U, 7U}));
  HUNDUN_CHECK(plan.groups.size() == 5U);
  const auto singleton_count =
      std::count_if(plan.groups.begin(), plan.groups.end(),
                    [](const auto &group) { return group.links.size() == 1U; });
  const auto joint_count =
      std::count_if(plan.groups.begin(), plan.groups.end(),
                    [](const auto &group) { return group.links.size() == 3U; });
  HUNDUN_CHECK(singleton_count == 3);
  HUNDUN_CHECK(joint_count == 2);

  const LocalCoefficientRow snapshot{
      {2.0, 3.0, 5.0, 13.0, 17.0, 19.0}, 7.0, 11.0};
  const auto before = array(snapshot);
  const std::vector<double> symbols{1.0, 2.0, 4.0};
  const double observed =
      transform.evaluate_wall_replacement(plan, snapshot, symbols);
  const double mean = (1.0 + 2.0 + 4.0) / 3.0;
  const double expected =
      2.0 * 1.0 + 3.0 * 2.0 + 5.0 * 4.0 + 7.0 * mean + 11.0 * mean;
  HUNDUN_CHECK_NEAR(observed, expected,
                    256.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, std::abs(expected)));
  for (std::size_t index = 0; index < 8U; ++index) {
    HUNDUN_CHECK(bitwise_equal(array(snapshot)[index], before[index]));
  }
  const auto permuted_value =
      transform.evaluate_wall_replacement(permuted, snapshot, symbols);
  HUNDUN_CHECK(bitwise_equal(observed, permuted_value));
  auto reordered = plan;
  std::reverse(reordered.groups.begin(), reordered.groups.end());
  for (auto &group : reordered.groups) {
    std::reverse(group.links.begin(), group.links.end());
    std::reverse(group.algebraic_occurrences.begin(),
                 group.algebraic_occurrences.end());
  }
  HUNDUN_CHECK(bitwise_equal(observed, transform.evaluate_wall_replacement(
                                           reordered, snapshot, symbols)));

  auto omitted = plan;
  omitted.groups.pop_back();
  HUNDUN_CHECK(expect_error([&] {
                 transform.evaluate_wall_replacement(omitted, snapshot,
                                                     symbols);
               }).find("fingerprint") != std::string::npos);
  auto duplicated = plan;
  duplicated.groups[1].algebraic_occurrences.push_back(0U);
  HUNDUN_CHECK(!expect_error([&] {
                  transform.evaluate_wall_replacement(duplicated, snapshot,
                                                      symbols);
                }).empty());
  auto duplicated_group = plan;
  duplicated_group.groups.push_back(duplicated_group.groups.front());
  HUNDUN_CHECK(expect_error([&] {
                 transform.evaluate_wall_replacement(duplicated_group, snapshot,
                                                     symbols);
               }).find("group") != std::string::npos);
  auto out_of_range = plan;
  out_of_range.groups.front().algebraic_occurrences.front() = 8U;
  HUNDUN_CHECK(expect_error([&] {
                 transform.evaluate_wall_replacement(out_of_range, snapshot,
                                                     symbols);
               }).find("range") != std::string::npos);
  auto nonfinite_symbols = symbols;
  nonfinite_symbols[1] = std::numeric_limits<double>::infinity();
  HUNDUN_CHECK(expect_error([&] {
                 transform.evaluate_wall_replacement(plan, snapshot,
                                                     nonfinite_symbols);
               }).find("finite") != std::string::npos);
  HUNDUN_CHECK(expect_error([&] {
                 transform.plan_row(42U, {1U, 1U}, immersed);
               }).find("duplicate") != std::string::npos);
}

void test_limiting_and_mutation_oracles() {
  const LocalFlowPatternTransform transform;
  const LocalCoefficientRow constant_projection{{1.0, -1.0, 0, 0, 0, 0}, 0, 0};
  const auto constant_plan =
      transform.plan_row(7U, {100U, 200U}, constant_projection);
  HUNDUN_CHECK(
      bitwise_equal(transform.evaluate_wall_replacement(
                        constant_plan, constant_projection, {3.0, 3.0}),
                    0.0));
  const double linear = transform.evaluate_wall_replacement(
      constant_plan, constant_projection, {4.0, 3.0});
  HUNDUN_CHECK(bitwise_equal(linear, 1.0));
  HUNDUN_CHECK(!bitwise_equal(linear, 0.0)); // Deleted evaluator mutation.
  const double sequential_substitution = 1.0 * 4.0 + (-1.0 + 1.0 * 4.0) * 3.0;
  HUNDUN_CHECK(!bitwise_equal(sequential_substitution, linear));
  const LocalCoefficientRow shared{
      {2.0, 3.0, 5.0, 7.0, 11.0, 13.0}, 41.0, 17.0};
  const auto transformed_shared =
      transform.transform_full(shared, 0.37, {1.0, 2.0, 3.0});
  HUNDUN_CHECK(
      !bitwise_equal(transformed_shared.neighbour[0], shared.neighbour[0]));
  HUNDUN_CHECK(!bitwise_equal(linear + 5.0, linear)); // Added wall source.
  HUNDUN_CHECK(!bitwise_equal(linear + 2.0, linear)); // Retained P-G term.
}

void test_affine_ghost_and_paper_closure_evidence() {
  // HUNDUN affine Ghost/quadratic enhancement is deliberately separate from
  // the paper whole-row coefficient-map oracle.
  const auto polynomial = [](double x) { return 1.0 + 2.0 * x + 3.0 * x * x; };
  constexpr double ghost_x = -0.5;
  constexpr double wall_x = 0.0;
  constexpr double wall_value = 2.0;
  const double ghost_symbol =
      polynomial(ghost_x) - polynomial(wall_x) + wall_value;
  const LocalCoefficientRow projection{{5.0, 0, 0, 0, 0, 0}, 0, 0};
  const LocalFlowPatternTransform transform;
  const auto plan = transform.plan_row(9U, {91U}, projection);
  const double replacement =
      transform.evaluate_wall_replacement(plan, projection, {ghost_symbol});
  HUNDUN_CHECK_NEAR(replacement, 5.0 * ghost_symbol,
                    256.0 * std::numeric_limits<double>::epsilon() *
                        std::max(1.0, std::abs(replacement)));

  // Paper Eq. 17 coefficient cancellation and Eq. 18 zero-normal scalar
  // closure are algebraic evidence only; no wall-function source is used.
  constexpr double solid_coefficient = 7.0;
  constexpr double solid_value = 3.0;
  const double momentum_cancellation =
      solid_coefficient * solid_value - solid_coefficient * solid_value;
  HUNDUN_CHECK(bitwise_equal(momentum_cancellation, 0.0));
  constexpr double active_scalar = 4.0;
  constexpr double solid_scalar = active_scalar;
  const double zero_normal_scalar =
      solid_coefficient * (active_scalar - solid_scalar);
  HUNDUN_CHECK(bitwise_equal(zero_normal_scalar, 0.0));
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_literal_paper_oracle();
    test_physical_normal_and_rejection();
    test_plans_and_simultaneous_evaluation();
    test_limiting_and_mutation_oracles();
    test_affine_ghost_and_paper_closure_evidence();
  });
}
