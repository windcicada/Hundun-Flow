// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/les_wale_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>

namespace {

double oracle(const std::array<double, 9> &g, double coefficient,
              double delta) {
  const auto at = [](const auto &a, std::size_t i, std::size_t j) {
    return a[3U * i + j];
  };
  std::array<double, 9> s{};
  std::array<double, 9> gg{};
  for (std::size_t i = 0U; i < 3U; ++i)
    for (std::size_t j = 0U; j < 3U; ++j) {
      s[3U * i + j] = 0.5 * (at(g, i, j) + at(g, j, i));
      for (std::size_t k = 0U; k < 3U; ++k)
        gg[3U * i + j] += at(g, i, k) * at(g, k, j);
    }
  const double trace = gg[0] + gg[4] + gg[8];
  double ss{};
  double sd2{};
  for (std::size_t i = 0U; i < 3U; ++i)
    for (std::size_t j = 0U; j < 3U; ++j) {
      const double sd = 0.5 * (at(gg, i, j) + at(gg, j, i)) -
                        (i == j ? trace / 3.0 : 0.0);
      ss += s[3U * i + j] * s[3U * i + j];
      sd2 += sd * sd;
    }
  if (ss == 0.0 && sd2 == 0.0)
    return 0.0;
  return std::pow(coefficient * delta, 2) * std::pow(sd2, 1.5) /
         (std::pow(ss, 2.5) + std::pow(sd2, 1.25));
}

bool near(double left, double right, double tolerance = 2.0e-13) {
  return std::abs(left - right) <=
         tolerance * std::max({1.0, std::abs(left), std::abs(right)});
}

std::array<double, 9> rotate_z_90(const std::array<double, 9> &g) {
  const std::array<double, 9> q{0.0, -1.0, 0.0, 1.0, 0.0,
                                0.0, 0.0,  0.0, 1.0};
  const auto at = [](const auto &a, std::size_t i, std::size_t j) {
    return a[3U * i + j];
  };
  std::array<double, 9> result{};
  for (std::size_t i = 0U; i < 3U; ++i)
    for (std::size_t j = 0U; j < 3U; ++j)
      for (std::size_t k = 0U; k < 3U; ++k)
        for (std::size_t l = 0U; l < 3U; ++l)
          result[3U * i + j] +=
              at(q, i, k) * at(g, k, l) * at(q, j, l);
  return result;
}

} // namespace

int main() {
  return hundun::test::run([] {
    constexpr double cw = 0.5;
    constexpr double delta = 0.125;
    const std::array<double, 9> tensor{0.2, -1.1, 0.7, 0.4, 0.3,
                                       0.9, -0.2, 0.8, -0.5};
    const double expected = oracle(tensor, cw, delta);
    const double actual = hundun::les::test::wale_kinematic_viscosity_for_test(
        tensor, cw, delta);
    HUNDUN_CHECK(near(actual, expected));
    HUNDUN_CHECK(near(
        hundun::les::test::wale_kinematic_viscosity_for_test(
            rotate_z_90(tensor), cw, delta),
        actual));

    std::array<double, 9> scaled = tensor;
    for (double &value : scaled)
      value *= 7.0;
    HUNDUN_CHECK(near(
        hundun::les::test::wale_kinematic_viscosity_for_test(scaled, cw,
                                                             delta),
        7.0 * actual));
    HUNDUN_CHECK(near(
        hundun::les::test::wale_kinematic_viscosity_for_test(tensor, cw,
                                                             3.0 * delta),
        9.0 * actual));

    const std::array<double, 9> zero{};
    const double exact_zero =
        hundun::les::test::wale_kinematic_viscosity_for_test(zero, cw, delta);
    HUNDUN_CHECK(exact_zero == 0.0);
    HUNDUN_CHECK(!std::signbit(exact_zero));

    const auto near_wall = [&](double y) {
      return hundun::les::test::wale_kinematic_viscosity_for_test(
          {0.0, 1.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, 0.0}, cw, delta);
    };
    const double y0 = 1.0e-4;
    const double slope = std::log(near_wall(2.0 * y0) / near_wall(y0)) /
                         std::log(2.0);
    HUNDUN_CHECK(slope >= 2.9 && slope <= 3.1);

    std::array<double, 9> non_finite = tensor;
    non_finite[4] = std::numeric_limits<double>::infinity();
    bool rejected = false;
    try {
      static_cast<void>(
          hundun::les::test::wale_kinematic_viscosity_for_test(non_finite, cw,
                                                               delta));
    } catch (...) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  });
}
