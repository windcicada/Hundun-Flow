// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cmath>
#include <cstddef>

namespace hundun::test::stage3 {

struct Dual3 final {
  double value{};
  std::array<double, 3> gradient{};
  std::array<std::array<double, 3>, 3> hessian{};

  constexpr Dual3() noexcept = default;
  constexpr Dual3(double supplied) noexcept : value(supplied) {}

  static Dual3 variable(double supplied, std::size_t axis) {
    Dual3 result(supplied);
    result.gradient.at(axis) = 1.0;
    return result;
  }
};

inline Dual3 operator+(const Dual3 &left, const Dual3 &right) noexcept {
  Dual3 result(left.value + right.value);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] = left.gradient[i] + right.gradient[i];
    for (std::size_t j = 0U; j < 3U; ++j)
      result.hessian[i][j] = left.hessian[i][j] + right.hessian[i][j];
  }
  return result;
}

inline Dual3 operator-(const Dual3 &left, const Dual3 &right) noexcept {
  Dual3 result(left.value - right.value);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] = left.gradient[i] - right.gradient[i];
    for (std::size_t j = 0U; j < 3U; ++j)
      result.hessian[i][j] = left.hessian[i][j] - right.hessian[i][j];
  }
  return result;
}

inline Dual3 operator-(const Dual3 &value) noexcept {
  Dual3 result(-value.value);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] = -value.gradient[i];
    for (std::size_t j = 0U; j < 3U; ++j)
      result.hessian[i][j] = -value.hessian[i][j];
  }
  return result;
}

inline Dual3 operator*(const Dual3 &left, const Dual3 &right) noexcept {
  Dual3 result(left.value * right.value);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] =
        left.gradient[i] * right.value + left.value * right.gradient[i];
    for (std::size_t j = 0U; j < 3U; ++j) {
      result.hessian[i][j] = left.hessian[i][j] * right.value +
                             left.gradient[i] * right.gradient[j] +
                             left.gradient[j] * right.gradient[i] +
                             left.value * right.hessian[i][j];
    }
  }
  return result;
}

inline Dual3 reciprocal(const Dual3 &input) noexcept {
  const double inverse = 1.0 / input.value;
  const double inverse2 = inverse * inverse;
  const double inverse3 = inverse2 * inverse;
  Dual3 result(inverse);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] = -input.gradient[i] * inverse2;
    for (std::size_t j = 0U; j < 3U; ++j) {
      result.hessian[i][j] =
          2.0 * input.gradient[i] * input.gradient[j] * inverse3 -
          input.hessian[i][j] * inverse2;
    }
  }
  return result;
}

inline Dual3 operator/(const Dual3 &left, const Dual3 &right) noexcept {
  return left * reciprocal(right);
}

inline Dual3 sin(const Dual3 &input) noexcept {
  const double sine = std::sin(input.value);
  const double cosine = std::cos(input.value);
  Dual3 result(sine);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] = cosine * input.gradient[i];
    for (std::size_t j = 0U; j < 3U; ++j) {
      result.hessian[i][j] = cosine * input.hessian[i][j] -
                             sine * input.gradient[i] * input.gradient[j];
    }
  }
  return result;
}

inline Dual3 cos(const Dual3 &input) noexcept {
  const double sine = std::sin(input.value);
  const double cosine = std::cos(input.value);
  Dual3 result(cosine);
  for (std::size_t i = 0U; i < 3U; ++i) {
    result.gradient[i] = -sine * input.gradient[i];
    for (std::size_t j = 0U; j < 3U; ++j) {
      result.hessian[i][j] = -sine * input.hessian[i][j] -
                             cosine * input.gradient[i] * input.gradient[j];
    }
  }
  return result;
}

} // namespace hundun::test::stage3
