// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace hundun::test {
inline void check(bool condition, const char* expression,
                  const char* file, int line) {
  if (!condition) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             " check failed: " + expression);
  }
}
inline void check_near(double actual, double expected, double tolerance,
                       const char* file, int line) {
  if (std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             " numerical check failed");
  }
}
template <class Function>
int run(Function&& function) {
  try {
    function();
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
}
#define HUNDUN_CHECK(expr) \
  ::hundun::test::check(static_cast<bool>(expr), #expr, __FILE__, __LINE__)
#define HUNDUN_CHECK_NEAR(actual, expected, tolerance) \
  ::hundun::test::check_near((actual), (expected), (tolerance), __FILE__, __LINE__)
