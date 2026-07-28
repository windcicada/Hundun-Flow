// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace hundun::application::detail {

inline std::uint64_t performance_crc64_ecma(std::string_view input) noexcept {
  constexpr std::uint64_t polynomial = UINT64_C(0x42f0e1eba9ea3693);
  std::uint64_t crc = 0U;
  for (const char character : input) {
    const auto byte = static_cast<unsigned char>(character);
    crc ^= static_cast<std::uint64_t>(byte) << 56U;
    for (unsigned bit = 0; bit < 8U; ++bit)
      crc = (crc & (UINT64_C(1) << 63U)) != 0U
                ? (crc << 1U) ^ polynomial
                : crc << 1U;
  }
  return crc;
}

inline std::string tagged_performance_crc64(std::string_view tag,
                                            std::string_view input) {
  constexpr char digits[] = "0123456789abcdef";
  const std::uint64_t crc = performance_crc64_ecma(input);
  std::string result(tag);
  result.reserve(tag.size() + 17U);
  result.push_back(':');
  for (int shift = 60; shift >= 0; shift -= 4)
    result.push_back(digits[(crc >> static_cast<unsigned>(shift)) & 0x0fU]);
  return result;
}

}  // namespace hundun::application::detail
