// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include <string_view>

namespace hundun::diagnostics::detail {

std::uint64_t structured_crc64_ecma_for_test(std::string_view) noexcept;

} // namespace hundun::diagnostics::detail

namespace hundun::diagnostics::test {

inline std::uint64_t crc64_ecma(std::string_view bytes) noexcept {
  return detail::structured_crc64_ecma_for_test(bytes);
}

} // namespace hundun::diagnostics::test
