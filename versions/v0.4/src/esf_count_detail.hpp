// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <cstddef>
namespace hundun::v04::esf {
// Upper storage bound from the native Restart field budget. Actual admission
// also accounts for mean/scalar/cache fields and the per-rank memory budget.
inline constexpr std::size_t maximum_fields = 64;
inline constexpr bool valid_field_count(std::size_t n) noexcept {
  return n >= 2 && n <= maximum_fields && n % 2 == 0;
}
} // namespace hundun::v04::esf
