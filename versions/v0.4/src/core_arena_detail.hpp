// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

inline constexpr std::size_t kCacheLineBytes = 64U;
inline constexpr std::size_t kDoublesPerCacheLine =
    kCacheLineBytes / sizeof(double);

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& out) noexcept;
bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& out) noexcept;
bool checked_align(std::size_t value, std::size_t alignment,
                   std::size_t& out) noexcept;
std::uint64_t issue_identity() noexcept;

}  // namespace hundun::v04::detail
