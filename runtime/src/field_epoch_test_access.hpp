// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>

namespace hundun::runtime {

class FieldStorage;

namespace detail {

struct FieldEpochTestAccess final {
  static std::uint64_t generation(const FieldStorage &storage) noexcept;
  static void force_generation(FieldStorage &storage,
                               std::uint64_t generation) noexcept;
};

}  // namespace detail
}  // namespace hundun::runtime
