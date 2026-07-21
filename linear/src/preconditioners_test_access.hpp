// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/execution/execution.hpp"
#include "hundun/linear/preconditioners.hpp"

#include <array>
#include <cstdint>

namespace hundun::linear::test {

struct JacobiStorageSnapshot final {
  std::array<execution::AllocationIdentity, 2> allocation_identities{};
  std::uint64_t revision{};
  bool cache_valid{};
};

class PreconditionerTestAccess final {
 public:
  static JacobiStorageSnapshot
  jacobi_storage(const JacobiPreconditioner& preconditioner) noexcept;
};

}  // namespace hundun::linear::test
