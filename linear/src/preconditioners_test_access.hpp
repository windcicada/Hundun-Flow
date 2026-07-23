// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/execution/execution.hpp"
#include "hundun/linear/preconditioners.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::linear::test {

struct JacobiStorageSnapshot final {
  std::array<execution::AllocationIdentity, 2> allocation_identities{};
  std::array<std::size_t, 2> byte_sizes{};
  std::vector<double> cached_inverse;
  std::uint64_t revision{};
  bool cache_valid{};
};

class PreconditionerTestAccess final {
 public:
  static JacobiStorageSnapshot
  jacobi_storage(const JacobiPreconditioner& preconditioner);
  static void arm_fail_next_cold_update_before_publication() noexcept;
  static void reset_cold_update_fault() noexcept;
};

}  // namespace hundun::linear::test
