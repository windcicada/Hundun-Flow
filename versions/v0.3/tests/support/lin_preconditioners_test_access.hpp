// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#ifndef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
#error "Preconditioner test access is available only in tests-on builds"
#endif

#include "hundun/exec_execution.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "src/lin_preconditioners_detail.hpp"

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
  jacobi_storage(const JacobiPreconditioner& preconditioner) {
    JacobiStorageSnapshot snapshot;
    if (!preconditioner.state_) {
      return snapshot;
    }
    const auto& state = *preconditioner.state_;
    snapshot.revision = state.revision;
    snapshot.cache_valid = state.cache_valid;
    detail::jacobi_storage_raw(
        state.inverse_diagonal, state.staging_inverse_diagonal,
        state.domain.owned_count(), snapshot.allocation_identities,
        snapshot.byte_sizes, snapshot.cached_inverse);
    return snapshot;
  }

  static void arm_fail_next_cold_update_before_publication() noexcept {
    detail::arm_fail_next_jacobi_cold_update_before_publication_raw();
  }

  static void reset_cold_update_fault() noexcept {
    detail::reset_jacobi_cold_update_fault_raw();
  }
};

}  // namespace hundun::linear::test
