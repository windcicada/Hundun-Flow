// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/lin_preconditioners.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace hundun::linear {

struct IdentityPreconditioner::State final {
  explicit State(execution::ExecutionContext& context_value) noexcept
      : context(&context_value) {}

  execution::ExecutionContext* context;
  const LinearOperator* linear_operator{};
  std::uint64_t revision{};
  VectorLayout domain;
  VectorLayout range;
  bool cache_valid{false};
};

struct JacobiPreconditioner::State final {
  explicit State(execution::ExecutionContext& context_value) noexcept
      : context(&context_value) {}

  execution::ExecutionContext* context;
  const LinearOperator* linear_operator{};
  std::uint64_t revision{};
  VectorLayout domain;
  VectorLayout range;
  std::optional<execution::Buffer> inverse_diagonal;
  std::optional<execution::Buffer> staging_inverse_diagonal;
  bool cache_valid{false};
};

#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
namespace detail {

void jacobi_storage_raw(
    const std::optional<execution::Buffer>& inverse_diagonal,
    const std::optional<execution::Buffer>& staging_inverse_diagonal,
    std::size_t owned_count,
    std::array<execution::AllocationIdentity, 2>& allocation_identities,
    std::array<std::size_t, 2>& byte_sizes,
    std::vector<double>& cached_inverse);
void arm_fail_next_jacobi_cold_update_before_publication_raw() noexcept;
void reset_jacobi_cold_update_fault_raw() noexcept;

}  // namespace detail
#endif

}  // namespace hundun::linear
