// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace hundun::v04::detail {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)

struct HypreLifecycleSnapshot {
  std::uint64_t matrix_creates{};
  std::uint64_t matrix_destroys{};
  std::uint64_t solver_creates{};
  std::uint64_t solver_destroys{};
  std::uint64_t solver_setups{};
};

enum class HypreFailurePoint : std::uint8_t {
  none,
  native_objects,
  matrix_fill,
  solver_setup,
};

HypreLifecycleSnapshot hypre_lifecycle_snapshot_for_test() noexcept;
void reset_hypre_lifecycle_for_test() noexcept;
void set_hypre_failure_for_test(HypreFailurePoint point) noexcept;
void clear_hypre_failure_for_test() noexcept;

#endif

}  // namespace hundun::v04::detail
