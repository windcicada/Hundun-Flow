// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <optional>

namespace hundun::finite_volume::detail {

enum class PoissonApplyPhase {
  halo_begin,
  local_rows,
  halo_wait,
  partition_rows
};

struct PoissonTestOptions final {
  std::optional<PoissonApplyPhase> injected_apply_failure{};
  bool force_revision_wrap{};
  bool observe_apply{};
  bool inject_reentrant_apply{};
};

struct PoissonApplyTrace final {
  std::array<PoissonApplyPhase, 4> phases{};
  std::size_t count{};
};

void set_poisson_test_options(PoissonTestOptions options) noexcept;
void reset_poisson_apply_trace() noexcept;
PoissonApplyTrace poisson_apply_trace() noexcept;
void record_poisson_apply_phase(PoissonApplyPhase phase);
void record_poisson_cleanup_wait() noexcept;
bool consume_force_revision_wrap() noexcept;
bool consume_injected_reentrant_apply() noexcept;

}  // namespace hundun::finite_volume::detail
