// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_stage3_performance.hpp"

#include <array>
#include <cstdint>
#include <type_traits>

static_assert(std::is_final_v<hundun::diagnostics::Stage3PerformanceCounters>);
static_assert(hundun::diagnostics::kStage3PerformanceCounterIds.size() == 17U);
static_assert(std::is_same_v<
              decltype(hundun::diagnostics::Stage3PerformanceCounters::
                           checkpoint_logical_io_bytes),
              std::uint64_t>);

int main() { return 0; }
