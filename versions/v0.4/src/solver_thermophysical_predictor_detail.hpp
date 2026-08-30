// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace hundun::v04::detail {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)

// Arm this on every participating rank.  The next second-order predictor
// call reports the same collective low-BDF-base failure as the physical
// admissibility check, then clears the one-shot state before a BE fallback.
void arm_low_bdf_source_base_failure_once_for_test() noexcept;
void clear_low_bdf_source_base_failure_for_test() noexcept;

#endif

}  // namespace hundun::v04::detail
