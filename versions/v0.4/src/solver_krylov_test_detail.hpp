// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_linear.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace hundun::v04::detail {

inline constexpr std::size_t kMaximumFgmresBasisUpdateCountForTest =
    (std::numeric_limits<std::uint8_t>::max() - 8U) / 2U;
static_assert(kMaximumFgmresBasisUpdateCountForTest == 123U);

// This declaration is private to the test-core build.  It is intentionally
// absent from every installed/public header and is not compiled into the
// production core.
Status fused_krylov_basis_update_for_test(
    FieldView destination, const ConstFieldView* sources,
    const double* scales, std::size_t count) noexcept;

bool krylov_basis_update_inputs_disjoint_for_test(
    FieldView destination, const ConstFieldView* sources,
    std::size_t count) noexcept;

Status single_reduction_fgmres_norm_for_test(
    Span<const double> projections, double bar_norm, double& t,
    double& next_norm) noexcept;

void force_single_reduction_fgmres_breakdown_for_test(
    std::size_t projection_count) noexcept;

}  // namespace hundun::v04::detail
#endif
