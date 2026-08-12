// SPDX-License-Identifier: Apache-2.0

#include "hundun/les_wale.hpp"

#include <type_traits>

using Model = hundun::les::WaleModel;
using Coefficients = hundun::les::WaleAttemptCoefficients;

static_assert(std::is_final_v<Model>);
static_assert(std::is_nothrow_move_constructible_v<Model>);
static_assert(!std::is_copy_constructible_v<Model>);
static_assert(std::is_final_v<Coefficients>);
static_assert(std::is_nothrow_move_constructible_v<Coefficients>);
static_assert(!std::is_copy_constructible_v<Coefficients>);
static_assert(std::is_same_v<
              decltype(std::declval<const Coefficients &>().identity()),
              hundun::les::WaleCoefficientIdentity>);
static_assert(std::is_same_v<
              decltype(std::declval<const Model &>().performance_counters()),
              hundun::diagnostics::Stage3PerformanceCounters>);

int main() { return 0; }
