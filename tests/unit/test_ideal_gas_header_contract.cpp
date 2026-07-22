// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
#include "hundun/flow/ideal_gas_closure.hpp"
#include "hundun/flow/ideal_gas_piso.hpp"

#include <cstdint>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

using Closure = hundun::flow::IdealGasClosure;
using ClosureReport = hundun::flow::IdealGasClosureReport;
using ClosureSource = hundun::flow::IdealGasClosureDiagnosticSource;
using Flow = hundun::flow::FixedStepIdealGasFlow;
using StepReport = hundun::flow::IdealGasStepAttemptReport;

static_assert(std::is_final_v<hundun::flow::IdealGasClosureSpec>);
static_assert(std::is_final_v<hundun::flow::IdealGasClosureState>);
static_assert(std::is_nothrow_destructible_v<Closure>);
static_assert(std::is_nothrow_move_constructible_v<Closure>);
static_assert(!std::is_copy_constructible_v<Closure>);
static_assert(!std::is_copy_assignable_v<Closure>);
static_assert(!std::is_move_assignable_v<Closure>);
static_assert(std::is_copy_constructible_v<ClosureReport>);
static_assert(std::is_nothrow_move_constructible_v<ClosureReport>);
static_assert(std::is_copy_constructible_v<StepReport>);
static_assert(std::is_nothrow_move_constructible_v<StepReport>);
static_assert(std::is_nothrow_destructible_v<ClosureSource>);
static_assert(std::is_nothrow_move_constructible_v<ClosureSource>);
static_assert(!std::is_copy_constructible_v<ClosureSource>);
static_assert(!std::is_copy_assignable_v<ClosureSource>);
static_assert(!std::is_move_assignable_v<ClosureSource>);
static_assert(std::is_nothrow_destructible_v<Flow>);
static_assert(std::is_nothrow_move_constructible_v<Flow>);
static_assert(!std::is_copy_constructible_v<Flow>);
static_assert(!std::is_move_assignable_v<Flow>);

static_assert(static_cast<std::uint8_t>(
                  hundun::flow::IdealGasPressureMode::closed_dynamic) == 0U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::IdealGasPressureMode::open_fixed) == 1U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::StepFailureReason::collective_operation) ==
              12U);
static_assert(static_cast<std::uint8_t>(
                  hundun::flow::StepFailureReason::density_closure_failure) ==
              13U);
static_assert(std::is_same_v<decltype(std::declval<const Closure &>().state()),
                             hundun::flow::IdealGasClosureState>);
static_assert(std::is_same_v<decltype(std::declval<const ClosureReport &>()
                                          .candidate_pressure_available()),
                             bool>);
static_assert(std::is_same_v<decltype(std::declval<const ClosureSource &>()
                                          .fingerprint_field_id(0U)),
                             std::string_view>);
static_assert(
    std::is_same_v<decltype(std::declval<const StepReport &>().flow()),
                   const hundun::flow::MaterialDensityStepAttemptReport &>);

} // namespace

int main() { return 0; }
