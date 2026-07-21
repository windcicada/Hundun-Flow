// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/flow/flow_state.hpp"

#include <cstdint>
#include <type_traits>

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "Task 18 tests-on targets must share the test-access class definition"
#endif

using PressureCorrect = hundun::flow::PressureCorrectionReport (
    hundun::flow::PisoCoupler::*)(
    hundun::flow::FlowState &, double,
    const hundun::runtime::FieldView<const double> &,
    const hundun::linear::SolveControl &) const;

static_assert(!std::is_copy_constructible_v<hundun::flow::FlowState>);
static_assert(!std::is_copy_constructible_v<hundun::flow::PisoCoupler>);
static_assert(
    !std::is_copy_constructible_v<hundun::flow::FixedStepConstantDensityFlow>);
static_assert(std::is_enum_v<hundun::flow::PressureCorrectionDisposition>);
static_assert(std::is_same_v<
              std::underlying_type_t<
                  hundun::flow::PressureCorrectionDisposition>,
              std::uint8_t>);
static_assert(std::is_final_v<hundun::flow::PressureCorrectionReport>);
static_assert(std::is_same_v<decltype(&hundun::flow::PisoCoupler::correct),
                             PressureCorrect>);
static_assert(std::is_same_v<
              decltype(hundun::flow::PressureCorrectionReport::disposition),
              hundun::flow::PressureCorrectionDisposition>);
static_assert(std::is_same_v<
              decltype(hundun::flow::PressureCorrectionReport::reason),
              hundun::flow::StepFailureReason>);
static_assert(std::is_same_v<
              decltype(hundun::flow::PressureCorrectionReport::
                           lowest_failing_rank),
              int>);

int main() { return 0; }
