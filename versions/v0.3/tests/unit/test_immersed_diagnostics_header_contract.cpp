// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_immersed_module.hpp"

#include <cstdint>
#include <type_traits>
#include <utility>

using Source = hundun::flow::ImmersedFlowDiagnosticSource;
using Flow = hundun::flow::FixedStepImmersedFlow;

static_assert(std::is_final_v<Source>);
static_assert(std::is_copy_constructible_v<Source>);
static_assert(std::is_same_v<
              decltype(std::declval<const Flow &>().diagnostic_source(
                  std::declval<const hundun::flow::FlowState &>(),
                  std::declval<const hundun::flow::ImmersedFlowStepAttemptReport &>())),
              Source>);
static_assert(std::is_same_v<
              decltype(std::declval<const Source &>().committed_step()),
              std::uint64_t>);
static_assert(std::is_same_v<
              decltype(std::declval<const Source &>()
                           .maximum_wall_penetration_m_per_s()),
              double>);

int main() { return 0; }
