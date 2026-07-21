// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/flow/flow_state.hpp"

#include <type_traits>

static_assert(!std::is_copy_constructible_v<hundun::flow::FlowState>);
static_assert(!std::is_copy_constructible_v<hundun::flow::PisoCoupler>);
static_assert(
    !std::is_copy_constructible_v<hundun::flow::FixedStepConstantDensityFlow>);

int main() { return 0; }
