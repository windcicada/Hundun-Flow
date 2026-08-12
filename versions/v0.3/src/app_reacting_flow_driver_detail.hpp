// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case_v4.hpp"

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::application {

struct ReactingDriverPlan final {
  bool closed_pressure{};
  bool immersed_boundary{};
  bool wale{};
  std::uint32_t backend_runtime_count{};
  std::uint32_t workspace_pool_count{};
  std::vector<std::string> operator_order;
};

ReactingDriverPlan
plan_reacting_flow_case(const config::ResolvedReactingCaseV4 &);

int run_reacting_flow_case(const config::ResolvedReactingCaseV4 &, MPI_Comm);

} // namespace hundun::application
