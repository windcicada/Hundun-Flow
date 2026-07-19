// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/resolved_case.hpp"

#include <filesystem>
#include <utility>

namespace hundun::application::detail {

template <class RootPhase, class ConfigPhase, class RankAndDispatchPhase>
int dispatch_in_root_config_rank_order(
    RootPhase&& root_phase, ConfigPhase&& config_phase,
    RankAndDispatchPhase&& rank_and_dispatch_phase) {
  const std::filesystem::path authoritative_root =
      std::forward<RootPhase>(root_phase)();
  const config::ResolvedCase resolved =
      std::forward<ConfigPhase>(config_phase)();
  return std::forward<RankAndDispatchPhase>(rank_and_dispatch_phase)(
      resolved, authoritative_root);
}

}  // namespace hundun::application::detail
