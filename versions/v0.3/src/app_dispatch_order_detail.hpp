// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case.hpp"

#include <filesystem>
#include <utility>

namespace hundun::application::detail {

template <class RootPhase, class ConfigPhase, class RankAndDispatchPhase>
int dispatch_in_root_config_rank_order(
    RootPhase&& root_phase, ConfigPhase&& config_phase,
    RankAndDispatchPhase&& rank_and_dispatch_phase) {
  const std::filesystem::path authoritative_root =
      std::forward<RootPhase>(root_phase)();
  const auto resolved = std::forward<ConfigPhase>(config_phase)();
  return std::forward<RankAndDispatchPhase>(rank_and_dispatch_phase)(
      resolved, authoritative_root);
}

}  // namespace hundun::application::detail
