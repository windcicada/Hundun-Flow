// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <variant>

namespace hundun::config {

enum class ImmersedBoundaryModel : std::uint8_t {
  none,
  local_flow_pattern_ghost_cell
};

enum class ImmersedFluidSide : std::uint8_t { inside, outside };

enum class LesModel : std::uint8_t { none, wale };

struct StlGeometryConfig final {
  std::filesystem::path file;
  double length_scale_to_m{};
  ImmersedFluidSide fluid_side{ImmersedFluidSide::outside};
};

struct StaticImmersedWallConfig final {
  runtime::Real3 velocity_m_per_s{};
};

struct ImmersedBoundaryConfig final {
  ImmersedBoundaryModel model{ImmersedBoundaryModel::none};
  std::optional<StlGeometryConfig> geometry;
  std::optional<StaticImmersedWallConfig> wall;
};

struct WaleConfig final {
  double coefficient{};
  double turbulent_prandtl{};
  double turbulent_schmidt{};
};

struct LesConfig final {
  LesModel model{LesModel::none};
  std::optional<WaleConfig> wale;
};

struct ImmersedFlowCaseConfig final {
  int schema_version{3};
  FlowCaseConfig common_flow;
  ImmersedBoundaryConfig immersed_boundary;
  LesConfig les;
};

using ResolvedCaseV3 =
    std::variant<CaseConfig, FlowCaseConfig, ImmersedFlowCaseConfig>;

} // namespace hundun::config
