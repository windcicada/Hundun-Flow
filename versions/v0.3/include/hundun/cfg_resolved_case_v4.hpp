// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case_v3.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace hundun::config {

enum class PressureConstraintMode : std::uint8_t {
  open_fixed_p0,
  closed,
  partially_closed
};

enum class ReactingThermalMode : std::uint8_t { adiabatic, isothermal };

struct ReactingMechanismConfig final {
  std::filesystem::path file;
  std::string sha256;
  std::string phase;
};

struct ChemistrySolverConfig final {
  double relative_tolerance{};
  double absolute_tolerance{};
  int maximum_internal_steps{};
};

struct ReactingThermalBoundaryConfig final {
  ReactingThermalMode mode{ReactingThermalMode::adiabatic};
  std::optional<double> temperature_k;
};

struct ReactingBoundaryConfig final {
  bool non_catalytic_impermeable{};
  std::optional<ReactingThermalBoundaryConfig> thermal;
};

struct ResolvedReactingCaseV4 final {
  int schema_version{4};
  FlowCaseConfig common_flow;
  ImmersedBoundaryConfig immersed_boundary;
  LesConfig les;
  ReactingMechanismConfig mechanism;
  ChemistrySolverConfig chemistry;
  double initial_p0_pa{};
  double initial_temperature_k{};
  std::vector<std::string> species_names;
  std::vector<double> initial_mass_fractions;
  PressureConstraintMode pressure_mode{PressureConstraintMode::open_fixed_p0};
  std::array<ReactingBoundaryConfig, 6> boundary_reacting;
  std::uint64_t composition_fingerprint{};
};

} // namespace hundun::config
