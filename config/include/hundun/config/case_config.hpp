// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/types.hpp"

#include <array>
#include <filesystem>
#include <optional>
#include <string>

namespace hundun::config {

struct MeshConfig {
  runtime::Int3 cells;
  runtime::Real3 origin_m;
  runtime::Real3 length_m;
  std::array<bool, 3> periodic;
};

struct TimeConfig {
  double dt_s;
  int steps;
};

struct TransportConfig {
  runtime::Real3 velocity_m_per_s;
  double diffusivity_m2_per_s;
};

struct OutputConfig {
  std::filesystem::path directory;
  int write_interval;
  int restart_interval;
};

struct CaseConfig {
  int schema_version;
  std::string case_name;
  std::optional<int> expected_ranks;
  std::optional<runtime::Int3> process_grid;
  MeshConfig mesh;
  TimeConfig time;
  TransportConfig transport;
  std::string initial_condition;
  OutputConfig output;
};

}  // namespace hundun::config
