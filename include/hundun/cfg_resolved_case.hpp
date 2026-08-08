// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_case_config.hpp"

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace hundun::config {

enum class SimulationType { passive_scalar, variable_density_flow };
enum class DensityModel { constant, material, ideal_gas };
enum class MeshMapping { uniform_box, analytic_warped_box };
enum class TimeMode { fixed, adaptive };
enum class PatchName { x_min, x_max, y_min, y_max, z_min, z_max };
enum class BoundaryType {
  periodic,
  no_slip_wall,
  symmetry,
  velocity_inlet,
  pressure_outlet
};
enum class InletThermalAuthority { temperature, enthalpy };

struct FlowResourcesConfig {
  std::optional<int> expected_ranks;
  std::optional<runtime::Int3> process_grid;
};

struct FlowMeshConfig {
  runtime::Int3 cells;
  runtime::Real3 origin_m;
  runtime::Real3 length_m;
  MeshMapping mapping;
  std::optional<runtime::Real3> warp_amplitude;
};

struct FlowTimeConfig {
  TimeMode mode;
  int steps;
  double initial_dt_s;
  double min_dt_s;
  double max_dt_s;
  double cfl_target;
  double diffusion_number_target;
  double growth_factor;
  double retry_factor;
  int max_retries;
};

struct FlowPhysicsConfig {
  double rho_ref_kg_per_m3;
  double dynamic_viscosity_pa_s;
  double inlet_consistency_rtol;
  std::optional<double> cp_J_per_kg_K;
  std::optional<double> gas_constant_J_per_kg_K;
  std::optional<double> thermodynamic_pressure_pa;
};

struct FlowScalarConfig {
  std::string name;
  double diffusivity_m2_per_s;
};

struct InletScalarValue {
  std::string name;
  double value;
};

struct FlowBoundaryConfig {
  PatchName patch;
  BoundaryType type;
  std::optional<runtime::Real3> velocity_m_per_s;
  std::optional<InletThermalAuthority> thermal_authority;
  std::optional<double> temperature_K;
  std::optional<double> enthalpy_J_per_kg;
  std::optional<double> density_kg_per_m3;
  std::optional<std::vector<InletScalarValue>> scalar_values;
  std::optional<double> pressure_perturbation_pa;
};

struct FlowRestartConfig {
  bool read;
  std::optional<std::filesystem::path> read_directory;
  std::filesystem::path write_directory;
  int write_interval;
};

struct FlowDiagnosticsConfig {
  std::filesystem::path directory;
  int write_interval;
  bool write_mesh;
};

struct FlowPerformanceConfig {
  bool enabled;
  std::filesystem::path directory;
  int warmup_steps;
  int measured_steps;
  int repetitions;
};

struct FlowCaseConfig {
  int schema_version;
  std::string case_name;
  SimulationType simulation_type;
  DensityModel density_model;
  FlowResourcesConfig resources;
  FlowMeshConfig mesh;
  FlowTimeConfig time;
  FlowPhysicsConfig physics;
  std::vector<FlowScalarConfig> scalars;
  std::array<FlowBoundaryConfig, 6> boundaries;
  FlowRestartConfig restart;
  FlowDiagnosticsConfig diagnostics;
  FlowPerformanceConfig performance;
};

using ResolvedCase = std::variant<CaseConfig, FlowCaseConfig>;

}  // namespace hundun::config
