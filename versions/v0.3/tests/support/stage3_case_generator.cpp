// SPDX-License-Identifier: Apache-2.0

#include "tests/support/stage3_performance_evidence.hpp"

#include <array>

namespace hundun::test {

runtime::Int3 stage3_performance_process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw std::invalid_argument("unsupported Stage 3 performance rank count");
}

config::FlowCaseConfig stage3_performance_case(int ranks, int cells) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-performance";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = stage3_performance_process_grid(ranks);
  result.mesh.cells = {cells, cells, cells};
  result.mesh.origin_m = {-0.5, -0.5, -0.5};
  result.mesh.length_m = {2.0, 2.0, 2.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = 0.01;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  result.time.mode = config::TimeMode::fixed;
  result.time.steps = 1;
  result.time.initial_dt_s = 0.01;
  result.time.min_dt_s = 0.01;
  result.time.max_dt_s = 0.01;
  result.time.cfl_target = 0.5;
  result.time.diffusion_number_target = 0.25;
  result.time.growth_factor = 1.25;
  result.time.retry_factor = 0.5;
  result.time.max_retries = 8;
  result.restart.read = false;
  result.restart.write_directory = "checkpoint";
  result.restart.write_interval = 1;
  result.diagnostics.directory = "diagnostics";
  result.diagnostics.write_interval = 1;
  result.performance.enabled = false;
  result.performance.directory = "performance";
  result.performance.warmup_steps = 1;
  result.performance.measured_steps = 1;
  result.performance.repetitions = 1;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t index = 0U; index < names.size(); ++index) {
    result.boundaries[index].patch = names[index];
    result.boundaries[index].type = config::BoundaryType::periodic;
  }
  return result;
}

std::vector<StlFixtureTriangle> stage3_performance_body() {
  auto coarse = outward_cube();
  for (auto &triangle : coarse)
    for (auto &vertex : triangle.vertices) {
      vertex.x = 0.30 + 0.40 * vertex.x;
      vertex.y = 0.30 + 0.40 * vertex.y;
      vertex.z = 0.30 + 0.40 * vertex.z;
    }
  const auto midpoint = [](runtime::Real3 left, runtime::Real3 right) {
    return runtime::Real3{0.5 * (left.x + right.x),
                          0.5 * (left.y + right.y),
                          0.5 * (left.z + right.z)};
  };
  std::vector<StlFixtureTriangle> result;
  result.reserve(4U * coarse.size());
  for (const auto &triangle : coarse) {
    const auto ab = midpoint(triangle.vertices[0], triangle.vertices[1]);
    const auto bc = midpoint(triangle.vertices[1], triangle.vertices[2]);
    const auto ca = midpoint(triangle.vertices[2], triangle.vertices[0]);
    result.push_back({triangle.file_normal, {triangle.vertices[0], ab, ca}});
    result.push_back({triangle.file_normal, {ab, triangle.vertices[1], bc}});
    result.push_back({triangle.file_normal, {ca, bc, triangle.vertices[2]}});
    result.push_back({triangle.file_normal, {ab, bc, ca}});
  }
  return result;
}

} // namespace hundun::test
