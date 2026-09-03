// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_mesh.hpp"

#include <mpi.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <type_traits>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

CartesianMeshSpec uniform_mesh(Int3 cells = {8, 6, 4}) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-1.0, 2.0, 0.5};
  mesh.upper = {3.0, 5.0, 2.5};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {1.0e-6, 1.0e-6, 1.0e-6};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

bool same(Int3 left, Int3 right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool test_uniform_metrics() {
  static_assert(!std::is_copy_constructible_v<AxisMetrics>);
  static_assert(!std::is_copy_assignable_v<AxisMetrics>);
  static_assert(std::is_nothrow_move_constructible_v<AxisMetrics>);
  static_assert(std::is_nothrow_move_assignable_v<AxisMetrics>);
  static_assert(!std::is_copy_constructible_v<CartesianGeometryPlan>);
  static_assert(!std::is_copy_assignable_v<CartesianGeometryPlan>);
  static_assert(std::is_nothrow_move_constructible_v<CartesianGeometryPlan>);
  static_assert(std::is_nothrow_move_assignable_v<CartesianGeometryPlan>);
  const CartesianMeshSpec mesh = uniform_mesh();
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  const Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh, GeometryBudget{1024U, 128U}, geometry, patch);
  bool passed = expect(static_cast<bool>(status), "uniform geometry compiles");
  passed &= expect(same(geometry.global_cells(), mesh.exact_cells),
                   "exact global cells survive planning");
  passed &= expect(geometry.kind() == GeometryKind::uniform &&
                       geometry.topology_revision() == 1U &&
                       geometry.fingerprint() != 0U,
                   "uniform plan publishes immutable identity");
  passed &= expect(geometry.x().faces().size == 9U &&
                       geometry.x().faces().data[0] == -1.0 &&
                       geometry.x().faces().data[8] == 3.0,
                   "uniform x faces preserve exact endpoints");
  passed &= expect(geometry.x().uniform() &&
                       geometry.x().uniform_width() == 0.5 &&
                       geometry.x().uniform_inverse_width() == 2.0,
                   "uniform fast metrics are cached");
  passed &= expect(geometry.x().centres().size == 8U &&
                       geometry.x().widths().size == 8U &&
                       geometry.x().inverse_widths().size == 8U &&
                       geometry.x().inverse_centre_distances().size == 7U,
                   "all one-dimensional metrics are cached once");
  passed &= expect(same(patch.begin, {0, 0, 0}) &&
                       same(patch.cells, mesh.exact_cells) &&
                       same(patch.process_grid, {1, 1, 1}),
                   "one rank owns one complete contiguous patch");

  CartesianGeometryPlan repeated;
  MeshPatch repeated_patch;
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, mesh, GeometryBudget{1024U, 128U},
                       repeated, repeated_patch)),
                   "same uniform plan compiles repeatedly");
  passed &= expect(repeated.fingerprint() == geometry.fingerprint(),
                   "uniform coordinates are deterministic");
  return passed;
}

bool test_coast_runtime_axes_metrics() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::coast_runtime_axes_v1;
  mesh.axes_file = "axes.dat";
  mesh.coast_runtime_faces = {
      std::vector<double>{-0.1, -0.025, 0.05},
      std::vector<double>{-0.02, 0.0, 0.03},
      std::vector<double>{0.0, 0.06}};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {2, 2, 1};
  mesh.lower = {-0.1, -0.02, 0.0};
  mesh.upper = {0.05, 0.03, 0.06};
  mesh.limits = {32U, 1U << 20U};

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  const Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch);
  bool passed = expect(static_cast<bool>(status),
                       "COAST runtime axes compile directly");
  passed &= expect(same(geometry.global_cells(), {2, 2, 1}) &&
                       geometry.kind() ==
                           GeometryKind::coast_runtime_axes_v1,
                   "COAST runtime axes preserve topology");
  passed &= expect(
      status && geometry.x().faces().size == 3U &&
          geometry.x().faces().data[1] ==
              static_cast<double>(static_cast<float>(-0.025)),
      "COAST runtime axes use the solver's float32 effective coordinate");
  const PlanFingerprint retained = geometry.fingerprint();

  CartesianMeshSpec collapsed = mesh;
  collapsed.coast_runtime_faces[0U][1U] = -0.1 + 1.0e-12;
  const Status collapsed_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, collapsed, GeometryBudget{}, geometry, patch);
  passed &= expect(collapsed_status.code == StatusCode::invalid_plan &&
                       geometry.fingerprint() == retained,
                   "float32-collapsed COAST faces are rejected atomically");

  CartesianMeshSpec mismatched = mesh;
  mismatched.exact_cells.x = 3;
  const Status mismatched_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mismatched, GeometryBudget{}, geometry, patch);
  passed &= expect(mismatched_status.code == StatusCode::invalid_plan &&
                       geometry.fingerprint() == retained,
                   "COAST face counts must match declared cells");

  CartesianMeshSpec changed = mesh;
  changed.coast_runtime_faces[0U][1U] = -0.02;
  CartesianGeometryPlan changed_geometry;
  MeshPatch changed_patch;
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, changed, GeometryBudget{},
                       changed_geometry, changed_patch)) &&
                       changed_geometry.fingerprint() != retained,
                   "effective COAST faces bind geometry identity");
  return passed;
}

bool test_hard_limits_and_publication() {
  bool passed = true;
  CartesianGeometryPlan retained;
  MeshPatch retained_patch;
  CartesianMeshSpec valid = uniform_mesh();
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, valid, GeometryBudget{0U, 1U}, retained,
                       retained_patch)),
                   "baseline plan exists before rejected replacements");
  const auto retained_fingerprint = retained.fingerprint();

  CartesianMeshSpec cell_limited = valid;
  cell_limited.limits.max_global_cells = 100U;
  const Status cell_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, cell_limited, GeometryBudget{0U, 1U}, retained,
      retained_patch);
  passed &= expect(cell_status.code == StatusCode::invalid_plan,
                   "global cell hard limit rejects the plan");
  passed &= expect(retained.fingerprint() == retained_fingerprint,
                   "rejected cell limit preserves published geometry");

  CartesianMeshSpec memory_limited = valid;
  memory_limited.limits.max_memory_bytes_per_rank = 100U;
  const Status memory_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, memory_limited, GeometryBudget{64U, 128U}, retained,
      retained_patch);
  passed &= expect(memory_status.code == StatusCode::invalid_plan,
                   "conservative local memory limit rejects the plan");

  CartesianMeshSpec preallocation_limited = valid;
  preallocation_limited.exact_cells = {1000000, 1, 1};
  preallocation_limited.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  preallocation_limited.limits.max_global_cells = 1000000U;
  preallocation_limited.limits.max_memory_bytes_per_rank = 4096U;
  const Status preallocation_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, preallocation_limited, GeometryBudget{}, retained,
      retained_patch);
  passed &= expect(preallocation_status.code == StatusCode::invalid_plan,
                   "exact geometry peak memory rejects before face allocation");

  CartesianMeshSpec invalid = valid;
  invalid.has_exact_cells = false;
  const Status invalid_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, invalid, GeometryBudget{}, retained, retained_patch);
  passed &= expect(invalid_status.code == StatusCode::invalid_plan,
                   "uniform input without exact cells is rejected");

  const Status communicator_status = CartesianGeometryCompiler::compile(
      MPI_COMM_NULL, valid, GeometryBudget{}, retained, retained_patch);
  passed &= expect(communicator_status.code == StatusCode::invalid_plan,
                   "null communicator is rejected");
  return passed;
}

bool test_cpu_tiles() {
  const MeshPatch patch{{5, 7, 11}, {10, 7, 5}, {1, 1, 1}, {0, 0, 0}};
  std::vector<CpuTile> tiles;
  bool passed = expect(static_cast<bool>(build_cpu_tiles(
                           patch, Int3{4, 3, 2}, tiles)),
                       "rank-local CPU tiles compile");
  passed &= expect(tiles.size() == 27U, "edge-clipped tile count is exact");
  passed &= expect(same(tiles.front().begin, patch.begin) &&
                       same(tiles.front().cells, {4, 3, 2}),
                   "first tile preserves global patch coordinates");
  passed &= expect(same(tiles.back().begin, {13, 13, 15}) &&
                       same(tiles.back().cells, {2, 1, 1}),
                   "last tile is clipped without changing ownership");
  const std::vector<CpuTile> retained = tiles;
  const Status invalid = build_cpu_tiles(patch, Int3{0, 3, 2}, tiles);
  passed &= expect(invalid.code == StatusCode::invalid_plan &&
                       tiles.size() == retained.size(),
                   "invalid tile shape preserves published schedule");
  return passed;
}

bool test_extreme_decomposition_arithmetic() {
  const MeshPatch patch{{0, 0, 0},
                        {std::numeric_limits<std::int32_t>::max(), 1, 1},
                        {std::numeric_limits<std::int32_t>::max(), 1, 1},
                        {0, 0, 0}};
  std::vector<CpuTile> tiles;
  bool passed = expect(static_cast<bool>(build_cpu_tiles(
                           patch,
                           Int3{std::numeric_limits<std::int32_t>::max(), 1, 1},
                           tiles)),
                       "extreme int32 extent avoids ceil-division overflow");
  passed &= expect(tiles.size() == 1U &&
                       tiles.front().cells.x ==
                           std::numeric_limits<std::int32_t>::max(),
                   "extreme extent produces one exact tile");
  const MeshPatch overflowing_begin{
      {std::numeric_limits<std::int32_t>::max(), 0, 0},
      {2, 1, 1}, {1, 1, 1}, {0, 0, 0}};
  const std::size_t retained_size = tiles.size();
  const Status overflow_status =
      build_cpu_tiles(overflowing_begin, Int3{1, 1, 1}, tiles);
  passed &= expect(overflow_status.code == StatusCode::invalid_plan &&
                       tiles.size() == retained_size,
                   "tile global begin overflow is rejected atomically");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_uniform_metrics();
  passed &= test_coast_runtime_axes_metrics();
  passed &= test_hard_limits_and_publication();
  passed &= test_cpu_tiles();
  passed &= test_extreme_decomposition_arithmetic();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 Cartesian geometry tests passed\n";
  return 0;
}
