// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_mesh.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

CartesianMeshSpec focus_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::tensor_stretched;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {4.0, 2.0, 1.0};
  mesh.has_base_spacing = true;
  mesh.base_spacing = {0.25, 0.25, 0.25};
  mesh.minimum_spacing = {0.025, 0.025, 0.025};
  mesh.max_growth_ratio = 1.25;
  mesh.focus_regions.push_back(
      {{1.0, 0.5, 0.25}, {3.0, 1.5, 0.75}, {0.1, 0.1, 0.1}});
  mesh.limits = {100000000U, 1U << 30U};
  return mesh;
}

bool check_axis(const AxisMetrics& axis, double lower, double upper,
                double minimum, double maximum, double growth,
                double focus_lower, double focus_upper,
                double focus_target) {
  bool passed = true;
  const auto faces = axis.faces();
  const auto centres = axis.centres();
  const auto widths = axis.widths();
  passed &= expect(faces.size == widths.size + 1U &&
                       centres.size == widths.size,
                   "stretched axis metric extents agree");
  passed &= expect(faces.data[0] == lower &&
                       faces.data[faces.size - 1U] == upper,
                   "stretched axis preserves exact endpoints");
  for (std::size_t index = 0U; index < widths.size; ++index) {
    passed &= expect(faces.data[index] < faces.data[index + 1U],
                     "stretched faces are strictly monotone");
    passed &= expect(widths.data[index] >= minimum * (1.0 - 1.0e-12) &&
                         widths.data[index] <= maximum * (1.0 + 1.0e-12),
                     "stretched width satisfies minimum and base caps");
    if (centres.data[index] >= focus_lower &&
        centres.data[index] <= focus_upper) {
      passed &= expect(widths.data[index] <= focus_target * (1.0 + 1.0e-11),
                       "focus target is attained at focused centres");
    }
    if (index > 0U) {
      const double ratio = std::max(widths.data[index] / widths.data[index - 1U],
                                    widths.data[index - 1U] / widths.data[index]);
      passed &= expect(ratio <= growth * (1.0 + 1.0e-11),
                       "adjacent stretched cells obey growth hard limit");
    }
  }
  return passed;
}

bool test_automatic_focus_grid() {
  CartesianMeshSpec mesh = focus_mesh();
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  const Status compile_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch);
  bool passed = expect(static_cast<bool>(compile_status),
                       "automatic tensor-focus geometry compiles");
  if (!passed) {
    return false;
  }
  passed &= check_axis(geometry.x(), 0.0, 4.0, 0.025, 0.25, 1.25,
                       1.0, 3.0, 0.1);
  passed &= check_axis(geometry.y(), 0.0, 2.0, 0.025, 0.25, 1.25,
                       0.5, 1.5, 0.1);
  passed &= check_axis(geometry.z(), 0.0, 1.0, 0.025, 0.25, 1.25,
                       0.25, 0.75, 0.1);
  passed &= expect(!geometry.x().uniform(),
                   "focus plan does not advertise uniform fast path");

  CartesianGeometryPlan repeated;
  MeshPatch repeated_patch;
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, mesh, GeometryBudget{}, repeated,
                       repeated_patch)),
                   "focus plan compiles repeatedly");
  passed &= expect(repeated.fingerprint() == geometry.fingerprint(),
                   "focus coordinate generation is deterministic");
  return passed;
}

bool test_exact_cell_contract() {
  CartesianMeshSpec automatic = focus_mesh();
  CartesianGeometryPlan generated;
  MeshPatch generated_patch;
  const Status automatic_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, automatic, GeometryBudget{}, generated, generated_patch);
  bool passed = expect(static_cast<bool>(automatic_status),
                       "automatic focus count is available");
  if (!passed) {
    return false;
  }

  CartesianMeshSpec exact = automatic;
  exact.has_exact_cells = true;
  exact.exact_cells = generated.global_cells();
  CartesianGeometryPlan matched;
  MeshPatch matched_patch;
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, exact, GeometryBudget{}, matched,
                       matched_patch)),
                   "sufficient exact focus counts compile");
  passed &= expect(matched.global_cells().x == exact.exact_cells.x &&
                       matched.global_cells().y == exact.exact_cells.y &&
                       matched.global_cells().z == exact.exact_cells.z,
                   "exact cells are honored in all axes");

  exact.exact_cells.x = std::max(1, generated.global_cells().x / 2);
  CartesianGeometryPlan ignored;
  MeshPatch ignored_patch;
  const Status too_coarse = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, exact, GeometryBudget{}, ignored, ignored_patch);
  passed &= expect(too_coarse.code == StatusCode::invalid_plan,
                   "exact N below monitor mass is rejected");
  return passed;
}

bool test_bounded_generation() {
  CartesianMeshSpec excessive = focus_mesh();
  excessive.upper.x = 1.0e9;
  excessive.focus_regions.clear();
  excessive.base_spacing.x = 1.0e-6;
  excessive.minimum_spacing.x = 1.0e-6;
  CartesianGeometryPlan ignored;
  MeshPatch ignored_patch;
  const Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, excessive, GeometryBudget{}, ignored, ignored_patch);
  return expect(status.code == StatusCode::invalid_plan,
                "metric sampling is hard-capped before large allocation");
}

bool test_focus_peak_memory_gate() {
  CartesianMeshSpec limited = focus_mesh();
  limited.limits.max_memory_bytes_per_rank = 2048U;
  CartesianGeometryPlan ignored;
  MeshPatch ignored_patch;
  const Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, limited, GeometryBudget{}, ignored, ignored_patch);
  return expect(status.code == StatusCode::invalid_plan,
                "focus samples and scratch reject before peak allocation");
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_automatic_focus_grid();
  passed &= test_exact_cell_contract();
  passed &= test_bounded_generation();
  passed &= test_focus_peak_memory_gate();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 tensor-focus tests passed\n";
  return 0;
}
