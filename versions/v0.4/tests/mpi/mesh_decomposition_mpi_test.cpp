// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_mesh.hpp"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-5.0, -2.0, 0.0};
  mesh.upper = {7.0, 8.0, 3.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {13, 10, 3};
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {100000U, 1U << 30U};
  return mesh;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  CartesianMeshSpec local_mesh = mesh_spec();
  if (rank != 0) {
    local_mesh.kind = GeometryKind::tensor_stretched;
    local_mesh.lower = {-900.0, -800.0, -700.0};
    local_mesh.upper = {-600.0, -500.0, -400.0};
  }
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  const Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_WORLD, local_mesh, GeometryBudget{4096U, 64U}, geometry, patch);
  bool passed = expect(static_cast<bool>(status), rank,
                       "collective geometry/decomposition compiles");

  std::uint64_t minimum_fingerprint = geometry.fingerprint();
  std::uint64_t maximum_fingerprint = geometry.fingerprint();
  MPI_Allreduce(MPI_IN_PLACE, &minimum_fingerprint, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &maximum_fingerprint, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(minimum_fingerprint == maximum_fingerprint, rank,
                   "all ranks receive identical global geometry");
  passed &= expect(geometry.kind() == GeometryKind::uniform &&
                       geometry.lower().x == -5.0 &&
                       geometry.lower().y == -2.0 &&
                       geometry.lower().z == 0.0 &&
                       geometry.upper().x == 7.0 &&
                       geometry.upper().y == 8.0 &&
                       geometry.upper().z == 3.0,
                   rank,
                   "rank-zero geometry metadata overrides poisoned non-root input");

  CartesianGeometryPlan retained = std::move(geometry);
  MeshPatch retained_patch = patch;
  const auto retained_fingerprint = retained.fingerprint();
  const GeometryBudget asymmetric_budget{
      rank == size - 1 ? (1U << 30U) : 4096U, 64U};
  const Status rejected = CartesianGeometryCompiler::compile(
      MPI_COMM_WORLD, mesh_spec(), asymmetric_budget, retained, retained_patch);
  int rejected_code = static_cast<int>(rejected.code);
  int minimum_rejected_code = rejected_code;
  int maximum_rejected_code = rejected_code;
  MPI_Allreduce(MPI_IN_PLACE, &minimum_rejected_code, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &maximum_rejected_code, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       minimum_rejected_code == maximum_rejected_code,
                   rank, "asymmetric local memory failure reaches consensus");
  passed &= expect(retained.fingerprint() == retained_fingerprint &&
                       retained_patch.begin.x == patch.begin.x &&
                       retained_patch.begin.y == patch.begin.y &&
                       retained_patch.begin.z == patch.begin.z,
                   rank, "collective rejection preserves published outputs");

  const std::array<std::int32_t, 12U> local{
      patch.begin.x,       patch.begin.y,       patch.begin.z,
      patch.cells.x,       patch.cells.y,       patch.cells.z,
      patch.process_grid.x, patch.process_grid.y, patch.process_grid.z,
      patch.process_coord.x, patch.process_coord.y, patch.process_coord.z};
  std::vector<std::int32_t> gathered;
  if (rank == 0) {
    gathered.resize(static_cast<std::size_t>(size) * local.size());
  }
  MPI_Gather(local.data(), static_cast<int>(local.size()), MPI_INT32_T,
             gathered.data(), static_cast<int>(local.size()), MPI_INT32_T, 0,
             MPI_COMM_WORLD);
  if (rank == 0) {
    std::vector<std::uint8_t> ownership(13U * 10U * 3U, 0U);
    for (int owner = 0; owner < size; ++owner) {
      const auto* record = gathered.data() +
                           static_cast<std::size_t>(owner) * local.size();
      passed &= expect(record[6] * record[7] * record[8] == size, rank,
                       "factor enumeration produces the exact rank product");
      passed &= expect(record[9] >= 0 && record[9] < record[6] &&
                           record[10] >= 0 && record[10] < record[7] &&
                           record[11] >= 0 && record[11] < record[8],
                       rank, "rank has a valid deterministic process coordinate");
      for (std::int32_t z = 0; z < record[5]; ++z) {
        for (std::int32_t y = 0; y < record[4]; ++y) {
          for (std::int32_t x = 0; x < record[3]; ++x) {
            const std::size_t index =
                static_cast<std::size_t>(record[0] + x) +
                13U * (static_cast<std::size_t>(record[1] + y) +
                       10U * static_cast<std::size_t>(record[2] + z));
            ++ownership[index];
          }
        }
      }
    }
    for (const std::uint8_t owners : ownership) {
      passed &= expect(owners == 1U, rank,
                       "q/r patches cover every global cell exactly once");
    }
  }

  int local_passed = passed ? 1 : 0;
  int all_passed = 0;
  MPI_Allreduce(&local_passed, &all_passed, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Finalize();
  if (all_passed == 0) {
    return 1;
  }
  if (rank == 0) {
    std::cout << "v0.4 MPI decomposition tests passed\n";
  }
  return 0;
}
