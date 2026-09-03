// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_mesh.hpp"

#include "mesh_stl_scan_detail.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr StlScanBudget kBudget{UINT64_C(268435456), UINT64_C(536870912),
                                UINT64_C(4000000), UINT64_C(10000), 4U};
constexpr std::uint64_t kExpectedGlobalCrc = UINT64_C(6874342017398944037);

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-2.0, -2.0, -2.0};
  mesh.upper = {2.0, 2.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {16, 16, 16};
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(1000000), UINT64_C(1073741824)};
  return mesh;
}

std::uint64_t hash_byte(std::uint64_t hash, std::uint8_t value) noexcept {
  hash ^= value;
  return hash * UINT64_C(1099511628211);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS || size <= 0) {
    static_cast<void>(MPI_Abort(MPI_COMM_WORLD, 3));
    return 3;
  }
  const char* data_root = std::getenv("HUNDUN_V04_TEST_DATA");
  bool passed = data_root != nullptr;

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (passed) {
    passed = static_cast<bool>(CartesianGeometryCompiler::compile(
        MPI_COMM_WORLD, mesh_spec(), GeometryBudget{}, geometry, patch));
  }
  StlScanPlan plan;
  detail::reset_stl_open_count_for_test();
  if (passed) {
    const std::filesystem::path local_root =
        rank == 0 ? std::filesystem::path{data_root}
                  : std::filesystem::path{"/definitely/unreadable/nonroot"};
    passed = static_cast<bool>(StlScanCompiler::compile(
        MPI_COMM_WORLD, local_root, std::filesystem::path{"cube_binary.stl"},
        geometry, patch, CartesianAxis::y, kBudget, plan));
  }
  std::uint64_t global_open_count = detail::stl_open_count_for_test();
  if (MPI_Allreduce(MPI_IN_PLACE, &global_open_count, 1, MPI_UINT64_T, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    passed = false;
  }
  passed = passed && global_open_count == 1U;
  std::vector<std::uint8_t> local;
  if (passed) {
    const std::size_t count = static_cast<std::size_t>(patch.cells.x) *
                              static_cast<std::size_t>(patch.cells.y) *
                              static_cast<std::size_t>(patch.cells.z);
    local.resize(count);
    passed = static_cast<bool>(
        plan.classify(geometry, patch, {local.data(), local.size()}));
  }

  std::vector<std::uint8_t> global(16U * 16U * 16U, 0xffU);
  if (passed) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        for (std::int32_t x = 0; x < patch.cells.x; ++x) {
          const std::size_t local_index =
              static_cast<std::size_t>(x) +
              static_cast<std::size_t>(patch.cells.x) *
                  (static_cast<std::size_t>(y) +
                   static_cast<std::size_t>(patch.cells.y) *
                       static_cast<std::size_t>(z));
          const std::size_t global_index =
              static_cast<std::size_t>(patch.begin.x + x) + 16U *
                  (static_cast<std::size_t>(patch.begin.y + y) + 16U *
                       static_cast<std::size_t>(patch.begin.z + z));
          global[global_index] = local[local_index];
        }
      }
    }
  }
  if (MPI_Allreduce(MPI_IN_PLACE, global.data(), static_cast<int>(global.size()),
                    MPI_UINT8_T, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS) {
    passed = false;
  }
  std::uint64_t hash = UINT64_C(14695981039346656037);
  std::size_t solids = 0U;
  for (const std::uint8_t value : global) {
    hash = hash_byte(hash, value);
    solids += value == static_cast<std::uint8_t>(RegionFlag::solid) ? 1U : 0U;
  }
  std::uint64_t minimum_hash = hash;
  std::uint64_t maximum_hash = hash;
  if (MPI_Allreduce(MPI_IN_PLACE, &minimum_hash, 1, MPI_UINT64_T, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, &maximum_hash, 1, MPI_UINT64_T, MPI_MAX,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    passed = false;
  }
  passed = passed && minimum_hash == maximum_hash &&
           hash == kExpectedGlobalCrc && solids == 512U;
  int accepted = passed ? 1 : 0;
  if (MPI_Allreduce(MPI_IN_PLACE, &accepted, 1, MPI_INT, MPI_MIN,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    accepted = 0;
  }
  if (MPI_Finalize() != MPI_SUCCESS) {
    return 4;
  }
  if (accepted == 0) {
    return 1;
  }
  if (rank == 0) {
    std::cout << "v0.4 MPI STL classification CRC passed: " << hash << '\n';
  }
  return 0;
}
