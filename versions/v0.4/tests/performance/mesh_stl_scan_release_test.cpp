// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_mesh.hpp"

#include "mesh_stl_scan_detail.hpp"

#include <mpi.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

namespace allocation_observer {
std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void* allocate(std::size_t size) {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
  if (void* pointer = std::malloc(size == 0U ? 1U : size)) {
    return pointer;
  }
  throw std::bad_alloc{};
}

struct Guard {
  Guard() {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_relaxed);
  }
  ~Guard() { enabled.store(false, std::memory_order_relaxed); }
};
}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

constexpr StlScanBudget kBudget{UINT64_C(268435456), UINT64_C(536870912),
                                UINT64_C(4000000), UINT64_C(10000), 4U};
constexpr std::uint64_t kExpectedRegionCrc = UINT64_C(6874342017398944037);
constexpr auto kCatastrophicRegressionLimit = std::chrono::seconds(30);

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

std::uint64_t region_crc(const std::vector<std::uint8_t>& region) noexcept {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const std::uint8_t value : region) {
    hash ^= value;
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = true;
  const char* data_root = std::getenv("HUNDUN_V04_TEST_DATA");
  passed = data_root != nullptr;

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (passed) {
    passed = static_cast<bool>(CartesianGeometryCompiler::compile(
        MPI_COMM_SELF, mesh_spec(), GeometryBudget{}, geometry, patch));
  }

  detail::reset_stl_open_count_for_test();
  StlScanPlan plan;
  const auto compile_begin = std::chrono::steady_clock::now();
  if (passed) {
    passed = static_cast<bool>(StlScanCompiler::compile(
        MPI_COMM_SELF, data_root, std::filesystem::path{"cube_binary.stl"},
        geometry, patch, CartesianAxis::y, kBudget, plan));
  }
  const auto compile_elapsed = std::chrono::steady_clock::now() - compile_begin;
  passed = passed && detail::stl_open_count_for_test() == 1U &&
           compile_elapsed < kCatastrophicRegressionLimit;

  std::vector<std::uint8_t> region(16U * 16U * 16U);
  Status classification{};
  std::size_t allocations = std::numeric_limits<std::size_t>::max();
  if (passed) {
    allocation_observer::Guard guard;
    classification =
        plan.classify(geometry, patch, {region.data(), region.size()});
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed = passed && static_cast<bool>(classification) && allocations == 0U &&
           region_crc(region) == kExpectedRegionCrc;

  if (MPI_Finalize() != MPI_SUCCESS) {
    return 3;
  }
  if (!passed) {
    std::cerr << "v0.4 Release STL timing/load/allocation/CRC gate failed\n";
    return 1;
  }
  std::cout << "v0.4 Release STL gate passed in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   compile_elapsed)
                   .count()
            << " ms\n";
  return 0;
}
