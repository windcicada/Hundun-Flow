// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_mesh.hpp"

#include "mesh_stl_scan_detail.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <optional>
#include <string_view>
#include <type_traits>
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
                                UINT64_C(4000000), UINT64_C(10000), 1U};

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

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

bool compile_geometry(CartesianGeometryPlan& geometry, MeshPatch& patch) {
  return static_cast<bool>(CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh_spec(), GeometryBudget{}, geometry, patch));
}

std::array<TriangleInput, 12U> cube_triangles() {
  const Real3 p000{-1.0, -1.0, -1.0};
  const Real3 p001{-1.0, -1.0, 1.0};
  const Real3 p010{-1.0, 1.0, -1.0};
  const Real3 p011{-1.0, 1.0, 1.0};
  const Real3 p100{1.0, -1.0, -1.0};
  const Real3 p101{1.0, -1.0, 1.0};
  const Real3 p110{1.0, 1.0, -1.0};
  const Real3 p111{1.0, 1.0, 1.0};
  return {TriangleInput{p000, p001, p011},
          TriangleInput{p000, p011, p010},
          TriangleInput{p100, p110, p111},
          TriangleInput{p100, p111, p101},
          TriangleInput{p000, p100, p101},
          TriangleInput{p000, p101, p001},
          TriangleInput{p010, p011, p111},
          TriangleInput{p010, p111, p110},
          TriangleInput{p000, p010, p110},
          TriangleInput{p000, p110, p100},
          TriangleInput{p001, p101, p111},
          TriangleInput{p001, p111, p011}};
}

std::uint64_t region_hash(const std::vector<std::uint8_t>& region) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const std::uint8_t byte : region) {
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

std::size_t solid_count(const std::vector<std::uint8_t>& region) {
  std::size_t count = 0U;
  for (const std::uint8_t value : region) {
    count += value == static_cast<std::uint8_t>(RegionFlag::solid) ? 1U : 0U;
  }
  return count;
}

bool test_axis_permutations_and_order() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(geometry, patch),
                       "uniform scan geometry compiles");
  const auto cube = cube_triangles();
  std::uint64_t expected_hash = 0U;
  for (const CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    StlScanPlan plan;
    const Status status = StlScanCompiler::compile_triangles(
        geometry, patch, {cube.data(), cube.size()}, axis, kBudget, plan);
    passed &= expect(static_cast<bool>(status),
                     "closed cube compiles for each scan-axis permutation");
    std::vector<std::uint8_t> region(16U * 16U * 16U, 99U);
    passed &= expect(static_cast<bool>(plan.classify(
                         geometry, patch, {region.data(), region.size()})),
                     "precomputed cube spans classify");
    const std::uint64_t hash = region_hash(region);
    if (expected_hash == 0U) {
      expected_hash = hash;
    }
    passed &= expect(hash == expected_hash,
                     "scan-axis permutations produce identical region CRC");
    passed &= expect(solid_count(region) == 8U * 8U * 8U,
                     "cube interior cell count is exact");
  }

  std::array<TriangleInput, 12U> reversed{};
  for (std::size_t index = 0U; index < cube.size(); ++index) {
    const TriangleInput& source = cube[cube.size() - index - 1U];
    reversed[index] = source;
  }
  StlScanPlan reversed_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {reversed.data(), reversed.size()},
                       CartesianAxis::y, kBudget, reversed_plan)),
                   "reversed input triangle order compiles");
  std::vector<std::uint8_t> reversed_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(reversed_plan.classify(
                       geometry, patch,
                       {reversed_region.data(), reversed_region.size()})) &&
                       region_hash(reversed_region) == expected_hash,
                   "input triangle order leaves the region CRC unchanged");
  StlScanPlan original_order_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, kBudget, original_order_plan)) &&
                       reversed_plan.triangles().fingerprint() ==
                           original_order_plan.triangles().fingerprint(),
                   "TriangleSoA identity ignores only input triangle order");

  std::array<TriangleInput, 12U> reversed_winding = cube;
  for (TriangleInput& triangle : reversed_winding) {
    std::swap(triangle.b, triangle.c);
  }
  StlScanPlan reversed_winding_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch,
                       {reversed_winding.data(), reversed_winding.size()},
                       CartesianAxis::y, kBudget, reversed_winding_plan)) &&
                       reversed_winding_plan.triangles().fingerprint() !=
                           original_order_plan.triangles().fingerprint(),
                   "TriangleSoA identity preserves winding and normal orientation");
  std::vector<std::uint8_t> winding_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(reversed_winding_plan.classify(
                       geometry, patch,
                       {winding_region.data(), winding_region.size()})) &&
                       region_hash(winding_region) == expected_hash,
                   "global winding reversal preserves parity classification");
  return passed;
}

bool test_duplicate_edges_open_surface_and_publication() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(geometry, patch),
                       "duplicate/open test geometry compiles");
  auto cube = cube_triangles();
  StlScanPlan baseline_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, kBudget, baseline_plan)),
                   "baseline cube plan compiles");
  std::vector<std::uint8_t> baseline_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(baseline_plan.classify(
                       geometry, patch,
                       {baseline_region.data(), baseline_region.size()})),
                   "baseline cube classifies");
  std::vector<TriangleInput> duplicated(cube.begin(), cube.end());
  // cube[4] is a y-normal facet and intersects interior y scan lines.
  duplicated.push_back(cube[4]);
  StlScanPlan duplicate_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch,
                       {duplicated.data(), duplicated.size()}, CartesianAxis::y,
                       kBudget, duplicate_plan)),
                   "same-coordinate same-normal duplicate facets deduplicate");
  std::vector<std::uint8_t> duplicate_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(duplicate_plan.classify(
                       geometry, patch,
                       {duplicate_region.data(), duplicate_region.size()})) &&
                       region_hash(duplicate_region) ==
                           region_hash(baseline_region),
                   "intersecting same-coordinate same-sign duplicate leaves parity unchanged");

  StlScanPlan retained = std::move(duplicate_plan);
  const auto retained_fingerprint = retained.fingerprint();

  const TriangleInput coincident{{-10.0, 0.0, -10.0},
                                 {10.0, 0.0, -10.0},
                                 {0.0, 0.0, 20.0}};
  const std::array<TriangleInput, 2U> opposite{
      coincident,
      TriangleInput{coincident.a, coincident.c, coincident.b}};
  StlScanPlan opposite_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {opposite.data(), opposite.size()},
                       CartesianAxis::y, kBudget, opposite_plan)),
                   "same-coordinate opposite-normal intersections are retained as an even pair");
  std::vector<std::uint8_t> opposite_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(opposite_plan.classify(
                       geometry, patch,
                       {opposite_region.data(), opposite_region.size()})) &&
                       solid_count(opposite_region) == 0U,
                   "opposite-normal coincident pair forms only zero-width spans");
  const std::array<TriangleInput, 2U> open{cube[4], cube[5]};
  const Status open_status = StlScanCompiler::compile_triangles(
      geometry, patch, {open.data(), open.size()}, CartesianAxis::y, kBudget,
      retained);
  passed &= expect(open_status.code == StatusCode::invalid_plan,
                   "odd-parity open surface is rejected");
  passed &= expect(retained.fingerprint() == retained_fingerprint,
                   "rejected open surface preserves the published plan");

  const std::array<TriangleInput, 1U> degenerate{
      TriangleInput{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                    {2.0, 0.0, 0.0}}};
  passed &= expect(StlScanCompiler::compile_triangles(
                       geometry, patch,
                       {degenerate.data(), degenerate.size()}, CartesianAxis::x,
                       kBudget, retained)
                           .code == StatusCode::invalid_plan,
                   "degenerate triangle is rejected");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const std::array<TriangleInput, 1U> nonfinite{
      TriangleInput{{nan, 0.0, 0.0}, {1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0}}};
  passed &= expect(StlScanCompiler::compile_triangles(
                       geometry, patch, {nonfinite.data(), nonfinite.size()},
                       CartesianAxis::z, kBudget, retained)
                           .code == StatusCode::invalid_plan,
                   "non-finite triangle is rejected");
  return passed;
}

bool test_inclusive_barycentric_and_transverse_reject() {
  CartesianMeshSpec mesh = mesh_spec();
  mesh.lower = {-1.5, -1.5, -1.5};
  mesh.upper = {1.5, 1.5, 1.5};
  mesh.exact_cells = {3, 3, 3};
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_SELF, mesh, GeometryBudget{}, geometry,
                           patch)),
                       "inclusive barycentric geometry compiles");
  const std::array<TriangleInput, 2U> slab{
      TriangleInput{{-1.0, -0.5, -1.0}, {1.0, -0.5, -1.0},
                    {-1.0, -0.5, 1.0}},
      TriangleInput{{-1.0, 0.5, -1.0}, {-1.0, 0.5, 1.0},
                    {1.0, 0.5, -1.0}}};
  StlScanPlan plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {slab.data(), slab.size()},
                       CartesianAxis::y, kBudget, plan)),
                   "parallel triangular slab compiles");
  std::vector<std::uint8_t> region(27U);
  passed &= expect(static_cast<bool>(
                       plan.classify(geometry, patch,
                                     {region.data(), region.size()})),
                   "triangular slab classifies");
  const auto at = [&](std::size_t x, std::size_t y, std::size_t z) {
    return region[x + 3U * (y + 3U * z)];
  };
  passed &= expect(at(0U, 1U, 0U) ==
                           static_cast<std::uint8_t>(RegionFlag::solid),
                     "scan through an exact triangle vertex is inclusive");
  passed &= expect(at(1U, 1U, 0U) ==
                           static_cast<std::uint8_t>(RegionFlag::solid) &&
                       at(0U, 1U, 1U) ==
                           static_cast<std::uint8_t>(RegionFlag::solid),
                   "scan through exact triangle edges is inclusive");
  passed &= expect(at(2U, 1U, 2U) ==
                           static_cast<std::uint8_t>(RegionFlag::fluid),
                   "transverse AABB/barycentric reject leaves outside line fluid");
  return passed;
}

bool test_hot_classify_and_empty_file_plan() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(geometry, patch),
                       "hot-path test geometry compiles");
  const auto cube = cube_triangles();
  StlScanPlan plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, kBudget, plan)),
                   "hot-path cube plan compiles");
  std::vector<std::uint8_t> region(16U * 16U * 16U);
  Status classify_status{};
  std::size_t allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    classify_status =
        plan.classify(geometry, patch, {region.data(), region.size()});
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(classify_status) && allocations == 0U,
                   "hot classify performs zero dynamic allocations");

  StlScanPlan empty;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile(
                       MPI_COMM_SELF, ".", std::nullopt, geometry, patch,
                       CartesianAxis::y, kBudget, empty)),
                   "missing optional STL compiles an empty plan");
  std::fill(region.begin(), region.end(), 0U);
  passed &= expect(static_cast<bool>(empty.classify(
                       geometry, patch, {region.data(), region.size()})) &&
                       solid_count(region) == 0U,
                   "empty STL plan classifies every cell as fluid");
  return passed;
}

bool test_budget_and_public_plan_contract() {
  static_assert(!std::is_copy_constructible_v<TriangleSoA>);
  static_assert(std::is_nothrow_move_constructible_v<TriangleSoA>);
  static_assert(!std::is_copy_constructible_v<StlScanPlan>);
  static_assert(std::is_nothrow_move_constructible_v<StlScanPlan>);

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(geometry, patch),
                       "budget test geometry compiles");
  const auto cube = cube_triangles();
  StlScanPlan retained;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, kBudget, retained)),
                   "baseline budgeted plan compiles");
  const PlanFingerprint fingerprint = retained.fingerprint();
  const StlScanBudget tiny_memory{1024U, 2048U, 1000U, 100U, 1U};
  passed &= expect(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, tiny_memory, retained)
                           .code == StatusCode::invalid_plan &&
                       retained.fingerprint() == fingerprint,
                   "persistent/build peak budget rejects before publication");
  StlScanBudget tiny_bins = kBudget;
  tiny_bins.max_bin_references = 1U;
  passed &= expect(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, tiny_bins, retained)
                           .code == StatusCode::invalid_plan &&
                       retained.fingerprint() == fingerprint,
                   "bin reference overflow rejects without brute-force fallback");
  StlScanBudget tiny_events = kBudget;
  tiny_events.max_events_per_line = 1U;
  passed &= expect(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, tiny_events, retained)
                           .code == StatusCode::invalid_plan &&
                       retained.fingerprint() == fingerprint,
                   "per-line scratch bound rejects before vector growth");
  // For this 12-triangle, 256-line plan the prior estimate was below 54,000
  // bytes, while charging the final and worker span stores during merge puts
  // the conservative peak above it.  This catches omission of either copy.
  const StlScanBudget merge_peak_boundary{20000U, 54000U, 1000U, 100U, 4U};
  passed &= expect(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, merge_peak_boundary, retained)
                           .code == StatusCode::invalid_plan &&
                       retained.fingerprint() == fingerprint,
                   "merge peak accounts for coexisting worker and final span storage");
  passed &= expect(retained.triangles().ax().size == cube.size() &&
                       retained.triangles().nx().size == cube.size() &&
                       retained.triangles().max_z().size == cube.size(),
                   "TriangleSoA exposes immutable public spans");

  StlScanBudget parallel_budget = kBudget;
  parallel_budget.worker_threads = 4U;
  StlScanPlan parallel_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {cube.data(), cube.size()},
                       CartesianAxis::y, parallel_budget, parallel_plan)) &&
                       parallel_plan.fingerprint() == retained.fingerprint(),
                   "fixed one-thread and four-thread chunks produce bitwise-identical plans");
  std::vector<std::uint8_t> serial_region(16U * 16U * 16U);
  std::vector<std::uint8_t> parallel_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(retained.classify(
                       geometry, patch,
                       {serial_region.data(), serial_region.size()})) &&
                       static_cast<bool>(parallel_plan.classify(
                           geometry, patch,
                           {parallel_region.data(), parallel_region.size()})) &&
                       region_hash(serial_region) == region_hash(parallel_region),
                   "fixed thread chunks preserve classification CRC");

  const std::size_t maximum = std::numeric_limits<std::size_t>::max();
  const std::array<detail::LineChunk, 3U> extreme_chunks{
      detail::fixed_line_chunk(maximum, 0U, 3U),
      detail::fixed_line_chunk(maximum, 1U, 3U),
      detail::fixed_line_chunk(maximum, 2U, 3U)};
  const std::size_t quotient = maximum / 3U;
  const std::size_t remainder = maximum % 3U;
  passed &= expect(extreme_chunks[0].begin == 0U &&
                       extreme_chunks[0].end ==
                           quotient + (remainder > 0U ? 1U : 0U) &&
                       extreme_chunks[1].begin == extreme_chunks[0].end &&
                       extreme_chunks[1].end - extreme_chunks[1].begin ==
                           quotient + (remainder > 1U ? 1U : 0U) &&
                       extreme_chunks[2].begin == extreme_chunks[1].end &&
                       extreme_chunks[2].end - extreme_chunks[2].begin ==
                           quotient + (remainder > 2U ? 1U : 0U) &&
                       extreme_chunks[2].end == maximum,
                   "fixed chunk arithmetic covers SIZE_MAX without overflow");

  const PlanFingerprint before_launch_failure = retained.fingerprint();
  detail::fail_thread_launch_after_for_test(1U);
  const Status launch_failure = StlScanCompiler::compile_triangles(
      geometry, patch, {cube.data(), cube.size()}, CartesianAxis::y,
      parallel_budget, retained);
  detail::reset_thread_launch_failure_for_test();
  passed &= expect(launch_failure.code == StatusCode::allocation_failure &&
                       retained.fingerprint() == before_launch_failure,
                   "partial thread launch failure joins workers and preserves publication");
  return passed;
}

bool test_file_formats() {
  const char* data_root = std::getenv("HUNDUN_V04_TEST_DATA");
  if (data_root == nullptr) {
    return expect(false, "HUNDUN_V04_TEST_DATA is configured");
  }
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(geometry, patch),
                       "file-format test geometry compiles");
  StlScanPlan ascii_plan;
  StlScanPlan binary_plan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile(
                       MPI_COMM_SELF, data_root,
                       std::filesystem::path{"cylinder_ascii.stl"}, geometry,
                       patch, CartesianAxis::z, kBudget, ascii_plan)),
                   "ASCII STL loads once and compiles");
  passed &= expect(static_cast<bool>(StlScanCompiler::compile(
                       MPI_COMM_SELF, data_root,
                       std::filesystem::path{"cube_binary.stl"}, geometry,
                       patch, CartesianAxis::y, kBudget, binary_plan)),
                   "binary STL loads once and compiles");
  passed &= expect(ascii_plan.triangle_count() > 0U &&
                       binary_plan.triangle_count() == 12U,
                   "ASCII and binary parsers publish triangle SoA counts");
  std::vector<std::uint8_t> ascii_region(16U * 16U * 16U);
  std::vector<std::uint8_t> binary_region(16U * 16U * 16U);
  passed &= expect(static_cast<bool>(ascii_plan.classify(
                       geometry, patch,
                       {ascii_region.data(), ascii_region.size()})) &&
                       static_cast<bool>(binary_plan.classify(
                           geometry, patch,
                           {binary_region.data(), binary_region.size()})) &&
                       region_hash(ascii_region) == region_hash(binary_region) &&
                       ascii_plan.triangles().fingerprint() ==
                           binary_plan.triangles().fingerprint(),
                   "ASCII/binary encodings of the same cube have identical geometry and region identity");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_axis_permutations_and_order();
  passed &= test_duplicate_edges_open_surface_and_publication();
  passed &= test_inclusive_barycentric_and_transverse_reject();
  passed &= test_hot_classify_and_empty_file_plan();
  passed &= test_budget_and_public_plan_contract();
  passed &= test_file_formats();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 STL scan tests passed\n";
  return 0;
}
