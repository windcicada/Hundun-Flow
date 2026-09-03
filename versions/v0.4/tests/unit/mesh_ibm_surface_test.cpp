// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <utility>
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
void operator delete(void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
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

bool near(double lhs, double rhs, double tolerance = 2.0e-13) {
  return std::isfinite(lhs) && std::isfinite(rhs) &&
         std::abs(lhs - rhs) <= tolerance;
}

bool near(Real3 lhs, Real3 rhs, double tolerance = 2.0e-13) {
  return near(lhs.x, rhs.x, tolerance) &&
         near(lhs.y, rhs.y, tolerance) &&
         near(lhs.z, rhs.z, tolerance);
}

Real3 add(Real3 lhs, Real3 rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 multiply(double scalar, Real3 value) {
  return {scalar * value.x, scalar * value.y, scalar * value.z};
}

double dot(Real3 lhs, Real3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-3.0, -3.0, -3.0};
  mesh.upper = {3.0, 3.0, 3.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {24, 24, 24};
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(1000000), UINT64_C(1073741824)};
  return mesh;
}

bool compile_geometry(CartesianGeometryPlan& geometry, MeshPatch& patch) {
  return static_cast<bool>(CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh_spec(), GeometryBudget{}, geometry, patch));
}

void append_box(std::vector<TriangleInput>& triangles, Real3 lower,
                Real3 upper) {
  const Real3 p000{lower.x, lower.y, lower.z};
  const Real3 p001{lower.x, lower.y, upper.z};
  const Real3 p010{lower.x, upper.y, lower.z};
  const Real3 p011{lower.x, upper.y, upper.z};
  const Real3 p100{upper.x, lower.y, lower.z};
  const Real3 p101{upper.x, lower.y, upper.z};
  const Real3 p110{upper.x, upper.y, lower.z};
  const Real3 p111{upper.x, upper.y, upper.z};
  const std::array<TriangleInput, 12U> box{
      TriangleInput{p000, p001, p011},
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
  triangles.insert(triangles.end(), box.begin(), box.end());
}

std::vector<TriangleInput> cube_triangles() {
  std::vector<TriangleInput> triangles;
  append_box(triangles, {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
  return triangles;
}

bool compile_scan(const CartesianGeometryPlan& geometry,
                  const MeshPatch& patch,
                  const std::vector<TriangleInput>& triangles,
                  StlScanPlan& scan) {
  return static_cast<bool>(StlScanCompiler::compile_triangles(
      geometry, patch, {triangles.data(), triangles.size()},
      CartesianAxis::y, kBudget, scan));
}

bool same_triangle(const SurfaceTriangle& lhs,
                   const SurfaceTriangle& rhs) {
  if (lhs.id != rhs.id || !near(lhs.area, rhs.area, 0.0) ||
      !near(lhs.centroid, rhs.centroid, 0.0) ||
      !near(lhs.geometric_outward_normal,
            rhs.geometric_outward_normal, 0.0)) {
    return false;
  }
  for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
    if (!near(lhs.vertices[vertex], rhs.vertices[vertex], 0.0)) {
      return false;
    }
  }
  return true;
}

bool same_plan(const ImmersedSurfacePlan& lhs,
               const ImmersedSurfacePlan& rhs) {
  if (lhs.fingerprint() != rhs.fingerprint() ||
      lhs.triangles().size != rhs.triangles().size ||
      !near(lhs.closed_volume(), rhs.closed_volume(), 0.0) ||
      !near(lhs.bounding_box_min(), rhs.bounding_box_min(), 0.0) ||
      !near(lhs.bounding_box_max(), rhs.bounding_box_max(), 0.0)) {
    return false;
  }
  for (std::size_t index = 0U; index < lhs.triangles().size; ++index) {
    if (!same_triangle(lhs.triangles().data[index],
                       rhs.triangles().data[index])) {
      return false;
    }
  }
  return true;
}

bool test_canonical_closed_surface_and_quadrature_basis(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    ImmersedSurfacePlan& baseline) {
  const std::vector<TriangleInput> cube = cube_triangles();
  StlScanPlan scan;
  bool passed = expect(compile_scan(geometry, patch, cube, scan),
                       "closed cube scan compiles");
  passed &= expect(static_cast<bool>(
                       ImmersedSurfaceCompiler::compile(scan, baseline)),
                   "closed cube surface compiles");
  passed &= expect(baseline.triangles().size == 12U,
                   "cube publishes twelve canonical triangles");
  passed &= expect(near(baseline.bounding_box_min(), {-1.0, -1.0, -1.0}) &&
                       near(baseline.bounding_box_max(), {1.0, 1.0, 1.0}),
                   "cube publishes its exact bounding box");
  passed &= expect(near(baseline.closed_volume(), 8.0),
                   "cube publishes positive closed volume");

  double total_area = 0.0;
  for (std::size_t index = 0U; index < baseline.triangles().size; ++index) {
    const SurfaceTriangle& triangle = baseline.triangles().data[index];
    passed &= expect(triangle.id == index,
                     "canonical triangle ids are dense and ordered");
    passed &= expect(near(triangle.area, 2.0),
                     "each cube triangle has exact area");
    passed &= expect(dot(triangle.geometric_outward_normal,
                         triangle.centroid) > 0.0,
                     "canonical normals point out of the closed cube");
    total_area += triangle.area;

    const std::array<std::array<double, 3U>, 3U> barycentric{
        std::array<double, 3U>{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
        std::array<double, 3U>{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
        std::array<double, 3U>{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}};
    Real3 first_moment{};
    double x_squared_integral = 0.0;
    for (const auto& lambda : barycentric) {
      Real3 point{};
      for (std::size_t vertex = 0U; vertex < 3U; ++vertex) {
        point = add(point,
                    multiply(lambda[vertex], triangle.vertices[vertex]));
      }
      first_moment =
          add(first_moment, multiply(triangle.area / 3.0, point));
      x_squared_integral += triangle.area / 3.0 * point.x * point.x;
    }
    const double x0 = triangle.vertices[0].x;
    const double x1 = triangle.vertices[1].x;
    const double x2 = triangle.vertices[2].x;
    const double exact_x_squared =
        triangle.area / 6.0 *
        (x0 * x0 + x1 * x1 + x2 * x2 + x0 * x1 + x0 * x2 +
         x1 * x2);
    passed &= expect(near(first_moment,
                          multiply(triangle.area, triangle.centroid)),
                     "three-point triangle rule reproduces first moments");
    passed &= expect(near(x_squared_integral, exact_x_squared),
                     "three-point triangle rule reproduces quadratics");
  }
  passed &= expect(near(total_area, 24.0),
                   "canonical cube surface area is exact");
  return passed;
}

bool test_order_winding_and_cycle_invariance(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    const ImmersedSurfacePlan& baseline) {
  bool passed = true;
  std::vector<TriangleInput> reordered = cube_triangles();
  std::reverse(reordered.begin(), reordered.end());
  StlScanPlan reordered_scan;
  ImmersedSurfacePlan reordered_surface;
  passed &= expect(compile_scan(geometry, patch, reordered, reordered_scan) &&
                       static_cast<bool>(ImmersedSurfaceCompiler::compile(
                           reordered_scan, reordered_surface)) &&
                       same_plan(baseline, reordered_surface),
                   "facet input order cannot change canonical surface bytes");

  std::vector<TriangleInput> cycled = cube_triangles();
  for (TriangleInput& triangle : cycled) {
    const TriangleInput source = triangle;
    triangle = {source.b, source.c, source.a};
  }
  StlScanPlan cycled_scan;
  ImmersedSurfacePlan cycled_surface;
  passed &= expect(compile_scan(geometry, patch, cycled, cycled_scan) &&
                       static_cast<bool>(ImmersedSurfaceCompiler::compile(
                           cycled_scan, cycled_surface)) &&
                       same_plan(baseline, cycled_surface),
                   "cyclic facet vertices cannot change canonical surface");

  std::vector<TriangleInput> inward = cube_triangles();
  for (TriangleInput& triangle : inward) {
    std::swap(triangle.b, triangle.c);
  }
  StlScanPlan inward_scan;
  ImmersedSurfacePlan inward_surface;
  passed &= expect(compile_scan(geometry, patch, inward, inward_scan) &&
                       static_cast<bool>(ImmersedSurfaceCompiler::compile(
                           inward_scan, inward_surface)) &&
                       same_plan(baseline, inward_surface),
                   "whole-surface inward winding is canonicalized outward");
  return passed;
}

bool test_disconnected_components_use_local_volume_scale(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch) {
  std::vector<TriangleInput> triangles;
  constexpr double width = 1.0e-4;
  append_box(triangles, {-2.0, -2.0, -2.0},
             {-2.0 + width, -2.0 + width, -2.0 + width});
  append_box(triangles, {2.0 - width, 2.0 - width, 2.0 - width},
             {2.0, 2.0, 2.0});
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  const bool compiled = compile_scan(geometry, patch, triangles, scan) &&
                        static_cast<bool>(
                            ImmersedSurfaceCompiler::compile(scan, surface));
  bool passed = expect(
      compiled,
      "distant tiny closed components use component-local volume scales");
  if (compiled) {
    passed &= expect(near(surface.closed_volume(), 2.0e-12, 2.0e-23),
                     "disconnected component volumes are summed positively");
  }
  return passed;
}

bool test_shared_stl_vertices_survive_edge_storage_roundtrip(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch) {
  // These are exactly representable binary32 STL coordinates promoted to
  // binary64.  Reconstructing b/c as a + (b/c - a) changes the near-zero x
  // coordinate, so topology must consume the preserved source vertices.
  const Real3 a{0.03270156309008598, 0.4989294707775116, -1.0};
  const Real3 b{3.0616171314629196e-17, 0.5, -1.0};
  const Real3 c{3.0616171314629196e-17, 0.5, 0.0};
  const Real3 d{0.1, 0.4, -0.5};
  const std::vector<TriangleInput> tetrahedron{
      {a, b, c}, {a, d, b}, {a, c, d}, {b, d, c}};
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  const bool compiled = compile_scan(geometry, patch, tetrahedron, scan) &&
                        static_cast<bool>(
                            ImmersedSurfaceCompiler::compile(scan, surface));
  bool passed = expect(
      compiled,
      "shared STL vertices remain bit-identical through the scan plan");
  if (compiled) {
    passed &= expect(surface.triangles().size == tetrahedron.size(),
                     "roundtrip fixture publishes all closed facets");
  }
  return passed;
}

bool test_topology_rejection_and_atomic_publication(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    ImmersedSurfacePlan& baseline) {
  bool passed = true;
  const PlanFingerprint fingerprint = baseline.fingerprint();
  const double volume = baseline.closed_volume();

  std::vector<TriangleInput> nested;
  append_box(nested, {-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
  append_box(nested, {-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5});
  StlScanPlan nested_scan;
  passed &= expect(compile_scan(geometry, patch, nested, nested_scan),
                   "scan stage accepts concentric closed-shell fixture");
  const Status nested_status =
      ImmersedSurfaceCompiler::compile(nested_scan, baseline);
  passed &= expect(!nested_status &&
                       nested_status.code == StatusCode::invalid_plan,
                   "surface stage rejects nested closed components");
  passed &= expect(baseline.fingerprint() == fingerprint &&
                       near(baseline.closed_volume(), volume, 0.0),
                   "nested-component failure leaves published plan unchanged");

  std::vector<TriangleInput> local_flip = cube_triangles();
  std::swap(local_flip[0].b, local_flip[0].c);
  StlScanPlan local_flip_scan;
  passed &= expect(compile_scan(geometry, patch, local_flip, local_flip_scan),
                   "scan stage accepts local winding mutation fixture");
  const Status local_flip_status =
      ImmersedSurfaceCompiler::compile(local_flip_scan, baseline);
  passed &= expect(!local_flip_status &&
                       local_flip_status.code == StatusCode::invalid_plan,
                   "surface stage rejects a locally flipped facet");
  passed &= expect(baseline.fingerprint() == fingerprint &&
                       near(baseline.closed_volume(), volume, 0.0),
                   "local flip failure leaves published plan unchanged");

  std::vector<TriangleInput> open = cube_triangles();
  open.erase(open.begin());
  StlScanPlan open_scan;
  passed &= expect(compile_scan(geometry, patch, open, open_scan),
                   "scan stage accepts non-scan-normal open fixture");
  const Status open_status = ImmersedSurfaceCompiler::compile(open_scan,
                                                               baseline);
  passed &= expect(!open_status && open_status.code == StatusCode::invalid_plan,
                   "surface stage rejects an open edge");

  std::vector<TriangleInput> duplicate = cube_triangles();
  duplicate.push_back(duplicate.front());
  StlScanPlan duplicate_scan;
  passed &= expect(compile_scan(geometry, patch, duplicate, duplicate_scan),
                   "scan stage accepts non-scan-normal duplicate fixture");
  const Status duplicate_status =
      ImmersedSurfaceCompiler::compile(duplicate_scan, baseline);
  passed &= expect(!duplicate_status &&
                       duplicate_status.code == StatusCode::invalid_plan,
                   "surface stage rejects duplicate/non-manifold facets");

  std::vector<TriangleInput> bow_tie;
  append_box(bow_tie, {-1.0, -1.0, -1.0}, {0.0, 0.0, 0.0});
  append_box(bow_tie, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  StlScanPlan bow_tie_scan;
  passed &= expect(compile_scan(geometry, patch, bow_tie, bow_tie_scan),
                   "scan stage accepts vertex-touching closed shells");
  const Status bow_tie_status =
      ImmersedSurfaceCompiler::compile(bow_tie_scan, baseline);
  passed &= expect(!bow_tie_status &&
                       bow_tie_status.code == StatusCode::invalid_plan,
                   "surface stage rejects a non-manifold bow-tie vertex");
  passed &= expect(baseline.fingerprint() == fingerprint &&
                       near(baseline.closed_volume(), volume, 0.0),
                   "all topology failures are publication-atomic");
  return passed;
}

bool test_queries_and_no_hot_allocations(const ImmersedSurfacePlan& surface) {
  bool passed = true;
  ClosestSurfacePoint closest{};
  passed &= expect(static_cast<bool>(surface.closest_point(
                       {2.0, 0.25, 0.10}, closest)) &&
                       near(closest.point, {1.0, 0.25, 0.10}) &&
                       near(closest.squared_distance, 1.0) &&
                       near(closest.geometric_outward_normal,
                            {1.0, 0.0, 0.0}),
                   "closest-point query finds the analytical cube face");
  passed &= expect(static_cast<bool>(surface.closest_point(
                       {0.80, 0.20, 0.10}, closest)) &&
                       near(closest.point, {1.0, 0.20, 0.10}) &&
                       near(closest.squared_distance, 0.04),
                   "closest-point query works inside the closed body");
  passed &= expect(static_cast<bool>(surface.closest_point(
                       {0.0, 0.0, 0.0}, closest)) &&
                       closest.triangle == 0U &&
                       near(closest.squared_distance, 1.0),
                   "equal-distance closest-point ties choose minimum id");

  SurfaceSegmentIntersection hit{};
  passed &= expect(static_cast<bool>(surface.first_segment_intersection(
                       {-2.0, 0.20, 0.10}, {2.0, 0.20, 0.10}, hit)) &&
                       hit.triangle != kInvalidSurfaceTriangle &&
                       near(hit.point, {-1.0, 0.20, 0.10}) &&
                       near(hit.segment_fraction, 0.25),
                   "segment query publishes the first forward intersection");
  passed &= expect(static_cast<bool>(surface.first_segment_intersection(
                       {2.0, 0.20, 0.10}, {-2.0, 0.20, 0.10}, hit)) &&
                       near(hit.point, {1.0, 0.20, 0.10}) &&
                       near(hit.segment_fraction, 0.25),
                   "segment direction changes which surface is first");

  hit = {17U, {9.0, 9.0, 9.0}, 0.75};
  passed &= expect(static_cast<bool>(surface.first_segment_intersection(
                       {-2.0, 2.0, 0.0}, {2.0, 2.0, 0.0}, hit)) &&
                       hit.triangle == kInvalidSurfaceTriangle &&
                       near(hit.point, {0.0, 0.0, 0.0}, 0.0) &&
                       near(hit.segment_fraction, 0.0, 0.0),
                   "valid no-hit segment resets the output deterministically");

  ClosestSurfacePoint failure_sentinel{7U, {8.0, 8.0, 8.0},
                                       {3.0, 4.0, 5.0}, 6.0};
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const Status non_finite =
      surface.closest_point({nan, 0.0, 0.0}, failure_sentinel);
  passed &= expect(!non_finite && failure_sentinel.triangle == 7U &&
                       near(failure_sentinel.point, {8.0, 8.0, 8.0}, 0.0),
                   "invalid closest-point query is output-atomic");
  hit = {19U, {7.0, 7.0, 7.0}, 0.5};
  const Status zero_length =
      surface.first_segment_intersection({0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}, hit);
  passed &= expect(!zero_length && hit.triangle == 19U &&
                       near(hit.point, {7.0, 7.0, 7.0}, 0.0),
                   "zero-length segment failure is output-atomic");

  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < 1000U; ++iteration) {
      ClosestSurfacePoint local_closest{};
      SurfaceSegmentIntersection local_hit{};
      const double offset = static_cast<double>(iteration % 17U) * 0.01;
      if (!surface.closest_point({2.0, offset, -offset}, local_closest) ||
          !surface.first_segment_intersection(
              {-2.0, offset, -offset}, {2.0, offset, -offset}, local_hit)) {
        passed = false;
        break;
      }
    }
  }
  passed &= expect(allocation_observer::count.load(std::memory_order_relaxed) ==
                       0U,
                   "closest/segment BVH queries allocate no heap storage");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(compile_geometry(geometry, patch),
                       "surface test geometry compiles");
  ImmersedSurfacePlan baseline;
  if (passed) {
    passed &= test_canonical_closed_surface_and_quadrature_basis(
        geometry, patch, baseline);
  }
  if (baseline.fingerprint() != 0U) {
    passed &= test_order_winding_and_cycle_invariance(geometry, patch,
                                                       baseline);
    passed &= test_disconnected_components_use_local_volume_scale(geometry,
                                                                   patch);
    passed &= test_shared_stl_vertices_survive_edge_storage_roundtrip(geometry,
                                                                      patch);
    passed &= test_topology_rejection_and_atomic_publication(geometry, patch,
                                                              baseline);
    passed &= test_queries_and_no_hot_allocations(baseline);
  }
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
