// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr StlScanBudget kScanBudget{UINT64_C(268435456),
                                    UINT64_C(536870912),
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
  mesh.lower = {-1.5, -1.5, -1.5};
  mesh.upper = {1.5, 1.5, 1.5};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {24, 24, 24};
  mesh.minimum_spacing = {1.0e-9, 1.0e-9, 1.0e-9};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(1000000), UINT64_C(1073741824)};
  return mesh;
}

std::array<TriangleInput, 12U> cube() {
  const Real3 a{-0.5, -0.5, -0.5};
  const Real3 b{-0.5, -0.5, 0.5};
  const Real3 c{-0.5, 0.5, -0.5};
  const Real3 d{-0.5, 0.5, 0.5};
  const Real3 e{0.5, -0.5, -0.5};
  const Real3 f{0.5, -0.5, 0.5};
  const Real3 g{0.5, 0.5, -0.5};
  const Real3 h{0.5, 0.5, 0.5};
  return {TriangleInput{a, d, b}, TriangleInput{a, c, d},
          TriangleInput{e, f, h}, TriangleInput{e, h, g},
          TriangleInput{a, b, f}, TriangleInput{a, f, e},
          TriangleInput{c, g, h}, TriangleInput{c, h, d},
          TriangleInput{a, e, g}, TriangleInput{a, g, c},
          TriangleInput{b, d, h}, TriangleInput{b, h, f}};
}

std::array<TriangleInput, 12U> box(Real3 lower, Real3 upper) {
  const Real3 a{lower.x, lower.y, lower.z};
  const Real3 b{lower.x, lower.y, upper.z};
  const Real3 c{lower.x, upper.y, lower.z};
  const Real3 d{lower.x, upper.y, upper.z};
  const Real3 e{upper.x, lower.y, lower.z};
  const Real3 f{upper.x, lower.y, upper.z};
  const Real3 g{upper.x, upper.y, lower.z};
  const Real3 h{upper.x, upper.y, upper.z};
  return {TriangleInput{a, d, b}, TriangleInput{a, c, d},
          TriangleInput{e, f, h}, TriangleInput{e, h, g},
          TriangleInput{a, b, f}, TriangleInput{a, f, e},
          TriangleInput{c, g, h}, TriangleInput{c, h, d},
          TriangleInput{a, e, g}, TriangleInput{a, g, c},
          TriangleInput{b, d, h}, TriangleInput{b, h, f}};
}

std::array<TriangleInput, 4U> oblique_tetrahedron() {
  const Real3 a{-1.0, -1.0, -1.0};
  const Real3 b{1.0, -1.0, -1.0};
  const Real3 c{-1.0, 1.0, -1.0};
  const Real3 d{-1.0, -1.0, 1.0};
  return {{{a, c, b}, {a, b, d}, {a, d, c}, {b, c, d}}};
}

std::uint64_t sum_u64(std::uint64_t local) {
  std::uint64_t global{};
  MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

bool test_dual_interface_metric(const CartesianGeometryPlan& geometry,
                                const MeshPatch& patch) {
  const auto triangles = oblique_tetrahedron();
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  EBTopology topology;
  ImmersedPlanLimits limits;
  bool passed = expect(
      static_cast<bool>(StlScanCompiler::compile_triangles(
          geometry, patch, {triangles.data(), triangles.size()},
          CartesianAxis::y, kScanBudget, scan)) &&
          static_cast<bool>(ImmersedSurfaceCompiler::compile(scan, surface)) &&
          static_cast<bool>(EBTopologyCompiler::compile(
              MPI_COMM_WORLD, geometry, patch, scan, surface,
              ImmersedFluidSide::outside, limits, topology)),
      "oblique dual-metric fixture compiles");
  if (!passed) return false;

  const IbmInterfaceMetricPlan& metric = topology.interface_metric();
  const Span<const ImmersedLink> links = topology.links();
  const Span<const IbmInterfaceLinkMetric> physical = metric.links();
  passed &= expect(metric.fingerprint() != 0U && physical.size == links.size,
                   "every Cartesian link binds one physical metric");
  double local_control = 0.0;
  double local_area = 0.0;
  std::array<double, 3U> local_normal{};
  std::array<double, 3U> local_first_moment{};
  std::array<double, 9U> local_normal_moment{};
  std::array<double, 9U> local_normal_second_moment{};
  for (std::size_t index = 0U; index < links.size && index < physical.size;
       ++index) {
    passed &= expect(physical.data[index].global_link ==
                             links.data[index].global_link &&
                         physical.data[index].physical_quadrature_area > 0.0 &&
                         links.data[index].cartesian_control_face_area > 0.0,
                     "link metric preserves explicit control/physical identity");
    local_control += links.data[index].cartesian_control_face_area;
    local_area += physical.data[index].physical_quadrature_area;
    local_normal[0U] += physical.data[index].physical_area_vector.x;
    local_normal[1U] += physical.data[index].physical_area_vector.y;
    local_normal[2U] += physical.data[index].physical_area_vector.z;
    local_first_moment[0U] +=
        physical.data[index].physical_first_moment.x;
    local_first_moment[1U] +=
        physical.data[index].physical_first_moment.y;
    local_first_moment[2U] +=
        physical.data[index].physical_first_moment.z;
    for (std::size_t entry = 0U; entry < local_normal_moment.size(); ++entry) {
      local_normal_moment[entry] +=
          physical.data[index].normal_first_moment[entry];
      local_normal_second_moment[entry] +=
          physical.data[index].normal_second_moment[entry];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &local_control, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &local_area, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, local_normal.data(),
                static_cast<int>(local_normal.size()), MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, local_first_moment.data(),
                static_cast<int>(local_first_moment.size()), MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, local_normal_moment.data(),
                static_cast<int>(local_normal_moment.size()), MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, local_normal_second_moment.data(),
                static_cast<int>(local_normal_second_moment.size()),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double expected_area = 6.0 + 2.0 * std::sqrt(3.0);
  passed &= expect(std::abs(local_area - expected_area) <= 2.0e-12,
                   "physical link quadrature closes oblique STL area");
  passed &= expect(std::abs(local_control - local_area) > 0.5,
                   "Cartesian control area is not substituted for oblique STL area");
  passed &= expect(std::abs(local_normal[0U]) <= 2.0e-12 &&
                       std::abs(local_normal[1U]) <= 2.0e-12 &&
                       std::abs(local_normal[2U]) <= 2.0e-12,
                   "closed-surface uniform pressure has zero integrated normal");
  const double expected_first_moment =
      -(10.0 + 2.0 * std::sqrt(3.0)) / 3.0;
  for (const double value : local_first_moment)
    passed &= expect(std::abs(value - expected_first_moment) <= 3.0e-12,
                     "physical surface first moment closes analytically");
  constexpr double volume = 4.0 / 3.0;
  for (std::size_t row = 0U; row < 3U; ++row)
    for (std::size_t column = 0U; column < 3U; ++column)
      passed &= expect(
          std::abs(local_normal_moment[3U * row + column] -
                   (row == column ? volume : 0.0)) <= 3.0e-12,
          "first normal moment integrates linear traction exactly");
  const double oblique_normal_moment = 2.0 * std::sqrt(3.0) / 3.0;
  for (std::size_t row = 0U; row < 3U; ++row)
    for (std::size_t column = 0U; column < 3U; ++column)
      passed &= expect(
          std::abs(local_normal_second_moment[3U * row + column] -
                   (row == column ? 2.0 + oblique_normal_moment
                                  : oblique_normal_moment)) <= 3.0e-12,
          "physical normal second moment closes analytically");
  const IbmInterfaceMetricConservation& conservation =
      metric.conservation();
  passed &= expect(
      std::abs(conservation.cartesian_control_area - local_control) <=
              2.0e-12 &&
          std::abs(conservation.physical_quadrature_area - local_area) <=
              2.0e-12 &&
          std::abs(conservation.physical_area_vector.x - local_normal[0U]) <=
              2.0e-12 &&
          std::abs(conservation.physical_area_vector.y - local_normal[1U]) <=
              2.0e-12 &&
          std::abs(conservation.physical_area_vector.z - local_normal[2U]) <=
              2.0e-12 &&
          std::abs(conservation.physical_first_moment.x -
                   local_first_moment[0U]) <= 3.0e-12 &&
          std::abs(conservation.physical_first_moment.y -
                   local_first_moment[1U]) <= 3.0e-12 &&
          std::abs(conservation.physical_first_moment.z -
                   local_first_moment[2U]) <= 3.0e-12,
      "published conservation equals the link partition-of-unity sums");
  for (std::size_t entry = 0U; entry < local_normal_moment.size(); ++entry) {
    passed &= expect(
        std::abs(conservation.normal_first_moment[entry] -
                 local_normal_moment[entry]) <= 3.0e-12 &&
            std::abs(conservation.normal_second_moment[entry] -
                     local_normal_second_moment[entry]) <= 3.0e-12,
        "published tensor conservation equals the link partition sums");
  }
  const IbmInterfaceMetricResources resources = metric.resources();
  passed &= expect(resources.persistent_bytes_per_rank >=
                           physical.size * sizeof(IbmInterfaceLinkMetric) &&
                       resources.peak_bytes_per_rank >=
                           resources.persistent_bytes_per_rank &&
                       resources.collective_doubles_per_rank > 0U,
                   "metric plan publishes cold memory/collective budget");
  IbmInterfaceMetricPlan independently_compiled;
  passed &= expect(
      static_cast<bool>(IbmInterfaceMetricCompiler::compile(
          MPI_COMM_WORLD, geometry, patch, surface, topology, limits,
          independently_compiled)) &&
          independently_compiled.physical_fingerprint() ==
              metric.physical_fingerprint() &&
          independently_compiled.fingerprint() == metric.fingerprint(),
      "standalone metric compiler reproduces the sealed topology authority");
  const PlanFingerprint retained_metric =
      independently_compiled.fingerprint();
  const std::size_t retained_links = independently_compiled.links().size;
  ImmersedPlanLimits exhausted = limits;
  exhausted.maximum_persistent_bytes_per_rank = 1U;
  exhausted.maximum_peak_bytes_per_rank = 1U;
  const Status exhausted_status = IbmInterfaceMetricCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, surface, topology, exhausted,
      independently_compiled);
  passed &= expect(exhausted_status.code == StatusCode::invalid_plan &&
                       independently_compiled.fingerprint() ==
                           retained_metric &&
                       independently_compiled.links().size == retained_links,
                   "metric budget rejection is collective and atomic");
  return passed;
}

bool test_domain_clipped_surface_quadrature(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch) {
  const auto triangles = box({-0.9, -0.9, -2.0}, {-0.4, -0.4, 2.0});
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  EBTopology topology;
  BoundaryStencilPlan boundary;
  SurfaceQuadraturePlan quadrature;
  ImmersedPlanLimits limits;
  ImmersedDomainBoundaryPolicy policy;
  policy.allow_one_sided_quadratic[4U] = true;
  policy.allow_one_sided_quadratic[5U] = true;
  bool passed = expect(
      static_cast<bool>(StlScanCompiler::compile_triangles(
          geometry, patch, {triangles.data(), triangles.size()},
          CartesianAxis::y, kScanBudget, scan)) &&
          static_cast<bool>(ImmersedSurfaceCompiler::compile(scan, surface)) &&
          static_cast<bool>(EBTopologyCompiler::compile(
              MPI_COMM_WORLD, geometry, patch, scan, surface,
              ImmersedFluidSide::outside, limits, topology)) &&
          static_cast<bool>(BoundaryStencilCompiler::compile(
              MPI_COMM_WORLD, geometry, patch, surface, topology, policy,
              limits, boundary)),
      "EB/external-boundary intersection keeps full quadratic link plans");
  if (!passed) return false;
  const Span<const QuadraticStencilGroup> groups =
      boundary.reconstruction().groups();
  bool exercised_one_sided = false;
  for (std::size_t index = 0U; index < groups.size; ++index) {
    const QuadraticStencilQuality quality = groups.data[index].quality;
    exercised_one_sided |= quality.required_quadrant_mask != 0x0fU;
    passed &= expect(
        quality.required_quadrant_mask != 0U &&
            (quality.quadrant_mask & quality.required_quadrant_mask) ==
                quality.required_quadrant_mask &&
            quality.normal_band_count >= 3U && quality.rank == 10U &&
            quality.condition_estimate <= 1.0e8,
        "boundary-intersection stencil retains coverage/rank/condition gates");
  }
  const std::uint64_t one_sided_count =
      sum_u64(exercised_one_sided ? 1U : 0U);
  passed &= expect(one_sided_count > 0U,
                   "fixture exercises the distinct one-sided certificate");
  const RemoteDonorFieldSpec donor_field{1U, 1U};
  RemoteDonorExchangePlan sparse_donors;
  const Status sparse_status = RemoteDonorExchangePlan::analyze(
      MPI_COMM_WORLD, geometry.global_cells(), patch,
      boundary.reconstruction(), {&donor_field, 1U}, 171U, sparse_donors);
  passed &= expect(static_cast<bool>(sparse_status),
                   "ranks with no local IBM rows join sparse donor analysis");
  const Status quadrature_status = SurfaceQuadratureCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, surface, topology, policy, limits,
      quadrature);
  passed &= expect(static_cast<bool>(quadrature_status),
                   "domain-clipped surface quadrature compiles");
  double local_area = 0.0;
  const Span<const SurfaceQuadraturePoint> points = quadrature.local_points();
  for (std::size_t index = 0U; index < points.size; ++index) {
    const SurfaceQuadraturePoint& point = points.data[index];
    local_area += point.weight;
    passed &= expect(point.position.x >= -1.5 && point.position.x <= 1.5 &&
                         point.position.y >= -1.5 && point.position.y <= 1.5 &&
                         point.position.z >= -1.5 && point.position.z <= 1.5,
                     "no quadrature point lies outside the Cartesian domain");
  }
  double global_area = 0.0;
  MPI_Allreduce(&local_area, &global_area, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  passed &= expect(std::abs(global_area - 6.0) <= 2.0e-12,
                   "clipped extrusion integrates only in-domain side area");
  return passed;
}

struct PeriodicFixture {
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  EBTopology topology;
};

Status compile_periodic_fixture(const CartesianGeometryPlan& geometry,
                                const MeshPatch& patch, Real3 lower,
                                Real3 upper, PeriodicFixture& fixture) {
  const auto triangles = box(lower, upper);
  const Status scan_status = StlScanCompiler::compile_triangles(
      geometry, patch, {triangles.data(), triangles.size()}, CartesianAxis::y,
      kScanBudget, fixture.scan);
  if (!scan_status) return scan_status;
  const Status surface_status =
      ImmersedSurfaceCompiler::compile(fixture.scan, fixture.surface);
  if (!surface_status) return surface_status;
  ImmersedPlanLimits limits;
  return EBTopologyCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, fixture.scan, fixture.surface,
      ImmersedFluidSide::outside, limits, fixture.topology);
}

// The frozen periodic-IBM successor keeps the old rejection without explicit
// authority and exercises parameterized image/corner paths.
bool test_periodic_image_variant(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    Real3 lower, Real3 upper, const ImmersedDomainBoundaryPolicy& policy,
    bool check_no_periodic, bool require_lower_z_alias,
    bool require_double_axis_alias) {
  PeriodicFixture fixture;
  BoundaryStencilPlan boundary;
  ImmersedPlanLimits limits;
  const Status topology_status =
      compile_periodic_fixture(geometry, patch, lower, upper, fixture);
  const ImmersedSurfacePlan& surface = fixture.surface;
  const EBTopology& topology = fixture.topology;
  bool passed = true;
  if (check_no_periodic) {
    BoundaryStencilPlan no_periodic;
    ImmersedDomainBoundaryPolicy no_policy;
    const Status no_periodic_status =
        topology_status ? BoundaryStencilCompiler::compile(
                              MPI_COMM_WORLD, geometry, patch, surface,
                              topology, no_policy, limits, no_periodic)
                        : topology_status;
    passed &= expect(
        !no_periodic_status &&
            no_periodic_status.code == StatusCode::invalid_plan &&
            no_periodic_status.detail == 13406U,
        "same periodic-image fixture rejects without explicit policy (RED)");
  }
  const Status boundary_status =
      topology_status ? BoundaryStencilCompiler::compile(
                            MPI_COMM_WORLD, geometry, patch, surface, topology,
                            policy, limits, boundary)
                      : topology_status;
  passed &= expect(static_cast<bool>(boundary_status),
                   "periodic z image donor plan compiles (GREEN)");
  if (!boundary_status) return passed;
  passed &= expect(boundary.reconstruction().periodic_axis(CartesianAxis::z),
                   "periodic z authority is published on the reconstruction");
  const bool y_periodic = policy.allow_periodic_images[2U] &&
                          policy.allow_periodic_images[3U];
  if (y_periodic) {
    passed &= expect(
        boundary.reconstruction().periodic_axis(CartesianAxis::y),
        "periodic y authority is published on the reconstruction");
  }

  const Span<const Int3> donor_indices =
      boundary.reconstruction().donor_local_indices();
  const Span<const GlobalCellId> donor_globals =
      boundary.reconstruction().donor_global_cells();
  const std::int32_t gx = geometry.global_cells().x;
  const std::int32_t gy = geometry.global_cells().y;
  const std::int32_t gz = geometry.global_cells().z;
  const std::uint64_t plane = static_cast<std::uint64_t>(gx) * gy;
  std::uint64_t local_aliases = 0U;
  std::uint64_t local_lower_z_aliases = 0U;
  std::uint64_t local_upper_z_aliases = 0U;
  std::uint64_t local_double_aliases = 0U;
  const Span<const QuadraticStencilGroup> groups =
      boundary.reconstruction().groups();
  passed &= expect(donor_indices.size == donor_globals.size,
                   "periodic donor index/global spans have equal size");
  for (std::size_t group_index = 0U; group_index < groups.size;
       ++group_index) {
    const QuadraticStencilGroup& group = groups.data[group_index];
    const std::size_t begin = group.donor_begin;
    const std::size_t count = group.quality.donor_count;
    const bool in_range = begin <= donor_globals.size &&
                          count <= donor_globals.size - begin &&
                          count <= donor_indices.size -
                                       std::min(begin, donor_indices.size);
    passed &= expect(in_range,
                     "periodic stencil donor range stays within published spans");
    if (!in_range) continue;
    for (std::size_t first = 0U; first < count; ++first) {
      for (std::size_t second = 0U; second < first; ++second) {
        passed &= expect(
            donor_globals.data[begin + first] !=
                donor_globals.data[begin + second],
            "canonical donor IDs are unique within each periodic stencil");
      }
    }
  }
  for (std::size_t index = 0U; index < donor_indices.size; ++index) {
    const Int3 local = donor_indices.data[index];
    const Int3 raw{patch.begin.x + local.x, patch.begin.y + local.y,
                   patch.begin.z + local.z};
    const bool y_alias = raw.y < 0 || raw.y >= gy;
    const bool z_alias = raw.z < 0 || raw.z >= gz;
    const bool alias = y_alias || z_alias;
    if (alias) {
      ++local_aliases;
      const GlobalCellId id = donor_globals.data[index];
      const std::uint64_t remainder = id % plane;
      const std::int32_t canonical_y =
          static_cast<std::int32_t>(remainder / gx);
      const std::int32_t canonical_z =
          static_cast<std::int32_t>(id / plane);
      if (y_alias) {
        const std::int32_t expected_y = raw.y < 0 ? raw.y + gy : raw.y - gy;
        passed &= expect(
            canonical_y == expected_y,
            "periodic donor stores canonical y ID beside raw alias");
      }
      if (z_alias) {
        if (raw.z < 0)
          ++local_lower_z_aliases;
        else
          ++local_upper_z_aliases;
        const std::int32_t expected_z = raw.z < 0 ? raw.z + gz : raw.z - gz;
        passed &= expect(
            canonical_z == expected_z,
            "periodic donor stores canonical z ID beside raw alias");
      }
      if (y_alias && z_alias) ++local_double_aliases;
    }
  }
  const std::uint64_t global_aliases = sum_u64(local_aliases);
  passed &= expect(global_aliases > 0U,
                   "periodic fixture exercises an out-of-domain donor alias");
  const std::uint64_t global_lower_z_aliases =
      sum_u64(local_lower_z_aliases);
  const std::uint64_t global_upper_z_aliases =
      sum_u64(local_upper_z_aliases);
  if (require_lower_z_alias) {
    passed &= expect(global_lower_z_aliases > 0U,
                     "z-min fixture exercises lower periodic aliases");
  } else {
    passed &= expect(global_upper_z_aliases > 0U,
                     "z-max fixture exercises upper periodic aliases");
  }
  if (require_double_axis_alias) {
    passed &= expect(sum_u64(local_double_aliases) > 0U,
                     "corner fixture exercises simultaneous y/z aliases");
  }

  constexpr std::int32_t ghost = 4;
  const std::size_t nx = static_cast<std::size_t>(patch.cells.x + 2 * ghost);
  const std::size_t ny = static_cast<std::size_t>(patch.cells.y + 2 * ghost);
  const std::size_t nz = static_cast<std::size_t>(patch.cells.z + 2 * ghost);
  const std::size_t component_stride = nx * ny * nz;
  std::vector<double> storage(component_stride, -777.0);
  FieldView donor_field;
  donor_field.base = storage.data() + ghost + nx * ghost + nx * ny * ghost;
  donor_field.interior = patch.cells;
  donor_field.ghosts = {ghost, ghost, ghost};
  donor_field.components = 1U;
  donor_field.stride_y = nx;
  donor_field.stride_z = nx * ny;
  donor_field.component_stride = component_stride;
  donor_field.field = 778U;
  donor_field.revision = 20U;
  donor_field.storage_identity = 7781U;
  donor_field.revision_domain = 7782U;
  for (std::int32_t z = 0; z < patch.cells.z; ++z)
    for (std::int32_t y = 0; y < patch.cells.y; ++y)
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 local{x, y, z};
        const Int3 global{patch.begin.x + x, patch.begin.y + y,
                          patch.begin.z + z};
        donor_field.unchecked(local, 0U) =
            1.0 + global.x + 10.0 * global.y + 100.0 * global.z;
      }
  const RemoteDonorFieldSpec donor_field_spec{donor_field.field, 1U};
  RemoteDonorExchangePlan donor_exchange;
  const Status exchange_compile = RemoteDonorExchangePlan::compile(
      MPI_COMM_WORLD, geometry.global_cells(), patch,
      boundary.reconstruction(), {&donor_field_spec, 1U}, 178U,
      donor_exchange);
  passed &= expect(static_cast<bool>(exchange_compile) && donor_exchange.ready(),
                   "periodic donor exchange binds self/remote aliases");
  std::array<FieldView, 1U> exchange_fields{{donor_field}};
  if (exchange_compile) {
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::array<FieldView, 1U> invalid_fields = exchange_fields;
    if (rank == size - 1) ++invalid_fields[0U].field;
    const RemoteDonorExchangeCounters counters_before =
        donor_exchange.runtime_counters();
    Status preflight = donor_exchange.preflight_exchange(
        178U, {invalid_fields.data(), invalid_fields.size()});
    const int local_preflight_ok = preflight ? 1 : 0;
    int global_preflight_ok = 0;
    MPI_Allreduce(&local_preflight_ok, &global_preflight_ok, 1, MPI_INT,
                  MPI_MIN, MPI_COMM_WORLD);
    const RemoteDonorExchangeCounters counters_after_preflight =
        donor_exchange.runtime_counters();
    passed &= expect(
        global_preflight_ok == 0 &&
            (rank != size - 1 ||
             preflight.code == StatusCode::invalid_plan) &&
            counters_after_preflight.exchange_calls ==
                counters_before.exchange_calls &&
            counters_after_preflight.peer_messages ==
                counters_before.peer_messages &&
            counters_after_preflight.bytes == counters_before.bytes,
        "one-rank invalid donor view is rejected by a no-write preflight before persistent communication");
    const Status collective_rejection = donor_exchange.exchange(
        178U, {invalid_fields.data(), invalid_fields.size()});
    const std::array<std::uint64_t, 2U> local_rejection{{
        static_cast<std::uint64_t>(collective_rejection.code),
        collective_rejection.detail}};
    std::array<std::uint64_t, 2U> minimum_rejection{};
    std::array<std::uint64_t, 2U> maximum_rejection{};
    MPI_Allreduce(local_rejection.data(), minimum_rejection.data(), 2,
                  MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(local_rejection.data(), maximum_rejection.data(), 2,
                  MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    const RemoteDonorExchangeCounters counters_after_rejection =
        donor_exchange.runtime_counters();
    passed &= expect(
        minimum_rejection == maximum_rejection &&
            collective_rejection.code == StatusCode::invalid_plan &&
            counters_after_rejection.exchange_calls ==
                counters_before.exchange_calls &&
            counters_after_rejection.peer_messages ==
                counters_before.peer_messages &&
            counters_after_rejection.bytes == counters_before.bytes,
        "one-rank invalid donor view is globalized inside exchange before persistent communication");
    preflight = donor_exchange.preflight_exchange(
        178U, {exchange_fields.data(), exchange_fields.size()});
    passed &= expect(static_cast<bool>(donor_exchange.exchange(
                             178U, {exchange_fields.data(),
                                    exchange_fields.size()})) &&
                         static_cast<bool>(preflight),
                     "periodic donor exchange completes through persistent plan");
    for (std::size_t index = 0U; index < donor_indices.size; ++index) {
      const Int3 local = donor_indices.data[index];
      const bool ghost_target = local.x < 0 || local.y < 0 || local.z < 0 ||
                                local.x >= patch.cells.x ||
                                local.y >= patch.cells.y ||
                                local.z >= patch.cells.z;
      if (!ghost_target) continue;
      const GlobalCellId id = donor_globals.data[index];
      const std::uint64_t remainder = id % plane;
      const double expected =
          1.0 + static_cast<double>(remainder % gx) +
          10.0 * static_cast<double>(remainder / gx) +
          100.0 * static_cast<double>(id / plane);
      passed &= expect(exchange_fields[0U].unchecked(local, 0U) == expected,
                       "periodic raw ghost target receives canonical self/remote value");
    }
  }
  return passed;
}

bool test_periodic_image_donor_red(const CartesianGeometryPlan& geometry,
                                   const MeshPatch& patch) {
  ImmersedDomainBoundaryPolicy policy;
  policy.allow_periodic_images[4U] = true;
  policy.allow_periodic_images[5U] = true;
  return test_periodic_image_variant(
      geometry, patch, {-0.4, -0.4, 1.25}, {0.4, 0.4, 1.75}, policy, true,
      false, false);
}

bool test_periodic_image_donor_zmin(const CartesianGeometryPlan& geometry,
                                    const MeshPatch& patch) {
  ImmersedDomainBoundaryPolicy policy;
  policy.allow_periodic_images[4U] = true;
  policy.allow_periodic_images[5U] = true;
  return test_periodic_image_variant(
      geometry, patch, {-0.4, -0.4, -1.75}, {0.4, 0.4, -1.25}, policy, false,
      true, false);
}

bool test_periodic_image_donor_corner(const CartesianGeometryPlan& geometry,
                                      const MeshPatch& patch) {
  ImmersedDomainBoundaryPolicy policy;
  policy.allow_periodic_images[2U] = true;
  policy.allow_periodic_images[3U] = true;
  policy.allow_periodic_images[4U] = true;
  policy.allow_periodic_images[5U] = true;
  return test_periodic_image_variant(
      geometry, patch, {-0.4, 1.25, 1.25}, {0.4, 1.75, 1.75}, policy, false,
      false, true);
}

bool test_periodic_policy_rejections() {
  CartesianMeshSpec tiny_mesh = mesh_spec();
  tiny_mesh.exact_cells = {24, 24, 8};
  CartesianGeometryPlan tiny_geometry;
  MeshPatch tiny_patch;
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, tiny_mesh, GeometryBudget{}, tiny_geometry,
          tiny_patch)),
      "small periodic extent geometry compiles");
  PeriodicFixture fixture;
  const Status fixture_status = compile_periodic_fixture(
      tiny_geometry, tiny_patch, {-0.4, -0.4, 1.25}, {0.4, 0.4, 1.75},
      fixture);
  passed &= expect(static_cast<bool>(fixture_status),
                   "small periodic extent fixture compiles");
  if (!fixture_status) return passed;

  ImmersedPlanLimits limits;
  BoundaryStencilPlan rejected;
  ImmersedDomainBoundaryPolicy half_pair;
  half_pair.allow_periodic_images[4U] = true;
  const Status half_pair_status = BoundaryStencilCompiler::compile(
      MPI_COMM_WORLD, tiny_geometry, tiny_patch, fixture.surface,
      fixture.topology, half_pair, limits, rejected);
  passed &= expect(half_pair_status.code == StatusCode::invalid_plan &&
                       half_pair_status.detail == 13401U,
                   "half-periodic pair fails closed");

  ImmersedDomainBoundaryPolicy mixed_face;
  mixed_face.allow_periodic_images[4U] = true;
  mixed_face.allow_periodic_images[5U] = true;
  mixed_face.allow_one_sided_quadratic[4U] = true;
  const Status mixed_face_status = BoundaryStencilCompiler::compile(
      MPI_COMM_WORLD, tiny_geometry, tiny_patch, fixture.surface,
      fixture.topology, mixed_face, limits, rejected);
  passed &= expect(mixed_face_status.code == StatusCode::invalid_plan &&
                       mixed_face_status.detail == 13401U,
                   "periodic and one-sided authority on one face fails closed");

  ImmersedDomainBoundaryPolicy short_extent;
  short_extent.allow_periodic_images[4U] = true;
  short_extent.allow_periodic_images[5U] = true;
  const Status short_extent_status = BoundaryStencilCompiler::compile(
      MPI_COMM_WORLD, tiny_geometry, tiny_patch, fixture.surface,
      fixture.topology, short_extent, limits, rejected);
  passed &= expect(short_extent_status.code == StatusCode::invalid_plan &&
                       short_extent_status.detail == 13401U,
                   "periodic extent at twice reach fails closed");
  return passed;
}

bool run() {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_WORLD, mesh_spec(), GeometryBudget{},
                           geometry, patch)),
                       "geometry compiles");
  const auto triangles = cube();
  StlScanPlan scan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {triangles.data(), triangles.size()},
                       CartesianAxis::y, kScanBudget, scan)),
                   "scan compiles");
  ImmersedSurfacePlan surface;
  passed &= expect(static_cast<bool>(
                       ImmersedSurfaceCompiler::compile(scan, surface)),
                   "canonical surface compiles");
  const Span<const std::array<std::uint32_t, 3U>> neighbours =
      surface.triangle_neighbours();
  bool reciprocal_neighbours = neighbours.size == surface.triangles().size;
  for (std::size_t triangle = 0U;
       reciprocal_neighbours && triangle < neighbours.size; ++triangle) {
    for (const std::uint32_t adjacent : neighbours.data[triangle]) {
      reciprocal_neighbours &= adjacent < neighbours.size;
      if (adjacent < neighbours.size) {
        const auto other = neighbours.data[adjacent];
        reciprocal_neighbours &=
            std::find(other.begin(), other.end(), triangle) != other.end();
      }
    }
  }
  passed &= expect(reciprocal_neighbours,
                   "canonical surface persists reciprocal triangle adjacency");

  ImmersedPlanLimits limits;
  EBTopology topology;
  passed &= expect(static_cast<bool>(EBTopologyCompiler::compile(
                       MPI_COMM_WORLD, geometry, patch, scan, surface,
                       ImmersedFluidSide::outside, limits, topology)),
                   "outside-fluid topology compiles");
  passed &= expect(!topology.links().size || topology.fingerprint() != 0U,
                   "topology publishes nonzero identity");

  BoundaryStencilPlan boundary;
  const Status boundary_status = BoundaryStencilCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, surface, topology, limits, boundary);
  if (!boundary_status) {
    std::cerr << "boundary status code="
              << static_cast<unsigned>(boundary_status.code)
              << " detail=" << boundary_status.detail << '\n';
  }
  passed &= expect(static_cast<bool>(boundary_status),
                   "strict quadratic boundary plan compiles");
  passed &= expect(boundary.links().size == topology.links().size,
                   "each local topology link owns one boundary stencil");
  passed &= expect(boundary.maximum_halo_reach() <= 4U &&
                       boundary.fingerprint() != 0U &&
                       boundary.reconstruction().fingerprint() != 0U,
                   "boundary plan respects halo and identity contracts");
  const Span<const QuadraticStencilGroup> groups =
      boundary.reconstruction().groups();
  for (std::size_t index = 0U; index < groups.size; ++index) {
    const QuadraticStencilGroup& group = groups.data[index];
    passed &= expect(group.quality.rank == 10U &&
                         group.quality.donor_count >= 14U &&
                         group.quality.donor_count <= 32U &&
                         group.quality.normal_band_count >= 3U &&
                         group.quality.quadrant_mask == 0x0fU &&
                         group.quality.required_quadrant_mask == 0x0fU &&
                         group.quality.reach <= 4U &&
                         group.quality.condition_estimate <= 1.0e8,
                     "every boundary group satisfies the strict contract");
  }
  std::uint64_t local_positive_ghost_donors = 0U;
  const Span<const BoundaryStencilLink> boundary_links = boundary.links();
  const Span<const ImmersedLink> topology_links = topology.links();
  const Span<const Int3> all_donor_indices =
      boundary.reconstruction().donor_local_indices();
  for (std::size_t index = 0U; index < boundary_links.size; ++index) {
    const BoundaryStencilLink& binding = boundary_links.data[index];
    const ImmersedLink& link = topology_links.data[binding.topology_link];
    const QuadraticStencilGroup& group =
        groups.data[binding.reconstruction_group];
    for (std::size_t donor = 0U; donor < group.quality.donor_count; ++donor) {
      const Int3 local =
          all_donor_indices.data[group.donor_begin + donor];
      const Int3 global{patch.begin.x + local.x, patch.begin.y + local.y,
                        patch.begin.z + local.z};
      const Real3 centre{
          geometry.x().centres().data[static_cast<std::size_t>(global.x)],
          geometry.y().centres().data[static_cast<std::size_t>(global.y)],
          geometry.z().centres().data[static_cast<std::size_t>(global.z)]};
      const Real3 offset{centre.x - link.wall_point.x,
                         centre.y - link.wall_point.y,
                         centre.z - link.wall_point.z};
      const double normal_coordinate =
          offset.x * link.solid_to_fluid_normal.x +
          offset.y * link.solid_to_fluid_normal.y +
          offset.z * link.solid_to_fluid_normal.z;
      passed &= expect(normal_coordinate > 0.0,
                       "every production boundary donor is positive-normal");
      const bool ghost = local.x < 0 || local.y < 0 || local.z < 0 ||
                         local.x >= patch.cells.x ||
                         local.y >= patch.cells.y ||
                         local.z >= patch.cells.z;
      if (ghost && normal_coordinate > 0.0) {
        ++local_positive_ghost_donors;
      }
    }
  }
  const std::uint64_t global_positive_ghost_donors =
      sum_u64(local_positive_ghost_donors);
  int stencil_communicator_size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &stencil_communicator_size);
  if (stencil_communicator_size > 1) {
    passed &= expect(
        global_positive_ghost_donors > 0U,
        "distributed production fixture exercises positive-normal ghost donors");
  }
  {
    constexpr std::int32_t ghost = 4;
    const std::size_t nx = static_cast<std::size_t>(patch.cells.x + 2 * ghost);
    const std::size_t ny = static_cast<std::size_t>(patch.cells.y + 2 * ghost);
    const std::size_t nz = static_cast<std::size_t>(patch.cells.z + 2 * ghost);
    const std::size_t component_stride = nx * ny * nz;
    std::vector<double> storage(2U * component_stride, -999.0);
    FieldView donor_field;
    donor_field.base = storage.data() + ghost + nx * ghost + nx * ny * ghost;
    donor_field.interior = patch.cells;
    donor_field.ghosts = {ghost, ghost, ghost};
    donor_field.components = 2U;
    donor_field.stride_y = nx;
    donor_field.stride_z = nx * ny;
    donor_field.component_stride = component_stride;
    donor_field.field = 777U;
    donor_field.revision = 19U;
    donor_field.storage_identity = 7771U;
    donor_field.revision_domain = 7772U;
    for (std::int32_t z = 0; z < patch.cells.z; ++z)
      for (std::int32_t y = 0; y < patch.cells.y; ++y)
        for (std::int32_t x = 0; x < patch.cells.x; ++x) {
          const Int3 local{x, y, z};
          const Int3 global{patch.begin.x + x, patch.begin.y + y,
                            patch.begin.z + z};
          const double value = 1.0 + global.x + 10.0 * global.y +
                               100.0 * global.z;
          donor_field.unchecked(local, 0U) = value;
          donor_field.unchecked(local, 1U) = -value;
        }
    const std::array<RemoteDonorFieldSpec, 1U> donor_fields{{
        {donor_field.field, donor_field.components}}};
    RemoteDonorExchangePlan donor_exchange;
    passed &= expect(
        static_cast<bool>(RemoteDonorExchangePlan::analyze(
            MPI_COMM_WORLD, geometry.global_cells(), patch,
            boundary.reconstruction(),
            {donor_fields.data(), donor_fields.size()}, 177U,
            donor_exchange)) &&
            !donor_exchange.ready() && donor_exchange.fingerprint() != 0U,
        "compact remote-donor metadata is analyzed before buffer binding");
    std::array<FieldView, 1U> exchange_fields{{donor_field}};
    passed &= expect(
        donor_exchange.exchange(
            177U, {exchange_fields.data(), exchange_fields.size()})
                .code == StatusCode::invalid_plan &&
            static_cast<bool>(donor_exchange.bind(MPI_COMM_WORLD)) &&
            donor_exchange.ready(),
        "remote-donor buffers and persistent requests bind after analysis");
    const Status wrong_stage = donor_exchange.exchange(
        176U, {exchange_fields.data(), exchange_fields.size()});
    passed &= expect(wrong_stage.code == StatusCode::invalid_plan,
                     "remote-donor exchange rejects the wrong hot stage");
    passed &= expect(
        static_cast<bool>(donor_exchange.exchange(
            177U, {exchange_fields.data(), exchange_fields.size()})),
        "compact remote-donor exchange completes without a grown-box gather");
    const auto donor_globals =
        boundary.reconstruction().donor_global_cells();
    for (std::size_t index = 0U; index < all_donor_indices.size; ++index) {
      const Int3 local = all_donor_indices.data[index];
      const GlobalCellId id = donor_globals.data[index];
      const std::uint64_t gx =
          static_cast<std::uint64_t>(geometry.global_cells().x);
      const std::uint64_t gy =
          static_cast<std::uint64_t>(geometry.global_cells().y);
      const std::uint64_t plane = gx * gy;
      const double global_z = static_cast<double>(id / plane);
      const std::uint64_t remainder = id % plane;
      const double global_y = static_cast<double>(remainder / gx);
      const double global_x = static_cast<double>(remainder % gx);
      const double expected =
          1.0 + global_x + 10.0 * global_y + 100.0 * global_z;
      passed &= expect(
          exchange_fields[0U].unchecked(local, 0U) == expected &&
              exchange_fields[0U].unchecked(local, 1U) == -expected,
          "every local, face, edge, and corner donor equals the global oracle");
    }
    const RemoteDonorExchangeStats donor_stats = donor_exchange.stats();
    passed &= expect(
        donor_stats.bytes_per_exchange ==
                2U * sizeof(double) *
                    (donor_stats.received_cells + donor_stats.supplied_cells) &&
            (stencil_communicator_size == 1 ||
             donor_stats.peer_messages != 0U),
        "compact donor resource counters report only requested values");
  }
  if (boundary.reconstruction().rows().size != 0U &&
      boundary.reconstruction().donor_local_indices().size != 0U) {
    const Span<const Int3> donor_indices =
        boundary.reconstruction().donor_local_indices();
    const QuadraticStencilGroup& group = groups.data[0U];
    std::vector<double> storage(24U * 24U * 24U, 1.0);
    ConstFieldView field;
    field.base = storage.data();
    field.interior = patch.cells;
    field.ghosts = {0, 0, 0};
    field.components = 1U;
    field.stride_y = 24U;
    field.stride_z = 24U * 24U;
    field.component_stride = storage.size();
    field.field = 1U;
    field.revision = 1U;
    field.storage_identity = 1U;
    field.revision_domain = 1U;
    bool all_owned = true;
    for (std::size_t donor = 0U; donor < group.quality.donor_count; ++donor) {
      const Int3 index = donor_indices.data[group.donor_begin + donor];
      all_owned = all_owned && index.x >= 0 && index.y >= 0 && index.z >= 0 &&
                  index.x < patch.cells.x && index.y < patch.cells.y &&
                  index.z < patch.cells.z;
    }
    if (all_owned) {
      double evaluated = -1.0;
      passed &= expect(static_cast<bool>(evaluate_quadratic_row(
                           boundary.reconstruction(), group.row_begin, field,
                           0U, 1.0, 0.0, evaluated)) &&
                           std::isfinite(evaluated),
                       "compiled boundary reconstruction is executable");
    }
  }

  SurfaceQuadraturePlan quadrature;
  const Status quadrature_status = SurfaceQuadratureCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, surface, topology, limits, quadrature);
  if (!quadrature_status) {
    std::cerr << "quadrature status code="
              << static_cast<unsigned>(quadrature_status.code)
              << " detail=" << quadrature_status.detail << '\n';
  }
  passed &= expect(static_cast<bool>(quadrature_status),
                   "surface quadrature plan compiles");
  passed &= expect(quadrature.reconstruction().fingerprint() != 0U,
                   "quadrature reconstruction publishes executable identity");
  passed &= expect(quadrature.global_point_count() ==
                       3U * surface.triangles().size &&
                       sum_u64(quadrature.local_points().size) ==
                           quadrature.global_point_count(),
                   "three-point quadrature has unique global ownership");
  double local_weight = 0.0;
  bool exercised_pressure_row = false;
  const Span<const SurfaceQuadraturePoint> local_points =
      quadrature.local_points();
  for (std::size_t index = 0U; index < local_points.size; ++index) {
    const SurfaceQuadraturePoint& point = local_points.data[index];
    passed &= expect(point.owner_rank >= 0 && point.weight > 0.0 &&
                         point.reconstruction_group != kInvalidIbmIndex,
                     "local quadrature point has owner and reconstruction");
    if (!exercised_pressure_row) {
      const QuadraticStencilPlan& reconstruction = quadrature.reconstruction();
      const QuadraticStencilGroup& group =
          reconstruction.groups().data[point.reconstruction_group];
      const Span<const Int3> donor_indices =
          reconstruction.donor_local_indices();
      bool addressable = true;
      for (std::size_t donor = 0U; donor < group.quality.donor_count; ++donor) {
        const Int3 index = donor_indices.data[group.donor_begin + donor];
        addressable = addressable && index.x >= -4 && index.y >= -4 &&
                      index.z >= -4 && index.x < patch.cells.x + 4 &&
                      index.y < patch.cells.y + 4 &&
                      index.z < patch.cells.z + 4;
      }
      if (addressable) {
        constexpr std::int32_t ghost = 4;
        const std::size_t nx = static_cast<std::size_t>(patch.cells.x + 8);
        const std::size_t ny = static_cast<std::size_t>(patch.cells.y + 8);
        const std::size_t nz = static_cast<std::size_t>(patch.cells.z + 8);
        std::vector<double> pressure_storage(nx * ny * nz, 2.75);
        ConstFieldView pressure;
        pressure.base = pressure_storage.data() + ghost + nx * ghost +
                        nx * ny * ghost;
        pressure.interior = patch.cells;
        pressure.ghosts = {ghost, ghost, ghost};
        pressure.components = 1U;
        pressure.stride_y = nx;
        pressure.stride_z = nx * ny;
        pressure.component_stride = pressure_storage.size();
        pressure.field = 2U;
        pressure.revision = 1U;
        pressure.storage_identity = 2U;
        pressure.revision_domain = 1U;
        double first = 0.0;
        double mutated = 0.0;
        passed &= expect(
            static_cast<bool>(evaluate_quadratic_row(
                reconstruction, point.wall_value_row, pressure, 0U, -101.0,
                9.0, first)) &&
                static_cast<bool>(evaluate_quadratic_row(
                    reconstruction, point.wall_value_row, pressure, 0U,
                    407.0, -13.0, mutated)) &&
                std::abs(first - 2.75) <= 1.0e-11 &&
                std::abs(mutated - first) <= 1.0e-13,
            "surface pressure row reconstructs donors and ignores wall inputs");
        exercised_pressure_row = true;
      }
    }
    local_weight += point.weight;
  }
  double global_weight = 0.0;
  MPI_Allreduce(&local_weight, &global_weight, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  passed &= expect(std::abs(global_weight - 6.0) <= 1.0e-12,
                   "quadrature weights close cube area");
  passed &= expect(
      std::abs(topology.interface_metric()
                   .conservation()
                   .physical_quadrature_area -
               global_weight) <= 1.0e-12,
      "link physical metric and surface quadrature share one STL area");
  const std::uint64_t global_pressure_rows =
      sum_u64(exercised_pressure_row ? 1U : 0U);
  passed &= expect(global_pressure_rows > 0U,
                   "at least one rank executes the surface pressure mutation RED");

  std::uint64_t maximum_local_points = local_points.size;
  std::uint64_t minimum_local_points = local_points.size;
  MPI_Allreduce(MPI_IN_PLACE, &maximum_local_points, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &minimum_local_points, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  int communicator_size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &communicator_size);
  passed &= expect(maximum_local_points <= quadrature.global_point_count() &&
                       minimum_local_points <= maximum_local_points,
                   "quadrature reports a bounded deterministic local distribution");
  if (communicator_size > 1) {
    passed &= expect(maximum_local_points < quadrature.global_point_count(),
                     "distributed ranks own less than the replicated global plan");
  }

  ImmersedPlanLimits tight_quadrature_limits = limits;
  tight_quadrature_limits.maximum_local_quadrature_points =
      maximum_local_points;
  SurfaceQuadraturePlan tight_quadrature;
  passed &= expect(
      static_cast<bool>(SurfaceQuadratureCompiler::compile(
          MPI_COMM_WORLD, geometry, patch, surface, topology,
          tight_quadrature_limits, tight_quadrature)) &&
          tight_quadrature.local_points().size == local_points.size &&
          tight_quadrature.global_point_count() ==
              quadrature.global_point_count(),
      "exact maximum-local budget compiles without reserving the global plan");

  if (maximum_local_points > 0U) {
    const PlanFingerprint retained_layout =
        tight_quadrature.local_layout_fingerprint();
    const PlanFingerprint retained_physical =
        tight_quadrature.physical_fingerprint();
    const std::size_t retained_local_count =
        tight_quadrature.local_points().size;
    ImmersedPlanLimits insufficient_quadrature_limits =
        tight_quadrature_limits;
    insufficient_quadrature_limits.maximum_local_quadrature_points =
        maximum_local_points - 1U;
    const Status insufficient = SurfaceQuadratureCompiler::compile(
        MPI_COMM_WORLD, geometry, patch, surface, topology,
        insufficient_quadrature_limits, tight_quadrature);
    passed &= expect(
        insufficient.code == StatusCode::invalid_plan &&
            tight_quadrature.local_layout_fingerprint() == retained_layout &&
            tight_quadrature.physical_fingerprint() == retained_physical &&
            tight_quadrature.local_points().size == retained_local_count,
        "under-local quadrature budget rejects collectively and atomically");
  }

  const PlanFingerprint retained = boundary.fingerprint();
  ImmersedPlanLimits rejected_limits = limits;
  rejected_limits.stencil.maximum_donors = 13U;
  const Status rejected = BoundaryStencilCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, surface, topology, rejected_limits,
      boundary);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       boundary.fingerprint() == retained,
                   "invalid donor policy rejects without publishing");
  passed &= test_dual_interface_metric(geometry, patch);
  passed &= test_domain_clipped_surface_quadrature(geometry, patch);
  passed &= test_periodic_image_donor_red(geometry, patch);
  passed &= test_periodic_image_donor_zmin(geometry, patch);
  passed &= test_periodic_image_donor_corner(geometry, patch);
  passed &= test_periodic_policy_rejections();
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool local = run();
  const int value = local ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return global == 1 ? 0 : 1;
}
