// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr Int3 kGlobal{16, 16, 16};
constexpr PlanFingerprint kMarkerSource = UINT64_C(0x6a7f31d94bc285e1);
constexpr std::size_t kCellCount = 16U * 16U * 16U;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

std::uint64_t sum(std::uint64_t local) {
  std::uint64_t global = 0U;
  MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

std::size_t flat(Int3 shape, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(cell.z));
}

GlobalCellId global_id(Int3 cell) noexcept {
  return static_cast<GlobalCellId>(flat(kGlobal, cell));
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-0.008, -0.008, -0.008};
  mesh.upper = {0.008, 0.008, 0.008};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kGlobal;
  mesh.minimum_spacing = {1.0e-12, 1.0e-12, 1.0e-12};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(1000000), UINT64_C(1073741824)};
  return mesh;
}

std::vector<std::uint8_t> plane_marker() {
  std::vector<std::uint8_t> marker(kCellCount, 0U);
  for (std::int32_t z = 0; z < kGlobal.z; ++z)
    for (std::int32_t y = 0; y < kGlobal.y; ++y)
      for (std::int32_t x = 8; x < kGlobal.x; ++x)
        marker[flat(kGlobal, {x, y, z})] = 1U;
  return marker;
}

std::vector<std::uint8_t> corner_marker() {
  std::vector<std::uint8_t> marker(kCellCount, 0U);
  for (std::int32_t z = 8; z < kGlobal.z; ++z)
    for (std::int32_t y = 8; y < kGlobal.y; ++y)
      for (std::int32_t x = 8; x < kGlobal.x; ++x)
        marker[flat(kGlobal, {x, y, z})] = 1U;
  return marker;
}

std::vector<std::uint8_t> thin_sheet_marker() {
  std::vector<std::uint8_t> marker(kCellCount, 0U);
  // An internal one-cell-thick fluid sheet has many tangential donors but
  // exactly one positive-normal z band on both faces.  It is the minimal
  // deterministic form of the GTMC 1304 defect; a six-neighbour aggregate
  // normal would cancel here and therefore cannot define the wall geometry.
  for (std::int32_t y = 0; y < kGlobal.y; ++y)
    for (std::int32_t x = 0; x < kGlobal.x; ++x)
      marker[flat(kGlobal, {x, y, 8})] = 1U;
  return marker;
}

std::vector<std::uint8_t> selection_edge_marker() {
  std::vector<std::uint8_t> marker(kCellCount, 0U);
  // At the x-negative wall along y=8, more than 32 nearer fluid donors lie in
  // the y=8 plane, while equally valid reach-four donors also exist at y<8.
  // Candidate truncation must retain that available tangential diversity.
  for (std::int32_t z = 0; z < kGlobal.z; ++z)
    for (std::int32_t y = 0; y <= 8; ++y)
      for (std::int32_t x = 8; x < kGlobal.x; ++x)
        marker[flat(kGlobal, {x, y, z})] = 1U;
  return marker;
}

Real3 expected_normal(ImmersedFaceDirection direction) noexcept {
  switch (direction) {
    case ImmersedFaceDirection::x_negative:
      return {1.0, 0.0, 0.0};
    case ImmersedFaceDirection::x_positive:
      return {-1.0, 0.0, 0.0};
    case ImmersedFaceDirection::y_negative:
      return {0.0, 1.0, 0.0};
    case ImmersedFaceDirection::y_positive:
      return {0.0, -1.0, 0.0};
    case ImmersedFaceDirection::z_negative:
      return {0.0, 0.0, 1.0};
    case ImmersedFaceDirection::z_positive:
      return {0.0, 0.0, -1.0};
  }
  return {};
}

bool same(Real3 left, Real3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool compile_and_check(const CartesianGeometryPlan& geometry,
                       const MeshPatch& patch,
                       const std::vector<std::uint8_t>& marker,
                       std::uint64_t expected_links, bool corner, int rank,
                       EBTopology& topology, BoundaryStencilPlan& boundary,
                       SurfaceQuadraturePlan& quadrature) {
  ImmersedPlanLimits limits;
  const Status status = ImportedIbmCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, {marker.data(), marker.size()},
      kMarkerSource, limits, topology, boundary, quadrature);
  bool passed = expect(static_cast<bool>(status), rank,
                       corner ? "corner marker compiles"
                              : "plane marker compiles");
  if (!passed) {
    std::cerr << "rank " << rank << " status="
              << static_cast<unsigned>(status.code)
              << " detail=" << status.detail << '\n';
    return false;
  }
  const Span<const ImmersedLink> links = topology.links();
  const Span<const IbmInterfaceLinkMetric> metrics =
      topology.interface_metric().links();
  const Span<const BoundaryStencilLink> boundary_links = boundary.links();
  const Span<const SurfaceQuadraturePoint> points = quadrature.local_points();
  passed &= expect(sum(links.size) == expected_links, rank,
                   "imported marker has exact global link count");
  passed &= expect(metrics.size == links.size &&
                       boundary_links.size == links.size &&
                       points.size == links.size &&
                       quadrature.global_point_count() == expected_links,
                   rank, "topology, metric, stencil and quadrature are bijective");
  passed &= expect(topology.region_halo_width() == 4U &&
                       topology.fingerprint() != 0U &&
                       boundary.fingerprint() != 0U &&
                       quadrature.physical_fingerprint() != 0U,
                   rank, "imported plans publish sealed identities");

  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 global_cell{patch.begin.x + x, patch.begin.y + y,
                               patch.begin.z + z};
        passed &= expect(
            topology.region().data[flat(patch.cells, {x, y, z})] ==
                marker[flat(kGlobal, global_cell)],
            rank, "owned region is bit-identical to imported marker");
      }
    }
  }

  std::uint64_t local_corner_links = 0U;
  for (std::size_t index = 0U; index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    const IbmInterfaceLinkMetric& metric = metrics.data[index];
    const SurfaceQuadraturePoint& point = points.data[index];
    const BoundaryStencilLink& stencil = boundary_links.data[index];
    const Real3 normal = expected_normal(link.direction);
    passed &= expect(link.global_link == 6U * link.fluid_cell +
                                           static_cast<std::uint8_t>(
                                               link.direction) &&
                         link.triangle == link.global_link &&
                         same(link.solid_to_fluid_normal, normal) &&
                         link.cartesian_control_face_area > 0.0,
                     rank, "link has canonical id, axis normal and area");
    passed &= expect(metric.global_link == link.global_link &&
                         metric.source_triangle == link.triangle &&
                         metric.physical_quadrature_area ==
                             link.cartesian_control_face_area &&
                         metric.physical_area_vector.x ==
                             link.cartesian_control_face_area * normal.x &&
                         metric.physical_area_vector.y ==
                             link.cartesian_control_face_area * normal.y &&
                         metric.physical_area_vector.z ==
                             link.cartesian_control_face_area * normal.z,
                     rank, "Cartesian face metric is exact");
    passed &= expect(point.triangle == link.triangle &&
                         point.position.x == link.wall_point.x &&
                         point.position.y == link.wall_point.y &&
                         point.position.z == link.wall_point.z &&
                         same(point.solid_to_fluid_normal, normal) &&
                         point.weight == link.cartesian_control_face_area &&
                         point.owner_cell == link.fluid_cell &&
                         stencil.topology_link == index &&
                         stencil.wall_value_row == point.wall_value_row &&
                         stencil.wall_normal_gradient_row ==
                             point.wall_normal_gradient_row,
                     rank, "quadrature reuses certified boundary wall rows");
    local_corner_links +=
        corner && link.fluid_cell == global_id({8, 8, 8}) ? 1U : 0U;
  }
  if (corner) {
    passed &= expect(sum(local_corner_links) == 3U, rank,
                     "concave Cartesian corner owns three directional links");
  }

  const double area = topology.interface_metric()
                          .conservation()
                          .physical_quadrature_area;
  const double expected_area = static_cast<double>(expected_links) * 1.0e-6;
  passed &= expect(std::abs(area - expected_area) <= 1.0e-16, rank,
                   "global metric area equals exact Cartesian face sum");
  return passed;
}

bool test_thin_sheet(const CartesianGeometryPlan& geometry,
                     const MeshPatch& patch, int rank, EBTopology& topology,
                     BoundaryStencilPlan& boundary,
                     SurfaceQuadraturePlan& quadrature) {
  const std::vector<std::uint8_t> marker = thin_sheet_marker();
  const PlanFingerprint retained_topology = topology.fingerprint();
  const PlanFingerprint retained_boundary = boundary.fingerprint();
  const PlanFingerprint retained_quadrature =
      quadrature.local_layout_fingerprint();
  ImmersedPlanLimits strict;
  const Status rejected = ImportedIbmCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, {marker.data(), marker.size()},
      kMarkerSource, strict, topology, boundary, quadrature);
  if (rejected.detail != 1304U) {
    std::cerr << "rank " << rank << " thin strict status="
              << static_cast<unsigned>(rejected.code)
              << " detail=" << rejected.detail << '\n';
  }
  bool passed = expect(
      rejected.code == StatusCode::invalid_plan && rejected.detail == 1304U &&
          topology.fingerprint() == retained_topology &&
          boundary.fingerprint() == retained_boundary &&
          quadrature.local_layout_fingerprint() == retained_quadrature,
      rank, "strict reconstruction rejects one-normal-band sheet atomically");

  ImmersedPlanLimits adaptive;
  adaptive.stencil.policy = IbmReconstructionPolicy::adaptive_order;
  const Status compiled = ImportedIbmCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, {marker.data(), marker.size()},
      kMarkerSource, adaptive, topology, boundary, quadrature);
  passed &= expect(static_cast<bool>(compiled), rank,
                   "adaptive imported marker certifies two-point thin links");
  if (!compiled) {
    std::cerr << "rank " << rank << " thin status="
              << static_cast<unsigned>(compiled.code)
              << " detail=" << compiled.detail << '\n';
    return false;
  }

  const Span<const ImmersedLink> links = topology.links();
  const QuadraticStencilPlan& reconstruction = boundary.reconstruction();
  const Span<const QuadraticStencilGroup> groups = reconstruction.groups();
  const Span<const QuadraticAffineRow> rows = reconstruction.rows();
  const Span<const GlobalCellId> donors =
      reconstruction.donor_global_cells();
  const Span<const double> weights = reconstruction.weights();
  passed &= expect(sum(links.size) == 512U && groups.size == links.size, rank,
                   "thin-sheet fixture has exact link/group count");
  std::uint64_t local_two_point = 0U;
  for (std::size_t index = 0U; index < groups.size; ++index) {
    const QuadraticStencilGroup& group = groups.data[index];
    if (group.quality.donor_count != 1U) continue;
    ++local_two_point;
    const bool storage = group.row_count == 4U && rows.size >= 4U &&
                         group.donor_begin < donors.size &&
                         group.row_begin <= rows.size - 4U;
    passed &= expect(storage &&
                         group.quality.order ==
                             IbmReconstructionOrder::linear &&
                         group.quality.normal_band_count == 1U &&
                         group.quality.fallback_reason ==
                             IbmReconstructionFallbackReason::
                                 quadratic_coverage &&
                         donors.data[group.donor_begin] ==
                             links.data[index].fluid_cell,
                     rank,
                     "thin closure is explicit one-dimensional linear fallback");
    if (!storage) continue;
    const QuadraticAffineRow& dirichlet = rows.data[group.row_begin];
    const QuadraticAffineRow& zero_normal = rows.data[group.row_begin + 1U];
    const QuadraticAffineRow& wall_value = rows.data[group.row_begin + 2U];
    const QuadraticAffineRow& wall_derivative =
        rows.data[group.row_begin + 3U];
    passed &= expect(
        dirichlet.weight_begin < weights.size &&
            zero_normal.weight_begin < weights.size &&
            wall_value.weight_begin < weights.size &&
            wall_derivative.weight_begin < weights.size &&
            std::abs(weights.data[dirichlet.weight_begin] + 1.0) <= 1.0e-13 &&
            std::abs(dirichlet.wall_value_weight - 2.0) <= 1.0e-13 &&
            std::abs(weights.data[zero_normal.weight_begin] - 1.0) <=
                1.0e-13 &&
            std::abs(zero_normal.wall_normal_gradient_weight_m + 0.001) <=
                1.0e-15 &&
            std::abs(weights.data[wall_value.weight_begin] - 1.0) <=
                1.0e-13 &&
            std::abs(weights.data[wall_derivative.weight_begin] - 2000.0) <=
                1.0e-9 &&
            std::abs(wall_derivative.wall_value_weight + 2000.0) <= 1.0e-9,
        rank, "two-point rows reproduce exact uniform half-cell formulas");
  }
  const IbmReconstructionAudit& audit = reconstruction.audit();
  passed &= expect(sum(local_two_point) == 512U &&
                       sum(audit.group_count) == 512U &&
                       sum(audit.linear_groups) == 512U &&
                       sum(audit.coverage_fallback_groups) == 512U,
                   rank, "thin-link count and adaptive audit are explicit");
  return passed;
}

bool test_selection_edge(int rank, EBTopology& topology,
                         BoundaryStencilPlan& boundary,
                         SurfaceQuadraturePlan& quadrature) {
  CartesianMeshSpec stretched = mesh_spec();
  stretched.lower.y = -0.08;
  stretched.upper.y = 0.08;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_WORLD, stretched, {}, geometry, patch)),
                       rank, "stretched selection-edge geometry compiles");
  if (!passed) return false;
  const std::vector<std::uint8_t> marker = selection_edge_marker();
  ImmersedPlanLimits limits;
  limits.stencil.policy = IbmReconstructionPolicy::adaptive_order;
  const Status status = ImportedIbmCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, {marker.data(), marker.size()},
      kMarkerSource, limits, topology, boundary, quadrature);
  passed &= expect(static_cast<bool>(status), rank,
                   "selection edge retains available tangential donors");
  if (!status) {
    std::cerr << "rank " << rank << " selection edge status="
              << static_cast<unsigned>(status.code)
              << " detail=" << status.detail << '\n';
    return false;
  }
  const Span<const ImmersedLink> links = topology.links();
  const Span<const QuadraticStencilGroup> groups =
      boundary.reconstruction().groups();
  const Span<const GlobalCellId> donors =
      boundary.reconstruction().donor_global_cells();
  const GlobalCellId target =
      6U * global_id({8, 8, 8}) +
      static_cast<std::uint8_t>(ImmersedFaceDirection::x_negative);
  std::uint64_t local_target = 0U;
  for (std::size_t index = 0U; index < links.size; ++index) {
    if (links.data[index].global_link != target) continue;
    ++local_target;
    const QuadraticStencilGroup& group = groups.data[index];
    std::int32_t minimum_y = kGlobal.y;
    std::int32_t maximum_y = -1;
    for (std::size_t donor = 0U; donor < group.quality.donor_count; ++donor) {
      const GlobalCellId id = donors.data[group.donor_begin + donor];
      const std::int32_t y = static_cast<std::int32_t>(
          (id / static_cast<GlobalCellId>(kGlobal.x)) %
          static_cast<GlobalCellId>(kGlobal.y));
      minimum_y = std::min(minimum_y, y);
      maximum_y = std::max(maximum_y, y);
    }
    passed &= expect(group.quality.donor_count >= 4U && minimum_y < maximum_y,
                     rank,
                     "edge reconstruction uses more than one tangential band");
  }
  passed &= expect(sum(links.size) == 272U && sum(local_target) == 1U, rank,
                   "selection-edge fixture has exact topology and target");
  return passed;
}

bool run(int rank) {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_WORLD, mesh_spec(), {}, geometry, patch)),
                       rank, "Cartesian fixture compiles");
  if (!passed) return false;

  EBTopology topology;
  BoundaryStencilPlan boundary;
  SurfaceQuadraturePlan quadrature;
  const std::vector<std::uint8_t> plane = plane_marker();
  passed &= compile_and_check(geometry, patch, plane, 256U, false, rank,
                              topology, boundary, quadrature);

  const std::vector<std::uint8_t> corner = corner_marker();
  passed &= compile_and_check(geometry, patch, corner, 192U, true, rank,
                              topology, boundary, quadrature);

  passed &= test_thin_sheet(geometry, patch, rank, topology, boundary,
                            quadrature);

  passed &= test_selection_edge(rank, topology, boundary, quadrature);

  const PlanFingerprint retained_topology = topology.fingerprint();
  const PlanFingerprint retained_boundary = boundary.fingerprint();
  const PlanFingerprint retained_quadrature =
      quadrature.local_layout_fingerprint();
  std::vector<std::uint8_t> invalid = corner;
  if (rank == 0) invalid.front() = 2U;
  ImmersedPlanLimits limits;
  const Status rejected = ImportedIbmCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, {invalid.data(), invalid.size()},
      kMarkerSource, limits, topology, boundary, quadrature);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       topology.fingerprint() == retained_topology &&
                       boundary.fingerprint() == retained_boundary &&
                       quadrature.local_layout_fingerprint() ==
                           retained_quadrature,
                   rank, "invalid marker rejects collectively and atomically");
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const bool passed = run(rank);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
