// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

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

std::uint64_t sum_u64(std::uint64_t local) {
  std::uint64_t global{};
  MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  return global;
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
                         group.quality.reach <= 4U &&
                         group.quality.condition_estimate <= 1.0e8,
                     "every boundary group satisfies the strict contract");
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
