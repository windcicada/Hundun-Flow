// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr Int3 kGlobal{8, 8, 8};
constexpr StlScanBudget kScanBudget{UINT64_C(67108864), UINT64_C(134217728),
                                    UINT64_C(1000000), UINT64_C(10000), 1U};

static_assert(!std::is_copy_constructible_v<EBTopology>);
static_assert(std::is_nothrow_move_constructible_v<EBTopology>);

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

std::uint64_t packed(Status status) noexcept {
  return (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
}

bool identical(std::uint64_t value) {
  std::uint64_t minimum = value;
  std::uint64_t maximum = value;
  return MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-2.0, -2.0, -2.0};
  mesh.upper = {2.0, 2.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kGlobal;
  mesh.minimum_spacing = {1.0e-10, 1.0e-10, 1.0e-10};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(1000000), UINT64_C(1073741824)};
  return mesh;
}

std::array<TriangleInput, 12U> cube_triangles() {
  constexpr double l = -1.0;
  constexpr double h = 1.0;
  return {{
      {{l, l, l}, {l, l, h}, {l, h, h}},
      {{l, l, l}, {l, h, h}, {l, h, l}},
      {{h, l, l}, {h, h, l}, {h, h, h}},
      {{h, l, l}, {h, h, h}, {h, l, h}},
      {{l, l, l}, {h, l, l}, {h, l, h}},
      {{l, l, l}, {h, l, h}, {l, l, h}},
      {{l, h, l}, {l, h, h}, {h, h, h}},
      {{l, h, l}, {h, h, h}, {h, h, l}},
      {{l, l, l}, {l, h, l}, {h, h, l}},
      {{l, l, l}, {h, h, l}, {h, l, l}},
      {{l, l, h}, {h, l, h}, {h, h, h}},
      {{l, l, h}, {h, h, h}, {l, h, h}},
  }};
}

std::array<TriangleInput, 4U> tetrahedron_triangles() {
  const Real3 a{-1.0, -1.0, -1.0};
  const Real3 b{1.0, -1.0, -1.0};
  const Real3 c{-1.0, 1.0, -1.0};
  const Real3 d{-1.0, -1.0, 1.0};
  return {{{a, c, b}, {a, b, d}, {a, d, c}, {b, c, d}}};
}

std::size_t flat(Int3 shape, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(cell.z));
}

bool inside(Int3 cell) noexcept {
  return cell.x >= 0 && cell.x < kGlobal.x && cell.y >= 0 &&
         cell.y < kGlobal.y && cell.z >= 0 && cell.z < kGlobal.z;
}

bool raw_inside_cube(Int3 cell) noexcept {
  return cell.x >= 2 && cell.x <= 5 && cell.y >= 2 && cell.y <= 5 &&
         cell.z >= 2 && cell.z <= 5;
}

bool is_fluid(ImmersedFluidSide side, Int3 cell) noexcept {
  return side == ImmersedFluidSide::inside ? raw_inside_cube(cell)
                                            : !raw_inside_cube(cell);
}

constexpr std::array<Int3, 6U> kDelta{{
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
    {0, 1, 0},  {0, 0, -1}, {0, 0, 1},
}};

std::vector<std::uint64_t> expected_link_ids(ImmersedFluidSide side) {
  std::vector<std::uint64_t> result(
      static_cast<std::size_t>(kGlobal.x * kGlobal.y * kGlobal.z) * 6U,
      std::numeric_limits<std::uint64_t>::max());
  for (std::int32_t z = 0; z < kGlobal.z; ++z) {
    for (std::int32_t y = 0; y < kGlobal.y; ++y) {
      for (std::int32_t x = 0; x < kGlobal.x; ++x) {
        const Int3 fluid{x, y, z};
        if (!is_fluid(side, fluid)) {
          continue;
        }
        for (std::size_t direction = 0U; direction < kDelta.size();
             ++direction) {
          const Int3 delta = kDelta[direction];
          const Int3 solid{fluid.x + delta.x, fluid.y + delta.y,
                           fluid.z + delta.z};
          if (inside(solid) && !is_fluid(side, solid)) {
            const std::uint64_t key =
                6U * flat(kGlobal, fluid) + direction;
            result[static_cast<std::size_t>(key)] = key;
          }
        }
      }
    }
  }
  return result;
}

Real3 centre(const CartesianGeometryPlan& geometry, Int3 cell) {
  return {geometry.x().centres().data[static_cast<std::size_t>(cell.x)],
          geometry.y().centres().data[static_cast<std::size_t>(cell.y)],
          geometry.z().centres().data[static_cast<std::size_t>(cell.z)]};
}

bool close(double left, double right, double tolerance = 1.0e-13) noexcept {
  return std::abs(left - right) <= tolerance;
}

bool validate_topology(const CartesianGeometryPlan& geometry,
                       const MeshPatch& patch, ImmersedFluidSide side,
                       const EBTopology& topology, int rank) {
  bool passed = expect(topology.fluid_side() == side &&
                           topology.region_halo_width() == 4U &&
                           topology.geometry_revision() ==
                               geometry.topology_revision() &&
                           topology.fingerprint() != 0U &&
                           topology.lowest_failing_rank() == -1,
                       rank, "topology publishes identity and four-cell region halo");
  const std::size_t local_count = static_cast<std::size_t>(patch.cells.x) *
                                  static_cast<std::size_t>(patch.cells.y) *
                                  static_cast<std::size_t>(patch.cells.z);
  passed &= expect(topology.region().size == local_count, rank,
                   "topology persists one region byte per owned cell");

  std::uint64_t local_fluid = 0U;
  std::vector<std::uint32_t> expected_interfaces;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 local{x, y, z};
        const Int3 global{patch.begin.x + x, patch.begin.y + y,
                          patch.begin.z + z};
        const std::size_t index = flat(patch.cells, local);
        const bool fluid = is_fluid(side, global);
        local_fluid += fluid ? 1U : 0U;
        passed &= expect(topology.region().data[index] ==
                             static_cast<std::uint8_t>(
                                 fluid ? RegionFlag::fluid : RegionFlag::solid),
                         rank, "owned region follows selected fluid side");
        passed &= expect(topology.is_fluid_global(global) == fluid, rank,
                         "public region query returns each owned cell");
        bool interface_cell = false;
        if (fluid) {
          for (const Int3 delta : kDelta) {
            const Int3 neighbor{global.x + delta.x, global.y + delta.y,
                                global.z + delta.z};
            interface_cell |= inside(neighbor) && !is_fluid(side, neighbor);
          }
        }
        if (interface_cell) {
          expected_interfaces.push_back(static_cast<std::uint32_t>(index));
        }
      }
    }
  }
  std::uint64_t global_fluid = local_fluid;
  passed &= expect(MPI_Allreduce(MPI_IN_PLACE, &global_fluid, 1, MPI_UINT64_T,
                                 MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS &&
                       global_fluid ==
                           (side == ImmersedFluidSide::inside ? 64U : 448U),
                   rank, "selected fluid side has the hand-counted volume");
  passed &= expect(topology.interface_cells().size ==
                           expected_interfaces.size() &&
                       std::equal(expected_interfaces.begin(),
                                  expected_interfaces.end(),
                                  topology.interface_cells().data),
                   rank, "interface cells are unique local-flat sorted owners");
  for (std::int32_t z = std::max(0, patch.begin.z - 4);
       z < std::min(kGlobal.z, patch.begin.z + patch.cells.z + 4); ++z) {
    for (std::int32_t y = std::max(0, patch.begin.y - 4);
         y < std::min(kGlobal.y, patch.begin.y + patch.cells.y + 4); ++y) {
      for (std::int32_t x = std::max(0, patch.begin.x - 4);
           x < std::min(kGlobal.x, patch.begin.x + patch.cells.x + 4); ++x) {
        const Int3 global{x, y, z};
        passed &= expect(topology.is_fluid_global(global) ==
                             is_fluid(side, global),
                         rank,
                         "public region query exposes the four-cell halo");
      }
    }
  }

  const std::vector<std::uint64_t> expected_ids = expected_link_ids(side);
  std::vector<unsigned int> global_id_owners(
      static_cast<std::size_t>(kGlobal.x * kGlobal.y * kGlobal.z) * 6U, 0U);
  const Span<const ImmersedLink> links = topology.links();
  std::uint64_t local_links = links.size;
  for (std::size_t link_index = 0U; link_index < links.size; ++link_index) {
    const ImmersedLink& link = links.data[link_index];
    const std::size_t direction = static_cast<std::size_t>(link.direction);
    const Int3 fluid = link.fluid_global_index;
    const Int3 delta = kDelta[direction];
    const Int3 solid{fluid.x + delta.x, fluid.y + delta.y,
                     fluid.z + delta.z};
    const std::uint64_t expected_id =
        expected_ids[6U * flat(kGlobal, fluid) + direction];
    passed &= expect(direction < kDelta.size() && inside(fluid) &&
                         inside(solid) && is_fluid(side, fluid) &&
                         !is_fluid(side, solid) &&
                         link.solid_global_index.x == solid.x &&
                         link.solid_global_index.y == solid.y &&
                         link.solid_global_index.z == solid.z &&
                         link.fluid_cell == flat(kGlobal, fluid) &&
                         link.solid_cell == flat(kGlobal, solid) &&
                         link.global_link == expected_id &&
                         expected_id < global_id_owners.size(),
                     rank, "link ownership, direction, cells, and global id are exact");
    if (expected_id < global_id_owners.size()) {
      ++global_id_owners[static_cast<std::size_t>(expected_id)];
    }
    passed &= expect(
        link.fluid_local_index.x == fluid.x - patch.begin.x &&
            link.fluid_local_index.y == fluid.y - patch.begin.y &&
            link.fluid_local_index.z == fluid.z - patch.begin.z &&
            link.solid_local_index.x == solid.x - patch.begin.x &&
            link.solid_local_index.y == solid.y - patch.begin.y &&
            link.solid_local_index.z == solid.z - patch.begin.z,
        rank, "link local indices address owned fluid and four-cell halo solid");
    const Real3 fluid_centre = centre(geometry, fluid);
    const Real3 solid_centre = centre(geometry, solid);
    const Real3 midpoint{0.5 * (fluid_centre.x + solid_centre.x),
                         0.5 * (fluid_centre.y + solid_centre.y),
                         0.5 * (fluid_centre.z + solid_centre.z)};
    passed &= expect(link.triangle < 12U && close(link.wall_point.x, midpoint.x) &&
                         close(link.wall_point.y, midpoint.y) &&
                         close(link.wall_point.z, midpoint.z) &&
                         close(link.solid_to_fluid_normal.x,
                               -static_cast<double>(delta.x)) &&
                         close(link.solid_to_fluid_normal.y,
                               -static_cast<double>(delta.y)) &&
                         close(link.solid_to_fluid_normal.z,
                               -static_cast<double>(delta.z)),
                     rank, "segment wall point and solid-to-fluid normal are exact");
    constexpr double link_face_area = 0.25;
    passed &= expect(
        close(link.surface_patch_centroid.x, link.wall_point.x) &&
            close(link.surface_patch_centroid.y, link.wall_point.y) &&
            close(link.surface_patch_centroid.z, link.wall_point.z) &&
            close(link.surface_measure_vector.x,
                  link_face_area * link.solid_to_fluid_normal.x) &&
            close(link.surface_measure_vector.y,
                  link_face_area * link.solid_to_fluid_normal.y) &&
            close(link.surface_measure_vector.z,
                  link_face_area * link.solid_to_fluid_normal.z),
        rank,
        "link control-face measure and centroid are deterministic geometry");
  }
  passed &= expect(MPI_Allreduce(MPI_IN_PLACE, &local_links, 1, MPI_UINT64_T,
                                 MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS &&
                       local_links == 96U,
                   rank, "cube owns exactly 96 physical-domain fluid-solid links");
  bool exact_key_owners =
      MPI_Allreduce(MPI_IN_PLACE, global_id_owners.data(),
                    static_cast<int>(global_id_owners.size()), MPI_UNSIGNED,
                    MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS;
  for (std::size_t key = 0U; key < global_id_owners.size(); ++key) {
    const bool expected =
        expected_ids[key] != std::numeric_limits<std::uint64_t>::max();
    exact_key_owners &= global_id_owners[key] == (expected ? 1U : 0U);
  }
  passed &= expect(exact_key_owners, rank,
                   "global link keys are unique and partition independent");
  passed &= expect(identical(topology.fingerprint()), rank,
                   "topology fingerprint is identical on every rank");
  return all_true(passed);
}

bool validate_partition_independence(
    const CartesianGeometryPlan& geometry, const StlScanPlan& scan,
    const ImmersedSurfacePlan& surface, ImmersedPlanLimits limits,
    ImmersedFluidSide side, PlanFingerprint distributed_fingerprint,
    int rank) {
  MeshPatch full_patch{{0, 0, 0}, kGlobal, {1, 1, 1}, {0, 0, 0}};
  StlScanPlan full_scan;
  const auto triangles = cube_triangles();
  bool passed = expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                           geometry, full_patch,
                           {triangles.data(), triangles.size()},
                           CartesianAxis::y, kScanBudget, full_scan)),
                       rank, "full-domain scan compiles on MPI_COMM_SELF");
  ImmersedSurfacePlan full_surface;
  passed &= expect(static_cast<bool>(ImmersedSurfaceCompiler::compile(
                       full_scan, full_surface)) &&
                       full_surface.fingerprint() == surface.fingerprint(),
                   rank, "full-domain canonical surface is unchanged");
  EBTopology full;
  limits.maximum_persistent_bytes_per_rank = UINT64_C(1048576);
  limits.maximum_peak_bytes_per_rank = UINT64_C(2097152);
  passed &= expect(static_cast<bool>(EBTopologyCompiler::compile(
                       MPI_COMM_SELF, geometry, full_patch, full_scan,
                       full_surface, side, limits, full)) &&
                       full.fingerprint() == distributed_fingerprint,
                   rank,
                   "distributed and full-domain topology fingerprints agree");
  const Int3 halo_probe{std::min(kGlobal.x - 1,
                                full_patch.begin.x + full_patch.cells.x / 2),
                       3, 3};
  passed &= expect(full.is_fluid_global(halo_probe) ==
                       is_fluid(side, halo_probe) &&
                       !full.is_fluid_global({-1, 3, 3}),
                   rank, "public region query rejects physical-domain exterior");
  (void)scan;
  return all_true(passed);
}

bool validate_oblique_surface_normal(const CartesianGeometryPlan& geometry,
                                     const MeshPatch& patch,
                                     ImmersedPlanLimits limits, int rank) {
  const auto triangles = tetrahedron_triangles();
  StlScanPlan scan;
  bool passed = expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                           geometry, patch,
                           {triangles.data(), triangles.size()},
                           CartesianAxis::y, kScanBudget, scan)),
                       rank, "oblique tetrahedron scan compiles");
  ImmersedSurfacePlan surface;
  passed &= expect(static_cast<bool>(ImmersedSurfaceCompiler::compile(
                       scan, surface)),
                   rank, "oblique tetrahedron surface compiles");
  EBTopology topology;
  passed &= expect(static_cast<bool>(EBTopologyCompiler::compile(
                       MPI_COMM_WORLD, geometry, patch, scan, surface,
                       ImmersedFluidSide::outside, limits, topology)),
                   rank, "oblique tetrahedron topology compiles");
  if (!all_true(passed)) {
    return false;
  }
  const Span<const SurfaceTriangle> canonical = surface.triangles();
  const Span<const ImmersedLink> links = topology.links();
  std::uint64_t local_oblique = 0U;
  for (std::size_t index = 0U; index < links.size; ++index) {
    const ImmersedLink& link = links.data[index];
    if (link.triangle >= canonical.size) {
      passed &= expect(false, rank, "topology link references a canonical triangle");
      continue;
    }
    const Real3 normal =
        canonical.data[static_cast<std::size_t>(link.triangle)]
            .geometric_outward_normal;
    const bool oblique = std::abs(normal.x) > 0.5 &&
                         std::abs(normal.y) > 0.5 &&
                         std::abs(normal.z) > 0.5;
    if (oblique) {
      ++local_oblique;
      passed &= expect(close(link.solid_to_fluid_normal.x, normal.x) &&
                           close(link.solid_to_fluid_normal.y, normal.y) &&
                           close(link.solid_to_fluid_normal.z, normal.z),
                       rank,
                       "outside-fluid link uses canonical oblique surface normal");
    }
  }
  std::uint64_t global_oblique = local_oblique;
  passed &= expect(MPI_Allreduce(MPI_IN_PLACE, &global_oblique, 1,
                                 MPI_UINT64_T, MPI_SUM,
                                 MPI_COMM_WORLD) == MPI_SUCCESS &&
                       global_oblique > 0U,
                   rank, "oblique surface owns at least one Cartesian link");
  return all_true(passed);
}

bool validate_narrow_partition_forwarding(
    const CartesianGeometryPlan& geometry, ImmersedPlanLimits limits,
    PlanFingerprint reference_fingerprint, int rank, int size) {
  if (size != 4) {
    return true;
  }
  const MeshPatch narrow_patch{{2 * rank, 0, 0}, {2, 8, 8}, {4, 1, 1},
                               {rank, 0, 0}};
  const auto triangles = cube_triangles();
  StlScanPlan narrow_scan;
  bool passed = expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                           geometry, narrow_patch,
                           {triangles.data(), triangles.size()},
                           CartesianAxis::y, kScanBudget, narrow_scan)),
                       rank, "two-cell-wide partition scan compiles");
  ImmersedSurfacePlan narrow_surface;
  passed &= expect(static_cast<bool>(ImmersedSurfaceCompiler::compile(
                       narrow_scan, narrow_surface)),
                   rank, "two-cell-wide partition surface compiles");
  limits.maximum_persistent_bytes_per_rank = UINT64_C(1048576);
  limits.maximum_peak_bytes_per_rank = UINT64_C(2097152);
  EBTopology narrow;
  passed &= expect(static_cast<bool>(EBTopologyCompiler::compile(
                       MPI_COMM_WORLD, geometry, narrow_patch, narrow_scan,
                       narrow_surface, ImmersedFluidSide::outside, limits,
                       narrow)),
                   rank, "two-cell-wide partition topology compiles");
  if (!all_true(passed)) {
    return false;
  }
  for (std::int32_t z = 0; z < kGlobal.z; ++z) {
    for (std::int32_t y = 0; y < kGlobal.y; ++y) {
      for (std::int32_t x = std::max(0, narrow_patch.begin.x - 4);
           x < std::min(kGlobal.x,
                        narrow_patch.begin.x + narrow_patch.cells.x + 4);
           ++x) {
        passed &= expect(narrow.is_fluid_global({x, y, z}) ==
                             is_fluid(ImmersedFluidSide::outside, {x, y, z}),
                         rank,
                         "four-cell halo forwards across multiple neighbors");
      }
    }
  }
  passed &= expect(narrow.fingerprint() == reference_fingerprint &&
                       identical(narrow.fingerprint()),
                   rank,
                   "narrow and balanced decompositions publish one topology identity");
  return all_true(passed);
}

bool run(int rank, int size) {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_WORLD, mesh_spec(), GeometryBudget{},
                           geometry, patch)),
                       rank, "topology geometry compiles");
  const auto triangles = cube_triangles();
  StlScanPlan scan;
  passed &= expect(static_cast<bool>(StlScanCompiler::compile_triangles(
                       geometry, patch, {triangles.data(), triangles.size()},
                       CartesianAxis::y, kScanBudget, scan)),
                   rank, "local cube parity scan compiles");
  ImmersedSurfacePlan surface;
  passed &= expect(static_cast<bool>(ImmersedSurfaceCompiler::compile(
                       scan, surface)),
                   rank, "canonical cube surface compiles");
  if (!all_true(passed)) {
    return false;
  }

  ImmersedPlanLimits limits;
  std::uint64_t maximum_local_region =
      static_cast<std::uint64_t>(patch.cells.x) * patch.cells.y * patch.cells.z;
  passed &= expect(MPI_Allreduce(MPI_IN_PLACE, &maximum_local_region, 1,
                                 MPI_UINT64_T, MPI_MAX,
                                 MPI_COMM_WORLD) == MPI_SUCCESS,
                   rank, "maximum owned-region size is collective");
  limits.maximum_persistent_bytes_per_rank = UINT64_C(1048576);
  limits.maximum_peak_bytes_per_rank = UINT64_C(2097152);
  limits.maximum_local_links = 1024U;
  EBTopology outside;
  const Status outside_status = EBTopologyCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, scan, surface,
      ImmersedFluidSide::outside, limits, outside);
  passed &= expect(static_cast<bool>(outside_status), rank,
                   "outside-fluid topology compiles");
  if (all_true(passed)) {
    passed &= validate_topology(geometry, patch, ImmersedFluidSide::outside,
                                outside, rank);
    passed &= validate_partition_independence(
        geometry, scan, surface, limits, ImmersedFluidSide::outside,
        outside.fingerprint(), rank);
    passed &= validate_narrow_partition_forwarding(
        geometry, limits, outside.fingerprint(), rank, size);
  }

  EBTopology inside;
  const Status inside_status = EBTopologyCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, scan, surface,
      ImmersedFluidSide::inside, limits, inside);
  passed &= expect(static_cast<bool>(inside_status), rank,
                   "inside-fluid topology compiles");
  if (all_true(passed)) {
    passed &= validate_topology(geometry, patch, ImmersedFluidSide::inside,
                                inside, rank);
    passed &= validate_partition_independence(
        geometry, scan, surface, limits, ImmersedFluidSide::inside,
        inside.fingerprint(), rank);
    passed &= expect(outside.fingerprint() != inside.fingerprint(), rank,
                     "fluid-side selection changes topology identity");
  }
  passed &= validate_oblique_surface_normal(geometry, patch, limits, rank);

  const PlanFingerprint retained_fingerprint = outside.fingerprint();
  const std::size_t retained_region = outside.region().size;
  const std::size_t retained_links = outside.links().size;
  ImmersedPlanLimits asymmetric = limits;
  if (rank == size - 1) {
    const std::uint64_t local_region =
        static_cast<std::uint64_t>(patch.cells.x) * patch.cells.y *
        patch.cells.z;
    asymmetric.maximum_persistent_bytes_per_rank = local_region - 1U;
    asymmetric.maximum_peak_bytes_per_rank = local_region - 1U;
  }
  const Status rejected = EBTopologyCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, scan, surface,
      ImmersedFluidSide::outside, asymmetric, outside);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       identical(packed(rejected)) &&
                       outside.lowest_failing_rank() == size - 1 &&
                       outside.fingerprint() == retained_fingerprint &&
                       outside.region().size == retained_region &&
                       outside.links().size == retained_links,
                   rank,
                   "rank-local owned-region gate rejects collectively without publishing");
  if (size > 1) {
    ImmersedPlanLimits mismatched = limits;
    if (rank == size - 1) {
      ++mismatched.maximum_peak_bytes_per_rank;
    }
    const Status mismatch = EBTopologyCompiler::compile(
        MPI_COMM_WORLD, geometry, patch, scan, surface,
        ImmersedFluidSide::outside, mismatched, outside);
    passed &= expect(mismatch.code == StatusCode::invalid_plan &&
                         identical(packed(mismatch)) &&
                         outside.lowest_failing_rank() == size - 1 &&
                         outside.fingerprint() == retained_fingerprint,
                     rank,
                     "rank-local topology budget mismatch rejects atomically");

    MeshPatch wrong_owner = patch;
    if (rank == size - 1) {
      wrong_owner.process_coord = {};
    }
    const Status ownership = EBTopologyCompiler::compile(
        MPI_COMM_WORLD, geometry, wrong_owner, scan, surface,
        ImmersedFluidSide::outside, limits, outside);
    passed &= expect(ownership.code == StatusCode::invalid_plan &&
                         identical(packed(ownership)) &&
                         outside.lowest_failing_rank() == size - 1 &&
                         outside.fingerprint() == retained_fingerprint,
                     rank,
                     "rank-local patch ownership mismatch rejects atomically");
  }
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  bool passed = expect(size == 1 || size == 2 || size == 4, rank,
                       "EB topology MPI test runs at 1, 2, or 4 ranks");
  passed &= run(rank, size);
  passed = all_true(passed);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
