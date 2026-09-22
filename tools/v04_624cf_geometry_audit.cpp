// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string_view>

namespace {

using namespace hundun::v04;

constexpr std::uint64_t kExpectedFluidCells = UINT64_C(325846);
constexpr std::uint64_t kExpectedLinks = UINT64_C(87767);
constexpr double kExpectedFluidVolume = 0.0019840314423936137;
constexpr double kExpectedInletArea = 0.00099071213268392584;
constexpr double kExpectedOutletArea = 0.0030952248751428104;
constexpr double kExpectedSurfaceVolume = 0.0020246215943216885;

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.1381004, -0.05400001};
  mesh.upper = {0.3245059, 0.3088898, 0.05400001};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {160, 96, 64};
  mesh.minimum_spacing = {1.0e-12, 1.0e-12, 1.0e-12};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {UINT64_C(2000000), UINT64_C(8589934592)};
  return mesh;
}

bool near(double actual, double expected) noexcept {
  return std::abs(actual - expected) <=
         5.0e-13 * std::max({1.0, std::abs(actual), std::abs(expected)});
}

void print_status(std::string_view stage, Status status, int rank) {
  if (rank == 0) {
    std::cerr << stage << " failed: " << status_message(status) << " ("
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << ")\n";
  }
}

std::size_t flat(Int3 cells, std::int32_t x, std::int32_t y,
                 std::int32_t z) noexcept {
  return static_cast<std::size_t>(x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(z));
}

int run(const std::filesystem::path& stl, int rank) {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  Status status = CartesianGeometryCompiler::compile(
      MPI_COMM_WORLD, mesh_spec(), GeometryBudget{}, geometry, patch);
  if (!status) {
    print_status("Cartesian geometry", status, rank);
    return 3;
  }

  const StlScanBudget scan_budget{
      UINT64_C(2147483648), UINT64_C(4294967296),
      UINT64_C(33554432), UINT64_C(1000000), 1U};
  StlScanPlan scan;
  status = StlScanCompiler::compile(
      MPI_COMM_WORLD, stl.parent_path(),
      std::optional<std::filesystem::path>{stl.filename()}, geometry, patch,
      CartesianAxis::y, scan_budget, scan);
  if (!status) {
    print_status("STL scan", status, rank);
    return 4;
  }

  ImmersedSurfacePlan surface;
  status = ImmersedSurfaceCompiler::compile(scan, surface);
  if (!status) {
    print_status("immersed surface", status, rank);
    return 5;
  }

  ImmersedPlanLimits limits;
  limits.maximum_persistent_bytes_per_rank = UINT64_C(4294967296);
  limits.maximum_peak_bytes_per_rank = UINT64_C(8589934592);
  EBTopology topology;
  status = EBTopologyCompiler::compile(
      MPI_COMM_WORLD, geometry, patch, scan, surface,
      ImmersedFluidSide::inside, limits, topology);
  if (!status) {
    print_status("EB topology", status, rank);
    return 6;
  }

  const Span<const std::uint8_t> region = topology.region();
  std::uint64_t local_fluid_cells = 0U;
  std::uint64_t local_inlet_cells = 0U;
  std::uint64_t local_outlet_cells = 0U;
  double local_fluid_volume = 0.0;
  double local_inlet_area = 0.0;
  double local_outlet_area = 0.0;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    const std::int32_t gz = patch.begin.z + z;
    const double dz = geometry.z().widths().data[static_cast<std::size_t>(gz)];
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      const std::int32_t gy = patch.begin.y + y;
      const double dy =
          geometry.y().widths().data[static_cast<std::size_t>(gy)];
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        if (region.data[flat(patch.cells, x, y, z)] !=
            static_cast<std::uint8_t>(RegionFlag::fluid)) {
          continue;
        }
        const std::int32_t gx = patch.begin.x + x;
        const double dx =
            geometry.x().widths().data[static_cast<std::size_t>(gx)];
        ++local_fluid_cells;
        local_fluid_volume += dx * dy * dz;
        if (gx == 0) {
          ++local_inlet_cells;
          local_inlet_area += dy * dz;
        }
        if (gx + 1 == geometry.global_cells().x) {
          ++local_outlet_cells;
          local_outlet_area += dy * dz;
        }
      }
    }
  }

  std::uint64_t counts[4U]{local_fluid_cells, local_inlet_cells,
                           local_outlet_cells,
                           static_cast<std::uint64_t>(topology.links().size)};
  double measures[3U]{local_fluid_volume, local_inlet_area,
                      local_outlet_area};
  MPI_Allreduce(MPI_IN_PLACE, counts, 4, MPI_UINT64_T, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, measures, 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  const bool passed =
      counts[0] == kExpectedFluidCells && counts[3] == kExpectedLinks &&
      near(measures[0], kExpectedFluidVolume) &&
      near(measures[1], kExpectedInletArea) &&
      near(measures[2], kExpectedOutletArea) &&
      near(surface.closed_volume(), kExpectedSurfaceVolume);
  if (rank == 0) {
    std::cout << std::setprecision(17)
              << "{\n"
              << "  \"status\": \"" << (passed ? "pass" : "mismatch")
              << "\",\n"
              << "  \"triangles\": " << scan.triangle_count() << ",\n"
              << "  \"surface_closed_volume_m3\": "
              << surface.closed_volume() << ",\n"
              << "  \"fluid_cells\": " << counts[0] << ",\n"
              << "  \"fluid_volume_m3\": " << measures[0] << ",\n"
              << "  \"inlet_cells\": " << counts[1] << ",\n"
              << "  \"inlet_area_m2\": " << measures[1] << ",\n"
              << "  \"outlet_cells\": " << counts[2] << ",\n"
              << "  \"outlet_area_m2\": " << measures[2] << ",\n"
              << "  \"immersed_links\": " << counts[3] << "\n"
              << "}\n";
  }
  return passed ? 0 : 7;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int result = 0;
  if (argc != 2) {
    if (rank == 0) {
      std::cerr << "usage: v04_624cf_geometry_audit /path/to/cf1.stl\n";
    }
    result = 2;
  } else {
    result = run(std::filesystem::absolute(argv[1]).lexically_normal(), rank);
  }
  MPI_Finalize();
  return result;
}
