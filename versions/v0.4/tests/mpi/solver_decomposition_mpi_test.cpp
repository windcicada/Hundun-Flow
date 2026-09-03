// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_execution.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
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

bool all_true(bool local, MPI_Comm communicator) {
  const int value = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         result != 0;
}

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-1.0, -0.75, 0.25};
  mesh.upper = {2.0, 1.25, 1.75};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {17, 11, 7};
  mesh.minimum_spacing = {1.0e-8, 1.0e-8, 1.0e-8};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {20000U, 1U << 28U};
  return mesh;
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec value;
  value.flow_kind = BoundaryKind::no_slip_wall;
  value.thermal_kind = BoundaryKind::adiabatic_wall;
  value.mach_limit = 0.95;
  return value;
}

ValidatedModel model_spec() {
  ValidatedModel value;
  value.fingerprint = 0x9382574U;
  value.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : value.boundaries) {
    face = wall();
  }
  value.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
  value.boundaries[0U].thermal_kind = BoundaryKind::none;
  value.boundaries[0U].velocity = {0.4, 0.0, 0.0};
  value.boundaries[0U].temperature = 300.0;
  value.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
  value.boundaries[1U].thermal_kind = BoundaryKind::none;
  value.boundaries[1U].pressure = 101325.0;
  value.schemes = SchemeSpec{};
  value.time = TimeControlSpec{};
  return value;
}

struct OwnedField {
  std::vector<double> bytes;
  FieldView view{};
};

OwnedField make_field(FieldId field, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision) {
  OwnedField result;
  const std::size_t stride_y =
      static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t stride_z =
      stride_y * static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t component_stride =
      stride_z * static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.bytes.assign(component_stride * components, -9.87654321e99);
  result.view.base =
      result.bytes.data() + ghosts + static_cast<std::size_t>(ghosts) * stride_y +
      static_cast<std::size_t>(ghosts) * stride_z;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = stride_y;
  result.view.stride_z = stride_z;
  result.view.component_stride = component_stride;
  result.view.field = field;
  result.view.revision = revision;
  result.view.storage_identity = static_cast<StorageIdentity>(field) + 101U;
  result.view.revision_domain = 5001U;
  return result;
}

double rho_oracle(Int3 global) {
  return 1.0 + 0.015 * global.x + 0.009 * global.y + 0.007 * global.z;
}

Real3 velocity_oracle(Int3 global) {
  return {0.4 + 0.006 * global.x, -0.23 + 0.004 * global.y,
          0.17 - 0.003 * global.z};
}

void fill_owned(OwnedField& rho, OwnedField& velocity,
                const MeshPatch& patch) {
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 local{x, y, z};
        const Int3 global{patch.begin.x + x, patch.begin.y + y,
                          patch.begin.z + z};
        rho.view.unchecked(local, 0U) = rho_oracle(global);
        const Real3 u = velocity_oracle(global);
        velocity.view.unchecked(local, 0U) = u.x;
        velocity.view.unchecked(local, 1U) = u.y;
        velocity.view.unchecked(local, 2U) = u.z;
      }
    }
  }
}

void fill_physical_ghosts(OwnedField& rho, OwnedField& velocity,
                          const CartesianGeometryPlan& geometry,
                          const MeshPatch& patch) {
  constexpr std::int32_t reach = 2;
  const Int3 global_cells = geometry.global_cells();
  for (std::int32_t z = -reach; z < patch.cells.z + reach; ++z) {
    for (std::int32_t y = -reach; y < patch.cells.y + reach; ++y) {
      for (std::int32_t x = -reach; x < patch.cells.x + reach; ++x) {
        const bool owned = x >= 0 && x < patch.cells.x && y >= 0 &&
                           y < patch.cells.y && z >= 0 && z < patch.cells.z;
        if (owned) {
          continue;
        }
        const Int3 global{patch.begin.x + x, patch.begin.y + y,
                          patch.begin.z + z};
        const bool physical = global.x < 0 || global.y < 0 || global.z < 0 ||
                              global.x >= global_cells.x ||
                              global.y >= global_cells.y ||
                              global.z >= global_cells.z;
        if (!physical) {
          continue;
        }
        rho.view.unchecked({x, y, z}, 0U) = rho_oracle(global);
        const Real3 u = velocity_oracle(global);
        velocity.view.unchecked({x, y, z}, 0U) = u.x;
        velocity.view.unchecked({x, y, z}, 1U) = u.y;
        velocity.view.unchecked({x, y, z}, 2U) = u.z;
      }
    }
  }
}

KernelBox interior_box(Int3 cells, std::uint8_t reach) {
  const auto width = static_cast<std::int32_t>(reach);
  return {{width, width, width},
          {std::max(0, cells.x - 2 * width),
           std::max(0, cells.y - 2 * width),
           std::max(0, cells.z - 2 * width)}};
}

bool on_shell(Int3 index, Int3 cells, std::uint8_t reach) {
  const auto width = static_cast<std::int32_t>(reach);
  return index.x < width || index.y < width || index.z < width ||
         index.x >= cells.x - width || index.y >= cells.y - width ||
         index.z >= cells.z - width;
}

std::uint64_t owned_face_count(KernelBox box) {
  const auto cx = static_cast<std::uint64_t>(box.cells.x);
  const auto cy = static_cast<std::uint64_t>(box.cells.y);
  const auto cz = static_cast<std::uint64_t>(box.cells.z);
  return (cx + (box.begin.x == 0 ? 1U : 0U)) * cy * cz +
         cx * (cy + (box.begin.y == 0 ? 1U : 0U)) * cz +
         cx * cy * (cz + (box.begin.z == 0 ? 1U : 0U));
}

bool compare_rank_interfaces(const MeshPatch& patch, ConstFaceFluxView flux,
                             MPI_Comm communicator) {
  int size = 0;
  MPI_Comm_size(communicator, &size);
  std::vector<std::uint64_t> local;
  const auto append = [&](CartesianAxis axis, Int3 global_face,
                          double value) {
    local.push_back(static_cast<std::uint64_t>(axis));
    local.push_back(static_cast<std::uint64_t>(global_face.x));
    local.push_back(static_cast<std::uint64_t>(global_face.y));
    local.push_back(static_cast<std::uint64_t>(global_face.z));
    local.push_back(bits(value));
  };
  if (patch.process_coord.x > 0) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        append(CartesianAxis::x,
               {patch.begin.x, patch.begin.y + y, patch.begin.z + z},
               flux.x.unchecked({0, y, z}));
      }
    }
  }
  if (patch.process_coord.x + 1 < patch.process_grid.x) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        append(CartesianAxis::x,
               {patch.begin.x + patch.cells.x, patch.begin.y + y,
                patch.begin.z + z},
               flux.x.unchecked({patch.cells.x, y, z}));
      }
    }
  }
  if (patch.process_coord.y > 0) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        append(CartesianAxis::y,
               {patch.begin.x + x, patch.begin.y, patch.begin.z + z},
               flux.y.unchecked({x, 0, z}));
      }
    }
  }
  if (patch.process_coord.y + 1 < patch.process_grid.y) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        append(CartesianAxis::y,
               {patch.begin.x + x, patch.begin.y + patch.cells.y,
                patch.begin.z + z},
               flux.y.unchecked({x, patch.cells.y, z}));
      }
    }
  }
  if (patch.process_coord.z > 0) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        append(CartesianAxis::z,
               {patch.begin.x + x, patch.begin.y + y, patch.begin.z},
               flux.z.unchecked({x, y, 0}));
      }
    }
  }
  if (patch.process_coord.z + 1 < patch.process_grid.z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        append(CartesianAxis::z,
               {patch.begin.x + x, patch.begin.y + y,
                patch.begin.z + patch.cells.z},
               flux.z.unchecked({x, y, patch.cells.z}));
      }
    }
  }

  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(size));
  if (MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                    communicator) != MPI_SUCCESS) {
    return false;
  }
  std::vector<int> displacements(static_cast<std::size_t>(size));
  int total = 0;
  for (int rank = 0; rank < size; ++rank) {
    displacements[static_cast<std::size_t>(rank)] = total;
    total += counts[static_cast<std::size_t>(rank)];
  }
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(total));
  if (MPI_Allgatherv(local.data(), local_count, MPI_UINT64_T, gathered.data(),
                     counts.data(), displacements.data(), MPI_UINT64_T,
                     communicator) != MPI_SUCCESS) {
    return false;
  }
  // Compare complete records rather than the sorted words above.  Every
  // internal interface record must occur exactly twice with an identical bit
  // pattern; physical faces were not inserted.
  for (std::size_t first = 0U; first < static_cast<std::size_t>(total);
       first += 5U) {
    const std::uint64_t* const record = gathered.data() + first;
    int matches = 0;
    for (std::size_t candidate = 0U;
         candidate < static_cast<std::size_t>(total); candidate += 5U) {
      bool same = true;
      for (std::size_t word = 0U; word < 5U; ++word) {
        same = same && gathered[candidate + word] == record[word];
      }
      matches += same ? 1 : 0;
    }
    if (matches != 2) {
      return false;
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_WORLD, mesh_spec(), GeometryBudget{},
                           geometry, patch)),
                       rank, "non-divisible Cartesian geometry compiles");
  constexpr FieldId kRho = 0U;
  constexpr FieldId kVelocity = 1U;
  constexpr FieldId kDivergence = 2U;
  OwnedField rho = make_field(kRho, patch.cells, 1U, 2U, 11U);
  OwnedField velocity = make_field(kVelocity, patch.cells, 3U, 2U, 12U);
  OwnedField divergence = make_field(kDivergence, patch.cells, 1U, 0U, 13U);
  fill_owned(rho, velocity, patch);

  const std::array halo_specs{HaloFieldSpec{kRho, 2U, 1U},
                              HaloFieldSpec{kVelocity, 2U, 3U}};
  HaloEngine halo;
  passed &= expect(static_cast<bool>(halo.reserve(
                       MPI_COMM_WORLD, patch,
                       Span<const HaloFieldSpec>{halo_specs.data(),
                                                 halo_specs.size()})),
                   rank, "cell halo reserves before kernel execution");

  CartesianKernelPlan plan;
  FieldRegistry registry;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  passed &= expect(static_cast<bool>(BoundaryCompiler::compile(
                       MPI_COMM_WORLD, model_spec(), geometry, patch, registry,
                       boundary, schemes, time)),
                   rank, "boundary and scheme plans compile");
  passed &= expect(static_cast<bool>(CartesianKernelPlan::compile(
                       schemes, geometry, patch, boundary, plan)),
                   rank, "kernel plan compiles with global metric offsets");
  FaceFluxStorage storage;
  passed &= expect(static_cast<bool>(
                       FaceFluxStorage::allocate_workspace(patch.cells, 1U,
                                                           storage)),
                   rank, "rank-local face storage allocates");
  FaceFluxView flux;
  passed &= expect(static_cast<bool>(storage.workspace_view(0U, 21U, flux)), rank,
                   "pending face storage view is available");

  std::array<FieldView, 2U> halo_fields{rho.view, velocity.view};
  const std::array<ConstFieldView, 2U> reads{as_const(rho.view),
                                            as_const(velocity.view)};
  KernelCounters counters;
  KernelInvocation invocation{
      Span<const ConstFieldView>{reads.data(), reads.size()},
      {}, interior_box(patch.cells, 2U), 0U, 0U, 1U, 0U, &counters};

  HaloTicket ticket;
  passed &= expect(static_cast<bool>(halo.begin(
                       91U,
                       Span<const FieldView>{halo_fields.data(),
                                             halo_fields.size()},
                       ticket)),
                   rank, "halo begins before interior kernel");
  passed &= expect(static_cast<bool>(reconstruct_mass_flux(
                       plan, invocation, flux)),
                   rank, "interior face reconstruction overlaps halo");
  passed &= expect(static_cast<bool>(halo.finish(
                       ticket, Span<FieldView>{halo_fields.data(),
                                               halo_fields.size()})),
                   rank, "halo finishes before shell kernel");
  rho.view = halo_fields[0U];
  velocity.view = halo_fields[1U];
  fill_physical_ghosts(rho, velocity, geometry, patch);
  const std::array<ConstFieldView, 2U> shell_reads{as_const(rho.view),
                                                  as_const(velocity.view)};
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        if (!on_shell({x, y, z}, patch.cells, 2U)) {
          continue;
        }
        KernelInvocation shell{
            Span<const ConstFieldView>{shell_reads.data(), shell_reads.size()},
            {}, {{x, y, z}, {1, 1, 1}}, 0U, 0U, 1U, 0U, &counters};
        passed &= static_cast<bool>(reconstruct_mass_flux(plan, shell, flux));
      }
    }
  }

  const ConstFaceFluxView const_flux = as_const(flux);
  const std::uint64_t local_cells =
      static_cast<std::uint64_t>(patch.cells.x) *
      static_cast<std::uint64_t>(patch.cells.y) *
      static_cast<std::uint64_t>(patch.cells.z);
  const std::uint64_t local_faces =
      owned_face_count({{0, 0, 0}, patch.cells});
  passed &= expect(counters.cells == local_cells &&
                       counters.faces == local_faces &&
                       counters.logical_bytes_read ==
                           4U * local_faces * sizeof(double) &&
                       counters.logical_bytes_written ==
                           local_faces * sizeof(double),
                   rank,
                   "interior plus shell reconstruct every owned face exactly once");
  passed &= expect(compare_rank_interfaces(patch, const_flux, MPI_COMM_WORLD),
                   rank,
                   "duplicated rank-interface faces are bitwise canonical");

  const std::array<FieldView, 1U> writes{divergence.view};
  KernelInvocation divergence_call{
      {}, Span<const FieldView>{writes.data(), writes.size()},
      {{0, 0, 0}, patch.cells}, 0U, 0U, 1U, flux.revision, &counters};
  passed &= expect(static_cast<bool>(cartesian_provisional_face_divergence(
                       plan, const_flux, divergence_call)),
                   rank, "face divergence consumes the exact flux revision");
  passed &= expect(counters.cells == 2U * local_cells &&
                       counters.faces == local_faces + 6U * local_cells &&
                       counters.logical_bytes_read ==
                           (4U * local_faces + 6U * local_cells) *
                               sizeof(double) &&
                       counters.logical_bytes_written ==
                           (local_faces + local_cells) * sizeof(double),
                   rank,
                   "decomposed reconstruction and divergence counters are additive");

  long double local_integral = 0.0L;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    const double dz = geometry.z().widths().data[patch.begin.z + z];
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      const double dy = geometry.y().widths().data[patch.begin.y + y];
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const double dx = geometry.x().widths().data[patch.begin.x + x];
        local_integral +=
            static_cast<long double>(divergence.view.unchecked({x, y, z}, 0U)) *
            dx * dy * dz;
      }
    }
  }
  long double global_integral = 0.0L;
  passed &= MPI_Allreduce(&local_integral, &global_integral, 1,
                          MPI_LONG_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD) == MPI_SUCCESS;

  long double local_boundary = 0.0L;
  if (patch.begin.x == 0) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        local_boundary -= const_flux.x.unchecked({0, y, z});
      }
    }
  }
  if (patch.begin.x + patch.cells.x == geometry.global_cells().x) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y; ++y) {
        local_boundary +=
            const_flux.x.unchecked({patch.cells.x, y, z});
      }
    }
  }
  if (patch.begin.y == 0) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        local_boundary -= const_flux.y.unchecked({x, 0, z});
      }
    }
  }
  if (patch.begin.y + patch.cells.y == geometry.global_cells().y) {
    for (std::int32_t z = 0; z < patch.cells.z; ++z) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        local_boundary +=
            const_flux.y.unchecked({x, patch.cells.y, z});
      }
    }
  }
  if (patch.begin.z == 0) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        local_boundary -= const_flux.z.unchecked({x, y, 0});
      }
    }
  }
  if (patch.begin.z + patch.cells.z == geometry.global_cells().z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        local_boundary +=
            const_flux.z.unchecked({x, y, patch.cells.z});
      }
    }
  }
  long double global_boundary = 0.0L;
  passed &= MPI_Allreduce(&local_boundary, &global_boundary, 1,
                          MPI_LONG_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD) == MPI_SUCCESS;
  const long double scale =
      std::max({1.0L, std::abs(global_integral), std::abs(global_boundary)});
  if (std::abs(global_integral - global_boundary) >
      512.0L * std::numeric_limits<double>::epsilon() * scale) {
    std::cerr << "rank " << rank << " integrated="
              << static_cast<double>(global_integral) << " boundary="
              << static_cast<double>(global_boundary) << " delta="
              << static_cast<double>(global_integral - global_boundary)
              << '\n';
  }
  passed &= expect(std::abs(global_integral - global_boundary) <=
                       512.0L * std::numeric_limits<double>::epsilon() * scale,
                   rank,
                   "integrated divergence equals physical boundary mass rate");

  passed = all_true(passed, MPI_COMM_WORLD);
  if (rank == 0 && passed) {
    std::cout << "v0.4 Cartesian decomposition kernel tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
