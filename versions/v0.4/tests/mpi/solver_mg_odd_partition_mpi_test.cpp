// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {
using namespace hundun::v04;

constexpr Int3 kGlobal{29, 17, 9};
constexpr std::size_t kGuardDoubles = 16U;
constexpr double kGuard = -9.87654321012345e107;
constexpr FieldId kWorkspaceField = 30U;

bool expect(bool condition, int rank, std::string_view message) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << message << '\n';
  }
  return condition;
}

bool all_true(bool local) {
  const int value = local ? 1 : 0;
  int result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         result != 0;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kGlobal;
  mesh.minimum_spacing = {1.0e-12, 1.0e-12, 1.0e-12};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
  std::size_t raw_begin{};
  std::size_t raw_values{};

  bool guards_intact() const noexcept {
    return std::all_of(storage.begin(),
                       storage.begin() + static_cast<std::ptrdiff_t>(raw_begin),
                       [](double value) { return value == kGuard; }) &&
           std::all_of(storage.begin() +
                           static_cast<std::ptrdiff_t>(raw_begin + raw_values),
                       storage.end(),
                       [](double value) { return value == kGuard; });
  }
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t ghosts,
                      StorageIdentity identity) {
  OwnedField result;
  const std::size_t nx =
      static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny =
      static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz =
      static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.raw_begin = kGuardDoubles;
  result.raw_values = nx * ny * nz;
  result.storage.assign(result.raw_values + 2U * kGuardDoubles, kGuard);
  std::fill_n(result.storage.data() + result.raw_begin, result.raw_values,
              0.0);
  result.view.base = result.storage.data() + result.raw_begin + ghosts +
                     static_cast<std::size_t>(ghosts) * (nx + nx * ny);
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = 1U;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = result.raw_values;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = identity;
  result.view.revision_domain = 701U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  ConstFaceFieldView view{};
};

OwnedFace make_face(CartesianAxis axis, Int3 cells,
                    StorageIdentity identity) {
  OwnedFace result;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) {
    ++extents.x;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
  } else {
    ++extents.z;
  }
  const std::size_t stride_y = static_cast<std::size_t>(extents.x);
  const std::size_t stride_z =
      stride_y * static_cast<std::size_t>(extents.y);
  result.storage.assign(stride_z * static_cast<std::size_t>(extents.z), 1.0);
  result.view = {result.storage.data(), extents, stride_y, stride_z, axis,
                 identity, 702U};
  return result;
}

std::int32_t selected(Int3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

bool mapped_levels(const MgWorkspaceRequirements& requirements) noexcept {
  if (requirements.level_count < 3U) {
    return false;
  }
  for (std::size_t level = 0U; level + 1U < requirements.level_count;
       ++level) {
    const auto& fine = requirements.levels[level];
    const auto& coarse = requirements.levels[level + 1U];
    for (int axis = 0; axis < 3; ++axis) {
      const std::uint8_t mask = static_cast<std::uint8_t>(1U << axis);
      const std::int32_t begin = selected(fine.patch.begin, axis);
      const std::int32_t end = begin + selected(fine.patch.cells, axis);
      const std::int32_t expected_begin =
          (fine.coarsen_axis_mask & mask) != 0U ? (begin + 1) / 2 : begin;
      const std::int32_t expected_end =
          (fine.coarsen_axis_mask & mask) != 0U ? (end + 1) / 2 : end;
      const std::int32_t expected_global =
          (fine.coarsen_axis_mask & mask) != 0U
              ? (selected(fine.global_shape, axis) + 1) / 2
              : selected(fine.global_shape, axis);
      if (selected(coarse.patch.begin, axis) != expected_begin ||
          selected(coarse.patch.cells, axis) !=
              expected_end - expected_begin ||
          selected(coarse.global_shape, axis) != expected_global) {
        return false;
      }
    }
  }
  return true;
}

double local_squared(ConstFieldView field) noexcept {
  double sum = 0.0;
  for (std::int32_t k = 0; k < field.interior.z; ++k) {
    for (std::int32_t j = 0; j < field.interior.y; ++j) {
      for (std::int32_t i = 0; i < field.interior.x; ++i) {
        const double value = field.unchecked({i, j, k}, 0U);
        sum += value * value;
      }
    }
  }
  return sum;
}

double global_norm(ConstFieldView field) noexcept {
  const double local = local_squared(field);
  double global = 0.0;
  if (MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::sqrt(global);
}

bool finite(ConstFieldView field) noexcept {
  for (std::int32_t k = 0; k < field.interior.z; ++k) {
    for (std::int32_t j = 0; j < field.interior.y; ++j) {
      for (std::int32_t i = 0; i < field.interior.x; ++i) {
        if (!std::isfinite(field.unchecked({i, j, k}, 0U))) {
          return false;
        }
      }
    }
  }
  return true;
}

void mark_workspace_padding(const MgWorkspaceRequirements& requirements,
                            OwnedField& arena) noexcept {
  double* const raw = arena.storage.data() + arena.raw_begin;
  for (std::size_t level = 0U; level < requirements.level_count; ++level) {
    const MgWorkspaceLevelRequirements& selected_level =
        requirements.levels[level];
    const std::size_t valid_x =
        static_cast<std::size_t>(selected_level.patch.cells.x) +
        2U * requirements.ghost_width;
    const std::size_t padded_y =
        static_cast<std::size_t>(selected_level.patch.cells.y) +
        2U * requirements.ghost_width;
    const std::size_t padded_z =
        static_cast<std::size_t>(selected_level.patch.cells.z) +
        2U * requirements.ghost_width;
    for (std::size_t slot = 0U; slot < requirements.slots_per_level; ++slot) {
      double* const slot_base =
          raw + selected_level.offset_doubles +
          slot * selected_level.slot_stride_doubles;
      for (std::size_t z = 0U; z < padded_z; ++z) {
        for (std::size_t y = 0U; y < padded_y; ++y) {
          double* const row = slot_base + z * selected_level.stride_z +
                              y * selected_level.stride_y;
          std::fill(row + valid_x, row + selected_level.stride_y, kGuard);
        }
      }
    }
  }
}

bool workspace_padding_intact(
    const MgWorkspaceRequirements& requirements,
    const OwnedField& arena) noexcept {
  const double* const raw = arena.storage.data() + arena.raw_begin;
  for (std::size_t level = 0U; level < requirements.level_count; ++level) {
    const MgWorkspaceLevelRequirements& selected_level =
        requirements.levels[level];
    const std::size_t valid_x =
        static_cast<std::size_t>(selected_level.patch.cells.x) +
        2U * requirements.ghost_width;
    const std::size_t padded_y =
        static_cast<std::size_t>(selected_level.patch.cells.y) +
        2U * requirements.ghost_width;
    const std::size_t padded_z =
        static_cast<std::size_t>(selected_level.patch.cells.z) +
        2U * requirements.ghost_width;
    for (std::size_t slot = 0U; slot < requirements.slots_per_level; ++slot) {
      const double* const slot_base =
          raw + selected_level.offset_doubles +
          slot * selected_level.slot_stride_doubles;
      for (std::size_t z = 0U; z < padded_z; ++z) {
        for (std::size_t y = 0U; y < padded_y; ++y) {
          const double* const row =
              slot_base + z * selected_level.stride_z +
              y * selected_level.stride_y;
          if (!std::all_of(row + valid_x, row + selected_level.stride_y,
                           [](double value) { return value == kGuard; })) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

double neighbor(ConstFieldView field, const MeshPatch& patch, Int3 cell,
                int axis, int direction) noexcept {
  Int3 selected_cell = cell;
  std::int32_t* index = axis == 0
                            ? &selected_cell.x
                            : (axis == 1 ? &selected_cell.y
                                         : &selected_cell.z);
  const std::int32_t local_extent = selected(field.interior, axis);
  *index += direction;
  if (*index >= 0 && *index < local_extent) {
    return field.unchecked(selected_cell, 0U);
  }
  const std::int32_t begin = selected(patch.begin, axis);
  const std::int32_t global_extent = selected(kGlobal, axis);
  const bool physical = direction < 0 ? begin == 0
                                      : begin + local_extent == global_extent;
  return physical ? 0.0 : field.unchecked(selected_cell, 0U);
}

void independent_residual(ConstFieldView old_residual,
                          ConstFieldView correction, const MeshPatch& patch,
                          FieldView next_residual) noexcept {
  const Int3 cells = old_residual.interior;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        double applied = 6.0 * correction.unchecked(cell, 0U);
        applied -= neighbor(correction, patch, cell, 0, -1);
        applied -= neighbor(correction, patch, cell, 0, 1);
        applied -= neighbor(correction, patch, cell, 1, -1);
        applied -= neighbor(correction, patch, cell, 1, 1);
        applied -= neighbor(correction, patch, cell, 2, -1);
        applied -= neighbor(correction, patch, cell, 2, 1);
        next_residual.unchecked(cell, 0U) =
            old_residual.unchecked(cell, 0U) - applied;
      }
    }
  }
}

bool test_odd_partition_transfer(int rank) {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, mesh_spec(), GeometryBudget{}, geometry, patch)),
      rank, "odd global Cartesian geometry compiles");
  if (!all_true(passed)) {
    return false;
  }

  NativeCartesianMgSpec spec;
  spec.communicator = MPI_COMM_WORLD;
  spec.geometry = &geometry;
  spec.patch = patch;
  spec.identity = {101U, 102U, 103U, 104U, 105U};
  spec.coefficients = {201U, 202U, 0.0};
  spec.policy.pre_sweeps = 2U;
  spec.policy.post_sweeps = 2U;
  spec.policy.maximum_levels = 16U;
  spec.policy.coarse_sweeps = 32U;

  MgWorkspaceRequirements requirements{};
  passed &= expect(static_cast<bool>(make_mg_workspace_requirements(
                       MPI_COMM_WORLD, geometry, patch, spec.policy, 301U,
                       requirements)) &&
                       mapped_levels(requirements),
                   rank, "odd patches retain canonical ceil ownership on every level");

  bool local_odd_interface = false;
  for (int axis = 0; axis < 3; ++axis) {
    const std::int32_t begin = selected(patch.begin, axis);
    const std::int32_t end = begin + selected(patch.cells, axis);
    const std::int32_t extent = selected(kGlobal, axis);
    local_odd_interface |= (begin > 0 && (begin & 1) != 0) ||
                           (end < extent && (end & 1) != 0);
  }
  const int local_odd = local_odd_interface ? 1 : 0;
  int global_odd = 0;
  MPI_Allreduce(&local_odd, &global_odd, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(global_odd != 0, rank,
                   "decomposition exercises an odd internal fine boundary");
  if (!all_true(passed)) {
    return false;
  }

  OwnedField diagonal = make_field(1U, patch.cells, 0U, 401U);
  OwnedField residual = make_field(2U, patch.cells, 0U, 402U);
  OwnedField correction = make_field(3U, patch.cells, 1U, 403U);
  OwnedField next_residual = make_field(4U, patch.cells, 0U, 404U);
  OwnedField arena =
      make_field(kWorkspaceField, requirements.arena_shape, 0U, 405U);
  mark_workspace_padding(requirements, arena);
  OwnedFace x = make_face(CartesianAxis::x, patch.cells, 411U);
  OwnedFace y = make_face(CartesianAxis::y, patch.cells, 412U);
  OwnedFace z = make_face(CartesianAxis::z, patch.cells, 413U);
  for (std::int32_t k = 0; k < patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < patch.cells.x; ++i) {
        const Int3 global{patch.begin.x + i, patch.begin.y + j,
                          patch.begin.z + k};
        diagonal.view.unchecked({i, j, k}, 0U) = 6.0;
        residual.view.unchecked({i, j, k}, 0U) =
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(global.x) + 0.5) / kGlobal.x) *
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(global.y) + 0.5) / kGlobal.y) *
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(global.z) + 0.5) / kGlobal.z);
      }
    }
  }

  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine finest_halo;
  HaloEngine oracle_halo;
  const std::array<HaloFieldSpec, 1U> workspace_fields{{
      {kWorkspaceField, 1U, 1U}}};
  const std::array<HaloFieldSpec, 1U> oracle_fields{{{3U, 1U, 1U}}};
  passed &= expect(
      static_cast<bool>(MgWorkspace::bind(requirements, arena.view,
                                          workspace)) &&
          static_cast<bool>(ReductionEngine::compile(
              MPI_COMM_WORLD, ReductionMode::mpi_allreduce, 4U,
              reductions)) &&
          static_cast<bool>(finest_halo.reserve(
              MPI_COMM_WORLD, requirements.levels[0U].patch,
              {workspace_fields.data(), workspace_fields.size()})) &&
          static_cast<bool>(oracle_halo.reserve(
              MPI_COMM_WORLD, patch,
              {oracle_fields.data(), oracle_fields.size()})),
      rank, "persistent workspace, reduction, and oracle halos reserve");

  std::vector<HaloEngine> coarse_halos(requirements.level_count - 1U);
  std::vector<HaloEngine*> coarse_halo_pointers(coarse_halos.size());
  for (std::size_t level = 1U; level < requirements.level_count; ++level) {
    passed &= expect(static_cast<bool>(coarse_halos[level - 1U].reserve(
                         MPI_COMM_WORLD, requirements.levels[level].patch,
                         {workspace_fields.data(), workspace_fields.size()})),
                     rank, "coarse persistent halo reserves");
    coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
  }
  if (!all_true(passed)) {
    return false;
  }

  NativeCartesianMgPlan plan;
  const MgRuntimeServices services{
      &finest_halo, &reductions, &workspace,
      {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          spec, services,
          MgCoefficientViews{as_const(diagonal.view), x.view, y.view, z.view},
          plan)) &&
          plan.level_count() >= 3U,
      rank, "odd-partition multilevel MG plan compiles");
  if (!all_true(passed)) {
    return false;
  }

  const double initial = global_norm(as_const(residual.view));
  double previous = initial;
  for (std::uint32_t cycle = 0U; cycle < 3U; ++cycle) {
    const Status applied =
        plan.apply(as_const(residual.view), correction.view, cycle);
    HaloTicket ticket;
    std::array<FieldView, 1U> fields{correction.view};
    Status exchanged = oracle_halo.begin(
        static_cast<StageId>(1400U + cycle),
        {fields.data(), fields.size()}, ticket);
    if (exchanged) {
      exchanged = oracle_halo.finish(ticket, {fields.data(), fields.size()});
    }
    if (applied && exchanged) {
      independent_residual(as_const(residual.view), as_const(fields[0U]),
                           patch, next_residual.view);
    }
    const double current = global_norm(as_const(next_residual.view));
    if (!applied || !exchanged || !std::isfinite(current) ||
        current > previous * (1.0 + 1.0e-10)) {
      std::cerr << "rank " << rank << " odd cycle " << cycle
                << " residual " << previous << " -> " << current
                << " apply=" << static_cast<unsigned>(applied.code) << ':'
                << applied.detail << " halo="
                << static_cast<unsigned>(exchanged.code) << ':'
                << exchanged.detail << '\n';
    }
    passed &= expect(static_cast<bool>(applied) &&
                         static_cast<bool>(exchanged) &&
                         finite(as_const(correction.view)) &&
                         finite(as_const(next_residual.view)) &&
                         std::isfinite(current) &&
                         current <= previous * (1.0 + 1.0e-10),
                     rank,
                     "odd-boundary V-cycle has finite monotone independent residual");
    std::swap(residual, next_residual);
    previous = current;
  }
  passed &= expect(previous < 0.80 * initial, rank,
                   "three odd-boundary V-cycles materially reduce true residual");
  passed &= expect(arena.guards_intact() && diagonal.guards_intact() &&
                       residual.guards_intact() &&
                       next_residual.guards_intact() &&
                       correction.guards_intact() &&
                       workspace_padding_intact(requirements, arena),
                   rank,
                   "caller guards and every MG slot row padding remain exact");
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
  bool passed = expect(size == 2 || size == 4, rank,
                       "odd-partition RED runs at 2 or 4 ranks");
  passed &= test_odd_partition_transfer(rank);
  passed = all_true(passed);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
