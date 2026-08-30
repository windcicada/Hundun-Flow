// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"
#include "solver_mg_detail.hpp"

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

constexpr Int3 kCells{64, 12, 8};
constexpr FieldId kLineField = 30U;
constexpr FieldId kPointField = 31U;

bool expect(bool condition, std::string_view message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
  }
  return condition;
}

CartesianMeshSpec anisotropic_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {0.05, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kCells;
  mesh.minimum_spacing = {1.0e-12, 1.0e-12, 1.0e-12};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField field(FieldId id, Int3 cells, std::uint8_t ghosts,
                 StorageIdentity identity) {
  OwnedField result;
  const std::size_t nx =
      static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny =
      static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz =
      static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz, 0.0);
  result.view.base = result.storage.data() + ghosts +
                     static_cast<std::size_t>(ghosts) * (nx + nx * ny);
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = 1U;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = identity;
  result.view.revision_domain = 801U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  ConstFaceFieldView view{};
};

OwnedFace face(CartesianAxis axis, Int3 cells, double coefficient,
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
  result.storage.assign(stride_z * static_cast<std::size_t>(extents.z),
                        coefficient);
  result.view = {result.storage.data(), extents, stride_y, stride_z, axis,
                 identity, 802U};
  return result;
}

struct Runtime {
  MgWorkspaceRequirements requirements{};
  OwnedField arena;
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine finest;
  std::vector<HaloEngine> coarse;
  std::vector<HaloEngine*> coarse_pointers;

  bool create(const CartesianGeometryPlan& geometry, MeshPatch patch,
              const MgHierarchyPolicy& policy, FieldId workspace_field,
              RevisionToken execution_revision,
              bool force_point_smoother) {
    if (!make_mg_workspace_requirements(MPI_COMM_SELF, geometry, patch,
                                        policy, execution_revision,
                                        requirements)) {
      return false;
    }
    if (force_point_smoother) {
      for (std::size_t level = 0U; level < requirements.level_count; ++level) {
        requirements.levels[level].line_axis_mask = 0U;
      }
      requirements.fingerprint += 1U;
      requirements.collective_fingerprint += 1U;
    }
    arena = field(workspace_field, requirements.arena_shape, 0U,
                  900U + workspace_field);
    if (!MgWorkspace::bind(requirements, arena.view, workspace) ||
        !ReductionEngine::compile(MPI_COMM_SELF,
                                  ReductionMode::mpi_allreduce, 4U,
                                  reductions)) {
      return false;
    }
    const std::array<HaloFieldSpec, 1U> fields{{
        {workspace_field, 1U, 1U}}};
    if (!finest.reserve(MPI_COMM_SELF, requirements.levels[0U].patch,
                        {fields.data(), fields.size()})) {
      return false;
    }
    coarse.resize(requirements.level_count - 1U);
    coarse_pointers.resize(coarse.size());
    for (std::size_t level = 1U; level < requirements.level_count; ++level) {
      if (!coarse[level - 1U].reserve(
              MPI_COMM_SELF, requirements.levels[level].patch,
              {fields.data(), fields.size()})) {
        return false;
      }
      coarse_pointers[level - 1U] = &coarse[level - 1U];
    }
    return true;
  }

  MgRuntimeServices services() noexcept {
    return {&finest, &reductions, &workspace,
            {coarse_pointers.data(), coarse_pointers.size()}};
  }
};

double norm(ConstFieldView field) noexcept {
  double sum = 0.0;
  for (std::int32_t k = 0; k < field.interior.z; ++k) {
    for (std::int32_t j = 0; j < field.interior.y; ++j) {
      for (std::int32_t i = 0; i < field.interior.x; ++i) {
        const double value = field.unchecked({i, j, k}, 0U);
        sum += value * value;
      }
    }
  }
  return std::sqrt(sum);
}

void apply_operator(ConstFieldView input, FieldView output,
                    double ax, double ay, double az) noexcept {
  const Int3 cells = input.interior;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        double value = 2.0 * (ax + ay + az) * input.unchecked(cell, 0U);
        value -= ax * (i == 0 ? 0.0
                              : input.unchecked({i - 1, j, k}, 0U));
        value -= ax * (i + 1 == cells.x
                           ? 0.0
                           : input.unchecked({i + 1, j, k}, 0U));
        value -= ay * (j == 0 ? 0.0
                              : input.unchecked({i, j - 1, k}, 0U));
        value -= ay * (j + 1 == cells.y
                           ? 0.0
                           : input.unchecked({i, j + 1, k}, 0U));
        value -= az * (k == 0 ? 0.0
                              : input.unchecked({i, j, k - 1}, 0U));
        value -= az * (k + 1 == cells.z
                           ? 0.0
                           : input.unchecked({i, j, k + 1}, 0U));
        output.unchecked(cell, 0U) = value;
      }
    }
  }
}

double post_cycle_residual(const NativeCartesianMgPlan& plan,
                           ConstFieldView rhs, FieldView correction,
                           FieldView applied, double ax, double ay,
                           double az) noexcept {
  const Status status = const_cast<NativeCartesianMgPlan&>(plan).apply(
      rhs, correction, 0U);
  if (!status) {
    return std::numeric_limits<double>::infinity();
  }
  apply_operator(as_const(correction), applied, ax, ay, az);
  for (std::int32_t k = 0; k < rhs.interior.z; ++k) {
    for (std::int32_t j = 0; j < rhs.interior.y; ++j) {
      for (std::int32_t i = 0; i < rhs.interior.x; ++i) {
        applied.unchecked({i, j, k}, 0U) =
            rhs.unchecked({i, j, k}, 0U) -
            applied.unchecked({i, j, k}, 0U);
      }
    }
  }
  return norm(as_const(applied));
}

bool test_line_relaxation_effect() {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, anisotropic_mesh(), GeometryBudget{}, geometry,
          patch)),
      "uniform but strongly anisotropic geometry compiles");
  if (!passed) {
    return false;
  }

  const double dx = geometry.x().uniform_width();
  const double dy = geometry.y().uniform_width();
  const double dz = geometry.z().uniform_width();
  const double ax = 1.0 / (dx * dx);
  const double ay = 1.0 / (dy * dy);
  const double az = 1.0 / (dz * dz);
  OwnedField diagonal = field(1U, patch.cells, 0U, 101U);
  OwnedField rhs = field(2U, patch.cells, 0U, 102U);
  OwnedField line_correction = field(3U, patch.cells, 0U, 103U);
  OwnedField point_correction = field(4U, patch.cells, 0U, 104U);
  OwnedField applied = field(5U, patch.cells, 0U, 105U);
  for (std::int32_t k = 0; k < patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < patch.cells.x; ++i) {
        diagonal.view.unchecked({i, j, k}, 0U) =
            2.0 * (ax + ay + az);
        rhs.view.unchecked({i, j, k}, 0U) =
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(i) + 0.5) / kCells.x) *
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(j) + 0.5) / kCells.y) *
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(k) + 0.5) / kCells.z);
      }
    }
  }
  OwnedFace x = face(CartesianAxis::x, patch.cells, ax, 111U);
  OwnedFace y = face(CartesianAxis::y, patch.cells, ay, 112U);
  OwnedFace z = face(CartesianAxis::z, patch.cells, az, 113U);
  const MgCoefficientViews coefficients{as_const(diagonal.view), x.view,
                                        y.view, z.view};

  NativeCartesianMgSpec line_spec;
  line_spec.communicator = MPI_COMM_SELF;
  line_spec.geometry = &geometry;
  line_spec.patch = patch;
  line_spec.identity = {201U, 202U, 203U, 204U, 205U};
  line_spec.coefficients = {301U, 302U, 0.0};
  line_spec.policy.anisotropy_threshold = 4.0;
  line_spec.policy.pre_sweeps = 1U;
  line_spec.policy.post_sweeps = 1U;
  line_spec.policy.maximum_levels = 12U;
  line_spec.policy.coarse_sweeps = 16U;
  line_spec.policy.point_smoother = MgPointSmootherKind::chebyshev_jacobi;
  line_spec.operator_class =
      MgOperatorClass::symmetric_diagonally_dominant_m_matrix;

  Runtime line_runtime;
  Runtime point_runtime;
  passed &= expect(line_runtime.create(geometry, patch, line_spec.policy,
                                       kLineField, 401U, false) &&
                       point_runtime.create(geometry, patch, line_spec.policy,
                                            kPointField, 402U, true),
                   "line and point-only runtime workspaces compile");
  if (!passed) {
    return false;
  }
  NativeCartesianMgSpec point_spec = line_spec;
  point_spec.identity = {211U, 212U, 213U, 214U, 215U};
  NativeCartesianMgPlan line_plan;
  NativeCartesianMgPlan point_plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          line_spec, line_runtime.services(), coefficients, line_plan)) &&
          static_cast<bool>(NativeCartesianMgPlan::compile(
              point_spec, point_runtime.services(), coefficients,
              point_plan)) &&
          line_plan.line_axis_mask() == 1U &&
          point_plan.line_axis_mask() == 0U,
      "anisotropic plan selects x lines while control forces point smoothing");
  if (!passed) {
    return false;
  }

  const double initial = norm(as_const(rhs.view));
  const double line_residual = post_cycle_residual(
      line_plan, as_const(rhs.view), line_correction.view, applied.view,
      ax, ay, az);
  const double point_residual = post_cycle_residual(
      point_plan, as_const(rhs.view), point_correction.view, applied.view,
      ax, ay, az);
  const detail::MgMatrixWorkCounters line_work =
      detail::mg_matrix_work_counters_for_test(line_plan);
  const detail::MgMatrixWorkCounters point_work =
      detail::mg_matrix_work_counters_for_test(point_plan);
  passed &= expect(line_work.chebyshev_stages == 0U &&
                       point_work.chebyshev_stages > 0U &&
                       point_work.chebyshev_exchange_actions ==
                           point_work.chebyshev_stages +
                               point_work.chebyshev_retained_final_defect_actions &&
                       point_work.chebyshev_defect_actions ==
                           point_work.chebyshev_exchange_actions &&
                       point_work.chebyshev_retained_final_defect_actions > 0U,
                   "line levels retain red/black while point levels use certified Chebyshev");
  passed &= expect(std::isfinite(line_residual) &&
                       std::isfinite(point_residual) &&
                       line_residual <= 0.10 * initial,
                   "line smoother produces a finite strongly residual-reducing cycle");
  passed &= expect(line_residual <= 0.70 * point_residual,
                   "actual line relaxation materially outperforms point-only smoothing");
  if (!(std::isfinite(line_residual) && std::isfinite(point_residual) &&
        line_residual <= 0.10 * initial &&
        line_residual <= 0.70 * point_residual)) {
    std::cerr << "line/point residuals: " << line_residual << " / "
              << point_residual << " (initial " << initial << ")\n";
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_line_relaxation_effect();
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
