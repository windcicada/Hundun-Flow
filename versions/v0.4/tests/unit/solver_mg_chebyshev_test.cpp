// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_parallel.hpp"
#include "solver_mg_detail.hpp"

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

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool same_counters(const MgPlanCounters& left,
                   const MgPlanCounters& right) noexcept {
  return left.symbolic_builds == right.symbolic_builds &&
         left.numeric_refreshes == right.numeric_refreshes &&
         left.hierarchy_rebuilds == right.hierarchy_rebuilds &&
         left.applications == right.applications &&
         left.blocking_collectives == right.blocking_collectives &&
         left.collective_logical_bytes == right.collective_logical_bytes;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField field(Int3 cells, std::uint8_t ghosts, FieldId id,
                 StorageIdentity storage_identity) {
  OwnedField result;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz, 0.0);
  result.view.base = result.storage.data() + ghosts + ghosts * nx +
                     ghosts * nx * ny;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = 1U;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = storage_identity;
  result.view.revision_domain = 900U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  FaceFieldView view{};
};

OwnedFace face(Int3 extents, CartesianAxis axis,
               StorageIdentity storage_identity, double value = 1.0) {
  OwnedFace result;
  const std::size_t count = static_cast<std::size_t>(extents.x) *
                            static_cast<std::size_t>(extents.y) *
                            static_cast<std::size_t>(extents.z);
  result.storage.assign(count, value);
  result.view.base = result.storage.data();
  result.view.extents = extents;
  result.view.stride_y = static_cast<std::size_t>(extents.x);
  result.view.stride_z = result.view.stride_y *
                         static_cast<std::size_t>(extents.y);
  result.view.axis = axis;
  result.view.storage_identity = storage_identity;
  result.view.revision_domain = 901U;
  return result;
}

ConstFaceFieldView as_const_face(const FaceFieldView& view) noexcept {
  return {view.base, view.extents, view.stride_y, view.stride_z, view.axis,
          view.storage_identity, view.revision_domain};
}

CartesianMeshSpec mesh_spec(Int3 cells = {8, 8, 8}) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {1.0e-8, 1.0e-8, 1.0e-8};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 28U};
  return mesh;
}

struct Fixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  OwnedField diagonal;
  OwnedFace x_faces;
  OwnedFace y_faces;
  OwnedFace z_faces;
  OwnedField rhs;
  OwnedField correction;
  OwnedField vectors;
  std::vector<std::uint8_t> activity_cells;
  std::vector<std::uint8_t> activity_x_faces;
  std::vector<std::uint8_t> activity_y_faces;
  std::vector<std::uint8_t> activity_z_faces;
  MgWorkspaceRequirements workspace_requirements{};
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine finest_halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  NativeCartesianMgSpec spec;
  NativeCartesianMgPlan plan;

  MgCoefficientViews coefficients() const noexcept {
    return {as_const(diagonal.view), as_const_face(x_faces.view),
            as_const_face(y_faces.view), as_const_face(z_faces.view)};
  }

  MgRuntimeServices services() noexcept {
    return {&finest_halo, &reductions, &workspace,
            {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
  }

  bool initialize(MgPointSmootherKind smoother, MgOperatorClass operator_class,
                  double diagonal_value = 7.0, double face_value = 1.0,
                  std::uint8_t stages = 2U,
                  MgBoundarySet boundaries = {},
                  MgNullSpace null_space = MgNullSpace::none,
                  bool with_inactive_cell = false,
                  Int3 mesh_cells = {8, 8, 8},
                  bool variable_coefficients = false,
                  double anisotropy_threshold = 4.0,
                  MgCycleKind cycle = MgCycleKind::v_cycle,
                  std::uint8_t maximum_levels = 2U) {
    if (!CartesianGeometryCompiler::compile(
            MPI_COMM_SELF, mesh_spec(mesh_cells),
            GeometryBudget{}, geometry, patch)) {
      return false;
    }
    const Int3 cells = patch.cells;
    diagonal = field(cells, 0U, 10U, 100U);
    std::fill(diagonal.storage.begin(), diagonal.storage.end(),
              diagonal_value);
    x_faces = face({cells.x + 1, cells.y, cells.z}, CartesianAxis::x, 101U,
                   face_value);
    y_faces = face({cells.x, cells.y + 1, cells.z}, CartesianAxis::y, 102U,
                   face_value);
    z_faces = face({cells.x, cells.y, cells.z + 1}, CartesianAxis::z, 103U,
                   face_value);
    if (variable_coefficients) {
      for (std::int32_t k = 0; k < cells.z; ++k) {
        for (std::int32_t j = 0; j < cells.y; ++j) {
          for (std::int32_t i = 0; i < cells.x; ++i) {
            diagonal.view.unchecked({i, j, k}, 0U) =
                diagonal_value + 0.03 * static_cast<double>(i + 1) +
                0.02 * static_cast<double>(j + 1) +
                0.01 * static_cast<double>(k + 1);
          }
        }
      }
      for (std::int32_t k = 0; k < cells.z; ++k) {
        for (std::int32_t j = 0; j < cells.y; ++j) {
          for (std::int32_t i = 0; i <= cells.x; ++i) {
            x_faces.view.unchecked({i, j, k}) =
                face_value + 0.005 * static_cast<double>(i + 1) +
                0.003 * static_cast<double>(j + 1) +
                0.001 * static_cast<double>(k + 1);
          }
        }
      }
      for (std::int32_t k = 0; k < cells.z; ++k) {
        for (std::int32_t j = 0; j <= cells.y; ++j) {
          for (std::int32_t i = 0; i < cells.x; ++i) {
            y_faces.view.unchecked({i, j, k}) =
                face_value + 0.004 * static_cast<double>(i + 1) +
                0.002 * static_cast<double>(j + 1) +
                0.001 * static_cast<double>(k + 1);
          }
        }
      }
      for (std::int32_t k = 0; k <= cells.z; ++k) {
        for (std::int32_t j = 0; j < cells.y; ++j) {
          for (std::int32_t i = 0; i < cells.x; ++i) {
            z_faces.view.unchecked({i, j, k}) =
                face_value + 0.003 * static_cast<double>(i + 1) +
                0.002 * static_cast<double>(j + 1) +
                0.004 * static_cast<double>(k + 1);
          }
        }
      }
    }
    rhs = field(cells, 1U, 20U, 110U);
    correction = field(cells, 1U, 21U, 111U);

    spec = {};
    spec.communicator = MPI_COMM_SELF;
    spec.geometry = &geometry;
    spec.patch = patch;
    spec.boundaries = boundaries;
    spec.null_space = null_space;
    spec.operator_class = operator_class;
    spec.policy.anisotropy_threshold = anisotropy_threshold;
    spec.policy.pre_sweeps = stages;
    spec.policy.post_sweeps = stages;
    spec.policy.maximum_levels = maximum_levels;
    spec.policy.coarse_sweeps = 2U;
    spec.policy.point_smoother = smoother;
    spec.policy.cycle = cycle;
    spec.policy.chebyshev_lower_spectrum_fraction = 0.3;
    spec.identity = {301U, 302U, 303U, 304U, 305U};
    spec.coefficients = {1U, 401U, 0.0};
    if (with_inactive_cell) {
      const std::size_t cells_count = static_cast<std::size_t>(cells.x) *
                                      static_cast<std::size_t>(cells.y) *
                                      static_cast<std::size_t>(cells.z);
      const Int3 xe{cells.x + 1, cells.y, cells.z};
      const Int3 ye{cells.x, cells.y + 1, cells.z};
      const Int3 ze{cells.x, cells.y, cells.z + 1};
      const auto count = [](Int3 shape) noexcept {
        return static_cast<std::size_t>(shape.x) *
               static_cast<std::size_t>(shape.y) *
               static_cast<std::size_t>(shape.z);
      };
      activity_cells.assign(cells_count, 1U);
      activity_x_faces.assign(count(xe), 1U);
      activity_y_faces.assign(count(ye), 1U);
      activity_z_faces.assign(count(ze), 1U);
      activity_cells[0U] = 0U;
      activity_x_faces[0U] = 0U;
      activity_x_faces[1U] = 0U;
      activity_y_faces[0U] = 0U;
      activity_y_faces[static_cast<std::size_t>(cells.x)] = 0U;
      activity_z_faces[0U] = 0U;
      activity_z_faces[static_cast<std::size_t>(cells.x) *
                       static_cast<std::size_t>(cells.y)] = 0U;
      spec.activity = {
          {activity_cells.data(), activity_cells.size()},
          {activity_x_faces.data(), activity_x_faces.size()},
          {activity_y_faces.data(), activity_y_faces.size()},
          {activity_z_faces.data(), activity_z_faces.size()}, 601U, 602U};
    }

    // Smoother choice does not change layout.  Cycle kind is nevertheless a
    // workspace semantic identity and must match the compiled plan.
    MgHierarchyPolicy workspace_policy = spec.policy;
    workspace_policy.point_smoother = MgPointSmootherKind::red_black;
    const Status workspace_status = make_mg_workspace_requirements(
        MPI_COMM_SELF, geometry, patch, workspace_policy, 81U,
        workspace_requirements);
    if (!workspace_status) {
      return false;
    }
    vectors = field(workspace_requirements.arena_shape, 1U, 30U, 120U);
    const Status bind_status =
        MgWorkspace::bind(workspace_requirements, vectors.view, workspace);
    const Status reduction_status = ReductionEngine::compile(
        MPI_COMM_SELF, ReductionMode::mpi_allreduce, 4U, reductions);
    if (!bind_status || !reduction_status) {
      return false;
    }
    const std::array<HaloFieldSpec, 1U> fields{{{30U, 1U, 1U}}};
    const HaloTopology topology{
        boundaries.x_min == MgBoundaryKind::periodic,
        boundaries.y_min == MgBoundaryKind::periodic,
        boundaries.z_min == MgBoundaryKind::periodic};
    if (!finest_halo.reserve(MPI_COMM_SELF,
                             workspace_requirements.levels[0U].patch,
                             {fields.data(), fields.size()}, topology)) {
      return false;
    }
    coarse_halos.resize(workspace_requirements.level_count - 1U);
    coarse_halo_pointers.resize(coarse_halos.size());
    for (std::size_t level = 1U;
         level < workspace_requirements.level_count; ++level) {
      if (!coarse_halos[level - 1U].reserve(
              MPI_COMM_SELF, workspace_requirements.levels[level].patch,
              {fields.data(), fields.size()}, topology)) {
        return false;
      }
      coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
    }
    const Status compile_status = NativeCartesianMgPlan::compile(
        spec, services(), coefficients(), plan);
    return static_cast<bool>(compile_status);
  }

  void fill_rhs(bool checkerboard = false) noexcept {
    for (std::int32_t k = 0; k < patch.cells.z; ++k) {
      for (std::int32_t j = 0; j < patch.cells.y; ++j) {
        for (std::int32_t i = 0; i < patch.cells.x; ++i) {
          rhs.view.unchecked({i, j, k}, 0U) =
              checkerboard
                  ? (((i + j + k) & 1) == 0 ? 1.0 : -1.0)
                  : (std::sin(0.23 * static_cast<double>(i + 1)) +
                     0.7 * std::cos(0.31 * static_cast<double>(j + 1)) -
                     0.2 * static_cast<double>(k));
        }
      }
    }
  }
};

double apply_residual_norm(const Fixture& fixture,
                           ConstFieldView solution) noexcept {
  double sum = 0.0;
  for (std::int32_t k = 0; k < fixture.patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < fixture.patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < fixture.patch.cells.x; ++i) {
        const Int3 cell{i, j, k};
        double value = fixture.diagonal.view.unchecked(cell, 0U) *
                       solution.unchecked(cell, 0U);
        const auto neighbor = [&](CartesianAxis axis, int direction) {
          Int3 selected = cell;
          std::int32_t* coordinate = axis == CartesianAxis::x
                                          ? &selected.x
                                          : (axis == CartesianAxis::y
                                                 ? &selected.y
                                                 : &selected.z);
          *coordinate += direction;
          const std::int32_t extent =
              axis == CartesianAxis::x
                  ? fixture.patch.cells.x
                  : (axis == CartesianAxis::y ? fixture.patch.cells.y
                                               : fixture.patch.cells.z);
          return *coordinate >= 0 && *coordinate < extent
                     ? solution.unchecked(selected, 0U)
                     : 0.0;
        };
        value -= fixture.x_faces.view.unchecked({i, j, k}) *
                 neighbor(CartesianAxis::x, -1);
        value -= fixture.x_faces.view.unchecked({i + 1, j, k}) *
                 neighbor(CartesianAxis::x, 1);
        value -= fixture.y_faces.view.unchecked({i, j, k}) *
                 neighbor(CartesianAxis::y, -1);
        value -= fixture.y_faces.view.unchecked({i, j + 1, k}) *
                 neighbor(CartesianAxis::y, 1);
        value -= fixture.z_faces.view.unchecked({i, j, k}) *
                 neighbor(CartesianAxis::z, -1);
        value -= fixture.z_faces.view.unchecked({i, j, k + 1}) *
                 neighbor(CartesianAxis::z, 1);
        const double residual = fixture.rhs.view.unchecked(cell, 0U) - value;
        sum += residual * residual;
      }
    }
  }
  return std::sqrt(sum);
}

std::vector<double> dense_action(const std::vector<double>& x,
                                 double diagonal) {
  std::vector<double> result(x.size(), 0.0);
  for (std::size_t i = 0U; i < x.size(); ++i) {
    result[i] = diagonal * x[i];
    if (i != 0U) result[i] -= x[i - 1U];
    if (i + 1U != x.size()) result[i] -= x[i + 1U];
  }
  return result;
}

void dense_red_black_smooth(std::vector<double>& x,
                            const std::vector<double>& rhs,
                            double diagonal, std::uint32_t sweeps,
                            bool reverse) {
  constexpr double omega = 0.72;
  for (std::uint32_t sweep = 0U; sweep < sweeps; ++sweep) {
    const int first = reverse ? 1 : 0;
    for (int pass = 0; pass < 2; ++pass) {
      const int color = (first + pass) & 1;
      for (std::size_t order = 0U; order < x.size(); ++order) {
        const std::size_t i = reverse ? x.size() - 1U - order : order;
        if (static_cast<int>(i & 1U) != color) continue;
        double action = diagonal * x[i];
        if (i != 0U) action -= x[i - 1U];
        if (i + 1U != x.size()) action -= x[i + 1U];
        x[i] += omega * (rhs[i] - action) / diagonal;
      }
    }
  }
}

std::vector<double> independent_two_level_cycle(
    const std::vector<double>& rhs, bool second_coarse_call) {
  std::vector<double> fine(rhs.size(), 0.0);
  dense_red_black_smooth(fine, rhs, 7.0, 2U, false);
  const std::vector<double> action = dense_action(fine, 7.0);
  std::vector<double> coarse_rhs(rhs.size() / 2U, 0.0);
  for (std::size_t coarse = 0U; coarse < coarse_rhs.size(); ++coarse) {
    coarse_rhs[coarse] = rhs[2U * coarse] - action[2U * coarse] +
                         rhs[2U * coarse + 1U] - action[2U * coarse + 1U];
  }
  // Summing two fine finite-volume rows gives reaction 2, x faces 1/1,
  // and the two y/z boundary-face pairs 2/2, hence diagonal 12.
  std::vector<double> coarse(coarse_rhs.size(), 0.0);
  dense_red_black_smooth(coarse, coarse_rhs, 12.0, 2U, false);
  if (second_coarse_call) {
    dense_red_black_smooth(coarse, coarse_rhs, 12.0, 2U, false);
  }
  for (std::size_t i = 0U; i < fine.size(); ++i) {
    const std::size_t parent = i / 2U;
    double correction = 0.0;
    if ((i & 1U) == 0U) {
      correction = parent == 0U
                       ? coarse[0U]
                       : 0.25 * coarse[parent - 1U] + 0.75 * coarse[parent];
    } else {
      correction = parent + 1U == coarse.size()
                       ? coarse.back()
                       : 0.75 * coarse[parent] + 0.25 * coarse[parent + 1U];
    }
    fine[i] += correction;
  }
  dense_red_black_smooth(fine, rhs, 7.0, 2U, true);
  // Native-MG publishes the scalar-minimal preconditioner image alpha*z,
  // where alpha minimizes ||r-alpha*A*z||.  This projection is part of the
  // public apply map, so the independent composition includes it explicitly.
  const std::vector<double> fine_action = dense_action(fine, 7.0);
  double numerator = 0.0;
  double denominator = 0.0;
  double residual_squared = 0.0;
  for (std::size_t i = 0U; i < fine.size(); ++i) {
    residual_squared += rhs[i] * rhs[i];
    numerator += rhs[i] * fine_action[i];
    denominator += fine_action[i] * fine_action[i];
  }
  const double alpha = denominator >
                               std::numeric_limits<double>::epsilon() *
                                   std::max(1.0, residual_squared)
                           ? numerator / denominator
                           : 0.0;
  for (double& value : fine) value *= alpha;
  return fine;
}

std::uint64_t packed(Status status) noexcept {
  return (static_cast<std::uint64_t>(status.code) << 32U) |
         static_cast<std::uint64_t>(status.detail);
}

bool bitwise_same(const OwnedField& left, const OwnedField& right) noexcept {
  return left.storage.size() == right.storage.size() &&
         std::memcmp(left.storage.data(), right.storage.data(),
                     left.storage.size() * sizeof(double)) == 0;
}

bool bitwise_same(double left, double right) noexcept {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

bool test_public_contract() {
  MgHierarchyPolicy policy{};
  NativeCartesianMgSpec spec{};
  return expect(policy.point_smoother == MgPointSmootherKind::red_black &&
                    policy.cycle == MgCycleKind::v_cycle &&
                    policy.chebyshev_lower_spectrum_fraction == 0.3 &&
                    spec.operator_class == MgOperatorClass::general,
                "public defaults retain V-cycle, red/black, and general operator");
}

bool test_certified_apply_and_hot_schedule() {
  Fixture fixture;
  bool passed = expect(
      fixture.initialize(
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix),
      "certified Chebyshev fixture compiles from default workspace shape");
  if (!passed) return false;
  fixture.fill_rhs();
  const double initial =
      apply_residual_norm(fixture, as_const(fixture.correction.view));
  const Status status = fixture.plan.apply(as_const(fixture.rhs.view),
                                           fixture.correction.view, 0U);
  const detail::MgMatrixWorkCounters work =
      detail::mg_matrix_work_counters_for_test(fixture.plan);
  const double final_residual =
      apply_residual_norm(fixture, as_const(fixture.correction.view));
  passed &= expect(static_cast<bool>(status) && std::isfinite(final_residual) &&
                       final_residual < initial,
                   "certified Chebyshev V-cycle returns finite residual reduction");
  const std::size_t cells = static_cast<std::size_t>(fixture.patch.cells.x) *
                            static_cast<std::size_t>(fixture.patch.cells.y) *
                            static_cast<std::size_t>(fixture.patch.cells.z);
  passed &= expect(
      fixture.plan.line_axis_mask() == 0U && work.chebyshev_stages == 4U &&
          work.chebyshev_exchange_actions == 5U &&
          work.chebyshev_defect_actions == 5U &&
          work.chebyshev_retained_final_defect_actions == 1U &&
          work.chebyshev_updates == 4U * cells,
      "compiled point route records exact degree-two pre/post hot schedule");
  return passed;
}

bool test_fixed_f_cycle_schedule_and_identity() {
  Fixture v_fixture;
  Fixture f_fixture;
  bool passed = expect(
      v_fixture.initialize(
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 7.0, 1.0,
          2U, {}, MgNullSpace::none, false, {8, 8, 8}, false, 4.0,
          MgCycleKind::v_cycle, 32U) &&
          f_fixture.initialize(
              MgPointSmootherKind::chebyshev_jacobi,
              MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 7.0,
              1.0, 2U, {}, MgNullSpace::none, false, {8, 8, 8}, false,
              4.0, MgCycleKind::f_cycle, 32U),
      "V/F identity fixtures compile with three levels");
  if (!passed) return false;
  passed &= expect(
      v_fixture.workspace_requirements.level_count == 3U &&
          f_fixture.workspace_requirements.level_count == 3U &&
          v_fixture.workspace_requirements.fingerprint !=
          f_fixture.workspace_requirements.fingerprint &&
          v_fixture.plan.symbolic_fingerprint() !=
              f_fixture.plan.symbolic_fingerprint() &&
          v_fixture.plan.certificate().collective_fingerprint !=
              f_fixture.plan.certificate().collective_fingerprint,
      "cycle kind mutates workspace, public symbolic, and prepared preconditioner identities");
  f_fixture.fill_rhs();
  const Status status = f_fixture.plan.apply(as_const(f_fixture.rhs.view),
                                             f_fixture.correction.view, 0U);
  const detail::MgMatrixWorkCounters work =
      detail::mg_matrix_work_counters_for_test(f_fixture.plan);
  passed &= expect(
      static_cast<bool>(status) && work.cycle_level_calls[0U] == 1U &&
          work.cycle_level_calls[1U] == 2U &&
          work.cycle_level_calls[2U] == 3U &&
          work.cycle_restrictions == 3U &&
          work.cycle_prolongations == 3U,
      "fixed kappa-two cycle executes exact 1/2/3 level-call schedule");

  NativeCartesianMgSpec invalid_spec = v_fixture.spec;
  invalid_spec.policy.cycle = static_cast<MgCycleKind>(255U);
  NativeCartesianMgPlan invalid_plan;
  const Status invalid_compile = NativeCartesianMgPlan::compile(
      invalid_spec, v_fixture.services(), v_fixture.coefficients(),
      invalid_plan);
  MgWorkspaceRequirements invalid_requirements;
  const Status invalid_workspace = make_mg_workspace_requirements(
      MPI_COMM_SELF, v_fixture.geometry, v_fixture.patch,
      invalid_spec.policy, 82U, invalid_requirements);
  passed &= expect(invalid_compile.code == StatusCode::invalid_plan &&
                       invalid_workspace.code == StatusCode::invalid_plan,
                   "unknown cycle enum is rejected by workspace and plan");

  bool strict_improvement = false;
  bool all_no_worse = true;
  for (std::size_t levels = 2U; levels <= 7U; ++levels) {
    const std::int32_t extent =
        static_cast<std::int32_t>(3U << (levels - 1U));
    Fixture v_schedule;
    Fixture f_schedule;
    passed &= expect(
        v_schedule.initialize(
            MgPointSmootherKind::chebyshev_jacobi,
            MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 7.0,
            1.0, 2U, {}, MgNullSpace::none, false, {extent, 1, 1}, false,
            1.0e9, MgCycleKind::v_cycle, 32U) &&
            f_schedule.initialize(
                MgPointSmootherKind::chebyshev_jacobi,
                MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 7.0,
                1.0, 2U, {}, MgNullSpace::none, false, {extent, 1, 1}, false,
                1.0e9, MgCycleKind::f_cycle, 32U) &&
            v_schedule.plan.level_count() == levels &&
            f_schedule.plan.level_count() == levels,
        "generic V/F fixtures compile for two through seven levels");
    if (!passed) break;
    v_schedule.fill_rhs();
    f_schedule.fill_rhs();
    const bool same_seed_and_rhs =
        bitwise_same(v_schedule.rhs, f_schedule.rhs) &&
        bitwise_same(v_schedule.correction, f_schedule.correction);
    const Status v_status = v_schedule.plan.apply(
        as_const(v_schedule.rhs.view), v_schedule.correction.view, 0U);
    const Status f_status = f_schedule.plan.apply(
        as_const(f_schedule.rhs.view), f_schedule.correction.view, 0U);
    const detail::MgMatrixWorkCounters v_work =
        detail::mg_matrix_work_counters_for_test(v_schedule.plan);
    const detail::MgMatrixWorkCounters f_work =
        detail::mg_matrix_work_counters_for_test(f_schedule.plan);
    bool calls_exact = true;
    for (std::size_t level = 0U; level < detail::kMgMaximumLevels; ++level) {
      const std::uint64_t expected_v = level < levels ? 1U : 0U;
      const std::uint64_t expected_f =
          level < levels ? static_cast<std::uint64_t>(level + 1U) : 0U;
      calls_exact &= v_work.cycle_level_calls[level] == expected_v &&
                     f_work.cycle_level_calls[level] == expected_f;
    }
    const std::uint64_t v_transfers = levels - 1U;
    const std::uint64_t f_transfers = levels * (levels - 1U) / 2U;
    const double v_residual = apply_residual_norm(
        v_schedule, as_const(v_schedule.correction.view));
    const double f_residual = apply_residual_norm(
        f_schedule, as_const(f_schedule.correction.view));
    if (levels >= 3U) {
      strict_improvement |= f_residual < v_residual;
      all_no_worse &= f_residual <= v_residual;
    }
    const bool route_valid =
        static_cast<bool>(v_status) && static_cast<bool>(f_status) &&
            calls_exact && v_work.cycle_restrictions == v_transfers &&
            v_work.cycle_prolongations == v_transfers &&
            f_work.cycle_restrictions == f_transfers &&
            f_work.cycle_prolongations == f_transfers &&
            same_seed_and_rhs && std::isfinite(v_residual) &&
            std::isfinite(f_residual);
    passed &= expect(
        route_valid,
        "generic V=ones/F=1..n schedule is exact with finite outputs");
  }
  passed &= expect(all_no_worse && strict_improvement,
                   "fixed F-cycle is no worse and strictly improves an independent true-residual fixture");
  return passed;
}

bool test_independent_fixed_cycle_composition() {
  Fixture v_fixture;
  Fixture f_fixture;
  bool passed = expect(
      v_fixture.initialize(MgPointSmootherKind::red_black,
                           MgOperatorClass::general, 7.0, 1.0, 2U, {},
                           MgNullSpace::none, false, {6, 1, 1}, false, 1.0e9,
                           MgCycleKind::v_cycle, 2U) &&
          f_fixture.initialize(MgPointSmootherKind::red_black,
                               MgOperatorClass::general, 7.0, 1.0, 2U, {},
                               MgNullSpace::none, false, {6, 1, 1}, false,
                               1.0e9, MgCycleKind::f_cycle, 2U),
      "independent dense-composition fixtures compile");
  if (!passed) return false;
  v_fixture.fill_rhs();
  f_fixture.fill_rhs();
  std::vector<double> rhs(6U, 0.0);
  for (std::size_t i = 0U; i < rhs.size(); ++i) {
    rhs[i] = v_fixture.rhs.view.unchecked(
        {static_cast<std::int32_t>(i), 0, 0}, 0U);
  }
  const std::vector<double> expected_v =
      independent_two_level_cycle(rhs, false);
  const std::vector<double> expected_f =
      independent_two_level_cycle(rhs, true);
  const Status v_status = v_fixture.plan.apply(
      as_const(v_fixture.rhs.view), v_fixture.correction.view, 0U);
  const Status f_status = f_fixture.plan.apply(
      as_const(f_fixture.rhs.view), f_fixture.correction.view, 0U);
  double maximum_v_error = 0.0;
  double maximum_f_error = 0.0;
  double maximum_vf_difference = 0.0;
  for (std::size_t i = 0U; i < rhs.size(); ++i) {
    const Int3 cell{static_cast<std::int32_t>(i), 0, 0};
    const double actual_v = v_fixture.correction.view.unchecked(cell, 0U);
    const double actual_f = f_fixture.correction.view.unchecked(cell, 0U);
    maximum_v_error =
        std::max(maximum_v_error, std::abs(actual_v - expected_v[i]));
    maximum_f_error =
        std::max(maximum_f_error, std::abs(actual_f - expected_f[i]));
    maximum_vf_difference =
        std::max(maximum_vf_difference, std::abs(actual_v - actual_f));
  }
  const double tolerance = 32.0 * std::numeric_limits<double>::epsilon();
  passed &= expect(
      static_cast<bool>(v_status) && static_cast<bool>(f_status) &&
          detail::mg_level_diagonal_for_test(v_fixture.plan, 1U, 0U) ==
              12.0 &&
          maximum_v_error <= tolerance && maximum_f_error <= tolerance &&
          maximum_vf_difference > tolerance,
      "production V/F agree with independent pre-R-two-coarse-P-post composition and differ when the second call is active");
  return passed;
}

bool test_certificate_rejection_is_transactional() {
  Fixture fixture;
  if (!expect(
          fixture.initialize(
              MgPointSmootherKind::chebyshev_jacobi,
              MgOperatorClass::symmetric_diagonally_dominant_m_matrix),
          "transaction fixture compiles")) {
    return false;
  }
  fixture.fill_rhs();
  const Status baseline_apply = fixture.plan.apply(
      as_const(fixture.rhs.view), fixture.correction.view, 0U);
  if (!expect(static_cast<bool>(baseline_apply),
              "transaction fixture establishes an active numerical state")) {
    return false;
  }
  const std::vector<double> correction_before = fixture.correction.storage;
  const PlanFingerprint numeric = fixture.plan.numeric_fingerprint();
  const PlanFingerprint hierarchy = fixture.plan.hierarchy_fingerprint();
  const RevisionToken generation = fixture.plan.generation();
  const MgPlanCounters counters = fixture.plan.counters();
  const double initial = fixture.plan.last_cycle_initial_residual();
  const double final = fixture.plan.last_cycle_final_residual();
  const std::uintptr_t hierarchy_storage =
      fixture.plan.hierarchy_storage_address();
  const std::uintptr_t workspace_storage =
      fixture.plan.workspace_storage_address();
  const std::size_t workspace_doubles = fixture.plan.workspace_doubles();
  const std::size_t level_count = fixture.plan.level_count();
  const std::uint8_t line_axis_mask = fixture.plan.line_axis_mask();
  std::fill(fixture.diagonal.storage.begin(), fixture.diagonal.storage.end(),
            5.0);
  LinearIdentity next = fixture.spec.identity;
  next.numeric = 902U;
  next.hierarchy = 903U;
  next.fingerprint = 904U;
  const Status rejected = fixture.plan.update_coefficients(
      next, MgCoefficientIdentity{2U, 402U, 0.5}, fixture.coefficients());
  bool passed = expect(
      rejected.code == StatusCode::numerical_failure &&
          fixture.plan.numeric_fingerprint() == numeric &&
          fixture.plan.hierarchy_fingerprint() == hierarchy &&
          fixture.plan.generation() == generation &&
          same_counters(fixture.plan.counters(), counters) &&
          fixture.plan.last_cycle_initial_residual() == initial &&
          fixture.plan.last_cycle_final_residual() == final &&
          fixture.plan.hierarchy_storage_address() == hierarchy_storage &&
          fixture.plan.workspace_storage_address() == workspace_storage &&
          fixture.plan.workspace_doubles() == workspace_doubles &&
          fixture.plan.level_count() == level_count &&
          fixture.plan.line_axis_mask() == line_axis_mask,
      "non-dominant certified refresh rejects with complete active state and counter preservation");
  if (!passed) return false;
  const Status after_failure_apply = fixture.plan.apply(
      as_const(fixture.rhs.view), fixture.correction.view, 1U);
  bool output_unchanged = static_cast<bool>(after_failure_apply);
  for (std::size_t index = 0U; index < fixture.correction.storage.size();
       ++index) {
    output_unchanged = output_unchanged &&
                       fixture.correction.storage[index] == correction_before[index];
  }
  return expect(output_unchanged,
                "failed certified refresh leaves the active hierarchy numerically unchanged");
}

double scalar_chebyshev_step(double rhs, double diagonal,
                             double operator_eigenvalue,
                             std::uint8_t stages,
                             double lower_fraction) noexcept {
  const double beta = 2.0;
  const double alpha = lower_fraction * beta;
  const double theta = (beta + alpha) * 0.5;
  const double delta = (beta - alpha) * 0.5;
  const double sigma = theta / delta;
  double solution = 0.0;
  double direction = 0.0;
  double previous_rho = 1.0 / sigma;
  for (std::uint8_t stage = 0U; stage < stages; ++stage) {
    const double defect = rhs - operator_eigenvalue * solution;
    if (stage == 0U) {
      direction = defect / (theta * diagonal);
    } else {
      const double rho = 1.0 / (2.0 * sigma - previous_rho);
      direction = rho * previous_rho * direction +
                  (2.0 * rho / delta) * defect / diagonal;
      previous_rho = rho;
    }
    solution += direction;
  }
  return solution;
}

double scalar_two_pass_chebyshev(double rhs, double diagonal,
                                 double operator_eigenvalue,
                                 std::uint8_t stages,
                                 double lower_fraction) noexcept {
  const double pre = scalar_chebyshev_step(
      rhs, diagonal, operator_eigenvalue, stages, lower_fraction);
  const double post = scalar_chebyshev_step(
      rhs - operator_eigenvalue * pre, diagonal, operator_eigenvalue, stages,
      lower_fraction);
  return pre + post;
}

bool test_degree_one_through_four_oracle() {
  bool passed = true;
  for (std::uint8_t degree = 1U; degree <= 4U; ++degree) {
    Fixture fixture;
    passed &= expect(
        fixture.initialize(
            MgPointSmootherKind::chebyshev_jacobi,
            MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 2.0,
            0.0, degree),
        "diagonal degree oracle fixture compiles");
    if (!passed) continue;
    fixture.fill_rhs(true);
    detail::set_mg_skip_final_projection_for_test(fixture.plan, true);
    const Status status = fixture.plan.apply(as_const(fixture.rhs.view),
                                             fixture.correction.view, 0U);
    const detail::MgMatrixWorkCounters work =
        detail::mg_matrix_work_counters_for_test(fixture.plan);
    const double expected = scalar_two_pass_chebyshev(
        1.0, 2.0, 2.0, degree,
        fixture.spec.policy.chebyshev_lower_spectrum_fraction);
    bool values_match = true;
    for (std::int32_t k = 0; k < fixture.patch.cells.z; ++k) {
      for (std::int32_t j = 0; j < fixture.patch.cells.y; ++j) {
        for (std::int32_t i = 0; i < fixture.patch.cells.x; ++i) {
          const double value =
              fixture.correction.view.unchecked({i, j, k}, 0U);
          const double sign = ((i + j + k) & 1) == 0 ? 1.0 : -1.0;
          values_match =
              values_match &&
              std::abs(value - sign * expected) < 2.0e-12;
        }
      }
    }
    const std::size_t cells = static_cast<std::size_t>(fixture.patch.cells.x) *
                              static_cast<std::size_t>(fixture.patch.cells.y) *
                              static_cast<std::size_t>(fixture.patch.cells.z);
    passed &= expect(static_cast<bool>(status) && values_match &&
                         work.chebyshev_stages == 2U * degree &&
                         work.chebyshev_exchange_actions == 2U * degree + 1U &&
                         work.chebyshev_defect_actions == 2U * degree + 1U &&
                         work.chebyshev_updates == 2U * degree * cells,
                     "independent degree-one-through-four recurrence oracle");
  }
  return passed;
}

bool test_streamed_chebyshev_matches_legacy_degrees_and_shapes() {
  struct Case {
    Int3 cells;
    MgBoundarySet boundaries;
    MgNullSpace null_space;
    double diagonal;
    double face;
    bool variable;
    double anisotropy_threshold;
  };
  const std::array<Case, 5U> cases{{
      {{8, 8, 8}, {}, MgNullSpace::none, 7.0, 1.0, false, 4.0},
      {{6, 4, 5},
       {MgBoundaryKind::periodic, MgBoundaryKind::periodic,
        MgBoundaryKind::neumann, MgBoundaryKind::dirichlet,
        MgBoundaryKind::dirichlet, MgBoundaryKind::neumann},
       MgNullSpace::none, 8.0, 0.3, true, 4.0},
      {{1, 6, 3}, {}, MgNullSpace::none, 4.0, 0.2, true, 100.0},
      {{2, 6, 4},
       {MgBoundaryKind::neumann, MgBoundaryKind::neumann,
        MgBoundaryKind::neumann, MgBoundaryKind::neumann,
        MgBoundaryKind::neumann, MgBoundaryKind::neumann},
       MgNullSpace::none, 8.0, 0.4, true, 4.0},
      {{3, 6, 6},
       {MgBoundaryKind::periodic, MgBoundaryKind::periodic,
        MgBoundaryKind::periodic, MgBoundaryKind::periodic,
        MgBoundaryKind::periodic, MgBoundaryKind::periodic},
       MgNullSpace::constant, 6.0, 1.0, false, 4.0},
  }};

  bool passed = true;
  for (const Case& selected : cases) {
    for (std::uint8_t degree = 1U; degree <= 4U; ++degree) {
      Fixture streamed;
      Fixture legacy;
      const bool fixtures_compiled =
          streamed.initialize(
              MgPointSmootherKind::chebyshev_jacobi,
              MgOperatorClass::symmetric_diagonally_dominant_m_matrix,
              selected.diagonal, selected.face, degree, selected.boundaries,
              selected.null_space, false, selected.cells, selected.variable,
              selected.anisotropy_threshold) &&
          legacy.initialize(
              MgPointSmootherKind::chebyshev_jacobi,
              MgOperatorClass::symmetric_diagonally_dominant_m_matrix,
              selected.diagonal, selected.face, degree,
              selected.boundaries, selected.null_space, false,
              selected.cells, selected.variable,
              selected.anisotropy_threshold);
      passed &= expect(fixtures_compiled,
                       "streamed/legacy Chebyshev fixtures compile");
      if (!fixtures_compiled) continue;
      streamed.fill_rhs();
      legacy.fill_rhs();
      detail::set_mg_reference_chebyshev_lifecycle_for_test(streamed.plan,
                                                             false);
      detail::set_mg_reference_chebyshev_lifecycle_for_test(legacy.plan, true);
      const Status streamed_status = streamed.plan.apply(
          as_const(streamed.rhs.view), streamed.correction.view, 0U);
      const Status legacy_status = legacy.plan.apply(
          as_const(legacy.rhs.view), legacy.correction.view, 0U);
      const detail::MgMatrixWorkCounters streamed_work =
          detail::mg_matrix_work_counters_for_test(streamed.plan);
      const detail::MgMatrixWorkCounters legacy_work =
          detail::mg_matrix_work_counters_for_test(legacy.plan);
      const std::uint64_t cells =
          static_cast<std::uint64_t>(selected.cells.x) *
          static_cast<std::uint64_t>(selected.cells.y) *
          static_cast<std::uint64_t>(selected.cells.z);
      const std::uint64_t stages = 2U * degree;
      const std::uint64_t defects = stages + 1U;
      const bool common_work =
          streamed_work.full_actions == legacy_work.full_actions &&
          streamed_work.full_cell_visits == legacy_work.full_cell_visits &&
          streamed_work.residual_writes == legacy_work.residual_writes &&
          streamed_work.retained_final_full_actions ==
              legacy_work.retained_final_full_actions &&
          streamed_work.retained_final_defect_actions ==
              legacy_work.retained_final_defect_actions &&
          streamed_work.defect_cell_visits == legacy_work.defect_cell_visits &&
          streamed_work.residual_finish_actions ==
              legacy_work.residual_finish_actions &&
          streamed_work.residual_finish_cell_visits ==
              legacy_work.residual_finish_cell_visits &&
          streamed_work.fused_color_actions == legacy_work.fused_color_actions &&
          streamed_work.fused_cell_visits == legacy_work.fused_cell_visits &&
          streamed_work.fused_updates == legacy_work.fused_updates &&
          streamed_work.separated_color_actions ==
              legacy_work.separated_color_actions &&
          streamed_work.separated_color_updates ==
              legacy_work.separated_color_updates &&
          streamed_work.chebyshev_stages == stages &&
          streamed_work.chebyshev_stages == legacy_work.chebyshev_stages &&
          streamed_work.chebyshev_exchange_actions == stages + 1U &&
          streamed_work.chebyshev_exchange_actions ==
              legacy_work.chebyshev_exchange_actions &&
          streamed_work.chebyshev_defect_actions == defects &&
          streamed_work.chebyshev_defect_actions ==
              legacy_work.chebyshev_defect_actions &&
          streamed_work.chebyshev_retained_final_defect_actions == 1U &&
          streamed_work.chebyshev_retained_final_defect_actions ==
              legacy_work.chebyshev_retained_final_defect_actions &&
          streamed_work.chebyshev_stencil_evaluations == defects * cells &&
          streamed_work.chebyshev_stencil_evaluations ==
              legacy_work.chebyshev_stencil_evaluations &&
          streamed_work.chebyshev_updates == stages * cells &&
          streamed_work.chebyshev_updates == legacy_work.chebyshev_updates;
      const bool declared_store_difference =
          streamed_work.defect_writes == cells &&
          legacy_work.defect_writes == defects * cells &&
          streamed_work.defect_writes + stages * cells ==
              legacy_work.defect_writes &&
          streamed_work.chebyshev_elided_intermediate_residual_publications ==
              stages * cells &&
          legacy_work.chebyshev_elided_intermediate_residual_publications == 0U;
      const bool copyback_expected =
          streamed_work.chebyshev_copyback_cells ==
              ((degree & 1U) != 0U ? 2U * cells : 0U) &&
          legacy_work.chebyshev_copyback_cells == 0U;
      passed &= expect(
          static_cast<bool>(streamed_status) &&
              packed(streamed_status) == packed(legacy_status) &&
              bitwise_same(streamed.correction, legacy.correction) &&
              bitwise_same(streamed.plan.last_cycle_initial_residual(),
                           legacy.plan.last_cycle_initial_residual()) &&
              bitwise_same(streamed.plan.last_cycle_final_residual(),
                           legacy.plan.last_cycle_final_residual()),
          "streamed Chebyshev correction/residual/status is bitwise legacy-equivalent");
      passed &= expect(common_work && declared_store_difference &&
                           copyback_expected,
                       "streamed Chebyshev preserves work and declares eliminated stores");
    }
  }
  return passed;
}

bool test_streamed_chebyshev_nonfinite_failure_is_unpublished() {
  bool passed = true;
  for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
    Fixture streamed;
    Fixture legacy;
    const MgBoundarySet mixed{
        MgBoundaryKind::periodic, MgBoundaryKind::periodic,
        MgBoundaryKind::neumann, MgBoundaryKind::dirichlet,
        MgBoundaryKind::dirichlet, MgBoundaryKind::neumann};
    const bool fixtures_compiled =
        streamed.initialize(
            MgPointSmootherKind::chebyshev_jacobi,
            MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 8.0, 0.3,
            3U, mixed, MgNullSpace::none, false, {6, 4, 5}, true, 4.0) &&
        legacy.initialize(
            MgPointSmootherKind::chebyshev_jacobi,
            MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 8.0, 0.3,
            3U, mixed, MgNullSpace::none, false, {6, 4, 5}, true, 4.0);
    passed &= expect(fixtures_compiled,
                     "nonfinite streamed/legacy fixtures compile");
    if (!fixtures_compiled) continue;
    streamed.fill_rhs();
    legacy.fill_rhs();
    streamed.rhs.view.unchecked({1, 1, 1}, 0U) = invalid;
    legacy.rhs.view.unchecked({1, 1, 1}, 0U) = invalid;
    std::fill(streamed.correction.storage.begin(),
              streamed.correction.storage.end(), -13.0);
    std::fill(legacy.correction.storage.begin(), legacy.correction.storage.end(),
              -13.0);
    const std::vector<double> streamed_before = streamed.correction.storage;
    const std::vector<double> legacy_before = legacy.correction.storage;
    const MgPlanCounters streamed_counters_before = streamed.plan.counters();
    const MgPlanCounters legacy_counters_before = legacy.plan.counters();
    detail::set_mg_reference_chebyshev_lifecycle_for_test(streamed.plan, false);
    detail::set_mg_reference_chebyshev_lifecycle_for_test(legacy.plan, true);
    const Status streamed_status = streamed.plan.apply(
        as_const(streamed.rhs.view), streamed.correction.view, 0U);
    const Status legacy_status = legacy.plan.apply(
        as_const(legacy.rhs.view), legacy.correction.view, 0U);
    const auto unchanged = [](const OwnedField& field,
                              const std::vector<double>& before) noexcept {
      return field.storage.size() == before.size() &&
             std::memcmp(field.storage.data(), before.data(),
                         before.size() * sizeof(double)) == 0;
    };
    passed &= expect(
        streamed_status.code == StatusCode::numerical_failure &&
            packed(streamed_status) == packed(legacy_status) &&
            unchanged(streamed.correction, streamed_before) &&
            unchanged(legacy.correction, legacy_before) &&
            same_counters(streamed.plan.counters(), streamed_counters_before) &&
            same_counters(legacy.plan.counters(), legacy_counters_before),
        "streamed/legacy NaN and Inf failures match without public correction publication");
  }
  return passed;
}

bool test_non_cubic_seven_point_oracle() {
  const MgBoundarySet periodic{MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic};
  Fixture fixture;
  bool passed = expect(
      fixture.initialize(
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 8.0, 1.0,
          3U, periodic, MgNullSpace::none, false, {6, 4, 6}),
      "non-cubic seven-point certified oracle fixture compiles");
  if (!passed) return false;
  fixture.fill_rhs(true);
  detail::set_mg_skip_final_projection_for_test(fixture.plan, true);
  const Status status = fixture.plan.apply(as_const(fixture.rhs.view),
                                           fixture.correction.view, 0U);
  const detail::MgMatrixWorkCounters work =
      detail::mg_matrix_work_counters_for_test(fixture.plan);
  const double expected = scalar_two_pass_chebyshev(
      1.0, 8.0, 14.0, 3U,
      fixture.spec.policy.chebyshev_lower_spectrum_fraction);
  bool values_match = true;
  for (std::int32_t k = 0; k < fixture.patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < fixture.patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < fixture.patch.cells.x; ++i) {
        const double sign = ((i + j + k) & 1) == 0 ? 1.0 : -1.0;
        values_match = values_match &&
                       std::abs(fixture.correction.view.unchecked(
                                    {i, j, k}, 0U) -
                                sign * expected) < 3.0e-12;
      }
    }
  }
  const std::size_t cells = static_cast<std::size_t>(fixture.patch.cells.x) *
                            static_cast<std::size_t>(fixture.patch.cells.y) *
                            static_cast<std::size_t>(fixture.patch.cells.z);
  passed &= expect(
      static_cast<bool>(status) && values_match &&
          work.chebyshev_stages == 6U &&
          work.chebyshev_exchange_actions == 7U &&
          work.chebyshev_defect_actions == 7U &&
          work.chebyshev_retained_final_defect_actions == 1U &&
          work.chebyshev_updates == 6U * cells,
      "independent degree-three recurrence oracle covers a non-cubic seven-point periodic system");
  return passed;
}

bool test_public_identity_and_configuration_rejection() {
  Fixture red;
  Fixture red_certified;
  Fixture chebyshev;
  bool passed = expect(
      red.initialize(MgPointSmootherKind::red_black,
                     MgOperatorClass::general) &&
          red_certified.initialize(
              MgPointSmootherKind::red_black,
              MgOperatorClass::symmetric_diagonally_dominant_m_matrix) &&
          chebyshev.initialize(
              MgPointSmootherKind::chebyshev_jacobi,
              MgOperatorClass::symmetric_diagonally_dominant_m_matrix),
      "red and certified configurations compile");
  if (!passed) return false;
  passed &= expect(red.plan.symbolic_fingerprint() !=
                       chebyshev.plan.symbolic_fingerprint(),
                   "smoother kind mutates symbolic identity");
  passed &= expect(red.plan.symbolic_fingerprint() !=
                       red_certified.plan.symbolic_fingerprint(),
                   "operator class mutates symbolic identity");
  Fixture changed_fraction;
  passed &= expect(
      changed_fraction.initialize(
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix),
      "fraction identity fixture compiles");
  if (passed) {
    changed_fraction.spec.policy.chebyshev_lower_spectrum_fraction = 0.4;
    NativeCartesianMgPlan changed_fraction_plan;
    const Status changed_fraction_status = NativeCartesianMgPlan::compile(
        changed_fraction.spec, changed_fraction.services(),
        changed_fraction.coefficients(), changed_fraction_plan);
    passed &= expect(
        static_cast<bool>(changed_fraction_status) &&
            changed_fraction_plan.symbolic_fingerprint() !=
                changed_fraction.plan.symbolic_fingerprint(),
        "Chebyshev lower-spectrum fraction mutates symbolic identity");
  }
  Fixture invalid;
  passed &= expect(!invalid.initialize(MgPointSmootherKind::chebyshev_jacobi,
                                       MgOperatorClass::general),
                    "Chebyshev without the explicit operator certificate is rejected");
  return passed;
}

bool test_boundary_nullspace_and_activity_routes() {
  bool passed = true;
  const MgBoundarySet periodic{MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic,
                               MgBoundaryKind::periodic};
  for (const MgBoundarySet boundaries :
       {MgBoundarySet{},
        MgBoundarySet{MgBoundaryKind::neumann, MgBoundaryKind::neumann,
                      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet,
                      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet},
        periodic}) {
    Fixture fixture;
    passed &= expect(
        fixture.initialize(
            MgPointSmootherKind::chebyshev_jacobi,
            MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 7.0, 1.0,
            2U, boundaries),
        "certified smoother compiles for Dirichlet/Neumann/periodic faces");
    if (!passed) continue;
    fixture.fill_rhs();
    passed &= expect(static_cast<bool>(fixture.plan.apply(
                         as_const(fixture.rhs.view), fixture.correction.view,
                         0U)),
                     "certified smoother applies on mixed topology");
  }

  Fixture nullspace;
  passed &= expect(
      nullspace.initialize(
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 6.0, 1.0,
          2U, periodic, MgNullSpace::constant),
      "certified weakly dominant constant-null-space hierarchy compiles");
  if (passed) {
    nullspace.fill_rhs();
    passed &= expect(static_cast<bool>(nullspace.plan.apply(
                         as_const(nullspace.rhs.view),
                         nullspace.correction.view, 0U)),
                     "constant-null-space certified apply remains finite");
  }

  Fixture activity;
  passed &= expect(
      activity.initialize(
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix, 7.0, 1.0,
          2U, MgBoundarySet{}, MgNullSpace::none, true),
      "certified activity map compiles with an inactive identity row");
  if (passed) {
    activity.fill_rhs();
    passed &= expect(static_cast<bool>(activity.plan.apply(
                         as_const(activity.rhs.view), activity.correction.view,
                         0U)) &&
                         activity.correction.view.unchecked({0, 0, 0}, 0U) ==
                             0.0,
                     "inactive identity row remains zero in certified apply");
  }
  return passed;
}

bool test_numeric_input_rejection() {
  bool passed = true;
  for (const double invalid_diagonal : {0.0, -1.0,
                                        std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::infinity()}) {
    Fixture fixture;
    passed &= expect(
        fixture.initialize(MgPointSmootherKind::red_black,
                           MgOperatorClass::general),
        "numeric rejection baseline compiles");
    if (!passed) continue;
    fixture.diagonal.storage[0U] = invalid_diagonal;
    NativeCartesianMgPlan candidate;
    const Status status = NativeCartesianMgPlan::compile(
        fixture.spec, fixture.services(), fixture.coefficients(), candidate);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "zero/negative/NaN/Inf diagonal is rejected before publication");
  }

  Fixture negative_face;
  passed &= expect(negative_face.initialize(MgPointSmootherKind::red_black,
                                            MgOperatorClass::general),
                   "negative face rejection baseline compiles");
  if (passed) {
    negative_face.x_faces.storage[0U] = -1.0;
    NativeCartesianMgPlan candidate;
    const Status status = NativeCartesianMgPlan::compile(
        negative_face.spec, negative_face.services(),
        negative_face.coefficients(), candidate);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "negative face coefficient is rejected before publication");
  }

  Fixture invalid_fraction;
  passed &= expect(invalid_fraction.initialize(
                       MgPointSmootherKind::chebyshev_jacobi,
                       MgOperatorClass::
                           symmetric_diagonally_dominant_m_matrix),
                   "fraction validation baseline compiles");
  if (passed) {
    invalid_fraction.spec.policy.chebyshev_lower_spectrum_fraction = 1.0;
    NativeCartesianMgPlan candidate;
    const Status status = NativeCartesianMgPlan::compile(
        invalid_fraction.spec, invalid_fraction.services(),
        invalid_fraction.coefficients(), candidate);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "Chebyshev lower-spectrum fraction must be strictly inside (0,1)");
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_public_contract();
  passed &= test_certified_apply_and_hot_schedule();
  passed &= test_fixed_f_cycle_schedule_and_identity();
  passed &= test_independent_fixed_cycle_composition();
  passed &= test_certificate_rejection_is_transactional();
  passed &= test_degree_one_through_four_oracle();
  passed &= test_streamed_chebyshev_matches_legacy_degrees_and_shapes();
  passed &= test_streamed_chebyshev_nonfinite_failure_is_unpublished();
  passed &= test_non_cubic_seven_point_oracle();
  passed &= test_public_identity_and_configuration_rejection();
  passed &= test_boundary_nullspace_and_activity_routes();
  passed &= test_numeric_input_rejection();
  MPI_Finalize();
  return passed ? 0 : 1;
}
