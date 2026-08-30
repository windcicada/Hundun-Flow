// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "solver_hypre_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

using namespace hundun::v04;

static_assert(!std::is_copy_constructible_v<HypreStructAdapter>);
static_assert(!std::is_copy_assignable_v<HypreStructAdapter>);
static_assert(std::is_nothrow_move_constructible_v<HypreStructAdapter>);
static_assert(std::is_nothrow_move_assignable_v<HypreStructAdapter>);
static_assert(std::is_base_of_v<LinearPreconditioner, HypreStructAdapter>);

// If a HYPRE header leaks through the installed HUNDUN header, its public
// include guards make this HYPRE-off build fail here.
#if defined(HYPRE_H) || defined(HYPRE_STRUCT_LS_HEADER) || \
    defined(HYPRE_STRUCT_MV_HEADER) || defined(HYPRE_UTILITIES_HEADER)
#error "v04_linear.hpp must not expose HYPRE headers or types"
#endif

#if defined(HUNDUN_V04_TEST_HAVE_HYPRE)
#include <HYPRE.h>
#include <HYPRE_struct_mv.h>
#endif

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

#if defined(HUNDUN_V04_TEST_HAVE_HYPRE)
bool external_hypre_grid_round_trip() noexcept {
  HYPRE_StructGrid grid = nullptr;
  HYPRE_Int lower[3]{0, 0, 0};
  HYPRE_Int upper[3]{1, 1, 1};
  bool valid = HYPRE_StructGridCreate(MPI_COMM_SELF, 3, &grid) == 0;
  valid = valid && grid != nullptr &&
          HYPRE_StructGridSetExtents(grid, lower, upper) == 0 &&
          HYPRE_StructGridAssemble(grid) == 0;
  if (grid != nullptr) {
    valid = HYPRE_StructGridDestroy(grid) == 0 && valid;
  }
  if (!valid) {
    (void)HYPRE_ClearAllErrors();
  }
  return valid;
}
#endif

bool empty(LinearPreconditionerCertificate certificate) noexcept {
  return certificate.identity.symbolic == 0U &&
         certificate.identity.numeric == 0U &&
         certificate.identity.hierarchy == 0U &&
         certificate.identity.workspace == 0U &&
         certificate.identity.fingerprint == 0U &&
         certificate.collective_fingerprint == 0U;
}

bool same(LinearIdentity left, LinearIdentity right) noexcept {
  return left.symbolic == right.symbolic && left.numeric == right.numeric &&
         left.hierarchy == right.hierarchy &&
         left.workspace == right.workspace &&
         left.fingerprint == right.fingerprint;
}

bool all_true(bool local) noexcept {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

std::uint64_t packed(Status status) noexcept {
  return (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
}

bool identical(Status status) noexcept {
  const std::uint64_t local = packed(status);
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&local, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&local, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

bool test_unavailable_adapter_is_atomic() {
  bool passed = expect(!HypreStructAdapter::available(),
                       "HYPRE-off build reports the adapter unavailable");
  HypreStructAdapter adapter;
  const PlanFingerprint fingerprint = adapter.fingerprint();
  const PlanFingerprint numeric = adapter.numeric_fingerprint();
  const std::uintptr_t handles = adapter.native_handle_storage_address();
  const LinearPreconditionerCertificate certificate = adapter.certificate();

  NativeCartesianMgSpec spec;
  spec.communicator = MPI_COMM_SELF;
  spec.patch.cells = {4, 3, 2};
  spec.patch.process_grid = {1, 1, 1};
  spec.identity = {11U, 12U, 13U, 14U, 15U};
  MgRuntimeServices services{};
  const Status compiled =
      HypreStructAdapter::compile(spec, services, adapter);
  passed &= expect(compiled.code == StatusCode::invalid_plan,
                   "HYPRE-off compile returns an explicit unavailable status");
  passed &= expect(adapter.fingerprint() == fingerprint &&
                       adapter.numeric_fingerprint() == numeric &&
                       adapter.native_handle_storage_address() == handles &&
                       empty(adapter.certificate()) && empty(certificate),
                   "unavailable compile preserves the adapter atomically");

  const Status updated = adapter.update_coefficients(
      spec.identity, MgCoefficientIdentity{21U, 22U, 0.0},
      MgCoefficientViews{});
  passed &= expect(updated.code == StatusCode::invalid_plan &&
                       adapter.fingerprint() == fingerprint &&
                       adapter.numeric_fingerprint() == numeric &&
                       adapter.native_handle_storage_address() == handles,
                   "unavailable coefficient update publishes no native state");

  double residual_value = 1.0;
  double correction_value = -7.0;
  ConstFieldView residual;
  residual.base = &residual_value;
  residual.interior = {1, 1, 1};
  residual.components = 1U;
  residual.stride_y = 1U;
  residual.stride_z = 1U;
  residual.component_stride = 1U;
  residual.field = 1U;
  residual.revision = 1U;
  residual.storage_identity = 1U;
  residual.revision_domain = 1U;
  FieldView correction;
  correction.base = &correction_value;
  correction.interior = {1, 1, 1};
  correction.components = 1U;
  correction.stride_y = 1U;
  correction.stride_z = 1U;
  correction.component_stride = 1U;
  correction.field = 2U;
  correction.revision = 2U;
  correction.storage_identity = 2U;
  correction.revision_domain = 1U;
  const Status applied = adapter.apply(residual, correction, 0U);
  passed &= expect(applied.code == StatusCode::invalid_plan &&
                       correction_value == -7.0,
                   "unavailable apply leaves caller output unchanged");
  return passed;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField field(Int3 cells, FieldId id, StorageIdentity storage) {
  OwnedField result;
  const std::size_t count = static_cast<std::size_t>(cells.x) *
                            static_cast<std::size_t>(cells.y) *
                            static_cast<std::size_t>(cells.z);
  result.storage.assign(count, 0.0);
  result.view.base = result.storage.data();
  result.view.interior = cells;
  result.view.components = 1U;
  result.view.stride_y = static_cast<std::size_t>(cells.x);
  result.view.stride_z = result.view.stride_y *
                         static_cast<std::size_t>(cells.y);
  result.view.component_stride = count;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = storage;
  result.view.revision_domain = 90U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  ConstFaceFieldView view{};
};

OwnedFace face(Int3 extents, CartesianAxis axis, StorageIdentity storage) {
  OwnedFace result;
  const std::size_t count = static_cast<std::size_t>(extents.x) *
                            static_cast<std::size_t>(extents.y) *
                            static_cast<std::size_t>(extents.z);
  result.storage.assign(count, 1.0);
  result.view.base = result.storage.data();
  result.view.extents = extents;
  result.view.stride_y = static_cast<std::size_t>(extents.x);
  result.view.stride_z = result.view.stride_y *
                         static_cast<std::size_t>(extents.y);
  result.view.axis = axis;
  result.view.storage_identity = storage;
  result.view.revision_domain = 91U;
  return result;
}

CartesianMeshSpec mesh_spec(Int3 cells) {
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

double neighbor(ConstFieldView input, Int3 cell, CartesianAxis axis,
                int direction) noexcept {
  Int3 selected = cell;
  std::int32_t* coordinate = axis == CartesianAxis::x
                                 ? &selected.x
                                 : (axis == CartesianAxis::y ? &selected.y
                                                             : &selected.z);
  const std::int32_t extent = axis == CartesianAxis::x
                                  ? input.interior.x
                                  : (axis == CartesianAxis::y
                                         ? input.interior.y
                                         : input.interior.z);
  *coordinate += direction;
  return *coordinate >= 0 && *coordinate < extent
             ? input.unchecked(selected, 0U)
             : 0.0;
}

void apply_operator(ConstFieldView input, FieldView output,
                    ConstFieldView diagonal, ConstFaceFieldView x,
                    ConstFaceFieldView y, ConstFaceFieldView z) noexcept {
  const Int3 cells = input.interior;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        double value = diagonal.unchecked(cell, 0U) *
                       input.unchecked(cell, 0U);
        value -= x.unchecked({i, j, k}) *
                 neighbor(input, cell, CartesianAxis::x, -1);
        value -= x.unchecked({i + 1, j, k}) *
                 neighbor(input, cell, CartesianAxis::x, 1);
        value -= y.unchecked({i, j, k}) *
                 neighbor(input, cell, CartesianAxis::y, -1);
        value -= y.unchecked({i, j + 1, k}) *
                 neighbor(input, cell, CartesianAxis::y, 1);
        value -= z.unchecked({i, j, k}) *
                 neighbor(input, cell, CartesianAxis::z, -1);
        value -= z.unchecked({i, j, k + 1}) *
                 neighbor(input, cell, CartesianAxis::z, 1);
        output.unchecked(cell, 0U) = value;
      }
    }
  }
}

double residual_norm(ConstFieldView rhs, ConstFieldView solution,
                     OwnedField& applied, ConstFieldView diagonal,
                     ConstFaceFieldView x, ConstFaceFieldView y,
                     ConstFaceFieldView z) noexcept {
  apply_operator(solution, applied.view, diagonal, x, y, z);
  double sum = 0.0;
  const Int3 cells = rhs.interior;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double value = rhs.unchecked({i, j, k}, 0U) -
                             applied.view.unchecked({i, j, k}, 0U);
        sum += value * value;
      }
    }
  }
  return std::sqrt(sum);
}

bool test_available_adapter_reuses_native_handles() {
  detail::reset_hypre_lifecycle_for_test();
  bool passed = expect(HypreStructAdapter::available(),
                       "HYPRE-on build reports the adapter available");
  constexpr Int3 cells{12, 10, 8};
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  if (!expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                  MPI_COMM_SELF, mesh_spec(cells), GeometryBudget{}, geometry,
                  patch)),
              "HYPRE fixture geometry compiles")) {
    return false;
  }
  ReductionEngine reductions;
  if (!expect(static_cast<bool>(ReductionEngine::compile(
                  MPI_COMM_SELF, ReductionMode::mpi_allreduce, 4U,
                  reductions)),
              "HYPRE fixture reductions compile")) {
    return false;
  }
  OwnedField diagonal = field(cells, 1U, 100U);
  OwnedFace x = face({cells.x + 1, cells.y, cells.z}, CartesianAxis::x,
                     101U);
  OwnedFace y = face({cells.x, cells.y + 1, cells.z}, CartesianAxis::y,
                     102U);
  OwnedFace z = face({cells.x, cells.y, cells.z + 1}, CartesianAxis::z,
                     103U);
  std::fill(diagonal.storage.begin(), diagonal.storage.end(), 6.0);
  OwnedField exact = field(cells, 2U, 110U);
  OwnedField rhs = field(cells, 3U, 111U);
  OwnedField correction = field(cells, 4U, 112U);
  OwnedField solution = field(cells, 5U, 113U);
  OwnedField applied = field(cells, 6U, 114U);
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        exact.view.unchecked({i, j, k}, 0U) =
            std::sin(0.17 * static_cast<double>(i + 1)) *
            std::sin(0.21 * static_cast<double>(j + 1)) *
            std::sin(0.25 * static_cast<double>(k + 1));
      }
    }
  }
  apply_operator(as_const(exact.view), rhs.view, as_const(diagonal.view),
                 x.view, y.view, z.view);
  NativeCartesianMgSpec spec;
  spec.communicator = MPI_COMM_SELF;
  spec.geometry = &geometry;
  spec.patch = patch;
  spec.boundaries = {};
  spec.null_space = MgNullSpace::none;
  spec.policy.pre_sweeps = 2U;
  spec.policy.post_sweeps = 2U;
  spec.policy.maximum_levels = 16U;
  spec.identity = {11U, 12U, 13U, 14U, 15U};
  spec.coefficients = {1U, 21U, 0.0};
  MgRuntimeServices services{};
  services.reductions = &reductions;

  NativeCartesianMgSpec chebyshev_spec = spec;
  chebyshev_spec.operator_class =
      MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
  chebyshev_spec.policy.point_smoother =
      MgPointSmootherKind::chebyshev_jacobi;
  HypreStructAdapter chebyshev_rejected;
  const Status chebyshev_status = HypreStructAdapter::compile(
      chebyshev_spec, services, chebyshev_rejected);
  passed &= expect(
      chebyshev_status.code == StatusCode::invalid_plan &&
          chebyshev_rejected.fingerprint() == 0U,
      "isolated HYPRE adapter explicitly rejects certified Chebyshev policy");

  NativeCartesianMgSpec f_cycle_spec = spec;
  f_cycle_spec.policy.cycle = MgCycleKind::f_cycle;
  HypreStructAdapter f_cycle_rejected;
  const Status f_cycle_status = HypreStructAdapter::compile(
      f_cycle_spec, services, f_cycle_rejected);
  passed &= expect(
      f_cycle_status.code == StatusCode::invalid_plan &&
          f_cycle_rejected.fingerprint() == 0U,
      "isolated HYPRE adapter explicitly rejects native F-cycle policy");

#if defined(HUNDUN_V04_TEST_HAVE_HYPRE)
  detail::set_hypre_failure_for_test(
      detail::HypreFailurePoint::native_objects);
  HypreStructAdapter rejected_runtime;
  const Status runtime_failure =
      HypreStructAdapter::compile(spec, services, rejected_runtime);
  detail::clear_hypre_failure_for_test();
  passed &= expect(
      !runtime_failure && rejected_runtime.fingerprint() == 0U &&
          HYPRE_Initialized() != 0 && external_hypre_grid_round_trip(),
      "failed adapter construction preserves initialized HYPRE for external users");
#endif

  HypreStructAdapter adapter;
  passed &= expect(static_cast<bool>(
                       HypreStructAdapter::compile(spec, services, adapter)),
                   "HYPRE adapter compiles persistent topology objects");
  passed &= expect(empty(adapter.certificate()),
                   "adapter remains unpublished until numeric setup");
  const MgCoefficientViews coefficients{as_const(diagonal.view), x.view,
                                         y.view, z.view};
  passed &= expect(static_cast<bool>(adapter.update_coefficients(
                       spec.identity, spec.coefficients, coefficients)),
                   "HYPRE adapter assembles and sets up the first matrix");
  const std::uintptr_t handle_address =
      adapter.native_handle_storage_address();
  const PlanFingerprint numeric = adapter.numeric_fingerprint();
  passed &= expect(handle_address != 0U && numeric != 0U &&
                       !empty(adapter.certificate()) &&
                       adapter.certificate().preconditioner_class ==
                           LinearPreconditionerClass::fixed_general &&
                       adapter.certificate().status_scope ==
                           LinearPreconditionerStatusScope::rank_local &&
                       adapter.certificate().apply_lifecycle ==
                           LinearPreconditionerApplyLifecycle::per_call_checked,
                   "numeric setup publishes rank-local certificate and native handles");
  std::fill(correction.storage.begin(), correction.storage.end(), -91.0);
  passed &= expect(static_cast<bool>(adapter.apply(
                       as_const(rhs.view), correction.view, 0U)),
                   "one zero-guess PFMG cycle succeeds");
  bool finite = true;
  for (double value : correction.storage) {
    finite = finite && std::isfinite(value);
  }
  const double before = residual_norm(
      as_const(rhs.view), as_const(solution.view), applied,
      as_const(diagonal.view), x.view, y.view, z.view);
  solution.storage = correction.storage;
  solution.view.base = solution.storage.data();
  const double after = residual_norm(
      as_const(rhs.view), as_const(solution.view), applied,
      as_const(diagonal.view), x.view, y.view, z.view);
  passed &= expect(finite && std::isfinite(after) && after < before,
                   "one PFMG cycle returns finite residual-reducing output");

  passed &= expect(static_cast<bool>(adapter.update_coefficients(
                       spec.identity, spec.coefficients, coefficients)) &&
                       adapter.native_handle_storage_address() ==
                           handle_address &&
                       adapter.numeric_fingerprint() == numeric,
                   "unchanged identity performs no rebuild and reuses handles");
  LinearIdentity next_identity = spec.identity;
  MgCoefficientIdentity next_coefficient{};
  for (std::uint64_t update = 0U; update < 8U; ++update) {
    next_identity.numeric = 32U + update;
    next_identity.hierarchy = 42U + update;
    next_identity.fingerprint = 52U + update;
    next_coefficient = {2U + update, 22U + update, 0.1};
    std::fill(diagonal.storage.begin(), diagonal.storage.end(),
              6.5 + 0.125 * static_cast<double>(update));
    passed &= expect(
        static_cast<bool>(adapter.update_coefficients(
            next_identity, next_coefficient, coefficients)) &&
            adapter.native_handle_storage_address() == handle_address &&
            adapter.numeric_fingerprint() != numeric &&
            adapter.certificate().identity.numeric == next_identity.numeric,
        "eight changed identities refresh fresh numeric candidates while retaining topology handles");

    std::fill(correction.storage.begin(), correction.storage.end(), -81.0);
    passed &= expect(static_cast<bool>(adapter.apply(
                         as_const(rhs.view), correction.view,
                         static_cast<std::uint32_t>(update + 1U))),
                     "every changed numeric candidate remains applicable");
    bool update_finite = true;
    for (double value : correction.storage) {
      update_finite = update_finite && std::isfinite(value);
    }
    std::fill(solution.storage.begin(), solution.storage.end(), 0.0);
    const double update_before = residual_norm(
        as_const(rhs.view), as_const(solution.view), applied,
        as_const(diagonal.view), x.view, y.view, z.view);
    solution.storage = correction.storage;
    solution.view.base = solution.storage.data();
    const double update_after = residual_norm(
        as_const(rhs.view), as_const(solution.view), applied,
        as_const(diagonal.view), x.view, y.view, z.view);
    passed &= expect(update_finite && std::isfinite(update_after) &&
                         update_after < update_before,
                     "every changed numeric candidate reduces the residual");
  }
  const LinearPreconditionerCertificate accepted = adapter.certificate();
  const PlanFingerprint accepted_numeric = adapter.numeric_fingerprint();
  const std::uintptr_t accepted_handles =
      adapter.native_handle_storage_address();
  const MgCoefficientViews invalid_coefficients{
      ConstFieldView{}, x.view, y.view, z.view};
  LinearIdentity rejected_identity = next_identity;
  rejected_identity.numeric = 42U;
  rejected_identity.hierarchy = 43U;
  rejected_identity.fingerprint = 45U;
  const Status rejected = adapter.update_coefficients(
      rejected_identity, MgCoefficientIdentity{3U, 23U, 0.1},
      invalid_coefficients);
  passed &= expect(!rejected &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.certificate().identity.numeric ==
                           accepted.identity.numeric,
                   "rejected numeric update leaves the published native state unchanged");
  std::fill(correction.storage.begin(), correction.storage.end(), -73.0);
  passed &= expect(static_cast<bool>(adapter.apply(
                       as_const(rhs.view), correction.view, 1U)),
                   "adapter remains usable after a rejected numeric update");

  LinearIdentity setup_failure_identity = next_identity;
  setup_failure_identity.numeric = 62U;
  setup_failure_identity.hierarchy = 63U;
  setup_failure_identity.fingerprint = 65U;
  const MgCoefficientIdentity setup_failure_coefficient{12U, 32U, 0.1};
  std::fill(diagonal.storage.begin(), diagonal.storage.end(), 7.75);
  detail::set_hypre_failure_for_test(
      detail::HypreFailurePoint::solver_setup);
  const Status setup_failure = adapter.update_coefficients(
      setup_failure_identity, setup_failure_coefficient, coefficients);
  detail::clear_hypre_failure_for_test();
  passed &= expect(!setup_failure &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.certificate().identity.numeric ==
                           accepted.identity.numeric,
                   "failed fresh setup preserves the published active pair");
  std::fill(correction.storage.begin(), correction.storage.end(), -63.0);
  passed &= expect(static_cast<bool>(adapter.apply(
                       as_const(rhs.view), correction.view, 10U)),
                   "active pair remains usable after a failed fresh setup");
  passed &= expect(static_cast<bool>(adapter.update_coefficients(
                       setup_failure_identity, setup_failure_coefficient,
                       coefficients)) &&
                       adapter.certificate().identity.numeric ==
                           setup_failure_identity.numeric,
                   "failed fresh setup can be retried with a new pair");

  const detail::HypreLifecycleSnapshot before_fill_failure =
      detail::hypre_lifecycle_snapshot_for_test();
  LinearIdentity fill_failure_identity = setup_failure_identity;
  fill_failure_identity.numeric = 72U;
  fill_failure_identity.hierarchy = 73U;
  fill_failure_identity.fingerprint = 75U;
  detail::set_hypre_failure_for_test(detail::HypreFailurePoint::matrix_fill);
  const Status fill_failure = adapter.update_coefficients(
      fill_failure_identity, MgCoefficientIdentity{13U, 33U, 0.1},
      coefficients);
  detail::clear_hypre_failure_for_test();
  const detail::HypreLifecycleSnapshot after_fill_failure =
      detail::hypre_lifecycle_snapshot_for_test();
  passed &= expect(
      !fill_failure && adapter.certificate().identity.numeric ==
                           setup_failure_identity.numeric &&
          after_fill_failure.matrix_creates ==
              before_fill_failure.matrix_creates + 1U &&
          after_fill_failure.matrix_destroys ==
              before_fill_failure.matrix_destroys + 1U &&
          after_fill_failure.solver_creates ==
              before_fill_failure.solver_creates + 1U &&
          after_fill_failure.solver_destroys ==
              before_fill_failure.solver_destroys + 1U &&
          after_fill_failure.solver_setups ==
              before_fill_failure.solver_setups,
      "failed matrix fill destroys its fresh pair before setup and preserves active state");

  const detail::HypreLifecycleSnapshot lifecycle =
      detail::hypre_lifecycle_snapshot_for_test();
  const bool matrix_pair_live =
      lifecycle.matrix_creates >= lifecycle.matrix_destroys &&
      lifecycle.matrix_creates - lifecycle.matrix_destroys == 2U;
  const bool solver_pair_live =
      lifecycle.solver_creates >= lifecycle.solver_destroys &&
      lifecycle.solver_creates - lifecycle.solver_destroys == 2U;
  passed &= expect(lifecycle.matrix_creates == lifecycle.solver_creates &&
                       lifecycle.solver_creates ==
                           lifecycle.solver_setups + 1U &&
                       lifecycle.solver_setups >= 11U && matrix_pair_live &&
                       solver_pair_live,
                   "every setup attempt owns a fresh pair, the pre-setup failure is reclaimed, and only two pairs remain live");
  adapter = HypreStructAdapter{};
  const detail::HypreLifecycleSnapshot released =
      detail::hypre_lifecycle_snapshot_for_test();
  passed &= expect(released.matrix_creates == released.matrix_destroys &&
                       released.solver_creates == released.solver_destroys,
                   "adapter release balances every native pair lifetime");
#if defined(HUNDUN_V04_TEST_HAVE_HYPRE)
  const bool runtime_survives_release = HYPRE_Initialized() != 0;
  const bool external_round_trip =
      runtime_survives_release && external_hypre_grid_round_trip();
  passed &= expect(
      runtime_survives_release && external_round_trip,
      "adapter release preserves process HYPRE runtime for external users");
#endif
  return passed;
}

bool test_reversed_reduction_communicator_is_rejected(int rank, int size) {
  if (!HypreStructAdapter::available()) {
    return true;
  }
  const Int3 global_cells{4 * size, 8, 6};
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, mesh_spec(global_cells), GeometryBudget{}, geometry,
          patch)),
      "MPI HYPRE fixture geometry compiles");

  ReductionEngine reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_WORLD, ReductionMode::mpi_allreduce, 4U,
                       reductions)),
                   "MPI HYPRE fixture reductions compile");
  MPI_Comm reversed = MPI_COMM_NULL;
  passed &= expect(MPI_Comm_split(MPI_COMM_WORLD, 0, size - 1 - rank,
                                  &reversed) == MPI_SUCCESS,
                   "reversed communicator compiles");
  ReductionEngine reversed_reductions;
  if (reversed != MPI_COMM_NULL) {
    passed &= expect(static_cast<bool>(ReductionEngine::compile(
                         reversed, ReductionMode::mpi_allreduce, 4U,
                         reversed_reductions)),
                     "reversed reduction engine compiles");
  }
  if (!all_true(passed)) {
    if (reversed != MPI_COMM_NULL) {
      (void)MPI_Comm_free(&reversed);
    }
    return false;
  }

  OwnedField diagonal = field(patch.cells, 20U, 201U);
  OwnedFace x = face({patch.cells.x + 1, patch.cells.y, patch.cells.z},
                     CartesianAxis::x, 202U);
  OwnedFace y = face({patch.cells.x, patch.cells.y + 1, patch.cells.z},
                     CartesianAxis::y, 203U);
  OwnedFace z = face({patch.cells.x, patch.cells.y, patch.cells.z + 1},
                     CartesianAxis::z, 204U);
  std::fill(diagonal.storage.begin(), diagonal.storage.end(), 7.0);
  const MgCoefficientViews coefficients{as_const(diagonal.view), x.view,
                                         y.view, z.view};
  NativeCartesianMgSpec spec;
  spec.communicator = MPI_COMM_WORLD;
  spec.geometry = &geometry;
  spec.patch = patch;
  spec.policy.pre_sweeps = 1U;
  spec.policy.post_sweeps = 1U;
  spec.policy.maximum_levels = 16U;
  spec.identity = {101U, 102U, 103U, 104U, 105U};
  spec.coefficients = {201U, 202U, 0.0};

  MgRuntimeServices accepted_services{};
  accepted_services.reductions = &reductions;
  HypreStructAdapter adapter;
  passed &= expect(static_cast<bool>(HypreStructAdapter::compile(
                       spec, accepted_services, adapter)) &&
                       static_cast<bool>(adapter.update_coefficients(
                           spec.identity, spec.coefficients, coefficients)),
                   "accepted MPI HYPRE adapter publishes numeric state");
  const PlanFingerprint accepted_fingerprint = adapter.fingerprint();
  const PlanFingerprint accepted_numeric = adapter.numeric_fingerprint();
  const std::uintptr_t accepted_handles =
      adapter.native_handle_storage_address();
  const LinearPreconditionerCertificate accepted = adapter.certificate();

  MgRuntimeServices asymmetric_services{};
  asymmetric_services.reductions = rank == 0 ? nullptr : &reductions;
  const Status asymmetric_rejected =
      HypreStructAdapter::compile(spec, asymmetric_services, adapter);
  passed &= expect(!asymmetric_rejected &&
                       asymmetric_rejected.code == StatusCode::invalid_plan &&
                       identical(asymmetric_rejected) &&
                       adapter.fingerprint() == accepted_fingerprint &&
                       adapter.numeric_fingerprint() == accepted_numeric,
                   "spec communicator propagates a rank-local null reduction service without publishing state");

  MgRuntimeServices rejected_services{};
  rejected_services.reductions = &reversed_reductions;
  const Status rejected =
      HypreStructAdapter::compile(spec, rejected_services, adapter);
  passed &= expect(!rejected && rejected.code == StatusCode::invalid_plan &&
                       identical(rejected),
                   "MPI_SIMILAR reduction communicator is collectively rejected");
  passed &= expect(adapter.fingerprint() == accepted_fingerprint &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       same(adapter.certificate().identity,
                            accepted.identity),
                   "rejected compile preserves the accepted adapter");

  OwnedField residual = field(patch.cells, 21U, 211U);
  OwnedField correction = field(patch.cells, 22U, 212U);
  std::fill(residual.storage.begin(), residual.storage.end(), 1.0);
  std::fill(correction.storage.begin(), correction.storage.end(), -9.0);
  passed &= expect(static_cast<bool>(adapter.apply(
                       as_const(residual.view), correction.view, 0U)),
                   "accepted adapter remains reusable after rejection");
  for (double value : correction.storage) {
    passed &= std::isfinite(value);
  }

  const double reduction_local = static_cast<double>(rank + 1);
  double reduction_global = 0.0;
  const Status reduced = reversed_reductions.checked_sum(
      {&reduction_local, 1U}, {&reduction_global, 1U});
  const double expected =
      0.5 * static_cast<double>(size) * static_cast<double>(size + 1);
  passed &= expect(static_cast<bool>(reduced) && reduction_global == expected,
                   "rejected reduction engine remains reusable");

  LinearIdentity mutation_identity = spec.identity;
  mutation_identity.numeric = 302U;
  mutation_identity.hierarchy = 303U;
  mutation_identity.fingerprint = 305U;
  const MgCoefficientIdentity mutation_coefficient{301U, 302U, 0.1};

  ReductionEngine moved_reductions;
  if (rank == size - 1) {
    moved_reductions = std::move(reductions);
  }
  const Status moved_update = adapter.update_coefficients(
      mutation_identity, mutation_coefficient, coefficients);
  passed &= expect(!moved_update &&
                       moved_update.code == StatusCode::invalid_plan &&
                       identical(moved_update) &&
                       adapter.fingerprint() == accepted_fingerprint &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       same(adapter.certificate().identity,
                            accepted.identity),
                   "rank-local moved reduction service is rejected collectively before numeric publication");
  std::fill(correction.storage.begin(), correction.storage.end(), -17.0);
  const Status moved_apply =
      adapter.apply(as_const(residual.view), correction.view, 1U);
  const bool moved_output_unchanged = std::all_of(
      correction.storage.begin(), correction.storage.end(),
      [](double value) noexcept { return value == -17.0; });
  passed &= expect(!moved_apply &&
                       moved_apply.code == StatusCode::invalid_plan &&
                       identical(moved_apply) && moved_output_unchanged &&
                       adapter.fingerprint() == accepted_fingerprint &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       same(adapter.certificate().identity,
                            accepted.identity),
                   "rank-local moved reduction service is rejected collectively before HYPRE apply");
  if (rank == size - 1) {
    reductions = std::move(moved_reductions);
  }
  std::fill(correction.storage.begin(), correction.storage.end(), -15.0);
  passed &= expect(static_cast<bool>(adapter.apply(
                       as_const(residual.view), correction.view, 2U)),
                   "restoring the same reduction implementation restores adapter usability");

  ReductionEngine self_reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_SELF, ReductionMode::mpi_allreduce, 4U,
                       self_reductions)),
                   "replacement self reduction engine compiles on every rank");
  reductions = std::move(self_reductions);
  const Status rebound_update = adapter.update_coefficients(
      mutation_identity, mutation_coefficient, coefficients);
  passed &= expect(!rebound_update &&
                       rebound_update.code == StatusCode::invalid_plan &&
                       identical(rebound_update) &&
                       adapter.fingerprint() == accepted_fingerprint &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       same(adapter.certificate().identity,
                            accepted.identity),
                   "all-rank self-recompiled reduction service is rejected before numeric publication");
  std::fill(correction.storage.begin(), correction.storage.end(), -13.0);
  const Status rebound_apply =
      adapter.apply(as_const(residual.view), correction.view, 3U);
  const bool rebound_output_unchanged = std::all_of(
      correction.storage.begin(), correction.storage.end(),
      [](double value) noexcept { return value == -13.0; });
  passed &= expect(!rebound_apply &&
                       rebound_apply.code == StatusCode::invalid_plan &&
                       identical(rebound_apply) &&
                       rebound_output_unchanged &&
                       adapter.fingerprint() == accepted_fingerprint &&
                       adapter.numeric_fingerprint() == accepted_numeric &&
                       adapter.native_handle_storage_address() ==
                           accepted_handles &&
                       same(adapter.certificate().identity,
                            accepted.identity),
                   "all-rank self-recompiled reduction service is rejected before HYPRE apply");

  if (reversed != MPI_COMM_NULL) {
    passed &= expect(MPI_Comm_free(&reversed) == MPI_SUCCESS,
                     "reversed communicator releases");
  }
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 1;
  const bool mpi_ready = MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS &&
                         MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS;
  const bool passed =
      mpi_ready && size > 1
          ? test_reversed_reduction_communicator_is_rejected(rank, size)
          : (HypreStructAdapter::available()
                 ? test_available_adapter_reuses_native_handles()
                 : test_unavailable_adapter_is_atomic());
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
