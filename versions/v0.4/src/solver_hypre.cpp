// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "solver_hypre_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <utility>
#include <vector>

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
#include <atomic>
#endif

#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
#include <HYPRE.h>
#include <HYPRE_struct_ls.h>
#include <HYPRE_struct_mv.h>

#include <mutex>
#endif

namespace hundun::v04 {
namespace {

#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
constexpr std::uint32_t kHyprePlan = 7202U;
constexpr std::uint32_t kHypreCollective = 7203U;
constexpr std::uint32_t kHypreCoefficient = 7204U;
constexpr std::uint32_t kHypreApply = 7205U;
constexpr std::uint32_t kHypreRuntime = 7206U;
#else
constexpr std::uint32_t kHypreUnavailable = 7201U;
#endif

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
struct HypreLifecycleCounters {
  std::atomic<std::uint64_t> matrix_creates{};
  std::atomic<std::uint64_t> matrix_destroys{};
  std::atomic<std::uint64_t> solver_creates{};
  std::atomic<std::uint64_t> solver_destroys{};
  std::atomic<std::uint64_t> solver_setups{};
};

HypreLifecycleCounters& hypre_lifecycle_counters() noexcept {
  static HypreLifecycleCounters counters;
  return counters;
}

std::atomic<detail::HypreFailurePoint> g_hypre_failure_point{
    detail::HypreFailurePoint::none};
#endif

#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE && \
    defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
HYPRE_Int create_matrix(MPI_Comm communicator, HYPRE_StructGrid grid,
                        HYPRE_StructStencil stencil,
                        HYPRE_StructMatrix* matrix) noexcept {
  const HYPRE_Int result =
      HYPRE_StructMatrixCreate(communicator, grid, stencil, matrix);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (result == 0) {
    hypre_lifecycle_counters().matrix_creates.fetch_add(
        1U, std::memory_order_relaxed);
  }
#endif
  return result;
}

HYPRE_Int destroy_matrix(HYPRE_StructMatrix matrix) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  hypre_lifecycle_counters().matrix_destroys.fetch_add(
      1U, std::memory_order_relaxed);
#endif
  return HYPRE_StructMatrixDestroy(matrix);
}

HYPRE_Int fill_matrix(HYPRE_StructMatrix matrix, HYPRE_Int* lower,
                      HYPRE_Int* upper, HYPRE_Int entry_count,
                      HYPRE_Int* entries, HYPRE_Complex* values) noexcept {
  if (g_hypre_failure_point.load(std::memory_order_relaxed) ==
      detail::HypreFailurePoint::matrix_fill) {
    return -1;
  }
  return HYPRE_StructMatrixSetBoxValues(matrix, lower, upper, entry_count,
                                        entries, values);
}

HYPRE_Int create_solver(MPI_Comm communicator,
                        HYPRE_StructSolver* solver) noexcept {
  const HYPRE_Int result = HYPRE_StructPFMGCreate(communicator, solver);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (result == 0) {
    hypre_lifecycle_counters().solver_creates.fetch_add(
        1U, std::memory_order_relaxed);
  }
#endif
  return result;
}

HYPRE_Int destroy_solver(HYPRE_StructSolver solver) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  hypre_lifecycle_counters().solver_destroys.fetch_add(
      1U, std::memory_order_relaxed);
#endif
  return HYPRE_StructPFMGDestroy(solver);
}

HYPRE_Int setup_solver(HYPRE_StructSolver solver,
                       HYPRE_StructMatrix matrix, HYPRE_StructVector rhs,
                       HYPRE_StructVector solution) noexcept {
  hypre_lifecycle_counters().solver_setups.fetch_add(
      1U, std::memory_order_relaxed);
  if (g_hypre_failure_point.load(std::memory_order_relaxed) ==
      detail::HypreFailurePoint::solver_setup) {
    return -1;
  }
  return HYPRE_StructPFMGSetup(solver, matrix, rhs, solution);
}
#elif defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
#define create_matrix HYPRE_StructMatrixCreate
#define destroy_matrix HYPRE_StructMatrixDestroy
#define create_solver HYPRE_StructPFMGCreate
#define destroy_solver HYPRE_StructPFMGDestroy
#define setup_solver HYPRE_StructPFMGSetup
#define fill_matrix HYPRE_StructMatrixSetBoxValues
#endif

#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
  return hash;
}

PlanFingerprint finish(std::uint64_t hash) noexcept {
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same_identity(LinearIdentity left, LinearIdentity right) noexcept {
  return left.symbolic == right.symbolic && left.numeric == right.numeric &&
         left.hierarchy == right.hierarchy &&
         left.workspace == right.workspace &&
         left.fingerprint == right.fingerprint;
}

bool valid_identity(LinearIdentity identity) noexcept {
  return identity.symbolic != 0U && identity.numeric != 0U &&
         identity.hierarchy != 0U && identity.workspace != 0U &&
         identity.fingerprint != 0U;
}

bool valid_pair(MgBoundaryKind minimum, MgBoundaryKind maximum) noexcept {
  return (minimum == MgBoundaryKind::periodic) ==
         (maximum == MgBoundaryKind::periodic);
}

bool valid_policy(MgHierarchyPolicy policy) noexcept {
  return std::isfinite(policy.anisotropy_threshold) &&
         policy.anisotropy_threshold > 1.0 &&
         std::isfinite(policy.coefficient_change_rebuild_ratio) &&
         policy.coefficient_change_rebuild_ratio >= 0.0 &&
         policy.pre_sweeps != 0U && policy.post_sweeps != 0U &&
         policy.maximum_levels >= 2U &&
         policy.maximum_levels <= kMgMaximumLevels &&
         policy.minimum_coarse_extent >= 2U;
}

bool positive_shape(Int3 shape) noexcept {
  return shape.x > 0 && shape.y > 0 && shape.z > 0;
}

std::size_t cell_count(Int3 shape) noexcept {
  if (!positive_shape(shape)) {
    return 0U;
  }
  const std::size_t x = static_cast<std::size_t>(shape.x);
  const std::size_t y = static_cast<std::size_t>(shape.y);
  const std::size_t z = static_cast<std::size_t>(shape.z);
  if (x > std::numeric_limits<std::size_t>::max() / y ||
      x * y > std::numeric_limits<std::size_t>::max() / z) {
    return 0U;
  }
  return x * y * z;
}

bool valid_scalar(ConstFieldView field, Int3 shape) noexcept {
  return field.base != nullptr && field.components == 1U &&
         same_shape(field.interior, shape) && field.ghosts.x >= 0 &&
         field.ghosts.y >= 0 && field.ghosts.z >= 0 &&
         field.stride_y >= static_cast<std::size_t>(shape.x) &&
         field.stride_z >=
             field.stride_y * static_cast<std::size_t>(shape.y) &&
         field.component_stride >=
             field.stride_z * static_cast<std::size_t>(shape.z) &&
         field.storage_identity != 0U && field.revision_domain != 0U;
}

bool valid_face(ConstFaceFieldView field, CartesianAxis axis,
                Int3 cells) noexcept {
  Int3 expected = cells;
  if (axis == CartesianAxis::x) {
    ++expected.x;
  } else if (axis == CartesianAxis::y) {
    ++expected.y;
  } else {
    ++expected.z;
  }
  return field.base != nullptr && field.axis == axis &&
         same_shape(field.extents, expected) &&
         field.stride_y >= static_cast<std::size_t>(expected.x) &&
         field.stride_z >=
             field.stride_y * static_cast<std::size_t>(expected.y) &&
         field.storage_identity != 0U && field.revision_domain != 0U;
}

bool valid_coefficient_views(MgCoefficientViews coefficients,
                             Int3 cells) noexcept {
  if (!valid_scalar(coefficients.diagonal, cells) ||
      !valid_face(coefficients.x, CartesianAxis::x, cells) ||
      !valid_face(coefficients.y, CartesianAxis::y, cells) ||
      !valid_face(coefficients.z, CartesianAxis::z, cells)) {
    return false;
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double diagonal =
            coefficients.diagonal.unchecked({i, j, k}, 0U);
        if (!(diagonal > 0.0) || !std::isfinite(diagonal)) {
          return false;
        }
      }
    }
  }
  const auto finite_faces = [](ConstFaceFieldView field) noexcept {
    for (std::int32_t k = 0; k < field.extents.z; ++k) {
      for (std::int32_t j = 0; j < field.extents.y; ++j) {
        for (std::int32_t i = 0; i < field.extents.x; ++i) {
          const double value = field.unchecked({i, j, k});
          if (value < 0.0 || !std::isfinite(value)) {
            return false;
          }
        }
      }
    }
    return true;
  };
  return finite_faces(coefficients.x) && finite_faces(coefficients.y) &&
         finite_faces(coefficients.z);
}

bool mpi_live() noexcept {
  int initialized = 0;
  int finalized = 0;
  (void)MPI_Initialized(&initialized);
  if (initialized != 0) {
    (void)MPI_Finalized(&finalized);
  }
  return initialized != 0 && finalized == 0;
}

struct HypreRuntimeState {
  std::mutex mutex;
};

HypreRuntimeState& runtime_state() noexcept {
  static HypreRuntimeState state;
  return state;
}

Status ensure_runtime() noexcept {
  HypreRuntimeState& state = runtime_state();
  std::lock_guard<std::mutex> lock(state.mutex);
  if (HYPRE_Initialized() == 0) {
    if (HYPRE_Initialize() != 0) {
      (void)HYPRE_ClearAllErrors();
      return {StatusCode::invalid_plan, kHypreRuntime};
    }
  }
  // HYPRE is process-wide and may be shared with callers outside HUNDUN.
  // Adapter destruction owns only adapter-native objects; it must never
  // terminate the shared runtime.  The embedding application retains final
  // process shutdown authority after HUNDUN's initialize-on-first-use path.
  return {};
}

Status hypre_status(HYPRE_Int code, std::uint32_t detail) noexcept {
  if (code == 0) {
    return {};
  }
  (void)HYPRE_ClearAllErrors();
  return {StatusCode::numerical_failure, detail};
}

bool valid_status_code(StatusCode code) noexcept {
  return static_cast<std::uint16_t>(code) <=
         static_cast<std::uint16_t>(StatusCode::io_failure);
}

Status raw_consensus(MPI_Comm communicator, Status local) noexcept {
  if (communicator == MPI_COMM_NULL) {
    return local ? Status{StatusCode::invalid_plan, kHyprePlan} : local;
  }
  if (!valid_status_code(local.code)) {
    local = {StatusCode::invalid_plan, kHyprePlan};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kHypreCollective};
  }
  const int candidate = local ? size : rank;
  int selected = size;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kHypreCollective};
  }
  if (selected == size) {
    return {};
  }
  std::uint64_t wire[2]{};
  if (rank == selected) {
    wire[0] = static_cast<std::uint64_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire, 2, MPI_UINT64_T, selected, communicator) !=
          MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    return {StatusCode::mpi_failure, kHypreCollective};
  }
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

PlanFingerprint structural_fingerprint(
    const NativeCartesianMgSpec& spec) noexcept {
  const Int3 global = spec.geometry->global_cells();
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, spec.geometry->fingerprint());
  hash = mix(hash, static_cast<std::uint64_t>(global.x));
  hash = mix(hash, static_cast<std::uint64_t>(global.y));
  hash = mix(hash, static_cast<std::uint64_t>(global.z));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.null_space));
  hash = mix(hash, spec.identity.symbolic);
  hash = mix(hash, spec.identity.workspace);
  hash = mix(hash, spec.policy.pre_sweeps);
  hash = mix(hash, spec.policy.post_sweeps);
  hash = mix(hash, spec.policy.maximum_levels);
  hash = mix(hash, spec.policy.minimum_coarse_extent);
  return finish(hash);
}

PlanFingerprint make_numeric_fingerprint(
    LinearIdentity identity, MgCoefficientIdentity coefficient,
    PlanFingerprint structural) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, structural);
  hash = mix(hash, identity.numeric);
  hash = mix(hash, identity.hierarchy);
  hash = mix(hash, coefficient.revision);
  hash = mix(hash, coefficient.fingerprint);
  return finish(hash);
}

PlanFingerprint update_contract(LinearIdentity identity,
                                MgCoefficientIdentity coefficient) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, identity.symbolic);
  hash = mix(hash, identity.numeric);
  hash = mix(hash, identity.hierarchy);
  hash = mix(hash, identity.workspace);
  hash = mix(hash, identity.fingerprint);
  hash = mix(hash, coefficient.revision);
  hash = mix(hash, coefficient.fingerprint);
  return finish(hash);
}

Status validate_compile(const NativeCartesianMgSpec& spec,
                        MgRuntimeServices services) noexcept {
  if (!mpi_live() || spec.communicator == MPI_COMM_NULL ||
      spec.geometry == nullptr || services.reductions == nullptr ||
      services.reductions->capacity() < 2U ||
      !positive_shape(spec.patch.cells) ||
      !valid_identity(spec.identity) || spec.coefficients.revision == 0U ||
      spec.coefficients.fingerprint == 0U ||
      !std::isfinite(spec.coefficients.maximum_relative_change) ||
      spec.coefficients.maximum_relative_change < 0.0 ||
      !valid_policy(spec.policy) ||
      !valid_pair(spec.boundaries.x_min, spec.boundaries.x_max) ||
      !valid_pair(spec.boundaries.y_min, spec.boundaries.y_max) ||
      !valid_pair(spec.boundaries.z_min, spec.boundaries.z_max)) {
    return {StatusCode::invalid_plan, kHyprePlan};
  }
  const Int3 global = spec.geometry->global_cells();
  const auto within = [](std::int32_t begin, std::int32_t count,
                         std::int32_t extent) noexcept {
    return begin >= 0 && count > 0 && begin <= extent &&
           count <= extent - begin;
  };
  return within(spec.patch.begin.x, spec.patch.cells.x, global.x) &&
                 within(spec.patch.begin.y, spec.patch.cells.y, global.y) &&
                 within(spec.patch.begin.z, spec.patch.cells.z, global.z)
             ? Status{}
             : Status{StatusCode::invalid_plan, kHyprePlan};
}
#endif

}  // namespace

struct HypreStructAdapter::Impl {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  MPI_Comm communicator{MPI_COMM_NULL};
  NativeCartesianMgSpec spec{};
  MgRuntimeServices services{};
  ReductionEngine* reductions_object{};
  std::uintptr_t reductions_identity{};
  HYPRE_StructGrid grid{};
  HYPRE_StructStencil stencil{};
  HYPRE_StructMatrix matrices[2]{};
  HYPRE_StructVector rhs{};
  HYPRE_StructVector solution{};
  HYPRE_StructSolver solvers[2]{};
  HYPRE_Int lower[3]{};
  HYPRE_Int upper[3]{};
  HYPRE_Int entries[7]{0, 1, 2, 3, 4, 5, 6};
  std::vector<HYPRE_Complex> matrix_values;
  std::vector<HYPRE_Complex> rhs_values;
  std::vector<HYPRE_Complex> solution_values;
  MgCoefficientIdentity coefficient_identity{};
  PlanFingerprint structural{};
  PlanFingerprint numeric{};
  std::uint8_t active_slot{};
  bool pair_setup[2]{};
  bool ready{};
#endif
};

#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
namespace {

template <class Implementation>
Status validate_borrowed_services(
    const Implementation& implementation) noexcept {
  if (implementation.services.reductions == nullptr ||
      implementation.services.reductions != implementation.reductions_object ||
      implementation.services.reductions->instance_identity() !=
          implementation.reductions_identity) {
    return {StatusCode::invalid_plan, kHyprePlan};
  }
  return implementation.services.reductions->validate_communicator(
      implementation.communicator);
}

template <class Implementation>
Status preflight_borrowed_services(
    const Implementation& implementation) noexcept {
  return raw_consensus(implementation.communicator,
                       validate_borrowed_services(implementation));
}

template <class Implementation>
void destroy_impl(Implementation* implementation) noexcept {
  if (implementation == nullptr) {
    return;
  }
  if (mpi_live()) {
    for (HYPRE_StructSolver& solver : implementation->solvers) {
      if (solver != nullptr) {
        (void)destroy_solver(solver);
      }
    }
    if (implementation->solution != nullptr) {
      (void)HYPRE_StructVectorDestroy(implementation->solution);
    }
    if (implementation->rhs != nullptr) {
      (void)HYPRE_StructVectorDestroy(implementation->rhs);
    }
    for (HYPRE_StructMatrix& matrix : implementation->matrices) {
      if (matrix != nullptr) {
        (void)destroy_matrix(matrix);
      }
    }
    if (implementation->stencil != nullptr) {
      (void)HYPRE_StructStencilDestroy(implementation->stencil);
    }
    if (implementation->grid != nullptr) {
      (void)HYPRE_StructGridDestroy(implementation->grid);
    }
    if (implementation->communicator != MPI_COMM_NULL) {
      (void)MPI_Comm_free(&implementation->communicator);
    }
  }
  delete implementation;
}

template <class Implementation>
Status consensus(Implementation& implementation, Status local) noexcept {
  return implementation.services.reductions->consensus(local);
}

void destroy_native_pair(HYPRE_StructMatrix& matrix,
                         HYPRE_StructSolver& solver) noexcept {
  if (solver != nullptr) {
    (void)destroy_solver(solver);
    solver = nullptr;
  }
  if (matrix != nullptr) {
    (void)destroy_matrix(matrix);
    matrix = nullptr;
  }
}

template <class Implementation>
Status create_native_pair(Implementation& implementation,
                          HYPRE_StructMatrix& matrix,
                          HYPRE_StructSolver& solver) noexcept {
  matrix = nullptr;
  solver = nullptr;
  Status local = hypre_status(
      create_matrix(implementation.communicator, implementation.grid,
                    implementation.stencil, &matrix),
      kHypreRuntime);
  if (local) {
    local = hypre_status(HYPRE_StructMatrixSetSymmetric(matrix, 0),
                         kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructMatrixInitialize(matrix), kHypreRuntime);
  }
  if (local) {
    local = hypre_status(create_solver(implementation.communicator, &solver),
                         kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructPFMGSetTol(solver, 0.0), kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructPFMGSetMaxIter(solver, 1), kHypreRuntime);
  }
  if (local) {
    local = hypre_status(
        HYPRE_StructPFMGSetMaxLevels(
            solver,
            static_cast<HYPRE_Int>(implementation.spec.policy.maximum_levels)),
        kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructPFMGSetZeroGuess(solver), kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructPFMGSetRelaxType(solver, 2),
                         kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructPFMGSetRAPType(solver, 0), kHypreRuntime);
  }
  if (local) {
    local = hypre_status(
        HYPRE_StructPFMGSetNumPreRelax(
            solver,
            static_cast<HYPRE_Int>(implementation.spec.policy.pre_sweeps)),
        kHypreRuntime);
  }
  if (local) {
    local = hypre_status(
        HYPRE_StructPFMGSetNumPostRelax(
            solver,
            static_cast<HYPRE_Int>(implementation.spec.policy.post_sweeps)),
        kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructPFMGSetLogging(solver, 0), kHypreRuntime);
  }
  return local;
}

template <class Implementation>
Status make_native_objects(Implementation& implementation) noexcept {
  Status local = ensure_runtime();
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (local &&
      g_hypre_failure_point.load(std::memory_order_relaxed) ==
          detail::HypreFailurePoint::native_objects) {
    local = {StatusCode::invalid_plan, kHypreRuntime};
  }
#endif
  Status agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  (void)HYPRE_ClearAllErrors();
  local = hypre_status(
      HYPRE_StructGridCreate(implementation.communicator, 3,
                             &implementation.grid),
      kHypreRuntime);
  if (local) {
    local = hypre_status(HYPRE_StructGridSetExtents(
                             implementation.grid, implementation.lower,
                             implementation.upper),
                         kHypreRuntime);
  }
  const Int3 global = implementation.spec.geometry->global_cells();
  HYPRE_Int periodic[3]{
      implementation.spec.boundaries.x_min == MgBoundaryKind::periodic
          ? static_cast<HYPRE_Int>(global.x)
          : 0,
      implementation.spec.boundaries.y_min == MgBoundaryKind::periodic
          ? static_cast<HYPRE_Int>(global.y)
          : 0,
      implementation.spec.boundaries.z_min == MgBoundaryKind::periodic
          ? static_cast<HYPRE_Int>(global.z)
          : 0};
  if (local && (periodic[0] != 0 || periodic[1] != 0 || periodic[2] != 0)) {
    local = hypre_status(
        HYPRE_StructGridSetPeriodic(implementation.grid, periodic),
        kHypreRuntime);
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructGridAssemble(implementation.grid),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }

  local = hypre_status(HYPRE_StructStencilCreate(3, 7,
                                                  &implementation.stencil),
                       kHypreRuntime);
  HYPRE_Int offsets[7][3]{{0, 0, 0},  {-1, 0, 0}, {1, 0, 0},
                          {0, -1, 0}, {0, 1, 0},  {0, 0, -1},
                          {0, 0, 1}};
  for (HYPRE_Int entry = 0; entry < 7 && local; ++entry) {
    local = hypre_status(HYPRE_StructStencilSetEntry(
                             implementation.stencil, entry, offsets[entry]),
                         kHypreRuntime);
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }

  for (std::size_t slot = 0U; slot < 2U && local; ++slot) {
    local = create_native_pair(implementation, implementation.matrices[slot],
                               implementation.solvers[slot]);
  }
  if (local) {
    local = hypre_status(HYPRE_StructVectorCreate(
                             implementation.communicator, implementation.grid,
                             &implementation.rhs),
                         kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructVectorCreate(
                             implementation.communicator, implementation.grid,
                             &implementation.solution),
                         kHypreRuntime);
  }
  if (local) {
    local = hypre_status(HYPRE_StructVectorInitialize(implementation.rhs),
                         kHypreRuntime);
  }
  if (local) {
    local = hypre_status(
        HYPRE_StructVectorInitialize(implementation.solution), kHypreRuntime);
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }

  std::fill(implementation.rhs_values.begin(),
            implementation.rhs_values.end(), HYPRE_Complex{0.0});
  local = hypre_status(HYPRE_StructVectorSetBoxValues(
                           implementation.rhs, implementation.lower,
                           implementation.upper,
                           implementation.rhs_values.data()),
                       kHypreRuntime);
  if (local) {
    local = hypre_status(HYPRE_StructVectorSetBoxValues(
                             implementation.solution, implementation.lower,
                             implementation.upper,
                             implementation.rhs_values.data()),
                         kHypreRuntime);
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructVectorAssemble(implementation.rhs),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructVectorAssemble(implementation.solution),
                       kHypreCollective);
  return consensus(implementation, local);
}

MgBoundaryKind minimum_boundary(const MgBoundarySet& boundaries,
                                CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? boundaries.x_min
             : (axis == CartesianAxis::y ? boundaries.y_min
                                         : boundaries.z_min);
}

MgBoundaryKind maximum_boundary(const MgBoundarySet& boundaries,
                                CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? boundaries.x_max
             : (axis == CartesianAxis::y ? boundaries.y_max
                                         : boundaries.z_max);
}

template <class Implementation>
Status fill_matrix_values(Implementation& implementation,
                          MgCoefficientViews coefficients,
                          double& local_anchor) noexcept {
  const Int3 cells = implementation.spec.patch.cells;
  const Int3 begin = implementation.spec.patch.begin;
  const Int3 global = implementation.spec.geometry->global_cells();
  const MgBoundarySet boundaries = implementation.spec.boundaries;
  local_anchor = 0.0;
  std::size_t cursor = 0U;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        double diagonal = coefficients.diagonal.unchecked({i, j, k}, 0U);
        double off_diagonal_sum = 0.0;
        double values[6]{coefficients.x.unchecked({i, j, k}),
                         coefficients.x.unchecked({i + 1, j, k}),
                         coefficients.y.unchecked({i, j, k}),
                         coefficients.y.unchecked({i, j + 1, k}),
                         coefficients.z.unchecked({i, j, k}),
                         coefficients.z.unchecked({i, j, k + 1})};
        const std::int32_t coordinate[3]{begin.x + i, begin.y + j,
                                         begin.z + k};
        const std::int32_t extent[3]{global.x, global.y, global.z};
        for (std::size_t axis = 0U; axis < 3U; ++axis) {
          const CartesianAxis selected =
              axis == 0U ? CartesianAxis::x
                         : (axis == 1U ? CartesianAxis::y : CartesianAxis::z);
          const bool at_minimum = coordinate[axis] == 0;
          const bool at_maximum = coordinate[axis] == extent[axis] - 1;
          if (at_minimum &&
              minimum_boundary(boundaries, selected) !=
                  MgBoundaryKind::periodic) {
            if (minimum_boundary(boundaries, selected) ==
                MgBoundaryKind::neumann) {
              diagonal -= values[2U * axis];
            }
            values[2U * axis] = 0.0;
          }
          if (at_maximum &&
              maximum_boundary(boundaries, selected) !=
                  MgBoundaryKind::periodic) {
            if (maximum_boundary(boundaries, selected) ==
                MgBoundaryKind::neumann) {
              diagonal -= values[2U * axis + 1U];
            }
            values[2U * axis + 1U] = 0.0;
          }
        }
        for (double value : values) {
          off_diagonal_sum += value;
        }
        const double tolerance =
            64.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(diagonal), off_diagonal_sum});
        if (!(diagonal > 0.0) || !std::isfinite(diagonal) ||
            diagonal + tolerance < off_diagonal_sum) {
          return {StatusCode::numerical_failure, kHypreCoefficient};
        }
        local_anchor =
            std::max(local_anchor, std::max(0.0, diagonal - off_diagonal_sum));
        implementation.matrix_values[cursor++] =
            static_cast<HYPRE_Complex>(diagonal);
        for (double value : values) {
          implementation.matrix_values[cursor++] =
              static_cast<HYPRE_Complex>(-value);
        }
      }
    }
  }
  return cursor == implementation.matrix_values.size()
             ? Status{}
             : Status{StatusCode::invalid_plan, kHypreCoefficient};
}

template <class Implementation>
Status project_values(Implementation& implementation,
                      std::vector<HYPRE_Complex>& values) noexcept {
  const MeshPatch patch = implementation.spec.patch;
  const AxisMetrics& x = implementation.spec.geometry->x();
  const AxisMetrics& y = implementation.spec.geometry->y();
  const AxisMetrics& z = implementation.spec.geometry->z();
  const Span<const double> dx = x.widths();
  const Span<const double> dy = y.widths();
  const Span<const double> dz = z.widths();
  double local[2]{};
  std::size_t cursor = 0U;
  for (std::int32_t k = 0; k < patch.cells.z; ++k) {
    const std::size_t gk = static_cast<std::size_t>(patch.begin.z + k);
    for (std::int32_t j = 0; j < patch.cells.y; ++j) {
      const std::size_t gj = static_cast<std::size_t>(patch.begin.y + j);
      for (std::int32_t i = 0; i < patch.cells.x; ++i) {
        const std::size_t gi = static_cast<std::size_t>(patch.begin.x + i);
        const double volume = dx.data[gi] * dy.data[gj] * dz.data[gk];
        local[0] += static_cast<double>(values[cursor++]) * volume;
        local[1] += volume;
      }
    }
  }
  double global[2]{};
  Status status = implementation.services.reductions->checked_sum(
      {local, 2U}, {global, 2U});
  if (!status || !(global[1] > 0.0) || !std::isfinite(global[0]) ||
      !std::isfinite(global[1])) {
    return status ? Status{StatusCode::numerical_failure, kHypreApply}
                  : status;
  }
  const double mean = global[0] / global[1];
  for (HYPRE_Complex& value : values) {
    value = static_cast<HYPRE_Complex>(static_cast<double>(value) - mean);
  }
  return {};
}

}  // namespace
#endif

HypreStructAdapter::~HypreStructAdapter() noexcept { release(); }

HypreStructAdapter::HypreStructAdapter(HypreStructAdapter&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

HypreStructAdapter& HypreStructAdapter::operator=(
    HypreStructAdapter&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void HypreStructAdapter::release() noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  destroy_impl(std::exchange(implementation_, nullptr));
#else
  delete std::exchange(implementation_, nullptr);
#endif
}

bool HypreStructAdapter::available() noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  return true;
#else
  return false;
#endif
}

Status HypreStructAdapter::compile(const NativeCartesianMgSpec& spec,
                                   MgRuntimeServices services,
                                   HypreStructAdapter& out) noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  Status local = validate_compile(spec, services);
  if (services.reductions != nullptr && local) {
    local = services.reductions->validate_communicator(spec.communicator);
  } else if (services.reductions == nullptr && local) {
    local = {StatusCode::invalid_plan, kHyprePlan};
  }
  Status agreed = raw_consensus(spec.communicator, local);
  if (!agreed) {
    return agreed;
  }
  agreed = services.reductions->consensus({});
  if (!agreed) {
    return agreed;
  }
  const PlanFingerprint structural = structural_fingerprint(spec);
  agreed = services.reductions->consensus_contract(structural);
  if (!agreed) {
    return agreed;
  }

  Impl* candidate = new (std::nothrow) Impl;
  local = candidate == nullptr ? Status{StatusCode::allocation_failure, 0U}
                               : Status{};
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    delete candidate;
    return agreed;
  }
  candidate->spec = spec;
  candidate->services = services;
  candidate->reductions_object = services.reductions;
  candidate->reductions_identity = services.reductions->instance_identity();
  candidate->structural = structural;
  const Int3 begin = spec.patch.begin;
  const Int3 cells = spec.patch.cells;
  candidate->lower[0] = static_cast<HYPRE_Int>(begin.x);
  candidate->lower[1] = static_cast<HYPRE_Int>(begin.y);
  candidate->lower[2] = static_cast<HYPRE_Int>(begin.z);
  candidate->upper[0] = static_cast<HYPRE_Int>(begin.x + cells.x - 1);
  candidate->upper[1] = static_cast<HYPRE_Int>(begin.y + cells.y - 1);
  candidate->upper[2] = static_cast<HYPRE_Int>(begin.z + cells.z - 1);
  local = MPI_Comm_dup(spec.communicator, &candidate->communicator) ==
                  MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure, kHypreCollective};
  if (local) {
    local = MPI_Comm_set_errhandler(candidate->communicator,
                                    MPI_ERRORS_RETURN) == MPI_SUCCESS
                ? Status{}
                : Status{StatusCode::mpi_failure, kHypreCollective};
  }
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    destroy_impl(candidate);
    return agreed;
  }
  const std::size_t cells_count = cell_count(cells);
  try {
    candidate->matrix_values.assign(cells_count * 7U, HYPRE_Complex{0.0});
    candidate->rhs_values.assign(cells_count, HYPRE_Complex{0.0});
    candidate->solution_values.assign(cells_count, HYPRE_Complex{0.0});
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kHyprePlan};
  }
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    destroy_impl(candidate);
    return agreed;
  }
  agreed = make_native_objects(*candidate);
  if (!agreed) {
    destroy_impl(candidate);
    return agreed;
  }
  out.release();
  out.implementation_ = candidate;
  return {};
#else
  (void)spec;
  (void)services;
  (void)out;
  return {StatusCode::invalid_plan, kHypreUnavailable};
#endif
}

Status HypreStructAdapter::update_coefficients(
    LinearIdentity next_identity, MgCoefficientIdentity identity,
    MgCoefficientViews coefficients) noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kHyprePlan};
  }
  Impl& implementation = *implementation_;
  Status agreed = preflight_borrowed_services(implementation);
  if (!agreed) {
    return agreed;
  }
  const bool locally_unchanged =
      implementation.ready &&
      same_identity(next_identity, implementation.spec.identity) &&
      identity.revision == implementation.coefficient_identity.revision &&
      identity.fingerprint == implementation.coefficient_identity.fingerprint;
  Status local =
      !valid_identity(next_identity) ||
              next_identity.symbolic != implementation.spec.identity.symbolic ||
              next_identity.workspace !=
                  implementation.spec.identity.workspace ||
              identity.revision == 0U || identity.fingerprint == 0U ||
              !std::isfinite(identity.maximum_relative_change) ||
              identity.maximum_relative_change < 0.0
          ? Status{StatusCode::numerical_failure, kHypreCoefficient}
          : Status{};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  agreed = implementation.services.reductions->consensus_contract(
      update_contract(next_identity, identity));
  if (!agreed) {
    return agreed;
  }
  double unchanged_local = locally_unchanged ? 1.0 : 0.0;
  double unchanged_global = 0.0;
  agreed = implementation.services.reductions->checked_sum(
      {&unchanged_local, 1U}, {&unchanged_global, 1U});
  if (!agreed) {
    return agreed;
  }
  int communicator_size = 0;
  local = MPI_Comm_size(implementation.communicator, &communicator_size) ==
                  MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure, kHypreCollective};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  if (unchanged_global == static_cast<double>(communicator_size)) {
    return {};
  }
  local = valid_coefficient_views(coefficients,
                                  implementation.spec.patch.cells)
              ? Status{}
              : Status{StatusCode::numerical_failure, kHypreCoefficient};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  double local_change = identity.maximum_relative_change;
  double global_change = 0.0;
  agreed = implementation.services.reductions->checked_max(
      {&local_change, 1U}, {&global_change, 1U});
  if (!agreed || !std::isfinite(global_change)) {
    return agreed ? Status{StatusCode::numerical_failure, kHypreCoefficient}
                  : agreed;
  }
  identity.maximum_relative_change = global_change;
  double local_anchor = 0.0;
  local = fill_matrix_values(implementation, coefficients, local_anchor);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  double global_anchor = 0.0;
  agreed = implementation.services.reductions->checked_max(
      {&local_anchor, 1U}, {&global_anchor, 1U});
  if (!agreed) {
    return agreed;
  }
  const double anchor_tolerance =
      64.0 * std::numeric_limits<double>::epsilon();
  local = implementation.spec.null_space == MgNullSpace::none
              ? (global_anchor > anchor_tolerance
                     ? Status{}
                     : Status{StatusCode::numerical_failure,
                              kHypreCoefficient})
              : (global_anchor <= anchor_tolerance
                     ? Status{}
                     : Status{StatusCode::numerical_failure,
                              kHypreCoefficient});
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }

  const std::uint8_t inactive =
      implementation.ready
          ? static_cast<std::uint8_t>(1U - implementation.active_slot)
          : implementation.active_slot;
  const bool use_reserved_pair = !implementation.pair_setup[inactive] &&
                                 implementation.matrices[inactive] != nullptr &&
                                 implementation.solvers[inactive] != nullptr;
  HYPRE_StructMatrix candidate_matrix =
      use_reserved_pair ? implementation.matrices[inactive] : nullptr;
  HYPRE_StructSolver candidate_solver =
      use_reserved_pair ? implementation.solvers[inactive] : nullptr;
  if (!use_reserved_pair) {
    local = create_native_pair(implementation, candidate_matrix,
                               candidate_solver);
    agreed = consensus(implementation, local);
    if (!agreed) {
      destroy_native_pair(candidate_matrix, candidate_solver);
      return agreed;
    }
  }
  const auto discard_candidate = [&]() noexcept {
    destroy_native_pair(candidate_matrix, candidate_solver);
    if (use_reserved_pair) {
      implementation.matrices[inactive] = nullptr;
      implementation.solvers[inactive] = nullptr;
    }
  };
  (void)HYPRE_ClearAllErrors();
  local = hypre_status(fill_matrix(candidate_matrix, implementation.lower,
                                   implementation.upper, 7,
                                   implementation.entries,
                                   implementation.matrix_values.data()),
                       kHypreCoefficient);
  agreed = consensus(implementation, local);
  if (!agreed) {
    discard_candidate();
    return agreed;
  }
  local = hypre_status(HYPRE_StructMatrixAssemble(candidate_matrix),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    discard_candidate();
    return agreed;
  }
  local = hypre_status(setup_solver(candidate_solver, candidate_matrix,
                                    implementation.rhs,
                                    implementation.solution),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    discard_candidate();
    return agreed;
  }
  if (!use_reserved_pair) {
    destroy_native_pair(implementation.matrices[inactive],
                        implementation.solvers[inactive]);
    implementation.matrices[inactive] = candidate_matrix;
    implementation.solvers[inactive] = candidate_solver;
  }
  implementation.pair_setup[inactive] = true;
  implementation.spec.identity = next_identity;
  implementation.spec.coefficients = identity;
  implementation.coefficient_identity = identity;
  implementation.numeric = make_numeric_fingerprint(
      next_identity, identity, implementation.structural);
  implementation.active_slot = inactive;
  implementation.ready = true;
  return {};
#else
  (void)next_identity;
  (void)identity;
  (void)coefficients;
  return {StatusCode::invalid_plan, kHypreUnavailable};
#endif
}

LinearPreconditionerCertificate HypreStructAdapter::certificate()
    const noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  if (implementation_ == nullptr || !implementation_->ready) {
    return {};
  }
  return {implementation_->spec.identity, implementation_->structural,
          LinearPreconditionerClass::fixed_general};
#else
  return {};
#endif
}

Status HypreStructAdapter::apply(ConstFieldView residual, FieldView correction,
                                 std::uint32_t iteration) noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  (void)iteration;
  if (implementation_ == nullptr || !implementation_->ready) {
    return {StatusCode::invalid_plan, kHypreApply};
  }
  Impl& implementation = *implementation_;
  Status agreed = preflight_borrowed_services(implementation);
  if (!agreed) {
    return agreed;
  }
  Status local =
      !valid_scalar(residual, implementation.spec.patch.cells) ||
              !valid_scalar(as_const(correction),
                            implementation.spec.patch.cells) ||
              residual.base == correction.base ||
              residual.storage_identity == correction.storage_identity
          ? Status{StatusCode::invalid_plan, kHypreApply}
          : Status{};
  std::size_t cursor = 0U;
  if (local) {
    const Int3 cells = implementation.spec.patch.cells;
    for (std::int32_t k = 0; k < cells.z && local; ++k) {
      for (std::int32_t j = 0; j < cells.y && local; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const double value = residual.unchecked({i, j, k}, 0U);
          if (!std::isfinite(value)) {
            local = {StatusCode::numerical_failure, kHypreApply};
            break;
          }
          implementation.rhs_values[cursor++] =
              static_cast<HYPRE_Complex>(value);
        }
      }
    }
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  if (implementation.spec.null_space == MgNullSpace::constant) {
    agreed = project_values(implementation, implementation.rhs_values);
    if (!agreed) {
      return agreed;
    }
  }
  std::fill(implementation.solution_values.begin(),
            implementation.solution_values.end(), HYPRE_Complex{0.0});
  (void)HYPRE_ClearAllErrors();
  local = hypre_status(HYPRE_StructVectorSetBoxValues(
                           implementation.rhs, implementation.lower,
                           implementation.upper,
                           implementation.rhs_values.data()),
                       kHypreApply);
  if (local) {
    local = hypre_status(HYPRE_StructVectorSetBoxValues(
                             implementation.solution, implementation.lower,
                             implementation.upper,
                             implementation.solution_values.data()),
                         kHypreApply);
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructVectorAssemble(implementation.rhs),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructVectorAssemble(implementation.solution),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructPFMGSolve(
                           implementation.solvers[implementation.active_slot],
                           implementation.matrices[implementation.active_slot],
                           implementation.rhs, implementation.solution),
                       kHypreCollective);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  local = hypre_status(HYPRE_StructVectorGetBoxValues(
                           implementation.solution, implementation.lower,
                           implementation.upper,
                           implementation.solution_values.data()),
                       kHypreApply);
  if (local) {
    for (HYPRE_Complex value : implementation.solution_values) {
      if (!std::isfinite(static_cast<double>(value))) {
        local = {StatusCode::numerical_failure, kHypreApply};
        break;
      }
    }
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  if (implementation.spec.null_space == MgNullSpace::constant) {
    agreed = project_values(implementation, implementation.solution_values);
    if (!agreed) {
      return agreed;
    }
  }
  cursor = 0U;
  const Int3 cells = implementation.spec.patch.cells;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        correction.unchecked({i, j, k}, 0U) =
            static_cast<double>(implementation.solution_values[cursor++]);
      }
    }
  }
  return {};
#else
  (void)residual;
  (void)correction;
  (void)iteration;
  return {StatusCode::invalid_plan, kHypreUnavailable};
#endif
}

PlanFingerprint HypreStructAdapter::numeric_fingerprint() const noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  return implementation_ == nullptr ? 0U : implementation_->numeric;
#else
  return 0U;
#endif
}

PlanFingerprint HypreStructAdapter::fingerprint() const noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  return implementation_ == nullptr ? 0U : implementation_->structural;
#else
  return 0U;
#endif
}

std::uintptr_t HypreStructAdapter::native_handle_storage_address()
    const noexcept {
#if defined(HUNDUN_V04_HAVE_HYPRE) && HUNDUN_V04_HAVE_HYPRE
  return implementation_ == nullptr
             ? 0U
             : reinterpret_cast<std::uintptr_t>(&implementation_->grid);
#else
  return 0U;
#endif
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

HypreLifecycleSnapshot hypre_lifecycle_snapshot_for_test() noexcept {
  const HypreLifecycleCounters& counters = hypre_lifecycle_counters();
  return {counters.matrix_creates.load(std::memory_order_relaxed),
          counters.matrix_destroys.load(std::memory_order_relaxed),
          counters.solver_creates.load(std::memory_order_relaxed),
          counters.solver_destroys.load(std::memory_order_relaxed),
          counters.solver_setups.load(std::memory_order_relaxed)};
}

void reset_hypre_lifecycle_for_test() noexcept {
  HypreLifecycleCounters& counters = hypre_lifecycle_counters();
  counters.matrix_creates.store(0U, std::memory_order_relaxed);
  counters.matrix_destroys.store(0U, std::memory_order_relaxed);
  counters.solver_creates.store(0U, std::memory_order_relaxed);
  counters.solver_destroys.store(0U, std::memory_order_relaxed);
  counters.solver_setups.store(0U, std::memory_order_relaxed);
  clear_hypre_failure_for_test();
}

void set_hypre_failure_for_test(HypreFailurePoint point) noexcept {
  g_hypre_failure_point.store(point, std::memory_order_relaxed);
}

void clear_hypre_failure_for_test() noexcept {
  g_hypre_failure_point.store(HypreFailurePoint::none,
                              std::memory_order_relaxed);
}

}  // namespace detail
#endif

}  // namespace hundun::v04
