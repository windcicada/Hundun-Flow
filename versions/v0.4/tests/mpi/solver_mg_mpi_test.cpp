// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_parallel.hpp"
#include "parallel_halo_detail.hpp"
#include "solver_mg_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <type_traits>
#include <vector>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  void* const result = std::malloc(size == 0U ? 1U : size);
  if (result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* result = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&result, alignment, requested) != 0 || result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

constexpr RevisionToken kExecutionRevision = 12001U;
constexpr FieldId kDiagonal = 20U;
constexpr FieldId kResidual = 21U;
constexpr FieldId kCorrection = 22U;
constexpr FieldId kWorkspaceVectors = 30U;

static_assert(!std::is_copy_constructible_v<NativeCartesianMgPlan>);
static_assert(!std::is_copy_assignable_v<NativeCartesianMgPlan>);
static_assert(std::is_nothrow_move_constructible_v<NativeCartesianMgPlan>);
static_assert(std::is_nothrow_move_assignable_v<NativeCartesianMgPlan>);

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local, MPI_Comm communicator = MPI_COMM_WORLD) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         output != 0;
}

std::uint64_t packed(Status status) noexcept {
  return (static_cast<std::uint64_t>(status.code) << 32U) |
         static_cast<std::uint64_t>(status.detail);
}

bool identical(std::uint64_t value,
               MPI_Comm communicator = MPI_COMM_WORLD) noexcept {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(MgPlanCounters left, MgPlanCounters right) noexcept {
  return left.symbolic_builds == right.symbolic_builds &&
         left.numeric_refreshes == right.numeric_refreshes &&
         left.hierarchy_rebuilds == right.hierarchy_rebuilds &&
         left.applications == right.applications &&
         left.blocking_collectives == right.blocking_collectives &&
         left.collective_logical_bytes == right.collective_logical_bytes &&
         left.point_to_point_messages == right.point_to_point_messages &&
         left.point_to_point_bytes == right.point_to_point_bytes;
}

bool same(detail::MgMatrixWorkCounters left,
          detail::MgMatrixWorkCounters right) noexcept {
  return left.cycle_level_calls == right.cycle_level_calls &&
         left.cycle_restrictions == right.cycle_restrictions &&
         left.cycle_prolongations == right.cycle_prolongations &&
         left.full_actions == right.full_actions &&
         left.full_cell_visits == right.full_cell_visits &&
         left.residual_writes == right.residual_writes &&
         left.retained_final_full_actions ==
             right.retained_final_full_actions &&
         left.retained_final_defect_actions ==
             right.retained_final_defect_actions &&
         left.defect_cell_visits == right.defect_cell_visits &&
         left.defect_writes == right.defect_writes &&
         left.residual_finish_actions == right.residual_finish_actions &&
         left.residual_finish_cell_visits ==
             right.residual_finish_cell_visits &&
         left.fused_color_actions == right.fused_color_actions &&
         left.fused_cell_visits == right.fused_cell_visits &&
         left.fused_updates == right.fused_updates &&
         left.separated_color_actions == right.separated_color_actions &&
         left.separated_color_updates == right.separated_color_updates &&
         left.chebyshev_stages == right.chebyshev_stages &&
         left.chebyshev_exchange_actions ==
             right.chebyshev_exchange_actions &&
         left.chebyshev_defect_actions == right.chebyshev_defect_actions &&
         left.chebyshev_retained_final_defect_actions ==
             right.chebyshev_retained_final_defect_actions &&
         left.chebyshev_stencil_evaluations ==
             right.chebyshev_stencil_evaluations &&
         left.chebyshev_updates == right.chebyshev_updates &&
         left.chebyshev_elided_intermediate_residual_publications ==
             right.chebyshev_elided_intermediate_residual_publications &&
         left.chebyshev_copyback_cells == right.chebyshev_copyback_cells;
}

bool same_chebyshev_scientific_work(
    detail::MgMatrixWorkCounters left,
    detail::MgMatrixWorkCounters right) noexcept {
  // These are the only intentionally different physical lifecycle fields.
  // Normalize them solely for the Chebyshev oracle; the generic comparator
  // above remains strict for every other test.
  left.defect_writes = 0U;
  right.defect_writes = 0U;
  left.chebyshev_elided_intermediate_residual_publications = 0U;
  right.chebyshev_elided_intermediate_residual_publications = 0U;
  left.chebyshev_copyback_cells = 0U;
  right.chebyshev_copyback_cells = 0U;
  return same(left, right);
}

bool exact_cycle_schedule(const detail::MgMatrixWorkCounters& work,
                          std::size_t levels, MgCycleKind cycle) noexcept {
  std::uint64_t transfers = 0U;
  for (std::size_t level = 0U; level < detail::kMgMaximumLevels; ++level) {
    const std::uint64_t expected =
        level < levels
            ? (cycle == MgCycleKind::f_cycle
                   ? static_cast<std::uint64_t>(level + 1U)
                   : 1U)
            : 0U;
    if (work.cycle_level_calls[level] != expected) return false;
    if (level + 1U < levels) transfers += expected;
  }
  return work.cycle_restrictions == transfers &&
         work.cycle_prolongations == transfers;
}

bool greater_or_equal(LinearReductionCounters after,
                      LinearReductionCounters before) noexcept {
  return after.calls >= before.calls && after.scalars >= before.scalars &&
         after.logical_bytes >= before.logical_bytes &&
         after.tree_messages >= before.tree_messages;
}

LinearReductionCounters difference(LinearReductionCounters after,
                                   LinearReductionCounters before) noexcept {
  return {after.calls - before.calls, after.scalars - before.scalars,
          after.logical_bytes - before.logical_bytes,
          after.tree_messages - before.tree_messages};
}

HaloRuntimeCounters difference(HaloRuntimeCounters after,
                               HaloRuntimeCounters before) noexcept {
  return {after.begin_calls - before.begin_calls,
          after.finish_calls - before.finish_calls,
          after.messages_started - before.messages_started,
          after.bytes_packed - before.bytes_packed,
          after.bytes_unpacked - before.bytes_unpacked,
          after.control_consensus_calls - before.control_consensus_calls};
}

bool same(LinearReductionCounters left,
          LinearReductionCounters right) noexcept {
  return left.calls == right.calls && left.scalars == right.scalars &&
         left.logical_bytes == right.logical_bytes &&
         left.tree_messages == right.tree_messages;
}

bool same(HaloRuntimeCounters left, HaloRuntimeCounters right) noexcept {
  return left.begin_calls == right.begin_calls &&
         left.finish_calls == right.finish_calls &&
         left.messages_started == right.messages_started &&
         left.bytes_packed == right.bytes_packed &&
         left.bytes_unpacked == right.bytes_unpacked;
}

bool bitwise_same(double left, double right) noexcept {
  return std::memcmp(&left, &right, sizeof(left)) == 0;
}

CartesianMeshSpec mesh_spec(Int3 cells = {24, 12, 8}) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {1.0e-10, 1.0e-10, 1.0e-10};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

LinearIdentity identity(std::uint64_t salt = 0U) noexcept {
  return {101U + salt, 102U + salt, 103U + salt, 104U + salt,
          105U + salt};
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t ghosts,
                      std::uint8_t components, RevisionToken revision,
                      StorageIdentity storage,
                      RevisionDomainIdentity domain) {
  OwnedField result;
  const std::size_t stride_y =
      static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t stride_z =
      stride_y * static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t component_stride =
      stride_z * static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(component_stride * components, 0.0);
  result.view.base =
      result.storage.data() + ghosts + static_cast<std::size_t>(ghosts) *
                                             (stride_y + stride_z);
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = stride_y;
  result.view.stride_z = stride_z;
  result.view.component_stride = component_stride;
  result.view.field = id;
  result.view.revision = revision;
  result.view.storage_identity = storage;
  result.view.revision_domain = domain;
  return result;
}

struct OwnedFaceField {
  std::vector<double> storage;
  FaceFieldView view{};
};

OwnedFaceField make_face(CartesianAxis axis, Int3 cells,
                         StorageIdentity storage,
                         RevisionDomainIdentity domain) {
  OwnedFaceField result;
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
                 storage, domain};
  return result;
}

void fill(OwnedField& field, double value) noexcept {
  for (std::int32_t z = 0; z < field.view.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.view.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.view.interior.x; ++x) {
        field.view.unchecked({x, y, z}, 0U) = value;
      }
    }
  }
}

bool finite(const FieldView& field) noexcept {
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        if (!std::isfinite(field.unchecked({x, y, z}, 0U))) {
          return false;
        }
      }
    }
  }
  return true;
}

std::uint64_t checksum(const FieldView& field) noexcept {
  std::uint64_t result = UINT64_C(1469598103934665603);
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        std::uint64_t bits = 0U;
        const double value = field.unchecked({x, y, z}, 0U);
        std::memcpy(&bits, &value, sizeof(bits));
        result ^= bits;
        result *= UINT64_C(1099511628211);
      }
    }
  }
  return result;
}

struct Fixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  OwnedField diagonal;
  OwnedFaceField x;
  OwnedFaceField y;
  OwnedFaceField z;
  OwnedField residual;
  OwnedField correction;
  OwnedField vectors;
  OwnedField solver_vectors;
  OwnedField solver_scalars;
  MgWorkspaceRequirements requirements{};
  MgWorkspace workspace;
  LinearWorkspaceRequirements solver_requirements{};
  SolverWorkspace solver_workspace;
  ReductionEngine reductions;
  HaloEngine halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  double anisotropy_threshold{4.0};
  std::uint8_t maximum_levels{32U};
  MgCycleKind cycle{MgCycleKind::v_cycle};
};

bool initialize(Fixture& fixture, HaloTopology topology = {},
                double anisotropy_threshold = 4.0,
                Int3 exact_cells = {24, 12, 8},
                std::uint8_t maximum_levels = 32U,
                MgCycleKind cycle = MgCycleKind::v_cycle) {
  fixture.anisotropy_threshold = anisotropy_threshold;
  fixture.maximum_levels = maximum_levels;
  fixture.cycle = cycle;
  if (!CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, mesh_spec(exact_cells), GeometryBudget{},
          fixture.geometry,
          fixture.patch)) {
    return false;
  }
  constexpr StorageIdentity field_storage = 31001U;
  constexpr RevisionDomainIdentity domain = 32001U;
  fixture.diagonal = make_field(kDiagonal, fixture.patch.cells, 0U, 1U, 11U,
                                field_storage, domain);
  fixture.x = make_face(CartesianAxis::x, fixture.patch.cells, field_storage,
                        domain);
  fixture.y = make_face(CartesianAxis::y, fixture.patch.cells, field_storage,
                        domain);
  fixture.z = make_face(CartesianAxis::z, fixture.patch.cells, field_storage,
                        domain);
  fixture.residual = make_field(kResidual, fixture.patch.cells, 0U, 1U, 12U,
                                field_storage + 1U, domain);
  fixture.correction =
      make_field(kCorrection, fixture.patch.cells, 0U, 1U, 13U,
                 field_storage + 2U, domain);
  fill(fixture.diagonal, 6.0);
  fill(fixture.residual, 1.0);
  fill(fixture.correction, -7.0);

  NativeCartesianMgSpec workspace_spec;
  workspace_spec.policy.anisotropy_threshold = anisotropy_threshold;
  workspace_spec.policy.coefficient_change_rebuild_ratio = 0.25;
  workspace_spec.policy.maximum_levels = maximum_levels;
  workspace_spec.policy.cycle = cycle;
  if (!make_mg_workspace_requirements(
          MPI_COMM_WORLD, fixture.geometry, fixture.patch,
          workspace_spec.policy, kExecutionRevision, fixture.requirements)) {
    return false;
  }
  constexpr StorageIdentity workspace_storage = 33001U;
  fixture.vectors = make_field(
      kWorkspaceVectors, fixture.requirements.arena_shape, 0U, 1U, 15U,
      workspace_storage, domain);
  const std::array halo_specs{
      HaloFieldSpec{kWorkspaceVectors, 1U, 1U}};
  if (!MgWorkspace::bind(fixture.requirements, fixture.vectors.view,
                         fixture.workspace) ||
      !ReductionEngine::compile(MPI_COMM_WORLD,
                                ReductionMode::mpi_allreduce, 4U,
                                fixture.reductions) ||
      !fixture.halo.reserve(
          MPI_COMM_WORLD, fixture.requirements.levels[0U].patch,
          Span<const HaloFieldSpec>{halo_specs.data(), halo_specs.size()},
          topology)) {
    return false;
  }
  if (!make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, fixture.patch.cells, 0U, 4U,
          ReductionMode::mpi_allreduce, kExecutionRevision,
          fixture.solver_requirements)) {
    return false;
  }
  fixture.solver_vectors = make_field(
      kWorkspaceVectors + 1U, fixture.patch.cells, 0U,
      fixture.solver_requirements.vector_slots, 16U, workspace_storage + 1U,
      domain);
  fixture.solver_scalars = make_field(
      kWorkspaceVectors + 2U,
      {static_cast<std::int32_t>(fixture.solver_requirements.scalar_doubles),
       1, 1},
      0U, 1U, 17U, workspace_storage + 1U, domain);
  if (!SolverWorkspace::bind(fixture.solver_requirements,
                              fixture.solver_vectors.view,
                              fixture.solver_scalars.view,
                              fixture.solver_workspace)) {
    return false;
  }
  fixture.coarse_halos.resize(fixture.requirements.level_count - 1U);
  fixture.coarse_halo_pointers.resize(fixture.coarse_halos.size());
  for (std::size_t level = 1U; level < fixture.requirements.level_count;
       ++level) {
    if (!fixture.coarse_halos[level - 1U].reserve(
            MPI_COMM_WORLD, fixture.requirements.levels[level].patch,
            Span<const HaloFieldSpec>{halo_specs.data(), halo_specs.size()},
            topology)) {
      return false;
    }
    fixture.coarse_halo_pointers[level - 1U] =
        &fixture.coarse_halos[level - 1U];
  }
  return true;
}

MgCoefficientViews coefficient_views(const Fixture& fixture) noexcept {
  return {as_const(fixture.diagonal.view),
          ConstFaceFieldView{fixture.x.view.base, fixture.x.view.extents,
                             fixture.x.view.stride_y, fixture.x.view.stride_z,
                             fixture.x.view.axis,
                             fixture.x.view.storage_identity,
                             fixture.x.view.revision_domain},
          ConstFaceFieldView{fixture.y.view.base, fixture.y.view.extents,
                             fixture.y.view.stride_y, fixture.y.view.stride_z,
                             fixture.y.view.axis,
                             fixture.y.view.storage_identity,
                             fixture.y.view.revision_domain},
          ConstFaceFieldView{fixture.z.view.base, fixture.z.view.extents,
                             fixture.z.view.stride_y, fixture.z.view.stride_z,
                             fixture.z.view.axis,
                             fixture.z.view.storage_identity,
                             fixture.z.view.revision_domain}};
}

NativeCartesianMgSpec mg_spec(
    const Fixture& fixture, std::uint64_t salt = 0U,
    MgPointSmootherKind smoother = MgPointSmootherKind::red_black,
    MgOperatorClass operator_class = MgOperatorClass::general) noexcept {
  NativeCartesianMgSpec result;
  result.communicator = MPI_COMM_WORLD;
  result.geometry = &fixture.geometry;
  result.patch = fixture.patch;
  result.identity = identity(salt);
  result.coefficients = {41U + salt, 51U + salt, 0.0};
  result.policy.anisotropy_threshold = fixture.anisotropy_threshold;
  result.policy.coefficient_change_rebuild_ratio = 0.25;
  result.policy.maximum_levels = fixture.maximum_levels;
  result.policy.point_smoother = smoother;
  result.policy.cycle = fixture.cycle;
  result.operator_class = operator_class;
  return result;
}

MgRuntimeServices services(Fixture& fixture) noexcept {
  return {&fixture.halo, &fixture.reductions, &fixture.workspace,
          {fixture.coarse_halo_pointers.data(),
           fixture.coarse_halo_pointers.size()}};
}

HaloRuntimeCounters halo_counters(const Fixture& fixture) noexcept {
  HaloRuntimeCounters result = fixture.halo.runtime_counters();
  for (const HaloEngine& halo : fixture.coarse_halos) {
    const HaloRuntimeCounters level = halo.runtime_counters();
    result.begin_calls += level.begin_calls;
    result.finish_calls += level.finish_calls;
    result.messages_started += level.messages_started;
    result.bytes_packed += level.bytes_packed;
    result.bytes_unpacked += level.bytes_unpacked;
    result.control_consensus_calls += level.control_consensus_calls;
  }
  return result;
}

MgBoundarySet boundaries(HaloTopology topology) noexcept;
void fill_variable_point_fixture(Fixture& fixture) noexcept;

bool test_periodic_prolongation_wraps_coarse_z_endpoint(int rank, int size) {
  constexpr Int3 cells{8, 8, 32};
  constexpr HaloTopology periodic_z{false, false, true};
  Fixture periodic_fixture;
  Fixture physical_fixture;
  bool passed =
      expect(initialize(periodic_fixture, periodic_z, 100.0, cells, 2U) &&
                 initialize(physical_fixture, {}, 100.0, cells, 2U),
             rank, "periodic/physical prolongation fixtures initialize");
  if (!all_true(passed)) return false;

  NativeCartesianMgSpec periodic_spec = mg_spec(periodic_fixture, 905U);
  periodic_spec.boundaries = boundaries(periodic_z);
  NativeCartesianMgSpec physical_spec = mg_spec(physical_fixture, 906U);
  physical_spec.boundaries.z_min = MgBoundaryKind::neumann;
  physical_spec.boundaries.z_max = MgBoundaryKind::neumann;
  NativeCartesianMgPlan periodic_plan;
  NativeCartesianMgPlan physical_plan;
  const Status periodic_compile = NativeCartesianMgPlan::compile(
      periodic_spec, services(periodic_fixture),
      coefficient_views(periodic_fixture), periodic_plan);
  const Status physical_compile = NativeCartesianMgPlan::compile(
      physical_spec, services(physical_fixture),
      coefficient_views(physical_fixture), physical_plan);
  const bool fixture_contract =
      static_cast<bool>(periodic_compile) &&
      static_cast<bool>(physical_compile) &&
      periodic_plan.level_count() == 2U && physical_plan.level_count() == 2U &&
      (size == 1 ? periodic_fixture.patch.process_grid.z == 1
                 : periodic_fixture.patch.process_grid.z > 1) &&
      (size == 1 ? physical_fixture.patch.process_grid.z == 1
                 : physical_fixture.patch.process_grid.z > 1) &&
      periodic_fixture.requirements.levels[1U].global_shape.z == 16 &&
      physical_fixture.requirements.levels[1U].global_shape.z == 16;
  if (!fixture_contract) {
    std::cerr << "rank " << rank << " prolongation fixture diagnostics: "
              << "periodic_compile=" << packed(periodic_compile)
              << " physical_compile=" << packed(physical_compile)
              << " periodic_levels=" << periodic_plan.level_count()
              << " physical_levels=" << physical_plan.level_count()
              << " periodic_grid_z=" << periodic_fixture.patch.process_grid.z
              << " physical_grid_z=" << physical_fixture.patch.process_grid.z
              << " periodic_coarse_z="
              << periodic_fixture.requirements.levels[1U].global_shape.z
              << " physical_coarse_z="
              << physical_fixture.requirements.levels[1U].global_shape.z
              << '\n';
  }
  passed &= expect(
      fixture_contract, rank,
      "prolongation fixture puts periodic z endpoint across the MPI partition");
  if (!all_true(passed)) return false;

  const auto seed_coarse_and_clear_fine = [](Fixture& fixture) noexcept {
    FieldView coarse = fixture.workspace.level(1U, MgWorkspaceSlot::solution);
    const MeshPatch coarse_patch = fixture.requirements.levels[1U].patch;
    for (std::int32_t k = 0; k < coarse.interior.z; ++k) {
      const double value = static_cast<double>(coarse_patch.begin.z + k);
      for (std::int32_t j = 0; j < coarse.interior.y; ++j) {
        for (std::int32_t i = 0; i < coarse.interior.x; ++i) {
          coarse.unchecked({i, j, k}, 0U) = value;
        }
      }
    }
    FieldView fine = fixture.workspace.level(0U, MgWorkspaceSlot::solution);
    for (std::int32_t k = 0; k < fine.interior.z; ++k) {
      for (std::int32_t j = 0; j < fine.interior.y; ++j) {
        for (std::int32_t i = 0; i < fine.interior.x; ++i) {
          fine.unchecked({i, j, k}, 0U) = 0.0;
        }
      }
    }
  };
  seed_coarse_and_clear_fine(periodic_fixture);
  seed_coarse_and_clear_fine(physical_fixture);

  const Status periodic_status =
      detail::mg_prolongate_add_for_test(periodic_plan, 0U, 12050U);
  const Status physical_status =
      detail::mg_prolongate_add_for_test(physical_plan, 0U, 12051U);
  passed &= expect(static_cast<bool>(periodic_status) &&
                       static_cast<bool>(physical_status),
                   rank, "production prolongation succeeds");
  if (!all_true(passed)) return false;

  const auto endpoint_values_are = [](const Fixture& fixture,
                                      double expected_minimum,
                                      double expected_maximum) noexcept {
    const FieldView fine =
        fixture.workspace.level(0U, MgWorkspaceSlot::solution);
    bool local = true;
    for (std::int32_t k = 0; k < fine.interior.z; ++k) {
      const std::int32_t global_z = fixture.patch.begin.z + k;
      if (global_z != 0 && global_z != cells.z - 1) continue;
      const double expected =
          global_z == 0 ? expected_minimum : expected_maximum;
      for (std::int32_t j = 0; j < fine.interior.y; ++j) {
        for (std::int32_t i = 0; i < fine.interior.x; ++i) {
          local = local && fine.unchecked({i, j, k}, 0U) == expected;
        }
      }
    }
    return local;
  };
  // With q_c(k)=k on 16 coarse cells, cell-centred linear interpolation at
  // periodic fine endpoints is 3/4*q_c(0)+1/4*q_c(15)=3.75 and
  // 3/4*q_c(15)+1/4*q_c(0)=11.25.  A physical symmetry boundary clamps.
  passed &= expect(endpoint_values_are(periodic_fixture, 3.75, 11.25), rank,
                   "periodic z prolongation consumes wrapped coarse neighbor");
  passed &= expect(endpoint_values_are(physical_fixture, 0.0, 15.0), rank,
                   "physical z prolongation keeps endpoint clamp");
  return all_true(passed);
}

bool test_cross_partition_periodic_prolongation_uses_complete_donors(
    int rank, int size) {
  constexpr Int3 cells{32, 8, 32};
  constexpr HaloTopology periodic_z{false, false, true};
  constexpr double ghost_sentinel = -777.0;
  Fixture periodic_fixture;
  Fixture physical_fixture;
  bool passed = expect(
      initialize(periodic_fixture, periodic_z, 100.0, cells, 2U) &&
          initialize(physical_fixture, {}, 100.0, cells, 2U),
      rank, "cross-partition prolongation fixtures initialize");
  if (!all_true(passed)) return false;

  NativeCartesianMgSpec periodic_spec = mg_spec(periodic_fixture, 907U);
  periodic_spec.boundaries = boundaries(periodic_z);
  NativeCartesianMgSpec physical_spec = mg_spec(physical_fixture, 908U);
  physical_spec.boundaries.z_min = MgBoundaryKind::neumann;
  physical_spec.boundaries.z_max = MgBoundaryKind::neumann;
  NativeCartesianMgPlan periodic_plan;
  NativeCartesianMgPlan physical_plan;
  const Status periodic_compile = NativeCartesianMgPlan::compile(
      periodic_spec, services(periodic_fixture),
      coefficient_views(periodic_fixture), periodic_plan);
  const Status physical_compile = NativeCartesianMgPlan::compile(
      physical_spec, services(physical_fixture),
      coefficient_views(physical_fixture), physical_plan);
  const Int3 expected_grid = size == 4   ? Int3{2, 1, 2}
                             : size == 2 ? Int3{1, 1, 2}
                                       : Int3{1, 1, 1};
  const auto same_grid = [expected_grid](Int3 grid) noexcept {
    return grid.x == expected_grid.x && grid.y == expected_grid.y &&
           grid.z == expected_grid.z;
  };
  passed &= expect(
      static_cast<bool>(periodic_compile) &&
          static_cast<bool>(physical_compile) &&
          periodic_plan.level_count() == 2U &&
          physical_plan.level_count() == 2U &&
          same_grid(periodic_fixture.patch.process_grid) &&
          same_grid(physical_fixture.patch.process_grid) &&
          periodic_fixture.requirements.levels[1U].global_shape.x == 16 &&
          periodic_fixture.requirements.levels[1U].global_shape.z == 16 &&
          physical_fixture.requirements.levels[1U].global_shape.x == 16 &&
          physical_fixture.requirements.levels[1U].global_shape.z == 16,
      rank,
      "four-rank cross-partition fixture decomposes exactly two-by-one-by-two");
  if (!all_true(passed)) return false;

  const auto coarse_value = [](std::int32_t global_x,
                               std::int32_t global_z) noexcept {
    const double x = static_cast<double>(global_x);
    const double z = static_cast<double>(global_z);
    return 100.0 * x * z + x + z;
  };
  const auto seed = [coarse_value](Fixture& fixture) noexcept {
    FieldView coarse =
        fixture.workspace.level(1U, MgWorkspaceSlot::solution);
    for (std::int32_t k = -1; k <= coarse.interior.z; ++k) {
      for (std::int32_t j = -1; j <= coarse.interior.y; ++j) {
        for (std::int32_t i = -1; i <= coarse.interior.x; ++i) {
          coarse.unchecked({i, j, k}, 0U) = ghost_sentinel;
        }
      }
    }
    const MeshPatch coarse_patch = fixture.requirements.levels[1U].patch;
    for (std::int32_t k = 0; k < coarse.interior.z; ++k) {
      for (std::int32_t j = 0; j < coarse.interior.y; ++j) {
        for (std::int32_t i = 0; i < coarse.interior.x; ++i) {
          coarse.unchecked({i, j, k}, 0U) = coarse_value(
              coarse_patch.begin.x + i, coarse_patch.begin.z + k);
        }
      }
    }
    FieldView fine = fixture.workspace.level(0U, MgWorkspaceSlot::solution);
    for (std::int32_t k = 0; k < fine.interior.z; ++k) {
      for (std::int32_t j = 0; j < fine.interior.y; ++j) {
        for (std::int32_t i = 0; i < fine.interior.x; ++i) {
          fine.unchecked({i, j, k}, 0U) = 0.0;
        }
      }
    }
  };
  seed(periodic_fixture);
  seed(physical_fixture);

  const MgPlanCounters periodic_before = periodic_plan.counters();
  const Status periodic_status = detail::mg_prolongate_add_for_test(
      periodic_plan, 0U, 12052U);
  const Status physical_status = detail::mg_prolongate_add_for_test(
      physical_plan, 0U, 12053U);
  passed &= expect(static_cast<bool>(periodic_status) &&
                       static_cast<bool>(physical_status),
                   rank, "cross-partition production prolongation succeeds");
  if (!all_true(passed)) return false;
  const MgPlanCounters periodic_after = periodic_plan.counters();
  const std::uint64_t expected_messages = size == 4 ? 2U : 0U;
  const std::uint64_t expected_bytes = size == 4 ? 1280U : 0U;
  passed &= expect(
      periodic_after.point_to_point_messages -
                  periodic_before.point_to_point_messages ==
              expected_messages &&
          periodic_after.point_to_point_bytes -
                  periodic_before.point_to_point_bytes ==
              expected_bytes,
      rank, "staged prolongation publishes exact point-to-point resources");

  bool local_corner_complete = true;
  bool local_witness_correct = true;
  if (size == 4 && periodic_fixture.patch.process_coord.x == 1 &&
      periodic_fixture.patch.process_coord.z == 0) {
    const FieldView coarse =
        periodic_fixture.workspace.level(1U, MgWorkspaceSlot::solution);
    const FieldView fine =
        periodic_fixture.workspace.level(0U, MgWorkspaceSlot::solution);
    const double diagonal_corner = coarse.unchecked({-1, 0, -1}, 0U);
    const double witness = fine.unchecked({0, 0, 0}, 0U);
    local_corner_complete = diagonal_corner == 10522.0;
    local_witness_correct = witness == 2917.75;
    if (!local_witness_correct) {
      std::cerr << "rank " << rank
                << " cross-partition periodic prolongation witness="
                << witness << " expected=2917.75 corner=" << diagonal_corner
                << '\n';
    }
  }
  passed &= expect(local_corner_complete, rank,
                   "diagonal donor authority fills the x-z corner ghost");
  passed &= expect(local_witness_correct, rank,
                   "cross-partition periodic endpoint ignores stale corner ghost");

  const auto interpolated_coordinate = [](std::int32_t fine,
                                          std::int32_t coarse_extent,
                                          bool periodic) noexcept {
    if (fine == 0) {
      return periodic ? 0.25 * static_cast<double>(coarse_extent - 1) : 0.0;
    }
    if (fine == 2 * coarse_extent - 1) {
      return periodic ? 0.75 * static_cast<double>(coarse_extent - 1)
                      : static_cast<double>(coarse_extent - 1);
    }
    return 0.5 * static_cast<double>(fine) - 0.25;
  };
  const auto matches_analytic = [interpolated_coordinate](
                                    const Fixture& fixture,
                                    bool periodic) noexcept {
    const FieldView fine =
        fixture.workspace.level(0U, MgWorkspaceSlot::solution);
    bool local = true;
    for (std::int32_t k = 0; k < fine.interior.z; ++k) {
      const std::int32_t global_z = fixture.patch.begin.z + k;
      const double z = interpolated_coordinate(global_z, 16, periodic);
      for (std::int32_t j = 0; j < fine.interior.y; ++j) {
        for (std::int32_t i = 0; i < fine.interior.x; ++i) {
          const std::int32_t global_x = fixture.patch.begin.x + i;
          const double x = interpolated_coordinate(global_x, 16, false);
          const double expected = 100.0 * x * z + x + z;
          const double actual = fine.unchecked({i, j, k}, 0U);
          local = local &&
                  std::abs(actual - expected) <=
                      1.0e-13 * std::max(1.0, std::abs(expected));
        }
      }
    }
    return local;
  };
  passed &= expect(matches_analytic(periodic_fixture, true), rank,
                   "periodic cross term is decomposition-independent");
  passed &= expect(matches_analytic(physical_fixture, false), rank,
                   "physical z endpoint remains clamped with x partitioning");

  const auto conservative = [](const Fixture& fixture) noexcept {
    const FieldView coarse =
        fixture.workspace.level(1U, MgWorkspaceSlot::solution);
    const FieldView fine =
        fixture.workspace.level(0U, MgWorkspaceSlot::solution);
    std::array<double, 2U> local{};
    for (std::int32_t k = 0; k < coarse.interior.z; ++k)
      for (std::int32_t j = 0; j < coarse.interior.y; ++j)
        for (std::int32_t i = 0; i < coarse.interior.x; ++i)
          local[0U] += coarse.unchecked({i, j, k}, 0U);
    for (std::int32_t k = 0; k < fine.interior.z; ++k)
      for (std::int32_t j = 0; j < fine.interior.y; ++j)
        for (std::int32_t i = 0; i < fine.interior.x; ++i)
          local[1U] += fine.unchecked({i, j, k}, 0U);
    std::array<double, 2U> global{};
    if (MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS) {
      return false;
    }
    return std::abs(global[1U] - 8.0 * global[0U]) <=
           1.0e-13 * std::max(1.0, std::abs(global[1U]));
  };
  passed &= expect(conservative(periodic_fixture) &&
                       conservative(physical_fixture),
                   rank,
                   "combined prolongation preserves the global integrated sum");
  return all_true(passed);
}

bool test_fused_point_smoother_reference_oracle(int rank,
                                                bool zero_pre_sweeps = false) {
  // Odd local extents exercise both color parities without creating a
  // periodic wrap.  The seven-point stencil remains bipartite here.
  Fixture optimized_fixture;
  Fixture reference_fixture;
  bool passed = expect(
      initialize(optimized_fixture, {}, 4.0, {25, 11, 9}) &&
          initialize(reference_fixture, {}, 4.0, {25, 11, 9}),
      rank, "fused/reference MG fixtures initialize");
  if (!all_true(passed)) return false;

  NativeCartesianMgPlan optimized;
  NativeCartesianMgPlan reference;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(optimized_fixture, 910U), services(optimized_fixture),
          coefficient_views(optimized_fixture), optimized)) &&
          static_cast<bool>(NativeCartesianMgPlan::compile(
              mg_spec(reference_fixture, 910U), services(reference_fixture),
              coefficient_views(reference_fixture), reference)) &&
          optimized.line_axis_mask() == 0U && reference.line_axis_mask() == 0U,
      rank, "fused and separated point plans compile");
  if (!all_true(passed)) return false;

  detail::set_mg_reference_point_actions_for_test(optimized, false);
  detail::set_mg_reference_point_actions_for_test(reference, true);
  if (zero_pre_sweeps) {
    // The public policy rejects zero pre-sweeps.  This test-only override
    // reaches the already compiled cycle to exercise the internal fail-safe
    // branch without changing that production contract.
    detail::set_mg_pre_sweeps_for_test(optimized, 0U);
    detail::set_mg_pre_sweeps_for_test(reference, 0U);
  }
  const HaloRuntimeCounters optimized_halo_before =
      halo_counters(optimized_fixture);
  const HaloRuntimeCounters reference_halo_before =
      halo_counters(reference_fixture);
  const LinearReductionCounters optimized_reductions_before =
      optimized_fixture.reductions.counters();
  const LinearReductionCounters reference_reductions_before =
      reference_fixture.reductions.counters();

  const Status optimized_status = optimized.apply(
      as_const(optimized_fixture.residual.view),
      optimized_fixture.correction.view, 0U);
  const Status reference_status = reference.apply(
      as_const(reference_fixture.residual.view),
      reference_fixture.correction.view, 0U);
  const detail::MgMatrixWorkCounters optimized_work =
      detail::mg_matrix_work_counters_for_test(optimized);
  const detail::MgMatrixWorkCounters reference_work =
      detail::mg_matrix_work_counters_for_test(reference);
  const HaloRuntimeCounters optimized_halo = difference(
      halo_counters(optimized_fixture), optimized_halo_before);
  const HaloRuntimeCounters reference_halo = difference(
      halo_counters(reference_fixture), reference_halo_before);
  const LinearReductionCounters optimized_reductions = difference(
      optimized_fixture.reductions.counters(), optimized_reductions_before);
  const LinearReductionCounters reference_reductions = difference(
      reference_fixture.reductions.counters(), reference_reductions_before);
  const bool output_same =
      optimized_fixture.correction.storage.size() ==
          reference_fixture.correction.storage.size() &&
      std::memcmp(optimized_fixture.correction.storage.data(),
                  reference_fixture.correction.storage.data(),
                  optimized_fixture.correction.storage.size() * sizeof(double)) ==
          0;

  passed &= expect(
      static_cast<bool>(optimized_status) &&
          packed(optimized_status) == packed(reference_status) && output_same &&
          bitwise_same(optimized.last_cycle_initial_residual(),
                       reference.last_cycle_initial_residual()) &&
          bitwise_same(optimized.last_cycle_final_residual(),
                       reference.last_cycle_final_residual()),
      rank, "fused color updates are bitwise identical to separated updates");
  passed &= expect(
      same(optimized.counters(), reference.counters()) &&
          same(optimized_halo, reference_halo) &&
          same(optimized_reductions, reference_reductions),
      rank, "fused updates preserve plan, halo, and reduction schedules");
  const bool work_contract =
      optimized_work.fused_color_actions > 0U &&
      optimized_work.fused_cell_visits == optimized_work.fused_updates &&
      optimized_work.full_actions > 0U &&
      optimized_work.retained_final_defect_actions > 0U &&
      optimized_work.defect_cell_visits == optimized_work.defect_writes &&
      optimized_work.residual_finish_actions <
          reference_work.residual_finish_actions &&
      reference_work.fused_color_actions == 0U &&
      reference_work.separated_color_actions > 0U &&
      reference_work.separated_color_updates > 0U &&
      reference_work.retained_final_defect_actions == 0U &&
      (zero_pre_sweeps ? reference_work.retained_final_full_actions == 0U
                       : reference_work.retained_final_full_actions > 0U) &&
          reference_work.full_actions > optimized_work.full_actions &&
          reference_work.residual_writes > optimized_work.residual_writes &&
          reference_work.full_cell_visits > optimized_work.fused_cell_visits;
  passed &= expect(work_contract, rank,
                   "fused path retains direct defect and removes redundant finish pass");
  if (zero_pre_sweeps) {
    passed &= expect(
        optimized_work.retained_final_defect_actions > 0U &&
            optimized_work.defect_cell_visits > 0U &&
            optimized_work.defect_cell_visits == optimized_work.defect_writes,
        rank, "zero pre-sweeps still forms rhs minus A*x directly");
  }
  return all_true(passed);
}

bool test_odd_periodic_point_smoother_fallback(int rank) {
  constexpr HaloTopology periodic_x{true, false, false};
  Fixture fixture;
  bool passed = expect(initialize(fixture, periodic_x, 4.0, {25, 12, 8}), rank,
                       "odd-periodic MG fixture initializes");
  if (!all_true(passed)) return false;
  NativeCartesianMgSpec spec = mg_spec(fixture, 920U);
  spec.boundaries = boundaries(periodic_x);
  spec.policy.maximum_levels = 2U;
  NativeCartesianMgPlan plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          spec, services(fixture), coefficient_views(fixture), plan)) &&
          plan.line_axis_mask() == 0U,
      rank, "odd-periodic point plan compiles without line relaxation");
  if (!all_true(passed)) return false;

  detail::set_mg_reference_point_actions_for_test(plan, false);
  const Status status = plan.apply(as_const(fixture.residual.view),
                                   fixture.correction.view, 0U);
  const detail::MgMatrixWorkCounters work =
      detail::mg_matrix_work_counters_for_test(plan);
  passed &= expect(static_cast<bool>(status) && finite(fixture.correction.view),
                   rank, "odd-periodic point apply succeeds");
  passed &= expect(work.separated_color_actions > 0U &&
                       work.separated_color_updates > 0U &&
                       work.full_actions > 0U &&
                       work.retained_final_full_actions == 0U &&
                       work.retained_final_defect_actions > 0U &&
                       work.defect_cell_visits == work.defect_writes,
                   rank,
                   "odd periodic extent uses separated updates with retained defect");
  return all_true(passed);
}

bool test_explicit_chebyshev_distributed_routes(int rank) {
  constexpr HaloTopology periodic_x{true, false, false};
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const MgOperatorClass certified =
      MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
  Fixture odd_periodic;
  bool passed = expect(
      initialize(odd_periodic, periodic_x, 4.0, {25, 12, 8}, 2U), rank,
      "explicit Chebyshev odd-periodic fixture initializes");
  if (!all_true(passed)) return false;
  NativeCartesianMgSpec odd_spec = mg_spec(
      odd_periodic, 925U, MgPointSmootherKind::chebyshev_jacobi, certified);
  odd_spec.boundaries = boundaries(periodic_x);
  NativeCartesianMgPlan odd_plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          odd_spec, services(odd_periodic), coefficient_views(odd_periodic),
          odd_plan)) &&
          odd_plan.line_axis_mask() == 0U && odd_plan.level_count() == 2U &&
          detail::mg_replicated_coarse_enabled_for_test(odd_plan) ==
              (size > 1),
      rank,
      "explicit Chebyshev compiles on odd-periodic point levels while retaining the distributed/replicated coarse choice");
  if (!all_true(passed)) return false;
  std::size_t allocations = 0U;
  Status odd_status{};
  {
    allocation_observer::Guard guard;
    odd_status = odd_plan.apply(as_const(odd_periodic.residual.view),
                                odd_periodic.correction.view, 0U);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const detail::MgMatrixWorkCounters odd_work =
      detail::mg_matrix_work_counters_for_test(odd_plan);
  const std::size_t local_cells =
      static_cast<std::size_t>(odd_periodic.patch.cells.x) *
      static_cast<std::size_t>(odd_periodic.patch.cells.y) *
      static_cast<std::size_t>(odd_periodic.patch.cells.z);
  passed &= expect(
      static_cast<bool>(odd_status) && finite(odd_periodic.correction.view) &&
          allocations == 0U && odd_work.chebyshev_stages == 4U &&
          odd_work.chebyshev_exchange_actions == 5U &&
          odd_work.chebyshev_defect_actions == 5U &&
          odd_work.chebyshev_retained_final_defect_actions == 1U &&
          odd_work.chebyshev_updates == 4U * local_cells,
      rank,
      "odd-periodic distributed Chebyshev keeps exact pre/post work and zero hot allocation");

  Fixture nullspace;
  passed &= expect(initialize(nullspace, {}, 4.0, {25, 12, 8}, 2U), rank,
                   "explicit Chebyshev null-space fixture initializes");
  if (!all_true(passed)) return false;
  NativeCartesianMgSpec null_spec = mg_spec(
      nullspace, 926U, MgPointSmootherKind::chebyshev_jacobi, certified);
  null_spec.boundaries = {
      MgBoundaryKind::neumann, MgBoundaryKind::neumann,
      MgBoundaryKind::neumann, MgBoundaryKind::neumann,
      MgBoundaryKind::neumann, MgBoundaryKind::neumann};
  null_spec.null_space = MgNullSpace::constant;
  NativeCartesianMgPlan null_plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          null_spec, services(nullspace), coefficient_views(nullspace),
          null_plan)) &&
          null_plan.line_axis_mask() == 0U && null_plan.level_count() == 2U,
      rank, "explicit Chebyshev compiles with a distributed constant null space");
  if (!all_true(passed)) return false;
  for (std::int32_t k = 0; k < nullspace.patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < nullspace.patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < nullspace.patch.cells.x; ++i) {
        nullspace.residual.view.unchecked({i, j, k}, 0U) =
            static_cast<double>(i + 2 * j + 3 * k);
      }
    }
  }
  allocations = 0U;
  Status null_status{};
  {
    allocation_observer::Guard guard;
    null_status = null_plan.apply(as_const(nullspace.residual.view),
                                  nullspace.correction.view, 0U);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const detail::MgMatrixWorkCounters null_work =
      detail::mg_matrix_work_counters_for_test(null_plan);
  passed &= expect(
      static_cast<bool>(null_status) && finite(nullspace.correction.view) &&
          allocations == 0U && null_work.chebyshev_stages == 4U &&
          null_work.chebyshev_exchange_actions == 5U &&
          null_work.chebyshev_defect_actions == 5U &&
          null_work.chebyshev_retained_final_defect_actions == 1U,
      rank,
      "distributed constant-null-space Chebyshev preserves exact lifecycle and zero hot allocation");
  return all_true(passed);
}

bool test_replicated_coarse_reference_oracle(int rank, int size) {
  Fixture replicated_fixture;
  Fixture reference_fixture;
  bool passed = expect(
      initialize(replicated_fixture, {}, 4.0, {24, 12, 8}) &&
          initialize(reference_fixture, {}, 4.0, {24, 12, 8}), rank,
      "replicated coarse oracle fixtures initialize");
  if (!all_true(passed)) return false;

  NativeCartesianMgPlan replicated;
  NativeCartesianMgPlan reference;
  MgPlanCounters replicated_external{};
  MgPlanCounters reference_external{};
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(replicated_fixture, 930U), services(replicated_fixture),
          coefficient_views(replicated_fixture), replicated,
          &replicated_external)) &&
          static_cast<bool>(NativeCartesianMgPlan::compile(
              mg_spec(reference_fixture, 930U), services(reference_fixture),
              coefficient_views(reference_fixture), reference,
              &reference_external)),
      rank, "replicated coarse oracle plans compile");
  if (!all_true(passed)) return false;

  detail::set_mg_replicated_coarse_mode_for_test(replicated, 2U);
  detail::set_mg_replicated_coarse_mode_for_test(reference, 1U);
  const bool replicated_enabled =
      detail::mg_replicated_coarse_enabled_for_test(replicated);
  const bool reference_enabled =
      detail::mg_replicated_coarse_enabled_for_test(reference);
  passed &= expect(replicated_enabled == (size > 1) && !reference_enabled,
                   rank, "replicated mode obeys communicator eligibility");
  if (!all_true(passed)) return false;

  const MgPlanCounters replicated_before = replicated.counters();
  const MgPlanCounters reference_before = reference.counters();
  const Status replicated_status = replicated.apply(
      as_const(replicated_fixture.residual.view),
      replicated_fixture.correction.view, 0U);
  const Status reference_status = reference.apply(
      as_const(reference_fixture.residual.view), reference_fixture.correction.view,
      0U);
  const bool output_same =
      replicated_fixture.correction.storage.size() ==
          reference_fixture.correction.storage.size() &&
      std::memcmp(replicated_fixture.correction.storage.data(),
                  reference_fixture.correction.storage.data(),
                  replicated_fixture.correction.storage.size() * sizeof(double)) ==
          0;
  passed &= expect(static_cast<bool>(replicated_status) &&
                       static_cast<bool>(reference_status) &&
                       packed(replicated_status) == packed(reference_status) &&
                       output_same &&
                       bitwise_same(replicated.last_cycle_initial_residual(),
                                    reference.last_cycle_initial_residual()) &&
                       bitwise_same(replicated.last_cycle_final_residual(),
                                    reference.last_cycle_final_residual()),
                   rank, "replicated coarse correction is bitwise reference-equivalent");
  const MgPlanCounters replicated_after = replicated.counters();
  const MgPlanCounters reference_after = reference.counters();
  if (replicated_enabled) {
    MgLevelView coarse{};
    passed &= expect(static_cast<bool>(replicated.level(
                         replicated.level_count() - 1U, coarse)),
                     rank, "replicated coarse level is queryable");
    const std::size_t global_cells =
        static_cast<std::size_t>(coarse.global_shape.x) *
        static_cast<std::size_t>(coarse.global_shape.y) *
        static_cast<std::size_t>(coarse.global_shape.z);
    const std::uint64_t expected_compile_bytes =
        static_cast<std::uint64_t>(6U * static_cast<std::size_t>(size) *
                                   sizeof(std::int32_t)) +
        static_cast<std::uint64_t>(7U * global_cells * sizeof(double)) +
        0U;
    const std::uint64_t expected_bytes =
        expected_compile_bytes +
        static_cast<std::uint64_t>(global_cells * sizeof(double));
    passed &= expect(replicated_external.blocking_collectives == 2U &&
                         replicated_external.collective_logical_bytes ==
                             expected_compile_bytes &&
                         replicated_before.blocking_collectives == 2U &&
                         replicated_after.blocking_collectives == 3U &&
                         replicated_after.collective_logical_bytes ==
                             expected_bytes,
                     rank, "replicated descriptor/operator/RHS accounting is exact");
  } else {
    passed &= expect(replicated_external.blocking_collectives == 0U &&
                         replicated_before.blocking_collectives == 0U &&
                         replicated_after.blocking_collectives == 0U,
                     rank, "single-rank coarse route has no replication collectives");
  }
  passed &= expect(reference_external.blocking_collectives ==
                       reference_before.blocking_collectives &&
                       reference_after.blocking_collectives ==
                           reference_before.blocking_collectives,
                   rank, "reference coarse route has no direct MG collectives");
  return all_true(passed);
}

struct ReplicatedCoarseCase {
  Int3 cells{24, 12, 8};
  HaloTopology topology{};
  MgBoundarySet boundaries{};
  MgNullSpace null_space{MgNullSpace::none};
  MgCycleKind cycle{MgCycleKind::v_cycle};
  bool variable{};
  bool activity{};
};

void make_activity(Fixture& fixture, std::vector<std::uint8_t>& cells,
                   std::vector<std::uint8_t>& x_faces,
                   std::vector<std::uint8_t>& y_faces,
                   std::vector<std::uint8_t>& z_faces) {
  const Int3 local = fixture.patch.cells;
  cells.assign(static_cast<std::size_t>(local.x) * local.y * local.z, 1U);
  x_faces.assign(static_cast<std::size_t>(local.x + 1) * local.y * local.z,
                 1U);
  y_faces.assign(static_cast<std::size_t>(local.x) * (local.y + 1) * local.z,
                 1U);
  z_faces.assign(static_cast<std::size_t>(local.x) * local.y * (local.z + 1),
                 1U);
  cells[0U] = 0U;
  x_faces[0U] = 0U;
  x_faces[1U] = 0U;
  y_faces[0U] = 0U;
  y_faces[static_cast<std::size_t>(local.x)] = 0U;
  z_faces[0U] = 0U;
  z_faces[static_cast<std::size_t>(local.x) * local.y] = 0U;
}

bool compare_replicated_coarse_case(int rank, int size,
                                    const ReplicatedCoarseCase& options,
                                    std::uint64_t salt) {
  Fixture replicated_fixture;
  Fixture reference_fixture;
  bool passed = expect(
      initialize(replicated_fixture, options.topology, 4.0, options.cells,
                 32U, options.cycle) &&
          initialize(reference_fixture, options.topology, 4.0, options.cells,
                     32U, options.cycle),
      rank, "replicated coarse case fixtures initialize");
  if (!all_true(passed)) return false;

  for (std::int32_t k = 0; k < replicated_fixture.patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < replicated_fixture.patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < replicated_fixture.patch.cells.x; ++i) {
        const double value =
            1.0 + 0.01 * static_cast<double>(
                             replicated_fixture.patch.begin.x + i +
                             3 * (replicated_fixture.patch.begin.y + j) +
                             5 * (replicated_fixture.patch.begin.z + k));
        replicated_fixture.residual.view.unchecked({i, j, k}, 0U) = value;
        reference_fixture.residual.view.unchecked({i, j, k}, 0U) = value;
      }
    }
  }
  if (options.variable) {
    fill_variable_point_fixture(replicated_fixture);
    fill_variable_point_fixture(reference_fixture);
  }

  std::vector<std::uint8_t> replicated_cells;
  std::vector<std::uint8_t> replicated_x_faces;
  std::vector<std::uint8_t> replicated_y_faces;
  std::vector<std::uint8_t> replicated_z_faces;
  std::vector<std::uint8_t> reference_cells;
  std::vector<std::uint8_t> reference_x_faces;
  std::vector<std::uint8_t> reference_y_faces;
  std::vector<std::uint8_t> reference_z_faces;
  if (options.activity) {
    make_activity(replicated_fixture, replicated_cells, replicated_x_faces,
                  replicated_y_faces, replicated_z_faces);
    make_activity(reference_fixture, reference_cells, reference_x_faces,
                  reference_y_faces, reference_z_faces);
  }
  NativeCartesianMgSpec replicated_spec = mg_spec(replicated_fixture, salt);
  NativeCartesianMgSpec reference_spec = mg_spec(reference_fixture, salt);
  replicated_spec.boundaries = options.boundaries;
  reference_spec.boundaries = options.boundaries;
  replicated_spec.null_space = options.null_space;
  reference_spec.null_space = options.null_space;
  if (options.activity) {
    replicated_spec.activity = {
        {replicated_cells.data(), replicated_cells.size()},
        {replicated_x_faces.data(), replicated_x_faces.size()},
        {replicated_y_faces.data(), replicated_y_faces.size()},
        {replicated_z_faces.data(), replicated_z_faces.size()}, 7811U, 8811U};
    reference_spec.activity = {
        {reference_cells.data(), reference_cells.size()},
        {reference_x_faces.data(), reference_x_faces.size()},
        {reference_y_faces.data(), reference_y_faces.size()},
        {reference_z_faces.data(), reference_z_faces.size()}, 7811U, 8811U};
  }
  NativeCartesianMgPlan replicated;
  NativeCartesianMgPlan reference;
  MgPlanCounters replicated_external{};
  MgPlanCounters reference_external{};
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          replicated_spec, services(replicated_fixture),
          coefficient_views(replicated_fixture), replicated,
          &replicated_external)) &&
          static_cast<bool>(NativeCartesianMgPlan::compile(
              reference_spec, services(reference_fixture),
              coefficient_views(reference_fixture), reference,
              &reference_external)),
      rank, "replicated coarse case plans compile");
  if (!all_true(passed)) return false;
  detail::set_mg_replicated_coarse_mode_for_test(replicated, 2U);
  detail::set_mg_replicated_coarse_mode_for_test(reference, 1U);
  const bool replicated_enabled =
      detail::mg_replicated_coarse_enabled_for_test(replicated);
  passed &= expect(replicated_enabled == (size > 1), rank,
                   "replicated case meets automatic eligibility");
  if (!all_true(passed)) return false;

  std::size_t allocations = 0U;
  Status replicated_status{};
  {
    allocation_observer::Guard guard;
    replicated_status = replicated.apply(
        as_const(replicated_fixture.residual.view),
        replicated_fixture.correction.view, 0U);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const Status reference_status = reference.apply(
      as_const(reference_fixture.residual.view), reference_fixture.correction.view,
      0U);
  const detail::MgMatrixWorkCounters replicated_work =
      detail::mg_matrix_work_counters_for_test(replicated);
  const detail::MgMatrixWorkCounters reference_work =
      detail::mg_matrix_work_counters_for_test(reference);
  const bool output_same =
      replicated_fixture.correction.storage.size() ==
          reference_fixture.correction.storage.size() &&
      std::memcmp(replicated_fixture.correction.storage.data(),
                  reference_fixture.correction.storage.data(),
                  replicated_fixture.correction.storage.size() * sizeof(double)) ==
          0;
  passed &= expect(allocations == 0U, rank,
                   "replicated coarse apply has zero hot allocation");
  passed &= expect(
      exact_cycle_schedule(replicated_work, replicated.level_count(),
                           options.cycle) &&
          exact_cycle_schedule(reference_work, reference.level_count(),
                               options.cycle) &&
          replicated_work.cycle_level_calls ==
              reference_work.cycle_level_calls &&
          replicated_work.cycle_restrictions ==
              reference_work.cycle_restrictions &&
          replicated_work.cycle_prolongations ==
              reference_work.cycle_prolongations,
      rank, "fixed multigrid cycle follows its exact level-call schedule");
  passed &= expect(static_cast<bool>(replicated_status) &&
                       static_cast<bool>(reference_status) &&
                       packed(replicated_status) == packed(reference_status) &&
                       output_same &&
                       bitwise_same(replicated.last_cycle_initial_residual(),
                                    reference.last_cycle_initial_residual()) &&
                       bitwise_same(replicated.last_cycle_final_residual(),
                                    reference.last_cycle_final_residual()),
                   rank, "replicated case matches reference correction/residuals");
  if (!all_true(passed)) return false;

  const MgPlanCounters replicated_before_refresh = replicated.counters();
  const MgPlanCounters reference_before_refresh = reference.counters();
  for (std::int32_t k = 0; k < replicated_fixture.patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < replicated_fixture.patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < replicated_fixture.patch.cells.x; ++i) {
        replicated_fixture.diagonal.view.unchecked({i, j, k}, 0U) += 0.1;
        reference_fixture.diagonal.view.unchecked({i, j, k}, 0U) += 0.1;
      }
    }
  }
  LinearIdentity next = replicated_spec.identity;
  ++next.numeric;
  ++next.hierarchy;
  ++next.fingerprint;
  const MgCoefficientIdentity next_coefficients{salt + 101U, salt + 201U,
                                                0.1};
  const Status replicated_refresh = replicated.update_coefficients(
      next, next_coefficients, coefficient_views(replicated_fixture),
      &replicated_external);
  const Status reference_refresh = reference.update_coefficients(
      next, next_coefficients, coefficient_views(reference_fixture),
      &reference_external);
  passed &= expect(static_cast<bool>(replicated_refresh) &&
                       static_cast<bool>(reference_refresh) &&
                       replicated.counters().numeric_refreshes ==
                           replicated_before_refresh.numeric_refreshes + 1U &&
                       reference.counters().numeric_refreshes ==
                           reference_before_refresh.numeric_refreshes + 1U,
                   rank, "replicated operator refresh publishes transactionally");
  if (replicated_enabled) {
    passed &= expect(
        replicated.counters().blocking_collectives ==
                replicated_before_refresh.blocking_collectives + 1U &&
            reference.counters().blocking_collectives ==
                reference_before_refresh.blocking_collectives + 1U,
        rank, "replicated refresh accounts one operator gather");
  }
  if (!all_true(passed)) return false;

  const MgPlanCounters before_reuse = replicated.counters();
  const Status reused = replicated.update_coefficients(
      next, next_coefficients, coefficient_views(replicated_fixture),
      &replicated_external);
  passed &= expect(static_cast<bool>(reused) &&
                       same(replicated.counters(), before_reuse),
                   rank, "unchanged coefficient refresh reuses active operator");
  if (!all_true(passed)) return false;

  const PlanFingerprint numeric_before_failure =
      replicated.numeric_fingerprint();
  const PlanFingerprint hierarchy_before_failure =
      replicated.hierarchy_fingerprint();
  const RevisionToken generation_before_failure = replicated.generation();
  const MgPlanCounters counters_before_failure = replicated.counters();
  if (rank == size - 1) {
    replicated_fixture.diagonal.view.unchecked({0, 0, 0}, 0U) =
        std::numeric_limits<double>::quiet_NaN();
    reference_fixture.diagonal.view.unchecked({0, 0, 0}, 0U) =
        std::numeric_limits<double>::quiet_NaN();
  }
  LinearIdentity failed_identity = next;
  ++failed_identity.numeric;
  ++failed_identity.hierarchy;
  ++failed_identity.fingerprint;
  const Status replicated_failed = replicated.update_coefficients(
      failed_identity, MgCoefficientIdentity{salt + 301U, salt + 401U, 0.2},
      coefficient_views(replicated_fixture), &replicated_external);
  const Status reference_failed = reference.update_coefficients(
      failed_identity, MgCoefficientIdentity{salt + 301U, salt + 401U, 0.2},
      coefficient_views(reference_fixture), &reference_external);
  passed &= expect(replicated_failed.code == StatusCode::numerical_failure &&
                       reference_failed.code == StatusCode::numerical_failure &&
                       identical(packed(replicated_failed)) &&
                       replicated.numeric_fingerprint() == numeric_before_failure &&
                       replicated.hierarchy_fingerprint() == hierarchy_before_failure &&
                       replicated.generation() == generation_before_failure &&
                       same(replicated.counters(), counters_before_failure),
                   rank, "failed refresh preserves active replicated operator");
  if (rank == size - 1) {
    replicated_fixture.diagonal.view.unchecked({0, 0, 0}, 0U) = 6.1;
    reference_fixture.diagonal.view.unchecked({0, 0, 0}, 0U) = 6.1;
  }
  if (!all_true(passed)) return false;
  fill(replicated_fixture.correction, -17.0);
  fill(reference_fixture.correction, -17.0);
  const Status replicated_after_failure = replicated.apply(
      as_const(replicated_fixture.residual.view),
      replicated_fixture.correction.view, 0U);
  const Status reference_after_failure = reference.apply(
      as_const(reference_fixture.residual.view), reference_fixture.correction.view,
      0U);
  const bool correction_after_failure_same =
      replicated_fixture.correction.storage.size() ==
          reference_fixture.correction.storage.size() &&
      std::memcmp(replicated_fixture.correction.storage.data(),
                  reference_fixture.correction.storage.data(),
                  replicated_fixture.correction.storage.size() * sizeof(double)) ==
          0;
  passed &= expect(static_cast<bool>(replicated_after_failure) &&
                       static_cast<bool>(reference_after_failure) &&
                       correction_after_failure_same,
                   rank, "post-failure apply uses unchanged active operator");
  return all_true(passed);
}

bool test_replicated_coarse_over_threshold_fallback(int rank, int size) {
  Fixture fixture;
  bool passed = expect(initialize(fixture, {}, 4.0, {64, 64, 64}, 2U), rank,
                       "over-threshold coarse fixture initializes");
  if (!all_true(passed)) return false;
  NativeCartesianMgPlan plan;
  MgPlanCounters counters{};
  passed &= expect(static_cast<bool>(NativeCartesianMgPlan::compile(
                       mg_spec(fixture, 960U), services(fixture),
                       coefficient_views(fixture), plan, &counters)),
                   rank, "over-threshold coarse plan compiles");
  if (!all_true(passed)) return false;
  detail::set_mg_replicated_coarse_mode_for_test(plan, 2U);
  passed &= expect(!detail::mg_replicated_coarse_enabled_for_test(plan), rank,
                   "over-threshold coarsest level falls back to reference route");
  passed &= expect(counters.blocking_collectives == 0U &&
                       counters.collective_logical_bytes == 0U,
                   rank, "fallback has no replicated collective accounting");
  return all_true(passed);
}

bool test_collective_compile_and_ownership(int rank, int size) {
  Fixture fixture;
  bool passed = expect(initialize(fixture), rank, "MG fixture initializes");
  if (!all_true(passed)) {
    return false;
  }
  NativeCartesianMgPlan plan;
  MgPlanCounters counters{};
  const Status compiled = NativeCartesianMgPlan::compile(
      mg_spec(fixture), services(fixture), coefficient_views(fixture), plan,
      &counters);
  passed &= expect(static_cast<bool>(compiled) && plan.level_count() >= 2U &&
                       plan.symbolic_fingerprint() != 0U &&
                       plan.numeric_fingerprint() != 0U &&
                       plan.hierarchy_fingerprint() != 0U &&
                       plan.generation() != 0U &&
                       plan.hierarchy_storage_address() != 0U &&
                       plan.workspace_storage_address() ==
                           fixture.workspace.storage_address() &&
                       counters.symbolic_builds == 1U &&
                       counters.hierarchy_rebuilds == 1U,
                   rank, "collective MG plan publishes one persistent hierarchy");
  MgLevelView finest{};
  passed &= expect(static_cast<bool>(plan.level(0U, finest)) &&
                       same(finest.global_shape,
                            fixture.geometry.global_cells()) &&
                       same(finest.local_shape, fixture.patch.cells),
                   rank, "finest MG level preserves deterministic patch ownership");
  std::array<std::int64_t, 6U> local{
      fixture.patch.begin.x, fixture.patch.begin.y, fixture.patch.begin.z,
      fixture.patch.cells.x, fixture.patch.cells.y, fixture.patch.cells.z};
  std::vector<std::int64_t> ownership;
  if (rank == 0) {
    ownership.resize(static_cast<std::size_t>(size) * local.size());
  }
  MPI_Gather(local.data(), static_cast<int>(local.size()), MPI_INT64_T,
             ownership.data(), static_cast<int>(local.size()), MPI_INT64_T, 0,
             MPI_COMM_WORLD);
  if (rank == 0) {
    const Int3 global = fixture.geometry.global_cells();
    std::vector<std::uint8_t> owners(
        static_cast<std::size_t>(global.x) *
            static_cast<std::size_t>(global.y) *
            static_cast<std::size_t>(global.z),
        0U);
    for (int owner = 0; owner < size; ++owner) {
      const std::int64_t* const record =
          ownership.data() + static_cast<std::size_t>(owner) * local.size();
      for (std::int64_t z = record[2]; z < record[2] + record[5]; ++z) {
        for (std::int64_t y = record[1]; y < record[1] + record[4]; ++y) {
          for (std::int64_t x = record[0]; x < record[0] + record[3]; ++x) {
            const std::size_t index =
                static_cast<std::size_t>(x) +
                static_cast<std::size_t>(global.x) *
                    (static_cast<std::size_t>(y) +
                     static_cast<std::size_t>(global.y) *
                         static_cast<std::size_t>(z));
            if (index < owners.size() && owners[index] !=
                                             std::numeric_limits<std::uint8_t>::max()) {
              ++owners[index];
            }
          }
        }
      }
    }
    passed &= expect(std::all_of(owners.begin(), owners.end(),
                                 [](std::uint8_t owner_count) {
                                   return owner_count == 1U;
                                 }),
                     rank, "MG ownership covers every global cell exactly once");
  }
  passed &= expect(identical(plan.symbolic_fingerprint()) &&
                       identical(plan.hierarchy_fingerprint()),
                   rank, "collective MG structural identities agree on all ranks");
  return all_true(passed);
}

MgBoundarySet boundaries(HaloTopology topology) noexcept {
  const auto kind = [](bool periodic) noexcept {
    return periodic ? MgBoundaryKind::periodic : MgBoundaryKind::dirichlet;
  };
  return {
      kind(topology.periodic_x), kind(topology.periodic_x),
      kind(topology.periodic_y), kind(topology.periodic_y),
      kind(topology.periodic_z), kind(topology.periodic_z),
  };
}

void fill_variable_point_fixture(Fixture& fixture) noexcept {
  const Int3 cells = fixture.patch.cells;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double coordinate = static_cast<double>(i + 2 * j + 3 * k);
        fixture.diagonal.view.unchecked({i, j, k}, 0U) =
            6.0 + 0.01 * coordinate;
        fixture.residual.view.unchecked({i, j, k}, 0U) =
            1.0 + 0.02 * coordinate;
      }
    }
  }
  const auto fill_face = [](OwnedFaceField& face) noexcept {
    for (std::int32_t k = 0; k < face.view.extents.z; ++k) {
      for (std::int32_t j = 0; j < face.view.extents.y; ++j) {
        for (std::int32_t i = 0; i < face.view.extents.x; ++i) {
          face.view.unchecked({i, j, k}) =
              0.4 + 0.01 * static_cast<double>(i + 2 * j + 3 * k);
        }
      }
    }
  };
  fill_face(fixture.x);
  fill_face(fixture.y);
  fill_face(fixture.z);
}

bool test_point_row_boundary_and_small_x_oracle(int rank) {
  const MgBoundarySet dirichlet{
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet,
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet,
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet};
  const MgBoundarySet neumann{
      MgBoundaryKind::neumann, MgBoundaryKind::neumann,
      MgBoundaryKind::neumann, MgBoundaryKind::neumann,
      MgBoundaryKind::neumann, MgBoundaryKind::neumann};
  const MgBoundarySet periodic{
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::periodic, MgBoundaryKind::periodic};
  const MgBoundarySet mixed{
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::neumann, MgBoundaryKind::dirichlet,
      MgBoundaryKind::dirichlet, MgBoundaryKind::neumann};
  struct Case {
    Int3 cells;
    HaloTopology topology;
    MgBoundarySet boundary;
    bool variable;
    std::uint64_t salt;
  };
  const std::array<Case, 4U> cases{{
      {{1, 5, 3}, {}, dirichlet, false, 940U},
      {{2, 5, 3}, {}, neumann, true, 941U},
      {{3, 6, 6}, {true, true, true}, periodic, true, 942U},
      {{3, 6, 6}, {true, false, false}, mixed, true, 943U},
  }};

  bool passed = true;
  for (const Case& selected : cases) {
    Fixture optimized_fixture;
    Fixture reference_fixture;
    passed &= expect(
        initialize(optimized_fixture, selected.topology, 100.0,
                   selected.cells, 2U, MgCycleKind::f_cycle) &&
            initialize(reference_fixture, selected.topology, 100.0,
                       selected.cells, 2U, MgCycleKind::f_cycle),
        rank, "point-row boundary/small-x fixtures initialize");
    if (!all_true(passed)) return false;
    if (selected.variable) {
      fill_variable_point_fixture(optimized_fixture);
      fill_variable_point_fixture(reference_fixture);
    }
    NativeCartesianMgSpec optimized_spec =
        mg_spec(optimized_fixture, selected.salt);
    NativeCartesianMgSpec reference_spec =
        mg_spec(reference_fixture, selected.salt);
    optimized_spec.boundaries = selected.boundary;
    reference_spec.boundaries = selected.boundary;
    NativeCartesianMgPlan optimized;
    NativeCartesianMgPlan reference;
    passed &= expect(
        static_cast<bool>(NativeCartesianMgPlan::compile(
            optimized_spec, services(optimized_fixture),
            coefficient_views(optimized_fixture), optimized)) &&
            static_cast<bool>(NativeCartesianMgPlan::compile(
                reference_spec, services(reference_fixture),
                coefficient_views(reference_fixture), reference)),
        rank, "point-row boundary/small-x plans compile");
    if (!all_true(passed)) return false;
    detail::set_mg_reference_point_actions_for_test(optimized, false);
    detail::set_mg_reference_point_actions_for_test(reference, true);
    const Status optimized_status = optimized.apply(
        as_const(optimized_fixture.residual.view),
        optimized_fixture.correction.view, 0U);
    const Status reference_status = reference.apply(
        as_const(reference_fixture.residual.view),
        reference_fixture.correction.view, 0U);
    const bool output_same =
        optimized_fixture.correction.storage.size() ==
            reference_fixture.correction.storage.size() &&
        std::memcmp(optimized_fixture.correction.storage.data(),
                    reference_fixture.correction.storage.data(),
                    optimized_fixture.correction.storage.size() * sizeof(double)) ==
            0;
    passed &= expect(
        static_cast<bool>(optimized_status) &&
            packed(optimized_status) == packed(reference_status) && output_same &&
            bitwise_same(optimized.last_cycle_initial_residual(),
                         reference.last_cycle_initial_residual()) &&
            bitwise_same(optimized.last_cycle_final_residual(),
                         reference.last_cycle_final_residual()),
        rank, "point-row boundary/small-x optimized path is bitwise oracle-equivalent");
  }
  return all_true(passed);
}

bool test_chebyshev_point_row_reference_oracle(int rank) {
  constexpr HaloTopology periodic_x{true, false, false};
  constexpr Int3 cells{25, 12, 8};
  const MgBoundarySet mixed{
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::neumann, MgBoundaryKind::dirichlet,
      MgBoundaryKind::dirichlet, MgBoundaryKind::neumann};
  bool passed = true;
  for (std::uint8_t degree = 1U; degree <= 4U; ++degree) {
    Fixture optimized_fixture;
    Fixture reference_fixture;
    passed &= expect(
        initialize(optimized_fixture, periodic_x, 4.0, cells, 2U,
                   MgCycleKind::f_cycle) &&
            initialize(reference_fixture, periodic_x, 4.0, cells, 2U,
                       MgCycleKind::f_cycle), rank,
        "Chebyshev point-row oracle fixtures initialize for degree 1-4");
    if (!all_true(passed)) return false;
    fill_variable_point_fixture(optimized_fixture);
    fill_variable_point_fixture(reference_fixture);

    NativeCartesianMgSpec optimized_spec = mg_spec(
        optimized_fixture, 944U, MgPointSmootherKind::chebyshev_jacobi,
        MgOperatorClass::symmetric_diagonally_dominant_m_matrix);
    NativeCartesianMgSpec reference_spec = mg_spec(
        reference_fixture, 944U, MgPointSmootherKind::chebyshev_jacobi,
        MgOperatorClass::symmetric_diagonally_dominant_m_matrix);
    optimized_spec.boundaries = mixed;
    reference_spec.boundaries = mixed;
    optimized_spec.policy.pre_sweeps = degree;
    optimized_spec.policy.post_sweeps = degree;
    reference_spec.policy.pre_sweeps = degree;
    reference_spec.policy.post_sweeps = degree;

    NativeCartesianMgPlan optimized;
    NativeCartesianMgPlan reference;
    passed &= expect(
        static_cast<bool>(NativeCartesianMgPlan::compile(
            optimized_spec, services(optimized_fixture),
            coefficient_views(optimized_fixture), optimized)) &&
            static_cast<bool>(NativeCartesianMgPlan::compile(
                reference_spec, services(reference_fixture),
                coefficient_views(reference_fixture), reference)) &&
            optimized.line_axis_mask() == 0U && reference.line_axis_mask() == 0U,
        rank, "Chebyshev point-row oracle plans compile for degree 1-4");
    if (!all_true(passed)) return false;
    detail::set_mg_reference_point_actions_for_test(optimized, false);
    detail::set_mg_reference_point_actions_for_test(reference, false);
    detail::set_mg_reference_chebyshev_lifecycle_for_test(reference, true);

    const HaloRuntimeCounters optimized_halo_before =
        halo_counters(optimized_fixture);
    const HaloRuntimeCounters reference_halo_before =
        halo_counters(reference_fixture);
    const LinearReductionCounters optimized_reductions_before =
        optimized_fixture.reductions.counters();
    const LinearReductionCounters reference_reductions_before =
        reference_fixture.reductions.counters();
    const MgPlanCounters optimized_plan_before = optimized.counters();
    const MgPlanCounters reference_plan_before = reference.counters();

    const Status optimized_status = optimized.apply(
        as_const(optimized_fixture.residual.view),
        optimized_fixture.correction.view, 0U);
    const Status reference_status = reference.apply(
        as_const(reference_fixture.residual.view),
        reference_fixture.correction.view, 0U);
    const detail::MgMatrixWorkCounters optimized_work =
        detail::mg_matrix_work_counters_for_test(optimized);
    const detail::MgMatrixWorkCounters reference_work =
        detail::mg_matrix_work_counters_for_test(reference);
    const HaloRuntimeCounters optimized_halo = difference(
        halo_counters(optimized_fixture), optimized_halo_before);
    const HaloRuntimeCounters reference_halo = difference(
        halo_counters(reference_fixture), reference_halo_before);
    const LinearReductionCounters optimized_reductions = difference(
        optimized_fixture.reductions.counters(), optimized_reductions_before);
    const LinearReductionCounters reference_reductions = difference(
        reference_fixture.reductions.counters(), reference_reductions_before);
    const MgPlanCounters optimized_plan_after = optimized.counters();
    const MgPlanCounters reference_plan_after = reference.counters();
    const Int3 local_shape = optimized_fixture.patch.cells;
    const std::uint64_t local_cells =
        static_cast<std::uint64_t>(local_shape.x) *
        static_cast<std::uint64_t>(local_shape.y) *
        static_cast<std::uint64_t>(local_shape.z);
    const std::uint64_t stages = 2U * static_cast<std::uint64_t>(degree);
    const std::uint64_t defects = stages + 1U;
    const bool output_same =
        optimized_fixture.correction.storage.size() ==
            reference_fixture.correction.storage.size() &&
        std::memcmp(optimized_fixture.correction.storage.data(),
                    reference_fixture.correction.storage.data(),
                    optimized_fixture.correction.storage.size() * sizeof(double)) ==
            0;
    passed &= expect(
        static_cast<bool>(optimized_status) &&
            packed(optimized_status) == packed(reference_status) && output_same &&
            bitwise_same(optimized.last_cycle_initial_residual(),
                         reference.last_cycle_initial_residual()) &&
            bitwise_same(optimized.last_cycle_final_residual(),
                         reference.last_cycle_final_residual()),
        rank, "Chebyshev point-row correction/residual/status are bitwise oracle-equivalent");
    passed &= expect(
        same_chebyshev_scientific_work(optimized_work, reference_work) &&
            same(optimized_halo, reference_halo) &&
            same(optimized_reductions, reference_reductions) &&
            same(optimized_plan_after, reference_plan_after) &&
            same(optimized_plan_before, reference_plan_before),
        rank, "Chebyshev point-row work, halo, reduction, and plan schedules are identical");
    passed &= expect(
        optimized_work.chebyshev_stages == stages &&
            optimized_work.chebyshev_stages == reference_work.chebyshev_stages &&
            optimized_work.chebyshev_exchange_actions == stages + 1U &&
            optimized_work.chebyshev_exchange_actions ==
                reference_work.chebyshev_exchange_actions &&
            optimized_work.chebyshev_defect_actions == defects &&
            optimized_work.chebyshev_defect_actions ==
                reference_work.chebyshev_defect_actions &&
            optimized_work.chebyshev_updates == stages * local_cells &&
            optimized_work.chebyshev_updates == reference_work.chebyshev_updates &&
            optimized_work.chebyshev_stencil_evaluations == defects * local_cells &&
            optimized_work.chebyshev_stencil_evaluations ==
                reference_work.chebyshev_stencil_evaluations &&
            optimized_work.defect_writes == local_cells &&
            reference_work.defect_writes == defects * local_cells &&
            optimized_work.chebyshev_elided_intermediate_residual_publications ==
                stages * local_cells &&
            reference_work.chebyshev_elided_intermediate_residual_publications ==
                0U &&
            optimized_work.chebyshev_copyback_cells ==
                ((degree & 1U) != 0U ? 2U * local_cells : 0U) &&
            reference_work.chebyshev_copyback_cells == 0U,
        rank, "Chebyshev point-row oracle executes degrees 1-4 with declared lifecycle differences");
  }
  return all_true(passed);
}

bool test_boundary_topology_contract_is_bidirectional(int rank) {
  constexpr HaloTopology periodic_halo{true, false, true};
  Fixture periodic_fixture;
  bool passed = expect(initialize(periodic_fixture, periodic_halo), rank,
                       "periodic MG topology fixture initializes");
  if (!all_true(passed)) {
    return false;
  }

  NativeCartesianMgSpec matching = mg_spec(periodic_fixture, 700U);
  matching.boundaries = boundaries(periodic_halo);
  NativeCartesianMgPlan matching_plan;
  const Status matching_status = NativeCartesianMgPlan::compile(
      matching, services(periodic_fixture), coefficient_views(periodic_fixture),
      matching_plan);
  NativeCartesianMgPlan periodic_to_nonperiodic_plan;
  const Status periodic_to_nonperiodic = NativeCartesianMgPlan::compile(
      mg_spec(periodic_fixture, 710U), services(periodic_fixture),
      coefficient_views(periodic_fixture), periodic_to_nonperiodic_plan);

  Fixture nonperiodic_fixture;
  passed &= expect(initialize(nonperiodic_fixture), rank,
                   "nonperiodic MG topology fixture initializes");
  if (!all_true(passed)) {
    return false;
  }
  NativeCartesianMgSpec expects_periodic =
      mg_spec(nonperiodic_fixture, 720U);
  expects_periodic.boundaries = boundaries(periodic_halo);
  NativeCartesianMgPlan nonperiodic_to_periodic_plan;
  const Status nonperiodic_to_periodic = NativeCartesianMgPlan::compile(
      expects_periodic, services(nonperiodic_fixture),
      coefficient_views(nonperiodic_fixture), nonperiodic_to_periodic_plan);

  passed &= expect(
      static_cast<bool>(matching_status) &&
          periodic_to_nonperiodic.code == StatusCode::invalid_plan &&
          nonperiodic_to_periodic.code == StatusCode::invalid_plan &&
          identical(packed(periodic_to_nonperiodic)) &&
          identical(packed(nonperiodic_to_periodic)),
      rank,
      "MG maps boundary pairs into halo topology and rejects both mismatch directions");
  return all_true(passed);
}

bool test_prepared_apply_equivalence(int rank, int size) {
  Fixture direct_fixture;
  Fixture prepared_fixture;
  bool passed = expect(initialize(direct_fixture), rank,
                       "Native-MG direct equivalence fixture initializes");
  passed &= expect(initialize(prepared_fixture), rank,
                   "Native-MG prepared equivalence fixture initializes");
  if (!all_true(passed)) return false;

  NativeCartesianMgPlan direct;
  NativeCartesianMgPlan prepared;
  MgPlanCounters direct_external{};
  MgPlanCounters prepared_external{};
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(direct_fixture), services(direct_fixture),
          coefficient_views(direct_fixture), direct, &direct_external)),
      rank, "Native-MG direct equivalence plan compiles");
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(prepared_fixture), services(prepared_fixture),
          coefficient_views(prepared_fixture), prepared,
          &prepared_external)),
      rank, "Native-MG prepared equivalence plan compiles");
  if (!all_true(passed)) return false;
  passed &= expect(
      prepared.certificate().apply_lifecycle ==
              LinearPreconditionerApplyLifecycle::prepared_batch &&
          prepared.certificate().status_scope ==
              LinearPreconditionerStatusScope::collective,
      rank, "Native-MG certificate claims only the collective prepared lifecycle");

  const LinearReductionCounters direct_before =
      direct_fixture.reductions.counters();
  const LinearReductionCounters prepared_before =
      prepared_fixture.reductions.counters();
  const HaloRuntimeCounters direct_halo_before =
      halo_counters(direct_fixture);
  const HaloRuntimeCounters prepared_halo_before =
      halo_counters(prepared_fixture);
  const MgPlanCounters direct_plan_before = direct.counters();
  const MgPlanCounters prepared_plan_before = prepared.counters();
  const Status direct_status = direct.apply(
      as_const(direct_fixture.residual.view), direct_fixture.correction.view,
      0U);

  const LinearPreconditionerBatchDescriptor descriptor{
      &prepared_fixture.solver_workspace, prepared_fixture.patch.cells, 1U,
      6U, 4U, 2U};
  LinearPreconditionerBatchTicket ticket;
  std::size_t allocations = 0U;
  Status cold_status{};
  Status prepared_status{};
  {
    allocation_observer::Guard guard;
    cold_status = prepared.prepare_batch(descriptor, ticket);
    if (cold_status) {
      prepared_status = prepared.apply_prepared(
          as_const(prepared_fixture.residual.view),
          prepared_fixture.correction.view, 0U, ticket);
    }
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const LinearReductionCounters direct_after =
      direct_fixture.reductions.counters();
  const LinearReductionCounters prepared_after =
      prepared_fixture.reductions.counters();
  const HaloRuntimeCounters direct_halo_after =
      halo_counters(direct_fixture);
  const HaloRuntimeCounters prepared_halo_after =
      halo_counters(prepared_fixture);
  const MgPlanCounters direct_plan_after = direct.counters();
  const MgPlanCounters prepared_plan_after = prepared.counters();
  const HaloRuntimeCounters direct_halo_delta =
      difference(direct_halo_after, direct_halo_before);
  const HaloRuntimeCounters prepared_halo_delta =
      difference(prepared_halo_after, prepared_halo_before);
  const std::uint64_t direct_checksum = checksum(direct_fixture.correction.view);
  const std::uint64_t prepared_checksum =
      checksum(prepared_fixture.correction.view);
  const bool direct_work_same =
      direct_plan_after.symbolic_builds == prepared_plan_after.symbolic_builds &&
      direct_plan_after.numeric_refreshes ==
          prepared_plan_after.numeric_refreshes &&
      direct_plan_after.hierarchy_rebuilds ==
          prepared_plan_after.hierarchy_rebuilds &&
      direct_plan_after.applications == prepared_plan_after.applications;
  const bool direct_reduction_delta =
      direct_after.calls - direct_before.calls ==
              prepared_after.calls - prepared_before.calls &&
      direct_after.blocking_operations - direct_before.blocking_operations ==
          prepared_after.blocking_operations -
                  prepared_before.blocking_operations +
              5U;
  const bool exact_equivalence =
      static_cast<bool>(direct_status) && static_cast<bool>(cold_status) &&
      static_cast<bool>(prepared_status) && allocations == 0U &&
      direct_plan_after.applications == direct_plan_before.applications + 1U &&
      prepared_plan_after.applications ==
          prepared_plan_before.applications + 1U &&
      same(direct_external, prepared_external) && direct_work_same &&
      direct_reduction_delta && direct_checksum == prepared_checksum &&
      same(direct_halo_delta, prepared_halo_delta) &&
      bitwise_same(direct.last_cycle_initial_residual(),
                   prepared.last_cycle_initial_residual()) &&
      bitwise_same(direct.last_cycle_final_residual(),
                   prepared.last_cycle_final_residual()) &&
      finite(direct_fixture.correction.view) &&
      finite(prepared_fixture.correction.view);
  passed &= expect(
      exact_equivalence, rank,
      "Native-MG direct/prepared paths have exact work and numerical identity with five fewer blocking collectives");
  passed &= expect(
      direct_halo_delta.control_consensus_calls ==
              4U * direct_halo_delta.begin_calls &&
          prepared_halo_delta.control_consensus_calls == 0U &&
          prepared_halo_delta.begin_calls > 0U,
      rank,
      "prepared Native-MG performs no per-Halo control consensus");

  fill(prepared_fixture.correction, -71.0);
  const std::uint64_t failure_checksum =
      checksum(prepared_fixture.correction.view);
  const HaloRuntimeCounters failure_halo_before =
      halo_counters(prepared_fixture);
  detail::set_halo_failure_for_test(detail::HaloFailurePoint::unpack,
                                    size - 1);
  Status deferred_failure{};
  std::size_t failure_allocations = 0U;
  {
    allocation_observer::Guard guard;
    deferred_failure = prepared.apply_prepared(
        as_const(prepared_fixture.residual.view),
        prepared_fixture.correction.view, 1U, ticket);
    failure_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  detail::clear_halo_failure_for_test();
  const HaloRuntimeCounters failure_halo = difference(
      halo_counters(prepared_fixture), failure_halo_before);
  bool ghost_unpublished =
      prepared_fixture.halo.ghost_revision(kWorkspaceVectors) == 0U;
  for (const HaloEngine& halo : prepared_fixture.coarse_halos) {
    ghost_unpublished &= halo.ghost_revision(kWorkspaceVectors) == 0U;
  }
  passed &= expect(
      deferred_failure.code == StatusCode::invalid_plan &&
          deferred_failure.detail == detail::halo_detail_unpack_failure &&
          identical(packed(deferred_failure)) &&
          prepared.lowest_failing_rank() == size - 1 &&
          checksum(prepared_fixture.correction.view) == failure_checksum &&
          failure_allocations == 0U && ghost_unpublished &&
          failure_halo.begin_calls == prepared_halo_delta.begin_calls &&
          failure_halo.finish_calls == prepared_halo_delta.finish_calls &&
          failure_halo.control_consensus_calls == 0U,
      rank,
      "prepared Native-MG defers a rank-local Halo failure through the fixed schedule and publishes no ghost certificate or correction");

  LinearPreconditionerBatchTicket invalid_ticket;
  const MgPlanCounters before_invalid = prepared.counters();
  const Status invalid = prepared.apply_prepared(
      as_const(prepared_fixture.residual.view), prepared_fixture.correction.view,
      0U, invalid_ticket);
  const Status over_boundary = prepared.apply_prepared(
      as_const(prepared_fixture.residual.view), prepared_fixture.correction.view,
      2U, ticket);
  passed &= expect(
      invalid.code == StatusCode::invalid_plan &&
          over_boundary.code == StatusCode::invalid_plan &&
          same(prepared.counters(), before_invalid),
      rank,
      "Native-MG rejects invalid/exhausted tickets and the application-counter boundary before hot work");

  MgPlanCounters overflow_counters{};
  overflow_counters.applications =
      std::numeric_limits<std::uint64_t>::max();
  detail::set_mg_runtime_counters_for_test(
      prepared, overflow_counters, prepared.generation());
  const MgPlanCounters overflow_before = prepared.counters();
  LinearPreconditionerBatchTicket overflow_ticket;
  const Status overflow_cold =
      prepared.prepare_batch(descriptor, overflow_ticket);
  passed &= expect(overflow_cold.code == StatusCode::invalid_plan &&
                       same(prepared.counters(), overflow_before),
                   rank,
                   "Native-MG cold prepare rejects application-counter overflow before hot work");
  return all_true(passed);
}

bool test_progress_and_hot_reuse(int rank) {
  Fixture fixture;
  bool passed = expect(initialize(fixture, {}, 4.0, {24, 12, 8}, 32U,
                                  MgCycleKind::f_cycle), rank,
                       "MG hot reuse fixture initializes");
  if (!all_true(passed)) {
    return false;
  }
  NativeCartesianMgPlan plan;
  MgPlanCounters external{};
  passed &= expect(static_cast<bool>(NativeCartesianMgPlan::compile(
                       mg_spec(fixture), services(fixture),
                       coefficient_views(fixture), plan, &external)),
                   rank, "MG hot reuse plan compiles");
  if (!all_true(passed)) {
    return false;
  }
  const std::uintptr_t hierarchy = plan.hierarchy_storage_address();
  const std::uintptr_t workspace = plan.workspace_storage_address();
  const HaloRuntimeCounters halo_before = fixture.halo.runtime_counters();
  const LinearReductionCounters reductions_before =
      fixture.reductions.counters();
  Status hot_status{};
  std::size_t allocations = 0U;
  {
    allocation_observer::Guard guard;
    for (std::uint32_t repetition = 0U; repetition < 100U; ++repetition) {
      const Status applied = plan.apply(as_const(fixture.residual.view),
                                        fixture.correction.view, repetition);
      if (!applied && hot_status) {
        hot_status = applied;
      }
    }
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const HaloRuntimeCounters halo_after = fixture.halo.runtime_counters();
  const LinearReductionCounters reductions_after =
      fixture.reductions.counters();
  const MgPlanCounters after = plan.counters();
  passed &= expect(static_cast<bool>(hot_status) && allocations == 0U &&
                       finite(fixture.correction.view) &&
                       after.applications == 100U &&
                       plan.hierarchy_storage_address() == hierarchy &&
                       plan.workspace_storage_address() == workspace,
                   rank, "100 MG applications reuse stable storage with zero allocation");
  passed &= expect(halo_after.begin_calls > halo_before.begin_calls &&
                       halo_after.finish_calls > halo_before.finish_calls &&
                       greater_or_equal(reductions_after, reductions_before) &&
                       reductions_after.calls > reductions_before.calls,
                   rank, "MG applications make both halo and reduction progress");
  return all_true(passed);
}

bool test_collective_failure_is_transactional(int rank, int size) {
  Fixture fixture;
  bool passed = expect(initialize(fixture), rank,
                       "MG transactional fixture initializes");
  if (!all_true(passed)) {
    return false;
  }
  NativeCartesianMgPlan plan;
  MgPlanCounters counters{};
  passed &= expect(static_cast<bool>(NativeCartesianMgPlan::compile(
                       mg_spec(fixture), services(fixture),
                       coefficient_views(fixture), plan, &counters)),
                   rank, "baseline MG plan compiles");
  if (!all_true(passed)) {
    return false;
  }

  const PlanFingerprint symbolic = plan.symbolic_fingerprint();
  const PlanFingerprint numeric = plan.numeric_fingerprint();
  const PlanFingerprint hierarchy = plan.hierarchy_fingerprint();
  const RevisionToken generation = plan.generation();
  const std::uintptr_t hierarchy_address = plan.hierarchy_storage_address();
  const std::uintptr_t workspace_address = plan.workspace_storage_address();
  const MgPlanCounters before = counters;
  NativeCartesianMgSpec mismatched = mg_spec(fixture);
  if (rank == size - 1) {
    mismatched.boundaries.x_min = MgBoundaryKind::periodic;
  }
  const Status rejected = NativeCartesianMgPlan::compile(
      mismatched, services(fixture), coefficient_views(fixture), plan,
      &counters);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       identical(packed(rejected)) &&
                       plan.lowest_failing_rank() == size - 1 &&
                       plan.symbolic_fingerprint() == symbolic &&
                       plan.numeric_fingerprint() == numeric &&
                       plan.hierarchy_fingerprint() == hierarchy &&
                       plan.generation() == generation &&
                       plan.hierarchy_storage_address() == hierarchy_address &&
                       plan.workspace_storage_address() == workspace_address &&
                       same(counters, before),
                   rank, "rank-local compile mismatch rejects collectively and atomically");

  MgCoefficientViews poisoned = coefficient_views(fixture);
  if (rank == size - 1) {
    fixture.diagonal.view.unchecked({0, 0, 0}, 0U) =
        std::numeric_limits<double>::quiet_NaN();
    poisoned.diagonal = as_const(fixture.diagonal.view);
  }
  const Status update = plan.update_coefficients(
      mg_spec(fixture).identity, MgCoefficientIdentity{42U, 52U, 0.1},
      poisoned, &counters);
  passed &= expect(update.code == StatusCode::numerical_failure &&
                       identical(packed(update)) &&
                       plan.lowest_failing_rank() == size - 1 &&
                       plan.numeric_fingerprint() == numeric &&
                       plan.hierarchy_fingerprint() == hierarchy &&
                       plan.generation() == generation &&
                       plan.hierarchy_storage_address() == hierarchy_address &&
                       same(counters, before),
                   rank, "rank-local coefficient failure preserves published hierarchy");

  fill(fixture.correction, -19.0);
  const std::uint64_t correction_before = checksum(fixture.correction.view);
  FieldView bad_output = fixture.correction.view;
  if (rank == size - 1) {
    bad_output.interior.x += 1;
  }
  const Status application =
      plan.apply(as_const(fixture.residual.view), bad_output, 0U);
  passed &= expect(application.code == StatusCode::invalid_plan &&
                       identical(packed(application)) &&
                       plan.lowest_failing_rank() == size - 1 &&
                       checksum(fixture.correction.view) == correction_before &&
                       plan.hierarchy_storage_address() == hierarchy_address &&
                       plan.workspace_storage_address() == workspace_address,
                   rank, "rank-local apply failure leaves caller correction unchanged");
  return all_true(passed);
}

bool test_default_v_cycle_is_explicit_v_cycle(int rank) {
  Fixture default_fixture;
  Fixture explicit_fixture;
  bool passed = expect(initialize(default_fixture) &&
                           initialize(explicit_fixture, {}, 4.0,
                                      {24, 12, 8}, 32U,
                                      MgCycleKind::v_cycle),
                       rank,
                       "default and explicit V-cycle fixtures initialize");
  if (!all_true(passed)) return false;
  NativeCartesianMgPlan default_plan;
  NativeCartesianMgPlan explicit_plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(default_fixture, 1029U), services(default_fixture),
          coefficient_views(default_fixture), default_plan)) &&
          static_cast<bool>(NativeCartesianMgPlan::compile(
              mg_spec(explicit_fixture, 1029U), services(explicit_fixture),
              coefficient_views(explicit_fixture), explicit_plan)),
      rank, "default and explicit V-cycle plans compile");
  if (!all_true(passed)) return false;
  const HaloRuntimeCounters default_halo_before =
      halo_counters(default_fixture);
  const HaloRuntimeCounters explicit_halo_before =
      halo_counters(explicit_fixture);
  const Status default_status = default_plan.apply(
      as_const(default_fixture.residual.view),
      default_fixture.correction.view, 0U);
  const Status explicit_status = explicit_plan.apply(
      as_const(explicit_fixture.residual.view),
      explicit_fixture.correction.view, 0U);
  const detail::MgMatrixWorkCounters default_work =
      detail::mg_matrix_work_counters_for_test(default_plan);
  const detail::MgMatrixWorkCounters explicit_work =
      detail::mg_matrix_work_counters_for_test(explicit_plan);
  const bool correction_same =
      default_fixture.correction.storage.size() ==
          explicit_fixture.correction.storage.size() &&
      std::memcmp(default_fixture.correction.storage.data(),
                  explicit_fixture.correction.storage.data(),
                  default_fixture.correction.storage.size() *
                      sizeof(double)) == 0;
  passed &= expect(
      static_cast<bool>(default_status) &&
          packed(default_status) == packed(explicit_status) &&
          default_fixture.requirements.fingerprint ==
              explicit_fixture.requirements.fingerprint &&
          default_plan.symbolic_fingerprint() ==
              explicit_plan.symbolic_fingerprint() &&
          default_plan.certificate().collective_fingerprint ==
              explicit_plan.certificate().collective_fingerprint &&
          correction_same && same(default_work, explicit_work) &&
          same(difference(halo_counters(default_fixture),
                          default_halo_before),
               difference(halo_counters(explicit_fixture),
                          explicit_halo_before)) &&
          bitwise_same(default_plan.last_cycle_initial_residual(),
                       explicit_plan.last_cycle_initial_residual()) &&
          bitwise_same(default_plan.last_cycle_final_residual(),
                       explicit_plan.last_cycle_final_residual()),
      rank,
      "default V-cycle is bitwise, work, Halo-stage, workspace, public, and prepared identity equivalent to explicit V");
  return all_true(passed);
}

bool test_fixed_f_cycle_phase_failures_and_retry(int rank, int size) {
  Fixture fixture;
  bool passed = expect(
      initialize(fixture, {}, 4.0, {24, 12, 8}, 2U,
                 MgCycleKind::f_cycle),
      rank, "fixed F-cycle phase-failure fixture initializes");
  if (!all_true(passed)) return false;
  NativeCartesianMgPlan plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(fixture, 1030U), services(fixture),
          coefficient_views(fixture), plan)) &&
          plan.level_count() == 2U,
      rank, "fixed F-cycle phase-failure plan compiles with two levels");
  if (!all_true(passed)) return false;

  struct Case {
    detail::MgCycleFailurePhase phase;
    std::size_t level;
    bool deferred;
    std::uint64_t expected_terminal_calls;
  };
  const std::array<Case, 8U> cases{{
      {detail::MgCycleFailurePhase::pre_smooth, 0U, false, 0U},
      {detail::MgCycleFailurePhase::restriction, 0U, false, 0U},
      {detail::MgCycleFailurePhase::first_coarse, 0U, false, 0U},
      {detail::MgCycleFailurePhase::second_coarse, 0U, false, 1U},
      {detail::MgCycleFailurePhase::prolongation, 0U, false, 2U},
      {detail::MgCycleFailurePhase::post_smooth, 0U, false, 2U},
      {detail::MgCycleFailurePhase::terminal, 1U, false, 1U},
      {detail::MgCycleFailurePhase::second_coarse, 0U, true, 2U},
  }};
  for (const Case& selected : cases) {
    fill(fixture.correction, -37.0);
    const std::uint64_t correction_before = checksum(fixture.correction.view);
    const detail::MgMatrixWorkCounters work_before =
        detail::mg_matrix_work_counters_for_test(plan);
    detail::set_mg_cycle_failure_for_test(
        plan, selected.phase, selected.level, size - 1, selected.deferred);
    Status failed{};
    std::size_t allocations = 0U;
    {
      allocation_observer::Guard guard;
      failed = plan.apply(as_const(fixture.residual.view),
                          fixture.correction.view, 0U);
      allocations = allocation_observer::count.load(std::memory_order_relaxed);
    }
    const detail::MgMatrixWorkCounters work_after =
        detail::mg_matrix_work_counters_for_test(plan);
    const StatusCode expected = selected.deferred
                                    ? StatusCode::numerical_failure
                                    : StatusCode::invalid_plan;
    const bool failure_exact =
        failed.code == expected && identical(packed(failed)) &&
            plan.lowest_failing_rank() == size - 1 && allocations == 0U &&
            checksum(fixture.correction.view) == correction_before &&
            work_after.cycle_level_calls[0U] ==
                work_before.cycle_level_calls[0U] + 1U &&
            work_after.cycle_level_calls[1U] ==
                work_before.cycle_level_calls[1U] +
                    selected.expected_terminal_calls;
    if (!failure_exact) {
      std::cerr << "rank " << rank << " F-fault phase="
                << static_cast<unsigned>(selected.phase) << " deferred="
                << selected.deferred << " status="
                << static_cast<unsigned>(failed.code) << '/' << failed.detail
                << " lowest=" << plan.lowest_failing_rank()
                << " alloc=" << allocations << " checksum="
                << (checksum(fixture.correction.view) == correction_before)
                << " calls="
                << work_after.cycle_level_calls[0U] -
                       work_before.cycle_level_calls[0U]
                << '/'
                << work_after.cycle_level_calls[1U] -
                       work_before.cycle_level_calls[1U]
                << " expected=1/" << selected.expected_terminal_calls
                << '\n';
    }
    passed &= expect(
        failure_exact,
        rank,
        "fixed F-cycle phase failure is collective, transactional, and follows exact coarse-call chronology");
    if (!all_true(passed)) return false;

    const detail::MgMatrixWorkCounters retry_before =
        detail::mg_matrix_work_counters_for_test(plan);
    const Status retried = plan.apply(as_const(fixture.residual.view),
                                      fixture.correction.view, 1U);
    const detail::MgMatrixWorkCounters retry_after =
        detail::mg_matrix_work_counters_for_test(plan);
    passed &= expect(
        static_cast<bool>(retried) &&
            retry_after.cycle_level_calls[0U] ==
                retry_before.cycle_level_calls[0U] + 1U &&
            retry_after.cycle_level_calls[1U] ==
                retry_before.cycle_level_calls[1U] + 2U,
        rank,
        "fixed F-cycle consumes one phase fault and retries on the same resources without a clear");
    if (!all_true(passed)) return false;
  }
  return all_true(passed);
}

bool test_fixed_cycle_generic_scientific_work(int rank, int size) {
  bool passed = true;
  for (std::size_t levels = 2U; levels <= 7U; ++levels) {
    const std::int32_t extent = static_cast<std::int32_t>(
        (3U * static_cast<std::size_t>(size)) << (levels - 1U));
    for (const MgCycleKind cycle :
         {MgCycleKind::v_cycle, MgCycleKind::f_cycle}) {
      Fixture fixture;
      passed &= expect(
          initialize(fixture, {}, 1.0e9, {extent, 1, 1},
                     static_cast<std::uint8_t>(levels), cycle),
          rank,
          "generic fixed-cycle scientific-work fixture initializes");
      if (!all_true(passed)) return false;
      NativeCartesianMgSpec spec = mg_spec(
          fixture, 1032U + static_cast<std::uint64_t>(levels),
          MgPointSmootherKind::chebyshev_jacobi,
          MgOperatorClass::symmetric_diagonally_dominant_m_matrix);
      const std::uint8_t pre_sweeps =
          cycle == MgCycleKind::f_cycle ? 1U : 3U;
      const std::uint8_t post_sweeps =
          cycle == MgCycleKind::f_cycle ? 2U : 3U;
      spec.policy.pre_sweeps = pre_sweeps;
      spec.policy.post_sweeps = post_sweeps;
      spec.policy.coarse_sweeps = 24U;
      NativeCartesianMgPlan plan;
      passed &= expect(
          static_cast<bool>(NativeCartesianMgPlan::compile(
              spec, services(fixture), coefficient_views(fixture), plan)) &&
              plan.level_count() == levels && plan.line_axis_mask() == 0U,
          rank,
          "generic fixed-cycle scientific-work plan compiles for two through seven levels");
      if (!all_true(passed)) return false;
      if (cycle == MgCycleKind::f_cycle && levels == 7U) {
        NativeCartesianMgSpec predecessor_spec = spec;
        predecessor_spec.policy.pre_sweeps = 3U;
        predecessor_spec.policy.post_sweeps = 3U;
        NativeCartesianMgPlan predecessor;
        passed &= expect(
            static_cast<bool>(NativeCartesianMgPlan::compile(
                predecessor_spec, services(fixture),
                coefficient_views(fixture), predecessor)) &&
                identical(plan.symbolic_fingerprint()) &&
                identical(predecessor.symbolic_fingerprint()) &&
                plan.symbolic_fingerprint() !=
                    predecessor.symbolic_fingerprint() &&
                identical(plan.certificate().collective_fingerprint) &&
                identical(predecessor.certificate().collective_fingerprint) &&
                plan.certificate().collective_fingerprint !=
                    predecessor.certificate().collective_fingerprint &&
                plan.workspace_storage_address() ==
                    predecessor.workspace_storage_address(),
            rank,
            "seven-level F/1/2 and F/3/3 have collective-distinct symbolic and prepared identities on the same workspace");
        if (!all_true(passed)) return false;
      }
      // Keep this formula oracle on the distributed terminal path so every
      // synthetic terminal stage is backed by a real HaloEngine call.  The
      // separate replicated/distributed equivalence matrix covers the
      // all-gather terminal route and its exact recursive call counts.
      detail::set_mg_replicated_coarse_mode_for_test(plan, 1U);

      const HaloRuntimeCounters halo_before = halo_counters(fixture);
      Status applied{};
      std::size_t allocations = 0U;
      {
        allocation_observer::Guard guard;
        applied = plan.apply(as_const(fixture.residual.view),
                             fixture.correction.view, 0U);
        allocations =
            allocation_observer::count.load(std::memory_order_relaxed);
      }
      const detail::MgMatrixWorkCounters work =
          detail::mg_matrix_work_counters_for_test(plan);
      const HaloRuntimeCounters halo =
          difference(halo_counters(fixture), halo_before);
      const std::uint64_t terminal_calls =
          cycle == MgCycleKind::f_cycle
              ? static_cast<std::uint64_t>(levels)
              : 1U;
      const std::uint64_t nonterminal_calls =
          cycle == MgCycleKind::f_cycle
              ? static_cast<std::uint64_t>(levels * (levels - 1U) / 2U)
              : static_cast<std::uint64_t>(levels - 1U);
      // Every nonterminal visit owns the selected pre+post Chebyshev stages,
      // one additional retained-defect exchange, one restriction exchange, and one
      // prolongation exchange.  Each 24-sweep terminal red/black visit owns
      // an initial exchange plus 2*24-1 post-color exchanges.  The complete
      // apply owns one final finest-level true-residual exchange.
      const std::uint64_t expected_pre_stages =
          static_cast<std::uint64_t>(pre_sweeps) * nonterminal_calls;
      const std::uint64_t expected_post_stages =
          static_cast<std::uint64_t>(post_sweeps) * nonterminal_calls;
      const std::uint64_t expected_chebyshev_stages =
          expected_pre_stages + expected_post_stages;
      const std::uint64_t expected_chebyshev_exchanges =
          expected_chebyshev_stages + nonterminal_calls;
      const std::uint64_t expected_halo_stages =
          expected_chebyshev_exchanges + 2U * nonterminal_calls +
          48U * terminal_calls + 1U;
      bool exact_level_stage_vector = true;
      for (std::size_t level = 0U; level + 1U < levels; ++level) {
        const std::uint64_t expected_calls =
            cycle == MgCycleKind::f_cycle
                ? static_cast<std::uint64_t>(level + 1U)
                : 1U;
        const std::uint64_t observed_level_stages =
            static_cast<std::uint64_t>(pre_sweeps + post_sweeps) *
            work.cycle_level_calls[level];
        exact_level_stage_vector &=
            observed_level_stages ==
            static_cast<std::uint64_t>(pre_sweeps + post_sweeps) *
                expected_calls;
      }
      const bool exact_asymmetric_seven_level_f =
          cycle != MgCycleKind::f_cycle || levels != 7U ||
          (expected_pre_stages == 21U && expected_post_stages == 42U &&
           expected_chebyshev_stages == 63U);
      const bool exact =
          static_cast<bool>(applied) && allocations == 0U &&
              exact_cycle_schedule(work, levels, cycle) &&
              exact_level_stage_vector &&
              exact_asymmetric_seven_level_f &&
              work.chebyshev_stages == expected_chebyshev_stages &&
              work.chebyshev_exchange_actions ==
                  expected_chebyshev_exchanges &&
              work.chebyshev_defect_actions ==
                  expected_chebyshev_exchanges &&
              work.chebyshev_retained_final_defect_actions ==
                  nonterminal_calls &&
              work.retained_final_defect_actions == nonterminal_calls &&
              halo.begin_calls == expected_halo_stages &&
              halo.finish_calls == expected_halo_stages;
      if (!exact && rank == 0) {
        std::cerr << "generic-work levels=" << levels << " cycle="
                  << static_cast<unsigned>(cycle) << " status="
                  << static_cast<unsigned>(applied.code) << '/'
                  << applied.detail << " alloc=" << allocations
                  << " stages=" << work.chebyshev_stages << '/'
                  << expected_chebyshev_stages << " exchanges="
                  << work.chebyshev_exchange_actions << '/'
                  << expected_chebyshev_exchanges << " defects="
                  << work.chebyshev_defect_actions << '/'
                  << expected_chebyshev_exchanges << " retained="
                  << work.chebyshev_retained_final_defect_actions << '/'
                  << work.retained_final_defect_actions << '/'
                  << nonterminal_calls << " halo=" << halo.begin_calls << '/'
                  << halo.finish_calls << '/' << expected_halo_stages << '\n';
      }
      passed &= expect(
          exact,
          rank,
          "generic V/3/3 and F/1/2 schedules match derived pre/post, retained-defect, transfer, terminal, Halo-stage, and zero-allocation work");
      if (!all_true(passed)) return false;
    }
  }
  return all_true(passed);
}

bool test_fixed_f_cycle_real_halo_failure_provenance(int rank, int size) {
  Fixture fixture;
  bool passed = expect(
      initialize(fixture, {}, 4.0, {24, 12, 8}, 2U,
                 MgCycleKind::f_cycle),
      rank, "fixed F-cycle real-Halo failure fixture initializes");
  if (!all_true(passed)) return false;
  NativeCartesianMgPlan plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          mg_spec(fixture, 1031U), services(fixture),
          coefficient_views(fixture), plan)) &&
          plan.level_count() == 2U,
      rank, "fixed F-cycle real-Halo failure plan compiles");
  if (!all_true(passed)) return false;

  struct Case {
    detail::HaloFailurePoint point;
    std::uint32_t expected_detail;
  };
  const std::array<Case, 2U> cases{{
      {detail::HaloFailurePoint::start, detail::halo_detail_start_failure},
      {detail::HaloFailurePoint::completion,
       detail::halo_detail_completion_failure},
  }};
  for (const Case& selected : cases) {
    fill(fixture.correction, -41.0);
    const std::uint64_t correction_before = checksum(fixture.correction.view);
    detail::set_halo_failure_for_test(selected.point, size - 1);
    Status failed{};
    std::size_t allocations = 0U;
    {
      allocation_observer::Guard guard;
      failed = plan.apply(as_const(fixture.residual.view),
                          fixture.correction.view, 0U);
      allocations = allocation_observer::count.load(std::memory_order_relaxed);
    }
    const int halo_lowest_failing_rank =
        fixture.halo.lowest_failing_rank();
    detail::clear_halo_failure_for_test();
    passed &= expect(
        failed.code == StatusCode::mpi_failure &&
            failed.detail == selected.expected_detail &&
            identical(packed(failed)) &&
            plan.lowest_failing_rank() == size - 1 &&
            halo_lowest_failing_rank == size - 1 && allocations == 0U &&
            checksum(fixture.correction.view) == correction_before,
        rank,
        "fixed F-cycle preserves real Halo start/completion provenance and caller output");
    if (!all_true(passed)) return false;

    const Status retried = plan.apply(as_const(fixture.residual.view),
                                      fixture.correction.view, 1U);
    passed &= expect(static_cast<bool>(retried) &&
                         plan.lowest_failing_rank() == -1 &&
                         finite(fixture.correction.view),
                     rank,
                     "fixed F-cycle retries after a real Halo failure on the same resources");
    if (!all_true(passed)) return false;
  }
  return all_true(passed);
}

bool test_rank_local_chebyshev_nonfinite_is_collective(int rank, int size) {
  Fixture fixture;
  bool passed = expect(initialize(fixture, {}, 4.0, {24, 12, 8}, 32U,
                                  MgCycleKind::f_cycle), rank,
                       "Chebyshev nonfinite RED fixture initializes");
  if (!all_true(passed)) {
    return false;
  }
  const MgOperatorClass certified =
      MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
  NativeCartesianMgSpec spec = mg_spec(
      fixture, 895U, MgPointSmootherKind::chebyshev_jacobi, certified);
  NativeCartesianMgPlan plan;
  passed &= expect(static_cast<bool>(NativeCartesianMgPlan::compile(
                       spec, services(fixture), coefficient_views(fixture),
                       plan)),
                   rank, "Chebyshev nonfinite RED plan compiles");
  if (!all_true(passed)) {
    return false;
  }

  const HaloRuntimeCounters healthy_halo_before =
      halo_counters(fixture);
  const Status healthy = plan.apply(as_const(fixture.residual.view),
                                    fixture.correction.view, 0U);
  const HaloRuntimeCounters healthy_halo = difference(
      halo_counters(fixture), healthy_halo_before);
  passed &= expect(static_cast<bool>(healthy), rank,
                   "Chebyshev nonfinite RED baseline apply succeeds");
  if (!all_true(passed)) {
    return false;
  }

  detail::set_mg_force_chebyshev_invalid_for_test(plan,
                                                   rank == size - 1);
  fill(fixture.correction, -31.0);
  const std::uint64_t correction_before = checksum(fixture.correction.view);
  const HaloRuntimeCounters failure_halo_before = halo_counters(fixture);
  Status applied{};
  std::size_t allocations = 0U;
  {
    allocation_observer::Guard guard;
    applied = plan.apply(as_const(fixture.residual.view),
                         fixture.correction.view, 1U);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const HaloRuntimeCounters failure_halo = difference(
      halo_counters(fixture), failure_halo_before);
  passed &= expect(applied.code == StatusCode::numerical_failure &&
                       identical(packed(applied)) &&
                       plan.lowest_failing_rank() == size - 1,
                   rank,
                   "rank-local Chebyshev nonfinite returns collectively");
  passed &= expect(checksum(fixture.correction.view) == correction_before,
                   rank,
                   "rank-local Chebyshev nonfinite publishes no output");
  passed &= expect(allocations == 0U, rank,
                   "rank-local Chebyshev nonfinite allocates no hot storage");
  passed &= expect(
      failure_halo.begin_calls + 1U == healthy_halo.begin_calls &&
          failure_halo.finish_calls + 1U == healthy_halo.finish_calls,
      rank,
      "rank-local Chebyshev nonfinite completes the fixed F-cycle halo schedule");
  return all_true(passed);
}

bool test_rank_local_line_pivot_failure_is_collective(int rank, int size) {
  Fixture point_fixture;
  bool passed = expect(initialize(point_fixture, {}, 4.0, {24, 12, 8}, 32U,
                                  MgCycleKind::f_cycle), rank,
                       "point-smoother reduction baseline initializes");
  if (!all_true(passed)) {
    return false;
  }
  NativeCartesianMgPlan point_plan;
  passed &= expect(static_cast<bool>(NativeCartesianMgPlan::compile(
                       mg_spec(point_fixture, 890U), services(point_fixture),
                       coefficient_views(point_fixture), point_plan)),
                   rank, "point-smoother reduction baseline compiles");
  if (!all_true(passed)) {
    return false;
  }
  const LinearReductionCounters point_reductions_before =
      point_fixture.reductions.counters();
  const Status point_applied = point_plan.apply(
      as_const(point_fixture.residual.view), point_fixture.correction.view, 0U);
  const LinearReductionCounters point_reductions = difference(
      point_fixture.reductions.counters(), point_reductions_before);

  Fixture fixture;
  passed &= expect(initialize(fixture, {}, 2.0, {24, 12, 8}, 32U,
                              MgCycleKind::f_cycle), rank,
                   "line-pivot fixture initializes with x line relaxation");
  if (!all_true(passed)) {
    return false;
  }
  NativeCartesianMgPlan plan;
  passed &= expect(static_cast<bool>(NativeCartesianMgPlan::compile(
                       mg_spec(fixture, 900U), services(fixture),
                       coefficient_views(fixture), plan)),
                   rank, "line-pivot MG plan compiles legal coefficients");
  MgLevelView finest{};
  passed &= expect(static_cast<bool>(plan.level(0U, finest)) &&
                       finest.line_axis_mask == 1U &&
                       fixture.patch.cells.x >= 2,
                   rank, "mutation reaches the Thomas x-line smoother");
  if (!all_true(passed)) {
    return false;
  }

  const HaloRuntimeCounters healthy_halo_before =
      fixture.halo.runtime_counters();
  const LinearReductionCounters healthy_reductions_before =
      fixture.reductions.counters();
  const Status healthy_applied = plan.apply(
      as_const(fixture.residual.view), fixture.correction.view, 0U);
  const HaloRuntimeCounters healthy_halo = difference(
      fixture.halo.runtime_counters(), healthy_halo_before);
  const LinearReductionCounters healthy_reductions = difference(
      fixture.reductions.counters(), healthy_reductions_before);
  passed &= expect(static_cast<bool>(point_applied) &&
                       static_cast<bool>(healthy_applied) &&
                       same(healthy_reductions, point_reductions),
                   rank,
                   "line smoothing adds no hot-path reduction to MG apply");
  if (!all_true(passed)) {
    return false;
  }

  if (rank == size - 1) {
    // Positive, finite coefficients remain compile-legal, but the first two
    // rows form the Thomas pivot 1 - (-1)*(-1) == 0 on one rank only.
    fixture.diagonal.view.unchecked({0, 0, 0}, 0U) = 1.0;
    fixture.diagonal.view.unchecked({1, 0, 0}, 0U) = 1.0;
  }
  LinearIdentity next = mg_spec(fixture, 900U).identity;
  next.numeric += 1U;
  next.hierarchy += 1U;
  next.fingerprint += 1U;
  const Status updated = plan.update_coefficients(
      next, MgCoefficientIdentity{942U, 952U, 0.1},
      coefficient_views(fixture), nullptr);
  passed &= expect(static_cast<bool>(updated), rank,
                   "rank-local pivot mutation remains a legal numeric refresh");
  if (!all_true(passed)) {
    return false;
  }

  fill(fixture.correction, -29.0);
  const std::uint64_t correction_before = checksum(fixture.correction.view);
  const HaloRuntimeCounters failure_halo_before =
      fixture.halo.runtime_counters();
  Status applied{};
  std::size_t allocations = 0U;
  {
    allocation_observer::Guard guard;
    applied = plan.apply(as_const(fixture.residual.view),
                         fixture.correction.view, 0U);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const HaloRuntimeCounters failure_halo = difference(
      fixture.halo.runtime_counters(), failure_halo_before);
  passed &= expect(applied.code == StatusCode::numerical_failure &&
                       identical(packed(applied)),
                   rank,
                   "rank-local Thomas pivot failure returns collectively");
  passed &= expect(checksum(fixture.correction.view) == correction_before,
                   rank, "rank-local Thomas pivot failure publishes no output");
  passed &= expect(allocations == 0U, rank,
                   "rank-local Thomas pivot failure allocates no hot storage");
  passed &= expect(
      failure_halo.begin_calls + 1U == healthy_halo.begin_calls &&
          failure_halo.finish_calls + 1U == healthy_halo.finish_calls,
      rank,
      "rank-local Thomas pivot failure completes the fixed F-cycle halo schedule");
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
                       "MG MPI RED runs at 1, 2, or 4 ranks");
  passed &= test_periodic_prolongation_wraps_coarse_z_endpoint(rank, size);
  passed &= test_cross_partition_periodic_prolongation_uses_complete_donors(
      rank, size);
  passed &= test_fused_point_smoother_reference_oracle(rank);
  passed &= test_fused_point_smoother_reference_oracle(rank, true);
  passed &= test_odd_periodic_point_smoother_fallback(rank);
  passed &= test_explicit_chebyshev_distributed_routes(rank);
  passed &= test_replicated_coarse_reference_oracle(rank, size);
  ReplicatedCoarseCase odd_case;
  odd_case.cells = {25, 11, 9};
  odd_case.cycle = MgCycleKind::f_cycle;
  passed &= compare_replicated_coarse_case(rank, size, odd_case, 970U);
  ReplicatedCoarseCase odd_periodic_case;
  odd_periodic_case.cells = {25, 12, 8};
  odd_periodic_case.topology = {true, false, false};
  odd_periodic_case.boundaries = {
      MgBoundaryKind::periodic, MgBoundaryKind::periodic,
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet,
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet};
  odd_periodic_case.cycle = MgCycleKind::f_cycle;
  passed &= compare_replicated_coarse_case(rank, size, odd_periodic_case, 980U);
  ReplicatedCoarseCase neumann_case;
  neumann_case.boundaries = {
      MgBoundaryKind::neumann, MgBoundaryKind::neumann,
      MgBoundaryKind::neumann, MgBoundaryKind::neumann,
      MgBoundaryKind::neumann, MgBoundaryKind::neumann};
  neumann_case.cycle = MgCycleKind::f_cycle;
  passed &= compare_replicated_coarse_case(rank, size, neumann_case, 990U);
  ReplicatedCoarseCase null_space_case;
  null_space_case.boundaries = neumann_case.boundaries;
  null_space_case.null_space = MgNullSpace::constant;
  null_space_case.cycle = MgCycleKind::f_cycle;
  passed &= compare_replicated_coarse_case(rank, size, null_space_case, 1000U);
  ReplicatedCoarseCase activity_case;
  activity_case.activity = true;
  activity_case.cycle = MgCycleKind::f_cycle;
  passed &= compare_replicated_coarse_case(rank, size, activity_case, 1010U);
  ReplicatedCoarseCase variable_case;
  variable_case.variable = true;
  variable_case.cycle = MgCycleKind::f_cycle;
  passed &= compare_replicated_coarse_case(rank, size, variable_case, 1020U);
  passed &= test_replicated_coarse_over_threshold_fallback(rank, size);
  passed &= test_collective_compile_and_ownership(rank, size);
  passed &= test_point_row_boundary_and_small_x_oracle(rank);
  passed &= test_chebyshev_point_row_reference_oracle(rank);
  passed &= test_boundary_topology_contract_is_bidirectional(rank);
  passed &= test_prepared_apply_equivalence(rank, size);
  passed &= test_progress_and_hot_reuse(rank);
  passed &= test_collective_failure_is_transactional(rank, size);
  passed &= test_default_v_cycle_is_explicit_v_cycle(rank);
  passed &= test_fixed_cycle_generic_scientific_work(rank, size);
  passed &= test_fixed_f_cycle_phase_failures_and_retry(rank, size);
  passed &= test_fixed_f_cycle_real_halo_failure_provenance(rank, size);
  if (size > 1) {
    passed &=
        test_rank_local_chebyshev_nonfinite_is_collective(rank, size);
    passed &= test_rank_local_line_pivot_failure_is_collective(rank, size);
  }
  passed = all_true(passed);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
