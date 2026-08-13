// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_parallel.hpp"

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
         left.applications == right.applications;
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
          after.bytes_unpacked - before.bytes_unpacked};
}

bool same(LinearReductionCounters left,
          LinearReductionCounters right) noexcept {
  return left.calls == right.calls && left.scalars == right.scalars &&
         left.logical_bytes == right.logical_bytes &&
         left.tree_messages == right.tree_messages;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {24, 12, 8};
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
  MgWorkspaceRequirements requirements{};
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  double anisotropy_threshold{4.0};
};

bool initialize(Fixture& fixture, HaloTopology topology = {},
                double anisotropy_threshold = 4.0) {
  fixture.anisotropy_threshold = anisotropy_threshold;
  if (!CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, mesh_spec(), GeometryBudget{}, fixture.geometry,
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

NativeCartesianMgSpec mg_spec(const Fixture& fixture,
                              std::uint64_t salt = 0U) noexcept {
  NativeCartesianMgSpec result;
  result.communicator = MPI_COMM_WORLD;
  result.geometry = &fixture.geometry;
  result.patch = fixture.patch;
  result.identity = identity(salt);
  result.coefficients = {41U + salt, 51U + salt, 0.0};
  result.policy.anisotropy_threshold = fixture.anisotropy_threshold;
  result.policy.coefficient_change_rebuild_ratio = 0.25;
  return result;
}

MgRuntimeServices services(Fixture& fixture) noexcept {
  return {&fixture.halo, &fixture.reductions, &fixture.workspace,
          {fixture.coarse_halo_pointers.data(),
           fixture.coarse_halo_pointers.size()}};
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

bool test_progress_and_hot_reuse(int rank) {
  Fixture fixture;
  bool passed = expect(initialize(fixture), rank,
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

bool test_rank_local_line_pivot_failure_is_collective(int rank, int size) {
  Fixture point_fixture;
  bool passed = expect(initialize(point_fixture), rank,
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
  passed &= expect(initialize(fixture, {}, 2.0), rank,
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
      "rank-local Thomas pivot failure completes the V-cycle halo schedule");
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
  passed &= test_collective_compile_and_ownership(rank, size);
  passed &= test_boundary_topology_contract_is_bidirectional(rank);
  passed &= test_progress_and_hot_reuse(rank);
  passed &= test_collective_failure_is_transactional(rank, size);
  if (size > 1) {
    passed &= test_rank_local_line_pivot_failure_is_collective(rank, size);
  }
  passed = all_true(passed);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
