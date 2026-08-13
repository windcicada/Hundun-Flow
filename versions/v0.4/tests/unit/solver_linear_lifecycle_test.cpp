// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"
#include "hundun/v04_field.hpp"
#include "hundun/v04_linear.hpp"

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>

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

constexpr Int3 kMaximumShape{9, 5, 3};
constexpr std::uint8_t kGhostWidth = 1U;
constexpr std::uint32_t kRestart = 8U;
constexpr RevisionToken kExecutionRevision = 41U;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(LinearLifecycleCounters left,
          LinearLifecycleCounters right) noexcept {
  return left.symbolic_builds == right.symbolic_builds &&
         left.numeric_refills == right.numeric_refills &&
         left.hierarchy_refreshes == right.hierarchy_refreshes &&
         left.hierarchy_rebuilds == right.hierarchy_rebuilds &&
         left.workspace_bindings == right.workspace_bindings &&
         left.workspace_replacements == right.workspace_replacements;
}

SymbolicSpec symbolic_spec(RevisionToken topology = 11U) noexcept {
  return SymbolicSpec{17U,
                      LinearLocation::cell,
                      topology,
                      12U,
                      13U,
                      14U,
                      15U,
                      LinearBackend::native_cartesian};
}

CoefficientRevisions coefficient_revisions(
    RevisionToken diagonal = 21U) noexcept {
  return CoefficientRevisions{diagonal, 22U, 23U, 24U, 25U, 26U};
}

HierarchyPolicyIdentity hierarchy_policy(
    RevisionToken epoch = 33U) noexcept {
  return HierarchyPolicyIdentity{31U, 32U, epoch,
                                 LinearBackend::native_cartesian};
}

bool make_requirements(LinearWorkspaceRequirements& out) {
  return static_cast<bool>(make_linear_workspace_requirements(
      LinearAlgorithm::fgmres, kMaximumShape, kGhostWidth, kRestart,
      ReductionMode::mpi_allreduce, kExecutionRevision, out));
}

bool test_workspace_requirement_contracts() {
  bool passed = true;
  LinearWorkspaceRequirements pcg{};
  passed &= expect(static_cast<bool>(make_linear_workspace_requirements(
                       LinearAlgorithm::pcg, kMaximumShape, kGhostWidth, 0U,
                       ReductionMode::mpi_allreduce, kExecutionRevision, pcg)),
                   "PCG workspace requirements compile");
  passed &= expect(pcg.vector_slots == 5U && pcg.reduction_capacity >= 2U &&
                       same(pcg.maximum_shape, kMaximumShape) &&
                       pcg.maximum_restart == 0U,
                   "PCG reserves exactly five persistent vectors");

  LinearWorkspaceRequirements bicgstab{};
  passed &= expect(static_cast<bool>(make_linear_workspace_requirements(
                       LinearAlgorithm::bicgstab, kMaximumShape, kGhostWidth,
                       0U, ReductionMode::reproducible_tree,
                       kExecutionRevision, bicgstab)),
                   "BiCGSTAB workspace requirements compile");
  passed &= expect(bicgstab.vector_slots == 10U &&
                       bicgstab.reduction_capacity >= 2U &&
                       bicgstab.maximum_restart == 0U,
                   "BiCGSTAB reserves exactly ten persistent vectors");

  LinearWorkspaceRequirements fgmres{};
  passed &= expect(make_requirements(fgmres),
                   "FGMRES workspace requirements compile");
  passed &= expect(fgmres.vector_slots == 2U * kRestart + 3U &&
                       fgmres.maximum_restart == kRestart &&
                       fgmres.scalar_doubles > 0U &&
                       fgmres.reduction_capacity >= kRestart + 1U,
                   "FGMRES reserves x, V, Z, and scalar recurrence storage");

  const LinearWorkspaceRequirements preserved = pcg;
  const Status pcg_restart = make_linear_workspace_requirements(
      LinearAlgorithm::pcg, kMaximumShape, kGhostWidth, 1U,
      ReductionMode::mpi_allreduce, kExecutionRevision, pcg);
  passed &= expect(pcg_restart.code == StatusCode::invalid_plan &&
                       pcg.vector_slots == preserved.vector_slots &&
                       same(pcg.maximum_shape, preserved.maximum_shape),
                   "PCG rejects a restart dimension atomically");

  LinearWorkspaceRequirements ignored{};
  passed &= expect(
      make_linear_workspace_requirements(
          LinearAlgorithm::bicgstab, kMaximumShape, kGhostWidth, 1U,
          ReductionMode::mpi_allreduce, kExecutionRevision, ignored)
              .code == StatusCode::invalid_plan,
      "BiCGSTAB rejects a restart dimension");
  passed &= expect(
      make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, kMaximumShape, kGhostWidth, 0U,
          ReductionMode::mpi_allreduce, kExecutionRevision, ignored)
              .code == StatusCode::invalid_plan,
      "FGMRES rejects a zero restart dimension");
  passed &= expect(
      make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, kMaximumShape, kGhostWidth, 127U,
          ReductionMode::mpi_allreduce, kExecutionRevision, ignored)
              .code == StatusCode::invalid_plan,
      "FGMRES rejects a vector-slot count beyond uint8 field capacity");
  passed &= expect(
      make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, Int3{0, 5, 3}, kGhostWidth, kRestart,
          ReductionMode::mpi_allreduce, kExecutionRevision, ignored)
              .code == StatusCode::invalid_plan,
      "workspace planning rejects non-positive maximum shape");
  passed &= expect(
      make_linear_workspace_requirements(
          LinearAlgorithm::fgmres,
          Int3{std::numeric_limits<std::int32_t>::max(),
               std::numeric_limits<std::int32_t>::max(),
               std::numeric_limits<std::int32_t>::max()},
          kGhostWidth, 126U, ReductionMode::mpi_allreduce,
          kExecutionRevision, ignored)
              .code == StatusCode::invalid_plan,
      "workspace planning rejects checked-size overflow");

  LinearWorkspaceRequirements maximum{};
  passed &= expect(
      static_cast<bool>(make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, Int3{1, 1, 1}, 0U, 126U,
          ReductionMode::mpi_allreduce, kExecutionRevision, maximum)) &&
          maximum.vector_slots ==
              std::numeric_limits<std::uint8_t>::max(),
      "the maximum legal FGMRES restart uses all 255 uint8 components");
  FieldRegistry maximum_registry;
  LinearWorkspaceFieldIds maximum_ids{};
  FieldSchema maximum_schema;
  passed &= expect(
      static_cast<bool>(register_linear_workspace(
          maximum_registry, "maximum_linear", maximum, maximum_ids)) &&
          static_cast<bool>(maximum_registry.freeze(maximum_schema)) &&
          maximum_schema[maximum_ids.vector_bundle].components ==
              std::numeric_limits<std::uint8_t>::max(),
      "the field registry preserves the full uint8 component capacity");

  FieldRegistry collision_registry;
  FieldId collision{};
  passed &= expect(
      static_cast<bool>(collision_registry.declare_field(
          "atomic_linear.scalars", 1U, 0U, collision)),
      "atomic registration collision fixture declares");
  const PlanFingerprint registry_before = collision_registry.fingerprint();
  LinearWorkspaceFieldIds preserved_ids{91U, 92U};
  const Status collision_status = register_linear_workspace(
      collision_registry, "atomic_linear", pcg, preserved_ids);
  passed &= expect(
      collision_status.code == StatusCode::invalid_plan &&
          collision_registry.fingerprint() == registry_before &&
          preserved_ids.vector_bundle == 91U &&
          preserved_ids.scalar_buffer == 92U,
      "two-field workspace registration rejects atomically after a late collision");
  return passed;
}

bool test_runtime_lifetime_views() {
  FieldRegistry registry;
  FieldId state{};
  FieldId pending{};
  FieldId persistent{};
  FieldId scratch{};
  FieldSchema schema;
  bool passed = expect(
      static_cast<bool>(registry.declare_field("state", 1U, 0U, state)) &&
          static_cast<bool>(registry.declare_field("pending", 1U, 0U,
                                                   pending)) &&
          static_cast<bool>(registry.declare_field("persistent", 1U, 0U,
                                                   persistent)) &&
          static_cast<bool>(registry.declare_field("scratch", 1U, 0U,
                                                   scratch)) &&
          static_cast<bool>(registry.freeze(schema)),
      "runtime lifetime schema freezes");
  if (!passed) {
    return false;
  }
  constexpr Int3 shape{3, 2, 1};
  const std::array requests{
      ArenaFieldRequest{state, shape, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{pending, shape, FieldPlacement{0U},
                        FieldLifetime::pending_cache},
      ArenaFieldRequest{persistent, shape, FieldPlacement{0U},
                        FieldLifetime::persistent_workspace},
      ArenaFieldRequest{scratch, shape, FieldPlacement{0U},
                        FieldLifetime::step_scratch},
  };
  ArenaLayout layout;
  StateLayers layers;
  passed &= expect(
      static_cast<bool>(ArenaLayout::compile(
          schema, Span<const ArenaFieldRequest>{requests.data(),
                                                requests.size()},
          layout)) &&
          static_cast<bool>(StateLayers::allocate(layout, layers)),
      "all runtime lifetime classes allocate in one arena");
  if (!passed) {
    return false;
  }

  const std::array runtime_ids{pending, persistent, scratch};
  constexpr std::array runtime_lifetimes{
      FieldLifetime::pending_cache, FieldLifetime::persistent_workspace,
      FieldLifetime::step_scratch};
  RevisionToken prior{};
  for (std::size_t index = 0U; index < runtime_ids.size(); ++index) {
    FieldView view;
    const RevisionToken revision =
        layers.runtime_revision(runtime_lifetimes[index], runtime_ids[index]);
    passed &= expect(
        revision != 0U && revision != prior &&
            static_cast<bool>(layers.runtime_view(
                runtime_lifetimes[index], runtime_ids[index], view)) &&
            view.field == runtime_ids[index] && view.replica == 0U &&
            view.revision == revision,
        "each non-state lifetime publishes one distinct single-replica view");
    prior = revision;
  }
  passed &= expect(layers.counters().aligned_payload_allocations == 1U,
                   "runtime fields add no second full-field allocation");

  FieldView sentinel;
  sentinel.base = reinterpret_cast<double*>(std::uintptr_t{1U});
  sentinel.field = std::numeric_limits<FieldId>::max();
  const Status wrong_runtime = layers.runtime_view(
      FieldLifetime::step_scratch, persistent, sentinel);
  passed &= expect(
      wrong_runtime.code == StatusCode::invalid_plan &&
          layers.runtime_revision(FieldLifetime::step_scratch, persistent) ==
              0U &&
          sentinel.base == reinterpret_cast<double*>(std::uintptr_t{1U}) &&
          sentinel.field == std::numeric_limits<FieldId>::max(),
      "runtime view rejects a lifetime mismatch without publishing output");
  FieldView wrong_state;
  passed &= expect(
      layers.view(StateRole::accepted_n, persistent, wrong_state).code ==
              StatusCode::invalid_plan &&
          layers.runtime_view(FieldLifetime::state_layer, state, wrong_state)
                  .code == StatusCode::invalid_plan,
      "state and runtime view APIs reject fields from the other lifetime class");
  return passed;
}

bool test_symbolic_and_numeric_lifecycle() {
  LinearLifecycleCounters counters{};
  SymbolicPlan symbolic;
  bool passed = expect(
      static_cast<bool>(SymbolicPlan::compile(symbolic_spec(), symbolic,
                                               &counters)),
      "initial symbolic structure compiles");
  const PlanFingerprint first_symbolic = symbolic.fingerprint();
  const RevisionToken first_symbolic_generation = symbolic.generation();
  passed &= expect(symbolic.valid() && first_symbolic != 0U &&
                       first_symbolic_generation != 0U &&
                       counters.symbolic_builds == 1U,
                   "initial symbolic build publishes one structural identity");

  passed &= expect(
      static_cast<bool>(SymbolicPlan::compile(symbolic_spec(), symbolic,
                                               &counters)) &&
          symbolic.fingerprint() == first_symbolic &&
          symbolic.generation() == first_symbolic_generation &&
          counters.symbolic_builds == 1U,
      "replaying an identical symbolic specification is an exact no-op");

  NumericState numeric;
  passed &= expect(static_cast<bool>(numeric.refill(
                       symbolic, coefficient_revisions(), 501U, &counters)),
                   "initial numeric coefficients publish");
  const PlanFingerprint first_numeric = numeric.fingerprint();
  const RevisionToken first_numeric_generation = numeric.generation();
  passed &= expect(numeric.valid_for(symbolic) && first_numeric != 0U &&
                       first_numeric_generation != 0U &&
                       counters.numeric_refills == 1U,
                   "numeric state records its symbolic authority");

  passed &= expect(
      static_cast<bool>(numeric.refill(symbolic, coefficient_revisions(),
                                       501U, &counters)) &&
          numeric.fingerprint() == first_numeric &&
          numeric.generation() == first_numeric_generation &&
          counters.numeric_refills == 1U,
      "identical coefficient revisions and content are an exact no-op");

  const LinearLifecycleCounters before_mutation = counters;
  const Status content_mutation = numeric.refill(
      symbolic, coefficient_revisions(), 502U, &counters);
  passed &= expect(content_mutation.code == StatusCode::invalid_plan &&
                       numeric.fingerprint() == first_numeric &&
                       numeric.generation() == first_numeric_generation &&
                       same(counters, before_mutation),
                   "changed numeric content under unchanged revisions is rejected atomically");

  passed &= expect(static_cast<bool>(numeric.refill(
                       symbolic, coefficient_revisions(27U), 502U,
                       &counters)) &&
                       numeric.fingerprint() != first_numeric &&
                       numeric.generation() > first_numeric_generation &&
                       counters.numeric_refills == 2U &&
                       counters.symbolic_builds == 1U,
                   "coefficient-only change performs exactly one refill without rebuilding symbolic state");

  HierarchyState hierarchy;
  passed &= expect(static_cast<bool>(hierarchy.update(
                       symbolic, numeric, hierarchy_policy(),
                       HierarchyUpdate::rebuild, 701U, &counters)),
                   "initial hierarchy build publishes");
  const PlanFingerprint first_hierarchy = hierarchy.fingerprint();
  const RevisionToken first_hierarchy_generation = hierarchy.generation();
  passed &= expect(hierarchy.valid_for(symbolic, numeric) &&
                       first_hierarchy != 0U &&
                       counters.hierarchy_rebuilds == 1U &&
                       counters.hierarchy_refreshes == 0U,
                   "initial hierarchy records symbolic and numeric identities");

  passed &= expect(
      static_cast<bool>(hierarchy.update(
          symbolic, numeric, hierarchy_policy(), HierarchyUpdate::retain,
          701U, &counters)) &&
          hierarchy.fingerprint() == first_hierarchy &&
          hierarchy.generation() == first_hierarchy_generation &&
          counters.hierarchy_rebuilds == 1U,
      "retain is a no-op while all hierarchy dependencies are unchanged");

  passed &= expect(static_cast<bool>(numeric.refill(
                       symbolic, coefficient_revisions(28U), 503U,
                       &counters)),
                   "a second coefficient revision refills numeric state");
  const PlanFingerprint stale_hierarchy = hierarchy.fingerprint();
  const RevisionToken stale_hierarchy_generation = hierarchy.generation();
  const LinearLifecycleCounters before_illegal_retain = counters;
  const Status illegal_retain = hierarchy.update(
      symbolic, numeric, hierarchy_policy(), HierarchyUpdate::retain, 701U,
      &counters);
  passed &= expect(illegal_retain.code == StatusCode::invalid_plan &&
                       hierarchy.fingerprint() == stale_hierarchy &&
                       hierarchy.generation() == stale_hierarchy_generation &&
                       same(counters, before_illegal_retain),
                   "retain rejects a changed numeric identity atomically");

  passed &= expect(static_cast<bool>(hierarchy.update(
                       symbolic, numeric, hierarchy_policy(),
                       HierarchyUpdate::refresh_numeric, 702U, &counters)) &&
                       hierarchy.valid_for(symbolic, numeric) &&
                       counters.hierarchy_refreshes == 1U &&
                       counters.hierarchy_rebuilds == 1U,
                   "below-policy coefficient change refreshes hierarchy numeric data only");
  const PlanFingerprint refreshed_hierarchy = hierarchy.fingerprint();
  const RevisionToken refreshed_generation = hierarchy.generation();

  const LinearLifecycleCounters before_policy_mismatch = counters;
  const Status policy_mismatch = hierarchy.update(
      symbolic, numeric, hierarchy_policy(34U),
      HierarchyUpdate::refresh_numeric, 703U, &counters);
  passed &= expect(policy_mismatch.code == StatusCode::invalid_plan &&
                       hierarchy.fingerprint() == refreshed_hierarchy &&
                       hierarchy.generation() == refreshed_generation &&
                       same(counters, before_policy_mismatch),
                   "hierarchy policy change cannot masquerade as numeric refresh");

  passed &= expect(static_cast<bool>(hierarchy.update(
                       symbolic, numeric, hierarchy_policy(34U),
                       HierarchyUpdate::rebuild, 703U, &counters)) &&
                       counters.hierarchy_rebuilds == 2U &&
                       counters.hierarchy_refreshes == 1U,
                   "policy epoch change performs exactly one hierarchy rebuild");
  const PlanFingerprint rebuilt_hierarchy = hierarchy.fingerprint();
  const RevisionToken rebuilt_generation = hierarchy.generation();
  const LinearLifecycleCounters before_hierarchy_mutation = counters;
  const Status hierarchy_mutation = hierarchy.update(
      symbolic, numeric, hierarchy_policy(34U), HierarchyUpdate::retain, 704U,
      &counters);
  passed &= expect(hierarchy_mutation.code == StatusCode::invalid_plan &&
                       hierarchy.fingerprint() == rebuilt_hierarchy &&
                       hierarchy.generation() == rebuilt_generation &&
                       same(counters, before_hierarchy_mutation),
                   "changed hierarchy content under unchanged identity is rejected atomically");

  passed &= expect(
      static_cast<bool>(SymbolicPlan::compile(symbolic_spec(16U), symbolic,
                                               &counters)) &&
          symbolic.fingerprint() != first_symbolic &&
          symbolic.generation() > first_symbolic_generation &&
          counters.symbolic_builds == 2U,
      "topology change replaces symbolic structure exactly once");
  passed &= expect(!numeric.valid_for(symbolic) &&
                       !hierarchy.valid_for(symbolic, numeric),
                   "structural replacement makes prior numeric and hierarchy states stale");
  return passed;
}

struct WorkspaceFixture {
  LinearWorkspaceRequirements requirements;
  LinearWorkspaceFieldIds ids;
  FieldId anchor{};
  FieldSchema schema;
  ArenaLayout layout;
  StateLayers layers;
  FieldView vector_bundle{};
  FieldView scalar_buffer{};
};

bool make_workspace_fixture(WorkspaceFixture& fixture) {
  if (!make_requirements(fixture.requirements)) {
    return false;
  }
  FieldRegistry registry;
  if (!registry.declare_field("linear_anchor", 1U, 0U, fixture.anchor) ||
      !register_linear_workspace(registry, "pressure_linear",
                                 fixture.requirements, fixture.ids) ||
      !registry.freeze(fixture.schema)) {
    return false;
  }
  LinearWorkspaceFieldIds late;
  if (register_linear_workspace(registry, "late_linear", fixture.requirements,
                                late)
          .code != StatusCode::invalid_plan) {
    return false;
  }

  if (fixture.requirements.scalar_doubles == 0U ||
      fixture.requirements.scalar_doubles >
          static_cast<std::size_t>(
              std::numeric_limits<std::int32_t>::max())) {
    return false;
  }
  const std::array requests{
      ArenaFieldRequest{fixture.anchor, kMaximumShape, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{fixture.ids.vector_bundle, kMaximumShape,
                        FieldPlacement{0U},
                        FieldLifetime::persistent_workspace},
      ArenaFieldRequest{
          fixture.ids.scalar_buffer,
          Int3{static_cast<std::int32_t>(fixture.requirements.scalar_doubles),
               1, 1},
          FieldPlacement{0U}, FieldLifetime::persistent_workspace},
  };
  if (!ArenaLayout::compile(
          fixture.schema,
          Span<const ArenaFieldRequest>{requests.data(), requests.size()},
          fixture.layout) ||
      !StateLayers::allocate(fixture.layout, fixture.layers) ||
      !fixture.layers.runtime_view(FieldLifetime::persistent_workspace,
                                   fixture.ids.vector_bundle,
                                   fixture.vector_bundle) ||
      !fixture.layers.runtime_view(FieldLifetime::persistent_workspace,
                                   fixture.ids.scalar_buffer,
                                   fixture.scalar_buffer)) {
    return false;
  }
  return true;
}

struct WorkspaceSnapshot {
  PlanFingerprint fingerprint{};
  RevisionToken binding{};
  std::uintptr_t vectors{};
  std::uintptr_t scalars{};
};

WorkspaceSnapshot snapshot(const SolverWorkspace& workspace) noexcept {
  return WorkspaceSnapshot{workspace.fingerprint(),
                           workspace.binding_identity(),
                           workspace.vector_storage_address(),
                           workspace.scalar_storage_address()};
}

bool same(WorkspaceSnapshot left, WorkspaceSnapshot right) noexcept {
  return left.fingerprint == right.fingerprint &&
         left.binding == right.binding && left.vectors == right.vectors &&
         left.scalars == right.scalars;
}

bool test_workspace_internal_overlap_is_rejected_atomically() {
  WorkspaceFixture fixture;
  bool passed = expect(
      make_workspace_fixture(fixture),
      "workspace-overlap fixture allocates two non-overlapping arena fields");
  if (!passed) {
    return false;
  }

  FieldView exact_overlap = fixture.scalar_buffer;
  exact_overlap.base = fixture.vector_bundle.base;
  SolverWorkspace exact_workspace;
  LinearLifecycleCounters exact_counters{};
  exact_counters.workspace_bindings = 7U;
  const WorkspaceSnapshot exact_before = snapshot(exact_workspace);
  const LinearLifecycleCounters exact_counters_before = exact_counters;
  const Status exact = SolverWorkspace::bind(
      fixture.requirements, fixture.vector_bundle, exact_overlap,
      exact_workspace, &exact_counters);
  passed &= expect(exact.code == StatusCode::invalid_plan &&
                       same(snapshot(exact_workspace), exact_before) &&
                       same(exact_counters, exact_counters_before),
                   "an exactly overlapping scalar/vector binding is rejected atomically");

  FieldView shifted_overlap = fixture.scalar_buffer;
  shifted_overlap.base = fixture.vector_bundle.base + 1U;
  SolverWorkspace shifted_workspace;
  LinearLifecycleCounters shifted_counters{};
  shifted_counters.workspace_bindings = 11U;
  const WorkspaceSnapshot shifted_before = snapshot(shifted_workspace);
  const LinearLifecycleCounters shifted_counters_before = shifted_counters;
  const Status shifted = SolverWorkspace::bind(
      fixture.requirements, fixture.vector_bundle, shifted_overlap,
      shifted_workspace, &shifted_counters);
  passed &= expect(shifted.code == StatusCode::invalid_plan &&
                       same(snapshot(shifted_workspace), shifted_before) &&
                       same(shifted_counters, shifted_counters_before),
                   "a shifted partial scalar/vector overlap is rejected atomically");

  SolverWorkspace valid_workspace;
  passed &= expect(
      static_cast<bool>(SolverWorkspace::bind(
          fixture.requirements, fixture.vector_bundle, fixture.scalar_buffer,
          valid_workspace)) &&
          valid_workspace.vector_storage_address() ==
              reinterpret_cast<std::uintptr_t>(fixture.vector_bundle.base) &&
          valid_workspace.scalar_storage_address() ==
              reinterpret_cast<std::uintptr_t>(fixture.scalar_buffer.base),
      "non-overlapping workspace fields in the same arena remain valid");
  return passed;
}

bool test_arena_workspace_binding_and_borrowing() {
  WorkspaceFixture fixture;
  bool passed = expect(make_workspace_fixture(fixture),
                       "workspace fields register before freeze and bind from the arena");
  if (!passed) {
    return false;
  }
  passed &= expect(fixture.vector_bundle.components ==
                           fixture.requirements.vector_slots &&
                       fixture.vector_bundle.ghosts.x >= kGhostWidth &&
                       fixture.scalar_buffer.interior.x >=
                           static_cast<std::int32_t>(
                               fixture.requirements.scalar_doubles),
                   "registered field shapes provide the exact planned capacities");
  passed &= expect(
      fixture.layers.runtime_revision(FieldLifetime::persistent_workspace,
                                      fixture.ids.vector_bundle) ==
              fixture.vector_bundle.revision &&
          fixture.layers.runtime_revision(FieldLifetime::persistent_workspace,
                                          fixture.ids.scalar_buffer) ==
              fixture.scalar_buffer.revision,
      "StateLayers publishes one runtime revision per workspace field");

  const StateLayers& const_layers = fixture.layers;
  ConstFieldView const_vectors;
  passed &= expect(
      static_cast<bool>(const_layers.runtime_view(
          FieldLifetime::persistent_workspace, fixture.ids.vector_bundle,
          const_vectors)) &&
          const_vectors.base == fixture.vector_bundle.base &&
          const_vectors.revision == fixture.vector_bundle.revision,
      "const runtime workspace view borrows the same single replica");

  LinearLifecycleCounters counters;
  SolverWorkspace workspace;
  passed &= expect(
      static_cast<bool>(SolverWorkspace::bind(
          fixture.requirements, fixture.vector_bundle, fixture.scalar_buffer,
          workspace, &counters)),
      "solver workspace binds non-owning arena views");
  const WorkspaceSnapshot initial = snapshot(workspace);
  passed &= expect(workspace.valid_for(fixture.requirements) &&
                       initial.fingerprint != 0U && initial.binding != 0U &&
                       initial.vectors == reinterpret_cast<std::uintptr_t>(
                                              fixture.vector_bundle.base) &&
                       initial.scalars == reinterpret_cast<std::uintptr_t>(
                                              fixture.scalar_buffer.base) &&
                       counters.workspace_bindings == 1U &&
                       counters.workspace_replacements == 0U,
                   "initial binding records stable external storage without replacement");

  passed &= expect(
      static_cast<bool>(SolverWorkspace::bind(
          fixture.requirements, fixture.vector_bundle, fixture.scalar_buffer,
          workspace, &counters)) &&
          same(snapshot(workspace), initial) &&
          counters.workspace_bindings == 1U &&
          counters.workspace_replacements == 0U,
      "revalidating the same workspace is an exact no-op");

  SymbolicPlan topology;
  passed &= expect(static_cast<bool>(
                       SymbolicPlan::compile(symbolic_spec(), topology,
                                             &counters)) &&
                       static_cast<bool>(SymbolicPlan::compile(
                           symbolic_spec(16U), topology, &counters)) &&
                       static_cast<bool>(SolverWorkspace::bind(
                           fixture.requirements, fixture.vector_bundle,
                           fixture.scalar_buffer, workspace, &counters)) &&
                       same(snapshot(workspace), initial) &&
                       counters.workspace_replacements == 0U,
                   "topology replacement with compatible capacity preserves workspace addresses");

  LinearWorkspaceRequirements larger{};
  passed &= expect(static_cast<bool>(make_linear_workspace_requirements(
                       LinearAlgorithm::fgmres, kMaximumShape, kGhostWidth,
                       kRestart + 1U, ReductionMode::mpi_allreduce,
                       kExecutionRevision, larger)),
                   "a larger valid workspace requirement compiles");
  const LinearLifecycleCounters before_capacity = counters;
  const Status insufficient = SolverWorkspace::bind(
      larger, fixture.vector_bundle, fixture.scalar_buffer, workspace,
      &counters);
  passed &= expect(insufficient.code == StatusCode::invalid_plan &&
                       same(snapshot(workspace), initial) &&
                       same(counters, before_capacity),
                   "capacity overflow rejects binding atomically without moving storage");

  WorkspaceFixture foreign;
  if (!expect(make_workspace_fixture(foreign),
              "an independent arena workspace fixture allocates")) {
    return false;
  }
  const LinearLifecycleCounters before_foreign = counters;
  const Status mixed_storage = SolverWorkspace::bind(
      fixture.requirements, foreign.vector_bundle, fixture.scalar_buffer,
      workspace, &counters);
  passed &= expect(mixed_storage.code == StatusCode::invalid_plan &&
                       same(snapshot(workspace), initial) &&
                       same(counters, before_foreign),
                   "cross-storage workspace views are rejected atomically");

  FieldView stale_vectors = fixture.vector_bundle;
  ++stale_vectors.revision;
  const LinearLifecycleCounters before_stale = counters;
  const Status stale = SolverWorkspace::bind(
      fixture.requirements, stale_vectors, fixture.scalar_buffer, workspace,
      &counters);
  passed &= expect(stale.code == StatusCode::invalid_plan &&
                       same(snapshot(workspace), initial) &&
                       same(counters, before_stale),
                   "a stale rebind view is rejected without changing the live workspace");

  FieldView first = workspace.vector(0U, Int3{7, 4, 2});
  FieldView second = workspace.vector(1U, Int3{7, 4, 2});
  FieldView last = workspace.vector(
      static_cast<std::uint8_t>(fixture.requirements.vector_slots - 1U),
      kMaximumShape);
  passed &= expect(first.base == fixture.vector_bundle.base &&
                       first.components == 1U &&
                       first.revision != 0U && second.revision != 0U &&
                       first.revision != second.revision &&
                       same(first.interior, Int3{7, 4, 2}) &&
                       last.base ==
                           fixture.vector_bundle.base +
                               (fixture.requirements.vector_slots - 1U) *
                                   fixture.vector_bundle.component_stride,
                   "vector borrowing selects one persistent SoA slot without allocation");
  const RevisionToken first_before_revise = first.revision;
  const RevisionToken second_before_revise = second.revision;
  const RevisionToken last_before_revise = last.revision;
  passed &= expect(
      static_cast<bool>(workspace.revise_vector(0U)) &&
          workspace.vector(0U, Int3{7, 4, 2}).revision >
              first_before_revise &&
          workspace.vector(0U, Int3{7, 4, 2}).revision >
              last_before_revise &&
          workspace.vector(1U, Int3{7, 4, 2}).revision ==
              second_before_revise,
      "explicit vector publication advances only the written slot revision");
  passed &= expect(workspace.vector(fixture.requirements.vector_slots,
                                    kMaximumShape)
                           .base == nullptr &&
                       workspace.vector(0U, Int3{kMaximumShape.x + 1,
                                                kMaximumShape.y,
                                                kMaximumShape.z})
                               .base == nullptr,
                   "vector borrowing rejects slot and active-shape overflow");

  Span<double> scalars = workspace.scalars(3U, 7U);
  passed &= expect(scalars.data == fixture.scalar_buffer.base + 3U &&
                       scalars.size == 7U,
                   "scalar recurrence borrowing returns one bounded contiguous span");
  passed &= expect(
      workspace
              .scalars(fixture.requirements.scalar_doubles - 1U, 2U)
              .data == nullptr,
      "scalar recurrence borrowing rejects offset-plus-count overflow");
  return passed;
}

bool build_consistent_lifecycle(SymbolicPlan& symbolic, NumericState& numeric,
                                HierarchyState& hierarchy,
                                LinearLifecycleCounters& counters) {
  return SymbolicPlan::compile(symbolic_spec(), symbolic, &counters) &&
         numeric.refill(symbolic, coefficient_revisions(), 801U, &counters) &&
         hierarchy.update(symbolic, numeric, hierarchy_policy(),
                          HierarchyUpdate::rebuild, 802U, &counters);
}

bool test_composed_identity_and_hot_noops() {
  WorkspaceFixture fixture;
  if (!make_workspace_fixture(fixture)) {
    return expect(false, "composed identity workspace fixture initializes");
  }
  LinearLifecycleCounters counters{};
  SymbolicPlan symbolic;
  NumericState numeric;
  HierarchyState hierarchy;
  SolverWorkspace workspace;
  bool passed = expect(
      build_consistent_lifecycle(symbolic, numeric, hierarchy, counters) &&
          SolverWorkspace::bind(fixture.requirements, fixture.vector_bundle,
                                fixture.scalar_buffer, workspace, &counters),
      "consistent four-layer lifecycle initializes");
  if (!passed) {
    return false;
  }

  const LinearIdentity valid =
      compose_linear_identity(symbolic, numeric, hierarchy, workspace);
  passed &= expect(valid.fingerprint != 0U,
                   "four mutually compatible layers compose a nonzero identity");

  passed &= expect(static_cast<bool>(numeric.refill(
                       symbolic, coefficient_revisions(29U), 803U,
                       &counters)),
                   "numeric mutation for stale compose test publishes");
  const LinearIdentity stale =
      compose_linear_identity(symbolic, numeric, hierarchy, workspace);
  passed &= expect(stale.fingerprint == 0U,
                   "compose rejects a hierarchy stale against current numeric state");
  passed &= expect(static_cast<bool>(hierarchy.update(
                       symbolic, numeric, hierarchy_policy(),
                       HierarchyUpdate::refresh_numeric, 804U, &counters)),
                   "hierarchy refresh restores a consistent identity");
  const LinearIdentity refreshed =
      compose_linear_identity(symbolic, numeric, hierarchy, workspace);
  passed &= expect(refreshed.fingerprint != 0U &&
                       refreshed.fingerprint != valid.fingerprint,
                   "refreshed compatible layers compose a new identity");

  const LinearLifecycleCounters before_hot = counters;
  bool hot_ok = true;
  std::size_t hot_allocations = 0U;
  {
    allocation_observer::Guard guard;
    for (std::size_t repetition = 0U; repetition < 100U; ++repetition) {
      hot_ok = static_cast<bool>(SymbolicPlan::compile(
                   symbolic_spec(), symbolic, &counters)) &&
               hot_ok;
      hot_ok = static_cast<bool>(numeric.refill(
                   symbolic, coefficient_revisions(29U), 803U, &counters)) &&
               hot_ok;
      hot_ok = static_cast<bool>(hierarchy.update(
                   symbolic, numeric, hierarchy_policy(),
                   HierarchyUpdate::retain, 804U, &counters)) &&
               hot_ok;
      hot_ok = static_cast<bool>(SolverWorkspace::bind(
                   fixture.requirements, fixture.vector_bundle,
                   fixture.scalar_buffer, workspace, &counters)) &&
               hot_ok;
      FieldView vector = workspace.vector(
          static_cast<std::uint8_t>(repetition %
                                    fixture.requirements.vector_slots),
          Int3{7, 4, 2});
      Span<double> scalar = workspace.scalars(2U, 4U);
      hot_ok = vector.base != nullptr && scalar.data != nullptr && hot_ok;
      vector.unchecked(Int3{0, 0, 0}, 0U) =
          static_cast<double>(repetition);
      scalar.data[0U] = static_cast<double>(repetition);
      hot_ok = compose_linear_identity(symbolic, numeric, hierarchy, workspace)
                       .fingerprint != 0U &&
               hot_ok;
    }
    hot_allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(hot_ok, "100 repeated lifecycle no-ops and borrows succeed");
  passed &= expect(hot_allocations == 0U,
                   "100 repeated lifecycle no-ops and borrows allocate zero heap bytes");
  passed &= expect(same(counters, before_hot),
                   "100 repeated no-ops leave every lifecycle event counter unchanged");
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  passed &= test_workspace_requirement_contracts();
  passed &= test_runtime_lifetime_views();
  passed &= test_symbolic_and_numeric_lifecycle();
  passed &= test_workspace_internal_overlap_is_rejected_atomically();
  passed &= test_arena_workspace_binding_and_borrowing();
  passed &= test_composed_identity_and_hot_noops();
  return passed ? 0 : 1;
}
