// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include "solver_mg_detail.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <vector>

namespace allocation_observer {
std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};
void* allocate(std::size_t bytes) {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
  void* const result = std::malloc(bytes == 0U ? 1U : bytes);
  if (result == nullptr) {
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
};
}  // namespace allocation_observer

void* operator new(std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new[](std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  void* result = nullptr;
  if (allocation_observer::enabled.load(std::memory_order_relaxed)) {
    allocation_observer::count.fetch_add(1U, std::memory_order_relaxed);
  }
  if (posix_memalign(&result, static_cast<std::size_t>(alignment),
                     bytes == 0U ? static_cast<std::size_t>(alignment)
                                 : bytes) != 0 ||
      result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
  return ::operator new(bytes, alignment);
}
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(bytes);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t bytes, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(bytes);
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

constexpr Int3 kCells{16, 12, 8};

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
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

CartesianMeshSpec uniform_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kCells;
  mesh.minimum_spacing = {1.0e-8, 1.0e-8, 1.0e-8};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

CartesianMeshSpec stretched_x_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::tensor_stretched;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {4.0, 1.0, 1.0};
  mesh.has_exact_cells = false;
  mesh.has_base_spacing = true;
  mesh.base_spacing = {0.25, 1.0 / 12.0, 0.125};
  mesh.minimum_spacing = {0.025, 1.0 / 12.0, 0.125};
  mesh.max_growth_ratio = 1.25;
  mesh.focus_regions.push_back(
      {{1.0, 0.0, 0.0}, {3.0, 1.0, 1.0}, {0.05, 1.0 / 12.0, 0.125}});
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField field(Int3 cells, std::uint8_t components, std::uint8_t ghosts,
                 FieldId id, StorageIdentity storage) {
  OwnedField result;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz * components, 0.0);
  result.view.base = result.storage.data() + ghosts + ghosts * nx +
                     ghosts * nx * ny;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = storage;
  result.view.revision_domain = 900U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  ConstFaceFieldView view{};
};

OwnedFace face(Int3 extents, CartesianAxis axis, double value) {
  OwnedFace result;
  const auto nx = static_cast<std::size_t>(extents.x);
  const auto ny = static_cast<std::size_t>(extents.y);
  const auto nz = static_cast<std::size_t>(extents.z);
  result.storage.assign(nx * ny * nz, value);
  result.view.base = result.storage.data();
  result.view.extents = extents;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.axis = axis;
  result.view.storage_identity = 800U + static_cast<std::uint8_t>(axis);
  result.view.revision_domain = 901U;
  return result;
}

struct Fixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  OwnedField diagonal;
  OwnedFace x;
  OwnedFace y;
  OwnedFace z;
  OwnedField vectors;
  MgWorkspaceRequirements workspace_requirements{};
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  NativeCartesianMgSpec spec{};
  MgRuntimeServices services{};

  bool create(const CartesianMeshSpec& mesh) {
    if (!CartesianGeometryCompiler::compile(
            MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch)) {
      return false;
    }
    const Int3 cells = patch.cells;
    diagonal = field(cells, 1U, 0U, 1U, 700U);
    x = face({cells.x + 1, cells.y, cells.z}, CartesianAxis::x, 1.0);
    y = face({cells.x, cells.y + 1, cells.z}, CartesianAxis::y, 1.0);
    z = face({cells.x, cells.y, cells.z + 1}, CartesianAxis::z, 1.0);
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          diagonal.view.unchecked({i, j, k}, 0U) = 6.0;
        }
      }
    }
    spec.communicator = MPI_COMM_SELF;
    spec.geometry = &geometry;
    spec.patch = patch;
    spec.boundaries = {};
    spec.policy.anisotropy_threshold = 4.0;
    spec.policy.coefficient_change_rebuild_ratio = 0.25;
    spec.policy.pre_sweeps = 1U;
    spec.policy.post_sweeps = 1U;
    spec.policy.maximum_levels = 16U;
    spec.policy.coarse_sweeps = 20U;
    spec.policy.minimum_coarse_extent = 3U;
    spec.policy.line_relaxation_maximum_extent = 4096U;
    spec.identity = {101U, 102U, 103U, 104U, 105U};
    spec.coefficients = {1U, 201U, 0.0};
    if (!make_mg_workspace_requirements(MPI_COMM_SELF, geometry, patch,
                                        spec.policy, 71U,
                                        workspace_requirements)) {
      return false;
    }
    vectors = field(workspace_requirements.arena_shape, 1U, 0U, 10U, 500U);
    if (!MgWorkspace::bind(workspace_requirements, vectors.view, workspace) ||
        !ReductionEngine::compile(MPI_COMM_SELF,
                                  ReductionMode::mpi_allreduce, 4U,
                                  reductions)) {
      return false;
    }
    const std::array<HaloFieldSpec, 1U> fields{{{10U, 1U, 1U}}};
    if (!halo.reserve(MPI_COMM_SELF,
                      workspace_requirements.levels[0U].patch,
                      {fields.data(), fields.size()})) {
      return false;
    }
    coarse_halos.resize(workspace_requirements.level_count - 1U);
    coarse_halo_pointers.resize(coarse_halos.size());
    for (std::size_t level = 1U;
         level < workspace_requirements.level_count; ++level) {
      if (!coarse_halos[level - 1U].reserve(
              MPI_COMM_SELF, workspace_requirements.levels[level].patch,
              {fields.data(), fields.size()})) {
        return false;
      }
      coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
    }
    services = {&halo, &reductions, &workspace,
                {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
    return true;
  }

  MgCoefficientViews coefficients() const noexcept {
    return {as_const(diagonal.view), x.view, y.view, z.view};
  }
};

bool test_coarsening_selection() {
  Fixture isotropic;
  bool passed = expect(isotropic.create(uniform_mesh()),
                       "isotropic fixture compiles");
  NativeCartesianMgPlan full;
  MgPlanCounters full_counters{};
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          isotropic.spec, isotropic.services, isotropic.coefficients(), full,
          &full_counters)),
      "isotropic native hierarchy compiles");
  passed &= expect(full.finest_coarsening() == CoarseningKind::full_xyz &&
                       full.line_axis_mask() == 0U && full.level_count() > 1U,
                   "near-isotropic metrics select full XYZ coarsening");
  MgLevelView level{};
  passed &= expect(static_cast<bool>(full.level(0U, level)) &&
                       level.coarsening == CoarseningKind::full_xyz &&
                       level.line_axis_mask == 0U,
                   "published finest-level strategy is full coarsening");

  Fixture stretched;
  passed &= expect(stretched.create(stretched_x_mesh()),
                   "x-stretched fixture compiles");
  NativeCartesianMgPlan semi;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          stretched.spec, stretched.services, stretched.coefficients(), semi)),
      "stretched native hierarchy compiles");
  constexpr std::uint8_t kXMask = 1U;
  passed &= expect(semi.finest_coarsening() == CoarseningKind::semi_yz &&
                       semi.line_axis_mask() == kXMask,
                   "strong x stretching stays uncoarsened and selects x lines");
  return passed;
}

bool test_reuse_policy_and_hot_lifetime() {
  Fixture fixture;
  if (!expect(fixture.create(uniform_mesh()), "reuse fixture compiles")) {
    return false;
  }
  NativeCartesianMgPlan plan;
  MgPlanCounters external{};
  bool passed = expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          fixture.spec, fixture.services, fixture.coefficients(), plan,
          &external)),
      "baseline hierarchy compiles");
  const auto hierarchy_address = plan.hierarchy_storage_address();
  const auto workspace_address = plan.workspace_storage_address();
  const auto symbolic = plan.symbolic_fingerprint();
  const auto numeric = plan.numeric_fingerprint();
  const auto hierarchy = plan.hierarchy_fingerprint();
  const auto generation = plan.generation();
  const MgPlanCounters baseline = plan.counters();
  passed &= expect(hierarchy_address != 0U && workspace_address != 0U &&
                       symbolic != 0U && hierarchy != 0U && generation != 0U &&
                       baseline.symbolic_builds == 1U &&
                       baseline.hierarchy_rebuilds == 1U,
                   "cold compile publishes one symbolic hierarchy and storage");

  passed &= expect(
      static_cast<bool>(plan.update_coefficients(
          fixture.spec.identity, fixture.spec.coefficients,
          fixture.coefficients(), &external)) &&
          same(plan.counters(), baseline) &&
          plan.hierarchy_storage_address() == hierarchy_address &&
          plan.workspace_storage_address() == workspace_address,
      "unchanged coefficient identity performs no lifecycle work");

  const MgCoefficientIdentity below{2U, 202U, 0.24};
  LinearIdentity below_identity = fixture.spec.identity;
  below_identity.numeric = 112U;
  below_identity.hierarchy = 113U;
  below_identity.fingerprint = 115U;
  const Status below_status = plan.update_coefficients(
      below_identity, below, fixture.coefficients(), &external);
  const bool below_ok =
      static_cast<bool>(below_status) &&
          plan.counters().numeric_refreshes == baseline.numeric_refreshes + 1U &&
          plan.counters().hierarchy_rebuilds == baseline.hierarchy_rebuilds &&
          plan.symbolic_fingerprint() == symbolic &&
          plan.numeric_fingerprint() != numeric &&
          plan.hierarchy_fingerprint() == hierarchy &&
          plan.generation() == generation &&
          plan.hierarchy_storage_address() == hierarchy_address &&
          plan.workspace_storage_address() == workspace_address;
  passed &= expect(
      below_ok,
      "below-threshold refresh preserves hierarchy identity and generation");

  const auto below_numeric = plan.numeric_fingerprint();
  const auto below_hierarchy = plan.hierarchy_fingerprint();
  const auto below_generation = plan.generation();
  const MgCoefficientIdentity above{3U, 203U, 0.26};
  LinearIdentity above_identity = below_identity;
  above_identity.numeric = 122U;
  above_identity.hierarchy = 123U;
  above_identity.fingerprint = 125U;
  passed &= expect(
      static_cast<bool>(plan.update_coefficients(
          above_identity, above, fixture.coefficients(), &external)) &&
          plan.counters().numeric_refreshes == baseline.numeric_refreshes + 2U &&
          plan.counters().hierarchy_rebuilds == baseline.hierarchy_rebuilds + 1U &&
          plan.symbolic_fingerprint() == symbolic &&
          plan.numeric_fingerprint() != below_numeric &&
          plan.hierarchy_fingerprint() != below_hierarchy &&
          plan.generation() > below_generation &&
          plan.hierarchy_storage_address() == hierarchy_address &&
          plan.workspace_storage_address() == workspace_address,
      "above-threshold change rebuilds exactly one hierarchy");

  const MgPlanCounters before_hot = plan.counters();
  const auto hot_hierarchy_address = plan.hierarchy_storage_address();
  const auto hot_workspace_address = plan.workspace_storage_address();
  {
    allocation_observer::Guard guard;
    for (std::uint32_t repeat = 0U; repeat < 100U; ++repeat) {
      if (!plan.update_coefficients(above_identity, above,
                                    fixture.coefficients(), nullptr)) {
        passed = false;
      }
    }
  }
  passed &= expect(allocation_observer::count.load(std::memory_order_relaxed) ==
                           0U &&
                       same(plan.counters(), before_hot) &&
                       plan.hierarchy_storage_address() ==
                           hot_hierarchy_address &&
                       plan.workspace_storage_address() == hot_workspace_address,
                   "100 unchanged hot updates allocate nothing and preserve addresses");
  return passed;
}

bool test_level_halo_pointer_table_is_plan_owned() {
  Fixture fixture;
  if (!expect(fixture.create(uniform_mesh()),
              "pointer-table fixture compiles")) {
    return false;
  }
  NativeCartesianMgPlan plan;
  bool passed = expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          fixture.spec, fixture.services, fixture.coefficients(), plan)),
      "plan compiles before caller mutates its pointer table");
  if (!passed || fixture.coarse_halo_pointers.empty()) {
    return passed;
  }

  for (HaloEngine*& pointer : fixture.coarse_halo_pointers) {
    pointer = &fixture.halo;
  }
  OwnedField residual = field(kCells, 1U, 0U, 20U, 600U);
  OwnedField correction = field(kCells, 1U, 0U, 21U, 601U);
  for (double& value : residual.storage) {
    value = 1.0;
  }
  passed &= expect(
      static_cast<bool>(plan.apply(as_const(residual.view), correction.view,
                                   0U)),
      "caller pointer-table mutation cannot alter compiled plan services");
  return passed;
}

bool test_borrowed_service_rebinding_is_rejected() {
  bool passed = true;
  {
    Fixture fixture;
    NativeCartesianMgPlan plan;
    passed &= expect(fixture.create(uniform_mesh()) &&
                         static_cast<bool>(NativeCartesianMgPlan::compile(
                             fixture.spec, fixture.services,
                             fixture.coefficients(), plan)),
                     "workspace-rebind plan compiles");
    const PlanFingerprint numeric = plan.numeric_fingerprint();
    OwnedField alternate = field(fixture.workspace_requirements.arena_shape,
                                 1U, 0U, 10U, 650U);
    MgWorkspace replacement;
    passed &= expect(static_cast<bool>(MgWorkspace::bind(
                         fixture.workspace_requirements, alternate.view,
                         replacement)),
                     "replacement workspace binds");
    fixture.workspace = std::move(replacement);
    LinearIdentity next = fixture.spec.identity;
    next.numeric = 412U;
    next.hierarchy = 413U;
    next.fingerprint = 415U;
    const Status rejected = plan.update_coefficients(
        next, MgCoefficientIdentity{2U, 502U, 0.10},
        fixture.coefficients(), nullptr);
    passed &= expect(rejected.code == StatusCode::invalid_plan &&
                         plan.numeric_fingerprint() == numeric,
                     "rebound workspace is rejected before numeric mutation");
  }
  {
    Fixture fixture;
    NativeCartesianMgPlan plan;
    passed &= expect(fixture.create(uniform_mesh()) &&
                         static_cast<bool>(NativeCartesianMgPlan::compile(
                             fixture.spec, fixture.services,
                             fixture.coefficients(), plan)),
                     "reduction-rebind plan compiles");
    ReductionEngine replacement;
    passed &= expect(static_cast<bool>(ReductionEngine::compile(
                         MPI_COMM_SELF, ReductionMode::mpi_allreduce, 4U,
                         replacement)),
                     "replacement reductions compile");
    fixture.reductions = std::move(replacement);
    LinearIdentity next = fixture.spec.identity;
    next.numeric = 512U;
    next.hierarchy = 513U;
    next.fingerprint = 515U;
    const Status rejected = plan.update_coefficients(
        next, MgCoefficientIdentity{2U, 602U, 0.10},
        fixture.coefficients(), nullptr);
    passed &= expect(rejected.code == StatusCode::invalid_plan,
                     "rebound reduction engine is rejected pre-use");
  }
  {
    Fixture fixture;
    NativeCartesianMgPlan plan;
    passed &= expect(fixture.create(uniform_mesh()) &&
                         static_cast<bool>(NativeCartesianMgPlan::compile(
                             fixture.spec, fixture.services,
                             fixture.coefficients(), plan)),
                     "halo-move plan compiles");
    HaloEngine moved = std::move(fixture.halo);
    LinearIdentity next = fixture.spec.identity;
    next.numeric = 612U;
    next.hierarchy = 613U;
    next.fingerprint = 615U;
    const Status rejected = plan.update_coefficients(
        next, MgCoefficientIdentity{2U, 702U, 0.10},
        fixture.coefficients(), nullptr);
    passed &= expect(rejected.code == StatusCode::invalid_plan && moved.ready(),
                     "moved finest halo service is rejected pre-use");
  }
  {
    Fixture fixture;
    NativeCartesianMgPlan plan;
    passed &= expect(fixture.create(uniform_mesh()) &&
                         static_cast<bool>(NativeCartesianMgPlan::compile(
                             fixture.spec, fixture.services,
                             fixture.coefficients(), plan)) &&
                         !fixture.coarse_halos.empty(),
                     "coarse-halo-move plan compiles");
    if (!fixture.coarse_halos.empty()) {
      HaloEngine moved = std::move(fixture.coarse_halos[0U]);
      LinearIdentity next = fixture.spec.identity;
      next.numeric = 712U;
      next.hierarchy = 713U;
      next.fingerprint = 715U;
      const Status rejected = plan.update_coefficients(
          next, MgCoefficientIdentity{2U, 802U, 0.10},
          fixture.coefficients(), nullptr);
      passed &= expect(rejected.code == StatusCode::invalid_plan && moved.ready(),
                       "moved coarse halo service is rejected pre-use");
    }
  }
  return passed;
}

bool test_rebuild_generation_overflow_is_atomic() {
  Fixture fixture;
  if (!expect(fixture.create(uniform_mesh()), "overflow fixture compiles")) {
    return false;
  }
  NativeCartesianMgPlan plan;
  bool passed = expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          fixture.spec, fixture.services, fixture.coefficients(), plan)),
      "plan compiles before injected generation overflow");
  if (!passed) {
    return false;
  }
  const MgPlanCounters saturated = plan.counters();
  detail::set_mg_runtime_counters_for_test(
      plan, saturated, std::numeric_limits<RevisionToken>::max());
  const PlanFingerprint numeric_before = plan.numeric_fingerprint();
  const PlanFingerprint hierarchy_before = plan.hierarchy_fingerprint();
  const RevisionToken generation_before = plan.generation();
  LinearIdentity next = fixture.spec.identity;
  next.numeric = 312U;
  next.hierarchy = 313U;
  next.fingerprint = 315U;
  const Status rejected = plan.update_coefficients(
      next, MgCoefficientIdentity{2U, 402U, 0.26},
      fixture.coefficients(), nullptr);
  passed &= expect(
      rejected.code == StatusCode::invalid_plan &&
          plan.numeric_fingerprint() == numeric_before &&
          plan.hierarchy_fingerprint() == hierarchy_before &&
          plan.generation() == generation_before &&
          same(plan.counters(), saturated),
      "saturated rebuild counter/generation rejects before state publication");
  return passed;
}

bool test_numeric_counter_overflow_is_atomic() {
  Fixture fixture;
  if (!expect(fixture.create(uniform_mesh()),
              "counter-overflow fixture compiles")) {
    return false;
  }
  NativeCartesianMgPlan plan;
  bool passed = expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          fixture.spec, fixture.services, fixture.coefficients(), plan)),
      "plan compiles before injected counter overflow");
  if (!passed) {
    return false;
  }
  MgPlanCounters saturated = plan.counters();
  saturated.numeric_refreshes =
      std::numeric_limits<std::uint64_t>::max();
  detail::set_mg_runtime_counters_for_test(plan, saturated, plan.generation());
  const PlanFingerprint numeric_before = plan.numeric_fingerprint();
  LinearIdentity next = fixture.spec.identity;
  next.numeric = 812U;
  next.hierarchy = 813U;
  next.fingerprint = 815U;
  const Status rejected = plan.update_coefficients(
      next, MgCoefficientIdentity{2U, 902U, 0.10},
      fixture.coefficients(), nullptr);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       plan.numeric_fingerprint() == numeric_before &&
                       same(plan.counters(), saturated),
                   "saturated numeric counter rejects atomically");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_coarsening_selection();
  passed &= test_reuse_policy_and_hot_lifetime();
  passed &= test_level_halo_pointer_table_is_plan_owned();
  passed &= test_borrowed_service_rebinding_is_rejected();
  passed &= test_rebuild_generation_overflow_is_atomic();
  passed &= test_numeric_counter_overflow_is_atomic();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 native Cartesian MG reuse tests passed\n";
  return 0;
}
