// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <array>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
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

bool all_true(bool local) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {16, 8, 4};
  mesh.minimum_spacing = {1.0e-10, 1.0e-10, 1.0e-10};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 30U};
  return mesh;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t ghosts,
                      std::uint8_t components, StorageIdentity storage,
                      RevisionDomainIdentity domain) {
  OwnedField result;
  const std::size_t stride_y =
      static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t stride_z =
      stride_y * static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t component_stride =
      stride_z * static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(component_stride * components, 0.0);
  result.view.base = result.storage.data() + ghosts +
                     static_cast<std::size_t>(ghosts) *
                         (stride_y + stride_z);
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = stride_y;
  result.view.stride_z = stride_z;
  result.view.component_stride = component_stride;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = storage;
  result.view.revision_domain = domain;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  ConstFaceFieldView view{};
};

OwnedFace make_face(CartesianAxis axis, Int3 cells,
                    StorageIdentity storage,
                    RevisionDomainIdentity domain) {
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
                 storage, domain};
  return result;
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

struct Fixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  OwnedField diagonal;
  OwnedFace x;
  OwnedFace y;
  OwnedFace z;
  OwnedField residual;
  OwnedField correction;
  OwnedField vectors;
  MgWorkspaceRequirements requirements{};
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  NativeCartesianMgSpec spec{};
  NativeCartesianMgPlan plan;

  bool create() {
    if (!CartesianGeometryCompiler::compile(
            MPI_COMM_WORLD, mesh_spec(), GeometryBudget{}, geometry, patch)) {
      return false;
    }
    constexpr RevisionDomainIdentity domain = 701U;
    diagonal = make_field(1U, patch.cells, 0U, 1U, 710U, domain);
    x = make_face(CartesianAxis::x, patch.cells, 711U, domain);
    y = make_face(CartesianAxis::y, patch.cells, 712U, domain);
    z = make_face(CartesianAxis::z, patch.cells, 713U, domain);
    residual = make_field(2U, patch.cells, 0U, 1U, 714U, domain);
    correction = make_field(3U, patch.cells, 0U, 1U, 715U, domain);
    for (std::int32_t k = 0; k < patch.cells.z; ++k) {
      for (std::int32_t j = 0; j < patch.cells.y; ++j) {
        for (std::int32_t i = 0; i < patch.cells.x; ++i) {
          diagonal.view.unchecked({i, j, k}, 0U) = 6.0;
          residual.view.unchecked({i, j, k}, 0U) = 1.0;
          correction.view.unchecked({i, j, k}, 0U) = -19.0;
        }
      }
    }
    spec.communicator = MPI_COMM_WORLD;
    spec.geometry = &geometry;
    spec.patch = patch;
    spec.identity = {11U, 12U, 13U, 14U, 15U};
    spec.coefficients = {21U, 22U, 0.0};
    if (!make_mg_workspace_requirements(MPI_COMM_WORLD, geometry, patch,
                                        spec.policy, 801U, requirements)) {
      return false;
    }
    vectors = make_field(30U, requirements.arena_shape, 0U, 1U, 720U,
                         domain);
    if (!MgWorkspace::bind(requirements, vectors.view, workspace) ||
        !ReductionEngine::compile(MPI_COMM_WORLD,
                                  ReductionMode::mpi_allreduce, 4U,
                                  reductions)) {
      return false;
    }
    const std::array<HaloFieldSpec, 1U> fields{{{30U, 1U, 1U}}};
    if (!halo.reserve(MPI_COMM_WORLD, requirements.levels[0U].patch,
                      {fields.data(), fields.size()})) {
      return false;
    }
    coarse_halos.resize(requirements.level_count - 1U);
    coarse_halo_pointers.resize(coarse_halos.size());
    for (std::size_t level = 1U; level < requirements.level_count; ++level) {
      if (!coarse_halos[level - 1U].reserve(
              MPI_COMM_WORLD, requirements.levels[level].patch,
              {fields.data(), fields.size()})) {
        return false;
      }
      coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
    }
    return static_cast<bool>(NativeCartesianMgPlan::compile(
        spec,
        MgRuntimeServices{
            &halo, &reductions, &workspace,
            {coarse_halo_pointers.data(), coarse_halo_pointers.size()}},
        coefficients(), plan));
  }

  MgCoefficientViews coefficients() const noexcept {
    return {as_const(diagonal.view), x.view, y.view, z.view};
  }
};

bool test_divergent_update_contract(int rank, int size) {
  Fixture fixture;
  bool passed = expect(fixture.create(), rank, "baseline MG plan compiles");
  if (!all_true(passed)) {
    return false;
  }
  const PlanFingerprint numeric = fixture.plan.numeric_fingerprint();
  const PlanFingerprint hierarchy = fixture.plan.hierarchy_fingerprint();
  const RevisionToken generation = fixture.plan.generation();
  const auto hierarchy_address = fixture.plan.hierarchy_storage_address();
  const auto workspace_address = fixture.plan.workspace_storage_address();
  const MgPlanCounters before = fixture.plan.counters();
  const std::vector<double> output_before = fixture.correction.storage;

  // Rank zero selects the unchanged fast path; every other rank selects the
  // changed path. A correct implementation must compare the collective
  // identity before any rank returns or mutates hierarchy state.
  const MgCoefficientIdentity identity =
      rank == 0 ? fixture.spec.coefficients
                : MgCoefficientIdentity{23U, 24U, 0.10};
  const Status status = fixture.plan.update_coefficients(
      fixture.spec.identity, identity, fixture.coefficients(), nullptr);

  std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  const int expected_lowest_mismatch = size > 1 ? 1 : -1;
  passed &= expect(
      status.code == StatusCode::invalid_plan && minimum == maximum &&
          fixture.plan.lowest_failing_rank() == expected_lowest_mismatch &&
          fixture.plan.numeric_fingerprint() == numeric &&
          fixture.plan.hierarchy_fingerprint() == hierarchy &&
          fixture.plan.generation() == generation &&
          fixture.plan.hierarchy_storage_address() == hierarchy_address &&
          fixture.plan.workspace_storage_address() == workspace_address &&
          same(fixture.plan.counters(), before) &&
          fixture.correction.storage == output_before,
      rank,
      "unchanged/changed rank divergence rejects collectively and atomically");
  return all_true(passed);
}

bool test_reversed_runtime_communicators(int rank, int size) {
  Fixture fixture;
  bool passed = expect(fixture.create(), rank,
                       "baseline plan exists before communicator RED");
  if (!all_true(passed)) return false;

  MPI_Comm reversed = MPI_COMM_NULL;
  const int split_status = MPI_Comm_split(
      MPI_COMM_WORLD, 0, size - 1 - rank, &reversed);
  passed &= expect(split_status == MPI_SUCCESS && reversed != MPI_COMM_NULL,
                   rank, "reversed communicator constructs");
  if (!all_true(passed)) {
    if (reversed != MPI_COMM_NULL) MPI_Comm_free(&reversed);
    return false;
  }

  ReductionEngine reversed_reductions;
  bool setup = static_cast<bool>(ReductionEngine::compile(
      reversed, ReductionMode::mpi_allreduce, 4U, reversed_reductions));
  passed &= expect(setup, rank,
                   "reversed reduction service independently compiles");
  if (!all_true(passed)) {
    MPI_Comm_free(&reversed);
    return false;
  }

  const PlanFingerprint numeric = fixture.plan.numeric_fingerprint();
  const PlanFingerprint hierarchy = fixture.plan.hierarchy_fingerprint();
  const RevisionToken generation = fixture.plan.generation();
  const auto hierarchy_address = fixture.plan.hierarchy_storage_address();
  const auto workspace_address = fixture.plan.workspace_storage_address();
  const MgPlanCounters before = fixture.plan.counters();
  const std::vector<double> output_before = fixture.correction.storage;
  const Status status = NativeCartesianMgPlan::compile(
      fixture.spec,
      MgRuntimeServices{
          &fixture.halo, &reversed_reductions, &fixture.workspace,
          {fixture.coarse_halo_pointers.data(),
           fixture.coarse_halo_pointers.size()}},
      fixture.coefficients(), fixture.plan);
  std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(
      status.code == StatusCode::invalid_plan && minimum == maximum &&
          fixture.plan.numeric_fingerprint() == numeric &&
          fixture.plan.hierarchy_fingerprint() == hierarchy &&
          fixture.plan.generation() == generation &&
          fixture.plan.hierarchy_storage_address() == hierarchy_address &&
          fixture.plan.workspace_storage_address() == workspace_address &&
          same(fixture.plan.counters(), before) &&
          fixture.correction.storage == output_before,
      rank,
      "reversed reductions/halos reject collectively before mutation");
  MPI_Comm_free(&reversed);
  return all_true(passed);
}

bool test_rank_local_null_reduction_preflight(int rank, int size) {
  Fixture fixture;
  bool passed = expect(fixture.create(), rank,
                       "baseline plan exists before null-reduction RED");
  if (!all_true(passed)) return false;

  const PlanFingerprint symbolic = fixture.plan.symbolic_fingerprint();
  const PlanFingerprint numeric = fixture.plan.numeric_fingerprint();
  const PlanFingerprint hierarchy = fixture.plan.hierarchy_fingerprint();
  const RevisionToken generation = fixture.plan.generation();
  const auto hierarchy_address = fixture.plan.hierarchy_storage_address();
  const auto workspace_address = fixture.plan.workspace_storage_address();
  const MgPlanCounters plan_before = fixture.plan.counters();
  MgPlanCounters external{31U, 32U, 33U, 34U};
  const MgPlanCounters external_before = external;

  MgRuntimeServices runtime{
      &fixture.halo,
      rank == size - 1 ? nullptr : &fixture.reductions,
      &fixture.workspace,
      {fixture.coarse_halo_pointers.data(),
       fixture.coarse_halo_pointers.size()}};
  const Status status = NativeCartesianMgPlan::compile(
      fixture.spec, runtime, fixture.coefficients(), fixture.plan, &external);
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);

  passed &= expect(
      status.code == StatusCode::invalid_plan && minimum == maximum &&
          fixture.plan.lowest_failing_rank() == size - 1 &&
          fixture.plan.symbolic_fingerprint() == symbolic &&
          fixture.plan.numeric_fingerprint() == numeric &&
          fixture.plan.hierarchy_fingerprint() == hierarchy &&
          fixture.plan.generation() == generation &&
          fixture.plan.hierarchy_storage_address() == hierarchy_address &&
          fixture.plan.workspace_storage_address() == workspace_address &&
          same(fixture.plan.counters(), plan_before) &&
          same(external, external_before),
      rank,
      "rank-local null reduction rejects collectively before publication");
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
                       "update contract RED runs at 2 or 4 ranks");
  passed &= test_divergent_update_contract(rank, size);
  passed &= test_reversed_runtime_communicators(rank, size);
  passed &= test_rank_local_null_reduction_preflight(rank, size);
  passed = all_true(passed);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
