// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {
using namespace hundun::v04;

constexpr Int3 kCells{24, 16, 8};

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kCells;
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

OwnedFace make_face(CartesianAxis axis, StorageIdentity storage) {
  OwnedFace result;
  Int3 extents = kCells;
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
                 storage, 501U};
  return result;
}

struct Runtime {
  OwnedField vectors;
  MgWorkspaceRequirements requirements{};
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;

  bool create(const CartesianGeometryPlan& geometry, const MeshPatch& patch,
              const MgHierarchyPolicy& policy, StorageIdentity storage,
              RevisionDomainIdentity domain) {
    if (!make_mg_workspace_requirements(MPI_COMM_SELF, geometry, patch, policy,
                                        601U + storage, requirements)) {
      return false;
    }
    vectors = make_field(static_cast<FieldId>(30U + storage),
                         requirements.arena_shape, 0U, 1U, storage, domain);
    const std::array<HaloFieldSpec, 1U> halo_fields{{
        {vectors.view.field, 1U, 1U}}};
    if (!MgWorkspace::bind(requirements, vectors.view, workspace) ||
        !ReductionEngine::compile(MPI_COMM_SELF,
                                  ReductionMode::mpi_allreduce, 4U,
                                  reductions) ||
        !halo.reserve(MPI_COMM_SELF, requirements.levels[0U].patch,
                      {halo_fields.data(), halo_fields.size()})) {
      return false;
    }
    coarse_halos.resize(requirements.level_count - 1U);
    coarse_halo_pointers.resize(coarse_halos.size());
    for (std::size_t level = 1U; level < requirements.level_count; ++level) {
      if (!coarse_halos[level - 1U].reserve(
              MPI_COMM_SELF, requirements.levels[level].patch,
              {halo_fields.data(), halo_fields.size()})) {
        return false;
      }
      coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
    }
    return true;
  }

  MgRuntimeServices services() noexcept {
    return {&halo, &reductions, &workspace,
            {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
  }
};

double norm(ConstFieldView field) noexcept {
  double sum = 0.0;
  for (std::int32_t z = 0; z < field.interior.z; ++z) {
    for (std::int32_t y = 0; y < field.interior.y; ++y) {
      for (std::int32_t x = 0; x < field.interior.x; ++x) {
        const double value = field.unchecked({x, y, z}, 0U);
        sum += value * value;
      }
    }
  }
  return std::sqrt(sum);
}

double difference(ConstFieldView left, ConstFieldView right) noexcept {
  double sum = 0.0;
  for (std::int32_t z = 0; z < left.interior.z; ++z) {
    for (std::int32_t y = 0; y < left.interior.y; ++y) {
      for (std::int32_t x = 0; x < left.interior.x; ++x) {
        const double value = left.unchecked({x, y, z}, 0U) -
                             right.unchecked({x, y, z}, 0U);
        sum += value * value;
      }
    }
  }
  return std::sqrt(sum);
}

bool test_coarse_policy_changes_real_correction() {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh_spec(), GeometryBudget{}, geometry, patch)),
      "coarse-use geometry compiles");
  if (!passed) {
    return false;
  }
  constexpr RevisionDomainIdentity domain = 501U;
  OwnedField diagonal = make_field(1U, patch.cells, 0U, 1U, 510U, domain);
  OwnedField residual = make_field(2U, patch.cells, 0U, 1U, 511U, domain);
  OwnedField shallow_correction =
      make_field(3U, patch.cells, 0U, 1U, 512U, domain);
  OwnedField deep_correction =
      make_field(4U, patch.cells, 0U, 1U, 513U, domain);
  for (std::int32_t k = 0; k < patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < patch.cells.x; ++i) {
        diagonal.view.unchecked({i, j, k}, 0U) = 6.0;
        // A smooth, low-frequency residual is deliberately chosen: local
        // smoothing treats both plans identically while coarse correction
        // strength depends strongly on actual transfer/coarse solve use.
        residual.view.unchecked({i, j, k}, 0U) =
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(i) + 0.5) / kCells.x) *
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(j) + 0.5) / kCells.y) *
            std::sin(3.14159265358979323846 *
                     (static_cast<double>(k) + 0.5) / kCells.z);
      }
    }
  }
  OwnedFace x = make_face(CartesianAxis::x, 520U);
  OwnedFace y = make_face(CartesianAxis::y, 521U);
  OwnedFace z = make_face(CartesianAxis::z, 522U);
  const MgCoefficientViews coefficients{as_const(diagonal.view), x.view,
                                        y.view, z.view};

  NativeCartesianMgSpec shallow_spec;
  shallow_spec.communicator = MPI_COMM_SELF;
  shallow_spec.geometry = &geometry;
  shallow_spec.patch = patch;
  shallow_spec.identity = {11U, 12U, 13U, 14U, 15U};
  shallow_spec.coefficients = {21U, 22U, 0.0};
  shallow_spec.policy.pre_sweeps = 1U;
  shallow_spec.policy.post_sweeps = 1U;
  shallow_spec.policy.maximum_levels = 2U;
  shallow_spec.policy.coarse_sweeps = 1U;

  NativeCartesianMgSpec deep_spec = shallow_spec;
  deep_spec.policy.maximum_levels = 16U;
  deep_spec.policy.coarse_sweeps = 48U;

  Runtime shallow_runtime;
  Runtime deep_runtime;
  passed &= expect(
      shallow_runtime.create(geometry, patch, shallow_spec.policy, 530U,
                             domain) &&
          deep_runtime.create(geometry, patch, deep_spec.policy, 540U, domain),
      "independent runtime services compile");
  if (!passed) {
    return false;
  }

  NativeCartesianMgPlan shallow;
  NativeCartesianMgPlan deep;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          shallow_spec,
          shallow_runtime.services(),
          coefficients, shallow)) &&
          static_cast<bool>(NativeCartesianMgPlan::compile(
              deep_spec,
              deep_runtime.services(),
              coefficients, deep)) &&
          shallow.level_count() == 2U && deep.level_count() > 2U,
      "shallow and deep hierarchies differ only in coarse policy");
  if (!passed) {
    return false;
  }

  passed &= expect(
      static_cast<bool>(shallow.apply(as_const(residual.view),
                                      shallow_correction.view, 0U)) &&
          static_cast<bool>(deep.apply(as_const(residual.view),
                                       deep_correction.view, 0U)),
      "both coarse policies produce finite corrections");
  const double correction_difference =
      difference(as_const(shallow_correction.view),
                 as_const(deep_correction.view));
  const double deep_norm = norm(as_const(deep_correction.view));
  passed &= expect(
      std::isfinite(correction_difference) &&
          correction_difference > 1.0e-8 * std::max(1.0, deep_norm) &&
          deep.last_cycle_final_residual() <
              shallow.last_cycle_final_residual(),
      "coarse depth/sweeps materially change correction and improve true cycle residual");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_coarse_policy_changes_real_correction();
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
