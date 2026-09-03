// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string_view>
#include <vector>

namespace {
using namespace hundun::v04;

constexpr Int3 kCells{12, 8, 6};

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

OwnedFace make_face(CartesianAxis axis, Int3 cells,
                    StorageIdentity storage) {
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
                 storage, 80U};
  return result;
}

std::uint64_t checksum(const std::vector<double>& values) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  for (const double value : values) {
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    hash ^= bits;
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

bool test_mg_never_borrows_krylov_slots() {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh_spec(), GeometryBudget{}, geometry, patch)),
      "isolation geometry compiles");
  if (!passed) {
    return false;
  }

  LinearWorkspaceRequirements requirements{};
  passed &= expect(static_cast<bool>(make_linear_workspace_requirements(
                       LinearAlgorithm::fgmres, patch.cells, 1U, 8U,
                       ReductionMode::mpi_allreduce, 701U, requirements)),
                   "independent Krylov workspace requirements compile");
  OwnedField krylov_vectors =
      make_field(30U, patch.cells, 1U, requirements.vector_slots, 100U, 200U);
  OwnedField krylov_scalars = make_field(
      31U,
      {static_cast<std::int32_t>(requirements.scalar_doubles), 1, 1}, 0U,
      1U, 100U, 200U);
  for (std::size_t index = 0U; index < krylov_vectors.storage.size(); ++index) {
    krylov_vectors.storage[index] =
        1000.0 + static_cast<double>(index) * 0.125;
  }
  for (std::size_t index = 0U; index < krylov_scalars.storage.size(); ++index) {
    krylov_scalars.storage[index] =
        -2000.0 - static_cast<double>(index) * 0.25;
  }
  const std::uint64_t vector_before = checksum(krylov_vectors.storage);
  const std::uint64_t scalar_before = checksum(krylov_scalars.storage);
  SolverWorkspace krylov_workspace;
  passed &= expect(static_cast<bool>(SolverWorkspace::bind(
                       requirements, krylov_vectors.view, krylov_scalars.view,
                       krylov_workspace)),
                   "Krylov workspace binds before MG compile");

  OwnedField diagonal = make_field(1U, patch.cells, 0U, 1U, 110U, 210U);
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        diagonal.view.unchecked({x, y, z}, 0U) = 6.0;
      }
    }
  }
  OwnedFace x = make_face(CartesianAxis::x, patch.cells, 120U);
  OwnedFace y = make_face(CartesianAxis::y, patch.cells, 121U);
  OwnedFace z = make_face(CartesianAxis::z, patch.cells, 122U);
  OwnedField residual = make_field(2U, patch.cells, 0U, 1U, 111U, 211U);
  OwnedField correction = make_field(3U, patch.cells, 0U, 1U, 112U, 212U);
  for (std::int32_t k = 0; k < patch.cells.z; ++k) {
    for (std::int32_t j = 0; j < patch.cells.y; ++j) {
      for (std::int32_t i = 0; i < patch.cells.x; ++i) {
        residual.view.unchecked({i, j, k}, 0U) =
            1.0 + 0.01 * static_cast<double>(i + 2 * j + 3 * k);
        correction.view.unchecked({i, j, k}, 0U) = -9.0;
      }
    }
  }

  ReductionEngine reductions;
  HaloEngine halo;
  NativeCartesianMgSpec spec;
  spec.communicator = MPI_COMM_SELF;
  spec.geometry = &geometry;
  spec.patch = patch;
  spec.identity = {11U, 12U, 13U, krylov_workspace.fingerprint(), 15U};
  spec.coefficients = {21U, 22U, 0.0};
  MgWorkspaceRequirements mg_requirements{};
  passed &= expect(static_cast<bool>(make_mg_workspace_requirements(
                       MPI_COMM_SELF, geometry, patch, spec.policy, 702U,
                       mg_requirements)),
                   "independent MG workspace requirements compile");
  OwnedField mg_vectors = make_field(32U, mg_requirements.arena_shape, 0U,
                                     1U, 101U, 201U);
  MgWorkspace mg_workspace;
  passed &= expect(static_cast<bool>(MgWorkspace::bind(
                       mg_requirements, mg_vectors.view, mg_workspace)),
                   "independent MG workspace binds");
  const std::array<HaloFieldSpec, 1U> halo_fields{{{32U, 1U, 1U}}};
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_SELF, ReductionMode::mpi_allreduce,
          requirements.reduction_capacity, reductions)) &&
          static_cast<bool>(halo.reserve(
              MPI_COMM_SELF, mg_requirements.levels[0U].patch,
              {halo_fields.data(), halo_fields.size()})),
      "MG runtime services compile around the independent Krylov workspace");
  std::vector<HaloEngine> coarse_halos(mg_requirements.level_count - 1U);
  std::vector<HaloEngine*> coarse_halo_pointers(coarse_halos.size());
  for (std::size_t level = 1U; level < mg_requirements.level_count; ++level) {
    passed &= expect(static_cast<bool>(coarse_halos[level - 1U].reserve(
                         MPI_COMM_SELF, mg_requirements.levels[level].patch,
                         {halo_fields.data(), halo_fields.size()})),
                     "coarse MG halo reserves");
    coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
  }
  NativeCartesianMgPlan plan;
  passed &= expect(
      static_cast<bool>(NativeCartesianMgPlan::compile(
          spec,
          MgRuntimeServices{
              &halo, &reductions, &mg_workspace,
              {coarse_halo_pointers.data(), coarse_halo_pointers.size()}},
          MgCoefficientViews{as_const(diagonal.view), x.view, y.view, z.view},
          plan)) &&
          plan.workspace_storage_address() != 0U &&
          plan.workspace_storage_address() !=
              krylov_workspace.vector_storage_address() &&
          plan.hierarchy_storage_address() != 0U,
      "MG owns persistent workspace distinct from every Krylov slot");
  passed &= expect(checksum(krylov_vectors.storage) == vector_before &&
                       checksum(krylov_scalars.storage) == scalar_before,
                   "MG compile leaves all Krylov vector/scalar sentinels exact");

  passed &= expect(static_cast<bool>(plan.apply(as_const(residual.view),
                                                correction.view, 0U)),
                   "MG applies using only its own persistent workspace");
  passed &= expect(checksum(krylov_vectors.storage) == vector_before &&
                       checksum(krylov_scalars.storage) == scalar_before,
                   "MG apply leaves every Krylov slot and scalar untouched");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_mg_never_borrows_krylov_slots();
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
