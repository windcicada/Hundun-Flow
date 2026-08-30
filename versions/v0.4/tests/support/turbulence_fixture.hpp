// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04::test {

struct TurbulenceOwnedField {
  std::vector<double> storage;
  FieldView view{};
};

inline TurbulenceOwnedField make_turbulence_field(
    FieldId field, Int3 cells, RevisionToken revision,
    StorageIdentity storage_identity) {
  TurbulenceOwnedField result;
  const std::size_t count = static_cast<std::size_t>(cells.x) * cells.y * cells.z;
  result.storage.assign(count, 0.0);
  result.view = {result.storage.data(), cells, {0, 0, 0}, 1U,
                 static_cast<std::size_t>(cells.x),
                 static_cast<std::size_t>(cells.x) * cells.y, count, 0U,
                 field, revision, storage_identity, 26001U};
  return result;
}

inline void fill(TurbulenceOwnedField& field, double value) {
  std::fill(field.storage.begin(), field.storage.end(), value);
}

inline CartesianMeshSpec turbulence_mesh(GeometryKind kind,
                                         std::int32_t cells = 4) {
  CartesianMeshSpec mesh;
  mesh.kind = kind;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {cells, cells, cells};
  mesh.minimum_spacing = {0.96 / cells, 0.96 / cells, 0.96 / cells};
  mesh.max_growth_ratio =
      kind == GeometryKind::uniform ? 1.0 : 1.0 + 0.8 / cells;
  if (kind == GeometryKind::tensor_stretched) {
    mesh.has_base_spacing = true;
    mesh.base_spacing = {1.04 / cells, 1.04 / cells, 1.04 / cells};
    mesh.focus_regions.push_back({{0.3, 0.3, 0.3},
                                  {0.7, 0.7, 0.7},
                                  {0.98 / cells, 0.98 / cells,
                                   0.98 / cells}});
  }
  mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(cells) * cells * cells;
  mesh.limits.max_memory_bytes_per_rank = 1U << 27U;
  return mesh;
}

struct TurbulenceFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  ContributionRegistry contributions;
  TurbulencePlan plan;
  TurbulenceOwnedField density;
  TurbulenceOwnedField molecular;
  TurbulenceOwnedField effective;
  std::vector<VelocityGradient> gradients;

  bool initialize(const TurbulencePlanSpec& spec,
                  GeometryKind kind = GeometryKind::uniform,
                  MPI_Comm communicator = MPI_COMM_SELF,
                  std::int32_t global_cells = 4) {
    if (!CartesianGeometryCompiler::compile(
            communicator, turbulence_mesh(kind, global_cells), {}, geometry,
            patch)) {
      return false;
    }
    const std::array<FieldId, 3U> fields{0U, 1U, 2U};
    if (!contributions.configure({fields.data(), fields.size()}) ||
        !TurbulencePlan::compile(communicator, spec, geometry, patch, 2U,
                                 41U, contributions, plan) ||
        !contributions.freeze()) {
      return false;
    }
    density = make_turbulence_field(0U, patch.cells, 101U, 201U);
    molecular = make_turbulence_field(1U, patch.cells, 102U, 202U);
    effective = make_turbulence_field(2U, patch.cells, 103U, 203U);
    fill(density, 1.2);
    fill(molecular, 1.8e-5);
    gradients.resize(static_cast<std::size_t>(patch.cells.x) * patch.cells.y *
                     patch.cells.z);
    return true;
  }

  TurbulenceUpdateInput input(RevisionToken gradient_revision = 104U) const {
    return {as_const(density.view), as_const(molecular.view),
            {gradients.data(), gradients.size()}, gradient_revision};
  }
};

}  // namespace hundun::v04::test
