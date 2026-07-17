// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/uniform_structured_mesh.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/structured_decomposition.hpp"

#include <cmath>

namespace hundun::mesh {
namespace {

bool same(runtime::Int3 lhs, runtime::Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool positive(runtime::Int3 value) noexcept {
  return value.x > 0 && value.y > 0 && value.z > 0;
}

bool positive(runtime::Real3 value) noexcept {
  return value.x > 0.0 && value.y > 0.0 && value.z > 0.0;
}

}  // namespace

UniformStructuredMesh::UniformStructuredMesh(
    runtime::Int3 global_cells, runtime::Real3 origin_m,
    runtime::Real3 length_m,
    const runtime::StructuredDecomposition& decomposition) {
  if (!positive(global_cells)) {
    throw runtime::Error("mesh global cell counts must be positive");
  }
  if (!same(global_cells, decomposition.global_extent())) {
    throw runtime::Error(
        "mesh global cell counts must match the decomposition extent");
  }
  if (!finite(origin_m)) {
    throw runtime::Error("mesh origin components must be finite");
  }
  if (!finite(length_m) || !positive(length_m)) {
    throw runtime::Error("mesh lengths must be finite and positive");
  }

  const runtime::Real3 spacing{
      length_m.x / static_cast<double>(global_cells.x),
      length_m.y / static_cast<double>(global_cells.y),
      length_m.z / static_cast<double>(global_cells.z)};
  if (!finite(spacing) || !positive(spacing)) {
    throw runtime::Error("mesh spacing must be finite and positive");
  }

  const double cell_volume = spacing.x * spacing.y * spacing.z;
  if (!std::isfinite(cell_volume) || cell_volume <= 0.0) {
    throw runtime::Error("mesh cell volume must be finite and positive");
  }

  const runtime::Real3 upper_endpoint{origin_m.x + length_m.x,
                                      origin_m.y + length_m.y,
                                      origin_m.z + length_m.z};
  if (!finite(upper_endpoint)) {
    throw runtime::Error("mesh upper endpoints must be finite");
  }

  origin_m_ = origin_m;
  spacing_m_ = spacing;
  cell_volume_m3_ = cell_volume;
  local_extent_ = decomposition.local_extent();
  owned_global_box_ = decomposition.owned_box();
}

runtime::Real3 UniformStructuredMesh::spacing_m() const noexcept {
  return spacing_m_;
}

runtime::Real3 UniformStructuredMesh::cell_center(
    runtime::Int3 local_cell) const {
  if (local_cell.x < 0 || local_cell.x >= local_extent_.x ||
      local_cell.y < 0 || local_cell.y >= local_extent_.y ||
      local_cell.z < 0 || local_cell.z >= local_extent_.z) {
    throw runtime::Error("local cell is outside the mesh owned box");
  }

  const runtime::Int3 global_cell{owned_global_box_.begin.x + local_cell.x,
                                  owned_global_box_.begin.y + local_cell.y,
                                  owned_global_box_.begin.z + local_cell.z};
  return runtime::Real3{
      origin_m_.x + (static_cast<double>(global_cell.x) + 0.5) * spacing_m_.x,
      origin_m_.y + (static_cast<double>(global_cell.y) + 0.5) * spacing_m_.y,
      origin_m_.z + (static_cast<double>(global_cell.z) + 0.5) * spacing_m_.z};
}

double UniformStructuredMesh::cell_volume_m3() const noexcept {
  return cell_volume_m3_;
}

runtime::Int3 UniformStructuredMesh::local_extent() const noexcept {
  return local_extent_;
}

runtime::Box3 UniformStructuredMesh::owned_global_box() const noexcept {
  return owned_global_box_;
}

}  // namespace hundun::mesh
