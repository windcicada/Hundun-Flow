// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

namespace hundun::diagnostics {

struct MeshDiagnosticVertexV2 final {
  std::uint64_t global_id{};
  runtime::Int3 logical{};
  runtime::Real3 position_m{};
};

struct MeshDiagnosticCellV2 final {
  std::uint64_t global_id{};
  runtime::Real3 centre_m{};
  double volume_m3{};
  double minimum_jacobian_m3{};
  runtime::Real3 closure_m2{};
};

struct MeshDiagnosticFaceV2 final {
  std::uint64_t global_id{};
  mesh::FaceAxis axis{};
  runtime::Int3 logical{};
  std::uint64_t owner_global_cell{};
  std::optional<std::uint64_t> neighbour_global_cell;
  std::optional<std::uint32_t> patch_id;
  std::optional<std::uint64_t> periodic_pair;
  std::array<std::uint64_t, 4> vertex_ids{};
  runtime::Real3 centre_m{};
  runtime::Real3 owner_area_vector_m2{};
  std::optional<runtime::Real3> neighbour_area_vector_m2;
  double area_m2{};
  double skewness{};
  double non_orthogonality_degrees{};
};

struct MeshDiagnosticV2 final {
  int rank{};
  int rank_count{};
  runtime::Int3 global_extent{};
  runtime::Box3 owned_box{};
  mesh::MappingKind mapping_kind{mesh::MappingKind::uniform_box};
  runtime::Real3 origin_m{};
  runtime::Real3 length_m{};
  std::vector<MeshDiagnosticVertexV2> vertices;
  std::vector<MeshDiagnosticCellV2> cells;
  std::vector<MeshDiagnosticFaceV2> faces;
  double owned_volume_sum{};
  double maximum_cell_closure_norm{};
  std::uint64_t reciprocal_check_count{};
  std::uint64_t reciprocal_failure_count{};
};

std::uint64_t mesh_diagnostic_global_vertex_id(
    runtime::Int3 global_extent, runtime::Int3 logical_vertex);

MeshDiagnosticV2 make_mesh_diagnostic_v2(
    int rank, int rank_count, const mesh::MeshTopology& topology,
    const mesh::MeshGeometry& geometry);

std::vector<std::byte> encode_mesh_diagnostic_v2(
    const MeshDiagnosticV2& value);
MeshDiagnosticV2 decode_mesh_diagnostic_v2(
    const std::vector<std::byte>& bytes);

void write_mesh_diagnostic_v2_file(const std::filesystem::path& path,
                                   const MeshDiagnosticV2& value);
MeshDiagnosticV2 read_mesh_diagnostic_v2_file(
    const std::filesystem::path& path);

void write_mesh_diagnostic_v2(
    const runtime::MpiContext& mpi, const std::filesystem::path& directory,
    const MeshDiagnosticV2& value);

}  // namespace hundun::diagnostics
