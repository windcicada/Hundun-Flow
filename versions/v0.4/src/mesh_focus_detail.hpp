// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_case.hpp"
#include "hundun/v04_mesh.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace hundun::v04::detail {

struct AxisPayload {
  std::vector<double> faces;
  std::vector<double> centres;
  std::vector<double> widths;
  std::vector<double> inverse_widths;
  std::vector<double> inverse_centre_distances;
  bool uniform{};
  double uniform_width{};
  double uniform_inverse_width{};
};

Status generate_cartesian_faces(
    const CartesianMeshSpec& mesh, std::uint64_t max_peak_payload_bytes,
    std::array<std::vector<double>, 3U>& faces) noexcept;

Status finish_axis_metrics(std::vector<double> faces,
                           AxisPayload& out) noexcept;

Status make_mesh_patch(std::int32_t rank, std::int32_t rank_count,
                       Int3 global_cells, MeshPatch& out) noexcept;

std::uint64_t maximum_patch_cells(Int3 global_cells,
                                  Int3 process_grid) noexcept;

}  // namespace hundun::v04::detail
