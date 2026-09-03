// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_view.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace hundun::mesh {
class MeshGeometry;
}

namespace hundun::immersed {

struct WeightedDonor final {
  mesh::GlobalCellId global_cell{};
  double weight{};
};

struct ReconstructionQuality final {
  std::uint32_t rank{};
  double condition_estimate{};
  std::uint32_t halo_reach{};
  std::uint64_t pivot_fingerprint{};
};

namespace detail {
struct QuadraticReconstructionStorage;
class QuadraticReconstructionWeights;
} // namespace detail

class QuadraticReconstruction final {
public:
  static QuadraticReconstruction
  create(runtime::Real3 origin_m, runtime::Real3 normal,
         runtime::Real3 tangent1, runtime::Real3 tangent2, double scale_m,
         runtime::Int3 stencil_anchor_global_cell,
         const std::vector<runtime::Int3> &donor_global_cells,
         const mesh::MeshTopology &, const mesh::MeshGeometry &);

  double value(runtime::Real3 point_m,
               const runtime::FieldView<const double> &field,
               std::size_t component) const;
  runtime::Real3 gradient(runtime::Real3 point_m,
                          const runtime::FieldView<const double> &field,
                          std::size_t component) const;
  const ReconstructionQuality &quality() const noexcept;

private:
  friend class detail::QuadraticReconstructionWeights;

  explicit QuadraticReconstruction(
      std::shared_ptr<const detail::QuadraticReconstructionStorage> storage)
      : storage_(std::move(storage)) {}

  std::shared_ptr<const detail::QuadraticReconstructionStorage> storage_;
};

} // namespace hundun::immersed
