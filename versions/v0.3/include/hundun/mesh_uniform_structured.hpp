// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_types.hpp"

namespace hundun::runtime {
class StructuredDecomposition;
}

namespace hundun::mesh {

class UniformStructuredMesh final {
 public:
  UniformStructuredMesh(runtime::Int3 global_cells,
                        runtime::Real3 origin_m,
                        runtime::Real3 length_m,
                        const runtime::StructuredDecomposition&);
  runtime::Real3 spacing_m() const noexcept;
  runtime::Real3 cell_center(runtime::Int3 local_cell) const;
  double cell_volume_m3() const noexcept;
  runtime::Int3 local_extent() const noexcept;
  runtime::Box3 owned_global_box() const noexcept;

 private:
  runtime::Real3 origin_m_{};
  runtime::Real3 spacing_m_{};
  double cell_volume_m3_{};
  runtime::Int3 local_extent_{};
  runtime::Box3 owned_global_box_{};
};

}  // namespace hundun::mesh
