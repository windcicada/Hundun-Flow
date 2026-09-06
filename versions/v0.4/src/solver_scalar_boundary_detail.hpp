// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/v04_boundary.hpp"

namespace hundun::v04::detail {

// A physical ghost is a mirror/gradient extension, not an upstream control
// volume. Recover the actual boundary face value. MPI/periodic neighbours
// remain genuine donor cells and must not be averaged with the owner.
inline double scalar_upwind_donor(const BoundaryPlan& boundary,
                                  ConstFieldView value, Int3 donor,
                                  ConstFieldView boundary_velocity = {}) noexcept {
  const int coordinates[3]{donor.x,donor.y,donor.z};
  const int extents[3]{value.interior.x,value.interior.y,value.interior.z};
  for (int axis=0;axis<3;++axis) {
    const int n=coordinates[axis];
    if(n>=0 && n<extents[axis]) continue;
    const bool high=n>=extents[axis];
    const BoundaryFacePlan* face=nullptr;
    if(boundary.face(static_cast<CartesianFace>(2*axis+(high?1:0)),face) &&
       face!=nullptr && face->local_owner && !face->periodic) {
      Int3 owner=donor;
      (axis==0?owner.x:axis==1?owner.y:owner.z)=high?extents[axis]-1:0;
      // A conditional outlet's resolved outflow target was captured before
      // the solve. Its zero-normal-gradient relation follows the current
      // unknown, not that stale target. The branch stays tied to final U.
      if (boundary_velocity.base != nullptr &&
          face->flow_kind == BoundaryKind::pressure_outlet &&
          boundary.allow_backflow().data[face->flow_parameter] != 0U &&
          (high ? 1.0 : -1.0) * boundary_velocity.unchecked(owner, axis) >= 0.0)
        return value.unchecked(owner,0U);
      return 0.5*(value.unchecked(donor,0U)+value.unchecked(owner,0U));
    }
  }
  return value.unchecked(donor,0U);
}

}  // namespace hundun::v04::detail
