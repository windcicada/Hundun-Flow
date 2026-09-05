// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once

#include "hundun/v04_boundary.hpp"

namespace hundun::v04::detail {
// Local-only identity calculation for a storage-schema migration. Reuses the
// BoundaryCompiler's exact semantic hash; does not build boundary buffers.
Status boundary_identity_for_registry(const ValidatedModel& model,
                                      const BoundaryPlan& boundary,
                                      PlanFingerprint registry,
                                      PlanFingerprint& out) noexcept;
}  // namespace hundun::v04::detail
