// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#pragma once
#include "hundun/v04_app.hpp"

namespace hundun::v04::detail {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void arm_application_local_allocation_failure_for_test(
    ApplicationFailurePhase phase, int rank) noexcept;
#endif
}  // namespace hundun::v04::detail
