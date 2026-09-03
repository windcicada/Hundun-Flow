// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstddef>

namespace hundun::linear::detail {

std::size_t checked_cg_workspace_bytes(std::size_t owned_count);

}  // namespace hundun::linear::detail
