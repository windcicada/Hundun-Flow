// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>

namespace hundun::linear::detail {

std::size_t checked_bicgstab_workspace_bytes(std::size_t owned_count);

}  // namespace hundun::linear::detail
