// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>

namespace hundun::flow::test {

// Private test seam. This header is not part of the installed/public ABI.
class MomentumPredictorTestAccess final {
public:
  static void reset_collective_selection_calls() noexcept;
  static std::size_t collective_selection_calls() noexcept;
};

} // namespace hundun::flow::test
