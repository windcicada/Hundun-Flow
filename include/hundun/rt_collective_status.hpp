// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_mpi_context.hpp"

#include <optional>
#include <string>
#include <string_view>

namespace hundun::runtime {

struct CollectiveStatus {
  bool ok;
  int failing_rank;
  std::string message;
};

CollectiveStatus collective_status(const MpiContext& context, bool local_ok,
                                   std::string_view local_message);

void require_expected_ranks(const MpiContext& context,
                            std::optional<int> expected_ranks);

}  // namespace hundun::runtime
