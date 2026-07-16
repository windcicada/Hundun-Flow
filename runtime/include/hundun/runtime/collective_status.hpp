// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <mpi.h>

#include <optional>
#include <string>
#include <string_view>

namespace hundun::runtime {

struct CollectiveStatus {
  bool ok;
  int failing_rank;
  std::string message;
};

CollectiveStatus collective_status(MPI_Comm comm, bool local_ok,
                                   std::string_view local_message);

void require_expected_ranks(MPI_Comm comm,
                            std::optional<int> expected_ranks);

}  // namespace hundun::runtime
