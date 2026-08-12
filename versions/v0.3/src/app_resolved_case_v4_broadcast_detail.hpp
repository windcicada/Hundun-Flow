// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_resolved_case_v4.hpp"

#include <mpi.h>

namespace hundun::config {

ResolvedReactingCaseV4
broadcast_resolved_reacting_case_v4(MPI_Comm comm, int root,
                                    const ResolvedReactingCaseV4 *root_case);

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
namespace detail {

enum class ResolvedReactingCaseV4BroadcastFault {
  none,
  reverse_species,
  rank_local_error
};

ResolvedReactingCaseV4 broadcast_resolved_reacting_case_v4_with_fault(
    MPI_Comm comm, int root, const ResolvedReactingCaseV4 *root_case,
    ResolvedReactingCaseV4BroadcastFault fault);

} // namespace detail
#endif

} // namespace hundun::config
