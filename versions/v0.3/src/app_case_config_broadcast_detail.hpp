// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/cfg_case_config.hpp"
#include "hundun/cfg_resolved_case.hpp"

#include <mpi.h>

namespace hundun::application {

config::CaseConfig broadcast_case_config(MPI_Comm comm, int root,
                                         const config::CaseConfig *root_config);

config::ResolvedCase broadcast_resolved_case(
    MPI_Comm comm, int root, const config::ResolvedCase *root_case);

} // namespace hundun::application
