// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/case_config.hpp"

#include <mpi.h>

namespace hundun::application {

config::CaseConfig broadcast_case_config(MPI_Comm comm, int root,
                                         const config::CaseConfig *root_config);

} // namespace hundun::application
