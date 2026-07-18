// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/case_config.hpp"

#include <mpi.h>

namespace hundun::application::detail {

enum class CaseConfigPathTransfer {
  restart_read_directory,
  restart_write_directory,
  output_directory
};

enum class CaseConfigPathFailurePhase { preparation, reconstruction };

struct CaseConfigPathFailureInjection final {
  CaseConfigPathTransfer transfer;
  CaseConfigPathFailurePhase phase;
  int rank;
};

config::CaseConfig broadcast_case_config_with_path_failure(
    MPI_Comm comm, int root, const config::CaseConfig *root_config,
    CaseConfigPathFailureInjection injection);

} // namespace hundun::application::detail
