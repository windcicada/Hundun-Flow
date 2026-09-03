// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_case_config.hpp"

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

config::CaseConfig broadcast_case_config_with_path_failure_raw(
    MPI_Comm comm, int root, const config::CaseConfig *root_config,
    int transfer, int phase, int rank);

inline config::CaseConfig broadcast_case_config_with_path_failure(
    MPI_Comm comm, int root, const config::CaseConfig *root_config,
    CaseConfigPathFailureInjection injection) {
  return broadcast_case_config_with_path_failure_raw(
      comm, root, root_config, static_cast<int>(injection.transfer),
      static_cast<int>(injection.phase), injection.rank);
}

} // namespace hundun::application::detail
