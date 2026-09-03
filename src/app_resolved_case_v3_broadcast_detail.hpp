// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case_v3.hpp"

#include <mpi.h>

namespace hundun::config {

ResolvedCaseV3 broadcast_resolved_case_v3(MPI_Comm comm, int root,
                                          const ResolvedCaseV3 *root_case);

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
namespace detail {

enum class ResolvedCaseV3BroadcastFault {
  none,
  version_tag,
  variant_tag,
  truncate_payload,
  mpi_operation
};

ResolvedCaseV3
broadcast_resolved_case_v3_with_fault(MPI_Comm comm, int root,
                                      const ResolvedCaseV3 *root_case,
                                      ResolvedCaseV3BroadcastFault fault);

} // namespace detail
#endif

} // namespace hundun::config
