// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_flow.hpp"
#include "hundun/v04_io.hpp"

#include <mpi.h>

namespace hundun::v04::detail {

Status runtime_committed_cfl(
    MPI_Comm communicator,
    const CommittedConvectiveCflCertificate& certificate,
    RuntimeCommittedConvectiveCflAudit& runtime) noexcept;

Status runtime_advective_cfl(
    MPI_Comm communicator,
    const MomentumAdvectiveCflCertificate& certificate,
    RuntimeAdvectiveCflAudit& runtime) noexcept;

}  // namespace hundun::v04::detail
