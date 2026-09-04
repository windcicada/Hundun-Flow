// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_flow.hpp"
#include "hundun/v04_io.hpp"

#include <mpi.h>

#include <string>

namespace hundun::v04::detail {

Status runtime_committed_cfl(
    MPI_Comm communicator,
    const CommittedConvectiveCflCertificate& certificate,
    RuntimeCommittedConvectiveCflAudit& runtime) noexcept;

Status runtime_advective_cfl(
    MPI_Comm communicator,
    const MomentumAdvectiveCflCertificate& certificate,
    RuntimeAdvectiveCflAudit& runtime) noexcept;

// Validate and serialize one globally identical runtime record without
// imposing a file-lifecycle policy. Long-running observers use this seam to
// keep their append descriptor open and synchronize it at checkpoints.
Status runtime_encode_evidence_line(
    MPI_Comm communicator, const IoServicePlan& services,
    const RuntimeEvidenceRecord& record, std::string& line) noexcept;

}  // namespace hundun::v04::detail
