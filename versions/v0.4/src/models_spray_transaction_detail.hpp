// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "models_spray_parcel_detail.hpp"
#include "hundun/v04_execution.hpp"

namespace hundun::v04::spray::detail {

struct SprayAttemptFinishReport {
  Status status{};
  bool committed{};
  int lowest_failing_rank{-1};
};

// Synchronous finish over the existing field transaction and parcel participants.
// The caller has already staged all physical increments, injection, deletion and
// migrated owners, and passes every active injector exactly once. There are no
// callbacks or mutations between preflight, collective prepare and publication.
// Model source admission, source algebra and the driver schedule remain separate
// responsibilities; this seam cannot certify that those have run.
Status finish_spray_attempt(MPI_Comm communicator, AttemptTransaction& gas,
    ParcelContainer& parcels, Span<DeterministicInjector* const> injectors,
    Status local_status, SprayAttemptFinishReport& report) noexcept;

}  // namespace hundun::v04::spray::detail
