// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_transaction_detail.hpp"

namespace hundun::v04::spray::detail {
namespace {
bool retryable(Status status) noexcept {
  return status.code == StatusCode::numerical_failure ||
         status.code == StatusCode::rejected_step;
}
void include_failure(Status candidate, Status& local) noexcept {
  if (!candidate && (local || (retryable(local) && !retryable(candidate)))) local = candidate;
}
}

Status finish_spray_attempt(MPI_Comm communicator, AttemptTransaction& gas,
    ParcelContainer& parcels, Span<DeterministicInjector* const> injectors,
    Status local, SprayAttemptFinishReport& report) noexcept {
  report = {};
  const bool valid_span = injectors.size == 0U || injectors.data != nullptr;
  include_failure(parcels.preflight_commit(), local);
  if (!valid_span) include_failure({StatusCode::invalid_plan, 26200U}, local);
  if (valid_span) for (std::size_t i = 0; i < injectors.size; ++i) {
    if (injectors.data[i] == nullptr) {
      include_failure({StatusCode::invalid_plan, 26201U}, local);
      continue;
    }
    for (std::size_t j = 0; j < i; ++j)
      if (injectors.data[j] == injectors.data[i])
        include_failure({StatusCode::invalid_plan, 26202U}, local);
    include_failure(injectors.data[i]->preflight_commit(), local);
  }
  PreparedAttemptFinish prepared;
  const Status status = gas.collective_prepare(communicator, local, prepared);
  if (status && prepared.decision() == AttemptFinishDecision::accept) {
    // All participants are still exactly those validated above. Publication is
    // allocation-free and cannot return a later failure on only one rank.
    gas.commit_accept(prepared);
    parcels.publish_preflighted_trial();
    for (std::size_t i = 0; i < injectors.size; ++i)
      injectors.data[i]->publish_preflighted_trial();
    report.committed = true;
    return {};
  }
  report.status = status ? prepared.outcome() : status;
  report.lowest_failing_rank = status ? prepared.lowest_failing_rank() : -1;
  if (status) gas.commit_reject(prepared);
  if (parcels.trial_active()) (void)parcels.rollback_trial();
  if (valid_span) for (std::size_t i = 0; i < injectors.size; ++i)
    if (injectors.data[i] != nullptr && injectors.data[i]->trial_active())
      (void)injectors.data[i]->rollback_trial();
  return report.status;
}

}  // namespace hundun::v04::spray::detail
