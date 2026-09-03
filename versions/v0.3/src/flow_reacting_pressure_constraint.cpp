// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_reacting_pressure_constraint_detail.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace hundun::flow::detail {

ReactingPressureState::ReactingPressureState(double history_p0_pa,
                                             double committed_p0_pa)
    : history_p0_pa_(history_p0_pa), committed_p0_pa_(committed_p0_pa),
      trial_p0_pa_(committed_p0_pa) {
  if (!std::isfinite(history_p0_pa_) || history_p0_pa_ <= 0.0 ||
      !std::isfinite(committed_p0_pa_) || committed_p0_pa_ <= 0.0)
    throw std::invalid_argument("reacting pressure state is invalid");
}
double ReactingPressureState::history_p0_pa() const noexcept {
  return history_p0_pa_;
}
double ReactingPressureState::committed_p0_pa() const noexcept {
  return committed_p0_pa_;
}
double ReactingPressureState::trial_p0_pa() const noexcept {
  return trial_p0_pa_;
}
bool ReactingPressureState::attempt_active() const noexcept {
  return attempt_active_;
}
void ReactingPressureState::begin_attempt() {
  if (attempt_active_)
    throw std::logic_error("reacting pressure attempt is already active");
  trial_p0_pa_ = committed_p0_pa_;
  attempt_active_ = true;
}
void ReactingPressureState::set_trial_p0_pa(double value) {
  if (!attempt_active_ || !std::isfinite(value) || value <= 0.0)
    throw std::invalid_argument("reacting trial pressure is invalid");
  trial_p0_pa_ = value;
}
bool ReactingPressureState::commit(const runtime::CollectiveStatus &status) {
  if (!attempt_active_)
    throw std::logic_error("reacting pressure commit requires an attempt");
  if (!status.ok) {
    rollback();
    return false;
  }
  history_p0_pa_ = std::exchange(committed_p0_pa_, trial_p0_pa_);
  trial_p0_pa_ = committed_p0_pa_;
  attempt_active_ = false;
  return true;
}
void ReactingPressureState::rollback() {
  if (!attempt_active_)
    throw std::logic_error("reacting pressure rollback requires an attempt");
  trial_p0_pa_ = committed_p0_pa_;
  attempt_active_ = false;
}

ReactingPressureConstraintReport attempt_reacting_pressure_constraint(
    ReactingPressureState &state, const runtime::MpiContext &mpi,
    const ReactingPressureConstraintRequest &request) {
  if (!state.attempt_active() || !std::isfinite(request.duration_s) ||
      request.duration_s <= 0.0)
    throw std::invalid_argument("reacting pressure constraint is invalid");
  ReactingPressureConstraintReport report;
  report.domain = request.domain;
  report.previous_p0_pa = state.committed_p0_pa();
  report.predictor_p0_pa = report.previous_p0_pa;
  report.corrected_p0_pa = report.previous_p0_pa;
  if (request.domain == ReactingPressureDomain::open) {
    state.set_trial_p0_pa(report.previous_p0_pa);
    return report;
  }
  double global[4]{request.local_pressure_capacity_m3_per_pa,
                   request.local_predictor_expansion_m3,
                   request.local_corrector_expansion_m3,
                   request.local_boundary_volume_outflow_m3};
  for (double value : global)
    if (!std::isfinite(value))
      throw std::invalid_argument("reacting pressure integral is non-finite");
  mpi.allreduce_fp64_in_place(global, 4U,
                              runtime::Fp64ReductionOperation::sum);
  report.global_reduction_count = 1U;
  if (global[0] <= 0.0)
    throw std::invalid_argument("reacting pressure capacity is invalid");
  report.integrated_expansion_m3 = 0.5 * (global[1] + global[2]);
  report.integrated_boundary_outflow_m3 = global[3];
  const double predictor_delta = (global[1] - global[3]) / global[0];
  const double corrected_delta =
      (report.integrated_expansion_m3 - global[3]) / global[0];
  report.predictor_p0_pa = report.previous_p0_pa + predictor_delta;
  report.corrected_p0_pa = report.previous_p0_pa + corrected_delta;
  report.midpoint_dp0_dt_pa_per_s = corrected_delta / request.duration_s;
  state.set_trial_p0_pa(report.corrected_p0_pa);
  return report;
}

} // namespace hundun::flow::detail
