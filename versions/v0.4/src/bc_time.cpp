// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTimeInput = 451U;
constexpr std::uint32_t kTimeLimitValue = 452U;
constexpr std::uint32_t kTimeState = 453U;
constexpr std::uint32_t kTimeMinimum = 454U;
constexpr std::uint32_t kTimeRetry = 455U;
constexpr std::uint32_t kTimeCollective = 456U;
constexpr std::uint32_t kTimeGeneration = 457U;
constexpr std::uint32_t kTimeIdentity = 458U;

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool positive_or_infinity(double value) noexcept {
  return (std::isfinite(value) && value > 0.0) ||
         value == std::numeric_limits<double>::infinity();
}

bool valid_status_code(StatusCode code) noexcept {
  return static_cast<std::uint16_t>(code) <=
         static_cast<std::uint16_t>(StatusCode::io_failure);
}

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
  return hash;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

Status collective_status(MPI_Comm communicator, Status local,
                         int& lowest) noexcept {
  if (communicator == MPI_COMM_NULL) {
    lowest = -1;
    return {StatusCode::invalid_plan, kTimeCollective};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    lowest = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  if (!valid_status_code(local.code)) {
    local = {StatusCode::invalid_plan, kTimeState};
  }
  const std::uint64_t candidate =
      local ? std::numeric_limits<std::uint64_t>::max()
            : static_cast<std::uint64_t>(rank);
  std::uint64_t selected = std::numeric_limits<std::uint64_t>::max();
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  if (selected == std::numeric_limits<std::uint64_t>::max()) {
    lowest = -1;
    return {};
  }
  if (selected >= static_cast<std::uint64_t>(size)) {
    lowest = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  const int source = static_cast<int>(selected);
  std::uint32_t wire[2]{};
  if (rank == source) {
    wire[0] = static_cast<std::uint32_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire, 2, MPI_UINT32_T, source, communicator) != MPI_SUCCESS ||
      wire[0] > static_cast<std::uint32_t>(StatusCode::io_failure)) {
    lowest = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  lowest = source;
  return {static_cast<StatusCode>(wire[0]), wire[1]};
}

Status collective_finish_decision(MPI_Comm communicator, Status local,
                                  std::uint32_t priority,
                                  int& lowest,
                                  std::uint32_t& selected_priority) noexcept {
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    lowest = -1;
    selected_priority = std::numeric_limits<std::uint32_t>::max();
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  constexpr std::uint64_t kNoDecision =
      std::numeric_limits<std::uint64_t>::max();
  const std::uint64_t candidate =
      local.code == StatusCode::ok
          ? kNoDecision
          : (static_cast<std::uint64_t>(priority) << 32U) |
                static_cast<std::uint32_t>(rank);
  std::uint64_t selected = kNoDecision;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    selected_priority = std::numeric_limits<std::uint32_t>::max();
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  if (selected == kNoDecision) {
    lowest = -1;
    selected_priority = std::numeric_limits<std::uint32_t>::max();
    return {};
  }
  selected_priority = static_cast<std::uint32_t>(selected >> 32U);
  const std::uint32_t selected_rank =
      static_cast<std::uint32_t>(selected & UINT64_C(0xffffffff));
  if (selected_priority > 2U ||
      selected_rank >= static_cast<std::uint32_t>(size)) {
    lowest = -1;
    selected_priority = std::numeric_limits<std::uint32_t>::max();
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  const int source = static_cast<int>(selected_rank);
  std::uint32_t wire[2]{};
  if (rank == source) {
    wire[0] = static_cast<std::uint32_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire, 2, MPI_UINT32_T, source, communicator) != MPI_SUCCESS ||
      wire[0] > static_cast<std::uint32_t>(StatusCode::io_failure)) {
    lowest = -1;
    selected_priority = std::numeric_limits<std::uint32_t>::max();
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  lowest = source;
  return {static_cast<StatusCode>(wire[0]), wire[1]};
}

bool same_double(double left, double right) noexcept {
  return double_bits(left) == double_bits(right);
}

std::array<std::uint64_t, 10U> step_wire(const StepTime& step) noexcept {
  return {double_bits(step.time),
          double_bits(step.dt),
          step.accepted_step,
          static_cast<std::uint64_t>(step.attempt),
          static_cast<std::uint64_t>(step.origin),
          double_bits(step.bdf.a0),
          double_bits(step.bdf.a1),
          double_bits(step.bdf.a2),
          static_cast<std::uint64_t>(step.bdf.order),
          step.generation};
}

bool decode_step_wire(const std::array<std::uint64_t, 10U>& wire,
                      StepTime& step) noexcept {
  if (wire[3] > std::numeric_limits<std::uint32_t>::max() ||
      wire[4] > static_cast<std::uint64_t>(StepOrigin::retry) ||
      wire[8] > std::numeric_limits<std::uint8_t>::max()) {
    return false;
  }
  std::memcpy(&step.time, &wire[0], sizeof(step.time));
  std::memcpy(&step.dt, &wire[1], sizeof(step.dt));
  step.accepted_step = wire[2];
  step.attempt = static_cast<std::uint32_t>(wire[3]);
  step.origin = static_cast<StepOrigin>(wire[4]);
  std::memcpy(&step.bdf.a0, &wire[5], sizeof(step.bdf.a0));
  std::memcpy(&step.bdf.a1, &wire[6], sizeof(step.bdf.a1));
  std::memcpy(&step.bdf.a2, &wire[7], sizeof(step.bdf.a2));
  step.bdf.order = static_cast<std::uint8_t>(wire[8]);
  step.generation = wire[9];
  return true;
}

}  // namespace

Status TimeSchemePlan::compile(const TimeControlSpec& spec,
                               TimeSchemePlan& out) noexcept {
  const bool valid_control =
      static_cast<std::uint8_t>(spec.control) <=
      static_cast<std::uint8_t>(TimeControlKind::adaptive_acoustic);
  const bool valid_scheme =
      static_cast<std::uint8_t>(spec.scheme) <=
      static_cast<std::uint8_t>(TimeScheme::variable_bdf2);
  if (!valid_control || !valid_scheme || !finite_positive(spec.initial_dt) ||
      !finite_positive(spec.minimum_dt) ||
      !finite_positive(spec.maximum_dt) ||
      spec.minimum_dt > spec.initial_dt ||
      spec.initial_dt > spec.maximum_dt ||
      !finite_positive(spec.convective_cfl) ||
      !finite_positive(spec.viscous_cfl) ||
      !finite_positive(spec.thermal_cfl) ||
      !finite_positive(spec.species_cfl) ||
      !finite_positive(spec.acoustic_cfl) ||
      !std::isfinite(spec.maximum_growth) || spec.maximum_growth < 1.0 ||
      !std::isfinite(spec.retry_factor) || !(spec.retry_factor > 0.0) ||
      !(spec.retry_factor < 1.0) || spec.maximum_retries == 0U ||
      !finite_positive(spec.minimum_bdf_ratio) ||
      !finite_positive(spec.maximum_bdf_ratio) ||
      spec.minimum_bdf_ratio > spec.maximum_bdf_ratio) {
    return {StatusCode::invalid_plan, kTimeInput};
  }

  TimeSchemePlan candidate;
  candidate.spec_ = spec;
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = hash_mix(hash, static_cast<std::uint8_t>(spec.control));
  hash = hash_mix(hash, static_cast<std::uint8_t>(spec.scheme));
  hash = hash_mix(hash, double_bits(spec.initial_dt));
  hash = hash_mix(hash, double_bits(spec.minimum_dt));
  hash = hash_mix(hash, double_bits(spec.maximum_dt));
  hash = hash_mix(hash, double_bits(spec.convective_cfl));
  hash = hash_mix(hash, double_bits(spec.viscous_cfl));
  hash = hash_mix(hash, double_bits(spec.thermal_cfl));
  hash = hash_mix(hash, double_bits(spec.species_cfl));
  hash = hash_mix(hash, double_bits(spec.acoustic_cfl));
  hash = hash_mix(hash, double_bits(spec.maximum_growth));
  hash = hash_mix(hash, double_bits(spec.retry_factor));
  hash = hash_mix(hash, spec.maximum_retries);
  hash = hash_mix(hash, double_bits(spec.minimum_bdf_ratio));
  hash = hash_mix(hash, double_bits(spec.maximum_bdf_ratio));
  candidate.fingerprint_ = hash == 0U ? 1U : hash;
  out = candidate;
  return {};
}

Status TimeSchemePlan::local_candidate(LocalTimeLimits limits, double& dt,
                                       TimeLimit& active_limit) const noexcept {
  if (fingerprint_ == 0U) {
    return {StatusCode::invalid_plan, kTimeState};
  }
  if (spec_.control == TimeControlKind::fixed) {
    dt = spec_.initial_dt;
    active_limit = TimeLimit::fixed;
    return {};
  }

  // Each input is a unit-CFL stability scale (inverse local spectral rate).
  // The compiled targets below turn those scales into actual dt upper bounds.
  const std::array<double, 5U> scales{
      limits.convective, limits.viscous, limits.thermal, limits.species,
      limits.acoustic};
  const std::array<double, 5U> targets{
      spec_.convective_cfl, spec_.viscous_cfl, spec_.thermal_cfl,
      spec_.species_cfl, spec_.acoustic_cfl};
  const std::array<TimeLimit, 5U> kinds{
      TimeLimit::convective, TimeLimit::viscous, TimeLimit::thermal,
      TimeLimit::species, TimeLimit::acoustic};
  const std::size_t count =
      spec_.control == TimeControlKind::adaptive_acoustic ? 5U : 4U;
  double candidate = spec_.maximum_dt;
  TimeLimit candidate_kind = TimeLimit::growth;
  constexpr double tie_scale = 32.0 * std::numeric_limits<double>::epsilon();
  for (std::size_t index = 0U; index < count; ++index) {
    if (!positive_or_infinity(scales[index])) {
      return {StatusCode::numerical_failure, kTimeLimitValue};
    }
    const double limited = targets[index] * scales[index];
    if (!(limited > 0.0) || std::isnan(limited)) {
      return {StatusCode::numerical_failure, kTimeLimitValue};
    }
    const double tolerance =
        tie_scale * std::max({1.0, std::abs(limited), std::abs(candidate)});
    if (limited < candidate) {
      if (limited < candidate - tolerance) {
        candidate_kind = kinds[index];
      }
      candidate = limited;
    }
  }
  dt = candidate;
  active_limit = candidate_kind;
  return {};
}

Status TimeControllerState::start(const TimeSchemePlan& plan,
                                  double start_time,
                                  TimeControllerState& out) noexcept {
  if (plan.fingerprint() == 0U || !std::isfinite(start_time)) {
    return {StatusCode::invalid_plan, kTimeInput};
  }
  TimeControllerState candidate;
  candidate.plan_ = plan;
  candidate.time_ = start_time;
  candidate.next_origin_ = StepOrigin::fresh_start;
  candidate.force_backward_euler_ = true;
  out = candidate;
  return {};
}

Status TimeControllerState::restart(const TimeSchemePlan& plan, double time,
                                    double last_accepted_dt,
                                    std::uint64_t accepted_step,
                                    TimeControllerState& out) noexcept {
  const double first_candidate_time = time + last_accepted_dt;
  if (plan.fingerprint() == 0U || !std::isfinite(time) ||
      !finite_positive(last_accepted_dt) ||
      last_accepted_dt < plan.spec().minimum_dt ||
      last_accepted_dt > plan.spec().maximum_dt ||
      accepted_step == std::numeric_limits<std::uint64_t>::max() ||
      !std::isfinite(first_candidate_time) ||
      !(first_candidate_time > time)) {
    return {StatusCode::invalid_plan, kTimeInput};
  }
  TimeControllerState candidate;
  candidate.plan_ = plan;
  candidate.time_ = time;
  candidate.last_accepted_dt_ = last_accepted_dt;
  candidate.accepted_step_ = accepted_step;
  candidate.next_origin_ = StepOrigin::restart;
  candidate.force_backward_euler_ = true;
  out = candidate;
  return {};
}

Status TimeControllerState::coefficients(
    double dt, bool force_backward_euler,
    BdfCoefficients& out) const noexcept {
  if (!finite_positive(dt)) {
    return {StatusCode::numerical_failure, kTimeLimitValue};
  }
  BdfCoefficients candidate{1.0 / dt, -1.0 / dt, 0.0, 1U};
  if (!force_backward_euler &&
      plan_.spec().scheme == TimeScheme::variable_bdf2 &&
      finite_positive(last_accepted_dt_)) {
    const double ratio = dt / last_accepted_dt_;
    if (std::isfinite(ratio) &&
        ratio >= plan_.spec().minimum_bdf_ratio &&
        ratio <= plan_.spec().maximum_bdf_ratio) {
      candidate.a0 = (1.0 + 2.0 * ratio) / ((1.0 + ratio) * dt);
      candidate.a1 = -(1.0 + ratio) / dt;
      candidate.a2 = ratio * ratio / ((1.0 + ratio) * dt);
      candidate.order = 2U;
    }
  }
  if (!std::isfinite(candidate.a0) || !std::isfinite(candidate.a1) ||
      !std::isfinite(candidate.a2)) {
    return {StatusCode::numerical_failure, kTimeLimitValue};
  }
  out = candidate;
  return {};
}

Status TimeControllerState::propose(MPI_Comm communicator,
                                    LocalTimeLimits limits,
                                    StepTime& out) noexcept {
  if (communicator == MPI_COMM_NULL) {
    lowest_failing_rank_ = -1;
    return {StatusCode::invalid_plan, kTimeCollective};
  }
  // Rank zero is the transaction authority for the immutable plan and the
  // accepted controller state.  Compare their exact wire representation
  // before using rank-local limits, so a valid-but-different rank cannot
  // publish a divergent StepTime.  The fixed-size stack wire keeps propose's
  // hot path allocation-free.
  const std::array<std::uint64_t, 19U> identity{
      plan_.fingerprint(),
      double_bits(time_),
      double_bits(last_accepted_dt_),
      accepted_step_,
      next_generation_,
      static_cast<std::uint64_t>(retry_count_),
      static_cast<std::uint64_t>(next_origin_),
      static_cast<std::uint64_t>(force_backward_euler_),
      static_cast<std::uint64_t>(active_),
      double_bits(active_step_.time),
      double_bits(active_step_.dt),
      active_step_.accepted_step,
      static_cast<std::uint64_t>(active_step_.attempt),
      static_cast<std::uint64_t>(active_step_.origin),
      double_bits(active_step_.bdf.a0),
      double_bits(active_step_.bdf.a1),
      double_bits(active_step_.bdf.a2),
      static_cast<std::uint64_t>(active_step_.bdf.order),
      active_step_.generation};
  std::array<std::uint64_t, identity.size()> authority = identity;
  if (MPI_Bcast(authority.data(), static_cast<int>(authority.size()),
                MPI_UINT64_T, 0, communicator) != MPI_SUCCESS) {
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }

  Status local = authority != identity
                     ? Status{StatusCode::invalid_plan, kTimeIdentity}
                     : (active_
                            ? Status{StatusCode::invalid_plan, kTimeState}
                            : Status{});
  double local_dt = 0.0;
  TimeLimit local_limit = TimeLimit::fixed;
  if (local) {
    local = plan_.local_candidate(limits, local_dt, local_limit);
  }
  int lowest = -1;
  Status status = collective_status(communicator, local, lowest);
  lowest_failing_rank_ = lowest;
  if (!status) {
    return status;
  }

  double global_dt = 0.0;
  if (MPI_Allreduce(&local_dt, &global_dt, 1, MPI_DOUBLE, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  static_cast<void>(local_limit);

  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS) {
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  StepTime candidate;
  Status authority_status;
  if (rank == 0) {
    const TimeControlSpec& spec = plan_.spec();
    double proposed_dt = global_dt;
    if (last_accepted_dt_ == 0.0) {
      if (spec.initial_dt < proposed_dt) {
        proposed_dt = spec.initial_dt;
      }
    } else {
      const double growth = last_accepted_dt_ * spec.maximum_growth;
      if (growth < proposed_dt) {
        proposed_dt = growth;
      }
    }
    authority_status = finite_positive(proposed_dt)
                           ? (proposed_dt < spec.minimum_dt
                                  ? Status{StatusCode::rejected_step,
                                           kTimeMinimum}
                                  : Status{})
                           : Status{StatusCode::numerical_failure,
                                    kTimeLimitValue};
    if (authority_status) {
      candidate.time = time_;
      candidate.dt = proposed_dt;
      candidate.accepted_step = accepted_step_;
      candidate.attempt = retry_count_;
      candidate.origin = next_origin_;
      candidate.generation = next_generation_;
      authority_status =
          coefficients(candidate.dt, force_backward_euler_, candidate.bdf);
    }
  }

  std::array<std::uint64_t, 12U> publication{};
  if (rank == 0) {
    publication[0] = static_cast<std::uint64_t>(authority_status.code);
    publication[1] = authority_status.detail;
    const auto candidate_wire = step_wire(candidate);
    std::copy(candidate_wire.begin(), candidate_wire.end(),
              publication.begin() + 2);
  }
  if (MPI_Bcast(publication.data(), static_cast<int>(publication.size()),
                MPI_UINT64_T, 0, communicator) != MPI_SUCCESS) {
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  if (publication[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      publication[1] > std::numeric_limits<std::uint32_t>::max()) {
    local = {StatusCode::invalid_plan, kTimeState};
  } else {
    authority_status = {static_cast<StatusCode>(publication[0]),
                        static_cast<std::uint32_t>(publication[1])};
    std::array<std::uint64_t, 10U> candidate_wire{};
    std::copy(publication.begin() + 2, publication.end(),
              candidate_wire.begin());
    local = decode_step_wire(candidate_wire, candidate)
                ? Status{}
                : Status{StatusCode::invalid_plan, kTimeState};
  }
  if (local && authority_status.code == StatusCode::ok) {
    const bool valid_candidate =
        same_double(candidate.time, time_) && finite_positive(candidate.dt) &&
        candidate.dt >= plan_.spec().minimum_dt &&
        candidate.dt <= plan_.spec().maximum_dt &&
        candidate.accepted_step == accepted_step_ &&
        candidate.attempt == retry_count_ && candidate.origin == next_origin_ &&
        std::isfinite(candidate.bdf.a0) && std::isfinite(candidate.bdf.a1) &&
        std::isfinite(candidate.bdf.a2) &&
        candidate.bdf.order >= 1U && candidate.bdf.order <= 2U &&
        candidate.generation == next_generation_;
    if (!valid_candidate) {
      local = {StatusCode::invalid_plan, kTimeState};
    }
  }
  if (local && authority_status.code != StatusCode::ok) {
    local = authority_status;
  }
  status = collective_status(communicator, local, lowest);
  lowest_failing_rank_ = lowest;
  if (!status) {
    return status;
  }
  active_step_ = candidate;
  active_ = true;
  out = candidate;
  return {};
}

bool TimeControllerState::matches_active(const StepTime& step) const noexcept {
  return active_ && step.generation == active_step_.generation &&
         step.accepted_step == active_step_.accepted_step &&
         step.attempt == active_step_.attempt &&
         step.origin == active_step_.origin &&
         same_double(step.time, active_step_.time) &&
         same_double(step.dt, active_step_.dt) &&
         step.bdf.order == active_step_.bdf.order &&
         same_double(step.bdf.a0, active_step_.bdf.a0) &&
         same_double(step.bdf.a1, active_step_.bdf.a1) &&
         same_double(step.bdf.a2, active_step_.bdf.a2);
}

std::array<std::uint64_t, 29U> TimeControllerState::finish_identity(
    const StepTime& step) const noexcept {
  return {plan_.fingerprint(),
          double_bits(time_),
          double_bits(last_accepted_dt_),
          accepted_step_,
          next_generation_,
          static_cast<std::uint64_t>(retry_count_),
          static_cast<std::uint64_t>(next_origin_),
          static_cast<std::uint64_t>(force_backward_euler_),
          static_cast<std::uint64_t>(active_),
          double_bits(active_step_.time),
          double_bits(active_step_.dt),
          active_step_.accepted_step,
          static_cast<std::uint64_t>(active_step_.attempt),
          static_cast<std::uint64_t>(active_step_.origin),
          double_bits(active_step_.bdf.a0),
          double_bits(active_step_.bdf.a1),
          double_bits(active_step_.bdf.a2),
          static_cast<std::uint64_t>(active_step_.bdf.order),
          active_step_.generation,
          double_bits(step.time),
          double_bits(step.dt),
          step.accepted_step,
          static_cast<std::uint64_t>(step.attempt),
          static_cast<std::uint64_t>(step.origin),
          double_bits(step.bdf.a0),
          double_bits(step.bdf.a1),
          double_bits(step.bdf.a2),
          static_cast<std::uint64_t>(step.bdf.order),
          step.generation};
}

bool TimeControllerState::matches_prepared(
    const PreparedTimeFinish& prepared) const noexcept {
  return prepared.valid() && prepared.commit_ready_ &&
         prepared.owner_ == this &&
         prepared.identity_ == finish_identity(prepared.proposal_) &&
         matches_active(prepared.proposal_);
}

Status TimeControllerState::prepare_finish(
    MPI_Comm communicator, const StepTime& step, Status local_outcome,
    PreparedTimeFinish& out) noexcept {
  assert(!out.valid());
  if (out.valid()) {
    return {StatusCode::invalid_plan, kTimeState};
  }
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kTimeCollective};
  }

  // Finalization is one collective transaction.  Rank zero publishes the
  // accepted-state and proposal-ticket identity; every rank checks the exact
  // bits before an outcome is reduced or any local state is mutated.
  const std::array<std::uint64_t, 29U> identity = finish_identity(step);
  std::array<std::uint64_t, identity.size()> authority = identity;
  if (MPI_Bcast(authority.data(), static_cast<int>(authority.size()),
                MPI_UINT64_T, 0, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kTimeCollective};
  }

  const bool valid_outcome = valid_status_code(local_outcome.code);
  const bool retryable =
      local_outcome.code == StatusCode::numerical_failure ||
      local_outcome.code == StatusCode::rejected_step;
  const bool local_identity_failure =
      authority != identity || !matches_active(step) || !valid_outcome;
  Status local = authority != identity
                     ? Status{StatusCode::invalid_plan, kTimeIdentity}
                 : !matches_active(step)
                     ? Status{StatusCode::invalid_plan, kTimeGeneration}
                 : !valid_outcome
                     ? Status{StatusCode::invalid_plan, kTimeState}
                 : local_outcome;
  const std::uint32_t priority =
      local_identity_failure
          ? 0U
          : (retryable ? 2U
                       : (local_outcome.code == StatusCode::ok ? 2U : 1U));
  int lowest = -1;
  std::uint32_t selected_priority =
      std::numeric_limits<std::uint32_t>::max();
  const Status decision = collective_finish_decision(
      communicator, local, priority, lowest, selected_priority);
  if (decision.code == StatusCode::mpi_failure) {
    return decision;
  }
  const int outcome_rank = lowest;

  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  const bool accept = outcome_rank < 0;
  const bool retry = decision.code == StatusCode::numerical_failure ||
                     decision.code == StatusCode::rejected_step;
  StepTime candidate;
  Status authority_status = decision;
  if (rank == 0 && accept) {
    authority_status =
        accepted_step_ != std::numeric_limits<std::uint64_t>::max() &&
                next_generation_ != std::numeric_limits<std::uint64_t>::max()
            ? Status{}
            : Status{StatusCode::invalid_plan, kTimeState};
    if (authority_status) {
      candidate = step;
      candidate.time = time_ + step.dt;
      candidate.accepted_step = accepted_step_ + 1U;
      candidate.generation = next_generation_ + 1U;
      if (!std::isfinite(candidate.time) || !(candidate.time > time_)) {
        authority_status = {StatusCode::invalid_plan, kTimeState};
      }
    }
  } else if (rank == 0 && retry) {
    authority_status =
        retry_count_ < plan_.spec().maximum_retries &&
                next_generation_ != std::numeric_limits<std::uint64_t>::max()
            ? Status{}
            : Status{StatusCode::rejected_step, kTimeRetry};
    const double retry_dt = step.dt * plan_.spec().retry_factor;
    if (authority_status &&
        (!finite_positive(retry_dt) || retry_dt < plan_.spec().minimum_dt)) {
      authority_status = {StatusCode::rejected_step, kTimeMinimum};
    }
    if (authority_status) {
      candidate.time = time_;
      candidate.dt = retry_dt;
      candidate.accepted_step = accepted_step_;
      candidate.attempt = retry_count_ + 1U;
      candidate.origin = StepOrigin::retry;
      candidate.generation = next_generation_ + 1U;
      authority_status = coefficients(retry_dt, true, candidate.bdf);
    }
  }

  std::array<std::uint64_t, 13U> publication{};
  if (rank == 0) {
    publication[0] = static_cast<std::uint64_t>(authority_status.code);
    publication[1] = authority_status.detail;
    publication[2] = accept ? 1U : (retry ? 2U : 0U);
    const auto candidate_wire = step_wire(candidate);
    std::copy(candidate_wire.begin(), candidate_wire.end(),
              publication.begin() + 3);
  }
  if (MPI_Bcast(publication.data(), static_cast<int>(publication.size()),
                MPI_UINT64_T, 0, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kTimeCollective};
  }
  if (publication[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      publication[1] > std::numeric_limits<std::uint32_t>::max() ||
      publication[2] > 2U) {
    local = {StatusCode::invalid_plan, kTimeState};
  } else {
    authority_status = {static_cast<StatusCode>(publication[0]),
                        static_cast<std::uint32_t>(publication[1])};
    std::array<std::uint64_t, 10U> candidate_wire{};
    std::copy(publication.begin() + 3, publication.end(),
              candidate_wire.begin());
    local = decode_step_wire(candidate_wire, candidate)
                ? Status{}
                : Status{StatusCode::invalid_plan, kTimeState};
  }
  if (local && authority_status.code == StatusCode::ok &&
      publication[2] == 1U) {
    const bool valid_accept =
        std::isfinite(candidate.time) && candidate.time > time_ &&
        same_double(candidate.dt, step.dt) &&
        candidate.accepted_step == accepted_step_ + 1U &&
        candidate.attempt == step.attempt && candidate.origin == step.origin &&
        same_double(candidate.bdf.a0, step.bdf.a0) &&
        same_double(candidate.bdf.a1, step.bdf.a1) &&
        same_double(candidate.bdf.a2, step.bdf.a2) &&
        candidate.bdf.order == step.bdf.order &&
        candidate.generation == next_generation_ + 1U;
    if (!valid_accept) {
      local = {StatusCode::invalid_plan, kTimeState};
    }
  }
  if (local && authority_status.code == StatusCode::ok &&
      publication[2] == 2U) {
    const bool valid_retry =
        same_double(candidate.time, time_) && finite_positive(candidate.dt) &&
        candidate.dt >= plan_.spec().minimum_dt &&
        candidate.accepted_step == accepted_step_ &&
        candidate.attempt == retry_count_ + 1U &&
        candidate.origin == StepOrigin::retry &&
        std::isfinite(candidate.bdf.a0) &&
        std::isfinite(candidate.bdf.a1) && candidate.bdf.a0 > 0.0 &&
        candidate.bdf.a1 < 0.0 &&
        candidate.bdf.a2 == 0.0 && candidate.bdf.order == 1U &&
        candidate.generation == next_generation_ + 1U;
    if (!valid_retry) {
      local = {StatusCode::invalid_plan, kTimeState};
    }
  }
  if (local && authority_status.code != StatusCode::ok) {
    local = authority_status;
  }
  Status publish_status = collective_status(communicator, local, lowest);
  if (!publish_status) {
    if (publish_status.code == StatusCode::mpi_failure) {
      return publish_status;
    }
    out.owner_ = this;
    out.identity_ = identity;
    out.proposal_ = step;
    out.candidate_ = candidate;
    out.outcome_ = publish_status;
    out.lowest_failing_rank_ =
        authority_status.code != StatusCode::ok
            ? (outcome_rank >= 0 ? outcome_rank : 0)
            : lowest;
    out.decision_ = TimeFinishDecision::fatal;
    out.commit_ready_ =
        selected_priority != 0U && authority_status.code != StatusCode::ok;
    out.valid_ = true;
    out.consumed_ = false;
    return {};
  }

  if (publication[2] == 0U) {
    out.owner_ = this;
    out.identity_ = identity;
    out.proposal_ = step;
    out.candidate_ = candidate;
    out.outcome_ = decision;
    out.lowest_failing_rank_ = outcome_rank;
    out.decision_ = TimeFinishDecision::fatal;
    out.commit_ready_ = selected_priority != 0U;
    out.valid_ = true;
    out.consumed_ = false;
    return {};
  }
  if (publication[2] == 1U) {
    out.owner_ = this;
    out.identity_ = identity;
    out.proposal_ = step;
    out.candidate_ = candidate;
    out.outcome_ = {};
    out.lowest_failing_rank_ = -1;
    out.decision_ = TimeFinishDecision::accept;
    out.commit_ready_ = true;
    out.valid_ = true;
    out.consumed_ = false;
    return {};
  }
  out.owner_ = this;
  out.identity_ = identity;
  out.proposal_ = step;
  out.candidate_ = candidate;
  out.outcome_ = decision;
  out.lowest_failing_rank_ = outcome_rank;
  out.decision_ = TimeFinishDecision::retry;
  out.commit_ready_ = true;
  out.valid_ = true;
  out.consumed_ = false;
  return {};
}

void TimeControllerState::commit_accept(
    PreparedTimeFinish& prepared) noexcept {
  const bool valid = matches_prepared(prepared) &&
                     prepared.decision_ == TimeFinishDecision::accept;
  assert(valid);
  if (!valid) return;
  const StepTime proposal = prepared.proposal_;
  const StepTime candidate = prepared.candidate_;
  prepared.consumed_ = true;
  prepared.valid_ = false;
  time_ = candidate.time;
  last_accepted_dt_ = proposal.dt;
  accepted_step_ = candidate.accepted_step;
  next_generation_ = candidate.generation;
  retry_count_ = 0U;
  lowest_failing_rank_ = -1;
  next_origin_ = StepOrigin::accepted;
  force_backward_euler_ = false;
  active_ = false;
  active_step_ = {};
}

void TimeControllerState::commit_retry(PreparedTimeFinish& prepared,
                                       StepTime& next) noexcept {
  const bool valid = matches_prepared(prepared) &&
                     prepared.decision_ == TimeFinishDecision::retry;
  assert(valid);
  if (!valid) return;
  const StepTime candidate = prepared.candidate_;
  const int lowest = prepared.lowest_failing_rank_;
  prepared.consumed_ = true;
  prepared.valid_ = false;
  retry_count_ = candidate.attempt;
  next_generation_ = candidate.generation;
  lowest_failing_rank_ = lowest;
  next_origin_ = StepOrigin::retry;
  force_backward_euler_ = true;
  active_step_ = candidate;
  active_ = true;
  next = candidate;
}

void TimeControllerState::commit_fatal(
    PreparedTimeFinish& prepared) noexcept {
  const bool valid = matches_prepared(prepared) &&
                     prepared.decision_ == TimeFinishDecision::fatal;
  // Unlike accept/retry, fatal finalization is also the recovery seam after a
  // caller has already rolled back its attempt transaction.  Treat an
  // unrelated or replayed certificate as a no-op so this cleanup path cannot
  // introduce a second failure while preserving the valid ticket's one-shot
  // ownership contract.
  if (!valid) return;
  const int lowest = prepared.lowest_failing_rank_;
  const std::uint64_t retired_generation = prepared.proposal_.generation;
  prepared.consumed_ = true;
  prepared.valid_ = false;
  if (retired_generation != std::numeric_limits<std::uint64_t>::max()) {
    next_generation_ = retired_generation + 1U;
  }
  retry_count_ = 0U;
  lowest_failing_rank_ = lowest;
  next_origin_ = StepOrigin::retry;
  force_backward_euler_ = true;
  active_step_ = {};
  active_ = false;
}

Status TimeControllerState::finish(MPI_Comm communicator,
                                   const StepTime& step,
                                   Status local_outcome,
                                   StepTime& next) noexcept {
  PreparedTimeFinish prepared;
  const Status status =
      prepare_finish(communicator, step, local_outcome, prepared);
  if (!status) {
    lowest_failing_rank_ = -1;
    return status;
  }
  if (prepared.decision_ == TimeFinishDecision::accept) {
    commit_accept(prepared);
    return {};
  }
  if (prepared.decision_ == TimeFinishDecision::retry) {
    commit_retry(prepared, next);
    return {};
  }
  const Status outcome = prepared.outcome_;
  commit_fatal(prepared);
  return outcome;
}

}  // namespace hundun::v04
