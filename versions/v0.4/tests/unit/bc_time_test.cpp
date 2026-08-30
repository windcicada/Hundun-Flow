// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cfenv>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <type_traits>
#include <utility>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  if (void* pointer = std::malloc(size == 0U ? 1U : size)) {
    return pointer;
  }
  throw std::bad_alloc{};
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* pointer = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&pointer, alignment, requested) == 0 &&
      pointer != nullptr) {
    return pointer;
  }
  throw std::bad_alloc{};
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}

void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}

void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}

void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double left, double right, double tolerance = 1.0e-12) {
  return std::abs(left - right) <=
         tolerance * std::max({1.0, std::abs(left), std::abs(right)});
}

std::uint64_t double_bits(double value) {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::array<std::uint64_t, 10U> step_wire(const StepTime& step) {
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

bool same_step(const StepTime& left, const StepTime& right) {
  return step_wire(left) == step_wire(right);
}

bool collective_status_identical(Status status) {
  const std::uint64_t packed =
      static_cast<std::uint64_t>(status.detail) << 16U |
      static_cast<std::uint64_t>(status.code);
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

bool collective_step_identical(const StepTime& step) {
  const auto local = step_wire(step);
  std::array<std::uint64_t, local.size()> minimum{};
  std::array<std::uint64_t, local.size()> maximum{};
  return MPI_Allreduce(local.data(), minimum.data(),
                       static_cast<int>(local.size()), MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(local.data(), maximum.data(),
                       static_cast<int>(local.size()), MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

StepTime sentinel_step() {
  StepTime step;
  step.time = -17.25;
  step.dt = 9.5;
  step.accepted_step = 91U;
  step.attempt = 7U;
  step.origin = StepOrigin::retry;
  step.bdf = {-3.0, 4.0, -5.0, 2U};
  step.generation = 1234U;
  return step;
}

Status accept_self(TimeControllerState& state, const StepTime& step) {
  StepTime next = sentinel_step();
  return state.finish(MPI_COMM_SELF, step, Status{}, next);
}

Status retry_self(TimeControllerState& state, const StepTime& step,
                  Status cause, StepTime& retry) {
  return state.finish(MPI_COMM_SELF, step, cause, retry);
}

TimeControlSpec valid_spec(TimeControlKind control =
                               TimeControlKind::adaptive_flow,
                           TimeScheme scheme = TimeScheme::variable_bdf2) {
  TimeControlSpec spec;
  spec.control = control;
  spec.scheme = scheme;
  spec.initial_dt = 0.1;
  spec.minimum_dt = 1.0e-6;
  spec.maximum_dt = 10.0;
  spec.convective_cfl = 0.8;
  spec.viscous_cfl = 0.5;
  spec.thermal_cfl = 0.4;
  spec.species_cfl = 0.3;
  spec.acoustic_cfl = 0.2;
  spec.maximum_growth = 1.25;
  spec.retry_factor = 0.5;
  spec.maximum_retries = 3U;
  spec.minimum_bdf_ratio = 0.2;
  spec.maximum_bdf_ratio = 5.0;
  return spec;
}

LocalTimeLimits loose_limits() {
  const double infinity = std::numeric_limits<double>::infinity();
  return {infinity, infinity, infinity, infinity, infinity};
}

bool test_compile_and_local_limits() {
  bool passed = true;
  TimeSchemePlan retained;
  TimeControlSpec spec = valid_spec();
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(spec, retained)) &&
                       retained.fingerprint() != 0U,
                   "valid adaptive-flow time plan compiles");
  const PlanFingerprint retained_fingerprint = retained.fingerprint();

  const std::array<double TimeControlSpec::*, 12U> positive_members{
      &TimeControlSpec::initial_dt,
      &TimeControlSpec::minimum_dt,
      &TimeControlSpec::maximum_dt,
      &TimeControlSpec::convective_cfl,
      &TimeControlSpec::viscous_cfl,
      &TimeControlSpec::thermal_cfl,
      &TimeControlSpec::species_cfl,
      &TimeControlSpec::acoustic_cfl,
      &TimeControlSpec::maximum_growth,
      &TimeControlSpec::retry_factor,
      &TimeControlSpec::minimum_bdf_ratio,
      &TimeControlSpec::maximum_bdf_ratio};
  for (double TimeControlSpec::*member : positive_members) {
    TimeControlSpec invalid = spec;
    invalid.*member = std::numeric_limits<double>::quiet_NaN();
    passed &= expect(TimeSchemePlan::compile(invalid, retained).code ==
                             StatusCode::invalid_plan &&
                         retained.fingerprint() == retained_fingerprint,
                     "non-finite time parameter rejects atomically");
  }
  TimeControlSpec invalid = spec;
  invalid.minimum_dt = invalid.initial_dt * 2.0;
  passed &= expect(TimeSchemePlan::compile(invalid, retained).code ==
                           StatusCode::invalid_plan,
                   "minimum dt cannot exceed initial dt");
  invalid = spec;
  invalid.maximum_growth = 0.99;
  passed &= expect(TimeSchemePlan::compile(invalid, retained).code ==
                           StatusCode::invalid_plan,
                   "growth below one is rejected");
  invalid = spec;
  invalid.retry_factor = 1.0;
  passed &= expect(TimeSchemePlan::compile(invalid, retained).code ==
                           StatusCode::invalid_plan,
                   "retry factor must strictly reduce dt");
  invalid = spec;
  invalid.maximum_retries = 0U;
  passed &= expect(TimeSchemePlan::compile(invalid, retained).code ==
                           StatusCode::invalid_plan,
                   "zero retry allowance is rejected");
  invalid = spec;
  invalid.minimum_bdf_ratio = 2.0;
  invalid.maximum_bdf_ratio = 1.0;
  passed &= expect(TimeSchemePlan::compile(invalid, retained).code ==
                           StatusCode::invalid_plan,
                   "reversed BDF ratio interval is rejected");

  TimeSchemePlan flow;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(spec, flow)) &&
                       !flow.acoustic_hard_limit(),
                   "adaptive flow excludes acoustic hard limit");
  LocalTimeLimits limits{10.0, 10.0, 10.0, 10.0, 1.0e-30};
  double dt = 0.0;
  TimeLimit active = TimeLimit::fixed;
  passed &= expect(static_cast<bool>(flow.local_candidate(limits, dt, active)) &&
                       close(dt, 3.0) && active == TimeLimit::species,
                   "flow candidate applies CFL multipliers and ignores sound");

  TimeControlSpec acoustic_spec = spec;
  acoustic_spec.control = TimeControlKind::adaptive_acoustic;
  TimeSchemePlan acoustic;
  passed &= expect(static_cast<bool>(
                       TimeSchemePlan::compile(acoustic_spec, acoustic)) &&
                       acoustic.acoustic_hard_limit() &&
                       static_cast<bool>(
                           acoustic.local_candidate(limits, dt, active)) &&
                       close(dt, 2.0e-31) &&
                       active == TimeLimit::acoustic,
                   "adaptive acoustic includes the acoustic hard limit");

  TimeControlSpec fixed_spec = spec;
  fixed_spec.control = TimeControlKind::fixed;
  TimeSchemePlan fixed;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(fixed_spec, fixed)) &&
                       static_cast<bool>(fixed.local_candidate(
                           LocalTimeLimits{}, dt, active)) &&
                       close(dt, fixed_spec.initial_dt) &&
                       active == TimeLimit::fixed,
                   "fixed control does not inspect physical limits");

  limits = {1.0, 1.6, 2.0, 8.0 / 3.0, 5.0};
  passed &= expect(static_cast<bool>(flow.local_candidate(limits, dt, active)) &&
                       close(dt, 0.8) && active == TimeLimit::convective,
                   "equal physical candidates use deterministic enum priority");
  const double just_below_max =
      std::nextafter(spec.maximum_dt, 0.0) / spec.convective_cfl;
  limits = loose_limits();
  limits.convective = just_below_max;
  passed &= expect(static_cast<bool>(flow.local_candidate(limits, dt, active)) &&
                       dt < spec.maximum_dt &&
                       dt == spec.convective_cfl * just_below_max,
                   "candidate never rounds above a stability upper bound");

  const std::array<TimeLimit, 4U> expected_kinds{
      TimeLimit::convective, TimeLimit::viscous, TimeLimit::thermal,
      TimeLimit::species};
  for (std::size_t selected = 0U; selected < expected_kinds.size(); ++selected) {
    limits = {100.0, 100.0, 100.0, 100.0,
              std::numeric_limits<double>::infinity()};
    double* const values[4]{&limits.convective, &limits.viscous,
                            &limits.thermal, &limits.species};
    *values[selected] = 1.0;
    passed &= expect(static_cast<bool>(flow.local_candidate(limits, dt, active)) &&
                         active == expected_kinds[selected],
                     "each flow stability family can own the minimum");
  }
  limits.convective = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(flow.local_candidate(limits, dt, active).code ==
                           StatusCode::numerical_failure,
                   "invalid active physical limit is numerical failure");
  limits = loose_limits();
  limits.acoustic = std::numeric_limits<double>::quiet_NaN();
  passed &= expect(static_cast<bool>(flow.local_candidate(limits, dt, active)),
                   "adaptive flow never validates unused acoustic limit");
  return passed;
}

bool test_state_and_bdf() {
  TimeSchemePlan plan;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(
                           valid_spec(), plan)),
                       "BDF state plan compiles");
  TimeControllerState state;
  passed &= expect(static_cast<bool>(
                       TimeControllerState::start(plan, 2.0, state)),
                   "fresh controller starts");
  LocalTimeLimits limits = loose_limits();
  limits.convective = 1.0;
  StepTime first;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_SELF, limits,
                                                   first)) &&
                       first.origin == StepOrigin::fresh_start &&
                       first.accepted_step == 0U && first.attempt == 0U &&
                       close(first.time, 2.0) && close(first.dt, 0.1) &&
                       first.bdf.order == 1U &&
                       close(first.bdf.a0, 10.0) &&
                       close(first.bdf.a1, -10.0) && first.bdf.a2 == 0.0,
                   "fresh proposal uses capped initial dt and BE");
  StepTime ignored;
  passed &= expect(state.propose(MPI_COMM_SELF, limits, ignored).code ==
                           StatusCode::invalid_plan,
                   "only one proposal may be active");
  passed &= expect(static_cast<bool>(accept_self(state, first)) &&
                       close(state.time(), 2.1) &&
                       state.accepted_step() == 1U &&
                       close(state.last_accepted_dt(), 0.1),
                   "accept advances authoritative time once");
  passed &= expect(accept_self(state, first).code == StatusCode::invalid_plan,
                   "stale proposal cannot be accepted twice");

  StepTime second;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_SELF, limits,
                                                   second)) &&
                       second.origin == StepOrigin::accepted &&
                       close(second.dt, 0.125) && second.bdf.order == 2U,
                   "second step obeys growth and enables variable BDF2");
  const double ratio = second.dt / first.dt;
  passed &= expect(close(second.bdf.a0,
                         (1.0 + 2.0 * ratio) /
                             ((1.0 + ratio) * second.dt)) &&
                       close(second.bdf.a1,
                             -(1.0 + ratio) / second.dt) &&
                       close(second.bdf.a2,
                             ratio * ratio /
                                 ((1.0 + ratio) * second.dt)) &&
                       close(second.bdf.a0 + second.bdf.a1 + second.bdf.a2,
                             0.0),
                   "variable-step BDF2 coefficients are analytic");

  StepTime retry;
  passed &= expect(static_cast<bool>(retry_self(
                       state,
                       second, Status{StatusCode::numerical_failure, 9U},
                       retry)) &&
                       retry.origin == StepOrigin::retry && retry.attempt == 1U &&
                       close(retry.time, state.time()) &&
                       close(retry.dt, second.dt * 0.5) &&
                       retry.bdf.order == 1U && state.retry_count() == 1U &&
                       state.has_active_proposal(),
                   "retry reduces dt, preserves time, and forces BE");
  passed &= expect(retry_self(state, second,
                              Status{StatusCode::numerical_failure, 9U},
                              ignored)
                           .code == StatusCode::invalid_plan,
                   "stale rejected proposal cannot mutate retry state");
  passed &= expect(static_cast<bool>(accept_self(state, retry)) &&
                       state.retry_count() == 0U,
                   "accepted retry clears retry counter");

  StepTime after_retry;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_SELF, limits,
                                                   after_retry)) &&
                       after_retry.bdf.order == 2U,
                   "new step after accepted retry may use BDF2 history");
  passed &= expect(static_cast<bool>(accept_self(state, after_retry)),
                   "post-retry proposal accepts");

  TimeControllerState restarted;
  passed &= expect(static_cast<bool>(TimeControllerState::restart(
                       plan, state.time(), state.last_accepted_dt(),
                       state.accepted_step(), restarted)),
                   "restart restores scalar accepted metadata");
  StepTime restart_first;
  passed &= expect(static_cast<bool>(restarted.propose(
                       MPI_COMM_SELF, limits, restart_first)) &&
                       restart_first.origin == StepOrigin::restart &&
                       restart_first.bdf.order == 1U,
                   "restart first step always uses BE");
  TimeControllerState retained_restart = restarted;
  passed &= expect(TimeControllerState::restart(
                       plan, state.time(), plan.spec().maximum_dt * 2.0,
                       state.accepted_step(), retained_restart)
                           .code == StatusCode::invalid_plan &&
                       retained_restart.accepted_step() ==
                           restarted.accepted_step(),
                   "restart rejects out-of-plan dt atomically");
  passed &= expect(
      TimeControllerState::restart(
          plan, state.time(), state.last_accepted_dt(),
          std::numeric_limits<std::uint64_t>::max(), retained_restart)
              .code == StatusCode::invalid_plan &&
          retained_restart.accepted_step() == restarted.accepted_step(),
      "restart rejects an accepted-step counter that cannot advance");
  passed &= expect(
      TimeControllerState::restart(plan,
                                   std::numeric_limits<double>::max(),
                                   state.last_accepted_dt(),
                                   state.accepted_step(), retained_restart)
                  .code == StatusCode::invalid_plan &&
          retained_restart.accepted_step() == restarted.accepted_step(),
      "restart rejects a time whose next accepted dt cannot advance it");

  TimeControlSpec narrow = valid_spec();
  narrow.minimum_bdf_ratio = 0.9;
  narrow.maximum_bdf_ratio = 1.1;
  TimeSchemePlan narrow_plan;
  TimeControllerState narrow_state;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(
                       narrow, narrow_plan)) &&
                       static_cast<bool>(TimeControllerState::start(
                           narrow_plan, 0.0, narrow_state)),
                   "narrow ratio plan starts");
  StepTime narrow_first;
  StepTime narrow_second;
  passed &= expect(static_cast<bool>(narrow_state.propose(
                       MPI_COMM_SELF, limits, narrow_first)) &&
                       static_cast<bool>(accept_self(narrow_state,
                                                     narrow_first)) &&
                       static_cast<bool>(narrow_state.propose(
                           MPI_COMM_SELF, limits, narrow_second)) &&
                       narrow_second.bdf.order == 1U,
                   "out-of-range dt ratio safely falls back to BE");

  TimeSchemePlan be_plan;
  TimeControllerState be_state;
  TimeControlSpec be_spec = valid_spec(TimeControlKind::adaptive_flow,
                                       TimeScheme::backward_euler);
  StepTime be_first;
  StepTime be_second;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(be_spec, be_plan)) &&
                       static_cast<bool>(TimeControllerState::start(
                           be_plan, 0.0, be_state)) &&
                       static_cast<bool>(be_state.propose(
                           MPI_COMM_SELF, limits, be_first)) &&
                       static_cast<bool>(accept_self(be_state, be_first)) &&
                       static_cast<bool>(be_state.propose(
                           MPI_COMM_SELF, limits, be_second)) &&
                       be_second.bdf.order == 1U,
                   "backward-Euler scheme never upgrades to BDF2");
  return passed;
}

bool test_prepare_then_commit_time_finish() {
  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(
                           valid_spec(TimeControlKind::fixed), plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 4.0, state)),
                       "two-phase time controller starts");
  StepTime proposal;
  passed &= expect(static_cast<bool>(state.propose(
                       MPI_COMM_WORLD, LocalTimeLimits{}, proposal)),
                   "two-phase accept proposal publishes");

  const double time_before = state.time();
  const double dt_before = state.last_accepted_dt();
  const std::uint64_t step_before = state.accepted_step();
  const std::uint32_t retry_before = state.retry_count();
  PreparedTimeFinish accept;
  const Status accept_prepared =
      state.prepare_finish(MPI_COMM_WORLD, proposal, Status{}, accept);
  passed &= expect(
      static_cast<bool>(accept_prepared) && accept.valid() &&
          accept.decision() == TimeFinishDecision::accept &&
          static_cast<bool>(accept.outcome()) && accept.lowest_failing_rank() == -1 &&
          close(accept.candidate().time, proposal.time + proposal.dt),
      "accept preparation returns one bound accept certificate");
  passed &= expect(
      state.time() == time_before && state.last_accepted_dt() == dt_before &&
          state.accepted_step() == step_before &&
          state.retry_count() == retry_before && state.has_active_proposal(),
      "accept preparation leaves every controller state value unchanged");
  state.commit_accept(accept);
  passed &= expect(!accept.valid() && close(state.time(), time_before + proposal.dt) &&
                       close(state.last_accepted_dt(), proposal.dt) &&
                       state.accepted_step() == step_before + 1U &&
                       state.retry_count() == 0U &&
                       !state.has_active_proposal(),
                   "accept commit alone advances authoritative time");

  StepTime retry_proposal;
  passed &= expect(static_cast<bool>(state.propose(
                       MPI_COMM_WORLD, LocalTimeLimits{}, retry_proposal)),
                   "two-phase retry proposal publishes");
  const double retry_time_before = state.time();
  const std::uint64_t retry_step_before = state.accepted_step();
  const int rank = [] {
    int value = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &value);
    return value;
  }();
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int failing_rank = size == 1 ? 0 : size - 1;
  const Status local = rank == failing_rank
                           ? Status{StatusCode::numerical_failure, 991U}
                           : Status{};
  PreparedTimeFinish retry;
  const Status retry_prepared =
      state.prepare_finish(MPI_COMM_WORLD, retry_proposal, local, retry);
  passed &= expect(
      static_cast<bool>(retry_prepared) && retry.valid() &&
          retry.decision() == TimeFinishDecision::retry &&
          retry.outcome().code == StatusCode::numerical_failure &&
          retry.outcome().detail == 991U &&
          retry.lowest_failing_rank() == failing_rank &&
          collective_status_identical(retry.outcome()) &&
          collective_step_identical(retry.candidate()),
      "one-rank numerical failure prepares one identical retry decision");
  passed &= expect(state.time() == retry_time_before &&
                       state.accepted_step() == retry_step_before &&
                       state.retry_count() == 0U &&
                       state.has_active_proposal(),
                   "retry preparation leaves the active proposal unchanged");
  const StepTime expected_retry = retry.candidate();
  StepTime next = sentinel_step();
  state.commit_retry(retry, next);
  passed &= expect(!retry.valid() && state.time() == retry_time_before &&
                       state.accepted_step() == retry_step_before &&
                       state.retry_count() == 1U &&
                       state.has_active_proposal() &&
                       same_step(next, expected_retry) &&
                       next.origin == StepOrigin::retry,
                   "retry commit alone publishes the reduced active proposal");

  PreparedTimeFinish fatal;
  const Status fatal_prepared = state.prepare_finish(
      MPI_COMM_WORLD, next,
      rank == failing_rank ? Status{StatusCode::io_failure, 992U} : Status{},
      fatal);
  passed &= expect(static_cast<bool>(fatal_prepared) && fatal.valid() &&
                       fatal.decision() == TimeFinishDecision::fatal &&
                       fatal.outcome().code == StatusCode::io_failure &&
                       fatal.lowest_failing_rank() == failing_rank &&
                       state.has_active_proposal() &&
                       state.retry_count() == 1U,
                   "fatal preparation reports failure without mutation");
  return passed;
}

bool test_retry_exhaustion_and_minimum() {
  TimeControlSpec spec = valid_spec();
  spec.minimum_dt = 0.04;
  spec.initial_dt = 0.1;
  spec.maximum_retries = 1U;
  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(spec, plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 0.0, state)),
                       "retry-bound controller starts");
  LocalTimeLimits limits = loose_limits();
  limits.convective = 1.0;
  StepTime proposal;
  StepTime retry;
  passed &= expect(static_cast<bool>(
                       state.propose(MPI_COMM_SELF, limits, proposal)) &&
                       static_cast<bool>(retry_self(
                           state,
                           proposal,
                           Status{StatusCode::rejected_step, 77U}, retry)),
                   "one configured retry is produced");
  PreparedTimeFinish nonretry;
  const Status nonretry_prepared = state.prepare_finish(
      MPI_COMM_SELF, retry, Status{StatusCode::io_failure, 92U}, nonretry);
  passed &= expect(static_cast<bool>(nonretry_prepared) && nonretry.valid() &&
                       nonretry.decision() == TimeFinishDecision::fatal &&
                       nonretry.outcome().code == StatusCode::io_failure &&
                       state.has_active_proposal() &&
                       state.retry_count() == 1U,
                   "fatal preparation remains non-mutating until explicit commit");
  StepTime retained = retry;
  const Status exhausted = retry_self(
      state,
      retry, Status{StatusCode::rejected_step, 78U}, retained);
  passed &= expect(exhausted.code == StatusCode::rejected_step &&
                       close(state.time(), 0.0) &&
                       state.accepted_step() == 0U &&
                       !state.has_active_proposal() &&
                       state.retry_count() == 0U &&
                       retained.generation == retry.generation,
                   "one-phase retry exhaustion retires the failed proposal");
  StepTime recovery;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_SELF, limits,
                                                   recovery)) &&
                       recovery.origin == StepOrigin::retry &&
                       recovery.attempt == 0U &&
                       recovery.bdf.order == 1U &&
                       recovery.generation > retry.generation,
                   "one-phase fatal finish remains reusable through a fresh BE proposal");

  TimeControlSpec below = valid_spec();
  below.initial_dt = 0.05;
  below.minimum_dt = 0.04;
  below.retry_factor = 0.5;
  TimeSchemePlan below_plan;
  TimeControllerState below_state;
  StepTime below_proposal;
  StepTime below_retry;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(
                       below, below_plan)) &&
                       static_cast<bool>(TimeControllerState::start(
                           below_plan, 0.0, below_state)) &&
                       static_cast<bool>(below_state.propose(
                           MPI_COMM_SELF, limits, below_proposal)) &&
                       retry_self(
                           below_state,
                           below_proposal,
                           Status{StatusCode::rejected_step, 79U}, below_retry)
                               .code == StatusCode::rejected_step &&
                       !below_state.has_active_proposal() &&
                       below_state.retry_count() == 0U,
                   "retry below minimum dt retires the proposal without unsafe clamping");
  return passed;
}

bool test_fatal_consume_arms_collective_be_recovery() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int failing_rank = size - 1;

  TimeControlSpec spec = valid_spec(TimeControlKind::fixed);
  spec.maximum_retries = 1U;
  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(spec, plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 3.0, state)),
                       "fatal-recovery controller starts");

  StepTime accepted_proposal;
  StepTime ignored = sentinel_step();
  passed &= expect(static_cast<bool>(state.propose(
                       MPI_COMM_WORLD, LocalTimeLimits{}, accepted_proposal)) &&
                       static_cast<bool>(state.finish(
                           MPI_COMM_WORLD, accepted_proposal, Status{},
                           ignored)),
                   "fatal-recovery fixture establishes accepted BDF history");
  const double accepted_time = state.time();
  const double accepted_dt = state.last_accepted_dt();
  const std::uint64_t accepted_step = state.accepted_step();

  StepTime bdf2_proposal;
  StepTime retry;
  const Status first_failure =
      rank == failing_rank
          ? Status{StatusCode::numerical_failure, 881U}
          : Status{};
  passed &= expect(static_cast<bool>(state.propose(
                       MPI_COMM_WORLD, LocalTimeLimits{}, bdf2_proposal)) &&
                       bdf2_proposal.bdf.order == 2U &&
                       static_cast<bool>(state.finish(
                           MPI_COMM_WORLD, bdf2_proposal, first_failure,
                           retry)) &&
                       retry.bdf.order == 1U && retry.attempt == 1U,
                   "one failed BDF2 proposal publishes its allowed BE retry");

  const Status exhausted_failure =
      rank == failing_rank ? Status{StatusCode::rejected_step, 882U}
                           : Status{};
  PreparedTimeFinish fatal;
  const Status prepared = state.prepare_finish(
      MPI_COMM_WORLD, retry, exhausted_failure, fatal);
  passed &= expect(static_cast<bool>(prepared) && fatal.valid() &&
                       fatal.decision() == TimeFinishDecision::fatal &&
                       fatal.outcome().code == StatusCode::rejected_step &&
                       fatal.lowest_failing_rank() == failing_rank &&
                       collective_status_identical(fatal.outcome()) &&
                       state.has_active_proposal() &&
                       state.retry_count() == 1U,
                   "retry exhaustion prepares one collective fatal ticket without mutation");

  state.commit_fatal(fatal);
  passed &= expect(!fatal.valid() && !state.has_active_proposal() &&
                       state.retry_count() == 0U &&
                       state.time() == accepted_time &&
                       state.last_accepted_dt() == accepted_dt &&
                       state.accepted_step() == accepted_step,
                   "fatal commit retires only attempt-local time state");

  StepTime recovery;
  passed &= expect(static_cast<bool>(state.propose(
                       MPI_COMM_WORLD, LocalTimeLimits{}, recovery)) &&
                       recovery.time == accepted_time &&
                       recovery.accepted_step == accepted_step &&
                       recovery.attempt == 0U &&
                       recovery.origin == StepOrigin::retry &&
                       recovery.bdf.order == 1U &&
                       recovery.generation > retry.generation &&
                       collective_step_identical(recovery),
                   "the next collective proposal is a fresh non-replayed BE recovery");
  return passed;
}

bool test_fatal_ticket_is_owner_bound_and_one_shot() {
  TimeSchemePlan plan;
  TimeControllerState owner;
  TimeControllerState foreign_owner;
  bool passed = expect(
      static_cast<bool>(TimeSchemePlan::compile(
          valid_spec(TimeControlKind::fixed), plan)) &&
          static_cast<bool>(TimeControllerState::start(plan, 5.0, owner)) &&
          static_cast<bool>(
              TimeControllerState::start(plan, 17.0, foreign_owner)),
      "fatal-ticket owner fixtures start");

  StepTime proposal;
  StepTime foreign_proposal;
  passed &= expect(static_cast<bool>(owner.propose(
                       MPI_COMM_SELF, LocalTimeLimits{}, proposal)) &&
                       static_cast<bool>(foreign_owner.propose(
                           MPI_COMM_SELF, LocalTimeLimits{},
                           foreign_proposal)),
                   "fatal-ticket owners publish independent proposals");

  PreparedTimeFinish fatal;
  PreparedTimeFinish stale;
  passed &= expect(static_cast<bool>(owner.prepare_finish(
                       MPI_COMM_SELF, proposal,
                       Status{StatusCode::io_failure, 883U}, fatal)) &&
                       static_cast<bool>(owner.prepare_finish(
                           MPI_COMM_SELF, proposal,
                           Status{StatusCode::io_failure, 883U}, stale)) &&
                       fatal.valid() && stale.valid() &&
                       fatal.decision() == TimeFinishDecision::fatal &&
                       stale.decision() == TimeFinishDecision::fatal,
                   "one active proposal can issue equivalent bound fatal certificates");

  foreign_owner.commit_fatal(fatal);
  passed &= expect(fatal.valid() && owner.has_active_proposal() &&
                       foreign_owner.has_active_proposal() &&
                       owner.time() == 5.0 && foreign_owner.time() == 17.0,
                   "a foreign controller cannot consume or apply a fatal certificate");

  owner.commit_fatal(fatal);
  passed &= expect(!fatal.valid() && !owner.has_active_proposal() &&
                       owner.time() == 5.0 &&
                       owner.last_accepted_dt() == 0.0 &&
                       owner.accepted_step() == 0U &&
                       owner.retry_count() == 0U,
                   "the owning controller consumes its fatal certificate exactly once");

  StepTime recovery;
  passed &= expect(static_cast<bool>(owner.propose(
                       MPI_COMM_SELF, LocalTimeLimits{}, recovery)) &&
                       recovery.origin == StepOrigin::retry &&
                       recovery.bdf.order == 1U &&
                       recovery.generation > proposal.generation,
                   "the consumed fatal ticket cannot alias its recovery proposal");
  owner.commit_fatal(fatal);
  passed &= expect(owner.has_active_proposal() && !fatal.valid(),
                   "replaying a consumed fatal certificate is a no-op");
  owner.commit_fatal(stale);
  passed &= expect(owner.has_active_proposal() && stale.valid(),
                   "an unconsumed but stale fatal certificate is a no-op");

  StepTime ignored = sentinel_step();
  passed &= expect(static_cast<bool>(owner.finish(
                       MPI_COMM_SELF, recovery, Status{}, ignored)) &&
                       owner.accepted_step() == 1U &&
                       static_cast<bool>(foreign_owner.finish(
                           MPI_COMM_SELF, foreign_proposal, Status{}, ignored)),
                   "invalid fatal commits leave both owners' live proposals usable");
  return passed;
}

bool test_fixed_retry_growth_recovery() {
  TimeControlSpec spec = valid_spec(TimeControlKind::fixed);
  spec.initial_dt = 0.1;
  spec.maximum_growth = 1.25;
  spec.retry_factor = 0.5;
  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(spec, plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 0.0, state)),
                       "fixed recovery controller starts");
  StepTime first;
  StepTime retry;
  const LocalTimeLimits limits{};
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_SELF, limits,
                                                   first)) &&
                       static_cast<bool>(retry_self(
                           state,
                           first, Status{StatusCode::rejected_step, 91U},
                           retry)) &&
                       close(retry.dt, 0.05) &&
                       static_cast<bool>(accept_self(state, retry)),
                   "fixed step accepts a reduced retry");
  StepTime recovering;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_SELF, limits,
                                                   recovering)) &&
                       close(recovering.dt, 0.0625),
                   "fixed control recovers through growth cap");
  return passed;
}

bool test_world_collective() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(
                           valid_spec(), plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 0.0, state)),
                       "collective controller starts");
  LocalTimeLimits limits = loose_limits();
  limits.convective = static_cast<double>(rank + 1);
  StepTime proposal;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_WORLD, limits,
                                                   proposal)) &&
                       close(proposal.dt, 0.1) &&
                       collective_step_identical(proposal),
                   "collective proposal is bitwise identical in every field");
  StepTime accepted_next = sentinel_step();
  passed &= expect(static_cast<bool>(state.finish(
                       MPI_COMM_WORLD, proposal, Status{}, accepted_next)),
                   "collective proposal accepts after consensus");

  limits.convective = rank == 0
                          ? std::numeric_limits<double>::quiet_NaN()
                          : 1.0;
  StepTime retained = proposal;
  const Status failure = state.propose(MPI_COMM_WORLD, limits, retained);
  passed &= expect(failure.code == StatusCode::numerical_failure &&
                       collective_status_identical(failure) &&
                       state.lowest_failing_rank() == 0 &&
                       !state.has_active_proposal() &&
                       retained.generation == proposal.generation,
                   "one-rank invalid limit returns one failure and no publish");
  return passed;
}

bool test_collective_finish_consensus() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    return true;
  }

  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(
                           valid_spec(), plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 3.0, state)),
                       "collective-finish controller starts");
  LocalTimeLimits limits = loose_limits();
  limits.convective = 1.0;
  StepTime proposal;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_WORLD, limits,
                                                   proposal)),
                   "collective-finish proposal publishes");
  const Status local_retry =
      rank == 1 ? Status{StatusCode::numerical_failure, 701U} : Status{};
  StepTime retry = sentinel_step();
  const Status retry_status =
      state.finish(MPI_COMM_WORLD, proposal, local_retry, retry);
  passed &= expect(static_cast<bool>(retry_status) &&
                       state.lowest_failing_rank() == 1 &&
                       state.has_active_proposal() &&
                       state.retry_count() == 1U &&
                       retry.origin == StepOrigin::retry &&
                       retry.attempt == 1U &&
                       close(retry.dt, proposal.dt * 0.5) &&
                       collective_step_identical(retry),
                   "lowest failing rank makes every rank publish one retry");

  const double retained_time = state.time();
  const std::uint64_t retained_step = state.accepted_step();
  const std::uint32_t retained_retry_count = state.retry_count();
  const StepTime retained_active = retry;
  const Status local_nonretry =
      rank == 1 ? Status{StatusCode::io_failure, 702U} : Status{};
  PreparedTimeFinish nonretry;
  const Status nonretry_status = state.prepare_finish(
      MPI_COMM_WORLD, retry, local_nonretry, nonretry);
  passed &= expect(static_cast<bool>(nonretry_status) && nonretry.valid() &&
                       nonretry.decision() == TimeFinishDecision::fatal &&
                       nonretry.outcome().code == StatusCode::io_failure &&
                       nonretry.outcome().detail == 702U &&
                       collective_status_identical(nonretry.outcome()) &&
                       state.lowest_failing_rank() == 1 &&
                       state.time() == retained_time &&
                       state.accepted_step() == retained_step &&
                       state.retry_count() == retained_retry_count &&
                       state.has_active_proposal(),
                   "non-retryable outcome prepares collectively without mutation");

  StepTime mismatched = retained_active;
  if (rank == 1) {
    ++mismatched.generation;
  }
  StepTime mismatch_output = sentinel_step();
  const Status mismatch_status =
      state.finish(MPI_COMM_WORLD, mismatched, Status{}, mismatch_output);
  passed &= expect(mismatch_status.code == StatusCode::invalid_plan &&
                       collective_status_identical(mismatch_status) &&
                       state.lowest_failing_rank() == 1 &&
                       state.time() == retained_time &&
                       state.accepted_step() == retained_step &&
                       state.retry_count() == retained_retry_count &&
                       state.has_active_proposal() &&
                       same_step(mismatch_output, sentinel_step()),
                   "ticket mismatch fails collectively without mutation");

  state.commit_fatal(nonretry);
  StepTime recovery;
  passed &= expect(!nonretry.valid() && !state.has_active_proposal() &&
                       state.retry_count() == 0U &&
                       static_cast<bool>(state.propose(
                           MPI_COMM_WORLD, limits, recovery)) &&
                       recovery.origin == StepOrigin::retry &&
                       recovery.attempt == 0U &&
                       recovery.bdf.order == 1U &&
                       collective_step_identical(recovery),
                   "collective fatal commit publishes one reusable BE recovery state");
  return passed;
}

bool test_collective_finish_rounding_authority() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    return true;
  }

  const int saved_rounding = std::fegetround();
  const int requested_rounding = rank == 0 ? FE_DOWNWARD : FE_UPWARD;
  bool passed = expect(saved_rounding != -1,
                       "rank-local rounding mode is readable");

  TimeControlSpec spec = valid_spec(TimeControlKind::fixed);
  spec.initial_dt = 0.1;
  spec.minimum_dt = 1.0e-9;
  spec.retry_factor = 0.3;
  TimeSchemePlan plan;
  TimeControllerState state;
  StepTime first;
  StepTime ignored = sentinel_step();
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(spec, plan)) &&
                       static_cast<bool>(TimeControllerState::start(
                           plan, 0.1, state)) &&
                       static_cast<bool>(state.propose(
                           MPI_COMM_WORLD, LocalTimeLimits{}, first)),
                   "mixed-rounding accept proposal publishes");
  passed &= expect(std::fesetround(requested_rounding) == 0,
                   "rank-local rounding modes are established for accept");
  passed &= expect(static_cast<bool>(state.finish(
                       MPI_COMM_WORLD, first, Status{}, ignored)),
                   "all-ok finish accepts under mixed rounding modes");
  const std::uint64_t accepted_time_bits = double_bits(state.time());
  std::uint64_t accepted_min = 0U;
  std::uint64_t accepted_max = 0U;
  MPI_Allreduce(&accepted_time_bits, &accepted_min, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&accepted_time_bits, &accepted_max, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(accepted_min == accepted_max,
                   "accepted time is bitwise identical under mixed rounding");

  StepTime second;
  passed &= expect(std::fesetround(FE_TONEAREST) == 0,
                   "proposal rounding mode is normalized");
  passed &= expect(static_cast<bool>(state.propose(
                       MPI_COMM_WORLD, LocalTimeLimits{}, second)),
                   "retry authority proposal publishes");
  passed &= expect(std::fesetround(requested_rounding) == 0,
                   "rank-local rounding modes are established for retry");
  StepTime retry = sentinel_step();
  const Status local_outcome =
      rank == 1 ? Status{StatusCode::numerical_failure, 721U} : Status{};
  const Status retry_status =
      state.finish(MPI_COMM_WORLD, second, local_outcome, retry);
  passed &= expect(static_cast<bool>(retry_status) &&
                       collective_status_identical(retry_status) &&
                       collective_step_identical(retry) &&
                       state.lowest_failing_rank() == 1,
                   "retry ticket and status are bitwise collective under mixed rounding");

  if (saved_rounding != -1) {
    passed &= expect(std::fesetround(saved_rounding) == 0,
                     "rank-local rounding mode is restored");
  }
  return passed;
}

bool test_collective_propose_rounding_authority() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    return true;
  }

  const int saved_rounding = std::fegetround();
  const int requested_rounding = rank == 0 ? FE_DOWNWARD : FE_UPWARD;
  bool passed = expect(saved_rounding != -1,
                       "propose rounding mode is readable");
  TimeControlSpec spec = valid_spec(TimeControlKind::fixed);
  spec.initial_dt = 0.1;
  spec.minimum_dt = 1.0e-9;
  TimeSchemePlan plan;
  TimeControllerState state;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(spec, plan)) &&
                       static_cast<bool>(TimeControllerState::start(
                           plan, 0.1, state)) &&
                       std::fesetround(requested_rounding) == 0,
                   "mixed-rounding propose controller starts");
  StepTime proposal = sentinel_step();
  const Status proposed =
      state.propose(MPI_COMM_WORLD, LocalTimeLimits{}, proposal);
  StepTime next = sentinel_step();
  const Status finished =
      state.finish(MPI_COMM_WORLD, proposal, Status{}, next);
  passed &= expect(static_cast<bool>(proposed) &&
                       collective_status_identical(proposed) &&
                       collective_step_identical(proposal) &&
                       static_cast<bool>(finished) &&
                       collective_status_identical(finished),
                   "mixed-rounding propose publishes one bitwise ticket and finishes");
  if (saved_rounding != -1) {
    passed &= expect(std::fesetround(saved_rounding) == 0,
                     "propose rounding mode is restored");
  }
  return passed;
}

bool test_collective_plan_and_state_identity() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size < 2) {
    return true;
  }

  bool passed = true;
  LocalTimeLimits limits = loose_limits();
  limits.convective = 1.0;

  TimeControlSpec differing_spec = valid_spec();
  if (rank != 0) {
    differing_spec.maximum_growth = 1.5;
  }
  TimeSchemePlan differing_plan;
  TimeControllerState differing_plan_state;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(
                       differing_spec, differing_plan)) &&
                       static_cast<bool>(TimeControllerState::start(
                           differing_plan, 0.0, differing_plan_state)),
                   "rank-local differing plans are individually valid");
  const StepTime plan_sentinel = sentinel_step();
  StepTime plan_output = plan_sentinel;
  const Status plan_status =
      differing_plan_state.propose(MPI_COMM_WORLD, limits, plan_output);
  passed &= expect(plan_status.code == StatusCode::invalid_plan &&
                       collective_status_identical(plan_status) &&
                       differing_plan_state.lowest_failing_rank() == 1 &&
                       !differing_plan_state.has_active_proposal() &&
                       differing_plan_state.time() == 0.0 &&
                       differing_plan_state.accepted_step() == 0U &&
                       same_step(plan_output, plan_sentinel),
                   "plan fingerprint mismatch rejects collectively and atomically");

  TimeSchemePlan common_plan;
  passed &= expect(static_cast<bool>(TimeSchemePlan::compile(
                       valid_spec(), common_plan)),
                   "common identity-test plan compiles");
  TimeControllerState differing_start_state;
  const double start_time = rank == 1 ? 1.0 : 0.0;
  passed &= expect(static_cast<bool>(TimeControllerState::start(
                       common_plan, start_time, differing_start_state)),
                   "rank-local differing start states are individually valid");
  const StepTime start_sentinel = sentinel_step();
  StepTime start_output = start_sentinel;
  const Status start_status =
      differing_start_state.propose(MPI_COMM_WORLD, limits, start_output);
  passed &= expect(start_status.code == StatusCode::invalid_plan &&
                       collective_status_identical(start_status) &&
                       differing_start_state.lowest_failing_rank() == 1 &&
                       !differing_start_state.has_active_proposal() &&
                       differing_start_state.time() == start_time &&
                       differing_start_state.accepted_step() == 0U &&
                       same_step(start_output, start_sentinel),
                   "start-state mismatch rejects collectively and atomically");

  TimeControllerState differing_restart_state;
  const double restart_dt = rank == 1 ? 0.2 : 0.1;
  const std::uint64_t restart_step = rank == 1 ? 12U : 11U;
  passed &= expect(static_cast<bool>(TimeControllerState::restart(
                       common_plan, 5.0, restart_dt, restart_step,
                       differing_restart_state)),
                   "rank-local differing restart states are individually valid");
  const StepTime restart_sentinel = sentinel_step();
  StepTime restart_output = restart_sentinel;
  const Status restart_status =
      differing_restart_state.propose(MPI_COMM_WORLD, limits, restart_output);
  passed &= expect(restart_status.code == StatusCode::invalid_plan &&
                       collective_status_identical(restart_status) &&
                       differing_restart_state.lowest_failing_rank() == 1 &&
                       !differing_restart_state.has_active_proposal() &&
                       differing_restart_state.time() == 5.0 &&
                       differing_restart_state.last_accepted_dt() ==
                           restart_dt &&
                       differing_restart_state.accepted_step() ==
                           restart_step &&
                       same_step(restart_output, restart_sentinel),
                   "restart-state mismatch rejects collectively and atomically");
  return passed;
}

bool test_collective_hot_path_no_allocations() {
  TimeSchemePlan plan;
  TimeControllerState state;
  bool passed = expect(static_cast<bool>(TimeSchemePlan::compile(
                           valid_spec(TimeControlKind::fixed), plan)) &&
                           static_cast<bool>(TimeControllerState::start(
                               plan, 0.0, state)),
                       "allocation-test controller starts");
  const LocalTimeLimits limits{};
  StepTime warmup;
  const StepTime warmup_next = sentinel_step();
  StepTime accepted_output = warmup_next;
  passed &= expect(static_cast<bool>(state.propose(MPI_COMM_WORLD, limits,
                                                   warmup)) &&
                       static_cast<bool>(state.finish(
                           MPI_COMM_WORLD, warmup, Status{},
                           accepted_output)) &&
                       same_step(accepted_output, warmup_next),
                   "allocation-test collective is warmed up without writing accept output");

  constexpr std::size_t repetitions = 16U;
  std::size_t local_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < repetitions; ++iteration) {
      StepTime proposal;
      const Status proposed = state.propose(MPI_COMM_WORLD, limits, proposal);
      StepTime next = sentinel_step();
      const Status accepted = state.finish(
          MPI_COMM_WORLD, proposal, Status{}, next);
      passed &= static_cast<bool>(proposed) && static_cast<bool>(accepted) &&
                same_step(next, sentinel_step());
    }
    StepTime proposal;
    StepTime retry;
    const Status proposed = state.propose(MPI_COMM_WORLD, limits, proposal);
    const Status retried = state.finish(
        MPI_COMM_WORLD, proposal,
        Status{StatusCode::numerical_failure, 811U}, retry);
    StepTime accept_retry_output = sentinel_step();
    const Status accepted_retry = state.finish(
        MPI_COMM_WORLD, retry, Status{}, accept_retry_output);
    passed &= static_cast<bool>(proposed) && static_cast<bool>(retried) &&
              retry.origin == StepOrigin::retry &&
              static_cast<bool>(accepted_retry) &&
              same_step(accept_retry_output, sentinel_step());
    local_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  std::uint64_t global_allocations = 0U;
  const std::uint64_t local_count = local_allocations;
  MPI_Allreduce(&local_count, &global_allocations, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(global_allocations == 0U,
                   "repeated collective propose and finish perform zero C++ heap allocations");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  static_assert(!std::is_copy_constructible_v<PreparedTimeFinish> &&
                !std::is_move_constructible_v<PreparedTimeFinish>);
  static_assert(noexcept(std::declval<TimeControllerState&>().prepare_finish(
      std::declval<MPI_Comm>(), std::declval<const StepTime&>(),
      std::declval<Status>(), std::declval<PreparedTimeFinish&>())));
  static_assert(noexcept(std::declval<TimeControllerState&>().commit_accept(
      std::declval<PreparedTimeFinish&>())));
  static_assert(noexcept(std::declval<TimeControllerState&>().commit_retry(
      std::declval<PreparedTimeFinish&>(), std::declval<StepTime&>())));
  static_assert(noexcept(std::declval<TimeControllerState&>().commit_fatal(
      std::declval<PreparedTimeFinish&>())));
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_compile_and_local_limits();
  passed &= test_state_and_bdf();
  passed &= test_prepare_then_commit_time_finish();
  passed &= test_retry_exhaustion_and_minimum();
  passed &= test_fatal_consume_arms_collective_be_recovery();
  passed &= test_fatal_ticket_is_owner_bound_and_one_shot();
  passed &= test_fixed_retry_growth_recovery();
  passed &= test_world_collective();
  passed &= test_collective_finish_consensus();
  passed &= test_collective_finish_rounding_authority();
  passed &= test_collective_propose_rounding_authority();
  passed &= test_collective_plan_and_state_identity();
  passed &= test_collective_hot_path_no_allocations();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 time-control tests passed\n";
  return 0;
}
