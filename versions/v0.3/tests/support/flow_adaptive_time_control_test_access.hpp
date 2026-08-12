// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flow_adaptive_time_control_detail.hpp"

namespace hundun::flow::detail {

using TimeControlAttemptObserverRaw =
    void (*)(void *, const MomentumTimeStencil &, const FlowState &);
void adaptive_set_post_return_mutation_raw(std::uint8_t, int) noexcept;
void adaptive_set_preflight_fault_raw(std::uint8_t, int) noexcept;
void adaptive_set_recoverable_failures_raw(std::uint32_t, int) noexcept;
void adaptive_set_recoverable_failure_reason_raw(StepFailureReason) noexcept;
void adaptive_set_post_return_iteration_value_raw(std::uint64_t) noexcept;
void adaptive_reset_faults_raw() noexcept;
std::uint8_t adaptive_preflight_category_raw(
    const TimeAdvanceReport &) noexcept;
void adaptive_trusted_tail_observation_raw(
    std::uint64_t &, std::uint64_t &, std::uint64_t &,
    std::uint64_t &) noexcept;
void adaptive_set_raw_fault(std::size_t, int) noexcept;
std::size_t adaptive_raw_operation_count() noexcept;
void adaptive_raw_fault_observation(std::size_t &, std::size_t &,
                                    std::size_t &, int &) noexcept;
void adaptive_set_attempt_observer_raw(TimeControlAttemptObserverRaw,
                                       void *) noexcept;
void adaptive_exercise_trusted_tail_attempt_observer_raw(
    const FlowState &, const MomentumTimeStencil &);
void adaptive_set_active_raw(Bdf2RetryController &, bool) noexcept;
void adaptive_stability_rates_raw(Bdf2RetryController &, const FlowState &,
                                  double, bool, double, double &, double &);
void adaptive_stability_reduction_observation_raw(
    std::uint64_t &, std::uint64_t &) noexcept;

} // namespace hundun::flow::detail

namespace hundun::flow::test {

enum class TimeControlPostReturnMutation : std::uint8_t {
  none,
  attempted_dt,
  reason,
  momentum_x_iterations,
  momentum_y_iterations,
  momentum_z_iterations,
  pressure_one_iterations,
  pressure_two_iterations
};

enum class TimeControlPreflightFault : std::uint8_t {
  none,
  layout,
  capability,
  preparation
};

struct TimeControlTrustedTailObservation final {
  std::uint64_t allocation_attempts{};
  std::uint64_t controller_collectives{};
  std::uint64_t field_state_traversals{};
  std::uint64_t callbacks_or_sinks{};
};

struct TimeControlRawFaultObservation final {
  std::size_t raw_operations{};
  std::size_t local_origins{};
  std::size_t fault_ordinal{};
  int requested_rank{-1};
};

using TimeControlAttemptObserver =
    void (*)(void *, const MomentumTimeStencil &, const FlowState &);

class AdaptiveTimeControlTestAccess final {
public:
  static TimeAdvanceReport report() noexcept {
    return TimeAdvanceReport(TimeAdvanceReport::ConstantReportTag{});
  }
  static std::uint64_t seal(const config::FlowTimeConfig &config,
                            config::DensityModel model,
                            const TimeControlState &state) noexcept {
    return detail::TimeControlStateCodec::seal(config, model, state);
  }
  static bool valid_config(const config::FlowTimeConfig &config) noexcept {
    return detail::time_control_config_valid(config);
  }
  static bool valid_model(config::DensityModel model) noexcept {
    return detail::density_model_valid(model);
  }
  static bool valid_control(const linear::SolveControl &control) noexcept {
    return detail::solve_control_valid(control);
  }
  static void set_post_return_mutation(TimeControlPostReturnMutation,
                                       int rank) noexcept;
  static void set_preflight_fault(TimeControlPreflightFault,
                                  int rank) noexcept;
  static void set_recoverable_failures(std::uint32_t count,
                                       int rank) noexcept;
  static void
  set_recoverable_failure_reason(StepFailureReason reason) noexcept;
  static void
  set_post_return_iteration_value(std::uint64_t value) noexcept;
  static void reset_faults() noexcept;
  static std::uint8_t
  preflight_category(const TimeAdvanceReport &) noexcept;
  static TimeControlTrustedTailObservation
  trusted_tail_observation() noexcept;
  static void set_raw_fault(std::size_t ordinal, int rank) noexcept;
  static std::size_t raw_operation_count() noexcept;
  static TimeControlRawFaultObservation raw_fault_observation() noexcept;
  static void set_attempt_observer(TimeControlAttemptObserver,
                                   void *) noexcept;
  static void exercise_trusted_tail_attempt_observer(
      const FlowState &, const MomentumTimeStencil &);
  static void set_active(Bdf2RetryController &, bool) noexcept;
  static std::array<double, 2>
  stability_rates(Bdf2RetryController &, const FlowState &,
                  double density_constant, bool variable_density,
                  double gamma);
  static std::array<std::uint64_t, 2>
  stability_reduction_observation() noexcept;
};

} // namespace hundun::flow::test

inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    set_post_return_mutation(TimeControlPostReturnMutation mutation,
                             int rank) noexcept {
  detail::adaptive_set_post_return_mutation_raw(
      static_cast<std::uint8_t>(mutation), rank);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    set_preflight_fault(TimeControlPreflightFault fault, int rank) noexcept {
  detail::adaptive_set_preflight_fault_raw(static_cast<std::uint8_t>(fault),
                                           rank);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    set_recoverable_failures(std::uint32_t count, int rank) noexcept {
  detail::adaptive_set_recoverable_failures_raw(count, rank);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    set_recoverable_failure_reason(StepFailureReason reason) noexcept {
  detail::adaptive_set_recoverable_failure_reason_raw(reason);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    set_post_return_iteration_value(std::uint64_t value) noexcept {
  detail::adaptive_set_post_return_iteration_value_raw(value);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    reset_faults() noexcept {
  detail::adaptive_reset_faults_raw();
}
inline std::uint8_t hundun::flow::test::AdaptiveTimeControlTestAccess::
    preflight_category(const TimeAdvanceReport &report) noexcept {
  return detail::adaptive_preflight_category_raw(report);
}
inline hundun::flow::test::TimeControlTrustedTailObservation
hundun::flow::test::AdaptiveTimeControlTestAccess::trusted_tail_observation()
    noexcept {
  TimeControlTrustedTailObservation result;
  detail::adaptive_trusted_tail_observation_raw(
      result.allocation_attempts, result.controller_collectives,
      result.field_state_traversals, result.callbacks_or_sinks);
  return result;
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::set_raw_fault(
    std::size_t ordinal, int rank) noexcept {
  detail::adaptive_set_raw_fault(ordinal, rank);
}
inline std::size_t hundun::flow::test::AdaptiveTimeControlTestAccess::
    raw_operation_count() noexcept {
  return detail::adaptive_raw_operation_count();
}
inline hundun::flow::test::TimeControlRawFaultObservation
hundun::flow::test::AdaptiveTimeControlTestAccess::raw_fault_observation()
    noexcept {
  TimeControlRawFaultObservation result;
  detail::adaptive_raw_fault_observation(
      result.raw_operations, result.local_origins, result.fault_ordinal,
      result.requested_rank);
  return result;
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    set_attempt_observer(TimeControlAttemptObserver observer,
                         void *context) noexcept {
  detail::adaptive_set_attempt_observer_raw(observer, context);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::
    exercise_trusted_tail_attempt_observer(
        const FlowState &state, const MomentumTimeStencil &stencil) {
  detail::adaptive_exercise_trusted_tail_attempt_observer_raw(state, stencil);
}
inline void hundun::flow::test::AdaptiveTimeControlTestAccess::set_active(
    Bdf2RetryController &controller, bool active) noexcept {
  detail::adaptive_set_active_raw(controller, active);
}
inline std::array<double, 2>
hundun::flow::test::AdaptiveTimeControlTestAccess::stability_rates(
    Bdf2RetryController &controller, const FlowState &state,
    double density_constant, bool variable_density, double gamma) {
  std::array<double, 2> result{};
  detail::adaptive_stability_rates_raw(
      controller, state, density_constant, variable_density, gamma,
      result[0], result[1]);
  return result;
}
inline std::array<std::uint64_t, 2>
hundun::flow::test::AdaptiveTimeControlTestAccess::
    stability_reduction_observation() noexcept {
  std::array<std::uint64_t, 2> result{};
  detail::adaptive_stability_reduction_observation_raw(result[0], result[1]);
  return result;
}
