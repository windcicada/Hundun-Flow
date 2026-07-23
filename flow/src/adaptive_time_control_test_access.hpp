// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "adaptive_time_control_detail.hpp"

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
  static void reset_faults() noexcept;
  static std::uint8_t
  preflight_category(const TimeAdvanceReport &) noexcept;
  static std::uint64_t post_commit_observation_count() noexcept;
};

} // namespace hundun::flow::test
