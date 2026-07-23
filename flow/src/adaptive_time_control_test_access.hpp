// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "adaptive_time_control_detail.hpp"

namespace hundun::flow::test {

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
};

} // namespace hundun::flow::test
