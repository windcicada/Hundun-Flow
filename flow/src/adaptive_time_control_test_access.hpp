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
};

} // namespace hundun::flow::test
