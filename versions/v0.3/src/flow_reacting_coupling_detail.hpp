// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "flow_reacting_transaction_detail.hpp"

namespace hundun::flow {

struct ReactingCouplingSchedule final {
  std::array<double, 2> chemistry_duration_fractions{0.5, 0.5};
  std::uint32_t chemistry_call_count{2U};
  std::uint32_t scalar_transport_count{1U};
  std::uint32_t pressure_corrector_count{2U};
  bool piso2_consumes_post_chemistry2{true};
  bool requires_predictor_to_final_delta_flux{true};

  static constexpr ReactingCouplingSchedule second_order() noexcept {
    return {};
  }
};

struct ReactingStepReport final {
  std::uint32_t chemistry_call_count{};
  std::uint32_t scalar_transport_count{};
  std::uint32_t pressure_corrector_count{};
  std::uint64_t post_chemistry2_epoch{};
  std::uint64_t piso2_consumed_epoch{};
  bool predictor_to_final_delta_flux_applied{};
  bool accepted{};
  std::string failure_stage;
  double p0_before_pa{};
  double p0_after_pa{};
  std::vector<double> integrated_species_delta_kg_per_m3;
  std::vector<double> integrated_enthalpy_delta_j_per_m3;
};

struct ReactingOperatorUpdate final {
  bool ok{true};
  std::string message;
  std::vector<double> species_delta_kg_per_m3;
  std::vector<double> enthalpy_delta_j_per_m3;
};

struct OpenReactingStepOperators final {
  std::function<ReactingOperatorUpdate(
      std::uint32_t, double, double, const detail::ReactingLayerState &)>
      chemistry;
  std::function<ReactingOperatorUpdate(
      double, double, const detail::ReactingLayerState &)>
      scalar_transport;
  std::function<bool(std::uint32_t, double, std::uint64_t,
                     const detail::ReactingLayerState &,
                     const std::vector<double> &, const std::vector<double> &,
                     std::string &)>
      pressure_corrector;
  std::function<runtime::CollectiveStatus(bool, const std::string &)>
      collective_validation;
};

ReactingStepReport attempt_open_reacting_step(
    detail::ReactingAttemptState &, double start_time_s, double duration_s,
    const OpenReactingStepOperators &);

enum class ReactingCouplingEventKind : std::uint8_t {
  chemistry,
  scalar_transport,
  piso
};

struct ReactingCouplingEvent final {
  ReactingCouplingEventKind kind{};
  double start_time_s{};
  double duration_s{};
  std::uint64_t state_epoch{};
};

inline bool validate_reacting_coupling_trace(
    const ReactingCouplingSchedule &schedule,
    const std::vector<ReactingCouplingEvent> &trace, double start_time_s,
    double duration_s) noexcept {
  if (!std::isfinite(start_time_s) || !std::isfinite(duration_s) ||
      duration_s <= 0.0 || trace.size() != 5U ||
      schedule.chemistry_call_count != 2U ||
      schedule.scalar_transport_count != 1U ||
      schedule.pressure_corrector_count != 2U ||
      schedule.chemistry_duration_fractions !=
          std::array<double, 2>{0.5, 0.5} ||
      !schedule.piso2_consumes_post_chemistry2 ||
      !schedule.requires_predictor_to_final_delta_flux) {
    return false;
  }
  const double half = 0.5 * duration_s;
  const double midpoint = start_time_s + half;
  const double end = start_time_s + duration_s;
  const double tolerance =
      1.0e-12 * std::max(1.0, std::abs(end));
  const auto same = [&](double left, double right) {
    return std::abs(left - right) <= tolerance;
  };
  return trace[0].kind == ReactingCouplingEventKind::chemistry &&
         same(trace[0].start_time_s, start_time_s) &&
         same(trace[0].duration_s, half) && trace[0].state_epoch == 1U &&
         trace[1].kind == ReactingCouplingEventKind::scalar_transport &&
         same(trace[1].start_time_s, midpoint) &&
         same(trace[1].duration_s, duration_s) &&
         trace[1].state_epoch == 1U &&
         trace[2].kind == ReactingCouplingEventKind::piso &&
         trace[2].state_epoch == 1U &&
         trace[3].kind == ReactingCouplingEventKind::chemistry &&
         same(trace[3].start_time_s, midpoint) &&
         same(trace[3].duration_s, half) && trace[3].state_epoch == 2U &&
         trace[4].kind == ReactingCouplingEventKind::piso &&
         same(trace[4].start_time_s, end) &&
         trace[4].state_epoch == trace[3].state_epoch;
}

inline bool validate_reacting_step_report(
    const ReactingCouplingSchedule &schedule,
    const ReactingStepReport &report) noexcept {
  return report.accepted && report.p0_before_pa == report.p0_after_pa &&
         report.chemistry_call_count == schedule.chemistry_call_count &&
         report.scalar_transport_count == schedule.scalar_transport_count &&
         report.pressure_corrector_count ==
             schedule.pressure_corrector_count &&
         report.post_chemistry2_epoch != 0U &&
         report.piso2_consumed_epoch == report.post_chemistry2_epoch &&
         report.predictor_to_final_delta_flux_applied ==
             schedule.requires_predictor_to_final_delta_flux;
}

} // namespace hundun::flow
