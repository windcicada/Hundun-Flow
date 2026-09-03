// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_adaptive_time_control_test_access.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/rt_error.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <cstdint>
#include <cstring>
#include <string_view>
#include <array>
#include <cmath>
#include <limits>

namespace {

void add_byte(std::uint64_t &hash, std::uint8_t value) {
  hash = (hash ^ value) * 1099511628211ULL;
}
template <class T> void add_le(std::uint64_t &hash, T value) {
  for (std::size_t i = 0; i < sizeof(T); ++i)
    add_byte(hash, static_cast<std::uint8_t>(value >> (8U * i)));
}
void add_double(std::uint64_t &hash, double value) {
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  add_le(hash, bits);
}

std::uint64_t independent_seal(
    const hundun::config::FlowTimeConfig &config,
    hundun::config::DensityModel model,
    const hundun::flow::TimeControlState &state) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (char value : std::string_view("hundun-time-control-state-seal-v1"))
    add_byte(hash, static_cast<std::uint8_t>(value));
  add_le(hash, static_cast<std::uint32_t>(config.mode));
  add_double(hash, config.initial_dt_s);
  add_double(hash, config.min_dt_s);
  add_double(hash, config.max_dt_s);
  add_double(hash, config.cfl_target);
  add_double(hash, config.diffusion_number_target);
  add_double(hash, config.growth_factor);
  add_double(hash, config.retry_factor);
  add_le(hash, static_cast<std::uint32_t>(config.max_retries));
  add_le(hash, static_cast<std::uint32_t>(model));
  add_le(hash, state.schema_version);
  add_le(hash, state.accepted_step);
  add_double(hash, state.proposed_next_dt_s);
  add_double(hash, state.last_accepted_dt_s);
  add_byte(hash, static_cast<std::uint8_t>(state.last_accepted_order));
  add_byte(hash, state.history_ready ? 1U : 0U);
  add_byte(hash, state.last_all_linear_solves_within_half_limit ? 1U : 0U);
  add_double(hash, state.last_convective_rate_per_s);
  add_double(hash, state.last_diffusive_rate_per_s);
  add_byte(hash, state.last_stability_metrics_available ? 1U : 0U);
  add_le(hash, state.last_retry_count);
  add_le(hash, state.revision);
  return hash;
}

void run() {
  HUNDUN_CHECK(
      hundun::test::time_control_state_equality_oracle_is_mutation_sensitive());
  HUNDUN_CHECK(
      hundun::test::adaptive_flow_state_equality_oracle_is_mutation_sensitive());
  hundun::config::FlowTimeConfig config{
      hundun::config::TimeMode::adaptive, 10, 0.1, 0.0125, 0.2,
      0.5, 0.25, 1.25, 0.5, 8};
  hundun::flow::TimeControlState state;
  state.accepted_step = (std::uint64_t{1} << 53U) + 1U;
  state.proposed_next_dt_s = 0.125;
  state.last_accepted_dt_s = 0.1;
  state.last_accepted_order = hundun::flow::MomentumTimeOrder::bdf2;
  state.history_ready = true;
  state.last_all_linear_solves_within_half_limit = true;
  state.last_convective_rate_per_s = 2.0;
  state.last_diffusive_rate_per_s = 1.0;
  state.last_stability_metrics_available = true;
  state.last_retry_count = 3U;
  state.revision = state.accepted_step;
  const auto expected = independent_seal(
      config, hundun::config::DensityModel::material, state);
  HUNDUN_CHECK(expected == 6963059103000454808ULL);
  HUNDUN_CHECK(hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
                   config, hundun::config::DensityModel::material, state) ==
               expected);
  auto changed = state;
  changed.last_retry_count = 4U;
  HUNDUN_CHECK(hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
                   config, hundun::config::DensityModel::material, changed) !=
               expected);
  const auto config_changes_seal = [&](auto mutate) {
    auto candidate = config;
    mutate(candidate);
    return hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
               candidate, hundun::config::DensityModel::material, state) !=
           expected;
  };
  HUNDUN_CHECK(config_changes_seal(
      [](auto &value) { value.mode = hundun::config::TimeMode::fixed; }));
  {
    auto different_horizon = config;
    ++different_horizon.steps;
    HUNDUN_CHECK(
        hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
            different_horizon, hundun::config::DensityModel::material,
            state) == expected);
  }
  HUNDUN_CHECK(config_changes_seal(
      [](auto &value) { value.initial_dt_s = 0.11; }));
  HUNDUN_CHECK(
      config_changes_seal([](auto &value) { value.min_dt_s = 0.01; }));
  HUNDUN_CHECK(
      config_changes_seal([](auto &value) { value.max_dt_s = 0.19; }));
  HUNDUN_CHECK(
      config_changes_seal([](auto &value) { value.cfl_target = 0.4; }));
  HUNDUN_CHECK(config_changes_seal(
      [](auto &value) { value.diffusion_number_target = 0.2; }));
  HUNDUN_CHECK(
      config_changes_seal([](auto &value) { value.growth_factor = 1.2; }));
  HUNDUN_CHECK(
      config_changes_seal([](auto &value) { value.retry_factor = 0.4; }));
  HUNDUN_CHECK(
      config_changes_seal([](auto &value) { --value.max_retries; }));
  HUNDUN_CHECK(hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
                   config, hundun::config::DensityModel::ideal_gas, state) !=
               expected);

  auto report =
      hundun::flow::test::AdaptiveTimeControlTestAccess::report();
  HUNDUN_CHECK(report.attempt_count() == 0U);
  auto moved = std::move(report);
  HUNDUN_CHECK(report.attempt_count() == 0U);
  HUNDUN_CHECK(report.reason() ==
               hundun::flow::StepFailureReason::invalid_input);
  HUNDUN_CHECK(!report.final_attempt_available());
  bool rejected = false;
  try {
    static_cast<void>(report.final_attempt());
  } catch (const hundun::runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(moved.attempt_count() == 0U);

  auto bad_config = config;
  bad_config.mode = static_cast<hundun::config::TimeMode>(255);
  HUNDUN_CHECK(
      !hundun::flow::test::AdaptiveTimeControlTestAccess::valid_config(
          bad_config));
  auto bad_model = static_cast<hundun::config::DensityModel>(255);
  HUNDUN_CHECK(
      !hundun::flow::test::AdaptiveTimeControlTestAccess::valid_model(
          bad_model));
  hundun::linear::SolveControl bad_control{};
  bad_control.residual_recompute_interval = 0U;
  HUNDUN_CHECK(
      !hundun::flow::test::AdaptiveTimeControlTestAccess::valid_control(
          bad_control));

  for (const double ratio : std::array<double, 4>{0.5, 1.0, 1.25, 2.0}) {
    const double previous = 0.08;
    const double current = ratio * previous;
    const auto stencil = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::bdf2, current, previous);
    HUNDUN_CHECK(stencil.order == hundun::flow::MomentumTimeOrder::bdf2);
    HUNDUN_CHECK(stencil.alpha0 == (1.0 + 2.0 * ratio) / (1.0 + ratio));
    HUNDUN_CHECK(stencil.alpha1 == -(1.0 + ratio));
    HUNDUN_CHECK(stencil.alpha2 ==
                 (ratio * ratio) / (1.0 + ratio));
  }
  const auto be = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
  HUNDUN_CHECK(be.alpha0 == 1.0);
  HUNDUN_CHECK(be.alpha1 == -1.0);
  HUNDUN_CHECK(be.alpha2 == 0.0);

  const auto rejects_stencil = [](auto order, double current,
                                  double previous) {
    try {
      static_cast<void>(hundun::flow::make_momentum_time_stencil(
          order, current, previous));
    } catch (const hundun::runtime::Error &) {
      return true;
    }
    return false;
  };
  const auto infinity = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  HUNDUN_CHECK(rejects_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, 0.0, 0.0));
  HUNDUN_CHECK(rejects_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, -0.01, 0.0));
  HUNDUN_CHECK(rejects_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, infinity, 0.0));
  HUNDUN_CHECK(rejects_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, nan, 0.0));
  HUNDUN_CHECK(rejects_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                               0.01, 0.0));
  HUNDUN_CHECK(rejects_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                               0.01, -0.01));
  HUNDUN_CHECK(rejects_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                               0.01, infinity));
  HUNDUN_CHECK(rejects_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                               0.01, nan));
  HUNDUN_CHECK(rejects_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                               0.049, 0.1));
  HUNDUN_CHECK(rejects_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                               0.201, 0.1));
  HUNDUN_CHECK(rejects_stencil(
      static_cast<hundun::flow::MomentumTimeOrder>(255), 0.01, 0.01));
}

} // namespace

int main() { return hundun::test::run(run); }
