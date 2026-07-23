// SPDX-License-Identifier: Apache-2.0

#include "flow/src/adaptive_time_control_test_access.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/runtime/error.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <cstdint>
#include <cstring>
#include <string_view>

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
  HUNDUN_CHECK(hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
                   config, hundun::config::DensityModel::material, state) ==
               expected);
  auto changed = state;
  changed.last_retry_count = 4U;
  HUNDUN_CHECK(hundun::flow::test::AdaptiveTimeControlTestAccess::seal(
                   config, hundun::config::DensityModel::material, changed) !=
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
}

} // namespace

int main() { return hundun::test::run(run); }
