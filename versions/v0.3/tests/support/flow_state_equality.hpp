// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_state.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/flow_ideal_gas_closure.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <optional>
#include <vector>

namespace hundun::test {

inline std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t bits{};
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline bool
fp64_vector_bitwise_equal(const std::vector<double> &left,
                          const std::vector<double> &right) noexcept {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (fp64_bits(left[index]) != fp64_bits(right[index]))
      return false;
  }
  return true;
}

inline bool
flow_layer_values_bitwise_equal(const flow::FlowLayerValues &left,
                                const flow::FlowLayerValues &right) noexcept {
  if (!fp64_vector_bitwise_equal(left.density, right.density) ||
      !fp64_vector_bitwise_equal(left.velocity, right.velocity) ||
      !fp64_vector_bitwise_equal(left.mechanical_pressure,
                                 right.mechanical_pressure) ||
      !fp64_vector_bitwise_equal(left.face_velocity, right.face_velocity) ||
      !fp64_vector_bitwise_equal(left.face_mass_flux, right.face_mass_flux) ||
      left.transported_cell_fields.size() !=
          right.transported_cell_fields.size())
    return false;
  for (std::size_t field = 0; field < left.transported_cell_fields.size();
       ++field) {
    if (!fp64_vector_bitwise_equal(left.transported_cell_fields[field],
                                   right.transported_cell_fields[field]))
      return false;
  }
  return true;
}

inline bool accepted_step_metadata_bitwise_equal(
    const flow::AcceptedStepMetadata &left,
    const flow::AcceptedStepMetadata &right) noexcept {
  return left.step == right.step && left.order == right.order &&
         fp64_bits(left.time_s) == fp64_bits(right.time_s) &&
         fp64_bits(left.dt_s) == fp64_bits(right.dt_s) &&
         fp64_bits(left.previous_dt_s) == fp64_bits(right.previous_dt_s);
}

struct IdealGasStateSnapshot final {
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  flow::FlowLayerValues trial;
  flow::AcceptedStepMetadata metadata;
  flow::IdealGasClosureState closure;
};

inline bool ideal_gas_closure_state_bitwise_equal(
    const flow::IdealGasClosureState &left,
    const flow::IdealGasClosureState &right) noexcept {
  return left.mode == right.mode && left.revision == right.revision &&
         fp64_bits(left.thermodynamic_pressure_pa) ==
             fp64_bits(right.thermodynamic_pressure_pa) &&
         left.target_mass_kg.has_value() == right.target_mass_kg.has_value() &&
         (!left.target_mass_kg ||
          fp64_bits(*left.target_mass_kg) == fp64_bits(*right.target_mass_kg));
}

inline bool ideal_gas_state_bitwise_equal(
    const IdealGasStateSnapshot &left,
    const IdealGasStateSnapshot &right) noexcept {
  return flow_layer_values_bitwise_equal(left.history, right.history) &&
         flow_layer_values_bitwise_equal(left.committed, right.committed) &&
         flow_layer_values_bitwise_equal(left.trial, right.trial) &&
         accepted_step_metadata_bitwise_equal(left.metadata, right.metadata) &&
         ideal_gas_closure_state_bitwise_equal(left.closure, right.closure);
}

inline bool ideal_gas_state_equality_oracle_is_mutation_sensitive() {
  IdealGasStateSnapshot baseline;
  baseline.history.density = {1.0};
  baseline.committed.density = {1.0};
  baseline.trial.density = {1.0};
  baseline.committed.transported_cell_fields = {{300000.0}};
  baseline.closure = {flow::IdealGasPressureMode::closed_dynamic, 101325.0,
                      1.0, 3U};
  const auto exact = baseline;
  auto ordinary = baseline;
  ordinary.committed.density.front() = 2.0;
  auto nested = baseline;
  nested.committed.transported_cell_fields.front().front() = 300001.0;
  auto pressure = baseline;
  pressure.closure.thermodynamic_pressure_pa = -101325.0;
  auto revision = baseline;
  ++revision.closure.revision;
  return ideal_gas_state_bitwise_equal(baseline, exact) &&
         !ideal_gas_state_bitwise_equal(baseline, ordinary) &&
         !ideal_gas_state_bitwise_equal(baseline, nested) &&
         !ideal_gas_state_bitwise_equal(baseline, pressure) &&
         !ideal_gas_state_bitwise_equal(baseline, revision);
}

inline bool flow_state_equality_oracle_is_mutation_sensitive() {
  flow::FlowLayerValues baseline;
  baseline.density = {0.0, 1.0};
  baseline.velocity = {2.0, 3.0, 4.0};
  baseline.mechanical_pressure = {5.0};
  baseline.face_velocity = {6.0, 7.0, 8.0};
  baseline.face_mass_flux = {9.0};
  baseline.transported_cell_fields = {{0.0, 10.0}, {11.0}};
  const auto exact = baseline;
  auto ordinary = baseline;
  ordinary.density.front() = -0.0;
  auto nested = baseline;
  nested.transported_cell_fields.front().front() = -0.0;
  return flow_layer_values_bitwise_equal(baseline, exact) &&
         !flow_layer_values_bitwise_equal(baseline, ordinary) &&
         !flow_layer_values_bitwise_equal(baseline, nested);
}

inline bool time_control_state_bitwise_equal(
    const flow::TimeControlState &left,
    const flow::TimeControlState &right) noexcept {
  return left.schema_version == right.schema_version &&
         left.accepted_step == right.accepted_step &&
         fp64_bits(left.proposed_next_dt_s) ==
             fp64_bits(right.proposed_next_dt_s) &&
         fp64_bits(left.last_accepted_dt_s) ==
             fp64_bits(right.last_accepted_dt_s) &&
         left.last_accepted_order == right.last_accepted_order &&
         left.history_ready == right.history_ready &&
         left.last_all_linear_solves_within_half_limit ==
             right.last_all_linear_solves_within_half_limit &&
         fp64_bits(left.last_convective_rate_per_s) ==
             fp64_bits(right.last_convective_rate_per_s) &&
         fp64_bits(left.last_diffusive_rate_per_s) ==
             fp64_bits(right.last_diffusive_rate_per_s) &&
         left.last_stability_metrics_available ==
             right.last_stability_metrics_available &&
         left.last_retry_count == right.last_retry_count &&
         left.revision == right.revision &&
         left.state_seal == right.state_seal;
}

inline bool time_control_state_equality_oracle_is_mutation_sensitive() {
  flow::TimeControlState baseline;
  baseline.accepted_step = 2U;
  baseline.proposed_next_dt_s = 0.1;
  baseline.last_accepted_dt_s = 0.05;
  baseline.last_accepted_order = flow::MomentumTimeOrder::bdf2;
  baseline.history_ready = true;
  baseline.last_all_linear_solves_within_half_limit = true;
  baseline.last_convective_rate_per_s = 2.0;
  baseline.last_diffusive_rate_per_s = 3.0;
  baseline.last_stability_metrics_available = true;
  baseline.last_retry_count = 2U;
  baseline.revision = 2U;
  baseline.state_seal = 4U;
  const auto exact = baseline;
  const auto changed = [&](auto mutate) {
    auto candidate = baseline;
    mutate(candidate);
    return !time_control_state_bitwise_equal(baseline, candidate);
  };
  return time_control_state_bitwise_equal(baseline, exact) &&
         changed([](auto &s) { ++s.schema_version; }) &&
         changed([](auto &s) { ++s.accepted_step; }) &&
         changed([](auto &s) { s.proposed_next_dt_s = -0.1; }) &&
         changed([](auto &s) { s.last_accepted_dt_s = -0.05; }) &&
         changed([](auto &s) {
           s.last_accepted_order = flow::MomentumTimeOrder::backward_euler;
         }) &&
         changed([](auto &s) { s.history_ready = false; }) &&
         changed([](auto &s) {
           s.last_all_linear_solves_within_half_limit = false;
         }) &&
         changed([](auto &s) { s.last_convective_rate_per_s = -2.0; }) &&
         changed([](auto &s) { s.last_diffusive_rate_per_s = -3.0; }) &&
         changed([](auto &s) {
           s.last_stability_metrics_available = false;
         }) &&
         changed([](auto &s) { ++s.last_retry_count; }) &&
         changed([](auto &s) { ++s.revision; }) &&
         changed([](auto &s) { ++s.state_seal; });
}

struct CheckpointV2StateSnapshot final {
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  flow::FlowLayerValues trial;
  flow::AcceptedStepMetadata metadata;
  flow::TimeControlState controller;
  std::optional<flow::IdealGasClosureState> closure;
};

inline bool checkpoint_v2_state_bitwise_equal(
    const CheckpointV2StateSnapshot &left,
    const CheckpointV2StateSnapshot &right) noexcept {
  return flow_layer_values_bitwise_equal(left.history, right.history) &&
         flow_layer_values_bitwise_equal(left.committed, right.committed) &&
         flow_layer_values_bitwise_equal(left.trial, right.trial) &&
         accepted_step_metadata_bitwise_equal(left.metadata, right.metadata) &&
         time_control_state_bitwise_equal(left.controller, right.controller) &&
         left.closure.has_value() == right.closure.has_value() &&
         (!left.closure ||
          ideal_gas_closure_state_bitwise_equal(*left.closure,
                                                *right.closure));
}

inline bool checkpoint_v2_state_equality_oracle_is_mutation_sensitive() {
  CheckpointV2StateSnapshot baseline;
  baseline.history.density = {1.0};
  baseline.committed.density = {2.0};
  baseline.trial.density = {2.0};
  baseline.committed.transported_cell_fields = {{3.0}, {4.0}};
  baseline.controller.proposed_next_dt_s = 0.1;
  baseline.closure = flow::IdealGasClosureState{
      flow::IdealGasPressureMode::closed_dynamic, 101325.0, 1.0, 3U};
  const auto exact = baseline;
  auto ordinary = baseline;
  ordinary.committed.density.front() = -2.0;
  auto nested = baseline;
  nested.committed.transported_cell_fields.front().front() = -3.0;
  auto metadata = baseline;
  metadata.metadata.time_s = -0.0;
  auto controller = baseline;
  controller.controller.proposed_next_dt_s = -0.1;
  auto closure = baseline;
  ++closure.closure->revision;
  return checkpoint_v2_state_bitwise_equal(baseline, exact) &&
         !checkpoint_v2_state_bitwise_equal(baseline, ordinary) &&
         !checkpoint_v2_state_bitwise_equal(baseline, nested) &&
         !checkpoint_v2_state_bitwise_equal(baseline, metadata) &&
         !checkpoint_v2_state_bitwise_equal(baseline, controller) &&
         !checkpoint_v2_state_bitwise_equal(baseline, closure);
}

struct AdaptiveFlowStateSnapshot final {
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  flow::FlowLayerValues trial;
  flow::AcceptedStepMetadata metadata;
  flow::TimeControlState controller;
  flow::IdealGasClosureState closure;
  std::uint64_t committed_allocation_identity{};
  std::uint64_t history_allocation_identity{};
  std::uint64_t trial_allocation_identity{};
  std::uint64_t committed_generation{};
  std::uint64_t history_generation{};
  std::uint64_t trial_generation{};
  std::uint64_t attempt_identity{};
  std::uint64_t diagnostic_identity{};
  std::uint64_t accepted_cache_identity{};
  std::size_t accepted_cache_capacity{};
  std::uint64_t accepted_cache_revision{};
  std::array<std::uint64_t, 4> operation_counters{};
};

inline bool adaptive_flow_state_bitwise_equal(
    const AdaptiveFlowStateSnapshot &left,
    const AdaptiveFlowStateSnapshot &right) noexcept {
  return flow_layer_values_bitwise_equal(left.history, right.history) &&
         flow_layer_values_bitwise_equal(left.committed, right.committed) &&
         flow_layer_values_bitwise_equal(left.trial, right.trial) &&
         accepted_step_metadata_bitwise_equal(left.metadata, right.metadata) &&
         time_control_state_bitwise_equal(left.controller, right.controller) &&
         ideal_gas_closure_state_bitwise_equal(left.closure, right.closure) &&
         left.committed_allocation_identity ==
             right.committed_allocation_identity &&
         left.history_allocation_identity == right.history_allocation_identity &&
         left.trial_allocation_identity == right.trial_allocation_identity &&
         left.committed_generation == right.committed_generation &&
         left.history_generation == right.history_generation &&
         left.trial_generation == right.trial_generation &&
         left.attempt_identity == right.attempt_identity &&
         left.diagnostic_identity == right.diagnostic_identity &&
         left.accepted_cache_identity == right.accepted_cache_identity &&
         left.accepted_cache_capacity == right.accepted_cache_capacity &&
         left.accepted_cache_revision == right.accepted_cache_revision &&
         left.operation_counters == right.operation_counters;
}

inline bool adaptive_flow_state_equality_oracle_is_mutation_sensitive() {
  AdaptiveFlowStateSnapshot baseline;
  const auto populate_layer = [](flow::FlowLayerValues &layer) {
    layer.density = {1.0};
    layer.velocity = {2.0, 3.0, 4.0};
    layer.mechanical_pressure = {5.0};
    layer.face_velocity = {6.0, 7.0, 8.0};
    layer.face_mass_flux = {9.0};
    layer.transported_cell_fields = {{10.0}, {11.0}};
  };
  populate_layer(baseline.history);
  populate_layer(baseline.committed);
  populate_layer(baseline.trial);
  baseline.metadata = {2U, 0.1, 0.05, 0.05,
                       flow::MomentumTimeOrder::bdf2};
  baseline.controller.proposed_next_dt_s = 0.1;
  baseline.controller.last_accepted_dt_s = 0.05;
  baseline.controller.accepted_step = 2U;
  baseline.controller.last_accepted_order = flow::MomentumTimeOrder::bdf2;
  baseline.controller.history_ready = true;
  baseline.controller.last_all_linear_solves_within_half_limit = true;
  baseline.controller.last_convective_rate_per_s = 2.0;
  baseline.controller.last_diffusive_rate_per_s = 3.0;
  baseline.controller.last_stability_metrics_available = true;
  baseline.controller.last_retry_count = 2U;
  baseline.controller.revision = 2U;
  baseline.controller.state_seal = 1U;
  baseline.closure = {flow::IdealGasPressureMode::closed_dynamic, 101325.0,
                      1.0, 2U};
  baseline.committed_allocation_identity = 10U;
  baseline.history_allocation_identity = 11U;
  baseline.trial_allocation_identity = 12U;
  baseline.committed_generation = 20U;
  baseline.history_generation = 21U;
  baseline.trial_generation = 22U;
  baseline.attempt_identity = 30U;
  baseline.diagnostic_identity = 31U;
  baseline.accepted_cache_identity = 40U;
  baseline.accepted_cache_capacity = 41U;
  baseline.accepted_cache_revision = 42U;
  baseline.operation_counters = {50U, 51U, 52U, 53U};
  const auto exact = baseline;
  const auto changed = [&](auto mutate) {
    auto candidate = baseline;
    mutate(candidate);
    return !adaptive_flow_state_bitwise_equal(baseline, candidate);
  };
  return adaptive_flow_state_bitwise_equal(baseline, exact) &&
         changed([](auto &s) { s.history.density[0] = -1.0; }) &&
         changed([](auto &s) { s.history.velocity[0] = -2.0; }) &&
         changed([](auto &s) {
           s.history.mechanical_pressure[0] = -5.0;
         }) &&
         changed([](auto &s) { s.history.face_velocity[0] = -6.0; }) &&
         changed([](auto &s) { s.history.face_mass_flux[0] = -9.0; }) &&
         changed([](auto &s) { s.committed.density[0] = -1.0; }) &&
         changed([](auto &s) { s.committed.velocity[0] = -2.0; }) &&
         changed([](auto &s) {
           s.committed.mechanical_pressure[0] = -5.0;
         }) &&
         changed([](auto &s) { s.committed.face_velocity[0] = -6.0; }) &&
         changed([](auto &s) {
           s.committed.face_mass_flux[0] = -9.0;
         }) &&
         changed([](auto &s) { s.trial.density[0] = -1.0; }) &&
         changed([](auto &s) { s.trial.velocity[0] = -2.0; }) &&
         changed([](auto &s) {
           s.trial.mechanical_pressure[0] = -5.0;
         }) &&
         changed([](auto &s) { s.trial.face_velocity[0] = -6.0; }) &&
         changed([](auto &s) { s.trial.face_mass_flux[0] = -9.0; }) &&
         changed([](auto &s) {
           s.history.transported_cell_fields[1][0] = -3.0;
         }) &&
         changed([](auto &s) {
           s.committed.transported_cell_fields[1][0] = -3.0;
         }) &&
         changed([](auto &s) {
           s.trial.transported_cell_fields[1][0] = -3.0;
         }) &&
         changed([](auto &s) {
           s.committed.transported_cell_fields[0].push_back(12.0);
         }) &&
         changed([](auto &s) {
           s.history.transported_cell_fields.push_back({4.0});
         }) &&
         changed([](auto &s) {
           s.committed.transported_cell_fields.push_back({4.0});
         }) &&
         changed([](auto &s) {
           s.trial.transported_cell_fields.push_back({4.0});
         }) &&
         changed([](auto &s) { ++s.metadata.step; }) &&
         changed([](auto &s) { s.metadata.time_s = 1.0; }) &&
         changed([](auto &s) { s.metadata.dt_s = 0.025; }) &&
         changed([](auto &s) { s.metadata.previous_dt_s = 0.025; }) &&
         changed([](auto &s) {
           s.metadata.order = flow::MomentumTimeOrder::backward_euler;
         }) &&
         changed([](auto &s) { ++s.controller.schema_version; }) &&
         changed([](auto &s) { ++s.controller.accepted_step; }) &&
         changed([](auto &s) { s.controller.proposed_next_dt_s = -0.1; }) &&
         changed([](auto &s) { s.controller.last_accepted_dt_s = -0.05; }) &&
         changed([](auto &s) {
           s.controller.last_accepted_order =
               flow::MomentumTimeOrder::backward_euler;
         }) &&
         changed([](auto &s) { s.controller.history_ready = false; }) &&
         changed([](auto &s) {
           s.controller.last_all_linear_solves_within_half_limit = false;
         }) &&
         changed([](auto &s) {
           s.controller.last_convective_rate_per_s = 1.0;
         }) &&
         changed([](auto &s) {
           s.controller.last_diffusive_rate_per_s = 1.0;
         }) &&
         changed([](auto &s) {
           s.controller.last_stability_metrics_available = false;
         }) &&
         changed([](auto &s) { ++s.controller.last_retry_count; }) &&
         changed([](auto &s) { ++s.controller.revision; }) &&
         changed([](auto &s) { ++s.controller.state_seal; }) &&
         changed([](auto &s) {
           s.closure.mode = flow::IdealGasPressureMode::open_fixed;
         }) &&
         changed([](auto &s) {
           s.closure.thermodynamic_pressure_pa = 101326.0;
         }) &&
         changed([](auto &s) { s.closure.target_mass_kg.reset(); }) &&
         changed([](auto &s) { s.closure.target_mass_kg = 2.0; }) &&
         changed([](auto &s) { ++s.closure.revision; }) &&
         changed([](auto &s) { ++s.committed_allocation_identity; }) &&
         changed([](auto &s) { ++s.history_allocation_identity; }) &&
         changed([](auto &s) { ++s.trial_allocation_identity; }) &&
         changed([](auto &s) { ++s.committed_generation; }) &&
         changed([](auto &s) { ++s.history_generation; }) &&
         changed([](auto &s) { ++s.trial_generation; }) &&
         changed([](auto &s) { ++s.attempt_identity; }) &&
         changed([](auto &s) { ++s.diagnostic_identity; }) &&
         changed([](auto &s) { ++s.accepted_cache_identity; }) &&
         changed([](auto &s) { ++s.accepted_cache_capacity; }) &&
         changed([](auto &s) { ++s.accepted_cache_revision; }) &&
         changed([](auto &s) { ++s.operation_counters[0]; }) &&
         changed([](auto &s) { ++s.operation_counters[1]; }) &&
         changed([](auto &s) { ++s.operation_counters[2]; }) &&
         changed([](auto &s) { ++s.operation_counters[3]; });
}

} // namespace hundun::test
