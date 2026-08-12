// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_state.hpp"
#include "tests/support/flow_state_equality.hpp"

#include <cmath>
#include <limits>

namespace hundun::test {

struct Stage3StateSnapshot final {
  flow::FlowLayerValues history;
  flow::FlowLayerValues committed;
  flow::FlowLayerValues trial;
  flow::AcceptedStepMetadata metadata;
};

inline Stage3StateSnapshot snapshot_stage3_state(const flow::FlowState &state) {
  return {state.snapshot(flow::FlowLayer::history),
          state.snapshot(flow::FlowLayer::committed),
          state.snapshot(flow::FlowLayer::trial), state.metadata()};
}

inline bool
stage3_state_bitwise_equal(const Stage3StateSnapshot &left,
                           const Stage3StateSnapshot &right) noexcept {
  return flow_layer_values_bitwise_equal(left.history, right.history) &&
         flow_layer_values_bitwise_equal(left.committed, right.committed) &&
         flow_layer_values_bitwise_equal(left.trial, right.trial) &&
         accepted_step_metadata_bitwise_equal(left.metadata, right.metadata);
}

inline bool stage3_state_equality_oracle_is_mutation_sensitive() {
  Stage3StateSnapshot baseline;
  const auto populate = [](flow::FlowLayerValues &layer) {
    layer.density = {1.0, 0.0};
    layer.velocity = {2.0, 3.0, 4.0, 0.0, 0.0, 0.0};
    layer.mechanical_pressure = {5.0, 0.0};
    layer.face_velocity = {6.0, 7.0, 8.0};
    layer.face_mass_flux = {9.0};
    layer.transported_cell_fields = {{10.0, 0.0}, {11.0, 0.0}};
  };
  populate(baseline.history);
  populate(baseline.committed);
  populate(baseline.trial);
  baseline.metadata = {2U, 0.1, 0.05, 0.025, flow::MomentumTimeOrder::bdf2};

  const auto exact = baseline;
  const auto changed = [&](auto mutate) {
    auto candidate = baseline;
    mutate(candidate);
    return !stage3_state_bitwise_equal(baseline, candidate);
  };
  return stage3_state_bitwise_equal(baseline, exact) &&
         changed([](auto &state) {
           state.committed.velocity[0] =
               std::nextafter(state.committed.velocity[0],
                              std::numeric_limits<double>::infinity());
         }) &&
         changed([](auto &state) {
           state.trial.transported_cell_fields[1][0] =
               std::nextafter(state.trial.transported_cell_fields[1][0],
                              std::numeric_limits<double>::infinity());
         }) &&
         changed([](auto &state) {
           state.metadata.dt_s = std::nextafter(
               state.metadata.dt_s, std::numeric_limits<double>::infinity());
         }) &&
         changed([](auto &state) { state.history.density[1] = -0.0; }) &&
         changed([](auto &state) {
           state.committed.face_mass_flux[0] =
               std::nextafter(state.committed.face_mass_flux[0],
                              std::numeric_limits<double>::infinity());
         }) &&
         changed([](auto &state) {
           state.committed.transported_cell_fields[0].push_back(12.0);
         }) &&
         changed([](auto &state) {
           state.committed.transported_cell_fields.push_back({12.0});
         });
}

} // namespace hundun::test
