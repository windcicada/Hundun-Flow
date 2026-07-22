// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/flow_state.hpp"
#include "hundun/flow/ideal_gas_closure.hpp"

#include <cstdint>
#include <cstring>
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

} // namespace hundun::test
