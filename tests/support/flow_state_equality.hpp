// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/flow_state.hpp"

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
