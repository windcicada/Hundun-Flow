// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/flow_state.hpp"
#include "hundun/flow_ideal_gas_closure.hpp"

#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

namespace hundun::application::detail {

inline std::uint64_t performance_fp64_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

inline void append_performance_u64(std::string& output,
                                   std::uint64_t value) {
  for (int shift = 56; shift >= 0; shift -= 8)
    output.push_back(static_cast<char>(
        (value >> static_cast<unsigned>(shift)) & UINT64_C(0xff)));
}

inline void append_performance_flow_layer(
    std::string& output, const flow::FlowLayerValues& values) {
  const auto append_vector = [&](const std::vector<double>& vector) {
    append_performance_u64(output,
                           static_cast<std::uint64_t>(vector.size()));
    for (const double value : vector)
      append_performance_u64(output, performance_fp64_bits(value));
  };
  append_vector(values.density);
  append_vector(values.velocity);
  append_vector(values.mechanical_pressure);
  append_vector(values.face_velocity);
  append_vector(values.face_mass_flux);
  append_performance_u64(
      output,
      static_cast<std::uint64_t>(values.transported_cell_fields.size()));
  for (const auto& nested : values.transported_cell_fields)
    append_vector(nested);
}

inline std::string performance_flow_layer_encoding(
    const flow::FlowLayerValues& values) {
  std::string result;
  append_performance_flow_layer(result, values);
  return result;
}

inline void append_performance_closure_state(
    std::string& output,
    const std::optional<flow::IdealGasClosureState>& closure) {
  append_performance_u64(output, closure.has_value() ? 1U : 0U);
  if (!closure.has_value())
    return;
  append_performance_u64(
      output, static_cast<std::uint64_t>(closure->mode));
  append_performance_u64(
      output, performance_fp64_bits(closure->thermodynamic_pressure_pa));
  append_performance_u64(
      output, closure->target_mass_kg.has_value() ? 1U : 0U);
  if (closure->target_mass_kg.has_value())
    append_performance_u64(
        output, performance_fp64_bits(*closure->target_mass_kg));
  append_performance_u64(output, closure->revision);
}

inline std::string performance_closure_state_encoding(
    const std::optional<flow::IdealGasClosureState>& closure) {
  std::string result;
  append_performance_closure_state(result, closure);
  return result;
}

}  // namespace hundun::application::detail
