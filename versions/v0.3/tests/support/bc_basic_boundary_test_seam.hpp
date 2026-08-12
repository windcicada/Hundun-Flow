// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/rt_error.hpp"

#include <exception>
#include <string_view>

namespace hundun::boundary::detail {

enum class OutletTopologyObservationForTest {
  unchanged,
  local_id_out_of_range,
  ownership_is_ghost,
  patch_id_is_not_outlet,
  outlet_patch_does_not_contain,
  global_id_differs,
  owned_cardinality_differs
};

std::string_view fixed_preflight_message(const runtime::Error &error) noexcept;
std::string_view fixed_preflight_message(const std::exception &error) noexcept;
std::string_view fixed_unknown_preflight_message() noexcept;

void set_next_outlet_topology_observation_raw(int observation);

inline void set_next_outlet_topology_observation_for_test(
    OutletTopologyObservationForTest observation) {
  set_next_outlet_topology_observation_raw(static_cast<int>(observation));
}

} // namespace hundun::boundary::detail
