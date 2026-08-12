// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/chem_composition.hpp"

#include <cstdint>
#include <vector>

namespace hundun::chemistry {

struct ChemistryIntervalRequest final {
  ThermochemicalPoint state;
  double start_time_s{};
  double duration_s{};
};

enum class ChemistryStatus : std::uint32_t {
  success = 0U,
  invalid_input = 1U,
  composition_mismatch = 2U,
  state_inversion_failure = 3U,
  integration_failure = 4U,
  non_finite_output = 5U,
  conservation_failure = 6U,
  workspace_failure = 7U
};

struct ChemistryIntervalReport final {
  ThermochemicalPoint final_state;
  std::vector<double> integrated_rho_y_delta_kg_per_m3;
  ChemistryStatus status{ChemistryStatus::invalid_input};
  double completed_duration_s{};
  std::uint32_t internal_step_count{};

  [[nodiscard]] bool succeeded() const noexcept {
    return status == ChemistryStatus::success;
  }
};

} // namespace hundun::chemistry
