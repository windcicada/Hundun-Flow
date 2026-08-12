// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_checkpoint_v4.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::flow::detail {

struct PreparedCheckpointV4Section final {
  CheckpointSectionDescriptor descriptor;
  std::vector<std::uint8_t> local_payload;
};

inline constexpr CheckpointSectionId kReactingCompositionSection = 4101U;
inline constexpr CheckpointSectionId kReactingStateSection = 4102U;
inline constexpr CheckpointSectionId kReactingPressureSection = 4103U;
inline constexpr CheckpointSectionId kReactingTimeSection = 4104U;
inline constexpr CheckpointSectionId kReactingBackendSection = 4105U;

struct ReactingCheckpointV4Data final {
  std::uint64_t composition_fingerprint{};
  std::vector<std::string> species_names;
  std::vector<double> history_rho_y_kg_per_m3;
  std::vector<double> committed_rho_y_kg_per_m3;
  std::vector<double> history_rho_h_tc_j_per_m3;
  std::vector<double> committed_rho_h_tc_j_per_m3;
  double history_p0_pa{};
  double committed_p0_pa{};
  std::uint64_t step{};
  double time_s{};
  double previous_dt_s{};
  std::uint32_t bdf_order{};
  std::string backend_id;
  std::string mechanism_sha256;
  std::string mechanism_phase;
};

struct EncodedCheckpointV4Section final {
  CheckpointSectionDescriptor descriptor;
  std::vector<std::uint8_t> payload;
  std::uint64_t crc64{};
};

std::vector<EncodedCheckpointV4Section>
encode_reacting_checkpoint_v4(const ReactingCheckpointV4Data &);

bool restore_reacting_checkpoint_v4(
    const std::vector<EncodedCheckpointV4Section> &,
    std::uint64_t expected_composition_fingerprint,
    const std::string &expected_mechanism_sha256,
    ReactingCheckpointV4Data &publish_target, std::string &message) noexcept;

} // namespace hundun::flow::detail
