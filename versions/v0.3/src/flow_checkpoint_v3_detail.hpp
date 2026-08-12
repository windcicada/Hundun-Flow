// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/flow_ideal_gas_closure.hpp"
#include "hundun/flow_state.hpp"
#include "hundun/rt_types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace hundun::flow::detail {

struct CheckpointV3AuthorityGradient final {
  std::uint64_t link{};
  double value{};
};

struct CheckpointV3RankAuthority final {
  std::int32_t rank{};
  runtime::Box3 owned_box{};
  std::array<std::uint64_t, 9> fingerprints{};
  std::string state_filename;
  std::uint64_t state_logical_bytes{};
  std::uint64_t state_actual_bytes{};
  std::uint64_t state_crc64{};
  bool history_available{};
  bool committed_available{};
  std::vector<CheckpointV3AuthorityGradient> history;
  std::vector<CheckpointV3AuthorityGradient> committed;
};

enum class CheckpointV3WaleTransientField : std::uint8_t {
  nu_t_m2_per_s = 1,
  mu_sgs_pa_s = 2,
  mu_eff_pa_s = 3
};

struct CheckpointV3WaleIdentity final {
  double coefficient{};
  double turbulent_prandtl{};
  double turbulent_schmidt{};
  std::uint64_t numerical_config_crc64{};
  std::uint32_t transient_schema_version{1U};
  std::array<CheckpointV3WaleTransientField, 3> transient_fields{
      CheckpointV3WaleTransientField::nu_t_m2_per_s,
      CheckpointV3WaleTransientField::mu_sgs_pa_s,
      CheckpointV3WaleTransientField::mu_eff_pa_s};
};

struct CheckpointV3Manifest final {
  CheckpointV3Presence presence{CheckpointV3Presence::constant_static_ibm};
  std::int32_t rank_count{};
  runtime::Int3 process_grid{};
  std::uint64_t payload_report_fingerprint{};
  std::uint64_t payload_manifest_crc64{};
  AcceptedStepMetadata metadata;
  CheckpointV3ControlState control;
  std::vector<CheckpointV3RankAuthority> ranks;
  std::uint32_t ibm_section_count{};
  std::uint64_t ibm_section_bytes{};
  std::uint32_t wale_section_count{};
  std::uint64_t wale_section_bytes{};
  std::optional<CheckpointV3WaleIdentity> wale;
  std::uint32_t ideal_gas_section_count{};
  std::uint64_t ideal_gas_section_bytes{};
  std::optional<IdealGasClosureState> ideal_gas;
};

std::vector<std::uint8_t>
encode_checkpoint_v3_manifest(const CheckpointV3Manifest &);
CheckpointV3Manifest
decode_checkpoint_v3_manifest(const std::vector<std::uint8_t> &);
bool checkpoint_v3_manifest_equal(const CheckpointV3Manifest &,
                                  const CheckpointV3Manifest &) noexcept;
std::uint64_t
checkpoint_v3_manifest_crc64(const std::vector<std::uint8_t> &) noexcept;
std::size_t checkpoint_v3_manifest_presence_offset() noexcept;

class ImmersedFlowCheckpointPreparedRestore final {
public:
  ~ImmersedFlowCheckpointPreparedRestore() noexcept;
  ImmersedFlowCheckpointPreparedRestore(
      ImmersedFlowCheckpointPreparedRestore &&) noexcept;
  ImmersedFlowCheckpointPreparedRestore &
  operator=(ImmersedFlowCheckpointPreparedRestore &&) = delete;
  ImmersedFlowCheckpointPreparedRestore(
      const ImmersedFlowCheckpointPreparedRestore &) = delete;
  ImmersedFlowCheckpointPreparedRestore &
  operator=(const ImmersedFlowCheckpointPreparedRestore &) = delete;

private:
  struct Impl;
  explicit ImmersedFlowCheckpointPreparedRestore(
      std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend struct ImmersedFlowCheckpointAccess;
};

struct ImmersedFlowCheckpointAccess final {
  static CheckpointV3RankAuthority
  snapshot(const FixedStepImmersedFlow &, std::int32_t rank,
           runtime::Box3 owned_box);
  static ImmersedFlowCheckpointPreparedRestore
  prepare_restore(const FixedStepImmersedFlow &,
                  const CheckpointV3RankAuthority &);
  static void publish_restore(FixedStepImmersedFlow &,
                              ImmersedFlowCheckpointPreparedRestore &&,
                              const FlowState &) noexcept;
};

} // namespace hundun::flow::detail
