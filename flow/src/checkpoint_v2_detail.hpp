// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/checkpoint_v2.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"

namespace hundun::flow {

struct FlowState::Impl final {
  Impl(const runtime::FieldRegistry &, runtime::FieldLayoutSet, FlowFieldIds,
       AcceptedStepMetadata);

  const runtime::FieldRegistry *registry;
  runtime::FieldLayoutSet layout;
  FlowFieldIds fields;
  AcceptedStepMetadata metadata;
  runtime::FieldAccessPlan access;
  runtime::FieldStorage history;
  runtime::FieldStorage committed;
  runtime::FieldStorage trial;
  runtime::FieldStorage rollback_snapshot;
  bool rollback_snapshot_valid{};
  bool attempt_active{};
  std::uint64_t attempt_identity{};
  std::uint64_t diagnostic_mutation_identity{1U};
  bool commit_prepared{};
  AcceptedStepMetadata prepared_metadata{};
};

namespace detail {

bool validate_ideal_gas_restore_state(
    const runtime::MpiContext &, const mesh::MeshTopology &,
    const mesh::MeshGeometry &, const boundary::BoundaryRegistry &,
    double cp_J_per_kg_K, double gas_constant_J_per_kg_K,
    double configured_pressure_pa, const FlowLayerValues &history,
    const FlowLayerValues &committed, const IdealGasClosureState &,
    std::uint64_t &collective_count);

struct CheckpointV2ReportValues final {
  CheckpointV2Operation operation{CheckpointV2Operation::write};
  CheckpointV2Disposition disposition{CheckpointV2Disposition::failed};
  CheckpointV2FailureReason reason{CheckpointV2FailureReason::invalid_input};
  CheckpointV2Phase phase{CheckpointV2Phase::preflight};
  int rank{};
  int lowest_failing_rank{-1};
  std::uint64_t step{};
  double time_s{};
  std::uint64_t local_logical_bytes{};
  std::uint64_t local_actual_bytes{};
  std::uint64_t global_logical_bytes{};
  std::uint64_t global_actual_bytes{};
  std::uint64_t local_crc64{};
  std::uint64_t manifest_crc64{};
  std::uint64_t file_count{};
  std::uint64_t crc_check_count{};
  std::uint64_t collective_count{};
  CheckpointV2CheckStatus rank_crc{CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus manifest_crc{CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus exact_size_eof{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus fingerprint{CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus partition{CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus transaction_entry{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus publication{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus rollback{CheckpointV2CheckStatus::not_checked};
};

struct FlowStateCheckpointAccess final {
  static bool live(const FlowState &) noexcept;
  static const runtime::FieldRegistry &registry(const FlowState &);
  static runtime::FieldLayoutSet layout(const FlowState &);
  static bool attempt_active(const FlowState &) noexcept;
  static bool diagnostic_identity_can_advance(const FlowState &) noexcept;
  static bool read_transaction_ready(FlowState &) noexcept;
  static void enter_read_transaction(FlowState &);
  static std::uint64_t diagnostic_identity(const FlowState &) noexcept;
  static FlowState prepare_replacement(const FlowState &,
                                       const FlowLayerValues &history,
                                       const FlowLayerValues &committed,
                                       AcceptedStepMetadata);
  static void publish_replacement(FlowState &, FlowState &&) noexcept;
};

struct CheckpointV2Access final {
  static CheckpointV2Report make(CheckpointV2ReportValues) noexcept;
  static CheckpointV2ReadResult
  make_read(CheckpointV2Report, TimeControlState,
            std::optional<IdealGasClosureState>, bool restored) noexcept;
  static CheckpointV2Report failed(CheckpointV2Operation, int,
                                   CheckpointV2FailureReason,
                                   CheckpointV2Phase) noexcept;
  static CheckpointV2ReadResult failed_read(int, CheckpointV2FailureReason,
                                            CheckpointV2Phase);
  static void authenticate(CheckpointV2Report &) noexcept;
};

} // namespace detail
} // namespace hundun::flow

struct hundun::flow::CheckpointV2DiagnosticSource::Impl final {
  explicit Impl(hundun::flow::CheckpointV2Report supplied)
      : report(std::move(supplied)) {}
  hundun::flow::CheckpointV2Report report;
};
