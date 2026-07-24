// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/config/resolved_case.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/flow/ideal_gas_closure.hpp"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>

namespace hundun::boundary {
class BoundaryRegistry;
}
namespace hundun::mesh {
class MeshGeometry;
class MeshTopology;
}
namespace hundun::runtime {
class MpiContext;
class StructuredDecomposition;
}
namespace hundun::diagnostics::detail {
struct CheckpointV2Adapter;
}

namespace hundun::flow {

namespace detail {
struct CheckpointV2Access;
}

enum class CheckpointV2Operation : std::uint8_t { write = 0, read = 1 };
enum class CheckpointV2Disposition : std::uint8_t {
  completed = 0,
  failed = 1
};
enum class CheckpointV2FailureReason : std::uint8_t {
  none = 0,
  invalid_input = 1,
  layout = 2,
  state = 3,
  file_integrity = 4,
  filesystem = 5
};
enum class CheckpointV2CheckStatus : std::uint8_t {
  not_checked = 0,
  passed = 1,
  failed = 2
};
enum class CheckpointV2Phase : std::uint8_t {
  none = 0,
  preflight = 1,
  transaction_entry = 2,
  rank_payload = 3,
  rank_temporary_file = 4,
  rank_publish = 5,
  manifest = 6,
  completed_marker = 7,
  marker_read = 8,
  manifest_read = 9,
  rank_read = 10,
  restore_prepare = 11,
  restore_publish = 12
};

class CheckpointV2DiagnosticSource;

class CheckpointV2Report final {
public:
  CheckpointV2Report(const CheckpointV2Report &) = default;
  CheckpointV2Report(CheckpointV2Report &&) noexcept = default;
  CheckpointV2Report &operator=(const CheckpointV2Report &) = default;
  CheckpointV2Report &operator=(CheckpointV2Report &&) noexcept = default;
  CheckpointV2Operation operation() const noexcept;
  CheckpointV2Disposition disposition() const noexcept;
  CheckpointV2FailureReason reason() const noexcept;
  CheckpointV2Phase phase() const noexcept;
  int rank() const noexcept;
  int lowest_failing_rank() const noexcept;
  std::uint64_t step() const noexcept;
  double time_s() const noexcept;
  std::uint64_t local_logical_bytes() const noexcept;
  std::uint64_t local_actual_bytes() const noexcept;
  std::uint64_t global_logical_bytes() const noexcept;
  std::uint64_t global_actual_bytes() const noexcept;
  std::uint64_t local_crc64() const noexcept;
  std::uint64_t manifest_crc64() const noexcept;
  std::uint64_t file_count() const noexcept;
  std::uint64_t crc_check_count() const noexcept;
  std::uint64_t collective_count() const noexcept;
  CheckpointV2CheckStatus rank_crc_status() const noexcept;
  CheckpointV2CheckStatus manifest_crc_status() const noexcept;
  CheckpointV2CheckStatus exact_size_and_eof_status() const noexcept;
  CheckpointV2CheckStatus fingerprint_status() const noexcept;
  CheckpointV2CheckStatus partition_status() const noexcept;
  CheckpointV2CheckStatus transaction_entry_status() const noexcept;
  CheckpointV2CheckStatus publication_status() const noexcept;
  CheckpointV2CheckStatus rollback_status() const noexcept;
  std::uint64_t semantic_fingerprint() const noexcept;

private:
  CheckpointV2Report() = default;
  CheckpointV2Operation operation_{CheckpointV2Operation::write};
  CheckpointV2Disposition disposition_{CheckpointV2Disposition::failed};
  CheckpointV2FailureReason reason_{CheckpointV2FailureReason::invalid_input};
  CheckpointV2Phase phase_{CheckpointV2Phase::preflight};
  int rank_{};
  int lowest_failing_rank_{-1};
  std::uint64_t step_{};
  double time_s_{};
  std::uint64_t local_logical_bytes_{};
  std::uint64_t local_actual_bytes_{};
  std::uint64_t global_logical_bytes_{};
  std::uint64_t global_actual_bytes_{};
  std::uint64_t local_crc64_{};
  std::uint64_t manifest_crc64_{};
  std::uint64_t file_count_{};
  std::uint64_t crc_check_count_{};
  std::uint64_t collective_count_{};
  CheckpointV2CheckStatus rank_crc_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus manifest_crc_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus exact_size_and_eof_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus fingerprint_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus partition_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus transaction_entry_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus publication_status_{
      CheckpointV2CheckStatus::not_checked};
  CheckpointV2CheckStatus rollback_status_{
      CheckpointV2CheckStatus::not_checked};
  std::uint64_t semantic_fingerprint_{};
  friend struct detail::CheckpointV2Access;
  friend CheckpointV2DiagnosticSource
  checkpoint_v2_diagnostic_source(const CheckpointV2Report &);
  friend struct ::hundun::diagnostics::detail::CheckpointV2Adapter;
};

class CheckpointV2ReadResult final {
public:
  CheckpointV2ReadResult(const CheckpointV2ReadResult &) = default;
  CheckpointV2ReadResult(CheckpointV2ReadResult &&) noexcept = default;
  CheckpointV2ReadResult &operator=(const CheckpointV2ReadResult &) = default;
  CheckpointV2ReadResult &
  operator=(CheckpointV2ReadResult &&) noexcept = default;
  const CheckpointV2Report &report() const noexcept;
  bool restored() const noexcept;
  const TimeControlState &time_control_state() const;
  bool ideal_gas_closure_state_available() const noexcept;
  const IdealGasClosureState &ideal_gas_closure_state() const;

private:
  explicit CheckpointV2ReadResult(CheckpointV2Report) noexcept;
  CheckpointV2Report report_;
  bool restored_{};
  TimeControlState time_control_;
  std::optional<IdealGasClosureState> closure_;
  friend struct detail::CheckpointV2Access;
};

class CheckpointV2DiagnosticSource final {
public:
  ~CheckpointV2DiagnosticSource() noexcept;
  CheckpointV2DiagnosticSource(CheckpointV2DiagnosticSource &&) noexcept;
  CheckpointV2DiagnosticSource &
  operator=(CheckpointV2DiagnosticSource &&) = delete;
  CheckpointV2DiagnosticSource(const CheckpointV2DiagnosticSource &) = delete;
  CheckpointV2DiagnosticSource &
  operator=(const CheckpointV2DiagnosticSource &) = delete;

private:
  struct Impl;
  explicit CheckpointV2DiagnosticSource(std::unique_ptr<Impl>) noexcept;
  std::unique_ptr<Impl> impl_;
  friend CheckpointV2DiagnosticSource
  checkpoint_v2_diagnostic_source(const CheckpointV2Report &);
  friend struct ::hundun::diagnostics::detail::CheckpointV2Adapter;
};

CheckpointV2DiagnosticSource
checkpoint_v2_diagnostic_source(const CheckpointV2Report &);

CheckpointV2Report write_checkpoint_v2(
    const runtime::MpiContext &, const runtime::StructuredDecomposition &,
    const mesh::MeshTopology &, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &, const config::FlowCaseConfig &,
    const FlowState &, const TimeControlState &,
    std::optional<IdealGasClosureState>,
    const std::filesystem::path &checkpoint_directory);

CheckpointV2ReadResult read_checkpoint_v2(
    const runtime::MpiContext &, const runtime::StructuredDecomposition &,
    const mesh::MeshTopology &, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &, const config::FlowCaseConfig &,
    FlowState &, const std::filesystem::path &checkpoint_directory);

} // namespace hundun::flow
