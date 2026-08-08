// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_checkpoint_v2.hpp"

#include "flow_checkpoint_v2_detail.hpp"
#include "rt_checkpoint_v2_protocol_detail.hpp"
#include "hundun/rt_mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>

namespace hundun::diagnostics::detail {
struct CheckpointV2Adapter final {
  static const flow::CheckpointV2Report *
  report_if_ready(const flow::CheckpointV2DiagnosticSource &source) noexcept {
    return source.impl_ ? &source.impl_->report : nullptr;
  }

  static const flow::CheckpointV2Report &
  report(const flow::CheckpointV2DiagnosticSource &source) {
    const auto *value = report_if_ready(source);
    if (value == nullptr)
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::invalid_input,
          "checkpoint-v2.diagnostics.stale-source", -1,
          "Checkpoint v2 diagnostic source has been moved from");
    return *value;
  }
};
} // namespace hundun::diagnostics::detail

namespace hundun::diagnostics {
namespace {

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
constexpr std::uint8_t kNoInjectedFault = 0U;
constexpr std::uint8_t kRawMpiInjectedFault = 1U;
constexpr std::uint8_t kRequestPreparationInjectedFault = 2U;
constexpr std::uint8_t kLocalRecordPreparationInjectedFault = 3U;
constexpr std::uint8_t kPostGatherPreparationInjectedFault = 4U;
constexpr std::uint8_t kFinalRecordPreparationInjectedFault = 5U;

std::uint8_t injected_fault{kNoInjectedFault};
int injected_fault_rank{-1};
std::uint64_t diagnostic_collective_calls{};

void before_collective(const runtime::MpiContext &mpi,
                       std::string_view operation) {
  ++diagnostic_collective_calls;
  if (injected_fault == kRawMpiInjectedFault &&
      (injected_fault_rank < 0 || injected_fault_rank == mpi.rank())) {
    injected_fault = kNoInjectedFault;
    runtime::check_mpi_result(MPI_ERR_OTHER, operation);
  }
}
#else
void before_collective(const runtime::MpiContext &, std::string_view) {}
#endif

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
bool consume_preparation_fault(
    const runtime::MpiContext &mpi, std::uint8_t phase) noexcept {
  if (injected_fault == phase &&
      (injected_fault_rank < 0 || injected_fault_rank == mpi.rank())) {
    injected_fault = kNoInjectedFault;
    return true;
  }
  return false;
}
#endif

void converge_preparation(const runtime::MpiContext &mpi,
                          bool local_ready) {
  const int candidate = local_ready ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  before_collective(mpi, "MPI_Allreduce(Checkpoint diagnostic preparation)");
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostic preparation)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::collective_operation,
        "checkpoint-v2.diagnostics.preparation", lowest,
        "Checkpoint v2 diagnostic collective preparation failed");
}

const flow::CheckpointV2Report &require_collective_source_readiness(
    const flow::CheckpointV2DiagnosticSource &source,
    const runtime::MpiContext &mpi) {
  const auto *report = detail::CheckpointV2Adapter::report_if_ready(source);
  const int candidate = report == nullptr ? mpi.rank() : mpi.size();
  int lowest = mpi.size();
  before_collective(mpi, "MPI_Allreduce(Checkpoint diagnostic source)");
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostic source)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input,
        "checkpoint-v2.diagnostics.stale-source", lowest,
        "Checkpoint v2 diagnostic source has been moved from");
  return *report;
}

constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

constexpr std::array<std::string_view, 26> kFingerprintIds{
    "checkpoint.collective-count",
    "checkpoint.crc-check-count",
    "checkpoint.disposition",
    "checkpoint.exact-size-eof-status",
    "checkpoint.file-count",
    "checkpoint.fingerprint-status",
    "checkpoint.global-actual-bytes",
    "checkpoint.global-logical-bytes",
    "checkpoint.local-actual-bytes",
    "checkpoint.local-crc64",
    "checkpoint.local-logical-bytes",
    "checkpoint.lowest-failing-rank",
    "checkpoint.manifest-crc-status",
    "checkpoint.manifest-crc64",
    "checkpoint.operation",
    "checkpoint.partition-status",
    "checkpoint.phase",
    "checkpoint.publication-status",
    "checkpoint.rank",
    "checkpoint.rank-crc-status",
    "checkpoint.reason",
    "checkpoint.rollback-status",
    "checkpoint.semantic-fingerprint",
    "checkpoint.step",
    "checkpoint.time",
    "checkpoint.transaction-entry-status"};

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::string phase_name(flow::CheckpointV2Phase phase) {
  switch (phase) {
  case flow::CheckpointV2Phase::none:
    return "none";
  case flow::CheckpointV2Phase::preflight:
    return "preflight";
  case flow::CheckpointV2Phase::transaction_entry:
    return "transaction-entry";
  case flow::CheckpointV2Phase::rank_payload:
    return "rank-payload";
  case flow::CheckpointV2Phase::rank_temporary_file:
    return "rank-temporary-file";
  case flow::CheckpointV2Phase::rank_publish:
    return "rank-publish";
  case flow::CheckpointV2Phase::manifest:
    return "manifest";
  case flow::CheckpointV2Phase::completed_marker:
    return "completed-marker";
  case flow::CheckpointV2Phase::marker_read:
    return "marker-read";
  case flow::CheckpointV2Phase::manifest_read:
    return "manifest-read";
  case flow::CheckpointV2Phase::rank_read:
    return "rank-read";
  case flow::CheckpointV2Phase::restore_prepare:
    return "restore-prepare";
  case flow::CheckpointV2Phase::restore_publish:
    return "restore-publish";
  }
  return "none";
}

std::string operation_name(flow::CheckpointV2Operation operation) {
  return operation == flow::CheckpointV2Operation::write ? "write" : "read";
}
std::string reason_name(flow::CheckpointV2FailureReason reason) {
  switch (reason) {
  case flow::CheckpointV2FailureReason::none:
    return "none";
  case flow::CheckpointV2FailureReason::invalid_input:
    return "invalid-input";
  case flow::CheckpointV2FailureReason::layout:
    return "layout";
  case flow::CheckpointV2FailureReason::state:
    return "state";
  case flow::CheckpointV2FailureReason::file_integrity:
    return "file-integrity";
  case flow::CheckpointV2FailureReason::filesystem:
    return "filesystem";
  }
  return "invalid-input";
}

DiagnosticFailureClass failure_class(flow::CheckpointV2FailureReason reason) {
  switch (reason) {
  case flow::CheckpointV2FailureReason::none:
    return DiagnosticFailureClass::none;
  case flow::CheckpointV2FailureReason::invalid_input:
  case flow::CheckpointV2FailureReason::state:
    return DiagnosticFailureClass::invalid_input;
  case flow::CheckpointV2FailureReason::layout:
    return DiagnosticFailureClass::layout;
  case flow::CheckpointV2FailureReason::file_integrity:
  case flow::CheckpointV2FailureReason::filesystem:
    return DiagnosticFailureClass::file_integrity;
  }
  return DiagnosticFailureClass::invalid_input;
}

std::string hex64(std::uint64_t value) {
  std::ostringstream output;
  output << std::hex << std::setfill('0') << std::setw(16) << value;
  return output.str();
}

void require_request(const flow::CheckpointV2Report &report,
                     const DiagnosticRequest &request,
                     DiagnosticScope expected) {
  if (request.level == DiagnosticLevel::bounded_state_sample ||
      !request.selected_fields.empty() || request.sample_budget != 0U)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "checkpoint-v2.diagnostics.capability", -1,
        "Checkpoint v2 diagnostic request uses an unsupported capability");
  try {
    validate(request);
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "checkpoint-v2.diagnostics.frame", -1,
        "Checkpoint v2 diagnostic request is invalid");
  }
  if (request.scope != expected)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "checkpoint-v2.diagnostics.capability", -1,
        "Checkpoint v2 diagnostic request uses an unsupported capability");
  if (request.frame.rank != report.rank() ||
      request.frame.step != report.step() ||
      bits(request.frame.time_s) != bits(report.time_s()) ||
      request.frame.phase != phase_name(report.phase()))
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "checkpoint-v2.diagnostics.frame", -1,
        "Checkpoint v2 diagnostic frame does not match the report");
}

enum class CollectiveRequestStatus : int {
  valid = 0,
  capability = 1,
  invalid_request = 2
};

bool valid_request_phase(std::string_view phase) noexcept {
  if (phase.empty() || phase.size() > 128U)
    return false;
  return std::all_of(phase.begin(), phase.end(), [](char character) {
    return (character >= 'a' && character <= 'z') ||
           (character >= '0' && character <= '9') || character == '_' ||
           character == '.' || character == '-';
  });
}

CollectiveRequestStatus
classify_collective_request(const flow::CheckpointV2Report &report,
                            const DiagnosticRequest &request) noexcept {
  const auto level = static_cast<std::uint8_t>(request.level);
  const auto scope = static_cast<std::uint8_t>(request.scope);
  if (level >
          static_cast<std::uint8_t>(DiagnosticLevel::bounded_state_sample) ||
      scope > static_cast<std::uint8_t>(DiagnosticScope::collective))
    return CollectiveRequestStatus::invalid_request;
  if (request.level == DiagnosticLevel::bounded_state_sample ||
      !request.selected_fields.empty() || request.sample_budget != 0U ||
      request.scope != DiagnosticScope::collective)
    return CollectiveRequestStatus::capability;
  if (request.frame.rank < 0 || !std::isfinite(request.frame.time_s) ||
      request.frame.time_s < 0.0 || !valid_request_phase(request.frame.phase) ||
      request.frame.rank != report.rank() ||
      request.frame.step != report.step() ||
      bits(request.frame.time_s) != bits(report.time_s()) ||
      request.frame.phase != phase_name(report.phase()))
    return CollectiveRequestStatus::invalid_request;
  return CollectiveRequestStatus::valid;
}

void require_collective_request_classification(
    const runtime::MpiContext &mpi, CollectiveRequestStatus local_status) {
  const int local = static_cast<int>(local_status);
  int lowest_status{};
  int highest_status{};
  before_collective(mpi,
                    "MPI_Allreduce(Checkpoint diagnostic request minimum)");
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &lowest_status, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostic request minimum)");
  before_collective(mpi,
                    "MPI_Allreduce(Checkpoint diagnostic request maximum)");
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &highest_status, 1, MPI_INT, MPI_MAX, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostic request maximum)");
  if (lowest_status == highest_status) {
    if (local_status == CollectiveRequestStatus::valid)
      return;
    if (local_status == CollectiveRequestStatus::capability)
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::capability,
          "checkpoint-v2.diagnostics.capability", -1,
          "Checkpoint v2 diagnostic request uses an unsupported capability");
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "checkpoint-v2.diagnostics.frame", -1,
        "Checkpoint v2 diagnostic request is invalid");
  }
  const int candidate =
      local_status == CollectiveRequestStatus::valid ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  before_collective(
      mpi, "MPI_Allreduce(Checkpoint diagnostic request disagreement)");
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostic request disagreement)");
  throw DiagnosticCollectionError(
      DiagnosticFailureClass::collective_operation,
      "checkpoint-v2.diagnostics.collective-agreement", lowest,
      "Checkpoint v2 collective diagnostic request disagrees");
}

std::uint64_t collective_signature(const flow::CheckpointV2Report &report,
    const DiagnosticRequest &request) {
  runtime::checkpoint_v2::Encoder encoder;
  constexpr std::string_view domain =
      "hundun.checkpoint-v2.diagnostic-agreement.v1";
  encoder.raw(domain.data(), domain.size());
  encoder.u8(0U);
  encoder.u32(1U);
  encoder.u8(static_cast<std::uint8_t>(request.level));
  encoder.u8(static_cast<std::uint8_t>(request.scope));
  encoder.u64(request.frame.step);
  encoder.f64(request.frame.time_s);
  encoder.string(std::string(request.frame.phase));
  encoder.u64(request.sample_budget);
  encoder.u8(static_cast<std::uint8_t>(report.operation()));
  encoder.u8(static_cast<std::uint8_t>(report.disposition()));
  encoder.u8(static_cast<std::uint8_t>(report.reason()));
  encoder.u8(static_cast<std::uint8_t>(report.phase()));
  encoder.i32(report.lowest_failing_rank());
  encoder.u64(report.step());
  encoder.f64(report.time_s());
  encoder.u64(report.global_logical_bytes());
  encoder.u64(report.global_actual_bytes());
  encoder.u64(report.manifest_crc64());
  encoder.u64(report.file_count());
  encoder.u64(report.crc_check_count());
  encoder.u64(report.collective_count());
  encoder.u8(static_cast<std::uint8_t>(report.manifest_crc_status()));
  encoder.u8(static_cast<std::uint8_t>(report.exact_size_and_eof_status()));
  encoder.u8(static_cast<std::uint8_t>(report.fingerprint_status()));
  encoder.u8(static_cast<std::uint8_t>(report.partition_status()));
  encoder.u8(static_cast<std::uint8_t>(report.transaction_entry_status()));
  encoder.u8(static_cast<std::uint8_t>(report.publication_status()));
  encoder.u8(static_cast<std::uint8_t>(report.rollback_status()));
  return runtime::checkpoint_v2::crc64_ecma(encoder.bytes().data(),
                                             encoder.bytes().size());
}

void require_collective_agreement(const flow::CheckpointV2Report &report,
                                  const runtime::MpiContext &mpi,
                                  const DiagnosticRequest &request) {
  require_collective_request_classification(
      mpi, classify_collective_request(report, request));
  bool local_valid = true;
  std::uint64_t signature{};
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(mpi, kRequestPreparationInjectedFault))
      throw std::bad_alloc();
#endif
    local_valid = mpi.rank() == report.rank();
    signature = collective_signature(report, request);
  } catch (...) {
    local_valid = false;
  }
  converge_preparation(mpi, local_valid);
  std::uint64_t root_signature = signature;
  before_collective(mpi, "MPI_Bcast(Checkpoint diagnostic agreement)");
  runtime::check_mpi_result(
      MPI_Bcast(&root_signature, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(Checkpoint diagnostic agreement)");
  const int candidate =
      local_valid && signature == root_signature ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  before_collective(mpi, "MPI_Allreduce(Checkpoint diagnostic agreement)");
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostic agreement)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::collective_operation,
        "checkpoint-v2.diagnostics.collective-agreement", lowest,
        "Checkpoint v2 collective diagnostic request or report disagrees");
}

void add_u64(DiagnosticFingerprintAccumulator &accumulator,
             std::string_view id, int rank, std::uint64_t value) {
  accumulator.add(id, static_cast<std::uint64_t>(rank), 0U,
                  describe_fp64(static_cast<double>(
                      static_cast<std::uint32_t>(value))));
  accumulator.add(id, static_cast<std::uint64_t>(rank), 1U,
                  describe_fp64(static_cast<double>(
                      static_cast<std::uint32_t>(value >> 32U))));
}
void add_count(DiagnosticFingerprintAccumulator &accumulator,
               std::string_view id, int rank, std::int64_t value) {
  accumulator.add(id, static_cast<std::uint64_t>(rank), 0U,
                  describe_fp64(static_cast<double>(value)));
}

DiagnosticFingerprintAccumulator
fingerprint(const flow::CheckpointV2Report &report) {
  DiagnosticFingerprintAccumulator result;
  add_u64(result, kFingerprintIds[0], report.rank(),
          report.collective_count());
  add_u64(result, kFingerprintIds[1], report.rank(),
          report.crc_check_count());
  add_count(result, kFingerprintIds[2], report.rank(),
            static_cast<std::uint8_t>(report.disposition()));
  add_count(result, kFingerprintIds[3], report.rank(),
            static_cast<std::uint8_t>(report.exact_size_and_eof_status()));
  add_u64(result, kFingerprintIds[4], report.rank(), report.file_count());
  add_count(result, kFingerprintIds[5], report.rank(),
            static_cast<std::uint8_t>(report.fingerprint_status()));
  add_u64(result, kFingerprintIds[6], report.rank(),
          report.global_actual_bytes());
  add_u64(result, kFingerprintIds[7], report.rank(),
          report.global_logical_bytes());
  add_u64(result, kFingerprintIds[8], report.rank(),
          report.local_actual_bytes());
  add_u64(result, kFingerprintIds[9], report.rank(), report.local_crc64());
  add_u64(result, kFingerprintIds[10], report.rank(),
          report.local_logical_bytes());
  add_count(result, kFingerprintIds[11], report.rank(),
            report.lowest_failing_rank());
  add_count(result, kFingerprintIds[12], report.rank(),
            static_cast<std::uint8_t>(report.manifest_crc_status()));
  add_u64(result, kFingerprintIds[13], report.rank(),
          report.manifest_crc64());
  add_count(result, kFingerprintIds[14], report.rank(),
            static_cast<std::uint8_t>(report.operation()));
  add_count(result, kFingerprintIds[15], report.rank(),
            static_cast<std::uint8_t>(report.partition_status()));
  add_count(result, kFingerprintIds[16], report.rank(),
            static_cast<std::uint8_t>(report.phase()));
  add_count(result, kFingerprintIds[17], report.rank(),
            static_cast<std::uint8_t>(report.publication_status()));
  add_count(result, kFingerprintIds[18], report.rank(),
            report.rank());
  add_count(result, kFingerprintIds[19], report.rank(),
            static_cast<std::uint8_t>(report.rank_crc_status()));
  add_count(result, kFingerprintIds[20], report.rank(),
            static_cast<std::uint8_t>(report.reason()));
  add_count(result, kFingerprintIds[21], report.rank(),
            static_cast<std::uint8_t>(report.rollback_status()));
  add_u64(result, kFingerprintIds[22], report.rank(),
          report.semantic_fingerprint());
  add_u64(result, kFingerprintIds[23], report.rank(), report.step());
  result.add(kFingerprintIds[24], static_cast<std::uint64_t>(report.rank()),
             0U, describe_fp64(report.time_s()));
  add_count(result, kFingerprintIds[25], report.rank(),
            static_cast<std::uint8_t>(report.transaction_entry_status()));
  return result;
}

DiagnosticInvariant invariant(std::string id,
                              flow::CheckpointV2CheckStatus status) {
  const double observed =
      status == flow::CheckpointV2CheckStatus::passed
          ? 1.0
          : status == flow::CheckpointV2CheckStatus::failed ? 0.0 : -1.0;
  return {std::move(id), "1", describe_fp64(observed), describe_fp64(1.0),
          InvariantRelation::equal,
          status == flow::CheckpointV2CheckStatus::passed};
}

DiagnosticRecord build_record(const flow::CheckpointV2Report &report,
                              const DiagnosticRequest &request) {
  DiagnosticRecord record;
  record.schema_version = 1U;
  record.module_kind = DiagnosticModuleKind::checkpoint;
  record.module_id = "checkpoint-v2";
  record.instance_id = "checkpoint-v2";
  record.level = request.level;
  record.scope = request.scope;
  record.rank = report.rank();
  record.step = report.step();
  record.time_s = describe_fp64(report.time_s());
  record.phase = phase_name(report.phase());
  record.status =
      report.disposition() == flow::CheckpointV2Disposition::completed
          ? DiagnosticStatus::ok
          : DiagnosticStatus::failed;
  record.failure.classification = failure_class(report.reason());
  record.failure.code =
      report.reason() == flow::CheckpointV2FailureReason::none
          ? "none"
          : "checkpoint-v2." + operation_name(report.operation()) + "." +
                phase_name(report.phase()) + "." + reason_name(report.reason());
  record.failure.lowest_failing_rank =
      request.scope == DiagnosticScope::collective
          ? report.lowest_failing_rank()
          : -1;
  if (request.level == DiagnosticLevel::summary) {
    const auto metric = [&](std::string id, std::string unit, double item) {
      record.metrics.push_back(
          {std::move(id), DiagnosticMetricKind::state_summary,
           std::move(unit), describe_fp64(item)});
    };
    metric("checkpoint.disposition", "count",
           static_cast<double>(static_cast<std::uint8_t>(
               report.disposition())));
    metric("checkpoint.global-actual-bytes-high", "byte",
           static_cast<double>(
               static_cast<std::uint32_t>(report.global_actual_bytes() >> 32)));
    metric("checkpoint.global-actual-bytes-low", "byte",
           static_cast<double>(
               static_cast<std::uint32_t>(report.global_actual_bytes())));
    metric("checkpoint.global-logical-bytes-high", "byte",
           static_cast<double>(static_cast<std::uint32_t>(
               report.global_logical_bytes() >> 32)));
    metric("checkpoint.global-logical-bytes-low", "byte",
           static_cast<double>(
               static_cast<std::uint32_t>(report.global_logical_bytes())));
    metric("checkpoint.local-actual-bytes-high", "byte",
           static_cast<double>(
               static_cast<std::uint32_t>(report.local_actual_bytes() >> 32)));
    metric("checkpoint.local-actual-bytes-low", "byte",
           static_cast<double>(
               static_cast<std::uint32_t>(report.local_actual_bytes())));
    metric("checkpoint.local-logical-bytes-high", "byte",
           static_cast<double>(static_cast<std::uint32_t>(
               report.local_logical_bytes() >> 32)));
    metric("checkpoint.local-logical-bytes-low", "byte",
           static_cast<double>(
               static_cast<std::uint32_t>(report.local_logical_bytes())));
    metric("checkpoint.operation", "count",
           static_cast<double>(
               static_cast<std::uint8_t>(report.operation())));
    metric("checkpoint.reason", "count",
           static_cast<double>(static_cast<std::uint8_t>(report.reason())));
  } else if (request.level == DiagnosticLevel::invariants) {
    record.invariants = {
        invariant("checkpoint.exact-size-eof",
                  report.exact_size_and_eof_status()),
        invariant("checkpoint.fingerprint", report.fingerprint_status()),
        invariant("checkpoint.manifest-crc", report.manifest_crc_status()),
        invariant("checkpoint.partition", report.partition_status()),
        invariant("checkpoint.publication", report.publication_status()),
        invariant("checkpoint.rank-crc", report.rank_crc_status()),
        invariant("checkpoint.rollback", report.rollback_status()),
        invariant("checkpoint.transaction-entry",
                  report.transaction_entry_status())};
  } else if (request.level == DiagnosticLevel::counters) {
    record.counters = {
        {"checkpoint.collective-count", "count", report.collective_count()},
        {"checkpoint.crc-check-count", "count", report.crc_check_count()},
        {"checkpoint.file-count", "count", report.file_count()},
        {"checkpoint.global-actual-bytes", "byte",
         report.global_actual_bytes()},
        {"checkpoint.global-logical-bytes", "byte",
         report.global_logical_bytes()},
        {"checkpoint.local-actual-bytes", "byte", report.local_actual_bytes()},
        {"checkpoint.local-crc64-high", "count", report.local_crc64() >> 32},
        {"checkpoint.local-crc64-low", "count",
         static_cast<std::uint32_t>(report.local_crc64())},
        {"checkpoint.local-logical-bytes", "byte",
         report.local_logical_bytes()},
        {"checkpoint.manifest-crc64-high", "count",
         report.manifest_crc64() >> 32},
        {"checkpoint.manifest-crc64-low", "count",
         static_cast<std::uint32_t>(report.manifest_crc64())}};
  } else {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "checkpoint-v2.diagnostics.level", -1,
        "Checkpoint v2 diagnostic sampling is unsupported");
  }
  const auto state_fingerprint = fingerprint(report).finish();
  record.identities.push_back(
      {"checkpoint-v2", hex64(report.semantic_fingerprint()), std::nullopt,
       std::nullopt, std::nullopt});
  record.state_fingerprint = state_fingerprint;
  return record;
}

void submit_validated(DiagnosticSink &sink, const DiagnosticRecord &record) {
  try {
    sink.submit(record);
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::sink_failure,
        "diagnostics.sink.submit", -1,
        "Checkpoint v2 diagnostic sink rejected a record");
  }
}

} // namespace

DiagnosticDescriptor describe_checkpoint_v2_diagnostics() noexcept {
  return {1U, DiagnosticModuleKind::checkpoint, "checkpoint-v2",
          "checkpoint-v2", kCapabilities};
}
DiagnosticDescriptor
describe_diagnostics(const flow::CheckpointV2DiagnosticSource &) noexcept {
  return describe_checkpoint_v2_diagnostics();
}
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::CheckpointV2DiagnosticSource &) {
  return {kFingerprintIds.begin(), kFingerprintIds.end()};
}

void collect_diagnostics(const flow::CheckpointV2DiagnosticSource &source,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const auto &report = detail::CheckpointV2Adapter::report(source);
  require_request(report, request, DiagnosticScope::local);
  auto record = build_record(report, request);
  validate(record);
  submit_validated(sink, record);
}

void collect_diagnostics(const flow::CheckpointV2DiagnosticSource &source,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const auto &report = require_collective_source_readiness(source, mpi);
  require_collective_agreement(report, mpi, request);
  const std::array<std::uint64_t, 4> local{
      report.local_logical_bytes(), report.local_actual_bytes(),
      report.local_crc64(),
      static_cast<std::uint64_t>(report.rank_crc_status())};
  DiagnosticRecord record;
  std::vector<std::uint64_t> gathered;
  bool local_ready = true;
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(mpi,
                                  kLocalRecordPreparationInjectedFault))
      throw std::bad_alloc();
#endif
    record = build_record(report, request);
    gathered.resize(static_cast<std::size_t>(mpi.size()) * local.size());
  } catch (...) {
    local_ready = false;
  }
  converge_preparation(mpi, local_ready);
  before_collective(mpi, "MPI_Allgather(Checkpoint diagnostic local fields)");
  runtime::check_mpi_result(
      MPI_Allgather(local.data(), static_cast<int>(local.size()),
                    MPI_UINT64_T, gathered.data(),
                    static_cast<int>(local.size()), MPI_UINT64_T,
                    mpi.comm()),
      "MPI_Allgather(Checkpoint diagnostic local fields)");
  std::uint64_t logical_sum{};
  std::uint64_t actual_sum{};
  flow::CheckpointV2CheckStatus rank_crc{
      flow::CheckpointV2CheckStatus::passed};
  std::array<std::uint64_t, 2> combined{};
  local_ready = true;
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(mpi,
                                  kPostGatherPreparationInjectedFault))
      throw std::bad_alloc();
#endif
    runtime::checkpoint_v2::Encoder crc_encoder;
    constexpr std::string_view crc_domain =
        "hundun.checkpoint-v2.diagnostic-rank-crcs.v1";
    crc_encoder.raw(crc_domain.data(), crc_domain.size());
    crc_encoder.u8(0U);
    crc_encoder.u32(1U);
    crc_encoder.u32(static_cast<std::uint32_t>(mpi.size()));
    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto offset = static_cast<std::size_t>(rank) * local.size();
      if (logical_sum >
              std::numeric_limits<std::uint64_t>::max() - gathered[offset] ||
          actual_sum >
              std::numeric_limits<std::uint64_t>::max() -
                  gathered[offset + 1U])
        throw std::overflow_error(
            "Checkpoint v2 diagnostic byte aggregation overflowed");
      logical_sum += gathered[offset];
      actual_sum += gathered[offset + 1U];
      crc_encoder.i32(rank);
      crc_encoder.u64(gathered[offset + 2U]);
      const auto status = static_cast<flow::CheckpointV2CheckStatus>(
          gathered[offset + 3U]);
      if (status == flow::CheckpointV2CheckStatus::failed)
        rank_crc = flow::CheckpointV2CheckStatus::failed;
      else if (status == flow::CheckpointV2CheckStatus::not_checked &&
               rank_crc == flow::CheckpointV2CheckStatus::passed)
        rank_crc = flow::CheckpointV2CheckStatus::not_checked;
    }
    const auto rank_crc_digest = runtime::checkpoint_v2::crc64_ecma(
        crc_encoder.bytes().data(), crc_encoder.bytes().size());
    if (request.level == DiagnosticLevel::summary) {
      for (auto &metric : record.metrics) {
        if (metric.id == "checkpoint.local-actual-bytes-high")
          metric.value = describe_fp64(static_cast<double>(
              static_cast<std::uint32_t>(actual_sum >> 32U)));
        else if (metric.id == "checkpoint.local-actual-bytes-low")
          metric.value = describe_fp64(
              static_cast<double>(static_cast<std::uint32_t>(actual_sum)));
        else if (metric.id == "checkpoint.local-logical-bytes-high")
          metric.value = describe_fp64(static_cast<double>(
              static_cast<std::uint32_t>(logical_sum >> 32U)));
        else if (metric.id == "checkpoint.local-logical-bytes-low")
          metric.value = describe_fp64(
              static_cast<double>(static_cast<std::uint32_t>(logical_sum)));
      }
    } else if (request.level == DiagnosticLevel::invariants) {
      for (auto &item : record.invariants)
        if (item.id == "checkpoint.rank-crc")
          item = invariant("checkpoint.rank-crc", rank_crc);
    } else if (request.level == DiagnosticLevel::counters) {
      for (auto &counter : record.counters) {
        if (counter.id == "checkpoint.local-actual-bytes")
          counter.value = actual_sum;
        else if (counter.id == "checkpoint.local-logical-bytes")
          counter.value = logical_sum;
        else if (counter.id == "checkpoint.local-crc64-high")
          counter.value = rank_crc_digest >> 32U;
        else if (counter.id == "checkpoint.local-crc64-low")
          counter.value = static_cast<std::uint32_t>(rank_crc_digest);
      }
    }
    const auto parts = fingerprint(report).parts();
    combined = {parts.xor64, parts.sum64};
  } catch (...) {
    local_ready = false;
  }
  converge_preparation(mpi, local_ready);
  before_collective(mpi, "MPI_Allreduce(Checkpoint diagnostics XOR)");
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, combined.data(), 1, MPI_UINT64_T, MPI_BXOR,
                    mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostics XOR)");
  before_collective(mpi, "MPI_Allreduce(Checkpoint diagnostics sum)");
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, combined.data() + 1, 1, MPI_UINT64_T,
                    MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Checkpoint diagnostics sum)");
  local_ready = true;
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(mpi,
                                  kFinalRecordPreparationInjectedFault))
      throw std::bad_alloc();
#endif
    DiagnosticFingerprintAccumulator accumulator;
    accumulator.combine({combined[0], combined[1]});
    record.state_fingerprint = accumulator.finish();
    record.identities.front().layout_fingerprint =
        record.state_fingerprint.hex;
    validate(record);
  } catch (...) {
    local_ready = false;
  }
  converge_preparation(mpi, local_ready);
  submit_validated(sink, record);
}

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
namespace detail {
void checkpoint_v2_diagnostics_set_fault_for_test(std::uint8_t fault,
                                                  int rank) noexcept {
  injected_fault = fault;
  injected_fault_rank = rank;
}
void checkpoint_v2_diagnostics_reset_for_test() noexcept {
  injected_fault = kNoInjectedFault;
  injected_fault_rank = -1;
  diagnostic_collective_calls = 0U;
}
std::uint64_t checkpoint_v2_diagnostics_collective_calls_for_test() noexcept {
  return diagnostic_collective_calls;
}
} // namespace detail
#endif

} // namespace hundun::diagnostics
