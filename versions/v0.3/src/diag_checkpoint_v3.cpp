// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_checkpoint_v3.hpp"

#include "hundun/rt_mpi_operation_error.hpp"

#include <mpi.h>

#include <array>
#include <cstring>
#include <optional>
#include <string>

namespace hundun::diagnostics {
namespace {

constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::string phase_name(flow::CheckpointV3Phase phase) {
  switch (phase) {
  case flow::CheckpointV3Phase::preflight:
    return "preflight";
  case flow::CheckpointV3Phase::rank_payload:
    return "rank-payload";
  case flow::CheckpointV3Phase::manifest:
    return "manifest";
  case flow::CheckpointV3Phase::completed_marker:
    return "completed-marker";
  case flow::CheckpointV3Phase::restore_prepare:
    return "restore-prepare";
  case flow::CheckpointV3Phase::restore_publish:
    return "restore-publish";
  }
  return "preflight";
}

DiagnosticFailureClass failure_class(
    flow::CheckpointV3FailureReason reason) noexcept {
  switch (reason) {
  case flow::CheckpointV3FailureReason::none:
    return DiagnosticFailureClass::none;
  case flow::CheckpointV3FailureReason::invalid_input:
  case flow::CheckpointV3FailureReason::state:
    return DiagnosticFailureClass::invalid_input;
  case flow::CheckpointV3FailureReason::layout:
  case flow::CheckpointV3FailureReason::fingerprint:
  case flow::CheckpointV3FailureReason::presence:
    return DiagnosticFailureClass::layout;
  case flow::CheckpointV3FailureReason::file_integrity:
  case flow::CheckpointV3FailureReason::filesystem:
    return DiagnosticFailureClass::file_integrity;
  case flow::CheckpointV3FailureReason::collective_operation:
    return DiagnosticFailureClass::collective_operation;
  }
  return DiagnosticFailureClass::invalid_input;
}

double check_value(flow::CheckpointV3CheckStatus status) noexcept {
  return status == flow::CheckpointV3CheckStatus::passed
             ? 1.0
             : status == flow::CheckpointV3CheckStatus::failed ? 0.0 : -1.0;
}

DiagnosticInvariant check_invariant(std::string id,
                                    flow::CheckpointV3CheckStatus status) {
  DiagnosticInvariant result{std::move(id), "1",
                             describe_fp64(check_value(status)),
                             describe_fp64(1.0), InvariantRelation::equal,
                             status == flow::CheckpointV3CheckStatus::passed};
  return result;
}

void require_request(const flow::CheckpointV3Report &report,
                     const DiagnosticRequest &request,
                     DiagnosticScope expected) {
  try {
    validate(request, describe_diagnostics(report));
    if (request.scope != expected ||
        request.frame.step != report.step() ||
        bits(request.frame.time_s) != bits(report.time_s()) ||
        request.frame.phase != phase_name(report.phase()) ||
        !request.selected_fields.empty() || request.sample_budget != 0U ||
        request.level == DiagnosticLevel::bounded_state_sample)
      throw std::runtime_error("checkpoint v3 request mismatch");
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "checkpoint-v3.diagnostics.request", -1,
        "Checkpoint v3 diagnostic request is invalid");
  }
}

DiagnosticRecord build_record(const flow::CheckpointV3Report &report,
                              const DiagnosticRequest &request) {
  const bool completed =
      report.disposition() == flow::CheckpointV3Disposition::completed;
  DiagnosticRecord record;
  record.schema_version = kDiagnosticRecordSchemaV1;
  record.module_kind = DiagnosticModuleKind::checkpoint;
  record.module_id = "checkpoint-v3";
  record.instance_id = "checkpoint-v3";
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = report.step();
  record.time_s = describe_fp64(report.time_s());
  record.phase = phase_name(report.phase());
  record.status = completed ? DiagnosticStatus::ok : DiagnosticStatus::failed;
  record.failure = {
      failure_class(report.reason()),
      report.reason() == flow::CheckpointV3FailureReason::none
          ? "none"
          : "checkpoint-v3.failure." +
                std::to_string(static_cast<unsigned>(report.reason())),
      request.scope == DiagnosticScope::collective && !completed
          ? report.lowest_failing_rank()
          : -1};
  record.identities.push_back(
      {"checkpoint-v3", std::nullopt, report.manifest_crc64(), std::nullopt,
       std::nullopt});
  DiagnosticFingerprintAccumulator fingerprint;
  fingerprint.add("checkpoint-crc", report.step(), 0U,
                  describe_fp64(
                      static_cast<double>(report.manifest_crc64() >> 32U)));
  fingerprint.add(
      "checkpoint-crc", report.step(), 1U,
      describe_fp64(static_cast<double>(
          static_cast<std::uint32_t>(report.manifest_crc64()))));
  fingerprint.add("checkpoint-status", report.step(), 0U,
                  describe_fp64(static_cast<double>(
                      static_cast<unsigned>(report.disposition()))));
  record.state_fingerprint = fingerprint.finish();

  if (request.level == DiagnosticLevel::summary) {
    record.metrics = {
        {"checkpoint.disposition", DiagnosticMetricKind::state_summary,
         "count",
         describe_fp64(static_cast<double>(
             static_cast<unsigned>(report.disposition())))},
        {"checkpoint.operation", DiagnosticMetricKind::state_summary, "count",
         describe_fp64(static_cast<double>(
             static_cast<unsigned>(report.operation())))},
        {"checkpoint.presence", DiagnosticMetricKind::state_summary, "count",
         describe_fp64(static_cast<double>(
             static_cast<unsigned>(report.presence())))},
        {"checkpoint.reason", DiagnosticMetricKind::state_summary, "count",
         describe_fp64(static_cast<double>(
             static_cast<unsigned>(report.reason())))},
    };
  } else if (request.level == DiagnosticLevel::invariants) {
    record.invariants = {
        check_invariant("checkpoint.crc", report.crc_status()),
        check_invariant("checkpoint.fingerprint", report.fingerprint_status()),
        check_invariant("checkpoint.partition", report.partition_status()),
        check_invariant("checkpoint.rollback", report.rollback_status()),
    };
  } else if (request.level == DiagnosticLevel::counters) {
    record.counters = {
        {"checkpoint.manifest-crc64-high", "count",
         report.manifest_crc64() >> 32U},
        {"checkpoint.manifest-crc64-low", "count",
         static_cast<std::uint32_t>(report.manifest_crc64())},
    };
  } else {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "checkpoint-v3.diagnostics.level", -1,
        "Checkpoint v3 diagnostic sampling is unsupported");
  }
  return record;
}

} // namespace

DiagnosticDescriptor
describe_diagnostics(const flow::CheckpointV3Report &) noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::checkpoint,
          "checkpoint-v3", "checkpoint-v3", kCapabilities};
}

std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const flow::CheckpointV3Report &) {
  return {"checkpoint-crc", "checkpoint-status"};
}

std::vector<DiagnosticModuleKind>
stage3_added_provider_inventory(flow::CheckpointV3Presence presence) {
  const auto value = static_cast<std::uint8_t>(presence);
  const bool ibm = value == 1U || value == 3U || value == 4U || value == 6U ||
                   value == 7U || value == 9U;
  const bool wale = value == 2U || value == 3U || value == 5U || value == 6U ||
                    value == 8U || value == 9U;
  std::vector<DiagnosticModuleKind> result;
  if (ibm) {
    result = {DiagnosticModuleKind::immersed_surface,
              DiagnosticModuleKind::ghost_stencil,
              DiagnosticModuleKind::local_flow_pattern,
              DiagnosticModuleKind::wall_force};
  }
  if (wale)
    result.push_back(DiagnosticModuleKind::les);
  return result;
}

std::vector<DiagnosticModuleKind>
stage3_added_provider_inventory(const flow::CheckpointV3Report &report) {
  return stage3_added_provider_inventory(report.presence());
}

void collect_diagnostics(const flow::CheckpointV3Report &report,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  require_request(report, request, DiagnosticScope::local);
  const auto record = build_record(report, request);
  validate(record, describe_diagnostics(report), request);
  try {
    sink.submit(record);
  } catch (...) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", -1,
                                    "Checkpoint v3 diagnostic sink failed");
  }
}

void collect_diagnostics(const flow::CheckpointV3Report &report,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  bool ready = true;
  try {
    require_request(report, request, DiagnosticScope::collective);
    if (report.disposition() == flow::CheckpointV3Disposition::failed &&
        report.lowest_failing_rank() < 0)
      ready = false;
  } catch (...) {
    ready = false;
  }
  const int candidate = ready ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(checkpoint v3 diagnostic preflight)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "checkpoint-v3.diagnostics.preflight", lowest,
        "Collective Checkpoint v3 diagnostic preflight failed");

  const std::array<std::uint64_t, 14> local_identity{
      static_cast<std::uint64_t>(request.level), report.step(),
      bits(report.time_s()), static_cast<std::uint64_t>(report.operation()),
      static_cast<std::uint64_t>(report.disposition()),
      static_cast<std::uint64_t>(report.reason()),
      static_cast<std::uint64_t>(report.phase()),
      static_cast<std::uint64_t>(report.presence()),
      static_cast<std::uint64_t>(report.lowest_failing_rank() + 1),
      report.manifest_crc64(),
      static_cast<std::uint64_t>(report.crc_status()),
      static_cast<std::uint64_t>(report.fingerprint_status()),
      static_cast<std::uint64_t>(report.partition_status()),
      static_cast<std::uint64_t>(report.rollback_status())};
  std::array<std::uint64_t, 14> minimum{};
  std::array<std::uint64_t, 14> maximum{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_identity.data(), minimum.data(),
                    static_cast<int>(local_identity.size()), MPI_UINT64_T,
                    MPI_MIN, mpi.comm()),
      "MPI_Allreduce(checkpoint v3 diagnostic identity minimum)");
  runtime::check_mpi_result(
      MPI_Allreduce(local_identity.data(), maximum.data(),
                    static_cast<int>(local_identity.size()), MPI_UINT64_T,
                    MPI_MAX, mpi.comm()),
      "MPI_Allreduce(checkpoint v3 diagnostic identity maximum)");
  if (minimum != maximum)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "checkpoint-v3.diagnostics.report-agreement", 0,
        "Collective Checkpoint v3 reports disagree");

  std::optional<DiagnosticRecord> record;
  bool record_failed = false;
  try {
    record.emplace(build_record(report, request));
    validate(*record, describe_diagnostics(report), request);
  } catch (...) {
    record_failed = true;
  }
  const int record_candidate = record_failed ? mpi.rank() : mpi.size();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&record_candidate, &lowest, 1, MPI_INT, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(checkpoint v3 diagnostic record)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::layout,
                                    "checkpoint-v3.diagnostics.record",
                                    lowest,
                                    "Collective Checkpoint v3 record failed");
  bool sink_failed = false;
  try {
    sink.submit(*record);
  } catch (...) {
    sink_failed = true;
  }
  const int sink_candidate = sink_failed ? mpi.rank() : mpi.size();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&sink_candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(checkpoint v3 diagnostic sink)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", lowest,
                                    "Collective diagnostic sink failed");
}

} // namespace hundun::diagnostics
