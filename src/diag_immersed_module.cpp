// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_immersed_module.hpp"

#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <optional>
#include <string>
#include <type_traits>

namespace hundun::diagnostics {
namespace {

constexpr std::string_view kModuleId = "hundun.flow.immersed";
constexpr std::string_view kInstanceId = "primary";
constexpr std::string_view kPhase = "immersed-flow.attempt-result";
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

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  constexpr std::uint64_t prime = UINT64_C(1099511628211);
  for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
    hash ^= (value >> (byte * 8U)) & UINT64_C(0xff);
    hash *= prime;
  }
  return hash;
}

std::uint64_t force_hash(
    const std::optional<flow::ForceAttemptReport> &force) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, force.has_value() ? 1U : 0U);
  if (!force.has_value())
    return hash;
  const auto add = [&](const immersed::ForceComponents &components) {
    for (const auto value : {components.pressure_N, components.total_N,
                             components.viscous_N}) {
      hash = mix(hash, bits(value.x));
      hash = mix(hash, bits(value.y));
      hash = mix(hash, bits(value.z));
    }
  };
  add(force->operator_force);
  add(force->budget_reaction);
  add(force->surface_traction);
  add(force->consistency);
  return hash;
}

const flow::StepAttemptReport &base(
    const flow::ImmersedFlowStepAttemptReport &report) {
  return std::visit(
      [](const auto &candidate) -> const flow::StepAttemptReport & {
        using Candidate = std::decay_t<decltype(candidate)>;
        if constexpr (std::is_same_v<Candidate, flow::StepAttemptReport>)
          return candidate;
        else if constexpr (std::is_same_v<
                               Candidate,
                               flow::MaterialDensityStepAttemptReport>)
          return candidate.flow();
        else
          return candidate.flow().flow();
      },
      report.base);
}

DiagnosticFailureClass failure_class(flow::StepFailureReason reason) noexcept {
  switch (reason) {
  case flow::StepFailureReason::none:
    return DiagnosticFailureClass::none;
  case flow::StepFailureReason::invalid_input:
    return DiagnosticFailureClass::invalid_input;
  case flow::StepFailureReason::non_finite_trial:
    return DiagnosticFailureClass::non_finite_state;
  case flow::StepFailureReason::momentum_linear_solve:
  case flow::StepFailureReason::pressure_linear_solve:
    return DiagnosticFailureClass::non_convergence;
  case flow::StepFailureReason::final_momentum_residual:
  case flow::StepFailureReason::final_transport_residual:
  case flow::StepFailureReason::final_continuity_residual:
  case flow::StepFailureReason::final_pressure_residual:
    return DiagnosticFailureClass::numerical_breakdown;
  case flow::StepFailureReason::final_conservation_defect:
    return DiagnosticFailureClass::conservation;
  case flow::StepFailureReason::boundary_backflow:
    return DiagnosticFailureClass::boundary;
  case flow::StepFailureReason::collective_operation:
    return DiagnosticFailureClass::collective_operation;
  case flow::StepFailureReason::transport_failure:
  case flow::StepFailureReason::density_closure_failure:
    return DiagnosticFailureClass::numerical_breakdown;
  }
  return DiagnosticFailureClass::invalid_input;
}

std::string failure_code(flow::StepFailureReason reason) {
  return reason == flow::StepFailureReason::none
             ? "none"
             : "immersed-flow.failure." +
                   std::to_string(static_cast<unsigned>(reason));
}

struct Snapshot final {
  double wall_maximum{};
  double wall_mean{};
  std::uint64_t classified_cells{};
  std::uint64_t active_cells{};
  std::uint64_t immersed_links{};
  std::uint64_t donor_references{};
  std::uint64_t wall_points{};
  runtime::Fp64ReductionCounters reduction;
  runtime::HaloPerformanceCounters halo;
  execution::AllocationCounters allocation;
};

Snapshot snapshot(const flow::ImmersedFlowDiagnosticSource &source) {
  return {source.maximum_wall_penetration_m_per_s(),
          source.mean_wall_penetration_m_per_s(),
          source.classified_cell_count(),
          source.active_cell_count(),
          source.immersed_link_count(),
          source.donor_reference_count(),
          source.wall_quadrature_point_count(),
          source.reduction_counters(),
          source.halo_counters(),
          source.allocation_counters()};
}

void require_request(const flow::ImmersedFlowDiagnosticSource &source,
                     const DiagnosticRequest &request,
                     DiagnosticScope expected) {
  try {
    if (request.scope != expected || request.frame.rank != source.rank() ||
        request.frame.step != source.committed_step() ||
        bits(request.frame.time_s) != bits(source.committed_time_s()) ||
        request.frame.phase != kPhase || !request.selected_fields.empty() ||
        request.sample_budget != 0U ||
        request.level == DiagnosticLevel::bounded_state_sample)
      throw runtime::Error("immersed diagnostic request mismatch");
    validate(request, describe_diagnostics(source));
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "immersed-flow.diagnostics.request", -1,
        "Immersed-flow diagnostic request is invalid");
  }
}

DiagnosticStateFingerprint fingerprint(
    const flow::ImmersedFlowDiagnosticSource &source,
    const Snapshot &observation) {
  DiagnosticFingerprintAccumulator accumulator;
  const auto &report = base(source.report());
  accumulator.add("continuity", source.committed_step(), 0U,
                  describe_fp64(report.final_continuity_normalized_l2));
  for (std::uint32_t component = 0U; component < 3U; ++component)
    accumulator.add("momentum", source.committed_step(), component,
                    describe_fp64(
                        report.final_momentum_normalized_l2[component]));
  accumulator.add("wall-maximum", source.committed_step(), 0U,
                  describe_fp64(observation.wall_maximum));
  accumulator.add("wall-mean", source.committed_step(), 0U,
                  describe_fp64(observation.wall_mean));
  return accumulator.finish();
}

void add_force_metrics(DiagnosticRecord &record, std::string_view authority,
                       const immersed::ForceComponents &force) {
  constexpr std::array<char, 3> axes{'x', 'y', 'z'};
  const std::array<std::pair<std::string_view, runtime::Real3>, 3> parts{{
      {"pressure", force.pressure_N},
      {"total", force.total_N},
      {"viscous", force.viscous_N},
  }};
  for (const auto &[part, value] : parts) {
    const std::array<double, 3> components{value.x, value.y, value.z};
    for (std::size_t component = 0U; component < components.size(); ++component)
      record.metrics.push_back(
          {"force." + std::string(authority) + "." + std::string(part) +
               "." + axes[component],
           DiagnosticMetricKind::conservation, "N",
           describe_fp64(components[component])});
  }
}

DiagnosticRecord build_record(const flow::ImmersedFlowDiagnosticSource &source,
                              const DiagnosticRequest &request,
                              const Snapshot &observation) {
  const auto &attempt = source.report();
  const auto &report = base(attempt);
  const bool committed =
      report.disposition == flow::StepAttemptDisposition::committed;
  DiagnosticRecord record;
  record.schema_version = kDiagnosticRecordSchemaV1;
  record.module_kind = DiagnosticModuleKind::piso;
  record.module_id = std::string(kModuleId);
  record.instance_id = std::string(kInstanceId);
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = request.frame.step;
  record.time_s = describe_fp64(request.frame.time_s);
  record.phase = std::string(kPhase);
  record.status = committed ? DiagnosticStatus::ok : DiagnosticStatus::failed;
  record.failure = {failure_class(report.reason), failure_code(report.reason),
                    request.scope == DiagnosticScope::collective && !committed
                        ? report.lowest_failing_rank
                        : -1};
  record.identities.push_back(
      {"immersed-authority", std::nullopt, std::nullopt,
       source.committed_step(), std::nullopt});
  record.state_fingerprint = fingerprint(source, observation);

  if (request.level == DiagnosticLevel::summary) {
    record.metrics.push_back(
        {"continuity.final-normalized-l2", DiagnosticMetricKind::residual, "1",
         committed ? describe_fp64(report.final_continuity_normalized_l2)
                   : DiagnosticFp64{}});
    constexpr std::array<char, 3> axes{'x', 'y', 'z'};
    for (std::size_t component = 0U; component < axes.size(); ++component)
      record.metrics.push_back(
          {std::string("momentum.final-normalized-l2.") + axes[component],
           DiagnosticMetricKind::residual, "1",
           committed
               ? describe_fp64(report.final_momentum_normalized_l2[component])
               : DiagnosticFp64{}});
    record.metrics.push_back(
        {"wall-penetration.maximum", DiagnosticMetricKind::residual, "m/s",
         describe_fp64(observation.wall_maximum)});
    record.metrics.push_back(
        {"wall-penetration.mean", DiagnosticMetricKind::residual, "m/s",
         describe_fp64(observation.wall_mean)});
    if (attempt.force.has_value()) {
      add_force_metrics(record, "budget-reaction",
                        attempt.force->budget_reaction);
      add_force_metrics(record, "consistency", attempt.force->consistency);
      add_force_metrics(record, "operator", attempt.force->operator_force);
      add_force_metrics(record, "surface-traction",
                        attempt.force->surface_traction);
    }
  } else if (request.level == DiagnosticLevel::invariants) {
    const auto add = [&](std::string id, double observed, double limit,
                         InvariantRelation relation) {
      DiagnosticInvariant value{std::move(id), "1", describe_fp64(observed),
                                describe_fp64(limit), relation, false};
      value.passed = evaluate_invariant(value);
      record.invariants.push_back(std::move(value));
    };
    add("continuity.finite", report.final_continuity_normalized_l2, 0.0,
        InvariantRelation::finite);
    add("pressure.correctors",
        static_cast<double>(report.pressure_corrector_count), 2.0,
        InvariantRelation::equal);
    add("wall-penetration.finite", observation.wall_maximum, 0.0,
        InvariantRelation::finite);
  } else if (request.level == DiagnosticLevel::counters) {
    const auto &reduction = observation.reduction;
    const auto &halo = observation.halo;
    const auto &allocation = observation.allocation;
    record.counters = {
        {"allocation.allocated-bytes", "byte", allocation.allocated_bytes},
        {"allocation.events", "count", allocation.allocation_events},
        {"allocation.live-bytes", "byte", allocation.live_bytes},
        {"allocation.peak-live-bytes", "byte", allocation.peak_live_bytes},
        {"construction.active-cells", "count", observation.active_cells},
        {"construction.classified-cells", "count",
         observation.classified_cells},
        {"construction.donor-references", "count",
         observation.donor_references},
        {"construction.immersed-links", "count",
         observation.immersed_links},
        {"construction.wall-quadrature-points", "count",
         observation.wall_points},
        {"halo.exchanges", "count", halo.completed_exchanges},
        {"halo.messages", "count",
         halo.send_messages + halo.receive_messages},
        {"halo.payload-bytes", "byte",
         halo.send_payload_bytes + halo.receive_payload_bytes},
        {"pressure.correctors", "count", report.pressure_corrector_count},
        {"reduction.calls", "count", reduction.collective_calls},
        {"reduction.logical-payload-bytes", "byte",
         reduction.logical_payload_bytes},
        {"reduction.scalars", "count", reduction.reduced_scalars},
    };
  } else {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "immersed-flow.diagnostics.level", -1,
        "Immersed-flow diagnostic sampling is unsupported");
  }
  std::sort(record.metrics.begin(), record.metrics.end(),
            [](const auto &left, const auto &right) {
              return left.id < right.id;
            });
  std::sort(record.invariants.begin(), record.invariants.end(),
            [](const auto &left, const auto &right) {
              return left.id < right.id;
            });
  return record;
}

void submit(DiagnosticSink &sink, const DiagnosticRecord &record,
            int failing_rank = -1) {
  try {
    sink.submit(record);
  } catch (...) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", failing_rank,
                                    "Diagnostic sink rejected a record");
  }
}

} // namespace

DiagnosticDescriptor describe_diagnostics(
    const flow::ImmersedFlowDiagnosticSource &) noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::piso, kModuleId,
          kInstanceId, kCapabilities};
}

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::ImmersedFlowDiagnosticSource &) {
  return {"continuity", "momentum", "wall-maximum", "wall-mean"};
}

void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &source,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  require_request(source, request, DiagnosticScope::local);
  const auto record = build_record(source, request, snapshot(source));
  validate(record, describe_diagnostics(source), request);
  submit(sink, record);
}

void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &source,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  bool ready = true;
  try {
    require_request(source, request, DiagnosticScope::collective);
  } catch (...) {
    ready = false;
  }
  const int candidate = ready ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(immersed diagnostic preflight)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "immersed-flow.diagnostics.preflight", lowest,
        "Collective immersed-flow diagnostic preflight failed");

  const auto &report = base(source.report());
  std::array<std::uint64_t, 11> local_identity{
      source.committed_step(), bits(source.committed_time_s()),
      static_cast<std::uint64_t>(report.disposition),
      static_cast<std::uint64_t>(report.reason),
      static_cast<std::uint64_t>(report.lowest_failing_rank + 1),
      report.pressure_corrector_count,
      bits(report.final_continuity_normalized_l2),
      bits(report.final_momentum_normalized_l2[0]),
      bits(report.final_momentum_normalized_l2[1]),
      bits(report.final_momentum_normalized_l2[2]),
      force_hash(source.report().force)};
  std::array<std::uint64_t, 11> minimum{};
  std::array<std::uint64_t, 11> maximum{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_identity.data(), minimum.data(),
                    static_cast<int>(local_identity.size()), MPI_UINT64_T,
                    MPI_MIN, mpi.comm()),
      "MPI_Allreduce(immersed diagnostic identity minimum)");
  runtime::check_mpi_result(
      MPI_Allreduce(local_identity.data(), maximum.data(),
                    static_cast<int>(local_identity.size()), MPI_UINT64_T,
                    MPI_MAX, mpi.comm()),
      "MPI_Allreduce(immersed diagnostic identity maximum)");
  if (minimum != maximum)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "immersed-flow.diagnostics.report-agreement", 0,
        "Collective immersed-flow reports disagree");

  const auto local_snapshot = snapshot(source);
  const std::array<std::uint64_t, 16> local_counts{
      local_snapshot.classified_cells,
      local_snapshot.active_cells,
      local_snapshot.immersed_links,
      local_snapshot.donor_references,
      local_snapshot.wall_points,
      local_snapshot.reduction.collective_calls,
      local_snapshot.reduction.reduced_scalars,
      local_snapshot.reduction.logical_payload_bytes,
      local_snapshot.halo.completed_exchanges,
      local_snapshot.halo.send_messages + local_snapshot.halo.receive_messages,
      local_snapshot.halo.send_payload_bytes +
          local_snapshot.halo.receive_payload_bytes,
      local_snapshot.allocation.allocation_events,
      local_snapshot.allocation.allocated_bytes,
      local_snapshot.allocation.live_bytes,
      local_snapshot.allocation.peak_live_bytes,
      local_snapshot.immersed_links};
  std::array<std::uint64_t, 16> global_counts{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_counts.data(), global_counts.data(),
                    static_cast<int>(local_counts.size()), MPI_UINT64_T,
                    MPI_SUM, mpi.comm()),
      "MPI_Allreduce(immersed diagnostic counters)");
  double local_wall[2]{
      local_snapshot.wall_maximum,
      local_snapshot.wall_mean *
          static_cast<double>(local_snapshot.immersed_links)};
  double global_wall_maximum{};
  double global_wall_sum{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_wall, &global_wall_maximum, 1, MPI_DOUBLE, MPI_MAX,
                    mpi.comm()),
      "MPI_Allreduce(immersed diagnostic wall maximum)");
  runtime::check_mpi_result(
      MPI_Allreduce(local_wall + 1, &global_wall_sum, 1, MPI_DOUBLE, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(immersed diagnostic wall sum)");
  Snapshot global_snapshot = local_snapshot;
  global_snapshot.wall_maximum = global_wall_maximum;
  global_snapshot.wall_mean =
      global_counts[15] == 0U
          ? 0.0
          : global_wall_sum / static_cast<double>(global_counts[15]);
  global_snapshot.classified_cells = global_counts[0];
  global_snapshot.active_cells = global_counts[1];
  global_snapshot.immersed_links = global_counts[2];
  global_snapshot.donor_references = global_counts[3];
  global_snapshot.wall_points = global_counts[4];
  global_snapshot.reduction = {global_counts[5], global_counts[6],
                               global_counts[7]};
  global_snapshot.halo.completed_exchanges = global_counts[8];
  global_snapshot.halo.send_messages = global_counts[9];
  global_snapshot.halo.receive_messages = 0U;
  global_snapshot.halo.send_payload_bytes = global_counts[10];
  global_snapshot.halo.receive_payload_bytes = 0U;
  global_snapshot.allocation.allocation_events = global_counts[11];
  global_snapshot.allocation.allocated_bytes = global_counts[12];
  global_snapshot.allocation.live_bytes = global_counts[13];
  global_snapshot.allocation.peak_live_bytes = global_counts[14];

  auto record = build_record(source, request, global_snapshot);
  bool valid = true;
  try {
    validate(record, describe_diagnostics(source), request);
  } catch (...) {
    valid = false;
  }
  const int record_candidate = valid ? mpi.size() : mpi.rank();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&record_candidate, &lowest, 1, MPI_INT, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(immersed diagnostic record)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::layout,
                                    "immersed-flow.diagnostics.record",
                                    lowest,
                                    "Collective immersed-flow record failed");
  bool sink_failed = false;
  try {
    sink.submit(record);
  } catch (...) {
    sink_failed = true;
  }
  const int sink_candidate = sink_failed ? mpi.rank() : mpi.size();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&sink_candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(immersed diagnostic sink)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", lowest,
                                    "Collective diagnostic sink failed");
}

namespace {

constexpr std::string_view kWallForceModuleId =
    "hundun.immersed.wall-force";
constexpr std::string_view kWallForcePhase = "immersed-flow.wall-force";
constexpr DiagnosticCapabilityFlags kWallForceCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

bool same_real3(runtime::Real3 left, runtime::Real3 right) noexcept {
  return bits(left.x) == bits(right.x) && bits(left.y) == bits(right.y) &&
         bits(left.z) == bits(right.z);
}

bool same_force(const immersed::ForceComponents &left,
                const immersed::ForceComponents &right) noexcept {
  return same_real3(left.pressure_N, right.pressure_N) &&
         same_real3(left.total_N, right.total_N) &&
         same_real3(left.viscous_N, right.viscous_N);
}

void require_wall_force_source(
    const flow::ImmersedFlowDiagnosticSource &source) {
  try {
    const auto &report = base(source.report());
    if (report.disposition != flow::StepAttemptDisposition::committed ||
        !source.report().force.has_value() ||
        !source.wall_force_available() || source.snapshot_seal() == 0U) {
      throw runtime::Error("wall-force source is unavailable");
    }
    const auto &sample = source.wall_force_sample();
    if (sample.lowest_failing_rank != -1 ||
        !same_force(source.report().force->surface_traction,
                    sample.surface_traction)) {
      throw runtime::Error("wall-force source is unauthenticated");
    }
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::unavailable,
        "stage3.wall-force.diagnostics.unavailable", -1,
        "Accepted wall-force diagnostic source is unavailable");
  }
}

void require_wall_force_request(
    const flow::ImmersedFlowDiagnosticSource &source,
    const DiagnosticRequest &request, DiagnosticScope scope) {
  try {
    require_wall_force_source(source);
    const DiagnosticDescriptor descriptor{
        kDiagnosticRecordSchemaV1, DiagnosticModuleKind::wall_force,
        kWallForceModuleId, kInstanceId, kWallForceCapabilities};
    validate(request, descriptor);
    if (request.scope != scope || request.frame.rank != source.rank() ||
        request.frame.step != source.committed_step() ||
        bits(request.frame.time_s) != bits(source.committed_time_s()) ||
        request.frame.phase != kWallForcePhase ||
        !request.selected_fields.empty() || request.sample_budget != 0U ||
        (request.level != DiagnosticLevel::summary &&
         request.level != DiagnosticLevel::counters)) {
      throw runtime::Error("wall-force request mismatch");
    }
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (...) {
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "stage3.wall-force.diagnostics.request", -1,
        "Stage 3 wall-force diagnostic request is invalid");
  }
}

void add_fingerprint_real3(DiagnosticFingerprintAccumulator &accumulator,
                           std::string_view field, std::uint32_t offset,
                           runtime::Real3 value) {
  accumulator.add(field, 0U, offset, describe_fp64(value.x));
  accumulator.add(field, 0U, offset + 1U, describe_fp64(value.y));
  accumulator.add(field, 0U, offset + 2U, describe_fp64(value.z));
}

void add_fingerprint_force(DiagnosticFingerprintAccumulator &accumulator,
                           std::string_view field,
                           const immersed::ForceComponents &force) {
  add_fingerprint_real3(accumulator, field, 0U, force.pressure_N);
  add_fingerprint_real3(accumulator, field, 3U, force.total_N);
  add_fingerprint_real3(accumulator, field, 6U, force.viscous_N);
}

DiagnosticStateFingerprint wall_force_fingerprint(
    const flow::ImmersedFlowDiagnosticSource &source) {
  DiagnosticFingerprintAccumulator accumulator;
  const auto &force = *source.report().force;
  const auto &sample = source.wall_force_sample();
  add_fingerprint_force(accumulator, "force.budget-reaction",
                        force.budget_reaction);
  add_fingerprint_force(accumulator, "force.consistency", force.consistency);
  add_fingerprint_force(accumulator, "force.operator", force.operator_force);
  add_fingerprint_force(accumulator, "force.surface-traction",
                        force.surface_traction);
  add_fingerprint_real3(accumulator, "moment", 0U,
                        sample.moment_about_global_origin.pressure_N_m);
  add_fingerprint_real3(accumulator, "moment", 3U,
                        sample.moment_about_global_origin.total_N_m);
  add_fingerprint_real3(accumulator, "moment", 6U,
                        sample.moment_about_global_origin.viscous_N_m);
  add_fingerprint_real3(accumulator, "area-closure", 0U,
                        sample.area_vector_closure_m2);
  accumulator.add("point-count", 0U, 0U,
                  describe_fp64(static_cast<double>(
                      sample.quadrature_point_count >> 32U)));
  accumulator.add("point-count", 0U, 1U,
                  describe_fp64(static_cast<double>(
                      static_cast<std::uint32_t>(
                          sample.quadrature_point_count))));
  accumulator.add("lowest-rank", 0U, 0U,
                  describe_fp64(
                      static_cast<double>(sample.lowest_failing_rank + 1)));
  accumulator.add("density-variant", 0U, 0U,
                  describe_fp64(static_cast<double>(source.density_model())));
  accumulator.add("wale", 0U, 0U,
                  describe_fp64(source.report().wale.has_value() ? 1.0 : 0.0));
  return accumulator.finish();
}

std::uint64_t wall_force_sample_hash(
    const immersed::WallForceSample &sample) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto add = [&](runtime::Real3 value) {
    hash = mix(hash, bits(value.x));
    hash = mix(hash, bits(value.y));
    hash = mix(hash, bits(value.z));
  };
  add(sample.surface_traction.pressure_N);
  add(sample.surface_traction.total_N);
  add(sample.surface_traction.viscous_N);
  add(sample.moment_about_global_origin.pressure_N_m);
  add(sample.moment_about_global_origin.total_N_m);
  add(sample.moment_about_global_origin.viscous_N_m);
  add(sample.area_vector_closure_m2);
  hash = mix(hash, sample.quadrature_point_count);
  hash = mix(hash,
             static_cast<std::uint64_t>(sample.lowest_failing_rank + 1));
  return hash == 0U ? 1U : hash;
}

void add_moment_metrics(DiagnosticRecord &record,
                        const immersed::MomentComponents &moment) {
  constexpr std::array<char, 3> axes{'x', 'y', 'z'};
  const std::array<std::pair<std::string_view, runtime::Real3>, 3> parts{{
      {"pressure", moment.pressure_N_m},
      {"total", moment.total_N_m},
      {"viscous", moment.viscous_N_m},
  }};
  for (const auto &[part, value] : parts) {
    const std::array<double, 3> components{value.x, value.y, value.z};
    for (std::size_t component = 0U; component < components.size(); ++component)
      record.metrics.push_back(
          {"moment." + std::string(part) + "." + axes[component],
           DiagnosticMetricKind::conservation, "N*m",
           describe_fp64(components[component])});
  }
}

DiagnosticRecord build_wall_force_record(
    const flow::ImmersedFlowDiagnosticSource &source,
    const DiagnosticRequest &request,
    std::optional<std::uint64_t> collective_seal = std::nullopt) {
  const auto &force = *source.report().force;
  const auto &sample = source.wall_force_sample();
  DiagnosticRecord record;
  record.schema_version = kDiagnosticRecordSchemaV1;
  record.module_kind = DiagnosticModuleKind::wall_force;
  record.module_id = std::string(kWallForceModuleId);
  record.instance_id = std::string(kInstanceId);
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = source.committed_step();
  record.time_s = describe_fp64(source.committed_time_s());
  record.phase = std::string(kWallForcePhase);
  record.status = DiagnosticStatus::ok;
  record.identities.push_back(
      {"wall-force-attempt", std::nullopt,
       collective_seal.value_or(source.snapshot_seal()),
       source.committed_step(), std::nullopt});
  record.state_fingerprint = wall_force_fingerprint(source);
  if (request.level == DiagnosticLevel::summary) {
    add_force_metrics(record, "budget-reaction", force.budget_reaction);
    add_force_metrics(record, "consistency", force.consistency);
    add_force_metrics(record, "operator", force.operator_force);
    add_force_metrics(record, "surface-traction", force.surface_traction);
    add_moment_metrics(record, sample.moment_about_global_origin);
    constexpr std::array<char, 3> axes{'x', 'y', 'z'};
    const std::array<double, 3> closure{sample.area_vector_closure_m2.x,
                                        sample.area_vector_closure_m2.y,
                                        sample.area_vector_closure_m2.z};
    for (std::size_t component = 0U; component < closure.size(); ++component)
      record.metrics.push_back(
          {std::string("surface-area-vector-closure.") + axes[component],
           DiagnosticMetricKind::conservation, "m2",
           describe_fp64(closure[component])});
  } else {
    record.counters = {
        {"density-variant", "count",
         static_cast<std::uint64_t>(source.density_model())},
        {"lowest-failing-rank", "count",
         static_cast<std::uint64_t>(sample.lowest_failing_rank + 1)},
        {"quadrature-points", "count", sample.quadrature_point_count},
        {"snapshot-seal-high", "count", source.snapshot_seal() >> 32U},
        {"snapshot-seal-low", "count",
         static_cast<std::uint32_t>(source.snapshot_seal())},
        {"wale-available", "count",
         source.report().wale.has_value() ? 1U : 0U},
    };
  }
  std::sort(record.metrics.begin(), record.metrics.end(),
            [](const auto &left, const auto &right) {
              return left.id < right.id;
            });
  return record;
}

} // namespace

DiagnosticDescriptor describe_diagnostics(
    const flow::ImmersedFlowDiagnosticSource &source,
    DiagnosticModuleKind kind) {
  if (kind != DiagnosticModuleKind::wall_force)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "stage3.immersed-attempt.diagnostics.kind", -1,
        "Diagnostic kind is not an immersed attempt provider");
  require_wall_force_source(source);
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::wall_force,
          kWallForceModuleId, kInstanceId, kWallForceCapabilities};
}

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::ImmersedFlowDiagnosticSource &source,
    DiagnosticModuleKind kind) {
  static_cast<void>(describe_diagnostics(source, kind));
  return {"area-closure", "density-variant", "force.budget-reaction",
          "force.consistency", "force.operator", "force.surface-traction",
          "lowest-rank", "moment", "point-count", "wale"};
}

void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &source,
                         DiagnosticModuleKind kind,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  if (kind != DiagnosticModuleKind::wall_force)
    static_cast<void>(describe_diagnostics(source, kind));
  require_wall_force_request(source, request, DiagnosticScope::local);
  const auto record = build_wall_force_record(source, request);
  validate(record, describe_diagnostics(source, kind), request);
  submit(sink, record);
}

void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &source,
                         DiagnosticModuleKind kind,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  bool ready = kind == DiagnosticModuleKind::wall_force;
  try {
    require_wall_force_request(source, request, DiagnosticScope::collective);
  } catch (...) {
    ready = false;
  }
  const int candidate = ready ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(Stage 3 wall-force diagnostic preflight)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "stage3.wall-force.diagnostics.preflight", lowest,
        "Collective Stage 3 wall-force diagnostic preflight failed");

  const std::array<std::uint64_t, 5> identity{
      source.committed_step(), bits(source.committed_time_s()),
      static_cast<std::uint64_t>(source.density_model()),
      force_hash(source.report().force),
      wall_force_sample_hash(source.wall_force_sample())};
  std::array<std::uint64_t, 5> minimum{};
  std::array<std::uint64_t, 5> maximum{};
  runtime::check_mpi_result(
      MPI_Allreduce(identity.data(), minimum.data(),
                    static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(Stage 3 wall-force identity minimum)");
  runtime::check_mpi_result(
      MPI_Allreduce(identity.data(), maximum.data(),
                    static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MAX,
                    mpi.comm()),
      "MPI_Allreduce(Stage 3 wall-force identity maximum)");
  if (minimum != maximum)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "stage3.wall-force.diagnostics.agreement", 0,
        "Collective Stage 3 wall-force sources disagree");

  const std::array<std::uint64_t, 2> local_seal{
      source.snapshot_seal(), source.snapshot_seal()};
  std::array<std::uint64_t, 2> combined_seal{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_seal.data(), combined_seal.data(), 1, MPI_UINT64_T,
                    MPI_BXOR, mpi.comm()),
      "MPI_Allreduce(Stage 3 wall-force seal xor)");
  runtime::check_mpi_result(
      MPI_Allreduce(local_seal.data() + 1, combined_seal.data() + 1, 1,
                    MPI_UINT64_T, MPI_SUM, mpi.comm()),
      "MPI_Allreduce(Stage 3 wall-force seal sum)");
  std::uint64_t collective_seal = UINT64_C(1469598103934665603);
  collective_seal = mix(collective_seal, combined_seal[0]);
  collective_seal = mix(collective_seal, combined_seal[1]);
  collective_seal = mix(collective_seal,
                        static_cast<std::uint64_t>(mpi.size()));
  if (collective_seal == 0U)
    collective_seal = 1U;

  std::optional<DiagnosticRecord> record;
  bool record_failed = false;
  try {
    record.emplace(
        build_wall_force_record(source, request, collective_seal));
    validate(*record, describe_diagnostics(source, kind), request);
  } catch (...) {
    record_failed = true;
  }
  const int record_candidate = record_failed ? mpi.rank() : mpi.size();
  lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&record_candidate, &lowest, 1, MPI_INT, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(Stage 3 wall-force record)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::layout,
                                    "stage3.wall-force.diagnostics.record",
                                    lowest,
                                    "Collective wall-force record failed");
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
      "MPI_Allreduce(Stage 3 wall-force sink)");
  if (lowest != mpi.size())
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", lowest,
                                    "Collective diagnostic sink failed");
}

} // namespace hundun::diagnostics
