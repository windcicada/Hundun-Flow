// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/stage2_module_diagnostics.hpp"

#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <utility>

namespace hundun::diagnostics {
namespace {

constexpr DiagnosticCapabilityFlags kSummary =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary);
constexpr DiagnosticCapabilityFlags kSummaryCounters =
    kSummary |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters);
constexpr DiagnosticCapabilityFlags kSummaryInvariants =
    kSummary |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants);

DiagnosticDescriptor descriptor(DiagnosticModuleKind kind,
                                std::string_view module_id,
                                DiagnosticCapabilityFlags capabilities) noexcept {
  return {kDiagnosticRecordSchemaV1, kind, module_id, module_id, capabilities};
}

DiagnosticRecord base_record(const DiagnosticDescriptor& descriptor,
                             const DiagnosticRequest& request) {
  DiagnosticRecord record;
  record.module_kind = descriptor.module_kind;
  record.module_id = std::string(descriptor.module_id);
  record.instance_id = std::string(descriptor.instance_id);
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = request.frame.step;
  record.time_s = describe_fp64(request.frame.time_s);
  record.phase = std::string(request.frame.phase);
  return record;
}

DiagnosticMetric metric(std::string id, std::string unit, double value) {
  return {std::move(id), DiagnosticMetricKind::state_summary, std::move(unit),
          describe_fp64(value)};
}

DiagnosticCounter counter(std::string id, std::uint64_t value,
                          std::string unit = "count") {
  return {std::move(id), std::move(unit), value};
}

DiagnosticInvariant finite_invariant(std::string id, std::string unit,
                                     double value) {
  DiagnosticInvariant result{std::move(id),
                             std::move(unit),
                             describe_fp64(value),
                             {},
                             InvariantRelation::finite,
                             false};
  result.passed = evaluate_invariant(result);
  return result;
}

DiagnosticInvariant positive_invariant(std::string id, std::string unit,
                                       double value) {
  DiagnosticInvariant result{std::move(id),
                             std::move(unit),
                             describe_fp64(value),
                             {},
                             InvariantRelation::positive,
                             false};
  result.passed = evaluate_invariant(result);
  return result;
}

class FingerprintBuilder final {
 public:
  void add_fp64(std::string_view id, double value) {
    observe(id);
    accumulator_.add(id, ordinal_++, 0U, describe_fp64(value));
  }
  void add_u64(std::string_view id, std::uint64_t value) {
    observe(id);
    accumulator_.add(id, ordinal_, 0U,
                     describe_fp64(static_cast<double>(
                         static_cast<std::uint32_t>(value >> 32U))));
    accumulator_.add(id, ordinal_++, 1U,
                     describe_fp64(static_cast<double>(
                         static_cast<std::uint32_t>(value))));
  }
  void add_i64(std::string_view id, std::int64_t value) {
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    add_u64(id, bits);
  }
  void add_text(std::string_view id, std::string_view length_id,
                std::string_view value) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (char raw : value) {
      hash ^= static_cast<unsigned char>(raw);
      hash *= 1099511628211ULL;
    }
    add_u64(id, hash);
    add_u64(length_id, static_cast<std::uint64_t>(value.size()));
  }
  DiagnosticStateFingerprint finish() const { return accumulator_.finish(); }
  DiagnosticStateFingerprint finish(
      const std::vector<std::string_view>& advertised) const {
    auto published = advertised;
    auto observed = observed_;
    const auto normalize = [](auto& ids) {
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    };
    normalize(published);
    normalize(observed);
    if (published != advertised || observed != published)
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::invalid_input,
          "diagnostics.fingerprint.fields", -1,
          "diagnostic fingerprint field IDs do not match hashed fields");
    return finish();
  }
  DiagnosticFingerprintParts parts() const noexcept {
    return accumulator_.parts();
  }

 private:
  void observe(std::string_view id) { observed_.push_back(id); }

  DiagnosticFingerprintAccumulator accumulator_;
  std::vector<std::string_view> observed_;
  std::uint64_t ordinal_{};
};

std::string request_agreement_bytes(const DiagnosticRequest& request) {
  std::string result;
  const auto append_u64 = [&](std::uint64_t value) {
    for (unsigned byte = 0; byte < 8U; ++byte)
      result.push_back(static_cast<char>((value >> (8U * byte)) & 0xffU));
  };
  append_u64(static_cast<std::uint64_t>(request.level));
  append_u64(static_cast<std::uint64_t>(request.scope));
  append_u64(request.frame.step);
  std::uint64_t time_bits{};
  std::memcpy(&time_bits, &request.frame.time_s, sizeof(time_bits));
  append_u64(time_bits);
  append_u64(request.sample_budget);
  append_u64(request.frame.phase.size());
  result.append(request.frame.phase);
  append_u64(request.selected_fields.size());
  for (const auto field : request.selected_fields) {
    append_u64(field.size());
    result.append(field);
  }
  return result;
}

void require_exact_bytes_agreement(const runtime::MpiContext& mpi,
                                   std::string_view local,
                                   std::string_view subject) {
  std::uint64_t root_size =
      mpi.rank() == 0 ? static_cast<std::uint64_t>(local.size()) : 0U;
  runtime::check_mpi_result(
      MPI_Bcast(&root_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast diagnostic agreement size");
  bool allocation_ok =
      root_size <= static_cast<std::uint64_t>(
                       std::numeric_limits<std::size_t>::max()) &&
      root_size <= static_cast<std::uint64_t>(
                       std::numeric_limits<int>::max());
  std::string root;
  if (allocation_ok) {
    try {
      root.resize(static_cast<std::size_t>(root_size));
    } catch (...) {
      allocation_ok = false;
    }
  }
  auto status = runtime::collective_status(
      mpi, allocation_ok,
      std::string(subject) + " agreement allocation failed");
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::collective_operation,
        "diagnostics.collective.agreement-allocation",
        status.failing_rank, status.message);
  if (mpi.rank() == 0)
    std::copy(local.begin(), local.end(), root.begin());
  runtime::check_mpi_result(
      MPI_Bcast(root.data(), static_cast<int>(root_size), MPI_BYTE, 0,
                mpi.comm()),
      "MPI_Bcast diagnostic agreement bytes");
  status = runtime::collective_status(
      mpi, local == root, std::string(subject) + " differs across ranks");
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "diagnostics.collective.agreement", status.failing_rank,
        status.message);
}

template <class Fingerprint>
DiagnosticStateFingerprint authoritative_fingerprint(
    const std::vector<std::string_view>& field_ids, Fingerprint&& fingerprint) {
  FingerprintBuilder builder;
  std::forward<Fingerprint>(fingerprint)(builder);
  return builder.finish(field_ids);
}

void add_solve_fingerprint(FingerprintBuilder& fp, std::string_view prefix,
                           const linear::SolveReport& solve) {
  fp.add_text("solve.prefix", "solve.prefix.length", prefix);
  fp.add_u64("solve.reason",
             static_cast<std::uint64_t>(solve.reason));
  fp.add_u64("solve.iterations", solve.iterations);
  fp.add_fp64("solve.initial_residual", solve.initial_residual);
  fp.add_fp64("solve.recursive_residual", solve.recursive_residual);
  fp.add_fp64("solve.final_residual", solve.final_residual);
  fp.add_u64("solve.matvec_count", solve.matvec_count);
  fp.add_u64("solve.preconditioner_apply_count",
             solve.preconditioner_apply_count);
  fp.add_u64("solve.global_reduction_count",
             solve.global_reduction_count);
  fp.add_i64("solve.lowest_failing_rank", solve.lowest_failing_rank);
}

template <class Fingerprint, class Fill>
void collect_common(const DiagnosticDescriptor& descriptor,
                    const DiagnosticRequest& request, DiagnosticSink& sink,
                    const std::vector<std::string_view>& field_ids,
                    Fingerprint&& fingerprint,
                    Fill&& fill) {
  validate(request, descriptor);
  DiagnosticRecord record = base_record(descriptor, request);
  std::forward<Fill>(fill)(record);
  std::sort(record.invariants.begin(), record.invariants.end(),
            [](const auto& left, const auto& right) {
              return left.id < right.id;
            });
  std::sort(record.metrics.begin(), record.metrics.end(),
            [](const auto& left, const auto& right) {
              return left.id < right.id;
            });
  std::sort(record.counters.begin(), record.counters.end(),
            [](const auto& left, const auto& right) {
              return left.id < right.id;
            });
  std::sort(record.identities.begin(), record.identities.end(),
            [](const auto& left, const auto& right) {
              return left.subject_id < right.subject_id;
            });
  record.state_fingerprint =
      authoritative_fingerprint(field_ids,
                                std::forward<Fingerprint>(fingerprint));
  validate(record, descriptor, request);
  try {
    sink.submit(record);
  } catch (const DiagnosticCollectionError&) {
    throw;
  } catch (const std::exception& error) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", -1,
                                    error.what());
  }
}

void require_summary(const DiagnosticRequest& request) {
  if (request.level != DiagnosticLevel::summary)
    throw DiagnosticCollectionError(DiagnosticFailureClass::capability,
                                    "diagnostics.request.capability", -1,
                                    "unsupported diagnostic request");
}

template <class Integer>
double as_double(Integer value) {
  return static_cast<double>(value);
}

std::uint64_t checked_u64(std::size_t value) {
  return static_cast<std::uint64_t>(value);
}

std::string failure_code(flow::StepFailureReason reason) {
  switch (reason) {
  case flow::StepFailureReason::none:
    return "none";
  case flow::StepFailureReason::invalid_input:
    return "flow.invalid-input";
  case flow::StepFailureReason::momentum_linear_solve:
    return "flow.momentum-linear-solve";
  case flow::StepFailureReason::pressure_linear_solve:
    return "flow.pressure-linear-solve";
  case flow::StepFailureReason::non_finite_trial:
    return "flow.non-finite-trial";
  case flow::StepFailureReason::boundary_backflow:
    return "flow.boundary-backflow";
  case flow::StepFailureReason::transport_failure:
    return "flow.transport-failure";
  case flow::StepFailureReason::final_momentum_residual:
    return "flow.final-momentum-residual";
  case flow::StepFailureReason::final_transport_residual:
    return "flow.final-transport-residual";
  case flow::StepFailureReason::final_conservation_defect:
    return "flow.final-conservation-defect";
  case flow::StepFailureReason::final_continuity_residual:
    return "flow.final-continuity-residual";
  case flow::StepFailureReason::final_pressure_residual:
    return "flow.final-pressure-residual";
  case flow::StepFailureReason::collective_operation:
    return "flow.collective-operation";
  case flow::StepFailureReason::density_closure_failure:
    return "flow.density-closure-failure";
  }
  return "flow.invalid-input";
}

DiagnosticFailureClass failure_class(flow::StepFailureReason reason) {
  switch (reason) {
  case flow::StepFailureReason::none:
    return DiagnosticFailureClass::none;
  case flow::StepFailureReason::invalid_input:
    return DiagnosticFailureClass::invalid_input;
  case flow::StepFailureReason::non_finite_trial:
    return DiagnosticFailureClass::non_finite_state;
  case flow::StepFailureReason::boundary_backflow:
    return DiagnosticFailureClass::boundary;
  case flow::StepFailureReason::collective_operation:
    return DiagnosticFailureClass::collective_operation;
  case flow::StepFailureReason::momentum_linear_solve:
  case flow::StepFailureReason::pressure_linear_solve:
  case flow::StepFailureReason::transport_failure:
  case flow::StepFailureReason::density_closure_failure:
    return DiagnosticFailureClass::numerical_breakdown;
  case flow::StepFailureReason::final_momentum_residual:
  case flow::StepFailureReason::final_transport_residual:
  case flow::StepFailureReason::final_continuity_residual:
  case flow::StepFailureReason::final_pressure_residual:
    return DiagnosticFailureClass::non_convergence;
  case flow::StepFailureReason::final_conservation_defect:
    return DiagnosticFailureClass::conservation;
  }
  return DiagnosticFailureClass::invalid_input;
}

void set_failure(DiagnosticRecord& record, flow::StepFailureReason reason,
                 int lowest_rank) {
  if (reason == flow::StepFailureReason::none)
    return;
  record.status = DiagnosticStatus::failed;
  record.failure = {failure_class(reason), failure_code(reason),
                    record.scope == DiagnosticScope::collective
                        ? std::max(0, lowest_rank)
                        : -1};
}

}  // namespace

DiagnosticDescriptor describe_diagnostics(const runtime::MpiContext&) noexcept {
  return descriptor(DiagnosticModuleKind::mpi, "hundun.runtime.mpi_context",
                    kSummaryCounters |
                        static_cast<DiagnosticCapabilityFlags>(
                            DiagnosticCapability::collective));
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const runtime::MpiContext&) {
  return {"collective_calls", "logical_payload_bytes", "rank",
          "reduced_scalars", "size", "thread_level"};
}
void collect_diagnostics(const runtime::MpiContext& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fingerprint) {
                   const auto values = source.fp64_reduction_counters();
                   fingerprint.add_i64("rank", source.rank());
                   fingerprint.add_i64("size", source.size());
                   fingerprint.add_i64("thread_level", source.thread_level());
                   fingerprint.add_u64("collective_calls",
                                       values.collective_calls);
                   fingerprint.add_u64("logical_payload_bytes",
                                       values.logical_payload_bytes);
                   fingerprint.add_u64("reduced_scalars",
                                       values.reduced_scalars);
                 },
                 [&](DiagnosticRecord& record) {
                   if (request.scope != DiagnosticScope::local)
                     throw DiagnosticCollectionError(
                         DiagnosticFailureClass::capability,
                         "diagnostics.request.collective-context", -1,
                         "collective context is required");
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("rank", "count", as_double(source.rank())),
                         metric("size", "count", as_double(source.size())),
                         metric("thread_level", "count",
                                as_double(source.thread_level()))};
                   } else {
                     const auto values = source.fp64_reduction_counters();
                     record.counters = {
                         counter("collective_calls", values.collective_calls),
                         counter("logical_payload_bytes",
                                 values.logical_payload_bytes, "byte"),
                         counter("reduced_scalars", values.reduced_scalars)};
                   }
                 });
}
void collect_diagnostics(const runtime::MpiContext& source,
                         const runtime::MpiContext& collective_context,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  bool local_ok = true;
  std::string local_message;
  try {
    validate(request, describe_diagnostics(source));
    if (request.scope != DiagnosticScope::collective)
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::invalid_request,
          "diagnostics.request.scope", -1,
          "collective overload requires collective scope");
    if (request.frame.rank != collective_context.rank() ||
        source.rank() != collective_context.rank() ||
        source.size() != collective_context.size())
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::invalid_request,
          "diagnostics.request.rank", -1,
          "collective diagnostic rank/provider mismatch");
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  }
  auto status =
      runtime::collective_status(collective_context, local_ok, local_message);
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "diagnostics.collective.preflight", status.failing_rank,
        status.message);
  std::string request_bytes;
  try {
    request_bytes = request_agreement_bytes(request);
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  }
  status =
      runtime::collective_status(collective_context, local_ok, local_message);
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "diagnostics.collective.request", status.failing_rank,
        status.message);
  require_exact_bytes_agreement(collective_context, request_bytes,
                                "diagnostic request");

  const std::array<int, 3> local{source.rank(), source.size(),
                                 source.thread_level()};
  std::array<int, 3> minimum = local;
  std::array<int, 3> maximum = local;
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), minimum.data(), 3, MPI_INT, MPI_MIN,
                    collective_context.comm()),
      "MPI_Allreduce MPI diagnostic minimum");
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), maximum.data(), 3, MPI_INT, MPI_MAX,
                    collective_context.comm()),
      "MPI_Allreduce MPI diagnostic maximum");
  if (minimum[1] != maximum[1] || minimum[2] != maximum[2])
    local_ok = false;
  status =
      runtime::collective_status(collective_context, local_ok,
                                 "MPI diagnostic provider differs");
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input,
        "diagnostics.collective.provider", status.failing_rank,
        status.message);

  const auto counters = source.fp64_reduction_counters();
  std::array<std::uint64_t, 3> local_counters{
      counters.collective_calls, counters.logical_payload_bytes,
      counters.reduced_scalars};
  std::array<std::uint64_t, 3> global_counters{};
  runtime::check_mpi_result(
      MPI_Allreduce(local_counters.data(), global_counters.data(), 3,
                    MPI_UINT64_T, MPI_SUM, collective_context.comm()),
      "MPI_Allreduce MPI diagnostic counters");

  FingerprintBuilder local_fingerprint;
  local_fingerprint.add_i64("rank", source.rank());
  local_fingerprint.add_i64("size", source.size());
  local_fingerprint.add_i64("thread_level", source.thread_level());
  local_fingerprint.add_u64("collective_calls",
                            counters.collective_calls);
  local_fingerprint.add_u64("logical_payload_bytes",
                            counters.logical_payload_bytes);
  local_fingerprint.add_u64("reduced_scalars",
                            counters.reduced_scalars);
  static_cast<void>(
      local_fingerprint.finish(diagnostic_fingerprint_field_ids(source)));
  const auto local_parts = local_fingerprint.parts();
  std::array<std::uint64_t, 2> combined_parts{};
  runtime::check_mpi_result(
      MPI_Allreduce(&local_parts.xor64, &combined_parts[0], 1,
                    MPI_UINT64_T, MPI_BXOR, collective_context.comm()),
      "MPI_Allreduce MPI diagnostic fingerprint xor");
  runtime::check_mpi_result(
      MPI_Allreduce(&local_parts.sum64, &combined_parts[1], 1,
                    MPI_UINT64_T, MPI_SUM, collective_context.comm()),
      "MPI_Allreduce MPI diagnostic fingerprint sum");
  DiagnosticFingerprintAccumulator combined_fingerprint;
  combined_fingerprint.combine(
      {combined_parts[0], combined_parts[1]});

  DiagnosticRecord record;
  try {
    record = base_record(describe_diagnostics(source), request);
    if (request.level == DiagnosticLevel::summary) {
      record.metrics = {
          metric("rank_max", "count", as_double(maximum[0])),
          metric("rank_min", "count", as_double(minimum[0])),
          metric("size", "count", as_double(maximum[1])),
          metric("thread_level", "count", as_double(maximum[2]))};
    } else {
      record.counters = {
          counter("collective_calls", global_counters[0]),
          counter("logical_payload_bytes", global_counters[1], "byte"),
          counter("reduced_scalars", global_counters[2])};
    }
    record.state_fingerprint = combined_fingerprint.finish();
    validate(record, describe_diagnostics(source), request);
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  }
  status =
      runtime::collective_status(collective_context, local_ok, local_message);
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input,
        "diagnostics.collective.record", status.failing_rank,
        status.message);

  local_ok = true;
  local_message.clear();
  try {
    sink.submit(record);
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  }
  status =
      runtime::collective_status(collective_context, local_ok, local_message);
  if (!status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::sink_failure,
        "diagnostics.sink.submit", status.failing_rank, status.message);
}

#define HUNDUN_SIMPLE_SUMMARY_ADAPTER(Type, Kind, Id, Fields, Fingerprint,  \
                                      Body)                                 \
  DiagnosticDescriptor describe_diagnostics(const Type&) noexcept {         \
    return descriptor(DiagnosticModuleKind::Kind, Id, kSummary);            \
  }                                                                         \
  std::vector<std::string_view> diagnostic_fingerprint_field_ids(           \
      const Type&) {                                                        \
    return Fields;                                                          \
  }                                                                         \
  void collect_diagnostics(const Type& source,                              \
                           const DiagnosticRequest& request,                \
                           DiagnosticSink& sink) {                           \
    collect_common(describe_diagnostics(source), request, sink,             \
                   diagnostic_fingerprint_field_ids(source),                \
                   [&](FingerprintBuilder& fp) { Fingerprint; },             \
                   [&](DiagnosticRecord& record) {                           \
                     require_summary(request);                              \
                     record.metrics = Body;                                 \
                   });                                                      \
  }

DiagnosticDescriptor describe_diagnostics(
    const runtime::StructuredDecomposition&) noexcept {
  return descriptor(DiagnosticModuleKind::runtime,
                    "hundun.runtime.structured_decomposition", kSummary);
}
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const runtime::StructuredDecomposition&) {
  return {"decomposition", "periodicity"};
}
void collect_diagnostics(const runtime::StructuredDecomposition& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  const auto global = source.global_extent();
  const auto grid = source.process_grid();
  const auto coordinates = source.process_coordinates();
  const auto box = source.owned_box();
  const auto local = source.local_extent();
  const auto periodic = source.periodic();
  collect_common(
      describe_diagnostics(source), request, sink,
      diagnostic_fingerprint_field_ids(source),
      [&](FingerprintBuilder& fp) {
        for (const auto value :
             {global.x, global.y, global.z, grid.x, grid.y, grid.z,
              coordinates.x, coordinates.y, coordinates.z, box.begin.x,
              box.begin.y, box.begin.z, box.end.x, box.end.y, box.end.z,
              local.x, local.y, local.z})
          fp.add_i64("decomposition", value);
        for (bool value : periodic)
          fp.add_u64("periodicity", value ? 1U : 0U);
      },
      [&](DiagnosticRecord& record) {
        require_summary(request);
        const auto product = [](runtime::Int3 value) {
          return static_cast<std::uint64_t>(value.x) *
                 static_cast<std::uint64_t>(value.y) *
                 static_cast<std::uint64_t>(value.z);
        };
        record.metrics = {
            metric("global_cell_count", "count", as_double(product(global))),
            metric("local_cell_count", "count", as_double(product(local))),
            metric("process_grid_x", "count", as_double(grid.x)),
            metric("process_grid_y", "count", as_double(grid.y)),
            metric("process_grid_z", "count", as_double(grid.z)),
            metric("process_coordinate_x", "count", as_double(coordinates.x)),
            metric("process_coordinate_y", "count", as_double(coordinates.y)),
            metric("process_coordinate_z", "count", as_double(coordinates.z)),
            metric("owned_begin_x", "count", as_double(box.begin.x)),
            metric("owned_begin_y", "count", as_double(box.begin.y)),
            metric("owned_begin_z", "count", as_double(box.begin.z)),
            metric("owned_end_x", "count", as_double(box.end.x)),
            metric("owned_end_y", "count", as_double(box.end.y)),
            metric("owned_end_z", "count", as_double(box.end.z)),
            metric("periodic_x", "count", periodic[0] ? 1.0 : 0.0),
            metric("periodic_y", "count", periodic[1] ? 1.0 : 0.0),
            metric("periodic_z", "count", periodic[2] ? 1.0 : 0.0)};
      });
}

DiagnosticDescriptor
describe_diagnostics(const runtime::ExchangePlan&) noexcept {
  return descriptor(DiagnosticModuleKind::halo, "hundun.runtime.halo",
                    kSummary);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const runtime::ExchangePlan&) {
  return {"ghost_width", "regions"};
}
void collect_diagnostics(const runtime::ExchangePlan& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(
      describe_diagnostics(source), request, sink,
      diagnostic_fingerprint_field_ids(source),
      [&](FingerprintBuilder& fp) {
        fp.add_i64("ghost_width", source.ghost_width());
        fp.add_u64("regions", source.regions().size());
        for (const auto& region : source.regions()) {
          for (const auto value :
               {region.offset.x, region.offset.y, region.offset.z,
                region.neighbor_rank, region.send_box.begin.x,
                region.send_box.begin.y, region.send_box.begin.z,
                region.send_box.end.x, region.send_box.end.y,
                region.send_box.end.z, region.receive_box.begin.x,
                region.receive_box.begin.y, region.receive_box.begin.z,
                region.receive_box.end.x, region.receive_box.end.y,
                region.receive_box.end.z})
            fp.add_i64("regions", value);
        }
      },
      [&](DiagnosticRecord& record) {
        require_summary(request);
        record.metrics = {
            metric("ghost_width", "count", as_double(source.ghost_width())),
            metric("region_count", "count",
                   as_double(source.regions().size()))};
        for (std::size_t index = 0; index < source.regions().size(); ++index) {
          const auto& region = source.regions()[index];
          const auto prefix = "region." + std::to_string(index) + ".";
          record.metrics.push_back(metric(prefix + "neighbor_rank", "count",
                                          as_double(region.neighbor_rank)));
          const std::array<std::pair<std::string_view, int>, 15> values{{
              {"offset_x", region.offset.x},
              {"offset_y", region.offset.y},
              {"offset_z", region.offset.z},
              {"send_begin_x", region.send_box.begin.x},
              {"send_begin_y", region.send_box.begin.y},
              {"send_begin_z", region.send_box.begin.z},
              {"send_end_x", region.send_box.end.x},
              {"send_end_y", region.send_box.end.y},
              {"send_end_z", region.send_box.end.z},
              {"receive_begin_x", region.receive_box.begin.x},
              {"receive_begin_y", region.receive_box.begin.y},
              {"receive_begin_z", region.receive_box.begin.z},
              {"receive_end_x", region.receive_box.end.x},
              {"receive_end_y", region.receive_box.end.y},
              {"receive_end_z", region.receive_box.end.z},
          }};
          for (const auto& [name, value] : values)
            record.metrics.push_back(
                metric(prefix + std::string(name), "count", as_double(value)));
        }
      });
}

DiagnosticDescriptor
describe_diagnostics(const mesh::MeshTopology&) noexcept {
  return descriptor(DiagnosticModuleKind::mesh_topology,
                    "hundun.mesh.topology", kSummary);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const mesh::MeshTopology&) {
  return {"cells.global_id", "cells.ownership", "cells.x", "cells.y",
          "cells.z", "faces.axis", "faces.global_id", "faces.neighbour",
          "faces.neighbour.present", "faces.owner", "faces.ownership",
          "faces.pair", "faces.pair.present", "faces.patch",
          "faces.patch.present", "faces.x", "faces.y", "faces.z",
          "mesh_metadata", "patches.face", "patches.id", "patches.name",
          "patches.name.length", "patches.paired",
          "patches.paired.present", "patches.pairing"};
}
void collect_diagnostics(const mesh::MeshTopology& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(
      describe_diagnostics(source), request, sink,
      diagnostic_fingerprint_field_ids(source),
      [&](FingerprintBuilder& fp) {
        const auto extent = source.global_extent();
        const auto box = source.owned_global_box();
        for (const auto value :
             {extent.x, extent.y, extent.z, box.begin.x, box.begin.y,
              box.begin.z, box.end.x, box.end.y, box.end.z})
          fp.add_i64("mesh_metadata", value);
        for (mesh::LocalCellId cell = 0; cell < source.local_cell_count();
             ++cell) {
          const auto logical = source.global_cell(cell);
          fp.add_u64("cells.global_id", source.global_cell_id(cell));
          fp.add_i64("cells.x", logical.x);
          fp.add_i64("cells.y", logical.y);
          fp.add_i64("cells.z", logical.z);
          fp.add_u64("cells.ownership",
                     static_cast<std::uint64_t>(source.cell_ownership(cell)));
        }
        for (mesh::LocalFaceId face = 0; face < source.local_face_count();
             ++face) {
          const auto logical = source.logical_face(face);
          fp.add_u64("faces.global_id", source.global_face_id(face));
          fp.add_u64("faces.axis", static_cast<std::uint64_t>(logical.axis));
          fp.add_i64("faces.x", logical.coordinate.x);
          fp.add_i64("faces.y", logical.coordinate.y);
          fp.add_i64("faces.z", logical.coordinate.z);
          fp.add_u64("faces.ownership",
                     static_cast<std::uint64_t>(
                         source.face_ownership(face)));
          fp.add_u64("faces.owner",
                     source.global_cell_id(source.owner(face)));
          const auto neighbour = source.neighbour(face);
          fp.add_u64("faces.neighbour.present", neighbour ? 1U : 0U);
          fp.add_u64("faces.neighbour",
                     neighbour ? source.global_cell_id(*neighbour)
                               : std::numeric_limits<std::uint64_t>::max());
          const auto patch = source.patch_id(face);
          fp.add_u64("faces.patch.present", patch ? 1U : 0U);
          fp.add_u64("faces.patch",
                     patch.value_or(
                         std::numeric_limits<std::uint32_t>::max()));
          const auto pair = source.periodic_pair(face);
          fp.add_u64("faces.pair.present", pair ? 1U : 0U);
          fp.add_u64("faces.pair",
                     pair.value_or(
                         std::numeric_limits<std::uint64_t>::max()));
        }
        for (const auto& patch : source.patches()) {
          fp.add_u64("patches.id", patch.stable_id());
          fp.add_text("patches.name", "patches.name.length", patch.name());
          fp.add_u64("patches.pairing",
                     static_cast<std::uint64_t>(patch.pairing_kind()));
          const auto paired = patch.paired_patch_id();
          fp.add_u64("patches.paired.present", paired ? 1U : 0U);
          fp.add_u64("patches.paired",
                     paired.value_or(
                         std::numeric_limits<std::uint32_t>::max()));
          fp.add_u64("patches.face", patch.local_faces().size());
          for (const auto face : patch.local_faces())
            fp.add_u64("patches.face", source.global_face_id(face));
        }
      },
      [&](DiagnosticRecord& record) {
        require_summary(request);
        record.metrics = {
            metric("global_cell_count", "count",
                   as_double(source.global_cell_count())),
            metric("global_face_count", "count",
                   as_double(source.global_face_count())),
            metric("local_cell_count", "count",
                   as_double(source.local_cell_count())),
            metric("local_face_count", "count",
                   as_double(source.local_face_count())),
            metric("owned_cell_count", "count",
                   as_double(source.owned_cell_count())),
            metric("ghost_cell_count", "count",
                   as_double(source.ghost_cell_count())),
            metric("owned_face_count", "count",
                   as_double(source.owned_face_count())),
            metric("ghost_face_count", "count",
                   as_double(source.ghost_face_count()))};
        for (const auto& patch : source.patches()) {
          record.identities.push_back(
              {"patch." + std::to_string(patch.stable_id()),
               std::string(patch.name()), std::nullopt, std::nullopt,
               std::nullopt});
          record.metrics.push_back(
              metric("patch." + std::to_string(patch.stable_id()) +
                         ".face_count",
                     "count", as_double(patch.local_faces().size())));
          record.metrics.push_back(
              metric("patch." + std::to_string(patch.stable_id()) +
                         ".pairing_kind",
                     "count", as_double(static_cast<int>(
                                  patch.pairing_kind()))));
        }
      });
}

DiagnosticDescriptor
describe_diagnostics(const FieldLayoutDiagnosticSource&) noexcept {
  return descriptor(DiagnosticModuleKind::field, "hundun.runtime.field_layout",
                    kSummary);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const FieldLayoutDiagnosticSource&) {
  return {"field.components", "field.conservative", "field.ghost_width",
          "field.id", "field.name", "field.name.length", "field.output",
          "field.owner", "field.owner.length", "field.restart",
          "field.scalar_type", "field.space", "field.unit",
          "field.unit.length", "layout.cell_extent", "layout.face_count",
          "role.field", "role.name", "role.name.length"};
}
void collect_diagnostics(const FieldLayoutDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  if (source.registry == nullptr || !source.registry->frozen())
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input, "field-layout.invalid", -1,
        "field layout source is not frozen");
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fp) {
                   for (const auto value :
                        {source.layout.cell_interior_extent.x,
                         source.layout.cell_interior_extent.y,
                         source.layout.cell_interior_extent.z})
                     fp.add_i64("layout.cell_extent", value);
                   fp.add_u64("layout.face_count", source.layout.face_count);
                   for (runtime::FieldId id = 0U;
                        id < source.registry->size(); ++id) {
                     const auto& field = source.registry->descriptor(id);
                     fp.add_u64("field.id", id);
                     fp.add_text("field.name", "field.name.length",
                                 field.name);
                     fp.add_text("field.unit", "field.unit.length",
                                 field.unit);
                     fp.add_text("field.owner", "field.owner.length",
                                 field.owner);
                     fp.add_u64("field.space",
                                static_cast<std::uint64_t>(field.space));
                     fp.add_u64("field.scalar_type",
                                static_cast<std::uint64_t>(field.scalar_type));
                     fp.add_u64("field.components", field.components);
                     fp.add_i64("field.ghost_width", field.ghost_width);
                     fp.add_u64("field.conservative",
                                field.conservative ? 1U : 0U);
                     fp.add_u64("field.restart",
                                static_cast<std::uint64_t>(field.restart));
                     fp.add_u64("field.output",
                                static_cast<std::uint64_t>(field.output));
                   }
                   for (const auto& role : source.roles) {
                     fp.add_text("role.name", "role.name.length", role.role);
                     fp.add_u64("role.field", role.field);
                   }
                 },
                 [&](DiagnosticRecord& record) {
                   require_summary(request);
                   record.metrics = {
                       metric("cell_count", "count",
                              as_double(static_cast<std::uint64_t>(
                                  source.layout.cell_interior_extent.x) *
                                        static_cast<std::uint64_t>(
                                            source.layout.cell_interior_extent.y) *
                                        static_cast<std::uint64_t>(
                                            source.layout.cell_interior_extent.z))),
                       metric("face_count", "count",
                              as_double(source.layout.face_count)),
                       metric("field_count", "count",
                              as_double(source.registry->size()))};
                   for (runtime::FieldId id = 0U;
                        id < source.registry->size(); ++id) {
                     const auto& field = source.registry->descriptor(id);
                     const auto prefix = "field." + std::to_string(id) + ".";
                     const auto descriptor_fingerprint =
                         authoritative_fingerprint(
                             {"components", "conservative", "ghost_width",
                              "name", "name.length", "output", "owner",
                              "owner.length", "restart", "scalar_type",
                              "space", "unit", "unit.length"},
                             [&](FingerprintBuilder& descriptor_fp) {
                               descriptor_fp.add_text("name", "name.length",
                                                      field.name);
                               descriptor_fp.add_text("unit", "unit.length",
                                                      field.unit);
                               descriptor_fp.add_text("owner", "owner.length",
                                                      field.owner);
                               descriptor_fp.add_u64(
                                   "space",
                                   static_cast<std::uint64_t>(field.space));
                               descriptor_fp.add_u64(
                                   "scalar_type",
                                   static_cast<std::uint64_t>(
                                       field.scalar_type));
                               descriptor_fp.add_u64("components",
                                                     field.components);
                               descriptor_fp.add_i64("ghost_width",
                                                     field.ghost_width);
                               descriptor_fp.add_u64(
                                   "conservative",
                                   field.conservative ? 1U : 0U);
                               descriptor_fp.add_u64(
                                   "restart",
                                   static_cast<std::uint64_t>(
                                       field.restart));
                               descriptor_fp.add_u64(
                                   "output",
                                   static_cast<std::uint64_t>(
                                       field.output));
                             });
                     record.identities.push_back(
                         {prefix + "descriptor",
                          descriptor_fingerprint.hex,
                          std::nullopt, std::nullopt, std::nullopt});
                     record.metrics.push_back(
                         metric(prefix + "space", "count",
                                as_double(static_cast<int>(field.space))));
                     record.metrics.push_back(
                         metric(prefix + "scalar_type", "count",
                                as_double(
                                    static_cast<int>(field.scalar_type))));
                     record.metrics.push_back(
                         metric(prefix + "components", "count",
                                as_double(field.components)));
                     record.metrics.push_back(
                         metric(prefix + "ghost_width", "count",
                                as_double(field.ghost_width)));
                     record.metrics.push_back(
                         metric(prefix + "conservative", "count",
                                field.conservative ? 1.0 : 0.0));
                     record.metrics.push_back(
                         metric(prefix + "restart", "count",
                                as_double(static_cast<int>(field.restart))));
                     record.metrics.push_back(
                         metric(prefix + "output", "count",
                                as_double(static_cast<int>(field.output))));
                   }
                   for (const auto& role : source.roles) {
                     if (role.role.empty() ||
                         role.field >= source.registry->size())
                       throw DiagnosticCollectionError(
                           DiagnosticFailureClass::invalid_input,
                           "field-layout.role", -1,
                           "field layout role is invalid");
                     record.identities.push_back(
                         {"role." + role.role, std::nullopt, role.field,
                          std::nullopt, std::nullopt});
                   }
                 });
}

DiagnosticDescriptor
describe_diagnostics(const mesh::MeshGeometry&) noexcept {
  return descriptor(DiagnosticModuleKind::mesh_geometry,
                    "hundun.mesh.geometry", kSummaryInvariants);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const mesh::MeshGeometry&) {
  return {"cells", "faces", "faces.local_count", "mapping.extent",
          "mapping.kind", "mapping.metric"};
}
void collect_diagnostics(const mesh::MeshGeometry& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fp) {
                   const auto extent = source.global_extent();
                   const auto box = source.owned_global_box();
                   const auto origin = source.origin_m();
                   const auto length = source.length_m();
                   fp.add_u64("mapping.kind",
                              static_cast<std::uint64_t>(
                                  source.mapping_kind()));
                   for (const auto value :
                        {extent.x, extent.y, extent.z, box.begin.x,
                         box.begin.y, box.begin.z, box.end.x, box.end.y,
                         box.end.z})
                     fp.add_i64("mapping.extent", value);
                   for (const auto value :
                        {origin.x, origin.y, origin.z, length.x, length.y,
                         length.z})
                     fp.add_fp64("mapping.metric", value);
                   const auto cell_count =
                       static_cast<std::size_t>(box.end.x - box.begin.x) *
                       static_cast<std::size_t>(box.end.y - box.begin.y) *
                       static_cast<std::size_t>(box.end.z - box.begin.z);
                   for (std::size_t cell = 0; cell < cell_count; ++cell) {
                     const auto centre = source.cell_center_m(cell);
                     const auto closure = source.cell_closure_m2(cell);
                     for (const auto value :
                          {centre.x, centre.y, centre.z,
                           source.cell_volume_m3(cell),
                           source.minimum_jacobian_determinant_m3(cell),
                           closure.x, closure.y, closure.z})
                       fp.add_fp64("cells", value);
                   }
                   fp.add_u64("faces.local_count",
                              source.local_face_count());
                   for (std::size_t face = 0;
                        face < source.local_face_count(); ++face) {
                     const auto centre = source.face_center_m(face);
                     const auto displacement =
                         source.face_displacement_m(face);
                     const auto owner = source.face_area_vector_m2(
                         face, mesh::FaceSide::owner);
                     for (const auto value :
                          {centre.x, centre.y, centre.z, displacement.x,
                           displacement.y, displacement.z, owner.x, owner.y,
                           owner.z, source.face_area_m2(face),
                           source.face_skewness(face),
                           source.face_non_orthogonality_degrees(face)})
                       fp.add_fp64("faces", value);
                   }
                 },
                 [&](DiagnosticRecord& record) {
                   const auto box = source.owned_global_box();
                   const std::size_t count =
                       static_cast<std::size_t>(box.end.x - box.begin.x) *
                       static_cast<std::size_t>(box.end.y - box.begin.y) *
                       static_cast<std::size_t>(box.end.z - box.begin.z);
                   double volume_sum = 0.0;
                   double minimum_jacobian =
                       std::numeric_limits<double>::infinity();
                   double maximum_closure = 0.0;
                   double maximum_face_area = 0.0;
                   double maximum_skewness = 0.0;
                   double maximum_non_orthogonality = 0.0;
                   for (std::size_t cell = 0; cell < count; ++cell) {
                     volume_sum += source.cell_volume_m3(cell);
                     minimum_jacobian =
                         std::min(minimum_jacobian,
                                  source.minimum_jacobian_determinant_m3(cell));
                     const auto closure = source.cell_closure_m2(cell);
                     maximum_closure =
                         std::max(maximum_closure,
                                  std::hypot(closure.x,
                                             std::hypot(closure.y, closure.z)));
                   }
                   for (std::size_t face = 0;
                        face < source.local_face_count(); ++face) {
                     maximum_face_area =
                         std::max(maximum_face_area,
                                  source.face_area_m2(face));
                     maximum_skewness =
                         std::max(maximum_skewness,
                                  source.face_skewness(face));
                     maximum_non_orthogonality = std::max(
                         maximum_non_orthogonality,
                         source.face_non_orthogonality_degrees(face));
                   }
                   const auto origin = source.origin_m();
                   const auto length = source.length_m();
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("mapping_kind", "count",
                                as_double(static_cast<int>(
                                    source.mapping_kind()))),
                         metric("origin_x", "m", origin.x),
                         metric("origin_y", "m", origin.y),
                         metric("origin_z", "m", origin.z),
                         metric("length_x", "m", length.x),
                         metric("length_y", "m", length.y),
                         metric("length_z", "m", length.z),
                         metric("local_face_count", "count",
                                as_double(source.local_face_count())),
                         metric("cell_volume_sum", "m3", volume_sum),
                         metric("maximum_closure_norm", "m2",
                                maximum_closure),
                         metric("minimum_jacobian", "m3",
                                minimum_jacobian),
                         metric("maximum_face_area", "m2",
                                maximum_face_area),
                         metric("maximum_skewness", "1",
                                maximum_skewness),
                         metric("maximum_non_orthogonality", "degree",
                                maximum_non_orthogonality)};
                   } else {
                     record.invariants = {
                         finite_invariant("cell_closure_finite", "m2",
                                          maximum_closure),
                         positive_invariant("cell_volume_sum_positive", "m3",
                                            volume_sum),
                         positive_invariant("jacobian_positive", "m3",
                                            minimum_jacobian)};
                   }
                 });
}

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    execution::ExecutionContext, execution, "hundun.execution.context",
    (std::vector<std::string_view>{"backend_identity", "backend_name",
                                   "backend_name.length", "capabilities",
                                   "ordered", "space"}),
    ([&] {
      fp.add_text("backend_name", "backend_name.length",
                  source.backend_name());
      fp.add_u64("backend_identity", source.backend_identity());
      fp.add_u64("space", static_cast<std::uint64_t>(source.space()));
      fp.add_u64("ordered", source.ordered() ? 1U : 0U);
      for (const auto capability :
           {execution::ExecutionCapability::buffer_allocation,
            execution::ExecutionCapability::host_access,
            execution::ExecutionCapability::transfer,
            execution::ExecutionCapability::asynchronous_event})
        fp.add_u64("capabilities",
                   source.supports(capability) ? 1U : 0U);
    }()),
    (std::vector<DiagnosticMetric>{
        metric("backend_identity", "count",
               as_double(source.backend_identity())),
        metric("ordered", "count", source.ordered() ? 1.0 : 0.0),
        metric("space", "count",
               source.space() == execution::ExecutionSpace::host ? 0.0
                                                                 : 1.0),
        metric("capability_buffer_allocation", "count",
               source.supports(
                   execution::ExecutionCapability::buffer_allocation)
                   ? 1.0
                   : 0.0),
        metric("capability_host_access", "count",
               source.supports(execution::ExecutionCapability::host_access)
                   ? 1.0
                   : 0.0),
        metric("capability_transfer", "count",
               source.supports(execution::ExecutionCapability::transfer)
                   ? 1.0
                   : 0.0),
        metric("capability_asynchronous_event", "count",
               source.supports(
                   execution::ExecutionCapability::asynchronous_event)
                   ? 1.0
                   : 0.0)}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    execution::Buffer, execution, "hundun.execution.buffer",
    (std::vector<std::string_view>{"allocation_identity", "backend_identity",
                                   "byte_size", "epoch", "space"}),
    fp.add_u64("allocation_identity", source.allocation_identity());
    fp.add_u64("byte_size", source.byte_size());
    fp.add_u64("epoch", source.epoch());
    fp.add_u64("backend_identity", source.backend_identity());
    fp.add_u64("space", static_cast<std::uint64_t>(source.space())),
    (std::vector<DiagnosticMetric>{
        metric("allocation_identity", "count",
               as_double(source.allocation_identity())),
        metric("byte_size", "byte", as_double(source.byte_size())),
        metric("epoch", "count", as_double(source.epoch())),
        metric("backend_identity", "count",
               as_double(source.backend_identity())),
        metric("space", "count",
               as_double(static_cast<int>(source.space())))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    execution::VectorView<const double>, execution,
    "hundun.execution.vector_view",
    (std::vector<std::string_view>{"allocation_identity", "backend_identity",
                                   "element_count", "epoch", "offset_bytes",
                                   "scalar_format", "space", "stride",
                                   "writable"}),
    fp.add_u64("allocation_identity", source.allocation_identity());
    fp.add_u64("element_count", source.size());
    fp.add_u64("epoch", source.epoch());
    fp.add_u64("offset_bytes", source.offset_bytes());
    fp.add_u64("stride", source.stride());
    fp.add_u64("backend_identity", source.backend_identity());
    fp.add_u64("space", static_cast<std::uint64_t>(source.space()));
    fp.add_u64("scalar_format",
               static_cast<std::uint64_t>(source.scalar_format()));
    fp.add_u64("writable", source.writable() ? 1U : 0U),
    (std::vector<DiagnosticMetric>{
        metric("allocation_identity", "count",
               as_double(source.allocation_identity())),
        metric("element_count", "count", as_double(source.size())),
        metric("epoch", "count", as_double(source.epoch())),
        metric("offset_bytes", "byte", as_double(source.offset_bytes())),
        metric("stride", "count", as_double(source.stride())),
        metric("backend_identity", "count",
               as_double(source.backend_identity())),
        metric("space", "count",
               as_double(static_cast<int>(source.space()))),
        metric("scalar_format", "count",
               as_double(static_cast<int>(source.scalar_format()))),
        metric("writable", "count", source.writable() ? 1.0 : 0.0)}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    linear::GhostedVector, execution, "hundun.linear.ghosted_vector",
    (std::vector<std::string_view>{"allocation_identity", "backend_identity",
                                   "epoch", "ghost_count", "layout.global_id",
                                   "local_count", "owned_count", "space"}),
    fp.add_u64("owned_count", source.owned_count());
    fp.add_u64("ghost_count", source.ghost_count());
    fp.add_u64("local_count", source.local_count());
    for (const auto id : source.layout().global_ids())
      fp.add_u64("layout.global_id", id);
    fp.add_u64("allocation_identity", source.allocation_identity());
    fp.add_u64("epoch", source.epoch());
    fp.add_u64("backend_identity", source.backend_identity());
    fp.add_u64("space", static_cast<std::uint64_t>(source.space())),
    (std::vector<DiagnosticMetric>{
        metric("allocation_identity", "count",
               as_double(source.allocation_identity())),
        metric("epoch", "count", as_double(source.epoch())),
        metric("ghost_count", "count", as_double(source.ghost_count())),
        metric("local_count", "count", as_double(source.local_count())),
        metric("owned_count", "count", as_double(source.owned_count())),
        metric("backend_identity", "count",
               as_double(source.backend_identity())),
        metric("space", "count",
               as_double(static_cast<int>(source.space())))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    linear::GhostedVectorHalo, halo, "hundun.linear.ghosted_vector_halo",
    (std::vector<std::string_view>{"ghost_count", "owned_count", "path",
                                   "receive_value_count",
                                   "send_value_count"}),
    fp.add_u64("path", static_cast<std::uint64_t>(source.path()));
    fp.add_u64("owned_count", source.owned_count());
    fp.add_u64("ghost_count", source.ghost_count());
    fp.add_u64("send_value_count", source.send_value_count());
    fp.add_u64("receive_value_count", source.receive_value_count()),
    (std::vector<DiagnosticMetric>{
        metric("ghost_count", "count", as_double(source.ghost_count())),
        metric("owned_count", "count", as_double(source.owned_count())),
        metric("path", "count", as_double(static_cast<int>(source.path()))),
        metric("receive_value_count", "count",
               as_double(source.receive_value_count())),
        metric("send_value_count", "count",
               as_double(source.send_value_count()))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    linear::LinearOperator, linear_operator, "hundun.linear.operator",
    (std::vector<std::string_view>{
        "context.backend_identity", "context.backend_name",
        "context.backend_name.length", "context.space", "diagonal_available",
        "domain.ghost", "domain.global_id", "domain.local", "domain.owned",
        "range.ghost", "range.global_id", "range.local", "range.owned",
        "revision"}),
    ([&] {
      const auto domain = source.domain_layout();
      const auto range = source.range_layout();
      fp.add_u64("domain.owned", domain.owned_count());
      fp.add_u64("domain.ghost", domain.ghost_count());
      fp.add_u64("domain.local", domain.local_count());
      for (const auto id : domain.global_ids())
        fp.add_u64("domain.global_id", id);
      fp.add_u64("range.owned", range.owned_count());
      fp.add_u64("range.ghost", range.ghost_count());
      fp.add_u64("range.local", range.local_count());
      for (const auto id : range.global_ids())
        fp.add_u64("range.global_id", id);
      fp.add_text("context.backend_name", "context.backend_name.length",
                  source.context().backend_name());
      fp.add_u64("context.backend_identity",
                 source.context().backend_identity());
      fp.add_u64("context.space",
                 static_cast<std::uint64_t>(source.context().space()));
      fp.add_u64("revision", source.revision());
      fp.add_u64("diagonal_available",
                 source.has_diagonal() ? 1U : 0U);
    }()),
    (std::vector<DiagnosticMetric>{
        metric("diagonal_available", "count",
               source.has_diagonal() ? 1.0 : 0.0),
        metric("domain_owned_count", "count",
               as_double(source.domain_layout().owned_count())),
        metric("domain_ghost_count", "count",
               as_double(source.domain_layout().ghost_count())),
        metric("domain_local_count", "count",
               as_double(source.domain_layout().local_count())),
        metric("range_owned_count", "count",
               as_double(source.range_layout().owned_count())),
        metric("range_ghost_count", "count",
               as_double(source.range_layout().ghost_count())),
        metric("range_local_count", "count",
               as_double(source.range_layout().local_count())),
        metric("revision", "count", as_double(source.revision())),
        metric("context_backend_identity", "count",
               as_double(source.context().backend_identity())),
        metric("context_space", "count",
               as_double(static_cast<int>(source.context().space())))}))

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    finite_volume::MatrixFreePoissonOperator, linear_operator,
    "hundun.finite_volume.poisson",
    (std::vector<std::string_view>{
        "constraint_mode", "context.backend_identity", "context.backend_name",
        "context.backend_name.length", "context.space", "diagonal_available",
        "domain.global_id", "domain.owned", "pressure_reference_patch",
        "pressure_reference_patch.present", "range.global_id", "range.owned",
        "revision", "solver_family"}),
    ([&] {
      const auto domain = source.domain_layout();
      const auto range = source.range_layout();
      for (const auto id : domain.global_ids())
        fp.add_u64("domain.global_id", id);
      for (const auto id : range.global_ids())
        fp.add_u64("range.global_id", id);
      fp.add_u64("domain.owned", domain.owned_count());
      fp.add_u64("range.owned", range.owned_count());
      fp.add_text("context.backend_name", "context.backend_name.length",
                  source.context().backend_name());
      fp.add_u64("context.backend_identity",
                 source.context().backend_identity());
      fp.add_u64("context.space",
                 static_cast<std::uint64_t>(source.context().space()));
      fp.add_u64("revision", source.revision());
      fp.add_u64("diagonal_available",
                 source.has_diagonal() ? 1U : 0U);
      fp.add_u64("constraint_mode",
                 static_cast<std::uint64_t>(source.constraint_mode()));
      const auto reference = source.pressure_reference_patch_id();
      fp.add_u64("pressure_reference_patch.present",
                 reference ? 1U : 0U);
      fp.add_u64("pressure_reference_patch",
                 reference.value_or(
                     std::numeric_limits<std::uint32_t>::max()));
      fp.add_u64("solver_family",
                 static_cast<std::uint64_t>(source.solver_family()));
    }()),
    (std::vector<DiagnosticMetric>{
        metric("constraint_mode", "count",
               as_double(static_cast<int>(source.constraint_mode()))),
        metric("diagonal_available", "count",
               source.has_diagonal() ? 1.0 : 0.0),
        metric("revision", "count", as_double(source.revision())),
        metric("solver_family", "count",
               as_double(static_cast<int>(source.solver_family()))),
        metric("pressure_reference_present", "count",
               source.pressure_reference_patch_id() ? 1.0 : 0.0),
        metric("pressure_reference_patch", "count",
               as_double(source.pressure_reference_patch_id().value_or(
                   0U)))}))

DiagnosticDescriptor
describe_diagnostics(const LinearSolveDiagnosticSource& source) noexcept {
  return {kDiagnosticRecordSchemaV1,
          DiagnosticModuleKind::linear_solver,
          "hundun.linear.solve",
          source.instance_id,
          kSummaryCounters};
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const LinearSolveDiagnosticSource&) {
  return {"final_residual", "global_reduction_count", "initial_residual",
          "instance_id", "instance_id.length", "iterations",
          "lowest_failing_rank", "matvec_count",
          "preconditioner_apply_count", "recursive_residual", "termination"};
}
void collect_diagnostics(const LinearSolveDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  if (source.report == nullptr || source.instance_id.empty())
    throw DiagnosticCollectionError(DiagnosticFailureClass::invalid_input,
                                    "linear-solve.invalid", -1,
                                    "linear solve source is invalid");
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fp) {
                   fp.add_text("instance_id", "instance_id.length",
                               source.instance_id);
                   fp.add_u64("termination",
                              static_cast<std::uint64_t>(
                                  source.report->reason));
                   fp.add_u64("iterations", source.report->iterations);
                   fp.add_fp64("initial_residual",
                               source.report->initial_residual);
                   fp.add_fp64("recursive_residual",
                               source.report->recursive_residual);
                   fp.add_fp64("final_residual",
                               source.report->final_residual);
                   fp.add_u64("matvec_count",
                              source.report->matvec_count);
                   fp.add_u64("preconditioner_apply_count",
                              source.report->preconditioner_apply_count);
                   fp.add_u64("global_reduction_count",
                              source.report->global_reduction_count);
                   fp.add_i64("lowest_failing_rank",
                              source.report->lowest_failing_rank);
                 },
                 [&](DiagnosticRecord& record) {
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("iterations", "count",
                                as_double(source.report->iterations)),
                         metric("termination", "count",
                                as_double(static_cast<int>(
                                    source.report->reason))),
                         metric("lowest_failing_rank", "count",
                                as_double(
                                    source.report->lowest_failing_rank))};
                   } else {
                     record.counters = {
                         counter("global_reduction_count",
                                 source.report->global_reduction_count),
                         counter("matvec_count",
                                 source.report->matvec_count),
                         counter("preconditioner_apply_count",
                                 source.report->preconditioner_apply_count)};
                   }
                 });
}

HUNDUN_SIMPLE_SUMMARY_ADAPTER(
    SharedFluxDiagnosticSource, finite_volume,
    "hundun.finite_volume.shared_flux",
    (std::vector<std::string_view>{"face_count", "field_id", "final_flux"}),
    fp.add_u64("field_id", source.field);
    fp.add_u64("face_count", source.face_count);
    fp.add_u64("final_flux", source.final_flux ? 1U : 0U),
    (std::vector<DiagnosticMetric>{
        metric("face_count", "count", as_double(source.face_count)),
        metric("field_id", "count",
               as_double(static_cast<std::uint64_t>(source.field))),
        metric("final_flux", "count", source.final_flux ? 1.0 : 0.0)}))

DiagnosticDescriptor
describe_diagnostics(const boundary::BoundaryRegistry&) noexcept {
  return descriptor(DiagnosticModuleKind::boundary,
                    "hundun.boundary.registry", kSummaryCounters);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const boundary::BoundaryRegistry& source) {
  std::vector<std::string_view> ids{
      "inlet_patch.present", "open_domain", "outlet_patch.present",
      "patch.density_rule", "patch.enthalpy_rule", "patch.id",
      "patch.inlet.present", "patch.kind", "patch.mass_flux_rule",
      "patch.name", "patch.name.length", "patch.paired.present",
      "patch.pressure.present", "patch.pressure_rule", "patch.scalar_rule",
      "patch.velocity_rule"};
  if (source.scalar_count() != 0U)
    ids.insert(ids.end(), {"scalar", "scalar.length"});
  if (source.velocity_inlet_patch_id())
    ids.push_back("inlet_patch");
  if (source.pressure_outlet_patch_id())
    ids.push_back("outlet_patch");
  bool has_pair = false;
  bool has_inlet = false;
  bool has_inlet_scalar = false;
  bool has_inlet_temperature = false;
  bool has_pressure = false;
  for (std::uint32_t patch = 0; patch < 6U; ++patch) {
    const auto& descriptor = source.patch(patch);
    has_pair = has_pair || descriptor.paired_patch_id().has_value();
    has_pressure =
        has_pressure || descriptor.pressure_value_pa().has_value();
    if (const auto& inlet = descriptor.inlet_state()) {
      has_inlet = true;
      has_inlet_scalar =
          has_inlet_scalar || !inlet->scalar_values.empty();
      has_inlet_temperature =
          has_inlet_temperature || inlet->temperature_K.has_value();
    }
  }
  if (has_pair)
    ids.push_back("patch.paired");
  if (has_inlet)
    ids.insert(ids.end(),
               {"patch.inlet", "patch.inlet.temperature.present"});
  if (has_inlet_scalar)
    ids.push_back("patch.inlet.scalar");
  if (has_inlet_temperature)
    ids.push_back("patch.inlet.temperature");
  if (has_pressure)
    ids.push_back("patch.pressure");
  std::sort(ids.begin(), ids.end());
  return ids;
}
void collect_diagnostics(const boundary::BoundaryRegistry& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fp) {
                   fp.add_u64("open_domain",
                              source.open_domain() ? 1U : 0U);
                   for (std::size_t scalar = 0;
                        scalar < source.scalar_count(); ++scalar)
                     fp.add_text("scalar", "scalar.length",
                                 source.scalar_name(scalar));
                   const auto inlet_patch =
                       source.velocity_inlet_patch_id();
                   const auto outlet_patch =
                       source.pressure_outlet_patch_id();
                   fp.add_u64("inlet_patch.present",
                              inlet_patch ? 1U : 0U);
                   if (inlet_patch)
                     fp.add_u64("inlet_patch", *inlet_patch);
                   fp.add_u64("outlet_patch.present",
                              outlet_patch ? 1U : 0U);
                   if (outlet_patch)
                     fp.add_u64("outlet_patch", *outlet_patch);
                   for (std::uint32_t patch = 0; patch < 6U; ++patch) {
                     const auto& descriptor = source.patch(patch);
                     fp.add_u64("patch.id", descriptor.stable_id());
                     fp.add_text("patch.name", "patch.name.length",
                                 descriptor.name());
                     fp.add_u64("patch.kind",
                                static_cast<std::uint64_t>(
                                    descriptor.kind()));
                     fp.add_u64("patch.velocity_rule",
                                static_cast<std::uint64_t>(
                                    descriptor.velocity_rule()));
                     fp.add_u64("patch.pressure_rule",
                                static_cast<std::uint64_t>(
                                    descriptor.pressure_rule()));
                     fp.add_u64("patch.density_rule",
                                static_cast<std::uint64_t>(
                                    descriptor.density_rule()));
                     fp.add_u64("patch.enthalpy_rule",
                                static_cast<std::uint64_t>(
                                    descriptor.enthalpy_rule()));
                     fp.add_u64("patch.scalar_rule",
                                static_cast<std::uint64_t>(
                                    descriptor.scalar_rule()));
                     fp.add_u64("patch.mass_flux_rule",
                                static_cast<std::uint64_t>(
                                    descriptor.mass_flux_rule()));
                     const auto paired = descriptor.paired_patch_id();
                     fp.add_u64("patch.paired.present",
                                paired ? 1U : 0U);
                     if (paired)
                       fp.add_u64("patch.paired", *paired);
                     const auto& inlet = descriptor.inlet_state();
                     fp.add_u64("patch.inlet.present",
                                inlet ? 1U : 0U);
                     if (inlet) {
                       for (const auto value :
                            {inlet->velocity_m_per_s.x,
                             inlet->velocity_m_per_s.y,
                             inlet->velocity_m_per_s.z,
                             inlet->density_kg_per_m3,
                             inlet->enthalpy_J_per_kg})
                         fp.add_fp64("patch.inlet", value);
                       fp.add_u64("patch.inlet.temperature.present",
                                  inlet->temperature_K ? 1U : 0U);
                       if (inlet->temperature_K)
                         fp.add_fp64("patch.inlet.temperature",
                                     *inlet->temperature_K);
                       for (double scalar : inlet->scalar_values)
                         fp.add_fp64("patch.inlet.scalar", scalar);
                     }
                     const auto pressure =
                         descriptor.pressure_value_pa();
                     fp.add_u64("patch.pressure.present",
                                pressure ? 1U : 0U);
                     if (pressure)
                       fp.add_fp64("patch.pressure", *pressure);
                   }
                 },
                 [&](DiagnosticRecord& record) {
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("open_domain", "count",
                                source.open_domain() ? 1.0 : 0.0),
                         metric("scalar_count", "count",
                                as_double(source.scalar_count())),
                         metric("velocity_inlet_present", "count",
                                source.velocity_inlet_patch_id() ? 1.0
                                                                 : 0.0),
                         metric("velocity_inlet_patch", "count",
                                as_double(
                                    source.velocity_inlet_patch_id()
                                        .value_or(0U))),
                         metric("pressure_outlet_present", "count",
                                source.pressure_outlet_patch_id() ? 1.0
                                                                  : 0.0),
                         metric("pressure_outlet_patch", "count",
                                as_double(
                                    source.pressure_outlet_patch_id()
                                        .value_or(0U)))};
                     for (std::size_t scalar = 0;
                          scalar < source.scalar_count(); ++scalar)
                       record.identities.push_back(
                           {"scalar." + std::to_string(scalar),
                            std::string(source.scalar_name(scalar)),
                            std::nullopt, std::nullopt, std::nullopt});
                     for (std::uint32_t patch = 0; patch < 6U; ++patch) {
                       const auto& descriptor = source.patch(patch);
                       const auto prefix =
                           "patch." + std::to_string(patch) + ".";
                       record.identities.push_back(
                           {prefix + "descriptor",
                            std::string(descriptor.name()), std::nullopt,
                            std::nullopt, std::nullopt});
                       record.metrics.push_back(metric(
                           prefix + "kind", "count",
                           as_double(static_cast<int>(descriptor.kind()))));
                       record.metrics.push_back(metric(
                           prefix + "velocity_rule", "count",
                           as_double(static_cast<int>(
                               descriptor.velocity_rule()))));
                       record.metrics.push_back(metric(
                           prefix + "pressure_rule", "count",
                           as_double(static_cast<int>(
                               descriptor.pressure_rule()))));
                       record.metrics.push_back(metric(
                           prefix + "density_rule", "count",
                           as_double(static_cast<int>(
                               descriptor.density_rule()))));
                       record.metrics.push_back(metric(
                           prefix + "enthalpy_rule", "count",
                           as_double(static_cast<int>(
                               descriptor.enthalpy_rule()))));
                       record.metrics.push_back(metric(
                           prefix + "scalar_rule", "count",
                           as_double(static_cast<int>(
                               descriptor.scalar_rule()))));
                       record.metrics.push_back(metric(
                           prefix + "mass_flux_rule", "count",
                           as_double(static_cast<int>(
                               descriptor.mass_flux_rule()))));
                       if (const auto& inlet = descriptor.inlet_state()) {
                         record.metrics.push_back(metric(
                             prefix + "inlet_velocity_x", "m/s",
                             inlet->velocity_m_per_s.x));
                         record.metrics.push_back(metric(
                             prefix + "inlet_velocity_y", "m/s",
                             inlet->velocity_m_per_s.y));
                         record.metrics.push_back(metric(
                             prefix + "inlet_velocity_z", "m/s",
                             inlet->velocity_m_per_s.z));
                         record.metrics.push_back(metric(
                             prefix + "inlet_density", "kg/m3",
                             inlet->density_kg_per_m3));
                         record.metrics.push_back(metric(
                             prefix + "inlet_enthalpy", "J/kg",
                             inlet->enthalpy_J_per_kg));
                         if (inlet->temperature_K)
                           record.metrics.push_back(metric(
                               prefix + "inlet_temperature", "K",
                               *inlet->temperature_K));
                         for (std::size_t scalar = 0;
                              scalar < inlet->scalar_values.size(); ++scalar)
                           record.metrics.push_back(metric(
                               prefix + "inlet_scalar." +
                                   std::to_string(scalar),
                               "1", inlet->scalar_values[scalar]));
                       }
                       if (const auto pressure =
                               descriptor.pressure_value_pa())
                         record.metrics.push_back(metric(
                             prefix + "pressure", "Pa", *pressure));
                     }
                   } else {
                     std::uint64_t periodic = 0;
                     std::uint64_t walls = 0;
                     for (std::uint32_t patch = 0; patch < 6U; ++patch) {
                       const auto kind = source.patch(patch).kind();
                       periodic += kind == boundary::BoundaryKind::periodic;
                       walls += kind == boundary::BoundaryKind::no_slip_wall;
                     }
                     record.counters = {
                         counter("no_slip_wall_count", walls),
                         counter("patch_count", 6U),
                         counter("periodic_patch_count", periodic)};
                   }
                 });
}

DiagnosticDescriptor describe_diagnostics(
    const ConstantDensityPisoDiagnosticSource&) noexcept {
  return descriptor(DiagnosticModuleKind::piso,
                    "hundun.flow.constant_density_piso",
                    kSummaryCounters);
}
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const ConstantDensityPisoDiagnosticSource& source) {
  std::vector<std::string_view> ids{
      "attempted_dt", "backflow.present", "continuity", "disposition",
      "lowest_failing_rank", "mass_conservation", "momentum_conservation",
      "momentum_residual", "pressure", "pressure_corrector_count", "reason",
      "solve.final_residual", "solve.global_reduction_count",
      "solve.initial_residual", "solve.iterations",
      "solve.lowest_failing_rank", "solve.matvec_count",
      "solve.preconditioner_apply_count", "solve.prefix",
      "solve.prefix.length", "solve.reason", "solve.recursive_residual",
      "suggested_dt", "transport_conservation", "transport_residual"};
  if (source.report != nullptr && source.report->final_backflow_evidence) {
    ids.insert(ids.end(),
               {"backflow.face", "backflow.minimum_flux", "backflow.patch",
                "backflow.rank", "backflow.step", "backflow.time"});
    std::sort(ids.begin(), ids.end());
  }
  return ids;
}
void collect_diagnostics(const ConstantDensityPisoDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  if (source.report == nullptr ||
      ((source.report->disposition == flow::StepAttemptDisposition::committed) !=
       (source.report->reason == flow::StepFailureReason::none)) ||
      source.report->final_transport_normalized_l2.size() !=
          source.report->final_transport_relative_conservation_defect.size() ||
      (source.report->disposition ==
           flow::StepAttemptDisposition::committed &&
       source.report->pressure_corrector_count != 2U))
    throw DiagnosticCollectionError(DiagnosticFailureClass::invalid_input,
                                    "constant-piso.invalid", -1,
                                    "constant density PISO source is invalid");
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fp) {
                   const auto& report = *source.report;
                   fp.add_u64("disposition",
                              static_cast<std::uint64_t>(
                                  report.disposition));
                   fp.add_u64("reason",
                              static_cast<std::uint64_t>(report.reason));
                   fp.add_i64("lowest_failing_rank",
                              report.lowest_failing_rank);
                   fp.add_u64("pressure_corrector_count",
                              report.pressure_corrector_count);
                   fp.add_fp64("attempted_dt", report.attempted_dt_s);
                   fp.add_fp64("suggested_dt", report.suggested_dt_s);
                   for (std::size_t component = 0; component < 3U;
                        ++component)
                     add_solve_fingerprint(
                         fp, "momentum." + std::to_string(component),
                         report.momentum.components[component]);
                   for (std::size_t corrector = 0; corrector < 2U;
                        ++corrector)
                     add_solve_fingerprint(
                         fp, "pressure." + std::to_string(corrector),
                         report.pressure[corrector]);
                   fp.add_fp64("continuity",
                               report.final_continuity_normalized_l2);
                   fp.add_fp64("pressure",
                               report.final_pressure_residual_l2);
                   for (double value :
                        report.final_momentum_normalized_l2)
                     fp.add_fp64("momentum_residual", value);
                   for (double value :
                        report.final_transport_normalized_l2)
                     fp.add_fp64("transport_residual", value);
                   fp.add_fp64(
                       "mass_conservation",
                       report.final_mass_relative_conservation_defect);
                   for (double value :
                        report.final_momentum_relative_conservation_defect)
                     fp.add_fp64("momentum_conservation", value);
                   for (double value :
                        report.final_transport_relative_conservation_defect)
                     fp.add_fp64("transport_conservation", value);
                   fp.add_u64("backflow.present",
                              report.final_backflow_evidence ? 1U : 0U);
                   if (report.final_backflow_evidence) {
                     const auto& evidence =
                         *report.final_backflow_evidence;
                     fp.add_u64("backflow.patch", evidence.patch_id);
                     fp.add_u64("backflow.step", evidence.step);
                     fp.add_fp64("backflow.time", evidence.time_s);
                     fp.add_fp64(
                         "backflow.minimum_flux",
                         evidence.minimum_outward_mass_flux_kg_per_s);
                     fp.add_u64("backflow.face",
                                evidence.global_face_id);
                     fp.add_i64("backflow.rank",
                                evidence.lowest_failing_rank);
                   }
                 },
                 [&](DiagnosticRecord& record) {
                   const auto& report = *source.report;
                   set_failure(record, report.reason,
                               report.lowest_failing_rank);
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("disposition", "count",
                                as_double(static_cast<int>(
                                    report.disposition))),
                         metric("reason", "count",
                                as_double(static_cast<int>(report.reason))),
                         metric("lowest_failing_rank", "count",
                                as_double(report.lowest_failing_rank)),
                         metric("attempted_dt", "s",
                                report.attempted_dt_s),
                         metric("suggested_dt", "s",
                                report.suggested_dt_s),
                         metric("continuity_normalized_l2", "1",
                                report.final_continuity_normalized_l2),
                         metric("pressure_residual_l2", "1",
                                report.final_pressure_residual_l2)};
                     for (std::size_t component = 0; component < 3U;
                          ++component) {
                       record.metrics.push_back(metric(
                           "momentum." + std::to_string(component) +
                               ".normalized_l2",
                           "1",
                           report.final_momentum_normalized_l2[component]));
                       record.metrics.push_back(metric(
                           "momentum." + std::to_string(component) +
                               ".conservation_defect",
                           "1",
                           report
                               .final_momentum_relative_conservation_defect
                                   [component]));
                       record.metrics.push_back(metric(
                           "momentum." + std::to_string(component) +
                               ".solve_reason",
                           "count",
                           as_double(static_cast<int>(
                               report.momentum.components[component]
                                   .reason))));
                       record.metrics.push_back(metric(
                           "momentum." + std::to_string(component) +
                               ".solve_iterations",
                           "count",
                           as_double(
                               report.momentum.components[component]
                                   .iterations)));
                     }
                     for (std::size_t field = 0;
                          field <
                          report.final_transport_normalized_l2.size();
                          ++field) {
                       record.metrics.push_back(metric(
                           "transport." + std::to_string(field) +
                               ".normalized_l2",
                           "1",
                           report.final_transport_normalized_l2[field]));
                       record.metrics.push_back(metric(
                           "transport." + std::to_string(field) +
                               ".conservation_defect",
                           "1",
                           report
                               .final_transport_relative_conservation_defect
                                   [field]));
                     }
                     record.metrics.push_back(metric(
                         "mass_conservation_defect", "1",
                         report.final_mass_relative_conservation_defect));
                     for (std::size_t corrector = 0; corrector < 2U;
                          ++corrector) {
                       record.metrics.push_back(metric(
                           "pressure." + std::to_string(corrector) +
                               ".solve_reason",
                           "count",
                           as_double(static_cast<int>(
                               report.pressure[corrector].reason))));
                       record.metrics.push_back(metric(
                           "pressure." + std::to_string(corrector) +
                               ".solve_iterations",
                           "count",
                           as_double(
                               report.pressure[corrector].iterations)));
                     }
                     if (report.final_backflow_evidence) {
                       const auto& evidence =
                           *report.final_backflow_evidence;
                       record.metrics.push_back(metric(
                           "backflow_patch", "count",
                           as_double(evidence.patch_id)));
                       record.metrics.push_back(metric(
                           "backflow_step", "count",
                           as_double(evidence.step)));
                       record.metrics.push_back(metric(
                           "backflow_time", "s", evidence.time_s));
                       record.metrics.push_back(metric(
                           "backflow_minimum_flux", "kg/s",
                           evidence.minimum_outward_mass_flux_kg_per_s));
                       record.metrics.push_back(metric(
                           "backflow_global_face", "count",
                           as_double(evidence.global_face_id)));
                     }
                   } else {
                     record.counters = {
                         counter("pressure_corrector_count",
                                 report.pressure_corrector_count),
                         counter("transported_field_count",
                                 checked_u64(
                                     report.final_transport_normalized_l2
                                         .size()))};
                     for (std::size_t component = 0; component < 3U;
                          ++component) {
                       const auto& solve =
                           report.momentum.components[component];
                       const auto prefix =
                           "momentum." + std::to_string(component) + ".";
                       record.counters.push_back(
                           counter(prefix + "iterations", solve.iterations));
                       record.counters.push_back(
                           counter(prefix + "matvec_count",
                                   solve.matvec_count));
                       record.counters.push_back(counter(
                           prefix + "preconditioner_apply_count",
                           solve.preconditioner_apply_count));
                       record.counters.push_back(counter(
                           prefix + "global_reduction_count",
                           solve.global_reduction_count));
                     }
                     for (std::size_t corrector = 0; corrector < 2U;
                          ++corrector) {
                       const auto& solve = report.pressure[corrector];
                       const auto prefix =
                           "pressure." + std::to_string(corrector) + ".";
                       record.counters.push_back(
                           counter(prefix + "iterations", solve.iterations));
                       record.counters.push_back(
                           counter(prefix + "matvec_count",
                                   solve.matvec_count));
                       record.counters.push_back(counter(
                           prefix + "preconditioner_apply_count",
                           solve.preconditioner_apply_count));
                       record.counters.push_back(counter(
                           prefix + "global_reduction_count",
                           solve.global_reduction_count));
                     }
                   }
                 });
}

DiagnosticDescriptor
describe_diagnostics(const FlowDriverDiagnosticSource&) noexcept {
  return descriptor(DiagnosticModuleKind::flow_driver,
                    "hundun.application.flow_driver", kSummaryCounters);
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const FlowDriverDiagnosticSource&) {
  return {"attempt_count", "density_model", "disposition",
          "lowest_failing_rank", "reason", "step", "time"};
}
void collect_diagnostics(const FlowDriverDiagnosticSource& source,
                         const DiagnosticRequest& request,
                         DiagnosticSink& sink) {
  if ((source.disposition == flow::TimeAdvanceDisposition::committed) !=
      (source.reason == flow::StepFailureReason::none))
    throw DiagnosticCollectionError(DiagnosticFailureClass::invalid_input,
                                    "flow-driver.invalid", -1,
                                    "flow driver source is inconsistent");
  collect_common(describe_diagnostics(source), request, sink,
                 diagnostic_fingerprint_field_ids(source),
                 [&](FingerprintBuilder& fp) {
                   fp.add_u64("density_model",
                              static_cast<std::uint64_t>(
                                  source.density_model));
                   fp.add_u64("step", source.step);
                   fp.add_fp64("time", source.time_s);
                   fp.add_u64("attempt_count", source.attempt_count);
                   fp.add_u64("disposition",
                              static_cast<std::uint64_t>(
                                  source.disposition));
                   fp.add_u64("reason",
                              static_cast<std::uint64_t>(source.reason));
                   fp.add_i64("lowest_failing_rank",
                              source.lowest_failing_rank);
                 },
                 [&](DiagnosticRecord& record) {
                   set_failure(record, source.reason,
                               source.lowest_failing_rank);
                   if (request.level == DiagnosticLevel::summary) {
                     record.metrics = {
                         metric("density_model", "count",
                                as_double(static_cast<int>(
                                    source.density_model))),
                         metric("time", "s", source.time_s),
                         metric("disposition", "count",
                                as_double(static_cast<int>(
                                    source.disposition))),
                         metric("reason", "count",
                                as_double(static_cast<int>(source.reason))),
                         metric("lowest_failing_rank", "count",
                                as_double(source.lowest_failing_rank))};
                   } else {
                     record.counters = {
                         counter("attempt_count",
                                 checked_u64(source.attempt_count)),
                         counter("step", source.step)};
                   }
                 });
}

#undef HUNDUN_SIMPLE_SUMMARY_ADAPTER

}  // namespace hundun::diagnostics
