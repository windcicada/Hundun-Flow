// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/material_density_transport_diagnostics.hpp"

#include "hundun/flow/material_density_transport.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#include "material_density_transport_diagnostics_test_access.hpp"
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace hundun::diagnostics {
namespace {

constexpr std::string_view kModuleId = "hundun.flow.material_density_transport";
constexpr std::string_view kInstanceId = "primary";
constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(
        DiagnosticCapability::bounded_state_sample) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
test::MaterialDiagnosticWorkCounts diagnostic_work_counts;
struct MaterialDiagnosticTestControl final {
  std::optional<std::pair<std::size_t, std::uint64_t>> global_id_override;
  bool provider_failure{};
  bool record_failure{};
  bool request_size_overflow{};
  bool sample_wire_overflow{};
  int rank{-1};
};
MaterialDiagnosticTestControl diagnostic_test_control;
MaterialDiagnosticTestControl take_test_control(int rank) noexcept {
  auto result = std::exchange(diagnostic_test_control, {});
  if (result.rank >= 0 && result.rank != rank)
    return {};
  return result;
}
#define HUNDUN_DIAGNOSTIC_WORK(member) (++diagnostic_work_counts.member)
#else
#define HUNDUN_DIAGNOSTIC_WORK(member) ((void)0)
struct MaterialDiagnosticTestControl final {};
MaterialDiagnosticTestControl take_test_control(int) noexcept { return {}; }
#endif

[[noreturn]] void collection_error(DiagnosticFailureClass classification,
                                   std::string code, int rank,
                                   std::string message) {
  throw DiagnosticCollectionError(classification, std::move(code), rank,
                                  std::move(message));
}

DiagnosticFailureClass
failure_class(flow::MaterialTransportFailureReason reason) noexcept {
  using R = flow::MaterialTransportFailureReason;
  switch (reason) {
  case R::none:
    return DiagnosticFailureClass::none;
  case R::invalid_input:
    return DiagnosticFailureClass::invalid_input;
  case R::non_finite_state:
    return DiagnosticFailureClass::non_finite_state;
  case R::non_positive_density:
    return DiagnosticFailureClass::non_positive_state;
  case R::final_density_residual:
  case R::final_transport_residual:
    return DiagnosticFailureClass::non_convergence;
  case R::final_conservation_defect:
    return DiagnosticFailureClass::conservation;
  case R::collective_operation:
    return DiagnosticFailureClass::collective_operation;
  }
  return DiagnosticFailureClass::invalid_input;
}

std::string failure_code(flow::MaterialTransportFailureReason reason) {
  using R = flow::MaterialTransportFailureReason;
  switch (reason) {
  case R::none:
    return "none";
  case R::invalid_input:
    return "material.invalid-input";
  case R::non_finite_state:
    return "material.non-finite-state";
  case R::non_positive_density:
    return "material.non-positive-density";
  case R::final_density_residual:
    return "material.density-residual";
  case R::final_transport_residual:
    return "material.transport-residual";
  case R::final_conservation_defect:
    return "material.conservation";
  case R::collective_operation:
    return "material.collective-operation";
  }
  return "material.invalid-input";
}

struct LocalObservation final {
  double density_min{std::numeric_limits<double>::infinity()};
  double density_max{-std::numeric_limits<double>::infinity()};
  double total_mass{};
  std::vector<double> transported_totals;
  DiagnosticFingerprintParts fingerprint;
  std::uint64_t eligible_count{};
  std::vector<DiagnosticSample> samples;
};

bool selected(const DiagnosticRequest &request, std::string_view field) {
  return request.selected_fields.empty() ||
         std::binary_search(request.selected_fields.begin(),
                            request.selected_fields.end(), field);
}

void require_provider_fields(
    const flow::MaterialDensityDiagnosticSource &source,
    const DiagnosticRequest &request) {
  const auto ids = diagnostic_fingerprint_field_ids(source);
  for (const auto requested : request.selected_fields) {
    if (!std::binary_search(ids.begin(), ids.end(), requested))
      collection_error(DiagnosticFailureClass::invalid_request,
                       "material.diagnostics.unknown-field", -1,
                       "unknown material diagnostic field");
  }
}

LocalObservation
observe_local(const flow::MaterialDensityDiagnosticSource &source,
              const DiagnosticRequest &request,
              const MaterialDiagnosticTestControl &control) {
  static_cast<void>(control);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (control.provider_failure)
    collection_error(DiagnosticFailureClass::layout,
                     "material.diagnostics.injected-provider", -1,
                     "injected material diagnostic provider failure");
#endif
  const auto ids = diagnostic_fingerprint_field_ids(source);
  LocalObservation result;
  const bool summary = request.level == DiagnosticLevel::summary;
  const bool invariants = request.level == DiagnosticLevel::invariants;
  if (summary)
    result.transported_totals.assign(ids.size() - 1U, 0.0);
  DiagnosticFingerprintAccumulator fingerprint;
  for (std::size_t field = 0; field < ids.size(); ++field) {
    for (std::size_t cell = 0; cell < source.owned_cell_count(); ++cell) {
      auto global = source.global_cell_id(cell);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
      if (control.global_id_override &&
          control.global_id_override->first == cell)
        global = control.global_id_override->second;
#endif
      const double value = source.field_value(field, cell);
      HUNDUN_DIAGNOSTIC_WORK(field_reads);
      const auto represented = describe_fp64(value);
      fingerprint.add(ids[field], global, 0U, represented);
      HUNDUN_DIAGNOSTIC_WORK(fingerprint_items);
      if (field == 0U && (summary || invariants)) {
        HUNDUN_DIAGNOSTIC_WORK(summary_accumulations);
        if (std::isnan(value)) {
          if (!std::isnan(result.density_min)) {
            result.density_min = value;
            result.density_max = value;
          }
        } else if (!std::isnan(result.density_min)) {
          result.density_min = std::min(result.density_min, value);
          result.density_max = std::max(result.density_max, value);
        }
        if (summary) {
          const double volume = source.cell_volume_m3(cell);
          HUNDUN_DIAGNOSTIC_WORK(volume_reads);
          result.total_mass += volume * value;
        }
      } else if (summary) {
        const double volume = source.cell_volume_m3(cell);
        HUNDUN_DIAGNOSTIC_WORK(volume_reads);
        HUNDUN_DIAGNOSTIC_WORK(summary_accumulations);
        result.transported_totals[field - 1U] += volume * value;
      }
      if (request.level == DiagnosticLevel::bounded_state_sample &&
          selected(request, ids[field])) {
        if (result.eligible_count == std::numeric_limits<std::uint64_t>::max())
          collection_error(DiagnosticFailureClass::layout,
                           "material.diagnostics.sample-count-overflow", -1,
                           "material sample count overflows");
        ++result.eligible_count;
        HUNDUN_DIAGNOSTIC_WORK(sample_candidates);
        if (result.samples.size() < request.sample_budget) {
          result.samples.push_back({std::string(ids[field]), global, 0U,
                                    std::string(source.field_unit(field)),
                                    represented});
        }
      }
    }
  }
  result.fingerprint = fingerprint.parts();
  return result;
}

DiagnosticStateFingerprint finish_fingerprint(DiagnosticFingerprintParts p) {
  DiagnosticFingerprintAccumulator accumulator;
  accumulator.combine(p);
  return accumulator.finish();
}

std::string indexed_id(std::string_view prefix, std::size_t index) {
  std::string result(prefix);
  std::array<char, 20> digits{};
  for (std::size_t place = digits.size(); place != 0U; --place) {
    digits[place - 1U] = static_cast<char>('0' + index % 10U);
    index /= 10U;
  }
  result.append(digits.data(), digits.size());
  return result;
}

DiagnosticRecord
build_record(const flow::MaterialDensityDiagnosticSource &source,
             const DiagnosticRequest &request,
             const LocalObservation &observation, DiagnosticScope scope,
             int failure_rank) {
  const auto descriptor = describe_diagnostics(source);
  const auto &report = source.report();
  DiagnosticRecord record;
  record.schema_version = descriptor.schema_version;
  record.module_kind = descriptor.module_kind;
  record.module_id = std::string(descriptor.module_id);
  record.instance_id = std::string(descriptor.instance_id);
  record.level = request.level;
  record.scope = scope;
  record.rank = request.frame.rank;
  record.step = request.frame.step;
  record.time_s = describe_fp64(request.frame.time_s);
  record.phase = std::string(request.frame.phase);
  const bool ok =
      report.disposition() == flow::MaterialTransportDisposition::finalized;
  record.status = ok ? DiagnosticStatus::ok : DiagnosticStatus::failed;
  record.failure =
      ok ? DiagnosticFailure{}
         : DiagnosticFailure{failure_class(report.reason()),
                             failure_code(report.reason()), failure_rank};
  record.identities.push_back({"flow_state.attempt", std::nullopt,
                               report.attempt_identity(), std::nullopt,
                               std::nullopt});
  record.identities.push_back({"material.finalization", std::nullopt,
                               report.finalization_identity(), std::nullopt,
                               std::nullopt});
  record.state_fingerprint = finish_fingerprint(observation.fingerprint);

  if (request.level == DiagnosticLevel::summary) {
    record.metrics.push_back({"density.maximum",
                              DiagnosticMetricKind::state_summary, "kg/m3",
                              describe_fp64(observation.density_max)});
    record.metrics.push_back({"density.minimum",
                              DiagnosticMetricKind::state_summary, "kg/m3",
                              describe_fp64(observation.density_min)});
    record.metrics.push_back(
        {"density.residual", DiagnosticMetricKind::residual, "1",
         report.density_residual_available()
             ? describe_fp64(report.density_normalized_l2())
             : DiagnosticFp64{}});
    record.metrics.push_back(
        {"mass.relative-defect", DiagnosticMetricKind::conservation, "1",
         report.mass_conservation_available()
             ? describe_fp64(report.mass_relative_conservation_defect())
             : DiagnosticFp64{}});
    record.metrics.push_back({"mass.total", DiagnosticMetricKind::conservation,
                              "kg", describe_fp64(observation.total_mass)});
    for (std::size_t field = 0; field < report.transport_normalized_l2().size();
         ++field) {
      const std::string root = indexed_id("transport.s", field);
      record.metrics.push_back(
          {root + ".relative-defect", DiagnosticMetricKind::conservation, "1",
           report.transport_conservation_availability()[field] != 0U
               ? describe_fp64(
                     report.transport_relative_conservation_defect()[field])
               : DiagnosticFp64{}});
      record.metrics.push_back(
          {root + ".residual", DiagnosticMetricKind::residual, "1",
           report.transport_residual_availability()[field] != 0U
               ? describe_fp64(report.transport_normalized_l2()[field])
               : DiagnosticFp64{}});
      record.metrics.push_back(
          {root + ".total", DiagnosticMetricKind::conservation,
           field == 0U ? "J" : "kg",
           describe_fp64(observation.transported_totals[field])});
    }
  } else if (request.level == DiagnosticLevel::invariants) {
    auto add = [&](std::string id, std::string unit, double observed,
                   double limit, InvariantRelation relation) {
      DiagnosticInvariant invariant{std::move(id),
                                    std::move(unit),
                                    describe_fp64(observed),
                                    relation == InvariantRelation::positive
                                        ? DiagnosticFp64{}
                                        : describe_fp64(limit),
                                    relation,
                                    false};
      invariant.passed = evaluate_invariant(invariant);
      record.invariants.push_back(std::move(invariant));
    };
    add("density.positive", "kg/m3", observation.density_min, 0.0,
        InvariantRelation::positive);
    auto add_availability = [&](std::string id, bool available) {
      add(std::move(id), "1", available ? 1.0 : 0.0, 1.0,
          InvariantRelation::equal);
    };
    if (report.density_residual_available())
      add("density.residual", "1", report.density_normalized_l2(), 1.0e-10,
          InvariantRelation::less_equal);
    add_availability("density.residual.available",
                     report.density_residual_available());
    if (report.mass_conservation_available())
      add("mass.relative-defect", "1",
          report.mass_relative_conservation_defect(), 5.0e-11,
          InvariantRelation::less_equal);
    add_availability("mass.relative-defect.available",
                     report.mass_conservation_available());
    for (std::size_t field = 0; field < report.transport_normalized_l2().size();
         ++field) {
      const std::string root = indexed_id("transport.s", field);
      const bool conservation_available =
          report.transport_conservation_availability()[field] != 0U;
      if (conservation_available)
        add(root + ".relative-defect", "1",
            report.transport_relative_conservation_defect()[field], 5.0e-11,
            InvariantRelation::less_equal);
      add_availability(root + ".relative-defect.available",
                       conservation_available);
      const bool residual_available =
          report.transport_residual_availability()[field] != 0U;
      if (residual_available)
        add(root + ".residual", "1", report.transport_normalized_l2()[field],
            1.0e-9, InvariantRelation::less_equal);
      add_availability(root + ".residual.available", residual_available);
    }
  } else if (request.level == DiagnosticLevel::counters) {
    record.counters.push_back(
        {"attempt.identity", "count", report.attempt_identity()});
    record.counters.push_back(
        {"finalization.identity", "count", report.finalization_identity()});
    record.counters.push_back(
        {"flux.provenance", "count",
         static_cast<std::uint64_t>(report.flux_provenance())});
    record.counters.push_back(
        {"owned.cells", "count",
         static_cast<std::uint64_t>(source.owned_cell_count())});
    record.counters.push_back(
        {"transport.fields", "count",
         static_cast<std::uint64_t>(report.transport_normalized_l2().size())});
  } else {
    record.sample_budget = request.sample_budget;
    record.eligible_sample_count = observation.eligible_count;
    record.samples_truncated =
        observation.eligible_count > observation.samples.size();
    record.samples = observation.samples;
  }
  return record;
}

template <class Function>
void submit_local(Function &&build, DiagnosticSink &sink) {
  DiagnosticRecord record = build();
  try {
    sink.submit(record);
  } catch (const std::exception &) {
    collection_error(DiagnosticFailureClass::sink_failure,
                     "diagnostics.sink.submit", -1,
                     "diagnostic sink submission failed");
  } catch (...) {
    collection_error(DiagnosticFailureClass::sink_failure,
                     "diagnostics.sink.submit", -1,
                     "diagnostic sink submission failed");
  }
}

int first_failing_rank(const runtime::MpiContext &mpi, bool failed,
                       std::string_view operation) {
  int local = failed ? mpi.rank() : mpi.size();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      operation);
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  return lowest == mpi.size() ? -1 : lowest;
}

std::string broadcast_failure_code(const runtime::MpiContext &mpi, int root,
                                   std::string local_code) {
  std::uint64_t size =
      mpi.rank() == root ? static_cast<std::uint64_t>(local_code.size()) : 0U;
  runtime::check_mpi_result(MPI_Bcast(&size, 1, MPI_UINT64_T, root, mpi.comm()),
                            "MPI_Bcast(material diagnostic failure code size)");
  if (size > 256U)
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.failure-code-size", root,
                     "material diagnostic failure code is too large");
  std::array<char, 256> bytes{};
  if (mpi.rank() == root)
    std::copy(local_code.begin(), local_code.end(), bytes.begin());
  runtime::check_mpi_result(MPI_Bcast(bytes.data(), static_cast<int>(size),
                                      MPI_CHAR, root, mpi.comm()),
                            "MPI_Bcast(material diagnostic failure code)");
  return std::string(bytes.data(), static_cast<std::size_t>(size));
}

std::vector<unsigned char> request_key(const DiagnosticRequest &request) {
  std::vector<unsigned char> result;
  auto append = [&](const void *data, std::size_t size) {
    const auto *bytes = static_cast<const unsigned char *>(data);
    result.insert(result.end(), bytes, bytes + size);
  };
  const auto level = static_cast<std::uint8_t>(request.level);
  const auto scope = static_cast<std::uint8_t>(request.scope);
  append(&level, sizeof(level));
  append(&scope, sizeof(scope));
  append(&request.frame.step, sizeof(request.frame.step));
  const auto time = describe_fp64(request.frame.time_s).bits;
  append(&time, sizeof(time));
  const auto append_text = [&](std::string_view text) {
    const auto size = static_cast<std::uint64_t>(text.size());
    append(&size, sizeof(size));
    append(text.data(), text.size());
  };
  append_text(request.frame.phase);
  const auto fields =
      static_cast<std::uint64_t>(request.selected_fields.size());
  append(&fields, sizeof(fields));
  for (const auto field : request.selected_fields)
    append_text(field);
  const auto budget = static_cast<std::uint64_t>(request.sample_budget);
  append(&budget, sizeof(budget));
  return result;
}

void require_request_agreement(const runtime::MpiContext &mpi,
                               const DiagnosticRequest &request,
                               const MaterialDiagnosticTestControl &control) {
  constexpr std::size_t limit =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t serialized_size = 2U + sizeof(request.frame.step) +
                                sizeof(std::uint64_t) +
                                2U * sizeof(std::uint64_t);
  const auto add_text_size = [&](std::string_view value) {
    if (serialized_size > limit - sizeof(std::uint64_t) ||
        value.size() > limit - sizeof(std::uint64_t) - serialized_size)
      return false;
    serialized_size += sizeof(std::uint64_t) + value.size();
    return true;
  };
  bool local_size_valid = add_text_size(request.frame.phase);
  for (const auto field : request.selected_fields)
    local_size_valid = local_size_valid && add_text_size(field);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (control.request_size_overflow)
    local_size_valid = false;
#else
  static_cast<void>(control);
#endif
  const int size_failure_rank = first_failing_rank(
      mpi, !local_size_valid,
      "MPI_Allreduce(material diagnostic request size validity)");
  if (size_failure_rank >= 0)
    collection_error(DiagnosticFailureClass::invalid_request,
                     "material.diagnostics.request-size", size_failure_rank,
                     "diagnostic request is too large");
  std::vector<unsigned char> key;
  bool key_failed = false;
  try {
    key = request_key(request);
  } catch (...) {
    key_failed = true;
  }
  const int key_failure_rank = first_failing_rank(
      mpi, key_failed,
      "MPI_Allreduce(material diagnostic request serialization)");
  if (key_failure_rank >= 0)
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.request-serialization",
                     key_failure_rank,
                     "diagnostic request serialization failed");
  std::uint64_t reference_size = static_cast<std::uint64_t>(key.size());
  runtime::check_mpi_result(
      MPI_Bcast(&reference_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(material diagnostic request size)");
  if (reference_size >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    collection_error(DiagnosticFailureClass::invalid_request,
                     "material.diagnostics.request-size", 0,
                     "diagnostic request is too large");
  std::vector<unsigned char> reference;
  bool reference_failed = false;
  try {
    reference.resize(static_cast<std::size_t>(reference_size));
    if (mpi.rank() == 0)
      reference = key;
  } catch (...) {
    reference_failed = true;
  }
  const int reference_failure_rank = first_failing_rank(
      mpi, reference_failed,
      "MPI_Allreduce(material diagnostic request buffer allocation)");
  if (reference_failure_rank >= 0)
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.request-buffer",
                     reference_failure_rank,
                     "diagnostic request buffer allocation failed");
  runtime::check_mpi_result(MPI_Bcast(reference.data(),
                                      static_cast<int>(reference.size()),
                                      MPI_BYTE, 0, mpi.comm()),
                            "MPI_Bcast(material diagnostic request)");
  const int lowest = first_failing_rank(
      mpi, key != reference, "MPI_Allreduce(material request agreement)");
  if (lowest >= 0)
    collection_error(DiagnosticFailureClass::invalid_request,
                     "material.diagnostics.request-agreement", lowest,
                     "collective diagnostic requests differ");
}

struct SampleWire final {
  std::array<char, 129> field{};
  std::uint64_t global{};
  std::uint32_t component{};
  std::uint64_t bits{};
  std::array<char, 8> unit{};
};

LocalObservation
collect_global_observation(const runtime::MpiContext &mpi,
                           const flow::MaterialDensityDiagnosticSource &source,
                           const DiagnosticRequest &request,
                           LocalObservation local,
                           const MaterialDiagnosticTestControl &control) {
  static_cast<void>(control);
  LocalObservation global;
  if (request.level == DiagnosticLevel::summary ||
      request.level == DiagnosticLevel::invariants) {
    const std::size_t scalar_count =
        request.level == DiagnosticLevel::summary ? 3U : 2U;
    std::array<double, 3> scalars{local.density_min, local.density_max,
                                  local.total_mass};
    std::vector<double> gathered(static_cast<std::size_t>(mpi.size()) *
                                 scalar_count);
    runtime::check_mpi_result(
        MPI_Allgather(scalars.data(), static_cast<int>(scalar_count),
                      MPI_DOUBLE, gathered.data(),
                      static_cast<int>(scalar_count), MPI_DOUBLE, mpi.comm()),
        "MPI_Allgather(material diagnostic summaries)");
    HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
    global.density_min = gathered[0];
    global.density_max = gathered[1];
    global.total_mass = 0.0;
    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto offset = static_cast<std::size_t>(rank) * scalar_count;
      if (std::isnan(global.density_min) || std::isnan(gathered[offset])) {
        if (!std::isnan(global.density_min)) {
          global.density_min = gathered[offset];
          global.density_max = gathered[offset + 1U];
        }
      } else {
        global.density_min = std::min(global.density_min, gathered[offset]);
        global.density_max =
            std::max(global.density_max, gathered[offset + 1U]);
      }
      if (request.level == DiagnosticLevel::summary)
        global.total_mass += gathered[offset + 2U];
    }
    if (request.level == DiagnosticLevel::summary) {
      global.transported_totals.resize(local.transported_totals.size());
      if (!global.transported_totals.empty()) {
        runtime::check_mpi_result(
            MPI_Allreduce(local.transported_totals.data(),
                          global.transported_totals.data(),
                          static_cast<int>(global.transported_totals.size()),
                          MPI_DOUBLE, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(material diagnostic transport totals)");
        HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
      }
    }
  }
  std::array<std::uint64_t, 2> local_parts{local.fingerprint.xor64,
                                           local.fingerprint.sum64};
  std::uint64_t global_xor{};
  std::uint64_t global_sum{};
  runtime::check_mpi_result(
      MPI_Allreduce(&local_parts[0], &global_xor, 1, MPI_UINT64_T, MPI_BXOR,
                    mpi.comm()),
      "MPI_Allreduce(material diagnostic fingerprint xor)");
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  runtime::check_mpi_result(
      MPI_Allreduce(&local_parts[1], &global_sum, 1, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(material diagnostic fingerprint sum)");
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  global.fingerprint = {global_xor, global_sum};

  if (request.level != DiagnosticLevel::bounded_state_sample)
    return global;

  std::uint64_t local_max_id{};
  for (std::size_t cell = 0; cell < source.owned_cell_count(); ++cell)
    local_max_id = std::max(local_max_id, source.global_cell_id(cell));
  std::uint64_t global_max_id{};
  runtime::check_mpi_result(
      MPI_Allreduce(&local_max_id, &global_max_id, 1, MPI_UINT64_T, MPI_MAX,
                    mpi.comm()),
      "MPI_Allreduce(material diagnostic maximum global ID)");
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  constexpr std::uint64_t chunk_width = 256U;
  for (std::uint64_t begin = 0U;;) {
    const bool final_chunk = global_max_id - begin < chunk_width;
    const std::uint64_t end = final_chunk ? global_max_id : begin + chunk_width;
    std::vector<std::uint64_t> local_ids;
    local_ids.reserve(static_cast<std::size_t>(end - begin) +
                      (final_chunk ? 1U : 0U));
    for (std::size_t cell = 0; cell < source.owned_cell_count(); ++cell) {
      auto id = source.global_cell_id(cell);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
      if (control.global_id_override &&
          control.global_id_override->first == cell)
        id = control.global_id_override->second;
#endif
      if (id >= begin && (final_chunk ? id <= end : id < end))
        local_ids.push_back(id);
    }
    std::sort(local_ids.begin(), local_ids.end());
    const bool local_duplicate =
        std::adjacent_find(local_ids.begin(), local_ids.end()) !=
        local_ids.end();
    const int duplicate_rank = first_failing_rank(
        mpi, local_duplicate,
        "MPI_Allreduce(material diagnostic local ownership uniqueness)");
    if (duplicate_rank >= 0)
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.duplicate-owned-sample",
                       duplicate_rank,
                       "duplicate owned material diagnostic tuple");
    const int local_count = static_cast<int>(local_ids.size());
    std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
    runtime::check_mpi_result(
        MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT,
                      mpi.comm()),
        "MPI_Allgather(material diagnostic ownership counts)");
    HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
    std::vector<int> displacements(counts.size());
    int total{};
    for (std::size_t rank = 0; rank < counts.size(); ++rank) {
      displacements[rank] = total;
      total += counts[rank];
    }
    std::vector<std::uint64_t> all_ids(static_cast<std::size_t>(total));
    runtime::check_mpi_result(
        MPI_Allgatherv(local_ids.data(), local_count, MPI_UINT64_T,
                       all_ids.data(), counts.data(), displacements.data(),
                       MPI_UINT64_T, mpi.comm()),
        "MPI_Allgatherv(material diagnostic ownership IDs)");
    HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
    std::sort(all_ids.begin(), all_ids.end());
    if (std::adjacent_find(all_ids.begin(), all_ids.end()) != all_ids.end())
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.duplicate-owned-sample", 0,
                       "duplicate owned material diagnostic tuple");
    if (final_chunk)
      break;
    begin = end;
  }

  std::vector<std::uint64_t> eligible(static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(
      MPI_Allgather(&local.eligible_count, 1, MPI_UINT64_T, eligible.data(), 1,
                    MPI_UINT64_T, mpi.comm()),
      "MPI_Allgather(material diagnostic eligible counts)");
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  for (const auto count : eligible) {
    if (count >
        std::numeric_limits<std::uint64_t>::max() - global.eligible_count)
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.global-count-overflow", 0,
                       "global material sample count overflows");
    global.eligible_count += count;
  }
  std::vector<SampleWire> wire(local.samples.size());
  for (std::size_t i = 0; i < local.samples.size(); ++i) {
    std::copy(local.samples[i].field_id.begin(),
              local.samples[i].field_id.end(), wire[i].field.begin());
    std::copy(local.samples[i].unit.begin(), local.samples[i].unit.end(),
              wire[i].unit.begin());
    wire[i].global = local.samples[i].global_id;
    wire[i].component = local.samples[i].component;
    wire[i].bits = local.samples[i].value.bits;
  }
  std::uint64_t local_bytes =
      static_cast<std::uint64_t>(wire.size() * sizeof(SampleWire));
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (control.sample_wire_overflow)
    local_bytes =
        static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1U;
#endif
  std::vector<std::uint64_t> byte_counts(static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(MPI_Allgather(&local_bytes, 1, MPI_UINT64_T,
                                          byte_counts.data(), 1, MPI_UINT64_T,
                                          mpi.comm()),
                            "MPI_Allgather(material diagnostic sample sizes)");
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  std::vector<int> counts(byte_counts.size());
  std::vector<int> displacements(byte_counts.size());
  std::size_t total_bytes = 0U;
  for (std::size_t rank = 0; rank < byte_counts.size(); ++rank) {
    if (byte_counts[rank] >
            static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        total_bytes >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                static_cast<std::size_t>(byte_counts[rank]))
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.sample-wire-size", 0,
                       "material diagnostic sample exchange is too large");
    counts[rank] = static_cast<int>(byte_counts[rank]);
    displacements[rank] = static_cast<int>(total_bytes);
    total_bytes += static_cast<std::size_t>(byte_counts[rank]);
  }
  std::vector<unsigned char> all_bytes(total_bytes);
  runtime::check_mpi_result(
      MPI_Allgatherv(wire.data(), static_cast<int>(local_bytes), MPI_BYTE,
                     all_bytes.data(), counts.data(), displacements.data(),
                     MPI_BYTE, mpi.comm()),
      "MPI_Allgatherv(material diagnostic samples)");
  HUNDUN_DIAGNOSTIC_WORK(raw_collectives);
  const auto total_items = total_bytes / sizeof(SampleWire);
  std::vector<SampleWire> all(total_items);
  if (total_bytes != 0U)
    std::memcpy(all.data(), all_bytes.data(), total_bytes);
  std::sort(all.begin(), all.end(), [](const auto &left, const auto &right) {
    return std::tie(left.field, left.global, left.component) <
           std::tie(right.field, right.global, right.component);
  });
  for (std::size_t i = 1; i < all.size(); ++i) {
    if (std::tie(all[i - 1U].field, all[i - 1U].global,
                 all[i - 1U].component) ==
        std::tie(all[i].field, all[i].global, all[i].component))
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.duplicate-owned-sample", 0,
                       "duplicate owned material diagnostic tuple");
  }
  const std::size_t retain = std::min(request.sample_budget, all.size());
  global.samples.reserve(retain);
  for (std::size_t i = 0; i < retain; ++i) {
    const std::string field(all[i].field.data());
    const std::string unit(all[i].unit.data());
    double value{};
    std::memcpy(&value, &all[i].bits, sizeof(value));
    global.samples.push_back(
        {field, all[i].global, all[i].component, unit, describe_fp64(value)});
  }
  return global;
}

} // namespace

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
void test::MaterialDensityTransportDiagnosticsTestAccess::reset() noexcept {
  diagnostic_work_counts = {};
  diagnostic_test_control = {};
}
test::MaterialDiagnosticWorkCounts
test::MaterialDensityTransportDiagnosticsTestAccess::work_counts() noexcept {
  return diagnostic_work_counts;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::override_global_id(
    std::size_t local_cell, std::uint64_t global_id, int rank) noexcept {
  diagnostic_test_control.global_id_override = {{local_cell, global_id}};
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_provider_failure(int rank) noexcept {
  diagnostic_test_control.provider_failure = true;
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::inject_record_failure(
    int rank) noexcept {
  diagnostic_test_control.record_failure = true;
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_request_size_overflow(int rank) noexcept {
  diagnostic_test_control.request_size_overflow = true;
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_sample_wire_overflow(int rank) noexcept {
  diagnostic_test_control.sample_wire_overflow = true;
  diagnostic_test_control.rank = rank;
}
#endif

DiagnosticDescriptor
describe_diagnostics(const flow::MaterialDensityDiagnosticSource &) noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::density_transport,
          kModuleId, kInstanceId, kCapabilities};
}

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::MaterialDensityDiagnosticSource &source) {
  std::vector<std::string_view> result;
  result.reserve(source.fingerprint_field_count());
  for (std::size_t field = 0; field < source.fingerprint_field_count(); ++field)
    result.push_back(source.fingerprint_field_id(field));
  if (!std::is_sorted(result.begin(), result.end()) ||
      std::adjacent_find(result.begin(), result.end()) != result.end())
    collection_error(DiagnosticFailureClass::layout,
                     "material.diagnostics.field-order", -1,
                     "material diagnostic field IDs are not canonical");
  return result;
}

void collect_diagnostics(const flow::MaterialDensityDiagnosticSource &source,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const auto control = take_test_control(request.frame.rank);
  DiagnosticDescriptor descriptor;
  LocalObservation observation;
  try {
    descriptor = describe_diagnostics(source);
    validate(descriptor);
    validate(request, descriptor);
    if (request.scope != DiagnosticScope::local)
      collection_error(DiagnosticFailureClass::invalid_request,
                       "material.diagnostics.local-scope", -1,
                       "local adapter requires local scope");
    if (request.frame.phase != "material.finalized-trial")
      collection_error(DiagnosticFailureClass::invalid_request,
                       "material.diagnostics.phase", -1,
                       "material diagnostic phase is invalid");
    require_provider_fields(source, request);
    observation = observe_local(source, request, control);
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (const std::exception &) {
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.provider", -1,
                     "material diagnostic provider failed");
  } catch (...) {
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.provider", -1,
                     "material diagnostic provider failed");
  }
  DiagnosticRecord record;
  try {
    record =
        build_record(source, request, observation, DiagnosticScope::local, -1);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (control.record_failure)
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.injected-record", -1,
                       "injected material diagnostic record failure");
#endif
    validate(record);
    validate(record, descriptor, request);
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (const std::exception &) {
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.record", -1,
                     "material diagnostic record failed validation");
  } catch (...) {
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.record", -1,
                     "material diagnostic record failed validation");
  }
  submit_local([&] { return std::move(record); }, sink);
}

void collect_diagnostics(const flow::MaterialDensityDiagnosticSource &source,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const auto control = take_test_control(mpi.rank());
  const auto descriptor = describe_diagnostics(source);
  bool preflight_failed = false;
  DiagnosticFailureClass classification = DiagnosticFailureClass::invalid_input;
  std::string code = "material.diagnostics.preflight";
  try {
    validate(descriptor);
    validate(request, descriptor);
    if (request.scope != DiagnosticScope::collective ||
        request.frame.rank != mpi.rank())
      collection_error(DiagnosticFailureClass::invalid_request,
                       "material.diagnostics.collective-frame", -1,
                       "collective diagnostic rank or scope is invalid");
    if (request.frame.phase != "material.finalized-trial")
      collection_error(DiagnosticFailureClass::invalid_request,
                       "material.diagnostics.phase", -1,
                       "material diagnostic phase is invalid");
    require_provider_fields(source, request);
  } catch (const DiagnosticCollectionError &error) {
    preflight_failed = true;
    classification = error.classification();
    code = std::string(error.code());
  } catch (const std::exception &) {
    preflight_failed = true;
  } catch (...) {
    preflight_failed = true;
  }
  const int preflight_rank = first_failing_rank(
      mpi, preflight_failed, "MPI_Allreduce(material diagnostic preflight)");
  if (preflight_rank >= 0) {
    int class_value =
        mpi.rank() == preflight_rank ? static_cast<int>(classification) : 0;
    runtime::check_mpi_result(
        MPI_Bcast(&class_value, 1, MPI_INT, preflight_rank, mpi.comm()),
        "MPI_Bcast(material diagnostic preflight class)");
    code = broadcast_failure_code(mpi, preflight_rank, std::move(code));
    collection_error(static_cast<DiagnosticFailureClass>(class_value),
                     std::move(code), preflight_rank,
                     "material diagnostic preflight failed");
  }
  require_request_agreement(mpi, request, control);

  LocalObservation local;
  bool provider_failed = false;
  DiagnosticFailureClass provider_class = DiagnosticFailureClass::invalid_input;
  std::string provider_code = "material.diagnostics.provider";
  try {
    local = observe_local(source, request, control);
  } catch (const DiagnosticCollectionError &error) {
    provider_failed = true;
    provider_class = error.classification();
    provider_code = std::string(error.code());
  } catch (const std::exception &) {
    provider_failed = true;
  } catch (...) {
    provider_failed = true;
  }
  const int provider_rank = first_failing_rank(
      mpi, provider_failed, "MPI_Allreduce(material diagnostic provider)");
  if (provider_rank >= 0) {
    int class_value =
        mpi.rank() == provider_rank ? static_cast<int>(provider_class) : 0;
    runtime::check_mpi_result(
        MPI_Bcast(&class_value, 1, MPI_INT, provider_rank, mpi.comm()),
        "MPI_Bcast(material diagnostic provider class)");
    provider_code =
        broadcast_failure_code(mpi, provider_rank, std::move(provider_code));
    collection_error(static_cast<DiagnosticFailureClass>(class_value),
                     std::move(provider_code), provider_rank,
                     "material diagnostic provider failed");
  }

  LocalObservation global = collect_global_observation(
      mpi, source, request, std::move(local), control);
  DiagnosticRecord record;
  bool validation_failed = false;
  DiagnosticFailureClass validation_class =
      DiagnosticFailureClass::invalid_input;
  std::string validation_code = "material.diagnostics.record";
  try {
    record = build_record(source, request, global, DiagnosticScope::collective,
                          source.report().lowest_failing_rank());
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (control.record_failure)
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.injected-record", -1,
                       "injected material diagnostic record failure");
#endif
    validate(record);
    validate(record, descriptor, request);
  } catch (const DiagnosticCollectionError &error) {
    validation_failed = true;
    validation_class = error.classification();
    validation_code = std::string(error.code());
  } catch (const std::exception &) {
    validation_failed = true;
  } catch (...) {
    validation_failed = true;
  }
  const int validation_rank = first_failing_rank(
      mpi, validation_failed, "MPI_Allreduce(material record validation)");
  if (validation_rank >= 0) {
    int class_value =
        mpi.rank() == validation_rank ? static_cast<int>(validation_class) : 0;
    runtime::check_mpi_result(
        MPI_Bcast(&class_value, 1, MPI_INT, validation_rank, mpi.comm()),
        "MPI_Bcast(material diagnostic record class)");
    validation_code = broadcast_failure_code(mpi, validation_rank,
                                             std::move(validation_code));
    collection_error(static_cast<DiagnosticFailureClass>(class_value),
                     std::move(validation_code), validation_rank,
                     "material diagnostic record failed validation");
  }

  bool sink_failed = false;
  try {
    sink.submit(record);
  } catch (...) {
    sink_failed = true;
  }
  const int sink_rank = first_failing_rank(
      mpi, sink_failed, "MPI_Allreduce(material diagnostic sink)");
  if (sink_rank >= 0)
    collection_error(DiagnosticFailureClass::sink_failure,
                     "diagnostics.sink.submit", sink_rank,
                     "material diagnostic sink submission failed");
}

} // namespace hundun::diagnostics

#undef HUNDUN_DIAGNOSTIC_WORK
