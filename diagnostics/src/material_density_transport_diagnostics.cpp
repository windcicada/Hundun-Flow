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
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
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

enum class AggregationPreparationPoint : std::uint8_t {
  summary_gather,
  transport_totals,
  ownership_box_proof,
  eligible_counts,
  local_sample_wire_and_size_counts,
  sample_exchange_buffers,
  decoded_and_retained_samples,
  provider_key,
  provider_reference
};

enum class RawCollectivePoint : std::uint8_t {
  none,
  other,
  preparation_summary_gather,
  preparation_transport_totals,
  preparation_ownership_box_proof,
  preparation_eligible_counts,
  preparation_local_sample_wire_and_size_counts,
  preparation_sample_exchange_buffers,
  preparation_decoded_and_retained_samples,
  sample_size_exchange,
  preparation_provider_key,
  preparation_provider_reference
};

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
test::MaterialDiagnosticWorkCounts diagnostic_work_counts;
struct MaterialDiagnosticTestControl final {
  std::optional<std::pair<std::size_t, std::uint64_t>> global_id_override;
  std::optional<runtime::Box3> owned_box_override;
  bool transport_total_count_overflow{};
  bool provider_failure{};
  bool record_failure{};
  bool request_size_overflow{};
  bool sample_wire_overflow{};
  bool sample_wire_prefix_overflow{};
  bool provider_key_value_difference{};
  std::optional<AggregationPreparationPoint> allocation_failure;
  std::optional<std::vector<unsigned char>> request_wire_override;
  std::optional<std::vector<unsigned char>> sample_wire_override;
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

void record_raw_collective(
    RawCollectivePoint point = RawCollectivePoint::other) noexcept {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  ++diagnostic_work_counts.raw_collectives;
  diagnostic_work_counts.last_raw_collective =
      static_cast<test::MaterialDiagnosticRawCollectivePoint>(point);
#else
  static_cast<void>(point);
#endif
}

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
  const std::string layout_fingerprint(
      scope == DiagnosticScope::local
          ? source.owned_cell_layout_fingerprint()
          : source.global_cell_layout_fingerprint());
  const auto field_ids = diagnostic_fingerprint_field_ids(source);
  record.identities.reserve(field_ids.size() + 3U);
  for (const auto field_id : field_ids) {
    record.identities.push_back({"field." + std::string(field_id),
                                 layout_fingerprint, std::nullopt, std::nullopt,
                                 std::nullopt});
  }
  record.identities.push_back({"flow_state.attempt", std::nullopt,
                               report.attempt_identity(), std::nullopt,
                               std::nullopt});
  record.identities.push_back({"layout.cells", layout_fingerprint, std::nullopt,
                               std::nullopt, std::nullopt});
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
                       std::string_view operation,
                       RawCollectivePoint point = RawCollectivePoint::other) {
  int local = failed ? mpi.rank() : mpi.size();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      operation);
  record_raw_collective(point);
  return lowest == mpi.size() ? -1 : lowest;
}

RawCollectivePoint
preparation_trace(AggregationPreparationPoint point) noexcept {
  switch (point) {
  case AggregationPreparationPoint::summary_gather:
    return RawCollectivePoint::preparation_summary_gather;
  case AggregationPreparationPoint::transport_totals:
    return RawCollectivePoint::preparation_transport_totals;
  case AggregationPreparationPoint::ownership_box_proof:
    return RawCollectivePoint::preparation_ownership_box_proof;
  case AggregationPreparationPoint::eligible_counts:
    return RawCollectivePoint::preparation_eligible_counts;
  case AggregationPreparationPoint::local_sample_wire_and_size_counts:
    return RawCollectivePoint::preparation_local_sample_wire_and_size_counts;
  case AggregationPreparationPoint::sample_exchange_buffers:
    return RawCollectivePoint::preparation_sample_exchange_buffers;
  case AggregationPreparationPoint::decoded_and_retained_samples:
    return RawCollectivePoint::preparation_decoded_and_retained_samples;
  case AggregationPreparationPoint::provider_key:
    return RawCollectivePoint::preparation_provider_key;
  case AggregationPreparationPoint::provider_reference:
    return RawCollectivePoint::preparation_provider_reference;
  }
  return RawCollectivePoint::other;
}

template <class Function>
void prepare_aggregation(const runtime::MpiContext &mpi,
                         AggregationPreparationPoint point,
                         const MaterialDiagnosticTestControl &control,
                         Function &&prepare) {
  bool failed = false;
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (control.allocation_failure && *control.allocation_failure == point)
      throw std::bad_alloc{};
#else
    static_cast<void>(control);
#endif
    prepare();
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    failed = true;
  }
  const int rank = first_failing_rank(
      mpi, failed, "MPI_Allreduce(material diagnostic aggregation preparation)",
      preparation_trace(point));
  if (rank >= 0)
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.aggregation-preparation", rank,
                     "material diagnostic aggregation preparation failed");
}

template <class Function>
void prepare_provider_agreement(const runtime::MpiContext &mpi,
                                AggregationPreparationPoint point,
                                const MaterialDiagnosticTestControl &control,
                                Function &&prepare) {
  bool failed = false;
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (control.allocation_failure && *control.allocation_failure == point)
      throw std::bad_alloc{};
#else
    static_cast<void>(control);
#endif
    prepare();
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    failed = true;
  }
  const int rank = first_failing_rank(
      mpi, failed, "MPI_Allreduce(material diagnostic provider preparation)",
      preparation_trace(point));
  if (rank >= 0)
    collection_error(
        DiagnosticFailureClass::invalid_input,
        "material.diagnostics.provider-agreement-preparation", rank,
        "material diagnostic provider agreement preparation failed");
}

std::string broadcast_failure_code(const runtime::MpiContext &mpi, int root,
                                   std::string local_code) {
  std::uint64_t size =
      mpi.rank() == root ? static_cast<std::uint64_t>(local_code.size()) : 0U;
  runtime::check_mpi_result(MPI_Bcast(&size, 1, MPI_UINT64_T, root, mpi.comm()),
                            "MPI_Bcast(material diagnostic failure code size)");
  record_raw_collective();
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
  record_raw_collective();
  return std::string(bytes.data(), static_cast<std::size_t>(size));
}

class WireError final : public std::runtime_error {
public:
  explicit WireError(const char *message) : std::runtime_error(message) {}
};

class ByteWriter final {
public:
  void append_u8(std::uint8_t value) { bytes_.push_back(value); }

  void append_u16(std::uint16_t value) {
    for (unsigned shift = 0U; shift < 16U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }

  void append_u32(std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }

  void append_u64(std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }

  void append_text(std::string_view value) {
    append_u64(static_cast<std::uint64_t>(value.size()));
    bytes_.insert(bytes_.end(), value.begin(), value.end());
  }

  std::vector<unsigned char> finish() && { return std::move(bytes_); }

private:
  std::vector<unsigned char> bytes_;
};

class ByteReader final {
public:
  explicit ByteReader(const std::vector<unsigned char> &bytes) noexcept
      : bytes_(&bytes) {}

  std::uint8_t read_u8() { return read_integer<std::uint8_t>(); }
  std::uint16_t read_u16() { return read_integer<std::uint16_t>(); }
  std::uint32_t read_u32() { return read_integer<std::uint32_t>(); }
  std::uint64_t read_u64() { return read_integer<std::uint64_t>(); }

  std::string read_text(std::size_t maximum) {
    const auto encoded_size = read_u64();
    if (encoded_size > static_cast<std::uint64_t>(maximum) ||
        encoded_size > static_cast<std::uint64_t>(remaining()))
      throw WireError("diagnostic wire text length is invalid");
    const auto size = static_cast<std::size_t>(encoded_size);
    const auto first = bytes_->begin() + static_cast<std::ptrdiff_t>(offset_);
    offset_ += size;
    return std::string(first, first + static_cast<std::ptrdiff_t>(size));
  }

  std::size_t remaining() const noexcept { return bytes_->size() - offset_; }
  void require_end() const {
    if (remaining() != 0U)
      throw WireError("diagnostic wire has trailing bytes");
  }

private:
  template <class Integer> Integer read_integer() {
    if (sizeof(Integer) > remaining())
      throw WireError("diagnostic wire integer is truncated");
    std::uint64_t value{};
    for (std::size_t byte = 0; byte < sizeof(Integer); ++byte)
      value |= static_cast<std::uint64_t>((*bytes_)[offset_ + byte])
               << (8U * byte);
    offset_ += sizeof(Integer);
    return static_cast<Integer>(value);
  }

  const std::vector<unsigned char> *bytes_{};
  std::size_t offset_{};
};

constexpr std::uint32_t kRequestWireSchemaV1 = 1U;
constexpr std::uint32_t kSampleWireSchemaV1 = 1U;

bool valid_wire_id(std::string_view value) noexcept {
  if (value.empty() || value.size() > 128U)
    return false;
  return std::all_of(value.begin(), value.end(), [](char character) {
    return (character >= 'a' && character <= 'z') ||
           (character >= '0' && character <= '9') || character == '_' ||
           character == '.' || character == '-';
  });
}

bool valid_wire_unit(std::string_view unit) noexcept {
  constexpr std::array<std::string_view, 19> units{
      "1",    "s",     "m",      "m2",     "m3",   "m/s", "m2/s",
      "kg",   "kg/m3", "kg/s",   "kg*m/s", "Pa",   "J",   "J/m3",
      "J/kg", "K",     "degree", "byte",   "count"};
  return std::find(units.begin(), units.end(), unit) != units.end();
}

DiagnosticValueStatus value_status_for_bits(std::uint64_t bits) noexcept {
  const std::uint64_t exponent = (bits >> 52U) & UINT64_C(0x7ff);
  const std::uint64_t fraction = bits & UINT64_C(0x000fffffffffffff);
  if (exponent != UINT64_C(0x7ff))
    return DiagnosticValueStatus::finite;
  if (fraction == 0U)
    return (bits >> 63U) == 0U ? DiagnosticValueStatus::positive_infinity
                               : DiagnosticValueStatus::negative_infinity;
  return (fraction & UINT64_C(0x0008000000000000)) != 0U
             ? DiagnosticValueStatus::quiet_nan
             : DiagnosticValueStatus::signaling_nan;
}

void append_fp64(ByteWriter &writer, DiagnosticFp64 value) {
  writer.append_u8(static_cast<std::uint8_t>(value.status));
  writer.append_u64(value.bits);
}

DiagnosticFp64 read_fp64(ByteReader &reader, bool allow_non_finite) {
  const auto status_value = reader.read_u8();
  if (status_value >
      static_cast<std::uint8_t>(DiagnosticValueStatus::signaling_nan))
    throw WireError("diagnostic wire FP64 status is invalid");
  const auto status = static_cast<DiagnosticValueStatus>(status_value);
  const auto bits = reader.read_u64();
  if (value_status_for_bits(bits) != status ||
      (!allow_non_finite && status != DiagnosticValueStatus::finite))
    throw WireError("diagnostic wire FP64 representation is inconsistent");
  return {status, bits};
}

struct DecodedRequestWire final {
  DiagnosticLevel level{};
  DiagnosticScope scope{};
  std::uint64_t step{};
  DiagnosticFp64 time;
  std::string phase;
  std::vector<std::string> selected_fields;
  std::uint64_t sample_budget{};
};

std::vector<unsigned char>
encode_request_wire(const DecodedRequestWire &request) {
  ByteWriter writer;
  writer.append_u32(kRequestWireSchemaV1);
  writer.append_u8(static_cast<std::uint8_t>(request.level));
  writer.append_u8(static_cast<std::uint8_t>(request.scope));
  writer.append_u64(request.step);
  append_fp64(writer, request.time);
  writer.append_text(request.phase);
  writer.append_u64(static_cast<std::uint64_t>(request.selected_fields.size()));
  for (const auto &field : request.selected_fields)
    writer.append_text(field);
  writer.append_u64(request.sample_budget);
  return std::move(writer).finish();
}

std::vector<unsigned char> request_key(const DiagnosticRequest &request) {
  DecodedRequestWire wire;
  wire.level = request.level;
  wire.scope = request.scope;
  wire.step = request.frame.step;
  wire.time = describe_fp64(request.frame.time_s);
  wire.phase = request.frame.phase;
  wire.selected_fields.reserve(request.selected_fields.size());
  for (const auto field : request.selected_fields)
    wire.selected_fields.emplace_back(field);
  wire.sample_budget = static_cast<std::uint64_t>(request.sample_budget);
  return encode_request_wire(wire);
}

DecodedRequestWire
decode_request_wire(const std::vector<unsigned char> &bytes) {
  ByteReader reader(bytes);
  if (reader.read_u32() != kRequestWireSchemaV1)
    throw WireError("diagnostic request wire schema is invalid");
  const auto level = reader.read_u8();
  const auto scope = reader.read_u8();
  if (level >
          static_cast<std::uint8_t>(DiagnosticLevel::bounded_state_sample) ||
      scope > static_cast<std::uint8_t>(DiagnosticScope::collective))
    throw WireError("diagnostic request wire enum is invalid");
  DecodedRequestWire result;
  result.level = static_cast<DiagnosticLevel>(level);
  result.scope = static_cast<DiagnosticScope>(scope);
  result.step = reader.read_u64();
  result.time = read_fp64(reader, false);
  if ((result.time.bits >> 63U) != 0U &&
      result.time.bits != UINT64_C(0x8000000000000000))
    throw WireError("diagnostic request wire time is negative");
  result.phase = reader.read_text(128U);
  if (!valid_wire_id(result.phase))
    throw WireError("diagnostic request wire phase is invalid");
  const auto field_count = reader.read_u64();
  if (field_count > static_cast<std::uint64_t>(reader.remaining() / 9U))
    throw WireError("diagnostic request wire field count is invalid");
  result.selected_fields.reserve(static_cast<std::size_t>(field_count));
  for (std::uint64_t field = 0; field < field_count; ++field) {
    auto value = reader.read_text(128U);
    if (!valid_wire_id(value))
      throw WireError("diagnostic request wire field is invalid");
    if (!result.selected_fields.empty() &&
        !(result.selected_fields.back() < value))
      throw WireError("diagnostic request wire fields are not canonical");
    result.selected_fields.push_back(std::move(value));
  }
  result.sample_budget = reader.read_u64();
  if (result.sample_budget > static_cast<std::uint64_t>(kMaximumStateSamplesV1))
    throw WireError("diagnostic request wire budget is invalid");
  if ((result.level == DiagnosticLevel::bounded_state_sample &&
       result.sample_budget == 0U) ||
      (result.level != DiagnosticLevel::bounded_state_sample &&
       (result.sample_budget != 0U || !result.selected_fields.empty())))
    throw WireError("diagnostic request wire selection is invalid");
  reader.require_end();
  return result;
}

std::vector<unsigned char>
encode_sample_wire(const std::vector<DiagnosticSample> &samples) {
  ByteWriter writer;
  writer.append_u32(kSampleWireSchemaV1);
  writer.append_u64(static_cast<std::uint64_t>(samples.size()));
  for (const auto &sample : samples) {
    writer.append_text(sample.field_id);
    writer.append_u64(sample.global_id);
    writer.append_u32(sample.component);
    writer.append_text(sample.unit);
    append_fp64(writer, sample.value);
  }
  return std::move(writer).finish();
}

std::vector<DiagnosticSample>
decode_sample_wire(const std::vector<unsigned char> &bytes) {
  ByteReader reader(bytes);
  if (reader.read_u32() != kSampleWireSchemaV1)
    throw WireError("diagnostic sample wire schema is invalid");
  const auto count = reader.read_u64();
  if (count > static_cast<std::uint64_t>(kMaximumStateSamplesV1) ||
      count > static_cast<std::uint64_t>(reader.remaining() / 30U))
    throw WireError("diagnostic sample wire count is invalid");
  std::vector<DiagnosticSample> result;
  result.reserve(static_cast<std::size_t>(count));
  for (std::uint64_t item = 0; item < count; ++item) {
    auto field = reader.read_text(128U);
    const auto global = reader.read_u64();
    const auto component = reader.read_u32();
    auto unit = reader.read_text(7U);
    const auto value = read_fp64(reader, true);
    if (!valid_wire_id(field) || !valid_wire_unit(unit))
      throw WireError("diagnostic sample wire text is invalid");
    result.push_back(
        {std::move(field), global, component, std::move(unit), value});
  }
  reader.require_end();
  return result;
}

std::vector<unsigned char>
provider_key(const flow::MaterialDensityDiagnosticSource &source,
             const MaterialDiagnosticTestControl &control) {
  constexpr std::uint32_t provider_key_schema = 1U;
  const auto descriptor = describe_diagnostics(source);
  const auto field_ids = diagnostic_fingerprint_field_ids(source);
  const auto &report = source.report();
  ByteWriter builder;
  builder.append_u32(provider_key_schema);
  builder.append_u32(descriptor.schema_version);
  builder.append_u16(static_cast<std::uint16_t>(descriptor.module_kind));
  builder.append_text(descriptor.module_id);
  builder.append_text(descriptor.instance_id);
  builder.append_u32(descriptor.capabilities);
  builder.append_u64(static_cast<std::uint64_t>(field_ids.size()));
  for (std::size_t field = 0; field < field_ids.size(); ++field) {
    builder.append_text(field_ids[field]);
    builder.append_text(source.field_unit(field));
  }
  builder.append_text(source.global_cell_layout_fingerprint());
  builder.append_u8(static_cast<std::uint8_t>(report.disposition()));
  builder.append_u8(static_cast<std::uint8_t>(report.reason()));
  builder.append_u32(static_cast<std::uint32_t>(report.lowest_failing_rank()));
  builder.append_u64(report.attempt_identity());
  builder.append_u64(report.finalization_identity());
  builder.append_u8(static_cast<std::uint8_t>(report.flux_provenance()));
  builder.append_u8(report.density_residual_available() ? 1U : 0U);
  builder.append_u64(describe_fp64(report.density_normalized_l2()).bits);
  builder.append_u8(report.mass_conservation_available() ? 1U : 0U);
  builder.append_u64(
      describe_fp64(report.mass_relative_conservation_defect()).bits);
  builder.append_u64(static_cast<std::uint64_t>(
      report.transport_residual_availability().size()));
  for (const auto available : report.transport_residual_availability())
    builder.append_u8(available);
  builder.append_u64(
      static_cast<std::uint64_t>(report.transport_normalized_l2().size()));
  for (std::size_t field = 0; field < report.transport_normalized_l2().size();
       ++field) {
    auto bits = describe_fp64(report.transport_normalized_l2()[field]).bits;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (field == 0U && control.provider_key_value_difference)
      bits ^= UINT64_C(1);
#else
    static_cast<void>(control);
#endif
    builder.append_u64(bits);
  }
  builder.append_u64(static_cast<std::uint64_t>(
      report.transport_conservation_availability().size()));
  for (const auto available : report.transport_conservation_availability())
    builder.append_u8(available);
  builder.append_u64(static_cast<std::uint64_t>(
      report.transport_relative_conservation_defect().size()));
  for (const auto value : report.transport_relative_conservation_defect())
    builder.append_u64(describe_fp64(value).bits);
  return std::move(builder).finish();
}

void require_provider_agreement(
    const runtime::MpiContext &mpi,
    const flow::MaterialDensityDiagnosticSource &source,
    const MaterialDiagnosticTestControl &control) {
  std::vector<unsigned char> key;
  prepare_provider_agreement(mpi, AggregationPreparationPoint::provider_key,
                             control,
                             [&] { key = provider_key(source, control); });

  std::uint64_t reference_size =
      mpi.rank() == 0 ? static_cast<std::uint64_t>(key.size()) : 0U;
  runtime::check_mpi_result(
      MPI_Bcast(&reference_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(material diagnostic provider key size)");
  record_raw_collective();
  if (reference_size >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    collection_error(DiagnosticFailureClass::invalid_input,
                     "material.diagnostics.provider-agreement-preparation", 0,
                     "material diagnostic provider agreement key is too large");

  std::vector<unsigned char> reference;
  prepare_provider_agreement(
      mpi, AggregationPreparationPoint::provider_reference, control, [&] {
        reference.resize(static_cast<std::size_t>(reference_size));
        if (mpi.rank() == 0)
          std::copy(key.begin(), key.end(), reference.begin());
      });
  runtime::check_mpi_result(MPI_Bcast(reference.data(),
                                      static_cast<int>(reference.size()),
                                      MPI_BYTE, 0, mpi.comm()),
                            "MPI_Bcast(material diagnostic provider key)");
  record_raw_collective();
  const int lowest = first_failing_rank(
      mpi, key != reference,
      "MPI_Allreduce(material diagnostic provider agreement)");
  if (lowest >= 0)
    collection_error(DiagnosticFailureClass::collective_operation,
                     "material.diagnostics.provider-agreement", lowest,
                     "collective diagnostic providers differ");
}

void require_request_agreement(const runtime::MpiContext &mpi,
                               const DiagnosticRequest &request,
                               const MaterialDiagnosticTestControl &control) {
  constexpr std::size_t limit =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t serialized_size =
      sizeof(std::uint32_t) + 2U * sizeof(std::uint8_t) +
      sizeof(std::uint64_t) + sizeof(std::uint8_t) + 3U * sizeof(std::uint64_t);
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
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (control.request_wire_override)
      key = *control.request_wire_override;
#endif
    static_cast<void>(decode_request_wire(key));
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
  record_raw_collective();
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
  record_raw_collective();
  const int lowest = first_failing_rank(
      mpi, key != reference, "MPI_Allreduce(material request agreement)");
  if (lowest >= 0)
    collection_error(DiagnosticFailureClass::invalid_request,
                     "material.diagnostics.request-agreement", lowest,
                     "collective diagnostic requests differ");
}

bool checked_product(std::uint64_t left, std::uint64_t right,
                     std::uint64_t &result) noexcept {
  if (right != 0U && left > std::numeric_limits<std::uint64_t>::max() / right)
    return false;
  result = left * right;
  return true;
}

bool checked_box_volume(runtime::Box3 box, std::uint64_t &result) noexcept {
  if (box.begin.x < 0 || box.begin.y < 0 || box.begin.z < 0 ||
      box.end.x <= box.begin.x || box.end.y <= box.begin.y ||
      box.end.z <= box.begin.z)
    return false;
  std::uint64_t xy{};
  return checked_product(static_cast<std::uint64_t>(box.end.x - box.begin.x),
                         static_cast<std::uint64_t>(box.end.y - box.begin.y),
                         xy) &&
         checked_product(
             xy, static_cast<std::uint64_t>(box.end.z - box.begin.z), result);
}

bool boxes_overlap(runtime::Box3 left, runtime::Box3 right) noexcept {
  return left.begin.x < right.end.x && right.begin.x < left.end.x &&
         left.begin.y < right.end.y && right.begin.y < left.end.y &&
         left.begin.z < right.end.z && right.begin.z < left.end.z;
}

void require_complete_ownership(
    const runtime::MpiContext &mpi,
    const flow::MaterialDensityDiagnosticSource &source,
    const MaterialDiagnosticTestControl &control) {
  const auto extent = source.global_cell_extent();
  auto box = source.owned_global_box();
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (control.owned_box_override)
    box = *control.owned_box_override;
#endif
  std::uint64_t box_volume{};
  const bool extent_valid = extent.x > 0 && extent.y > 0 && extent.z > 0;
  const bool box_valid =
      extent_valid && checked_box_volume(box, box_volume) &&
      box.end.x <= extent.x && box.end.y <= extent.y && box.end.z <= extent.z &&
      box_volume == static_cast<std::uint64_t>(source.owned_cell_count());
  bool ids_valid = true;
  for (std::size_t cell = 0; cell < source.owned_cell_count(); ++cell) {
    auto observed = source.global_cell_id(cell);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (control.global_id_override && control.global_id_override->first == cell)
      observed = control.global_id_override->second;
#endif
    if (observed != source.global_cell_id(cell)) {
      ids_valid = false;
      break;
    }
  }
  const int local_failure_rank = first_failing_rank(
      mpi, !box_valid || !ids_valid,
      "MPI_Allreduce(material diagnostic local ownership proof)");
  if (local_failure_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "material.diagnostics.duplicate-owned-sample",
                     local_failure_rank,
                     "material diagnostic ownership proof failed");

  std::array<std::uint64_t, 7> local{static_cast<std::uint64_t>(box.begin.x),
                                     static_cast<std::uint64_t>(box.begin.y),
                                     static_cast<std::uint64_t>(box.begin.z),
                                     static_cast<std::uint64_t>(box.end.x),
                                     static_cast<std::uint64_t>(box.end.y),
                                     static_cast<std::uint64_t>(box.end.z),
                                     box_volume};
  std::vector<std::uint64_t> gathered;
  std::vector<runtime::Box3> boxes;
  prepare_aggregation(
      mpi, AggregationPreparationPoint::ownership_box_proof, control, [&] {
        gathered.resize(static_cast<std::size_t>(mpi.size()) * local.size());
        boxes.resize(static_cast<std::size_t>(mpi.size()));
      });
  runtime::check_mpi_result(
      MPI_Allgather(local.data(), static_cast<int>(local.size()), MPI_UINT64_T,
                    gathered.data(), static_cast<int>(local.size()),
                    MPI_UINT64_T, mpi.comm()),
      "MPI_Allgather(material diagnostic ownership boxes)");
  record_raw_collective();

  std::uint64_t covered{};
  int invalid_rank = -1;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const auto offset = static_cast<std::size_t>(rank) * local.size();
    bool representable = true;
    for (std::size_t coordinate = 0; coordinate < 6U; ++coordinate)
      representable = representable && gathered[offset + coordinate] <=
                                           static_cast<std::uint64_t>(
                                               std::numeric_limits<int>::max());
    if (!representable) {
      invalid_rank = rank;
      break;
    }
    auto &candidate = boxes[static_cast<std::size_t>(rank)];
    candidate = {{static_cast<int>(gathered[offset]),
                  static_cast<int>(gathered[offset + 1U]),
                  static_cast<int>(gathered[offset + 2U])},
                 {static_cast<int>(gathered[offset + 3U]),
                  static_cast<int>(gathered[offset + 4U]),
                  static_cast<int>(gathered[offset + 5U])}};
    std::uint64_t candidate_volume{};
    if (!checked_box_volume(candidate, candidate_volume) ||
        candidate.end.x > extent.x || candidate.end.y > extent.y ||
        candidate.end.z > extent.z ||
        candidate_volume != gathered[offset + 6U] ||
        gathered[offset + 6U] >
            std::numeric_limits<std::uint64_t>::max() - covered) {
      invalid_rank = rank;
      break;
    }
    for (int earlier = 0; earlier < rank; ++earlier) {
      if (boxes_overlap(boxes[static_cast<std::size_t>(earlier)], candidate)) {
        invalid_rank = rank;
        break;
      }
    }
    if (invalid_rank >= 0)
      break;
    covered += gathered[offset + 6U];
  }
  std::uint64_t xy{};
  std::uint64_t global_volume{};
  const bool global_volume_valid =
      extent_valid &&
      checked_product(static_cast<std::uint64_t>(extent.x),
                      static_cast<std::uint64_t>(extent.y), xy) &&
      checked_product(xy, static_cast<std::uint64_t>(extent.z), global_volume);
  if (invalid_rank < 0 && (!global_volume_valid || covered != global_volume))
    invalid_rank = 0;
  if (invalid_rank >= 0)
    collection_error(
        DiagnosticFailureClass::layout,
        "material.diagnostics.duplicate-owned-sample", invalid_rank,
        "material diagnostic ownership boxes do not cover the layout");
}

LocalObservation
collect_global_observation(const runtime::MpiContext &mpi,
                           const flow::MaterialDensityDiagnosticSource &source,
                           const DiagnosticRequest &request,
                           LocalObservation local,
                           const MaterialDiagnosticTestControl &control) {
  LocalObservation global;
  if (request.level == DiagnosticLevel::summary ||
      request.level == DiagnosticLevel::invariants) {
    const std::size_t scalar_count =
        request.level == DiagnosticLevel::summary ? 3U : 2U;
    std::array<double, 3> scalars{local.density_min, local.density_max,
                                  local.total_mass};
    std::vector<double> gathered;
    prepare_aggregation(
        mpi, AggregationPreparationPoint::summary_gather, control, [&] {
          gathered.resize(static_cast<std::size_t>(mpi.size()) * scalar_count);
        });
    runtime::check_mpi_result(
        MPI_Allgather(scalars.data(), static_cast<int>(scalar_count),
                      MPI_DOUBLE, gathered.data(),
                      static_cast<int>(scalar_count), MPI_DOUBLE, mpi.comm()),
        "MPI_Allgather(material diagnostic summaries)");
    record_raw_collective();
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
      prepare_aggregation(
          mpi, AggregationPreparationPoint::transport_totals, control, [&] {
            global.transported_totals.resize(local.transported_totals.size());
          });
      if (!global.transported_totals.empty()) {
        std::uint64_t reported_count =
            static_cast<std::uint64_t>(global.transported_totals.size());
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
        if (control.transport_total_count_overflow)
          reported_count =
              static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1U;
#endif
        const int count_failure_rank = first_failing_rank(
            mpi,
            reported_count >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
            "MPI_Allreduce(material diagnostic transport total count)");
        if (count_failure_rank >= 0)
          collection_error(
              DiagnosticFailureClass::invalid_input,
              "material.diagnostics.transport-total-count", count_failure_rank,
              "material diagnostic transport total count is too large");
        runtime::check_mpi_result(
            MPI_Allreduce(local.transported_totals.data(),
                          global.transported_totals.data(),
                          static_cast<int>(reported_count), MPI_DOUBLE, MPI_SUM,
                          mpi.comm()),
            "MPI_Allreduce(material diagnostic transport totals)");
        record_raw_collective();
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
  record_raw_collective();
  runtime::check_mpi_result(
      MPI_Allreduce(&local_parts[1], &global_sum, 1, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(material diagnostic fingerprint sum)");
  record_raw_collective();
  global.fingerprint = {global_xor, global_sum};

  if (request.level != DiagnosticLevel::bounded_state_sample)
    return global;

  require_complete_ownership(mpi, source, control);

  std::vector<std::uint64_t> eligible;
  prepare_aggregation(
      mpi, AggregationPreparationPoint::eligible_counts, control,
      [&] { eligible.resize(static_cast<std::size_t>(mpi.size())); });
  runtime::check_mpi_result(
      MPI_Allgather(&local.eligible_count, 1, MPI_UINT64_T, eligible.data(), 1,
                    MPI_UINT64_T, mpi.comm()),
      "MPI_Allgather(material diagnostic eligible counts)");
  record_raw_collective();
  for (const auto count : eligible) {
    if (count >
        std::numeric_limits<std::uint64_t>::max() - global.eligible_count)
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.global-count-overflow", 0,
                       "global material sample count overflows");
    global.eligible_count += count;
  }
  std::vector<unsigned char> wire;
  std::vector<std::uint64_t> byte_counts;
  std::uint64_t local_bytes{};
  prepare_aggregation(
      mpi, AggregationPreparationPoint::local_sample_wire_and_size_counts,
      control, [&] {
        wire = encode_sample_wire(local.samples);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
        if (control.sample_wire_override)
          wire = *control.sample_wire_override;
#endif
        local_bytes = static_cast<std::uint64_t>(wire.size());
        byte_counts.resize(static_cast<std::size_t>(mpi.size()));
      });
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (control.sample_wire_overflow)
    local_bytes =
        static_cast<std::uint64_t>(std::numeric_limits<int>::max()) + 1U;
  if (control.sample_wire_prefix_overflow) {
    if (mpi.size() == 1)
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.sample-wire-size", 0,
                       "material diagnostic sample exchange is too large");
    local_bytes =
        mpi.rank() == 0
            ? static_cast<std::uint64_t>(std::numeric_limits<int>::max())
            : (mpi.rank() == 1 ? 1U : 0U);
  }
#endif
  runtime::check_mpi_result(MPI_Allgather(&local_bytes, 1, MPI_UINT64_T,
                                          byte_counts.data(), 1, MPI_UINT64_T,
                                          mpi.comm()),
                            "MPI_Allgather(material diagnostic sample sizes)");
  record_raw_collective(RawCollectivePoint::sample_size_exchange);
  std::size_t total_bytes = 0U;
  for (std::size_t rank = 0; rank < byte_counts.size(); ++rank) {
    if (byte_counts[rank] >
            static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        total_bytes >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                static_cast<std::size_t>(byte_counts[rank]))
      collection_error(DiagnosticFailureClass::layout,
                       "material.diagnostics.sample-wire-size",
                       static_cast<int>(rank),
                       "material diagnostic sample exchange is too large");
    total_bytes += static_cast<std::size_t>(byte_counts[rank]);
  }
  std::vector<int> counts;
  std::vector<int> displacements;
  std::vector<unsigned char> all_bytes;
  prepare_aggregation(
      mpi, AggregationPreparationPoint::sample_exchange_buffers, control, [&] {
        counts.resize(byte_counts.size());
        displacements.resize(byte_counts.size());
        std::size_t displacement{};
        for (std::size_t rank = 0; rank < byte_counts.size(); ++rank) {
          counts[rank] = static_cast<int>(byte_counts[rank]);
          displacements[rank] = static_cast<int>(displacement);
          displacement += static_cast<std::size_t>(byte_counts[rank]);
        }
        all_bytes.resize(total_bytes);
      });
  runtime::check_mpi_result(
      MPI_Allgatherv(wire.data(), static_cast<int>(local_bytes), MPI_BYTE,
                     all_bytes.data(), counts.data(), displacements.data(),
                     MPI_BYTE, mpi.comm()),
      "MPI_Allgatherv(material diagnostic samples)");
  record_raw_collective();
  std::vector<DiagnosticSample> all;
  bool duplicate = false;
  int malformed_rank = -1;
  prepare_aggregation(
      mpi, AggregationPreparationPoint::decoded_and_retained_samples, control,
      [&] {
        for (std::size_t rank = 0; rank < counts.size(); ++rank) {
          const auto begin = all_bytes.begin() + displacements[rank];
          const auto end = begin + counts[rank];
          std::vector<unsigned char> segment(begin, end);
          try {
            auto decoded = decode_sample_wire(segment);
            all.insert(all.end(), std::make_move_iterator(decoded.begin()),
                       std::make_move_iterator(decoded.end()));
          } catch (const WireError &) {
            malformed_rank = static_cast<int>(rank);
            return;
          }
        }
        std::sort(
            all.begin(), all.end(), [](const auto &left, const auto &right) {
              return std::tie(left.field_id, left.global_id, left.component) <
                     std::tie(right.field_id, right.global_id, right.component);
            });
        for (std::size_t i = 1; i < all.size(); ++i) {
          if (std::tie(all[i - 1U].field_id, all[i - 1U].global_id,
                       all[i - 1U].component) ==
              std::tie(all[i].field_id, all[i].global_id, all[i].component))
            duplicate = true;
        }
        const std::size_t retain = std::min(request.sample_budget, all.size());
        global.samples.reserve(retain);
        for (std::size_t i = 0; i < retain; ++i)
          global.samples.push_back(std::move(all[i]));
      });
  if (malformed_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "material.diagnostics.sample-wire-malformed",
                     malformed_rank,
                     "material diagnostic sample segment is malformed");
  if (duplicate)
    collection_error(DiagnosticFailureClass::layout,
                     "material.diagnostics.duplicate-owned-sample", 0,
                     "duplicate owned material diagnostic tuple");
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
    override_owned_global_box(runtime::Box3 box, int rank) noexcept {
  diagnostic_test_control.owned_box_override = box;
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_transport_total_count_overflow(int rank) noexcept {
  diagnostic_test_control.transport_total_count_overflow = true;
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
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_allocation_failure(MaterialDiagnosticAllocationPoint point,
                              int rank) noexcept {
  diagnostic_test_control.allocation_failure =
      static_cast<AggregationPreparationPoint>(point);
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_sample_wire_prefix_overflow() noexcept {
  diagnostic_test_control.sample_wire_prefix_overflow = true;
  diagnostic_test_control.rank = -1;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    inject_provider_key_value_difference(int rank) noexcept {
  diagnostic_test_control.provider_key_value_difference = true;
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    override_request_wire_bytes(const unsigned char *bytes, std::size_t size,
                                int rank) {
  diagnostic_test_control.request_wire_override =
      std::vector<unsigned char>(bytes, bytes + size);
  diagnostic_test_control.rank = rank;
}
void test::MaterialDensityTransportDiagnosticsTestAccess::
    override_sample_wire_bytes(const unsigned char *bytes, std::size_t size,
                               int rank) {
  diagnostic_test_control.sample_wire_override =
      std::vector<unsigned char>(bytes, bytes + size);
  diagnostic_test_control.rank = rank;
}
std::vector<unsigned char>
test::MaterialDensityTransportDiagnosticsTestAccess::request_wire_bytes(
    const DiagnosticRequest &request) {
  return request_key(request);
}
bool test::MaterialDensityTransportDiagnosticsTestAccess::request_wire_is_valid(
    const std::vector<unsigned char> &bytes) noexcept {
  try {
    static_cast<void>(decode_request_wire(bytes));
    return true;
  } catch (...) {
    return false;
  }
}
std::vector<unsigned char>
test::MaterialDensityTransportDiagnosticsTestAccess::request_wire_round_trip(
    const std::vector<unsigned char> &bytes) {
  return encode_request_wire(decode_request_wire(bytes));
}
std::vector<unsigned char>
test::MaterialDensityTransportDiagnosticsTestAccess::sample_wire_bytes(
    const DiagnosticSample &sample) {
  return encode_sample_wire({sample});
}
bool test::MaterialDensityTransportDiagnosticsTestAccess::sample_wire_is_valid(
    const std::vector<unsigned char> &bytes) noexcept {
  try {
    static_cast<void>(decode_sample_wire(bytes));
    return true;
  } catch (...) {
    return false;
  }
}
DiagnosticSample
test::MaterialDensityTransportDiagnosticsTestAccess::decode_single_sample_wire(
    const std::vector<unsigned char> &bytes) {
  auto decoded = decode_sample_wire(bytes);
  if (decoded.size() != 1U)
    throw std::invalid_argument("sample wire does not contain one sample");
  return std::move(decoded.front());
}
std::vector<unsigned char>
test::MaterialDensityTransportDiagnosticsTestAccess::provider_key_bytes(
    const flow::MaterialDensityDiagnosticSource &source) {
  return provider_key(source, {});
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
    record_raw_collective();
    code = broadcast_failure_code(mpi, preflight_rank, std::move(code));
    collection_error(static_cast<DiagnosticFailureClass>(class_value),
                     std::move(code), preflight_rank,
                     "material diagnostic preflight failed");
  }
  require_request_agreement(mpi, request, control);
  require_provider_agreement(mpi, source, control);

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
    record_raw_collective();
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
    record_raw_collective();
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
