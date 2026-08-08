// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_material_density_piso.hpp"

#include "hundun/flow_material_density_piso.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <charconv>
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

constexpr std::string_view kModuleId =
    "hundun.flow.fixed_step_material_density";
constexpr std::string_view kInstanceId = "primary";
constexpr std::string_view kPhase = "material-density.attempt-result";
constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(
        DiagnosticCapability::bounded_state_sample) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
constexpr std::uint8_t kNoInjectedFault = 0U;
constexpr std::uint8_t kProviderAgreementInjectedFault = 1U;
constexpr std::uint8_t kCellExactCoverInjectedFault = 2U;
constexpr std::uint8_t kFaceExactCoverInjectedFault = 3U;
constexpr std::uint8_t kSampleSendPreparationInjectedFault = 4U;
constexpr std::uint8_t kSampleReceivePreparationInjectedFault = 5U;
constexpr std::uint8_t kSampleWireMalformedInjectedFault = 6U;
constexpr std::uint8_t kRecordValidationInjectedFault = 7U;

std::uint8_t injected_fault{kNoInjectedFault};

bool fault(std::uint8_t value) noexcept {
  return injected_fault == value;
}
#endif

[[noreturn]] void collection_error(DiagnosticFailureClass classification,
                                   std::string code, int rank,
                                   std::string message) {
  throw DiagnosticCollectionError(classification, std::move(code), rank,
                                  std::move(message));
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

class CompensatedSum final {
public:
  void add(double value) noexcept {
    const double next = sum_ + value;
    correction_ += std::abs(sum_) >= std::abs(value)
                       ? (sum_ - next) + value
                       : (value - next) + sum_;
    sum_ = next;
  }
  double value() const noexcept { return sum_ + correction_; }

private:
  double sum_{};
  double correction_{};
};

struct Observation final {
  double density_min{std::numeric_limits<double>::infinity()};
  double density_max{-std::numeric_limits<double>::infinity()};
  double total_mass{};
  std::array<double, 3> total_momentum{};
  DiagnosticFingerprintParts fingerprint{};
  std::uint64_t eligible{};
  std::vector<DiagnosticSample> samples;
};

class AgreementWriter final {
public:
  void u8(std::uint8_t value) { bytes_.push_back(value); }
  void u32(std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }
  void u64(std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }
  void text(std::string_view value) {
    u64(static_cast<std::uint64_t>(value.size()));
    bytes_.insert(bytes_.end(), value.begin(), value.end());
  }
  std::vector<unsigned char> finish() && { return std::move(bytes_); }

private:
  std::vector<unsigned char> bytes_;
};

constexpr std::uint32_t kSampleWireSchemaV1 = 1U;
constexpr std::size_t kSampleWireHeaderBytes = 12U;
constexpr std::size_t kSampleWireItemBytes = 25U;
static_assert(sizeof(std::uint8_t) == 1U);
static_assert(sizeof(std::uint32_t) == 4U);
static_assert(sizeof(std::uint64_t) == 8U);
constexpr std::array<std::string_view, 5> kSampleFieldIds{
    "face_mass_flux", "face_velocity", "pi", "rho", "velocity"};
constexpr std::array<std::uint32_t, 5> kSampleFieldComponents{1U, 3U, 1U,
                                                              1U, 3U};
constexpr std::array<std::string_view, 5> kSampleFieldUnits{
    "kg/s", "m/s", "Pa", "kg/m3", "m/s"};

struct SampleWireItem final {
  std::uint32_t field{};
  std::uint64_t global_id{};
  std::uint32_t component{};
  std::uint64_t value_bits{};
  DiagnosticValueStatus status{DiagnosticValueStatus::finite};
};

DiagnosticValueStatus status_for_bits(std::uint64_t value) noexcept {
  const auto exponent = (value >> 52U) & UINT64_C(0x7ff);
  const auto fraction = value & UINT64_C(0x000fffffffffffff);
  if (exponent != UINT64_C(0x7ff))
    return DiagnosticValueStatus::finite;
  if (fraction == 0U)
    return (value >> 63U) == 0U
               ? DiagnosticValueStatus::positive_infinity
               : DiagnosticValueStatus::negative_infinity;
  return (fraction & UINT64_C(0x0008000000000000)) != 0U
             ? DiagnosticValueStatus::quiet_nan
             : DiagnosticValueStatus::signaling_nan;
}

void validate_sample_wire_item(const SampleWireItem &item) {
  if (item.field >= kSampleFieldIds.size() ||
      item.component >= kSampleFieldComponents[item.field] ||
      static_cast<std::uint8_t>(item.status) >
          static_cast<std::uint8_t>(DiagnosticValueStatus::signaling_nan) ||
      status_for_bits(item.value_bits) != item.status)
    throw runtime::Error(
        "material flow diagnostic sample representation is invalid");
}

class SampleWireWriter final {
public:
  explicit SampleWireWriter(std::size_t capacity) { bytes_.reserve(capacity); }

  void u8(std::uint8_t value) { bytes_.push_back(value); }
  void u32(std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }
  void u64(std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      bytes_.push_back(static_cast<unsigned char>(value >> shift));
  }
  std::vector<unsigned char> finish() && { return std::move(bytes_); }

private:
  std::vector<unsigned char> bytes_;
};

class SampleWireReader final {
public:
  SampleWireReader(const unsigned char *data, std::size_t size) noexcept
      : data_(data), size_(size) {}

  std::uint8_t u8() { return integer<std::uint8_t>(); }
  std::uint32_t u32() { return integer<std::uint32_t>(); }
  std::uint64_t u64() { return integer<std::uint64_t>(); }
  std::size_t remaining() const noexcept { return size_ - offset_; }
  void require_end() const {
    if (remaining() != 0U)
      throw runtime::Error(
          "material flow diagnostic sample wire has trailing bytes");
  }

private:
  template <class Integer> Integer integer() {
    if (sizeof(Integer) > remaining())
      throw runtime::Error(
          "material flow diagnostic sample wire is truncated");
    std::uint64_t value{};
    for (std::size_t byte = 0; byte < sizeof(Integer); ++byte)
      value |= static_cast<std::uint64_t>(data_[offset_ + byte])
               << (8U * byte);
    offset_ += sizeof(Integer);
    return static_cast<Integer>(value);
  }

  const unsigned char *data_{};
  std::size_t size_{};
  std::size_t offset_{};
};

std::vector<unsigned char>
encode_sample_payload(const std::vector<SampleWireItem> &items) {
  if (items.size() > kMaximumStateSamplesV1)
    throw runtime::Error(
        "material flow diagnostic sample wire count is invalid");
  const std::size_t encoded_size =
      kSampleWireHeaderBytes + items.size() * kSampleWireItemBytes;
  SampleWireWriter writer(encoded_size);
  writer.u32(kSampleWireSchemaV1);
  writer.u64(static_cast<std::uint64_t>(items.size()));
  for (const auto &item : items) {
    validate_sample_wire_item(item);
    writer.u32(item.field);
    writer.u64(item.global_id);
    writer.u32(item.component);
    writer.u64(item.value_bits);
    writer.u8(static_cast<std::uint8_t>(item.status));
  }
  return std::move(writer).finish();
}

void decode_sample_payload(const unsigned char *data, std::size_t size,
                           std::vector<SampleWireItem> &output) {
  const std::size_t original_size = output.size();
  try {
    SampleWireReader reader(data, size);
    if (reader.u32() != kSampleWireSchemaV1)
      throw runtime::Error(
          "material flow diagnostic sample wire schema is invalid");
    const auto count = reader.u64();
    if (count > static_cast<std::uint64_t>(kMaximumStateSamplesV1) ||
        count > static_cast<std::uint64_t>(reader.remaining() /
                                           kSampleWireItemBytes) ||
        count > static_cast<std::uint64_t>(output.capacity() - output.size()))
      throw runtime::Error(
          "material flow diagnostic sample wire count is invalid");
    for (std::uint64_t index = 0; index < count; ++index) {
      SampleWireItem item;
      item.field = reader.u32();
      item.global_id = reader.u64();
      item.component = reader.u32();
      item.value_bits = reader.u64();
      item.status = static_cast<DiagnosticValueStatus>(reader.u8());
      validate_sample_wire_item(item);
      output.push_back(item);
    }
    reader.require_end();
  } catch (...) {
    output.resize(original_size);
    throw;
  }
}

bool selected(const DiagnosticRequest &request, std::string_view id) {
  return request.selected_fields.empty() ||
         std::binary_search(request.selected_fields.begin(),
                            request.selected_fields.end(), id);
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

void require_request(
  const flow::MaterialDensityFlowDiagnosticSource &source,
    const DiagnosticRequest &request, DiagnosticScope expected_scope) {
  bool source_valid = true;
  int rank{};
  std::uint64_t step{};
  double time{};
  try {
    rank = source.relative_rank();
    step = source.committed_step();
    time = source.committed_time_s();
  } catch (const runtime::Error &error) {
    source_valid = false;
    if (std::string_view(error.what()) ==
        "material flow diagnostic source is stale")
      collection_error(DiagnosticFailureClass::invalid_input,
                       "flow.diagnostics.stale-source", -1,
                       "material flow diagnostic source is stale");
  }
  if (!source_valid)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.layout", -1,
                     "material flow diagnostic provider is invalid");
  if (request.scope != expected_scope || request.frame.rank != rank ||
      request.frame.step != step || bits(request.frame.time_s) != bits(time) ||
      request.frame.phase != kPhase)
    collection_error(DiagnosticFailureClass::invalid_request,
                     "flow.diagnostics.frame", -1,
                     "material flow diagnostic frame does not match");
  const auto ids = diagnostic_fingerprint_field_ids(source);
  for (const auto field : request.selected_fields)
    if (!std::binary_search(ids.begin(), ids.end(), field))
      collection_error(DiagnosticFailureClass::invalid_request,
                       "flow.diagnostics.unknown-field", -1,
                       "unknown material flow diagnostic field");
  try {
    validate(request, describe_diagnostics(source));
  } catch (const std::exception &) {
    collection_error(DiagnosticFailureClass::invalid_request,
                     "flow.diagnostics.frame", -1,
                     "material flow diagnostic request is invalid");
  }
}

Observation observe(const flow::MaterialDensityFlowDiagnosticSource &source,
                    const DiagnosticRequest &request) {
  Observation result;
  const auto ids = diagnostic_fingerprint_field_ids(source);
  std::vector<double> density;
  if (request.level == DiagnosticLevel::summary ||
      request.level == DiagnosticLevel::invariants)
    density.resize(source.owned_cell_count());
  DiagnosticFingerprintAccumulator fingerprint;
  for (std::size_t field = 0; field < ids.size(); ++field) {
    const std::size_t items = source.field_item_count(field);
    const std::size_t components = source.field_component_count(field);
    for (std::size_t item = 0; item < items; ++item) {
      const auto global = source.field_global_id(field, item);
      for (std::size_t component = 0; component < components; ++component) {
        const double value = source.field_value(field, item, component);
        const auto represented = describe_fp64(value);
        fingerprint.add(ids[field], global,
                        static_cast<std::uint32_t>(component), represented);
        if (request.level == DiagnosticLevel::bounded_state_sample &&
            selected(request, ids[field])) {
          if (result.eligible == std::numeric_limits<std::uint64_t>::max())
            collection_error(DiagnosticFailureClass::layout,
                             "flow.diagnostics.layout", -1,
                             "material flow diagnostic sample count wraps");
          ++result.eligible;
          if (result.samples.size() < request.sample_budget)
            result.samples.push_back(
                {std::string(ids[field]), global,
                 static_cast<std::uint32_t>(component),
                 std::string(source.field_unit(field)), represented});
        }
        if (field == 3U && component == 0U && !density.empty()) {
          density[item] = value;
          result.density_min = std::min(result.density_min, value);
          result.density_max = std::max(result.density_max, value);
        }
      }
    }
  }
  if (request.level == DiagnosticLevel::summary) {
    CompensatedSum mass;
    std::array<CompensatedSum, 3> momentum;
    for (std::size_t cell = 0; cell < density.size(); ++cell) {
      const double volume = source.cell_volume_m3(cell);
      mass.add(volume * density[cell]);
      for (std::size_t component = 0; component < 3U; ++component)
        momentum[component].add(
            volume * density[cell] * source.field_value(4U, cell, component));
    }
    result.total_mass = mass.value();
    for (std::size_t component = 0; component < 3U; ++component)
      result.total_momentum[component] = momentum[component].value();
  }
  result.fingerprint = fingerprint.parts();
  return result;
}

DiagnosticFailure mapped_failure(const flow::MaterialDensityStepAttemptReport &r,
                                 int rank) {
  using R = flow::StepFailureReason;
  switch (r.flow().reason) {
  case R::none:
    return {};
  case R::invalid_input:
    return {DiagnosticFailureClass::invalid_input, "flow.invalid-input", rank};
  case R::momentum_linear_solve:
    return {DiagnosticFailureClass::non_convergence,
            "flow.momentum-linear-solve", rank};
  case R::pressure_linear_solve:
    return {DiagnosticFailureClass::non_convergence,
            "flow.pressure-linear-solve", rank};
  case R::non_finite_trial:
    return {DiagnosticFailureClass::non_finite_state,
            "flow.non-finite-trial", rank};
  case R::boundary_backflow:
    return {DiagnosticFailureClass::boundary, "flow.boundary-backflow", rank};
  case R::transport_failure:
    if (r.material_failure_reason() ==
        flow::MaterialTransportFailureReason::non_positive_density)
      return {DiagnosticFailureClass::non_positive_state,
              "flow.transport-non-positive-density", rank};
    return {DiagnosticFailureClass::non_finite_state,
            "flow.transport-non-finite-state", rank};
  case R::final_momentum_residual:
    return {DiagnosticFailureClass::non_convergence,
            "flow.final-momentum-residual", rank};
  case R::final_transport_residual:
    return {DiagnosticFailureClass::non_convergence,
            "flow.final-transport-residual", rank};
  case R::final_conservation_defect:
    return {DiagnosticFailureClass::conservation,
            "flow.final-conservation-defect", rank};
  case R::final_continuity_residual:
    return {DiagnosticFailureClass::non_convergence,
            "flow.final-continuity-residual", rank};
  case R::final_pressure_residual:
    return {DiagnosticFailureClass::non_convergence,
            "flow.final-pressure-residual", rank};
  case R::collective_operation:
    return {DiagnosticFailureClass::collective_operation,
            "flow.collective-operation", rank};
  case R::density_closure_failure:
    collection_error(DiagnosticFailureClass::capability,
                     "closure.diagnostics.capability", rank,
                     "material diagnostics cannot represent closure failure");
  }
  return {DiagnosticFailureClass::invalid_input, "flow.invalid-input", rank};
}

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.layout", -1,
                     "material flow diagnostic counter wraps");
  return left + right;
}

DiagnosticRecord build_record(
    const flow::MaterialDensityFlowDiagnosticSource &source,
    const DiagnosticRequest &request, const Observation &observation,
    DiagnosticScope scope) {
  const auto descriptor = describe_diagnostics(source);
  const auto &parent = source.report();
  const auto &flow_report = parent.flow();
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
  record.phase = std::string(kPhase);
  const bool committed = flow_report.disposition ==
                         flow::StepAttemptDisposition::committed;
  record.status = committed ? DiagnosticStatus::ok : DiagnosticStatus::failed;
  record.failure = mapped_failure(
      parent, scope == DiagnosticScope::collective && !committed
                  ? flow_report.lowest_failing_rank
                  : -1);
  const auto ids = diagnostic_fingerprint_field_ids(source);
  for (std::size_t field = 0; field < ids.size(); ++field) {
    const std::string layout(
        source.field_entity(field) == flow::MaterialDensityDiagnosticEntity::face
            ? (scope == DiagnosticScope::local
                   ? source.owned_face_layout_fingerprint()
                   : source.global_face_layout_fingerprint())
            : (scope == DiagnosticScope::local
                   ? source.owned_cell_layout_fingerprint()
                   : source.global_cell_layout_fingerprint()));
    record.identities.push_back({"field." + std::string(ids[field]), layout,
                                 std::nullopt, std::nullopt, std::nullopt});
  }
  record.identities.push_back({"flow.attempt", std::nullopt,
                               parent.attempt_identity(), std::nullopt,
                               std::nullopt});
  record.identities.push_back(
      {"layout.cells",
       std::string(scope == DiagnosticScope::local
                       ? source.owned_cell_layout_fingerprint()
                       : source.global_cell_layout_fingerprint()),
       std::nullopt, std::nullopt, std::nullopt});
  record.identities.push_back(
      {"layout.faces",
       std::string(scope == DiagnosticScope::local
                       ? source.owned_face_layout_fingerprint()
                       : source.global_face_layout_fingerprint()),
       std::nullopt, std::nullopt, std::nullopt});
  std::optional<std::uint64_t> finalization;
  if (parent.material_report_available())
    finalization = parent.material_report().finalization_identity();
  record.identities.push_back({"material.finalization", std::nullopt,
                               finalization, std::nullopt, std::nullopt});
  DiagnosticFingerprintAccumulator fingerprint;
  fingerprint.combine(observation.fingerprint);
  record.state_fingerprint = fingerprint.finish();

  const auto unavailable = [] { return DiagnosticFp64{}; };
  const auto material = parent.material_report_available()
                            ? &parent.material_report()
                            : nullptr;
  if (request.level == DiagnosticLevel::summary) {
    record.metrics.push_back({"attempt.dt", DiagnosticMetricKind::state_summary,
                              "s", describe_fp64(flow_report.attempted_dt_s)});
    record.metrics.push_back(
        {"continuity.residual", DiagnosticMetricKind::residual, "1",
         parent.final_continuity_residual_available()
             ? describe_fp64(flow_report.final_continuity_normalized_l2)
             : unavailable()});
    record.metrics.push_back({"density.maximum",
                              DiagnosticMetricKind::state_summary, "kg/m3",
                              describe_fp64(observation.density_max)});
    record.metrics.push_back({"density.minimum",
                              DiagnosticMetricKind::state_summary, "kg/m3",
                              describe_fp64(observation.density_min)});
    record.metrics.push_back(
        {"mass.relative-defect", DiagnosticMetricKind::conservation, "1",
         parent.mass_conservation_available()
             ? describe_fp64(flow_report.final_mass_relative_conservation_defect)
             : unavailable()});
    record.metrics.push_back({"mass.total", DiagnosticMetricKind::conservation,
                              "kg", describe_fp64(observation.total_mass)});
    record.metrics.push_back(
        {"material.density.residual", DiagnosticMetricKind::residual, "1",
         material != nullptr && material->density_residual_available()
             ? describe_fp64(material->density_normalized_l2())
             : unavailable()});
    record.metrics.push_back(
        {"material.mass.relative-defect", DiagnosticMetricKind::conservation,
         "1", material != nullptr && material->mass_conservation_available()
                  ? describe_fp64(material->mass_relative_conservation_defect())
                  : unavailable()});
    constexpr std::array<char, 3> components{'x', 'y', 'z'};
    for (std::size_t c = 0; c < 3U; ++c) {
      const std::string root = std::string("momentum.") + components[c];
      record.metrics.push_back(
          {root + ".relative-defect", DiagnosticMetricKind::conservation, "1",
           parent.momentum_conservation_availability()[c] != 0U
               ? describe_fp64(
                     flow_report.final_momentum_relative_conservation_defect[c])
               : unavailable()});
      record.metrics.push_back(
          {root + ".residual", DiagnosticMetricKind::residual, "1",
           parent.final_momentum_residual_availability()[c] != 0U
               ? describe_fp64(flow_report.final_momentum_normalized_l2[c])
               : unavailable()});
      record.metrics.push_back(
          {root + ".total", DiagnosticMetricKind::conservation, "kg*m/s",
           describe_fp64(observation.total_momentum[c])});
    }
    record.metrics.push_back(
        {"pressure.residual", DiagnosticMetricKind::residual, "1",
         parent.final_pressure_residual_available()
             ? describe_fp64(parent.final_pressure_normalized_residual())
             : unavailable()});
    for (std::size_t index = 0; index < parent.material_field_count(); ++index) {
      const std::string root = indexed_id("transport.s", index);
      const bool conservation =
          material != nullptr &&
          index < material->transport_conservation_availability().size() &&
          material->transport_conservation_availability()[index] != 0U;
      const bool residual =
          material != nullptr &&
          index < material->transport_residual_availability().size() &&
          material->transport_residual_availability()[index] != 0U;
      record.metrics.push_back(
          {root + ".relative-defect", DiagnosticMetricKind::conservation, "1",
           conservation
               ? describe_fp64(
                     material->transport_relative_conservation_defect()[index])
               : unavailable()});
      record.metrics.push_back(
          {root + ".residual", DiagnosticMetricKind::residual, "1",
           residual ? describe_fp64(material->transport_normalized_l2()[index])
                    : unavailable()});
    }
  } else if (request.level == DiagnosticLevel::invariants) {
    auto add = [&](std::string id, std::string unit, double observed,
                   double limit, InvariantRelation relation) {
      DiagnosticInvariant item{std::move(id), std::move(unit),
                               describe_fp64(observed),
                               relation == InvariantRelation::positive
                                   ? DiagnosticFp64{}
                                   : describe_fp64(limit),
                               relation, false};
      item.passed = evaluate_invariant(item);
      record.invariants.push_back(std::move(item));
    };
    auto available = [&](std::string id, bool value) {
      add(std::move(id), "1", value ? 1.0 : 0.0, 1.0,
          InvariantRelation::equal);
    };
    if (parent.final_continuity_residual_available())
      add("continuity.residual", "1",
          flow_report.final_continuity_normalized_l2, 1.0e-10,
          InvariantRelation::less_equal);
    available("continuity.residual.available",
              parent.final_continuity_residual_available());
    add("density.positive", "kg/m3", observation.density_min, 0.0,
        InvariantRelation::positive);
    add("flux.final-corrected", "1",
        parent.flux_provenance() == flow::MaterialFluxProvenance::final_corrected
            ? 1.0
            : 0.0,
        1.0, InvariantRelation::equal);
    if (parent.mass_conservation_available())
      add("mass.relative-defect", "1",
          flow_report.final_mass_relative_conservation_defect, 5.0e-11,
          InvariantRelation::less_equal);
    available("mass.relative-defect.available",
              parent.mass_conservation_available());
    if (material != nullptr && material->density_residual_available())
      add("material.density.residual", "1",
          material->density_normalized_l2(), 1.0e-10,
          InvariantRelation::less_equal);
    available("material.density.residual.available",
              material != nullptr && material->density_residual_available());
    add("material.finalized", "1",
        material != nullptr && material->disposition() ==
                                   flow::MaterialTransportDisposition::finalized
            ? 1.0
            : 0.0,
        1.0, InvariantRelation::equal);
    if (material != nullptr && material->mass_conservation_available())
      add("material.mass.relative-defect", "1",
          material->mass_relative_conservation_defect(), 5.0e-11,
          InvariantRelation::less_equal);
    available("material.mass.relative-defect.available",
              material != nullptr && material->mass_conservation_available());
    constexpr std::array<char, 3> components{'x', 'y', 'z'};
    for (std::size_t c = 0; c < 3U; ++c) {
      const std::string root = std::string("momentum.") + components[c];
      const bool conservation =
          parent.momentum_conservation_availability()[c] != 0U;
      if (conservation)
        add(root + ".relative-defect", "1",
            flow_report.final_momentum_relative_conservation_defect[c],
            5.0e-11, InvariantRelation::less_equal);
      available(root + ".relative-defect.available", conservation);
      const bool residual =
          parent.final_momentum_residual_availability()[c] != 0U;
      if (residual)
        add(root + ".residual", "1",
            flow_report.final_momentum_normalized_l2[c], 1.0e-9,
            InvariantRelation::less_equal);
      available(root + ".residual.available", residual);
    }
    add("pressure.correctors", "count",
        static_cast<double>(flow_report.pressure_corrector_count), 2.0,
        InvariantRelation::equal);
    if (parent.final_pressure_residual_available())
      add("pressure.residual", "1",
          parent.final_pressure_normalized_residual(), 1.0,
          InvariantRelation::less_equal);
    available("pressure.residual.available",
              parent.final_pressure_residual_available());
    for (std::size_t index = 0; index < parent.material_field_count(); ++index) {
      const std::string root = indexed_id("transport.s", index);
      const bool conservation =
          material != nullptr &&
          index < material->transport_conservation_availability().size() &&
          material->transport_conservation_availability()[index] != 0U;
      if (conservation)
        add(root + ".relative-defect", "1",
            material->transport_relative_conservation_defect()[index],
            5.0e-11, InvariantRelation::less_equal);
      available(root + ".relative-defect.available", conservation);
      const bool residual =
          material != nullptr &&
          index < material->transport_residual_availability().size() &&
          material->transport_residual_availability()[index] != 0U;
      if (residual)
        add(root + ".residual", "1",
            material->transport_normalized_l2()[index], 1.0e-9,
            InvariantRelation::less_equal);
      available(root + ".residual.available", residual);
    }
  } else if (request.level == DiagnosticLevel::counters) {
    auto sum_reports = [&](auto getter) {
      std::uint64_t value{};
      for (const auto &item : flow_report.momentum.components)
        value = checked_add(value, getter(item));
      return value;
    };
    auto sum_pressure = [&](auto getter) {
      std::uint64_t value{};
      for (const auto &item : flow_report.pressure)
        value = checked_add(value, getter(item));
      return value;
    };
    record.counters.push_back(
        {"attempt.identity", "count", parent.attempt_identity()});
    record.counters.push_back(
        {"material.fields", "count", parent.material_field_count()});
    record.counters.push_back(
        {"material.finalization.identity", "count",
         material == nullptr ? 0U : material->finalization_identity()});
    record.counters.push_back(
        {"pressure.correctors", "count", flow_report.pressure_corrector_count});
    record.counters.push_back(
        {"solver.momentum.iterations", "count",
         sum_reports([](const auto &x) { return x.iterations; })});
    record.counters.push_back(
        {"solver.momentum.matvec", "count",
         sum_reports([](const auto &x) { return x.matvec_count; })});
    record.counters.push_back(
        {"solver.momentum.preconditioner", "count",
         sum_reports([](const auto &x) { return x.preconditioner_apply_count; })});
    record.counters.push_back(
        {"solver.momentum.reductions", "count",
         sum_reports([](const auto &x) { return x.global_reduction_count; })});
    record.counters.push_back(
        {"solver.pressure.iterations", "count",
         sum_pressure([](const auto &x) { return x.iterations; })});
    record.counters.push_back(
        {"solver.pressure.matvec", "count",
         sum_pressure([](const auto &x) { return x.matvec_count; })});
    record.counters.push_back(
        {"solver.pressure.preconditioner", "count",
         sum_pressure([](const auto &x) { return x.preconditioner_apply_count; })});
    record.counters.push_back(
        {"solver.pressure.reductions", "count",
         sum_pressure([](const auto &x) { return x.global_reduction_count; })});
  } else {
    record.sample_budget = request.sample_budget;
    record.eligible_sample_count = observation.eligible;
    record.samples_truncated = observation.eligible > observation.samples.size();
    record.samples = observation.samples;
  }
  return record;
}

void submit_local(DiagnosticRecord record, DiagnosticSink &sink) {
  try {
    sink.submit(record);
  } catch (...) {
    collection_error(DiagnosticFailureClass::sink_failure,
                     "diagnostics.sink.submit", -1,
                     "diagnostic sink submission failed");
  }
}

int lowest_rank(const runtime::MpiContext &mpi, bool failure,
                std::string_view operation) {
  int local = failure ? mpi.rank() : mpi.size();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      operation);
  return lowest == mpi.size() ? -1 : lowest;
}

std::vector<unsigned char> request_agreement_key(
    const DiagnosticRequest &request) {
  AgreementWriter writer;
  writer.u8(static_cast<std::uint8_t>(request.level));
  writer.u8(static_cast<std::uint8_t>(request.scope));
  writer.u64(request.frame.step);
  writer.u64(bits(request.frame.time_s));
  writer.text(request.frame.phase);
  writer.u64(static_cast<std::uint64_t>(request.selected_fields.size()));
  for (const auto field : request.selected_fields)
    writer.text(field);
  writer.u64(static_cast<std::uint64_t>(request.sample_budget));
  return std::move(writer).finish();
}

std::vector<unsigned char> provider_agreement_key(
    const flow::MaterialDensityFlowDiagnosticSource &source) {
  AgreementWriter writer;
  const auto descriptor = describe_diagnostics(source);
  writer.u32(descriptor.schema_version);
  writer.u32(static_cast<std::uint32_t>(descriptor.module_kind));
  writer.text(descriptor.module_id);
  writer.text(descriptor.instance_id);
  writer.u32(descriptor.capabilities);
  writer.text(source.global_cell_layout_fingerprint());
  writer.text(source.global_face_layout_fingerprint());
  const auto extent = source.global_cell_extent();
  writer.u32(static_cast<std::uint32_t>(extent.x));
  writer.u32(static_cast<std::uint32_t>(extent.y));
  writer.u32(static_cast<std::uint32_t>(extent.z));
  writer.u64(static_cast<std::uint64_t>(source.fingerprint_field_count()));
  for (std::size_t field = 0; field < source.fingerprint_field_count(); ++field) {
    writer.text(source.fingerprint_field_id(field));
    writer.text(source.field_unit(field));
    writer.u8(static_cast<std::uint8_t>(source.field_entity(field)));
    writer.u64(static_cast<std::uint64_t>(source.field_component_count(field)));
  }

  const auto &parent = source.report();
  const auto &report = parent.flow();
  writer.u8(static_cast<std::uint8_t>(report.disposition));
  writer.u8(static_cast<std::uint8_t>(report.reason));
  writer.u32(static_cast<std::uint32_t>(report.lowest_failing_rank));
  writer.u64(report.pressure_corrector_count);
  writer.u64(bits(report.attempted_dt_s));
  writer.u64(bits(report.suggested_dt_s));
  const auto append_solve = [&](const linear::SolveReport &solve) {
    writer.u8(static_cast<std::uint8_t>(solve.reason));
    writer.u64(solve.iterations);
    writer.u64(bits(solve.initial_residual));
    writer.u64(bits(solve.recursive_residual));
    writer.u64(bits(solve.final_residual));
    writer.u64(solve.matvec_count);
    writer.u64(solve.preconditioner_apply_count);
    writer.u64(solve.global_reduction_count);
    writer.u32(static_cast<std::uint32_t>(solve.lowest_failing_rank));
  };
  for (const auto &solve : report.momentum.components)
    append_solve(solve);
  for (const auto &solve : report.pressure)
    append_solve(solve);
  writer.u64(bits(report.final_continuity_normalized_l2));
  writer.u64(bits(report.final_pressure_residual_l2));
  for (const auto value : report.final_momentum_normalized_l2)
    writer.u64(bits(value));
  writer.u64(static_cast<std::uint64_t>(
      report.final_transport_normalized_l2.size()));
  for (const auto value : report.final_transport_normalized_l2)
    writer.u64(bits(value));
  writer.u64(bits(report.final_mass_relative_conservation_defect));
  for (const auto value : report.final_momentum_relative_conservation_defect)
    writer.u64(bits(value));
  writer.u64(static_cast<std::uint64_t>(
      report.final_transport_relative_conservation_defect.size()));
  for (const auto value : report.final_transport_relative_conservation_defect)
    writer.u64(bits(value));

  writer.u8(parent.material_report_available() ? 1U : 0U);
  writer.u8(static_cast<std::uint8_t>(parent.material_failure_reason()));
  writer.u64(parent.material_field_count());
  writer.u32(parent.shared_face_mass_flux_field());
  writer.u8(static_cast<std::uint8_t>(parent.flux_provenance()));
  writer.u64(parent.attempt_identity());
  writer.u8(parent.final_continuity_residual_available() ? 1U : 0U);
  writer.u8(parent.final_pressure_residual_available() ? 1U : 0U);
  writer.u64(bits(parent.final_pressure_normalized_residual()));
  for (const auto available : parent.final_momentum_residual_availability())
    writer.u8(available);
  writer.u8(parent.mass_conservation_available() ? 1U : 0U);
  for (const auto available : parent.momentum_conservation_availability())
    writer.u8(available);
  if (parent.material_report_available()) {
    const auto &material = parent.material_report();
    writer.u8(static_cast<std::uint8_t>(material.disposition()));
    writer.u8(static_cast<std::uint8_t>(material.reason()));
    writer.u32(static_cast<std::uint32_t>(material.lowest_failing_rank()));
    writer.u8(static_cast<std::uint8_t>(material.flux_provenance()));
    writer.u64(material.attempt_identity());
    writer.u64(material.finalization_identity());
    writer.u32(material.shared_face_mass_flux_field());
    writer.u8(material.density_residual_available() ? 1U : 0U);
    writer.u64(bits(material.density_normalized_l2()));
    writer.u64(static_cast<std::uint64_t>(
        material.transport_residual_availability().size()));
    for (const auto available : material.transport_residual_availability())
      writer.u8(available);
    writer.u64(static_cast<std::uint64_t>(
        material.transport_normalized_l2().size()));
    for (const auto value : material.transport_normalized_l2())
      writer.u64(bits(value));
    writer.u8(material.mass_conservation_available() ? 1U : 0U);
    writer.u64(bits(material.mass_relative_conservation_defect()));
    writer.u64(static_cast<std::uint64_t>(
        material.transport_conservation_availability().size()));
    for (const auto available : material.transport_conservation_availability())
      writer.u8(available);
    writer.u64(static_cast<std::uint64_t>(
        material.transport_relative_conservation_defect().size()));
    for (const auto value : material.transport_relative_conservation_defect())
      writer.u64(bits(value));
  }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(kProviderAgreementInjectedFault))
    writer.u8(0xffU);
#endif
  return std::move(writer).finish();
}

void require_byte_agreement(const runtime::MpiContext &mpi,
                            const std::vector<unsigned char> &local,
                            DiagnosticFailureClass classification,
                            std::string code, std::string_view label) {
  std::uint64_t reference_size =
      mpi.rank() == 0 ? static_cast<std::uint64_t>(local.size()) : 0U;
  runtime::check_mpi_result(
      MPI_Bcast(&reference_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(flow diagnostic agreement size)");
  bool preparation_failed =
      reference_size > static_cast<std::uint64_t>(
                           std::numeric_limits<int>::max());
  std::vector<unsigned char> reference;
  try {
    if (!preparation_failed)
      reference.resize(static_cast<std::size_t>(reference_size));
  } catch (...) {
    preparation_failed = true;
  }
  const int preparation_rank = lowest_rank(
      mpi, preparation_failed,
      "MPI_Allreduce(flow diagnostic agreement preparation)");
  if (preparation_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.agreement-preparation",
                     preparation_rank,
                     "collective material flow diagnostic agreement preparation failed");
  if (mpi.rank() == 0)
    std::copy(local.begin(), local.end(), reference.begin());
  runtime::check_mpi_result(
      MPI_Bcast(reference.data(), static_cast<int>(reference.size()), MPI_BYTE,
                0, mpi.comm()),
      "MPI_Bcast(flow diagnostic agreement bytes)");
  const int mismatch = lowest_rank(
      mpi, local != reference, "MPI_Allreduce(flow diagnostic agreement)");
  if (mismatch >= 0)
    collection_error(classification, std::move(code), mismatch,
                     std::string(label) + " differ across ranks");
}

std::size_t global_face_count(
    const flow::MaterialDensityFlowDiagnosticSource &source) {
  constexpr std::string_view prefix = "face.f64.global.";
  const auto fingerprint = source.global_face_layout_fingerprint();
  if (fingerprint.substr(0U, prefix.size()) != prefix)
    throw runtime::Error("material flow global face layout is invalid");
  std::size_t result{};
  const auto tail = fingerprint.substr(prefix.size());
  const auto parsed = std::from_chars(tail.data(), tail.data() + tail.size(),
                                      result);
  if (parsed.ec != std::errc{} || parsed.ptr != tail.data() + tail.size())
    throw runtime::Error("material flow global face count is invalid");
  return result;
}

void require_exact_cover(
    const flow::MaterialDensityFlowDiagnosticSource &source,
    const runtime::MpiContext &mpi) {
  bool preparation_failed = false;
  std::vector<int> cells;
  std::vector<int> faces;
  try {
    const auto extent = source.global_cell_extent();
    if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0)
      throw runtime::Error("material flow global cell extent is invalid");
    const std::uint64_t cell_count =
        static_cast<std::uint64_t>(extent.x) *
        static_cast<std::uint64_t>(extent.y) *
        static_cast<std::uint64_t>(extent.z);
    const auto face_count = global_face_count(source);
    if (cell_count > static_cast<std::uint64_t>(
                         std::numeric_limits<int>::max()) ||
        face_count > static_cast<std::size_t>(
                         std::numeric_limits<int>::max()))
      throw runtime::Error("material flow diagnostic exact cover is too large");
    cells.assign(static_cast<std::size_t>(cell_count), 0);
    faces.assign(face_count, 0);
    for (std::size_t item = 0; item < source.owned_cell_count(); ++item) {
      const auto global = source.field_global_id(3U, item);
      if (global >= cells.size() || ++cells[global] != 1)
        throw runtime::Error("material flow owned cell cover is invalid");
    }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (fault(kCellExactCoverInjectedFault))
      throw runtime::Error("injected material cell cover failure");
#endif
    for (std::size_t item = 0; item < source.canonical_owned_face_count();
         ++item) {
      const auto global = source.field_global_id(0U, item);
      if (global >= faces.size() || ++faces[global] != 1)
        throw runtime::Error("material flow owned face cover is invalid");
    }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (fault(kFaceExactCoverInjectedFault))
      throw runtime::Error("injected material face cover failure");
#endif
  } catch (...) {
    preparation_failed = true;
  }
  const int preparation_rank = lowest_rank(
      mpi, preparation_failed,
      "MPI_Allreduce(flow diagnostic cover preparation)");
  if (preparation_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.exact-cover", preparation_rank,
                     "material flow diagnostic exact cover is invalid");
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, cells.data(), static_cast<int>(cells.size()),
                    MPI_INT, MPI_SUM, mpi.comm()),
      "MPI_Allreduce(flow diagnostic cell cover)");
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, faces.data(), static_cast<int>(faces.size()),
                    MPI_INT, MPI_SUM, mpi.comm()),
      "MPI_Allreduce(flow diagnostic face cover)");
  const bool invalid =
      !std::all_of(cells.begin(), cells.end(), [](int value) {
        return value == 1;
      }) ||
      !std::all_of(faces.begin(), faces.end(), [](int value) {
        return value == 1;
      });
  const int invalid_rank = lowest_rank(
      mpi, invalid, "MPI_Allreduce(flow diagnostic cover validation)");
  if (invalid_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.exact-cover", invalid_rank,
                     "material flow diagnostic ownership cover is incomplete");
}

void collective_preflight(
    const flow::MaterialDensityFlowDiagnosticSource &source,
    const runtime::MpiContext &mpi, const DiagnosticRequest &request) {
  int local_kind = 0;
  try {
    require_request(source, request, DiagnosticScope::collective);
    if (source.relative_rank() != mpi.rank())
      local_kind = 1;
  } catch (const DiagnosticCollectionError &error) {
    local_kind = error.code() == "flow.diagnostics.stale-source"
                     ? 2
                     : (error.code() == "flow.diagnostics.unknown-field" ? 4
                                                                          : 1);
  } catch (...) {
    local_kind = 3;
  }
  const int rank = lowest_rank(mpi, local_kind != 0,
                               "MPI_Allreduce(flow diagnostic preflight)");
  if (rank >= 0) {
    int selected = mpi.rank() == rank ? local_kind : 0;
    runtime::check_mpi_result(
        MPI_Bcast(&selected, 1, MPI_INT, rank, mpi.comm()),
        "MPI_Bcast(flow diagnostic preflight)");
    if (selected == 2)
      collection_error(DiagnosticFailureClass::invalid_input,
                       "flow.diagnostics.stale-source", rank,
                       "material flow diagnostic source is stale");
    if (selected == 4)
      collection_error(DiagnosticFailureClass::invalid_request,
                       "flow.diagnostics.unknown-field", rank,
                       "unknown collective material flow diagnostic field");
    collection_error(selected == 3 ? DiagnosticFailureClass::layout
                                   : DiagnosticFailureClass::invalid_request,
                     selected == 3 ? "flow.diagnostics.layout"
                                   : "flow.diagnostics.frame",
                     rank, "collective material flow diagnostic preflight failed");
  }
  std::vector<unsigned char> request_key;
  std::vector<unsigned char> provider_key;
  bool preparation_failed = false;
  try {
    request_key = request_agreement_key(request);
    provider_key = provider_agreement_key(source);
  } catch (...) {
    preparation_failed = true;
  }
  const int preparation_rank = lowest_rank(
      mpi, preparation_failed,
      "MPI_Allreduce(flow diagnostic preflight preparation)");
  if (preparation_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.agreement-preparation",
                     preparation_rank,
                     "collective material flow diagnostic key preparation failed");
  require_byte_agreement(mpi, request_key,
                         DiagnosticFailureClass::invalid_request,
                         "flow.diagnostics.frame",
                         "collective material flow diagnostic requests");
  require_byte_agreement(mpi, provider_key, DiagnosticFailureClass::layout,
                         "flow.diagnostics.provider-agreement",
                         "collective material flow diagnostic providers");
  require_exact_cover(source, mpi);
}

Observation aggregate(const runtime::MpiContext &mpi, Observation local,
                      DiagnosticLevel level, std::size_t sample_budget) {
  std::array<double, 5> values{local.total_mass, local.total_momentum[0],
                               local.total_momentum[1],
                               local.total_momentum[2], local.density_min};
  mpi.allreduce_fp64_in_place(values.data(), 4U,
                              runtime::Fp64ReductionOperation::sum);
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, values.data() + 4U, 1, MPI_DOUBLE, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce(flow diagnostic density minimum)");
  double maximum = local.density_max;
  mpi.allreduce_fp64_in_place(&maximum, 1U,
                              runtime::Fp64ReductionOperation::maximum);
  local.total_mass = values[0];
  for (std::size_t c = 0; c < 3U; ++c)
    local.total_momentum[c] = values[c + 1U];
  local.density_min = values[4];
  local.density_max = maximum;
  std::uint64_t fp[2]{local.fingerprint.xor64, local.fingerprint.sum64};
  runtime::check_mpi_result(MPI_Allreduce(MPI_IN_PLACE, fp, 1, MPI_UINT64_T,
                                          MPI_BXOR, mpi.comm()),
                            "MPI_Allreduce(flow diagnostic fingerprint xor)");
  runtime::check_mpi_result(MPI_Allreduce(MPI_IN_PLACE, fp + 1, 1,
                                          MPI_UINT64_T, MPI_SUM, mpi.comm()),
                            "MPI_Allreduce(flow diagnostic fingerprint sum)");
  local.fingerprint = {fp[0], fp[1]};
  if (level == DiagnosticLevel::bounded_state_sample) {
    runtime::check_mpi_result(MPI_Allreduce(MPI_IN_PLACE, &local.eligible, 1,
                                            MPI_UINT64_T, MPI_SUM, mpi.comm()),
                              "MPI_Allreduce(flow diagnostic sample count)");
    // Local prefixes are sufficient for the global prefix once merged: an
    // item beyond a rank's local budget already has that many smaller local
    // items ahead of it.
    std::vector<unsigned char> send;
    std::vector<int> counts;
    bool send_preparation_failed = false;
    int send_bytes{};
    try {
      std::vector<SampleWireItem> items;
      items.reserve(local.samples.size());
      for (const auto &sample : local.samples) {
        const auto it = std::find(kSampleFieldIds.begin(), kSampleFieldIds.end(),
                                  sample.field_id);
        if (it == kSampleFieldIds.end())
          throw runtime::Error("material flow diagnostic sample field is invalid");
        items.push_back(
            {static_cast<std::uint32_t>(it - kSampleFieldIds.begin()),
             sample.global_id, sample.component, sample.value.bits,
             sample.value.status});
      }
      send = encode_sample_payload(items);
      if (send.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw runtime::Error(
            "material flow diagnostic sample wire is too large");
      send_bytes = static_cast<int>(send.size());
      counts.resize(static_cast<std::size_t>(mpi.size()));
    } catch (...) {
      send_preparation_failed = true;
    }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    send_preparation_failed =
        send_preparation_failed ||
        fault(kSampleSendPreparationInjectedFault);
#endif
    const int send_preparation_rank = lowest_rank(
        mpi, send_preparation_failed,
        "MPI_Allreduce(flow diagnostic sample preparation)");
    if (send_preparation_rank >= 0)
      collection_error(DiagnosticFailureClass::layout,
                       "flow.diagnostics.sample-preparation",
                       send_preparation_rank,
                       "material flow diagnostic sample preparation failed");
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (fault(kSampleWireMalformedInjectedFault))
      send.front() ^= 0x80U;
#endif
    runtime::check_mpi_result(MPI_Allgather(&send_bytes, 1, MPI_INT,
                                            counts.data(), 1, MPI_INT,
                                            mpi.comm()),
                              "MPI_Allgather(flow diagnostic sample sizes)");
    std::vector<int> offsets;
    std::vector<unsigned char> received;
    std::vector<SampleWireItem> merged;
    bool receive_preparation_failed = false;
    int total{};
    try {
      offsets.resize(counts.size());
      std::size_t maximum_items{};
      for (std::size_t rank = 0; rank < counts.size(); ++rank) {
        offsets[rank] = total;
        if (counts[rank] < 0 ||
            counts[rank] > std::numeric_limits<int>::max() - total)
          throw runtime::Error(
              "material flow diagnostic sample exchange wraps");
        const auto payload_size = static_cast<std::size_t>(counts[rank]);
        const std::size_t possible_items =
            payload_size < kSampleWireHeaderBytes
                ? 0U
                : std::min(kMaximumStateSamplesV1,
                           (payload_size - kSampleWireHeaderBytes) /
                               kSampleWireItemBytes);
        if (possible_items >
            std::numeric_limits<std::size_t>::max() - maximum_items)
          throw runtime::Error(
              "material flow diagnostic sample capacity wraps");
        maximum_items += possible_items;
        total += counts[rank];
      }
      received.resize(static_cast<std::size_t>(total));
      merged.reserve(maximum_items);
      local.samples.reserve(sample_budget);
    } catch (...) {
      receive_preparation_failed = true;
    }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    receive_preparation_failed = receive_preparation_failed ||
        fault(kSampleReceivePreparationInjectedFault);
#endif
    const int receive_preparation_rank = lowest_rank(
        mpi, receive_preparation_failed,
        "MPI_Allreduce(flow diagnostic sample receive preparation)");
    if (receive_preparation_rank >= 0)
      collection_error(
          DiagnosticFailureClass::layout,
          "flow.diagnostics.sample-preparation", receive_preparation_rank,
          "material flow diagnostic sample receive preparation failed");
    runtime::check_mpi_result(
        MPI_Allgatherv(send.data(), send_bytes, MPI_BYTE, received.data(),
                       counts.data(), offsets.data(), MPI_BYTE, mpi.comm()),
        "MPI_Allgatherv(flow diagnostic samples)");
    int wire_source = mpi.size();
    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto rank_index = static_cast<std::size_t>(rank);
      try {
        decode_sample_payload(
            received.data() + static_cast<std::size_t>(offsets[rank_index]),
            static_cast<std::size_t>(counts[rank_index]), merged);
      } catch (...) {
        wire_source = rank;
        break;
      }
    }
    int agreed_wire_source = mpi.size();
    runtime::check_mpi_result(
        MPI_Allreduce(&wire_source, &agreed_wire_source, 1, MPI_INT, MPI_MIN,
                      mpi.comm()),
        "MPI_Allreduce(flow diagnostic sample validation)");
    if (agreed_wire_source < mpi.size())
      collection_error(DiagnosticFailureClass::layout,
                       "flow.diagnostics.sample-wire", agreed_wire_source,
                       "material flow diagnostic sample wire is invalid");
    std::sort(merged.begin(), merged.end(),
              [](const SampleWireItem &a, const SampleWireItem &b) {
      return std::tie(a.field, a.global_id, a.component) <
             std::tie(b.field, b.global_id, b.component);
    });
    const bool duplicate = std::adjacent_find(
                               merged.begin(), merged.end(),
                               [](const SampleWireItem &a,
                                  const SampleWireItem &b) {
                                 return std::tie(a.field, a.global_id,
                                                 a.component) ==
                                        std::tie(b.field, b.global_id,
                                                 b.component);
                               }) != merged.end();
    const int duplicate_rank = lowest_rank(
        mpi, duplicate, "MPI_Allreduce(flow diagnostic sample uniqueness)");
    if (duplicate_rank >= 0)
      collection_error(DiagnosticFailureClass::layout,
                       "flow.diagnostics.sample-wire", duplicate_rank,
                       "material flow diagnostic samples are duplicated");
    local.samples.clear();
    const std::size_t retained = std::min(merged.size(), sample_budget);
    local.samples.reserve(retained);
    for (std::size_t index = 0; index < retained; ++index) {
      const auto &item = merged[index];
      local.samples.push_back(
          {std::string(kSampleFieldIds[item.field]), item.global_id,
           item.component, std::string(kSampleFieldUnits[item.field]),
           {item.status, item.value_bits}});
    }
  }
  return local;
}

} // namespace

DiagnosticDescriptor describe_diagnostics(
    const flow::MaterialDensityFlowDiagnosticSource &) noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::piso, kModuleId,
          kInstanceId, kCapabilities};
}

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::MaterialDensityFlowDiagnosticSource &source) {
  std::vector<std::string_view> result;
  result.reserve(source.fingerprint_field_count());
  for (std::size_t field = 0; field < source.fingerprint_field_count(); ++field)
    result.push_back(source.fingerprint_field_id(field));
  return result;
}

void collect_diagnostics(
    const flow::MaterialDensityFlowDiagnosticSource &source,
    const DiagnosticRequest &request, DiagnosticSink &sink) {
  require_request(source, request, DiagnosticScope::local);
  Observation observation;
  try {
    observation = observe(source, request);
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (const runtime::Error &error) {
    if (std::string_view(error.what()) ==
        "material flow diagnostic source is stale")
      collection_error(DiagnosticFailureClass::invalid_input,
                       "flow.diagnostics.stale-source", -1,
                       "material flow diagnostic source is stale");
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.layout", -1,
                     "material flow diagnostic provider failed");
  }
  auto record = build_record(source, request, observation,
                             DiagnosticScope::local);
  validate(record, describe_diagnostics(source), request);
  submit_local(std::move(record), sink);
}

void collect_diagnostics(
    const flow::MaterialDensityFlowDiagnosticSource &source,
    const runtime::MpiContext &mpi, const DiagnosticRequest &request,
    DiagnosticSink &sink) {
  collective_preflight(source, mpi, request);
  Observation local;
  bool provider_failed = false;
  try {
    local = observe(source, request);
  } catch (...) {
    provider_failed = true;
  }
  const int provider_rank = lowest_rank(
      mpi, provider_failed, "MPI_Allreduce(flow diagnostic provider)");
  if (provider_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.layout", provider_rank,
                     "collective material flow diagnostic provider failed");
  std::optional<Observation> global;
  bool aggregation_failed = false;
  try {
    global.emplace(
        aggregate(mpi, std::move(local), request.level, request.sample_budget));
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    aggregation_failed = true;
  }
  const int aggregation_rank = lowest_rank(
      mpi, aggregation_failed,
      "MPI_Allreduce(flow diagnostic aggregation completion)");
  if (aggregation_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.aggregation", aggregation_rank,
                     "collective material flow diagnostic aggregation failed");

  std::optional<DiagnosticRecord> prepared_record;
  bool record_failed = false;
  try {
    prepared_record.emplace(build_record(source, request, *global,
                                         DiagnosticScope::collective));
    validate(*prepared_record, describe_diagnostics(source), request);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (fault(kRecordValidationInjectedFault))
      throw runtime::Error("injected material diagnostic record failure");
#endif
  } catch (...) {
    record_failed = true;
  }
  const int record_rank = lowest_rank(
      mpi, record_failed, "MPI_Allreduce(flow diagnostic record validation)");
  if (record_rank >= 0)
    collection_error(DiagnosticFailureClass::layout,
                     "flow.diagnostics.record", record_rank,
                     "collective material flow diagnostic record is invalid");
  bool sink_failed = false;
  try {
    sink.submit(*prepared_record);
  } catch (...) {
    sink_failed = true;
  }
  const int sink_rank = lowest_rank(
      mpi, sink_failed, "MPI_Allreduce(flow diagnostic sink)");
  if (sink_rank >= 0)
    collection_error(DiagnosticFailureClass::sink_failure,
                     "diagnostics.sink.submit", sink_rank,
                     "collective diagnostic sink submission failed");
}

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
namespace detail {
void material_density_piso_set_fault_for_test(std::uint8_t value) noexcept {
  injected_fault = value;
}

void material_density_piso_reset_for_test() noexcept {
  injected_fault = kNoInjectedFault;
}

std::vector<unsigned char>
material_density_piso_encode_sample_wire_for_test(
    const std::vector<std::array<std::uint64_t, 5>> &items) {
  std::vector<SampleWireItem> encoded_items;
  encoded_items.reserve(items.size());
  for (const auto &item : items)
    encoded_items.push_back(
        {static_cast<std::uint32_t>(item[0]), item[1],
         static_cast<std::uint32_t>(item[2]), item[3],
         static_cast<DiagnosticValueStatus>(item[4])});
  return encode_sample_payload(encoded_items);
}

std::vector<std::array<std::uint64_t, 5>>
material_density_piso_decode_sample_wire_for_test(
    const std::vector<unsigned char> &bytes) {
  std::vector<SampleWireItem> decoded_items;
  decoded_items.reserve(kMaximumStateSamplesV1);
  decode_sample_payload(bytes.data(), bytes.size(), decoded_items);
  std::vector<std::array<std::uint64_t, 5>> result;
  result.reserve(decoded_items.size());
  for (const auto &item : decoded_items)
    result.push_back(
        {item.field, item.global_id, item.component, item.value_bits,
         static_cast<std::uint8_t>(item.status)});
  return result;
}
} // namespace detail
#endif

} // namespace hundun::diagnostics
