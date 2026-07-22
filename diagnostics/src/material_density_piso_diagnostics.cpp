// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/material_density_piso_diagnostics.hpp"

#include "hundun/flow/material_density_piso.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

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
  std::uint64_t selection_hash = UINT64_C(1469598103934665603);
  for (const auto field : request.selected_fields) {
    for (const char character : field) {
      const auto value = static_cast<unsigned char>(character);
      selection_hash ^= value;
      selection_hash *= UINT64_C(1099511628211);
    }
    selection_hash ^= 0xffU;
    selection_hash *= UINT64_C(1099511628211);
  }
  std::array<std::uint64_t, 6> key{
      static_cast<std::uint64_t>(request.level), request.frame.step,
      bits(request.frame.time_s), request.sample_budget,
      static_cast<std::uint64_t>(request.selected_fields.size()),
      selection_hash};
  auto root = key;
  runtime::check_mpi_result(
      MPI_Bcast(root.data(), static_cast<int>(root.size()), MPI_UINT64_T, 0,
                mpi.comm()),
      "MPI_Bcast(flow diagnostic request key)");
  const int mismatch = lowest_rank(
      mpi, key != root, "MPI_Allreduce(flow diagnostic request agreement)");
  if (mismatch >= 0)
    collection_error(DiagnosticFailureClass::invalid_request,
                     "flow.diagnostics.frame", mismatch,
                     "collective material flow diagnostic requests differ");
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
    struct Wire final {
      std::uint64_t global{};
      std::uint64_t value_bits{};
      std::uint32_t field{};
      std::uint32_t component{};
      std::uint32_t status{};
    };
    const std::array<std::string_view, 5> ids{
        "face_mass_flux", "face_velocity", "pi", "rho", "velocity"};
    std::vector<Wire> send;
    send.reserve(local.samples.size());
    for (const auto &sample : local.samples) {
      const auto it = std::find(ids.begin(), ids.end(), sample.field_id);
      send.push_back({sample.global_id, sample.value.bits,
                      static_cast<std::uint32_t>(it - ids.begin()),
                      sample.component,
                      static_cast<std::uint32_t>(sample.value.status)});
    }
    const std::size_t bytes_size = send.size() * sizeof(Wire);
    if (bytes_size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      collection_error(DiagnosticFailureClass::layout,
                       "flow.diagnostics.layout", -1,
                       "material flow diagnostic sample wire is too large");
    int send_bytes = static_cast<int>(bytes_size);
    std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
    runtime::check_mpi_result(MPI_Allgather(&send_bytes, 1, MPI_INT,
                                            counts.data(), 1, MPI_INT,
                                            mpi.comm()),
                              "MPI_Allgather(flow diagnostic sample sizes)");
    std::vector<int> offsets(counts.size());
    int total{};
    for (std::size_t rank = 0; rank < counts.size(); ++rank) {
      offsets[rank] = total;
      if (counts[rank] > std::numeric_limits<int>::max() - total)
        collection_error(DiagnosticFailureClass::layout,
                         "flow.diagnostics.layout", -1,
                         "material flow diagnostic sample exchange wraps");
      total += counts[rank];
    }
    std::vector<unsigned char> received(static_cast<std::size_t>(total));
    runtime::check_mpi_result(
        MPI_Allgatherv(send.data(), send_bytes, MPI_BYTE, received.data(),
                       counts.data(), offsets.data(), MPI_BYTE, mpi.comm()),
        "MPI_Allgatherv(flow diagnostic samples)");
    std::vector<Wire> merged(received.size() / sizeof(Wire));
    if (!received.empty())
      std::memcpy(merged.data(), received.data(), received.size());
    std::sort(merged.begin(), merged.end(), [](const Wire &a, const Wire &b) {
      return std::tie(a.field, a.global, a.component) <
             std::tie(b.field, b.global, b.component);
    });
    local.samples.clear();
    const std::size_t retained = std::min(merged.size(), sample_budget);
    local.samples.reserve(retained);
    for (std::size_t index = 0; index < retained; ++index) {
      const auto &item = merged[index];
      local.samples.push_back(
          {std::string(ids[item.field]), item.global, item.component,
           item.field == 0U   ? "kg/s"
           : item.field == 2U ? "Pa"
           : item.field == 3U ? "kg/m3"
                              : "m/s",
           {static_cast<DiagnosticValueStatus>(item.status), item.value_bits}});
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
  auto global =
      aggregate(mpi, std::move(local), request.level, request.sample_budget);
  auto record = build_record(source, request, global,
                             DiagnosticScope::collective);
  validate(record, describe_diagnostics(source), request);
  bool sink_failed = false;
  try {
    sink.submit(record);
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

} // namespace hundun::diagnostics
