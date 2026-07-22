// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"

#include "density_closure_detail.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#include "ideal_gas_closure_diagnostics_test_access.hpp"
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

namespace hundun::diagnostics {
namespace {

constexpr std::string_view kModuleId = "flow.ideal-gas-closure";
constexpr std::string_view kInstanceId = "primary";
constexpr std::string_view kPhase = "ideal-gas-closure.attempt-result";
constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(
        DiagnosticCapability::bounded_state_sample) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

constexpr std::array<std::string_view, 5> kSampleIds{"enthalpy", "p0", "rho",
                                                     "rho_h", "temperature"};

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
test::IdealGasClosureDiagnosticFault injected_fault{
    test::IdealGasClosureDiagnosticFault::none};
int injected_fault_rank{-1};
test::IdealGasClosureDiagnosticWork diagnostic_work{};

bool fault(test::IdealGasClosureDiagnosticFault value, int rank) noexcept {
  return injected_fault == value &&
         (injected_fault_rank < 0 || injected_fault_rank == rank);
}
#endif

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double relative_product_error(double a, double b, double c,
                              double expected) noexcept {
  const long double actual_value = static_cast<long double>(a) * b * c;
  const long double expected_value = expected;
  const long double denominator =
      std::max({std::abs(actual_value), std::abs(expected_value),
                std::numeric_limits<long double>::min()});
  const long double value =
      std::abs(actual_value - expected_value) / denominator;
  return std::isfinite(value) &&
                 value <= static_cast<long double>(
                              std::numeric_limits<double>::max())
             ? static_cast<double>(value)
             : std::numeric_limits<double>::infinity();
}

[[noreturn]] void fail(DiagnosticFailureClass classification, std::string code,
                       int rank, std::string message) {
  throw DiagnosticCollectionError(classification, std::move(code), rank,
                                  std::move(message));
}

void require_request(const flow::IdealGasClosureDiagnosticSource &source,
                     const DiagnosticRequest &request,
                     DiagnosticScope expected) {
  try {
    validate(request);
  } catch (const DiagnosticCollectionError &error) {
    if (error.code() ==
        std::string_view("diagnostics.request.selected-fields"))
      fail(DiagnosticFailureClass::invalid_request,
           "closure.diagnostics.selected-field", -1,
           "ideal-gas selected fields are not canonical");
    fail(DiagnosticFailureClass::invalid_request, "closure.diagnostics.frame",
         -1, "ideal-gas diagnostic request is invalid");
  } catch (...) {
    fail(DiagnosticFailureClass::invalid_request, "closure.diagnostics.frame",
         -1, "ideal-gas diagnostic request is invalid");
  }
  if (request.scope != expected ||
      request.frame.rank != source.relative_rank() ||
      request.frame.step != source.committed_step() ||
      bits(request.frame.time_s) != bits(source.committed_time_s()) ||
      request.frame.phase != kPhase)
    fail(DiagnosticFailureClass::invalid_request, "closure.diagnostics.frame",
         -1, "ideal-gas diagnostic frame does not match the source");
  if (request.sample_budget > kMaximumStateSamplesV1)
    fail(DiagnosticFailureClass::invalid_request, "closure.diagnostics.frame",
         -1, "ideal-gas diagnostic sample budget is invalid");
  if (!request.selected_fields.empty()) {
    if (!std::is_sorted(request.selected_fields.begin(),
                        request.selected_fields.end()) ||
        std::adjacent_find(request.selected_fields.begin(),
                           request.selected_fields.end()) !=
            request.selected_fields.end())
      fail(DiagnosticFailureClass::invalid_request,
           "closure.diagnostics.selected-field", -1,
           "ideal-gas selected fields are not canonical");
    for (const auto selected : request.selected_fields)
      if (!std::binary_search(kSampleIds.begin(), kSampleIds.end(), selected))
        fail(DiagnosticFailureClass::invalid_request,
             "closure.diagnostics.selected-field", -1,
             "ideal-gas selected field is unknown");
  }
}

void append_u64(std::vector<unsigned char> &bytes, std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U)
    bytes.push_back(static_cast<unsigned char>(value >> shift));
}

void append_text(std::vector<unsigned char> &bytes, std::string_view value) {
  append_u64(bytes, static_cast<std::uint64_t>(value.size()));
  bytes.insert(bytes.end(), value.begin(), value.end());
}

std::vector<unsigned char> request_key(const DiagnosticRequest &request) {
  std::vector<unsigned char> bytes;
  bytes.push_back(static_cast<unsigned char>(request.level));
  bytes.push_back(static_cast<unsigned char>(request.scope));
  append_u64(bytes, request.frame.step);
  append_u64(bytes, bits(request.frame.time_s));
  append_text(bytes, request.frame.phase);
  append_u64(bytes, static_cast<std::uint64_t>(request.selected_fields.size()));
  for (const auto field : request.selected_fields)
    append_text(bytes, field);
  append_u64(bytes, static_cast<std::uint64_t>(request.sample_budget));
  return bytes;
}

std::vector<unsigned char>
provider_key(const flow::IdealGasClosureDiagnosticSource &source) {
  std::vector<unsigned char> bytes;
  append_text(bytes, source.global_cell_layout_fingerprint());
  append_u64(bytes, source.committed_step());
  append_u64(bytes, bits(source.committed_time_s()));
  append_u64(
      bytes,
      bits(
          flow::detail::DensityClosureDiagnosticAccess::gas_constant_J_per_kg_K(
              source)));
  const auto closure = source.closure_state();
  bytes.push_back(static_cast<unsigned char>(closure.mode));
  append_u64(bytes, bits(closure.thermodynamic_pressure_pa));
  append_u64(bytes, closure.revision);
  bytes.push_back(closure.target_mass_kg.has_value() ? 1U : 0U);
  if (closure.target_mass_kg)
    append_u64(bytes, bits(*closure.target_mass_kg));
  const auto &outer = source.report();
  append_u64(bytes, outer.attempt_identity());
  const auto &parent = outer.flow();
  const auto &flow_report = parent.flow();
  bytes.push_back(static_cast<unsigned char>(flow_report.disposition));
  bytes.push_back(static_cast<unsigned char>(flow_report.reason));
  append_u64(bytes,
             static_cast<std::uint64_t>(flow_report.lowest_failing_rank + 1));
  append_u64(bytes, flow_report.pressure_corrector_count);
  append_u64(bytes, bits(flow_report.attempted_dt_s));
  append_u64(bytes, bits(flow_report.suggested_dt_s));
  const auto append_solve = [&](const linear::SolveReport &solve) {
    bytes.push_back(static_cast<unsigned char>(solve.reason));
    append_u64(bytes, solve.iterations);
    append_u64(bytes, bits(solve.initial_residual));
    append_u64(bytes, bits(solve.recursive_residual));
    append_u64(bytes, bits(solve.final_residual));
    append_u64(bytes, solve.matvec_count);
    append_u64(bytes, solve.preconditioner_apply_count);
    append_u64(bytes, solve.global_reduction_count);
    append_u64(bytes,
               static_cast<std::uint64_t>(solve.lowest_failing_rank + 1));
  };
  for (const auto &solve : flow_report.momentum.components)
    append_solve(solve);
  for (const auto &solve : flow_report.pressure)
    append_solve(solve);
  for (const double value :
       {flow_report.final_continuity_normalized_l2,
        flow_report.final_pressure_residual_l2,
        flow_report.final_mass_relative_conservation_defect})
    append_u64(bytes, bits(value));
  for (const double value : flow_report.final_momentum_normalized_l2)
    append_u64(bytes, bits(value));
  for (const double value : flow_report.final_transport_normalized_l2)
    append_u64(bytes, bits(value));
  for (const double value :
       flow_report.final_momentum_relative_conservation_defect)
    append_u64(bytes, bits(value));
  for (const double value :
       flow_report.final_transport_relative_conservation_defect)
    append_u64(bytes, bits(value));
  bytes.push_back(flow_report.final_backflow_evidence.has_value() ? 1U : 0U);
  append_u64(bytes, parent.material_field_count());
  append_u64(bytes, parent.shared_face_mass_flux_field());
  bytes.push_back(static_cast<unsigned char>(parent.flux_provenance()));
  bytes.push_back(static_cast<unsigned char>(parent.material_failure_reason()));
  bytes.push_back(parent.final_continuity_residual_available() ? 1U : 0U);
  bytes.push_back(parent.final_pressure_residual_available() ? 1U : 0U);
  append_u64(bytes, bits(parent.final_pressure_normalized_residual()));
  for (const auto value : parent.final_momentum_residual_availability())
    bytes.push_back(value);
  bytes.push_back(parent.mass_conservation_available() ? 1U : 0U);
  for (const auto value : parent.momentum_conservation_availability())
    bytes.push_back(value);
  bytes.push_back(parent.material_report_available() ? 1U : 0U);
  if (parent.material_report_available()) {
    const auto &material = parent.material_report();
    bytes.push_back(static_cast<unsigned char>(material.disposition()));
    bytes.push_back(static_cast<unsigned char>(material.reason()));
    append_u64(bytes,
               static_cast<std::uint64_t>(material.lowest_failing_rank() + 1));
    append_u64(bytes, material.attempt_identity());
    append_u64(bytes, material.finalization_identity());
    append_u64(bytes, material.shared_face_mass_flux_field());
    bytes.push_back(static_cast<unsigned char>(material.flux_provenance()));
    bytes.push_back(material.density_residual_available() ? 1U : 0U);
    append_u64(bytes, bits(material.density_normalized_l2()));
    bytes.push_back(material.mass_conservation_available() ? 1U : 0U);
    append_u64(bytes, bits(material.mass_relative_conservation_defect()));
    for (const auto value : material.transport_residual_availability())
      bytes.push_back(value);
    for (const double value : material.transport_normalized_l2())
      append_u64(bytes, bits(value));
    for (const auto value : material.transport_conservation_availability())
      bytes.push_back(value);
    for (const double value : material.transport_relative_conservation_defect())
      append_u64(bytes, bits(value));
    bytes.push_back(material.minimum_density_available() ? 1U : 0U);
    append_u64(bytes, bits(material.minimum_density_kg_per_m3()));
    append_u64(bytes, material.minimum_density_global_cell());
    append_u64(bytes,
               static_cast<std::uint64_t>(material.minimum_density_rank() + 1));
  }
  bytes.push_back(outer.closure_report_available() ? 1U : 0U);
  if (outer.closure_report_available()) {
    const auto &report = outer.closure_report();
    bytes.push_back(static_cast<unsigned char>(report.disposition()));
    bytes.push_back(static_cast<unsigned char>(report.reason()));
    bytes.push_back(static_cast<unsigned char>(report.stage()));
    append_u64(bytes,
               static_cast<std::uint64_t>(report.lowest_failing_rank() + 1));
    append_u64(bytes, report.evaluation_count());
    append_u64(bytes, report.collective_count());
    bytes.push_back(report.candidate_pressure_available() ? 1U : 0U);
    if (report.candidate_pressure_available())
      append_u64(bytes, bits(report.candidate_pressure_pa()));
    bytes.push_back(report.final_metrics_available() ? 1U : 0U);
    if (report.final_metrics_available()) {
      for (const double value :
           {report.actual_mass_kg(), report.temperature_min_K(),
            report.temperature_max_K(), report.enthalpy_min_J_per_kg(),
            report.enthalpy_max_J_per_kg(), report.density_min_kg_per_m3(),
            report.density_max_kg_per_m3(), report.rho_remap_normalized_l2(),
            report.rho_h_remap_normalized_l2(),
            report.rho_remap_relative_conservation_defect(),
            report.rho_h_remap_relative_conservation_defect(),
            report.enthalpy_temperature_max_relative_error(),
            report.eos_max_relative_error()})
        append_u64(bytes, bits(value));
    }
  }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(test::IdealGasClosureDiagnosticFault::provider_agreement,
            source.relative_rank()))
    bytes.push_back(0xffU);
#endif
  return bytes;
}

int lowest_rank(const runtime::MpiContext &mpi, bool failure,
                std::string_view operation) {
  const int local = failure ? mpi.rank() : mpi.size();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&local, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      operation);
  return lowest == mpi.size() ? -1 : lowest;
}

void require_byte_agreement(const runtime::MpiContext &mpi,
                            const std::vector<unsigned char> &local,
                            DiagnosticFailureClass classification,
                            const char *code, const char *message) {
  std::uint64_t size =
      mpi.rank() == 0 ? static_cast<std::uint64_t>(local.size()) : 0U;
  runtime::check_mpi_result(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
                            "MPI_Bcast(ideal-gas diagnostic agreement size)");
  if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.aggregation",
         mpi.rank(), "ideal-gas diagnostic agreement is too large");
  std::vector<unsigned char> reference(static_cast<std::size_t>(size));
  if (mpi.rank() == 0)
    reference = local;
  runtime::check_mpi_result(MPI_Bcast(reference.data(),
                                      static_cast<int>(reference.size()),
                                      MPI_BYTE, 0, mpi.comm()),
                            "MPI_Bcast(ideal-gas diagnostic agreement bytes)");
  const int mismatch =
      lowest_rank(mpi, local != reference,
                  "MPI_Allreduce(ideal-gas diagnostic agreement rank)");
  if (mismatch >= 0)
    fail(classification, code, mismatch, message);
}

void require_source_communicator(
    const flow::IdealGasClosureDiagnosticSource &source,
    const runtime::MpiContext &mpi) {
  const auto &bound = flow::detail::DensityClosureDiagnosticAccess::mpi(source);
  int comparison = MPI_UNEQUAL;
  runtime::check_mpi_result(
      MPI_Comm_compare(bound.comm(), mpi.comm(), &comparison),
      "MPI_Comm_compare(ideal-gas diagnostics)");
  if (comparison != MPI_IDENT || bound.rank() != mpi.rank() ||
      bound.size() != mpi.size())
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.communicator", mpi.rank(),
         "ideal-gas diagnostic communicator does not match the source");
}

void require_exact_cell_cover(
    const flow::IdealGasClosureDiagnosticSource &source,
    const runtime::MpiContext &mpi) {
  const auto session =
      flow::detail::DensityClosureDiagnosticAccess::acquire_committed(source);
  struct Wire final {
    std::int64_t begin[3];
    std::int64_t end[3];
    std::int64_t global[3];
    std::uint64_t count;
  };
  const auto box = session.owned_box;
  const auto global = session.global_extent;
  Wire local{{box.begin.x, box.begin.y, box.begin.z},
             {box.end.x, box.end.y, box.end.z},
             {global.x, global.y, global.z},
             source.owned_cell_count()};
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(test::IdealGasClosureDiagnosticFault::ownership_layout,
            mpi.rank()))
    ++local.count;
#endif
  std::vector<Wire> all(static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(
      MPI_Allgather(&local, static_cast<int>(sizeof(Wire)), MPI_BYTE,
                    all.data(), static_cast<int>(sizeof(Wire)), MPI_BYTE,
                    mpi.comm()),
      "MPI_Allgather(ideal-gas diagnostic ownership)");
  std::uint64_t covered{};
  for (std::size_t left = 0; left < all.size(); ++left) {
    std::uint64_t volume = 1U;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (all[left].global[axis] != local.global[axis] ||
          all[left].global[axis] <= 0 || all[left].begin[axis] < 0 ||
          all[left].end[axis] <= all[left].begin[axis] ||
          all[left].end[axis] > all[left].global[axis])
        fail(DiagnosticFailureClass::layout, "closure.diagnostics.ownership",
             static_cast<int>(left),
             "ideal-gas diagnostic ownership is invalid");
      const auto width = static_cast<std::uint64_t>(all[left].end[axis] -
                                                    all[left].begin[axis]);
      if (volume > std::numeric_limits<std::uint64_t>::max() / width)
        fail(DiagnosticFailureClass::layout, "closure.diagnostics.ownership",
             static_cast<int>(left),
             "ideal-gas diagnostic ownership overflows");
      volume *= width;
    }
    if (volume != all[left].count ||
        covered > std::numeric_limits<std::uint64_t>::max() - volume)
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.ownership",
           static_cast<int>(left),
           "ideal-gas diagnostic ownership count is invalid");
    covered += volume;
    for (std::size_t right = 0; right < left; ++right) {
      bool overlap = true;
      for (std::size_t axis = 0; axis < 3U; ++axis)
        overlap = overlap && all[left].begin[axis] < all[right].end[axis] &&
                  all[right].begin[axis] < all[left].end[axis];
      if (overlap)
        fail(DiagnosticFailureClass::layout, "closure.diagnostics.ownership",
             static_cast<int>(left),
             "ideal-gas diagnostic ownership overlaps");
    }
  }
  std::uint64_t global_count = 1U;
  for (const auto extent : local.global) {
    const auto width = static_cast<std::uint64_t>(extent);
    if (global_count > std::numeric_limits<std::uint64_t>::max() / width)
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.ownership",
           mpi.rank(), "ideal-gas diagnostic global count overflows");
    global_count *= width;
  }
  if (covered != global_count)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.ownership",
         mpi.rank(), "ideal-gas diagnostic ownership does not exactly cover");
}

struct CommittedObservation final {
  double mass{};
  double energy{};
  double rho_min{std::numeric_limits<double>::infinity()};
  double rho_max{-std::numeric_limits<double>::infinity()};
  double h_min{std::numeric_limits<double>::infinity()};
  double h_max{-std::numeric_limits<double>::infinity()};
  double t_min{std::numeric_limits<double>::infinity()};
  double t_max{-std::numeric_limits<double>::infinity()};
  double eos_max{};
  DiagnosticFingerprintParts fingerprint{};
  std::vector<double> rho;
  std::vector<double> rho_h;
};

CommittedObservation
observe(const flow::IdealGasClosureDiagnosticSource &source,
        DiagnosticLevel level) {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  ++diagnostic_work.observations;
  if (fault(test::IdealGasClosureDiagnosticFault::provider_agreement,
            source.relative_rank()))
    throw runtime::Error("injected ideal-gas diagnostic provider failure");
#endif
  CommittedObservation result;
  DiagnosticFingerprintAccumulator fingerprint;
  const auto closure = source.closure_state();
  const double p0 = closure.thermodynamic_pressure_pa;
  const double gas_constant =
      flow::detail::DensityClosureDiagnosticAccess::gas_constant_J_per_kg_K(
          source);
  const double cp =
      flow::detail::DensityClosureDiagnosticAccess::cp_J_per_kg_K(source);
  const auto session =
      flow::detail::DensityClosureDiagnosticAccess::acquire_committed(source);
  const auto extent = session.density.interior_extent();
  const std::size_t count = source.owned_cell_count();
  if (source.relative_rank() == 0) {
    fingerprint.add("p0", 0U, 0U, describe_fp64(p0));
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    ++diagnostic_work.fingerprint_items;
#endif
  }
  if (level == DiagnosticLevel::bounded_state_sample) {
    result.rho.resize(count);
    result.rho_h.resize(count);
  }
  for (std::size_t cell = 0; cell < count; ++cell) {
    const auto plane =
        static_cast<std::size_t>(extent.x) * static_cast<std::size_t>(extent.y);
    const int k = static_cast<int>(cell / plane);
    const auto within = cell % plane;
    const int j = static_cast<int>(within / static_cast<std::size_t>(extent.x));
    const int i = static_cast<int>(within % static_cast<std::size_t>(extent.x));
    const double rho = session.density(i, j, k, 0);
    const double rho_h = session.enthalpy_density(i, j, k, 0);
    const auto global_id = source.fingerprint_field_global_id(1U, cell);
    fingerprint.add("rho", global_id, 0U, describe_fp64(rho));
    fingerprint.add("rho_h", global_id, 0U, describe_fp64(rho_h));
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    diagnostic_work.fingerprint_items += 2U;
#endif
    if (level == DiagnosticLevel::bounded_state_sample) {
      result.rho[cell] = rho;
      result.rho_h[cell] = rho_h;
    }
    if (level != DiagnosticLevel::summary &&
        level != DiagnosticLevel::invariants)
      continue;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    ++diagnostic_work.summary_items;
#endif
    const double h = rho_h / rho;
    const double temperature = h / cp;
    const double volume = source.cell_volume_m3(cell);
    if (level == DiagnosticLevel::summary) {
      result.mass += volume * rho;
      result.energy += volume * rho_h;
    }
    result.rho_min = std::min(result.rho_min, rho);
    result.h_min = std::min(result.h_min, h);
    result.t_min = std::min(result.t_min, temperature);
    if (level == DiagnosticLevel::summary) {
      result.rho_max = std::max(result.rho_max, rho);
      result.h_max = std::max(result.h_max, h);
      result.t_max = std::max(result.t_max, temperature);
    }
    result.eos_max =
        std::max(result.eos_max,
                 relative_product_error(rho, gas_constant, temperature, p0));
  }
  result.fingerprint = fingerprint.parts();
  return result;
}

DiagnosticFailure mapped_failure(const flow::IdealGasStepAttemptReport &outer,
                                 DiagnosticStatus &status, int rank) {
  if (!outer.closure_report_available()) {
    status = DiagnosticStatus::unavailable;
    return {DiagnosticFailureClass::unavailable, "closure.not-evaluated", rank};
  }
  const auto &report = outer.closure_report();
  if (report.disposition() == flow::IdealGasClosureDisposition::closed) {
    status = outer.flow().flow().disposition ==
                     flow::StepAttemptDisposition::committed
                 ? DiagnosticStatus::ok
                 : DiagnosticStatus::warning;
    return {};
  }
  status = DiagnosticStatus::failed;
  using R = flow::IdealGasClosureFailureReason;
  switch (report.reason()) {
  case R::none:
    return {};
  case R::invalid_input:
    return {DiagnosticFailureClass::invalid_input, "closure.invalid-input",
            rank};
  case R::non_finite_enthalpy:
    return {DiagnosticFailureClass::non_finite_state,
            "closure.enthalpy-non-finite", rank};
  case R::non_positive_enthalpy:
    return {DiagnosticFailureClass::non_positive_state,
            "closure.enthalpy-non-positive", rank};
  case R::non_finite_temperature:
    return {DiagnosticFailureClass::non_finite_state,
            "closure.temperature-non-finite", rank};
  case R::non_positive_temperature:
    return {DiagnosticFailureClass::non_positive_state,
            "closure.temperature-non-positive", rank};
  case R::denominator_breakdown:
    return {DiagnosticFailureClass::numerical_breakdown,
            "closure.denominator-breakdown", rank};
  case R::non_finite_pressure:
    return {DiagnosticFailureClass::non_finite_state,
            "closure.pressure-non-finite", rank};
  case R::non_positive_pressure:
    return {DiagnosticFailureClass::non_positive_state,
            "closure.pressure-non-positive", rank};
  case R::non_finite_density:
    return {DiagnosticFailureClass::non_finite_state,
            "closure.density-non-finite", rank};
  case R::non_positive_density:
    return {DiagnosticFailureClass::non_positive_state,
            "closure.density-non-positive", rank};
  case R::eos_residual:
    return {DiagnosticFailureClass::non_convergence, "closure.eos-residual",
            rank};
  case R::remap_residual:
    return {DiagnosticFailureClass::non_convergence, "closure.remap-residual",
            rank};
  case R::mass_conservation:
    return {DiagnosticFailureClass::conservation, "closure.mass-conservation",
            rank};
  case R::enthalpy_conservation:
    return {DiagnosticFailureClass::conservation,
            "closure.enthalpy-conservation", rank};
  case R::collective_operation:
    return {DiagnosticFailureClass::collective_operation,
            "closure.collective-operation", rank};
  }
  return {DiagnosticFailureClass::invalid_input, "closure.invalid-input", rank};
}

DiagnosticMetric metric(std::string id, DiagnosticMetricKind kind,
                        std::string unit, double value) {
  return {std::move(id), kind, std::move(unit), describe_fp64(value)};
}

DiagnosticMetric unavailable_metric(std::string id, DiagnosticMetricKind kind,
                                    std::string unit) {
  return {std::move(id), kind, std::move(unit), {}};
}

DiagnosticInvariant make_invariant(std::string id, std::string unit,
                                   double observed, double limit,
                                   InvariantRelation relation) {
  DiagnosticInvariant item{std::move(id),
                           std::move(unit),
                           describe_fp64(observed),
                           relation == InvariantRelation::positive ||
                                   relation == InvariantRelation::finite
                               ? DiagnosticFp64{}
                               : describe_fp64(limit),
                           relation,
                           false};
  item.passed = evaluate_invariant(item);
  return item;
}

DiagnosticRecord
build_record(const flow::IdealGasClosureDiagnosticSource &source,
             const DiagnosticRequest &request,
             const CommittedObservation &observation, DiagnosticScope scope) {
  DiagnosticRecord record;
  const auto descriptor = describe_diagnostics(source);
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
  const auto &outer = source.report();
  int failure_rank = -1;
  if (scope == DiagnosticScope::collective) {
    if (outer.closure_report_available() &&
        outer.closure_report().disposition() !=
            flow::IdealGasClosureDisposition::closed)
      failure_rank = outer.closure_report().lowest_failing_rank();
    else if (!outer.closure_report_available())
      failure_rank = outer.flow().flow().lowest_failing_rank;
  }
  record.failure = mapped_failure(outer, record.status, failure_rank);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(test::IdealGasClosureDiagnosticFault::status_warning,
            source.relative_rank())) {
    record.status = DiagnosticStatus::warning;
    record.failure = {};
  } else if (fault(test::IdealGasClosureDiagnosticFault::status_failed,
                   source.relative_rank())) {
    record.status = DiagnosticStatus::failed;
    record.failure = {DiagnosticFailureClass::non_convergence,
                      "closure.injected-failure",
                      scope == DiagnosticScope::collective ? 0 : -1};
  } else if (fault(test::IdealGasClosureDiagnosticFault::status_unavailable,
                   source.relative_rank())) {
    record.status = DiagnosticStatus::unavailable;
    record.failure = {DiagnosticFailureClass::unavailable,
                      "closure.injected-unavailable",
                      scope == DiagnosticScope::collective ? 0 : -1};
  }
#endif
  const auto closure = source.closure_state();
  const auto layout = scope == DiagnosticScope::collective
                          ? source.global_cell_layout_fingerprint()
                          : source.owned_cell_layout_fingerprint();
  record.identities = {{"closure.attempt", std::nullopt,
                        outer.attempt_identity(), std::nullopt, std::nullopt},
                       {"closure.state", std::nullopt, closure.revision,
                        std::nullopt, std::nullopt},
                       {"field.p0",
                        std::string("replicated-scalar.rank0-owner.v1"),
                        std::nullopt, std::nullopt, std::nullopt},
                       {"field.rho", std::string(layout), std::nullopt,
                        std::nullopt, std::nullopt},
                       {"field.rho_h", std::string(layout), std::nullopt,
                        std::nullopt, std::nullopt},
                       {"layout.cells", std::string(layout), std::nullopt,
                        std::nullopt, std::nullopt}};
  DiagnosticFingerprintAccumulator fingerprint;
  fingerprint.combine(observation.fingerprint);
  record.state_fingerprint = fingerprint.finish();

  const flow::IdealGasClosureReport *report =
      outer.closure_report_available() ? &outer.closure_report() : nullptr;
  if (request.level == DiagnosticLevel::summary) {
    const auto final_value = [&](std::string id, DiagnosticMetricKind kind,
                                 std::string unit, auto getter) {
      return report != nullptr && report->final_metrics_available()
                 ? metric(std::move(id), kind, std::move(unit), getter(*report))
                 : unavailable_metric(std::move(id), kind, std::move(unit));
    };
    record.metrics = {
        final_value("closure.actual-mass", DiagnosticMetricKind::conservation,
                    "kg", [](const auto &r) { return r.actual_mass_kg(); }),
        report != nullptr && report->candidate_pressure_available()
            ? metric("closure.candidate-pressure",
                     DiagnosticMetricKind::state_summary, "Pa",
                     report->candidate_pressure_pa())
            : unavailable_metric("closure.candidate-pressure",
                                 DiagnosticMetricKind::state_summary, "Pa"),
        metric("closure.committed-mass", DiagnosticMetricKind::conservation,
               "kg", observation.mass),
        metric("closure.committed-pressure",
               DiagnosticMetricKind::state_summary, "Pa",
               closure.thermodynamic_pressure_pa),
        metric("closure.density-maximum", DiagnosticMetricKind::state_summary,
               "kg/m3", observation.rho_max),
        metric("closure.density-minimum", DiagnosticMetricKind::state_summary,
               "kg/m3", observation.rho_min),
        metric("closure.enthalpy-maximum", DiagnosticMetricKind::state_summary,
               "J/kg", observation.h_max),
        metric("closure.enthalpy-minimum", DiagnosticMetricKind::state_summary,
               "J/kg", observation.h_min),
        final_value("closure.enthalpy-temperature-error",
                    DiagnosticMetricKind::residual, "1",
                    [](const auto &r) {
                      return r.enthalpy_temperature_max_relative_error();
                    }),
        final_value("closure.eos-error", DiagnosticMetricKind::residual, "1",
                    [](const auto &r) { return r.eos_max_relative_error(); }),
        final_value("closure.rho-h-remap-conservation",
                    DiagnosticMetricKind::conservation, "1",
                    [](const auto &r) {
                      return r.rho_h_remap_relative_conservation_defect();
                    }),
        final_value(
            "closure.rho-h-remap-l2", DiagnosticMetricKind::residual, "1",
            [](const auto &r) { return r.rho_h_remap_normalized_l2(); }),
        final_value("closure.rho-remap-conservation",
                    DiagnosticMetricKind::conservation, "1",
                    [](const auto &r) {
                      return r.rho_remap_relative_conservation_defect();
                    }),
        final_value("closure.rho-remap-l2", DiagnosticMetricKind::residual, "1",
                    [](const auto &r) { return r.rho_remap_normalized_l2(); }),
        closure.target_mass_kg
            ? metric("closure.target-mass", DiagnosticMetricKind::state_summary,
                     "kg", *closure.target_mass_kg)
            : unavailable_metric("closure.target-mass",
                                 DiagnosticMetricKind::state_summary, "kg"),
        metric("closure.temperature-maximum",
               DiagnosticMetricKind::state_summary, "K", observation.t_max),
        metric("closure.temperature-minimum",
               DiagnosticMetricKind::state_summary, "K", observation.t_min)};
  } else if (request.level == DiagnosticLevel::invariants) {
    auto add = [&](std::string id, std::string unit, double value, double limit,
                   InvariantRelation relation) {
      record.invariants.push_back(make_invariant(std::move(id), std::move(unit),
                                                 value, limit, relation));
    };
    add("closure.committed-eos", "1", observation.eos_max, 1.0e-12,
        InvariantRelation::less_equal);
    add("closure.final-evidence-available", "1",
        report != nullptr && report->final_metrics_available() ? 1.0 : 0.0, 1.0,
        InvariantRelation::equal);
    add("closure.pressure-positive", "Pa", closure.thermodynamic_pressure_pa,
        0.0, InvariantRelation::positive);
    add("density.positive", "kg/m3", observation.rho_min, 0.0,
        InvariantRelation::positive);
    add("enthalpy.positive", "J/kg", observation.h_min, 0.0,
        InvariantRelation::positive);
    add("temperature.positive", "K", observation.t_min, 0.0,
        InvariantRelation::positive);
  } else if (request.level == DiagnosticLevel::counters) {
    record.counters = {{"closure.collectives", "count",
                        report == nullptr ? 0U : report->collective_count()},
                       {"closure.evaluations", "count",
                        report == nullptr ? 0U : report->evaluation_count()},
                       {"closure.revision", "count", closure.revision}};
  } else {
    const auto &selected = request.selected_fields.empty()
                               ? std::vector<std::string_view>(
                                     kSampleIds.begin(), kSampleIds.end())
                               : request.selected_fields;
    record.sample_budget = request.sample_budget;
    for (const auto id : selected) {
      const auto field = static_cast<std::size_t>(
          std::find(kSampleIds.begin(), kSampleIds.end(), id) -
          kSampleIds.begin());
      const auto eligible = source.sample_field_item_count(field);
      if (eligible > std::numeric_limits<std::uint64_t>::max() -
                         record.eligible_sample_count)
        fail(DiagnosticFailureClass::layout,
             "closure.diagnostics.sample-preparation", source.relative_rank(),
             "ideal-gas eligible sample count overflows");
      record.eligible_sample_count += eligible;
      for (std::size_t item = 0; item < source.sample_field_item_count(field) &&
                                 record.samples.size() < request.sample_budget;
           ++item) {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
        ++diagnostic_work.sample_items;
#endif
        double value = closure.thermodynamic_pressure_pa;
        if (field != 1U) {
          const double rho = observation.rho[item];
          const double rho_h = observation.rho_h[item];
          const double enthalpy = rho_h / rho;
          value = field == 0U   ? enthalpy
                  : field == 2U ? rho
                  : field == 3U
                      ? rho_h
                      : enthalpy /
                            flow::detail::DensityClosureDiagnosticAccess::
                                cp_J_per_kg_K(source);
        }
        record.samples.push_back(
            {std::string(id), source.sample_field_global_id(field, item), 0U,
             std::string(source.sample_field_unit(field)),
             describe_fp64(value)});
      }
    }
    std::sort(record.samples.begin(), record.samples.end(),
              [](const auto &a, const auto &b) {
                return std::tie(a.field_id, a.global_id, a.component) <
                       std::tie(b.field_id, b.global_id, b.component);
              });
    if (record.samples.size() > request.sample_budget)
      record.samples.resize(request.sample_budget);
    record.samples_truncated =
        record.eligible_sample_count > record.samples.size();
  }
  return record;
}

void submit(DiagnosticRecord record,
            const flow::IdealGasClosureDiagnosticSource &source,
            const DiagnosticRequest &request, DiagnosticSink &sink) {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(test::IdealGasClosureDiagnosticFault::record_validation,
            source.relative_rank()))
    record.schema_version = 0U;
#endif
  try {
    validate(record, describe_diagnostics(source), request);
  } catch (...) {
    fail(DiagnosticFailureClass::invalid_input, "closure.diagnostics.record",
         -1, "ideal-gas diagnostic record is invalid");
  }
  try {
    sink.submit(record);
  } catch (...) {
    fail(DiagnosticFailureClass::sink_failure, "diagnostics.sink.submit", -1,
         "diagnostic sink submission failed");
  }
}

void aggregate_samples(DiagnosticRecord &record,
                       const runtime::MpiContext &mpi) {
  if (record.level != DiagnosticLevel::bounded_state_sample)
    return;
  struct Wire final {
    std::uint64_t global_id{};
    std::uint64_t value_bits{};
    std::uint32_t field{};
    std::uint8_t status{};
  };
  std::vector<Wire> local;
  local.reserve(record.samples.size());
  for (const auto &sample : record.samples) {
    const auto found =
        std::find(kSampleIds.begin(), kSampleIds.end(), sample.field_id);
    if (found == kSampleIds.end())
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
           mpi.rank(), "ideal-gas diagnostic sample field is invalid");
    local.push_back({sample.global_id, sample.value.bits,
                     static_cast<std::uint32_t>(found - kSampleIds.begin()),
                     static_cast<std::uint8_t>(sample.value.status)});
  }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (!local.empty() &&
      fault(test::IdealGasClosureDiagnosticFault::sample_wire, mpi.rank()))
    local.push_back(local.front());
#endif
  if (local.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(Wire))
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", mpi.rank(),
         "ideal-gas diagnostic sample payload is too large");
  const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
  std::vector<int> counts(static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(MPI_Allgather(&local_bytes, 1, MPI_INT,
                                          counts.data(), 1, MPI_INT,
                                          mpi.comm()),
                            "MPI_Allgather(ideal-gas diagnostic sample sizes)");
  std::vector<int> displacements(counts.size());
  int total{};
  for (std::size_t rank = 0; rank < counts.size(); ++rank) {
    if (counts[rank] < 0 ||
        counts[rank] % static_cast<int>(sizeof(Wire)) != 0 ||
        counts[rank] > std::numeric_limits<int>::max() - total)
      fail(DiagnosticFailureClass::layout,
           "closure.diagnostics.sample-preparation", mpi.rank(),
           "ideal-gas diagnostic sample size is invalid");
    displacements[rank] = total;
    total += counts[rank];
  }
  std::vector<unsigned char> merged_bytes(static_cast<std::size_t>(total));
  runtime::check_mpi_result(
      MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, merged_bytes.data(),
                     counts.data(), displacements.data(), MPI_BYTE, mpi.comm()),
      "MPI_Allgatherv(ideal-gas diagnostic samples)");
  std::vector<Wire> merged(merged_bytes.size() / sizeof(Wire));
  if (!merged.empty())
    std::memcpy(merged.data(), merged_bytes.data(), merged_bytes.size());
  std::sort(merged.begin(), merged.end(),
            [](const Wire &left, const Wire &right) {
              return std::tie(left.field, left.global_id) <
                     std::tie(right.field, right.global_id);
            });
  if (std::adjacent_find(merged.begin(), merged.end(),
                         [](const Wire &left, const Wire &right) {
                           return left.field == right.field &&
                                  left.global_id == right.global_id;
                         }) != merged.end())
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
         injected_fault == test::IdealGasClosureDiagnosticFault::sample_wire
             ? (injected_fault_rank < 0 ? 0 : injected_fault_rank)
             : mpi.rank(),
#else
         mpi.rank(),
#endif
         "ideal-gas diagnostic samples are duplicated");
  std::vector<std::uint64_t> eligible_by_rank(
      static_cast<std::size_t>(mpi.size()));
  runtime::check_mpi_result(
      MPI_Allgather(&record.eligible_sample_count, 1, MPI_UINT64_T,
                    eligible_by_rank.data(), 1, MPI_UINT64_T, mpi.comm()),
      "MPI_Allgather(ideal-gas diagnostic eligible samples)");
  std::uint64_t eligible{};
  for (const auto value : eligible_by_rank) {
    if (value > std::numeric_limits<std::uint64_t>::max() - eligible)
      fail(DiagnosticFailureClass::layout,
           "closure.diagnostics.sample-preparation", mpi.rank(),
           "ideal-gas eligible sample count overflows");
    eligible += value;
  }
  record.eligible_sample_count = eligible;
  if (merged.size() > record.sample_budget)
    merged.resize(record.sample_budget);
  record.samples.clear();
  record.samples.reserve(merged.size());
  for (const auto &item : merged) {
    if (item.field >= kSampleIds.size() ||
        item.status >
            static_cast<std::uint8_t>(DiagnosticValueStatus::unavailable))
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
           mpi.rank(), "ideal-gas diagnostic sample wire is invalid");
    const std::size_t field = item.field;
    record.samples.push_back(
        {std::string(kSampleIds[field]),
         item.global_id,
         0U,
         std::string(field == 0U   ? "J/kg"
                     : field == 1U ? "Pa"
                     : field == 2U ? "kg/m3"
                     : field == 3U ? "J/m3"
                                   : "K"),
         {static_cast<DiagnosticValueStatus>(item.status), item.value_bits}});
  }
  record.samples_truncated =
      record.eligible_sample_count > record.samples.size();
}

void submit_collective(const DiagnosticRecord &record,
                       const flow::IdealGasClosureDiagnosticSource &source,
                       const DiagnosticRequest &request,
                       const runtime::MpiContext &mpi, DiagnosticSink &sink) {
  bool invalid = false;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  invalid = fault(test::IdealGasClosureDiagnosticFault::record_validation,
                  mpi.rank());
#endif
  try {
    if (!invalid)
      validate(record, describe_diagnostics(source), request);
  } catch (...) {
    invalid = true;
  }
  int rank = lowest_rank(
      mpi, invalid, "MPI_Allreduce(ideal-gas diagnostic record validation)");
  if (rank >= 0)
    fail(DiagnosticFailureClass::invalid_input, "closure.diagnostics.record",
         rank, "ideal-gas diagnostic record is invalid");
  bool sink_failed = false;
  try {
    sink.submit(record);
  } catch (...) {
    sink_failed = true;
  }
  rank = lowest_rank(mpi, sink_failed,
                     "MPI_Allreduce(ideal-gas diagnostic sink submission)");
  if (rank >= 0)
    fail(DiagnosticFailureClass::sink_failure, "diagnostics.sink.submit", rank,
         "collective diagnostic sink submission failed");
}

} // namespace

DiagnosticDescriptor describe_ideal_gas_closure_diagnostics() noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::density_closure,
          kModuleId, kInstanceId, kCapabilities};
}

DiagnosticDescriptor
describe_diagnostics(const flow::IdealGasClosureDiagnosticSource &) noexcept {
  return describe_ideal_gas_closure_diagnostics();
}

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::IdealGasClosureDiagnosticSource &source) {
  std::vector<std::string_view> result;
  result.reserve(source.fingerprint_field_count());
  for (std::size_t field = 0; field < source.fingerprint_field_count(); ++field)
    result.push_back(source.fingerprint_field_id(field));
  return result;
}

void collect_diagnostics(const flow::IdealGasClosureDiagnosticSource &source,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  try {
    require_request(source, request, DiagnosticScope::local);
    submit(build_record(source, request, observe(source, request.level),
                        DiagnosticScope::local),
           source, request, sink);
  } catch (const DiagnosticCollectionError &) {
    throw;
  } catch (const runtime::Error &) {
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.stale-source", -1,
         "ideal-gas diagnostic source is stale");
  }
}

void collect_diagnostics(const flow::IdealGasClosureDiagnosticSource &source,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) try {
  require_source_communicator(source, mpi);
  bool stale = false;
  bool request_invalid = false;
  bool selected_invalid = false;
  try {
    require_request(source, request, DiagnosticScope::collective);
  } catch (const DiagnosticCollectionError &error) {
    selected_invalid =
        error.code() == std::string_view("closure.diagnostics.selected-field");
    request_invalid = !selected_invalid;
  } catch (...) {
    stale = true;
  }
  int failed_rank = lowest_rank(
      mpi, stale, "MPI_Allreduce(ideal-gas diagnostic stale source)");
  if (failed_rank >= 0)
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.stale-source", failed_rank,
         "ideal-gas diagnostic source is stale");
  failed_rank = lowest_rank(
      mpi, selected_invalid,
      "MPI_Allreduce(ideal-gas diagnostic selected-field validation)");
  if (failed_rank >= 0)
    fail(DiagnosticFailureClass::invalid_request,
         "closure.diagnostics.selected-field", failed_rank,
         "ideal-gas diagnostic selected fields are invalid");
  failed_rank =
      lowest_rank(mpi, request_invalid,
                  "MPI_Allreduce(ideal-gas diagnostic request validation)");
  if (failed_rank >= 0)
    fail(DiagnosticFailureClass::invalid_request, "closure.diagnostics.frame",
         failed_rank, "ideal-gas diagnostic frame is invalid");
  require_byte_agreement(mpi, request_key(request),
                         DiagnosticFailureClass::invalid_request,
                         "closure.diagnostics.request-agreement",
                         "ideal-gas diagnostic requests disagree");
  require_byte_agreement(mpi, provider_key(source),
                         DiagnosticFailureClass::invalid_input,
                         "closure.diagnostics.provider-agreement",
                         "ideal-gas diagnostic providers disagree");
  require_exact_cell_cover(source, mpi);
  bool provider_failed = false;
  std::optional<CommittedObservation> observation;
  try {
    observation.emplace(observe(source, request.level));
  } catch (...) {
    provider_failed = true;
  }
  const int local_rank = provider_failed ? mpi.rank() : mpi.size();
  failed_rank = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&local_rank, &failed_rank, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(ideal-gas diagnostic provider)");
  if (failed_rank != mpi.size())
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.provider-agreement", failed_rank,
         "ideal-gas diagnostic providers disagree");
  if (request.level == DiagnosticLevel::summary ||
      request.level == DiagnosticLevel::invariants) {
    double sums[2]{observation->mass, observation->energy};
    if (request.level == DiagnosticLevel::summary)
      runtime::check_mpi_result(
          MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, mpi.comm()),
          "MPI_Allreduce(ideal-gas diagnostic sums)");
    double maxima[7]{-observation->rho_min, observation->rho_max,
                     -observation->h_min,   observation->h_max,
                     -observation->t_min,   observation->t_max,
                     observation->eos_max};
    runtime::check_mpi_result(
        MPI_Allreduce(MPI_IN_PLACE, maxima, 7, MPI_DOUBLE, MPI_MAX, mpi.comm()),
        "MPI_Allreduce(ideal-gas diagnostic extrema)");
    observation->mass = sums[0];
    observation->energy = sums[1];
    observation->rho_min = -maxima[0];
    observation->rho_max = maxima[1];
    observation->h_min = -maxima[2];
    observation->h_max = maxima[3];
    observation->t_min = -maxima[4];
    observation->t_max = maxima[5];
    observation->eos_max = maxima[6];
  }
  std::uint64_t fingerprint_xor = observation->fingerprint.xor64;
  std::uint64_t fingerprint_sum = observation->fingerprint.sum64;
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, &fingerprint_xor, 1, MPI_UINT64_T, MPI_BXOR,
                    mpi.comm()),
      "MPI_Allreduce(ideal-gas diagnostic fingerprint xor)");
  runtime::check_mpi_result(
      MPI_Allreduce(MPI_IN_PLACE, &fingerprint_sum, 1, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(ideal-gas diagnostic fingerprint sum)");
  observation->fingerprint = {fingerprint_xor, fingerprint_sum};
  auto record =
      build_record(source, request, *observation, DiagnosticScope::collective);
  aggregate_samples(record, mpi);
  submit_collective(record, source, request, mpi, sink);
} catch (const DiagnosticCollectionError &) {
  throw;
} catch (const runtime::MpiOperationError &) {
  fail(DiagnosticFailureClass::collective_operation,
       "closure.diagnostics.collective-operation", -1,
       "ideal-gas diagnostic collective operation failed");
} catch (const runtime::Error &) {
  fail(DiagnosticFailureClass::invalid_input,
       "closure.diagnostics.stale-source", -1,
       "ideal-gas diagnostic source is stale");
}

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
void test::IdealGasClosureDiagnosticTestAccess::set_fault(
    IdealGasClosureDiagnosticFault value, int rank) noexcept {
  injected_fault = value;
  injected_fault_rank = rank;
}

void test::IdealGasClosureDiagnosticTestAccess::reset() noexcept {
  injected_fault = IdealGasClosureDiagnosticFault::none;
  injected_fault_rank = -1;
  diagnostic_work = {};
}

test::IdealGasClosureDiagnosticWork
test::IdealGasClosureDiagnosticTestAccess::work() noexcept {
  return diagnostic_work;
}

DiagnosticInvariant
test::IdealGasClosureDiagnosticTestAccess::positive_invariant(std::string id,
                                                              std::string unit,
                                                              double observed) {
  return make_invariant(std::move(id), std::move(unit), observed, 0.0,
                        InvariantRelation::positive);
}
#endif

} // namespace hundun::diagnostics
