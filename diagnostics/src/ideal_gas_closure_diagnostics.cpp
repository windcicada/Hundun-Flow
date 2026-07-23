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
auto &test_state(const flow::IdealGasClosureDiagnosticSource &source) {
  return flow::detail::DensityClosureDiagnosticAccess::test_state(source);
}

bool fault(const flow::IdealGasClosureDiagnosticSource &source,
           test::IdealGasClosureDiagnosticFault value, int rank) noexcept {
  auto &state = test_state(source);
  const bool reached =
      state.fault == static_cast<std::uint8_t>(value) && !state.consumed;
  const bool selected =
      reached && (state.fault_rank < 0 || state.fault_rank == rank);
  if (reached)
    state.consumed = true;
  return selected;
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
            std::string_view("diagnostics.request.selected-fields") ||
        error.code() == std::string_view("diagnostics.selected-field.invalid"))
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
                            const char *code, const char *message,
                            bool force_oversized = false) {
  std::uint64_t size =
      mpi.rank() == 0
          ? (force_oversized
                 ? static_cast<std::uint64_t>(
                       std::numeric_limits<int>::max()) +
                       1U
                 : static_cast<std::uint64_t>(local.size()))
          : 0U;
  runtime::check_mpi_result(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
                            "MPI_Bcast(ideal-gas diagnostic agreement size)");
  if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.aggregation",
         0, "ideal-gas diagnostic agreement is too large");
  std::vector<unsigned char> reference;
  bool allocation_failed = false;
  try {
    reference.resize(static_cast<std::size_t>(size));
  } catch (...) {
    allocation_failed = true;
  }
  const int allocation_rank = lowest_rank(
      mpi, allocation_failed,
      "MPI_Allreduce(ideal-gas diagnostic agreement allocation)");
  if (allocation_rank >= 0)
    fail(classification, code, allocation_rank,
         "ideal-gas diagnostic agreement preparation failed");
  if (mpi.rank() == 0)
    std::copy(local.begin(), local.end(), reference.begin());
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
  runtime::check_mpi_result(MPI_Comm_compare(bound.comm(), mpi.comm(),
                                             &comparison),
                            "MPI_Comm_compare(ideal-gas diagnostics)");
  const bool mismatch = comparison != MPI_IDENT || bound.rank() != mpi.rank() ||
                        bound.size() != mpi.size();
  const int rank = lowest_rank(
      bound, mismatch,
      "MPI_Allreduce(ideal-gas diagnostic communicator validation)");
  if (rank >= 0)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout", rank,
         "ideal-gas diagnostic communicator does not match the source");
}

void require_exact_cell_cover(
    const flow::IdealGasClosureDiagnosticSource &source,
    const runtime::MpiContext &mpi) {
  struct Wire final {
    std::int64_t begin[3];
    std::int64_t end[3];
    std::int64_t global[3];
    std::uint64_t count;
  };
  std::optional<Wire> prepared_local;
  std::vector<Wire> all;
  bool preparation_failed = false;
  try {
    const auto session =
        flow::detail::DensityClosureDiagnosticAccess::acquire_committed(source);
    const auto box = session.owned_box;
    const auto global = session.global_extent;
    prepared_local.emplace(
        Wire{{box.begin.x, box.begin.y, box.begin.z},
             {box.end.x, box.end.y, box.end.z},
             {global.x, global.y, global.z}, source.owned_cell_count()});
    all.resize(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    preparation_failed = true;
  }
  const int preparation_rank = lowest_rank(
      mpi, preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic layout preparation)");
  if (preparation_rank >= 0)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout",
         preparation_rank, "ideal-gas diagnostic layout preparation failed");
  Wire local = *prepared_local;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(source, test::IdealGasClosureDiagnosticFault::ownership_layout,
            mpi.rank()))
    ++local.count;
  if (fault(source, test::IdealGasClosureDiagnosticFault::global_extent,
            mpi.rank()))
    ++local.global[0];
#endif
  runtime::check_mpi_result(
      MPI_Allgather(&local, static_cast<int>(sizeof(Wire)), MPI_BYTE,
                    all.data(), static_cast<int>(sizeof(Wire)), MPI_BYTE,
                    mpi.comm()),
      "MPI_Allgather(ideal-gas diagnostic ownership)");
  std::uint64_t covered{};
  const auto &canonical_global = all.front().global;
  for (std::size_t left = 0; left < all.size(); ++left) {
    std::uint64_t volume = 1U;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      if (all[left].global[axis] != canonical_global[axis] ||
          all[left].global[axis] <= 0 || all[left].begin[axis] < 0 ||
          all[left].end[axis] <= all[left].begin[axis] ||
          all[left].end[axis] > all[left].global[axis])
        fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout",
             static_cast<int>(left),
             "ideal-gas diagnostic ownership is invalid");
      const auto width = static_cast<std::uint64_t>(all[left].end[axis] -
                                                    all[left].begin[axis]);
      if (volume > std::numeric_limits<std::uint64_t>::max() / width)
        fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout",
             static_cast<int>(left),
             "ideal-gas diagnostic ownership overflows");
      volume *= width;
    }
    if (volume != all[left].count ||
        covered > std::numeric_limits<std::uint64_t>::max() - volume)
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout",
           static_cast<int>(left),
           "ideal-gas diagnostic ownership count is invalid");
    covered += volume;
    for (std::size_t right = 0; right < left; ++right) {
      bool overlap = true;
      for (std::size_t axis = 0; axis < 3U; ++axis)
        overlap = overlap && all[left].begin[axis] < all[right].end[axis] &&
                  all[right].begin[axis] < all[left].end[axis];
      if (overlap)
        fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout",
             static_cast<int>(left),
             "ideal-gas diagnostic ownership overlaps");
    }
  }
  std::uint64_t global_count = 1U;
  for (const auto extent : canonical_global) {
    const auto width = static_cast<std::uint64_t>(extent);
    if (global_count > std::numeric_limits<std::uint64_t>::max() / width)
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout", 0,
           "ideal-gas diagnostic global count overflows");
    global_count *= width;
  }
  if (covered != global_count)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout", 0,
         "ideal-gas diagnostic ownership does not exactly cover");
}

struct CommittedObservation final {
  double mass{};
  double rho_min{std::numeric_limits<double>::infinity()};
  double rho_max{-std::numeric_limits<double>::infinity()};
  double h_min{std::numeric_limits<double>::infinity()};
  double h_max{-std::numeric_limits<double>::infinity()};
  double t_min{std::numeric_limits<double>::infinity()};
  double t_max{-std::numeric_limits<double>::infinity()};
  double eos_max{};
  DiagnosticFingerprintParts fingerprint{};
  std::uint64_t eligible_sample_count{};
  std::vector<DiagnosticSample> retained_samples;
};

CommittedObservation
observe(const flow::IdealGasClosureDiagnosticSource &source,
        const DiagnosticRequest &request) {
  const auto level = request.level;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  auto &work = test_state(source);
  ++work.observations;
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
    ++work.fingerprint_items;
#endif
  }
  const auto cell_indices = [&](std::size_t cell) {
    const auto plane =
        static_cast<std::size_t>(extent.x) * static_cast<std::size_t>(extent.y);
    const int k = static_cast<int>(cell / plane);
    const auto within = cell % plane;
    const int j = static_cast<int>(within / static_cast<std::size_t>(extent.x));
    const int i = static_cast<int>(within % static_cast<std::size_t>(extent.x));
    return std::array<int, 3>{i, j, k};
  };
  for (std::size_t cell = 0; cell < count; ++cell) {
    const auto index = cell_indices(cell);
    const int i = index[0], j = index[1], k = index[2];
    const double rho = session.density(i, j, k, 0);
    const double rho_h = session.enthalpy_density(i, j, k, 0);
    const auto global_id = source.fingerprint_field_global_id(1U, cell);
    fingerprint.add("rho", global_id, 0U, describe_fp64(rho));
    fingerprint.add("rho_h", global_id, 0U, describe_fp64(rho_h));
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    work.fingerprint_items += 2U;
#endif
    if (level != DiagnosticLevel::summary &&
        level != DiagnosticLevel::invariants)
      continue;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (level == DiagnosticLevel::summary)
      ++work.summary_items;
    else
      ++work.invariant_items;
#endif
    const double h = rho_h / rho;
    const double temperature = h / cp;
    const double volume = source.cell_volume_m3(cell);
    if (level == DiagnosticLevel::summary) {
      result.mass += volume * rho;
    }
    result.rho_min = std::min(result.rho_min, rho);
    result.h_min = std::min(result.h_min, h);
    result.t_min = std::min(result.t_min, temperature);
    if (level == DiagnosticLevel::summary) {
      result.rho_max = std::max(result.rho_max, rho);
      result.h_max = std::max(result.h_max, h);
      result.t_max = std::max(result.t_max, temperature);
    }
    if (level == DiagnosticLevel::invariants)
      result.eos_max = std::max(
          result.eos_max,
          relative_product_error(rho, gas_constant, temperature, p0));
  }
  if (level == DiagnosticLevel::bounded_state_sample) {
    const auto selected = request.selected_fields.empty()
                              ? std::vector<std::string_view>(
                                    kSampleIds.begin(), kSampleIds.end())
                              : request.selected_fields;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    ++work.allocation_events;
#endif
    result.retained_samples.reserve(
        std::min(request.sample_budget, selected.size() * count + 1U));
    for (const auto id : selected) {
      const auto field = static_cast<std::size_t>(
          std::find(kSampleIds.begin(), kSampleIds.end(), id) -
          kSampleIds.begin());
      const auto eligible = source.sample_field_item_count(field);
      if (eligible > std::numeric_limits<std::uint64_t>::max() -
                         result.eligible_sample_count)
        fail(DiagnosticFailureClass::layout,
             "closure.diagnostics.sample-preparation", source.relative_rank(),
             "ideal-gas eligible sample count overflows");
      result.eligible_sample_count += eligible;
      for (std::size_t item = 0;
           item < eligible &&
           result.retained_samples.size() < request.sample_budget;
           ++item) {
        double value = p0;
        if (field != 1U) {
          const auto index = cell_indices(item);
          const double rho = session.density(index[0], index[1], index[2], 0);
          const double rho_h =
              session.enthalpy_density(index[0], index[1], index[2], 0);
          const double enthalpy = rho_h / rho;
          value = field == 0U   ? enthalpy
                  : field == 2U ? rho
                  : field == 3U ? rho_h
                                : enthalpy / cp;
        }
        result.retained_samples.push_back(
            {std::string(id), source.sample_field_global_id(field, item), 0U,
             std::string(source.sample_field_unit(field)), describe_fp64(value)});
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
        ++work.sample_items;
        ++work.retained_sample_items;
#endif
      }
    }
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
    record.sample_budget = request.sample_budget;
    record.eligible_sample_count = observation.eligible_sample_count;
    record.samples = observation.retained_samples;
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

void require_exact_record_contract(
    const DiagnosticRecord &record, const DiagnosticRecord &expected,
    const flow::IdealGasClosureDiagnosticSource &source,
    const DiagnosticRequest &request) {
  validate(record, describe_diagnostics(source), request);
  validate(expected, describe_diagnostics(source), request);
  if (to_canonical_json(record) != to_canonical_json(expected))
    throw runtime::Error("ideal-gas diagnostic record contract differs");
}

void submit(DiagnosticRecord record,
            const flow::IdealGasClosureDiagnosticSource &source,
            const DiagnosticRequest &request, DiagnosticSink &sink) {
  const auto expected = record;
#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  (void)source;
#endif
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(source, test::IdealGasClosureDiagnosticFault::record_validation,
            source.relative_rank()))
    record.schema_version = 0U;
#endif
  try {
    require_exact_record_contract(record, expected, source, request);
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
                       const flow::IdealGasClosureDiagnosticSource &source,
                       const runtime::MpiContext &mpi) {
#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  (void)source;
#endif
  if (record.level != DiagnosticLevel::bounded_state_sample)
    return;
  struct Wire final {
    std::uint64_t global_id{};
    std::uint64_t value_bits{};
    std::uint32_t field{};
    std::uint8_t status{};
    std::int32_t owner_rank{};
  };
  std::optional<std::vector<Wire>> prepared_local;
  std::optional<std::vector<int>> prepared_counts;
  bool wire_invalid = false;
  bool local_preparation_failed = false;
  try {
    prepared_local.emplace();
    prepared_local->reserve(record.samples.size());
    for (const auto &sample : record.samples) {
      const auto found =
          std::find(kSampleIds.begin(), kSampleIds.end(), sample.field_id);
      if (found == kSampleIds.end()) {
        wire_invalid = true;
        break;
      }
      prepared_local->push_back(
          {sample.global_id, sample.value.bits,
           static_cast<std::uint32_t>(found - kSampleIds.begin()),
           static_cast<std::uint8_t>(sample.value.status), mpi.rank()});
    }
    prepared_counts.emplace(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    local_preparation_failed = true;
  }
  int failure_rank = lowest_rank(
      mpi, wire_invalid,
      "MPI_Allreduce(ideal-gas diagnostic local sample validation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
         failure_rank, "ideal-gas diagnostic sample field is invalid");
  failure_rank = lowest_rank(
      mpi, local_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic local sample preparation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas diagnostic sample preparation failed");
  auto local = std::move(*prepared_local);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (!local.empty() && fault(source,
                              test::IdealGasClosureDiagnosticFault::sample_wire,
                              mpi.rank()))
    local.push_back(local.front());
#endif
  const bool payload_too_large =
      local.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(Wire);
  failure_rank = lowest_rank(
      mpi, payload_too_large,
      "MPI_Allreduce(ideal-gas diagnostic sample payload validation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas diagnostic sample payload is too large");
  const int local_bytes = static_cast<int>(local.size() * sizeof(Wire));
  auto counts = std::move(*prepared_counts);
  runtime::check_mpi_result(MPI_Allgather(&local_bytes, 1, MPI_INT,
                                          counts.data(), 1, MPI_INT,
                                          mpi.comm()),
                            "MPI_Allgather(ideal-gas diagnostic sample sizes)");
  std::optional<std::vector<int>> prepared_displacements;
  std::optional<std::vector<unsigned char>> prepared_merged_bytes;
  bool sizes_invalid = false;
  bool receive_preparation_failed = false;
  int total{};
  try {
    prepared_displacements.emplace(counts.size());
    for (std::size_t rank = 0; rank < counts.size(); ++rank) {
      if (counts[rank] < 0 ||
          counts[rank] % static_cast<int>(sizeof(Wire)) != 0 ||
          counts[rank] > std::numeric_limits<int>::max() - total) {
        sizes_invalid = true;
        break;
      }
      (*prepared_displacements)[rank] = total;
      total += counts[rank];
    }
    if (!sizes_invalid)
      prepared_merged_bytes.emplace(static_cast<std::size_t>(total));
  } catch (...) {
    receive_preparation_failed = true;
  }
  failure_rank = lowest_rank(
      mpi, sizes_invalid,
      "MPI_Allreduce(ideal-gas diagnostic sample size validation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas diagnostic sample size is invalid");
  failure_rank = lowest_rank(
      mpi, receive_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic sample receive preparation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas diagnostic sample preparation failed");
  auto displacements = std::move(*prepared_displacements);
  auto merged_bytes = std::move(*prepared_merged_bytes);
  runtime::check_mpi_result(
      MPI_Allgatherv(local.data(), local_bytes, MPI_BYTE, merged_bytes.data(),
                     counts.data(), displacements.data(), MPI_BYTE, mpi.comm()),
      "MPI_Allgatherv(ideal-gas diagnostic samples)");
  std::optional<std::vector<Wire>> prepared_merged;
  std::optional<std::vector<std::uint64_t>> prepared_eligible_by_rank;
  int duplicate_owner = -1;
  int invalid_owner = -1;
  bool merge_preparation_failed = false;
  try {
    prepared_merged.emplace(merged_bytes.size() / sizeof(Wire));
    if (!prepared_merged->empty())
      std::memcpy(prepared_merged->data(), merged_bytes.data(),
                  merged_bytes.size());
    std::sort(prepared_merged->begin(), prepared_merged->end(),
              [](const Wire &left, const Wire &right) {
                return std::tie(left.field, left.global_id) <
                       std::tie(right.field, right.global_id);
              });
    const auto duplicate =
        std::adjacent_find(prepared_merged->begin(), prepared_merged->end(),
                           [](const Wire &left, const Wire &right) {
                             return left.field == right.field &&
                                    left.global_id == right.global_id;
                           });
    if (duplicate != prepared_merged->end())
      duplicate_owner = duplicate->owner_rank;
    for (const auto &item : *prepared_merged) {
      if (item.field >= kSampleIds.size() ||
          item.status >
              static_cast<std::uint8_t>(DiagnosticValueStatus::unavailable)) {
        invalid_owner = item.owner_rank;
        break;
      }
    }
    prepared_eligible_by_rank.emplace(static_cast<std::size_t>(mpi.size()));
  } catch (...) {
    merge_preparation_failed = true;
  }
  failure_rank = lowest_rank(
      mpi, duplicate_owner >= 0 && mpi.rank() == duplicate_owner,
      "MPI_Allreduce(ideal-gas diagnostic duplicate sample owner)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
         failure_rank, "ideal-gas diagnostic samples are duplicated");
  failure_rank = lowest_rank(
      mpi, invalid_owner >= 0 && mpi.rank() == invalid_owner,
      "MPI_Allreduce(ideal-gas diagnostic invalid sample owner)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
         failure_rank, "ideal-gas diagnostic sample wire is invalid");
  failure_rank = lowest_rank(
      mpi, merge_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic sample merge preparation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas diagnostic sample preparation failed");
  auto merged = std::move(*prepared_merged);
  auto eligible_by_rank = std::move(*prepared_eligible_by_rank);
  runtime::check_mpi_result(
      MPI_Allgather(&record.eligible_sample_count, 1, MPI_UINT64_T,
                    eligible_by_rank.data(), 1, MPI_UINT64_T, mpi.comm()),
      "MPI_Allgather(ideal-gas diagnostic eligible samples)");
  std::uint64_t eligible{};
  bool eligible_invalid = false;
  for (const auto value : eligible_by_rank) {
    if (value > std::numeric_limits<std::uint64_t>::max() - eligible) {
      eligible_invalid = true;
      break;
    }
    eligible += value;
  }
  failure_rank = lowest_rank(
      mpi, eligible_invalid,
      "MPI_Allreduce(ideal-gas diagnostic eligible sample validation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas eligible sample count overflows");
  if (merged.size() > record.sample_budget)
    merged.resize(record.sample_budget);
  std::optional<std::vector<DiagnosticSample>> prepared_samples;
  bool record_preparation_failed = false;
  try {
    prepared_samples.emplace();
    prepared_samples->reserve(merged.size());
    for (const auto &item : merged) {
      const std::size_t field = item.field;
      prepared_samples->push_back(
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
  } catch (...) {
    record_preparation_failed = true;
  }
  failure_rank = lowest_rank(
      mpi, record_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic sample record preparation)");
  if (failure_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", failure_rank,
         "ideal-gas diagnostic sample preparation failed");
  record.eligible_sample_count = eligible;
  record.samples = std::move(*prepared_samples);
  record.samples_truncated =
      record.eligible_sample_count > record.samples.size();
}

void submit_collective(const DiagnosticRecord &record,
                       const flow::IdealGasClosureDiagnosticSource &source,
                       const DiagnosticRequest &request,
                       const runtime::MpiContext &mpi, DiagnosticSink &sink) {
  auto candidate = record;
  bool invalid = false;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  if (fault(source, test::IdealGasClosureDiagnosticFault::record_validation,
            mpi.rank()))
    candidate.schema_version = 0U;
#endif
  try {
    require_exact_record_contract(candidate, record, source, request);
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
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (fault(source, test::IdealGasClosureDiagnosticFault::stale_source,
              source.relative_rank()))
      fail(DiagnosticFailureClass::invalid_input,
           "closure.diagnostics.stale-source", -1,
           "ideal-gas diagnostic source is stale");
    if (fault(source, test::IdealGasClosureDiagnosticFault::capability,
              source.relative_rank()))
      fail(DiagnosticFailureClass::capability, "closure.diagnostics.capability",
           -1, "ideal-gas diagnostic capability is unavailable");
    if (fault(source, test::IdealGasClosureDiagnosticFault::provider_agreement,
              source.relative_rank()))
      fail(DiagnosticFailureClass::invalid_input,
           "closure.diagnostics.provider-agreement", -1,
           "ideal-gas diagnostic provider is inconsistent");
    if (fault(source, test::IdealGasClosureDiagnosticFault::ownership_layout,
              source.relative_rank()) ||
        fault(source, test::IdealGasClosureDiagnosticFault::global_extent,
              source.relative_rank()))
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.layout", -1,
           "ideal-gas diagnostic layout is inconsistent");
    if (fault(source, test::IdealGasClosureDiagnosticFault::request_preparation,
              source.relative_rank()))
      fail(DiagnosticFailureClass::invalid_request,
           "closure.diagnostics.request-agreement", -1,
           "ideal-gas diagnostic request preparation failed");
    if (fault(source, test::IdealGasClosureDiagnosticFault::sample_preparation,
              source.relative_rank()))
      fail(DiagnosticFailureClass::layout,
           "closure.diagnostics.sample-preparation", -1,
           "ideal-gas diagnostic sample preparation failed");
    if (fault(source, test::IdealGasClosureDiagnosticFault::sample_wire,
              source.relative_rank()))
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.sample-wire",
           -1, "ideal-gas diagnostic sample wire is malformed");
    if (fault(source, test::IdealGasClosureDiagnosticFault::aggregation,
              source.relative_rank()) ||
        fault(source, test::IdealGasClosureDiagnosticFault::oversized_agreement,
              source.relative_rank()))
      fail(DiagnosticFailureClass::layout, "closure.diagnostics.aggregation",
           -1, "ideal-gas diagnostic aggregation failed");
#endif
    require_request(source, request, DiagnosticScope::local);
    submit(build_record(source, request, observe(source, request),
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
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  auto &work = test_state(source);
  ++work.collective_calls;
  int injected_rank = lowest_rank(
      mpi,
      fault(source, test::IdealGasClosureDiagnosticFault::stale_source,
            mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic stale-source injection)");
  if (injected_rank >= 0)
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.stale-source", injected_rank,
         "ideal-gas diagnostic source is stale");
  injected_rank = lowest_rank(
      mpi,
      fault(source, test::IdealGasClosureDiagnosticFault::request_preparation,
            mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic request preparation)");
  if (injected_rank >= 0)
    fail(DiagnosticFailureClass::invalid_request,
         "closure.diagnostics.request-agreement", injected_rank,
         "ideal-gas diagnostic request preparation failed");
  injected_rank = lowest_rank(
      mpi, fault(source, test::IdealGasClosureDiagnosticFault::capability,
                 mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic capability)");
  if (injected_rank >= 0)
    fail(DiagnosticFailureClass::capability,
         "closure.diagnostics.capability", injected_rank,
         "ideal-gas diagnostic capability is unavailable");
  injected_rank = lowest_rank(
      mpi, fault(source, test::IdealGasClosureDiagnosticFault::raw_mpi,
                 mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic typed failure preparation)");
  if (injected_rank >= 0)
    runtime::check_mpi_result(MPI_ERR_OTHER,
                              "MPI_Task21(ideal-gas diagnostic typed fault)");
#endif
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
  std::optional<std::vector<unsigned char>> prepared_request;
  bool request_preparation_failed = false;
  try {
    prepared_request.emplace(request_key(request));
  } catch (...) {
    request_preparation_failed = true;
  }
  failed_rank = lowest_rank(
      mpi, request_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic request serialization)");
  if (failed_rank >= 0)
    fail(DiagnosticFailureClass::invalid_request,
         "closure.diagnostics.request-agreement", failed_rank,
         "ideal-gas diagnostic request preparation failed");
  bool force_oversized = false;
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  force_oversized =
      lowest_rank(mpi,
                  fault(source,
                        test::IdealGasClosureDiagnosticFault::oversized_agreement,
                        mpi.rank()),
                  "MPI_Allreduce(ideal-gas diagnostic oversized agreement)") >=
      0;
#endif
  require_byte_agreement(mpi, *prepared_request,
                         DiagnosticFailureClass::invalid_request,
                         "closure.diagnostics.request-agreement",
                         "ideal-gas diagnostic requests disagree",
                         force_oversized);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  injected_rank = lowest_rank(
      mpi,
      fault(source, test::IdealGasClosureDiagnosticFault::provider_agreement,
            mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic provider preparation)");
  if (injected_rank >= 0)
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.provider-agreement", injected_rank,
         "ideal-gas diagnostic provider preparation failed");
#endif
  std::optional<std::vector<unsigned char>> prepared_provider;
  bool provider_preparation_failed = false;
  try {
    prepared_provider.emplace(provider_key(source));
  } catch (...) {
    provider_preparation_failed = true;
  }
  failed_rank = lowest_rank(
      mpi, provider_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic provider serialization)");
  if (failed_rank >= 0)
    fail(DiagnosticFailureClass::invalid_input,
         "closure.diagnostics.provider-agreement", failed_rank,
         "ideal-gas diagnostic provider preparation failed");
  require_byte_agreement(mpi, *prepared_provider,
                         DiagnosticFailureClass::invalid_input,
                         "closure.diagnostics.provider-agreement",
                         "ideal-gas diagnostic providers disagree");
  require_exact_cell_cover(source, mpi);
  bool provider_failed = false;
  std::optional<CommittedObservation> observation;
  try {
    observation.emplace(observe(source, request));
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
    if (request.level == DiagnosticLevel::summary) {
      runtime::check_mpi_result(
          MPI_Allreduce(MPI_IN_PLACE, &observation->mass, 1, MPI_DOUBLE,
                        MPI_SUM, mpi.comm()),
          "MPI_Allreduce(ideal-gas diagnostic mass)");
      double extrema[6]{-observation->rho_min, observation->rho_max,
                        -observation->h_min,   observation->h_max,
                        -observation->t_min,   observation->t_max};
      runtime::check_mpi_result(
          MPI_Allreduce(MPI_IN_PLACE, extrema, 6, MPI_DOUBLE, MPI_MAX,
                        mpi.comm()),
          "MPI_Allreduce(ideal-gas diagnostic summary extrema)");
      observation->rho_min = -extrema[0];
      observation->rho_max = extrema[1];
      observation->h_min = -extrema[2];
      observation->h_max = extrema[3];
      observation->t_min = -extrema[4];
      observation->t_max = extrema[5];
    } else {
      double invariants[4]{-observation->rho_min, -observation->h_min,
                           -observation->t_min, observation->eos_max};
      runtime::check_mpi_result(
          MPI_Allreduce(MPI_IN_PLACE, invariants, 4, MPI_DOUBLE, MPI_MAX,
                        mpi.comm()),
          "MPI_Allreduce(ideal-gas diagnostic invariant operands)");
      observation->rho_min = -invariants[0];
      observation->h_min = -invariants[1];
      observation->t_min = -invariants[2];
      observation->eos_max = invariants[3];
    }
  }
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  injected_rank = lowest_rank(
      mpi,
      fault(source, test::IdealGasClosureDiagnosticFault::aggregation,
            mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic aggregation preparation)");
  if (injected_rank >= 0)
    fail(DiagnosticFailureClass::layout, "closure.diagnostics.aggregation",
         injected_rank, "ideal-gas diagnostic aggregation failed");
#endif
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
  std::optional<DiagnosticRecord> prepared_record;
  bool record_preparation_failed = false;
  try {
    prepared_record.emplace(build_record(source, request, *observation,
                                         DiagnosticScope::collective));
  } catch (...) {
    record_preparation_failed = true;
  }
  failed_rank = lowest_rank(
      mpi, record_preparation_failed,
      "MPI_Allreduce(ideal-gas diagnostic record preparation)");
  if (failed_rank >= 0)
    fail(DiagnosticFailureClass::invalid_input, "closure.diagnostics.record",
         failed_rank, "ideal-gas diagnostic record preparation failed");
  auto record = std::move(*prepared_record);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  injected_rank = lowest_rank(
      mpi,
      fault(source, test::IdealGasClosureDiagnosticFault::sample_preparation,
            mpi.rank()),
      "MPI_Allreduce(ideal-gas diagnostic sample preparation)");
  if (injected_rank >= 0)
    fail(DiagnosticFailureClass::layout,
         "closure.diagnostics.sample-preparation", injected_rank,
         "ideal-gas diagnostic sample preparation failed");
#endif
  aggregate_samples(record, source, mpi);
  submit_collective(record, source, request, mpi, sink);
} catch (const DiagnosticCollectionError &) {
  throw;
} catch (const runtime::MpiOperationError &) {
  throw;
} catch (const runtime::Error &) {
  fail(DiagnosticFailureClass::invalid_input,
       "closure.diagnostics.stale-source", -1,
       "ideal-gas diagnostic source is stale");
}

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
void test::IdealGasClosureDiagnosticTestAccess::set_fault(
    const flow::IdealGasClosureDiagnosticSource &source,
    IdealGasClosureDiagnosticFault value, int rank) noexcept {
  auto &state = test_state(source);
  state.fault = static_cast<std::uint8_t>(value);
  state.fault_rank = rank;
  state.consumed = false;
}

void test::IdealGasClosureDiagnosticTestAccess::reset(
    const flow::IdealGasClosureDiagnosticSource &source) noexcept {
  test_state(source) = {};
  test_state(source).fault_rank = -1;
}

test::IdealGasClosureDiagnosticWork
test::IdealGasClosureDiagnosticTestAccess::work(
    const flow::IdealGasClosureDiagnosticSource &source) noexcept {
  const auto &state = test_state(source);
  return {state.observations,
          state.fingerprint_items,
          state.summary_items,
          state.invariant_items,
          state.sample_items,
          state.retained_sample_items,
          state.allocation_events,
          state.full_field_copy_attempts,
          state.collective_calls};
}

DiagnosticInvariant
test::IdealGasClosureDiagnosticTestAccess::positive_invariant(std::string id,
                                                              std::string unit,
                                                              double observed) {
  return make_invariant(std::move(id), std::move(unit), observed, 0.0,
                        InvariantRelation::positive);
}

void test::IdealGasClosureDiagnosticTestAccess::submit_record_contract(
    const flow::IdealGasClosureDiagnosticSource &source,
    const DiagnosticRequest &request, const DiagnosticRecord &expected,
    const DiagnosticRecord &candidate, DiagnosticSink &sink) {
  try {
    require_exact_record_contract(candidate, expected, source, request);
  } catch (...) {
    fail(DiagnosticFailureClass::invalid_input, "closure.diagnostics.record",
         -1, "ideal-gas diagnostic record is invalid");
  }
  sink.submit(candidate);
}
#endif

} // namespace hundun::diagnostics
