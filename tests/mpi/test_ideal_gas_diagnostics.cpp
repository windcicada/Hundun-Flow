// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
#include "diagnostics/src/ideal_gas_closure_diagnostics_test_access.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

namespace {
void run_task21_diagnostic_matrix(
    const hundun::flow::IdealGasClosureDiagnosticSource &,
    const hundun::runtime::MpiContext &, const hundun::flow::FlowState &,
    const hundun::flow::IdealGasStepAttemptReport &);
void run_task21_status_matrix(
    const hundun::flow::IdealGasClosureDiagnosticSource &,
    const hundun::runtime::MpiContext &,
    hundun::diagnostics::DiagnosticStatus,
    hundun::diagnostics::DiagnosticFailureClass, std::string_view, int);
}

#define HUNDUN_TASK21_IDEAL_GAS_FIXTURE_ONLY 1
#define HUNDUN_TASK21_DIAGNOSTICS_EXHAUSTIVE 1
#define HUNDUN_TASK21_DIAGNOSTIC_MATRIX(source, mpi, state, report)           \
  run_task21_diagnostic_matrix(source, mpi, state, report)
#define HUNDUN_TASK21_DIAGNOSTIC_STATUS(source, mpi, status, classification,  \
                                        code, rank)                           \
  run_task21_status_matrix(source, mpi, status, classification, code, rank)
#include "tests/mpi/test_ideal_gas_piso.cpp"

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace {

namespace diagnostics = hundun::diagnostics;
using Access = diagnostics::test::IdealGasClosureDiagnosticTestAccess;
using Fault = diagnostics::test::IdealGasClosureDiagnosticFault;

std::string diagnostic_reference_mode;
std::filesystem::path diagnostic_reference_path;

class RecordingSink final : public diagnostics::DiagnosticSink {
public:
  void submit(const diagnostics::DiagnosticRecord &record) override {
    ++calls;
    if (throw_on_submit)
      throw std::runtime_error("injected ideal-gas diagnostic sink failure");
    records.push_back(record);
  }

  bool throw_on_submit{};
  std::size_t calls{};
  std::vector<diagnostics::DiagnosticRecord> records;
};

class FaultReset final {
public:
  explicit FaultReset(
      const hundun::flow::IdealGasClosureDiagnosticSource &source) noexcept
      : source_(&source) {}
  ~FaultReset() noexcept { Access::reset(*source_); }

private:
  const hundun::flow::IdealGasClosureDiagnosticSource *source_;
};

void write_diagnostic_reference(const std::string &canonical) {
  std::ofstream output(diagnostic_reference_path,
                       std::ios::binary | std::ios::trunc);
  if (!output)
    throw std::runtime_error("cannot create ideal-gas diagnostic reference");
  constexpr std::uint64_t magic = UINT64_C(0x48554e4449444941);
  const auto size = static_cast<std::uint64_t>(canonical.size());
  output.write(reinterpret_cast<const char *>(&magic), sizeof(magic));
  output.write(reinterpret_cast<const char *>(&size), sizeof(size));
  output.write(canonical.data(), static_cast<std::streamsize>(canonical.size()));
  if (!output)
    throw std::runtime_error("cannot write ideal-gas diagnostic reference");
}

std::string read_diagnostic_reference() {
  std::ifstream input(diagnostic_reference_path, std::ios::binary);
  if (!input)
    throw std::runtime_error("ideal-gas diagnostic reference is missing");
  std::uint64_t magic{}, size{};
  input.read(reinterpret_cast<char *>(&magic), sizeof(magic));
  input.read(reinterpret_cast<char *>(&size), sizeof(size));
  if (magic != UINT64_C(0x48554e4449444941) ||
      size > static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max()))
    throw std::runtime_error("ideal-gas diagnostic reference header differs");
  std::string canonical(static_cast<std::size_t>(size), '\0');
  input.read(canonical.data(), static_cast<std::streamsize>(size));
  if (!input || input.peek() != std::ifstream::traits_type::eof())
    throw std::runtime_error("ideal-gas diagnostic reference is malformed");
  return canonical;
}

struct ExactState final {
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues committed;
  hundun::flow::FlowLayerValues trial;
  hundun::flow::AcceptedStepMetadata metadata;
  hundun::flow::IdealGasClosureState closure;
  std::uint64_t source_generation{};
  hundun::runtime::Fp64ReductionCounters reductions;
};

ExactState capture(
    const hundun::flow::FlowState &state,
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    const hundun::runtime::MpiContext &mpi) {
  return {state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata(),
          source.closure_state(),
          hundun::flow::test::IdealGasClosureTestAccess::source_generation(
              source),
          mpi.fp64_reduction_counters()};
}

void require_pure(
    const ExactState &expected, const hundun::flow::FlowState &state,
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    const hundun::runtime::MpiContext &mpi) {
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.history, state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.committed, state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.trial, state.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      expected.metadata, state.metadata()));
  HUNDUN_CHECK(hundun::test::ideal_gas_closure_state_bitwise_equal(
      expected.closure, source.closure_state()));
  HUNDUN_CHECK(
      expected.source_generation ==
      hundun::flow::test::IdealGasClosureTestAccess::source_generation(source));
  const auto reductions = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(reductions.collective_calls ==
               expected.reductions.collective_calls);
  HUNDUN_CHECK(reductions.reduced_scalars ==
               expected.reductions.reduced_scalars);
  HUNDUN_CHECK(reductions.logical_payload_bytes ==
               expected.reductions.logical_payload_bytes);
}

diagnostics::DiagnosticRequest request_for(
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    diagnostics::DiagnosticLevel level, diagnostics::DiagnosticScope scope,
    std::size_t budget = 7U,
    std::vector<std::string_view> selected_fields = {}) {
  return {level,
          scope,
          {source.relative_rank(), source.committed_step(),
           source.committed_time_s(), "ideal-gas-closure.attempt-result"},
          std::move(selected_fields),
          level == diagnostics::DiagnosticLevel::bounded_state_sample ? budget
                                                                      : 0U};
}

void collect(const hundun::flow::IdealGasClosureDiagnosticSource &source,
             const hundun::runtime::MpiContext &mpi,
             const diagnostics::DiagnosticRequest &request,
             diagnostics::DiagnosticSink &sink) {
  if (request.scope == diagnostics::DiagnosticScope::local)
    diagnostics::collect_diagnostics(source, request, sink);
  else
    diagnostics::collect_diagnostics(source, mpi, request, sink);
}

template <class Items, class Expected>
void require_ids(const Items &items,
                 const Expected &expected) {
  HUNDUN_CHECK(items.size() == expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index)
    HUNDUN_CHECK(items[index].id == expected[index]);
}

void require_payload(const diagnostics::DiagnosticRecord &record) {
  struct MetricShape final {
    std::string_view id;
    diagnostics::DiagnosticMetricKind kind;
    std::string_view unit;
  };
  static const std::array metric_shapes{
      MetricShape{"closure.actual-mass",
                  diagnostics::DiagnosticMetricKind::conservation, "kg"},
      MetricShape{"closure.candidate-pressure",
                  diagnostics::DiagnosticMetricKind::state_summary, "Pa"},
      MetricShape{"closure.committed-mass",
                  diagnostics::DiagnosticMetricKind::conservation, "kg"},
      MetricShape{"closure.committed-pressure",
                  diagnostics::DiagnosticMetricKind::state_summary, "Pa"},
      MetricShape{"closure.density-maximum",
                  diagnostics::DiagnosticMetricKind::state_summary, "kg/m3"},
      MetricShape{"closure.density-minimum",
                  diagnostics::DiagnosticMetricKind::state_summary, "kg/m3"},
      MetricShape{"closure.enthalpy-maximum",
                  diagnostics::DiagnosticMetricKind::state_summary, "J/kg"},
      MetricShape{"closure.enthalpy-minimum",
                  diagnostics::DiagnosticMetricKind::state_summary, "J/kg"},
      MetricShape{"closure.enthalpy-temperature-error",
                  diagnostics::DiagnosticMetricKind::residual, "1"},
      MetricShape{"closure.eos-error",
                  diagnostics::DiagnosticMetricKind::residual, "1"},
      MetricShape{"closure.rho-h-remap-conservation",
                  diagnostics::DiagnosticMetricKind::conservation, "1"},
      MetricShape{"closure.rho-h-remap-l2",
                  diagnostics::DiagnosticMetricKind::residual, "1"},
      MetricShape{"closure.rho-remap-conservation",
                  diagnostics::DiagnosticMetricKind::conservation, "1"},
      MetricShape{"closure.rho-remap-l2",
                  diagnostics::DiagnosticMetricKind::residual, "1"},
      MetricShape{"closure.target-mass",
                  diagnostics::DiagnosticMetricKind::state_summary, "kg"},
      MetricShape{"closure.temperature-maximum",
                  diagnostics::DiagnosticMetricKind::state_summary, "K"},
      MetricShape{"closure.temperature-minimum",
                  diagnostics::DiagnosticMetricKind::state_summary, "K"}};
  struct InvariantShape final {
    std::string_view id;
    std::string_view unit;
    diagnostics::InvariantRelation relation;
    double limit;
    bool limit_available;
  };
  static const std::array invariant_shapes{
      InvariantShape{"closure.committed-eos", "1",
                     diagnostics::InvariantRelation::less_equal, 1.0e-12,
                     true},
      InvariantShape{"closure.final-evidence-available", "1",
                     diagnostics::InvariantRelation::equal, 1.0, true},
      InvariantShape{"closure.pressure-positive", "Pa",
                     diagnostics::InvariantRelation::positive, 0.0, false},
      InvariantShape{"density.positive", "kg/m3",
                     diagnostics::InvariantRelation::positive, 0.0, false},
      InvariantShape{"enthalpy.positive", "J/kg",
                     diagnostics::InvariantRelation::positive, 0.0, false},
      InvariantShape{"temperature.positive", "K",
                     diagnostics::InvariantRelation::positive, 0.0, false}};
  static const std::array<std::string_view, 3> counter_ids{
      "closure.collectives", "closure.evaluations", "closure.revision"};
  static const std::array<std::string_view, 6> identity_ids{
      "closure.attempt", "closure.state", "field.p0", "field.rho",
      "field.rho_h", "layout.cells"};

  HUNDUN_CHECK(record.identities.size() == identity_ids.size());
  for (std::size_t index = 0; index < identity_ids.size(); ++index) {
    HUNDUN_CHECK(record.identities[index].subject_id == identity_ids[index]);
    HUNDUN_CHECK(!record.identities[index].generation.has_value());
    HUNDUN_CHECK(!record.identities[index].allocation_identity.has_value());
  }
  HUNDUN_CHECK(!record.identities[0].layout_fingerprint.has_value());
  HUNDUN_CHECK(record.identities[0].revision.has_value());
  HUNDUN_CHECK(!record.identities[1].layout_fingerprint.has_value());
  HUNDUN_CHECK(record.identities[1].revision.has_value());
  HUNDUN_CHECK(record.identities[2].layout_fingerprint ==
               std::optional<std::string>{
                   "replicated-scalar.rank0-owner.v1"});
  HUNDUN_CHECK(!record.identities[2].revision.has_value());
  for (std::size_t index = 3U; index < identity_ids.size(); ++index) {
    HUNDUN_CHECK(record.identities[index].layout_fingerprint.has_value());
    HUNDUN_CHECK(!record.identities[index].layout_fingerprint->empty());
    HUNDUN_CHECK(!record.identities[index].revision.has_value());
  }
  HUNDUN_CHECK(record.identities[3].layout_fingerprint ==
               record.identities[4].layout_fingerprint);
  HUNDUN_CHECK(record.identities[4].layout_fingerprint ==
               record.identities[5].layout_fingerprint);
  HUNDUN_CHECK(record.time_s.status ==
               diagnostics::DiagnosticValueStatus::finite);
  if (record.level == diagnostics::DiagnosticLevel::summary) {
    HUNDUN_CHECK(record.metrics.size() == metric_shapes.size());
    for (std::size_t index = 0; index < metric_shapes.size(); ++index) {
      HUNDUN_CHECK(record.metrics[index].id == metric_shapes[index].id);
      HUNDUN_CHECK(record.metrics[index].kind == metric_shapes[index].kind);
      HUNDUN_CHECK(record.metrics[index].unit == metric_shapes[index].unit);
      HUNDUN_CHECK(record.metrics[index].value.status ==
                       diagnostics::DiagnosticValueStatus::finite ||
                   record.metrics[index].value.status ==
                       diagnostics::DiagnosticValueStatus::unavailable);
      if (record.metrics[index].value.status ==
          diagnostics::DiagnosticValueStatus::unavailable)
        HUNDUN_CHECK(record.metrics[index].value.bits == 0U);
    }
    HUNDUN_CHECK(record.invariants.empty());
    HUNDUN_CHECK(record.counters.empty());
    HUNDUN_CHECK(record.samples.empty());
  } else if (record.level == diagnostics::DiagnosticLevel::invariants) {
    HUNDUN_CHECK(record.invariants.size() == invariant_shapes.size());
    for (std::size_t index = 0; index < invariant_shapes.size(); ++index) {
      const auto &item = record.invariants[index];
      const auto &shape = invariant_shapes[index];
      HUNDUN_CHECK(item.id == shape.id);
      HUNDUN_CHECK(item.unit == shape.unit);
      HUNDUN_CHECK(item.relation == shape.relation);
      HUNDUN_CHECK(item.observed.status ==
                   diagnostics::DiagnosticValueStatus::finite);
      if (shape.limit_available) {
        HUNDUN_CHECK(item.limit.status ==
                     diagnostics::DiagnosticValueStatus::finite);
        HUNDUN_CHECK(item.limit.bits ==
                     diagnostics::describe_fp64(shape.limit).bits);
      } else {
        HUNDUN_CHECK(item.limit.status ==
                     diagnostics::DiagnosticValueStatus::unavailable);
        HUNDUN_CHECK(item.limit.bits == 0U);
      }
      HUNDUN_CHECK(item.passed == diagnostics::evaluate_invariant(item));
    }
    HUNDUN_CHECK(record.metrics.empty());
    HUNDUN_CHECK(record.counters.empty());
    HUNDUN_CHECK(record.samples.empty());
  } else if (record.level == diagnostics::DiagnosticLevel::counters) {
    require_ids(record.counters, counter_ids);
    for (const auto &counter : record.counters)
      HUNDUN_CHECK(counter.unit == "count");
    HUNDUN_CHECK(record.metrics.empty());
    HUNDUN_CHECK(record.invariants.empty());
    HUNDUN_CHECK(record.samples.empty());
  } else {
    HUNDUN_CHECK(record.metrics.empty());
    HUNDUN_CHECK(record.invariants.empty());
    HUNDUN_CHECK(record.counters.empty());
    HUNDUN_CHECK(record.samples.size() <= record.sample_budget);
    HUNDUN_CHECK(record.eligible_sample_count >= record.samples.size());
    HUNDUN_CHECK(std::is_sorted(
        record.samples.begin(), record.samples.end(),
        [](const auto &left, const auto &right) {
          return std::tie(left.field_id, left.global_id, left.component) <
                 std::tie(right.field_id, right.global_id, right.component);
        }));
    for (const auto &sample : record.samples) {
      HUNDUN_CHECK(sample.component == 0U);
      HUNDUN_CHECK(sample.value.status ==
                   diagnostics::DiagnosticValueStatus::finite);
      const std::string_view expected_unit =
          sample.field_id == "enthalpy"   ? "J/kg"
          : sample.field_id == "p0"       ? "Pa"
          : sample.field_id == "rho"      ? "kg/m3"
          : sample.field_id == "rho_h"    ? "J/m3"
          : sample.field_id == "temperature" ? "K"
                                             : "";
      HUNDUN_CHECK(!expected_unit.empty());
      HUNDUN_CHECK(sample.unit == expected_unit);
    }
  }
}

void require_success_payload_values(
    const diagnostics::DiagnosticRecord &record,
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    const hundun::runtime::MpiContext &mpi) {
  const auto &outer = source.report();
  const auto &closure_report = outer.closure_report();
  const auto closure_state = source.closure_state();
  HUNDUN_CHECK(record.identities[0].revision ==
               std::optional<std::uint64_t>{outer.attempt_identity()});
  HUNDUN_CHECK(record.identities[1].revision ==
               std::optional<std::uint64_t>{closure_state.revision});
  HUNDUN_CHECK(record.time_s.bits ==
               diagnostics::describe_fp64(source.committed_time_s()).bits);
  if (record.level == diagnostics::DiagnosticLevel::summary) {
    std::array<double, 17> expected{};
    expected[0] = closure_report.actual_mass_kg();
    expected[1] = closure_report.candidate_pressure_pa();
    double local_mass = 0.0;
    for (std::size_t cell = 0; cell < source.owned_cell_count(); ++cell)
      local_mass += source.sample_field_value(2U, cell) *
                    source.cell_volume_m3(cell);
    expected[2] = local_mass;
    if (record.scope == diagnostics::DiagnosticScope::collective)
      HUNDUN_CHECK(MPI_Allreduce(&local_mass, &expected[2], 1, MPI_DOUBLE,
                                 MPI_SUM, mpi.comm()) == MPI_SUCCESS);
    expected[3] = closure_state.thermodynamic_pressure_pa;
    expected[4] = source.sample_field_value(2U, 0U);
    expected[5] = expected[4];
    expected[6] = source.sample_field_value(0U, 0U);
    expected[7] = expected[6];
    expected[8] =
        closure_report.enthalpy_temperature_max_relative_error();
    expected[9] = closure_report.eos_max_relative_error();
    expected[10] =
        closure_report.rho_h_remap_relative_conservation_defect();
    expected[11] = closure_report.rho_h_remap_normalized_l2();
    expected[12] = closure_report.rho_remap_relative_conservation_defect();
    expected[13] = closure_report.rho_remap_normalized_l2();
    expected[14] = *closure_state.target_mass_kg;
    expected[15] = source.sample_field_value(4U, 0U);
    expected[16] = expected[15];
    for (std::size_t index = 0; index < expected.size(); ++index) {
      HUNDUN_CHECK(record.metrics[index].value.status ==
                   diagnostics::DiagnosticValueStatus::finite);
      HUNDUN_CHECK(record.metrics[index].value.bits ==
                   diagnostics::describe_fp64(expected[index]).bits);
    }
  } else if (record.level == diagnostics::DiagnosticLevel::invariants) {
    HUNDUN_CHECK(std::all_of(record.invariants.begin(),
                             record.invariants.end(),
                             [](const auto &item) { return item.passed; }));
  } else if (record.level == diagnostics::DiagnosticLevel::counters) {
    HUNDUN_CHECK(record.counters[0].value ==
                 closure_report.collective_count());
    HUNDUN_CHECK(record.counters[1].value ==
                 closure_report.evaluation_count());
    HUNDUN_CHECK(record.counters[2].value == closure_state.revision);
  } else {
    for (const auto &sample : record.samples) {
      const double expected =
          sample.field_id == "p0"
              ? closure_state.thermodynamic_pressure_pa
          : sample.field_id == "enthalpy"
              ? source.sample_field_value(0U, 0U)
          : sample.field_id == "rho"
              ? source.sample_field_value(2U, 0U)
          : sample.field_id == "rho_h"
              ? source.sample_field_value(3U, 0U)
              : source.sample_field_value(4U, 0U);
      HUNDUN_CHECK(sample.value.bits == diagnostics::describe_fp64(expected).bits);
    }
  }
}

void require_canonical_oracle_is_mutation_sensitive(
    const diagnostics::DiagnosticRecord &record) {
  const auto baseline = diagnostics::to_canonical_json(record);
  const auto prove = [&](auto mutation) {
    auto changed = record;
    mutation(changed);
    bool distinguished = false;
    try {
      distinguished = diagnostics::to_canonical_json(changed) != baseline;
    } catch (...) {
      distinguished = true;
    }
    HUNDUN_CHECK(distinguished);
  };
  prove([](auto &changed) { ++changed.schema_version; });
  prove([](auto &changed) { changed.module_id += ".mutated"; });
  prove([](auto &changed) { changed.instance_id += ".mutated"; });
  prove([](auto &changed) { ++changed.rank; });
  prove([](auto &changed) { ++changed.step; });
  prove([](auto &changed) { changed.time_s.bits ^= 1U; });
  prove([](auto &changed) { changed.phase += ".mutated"; });
  prove([](auto &changed) {
    changed.status = diagnostics::DiagnosticStatus::warning;
  });
  prove([](auto &changed) {
    changed.failure.classification =
        diagnostics::DiagnosticFailureClass::invalid_input;
  });
  prove([](auto &changed) { changed.failure.code += ".mutated"; });
  prove([](auto &changed) { ++changed.failure.lowest_failing_rank; });
  for (std::size_t index = 0; index < record.identities.size(); ++index)
    prove([&](auto &changed) {
      auto &identity = changed.identities[index];
      if (identity.revision)
        ++*identity.revision;
      else
        identity.layout_fingerprint->append(".mutated");
    });
  for (std::size_t index = 0; index < record.metrics.size(); ++index)
    prove([&](auto &changed) { changed.metrics[index].value.bits ^= 1U; });
  for (std::size_t index = 0; index < record.invariants.size(); ++index)
    prove([&](auto &changed) {
      changed.invariants[index].observed.bits ^= 1U;
    });
  for (std::size_t index = 0; index < record.counters.size(); ++index)
    prove([&](auto &changed) { ++changed.counters[index].value; });
  for (std::size_t index = 0; index < record.samples.size(); ++index)
    prove([&](auto &changed) { changed.samples[index].value.bits ^= 1U; });
}

void require_work(
    diagnostics::DiagnosticLevel level, std::size_t budget,
    diagnostics::DiagnosticScope scope,
    const hundun::flow::IdealGasClosureDiagnosticSource &source) {
  const auto work = Access::work(source);
  const auto cells = source.owned_cell_count();
  HUNDUN_CHECK(work.observations == 1U);
  HUNDUN_CHECK(work.fingerprint_items ==
               2U * cells + (source.relative_rank() == 0 ? 1U : 0U));
  HUNDUN_CHECK(work.summary_items ==
               (level == diagnostics::DiagnosticLevel::summary ? cells : 0U));
  HUNDUN_CHECK(work.invariant_items ==
               (level == diagnostics::DiagnosticLevel::invariants ? cells
                                                                  : 0U));
  HUNDUN_CHECK(work.sample_items == work.retained_sample_items);
  HUNDUN_CHECK(work.sample_items <= budget);
  HUNDUN_CHECK(work.full_field_copy_attempts == 0U);
  HUNDUN_CHECK(work.allocation_events ==
               (level == diagnostics::DiagnosticLevel::bounded_state_sample
                    ? 1U
                    : 0U));
  HUNDUN_CHECK(scope == diagnostics::DiagnosticScope::local
                   ? work.collective_calls == 0U
                   : work.collective_calls == 1U);
}

template <class Operation>
void require_error(Operation &&operation,
                   diagnostics::DiagnosticFailureClass classification,
                   std::string_view code, int rank, RecordingSink &sink,
                   std::size_t expected_calls = 0U,
                   std::size_t expected_records = 0U) {
  bool caught = false;
  try {
    operation();
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    caught = true;
    HUNDUN_CHECK(error.classification() == classification);
    HUNDUN_CHECK(error.code() == code);
    HUNDUN_CHECK(error.lowest_failing_rank() == rank);
  }
  HUNDUN_CHECK(caught);
  HUNDUN_CHECK(sink.calls == expected_calls);
  HUNDUN_CHECK(sink.records.size() == expected_records);
}

void run_task21_diagnostic_matrix(
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    const hundun::runtime::MpiContext &mpi,
    const hundun::flow::FlowState &state,
    const hundun::flow::IdealGasStepAttemptReport &report) {
  FaultReset reset(source);
  const auto before = capture(state, source, mpi);
  HUNDUN_CHECK(hundun::flow::test::IdealGasClosureTestAccess::
                   report_authenticated(report));

  for (const auto scope : {diagnostics::DiagnosticScope::local,
                           diagnostics::DiagnosticScope::collective}) {
    std::optional<std::string> fingerprint;
    for (std::uint8_t encoded = 0U; encoded < 4U; ++encoded) {
      const auto level = static_cast<diagnostics::DiagnosticLevel>(encoded);
      const auto request = request_for(source, level, scope);
      Access::reset(source);
      RecordingSink first;
      collect(source, mpi, request, first);
      HUNDUN_CHECK(first.calls == 1U && first.records.size() == 1U);
      const auto &record = first.records.front();
      HUNDUN_CHECK(record.schema_version ==
                   diagnostics::kDiagnosticRecordSchemaV1);
      HUNDUN_CHECK(record.module_kind ==
                   diagnostics::DiagnosticModuleKind::density_closure);
      HUNDUN_CHECK(record.module_id == "flow.ideal-gas-closure");
      HUNDUN_CHECK(record.instance_id == "primary");
      HUNDUN_CHECK(record.level == level && record.scope == scope);
      HUNDUN_CHECK(record.rank == mpi.rank());
      HUNDUN_CHECK(record.step == source.committed_step());
      HUNDUN_CHECK(record.phase == "ideal-gas-closure.attempt-result");
      HUNDUN_CHECK(record.status == diagnostics::DiagnosticStatus::ok);
      HUNDUN_CHECK(record.failure.classification ==
                   diagnostics::DiagnosticFailureClass::none);
      HUNDUN_CHECK(record.failure.code == "none");
      HUNDUN_CHECK(record.failure.lowest_failing_rank == -1);
      HUNDUN_CHECK(record.state_fingerprint.algorithm ==
                   diagnostics::kStateFingerprintAlgorithmV1);
      HUNDUN_CHECK(!record.state_fingerprint.hex.empty());
      if (!fingerprint)
        fingerprint = record.state_fingerprint.hex;
      HUNDUN_CHECK(record.state_fingerprint.hex == *fingerprint);
      require_payload(record);
      require_success_payload_values(record, source, mpi);
      diagnostics::validate(record, diagnostics::describe_diagnostics(source),
                            request);
      const auto canonical = diagnostics::to_canonical_json(record);
      HUNDUN_CHECK(!canonical.empty());
      require_canonical_oracle_is_mutation_sensitive(record);
      require_work(level, 7U, scope, source);
      require_pure(before, state, source, mpi);

      Access::reset(source);
      RecordingSink repeated;
      collect(source, mpi, request, repeated);
      HUNDUN_CHECK(repeated.calls == 1U && repeated.records.size() == 1U);
      HUNDUN_CHECK(diagnostics::to_canonical_json(repeated.records.front()) ==
                   canonical);
      require_work(level, 7U, scope, source);
      require_pure(before, state, source, mpi);
    }

    for (const std::size_t budget : {1U, 7U, 256U}) {
      auto sampled = request_for(
          source, diagnostics::DiagnosticLevel::bounded_state_sample, scope,
          budget, {"rho", "temperature"});
      Access::reset(source);
      RecordingSink sink;
      collect(source, mpi, sampled, sink);
      HUNDUN_CHECK(sink.calls == 1U && sink.records.size() == 1U);
      HUNDUN_CHECK(sink.records.front().sample_budget == budget);
      HUNDUN_CHECK(sink.records.front().samples.size() <= budget);
      for (const auto &sample : sink.records.front().samples)
        HUNDUN_CHECK(sample.field_id == "rho" ||
                     sample.field_id == "temperature");
      require_work(diagnostics::DiagnosticLevel::bounded_state_sample, budget,
                   scope, source);
      require_pure(before, state, source, mpi);
    }

    {
      const auto p0_only = request_for(
          source, diagnostics::DiagnosticLevel::bounded_state_sample, scope,
          1U, {"p0"});
      Access::reset(source);
      RecordingSink sink;
      collect(source, mpi, p0_only, sink);
      HUNDUN_CHECK(sink.calls == 1U && sink.records.size() == 1U);
      const auto &record = sink.records.front();
      const std::size_t expected =
          scope == diagnostics::DiagnosticScope::collective || mpi.rank() == 0
              ? 1U
              : 0U;
      HUNDUN_CHECK(record.eligible_sample_count == expected);
      HUNDUN_CHECK(record.samples.size() == expected);
      if (expected != 0U) {
        HUNDUN_CHECK(record.samples.front().field_id == "p0");
        HUNDUN_CHECK(record.samples.front().global_id == 0U);
        HUNDUN_CHECK(record.samples.front().value.bits ==
                     diagnostics::describe_fp64(
                         source.closure_state().thermodynamic_pressure_pa)
                         .bits);
      }
      require_work(diagnostics::DiagnosticLevel::bounded_state_sample, 1U,
                   scope, source);
      require_pure(before, state, source, mpi);
    }

  }

  if (!diagnostic_reference_mode.empty()) {
    const auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 256U);
    Access::reset(source);
    RecordingSink sink;
    collect(source, mpi, request, sink);
    HUNDUN_CHECK(sink.calls == 1U && sink.records.size() == 1U);
    auto reference_record = sink.records.front();
    reference_record.rank = 0;
    const auto canonical = diagnostics::to_canonical_json(reference_record);
    if (diagnostic_reference_mode == "--reference-write") {
      HUNDUN_CHECK(mpi.size() == 1);
      if (mpi.rank() == 0)
        write_diagnostic_reference(canonical);
    } else {
      HUNDUN_CHECK(diagnostic_reference_mode == "--reference-read");
      HUNDUN_CHECK(read_diagnostic_reference() == canonical);
    }
    require_pure(before, state, source, mpi);
  }

  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    ++request.frame.step;
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", -1, sink);
  }
  {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::local, 0U);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", -1, sink);
  }
  for (auto fields : {std::vector<std::string_view>{"unknown"},
                      std::vector<std::string_view>{"RHO"},
                      std::vector<std::string_view>{"rho", "rho"},
                      std::vector<std::string_view>{"temperature", "rho"}}) {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::local, 7U, std::move(fields));
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.selected-field", -1, sink);
  }
  {
    Access::set_fault(source, Fault::capability, mpi.rank());
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::capability,
        "closure.diagnostics.capability", -1, sink);
    Access::reset(source);
  }
  {
    Access::set_fault(source, Fault::record_validation, mpi.rank());
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "closure.diagnostics.record", -1, sink);
    Access::reset(source);
  }
  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::local);
    RecordingSink sink;
    sink.throw_on_submit = true;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::sink_failure,
        "diagnostics.sink.submit", -1, sink, 1U);
  }

  const int target = mpi.size() > 1 ? 1 : 0;
  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    if (mpi.rank() == target)
      ++request.frame.step;
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", target, sink);
  }
  {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 0U);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.frame", 0, sink);
  }
  {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 7U);
    if (mpi.rank() == target)
      request.selected_fields = {"unknown"};
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.selected-field", target, sink);
  }
  if (mpi.size() > 1) {
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 7U);
    if (mpi.rank() == target)
      request.selected_fields = {"rho"};
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_request,
        "closure.diagnostics.request-agreement", target, sink);
  }
  {
    Access::set_fault(source, Fault::provider_agreement, target);
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "closure.diagnostics.provider-agreement", target, sink);
    Access::reset(source);
  }
  for (const auto &[fault, classification, code] :
       {std::tuple{Fault::capability,
                   diagnostics::DiagnosticFailureClass::capability,
                   std::string_view("closure.diagnostics.capability")},
        std::tuple{Fault::request_preparation,
                   diagnostics::DiagnosticFailureClass::invalid_request,
                   std::string_view("closure.diagnostics.request-agreement")},
        std::tuple{Fault::aggregation,
                   diagnostics::DiagnosticFailureClass::layout,
                   std::string_view("closure.diagnostics.aggregation")},
        std::tuple{Fault::sample_preparation,
                   diagnostics::DiagnosticFailureClass::layout,
                   std::string_view("closure.diagnostics.sample-preparation")}}) {
    Access::set_fault(source, fault, target);
    auto request = request_for(
        source,
        fault == Fault::sample_preparation
            ? diagnostics::DiagnosticLevel::bounded_state_sample
            : diagnostics::DiagnosticLevel::summary,
        diagnostics::DiagnosticScope::collective, 7U);
    RecordingSink sink;
    require_error([&] { collect(source, mpi, request, sink); }, classification,
                  code, target, sink);
    Access::reset(source);
  }
  {
    Access::set_fault(source, Fault::raw_mpi, target);
    const auto request = request_for(
        source, diagnostics::DiagnosticLevel::summary,
        diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    bool typed = false;
    try {
      collect(source, mpi, request, sink);
    } catch (const hundun::runtime::MpiOperationError &error) {
      typed = std::string_view(error.what()).find("MPI_Task21") !=
              std::string_view::npos;
    }
    HUNDUN_CHECK(typed && sink.calls == 0U);
    Access::reset(source);
  }
  {
    Access::set_fault(source, Fault::oversized_agreement, target);
    const auto request = request_for(
        source, diagnostics::DiagnosticLevel::summary,
        diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::layout,
        "closure.diagnostics.aggregation", 0, sink);
    Access::reset(source);
  }
  {
    auto duplicate =
        hundun::runtime::MpiContext::duplicate(mpi.comm());
    const auto request = request_for(
        source, diagnostics::DiagnosticLevel::summary,
        diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    const auto &supplied = mpi.rank() == target ? duplicate : mpi;
    require_error(
        [&] { diagnostics::collect_diagnostics(source, supplied, request, sink); },
        diagnostics::DiagnosticFailureClass::layout,
        "closure.diagnostics.layout", target, sink);
  }
  {
    Access::set_fault(source, Fault::ownership_layout, target);
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::layout,
        "closure.diagnostics.layout", target, sink);
    Access::reset(source);
  }
  {
    Access::set_fault(source, Fault::sample_wire, target);
    auto request = request_for(
        source, diagnostics::DiagnosticLevel::bounded_state_sample,
        diagnostics::DiagnosticScope::collective, 7U);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::layout,
        "closure.diagnostics.sample-wire", target, sink);
    Access::reset(source);
  }
  {
    Access::set_fault(source, Fault::record_validation, target);
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::invalid_input,
        "closure.diagnostics.record", target, sink);
    Access::reset(source);
  }
  {
    auto request = request_for(source, diagnostics::DiagnosticLevel::summary,
                               diagnostics::DiagnosticScope::collective);
    RecordingSink sink;
    sink.throw_on_submit = mpi.rank() == target;
    require_error(
        [&] { collect(source, mpi, request, sink); },
        diagnostics::DiagnosticFailureClass::sink_failure,
        "diagnostics.sink.submit", target, sink, 1U,
        mpi.rank() == target ? 0U : 1U);
  }
  require_pure(before, state, source, mpi);
}

void run_task21_status_matrix(
    const hundun::flow::IdealGasClosureDiagnosticSource &source,
    const hundun::runtime::MpiContext &mpi,
    diagnostics::DiagnosticStatus expected_status,
    diagnostics::DiagnosticFailureClass expected_class,
    std::string_view expected_code, int expected_rank) {
  for (const auto scope : {diagnostics::DiagnosticScope::local,
                           diagnostics::DiagnosticScope::collective}) {
    for (std::uint8_t encoded = 0U; encoded < 4U; ++encoded) {
      const auto level = static_cast<diagnostics::DiagnosticLevel>(encoded);
      const auto request = request_for(source, level, scope, 1U);
      RecordingSink sink;
      collect(source, mpi, request, sink);
      HUNDUN_CHECK(sink.calls == 1U && sink.records.size() == 1U);
      const auto &record = sink.records.front();
      HUNDUN_CHECK(record.status == expected_status);
      HUNDUN_CHECK(record.failure.classification == expected_class);
      HUNDUN_CHECK(record.failure.code == expected_code);
      HUNDUN_CHECK(record.failure.lowest_failing_rank ==
                   (scope == diagnostics::DiagnosticScope::local
                        ? -1
                        : expected_status == diagnostics::DiagnosticStatus::warning
                              ? -1
                              : expected_rank));
      require_payload(record);
      if (level == diagnostics::DiagnosticLevel::summary) {
        const auto candidate = std::find_if(
            record.metrics.begin(), record.metrics.end(), [](const auto &item) {
              return item.id == "closure.candidate-pressure";
            });
        HUNDUN_CHECK(candidate != record.metrics.end());
        if (source.report().closure_report_available() &&
            source.report().closure_report().candidate_pressure_available()) {
          HUNDUN_CHECK(candidate->value.bits ==
                       diagnostics::describe_fp64(
                           source.report()
                               .closure_report()
                               .candidate_pressure_pa())
                           .bits);
        } else {
          HUNDUN_CHECK(candidate->value.status ==
                       diagnostics::DiagnosticValueStatus::unavailable);
        }
      } else if (level == diagnostics::DiagnosticLevel::counters) {
        HUNDUN_CHECK(record.counters.size() == 3U);
        const auto *closure_report =
            source.report().closure_report_available()
                ? &source.report().closure_report()
                : nullptr;
        HUNDUN_CHECK(record.counters[0].value ==
                     (closure_report == nullptr
                          ? 0U
                          : closure_report->collective_count()));
        HUNDUN_CHECK(record.counters[1].value ==
                     (closure_report == nullptr
                          ? 0U
                          : closure_report->evaluation_count()));
        HUNDUN_CHECK(record.counters[2].value ==
                     source.closure_state().revision);
      }
      diagnostics::validate(record, diagnostics::describe_diagnostics(source),
                            request);
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 1 || argc == 3);
    if (argc == 3) {
      diagnostic_reference_mode = argv[1];
      diagnostic_reference_path = argv[2];
      HUNDUN_CHECK(diagnostic_reference_mode == "--reference-write" ||
                   diagnostic_reference_mode == "--reference-read");
    }
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const auto descriptor =
        diagnostics::describe_ideal_gas_closure_diagnostics();
    HUNDUN_CHECK(descriptor.module_kind ==
                 diagnostics::DiagnosticModuleKind::density_closure);
    HUNDUN_CHECK(descriptor.module_id ==
                 std::string_view("flow.ideal-gas-closure"));
    HUNDUN_CHECK(descriptor.instance_id == std::string_view("primary"));
    const auto positive = Access::positive_invariant(
        "density.positive", "kg/m3", 1.0);
    HUNDUN_CHECK(positive.relation ==
                 diagnostics::InvariantRelation::positive);
    HUNDUN_CHECK(positive.limit.status ==
                 diagnostics::DiagnosticValueStatus::unavailable);
    HUNDUN_CHECK(positive.limit.bits == 0U);
    HUNDUN_CHECK(positive.passed);
    run(mpi, false, true, false);
  });
}
