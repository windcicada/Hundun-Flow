// SPDX-License-Identifier: Apache-2.0

#include "tests/support/diag_checkpoint_v2_test_access.hpp"
#include "src/flow_checkpoint_v2_detail.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/diag_checkpoint_v2.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_mpi_operation_error.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/checkpoint_v2_product_oracle.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace {

class OneRecordSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &value) override {
    record = value;
    ++calls;
  }
  hundun::diagnostics::DiagnosticRecord record;
  int calls{};
};

class DiagnosticErrorSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &) override {
    throw hundun::diagnostics::DiagnosticCollectionError(
        hundun::diagnostics::DiagnosticFailureClass::invalid_request,
        "test.sink.original", -1, "deliberate sink rejection");
  }
};

template <class Function> bool rejects(Function &&function) {
  try {
    function();
  } catch (const hundun::diagnostics::DiagnosticCollectionError &) {
    return true;
  }
  return false;
}

hundun::flow::CheckpointV2Report
report(int rank, hundun::flow::CheckpointV2CheckStatus rank_crc =
                     hundun::flow::CheckpointV2CheckStatus::passed) {
  hundun::flow::detail::CheckpointV2ReportValues values;
  values.operation = hundun::flow::CheckpointV2Operation::write;
  values.disposition = hundun::flow::CheckpointV2Disposition::completed;
  values.reason = hundun::flow::CheckpointV2FailureReason::none;
  values.phase = hundun::flow::CheckpointV2Phase::completed_marker;
  values.rank = rank;
  values.step = 4U;
  values.time_s = 0.125;
  values.local_logical_bytes = 100U + static_cast<std::uint64_t>(rank);
  values.local_actual_bytes = 200U + static_cast<std::uint64_t>(rank);
  values.global_logical_bytes = 404U;
  values.global_actual_bytes = 900U;
  values.local_crc64 = 0x12340000U + static_cast<std::uint64_t>(rank);
  values.manifest_crc64 = 0x5678U;
  values.file_count = 6U;
  values.crc_check_count = 6U;
  values.collective_count = 12U;
  values.rank_crc = rank_crc;
  values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::passed;
  values.fingerprint = hundun::flow::CheckpointV2CheckStatus::passed;
  values.partition = hundun::flow::CheckpointV2CheckStatus::passed;
  values.publication = hundun::flow::CheckpointV2CheckStatus::passed;
  return hundun::flow::detail::CheckpointV2Access::make(values);
}

hundun::diagnostics::DiagnosticRequest
request(int rank, hundun::diagnostics::DiagnosticLevel level,
        hundun::diagnostics::DiagnosticScope scope,
        std::string_view phase = "completed-marker") {
  return {level, scope, {rank, 4U, 0.125, phase}, {}, 0U};
}

std::string_view phase_name(hundun::flow::CheckpointV2Phase phase) {
  using Phase = hundun::flow::CheckpointV2Phase;
  switch (phase) {
  case Phase::none:
    return "none";
  case Phase::preflight:
    return "preflight";
  case Phase::transaction_entry:
    return "transaction-entry";
  case Phase::rank_payload:
    return "rank-payload";
  case Phase::rank_temporary_file:
    return "rank-temporary-file";
  case Phase::rank_publish:
    return "rank-publish";
  case Phase::manifest:
    return "manifest";
  case Phase::completed_marker:
    return "completed-marker";
  case Phase::marker_read:
    return "marker-read";
  case Phase::manifest_read:
    return "manifest-read";
  case Phase::rank_read:
    return "rank-read";
  case Phase::restore_prepare:
    return "restore-prepare";
  case Phase::restore_publish:
    return "restore-publish";
  }
  return "none";
}

hundun::flow::CheckpointV2Report
failure_report(int rank, hundun::flow::CheckpointV2FailureReason reason,
               hundun::flow::CheckpointV2Phase phase,
               int lowest_failing_rank = 0) {
  hundun::flow::detail::CheckpointV2ReportValues values;
  values.operation = static_cast<std::uint8_t>(phase) >=
                             static_cast<std::uint8_t>(
                                 hundun::flow::CheckpointV2Phase::marker_read)
                         ? hundun::flow::CheckpointV2Operation::read
                         : hundun::flow::CheckpointV2Operation::write;
  values.disposition = hundun::flow::CheckpointV2Disposition::failed;
  values.reason = reason;
  values.phase = phase;
  values.rank = rank;
  values.lowest_failing_rank = lowest_failing_rank;
  values.step = 4U;
  values.time_s = 0.125;
  values.local_logical_bytes = 10U + static_cast<std::uint64_t>(rank);
  values.local_actual_bytes = 20U + static_cast<std::uint64_t>(rank);
  values.local_crc64 =
      UINT64_C(0x1020304050607000) + static_cast<std::uint64_t>(rank);
  values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::failed;
  values.fingerprint = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.partition = hundun::flow::CheckpointV2CheckStatus::passed;
  values.transaction_entry = hundun::flow::CheckpointV2CheckStatus::passed;
  values.publication = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.rollback = hundun::flow::CheckpointV2CheckStatus::passed;
  values.rank_crc = rank == 0
                        ? hundun::flow::CheckpointV2CheckStatus::failed
                        : hundun::flow::CheckpointV2CheckStatus::not_checked;
  return hundun::flow::detail::CheckpointV2Access::make(values);
}

std::uint64_t independent_rank_crc_digest(int ranks) {
  std::vector<std::uint8_t> bytes;
  constexpr std::string_view domain =
      "hundun.checkpoint-v2.diagnostic-rank-crcs.v1";
  bytes.insert(bytes.end(), domain.begin(), domain.end());
  bytes.push_back(0U);
  const auto u32 = [&](std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      bytes.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  const auto i32 = [&](std::int32_t value) {
    std::uint32_t wire{};
    std::memcpy(&wire, &value, sizeof(wire));
    u32(wire);
  };
  const auto u64 = [&](std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      bytes.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  u32(1U);
  u32(static_cast<std::uint32_t>(ranks));
  for (int rank = 0; rank < ranks; ++rank) {
    i32(rank);
    u64(UINT64_C(0x1020304050607000) + static_cast<std::uint64_t>(rank));
  }
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  std::uint64_t crc{};
  for (const auto byte : bytes) {
    crc ^= static_cast<std::uint64_t>(byte) << 56U;
    for (unsigned bit = 0; bit < 8U; ++bit)
      crc = (crc & (UINT64_C(1) << 63U)) != 0U ? (crc << 1U) ^ polynomial
                                               : crc << 1U;
  }
  return crc;
}

hundun::runtime::FieldDescriptor
diagnostic_cell(const char *name, std::uint32_t components, bool conservative) {
  return {name,
          "1",
          "checkpoint-v2-diagnostics",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor diagnostic_face(const char *name,
                                                 std::uint32_t components) {
  auto result = diagnostic_cell(name, components, false);
  result.space = hundun::runtime::FunctionSpace::face_value;
  result.ghost_width = 0;
  return result;
}
hundun::config::FlowCaseConfig diagnostic_case(int ranks) {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "checkpoint-v2-diagnostics";
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = hundun::runtime::Int3{ranks, 1, 1};
  config.mesh.cells = {4, 2, 2};
  config.mesh.origin_m = {0.0, 0.0, 0.0};
  config.mesh.length_m = {1.0, 1.0, 1.0};
  config.mesh.mapping = hundun::config::MeshMapping::uniform_box;
  config.time = {hundun::config::TimeMode::fixed,
                 2,
                 0.01,
                 0.001,
                 0.1,
                 0.5,
                 0.25,
                 1.25,
                 0.5,
                 8};
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.dynamic_viscosity_pa_s = 0.01;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

struct ProductReports final {
  hundun::flow::CheckpointV2Report write;
  hundun::flow::CheckpointV2Report read;
  hundun::flow::CheckpointV2Report failure;
};

ProductReports product_failure_reports(const hundun::runtime::MpiContext &mpi) {
  auto config = diagnostic_case(mpi.size());
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, config.mesh.cells, {true, true, true},
      hundun::runtime::DecompositionOptions{
          hundun::runtime::Int3{mpi.size(), 1, 1}});
  hundun::mesh::MeshTopology topology(decomposition);
  hundun::mesh::MeshGeometry geometry(
      topology, hundun::mesh::UniformBoxMapping(config.mesh.origin_m,
                                                config.mesh.length_m));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(config, topology);
  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(diagnostic_cell("rho", 1U, true));
  fields.velocity = registry.declare_field(diagnostic_cell("u", 3U, false));
  fields.mechanical_pressure =
      registry.declare_field(diagnostic_cell("pi", 1U, false));
  fields.face_velocity = registry.declare_field(diagnostic_face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  registry.freeze();
  const hundun::flow::AcceptedStepMetadata metadata{
      0U, 0.0, config.time.initial_dt_s, 0.0,
      hundun::flow::MomentumTimeOrder::backward_euler};
  const auto make_state = [&] {
    return hundun::flow::FlowState::create(
        registry, {decomposition.local_extent(), topology.local_face_count()},
        fields, metadata);
  };
  hundun::flow::FlowLayerValues values;
  values.density.assign(topology.owned_cell_count(), 1.0);
  values.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
  values.mechanical_pressure.assign(topology.owned_cell_count(), -0.0);
  values.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  values.face_mass_flux.assign(topology.local_face_count(), 0.0);
  auto source = make_state();
  source.seed_accepted_layers(values, values);
  auto controller = hundun::flow::Bdf2RetryController::create(
      config.time, config.density_model, topology, geometry, mpi, source);
  const auto directory =
      std::filesystem::temp_directory_path() /
      ("hundun-task23-diagnostic-product-" + std::to_string(mpi.size()));
  if (mpi.rank() == 0)
    std::filesystem::remove_all(directory);
  mpi.barrier();
  auto written = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, source,
      controller.state(), std::nullopt, directory);
  auto destination = make_state();
  destination.seed_accepted_layers(values, values);
  auto restored = hundun::flow::read_checkpoint_v2(mpi, decomposition, topology,
                                                   geometry, boundaries, config,
                                                   destination, directory);
  auto duplicate = hundun::flow::write_checkpoint_v2(
      mpi, decomposition, topology, geometry, boundaries, config, source,
      controller.state(), std::nullopt, directory);
  HUNDUN_CHECK(written.disposition() ==
               hundun::flow::CheckpointV2Disposition::completed);
  HUNDUN_CHECK(restored.restored());
  HUNDUN_CHECK(duplicate.disposition() ==
               hundun::flow::CheckpointV2Disposition::failed);
  mpi.barrier();
  if (mpi.rank() == 0)
    std::filesystem::remove_all(directory);
  return {std::move(written), restored.report(), std::move(duplicate)};
}

double finite_value(const hundun::diagnostics::DiagnosticFp64 &value) {
  HUNDUN_CHECK(value.status ==
               hundun::diagnostics::DiagnosticValueStatus::finite);
  double result{};
  std::memcpy(&result, &value.bits, sizeof(result));
  return result;
}

hundun::diagnostics::DiagnosticStateFingerprint
expected_local_fingerprint(const hundun::flow::CheckpointV2Report &value) {
  constexpr std::array<std::string_view, 26> ids{
      "checkpoint.collective-count",
      "checkpoint.crc-check-count",
      "checkpoint.disposition",
      "checkpoint.exact-size-eof-status",
      "checkpoint.file-count",
      "checkpoint.fingerprint-status",
      "checkpoint.global-actual-bytes",
      "checkpoint.global-logical-bytes",
      "checkpoint.local-actual-bytes",
      "checkpoint.local-crc64",
      "checkpoint.local-logical-bytes",
      "checkpoint.lowest-failing-rank",
      "checkpoint.manifest-crc-status",
      "checkpoint.manifest-crc64",
      "checkpoint.operation",
      "checkpoint.partition-status",
      "checkpoint.phase",
      "checkpoint.publication-status",
      "checkpoint.rank",
      "checkpoint.rank-crc-status",
      "checkpoint.reason",
      "checkpoint.rollback-status",
      "checkpoint.semantic-fingerprint",
      "checkpoint.step",
      "checkpoint.time",
      "checkpoint.transaction-entry-status"};
  hundun::diagnostics::DiagnosticFingerprintAccumulator accumulator;
  const auto limb = [&](std::size_t index, std::uint64_t item) {
    accumulator.add(ids[index], static_cast<std::uint64_t>(value.rank()), 0U,
        hundun::diagnostics::describe_fp64(
            static_cast<double>(static_cast<std::uint32_t>(item))));
    accumulator.add(ids[index], static_cast<std::uint64_t>(value.rank()), 1U,
                    hundun::diagnostics::describe_fp64(static_cast<double>(
                        static_cast<std::uint32_t>(item >> 32U))));
  };
  const auto count = [&](std::size_t index, std::int64_t item) {
    accumulator.add(
        ids[index], static_cast<std::uint64_t>(value.rank()), 0U,
        hundun::diagnostics::describe_fp64(static_cast<double>(item)));
  };
  limb(0U, value.collective_count());
  limb(1U, value.crc_check_count());
  count(2U, static_cast<std::uint8_t>(value.disposition()));
  count(3U, static_cast<std::uint8_t>(value.exact_size_and_eof_status()));
  limb(4U, value.file_count());
  count(5U, static_cast<std::uint8_t>(value.fingerprint_status()));
  limb(6U, value.global_actual_bytes());
  limb(7U, value.global_logical_bytes());
  limb(8U, value.local_actual_bytes());
  limb(9U, value.local_crc64());
  limb(10U, value.local_logical_bytes());
  count(11U, value.lowest_failing_rank());
  count(12U, static_cast<std::uint8_t>(value.manifest_crc_status()));
  limb(13U, value.manifest_crc64());
  count(14U, static_cast<std::uint8_t>(value.operation()));
  count(15U, static_cast<std::uint8_t>(value.partition_status()));
  count(16U, static_cast<std::uint8_t>(value.phase()));
  count(17U, static_cast<std::uint8_t>(value.publication_status()));
  count(18U, value.rank());
  count(19U, static_cast<std::uint8_t>(value.rank_crc_status()));
  count(20U, static_cast<std::uint8_t>(value.reason()));
  count(21U, static_cast<std::uint8_t>(value.rollback_status()));
  limb(22U, value.semantic_fingerprint());
  limb(23U, value.step());
  accumulator.add(ids[24U], static_cast<std::uint64_t>(value.rank()), 0U,
                  hundun::diagnostics::describe_fp64(value.time_s()));
  count(25U, static_cast<std::uint8_t>(value.transaction_entry_status()));
  return accumulator.finish();
}

hundun::diagnostics::DiagnosticStateFingerprint
expected_collective_fingerprint(const hundun::flow::CheckpointV2Report &value,
                                const hundun::runtime::MpiContext &mpi) {
  const auto local = expected_local_fingerprint(value);
  HUNDUN_CHECK(local.hex.size() == 32U);
  const auto decode_limb = [&](std::size_t begin) {
    std::uint64_t result{};
    for (std::size_t index = begin; index < begin + 16U; ++index) {
      const char digit = local.hex[index];
      const auto nibble = digit >= '0' && digit <= '9'
                              ? static_cast<unsigned>(digit - '0')
                              : static_cast<unsigned>(digit - 'a' + 10);
      HUNDUN_CHECK(nibble < 16U);
      result = (result << 4U) | nibble;
    }
    return result;
  };
  const std::array<std::uint64_t, 2> parts{decode_limb(0U), decode_limb(16U)};
  std::vector<std::uint64_t> gathered(static_cast<std::size_t>(mpi.size()) *
                                      parts.size());
  hundun::runtime::check_mpi_result(
      MPI_Allgather(parts.data(), static_cast<int>(parts.size()), MPI_UINT64_T,
                    gathered.data(), static_cast<int>(parts.size()),
                    MPI_UINT64_T, mpi.comm()),
      "MPI_Allgather(Task23 expected diagnostic fingerprint)");
  hundun::diagnostics::DiagnosticFingerprintParts combined{};
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const auto offset = static_cast<std::size_t>(rank) * parts.size();
    combined.xor64 ^= gathered[offset];
    combined.sum64 += gathered[offset + 1U];
  }
  hundun::diagnostics::DiagnosticFingerprintAccumulator accumulator;
  accumulator.combine(combined);
  return accumulator.finish();
}

void require_fp64(double expected,
                  const hundun::diagnostics::DiagnosticFp64 &actual) {
  HUNDUN_CHECK(actual.status ==
               hundun::diagnostics::DiagnosticValueStatus::finite);
  std::uint64_t expected_bits{};
  std::memcpy(&expected_bits, &expected, sizeof(expected_bits));
  HUNDUN_CHECK(actual.bits == expected_bits);
}

hundun::diagnostics::DiagnosticFailureClass
expected_failure_class(hundun::flow::CheckpointV2FailureReason reason) {
  using Reason = hundun::flow::CheckpointV2FailureReason;
  switch (reason) {
  case Reason::none:
    return hundun::diagnostics::DiagnosticFailureClass::none;
  case Reason::invalid_input:
  case Reason::state:
    return hundun::diagnostics::DiagnosticFailureClass::invalid_input;
  case Reason::layout:
    return hundun::diagnostics::DiagnosticFailureClass::layout;
  case Reason::file_integrity:
  case Reason::filesystem:
    return hundun::diagnostics::DiagnosticFailureClass::file_integrity;
  }
  return hundun::diagnostics::DiagnosticFailureClass::invalid_input;
}

std::string
expected_failure_code(const hundun::flow::CheckpointV2Report &report) {
  using Reason = hundun::flow::CheckpointV2FailureReason;
  if (report.reason() == Reason::none)
    return "none";
  const auto operation =
      report.operation() == hundun::flow::CheckpointV2Operation::write ? "write"
                                                                       : "read";
  const auto reason = [&] {
    switch (report.reason()) {
    case Reason::none:
      return "none";
    case Reason::invalid_input:
      return "invalid-input";
    case Reason::layout:
      return "layout";
    case Reason::state:
      return "state";
    case Reason::file_integrity:
      return "file-integrity";
    case Reason::filesystem:
      return "filesystem";
    }
    return "invalid-input";
  }();
  return std::string("checkpoint-v2.") + operation + "." +
         std::string(phase_name(report.phase())) + "." + reason;
}

void require_exact_product_diagnostic_record(
    const hundun::flow::CheckpointV2Report &report,
    hundun::diagnostics::DiagnosticLevel level,
    const hundun::diagnostics::DiagnosticRecord &record) {
  hundun::test::checkpoint_v2_oracle::require_exact_product_diagnostic_record(
      report, level, record);
  using Check = hundun::flow::CheckpointV2CheckStatus;
  HUNDUN_CHECK(record.schema_version == 1U);
  HUNDUN_CHECK(record.module_kind ==
               hundun::diagnostics::DiagnosticModuleKind::checkpoint);
  HUNDUN_CHECK(record.module_id == "checkpoint-v2");
  HUNDUN_CHECK(record.instance_id == "checkpoint-v2");
  HUNDUN_CHECK(record.level == level);
  HUNDUN_CHECK(record.scope == hundun::diagnostics::DiagnosticScope::local);
  HUNDUN_CHECK(record.rank == report.rank());
  HUNDUN_CHECK(record.step == report.step());
  require_fp64(report.time_s(), record.time_s);
  HUNDUN_CHECK(record.phase == phase_name(report.phase()));
  HUNDUN_CHECK(
      record.status ==
      (report.disposition() == hundun::flow::CheckpointV2Disposition::completed
           ? hundun::diagnostics::DiagnosticStatus::ok
           : hundun::diagnostics::DiagnosticStatus::failed));
  HUNDUN_CHECK(record.failure.classification ==
               expected_failure_class(report.reason()));
  HUNDUN_CHECK(record.failure.code == expected_failure_code(report));
  HUNDUN_CHECK(record.failure.lowest_failing_rank == -1);
  HUNDUN_CHECK(record.identities.size() == 1U);
  HUNDUN_CHECK(record.identities[0].subject_id == "checkpoint-v2");
  const auto expected_identity = [&] {
    std::array<char, 17> value{};
    std::snprintf(
        value.data(), value.size(), "%016llx",
        static_cast<unsigned long long>(report.semantic_fingerprint()));
    return std::string(value.data());
  }();
  HUNDUN_CHECK(record.identities[0].layout_fingerprint.has_value());
  HUNDUN_CHECK(*record.identities[0].layout_fingerprint == expected_identity);
  HUNDUN_CHECK(!record.identities[0].revision);
  HUNDUN_CHECK(!record.identities[0].generation);
  HUNDUN_CHECK(!record.identities[0].allocation_identity);
  HUNDUN_CHECK(record.state_fingerprint.hex ==
               expected_local_fingerprint(report).hex);
  HUNDUN_CHECK(record.state_fingerprint.algorithm ==
               hundun::diagnostics::kStateFingerprintAlgorithmV1);
  HUNDUN_CHECK(record.sample_budget == 0U);
  HUNDUN_CHECK(record.eligible_sample_count == 0U);
  HUNDUN_CHECK(!record.samples_truncated);
  HUNDUN_CHECK(record.samples.empty());

  if (level == hundun::diagnostics::DiagnosticLevel::summary) {
    struct Expected final {
      std::string_view id;
      std::string_view unit;
      double value;
    };
    const std::array expected{
        Expected{"checkpoint.disposition", "count",
                 static_cast<double>(
                     static_cast<std::uint8_t>(report.disposition()))},
        Expected{"checkpoint.global-actual-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_actual_bytes() >> 32U))},
        Expected{"checkpoint.global-actual-bytes-low", "byte",
                 static_cast<double>(
                     static_cast<std::uint32_t>(report.global_actual_bytes()))},
        Expected{"checkpoint.global-logical-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_logical_bytes() >> 32U))},
        Expected{"checkpoint.global-logical-bytes-low", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_logical_bytes()))},
        Expected{"checkpoint.local-actual-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.local_actual_bytes() >> 32U))},
        Expected{"checkpoint.local-actual-bytes-low", "byte",
                 static_cast<double>(
                     static_cast<std::uint32_t>(report.local_actual_bytes()))},
        Expected{"checkpoint.local-logical-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.local_logical_bytes() >> 32U))},
        Expected{"checkpoint.local-logical-bytes-low", "byte",
                 static_cast<double>(
                     static_cast<std::uint32_t>(report.local_logical_bytes()))},
        Expected{
            "checkpoint.operation", "count",
            static_cast<double>(static_cast<std::uint8_t>(report.operation()))},
        Expected{
            "checkpoint.reason", "count",
            static_cast<double>(static_cast<std::uint8_t>(report.reason()))}};
    HUNDUN_CHECK(record.metrics.size() == expected.size());
    HUNDUN_CHECK(record.invariants.empty() && record.counters.empty());
    for (std::size_t index = 0; index < expected.size(); ++index) {
      HUNDUN_CHECK(record.metrics[index].id == expected[index].id);
      HUNDUN_CHECK(record.metrics[index].unit == expected[index].unit);
      HUNDUN_CHECK(record.metrics[index].kind ==
                   hundun::diagnostics::DiagnosticMetricKind::state_summary);
      require_fp64(expected[index].value, record.metrics[index].value);
    }
  } else if (level == hundun::diagnostics::DiagnosticLevel::invariants) {
    struct Expected final {
      std::string_view id;
      Check status;
    };
    const std::array expected{
        Expected{"checkpoint.exact-size-eof",
                 report.exact_size_and_eof_status()},
        Expected{"checkpoint.fingerprint", report.fingerprint_status()},
        Expected{"checkpoint.manifest-crc", report.manifest_crc_status()},
        Expected{"checkpoint.partition", report.partition_status()},
        Expected{"checkpoint.publication", report.publication_status()},
        Expected{"checkpoint.rank-crc", report.rank_crc_status()},
        Expected{"checkpoint.rollback", report.rollback_status()},
        Expected{"checkpoint.transaction-entry",
                 report.transaction_entry_status()}};
    HUNDUN_CHECK(record.invariants.size() == expected.size());
    HUNDUN_CHECK(record.metrics.empty() && record.counters.empty());
    for (std::size_t index = 0; index < expected.size(); ++index) {
      const auto observed = expected[index].status == Check::passed   ? 1.0
                            : expected[index].status == Check::failed ? 0.0
                                                                      : -1.0;
      HUNDUN_CHECK(record.invariants[index].id == expected[index].id);
      HUNDUN_CHECK(record.invariants[index].unit == "1");
      HUNDUN_CHECK(record.invariants[index].relation ==
                   hundun::diagnostics::InvariantRelation::equal);
      require_fp64(observed, record.invariants[index].observed);
      require_fp64(1.0, record.invariants[index].limit);
      HUNDUN_CHECK(record.invariants[index].passed ==
                   (expected[index].status == Check::passed));
    }
  } else {
    struct Expected final {
      std::string_view id;
      std::string_view unit;
      std::uint64_t value;
    };
    const std::array expected{
        Expected{"checkpoint.collective-count", "count",
                 report.collective_count()},
        Expected{"checkpoint.crc-check-count", "count",
                 report.crc_check_count()},
        Expected{"checkpoint.file-count", "count", report.file_count()},
        Expected{"checkpoint.global-actual-bytes", "byte",
                 report.global_actual_bytes()},
        Expected{"checkpoint.global-logical-bytes", "byte",
                 report.global_logical_bytes()},
        Expected{"checkpoint.local-actual-bytes", "byte",
                 report.local_actual_bytes()},
        Expected{"checkpoint.local-crc64-high", "count",
                 report.local_crc64() >> 32U},
        Expected{"checkpoint.local-crc64-low", "count",
                 static_cast<std::uint32_t>(report.local_crc64())},
        Expected{"checkpoint.local-logical-bytes", "byte",
                 report.local_logical_bytes()},
        Expected{"checkpoint.manifest-crc64-high", "count",
                 report.manifest_crc64() >> 32U},
        Expected{"checkpoint.manifest-crc64-low", "count",
                 static_cast<std::uint32_t>(report.manifest_crc64())}};
    HUNDUN_CHECK(record.counters.size() == expected.size());
    HUNDUN_CHECK(record.metrics.empty() && record.invariants.empty());
    for (std::size_t index = 0; index < expected.size(); ++index) {
      HUNDUN_CHECK(record.counters[index].id == expected[index].id);
      HUNDUN_CHECK(record.counters[index].unit == expected[index].unit);
      HUNDUN_CHECK(record.counters[index].value == expected[index].value);
    }
  }
}

bool diagnostic_oracle_is_mutation_sensitive(
    const hundun::flow::CheckpointV2Report &report,
    hundun::diagnostics::DiagnosticLevel level,
    const hundun::diagnostics::DiagnosticRecord &record) {
  const auto rejected = [&](auto mutate) {
    auto changed = record;
    mutate(changed);
    try {
      require_exact_product_diagnostic_record(report, level, changed);
    } catch (const std::exception &) {
      return true;
    }
    return false;
  };
  const bool common =
      rejected([](auto &item) { ++item.step; }) &&
      rejected([](auto &item) { item.module_id += "-changed"; }) &&
      rejected([](auto &item) { item.phase += "-changed"; }) &&
      rejected([](auto &item) {
        item.failure.classification =
            hundun::diagnostics::DiagnosticFailureClass::unavailable;
      }) &&
      rejected([](auto &item) { item.failure.code += "-changed"; }) &&
      rejected([](auto &item) { ++item.failure.lowest_failing_rank; }) &&
      rejected([](auto &item) {
        item.state_fingerprint.hex.front() =
            item.state_fingerprint.hex.front() == '0' ? '1' : '0';
      }) &&
      rejected(
          [](auto &item) { item.identities.front().subject_id += "-changed"; });
  if (level == hundun::diagnostics::DiagnosticLevel::summary)
    return common && rejected([](auto &item) {
             item.metrics.front().value.bits ^= UINT64_C(1);
           });
  if (level == hundun::diagnostics::DiagnosticLevel::invariants)
    return common && rejected([](auto &item) {
             item.invariants.front().observed.status =
                 hundun::diagnostics::DiagnosticValueStatus::unavailable;
           }) &&
           rejected([](auto &item) {
             item.invariants.front().passed = !item.invariants.front().passed;
         });
  return common && rejected([](auto &item) { ++item.counters.front().value; });
}

void run(const hundun::runtime::MpiContext &mpi, bool acceptance) {
  const auto local_report = report(mpi.rank());
  const auto expected_fingerprint = expected_local_fingerprint(local_report);
  auto source = hundun::flow::checkpoint_v2_diagnostic_source(local_report);
  const auto descriptor = hundun::diagnostics::describe_diagnostics(source);
  HUNDUN_CHECK(descriptor.module_kind ==
               hundun::diagnostics::DiagnosticModuleKind::checkpoint);
  HUNDUN_CHECK(!hundun::diagnostics::has_capability(
      descriptor.capabilities,
      hundun::diagnostics::DiagnosticCapability::bounded_state_sample));
  const auto ids =
      hundun::diagnostics::diagnostic_fingerprint_field_ids(source);
  constexpr std::array<std::string_view, 26> expected_ids{
      "checkpoint.collective-count",
      "checkpoint.crc-check-count",
      "checkpoint.disposition",
      "checkpoint.exact-size-eof-status",
      "checkpoint.file-count",
      "checkpoint.fingerprint-status",
      "checkpoint.global-actual-bytes",
      "checkpoint.global-logical-bytes",
      "checkpoint.local-actual-bytes",
      "checkpoint.local-crc64",
      "checkpoint.local-logical-bytes",
      "checkpoint.lowest-failing-rank",
      "checkpoint.manifest-crc-status",
      "checkpoint.manifest-crc64",
      "checkpoint.operation",
      "checkpoint.partition-status",
      "checkpoint.phase",
      "checkpoint.publication-status",
      "checkpoint.rank",
      "checkpoint.rank-crc-status",
      "checkpoint.reason",
      "checkpoint.rollback-status",
      "checkpoint.semantic-fingerprint",
      "checkpoint.step",
      "checkpoint.time",
      "checkpoint.transaction-entry-status"};
  HUNDUN_CHECK(ids.size() == 26U);
  HUNDUN_CHECK(std::is_sorted(ids.begin(), ids.end()));
  HUNDUN_CHECK(std::equal(ids.begin(), ids.end(), expected_ids.begin()));

  for (const auto level : {hundun::diagnostics::DiagnosticLevel::summary,
           hundun::diagnostics::DiagnosticLevel::invariants,
           hundun::diagnostics::DiagnosticLevel::counters}) {
    OneRecordSink first;
    const auto local =
        request(mpi.rank(), level, hundun::diagnostics::DiagnosticScope::local);
    hundun::diagnostics::collect_diagnostics(source, local, first);
    HUNDUN_CHECK(first.calls == 1);
    require_exact_product_diagnostic_record(local_report, level, first.record);
    HUNDUN_CHECK(diagnostic_oracle_is_mutation_sensitive(local_report, level,
                                                         first.record));
    HUNDUN_CHECK(first.record.scope ==
                 hundun::diagnostics::DiagnosticScope::local);
    HUNDUN_CHECK(first.record.state_fingerprint.hex ==
                 expected_fingerprint.hex);
    if (level == hundun::diagnostics::DiagnosticLevel::summary) {
      constexpr std::array<std::string_view, 11> expected{
          "checkpoint.disposition",
          "checkpoint.global-actual-bytes-high",
          "checkpoint.global-actual-bytes-low",
          "checkpoint.global-logical-bytes-high",
          "checkpoint.global-logical-bytes-low",
          "checkpoint.local-actual-bytes-high",
          "checkpoint.local-actual-bytes-low",
          "checkpoint.local-logical-bytes-high",
          "checkpoint.local-logical-bytes-low",
          "checkpoint.operation",
          "checkpoint.reason"};
      HUNDUN_CHECK(first.record.metrics.size() == expected.size());
      for (std::size_t index = 0; index < expected.size(); ++index)
        HUNDUN_CHECK(first.record.metrics[index].id == expected[index]);
    } else if (level == hundun::diagnostics::DiagnosticLevel::invariants) {
      constexpr std::array<std::string_view, 8> expected{
          "checkpoint.exact-size-eof", "checkpoint.fingerprint",
          "checkpoint.manifest-crc",   "checkpoint.partition",
          "checkpoint.publication",    "checkpoint.rank-crc",
          "checkpoint.rollback",       "checkpoint.transaction-entry"};
      HUNDUN_CHECK(first.record.invariants.size() == expected.size());
      for (std::size_t index = 0; index < expected.size(); ++index)
        HUNDUN_CHECK(first.record.invariants[index].id == expected[index]);
    } else {
      constexpr std::array<std::string_view, 11> expected{
          "checkpoint.collective-count",     "checkpoint.crc-check-count",
          "checkpoint.file-count",           "checkpoint.global-actual-bytes",
          "checkpoint.global-logical-bytes", "checkpoint.local-actual-bytes",
          "checkpoint.local-crc64-high",     "checkpoint.local-crc64-low",
          "checkpoint.local-logical-bytes",  "checkpoint.manifest-crc64-high",
          "checkpoint.manifest-crc64-low"};
      HUNDUN_CHECK(first.record.counters.size() == expected.size());
      for (std::size_t index = 0; index < expected.size(); ++index)
        HUNDUN_CHECK(first.record.counters[index].id == expected[index]);
    }
    const auto first_json =
        hundun::diagnostics::to_canonical_json(first.record);
    OneRecordSink second;
    hundun::diagnostics::collect_diagnostics(source, local, second);
    HUNDUN_CHECK(first_json ==
                 hundun::diagnostics::to_canonical_json(second.record));

    OneRecordSink collective;
    const auto global = request(
        mpi.rank(), level, hundun::diagnostics::DiagnosticScope::collective);
    hundun::diagnostics::collect_diagnostics(source, mpi, global, collective);
    HUNDUN_CHECK(collective.calls == 1);
    HUNDUN_CHECK(collective.record.scope ==
                 hundun::diagnostics::DiagnosticScope::collective);
    const auto expected_logical =
        static_cast<std::uint64_t>(mpi.size()) * 100U +
        static_cast<std::uint64_t>(mpi.size() * (mpi.size() - 1) / 2);
    const auto expected_actual =
        static_cast<std::uint64_t>(mpi.size()) * 200U +
        static_cast<std::uint64_t>(mpi.size() * (mpi.size() - 1) / 2);
    if (level == hundun::diagnostics::DiagnosticLevel::summary) {
      const auto find_metric = [&](std::string_view id) {
        return std::find_if(collective.record.metrics.begin(),
                            collective.record.metrics.end(),
                            [&](const auto &item) { return item.id == id; });
      };
      const auto logical = find_metric("checkpoint.local-logical-bytes-low");
      const auto actual = find_metric("checkpoint.local-actual-bytes-low");
      HUNDUN_CHECK(logical != collective.record.metrics.end());
      HUNDUN_CHECK(actual != collective.record.metrics.end());
      HUNDUN_CHECK(finite_value(logical->value) ==
                   static_cast<double>(expected_logical));
      HUNDUN_CHECK(finite_value(actual->value) ==
                   static_cast<double>(expected_actual));
    } else if (level == hundun::diagnostics::DiagnosticLevel::invariants) {
      const auto rank_crc = std::find_if(
          collective.record.invariants.begin(),
          collective.record.invariants.end(),
          [](const auto &item) { return item.id == "checkpoint.rank-crc"; });
      HUNDUN_CHECK(rank_crc != collective.record.invariants.end());
      HUNDUN_CHECK(rank_crc->passed);
    } else {
      const auto find_counter = [&](std::string_view id) {
        return std::find_if(collective.record.counters.begin(),
                            collective.record.counters.end(),
                            [&](const auto &item) { return item.id == id; });
      };
      const auto logical = find_counter("checkpoint.local-logical-bytes");
      const auto actual = find_counter("checkpoint.local-actual-bytes");
      HUNDUN_CHECK(logical != collective.record.counters.end());
      HUNDUN_CHECK(actual != collective.record.counters.end());
      HUNDUN_CHECK(logical->value == expected_logical);
      HUNDUN_CHECK(actual->value == expected_actual);
    }
    std::array<char, 33> root{};
    if (mpi.rank() == 0)
      std::copy(collective.record.state_fingerprint.hex.begin(),
                collective.record.state_fingerprint.hex.end(), root.begin());
    hundun::runtime::check_mpi_result(
        MPI_Bcast(root.data(), static_cast<int>(root.size()), MPI_CHAR, 0,
                  mpi.comm()),
        "MPI_Bcast(Checkpoint diagnostics fingerprint)");
    HUNDUN_CHECK(collective.record.state_fingerprint.hex ==
                 std::string(root.data()));
  }
  for (int capability_case = 0; capability_case < 2; ++capability_case) {
    auto unsupported =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local);
    if (capability_case == 0)
      unsupported.selected_fields = {"rho"};
    else
      unsupported.sample_budget = 1U;
    OneRecordSink rejected_sink;
    bool exact = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, unsupported,
                                               rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      exact = error.classification() ==
                  hundun::diagnostics::DiagnosticFailureClass::capability &&
              error.code() == "checkpoint-v2.diagnostics.capability";
    }
    HUNDUN_CHECK(exact);
    HUNDUN_CHECK(rejected_sink.calls == 0);
  }
  if (!acceptance)
    return;

  using DiagnosticAccess =
      hundun::diagnostics::test::CheckpointV2DiagnosticsTestAccess;
  using DiagnosticFault =
      hundun::diagnostics::test::CheckpointV2DiagnosticFault;
  if (mpi.size() > 1) {
    for (int stale_rank = 0; stale_rank < mpi.size(); ++stale_rank) {
      auto stale_source =
          hundun::flow::checkpoint_v2_diagnostic_source(local_report);
      std::unique_ptr<hundun::flow::CheckpointV2DiagnosticSource> retained;
      if (mpi.rank() == stale_rank)
        retained = std::make_unique<hundun::flow::CheckpointV2DiagnosticSource>(
            std::move(stale_source));
      OneRecordSink rejected_sink;
      bool exact = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            stale_source, mpi,
            request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                    hundun::diagnostics::DiagnosticScope::collective),
            rejected_sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        exact =
            error.classification() ==
                hundun::diagnostics::DiagnosticFailureClass::invalid_input &&
            error.code() == "checkpoint-v2.diagnostics.stale-source" &&
            error.lowest_failing_rank() == stale_rank;
      }
      HUNDUN_CHECK(exact);
      HUNDUN_CHECK(rejected_sink.calls == 0);
    }
  }

  enum class CollectiveRequestCase : std::uint8_t {
    wrong_scope,
    selected_field,
    sample_budget,
    bounded_level,
    invalid_level,
    invalid_scope,
    frame_rank,
    frame_step,
    frame_time,
    frame_phase
  };
  constexpr std::array collective_request_cases{
      CollectiveRequestCase::wrong_scope,
      CollectiveRequestCase::selected_field,
      CollectiveRequestCase::sample_budget,
      CollectiveRequestCase::bounded_level,
      CollectiveRequestCase::invalid_level,
      CollectiveRequestCase::invalid_scope,
      CollectiveRequestCase::frame_rank,
      CollectiveRequestCase::frame_step,
      CollectiveRequestCase::frame_time,
      CollectiveRequestCase::frame_phase};
  const auto mutate_collective_request =
      [](auto &candidate, CollectiveRequestCase request_case) {
        using Level = hundun::diagnostics::DiagnosticLevel;
        using Scope = hundun::diagnostics::DiagnosticScope;
        switch (request_case) {
        case CollectiveRequestCase::wrong_scope:
          candidate.scope = Scope::local;
          break;
        case CollectiveRequestCase::selected_field:
          candidate.selected_fields = {"rho"};
          break;
        case CollectiveRequestCase::sample_budget:
          candidate.sample_budget = 1U;
          break;
        case CollectiveRequestCase::bounded_level:
          candidate.level = Level::bounded_state_sample;
          break;
        case CollectiveRequestCase::invalid_level:
          candidate.level = static_cast<Level>(255U);
          break;
        case CollectiveRequestCase::invalid_scope:
          candidate.scope = static_cast<Scope>(255U);
          break;
        case CollectiveRequestCase::frame_rank:
          ++candidate.frame.rank;
          break;
        case CollectiveRequestCase::frame_step:
          ++candidate.frame.step;
          break;
        case CollectiveRequestCase::frame_time:
          candidate.frame.time_s = std::nextafter(
              candidate.frame.time_s, std::numeric_limits<double>::infinity());
          break;
        case CollectiveRequestCase::frame_phase:
          candidate.frame.phase = "wrong-phase";
          break;
        }
      };
  const auto is_capability_case = [](CollectiveRequestCase request_case) {
    return request_case == CollectiveRequestCase::wrong_scope ||
           request_case == CollectiveRequestCase::selected_field ||
           request_case == CollectiveRequestCase::sample_budget ||
           request_case == CollectiveRequestCase::bounded_level;
  };
  for (const auto request_case : collective_request_cases) {
    auto candidate =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    mutate_collective_request(candidate, request_case);
    OneRecordSink rejected_sink;
    bool exact = false;
    try {
      hundun::diagnostics::collect_diagnostics(source, mpi, candidate,
                                               rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      exact = error.classification() ==
                  (is_capability_case(request_case)
                       ? hundun::diagnostics::DiagnosticFailureClass::capability
                       : hundun::diagnostics::DiagnosticFailureClass::
                             invalid_request) &&
              error.code() == (is_capability_case(request_case)
                                   ? "checkpoint-v2.diagnostics.capability"
                                   : "checkpoint-v2.diagnostics.frame") &&
              error.lowest_failing_rank() == -1;
    }
    HUNDUN_CHECK(exact);
    HUNDUN_CHECK(rejected_sink.calls == 0);
  }
  if (mpi.size() > 1) {
    for (int failing_rank = 0; failing_rank < mpi.size(); ++failing_rank)
      for (const auto request_case : collective_request_cases) {
        auto candidate =
            request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                    hundun::diagnostics::DiagnosticScope::collective);
        if (mpi.rank() == failing_rank)
          mutate_collective_request(candidate, request_case);
        OneRecordSink rejected_sink;
        bool exact = false;
        try {
          hundun::diagnostics::collect_diagnostics(source, mpi, candidate,
                                                   rejected_sink);
        } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
          exact = error.classification() ==
                      hundun::diagnostics::DiagnosticFailureClass::
                          collective_operation &&
                  error.code() ==
                      "checkpoint-v2.diagnostics.collective-agreement" &&
                  error.lowest_failing_rank() == failing_rank;
        }
        HUNDUN_CHECK(exact);
        HUNDUN_CHECK(rejected_sink.calls == 0);
      }
  }

  DiagnosticAccess::reset();
  DiagnosticAccess::set_fault(DiagnosticFault::raw_mpi);
  {
    OneRecordSink local_while_collective_armed;
    hundun::diagnostics::collect_diagnostics(
        source,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local),
        local_while_collective_armed);
    HUNDUN_CHECK(local_while_collective_armed.calls == 1);
    HUNDUN_CHECK(DiagnosticAccess::work().collective_calls == 0U);
  }
  bool typed_mpi_error = false;
  try {
    OneRecordSink collective_fault;
    hundun::diagnostics::collect_diagnostics(
        source, mpi,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective),
        collective_fault);
  } catch (const hundun::runtime::MpiOperationError &) {
    typed_mpi_error = true;
  }
  HUNDUN_CHECK(typed_mpi_error);
  HUNDUN_CHECK(DiagnosticAccess::work().collective_calls == 1U);
  DiagnosticAccess::reset();

  for (const auto preparation_fault :
       {DiagnosticFault::request_preparation,
        DiagnosticFault::local_record_preparation,
        DiagnosticFault::post_gather_preparation,
        DiagnosticFault::final_record_preparation}) {
    for (int failure_rank = 0; failure_rank < mpi.size(); ++failure_rank) {
      DiagnosticAccess::reset();
      DiagnosticAccess::set_fault(preparation_fault, failure_rank);
      OneRecordSink rejected_sink;
      bool exact_failure = false;
      try {
        hundun::diagnostics::collect_diagnostics(
            source, mpi,
            request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                    hundun::diagnostics::DiagnosticScope::collective),
            rejected_sink);
      } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
        exact_failure =
            error.classification() ==
                hundun::diagnostics::DiagnosticFailureClass::
                    collective_operation &&
            error.code() == "checkpoint-v2.diagnostics.preparation" &&
            error.lowest_failing_rank() == failure_rank;
      }
      HUNDUN_CHECK(exact_failure);
      HUNDUN_CHECK(rejected_sink.calls == 0);
    }
  }
  DiagnosticAccess::reset();

  auto actual = product_failure_reports(mpi);
  const std::array actual_reports{actual.write, actual.read, actual.failure};
  for (const auto &actual_report : actual_reports) {
    auto actual_source =
        hundun::flow::checkpoint_v2_diagnostic_source(actual_report);
    for (const auto level : {hundun::diagnostics::DiagnosticLevel::summary,
                             hundun::diagnostics::DiagnosticLevel::invariants,
                             hundun::diagnostics::DiagnosticLevel::counters}) {
      const hundun::diagnostics::DiagnosticRequest local{
          level,
          hundun::diagnostics::DiagnosticScope::local,
          {actual_report.rank(), actual_report.step(), actual_report.time_s(),
           phase_name(actual_report.phase())},
          {},
          0U};
      OneRecordSink first;
      hundun::diagnostics::collect_diagnostics(actual_source, local, first);
      HUNDUN_CHECK(first.calls == 1);
      require_exact_product_diagnostic_record(actual_report, level,
                                              first.record);
      HUNDUN_CHECK(diagnostic_oracle_is_mutation_sensitive(actual_report, level,
                                                           first.record));
      HUNDUN_CHECK(first.record.rank == actual_report.rank());
      HUNDUN_CHECK(first.record.step == actual_report.step());
      HUNDUN_CHECK(first.record.phase == phase_name(actual_report.phase()));
      HUNDUN_CHECK(first.record.status ==
                   (actual_report.disposition() ==
                            hundun::flow::CheckpointV2Disposition::completed
                        ? hundun::diagnostics::DiagnosticStatus::ok
                        : hundun::diagnostics::DiagnosticStatus::failed));
      HUNDUN_CHECK(first.record.state_fingerprint.hex ==
                   expected_local_fingerprint(actual_report).hex);
      OneRecordSink repeated;
      hundun::diagnostics::collect_diagnostics(actual_source, local, repeated);
      HUNDUN_CHECK(hundun::diagnostics::to_canonical_json(first.record) ==
                   hundun::diagnostics::to_canonical_json(repeated.record));

      auto collective_request = local;
      collective_request.scope =
          hundun::diagnostics::DiagnosticScope::collective;
      OneRecordSink collective_record;
      hundun::diagnostics::collect_diagnostics(
          actual_source, mpi, collective_request, collective_record);
      HUNDUN_CHECK(collective_record.calls == 1);
      HUNDUN_CHECK(collective_record.record.failure.lowest_failing_rank ==
                   actual_report.lowest_failing_rank());
    }
  }

  auto lifetime_source =
      hundun::flow::checkpoint_v2_diagnostic_source(actual.write);
  const auto lifetime_request = hundun::diagnostics::DiagnosticRequest{
      hundun::diagnostics::DiagnosticLevel::summary,
      hundun::diagnostics::DiagnosticScope::local,
      {actual.write.rank(), actual.write.step(), actual.write.time_s(),
       phase_name(actual.write.phase())},
      {},
      0U};
  OneRecordSink lifetime_before;
  hundun::diagnostics::collect_diagnostics(lifetime_source, lifetime_request,
                                           lifetime_before);
  actual = {failure_report(mpi.rank(),
                           hundun::flow::CheckpointV2FailureReason::state,
                           hundun::flow::CheckpointV2Phase::restore_prepare),
            failure_report(mpi.rank(),
                           hundun::flow::CheckpointV2FailureReason::layout,
                           hundun::flow::CheckpointV2Phase::manifest_read),
            report(mpi.rank())};
  OneRecordSink lifetime_after;
  hundun::diagnostics::collect_diagnostics(lifetime_source, lifetime_request,
                                           lifetime_after);
  HUNDUN_CHECK(hundun::diagnostics::to_canonical_json(lifetime_before.record) ==
               hundun::diagnostics::to_canonical_json(lifetime_after.record));

  if (mpi.size() > 1) {
    auto disagreement =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    if (mpi.rank() == mpi.size() - 1)
      ++disagreement.frame.step;
    bool rejected = false;
    try {
      OneRecordSink rejected_sink;
      hundun::diagnostics::collect_diagnostics(source, mpi, disagreement,
                                               rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = error.classification() ==
              hundun::diagnostics::DiagnosticFailureClass::
                  collective_operation &&
          error.lowest_failing_rank() == mpi.size() - 1;
    }
    HUNDUN_CHECK(rejected);

    auto capability_disagreement =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    if (mpi.rank() == mpi.size() - 1)
      capability_disagreement.selected_fields = {"rho"};
    bool capability_rejected = false;
    try {
      OneRecordSink rejected_sink;
      hundun::diagnostics::collect_diagnostics(
          source, mpi, capability_disagreement, rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      capability_rejected =
          error.classification() ==
                  hundun::diagnostics::DiagnosticFailureClass::
                      collective_operation &&
          error.code() == "checkpoint-v2.diagnostics.collective-agreement" &&
          error.lowest_failing_rank() == mpi.size() - 1;
    }
    HUNDUN_CHECK(capability_rejected);
  }

  auto moved = std::move(source);
  HUNDUN_CHECK(rejects([&] {
    OneRecordSink sink;
    hundun::diagnostics::collect_diagnostics(
        source,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local),
        sink);
  }));
  OneRecordSink sink;
  hundun::diagnostics::collect_diagnostics(
      moved,
      request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
              hundun::diagnostics::DiagnosticScope::local),
      sink);
  HUNDUN_CHECK(sink.calls == 1);

  const auto require_request_error =
      [&](const hundun::diagnostics::DiagnosticRequest &candidate,
          hundun::diagnostics::DiagnosticFailureClass classification,
          std::string_view code) {
        OneRecordSink rejected_sink;
        bool exact = false;
        try {
          hundun::diagnostics::collect_diagnostics(moved, candidate,
                                                   rejected_sink);
        } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
          exact = error.classification() == classification &&
                  error.code() == code && error.lowest_failing_rank() == -1;
        }
        HUNDUN_CHECK(exact);
        HUNDUN_CHECK(rejected_sink.calls == 0);
      };
  {
    auto wrong_scope =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    require_request_error(
        wrong_scope, hundun::diagnostics::DiagnosticFailureClass::capability,
        "checkpoint-v2.diagnostics.capability");
  }
  {
    auto selected =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local);
    selected.selected_fields = {"rho"};
    require_request_error(
        selected, hundun::diagnostics::DiagnosticFailureClass::capability,
        "checkpoint-v2.diagnostics.capability");
  }
  {
    auto sampled =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local);
    sampled.sample_budget = 1U;
    require_request_error(
        sampled, hundun::diagnostics::DiagnosticFailureClass::capability,
        "checkpoint-v2.diagnostics.capability");
  }
  for (int capability_case = 0; capability_case < 2; ++capability_case) {
    auto unsupported =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    if (capability_case == 0)
      unsupported.selected_fields = {"rho"};
    else
      unsupported.sample_budget = 1U;
    OneRecordSink rejected_sink;
    bool exact = false;
    try {
      hundun::diagnostics::collect_diagnostics(moved, mpi, unsupported,
                                               rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      exact = error.classification() ==
                  hundun::diagnostics::DiagnosticFailureClass::capability &&
              error.code() == "checkpoint-v2.diagnostics.capability" &&
              error.lowest_failing_rank() == -1;
    }
    HUNDUN_CHECK(exact);
    HUNDUN_CHECK(rejected_sink.calls == 0);
  }
  for (int frame_case = 0; frame_case < 4; ++frame_case) {
    auto invalid =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local);
    if (frame_case == 0)
      ++invalid.frame.rank;
    else if (frame_case == 1)
      ++invalid.frame.step;
    else if (frame_case == 2)
      invalid.frame.time_s = std::nextafter(
          invalid.frame.time_s, std::numeric_limits<double>::infinity());
    else
      invalid.frame.phase = "wrong-phase";
    require_request_error(
        invalid, hundun::diagnostics::DiagnosticFailureClass::invalid_request,
        "checkpoint-v2.diagnostics.frame");
  }
  {
    auto invalid_level =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local);
    invalid_level.level =
        static_cast<hundun::diagnostics::DiagnosticLevel>(255U);
    require_request_error(
        invalid_level,
        hundun::diagnostics::DiagnosticFailureClass::invalid_request,
        "checkpoint-v2.diagnostics.frame");
  }
  {
    auto invalid_scope =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local);
    invalid_scope.scope =
        static_cast<hundun::diagnostics::DiagnosticScope>(255U);
    require_request_error(
        invalid_scope,
        hundun::diagnostics::DiagnosticFailureClass::invalid_request,
        "checkpoint-v2.diagnostics.frame");
  }
  {
    DiagnosticErrorSink rejected_sink;
    try {
      hundun::diagnostics::collect_diagnostics(
          moved,
          request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                  hundun::diagnostics::DiagnosticScope::local),
          rejected_sink);
      HUNDUN_CHECK(false);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      HUNDUN_CHECK(error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::sink_failure);
      HUNDUN_CHECK(error.code() == "diagnostics.sink.submit");
    }
  }

  using Reason = hundun::flow::CheckpointV2FailureReason;
  using Phase = hundun::flow::CheckpointV2Phase;
  struct FailureCase final {
    Reason reason;
    Phase phase;
    hundun::diagnostics::DiagnosticFailureClass classification;
    std::string_view suffix;
  };
  constexpr std::array failure_cases{
      FailureCase{Reason::invalid_input, Phase::preflight,
                  hundun::diagnostics::DiagnosticFailureClass::invalid_input,
                  "write.preflight.invalid-input"},
      FailureCase{Reason::layout, Phase::preflight,
                  hundun::diagnostics::DiagnosticFailureClass::layout,
                  "write.preflight.layout"},
      FailureCase{Reason::state, Phase::restore_prepare,
                  hundun::diagnostics::DiagnosticFailureClass::invalid_input,
                  "read.restore-prepare.state"},
      FailureCase{Reason::file_integrity, Phase::rank_read,
                  hundun::diagnostics::DiagnosticFailureClass::file_integrity,
                  "read.rank-read.file-integrity"},
      FailureCase{Reason::filesystem, Phase::manifest,
                  hundun::diagnostics::DiagnosticFailureClass::file_integrity,
                  "write.manifest.filesystem"}};
  for (const auto &item : failure_cases) {
    const auto failed = failure_report(mpi.rank(), item.reason, item.phase);
    auto failed_source = hundun::flow::checkpoint_v2_diagnostic_source(failed);
    OneRecordSink local_failure;
    hundun::diagnostics::collect_diagnostics(
        failed_source,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::invariants,
                hundun::diagnostics::DiagnosticScope::local,
                phase_name(item.phase)),
        local_failure);
    HUNDUN_CHECK(local_failure.record.status ==
                 hundun::diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(local_failure.record.failure.classification ==
                 item.classification);
    HUNDUN_CHECK(local_failure.record.failure.code ==
                 std::string("checkpoint-v2.") + std::string(item.suffix));
    const auto exact = std::find_if(
        local_failure.record.invariants.begin(),
        local_failure.record.invariants.end(), [](const auto &invariant) {
          return invariant.id == "checkpoint.exact-size-eof";
        });
    const auto manifest = std::find_if(
        local_failure.record.invariants.begin(),
        local_failure.record.invariants.end(), [](const auto &invariant) {
          return invariant.id == "checkpoint.manifest-crc";
        });
    HUNDUN_CHECK(exact != local_failure.record.invariants.end());
    HUNDUN_CHECK(manifest != local_failure.record.invariants.end());
    HUNDUN_CHECK(finite_value(exact->observed) == 0.0);
    HUNDUN_CHECK(finite_value(manifest->observed) == -1.0);
  }

  constexpr std::array all_reasons{Reason::invalid_input, Reason::layout,
                                   Reason::state, Reason::file_integrity,
                                   Reason::filesystem};
  constexpr std::array all_phases{
      Phase::preflight,        Phase::transaction_entry,
      Phase::rank_payload,     Phase::rank_temporary_file,
      Phase::rank_publish,     Phase::manifest,
      Phase::completed_marker, Phase::marker_read,
      Phase::manifest_read,    Phase::rank_read,
      Phase::restore_prepare,  Phase::restore_publish};
  for (int lowest = 0; lowest < mpi.size(); ++lowest)
    for (const auto reason : all_reasons)
      for (const auto phase : all_phases) {
        const auto exhaustive =
            failure_report(mpi.rank(), reason, phase, lowest);
        auto exhaustive_source =
            hundun::flow::checkpoint_v2_diagnostic_source(exhaustive);
        const hundun::diagnostics::DiagnosticRequest exhaustive_request{
            hundun::diagnostics::DiagnosticLevel::summary,
            hundun::diagnostics::DiagnosticScope::collective,
            {mpi.rank(), exhaustive.step(), exhaustive.time_s(),
             phase_name(phase)},
            {},
            0U};
        OneRecordSink exhaustive_sink;
        hundun::diagnostics::collect_diagnostics(
            exhaustive_source, mpi, exhaustive_request, exhaustive_sink);
        HUNDUN_CHECK(exhaustive_sink.calls == 1);
        HUNDUN_CHECK(exhaustive_sink.record.status ==
                     hundun::diagnostics::DiagnosticStatus::failed);
        HUNDUN_CHECK(exhaustive_sink.record.phase == phase_name(phase));
        HUNDUN_CHECK(exhaustive_sink.record.failure.lowest_failing_rank ==
                     lowest);
        HUNDUN_CHECK(exhaustive_sink.record.state_fingerprint.hex ==
                     expected_collective_fingerprint(exhaustive, mpi).hex);
      }

  const auto failed =
      failure_report(mpi.rank(), Reason::file_integrity, Phase::rank_read);
  auto failed_source = hundun::flow::checkpoint_v2_diagnostic_source(failed);
  OneRecordSink failed_invariants;
  hundun::diagnostics::collect_diagnostics(
      failed_source, mpi,
      request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::invariants,
              hundun::diagnostics::DiagnosticScope::collective, "rank-read"),
      failed_invariants);
  const auto rank_crc = std::find_if(
      failed_invariants.record.invariants.begin(),
      failed_invariants.record.invariants.end(), [](const auto &invariant) {
        return invariant.id == "checkpoint.rank-crc";
      });
  HUNDUN_CHECK(rank_crc != failed_invariants.record.invariants.end());
  HUNDUN_CHECK(!rank_crc->passed);
  HUNDUN_CHECK(finite_value(rank_crc->observed) == 0.0);

  for (int matrix_case = 0; matrix_case < (mpi.size() > 1 ? 2 : 1);
       ++matrix_case) {
    const auto local_status =
        matrix_case == 0 || mpi.rank() == mpi.size() - 1
            ? hundun::flow::CheckpointV2CheckStatus::not_checked
            : hundun::flow::CheckpointV2CheckStatus::passed;
    const auto matrix_report = report(mpi.rank(), local_status);
    auto matrix_source =
        hundun::flow::checkpoint_v2_diagnostic_source(matrix_report);
    OneRecordSink matrix_sink;
    hundun::diagnostics::collect_diagnostics(
        matrix_source, mpi,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::invariants,
                hundun::diagnostics::DiagnosticScope::collective),
        matrix_sink);
    const auto matrix_rank_crc = std::find_if(
        matrix_sink.record.invariants.begin(),
        matrix_sink.record.invariants.end(), [](const auto &invariant) {
          return invariant.id == "checkpoint.rank-crc";
        });
    HUNDUN_CHECK(matrix_rank_crc != matrix_sink.record.invariants.end());
    HUNDUN_CHECK(!matrix_rank_crc->passed);
    HUNDUN_CHECK(finite_value(matrix_rank_crc->observed) == -1.0);
  }

  OneRecordSink failed_counters;
  hundun::diagnostics::collect_diagnostics(
      failed_source, mpi,
      request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::counters,
              hundun::diagnostics::DiagnosticScope::collective, "rank-read"),
      failed_counters);
  const auto counter = [&](std::string_view id) {
    return std::find_if(failed_counters.record.counters.begin(),
                        failed_counters.record.counters.end(),
                        [&](const auto &value) { return value.id == id; });
  };
  const auto high = counter("checkpoint.local-crc64-high");
  const auto low = counter("checkpoint.local-crc64-low");
  HUNDUN_CHECK(high != failed_counters.record.counters.end());
  HUNDUN_CHECK(low != failed_counters.record.counters.end());
  const auto digest = independent_rank_crc_digest(mpi.size());
  HUNDUN_CHECK(high->value == (digest >> 32U));
  HUNDUN_CHECK(low->value == static_cast<std::uint32_t>(digest));

  auto unsupported = request(
      mpi.rank(), hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
      hundun::diagnostics::DiagnosticScope::local);
  unsupported.sample_budget = 1U;
  require_request_error(unsupported,
                        hundun::diagnostics::DiagnosticFailureClass::capability,
                        "checkpoint-v2.diagnostics.capability");
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  const std::string mode(argv[1]);
  if (mode != "fast" && mode != "acceptance")
    return 2;
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    run(mpi, mode == "acceptance");
  });
}
