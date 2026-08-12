// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_performance_artifact.hpp"
#include "hundun/diag_stage3_performance.hpp"
#include "tests/support/test_main.hpp"
#include "yyjson.h"

#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using hundun::diagnostics::Artifact;
using hundun::diagnostics::ArtifactMetadata;
using hundun::diagnostics::CellExtents;
using hundun::diagnostics::ComparisonMode;
using hundun::diagnostics::ComparisonResult;
using hundun::diagnostics::ComparisonStatus;
using hundun::diagnostics::CompatibilityMetadata;
using hundun::diagnostics::ProcessGrid;
using hundun::diagnostics::RawSample;

void check_finite_near(double actual, double expected, double tolerance,
                       const char *file, int line) {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      !std::isfinite(tolerance) || tolerance < 0.0 ||
      std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             " finite numerical check failed");
  }
}

#define HUNDUN_CHECK_FINITE_NEAR(actual, expected, tolerance)                  \
  check_finite_near((actual), (expected), (tolerance), __FILE__, __LINE__)

template <class Function> void expect_invalid(Function &&function) {
  bool rejected = false;
  try {
    std::forward<Function>(function)();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

CompatibilityMetadata compatible_metadata() {
  CompatibilityMetadata metadata;
  metadata.hardware_identity = "cpu-family-model-stepping";
  metadata.node_identity = "node-a";
  metadata.mpi_identity = "Open MPI 4.1.6";
  metadata.compiler_identity = "Clang";
  metadata.compiler_version = "15.0.7";
  metadata.compiler_flags = "-O3 -DNDEBUG";
  metadata.link_flags = "-fuse-ld=lld";
  metadata.build_type = "Release";
  metadata.cpu_affinity = "rank:core";
  metadata.rank_placement = "ppr:1:core";
  metadata.problem_fingerprint = "sha256:case";
  metadata.numerical_tolerance_contract = "sha256:tolerances";
  metadata.measurement_method = "MPI_Wtime:measured-phase";
  metadata.warmup_steps = 5;
  metadata.measured_steps = 20;
  metadata.repetitions = 2;
  metadata.execution_backend = "cpu_reference";
  metadata.ranks = 2;
  metadata.threads = 1;
  metadata.process_grid = ProcessGrid{2, 1, 1};
  metadata.global_owned_cell_extents = CellExtents{64, 64, 64};
  metadata.per_rank_owned_cell_extents = CellExtents{32, 64, 64};
  return metadata;
}

Artifact make_artifact() {
  Artifact artifact;
  artifact.metadata.commit = "0123456789abcdef";
  artifact.metadata.clean = true;
  artifact.metadata.compatibility = compatible_metadata();
  artifact.aggregation = hundun::diagnostics::aggregate_samples(
      {{1, 1, 80.0, 20}, {0, 0, 20.0, 20}, {1, 0, 60.0, 20}, {0, 1, 40.0, 20}},
      2, 2);
  artifact.correctness.passed = true;
  artifact.correctness.summary = "all numerical contracts passed";
  artifact.comparison = ComparisonResult{
      ComparisonMode::identical, ComparisonStatus::comparable, {}};
  artifact.counters.allocated_bytes["field_storage"] = 4096;
  artifact.counters.halo_payload_bytes["pressure"] = 512;
  artifact.counters.halo_messages["pressure"] = 12;
  artifact.counters.collectives["allreduce"] = 7;
  artifact.counters.collective_logical_payload_bytes["allreduce"] = 112;
  artifact.counters.matvec["pressure"] = 18;
  artifact.counters.preconditioner_applications["pressure"] = 17;
  artifact.counters.logical_io_bytes["checkpoint"] = 8192;
  return artifact;
}

Artifact make_stage3_artifact() {
  Artifact artifact = make_artifact();
  artifact.schema_version = 2;
  artifact.metadata.tree_fingerprint = "sha256:tree";
  artifact.metadata.binary_fingerprint = "sha256:hundun";
  artifact.metadata.profile = "ideal-gas-static-ibm-wale";
  artifact.metadata.geometry_fingerprint = "sha256:geometry";
  artifact.metadata.cpuset = "0-3";
  artifact.metadata.thread_budget = 4;
  hundun::diagnostics::Stage3PerformanceCounters counters;
  counters.init_surface_triangles = 12U;
  counters.step_wale_evaluations = 3U;
  artifact.counters.algorithmic_work =
      hundun::diagnostics::stage3_algorithmic_work_map(counters);
  return artifact;
}

yyjson_val *required(yyjson_val *object, const char *key) {
  yyjson_val *value = yyjson_obj_get(object, key);
  HUNDUN_CHECK(value != nullptr);
  return value;
}

std::string byte_string(std::initializer_list<unsigned int> bytes) {
  std::string value;
  value.reserve(bytes.size());
  for (const unsigned int byte : bytes) {
    HUNDUN_CHECK(byte <= 0xffU);
    value.push_back(static_cast<char>(byte));
  }
  return value;
}

void test_finite_aware_comparison() {
  bool rejected = false;
  try {
    HUNDUN_CHECK_FINITE_NEAR(std::numeric_limits<double>::quiet_NaN(), 0.0,
                             0.0);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_elapsed_phase_measurement() {
  std::vector<double> times{10.0, 10.25};
  std::size_t next = 0;
  bool phase_called = false;
  const double elapsed = hundun::diagnostics::measure_elapsed_phase(
      [&] { return times.at(next++); }, [&] { phase_called = true; });
  HUNDUN_CHECK(phase_called);
  HUNDUN_CHECK_FINITE_NEAR(elapsed, 0.25, 0.0);

  expect_invalid([] {
    std::vector<double> backwards{2.0, 1.0};
    std::size_t index = 0;
    static_cast<void>(hundun::diagnostics::measure_elapsed_phase(
        [&] { return backwards.at(index++); }, [] {}));
  });
  expect_invalid([] {
    static_cast<void>(hundun::diagnostics::measure_elapsed_phase(
        [] { return std::numeric_limits<double>::infinity(); }, [] {}));
  });
  expect_invalid([] {
    std::vector<double> times{1.0, std::numeric_limits<double>::quiet_NaN()};
    std::size_t index = 0;
    static_cast<void>(hundun::diagnostics::measure_elapsed_phase(
        [&] { return times.at(index++); }, [] {}));
  });
  expect_invalid([] {
    std::vector<double> times{-std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::max()};
    std::size_t index = 0;
    static_cast<void>(hundun::diagnostics::measure_elapsed_phase(
        [&] { return times.at(index++); }, [] {}));
  });
}

void test_aggregation_and_medians() {
  const auto changing_slowest =
      hundun::diagnostics::aggregate_samples({{1, 2, 10.0, 2},
                                              {0, 2, 5.0, 2},
                                              {1, 0, 8.0, 2},
                                              {0, 0, 4.0, 2},
                                              {1, 1, 5.0, 2},
                                              {0, 1, 6.0, 2}},
                                             2, 3);
  HUNDUN_CHECK(changing_slowest.raw_samples.size() == 6U);
  for (std::size_t index = 0; index < changing_slowest.raw_samples.size();
       ++index) {
    HUNDUN_CHECK(changing_slowest.raw_samples[index].repetition ==
                 static_cast<int>(index / 3U));
    HUNDUN_CHECK(changing_slowest.raw_samples[index].relative_rank ==
                 static_cast<int>(index % 3U));
  }
  HUNDUN_CHECK(changing_slowest.repetition_maxima.size() == 2U);
  HUNDUN_CHECK(changing_slowest.repetition_maxima[0].slowest_relative_rank ==
               1);
  HUNDUN_CHECK_FINITE_NEAR(changing_slowest.repetition_maxima[0].step_seconds,
                           3.0, 0.0);
  HUNDUN_CHECK(changing_slowest.repetition_maxima[1].slowest_relative_rank ==
               2);
  HUNDUN_CHECK_FINITE_NEAR(changing_slowest.repetition_maxima[1].step_seconds,
                           5.0, 0.0);
  HUNDUN_CHECK_FINITE_NEAR(changing_slowest.median_step_seconds, 4.0, 0.0);

  const auto odd = hundun::diagnostics::aggregate_samples(
      {{0, 0, 9.0, 1}, {1, 0, 1.0, 1}, {2, 0, 5.0, 1}}, 3, 1);
  HUNDUN_CHECK_FINITE_NEAR(odd.median_step_seconds, 5.0, 0.0);
  const auto even = hundun::diagnostics::aggregate_samples(
      {{0, 0, 1.0, 1}, {1, 0, 9.0, 1}, {2, 0, 3.0, 1}, {3, 0, 7.0, 1}}, 4, 1);
  HUNDUN_CHECK_FINITE_NEAR(even.median_step_seconds, 5.0, 0.0);
}

void test_invalid_aggregation_inputs() {
  expect_invalid([] {
    static_cast<void>(hundun::diagnostics::aggregate_samples({}, 0, 1));
  });
  expect_invalid([] {
    static_cast<void>(hundun::diagnostics::aggregate_samples({}, 1, 0));
  });
  expect_invalid([] {
    static_cast<void>(hundun::diagnostics::aggregate_samples(
        {{0, 0, 1.0, 1}, {0, 0, 2.0, 1}}, 1, 2));
  });
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::aggregate_samples({{0, 0, 1.0, 1}}, 1, 2));
  });
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::aggregate_samples({{-1, 0, 1.0, 1}}, 1, 1));
  });
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::aggregate_samples({{1, 0, 1.0, 1}}, 1, 1));
  });
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::aggregate_samples({{0, -1, 1.0, 1}}, 1, 1));
  });
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::aggregate_samples({{0, 1, 1.0, 1}}, 1, 1));
  });
  expect_invalid([] {
    static_cast<void>(hundun::diagnostics::aggregate_samples(
        {{0, 0, 1.0, 2}, {0, 1, 1.0, 3}}, 1, 2));
  });
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::aggregate_samples({{0, 0, 1.0, 0}}, 1, 1));
  });
  for (const double elapsed :
       {0.0, -1.0, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN()}) {
    expect_invalid([elapsed] {
      static_cast<void>(
          hundun::diagnostics::aggregate_samples({{0, 0, elapsed, 1}}, 1, 1));
    });
  }
}

void test_compatibility_keys_and_scaling_rules() {
  const CompatibilityMetadata baseline = compatible_metadata();
  const auto equal = hundun::diagnostics::compare_compatibility(
      baseline, baseline, ComparisonMode::identical);
  HUNDUN_CHECK(equal.status == ComparisonStatus::comparable);
  HUNDUN_CHECK(equal.reasons.empty());

  CompatibilityMetadata changed = baseline;
  changed.hardware_identity = "different";
  changed.node_identity = "different";
  changed.mpi_identity = "different";
  changed.compiler_identity = "different";
  changed.compiler_version = "different";
  changed.compiler_flags = "different";
  changed.link_flags = "different";
  changed.build_type = "different";
  changed.cpu_affinity = "different";
  changed.rank_placement = "different";
  changed.problem_fingerprint = "different";
  changed.numerical_tolerance_contract = "different";
  changed.measurement_method = "different";
  changed.warmup_steps = 6;
  changed.measured_steps = 21;
  changed.repetitions = 3;
  changed.execution_backend = "different";
  changed.ranks = 4;
  changed.threads = 2;
  changed.process_grid = ProcessGrid{4, 1, 1};
  changed.global_owned_cell_extents = CellExtents{128, 64, 64};
  changed.per_rank_owned_cell_extents = CellExtents{32, 32, 64};
  const auto mismatch = hundun::diagnostics::compare_compatibility(
      baseline, changed, ComparisonMode::identical);
  HUNDUN_CHECK(mismatch.status == ComparisonStatus::incomparable);
  const std::vector<std::string> expected{"build_type",
                                          "compiler_flags",
                                          "compiler_identity",
                                          "compiler_version",
                                          "cpu_affinity",
                                          "execution_backend",
                                          "global_owned_cell_extents",
                                          "hardware_identity",
                                          "link_flags",
                                          "measured_steps",
                                          "measurement_method",
                                          "mpi_identity",
                                          "node_identity",
                                          "numerical_tolerance_contract",
                                          "per_rank_owned_cell_extents",
                                          "problem_fingerprint",
                                          "process_grid",
                                          "rank_placement",
                                          "ranks",
                                          "repetitions",
                                          "threads",
                                          "warmup_steps"};
  HUNDUN_CHECK(mismatch.reasons == expected);

  CompatibilityMetadata scaled = baseline;
  scaled.ranks = 4;
  scaled.process_grid = ProcessGrid{4, 1, 1};
  scaled.per_rank_owned_cell_extents = CellExtents{16, 64, 64};
  auto result = hundun::diagnostics::compare_compatibility(
      baseline, scaled, ComparisonMode::strong_scaling);
  HUNDUN_CHECK(result.status == ComparisonStatus::comparable);
  HUNDUN_CHECK(result.reasons.empty());
  scaled.global_owned_cell_extents = CellExtents{128, 64, 64};
  result = hundun::diagnostics::compare_compatibility(
      baseline, scaled, ComparisonMode::strong_scaling);
  HUNDUN_CHECK(result.reasons ==
               std::vector<std::string>{"global_owned_cell_extents"});

  scaled = baseline;
  scaled.ranks = 4;
  scaled.process_grid = ProcessGrid{4, 1, 1};
  scaled.global_owned_cell_extents = CellExtents{128, 64, 64};
  result = hundun::diagnostics::compare_compatibility(
      baseline, scaled, ComparisonMode::weak_scaling);
  HUNDUN_CHECK(result.status == ComparisonStatus::comparable);
  scaled.per_rank_owned_cell_extents = CellExtents{64, 64, 64};
  result = hundun::diagnostics::compare_compatibility(
      baseline, scaled, ComparisonMode::weak_scaling);
  HUNDUN_CHECK(result.reasons ==
               std::vector<std::string>{"per_rank_owned_cell_extents"});
}

void test_scaling_formulas() {
  HUNDUN_CHECK_FINITE_NEAR(
      hundun::diagnostics::strong_scaling_speedup(8.0, 2.0), 4.0, 0.0);
  HUNDUN_CHECK_FINITE_NEAR(
      hundun::diagnostics::strong_scaling_efficiency(8.0, 2.0, 4), 1.0, 0.0);
  HUNDUN_CHECK_FINITE_NEAR(
      hundun::diagnostics::weak_scaling_efficiency(8.0, 10.0), 0.8, 0.0);
  for (const double invalid :
       {0.0, -1.0, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN()}) {
    expect_invalid([invalid] {
      static_cast<void>(
          hundun::diagnostics::strong_scaling_speedup(invalid, 1.0));
    });
    expect_invalid([invalid] {
      static_cast<void>(
          hundun::diagnostics::strong_scaling_speedup(1.0, invalid));
    });
    expect_invalid([invalid] {
      static_cast<void>(
          hundun::diagnostics::weak_scaling_efficiency(invalid, 1.0));
    });
    expect_invalid([invalid] {
      static_cast<void>(
          hundun::diagnostics::weak_scaling_efficiency(1.0, invalid));
    });
  }
  expect_invalid([] {
    static_cast<void>(
        hundun::diagnostics::strong_scaling_efficiency(1.0, 1.0, 0));
  });
}

void test_json_utf8_contract() {
  const std::string valid_multibyte = u8"节点🔥";
  Artifact artifact = make_artifact();
  artifact.metadata.compatibility.node_identity = valid_multibyte;
  artifact.counters.allocated_bytes[valid_multibyte] = 7U;
  const std::string encoded = hundun::diagnostics::to_json(artifact);
  yyjson_doc *document = yyjson_read(encoded.data(), encoded.size(), 0);
  HUNDUN_CHECK(document != nullptr);
  yyjson_val *root = yyjson_doc_get_root(document);
  HUNDUN_CHECK(std::string(yyjson_get_str(required(root, "node_identity"))) ==
               valid_multibyte);
  yyjson_val *allocated =
      required(required(root, "exact_counters"), "allocated_bytes");
  yyjson_val *multibyte_counter =
      yyjson_obj_get(allocated, valid_multibyte.c_str());
  HUNDUN_CHECK(multibyte_counter != nullptr);
  HUNDUN_CHECK(yyjson_get_uint(multibyte_counter) == 7U);
  yyjson_doc_free(document);

  const std::vector<std::string> malformed{
      byte_string({0x80U}),                      // Isolated continuation.
      byte_string({0xe2U, 0x82U}),               // Truncated sequence.
      byte_string({0xc0U, 0xafU}),               // Overlong encoding.
      byte_string({0xedU, 0xa0U, 0x80U}),        // UTF-16 surrogate.
      byte_string({0xf4U, 0x90U, 0x80U, 0x80U}), // Above U+10FFFF.
  };
  for (const std::string &invalid_utf8 : malformed) {
    artifact = make_artifact();
    artifact.metadata.compatibility.node_identity = invalid_utf8;
    expect_invalid(
        [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });

    artifact = make_artifact();
    artifact.counters.allocated_bytes[invalid_utf8] = 1U;
    expect_invalid(
        [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  }
}

void test_deterministic_json_contract() {
  Artifact artifact = make_artifact();
  artifact.metadata.commit = "commit\"with\ncharacters";
  const std::string first = hundun::diagnostics::to_json(artifact);
  const std::string second = hundun::diagnostics::to_json(artifact);
  HUNDUN_CHECK(first == second);

  yyjson_doc *document = yyjson_read(first.data(), first.size(), 0);
  HUNDUN_CHECK(document != nullptr);
  yyjson_val *root = yyjson_doc_get_root(document);
  HUNDUN_CHECK(yyjson_is_obj(root));
  HUNDUN_CHECK(yyjson_is_int(required(root, "schema_version")));
  HUNDUN_CHECK(yyjson_get_sint(required(root, "schema_version")) == 1);
  HUNDUN_CHECK(yyjson_is_str(required(root, "commit")));
  HUNDUN_CHECK(std::string(yyjson_get_str(required(root, "commit"))) ==
               artifact.metadata.commit);
  HUNDUN_CHECK(yyjson_is_obj(required(root, "working_tree")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "case_fingerprint")));
  HUNDUN_CHECK(yyjson_is_obj(required(root, "compiler")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "link_flags")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "build_type")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "mpi_implementation")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "node_identity")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "hardware_identity")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "cpu_affinity")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "rank_placement")));
  HUNDUN_CHECK(yyjson_is_int(required(root, "ranks")));
  HUNDUN_CHECK(yyjson_is_int(required(root, "threads")));
  HUNDUN_CHECK(yyjson_is_arr(required(root, "process_grid")));
  HUNDUN_CHECK(yyjson_is_arr(required(root, "global_owned_cell_extents")));
  HUNDUN_CHECK(yyjson_is_arr(required(root, "per_rank_owned_cell_extents")));
  HUNDUN_CHECK(yyjson_is_str(required(root, "execution_backend")));
  HUNDUN_CHECK(
      yyjson_is_str(required(root, "numerical_tolerance_fingerprint")));
  HUNDUN_CHECK(yyjson_is_obj(required(root, "measurement")));
  HUNDUN_CHECK(yyjson_is_obj(required(root, "correctness")));

  yyjson_val *raw = required(root, "raw_samples");
  HUNDUN_CHECK(yyjson_is_arr(raw));
  HUNDUN_CHECK(yyjson_arr_size(raw) == 4U);
  yyjson_val *first_sample = yyjson_arr_get(raw, 0);
  HUNDUN_CHECK(yyjson_is_int(required(first_sample, "repetition")));
  HUNDUN_CHECK(yyjson_is_int(required(first_sample, "relative_rank")));
  HUNDUN_CHECK(yyjson_is_num(required(first_sample, "elapsed_seconds")));
  HUNDUN_CHECK(yyjson_is_int(required(first_sample, "measured_steps")));

  yyjson_val *maxima = required(root, "repetition_maxima");
  HUNDUN_CHECK(yyjson_is_arr(maxima));
  HUNDUN_CHECK(yyjson_arr_size(maxima) == 2U);
  HUNDUN_CHECK_FINITE_NEAR(
      yyjson_get_real(required(yyjson_arr_get(maxima, 0), "step_seconds")), 2.0,
      0.0);
  HUNDUN_CHECK_FINITE_NEAR(
      yyjson_get_real(required(yyjson_arr_get(maxima, 1), "step_seconds")), 4.0,
      0.0);
  HUNDUN_CHECK_FINITE_NEAR(
      yyjson_get_real(required(root, "median_step_seconds")), 3.0, 0.0);
  HUNDUN_CHECK(yyjson_is_obj(required(root, "comparison")));
  yyjson_val *counters = required(root, "exact_counters");
  HUNDUN_CHECK(yyjson_is_obj(counters));
  HUNDUN_CHECK(yyjson_is_obj(required(counters, "allocated_bytes")));
  HUNDUN_CHECK(yyjson_is_obj(required(counters, "halo_payload_bytes")));
  HUNDUN_CHECK(yyjson_is_obj(required(counters, "halo_messages")));
  HUNDUN_CHECK(yyjson_is_obj(required(counters, "collectives")));
  HUNDUN_CHECK(
      yyjson_is_obj(required(counters, "collective_logical_payload_bytes")));
  HUNDUN_CHECK(yyjson_is_obj(required(counters, "matvec")));
  HUNDUN_CHECK(
      yyjson_is_obj(required(counters, "preconditioner_applications")));
  HUNDUN_CHECK(yyjson_is_obj(required(counters, "logical_io_bytes")));
  yyjson_doc_free(document);
}

void test_artifact_validation() {
  Artifact artifact = make_artifact();
  artifact.metadata.clean = false;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  artifact.metadata.dirty_summary = "diff-sha256:abc";
  HUNDUN_CHECK(!hundun::diagnostics::to_json(artifact).empty());
  artifact.metadata.clean = true;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });

  artifact = make_artifact();
  artifact.aggregation.raw_samples[0].elapsed_seconds =
      std::numeric_limits<double>::infinity();
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });

  artifact = make_artifact();
  artifact.metadata.compatibility.ranks = 4;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  artifact = make_artifact();
  artifact.aggregation.median_step_seconds = 99.0;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  artifact = make_artifact();
  artifact.metadata.compatibility.process_grid = ProcessGrid{1, 1, 1};
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  artifact = make_artifact();
  artifact.comparison.status = ComparisonStatus::incomparable;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
}

void test_stage3_schema_v2() {
  Artifact artifact = make_stage3_artifact();
  const auto encoded = hundun::diagnostics::to_json(artifact);
  yyjson_doc *document = yyjson_read(encoded.data(), encoded.size(), 0);
  HUNDUN_CHECK(document != nullptr);
  yyjson_val *root = yyjson_doc_get_root(document);
  HUNDUN_CHECK(yyjson_get_sint(required(root, "schema_version")) == 2);
  yyjson_val *identity = required(root, "stage3_identity");
  HUNDUN_CHECK(std::string(yyjson_get_str(
                   required(identity, "binary_fingerprint"))) ==
               artifact.metadata.binary_fingerprint);
  HUNDUN_CHECK(yyjson_get_sint(required(identity, "thread_budget")) == 4);
  yyjson_val *algorithmic =
      required(required(root, "exact_counters"), "algorithmic_work");
  HUNDUN_CHECK(yyjson_obj_size(algorithmic) == 17U);
  HUNDUN_CHECK(yyjson_get_uint(
                   required(algorithmic, "init.surface.triangles")) == 12U);
  yyjson_doc_free(document);

  artifact.counters.algorithmic_work.erase("step.wale.evaluations");
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  artifact = make_stage3_artifact();
  artifact.counters.algorithmic_work["unknown"] = 1U;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
  artifact = make_stage3_artifact();
  artifact.metadata.binary_fingerprint.clear();
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });

  artifact = make_stage3_artifact();
  Artifact changed = artifact;
  changed.metadata.tree_fingerprint = "sha256:different-tree";
  auto comparison = hundun::diagnostics::compare_artifact_metadata(
      artifact.metadata, changed.metadata, ComparisonMode::identical);
  HUNDUN_CHECK(comparison.status == ComparisonStatus::incomparable);
  HUNDUN_CHECK(comparison.reasons ==
               std::vector<std::string>{
                   "source.tree-fingerprint.mismatch"});
  changed = artifact;
  ++changed.metadata.thread_budget;
  comparison = hundun::diagnostics::compare_artifact_metadata(
      artifact.metadata, changed.metadata, ComparisonMode::identical);
  HUNDUN_CHECK(comparison.reasons ==
               std::vector<std::string>{
                   "platform.thread-budget.mismatch"});

  artifact = make_artifact();
  const auto schema_v1 = hundun::diagnostics::to_json(artifact);
  HUNDUN_CHECK(schema_v1.find("stage3_identity") == std::string::npos);
  HUNDUN_CHECK(schema_v1.find("algorithmic_work") == std::string::npos);
  artifact.counters.algorithmic_work["init.surface.triangles"] = 1U;
  expect_invalid(
      [&] { static_cast<void>(hundun::diagnostics::to_json(artifact)); });
}
} // namespace

int main() {
  return hundun::test::run([] {
    test_finite_aware_comparison();
    test_elapsed_phase_measurement();
    test_aggregation_and_medians();
    test_invalid_aggregation_inputs();
    test_compatibility_keys_and_scaling_rules();
    test_scaling_formulas();
    test_json_utf8_contract();
    test_deterministic_json_contract();
    test_artifact_validation();
    test_stage3_schema_v2();
  });
}
