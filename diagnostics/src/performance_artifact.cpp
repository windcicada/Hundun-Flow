// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/performance_artifact.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>

namespace hundun::diagnostics {
namespace {

bool equal(ProcessGrid left, ProcessGrid right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool equal(CellExtents left, CellExtents right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool equal(const RawSample &left, const RawSample &right) {
  return left.repetition == right.repetition &&
         left.relative_rank == right.relative_rank &&
         left.elapsed_seconds == right.elapsed_seconds &&
         left.measured_steps == right.measured_steps;
}

bool equal(const RepetitionMaximum &left, const RepetitionMaximum &right) {
  return left.repetition == right.repetition &&
         left.slowest_relative_rank == right.slowest_relative_rank &&
         left.step_seconds == right.step_seconds;
}

void require_nonempty(const std::string &value, std::string_view key) {
  if (value.empty()) {
    throw std::invalid_argument(std::string(key) + " must not be empty");
  }
}

std::uint64_t checked_grid_size(ProcessGrid grid) {
  if (grid.x <= 0 || grid.y <= 0 || grid.z <= 0) {
    throw std::invalid_argument("process_grid entries must be positive");
  }
  const auto x = static_cast<std::uint64_t>(grid.x);
  const auto y = static_cast<std::uint64_t>(grid.y);
  const auto z = static_cast<std::uint64_t>(grid.z);
  if (x > std::numeric_limits<std::uint64_t>::max() / y ||
      x * y > std::numeric_limits<std::uint64_t>::max() / z) {
    throw std::invalid_argument("process_grid product overflows");
  }
  return x * y * z;
}

void validate_extents(CellExtents extents, std::string_view key) {
  if (extents.x == 0U || extents.y == 0U || extents.z == 0U) {
    throw std::invalid_argument(std::string(key) + " entries must be positive");
  }
}

void validate_compatibility(const CompatibilityMetadata &metadata) {
  require_nonempty(metadata.hardware_identity, "hardware_identity");
  require_nonempty(metadata.node_identity, "node_identity");
  require_nonempty(metadata.mpi_identity, "mpi_identity");
  require_nonempty(metadata.compiler_identity, "compiler_identity");
  require_nonempty(metadata.compiler_version, "compiler_version");
  require_nonempty(metadata.build_type, "build_type");
  require_nonempty(metadata.cpu_affinity, "cpu_affinity");
  require_nonempty(metadata.rank_placement, "rank_placement");
  require_nonempty(metadata.problem_fingerprint, "problem_fingerprint");
  require_nonempty(metadata.numerical_tolerance_contract,
                   "numerical_tolerance_contract");
  require_nonempty(metadata.measurement_method, "measurement_method");
  require_nonempty(metadata.execution_backend, "execution_backend");
  if (metadata.measured_steps == 0U || metadata.repetitions <= 0 ||
      metadata.ranks <= 0 || metadata.threads <= 0) {
    throw std::invalid_argument(
        "measured_steps, repetitions, ranks, and threads must be positive");
  }
  if (checked_grid_size(metadata.process_grid) !=
      static_cast<std::uint64_t>(metadata.ranks)) {
    throw std::invalid_argument("process_grid product must equal ranks");
  }
  validate_extents(metadata.global_owned_cell_extents,
                   "global_owned_cell_extents");
  validate_extents(metadata.per_rank_owned_cell_extents,
                   "per_rank_owned_cell_extents");
}

void add_mismatch(std::vector<std::string> &reasons, bool mismatch,
                  const char *key) {
  if (mismatch) {
    reasons.emplace_back(key);
  }
}

void validate_positive_time(double value, std::string_view key) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw std::invalid_argument(std::string(key) +
                                " must be positive and finite");
  }
}

const char *mode_name(ComparisonMode mode) {
  switch (mode) {
  case ComparisonMode::identical:
    return "identical";
  case ComparisonMode::strong_scaling:
    return "strong_scaling";
  case ComparisonMode::weak_scaling:
    return "weak_scaling";
  }
  throw std::invalid_argument("invalid comparison mode");
}

const char *status_name(ComparisonStatus status) {
  switch (status) {
  case ComparisonStatus::comparable:
    return "comparable";
  case ComparisonStatus::incomparable:
    return "incomparable";
  }
  throw std::invalid_argument("invalid comparison status");
}

void append_json_string(std::string &output, std::string_view value) {
  constexpr char hexadecimal[] = "0123456789abcdef";
  output.push_back('"');
  for (const char source_character : value) {
    const auto character = static_cast<unsigned char>(source_character);
    switch (character) {
    case '"':
      output += "\\\"";
      break;
    case '\\':
      output += "\\\\";
      break;
    case '\b':
      output += "\\b";
      break;
    case '\f':
      output += "\\f";
      break;
    case '\n':
      output += "\\n";
      break;
    case '\r':
      output += "\\r";
      break;
    case '\t':
      output += "\\t";
      break;
    default:
      if (character < 0x20U) {
        output += "\\u00";
        output.push_back(hexadecimal[character >> 4U]);
        output.push_back(hexadecimal[character & 0x0fU]);
      } else {
        output.push_back(static_cast<char>(character));
      }
    }
  }
  output.push_back('"');
}

void append_double(std::string &output, double value) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument("JSON numeric values must be finite");
  }
  std::ostringstream stream;
  stream.imbue(std::locale::classic());
  stream << std::setprecision(std::numeric_limits<double>::max_digits10)
         << value;
  std::string encoded = stream.str();
  if (encoded.find_first_of(".eE") == std::string::npos) {
    encoded += ".0";
  }
  output += encoded;
}

void append_extents(std::string &output, CellExtents extents) {
  output += '[' + std::to_string(extents.x) + ',' + std::to_string(extents.y) +
            ',' + std::to_string(extents.z) + ']';
}

void append_grid(std::string &output, ProcessGrid grid) {
  output += '[' + std::to_string(grid.x) + ',' + std::to_string(grid.y) + ',' +
            std::to_string(grid.z) + ']';
}

void validate_counter_map(const CounterMap &counters) {
  for (const auto &[key, value] : counters) {
    static_cast<void>(value);
    require_nonempty(key, "counter key");
  }
}

void append_counter_map(std::string &output, const CounterMap &counters) {
  output.push_back('{');
  bool first = true;
  for (const auto &[key, value] : counters) {
    if (!first) {
      output.push_back(',');
    }
    first = false;
    append_json_string(output, key);
    output.push_back(':');
    output += std::to_string(value);
  }
  output.push_back('}');
}

void require_same_aggregation(const SampleAggregation &recorded,
                              const SampleAggregation &recomputed) {
  if (recorded.repetitions != recomputed.repetitions ||
      recorded.ranks != recomputed.ranks ||
      recorded.measured_steps != recomputed.measured_steps ||
      recorded.median_step_seconds != recomputed.median_step_seconds ||
      recorded.raw_samples.size() != recomputed.raw_samples.size() ||
      recorded.repetition_maxima.size() !=
          recomputed.repetition_maxima.size()) {
    throw std::invalid_argument("artifact aggregation is inconsistent");
  }
  for (std::size_t index = 0; index < recorded.raw_samples.size(); ++index) {
    if (!equal(recorded.raw_samples[index], recomputed.raw_samples[index])) {
      throw std::invalid_argument("artifact raw samples are inconsistent");
    }
  }
  for (std::size_t index = 0; index < recorded.repetition_maxima.size();
       ++index) {
    if (!equal(recorded.repetition_maxima[index],
               recomputed.repetition_maxima[index])) {
      throw std::invalid_argument(
          "artifact repetition maxima are inconsistent");
    }
  }
}

void validate_comparison(const ComparisonResult &comparison) {
  static_cast<void>(mode_name(comparison.mode));
  static_cast<void>(status_name(comparison.status));
  const bool sorted =
      std::is_sorted(comparison.reasons.begin(), comparison.reasons.end());
  const bool unique =
      std::adjacent_find(comparison.reasons.begin(),
                         comparison.reasons.end()) == comparison.reasons.end();
  for (const std::string &reason : comparison.reasons) {
    require_nonempty(reason, "comparison reason");
  }
  if (!sorted || !unique ||
      (comparison.status == ComparisonStatus::comparable &&
       !comparison.reasons.empty()) ||
      (comparison.status == ComparisonStatus::incomparable &&
       comparison.reasons.empty())) {
    throw std::invalid_argument(
        "comparison status and reasons are inconsistent");
  }
}

} // namespace

SampleAggregation aggregate_samples(std::vector<RawSample> samples,
                                    int repetitions, int ranks) {
  if (repetitions <= 0 || ranks <= 0) {
    throw std::invalid_argument("repetition and rank domains must be positive");
  }
  const auto repetition_count = static_cast<std::uint64_t>(repetitions);
  const auto rank_count = static_cast<std::uint64_t>(ranks);
  if (repetition_count >
      std::numeric_limits<std::uint64_t>::max() / rank_count) {
    throw std::invalid_argument("sample domain size overflows");
  }
  const std::uint64_t expected = repetition_count * rank_count;
  if (expected != static_cast<std::uint64_t>(samples.size())) {
    throw std::invalid_argument(
        "exactly one sample is required for each repetition/rank pair");
  }

  std::sort(samples.begin(), samples.end(),
            [](const RawSample &left, const RawSample &right) {
              return std::tie(left.repetition, left.relative_rank) <
                     std::tie(right.repetition, right.relative_rank);
            });
  std::uint64_t common_steps = 0;
  for (std::size_t index = 0; index < samples.size(); ++index) {
    const RawSample &sample = samples[index];
    if (sample.repetition < 0 || sample.repetition >= repetitions ||
        sample.relative_rank < 0 || sample.relative_rank >= ranks) {
      throw std::invalid_argument("sample lies outside its declared domain");
    }
    validate_positive_time(sample.elapsed_seconds, "elapsed_seconds");
    if (sample.measured_steps == 0U) {
      throw std::invalid_argument("measured_steps must be positive");
    }
    if (common_steps == 0U) {
      common_steps = sample.measured_steps;
    } else if (sample.measured_steps != common_steps) {
      throw std::invalid_argument(
          "measured_steps must be common to all samples");
    }
    const auto expected_repetition = static_cast<int>(index / rank_count);
    const auto expected_rank = static_cast<int>(index % rank_count);
    if (sample.repetition != expected_repetition ||
        sample.relative_rank != expected_rank) {
      throw std::invalid_argument(
          "duplicate or missing repetition/rank sample");
    }
  }

  SampleAggregation result;
  result.repetitions = repetitions;
  result.ranks = ranks;
  result.measured_steps = common_steps;
  result.raw_samples = std::move(samples);
  result.repetition_maxima.reserve(static_cast<std::size_t>(repetitions));
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    RepetitionMaximum maximum{repetition, 0, -1.0};
    for (int rank = 0; rank < ranks; ++rank) {
      const auto index = static_cast<std::size_t>(repetition) *
                             static_cast<std::size_t>(ranks) +
                         static_cast<std::size_t>(rank);
      const double step_seconds = result.raw_samples[index].elapsed_seconds /
                                  static_cast<double>(common_steps);
      validate_positive_time(step_seconds, "rank step time");
      if (step_seconds > maximum.step_seconds) {
        maximum.slowest_relative_rank = rank;
        maximum.step_seconds = step_seconds;
      }
    }
    result.repetition_maxima.push_back(maximum);
  }

  std::vector<double> ordered_maxima;
  ordered_maxima.reserve(result.repetition_maxima.size());
  for (const RepetitionMaximum &maximum : result.repetition_maxima) {
    ordered_maxima.push_back(maximum.step_seconds);
  }
  std::sort(ordered_maxima.begin(), ordered_maxima.end());
  const std::size_t middle = ordered_maxima.size() / 2U;
  if (ordered_maxima.size() % 2U == 1U) {
    result.median_step_seconds = ordered_maxima[middle];
  } else {
    result.median_step_seconds =
        ordered_maxima[middle - 1U] +
        (ordered_maxima[middle] - ordered_maxima[middle - 1U]) / 2.0;
  }
  validate_positive_time(result.median_step_seconds, "median step time");
  return result;
}

double measure_elapsed_phase(const ClockFunction &clock,
                             const PhaseFunction &phase) {
  if (!clock || !phase) {
    throw std::invalid_argument("clock and phase callbacks are required");
  }
  const double start = clock();
  if (!std::isfinite(start)) {
    throw std::invalid_argument("phase start time must be finite");
  }
  phase();
  const double finish = clock();
  if (!std::isfinite(finish) || finish < start) {
    throw std::invalid_argument(
        "phase finish time must be finite and not precede start");
  }
  const double elapsed = finish - start;
  if (!std::isfinite(elapsed)) {
    throw std::invalid_argument("phase elapsed interval must be finite");
  }
  return elapsed;
}

ComparisonResult compare_compatibility(const CompatibilityMetadata &baseline,
                                       const CompatibilityMetadata &candidate,
                                       ComparisonMode mode) {
  validate_compatibility(baseline);
  validate_compatibility(candidate);
  static_cast<void>(mode_name(mode));
  std::vector<std::string> reasons;
  add_mismatch(reasons, baseline.build_type != candidate.build_type,
               "build_type");
  add_mismatch(reasons, baseline.compiler_flags != candidate.compiler_flags,
               "compiler_flags");
  add_mismatch(reasons,
               baseline.compiler_identity != candidate.compiler_identity,
               "compiler_identity");
  add_mismatch(reasons, baseline.compiler_version != candidate.compiler_version,
               "compiler_version");
  add_mismatch(reasons, baseline.cpu_affinity != candidate.cpu_affinity,
               "cpu_affinity");
  add_mismatch(reasons,
               baseline.execution_backend != candidate.execution_backend,
               "execution_backend");
  add_mismatch(reasons,
               baseline.hardware_identity != candidate.hardware_identity,
               "hardware_identity");
  add_mismatch(reasons, baseline.link_flags != candidate.link_flags,
               "link_flags");
  add_mismatch(reasons, baseline.measured_steps != candidate.measured_steps,
               "measured_steps");
  add_mismatch(reasons,
               baseline.measurement_method != candidate.measurement_method,
               "measurement_method");
  add_mismatch(reasons, baseline.mpi_identity != candidate.mpi_identity,
               "mpi_identity");
  add_mismatch(reasons, baseline.node_identity != candidate.node_identity,
               "node_identity");
  add_mismatch(reasons,
               baseline.numerical_tolerance_contract !=
                   candidate.numerical_tolerance_contract,
               "numerical_tolerance_contract");
  add_mismatch(reasons,
               baseline.problem_fingerprint != candidate.problem_fingerprint,
               "problem_fingerprint");
  add_mismatch(reasons, baseline.rank_placement != candidate.rank_placement,
               "rank_placement");
  add_mismatch(reasons, baseline.repetitions != candidate.repetitions,
               "repetitions");
  add_mismatch(reasons, baseline.threads != candidate.threads, "threads");
  add_mismatch(reasons, baseline.warmup_steps != candidate.warmup_steps,
               "warmup_steps");

  if (mode == ComparisonMode::identical) {
    add_mismatch(reasons,
                 !equal(baseline.global_owned_cell_extents,
                        candidate.global_owned_cell_extents),
                 "global_owned_cell_extents");
    add_mismatch(reasons,
                 !equal(baseline.per_rank_owned_cell_extents,
                        candidate.per_rank_owned_cell_extents),
                 "per_rank_owned_cell_extents");
    add_mismatch(reasons, !equal(baseline.process_grid, candidate.process_grid),
                 "process_grid");
    add_mismatch(reasons, baseline.ranks != candidate.ranks, "ranks");
  } else if (mode == ComparisonMode::strong_scaling) {
    add_mismatch(reasons,
                 !equal(baseline.global_owned_cell_extents,
                        candidate.global_owned_cell_extents),
                 "global_owned_cell_extents");
  } else {
    add_mismatch(reasons,
                 !equal(baseline.per_rank_owned_cell_extents,
                        candidate.per_rank_owned_cell_extents),
                 "per_rank_owned_cell_extents");
  }

  std::sort(reasons.begin(), reasons.end());
  return ComparisonResult{mode,
                          reasons.empty() ? ComparisonStatus::comparable
                                          : ComparisonStatus::incomparable,
                          std::move(reasons)};
}

double strong_scaling_speedup(double single_rank_step_seconds,
                              double parallel_step_seconds) {
  validate_positive_time(single_rank_step_seconds, "single-rank step time");
  validate_positive_time(parallel_step_seconds, "parallel step time");
  const double value = single_rank_step_seconds / parallel_step_seconds;
  validate_positive_time(value, "strong-scaling speedup");
  return value;
}

double strong_scaling_efficiency(double single_rank_step_seconds,
                                 double parallel_step_seconds, int ranks) {
  if (ranks <= 0) {
    throw std::invalid_argument("ranks must be positive");
  }
  const double value =
      strong_scaling_speedup(single_rank_step_seconds, parallel_step_seconds) /
      static_cast<double>(ranks);
  validate_positive_time(value, "strong-scaling efficiency");
  return value;
}

double weak_scaling_efficiency(double single_rank_step_seconds,
                               double parallel_step_seconds) {
  return strong_scaling_speedup(single_rank_step_seconds,
                                parallel_step_seconds);
}

std::string to_json(const Artifact &artifact) {
  if (artifact.schema_version != 1) {
    throw std::invalid_argument(
        "performance artifact schema_version must be 1");
  }
  require_nonempty(artifact.metadata.commit, "commit");
  if ((artifact.metadata.clean && !artifact.metadata.dirty_summary.empty()) ||
      (!artifact.metadata.clean && artifact.metadata.dirty_summary.empty())) {
    throw std::invalid_argument("clean/dirty metadata is inconsistent");
  }
  validate_compatibility(artifact.metadata.compatibility);
  if (artifact.correctness.summary.empty()) {
    throw std::invalid_argument("correctness summary must not be empty");
  }
  if (artifact.aggregation.repetitions !=
          artifact.metadata.compatibility.repetitions ||
      artifact.aggregation.ranks != artifact.metadata.compatibility.ranks ||
      artifact.aggregation.measured_steps !=
          artifact.metadata.compatibility.measured_steps) {
    throw std::invalid_argument(
        "artifact measurement metadata does not match samples");
  }
  const SampleAggregation recomputed = aggregate_samples(
      artifact.aggregation.raw_samples, artifact.aggregation.repetitions,
      artifact.aggregation.ranks);
  require_same_aggregation(artifact.aggregation, recomputed);
  validate_comparison(artifact.comparison);
  for (const CounterMap *counters :
       {&artifact.counters.allocated_bytes,
        &artifact.counters.halo_payload_bytes, &artifact.counters.halo_messages,
        &artifact.counters.collectives,
        &artifact.counters.collective_logical_payload_bytes,
        &artifact.counters.matvec,
        &artifact.counters.preconditioner_applications,
        &artifact.counters.logical_io_bytes}) {
    validate_counter_map(*counters);
  }

  const CompatibilityMetadata &metadata = artifact.metadata.compatibility;
  std::string output;
  output.reserve(2048U + artifact.aggregation.raw_samples.size() * 96U);
  output += "{\"schema_version\":1,\"commit\":";
  append_json_string(output, artifact.metadata.commit);
  output += ",\"working_tree\":{\"clean\":";
  output += artifact.metadata.clean ? "true" : "false";
  output += ",\"dirty_summary\":";
  append_json_string(output, artifact.metadata.dirty_summary);
  output += "},\"case_fingerprint\":";
  append_json_string(output, metadata.problem_fingerprint);
  output += ",\"compiler\":{\"identity\":";
  append_json_string(output, metadata.compiler_identity);
  output += ",\"version\":";
  append_json_string(output, metadata.compiler_version);
  output += ",\"flags\":";
  append_json_string(output, metadata.compiler_flags);
  output += "},\"link_flags\":";
  append_json_string(output, metadata.link_flags);
  output += ",\"build_type\":";
  append_json_string(output, metadata.build_type);
  output += ",\"mpi_implementation\":";
  append_json_string(output, metadata.mpi_identity);
  output += ",\"node_identity\":";
  append_json_string(output, metadata.node_identity);
  output += ",\"hardware_identity\":";
  append_json_string(output, metadata.hardware_identity);
  output += ",\"cpu_affinity\":";
  append_json_string(output, metadata.cpu_affinity);
  output += ",\"rank_placement\":";
  append_json_string(output, metadata.rank_placement);
  output += ",\"ranks\":" + std::to_string(metadata.ranks);
  output += ",\"threads\":" + std::to_string(metadata.threads);
  output += ",\"process_grid\":";
  append_grid(output, metadata.process_grid);
  output += ",\"global_owned_cell_extents\":";
  append_extents(output, metadata.global_owned_cell_extents);
  output += ",\"per_rank_owned_cell_extents\":";
  append_extents(output, metadata.per_rank_owned_cell_extents);
  output += ",\"execution_backend\":";
  append_json_string(output, metadata.execution_backend);
  output += ",\"numerical_tolerance_fingerprint\":";
  append_json_string(output, metadata.numerical_tolerance_contract);
  output += ",\"measurement\":{\"method\":";
  append_json_string(output, metadata.measurement_method);
  output += ",\"warmup_steps\":" + std::to_string(metadata.warmup_steps);
  output += ",\"measured_steps\":" + std::to_string(metadata.measured_steps);
  output += ",\"repetitions\":" + std::to_string(metadata.repetitions) + '}';
  output += ",\"correctness\":{\"passed\":";
  output += artifact.correctness.passed ? "true" : "false";
  output += ",\"summary\":";
  append_json_string(output, artifact.correctness.summary);
  output += "},\"raw_samples\":[";
  bool first = true;
  for (const RawSample &sample : artifact.aggregation.raw_samples) {
    if (!first) {
      output.push_back(',');
    }
    first = false;
    output += "{\"repetition\":" + std::to_string(sample.repetition) +
              ",\"relative_rank\":" + std::to_string(sample.relative_rank) +
              ",\"elapsed_seconds\":";
    append_double(output, sample.elapsed_seconds);
    output +=
        ",\"measured_steps\":" + std::to_string(sample.measured_steps) + '}';
  }
  output += "],\"repetition_maxima\":[";
  first = true;
  for (const RepetitionMaximum &maximum :
       artifact.aggregation.repetition_maxima) {
    if (!first) {
      output.push_back(',');
    }
    first = false;
    output += "{\"repetition\":" + std::to_string(maximum.repetition) +
              ",\"slowest_relative_rank\":" +
              std::to_string(maximum.slowest_relative_rank) +
              ",\"step_seconds\":";
    append_double(output, maximum.step_seconds);
    output.push_back('}');
  }
  output += "],\"median_step_seconds\":";
  append_double(output, artifact.aggregation.median_step_seconds);
  output += ",\"comparison\":{\"mode\":";
  append_json_string(output, mode_name(artifact.comparison.mode));
  output += ",\"status\":";
  append_json_string(output, status_name(artifact.comparison.status));
  output += ",\"reasons\":[";
  first = true;
  for (const std::string &reason : artifact.comparison.reasons) {
    if (!first) {
      output.push_back(',');
    }
    first = false;
    append_json_string(output, reason);
  }
  output += "]},\"exact_counters\":{\"allocated_bytes\":";
  append_counter_map(output, artifact.counters.allocated_bytes);
  output += ",\"halo_payload_bytes\":";
  append_counter_map(output, artifact.counters.halo_payload_bytes);
  output += ",\"halo_messages\":";
  append_counter_map(output, artifact.counters.halo_messages);
  output += ",\"collectives\":";
  append_counter_map(output, artifact.counters.collectives);
  output += ",\"collective_logical_payload_bytes\":";
  append_counter_map(output,
                     artifact.counters.collective_logical_payload_bytes);
  output += ",\"matvec\":";
  append_counter_map(output, artifact.counters.matvec);
  output += ",\"preconditioner_applications\":";
  append_counter_map(output, artifact.counters.preconditioner_applications);
  output += ",\"logical_io_bytes\":";
  append_counter_map(output, artifact.counters.logical_io_bytes);
  output += "}}";
  return output;
}

} // namespace hundun::diagnostics
