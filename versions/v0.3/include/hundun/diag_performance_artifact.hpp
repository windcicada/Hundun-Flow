// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace hundun::diagnostics {

struct RawSample {
  int repetition{};
  int relative_rank{};
  double elapsed_seconds{};
  std::uint64_t measured_steps{};
};

struct RepetitionMaximum {
  int repetition{};
  int slowest_relative_rank{};
  double step_seconds{};
};

struct SampleAggregation {
  int repetitions{};
  int ranks{};
  std::uint64_t measured_steps{};
  std::vector<RawSample> raw_samples;
  std::vector<RepetitionMaximum> repetition_maxima;
  double median_step_seconds{};
};

SampleAggregation aggregate_samples(std::vector<RawSample> samples,
                                    int repetitions, int ranks);

using ClockFunction = std::function<double()>;
using PhaseFunction = std::function<void()>;

double measure_elapsed_phase(const ClockFunction &clock,
                             const PhaseFunction &phase);

struct ProcessGrid {
  int x{};
  int y{};
  int z{};
};

struct CellExtents {
  std::uint64_t x{};
  std::uint64_t y{};
  std::uint64_t z{};
};

struct CompatibilityMetadata {
  std::string hardware_identity;
  std::string node_identity;
  std::string mpi_identity;
  std::string compiler_identity;
  std::string compiler_version;
  std::string compiler_flags;
  std::string link_flags;
  std::string build_type;
  std::string cpu_affinity;
  std::string rank_placement;
  std::string problem_fingerprint;
  std::string numerical_tolerance_contract;
  std::string measurement_method;
  std::uint64_t warmup_steps{};
  std::uint64_t measured_steps{};
  int repetitions{};
  std::string execution_backend;
  int ranks{};
  int threads{};
  ProcessGrid process_grid;
  CellExtents global_owned_cell_extents;
  CellExtents per_rank_owned_cell_extents;
};

enum class ComparisonMode { identical, strong_scaling, weak_scaling };
enum class ComparisonStatus { comparable, incomparable };

struct ComparisonResult {
  ComparisonMode mode{ComparisonMode::identical};
  ComparisonStatus status{ComparisonStatus::comparable};
  std::vector<std::string> reasons;
};

ComparisonResult compare_compatibility(const CompatibilityMetadata &baseline,
                                       const CompatibilityMetadata &candidate,
                                       ComparisonMode mode);

double strong_scaling_speedup(double single_rank_step_seconds,
                              double parallel_step_seconds);
double strong_scaling_efficiency(double single_rank_step_seconds,
                                 double parallel_step_seconds, int ranks);
double weak_scaling_efficiency(double single_rank_step_seconds,
                               double parallel_step_seconds);

struct ArtifactMetadata {
  std::string commit;
  bool clean{true};
  std::string dirty_summary;
  CompatibilityMetadata compatibility;
  std::string tree_fingerprint;
  std::string binary_fingerprint;
  std::string profile;
  std::string geometry_fingerprint;
  std::string cpuset;
  int thread_budget{};
};

ComparisonResult compare_artifact_metadata(const ArtifactMetadata& baseline,
                                           const ArtifactMetadata& candidate,
                                           ComparisonMode mode);

struct CorrectnessResult {
  bool passed{};
  std::string summary;
};

using CounterMap = std::map<std::string, std::uint64_t>;

struct ExactCounterMaps {
  CounterMap allocated_bytes;
  CounterMap halo_payload_bytes;
  CounterMap halo_messages;
  CounterMap collectives;
  CounterMap collective_logical_payload_bytes;
  CounterMap matvec;
  CounterMap preconditioner_applications;
  CounterMap logical_io_bytes;
  CounterMap algorithmic_work;
};

struct Artifact {
  int schema_version{1};
  ArtifactMetadata metadata;
  CorrectnessResult correctness;
  SampleAggregation aggregation;
  ComparisonResult comparison;
  ExactCounterMaps counters;
};

std::string to_json(const Artifact &artifact);

} // namespace hundun::diagnostics
