// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_performance_build.hpp"
#include "hundun/diag_performance_artifact.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/stage3_performance_evidence.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using hundun::diagnostics::Stage3PerformanceCounters;

std::string required_environment(const char *name) {
  const char *value = std::getenv(name);
  if (value == nullptr || *value == '\0')
    throw std::runtime_error(std::string("missing required environment: ") +
                             name);
  return value;
}

std::string optional_environment(const char *name,
                                 std::string fallback) {
  const char *value = std::getenv(name);
  return value == nullptr || *value == '\0' ? std::move(fallback)
                                             : std::string(value);
}

std::string read_first_matching_line(const std::filesystem::path &path,
                                     std::string_view prefix) {
  std::ifstream input(path);
  std::string line;
  while (std::getline(input, line)) {
    if (line.compare(0U, prefix.size(), prefix) != 0)
      continue;
    const auto begin = line.find_first_not_of(" \t:", prefix.size());
    return begin == std::string::npos ? std::string{} : line.substr(begin);
  }
  return {};
}

std::string mpi_identity() {
  std::array<char, MPI_MAX_LIBRARY_VERSION_STRING> buffer{};
  int length = 0;
  if (MPI_Get_library_version(buffer.data(), &length) != MPI_SUCCESS ||
      length <= 0)
    throw std::runtime_error("unable to identify the MPI implementation");
  std::string result(buffer.data(), static_cast<std::size_t>(length));
  for (char &character : result)
    if (character == '\n' || character == '\r' || character == '\t')
      character = ' ';
  while (!result.empty() &&
         (result.back() == ' ' || result.back() == '\0'))
    result.pop_back();
  return result;
}

struct Placement final {
  std::string nodes;
  std::string ranks;
};

Placement placement(const hundun::runtime::MpiContext &mpi) {
  std::array<char, MPI_MAX_PROCESSOR_NAME> local{};
  int length = 0;
  if (MPI_Get_processor_name(local.data(), &length) != MPI_SUCCESS ||
      length <= 0 || length >= static_cast<int>(local.size()))
    throw std::runtime_error("unable to identify MPI processor placement");
  local[static_cast<std::size_t>(length)] = '\0';
  std::vector<char> gathered(local.size() * static_cast<std::size_t>(mpi.size()));
  if (MPI_Allgather(local.data(), static_cast<int>(local.size()), MPI_CHAR,
                    gathered.data(), static_cast<int>(local.size()), MPI_CHAR,
                    mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to gather MPI processor placement");
  std::set<std::string> nodes;
  std::ostringstream ranks;
  for (int rank = 0; rank < mpi.size(); ++rank) {
    const char *name = gathered.data() +
                       static_cast<std::size_t>(rank) * local.size();
    nodes.emplace(name);
    if (rank != 0)
      ranks << ',';
    ranks << rank << ':' << name;
  }
  std::ostringstream node_list;
  for (auto iterator = nodes.begin(); iterator != nodes.end(); ++iterator) {
    if (iterator != nodes.begin())
      node_list << ',';
    node_list << *iterator;
  }
  return {node_list.str(), ranks.str()};
}

std::string tagged_hex(std::string_view tag, std::uint64_t value) {
  std::ostringstream result;
  result << tag << ':' << std::hex << std::setfill('0') << std::setw(16)
         << value;
  return result.str();
}

int thread_budget() {
  const std::string text = optional_environment("HUNDUN_STAGE3_THREAD_BUDGET",
                                                "1");
  std::size_t consumed = 0U;
  const long value = std::stol(text, &consumed);
  if (consumed != text.size() || value <= 0 || value > INT_MAX)
    throw std::runtime_error("invalid HUNDUN_STAGE3_THREAD_BUDGET");
  return static_cast<int>(value);
}

Stage3PerformanceCounters artifact_counters(
    const hundun::test::Stage3PerformanceEvidence &evidence) {
  auto result = evidence.measured_delta;
  result.init_surface_triangles =
      evidence.after_initialization.init_surface_triangles;
  result.init_query_closest_calls =
      evidence.after_initialization.init_query_closest_calls;
  result.init_query_segment_calls =
      evidence.after_initialization.init_query_segment_calls;
  result.init_classification_cells =
      evidence.after_initialization.init_classification_cells;
  result.init_ghost_qr_plans =
      evidence.after_initialization.init_ghost_qr_plans;
  result.init_ghost_rejected_plans =
      evidence.after_initialization.init_ghost_rejected_plans;
  result.init_ghost_donor_references =
      evidence.after_initialization.init_ghost_donor_references;
  result.init_wall_points = evidence.after_initialization.init_wall_points;
  result.checkpoint_logical_io_bytes =
      evidence.checkpoint_delta.checkpoint_logical_io_bytes;
  return result;
}

void write_formal_artifact(
    const hundun::runtime::MpiContext &mpi,
    const hundun::test::Stage3PerformanceRun &run,
    const hundun::test::Stage3PerformanceEvidence &evidence) {
  std::vector<double> elapsed(static_cast<std::size_t>(mpi.size()));
  if (MPI_Allgather(&evidence.elapsed_seconds, 1, MPI_DOUBLE, elapsed.data(), 1,
                    MPI_DOUBLE, mpi.comm()) != MPI_SUCCESS)
    throw std::runtime_error("unable to gather Stage 3 performance samples");

  const auto grid = hundun::test::stage3_performance_process_grid(mpi.size());
  const auto where = placement(mpi);
  const std::string cpuset = optional_environment(
      "HUNDUN_STAGE3_CPUSET",
      read_first_matching_line("/proc/self/status", "Cpus_allowed_list"));
  const std::string hardware = optional_environment(
      "HUNDUN_STAGE3_HARDWARE_IDENTITY",
      read_first_matching_line("/proc/cpuinfo", "model name"));
  if (cpuset.empty() || hardware.empty())
    throw std::runtime_error("unable to identify Stage 3 performance host");

  hundun::diagnostics::CompatibilityMetadata compatibility;
  compatibility.hardware_identity = hardware;
  compatibility.node_identity = where.nodes;
  compatibility.mpi_identity = mpi_identity();
#if defined(__clang__)
  compatibility.compiler_identity = "clang";
#elif defined(__GNUC__)
  compatibility.compiler_identity = "gcc";
#else
  compatibility.compiler_identity = "unknown-cxx";
#endif
  compatibility.compiler_version = __VERSION__;
  compatibility.compiler_flags = std::string(
      hundun::application::detail::performance_compiler_flags);
  compatibility.link_flags =
      std::string(hundun::application::detail::performance_link_flags);
  compatibility.build_type =
      std::string(hundun::application::detail::performance_build_type);
  compatibility.cpu_affinity = cpuset;
  compatibility.rank_placement = where.ranks;
  compatibility.problem_fingerprint =
      "stage3-exact-counters-n" + std::to_string(run.cells) +
      "-constant-static-ibm-wale";
  compatibility.numerical_tolerance_contract = "stage3-task11-frozen-v2";
  compatibility.measurement_method = "mpi-wtime-v1";
  compatibility.warmup_steps = run.warmup_steps;
  compatibility.measured_steps = run.measured_steps;
  compatibility.repetitions = run.repetitions;
  compatibility.execution_backend = "cpu_reference";
  compatibility.ranks = mpi.size();
  compatibility.threads = thread_budget();
  compatibility.process_grid = {grid.x, grid.y, grid.z};
  compatibility.global_owned_cell_extents = {
      static_cast<std::uint64_t>(run.cells),
      static_cast<std::uint64_t>(run.cells),
      static_cast<std::uint64_t>(run.cells)};
  compatibility.per_rank_owned_cell_extents = {
      static_cast<std::uint64_t>(run.cells / grid.x),
      static_cast<std::uint64_t>(run.cells / grid.y),
      static_cast<std::uint64_t>(run.cells / grid.z)};

  hundun::diagnostics::Artifact artifact;
  artifact.schema_version = 2;
  artifact.metadata = {
      std::string(hundun::application::detail::performance_source_commit),
      hundun::application::detail::performance_source_clean,
      std::string(
          hundun::application::detail::performance_source_dirty_summary),
      std::move(compatibility),
      required_environment("HUNDUN_STAGE3_TREE_FINGERPRINT"),
      required_environment("HUNDUN_STAGE3_BINARY_FINGERPRINT"),
      "constant-static-ibm-wale",
      tagged_hex("surface-domain-crc64",
                 evidence.surface_fingerprint ^
                     evidence.classification_fingerprint),
      cpuset,
      thread_budget()};
  artifact.correctness = {true, "stage3-exact-counters:no-retry:pass"};
  std::vector<hundun::diagnostics::RawSample> samples;
  samples.reserve(elapsed.size());
  for (int rank = 0; rank < mpi.size(); ++rank)
    samples.push_back({0, rank, elapsed[static_cast<std::size_t>(rank)],
                       run.measured_steps});
  artifact.aggregation = hundun::diagnostics::aggregate_samples(
      std::move(samples), run.repetitions, mpi.size());
  artifact.comparison = hundun::diagnostics::compare_artifact_metadata(
      artifact.metadata, artifact.metadata,
      hundun::diagnostics::ComparisonMode::identical);
  artifact.counters.algorithmic_work =
      hundun::diagnostics::stage3_algorithmic_work_map(
          artifact_counters(evidence));
  const std::string json = hundun::diagnostics::to_json(artifact);

  if (mpi.rank() == 0) {
    const std::filesystem::path directory =
        required_environment("HUNDUN_STAGE3_EVIDENCE_DIR");
    std::filesystem::create_directories(directory);
    const auto final = directory /
                       ("stage3-exact-counters-n" +
                        std::to_string(run.cells) + "-r" +
                        std::to_string(mpi.size()) + ".v2.json");
    const auto temporary = final.string() + ".tmp";
    {
      std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
      output.exceptions(std::ios::failbit | std::ios::badbit);
      output << json << '\n';
    }
    std::filesystem::rename(temporary, final);
  }
  mpi.barrier();
}

void verify_evidence(
    const hundun::test::Stage3PerformanceRun &run,
    const hundun::test::Stage3PerformanceEvidence &evidence) {
  HUNDUN_CHECK(evidence.committed);
  HUNDUN_CHECK(evidence.pressure_corrector_count == 2U);
  HUNDUN_CHECK(evidence.after_initialization.init_surface_triangles > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_query_closest_calls > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_query_segment_calls > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_classification_cells > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_ghost_qr_plans > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_ghost_rejected_plans > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_ghost_donor_references > 0U);
  HUNDUN_CHECK(evidence.after_initialization.init_wall_points > 0U);
  HUNDUN_CHECK(evidence.measured_delta.step_ghost_constraints > 0U);
  HUNDUN_CHECK(evidence.measured_delta.step_lfp_transforms > 0U);
  HUNDUN_CHECK(evidence.measured_delta.step_immersed_rows > 0U);
  HUNDUN_CHECK(evidence.measured_delta.step_pressure_wall_constraints ==
               2U * run.measured_steps * evidence.immersed_link_count);
  HUNDUN_CHECK(evidence.measured_delta.step_wall_quadrature_evaluations ==
               run.measured_steps * evidence.local_wall_point_count);
  HUNDUN_CHECK(evidence.measured_delta.step_force_reductions ==
               3U * run.measured_steps);
  HUNDUN_CHECK(evidence.measured_delta.step_wale_gradient_cells ==
               run.measured_steps * evidence.owned_active_cell_count);
  HUNDUN_CHECK(evidence.measured_delta.step_wale_evaluations ==
               run.measured_steps);
  HUNDUN_CHECK(evidence.measured_delta.init_surface_triangles == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_query_closest_calls == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_query_segment_calls == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_classification_cells == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_ghost_qr_plans == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_ghost_rejected_plans == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_ghost_donor_references == 0U);
  HUNDUN_CHECK(evidence.measured_delta.init_wall_points == 0U);
  HUNDUN_CHECK(evidence.checkpoint_delta.checkpoint_logical_io_bytes > 0U);
  if (run.inject_failed_attempt) {
    HUNDUN_CHECK(evidence.failed_attempt_rolled_back);
    HUNDUN_CHECK(evidence.failed_attempt_delta.step_ghost_constraints > 0U);
    HUNDUN_CHECK(evidence.failed_attempt_delta.step_wale_evaluations == 1U);
  } else {
    HUNDUN_CHECK(!evidence.failed_attempt_rolled_back);
    HUNDUN_CHECK(
        hundun::diagnostics::stage3_algorithmic_work_map(
            evidence.failed_attempt_delta) ==
        hundun::diagnostics::stage3_algorithmic_work_map({}));
  }
}

hundun::test::Stage3PerformanceRun parse_run(int argc, char **argv) {
  if (argc == 1)
    return {};
  if (argc != 6 || std::string_view(argv[1]) != "formal")
    throw std::invalid_argument(
        "usage: test_stage3_performance [formal cells warmup measured repetitions]");
  hundun::test::Stage3PerformanceRun result;
  result.cells = std::stoi(argv[2]);
  result.warmup_steps = static_cast<std::uint64_t>(std::stoull(argv[3]));
  result.measured_steps = static_cast<std::uint64_t>(std::stoull(argv[4]));
  result.repetitions = std::stoi(argv[5]);
  result.inject_failed_attempt = false;
  if (result.cells != 24 || result.warmup_steps != 2U ||
      result.measured_steps != 3U || result.repetitions != 1)
    throw std::invalid_argument("formal Stage 3 performance selector is frozen");
  return result;
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const auto run = parse_run(argc, argv);
    const auto evidence =
        hundun::test::run_stage3_performance_evidence(mpi, run);
    verify_evidence(run, evidence);
    if (!run.inject_failed_attempt)
      write_formal_artifact(mpi, run, evidence);
  });
}
