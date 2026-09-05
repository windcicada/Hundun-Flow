// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>
#include <unistd.h>

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "hundun/v04_app.hpp"

namespace {
// Observe the real MPI resource boundary; all communication still uses MPI.
// Sampling must not create communicators in proportion to accepted steps.
bool observing = false;
std::array<MPI_Comm, 32U> live{};
std::uint64_t creations = 0U, releases = 0U, split_nanoseconds = 0U;
bool tracking_overflow = false;
using Clock = std::chrono::steady_clock;
}  // namespace

extern "C" int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key,
                                   MPI_Info info, MPI_Comm* out) {
  const auto begin = Clock::now();
  const int status = PMPI_Comm_split_type(comm, split_type, key, info, out);
  if (observing && status == MPI_SUCCESS && *out != MPI_COMM_NULL) {
    ++creations;
    split_nanoseconds += static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() -
                                                             begin)
            .count());
    for (auto& slot : live) {
      if (slot == MPI_COMM_NULL) {
        slot = *out;
        return status;
      }
    }
    tracking_overflow = true;
  }
  return status;
}

extern "C" int MPI_Comm_free(MPI_Comm* comm) {
  const MPI_Comm before = *comm;
  const int status = PMPI_Comm_free(comm);
  if (observing && status == MPI_SUCCESS) {
    for (auto& slot : live) {
      if (slot == before) {
        slot = MPI_COMM_NULL;
        ++releases;
        break;
      }
    }
  }
  return status;
}

int main(int argc, char** argv) {
  using namespace hundun::v04;
  namespace fs = std::filesystem;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<char, 256U> root_text{};
  int ready = 1;
  if (rank == 0) {
    const fs::path root =
        fs::temp_directory_path() /
        ("hundun-app-resources-" + std::to_string(::getpid()));
    root.string().copy(root_text.data(), root_text.size() - 1U);
    ready =
        ApplicationService::initialize_case_directory(root / "case") ? 1 : 0;
  }
  MPI_Bcast(root_text.data(), static_cast<int>(root_text.size()), MPI_CHAR, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(&ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
  const fs::path root{root_text.data()};
  bool passed = ready != 0;
  std::uint64_t one_step_creations = 0U;
  for (const std::uint64_t steps : {1U, 4U, 2U}) {
    if (!passed) break;
    ApplicationRunOptions options;
    options.case_root = root / "case";
    options.run_directory = root / ("steps-" + std::to_string(steps));
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    options.steps = steps;
    options.output_interval = 0U;
    options.restart_interval = 0U;
    const bool fail_evidence = steps == 2U;
    if (fail_evidence && rank == 0)
      fs::create_directories(options.run_directory / "evidence.jsonl");
    MPI_Barrier(MPI_COMM_WORLD);
    ApplicationRunReport report;
    live.fill(MPI_COMM_NULL);
    creations = releases = split_nanoseconds = 0U;
    tracking_overflow = false;
    observing = true;
    const auto begin = Clock::now();
    const Status status =
        ApplicationService::run(MPI_COMM_WORLD, options, report);
    const double seconds =
        std::chrono::duration<double>(Clock::now() - begin).count();
    observing = false;
    if (steps == 1U) one_step_creations = creations;
    const bool expected_outcome =
        fail_evidence
            ? !status && report.accepted_steps == 1U &&
                  report.failure_phase == ApplicationFailurePhase::evidence
            : status && report.accepted_steps == steps;
    int okay = expected_outcome && report.final_time > 0.0 &&
                       !tracking_overflow && creations == releases &&
                       creations <= one_step_creations
                   ? 1
                   : 0;
    if (rank == 0 && !fail_evidence) {
      std::ifstream file(options.run_directory / "evidence.jsonl");
      std::string line;
      std::uint64_t records = 0U;
      while (std::getline(file, line)) {
        ++records;
        okay &= line.find("\"max_rank_rss_bytes\":") != std::string::npos &&
                line.find("\"max_node_rss_bytes\":") != std::string::npos &&
                line.find("\"max_rank_rss_bytes\":0") == std::string::npos &&
                line.find("\"max_node_rss_bytes\":0") == std::string::npos;
      }
      okay &= records == steps;
      std::cout << "resource lifecycle steps=" << steps
                << " shared_created=" << creations
                << " shared_freed=" << releases
                << " split_ns=" << split_nanoseconds << " run_s=" << seconds
                << " evidence_records=" << records << '\n';
    }
    if (rank == 0 && fail_evidence)
      std::cout << "resource lifecycle output failure shared_created="
                << creations << " shared_freed=" << releases
                << " checked=" << okay << '\n';
    MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    passed = okay != 0;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) fs::remove_all(root);
  if (rank == 0 && !passed)
    std::cerr << "FAIL: resource sampling must retain RSS evidence without "
                 "per-step communicator creation or unreleased communicators\n";
  MPI_Finalize();
  return passed ? 0 : 1;
}
