// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#include "hundun/v04_app.hpp"

#include <mpi.h>
#include <unistd.h>
#include <array>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  using namespace hundun::v04;
  namespace fs = std::filesystem;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  const bool missing_case = argc > 1 && std::string(argv[1]) == "--missing-case";
  std::array<char, 256U> text{};
  if (rank == 0) {
    const auto root = fs::temp_directory_path() / ("hundun-app-control-" + std::to_string(::getpid()));
    root.string().copy(text.data(), text.size() - 1U);
    if (!missing_case && !ApplicationService::initialize_case_directory(root / "case"))
      MPI_Abort(MPI_COMM_WORLD, 2);
  }
  MPI_Bcast(text.data(), static_cast<int>(text.size()), MPI_CHAR, 0, MPI_COMM_WORLD);
  const fs::path root{text.data()};
  bool passed = ranks >= 2;
  for (int mode = 0; mode < 8; ++mode) {
    ApplicationRunOptions options;
    options.case_root = root / "case";
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    options.run_directory = root / ("output-" + std::to_string(mode));
    options.steps = 1U;
    options.output_interval = options.restart_interval = 0U;
    if (mode == 4 || mode == 5) options.restart_directory = root / "checkpoint";
    if (rank == ranks - 1) {
      if (mode == 0) options.steps = 2U;
      if (mode == 1) options.output_interval = 1U;
      if (mode == 2) options.restart_interval = 1U;
      if (mode == 3) options.restart_directory = root / "checkpoint";
      if (mode == 4) options.restart_history_policy = RestartHistoryPolicy::rebuild_method_history;
      if (mode == 5) options.restart_storage_compatibility = RestartStorageCompatibility::mg_bundle_ghost_v1;
      if (mode == 6) options.initial_state = DriverInitialState{};
      if (mode == 7) options.diagnostics_interval = 1U;
    }
    ApplicationRunReport report;
    const Status status = ApplicationService::run(MPI_COMM_WORLD, options, report);
    const std::uint64_t wire = (std::uint64_t(status.code) << 32U) | status.detail;
    auto low = wire, high = wire;
    MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &high, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    const bool okay = status.code == StatusCode::invalid_case && status.detail == 10506U &&
        low == high && report.failure_phase == ApplicationFailurePhase::input &&
        report.accepted_steps == 0U && report.attempts == 0U && !fs::exists(options.run_directory);
    passed &= okay;
    if (rank == 0) std::cout << "app_control mode=" << mode << " status=" << unsigned(status.code)
        << '/' << status.detail << " cold_rejection=" << okay << '\n';
  }
  if (!missing_case) {
    ApplicationRunOptions options;
    options.case_root = root / "case";
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    options.run_directory = root / "rank-local-time-limits";
    options.steps = 1U;
    options.output_interval = options.restart_interval = 0U;
    const double scale = 1.0 + rank;
    options.time_limits = {scale, scale, scale, scale, scale};
    ApplicationRunReport report;
    const Status status = ApplicationService::run(MPI_COMM_WORLD, options, report);
    passed &= status && report.accepted_steps == 1U;
  }
  int accepted = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &accepted, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0) fs::remove_all(root);
  MPI_Finalize();
  return accepted ? 0 : 1;
}
