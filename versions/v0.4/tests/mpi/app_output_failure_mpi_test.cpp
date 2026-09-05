// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>
#include <unistd.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <string>

#include "app_driver_detail.hpp"
#include "hundun/v04_app.hpp"

int main(int argc, char** argv) {
  using namespace hundun::v04;
  namespace fs = std::filesystem;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  std::array<char, 256U> root_text{};
  int ready = 1;
  if (rank == 0) {
    const fs::path root =
        fs::temp_directory_path() /
        ("hundun-app-output-failure-" + std::to_string(::getpid()));
    root.string().copy(root_text.data(), root_text.size() - 1U);
    ready =
        ApplicationService::initialize_case_directory(root / "case") ? 1 : 0;
  }
  MPI_Bcast(root_text.data(), static_cast<int>(root_text.size()), MPI_CHAR, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(&ready, 1, MPI_INT, 0, MPI_COMM_WORLD);
  const fs::path root{root_text.data()};
  bool passed = ready != 0;
  for (const auto phase :
       {ApplicationFailurePhase::visit, ApplicationFailurePhase::screen,
        ApplicationFailurePhase::monitor, ApplicationFailurePhase::restart,
        ApplicationFailurePhase::resources,
        ApplicationFailurePhase::evidence}) {
    if (!passed) break;
    ApplicationRunOptions options;
    options.case_root = root / "case";
    options.run_directory =
        root / ("phase-" + std::to_string(static_cast<int>(phase)));
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    options.steps = 1U;
    ApplicationRunReport report;
    detail::arm_application_local_allocation_failure_for_test(phase, ranks - 1);
    const Status status =
        ApplicationService::run(MPI_COMM_WORLD, options, report);
    const std::uint64_t packed =
        (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
    std::uint64_t low = 0U, high = 0U;
    MPI_Allreduce(&packed, &low, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&packed, &high, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    int okay = status.code == StatusCode::allocation_failure && low == high &&
                       report.failure.code == status.code &&
                       report.failure.detail == status.detail &&
                       report.failure_phase == phase &&
                       report.accepted_steps == 1U && report.final_time > 0.0 &&
                       report.failed_stage == 0U
                   ? 1
                   : 0;
    MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    passed = okay != 0;
    if (rank == 0)
      std::cout << "app local failure phase=" << static_cast<int>(phase)
                << " consistent=" << passed << '\n';
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) fs::remove_all(root);
  MPI_Finalize();
  return passed ? 0 : 1;
}
