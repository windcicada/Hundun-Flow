// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>
#include <unistd.h>

#include <array>
#include <filesystem>
#include <fstream>
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
  {
    ApplicationRunOptions invalid;
    ApplicationRunReport reused;
    reused.accepted_steps = 99U;
    reused.final_time = 7.0;
    reused.attempts = 4U;
    reused.failed_stage = 50U;
    reused.failure_phase = ApplicationFailurePhase::advance;
    const Status status =
        ApplicationService::run(MPI_COMM_WORLD, invalid, reused);
    passed = passed && !status && reused.failure.code == status.code &&
             reused.failure.detail == status.detail &&
             reused.accepted_steps == 0U && reused.final_time == 0.0 &&
             reused.attempts == 0U && reused.failed_stage == 0U &&
             reused.failure_phase == ApplicationFailurePhase::input;
    if (rank == 0)
      std::cout << "app reused early report reset=" << passed << '\n';
  }
  for (int target = 0; target<ranks; target += ranks> 1 ? ranks - 1 : 1) {
    ApplicationRunOptions invalid;
    invalid.case_root = root / "case";
    invalid.run_directory = root / "invalid-local-options";
    invalid.source_root = HUNDUN_V04_SOURCE_ROOT;
    invalid.steps = rank == target ? 0U : 1U;
    ApplicationRunReport report;
    const Status status =
        ApplicationService::run(MPI_COMM_WORLD, invalid, report);
    passed &= status.code == StatusCode::invalid_case &&
              report.failure_phase == ApplicationFailurePhase::input &&
              report.accepted_steps == 0U && report.attempts == 0U;
  }
  for (int target = 0; target<ranks; target += ranks> 1 ? ranks - 1 : 1) {
    for (const auto phase :
         {ApplicationFailurePhase::initialize, ApplicationFailurePhase::visit,
          ApplicationFailurePhase::screen, ApplicationFailurePhase::monitor,
          ApplicationFailurePhase::restart, ApplicationFailurePhase::resources,
          ApplicationFailurePhase::evidence}) {
      if (!passed) break;
      ApplicationRunOptions options;
      options.case_root = root / "case";
      options.run_directory =
          root / ("phase-" + std::to_string(static_cast<int>(phase)) +
                  "-rank-" + std::to_string(target));
      options.source_root = HUNDUN_V04_SOURCE_ROOT;
      options.steps = 1U;
      ApplicationRunReport report;
      detail::arm_application_local_allocation_failure_for_test(phase, target);
      const Status status =
          ApplicationService::run(MPI_COMM_WORLD, options, report);
      const std::uint64_t packed =
          (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
      std::uint64_t low = 0U, high = 0U;
      MPI_Allreduce(&packed, &low, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&packed, &high, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
      const bool setup = phase == ApplicationFailurePhase::initialize;
      int okay =
          status.code == StatusCode::allocation_failure && low == high &&
                  report.failure.code == status.code &&
                  report.failure.detail == status.detail &&
                  report.failure_phase == phase &&
                  report.accepted_steps == (setup ? 0U : 1U) &&
                  (setup ? report.final_time == 0.0 && report.attempts == 0U
                         : report.final_time > 0.0 &&
                               report.step_completion.outcome) &&
                  report.failed_stage == 0U
              ? 1
              : 0;
      MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      passed = okay != 0;
      if (rank == 0)
        std::cout << "app local failure phase=" << static_cast<int>(phase)
                  << " target=" << target << " consistent=" << passed << '\n';
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  {
    ApplicationRunOptions options;
    options.case_root = root / "case";
    options.run_directory = root / "real-io-failure";
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    options.steps = 1U;
    if (rank == 0) {
      fs::create_directories(options.run_directory);
      std::ofstream occupied(options.run_directory / "Visit");
      occupied << "an existing regular file blocks directory creation\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ApplicationRunReport report;
    const Status status =
        ApplicationService::run(MPI_COMM_WORLD, options, report);
    passed &=
        status.code == StatusCode::io_failure &&
        report.failure_phase == ApplicationFailurePhase::visit &&
        report.accepted_steps == 1U && report.final_time > 0.0 &&
        report.failed_stage == 0U && report.io_failure.valid &&
        report.io_failure.rank == 0 && report.io_failure.system_error != 0 &&
        report.io_failure.operation == IoFailureOperation::create_directory &&
        fs::path(report.io_failure.path.data()) ==
            options.run_directory / "Visit";
    if (rank == 0)
      std::cout << "app real IO errno=" << report.io_failure.system_error
                << " accepted_steps=" << report.accepted_steps << '\n';
  }
  {
    ApplicationRunOptions options;
    options.case_root = root / "case";
    options.run_directory = root / "timing-success";
    options.source_root = HUNDUN_V04_SOURCE_ROOT;
    options.steps = 1U;
    ApplicationRunReport report;
    const Status status =
        ApplicationService::run(MPI_COMM_WORLD, options, report);
    std::uint64_t accounted = 0U;
    for (auto value : report.local_phase_nanoseconds) accounted += value;
    passed &= status && report.accepted_steps == 1U &&
              accounted == report.local_run_nanoseconds && accounted > 0U;
    for (auto value : report.local_phase_nanoseconds) passed &= value > 0U;
    if (rank == 0)
      std::cout << "app local wall accounted=" << accounted
                << " run=" << report.local_run_nanoseconds << '\n';
  }
  int final_passed = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &final_passed, 1, MPI_INT, MPI_MIN,
                MPI_COMM_WORLD);
  passed = final_passed != 0;
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) fs::remove_all(root);
  MPI_Finalize();
  return passed ? 0 : 1;
}
