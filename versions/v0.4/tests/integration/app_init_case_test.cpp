// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"

#include <mpi.h>
#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

bool run() {
  const fs::path root = fs::temp_directory_path() /
                        ("hundun-v04-init-case-" +
                         std::to_string(::getpid()));
  std::error_code error;
  fs::remove_all(root, error);
  const fs::path case_root = root / "case";
  bool passed = static_cast<bool>(
      ApplicationService::initialize_case_directory(case_root));
  passed &= fs::is_regular_file(case_root / "case.json") &&
            fs::is_regular_file(case_root / "thermophysics.d");
  passed &= ApplicationService::initialize_case_directory(case_root).code ==
            StatusCode::invalid_case;
  CaseValidationReport report;
  passed &= static_cast<bool>(
      ApplicationService::validate(MPI_COMM_SELF, case_root, report));
  passed &= report.case_model != 0U && report.product != 0U &&
            report.summary.sealed && report.summary.pressure_correctors == 2U;
  passed &= static_cast<bool>(ApplicationService::validate_run_directories(
      case_root, root / "run", HUNDUN_V04_SOURCE_ROOT));
  passed &= ApplicationService::validate_run_directories(
                case_root, case_root / "run", HUNDUN_V04_SOURCE_ROOT)
                .code == StatusCode::invalid_case;
  passed &= ApplicationService::validate_run_directories(
                case_root,
                fs::path(HUNDUN_V04_SOURCE_ROOT) / "forbidden-run",
                HUNDUN_V04_SOURCE_ROOT)
                .code == StatusCode::invalid_case;

  ApplicationRunReport run_report;
  ApplicationRunOptions run_options;
  run_options.case_root = case_root;
  run_options.run_directory = root / "run";
  run_options.source_root = HUNDUN_V04_SOURCE_ROOT;
  run_options.steps = 2U;
  run_options.output_interval = 1U;
  run_options.restart_interval = 1U;
  passed &= static_cast<bool>(
      ApplicationService::run(MPI_COMM_SELF, run_options, run_report));
  passed &= run_report.case_model == report.case_model &&
            run_report.product == report.product &&
            run_report.accepted_steps == 2U && run_report.final_time > 0.0 &&
            fs::is_regular_file(root / "run" / "screen.log") &&
            fs::is_regular_file(root / "run" / "monitor.jsonl") &&
            fs::is_regular_file(root / "run" / "evidence.jsonl") &&
            fs::is_regular_file(
                root / "run" / "Visit" /
                "step-00000000000000000002.visit") &&
            fs::exists(root / "run" / "Restart" / "current");
  {
    std::ifstream evidence(root / "run" / "evidence.jsonl",
                           std::ios::binary);
    const std::string text{std::istreambuf_iterator<char>(evidence),
                           std::istreambuf_iterator<char>()};
    passed &= text.find("HUNDUN_V04_EVIDENCE_V8") != std::string::npos &&
              text.find("\"coupling\":\"PISO\"") != std::string::npos &&
              text.find("\"momentum_predictor_passes\":1") !=
                  std::string::npos &&
              text.find("\"candidate_identity\":{") != std::string::npos &&
              text.find("\"previous_committed_time\":") !=
                  std::string::npos &&
              text.find(
                  "\"pressure_solve_contract\":\"continuity_energy_coupled\"") !=
                  std::string::npos &&
              text.find(
                  "\"pressure_energy_refinement_termination\":\"component_residuals_converged\"") !=
                  std::string::npos &&
              text.find("\"pressure_energy_refinement\":[") !=
                  std::string::npos &&
              text.find("\"terminal_physical_audit\":{\"present\":true") !=
                  std::string::npos &&
              text.find("\"committed_convective_cfl\":{") !=
                  std::string::npos &&
              text.find("\"scheme\":\"common_face_afc_v3_owner\"") !=
                  std::string::npos &&
              text.find("\"advective_cfl\":{\"present\":true") !=
                  std::string::npos &&
              text.find("\"final_flux_revision\":0") == std::string::npos &&
              text.find("\"face_flux_revision\":0") == std::string::npos &&
              text.find("\"max_rank_rss_bytes\":0") == std::string::npos &&
              text.find("\"max_node_rss_bytes\":0") == std::string::npos &&
              text.find("\"stages\":[{\"id\":10") != std::string::npos &&
              text.find("\"stages\":[{\"id\":10,\"min_ns\":0") ==
                  std::string::npos;
  }
  ApplicationRunOptions resumed_options = run_options;
  resumed_options.run_directory = root / "run-resumed";
  resumed_options.restart_directory = root / "run" / "Restart";
  resumed_options.steps = 2U;
  ApplicationRunReport resumed_report;
  passed &= static_cast<bool>(ApplicationService::run(
      MPI_COMM_SELF, resumed_options, resumed_report));
  passed &= resumed_report.case_model == report.case_model &&
            resumed_report.product == report.product &&
            resumed_report.accepted_steps == 4U &&
            resumed_report.final_time > run_report.final_time &&
            fs::is_regular_file(
                root / "run-resumed" / "Visit" /
                "step-00000000000000000004.visit");
  {
    std::ifstream evidence(root / "run-resumed" / "evidence.jsonl",
                           std::ios::binary);
    const std::string text{std::istreambuf_iterator<char>(evidence),
                           std::istreambuf_iterator<char>()};
    passed &= text.find("HUNDUN_V04_EVIDENCE_V8") != std::string::npos &&
              text.find("\"coupling\":\"PISO\"") != std::string::npos &&
              text.find("\"momentum_predictor_passes\":1") !=
                  std::string::npos &&
              text.find("\"candidate_identity\":{") != std::string::npos &&
              text.find("\"step\":3") != std::string::npos &&
              text.find("\"bdf_order\":2") != std::string::npos &&
              text.find("\"restart_recovery\":true") ==
                  std::string::npos &&
              text.find("\"restart_recovery\":false") !=
                  std::string::npos &&
              text.find(
                  "\"pressure_solve_contract\":\"continuity_energy_coupled\"") !=
                  std::string::npos &&
              text.find(
                  "\"pressure_energy_refinement_termination\":\"component_residuals_converged\"") !=
                  std::string::npos &&
              text.find("\"terminal_physical_audit\":{\"present\":true") !=
                  std::string::npos &&
              text.find("\"committed_convective_cfl\":{") !=
                  std::string::npos &&
              text.find("\"scheme\":\"common_face_afc_v3_owner\"") !=
                  std::string::npos &&
              text.find("\"advective_cfl\":{\"present\":true") !=
                  std::string::npos &&
              text.find("\"face_flux_revision\":0") == std::string::npos &&
              text.find("\"final_flux_revision\":0") == std::string::npos;
  }
  ApplicationRunOptions benchmark_options = run_options;
  benchmark_options.run_directory = root / "run-no-serialized-output";
  benchmark_options.steps = 2U;
  benchmark_options.output_interval = 0U;
  benchmark_options.restart_interval = 0U;
  ApplicationRunReport benchmark_report;
  passed &= static_cast<bool>(ApplicationService::run(
      MPI_COMM_SELF, benchmark_options, benchmark_report));
  passed &= benchmark_report.accepted_steps == 2U &&
            fs::is_regular_file(benchmark_options.run_directory /
                                "evidence.jsonl") &&
            !fs::exists(benchmark_options.run_directory / "Visit") &&
            !fs::exists(benchmark_options.run_directory / "Restart") &&
            !fs::exists(benchmark_options.run_directory / "screen.log") &&
            !fs::exists(benchmark_options.run_directory / "monitor.jsonl");
  ApplicationRunReport unchanged_run;
  unchanged_run.product = UINT64_C(0xcafef00d);
  ApplicationRunOptions invalid_run = run_options;
  invalid_run.run_directory = case_root / "inside";
  passed &= !ApplicationService::run(MPI_COMM_SELF, invalid_run,
                                     unchanged_run) &&
            unchanged_run.product == UINT64_C(0xcafef00d);

  for (fs::recursive_directory_iterator iterator(HUNDUN_V04_SOURCE_ROOT, error),
       end;
       !error && iterator != end; iterator.increment(error)) {
    if (iterator->is_regular_file(error) &&
        iterator->path().filename() == "case.json") {
      passed = false;
      break;
    }
  }
  passed &= !error;
  CaseValidationReport unchanged;
  unchanged.product = UINT64_C(0xdeadbeef);
  passed &= !ApplicationService::validate(MPI_COMM_SELF, root / "missing",
                                          unchanged) &&
            unchanged.product == UINT64_C(0xdeadbeef);
  fs::remove_all(root, error);
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed = run();
  if (!passed) std::cerr << "init-case/product validation failure\n";
  MPI_Finalize();
  return passed ? 0 : 1;
}
