// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"

#include <mpi.h>
#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <set>

namespace {

namespace fs = std::filesystem;
using namespace hundun::v04;

bool run() {
  const auto source_cases = [] {
    std::set<fs::path> paths;
    for (const auto& entry : fs::recursive_directory_iterator(HUNDUN_V04_SOURCE_ROOT))
      if (entry.is_regular_file() && entry.path().filename() == "case.json")
        paths.insert(entry.path());
    return paths;
  };
  const auto original_source_cases = source_cases();
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
  passed &= static_cast<bool>(ApplicationService::validate_run_directories(
      case_root, root / "case-sibling", HUNDUN_V04_SOURCE_ROOT));
  passed &= ApplicationService::validate_run_directories(
                case_root, root / "case" / "nested", HUNDUN_V04_SOURCE_ROOT)
                .code == StatusCode::invalid_case;
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
  if (run_report.piso.pressure_solve_calls != 2U ||
      run_report.momentum_predictor_solve.predictor_passes != 1U ||
      run_report.pressure_energy_globalization.trajectory_count == 0U) {
    std::cerr << "FAIL: successful ApplicationService run preserves the last "
                 "accepted step diagnostics; pressure_solve_calls="
              << static_cast<unsigned>(run_report.piso.pressure_solve_calls)
              << '\n';
    passed = false;
  }
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
    std::ifstream monitor(root / "run" / "monitor.jsonl", std::ios::binary);
    const std::string text{std::istreambuf_iterator<char>(monitor),
                           std::istreambuf_iterator<char>()};
    passed &= text.find("\"candidate_baseline_evaluations\":2") !=
                  std::string::npos &&
              text.find("\"candidate_incomplete_evaluations\":0") !=
                  std::string::npos;
  }
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
  // Independent diagnostics never require full-field Visit output.
  {
    ApplicationRunOptions diagnostics = run_options;
    diagnostics.run_directory = root / "run-diagnostics";
    diagnostics.output_interval = diagnostics.restart_interval = 0U;
    diagnostics.diagnostics_interval = 1U;
    ApplicationRunReport diagnostic_report;
    passed &= static_cast<bool>(ApplicationService::run(
        MPI_COMM_SELF, diagnostics, diagnostic_report));
    std::ifstream input(diagnostics.run_directory / "diagnostics.jsonl");
    const std::string text{std::istreambuf_iterator<char>(input),
                           std::istreambuf_iterator<char>()};
    const bool ledger = std::count(text.begin(), text.end(), '\n') == 2 &&
        text.find("HUNDUN_V04_DEVELOPMENT_DIAGNOSTICS_V1") != std::string::npos &&
        text.find("\"mass_balance_defect_kg_s\":") != std::string::npos &&
        text.find("\"cumulative_energy_defect_J\":") != std::string::npos &&
        text.find("\"statistics_eligible\":false") != std::string::npos &&
        !fs::exists(diagnostics.run_directory / "Visit") &&
        !fs::exists(diagnostics.run_directory / "monitor.jsonl") &&
        !fs::exists(root / "run" / "diagnostics.jsonl");
    if (!ledger) std::cerr << "FAIL: optional committed conservation ledger without Visit\n";
    passed &= ledger;
  }
  // The ordinary application must distinguish an exact continuation from an
  // explicitly requested method recovery, even for an intact current image.
  ApplicationRunOptions recovery_options = resumed_options;
  recovery_options.run_directory = root / "run-method-recovery";
  recovery_options.restart_directory = resumed_options.run_directory / "Restart";
  recovery_options.restart_history_policy = RestartHistoryPolicy::rebuild_method_history;
  ApplicationRunReport recovery_report;
  const auto recovery_status = ApplicationService::run(
      MPI_COMM_SELF, recovery_options, recovery_report);
  std::ifstream recovery_evidence(recovery_options.run_directory / "evidence.jsonl");
  std::string recovery_first, recovery_second;
  std::getline(recovery_evidence, recovery_first);
  std::getline(recovery_evidence, recovery_second);
  const bool recovery_ok = recovery_status && recovery_report.accepted_steps == 6U &&
      recovery_first.find("\"bdf_order\":1") != std::string::npos &&
      recovery_first.find("\"restart_recovery\":true") != std::string::npos &&
      recovery_first.find("\"policy\":\"rebuild_method_history\"") != std::string::npos &&
      recovery_second.find("\"bdf_order\":2") != std::string::npos &&
      recovery_second.find("\"restart_recovery\":false") != std::string::npos;
  if (!recovery_ok) std::cerr << "FAIL: ApplicationService explicit method recovery must use BE then BDF2 and record its policy\n";
  passed &= recovery_ok;
  ApplicationRunOptions invalid_recovery = recovery_options;
  invalid_recovery.restart_directory.clear();
  invalid_recovery.run_directory = root / "recovery-without-source";
  ApplicationRunReport invalid_recovery_report;
  passed &= ApplicationService::run(MPI_COMM_SELF, invalid_recovery, invalid_recovery_report).code == StatusCode::invalid_case &&
      !fs::exists(invalid_recovery.run_directory);

  ApplicationRunOptions explicit_options = run_options;
  explicit_options.run_directory = root / "run-explicit-initial";
  explicit_options.steps = 1U;
  explicit_options.output_interval = 0U;
  explicit_options.initial_state = DriverInitialState{};
  explicit_options.initial_state->temperature = 350.0;
  ApplicationRunReport explicit_report;
  auto explicit_status = ApplicationService::run(MPI_COMM_SELF, explicit_options, explicit_report);
  ValidatedModel explicit_model;
  CompiledCasePlan explicit_plan;
  ProductDriver explicit_reader;
  RestartExpected explicit_expected;
  RestartImage explicit_image;
  if (explicit_status) explicit_status = CaseCompiler::load_and_compile(MPI_COMM_SELF, case_root, explicit_model);
  if (explicit_status) explicit_status = ProductCompiler::compile(MPI_COMM_SELF, explicit_model, case_root, explicit_plan);
  if (explicit_status) explicit_status = ProductDriver::create(MPI_COMM_SELF, std::move(explicit_plan), explicit_reader);
  if (explicit_status) explicit_status = explicit_reader.restart_expected(explicit_expected);
  if (explicit_status) explicit_status = RestartReader::load(MPI_COMM_SELF, explicit_options.run_directory / "Restart", explicit_expected, explicit_image);
  bool initial_ok = explicit_status && explicit_report.accepted_steps == 1U;
  bool have_enthalpy = false;
  for (const auto& f : explicit_image.fields) if (f.role == RestartFieldRole::enthalpy) {
    have_enthalpy = true;
    for (double h : f.values) initial_ok &= std::abs(h - 3.5*kUniversalGasConstant/28.96546*350.0) < 1e-6;
  }
  initial_ok &= have_enthalpy;
  if (!initial_ok) std::cerr << "FAIL: explicit application initial temperature must set the EOS-consistent initial enthalpy\n";
  passed &= initial_ok;
  explicit_options.restart_directory = run_options.run_directory / "Restart";
  explicit_options.run_directory = root / "ambiguous-initial-and-restart";
  passed &= ApplicationService::run(MPI_COMM_SELF, explicit_options, explicit_report).code == StatusCode::invalid_case &&
      !fs::exists(explicit_options.run_directory);

  const auto conflicting_case = root / "conflicting-initial-hints";
  passed &= static_cast<bool>(ApplicationService::initialize_case_directory(conflicting_case));
  {
    std::ifstream input(conflicting_case / "case.json");
    std::string definition{std::istreambuf_iterator<char>(input), {}};
    const auto temperature = definition.find("\"temperature\":300");
    passed &= temperature != std::string::npos;
    if (temperature != std::string::npos) definition.replace(temperature, 17U, "\"temperature\":350");
    std::ofstream output(conflicting_case / "case.json");
    output << definition;
  }
  ApplicationRunOptions conflicting_options = run_options;
  conflicting_options.case_root = conflicting_case;
  conflicting_options.run_directory = root / "run-conflicting-hints";
  ApplicationRunReport conflicting_report;
  const auto conflicting_status = ApplicationService::run(MPI_COMM_SELF, conflicting_options, conflicting_report);
  const bool conflict_rejected = conflicting_status.code == StatusCode::invalid_case &&
      conflicting_report.failure_phase == ApplicationFailurePhase::initialize &&
      conflicting_report.accepted_steps == 0U && !fs::exists(conflicting_options.run_directory);
  if (!conflict_rejected) std::cerr << "FAIL: conflicting boundary hints must not silently select the last temperature\n";
  passed &= conflict_rejected;
  conflicting_options.initial_state = DriverInitialState{};
  conflicting_options.run_directory = root / "run-explicit-over-hints";
  passed &= static_cast<bool>(ApplicationService::run(MPI_COMM_SELF, conflicting_options, conflicting_report));

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
  ApplicationRunOptions capped = benchmark_options;
  capped.run_directory = root / "cap";
  capped.time_limits.maximum_dt = 1.0e-4;
  ApplicationRunReport capped_report;
  passed &= static_cast<bool>(ApplicationService::run(MPI_COMM_SELF, capped, capped_report));
  passed &= capped_report.accepted_steps == 2U &&
            std::abs(capped_report.final_time - 2.0e-4) < 1.0e-18;
  ApplicationRunReport unchanged_run;
  unchanged_run.product = UINT64_C(0xcafef00d);
  ApplicationRunOptions invalid_run = run_options;
  invalid_run.run_directory = case_root / "inside";
  const Status invalid_status =
      ApplicationService::run(MPI_COMM_SELF, invalid_run, unchanged_run);
  // Run reports describe this invocation, including early input failures.
  // CaseValidationReport below intentionally retains its separate contract.
  passed &= !invalid_status && unchanged_run.product == 0U &&
            unchanged_run.accepted_steps == 0U &&
            unchanged_run.final_time == 0.0 &&
            unchanged_run.failed_stage == 0U && unchanged_run.attempts == 0U &&
            unchanged_run.failure_phase == ApplicationFailurePhase::input &&
            unchanged_run.failure.code == invalid_status.code &&
            unchanged_run.failure.detail == invalid_status.detail;

  for (const char* blocked : {"Visit", "Restart", "evidence.jsonl"}) {
    ApplicationRunOptions failure_options = run_options;
    failure_options.run_directory = root / (std::string("blocked-") + blocked);
    failure_options.steps = 1U;
    failure_options.output_interval = std::string(blocked) == "Visit" ? 1U : 0U;
    failure_options.restart_interval =
        std::string(blocked) == "Restart" ? 1U : 0U;
    fs::create_directories(failure_options.run_directory);
    const fs::path blocked_path = failure_options.run_directory / blocked;
    if (std::string(blocked) == "evidence.jsonl")
      fs::create_directory(blocked_path);
    else
      std::ofstream(blocked_path) << "test-owned output blocker\n";
    ApplicationRunReport failure_report;
    const Status failure =
        ApplicationService::run(MPI_COMM_SELF, failure_options, failure_report);
    const ApplicationFailurePhase expected_phase =
        std::string(blocked) == "Visit"     ? ApplicationFailurePhase::visit
        : std::string(blocked) == "Restart" ? ApplicationFailurePhase::restart
                                            : ApplicationFailurePhase::evidence;
    if (failure || failure_report.accepted_steps != 1U ||
        !(failure_report.final_time > 0.0) ||
        failure_report.failure.code != failure.code ||
        failure_report.failure.detail != failure.detail ||
        failure_report.failed_stage != 0U ||
        failure_report.failure_phase != expected_phase) {
      std::cerr << "FAIL: post-commit " << blocked
                << " failure retains accepted step/time and I/O status; "
                << "accepted_steps=" << failure_report.accepted_steps << '\n';
      passed = false;
    }
  }

  passed &= source_cases() == original_source_cases;
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
