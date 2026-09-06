// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_app.hpp"
#include "core_product_freeze_detail.hpp"
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <utility>

// Manual, read-only-checkpoint method probe; never registered as a passing
// CFD regression and never writes a new checkpoint or a production run.
int main(int argc, char** argv) {
  using namespace hundun::v04;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int result = 2;
  {
    ValidatedModel model;
    Status status = argc == 3 ? CaseCompiler::load_and_compile(MPI_COMM_WORLD, argv[1], model)
                             : Status{StatusCode::invalid_case, 1U};
    CompiledCasePlan plan;
    if (status) status = ProductCompiler::compile(MPI_COMM_WORLD, model, argv[1], plan);
    ProductDriver driver;
    if (status) status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
    RestartExpected expected;
    if (status) status = driver.restart_expected(expected);
    RestartImage image;
    if (status) status = RestartReader::load(MPI_COMM_WORLD, argv[2], expected, image);
    image.backward_euler_recovery = true;
    if (status) status = driver.initialize_restart(image);
    DriverStepReport report;
    if (status) {
      detail::arm_pressure_energy_candidate_globalization_once_for_test();
      status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, report);
    }
    detail::PressureEnergyCandidateGlobalizationDiagnostic diagnostic;
    const bool observed = detail::pressure_energy_candidate_globalization_diagnostic_for_test(diagnostic);
    if (rank == 0) {
      std::cout << std::setprecision(17) << "method-probe status="
                << static_cast<unsigned>(status.code) << '/' << status.detail
                << " accepted=" << report.accepted << " stage=" << report.failed_stage
                << " diagnostics=" << observed << '\n';
      if (observed) {
        std::cout << "linear-prediction C1="
                  << diagnostic.corrector_one_linear_predicted_normalized_continuity << ','
                  << diagnostic.corrector_one_linear_predicted_normalized_energy << " C2="
                  << diagnostic.corrector_two_linear_predicted_normalized_continuity << ','
                  << diagnostic.corrector_two_linear_predicted_normalized_energy << '\n';
        for (std::size_t slot = 0U; slot < diagnostic.linear_target_gap.size(); ++slot)
          std::cout << "linear-target-gap slot=" << slot << " C="
                    << diagnostic.linear_target_gap[slot][0U] << " E="
                    << diagnostic.linear_target_gap[slot][1U] << " observed="
                    << diagnostic.linear_target_observed[slot] << " worst-cell="
                    << diagnostic.linear_target_continuity_worst_cell[slot] << '\n';
      }
    }
    result = observed ? 0 : 2;
    detail::clear_pressure_energy_candidate_globalization_for_test();
  }
  MPI_Finalize();
  return result;
}
