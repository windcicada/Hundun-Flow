// SPDX-License-Identifier: Apache-2.0

#include "src/chem_workspace_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"

#include "tests/support/test_main.hpp"

#include <cmath>
#include <filesystem>
#include <memory>
#include <numeric>

namespace {

hundun::config::ResolvedReactingCaseV4
config(const std::filesystem::path &mechanism) {
  hundun::config::ResolvedReactingCaseV4 value;
  value.mechanism.file = mechanism;
  value.mechanism.sha256 =
      "c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee";
  value.mechanism.phase = "synthetic-gas";
  value.species_names = {"A", "B"};
  value.initial_mass_fractions = {0.9, 0.1};
  value.chemistry.relative_tolerance = 1.0e-10;
  value.chemistry.absolute_tolerance = 1.0e-18;
  value.chemistry.maximum_internal_steps = 5000;
  return value;
}

hundun::chemistry::ThermochemicalPoint
initial_point(const std::filesystem::path &mechanism) {
  auto solution =
      Cantera::newSolution(mechanism.string(), "synthetic-gas",
                           "mixture-averaged");
  solution->thermo()->setState_TPX(1000.0, Cantera::OneAtm,
                                   "A:0.9, B:0.1");
  std::vector<double> fractions(solution->thermo()->nSpecies());
  solution->thermo()->getMassFractions(fractions.data());
  return {Cantera::OneAtm, solution->thermo()->enthalpy_mass(), fractions};
}

std::unique_ptr<hundun::chemistry::CanteraBackend>
backend(const hundun::config::ResolvedReactingCaseV4 &resolved,
        std::shared_ptr<hundun::chemistry::CanteraBackendRuntime> &runtime,
        std::unique_ptr<hundun::chemistry::CanteraWorkspacePool> &pool) {
  runtime =
      std::make_shared<hundun::chemistry::CanteraBackendRuntime>(resolved);
  pool = std::make_unique<hundun::chemistry::CanteraWorkspacePool>(runtime, 1U);
  return hundun::chemistry::make_cantera_backend(resolved, *pool);
}

void test_interval(const std::filesystem::path &mechanism) {
  const auto resolved = config(mechanism);
  std::shared_ptr<hundun::chemistry::CanteraBackendRuntime> runtime;
  std::unique_ptr<hundun::chemistry::CanteraWorkspacePool> pool;
  auto chemistry = backend(resolved, runtime, pool);
  const auto initial = initial_point(mechanism);
  const auto full = chemistry->integrate({initial, 2.0, 1.0e-3});
  HUNDUN_CHECK(full.succeeded());
  HUNDUN_CHECK(full.completed_duration_s == 1.0e-3);
  HUNDUN_CHECK(full.internal_step_count > 0U);
  HUNDUN_CHECK(full.final_state.p0_pa == initial.p0_pa);
  HUNDUN_CHECK(full.final_state.h_tc_j_per_kg == initial.h_tc_j_per_kg);
  HUNDUN_CHECK(full.final_state.mass_fractions != initial.mass_fractions);
  HUNDUN_CHECK_NEAR(
      std::accumulate(full.integrated_rho_y_delta_kg_per_m3.begin(),
                      full.integrated_rho_y_delta_kg_per_m3.end(), 0.0),
      0.0, 1.0e-12);

  auto first = chemistry->integrate({initial, 2.0, 5.0e-4});
  auto second = chemistry->integrate({first.final_state, 2.0005, 5.0e-4});
  HUNDUN_CHECK(first.succeeded());
  HUNDUN_CHECK(second.succeeded());
  for (std::size_t index = 0; index < full.final_state.mass_fractions.size();
       ++index) {
    HUNDUN_CHECK_NEAR(second.final_state.mass_fractions[index],
                      full.final_state.mass_fractions[index], 1.0e-9);
  }

  auto invalid = initial;
  invalid.mass_fractions = {-0.1, 1.1};
  const auto failed = chemistry->integrate({invalid, 0.0, 1.0e-3});
  HUNDUN_CHECK(!failed.succeeded());
  HUNDUN_CHECK(failed.final_state.mass_fractions == invalid.mass_fractions);
  HUNDUN_CHECK(failed.integrated_rho_y_delta_kg_per_m3 ==
               std::vector<double>({0.0, 0.0}));
  HUNDUN_CHECK(failed.completed_duration_s == 0.0);
  HUNDUN_CHECK(failed.internal_step_count == 0U);
  for (double delta : failed.integrated_rho_y_delta_kg_per_m3) {
    HUNDUN_CHECK(delta == 0.0);
    HUNDUN_CHECK(!std::signbit(delta));
  }

  const auto identity = chemistry->integrate({initial, 3.0, 0.0});
  HUNDUN_CHECK(identity.succeeded());
  HUNDUN_CHECK(identity.final_state.mass_fractions ==
               initial.mass_fractions);

  auto impossible = initial;
  impossible.h_tc_j_per_kg = 1.0e30;
  const auto inversion = chemistry->integrate({impossible, 0.0, 1.0e-3});
  HUNDUN_CHECK(inversion.status ==
               hundun::chemistry::ChemistryStatus::state_inversion_failure);

  auto limited_config = resolved;
  limited_config.chemistry.maximum_internal_steps = 1;
  std::shared_ptr<hundun::chemistry::CanteraBackendRuntime> limited_runtime;
  std::unique_ptr<hundun::chemistry::CanteraWorkspacePool> limited_pool;
  auto limited = backend(limited_config, limited_runtime, limited_pool);
  const auto integration_failure =
      limited->integrate({initial, 0.0, 1.0});
  HUNDUN_CHECK(integration_failure.status ==
               hundun::chemistry::ChemistryStatus::integration_failure);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    test_interval(argv[1]);
  });
}
