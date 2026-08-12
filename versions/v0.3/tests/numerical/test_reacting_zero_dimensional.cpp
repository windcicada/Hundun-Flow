// SPDX-License-Identifier: Apache-2.0

#include "src/chem_workspace_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ReactorNet.h"

#include "tests/support/test_main.hpp"

#include <filesystem>
#include <memory>
#include <vector>

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

void test_zero_dimensional_oracle(const std::filesystem::path &mechanism) {
  const auto resolved = config(mechanism);
  auto runtime =
      std::make_shared<hundun::chemistry::CanteraBackendRuntime>(resolved);
  hundun::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto backend = hundun::chemistry::make_cantera_backend(resolved, pool);

  auto solution =
      Cantera::newSolution(mechanism.string(), "synthetic-gas",
                           "mixture-averaged");
  solution->thermo()->setState_TPX(1000.0, Cantera::OneAtm,
                                   "A:0.9, B:0.1");
  std::vector<double> initial_y(solution->thermo()->nSpecies());
  solution->thermo()->getMassFractions(initial_y.data());
  const double initial_density = solution->thermo()->density();
  const double initial_h = solution->thermo()->enthalpy_mass();
  const hundun::chemistry::ThermochemicalPoint initial{
      Cantera::OneAtm, initial_h, initial_y};

  auto reactor =
      std::make_shared<Cantera::IdealGasConstPressureReactor>(solution, false);
  reactor->syncState();
  Cantera::ReactorNet network(reactor);
  network.setTolerances(1.0e-10, 1.0e-18);
  network.setMaxSteps(5000);
  network.advance(1.0e-3);
  auto final_thermo = reactor->phase()->thermo();
  std::vector<double> expected_y(final_thermo->nSpecies());
  final_thermo->getMassFractions(expected_y.data());

  const auto actual = backend->integrate({initial, 0.0, 1.0e-3});
  HUNDUN_CHECK(actual.succeeded());
  for (std::size_t index = 0; index < expected_y.size(); ++index) {
    HUNDUN_CHECK_NEAR(actual.final_state.mass_fractions[index],
                      expected_y[index], 1.0e-11);
    HUNDUN_CHECK_NEAR(actual.integrated_rho_y_delta_kg_per_m3[index],
                      initial_density * (expected_y[index] - initial_y[index]),
                      1.0e-12);
  }
  const auto final_properties = backend->evaluate(actual.final_state);
  HUNDUN_CHECK_NEAR(final_properties.temperature_k,
                    final_thermo->temperature(), 1.0e-6);
  HUNDUN_CHECK(actual.final_state.h_tc_j_per_kg == initial_h);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    test_zero_dimensional_oracle(argv[1]);
  });
}
