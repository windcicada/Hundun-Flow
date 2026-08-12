// SPDX-License-Identifier: Apache-2.0

#include "src/chem_workspace_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"

#include "tests/support/test_main.hpp"

#include <filesystem>
#include <memory>
#include <stdexcept>
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

void compare_state(hundun::chemistry::CanteraBackend &backend,
                   Cantera::ThermoPhase &reference, double temperature,
                   double pressure, const std::vector<double> &fractions) {
  reference.setState_TPY(temperature, pressure, fractions.data());
  const hundun::chemistry::ThermochemicalPoint point{
      pressure, reference.enthalpy_mass(), fractions};
  const auto actual = backend.evaluate(point);
  HUNDUN_CHECK_NEAR(actual.temperature_k, reference.temperature(), 1.0e-9);
  HUNDUN_CHECK_NEAR(actual.density_kg_per_m3, reference.density(), 1.0e-12);
  HUNDUN_CHECK_NEAR(actual.cp_j_per_kg_k, reference.cp_mass(), 1.0e-9);
  HUNDUN_CHECK_NEAR(actual.mixture_molecular_weight_kg_per_kmol,
                    reference.meanMolecularWeight(), 1.0e-12);
  HUNDUN_CHECK_NEAR(point.p0_pa, reference.pressure(), 1.0e-10);
  HUNDUN_CHECK_NEAR(point.h_tc_j_per_kg, reference.enthalpy_mass(), 1.0e-8);
}

void test_inversion_and_invalid_inputs(const std::filesystem::path &mechanism) {
  const auto resolved = config(mechanism);
  auto runtime =
      std::make_shared<hundun::chemistry::CanteraBackendRuntime>(resolved);
  hundun::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto backend = hundun::chemistry::make_cantera_backend(resolved, pool);
  auto reference =
      Cantera::newSolution(mechanism.string(), "synthetic-gas",
                           "mixture-averaged")
          ->thermo();
  compare_state(*backend, *reference, 600.0, 101325.0, {0.9, 0.1});
  compare_state(*backend, *reference, 1400.0, 220000.0, {0.2, 0.8});

  for (const auto &point :
       std::vector<hundun::chemistry::ThermochemicalPoint>{
           {101325.0, 0.0, {-0.1, 1.1}},
           {101325.0, 0.0, {0.4, 0.5}},
           {101325.0, 0.0, {1.0}},
           {0.0, 0.0, {0.9, 0.1}}}) {
    bool rejected = false;
    try {
      static_cast<void>(backend->evaluate(point));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  bool inversion_failed = false;
  try {
    static_cast<void>(
        backend->evaluate({101325.0, 1.0e30, {0.9, 0.1}}));
  } catch (const std::runtime_error &) {
    inversion_failed = true;
  }
  HUNDUN_CHECK(inversion_failed);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    test_inversion_and_invalid_inputs(argv[1]);
  });
}
