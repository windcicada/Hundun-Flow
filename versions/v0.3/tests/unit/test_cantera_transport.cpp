// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/chem_workspace_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/Transport.h"

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
                   Cantera::Solution &reference, double temperature,
                   double pressure, const std::vector<double> &fractions) {
  auto thermo = reference.thermo();
  thermo->setState_TPY(temperature, pressure, fractions.data());
  const hundun::chemistry::ThermochemicalPoint point{
      pressure, thermo->enthalpy_mass(), fractions};
  const auto properties = backend.evaluate(point);
  const auto actual = backend.evaluate(point, properties);
  auto transport = reference.transport();
  std::vector<double> expected(thermo->nSpecies());
  transport->getMixDiffCoeffs(expected.data());
  HUNDUN_CHECK_NEAR(actual.viscosity_pa_s, transport->viscosity(), 1.0e-15);
  HUNDUN_CHECK_NEAR(actual.conductivity_w_per_m_k,
                    transport->thermalConductivity(), 1.0e-12);
  HUNDUN_CHECK(actual.mixture_diffusivity_m2_per_s.size() ==
               expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    HUNDUN_CHECK_NEAR(actual.mixture_diffusivity_m2_per_s[index],
                      expected[index], 1.0e-15);
    HUNDUN_CHECK(actual.mixture_diffusivity_m2_per_s[index] > 0.0);
  }
}

void test_transport(const std::filesystem::path &mechanism) {
  const auto resolved = config(mechanism);
  auto runtime =
      std::make_shared<hundun::chemistry::CanteraBackendRuntime>(resolved);
  hundun::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto backend = hundun::chemistry::make_cantera_backend(resolved, pool);
  auto reference =
      Cantera::newSolution(mechanism.string(), "synthetic-gas",
                           "mixture-averaged");
  compare_state(*backend, *reference, 600.0, 101325.0, {0.9, 0.1});
  compare_state(*backend, *reference, 1400.0, 220000.0, {0.2, 0.8});

  const hundun::chemistry::ThermochemicalPoint point{
      101325.0, reference->thermo()->enthalpy_mass(), {0.2, 0.8}};
  auto mismatch = backend->evaluate(point);
  mismatch.temperature_k += 1.0;
  bool rejected = false;
  try {
    static_cast<void>(backend->evaluate(point, mismatch));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    test_transport(argv[1]);
  });
}
