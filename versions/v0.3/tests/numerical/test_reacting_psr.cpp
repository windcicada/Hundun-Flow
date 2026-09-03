// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/chem_workspace_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"

#include "tests/support/test_main.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <vector>

namespace {

hundun::config::ResolvedReactingCaseV4
config(const std::filesystem::path &mechanism, int maximum_steps = 5000) {
  hundun::config::ResolvedReactingCaseV4 value;
  value.mechanism.file = mechanism;
  value.mechanism.sha256 =
      "c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee";
  value.mechanism.phase = "synthetic-gas";
  value.species_names = {"A", "B"};
  value.initial_mass_fractions = {0.9, 0.1};
  value.chemistry.relative_tolerance = 1.0e-10;
  value.chemistry.absolute_tolerance = 1.0e-18;
  value.chemistry.maximum_internal_steps = maximum_steps;
  return value;
}

hundun::chemistry::ThermochemicalPoint
inlet(const std::filesystem::path &mechanism) {
  auto solution =
      Cantera::newSolution(mechanism.string(), "synthetic-gas",
                           "mixture-averaged");
  solution->thermo()->setState_TPX(1000.0, Cantera::OneAtm,
                                   "A:0.9, B:0.1");
  std::vector<double> fractions(solution->thermo()->nSpecies());
  solution->thermo()->getMassFractions(fractions.data());
  return {Cantera::OneAtm, solution->thermo()->enthalpy_mass(), fractions};
}

void test_psr_residual_and_identity(const std::filesystem::path &mechanism) {
  const auto resolved = config(mechanism);
  auto runtime =
      std::make_shared<hundun::chemistry::CanteraBackendRuntime>(resolved);
  hundun::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto backend = hundun::chemistry::make_cantera_backend(resolved, pool);
  const std::uint64_t fingerprint = backend->composition().fingerprint;
  const auto feed = inlet(mechanism);
  auto state = feed;
  double residual = 1.0;
  for (std::size_t iteration = 0; iteration < 80U; ++iteration) {
    const auto chemistry =
        backend->integrate({state, 0.01 * static_cast<double>(iteration),
                            1.0e-4});
    HUNDUN_CHECK(chemistry.succeeded());
    const auto thermo = backend->evaluate(chemistry.final_state);
    const auto transport = backend->evaluate(chemistry.final_state, thermo);
    HUNDUN_CHECK(backend->composition().fingerprint == fingerprint);
    HUNDUN_CHECK(transport.mixture_diffusivity_m2_per_s.size() ==
                 backend->composition().species.size());
    residual = 0.0;
    for (std::size_t species = 0; species < state.mass_fractions.size();
         ++species) {
      const double next =
          0.5 * feed.mass_fractions[species] +
          0.5 * chemistry.final_state.mass_fractions[species];
      residual = std::max(residual,
                          std::abs(next - state.mass_fractions[species]));
      state.mass_fractions[species] = next;
    }
  }
  HUNDUN_CHECK(residual < 1.0e-12);

  const auto limited_resolved = config(mechanism, 1);
  auto limited_runtime =
      std::make_shared<hundun::chemistry::CanteraBackendRuntime>(
          limited_resolved);
  hundun::chemistry::CanteraWorkspacePool limited_pool(limited_runtime, 1U);
  auto limited =
      hundun::chemistry::make_cantera_backend(limited_resolved, limited_pool);
  const auto failure = limited->integrate({feed, 0.0, 1.0});
  HUNDUN_CHECK(failure.status ==
               hundun::chemistry::ChemistryStatus::integration_failure);
  HUNDUN_CHECK(limited->composition().fingerprint == fingerprint);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    test_psr_residual_and_identity(argv[1]);
  });
}
