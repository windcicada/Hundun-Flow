// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09

#include "hundun/v04_cantera.hpp"
#include "../support/chemistry_test_support.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <vector>

namespace {
using namespace hundun::v04;
constexpr std::size_t species_count = 7;
using Composition = std::array<double, species_count>;
struct State {
  unsigned id{};
  double pressure{}, enthalpy{};
  Composition composition{};
};

Composition advance(chemistry::CanteraBackend& gas, const State& state) {
  Composition result{}, delta{};
  const portable::GasAdvanceQuery query{
      {{1, 32002, 1}, gas.gas_identity().composition_fingerprint,
       portable::GasStateCoordinates::pressure_enthalpy,
       state.pressure, state.enthalpy, 0, state.composition.data(), species_count},
      0.021806953522829548, 6.5569904281136158e-7};
  portable::GasAdvanceOutput output{{}, result.data(), delta.data(), species_count};
  HUNDUN_CHECK(gas.advance_gas(query, output) == portable::Status::success);
  double total{};
  for (const auto value : result) {
    HUNDUN_CHECK(std::isfinite(value) && value >= 0 && value <= 1);
    total += value;
  }
  // Gas publication and ESF state admission share this closure contract.
  HUNDUN_CHECK_NEAR(total, 1., 2e-12);
  return result;
}

void check_pairs(const char* mechanism, const char* points, double tolerance) {
  chemistry::CanteraBackendConfig config;
  config.mechanism = {mechanism,
      "a1ca61b0b97847c59ea59e2c0373451c0d1c4cb6a46d3ddde21f7a1d64c05fa1", "jl4"};
  config.species_names = {"CH4", "O2", "CO2", "CO", "H2O", "H2", "N2"};
  config.chemistry = {tolerance, 1e-16, 100000};
  config.continuous_enthalpy = true;
  config.minimum_temperature = 273.15;
  config.maximum_temperature = 3500;
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(config);
  chemistry::CanteraWorkspacePool pool(runtime, 1);
  auto gas = chemistry::make_cantera_backend(config, pool);
  auto fine_config = config;
  fine_config.chemistry.relative_tolerance = 1e-12;
  auto fine_runtime = std::make_shared<chemistry::CanteraBackendRuntime>(fine_config);
  chemistry::CanteraWorkspacePool fine_pool(fine_runtime, 1);
  auto fine = chemistry::make_cantera_backend(fine_config, fine_pool);

  std::ifstream input(points);
  HUNDUN_CHECK(input.good());
  std::vector<State> states;
  State state;
  while (input >> state.id) {
    HUNDUN_CHECK(state.id == states.size());
    HUNDUN_CHECK(input >> state.pressure >> state.enthalpy);
    for (auto& value : state.composition) HUNDUN_CHECK(input >> value);
    states.push_back(state);
  }
  HUNDUN_CHECK(input.eof() && states.size() == 120);
  std::vector<Composition> answers;
  double pair_error{}, refinement_error{};
  for (const auto& value : states) {
    answers.push_back(advance(*gas, value));
    const auto reference = advance(*fine, value);
    for (std::size_t c = 0; c < species_count; ++c)
      refinement_error = std::max(refinement_error,
          std::abs(answers.back()[c] - reference[c]));
  }
  for (std::size_t i = 0; i < states.size(); i += 2) {
    HUNDUN_CHECK(states[i].pressure == states[i + 1].pressure);
    HUNDUN_CHECK(states[i].enthalpy == states[i + 1].enthalpy);
    for (std::size_t c = 0; c < species_count; ++c) {
      HUNDUN_CHECK_NEAR(states[i].composition[c], states[i + 1].composition[c], 1e-14);
      pair_error = std::max(pair_error, std::abs(answers[i][c] - answers[i + 1][c]));
    }
  }
  // Repeat actual queries after other cells and with independent workspaces.
  for (std::size_t i = states.size(); i-- > 0;) {
    HUNDUN_CHECK(advance(*gas, states[i]) == answers[i]);
    chemistry::CanteraWorkspacePool fresh_pool(runtime, 1);
    auto fresh = chemistry::make_cantera_backend(config, fresh_pool);
    HUNDUN_CHECK(advance(*fresh, states[i]) == answers[i]);
  }
  std::cout << std::setprecision(17) << "chemistry_pairs count=60 rtol=" << tolerance
            << " pair_max=" << pair_error << " refinement_max=" << refinement_error
            << " limit=1e-9 order=exact workspace=exact\n";
  HUNDUN_CHECK(pair_error < 1e-9);
  HUNDUN_CHECK(refinement_error < 1e-9);
}
} // namespace

int main(int argc, char** argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 3 || argc == 4);
    // The optional tolerance reproduces the original 1e-8 failing configuration.
    check_pairs(argv[1], argv[2], argc == 4 ? std::stod(argv[3]) : 1e-10);
  });
}
