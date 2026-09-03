// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/flow_reacting_coupling_detail.hpp"

#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <string>

namespace {

using State = std::array<double, 2>;

State chemistry(State value, double dt) {
  value[0] += dt * value[1];
  return value;
}

State transport(State value, double dt) {
  value[1] += dt * value[0];
  return value;
}

State strang(State value, double dt) {
  return chemistry(transport(chemistry(value, 0.5 * dt), dt), 0.5 * dt);
}

State exact(double time) {
  return {std::cosh(time) + 2.0 * std::sinh(time),
          std::sinh(time) + 2.0 * std::cosh(time)};
}

double error(State value, State expected) {
  return std::hypot(value[0] - expected[0], value[1] - expected[1]);
}

double slope(double coarse, double fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

void open_smoke(double chemistry_delta, double transport_delta) {
  constexpr std::size_t cells = 8U;
  hundun::flow::detail::ReactingLayerState initial;
  initial.p0_pa = 101325.0;
  initial.rho_y_kg_per_m3.resize(cells * 2U);
  initial.rho_h_tc_j_per_m3.assign(cells, 10.0);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    initial.rho_y_kg_per_m3[cell * 2U] = 0.7;
    initial.rho_y_kg_per_m3[cell * 2U + 1U] = 0.3;
  }
  auto history = initial;
  hundun::flow::detail::ReactingAttemptState state(
      cells, 2U, std::move(history), initial);
  hundun::flow::OpenReactingStepOperators operators;
  operators.chemistry =
      [&](std::uint32_t, double, double,
          const hundun::flow::detail::ReactingLayerState &) {
        hundun::flow::ReactingOperatorUpdate update;
        update.species_delta_kg_per_m3.resize(cells * 2U);
        for (std::size_t cell = 0; cell < cells; ++cell) {
          update.species_delta_kg_per_m3[cell * 2U] = -chemistry_delta;
          update.species_delta_kg_per_m3[cell * 2U + 1U] = chemistry_delta;
        }
        return update;
      };
  operators.scalar_transport =
      [&](double, double,
          const hundun::flow::detail::ReactingLayerState &) {
        hundun::flow::ReactingOperatorUpdate update;
        update.species_delta_kg_per_m3.resize(cells * 2U);
        update.enthalpy_delta_j_per_m3.assign(cells, transport_delta);
        for (std::size_t cell = 0; cell < cells; ++cell) {
          update.species_delta_kg_per_m3[cell * 2U] = -transport_delta;
          update.species_delta_kg_per_m3[cell * 2U + 1U] = transport_delta;
        }
        return update;
      };
  operators.pressure_corrector =
      [](std::uint32_t, double, std::uint64_t,
         const hundun::flow::detail::ReactingLayerState &,
         const std::vector<double> &, const std::vector<double> &,
         std::string &) { return true; };
  operators.collective_validation =
      [](bool ok, const std::string &message) {
        return hundun::runtime::CollectiveStatus{ok, ok ? -1 : 0, message};
      };
  const auto report = hundun::flow::attempt_open_reacting_step(
      state, 0.0, 0.01, operators);
  HUNDUN_CHECK(report.accepted);
  for (std::size_t index = 0; index < initial.rho_y_kg_per_m3.size(); ++index)
    HUNDUN_CHECK_NEAR(state.committed().rho_y_kg_per_m3[index] -
                          initial.rho_y_kg_per_m3[index],
                      report.integrated_species_delta_kg_per_m3[index],
                      1.0e-15);
  for (std::size_t cell = 0; cell < cells; ++cell)
    HUNDUN_CHECK_NEAR(state.committed().rho_h_tc_j_per_m3[cell] - 10.0,
                      report.integrated_enthalpy_delta_j_per_m3[cell],
                      1.0e-15);
}

} // namespace

int main() {
  return hundun::test::run([] {
    const State initial{1.0, 2.0};
    const double local_coarse = error(strang(initial, 0.1), exact(0.1));
    const double local_fine = error(strang(initial, 0.05), exact(0.05));
    HUNDUN_CHECK(slope(local_coarse, local_fine) > 2.9);

    auto global_error = [&](std::size_t steps) {
      State value = initial;
      const double dt = 1.0 / static_cast<double>(steps);
      for (std::size_t step = 0; step < steps; ++step) {
        value = strang(value, dt);
      }
      return error(value, exact(1.0));
    };
    HUNDUN_CHECK(slope(global_error(20U), global_error(40U)) > 1.95);

    const State ct = transport(chemistry(initial, 0.1), 0.1);
    const State tc = chemistry(transport(initial, 0.1), 0.1);
    HUNDUN_CHECK(error(ct, exact(0.1)) > local_coarse);
    HUNDUN_CHECK(error(tc, exact(0.1)) > local_coarse);
    const State wrong_half =
        chemistry(transport(chemistry(initial, 0.1), 0.1), 0.1);
    HUNDUN_CHECK(error(wrong_half, exact(0.1)) > local_coarse);

    open_smoke(0.01, 0.0);
    open_smoke(0.0, 0.01);
    open_smoke(0.01, 0.01);
  });
}
