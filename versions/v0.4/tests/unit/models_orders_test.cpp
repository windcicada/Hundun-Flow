// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_cantera.hpp"
#include "models_roundoff_detail.hpp"
#include "../support/chemistry_test_support.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    using namespace hundun::v04;
    {
      std::vector<double> valid{.4, .6}, exhausted{-1e-30, .4, .6};
      HUNDUN_CHECK(chemistry::detail::bound_chemistry_roundoff(valid, 1e-18));
      HUNDUN_CHECK(valid == std::vector<double>({.4, .6}));
      HUNDUN_CHECK(chemistry::detail::bound_chemistry_roundoff(exhausted, 1e-18));
      HUNDUN_CHECK(exhausted == std::vector<double>({0, .4, .6}));
      for (auto invalid : {std::vector<double>{-1e-17, .4, .6},
                           std::vector<double>{-5e-18, -5e-18, .4, .6},
                           std::vector<double>{-1e-10, .4, .6}}) {
        const auto before = invalid;
        HUNDUN_CHECK(!chemistry::detail::bound_chemistry_roundoff(invalid, 1e-18));
        HUNDUN_CHECK(invalid == before);
      }
      std::vector<double> excessive{-1e-10, .4, .6};
      HUNDUN_CHECK(!chemistry::detail::bound_chemistry_roundoff(excessive, 1e-6));
    }
    const bool reject = argc == 4 && std::string(argv[1]) == "--reject";
    HUNDUN_CHECK(argc == 2 || argc == 3 || reject);
    chemistry::CanteraBackendConfig config;
    config.mechanism = {reject ? argv[2] : argv[1],
        "15ebe7bb6d266e1dbcac745b8741cb631e1c29d3095f5396533326e9a1cdc85a",
        "synthetic-gas"};
    if (reject || argc == 3)
      config.mechanism.sha256 = argv[reject ? 3 : 2];
    config.species_names = {"A", "B"};
    config.chemistry = {1e-11, 1e-18, 20000};
    if (reject) {
      bool rejected = false;
      try {
        chemistry::CanteraBackendRuntime runtime(config);
      } catch (const std::exception &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
      return;
    }
    auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(config);
    chemistry::CanteraWorkspacePool pool(runtime, 2);
    auto gas = chemistry::make_cantera_backend(config, pool);
    auto other = chemistry::make_cantera_backend(config, pool);
    double y[]{0.5, 0.5}, diffusion[2]{}, h[2]{}, rates[2]{};
    portable::GasQuery q{{1, 1, 1}, gas->composition().fingerprint,
        portable::GasStateCoordinates::pressure_temperature, 100000, 0, 1000,
        y, 2};
    portable::GasQueryOutput out{{}, diffusion, h, rates, 2};
    const double mw = gas->gas_identity().molecular_weights_kg_per_kmol[0];
    for (double temperature : {500., 1500.})
      for (double pressure : {10000., 100000., 1000000.})
        for (double fraction : {0., 1e-12, .25, .75, 1.}) {
          q.temperature_k = temperature;
          q.pressure_pa = pressure;
          y[0] = fraction;
          y[1] = 1 - fraction;
          HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::success);
          const double c = out.sample.density_kg_per_m3 / mw;
          const double x = c * y[0] / .01;
          const double forward = x < 1
              ? std::pow(.01, .25) * (1.75 * x - .75 * x * x)
              : std::pow(c * y[0], .25);
          const double expected = .0002 * mw * forward /
              std::pow(std::max(.1, c * y[1]), .75);
          HUNDUN_CHECK_NEAR(rates[1], expected, 2e-13 * std::max(1., expected));
          HUNDUN_CHECK_NEAR(rates[0] + rates[1], 0., 1e-15);
        }
    // Equal thermo and molecular weights keep T and total concentration fixed.
    // In the quadratic tail with capped inhibitor, y' = -a*y + b*y*y.
    q.temperature_k = 1000;
    q.pressure_pa = 100000;
    y[0] = .1;
    y[1] = .9;
    HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::success);
    const double concentration = out.sample.density_kg_per_m3 / mw;
    HUNDUN_CHECK(concentration * y[0] < .01 && concentration < .1);
    const double a = .0002 * 1.75 * std::pow(.01, -.75) * std::pow(.1, -.75);
    const double b = .0002 * .75 * concentration * std::pow(.01, -1.75) *
                     std::pow(.1, -.75);
    const double expected = a * y[0] /
        (b * y[0] + (a - b * y[0]) * std::exp(4 * a));
    const double initial_h = out.sample.enthalpy_j_per_kg;
    q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
    q.enthalpy_j_per_kg = initial_h;
    double final[2]{}, delta[2]{};
    portable::GasAdvanceOutput advanced{{}, final, delta, 2};
    HUNDUN_CHECK(gas->advance_gas({q, 0, 4}, advanced) == portable::Status::success);
    HUNDUN_CHECK_NEAR(final[0], expected, 3e-10);
    HUNDUN_CHECK_NEAR(final[0] + final[1], 1., 1e-14);
    HUNDUN_CHECK_NEAR(advanced.final_sample.enthalpy_j_per_kg, initial_h, 1e-7);
    HUNDUN_CHECK(other->advance_gas({q, 0, 2}, advanced) == portable::Status::success);
    y[0] = final[0];
    y[1] = final[1];
    HUNDUN_CHECK(other->advance_gas({q, 2, 2}, advanced) == portable::Status::success);
    HUNDUN_CHECK_NEAR(final[0], expected, 3e-10);
  });
}
