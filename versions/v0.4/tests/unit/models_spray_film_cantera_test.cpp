// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "hundun/v04_cantera.hpp"
#include "models_spray_properties_detail.hpp"
#include "../support/chemistry_test_support.hpp"
#include <array>
#include <iomanip>
#include <string_view>

using namespace hundun::v04;
using namespace spray::detail;
int main(int argc, char** argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2 || argc == 3);
    chemistry::CanteraBackendConfig config;
    config.mechanism = {argv[1],
        "8b392d01bfc64baec3484c766cd02571ea3084f57e483c1635c8d0bfa5a6cad6", "film-gas"};
    config.species_names = {"O2", "N2", "C12H23"};
    config.chemistry = {1e-10, 1e-16, 10000};
    config.continuous_enthalpy = true;
    auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(config);
    chemistry::CanteraWorkspacePool pool(runtime, 1);
    auto gas = chemistry::make_cantera_backend(config, pool);
    LiquidAsset asset;
    asset.gas_identity = gas->gas_identity();
    asset.vapor_species_index = 2;
    asset.vapor_species_name = "C12H23";
    asset.vapor_molecular_weight_kg_per_kmol = 167.;
    asset.pack.minimum_temperature_k = 200.;
    asset.pack.maximum_temperature_k = 680.;
    const double y[]{.23, .76, .01};
    const double mw[]{32., 28., 167.}, a[]{3.5, 3.5, 30.}, b[]{.0001, .0002, .002};
    constexpr double r = 8314.46261815324;
    double far_h{};
    for (unsigned i=0; i<3; ++i) {
      HUNDUN_CHECK(asset.gas_identity.molecular_weights_kg_per_kmol[i] == mw[i]);
      far_h += y[i]*(r/mw[i]*(a[i]*900.+.5*b[i]*900.*900.)-(i==2 ? 1e6 : 0.));
    }
    KeroseneFilmInput input;
    input.expected_revision = input.transport_revision = {32004, 3, 1};
    input.far_gas = {input.expected_revision, asset.gas_identity.composition_fingerprint,
        portable::GasStateCoordinates::pressure_enthalpy, 101325., far_h, 0., y, 3};
    input.far_dynamic_viscosity_pa_s = 2e-5;
    KeroseneFilmWorkspace workspace(3);
    const auto sample = [&](double t, double p) {
      input.surface_temperature_k = t;
      input.far_gas.pressure_pa = p;
      const auto report = workspace.query(asset, *gas, input);
      HUNDUN_CHECK(report.available && report.revision == input.expected_revision);
      HUNDUN_CHECK_NEAR(report.far_temperature_k, 900., 1e-9);
      double cp{};
      for (unsigned i=0; i<3; ++i) cp += y[i]*r/mw[i]*(a[i]+b[i]*report.film_temperature_k);
      HUNDUN_CHECK_NEAR(report.gas_cp_j_per_kg_k, cp, 1e-9);
      HUNDUN_CHECK_NEAR(report.vapor_absolute_enthalpy_j_per_kg,
          r/167.*(30.*t+.001*t*t)-1e6, 1e-8);
      return report;
    };
    if (argc == 3 && std::string_view(argv[2]) == "--sample") {
      double t,p,unused;
      std::cout << std::setprecision(17);
      while (std::cin >> t >> p >> unused) {
        const auto s = sample(t,p);
        std::cout << s.film_temperature_k << ' ' << s.gas_cp_j_per_kg_k << ' '
            << s.gas_dynamic_viscosity_pa_s << ' ' << s.vapor_cp_j_per_kg_k << ' '
            << s.vapor_prandtl_number << ' ' << s.vapor_absolute_enthalpy_j_per_kg << ' '
            << s.surface_vapor_mass_fraction << '\n';
      }
      return;
    }
    std::array<KeroseneFilmReport, 6> forward{};
    const double temperatures[]{200., 350., 400., 477.95, 550., 650.};
    for (unsigned i=0; i<6; ++i) forward[i] = sample(temperatures[i],101325.);
    for (unsigned i=6; i-- > 0;) {
      const auto repeated = sample(temperatures[i],101325.);
      HUNDUN_CHECK(repeated.gas_cp_j_per_kg_k == forward[i].gas_cp_j_per_kg_k);
      HUNDUN_CHECK(repeated.vapor_absolute_enthalpy_j_per_kg == forward[i].vapor_absolute_enthalpy_j_per_kg);
    }
    // An intervening failed state keeps subsequent PH and pure-vapor queries valid.
    input.surface_temperature_k = 100.;
    HUNDUN_CHECK(!workspace.query(asset,*gas,input).available);
    HUNDUN_CHECK(sample(350.,101325.).gas_cp_j_per_kg_k == forward[1].gas_cp_j_per_kg_k);
    std::cout << "thick_film_real_backend PH_and_pure_vapor=passed order=exact recovery=passed\n";
  });
}
