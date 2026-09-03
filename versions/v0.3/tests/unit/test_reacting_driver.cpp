// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_reacting_flow_driver_detail.hpp"

#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <stdexcept>

namespace {

hundun::config::ResolvedReactingCaseV4 base_case() {
  hundun::config::ResolvedReactingCaseV4 value;
  value.common_flow.schema_version = 2;
  value.mechanism.file = "/tmp/mechanism.yaml";
  value.mechanism.sha256 = std::string(64U, 'a');
  value.mechanism.phase = "gas";
  value.initial_p0_pa = 101325.0;
  value.initial_temperature_k = 300.0;
  value.species_names = {"A", "B"};
  value.initial_mass_fractions = {0.25, 0.75};
  value.composition_fingerprint = 42U;
  for (auto &boundary : value.boundary_reacting)
    boundary.non_catalytic_impermeable = true;
  return value;
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    for (const auto pressure :
         {hundun::config::PressureConstraintMode::open_fixed_p0,
          hundun::config::PressureConstraintMode::closed}) {
      for (bool immersed : {false, true}) {
        for (bool wale : {false, true}) {
          auto resolved = base_case();
          resolved.pressure_mode = pressure;
          if (immersed) {
            resolved.immersed_boundary.model = hundun::config::
                ImmersedBoundaryModel::local_flow_pattern_ghost_cell;
            resolved.immersed_boundary.geometry =
                hundun::config::StlGeometryConfig{
                    "/tmp/body.stl", 1.0,
                    hundun::config::ImmersedFluidSide::outside};
          }
          if (wale) {
            resolved.les.model = hundun::config::LesModel::wale;
            resolved.les.wale = hundun::config::WaleConfig{0.5, 0.9, 0.7};
          }
          const auto plan =
              hundun::application::plan_reacting_flow_case(resolved);
          HUNDUN_CHECK(plan.backend_runtime_count == 1U);
          HUNDUN_CHECK(plan.workspace_pool_count == 1U);
          HUNDUN_CHECK(plan.immersed_boundary == immersed);
          HUNDUN_CHECK(plan.wale == wale);
          HUNDUN_CHECK(plan.operator_order[0] == "chemistry-half-1");
          HUNDUN_CHECK(plan.operator_order[4] == "piso-2");
          HUNDUN_CHECK(hundun::application::run_reacting_flow_case(
                           resolved, MPI_COMM_WORLD) == EXIT_SUCCESS);
        }
      }
    }

    for (int mutation = 0; mutation < 4; ++mutation) {
      auto invalid = base_case();
      if (mutation == 0)
        invalid.mechanism.file.clear();
      if (mutation == 1)
        invalid.species_names.push_back("C");
      if (mutation == 2) {
        invalid.boundary_reacting[0].non_catalytic_impermeable = false;
        invalid.boundary_reacting[0].thermal =
            hundun::config::ReactingThermalBoundaryConfig{};
      }
      if (mutation == 3)
        invalid.mechanism.sha256[0] = 'G';
      bool rejected = false;
      try {
        static_cast<void>(
            hundun::application::plan_reacting_flow_case(invalid));
      } catch (const std::invalid_argument &) {
        rejected = true;
      }
      HUNDUN_CHECK(rejected);
    }
  });
}
