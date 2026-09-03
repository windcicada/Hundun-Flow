// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "flow_reacting_boundary_detail.hpp"

#include "hundun/mesh_topology.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/chem_analytic_backend.hpp"
#include "tests/support/test_main.hpp"

#include <stdexcept>

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    auto decomposition = hundun::runtime::StructuredDecomposition::create(
        mpi, {12, 12, 12}, {false, false, false});
    hundun::mesh::MeshTopology topology(decomposition);
    const auto &patch = topology.patch(0U);
    auto backend = hundun::test::make_analytic_reacting_backend_for_tests();
    const std::vector<double> y{0.25, 0.75};

    const auto adiabatic = hundun::flow::detail::resolve_reacting_wall_boundary(
        patch, {}, 101325.0, y, *backend);
    HUNDUN_CHECK(adiabatic.patch_id == patch.stable_id());
    HUNDUN_CHECK(adiabatic.species_flux_kg_per_m2_s ==
                 std::vector<double>({0.0, 0.0}));
    HUNDUN_CHECK(adiabatic.heat_flux_w_per_m2 == 0.0);

    const auto isothermal =
        hundun::flow::detail::resolve_reacting_wall_boundary(
            patch,
            {false,
             hundun::flow::detail::ReactingWallThermalPolicy::isothermal,
             300.0},
            101325.0, y, *backend);
    HUNDUN_CHECK(isothermal.wall_h_tc_j_per_kg.has_value());
    const auto thermo = backend->evaluate(
        {101325.0, *isothermal.wall_h_tc_j_per_kg, y});
    HUNDUN_CHECK_NEAR(thermo.temperature_k, 300.0, 1.0e-10);

    const auto inlet = hundun::flow::detail::resolve_reacting_inlet(
        patch, {101325.0, y, 350.0, std::nullopt}, *backend);
    HUNDUN_CHECK_NEAR(inlet.temperature_k, 350.0, 1.0e-10);
    const auto redundant = hundun::flow::detail::resolve_reacting_inlet(
        patch,
        {101325.0, y, 350.0, inlet.thermochemical.h_tc_j_per_kg},
        *backend);
    HUNDUN_CHECK(redundant.thermochemical.h_tc_j_per_kg ==
                 inlet.thermochemical.h_tc_j_per_kg);

    const auto outlet =
        hundun::flow::detail::reacting_pressure_outlet_state(
            patch, 101325.0, 9000.0, inlet.thermochemical.h_tc_j_per_kg, y,
            *backend);
    HUNDUN_CHECK(outlet.p0_pa == 101325.0);
    HUNDUN_CHECK(outlet.p0_pa != 9000.0);

    bool rejected = false;
    try {
      static_cast<void>(hundun::flow::detail::resolve_reacting_wall_boundary(
          patch, {true, {}, std::nullopt}, 101325.0, y, *backend));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    rejected = false;
    try {
      static_cast<void>(hundun::flow::detail::resolve_reacting_inlet(
          patch,
          {101325.0, y, 360.0, inlet.thermochemical.h_tc_j_per_kg},
          *backend));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);

    for (std::size_t cell = 0; cell < 12U * 12U * 12U; ++cell) {
      const auto smoke =
          hundun::flow::detail::resolve_reacting_wall_boundary(
              patch, {}, 101325.0, y, *backend);
      HUNDUN_CHECK(smoke.heat_flux_w_per_m2 == 0.0);
    }
  });
}
