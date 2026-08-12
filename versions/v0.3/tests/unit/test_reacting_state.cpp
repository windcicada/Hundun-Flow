// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow_reacting_state.hpp"

#include "tests/support/chem_analytic_backend.hpp"
#include "tests/support/test_main.hpp"

#include <stdexcept>

int main() {
  return hundun::test::run([] {
    auto backend = hundun::test::make_analytic_reacting_backend_for_tests();
    hundun::runtime::FieldRegistry registry;
    const auto fields = hundun::flow::declare_reacting_fields(
        registry, backend->composition());
    HUNDUN_CHECK(fields.species_density.size() == 2U);
    HUNDUN_CHECK(fields.total_thermochemical_enthalpy !=
                 fields.species_density[0]);
    registry.freeze();

    auto state = hundun::flow::ReactingFlowState::create(
        registry, {2, 1, 1}, fields, backend->composition());
    state.initialize_uniform(*backend, 101325.0, 300.19, {0.5, 0.5});
    const auto committed =
        state.snapshot(hundun::flow::ReactingLayer::committed);
    HUNDUN_CHECK(committed.rho_y_kg_per_m3.size() == 2U);
    HUNDUN_CHECK(committed.rho_y_kg_per_m3[0].size() == 2U);
    HUNDUN_CHECK_NEAR(committed.rho_y_kg_per_m3[0][0], 1.0, 1.0e-12);
    HUNDUN_CHECK_NEAR(committed.rho_y_kg_per_m3[1][0], 1.0, 1.0e-12);
    HUNDUN_CHECK(committed.rho_h_tc_j_per_m3[0] != 300.19);

    const auto properties = state.thermodynamics(
        hundun::flow::ReactingLayer::committed, 0U, *backend);
    const auto cache_epoch = state.cache_epoch();
    HUNDUN_CHECK_NEAR(properties.temperature_k, 300.19, 1.0e-10);
    HUNDUN_CHECK(cache_epoch == state.writer_epoch());

    state.begin_attempt();
    state.write_trial_cell(0U, {0.8, 1.2},
                           committed.rho_h_tc_j_per_m3[0]);
    HUNDUN_CHECK(state.writer_epoch() != cache_epoch);
    HUNDUN_CHECK(state.cache_epoch() != state.writer_epoch());
    state.rollback_attempt();
    HUNDUN_CHECK(state.snapshot(hundun::flow::ReactingLayer::committed)
                     .rho_y_kg_per_m3 == committed.rho_y_kg_per_m3);

    state.begin_attempt();
    state.write_trial_cell(0U, {0.8, 1.2},
                           committed.rho_h_tc_j_per_m3[0]);
    state.commit_attempt();
    HUNDUN_CHECK(state.snapshot(hundun::flow::ReactingLayer::history)
                     .rho_y_kg_per_m3 == committed.rho_y_kg_per_m3);
    HUNDUN_CHECK_NEAR(
        state.snapshot(hundun::flow::ReactingLayer::committed)
            .rho_y_kg_per_m3[0][0],
        0.8, 0.0);

    bool rejected = false;
    try {
      state.begin_attempt();
      state.write_trial_cell(1U, {-1.0, 1.0}, 1.0);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    state.rollback_attempt();
  });
}
