// SPDX-License-Identifier: Apache-2.0

#include "src/flow_reacting_transaction_detail.hpp"
#include "tests/support/test_main.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace {

using hundun::flow::detail::ReactingAttemptState;
using hundun::flow::detail::ReactingLayerState;
using hundun::flow::detail::ReactingSourceIdentity;
using hundun::flow::detail::ReactingSourceKind;
using hundun::flow::detail::ReactingSourceTransaction;

ReactingLayerState layer(double offset) {
  return {{1.0 + offset, 2.0 + offset, 3.0 + offset, 4.0 + offset},
          {10.0 + offset, 20.0 + offset}, 101325.0 + offset};
}

template <class Function> void expect_rejected(Function &&function) {
  bool rejected = false;
  try {
    function();
  } catch (const std::exception &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_all_species_accumulation_provenance_and_publish_once() {
  ReactingAttemptState state(2U, 2U, layer(-1.0), layer(0.0));
  const auto old_committed = state.committed();
  state.begin_attempt();
  ReactingSourceTransaction transaction(state);
  const ReactingSourceIdentity chemistry{"chem-half-1",
                                          ReactingSourceKind::chemistry};
  const ReactingSourceIdentity transport{"scalar-transport",
                                          ReactingSourceKind::transport};

  transaction.add_species(chemistry, 0U, 0U, -0.25);
  transaction.add_species(chemistry, 0U, 1U, 0.25);
  transaction.add_species(chemistry, 1U, 0U, -0.5);
  transaction.add_species(chemistry, 1U, 1U, 0.5);
  transaction.add_enthalpy(transport, 0U, 5.0);
  transaction.add_enthalpy(transport, 1U, -5.0);

  HUNDUN_CHECK(transaction.species_delta_kg_per_m3() ==
               std::vector<double>({-0.25, 0.25, -0.5, 0.5}));
  HUNDUN_CHECK(transaction.enthalpy_delta_j_per_m3() ==
               std::vector<double>({5.0, -5.0}));
  HUNDUN_CHECK(transaction.records().size() == 6U);
  HUNDUN_CHECK(transaction.records()[0].source.id == "chem-half-1");
  HUNDUN_CHECK(transaction.records()[0].units == "kg/m^3");
  HUNDUN_CHECK(transaction.records()[4].units == "J/m^3");

  const hundun::runtime::CollectiveStatus accepted{true, -1, {}};
  HUNDUN_CHECK(transaction.commit(accepted));
  HUNDUN_CHECK(state.history() == old_committed);
  HUNDUN_CHECK(state.committed().rho_y_kg_per_m3 ==
               std::vector<double>({0.75, 2.25, 2.5, 4.5}));
  HUNDUN_CHECK(state.committed().rho_h_tc_j_per_m3 ==
               std::vector<double>({15.0, 15.0}));
  expect_rejected([&] { static_cast<void>(transaction.commit(accepted)); });
}

void test_invalid_sources_and_local_rollback() {
  ReactingAttemptState state(2U, 2U, layer(-1.0), layer(0.0));
  const auto original_history = state.history();
  const auto original_committed = state.committed();
  state.begin_attempt();
  ReactingSourceTransaction transaction(state);
  const ReactingSourceIdentity chemistry{"chemistry",
                                          ReactingSourceKind::chemistry};
  expect_rejected(
      [&] { transaction.add_enthalpy(chemistry, 0U, 1.0); });
  transaction.add_species(chemistry, 0U, 0U, -0.1);
  transaction.add_species(chemistry, 0U, 1U, 0.1);
  transaction.rollback();
  HUNDUN_CHECK(state.history() == original_history);
  HUNDUN_CHECK(state.committed() == original_committed);
  HUNDUN_CHECK(!state.attempt_active());

  state.begin_attempt();
  ReactingSourceTransaction invalid_mass(state);
  invalid_mass.add_species(chemistry, 0U, 0U, -0.1);
  const hundun::runtime::CollectiveStatus accepted{true, -1, {}};
  expect_rejected(
      [&] { static_cast<void>(invalid_mass.commit(accepted)); });
  HUNDUN_CHECK(state.committed() == original_committed);
  HUNDUN_CHECK(!state.attempt_active());
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_all_species_accumulation_provenance_and_publish_once();
    test_invalid_sources_and_local_rollback();
  });
}
