// SPDX-License-Identifier: Apache-2.0

#include "src/chem_workspace_detail.hpp"

#include "tests/support/test_main.hpp"

#include <filesystem>
#include <memory>
#include <stdexcept>

namespace {

std::filesystem::path mechanism_path;

hundun::config::ResolvedReactingCaseV4 config() {
  hundun::config::ResolvedReactingCaseV4 value;
  value.mechanism.file = mechanism_path;
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

void test_exact_identity_and_independent_lanes() {
  auto runtime = std::make_shared<hundun::chemistry::CanteraBackendRuntime>(
      config());
  HUNDUN_CHECK(runtime->composition().species.size() == 2U);
  HUNDUN_CHECK(runtime->composition().species[0].name == "A");
  HUNDUN_CHECK(runtime->composition().species[1].name == "B");
  HUNDUN_CHECK(runtime->mechanism_sha256() ==
               config().mechanism.sha256);

  hundun::chemistry::CanteraWorkspacePool pool(runtime, 2U);
  HUNDUN_CHECK(pool.workspace_count() == 2U);
  HUNDUN_CHECK(pool.workspaces_are_distinct());
  auto lane0 = hundun::chemistry::make_cantera_backend(config(), pool);
  auto lane1 = hundun::chemistry::make_cantera_backend(config(), pool);
  HUNDUN_CHECK(lane0->lane_index() == 0U);
  HUNDUN_CHECK(lane1->lane_index() == 1U);
  HUNDUN_CHECK(lane0->composition().fingerprint ==
               lane1->composition().fingerprint);

  bool missing = false;
  try {
    static_cast<void>(
        hundun::chemistry::make_cantera_backend(config(), pool));
  } catch (const std::runtime_error &) {
    missing = true;
  }
  HUNDUN_CHECK(missing);
}

void test_identity_mismatch_and_runtime_lifetime() {
  auto wrong_hash = config();
  wrong_hash.mechanism.sha256[0] = '0';
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::chemistry::CanteraBackendRuntime(wrong_hash));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto runtime = std::make_shared<hundun::chemistry::CanteraBackendRuntime>(
      config());
  hundun::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto changed = config();
  changed.species_names = {"B", "A"};
  rejected = false;
  try {
    static_cast<void>(
        hundun::chemistry::make_cantera_backend(changed, pool));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  runtime.reset();
  rejected = false;
  try {
    static_cast<void>(
        hundun::chemistry::make_cantera_backend(config(), pool));
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    mechanism_path = argv[1];
    test_exact_identity_and_independent_lanes();
    test_identity_mismatch_and_runtime_lifetime();
  });
}
