// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_cantera.hpp"

#include "../support/chemistry_test_support.hpp"

#include <filesystem>
#include <memory>
#include <limits>
#include <stdexcept>

namespace {

std::filesystem::path mechanism_path;

hundun::v04::chemistry::CanteraBackendConfig config() {
  hundun::v04::chemistry::CanteraBackendConfig value;
  value.mechanism.file = mechanism_path;
  value.mechanism.sha256 =
      "c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee";
  value.mechanism.phase = "synthetic-gas";
  value.species_names = {"A", "B"};
  value.chemistry.relative_tolerance = 1.0e-10;
  value.chemistry.absolute_tolerance = 1.0e-18;
  value.chemistry.maximum_internal_steps = 5000;
  return value;
}

void test_exact_identity_and_independent_lanes() {
  auto runtime = std::make_shared<hundun::v04::chemistry::CanteraBackendRuntime>(
      config());
  HUNDUN_CHECK(runtime->composition().species.size() == 2U);
  HUNDUN_CHECK(runtime->composition().species[0].name == "A");
  HUNDUN_CHECK(runtime->composition().species[1].name == "B");
  HUNDUN_CHECK(runtime->mechanism_sha256() ==
               config().mechanism.sha256);

  hundun::v04::chemistry::CanteraWorkspacePool pool(runtime, 2U);
  HUNDUN_CHECK(pool.workspace_count() == 2U);
  HUNDUN_CHECK(pool.workspaces_are_distinct());
  auto lane0 = hundun::v04::chemistry::make_cantera_backend(config(), pool);
  auto lane1 = hundun::v04::chemistry::make_cantera_backend(config(), pool);
  HUNDUN_CHECK(lane0->lane_index() == 0U);
  HUNDUN_CHECK(lane1->lane_index() == 1U);
  HUNDUN_CHECK(lane0->composition().fingerprint ==
               lane1->composition().fingerprint);

  bool missing = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::make_cantera_backend(config(), pool));
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
        hundun::v04::chemistry::CanteraBackendRuntime(wrong_hash));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto runtime = std::make_shared<hundun::v04::chemistry::CanteraBackendRuntime>(
      config());
  hundun::v04::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto changed = config();
  changed.species_names = {"B", "A"};
  rejected = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::make_cantera_backend(changed, pool));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  runtime.reset();
  rejected = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::make_cantera_backend(config(), pool));
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_controls_are_validated_without_product_schema() {
  auto runtime = std::make_shared<hundun::v04::chemistry::CanteraBackendRuntime>(
      config());
  hundun::v04::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  for (int mutation = 0; mutation < 4; ++mutation) {
    auto invalid = config();
    if (mutation == 0) invalid.chemistry.relative_tolerance = 0.0;
    if (mutation == 1) invalid.chemistry.absolute_tolerance = -1.0;
    if (mutation == 2) invalid.chemistry.maximum_internal_steps = 0;
    if (mutation == 3)
      invalid.chemistry.relative_tolerance =
          std::numeric_limits<double>::quiet_NaN();
    bool rejected = false;
    try {
      static_cast<void>(hundun::v04::chemistry::make_cantera_backend(invalid, pool));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  // Rejected controls must not consume a lane.
  auto valid = hundun::v04::chemistry::make_cantera_backend(config(), pool);
  HUNDUN_CHECK(valid->lane_index() == 0U);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    mechanism_path = argv[1];
    test_exact_identity_and_independent_lanes();
    test_identity_mismatch_and_runtime_lifetime();
    test_controls_are_validated_without_product_schema();
  });
}
