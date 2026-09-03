// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/flow_checkpoint_v4.hpp"
#include "flow_checkpoint_v4_detail.hpp"

#include "tests/support/test_main.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace {

using hundun::flow::CheckpointSectionPresence;
using hundun::flow::CheckpointSectionProvider;
using hundun::flow::CheckpointSectionRegistry;
using hundun::flow::CheckpointV4Manifest;
using hundun::flow::CheckpointV4Section;

bool valid_provider(const void *, std::string &) noexcept { return true; }
bool invalid_provider(const void *, std::string &message) noexcept {
  message = "invalid fixture";
  return false;
}
std::vector<std::uint8_t> write_local(const void *) { return {1U, 2U}; }

CheckpointSectionProvider provider(std::uint32_t id,
                                   CheckpointSectionPresence presence,
                                   bool valid = true) {
  return {{id, 1U, presence, "fixture-" + std::to_string(id)}, nullptr,
          valid ? &valid_provider : &invalid_provider, &write_local};
}

void test_atomic_registration_and_ordering() {
  CheckpointSectionRegistry registry;
  registry.register_providers(
      {provider(20U, CheckpointSectionPresence::optional),
       provider(10U, CheckpointSectionPresence::mandatory)});
  HUNDUN_CHECK(registry.providers().size() == 2U);
  HUNDUN_CHECK(registry.providers()[0].descriptor.id == 10U);
  HUNDUN_CHECK(registry.providers()[1].descriptor.id == 20U);

  bool rejected = false;
  try {
    registry.register_providers(
        {provider(30U, CheckpointSectionPresence::optional),
         provider(20U, CheckpointSectionPresence::mandatory)});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(registry.providers().size() == 2U);

  rejected = false;
  try {
    registry.register_provider(
        provider(30U, CheckpointSectionPresence::optional, false));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(!registry.find(30U));
}

void test_manifest_presence_and_optional_skip() {
  CheckpointSectionRegistry registry;
  registry.register_provider(
      provider(10U, CheckpointSectionPresence::mandatory));
  CheckpointV4Manifest manifest;
  manifest.sections = {
      CheckpointV4Section{10U, 1U, CheckpointSectionPresence::mandatory, 2U},
      CheckpointV4Section{99U, 1U, CheckpointSectionPresence::optional, 4U}};
  const auto compatibility =
      hundun::flow::validate_checkpoint_v4_manifest(manifest, registry);
  HUNDUN_CHECK(compatibility.restored_section_ids ==
               std::vector<std::uint32_t>{10U});
  HUNDUN_CHECK(compatibility.skipped_optional_section_ids ==
               std::vector<std::uint32_t>{99U});

  manifest.sections[1].presence = CheckpointSectionPresence::mandatory;
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::flow::validate_checkpoint_v4_manifest(manifest, registry));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  manifest.sections[1] = manifest.sections[0];
  rejected = false;
  try {
    static_cast<void>(
        hundun::flow::validate_checkpoint_v4_manifest(manifest, registry));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

hundun::flow::detail::ReactingCheckpointV4Data reacting_data() {
  return {77U,
          {"A", "B"},
          {0.7, 0.3, 0.6, 0.4},
          {0.65, 0.35, 0.55, 0.45},
          {10.0, 11.0},
          {12.0, 13.0},
          100000.0,
          101325.0,
          9U,
          0.09,
          0.01,
          2U,
          "cantera-3.2.0",
          std::string(64U, 'a'),
          "gas"};
}

bool same(const hundun::flow::detail::ReactingCheckpointV4Data &left,
          const hundun::flow::detail::ReactingCheckpointV4Data &right) {
  return left.composition_fingerprint == right.composition_fingerprint &&
         left.species_names == right.species_names &&
         left.history_rho_y_kg_per_m3 == right.history_rho_y_kg_per_m3 &&
         left.committed_rho_y_kg_per_m3 == right.committed_rho_y_kg_per_m3 &&
         left.history_rho_h_tc_j_per_m3 == right.history_rho_h_tc_j_per_m3 &&
         left.committed_rho_h_tc_j_per_m3 ==
             right.committed_rho_h_tc_j_per_m3 &&
         left.history_p0_pa == right.history_p0_pa &&
         left.committed_p0_pa == right.committed_p0_pa &&
         left.step == right.step && left.time_s == right.time_s &&
         left.previous_dt_s == right.previous_dt_s &&
         left.bdf_order == right.bdf_order &&
         left.backend_id == right.backend_id &&
         left.mechanism_sha256 == right.mechanism_sha256 &&
         left.mechanism_phase == right.mechanism_phase;
}

void test_reacting_codecs_validate_then_publish() {
  const auto expected = reacting_data();
  auto sections =
      hundun::flow::detail::encode_reacting_checkpoint_v4(expected);
  HUNDUN_CHECK(sections.size() == 5U);
  HUNDUN_CHECK(sections.front().descriptor.id ==
               hundun::flow::detail::kReactingCompositionSection);
  HUNDUN_CHECK(sections.back().descriptor.id ==
               hundun::flow::detail::kReactingBackendSection);
  auto restored = reacting_data();
  restored.step = 999U;
  std::string message;
  HUNDUN_CHECK(hundun::flow::detail::restore_reacting_checkpoint_v4(
      sections, expected.composition_fingerprint,
      expected.mechanism_sha256, restored, message));
  HUNDUN_CHECK(same(restored, expected));

  const auto unchanged = restored;
  sections[1].payload.back() ^= 1U;
  HUNDUN_CHECK(!hundun::flow::detail::restore_reacting_checkpoint_v4(
      sections, expected.composition_fingerprint,
      expected.mechanism_sha256, restored, message));
  HUNDUN_CHECK(same(restored, unchanged));
  sections = hundun::flow::detail::encode_reacting_checkpoint_v4(expected);
  HUNDUN_CHECK(!hundun::flow::detail::restore_reacting_checkpoint_v4(
      sections, expected.composition_fingerprint,
      std::string(64U, 'b'), restored, message));
  HUNDUN_CHECK(same(restored, unchanged));

  auto omitted_species = expected;
  omitted_species.committed_rho_y_kg_per_m3.pop_back();
  bool rejected = false;
  try {
    static_cast<void>(hundun::flow::detail::encode_reacting_checkpoint_v4(
        omitted_species));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_atomic_registration_and_ordering();
    test_manifest_presence_and_optional_skip();
    test_reacting_codecs_validate_then_publish();
  });
}
