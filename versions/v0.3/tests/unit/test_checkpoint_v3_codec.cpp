// SPDX-License-Identifier: Apache-2.0

#include "src/flow_checkpoint_v3_detail.hpp"

#include "hundun/rt_error.hpp"
#include "tests/support/test_main.hpp"

#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

namespace {

using hundun::flow::CheckpointV3Presence;
using hundun::flow::IdealGasClosureState;
using hundun::flow::IdealGasPressureMode;
using hundun::flow::detail::CheckpointV3AuthorityGradient;
using hundun::flow::detail::CheckpointV3Manifest;
using hundun::flow::detail::CheckpointV3WaleIdentity;
using hundun::flow::detail::CheckpointV3WaleTransientField;

constexpr std::size_t kProfile1GoldenSize = 558U;
constexpr std::uint64_t kProfile1GoldenCrc =
    UINT64_C(8909770059348032994);
constexpr std::uint64_t kIbmSectionBytes = 324U;
constexpr std::uint64_t kWaleSectionBytes = 47U;
constexpr std::uint64_t kIdealGasSectionBytes = 30U;

CheckpointV3Manifest fixture() {
  CheckpointV3Manifest value;
  value.presence = CheckpointV3Presence::constant_static_ibm;
  value.rank_count = 2;
  value.process_grid = {2, 1, 1};
  value.payload_report_fingerprint = 21U;
  value.payload_manifest_crc64 = 22U;
  value.metadata = {3U, 0.003, 0.001, 0.001,
                    hundun::flow::MomentumTimeOrder::bdf2};
  value.control.proposed_next_dt_s = 0.001;
  value.control.last_retry_count = 1U;
  value.ranks = {
      {0,
       {{0, 0, 0}, {6, 12, 12}},
       {11U, 12U, 13U, 14U, 15U, 16U, 17U, 18U, 19U},
       "rank-000000.v3.bin",
       100U,
       200U,
       30U,
       true,
       true,
       {{7U, -0.0}, {9U, 1.25}},
       {{7U, -2.5}, {9U, 3.75}}},
      {1,
       {{6, 0, 0}, {12, 12, 12}},
       {21U, 22U, 23U, 24U, 25U, 26U, 27U, 28U, 29U},
       "rank-000001.v3.bin",
       110U,
       210U,
       31U,
       true,
       true,
       {{17U, 0.5}, {19U, 1.5}},
       {{17U, -1.5}, {19U, 2.5}}},
  };
  return value;
}

CheckpointV3WaleIdentity wale_identity() {
  return {0.5,
          0.9,
          0.7,
          UINT64_C(0x8ed19352cd41af76),
          1U,
          {CheckpointV3WaleTransientField::nu_t_m2_per_s,
           CheckpointV3WaleTransientField::mu_sgs_pa_s,
           CheckpointV3WaleTransientField::mu_eff_pa_s}};
}

CheckpointV3Manifest wale_fixture(CheckpointV3Presence presence) {
  auto value = fixture();
  value.presence = presence;
  value.wale_section_count = 1U;
  value.wale_section_bytes = kWaleSectionBytes;
  value.wale = wale_identity();
  if (presence == CheckpointV3Presence::constant_body_fitted_wale) {
    for (auto &rank : value.ranks) {
      rank.fingerprints.fill(0U);
      rank.history_available = false;
      rank.committed_available = false;
      rank.history.clear();
      rank.committed.clear();
    }
  } else {
    value.ibm_section_count = 1U;
    value.ibm_section_bytes = kIbmSectionBytes;
  }
  return value;
}

bool profile_has_ibm(CheckpointV3Presence presence) noexcept {
  return presence == CheckpointV3Presence::material_static_ibm ||
         presence == CheckpointV3Presence::material_static_ibm_wale ||
         presence == CheckpointV3Presence::ideal_gas_static_ibm ||
         presence == CheckpointV3Presence::ideal_gas_static_ibm_wale;
}

bool profile_has_wale(CheckpointV3Presence presence) noexcept {
  return presence == CheckpointV3Presence::material_body_fitted_wale ||
         presence == CheckpointV3Presence::material_static_ibm_wale ||
         presence == CheckpointV3Presence::ideal_gas_body_fitted_wale ||
         presence == CheckpointV3Presence::ideal_gas_static_ibm_wale;
}

bool profile_has_ideal_gas(CheckpointV3Presence presence) noexcept {
  return presence == CheckpointV3Presence::ideal_gas_static_ibm ||
         presence == CheckpointV3Presence::ideal_gas_body_fitted_wale ||
         presence == CheckpointV3Presence::ideal_gas_static_ibm_wale;
}

CheckpointV3Manifest density_fixture(CheckpointV3Presence presence) {
  auto value = fixture();
  value.presence = presence;
  if (profile_has_ibm(presence)) {
    value.ibm_section_count = 1U;
    value.ibm_section_bytes = kIbmSectionBytes;
  } else {
    for (auto &rank : value.ranks) {
      rank.fingerprints.fill(0U);
      rank.history_available = false;
      rank.committed_available = false;
      rank.history.clear();
      rank.committed.clear();
    }
  }
  if (profile_has_wale(presence)) {
    value.wale_section_count = 1U;
    value.wale_section_bytes = kWaleSectionBytes;
    value.wale = wale_identity();
  }
  if (profile_has_ideal_gas(presence)) {
    value.ideal_gas_section_count = 1U;
    value.ideal_gas_section_bytes = kIdealGasSectionBytes;
    value.ideal_gas =
        presence == CheckpointV3Presence::ideal_gas_body_fitted_wale
            ? IdealGasClosureState{IdealGasPressureMode::open_fixed, 101325.0,
                                   std::nullopt, 19U}
            : IdealGasClosureState{IdealGasPressureMode::closed_dynamic,
                                   101111.0, 2.5, 23U};
  }
  return value;
}

void expect_rejected(std::vector<std::uint8_t> bytes) {
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::flow::detail::decode_checkpoint_v3_manifest(bytes));
  } catch (const hundun::runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void expect_encode_rejected(const CheckpointV3Manifest &manifest) {
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::flow::detail::encode_checkpoint_v3_manifest(manifest));
  } catch (const hundun::runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void check_profile_1_golden_bytes() {
  const auto bytes =
      hundun::flow::detail::encode_checkpoint_v3_manifest(fixture());
  HUNDUN_CHECK(bytes.size() == kProfile1GoldenSize);
  HUNDUN_CHECK(hundun::flow::detail::checkpoint_v3_manifest_crc64(bytes) ==
               kProfile1GoldenCrc);
}

void check_wale_profile_roundtrips() {
  for (const auto presence : {
           CheckpointV3Presence::constant_body_fitted_wale,
           CheckpointV3Presence::constant_static_ibm_wale,
       }) {
    const auto expected = wale_fixture(presence);
    std::optional<CheckpointV3Manifest> decoded;
    try {
      const auto bytes =
          hundun::flow::detail::encode_checkpoint_v3_manifest(expected);
      decoded = hundun::flow::detail::decode_checkpoint_v3_manifest(bytes);
    } catch (const hundun::runtime::Error &) {
    }
    HUNDUN_CHECK(decoded.has_value());
    HUNDUN_CHECK(hundun::flow::detail::checkpoint_v3_manifest_equal(
        expected, *decoded));
  }
}

void check_density_profile_roundtrips() {
  for (const auto presence : {
           CheckpointV3Presence::material_static_ibm,
           CheckpointV3Presence::material_body_fitted_wale,
           CheckpointV3Presence::material_static_ibm_wale,
           CheckpointV3Presence::ideal_gas_static_ibm,
           CheckpointV3Presence::ideal_gas_body_fitted_wale,
           CheckpointV3Presence::ideal_gas_static_ibm_wale,
       }) {
    const auto expected = density_fixture(presence);
    std::optional<CheckpointV3Manifest> decoded;
    try {
      const auto bytes =
          hundun::flow::detail::encode_checkpoint_v3_manifest(expected);
      decoded = hundun::flow::detail::decode_checkpoint_v3_manifest(bytes);
    } catch (const hundun::runtime::Error &) {
    }
    HUNDUN_CHECK(decoded.has_value());
    HUNDUN_CHECK(hundun::flow::detail::checkpoint_v3_manifest_equal(
        expected, *decoded));
  }
}

void check_roundtrip() {
  const auto expected = fixture();
  const auto bytes =
      hundun::flow::detail::encode_checkpoint_v3_manifest(expected);
  HUNDUN_CHECK(!bytes.empty());
  const auto decoded =
      hundun::flow::detail::decode_checkpoint_v3_manifest(bytes);
  HUNDUN_CHECK(
      hundun::flow::detail::checkpoint_v3_manifest_equal(expected, decoded));
  HUNDUN_CHECK(hundun::flow::detail::checkpoint_v3_manifest_crc64(bytes) !=
               0U);
  HUNDUN_CHECK(std::signbit(decoded.ranks.front().history.front().value));
}

void check_empty_local_authority_roundtrip() {
  auto expected = fixture();
  expected.ranks[1U].history.clear();
  expected.ranks[1U].committed.clear();
  const auto bytes =
      hundun::flow::detail::encode_checkpoint_v3_manifest(expected);
  const auto decoded =
      hundun::flow::detail::decode_checkpoint_v3_manifest(bytes);
  HUNDUN_CHECK(
      hundun::flow::detail::checkpoint_v3_manifest_equal(expected, decoded));
  HUNDUN_CHECK(decoded.ranks[1U].history_available);
  HUNDUN_CHECK(decoded.ranks[1U].committed_available);
  HUNDUN_CHECK(decoded.ranks[1U].history.empty());
  HUNDUN_CHECK(decoded.ranks[1U].committed.empty());
}

void check_zero_crc_roundtrip() {
  auto expected = fixture();
  expected.payload_report_fingerprint = 0U;
  expected.payload_manifest_crc64 = 0U;
  expected.ranks[0U].fingerprints[0U] = 0U;
  expected.ranks[0U].state_crc64 = 0U;
  const auto bytes =
      hundun::flow::detail::encode_checkpoint_v3_manifest(expected);
  const auto decoded =
      hundun::flow::detail::decode_checkpoint_v3_manifest(bytes);
  HUNDUN_CHECK(
      hundun::flow::detail::checkpoint_v3_manifest_equal(expected, decoded));
}

void check_semantic_mutations() {
  auto changed = fixture();
  ++changed.ranks[0U].fingerprints[4U];
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(
      fixture(), changed));
  changed = fixture();
  changed.ranks[0U].committed[0U].value =
      std::nextafter(changed.ranks[0U].committed[0U].value,
                     -std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(
      fixture(), changed));

  changed = fixture();
  changed.presence = static_cast<CheckpointV3Presence>(0U);
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::flow::detail::encode_checkpoint_v3_manifest(changed));
  } catch (const hundun::runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  for (const auto &mutate : {
           std::function<void(CheckpointV3Manifest &)>([](auto &value) {
             value.ranks[0U].history[1U].link =
                 value.ranks[0U].history[0U].link;
           }),
           std::function<void(CheckpointV3Manifest &)>([](auto &value) {
             value.ranks[0U].committed[0U].link = 10U;
           }),
           std::function<void(CheckpointV3Manifest &)>([](auto &value) {
             value.ranks[1U].rank = 0;
           }),
           std::function<void(CheckpointV3Manifest &)>([](auto &value) {
             value.control.proposed_next_dt_s =
                 std::numeric_limits<double>::quiet_NaN();
           }),
           std::function<void(CheckpointV3Manifest &)>([](auto &value) {
             value.metadata.order =
                 static_cast<hundun::flow::MomentumTimeOrder>(255U);
           }),
           std::function<void(CheckpointV3Manifest &)>([](auto &value) {
             value.metadata.step = 1U;
             value.metadata.time_s = value.metadata.dt_s;
             value.metadata.previous_dt_s = value.metadata.dt_s;
             value.metadata.order = hundun::flow::MomentumTimeOrder::bdf2;
           }),
       }) {
    changed = fixture();
    mutate(changed);
    bool invalid = false;
    try {
      static_cast<void>(
          hundun::flow::detail::encode_checkpoint_v3_manifest(changed));
    } catch (const hundun::runtime::Error &) {
      invalid = true;
    }
    HUNDUN_CHECK(invalid);
  }
}

void check_profile_section_mutations() {
  auto changed = fixture();
  changed.wale_section_count = 1U;
  changed.wale_section_bytes = kWaleSectionBytes;
  changed.wale = wale_identity();
  expect_encode_rejected(changed);

  changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  changed.ibm_section_count = 1U;
  changed.ibm_section_bytes = kIbmSectionBytes;
  expect_encode_rejected(changed);
  changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  changed.ibm_section_bytes = 1U;
  expect_encode_rejected(changed);
  changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  changed.ranks.front().fingerprints.front() = 1U;
  expect_encode_rejected(changed);

  changed = wale_fixture(CheckpointV3Presence::constant_static_ibm_wale);
  changed.ibm_section_count = 0U;
  changed.ibm_section_bytes = 0U;
  expect_encode_rejected(changed);
  changed = wale_fixture(CheckpointV3Presence::constant_static_ibm_wale);
  changed.wale.reset();
  changed.wale_section_count = 0U;
  changed.wale_section_bytes = 0U;
  expect_encode_rejected(changed);

  changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  changed.wale->transient_schema_version = 0U;
  expect_encode_rejected(changed);
  changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  changed.wale->transient_fields[2U] =
      CheckpointV3WaleTransientField::nu_t_m2_per_s;
  expect_encode_rejected(changed);

  for (const auto &mutate : {
           std::function<void(CheckpointV3WaleIdentity &)>([](auto &value) {
             value.coefficient = std::nextafter(1.0e-6, 0.0);
           }),
           std::function<void(CheckpointV3WaleIdentity &)>([](auto &value) {
             value.coefficient =
                 std::nextafter(1.0, std::numeric_limits<double>::infinity());
           }),
           std::function<void(CheckpointV3WaleIdentity &)>([](auto &value) {
             value.turbulent_prandtl = std::nextafter(0.1, 0.0);
           }),
           std::function<void(CheckpointV3WaleIdentity &)>([](auto &value) {
             value.turbulent_prandtl = std::nextafter(
                 10.0, std::numeric_limits<double>::infinity());
           }),
           std::function<void(CheckpointV3WaleIdentity &)>([](auto &value) {
             value.turbulent_schmidt = std::nextafter(0.1, 0.0);
           }),
           std::function<void(CheckpointV3WaleIdentity &)>([](auto &value) {
             value.turbulent_schmidt = std::nextafter(
                 10.0, std::numeric_limits<double>::infinity());
           }),
       }) {
    changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
    mutate(*changed.wale);
    expect_encode_rejected(changed);
  }

  changed = wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  changed.presence = static_cast<CheckpointV3Presence>(255U);
  expect_encode_rejected(changed);
}

void check_density_profile_section_mutations() {
  auto changed = density_fixture(CheckpointV3Presence::material_static_ibm);
  changed.wale_section_count = 1U;
  changed.wale_section_bytes = kWaleSectionBytes;
  changed.wale = wale_identity();
  expect_encode_rejected(changed);

  changed =
      density_fixture(CheckpointV3Presence::material_body_fitted_wale);
  changed.ibm_section_count = 1U;
  changed.ibm_section_bytes = kIbmSectionBytes;
  expect_encode_rejected(changed);

  changed = density_fixture(CheckpointV3Presence::ideal_gas_static_ibm);
  changed.ideal_gas_section_count = 0U;
  changed.ideal_gas_section_bytes = 0U;
  changed.ideal_gas.reset();
  expect_encode_rejected(changed);

  changed = density_fixture(CheckpointV3Presence::material_static_ibm);
  changed.ideal_gas_section_count = 1U;
  changed.ideal_gas_section_bytes = kIdealGasSectionBytes;
  changed.ideal_gas = {IdealGasPressureMode::closed_dynamic, 101325.0, 2.5,
                       1U};
  expect_encode_rejected(changed);

  for (const auto &mutate : {
           std::function<void(IdealGasClosureState &)>([](auto &state) {
             state.mode = static_cast<IdealGasPressureMode>(255U);
           }),
           std::function<void(IdealGasClosureState &)>([](auto &state) {
             state.thermodynamic_pressure_pa = 0.0;
           }),
           std::function<void(IdealGasClosureState &)>([](auto &state) {
             state.target_mass_kg.reset();
           }),
           std::function<void(IdealGasClosureState &)>([](auto &state) {
             state.target_mass_kg = -1.0;
           }),
           std::function<void(IdealGasClosureState &)>([](auto &state) {
             state.revision = std::numeric_limits<std::uint64_t>::max();
           }),
       }) {
    changed = density_fixture(CheckpointV3Presence::ideal_gas_static_ibm);
    mutate(*changed.ideal_gas);
    expect_encode_rejected(changed);
  }

  changed =
      density_fixture(CheckpointV3Presence::ideal_gas_body_fitted_wale);
  changed.ideal_gas->target_mass_kg = 2.5;
  expect_encode_rejected(changed);
}

void check_ideal_gas_identity_mutations() {
  const auto original =
      density_fixture(CheckpointV3Presence::ideal_gas_static_ibm_wale);
  auto changed = original;
  changed.ideal_gas->thermodynamic_pressure_pa = std::nextafter(
      changed.ideal_gas->thermodynamic_pressure_pa,
      std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(original,
                                                                    changed));
  changed = original;
  changed.ideal_gas->target_mass_kg = std::nextafter(
      *changed.ideal_gas->target_mass_kg,
      std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(original,
                                                                    changed));
  changed = original;
  ++changed.ideal_gas->revision;
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(original,
                                                                    changed));
}

void check_wale_transient_attempt_identity_is_rejected() {
  auto bytes = hundun::flow::detail::encode_checkpoint_v3_manifest(
      wale_fixture(CheckpointV3Presence::constant_body_fitted_wale));
  constexpr std::uint64_t forbidden_attempt_identity =
      UINT64_C(0xdecafbad12345678);
  for (std::size_t byte = 0U; byte < sizeof(forbidden_attempt_identity);
       ++byte) {
    bytes.push_back(static_cast<std::uint8_t>(
        forbidden_attempt_identity >> (8U * byte)));
  }
  expect_rejected(std::move(bytes));
}

void check_wale_identity_mutations() {
  const auto original =
      wale_fixture(CheckpointV3Presence::constant_body_fitted_wale);
  auto changed = original;
  changed.wale->coefficient = std::nextafter(
      changed.wale->coefficient, std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(original,
                                                                    changed));
  changed = original;
  ++changed.wale->numerical_config_crc64;
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(original,
                                                                    changed));
  changed = original;
  changed.wale->transient_fields[0U] =
      CheckpointV3WaleTransientField::mu_eff_pa_s;
  HUNDUN_CHECK(!hundun::flow::detail::checkpoint_v3_manifest_equal(original,
                                                                    changed));
}

void check_byte_mutations() {
  const auto bytes =
      hundun::flow::detail::encode_checkpoint_v3_manifest(fixture());
  auto truncated = bytes;
  truncated.pop_back();
  expect_rejected(std::move(truncated));
  auto trailing = bytes;
  trailing.push_back(0U);
  expect_rejected(std::move(trailing));
  auto wrong_magic = bytes;
  wrong_magic.front() ^= 0xffU;
  expect_rejected(std::move(wrong_magic));
  auto wrong_presence = bytes;
  HUNDUN_CHECK(
      hundun::flow::detail::checkpoint_v3_manifest_presence_offset() <
      wrong_presence.size());
  wrong_presence[hundun::flow::detail::checkpoint_v3_manifest_presence_offset()] =
      0U;
  expect_rejected(std::move(wrong_presence));
}

} // namespace

int main() {
  return hundun::test::run([] {
    check_profile_1_golden_bytes();
    check_wale_profile_roundtrips();
    check_density_profile_roundtrips();
    check_roundtrip();
    check_empty_local_authority_roundtrip();
    check_zero_crc_roundtrip();
    check_semantic_mutations();
    check_profile_section_mutations();
    check_density_profile_section_mutations();
    check_wale_transient_attempt_identity_is_rejected();
    check_wale_identity_mutations();
    check_ideal_gas_identity_mutations();
    check_byte_mutations();
  });
}
