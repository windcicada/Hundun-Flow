// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/chem_composition.hpp"
#include "hundun/chem_reports.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using hundun::chemistry::CompositionIdentity;
using hundun::chemistry::SpeciesIdentity;

CompositionIdentity two_species_identity() {
  CompositionIdentity identity;
  identity.element_names = {"C", "H"};
  identity.species = {
      SpeciesIdentity{"CH4", 16.043, {1, 4}},
      SpeciesIdentity{"H2", 2.016, {0, 2}},
  };
  identity.fingerprint =
      hundun::chemistry::composition_identity_fingerprint(identity);
  return identity;
}

template <class Function> void expect_invalid(Function &&function) {
  bool rejected = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_all_species_identity_and_ordered_fingerprint() {
  auto identity = two_species_identity();
  hundun::chemistry::validate_composition_identity(identity);
  HUNDUN_CHECK(identity.species.size() == 2U);
  HUNDUN_CHECK(identity.species[0].name == "CH4");
  HUNDUN_CHECK(identity.species[1].name == "H2");
  HUNDUN_CHECK(identity.species[0].element_counts ==
               std::vector<std::int32_t>({1, 4}));

  auto reordered = identity;
  std::swap(reordered.species[0], reordered.species[1]);
  HUNDUN_CHECK(
      hundun::chemistry::composition_identity_fingerprint(reordered) !=
      identity.fingerprint);

  auto missing_last_species = identity;
  missing_last_species.species.pop_back();
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(missing_last_species);
  });

  auto changed_elements = identity;
  changed_elements.species[0].element_counts[1] = 3;
  HUNDUN_CHECK(
      hundun::chemistry::composition_identity_fingerprint(changed_elements) !=
      identity.fingerprint);
}

void test_identity_validation_rejects_invalid_shape_and_values() {
  auto identity = two_species_identity();

  auto empty_elements = identity;
  empty_elements.element_names.clear();
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(empty_elements);
  });

  auto wrong_shape = identity;
  wrong_shape.species[1].element_counts.pop_back();
  wrong_shape.fingerprint =
      hundun::chemistry::composition_identity_fingerprint(wrong_shape);
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(wrong_shape);
  });

  auto negative_count = identity;
  negative_count.species[0].element_counts[0] = -1;
  negative_count.fingerprint =
      hundun::chemistry::composition_identity_fingerprint(negative_count);
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(negative_count);
  });

  auto zero_weight = identity;
  zero_weight.species[0].molecular_weight_kg_per_kmol = 0.0;
  zero_weight.fingerprint =
      hundun::chemistry::composition_identity_fingerprint(zero_weight);
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(zero_weight);
  });

  auto non_finite_weight = identity;
  non_finite_weight.species[0].molecular_weight_kg_per_kmol =
      std::numeric_limits<double>::infinity();
  non_finite_weight.fingerprint =
      hundun::chemistry::composition_identity_fingerprint(non_finite_weight);
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(non_finite_weight);
  });

  auto duplicate_species = identity;
  duplicate_species.species[1].name = duplicate_species.species[0].name;
  duplicate_species.fingerprint =
      hundun::chemistry::composition_identity_fingerprint(duplicate_species);
  expect_invalid([&] {
    hundun::chemistry::validate_composition_identity(duplicate_species);
  });
}

void test_total_thermochemical_enthalpy_uses_all_species() {
  const std::vector<double> mass_fractions{0.25, 0.75};
  const std::vector<double> sensible_enthalpies{100.0, 200.0};
  const std::vector<double> formation_enthalpies{-400.0, 20.0};
  const std::vector<double> total_species_enthalpies{
      sensible_enthalpies[0] + formation_enthalpies[0],
      sensible_enthalpies[1] + formation_enthalpies[1]};

  const double h_tc =
      hundun::chemistry::mixture_total_thermochemical_enthalpy_j_per_kg(
          mass_fractions, total_species_enthalpies);
  HUNDUN_CHECK_NEAR(h_tc, 90.0, 1.0e-14);
  HUNDUN_CHECK(h_tc !=
               hundun::chemistry::mixture_total_thermochemical_enthalpy_j_per_kg(
                   mass_fractions, sensible_enthalpies));

  hundun::chemistry::ThermochemicalPoint point;
  point.p0_pa = 101325.0;
  point.h_tc_j_per_kg = h_tc;
  point.mass_fractions = mass_fractions;
  HUNDUN_CHECK(point.mass_fractions.size() == 2U);
  HUNDUN_CHECK(point.h_tc_j_per_kg == 90.0);

  expect_invalid([&] {
    static_cast<void>(
        hundun::chemistry::mixture_total_thermochemical_enthalpy_j_per_kg(
            std::vector<double>{1.0}, total_species_enthalpies));
  });
}

void test_chemistry_status_and_report_contract() {
  using hundun::chemistry::ChemistryStatus;
  HUNDUN_CHECK(static_cast<std::uint32_t>(ChemistryStatus::success) == 0U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(ChemistryStatus::invalid_input) == 1U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(
                   ChemistryStatus::composition_mismatch) == 2U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(
                   ChemistryStatus::state_inversion_failure) == 3U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(
                   ChemistryStatus::integration_failure) == 4U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(
                   ChemistryStatus::non_finite_output) == 5U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(
                   ChemistryStatus::conservation_failure) == 6U);
  HUNDUN_CHECK(static_cast<std::uint32_t>(
                   ChemistryStatus::workspace_failure) == 7U);

  hundun::chemistry::ChemistryIntervalReport report;
  report.status = ChemistryStatus::success;
  report.integrated_rho_y_delta_kg_per_m3 = {-0.5, 0.5};
  HUNDUN_CHECK(report.succeeded());
  report.status = ChemistryStatus::integration_failure;
  HUNDUN_CHECK(!report.succeeded());
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_all_species_identity_and_ordered_fingerprint();
    test_identity_validation_rejects_invalid_shape_and_values();
    test_total_thermochemical_enthalpy_uses_all_species();
    test_chemistry_status_and_report_contract();
  });
}
