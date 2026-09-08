// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_breakup_detail.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

constexpr double kPi = 3.141592653589793238462643383279502884;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool near(double left, double right, double relative = 1.0e-12,
          double absolute = 1.0e-18) {
  return std::abs(left - right) <=
         absolute + relative * (std::abs(left) + std::abs(right));
}

SprayParcelState parent_parcel(double multiplicity = 100.0) {
  SprayParcelState parent;
  parent.id = {UINT64_C(0x123456789abcdef0),
               UINT64_C(0x0fedcba987654321)};
  parent.position_m = {1.0, 2.0, 3.0};
  parent.velocity_m_per_s = {20.0, -3.0, 4.0};
  parent.droplet_diameter_m = 1.0e-4;
  parent.droplet_mass_kg =
      800.0 * kPi / 6.0 * std::pow(parent.droplet_diameter_m, 3.0);
  parent.multiplicity = multiplicity;
  parent.temperature_k = 420.0;
  parent.liquid_material_fingerprint = 91U;
  parent.owner_global_cell = 17U;
  parent.age_s = 0.025;
  return parent;
}

TabBreakupReport breakup_trigger() {
  TabBreakupInput input;
  input.initial_deformation = 1.0;
  input.initial_deformation_rate_per_s = 2.0;
  input.relative_speed_m_per_s = 10.0;
  input.gas_density_kg_per_m3 = 1.0;
  input.liquid_density_kg_per_m3 = 800.0;
  input.liquid_viscosity_pa_s = 8.0e-4;
  input.surface_tension_n_per_m = 0.025;
  input.droplet_radius_m = 5.0e-5;
  input.breakup_threshold = 1.0;
  input.duration_s = 1.0e-4;
  input.coefficients = {2.0 / 3.0, 0.5, 8.0, 10.0};
  return evaluate_tab_breakup(input);
}

BreakupChildInput valid_input(double multiplicity = 100.0) {
  BreakupChildInput input;
  input.parent = parent_parcel(multiplicity);
  input.tab_trigger = breakup_trigger();
  input.accepted_step = 37U;
  input.breakup_ordinal = 9U;
  input.parameters.child_parcel_count = 4U;
  input.parameters.maximum_child_parcels = 8U;
  input.parameters.supplied_uniform_child_diameter_m = 5.0e-5;
  input.parameters.liquid_density_kg_per_m3 = 800.0;
  input.parameters.surface_tension_n_per_m = 0.025;
  input.parameters.liquid_absolute_thermochemical_enthalpy_j_per_kg =
      -2.0e5;
  input.parameters.supplied_deformation_energy_per_parent_droplet_j = 0.0;
  return input;
}

bool same_parcel(const SprayParcelState& left,
                 const SprayParcelState& right) {
  return left.id == right.id && left.position_m == right.position_m &&
         left.velocity_m_per_s == right.velocity_m_per_s &&
         left.droplet_mass_kg == right.droplet_mass_kg &&
         left.droplet_diameter_m == right.droplet_diameter_m &&
         left.multiplicity == right.multiplicity &&
         left.temperature_k == right.temperature_k &&
         left.liquid_material_fingerprint ==
             right.liquid_material_fingerprint &&
         left.owner_global_cell == right.owner_global_cell &&
         left.age_s == right.age_s;
}

bool canonical_candidate(const BreakupChildCandidate& candidate) {
  if (candidate.available || candidate.child_count != 0U) return false;
  for (const SprayParcelState& child : candidate.children) {
    if (!same_parcel(child, SprayParcelState{})) return false;
  }
  return true;
}

bool same_candidate(const BreakupChildCandidate& left,
                    const BreakupChildCandidate& right) {
  if (left.available != right.available ||
      left.child_count != right.child_count) {
    return false;
  }
  for (std::uint32_t index = 0U; index < left.child_count; ++index) {
    if (!same_parcel(left.children[index], right.children[index])) {
      return false;
    }
  }
  return true;
}

bool test_conservative_uniform_split() {
  const BreakupChildInput input = valid_input();
  const SprayParcelState parent_before = input.parent;
  const BreakupChildReport report =
      generate_supplied_diameter_children(input);

  bool passed = true;
  passed &= expect(report.succeeded() &&
                       report.model_id ==
                           "conservative_supplied_uniform_diameter_split_v1" &&
                       report.trigger_model_id ==
                           "tab_linear_oscillator_v1" &&
                       report.rng_model_id ==
                           "parcel_counter_splitmix64_v1:breakup_child" &&
                       report.parent_id == input.parent.id &&
                       report.accepted_step == input.accepted_step &&
                       report.breakup_ordinal == input.breakup_ordinal &&
                       report.candidate.child_count == 4U,
                   "report freezes supplied-diameter, TAB-trigger and RNG identities");

  const double expected_child_mass =
      input.parameters.liquid_density_kg_per_m3 * kPi / 6.0 *
      std::pow(input.parameters.supplied_uniform_child_diameter_m, 3.0);
  const double expected_total_multiplicity =
      input.parent.multiplicity *
      std::pow(input.parent.droplet_diameter_m /
                   input.parameters.supplied_uniform_child_diameter_m,
               3.0);
  const double expected_each_multiplicity =
      expected_total_multiplicity / 4.0;
  for (std::uint32_t index = 0U; index < 4U; ++index) {
    const SprayParcelState& child = report.candidate.children[index];
    passed &= expect(child.id != input.parent.id &&
                         near(child.droplet_mass_kg, expected_child_mass,
                              1.0e-14, 1.0e-24) &&
                         child.droplet_diameter_m ==
                             input.parameters.supplied_uniform_child_diameter_m &&
                         near(child.multiplicity,
                              expected_each_multiplicity) &&
                         child.position_m == input.parent.position_m &&
                         child.velocity_m_per_s == input.parent.velocity_m_per_s &&
                         child.temperature_k == input.parent.temperature_k &&
                         child.liquid_material_fingerprint ==
                             input.parent.liquid_material_fingerprint &&
                         child.owner_global_cell ==
                             input.parent.owner_global_cell &&
                         child.age_s == input.parent.age_s,
                     "child inherits parent material/state and supplied diameter");
    for (std::uint32_t other = 0U; other < index; ++other) {
      passed &= expect(child.id != report.candidate.children[other].id,
                       "every child ID is distinct");
    }
  }
  passed &= expect(same_parcel(input.parent, parent_before),
                   "pure generation never mutates or deletes the parent");

  const BreakupConservationReport& budget = report.conservation;
  passed &= expect(near(budget.mass_residual_kg, 0.0) &&
                       near(budget.momentum_residual_kg_m_per_s[0U], 0.0) &&
                       near(budget.momentum_residual_kg_m_per_s[1U], 0.0) &&
                       near(budget.momentum_residual_kg_m_per_s[2U], 0.0) &&
                       near(budget.bulk_kinetic_energy_residual_j, 0.0) &&
                       near(budget.thermochemical_enthalpy_residual_j, 0.0) &&
                       budget.parent_temperature_k == 420.0 &&
                       budget.minimum_child_temperature_k == 420.0 &&
                       budget.maximum_child_temperature_k == 420.0,
                   "mass, momentum, bulk kinetic, h_tc and temperature budgets close");
  passed &= expect(
      near(budget.children_total_multiplicity,
           expected_total_multiplicity) &&
          near(budget.multiplicity_residual, 0.0) &&
          budget.children_total_multiplicity > budget.parent_multiplicity,
      "multiplicity increases by the explicit volume ratio rather than being conserved");
  passed &= expect(
      budget.surface_energy_increase_j > 0.0 &&
          budget.supplied_deformation_energy_j == 0.0 &&
          budget.unassigned_deformation_energy_j < 0.0 &&
          budget.energy_disposition ==
              BreakupEnergyDisposition::supplied_deformation_energy_deficit,
      "surface/deformation energy is reported separately without a false total-energy claim");
  return passed;
}

bool test_breakup_rng_domain_and_retry_determinism() {
  const BreakupChildInput input = valid_input();
  const BreakupChildReport first =
      generate_supplied_diameter_children(input);
  const BreakupChildReport retry =
      generate_supplied_diameter_children(input);
  bool passed = true;
  passed &= expect(first.succeeded() && retry.succeeded() &&
                       same_candidate(first.candidate, retry.candidate),
                   "retry from the same accepted address is bitwise deterministic");

  ParcelRandomAddress address;
  address.seed = input.parent.id.high;
  address.accepted_step = input.accepted_step;
  address.stream_identity = input.parent.id.low;
  address.ordinal = input.breakup_ordinal;
  address.purpose = ParcelRandomPurpose::breakup_child;
  passed &= expect(
      first.candidate.children[0U].id.high ==
              parcel_random_u64(address, 0U) &&
          first.candidate.children[0U].id.low ==
              parcel_random_u64(address, 1U),
      "child ID uses the full parent ID, accepted step, breakup ordinal and breakup domain");
  address.purpose = ParcelRandomPurpose::injection_direction;
  passed &= expect(first.candidate.children[0U].id.high !=
                       parcel_random_u64(address, 0U),
                   "breakup IDs do not consume the injector RNG domain");

  BreakupChildInput next_step = input;
  ++next_step.accepted_step;
  BreakupChildInput next_ordinal = input;
  ++next_ordinal.breakup_ordinal;
  BreakupChildInput changed_high = input;
  ++changed_high.parent.id.high;
  BreakupChildInput changed_low = input;
  ++changed_low.parent.id.low;
  const BreakupChildReport step_report =
      generate_supplied_diameter_children(next_step);
  const BreakupChildReport ordinal_report =
      generate_supplied_diameter_children(next_ordinal);
  const BreakupChildReport high_report =
      generate_supplied_diameter_children(changed_high);
  const BreakupChildReport low_report =
      generate_supplied_diameter_children(changed_low);
  passed &= expect(
      step_report.succeeded() && ordinal_report.succeeded() &&
          high_report.succeeded() && low_report.succeeded() &&
          step_report.candidate.children[0U].id !=
              first.candidate.children[0U].id &&
          ordinal_report.candidate.children[0U].id !=
              first.candidate.children[0U].id &&
          high_report.candidate.children[0U].id !=
              first.candidate.children[0U].id &&
          low_report.candidate.children[0U].id !=
              first.candidate.children[0U].id,
      "every frozen address component changes the child ID domain");
  return passed;
}

bool test_energy_disposition_and_multiplicity_scaling() {
  BreakupChildInput deficit_input = valid_input();
  const BreakupChildReport deficit =
      generate_supplied_diameter_children(deficit_input);
  BreakupChildInput balanced_input = deficit_input;
  balanced_input.parameters.supplied_deformation_energy_per_parent_droplet_j =
      deficit.conservation.surface_energy_increase_j /
      deficit_input.parent.multiplicity;
  const BreakupChildReport balanced =
      generate_supplied_diameter_children(balanced_input);
  const BreakupChildReport doubled =
      generate_supplied_diameter_children(valid_input(200.0));

  bool passed = true;
  passed &= expect(
      deficit.succeeded() && balanced.succeeded() &&
          balanced.conservation.energy_disposition ==
              BreakupEnergyDisposition::balanced &&
          near(balanced.conservation.unassigned_deformation_energy_j, 0.0),
      "explicit deformation energy can balance surface creation without hidden conversion");
  passed &= expect(
      doubled.succeeded() &&
          near(doubled.conservation.parent_total_mass_kg,
               2.0 * deficit.conservation.parent_total_mass_kg) &&
          near(doubled.conservation.children_total_mass_kg,
               2.0 * deficit.conservation.children_total_mass_kg) &&
          near(doubled.conservation.children_total_multiplicity,
               2.0 * deficit.conservation.children_total_multiplicity) &&
          doubled.candidate.children[0U].id ==
              deficit.candidate.children[0U].id,
      "parcel multiplicity scales extensive budgets once but not deterministic IDs");
  return passed;
}

bool test_failures_publish_no_partial_children() {
  BreakupChildInput no_breakup = valid_input();
  no_breakup.tab_trigger.breakup_requested = false;
  const BreakupChildReport no_breakup_report =
      generate_supplied_diameter_children(no_breakup);

  BreakupChildInput too_large = valid_input();
  too_large.parameters.supplied_uniform_child_diameter_m =
      too_large.parent.droplet_diameter_m;
  const BreakupChildReport too_large_report =
      generate_supplied_diameter_children(too_large);

  BreakupChildInput overflow = valid_input();
  overflow.parameters.child_parcel_count = 9U;
  const BreakupChildReport overflow_report =
      generate_supplied_diameter_children(overflow);

  BreakupChildInput inconsistent = valid_input();
  inconsistent.parent.droplet_mass_kg *= 1.1;
  const BreakupChildReport inconsistent_report =
      generate_supplied_diameter_children(inconsistent);

  BreakupChildInput non_finite = valid_input();
  non_finite.parameters.surface_tension_n_per_m =
      std::numeric_limits<double>::quiet_NaN();
  const BreakupChildReport non_finite_report =
      generate_supplied_diameter_children(non_finite);

  bool passed = true;
  passed &= expect(
      no_breakup_report.status == BreakupChildStatus::breakup_not_requested &&
          canonical_candidate(no_breakup_report.candidate),
      "a non-triggering TAB report cannot stage children");
  passed &= expect(too_large_report.status ==
                           BreakupChildStatus::invalid_input &&
                       canonical_candidate(too_large_report.candidate),
                   "a supplied diameter that does not break the parent is rejected");
  passed &= expect(overflow_report.status ==
                           BreakupChildStatus::invalid_input &&
                       canonical_candidate(overflow_report.candidate),
                   "configured child bound is enforced before generation");
  passed &= expect(
      inconsistent_report.status ==
              BreakupChildStatus::inconsistent_parent_geometry &&
          canonical_candidate(inconsistent_report.candidate),
      "parent mass/diameter/density mismatch has no partial output");
  passed &= expect(non_finite_report.status ==
                           BreakupChildStatus::invalid_input &&
                       canonical_candidate(non_finite_report.candidate),
                   "non-finite input has canonical empty candidate and budget");
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  TabRepresentativeSplitInput tab;
  tab.parent = parent_parcel();
  tab.tab_trigger = breakup_trigger();
  tab.tab_trigger.candidate.deformation_rate_per_s = 0.0;
  tab.liquid_density_kg_per_m3 = 800.0;
  tab.surface_tension_n_per_m = 0.025;
  tab.liquid_absolute_thermochemical_enthalpy_j_per_kg = -2.0e5;
  const auto bag = generate_tab_representative_children(tab);
  passed &= expect(bag.succeeded() &&
      near(bag.representative_diameter_m, 4.2857142857142857e-5, 1e-13, 1e-20) &&
      near(bag.total_energy_residual_j, 0.0, 0.0, 1e-18),
      "TAB bag limit r/r32=7/3 closes surface and deformation energy");
  // Worked SI example: r=50 um, rho=800, sigma=.025, ydot=20000/s.
  // rho*r^3*ydot^2/(8*sigma)=.2; d32=3.94736842105263e-5 m;
  // transverse speed=.5 m/s, Ekin=5.235987755982989e-9 J (100 drops).
  tab.tab_trigger.candidate.deformation_rate_per_s = 20000.0;
  tab.child_parcel_count = 4;
  const auto moving = generate_tab_representative_children(tab);
  const auto retry = generate_tab_representative_children(tab);
  passed &= expect(moving.succeeded() &&
      near(moving.representative_diameter_m, 3.9473684210526316e-5, 1e-13, 1e-20) &&
      near(moving.transverse_speed_m_per_s, 0.5) &&
      near(moving.transverse_kinetic_energy_j, 5.235987755982989e-9, 1e-13, 1e-22) &&
      near(moving.total_energy_residual_j, 0, 0, 1e-18) &&
      same_candidate(moving.split.candidate, retry.split.candidate),
      "TAB rate response, dispersion energy and retry identity match independent values");
  Vector3 total_p{};
  for (std::uint32_t i=0; i<moving.split.candidate.child_count; ++i) {
    const auto& child=moving.split.candidate.children[i];
    for (std::size_t j=0;j<3;++j)
      total_p[j]+=child.droplet_mass_kg*child.multiplicity*child.velocity_m_per_s[j];
  }
  for (std::size_t j=0;j<3;++j)
    passed &= expect(near(total_p[j], tab.parent.droplet_mass_kg *
        tab.parent.multiplicity*tab.parent.velocity_m_per_s[j],1e-13,1e-22),
        "paired TAB child dispersion conserves momentum");
  tab.tab_trigger.candidate.deformation = 0.9;
  passed &= expect(!generate_tab_representative_children(tab).succeeded(),
      "TAB representative closure rejects non-unit trigger without child publication");
  tab.tab_trigger.candidate.deformation = 1.0;
  tab.child_parcel_count=3;
  passed &= expect(!generate_tab_representative_children(tab).split.candidate.available,
      "unpaired child count rejects before publication");
  passed &= test_conservative_uniform_split();
  passed &= test_breakup_rng_domain_and_retry_determinism();
  passed &= test_energy_disposition_and_multiplicity_scaling();
  passed &= test_failures_publish_no_partial_children();
  if (!passed) return 1;
  std::cout << "models_spray_breakup_test: PASS\n";
  return 0;
}
