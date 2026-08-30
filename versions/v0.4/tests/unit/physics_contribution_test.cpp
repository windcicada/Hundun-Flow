// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_physics.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

UnitDimension source_units(std::int8_t time_exponent = -1) {
  UnitDimension result;
  result.si_exponents = {1, 0, time_exponent, 0, 0, 0, 0};
  return result;
}

ContributionSpec contribution(FieldId conserved, StageId stage,
                              Span<const FieldId> reads,
                              FieldId explicit_source,
                              UnitDimension units = source_units()) {
  ContributionSpec result;
  result.conserved_quantity = conserved;
  result.units = units;
  result.stage = stage;
  result.reads = reads;
  result.explicit_source = explicit_source;
  result.capability = ContributionCapability::inert_source;
  return result;
}

bool same_contribution(const CompiledContribution& left,
                       const CompiledContribution& right) {
  return left.conserved_quantity == right.conserved_quantity &&
         left.units == right.units && left.stage == right.stage &&
         left.read_begin == right.read_begin &&
         left.read_count == right.read_count &&
         left.explicit_source == right.explicit_source &&
         left.implicit_diagonal == right.implicit_diagonal &&
         left.supplies_implicit_diagonal ==
             right.supplies_implicit_diagonal &&
         left.registration_ordinal == right.registration_ordinal;
}

bool same_frozen_registry(const ContributionRegistry& left,
                          const ContributionRegistry& right) {
  if (!left.frozen() || !right.frozen() ||
      left.fingerprint() != right.fingerprint() ||
      left.contributions().size != right.contributions().size ||
      left.reads().size != right.reads().size) {
    return false;
  }
  for (std::size_t index = 0U; index < left.contributions().size; ++index) {
    if (!same_contribution(left.contributions().data[index],
                           right.contributions().data[index])) {
      return false;
    }
  }
  for (std::size_t index = 0U; index < left.reads().size; ++index) {
    if (left.reads().data[index] != right.reads().data[index]) {
      return false;
    }
  }
  return true;
}

bool test_effective_viscosity_authority() {
  bool passed = true;
  ContributionRegistry registry;
  const std::array<FieldId, 2U> declared{0U, 1U};
  passed &= expect(static_cast<bool>(
                       registry.configure({declared.data(), declared.size()})),
                   "mu_eff registry configures");
  EffectiveViscosityAuthority authority;
  passed &= expect(!authority.claimed(), "mu_eff starts unclaimed");
  passed &= expect(static_cast<bool>(authority.claim(0U, 7U, registry)),
                   "the base mu_eff writer can claim field zero");
  passed &= expect(authority.claimed() && authority.output() == 0U &&
                       authority.stage() == 7U,
                   "a successful mu_eff claim publishes its authority");
  passed &= expect(!authority.claim(0U, 7U, registry),
                   "even an identical second mu_eff writer is rejected");
  passed &= expect(!authority.claim(1U, 8U, registry),
                   "a competing mu_eff writer is rejected");
  passed &= expect(authority.output() == 0U && authority.stage() == 7U,
                   "a rejected mu_eff claim preserves the first authority");

  EffectiveViscosityAuthority competing_instance;
  passed &= expect(!competing_instance.claim(1U, 8U, registry),
                   "a second authority instance cannot claim the registry");
  passed &= expect(registry.plan() == nullptr,
                   "mutable registration does not expose a runtime plan");
  passed &= expect(static_cast<bool>(registry.freeze()) &&
                       registry.fingerprint() != 0U,
                   "mu_eff authority participates in production freeze");
  const ContributionPlan* plan = registry.plan();
  passed &= expect(plan != nullptr && plan->valid() &&
                       plan->fingerprint() == registry.fingerprint() &&
                       plan->contributions().size ==
                           registry.contributions().size &&
                       plan->reads().size == registry.reads().size,
                   "freeze publishes one immutable ContributionPlan seam");
  return passed;
}

bool test_effective_viscosity_conflicts_with_contribution_outputs() {
  bool passed = true;
  const std::array<FieldId, 6U> declared{0U, 1U, 2U, 3U, 4U, 5U};
  const std::array<FieldId, 1U> reads{2U};

  ContributionRegistry authority_first;
  EffectiveViscosityAuthority authority;
  passed &= expect(
      static_cast<bool>(authority_first.configure(
          {declared.data(), declared.size()})) &&
          static_cast<bool>(authority.claim(4U, 7U, authority_first)),
      "mu_eff claims its output before contribution registration");
  ContributionSpec explicit_conflict =
      contribution(0U, 3U, {reads.data(), reads.size()}, 4U);
  passed &= expect(!authority_first.register_contribution(explicit_conflict),
                   "a contribution source cannot steal the mu_eff output");
  ContributionSpec diagonal_conflict =
      contribution(0U, 3U, {reads.data(), reads.size()}, 1U);
  diagonal_conflict.supplies_implicit_diagonal = true;
  diagonal_conflict.implicit_diagonal = 4U;
  passed &= expect(!authority_first.register_contribution(diagonal_conflict),
                   "a contribution diagonal cannot steal the mu_eff output");
  passed &= expect(authority_first.contributions().size == 0U,
                   "mu_eff output conflicts publish no contribution state");

  ContributionRegistry explicit_first;
  passed &= expect(static_cast<bool>(explicit_first.configure(
                       {declared.data(), declared.size()})) &&
                       static_cast<bool>(explicit_first.register_contribution(
                           contribution(0U, 3U,
                                        {reads.data(), reads.size()}, 4U))),
                   "a contribution claims an explicit output first");
  EffectiveViscosityAuthority explicit_competitor;
  passed &= expect(!explicit_competitor.claim(4U, 7U, explicit_first) &&
                       !explicit_competitor.claimed(),
                   "mu_eff cannot steal a contribution source output");

  ContributionRegistry diagonal_first;
  passed &= expect(static_cast<bool>(diagonal_first.configure(
                       {declared.data(), declared.size()})),
                   "diagonal-first registry configures");
  ContributionSpec diagonal_owner =
      contribution(0U, 3U, {reads.data(), reads.size()}, 1U);
  diagonal_owner.supplies_implicit_diagonal = true;
  diagonal_owner.implicit_diagonal = 4U;
  passed &= expect(
      static_cast<bool>(diagonal_first.register_contribution(diagonal_owner)),
      "a contribution claims an implicit output first");
  EffectiveViscosityAuthority diagonal_competitor;
  passed &= expect(!diagonal_competitor.claim(4U, 7U, diagonal_first) &&
                       !diagonal_competitor.claimed(),
                   "mu_eff cannot steal a contribution diagonal output");
  return passed;
}

bool test_configuration_is_fail_closed_and_transactional() {
  bool passed = true;
  ContributionRegistry registry;
  const std::array<FieldId, 1U> read{2U};
  const ContributionSpec spec =
      contribution(0U, 3U, {read.data(), read.size()}, 1U);
  passed &= expect(!registry.register_contribution(spec),
                   "registration before the declared-field catalog fails");
  passed &= expect(!registry.configure({nullptr, 0U}),
                   "an empty declared-field catalog is rejected");

  const std::array<FieldId, 4U> duplicate{0U, 2U, 1U, 2U};
  passed &= expect(!registry.configure({duplicate.data(), duplicate.size()}),
                   "duplicate declared fields are rejected");
  const std::array<FieldId, 3U> declared{2U, 0U, 1U};
  passed &= expect(static_cast<bool>(
                       registry.configure({declared.data(), declared.size()})),
                   "a failed configure leaves the registry configurable");
  passed &= expect(!registry.configure({declared.data(), declared.size()}),
                   "the declared-field catalog has one authority");
  passed &= expect(static_cast<bool>(registry.register_contribution(spec)),
                   "field zero remains a valid declared conserved quantity");
  return passed;
}

bool test_registration_validation() {
  bool passed = true;
  const std::array<FieldId, 8U> declared{0U, 1U, 2U, 3U,
                                        4U, 5U, 6U, 7U};
  ContributionRegistry registry;
  passed &= expect(static_cast<bool>(
                       registry.configure({declared.data(), declared.size()})),
                   "validation registry configures");
  const std::array<FieldId, 2U> valid_reads{2U, 3U};
  const auto valid = contribution(0U, 4U,
                                  {valid_reads.data(), valid_reads.size()},
                                  1U);

  ContributionSpec candidate = valid;
  candidate.stage = 0U;
  passed &= expect(!registry.register_contribution(candidate),
                   "stage zero cannot enter the frozen execution graph");
  candidate = valid;
  candidate.conserved_quantity = 20U;
  passed &= expect(!registry.register_contribution(candidate),
                   "an undeclared conserved field is rejected");
  candidate = valid;
  candidate.explicit_source = 20U;
  passed &= expect(!registry.register_contribution(candidate),
                   "an undeclared explicit source field is rejected");
  const std::array<FieldId, 2U> undeclared_reads{2U, 20U};
  candidate = contribution(0U, 4U,
                           {undeclared_reads.data(), undeclared_reads.size()},
                           1U);
  passed &= expect(!registry.register_contribution(candidate),
                   "an undeclared read is rejected");
  candidate = valid;
  candidate.supplies_implicit_diagonal = true;
  candidate.implicit_diagonal = 20U;
  passed &= expect(!registry.register_contribution(candidate),
                   "an undeclared implicit diagonal is rejected");
  candidate.supplies_implicit_diagonal = false;
  passed &= expect(static_cast<bool>(registry.register_contribution(candidate)),
                   "an inactive optional diagonal requires no field authority");

  candidate = valid;
  candidate.units = UnitDimension{};
  passed &= expect(!registry.register_contribution(candidate),
                   "a contribution without physical dimensions is rejected");
  candidate = valid;
  candidate.reads = {nullptr, 0U};
  passed &= expect(!registry.register_contribution(candidate),
                   "an absent read set is rejected");
  FieldId oversized_probe = 2U;
  candidate = valid;
  candidate.reads = {
      &oversized_probe,
      static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()) +
          1U};
  passed &= expect(!registry.register_contribution(candidate),
                   "an unrepresentable read set is rejected before access");
  const std::array<FieldId, 2U> duplicate_reads{2U, 2U};
  candidate = contribution(0U, 4U,
                           {duplicate_reads.data(), duplicate_reads.size()},
                           1U);
  passed &= expect(!registry.register_contribution(candidate),
                   "duplicate declared reads are rejected");

  candidate = valid;
  candidate.capability = ContributionCapability::chemistry;
  passed &= expect(!registry.register_contribution(candidate),
                   "chemistry contributions are outside v0.4 scope");
  candidate.capability = ContributionCapability::reacting;
  passed &= expect(!registry.register_contribution(candidate),
                   "reacting contributions are outside v0.4 scope");

  candidate = valid;
  candidate.explicit_source = candidate.conserved_quantity;
  passed &= expect(!registry.register_contribution(candidate),
                   "the conserved state cannot double as source output");
  const std::array<FieldId, 2U> source_read{1U, 2U};
  candidate = contribution(0U, 4U, {source_read.data(), source_read.size()},
                           1U);
  passed &= expect(!registry.register_contribution(candidate),
                   "an explicit source cannot also be a declared read");

  candidate = valid;
  candidate.supplies_implicit_diagonal = true;
  candidate.implicit_diagonal = candidate.conserved_quantity;
  passed &= expect(!registry.register_contribution(candidate),
                   "the conserved state cannot double as implicit output");
  candidate.implicit_diagonal = candidate.explicit_source;
  passed &= expect(!registry.register_contribution(candidate),
                   "explicit and implicit outputs must be distinct");
  candidate.implicit_diagonal = 2U;
  passed &= expect(!registry.register_contribution(candidate),
                   "an implicit output cannot also be a declared read");

  passed &= expect(registry.contributions().size == 1U &&
                       registry.reads().size == valid_reads.size(),
                   "rejected registrations publish no partial state");
  passed &= expect(static_cast<bool>(registry.freeze()),
                   "the valid survivor freezes");
  passed &= expect(registry.contributions().size == 1U &&
                       !registry.contributions().data[0].
                           supplies_implicit_diagonal &&
                       registry.contributions().data[0].implicit_diagonal ==
                           FieldId{},
                   "the inactive optional implicit field is canonicalized");
  return passed;
}

bool test_single_writers() {
  bool passed = true;
  const std::array<FieldId, 12U> declared{0U, 1U, 2U, 3U, 4U, 5U,
                                         6U, 7U, 8U, 9U, 10U, 11U};
  const std::array<FieldId, 1U> first_reads{2U};
  ContributionRegistry registry;
  passed &= expect(static_cast<bool>(
                       registry.configure({declared.data(), declared.size()})),
                   "writer registry configures");
  ContributionSpec first =
      contribution(0U, 2U, {first_reads.data(), first_reads.size()}, 1U);
  first.supplies_implicit_diagonal = true;
  first.implicit_diagonal = 3U;
  passed &= expect(static_cast<bool>(registry.register_contribution(first)),
                   "the first explicit and implicit writers register");

  const auto rejected_writer = [&](FieldId explicit_source,
                                   bool supplies_implicit,
                                   FieldId implicit) {
    const std::array<FieldId, 1U> reads{6U};
    ContributionSpec spec =
        contribution(4U, 5U, {reads.data(), reads.size()}, explicit_source);
    spec.supplies_implicit_diagonal = supplies_implicit;
    spec.implicit_diagonal = implicit;
    return !registry.register_contribution(spec);
  };
  passed &= expect(rejected_writer(1U, false, 0U),
                   "two explicit writers for one field are rejected");
  passed &= expect(rejected_writer(3U, false, 0U),
                   "an explicit writer cannot steal an implicit output");
  passed &= expect(rejected_writer(7U, true, 1U),
                   "an implicit writer cannot steal an explicit output");
  passed &= expect(rejected_writer(7U, true, 3U),
                   "two implicit writers for one field are rejected");
  passed &= expect(registry.contributions().size == 1U,
                   "writer conflicts consume no registrations");

  const std::array<FieldId, 1U> second_reads{6U};
  ContributionSpec second =
      contribution(0U, 5U, {second_reads.data(), second_reads.size()}, 7U);
  second.supplies_implicit_diagonal = true;
  second.implicit_diagonal = 8U;
  passed &= expect(static_cast<bool>(registry.register_contribution(second)),
                   "multiple contributions may target one equation with "
                   "distinct source authorities");
  passed &= expect(registry.contributions().data[1].registration_ordinal == 1U,
                   "failed writer claims do not consume stable ordinals");
  return passed;
}

bool populate_ordered_registry(ContributionRegistry& registry,
                               bool reverse_catalog = false,
                               UnitDimension first_units = source_units()) {
  const std::array<FieldId, 12U> ascending{0U, 1U, 2U, 3U, 4U, 5U,
                                          6U, 7U, 8U, 9U, 10U, 11U};
  const std::array<FieldId, 12U> descending{11U, 10U, 9U, 8U, 7U, 6U,
                                           5U, 4U, 3U, 2U, 1U, 0U};
  const auto& catalog = reverse_catalog ? descending : ascending;
  if (!registry.configure({catalog.data(), catalog.size()})) {
    return false;
  }

  const std::array<FieldId, 2U> reads_a{6U, 2U};
  ContributionSpec a =
      contribution(4U, 7U, {reads_a.data(), reads_a.size()}, 9U, first_units);
  a.supplies_implicit_diagonal = true;
  a.implicit_diagonal = 10U;
  if (!registry.register_contribution(a)) {
    return false;
  }

  const std::array<FieldId, 1U> reads_b{5U};
  if (!registry.register_contribution(
          contribution(4U, 2U, {reads_b.data(), reads_b.size()}, 8U))) {
    return false;
  }

  const std::array<FieldId, 1U> reads_c{7U};
  if (!registry.register_contribution(
          contribution(3U, 7U, {reads_c.data(), reads_c.size()}, 11U))) {
    return false;
  }
  return static_cast<bool>(registry.freeze());
}

bool test_deterministic_freeze_and_fingerprint() {
  bool passed = true;
  ContributionRegistry first;
  ContributionRegistry replay;
  passed &= expect(populate_ordered_registry(first),
                   "the reference contribution registry freezes");
  passed &= expect(populate_ordered_registry(replay, true),
                   "the replay contribution registry freezes");
  passed &= expect(same_frozen_registry(first, replay),
                   "declared-field input order cannot affect the frozen plan");
  passed &= expect(first.fingerprint() != 0U,
                   "a frozen contribution registry has an identity");

  const auto compiled = first.contributions();
  passed &= expect(compiled.size == 3U,
                   "every valid contribution is present after freeze");
  if (compiled.size == 3U) {
    passed &= expect(compiled.data[0].stage == 2U &&
                         compiled.data[0].conserved_quantity == 4U &&
                         compiled.data[0].registration_ordinal == 1U,
                     "stage is the first deterministic ordering key");
    passed &= expect(compiled.data[1].stage == 7U &&
                         compiled.data[1].conserved_quantity == 3U &&
                         compiled.data[1].registration_ordinal == 2U,
                     "conserved quantity is the second ordering key");
    passed &= expect(compiled.data[2].stage == 7U &&
                         compiled.data[2].conserved_quantity == 4U &&
                         compiled.data[2].registration_ordinal == 0U,
                     "the stable registration ordinal is retained");
    passed &= expect(compiled.data[2].read_count == 2U &&
                         first.reads().data[compiled.data[2].read_begin] == 2U &&
                         first.reads().data[compiled.data[2].read_begin + 1U] ==
                             6U,
                     "each compact read span is sorted and remapped with its "
                     "contribution");
  }

  const PlanFingerprint frozen_fingerprint = first.fingerprint();
  const std::size_t frozen_contributions = first.contributions().size;
  const std::size_t frozen_reads = first.reads().size;
  const std::array<FieldId, 1U> late_reads{0U};
  passed &= expect(!first.register_contribution(
                       contribution(1U, 9U,
                                    {late_reads.data(), late_reads.size()},
                                    2U)) &&
                       !first.freeze(),
                   "registration and a second freeze are rejected after "
                   "production freeze");
  const std::array<FieldId, 1U> late_catalog{0U};
  passed &= expect(!first.configure({late_catalog.data(), late_catalog.size()}),
                   "the declared-field catalog is immutable after freeze");
  passed &= expect(first.fingerprint() == frozen_fingerprint &&
                       first.contributions().size == frozen_contributions &&
                       first.reads().size == frozen_reads,
                   "post-freeze rejection preserves the frozen plan");

  ContributionRegistry changed_units;
  passed &= expect(populate_ordered_registry(
                       changed_units, false, source_units(-2)),
                   "the units mutation remains structurally valid");
  passed &= expect(changed_units.fingerprint() != replay.fingerprint(),
                   "unit dimensions participate in the plan fingerprint");
  return passed;
}

}  // namespace

int main() {
  bool passed = true;
  passed &= test_effective_viscosity_authority();
  passed &= test_effective_viscosity_conflicts_with_contribution_outputs();
  passed &= test_configuration_is_fail_closed_and_transactional();
  passed &= test_registration_validation();
  passed &= test_single_writers();
  passed &= test_deterministic_freeze_and_fingerprint();
  return passed ? 0 : 1;
}
