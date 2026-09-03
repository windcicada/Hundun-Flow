// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

bool has_dimensions(UnitDimension units) noexcept {
  return std::any_of(units.si_exponents.begin(), units.si_exponents.end(),
                     [](std::int8_t exponent) { return exponent != 0; });
}

bool contains(const std::vector<FieldId>& fields, FieldId field) noexcept {
  return std::binary_search(fields.begin(), fields.end(), field);
}

bool contribution_less(const CompiledContribution& left,
                       const CompiledContribution& right) noexcept {
  if (left.stage != right.stage) {
    return left.stage < right.stage;
  }
  if (left.conserved_quantity != right.conserved_quantity) {
    return left.conserved_quantity < right.conserved_quantity;
  }
  return left.registration_ordinal < right.registration_ordinal;
}

bool output_is_claimed(const std::vector<CompiledContribution>& contributions,
                       FieldId output) noexcept {
  return std::any_of(
      contributions.begin(), contributions.end(),
      [output](const CompiledContribution& contribution) {
        return contribution.explicit_source == output ||
               (contribution.supplies_implicit_diagonal &&
                contribution.implicit_diagonal == output);
      });
}

PlanFingerprint contribution_fingerprint(
    const std::vector<FieldId>& declared_fields,
    const std::vector<CompiledContribution>& contributions,
    const std::vector<FieldId>& reads, bool mu_eff_claimed,
    FieldId mu_eff_output, StageId mu_eff_stage) noexcept {
  std::uint64_t hash = kFnvOffset;
  hash = hash_mix(hash, declared_fields.size());
  for (const FieldId field : declared_fields) {
    hash = hash_mix(hash, field);
  }
  hash = hash_mix(hash, contributions.size());
  for (const CompiledContribution& contribution : contributions) {
    hash = hash_mix(hash, contribution.conserved_quantity);
    for (const std::int8_t exponent : contribution.units.si_exponents) {
      hash = hash_mix(hash, static_cast<std::uint8_t>(exponent));
    }
    hash = hash_mix(hash, contribution.stage);
    hash = hash_mix(hash, contribution.read_count);
    const std::size_t begin = contribution.read_begin;
    const std::size_t end = begin + contribution.read_count;
    if (end > reads.size()) {
      return 0U;
    }
    for (std::size_t index = begin; index < end; ++index) {
      hash = hash_mix(hash, reads[index]);
    }
    hash = hash_mix(hash, contribution.explicit_source);
    hash = hash_mix(hash, contribution.supplies_implicit_diagonal ? 1U : 0U);
    hash = hash_mix(hash, contribution.supplies_implicit_diagonal
                              ? contribution.implicit_diagonal
                              : 0U);
    hash = hash_mix(hash, contribution.registration_ordinal);
  }
  hash = hash_mix(hash, mu_eff_claimed ? 1U : 0U);
  if (mu_eff_claimed) {
    hash = hash_mix(hash, mu_eff_output);
    hash = hash_mix(hash, mu_eff_stage);
  }
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

}  // namespace

Status EffectiveViscosityAuthority::claim(
    FieldId output, StageId stage,
    ContributionRegistry& registry) noexcept {
  if (claimed_ || registry.frozen_ || registry.declared_fields_.empty() ||
      !contains(registry.declared_fields_, output) ||
      registry.effective_viscosity_claimed_ ||
      output_is_claimed(registry.contributions_, output)) {
    return {StatusCode::invalid_plan, 2U};
  }
  registry.effective_viscosity_claimed_ = true;
  registry.effective_viscosity_output_ = output;
  registry.effective_viscosity_stage_ = stage;
  output_ = output;
  stage_ = stage;
  claimed_ = true;
  return {};
}

Status ContributionRegistry::configure(
    Span<const FieldId> declared_fields) noexcept {
  if (frozen_ || !declared_fields_.empty() || !contributions_.empty() ||
      !reads_.empty() || next_ordinal_ != 0U || plan_.valid()) {
    return {StatusCode::invalid_plan, 1U};
  }
  if (declared_fields.data == nullptr || declared_fields.size == 0U) {
    return {StatusCode::invalid_plan, 2U};
  }
  if (declared_fields.size >
      static_cast<std::size_t>(std::numeric_limits<FieldId>::max()) + 1U) {
    return {StatusCode::invalid_plan, 3U};
  }
  try {
    std::vector<FieldId> candidate(declared_fields.data,
                                   declared_fields.data + declared_fields.size);
    std::sort(candidate.begin(), candidate.end());
    if (std::adjacent_find(candidate.begin(), candidate.end()) !=
        candidate.end()) {
      return {StatusCode::invalid_plan, 3U};
    }
    declared_fields_.swap(candidate);
  } catch (...) {
    return {StatusCode::allocation_failure, 0U};
  }
  return {};
}

Status ContributionRegistry::register_contribution(
    const ContributionSpec& spec) noexcept {
  if (frozen_ || declared_fields_.empty()) {
    return {StatusCode::invalid_plan, 4U};
  }
  if (spec.stage == 0U ||
      spec.capability != ContributionCapability::inert_source) {
    return {StatusCode::invalid_plan, 5U};
  }
  if (!has_dimensions(spec.units)) {
    return {StatusCode::invalid_plan, 6U};
  }
  if (spec.reads.data == nullptr || spec.reads.size == 0U ||
      spec.reads.size > std::numeric_limits<std::uint16_t>::max()) {
    return {StatusCode::invalid_plan, 8U};
  }
  if (!contains(declared_fields_, spec.conserved_quantity) ||
      !contains(declared_fields_, spec.explicit_source) ||
      (spec.supplies_implicit_diagonal &&
       !contains(declared_fields_, spec.implicit_diagonal))) {
    return {StatusCode::invalid_plan, 9U};
  }

  std::vector<FieldId> normalized_reads;
  try {
    normalized_reads.assign(spec.reads.data,
                            spec.reads.data + spec.reads.size);
    std::sort(normalized_reads.begin(), normalized_reads.end());
  } catch (...) {
    return {StatusCode::allocation_failure, 0U};
  }
  if (std::adjacent_find(normalized_reads.begin(), normalized_reads.end()) !=
      normalized_reads.end()) {
    return {StatusCode::invalid_plan, 10U};
  }
  if (std::any_of(normalized_reads.begin(), normalized_reads.end(),
                  [&](FieldId field) {
                    return !contains(declared_fields_, field);
                  })) {
    return {StatusCode::invalid_plan, 9U};
  }
  if (spec.explicit_source == spec.conserved_quantity ||
      contains(normalized_reads, spec.explicit_source) ||
      (spec.supplies_implicit_diagonal &&
       (spec.implicit_diagonal == spec.conserved_quantity ||
        spec.implicit_diagonal == spec.explicit_source ||
        contains(normalized_reads, spec.implicit_diagonal)))) {
    return {StatusCode::invalid_plan, 11U};
  }
  if (output_is_claimed(contributions_, spec.explicit_source) ||
      (spec.supplies_implicit_diagonal &&
       output_is_claimed(contributions_, spec.implicit_diagonal)) ||
      (effective_viscosity_claimed_ &&
       (spec.explicit_source == effective_viscosity_output_ ||
        (spec.supplies_implicit_diagonal &&
         spec.implicit_diagonal == effective_viscosity_output_)))) {
    return {StatusCode::invalid_plan, 12U};
  }
  if (next_ordinal_ == std::numeric_limits<std::uint32_t>::max() ||
      reads_.size() > std::numeric_limits<std::uint32_t>::max() -
                          normalized_reads.size()) {
    return {StatusCode::invalid_plan, 13U};
  }

  try {
    reads_.reserve(reads_.size() + normalized_reads.size());
    contributions_.reserve(contributions_.size() + 1U);
  } catch (...) {
    return {StatusCode::allocation_failure, 0U};
  }

  const auto read_begin = static_cast<std::uint32_t>(reads_.size());
  reads_.insert(reads_.end(), normalized_reads.begin(),
                normalized_reads.end());
  contributions_.push_back(CompiledContribution{
      spec.conserved_quantity,
      spec.units,
      spec.stage,
      read_begin,
      static_cast<std::uint16_t>(normalized_reads.size()),
      spec.explicit_source,
      spec.supplies_implicit_diagonal ? spec.implicit_diagonal : FieldId{},
      spec.supplies_implicit_diagonal,
      next_ordinal_});
  ++next_ordinal_;
  return {};
}

Status ContributionRegistry::freeze() noexcept {
  if (frozen_ || declared_fields_.empty()) {
    return {StatusCode::invalid_plan, 14U};
  }
  try {
    std::vector<CompiledContribution> ordered = contributions_;
    std::stable_sort(ordered.begin(), ordered.end(), contribution_less);

    std::vector<FieldId> ordered_reads;
    ordered_reads.reserve(reads_.size());
    for (CompiledContribution& contribution : ordered) {
      const std::size_t begin = contribution.read_begin;
      const std::size_t end = begin + contribution.read_count;
      if (end > reads_.size() ||
          ordered_reads.size() > std::numeric_limits<std::uint32_t>::max()) {
        return {StatusCode::invalid_plan, 15U};
      }
      contribution.read_begin =
          static_cast<std::uint32_t>(ordered_reads.size());
      ordered_reads.insert(ordered_reads.end(), reads_.begin() + begin,
                           reads_.begin() + end);
    }
    const PlanFingerprint candidate_fingerprint = contribution_fingerprint(
        declared_fields_, ordered, ordered_reads,
        effective_viscosity_claimed_, effective_viscosity_output_,
        effective_viscosity_stage_);
    if (candidate_fingerprint == 0U) {
      return {StatusCode::invalid_plan, 15U};
    }
    plan_.contributions_.swap(ordered);
    plan_.reads_.swap(ordered_reads);
    plan_.fingerprint_ = candidate_fingerprint;
    contributions_.clear();
    reads_.clear();
    frozen_ = true;
  } catch (...) {
    return {StatusCode::allocation_failure, 0U};
  }
  return {};
}

}  // namespace hundun::v04
