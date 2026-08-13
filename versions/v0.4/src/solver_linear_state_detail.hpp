// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_linear.hpp"

#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

constexpr std::uint64_t kLinearFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kLinearFnvPrime = UINT64_C(1099511628211);

inline std::uint64_t linear_hash_mix(std::uint64_t hash,
                                     std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kLinearFnvPrime;
  return hash;
}

inline PlanFingerprint finish_linear_hash(std::uint64_t hash) noexcept {
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

bool same_symbolic_spec(SymbolicSpec left, SymbolicSpec right) noexcept;
bool same_coefficients(CoefficientRevisions left,
                       CoefficientRevisions right) noexcept;
bool same_hierarchy_policy(HierarchyPolicyIdentity left,
                           HierarchyPolicyIdentity right) noexcept;
bool same_workspace_requirements(const LinearWorkspaceRequirements& left,
                                 const LinearWorkspaceRequirements& right)
    noexcept;
bool valid_workspace_requirements(
    const LinearWorkspaceRequirements& requirements) noexcept;
bool same_field_binding(FieldView left, FieldView right) noexcept;

PlanFingerprint symbolic_fingerprint(SymbolicSpec spec) noexcept;
PlanFingerprint numeric_fingerprint(PlanFingerprint symbolic,
                                    CoefficientRevisions revisions,
                                    PlanFingerprint content) noexcept;
PlanFingerprint hierarchy_fingerprint(
    PlanFingerprint symbolic, PlanFingerprint numeric,
    HierarchyPolicyIdentity policy, PlanFingerprint content) noexcept;
PlanFingerprint workspace_fingerprint(
    const LinearWorkspaceRequirements& requirements, FieldView vector_bundle,
    FieldView scalar_buffer) noexcept;

bool checked_linear_add(std::size_t left, std::size_t right,
                        std::size_t& out) noexcept;
bool checked_linear_multiply(std::size_t left, std::size_t right,
                             std::size_t& out) noexcept;
bool checked_counter_increment(std::uint64_t current,
                               std::uint64_t& out) noexcept;

}  // namespace hundun::v04::detail
