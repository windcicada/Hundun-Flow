// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include "core_arena_detail.hpp"

#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {

namespace {

constexpr std::uint32_t kRevisionInvalidField = 1U;
constexpr std::uint32_t kRevisionExhausted = 2U;

}  // namespace

Status RevisionSet::create(const FieldSchema& schema, RevisionSet& out) {
  try {
    RevisionSet candidate;
    candidate.tokens_.assign(schema.size(), RevisionToken{1U});
    candidate.identity_ = detail::issue_identity();
    if (candidate.identity_ == 0U) {
      return {StatusCode::invalid_plan, kRevisionExhausted};
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kRevisionInvalidField};
  }
}

RevisionToken RevisionSet::token(FieldId field) const noexcept {
  const std::size_t index = static_cast<std::size_t>(field);
  return index < tokens_.size() ? tokens_[index] : RevisionToken{0U};
}

Status RevisionSet::revise(FieldId field) noexcept {
  const std::size_t index = static_cast<std::size_t>(field);
  if (index >= tokens_.size()) {
    return {StatusCode::invalid_plan, kRevisionInvalidField};
  }
  if (tokens_[index] == std::numeric_limits<RevisionToken>::max()) {
    return {StatusCode::invalid_plan, kRevisionExhausted};
  }
  ++tokens_[index];
  return {};
}

}  // namespace hundun::v04
