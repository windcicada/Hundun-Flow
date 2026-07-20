// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/preconditioners.hpp"

#include "hundun/runtime/error.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <utility>

namespace hundun::linear {
namespace {

struct CacheCandidate final {
  const LinearOperator* linear_operator{};
  const execution::ExecutionContext* context{};
  std::uint64_t revision{};
  VectorLayout domain;
  VectorLayout range;
};

void validate_reference_context(
    const execution::ExecutionContext& context) {
  if (context.backend_identity() == 0U) {
    throw runtime::Error(
        "reference preconditioner context identity must be nonzero");
  }
  if (context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access)) {
    throw runtime::Error(
        "reference preconditioner requires a host-accessible context");
  }
  if (!context.supports(
          execution::ExecutionCapability::buffer_allocation)) {
    throw runtime::Error(
        "reference preconditioner context lacks buffer capability");
  }
  if (!context.supports(execution::ExecutionCapability::transfer)) {
    throw runtime::Error(
        "reference preconditioner context lacks transfer capability");
  }
}

CacheCandidate make_candidate(
    const LinearOperator& linear_operator,
    const execution::ExecutionContext& expected_context,
    std::uint64_t revision, const LinearOperator* cached_operator,
    std::uint64_t cached_revision, bool cache_valid) {
  if (&linear_operator.context() != &expected_context) {
    throw runtime::Error(
        "preconditioner operator context object does not match");
  }
  if (linear_operator.revision() != revision) {
    throw runtime::Error(
        "preconditioner update revision does not match the operator");
  }
  if (cache_valid && cached_operator == &linear_operator &&
      revision < cached_revision) {
    throw runtime::Error(
        "preconditioner operator revision cannot decrease");
  }

  CacheCandidate candidate;
  candidate.linear_operator = &linear_operator;
  candidate.context = &expected_context;
  candidate.revision = revision;
  candidate.domain = linear_operator.domain_layout();
  candidate.range = linear_operator.range_layout();
  if (candidate.domain != candidate.range) {
    throw runtime::Error(
        "preconditioner requires an exactly square operator layout");
  }
  return candidate;
}

bool cache_identity_matches(
    const CacheCandidate& candidate, const LinearOperator* cached_operator,
    const execution::ExecutionContext* cached_context,
    std::uint64_t cached_revision, const VectorLayout& cached_domain,
    const VectorLayout& cached_range, bool cache_valid) noexcept {
  return cache_valid && candidate.linear_operator == cached_operator &&
         candidate.context == cached_context &&
         candidate.revision == cached_revision &&
         candidate.domain == cached_domain &&
         candidate.range == cached_range;
}

void validate_live_cache(
    const LinearOperator* linear_operator,
    const execution::ExecutionContext* context, std::uint64_t revision,
    const VectorLayout& domain, const VectorLayout& range, bool cache_valid) {
  if (!cache_valid || linear_operator == nullptr || context == nullptr) {
    throw runtime::Error(
        "preconditioner apply requires a successful update");
  }
  if (&linear_operator->context() != context) {
    throw runtime::Error(
        "preconditioner cached operator context has changed");
  }
  if (linear_operator->revision() != revision) {
    throw runtime::Error(
        "preconditioner cached operator revision has changed");
  }
  if (linear_operator->domain_layout() != domain ||
      linear_operator->range_layout() != range) {
    throw runtime::Error(
        "preconditioner cached operator layout has changed");
  }
}

const double* validate_input(
    execution::VectorView<const double> input,
    const execution::ExecutionContext& context, std::size_t expected_size) {
  if (input.scalar_format() != execution::ScalarFormat::float64) {
    throw runtime::Error("preconditioner input must be float64");
  }
  if (input.size() != expected_size) {
    throw runtime::Error(
        "preconditioner input size does not match owned layout");
  }
  if (input.backend_identity() != context.backend_identity()) {
    throw runtime::Error(
        "preconditioner input backend identity does not match");
  }
  if (input.space() != context.space() ||
      input.space() != execution::ExecutionSpace::host) {
    throw runtime::Error(
        "preconditioner input execution space does not match host context");
  }
  return input.data();
}

double* validate_output(
    execution::VectorView<double> output,
    const execution::ExecutionContext& context, std::size_t expected_size) {
  if (!output.writable()) {
    throw runtime::Error("preconditioner output must be writable");
  }
  if (output.scalar_format() != execution::ScalarFormat::float64) {
    throw runtime::Error("preconditioner output must be float64");
  }
  if (output.size() != expected_size) {
    throw runtime::Error(
        "preconditioner output size does not match owned layout");
  }
  if (output.backend_identity() != context.backend_identity()) {
    throw runtime::Error(
        "preconditioner output backend identity does not match");
  }
  if (output.space() != context.space() ||
      output.space() != execution::ExecutionSpace::host) {
    throw runtime::Error(
        "preconditioner output execution space does not match host context");
  }
  return output.data();
}

bool exact_same_view(execution::VectorView<const double> input,
                     execution::VectorView<double> output) noexcept {
  return input.allocation_identity() == output.allocation_identity() &&
         input.offset_bytes() == output.offset_bytes() &&
         input.size() == output.size() && input.stride() == output.stride();
}

void reject_partial_alias(execution::VectorView<const double> input,
                          execution::VectorView<double> output) {
  if (input.allocation_identity() == output.allocation_identity() &&
      !exact_same_view(input, output)) {
    throw runtime::Error(
        "preconditioner output aliases input without exact view identity");
  }
}

std::size_t checked_vector_bytes(std::size_t count) {
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw runtime::Error("preconditioner vector byte size overflow");
  }
  return count * sizeof(double);
}

}  // namespace

struct IdentityPreconditioner::State final {
  explicit State(execution::ExecutionContext& context_value) noexcept
      : context(&context_value) {}

  execution::ExecutionContext* context;
  const LinearOperator* linear_operator{};
  std::uint64_t revision{};
  VectorLayout domain;
  VectorLayout range;
  bool cache_valid{false};
};

struct JacobiPreconditioner::State final {
  explicit State(execution::ExecutionContext& context_value) noexcept
      : context(&context_value) {}

  execution::ExecutionContext* context;
  const LinearOperator* linear_operator{};
  std::uint64_t revision{};
  VectorLayout domain;
  VectorLayout range;
  std::optional<execution::Buffer> inverse_diagonal;
  bool cache_valid{false};
};

IdentityPreconditioner::IdentityPreconditioner(
    execution::ExecutionContext& context) {
  validate_reference_context(context);
  try {
    state_ = std::make_unique<State>(context);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("identity preconditioner state allocation failed");
  }
}

IdentityPreconditioner::~IdentityPreconditioner() noexcept = default;

void IdentityPreconditioner::update(
    const LinearOperator& linear_operator, std::uint64_t revision) {
  const auto candidate = make_candidate(
      linear_operator, *state_->context, revision, state_->linear_operator,
      state_->revision, state_->cache_valid);
  state_->linear_operator = candidate.linear_operator;
  state_->revision = candidate.revision;
  state_->domain = candidate.domain;
  state_->range = candidate.range;
  state_->cache_valid = true;
}

execution::ExecutionEvent IdentityPreconditioner::apply(
    execution::VectorView<const double> residual,
    execution::VectorView<double> correction) const {
  validate_live_cache(state_->linear_operator, state_->context,
                      state_->revision, state_->domain, state_->range,
                      state_->cache_valid);
  const std::size_t count = state_->domain.owned_count();
  const double* input = validate_input(residual, *state_->context, count);
  double* output = validate_output(correction, *state_->context, count);
  reject_partial_alias(residual, correction);
  for (std::size_t index = 0; index < count; ++index) {
    if (!std::isfinite(input[index * residual.stride()])) {
      throw runtime::Error("identity preconditioner input must be finite");
    }
  }
  auto event = execution::ExecutionEvent::completed();
  if (!exact_same_view(residual, correction)) {
    for (std::size_t index = 0; index < count; ++index) {
      output[index * correction.stride()] =
          input[index * residual.stride()];
    }
  }
  return event;
}

JacobiPreconditioner::JacobiPreconditioner(
    execution::ExecutionContext& context) {
  validate_reference_context(context);
  try {
    state_ = std::make_unique<State>(context);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("Jacobi preconditioner state allocation failed");
  }
}

JacobiPreconditioner::~JacobiPreconditioner() noexcept = default;

void JacobiPreconditioner::update(
    const LinearOperator& linear_operator, std::uint64_t revision) {
  const auto candidate = make_candidate(
      linear_operator, *state_->context, revision, state_->linear_operator,
      state_->revision, state_->cache_valid);
  if (!linear_operator.has_diagonal()) {
    throw runtime::Error(
        "Jacobi preconditioner requires diagonal capability");
  }
  if (cache_identity_matches(
          candidate, state_->linear_operator, state_->context,
          state_->revision, state_->domain, state_->range,
          state_->cache_valid)) {
    return;
  }

  const std::size_t count = candidate.range.owned_count();
  execution::Buffer staging(*state_->context, checked_vector_bytes(count));
  auto staging_view = staging.view(0U, count);
  auto diagonal_event = linear_operator.diagonal(staging_view);
  diagonal_event.wait();

  for (std::size_t index = 0; index < count; ++index) {
    const double diagonal = staging_view[index];
    if (!std::isfinite(diagonal)) {
      throw runtime::Error("Jacobi diagonal value must be finite");
    }
    if (diagonal == 0.0) {
      throw runtime::Error("Jacobi diagonal value must be nonzero");
    }
    if (!std::isfinite(1.0 / diagonal)) {
      throw runtime::Error("Jacobi diagonal reciprocal must be finite");
    }
  }
  for (std::size_t index = 0; index < count; ++index) {
    staging_view[index] = 1.0 / staging_view[index];
  }

  std::unique_ptr<State> replacement;
  try {
    replacement = std::make_unique<State>(*state_->context);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("Jacobi preconditioner state allocation failed");
  }
  replacement->linear_operator = candidate.linear_operator;
  replacement->revision = candidate.revision;
  replacement->domain = candidate.domain;
  replacement->range = candidate.range;
  replacement->inverse_diagonal.emplace(std::move(staging));
  replacement->cache_valid = true;
  state_.swap(replacement);
}

execution::ExecutionEvent JacobiPreconditioner::apply(
    execution::VectorView<const double> residual,
    execution::VectorView<double> correction) const {
  validate_live_cache(state_->linear_operator, state_->context,
                      state_->revision, state_->domain, state_->range,
                      state_->cache_valid);
  if (!state_->inverse_diagonal.has_value()) {
    throw runtime::Error(
        "Jacobi preconditioner apply requires a cached diagonal");
  }
  const std::size_t count = state_->domain.owned_count();
  const double* input = validate_input(residual, *state_->context, count);
  double* output = validate_output(correction, *state_->context, count);
  reject_partial_alias(residual, correction);
  const auto inverse =
      static_cast<const execution::Buffer&>(*state_->inverse_diagonal)
          .view(0U, count);
  const double* inverse_values = inverse.data();
  for (std::size_t index = 0; index < count; ++index) {
    const double input_value = input[index * residual.stride()];
    const double result = input_value * inverse_values[index];
    if (!std::isfinite(input_value)) {
      throw runtime::Error("Jacobi preconditioner input must be finite");
    }
    if (!std::isfinite(result)) {
      throw runtime::Error("Jacobi preconditioner result must be finite");
    }
  }
  auto event = execution::ExecutionEvent::completed();
  for (std::size_t index = 0; index < count; ++index) {
    output[index * correction.stride()] =
        input[index * residual.stride()] * inverse_values[index];
  }
  return event;
}

}  // namespace hundun::linear
