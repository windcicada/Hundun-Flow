// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/preconditioners.hpp"

#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
#include "preconditioners_test_access.hpp"

#include <atomic>
#endif

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

#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
std::atomic<bool> fail_next_cold_update_before_publication{false};
#endif

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
  std::optional<execution::Buffer> staging_inverse_diagonal;
  bool cache_valid{false};
};

namespace {

void fill_inverse_diagonal(const LinearOperator& linear_operator,
                           execution::Buffer& destination,
                           std::size_t count) {
  auto staging_view = destination.view(0U, count);
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
}

bool reusable_jacobi_storage(
    bool cache_valid,
    const std::optional<execution::Buffer>& inverse_diagonal,
    const std::optional<execution::Buffer>& staging_inverse_diagonal,
    const execution::ExecutionContext* context, const VectorLayout& domain,
    const VectorLayout& range, const CacheCandidate& candidate,
    std::size_t byte_size) noexcept {
  return cache_valid && inverse_diagonal.has_value() &&
         staging_inverse_diagonal.has_value() &&
         candidate.context == context && candidate.domain == domain &&
         candidate.range == range &&
         inverse_diagonal->byte_size() == byte_size &&
         staging_inverse_diagonal->byte_size() == byte_size;
}

}  // namespace

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
  const std::size_t byte_size = checked_vector_bytes(count);
  if (reusable_jacobi_storage(
          state_->cache_valid, state_->inverse_diagonal,
          state_->staging_inverse_diagonal, state_->context, state_->domain,
          state_->range, candidate, byte_size)) {
    fill_inverse_diagonal(linear_operator,
                          *state_->staging_inverse_diagonal, count);
    state_->inverse_diagonal.swap(state_->staging_inverse_diagonal);
    state_->linear_operator = candidate.linear_operator;
    state_->revision = candidate.revision;
    return;
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
  replacement->inverse_diagonal.emplace(*state_->context, byte_size);
  replacement->staging_inverse_diagonal.emplace(*state_->context, byte_size);
  fill_inverse_diagonal(linear_operator, *replacement->inverse_diagonal,
                        count);
#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
  if (fail_next_cold_update_before_publication.exchange(
          false, std::memory_order_relaxed)) {
    throw runtime::Error(
        "injected Jacobi cold update failure before publication");
  }
#endif
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

#ifdef HUNDUN_LINEAR_ENABLE_TEST_ACCESS
test::JacobiStorageSnapshot test::PreconditionerTestAccess::jacobi_storage(
    const JacobiPreconditioner& preconditioner) {
  test::JacobiStorageSnapshot snapshot;
  if (!preconditioner.state_) {
    return snapshot;
  }
  const auto& state = *preconditioner.state_;
  snapshot.revision = state.revision;
  snapshot.cache_valid = state.cache_valid;
  if (state.inverse_diagonal.has_value()) {
    snapshot.allocation_identities[0] =
        state.inverse_diagonal->allocation_identity();
    snapshot.byte_sizes[0] = state.inverse_diagonal->byte_size();
    const auto count = state.domain.owned_count();
    const auto inverse = static_cast<const execution::Buffer&>(
                             *state.inverse_diagonal)
                             .view(0U, count);
    snapshot.cached_inverse.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
      snapshot.cached_inverse.push_back(inverse[index]);
    }
  }
  if (state.staging_inverse_diagonal.has_value()) {
    snapshot.allocation_identities[1] =
        state.staging_inverse_diagonal->allocation_identity();
    snapshot.byte_sizes[1] =
        state.staging_inverse_diagonal->byte_size();
  }
  return snapshot;
}

void test::PreconditionerTestAccess::
    arm_fail_next_cold_update_before_publication() noexcept {
  fail_next_cold_update_before_publication.store(
      true, std::memory_order_relaxed);
}

void test::PreconditionerTestAccess::reset_cold_update_fault() noexcept {
  fail_next_cold_update_before_publication.store(
      false, std::memory_order_relaxed);
}
#endif

}  // namespace hundun::linear
