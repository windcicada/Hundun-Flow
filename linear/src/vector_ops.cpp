// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/vector_ops.hpp"

#include "hundun/runtime/error.hpp"
#include "linear/src/vector_ops_detail.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>

namespace hundun::linear {
namespace detail {

double bad_dot_product_sentinel() noexcept {
  return std::numeric_limits<double>::infinity();
}

}  // namespace detail

namespace {

template <class T>
const double* validate_read_view(
    execution::VectorView<T> view,
    execution::BackendIdentity backend_identity) {
  static_assert(std::is_same_v<T, double> ||
                std::is_same_v<T, const double>);
  if (view.scalar_format() != execution::ScalarFormat::float64) {
    throw runtime::Error("VectorOps requires float64 vector views");
  }
  if (view.backend_identity() != backend_identity) {
    throw runtime::Error("VectorOps vector view backend does not match");
  }
  if (view.space() != execution::ExecutionSpace::host) {
    throw runtime::Error("VectorOps requires host-resident vector views");
  }
  return view.data();
}

double* validate_write_view(
    execution::VectorView<double> view,
    execution::BackendIdentity backend_identity) {
  if (!view.writable()) {
    throw runtime::Error("VectorOps output vector view is not writable");
  }
  const double* pointer = validate_read_view(view, backend_identity);
  return const_cast<double*>(pointer);
}

template <class Left, class Right>
bool same_view(execution::VectorView<Left> left,
               execution::VectorView<Right> right) noexcept {
  return left.allocation_identity() == right.allocation_identity() &&
         left.offset_bytes() == right.offset_bytes() &&
         left.size() == right.size() && left.stride() == right.stride();
}

template <class Left, class Right>
bool shares_allocation(execution::VectorView<Left> left,
                       execution::VectorView<Right> right) noexcept {
  return left.allocation_identity() == right.allocation_identity();
}

template <class Input>
void reject_nonidentical_output_alias(
    execution::VectorView<Input> input,
    execution::VectorView<double> output) {
  if (shares_allocation(input, output) && !same_view(input, output)) {
    throw runtime::Error(
        "VectorOps output aliases an input without identical view metadata");
  }
}

void require_finite(double value, const char* message) {
  if (!std::isfinite(value)) {
    throw runtime::Error(message);
  }
}

class CompensatedSum final {
 public:
  bool add(double value) noexcept {
    const double adjusted = value - correction_;
    const double updated = sum_ + adjusted;
    const double correction = (updated - sum_) - adjusted;
    if (!std::isfinite(adjusted) || !std::isfinite(updated) ||
        !std::isfinite(correction)) {
      return false;
    }
    sum_ = updated;
    correction_ = correction;
    return true;
  }

  double value() const noexcept { return sum_; }

 private:
  double sum_{0.0};
  double correction_{0.0};
};

}  // namespace

VectorOps::VectorOps(const execution::ExecutionContext& context)
    : backend_identity_(context.backend_identity()), space_(context.space()) {
  if (backend_identity_ == 0U) {
    throw runtime::Error("VectorOps context backend identity must be nonzero");
  }
  if (space_ != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access)) {
    throw runtime::Error("VectorOps requires a host-accessible context");
  }
}

void VectorOps::fill(execution::VectorView<double> output,
                     double value) const {
  require_finite(value, "VectorOps fill value must be finite");
  double* output_pointer = validate_write_view(output, backend_identity_);
  for (std::size_t index = 0; index < output.size(); ++index) {
    output_pointer[index * output.stride()] = value;
  }
}

void VectorOps::copy(execution::VectorView<const double> input,
                     execution::VectorView<double> output) const {
  const double* input_pointer = validate_read_view(input, backend_identity_);
  double* output_pointer = validate_write_view(output, backend_identity_);
  if (input.size() != output.size()) {
    throw runtime::Error("VectorOps copy input and output sizes differ");
  }
  reject_nonidentical_output_alias(input, output);
  for (std::size_t index = 0; index < input.size(); ++index) {
    require_finite(input_pointer[index * input.stride()],
                   "VectorOps copy input must be finite");
  }
  if (same_view(input, output)) {
    return;
  }
  for (std::size_t index = 0; index < input.size(); ++index) {
    output_pointer[index * output.stride()] =
        input_pointer[index * input.stride()];
  }
}

void VectorOps::scale(double alpha,
                      execution::VectorView<double> values) const {
  require_finite(alpha, "VectorOps scale coefficient must be finite");
  double* pointer = validate_write_view(values, backend_identity_);
  for (std::size_t index = 0; index < values.size(); ++index) {
    const double result = alpha * pointer[index * values.stride()];
    require_finite(pointer[index * values.stride()],
                   "VectorOps scale input must be finite");
    require_finite(result, "VectorOps scale result must be finite");
  }
  for (std::size_t index = 0; index < values.size(); ++index) {
    pointer[index * values.stride()] *= alpha;
  }
}

void VectorOps::axpy(double alpha,
                     execution::VectorView<const double> x,
                     execution::VectorView<double> y) const {
  require_finite(alpha, "VectorOps axpy coefficient must be finite");
  const double* x_pointer = validate_read_view(x, backend_identity_);
  double* y_pointer = validate_write_view(y, backend_identity_);
  if (x.size() != y.size()) {
    throw runtime::Error("VectorOps axpy input and output sizes differ");
  }
  reject_nonidentical_output_alias(x, y);
  for (std::size_t index = 0; index < x.size(); ++index) {
    const double x_value = x_pointer[index * x.stride()];
    const double y_value = y_pointer[index * y.stride()];
    require_finite(x_value, "VectorOps axpy input must be finite");
    require_finite(y_value, "VectorOps axpy output input must be finite");
    require_finite(alpha * x_value + y_value,
                   "VectorOps axpy result must be finite");
  }
  for (std::size_t index = 0; index < x.size(); ++index) {
    y_pointer[index * y.stride()] +=
        alpha * x_pointer[index * x.stride()];
  }
}

void VectorOps::linear_combination(
    double alpha, execution::VectorView<const double> x,
    double beta, execution::VectorView<const double> y,
    execution::VectorView<double> output) const {
  require_finite(alpha,
                 "VectorOps linear-combination alpha must be finite");
  require_finite(beta,
                 "VectorOps linear-combination beta must be finite");
  const double* x_pointer = validate_read_view(x, backend_identity_);
  const double* y_pointer = validate_read_view(y, backend_identity_);
  double* output_pointer = validate_write_view(output, backend_identity_);
  if (x.size() != y.size() || x.size() != output.size()) {
    throw runtime::Error(
        "VectorOps linear-combination input and output sizes differ");
  }
  reject_nonidentical_output_alias(x, output);
  reject_nonidentical_output_alias(y, output);
  for (std::size_t index = 0; index < x.size(); ++index) {
    const double x_value = x_pointer[index * x.stride()];
    const double y_value = y_pointer[index * y.stride()];
    require_finite(x_value,
                   "VectorOps linear-combination input must be finite");
    require_finite(y_value,
                   "VectorOps linear-combination input must be finite");
    require_finite(alpha * x_value + beta * y_value,
                   "VectorOps linear-combination result must be finite");
  }
  for (std::size_t index = 0; index < x.size(); ++index) {
    output_pointer[index * output.stride()] =
        alpha * x_pointer[index * x.stride()] +
        beta * y_pointer[index * y.stride()];
  }
}

double VectorOps::norm(execution::VectorView<const double> values,
                       const runtime::MpiContext& context) const {
  const double* pointer = nullptr;
  double local_scale = 0.0;
  bool has_local_source = false;
  try {
    pointer = validate_read_view(values, backend_identity_);
    for (std::size_t index = 0; index < values.size(); ++index) {
      const double value = pointer[index * values.stride()];
      if (!std::isfinite(value)) {
        has_local_source = true;
        local_scale = std::numeric_limits<double>::infinity();
        break;
      }
      local_scale = std::max(local_scale, std::abs(value));
    }
  } catch (const runtime::Error&) {
    has_local_source = true;
    local_scale = std::numeric_limits<double>::infinity();
  }

  context.allreduce_fp64_in_place(
      &local_scale, 1U, runtime::Fp64ReductionOperation::maximum);
  if (!std::isfinite(local_scale)) {
    throw detail::SynchronizedReductionError(
        "VectorOps norm input must be finite and valid", has_local_source);
  }

  CompensatedSum local_sum;
  bool local_sum_valid = true;
  if (local_scale != 0.0) {
    for (std::size_t index = 0; index < values.size(); ++index) {
      const double ratio = pointer[index * values.stride()] / local_scale;
      if (!std::isfinite(ratio) || !local_sum.add(ratio * ratio)) {
        has_local_source = true;
        local_sum_valid = false;
        break;
      }
    }
  }
  double global_sum = local_sum_valid
                          ? local_sum.value()
                          : std::numeric_limits<double>::infinity();
  context.allreduce_fp64_in_place(
      &global_sum, 1U, runtime::Fp64ReductionOperation::sum);
  if (!std::isfinite(global_sum) || global_sum < 0.0) {
    throw detail::SynchronizedReductionError(
        "VectorOps norm global scaled sum must be finite and nonnegative",
        has_local_source);
  }
  const double result = local_scale * std::sqrt(global_sum);
  if (!std::isfinite(result)) {
    throw detail::SynchronizedReductionError(
        "VectorOps norm result must be finite", false);
  }
  return result;
}

void VectorOps::dot_batch(
    const DotProductPair* pairs, std::size_t pair_count,
    execution::VectorView<double> results,
    const runtime::MpiContext& context) const {
  if (pair_count > static_cast<std::size_t>(INT_MAX)) {
    throw runtime::Error("VectorOps dot_batch pair count exceeds MPI INT_MAX");
  }
  if (pair_count != 0U && pairs == nullptr) {
    throw runtime::Error(
        "VectorOps dot_batch requires a non-null pair pointer");
  }
  if (results.size() != pair_count) {
    throw runtime::Error(
        "VectorOps dot_batch result size must equal the pair count");
  }
  if (results.stride() != 1U) {
    throw runtime::Error("VectorOps dot_batch result stride must be one");
  }
  double* result_pointer = validate_write_view(results, backend_identity_);
  if (pair_count == 0U) {
    return;
  }
  for (std::size_t pair = 0; pair < pair_count; ++pair) {
    if (shares_allocation(results, pairs[pair].left) ||
        shares_allocation(results, pairs[pair].right)) {
      throw runtime::Error(
          "VectorOps dot_batch result storage must be distinct from inputs");
    }
  }

  std::size_t first_local_source_pair = pair_count;
  for (std::size_t pair = 0; pair < pair_count; ++pair) {
    double local_result = detail::bad_dot_product_sentinel();
    try {
      const auto left = pairs[pair].left;
      const auto right = pairs[pair].right;
      const double* left_pointer =
          validate_read_view(left, backend_identity_);
      const double* right_pointer =
          validate_read_view(right, backend_identity_);
      if (left.size() != right.size()) {
        throw runtime::Error("VectorOps dot_batch pair sizes differ");
      }
      CompensatedSum local_sum;
      bool valid = true;
      for (std::size_t index = 0; index < left.size(); ++index) {
        const double left_value = left_pointer[index * left.stride()];
        const double right_value = right_pointer[index * right.stride()];
        const double product = left_value * right_value;
        if (!std::isfinite(left_value) || !std::isfinite(right_value) ||
            !std::isfinite(product) || !local_sum.add(product)) {
          valid = false;
          break;
        }
      }
      if (valid && std::isfinite(local_sum.value())) {
        local_result = local_sum.value();
      } else {
        first_local_source_pair =
            std::min(first_local_source_pair, pair);
      }
    } catch (const runtime::Error&) {
      first_local_source_pair = std::min(first_local_source_pair, pair);
      local_result = detail::bad_dot_product_sentinel();
    }
    result_pointer[pair] = local_result;
  }

  context.allreduce_fp64_in_place(
      result_pointer, pair_count, runtime::Fp64ReductionOperation::sum);
  for (std::size_t pair = 0; pair < pair_count; ++pair) {
    if (!std::isfinite(result_pointer[pair])) {
      throw detail::SynchronizedReductionError(
          "VectorOps dot_batch global result " + std::to_string(pair) +
              " must be finite",
          first_local_source_pair == pair);
    }
  }
}

}  // namespace hundun::linear
