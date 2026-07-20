// SPDX-License-Identifier: Apache-2.0

#include "linear/src/solve_collective_detail.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>

#include "hundun/runtime/error.hpp"

namespace hundun::linear::detail {
namespace {

constexpr std::uint64_t kFailureStride = 16U;
constexpr std::size_t kControlWordCount = 8U;

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t result = 0U;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void split_word(std::uint64_t word, double& low, double& high) noexcept {
  low = static_cast<double>(word & UINT64_C(0xffffffff));
  high = static_cast<double>(word >> 32U);
}

std::uint64_t join_word(double low, double high) {
  constexpr double kMaxHalf = 4294967295.0;
  if (!std::isfinite(low) || !std::isfinite(high) || low < 0.0 || high < 0.0 ||
      low > kMaxHalf || high > kMaxHalf || std::floor(low) != low ||
      std::floor(high) != high) {
    throw runtime::Error("solve control agreement produced corrupt halves");
  }
  return static_cast<std::uint64_t>(low) |
         (static_cast<std::uint64_t>(high) << 32U);
}

SolveTerminationReason reason_for(LocalSolveFailure failure) {
  switch (failure) {
    case LocalSolveFailure::invalid_control:
      return SolveTerminationReason::invalid_control;
    case LocalSolveFailure::collective_failure:
      return SolveTerminationReason::collective_failure;
    case LocalSolveFailure::non_finite_value:
      return SolveTerminationReason::non_finite_value;
    case LocalSolveFailure::none:
      break;
  }
  throw runtime::Error("solve failure protocol decoded an empty failure");
}

bool valid_control(const SolveControl& control) noexcept {
  return std::isfinite(control.atol) && control.atol >= 0.0 &&
         std::isfinite(control.rtol) && control.rtol >= 0.0 &&
         control.residual_recompute_interval != 0U;
}

template <class T>
void validate_view(execution::VectorView<T> view,
                   const execution::ExecutionContext& context,
                   std::size_t expected_size,
                   bool require_writable,
                   const char* solver_name,
                   const char* view_name) {
  if (view.scalar_format() != execution::ScalarFormat::float64) {
    throw runtime::Error(std::string(solver_name) + " " + view_name +
                         " view must be float64");
  }
  if (view.size() != expected_size) {
    throw runtime::Error(std::string(solver_name) + " " + view_name +
                         " size does not match the owned layout");
  }
  if (view.backend_identity() != context.backend_identity()) {
    throw runtime::Error(std::string(solver_name) + " " + view_name +
                         " backend identity does not match");
  }
  if (view.space() != context.space() ||
      view.space() != execution::ExecutionSpace::host) {
    throw runtime::Error(std::string(solver_name) + " " + view_name +
                         " execution space does not match host context");
  }
  if (require_writable && !view.writable()) {
    throw runtime::Error(std::string(solver_name) +
                         " solution view must be writable");
  }
  static_cast<void>(view.data());
}

}  // namespace

SolveCollectiveProtocol::SolveCollectiveProtocol(
    const runtime::MpiContext& context)
    : context_(&context),
      reductions_before_(context.fp64_reduction_counters().collective_calls) {}

FailureSelection SolveCollectiveProtocol::checkpoint(
    LocalSolveFailure local_failure) const {
  const int rank = context_->rank();
  const int size = context_->size();
  const auto code = static_cast<std::uint64_t>(local_failure);
  if (code >= kFailureStride || rank < 0 || rank >= size || size <= 0) {
    throw runtime::Error("solve failure protocol local state is invalid");
  }
  double encoded = 0.0;
  if (local_failure != LocalSolveFailure::none) {
    const std::uint64_t rank_priority = static_cast<std::uint64_t>(size - rank);
    if (rank_priority > (UINT64_C(1) << 53U) / kFailureStride) {
      throw runtime::Error("solve failure protocol rank encoding overflow");
    }
    encoded = static_cast<double>(rank_priority * kFailureStride + code);
  }
  context_->allreduce_fp64_in_place(&encoded, 1U,
                                    runtime::Fp64ReductionOperation::maximum);
  if (encoded == 0.0) {
    return {};
  }
  if (!std::isfinite(encoded) || encoded < 0.0 ||
      std::floor(encoded) != encoded ||
      encoded >= static_cast<double>(UINT64_C(1) << 53U)) {
    throw runtime::Error("solve failure protocol encoding is corrupt");
  }
  const auto integer = static_cast<std::uint64_t>(encoded);
  const auto decoded_code = integer % kFailureStride;
  const auto priority = integer / kFailureStride;
  if (decoded_code == 0U ||
      decoded_code >
          static_cast<std::uint64_t>(LocalSolveFailure::non_finite_value) ||
      priority == 0U || priority > static_cast<std::uint64_t>(size)) {
    throw runtime::Error("solve failure protocol encoding cannot be decoded");
  }
  const int decoded_rank = size - static_cast<int>(priority);
  const auto decoded_failure = static_cast<LocalSolveFailure>(decoded_code);
  return {reason_for(decoded_failure), decoded_rank, true};
}

FailureSelection SolveCollectiveProtocol::agree_control(
    const SolveControl& control) const {
  auto result =
      checkpoint(valid_control(control) ? LocalSolveFailure::none
                                        : LocalSolveFailure::invalid_control);
  if (result.failed) {
    return result;
  }

  std::array<double, kControlWordCount> root_words{};
  if (context_->rank() == 0) {
    split_word(double_bits(control.atol), root_words[0], root_words[1]);
    split_word(double_bits(control.rtol), root_words[2], root_words[3]);
    split_word(control.max_iterations, root_words[4], root_words[5]);
    split_word(control.residual_recompute_interval, root_words[6],
               root_words[7]);
  }
  context_->allreduce_fp64_in_place(root_words.data(), root_words.size(),
                                    runtime::Fp64ReductionOperation::sum);
  const bool matches =
      join_word(root_words[0], root_words[1]) == double_bits(control.atol) &&
      join_word(root_words[2], root_words[3]) == double_bits(control.rtol) &&
      join_word(root_words[4], root_words[5]) == control.max_iterations &&
      join_word(root_words[6], root_words[7]) ==
          control.residual_recompute_interval;
  return checkpoint(matches ? LocalSolveFailure::none
                            : LocalSolveFailure::collective_failure);
}

std::uint64_t SolveCollectiveProtocol::reduction_delta() const {
  const auto after = context_->fp64_reduction_counters().collective_calls;
  if (after < reductions_before_) {
    throw runtime::Error("solve global reduction counter decreased");
  }
  return after - reductions_before_;
}

void validate_solver_execution_context(
    const execution::ExecutionContext& context, const char* solver_name) {
  if (context.backend_identity() == 0U) {
    throw runtime::Error(std::string(solver_name) +
                         " execution context identity must be nonzero");
  }
  if (context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access)) {
    throw runtime::Error(std::string(solver_name) +
                         " requires a host-accessible execution context");
  }
  if (!context.supports(execution::ExecutionCapability::buffer_allocation)) {
    throw runtime::Error(std::string(solver_name) +
                         " execution context lacks buffer capability");
  }
  if (!context.supports(execution::ExecutionCapability::transfer)) {
    throw runtime::Error(std::string(solver_name) +
                         " execution context lacks transfer capability");
  }
}

void validate_solver_mpi_context(const runtime::MpiContext& context,
                                 const char* solver_name) {
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS ||
      MPI_Finalized(&finalized) != MPI_SUCCESS) {
    throw runtime::Error(std::string(solver_name) +
                         " could not query MPI context activity");
  }
  if (initialized == 0 || finalized != 0 || context.comm() == MPI_COMM_NULL ||
      context.size() <= 0 || context.rank() < 0 ||
      context.rank() >= context.size()) {
    throw runtime::Error(std::string(solver_name) +
                         " requires a live nonempty MPI context");
  }
}

std::size_t validate_solver_preflight(
    const LinearOperator& linear_operator,
    const execution::ExecutionContext& context,
    execution::VectorView<const double> b,
    execution::VectorView<double> x,
    const char* solver_name) {
  if (&linear_operator.context() != &context) {
    throw runtime::Error(std::string(solver_name) +
                         " operator context object does not match");
  }
  const VectorLayout domain = linear_operator.domain_layout();
  const VectorLayout range = linear_operator.range_layout();
  if (domain != range) {
    throw runtime::Error(std::string(solver_name) +
                         " requires an exactly square operator layout");
  }
  validate_view(b, context, domain.owned_count(), false, solver_name,
                "right-hand-side");
  validate_view(x, context, domain.owned_count(), true, solver_name,
                "solution");
  if (b.allocation_identity() == x.allocation_identity()) {
    throw runtime::Error(
        std::string(solver_name) +
        " right-hand-side and solution must not share an allocation");
  }
  return domain.owned_count();
}

bool finite_axpy_candidate(double alpha,
                           execution::VectorView<const double> input,
                           execution::VectorView<const double> current) {
  const double* input_pointer = input.data();
  const double* current_pointer = current.data();
  for (std::size_t index = 0; index < input.size(); ++index) {
    if (!std::isfinite(alpha * input_pointer[index * input.stride()] +
                       current_pointer[index * current.stride()])) {
      return false;
    }
  }
  return true;
}

bool finite_linear_candidate(double alpha,
                             execution::VectorView<const double> left,
                             double beta,
                             execution::VectorView<const double> right) {
  const double* left_pointer = left.data();
  const double* right_pointer = right.data();
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (!std::isfinite(alpha * left_pointer[index * left.stride()] +
                       beta * right_pointer[index * right.stride()])) {
      return false;
    }
  }
  return true;
}

bool finite_bicgstab_solution_candidate(
    execution::VectorView<const double> current,
    double alpha,
    execution::VectorView<const double> p_hat,
    double omega,
    execution::VectorView<const double> s_hat) {
  const double* current_pointer = current.data();
  const double* p_pointer = p_hat.data();
  const double* s_pointer = s_hat.data();
  for (std::size_t index = 0; index < current.size(); ++index) {
    const double candidate = current_pointer[index * current.stride()] +
                             alpha * p_pointer[index * p_hat.stride()] +
                             omega * s_pointer[index * s_hat.stride()];
    if (!std::isfinite(candidate)) {
      return false;
    }
  }
  return true;
}

SolveReport finish_report(const SolveCollectiveProtocol& protocol,
                          SolveReport report,
                          SolveTerminationReason reason,
                          int lowest_failing_rank) {
  report.reason = reason;
  report.lowest_failing_rank = lowest_failing_rank;
  report.global_reduction_count = protocol.reduction_delta();
  return report;
}

}  // namespace hundun::linear::detail
