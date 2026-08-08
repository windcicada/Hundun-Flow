// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "hundun/lin_system.hpp"
#include "hundun/rt_mpi_context.hpp"

namespace hundun::linear::detail {

enum class LocalSolveFailure : unsigned char {
  none = 0,
  invalid_control = 1,
  collective_failure = 2,
  non_finite_value = 3
};

struct FailureSelection final {
  SolveTerminationReason reason{SolveTerminationReason::collective_failure};
  int lowest_rank{-1};
  bool failed{false};
};

class SolveCollectiveProtocol final {
 public:
  explicit SolveCollectiveProtocol(const runtime::MpiContext& context);

  FailureSelection agree_control(const SolveControl& control) const;
  FailureSelection checkpoint(LocalSolveFailure local_failure) const;
  std::uint64_t reduction_delta() const;

 private:
  const runtime::MpiContext* context_;
  std::uint64_t reductions_before_;
};

void validate_solver_execution_context(
    const execution::ExecutionContext& context, const char* solver_name);
void validate_solver_mpi_context(const runtime::MpiContext& context,
                                 const char* solver_name);
std::size_t validate_solver_preflight(
    const LinearOperator& linear_operator,
    const execution::ExecutionContext& context,
    execution::VectorView<const double> b,
    execution::VectorView<double> x,
    const char* solver_name);

template <class T>
bool all_finite(execution::VectorView<T> values) {
  const double* pointer = values.data();
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(pointer[index * values.stride()])) {
      return false;
    }
  }
  return true;
}

bool finite_axpy_candidate(double alpha,
                           execution::VectorView<const double> input,
                           execution::VectorView<const double> current);
bool finite_linear_candidate(double alpha,
                             execution::VectorView<const double> left,
                             double beta,
                             execution::VectorView<const double> right);
bool finite_bicgstab_solution_candidate(
    execution::VectorView<const double> current,
    double alpha,
    execution::VectorView<const double> p_hat,
    double omega,
    execution::VectorView<const double> s_hat);

SolveReport finish_report(const SolveCollectiveProtocol& protocol,
                          SolveReport report,
                          SolveTerminationReason reason,
                          int lowest_failing_rank = -1);

}  // namespace hundun::linear::detail
