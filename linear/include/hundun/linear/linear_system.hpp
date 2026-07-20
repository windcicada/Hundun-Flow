// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/execution/execution.hpp"
#include "hundun/linear/ghosted_vector.hpp"

#include <cstdint>
#include <limits>

namespace hundun::linear {

enum class SolveTerminationReason {
  converged,
  zero_right_hand_side,
  maximum_iterations,
  numerical_breakdown,
  non_finite_value,
  invalid_control,
  collective_failure
};

struct SolveControl final {
  double atol{1.0e-12};
  double rtol{1.0e-10};
  std::uint64_t max_iterations{500U};
  std::uint64_t residual_recompute_interval{20U};
};

struct SolveReport final {
  SolveTerminationReason reason{SolveTerminationReason::invalid_control};
  std::uint64_t iterations{};
  double initial_residual{std::numeric_limits<double>::infinity()};
  double recursive_residual{std::numeric_limits<double>::infinity()};
  double final_residual{std::numeric_limits<double>::infinity()};
  std::uint64_t matvec_count{};
  std::uint64_t preconditioner_apply_count{};
  std::uint64_t global_reduction_count{};
  int lowest_failing_rank{-1};
};

class LinearOperator {
 public:
  virtual ~LinearOperator() noexcept = default;
  virtual VectorLayout domain_layout() const = 0;
  virtual VectorLayout range_layout() const = 0;
  virtual const execution::ExecutionContext& context() const = 0;
  virtual std::uint64_t revision() const = 0;
  virtual execution::ExecutionEvent apply(
      execution::VectorView<const double> x,
      execution::VectorView<double> y) const = 0;
  virtual bool has_diagonal() const = 0;
  virtual execution::ExecutionEvent diagonal(
      execution::VectorView<double> d) const = 0;
};

class Preconditioner {
 public:
  virtual ~Preconditioner() noexcept = default;
  virtual void update(const LinearOperator& linear_operator,
                      std::uint64_t revision) = 0;
  virtual execution::ExecutionEvent apply(
      execution::VectorView<const double> r,
      execution::VectorView<double> z) const = 0;
};

class LinearSolver {
 public:
  virtual ~LinearSolver() noexcept = default;
  virtual SolveReport solve(
      const LinearOperator& linear_operator,
      Preconditioner& preconditioner,
      execution::VectorView<const double> b,
      execution::VectorView<double> x,
      const SolveControl& control) const = 0;
};

}  // namespace hundun::linear
