// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <memory>

#include "hundun/lin_system.hpp"
#include "hundun/rt_mpi_context.hpp"

namespace hundun::linear {

class BiCGStabSolver final : public LinearSolver {
 public:
  BiCGStabSolver(execution::ExecutionContext& execution_context,
                 const runtime::MpiContext& mpi_context);
  ~BiCGStabSolver() noexcept override;
  BiCGStabSolver(const BiCGStabSolver&) = delete;
  BiCGStabSolver& operator=(const BiCGStabSolver&) = delete;

  SolveReport solve(const LinearOperator& linear_operator,
                    Preconditioner& preconditioner,
                    execution::VectorView<const double> b,
                    execution::VectorView<double> x,
                    const SolveControl& control) const override;

 private:
  struct State;
  std::unique_ptr<State> state_;
};

}  // namespace hundun::linear
