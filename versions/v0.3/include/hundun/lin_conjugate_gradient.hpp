// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/lin_system.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <memory>

namespace hundun::linear {

class ConjugateGradientSolver final : public LinearSolver {
 public:
  ConjugateGradientSolver(
      execution::ExecutionContext& execution_context,
      const runtime::MpiContext& mpi_context);
  ~ConjugateGradientSolver() noexcept override;
  ConjugateGradientSolver(const ConjugateGradientSolver&) = delete;
  ConjugateGradientSolver& operator=(
      const ConjugateGradientSolver&) = delete;

  SolveReport solve(
      const LinearOperator& linear_operator,
      Preconditioner& preconditioner,
      execution::VectorView<const double> b,
      execution::VectorView<double> x,
      const SolveControl& control) const override;

 private:
  struct State;
  std::unique_ptr<State> state_;
};

}  // namespace hundun::linear
