// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/lin_system.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <cstdint>
#include <memory>

namespace hundun::linear {

// Stores each preconditioned Arnoldi basis vector separately, so callers may
// use either a fixed right preconditioner (GMRES) or a varying approximate
// inverse (FGMRES) without changing the Krylov update.
class RestartedGmresSolver final : public LinearSolver {
public:
  RestartedGmresSolver(execution::ExecutionContext &execution_context,
                       const runtime::MpiContext &mpi_context,
                       std::uint32_t restart_length = 30U);
  ~RestartedGmresSolver() noexcept override;
  RestartedGmresSolver(const RestartedGmresSolver &) = delete;
  RestartedGmresSolver &operator=(const RestartedGmresSolver &) = delete;

  SolveReport solve(const LinearOperator &linear_operator,
                    Preconditioner &preconditioner,
                    execution::VectorView<const double> b,
                    execution::VectorView<double> x,
                    const SolveControl &control) const override;

private:
  struct State;
  std::unique_ptr<State> state_;
};

using RestartedFgmresSolver = RestartedGmresSolver;

} // namespace hundun::linear
