// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/conjugate_gradient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string>

#include "hundun/linear/vector_ops.hpp"
#include "hundun/runtime/error.hpp"
#include "linear/src/conjugate_gradient_detail.hpp"
#include "linear/src/solve_collective_detail.hpp"
#include "linear/src/vector_ops_detail.hpp"
#include "mpi_error.hpp"

namespace hundun::linear {
namespace {

constexpr std::size_t kWorkspaceVectorCount = 5U;

std::uint64_t checked_increment(std::uint64_t value, const char* name) {
  if (value == std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error(std::string("CG ") + name + " counter would overflow");
  }
  return value + 1U;
}

}  // namespace

namespace detail {

std::size_t checked_cg_workspace_bytes(std::size_t owned_count) {
  if (owned_count >
      std::numeric_limits<std::size_t>::max() / kWorkspaceVectorCount) {
    throw runtime::Error("CG workspace vector count overflow");
  }
  const std::size_t elements = owned_count * kWorkspaceVectorCount;
  if (elements > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw runtime::Error("CG workspace byte size overflow");
  }
  return elements * sizeof(double);
}

}  // namespace detail

struct ConjugateGradientSolver::State final {
  State(execution::ExecutionContext& execution_context_value,
        const runtime::MpiContext& mpi_context_value)
      : execution_context(&execution_context_value),
        mpi_context(&mpi_context_value),
        operations(execution_context_value) {}

  execution::ExecutionContext* execution_context;
  const runtime::MpiContext* mpi_context;
  VectorOps operations;
};

ConjugateGradientSolver::ConjugateGradientSolver(
    execution::ExecutionContext& execution_context,
    const runtime::MpiContext& mpi_context) {
  detail::validate_solver_execution_context(execution_context, "CG");
  detail::validate_solver_mpi_context(mpi_context, "CG");
  try {
    state_ = std::make_unique<State>(execution_context, mpi_context);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("CG solver state allocation failed");
  }
}

ConjugateGradientSolver::~ConjugateGradientSolver() noexcept = default;

SolveReport ConjugateGradientSolver::solve(
    const LinearOperator& linear_operator,
    Preconditioner& preconditioner,
    execution::VectorView<const double> b,
    execution::VectorView<double> x,
    const SolveControl& control) const {
  detail::validate_solver_mpi_context(*state_->mpi_context, "CG");
  const detail::SolveCollectiveProtocol protocol(*state_->mpi_context);
  SolveReport report;
  auto finish = [&](SolveTerminationReason reason, int rank = -1) {
    return detail::finish_report(protocol, report, reason, rank);
  };
  auto selected = protocol.agree_control(control);
  if (selected.failed) {
    return finish(selected.reason, selected.lowest_rank);
  }

  std::size_t owned_count = 0U;
  detail::LocalSolveFailure local_failure = detail::LocalSolveFailure::none;
  try {
    owned_count = detail::validate_solver_preflight(
        linear_operator, *state_->execution_context, b, x, "CG");
    if (!detail::all_finite(b) || !detail::all_finite(x)) {
      local_failure = detail::LocalSolveFailure::non_finite_value;
    }
  } catch (const runtime::detail::MpiOperationError&) {
    throw;
  } catch (...) {
    local_failure = detail::LocalSolveFailure::collective_failure;
  }
  selected = protocol.checkpoint(local_failure);
  if (selected.failed) {
    return finish(selected.reason, selected.lowest_rank);
  }

  std::array<std::optional<execution::Buffer>, kWorkspaceVectorCount> workspace;
  std::optional<execution::Buffer> scalar_buffer;
  try {
    const std::size_t workspace_bytes =
        detail::checked_cg_workspace_bytes(owned_count);
    const std::size_t vector_bytes = workspace_bytes / kWorkspaceVectorCount;
    for (auto& buffer : workspace) {
      buffer.emplace(*state_->execution_context, vector_bytes);
    }
    scalar_buffer.emplace(*state_->execution_context, sizeof(double));
  } catch (const runtime::detail::MpiOperationError&) {
    throw;
  } catch (...) {
    local_failure = detail::LocalSolveFailure::collective_failure;
  }
  selected = protocol.checkpoint(local_failure);
  if (selected.failed) {
    return finish(selected.reason, selected.lowest_rank);
  }

  auto r = workspace[0]->view(0U, owned_count);
  auto z = workspace[1]->view(0U, owned_count);
  auto p = workspace[2]->view(0U, owned_count);
  auto ap = workspace[3]->view(0U, owned_count);
  auto ax = workspace[4]->view(0U, owned_count);
  auto reduction = scalar_buffer->view(0U, 1U);

  auto phase_checkpoint = [&](detail::LocalSolveFailure failure) {
    selected = protocol.checkpoint(failure);
    return !selected.failed;
  };
  auto apply_operator = [&](execution::VectorView<const double> input,
                            execution::VectorView<double> output) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      report.matvec_count = checked_increment(report.matvec_count, "matvec");
      auto event = linear_operator.apply(input, output);
      event.wait();
      if (!detail::all_finite(output)) {
        failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (const runtime::detail::MpiOperationError&) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return phase_checkpoint(failure);
  };
  auto update_preconditioner = [&] {
    auto failure = detail::LocalSolveFailure::none;
    try {
      preconditioner.update(linear_operator, linear_operator.revision());
    } catch (const runtime::detail::MpiOperationError&) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return phase_checkpoint(failure);
  };
  auto apply_preconditioner = [&] {
    auto failure = detail::LocalSolveFailure::none;
    try {
      report.preconditioner_apply_count = checked_increment(
          report.preconditioner_apply_count, "preconditioner apply");
      auto event = preconditioner.apply(r, z);
      event.wait();
      if (!detail::all_finite(z)) {
        failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (const runtime::detail::MpiOperationError&) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return phase_checkpoint(failure);
  };
  auto vector_phase = [&](auto&& operation) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      operation();
    } catch (const runtime::detail::MpiOperationError&) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return phase_checkpoint(failure);
  };
  auto norm = [&](execution::VectorView<const double> values, double& result) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      result = state_->operations.norm(values, *state_->mpi_context);
      if (!std::isfinite(result)) {
        failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (const runtime::detail::MpiOperationError&) {
      throw;
    } catch (const detail::SynchronizedReductionError& error) {
      if (!phase_checkpoint(error.has_local_source()
                                ? detail::LocalSolveFailure::non_finite_value
                                : detail::LocalSolveFailure::none)) {
        return false;
      }
      selected = {SolveTerminationReason::non_finite_value, -1, true};
      return false;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return phase_checkpoint(failure);
  };
  auto dot = [&](execution::VectorView<const double> left,
                 execution::VectorView<const double> right, double& result) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      const DotProductPair pair{left, right};
      state_->operations.dot_batch(&pair, 1U, reduction, *state_->mpi_context);
      result = reduction[0];
      if (!std::isfinite(result)) {
        failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (const runtime::detail::MpiOperationError&) {
      throw;
    } catch (const detail::SynchronizedReductionError& error) {
      if (!phase_checkpoint(error.has_local_source()
                                ? detail::LocalSolveFailure::non_finite_value
                                : detail::LocalSolveFailure::none)) {
        return false;
      }
      selected = {SolveTerminationReason::non_finite_value, -1, true};
      return false;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return phase_checkpoint(failure);
  };

  auto independent_residual = [&](double& value) {
    if (!apply_operator(x, ax)) {
      return false;
    }
    auto failure = detail::LocalSolveFailure::none;
    try {
      if (!detail::finite_linear_candidate(1.0, b, -1.0, ax)) {
        failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    if (!phase_checkpoint(failure)) {
      return false;
    }
    if (!vector_phase([&] {
          state_->operations.linear_combination(1.0, b, -1.0, ax, r);
        })) {
      return false;
    }
    return norm(r, value);
  };

  double right_hand_side_norm = std::numeric_limits<double>::infinity();
  if (!norm(b, right_hand_side_norm)) {
    return finish(selected.reason, selected.lowest_rank);
  }
  const double relative_tolerance = control.rtol * right_hand_side_norm;
  const double tolerance = std::max(control.atol, relative_tolerance);
  if (!std::isfinite(relative_tolerance) || !std::isfinite(tolerance)) {
    return finish(SolveTerminationReason::non_finite_value);
  }
  if (!phase_checkpoint(detail::LocalSolveFailure::none)) {
    return finish(selected.reason, selected.lowest_rank);
  }

  double independent = std::numeric_limits<double>::infinity();
  if (!independent_residual(independent)) {
    return finish(selected.reason, selected.lowest_rank);
  }
  report.initial_residual = independent;
  report.recursive_residual = independent;
  report.final_residual = independent;
  bool independent_is_current = true;

  if (right_hand_side_norm == 0.0) {
    report.final_residual = std::numeric_limits<double>::infinity();
    independent_is_current = false;
    if (!vector_phase([&] { state_->operations.fill(x, 0.0); })) {
      return finish(selected.reason, selected.lowest_rank);
    }
    if (!independent_residual(independent)) {
      return finish(selected.reason, selected.lowest_rank);
    }
    independent_is_current = true;
    report.final_residual = independent;
    return finish(independent <= control.atol
                      ? SolveTerminationReason::zero_right_hand_side
                      : SolveTerminationReason::numerical_breakdown);
  }
  if (independent <= tolerance) {
    return finish(SolveTerminationReason::converged);
  }
  if (control.max_iterations == 0U) {
    return finish(SolveTerminationReason::maximum_iterations);
  }

  if (!update_preconditioner()) {
    return finish(selected.reason, selected.lowest_rank);
  }
  if (!apply_preconditioner()) {
    return finish(selected.reason, selected.lowest_rank);
  }
  double rho = 0.0;
  if (!dot(r, z, rho)) {
    return finish(selected.reason, selected.lowest_rank);
  }
  if (rho <= 0.0) {
    return finish(SolveTerminationReason::numerical_breakdown);
  }
  if (!vector_phase([&] { state_->operations.copy(z, p); })) {
    return finish(selected.reason, selected.lowest_rank);
  }

  auto refresh_final = [&] {
    if (!independent_is_current) {
      if (!independent_residual(independent)) {
        return false;
      }
      independent_is_current = true;
      report.final_residual = independent;
    }
    return true;
  };
  auto finish_failure_after_refresh =
      [&](detail::FailureSelection original_failure) {
        if (!refresh_final()) {
          return finish(selected.reason, selected.lowest_rank);
        }
        return finish(original_failure.reason, original_failure.lowest_rank);
      };

  while (report.iterations < control.max_iterations) {
    if (!apply_operator(p, ap)) {
      return finish(selected.reason, selected.lowest_rank);
    }
    double curvature = 0.0;
    if (!dot(p, ap, curvature)) {
      return finish_failure_after_refresh(selected);
    }
    if (curvature <= 0.0) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    const double alpha = rho / curvature;
    if (!std::isfinite(alpha)) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::non_finite_value);
    }
    local_failure = detail::LocalSolveFailure::none;
    try {
      if (!detail::finite_axpy_candidate(alpha, p, x) ||
          !detail::finite_axpy_candidate(-alpha, ap, r)) {
        local_failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (...) {
      local_failure = detail::LocalSolveFailure::collective_failure;
    }
    if (!phase_checkpoint(local_failure)) {
      return finish(selected.reason, selected.lowest_rank);
    }
    report.final_residual = std::numeric_limits<double>::infinity();
    independent_is_current = false;
    if (!vector_phase([&] { state_->operations.axpy(alpha, p, x); })) {
      return finish(selected.reason, selected.lowest_rank);
    }
    if (!vector_phase([&] { state_->operations.axpy(-alpha, ap, r); })) {
      return finish(selected.reason, selected.lowest_rank);
    }
    report.iterations = checked_increment(report.iterations, "iteration");
    if (!norm(r, report.recursive_residual)) {
      return finish(selected.reason, selected.lowest_rank);
    }

    const bool periodic =
        report.iterations % control.residual_recompute_interval == 0U;
    const bool candidate = report.recursive_residual <= tolerance;
    const bool exhausted = report.iterations == control.max_iterations;
    if (periodic || candidate || exhausted) {
      if (!independent_residual(independent)) {
        return finish(selected.reason, selected.lowest_rank);
      }
      independent_is_current = true;
      report.final_residual = independent;
      if (independent <= tolerance) {
        return finish(SolveTerminationReason::converged);
      }
      if (exhausted) {
        return finish(SolveTerminationReason::maximum_iterations);
      }
      report.recursive_residual = independent;
      if (!apply_preconditioner()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      if (!dot(r, z, rho)) {
        return finish(selected.reason, selected.lowest_rank);
      }
      if (rho <= 0.0) {
        return finish(SolveTerminationReason::numerical_breakdown);
      }
      if (!vector_phase([&] { state_->operations.copy(z, p); })) {
        return finish(selected.reason, selected.lowest_rank);
      }
      continue;
    }

    if (!apply_preconditioner()) {
      return finish_failure_after_refresh(selected);
    }
    double rho_new = 0.0;
    if (!dot(r, z, rho_new)) {
      return finish_failure_after_refresh(selected);
    }
    if (rho_new <= 0.0) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    const double beta = rho_new / rho;
    if (!std::isfinite(beta)) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::non_finite_value);
    }
    local_failure = detail::LocalSolveFailure::none;
    try {
      if (!detail::finite_linear_candidate(1.0, z, beta, p)) {
        local_failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (...) {
      local_failure = detail::LocalSolveFailure::collective_failure;
    }
    if (!phase_checkpoint(local_failure)) {
      return finish_failure_after_refresh(selected);
    }
    if (!vector_phase([&] {
          state_->operations.linear_combination(1.0, z, beta, p, p);
        })) {
      return finish_failure_after_refresh(selected);
    }
    rho = rho_new;
  }

  if (!refresh_final()) {
    return finish(selected.reason, selected.lowest_rank);
  }
  return finish(SolveTerminationReason::maximum_iterations);
}

}  // namespace hundun::linear
