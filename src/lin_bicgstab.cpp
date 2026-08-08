// SPDX-License-Identifier: Apache-2.0

#include "hundun/lin_bicgstab.hpp"

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

#include "hundun/lin_vector_ops.hpp"
#include "hundun/rt_error.hpp"
#include "lin_bicgstab_detail.hpp"
#include "lin_solve_collective_detail.hpp"
#include "lin_vector_ops_detail.hpp"
#include "rt_mpi_error_detail.hpp"

namespace hundun::linear {
namespace {

constexpr std::size_t kWorkspaceVectorCount = 9U;

std::uint64_t checked_increment(std::uint64_t value, const char* name) {
  if (value == std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error(std::string("BiCGStab ") + name +
                         " counter would overflow");
  }
  return value + 1U;
}

}  // namespace

namespace detail {

std::size_t checked_bicgstab_workspace_bytes(std::size_t owned_count) {
  if (owned_count >
      std::numeric_limits<std::size_t>::max() / kWorkspaceVectorCount) {
    throw runtime::Error("BiCGStab workspace vector count overflow");
  }
  const std::size_t elements = owned_count * kWorkspaceVectorCount;
  if (elements > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw runtime::Error("BiCGStab workspace byte size overflow");
  }
  return elements * sizeof(double);
}

}  // namespace detail

struct BiCGStabSolver::State final {
  State(execution::ExecutionContext& execution_context_value,
        const runtime::MpiContext& mpi_context_value)
      : execution_context(&execution_context_value),
        mpi_context(&mpi_context_value),
        operations(execution_context_value) {}

  execution::ExecutionContext* execution_context;
  const runtime::MpiContext* mpi_context;
  VectorOps operations;
};

BiCGStabSolver::BiCGStabSolver(execution::ExecutionContext& execution_context,
                               const runtime::MpiContext& mpi_context) {
  detail::validate_solver_execution_context(execution_context, "BiCGStab");
  detail::validate_solver_mpi_context(mpi_context, "BiCGStab");
  try {
    state_ = std::make_unique<State>(execution_context, mpi_context);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("BiCGStab solver state allocation failed");
  }
}

BiCGStabSolver::~BiCGStabSolver() noexcept = default;

SolveReport BiCGStabSolver::solve(const LinearOperator& linear_operator,
                                  Preconditioner& preconditioner,
                                  execution::VectorView<const double> b,
                                  execution::VectorView<double> x,
                                  const SolveControl& control) const {
  detail::validate_solver_mpi_context(*state_->mpi_context, "BiCGStab");
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
  auto local_failure = detail::LocalSolveFailure::none;
  try {
    owned_count = detail::validate_solver_preflight(
        linear_operator, *state_->execution_context, b, x, "BiCGStab");
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
        detail::checked_bicgstab_workspace_bytes(owned_count);
    const std::size_t vector_bytes = workspace_bytes / kWorkspaceVectorCount;
    for (auto& buffer : workspace) {
      buffer.emplace(*state_->execution_context, vector_bytes);
    }
    scalar_buffer.emplace(*state_->execution_context, 2U * sizeof(double));
  } catch (const runtime::detail::MpiOperationError&) {
    throw;
  } catch (...) {
    local_failure = detail::LocalSolveFailure::collective_failure;
  }
  selected = protocol.checkpoint(local_failure);
  if (selected.failed) {
    return finish(selected.reason, selected.lowest_rank);
  }

  const std::size_t n = owned_count;
  auto r = workspace[0]->view(0U, n);
  auto r_hat = workspace[1]->view(0U, n);
  auto p = workspace[2]->view(0U, n);
  auto v = workspace[3]->view(0U, n);
  auto s = workspace[4]->view(0U, n);
  auto t = workspace[5]->view(0U, n);
  auto p_hat = workspace[6]->view(0U, n);
  auto s_hat = workspace[7]->view(0U, n);
  auto ax = workspace[8]->view(0U, n);
  auto reduction_one = scalar_buffer->view(0U, 1U);
  auto reduction_two = scalar_buffer->view(0U, 2U);

  auto phase_checkpoint = [&](detail::LocalSolveFailure failure) {
    selected = protocol.checkpoint(failure);
    return !selected.failed;
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
  auto apply_preconditioner = [&](execution::VectorView<const double> input,
                                  execution::VectorView<double> output) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      report.preconditioner_apply_count = checked_increment(
          report.preconditioner_apply_count, "preconditioner apply");
      auto event = preconditioner.apply(input, output);
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
      state_->operations.dot_batch(&pair, 1U, reduction_one,
                                   *state_->mpi_context);
      result = reduction_one[0];
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
  auto dot_two = [&](double& ts, double& tt) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      const DotProductPair pairs[2]{{t, s}, {t, t}};
      state_->operations.dot_batch(pairs, 2U, reduction_two,
                                   *state_->mpi_context);
      ts = reduction_two[0];
      tt = reduction_two[1];
      if (!std::isfinite(ts) || !std::isfinite(tt)) {
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

  double b_norm = std::numeric_limits<double>::infinity();
  if (!norm(b, b_norm)) {
    return finish(selected.reason, selected.lowest_rank);
  }
  const double relative_tolerance = control.rtol * b_norm;
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

  if (b_norm == 0.0) {
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
  if (!vector_phase([&] {
        state_->operations.copy(r, r_hat);
        state_->operations.fill(p, 0.0);
        state_->operations.fill(v, 0.0);
      })) {
    return finish(selected.reason, selected.lowest_rank);
  }

  double rho_old = 1.0;
  double alpha = 1.0;
  double omega = 1.0;
  bool first_or_restart = true;
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
  auto restart = [&] {
    return vector_phase([&] {
      state_->operations.copy(r, r_hat);
      state_->operations.fill(p, 0.0);
      state_->operations.fill(v, 0.0);
      rho_old = 1.0;
      alpha = 1.0;
      omega = 1.0;
      first_or_restart = true;
    });
  };

  while (report.iterations < control.max_iterations) {
    double rho = 0.0;
    if (!dot(r_hat, r, rho)) {
      return finish_failure_after_refresh(selected);
    }
    if (rho == 0.0) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    if (first_or_restart) {
      if (!vector_phase([&] { state_->operations.copy(r, p); })) {
        return finish(selected.reason, selected.lowest_rank);
      }
    } else {
      const double beta = (rho / rho_old) * (alpha / omega);
      if (!std::isfinite(beta)) {
        if (!refresh_final()) {
          return finish(selected.reason, selected.lowest_rank);
        }
        return finish(SolveTerminationReason::non_finite_value);
      }
      local_failure = detail::LocalSolveFailure::none;
      try {
        const double* rp = r.data();
        const double* pp = p.data();
        const double* vp = v.data();
        for (std::size_t index = 0; index < n; ++index) {
          if (!std::isfinite(rp[index * r.stride()] +
                             beta * (pp[index * p.stride()] -
                                     omega * vp[index * v.stride()]))) {
            local_failure = detail::LocalSolveFailure::non_finite_value;
            break;
          }
        }
      } catch (...) {
        local_failure = detail::LocalSolveFailure::collective_failure;
      }
      if (!phase_checkpoint(local_failure)) {
        return finish_failure_after_refresh(selected);
      }
      if (!vector_phase([&] { state_->operations.axpy(-omega, v, p); }) ||
          !vector_phase([&] {
            state_->operations.linear_combination(1.0, r, beta, p, p);
          })) {
        return finish_failure_after_refresh(selected);
      }
    }

    if (!apply_preconditioner(p, p_hat) || !apply_operator(p_hat, v)) {
      return finish(selected.reason, selected.lowest_rank);
    }
    double denominator = 0.0;
    if (!dot(r_hat, v, denominator)) {
      return finish_failure_after_refresh(selected);
    }
    if (denominator == 0.0) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    alpha = rho / denominator;
    if (!std::isfinite(alpha)) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::non_finite_value);
    }
    local_failure = detail::LocalSolveFailure::none;
    try {
      if (!detail::finite_linear_candidate(1.0, r, -alpha, v)) {
        local_failure = detail::LocalSolveFailure::non_finite_value;
      }
    } catch (...) {
      local_failure = detail::LocalSolveFailure::collective_failure;
    }
    if (!phase_checkpoint(local_failure)) {
      return finish_failure_after_refresh(selected);
    }
    if (!vector_phase([&] {
          state_->operations.linear_combination(1.0, r, -alpha, v, s);
        })) {
      return finish(selected.reason, selected.lowest_rank);
    }
    double s_norm = 0.0;
    if (!norm(s, s_norm)) {
      return finish(selected.reason, selected.lowest_rank);
    }

    if (s_norm <= tolerance) {
      local_failure = detail::LocalSolveFailure::none;
      try {
        if (!detail::finite_axpy_candidate(alpha, p_hat, x)) {
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
      if (!vector_phase([&] { state_->operations.axpy(alpha, p_hat, x); })) {
        return finish(selected.reason, selected.lowest_rank);
      }
      report.iterations = checked_increment(report.iterations, "iteration");
      if (!independent_residual(independent)) {
        return finish(selected.reason, selected.lowest_rank);
      }
      independent_is_current = true;
      report.final_residual = independent;
      report.recursive_residual = independent;
      if (independent <= tolerance) {
        return finish(SolveTerminationReason::converged);
      }
      if (report.iterations == control.max_iterations) {
        return finish(SolveTerminationReason::maximum_iterations);
      }
      if (!restart()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      continue;
    }

    if (!apply_preconditioner(s, s_hat) || !apply_operator(s_hat, t)) {
      return finish(selected.reason, selected.lowest_rank);
    }
    double ts = 0.0;
    double tt = 0.0;
    if (!dot_two(ts, tt)) {
      return finish_failure_after_refresh(selected);
    }
    if (tt <= 0.0) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    omega = ts / tt;
    if (!std::isfinite(omega)) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::non_finite_value);
    }
    local_failure = detail::LocalSolveFailure::none;
    if (omega == 0.0) {
      if (!refresh_final()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    if (local_failure == detail::LocalSolveFailure::none) {
      try {
        if (!detail::finite_bicgstab_solution_candidate(x, alpha, p_hat, omega,
                                                        s_hat) ||
            !detail::finite_linear_candidate(1.0, s, -omega, t)) {
          local_failure = detail::LocalSolveFailure::non_finite_value;
        }
      } catch (...) {
        local_failure = detail::LocalSolveFailure::collective_failure;
      }
    }
    if (!phase_checkpoint(local_failure)) {
      return finish_failure_after_refresh(selected);
    }
    report.final_residual = std::numeric_limits<double>::infinity();
    independent_is_current = false;
    if (!vector_phase([&] { state_->operations.axpy(alpha, p_hat, x); })) {
      return finish(selected.reason, selected.lowest_rank);
    }
    if (!vector_phase([&] { state_->operations.axpy(omega, s_hat, x); }) ||
        !vector_phase([&] {
          state_->operations.linear_combination(1.0, s, -omega, t, r);
        })) {
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
      report.recursive_residual = independent;
      if (independent <= tolerance) {
        return finish(SolveTerminationReason::converged);
      }
      if (exhausted) {
        return finish(SolveTerminationReason::maximum_iterations);
      }
      if (!restart()) {
        return finish(selected.reason, selected.lowest_rank);
      }
      continue;
    }
    rho_old = rho;
    first_or_restart = false;
  }

  if (!refresh_final()) {
    return finish(selected.reason, selected.lowest_rank);
  }
  return finish(SolveTerminationReason::maximum_iterations);
}

}  // namespace hundun::linear
