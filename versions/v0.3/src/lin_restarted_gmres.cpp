// SPDX-License-Identifier: Apache-2.0

#include "hundun/lin_restarted_gmres.hpp"

#include "hundun/rt_error.hpp"
#include "lin_solve_collective_detail.hpp"
#include "rt_mpi_error_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string>
#include <vector>

namespace hundun::linear {
namespace {

std::uint64_t checked_increment(std::uint64_t value, const char *name) {
  if (value == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error(std::string("GMRES ") + name +
                         " counter would overflow");
  return value + 1U;
}

std::size_t checked_workspace_count(std::uint32_t restart_length) {
  const auto length = static_cast<std::size_t>(restart_length);
  if (length > (std::numeric_limits<std::size_t>::max() - 4U) / 2U)
    throw runtime::Error("GMRES restart length overflows workspace count");
  return 2U * length + 4U;
}

std::size_t checked_vector_bytes(std::size_t count) {
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(double))
    throw runtime::Error("GMRES vector byte size overflows");
  return count * sizeof(double);
}

} // namespace

struct RestartedGmresSolver::State final {
  State(execution::ExecutionContext &execution_context_value,
        const runtime::MpiContext &mpi_context_value,
        std::uint32_t restart_length_value)
      : execution_context(&execution_context_value),
        mpi_context(&mpi_context_value), restart_length(restart_length_value) {}

  execution::ExecutionContext *execution_context{};
  const runtime::MpiContext *mpi_context{};
  std::uint32_t restart_length{};
  std::vector<std::optional<execution::Buffer>> workspace;
  std::size_t workspace_owned_count{std::numeric_limits<std::size_t>::max()};
};

RestartedGmresSolver::RestartedGmresSolver(
    execution::ExecutionContext &execution_context,
    const runtime::MpiContext &mpi_context, std::uint32_t restart_length) {
  detail::validate_solver_execution_context(execution_context, "GMRES");
  detail::validate_solver_mpi_context(mpi_context, "GMRES");
  if (restart_length == 0U)
    throw runtime::Error("GMRES restart length must be positive");
  try {
    state_ = std::make_unique<State>(execution_context, mpi_context,
                                     restart_length);
  } catch (const std::bad_alloc &) {
    throw runtime::Error("GMRES solver state allocation failed");
  }
}

RestartedGmresSolver::~RestartedGmresSolver() noexcept = default;

SolveReport RestartedGmresSolver::solve(
    const LinearOperator &linear_operator, Preconditioner &preconditioner,
    execution::VectorView<const double> b, execution::VectorView<double> x,
    const SolveControl &control) const {
  detail::validate_solver_mpi_context(*state_->mpi_context, "GMRES");
  const detail::SolveCollectiveProtocol protocol(*state_->mpi_context);
  SolveReport report;
  auto finish = [&](SolveTerminationReason reason, int rank = -1) {
    return detail::finish_report(protocol, report, reason, rank);
  };
  auto selected = protocol.agree_control(control);
  if (selected.failed)
    return finish(selected.reason, selected.lowest_rank);

  double root_restart = state_->mpi_context->rank() == 0
                            ? static_cast<double>(state_->restart_length)
                            : 0.0;
  state_->mpi_context->allreduce_fp64_in_place(
      &root_restart, 1U, runtime::Fp64ReductionOperation::sum);
  selected = protocol.checkpoint(
      root_restart == static_cast<double>(state_->restart_length)
          ? detail::LocalSolveFailure::none
          : detail::LocalSolveFailure::collective_failure);
  if (selected.failed)
    return finish(selected.reason, selected.lowest_rank);

  std::size_t owned_count = 0U;
  auto local_failure = detail::LocalSolveFailure::none;
  try {
    owned_count = detail::validate_solver_preflight(
        linear_operator, *state_->execution_context, b, x, "GMRES");
    if (!detail::all_finite(b) || !detail::all_finite(x))
      local_failure = detail::LocalSolveFailure::non_finite_value;
  } catch (const runtime::detail::MpiOperationError &) {
    throw;
  } catch (...) {
    local_failure = detail::LocalSolveFailure::collective_failure;
  }
  selected = protocol.checkpoint(local_failure);
  if (selected.failed)
    return finish(selected.reason, selected.lowest_rank);

  const std::size_t restart = state_->restart_length;
  std::vector<std::optional<execution::Buffer>> replacement_workspace;
  std::vector<double> hessenberg;
  std::vector<double> cosine;
  std::vector<double> sine;
  std::vector<double> least_squares_rhs;
  std::vector<double> coefficients;
  std::vector<double> dots;
  const std::size_t required_workspace_count =
      checked_workspace_count(state_->restart_length);
  const bool replace_workspace =
      state_->workspace_owned_count != owned_count ||
      state_->workspace.size() != required_workspace_count ||
      std::any_of(state_->workspace.begin(), state_->workspace.end(),
                  [](const auto &buffer) { return !buffer.has_value(); });
  try {
    if (replace_workspace) {
      replacement_workspace.resize(required_workspace_count);
      const auto bytes = checked_vector_bytes(owned_count);
      for (auto &buffer : replacement_workspace)
        buffer.emplace(*state_->execution_context, bytes);
    }
    hessenberg.resize((restart + 1U) * restart, 0.0);
    cosine.resize(restart, 0.0);
    sine.resize(restart, 0.0);
    least_squares_rhs.resize(restart + 1U, 0.0);
    coefficients.resize(restart, 0.0);
    dots.resize(restart, 0.0);
  } catch (const runtime::detail::MpiOperationError &) {
    throw;
  } catch (...) {
    local_failure = detail::LocalSolveFailure::collective_failure;
  }
  selected = protocol.checkpoint(local_failure);
  if (selected.failed)
    return finish(selected.reason, selected.lowest_rank);

  if (replace_workspace) {
    state_->workspace = std::move(replacement_workspace);
    state_->workspace_owned_count = owned_count;
  }

  const std::size_t z_offset = restart + 1U;
  const std::size_t r_index = z_offset + restart;
  const std::size_t w_index = r_index + 1U;
  const std::size_t ax_index = w_index + 1U;
  const auto mutable_view = [&](std::size_t index) {
    return state_->workspace[index]->view(0U, owned_count);
  };
  const auto const_view = [&](std::size_t index) {
    return static_cast<const execution::Buffer &>(*state_->workspace[index])
        .view(0U, owned_count);
  };
  const auto h_index = [restart](std::size_t row, std::size_t column) {
    return column * (restart + 1U) + row;
  };
  const auto checkpoint = [&](detail::LocalSolveFailure failure) {
    selected = protocol.checkpoint(failure);
    return !selected.failed;
  };
  const auto validate_local_vector = [&](execution::VectorView<const double> v) {
    return detail::all_finite(v) ? detail::LocalSolveFailure::none
                                 : detail::LocalSolveFailure::non_finite_value;
  };
  const auto apply_operator = [&](execution::VectorView<const double> input,
                                  execution::VectorView<double> output) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      report.matvec_count = checked_increment(report.matvec_count, "matvec");
      linear_operator.apply(input, output).wait();
      failure = validate_local_vector(output);
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return checkpoint(failure);
  };
  const auto update_preconditioner = [&] {
    auto failure = detail::LocalSolveFailure::none;
    try {
      preconditioner.update(linear_operator, linear_operator.revision());
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return checkpoint(failure);
  };
  const auto apply_preconditioner = [&](execution::VectorView<const double> in,
                                        execution::VectorView<double> out) {
    auto failure = detail::LocalSolveFailure::none;
    try {
      report.preconditioner_apply_count = checked_increment(
          report.preconditioner_apply_count, "preconditioner apply");
      preconditioner.apply(in, out).wait();
      failure = validate_local_vector(out);
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (...) {
      failure = detail::LocalSolveFailure::collective_failure;
    }
    return checkpoint(failure);
  };
  const auto norm = [&](execution::VectorView<const double> values,
                        double &result) {
    double local = 0.0;
    auto failure = detail::LocalSolveFailure::none;
    for (std::size_t index = 0U; index < values.size(); ++index) {
      const double value = values[index];
      if (!std::isfinite(value) ||
          !std::isfinite(local + value * value)) {
        failure = detail::LocalSolveFailure::non_finite_value;
        break;
      }
      local += value * value;
    }
    if (!checkpoint(failure))
      return false;
    state_->mpi_context->allreduce_fp64_in_place(
        &local, 1U, runtime::Fp64ReductionOperation::sum);
    result = std::sqrt(local);
    return checkpoint(std::isfinite(result)
                          ? detail::LocalSolveFailure::none
                          : detail::LocalSolveFailure::non_finite_value);
  };
  const auto dot_batch = [&](std::size_t count,
                             execution::VectorView<const double> right) {
    std::fill(dots.data(), dots.data() + count, 0.0);
    auto failure = detail::LocalSolveFailure::none;
    for (std::size_t vector = 0U; vector < count; ++vector) {
      const auto left = const_view(vector);
      for (std::size_t index = 0U; index < right.size(); ++index) {
        const double candidate = dots[vector] + left[index] * right[index];
        if (!std::isfinite(candidate)) {
          failure = detail::LocalSolveFailure::non_finite_value;
          break;
        }
        dots[vector] = candidate;
      }
      if (failure != detail::LocalSolveFailure::none)
        break;
    }
    if (!checkpoint(failure))
      return false;
    state_->mpi_context->allreduce_fp64_in_place(
        dots.data(), count, runtime::Fp64ReductionOperation::sum);
    return checkpoint(std::all_of(dots.data(), dots.data() + count,
                                  [](double value) {
                                    return std::isfinite(value);
                                  })
                          ? detail::LocalSolveFailure::none
                          : detail::LocalSolveFailure::non_finite_value);
  };
  const auto true_residual = [&](double &value) {
    auto ax = mutable_view(ax_index);
    if (!apply_operator(x, ax))
      return false;
    auto r = mutable_view(r_index);
    auto failure = detail::LocalSolveFailure::none;
    for (std::size_t index = 0U; index < owned_count; ++index) {
      const double candidate = b[index] - ax[index];
      if (!std::isfinite(candidate)) {
        failure = detail::LocalSolveFailure::non_finite_value;
        break;
      }
      r[index] = candidate;
    }
    if (!checkpoint(failure))
      return false;
    return norm(const_view(r_index), value);
  };

  double rhs_norm = 0.0;
  if (!norm(b, rhs_norm))
    return finish(selected.reason, selected.lowest_rank);
  const double relative_tolerance = control.rtol * rhs_norm;
  const double tolerance = std::max(control.atol, relative_tolerance);
  if (!std::isfinite(relative_tolerance) || !std::isfinite(tolerance))
    return finish(SolveTerminationReason::non_finite_value);

  double independent = std::numeric_limits<double>::infinity();
  if (!true_residual(independent))
    return finish(selected.reason, selected.lowest_rank);
  report.initial_residual = independent;
  report.recursive_residual = independent;
  report.final_residual = independent;
  if (rhs_norm == 0.0) {
    for (std::size_t index = 0U; index < x.size(); ++index)
      x[index] = 0.0;
    if (!true_residual(independent))
      return finish(selected.reason, selected.lowest_rank);
    report.final_residual = independent;
    return finish(independent <= control.atol
                      ? SolveTerminationReason::zero_right_hand_side
                      : SolveTerminationReason::numerical_breakdown);
  }
  if (independent <= tolerance)
    return finish(SolveTerminationReason::converged);
  if (control.max_iterations == 0U)
    return finish(SolveTerminationReason::maximum_iterations);
  if (!update_preconditioner())
    return finish(selected.reason, selected.lowest_rank);

  while (report.iterations < control.max_iterations) {
    const double beta = independent;
    auto v0 = mutable_view(0U);
    for (std::size_t index = 0U; index < owned_count; ++index)
      v0[index] = const_view(r_index)[index] / beta;
    if (!checkpoint(validate_local_vector(v0)))
      return finish(selected.reason, selected.lowest_rank);
    std::fill(hessenberg.begin(), hessenberg.end(), 0.0);
    std::fill(cosine.begin(), cosine.end(), 0.0);
    std::fill(sine.begin(), sine.end(), 0.0);
    std::fill(least_squares_rhs.begin(), least_squares_rhs.end(), 0.0);
    least_squares_rhs[0] = beta;

    std::size_t cycle_iterations = 0U;
    bool happy_breakdown = false;
    const auto cycle_limit = std::min<std::uint64_t>(
        restart, control.max_iterations - report.iterations);
    for (std::size_t column = 0U; column < cycle_limit; ++column) {
      auto z = mutable_view(z_offset + column);
      if (!apply_preconditioner(const_view(column), z))
        return finish(selected.reason, selected.lowest_rank);
      auto w = mutable_view(w_index);
      if (!apply_operator(const_view(z_offset + column), w))
        return finish(selected.reason, selected.lowest_rank);

      for (std::size_t pass = 0U; pass < 2U; ++pass) {
        if (!dot_batch(column + 1U, const_view(w_index)))
          return finish(selected.reason, selected.lowest_rank);
        auto failure = detail::LocalSolveFailure::none;
        for (std::size_t row = 0U; row <= column; ++row) {
          hessenberg[h_index(row, column)] += dots[row];
          const auto basis = const_view(row);
          for (std::size_t index = 0U; index < owned_count; ++index) {
            const double candidate = w[index] - dots[row] * basis[index];
            if (!std::isfinite(candidate)) {
              failure = detail::LocalSolveFailure::non_finite_value;
              break;
            }
            w[index] = candidate;
          }
          if (failure != detail::LocalSolveFailure::none)
            break;
        }
        if (!checkpoint(failure))
          return finish(selected.reason, selected.lowest_rank);
      }
      double next_norm = 0.0;
      if (!norm(const_view(w_index), next_norm))
        return finish(selected.reason, selected.lowest_rank);
      hessenberg[h_index(column + 1U, column)] = next_norm;
      happy_breakdown = next_norm == 0.0;
      if (!happy_breakdown) {
        auto next_basis = mutable_view(column + 1U);
        for (std::size_t index = 0U; index < owned_count; ++index)
          next_basis[index] = w[index] / next_norm;
        if (!checkpoint(validate_local_vector(next_basis)))
          return finish(selected.reason, selected.lowest_rank);
      }

      for (std::size_t row = 0U; row < column; ++row) {
        const double upper = hessenberg[h_index(row, column)];
        const double lower = hessenberg[h_index(row + 1U, column)];
        hessenberg[h_index(row, column)] =
            cosine[row] * upper + sine[row] * lower;
        hessenberg[h_index(row + 1U, column)] =
            -sine[row] * upper + cosine[row] * lower;
      }
      const double diagonal = hessenberg[h_index(column, column)];
      const double subdiagonal =
          hessenberg[h_index(column + 1U, column)];
      const double magnitude = std::hypot(diagonal, subdiagonal);
      if (!(magnitude > 0.0) || !std::isfinite(magnitude))
        return finish(SolveTerminationReason::numerical_breakdown);
      cosine[column] = diagonal / magnitude;
      sine[column] = subdiagonal / magnitude;
      hessenberg[h_index(column, column)] = magnitude;
      hessenberg[h_index(column + 1U, column)] = 0.0;
      const double rhs_upper = least_squares_rhs[column];
      const double rhs_lower = least_squares_rhs[column + 1U];
      least_squares_rhs[column] =
          cosine[column] * rhs_upper + sine[column] * rhs_lower;
      least_squares_rhs[column + 1U] =
          -sine[column] * rhs_upper + cosine[column] * rhs_lower;
      if (!checkpoint(std::isfinite(least_squares_rhs[column]) &&
                              std::isfinite(least_squares_rhs[column + 1U])
                          ? detail::LocalSolveFailure::none
                          : detail::LocalSolveFailure::non_finite_value))
        return finish(selected.reason, selected.lowest_rank);
      report.iterations = checked_increment(report.iterations, "iteration");
      cycle_iterations = column + 1U;
      report.recursive_residual =
          std::abs(least_squares_rhs[column + 1U]);
      if (report.recursive_residual <= tolerance || happy_breakdown ||
          report.iterations == control.max_iterations)
        break;
    }

    std::fill(coefficients.begin(), coefficients.end(), 0.0);
    bool triangular_ok = cycle_iterations != 0U;
    for (std::size_t reverse = cycle_iterations; reverse-- > 0U;) {
      double value = least_squares_rhs[reverse];
      for (std::size_t column = reverse + 1U; column < cycle_iterations;
           ++column)
        value -= hessenberg[h_index(reverse, column)] * coefficients[column];
      const double diagonal = hessenberg[h_index(reverse, reverse)];
      if (diagonal == 0.0 || !std::isfinite(diagonal) ||
          !std::isfinite(value / diagonal)) {
        triangular_ok = false;
        break;
      }
      coefficients[reverse] = value / diagonal;
    }
    if (!checkpoint(triangular_ok ? detail::LocalSolveFailure::none
                                  : detail::LocalSolveFailure::non_finite_value))
      return finish(selected.reason, selected.lowest_rank);

    auto candidate = mutable_view(w_index);
    auto candidate_failure = detail::LocalSolveFailure::none;
    for (std::size_t index = 0U; index < owned_count; ++index) {
      double value = x[index];
      for (std::size_t column = 0U; column < cycle_iterations; ++column)
        value += coefficients[column] *
                 const_view(z_offset + column)[index];
      if (!std::isfinite(value)) {
        candidate_failure = detail::LocalSolveFailure::non_finite_value;
        break;
      }
      candidate[index] = value;
    }
    if (!checkpoint(candidate_failure))
      return finish(selected.reason, selected.lowest_rank);
    for (std::size_t index = 0U; index < owned_count; ++index)
      x[index] = candidate[index];
    if (!true_residual(independent))
      return finish(selected.reason, selected.lowest_rank);
    report.final_residual = independent;
    if (independent <= tolerance)
      return finish(SolveTerminationReason::converged);
    if (report.iterations == control.max_iterations)
      return finish(SolveTerminationReason::maximum_iterations);
    if (happy_breakdown)
      return finish(SolveTerminationReason::numerical_breakdown);
  }
  return finish(SolveTerminationReason::maximum_iterations);
}

} // namespace hundun::linear
