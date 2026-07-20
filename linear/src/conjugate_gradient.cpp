// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/conjugate_gradient.hpp"

#include "hundun/linear/vector_ops.hpp"
#include "hundun/runtime/error.hpp"
#include "linear/src/conjugate_gradient_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <string>

namespace hundun::linear {
namespace {

constexpr std::size_t kWorkspaceVectorCount = 5U;

void validate_execution_context(
    const execution::ExecutionContext& context) {
  if (context.backend_identity() == 0U) {
    throw runtime::Error("CG execution context identity must be nonzero");
  }
  if (context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access)) {
    throw runtime::Error("CG requires a host-accessible execution context");
  }
  if (!context.supports(
          execution::ExecutionCapability::buffer_allocation)) {
    throw runtime::Error("CG execution context lacks buffer capability");
  }
  if (!context.supports(execution::ExecutionCapability::transfer)) {
    throw runtime::Error("CG execution context lacks transfer capability");
  }
}

void validate_mpi_context(const runtime::MpiContext& context) {
  int initialized = 0;
  int finalized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS ||
      MPI_Finalized(&finalized) != MPI_SUCCESS) {
    throw runtime::Error("CG could not query MPI context activity");
  }
  if (initialized == 0 || finalized != 0) {
    throw runtime::Error("CG requires a live MPI context");
  }
  if (context.comm() == MPI_COMM_NULL || context.size() <= 0 ||
      context.rank() < 0 || context.rank() >= context.size()) {
    throw runtime::Error("CG requires a live nonempty MPI context");
  }
}

bool valid_control(const SolveControl& control) noexcept {
  return std::isfinite(control.atol) && control.atol >= 0.0 &&
         std::isfinite(control.rtol) && control.rtol >= 0.0 &&
         control.residual_recompute_interval != 0U;
}

template <class T>
void validate_view(execution::VectorView<T> view,
                   const execution::ExecutionContext& context,
                   std::size_t expected_size, bool require_writable,
                   const char* name) {
  if (view.scalar_format() != execution::ScalarFormat::float64) {
    throw runtime::Error(std::string("CG ") + name +
                         " view must be float64");
  }
  if (view.size() != expected_size) {
    throw runtime::Error(std::string("CG ") + name +
                         " size does not match the owned layout");
  }
  if (view.backend_identity() != context.backend_identity()) {
    throw runtime::Error(std::string("CG ") + name +
                         " backend identity does not match");
  }
  if (view.space() != context.space() ||
      view.space() != execution::ExecutionSpace::host) {
    throw runtime::Error(std::string("CG ") + name +
                         " execution space does not match host context");
  }
  if (require_writable && !view.writable()) {
    throw runtime::Error("CG solution view must be writable");
  }
  static_cast<void>(view.data());
}

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

bool finite_axpy_update(double alpha,
                        execution::VectorView<const double> input,
                        execution::VectorView<const double> current) {
  const double* input_pointer = input.data();
  const double* current_pointer = current.data();
  for (std::size_t index = 0; index < input.size(); ++index) {
    const double result =
        alpha * input_pointer[index * input.stride()] +
        current_pointer[index * current.stride()];
    if (!std::isfinite(result)) {
      return false;
    }
  }
  return true;
}

bool finite_linear_combination(
    double alpha, execution::VectorView<const double> left, double beta,
    execution::VectorView<const double> right) {
  const double* left_pointer = left.data();
  const double* right_pointer = right.data();
  for (std::size_t index = 0; index < left.size(); ++index) {
    const double result =
        alpha * left_pointer[index * left.stride()] +
        beta * right_pointer[index * right.stride()];
    if (!std::isfinite(result)) {
      return false;
    }
  }
  return true;
}

std::size_t validate_preflight(
    const LinearOperator& linear_operator,
    const execution::ExecutionContext& context,
    execution::VectorView<const double> b,
    execution::VectorView<double> x) {
  if (&linear_operator.context() != &context) {
    throw runtime::Error("CG operator context object does not match");
  }
  const VectorLayout domain = linear_operator.domain_layout();
  const VectorLayout range = linear_operator.range_layout();
  if (domain != range) {
    throw runtime::Error("CG requires an exactly square operator layout");
  }
  validate_view(b, context, domain.owned_count(), false, "right-hand-side");
  validate_view(x, context, domain.owned_count(), true, "solution");
  if (b.allocation_identity() == x.allocation_identity()) {
    throw runtime::Error(
        "CG right-hand-side and solution must not share an allocation");
  }
  return domain.owned_count();
}

}  // namespace

namespace detail {

std::size_t checked_cg_workspace_bytes(std::size_t owned_count) {
  if (owned_count >
      std::numeric_limits<std::size_t>::max() / kWorkspaceVectorCount) {
    throw runtime::Error("CG workspace vector count overflow");
  }
  const std::size_t elements = owned_count * kWorkspaceVectorCount;
  if (elements >
      std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw runtime::Error("CG workspace byte size overflow");
  }
  return elements * sizeof(double);
}

}  // namespace detail

namespace {

std::uint64_t checked_increment(std::uint64_t value, const char* name) {
  if (value == std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error(std::string("CG ") + name +
                         " counter would overflow");
  }
  return value + 1U;
}

}  // namespace

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
  validate_execution_context(execution_context);
  validate_mpi_context(mpi_context);
  try {
    state_ = std::make_unique<State>(execution_context, mpi_context);
  } catch (const std::bad_alloc&) {
    throw runtime::Error("CG solver state allocation failed");
  }
}

ConjugateGradientSolver::~ConjugateGradientSolver() noexcept = default;

SolveReport ConjugateGradientSolver::solve(
    const LinearOperator& linear_operator, Preconditioner& preconditioner,
    execution::VectorView<const double> b,
    execution::VectorView<double> x, const SolveControl& control) const {
  SolveReport report;
  if (!valid_control(control)) {
    return report;
  }

  validate_mpi_context(*state_->mpi_context);
  const std::size_t owned_count = validate_preflight(
      linear_operator, *state_->execution_context, b, x);
  const std::size_t workspace_bytes =
      detail::checked_cg_workspace_bytes(owned_count);
  const std::size_t vector_bytes = workspace_bytes / kWorkspaceVectorCount;
  execution::Buffer r_buffer(*state_->execution_context, vector_bytes);
  execution::Buffer z_buffer(*state_->execution_context, vector_bytes);
  execution::Buffer p_buffer(*state_->execution_context, vector_bytes);
  execution::Buffer ap_buffer(*state_->execution_context, vector_bytes);
  execution::Buffer ax_buffer(*state_->execution_context, vector_bytes);
  execution::Buffer scalar(*state_->execution_context, sizeof(double));
  auto r = r_buffer.view(0U, owned_count);
  auto z = z_buffer.view(0U, owned_count);
  auto p = p_buffer.view(0U, owned_count);
  auto ap = ap_buffer.view(0U, owned_count);
  auto ax = ax_buffer.view(0U, owned_count);
  auto reduction = scalar.view(0U, 1U);

  const std::uint64_t reductions_before =
      state_->mpi_context->fp64_reduction_counters().collective_calls;
  auto finish = [&](SolveTerminationReason reason) {
    const std::uint64_t reductions_after =
        state_->mpi_context->fp64_reduction_counters().collective_calls;
    if (reductions_after < reductions_before) {
      throw runtime::Error("CG global reduction counter decreased");
    }
    report.reason = reason;
    report.global_reduction_count = reductions_after - reductions_before;
    return report;
  };

  auto apply_operator = [&](execution::VectorView<const double> input,
                            execution::VectorView<double> output) {
    report.matvec_count = checked_increment(report.matvec_count, "matvec");
    auto event = linear_operator.apply(input, output);
    event.wait();
  };
  auto apply_preconditioner = [&] {
    report.preconditioner_apply_count = checked_increment(
        report.preconditioner_apply_count, "preconditioner apply");
    auto event = preconditioner.apply(r, z);
    event.wait();
    return all_finite(z);
  };
  auto dot = [&](execution::VectorView<const double> left,
                 execution::VectorView<const double> right) {
    const DotProductPair pair{left, right};
    state_->operations.dot_batch(&pair, 1U, reduction,
                                 *state_->mpi_context);
    return reduction[0];
  };
  auto independent_residual = [&] {
    apply_operator(x, ax);
    if (!all_finite(ax) || !finite_linear_combination(1.0, b, -1.0, ax)) {
      return std::numeric_limits<double>::infinity();
    }
    state_->operations.linear_combination(1.0, b, -1.0, ax, r);
    return state_->operations.norm(r, *state_->mpi_context);
  };

  const double right_hand_side_norm =
      state_->operations.norm(b, *state_->mpi_context);
  const double relative_tolerance = control.rtol * right_hand_side_norm;
  const double tolerance = std::max(control.atol, relative_tolerance);
  if (!std::isfinite(relative_tolerance) || !std::isfinite(tolerance)) {
    return finish(SolveTerminationReason::non_finite_value);
  }

  double independent = independent_residual();
  report.initial_residual = independent;
  report.recursive_residual = independent;
  report.final_residual = independent;
  bool independent_is_current = true;
  if (!std::isfinite(independent)) {
    return finish(SolveTerminationReason::non_finite_value);
  }

  if (right_hand_side_norm == 0.0) {
    state_->operations.fill(x, 0.0);
    independent_is_current = false;
    independent = independent_residual();
    independent_is_current = true;
    report.final_residual = independent;
    if (!std::isfinite(independent)) {
      return finish(SolveTerminationReason::non_finite_value);
    }
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

  preconditioner.update(linear_operator, linear_operator.revision());
  if (!apply_preconditioner()) {
    return finish(SolveTerminationReason::non_finite_value);
  }
  double rho = dot(r, z);
  if (!std::isfinite(rho)) {
    return finish(SolveTerminationReason::non_finite_value);
  }
  if (rho <= 0.0) {
    return finish(SolveTerminationReason::numerical_breakdown);
  }
  state_->operations.copy(z, p);

  auto refresh_final_residual = [&] {
    if (!independent_is_current) {
      independent = independent_residual();
      independent_is_current = true;
      report.final_residual = independent;
    }
    return std::isfinite(report.final_residual);
  };

  while (report.iterations < control.max_iterations) {
    apply_operator(p, ap);
    if (!all_finite(ap)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    const double curvature = dot(p, ap);
    if (!std::isfinite(curvature)) {
      return finish(SolveTerminationReason::non_finite_value);
    }
    if (curvature <= 0.0) {
      if (!refresh_final_residual()) {
        return finish(SolveTerminationReason::non_finite_value);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }

    const double alpha = rho / curvature;
    if (!std::isfinite(alpha)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    if (!finite_axpy_update(alpha, p, x)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    state_->operations.axpy(alpha, p, x);
    independent_is_current = false;
    if (!finite_axpy_update(-alpha, ap, r)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    state_->operations.axpy(-alpha, ap, r);
    report.iterations = checked_increment(report.iterations, "iteration");
    report.recursive_residual =
        state_->operations.norm(r, *state_->mpi_context);

    const bool periodic =
        report.iterations % control.residual_recompute_interval == 0U;
    const bool candidate = report.recursive_residual <= tolerance;
    const bool exhausted = report.iterations == control.max_iterations;
    if (periodic || candidate || exhausted) {
      independent = independent_residual();
      independent_is_current = true;
      report.final_residual = independent;
      if (!std::isfinite(independent)) {
        return finish(SolveTerminationReason::non_finite_value);
      }
      if (independent <= tolerance) {
        return finish(SolveTerminationReason::converged);
      }
      if (exhausted) {
        return finish(SolveTerminationReason::maximum_iterations);
      }

      report.recursive_residual = independent;
      if (!apply_preconditioner()) {
        return finish(SolveTerminationReason::non_finite_value);
      }
      rho = dot(r, z);
      if (!std::isfinite(rho)) {
        return finish(SolveTerminationReason::non_finite_value);
      }
      if (rho <= 0.0) {
        return finish(SolveTerminationReason::numerical_breakdown);
      }
      state_->operations.copy(z, p);
      continue;
    }

    if (!apply_preconditioner()) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    const double rho_new = dot(r, z);
    if (!std::isfinite(rho_new)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    if (rho_new <= 0.0) {
      if (!refresh_final_residual()) {
        return finish(SolveTerminationReason::non_finite_value);
      }
      return finish(SolveTerminationReason::numerical_breakdown);
    }
    const double beta = rho_new / rho;
    if (!std::isfinite(beta)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    if (!finite_linear_combination(1.0, z, beta, p)) {
      static_cast<void>(refresh_final_residual());
      return finish(SolveTerminationReason::non_finite_value);
    }
    state_->operations.linear_combination(1.0, z, beta, p, p);
    rho = rho_new;
  }

  if (!independent_is_current) {
    report.final_residual = independent_residual();
  }
  return finish(SolveTerminationReason::maximum_iterations);
}

}  // namespace hundun::linear
