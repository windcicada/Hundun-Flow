// SPDX-License-Identifier: Apache-2.0

#include "hundun/lin_system.hpp"
#include "hundun/lin_preconditioners.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace {

using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionEvent;
using hundun::execution::VectorView;
using hundun::linear::IdentityPreconditioner;
using hundun::linear::JacobiPreconditioner;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::SolveTerminationReason;
using hundun::linear::VectorLayout;

using OperatorApply = ExecutionEvent (LinearOperator::*)(
    VectorView<const double>, VectorView<double>) const;
using OperatorDiagonal =
    ExecutionEvent (LinearOperator::*)(VectorView<double>) const;
using PreconditionerUpdate =
    void (Preconditioner::*)(const LinearOperator&, std::uint64_t);
using PreconditionerApply = ExecutionEvent (Preconditioner::*)(
    VectorView<const double>, VectorView<double>) const;
using SolverSolve = SolveReport (LinearSolver::*)(
    const LinearOperator&, Preconditioner&, VectorView<const double>,
    VectorView<double>, const SolveControl&) const;

static_assert(std::is_same_v<decltype(static_cast<OperatorApply>(
                                 &LinearOperator::apply)),
                             OperatorApply>);
static_assert(std::is_same_v<decltype(static_cast<OperatorDiagonal>(
                                 &LinearOperator::diagonal)),
                             OperatorDiagonal>);
static_assert(std::is_same_v<decltype(static_cast<PreconditionerUpdate>(
                                 &Preconditioner::update)),
                             PreconditionerUpdate>);
static_assert(std::is_same_v<decltype(static_cast<PreconditionerApply>(
                                 &Preconditioner::apply)),
                             PreconditionerApply>);
static_assert(std::is_same_v<decltype(static_cast<SolverSolve>(
                                 &LinearSolver::solve)),
                             SolverSolve>);

static_assert(std::has_virtual_destructor_v<LinearOperator>);
static_assert(std::has_virtual_destructor_v<Preconditioner>);
static_assert(std::has_virtual_destructor_v<LinearSolver>);
static_assert(std::is_nothrow_destructible_v<LinearOperator>);
static_assert(std::is_nothrow_destructible_v<Preconditioner>);
static_assert(std::is_nothrow_destructible_v<LinearSolver>);
static_assert(std::is_abstract_v<LinearOperator>);
static_assert(std::is_abstract_v<Preconditioner>);
static_assert(std::is_abstract_v<LinearSolver>);

static_assert(std::is_final_v<SolveControl>);
static_assert(std::is_final_v<SolveReport>);
static_assert(std::is_final_v<IdentityPreconditioner>);
static_assert(std::is_final_v<JacobiPreconditioner>);
static_assert(std::is_base_of_v<Preconditioner, IdentityPreconditioner>);
static_assert(std::is_base_of_v<Preconditioner, JacobiPreconditioner>);
static_assert(std::is_constructible_v<IdentityPreconditioner,
                                      ExecutionContext&>);
static_assert(std::is_constructible_v<JacobiPreconditioner,
                                      ExecutionContext&>);
static_assert(!std::is_convertible_v<ExecutionContext&,
                                     IdentityPreconditioner>);
static_assert(!std::is_convertible_v<ExecutionContext&,
                                     JacobiPreconditioner>);
static_assert(!std::is_copy_constructible_v<IdentityPreconditioner>);
static_assert(!std::is_copy_assignable_v<IdentityPreconditioner>);
static_assert(!std::is_copy_constructible_v<JacobiPreconditioner>);
static_assert(!std::is_copy_assignable_v<JacobiPreconditioner>);

static_assert(std::is_same_v<decltype(&LinearOperator::domain_layout),
                             VectorLayout (LinearOperator::*)() const>);
static_assert(std::is_same_v<decltype(&LinearOperator::range_layout),
                             VectorLayout (LinearOperator::*)() const>);
static_assert(std::is_same_v<decltype(&LinearOperator::context),
                             const ExecutionContext&
                                 (LinearOperator::*)() const>);
static_assert(std::is_same_v<decltype(&LinearOperator::revision),
                             std::uint64_t (LinearOperator::*)() const>);
static_assert(std::is_same_v<decltype(&LinearOperator::has_diagonal),
                             bool (LinearOperator::*)() const>);
static_assert(static_cast<int>(SolveTerminationReason::converged) == 0);
static_assert(static_cast<int>(
                  SolveTerminationReason::zero_right_hand_side) == 1);
static_assert(static_cast<int>(SolveTerminationReason::maximum_iterations) ==
              2);
static_assert(static_cast<int>(SolveTerminationReason::numerical_breakdown) ==
              3);
static_assert(static_cast<int>(SolveTerminationReason::non_finite_value) == 4);
static_assert(static_cast<int>(SolveTerminationReason::invalid_control) == 5);
static_assert(static_cast<int>(SolveTerminationReason::collective_failure) ==
              6);

}  // namespace

int main() {
  const SolveControl control;
  if (control.atol != 1.0e-12 || control.rtol != 1.0e-10 ||
      control.max_iterations != 500U ||
      control.residual_recompute_interval != 20U) {
    return 1;
  }

  const SolveReport report;
  if (report.reason != SolveTerminationReason::invalid_control ||
      report.iterations != 0U ||
      !std::isinf(report.initial_residual) ||
      !std::isinf(report.recursive_residual) ||
      !std::isinf(report.final_residual) || report.matvec_count != 0U ||
      report.preconditioner_apply_count != 0U ||
      report.global_reduction_count != 0U ||
      report.lowest_failing_rank != -1) {
    return 1;
  }

  return static_cast<int>(SolveTerminationReason::collective_failure) != 6;
}
