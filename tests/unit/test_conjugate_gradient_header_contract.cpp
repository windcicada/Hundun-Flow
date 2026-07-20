// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/conjugate_gradient.hpp"

#include <memory>
#include <type_traits>

namespace {

using hundun::execution::ExecutionContext;
using hundun::execution::VectorView;
using hundun::linear::ConjugateGradientSolver;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::runtime::MpiContext;

using ConcreteSolve = SolveReport (ConjugateGradientSolver::*)(
    const LinearOperator&, Preconditioner&, VectorView<const double>,
    VectorView<double>, const SolveControl&) const;

static_assert(std::is_final_v<ConjugateGradientSolver>);
static_assert(std::is_base_of_v<LinearSolver, ConjugateGradientSolver>);
static_assert(std::has_virtual_destructor_v<ConjugateGradientSolver>);
static_assert(std::is_nothrow_destructible_v<ConjugateGradientSolver>);
static_assert(std::is_constructible_v<ConjugateGradientSolver,
                                      ExecutionContext&, const MpiContext&>);
static_assert(!std::is_copy_constructible_v<ConjugateGradientSolver>);
static_assert(!std::is_copy_assignable_v<ConjugateGradientSolver>);
static_assert(std::is_same_v<decltype(static_cast<ConcreteSolve>(
                                 &ConjugateGradientSolver::solve)),
                             ConcreteSolve>);

}  // namespace

int main() { return 0; }
