// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include <type_traits>

#include "hundun/lin_bicgstab.hpp"

namespace {

using hundun::execution::ExecutionContext;
using hundun::execution::VectorView;
using hundun::linear::BiCGStabSolver;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::runtime::MpiContext;

using ConcreteSolve =
    SolveReport (BiCGStabSolver::*)(const LinearOperator&,
                                    Preconditioner&,
                                    VectorView<const double>,
                                    VectorView<double>,
                                    const SolveControl&) const;

static_assert(std::is_final_v<BiCGStabSolver>);
static_assert(std::is_base_of_v<LinearSolver, BiCGStabSolver>);
static_assert(std::has_virtual_destructor_v<BiCGStabSolver>);
static_assert(std::is_nothrow_destructible_v<BiCGStabSolver>);
static_assert(std::is_constructible_v<BiCGStabSolver,
                                      ExecutionContext&,
                                      const MpiContext&>);
static_assert(!std::is_copy_constructible_v<BiCGStabSolver>);
static_assert(!std::is_copy_assignable_v<BiCGStabSolver>);
static_assert(
    std::is_same_v<decltype(static_cast<ConcreteSolve>(&BiCGStabSolver::solve)),
                   ConcreteSolve>);

}  // namespace

int main() { return 0; }
