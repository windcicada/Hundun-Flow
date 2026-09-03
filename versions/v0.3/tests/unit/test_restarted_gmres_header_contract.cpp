// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include <cstdint>
#include <type_traits>

#include "hundun/lin_restarted_gmres.hpp"

namespace {

using hundun::execution::ExecutionContext;
using hundun::execution::VectorView;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::RestartedFgmresSolver;
using hundun::linear::RestartedGmresSolver;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::runtime::MpiContext;

using ConcreteSolve =
    SolveReport (RestartedGmresSolver::*)(const LinearOperator &,
                                          Preconditioner &,
                                          VectorView<const double>,
                                          VectorView<double>,
                                          const SolveControl &) const;

static_assert(std::is_final_v<RestartedGmresSolver>);
static_assert(std::is_base_of_v<LinearSolver, RestartedGmresSolver>);
static_assert(std::is_same_v<RestartedFgmresSolver, RestartedGmresSolver>);
static_assert(std::has_virtual_destructor_v<RestartedGmresSolver>);
static_assert(std::is_nothrow_destructible_v<RestartedGmresSolver>);
static_assert(std::is_constructible_v<RestartedGmresSolver,
                                      ExecutionContext &, const MpiContext &,
                                      std::uint32_t>);
static_assert(!std::is_copy_constructible_v<RestartedGmresSolver>);
static_assert(!std::is_copy_assignable_v<RestartedGmresSolver>);
static_assert(std::is_same_v<
              decltype(static_cast<ConcreteSolve>(&RestartedGmresSolver::solve)),
              ConcreteSolve>);

} // namespace

int main() { return 0; }
