// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/exec_execution.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/lin_restarted_gmres.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/exec_execution_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace {

using namespace hundun;

class NonnormalOperator final : public linear::LinearOperator {
public:
  NonnormalOperator(execution::ExecutionContext &context,
                    linear::VectorLayout layout)
      : context_(&context), layout_(std::move(layout)) {}

  linear::VectorLayout domain_layout() const override { return layout_; }
  linear::VectorLayout range_layout() const override { return layout_; }
  const execution::ExecutionContext &context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return 1U; }
  execution::ExecutionEvent
  apply(execution::VectorView<const double> input,
        execution::VectorView<double> output) const override {
    output[0] = 4.0 * input[0] - 8.0 * input[1];
    output[1] = input[0] + 4.0 * input[1] - 8.0 * input[2];
    output[2] = input[1] + 4.0 * input[2];
    return execution::ExecutionEvent::completed();
  }
  bool has_diagonal() const override { return true; }
  execution::ExecutionEvent
  diagonal(execution::VectorView<double> output) const override {
    for (std::size_t index = 0U; index < output.size(); ++index)
      output[index] = 4.0;
    return execution::ExecutionEvent::completed();
  }

private:
  execution::ExecutionContext *context_{};
  linear::VectorLayout layout_;
};

linear::VectorLayout layout(int rank) {
  const auto first = static_cast<std::uint64_t>(rank) * 3U;
  return linear::VectorLayout(
      3U, std::vector<mesh::GlobalCellId>{first, first + 1U, first + 2U});
}

void run(const runtime::MpiContext &mpi) {
  execution::CpuReferenceContext execution;
  NonnormalOperator matrix(execution, layout(mpi.rank()));
  linear::JacobiPreconditioner preconditioner(execution);
  linear::RestartedGmresSolver solver(execution, mpi, 3U);
  execution::Buffer right_hand_side(execution, 3U * sizeof(double));
  execution::Buffer solution(execution, 3U * sizeof(double));
  const std::array<double, 3> exact{1.0, 2.0, -1.0};
  auto rhs = right_hand_side.view(0U, 3U);
  rhs[0] = -12.0;
  rhs[1] = 17.0;
  rhs[2] = -2.0;
  auto x = solution.view(0U, 3U);
  for (std::size_t index = 0U; index < x.size(); ++index)
    x[index] = 0.0;
  linear::SolveControl control;
  control.atol = 1.0e-13;
  control.rtol = 1.0e-12;
  control.max_iterations = 12U;
  const auto report = solver.solve(
      matrix, preconditioner,
      static_cast<const execution::Buffer &>(right_hand_side).view(0U, 3U),
      x, control);
  HUNDUN_CHECK(report.reason == linear::SolveTerminationReason::converged);
  HUNDUN_CHECK(report.iterations <= 3U);
  HUNDUN_CHECK(report.final_residual <=
               std::max(control.atol, control.rtol * report.initial_residual));
  for (std::size_t index = 0U; index < x.size(); ++index)
    HUNDUN_CHECK_NEAR(x[index], exact[index], 1.0e-11);

  for (std::size_t index = 0U; index < x.size(); ++index)
    x[index] = 0.0;
  const auto allocations_before_repeat = execution::allocation_counters();
  const auto repeated = solver.solve(
      matrix, preconditioner,
      static_cast<const execution::Buffer &>(right_hand_side).view(0U, 3U), x,
      control);
  const auto allocations_after_repeat = execution::allocation_counters();
  HUNDUN_CHECK(repeated.reason == linear::SolveTerminationReason::converged);
  HUNDUN_CHECK(allocations_after_repeat.allocation_events ==
               allocations_before_repeat.allocation_events);
  HUNDUN_CHECK(allocations_after_repeat.deallocation_events ==
               allocations_before_repeat.deallocation_events);
  HUNDUN_CHECK(allocations_after_repeat.allocated_bytes ==
               allocations_before_repeat.allocated_bytes);
  HUNDUN_CHECK(allocations_after_repeat.deallocated_bytes ==
               allocations_before_repeat.deallocated_bytes);
  HUNDUN_CHECK(allocations_after_repeat.live_bytes ==
               allocations_before_repeat.live_bytes);
  HUNDUN_CHECK(allocations_after_repeat.peak_live_bytes ==
               allocations_before_repeat.peak_live_bytes);

  constexpr std::size_t replacement_size = 4U;
  const auto replacement_first =
      static_cast<std::uint64_t>(mpi.rank()) * replacement_size + 1000U;
  linear::VectorLayout replacement_layout(
      replacement_size, std::vector<mesh::GlobalCellId>{
                            replacement_first, replacement_first + 1U,
                            replacement_first + 2U, replacement_first + 3U});
  NonnormalOperator replacement_matrix(execution,
                                       std::move(replacement_layout));
  execution::Buffer replacement_rhs(execution,
                                    replacement_size * sizeof(double));
  execution::Buffer replacement_solution(execution,
                                         replacement_size * sizeof(double));
  auto replacement_rhs_values = replacement_rhs.view(0U, replacement_size);
  auto replacement_x = replacement_solution.view(0U, replacement_size);
  for (std::size_t index = 0U; index < replacement_size; ++index) {
    replacement_rhs_values[index] = 1.0;
    replacement_x[index] = 7.0;
  }
  execution::test::ExecutionTestAccess::fail_next_allocation();
  const auto failed_resize =
      solver.solve(replacement_matrix, preconditioner,
                   static_cast<const execution::Buffer &>(replacement_rhs)
                       .view(0U, replacement_size),
                   replacement_x, control);
  HUNDUN_CHECK(failed_resize.reason ==
               linear::SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(failed_resize.lowest_failing_rank == 0);
  for (std::size_t index = 0U; index < replacement_size; ++index)
    HUNDUN_CHECK(replacement_x[index] == 7.0);

  for (std::size_t index = 0U; index < x.size(); ++index)
    x[index] = 0.0;
  const auto allocations_before_retry = execution::allocation_counters();
  const auto retry = solver.solve(
      matrix, preconditioner,
      static_cast<const execution::Buffer &>(right_hand_side).view(0U, 3U), x,
      control);
  const auto allocations_after_retry = execution::allocation_counters();
  HUNDUN_CHECK(retry.reason == linear::SolveTerminationReason::converged);
  HUNDUN_CHECK(allocations_after_retry.allocation_events ==
               allocations_before_retry.allocation_events);
  HUNDUN_CHECK(allocations_after_retry.deallocation_events ==
               allocations_before_retry.deallocation_events);
  HUNDUN_CHECK(allocations_after_retry.allocated_bytes ==
               allocations_before_retry.allocated_bytes);
  HUNDUN_CHECK(allocations_after_retry.deallocated_bytes ==
               allocations_before_retry.deallocated_bytes);

  linear::RestartedGmresSolver short_restart_solver(execution, mpi, 2U);
  auto short_restart_control = control;
  short_restart_control.max_iterations = 12U;
  for (std::size_t index = 0U; index < x.size(); ++index)
    x[index] = 0.0;
  const auto restarted = short_restart_solver.solve(
      matrix, preconditioner,
      static_cast<const execution::Buffer &>(right_hand_side).view(0U, 3U),
      x, short_restart_control);
  HUNDUN_CHECK(restarted.reason ==
               linear::SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(restarted.iterations > 2U);
  HUNDUN_CHECK(restarted.iterations == short_restart_control.max_iterations);
  HUNDUN_CHECK(restarted.final_residual < restarted.initial_residual);
  HUNDUN_CHECK_NEAR(restarted.recursive_residual,
                    restarted.final_residual, 1.0e-12);

  for (std::size_t index = 0U; index < x.size(); ++index) {
    rhs[index] = 0.0;
    x[index] = 7.0;
  }
  const auto zero = solver.solve(
      matrix, preconditioner,
      static_cast<const execution::Buffer &>(right_hand_side).view(0U, 3U),
      x, control);
  HUNDUN_CHECK(zero.reason ==
               linear::SolveTerminationReason::zero_right_hand_side);
  HUNDUN_CHECK(zero.final_residual == 0.0);
  for (std::size_t index = 0U; index < x.size(); ++index)
    HUNDUN_CHECK(x[index] == 0.0);

  bool rejected = false;
  try {
    linear::RestartedGmresSolver invalid(execution, mpi, 0U);
    static_cast<void>(invalid);
  } catch (const runtime::Error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

} // namespace

int main(int argc, char **argv) {
  runtime::MpiEnvironment environment(argc, argv);
  const auto mpi = runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return test::run([&] { run(mpi); });
}
