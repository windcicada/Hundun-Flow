// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/linear/vector_ops.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "execution/src/execution_test_access.hpp"
#include "linear/src/conjugate_gradient_detail.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using hundun::execution::BackendIdentity;
using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::execution::ExecutionCapability;
using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionEvent;
using hundun::execution::ExecutionSpace;
using hundun::execution::VectorView;
using hundun::execution::test::ExecutionTestAccess;
using hundun::linear::ConjugateGradientSolver;
using hundun::linear::IdentityPreconditioner;
using hundun::linear::JacobiPreconditioner;
using hundun::linear::LinearOperator;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::SolveTerminationReason;
using hundun::linear::VectorLayout;
using hundun::linear::VectorOps;
using hundun::runtime::Error;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
namespace allocation_probe = hundun::test::allocation_probe;

constexpr std::size_t kSystemSize = 63U;
constexpr std::size_t kRecordedCalls = 8U;

template <class Function>
std::string expect_error(Function&& function) {
  try {
    function();
  } catch (const Error& error) {
    return error.what();
  }
  throw std::runtime_error("expected hundun::runtime::Error");
}

template <class Function>
void expect_error_containing(Function&& function, const std::string& text) {
  const auto message = expect_error(std::forward<Function>(function));
  HUNDUN_CHECK(message.find(text) != std::string::npos);
}

VectorLayout make_layout(std::size_t count, int rank,
                         std::uint64_t offset = 0U) {
  std::vector<hundun::mesh::GlobalCellId> ids;
  ids.reserve(count);
  const std::uint64_t first =
      offset + static_cast<std::uint64_t>(rank) * UINT64_C(100000);
  for (std::size_t index = 0; index < count; ++index) {
    ids.push_back(first + static_cast<std::uint64_t>(index) + 1U);
  }
  return VectorLayout(count, std::move(ids));
}

void fill(VectorView<double> values, double value) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index] = value;
  }
}

class ConfigurableContext final : public ExecutionContext {
 public:
  ConfigurableContext(BackendIdentity identity, ExecutionSpace space,
                      bool buffer, bool host,
                      bool transfer) noexcept
      : identity_(identity),
        space_(space),
        buffer_(buffer),
        host_(host),
        transfer_(transfer) {}

  std::string_view backend_name() const noexcept override { return "test"; }
  BackendIdentity backend_identity() const noexcept override {
    return identity_;
  }
  ExecutionSpace space() const noexcept override { return space_; }
  bool ordered() const noexcept override { return true; }
  bool supports(ExecutionCapability capability) const noexcept override {
    switch (capability) {
      case ExecutionCapability::buffer_allocation:
        return buffer_;
      case ExecutionCapability::host_access:
        return host_;
      case ExecutionCapability::transfer:
        return transfer_;
      case ExecutionCapability::asynchronous_event:
        return false;
    }
    return false;
  }

 private:
  BackendIdentity identity_;
  ExecutionSpace space_;
  bool buffer_;
  bool host_;
  bool transfer_;
};

class TridiagonalOperator final : public LinearOperator {
 public:
  TridiagonalOperator(ExecutionContext& context, VectorLayout domain,
                      VectorLayout range, double diagonal,
                      double off_diagonal)
      : context_(&context),
        domain_(std::move(domain)),
        range_(std::move(range)),
        diagonal_(diagonal),
        off_diagonal_(off_diagonal),
        completed_(ExecutionEvent::completed()),
        failed_(ExecutionTestAccess::pending_event({})) {
    ExecutionTestAccess::complete_error(failed_, "CG operator event error");
  }

  VectorLayout domain_layout() const override { return domain_; }
  VectorLayout range_layout() const override { return range_; }
  const ExecutionContext& context() const override { return *context_; }
  std::uint64_t revision() const override { return revision_; }

  ExecutionEvent apply(VectorView<const double> input,
                       VectorView<double> output) const override {
    ++apply_calls_;
    if (fail_apply_call_ == apply_calls_) {
      return failed_;
    }
    if (input.size() != domain_.owned_count() ||
        output.size() != range_.owned_count()) {
      throw Error("test operator size mismatch");
    }
    if (apply_calls_ <= kRecordedCalls && input.size() <= kSystemSize) {
      recorded_sizes_[apply_calls_ - 1U] = input.size();
      for (std::size_t index = 0; index < input.size(); ++index) {
        recorded_inputs_[apply_calls_ - 1U][index] = input[index];
      }
    }
    if (tiny_curvature_call_ == apply_calls_) {
      fill(output, 0.0);
      if (output.size() != 0U) {
        output[0] = std::copysign(
            std::numeric_limits<double>::min(), input[0]);
      }
      return completed_;
    }
    for (std::size_t index = 0; index < output.size(); ++index) {
      double value = diagonal_ * input[index];
      if (index != 0U) {
        value += off_diagonal_ * input[index - 1U];
      }
      if (index + 1U != input.size()) {
        value += off_diagonal_ * input[index + 1U];
      }
      output[index] = value + bias_;
    }
    if (allocation_begin_call_ == apply_calls_) {
      allocation_probe::attempt_count = 0U;
      allocation_probe::count_attempts = true;
    }
    return completed_;
  }

  bool has_diagonal() const override { return true; }
  ExecutionEvent diagonal(VectorView<double> output) const override {
    for (std::size_t index = 0; index < output.size(); ++index) {
      output[index] = diagonal_;
    }
    return completed_;
  }

  void fail_apply_on(std::uint64_t call) noexcept { fail_apply_call_ = call; }
  void tiny_curvature_on(std::uint64_t call) noexcept {
    tiny_curvature_call_ = call;
  }
  void set_bias(double bias) noexcept { bias_ = bias; }
  void begin_allocation_observation_on(std::uint64_t begin_call) noexcept {
    allocation_begin_call_ = begin_call;
  }
  std::uint64_t apply_calls() const noexcept { return apply_calls_; }
  std::size_t recorded_size(std::size_t call) const noexcept {
    return recorded_sizes_[call - 1U];
  }
  double recorded(std::size_t call, std::size_t index) const noexcept {
    return recorded_inputs_[call - 1U][index];
  }
  double diagonal_value() const noexcept { return diagonal_; }
  double off_diagonal_value() const noexcept { return off_diagonal_; }

 private:
  ExecutionContext* context_;
  VectorLayout domain_;
  VectorLayout range_;
  double diagonal_;
  double off_diagonal_;
  double bias_{0.0};
  std::uint64_t revision_{1U};
  ExecutionEvent completed_;
  ExecutionEvent failed_;
  mutable std::uint64_t apply_calls_{0U};
  std::uint64_t fail_apply_call_{std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t tiny_curvature_call_{
      std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t allocation_begin_call_{
      std::numeric_limits<std::uint64_t>::max()};
  mutable std::array<std::array<double, kSystemSize>, kRecordedCalls>
      recorded_inputs_{};
  mutable std::array<std::size_t, kRecordedCalls> recorded_sizes_{};
};

class InstrumentedPreconditioner final : public Preconditioner {
 public:
  explicit InstrumentedPreconditioner(double multiplier = 1.0)
      : multiplier_(multiplier),
        completed_(ExecutionEvent::completed()),
        failed_(ExecutionTestAccess::pending_event({})) {
    ExecutionTestAccess::complete_error(failed_,
                                        "CG preconditioner event error");
  }

  void update(const LinearOperator&, std::uint64_t) override {
    ++update_calls_;
  }
  ExecutionEvent apply(VectorView<const double> residual,
                       VectorView<double> correction) const override {
    ++apply_calls_;
    for (std::size_t index = 0; index < residual.size(); ++index) {
      correction[index] = multiplier_ * residual[index];
    }
    return fail_apply_ ? failed_ : completed_;
  }

  void fail_apply(bool enabled) noexcept { fail_apply_ = enabled; }
  std::uint64_t update_calls() const noexcept { return update_calls_; }
  std::uint64_t apply_calls() const noexcept { return apply_calls_; }

 private:
  double multiplier_;
  ExecutionEvent completed_;
  ExecutionEvent failed_;
  std::uint64_t update_calls_{0U};
  mutable std::uint64_t apply_calls_{0U};
  bool fail_apply_{false};
};

void make_known_system(const TridiagonalOperator& linear_operator,
                       VectorView<double> exact,
                       VectorView<double> right_hand_side) {
  for (std::size_t index = 0; index < exact.size(); ++index) {
    exact[index] = 0.5 + std::sin(0.13 * static_cast<double>(index + 1U));
  }
  linear_operator.apply(exact, right_hand_side).wait();
}

double independent_residual_norm(
    const TridiagonalOperator& linear_operator, VectorView<const double> b,
    VectorView<const double> x, CpuReferenceContext& execution,
    const MpiContext& oracle_context) {
  Buffer residual_buffer(execution, b.size() * sizeof(double));
  auto residual = residual_buffer.view(0U, b.size());
  for (std::size_t index = 0; index < b.size(); ++index) {
    double value = linear_operator.diagonal_value() * x[index];
    if (index != 0U) {
      value += linear_operator.off_diagonal_value() * x[index - 1U];
    }
    if (index + 1U != b.size()) {
      value += linear_operator.off_diagonal_value() * x[index + 1U];
    }
    residual[index] = b[index] - value;
  }
  return VectorOps(execution).norm(residual, oracle_context);
}

void check_report_counter_oracles(
    const hundun::linear::SolveReport& report,
    const TridiagonalOperator& linear_operator,
    const InstrumentedPreconditioner& preconditioner,
    std::uint64_t operator_before, std::uint64_t preconditioner_before,
    std::uint64_t reductions_before, const MpiContext& context) {
  HUNDUN_CHECK(report.matvec_count ==
               linear_operator.apply_calls() - operator_before);
  HUNDUN_CHECK(report.preconditioner_apply_count ==
               preconditioner.apply_calls() - preconditioner_before);
  HUNDUN_CHECK(report.global_reduction_count ==
               context.fp64_reduction_counters().collective_calls -
                   reductions_before);
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void check_identical_report(const SolveReport& report,
                            const MpiContext& world) {
  constexpr std::size_t kWords = 9U;
  std::array<std::uint64_t, kWords> local{
      static_cast<std::uint64_t>(report.reason),
      report.iterations,
      bits(report.initial_residual),
      bits(report.recursive_residual),
      bits(report.final_residual),
      report.matvec_count,
      report.preconditioner_apply_count,
      report.global_reduction_count,
      static_cast<std::uint64_t>(
          static_cast<std::int64_t>(report.lowest_failing_rank))};
  std::array<std::array<std::uint64_t, kWords>, 4U> gathered{};
  HUNDUN_CHECK(world.size() <= static_cast<int>(gathered.size()));
  HUNDUN_CHECK(MPI_Allgather(local.data(), static_cast<int>(local.size()),
                             MPI_UINT64_T, gathered.data(),
                             static_cast<int>(local.size()), MPI_UINT64_T,
                             world.comm()) == MPI_SUCCESS);
  for (int rank = 1; rank < world.size(); ++rank) {
    HUNDUN_CHECK(gathered[static_cast<std::size_t>(rank)] == gathered[0]);
  }
}

void test_manufactured_spd(const MpiContext& world,
                           const MpiContext& oracle) {
  CpuReferenceContext execution;
  const auto layout = make_layout(kSystemSize, world.rank());
  TridiagonalOperator linear_operator(execution, layout, layout, 4.0, -1.0);
  Buffer exact_buffer(execution, kSystemSize * sizeof(double));
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto exact = exact_buffer.view(0U, kSystemSize);
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  make_known_system(linear_operator, exact, b);
  fill(x, 0.0);

  InstrumentedPreconditioner preconditioner;
  ConjugateGradientSolver solver(execution, world);
  const auto op_before = linear_operator.apply_calls();
  const auto preconditioner_before = preconditioner.apply_calls();
  const auto reductions_before =
      world.fp64_reduction_counters().collective_calls;
  const auto report = solver.solve(linear_operator, preconditioner, b, x, {});

  HUNDUN_CHECK(report.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(report.iterations > 0U);
  HUNDUN_CHECK(std::isfinite(report.initial_residual));
  HUNDUN_CHECK(std::isfinite(report.recursive_residual));
  HUNDUN_CHECK(std::isfinite(report.final_residual));
  HUNDUN_CHECK(report.final_residual <= report.initial_residual);
  HUNDUN_CHECK(preconditioner.update_calls() == 1U);
  check_report_counter_oracles(report, linear_operator, preconditioner,
                               op_before, preconditioner_before,
                               reductions_before, world);
  const double measured =
      independent_residual_norm(linear_operator, b, x, execution, oracle);
  const double oracle_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, measured);
  HUNDUN_CHECK(std::abs(report.final_residual - measured) <=
               oracle_tolerance);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    HUNDUN_CHECK(std::abs(x[index] - exact[index]) <= 2.0e-10);
  }
}

void test_accepted_preconditioners(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(kSystemSize, world.rank(), 2000000U);
  TridiagonalOperator linear_operator(execution, layout, layout, 4.0, -1.0);
  Buffer exact_buffer(execution, kSystemSize * sizeof(double));
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto exact = exact_buffer.view(0U, kSystemSize);
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  make_known_system(linear_operator, exact, b);
  ConjugateGradientSolver solver(execution, world);

  IdentityPreconditioner identity(execution);
  fill(x, 0.0);
  HUNDUN_CHECK(solver.solve(linear_operator, identity, b, x, {}).reason ==
               SolveTerminationReason::converged);

  JacobiPreconditioner jacobi(execution);
  fill(x, 0.0);
  HUNDUN_CHECK(solver.solve(linear_operator, jacobi, b, x, {}).reason ==
               SolveTerminationReason::converged);
}

void test_early_and_limit_exits(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(5U, world.rank(), 3000000U);
  TridiagonalOperator linear_operator(execution, layout, layout, 2.0, 0.0);
  ConjugateGradientSolver solver(execution, world);
  Buffer b_buffer(execution, 5U * sizeof(double));
  Buffer x_buffer(execution, 5U * sizeof(double));
  auto b = b_buffer.view(0U, 5U);
  auto x = x_buffer.view(0U, 5U);
  InstrumentedPreconditioner preconditioner;

  fill(b, 2.0);
  fill(x, 1.0);
  const auto exact = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(exact.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(exact.iterations == 0U);
  HUNDUN_CHECK(exact.initial_residual == 0.0);
  HUNDUN_CHECK(exact.recursive_residual == 0.0);
  HUNDUN_CHECK(exact.final_residual == 0.0);
  HUNDUN_CHECK(preconditioner.update_calls() == 0U);
  HUNDUN_CHECK(preconditioner.apply_calls() == 0U);

  fill(b, 0.0);
  fill(x, 7.0);
  const auto zero = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(zero.reason == SolveTerminationReason::zero_right_hand_side);
  HUNDUN_CHECK(zero.iterations == 0U);
  HUNDUN_CHECK(zero.matvec_count == 2U);
  HUNDUN_CHECK(preconditioner.update_calls() == 0U);
  HUNDUN_CHECK(preconditioner.apply_calls() == 0U);
  HUNDUN_CHECK(zero.initial_residual > 0.0);
  HUNDUN_CHECK(zero.recursive_residual == zero.initial_residual);
  HUNDUN_CHECK(zero.final_residual == 0.0);
  for (std::size_t index = 0; index < x.size(); ++index) {
    HUNDUN_CHECK(x[index] == 0.0);
  }

  TridiagonalOperator nonlinear_zero(execution, layout, layout, 2.0, 0.0);
  nonlinear_zero.set_bias(1.0);
  fill(b, 0.0);
  fill(x, 3.0);
  const auto rejected_zero =
      solver.solve(nonlinear_zero, preconditioner, b, x, {});
  HUNDUN_CHECK(rejected_zero.reason ==
               SolveTerminationReason::numerical_breakdown);
  HUNDUN_CHECK(rejected_zero.iterations == 0U);
  HUNDUN_CHECK(preconditioner.update_calls() == 0U);
  HUNDUN_CHECK(preconditioner.apply_calls() == 0U);
  HUNDUN_CHECK(rejected_zero.recursive_residual ==
               rejected_zero.initial_residual);
  HUNDUN_CHECK(rejected_zero.final_residual > 0.0);

  fill(b, 1.0);
  fill(x, 0.0);
  SolveControl no_iterations;
  no_iterations.max_iterations = 0U;
  const auto none =
      solver.solve(linear_operator, preconditioner, b, x, no_iterations);
  HUNDUN_CHECK(none.reason == SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(none.iterations == 0U);
  HUNDUN_CHECK(none.matvec_count == 1U);
  HUNDUN_CHECK(preconditioner.update_calls() == 0U);

  TridiagonalOperator coupled(execution, layout, layout, 4.0, -1.0);
  fill(b, 1.0);
  fill(x, 0.0);
  SolveControl one_iteration;
  one_iteration.atol = 0.0;
  one_iteration.rtol = 0.0;
  one_iteration.max_iterations = 1U;
  const auto limited =
      solver.solve(coupled, preconditioner, b, x, one_iteration);
  HUNDUN_CHECK(limited.reason == SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(limited.iterations == 1U);
  HUNDUN_CHECK(limited.matvec_count == 3U);
  HUNDUN_CHECK(std::isfinite(limited.final_residual));
}

void test_dynamic_spd_checks(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(4U, world.rank(), 4000000U);
  Buffer b_buffer(execution, 4U * sizeof(double));
  Buffer x_buffer(execution, 4U * sizeof(double));
  auto b = b_buffer.view(0U, 4U);
  auto x = x_buffer.view(0U, 4U);
  fill(b, 1.0);
  ConjugateGradientSolver solver(execution, world);

  for (double diagonal : {-1.0, 0.0}) {
    TridiagonalOperator linear_operator(execution, layout, layout, diagonal,
                                        0.0);
    InstrumentedPreconditioner preconditioner;
    fill(x, 0.0);
    const auto report =
        solver.solve(linear_operator, preconditioner, b, x, {});
    HUNDUN_CHECK(report.reason ==
                 SolveTerminationReason::numerical_breakdown);
    HUNDUN_CHECK(report.iterations == 0U);
  }

  TridiagonalOperator positive(execution, layout, layout, 1.0, 0.0);
  for (double multiplier : {-1.0, 0.0}) {
    InstrumentedPreconditioner preconditioner(multiplier);
    fill(x, 0.0);
    const auto report = solver.solve(positive, preconditioner, b, x, {});
    HUNDUN_CHECK(report.reason ==
                 SolveTerminationReason::numerical_breakdown);
    HUNDUN_CHECK(report.iterations == 0U);
    HUNDUN_CHECK(report.matvec_count == 1U);
  }
}

void test_nonfinite_curvature_refreshes_independent_residual(
    const MpiContext& world, const MpiContext& oracle) {
  CpuReferenceContext execution;
  const auto layout = make_layout(4U, world.rank(), 4500000U);
  TridiagonalOperator linear_operator(execution, layout, layout, 4.0, -1.0);
  Buffer b_buffer(execution, 4U * sizeof(double));
  Buffer x_buffer(execution, 4U * sizeof(double));
  auto b = b_buffer.view(0U, 4U);
  auto x = x_buffer.view(0U, 4U);
  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 10.0 * static_cast<double>(index + 1U);
  }
  fill(x, 0.0);
  linear_operator.tiny_curvature_on(3U);

  InstrumentedPreconditioner preconditioner;
  ConjugateGradientSolver solver(execution, world);
  SolveControl control;
  control.atol = 0.0;
  control.rtol = 0.0;
  control.max_iterations = 3U;

  const auto report =
      solver.solve(linear_operator, preconditioner, b, x, control);

  HUNDUN_CHECK(report.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(report.iterations == 1U);
  HUNDUN_CHECK(report.matvec_count == 4U);
  const double measured =
      independent_residual_norm(linear_operator, b, x, execution, oracle);
  const double oracle_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, measured);
  HUNDUN_CHECK(std::abs(report.final_residual - measured) <=
               oracle_tolerance);
}

void test_periodic_replacement_restart_and_allocations(
    const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(kSystemSize, world.rank(), 5000000U);
  TridiagonalOperator linear_operator(execution, layout, layout, 4.0, -1.0);
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    b[index] = index % 2U == 0U ? 1.0 : 2.0;
  }
  fill(x, 0.0);
  InstrumentedPreconditioner preconditioner;
  ConjugateGradientSolver solver(execution, world);
  SolveControl control;
  control.atol = 0.0;
  control.rtol = 0.0;
  control.max_iterations = 2U;
  control.residual_recompute_interval = 1U;
  linear_operator.begin_allocation_observation_on(2U);

  const auto reductions_before =
      world.fp64_reduction_counters().collective_calls;
  const auto report =
      solver.solve(linear_operator, preconditioner, b, x, control);
  const std::size_t iteration_allocation_attempts =
      allocation_probe::attempt_count;
  allocation_probe::count_attempts = false;
  HUNDUN_CHECK(report.reason == SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(report.iterations == 2U);
  HUNDUN_CHECK(report.matvec_count == 5U);
  HUNDUN_CHECK(report.preconditioner_apply_count == 2U);
  HUNDUN_CHECK(report.global_reduction_count == 54U);
  HUNDUN_CHECK(world.fp64_reduction_counters().collective_calls -
                   reductions_before ==
               54U);
  HUNDUN_CHECK(iteration_allocation_attempts == 0U);
  HUNDUN_CHECK(linear_operator.recorded_size(4U) == kSystemSize);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    double applied = 4.0 * linear_operator.recorded(3U, index);
    if (index != 0U) {
      applied -= linear_operator.recorded(3U, index - 1U);
    }
    if (index + 1U != kSystemSize) {
      applied -= linear_operator.recorded(3U, index + 1U);
    }
    const double expected_restart = b[index] - applied;
    HUNDUN_CHECK(std::abs(linear_operator.recorded(4U, index) -
                          expected_restart) <=
                 8.0 * std::numeric_limits<double>::epsilon());
  }
}

void test_control_preflight_and_events(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(3U, world.rank(), 6000000U);
  TridiagonalOperator linear_operator(execution, layout, layout, 2.0, 0.0);
  Buffer b_buffer(execution, 6U * sizeof(double));
  Buffer x_buffer(execution, 6U * sizeof(double));
  auto b = b_buffer.view(0U, 3U, 2U);
  auto x = x_buffer.view(0U, 3U, 2U);
  fill(b, 2.0);
  fill(x, 0.0);
  InstrumentedPreconditioner preconditioner;
  ConjugateGradientSolver solver(execution, world);

  const auto strided = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(strided.reason == SolveTerminationReason::converged);
  for (std::size_t index = 0; index < x.size(); ++index) {
    HUNDUN_CHECK_NEAR(x[index], 1.0, 1.0e-14);
  }

  const auto reductions_before =
      world.fp64_reduction_counters().collective_calls;
  const auto calls_before = linear_operator.apply_calls();
  const auto updates_before = preconditioner.update_calls();
  fill(x, 4.0);
  SolveControl invalid;
  invalid.residual_recompute_interval = 0U;
  const auto invalid_report =
      solver.solve(linear_operator, preconditioner, b, x, invalid);
  HUNDUN_CHECK(invalid_report.reason ==
               SolveTerminationReason::invalid_control);
  HUNDUN_CHECK(invalid_report.lowest_failing_rank == 0);
  HUNDUN_CHECK(invalid_report.global_reduction_count == 1U);
  HUNDUN_CHECK(linear_operator.apply_calls() == calls_before);
  HUNDUN_CHECK(preconditioner.update_calls() == updates_before);
  HUNDUN_CHECK(world.fp64_reduction_counters().collective_calls ==
               reductions_before + 1U);
  for (std::size_t index = 0; index < x.size(); ++index) {
    HUNDUN_CHECK(x[index] == 4.0);
  }

  invalid = {};
  invalid.atol = std::numeric_limits<double>::quiet_NaN();
  HUNDUN_CHECK(
      solver.solve(linear_operator, preconditioner, b, x, invalid).reason ==
      SolveTerminationReason::invalid_control);
  invalid = {};
  invalid.rtol = -1.0;
  HUNDUN_CHECK(
      solver.solve(linear_operator, preconditioner, b, x, invalid).reason ==
      SolveTerminationReason::invalid_control);

  invalid = {};
  invalid.atol = -1.0;
  ExecutionTestAccess::fail_next_allocation();
  HUNDUN_CHECK(
      solver.solve(linear_operator, preconditioner, b, x, invalid).reason ==
      SolveTerminationReason::invalid_control);
  expect_error_containing(
      [&] { Buffer allocation_would_have_failed(execution, sizeof(double)); },
      "allocation failed");

  Buffer shared(execution, 6U * sizeof(double));
  auto shared_x = shared.view(0U, 3U, 2U);
  fill(shared_x, 1.0);
  const auto shared_b = static_cast<VectorView<const double>>(shared_x);
  const auto before_shared = world.fp64_reduction_counters().collective_calls;
  const auto shared_report =
      solver.solve(linear_operator, preconditioner, shared_b, shared_x, {});
  HUNDUN_CHECK(shared_report.reason ==
               SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(shared_report.lowest_failing_rank == 0);
  HUNDUN_CHECK(world.fp64_reduction_counters().collective_calls ==
               before_shared + 4U);

  Buffer short_buffer(execution, 2U * sizeof(double));
  HUNDUN_CHECK(solver
                   .solve(linear_operator, preconditioner,
                          short_buffer.view(0U, 2U), x, {})
                   .reason == SolveTerminationReason::collective_failure);

  auto nonwritable_metadata = ExecutionTestAccess::metadata(x);
  nonwritable_metadata.writable = false;
  const auto nonwritable_x =
      ExecutionTestAccess::mutable_view(x_buffer, nonwritable_metadata);
  HUNDUN_CHECK(
      solver.solve(linear_operator, preconditioner, b, nonwritable_x, {})
          .reason == SolveTerminationReason::collective_failure);

  auto backend_metadata = ExecutionTestAccess::metadata(x);
  ++backend_metadata.backend_identity;
  const auto wrong_backend_x =
      ExecutionTestAccess::mutable_view(x_buffer, backend_metadata);
  HUNDUN_CHECK(
      solver.solve(linear_operator, preconditioner, b, wrong_backend_x, {})
          .reason == SolveTerminationReason::collective_failure);

  const auto other_range = make_layout(3U, world.rank(), 7000000U);
  TridiagonalOperator nonsquare(execution, layout, other_range, 1.0, 0.0);
  HUNDUN_CHECK(solver.solve(nonsquare, preconditioner, b, x, {}).reason ==
               SolveTerminationReason::collective_failure);

  CpuReferenceContext other_execution;
  TridiagonalOperator wrong_context(other_execution, layout, layout, 1.0, 0.0);
  HUNDUN_CHECK(solver.solve(wrong_context, preconditioner, b, x, {}).reason ==
               SolveTerminationReason::collective_failure);

  Buffer stale_buffer(execution, 3U * sizeof(double));
  auto stale = stale_buffer.view(0U, 3U);
  stale_buffer.reallocate(3U * sizeof(double));
  HUNDUN_CHECK(
      solver.solve(linear_operator, preconditioner, stale, x, {}).reason ==
      SolveTerminationReason::collective_failure);

  TridiagonalOperator failed_operator(execution, layout, layout, 1.0, 0.0);
  failed_operator.fail_apply_on(1U);
  HUNDUN_CHECK(solver.solve(failed_operator, preconditioner, b, x, {}).reason ==
               SolveTerminationReason::collective_failure);

  TridiagonalOperator identity(execution, layout, layout, 1.0, 0.0);
  InstrumentedPreconditioner failed_preconditioner;
  failed_preconditioner.fail_apply(true);
  fill(x, 0.0);
  HUNDUN_CHECK(solver.solve(identity, failed_preconditioner, b, x, {}).reason ==
               SolveTerminationReason::collective_failure);

  SolveControl overflow_tolerance;
  overflow_tolerance.rtol = std::numeric_limits<double>::max();
  fill(b, 2.0);
  fill(x, 0.0);
  const auto nonfinite =
      solver.solve(identity, preconditioner, b, x, overflow_tolerance);
  HUNDUN_CHECK(nonfinite.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(nonfinite.lowest_failing_rank == -1);
  HUNDUN_CHECK(nonfinite.matvec_count == 0U);
  HUNDUN_CHECK(nonfinite.global_reduction_count == 8U);
  check_identical_report(nonfinite, world);
}

void test_zero_owned_and_workspace_overflow(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(0U, world.rank(), 8000000U);
  TridiagonalOperator linear_operator(execution, layout, layout, 1.0, 0.0);
  Buffer b_buffer(execution, 0U);
  Buffer x_buffer(execution, 0U);
  InstrumentedPreconditioner preconditioner;
  ConjugateGradientSolver solver(execution, world);
  const auto report = solver.solve(linear_operator, preconditioner,
                                   b_buffer.view(0U, 0U),
                                   x_buffer.view(0U, 0U), {});
  HUNDUN_CHECK(report.reason == SolveTerminationReason::zero_right_hand_side);
  HUNDUN_CHECK(report.matvec_count == 2U);
  HUNDUN_CHECK(preconditioner.update_calls() == 0U);

  const std::size_t too_large =
      std::numeric_limits<std::size_t>::max() / 5U + 1U;
  expect_error_containing(
      [&] {
        static_cast<void>(
            hundun::linear::detail::checked_cg_workspace_bytes(too_large));
      },
      "overflow");
}

void test_constructor_capabilities(const MpiContext& world) {
  ConfigurableContext zero(0U, ExecutionSpace::host, true, true, true);
  expect_error_containing(
      [&] { ConjugateGradientSolver solver(zero, world); }, "identity");
  ConfigurableContext device(9U, ExecutionSpace::device, true, true, true);
  expect_error_containing(
      [&] { ConjugateGradientSolver solver(device, world); },
      "host-accessible");
  ConfigurableContext no_buffer(9U, ExecutionSpace::host, false, true, true);
  expect_error_containing(
      [&] { ConjugateGradientSolver solver(no_buffer, world); }, "buffer");
  ConfigurableContext no_transfer(9U, ExecutionSpace::host, true, true, false);
  expect_error_containing(
      [&] { ConjugateGradientSolver solver(no_transfer, world); },
      "transfer");

  auto moved_source = MpiContext::duplicate(world.comm());
  auto live_destination = std::move(moved_source);
  CpuReferenceContext execution;
  expect_error_containing(
      [&] { ConjugateGradientSolver solver(execution, moved_source); },
      "MPI context");
  HUNDUN_CHECK(live_destination.size() == world.size());
}

}  // namespace

int main(int argc, char** argv) {
  std::optional<MpiContext> finalized_context;
  int numerical_result = EXIT_SUCCESS;
  {
    MpiEnvironment environment(argc, argv);
    auto world = MpiContext::duplicate(MPI_COMM_WORLD);
    auto oracle = MpiContext::duplicate(MPI_COMM_WORLD);
    numerical_result = hundun::test::run([&] {
      test_manufactured_spd(world, oracle);
      test_accepted_preconditioners(world);
      test_early_and_limit_exits(world);
      test_dynamic_spd_checks(world);
      test_nonfinite_curvature_refreshes_independent_residual(world, oracle);
      test_periodic_replacement_restart_and_allocations(world);
      test_control_preflight_and_events(world);
      test_zero_owned_and_workspace_overflow(world);
      test_constructor_capabilities(world);
    });
    finalized_context.emplace(MpiContext::duplicate(MPI_COMM_WORLD));
  }
  if (numerical_result != EXIT_SUCCESS) {
    return numerical_result;
  }

  return hundun::test::run([&] {
    CpuReferenceContext execution;
    expect_error_containing(
        [&] {
          ConjugateGradientSolver solver(execution, *finalized_context);
        },
        "MPI context");
  });
}
