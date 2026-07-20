// SPDX-License-Identifier: Apache-2.0

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "execution/src/execution_test_access.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/linear/bicgstab.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "linear/src/bicgstab_detail.hpp"
#include "runtime/src/mpi_context_test_seam.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

namespace {

using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionEvent;
using hundun::execution::VectorView;
using hundun::execution::test::ExecutionTestAccess;
using hundun::linear::BiCGStabSolver;
using hundun::linear::ConjugateGradientSolver;
using hundun::linear::IdentityPreconditioner;
using hundun::linear::JacobiPreconditioner;
using hundun::linear::LinearOperator;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::SolveTerminationReason;
using hundun::linear::VectorLayout;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
namespace allocation_probe = hundun::test::allocation_probe;

constexpr std::size_t kSystemSize = 63U;

template <class Function>
void expect_runtime_error(Function&& function) {
  try {
    function();
  } catch (const hundun::runtime::Error&) {
    return;
  }
  throw std::runtime_error("expected hundun::runtime::Error");
}

VectorLayout make_layout(std::size_t count, int rank) {
  std::vector<hundun::mesh::GlobalCellId> ids;
  ids.reserve(count);
  const auto first = static_cast<std::uint64_t>(rank) * UINT64_C(100000);
  for (std::size_t index = 0; index < count; ++index) {
    ids.push_back(first + static_cast<std::uint64_t>(index) + 1U);
  }
  return VectorLayout(count, std::move(ids));
}

struct CandidateOverflowState final {
  int failure_rank{-1};
  std::uint64_t operator_call{std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t preconditioner_call{
      std::numeric_limits<std::uint64_t>::max()};
  std::optional<VectorView<double>> solution;
  std::optional<VectorView<double>> candidate;
};

class NonsymmetricOperator final : public LinearOperator {
 public:
  NonsymmetricOperator(ExecutionContext& context,
                       VectorLayout layout,
                       int rank = 0)
      : context_(&context),
        layout_(std::move(layout)),
        range_(layout_),
        rank_(rank),
        completed_(ExecutionEvent::completed()),
        failed_(ExecutionTestAccess::pending_event({})) {
    ExecutionTestAccess::complete_error(failed_,
                                        "BiCGStab operator event failure");
  }

  NonsymmetricOperator(ExecutionContext& context,
                       VectorLayout domain,
                       VectorLayout range,
                       int rank)
      : context_(&context),
        layout_(std::move(domain)),
        range_(std::move(range)),
        rank_(rank),
        completed_(ExecutionEvent::completed()),
        failed_(ExecutionTestAccess::pending_event({})) {
    ExecutionTestAccess::complete_error(failed_,
                                        "BiCGStab operator event failure");
  }

  VectorLayout domain_layout() const override { return layout_; }
  VectorLayout range_layout() const override { return range_; }
  const ExecutionContext& context() const override { return *context_; }
  std::uint64_t revision() const override { return 1U; }
  ExecutionEvent apply(VectorView<const double> input,
                       VectorView<double> output) const override {
    ++apply_calls_;
    if (apply_calls_ <= recorded_inputs_.size() &&
        input.size() <= kSystemSize) {
      recorded_sizes_[apply_calls_ - 1U] = input.size();
      for (std::size_t index = 0; index < input.size(); ++index) {
        recorded_inputs_[apply_calls_ - 1U][index] = input[index];
      }
    }
    if (rank_ == event_failure_rank_ && apply_calls_ == failure_call_) {
      return failed_;
    }
    if (apply_calls_ == zero_output_call_) {
      for (std::size_t index = 0; index < output.size(); ++index) {
        output[index] = 0.0;
      }
      return completed_;
    }
    if (apply_calls_ == orthogonal_output_call_) {
      for (std::size_t index = 0; index < output.size(); ++index) {
        output[index] = 0.0;
      }
      if (output.size() >= 2U) {
        output[0] = -input[1];
        output[1] = input[0];
      }
      return completed_;
    }
    if (apply_calls_ == tiny_output_call_) {
      for (std::size_t index = 0; index < output.size(); ++index) {
        output[index] = 0.0;
      }
      if (output.size() != 0U) {
        output[0] = std::copysign(std::numeric_limits<double>::min(), input[0]);
      }
      return completed_;
    }
    if (candidate_overflow_ != nullptr &&
        rank_ == candidate_overflow_->failure_rank &&
        candidate_overflow_->solution.has_value() &&
        input.allocation_identity() ==
            candidate_overflow_->solution->allocation_identity() &&
        (*candidate_overflow_->solution)[0] ==
            std::numeric_limits<double>::max()) {
      for (std::size_t index = 0; index < output.size(); ++index) {
        output[index] = 0.0;
      }
      return completed_;
    }
    for (std::size_t index = 0; index < output.size(); ++index) {
      double value = 4.0 * input[index];
      if (index != 0U) {
        value -= input[index - 1U];
      }
      if (index + 1U != output.size()) {
        value -= 2.0 * input[index + 1U];
      }
      output[index] = value + bias_;
    }
    if (candidate_overflow_ != nullptr &&
        rank_ == candidate_overflow_->failure_rank &&
        apply_calls_ == candidate_overflow_->operator_call &&
        output.size() != 0U) {
      HUNDUN_CHECK(candidate_overflow_->candidate.has_value());
      double local_ts = 0.0;
      for (std::size_t index = 0; index < output.size(); ++index) {
        local_ts += output[index] * input[index];
      }
      HUNDUN_CHECK(local_ts != 0.0);
      (*candidate_overflow_->solution)[0] =
          std::numeric_limits<double>::max();
      (*candidate_overflow_->candidate)[0] = std::copysign(
          std::numeric_limits<double>::max(), local_ts);
    }
    if (rank_ == nonfinite_output_rank_ && apply_calls_ == failure_call_ &&
        output.size() != 0U) {
      output[0] = std::numeric_limits<double>::quiet_NaN();
    }
    return completed_;
  }
  bool has_diagonal() const override { return true; }
  ExecutionEvent diagonal(VectorView<double> output) const override {
    for (std::size_t index = 0; index < output.size(); ++index) {
      output[index] = 4.0;
    }
    return completed_;
  }

  void fail_event_on(int rank, std::uint64_t call) noexcept {
    event_failure_rank_ = rank;
    failure_call_ = call;
  }
  void nonfinite_output_on(int rank, std::uint64_t call) noexcept {
    nonfinite_output_rank_ = rank;
    failure_call_ = call;
  }
  void zero_output_on(std::uint64_t call) noexcept { zero_output_call_ = call; }
  void orthogonal_output_on(std::uint64_t call) noexcept {
    orthogonal_output_call_ = call;
  }
  void tiny_output_on(std::uint64_t call) noexcept { tiny_output_call_ = call; }
  void set_bias(double bias) noexcept { bias_ = bias; }
  void overflow_solution_candidate_after_apply_on(
      int rank, std::uint64_t call, VectorView<double> solution,
      CandidateOverflowState& state) noexcept {
    state.failure_rank = rank;
    state.operator_call = call;
    state.solution = solution;
    candidate_overflow_ = &state;
  }
  std::uint64_t apply_calls() const noexcept { return apply_calls_; }
  double recorded(std::size_t call, std::size_t index) const noexcept {
    return recorded_inputs_[call - 1U][index];
  }
  std::size_t recorded_size(std::size_t call) const noexcept {
    return recorded_sizes_[call - 1U];
  }

 private:
  ExecutionContext* context_;
  VectorLayout layout_;
  VectorLayout range_;
  int rank_;
  ExecutionEvent completed_;
  ExecutionEvent failed_;
  mutable std::uint64_t apply_calls_{0U};
  int event_failure_rank_{-1};
  int nonfinite_output_rank_{-1};
  std::uint64_t failure_call_{std::numeric_limits<std::uint64_t>::max()};
  double bias_{0.0};
  std::uint64_t zero_output_call_{std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t orthogonal_output_call_{
      std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t tiny_output_call_{std::numeric_limits<std::uint64_t>::max()};
  CandidateOverflowState* candidate_overflow_{nullptr};
  mutable std::array<std::array<double, kSystemSize>, 8U> recorded_inputs_{};
  mutable std::array<std::size_t, 8U> recorded_sizes_{};
};

class InstrumentedPreconditioner final : public Preconditioner {
 public:
  explicit InstrumentedPreconditioner(int rank)
      : rank_(rank),
        completed_(ExecutionEvent::completed()),
        failed_(ExecutionTestAccess::pending_event({})) {
    ExecutionTestAccess::complete_error(
        failed_, "BiCGStab preconditioner event failure");
  }

  void update(const LinearOperator&, std::uint64_t) override {
    ++update_calls_;
    if (rank_ == update_failure_rank_) {
      throw hundun::runtime::Error("injected preconditioner update failure");
    }
    if (arm_allocation_probe_) {
      allocation_probe::attempt_count = 0U;
      allocation_probe::count_attempts = true;
    }
  }

  ExecutionEvent apply(VectorView<const double> input,
                       VectorView<double> output) const override {
    ++apply_calls_;
    for (std::size_t index = 0; index < output.size(); ++index) {
      output[index] = input[index];
    }
    if (candidate_overflow_ != nullptr &&
        rank_ == candidate_overflow_->failure_rank &&
        apply_calls_ == candidate_overflow_->preconditioner_call) {
      candidate_overflow_->candidate = output;
    }
    if (rank_ == nonfinite_output_rank_ && apply_calls_ == failure_call_ &&
        output.size() != 0U) {
      output[0] = std::numeric_limits<double>::infinity();
    }
    if (rank_ == event_failure_rank_ && apply_calls_ == failure_call_) {
      return failed_;
    }
    return completed_;
  }

  void fail_update_on(int rank) noexcept { update_failure_rank_ = rank; }
  void fail_event_on(int rank, std::uint64_t call) noexcept {
    event_failure_rank_ = rank;
    failure_call_ = call;
  }
  void nonfinite_output_on(int rank, std::uint64_t call) noexcept {
    nonfinite_output_rank_ = rank;
    failure_call_ = call;
  }
  void arm_allocation_probe_after_update() noexcept {
    arm_allocation_probe_ = true;
  }
  void capture_solution_candidate_on(std::uint64_t call,
                                     CandidateOverflowState& state) noexcept {
    state.preconditioner_call = call;
    candidate_overflow_ = &state;
  }
  std::uint64_t apply_calls() const noexcept { return apply_calls_; }
  std::uint64_t update_calls() const noexcept { return update_calls_; }

 private:
  int rank_;
  ExecutionEvent completed_;
  ExecutionEvent failed_;
  std::uint64_t update_calls_{0U};
  mutable std::uint64_t apply_calls_{0U};
  int update_failure_rank_{-1};
  int event_failure_rank_{-1};
  int nonfinite_output_rank_{-1};
  std::uint64_t failure_call_{std::numeric_limits<std::uint64_t>::max()};
  bool arm_allocation_probe_{false};
  CandidateOverflowState* candidate_overflow_{nullptr};
};

enum class ScriptedRecurrence { rho_zero, beta_overflow };

class ScriptedRecurrenceOperator final : public LinearOperator {
 public:
  ScriptedRecurrenceOperator(ExecutionContext& context,
                             VectorLayout layout,
                             ScriptedRecurrence script)
      : context_(&context),
        layout_(std::move(layout)),
        script_(script),
        completed_(ExecutionEvent::completed()) {}

  VectorLayout domain_layout() const override { return layout_; }
  VectorLayout range_layout() const override { return layout_; }
  const ExecutionContext& context() const override { return *context_; }
  std::uint64_t revision() const override { return 1U; }
  ExecutionEvent apply(VectorView<const double>,
                       VectorView<double> output) const override {
    ++apply_calls_;
    for (std::size_t index = 0; index < output.size(); ++index) {
      output[index] = 0.0;
    }
    if (apply_calls_ == 2U) {
      output[0] = 1.0;
    } else if (apply_calls_ == 3U) {
      if (script_ == ScriptedRecurrence::rho_zero) {
        output[0] = -1.0;
        output[1] = 1.0;
        output[2] = 1.0;
      } else {
        output[0] = 1.0;
        output[1] = 1.0;
        output[2] = 1.0e-154;
      }
    }
    return completed_;
  }
  bool has_diagonal() const override { return false; }
  ExecutionEvent diagonal(VectorView<double>) const override {
    throw hundun::runtime::Error("scripted operator has no diagonal");
  }
  std::uint64_t apply_calls() const noexcept { return apply_calls_; }

 private:
  ExecutionContext* context_;
  VectorLayout layout_;
  ScriptedRecurrence script_;
  ExecutionEvent completed_;
  mutable std::uint64_t apply_calls_{0U};
};

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

double global_l2(double local_sum, const MpiContext& world) {
  double global_sum = local_sum;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &global_sum, 1, MPI_DOUBLE, MPI_SUM,
                             world.comm()) == MPI_SUCCESS);
  return std::sqrt(global_sum);
}

void test_manufactured(const MpiContext& world) {
  CpuReferenceContext execution;
  NonsymmetricOperator linear_operator(execution,
                                       make_layout(kSystemSize, world.rank()));
  Buffer exact_buffer(execution, kSystemSize * sizeof(double));
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto exact = exact_buffer.view(0U, kSystemSize);
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    exact[index] = 0.5 + std::sin(0.13 * static_cast<double>(index + 1U));
    x[index] = 0.0;
  }
  linear_operator.apply(exact, b).wait();

  IdentityPreconditioner preconditioner(execution);
  BiCGStabSolver solver(execution, world);
  const auto report = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(report.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(report.lowest_failing_rank == -1);
  HUNDUN_CHECK(report.iterations > 0U);
  HUNDUN_CHECK(std::isfinite(report.final_residual));
  double local_b_squared = 0.0;
  double local_residual_squared = 0.0;
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    HUNDUN_CHECK_NEAR(x[index], exact[index], 2.0e-9);
    double applied = 4.0 * x[index];
    if (index != 0U) {
      applied -= x[index - 1U];
    }
    if (index + 1U != kSystemSize) {
      applied -= 2.0 * x[index + 1U];
    }
    const double residual = b[index] - applied;
    local_b_squared += b[index] * b[index];
    local_residual_squared += residual * residual;
  }
  const double b_norm = global_l2(local_b_squared, world);
  const double residual_norm = global_l2(local_residual_squared, world);
  HUNDUN_CHECK(report.final_residual <= std::max(1.0e-12, 1.0e-10 * b_norm));
  HUNDUN_CHECK(std::abs(report.final_residual - residual_norm) <=
               64.0 * std::numeric_limits<double>::epsilon() *
                   std::max(1.0, residual_norm));

  for (std::size_t index = 0; index < kSystemSize; ++index) {
    x[index] = 0.0;
  }
  JacobiPreconditioner jacobi(execution);
  const auto jacobi_report = solver.solve(linear_operator, jacobi, b, x, {});
  HUNDUN_CHECK(jacobi_report.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(jacobi_report.lowest_failing_rank == -1);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    HUNDUN_CHECK_NEAR(x[index], exact[index], 2.0e-9);
  }
}

void test_collective_control_and_input(const MpiContext& world) {
  CpuReferenceContext execution;
  NonsymmetricOperator linear_operator(execution,
                                       make_layout(3U, world.rank()));
  Buffer b_buffer(execution, 3U * sizeof(double));
  Buffer x_buffer(execution, 3U * sizeof(double));
  auto b = b_buffer.view(0U, 3U);
  auto x = x_buffer.view(0U, 3U);
  for (std::size_t index = 0; index < 3U; ++index) {
    b[index] = 1.0;
    x[index] = 7.0;
  }
  IdentityPreconditioner preconditioner(execution);
  BiCGStabSolver solver(execution, world);

  if (world.size() > 1) {
    SolveControl invalid;
    if (world.rank() == 1) {
      invalid.residual_recompute_interval = 0U;
    }
    const auto invalid_report =
        solver.solve(linear_operator, preconditioner, b, x, invalid);
    HUNDUN_CHECK(invalid_report.reason ==
                 SolveTerminationReason::invalid_control);
    HUNDUN_CHECK(invalid_report.lowest_failing_rank == 1);
    check_identical_report(invalid_report, world);
    for (std::size_t index = 0; index < 3U; ++index) {
      HUNDUN_CHECK(x[index] == 7.0);
    }

    SolveControl different;
    if (world.rank() == 1) {
      different.max_iterations = 7U;
    }
    const auto different_report =
        solver.solve(linear_operator, preconditioner, b, x, different);
    HUNDUN_CHECK(different_report.reason ==
                 SolveTerminationReason::collective_failure);
    HUNDUN_CHECK(different_report.lowest_failing_rank == 1);
    check_identical_report(different_report, world);

    b[0] = world.rank() == 0 ? std::numeric_limits<double>::quiet_NaN() : 1.0;
    const auto nonfinite_report =
        solver.solve(linear_operator, preconditioner, b, x, {});
    HUNDUN_CHECK(nonfinite_report.reason ==
                 SolveTerminationReason::non_finite_value);
    HUNDUN_CHECK(nonfinite_report.lowest_failing_rank == 0);
    check_identical_report(nonfinite_report, world);

    b[0] = 1.0;
    x[0] = world.rank() == 1 ? std::numeric_limits<double>::infinity() : 7.0;
    const auto nonfinite_x =
        solver.solve(linear_operator, preconditioner, b, x, {});
    HUNDUN_CHECK(nonfinite_x.reason ==
                 SolveTerminationReason::non_finite_value);
    HUNDUN_CHECK(nonfinite_x.lowest_failing_rank == 1);
    check_identical_report(nonfinite_x, world);
  }
}

void test_tolerance_overflow_is_global(const MpiContext& world) {
  CpuReferenceContext execution;
  NonsymmetricOperator linear_operator(execution,
                                       make_layout(3U, world.rank()));
  Buffer b_buffer(execution, 3U * sizeof(double));
  Buffer x_buffer(execution, 3U * sizeof(double));
  auto b = b_buffer.view(0U, 3U);
  auto x = x_buffer.view(0U, 3U);
  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 2.0;
    x[index] = 0.0;
  }
  SolveControl overflow;
  overflow.rtol = std::numeric_limits<double>::max();
  InstrumentedPreconditioner preconditioner(world.rank());
  BiCGStabSolver solver(execution, world);
  const auto report =
      solver.solve(linear_operator, preconditioner, b, x, overflow);
  HUNDUN_CHECK(report.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(report.lowest_failing_rank == -1);
  HUNDUN_CHECK(report.matvec_count == 0U);
  HUNDUN_CHECK(report.preconditioner_apply_count == 0U);
  HUNDUN_CHECK(report.global_reduction_count == 8U);
  check_identical_report(report, world);
}

void check_global_norm_overflow_report(const SolveReport& report,
                                       const MpiContext& world) {
  HUNDUN_CHECK(report.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(report.lowest_failing_rank == -1);
  HUNDUN_CHECK(report.iterations == 0U);
  HUNDUN_CHECK(report.matvec_count == 0U);
  HUNDUN_CHECK(report.preconditioner_apply_count == 0U);
  HUNDUN_CHECK(report.global_reduction_count == 8U);
  HUNDUN_CHECK(std::isinf(report.initial_residual));
  HUNDUN_CHECK(std::isinf(report.recursive_residual));
  HUNDUN_CHECK(std::isinf(report.final_residual));
  check_identical_report(report, world);
}

void check_local_dot_overflow_report(const SolveReport& report,
                                     const MpiContext& world,
                                     int expected_lowest_rank,
                                     std::uint64_t expected_preconditioners) {
  HUNDUN_CHECK(report.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(report.lowest_failing_rank == expected_lowest_rank);
  HUNDUN_CHECK(report.iterations == 0U);
  HUNDUN_CHECK(report.matvec_count == 1U);
  HUNDUN_CHECK(report.preconditioner_apply_count ==
               expected_preconditioners);
  HUNDUN_CHECK(report.global_reduction_count == 19U);
  HUNDUN_CHECK(std::isfinite(report.initial_residual));
  HUNDUN_CHECK(report.recursive_residual == report.initial_residual);
  HUNDUN_CHECK(report.final_residual == report.initial_residual);
  check_identical_report(report, world);
}

void test_synchronized_reduction_overflow_provenance(
    const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(kSystemSize, world.rank());
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  BiCGStabSolver bicgstab(execution, world);
  ConjugateGradientSolver cg(execution, world);

  const double norm_overflow_value =
      std::numeric_limits<double>::max() / 4.0;
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    b[index] = norm_overflow_value;
    x[index] = 0.0;
  }
  NonsymmetricOperator bicg_norm_operator(execution, layout, world.rank());
  InstrumentedPreconditioner bicg_norm_preconditioner(world.rank());
  const auto bicg_norm = bicgstab.solve(
      bicg_norm_operator, bicg_norm_preconditioner, b, x, {});
  check_global_norm_overflow_report(bicg_norm, world);
  HUNDUN_CHECK(bicg_norm_operator.apply_calls() == 0U);
  HUNDUN_CHECK(bicg_norm_preconditioner.update_calls() == 0U);

  NonsymmetricOperator cg_norm_operator(execution, layout, world.rank());
  InstrumentedPreconditioner cg_norm_preconditioner(world.rank());
  const auto cg_norm =
      cg.solve(cg_norm_operator, cg_norm_preconditioner, b, x, {});
  check_global_norm_overflow_report(cg_norm, world);
  HUNDUN_CHECK(cg_norm_operator.apply_calls() == 0U);
  HUNDUN_CHECK(cg_norm_preconditioner.update_calls() == 0U);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    HUNDUN_CHECK(x[index] == 0.0);
  }

  const int expected_lowest_rank = world.size() == 1 ? 0 : 1;
  const bool local_dot_overflow =
      world.size() == 1 || world.rank() == 1 ||
      (world.size() == 4 && world.rank() == 3);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    b[index] = local_dot_overflow ? 1.0e154 : 1.0;
    x[index] = 0.0;
  }
  NonsymmetricOperator bicg_dot_operator(execution, layout, world.rank());
  InstrumentedPreconditioner bicg_dot_preconditioner(world.rank());
  const auto bicg_dot = bicgstab.solve(
      bicg_dot_operator, bicg_dot_preconditioner, b, x, {});
  check_local_dot_overflow_report(bicg_dot, world, expected_lowest_rank, 0U);
  HUNDUN_CHECK(bicg_dot_preconditioner.update_calls() == 1U);

  NonsymmetricOperator cg_dot_operator(execution, layout, world.rank());
  InstrumentedPreconditioner cg_dot_preconditioner(world.rank());
  const auto cg_dot =
      cg.solve(cg_dot_operator, cg_dot_preconditioner, b, x, {});
  check_local_dot_overflow_report(cg_dot, world, expected_lowest_rank, 1U);
  HUNDUN_CHECK(cg_dot_preconditioner.update_calls() == 1U);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    HUNDUN_CHECK(x[index] == 0.0);
  }
}

void test_early_exits_and_breakdowns(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(5U, world.rank());
  NonsymmetricOperator linear_operator(execution, layout, world.rank());
  Buffer b_buffer(execution, 5U * sizeof(double));
  Buffer x_buffer(execution, 5U * sizeof(double));
  auto b = b_buffer.view(0U, 5U);
  auto x = x_buffer.view(0U, 5U);
  InstrumentedPreconditioner preconditioner(world.rank());
  BiCGStabSolver solver(execution, world);

  for (std::size_t index = 0; index < x.size(); ++index) {
    x[index] = 0.25 + static_cast<double>(index);
  }
  linear_operator.apply(x, b).wait();
  const auto exact = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(exact.reason == SolveTerminationReason::converged);
  HUNDUN_CHECK(exact.iterations == 0U);
  HUNDUN_CHECK(exact.final_residual == 0.0);
  HUNDUN_CHECK(exact.lowest_failing_rank == -1);
  HUNDUN_CHECK(preconditioner.update_calls() == 0U);

  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 0.0;
    x[index] = 7.0;
  }
  const auto zero = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(zero.reason == SolveTerminationReason::zero_right_hand_side);
  HUNDUN_CHECK(zero.iterations == 0U);
  HUNDUN_CHECK(zero.matvec_count == 2U);
  HUNDUN_CHECK(zero.final_residual == 0.0);
  for (std::size_t index = 0; index < x.size(); ++index) {
    HUNDUN_CHECK(x[index] == 0.0);
  }

  NonsymmetricOperator nonlinear_zero(execution, layout, world.rank());
  nonlinear_zero.set_bias(1.0);
  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 0.0;
    x[index] = 3.0;
  }
  const auto rejected_zero =
      solver.solve(nonlinear_zero, preconditioner, b, x, {});
  HUNDUN_CHECK(rejected_zero.reason ==
               SolveTerminationReason::numerical_breakdown);
  HUNDUN_CHECK(rejected_zero.iterations == 0U);
  HUNDUN_CHECK(rejected_zero.final_residual > 0.0);

  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 1.0;
    x[index] = 0.0;
  }
  SolveControl no_iterations;
  no_iterations.max_iterations = 0U;
  const auto none =
      solver.solve(linear_operator, preconditioner, b, x, no_iterations);
  HUNDUN_CHECK(none.reason == SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(none.iterations == 0U);
  HUNDUN_CHECK(none.matvec_count == 1U);

  NonsymmetricOperator zero_v(execution, layout, world.rank());
  zero_v.zero_output_on(2U);
  const auto denominator = solver.solve(zero_v, preconditioner, b, x, {});
  HUNDUN_CHECK(denominator.reason ==
               SolveTerminationReason::numerical_breakdown);
  HUNDUN_CHECK(denominator.iterations == 0U);
  HUNDUN_CHECK(std::isfinite(denominator.final_residual));

  NonsymmetricOperator zero_t(execution, layout, world.rank());
  zero_t.zero_output_on(3U);
  const auto tt = solver.solve(zero_t, preconditioner, b, x, {});
  HUNDUN_CHECK(tt.reason == SolveTerminationReason::numerical_breakdown);
  HUNDUN_CHECK(tt.iterations == 0U);
  HUNDUN_CHECK(std::isfinite(tt.final_residual));

  NonsymmetricOperator orthogonal_t(execution, layout, world.rank());
  orthogonal_t.orthogonal_output_on(3U);
  const auto omega = solver.solve(orthogonal_t, preconditioner, b, x, {});
  HUNDUN_CHECK(omega.reason == SolveTerminationReason::numerical_breakdown);
  HUNDUN_CHECK(omega.iterations == 0U);
  HUNDUN_CHECK(std::isfinite(omega.final_residual));

  NonsymmetricOperator tiny_v(execution, layout, world.rank());
  tiny_v.tiny_output_on(2U);
  const auto alpha = solver.solve(tiny_v, preconditioner, b, x, {});
  HUNDUN_CHECK(alpha.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(alpha.lowest_failing_rank == -1);
  HUNDUN_CHECK(std::isfinite(alpha.final_residual));

  const auto empty_layout = make_layout(0U, world.rank());
  NonsymmetricOperator empty_operator(execution, empty_layout, world.rank());
  Buffer empty_b(execution, 0U);
  Buffer empty_x(execution, 0U);
  InstrumentedPreconditioner empty_preconditioner(world.rank());
  const auto empty =
      solver.solve(empty_operator, empty_preconditioner, empty_b.view(0U, 0U),
                   empty_x.view(0U, 0U), {});
  HUNDUN_CHECK(empty.reason == SolveTerminationReason::zero_right_hand_side);
  HUNDUN_CHECK(empty.matvec_count == 2U);
}

void test_counter_batch_and_no_loop_allocation(const MpiContext& world) {
  CpuReferenceContext execution;
  NonsymmetricOperator linear_operator(
      execution, make_layout(kSystemSize, world.rank()), world.rank());
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    b[index] = index % 2U == 0U ? 1.0 : 2.0;
    x[index] = 0.0;
  }
  InstrumentedPreconditioner preconditioner(world.rank());
  preconditioner.arm_allocation_probe_after_update();
  BiCGStabSolver solver(execution, world);
  SolveControl control;
  control.atol = 0.0;
  control.rtol = 0.0;
  control.max_iterations = 1U;
  const auto reductions_before =
      world.fp64_reduction_counters().collective_calls;
  const auto report =
      solver.solve(linear_operator, preconditioner, b, x, control);
  const auto loop_allocations = allocation_probe::attempt_count;
  allocation_probe::count_attempts = false;

  HUNDUN_CHECK(report.reason == SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(report.iterations == 1U);
  HUNDUN_CHECK(report.matvec_count == 4U);
  HUNDUN_CHECK(report.preconditioner_apply_count == 2U);
  HUNDUN_CHECK(report.global_reduction_count == 46U);
  HUNDUN_CHECK(world.fp64_reduction_counters().collective_calls -
                   reductions_before ==
               46U);
  HUNDUN_CHECK(loop_allocations == 0U);
  HUNDUN_CHECK(std::isfinite(report.final_residual));
}

void test_rho_and_beta_mutations(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(3U, world.rank());
  Buffer b_buffer(execution, 3U * sizeof(double));
  Buffer x_buffer(execution, 3U * sizeof(double));
  auto b = b_buffer.view(0U, 3U);
  auto x = x_buffer.view(0U, 3U);
  InstrumentedPreconditioner preconditioner(world.rank());
  BiCGStabSolver solver(execution, world);
  SolveControl control;
  control.atol = 0.0;
  control.rtol = 0.0;
  control.max_iterations = 2U;

  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 0.0;
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  ScriptedRecurrenceOperator rho_operator(execution, layout,
                                          ScriptedRecurrence::rho_zero);
  const auto rho = solver.solve(rho_operator, preconditioner, b, x, control);
  HUNDUN_CHECK(rho.reason == SolveTerminationReason::numerical_breakdown);
  HUNDUN_CHECK(rho.iterations == 1U);
  HUNDUN_CHECK(rho.matvec_count == 4U);
  HUNDUN_CHECK(rho.lowest_failing_rank == -1);
  HUNDUN_CHECK(std::isfinite(rho.final_residual));

  b[0] = 1.0;
  b[1] = 1.0;
  b[2] = 1.0e-154;
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  ScriptedRecurrenceOperator beta_operator(execution, layout,
                                           ScriptedRecurrence::beta_overflow);
  const auto beta = solver.solve(beta_operator, preconditioner, b, x, control);
  HUNDUN_CHECK(beta.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(beta.iterations == 1U);
  HUNDUN_CHECK(beta.matvec_count == 4U);
  HUNDUN_CHECK(beta.lowest_failing_rank == -1);
  HUNDUN_CHECK(std::isfinite(beta.final_residual));
}

void test_periodic_replacement_restart(const MpiContext& world) {
  CpuReferenceContext execution;
  NonsymmetricOperator linear_operator(
      execution, make_layout(kSystemSize, world.rank()), world.rank());
  Buffer b_buffer(execution, kSystemSize * sizeof(double));
  Buffer x_buffer(execution, kSystemSize * sizeof(double));
  auto b = b_buffer.view(0U, kSystemSize);
  auto x = x_buffer.view(0U, kSystemSize);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    b[index] = index % 3U == 0U ? 2.0 : 1.0;
    x[index] = 0.0;
  }
  InstrumentedPreconditioner preconditioner(world.rank());
  BiCGStabSolver solver(execution, world);
  SolveControl control;
  control.atol = 0.0;
  control.rtol = 0.0;
  control.max_iterations = 2U;
  control.residual_recompute_interval = 1U;
  const auto report =
      solver.solve(linear_operator, preconditioner, b, x, control);

  HUNDUN_CHECK(report.reason == SolveTerminationReason::maximum_iterations);
  HUNDUN_CHECK(report.iterations == 2U);
  HUNDUN_CHECK(report.matvec_count == 7U);
  HUNDUN_CHECK(linear_operator.recorded_size(4U) == kSystemSize);
  HUNDUN_CHECK(linear_operator.recorded_size(5U) == kSystemSize);
  for (std::size_t index = 0; index < kSystemSize; ++index) {
    double applied = 4.0 * linear_operator.recorded(4U, index);
    if (index != 0U) {
      applied -= linear_operator.recorded(4U, index - 1U);
    }
    if (index + 1U != kSystemSize) {
      applied -= 2.0 * linear_operator.recorded(4U, index + 1U);
    }
    const double expected_restart = b[index] - applied;
    HUNDUN_CHECK(
        std::abs(linear_operator.recorded(5U, index) - expected_restart) <=
        16.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::abs(expected_restart)));
  }
}

void test_preflight_stride_overflow_and_mpi_error(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(3U, world.rank());
  NonsymmetricOperator linear_operator(execution, layout, world.rank());
  Buffer b_buffer(execution, 6U * sizeof(double));
  Buffer x_buffer(execution, 6U * sizeof(double));
  auto b = b_buffer.view(0U, 3U, 2U);
  auto x = x_buffer.view(0U, 3U, 2U);
  for (std::size_t index = 0; index < 3U; ++index) {
    b[index] = 2.0;
    x[index] = 0.0;
  }
  IdentityPreconditioner preconditioner(execution);
  BiCGStabSolver solver(execution, world);
  const auto strided = solver.solve(linear_operator, preconditioner, b, x, {});
  HUNDUN_CHECK(strided.reason == SolveTerminationReason::converged);
  for (std::size_t index = 0; index < 3U; ++index) {
    HUNDUN_CHECK(std::isfinite(x[index]));
  }

  Buffer shared(execution, 6U * sizeof(double));
  auto shared_x = shared.view(0U, 3U, 2U);
  const auto shared_b = static_cast<VectorView<const double>>(shared_x);
  const auto shared_report =
      solver.solve(linear_operator, preconditioner, shared_b, shared_x, {});
  HUNDUN_CHECK(shared_report.reason ==
               SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(shared_report.lowest_failing_rank == 0);

  Buffer stale_buffer(execution, 3U * sizeof(double));
  auto stale = stale_buffer.view(0U, 3U);
  stale_buffer.reallocate(3U * sizeof(double));
  const auto stale_report =
      solver.solve(linear_operator, preconditioner, stale, x, {});
  HUNDUN_CHECK(stale_report.reason ==
               SolveTerminationReason::collective_failure);

  const std::size_t too_large =
      std::numeric_limits<std::size_t>::max() / 9U + 1U;
  expect_runtime_error([&] {
    static_cast<void>(
        hundun::linear::detail::checked_bicgstab_workspace_bytes(too_large));
  });

  hundun::runtime::detail::inject_next_fp64_allreduce_result_for_test(
      MPI_ERR_OTHER);
  expect_runtime_error([&] {
    static_cast<void>(solver.solve(linear_operator, preconditioner, b, x, {}));
  });

  auto moved_source = MpiContext::duplicate(world.comm());
  auto live_destination = std::move(moved_source);
  expect_runtime_error(
      [&] { BiCGStabSolver moved_solver(execution, moved_source); });
  HUNDUN_CHECK(live_destination.size() == world.size());
}

void test_collective_phase_failures(const MpiContext& world) {
  if (world.size() < 2) {
    return;
  }
  CpuReferenceContext execution;
  const auto layout = make_layout(7U, world.rank());
  Buffer b_buffer(execution, 7U * sizeof(double));
  Buffer x_buffer(execution, 7U * sizeof(double));
  auto b = b_buffer.view(0U, 7U);
  auto x = x_buffer.view(0U, 7U);
  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 1.0 + static_cast<double>(index);
    x[index] = 0.0;
  }
  BiCGStabSolver solver(execution, world);
  ConjugateGradientSolver cg(execution, world);

  if (world.rank() == 1) {
    ExecutionTestAccess::fail_next_allocation();
  }
  NonsymmetricOperator allocation_operator(execution, layout, world.rank());
  InstrumentedPreconditioner allocation_preconditioner(world.rank());
  const auto allocation =
      solver.solve(allocation_operator, allocation_preconditioner, b, x, {});
  HUNDUN_CHECK(allocation.reason == SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(allocation.lowest_failing_rank == 1);
  HUNDUN_CHECK(allocation.matvec_count == 0U);
  check_identical_report(allocation, world);

  const auto other_range = make_layout(8U, world.rank());
  NonsymmetricOperator bad_layout(execution, layout,
                                  world.rank() == 1 ? other_range : layout,
                                  world.rank());
  InstrumentedPreconditioner layout_preconditioner(world.rank());
  const auto layout_report =
      solver.solve(bad_layout, layout_preconditioner, b, x, {});
  HUNDUN_CHECK(layout_report.reason ==
               SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(layout_report.lowest_failing_rank == 1);
  check_identical_report(layout_report, world);
  InstrumentedPreconditioner cg_layout_preconditioner(world.rank());
  const auto cg_layout =
      cg.solve(bad_layout, cg_layout_preconditioner, b, x, {});
  HUNDUN_CHECK(cg_layout.reason == SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(cg_layout.lowest_failing_rank == 1);
  check_identical_report(cg_layout, world);

  NonsymmetricOperator event_operator(execution, layout, world.rank());
  event_operator.fail_event_on(1, 1U);
  InstrumentedPreconditioner event_preconditioner(world.rank());
  const auto event =
      solver.solve(event_operator, event_preconditioner, b, x, {});
  HUNDUN_CHECK(event.reason == SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(event.lowest_failing_rank == 1);
  HUNDUN_CHECK(std::isinf(event.final_residual));
  check_identical_report(event, world);

  NonsymmetricOperator nonfinite_operator(execution, layout, world.rank());
  nonfinite_operator.nonfinite_output_on(1, 1U);
  InstrumentedPreconditioner nonfinite_preconditioner(world.rank());
  const auto nonfinite =
      solver.solve(nonfinite_operator, nonfinite_preconditioner, b, x, {});
  HUNDUN_CHECK(nonfinite.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(nonfinite.lowest_failing_rank == 1);
  check_identical_report(nonfinite, world);

  NonsymmetricOperator update_operator(execution, layout, world.rank());
  InstrumentedPreconditioner update_preconditioner(world.rank());
  update_preconditioner.fail_update_on(1);
  const auto update =
      solver.solve(update_operator, update_preconditioner, b, x, {});
  HUNDUN_CHECK(update.reason == SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(update.lowest_failing_rank == 1);
  check_identical_report(update, world);

  NonsymmetricOperator pre_event_operator(execution, layout, world.rank());
  InstrumentedPreconditioner pre_event(world.rank());
  pre_event.fail_event_on(1, 1U);
  const auto pre_event_report =
      solver.solve(pre_event_operator, pre_event, b, x, {});
  HUNDUN_CHECK(pre_event_report.reason ==
               SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(pre_event_report.lowest_failing_rank == 1);
  check_identical_report(pre_event_report, world);

  NonsymmetricOperator pre_nonfinite_operator(execution, layout, world.rank());
  InstrumentedPreconditioner pre_nonfinite(world.rank());
  pre_nonfinite.nonfinite_output_on(1, 1U);
  const auto pre_nonfinite_report =
      solver.solve(pre_nonfinite_operator, pre_nonfinite, b, x, {});
  HUNDUN_CHECK(pre_nonfinite_report.reason ==
               SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(pre_nonfinite_report.lowest_failing_rank == 1);
  check_identical_report(pre_nonfinite_report, world);

  NonsymmetricOperator late_event_operator(execution, layout, world.rank());
  late_event_operator.fail_event_on(1, 4U);
  InstrumentedPreconditioner late_event_preconditioner(world.rank());
  SolveControl late_control;
  late_control.atol = 0.0;
  late_control.rtol = 0.0;
  late_control.max_iterations = 3U;
  const auto late_event = solver.solve(
      late_event_operator, late_event_preconditioner, b, x, late_control);
  HUNDUN_CHECK(late_event.reason == SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(late_event.lowest_failing_rank == 1);
  HUNDUN_CHECK(late_event.iterations == 1U);
  HUNDUN_CHECK(std::isinf(late_event.final_residual));
  check_identical_report(late_event, world);

  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 1.0;
    x[index] = 0.0;
  }
  if (world.rank() == 1 || world.rank() == 3) {
    b[0] = std::numeric_limits<double>::quiet_NaN();
  }
  NonsymmetricOperator multiple_operator(execution, layout, world.rank());
  InstrumentedPreconditioner multiple_preconditioner(world.rank());
  const auto multiple =
      solver.solve(multiple_operator, multiple_preconditioner, b, x, {});
  HUNDUN_CHECK(multiple.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(multiple.lowest_failing_rank == 1);
  check_identical_report(multiple, world);

  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 1.0;
    x[index] = 0.0;
  }
  NonsymmetricOperator cg_operator(execution, layout, world.rank());
  InstrumentedPreconditioner cg_preconditioner(world.rank());
  if (world.rank() == 1) {
    b[0] = std::numeric_limits<double>::infinity();
  }
  const auto cg_nonfinite = cg.solve(cg_operator, cg_preconditioner, b, x, {});
  HUNDUN_CHECK(cg_nonfinite.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(cg_nonfinite.lowest_failing_rank == 1);
  check_identical_report(cg_nonfinite, world);

  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 1.0;
    x[index] = 0.0;
  }
  NonsymmetricOperator cg_event_operator(execution, layout, world.rank());
  cg_event_operator.fail_event_on(1, 1U);
  InstrumentedPreconditioner cg_event_preconditioner(world.rank());
  const auto cg_event =
      cg.solve(cg_event_operator, cg_event_preconditioner, b, x, {});
  HUNDUN_CHECK(cg_event.reason == SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(cg_event.lowest_failing_rank == 1);
  check_identical_report(cg_event, world);

  NonsymmetricOperator cg_pre_operator(execution, layout, world.rank());
  InstrumentedPreconditioner cg_preconditioner_event(world.rank());
  cg_preconditioner_event.fail_event_on(1, 1U);
  const auto cg_pre_event =
      cg.solve(cg_pre_operator, cg_preconditioner_event, b, x, {});
  HUNDUN_CHECK(cg_pre_event.reason ==
               SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(cg_pre_event.lowest_failing_rank == 1);
  check_identical_report(cg_pre_event, world);
}

void test_post_iteration_refresh_preserves_failure(const MpiContext& world) {
  CpuReferenceContext execution;
  const auto layout = make_layout(7U, world.rank());
  Buffer b_buffer(execution, 7U * sizeof(double));
  Buffer x_buffer(execution, 7U * sizeof(double));
  auto b = b_buffer.view(0U, 7U);
  auto x = x_buffer.view(0U, 7U);
  for (std::size_t index = 0; index < b.size(); ++index) {
    b[index] = 1.0 + static_cast<double>(index);
    x[index] = 0.0;
  }
  const int failure_rank = world.size() == 1 ? 0 : 1;
  SolveControl control;
  control.atol = 0.0;
  control.rtol = 0.0;
  control.max_iterations = 3U;

  NonsymmetricOperator cg_operator(execution, layout, world.rank());
  InstrumentedPreconditioner cg_preconditioner(world.rank());
  cg_preconditioner.fail_event_on(failure_rank, 2U);
  ConjugateGradientSolver cg(execution, world);
  const auto cg_report =
      cg.solve(cg_operator, cg_preconditioner, b, x, control);
  HUNDUN_CHECK(cg_report.reason ==
               SolveTerminationReason::collective_failure);
  HUNDUN_CHECK(cg_report.iterations == 1U);
  HUNDUN_CHECK(cg_report.lowest_failing_rank == failure_rank);
  HUNDUN_CHECK(std::isfinite(cg_report.final_residual));
  double cg_local_residual_squared = 0.0;
  for (std::size_t index = 0; index < b.size(); ++index) {
    double applied = 4.0 * x[index];
    if (index != 0U) {
      applied -= x[index - 1U];
    }
    if (index + 1U != x.size()) {
      applied -= 2.0 * x[index + 1U];
    }
    const double residual = b[index] - applied;
    cg_local_residual_squared += residual * residual;
  }
  const double cg_measured = global_l2(cg_local_residual_squared, world);
  HUNDUN_CHECK(
      std::abs(cg_report.final_residual - cg_measured) <=
      64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, cg_measured));
  check_identical_report(cg_report, world);

  for (std::size_t index = 0; index < x.size(); ++index) {
    x[index] = 0.0;
  }
  CandidateOverflowState candidate_overflow;
  NonsymmetricOperator bicg_operator(execution, layout, world.rank());
  bicg_operator.overflow_solution_candidate_after_apply_on(failure_rank, 5U,
                                                            x, candidate_overflow);
  InstrumentedPreconditioner bicg_preconditioner(world.rank());
  bicg_preconditioner.capture_solution_candidate_on(4U, candidate_overflow);
  BiCGStabSolver bicg(execution, world);
  const auto bicg_report =
      bicg.solve(bicg_operator, bicg_preconditioner, b, x, control);
  HUNDUN_CHECK(bicg_report.reason == SolveTerminationReason::non_finite_value);
  HUNDUN_CHECK(bicg_report.iterations == 1U);
  HUNDUN_CHECK(bicg_report.lowest_failing_rank == failure_rank);
  HUNDUN_CHECK(std::isfinite(bicg_report.final_residual));
  double bicg_local_residual_squared = 0.0;
  for (std::size_t index = 0; index < b.size(); ++index) {
    double applied = 0.0;
    if (world.rank() != failure_rank) {
      applied = 4.0 * x[index];
      if (index != 0U) {
        applied -= x[index - 1U];
      }
      if (index + 1U != x.size()) {
        applied -= 2.0 * x[index + 1U];
      }
    }
    const double residual = b[index] - applied;
    bicg_local_residual_squared += residual * residual;
  }
  const double bicg_measured = global_l2(bicg_local_residual_squared, world);
  HUNDUN_CHECK(
      std::abs(bicg_report.final_residual - bicg_measured) <=
      64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, bicg_measured));
  check_identical_report(bicg_report, world);
}

}  // namespace

int main(int argc, char** argv) {
  MpiEnvironment environment(argc, argv);
  auto world = MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    test_manufactured(world);
    test_collective_control_and_input(world);
    test_tolerance_overflow_is_global(world);
    test_synchronized_reduction_overflow_provenance(world);
    test_early_exits_and_breakdowns(world);
    test_counter_batch_and_no_loop_allocation(world);
    test_rho_and_beta_mutations(world);
    test_periodic_replacement_restart(world);
    test_preflight_stride_overflow_and_mpi_error(world);
    test_collective_phase_failures(world);
    test_post_iteration_refresh_preserves_failure(world);
  });
}
