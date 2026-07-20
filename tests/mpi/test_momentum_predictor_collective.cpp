// SPDX-License-Identifier: Apache-2.0

#include "execution/src/execution_test_access.hpp"
#include "flow/src/momentum_predictor_test_access.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/flow/momentum_predictor.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "runtime/src/mpi_error.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {

using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionEvent;
using hundun::execution::VectorView;
using hundun::execution::test::ExecutionTestAccess;
using hundun::flow::MomentumComponentEquation;
using hundun::flow::MomentumPredictor;
using hundun::flow::test::MomentumPredictorTestAccess;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::SolveTerminationReason;
using hundun::linear::VectorLayout;
using hundun::runtime::Error;
using hundun::runtime::Fp64ReductionOperation;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;

class NoopPreconditioner final : public Preconditioner {
public:
  void update(const LinearOperator &, std::uint64_t) override {}
  ExecutionEvent apply(VectorView<const double>,
                       VectorView<double>) const override {
    return ExecutionEvent::completed();
  }
};

class TestOperator final : public LinearOperator {
public:
  TestOperator(ExecutionContext &context, VectorLayout domain,
               VectorLayout range, std::vector<double> diagonal)
      : context_(&context), domain_(std::move(domain)),
        range_(std::move(range)), diagonal_(std::move(diagonal)) {}

  VectorLayout domain_layout() const override { return domain_; }
  VectorLayout range_layout() const override { return range_; }
  const ExecutionContext &context() const override { return *context_; }
  std::uint64_t revision() const override { return 1U; }
  ExecutionEvent apply(VectorView<const double> x,
                       VectorView<double> y) const override {
    for (std::size_t index = 0; index < x.size(); ++index)
      y[index] = x[index];
    return ExecutionEvent::completed();
  }
  bool has_diagonal() const override { return has_diagonal_; }
  ExecutionEvent diagonal(VectorView<double> output) const override {
    if (throw_mpi_operation_error_) {
      throw hundun::runtime::detail::MpiOperationError(
          "injected typed MPI operation failure");
    }
    if (!diagonal_error_.empty())
      throw Error(diagonal_error_);
    for (std::size_t index = 0; index < diagonal_.size(); ++index)
      output[index] = diagonal_[index];
    if (failed_event_) {
      auto event = ExecutionTestAccess::pending_event({});
      ExecutionTestAccess::complete_error(event,
                                          "injected diagonal event failure");
      return event;
    }
    return ExecutionEvent::completed();
  }
  void set_has_diagonal(bool value) noexcept { has_diagonal_ = value; }
  void set_failed_event(bool value) noexcept { failed_event_ = value; }
  void set_diagonal_error(std::string value) {
    diagonal_error_ = std::move(value);
  }
  void set_throw_mpi_operation_error(bool value) noexcept {
    throw_mpi_operation_error_ = value;
  }

private:
  ExecutionContext *context_;
  VectorLayout domain_;
  VectorLayout range_;
  std::vector<double> diagonal_;
  bool has_diagonal_{true};
  bool failed_event_{false};
  bool throw_mpi_operation_error_{false};
  std::string diagonal_error_;
};

class CountingSolver final : public LinearSolver {
public:
  SolveReport solve(const LinearOperator &, Preconditioner &,
                    VectorView<const double>, VectorView<double>,
                    const SolveControl &) const override {
    ++calls_;
    SolveReport report{};
    report.reason = SolveTerminationReason::converged;
    return report;
  }
  std::size_t calls() const noexcept { return calls_; }

private:
  mutable std::size_t calls_{};
};

struct Fixture final {
  static constexpr std::size_t count = 2U;

  explicit Fixture(const MpiContext &mpi)
      : layout(count, std::vector<hundun::mesh::GlobalCellId>{0U, 1U}),
        other_layout(count, std::vector<hundun::mesh::GlobalCellId>{1U, 0U}),
        op0(context, layout, layout, {2.0, 3.0}),
        op1(context, layout, layout, {4.0, 5.0}),
        op2(context, layout, layout, {6.0, 7.0}), rhs0(context, 16U),
        rhs1(context, 16U), rhs2(context, 16U), x0(context, 16U),
        x1(context, 16U), x2(context, 16U), d0(context, 16U), d1(context, 16U),
        d2(context, 16U) {
    static_cast<void>(mpi);
    reset();
  }

  std::array<MomentumComponentEquation, 3> equations() {
    return {{{&op0, &p0, static_cast<const Buffer &>(rhs0).view(0U, count),
              x0.view(0U, count), d0.view(0U, count)},
             {&op1, &p1, static_cast<const Buffer &>(rhs1).view(0U, count),
              x1.view(0U, count), d1.view(0U, count)},
             {&op2, &p2, static_cast<const Buffer &>(rhs2).view(0U, count),
              x2.view(0U, count), d2.view(0U, count)}}};
  }

  void reset() {
    for (Buffer *buffer : std::array<Buffer *, 3>{&x0, &x1, &x2}) {
      auto view = buffer->view(0U, count);
      for (std::size_t index = 0; index < count; ++index)
        view[index] = -11.0;
    }
    for (Buffer *buffer : std::array<Buffer *, 3>{&d0, &d1, &d2}) {
      auto view = buffer->view(0U, count);
      for (std::size_t index = 0; index < count; ++index)
        view[index] = -22.0;
    }
  }

  void check_sentinels() const {
    for (const Buffer *buffer : std::array<const Buffer *, 3>{&x0, &x1, &x2}) {
      const auto view = buffer->view(0U, count);
      for (std::size_t index = 0; index < count; ++index)
        HUNDUN_CHECK(view[index] == -11.0);
    }
    for (const Buffer *buffer : std::array<const Buffer *, 3>{&d0, &d1, &d2}) {
      const auto view = buffer->view(0U, count);
      for (std::size_t index = 0; index < count; ++index)
        HUNDUN_CHECK(view[index] == -22.0);
    }
  }

  CpuReferenceContext context;
  VectorLayout layout;
  VectorLayout other_layout;
  TestOperator op0;
  TestOperator op1;
  TestOperator op2;
  NoopPreconditioner p0;
  NoopPreconditioner p1;
  NoopPreconditioner p2;
  CountingSolver solver;
  Buffer rhs0;
  Buffer rhs1;
  Buffer rhs2;
  Buffer x0;
  Buffer x1;
  Buffer x2;
  Buffer d0;
  Buffer d1;
  Buffer d2;
};

void require_no_solver_calls(const MpiContext &mpi, const Fixture &fixture) {
  double calls = static_cast<double>(fixture.solver.calls());
  mpi.allreduce_fp64_in_place(&calls, 1U, Fp64ReductionOperation::sum);
  HUNDUN_CHECK(calls == 0.0);
  fixture.check_sentinels();
}

void check_null_pointer_is_collective(const MpiContext &mpi) {
  Fixture fixture(mpi);
  auto equations = fixture.equations();
  if (mpi.rank() == 1)
    equations[0].linear_operator = nullptr;
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(predictor.solve(mpi, equations, SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(message ==
               "momentum predictor preflight failed on rank 1: momentum "
               "equation has a null operator or preconditioner");
  require_no_solver_calls(mpi, fixture);
}

void check_lowest_failure_rank_wins(const MpiContext &mpi) {
  Fixture fixture(mpi);
  if (mpi.rank() == 1)
    fixture.op1.set_has_diagonal(false);
  if (mpi.rank() == 2) {
    fixture.op2 = TestOperator(fixture.context, fixture.layout, fixture.layout,
                               {6.0, 0.0});
  }
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(message ==
               "momentum predictor preflight failed on rank 1: momentum "
               "operator lacks diagonal capability");
  require_no_solver_calls(mpi, fixture);
}

void check_common_layout_mismatch_is_collective(const MpiContext &mpi) {
  Fixture fixture(mpi);
  if (mpi.rank() == 1) {
    fixture.op1 = TestOperator(fixture.context, fixture.other_layout,
                               fixture.other_layout, {4.0, 5.0});
  }
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(message ==
               "momentum predictor preflight failed on rank 1: momentum "
               "equation layouts are not identical");
  require_no_solver_calls(mpi, fixture);
}

void check_nonsquare_layout_is_collective(const MpiContext &mpi) {
  Fixture fixture(mpi);
  if (mpi.rank() == 1) {
    fixture.op1 = TestOperator(fixture.context, fixture.layout,
                               fixture.other_layout, {4.0, 5.0});
  }
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(message ==
               "momentum predictor preflight failed on rank 1: momentum "
               "operator is not square");
  require_no_solver_calls(mpi, fixture);
}

void check_event_failure_is_collective(const MpiContext &mpi) {
  Fixture fixture(mpi);
  if (mpi.rank() == 1)
    fixture.op1.set_failed_event(true);
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(message ==
               "momentum predictor preflight failed on rank 1: injected "
               "diagonal event failure");
  require_no_solver_calls(mpi, fixture);
}

void check_staging_allocation_failure_is_collective(const MpiContext &mpi) {
  Fixture fixture(mpi);
  if (mpi.rank() == 1)
    ExecutionTestAccess::fail_next_allocation();
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(message ==
               "momentum predictor preflight failed on rank 1: buffer "
               "allocation failed by the test seam");
  require_no_solver_calls(mpi, fixture);
}

void check_long_diagnostic_is_deterministically_bounded(const MpiContext &mpi) {
  Fixture fixture(mpi);
  if (mpi.rank() == 1)
    fixture.op1.set_diagonal_error(std::string(2048U, 'q'));
  MomentumPredictor predictor(fixture.solver);
  std::string message;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const Error &error) {
    message = error.what();
  }
  const std::string prefix = "momentum predictor preflight failed on rank 1: ";
  HUNDUN_CHECK(message.rfind(prefix, 0U) == 0U);
  HUNDUN_CHECK(message.size() == prefix.size() + 511U);
  for (std::size_t index = prefix.size(); index < message.size(); ++index)
    HUNDUN_CHECK(message[index] == 'q');
  require_no_solver_calls(mpi, fixture);
}

void check_typed_mpi_operation_error_is_rethrown(const MpiContext &mpi) {
  Fixture fixture(mpi);
  fixture.op0.set_throw_mpi_operation_error(true);
  MomentumPredictor predictor(fixture.solver);
  MomentumPredictorTestAccess::reset_collective_selection_calls();
  bool caught_typed = false;
  try {
    static_cast<void>(
        predictor.solve(mpi, fixture.equations(), SolveControl{}));
  } catch (const hundun::runtime::detail::MpiOperationError &error) {
    caught_typed =
        std::string(error.what()) == "injected typed MPI operation failure";
  }
  HUNDUN_CHECK(caught_typed);
  HUNDUN_CHECK(MomentumPredictorTestAccess::collective_selection_calls() == 0U);
  require_no_solver_calls(mpi, fixture);
}

} // namespace

int main(int argc, char **argv) {
  MpiEnvironment environment(argc, argv);
  MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    HUNDUN_CHECK(mpi.size() == 2 || mpi.size() == 4);
    check_null_pointer_is_collective(mpi);
    check_lowest_failure_rank_wins(mpi);
    check_common_layout_mismatch_is_collective(mpi);
    check_nonsquare_layout_is_collective(mpi);
    check_event_failure_is_collective(mpi);
    check_staging_allocation_failure_is_collective(mpi);
    check_long_diagnostic_is_deterministically_bounded(mpi);
    check_typed_mpi_operation_error_is_rethrown(mpi);
  });
}
