// SPDX-License-Identifier: Apache-2.0

#include "execution/src/execution_test_access.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/flow/momentum_predictor.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/runtime/error.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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
using hundun::flow::make_momentum_time_stencil;
using hundun::flow::MomentumComponentEquation;
using hundun::flow::MomentumPredictor;
using hundun::flow::MomentumTimeOrder;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::SolveTerminationReason;
using hundun::linear::VectorLayout;
using hundun::runtime::Error;

template <class Function> void expect_error(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &) {
    threw = true;
  }
  HUNDUN_CHECK(threw);
}

class NoopPreconditioner final : public Preconditioner {
public:
  void update(const LinearOperator &, std::uint64_t) override {}
  ExecutionEvent apply(VectorView<const double>,
                       VectorView<double>) const override {
    return ExecutionEvent::completed();
  }
};

class DiagonalOperator final : public LinearOperator {
public:
  DiagonalOperator(ExecutionContext &context, VectorLayout layout,
                   std::vector<double> diagonal)
      : context_(&context), layout_(std::move(layout)),
        diagonal_(std::move(diagonal)) {}

  VectorLayout domain_layout() const override { return layout_; }
  VectorLayout range_layout() const override { return layout_; }
  const ExecutionContext &context() const override { return *context_; }
  std::uint64_t revision() const override { return 7U; }
  ExecutionEvent apply(VectorView<const double> x,
                       VectorView<double> y) const override {
    for (std::size_t i = 0; i < x.size(); ++i) {
      y[i] = diagonal_[i] * x[i];
    }
    return ExecutionEvent::completed();
  }
  bool has_diagonal() const override { return has_diagonal_; }
  ExecutionEvent diagonal(VectorView<double> output) const override {
    ++diagonal_calls_;
    if (throw_from_diagonal_) {
      throw Error("injected diagonal failure");
    }
    for (std::size_t i = 0; i < diagonal_.size(); ++i) {
      output[i] = diagonal_[i];
    }
    if (failed_event_) {
      auto event = ExecutionTestAccess::pending_event({});
      ExecutionTestAccess::complete_error(event, "injected event failure");
      return event;
    }
    return ExecutionEvent::completed();
  }

  void set_has_diagonal(bool value) noexcept { has_diagonal_ = value; }
  void set_throw_from_diagonal(bool value) noexcept {
    throw_from_diagonal_ = value;
  }
  void set_failed_event(bool value) noexcept { failed_event_ = value; }
  std::size_t diagonal_calls() const noexcept { return diagonal_calls_; }

private:
  ExecutionContext *context_;
  VectorLayout layout_;
  std::vector<double> diagonal_;
  bool has_diagonal_{true};
  bool throw_from_diagonal_{false};
  bool failed_event_{false};
  mutable std::size_t diagonal_calls_{};
};

class RecordingSolver final : public LinearSolver {
public:
  std::array<const LinearOperator *, 3> expected_operators{};
  std::array<Preconditioner *, 3> expected_preconditioners{};
  std::array<SolveTerminationReason, 3> reasons{
      SolveTerminationReason::converged,
      SolveTerminationReason::zero_right_hand_side,
      SolveTerminationReason::converged};

  SolveReport solve(const LinearOperator &linear_operator,
                    Preconditioner &preconditioner,
                    VectorView<const double> rhs, VectorView<double> predictor,
                    const SolveControl &) const override {
    HUNDUN_CHECK(call_count_ < 3U);
    HUNDUN_CHECK(expected_operators[call_count_] == &linear_operator);
    HUNDUN_CHECK(expected_preconditioners[call_count_] == &preconditioner);
    HUNDUN_CHECK(rhs.size() == predictor.size());
    for (std::size_t i = 0; i < predictor.size(); ++i) {
      predictor[i] = 100.0 * static_cast<double>(call_count_ + 1U) + rhs[i];
    }
    SolveReport report;
    report.reason = reasons[call_count_];
    report.iterations = 10U + call_count_;
    report.initial_residual = 3.0;
    report.recursive_residual = 2.0;
    report.final_residual = 1.0;
    report.matvec_count = 20U + call_count_;
    ++call_count_;
    return report;
  }

  std::size_t call_count() const noexcept { return call_count_; }

private:
  mutable std::size_t call_count_{};
};

struct PredictorFixture final {
  static constexpr std::size_t count = 3U;

  PredictorFixture()
      : layout(count, std::vector<hundun::mesh::GlobalCellId>{0U, 1U, 2U}),
        op0(context, layout, {2.0, 3.0, 4.0}),
        op1(context, layout, {5.0, 6.0, 7.0}),
        op2(context, layout, {8.0, 9.0, 10.0}),
        rhs0(context, count * sizeof(double)),
        rhs1(context, count * sizeof(double)),
        rhs2(context, count * sizeof(double)),
        x0(context, count * sizeof(double)),
        x1(context, count * sizeof(double)),
        x2(context, count * sizeof(double)),
        d0(context, count * sizeof(double)),
        d1(context, count * sizeof(double)),
        d2(context, count * sizeof(double)) {
    operators = {&op0, &op1, &op2};
    preconditioners = {&p0, &p1, &p2};
    solver.expected_operators = operators;
    solver.expected_preconditioners = preconditioners;
    reset_outputs();
    for (std::size_t c = 0; c < 3U; ++c) {
      auto rhs = rhs_buffers()[c]->view(0U, count);
      for (std::size_t i = 0; i < count; ++i) {
        rhs[i] = static_cast<double>(10U * c + i + 1U);
      }
    }
  }

  std::array<Buffer *, 3> rhs_buffers() { return {&rhs0, &rhs1, &rhs2}; }
  std::array<Buffer *, 3> x_buffers() { return {&x0, &x1, &x2}; }
  std::array<Buffer *, 3> d_buffers() { return {&d0, &d1, &d2}; }

  std::array<MomentumComponentEquation, 3> equations() {
    const auto rhs = rhs_buffers();
    const auto x = x_buffers();
    const auto d = d_buffers();
    return {{{operators[0], preconditioners[0],
              static_cast<const Buffer &>(*rhs[0]).view(0U, count),
              x[0]->view(0U, count), d[0]->view(0U, count)},
             {operators[1], preconditioners[1],
              static_cast<const Buffer &>(*rhs[1]).view(0U, count),
              x[1]->view(0U, count), d[1]->view(0U, count)},
             {operators[2], preconditioners[2],
              static_cast<const Buffer &>(*rhs[2]).view(0U, count),
              x[2]->view(0U, count), d[2]->view(0U, count)}}};
  }

  void reset_outputs() {
    for (Buffer *buffer : x_buffers()) {
      auto view = buffer->view(0U, count);
      for (std::size_t i = 0; i < count; ++i)
        view[i] = -11.0;
    }
    for (Buffer *buffer : d_buffers()) {
      auto view = buffer->view(0U, count);
      for (std::size_t i = 0; i < count; ++i)
        view[i] = -22.0;
    }
  }

  void check_unchanged() {
    for (Buffer *buffer : x_buffers()) {
      const auto view = static_cast<const Buffer &>(*buffer).view(0U, count);
      for (std::size_t i = 0; i < count; ++i)
        HUNDUN_CHECK(view[i] == -11.0);
    }
    for (Buffer *buffer : d_buffers()) {
      const auto view = static_cast<const Buffer &>(*buffer).view(0U, count);
      for (std::size_t i = 0; i < count; ++i)
        HUNDUN_CHECK(view[i] == -22.0);
    }
  }

  CpuReferenceContext context;
  VectorLayout layout;
  DiagonalOperator op0;
  DiagonalOperator op1;
  DiagonalOperator op2;
  NoopPreconditioner p0;
  NoopPreconditioner p1;
  NoopPreconditioner p2;
  RecordingSolver solver;
  std::array<const LinearOperator *, 3> operators{};
  std::array<Preconditioner *, 3> preconditioners{};
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

void check_time_stencils() {
  const auto be = make_momentum_time_stencil(MomentumTimeOrder::backward_euler,
                                             0.125, 91.0);
  HUNDUN_CHECK(be.order == MomentumTimeOrder::backward_euler);
  HUNDUN_CHECK(be.dt_s == 0.125);
  HUNDUN_CHECK(be.previous_dt_s == 0.0);
  HUNDUN_CHECK(be.alpha0 == 1.0);
  HUNDUN_CHECK(be.alpha1 == -1.0);
  HUNDUN_CHECK(be.alpha2 == 0.0);

  for (const auto &[ratio, expected] :
       std::array<std::pair<double, std::array<double, 3>>, 3>{
           std::pair{0.5, std::array<double, 3>{4.0 / 3.0, -1.5, 1.0 / 6.0}},
           std::pair{1.0, std::array<double, 3>{1.5, -2.0, 0.5}},
           std::pair{2.0, std::array<double, 3>{5.0 / 3.0, -3.0, 4.0 / 3.0}}}) {
    const auto bdf =
        make_momentum_time_stencil(MomentumTimeOrder::bdf2, ratio * 2.0, 2.0);
    HUNDUN_CHECK(bdf.order == MomentumTimeOrder::bdf2);
    HUNDUN_CHECK(bdf.dt_s == ratio * 2.0);
    HUNDUN_CHECK(bdf.previous_dt_s == 2.0);
    HUNDUN_CHECK(bdf.alpha0 == expected[0]);
    HUNDUN_CHECK(bdf.alpha1 == expected[1]);
    HUNDUN_CHECK(bdf.alpha2 == expected[2]);
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (double invalid : {0.0, -1.0, nan, inf}) {
    expect_error([&] {
      static_cast<void>(make_momentum_time_stencil(
          MomentumTimeOrder::backward_euler, invalid, 1.0));
    });
  }
  expect_error([&] {
    static_cast<void>(make_momentum_time_stencil(
        static_cast<MomentumTimeOrder>(99), 1.0, 1.0));
  });
  for (const auto &pair : std::array<std::pair<double, double>, 6>{
           std::pair{1.0, 0.0}, std::pair{1.0, -1.0}, std::pair{1.0, nan},
           std::pair{1.0, inf}, std::pair{0.49, 1.0}, std::pair{2.01, 1.0}}) {
    expect_error([&] {
      static_cast<void>(make_momentum_time_stencil(MomentumTimeOrder::bdf2,
                                                   pair.first, pair.second));
    });
  }
}

void check_predictor_success_and_reports() {
  PredictorFixture fixture;
  MomentumPredictor predictor(fixture.solver);
  auto equations = fixture.equations();
  const auto report = predictor.solve(equations, SolveControl{});
  HUNDUN_CHECK(fixture.solver.call_count() == 3U);
  HUNDUN_CHECK(report.all_converged());
  for (std::size_t c = 0; c < 3U; ++c) {
    HUNDUN_CHECK(report.components[c].iterations == 10U + c);
    HUNDUN_CHECK(report.components[c].matvec_count == 20U + c);
    const auto diagonal = static_cast<const Buffer &>(*fixture.d_buffers()[c])
                              .view(0U, PredictorFixture::count);
    const auto x = static_cast<const Buffer &>(*fixture.x_buffers()[c])
                       .view(0U, PredictorFixture::count);
    for (std::size_t i = 0; i < PredictorFixture::count; ++i) {
      HUNDUN_CHECK(diagonal[i] ==
                   2.0 + 3.0 * static_cast<double>(c) + static_cast<double>(i));
      HUNDUN_CHECK(x[i] == 100.0 * static_cast<double>(c + 1U) +
                               static_cast<double>(10U * c + i + 1U));
    }
  }

  fixture.solver.reasons[1] = SolveTerminationReason::maximum_iterations;
  PredictorFixture failed;
  failed.solver.reasons[1] = SolveTerminationReason::maximum_iterations;
  MomentumPredictor failed_predictor(failed.solver);
  const auto failed_report =
      failed_predictor.solve(failed.equations(), SolveControl{});
  HUNDUN_CHECK(failed.solver.call_count() == 3U);
  HUNDUN_CHECK(!failed_report.all_converged());
  HUNDUN_CHECK(failed_report.components[1].reason ==
               SolveTerminationReason::maximum_iterations);
}

void check_predictor_preflight_transaction() {
  {
    PredictorFixture fixture;
    fixture.op1.set_has_diagonal(false);
    MomentumPredictor predictor(fixture.solver);
    expect_error([&] {
      static_cast<void>(predictor.solve(fixture.equations(), SolveControl{}));
    });
    HUNDUN_CHECK(fixture.solver.call_count() == 0U);
    HUNDUN_CHECK(fixture.op2.diagonal_calls() == 0U);
    fixture.check_unchanged();
  }
  {
    PredictorFixture fixture;
    fixture.op1.set_failed_event(true);
    MomentumPredictor predictor(fixture.solver);
    expect_error([&] {
      static_cast<void>(predictor.solve(fixture.equations(), SolveControl{}));
    });
    HUNDUN_CHECK(fixture.solver.call_count() == 0U);
    fixture.check_unchanged();
  }
  {
    PredictorFixture fixture;
    auto equations = fixture.equations();
    equations[1].linear_operator = nullptr;
    MomentumPredictor predictor(fixture.solver);
    expect_error(
        [&] { static_cast<void>(predictor.solve(equations, SolveControl{})); });
    fixture.check_unchanged();
  }
  {
    PredictorFixture fixture;
    auto equations = fixture.equations();
    equations[0].actual_diagonal = equations[0].predictor;
    MomentumPredictor predictor(fixture.solver);
    expect_error(
        [&] { static_cast<void>(predictor.solve(equations, SolveControl{})); });
    fixture.check_unchanged();
  }
  {
    PredictorFixture fixture;
    auto equations = fixture.equations();
    Buffer stale(fixture.context, PredictorFixture::count * sizeof(double));
    equations[0].predictor = stale.view(0U, PredictorFixture::count);
    stale.reallocate(PredictorFixture::count * sizeof(double));
    MomentumPredictor predictor(fixture.solver);
    expect_error(
        [&] { static_cast<void>(predictor.solve(equations, SolveControl{})); });
    fixture.check_unchanged();
  }
  for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    PredictorFixture fixture;
    fixture.op1 =
        DiagonalOperator(fixture.context, fixture.layout, {5.0, invalid, 7.0});
    fixture.operators[1] = &fixture.op1;
    fixture.solver.expected_operators = fixture.operators;
    MomentumPredictor predictor(fixture.solver);
    expect_error([&] {
      static_cast<void>(predictor.solve(fixture.equations(), SolveControl{}));
    });
    fixture.check_unchanged();
  }
}

} // namespace

int main() {
  return hundun::test::run([] {
    check_time_stencils();
    check_predictor_success_and_reports();
    check_predictor_preflight_transaction();
  });
}
