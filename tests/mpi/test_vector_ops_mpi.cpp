// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"
#include "hundun/linear/vector_ops.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "execution/src/execution_test_access.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>

namespace hundun::runtime::test {

class MpiContextTestAccess final {
 public:
  static void set_fp64_reduction_counters(
      MpiContext& context, Fp64ReductionCounters counters) noexcept {
    context.fp64_reduction_counters_ = counters;
  }
};

}  // namespace hundun::runtime::test

namespace {

using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::execution::VectorView;
using hundun::execution::test::ExecutionTestAccess;
using hundun::linear::DotProductPair;
using hundun::linear::VectorOps;
using hundun::runtime::Error;
using hundun::runtime::Fp64ReductionCounters;
using hundun::runtime::Fp64ReductionOperation;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::test::MpiContextTestAccess;
using hundun::test::allocation_probe::AllocationAttemptGuard;

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

void check_counter_delta(Fp64ReductionCounters before,
                         Fp64ReductionCounters after,
                         std::uint64_t calls, std::uint64_t scalars,
                         std::uint64_t bytes) {
  HUNDUN_CHECK(after.collective_calls - before.collective_calls == calls);
  HUNDUN_CHECK(after.reduced_scalars - before.reduced_scalars == scalars);
  HUNDUN_CHECK(after.logical_payload_bytes - before.logical_payload_bytes ==
               bytes);
}

void test_runtime_reductions(const MpiContext& world) {
  auto context = MpiContext::duplicate(world.comm());
  HUNDUN_CHECK(context.fp64_reduction_counters().collective_calls == 0U);

  double sum = static_cast<double>(context.rank() + 1);
  context.allreduce_fp64_in_place(&sum, 1U, Fp64ReductionOperation::sum);
  const double expected_sum =
      0.5 * static_cast<double>(context.size()) *
      static_cast<double>(context.size() + 1);
  HUNDUN_CHECK_NEAR(sum, expected_sum, 0.0);

  std::array<double, 2> maxima{
      -static_cast<double>(context.rank()),
      static_cast<double>(2 * context.rank() + 1)};
  context.allreduce_fp64_in_place(maxima.data(), maxima.size(),
                                 Fp64ReductionOperation::maximum);
  HUNDUN_CHECK_NEAR(maxima[0], 0.0, 0.0);
  HUNDUN_CHECK_NEAR(maxima[1],
                    static_cast<double>(2 * (context.size() - 1) + 1), 0.0);
  const auto after_nonempty = context.fp64_reduction_counters();
  HUNDUN_CHECK(after_nonempty.collective_calls == 2U);
  HUNDUN_CHECK(after_nonempty.reduced_scalars == 3U);
  HUNDUN_CHECK(after_nonempty.logical_payload_bytes == 3U * sizeof(double));

  context.allreduce_fp64_in_place(nullptr, 0U, Fp64ReductionOperation::sum);
  const auto after_empty = context.fp64_reduction_counters();
  check_counter_delta(after_nonempty, after_empty, 0U, 0U, 0U);

  expect_error_containing(
      [&] {
        context.allreduce_fp64_in_place(nullptr, 1U,
                                       Fp64ReductionOperation::sum);
      },
      "pointer");
  double dummy = 0.0;
  expect_error_containing(
      [&] {
        context.allreduce_fp64_in_place(
            &dummy, static_cast<std::size_t>(INT_MAX) + 1U,
            Fp64ReductionOperation::sum);
      },
      "INT_MAX");
  expect_error_containing(
      [&] {
        context.allreduce_fp64_in_place(
            &dummy, 1U, static_cast<Fp64ReductionOperation>(99));
      },
      "operation");
  check_counter_delta(after_empty, context.fp64_reduction_counters(), 0U, 0U,
                      0U);
}

void test_counter_overflow_and_moves(const MpiContext& world) {
  const auto maximum = std::numeric_limits<std::uint64_t>::max();
  double value = 1.0;

  auto calls = MpiContext::duplicate(world.comm());
  MpiContextTestAccess::set_fp64_reduction_counters(
      calls, Fp64ReductionCounters{maximum, 0U, 0U});
  expect_error_containing(
      [&] {
        calls.allreduce_fp64_in_place(&value, 1U,
                                      Fp64ReductionOperation::sum);
      },
      "counter");
  HUNDUN_CHECK(calls.fp64_reduction_counters().collective_calls == maximum);

  auto scalars = MpiContext::duplicate(world.comm());
  MpiContextTestAccess::set_fp64_reduction_counters(
      scalars, Fp64ReductionCounters{0U, maximum, 0U});
  expect_error_containing(
      [&] {
        scalars.allreduce_fp64_in_place(&value, 1U,
                                        Fp64ReductionOperation::sum);
      },
      "counter");
  HUNDUN_CHECK(scalars.fp64_reduction_counters().reduced_scalars == maximum);

  auto bytes = MpiContext::duplicate(world.comm());
  MpiContextTestAccess::set_fp64_reduction_counters(
      bytes, Fp64ReductionCounters{0U, 0U, maximum});
  expect_error_containing(
      [&] {
        bytes.allreduce_fp64_in_place(&value, 1U,
                                      Fp64ReductionOperation::sum);
      },
      "counter");
  HUNDUN_CHECK(bytes.fp64_reduction_counters().logical_payload_bytes ==
               maximum);

  auto source = MpiContext::duplicate(world.comm());
  source.allreduce_fp64_in_place(&value, 1U, Fp64ReductionOperation::sum);
  MpiContext moved(std::move(source));
  HUNDUN_CHECK(source.comm() == MPI_COMM_NULL);
  HUNDUN_CHECK(source.fp64_reduction_counters().collective_calls == 0U);
  HUNDUN_CHECK(moved.fp64_reduction_counters().collective_calls == 1U);
  expect_error_containing(
      [&] {
        source.allreduce_fp64_in_place(&value, 1U,
                                       Fp64ReductionOperation::sum);
      },
      "empty");

  auto replacement = MpiContext::duplicate(world.comm());
  std::array<double, 2> replacement_values{1.0, 2.0};
  replacement.allreduce_fp64_in_place(
      replacement_values.data(), replacement_values.size(),
      Fp64ReductionOperation::maximum);
  auto destination = MpiContext::duplicate(world.comm());
  destination.allreduce_fp64_in_place(&value, 1U,
                                      Fp64ReductionOperation::sum);
  destination = std::move(replacement);
  HUNDUN_CHECK(replacement.comm() == MPI_COMM_NULL);
  HUNDUN_CHECK(replacement.fp64_reduction_counters().collective_calls == 0U);
  HUNDUN_CHECK(destination.fp64_reduction_counters().collective_calls == 1U);
  HUNDUN_CHECK(destination.fp64_reduction_counters().reduced_scalars == 2U);
  HUNDUN_CHECK(destination.fp64_reduction_counters().logical_payload_bytes ==
               2U * sizeof(double));
}

void fill_norm_values(VectorView<double> values, int rank, double factor) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index] = factor *
                    static_cast<double>(rank + 1 + static_cast<int>(index));
  }
}

double expected_norm(int ranks, double factor) {
  double sum = 0.0;
  for (int rank = 0; rank < ranks; ++rank) {
    for (int index = 0; index < rank + 1; ++index) {
      const double value = static_cast<double>(rank + 1 + index);
      sum += value * value;
    }
  }
  return factor * std::sqrt(sum);
}

void test_global_norm(const MpiContext& world) {
  CpuReferenceContext execution;
  VectorOps operations(execution);
  const std::size_t local_count = static_cast<std::size_t>(world.rank() + 1);
  Buffer values_buffer(execution, local_count * sizeof(double));
  auto values = values_buffer.view(0U, local_count);
  fill_norm_values(values, world.rank(), 1.0e200);

  const auto before = world.fp64_reduction_counters();
  std::size_t attempts = 0U;
  double actual = 0.0;
  {
    AllocationAttemptGuard guard;
    actual = operations.norm(values, world);
    attempts = guard.attempts();
  }
  HUNDUN_CHECK(attempts == 0U);
  const double expected = expected_norm(world.size(), 1.0e200);
  HUNDUN_CHECK(std::abs(actual - expected) <=
               32.0 * std::numeric_limits<double>::epsilon() * expected);
  check_counter_delta(before, world.fp64_reduction_counters(), 2U, 2U,
                      2U * sizeof(double));

  fill_norm_values(values, world.rank(), 0.0);
  const auto before_zero = world.fp64_reduction_counters();
  HUNDUN_CHECK_NEAR(operations.norm(values, world), 0.0, 0.0);
  check_counter_delta(before_zero, world.fp64_reduction_counters(), 2U, 2U,
                      2U * sizeof(double));

  Buffer empty_buffer(execution, 0U);
  const auto before_empty = world.fp64_reduction_counters();
  HUNDUN_CHECK_NEAR(operations.norm(empty_buffer.view(0U, 0U), world), 0.0,
                    0.0);
  check_counter_delta(before_empty, world.fp64_reduction_counters(), 2U, 2U,
                      2U * sizeof(double));

  fill_norm_values(values, world.rank(), 1.0);
  if (world.rank() == world.size() - 1) {
    values[0] = std::numeric_limits<double>::quiet_NaN();
  }
  const auto before_nonfinite = world.fp64_reduction_counters();
  expect_error_containing([&] { static_cast<void>(operations.norm(values, world)); },
                          "finite");
  check_counter_delta(before_nonfinite, world.fp64_reduction_counters(), 1U,
                      1U, sizeof(double));

  fill_norm_values(values, world.rank(), 1.0);
  auto metadata = ExecutionTestAccess::metadata(values);
  if (world.rank() == world.size() - 1) {
    ++metadata.epoch;
  }
  const auto invalid = ExecutionTestAccess::const_view(values_buffer, metadata);
  const auto before_invalid = world.fp64_reduction_counters();
  expect_error_containing(
      [&] { static_cast<void>(operations.norm(invalid, world)); }, "finite");
  check_counter_delta(before_invalid, world.fp64_reduction_counters(), 1U, 1U,
                      sizeof(double));
}

double expected_dot(int ranks, int pair) {
  double sum = 0.0;
  for (int rank = 0; rank < ranks; ++rank) {
    const int count = rank + 4;
    for (int index = 0; index < count; ++index) {
      const double a = static_cast<double>(rank + 1 + index);
      const double b = static_cast<double>(index + 2);
      if (pair == 0) {
        sum += a * b;
      } else if (pair == 1) {
        sum += (2.0 * a) * (b - 1.0);
      } else {
        sum += (a - 0.5) * (-b);
      }
    }
  }
  return sum;
}

void fill_dot_inputs(VectorView<double> a, VectorView<double> b,
                     VectorView<double> c, VectorView<double> d,
                     VectorView<double> e, VectorView<double> f, int rank) {
  for (std::size_t index = 0; index < a.size(); ++index) {
    const double av =
        static_cast<double>(rank + 1 + static_cast<int>(index));
    const double bv = static_cast<double>(index + 2U);
    a[index] = av;
    b[index] = bv;
    c[index] = 2.0 * av;
    d[index] = bv - 1.0;
    e[index] = av - 0.5;
    f[index] = -bv;
  }
}

void test_dot_batch(const MpiContext& world) {
  CpuReferenceContext execution;
  VectorOps operations(execution);
  const std::size_t count = static_cast<std::size_t>(world.rank() + 4);
  Buffer a_buffer(execution, count * sizeof(double));
  Buffer b_buffer(execution, count * sizeof(double));
  Buffer c_buffer(execution, count * sizeof(double));
  Buffer d_buffer(execution, count * sizeof(double));
  Buffer e_buffer(execution, count * sizeof(double));
  Buffer f_buffer(execution, count * sizeof(double));
  Buffer results_buffer(execution, 8U * sizeof(double));
  auto a = a_buffer.view(0U, count);
  auto b = b_buffer.view(0U, count);
  auto c = c_buffer.view(0U, count);
  auto d = d_buffer.view(0U, count);
  auto e = e_buffer.view(0U, count);
  auto f = f_buffer.view(0U, count);
  auto results = results_buffer.view(0U, 3U);
  fill_dot_inputs(a, b, c, d, e, f, world.rank());
  const std::array<DotProductPair, 3> pairs{
      DotProductPair{a, b}, DotProductPair{c, d}, DotProductPair{e, f}};

  const auto before = world.fp64_reduction_counters();
  std::size_t attempts = 0U;
  {
    AllocationAttemptGuard guard;
    operations.dot_batch(pairs.data(), pairs.size(), results, world);
    attempts = guard.attempts();
  }
  HUNDUN_CHECK(attempts == 0U);
  for (int pair = 0; pair < 3; ++pair) {
    HUNDUN_CHECK_NEAR(results[static_cast<std::size_t>(pair)],
                      expected_dot(world.size(), pair), 1.0e-12);
  }
  check_counter_delta(before, world.fp64_reduction_counters(), 1U, 3U,
                      3U * sizeof(double));

  const auto before_empty = world.fp64_reduction_counters();
  operations.dot_batch(nullptr, 0U, results_buffer.view(0U, 0U), world);
  check_counter_delta(before_empty, world.fp64_reduction_counters(), 0U, 0U,
                      0U);

  auto check_preflight_rejection = [&](auto&& call, const std::string& text) {
    const auto before_failure = world.fp64_reduction_counters();
    expect_error_containing(std::forward<decltype(call)>(call), text);
    check_counter_delta(before_failure, world.fp64_reduction_counters(), 0U,
                        0U, 0U);
  };
  check_preflight_rejection(
      [&] { operations.dot_batch(nullptr, 1U, results_buffer.view(0U, 1U), world); },
      "pointer");
  check_preflight_rejection(
      [&] { operations.dot_batch(pairs.data(), 3U, results_buffer.view(0U, 2U), world); },
      "size");
  check_preflight_rejection(
      [&] { operations.dot_batch(pairs.data(), 3U, results_buffer.view(0U, 3U, 2U), world); },
      "stride");
  check_preflight_rejection(
      [&] { operations.dot_batch(pairs.data(), 3U, a_buffer.view(0U, 3U), world); },
      "distinct");
  check_preflight_rejection(
      [&] {
        operations.dot_batch(
            pairs.data(), static_cast<std::size_t>(INT_MAX) + 1U,
            results_buffer.view(0U, 0U), world);
      },
      "INT_MAX");

  auto bad_pairs = pairs;
  if (world.rank() == world.size() - 1) {
    bad_pairs[1] = DotProductPair{c, d_buffer.view(0U, count - 1U)};
  }
  const auto before_bad_size = world.fp64_reduction_counters();
  expect_error_containing(
      [&] {
        operations.dot_batch(bad_pairs.data(), bad_pairs.size(), results,
                             world);
      },
      "finite");
  check_counter_delta(before_bad_size, world.fp64_reduction_counters(), 1U,
                      3U, 3U * sizeof(double));

  fill_dot_inputs(a, b, c, d, e, f, world.rank());
  if (world.rank() == world.size() - 1) {
    e[0] = std::numeric_limits<double>::infinity();
  }
  const auto before_bad_value = world.fp64_reduction_counters();
  expect_error_containing(
      [&] { operations.dot_batch(pairs.data(), pairs.size(), results, world); },
      "finite");
  check_counter_delta(before_bad_value, world.fp64_reduction_counters(), 1U,
                      3U, 3U * sizeof(double));
}

}  // namespace

int main(int argc, char** argv) {
  MpiEnvironment environment(argc, argv);
  auto world = MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    test_runtime_reductions(world);
    test_counter_overflow_and_moves(world);
    test_global_norm(world);
    test_dot_batch(world);
  });
}
