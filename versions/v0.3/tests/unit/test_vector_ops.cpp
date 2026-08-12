// SPDX-License-Identifier: Apache-2.0

#include "hundun/exec_execution.hpp"
#include "hundun/lin_vector_ops.hpp"
#include "hundun/rt_error.hpp"
#include "tests/support/exec_execution_test_access.hpp"
#include "src/lin_vector_ops_detail.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
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
using hundun::execution::ExecutionSpace;
using hundun::execution::VectorView;
using hundun::execution::test::ExecutionTestAccess;
using hundun::linear::VectorOps;
using hundun::runtime::Error;
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
  const std::string message = expect_error(std::forward<Function>(function));
  HUNDUN_CHECK(message.find(text) != std::string::npos);
}

std::vector<double> read(VectorView<const double> view) {
  std::vector<double> values;
  values.reserve(view.size());
  for (std::size_t index = 0; index < view.size(); ++index) {
    values.push_back(view[index]);
  }
  return values;
}

void set(VectorView<double> view, const std::vector<double>& values) {
  HUNDUN_CHECK(view.size() == values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    view[index] = values[index];
  }
}

class TestContext final : public ExecutionContext {
 public:
  TestContext(BackendIdentity identity, ExecutionSpace space,
              bool host_access) noexcept
      : identity_(identity), space_(space), host_access_(host_access) {}

  std::string_view backend_name() const noexcept override { return "test"; }
  BackendIdentity backend_identity() const noexcept override {
    return identity_;
  }
  ExecutionSpace space() const noexcept override { return space_; }
  bool ordered() const noexcept override { return true; }
  bool supports(ExecutionCapability capability) const noexcept override {
    return capability == ExecutionCapability::host_access && host_access_;
  }

 private:
  BackendIdentity identity_;
  ExecutionSpace space_;
  bool host_access_;
};

void test_constructor_and_empty_operations() {
  CpuReferenceContext context;
  VectorOps operations(context);
  Buffer empty(context, 0U);
  auto output = empty.view(0U, 0U);
  const Buffer& constant_empty = empty;
  auto input = constant_empty.view(0U, 0U);

  operations.fill(output, 2.0);
  operations.copy(input, output);
  operations.scale(-3.0, output);
  operations.axpy(4.0, input, output);
  operations.linear_combination(2.0, input, -1.0, input, output);

  TestContext zero_identity(0U, ExecutionSpace::host, true);
  expect_error_containing([&] { VectorOps invalid(zero_identity); },
                          "identity");
  TestContext device(71U, ExecutionSpace::device, false);
  expect_error_containing([&] { VectorOps invalid(device); }, "host");
  TestContext inaccessible(72U, ExecutionSpace::host, false);
  expect_error_containing([&] { VectorOps invalid(inaccessible); }, "host");
}

void test_bad_dot_product_sentinel_policy() {
  const double sentinel =
      hundun::linear::detail::bad_dot_product_sentinel();
  HUNDUN_CHECK(std::isinf(sentinel));
  HUNDUN_CHECK(sentinel > 0.0);
}

void test_contiguous_strided_and_aliases() {
  CpuReferenceContext context;
  VectorOps operations(context);
  Buffer x_buffer(context, 12U * sizeof(double));
  Buffer y_buffer(context, 12U * sizeof(double));
  Buffer output_buffer(context, 12U * sizeof(double));
  auto x = x_buffer.view(1U, 4U, 2U);
  auto y = y_buffer.view(2U, 4U, 2U);
  auto output = output_buffer.view(3U, 4U, 2U);
  set(x, {1.0, 2.0, 3.0, 4.0});
  set(y, {5.0, 6.0, 7.0, 8.0});

  operations.fill(output, 9.0);
  HUNDUN_CHECK(read(output) == std::vector<double>({9.0, 9.0, 9.0, 9.0}));
  operations.copy(x, output);
  HUNDUN_CHECK(read(output) == std::vector<double>({1.0, 2.0, 3.0, 4.0}));
  operations.scale(-2.0, output);
  HUNDUN_CHECK(read(output) ==
               std::vector<double>({-2.0, -4.0, -6.0, -8.0}));
  operations.axpy(0.5, x, output);
  HUNDUN_CHECK(read(output) ==
               std::vector<double>({-1.5, -3.0, -4.5, -6.0}));
  operations.linear_combination(2.0, x, -1.0, y, output);
  HUNDUN_CHECK(read(output) ==
               std::vector<double>({-3.0, -2.0, -1.0, 0.0}));

  const auto before_copy_alias = read(x);
  operations.copy(x, x);
  HUNDUN_CHECK(read(x) == before_copy_alias);
  operations.axpy(2.0, x, x);
  HUNDUN_CHECK(read(x) == std::vector<double>({3.0, 6.0, 9.0, 12.0}));
  operations.linear_combination(2.0, x, -1.0, y, x);
  HUNDUN_CHECK(read(x) == std::vector<double>({1.0, 6.0, 11.0, 16.0}));
  operations.linear_combination(0.5, y, 0.25, y, y);
  HUNDUN_CHECK(read(y) == std::vector<double>({3.75, 4.5, 5.25, 6.0}));
}

void test_rejected_views_and_transactionality() {
  CpuReferenceContext context;
  VectorOps operations(context);
  Buffer shared(context, 8U * sizeof(double));
  Buffer other(context, 8U * sizeof(double));
  auto shared_all = shared.view(0U, 8U);
  auto other_all = other.view(0U, 8U);
  set(shared_all, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0});
  set(other_all, {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0});

  auto input = shared.view(0U, 4U);
  auto partial = shared.view(1U, 4U);
  const auto shared_before = read(shared_all);
  expect_error_containing([&] { operations.copy(input, partial); }, "alias");
  expect_error_containing([&] { operations.axpy(1.0, input, partial); },
                          "alias");
  expect_error_containing(
      [&] {
        operations.linear_combination(1.0, input, 1.0, other.view(0U, 4U),
                                      partial);
      },
      "alias");
  HUNDUN_CHECK(read(shared_all) == shared_before);

  const auto other_before = read(other_all);
  expect_error_containing(
      [&] { operations.copy(input, other.view(0U, 3U)); }, "size");
  HUNDUN_CHECK(read(other_all) == other_before);

  auto metadata = ExecutionTestAccess::metadata(other_all);
  metadata.writable = false;
  auto read_only_output = ExecutionTestAccess::mutable_view(other, metadata);
  expect_error_containing(
      [&] { operations.fill(read_only_output, 1.0); }, "writable");
  HUNDUN_CHECK(read(other_all) == other_before);

  metadata = ExecutionTestAccess::metadata(other_all);
  ++metadata.backend_identity;
  auto wrong_backend = ExecutionTestAccess::mutable_view(other, metadata);
  expect_error_containing([&] { operations.fill(wrong_backend, 1.0); },
                          "backend");
  HUNDUN_CHECK(read(other_all) == other_before);

  metadata = ExecutionTestAccess::metadata(other_all);
  metadata.space = ExecutionSpace::device;
  auto wrong_space = ExecutionTestAccess::mutable_view(other, metadata);
  expect_error_containing([&] { operations.fill(wrong_space, 1.0); },
                          "host");
  HUNDUN_CHECK(read(other_all) == other_before);

  auto stale = other.view(0U, 4U);
  other.reallocate(8U * sizeof(double));
  expect_error_containing([&] { operations.fill(stale, 1.0); }, "live");

  Buffer finite(context, 4U * sizeof(double));
  Buffer destination(context, 4U * sizeof(double));
  auto finite_values = finite.view(0U, 4U);
  auto destination_values = destination.view(0U, 4U);
  set(finite_values, {1.0, 2.0, 3.0, 4.0});
  set(destination_values, {8.0, 8.0, 8.0, 8.0});
  const auto destination_before = read(destination_values);

  expect_error_containing(
      [&] {
        operations.fill(destination_values,
                        std::numeric_limits<double>::infinity());
      },
      "finite");
  HUNDUN_CHECK(read(destination_values) == destination_before);
  finite_values[2] = std::numeric_limits<double>::quiet_NaN();
  expect_error_containing(
      [&] { operations.copy(finite_values, destination_values); }, "finite");
  HUNDUN_CHECK(read(destination_values) == destination_before);
  finite_values[2] = 3.0;
  expect_error_containing(
      [&] {
        operations.scale(std::numeric_limits<double>::infinity(),
                         destination_values);
      },
      "finite");
  HUNDUN_CHECK(read(destination_values) == destination_before);
  finite_values[0] = std::numeric_limits<double>::max();
  destination_values[0] = std::numeric_limits<double>::max();
  const auto overflow_before = read(destination_values);
  expect_error_containing(
      [&] { operations.scale(2.0, destination_values); }, "finite");
  HUNDUN_CHECK(read(destination_values) == overflow_before);
  expect_error_containing(
      [&] { operations.axpy(1.0, finite_values, destination_values); },
      "finite");
  HUNDUN_CHECK(read(destination_values) == overflow_before);
  expect_error_containing(
      [&] {
        operations.linear_combination(2.0, finite_values, 1.0,
                                      finite_values, destination_values);
      },
      "finite");
  HUNDUN_CHECK(read(destination_values) == overflow_before);
}

void test_successful_operations_do_not_allocate() {
  CpuReferenceContext context;
  VectorOps operations(context);
  Buffer x_buffer(context, 4U * sizeof(double));
  Buffer y_buffer(context, 4U * sizeof(double));
  Buffer output_buffer(context, 4U * sizeof(double));
  auto x = x_buffer.view(0U, 4U);
  auto y = y_buffer.view(0U, 4U);
  auto output = output_buffer.view(0U, 4U);
  set(x, {1.0, 2.0, 3.0, 4.0});
  set(y, {5.0, 6.0, 7.0, 8.0});

  std::size_t attempts = 0U;
  {
    AllocationAttemptGuard guard;
    operations.fill(output, 1.0);
    operations.copy(x, output);
    operations.scale(2.0, output);
    operations.axpy(-1.0, y, output);
    operations.linear_combination(2.0, x, 3.0, y, output);
    attempts = guard.attempts();
  }
  HUNDUN_CHECK(attempts == 0U);
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_bad_dot_product_sentinel_policy();
    test_constructor_and_empty_operations();
    test_contiguous_strided_and_aliases();
    test_rejected_views_and_transactionality();
    test_successful_operations_do_not_allocate();
  });
}
