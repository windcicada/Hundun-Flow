// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"
#include "hundun/linear/linear_system.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/runtime/error.hpp"
#include "execution/src/execution_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
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
using hundun::execution::test::TestViewMetadata;
using hundun::linear::IdentityPreconditioner;
using hundun::linear::JacobiPreconditioner;
using hundun::linear::LinearOperator;
using hundun::linear::LinearSolver;
using hundun::linear::Preconditioner;
using hundun::linear::SolveControl;
using hundun::linear::SolveReport;
using hundun::linear::VectorLayout;
using hundun::runtime::Error;

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

VectorLayout make_layout(std::size_t owned, std::size_t ghosts = 0U,
                         hundun::mesh::GlobalCellId first = 1U) {
  std::vector<hundun::mesh::GlobalCellId> ids;
  ids.reserve(owned + ghosts);
  for (std::size_t index = 0; index < owned + ghosts; ++index) {
    ids.push_back(first + static_cast<hundun::mesh::GlobalCellId>(index));
  }
  return VectorLayout(owned, std::move(ids));
}

void set_values(VectorView<double> view, const std::vector<double>& values) {
  HUNDUN_CHECK(view.size() == values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    view[index] = values[index];
  }
}

std::vector<double> read_values(VectorView<const double> view) {
  std::vector<double> result;
  result.reserve(view.size());
  for (std::size_t index = 0; index < view.size(); ++index) {
    result.push_back(view[index]);
  }
  return result;
}

bool bitwise_equal(const std::vector<double>& left,
                   const std::vector<double>& right) {
  return left.size() == right.size() &&
         std::memcmp(left.data(), right.data(),
                     left.size() * sizeof(double)) == 0;
}

class ConfigurableContext final : public ExecutionContext {
 public:
  ConfigurableContext(BackendIdentity identity, ExecutionSpace space,
                      bool buffer, bool host, bool transfer) noexcept
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

class FakeOperator final : public LinearOperator {
 public:
  FakeOperator(ExecutionContext& context, VectorLayout domain,
               VectorLayout range, std::uint64_t revision,
               std::vector<double> diagonal = {})
      : context_(&context),
        domain_(std::move(domain)),
        range_(std::move(range)),
        revision_(revision),
        diagonal_(std::move(diagonal)) {}

  ~FakeOperator() noexcept override {
    if (destroyed_ != nullptr) {
      *destroyed_ = true;
    }
  }

  VectorLayout domain_layout() const override { return domain_; }
  VectorLayout range_layout() const override { return range_; }
  const ExecutionContext& context() const override { return *context_; }
  std::uint64_t revision() const override { return revision_; }

  ExecutionEvent apply(VectorView<const double> x,
                       VectorView<double> y) const override {
    if (x.size() != domain_.owned_count() ||
        y.size() != range_.owned_count()) {
      throw Error("fake operator view size mismatch");
    }
    for (std::size_t index = 0; index < y.size(); ++index) {
      y[index] = x[index];
    }
    return ExecutionEvent::completed();
  }

  bool has_diagonal() const override { return has_diagonal_; }

  ExecutionEvent diagonal(VectorView<double> output) const override {
    ++diagonal_calls_;
    if (throw_on_diagonal_) {
      throw Error("fake diagonal throw");
    }
    if (output.size() != range_.owned_count() ||
        diagonal_.size() != range_.owned_count()) {
      throw Error("fake diagonal size mismatch");
    }
    for (std::size_t index = 0; index < diagonal_.size(); ++index) {
      output[index] = diagonal_[index];
    }
    if (fail_diagonal_event_) {
      auto event = ExecutionTestAccess::pending_event({});
      ExecutionTestAccess::complete_error(event, "fake diagonal event error");
      return event;
    }
    return ExecutionEvent::completed();
  }

  void set_revision(std::uint64_t revision) noexcept { revision_ = revision; }
  void set_layouts(VectorLayout domain, VectorLayout range) {
    domain_ = std::move(domain);
    range_ = std::move(range);
  }
  void set_context(ExecutionContext& context) noexcept { context_ = &context; }
  void set_diagonal(std::vector<double> diagonal) {
    diagonal_ = std::move(diagonal);
  }
  void set_has_diagonal(bool available) noexcept {
    has_diagonal_ = available;
  }
  void set_throw_on_diagonal(bool enabled) noexcept {
    throw_on_diagonal_ = enabled;
  }
  void set_fail_diagonal_event(bool enabled) noexcept {
    fail_diagonal_event_ = enabled;
  }
  std::size_t diagonal_calls() const noexcept { return diagonal_calls_; }
  void observe_destruction(bool& destroyed) noexcept { destroyed_ = &destroyed; }

 private:
  ExecutionContext* context_;
  VectorLayout domain_;
  VectorLayout range_;
  std::uint64_t revision_;
  std::vector<double> diagonal_;
  bool has_diagonal_{true};
  bool throw_on_diagonal_{false};
  bool fail_diagonal_event_{false};
  mutable std::size_t diagonal_calls_{0U};
  bool* destroyed_{nullptr};
};

class ReusedFakeOperatorStorage final {
 public:
  ReusedFakeOperatorStorage() = default;
  ~ReusedFakeOperatorStorage() noexcept { destroy(); }
  ReusedFakeOperatorStorage(const ReusedFakeOperatorStorage&) = delete;
  ReusedFakeOperatorStorage& operator=(const ReusedFakeOperatorStorage&) =
      delete;

  FakeOperator& emplace(ExecutionContext& context, VectorLayout domain,
                        VectorLayout range, std::uint64_t revision,
                        std::vector<double> diagonal) {
    if (live_) {
      throw std::logic_error("placement operator storage is already live");
    }
    auto* created = ::new (static_cast<void*>(storage_)) FakeOperator(
        context, std::move(domain), std::move(range), revision,
        std::move(diagonal));
    live_ = true;
    return *created;
  }

  void destroy() noexcept {
    if (live_) {
      current()->~FakeOperator();
      live_ = false;
    }
  }

 private:
  FakeOperator* current() noexcept {
    return std::launder(reinterpret_cast<FakeOperator*>(storage_));
  }

  alignas(FakeOperator) std::byte storage_[sizeof(FakeOperator)];
  bool live_{false};
};

class ProbePreconditioner final : public Preconditioner {
 public:
  explicit ProbePreconditioner(bool& destroyed) : destroyed_(destroyed) {}
  ~ProbePreconditioner() noexcept override { destroyed_ = true; }
  void update(const LinearOperator&, std::uint64_t) override {}
  ExecutionEvent apply(VectorView<const double>,
                       VectorView<double>) const override {
    return ExecutionEvent::completed();
  }

 private:
  bool& destroyed_;
};

class ProbeSolver final : public LinearSolver {
 public:
  explicit ProbeSolver(bool& destroyed) : destroyed_(destroyed) {}
  ~ProbeSolver() noexcept override { destroyed_ = true; }
  SolveReport solve(const LinearOperator&, Preconditioner&,
                    VectorView<const double>, VectorView<double>,
                    const SolveControl&) const override {
    return {};
  }

 private:
  bool& destroyed_;
};

void check_constructor_rejected(ExecutionContext& context,
                                const std::string& text) {
  expect_error_containing([&] { IdentityPreconditioner value(context); },
                          text);
  expect_error_containing([&] { JacobiPreconditioner value(context); }, text);
}

void test_completed_event_factory() {
  auto event = hundun::execution::ExecutionEvent::completed();
  HUNDUN_CHECK(event.ready());
  event.wait();
  event.wait();

  auto copied = event;
  HUNDUN_CHECK(copied.ready());
  copied.wait();

  auto moved = std::move(event);
  HUNDUN_CHECK(moved.ready());
  moved.wait();
  HUNDUN_CHECK(!event.ready());
  expect_error_containing([&] { event.wait(); }, "moved-from");
}

void test_completed_event_allocation_failure_is_one_shot() {
  CpuReferenceContext context;
  Buffer source(context, sizeof(double));
  Buffer destination(context, sizeof(double));
  source.view(0U, 1U)[0] = 3.0;
  destination.view(0U, 1U)[0] = -1.0;

  ExecutionTestAccess::fail_next_completed_event_allocation();
  hundun::execution::transfer(source.view(0U, 1U),
                              destination.view(0U, 1U), context)
      .wait();
  HUNDUN_CHECK(destination.view(0U, 1U)[0] == 3.0);
  expect_error_containing(
      [] { static_cast<void>(ExecutionEvent::completed()); },
      "completed event allocation");

  auto recovered = ExecutionEvent::completed();
  HUNDUN_CHECK(recovered.ready());
  recovered.wait();
}

void test_virtual_destruction() {
  CpuReferenceContext context;
  const auto layout = make_layout(1U);

  bool operator_destroyed = false;
  {
    auto concrete =
        std::make_unique<FakeOperator>(context, layout, layout, 1U,
                                       std::vector<double>{1.0});
    concrete->observe_destruction(operator_destroyed);
    std::unique_ptr<LinearOperator> base = std::move(concrete);
  }
  HUNDUN_CHECK(operator_destroyed);

  bool preconditioner_destroyed = false;
  {
    std::unique_ptr<Preconditioner> base =
        std::make_unique<ProbePreconditioner>(preconditioner_destroyed);
  }
  HUNDUN_CHECK(preconditioner_destroyed);

  bool solver_destroyed = false;
  {
    std::unique_ptr<LinearSolver> base =
        std::make_unique<ProbeSolver>(solver_destroyed);
  }
  HUNDUN_CHECK(solver_destroyed);
}

void test_context_and_update_validation() {
  CpuReferenceContext context;
  ConfigurableContext zero(0U, ExecutionSpace::host, true, true, true);
  ConfigurableContext device(91U, ExecutionSpace::device, true, false, true);
  ConfigurableContext no_buffer(92U, ExecutionSpace::host, false, true, true);
  ConfigurableContext no_host(93U, ExecutionSpace::host, true, false, true);
  ConfigurableContext no_transfer(94U, ExecutionSpace::host, true, true,
                                  false);
  check_constructor_rejected(zero, "identity");
  check_constructor_rejected(device, "host");
  check_constructor_rejected(no_buffer, "buffer");
  check_constructor_rejected(no_host, "host");
  check_constructor_rejected(no_transfer, "transfer");

  const auto layout = make_layout(2U, 1U);
  FakeOperator linear_operator(context, layout, layout, 5U, {2.0, 4.0});
  IdentityPreconditioner identity(context);

  Buffer input(context, 2U * sizeof(double));
  Buffer output(context, 2U * sizeof(double));
  set_values(input.view(0U, 2U), {1.0, 2.0});
  set_values(output.view(0U, 2U), {8.0, 8.0});
  const auto unchanged = read_values(output.view(0U, 2U));
  expect_error_containing(
      [&] { identity.apply(input.view(0U, 2U), output.view(0U, 2U)); },
      "update");
  HUNDUN_CHECK(read_values(output.view(0U, 2U)) == unchanged);

  expect_error_containing([&] { identity.update(linear_operator, 4U); },
                          "revision");
  FakeOperator nonsquare(context, layout, make_layout(2U, 1U, 50U), 1U,
                         {1.0, 1.0});
  expect_error_containing([&] { identity.update(nonsquare, 1U); }, "square");
  CpuReferenceContext same_backend_different_object;
  FakeOperator wrong_context(same_backend_different_object, layout, layout,
                             1U, {1.0, 1.0});
  expect_error_containing([&] { identity.update(wrong_context, 1U); },
                          "context");

  identity.update(linear_operator, 5U);
  linear_operator.set_revision(4U);
  expect_error_containing([&] { identity.update(linear_operator, 4U); },
                          "decrease");
  linear_operator.set_revision(5U);
  identity.apply(input.view(0U, 2U), output.view(0U, 2U)).wait();
}

void test_current_operator_revision_lineages() {
  CpuReferenceContext context;
  const auto layout = make_layout(2U);
  FakeOperator identity_first(context, layout, layout, 5U, {1.0, 1.0});
  FakeOperator identity_second(context, layout, layout, 2U, {1.0, 1.0});
  IdentityPreconditioner identity(context);
  identity.update(identity_first, 5U);
  identity.update(identity_second, 2U);
  identity_second.set_revision(1U);
  expect_error_containing([&] { identity.update(identity_second, 1U); },
                          "decrease");
  identity_second.set_revision(2U);

  Buffer input(context, 2U * sizeof(double));
  Buffer output(context, 2U * sizeof(double));
  set_values(input.view(0U, 2U), {4.0, -6.0});
  set_values(output.view(0U, 2U), {9.0, 9.0});
  identity.apply(input.view(0U, 2U), output.view(0U, 2U)).wait();
  HUNDUN_CHECK(read_values(output.view(0U, 2U)) ==
               std::vector<double>({4.0, -6.0}));

  FakeOperator jacobi_first(context, layout, layout, 8U, {2.0, 4.0});
  FakeOperator jacobi_second(context, layout, layout, 3U, {4.0, -2.0});
  JacobiPreconditioner jacobi(context);
  jacobi.update(jacobi_first, 8U);
  jacobi.update(jacobi_second, 3U);
  HUNDUN_CHECK(jacobi_second.diagonal_calls() == 1U);
  jacobi_second.set_revision(2U);
  expect_error_containing([&] { jacobi.update(jacobi_second, 2U); },
                          "decrease");
  HUNDUN_CHECK(jacobi_second.diagonal_calls() == 1U);
  jacobi_second.set_revision(3U);

  set_values(input.view(0U, 2U), {8.0, 6.0});
  set_values(output.view(0U, 2U), {-7.0, -7.0});
  jacobi.apply(input.view(0U, 2U), output.view(0U, 2U)).wait();
  HUNDUN_CHECK(read_values(output.view(0U, 2U)) ==
               std::vector<double>({2.0, -3.0}));
}

void test_reused_operator_address_starts_fresh_lineage() {
  CpuReferenceContext context;
  const auto layout = make_layout(2U);

  {
    ReusedFakeOperatorStorage reused_storage;
    FakeOperator distinct(context, layout, layout, 2U, {1.0, 1.0});
    IdentityPreconditioner identity(context);

    auto& first = reused_storage.emplace(context, layout, layout, 9U,
                                         {1.0, 1.0});
    identity.update(first, 9U);
    identity.update(distinct, 2U);
    reused_storage.destroy();

    auto& replacement = reused_storage.emplace(context, layout, layout, 1U,
                                               {1.0, 1.0});
    identity.update(replacement, 1U);

    Buffer input(context, 2U * sizeof(double));
    Buffer output(context, 2U * sizeof(double));
    set_values(input.view(0U, 2U), {5.0, -7.0});
    set_values(output.view(0U, 2U), {11.0, 11.0});
    identity.apply(input.view(0U, 2U), output.view(0U, 2U)).wait();
    HUNDUN_CHECK(read_values(output.view(0U, 2U)) ==
                 std::vector<double>({5.0, -7.0}));
  }

  {
    ReusedFakeOperatorStorage reused_storage;
    FakeOperator distinct(context, layout, layout, 3U, {4.0, 5.0});
    JacobiPreconditioner jacobi(context);

    auto& first = reused_storage.emplace(context, layout, layout, 10U,
                                         {2.0, 2.0});
    jacobi.update(first, 10U);
    jacobi.update(distinct, 3U);
    reused_storage.destroy();

    auto& replacement = reused_storage.emplace(context, layout, layout, 1U,
                                               {2.0, -4.0});
    jacobi.update(replacement, 1U);
    HUNDUN_CHECK(replacement.diagonal_calls() == 1U);

    Buffer residual(context, 2U * sizeof(double));
    Buffer correction(context, 2U * sizeof(double));
    set_values(residual.view(0U, 2U), {8.0, 12.0});
    set_values(correction.view(0U, 2U), {-9.0, -9.0});
    jacobi.apply(residual.view(0U, 2U), correction.view(0U, 2U)).wait();
    HUNDUN_CHECK(read_values(correction.view(0U, 2U)) ==
                 std::vector<double>({4.0, -3.0}));
  }
}

void test_completed_event_failure_preserves_apply_outputs() {
  CpuReferenceContext context;
  const auto layout = make_layout(3U);
  FakeOperator linear_operator(context, layout, layout, 4U,
                               {2.0, -4.0, 0.5});
  IdentityPreconditioner identity(context);
  JacobiPreconditioner jacobi(context);
  identity.update(linear_operator, 4U);
  jacobi.update(linear_operator, 4U);

  Buffer input(context, 3U * sizeof(double));
  Buffer output(context, 3U * sizeof(double));
  set_values(input.view(0U, 3U), {4.0, 8.0, -1.0});

  set_values(output.view(0U, 3U), {-0.0, 17.0, -23.0});
  const auto identity_before = read_values(output.view(0U, 3U));
  ExecutionTestAccess::fail_next_completed_event_allocation();
  expect_error_containing(
      [&] { identity.apply(input.view(0U, 3U), output.view(0U, 3U)); },
      "completed event allocation");
  HUNDUN_CHECK(bitwise_equal(read_values(output.view(0U, 3U)),
                             identity_before));
  ExecutionEvent::completed().wait();

  set_values(output.view(0U, 3U), {-0.0, -31.0, 47.0});
  const auto jacobi_before = read_values(output.view(0U, 3U));
  ExecutionTestAccess::fail_next_completed_event_allocation();
  expect_error_containing(
      [&] { jacobi.apply(input.view(0U, 3U), output.view(0U, 3U)); },
      "completed event allocation");
  HUNDUN_CHECK(bitwise_equal(read_values(output.view(0U, 3U)),
                             jacobi_before));
  ExecutionEvent::completed().wait();
}

void test_identity_apply_contract() {
  CpuReferenceContext context;
  CpuReferenceContext foreign_context;
  const auto layout = make_layout(3U, 2U);
  FakeOperator linear_operator(context, layout, layout, 7U,
                               {1.0, 1.0, 1.0});
  IdentityPreconditioner identity(context);
  identity.update(linear_operator, 7U);

  Buffer input(context, 10U * sizeof(double));
  Buffer output(context, 10U * sizeof(double));
  auto strided_input = input.view(1U, 3U, 2U);
  auto strided_output = output.view(2U, 3U, 2U);
  set_values(strided_input, {2.0, -3.0, 4.5});
  set_values(strided_output, {9.0, 9.0, 9.0});
  auto event = identity.apply(strided_input, strided_output);
  HUNDUN_CHECK(event.ready());
  event.wait();
  HUNDUN_CHECK(read_values(strided_output) ==
               std::vector<double>({2.0, -3.0, 4.5}));

  set_values(strided_input, {6.0, 7.0, 8.0});
  identity.apply(strided_input, strided_input).wait();
  HUNDUN_CHECK(read_values(strided_input) ==
               std::vector<double>({6.0, 7.0, 8.0}));

  Buffer shared(context, 6U * sizeof(double));
  set_values(shared.view(0U, 6U), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  const auto shared_before = read_values(shared.view(0U, 6U));
  expect_error_containing(
      [&] { identity.apply(shared.view(0U, 3U), shared.view(1U, 3U)); },
      "alias");
  HUNDUN_CHECK(read_values(shared.view(0U, 6U)) == shared_before);

  set_values(strided_output, {11.0, 12.0, 13.0});
  auto require_output_unchanged = [&](auto&& operation,
                                      const std::string& text) {
    const auto before = read_values(strided_output);
    expect_error_containing(std::forward<decltype(operation)>(operation),
                            text);
    HUNDUN_CHECK(read_values(strided_output) == before);
  };

  require_output_unchanged(
      [&] { identity.apply(input.view(0U, 2U), strided_output); }, "size");

  Buffer stale_owner(context, 3U * sizeof(double));
  auto stale = stale_owner.view(0U, 3U);
  stale_owner.reallocate(3U * sizeof(double));
  require_output_unchanged([&] { identity.apply(stale, strided_output); },
                           "live");

  auto metadata = ExecutionTestAccess::metadata(strided_input);
  ++metadata.backend_identity;
  auto foreign = ExecutionTestAccess::const_view(input, metadata);
  require_output_unchanged([&] { identity.apply(foreign, strided_output); },
                           "backend");

  auto output_metadata = ExecutionTestAccess::metadata(strided_output);
  output_metadata.writable = false;
  auto non_writable = ExecutionTestAccess::mutable_view(output,
                                                         output_metadata);
  const auto output_before_nonwritable = read_values(strided_output);
  expect_error_containing(
      [&] { identity.apply(strided_input, non_writable); }, "writable");
  HUNDUN_CHECK(read_values(strided_output) == output_before_nonwritable);

  strided_input[1] = std::numeric_limits<double>::quiet_NaN();
  require_output_unchanged(
      [&] { identity.apply(strided_input, strided_output); }, "finite");
  strided_input[1] = 7.0;

  linear_operator.set_revision(8U);
  require_output_unchanged(
      [&] { identity.apply(strided_input, strided_output); }, "revision");
  linear_operator.set_revision(7U);
  linear_operator.set_layouts(make_layout(3U, 2U, 100U),
                              make_layout(3U, 2U, 100U));
  require_output_unchanged(
      [&] { identity.apply(strided_input, strided_output); }, "layout");
  linear_operator.set_layouts(layout, layout);
  linear_operator.set_context(foreign_context);
  require_output_unchanged(
      [&] { identity.apply(strided_input, strided_output); }, "context");
  linear_operator.set_context(context);

  ExecutionTestAccess::fail_next_allocation();
  identity.apply(strided_input, strided_output).wait();
  expect_error_containing(
      [&] { Buffer allocation_probe(context, sizeof(double)); },
      "allocation");

  const VectorLayout zero_layout;
  FakeOperator zero_operator(context, zero_layout, zero_layout, 1U, {});
  identity.update(zero_operator, 1U);
  Buffer empty(context, 0U);
  identity.apply(empty.view(0U, 0U), empty.view(0U, 0U)).wait();
}

void test_jacobi_apply_and_cache() {
  CpuReferenceContext context;
  const auto layout = make_layout(3U, 2U);
  FakeOperator linear_operator(context, layout, layout, 10U,
                               {2.0, -4.0, 0.5});
  JacobiPreconditioner jacobi(context);
  jacobi.update(linear_operator, 10U);
  HUNDUN_CHECK(linear_operator.diagonal_calls() == 1U);
  jacobi.update(linear_operator, 10U);
  HUNDUN_CHECK(linear_operator.diagonal_calls() == 1U);

  Buffer residual(context, 10U * sizeof(double));
  Buffer correction(context, 10U * sizeof(double));
  auto r = residual.view(1U, 3U, 2U);
  auto z = correction.view(2U, 3U, 2U);
  set_values(r, {4.0, 8.0, -1.0});
  set_values(z, {99.0, 99.0, 99.0});
  auto event = jacobi.apply(r, z);
  HUNDUN_CHECK(event.ready());
  event.wait();
  HUNDUN_CHECK(read_values(z) == std::vector<double>({2.0, -2.0, -2.0}));

  set_values(r, {6.0, -8.0, 2.0});
  jacobi.apply(r, r).wait();
  HUNDUN_CHECK(read_values(r) == std::vector<double>({3.0, 2.0, 4.0}));

  linear_operator.set_revision(11U);
  linear_operator.set_diagonal({4.0, 2.0, -0.5});
  jacobi.update(linear_operator, 11U);
  HUNDUN_CHECK(linear_operator.diagonal_calls() == 2U);
  set_values(r, {8.0, 6.0, 1.0});
  jacobi.apply(r, z).wait();
  HUNDUN_CHECK(read_values(z) == std::vector<double>({2.0, 3.0, -2.0}));

  const auto changed_layout = make_layout(2U, 1U, 100U);
  linear_operator.set_layouts(changed_layout, changed_layout);
  linear_operator.set_diagonal({8.0, -2.0});
  jacobi.update(linear_operator, 11U);
  HUNDUN_CHECK(linear_operator.diagonal_calls() == 3U);

  FakeOperator distinct(context, changed_layout, changed_layout, 11U,
                        {2.0, 5.0});
  jacobi.update(distinct, 11U);
  HUNDUN_CHECK(distinct.diagonal_calls() == 1U);
  jacobi.update(distinct, 11U);
  HUNDUN_CHECK(distinct.diagonal_calls() == 1U);

  Buffer two_input(context, 2U * sizeof(double));
  Buffer two_output(context, 2U * sizeof(double));
  set_values(two_input.view(0U, 2U), {4.0, 15.0});
  set_values(two_output.view(0U, 2U), {7.0, 7.0});
  distinct.set_revision(12U);
  const auto before_stale = read_values(two_output.view(0U, 2U));
  expect_error_containing(
      [&] { jacobi.apply(two_input.view(0U, 2U), two_output.view(0U, 2U)); },
      "revision");
  HUNDUN_CHECK(read_values(two_output.view(0U, 2U)) == before_stale);
  distinct.set_revision(11U);

  ExecutionTestAccess::fail_next_allocation();
  jacobi.apply(two_input.view(0U, 2U), two_output.view(0U, 2U)).wait();
  expect_error_containing(
      [&] { Buffer allocation_probe(context, sizeof(double)); },
      "allocation");

  const VectorLayout zero_layout;
  FakeOperator zero_operator(context, zero_layout, zero_layout, 1U, {});
  JacobiPreconditioner zero(context);
  zero.update(zero_operator, 1U);
  HUNDUN_CHECK(zero_operator.diagonal_calls() == 1U);
  Buffer empty(context, 0U);
  zero.apply(empty.view(0U, 0U), empty.view(0U, 0U)).wait();
}

void test_jacobi_update_failures_are_transactional() {
  CpuReferenceContext context;
  const auto layout = make_layout(2U);
  FakeOperator good(context, layout, layout, 3U, {2.0, 4.0});
  JacobiPreconditioner jacobi(context);
  jacobi.update(good, 3U);

  Buffer residual(context, 2U * sizeof(double));
  Buffer correction(context, 2U * sizeof(double));
  set_values(residual.view(0U, 2U), {6.0, 8.0});
  auto prove_old_cache = [&] {
    set_values(correction.view(0U, 2U), {-9.0, -9.0});
    jacobi.apply(residual.view(0U, 2U), correction.view(0U, 2U)).wait();
    HUNDUN_CHECK(read_values(correction.view(0U, 2U)) ==
                 std::vector<double>({3.0, 2.0}));
  };
  prove_old_cache();

  auto failed_refresh = [&](std::vector<double> diagonal,
                            const std::string& text) {
    FakeOperator candidate(context, layout, layout, 4U,
                           std::move(diagonal));
    expect_error_containing([&] { jacobi.update(candidate, 4U); }, text);
    HUNDUN_CHECK(candidate.diagonal_calls() == 1U);
    prove_old_cache();
  };
  failed_refresh({0.0, 1.0}, "zero");
  failed_refresh({-0.0, 1.0}, "zero");
  failed_refresh({std::numeric_limits<double>::quiet_NaN(), 1.0}, "finite");
  failed_refresh({std::numeric_limits<double>::infinity(), 1.0}, "finite");
  failed_refresh({std::numeric_limits<double>::denorm_min(), 1.0},
                 "reciprocal");

  FakeOperator no_capability(context, layout, layout, 4U, {1.0, 1.0});
  no_capability.set_has_diagonal(false);
  expect_error_containing([&] { jacobi.update(no_capability, 4U); },
                          "diagonal");
  HUNDUN_CHECK(no_capability.diagonal_calls() == 0U);
  prove_old_cache();

  FakeOperator event_failure(context, layout, layout, 4U, {1.0, 1.0});
  event_failure.set_fail_diagonal_event(true);
  expect_error_containing([&] { jacobi.update(event_failure, 4U); },
                          "event error");
  HUNDUN_CHECK(event_failure.diagonal_calls() == 1U);
  prove_old_cache();

  FakeOperator throwing(context, layout, layout, 4U, {1.0, 1.0});
  throwing.set_throw_on_diagonal(true);
  expect_error_containing([&] { jacobi.update(throwing, 4U); }, "throw");
  HUNDUN_CHECK(throwing.diagonal_calls() == 1U);
  prove_old_cache();

  FakeOperator allocation_failure(context, layout, layout, 4U, {1.0, 1.0});
  ExecutionTestAccess::fail_next_allocation();
  expect_error_containing([&] { jacobi.update(allocation_failure, 4U); },
                          "allocation");
  HUNDUN_CHECK(allocation_failure.diagonal_calls() == 0U);
  prove_old_cache();

  CpuReferenceContext other_context;
  FakeOperator wrong_context(other_context, layout, layout, 4U, {1.0, 1.0});
  expect_error_containing([&] { jacobi.update(wrong_context, 4U); },
                          "context");
  prove_old_cache();
}

void test_jacobi_apply_failures_do_not_mutate() {
  CpuReferenceContext context;
  const auto layout = make_layout(2U);
  FakeOperator linear_operator(context, layout, layout, 1U,
                               {std::numeric_limits<double>::min(), -2.0});
  JacobiPreconditioner jacobi(context);
  jacobi.update(linear_operator, 1U);

  Buffer residual(context, 4U * sizeof(double));
  Buffer correction(context, 4U * sizeof(double));
  auto r = residual.view(0U, 2U);
  auto z = correction.view(0U, 2U);
  set_values(r, {8.0, 4.0});
  set_values(z, {5.0, 6.0});

  auto require_unchanged = [&](auto&& operation, const std::string& text) {
    const auto before = read_values(z);
    expect_error_containing(std::forward<decltype(operation)>(operation),
                            text);
    HUNDUN_CHECK(read_values(z) == before);
  };
  require_unchanged([&] { jacobi.apply(r, z); }, "finite");

  r[0] = 1.0;
  r[1] = std::numeric_limits<double>::infinity();
  require_unchanged([&] { jacobi.apply(r, z); }, "finite");
  r[1] = 2.0;

  require_unchanged(
      [&] { jacobi.apply(residual.view(0U, 1U), z); }, "size");

  Buffer shared(context, 3U * sizeof(double));
  set_values(shared.view(0U, 3U), {1.0, 2.0, 3.0});
  const auto shared_before = read_values(shared.view(0U, 3U));
  expect_error_containing(
      [&] { jacobi.apply(shared.view(0U, 2U), shared.view(1U, 2U)); },
      "alias");
  HUNDUN_CHECK(read_values(shared.view(0U, 3U)) == shared_before);

  auto metadata = ExecutionTestAccess::metadata(r);
  ++metadata.backend_identity;
  auto foreign = ExecutionTestAccess::const_view(residual, metadata);
  require_unchanged([&] { jacobi.apply(foreign, z); }, "backend");

  auto stale = residual.view(0U, 2U);
  residual.reallocate(4U * sizeof(double));
  require_unchanged([&] { jacobi.apply(stale, z); }, "live");
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_completed_event_factory();
    test_completed_event_allocation_failure_is_one_shot();
    test_virtual_destruction();
    test_context_and_update_validation();
    test_current_operator_revision_lineages();
    test_reused_operator_address_starts_fresh_lineage();
    test_completed_event_failure_preserves_apply_outputs();
    test_identity_apply_contract();
    test_jacobi_apply_and_cache();
    test_jacobi_update_failures_are_transactional();
    test_jacobi_apply_failures_do_not_mutate();
  });
}
