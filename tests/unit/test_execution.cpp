// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"

#include "execution/src/execution_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::execution::AllocationIdentity;
using hundun::execution::BackendIdentity;
using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::execution::ExecutionCapability;
using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionEvent;
using hundun::execution::ExecutionSpace;
using hundun::execution::ScalarFormat;
using hundun::execution::VectorView;
using hundun::execution::test::ExecutionTestAccess;
using hundun::execution::test::MetadataOnlyViewFixture;
using hundun::execution::test::TestViewMetadata;
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
  if (message.find(text) == std::string::npos) {
    throw std::runtime_error("expected error containing '" + text +
                             "', got '" + message + "'");
  }
}

void join_and_capture(std::thread& thread,
                      std::exception_ptr& cleanup_error) noexcept {
  if (!thread.joinable()) {
    return;
  }
  try {
    thread.join();
  } catch (...) {
    if (cleanup_error == nullptr) {
      cleanup_error = std::current_exception();
    }
  }
}

class ConfigurableContext final : public ExecutionContext {
 public:
  ConfigurableContext(BackendIdentity identity, ExecutionSpace space,
                      bool allocation, bool host_access, bool transfer,
                      bool asynchronous) noexcept
      : identity_(identity),
        space_(space),
        allocation_(allocation),
        host_access_(host_access),
        transfer_(transfer),
        asynchronous_(asynchronous) {}

  std::string_view backend_name() const noexcept override {
    return space_ == ExecutionSpace::device ? "device_test_double"
                                             : "host_test_double";
  }
  BackendIdentity backend_identity() const noexcept override {
    return identity_;
  }
  ExecutionSpace space() const noexcept override { return space_; }
  bool ordered() const noexcept override { return true; }
  bool supports(ExecutionCapability capability) const noexcept override {
    switch (capability) {
      case ExecutionCapability::buffer_allocation:
        return allocation_;
      case ExecutionCapability::host_access:
        return host_access_;
      case ExecutionCapability::transfer:
        return transfer_;
      case ExecutionCapability::asynchronous_event:
        return asynchronous_;
    }
    return false;
  }

 private:
  BackendIdentity identity_;
  ExecutionSpace space_;
  bool allocation_;
  bool host_access_;
  bool transfer_;
  bool asynchronous_;
};

void check_values(VectorView<const double> view,
                  const std::vector<double>& expected) {
  HUNDUN_CHECK(view.size() == expected.size());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    HUNDUN_CHECK_NEAR(view[index], expected[index], 0.0);
  }
}

std::vector<double> values(VectorView<const double> view) {
  std::vector<double> result;
  result.reserve(view.size());
  for (std::size_t index = 0; index < view.size(); ++index) {
    result.push_back(view[index]);
  }
  return result;
}

void test_public_traits_and_context() {
  static_assert(std::has_virtual_destructor_v<ExecutionContext>);
  static_assert(std::is_nothrow_destructible_v<ExecutionContext>);
  static_assert(std::is_constructible_v<VectorView<const double>,
                                        VectorView<double>>);
  static_assert(!std::is_constructible_v<VectorView<double>,
                                         VectorView<const double>>);
  static_assert(!std::is_copy_constructible_v<Buffer>);
  static_assert(std::is_nothrow_move_constructible_v<Buffer>);
  static_assert(std::is_nothrow_move_assignable_v<Buffer>);
  static_assert(std::is_copy_constructible_v<ExecutionEvent>);
  static_assert(std::is_nothrow_move_constructible_v<ExecutionEvent>);

  CpuReferenceContext first;
  CpuReferenceContext second;
  HUNDUN_CHECK(first.backend_name() == "cpu_reference");
  HUNDUN_CHECK(first.backend_identity() != 0);
  HUNDUN_CHECK(first.backend_identity() == second.backend_identity());
  HUNDUN_CHECK(first.space() == ExecutionSpace::host);
  HUNDUN_CHECK(first.ordered());
  HUNDUN_CHECK(first.supports(ExecutionCapability::buffer_allocation));
  HUNDUN_CHECK(first.supports(ExecutionCapability::host_access));
  HUNDUN_CHECK(first.supports(ExecutionCapability::transfer));
  HUNDUN_CHECK(!first.supports(ExecutionCapability::asynchronous_event));
  HUNDUN_CHECK(!first.supports(static_cast<ExecutionCapability>(999)));

  struct DestructionContext final : ExecutionContext {
    explicit DestructionContext(bool& destroyed) : destroyed_(destroyed) {}
    ~DestructionContext() noexcept override { destroyed_ = true; }
    std::string_view backend_name() const noexcept override { return "probe"; }
    BackendIdentity backend_identity() const noexcept override { return 7; }
    ExecutionSpace space() const noexcept override { return ExecutionSpace::host; }
    bool ordered() const noexcept override { return true; }
    bool supports(ExecutionCapability) const noexcept override { return false; }
    bool& destroyed_;
  };
  bool destroyed = false;
  std::unique_ptr<ExecutionContext> base =
      std::make_unique<DestructionContext>(destroyed);
  base.reset();
  HUNDUN_CHECK(destroyed);
}

void test_buffers_and_views() {
  CpuReferenceContext context;

  Buffer zero(context, 0);
  HUNDUN_CHECK(zero.byte_size() == 0);
  HUNDUN_CHECK(zero.allocation_identity() != 0);
  HUNDUN_CHECK(zero.epoch() != 0);
  HUNDUN_CHECK(zero.backend_identity() == context.backend_identity());
  HUNDUN_CHECK(zero.space() == ExecutionSpace::host);
  auto empty = zero.view(0, 0);
  HUNDUN_CHECK(empty.size() == 0);
  HUNDUN_CHECK(empty.stride() == 1);
  HUNDUN_CHECK(empty.data() == nullptr);
  HUNDUN_CHECK(empty.writable());
  HUNDUN_CHECK(empty.scalar_format() == ScalarFormat::float64);
  HUNDUN_CHECK(empty.offset_bytes() == 0);
  expect_error_containing([&] { static_cast<void>(zero.view(0, 0, 0)); },
                          "stride");

  Buffer first(context, 6 * sizeof(double));
  Buffer second(context, sizeof(double));
  HUNDUN_CHECK(first.allocation_identity() != second.allocation_identity());
  auto all = first.view(0, 6);
  HUNDUN_CHECK(reinterpret_cast<std::uintptr_t>(all.data()) % alignof(double) ==
               0);
  for (std::size_t index = 0; index < all.size(); ++index) {
    all[index] = static_cast<double>(index + 1);
  }
  const Buffer& constant_first = first;
  VectorView<const double> converted = all;
  HUNDUN_CHECK(converted.size() == all.size());
  HUNDUN_CHECK(converted.stride() == all.stride());
  HUNDUN_CHECK(converted.offset_bytes() == all.offset_bytes());
  HUNDUN_CHECK(converted.allocation_identity() == all.allocation_identity());
  HUNDUN_CHECK(converted.epoch() == all.epoch());
  HUNDUN_CHECK(!converted.writable());
  auto odds = constant_first.view(1, 3, 2);
  HUNDUN_CHECK(!odds.writable());
  HUNDUN_CHECK(odds.offset_bytes() == sizeof(double));
  HUNDUN_CHECK(odds.epoch() == first.epoch());
  HUNDUN_CHECK(odds.allocation_identity() == first.allocation_identity());
  HUNDUN_CHECK(odds.backend_identity() == context.backend_identity());
  HUNDUN_CHECK(odds.space() == ExecutionSpace::host);
  HUNDUN_CHECK(odds.data() != nullptr);
  check_values(odds, {2.0, 4.0, 6.0});

  auto at_end = first.view(6, 0, 3);
  HUNDUN_CHECK(at_end.data() == nullptr);
  expect_error_containing([&] { static_cast<void>(first.view(7, 0)); },
                          "range");
  expect_error_containing([&] { static_cast<void>(first.view(0, 7)); },
                          "range");
  expect_error_containing([&] { static_cast<void>(first.view(5, 2)); },
                          "range");
  expect_error_containing(
      [&] {
        static_cast<void>(first.view(
            std::numeric_limits<std::size_t>::max() / sizeof(double) + 1, 0));
      },
      "overflow");
  expect_error_containing(
      [&] {
        static_cast<void>(first.view(0, 2,
                                     std::numeric_limits<std::size_t>::max()));
      },
      "overflow");
  expect_error_containing([&] { static_cast<void>(all[all.size()]); },
                          "bounds");

  auto before_move = first.view(0, 2);
  const auto first_identity = first.allocation_identity();
  const auto first_epoch = first.epoch();
  Buffer moved(std::move(first));
  HUNDUN_CHECK(moved.allocation_identity() == first_identity);
  HUNDUN_CHECK(moved.epoch() == first_epoch);
  check_values(before_move, {1.0, 2.0});
  HUNDUN_CHECK(first.byte_size() == 0);
  HUNDUN_CHECK(first.allocation_identity() == 0);
  expect_error_containing([&] { static_cast<void>(first.view(0, 0)); },
                          "moved-from");

  Buffer destination(context, 2 * sizeof(double));
  auto destination_old = destination.view(0, 2);
  destination_old[0] = -1.0;
  auto source_view = moved.view(0, 2);
  destination = std::move(moved);
  expect_error_containing([&] { static_cast<void>(destination_old[0]); },
                          "live");
  check_values(source_view, {1.0, 2.0});
  HUNDUN_CHECK(destination.allocation_identity() == first_identity);
  expect_error_containing([&] { static_cast<void>(moved.view(0, 0)); },
                          "moved-from");

  auto stale = destination.view(0, 2);
  const auto old_identity = destination.allocation_identity();
  const auto old_epoch = destination.epoch();
  destination.reallocate(2 * sizeof(double));
  HUNDUN_CHECK(destination.allocation_identity() != old_identity);
  HUNDUN_CHECK(destination.epoch() == old_epoch + 1);
  HUNDUN_CHECK(destination.byte_size() == 2 * sizeof(double));
  expect_error_containing([&] { static_cast<void>(stale[0]); }, "live");

  auto current = destination.view(0, 2);
  current[0] = 31.0;
  current[1] = 41.0;
  const auto stable_identity = destination.allocation_identity();
  const auto stable_epoch = destination.epoch();
  ExecutionTestAccess::fail_next_allocation();
  expect_error_containing(
      [&] { destination.reallocate(4 * sizeof(double)); }, "allocation");
  HUNDUN_CHECK(destination.allocation_identity() == stable_identity);
  HUNDUN_CHECK(destination.epoch() == stable_epoch);
  check_values(destination.view(0, 2), {31.0, 41.0});
  expect_error_containing(
      [&] {
        destination.reallocate(std::numeric_limits<std::size_t>::max());
      },
      "byte size");
  HUNDUN_CHECK(destination.allocation_identity() == stable_identity);
  HUNDUN_CHECK(destination.epoch() == stable_epoch);
  check_values(destination.view(0, 2), {31.0, 41.0});

  ExecutionTestAccess::set_epoch(
      destination, std::numeric_limits<std::uint64_t>::max());
  expect_error_containing([&] { destination.reallocate(sizeof(double)); },
                          "epoch");
  HUNDUN_CHECK(destination.allocation_identity() == stable_identity);
  check_values(destination.view(0, 2), {31.0, 41.0});

  VectorView<double> dead = [&] {
    Buffer temporary(context, sizeof(double));
    auto view = temporary.view(0, 1);
    view[0] = 9.0;
    return view;
  }();
  expect_error_containing([&] { static_cast<void>(dead.data()); }, "live");

  const auto saved_next = ExecutionTestAccess::next_allocation_identity();
  ExecutionTestAccess::set_next_allocation_identity(
      std::numeric_limits<AllocationIdentity>::max());
  expect_error_containing(
      [&] {
        Buffer wrapped(context, 0);
        static_cast<void>(wrapped);
      },
      "identity");
  ExecutionTestAccess::set_next_allocation_identity(saved_next);

  ConfigurableContext zero_identity(0, ExecutionSpace::host, true, true, true,
                                    false);
  expect_error_containing(
      [&] {
        Buffer invalid(zero_identity, 0);
        static_cast<void>(invalid);
      },
      "identity");
  ConfigurableContext device(88, ExecutionSpace::device, true, false, true,
                             true);
  expect_error_containing(
      [&] {
        Buffer invalid(device, sizeof(double));
        static_cast<void>(invalid);
      },
      "host");
  HUNDUN_CHECK(device.backend_name() == "device_test_double");
  HUNDUN_CHECK(device.space() == ExecutionSpace::device);
  HUNDUN_CHECK(device.supports(ExecutionCapability::buffer_allocation));
  HUNDUN_CHECK(!device.supports(ExecutionCapability::host_access));
  HUNDUN_CHECK(device.supports(ExecutionCapability::transfer));
  HUNDUN_CHECK(device.supports(ExecutionCapability::asynchronous_event));
  HUNDUN_CHECK(!device.supports(static_cast<ExecutionCapability>(999)));
  ConfigurableContext no_allocation(89, ExecutionSpace::host, false, true,
                                    true, false);
  expect_error_containing(
      [&] {
        Buffer invalid(no_allocation, sizeof(double));
        static_cast<void>(invalid);
      },
      "allocation");
}

void test_transfers_and_validation() {
  CpuReferenceContext context;
  Buffer source(context, 8 * sizeof(double));
  Buffer destination(context, 8 * sizeof(double));
  auto source_values = source.view(0, 8);
  auto destination_values = destination.view(0, 8);
  for (std::size_t index = 0; index < 8; ++index) {
    source_values[index] = 10.0 + static_cast<double>(index);
    destination_values[index] = -5.0 - static_cast<double>(index);
  }
  const auto source_identity = source.allocation_identity();
  const auto source_epoch = source.epoch();
  const auto destination_identity = destination.allocation_identity();
  const auto destination_epoch = destination.epoch();

  auto contiguous_event = hundun::execution::transfer(
      VectorView<const double>(source.view(0, 4)), destination.view(0, 4),
      context);
  HUNDUN_CHECK(contiguous_event.ready());
  contiguous_event.wait();
  check_values(destination.view(0, 4), {10.0, 11.0, 12.0, 13.0});

  auto strided_event = hundun::execution::transfer(
      VectorView<const double>(source.view(1, 3, 2)),
      destination.view(1, 3, 2), context);
  strided_event.wait();
  check_values(destination.view(1, 3, 2), {11.0, 13.0, 15.0});

  auto zero_event = hundun::execution::transfer(
      VectorView<const double>(source.view(8, 0)), destination.view(8, 0),
      context);
  HUNDUN_CHECK(zero_event.ready());
  zero_event.wait();

  const auto source_before_self = values(source.view(0, 8));
  hundun::execution::transfer(VectorView<const double>(source.view(0, 4, 2)),
                              source.view(0, 4, 2), context)
      .wait();
  HUNDUN_CHECK(values(source.view(0, 8)) == source_before_self);

  Buffer same(context, 6 * sizeof(double));
  auto same_all = same.view(0, 6);
  for (std::size_t index = 0; index < 6; ++index) {
    same_all[index] = static_cast<double>(index + 20);
  }
  hundun::execution::transfer(VectorView<const double>(same.view(0, 2)),
                              same.view(4, 2), context)
      .wait();
  check_values(same.view(0, 6), {20.0, 21.0, 22.0, 23.0, 20.0, 21.0});

  HUNDUN_CHECK(source.allocation_identity() == source_identity);
  HUNDUN_CHECK(source.epoch() == source_epoch);
  HUNDUN_CHECK(destination.allocation_identity() == destination_identity);
  HUNDUN_CHECK(destination.epoch() == destination_epoch);

  auto require_unchanged_failure = [&](auto&& operation,
                                       const std::string& message) {
    const auto before = values(destination.view(0, 8));
    expect_error_containing(std::forward<decltype(operation)>(operation),
                            message);
    HUNDUN_CHECK(values(destination.view(0, 8)) == before);
  };

  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)),
            destination.view(0, 3), context));
      },
      "size");

  Buffer stale_owner(context, 2 * sizeof(double));
  auto stale_source = VectorView<const double>(stale_owner.view(0, 2));
  stale_owner.reallocate(2 * sizeof(double));
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            stale_source, destination.view(0, 2), context));
      },
      "live");

  VectorView<const double> dead_source = [&] {
    Buffer temporary(context, 2 * sizeof(double));
    return VectorView<const double>(temporary.view(0, 2));
  }();
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            dead_source, destination.view(0, 2), context));
      },
      "live");

  const auto same_before_overlap = values(same.view(0, 6));
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(same.view(0, 3)), same.view(1, 3),
            context));
      },
      "overlap");
  HUNDUN_CHECK(values(same.view(0, 6)) == same_before_overlap);
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(same.view(0, 2, 2)), same.view(1, 2, 2),
            context));
      },
      "overlap");

  ConfigurableContext wrong_identity(context.backend_identity() + 1,
                                     ExecutionSpace::host, true, true, true,
                                     false);
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)),
            destination.view(0, 2), wrong_identity));
      },
      "identity");
  ConfigurableContext zero_identity(0, ExecutionSpace::host, true, true, true,
                                    false);
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)),
            destination.view(0, 2), zero_identity));
      },
      "identity");
  ConfigurableContext compatible_by_identity(
      context.backend_identity(), ExecutionSpace::host, true, true, true,
      false);
  hundun::execution::transfer(
      VectorView<const double>(source.view(0, 2)), destination.view(0, 2),
      compatible_by_identity)
      .wait();
  check_values(destination.view(0, 2), {10.0, 11.0});
  ConfigurableContext no_transfer(context.backend_identity(),
                                  ExecutionSpace::host, true, true, false,
                                  false);
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)),
            destination.view(0, 2), no_transfer));
      },
      "transfer");
  ConfigurableContext no_host(context.backend_identity(), ExecutionSpace::host,
                              true, false, true, false);
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)),
            destination.view(0, 2), no_host));
      },
      "host");
  ConfigurableContext device_context(context.backend_identity(),
                                     ExecutionSpace::device, true, false, true,
                                     true);
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)),
            destination.view(0, 2), device_context));
      },
      "host");

  const TestViewMetadata valid_destination =
      ExecutionTestAccess::metadata(destination.view(0, 2));
  auto non_writable = ExecutionTestAccess::mutable_view(
      destination, TestViewMetadata{valid_destination.offset_bytes,
                                    valid_destination.element_count,
                                    valid_destination.stride_elements,
                                    valid_destination.allocation_identity,
                                    valid_destination.epoch,
                                    valid_destination.backend_identity,
                                    valid_destination.space,
                                    false});
  expect_error_containing(
      [&] { static_cast<void>(non_writable.data()); }, "writable");
  expect_error_containing(
      [&] { static_cast<void>(non_writable[0]); }, "writable");
  HUNDUN_CHECK(destination.view(0, 2)[0] == 10.0);
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            VectorView<const double>(source.view(0, 2)), non_writable,
            context));
      },
      "writable");

  constexpr BackendIdentity device_identity = 999;
  constexpr AllocationIdentity device_allocation = 7001;
  constexpr std::uint64_t device_epoch = 9;
  MetadataOnlyViewFixture device_fixture =
      ExecutionTestAccess::metadata_only_allocation(
          2 * sizeof(double), device_allocation, device_epoch,
          device_identity, ExecutionSpace::device);
  HUNDUN_CHECK(!ExecutionTestAccess::has_storage(device_fixture));
  auto device_source = ExecutionTestAccess::const_view(
      device_fixture, 0, 2, 1);
  auto device_destination = ExecutionTestAccess::mutable_view(
      device_fixture, 0, 2, 1, true);
  ConfigurableContext matching_device_context(
      device_identity, ExecutionSpace::device, true, false, true, true);
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            device_source, device_destination, matching_device_context));
      },
      "host context");

  constexpr AllocationIdentity cpu_device_allocation = 7002;
  MetadataOnlyViewFixture cpu_device_fixture =
      ExecutionTestAccess::metadata_only_allocation(
          2 * sizeof(double), cpu_device_allocation, device_epoch,
          context.backend_identity(), ExecutionSpace::device);
  HUNDUN_CHECK(!ExecutionTestAccess::has_storage(cpu_device_fixture));
  auto cpu_identity_device_source = ExecutionTestAccess::const_view(
      cpu_device_fixture, 0, 2, 1);
  expect_error_containing(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            cpu_identity_device_source, destination.view(0, 2), context));
      },
      "host vector views");

  auto misaligned = ExecutionTestAccess::const_view(
      source, TestViewMetadata{1,
                               1,
                               1,
                               source.allocation_identity(),
                               source.epoch(),
                               source.backend_identity(),
                               source.space(),
                               false});
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            misaligned, destination.view(0, 1), context));
      },
      "alignment");

  auto zero_stride = ExecutionTestAccess::const_view(
      source, TestViewMetadata{0,
                               1,
                               0,
                               source.allocation_identity(),
                               source.epoch(),
                               source.backend_identity(),
                               source.space(),
                               false});
  require_unchanged_failure(
      [&] {
        static_cast<void>(hundun::execution::transfer(
            zero_stride, destination.view(0, 1), context));
      },
      "stride");
}

void test_events_and_lifetime() {
  CpuReferenceContext context;
  Buffer source(context, sizeof(double));
  Buffer destination(context, sizeof(double));
  source.view(0, 1)[0] = 3.0;
  auto event = hundun::execution::transfer(
      VectorView<const double>(source.view(0, 1)), destination.view(0, 1),
      context);
  HUNDUN_CHECK(event.ready());
  event.wait();
  event.wait();
  ExecutionEvent copied = event;
  HUNDUN_CHECK(copied.ready());
  copied.wait();
  ExecutionEvent moved = std::move(event);
  HUNDUN_CHECK(moved.ready());
  HUNDUN_CHECK(!event.ready());
  expect_error_containing([&] { event.wait(); }, "moved-from");
  moved.wait();
  auto move_assigned = hundun::execution::transfer(
      VectorView<const double>(source.view(1, 0)), destination.view(1, 0),
      context);
  move_assigned = std::move(copied);
  HUNDUN_CHECK(move_assigned.ready());
  HUNDUN_CHECK(!copied.ready());
  expect_error_containing([&] { copied.wait(); }, "moved-from");

  auto retained_buffer =
      std::make_unique<Buffer>(context, 2 * sizeof(double));
  auto retained_view = retained_buffer->view(0, 2);
  retained_view[0] = 7.0;
  const auto observed = ExecutionTestAccess::observe_allocation(*retained_buffer);
  auto pending = ExecutionTestAccess::pending_event({retained_buffer.get()});
  ExecutionEvent pending_copy = pending;
  ExecutionEvent pending_waiter_copy_one = pending;
  ExecutionEvent pending_waiter_copy_two = pending;
  HUNDUN_CHECK(!pending.ready());
  HUNDUN_CHECK(!pending_copy.ready());
  HUNDUN_CHECK(!pending_waiter_copy_one.ready());
  HUNDUN_CHECK(!pending_waiter_copy_two.ready());
  std::atomic<bool> wait_one_returned{false};
  std::atomic<bool> wait_two_returned{false};
  std::exception_ptr wait_one_error;
  std::exception_ptr wait_two_error;
  std::thread waiter_one;
  std::thread waiter_two;
  bool ready_before_completion = false;
  bool wait_one_returned_before_completion = false;
  bool wait_two_returned_before_completion = false;
  bool retained_before_completion = false;
  std::string stale_view_message;
  try {
    waiter_one = std::thread([&] {
      try {
        pending_waiter_copy_one.wait();
      } catch (...) {
        wait_one_error = std::current_exception();
      }
      wait_one_returned.store(true, std::memory_order_release);
    });
    waiter_two = std::thread([&] {
      try {
        pending_waiter_copy_two.wait();
      } catch (...) {
        wait_two_error = std::current_exception();
      }
      wait_two_returned.store(true, std::memory_order_release);
    });
    ExecutionTestAccess::wait_until_waiter_count(pending_copy, 2);
    ready_before_completion = pending.ready();
    wait_one_returned_before_completion =
        wait_one_returned.load(std::memory_order_acquire);
    wait_two_returned_before_completion =
        wait_two_returned.load(std::memory_order_acquire);
    retained_buffer.reset();
    retained_before_completion = !observed.expired();
    try {
      static_cast<void>(retained_view[0]);
    } catch (const Error& error) {
      stale_view_message = error.what();
    }
    ExecutionTestAccess::complete_success(pending_copy);
    waiter_one.join();
    waiter_two.join();
  } catch (...) {
    const auto original_error = std::current_exception();
    std::exception_ptr cleanup_error;
    if (!pending_copy.ready()) {
      try {
        ExecutionTestAccess::complete_success(pending_copy);
      } catch (...) {
        cleanup_error = std::current_exception();
      }
    }
    join_and_capture(waiter_one, cleanup_error);
    join_and_capture(waiter_two, cleanup_error);
    if (cleanup_error != nullptr) {
      std::rethrow_exception(cleanup_error);
    }
    std::rethrow_exception(original_error);
  }
  HUNDUN_CHECK(!ready_before_completion);
  HUNDUN_CHECK(!wait_one_returned_before_completion);
  HUNDUN_CHECK(!wait_two_returned_before_completion);
  HUNDUN_CHECK(retained_before_completion);
  HUNDUN_CHECK(stale_view_message.find("live") != std::string::npos);
  HUNDUN_CHECK(wait_one_returned.load(std::memory_order_acquire));
  HUNDUN_CHECK(wait_two_returned.load(std::memory_order_acquire));
  HUNDUN_CHECK(wait_one_error == nullptr);
  HUNDUN_CHECK(wait_two_error == nullptr);
  HUNDUN_CHECK(pending.ready());
  HUNDUN_CHECK(pending_copy.ready());
  HUNDUN_CHECK(pending_waiter_copy_one.ready());
  HUNDUN_CHECK(pending_waiter_copy_two.ready());
  HUNDUN_CHECK(observed.expired());
  pending.wait();
  pending_copy.wait();
  pending_waiter_copy_one.wait();
  pending_waiter_copy_two.wait();

  auto failed = ExecutionTestAccess::pending_event({});
  HUNDUN_CHECK(!failed.ready());
  ExecutionEvent failed_copy = failed;
  std::atomic<bool> failed_wait_returned{false};
  std::string failed_wait_message;
  std::exception_ptr failed_wait_error;
  std::thread failed_waiter;
  bool failed_returned_before_completion = false;
  std::string normal_completion_message{"fixed numerical event error"};
  std::string cleanup_completion_message{"fixed numerical event error"};
  try {
    failed_waiter = std::thread([&] {
      try {
        failed_wait_message = expect_error([&] { failed_copy.wait(); });
      } catch (...) {
        failed_wait_error = std::current_exception();
      }
      failed_wait_returned.store(true, std::memory_order_release);
    });
    ExecutionTestAccess::wait_until_waiter_count(failed, 1);
    failed_returned_before_completion =
        failed_wait_returned.load(std::memory_order_acquire);
    ExecutionTestAccess::complete_error(
        failed, std::move(normal_completion_message));
    failed_waiter.join();
  } catch (...) {
    const auto original_error = std::current_exception();
    std::exception_ptr cleanup_error;
    if (!failed.ready()) {
      try {
        ExecutionTestAccess::complete_error(
            failed, std::move(cleanup_completion_message));
      } catch (...) {
        cleanup_error = std::current_exception();
      }
    }
    join_and_capture(failed_waiter, cleanup_error);
    if (cleanup_error != nullptr) {
      std::rethrow_exception(cleanup_error);
    }
    std::rethrow_exception(original_error);
  }
  HUNDUN_CHECK(!failed_returned_before_completion);
  HUNDUN_CHECK(failed_wait_returned.load(std::memory_order_acquire));
  if (failed_wait_error != nullptr) {
    std::rethrow_exception(failed_wait_error);
  }
  HUNDUN_CHECK(failed_wait_message == "fixed numerical event error");
  HUNDUN_CHECK(failed.ready());
  HUNDUN_CHECK(expect_error([&] { failed.wait(); }) ==
               "fixed numerical event error");
  HUNDUN_CHECK(expect_error([&] { failed.wait(); }) ==
               "fixed numerical event error");
  HUNDUN_CHECK(expect_error([&] { failed_copy.wait(); }) ==
               "fixed numerical event error");
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_public_traits_and_context();
    test_buffers_and_views();
    test_transfers_and_validation();
    test_events_and_lifetime();
  });
}
