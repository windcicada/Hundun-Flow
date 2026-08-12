// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/exec_execution.hpp"

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace hundun::execution::detail {

AllocationIdentity execution_next_allocation_identity_raw() noexcept;
void execution_set_allocation_counters_raw(AllocationCounters) noexcept;
void execution_set_next_allocation_identity_raw(AllocationIdentity) noexcept;
void execution_fail_next_allocation_raw() noexcept;
void execution_fail_next_completed_event_allocation_raw() noexcept;
void execution_set_control_epoch_raw(
    const std::shared_ptr<AllocationControl> &, std::uint64_t epoch);

std::shared_ptr<AllocationControl> execution_make_metadata_control_raw(
    std::size_t byte_size, AllocationIdentity allocation_identity,
    std::uint64_t epoch, BackendIdentity backend_identity,
    ExecutionSpace space);
bool execution_control_has_storage_raw(
    const std::shared_ptr<AllocationControl> &) noexcept;
bool execution_control_owner_live_raw(
    const std::shared_ptr<AllocationControl> &) noexcept;

std::shared_ptr<EventState> execution_make_pending_event_raw(
    std::vector<std::shared_ptr<AllocationControl>> retained_allocations);
void execution_complete_success_raw(const std::shared_ptr<EventState> &);
void execution_complete_error_raw(const std::shared_ptr<EventState> &,
                                  std::string message);
void execution_wait_until_waiter_count_raw(
    const std::shared_ptr<EventState> &, std::size_t expected);

} // namespace hundun::execution::detail

namespace hundun::execution::test {

class MetadataOnlyViewFixture final {
public:
  ~MetadataOnlyViewFixture() = default;
  MetadataOnlyViewFixture(MetadataOnlyViewFixture &&) noexcept = default;
  MetadataOnlyViewFixture &
  operator=(MetadataOnlyViewFixture &&) noexcept = default;
  MetadataOnlyViewFixture(const MetadataOnlyViewFixture &) = delete;
  MetadataOnlyViewFixture &
  operator=(const MetadataOnlyViewFixture &) = delete;

private:
  friend class ExecutionTestAccess;

  MetadataOnlyViewFixture(
      std::shared_ptr<detail::AllocationControl> control,
      AllocationIdentity allocation_identity,
      std::uint64_t epoch, BackendIdentity backend_identity,
      ExecutionSpace space) noexcept
      : control_(std::move(control)), allocation_identity_(allocation_identity), epoch_(epoch),
        backend_identity_(backend_identity), space_(space) {}

  std::shared_ptr<detail::AllocationControl> control_;
  AllocationIdentity allocation_identity_{};
  std::uint64_t epoch_{};
  BackendIdentity backend_identity_{};
  ExecutionSpace space_{ExecutionSpace::host};
};

struct TestViewMetadata {
  std::size_t offset_bytes;
  std::size_t element_count;
  std::size_t stride_elements;
  AllocationIdentity allocation_identity;
  std::uint64_t epoch;
  BackendIdentity backend_identity;
  ExecutionSpace space;
  bool writable;
};

class ExecutionTestAccess final {
public:
  static void set_allocation_counters_for_test(
      AllocationCounters counters) noexcept {
    detail::execution_set_allocation_counters_raw(counters);
  }

  static AllocationIdentity next_allocation_identity() noexcept {
    return detail::execution_next_allocation_identity_raw();
  }

  static void set_next_allocation_identity(AllocationIdentity next) noexcept {
    detail::execution_set_next_allocation_identity_raw(next);
  }

  static void fail_next_allocation() noexcept {
    detail::execution_fail_next_allocation_raw();
  }

  static void fail_next_completed_event_allocation() noexcept {
    detail::execution_fail_next_completed_event_allocation_raw();
  }

  static void set_epoch(Buffer &buffer, std::uint64_t epoch) {
    buffer.require_live();
    detail::execution_set_control_epoch_raw(buffer.control_, epoch);
  }

  static TestViewMetadata metadata(const VectorView<double> &view) noexcept {
    return TestViewMetadata{view.offset_bytes_,
                            view.element_count_,
                            view.stride_elements_,
                            view.allocation_identity_,
                            view.epoch_,
                            view.backend_identity_,
                            view.space_,
                            view.writable_};
  }

  static VectorView<double> mutable_view(Buffer &buffer,
                                         TestViewMetadata metadata) {
    buffer.require_live();
    return VectorView<double>(
        buffer.control_, metadata.offset_bytes, metadata.element_count,
        metadata.stride_elements, metadata.allocation_identity,
        metadata.epoch, metadata.backend_identity, metadata.space,
        metadata.writable);
  }

  static VectorView<const double> const_view(const Buffer &buffer,
                                             TestViewMetadata metadata) {
    buffer.require_live();
    return VectorView<const double>(
        buffer.control_, metadata.offset_bytes, metadata.element_count,
        metadata.stride_elements, metadata.allocation_identity,
        metadata.epoch, metadata.backend_identity, metadata.space, false);
  }

  static MetadataOnlyViewFixture metadata_only_allocation(
      std::size_t byte_size, AllocationIdentity allocation_identity,
      std::uint64_t epoch, BackendIdentity backend_identity,
      ExecutionSpace space) {
    auto control = detail::execution_make_metadata_control_raw(
        byte_size, allocation_identity, epoch, backend_identity, space);
    return MetadataOnlyViewFixture(std::move(control), allocation_identity, epoch,
                                   backend_identity, space);
  }

  static bool has_storage(
      const MetadataOnlyViewFixture &fixture) noexcept {
    return detail::execution_control_has_storage_raw(fixture.control_);
  }

  static VectorView<double> mutable_view(
      const MetadataOnlyViewFixture &fixture, std::size_t offset_bytes,
      std::size_t element_count, std::size_t stride_elements,
      bool writable) {
    if (!detail::execution_control_owner_live_raw(fixture.control_))
      throw runtime::Error("metadata-only allocation fixture is not live");
    return VectorView<double>(
        fixture.control_, offset_bytes, element_count, stride_elements,
        fixture.allocation_identity_, fixture.epoch_,
        fixture.backend_identity_, fixture.space_, writable);
  }

  static VectorView<const double> const_view(
      const MetadataOnlyViewFixture &fixture, std::size_t offset_bytes,
      std::size_t element_count, std::size_t stride_elements) {
    if (!detail::execution_control_owner_live_raw(fixture.control_))
      throw runtime::Error("metadata-only allocation fixture is not live");
    return VectorView<const double>(
        fixture.control_, offset_bytes, element_count, stride_elements,
        fixture.allocation_identity_, fixture.epoch_,
        fixture.backend_identity_, fixture.space_, false);
  }

  static std::weak_ptr<void> observe_allocation(const Buffer &buffer) {
    buffer.require_live();
    return buffer.control_;
  }

  static ExecutionEvent pending_event(
      std::initializer_list<Buffer *> retained_buffers) {
    std::vector<std::shared_ptr<detail::AllocationControl>> retained;
    retained.reserve(retained_buffers.size());
    for (Buffer *buffer : retained_buffers) {
      if (buffer == nullptr)
        throw runtime::Error("pending event cannot retain a null buffer");
      buffer->require_live();
      retained.push_back(buffer->control_);
    }
    return ExecutionEvent(
        detail::execution_make_pending_event_raw(std::move(retained)));
  }

  static void complete_success(ExecutionEvent &event) {
    detail::execution_complete_success_raw(event.state_);
  }

  static void complete_error(ExecutionEvent &event, std::string message) {
    detail::execution_complete_error_raw(event.state_, std::move(message));
  }

  static void wait_until_waiter_count(const ExecutionEvent &event,
                                      std::size_t expected) {
    detail::execution_wait_until_waiter_count_raw(event.state_, expected);
  }
};

} // namespace hundun::execution::test
