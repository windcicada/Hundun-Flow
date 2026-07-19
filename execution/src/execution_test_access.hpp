// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/execution/execution.hpp"

#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <string>

namespace hundun::execution::test {

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
  static AllocationIdentity next_allocation_identity() noexcept;
  static void set_next_allocation_identity(AllocationIdentity next) noexcept;
  static void fail_next_allocation() noexcept;
  static void set_epoch(Buffer& buffer, std::uint64_t epoch);

  static TestViewMetadata metadata(const VectorView<double>& view) noexcept;
  static VectorView<double> mutable_view(Buffer& buffer,
                                         TestViewMetadata metadata);
  static VectorView<const double> const_view(const Buffer& buffer,
                                             TestViewMetadata metadata);

  static std::weak_ptr<void> observe_allocation(const Buffer& buffer);
  static ExecutionEvent pending_event(
      std::initializer_list<Buffer*> retained_buffers);
  static void complete_success(ExecutionEvent& event);
  static void complete_error(ExecutionEvent& event, std::string message);
};

}  // namespace hundun::execution::test
