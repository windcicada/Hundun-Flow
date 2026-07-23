// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdlib>
#include <new>

extern "C" void hundun_task22_record_allocation_attempt() noexcept
    __attribute__((weak));

namespace hundun::test::allocation_probe {

inline thread_local bool count_attempts = false;
inline thread_local std::size_t attempt_count = 0U;
inline thread_local bool fail_next_attempt = false;

inline void record_attempt() noexcept {
  if (hundun_task22_record_allocation_attempt) {
    hundun_task22_record_allocation_attempt();
  }
  if (count_attempts) {
    ++attempt_count;
  }
}

inline void* allocate(std::size_t bytes) {
  record_attempt();
  if (fail_next_attempt) {
    fail_next_attempt = false;
    throw std::bad_alloc();
  }
  if (void* pointer = std::malloc(bytes == 0U ? 1U : bytes)) {
    return pointer;
  }
  throw std::bad_alloc();
}

inline void* allocate_aligned(std::size_t bytes, std::size_t alignment) {
  record_attempt();
  if (fail_next_attempt) {
    fail_next_attempt = false;
    throw std::bad_alloc();
  }
  void* pointer = nullptr;
  if (posix_memalign(&pointer, alignment, bytes == 0U ? 1U : bytes) == 0) {
    return pointer;
  }
  throw std::bad_alloc();
}

class AllocationAttemptGuard final {
 public:
  AllocationAttemptGuard() noexcept {
    attempt_count = 0U;
    count_attempts = true;
  }
  ~AllocationAttemptGuard() noexcept { count_attempts = false; }
  AllocationAttemptGuard(const AllocationAttemptGuard&) = delete;
  AllocationAttemptGuard& operator=(const AllocationAttemptGuard&) = delete;

  std::size_t attempts() const noexcept { return attempt_count; }
};

class FailNextAllocationGuard final {
 public:
  explicit FailNextAllocationGuard(bool enabled = true) noexcept
      : enabled_(enabled) {
    fail_next_attempt = enabled;
  }
  ~FailNextAllocationGuard() noexcept {
    if (enabled_) {
      fail_next_attempt = false;
    }
  }
  FailNextAllocationGuard(const FailNextAllocationGuard&) = delete;
  FailNextAllocationGuard& operator=(const FailNextAllocationGuard&) = delete;

 private:
  bool enabled_{};
};

}  // namespace hundun::test::allocation_probe
