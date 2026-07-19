// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdlib>
#include <new>

namespace hundun::test::allocation_probe {

inline thread_local bool count_attempts = false;
inline thread_local std::size_t attempt_count = 0U;

inline void record_attempt() noexcept {
  if (count_attempts) {
    ++attempt_count;
  }
}

inline void* allocate(std::size_t bytes) {
  record_attempt();
  if (void* pointer = std::malloc(bytes == 0U ? 1U : bytes)) {
    return pointer;
  }
  throw std::bad_alloc();
}

inline void* allocate_aligned(std::size_t bytes, std::size_t alignment) {
  record_attempt();
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

}  // namespace hundun::test::allocation_probe
