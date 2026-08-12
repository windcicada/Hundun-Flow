// SPDX-License-Identifier: Apache-2.0

#include "tests/support/allocation_attempt_guard.hpp"

void* operator new(std::size_t bytes) {
  return hundun::test::allocation_probe::allocate(bytes);
}
void* operator new[](std::size_t bytes) {
  return hundun::test::allocation_probe::allocate(bytes);
}
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
  try {
    return hundun::test::allocation_probe::allocate(bytes);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t bytes, const std::nothrow_t&) noexcept {
  try {
    return hundun::test::allocation_probe::allocate(bytes);
  } catch (...) {
    return nullptr;
  }
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  return hundun::test::allocation_probe::allocate_aligned(
      bytes, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
  return hundun::test::allocation_probe::allocate_aligned(
      bytes, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t bytes, std::align_val_t alignment,
                   const std::nothrow_t&) noexcept {
  try {
    return hundun::test::allocation_probe::allocate_aligned(
        bytes, static_cast<std::size_t>(alignment));
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t bytes, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
  try {
    return hundun::test::allocation_probe::allocate_aligned(
        bytes, static_cast<std::size_t>(alignment));
  } catch (...) {
    return nullptr;
  }
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t,
                     const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t,
                       const std::nothrow_t&) noexcept {
  std::free(pointer);
}
