// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
// Test-only LD_PRELOAD observer at POSIX/MPI boundaries. No production hooks.
#include <mpi.h>
#include <fcntl.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <execinfo.h>

namespace {
bool active = false, started = false, fired = false;
long calls = 0, failure = -1;
int evidence = -1, complete = -1;
std::uint64_t accepted_step = 0;
bool visit_mode() noexcept {
  return std::getenv("HUNDUN_TEST_VISIT_ALLOC") != nullptr;
}
int rank_value() noexcept {
  const char* value = std::getenv("OMPI_COMM_WORLD_RANK");
  if (!value) value = std::getenv("PMI_RANK");
  return value ? std::atoi(value) : -1;
}
bool target() noexcept {
  const char* value = std::getenv("HUNDUN_TEST_ALLOC_RANK");
  return value && rank_value() == std::atoi(value);
}
bool ends(const char* value, const char* suffix) noexcept {
  const auto n = std::strlen(value), m = std::strlen(suffix);
  return n >= m && std::strcmp(value + n - m, suffix) == 0;
}
void start() noexcept {
  if (started || !target()) return;
  const char* value = std::getenv("HUNDUN_TEST_ALLOC_INDEX");
  failure = value ? std::atol(value) : -1;
  active = started = true;
}
}

void* operator new(std::size_t bytes) {
  if (active && calls++ == failure) {
    fired = true;
    if (std::getenv("HUNDUN_TEST_ALLOC_TRACE")) {
      void* trace[24];
      const int count = ::backtrace(trace, 24);
      ::backtrace_symbols_fd(trace, count, STDERR_FILENO);
    }
    throw std::bad_alloc{};
  }
  if (void* value = std::malloc(bytes ? bytes : 1U)) return value;
  throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* value) noexcept { std::free(value); }
void operator delete[](void* value) noexcept { std::free(value); }
void operator delete(void* value, std::size_t) noexcept { std::free(value); }
void operator delete[](void* value, std::size_t) noexcept { std::free(value); }

extern "C" int open(const char* path, int flags, ...) {
  mode_t permissions = 0;
  if ((flags & O_CREAT) != 0) {
    va_list args;
    va_start(args, flags);
    permissions = va_arg(args, int);
    va_end(args);
  }
  const int fd = static_cast<int>(::syscall(SYS_openat, AT_FDCWD, path, flags, permissions));
  if (target() && fd >= 0) {
    if (ends(path, "/evidence.jsonl")) evidence = fd;
    if (ends(path, ".complete")) complete = fd;
  }
  return fd;
}
extern "C" int fsync(int fd) {
  const int result = static_cast<int>(::syscall(SYS_fsync, fd));
  if (result == 0 && fd == evidence) {
    if (visit_mode()) active = false;
    else start();
  }
  if (result == 0 && fd == complete) active = false;
  return result;
}
extern "C" int close(int fd) {
  if (fd == evidence) evidence = -1;
  if (fd == complete) complete = -1;
  return static_cast<int>(::syscall(SYS_close, fd));
}
extern "C" int MPI_Reduce(const void* send, void* receive, int count,
                          MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm) {
  if (!visit_mode() && rank_value() != 0) start();
  return PMPI_Reduce(send, receive, count, type, op, root, comm);
}
extern "C" void hundun_v04_observe_step(std::uint64_t step, int entering) {
  if (visit_mode() && entering == 0) {
    accepted_step = step;
    start();
  }
}
extern "C" int MPI_Finalize() {
  active = false;
  if (target())
    std::fprintf(stderr, "CHECKPOINT_ALLOC sites=%ld fired=%d started=%d\n",
                  calls, fired, started);
  if (visit_mode())
    std::fprintf(stderr, "VISIT_FINALIZE rank=%d accepted_step=%llu\n",
                 rank_value(), static_cast<unsigned long long>(accepted_step));
  return PMPI_Finalize();
}
