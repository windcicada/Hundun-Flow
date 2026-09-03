// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/sdk_plugin_loader.hpp"
#include "hundun/rt_error.hpp"
#include "tests/support/test_main.hpp"

#include <dlfcn.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstring>
#include <exception>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

namespace {

std::atomic<std::size_t> dlclose_calls{0};

enum class LookupCallState {
  disarmed,
  expect_prelookup_dlerror,
  expect_entry_dlsym,
  expect_postlookup_dlerror,
  complete,
  dlsym_before_prelookup_dlerror,
  invalid,
};

std::atomic<LookupCallState> lookup_call_state{LookupCallState::disarmed};

void observe_dlerror_call() noexcept {
  const auto state = lookup_call_state.load(std::memory_order_relaxed);
  if (state == LookupCallState::expect_prelookup_dlerror) {
    lookup_call_state.store(LookupCallState::expect_entry_dlsym,
                            std::memory_order_relaxed);
  } else if (state == LookupCallState::expect_postlookup_dlerror) {
    lookup_call_state.store(LookupCallState::complete,
                            std::memory_order_relaxed);
  } else if (state != LookupCallState::disarmed &&
             state != LookupCallState::dlsym_before_prelookup_dlerror &&
             state != LookupCallState::invalid) {
    lookup_call_state.store(LookupCallState::invalid,
                            std::memory_order_relaxed);
  }
}

void observe_dlsym_call(const char* symbol) noexcept {
  const auto state = lookup_call_state.load(std::memory_order_relaxed);
  if (state == LookupCallState::expect_entry_dlsym && symbol != nullptr &&
      std::strcmp(symbol, "hundun_plugin_entry_v1") == 0) {
    lookup_call_state.store(LookupCallState::expect_postlookup_dlerror,
                            std::memory_order_relaxed);
  } else if (state == LookupCallState::expect_prelookup_dlerror) {
    lookup_call_state.store(LookupCallState::dlsym_before_prelookup_dlerror,
                            std::memory_order_relaxed);
  } else if (state != LookupCallState::disarmed &&
             state != LookupCallState::dlsym_before_prelookup_dlerror &&
             state != LookupCallState::invalid) {
    lookup_call_state.store(LookupCallState::invalid,
                            std::memory_order_relaxed);
  }
}

}  // namespace

extern "C" int __real_dlclose(void* handle) noexcept;
extern "C" char* __real_dlerror() noexcept;
extern "C" void* __real_dlsym(void* handle, const char* symbol) noexcept;

extern "C" int __wrap_dlclose(void* handle) noexcept {
  dlclose_calls.fetch_add(1, std::memory_order_relaxed);
  return __real_dlclose(handle);
}

extern "C" char* __wrap_dlerror() noexcept {
  char* result = __real_dlerror();
  observe_dlerror_call();
  return result;
}

extern "C" void* __wrap_dlsym(void* handle, const char* symbol) noexcept {
  void* result = __real_dlsym(handle, symbol);
  observe_dlsym_call(symbol);
  return result;
}

namespace {

using hundun::sdk::PluginLoader;

constexpr std::string_view kExpectedName = "mock-plugin";
constexpr std::string_view kExpectedVersion = "1.2.3";
constexpr std::string_view kExpectedDescription =
    "metadata discovery fixture";
constexpr std::size_t kV1PrefixSize =
    offsetof(HundunPluginDescriptorV1, description) +
    sizeof(HundunPluginDescriptorV1::description);

std::size_t close_count() noexcept {
  return dlclose_calls.load(std::memory_order_relaxed);
}

void check_close_delta(std::size_t before, std::size_t expected,
                       const char* scenario) {
  const auto actual = close_count() - before;
  if (actual != expected) {
    throw std::runtime_error(std::string("dlclose count mismatch for ") +
                             scenario + ": expected " +
                             std::to_string(expected) + ", actual " +
                             std::to_string(actual));
  }
}

struct PluginPaths {
  std::filesystem::path valid;
  std::filesystem::path larger;
  std::filesystem::path overlap;
  std::filesystem::path short_object;
  std::filesystem::path range_only;
  std::filesystem::path no_overlap;
  std::filesystem::path reversed;
  std::filesystem::path capability;
  std::filesystem::path null_entry;
  std::filesystem::path null_name;
  std::filesystem::path empty_name;
  std::filesystem::path null_version;
  std::filesystem::path empty_version;
  std::filesystem::path null_description;
  std::filesystem::path empty_description;
  std::filesystem::path missing_entry;
};

void check_expected_metadata(const HundunPluginDescriptorV1& descriptor) {
  HUNDUN_CHECK(descriptor.abi_version_min == 1u);
  HUNDUN_CHECK(descriptor.abi_version_max >= 1u);
  HUNDUN_CHECK(descriptor.capability_flags == UINT64_C(0));
  HUNDUN_CHECK(std::string_view(descriptor.name) == kExpectedName);
  HUNDUN_CHECK(std::string_view(descriptor.version) == kExpectedVersion);
  HUNDUN_CHECK(std::string_view(descriptor.description) ==
               kExpectedDescription);
}

void expect_runtime_error(const std::filesystem::path& path) {
  bool threw = false;
  try {
    PluginLoader loader(path);
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string_view(error.what()).empty() == false);
  }
  HUNDUN_CHECK(threw);
}

void test_accepted_descriptors(const PluginPaths& paths) {
  const auto before = close_count();
  {
    PluginLoader exact(paths.valid);
    HUNDUN_CHECK(exact.descriptor().struct_size == kV1PrefixSize);
    HUNDUN_CHECK(exact.descriptor().abi_version_max == 1u);
    check_expected_metadata(exact.descriptor());

    PluginLoader larger(paths.larger);
    HUNDUN_CHECK(larger.descriptor().struct_size > kV1PrefixSize);
    check_expected_metadata(larger.descriptor());

    PluginLoader overlap(paths.overlap);
    HUNDUN_CHECK(overlap.descriptor().abi_version_min == 1u);
    HUNDUN_CHECK(overlap.descriptor().abi_version_max == 2u);
    check_expected_metadata(overlap.descriptor());
  }
  check_close_delta(before, 3, "accepted descriptors");
}

void test_rejected_descriptors(const PluginPaths& paths) {
  const std::array<std::filesystem::path, 13> rejected = {
      paths.short_object,
      paths.range_only,
      paths.no_overlap,
      paths.reversed,
      paths.capability,
      paths.null_entry,
      paths.null_name,
      paths.empty_name,
      paths.null_version,
      paths.empty_version,
      paths.null_description,
      paths.empty_description,
      paths.missing_entry};
  for (const auto& path : rejected) {
    const auto before = close_count();
    expect_runtime_error(path);
    check_close_delta(before, 1, "DSO-backed rejection");
  }

  const auto nonexistent =
      paths.valid.parent_path() / "hundun-plugin-does-not-exist.so";
  HUNDUN_CHECK(std::filesystem::exists(nonexistent) == false);
  const auto before = close_count();
  expect_runtime_error(nonexistent);
  check_close_delta(before, 0, "nonexistent-path rejection");
}

void test_independent_lifetimes(const std::filesystem::path& path) {
  const auto before = close_count();
  {
    PluginLoader first(path);
    PluginLoader second(path);
    check_expected_metadata(first.descriptor());
    check_expected_metadata(second.descriptor());
    check_close_delta(before, 0, "two live independent loaders");
  }
  check_close_delta(before, 2, "two destroyed independent loaders");
}

void test_move_lifetime(const std::filesystem::path& path) {
  const auto before = close_count();
  std::optional<PluginLoader> destination;
  {
    PluginLoader source(path);
    destination.emplace(std::move(source));
  }
  check_close_delta(before, 0, "destroyed move source");
  check_expected_metadata(destination->descriptor());
  destination.reset();
  check_close_delta(before, 1, "destroyed move destination");
}

void test_repeated_loads(const std::filesystem::path& path) {
  const auto before = close_count();
  for (int repetition = 0; repetition < 20; ++repetition) {
    PluginLoader repeated(path);
    check_expected_metadata(repeated.descriptor());
  }
  check_close_delta(before, 20, "sequential loads");
}

void test_lookup_call_order(const std::filesystem::path& path) {
  const auto before = close_count();
  const auto previous = lookup_call_state.exchange(
      LookupCallState::expect_prelookup_dlerror, std::memory_order_relaxed);
  HUNDUN_CHECK(previous == LookupCallState::disarmed);
  {
    PluginLoader loader(path);
    const auto observed = lookup_call_state.exchange(
        LookupCallState::disarmed, std::memory_order_relaxed);
    if (observed == LookupCallState::dlsym_before_prelookup_dlerror) {
      throw std::runtime_error(
          "plugin lookup call order began with dlsym instead of the "
          "pre-lookup dlerror");
    }
    if (observed != LookupCallState::complete) {
      throw std::runtime_error(
          "plugin lookup did not follow dlerror, entry dlsym, dlerror order");
    }
    check_expected_metadata(loader.descriptor());
  }
  check_close_delta(before, 1, "lookup-order observation load");
}

void test_parallel_loads(const std::filesystem::path& path) {
  constexpr std::size_t kThreadCount = 4;
  constexpr int kLoadsPerThread = 10;
  const auto before = close_count();
  std::array<std::exception_ptr, kThreadCount> errors{};
  std::array<std::thread, kThreadCount> threads;

  for (std::size_t index = 0; index < threads.size(); ++index) {
    threads[index] = std::thread([&, index] {
      try {
        for (int load = 0; load < kLoadsPerThread; ++load) {
          PluginLoader loader(path);
          check_expected_metadata(loader.descriptor());
        }
      } catch (...) {
        errors[index] = std::current_exception();
      }
    });
  }
  for (auto& thread : threads) {
    thread.join();
  }
  for (const auto& error : errors) {
    if (error) {
      std::rethrow_exception(error);
    }
  }
  check_close_delta(before, 40, "parallel loads");
}

PluginPaths parse_paths(int argc, char** argv) {
  constexpr int kExpectedArgumentCount = 17;
  HUNDUN_CHECK(argc == kExpectedArgumentCount);
  return PluginPaths{argv[1],  argv[2],  argv[3],  argv[4],
                     argv[5],  argv[6],  argv[7],  argv[8],
                     argv[9],  argv[10], argv[11], argv[12],
                     argv[13], argv[14], argv[15], argv[16]};
}

}  // namespace

int main(int argc, char** argv) {
  return hundun::test::run([&] {
    const auto paths = parse_paths(argc, argv);
    test_accepted_descriptors(paths);
    test_rejected_descriptors(paths);
    test_independent_lifetimes(paths.valid);
    test_move_lifetime(paths.valid);
    test_repeated_loads(paths.valid);
    test_lookup_call_order(paths.valid);
    test_parallel_loads(paths.valid);
  });
}
