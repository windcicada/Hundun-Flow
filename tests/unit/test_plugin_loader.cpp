// SPDX-License-Identifier: Apache-2.0

#include "hundun/sdk/plugin_loader.hpp"
#include "runtime/include/hundun/runtime/error.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <optional>
#include <string_view>
#include <thread>
#include <utility>

namespace {

using hundun::sdk::PluginLoader;

constexpr std::string_view kExpectedName = "mock-plugin";
constexpr std::string_view kExpectedVersion = "1.2.3";
constexpr std::string_view kExpectedDescription =
    "metadata discovery fixture";
constexpr std::size_t kV1PrefixSize =
    offsetof(HundunPluginDescriptorV1, description) +
    sizeof(HundunPluginDescriptorV1::description);

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
    expect_runtime_error(path);
  }

  const auto nonexistent =
      paths.valid.parent_path() / "hundun-plugin-does-not-exist.so";
  HUNDUN_CHECK(std::filesystem::exists(nonexistent) == false);
  expect_runtime_error(nonexistent);
}

void test_independent_lifetimes_and_move(const std::filesystem::path& path) {
  PluginLoader first(path);
  {
    PluginLoader second(path);
    check_expected_metadata(first.descriptor());
    check_expected_metadata(second.descriptor());
  }
  check_expected_metadata(first.descriptor());

  std::optional<PluginLoader> destination;
  {
    PluginLoader source(path);
    destination.emplace(std::move(source));
  }
  check_expected_metadata(destination->descriptor());
  destination.reset();

  for (int repetition = 0; repetition < 20; ++repetition) {
    PluginLoader repeated(path);
    check_expected_metadata(repeated.descriptor());
  }
}

void test_parallel_loads(const std::filesystem::path& path) {
  constexpr std::size_t kThreadCount = 4;
  constexpr int kLoadsPerThread = 10;
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
    test_independent_lifetimes_and_move(paths.valid);
    test_parallel_loads(paths.valid);
  });
}
