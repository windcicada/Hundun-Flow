// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/restart_binary.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "runtime/src/restart_detail.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

using hundun::runtime::Box3;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Error;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::RestartMetadata;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;
using hundun::runtime::detail::RestartFailureInjection;
using hundun::runtime::detail::RestartFailurePhase;

constexpr Int3 kGlobalExtent{4, 2, 2};
constexpr std::array<bool, 3> kPeriodic{true, false, true};
constexpr double kTimeScale = 0.0625;

template <class T>
T load_native(const std::vector<std::byte> &bytes, std::size_t offset) {
  HUNDUN_CHECK(offset <= bytes.size());
  HUNDUN_CHECK(sizeof(T) <= bytes.size() - offset);
  T value{};
  std::memcpy(&value, bytes.data() + offset, sizeof(T));
  return value;
}

template <class T>
void store_native(std::vector<std::byte> &bytes, std::size_t offset, T value) {
  HUNDUN_CHECK(offset <= bytes.size());
  HUNDUN_CHECK(sizeof(T) <= bytes.size() - offset);
  std::memcpy(bytes.data() + offset, &value, sizeof(T));
}

std::vector<std::byte> read_bytes(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  HUNDUN_CHECK(input.is_open());
  input.seekg(0, std::ios::end);
  const auto end = input.tellg();
  HUNDUN_CHECK(end >= 0);
  input.seekg(0, std::ios::beg);
  std::vector<std::byte> bytes(static_cast<std::size_t>(end));
  if (!bytes.empty()) {
    input.read(reinterpret_cast<char *>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
  }
  HUNDUN_CHECK(input.good() || input.eof());
  return bytes;
}

void write_bytes(const std::filesystem::path &path,
                 const std::vector<std::byte> &bytes) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  HUNDUN_CHECK(output.is_open());
  if (!bytes.empty()) {
    output.write(reinterpret_cast<const char *>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
  }
  output.flush();
  HUNDUN_CHECK(output.good());
  output.close();
  HUNDUN_CHECK(!output.fail());
}

std::uint64_t reference_crc(const std::vector<std::byte> &bytes) {
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  std::uint64_t crc = 0;
  for (const std::byte byte : bytes) {
    crc ^= static_cast<std::uint64_t>(std::to_integer<unsigned char>(byte))
           << 56U;
    for (int bit = 0; bit < 8; ++bit) {
      crc = (crc & (UINT64_C(1) << 63U)) != 0U ? (crc << 1U) ^ polynomial
                                               : crc << 1U;
    }
  }
  return crc;
}

void append_native(std::vector<std::byte> &bytes, const void *value,
                   std::size_t size) {
  const auto *begin = static_cast<const std::byte *>(value);
  bytes.insert(bytes.end(), begin, begin + size);
}

template <class T> void append_native(std::vector<std::byte> &bytes, T value) {
  append_native(bytes, &value, sizeof(T));
}

void rewrite_marker(const std::filesystem::path &directory) {
  const auto manifest = read_bytes(directory / "manifest.v1.bin");
  std::vector<std::byte> marker;
  const std::array<char, 8> magic{'H', 'U', 'N', 'D', 'C', 'M', 'P', '1'};
  append_native(marker, magic.data(), magic.size());
  append_native(marker, std::uint32_t{1});
  append_native(marker, UINT32_C(0x01020304));
  append_native(marker, static_cast<std::uint64_t>(manifest.size()));
  append_native(marker, reference_crc(manifest));
  write_bytes(directory / "COMPLETED", marker);
}

std::filesystem::path broadcast_root(const MpiContext &context) {
  std::string text;
  if (context.rank() == 0) {
    text = "/tmp/hundun-task10-restart-" +
           std::to_string(static_cast<long long>(::getpid()));
  }
  std::uint64_t length = static_cast<std::uint64_t>(text.size());
  HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, context.comm()) ==
               MPI_SUCCESS);
  if (context.rank() != 0) {
    text.resize(static_cast<std::size_t>(length));
  }
  HUNDUN_CHECK(length <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(length), MPI_BYTE, 0,
                         context.comm()) == MPI_SUCCESS);
  return text;
}

std::string rank_filename(int rank) {
  std::array<char, 64> buffer{};
  const int count =
      std::snprintf(buffer.data(), buffer.size(), "restart.rank%06d.bin", rank);
  HUNDUN_CHECK(count > 0);
  HUNDUN_CHECK(static_cast<std::size_t>(count) < buffer.size());
  return std::string(buffer.data(), static_cast<std::size_t>(count));
}

std::filesystem::path step_directory(const std::filesystem::path &root,
                                     std::int64_t step) {
  std::array<char, 64> buffer{};
  const int count = std::snprintf(buffer.data(), buffer.size(), "step%08lld",
                                  static_cast<long long>(step));
  HUNDUN_CHECK(count > 0);
  HUNDUN_CHECK(static_cast<std::size_t>(count) < buffer.size());
  return root / std::string(buffer.data(), static_cast<std::size_t>(count));
}

struct Fields final {
  FieldRegistry registry;
  FieldId real{};
  FieldId integer{};
  FieldId transient{};
  FieldStorage storage;

  Fields(Int3 extent, std::string real_name = "state")
      : real(registry.declare_field(FieldDescriptor{
            std::move(real_name), "1", "checkpoint-test",
            FunctionSpace::cell_average, ScalarType::float64, 2U, 1, true,
            RestartPolicy::persistent, OutputPolicy::selected})),
        integer(registry.declare_field(FieldDescriptor{
            "tag", "1", "checkpoint-test", FunctionSpace::cell_average,
            ScalarType::int32, 1U, 1, false, RestartPolicy::persistent,
            OutputPolicy::never})),
        transient(registry.declare_field(FieldDescriptor{
            "scratch", "1", "checkpoint-test", FunctionSpace::cell_average,
            ScalarType::float64, 1U, 1, false, RestartPolicy::transient,
            OutputPolicy::never})),
        storage(freeze(), extent) {}

private:
  FieldRegistry &freeze() {
    registry.freeze();
    return registry;
  }
};

enum class WriterMismatch { none, first_components, first_ghost_width };

struct WriterLayout final {
  FieldRegistry registry;
  FieldId real{};
  FieldId integer{};
  FieldId transient{};
  FieldStorage storage;

  WriterLayout(Int3 extent, WriterMismatch mismatch)
      : real(registry.declare_field(FieldDescriptor{
            "state", "1", "checkpoint-test", FunctionSpace::cell_average,
            ScalarType::float64,
            mismatch == WriterMismatch::first_components ? 3U : 2U,
            mismatch == WriterMismatch::first_ghost_width ? 2 : 1, true,
            RestartPolicy::persistent, OutputPolicy::selected})),
        integer(registry.declare_field(FieldDescriptor{
            "tag", "1", "checkpoint-test", FunctionSpace::cell_average,
            ScalarType::int32, 1U, 1, false, RestartPolicy::persistent,
            OutputPolicy::never})),
        transient(registry.declare_field(FieldDescriptor{
            "scratch", "1", "checkpoint-test", FunctionSpace::cell_average,
            ScalarType::float64, 1U, 1, false, RestartPolicy::transient,
            OutputPolicy::never})),
        storage(freeze(), extent) {}

private:
  FieldRegistry &freeze() {
    registry.freeze();
    return registry;
  }
};

enum class DestinationMismatch { none, first_scalar_type, later_scalar_type };

struct DestinationLayout final {
  DestinationMismatch mismatch;
  FieldRegistry registry;
  FieldId real{};
  FieldId integer{};
  FieldId transient{};
  FieldStorage storage;

  DestinationLayout(Int3 extent, DestinationMismatch mismatch_value)
      : mismatch(mismatch_value),
        real(registry.declare_field(FieldDescriptor{
            "state", "1", "checkpoint-test", FunctionSpace::cell_average,
            mismatch == DestinationMismatch::first_scalar_type
                ? ScalarType::int32
                : ScalarType::float64,
            2U, 1, true, RestartPolicy::persistent, OutputPolicy::selected})),
        integer(registry.declare_field(FieldDescriptor{
            "tag", "1", "checkpoint-test", FunctionSpace::cell_average,
            mismatch == DestinationMismatch::later_scalar_type
                ? ScalarType::float64
                : ScalarType::int32,
            1U, 1, false, RestartPolicy::persistent, OutputPolicy::never})),
        transient(registry.declare_field(FieldDescriptor{
            "scratch", "1", "checkpoint-test", FunctionSpace::cell_average,
            ScalarType::float64, 1U, 1, false, RestartPolicy::transient,
            OutputPolicy::never})),
        storage(freeze(), extent) {}

private:
  FieldRegistry &freeze() {
    registry.freeze();
    return registry;
  }
};

void fill_destination_sentinels(DestinationLayout &destination) {
  const Int3 extent = destination.storage.interior_extent();
  if (destination.mismatch == DestinationMismatch::first_scalar_type) {
    auto real = destination.storage.view<std::int32_t>(destination.real);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          real(i, j, k, 0) = 811;
          real(i, j, k, 1) = 812;
        }
      }
    }
  } else {
    auto real = destination.storage.view<double>(destination.real);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          real(i, j, k, 0) = 811.0;
          real(i, j, k, 1) = 812.0;
        }
      }
    }
  }

  if (destination.mismatch == DestinationMismatch::later_scalar_type) {
    auto integer = destination.storage.view<double>(destination.integer);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          integer(i, j, k, 0) = 821.0;
        }
      }
    }
  } else {
    auto integer = destination.storage.view<std::int32_t>(destination.integer);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          integer(i, j, k, 0) = 821;
        }
      }
    }
  }

  auto transient = destination.storage.view<double>(destination.transient);
  for (int k = -1; k <= extent.z; ++k) {
    for (int j = -1; j <= extent.y; ++j) {
      for (int i = -1; i <= extent.x; ++i) {
        transient(i, j, k, 0) = 831.0;
      }
    }
  }
}

void verify_destination_sentinels(const DestinationLayout &destination) {
  const Int3 extent = destination.storage.interior_extent();
  if (destination.mismatch == DestinationMismatch::first_scalar_type) {
    const auto real = destination.storage.view<std::int32_t>(destination.real);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          HUNDUN_CHECK(real(i, j, k, 0) == 811);
          HUNDUN_CHECK(real(i, j, k, 1) == 812);
        }
      }
    }
  } else {
    const auto real = destination.storage.view<double>(destination.real);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          HUNDUN_CHECK(real(i, j, k, 0) == 811.0);
          HUNDUN_CHECK(real(i, j, k, 1) == 812.0);
        }
      }
    }
  }

  if (destination.mismatch == DestinationMismatch::later_scalar_type) {
    const auto integer = destination.storage.view<double>(destination.integer);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          HUNDUN_CHECK(integer(i, j, k, 0) == 821.0);
        }
      }
    }
  } else {
    const auto integer =
        destination.storage.view<std::int32_t>(destination.integer);
    for (int k = -1; k <= extent.z; ++k) {
      for (int j = -1; j <= extent.y; ++j) {
        for (int i = -1; i <= extent.x; ++i) {
          HUNDUN_CHECK(integer(i, j, k, 0) == 821);
        }
      }
    }
  }

  const auto transient =
      destination.storage.view<double>(destination.transient);
  for (int k = -1; k <= extent.z; ++k) {
    for (int j = -1; j <= extent.y; ++j) {
      for (int i = -1; i <= extent.x; ++i) {
        HUNDUN_CHECK(transient(i, j, k, 0) == 831.0);
      }
    }
  }
}

void fill_source(Fields &fields, int rank) {
  const Int3 extent = fields.storage.interior_extent();
  auto real = fields.storage.view<double>(fields.real);
  auto integer = fields.storage.view<std::int32_t>(fields.integer);
  auto transient = fields.storage.view<double>(fields.transient);
  for (int k = -1; k <= extent.z; ++k) {
    for (int j = -1; j <= extent.y; ++j) {
      for (int i = -1; i <= extent.x; ++i) {
        transient(i, j, k, 0) = -3000.0;
        integer(i, j, k, 0) = -2000;
        real(i, j, k, 0) = -1000.0;
        real(i, j, k, 1) = -1000.0;
      }
    }
  }
  int ordinal = 0;
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        real(i, j, k, 0) = rank * 1000.0 + ordinal + 0.25;
        real(i, j, k, 1) = rank * 1000.0 + ordinal + 0.75;
        integer(i, j, k, 0) = rank * 1000 + ordinal;
        ++ordinal;
      }
    }
  }
}

void fill_sentinel(Fields &fields) {
  const Int3 extent = fields.storage.interior_extent();
  auto real = fields.storage.view<double>(fields.real);
  auto integer = fields.storage.view<std::int32_t>(fields.integer);
  auto transient = fields.storage.view<double>(fields.transient);
  for (int k = -1; k <= extent.z; ++k) {
    for (int j = -1; j <= extent.y; ++j) {
      for (int i = -1; i <= extent.x; ++i) {
        real(i, j, k, 0) = 701.0;
        real(i, j, k, 1) = 702.0;
        integer(i, j, k, 0) = 703;
        transient(i, j, k, 0) = 704.0;
      }
    }
  }
}

void verify_unchanged(const Fields &fields) {
  const Int3 extent = fields.storage.interior_extent();
  const auto real = fields.storage.view<double>(fields.real);
  const auto integer = fields.storage.view<std::int32_t>(fields.integer);
  const auto transient = fields.storage.view<double>(fields.transient);
  for (int k = -1; k <= extent.z; ++k) {
    for (int j = -1; j <= extent.y; ++j) {
      for (int i = -1; i <= extent.x; ++i) {
        HUNDUN_CHECK(real(i, j, k, 0) == 701.0);
        HUNDUN_CHECK(real(i, j, k, 1) == 702.0);
        HUNDUN_CHECK(integer(i, j, k, 0) == 703);
        HUNDUN_CHECK(transient(i, j, k, 0) == 704.0);
      }
    }
  }
}

void verify_restored(const Fields &fields, int rank) {
  const Int3 extent = fields.storage.interior_extent();
  const auto real = fields.storage.view<double>(fields.real);
  const auto integer = fields.storage.view<std::int32_t>(fields.integer);
  const auto transient = fields.storage.view<double>(fields.transient);
  int ordinal = 0;
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        HUNDUN_CHECK(real(i, j, k, 0) == rank * 1000.0 + ordinal + 0.25);
        HUNDUN_CHECK(real(i, j, k, 1) == rank * 1000.0 + ordinal + 0.75);
        HUNDUN_CHECK(integer(i, j, k, 0) == rank * 1000 + ordinal);
        HUNDUN_CHECK(transient(i, j, k, 0) == 704.0);
        ++ordinal;
      }
    }
  }
  for (const Int3 corner :
       {Int3{-1, -1, -1}, Int3{extent.x, extent.y, extent.z}}) {
    HUNDUN_CHECK(real(corner.x, corner.y, corner.z, 0) == 701.0);
    HUNDUN_CHECK(integer(corner.x, corner.y, corner.z, 0) == 703);
    HUNDUN_CHECK(transient(corner.x, corner.y, corner.z, 0) == 704.0);
  }
}

template <class Function>
std::string expect_same_collective_error(const MpiContext &context,
                                         Function &&function) {
  std::string message;
  bool caught = false;
  try {
    function();
  } catch (const Error &error) {
    message = error.what();
    caught = !message.empty();
  }
  HUNDUN_CHECK(caught);
  std::array<char, 512> local{};
  HUNDUN_CHECK(message.size() < local.size());
  std::copy(message.begin(), message.end(), local.begin());
  std::vector<char> all(local.size() *
                        static_cast<std::size_t>(context.size()));
  HUNDUN_CHECK(MPI_Allgather(local.data(), static_cast<int>(local.size()),
                             MPI_CHAR, all.data(),
                             static_cast<int>(local.size()), MPI_CHAR,
                             context.comm()) == MPI_SUCCESS);
  for (int rank = 0; rank < context.size(); ++rank) {
    const char *other =
        all.data() + static_cast<std::size_t>(rank) * local.size();
    HUNDUN_CHECK(message == other);
  }
  return message;
}

void validate_completed_files(const MpiContext &context,
                              const std::filesystem::path &directory) {
  if (context.rank() != 0) {
    return;
  }
  std::set<std::string> names;
  for (const auto &entry : std::filesystem::directory_iterator(directory)) {
    names.insert(entry.path().filename().string());
    HUNDUN_CHECK(entry.path().extension() != ".tmp");
  }
  HUNDUN_CHECK(names.count("manifest.v1.bin") == 1U);
  HUNDUN_CHECK(names.count("COMPLETED") == 1U);
  HUNDUN_CHECK(names.size() == static_cast<std::size_t>(context.size()) + 2U);
  for (int rank = 0; rank < context.size(); ++rank) {
    HUNDUN_CHECK(names.count(rank_filename(rank)) == 1U);
  }

  const auto marker = read_bytes(directory / "COMPLETED");
  HUNDUN_CHECK(marker.size() == 32U);
  HUNDUN_CHECK(std::memcmp(marker.data(), "HUNDCMP1", 8U) == 0);
  HUNDUN_CHECK(load_native<std::uint32_t>(marker, 8U) == 1U);
  HUNDUN_CHECK(load_native<std::uint32_t>(marker, 12U) == UINT32_C(0x01020304));
  const auto manifest = read_bytes(directory / "manifest.v1.bin");
  HUNDUN_CHECK(load_native<std::uint64_t>(marker, 16U) == manifest.size());
  HUNDUN_CHECK(load_native<std::uint64_t>(marker, 24U) ==
               reference_crc(manifest));
  HUNDUN_CHECK(std::memcmp(manifest.data(), "HUNDMAN1", 8U) == 0);
  HUNDUN_CHECK(load_native<std::int32_t>(manifest, 32U) == context.size());
  HUNDUN_CHECK(load_native<std::uint32_t>(manifest, 75U) ==
               static_cast<std::uint32_t>(context.size()));

  std::size_t cursor = 79U;
  for (int rank = 0; rank < context.size(); ++rank) {
    HUNDUN_CHECK(load_native<std::int32_t>(manifest, cursor) == rank);
    cursor += 28U;
    const auto name_size = load_native<std::uint32_t>(manifest, cursor);
    cursor += 4U;
    const std::string name(
        reinterpret_cast<const char *>(manifest.data() + cursor), name_size);
    HUNDUN_CHECK(name == rank_filename(rank));
    cursor += name_size;
    const auto byte_size = load_native<std::uint64_t>(manifest, cursor);
    cursor += 8U;
    const auto crc = load_native<std::uint64_t>(manifest, cursor);
    cursor += 8U;
    const auto rank_bytes = read_bytes(directory / name);
    HUNDUN_CHECK(byte_size == rank_bytes.size());
    HUNDUN_CHECK(crc == reference_crc(rank_bytes));
  }
  HUNDUN_CHECK(cursor == manifest.size());
}

void test_success(const MpiContext &context,
                  const StructuredDecomposition &decomposition,
                  const std::filesystem::path &root) {
  Fields source(decomposition.local_extent());
  fill_source(source, context.rank());
  const auto directory = step_directory(root, 10);
  hundun::runtime::write_restart_checkpoint(context, decomposition,
                                            source.registry, source.storage,
                                            directory, 10, 10.0 * kTimeScale);
  validate_completed_files(context, directory);
  context.barrier();

  Fields destination(decomposition.local_extent());
  fill_sentinel(destination);
  const RestartMetadata metadata = hundun::runtime::read_restart_checkpoint(
      context, decomposition, destination.registry, destination.storage,
      directory);
  HUNDUN_CHECK(metadata.step == 10);
  HUNDUN_CHECK(metadata.time_s == 10.0 * kTimeScale);
  verify_restored(destination, context.rank());
}

void test_preflight_rejection(const MpiContext &context,
                              const StructuredDecomposition &decomposition,
                              const std::filesystem::path &root) {
  Fields fields(decomposition.local_extent());
  fill_source(fields, context.rank());

  if (context.size() > 1) {
    const auto local_directory =
        step_directory(root, context.rank() == 0 ? 20 : 21);
    expect_same_collective_error(context, [&] {
      hundun::runtime::write_restart_checkpoint(
          context, decomposition, fields.registry, fields.storage,
          local_directory, context.rank() == 0 ? 20 : 21, 20.0 * kTimeScale);
    });
    context.barrier();
    if (context.rank() == 0) {
      HUNDUN_CHECK(!std::filesystem::exists(step_directory(root, 20)));
      HUNDUN_CHECK(!std::filesystem::exists(step_directory(root, 21)));
    }

    const auto step_path = step_directory(root, 22);
    expect_same_collective_error(context, [&] {
      hundun::runtime::write_restart_checkpoint(
          context, decomposition, fields.registry, fields.storage, step_path,
          context.rank() == 0 ? 22 : 23, 22.0 * kTimeScale);
    });
    if (context.rank() == 0) {
      HUNDUN_CHECK(!std::filesystem::exists(step_path));
    }

    const auto time_path = step_directory(root, 24);
    const double time = context.rank() == 0 ? 0.0 : -0.0;
    expect_same_collective_error(context, [&] {
      hundun::runtime::write_restart_checkpoint(context, decomposition,
                                                fields.registry, fields.storage,
                                                time_path, 24, time);
    });
    if (context.rank() == 0) {
      HUNDUN_CHECK(!std::filesystem::exists(time_path));
    }

    Fields schema_fields(decomposition.local_extent(),
                         context.rank() == 0 ? "state" : "other-state");
    fill_source(schema_fields, context.rank());
    const auto schema_path = step_directory(root, 25);
    expect_same_collective_error(context, [&] {
      hundun::runtime::write_restart_checkpoint(
          context, decomposition, schema_fields.registry, schema_fields.storage,
          schema_path, 25, 25.0 * kTimeScale);
    });
    if (context.rank() == 0) {
      HUNDUN_CHECK(!std::filesystem::exists(schema_path));
    }
  }

  const auto existing = step_directory(root, 26);
  if (context.rank() == 0) {
    std::filesystem::create_directories(existing);
    write_bytes(existing / "sentinel", {std::byte{7}});
  }
  context.barrier();
  expect_same_collective_error(context, [&] {
    hundun::runtime::write_restart_checkpoint(context, decomposition,
                                              fields.registry, fields.storage,
                                              existing, 26, 26.0 * kTimeScale);
  });
  if (context.rank() == 0) {
    HUNDUN_CHECK(std::distance(std::filesystem::directory_iterator(existing),
                               std::filesystem::directory_iterator{}) == 1);
  }
}

void test_writer_layout_preflight(const MpiContext &context,
                                  const StructuredDecomposition &decomposition,
                                  const std::filesystem::path &root) {
  Fields checkpoint_fields(decomposition.local_extent());
  const std::array<WriterMismatch, 2> mismatches{
      WriterMismatch::first_components, WriterMismatch::first_ghost_width};
  for (std::size_t index = 0; index < mismatches.size(); ++index) {
    const WriterMismatch local_mismatch =
        context.size() > 1 && context.rank() == 0 ? WriterMismatch::none
                                                  : mismatches[index];
    WriterLayout writer(decomposition.local_extent(), local_mismatch);
    const std::int64_t step = 53 + static_cast<std::int64_t>(index);
    const auto directory = step_directory(root, step);
    expect_same_collective_error(context, [&] {
      hundun::runtime::write_restart_checkpoint(
          context, decomposition, checkpoint_fields.registry, writer.storage,
          directory, step, static_cast<double>(step) * kTimeScale);
    });
    context.barrier();
    if (context.rank() == 0) {
      HUNDUN_CHECK(!std::filesystem::exists(directory));
    }
  }
}

void test_injected_failures(const MpiContext &context,
                            const StructuredDecomposition &decomposition,
                            const std::filesystem::path &root) {
  Fields fields(decomposition.local_extent());
  fill_source(fields, context.rank());
  const std::array<RestartFailureInjection, 8> injections{
      RestartFailureInjection{RestartFailurePhase::path_preparation,
                              context.size() > 1 ? 1 : 0},
      RestartFailureInjection{RestartFailurePhase::agreement_preparation,
                              context.size() > 1 ? 1 : 0},
      RestartFailureInjection{RestartFailurePhase::owned_box_preparation,
                              context.size() > 1 ? 1 : 0},
      RestartFailureInjection{RestartFailurePhase::filename_preparation,
                              context.size() > 1 ? 1 : 0},
      RestartFailureInjection{RestartFailurePhase::rank_file,
                              context.size() > 1 ? 1 : 0},
      RestartFailureInjection{RestartFailurePhase::record_gather_preparation,
                              0},
      RestartFailureInjection{RestartFailurePhase::manifest, 0},
      RestartFailureInjection{RestartFailurePhase::marker, 0}};
  for (std::size_t index = 0; index < injections.size(); ++index) {
    const std::int64_t step = 30 + static_cast<std::int64_t>(index);
    const auto directory = step_directory(root, step);
    expect_same_collective_error(context, [&] {
      hundun::runtime::detail::write_restart_checkpoint_with_failure(
          context, decomposition, fields.registry, fields.storage, directory,
          step, static_cast<double>(step) * kTimeScale, injections[index]);
    });
    context.barrier();
    if (context.rank() == 0) {
      if (injections[index].phase == RestartFailurePhase::path_preparation ||
          injections[index].phase ==
              RestartFailurePhase::agreement_preparation ||
          injections[index].phase ==
              RestartFailurePhase::owned_box_preparation) {
        HUNDUN_CHECK(!std::filesystem::exists(directory));
      } else {
        HUNDUN_CHECK(std::filesystem::exists(directory));
        HUNDUN_CHECK(!std::filesystem::exists(directory / "COMPLETED"));
      }
      if (injections[index].phase ==
          RestartFailurePhase::record_gather_preparation) {
        for (int rank = 0; rank < context.size(); ++rank) {
          HUNDUN_CHECK(
              std::filesystem::exists(directory / rank_filename(rank)));
        }
      }
    }
  }
}

template <class Mutation>
void expect_transactional_read_failure(
    const MpiContext &context, const StructuredDecomposition &decomposition,
    const std::filesystem::path &root, std::int64_t step, Mutation &&mutation) {
  Fields source(decomposition.local_extent());
  fill_source(source, context.rank());
  const auto directory = step_directory(root, step);
  hundun::runtime::write_restart_checkpoint(
      context, decomposition, source.registry, source.storage, directory, step,
      static_cast<double>(step) * kTimeScale);
  context.barrier();
  if (context.rank() == 0) {
    mutation(directory);
  }
  context.barrier();

  Fields destination(decomposition.local_extent());
  fill_sentinel(destination);
  expect_same_collective_error(context, [&] {
    static_cast<void>(hundun::runtime::read_restart_checkpoint(
        context, decomposition, destination.registry, destination.storage,
        directory));
  });
  verify_unchanged(destination);
}

void mutate_manifest_and_rebind(
    const std::filesystem::path &directory, std::size_t offset,
    const std::function<void(std::vector<std::byte> &, std::size_t)>
        &mutation) {
  auto manifest = read_bytes(directory / "manifest.v1.bin");
  mutation(manifest, offset);
  write_bytes(directory / "manifest.v1.bin", manifest);
  rewrite_marker(directory);
}

void test_transactional_failures(const MpiContext &context,
                                 const StructuredDecomposition &decomposition,
                                 const std::filesystem::path &root) {
  expect_transactional_read_failure(context, decomposition, root, 40,
                                    [](const std::filesystem::path &directory) {
                                      std::filesystem::remove(directory /
                                                              "COMPLETED");
                                    });
  expect_transactional_read_failure(
      context, decomposition, root, 41,
      [](const std::filesystem::path &directory) {
        auto marker = read_bytes(directory / "COMPLETED");
        marker[0] ^= std::byte{1};
        write_bytes(directory / "COMPLETED", marker);
      });
  expect_transactional_read_failure(context, decomposition, root, 42,
                                    [](const std::filesystem::path &directory) {
                                      std::filesystem::remove(
                                          directory / "manifest.v1.bin");
                                    });
  expect_transactional_read_failure(
      context, decomposition, root, 43,
      [](const std::filesystem::path &directory) {
        auto manifest = read_bytes(directory / "manifest.v1.bin");
        manifest[0] ^= std::byte{1};
        write_bytes(directory / "manifest.v1.bin", manifest);
        rewrite_marker(directory);
      });
  expect_transactional_read_failure(
      context, decomposition, root, 44,
      [](const std::filesystem::path &directory) {
        auto manifest = read_bytes(directory / "manifest.v1.bin");
        manifest.push_back(std::byte{0});
        write_bytes(directory / "manifest.v1.bin", manifest);
        rewrite_marker(directory);
      });
  expect_transactional_read_failure(
      context, decomposition, root, 45,
      [&](const std::filesystem::path &directory) {
        mutate_manifest_and_rebind(
            directory, 32U,
            [&](std::vector<std::byte> &bytes, std::size_t offset) {
              store_native<std::int32_t>(bytes, offset, context.size() + 1);
            });
      });
  expect_transactional_read_failure(
      context, decomposition, root, 46,
      [](const std::filesystem::path &directory) {
        mutate_manifest_and_rebind(
            directory, 51U,
            [](std::vector<std::byte> &bytes, std::size_t offset) {
              store_native<std::int32_t>(
                  bytes, offset, load_native<std::int32_t>(bytes, offset) + 1);
            });
      });
  expect_transactional_read_failure(
      context, decomposition, root, 47,
      [](const std::filesystem::path &directory) {
        mutate_manifest_and_rebind(
            directory, 83U,
            [](std::vector<std::byte> &bytes, std::size_t offset) {
              store_native<std::int32_t>(
                  bytes, offset, load_native<std::int32_t>(bytes, offset) + 1);
            });
      });
  expect_transactional_read_failure(
      context, decomposition, root, 48,
      [](const std::filesystem::path &directory) {
        mutate_manifest_and_rebind(
            directory, 67U,
            [](std::vector<std::byte> &bytes, std::size_t offset) {
              store_native<std::uint64_t>(
                  bytes, offset,
                  load_native<std::uint64_t>(bytes, offset) ^ UINT64_C(1));
            });
      });
  expect_transactional_read_failure(
      context, decomposition, root, 49,
      [&](const std::filesystem::path &directory) {
        std::filesystem::remove(directory / rank_filename(context.size() - 1));
      });
  expect_transactional_read_failure(
      context, decomposition, root, 50,
      [&](const std::filesystem::path &directory) {
        auto bytes = read_bytes(directory / rank_filename(context.size() - 1));
        bytes.back() ^= std::byte{1};
        write_bytes(directory / rank_filename(context.size() - 1), bytes);
      });
  expect_transactional_read_failure(
      context, decomposition, root, 51,
      [&](const std::filesystem::path &directory) {
        auto bytes = read_bytes(directory / rank_filename(context.size() - 1));
        bytes.push_back(std::byte{0});
        write_bytes(directory / rank_filename(context.size() - 1), bytes);
      });
}

void test_incompatible_destination_preflight(
    const MpiContext &context, const StructuredDecomposition &decomposition,
    const std::filesystem::path &root) {
  Fields source(decomposition.local_extent());
  fill_source(source, context.rank());
  const auto directory = step_directory(root, 52);
  hundun::runtime::write_restart_checkpoint(context, decomposition,
                                            source.registry, source.storage,
                                            directory, 52, 52.0 * kTimeScale);

  const DestinationMismatch mismatch =
      context.size() == 1
          ? DestinationMismatch::first_scalar_type
          : (context.rank() == 0 ? DestinationMismatch::none
                                 : DestinationMismatch::later_scalar_type);
  DestinationLayout destination(decomposition.local_extent(), mismatch);
  fill_destination_sentinels(destination);
  expect_same_collective_error(context, [&] {
    static_cast<void>(hundun::runtime::read_restart_checkpoint(
        context, decomposition, source.registry, destination.storage,
        directory));
  });
  verify_destination_sentinels(destination);
}

void run_full(const MpiContext &context,
              const StructuredDecomposition &decomposition,
              const std::filesystem::path &root) {
  test_success(context, decomposition, root);
  test_preflight_rejection(context, decomposition, root);
  test_writer_layout_preflight(context, decomposition, root);
  test_injected_failures(context, decomposition, root);
  test_transactional_failures(context, decomposition, root);
  test_incompatible_destination_preflight(context, decomposition, root);
}

void run_failure_repeat(const MpiContext &context,
                        const StructuredDecomposition &decomposition,
                        const std::filesystem::path &root) {
  test_injected_failures(context, decomposition, root);
  test_transactional_failures(context, decomposition, root);
}

int run_finalized(int argc, char **argv) {
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    auto context = MpiContext::duplicate(MPI_COMM_WORLD);
    auto decomposition = StructuredDecomposition::create(context, Int3{1, 1, 1},
                                                         {false, false, false});
    Fields fields(decomposition.local_extent());
    fill_source(fields, 0);
    const auto directory =
        std::filesystem::path{
            "/tmp/hundun-task10-finalized-" +
            std::to_string(static_cast<long long>(::getpid()))} /
        "step00000000";
    HUNDUN_CHECK(MPI_Finalize() == MPI_SUCCESS);
    result = hundun::test::run([&] {
      bool caught = false;
      try {
        hundun::runtime::write_restart_checkpoint(
            context, decomposition, fields.registry, fields.storage, directory,
            0, 0.0);
      } catch (const Error &error) {
        caught = std::strlen(error.what()) != 0U;
      }
      HUNDUN_CHECK(caught);
      HUNDUN_CHECK(!std::filesystem::exists(directory));
    });
  }
  return result;
}

} // namespace

int main(int argc, char **argv) {
  const std::string mode = argc > 1 ? argv[1] : "full";
  if (mode == "finalized") {
    return run_finalized(argc, argv);
  }

  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    auto context = MpiContext::duplicate(MPI_COMM_WORLD);
    const Int3 process_grid =
        context.size() == 1 ? Int3{1, 1, 1} : Int3{context.size(), 1, 1};
    auto decomposition = StructuredDecomposition::create(
        context, kGlobalExtent, kPeriodic, DecompositionOptions{process_grid});
    const auto root = broadcast_root(context);
    if (context.rank() == 0) {
      std::filesystem::remove_all(root);
    }
    context.barrier();
    result = hundun::test::run([&] {
      if (mode == "full" || mode == "evidence") {
        run_full(context, decomposition, root);
      } else if (mode == "failure") {
        run_failure_repeat(context, decomposition, root);
      } else {
        throw std::runtime_error("unknown restart checkpoint test mode");
      }
    });
    context.barrier();
    if (context.rank() == 0 && mode == "evidence") {
      std::cout << "EVIDENCE_ROOT=" << root.string() << '\n';
    }
    if (context.rank() == 0 && mode != "evidence") {
      std::filesystem::remove_all(root);
    }
    context.barrier();
  }
  return result;
}
