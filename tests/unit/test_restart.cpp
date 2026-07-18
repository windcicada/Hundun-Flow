// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/restart_binary.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "runtime/src/restart_detail.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::runtime::Box3;
using hundun::runtime::Error;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::OutputPolicy;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::detail::RestartRankMetadata;

constexpr Int3 kExtent{2, 2, 2};
constexpr Box3 kBox{Int3{3, 4, 5}, Int3{5, 6, 7}};
constexpr std::int64_t kStep = 17;
constexpr double kTime = 0.125;

template <class Function> void expect_error(Function &&function) {
  bool caught = false;
  try {
    function();
  } catch (const Error &error) {
    caught = std::strlen(error.what()) != 0U;
  }
  HUNDUN_CHECK(caught);
}

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

struct Fixture final {
  FieldRegistry registry;
  FieldId real{};
  FieldId integer{};
  FieldId transient{};
  FieldStorage storage;

  Fixture()
      : real(registry.declare_field(FieldDescriptor{
            u8"温度", "K", "restart-test", FunctionSpace::cell_average,
            ScalarType::float64, 2U, 1, true, RestartPolicy::persistent,
            OutputPolicy::selected})),
        integer(registry.declare_field(FieldDescriptor{
            "material", "1", "restart-test", FunctionSpace::cell_average,
            ScalarType::int32, 2U, 1, false, RestartPolicy::persistent,
            OutputPolicy::never})),
        transient(registry.declare_field(FieldDescriptor{
            "scratch", "1", "restart-test", FunctionSpace::cell_average,
            ScalarType::float64, 1U, 1, false, RestartPolicy::transient,
            OutputPolicy::never})),
        storage(freeze_and_return(), kExtent) {
    fill_source();
  }

private:
  FieldRegistry &freeze_and_return() {
    registry.freeze();
    return registry;
  }

  void fill_source() {
    auto real_view = storage.view<double>(real);
    auto int_view = storage.view<std::int32_t>(integer);
    auto transient_view = storage.view<double>(transient);
    for (int k = -1; k <= kExtent.z; ++k) {
      for (int j = -1; j <= kExtent.y; ++j) {
        for (int i = -1; i <= kExtent.x; ++i) {
          transient_view(i, j, k, 0) = -777.0;
          for (int component = 0; component < 2; ++component) {
            real_view(i, j, k, component) = -999.0;
            int_view(i, j, k, component) = -999;
          }
        }
      }
    }
    std::int32_t ordinal = 0;
    for (int k = 0; k < kExtent.z; ++k) {
      for (int j = 0; j < kExtent.y; ++j) {
        for (int i = 0; i < kExtent.x; ++i) {
          for (int component = 0; component < 2; ++component) {
            real_view(i, j, k, component) = 0.25 + static_cast<double>(ordinal);
            int_view(i, j, k, component) = 1000 + ordinal;
            ++ordinal;
          }
        }
      }
    }
  }
};

RestartRankMetadata metadata() {
  return RestartRankMetadata{0, 2, Int3{8, 9, 10}, kBox, kStep, kTime};
}

std::vector<std::byte> encoded_fixture(Fixture &fixture) {
  return hundun::runtime::detail::encode_restart_rank(
      metadata(), fixture.registry, fixture.storage);
}

void test_crc_standard_vector() {
  constexpr std::array<std::byte, 9> input{
      std::byte{'1'}, std::byte{'2'}, std::byte{'3'},
      std::byte{'4'}, std::byte{'5'}, std::byte{'6'},
      std::byte{'7'}, std::byte{'8'}, std::byte{'9'}};
  HUNDUN_CHECK(hundun::runtime::detail::crc64_ecma(
                   input.data(), input.size()) == UINT64_C(0x6C40DF5F0B497347));
}

void test_exact_bytes_and_round_trip() {
  Fixture source;
  const auto bytes = encoded_fixture(source);
  HUNDUN_CHECK(bytes.size() > 80U);
  const std::array<char, 8> magic{'H', 'U', 'N', 'D', 'U', 'N', 'R', '1'};
  HUNDUN_CHECK(std::memcmp(bytes.data(), magic.data(), magic.size()) == 0);
  HUNDUN_CHECK(load_native<std::uint32_t>(bytes, 8U) == 1U);
  HUNDUN_CHECK(load_native<std::uint32_t>(bytes, 12U) == UINT32_C(0x01020304));
  HUNDUN_CHECK(load_native<std::int32_t>(bytes, 16U) == 0);
  HUNDUN_CHECK(load_native<std::int32_t>(bytes, 20U) == 2);
  HUNDUN_CHECK(load_native<std::int64_t>(bytes, 60U) == kStep);
  HUNDUN_CHECK(load_native<double>(bytes, 68U) == kTime);
  HUNDUN_CHECK(load_native<std::uint32_t>(bytes, 76U) == 2U);

  std::size_t cursor = 80U;
  const std::uint32_t first_name_bytes =
      load_native<std::uint32_t>(bytes, cursor);
  cursor += sizeof(std::uint32_t);
  const std::string first_name(
      reinterpret_cast<const char *>(bytes.data() + cursor), first_name_bytes);
  HUNDUN_CHECK(first_name == u8"温度");
  cursor += first_name_bytes;
  HUNDUN_CHECK(load_native<std::uint32_t>(bytes, cursor) == 1U);
  cursor += sizeof(std::uint32_t);
  HUNDUN_CHECK(load_native<std::uint32_t>(bytes, cursor) == 2U);
  cursor += sizeof(std::uint32_t);
  HUNDUN_CHECK(load_native<std::uint64_t>(bytes, cursor) == 16U);
  cursor += sizeof(std::uint64_t);
  for (std::uint64_t index = 0; index < 16U; ++index) {
    HUNDUN_CHECK(load_native<double>(bytes, cursor) ==
                 0.25 + static_cast<double>(index));
    cursor += sizeof(double);
  }

  Fixture destination;
  auto real_view = destination.storage.view<double>(destination.real);
  auto int_view = destination.storage.view<std::int32_t>(destination.integer);
  auto transient_view = destination.storage.view<double>(destination.transient);
  for (int k = -1; k <= kExtent.z; ++k) {
    for (int j = -1; j <= kExtent.y; ++j) {
      for (int i = -1; i <= kExtent.x; ++i) {
        transient_view(i, j, k, 0) = 313.0;
        for (int component = 0; component < 2; ++component) {
          real_view(i, j, k, component) = 414.0;
          int_view(i, j, k, component) = 515;
        }
      }
    }
  }

  const auto staged = hundun::runtime::detail::decode_restart_rank(
      bytes, metadata(), destination.registry);
  const auto commit_plan = hundun::runtime::detail::make_restart_commit_plan(
      destination.registry, destination.storage, kExtent);
  hundun::runtime::detail::commit_restart_rank(staged, commit_plan);
  std::int32_t ordinal = 0;
  for (int k = 0; k < kExtent.z; ++k) {
    for (int j = 0; j < kExtent.y; ++j) {
      for (int i = 0; i < kExtent.x; ++i) {
        for (int component = 0; component < 2; ++component) {
          HUNDUN_CHECK(real_view(i, j, k, component) ==
                       0.25 + static_cast<double>(ordinal));
          HUNDUN_CHECK(int_view(i, j, k, component) == 1000 + ordinal);
          ++ordinal;
        }
      }
    }
  }
  for (const Int3 corner :
       {Int3{-1, -1, -1}, Int3{kExtent.x, kExtent.y, kExtent.z}}) {
    HUNDUN_CHECK(real_view(corner.x, corner.y, corner.z, 0) == 414.0);
    HUNDUN_CHECK(int_view(corner.x, corner.y, corner.z, 1) == 515);
    HUNDUN_CHECK(transient_view(corner.x, corner.y, corner.z, 0) == 313.0);
  }
  HUNDUN_CHECK(transient_view(0, 0, 0, 0) == 313.0);
}

void test_prefix_corruption() {
  Fixture fixture;
  const auto clean = encoded_fixture(fixture);
  const auto reject = [&](std::vector<std::byte> bytes) {
    expect_error([&] {
      static_cast<void>(hundun::runtime::detail::decode_restart_rank(
          bytes, metadata(), fixture.registry));
    });
  };

  auto bytes = clean;
  bytes[0] ^= std::byte{1};
  reject(bytes);
  bytes = clean;
  store_native<std::uint32_t>(bytes, 8U, 2U);
  reject(bytes);
  bytes = clean;
  store_native<std::uint32_t>(bytes, 12U, UINT32_C(0x04030201));
  reject(bytes);
  for (const std::size_t offset : {16U, 20U, 24U, 36U, 48U}) {
    bytes = clean;
    store_native<std::int32_t>(bytes, offset,
                               load_native<std::int32_t>(bytes, offset) + 1);
    reject(bytes);
  }
  bytes = clean;
  store_native<std::int32_t>(bytes, 36U,
                             std::numeric_limits<std::int32_t>::max());
  store_native<std::int32_t>(bytes, 48U, 1);
  reject(bytes);
  bytes = clean;
  store_native<std::int32_t>(bytes, 36U,
                             std::numeric_limits<std::int32_t>::min());
  store_native<std::int32_t>(bytes, 48U, -1);
  reject(bytes);
  bytes = clean;
  store_native<std::int64_t>(bytes, 60U, kStep + 1);
  reject(bytes);
  bytes = clean;
  store_native<double>(bytes, 68U, kTime + 1.0);
  reject(bytes);
}

void test_field_corruption_truncation_and_trailing() {
  Fixture fixture;
  const auto clean = encoded_fixture(fixture);
  const auto reject = [&](const std::vector<std::byte> &bytes) {
    expect_error([&] {
      static_cast<void>(hundun::runtime::detail::decode_restart_rank(
          bytes, metadata(), fixture.registry));
    });
  };

  const std::uint32_t name_bytes = load_native<std::uint32_t>(clean, 80U);
  const std::size_t name_offset = 84U;
  const std::size_t scalar_offset = name_offset + name_bytes;
  const std::size_t components_offset = scalar_offset + 4U;
  const std::size_t count_offset = components_offset + 4U;

  auto bytes = clean;
  bytes[name_offset] ^= std::byte{1};
  reject(bytes);
  bytes = clean;
  store_native<std::uint32_t>(bytes, scalar_offset, 2U);
  reject(bytes);
  bytes = clean;
  store_native<std::uint32_t>(bytes, scalar_offset, 99U);
  reject(bytes);
  bytes = clean;
  store_native<std::uint32_t>(bytes, components_offset, 3U);
  reject(bytes);
  bytes = clean;
  store_native<std::uint64_t>(bytes, count_offset, 15U);
  reject(bytes);
  bytes = clean;
  store_native<std::uint32_t>(bytes, 76U, 1U);
  reject(bytes);

  const std::array<std::size_t, 7> truncation_lengths{
      0U, 7U, 8U, 12U, 60U, 79U, clean.size() - 1U};
  for (const std::size_t length : truncation_lengths) {
    reject(std::vector<std::byte>(
        clean.begin(), clean.begin() + static_cast<std::ptrdiff_t>(length)));
  }
  bytes = clean;
  bytes.push_back(std::byte{0});
  reject(bytes);
}

void test_schema_and_checked_overflow() {
  Fixture fixture;
  const auto fingerprint =
      hundun::runtime::detail::restart_schema_fingerprint(fixture.registry);
  HUNDUN_CHECK(fingerprint != 0U);
  HUNDUN_CHECK(
      hundun::runtime::detail::checked_owned_value_count(kExtent, 2U) == 16U);
  expect_error([] {
    static_cast<void>(hundun::runtime::detail::checked_owned_value_count(
        Int3{std::numeric_limits<int>::max(), std::numeric_limits<int>::max(),
             std::numeric_limits<int>::max()},
        std::numeric_limits<std::uint32_t>::max()));
  });

  FieldRegistry registry;
  const auto id = registry.declare_field(
      FieldDescriptor{"mask", "1", "restart-test", FunctionSpace::cell_average,
                      ScalarType::uint8, 1U, 0, false,
                      RestartPolicy::persistent, OutputPolicy::never});
  static_cast<void>(id);
  registry.freeze();
  FieldStorage storage(registry, Int3{1, 1, 1});
  const RestartRankMetadata small{
      0, 1, Int3{1, 1, 1}, Box3{Int3{}, Int3{1, 1, 1}}, 0, 0.0};
  expect_error([&] {
    static_cast<void>(
        hundun::runtime::detail::encode_restart_rank(small, registry, storage));
  });
}

} // namespace

int main() {
  static_assert(std::is_same_v<decltype(hundun::runtime::RestartMetadata::step),
                               std::int64_t>);
  return hundun::test::run([] {
    test_crc_standard_vector();
    test_exact_bytes_and_round_trip();
    test_prefix_corruption();
    test_field_corruption_truncation_and_trailing();
    test_schema_and_checked_overflow();
  });
}
