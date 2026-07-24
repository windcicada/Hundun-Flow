// SPDX-License-Identifier: Apache-2.0

#include "checkpoint_v2_protocol.hpp"
#include "flow/src/checkpoint_v2_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace {

template <class Function> bool rejects(Function &&function) {
  try {
    function();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

std::uint64_t bits(double value) {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void test_crc_and_little_endian_codec() {
  using namespace hundun::runtime::checkpoint_v2;
  const std::string check = "123456789";
  HUNDUN_CHECK(crc64_ecma(check.data(), check.size()) ==
               UINT64_C(0x6C40DF5F0B497347));

  Encoder encoder;
  encoder.u8(0x7fU);
  encoder.u32(UINT32_C(0x01020304));
  encoder.i32(-2);
  encoder.u64(UINT64_C(0x0102030405060708));
  encoder.f64(-0.0);
  encoder.string("ok");
  const std::vector<std::uint8_t> expected{
      0x7f, 0x04, 0x03, 0x02, 0x01, 0xfe, 0xff, 0xff, 0xff, 0x08, 0x07,
      0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x80, 0x02, 0x00, 0x00, 0x00, 'o',  'k'};
  HUNDUN_CHECK(encoder.bytes() == expected);

  Decoder decoder(encoder.bytes());
  HUNDUN_CHECK(decoder.u8() == 0x7fU);
  HUNDUN_CHECK(decoder.u32() == UINT32_C(0x01020304));
  HUNDUN_CHECK(decoder.i32() == -2);
  HUNDUN_CHECK(decoder.u64() == UINT64_C(0x0102030405060708));
  HUNDUN_CHECK(bits(decoder.f64()) == bits(-0.0));
  HUNDUN_CHECK(decoder.string() == "ok");
  decoder.require_eof();
}

void test_rank_wrapper_literal_and_corruption() {
  using namespace hundun::runtime::checkpoint_v2;
  const std::vector<std::uint8_t> payload{0xaa, 0x55};
  const auto encoded = encode_rank_wrapper(3, 8, payload);
  const std::vector<std::uint8_t> expected{
      'H',  'F',  'C',  '2',  'R',  'N',  'K',  0x00, 0x02, 0x00, 0x00, 0x00,
      0x04, 0x03, 0x02, 0x01, 0x03, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00,
      0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xaa, 0x55};
  HUNDUN_CHECK(encoded == expected);
  const auto decoded = decode_rank_wrapper(encoded, payload.size());
  HUNDUN_CHECK(decoded.rank == 3);
  HUNDUN_CHECK(decoded.rank_count == 8);
  HUNDUN_CHECK(decoded.payload == payload);
  const auto expected_crc = crc64_ecma(encoded.data(), encoded.size());
  const auto is_authenticated = [&](const std::vector<std::uint8_t> &candidate,
                                    std::uint64_t crc) {
    return hundun::flow::test::checkpoint_v2_authenticate_rank_wrapper_for_test(
        candidate, crc, encoded.size(), 3, 8, payload.size());
  };
  HUNDUN_CHECK(is_authenticated(encoded, expected_crc));
  for (std::size_t index = 0; index < encoded.size(); ++index) {
    auto corrupted = encoded;
    corrupted[index] ^= 1U;
    HUNDUN_CHECK(!is_authenticated(corrupted, expected_crc));
  }

  for (std::size_t removed = 1; removed <= encoded.size(); ++removed) {
    auto truncated = encoded;
    truncated.resize(encoded.size() - removed);
    HUNDUN_CHECK(!is_authenticated(truncated, expected_crc));
  }
  auto trailing = encoded;
  trailing.push_back(0U);
  HUNDUN_CHECK(!is_authenticated(trailing, expected_crc));
  auto bad_magic = encoded;
  bad_magic.front() ^= 1U;
  HUNDUN_CHECK(!is_authenticated(bad_magic, expected_crc));
  auto rebuilt_identity = encoded;
  rebuilt_identity[16U] = 4U;
  HUNDUN_CHECK(
      !is_authenticated(rebuilt_identity, crc64_ecma(rebuilt_identity.data(),
                                                     rebuilt_identity.size())));
}

void test_codec_limits() {
  using namespace hundun::runtime::checkpoint_v2;
  Encoder encoder;
  encoder.string(std::string(4096U, 'x'));
  Decoder maximum_string_decoder(encoder.bytes());
  HUNDUN_CHECK(maximum_string_decoder.string() == std::string(4096U, 'x'));
  maximum_string_decoder.require_eof();
  HUNDUN_CHECK(rejects([&] {
    Encoder nul;
    nul.string(std::string("a\0b", 3U));
  }));
  encoder = Encoder{};
  HUNDUN_CHECK(rejects([&] { encoder.string(std::string(4097U, 'x')); }));
  HUNDUN_CHECK(rejects([&] { encoder.string(std::string("\xc0\xaf", 2)); }));
  HUNDUN_CHECK(
      rejects([&] { encoder.string(std::string("\xed\xa0\x80", 3)); }));

  Encoder invalid_bool;
  invalid_bool.u8(2U);
  Decoder bool_decoder(invalid_bool.bytes());
  HUNDUN_CHECK(rejects([&] { bool_decoder.boolean(); }));

  Encoder valid;
  valid.boolean(true);
  Decoder valid_decoder(valid.bytes());
  HUNDUN_CHECK(valid_decoder.boolean());
  valid_decoder.require_eof();

  Encoder invalid_utf8_wire;
  invalid_utf8_wire.u32(2U);
  invalid_utf8_wire.u8(0xc0U);
  invalid_utf8_wire.u8(0xafU);
  Decoder invalid_utf8_decoder(invalid_utf8_wire.bytes());
  HUNDUN_CHECK(rejects([&] { invalid_utf8_decoder.string(); }));

  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(checked_size(
        static_cast<std::uint64_t>(std::numeric_limits<std::ptrdiff_t>::max()) +
        1U));
  }));
  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(
        checked_product(std::numeric_limits<std::size_t>::max(), 2U));
  }));
  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(
        checked_sum_u64(std::numeric_limits<std::uint64_t>::max(), 1U));
  }));
}

void test_integer_and_binary64_edges() {
  using namespace hundun::runtime::checkpoint_v2;
  constexpr std::array<std::int32_t, 5> signed_values{
      std::numeric_limits<std::int32_t>::min(), -1, 0, 1,
      std::numeric_limits<std::int32_t>::max()};
  constexpr std::array<std::uint64_t, 8> binary64_values{
      UINT64_C(0x0000000000000000), UINT64_C(0x8000000000000000),
      UINT64_C(0x0000000000000001), UINT64_C(0x0010000000000000),
      UINT64_C(0x7fefffffffffffff), UINT64_C(0x7ff0000000000000),
      UINT64_C(0xfff0000000000000), UINT64_C(0x7ff8000000001234)};
  Encoder encoder;
  for (const auto value : signed_values)
    encoder.i32(value);
  for (const auto value : binary64_values) {
    double item{};
    std::memcpy(&item, &value, sizeof(item));
    encoder.f64(item);
  }
  Decoder decoder(encoder.bytes());
  for (const auto value : signed_values)
    HUNDUN_CHECK(decoder.i32() == value);
  for (const auto value : binary64_values)
    HUNDUN_CHECK(bits(decoder.f64()) == value);
  decoder.require_eof();
}

void test_manifest_and_completed_marker() {
  using namespace hundun::runtime::checkpoint_v2;
  Manifest manifest;
  manifest.rank_count = 1U;
  manifest.process_grid = {1, 1, 1};
  manifest.fingerprints = {1U, 2U, 3U, 4U, 5U};
  manifest.global_payload = {0x12U, 0x34U};
  manifest.ranks.push_back(
      {0, {0, 0, 0}, {2, 3, 4}, "rank-000000.v2.bin", 192U, 226U, 7U, 8U});
  const auto bytes = encode_manifest(manifest);
  std::vector<std::uint8_t> expected;
  const auto u32 = [&](std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      expected.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  const auto i32 = [&](std::int32_t value) {
    std::uint32_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    u32(bits);
  };
  const auto u64 = [&](std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      expected.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  const std::string magic("HFC2MAN\0", 8);
  expected.insert(expected.end(), magic.begin(), magic.end());
  u32(2U);
  u32(UINT32_C(0x01020304));
  u32(1U);
  i32(1);
  i32(1);
  i32(1);
  for (std::uint64_t value = 1U; value <= 5U; ++value)
    u64(value);
  u64(2U);
  expected.push_back(0x12U);
  expected.push_back(0x34U);
  u32(1U);
  i32(0);
  for (const auto value : {0, 0, 0, 2, 3, 4})
    i32(value);
  const std::string filename = "rank-000000.v2.bin";
  u32(static_cast<std::uint32_t>(filename.size()));
  expected.insert(expected.end(), filename.begin(), filename.end());
  for (const auto value : {192U, 226U, 7U, 8U})
    u64(value);
  HUNDUN_CHECK(bytes == expected);
  const auto decoded = decode_manifest(bytes, manifest.rank_count,
                                       manifest.global_payload.size());
  HUNDUN_CHECK(decoded.rank_count == manifest.rank_count);
  HUNDUN_CHECK(decoded.process_grid.x == manifest.process_grid.x);
  HUNDUN_CHECK(decoded.process_grid.y == manifest.process_grid.y);
  HUNDUN_CHECK(decoded.process_grid.z == manifest.process_grid.z);
  HUNDUN_CHECK(decoded.fingerprints == manifest.fingerprints);
  HUNDUN_CHECK(decoded.global_payload == manifest.global_payload);
  HUNDUN_CHECK(decoded.ranks.size() == 1U);
  HUNDUN_CHECK(decoded.ranks.front().filename == "rank-000000.v2.bin");
  const auto expected_manifest_crc = crc64_ecma(bytes.data(), bytes.size());
  const auto manifest_is_authenticated =
      [&](const std::vector<std::uint8_t> &candidate, std::uint64_t crc) {
        return hundun::flow::test::checkpoint_v2_authenticate_manifest_for_test(
            candidate, crc, bytes.size(), manifest);
      };
  HUNDUN_CHECK(manifest_is_authenticated(bytes, expected_manifest_crc));
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    auto corrupted = bytes;
    corrupted[index] ^= 1U;
    HUNDUN_CHECK(!manifest_is_authenticated(corrupted, expected_manifest_crc));
  }
  for (std::size_t size = 0; size < bytes.size(); ++size) {
    auto truncated = bytes;
    truncated.resize(size);
    HUNDUN_CHECK(!manifest_is_authenticated(truncated, expected_manifest_crc));
  }
  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(decode_manifest(bytes, manifest.rank_count + 1U,
                                      manifest.global_payload.size()));
  }));
  HUNDUN_CHECK(rejects([&] {
    static_cast<void>(decode_manifest(bytes, manifest.rank_count,
                                      manifest.global_payload.size() + 1U));
  }));
  auto manifest_trailing = bytes;
  manifest_trailing.push_back(0U);
  HUNDUN_CHECK(
      !manifest_is_authenticated(manifest_trailing, expected_manifest_crc));
  auto rebuilt_manifest = bytes;
  rebuilt_manifest[24U] ^= 1U;
  HUNDUN_CHECK(!manifest_is_authenticated(
      rebuilt_manifest,
      crc64_ecma(rebuilt_manifest.data(), rebuilt_manifest.size())));

  CompletedMarker marker{static_cast<std::uint64_t>(bytes.size()),
                         crc64_ecma(bytes.data(), bytes.size()), 99U};
  const auto marker_bytes = encode_completed_marker(marker);
  std::vector<std::uint8_t> expected_marker{'H',  'F',  'C',  '2',  'D',  'O',
                                            'N',  0x00, 0x02, 0x00, 0x00, 0x00,
                                            0x04, 0x03, 0x02, 0x01};
  const auto marker_u64 = [&](std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      expected_marker.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  marker_u64(marker.manifest_actual_size);
  marker_u64(marker.manifest_crc64);
  marker_u64(marker.common_fingerprint);
  HUNDUN_CHECK(marker_bytes == expected_marker);
  const auto decoded_marker = decode_completed_marker(marker_bytes);
  HUNDUN_CHECK(decoded_marker.manifest_actual_size ==
               marker.manifest_actual_size);
  HUNDUN_CHECK(decoded_marker.manifest_crc64 == marker.manifest_crc64);
  HUNDUN_CHECK(decoded_marker.common_fingerprint == marker.common_fingerprint);
  const auto marker_is_authenticated =
      [&](const std::vector<std::uint8_t> &candidate) {
        return hundun::flow::test::
            checkpoint_v2_authenticate_completed_marker_for_test(
                candidate, marker.manifest_actual_size, marker.manifest_crc64,
                marker.common_fingerprint);
  };
  HUNDUN_CHECK(marker_is_authenticated(marker_bytes));
  for (std::size_t index = 0; index < marker_bytes.size(); ++index) {
    auto corrupted = marker_bytes;
    corrupted[index] ^= 1U;
    HUNDUN_CHECK(!marker_is_authenticated(corrupted));
  }
  for (std::size_t size = 0; size < marker_bytes.size(); ++size) {
    auto truncated = marker_bytes;
    truncated.resize(size);
    HUNDUN_CHECK(!marker_is_authenticated(truncated));
  }
  auto marker_trailing = marker_bytes;
  marker_trailing.push_back(0U);
  HUNDUN_CHECK(!marker_is_authenticated(marker_trailing));
  HUNDUN_CHECK(
      !hundun::flow::test::checkpoint_v2_authenticate_completed_marker_for_test(
          marker_bytes, marker.manifest_actual_size + 1U, marker.manifest_crc64,
          marker.common_fingerprint));
  HUNDUN_CHECK(
      !hundun::flow::test::checkpoint_v2_authenticate_completed_marker_for_test(
          marker_bytes, marker.manifest_actual_size, marker.manifest_crc64 + 1U,
          marker.common_fingerprint));
  HUNDUN_CHECK(
      !hundun::flow::test::checkpoint_v2_authenticate_completed_marker_for_test(
          marker_bytes, marker.manifest_actual_size, marker.manifest_crc64,
          marker.common_fingerprint + 1U));
}

void test_exact_read_phase_failures() {
  using namespace hundun::runtime::checkpoint_v2;
  const auto path =
      std::filesystem::temp_directory_path() / "hundun-task23-exact-read.bin";
  {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    const std::array<char, 4> bytes{'a', '\0', '\0', '\0'};
    stream.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  }
  test::set_exact_read_fault(test::ExactReadFault::truncate_after_size);
  try {
    static_cast<void>(read_regular_file_exact(path, 4U));
    HUNDUN_CHECK(false);
  } catch (const NumericFileError &error) {
    HUNDUN_CHECK(error.failure() == NumericFileFailure::integrity);
  }

  {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    stream.write("abcd", 4);
  }
  test::set_exact_read_fault(test::ExactReadFault::read_failure);
  try {
    static_cast<void>(read_regular_file_exact(path, 4U));
    HUNDUN_CHECK(false);
  } catch (const NumericFileError &error) {
    HUNDUN_CHECK(error.failure() == NumericFileFailure::filesystem);
  }

  test::set_exact_read_fault(test::ExactReadFault::close_failure);
  try {
    static_cast<void>(read_regular_file_exact(path, 4U));
    HUNDUN_CHECK(false);
  } catch (const NumericFileError &error) {
    HUNDUN_CHECK(error.failure() == NumericFileFailure::filesystem);
  }
  test::set_exact_read_fault(test::ExactReadFault::none);
  std::filesystem::remove(path);
}

void test_numeric_file_phase_failures() {
  using namespace hundun::runtime::checkpoint_v2;
  using Point = test::NumericFilePoint;
  const auto root = std::filesystem::temp_directory_path() /
                    "hundun-task23-numeric-file-points";
  std::filesystem::remove_all(root);
  const auto expect_filesystem = [&](auto &&operation) {
    try {
      operation();
      HUNDUN_CHECK(false);
    } catch (const NumericFileError &error) {
      HUNDUN_CHECK(error.failure() == NumericFileFailure::filesystem);
    }
  };

  for (const auto point : {Point::directory_status, Point::parent_status,
                           Point::directory_create}) {
    std::filesystem::remove_all(root);
    test::set_numeric_file_fault(point);
    expect_filesystem([&] { create_directory_exclusive(root); });
    HUNDUN_CHECK(!std::filesystem::exists(root));
  }

  std::filesystem::create_directories(root);
  const std::vector<std::uint8_t> payload{'a', 'b', 'c', 'd'};
  for (const auto point :
       {Point::write_status, Point::write_open, Point::write_body,
        Point::write_flush, Point::write_close, Point::read_status,
        Point::read_size, Point::read_open, Point::read_body,
        Point::read_close}) {
    const auto path = root / "temporary.bin";
    std::filesystem::remove(path);
    test::set_numeric_file_fault(point);
    expect_filesystem(
        [&] { static_cast<void>(write_verified_temporary(path, payload)); });
    std::filesystem::remove(path);
  }

  for (const auto point : {Point::publish_status, Point::publish_rename}) {
    const auto temporary = root / "publish.tmp";
    const auto final_path = root / "publish.bin";
    std::filesystem::remove(temporary);
    std::filesystem::remove(final_path);
    {
      std::ofstream stream(temporary, std::ios::binary);
      stream.write("x", 1);
    }
    test::set_numeric_file_fault(point);
    expect_filesystem([&] { publish_no_overwrite(temporary, final_path); });
    HUNDUN_CHECK(!std::filesystem::exists(final_path));
  }

  {
    const auto item = root / "item";
    std::ofstream stream(item, std::ios::binary);
    stream.write("x", 1);
    stream.close();
    test::set_numeric_file_fault(Point::inventory_iteration);
    expect_filesystem(
        [&] { static_cast<void>(exact_directory_inventory(root, {"item"})); });
  }
  test::set_numeric_file_fault(Point::none);
  std::filesystem::remove_all(root);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  const std::string mode(argv[1]);
  if (mode != "fast" && mode != "acceptance")
    return 2;
  return hundun::test::run([&] {
    test_crc_and_little_endian_codec();
    test_codec_limits();
    test_exact_read_phase_failures();
    if (mode == "acceptance") {
      test_rank_wrapper_literal_and_corruption();
      test_integer_and_binary64_edges();
      test_manifest_and_completed_marker();
      test_numeric_file_phase_failures();
    }
  });
}
