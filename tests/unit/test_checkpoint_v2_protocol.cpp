// SPDX-License-Identifier: Apache-2.0

#include "checkpoint_v2_protocol.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

template <class Function>
bool rejects(Function &&function) {
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
      0x7f, 0x04, 0x03, 0x02, 0x01, 0xfe, 0xff, 0xff, 0xff,
      0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80,
      0x02, 0x00, 0x00, 0x00, 'o',  'k'};
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
  HUNDUN_CHECK(encoded.size() == 34U);
  HUNDUN_CHECK(std::string(encoded.begin(), encoded.begin() + 8) ==
               std::string("HFC2RNK\0", 8));
  const auto decoded = decode_rank_wrapper(encoded);
  HUNDUN_CHECK(decoded.rank == 3);
  HUNDUN_CHECK(decoded.rank_count == 8);
  HUNDUN_CHECK(decoded.payload == payload);

  for (std::size_t removed = 1; removed <= encoded.size(); ++removed) {
    auto truncated = encoded;
    truncated.resize(encoded.size() - removed);
    HUNDUN_CHECK(rejects([&] { decode_rank_wrapper(truncated); }));
  }
  auto trailing = encoded;
  trailing.push_back(0U);
  HUNDUN_CHECK(rejects([&] { decode_rank_wrapper(trailing); }));
  auto bad_magic = encoded;
  bad_magic.front() ^= 1U;
  HUNDUN_CHECK(rejects([&] { decode_rank_wrapper(bad_magic); }));
}

void test_codec_limits() {
  using namespace hundun::runtime::checkpoint_v2;
  Encoder encoder;
  HUNDUN_CHECK(
      rejects([&] { encoder.string(std::string(4097U, 'x')); }));

  Encoder invalid_bool;
  invalid_bool.u8(2U);
  Decoder bool_decoder(invalid_bool.bytes());
  HUNDUN_CHECK(rejects([&] { bool_decoder.boolean(); }));

  Encoder valid;
  valid.boolean(true);
  Decoder valid_decoder(valid.bytes());
  HUNDUN_CHECK(valid_decoder.boolean());
  valid_decoder.require_eof();
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
  const auto decoded = decode_manifest(bytes);
  HUNDUN_CHECK(decoded.rank_count == manifest.rank_count);
  HUNDUN_CHECK(decoded.process_grid.x == manifest.process_grid.x);
  HUNDUN_CHECK(decoded.process_grid.y == manifest.process_grid.y);
  HUNDUN_CHECK(decoded.process_grid.z == manifest.process_grid.z);
  HUNDUN_CHECK(decoded.fingerprints == manifest.fingerprints);
  HUNDUN_CHECK(decoded.global_payload == manifest.global_payload);
  HUNDUN_CHECK(decoded.ranks.size() == 1U);
  HUNDUN_CHECK(decoded.ranks.front().filename == "rank-000000.v2.bin");

  CompletedMarker marker{static_cast<std::uint64_t>(bytes.size()),
                         crc64_ecma(bytes.data(), bytes.size()), 99U};
  const auto marker_bytes = encode_completed_marker(marker);
  HUNDUN_CHECK(marker_bytes.size() == 40U);
  const auto decoded_marker = decode_completed_marker(marker_bytes);
  HUNDUN_CHECK(decoded_marker.manifest_actual_size ==
               marker.manifest_actual_size);
  HUNDUN_CHECK(decoded_marker.manifest_crc64 == marker.manifest_crc64);
  HUNDUN_CHECK(decoded_marker.common_fingerprint ==
               marker.common_fingerprint);
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_crc_and_little_endian_codec();
    test_rank_wrapper_literal_and_corruption();
    test_codec_limits();
    test_manifest_and_completed_marker();
  });
}
