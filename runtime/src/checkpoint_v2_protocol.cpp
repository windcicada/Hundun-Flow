// SPDX-License-Identifier: Apache-2.0

#include "checkpoint_v2_protocol.hpp"

#include "hundun/runtime/error.hpp"

#include <array>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <system_error>
#include <type_traits>

namespace hundun::runtime::checkpoint_v2 {
namespace {

constexpr std::array<std::uint8_t, 8> kRankMagic{
    'H', 'F', 'C', '2', 'R', 'N', 'K', '\0'};
constexpr std::array<std::uint8_t, 8> kManifestMagic{
    'H', 'F', 'C', '2', 'M', 'A', 'N', '\0'};
constexpr std::array<std::uint8_t, 8> kCompletedMagic{
    'H', 'F', 'C', '2', 'D', 'O', 'N', '\0'};
constexpr std::uint32_t kFormatVersion = 2U;
constexpr std::uint32_t kEndianMarker = UINT32_C(0x01020304);
constexpr std::size_t kMaximumStringBytes = 4096U;

[[noreturn]] void malformed(const char *message) { throw Error(message); }

std::size_t checked_size(std::uint64_t value) {
  if (value > static_cast<std::uint64_t>(
                  std::numeric_limits<std::size_t>::max()))
    malformed("Checkpoint v2 size exceeds this platform");
  return static_cast<std::size_t>(value);
}

} // namespace

std::uint64_t crc64_ecma(const void *data, std::size_t size) noexcept {
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  const auto *bytes = static_cast<const std::uint8_t *>(data);
  std::uint64_t crc = 0U;
  for (std::size_t index = 0; index < size; ++index) {
    crc ^= static_cast<std::uint64_t>(bytes[index]) << 56U;
    for (unsigned bit = 0; bit < 8U; ++bit)
      crc = (crc & (UINT64_C(1) << 63U)) != 0U
                ? (crc << 1U) ^ polynomial
                : crc << 1U;
  }
  return crc;
}

void Encoder::u8(std::uint8_t value) { bytes_.push_back(value); }
void Encoder::boolean(bool value) { u8(value ? 1U : 0U); }
void Encoder::u32(std::uint32_t value) {
  for (unsigned shift = 0U; shift < 32U; shift += 8U)
    u8(static_cast<std::uint8_t>(value >> shift));
}
void Encoder::i32(std::int32_t value) {
  std::uint32_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  u32(bits);
}
void Encoder::u64(std::uint64_t value) {
  for (unsigned shift = 0U; shift < 64U; shift += 8U)
    u8(static_cast<std::uint8_t>(value >> shift));
}
void Encoder::f64(double value) {
  static_assert(sizeof(value) == sizeof(std::uint64_t));
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  u64(bits);
}
void Encoder::string(const std::string &value) {
  if (value.size() > kMaximumStringBytes)
    malformed("Checkpoint v2 string exceeds 4096 bytes");
  if (value.find('\0') != std::string::npos)
    malformed("Checkpoint v2 string contains NUL");
  u32(static_cast<std::uint32_t>(value.size()));
  raw(value.data(), value.size());
}
void Encoder::raw(const void *data, std::size_t size) {
  if (size == 0U)
    return;
  if (data == nullptr)
    malformed("Checkpoint v2 encoder received null bytes");
  if (size > bytes_.max_size() - bytes_.size())
    malformed("Checkpoint v2 byte size overflows");
  const auto *first = static_cast<const std::uint8_t *>(data);
  bytes_.insert(bytes_.end(), first, first + size);
}
const std::vector<std::uint8_t> &Encoder::bytes() const noexcept {
  return bytes_;
}
std::vector<std::uint8_t> Encoder::take() && noexcept {
  return std::move(bytes_);
}

Decoder::Decoder(const std::vector<std::uint8_t> &bytes) : bytes_(&bytes) {}
std::vector<std::uint8_t> Decoder::raw(std::size_t count) {
  if (bytes_ == nullptr || count > bytes_->size() - offset_)
    malformed("Checkpoint v2 data is truncated");
  std::vector<std::uint8_t> result(bytes_->begin() +
                                       static_cast<std::ptrdiff_t>(offset_),
                                   bytes_->begin() + static_cast<std::ptrdiff_t>(
                                                         offset_ + count));
  offset_ += count;
  return result;
}
std::uint8_t Decoder::u8() {
  const auto result = raw(1U);
  return result.front();
}
bool Decoder::boolean() {
  const auto value = u8();
  if (value > 1U)
    malformed("Checkpoint v2 boolean is invalid");
  return value != 0U;
}
std::uint32_t Decoder::u32() {
  const auto value = raw(4U);
  std::uint32_t result{};
  for (unsigned index = 0U; index < 4U; ++index)
    result |= static_cast<std::uint32_t>(value[index]) << (8U * index);
  return result;
}
std::int32_t Decoder::i32() {
  const auto value = u32();
  std::int32_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}
std::uint64_t Decoder::u64() {
  const auto value = raw(8U);
  std::uint64_t result{};
  for (unsigned index = 0U; index < 8U; ++index)
    result |= static_cast<std::uint64_t>(value[index]) << (8U * index);
  return result;
}
double Decoder::f64() {
  const auto value = u64();
  double result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}
std::string Decoder::string() {
  const auto count = u32();
  if (count > kMaximumStringBytes)
    malformed("Checkpoint v2 string exceeds 4096 bytes");
  const auto value = raw(count);
  if (std::find(value.begin(), value.end(), std::uint8_t{}) != value.end())
    malformed("Checkpoint v2 string contains NUL");
  return std::string(value.begin(), value.end());
}
std::size_t Decoder::remaining() const noexcept {
  return bytes_ == nullptr ? 0U : bytes_->size() - offset_;
}
void Decoder::require_eof() const {
  if (remaining() != 0U)
    malformed("Checkpoint v2 data has trailing bytes");
}

std::vector<std::uint8_t>
encode_rank_wrapper(std::int32_t rank, std::int32_t rank_count,
                    const std::vector<std::uint8_t> &payload) {
  if (rank_count <= 0 || rank < 0 || rank >= rank_count)
    malformed("Checkpoint v2 rank wrapper identity is invalid");
  Encoder encoder;
  encoder.raw(kRankMagic.data(), kRankMagic.size());
  encoder.u32(kFormatVersion);
  encoder.u32(kEndianMarker);
  encoder.i32(rank);
  encoder.i32(rank_count);
  encoder.u64(static_cast<std::uint64_t>(payload.size()));
  encoder.raw(payload.data(), payload.size());
  return std::move(encoder).take();
}

RankWrapper decode_rank_wrapper(const std::vector<std::uint8_t> &bytes) {
  Decoder decoder(bytes);
  if (decoder.raw(kRankMagic.size()) !=
      std::vector<std::uint8_t>(kRankMagic.begin(), kRankMagic.end()))
    malformed("Checkpoint v2 rank magic is invalid");
  if (decoder.u32() != kFormatVersion)
    malformed("Checkpoint v2 rank version is invalid");
  if (decoder.u32() != kEndianMarker)
    malformed("Checkpoint v2 rank endian marker is invalid");
  RankWrapper result;
  result.rank = decoder.i32();
  result.rank_count = decoder.i32();
  if (result.rank_count <= 0 || result.rank < 0 ||
      result.rank >= result.rank_count)
    malformed("Checkpoint v2 rank identity is invalid");
  result.payload = decoder.raw(checked_size(decoder.u64()));
  decoder.require_eof();
  return result;
}

std::vector<std::uint8_t> encode_manifest(const Manifest &manifest) {
  if (manifest.rank_count == 0U || manifest.rank_count > 999999U ||
      manifest.ranks.size() != manifest.rank_count)
    malformed("Checkpoint v2 manifest rank count is invalid");
  Encoder encoder;
  encoder.raw(kManifestMagic.data(), kManifestMagic.size());
  encoder.u32(kFormatVersion);
  encoder.u32(kEndianMarker);
  encoder.u32(manifest.rank_count);
  encoder.i32(manifest.process_grid.x);
  encoder.i32(manifest.process_grid.y);
  encoder.i32(manifest.process_grid.z);
  for (const auto fingerprint : manifest.fingerprints)
    encoder.u64(fingerprint);
  encoder.u64(static_cast<std::uint64_t>(manifest.global_payload.size()));
  encoder.raw(manifest.global_payload.data(), manifest.global_payload.size());
  encoder.u32(static_cast<std::uint32_t>(manifest.ranks.size()));
  for (std::size_t index = 0; index < manifest.ranks.size(); ++index) {
    const auto &record = manifest.ranks[index];
    if (record.rank != static_cast<std::int32_t>(index) ||
        record.filename !=
            "rank-" +
                [&] {
                  std::string digits = std::to_string(index);
                  return std::string(6U - digits.size(), '0') + digits;
                }() +
                ".v2.bin")
      malformed("Checkpoint v2 manifest rank record is not canonical");
    encoder.i32(record.rank);
    encoder.i32(record.owned_box_begin.x);
    encoder.i32(record.owned_box_begin.y);
    encoder.i32(record.owned_box_begin.z);
    encoder.i32(record.owned_box_end.x);
    encoder.i32(record.owned_box_end.y);
    encoder.i32(record.owned_box_end.z);
    encoder.string(record.filename);
    encoder.u64(record.logical_byte_size);
    encoder.u64(record.actual_byte_size);
    encoder.u64(record.crc64);
    encoder.u64(record.local_layout_fingerprint);
  }
  return std::move(encoder).take();
}

Manifest decode_manifest(const std::vector<std::uint8_t> &bytes) {
  Decoder decoder(bytes);
  if (decoder.raw(kManifestMagic.size()) !=
      std::vector<std::uint8_t>(kManifestMagic.begin(), kManifestMagic.end()))
    malformed("Checkpoint v2 manifest magic is invalid");
  if (decoder.u32() != kFormatVersion || decoder.u32() != kEndianMarker)
    malformed("Checkpoint v2 manifest header is invalid");
  Manifest manifest;
  manifest.rank_count = decoder.u32();
  if (manifest.rank_count == 0U || manifest.rank_count > 999999U)
    malformed("Checkpoint v2 manifest rank count is invalid");
  manifest.process_grid = {decoder.i32(), decoder.i32(), decoder.i32()};
  for (auto &fingerprint : manifest.fingerprints)
    fingerprint = decoder.u64();
  manifest.global_payload = decoder.raw(checked_size(decoder.u64()));
  const auto record_count = decoder.u32();
  if (record_count != manifest.rank_count)
    malformed("Checkpoint v2 manifest record count is invalid");
  manifest.ranks.reserve(record_count);
  for (std::uint32_t index = 0; index < record_count; ++index) {
    ManifestRankRecord record;
    record.rank = decoder.i32();
    record.owned_box_begin = {decoder.i32(), decoder.i32(), decoder.i32()};
    record.owned_box_end = {decoder.i32(), decoder.i32(), decoder.i32()};
    record.filename = decoder.string();
    record.logical_byte_size = decoder.u64();
    record.actual_byte_size = decoder.u64();
    record.crc64 = decoder.u64();
    record.local_layout_fingerprint = decoder.u64();
    const std::string digits = std::to_string(index);
    const std::string expected =
        "rank-" + std::string(6U - digits.size(), '0') + digits + ".v2.bin";
    if (record.rank != static_cast<std::int32_t>(index) ||
        record.filename != expected)
      malformed("Checkpoint v2 manifest rank record is not canonical");
    manifest.ranks.push_back(std::move(record));
  }
  decoder.require_eof();
  return manifest;
}

std::vector<std::uint8_t>
encode_completed_marker(const CompletedMarker &marker) {
  Encoder encoder;
  encoder.raw(kCompletedMagic.data(), kCompletedMagic.size());
  encoder.u32(kFormatVersion);
  encoder.u32(kEndianMarker);
  encoder.u64(marker.manifest_actual_size);
  encoder.u64(marker.manifest_crc64);
  encoder.u64(marker.common_fingerprint);
  return std::move(encoder).take();
}

CompletedMarker
decode_completed_marker(const std::vector<std::uint8_t> &bytes) {
  Decoder decoder(bytes);
  if (decoder.raw(kCompletedMagic.size()) !=
      std::vector<std::uint8_t>(kCompletedMagic.begin(),
                                kCompletedMagic.end()))
    malformed("Checkpoint v2 completed magic is invalid");
  if (decoder.u32() != kFormatVersion || decoder.u32() != kEndianMarker)
    malformed("Checkpoint v2 completed header is invalid");
  CompletedMarker result;
  result.manifest_actual_size = decoder.u64();
  result.manifest_crc64 = decoder.u64();
  result.common_fingerprint = decoder.u64();
  decoder.require_eof();
  return result;
}

VerifiedFile write_verified_temporary(
    const std::filesystem::path &path,
    const std::vector<std::uint8_t> &bytes) {
  if (std::filesystem::exists(path))
    malformed("Checkpoint v2 temporary file already exists");
  {
    std::ofstream stream(path, std::ios::binary | std::ios::out);
    if (!stream)
      malformed("Checkpoint v2 temporary file could not be opened");
    if (!bytes.empty())
      stream.write(reinterpret_cast<const char *>(bytes.data()),
                   static_cast<std::streamsize>(bytes.size()));
    stream.flush();
    if (!stream)
      malformed("Checkpoint v2 temporary file write failed");
  }
  const auto reread =
      read_regular_file_exact(path, static_cast<std::uint64_t>(bytes.size()));
  if (reread != bytes)
    malformed("Checkpoint v2 temporary file verification failed");
  return {static_cast<std::uint64_t>(bytes.size()),
          crc64_ecma(bytes.data(), bytes.size())};
}

std::vector<std::uint8_t>
read_regular_file_exact(const std::filesystem::path &path,
                        std::uint64_t expected_size) {
  std::error_code error;
  const auto status = std::filesystem::symlink_status(path, error);
  if (error || status.type() != std::filesystem::file_type::regular)
    malformed("Checkpoint v2 entry is not a regular file");
  const auto actual = std::filesystem::file_size(path, error);
  if (error || actual != expected_size)
    malformed("Checkpoint v2 file size is invalid");
  std::vector<std::uint8_t> result(checked_size(actual));
  std::ifstream stream(path, std::ios::binary | std::ios::in);
  if (!stream)
    malformed("Checkpoint v2 file could not be opened");
  if (!result.empty())
    stream.read(reinterpret_cast<char *>(result.data()),
                static_cast<std::streamsize>(result.size()));
  if (!stream && !stream.eof())
    malformed("Checkpoint v2 file read failed");
  char trailing{};
  if (stream.get(trailing))
    malformed("Checkpoint v2 file has trailing bytes");
  return result;
}

void publish_no_overwrite(const std::filesystem::path &temporary,
                          const std::filesystem::path &final_path) {
  if (std::filesystem::exists(final_path) ||
      std::filesystem::is_symlink(final_path))
    malformed("Checkpoint v2 final file already exists");
  std::error_code error;
  std::filesystem::rename(temporary, final_path, error);
  if (error)
    malformed("Checkpoint v2 file publication failed");
}

} // namespace hundun::runtime::checkpoint_v2
