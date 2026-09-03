// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "rt_checkpoint_v2_protocol_detail.hpp"

#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_operation_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <system_error>
#include <type_traits>

namespace hundun::runtime::checkpoint_v2 {
namespace {

constexpr std::array<std::uint8_t, 8> kRankMagic{'H', 'F', 'C', '2',
                                                 'R', 'N', 'K', '\0'};
constexpr std::array<std::uint8_t, 8> kManifestMagic{'H', 'F', 'C', '2',
                                                     'M', 'A', 'N', '\0'};
constexpr std::array<std::uint8_t, 8> kCompletedMagic{'H', 'F', 'C', '2',
                                                      'D', 'O', 'N', '\0'};
constexpr std::uint32_t kFormatVersion = 2U;
constexpr std::uint32_t kEndianMarker = UINT32_C(0x01020304);
constexpr std::size_t kMaximumStringBytes = 4096U;

#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
test::ExactReadFault exact_read_fault{test::ExactReadFault::none};
std::uint32_t exact_read_fault_calls_before{};
test::PreparationPoint preparation_fault{test::PreparationPoint::none};
std::uint32_t preparation_fault_calls_before{};
test::NumericFilePoint numeric_file_fault{test::NumericFilePoint::none};
std::uint32_t numeric_file_fault_calls_before{};

test::ExactReadFault consume_exact_read_fault() noexcept {
  if (exact_read_fault_calls_before != 0U) {
    --exact_read_fault_calls_before;
    return test::ExactReadFault::none;
  }
  const auto result = exact_read_fault;
  exact_read_fault = test::ExactReadFault::none;
  return result;
}

bool consume_preparation_fault(test::PreparationPoint point) noexcept {
  if (preparation_fault != point)
    return false;
  if (preparation_fault_calls_before != 0U) {
    --preparation_fault_calls_before;
    return false;
  }
  preparation_fault = test::PreparationPoint::none;
  return true;
}

bool consume_numeric_file_fault(test::NumericFilePoint point) noexcept {
  if (numeric_file_fault != point)
    return false;
  if (numeric_file_fault_calls_before != 0U) {
    --numeric_file_fault_calls_before;
    return false;
  }
  numeric_file_fault = test::NumericFilePoint::none;
  return true;
}
#endif

[[noreturn]] void malformed(const char *message) { throw Error(message); }
[[noreturn]] void file_integrity(const char *message) {
  throw NumericFileError(NumericFileFailure::integrity, message);
}
[[noreturn]] void file_integrity_unchecked(const char *message) {
  throw NumericFileError(NumericFileFailure::integrity_unchecked, message);
}
[[noreturn]] void file_system(const char *message) {
  throw NumericFileError(NumericFileFailure::filesystem, message);
}

bool valid_utf8(const std::uint8_t *data, std::size_t size) noexcept {
  std::size_t index = 0U;
  while (index < size) {
    const auto first = data[index++];
    if (first <= 0x7fU)
      continue;
    std::uint32_t codepoint{};
    std::size_t continuation{};
    std::uint32_t minimum{};
    if (first >= 0xc2U && first <= 0xdfU) {
      codepoint = first & 0x1fU;
      continuation = 1U;
      minimum = 0x80U;
    } else if (first >= 0xe0U && first <= 0xefU) {
      codepoint = first & 0x0fU;
      continuation = 2U;
      minimum = 0x800U;
    } else if (first >= 0xf0U && first <= 0xf4U) {
      codepoint = first & 0x07U;
      continuation = 3U;
      minimum = 0x10000U;
    } else {
      return false;
    }
    if (continuation > size - index)
      return false;
    for (std::size_t offset = 0; offset < continuation; ++offset) {
      const auto next = data[index++];
      if ((next & 0xc0U) != 0x80U)
        return false;
      codepoint = (codepoint << 6U) | (next & 0x3fU);
    }
    if (codepoint < minimum || codepoint > 0x10ffffU ||
        (codepoint >= 0xd800U && codepoint <= 0xdfffU))
      return false;
  }
  return true;
}

} // namespace

std::size_t checked_size(std::uint64_t value) {
  if (value >
          static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
      value > static_cast<std::uint64_t>(
                  std::numeric_limits<std::ptrdiff_t>::max()))
    malformed("Checkpoint v2 size exceeds this platform");
  return static_cast<std::size_t>(value);
}

std::size_t checked_product(std::size_t left, std::size_t right) {
  if (right != 0U && left > std::numeric_limits<std::size_t>::max() / right)
    malformed("Checkpoint v2 size product overflows");
  return left * right;
}

std::uint64_t checked_sum_u64(std::uint64_t left, std::uint64_t right) {
  if (left > std::numeric_limits<std::uint64_t>::max() - right)
    malformed("Checkpoint v2 byte sum overflows");
  return left + right;
}

std::uint64_t crc64_ecma(const void *data, std::size_t size) noexcept {
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  const auto *bytes = static_cast<const std::uint8_t *>(data);
  std::uint64_t crc = 0U;
  for (std::size_t index = 0; index < size; ++index) {
    crc ^= static_cast<std::uint64_t>(bytes[index]) << 56U;
    for (unsigned bit = 0; bit < 8U; ++bit)
      crc = (crc & (UINT64_C(1) << 63U)) != 0U ? (crc << 1U) ^ polynomial
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
  if (!valid_utf8(reinterpret_cast<const std::uint8_t *>(value.data()),
                  value.size()))
    malformed("Checkpoint v2 string is not valid UTF-8");
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
  std::vector<std::uint8_t> result(
      bytes_->begin() + static_cast<std::ptrdiff_t>(offset_),
      bytes_->begin() + static_cast<std::ptrdiff_t>(offset_ + count));
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
  if (!valid_utf8(value.data(), value.size()))
    malformed("Checkpoint v2 string is not valid UTF-8");
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

RankWrapper decode_rank_wrapper(const std::vector<std::uint8_t> &bytes,
                                std::uint64_t expected_payload_size) {
  const auto expected_total = checked_sum_u64(32U, expected_payload_size);
  if (bytes.size() != checked_size(expected_total))
    malformed("Checkpoint v2 rank wrapper size is invalid");
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
  const auto declared_size = decoder.u64();
  if (declared_size != expected_payload_size)
    malformed("Checkpoint v2 rank payload size is invalid");
  result.payload = decoder.raw(checked_size(declared_size));
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
        record.filename != "rank-" + [&] {
                  std::string digits = std::to_string(index);
                  return std::string(6U - digits.size(), '0') + digits;
        }() + ".v2.bin")
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

Manifest decode_manifest(const std::vector<std::uint8_t> &bytes,
                         std::uint32_t expected_rank_count,
                         std::uint64_t expected_global_payload_size) {
  Decoder decoder(bytes);
  if (decoder.raw(kManifestMagic.size()) !=
      std::vector<std::uint8_t>(kManifestMagic.begin(), kManifestMagic.end()))
    malformed("Checkpoint v2 manifest magic is invalid");
  if (decoder.u32() != kFormatVersion || decoder.u32() != kEndianMarker)
    malformed("Checkpoint v2 manifest header is invalid");
  Manifest manifest;
  manifest.rank_count = decoder.u32();
  if (manifest.rank_count == 0U || manifest.rank_count > 999999U ||
      manifest.rank_count != expected_rank_count)
    malformed("Checkpoint v2 manifest rank count is invalid");
  manifest.process_grid = {decoder.i32(), decoder.i32(), decoder.i32()};
  for (auto &fingerprint : manifest.fingerprints)
    fingerprint = decoder.u64();
  const auto declared_global_size = decoder.u64();
  if (declared_global_size != expected_global_payload_size)
    malformed("Checkpoint v2 global payload size is invalid");
  manifest.global_payload = decoder.raw(checked_size(declared_global_size));
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
      std::vector<std::uint8_t>(kCompletedMagic.begin(), kCompletedMagic.end()))
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

VerifiedFile write_verified_temporary(const std::filesystem::path &path,
    const std::vector<std::uint8_t> &bytes) {
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::write_status))
    file_system("Checkpoint v2 injected temporary status failure");
#endif
  std::error_code status_error;
  const auto status = std::filesystem::symlink_status(path, status_error);
  if (status_error && status_error != std::errc::no_such_file_or_directory)
    file_system("Checkpoint v2 temporary file status failed");
  if (!status_error && status.type() != std::filesystem::file_type::not_found)
    file_system("Checkpoint v2 temporary file already exists");
  if (bytes.size() >
      static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
    malformed("Checkpoint v2 file exceeds stream size");
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::write_open))
    file_system("Checkpoint v2 injected temporary open failure");
#endif
  std::ofstream stream(path, std::ios::binary | std::ios::out);
  if (!stream)
    file_system("Checkpoint v2 temporary file could not be opened");
  if (!bytes.empty())
    stream.write(reinterpret_cast<const char *>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::write_body))
    file_system("Checkpoint v2 injected temporary write failure");
#endif
  if (!stream)
    file_system("Checkpoint v2 temporary file write failed");
  stream.flush();
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::write_flush))
    file_system("Checkpoint v2 injected temporary flush failure");
#endif
  if (!stream)
    file_system("Checkpoint v2 temporary file flush failed");
  stream.close();
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::write_close))
    file_system("Checkpoint v2 injected temporary close failure");
#endif
  if (stream.fail())
    file_system("Checkpoint v2 temporary file close failed");
  const auto reread =
      read_regular_file_exact(path, static_cast<std::uint64_t>(bytes.size()));
  if (reread != bytes)
    file_integrity("Checkpoint v2 temporary file verification failed");
  return {static_cast<std::uint64_t>(bytes.size()),
          crc64_ecma(bytes.data(), bytes.size())};
}

void create_directory_exclusive(const std::filesystem::path &path) {
  std::error_code error;
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::directory_status))
    file_system("Checkpoint v2 injected directory status failure");
#endif
  const auto target_status = std::filesystem::symlink_status(path, error);
  const bool absent =
      error == std::errc::no_such_file_or_directory ||
      target_status.type() == std::filesystem::file_type::not_found;
  if (!absent || std::filesystem::is_symlink(target_status))
    file_system("Checkpoint v2 target already exists");

  error.clear();
  auto parent = path.parent_path();
  if (parent.empty())
    parent = ".";
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::parent_status))
    file_system("Checkpoint v2 injected parent status failure");
#endif
  const auto parent_status = std::filesystem::symlink_status(parent, error);
  if (error || parent_status.type() != std::filesystem::file_type::directory)
    file_system("Checkpoint v2 parent is not a directory");
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::directory_create))
    file_system("Checkpoint v2 injected directory create failure");
#endif
  if (!std::filesystem::create_directory(path, error) || error)
    file_system("Checkpoint v2 directory creation failed");
}

std::vector<std::uint8_t>
read_regular_file_exact(const std::filesystem::path &path,
                        std::uint64_t expected_size) {
  std::error_code error;
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::read_status))
    file_system("Checkpoint v2 injected read status failure");
#endif
  const auto status = std::filesystem::symlink_status(path, error);
  if (error) {
    if (error == std::errc::no_such_file_or_directory)
      file_integrity_unchecked("Checkpoint v2 required file is missing");
    file_system("Checkpoint v2 file status failed");
  }
  if (status.type() != std::filesystem::file_type::regular)
    file_integrity_unchecked("Checkpoint v2 entry is not a regular file");
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::read_size))
    file_system("Checkpoint v2 injected read size failure");
#endif
  const auto actual = std::filesystem::file_size(path, error);
  if (error) {
    if (error == std::errc::no_such_file_or_directory)
      file_integrity_unchecked("Checkpoint v2 required file is missing");
    file_system("Checkpoint v2 file size query failed");
  }
  if (actual != expected_size)
    file_integrity("Checkpoint v2 file size is invalid");
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  const auto fault = consume_exact_read_fault();
  if (fault == test::ExactReadFault::truncate_after_size) {
    std::filesystem::resize_file(path, actual == 0U ? 0U : actual - 1U, error);
    if (error)
      file_system("Checkpoint v2 test truncation failed");
  }
#endif
  const auto allocation_size = checked_size(actual);
  if (allocation_size >
      static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
    malformed("Checkpoint v2 file exceeds stream size");
  std::vector<std::uint8_t> result(allocation_size);
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::read_open))
    file_system("Checkpoint v2 injected read open failure");
#endif
  std::ifstream stream(path, std::ios::binary | std::ios::in);
  if (!stream)
    file_system("Checkpoint v2 file could not be opened");
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (fault == test::ExactReadFault::read_failure)
    file_system("Checkpoint v2 injected read failure");
#endif
  if (!result.empty())
    stream.read(reinterpret_cast<char *>(result.data()),
                static_cast<std::streamsize>(result.size()));
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::read_body))
    file_system("Checkpoint v2 injected read body failure");
#endif
  if (stream.gcount() != static_cast<std::streamsize>(result.size()))
    file_integrity("Checkpoint v2 file became shorter while reading");
  if (!stream)
    file_system("Checkpoint v2 file read failed");
  char trailing{};
  if (stream.get(trailing))
    file_integrity("Checkpoint v2 file has trailing bytes");
  if (!stream.eof() && stream.fail())
    file_system("Checkpoint v2 file read failed");
  stream.clear();
  stream.close();
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::read_close))
    file_system("Checkpoint v2 injected read close failure");
#endif
  if (stream.fail())
    file_system("Checkpoint v2 file close failed");
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (fault == test::ExactReadFault::close_failure)
    file_system("Checkpoint v2 injected close failure");
#endif
  return result;
}

void publish_no_overwrite(const std::filesystem::path &temporary,
                          const std::filesystem::path &final_path) {
  std::error_code error;
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::publish_status))
    file_system("Checkpoint v2 injected publish status failure");
#endif
  const auto status = std::filesystem::symlink_status(final_path, error);
  if (error && error != std::errc::no_such_file_or_directory)
    file_system("Checkpoint v2 final file status failed");
  if (!error && status.type() != std::filesystem::file_type::not_found)
    file_system("Checkpoint v2 final file already exists");
  error.clear();
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::publish_rename))
    file_system("Checkpoint v2 injected publish rename failure");
#endif
  std::filesystem::rename(temporary, final_path, error);
  if (error)
    file_system("Checkpoint v2 file publication failed");
}

CollectiveResult converge_phase(const runtime::MpiContext &mpi, bool local_ok,
                                std::uint64_t &collective_count,
                                std::string_view operation) {
  const int candidate = local_ok ? mpi.size() : mpi.rank();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      operation);
  ++collective_count;
  return {lowest == mpi.size(), lowest == mpi.size() ? -1 : lowest};
}

bool opaque_bytes_agree(const runtime::MpiContext &mpi,
                        const std::vector<std::uint8_t> &bytes,
                        std::uint64_t &collective_count,
                        std::string_view operation) {
  return opaque_bytes_agreement(mpi, bytes, collective_count, operation).ok;
}

CollectiveResult opaque_bytes_agreement(const runtime::MpiContext &mpi,
                                        const std::vector<std::uint8_t> &bytes,
                                        std::uint64_t &collective_count,
                                        std::string_view operation) {
  constexpr std::uint64_t maximum = UINT64_C(1024) * UINT64_C(1024);
  std::uint64_t root_size =
      mpi.rank() == 0 ? static_cast<std::uint64_t>(bytes.size()) : 0U;
  runtime::check_mpi_result(
      MPI_Bcast(&root_size, 1, MPI_UINT64_T, 0, mpi.comm()), operation);
  ++collective_count;
  const bool size_valid =
      root_size <= maximum && bytes.size() <= checked_size(maximum);
  std::vector<std::uint8_t> root;
  bool prepared = true;
  try {
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(test::PreparationPoint::opaque_bytes_buffer))
      throw std::bad_alloc();
#endif
    root.resize(root_size <= maximum ? checked_size(root_size) : 0U);
    if (mpi.rank() == 0 && size_valid)
      root = bytes;
  } catch (...) {
    prepared = false;
  }
  const auto preparation =
      converge_phase(mpi, prepared, collective_count,
                     "MPI_Allreduce(Checkpoint opaque buffer preparation)");
  if (!preparation.ok)
    throw CollectivePreparationError(
        preparation.failing_rank,
        "Checkpoint v2 opaque buffer preparation failed");
  runtime::check_mpi_result(MPI_Bcast(root.data(),
                                      static_cast<int>(root.size()), MPI_BYTE,
                                      0, mpi.comm()),
      operation);
  ++collective_count;
  return converge_phase(mpi, size_valid && bytes == root, collective_count,
                        operation);
}

std::vector<std::uint64_t> allgather_u64(const runtime::MpiContext &mpi,
                                         const std::uint64_t *local,
                                         std::size_t count,
                                         std::uint64_t &collective_count,
              std::string_view operation) {
  if ((count != 0U && local == nullptr) ||
      count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    malformed("Checkpoint v2 opaque gather count is invalid");
  const auto total = checked_product(
      checked_size(static_cast<std::uint64_t>(mpi.size())), count);
  std::vector<std::uint64_t> result;
  bool prepared = true;
  try {
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(test::PreparationPoint::allgather_result))
      throw std::bad_alloc();
#endif
    result.resize(total);
  } catch (...) {
    prepared = false;
  }
  const auto preparation =
      converge_phase(mpi, prepared, collective_count,
                     "MPI_Allreduce(Checkpoint gather workspace preparation)");
  if (!preparation.ok)
    throw CollectivePreparationError(
        preparation.failing_rank,
        "Checkpoint v2 gather workspace preparation failed");
  runtime::check_mpi_result(
      MPI_Allgather(local, static_cast<int>(count), MPI_UINT64_T, result.data(),
                    static_cast<int>(count), MPI_UINT64_T, mpi.comm()),
      operation);
  ++collective_count;
  return result;
}

std::uint64_t allreduce_sum_u64(const runtime::MpiContext &mpi,
                                std::uint64_t local,
                                std::uint64_t &collective_count,
                                std::string_view operation) {
  std::vector<std::uint64_t> items;
  bool prepared = true;
  try {
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
    if (consume_preparation_fault(test::PreparationPoint::allreduce_workspace))
      throw std::bad_alloc();
#endif
    items.resize(checked_size(static_cast<std::uint64_t>(mpi.size())));
  } catch (...) {
    prepared = false;
  }
  const auto preparation = converge_phase(
      mpi, prepared, collective_count,
      "MPI_Allreduce(Checkpoint reduction workspace preparation)");
  if (!preparation.ok)
    throw CollectivePreparationError(
        preparation.failing_rank,
        "Checkpoint v2 reduction workspace preparation failed");
  runtime::check_mpi_result(MPI_Allgather(&local, 1, MPI_UINT64_T, items.data(),
                                          1, MPI_UINT64_T, mpi.comm()),
      operation);
  ++collective_count;
  std::uint64_t result{};
  for (const auto item : items)
    result = checked_sum_u64(result, item);
  return result;
}

bool exact_directory_inventory(const std::filesystem::path &directory,
                               const std::vector<std::string> &names) {
  std::error_code error;
  const auto status = std::filesystem::symlink_status(directory, error);
  if (error) {
    if (error == std::errc::no_such_file_or_directory)
      return false;
    file_system("Checkpoint v2 directory status failed");
  }
  if (status.type() != std::filesystem::file_type::directory)
    return false;
  const std::set<std::string> expected(names.begin(), names.end());
  if (expected.size() != names.size())
    return false;
  std::set<std::string> observed;
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  if (consume_numeric_file_fault(test::NumericFilePoint::inventory_iteration))
    file_system("Checkpoint v2 injected directory iteration failure");
#endif
  for (std::filesystem::directory_iterator iterator(directory, error), end;
       !error && iterator != end; iterator.increment(error)) {
    const auto entry_status = iterator->symlink_status(error);
    if (error)
      file_system("Checkpoint v2 directory entry status failed");
    if (entry_status.type() != std::filesystem::file_type::regular)
      return false;
    observed.insert(iterator->path().filename().generic_string());
  }
  if (error)
    file_system("Checkpoint v2 directory iteration failed");
  return observed == expected;
}

#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
namespace test {

void set_exact_read_fault(ExactReadFault fault,
                          std::uint32_t calls_before) noexcept {
  exact_read_fault = fault;
  exact_read_fault_calls_before = calls_before;
}

void set_preparation_fault(PreparationPoint point,
                           std::uint32_t calls_before) noexcept {
  preparation_fault = point;
  preparation_fault_calls_before = calls_before;
}

void set_numeric_file_fault(NumericFilePoint point,
                            std::uint32_t calls_before) noexcept {
  numeric_file_fault = point;
  numeric_file_fault_calls_before = calls_before;
}

} // namespace test
#endif

} // namespace hundun::runtime::checkpoint_v2
