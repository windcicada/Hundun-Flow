// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::runtime {
class MpiContext;
}

namespace hundun::runtime::checkpoint_v2 {

enum class NumericFileFailure : std::uint8_t {
  integrity_unchecked,
  integrity,
  filesystem
};

class NumericFileError final : public Error {
public:
  NumericFileError(NumericFileFailure failure, const std::string &message)
      : Error(message), failure_(failure) {}
  NumericFileFailure failure() const noexcept { return failure_; }

private:
  NumericFileFailure failure_;
};

std::uint64_t crc64_ecma(const void *data, std::size_t size) noexcept;
std::size_t checked_size(std::uint64_t value);
std::size_t checked_product(std::size_t left, std::size_t right);
std::uint64_t checked_sum_u64(std::uint64_t left, std::uint64_t right);

class Encoder final {
public:
  void u8(std::uint8_t);
  void boolean(bool);
  void u32(std::uint32_t);
  void i32(std::int32_t);
  void u64(std::uint64_t);
  void f64(double);
  void string(const std::string &);
  void raw(const void *, std::size_t);
  const std::vector<std::uint8_t> &bytes() const noexcept;
  std::vector<std::uint8_t> take() && noexcept;

private:
  std::vector<std::uint8_t> bytes_;
};

class Decoder final {
public:
  explicit Decoder(const std::vector<std::uint8_t> &);
  std::uint8_t u8();
  bool boolean();
  std::uint32_t u32();
  std::int32_t i32();
  std::uint64_t u64();
  double f64();
  std::string string();
  std::vector<std::uint8_t> raw(std::size_t);
  std::size_t remaining() const noexcept;
  void require_eof() const;

private:
  const std::vector<std::uint8_t> *bytes_{};
  std::size_t offset_{};
};

struct RankWrapper final {
  std::int32_t rank{};
  std::int32_t rank_count{};
  std::vector<std::uint8_t> payload;
};

std::vector<std::uint8_t>
encode_rank_wrapper(std::int32_t rank, std::int32_t rank_count,
                    const std::vector<std::uint8_t> &payload);
RankWrapper decode_rank_wrapper(const std::vector<std::uint8_t> &bytes,
                                std::uint64_t expected_payload_size);

struct ManifestRankRecord final {
  std::int32_t rank{};
  Int3 owned_box_begin{};
  Int3 owned_box_end{};
  std::string filename;
  std::uint64_t logical_byte_size{};
  std::uint64_t actual_byte_size{};
  std::uint64_t crc64{};
  std::uint64_t local_layout_fingerprint{};
};

struct Manifest final {
  std::uint32_t rank_count{};
  Int3 process_grid{};
  std::array<std::uint64_t, 5> fingerprints{};
  std::vector<std::uint8_t> global_payload;
  std::vector<ManifestRankRecord> ranks;
};

std::vector<std::uint8_t> encode_manifest(const Manifest &);
Manifest decode_manifest(const std::vector<std::uint8_t> &,
                         std::uint32_t expected_rank_count,
                         std::uint64_t expected_global_payload_size);

struct CompletedMarker final {
  std::uint64_t manifest_actual_size{};
  std::uint64_t manifest_crc64{};
  std::uint64_t common_fingerprint{};
};

std::vector<std::uint8_t> encode_completed_marker(const CompletedMarker &);
CompletedMarker
decode_completed_marker(const std::vector<std::uint8_t> &bytes);

struct VerifiedFile final {
  std::uint64_t actual_size{};
  std::uint64_t crc64{};
};

VerifiedFile write_verified_temporary(
    const std::filesystem::path &path,
    const std::vector<std::uint8_t> &bytes);
void create_directory_exclusive(const std::filesystem::path &path);
std::vector<std::uint8_t>
read_regular_file_exact(const std::filesystem::path &path,
                        std::uint64_t expected_size);
void publish_no_overwrite(const std::filesystem::path &temporary,
                          const std::filesystem::path &final_path);

struct CollectiveResult final {
  bool ok{};
  int failing_rank{-1};
};

CollectiveResult converge_phase(const runtime::MpiContext &, bool local_ok,
                                std::uint64_t &collective_count,
                                std::string_view operation);
CollectiveResult opaque_bytes_agreement(const runtime::MpiContext &,
                                        const std::vector<std::uint8_t> &,
                                        std::uint64_t &collective_count,
                                        std::string_view operation);
bool opaque_bytes_agree(const runtime::MpiContext &,
                        const std::vector<std::uint8_t> &,
                        std::uint64_t &collective_count,
                        std::string_view operation);
std::vector<std::uint64_t>
allgather_u64(const runtime::MpiContext &, const std::uint64_t *local,
              std::size_t count, std::uint64_t &collective_count,
              std::string_view operation);
std::uint64_t allreduce_sum_u64(const runtime::MpiContext &,
                                std::uint64_t local,
                                std::uint64_t &collective_count,
                                std::string_view operation);
bool exact_directory_inventory(const std::filesystem::path &,
                               const std::vector<std::string> &);

} // namespace hundun::runtime::checkpoint_v2
