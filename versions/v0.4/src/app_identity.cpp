// SPDX-License-Identifier: Apache-2.0

#include "app_identity_detail.hpp"

#include "app_identity_build.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <string>
#include <string_view>

namespace hundun::v04::detail {
namespace {

constexpr std::uint32_t kRuntimeIdentityFailure = 10505U;
constexpr std::string_view kIdentitySchema =
    "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V1";

constexpr std::array<std::uint32_t, 64U> kSha256RoundConstants{{
    UINT32_C(0x428a2f98), UINT32_C(0x71374491), UINT32_C(0xb5c0fbcf),
    UINT32_C(0xe9b5dba5), UINT32_C(0x3956c25b), UINT32_C(0x59f111f1),
    UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5), UINT32_C(0xd807aa98),
    UINT32_C(0x12835b01), UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
    UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe), UINT32_C(0x9bdc06a7),
    UINT32_C(0xc19bf174), UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786),
    UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc), UINT32_C(0x2de92c6f),
    UINT32_C(0x4a7484aa), UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
    UINT32_C(0x983e5152), UINT32_C(0xa831c66d), UINT32_C(0xb00327c8),
    UINT32_C(0xbf597fc7), UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147),
    UINT32_C(0x06ca6351), UINT32_C(0x14292967), UINT32_C(0x27b70a85),
    UINT32_C(0x2e1b2138), UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
    UINT32_C(0x650a7354), UINT32_C(0x766a0abb), UINT32_C(0x81c2c92e),
    UINT32_C(0x92722c85), UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b),
    UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3), UINT32_C(0xd192e819),
    UINT32_C(0xd6990624), UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
    UINT32_C(0x19a4c116), UINT32_C(0x1e376c08), UINT32_C(0x2748774c),
    UINT32_C(0x34b0bcb5), UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a),
    UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3), UINT32_C(0x748f82ee),
    UINT32_C(0x78a5636f), UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
    UINT32_C(0x90befffa), UINT32_C(0xa4506ceb), UINT32_C(0xbef9a3f7),
    UINT32_C(0xc67178f2),
}};

constexpr std::uint32_t rotate_right(std::uint32_t value,
                                     unsigned shift) noexcept {
  return (value >> shift) | (value << (32U - shift));
}

class Sha256 final {
 public:
  bool update(const unsigned char* bytes, std::size_t size) noexcept {
    if ((bytes == nullptr && size != 0U) ||
        size > std::numeric_limits<std::uint64_t>::max() - total_bytes_)
      return false;
    total_bytes_ += static_cast<std::uint64_t>(size);
    while (size != 0U) {
      const std::size_t available = block_.size() - block_size_;
      const std::size_t copied = size < available ? size : available;
      std::memcpy(block_.data() + block_size_, bytes, copied);
      block_size_ += copied;
      bytes += copied;
      size -= copied;
      if (block_size_ == block_.size()) {
        transform(block_.data());
        block_size_ = 0U;
      }
    }
    return true;
  }

  bool update(std::string_view text) noexcept {
    return update(reinterpret_cast<const unsigned char*>(text.data()),
                  text.size());
  }

  std::array<unsigned char, 32U> finish() noexcept {
    const std::uint64_t bit_count = total_bytes_ * UINT64_C(8);
    block_[block_size_++] = 0x80U;
    if (block_size_ > 56U) {
      while (block_size_ < block_.size()) block_[block_size_++] = 0U;
      transform(block_.data());
      block_size_ = 0U;
    }
    while (block_size_ < 56U) block_[block_size_++] = 0U;
    for (unsigned index = 0U; index < 8U; ++index) {
      block_[63U - index] = static_cast<unsigned char>(
          bit_count >> static_cast<unsigned>(index * 8U));
    }
    transform(block_.data());
    block_size_ = 0U;

    std::array<unsigned char, 32U> digest{};
    for (std::size_t word = 0U; word < state_.size(); ++word) {
      digest[word * 4U] =
          static_cast<unsigned char>(state_[word] >> 24U);
      digest[word * 4U + 1U] =
          static_cast<unsigned char>(state_[word] >> 16U);
      digest[word * 4U + 2U] =
          static_cast<unsigned char>(state_[word] >> 8U);
      digest[word * 4U + 3U] = static_cast<unsigned char>(state_[word]);
    }
    return digest;
  }

 private:
  void transform(const unsigned char* block) noexcept {
    std::array<std::uint32_t, 64U> schedule{};
    for (std::size_t index = 0U; index < 16U; ++index) {
      const std::size_t offset = index * 4U;
      schedule[index] =
          (static_cast<std::uint32_t>(block[offset]) << 24U) |
          (static_cast<std::uint32_t>(block[offset + 1U]) << 16U) |
          (static_cast<std::uint32_t>(block[offset + 2U]) << 8U) |
          static_cast<std::uint32_t>(block[offset + 3U]);
    }
    for (std::size_t index = 16U; index < schedule.size(); ++index) {
      const std::uint32_t previous15 = schedule[index - 15U];
      const std::uint32_t previous2 = schedule[index - 2U];
      const std::uint32_t sigma0 = rotate_right(previous15, 7U) ^
                                   rotate_right(previous15, 18U) ^
                                   (previous15 >> 3U);
      const std::uint32_t sigma1 = rotate_right(previous2, 17U) ^
                                   rotate_right(previous2, 19U) ^
                                   (previous2 >> 10U);
      schedule[index] = schedule[index - 16U] + sigma0 +
                        schedule[index - 7U] + sigma1;
    }

    std::uint32_t a = state_[0U];
    std::uint32_t b = state_[1U];
    std::uint32_t c = state_[2U];
    std::uint32_t d = state_[3U];
    std::uint32_t e = state_[4U];
    std::uint32_t f = state_[5U];
    std::uint32_t g = state_[6U];
    std::uint32_t h = state_[7U];
    for (std::size_t index = 0U; index < schedule.size(); ++index) {
      const std::uint32_t sum1 = rotate_right(e, 6U) ^ rotate_right(e, 11U) ^
                                 rotate_right(e, 25U);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t temporary1 =
          h + sum1 + choice + kSha256RoundConstants[index] + schedule[index];
      const std::uint32_t sum0 = rotate_right(a, 2U) ^ rotate_right(a, 13U) ^
                                 rotate_right(a, 22U);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temporary2 = sum0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temporary1;
      d = c;
      c = b;
      b = a;
      a = temporary1 + temporary2;
    }
    state_[0U] += a;
    state_[1U] += b;
    state_[2U] += c;
    state_[3U] += d;
    state_[4U] += e;
    state_[5U] += f;
    state_[6U] += g;
    state_[7U] += h;
  }

  std::array<std::uint32_t, 8U> state_{{
      UINT32_C(0x6a09e667), UINT32_C(0xbb67ae85), UINT32_C(0x3c6ef372),
      UINT32_C(0xa54ff53a), UINT32_C(0x510e527f), UINT32_C(0x9b05688c),
      UINT32_C(0x1f83d9ab), UINT32_C(0x5be0cd19),
  }};
  std::array<unsigned char, 64U> block_{};
  std::size_t block_size_{};
  std::uint64_t total_bytes_{};
};

using DigestText =
    std::array<char, kRuntimeSha256HexCharacters + 1U>;
using GitObjectText =
    std::array<char, kRuntimeGitObjectHexCharacters + 1U>;

DigestText hex_digest(const std::array<unsigned char, 32U>& digest) noexcept {
  constexpr char kHex[] = "0123456789abcdef";
  DigestText text{};
  for (std::size_t index = 0U; index < digest.size(); ++index) {
    text[index * 2U] = kHex[digest[index] >> 4U];
    text[index * 2U + 1U] = kHex[digest[index] & 0x0fU];
  }
  return text;
}

bool valid_hex(std::string_view value, std::size_t size) noexcept {
  if (value.size() != size) return false;
  for (char character : value)
    if (!((character >= '0' && character <= '9') ||
          (character >= 'a' && character <= 'f')))
      return false;
  return true;
}

bool valid_digest(std::string_view digest) noexcept {
  return valid_hex(digest, kRuntimeSha256HexCharacters);
}

bool valid_digest(const DigestText& digest) noexcept {
  return digest[kRuntimeSha256HexCharacters] == '\0' &&
         valid_digest(std::string_view{digest.data(),
                                       kRuntimeSha256HexCharacters});
}

bool valid_git_object(const GitObjectText& object) noexcept {
  return object[kRuntimeGitObjectHexCharacters] == '\0' &&
         valid_hex(std::string_view{object.data(),
                                    kRuntimeGitObjectHexCharacters},
                   kRuntimeGitObjectHexCharacters);
}

bool copy_digest(std::string_view source, DigestText& destination) noexcept {
  if (!valid_digest(source)) return false;
  std::copy(source.begin(), source.end(), destination.begin());
  destination[kRuntimeSha256HexCharacters] = '\0';
  return true;
}

bool copy_git_object(std::string_view source,
                     GitObjectText& destination) noexcept {
  if (!valid_hex(source, kRuntimeGitObjectHexCharacters))
    return false;
  std::copy(source.begin(), source.end(), destination.begin());
  destination[kRuntimeGitObjectHexCharacters] = '\0';
  return true;
}

bool sha256_file(const char* path, DigestText& out) noexcept {
  const int descriptor = ::open(path, O_RDONLY | O_CLOEXEC);
  if (descriptor < 0) return false;
  struct stat metadata {};
  bool valid = ::fstat(descriptor, &metadata) == 0 &&
               S_ISREG(metadata.st_mode) && metadata.st_size > 0;
  Sha256 hash;
  std::array<unsigned char, 1024U * 1024U> block{};
  while (valid) {
    const ssize_t received = ::read(descriptor, block.data(), block.size());
    if (received > 0) {
      valid = hash.update(block.data(), static_cast<std::size_t>(received));
    } else if (received == 0) {
      break;
    } else if (errno != EINTR) {
      valid = false;
    }
  }
  if (::close(descriptor) != 0) valid = false;
  if (valid) out = hex_digest(hash.finish());
  return valid;
}

std::string identity_payload(const RuntimeCandidateIdentity& identity) {
  std::string payload;
  payload.reserve(420U);
  payload.append("schema=");
  payload.append(kIdentitySchema);
  payload.append("\nhead=");
  payload.append(identity.head.data(), kRuntimeGitObjectHexCharacters);
  payload.append("\ntree=");
  payload.append(identity.tree.data(), kRuntimeGitObjectHexCharacters);
  payload.append("\nbuild_manifest_sha256=");
  payload.append(identity.build_manifest.data(), kRuntimeSha256HexCharacters);
  payload.append("\nexecutable_sha256=");
  payload.append(identity.executable.data(), kRuntimeSha256HexCharacters);
  payload.push_back('\n');
  return payload;
}

DigestText identity_digest(const RuntimeCandidateIdentity& identity) {
  const std::string payload = identity_payload(identity);
  Sha256 hash;
  if (!hash.update(payload)) return {};
  return hex_digest(hash.finish());
}

std::array<unsigned char,
           2U * (kRuntimeGitObjectHexCharacters + 1U) +
               3U * (kRuntimeSha256HexCharacters + 1U)>
identity_bytes(const RuntimeCandidateIdentity& identity) noexcept {
  std::array<unsigned char,
             2U * (kRuntimeGitObjectHexCharacters + 1U) +
                 3U * (kRuntimeSha256HexCharacters + 1U)>
      bytes{};
  std::size_t offset = 0U;
  const auto append = [&](const auto& field) noexcept {
    for (char value : field)
      bytes[offset++] = static_cast<unsigned char>(value);
  };
  append(identity.head);
  append(identity.tree);
  append(identity.build_manifest);
  append(identity.executable);
  append(identity.identity);
  return bytes;
}

}  // namespace

bool valid_runtime_candidate_identity(
    const RuntimeCandidateIdentity& identity) noexcept try {
  if (!valid_git_object(identity.head) || !valid_git_object(identity.tree) ||
      !valid_digest(identity.build_manifest) ||
      !valid_digest(identity.executable) || !valid_digest(identity.identity))
    return false;
  return identity.identity == identity_digest(identity);
} catch (...) {
  return false;
}

PlanFingerprint runtime_sha256_fingerprint(
    const DigestText& digest) noexcept {
  if (!valid_digest(digest)) return 0U;
  PlanFingerprint value = 0U;
  for (std::size_t index = 0U; index < 16U; ++index) {
    const char encoded = digest[index];
    const unsigned nibble = encoded <= '9'
                                ? static_cast<unsigned>(encoded - '0')
                                : static_cast<unsigned>(encoded - 'a' + 10);
    value = (value << 4U) | static_cast<PlanFingerprint>(nibble);
  }
  return value;
}

bool valid_runtime_sha256(const RuntimeSha256Digest& digest) noexcept {
  return valid_digest(digest);
}

bool runtime_sha256_bytes(Span<const std::uint8_t> bytes,
                          RuntimeSha256Digest& out) noexcept {
  if (bytes.data == nullptr && bytes.size != 0U) return false;
  Sha256 hash;
  if (!hash.update(bytes.data, bytes.size)) return false;
  out = hex_digest(hash.finish());
  return valid_digest(out);
}

Status runtime_candidate_identity(MPI_Comm communicator,
                                  RuntimeCandidateIdentity& out) noexcept try {
  if (communicator == MPI_COMM_NULL)
    return {StatusCode::invalid_plan, kRuntimeIdentityFailure};
  RuntimeCandidateIdentity local;
  bool valid = copy_git_object(identity_source_commit, local.head) &&
               copy_git_object(identity_source_tree, local.tree) &&
               copy_digest(identity_build_manifest_sha256,
                           local.build_manifest) &&
               sha256_file("/proc/self/exe", local.executable);
  if (valid) {
    local.identity = identity_digest(local);
    valid = valid_runtime_candidate_identity(local) &&
            runtime_sha256_fingerprint(local.build_manifest) != 0U &&
            runtime_sha256_fingerprint(local.executable) != 0U;
  }
  const int local_valid = valid ? 1 : 0;
  int all_valid = 0;
  if (MPI_Allreduce(&local_valid, &all_valid, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRuntimeIdentityFailure};
  const auto local_bytes = identity_bytes(local);
  auto minimum = local_bytes;
  auto maximum = local_bytes;
  if (MPI_Allreduce(MPI_IN_PLACE, minimum.data(),
                    static_cast<int>(minimum.size()), MPI_UNSIGNED_CHAR,
                    MPI_MIN, communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, maximum.data(),
                    static_cast<int>(maximum.size()), MPI_UNSIGNED_CHAR,
                    MPI_MAX, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRuntimeIdentityFailure};
  if (all_valid == 0 || minimum != maximum)
    return {StatusCode::invalid_plan, kRuntimeIdentityFailure};
  out = local;
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kRuntimeIdentityFailure};
} catch (...) {
  return {StatusCode::invalid_plan, kRuntimeIdentityFailure};
}

}  // namespace hundun::v04::detail
