// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_io.hpp"

#include <mpi.h>

namespace hundun::v04::detail {

inline constexpr std::string_view kRuntimeEvidenceSchema =
    "HUNDUN_V04_EVIDENCE_V8";
inline constexpr std::string_view kRuntimeCandidateIdentitySchema =
    "HUNDUN_V04_RUNTIME_CANDIDATE_IDENTITY_V2";

Status runtime_candidate_identity(MPI_Comm communicator,
                                  RuntimeCandidateIdentity& out) noexcept;
bool valid_runtime_candidate_identity(
    const RuntimeCandidateIdentity& identity) noexcept;
PlanFingerprint runtime_sha256_fingerprint(
    const std::array<char, kRuntimeSha256HexCharacters + 1U>& digest) noexcept;
bool valid_runtime_sha256(const RuntimeSha256Digest& digest) noexcept;
bool runtime_sha256_bytes(Span<const std::uint8_t> bytes,
                          RuntimeSha256Digest& out) noexcept;

}  // namespace hundun::v04::detail
