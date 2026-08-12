// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::diagnostics {

struct Stage3EvidenceRecord final {
  std::string row_id;
  std::string candidate_head;
  std::string tree_fingerprint;
  std::string diff_fingerprint;
  std::string build_role;
  std::string build_root;
  std::string cache_sha256;
  std::string binary_sha256;
  std::uint64_t binary_inode{};
  std::string compiler_identity;
  std::string libcxx_identity;
  std::string mpi_identity;
  std::string argv;
  std::string environment;
  std::string cpuset;
  int ranks{};
  std::string started_at;
  std::string ended_at;
  int exit_status{-1};
  double duration_seconds{};
  std::uint64_t peak_rss_kib{};
  std::string log_sha256;
  std::vector<std::string> artifact_sha256;
};

void validate_stage3_evidence_manifest(
    const std::vector<Stage3EvidenceRecord> &records);
std::string stage3_evidence_manifest_json(
    const std::vector<Stage3EvidenceRecord> &records);

} // namespace hundun::diagnostics
