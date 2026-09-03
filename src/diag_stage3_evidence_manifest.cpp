// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_stage3_evidence_manifest.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string_view>

namespace hundun::diagnostics {
namespace {

void required(const std::string &value, std::string_view name) {
  if (value.empty())
    throw std::invalid_argument(std::string(name) + " must not be empty");
}

void string(std::string &out, std::string_view value) {
  out.push_back('"');
  for (const char character : value) {
    if (character == '"' || character == '\\')
      out.push_back('\\');
    if (character == '\n')
      out += "\\n";
    else if (character == '\r')
      out += "\\r";
    else if (character == '\t')
      out += "\\t";
    else
      out.push_back(character);
  }
  out.push_back('"');
}

} // namespace

void validate_stage3_evidence_manifest(
    const std::vector<Stage3EvidenceRecord> &records) {
  if (records.empty())
    throw std::invalid_argument("Stage 3 evidence manifest is empty");
  std::set<std::string> rows;
  for (const auto &record : records) {
    required(record.row_id, "row_id");
    required(record.candidate_head, "candidate_head");
    required(record.tree_fingerprint, "tree_fingerprint");
    required(record.diff_fingerprint, "diff_fingerprint");
    required(record.build_role, "build_role");
    required(record.build_root, "build_root");
    required(record.cache_sha256, "cache_sha256");
    required(record.binary_sha256, "binary_sha256");
    required(record.compiler_identity, "compiler_identity");
    required(record.libcxx_identity, "libcxx_identity");
    required(record.mpi_identity, "mpi_identity");
    required(record.argv, "argv");
    required(record.environment, "environment");
    required(record.cpuset, "cpuset");
    required(record.started_at, "started_at");
    required(record.ended_at, "ended_at");
    required(record.log_sha256, "log_sha256");
    if (record.binary_inode == 0U || record.ranks <= 0 ||
        record.exit_status < 0 || !std::isfinite(record.duration_seconds) ||
        record.duration_seconds < 0.0 || record.artifact_sha256.empty())
      throw std::invalid_argument("Stage 3 evidence terminal fields are invalid");
    if (!rows.insert(record.row_id).second)
      throw std::invalid_argument("duplicate Stage 3 evidence row_id");
    for (const auto &artifact : record.artifact_sha256)
      required(artifact, "artifact_sha256");
  }
}

std::string stage3_evidence_manifest_json(
    const std::vector<Stage3EvidenceRecord> &records) {
  validate_stage3_evidence_manifest(records);
  std::string out = "{\"schema_version\":1,\"records\":[";
  bool first = true;
  for (const auto &record : records) {
    if (!first)
      out.push_back(',');
    first = false;
    out += "{\"row_id\":";
    string(out, record.row_id);
    out += ",\"candidate_head\":";
    string(out, record.candidate_head);
    out += ",\"tree_fingerprint\":";
    string(out, record.tree_fingerprint);
    out += ",\"diff_fingerprint\":";
    string(out, record.diff_fingerprint);
    out += ",\"build_role\":";
    string(out, record.build_role);
    out += ",\"build_root\":";
    string(out, record.build_root);
    out += ",\"cache_sha256\":";
    string(out, record.cache_sha256);
    out += ",\"binary_sha256\":";
    string(out, record.binary_sha256);
    out += ",\"binary_inode\":" + std::to_string(record.binary_inode);
    out += ",\"compiler_identity\":";
    string(out, record.compiler_identity);
    out += ",\"libcxx_identity\":";
    string(out, record.libcxx_identity);
    out += ",\"mpi_identity\":";
    string(out, record.mpi_identity);
    out += ",\"argv\":";
    string(out, record.argv);
    out += ",\"environment\":";
    string(out, record.environment);
    out += ",\"cpuset\":";
    string(out, record.cpuset);
    out += ",\"ranks\":" + std::to_string(record.ranks);
    out += ",\"started_at\":";
    string(out, record.started_at);
    out += ",\"ended_at\":";
    string(out, record.ended_at);
    out += ",\"exit_status\":" + std::to_string(record.exit_status);
    out += ",\"duration_seconds\":" + std::to_string(record.duration_seconds);
    out += ",\"peak_rss_kib\":" + std::to_string(record.peak_rss_kib);
    out += ",\"log_sha256\":";
    string(out, record.log_sha256);
    out += ",\"artifact_sha256\":[";
    bool first_artifact = true;
    for (const auto &artifact : record.artifact_sha256) {
      if (!first_artifact)
        out.push_back(',');
      first_artifact = false;
      string(out, artifact);
    }
    out += "]}";
  }
  out += "]}";
  return out;
}

} // namespace hundun::diagnostics
