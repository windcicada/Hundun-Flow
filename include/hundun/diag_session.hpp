// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace hundun::diagnostics {

class DiagnosticBatch final {
 public:
  std::size_t size() const noexcept;
  bool empty() const noexcept;
  std::string canonical_json_lines() const;

 private:
  friend class DiagnosticBatchSink;
  std::vector<DiagnosticRecord> records_;
};

class DiagnosticBatchSink final : public DiagnosticSink {
 public:
  DiagnosticBatchSink(DiagnosticDescriptor descriptor,
                      DiagnosticRequest request, DiagnosticBatch& batch);
  void submit(const DiagnosticRecord& record) override;

 private:
  DiagnosticDescriptor descriptor_;
  DiagnosticRequest request_;
  DiagnosticBatch* batch_{};
};

class DiagnosticSession final {
 public:
  DiagnosticSession(std::filesystem::path directory, int write_interval,
                    int rank);
  DiagnosticSession(DiagnosticSession&&) noexcept = default;
  DiagnosticSession& operator=(DiagnosticSession&&) noexcept = default;
  DiagnosticSession(const DiagnosticSession&) = delete;
  DiagnosticSession& operator=(const DiagnosticSession&) = delete;

  bool due(std::uint64_t step) const noexcept;
  const std::filesystem::path& directory() const noexcept;
  int write_interval() const noexcept;
  int rank() const noexcept;
  void publish(const runtime::MpiContext& mpi, std::uint64_t step,
               const DiagnosticBatch& batch) const;

 private:
  std::filesystem::path directory_;
  int write_interval_{};
  int rank_{};
};

}  // namespace hundun::diagnostics
