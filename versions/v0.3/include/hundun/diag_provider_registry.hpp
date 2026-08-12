// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_structured.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::diagnostics {

using DiagnosticProviderValidation =
    bool (*)(const void *, std::string &) noexcept;
using DiagnosticProviderFingerprint = std::uint64_t (*)(const void *) noexcept;
using LocalDiagnosticCallback =
    void (*)(const void *, const DiagnosticRequest &, DiagnosticSink &);

struct DiagnosticProvider final {
  DiagnosticModuleKind kind{};
  std::string provider_id;
  const void *source{};
  DiagnosticProviderValidation validate{};
  DiagnosticProviderFingerprint state_fingerprint{};
  LocalDiagnosticCallback collect_local{};
};

class DiagnosticProviderRegistry final {
public:
  void register_provider(DiagnosticProvider);
  void register_providers(std::vector<DiagnosticProvider>);
  const DiagnosticProvider *find(DiagnosticModuleKind) const noexcept;
  const std::vector<DiagnosticProvider> &providers() const noexcept;
  bool collect_local(DiagnosticModuleKind, const DiagnosticRequest &,
                     DiagnosticSink &) const;

private:
  std::vector<DiagnosticProvider> providers_;
};

} // namespace hundun::diagnostics
