// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_provider_registry.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace hundun::diagnostics {
namespace {

class BufferedSink final : public DiagnosticSink {
public:
  void submit(const DiagnosticRecord &record) override {
    records.push_back(record);
  }

  std::vector<DiagnosticRecord> records;
};

std::uint16_t kind_value(DiagnosticModuleKind kind) noexcept {
  return static_cast<std::uint16_t>(kind);
}

void validate_provider(const DiagnosticProvider &provider) {
  if (provider.provider_id.empty() || provider.validate == nullptr ||
      provider.state_fingerprint == nullptr ||
      provider.collect_local == nullptr) {
    throw std::invalid_argument("diagnostic provider is incomplete");
  }
  std::string message;
  if (!provider.validate(provider.source, message)) {
    throw std::invalid_argument(message.empty()
                                    ? "diagnostic provider validation failed"
                                    : std::move(message));
  }
}

} // namespace

void DiagnosticProviderRegistry::register_provider(DiagnosticProvider provider) {
  register_providers({std::move(provider)});
}

void DiagnosticProviderRegistry::register_providers(
    std::vector<DiagnosticProvider> providers) {
  std::vector<DiagnosticProvider> candidate = providers_;
  candidate.insert(candidate.end(), std::make_move_iterator(providers.begin()),
                   std::make_move_iterator(providers.end()));
  for (const auto &provider : candidate) {
    validate_provider(provider);
  }
  std::sort(candidate.begin(), candidate.end(),
            [](const auto &left, const auto &right) {
              return kind_value(left.kind) < kind_value(right.kind);
            });
  if (std::adjacent_find(candidate.begin(), candidate.end(),
                         [](const auto &left, const auto &right) {
                           return left.kind == right.kind;
                         }) != candidate.end()) {
    throw std::invalid_argument("duplicate diagnostic provider kind");
  }
  providers_.swap(candidate);
}

const DiagnosticProvider *
DiagnosticProviderRegistry::find(DiagnosticModuleKind kind) const noexcept {
  const auto found = std::lower_bound(
      providers_.begin(), providers_.end(), kind,
      [](const auto &provider, DiagnosticModuleKind value) {
        return kind_value(provider.kind) < kind_value(value);
      });
  return found != providers_.end() && found->kind == kind ? &*found : nullptr;
}

const std::vector<DiagnosticProvider> &
DiagnosticProviderRegistry::providers() const noexcept {
  return providers_;
}

bool DiagnosticProviderRegistry::collect_local(
    DiagnosticModuleKind kind, const DiagnosticRequest &request,
    DiagnosticSink &sink) const {
  if (request.scope != DiagnosticScope::local) {
    throw std::invalid_argument(
        "local diagnostic provider cannot execute a collective request");
  }
  const auto *provider = find(kind);
  if (provider == nullptr) {
    return false;
  }
  const std::uint64_t before =
      provider->state_fingerprint(provider->source);
  BufferedSink buffered;
  provider->collect_local(provider->source, request, buffered);
  const std::uint64_t after = provider->state_fingerprint(provider->source);
  if (before != after) {
    throw std::runtime_error("diagnostic provider mutated observed state");
  }
  for (const auto &record : buffered.records) {
    sink.submit(record);
  }
  return true;
}

} // namespace hundun::diagnostics
