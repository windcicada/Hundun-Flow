// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_provider_registry.hpp"

#include "tests/support/test_main.hpp"

#include <stdexcept>
#include <string>

namespace {

using hundun::diagnostics::DiagnosticModuleKind;
using hundun::diagnostics::DiagnosticProvider;
using hundun::diagnostics::DiagnosticProviderRegistry;
using hundun::diagnostics::DiagnosticRequest;
using hundun::diagnostics::DiagnosticScope;
using hundun::diagnostics::DiagnosticSink;

struct Source final {
  mutable std::uint64_t value{};
  mutable int calls{};
};

class Sink final : public DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &) override {
    ++submissions;
  }
  int submissions{};
};

bool validate_source(const void *source, std::string &) noexcept {
  return source != nullptr;
}
std::uint64_t fingerprint(const void *source) noexcept {
  return static_cast<const Source *>(source)->value;
}
void collect_stable(const void *source, const DiagnosticRequest &,
                    DiagnosticSink &sink) {
  ++static_cast<const Source *>(source)->calls;
  sink.submit({});
}
void collect_mutating(const void *source, const DiagnosticRequest &,
                      DiagnosticSink &sink) {
  sink.submit({});
  ++static_cast<const Source *>(source)->value;
}

DiagnosticProvider provider(DiagnosticModuleKind kind, Source &source,
                            hundun::diagnostics::LocalDiagnosticCallback cb) {
  return {kind, "fixture", &source, &validate_source, &fingerprint, cb};
}

void test_duplicate_kind_and_deterministic_order() {
  Source a;
  Source b;
  DiagnosticProviderRegistry registry;
  registry.register_providers(
      {provider(DiagnosticModuleKind::flow_driver, a, &collect_stable),
       provider(DiagnosticModuleKind::checkpoint, b, &collect_stable)});
  HUNDUN_CHECK(registry.providers()[0].kind ==
               DiagnosticModuleKind::checkpoint);
  HUNDUN_CHECK(registry.providers()[1].kind ==
               DiagnosticModuleKind::flow_driver);
  bool rejected = false;
  try {
    registry.register_provider(
        provider(DiagnosticModuleKind::checkpoint, a, &collect_stable));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(registry.providers().size() == 2U);
}

void test_absence_locality_and_mutation_guard() {
  Source source;
  DiagnosticProviderRegistry registry;
  registry.register_provider(
      provider(DiagnosticModuleKind::checkpoint, source, &collect_stable));
  DiagnosticRequest request;
  request.scope = DiagnosticScope::local;
  Sink sink;
  HUNDUN_CHECK(registry.collect_local(DiagnosticModuleKind::checkpoint,
                                      request, sink));
  HUNDUN_CHECK(source.calls == 1);
  HUNDUN_CHECK(sink.submissions == 1);
  HUNDUN_CHECK(!registry.collect_local(DiagnosticModuleKind::performance,
                                       request, sink));
  HUNDUN_CHECK(source.calls == 1);

  DiagnosticProviderRegistry mutating;
  mutating.register_provider(
      provider(DiagnosticModuleKind::checkpoint, source, &collect_mutating));
  bool rejected = false;
  try {
    static_cast<void>(mutating.collect_local(DiagnosticModuleKind::checkpoint,
                                             request, sink));
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(sink.submissions == 1);

  request.scope = DiagnosticScope::collective;
  rejected = false;
  try {
    static_cast<void>(registry.collect_local(DiagnosticModuleKind::checkpoint,
                                             request, sink));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  HUNDUN_CHECK(source.calls == 1);
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_duplicate_kind_and_deterministic_order();
    test_absence_locality_and_mutation_guard();
  });
}
