// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_combustion.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <new>

namespace {
int allocation_to_fail = -1;
}
void* operator new(std::size_t bytes) {
  if (allocation_to_fail == 0) {
    allocation_to_fail = -1;
    throw std::bad_alloc{};
  }
  if (allocation_to_fail > 0) --allocation_to_fail;
  if (void* result = std::malloc(bytes ? bytes : 1U)) return result;
  throw std::bad_alloc{};
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

using namespace hundun::v04::combustion;
class QueryOnlyBackend final : public ChemistryRepresentationAdapter {
 public:
  QueryOnlyBackend() {
    identity_.element_count = 1U;
    identity_.species = {{1.0, {1U}, 0.0}};
    identity_.fingerprint = chemistry_identity_fingerprint(identity_);
  }
  const ChemistryIdentity& identity() const noexcept override { return identity_; }
  ChemistryRepresentationKind representation() const noexcept override {
    return ChemistryRepresentationKind::direct_cantera;
  }
  ChemistryAdvanceReport integrate(const ChemistryAdvanceRequest&) noexcept override {
    ++calls;
    return {}; // Every injected allocation precedes the first backend call.
  }
  int calls{};
 private:
  ChemistryIdentity identity_;
};

int main() {
  // Make pre-fix noexcept termination observable without generating a core file.
  std::set_terminate([] { std::_Exit(90); });
  QueryOnlyBackend backend;
  CombustionClosureRequest request;
  request.composition_fingerprint = backend.identity().fingerprint;
  request.mean_state = {101325.0, 1.0, 0.0, {1.0}};
  request.duration_s = 1.0;
  request.accepted_step = 7U;
  request.source_state_revision = 11U;
  const CombustionClosureConfig config;
  for (int location = 0; location < 4; ++location) {
    allocation_to_fail = location;
    const auto report = evaluate_combustion_closure(request, config, backend, nullptr);
    allocation_to_fail = -1;
    if (report.status != CombustionClosureStatus::workspace_failure ||
        report.succeeded() || report.candidate.available || report.source_published ||
        !report.candidate.species_density_delta_kg_per_m3.empty() ||
        report.candidate.integrated_heat_release_j_per_m3 != 0.0 ||
        report.accepted_step != 7U || report.source_state_revision != 11U ||
        backend.calls != 0) {
      std::cerr << "closure-owned allocation failure must return empty candidate\n";
      return 1;
    }
  }
  const std::vector<ThermochemicalState> fields(4U, request.mean_state);
  EsfCommonSource source;
  source.composition_fingerprint = backend.identity().fingerprint;
  source.species_density_delta_kg_per_m3 = {0.0};
  for (int location = 0; location < 6; ++location) {
    allocation_to_fail = location;
    const auto report = apply_esf_common_source(
        request.mean_state, fields, backend.identity(), source);
    allocation_to_fail = -1;
    if (report.status != EsfCommonSourceStatus::workspace_failure ||
        report.available || !report.fields.empty() ||
        !report.mean_state.mass_fractions.empty()) {
      std::cerr << "allocation failure must withdraw every common-source field\n";
      return 1;
    }
  }
  return 0;
}
