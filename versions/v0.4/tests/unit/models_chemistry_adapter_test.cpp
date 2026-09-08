// SPDX-License-Identifier: Apache-2.0
#include "models_chemistry_adapter_detail.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
using namespace hundun::v04;
static bool deny_allocations = false;
void *operator new(std::size_t n) {
  if (deny_allocations)
    throw std::bad_alloc();
  if (auto *p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}
void operator delete(void *p) noexcept { std::free(p); }
class FailingProvider final : public portable::GasQueryProvider {
public:
  chemistry::detail::AnalyticIsomerBackend &backend;
  int calls{}, fail_at{2};
  bool stale{};
  explicit FailingProvider(chemistry::detail::AnalyticIsomerBackend &b)
      : backend(b) {}
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &o) noexcept override {
    ++calls;
    if (calls == fail_at) {
      o.sample = {};
      return portable::Status::provider_failure;
    }
    auto result = backend.query_gas(q, o);
    if (stale)
      ++o.sample.revision.input_revision;
    return result;
  }
};
bool check(bool ok, const char *message) {
  if (!ok)
    std::cerr << "FAIL: " << message << '\n';
  return ok;
}
int main() {
  chemistry::detail::AnalyticIsomerBackend backend;
  const double y[]{1.0, 0.0};
  portable::GasQuery q{{4, 9, 1},
                       backend.composition().fingerprint,
                       portable::GasStateCoordinates::pressure_temperature,
                       101325,
                       0,
                       300,
                       y,
                       2};
  double d[2]{}, h[2]{}, w[2]{};
  portable::GasQueryOutput out{{}, d, h, w, 2};
  // cp=1000; h_A(298.15)=100000, h_B(298.15)=0 J/kg.
  bool ok = check(backend.query_gas(q, out) == portable::Status::success,
                  "synthetic gas PT query succeeds");
  ok &= check(std::abs(out.sample.enthalpy_j_per_kg - 101850) < 1e-8,
              "PT enthalpy includes reference formation enthalpy");
  chemistry::detail::BackendAdapter adapter(
      backend, backend, backend.closure_identity(), q.revision);
  combustion::ChemistryAdvanceRequest advance;
  advance.state = {101325, out.sample.density_kg_per_m3, 101850, {1, 0}};
  advance.duration_s = 0.5;
  auto result = adapter.integrate(advance);
  ok &= check(result.succeeded(), "common adapter advances analytic backend");
  ok &= check(std::abs(result.final_state.mass_fractions[0] -
                       0.36787944117144233) < 1e-13,
              "first-order exact A decay at k dt=1");
  ok &=
      check(result.final_state.total_thermochemical_enthalpy_j_per_kg == 101850,
            "heat release is report only, no duplicate total enthalpy source");
  auto ph = q;
  ph.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  ph.enthalpy_j_per_kg = 101850;
  double final_y[2], delta_y[2];
  portable::GasAdvanceOutput bounded{{}, final_y, delta_y, 2};
  deny_allocations = true;
  const auto bounded_status = backend.advance_gas({ph, 0, 0.5}, bounded);
  deny_allocations = false;
  ok &= check(bounded_status == portable::Status::success &&
                  std::abs(final_y[0] - 0.36787944117144233) < 1e-13,
              "prepared analytic interval has no hot allocations");
  chemistry::detail::GasBatchWorkspace batch(2, 2);
  portable::GasQuery queries[]{q, q};
  auto batched = batch.query(backend, queries, 2, q.revision);
  ok &= check(batched.available && batch.results(),
              "prepared batch exposes complete result");
  queries[1].revision.input_revision++;
  batched = batch.query(backend, queries, 2, q.revision);
  ok &= check(!batched.available && batched.failure_index == 1 &&
                  batched.status == portable::Status::stale_revision &&
                  !batch.results(),
              "stale second input invalidates entire batch");
  queries[1] = q;
  deny_allocations = true;
  batched = batch.query(backend, queries, 2, q.revision);
  deny_allocations = false;
  ok &= check(batched.available, "hot batch does not allocate HUNDUN vectors");
  batched = batch.query(backend, queries, 3, q.revision);
  ok &= check(batched.status == portable::Status::capacity_exceeded &&
                  !batch.results(),
              "batch capacity rejects before reading excess input");
  FailingProvider failure(backend);
  batched = batch.query(failure, queries, 2, q.revision);
  ok &= check(!batched.available && batched.failure_index == 1 &&
                  !batch.results(),
              "second provider failure exposes no partial batch");
  failure.calls = 0;
  failure.fail_at = 0;
  failure.stale = true;
  chemistry::detail::BackendAdapter stale_adapter(
      backend, failure, backend.closure_identity(), q.revision);
  ok &= check(!stale_adapter.integrate(advance).succeeded(),
              "adapter rejects stale provider answer");
  combustion::CombustionClosureRequest closure;
  closure.composition_fingerprint = adapter.identity().fingerprint;
  closure.mean_state = advance.state;
  closure.duration_s = 0.01;
  combustion::CombustionClosureConfig cfg;
  auto finite =
      combustion::evaluate_combustion_closure(closure, cfg, adapter, &adapter);
  ok &= check(finite.succeeded() && finite.chemistry_call_count == 2,
              "analytic finite rate traverses common closure");
  cfg.tci_closure = combustion::TciClosureKind::pasr_algebraic_v1;
  cfg.mixing_time = {0.01, 1e-5, 1e-5, 1, 0.5};
  auto pasr =
      combustion::evaluate_combustion_closure(closure, cfg, adapter, &adapter);
  ok &= check(pasr.succeeded() && pasr.chemistry_call_count == 2 &&
                  pasr.timescales.kappa >= 0 && pasr.timescales.kappa <= 1,
              "analytic PaSR topology and bounds");
  ok &= check(std::abs(pasr.candidate.integrated_heat_release_j_per_m3 -
                       pasr.timescales.kappa *
                           finite.candidate.integrated_heat_release_j_per_m3) <
                  1e-10,
              "PaSR scales heat with identical reaction fraction");
  cfg.tci_closure = combustion::TciClosureKind::esf_tpdf;
  closure.stochastic_fields.assign(4, closure.mean_state);
  auto ensemble =
      combustion::evaluate_combustion_closure(closure, cfg, adapter, &adapter);
  ok &= check(ensemble.succeeded() && ensemble.chemistry_call_count == 8,
              "N4 ensemble topology 2N");
  chemistry::detail::AnalyticIsomerBackend reversed(2, true);
  chemistry::detail::BackendAdapter mapped(
      reversed, reversed, reversed.closure_identity(), q.revision);
  auto reordered = advance;
  reordered.state.mass_fractions = {0, 1};
  auto reversed_result = mapped.integrate(reordered);
  ok &= check(reversed_result.succeeded() &&
                  std::abs(reversed_result.final_state.mass_fractions[1] -
                           0.36787944117144233) < 1e-13,
              "reordered synthetic representation preserves physical result");
  bool rejected = false;
  try {
    chemistry::detail::BackendAdapter wrong(
        backend, reversed, backend.closure_identity(), q.revision);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  ok &= check(rejected, "unmapped species order is not silently accepted");
  chemistry::detail::AnalyticIsomerBackend different_rate(3, false);
  rejected = false;
  try {
    chemistry::detail::BackendAdapter wrong_rate(
        different_rate, backend, backend.closure_identity(), q.revision);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  ok &= check(rejected, "same composition does not authorize pairing different "
                        "chemical mechanism configurations");
  chemistry::detail::AnalyticIsomerBackend different_thermo(2, false, 1200);
  rejected = false;
  try {
    chemistry::detail::BackendAdapter wrong_asset(
        different_thermo, backend, backend.closure_identity(), q.revision);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  ok &= check(rejected, "same species and formation enthalpies do not "
                        "substitute for mechanism SHA identity");
  return ok ? 0 : 1;
}
