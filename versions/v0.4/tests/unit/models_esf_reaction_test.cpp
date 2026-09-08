// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#include "models_esf_detail.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
namespace {
std::size_t allocation_count{};
}
void *operator new(std::size_t n) {
  ++allocation_count;
  if (void *p = std::malloc(n ? n : 1))
    return p;
  throw std::bad_alloc();
}
void operator delete(void *p) noexcept { std::free(p); }
using namespace hundun::v04;
bool near(double actual, double expected, double absolute = 1e-14,
          double relative = 1e-13) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             absolute + relative * std::abs(expected);
}
struct Analytic final : portable::GasAdvanceProvider {
  const portable::GasIdentity *identity{};
  const portable::GasIdentity &gas_identity() const noexcept override {
    return *identity;
  }
  unsigned calls{}, fail_on{};
  bool bad_heat{}, bad_identity{}, bad_elements{};
  portable::Status
  advance_gas(const portable::GasAdvanceQuery &q,
              portable::GasAdvanceOutput &o) noexcept override {
    if (++calls == fail_on)
      return portable::Status::provider_failure;
    const double a = q.state.mass_fractions[0] * std::exp(-q.duration_s);
    o.final_mass_fractions[0] = a;
    o.final_mass_fractions[1] = 1 - a;
    o.integrated_species_density_delta_kg_per_m3[0] =
        a - q.state.mass_fractions[0];
    o.integrated_species_density_delta_kg_per_m3[1] =
        q.state.mass_fractions[0] - a;
    o.final_sample = {q.state.revision,
                      q.state.composition_fingerprint,
                      q.state.pressure_pa,
                      300,
                      1,
                      q.state.enthalpy_j_per_kg,
                      1000,
                      1e-5,
                      0.1};
    o.integrated_heat_release_j_per_m3 =
        -10 * o.integrated_species_density_delta_kg_per_m3[0];
    o.completed_duration_s = q.duration_s;
    o.internal_step_count = 1;
    if (bad_heat)
      o.integrated_heat_release_j_per_m3 += 1;
    if (bad_identity)
      ++o.final_sample.composition_fingerprint;
    if (bad_elements)
      o.integrated_species_density_delta_kg_per_m3[0] += 1;
    return portable::Status::success;
  }
};
int main() {
  combustion::ChemistryIdentity id{1, {{1, {1}, 10}, {1, {1}, 0}}, 42};
  portable::GasIdentity gas{
      "synthetic-sha",     "gas", {"A", "B"}, {"E"}, {1, 1}, {1, 1},
      "synthetic-total-h", 17,    42};
  double values[]{.2, .8, 100, .8, .2, 100, .2, .8, 100, .8, .2, 100},
      p[]{1e5, 1e5, 1e5, 1e5}, rho[]{1, 1, 1, 1};
  esf::detail::ReactionRequest q{
      {{0, 1, 1}, 42, 4, 2, values}, {0, 1, 1}, &id, &gas, p, rho, 0, 1};
  esf::detail::Workspace workspace(2);
  Analytic backend;
  backend.identity = &gas;
  auto r = workspace.react(q, backend);
  if (r.status != portable::Status::success || r.chemistry_call_count != 8 ||
      backend.calls != 8 || !near(r.candidate.values[0], 0.07357588823428846) ||
      r.candidate.values[2] != 100)
    return 1;
  if (!near(r.ensemble_heat_release_j_per_m3, 3.1606027941427884, 1e-12, 1e-13))
    return 2;
  backend.calls = 0;
  backend.fail_on = 3;
  auto failed = workspace.react(q, backend);
  if (failed.status != portable::Status::provider_failure ||
      failed.candidate.values || values[0] != 0.2 ||
      workspace.valid(r.candidate))
    return 3;
  backend.calls = 0;
  backend.fail_on = 0;
  q.accepted.fields = 2;
  if (workspace.react(q, backend).chemistry_call_count != 4)
    return 4;
  auto allocations = allocation_count;
  backend.calls = 0;
  auto repeat = workspace.react(q, backend);
  if (repeat.status != portable::Status::success ||
      allocation_count != allocations)
    return 5;
  backend.bad_heat = true;
  if (workspace.react(q, backend).status !=
      portable::Status::conservation_failure)
    return 6;
  backend.bad_heat = false;
  backend.bad_identity = true;
  if (workspace.react(q, backend).status != portable::Status::provider_failure)
    return 7;
  backend.bad_identity = false;
  backend.bad_elements = true;
  if (workspace.react(q, backend).status !=
      portable::Status::conservation_failure)
    return 8;
  backend.bad_elements = false;
  id.species[1].element_counts[0] = 2;
  gas.element_counts[1] = 2;
  if (workspace.react(q, backend).status !=
      portable::Status::conservation_failure)
    return 9;
  std::cout << "ESF persistent reaction PASS\n";
}
