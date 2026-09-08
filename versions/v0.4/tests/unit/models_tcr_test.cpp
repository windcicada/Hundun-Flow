// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#include "models_tcr_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04;
using namespace tcr::detail;
bool near(double actual, double expected, double absolute = 1e-14,
          double relative = 1e-13) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             absolute + relative * std::abs(expected);
}
int main() {
  // R -> 0: the admissible lower root is eta*R to first order. An
  // absolute tolerance near unity would hide catastrophic cancellation.
  const auto weak_reaction = algebra({0.25, 1e-20}, -1);
  if (weak_reaction.status != Status::success ||
      !near(weak_reaction.control, 2.5e-21, 0, 2e-14))
    return 26;
  // eta=1/4,r=1 -> D=1/4, roots 1/3 and 1, two admissible branches.
  const auto a = algebra({0.25, 1}, -1);
  if (a.status != Status::success || a.admissible_count != 2 ||
      a.signed_root != -0.5 || !near(a.control, 1.0 / 3, 1e-15, 1e-14))
    return 1;
  if (algebra({0.25, 1}).status != Status::branch_required)
    return 2;
  // eta=1/4,r=0 -> roots 0 and 4/3: only zero is physically admissible.
  const auto unique = algebra({0.25, 0});
  if (unique.status != Status::success || unique.admissible_count != 1 ||
      unique.control != 0)
    return 3;
  double y[]{0.2, 0.3, 0.5}, mw[]{2, 3, 5}, rates[]{2, 4};
  std::size_t reactants[]{0, 1};
  ReactantMappingInput mi{{0, 1, 1}, {0, 1, 1}, 42,    42, y, mw,   3,
                          reactants, 2,         rates, 2,  6, 1e-12};
  const auto mapped = ideal_gas_reactant_mole_fraction_v1(mi);
  if (mapped.status != Status::success ||
      !near(mapped.input.eta, 2.0 / 3, 1e-15, 1e-14) ||
      mapped.input.rate_ratio != 0.5)
    return 4;
  if (algebra({0.5, 2}, -1).status != Status::negative_discriminant ||
      algebra({1, 1}, -1).status != Status::weak_denominator ||
      algebra({0.5, 1}, -1).status != Status::fold_unresolved)
    return 5;
  History h;
  h.revision = {0, 1, 1};
  TrialRequest request;
  request.expected_revision = h.revision;
  request.mode = Mode::experimental;
  request.mapping = {Status::success, {0.25, 1}, 71};
  request.initialization_sign = -1;
  auto trial = prepare(h, request);
  auto retry = prepare(h, request);
  if (!trial.available || !trial.feedback ||
      trial.mixer_control != retry.mixer_control || h.initialized)
    return 6;
  auto accepted = accept(h, trial, {1, 2, 1});
  if (!accepted.available || !accepted.candidate.initialized ||
      accepted.candidate.branch_sign != -1)
    return 7;
  h = accepted.candidate;
  request.expected_revision = h.revision;
  request.initialization_sign = 0;
  if (!restore(h, h.revision).available)
    return 8;
  auto corrupted = h;
  corrupted.control = 0.99;
  if (restore(corrupted, h.revision).available)
    return 9;
  request.mapping.input = {0.25, 0.5};
  if (prepare(h, request).candidate.branch_sign != -1)
    return 10;
  request.initialization_sign = 1;
  if (prepare(h, request).status != Status::branch_crossing)
    return 11;
  request.initialization_sign = 0;
  request.mapping.input = {0.25, 1};
  request.fold = {true, h.revision, {0.25, 4.0 / 3}, 1};
  auto folded = prepare(h, request);
  if (!folded.available || folded.candidate.branch_sign != 1 ||
      folded.candidate.fold_count != 1)
    return 12;
  request.fold.at_fold = {0.25, 1.1};
  if (prepare(h, request).available)
    return 13;
  request.fold = {};
  request.mapping.input = {0.5, 2};
  if (prepare(h, request).available)
    return 14;
  request.mode = Mode::shadow;
  auto shadow = prepare(h, request);
  if (!shadow.available || shadow.feedback || shadow.mixer_control != 1 ||
      shadow.candidate.branch_sign != h.branch_sign)
    return 15;
  request.mode = Mode::validated;
  if (prepare(h, request).status != Status::scientific_evidence_unavailable)
    return 16;
  request.mode = Mode::off;
  if (!prepare(h, request).available || prepare(h, request).feedback)
    return 17;
  double xs[]{0, 1, 0, 1}, rs[]{2, 4, 2, 4};
  auto stats = statistics(xs, rs, 4);
  if (stats.status != Status::success || stats.scalar_mean != 0.5 ||
      stats.scalar_variance != 0.25 || stats.scalar_rate_covariance != 0.5)
    return 18;
  if (statistics(xs, rs, 3).status != Status::statistics_unavailable)
    return 19;
  mi.psr_progress_rate = 0;
  if (ideal_gas_reactant_mole_fraction_v1(mi).status !=
      Status::weak_denominator)
    return 20;
  mi.psr_progress_rate = -6;
  if (ideal_gas_reactant_mole_fraction_v1(mi).status == Status::success)
    return 21;
  // One shared mixer: off exactly reproduces IEM; experimental control=1/3
  // has exp(-cbrt(1/3)) variance ratio for dt=tau=1.
  combustion::ChemistryIdentity identity{1, {{1, {1}, 0}, {1, {1}, 0}}, 42};
  double fields[]{0.2, 0.8, 100, 0.8, 0.2, 300};
  esf::detail::Request eq;
  eq.accepted = {h.revision, 42, 2, 2, fields};
  eq.expected_revision = h.revision;
  eq.identity = &identity;
  eq.dt_s = 1;
  eq.mixing_time_s = 1;
  eq.random.accepted_step = h.revision.accepted_step;
  esf::detail::Workspace workspace(2);
  request.mapping.input = {0.25, 1};
  request.mode = Mode::experimental;
  auto combined = mix(h, request, eq, workspace);
  if (!combined.available || !combined.tcr.feedback ||
      combined.tcr.mixer_control != a.control ||
      !near(combined.esf.means[2], 200, 1e-12, 1e-13))
    return 22;
  const double actual = combined.esf.candidate.values[0];
  // Independently evaluated high-precision reference for .5-.3
  // exp(-cbrt(1/3)/2).
  if (!near(actual, 0.287890672504345570578038873963459))
    return 23;
  request.mode = Mode::shadow;
  if (!mix(h, request, eq, workspace).available)
    return 24;
  request.mode = Mode::experimental;
  request.mapping.input = {0.5, 2};
  auto failed = mix(h, request, eq, workspace);
  if (failed.available || failed.esf.candidate.values || h.branch_sign != -1)
    return 25;
  std::cout << "TCR focused PASS\n";
}
