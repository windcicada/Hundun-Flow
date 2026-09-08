// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09 HUNDUN equations explicitly ported from
// governance 8ffdf2b comb_tcr_algebra; independent value history, no inherited
// oracle authorization or checkpoint.
#include "models_tcr_detail.hpp"
#include <cmath>
#include <limits>
namespace hundun::v04::tcr::detail {
namespace {
constexpr double tolerance = 64 * std::numeric_limits<double>::epsilon();
bool close(double a, double b) noexcept {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= tolerance + tolerance * std::abs(b);
}
bool valid_sign(int s) noexcept { return s == -1 || s == 1; }
} // namespace
AlgebraReport algebra(AlgebraInput a, int sign) noexcept {
  AlgebraReport r;
  if (!std::isfinite(a.eta) || a.eta < 0 || a.eta > 1 ||
      !std::isfinite(a.rate_ratio) || a.rate_ratio < 0 ||
      (sign != 0 && !valid_sign(sign)))
    return r;
  const double d = 1 - 4 * a.eta * (1 - a.eta) * a.rate_ratio;
  r.discriminant = d;
  if (!std::isfinite(d))
    return r;
  if (d < 0) {
    r.status = Status::negative_discriminant;
    return r;
  }
  const double denominator = 2 * (1 - a.eta);
  if (std::abs(denominator) <= tolerance) {
    r.status = Status::weak_denominator;
    return r;
  }
  if (d <= tolerance) {
    r.status = Status::fold_unresolved;
    return r;
  }
  const double s = std::sqrt(d);
  // Rationalized lower root retains eta*R as R -> 0. This is algebraic
  // evaluation of the same equation, not clipping or an alternate model.
  r.roots = {2 * a.eta * a.rate_ratio / (1 + s), (1 + s) / denominator};
  for (unsigned i = 0; i < 2; ++i) {
    r.admissible[i] =
        std::isfinite(r.roots[i]) && r.roots[i] >= 0 && r.roots[i] <= 1;
    r.admissible_count += r.admissible[i] ? 1 : 0;
  }
  if (r.admissible_count == 0) {
    r.status = Status::no_admissible_root;
    return r;
  }
  if (sign == 0 && r.admissible_count == 2) {
    r.status = Status::branch_required;
    return r;
  }
  if (sign == 0)
    sign = r.admissible[0] ? -1 : 1;
  const unsigned index = sign < 0 ? 0 : 1;
  if (!r.admissible[index]) {
    r.status = Status::branch_crossing;
    return r;
  }
  r.status = Status::success;
  r.selected_sign = sign;
  r.signed_root = sign * s;
  r.control = r.roots[index];
  return r;
}
Statistics statistics(const double *scalar, const double *rates,
                      std::size_t n) noexcept {
  Statistics r;
  r.count = n;
  if (!scalar || !rates || (n != 2 && n != 4))
    return r;
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(scalar[i]) || !std::isfinite(rates[i]))
      return r;
    r.scalar_mean += scalar[i] / double(n);
    r.rate_mean += rates[i] / double(n);
  }
  for (std::size_t i = 0; i < n; ++i) {
    const double dx = scalar[i] - r.scalar_mean;
    r.scalar_variance += dx * dx / double(n);
    r.scalar_rate_covariance += dx * (rates[i] - r.rate_mean) / double(n);
  }
  if (!std::isfinite(r.scalar_mean) || !std::isfinite(r.rate_mean) ||
      !std::isfinite(r.scalar_variance) ||
      !std::isfinite(r.scalar_rate_covariance))
    return r;
  r.status = Status::success;
  r.minimal_pair_limited = n == 2;
  return r;
}
MappingReport map_statistics(const Statistics &s,
                             const StatisticsMappingProvider &p) noexcept {
  if (s.status != Status::success || !std::isfinite(s.scalar_mean) ||
      !std::isfinite(s.scalar_variance) || s.scalar_variance < 0 ||
      !std::isfinite(s.rate_mean) || !std::isfinite(s.scalar_rate_covariance) ||
      (s.count != 2 && s.count != 4))
    return {Status::statistics_unavailable};
  auto r = p.map(s);
  if (r.status != Status::success)
    return r;
  if (r.mapping_identity == 0 || !std::isfinite(r.input.eta) ||
      !std::isfinite(r.input.rate_ratio))
    return {Status::mapping_unavailable};
  return r;
}
MappingReport
ideal_gas_reactant_mole_fraction_v1(const ReactantMappingInput &q) noexcept {
  MappingReport r;
  if (q.revision != q.expected_revision || q.revision.algorithm_version != 1) {
    r.status = Status::stale_revision;
    return r;
  }
  if (q.composition_fingerprint == 0 ||
      q.composition_fingerprint != q.rate_composition_fingerprint ||
      !q.mean_mass_fractions || !q.molecular_weights_kg_per_kmol ||
      !q.reactant_indices || !q.field_progress_rates || q.species_count == 0 ||
      q.reactant_count == 0 || q.reactant_count > q.species_count ||
      (q.field_count != 2 && q.field_count != 4) ||
      !std::isfinite(q.psr_progress_rate) ||
      !std::isfinite(q.weak_rate_absolute_threshold) ||
      q.weak_rate_absolute_threshold <= 0)
    return r;
  if (std::abs(q.psr_progress_rate) <= q.weak_rate_absolute_threshold) {
    r.status = Status::weak_denominator;
    return r;
  }
  double sum_y = 0, moles = 0, reactant_moles = 0, mean_rate = 0;
  for (std::size_t s = 0; s < q.species_count; ++s) {
    const double y = q.mean_mass_fractions[s],
                 w = q.molecular_weights_kg_per_kmol[s];
    if (!std::isfinite(y) || y < 0 || y > 1 || !std::isfinite(w) || w <= 0)
      return r;
    sum_y += y;
    moles += y / w;
  }
  if (!close(sum_y, 1) || !std::isfinite(moles) || moles <= 0)
    return r;
  for (std::size_t i = 0; i < q.reactant_count; ++i) {
    const auto s = q.reactant_indices[i];
    if (s >= q.species_count)
      return r;
    for (std::size_t j = 0; j < i; ++j)
      if (s == q.reactant_indices[j])
        return r;
  }
  // Canonical species order also makes an all-species set yield eta=1
  // without a clamp, independent of the supplied index ordering.
  for (std::size_t s = 0; s < q.species_count; ++s)
    for (std::size_t i = 0; i < q.reactant_count; ++i)
      if (s == q.reactant_indices[i])
        reactant_moles +=
            q.mean_mass_fractions[s] / q.molecular_weights_kg_per_kmol[s];
  for (std::size_t f = 0; f < q.field_count; ++f) {
    if (!std::isfinite(q.field_progress_rates[f]))
      return r;
    mean_rate += q.field_progress_rates[f] / double(q.field_count);
  }
  if (!std::isfinite(mean_rate) ||
      (mean_rate != 0 &&
       std::signbit(mean_rate) != std::signbit(q.psr_progress_rate)))
    return r;
  r.input = {reactant_moles / moles, mean_rate / q.psr_progress_rate};
  if (!std::isfinite(r.input.eta) || r.input.eta < 0 || r.input.eta > 1 ||
      !std::isfinite(r.input.rate_ratio))
    return {};
  r.status = Status::success;
  r.mapping_identity = kReactantMoleFractionMappingIdentity;
  return r;
}
bool valid_history(const History &h) noexcept {
  if (h.revision.algorithm_version != 1)
    return false;
  if (!h.initialized)
    return h.branch_sign == 0 && h.signed_root == 0 && h.control == 0 &&
           h.mapping_identity == 0 && h.input.eta == 0 &&
           h.input.rate_ratio == 0 && h.initialization_sign == 0 &&
           h.fold_count == 0;
  if (!valid_sign(h.initialization_sign))
    return false;
  if (h.fold_count != 0 &&
      (algebra(h.last_fold).status != Status::fold_unresolved ||
       h.last_fold_base_revision.algorithm_version != 1 ||
       h.last_fold_base_revision.accepted_step > h.revision.accepted_step ||
       h.last_fold_base_revision.input_revision > h.revision.input_revision))
    return false;
  const auto a = algebra(h.input, h.branch_sign);
  return h.mapping_identity != 0 && a.status == Status::success &&
         close(h.signed_root, a.signed_root) && close(h.control, a.control);
}
Trial prepare(const History &h, const TrialRequest &q) noexcept {
  Trial r;
  r.base_revision = h.revision;
  r.mode = q.mode;
  r.request = q;
  if (h.revision != q.expected_revision) {
    r.status = Status::stale_revision;
    return r;
  }
  if (!valid_history(h)) {
    r.status = Status::invalid_input;
    return r;
  }
  if (q.mode != Mode::off && q.mode != Mode::shadow &&
      q.mode != Mode::experimental && q.mode != Mode::validated)
    return r;
  if (q.mode == Mode::validated) {
    r.status = Status::scientific_evidence_unavailable;
    return r;
  }
  if (q.mode == Mode::off) {
    r.status = Status::success;
    r.available = true;
    r.candidate = h;
    return r;
  }
  if (q.mapping.status != Status::success || q.mapping.mapping_identity == 0) {
    r.status = q.mapping.status == Status::success ? Status::mapping_unavailable
                                                   : q.mapping.status;
    if (q.mode == Mode::shadow) {
      r.available = true;
      r.candidate = h;
    }
    return r;
  }
  if (h.initialized && h.mapping_identity != q.mapping.mapping_identity) {
    r.status = Status::mapping_unavailable;
    return r;
  }
  int sign = h.initialized ? h.branch_sign : q.initialization_sign;
  if (h.initialized && q.initialization_sign != 0) {
    r.status = Status::branch_crossing;
    return r;
  }
  if (q.fold.supplied) {
    const auto fold = algebra(q.fold.at_fold);
    // Source-control continuation holds eta fixed; departure sign is explicit.
    if (!h.initialized || q.fold.base_revision != h.revision ||
        !valid_sign(q.fold.departure_sign) ||
        fold.status != Status::fold_unresolved ||
        !close(q.fold.at_fold.eta, h.input.eta) ||
        !close(q.mapping.input.eta, h.input.eta)) {
      r.status = Status::branch_crossing;
      return r;
    }
    sign = q.fold.departure_sign;
  }
  r.observed = algebra(q.mapping.input, sign);
  r.status = r.observed.status;
  if (r.status != Status::success) {
    if (q.mode == Mode::shadow) {
      r.available = true;
      r.candidate = h;
    }
    return r;
  }
  r.candidate = h;
  r.candidate.initialized = true;
  r.candidate.input = q.mapping.input;
  r.candidate.signed_root = r.observed.signed_root;
  r.candidate.control = r.observed.control;
  r.candidate.branch_sign = r.observed.selected_sign;
  r.candidate.mapping_identity = q.mapping.mapping_identity;
  if (!h.initialized)
    r.candidate.initialization_sign = r.observed.selected_sign;
  if (q.fold.supplied) {
    if (h.fold_count == std::numeric_limits<std::uint64_t>::max()) {
      r.status = Status::invalid_input;
      return r;
    }
    r.candidate.fold_count = h.fold_count + 1;
    r.candidate.last_fold = q.fold.at_fold;
    r.candidate.last_fold_base_revision = h.revision;
  }
  r.available = true;
  r.feedback = q.mode == Mode::experimental;
  r.mixer_control = r.feedback ? r.observed.control : 1;
  return r;
}
Trial accept(const History &h, const Trial &t,
             portable::Revision next) noexcept {
  Trial r;
  if (!t.available || !valid_history(h) || !valid_history(t.candidate) ||
      h.revision != t.base_revision || t.candidate.revision != h.revision ||
      next.algorithm_version != 1 ||
      h.revision.accepted_step == std::numeric_limits<std::uint64_t>::max() ||
      next.accepted_step != h.revision.accepted_step + 1 ||
      next.input_revision <= h.revision.input_revision) {
    r.status = Status::stale_revision;
    return r;
  }
  const auto checked = prepare(h, t.request);
  if (!checked.available ||
      checked.candidate.branch_sign != t.candidate.branch_sign ||
      checked.candidate.control != t.candidate.control ||
      checked.candidate.signed_root != t.candidate.signed_root ||
      checked.candidate.input.eta != t.candidate.input.eta ||
      checked.candidate.input.rate_ratio != t.candidate.input.rate_ratio ||
      checked.candidate.mapping_identity != t.candidate.mapping_identity ||
      checked.candidate.fold_count != t.candidate.fold_count) {
    r.status = Status::invalid_input;
    return r;
  }
  r = checked;
  r.candidate.revision = next;
  return r;
}
Trial restore(const History &snapshot, portable::Revision expected) noexcept {
  Trial r;
  if (snapshot.revision != expected) {
    r.status = Status::stale_revision;
    return r;
  }
  if (!valid_history(snapshot))
    return r;
  r.status = Status::success;
  r.available = true;
  r.base_revision = expected;
  r.candidate = snapshot;
  return r;
}
MixReport mix(const History &h, const TrialRequest &t, esf::detail::Request q,
              esf::detail::Workspace &workspace) noexcept {
  MixReport r;
  r.tcr = prepare(h, t);
  if (!r.tcr.available)
    return r;
  if (q.expected_revision != h.revision) {
    r.tcr.status = Status::stale_revision;
    r.tcr.available = false;
    return r;
  }
  q.tcr_control = r.tcr.mixer_control;
  r.esf = workspace.advance(q);
  r.available = r.esf.status == portable::Status::success;
  if (!r.available)
    r.tcr.available = false;
  return r;
}
} // namespace hundun::v04::tcr::detail
