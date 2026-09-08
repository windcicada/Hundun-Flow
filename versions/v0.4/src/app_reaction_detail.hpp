// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/v04_case.hpp"
#include <cmath>
#include <limits>

namespace hundun::v04::detail {
inline bool valid_reaction_spec(const ReactionSpec &r) {
  if (r.mode == ReactionMode::none)
    return r.mechanism_sha256.empty() && r.phase.empty() &&
           r.mechanism_file.empty() &&
           r.representation == ReactionSpec::Representation::external_provider;
  if (r.mode != ReactionMode::finite_rate_mean &&
      r.mode != ReactionMode::pasr_algebraic_v1 &&
      r.mode != ReactionMode::esf_tpdf)
    return false;
  if (r.mechanism_sha256.size() != 64 ||
      r.mechanism_sha256.find_first_not_of("0123456789abcdef") !=
          std::string::npos ||
      r.phase.empty() || r.phase.size() > 255 || r.phase.find('\0') != std::string::npos ||
      !std::isfinite(r.relative_tolerance) || r.relative_tolerance <= 0 ||
      r.relative_tolerance >= 1 || !std::isfinite(r.absolute_tolerance) ||
      r.absolute_tolerance <= 0 || r.maximum_internal_steps == 0 ||
      r.maximum_internal_steps >
          std::uint32_t(std::numeric_limits<int>::max()) ||
      !std::isfinite(r.mixing_c_z) || r.mixing_c_z <= 0 ||
      !std::isfinite(r.turbulent_schmidt) || r.turbulent_schmidt <= 0 ||
      !std::isfinite(r.analytic_rate_s) || r.analytic_rate_s < 0 ||
      (r.analytic_cp_j_per_kg_k != 1000 && r.analytic_cp_j_per_kg_k != 1200))
    return false;
  switch (r.representation) {
  case ReactionSpec::Representation::external_provider:
  case ReactionSpec::Representation::analytic_isomer:
    return r.mechanism_file.empty();
  case ReactionSpec::Representation::direct_cantera:
    // The case reader additionally opens a bounded direct-root regular file.
    return !r.mechanism_file.empty() && !r.mechanism_file.has_parent_path() &&
           r.mechanism_file.extension() == ".yaml" &&
           r.mechanism_file.native().size() <= 255;
  }
  return false;
}
} // namespace hundun::v04::detail
