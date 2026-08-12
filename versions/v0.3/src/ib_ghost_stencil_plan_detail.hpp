// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/ib_ghost_stencil_plan.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::immersed::detail {

struct ReconstructionFunctionalScore final {
  double amplification{};
  double condition_estimate{};
  std::size_t donor_count{};
  std::uint64_t pivot_fingerprint{};
};

bool functional_score_less(const ReconstructionFunctionalScore &,
                           const ReconstructionFunctionalScore &);
std::size_t select_minimum_functional_score(
    const std::vector<ReconstructionFunctionalScore> &);

int select_wall_quadrature_execution_owner(
    const std::vector<runtime::Box3> &owner_boxes,
    const std::vector<runtime::Int3> &point_donors,
    const std::vector<runtime::Int3> &pressure_authority_donors, int reach);

} // namespace hundun::immersed::detail
