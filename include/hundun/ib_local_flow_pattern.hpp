// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_surface.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_types.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace hundun::immersed {

struct LocalCoefficientRow final {
  std::array<double, 6> neighbour{};
  double diagonal{};
  double source{};
};

struct ReplacementGroup final {
  std::uint64_t stable_group_id{};
  std::vector<ImmersedLinkId> links;
  std::vector<std::uint32_t> algebraic_occurrences;
  std::uint64_t evaluator_fingerprint{};
};

struct RowReplacementPlan final {
  mesh::GlobalCellId active_cell{};
  std::vector<ReplacementGroup> groups;
  std::uint64_t fingerprint{};
};

class LocalFlowPatternTransform final {
public:
  LocalCoefficientRow
  transform_full(const LocalCoefficientRow &, double normal_scale,
                 runtime::Real3 solid_to_fluid_normal) const;
  RowReplacementPlan plan_row(mesh::GlobalCellId,
                              const std::vector<ImmersedLinkId> &,
                              const LocalCoefficientRow &) const;
  double evaluate_wall_replacement(
      const RowReplacementPlan &, const LocalCoefficientRow &immutable_snapshot,
      const std::vector<double> &link_local_symbols) const;
  std::uint64_t algorithm_fingerprint() const noexcept;
};

} // namespace hundun::immersed
