// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_local_flow_pattern.hpp"

#include <array>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::ImmersedLinkId;
using hundun::immersed::LocalCoefficientRow;
using hundun::immersed::LocalFlowPatternTransform;
using hundun::immersed::ReplacementGroup;
using hundun::immersed::RowReplacementPlan;
using hundun::mesh::GlobalCellId;
using hundun::runtime::Real3;

using TransformFull = LocalCoefficientRow (LocalFlowPatternTransform::*)(
    const LocalCoefficientRow &, double, Real3) const;
using PlanRow = RowReplacementPlan (LocalFlowPatternTransform::*)(
    GlobalCellId, const std::vector<ImmersedLinkId> &,
    const LocalCoefficientRow &) const;
using Evaluate = double (LocalFlowPatternTransform::*)(
    const RowReplacementPlan &, const LocalCoefficientRow &,
    const std::vector<double> &) const;

static_assert(std::is_final_v<LocalCoefficientRow>);
static_assert(std::is_same_v<decltype(LocalCoefficientRow::neighbour),
                             std::array<double, 6>>);
static_assert(std::is_same_v<decltype(LocalCoefficientRow::diagonal), double>);
static_assert(std::is_same_v<decltype(LocalCoefficientRow::source), double>);
static_assert(std::is_final_v<ReplacementGroup>);
static_assert(
    std::is_same_v<decltype(ReplacementGroup::stable_group_id), std::uint64_t>);
static_assert(std::is_same_v<decltype(ReplacementGroup::links),
                             std::vector<ImmersedLinkId>>);
static_assert(std::is_same_v<decltype(ReplacementGroup::algebraic_occurrences),
                             std::vector<std::uint32_t>>);
static_assert(std::is_final_v<RowReplacementPlan>);
static_assert(
    std::is_same_v<decltype(RowReplacementPlan::active_cell), GlobalCellId>);
static_assert(std::is_same_v<decltype(RowReplacementPlan::groups),
                             std::vector<ReplacementGroup>>);
static_assert(std::is_final_v<LocalFlowPatternTransform>);
static_assert(
    std::is_same_v<TransformFull,
                   decltype(&LocalFlowPatternTransform::transform_full)>);
static_assert(
    std::is_same_v<PlanRow, decltype(&LocalFlowPatternTransform::plan_row)>);
static_assert(std::is_same_v<
              Evaluate,
              decltype(&LocalFlowPatternTransform::evaluate_wall_replacement)>);
static_assert(noexcept(
    std::declval<const LocalFlowPatternTransform &>().algorithm_fingerprint()));

} // namespace

int main() { return 0; }
