// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_local_flow_pattern.hpp"

#include "hundun/rt_error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace hundun::immersed {
namespace {

constexpr std::uint64_t kAlgorithmFingerprint = 0x4c46505f56303031ULL;
constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr double kPi = 3.141592653589793238462643383279502884;

using runtime::Real3;

[[noreturn]] void fail(const std::string &message) {
  throw runtime::Error("local flow pattern: " + message);
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) {
  for (std::size_t byte = 0; byte < sizeof(value); ++byte) {
    hash ^= (value >> (byte * 8U)) & 0xffU;
    hash *= kFnvPrime;
  }
}

void hash_u32(std::uint64_t &hash, std::uint32_t value) {
  hash_u64(hash, static_cast<std::uint64_t>(value));
}

bool finite(Real3 value) {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

void require_finite(const LocalCoefficientRow &row) {
  for (const double coefficient : row.neighbour) {
    if (!std::isfinite(coefficient)) {
      fail("coefficient row must be finite");
    }
  }
  if (!std::isfinite(row.diagonal) || !std::isfinite(row.source)) {
    fail("coefficient row must be finite");
  }
}

Real3 cross(Real3 lhs, Real3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

Real3 normalized(Real3 value) {
  if (!finite(value)) {
    fail("normal must be finite");
  }
  const double scale =
      std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
  if (scale == 0.0) {
    fail("normal must be nonzero");
  }
  const Real3 scaled{value.x / scale, value.y / scale, value.z / scale};
  const double magnitude = std::sqrt(scaled.x * scaled.x + scaled.y * scaled.y +
                                     scaled.z * scaled.z);
  return {scaled.x / magnitude, scaled.y / magnitude, scaled.z / magnitude};
}

std::array<double, 3> rotation_factors(Real3 normal) {
  const std::array<double, 3> absolute{std::abs(normal.x), std::abs(normal.y),
                                       std::abs(normal.z)};
  std::size_t axis_index = 0U;
  if (absolute[1] < absolute[axis_index]) {
    axis_index = 1U;
  }
  if (absolute[2] < absolute[axis_index]) {
    axis_index = 2U;
  }
  std::array<Real3, 3> axes{Real3{1.0, 0.0, 0.0}, Real3{0.0, 1.0, 0.0},
                            Real3{0.0, 0.0, 1.0}};
  const Real3 tangent_one = normalized(cross(axes[axis_index], normal));
  const Real3 tangent_two = cross(normal, tangent_one);

  // Deterministic Z-Y-X decomposition of the right-handed local frame.
  const double pitch = std::asin(std::clamp(-tangent_one.z, -1.0, 1.0));
  const double roll = std::atan2(tangent_two.z, normal.z);
  const double yaw = std::atan2(tangent_one.y, tangent_one.x);
  const auto factor = [](double angle) {
    double folded = std::abs(std::remainder(angle, kPi));
    if (folded > 0.5 * kPi) {
      folded = kPi - folded;
    }
    return std::clamp(1.0 - 2.0 * folded / kPi, 0.0, 1.0);
  };
  return {factor(roll), factor(pitch), factor(yaw)};
}

void swap_direction_pairs(LocalCoefficientRow &row, Real3 normal) {
  if (normal.x < 0.0) {
    std::swap(row.neighbour[2], row.neighbour[3]);
  }
  if (normal.y < 0.0) {
    std::swap(row.neighbour[0], row.neighbour[1]);
  }
  if (normal.z < 0.0) {
    std::swap(row.neighbour[4], row.neighbour[5]);
  }
}

LocalCoefficientRow paper_transform_factors(const LocalCoefficientRow &input,
                                            double k0, double k1, double k2,
                                            double k3) {
  require_finite(input);
  for (const double factor : {k0, k1, k2, k3}) {
    if (!std::isfinite(factor) || factor < 0.0 || factor > 1.0) {
      fail("scale and rotation factors must be finite and in [0,1]");
    }
  }

  std::array<double, 7> level_one{};
  std::copy(input.neighbour.begin(), input.neighbour.end(), level_one.begin());
  level_one[6] = input.diagonal;

  std::array<double, 7> level_two{};
  level_two[0] = level_one[0];
  level_two[1] = level_one[0]; // Literal printed Eq. 16.
  level_two[2] = k1 * level_one[2] + (1.0 - k1) * level_one[5];
  level_two[3] = k1 * level_one[3] + (1.0 - k1) * level_one[4];
  level_two[4] = k1 * level_one[4] + (1.0 - k1) * level_one[2];
  level_two[5] = k1 * level_one[5] + (1.0 - k1) * level_one[3];
  level_two[6] = level_one[6];

  std::array<double, 7> level_three{};
  level_three[0] = k2 * level_two[0] + (1.0 - k2) * level_two[5];
  level_three[1] = k2 * level_two[1] + (1.0 - k2) * level_two[4];
  level_three[2] = level_two[2];
  level_three[3] = level_two[3];
  level_three[4] = k2 * level_two[4] + (1.0 - k2) * level_two[0];
  level_three[5] = k2 * level_two[5] + (1.0 - k2) * level_two[1];
  level_three[6] = level_two[6];

  std::array<double, 7> level_zero{};
  level_zero[0] = k3 * level_three[0] + (1.0 - k3) * level_three[2];
  level_zero[1] = k3 * level_three[1] + (1.0 - k3) * level_three[3];
  level_zero[2] = k3 * level_three[2] + (1.0 - k3) * level_three[1];
  level_zero[3] = k3 * level_three[3] + (1.0 - k3) * level_three[0];
  level_zero[4] = level_three[4];
  level_zero[5] = level_three[5];
  level_zero[6] = level_three[6];

  LocalCoefficientRow result{};
  double neighbour_sum = 0.0;
  for (std::size_t index = 0; index < result.neighbour.size(); ++index) {
    result.neighbour[index] = k0 * level_zero[index];
    neighbour_sum += level_zero[index];
  }
  result.diagonal = (1.0 - k0) * neighbour_sum + level_zero[6];
  result.source = input.source;
  require_finite(result);
  return result;
}

double occurrence_value(const LocalCoefficientRow &row,
                        std::uint32_t occurrence) {
  if (occurrence < 6U) {
    return row.neighbour[occurrence];
  }
  if (occurrence == 6U) {
    return row.diagonal;
  }
  if (occurrence == 7U) {
    return row.source;
  }
  fail("algebraic occurrence is out of range");
}

std::uint64_t evaluator_fingerprint(const ReplacementGroup &group) {
  std::uint64_t hash = kFnvOffset;
  hash_u64(hash, kAlgorithmFingerprint);
  hash_u64(hash, 0x4556414c55415445ULL);
  hash_u64(hash, static_cast<std::uint64_t>(group.links.size()));
  for (const ImmersedLinkId link : group.links) {
    hash_u64(hash, link);
  }
  hash_u64(hash,
           static_cast<std::uint64_t>(group.algebraic_occurrences.size()));
  for (const std::uint32_t occurrence : group.algebraic_occurrences) {
    hash_u32(hash, occurrence);
  }
  return hash;
}

std::uint64_t stable_group_id(const ReplacementGroup &group) {
  std::uint64_t hash = kFnvOffset;
  hash_u64(hash, 0x47524f55505f4944ULL);
  hash_u64(hash, evaluator_fingerprint(group));
  return hash;
}

void canonicalize_group(ReplacementGroup &group) {
  std::sort(group.links.begin(), group.links.end());
  std::sort(group.algebraic_occurrences.begin(),
            group.algebraic_occurrences.end());
}

bool group_less(const ReplacementGroup &lhs, const ReplacementGroup &rhs) {
  if (lhs.stable_group_id != rhs.stable_group_id) {
    return lhs.stable_group_id < rhs.stable_group_id;
  }
  if (lhs.links != rhs.links) {
    return lhs.links < rhs.links;
  }
  if (lhs.algebraic_occurrences != rhs.algebraic_occurrences) {
    return lhs.algebraic_occurrences < rhs.algebraic_occurrences;
  }
  return lhs.evaluator_fingerprint < rhs.evaluator_fingerprint;
}

std::uint64_t plan_fingerprint(mesh::GlobalCellId cell,
                               std::vector<ReplacementGroup> groups) {
  for (auto &group : groups) {
    canonicalize_group(group);
  }
  std::sort(groups.begin(), groups.end(), group_less);
  std::uint64_t hash = kFnvOffset;
  hash_u64(hash, kAlgorithmFingerprint);
  hash_u64(hash, cell);
  hash_u64(hash, static_cast<std::uint64_t>(groups.size()));
  for (const auto &group : groups) {
    hash_u64(hash, group.stable_group_id);
    hash_u64(hash, group.evaluator_fingerprint);
    hash_u64(hash, static_cast<std::uint64_t>(group.links.size()));
    for (const auto link : group.links) {
      hash_u64(hash, link);
    }
    hash_u64(hash,
             static_cast<std::uint64_t>(group.algebraic_occurrences.size()));
    for (const auto occurrence : group.algebraic_occurrences) {
      hash_u32(hash, occurrence);
    }
  }
  return hash;
}

ReplacementGroup make_group(std::vector<ImmersedLinkId> links) {
  ReplacementGroup group{};
  group.links = std::move(links);
  canonicalize_group(group);
  return group;
}

void finalize_group(ReplacementGroup &group) {
  canonicalize_group(group);
  group.evaluator_fingerprint = evaluator_fingerprint(group);
  group.stable_group_id = stable_group_id(group);
}

} // namespace

LocalCoefficientRow LocalFlowPatternTransform::transform_full(
    const LocalCoefficientRow &row, double normal_scale,
    runtime::Real3 solid_to_fluid_normal) const {
  if (!std::isfinite(normal_scale) || normal_scale < 0.0 ||
      normal_scale > 1.0) {
    fail("normal scale must be finite and in [0,1]");
  }
  const Real3 unit_normal = normalized(solid_to_fluid_normal);
  const Real3 positive_normal{std::abs(unit_normal.x), std::abs(unit_normal.y),
                              std::abs(unit_normal.z)};
  LocalCoefficientRow oriented = row;
  swap_direction_pairs(oriented, unit_normal);
  const auto factors = rotation_factors(positive_normal);
  auto result = paper_transform_factors(oriented, normal_scale, factors[0],
                                        factors[1], factors[2]);
  swap_direction_pairs(result, unit_normal);
  return result;
}

RowReplacementPlan LocalFlowPatternTransform::plan_row(
    mesh::GlobalCellId active_cell, const std::vector<ImmersedLinkId> &links,
    const LocalCoefficientRow &immersed_projection) const {
  require_finite(immersed_projection);
  std::vector<ImmersedLinkId> canonical_links = links;
  std::sort(canonical_links.begin(), canonical_links.end());
  if (std::adjacent_find(canonical_links.begin(), canonical_links.end()) !=
      canonical_links.end()) {
    fail("duplicate immersed link");
  }

  std::vector<std::uint32_t> nonzero_neighbours;
  for (std::uint32_t occurrence = 0U; occurrence < 6U; ++occurrence) {
    if (occurrence_value(immersed_projection, occurrence) != 0.0) {
      nonzero_neighbours.push_back(occurrence);
    }
  }
  const bool has_diagonal = immersed_projection.diagonal != 0.0;
  const bool has_source = immersed_projection.source != 0.0;
  const bool has_occurrence =
      !nonzero_neighbours.empty() || has_diagonal || has_source;
  if (has_occurrence && canonical_links.empty()) {
    fail("nonempty immersed row requires at least one link");
  }
  if (!has_occurrence && !canonical_links.empty()) {
    fail("every immersed link must be covered by an occurrence");
  }

  RowReplacementPlan plan{};
  plan.active_cell = active_cell;
  if (!has_occurrence) {
    plan.fingerprint = plan_fingerprint(active_cell, {});
    return plan;
  }

  std::vector<ReplacementGroup> singleton_groups;
  singleton_groups.reserve(canonical_links.size());
  for (const auto link : canonical_links) {
    singleton_groups.push_back(make_group({link}));
  }
  for (std::size_t index = 0; index < nonzero_neighbours.size(); ++index) {
    singleton_groups[index % singleton_groups.size()]
        .algebraic_occurrences.push_back(nonzero_neighbours[index]);
  }

  std::vector<ReplacementGroup> joint_groups;
  const auto add_local_occurrence = [&](std::uint32_t occurrence) {
    if (canonical_links.size() == 1U) {
      singleton_groups.front().algebraic_occurrences.push_back(occurrence);
    } else {
      auto group = make_group(canonical_links);
      group.algebraic_occurrences.push_back(occurrence);
      joint_groups.push_back(std::move(group));
    }
  };
  if (has_diagonal) {
    add_local_occurrence(6U);
  }
  if (has_source) {
    add_local_occurrence(7U);
  }

  for (auto &group : singleton_groups) {
    if (!group.algebraic_occurrences.empty()) {
      finalize_group(group);
      plan.groups.push_back(std::move(group));
    }
  }
  for (auto &group : joint_groups) {
    finalize_group(group);
    plan.groups.push_back(std::move(group));
  }
  std::sort(plan.groups.begin(), plan.groups.end(), group_less);

  std::vector<ImmersedLinkId> covered_links;
  for (const auto &group : plan.groups) {
    covered_links.insert(covered_links.end(), group.links.begin(),
                         group.links.end());
  }
  std::sort(covered_links.begin(), covered_links.end());
  covered_links.erase(std::unique(covered_links.begin(), covered_links.end()),
                      covered_links.end());
  if (covered_links != canonical_links) {
    fail("every immersed link must be covered by a replacement group");
  }
  plan.fingerprint = plan_fingerprint(active_cell, plan.groups);
  return plan;
}

double LocalFlowPatternTransform::evaluate_wall_replacement(
    const RowReplacementPlan &input_plan,
    const LocalCoefficientRow &immutable_snapshot,
    const std::vector<double> &link_local_symbols) const {
  require_finite(immutable_snapshot);
  for (const double symbol : link_local_symbols) {
    if (!std::isfinite(symbol)) {
      fail("link-local symbols must be finite");
    }
  }

  auto groups = input_plan.groups;
  std::vector<ImmersedLinkId> canonical_links;
  std::array<bool, 8> occurrence_seen{};
  for (auto &group : groups) {
    canonicalize_group(group);
    if (group.links.empty() || group.algebraic_occurrences.empty()) {
      fail("replacement group must be nonempty");
    }
    if (std::adjacent_find(group.links.begin(), group.links.end()) !=
        group.links.end()) {
      fail("replacement group contains a duplicate link");
    }
    if (std::adjacent_find(group.algebraic_occurrences.begin(),
                           group.algebraic_occurrences.end()) !=
        group.algebraic_occurrences.end()) {
      fail("replacement group contains a duplicate occurrence");
    }
    for (const auto occurrence : group.algebraic_occurrences) {
      if (occurrence >= occurrence_seen.size()) {
        fail("algebraic occurrence is out of range");
      }
      if (occurrence_seen[occurrence]) {
        fail("algebraic occurrence belongs to more than one group");
      }
      occurrence_seen[occurrence] = true;
    }
    if (group.evaluator_fingerprint != evaluator_fingerprint(group) ||
        group.stable_group_id != stable_group_id(group)) {
      fail("replacement group fingerprint is inconsistent");
    }
    canonical_links.insert(canonical_links.end(), group.links.begin(),
                           group.links.end());
  }
  std::sort(canonical_links.begin(), canonical_links.end());
  canonical_links.erase(
      std::unique(canonical_links.begin(), canonical_links.end()),
      canonical_links.end());
  if (canonical_links.size() != link_local_symbols.size()) {
    fail("link-local symbol layout is inconsistent");
  }
  if (input_plan.fingerprint !=
      plan_fingerprint(input_plan.active_cell, groups)) {
    fail("replacement plan fingerprint is inconsistent");
  }
  std::sort(groups.begin(), groups.end(), group_less);

  double result = 0.0;
  for (const auto &group : groups) {
    double symbol_sum = 0.0;
    for (const auto link : group.links) {
      const auto found = std::lower_bound(canonical_links.begin(),
                                          canonical_links.end(), link);
      if (found == canonical_links.end() || *found != link) {
        fail("replacement group link is not in the canonical layout");
      }
      const auto index =
          static_cast<std::size_t>(found - canonical_links.begin());
      symbol_sum += link_local_symbols[index];
    }
    const double group_symbol =
        symbol_sum / static_cast<double>(group.links.size());
    for (const auto occurrence : group.algebraic_occurrences) {
      result += occurrence_value(immutable_snapshot, occurrence) * group_symbol;
    }
  }
  if (!std::isfinite(result)) {
    fail("replacement evaluation produced a non-finite value");
  }
  return result;
}

std::uint64_t
LocalFlowPatternTransform::algorithm_fingerprint() const noexcept {
  return kAlgorithmFingerprint;
}

namespace detail {
#if defined(HUNDUN_IMMERSED_ENABLE_TEST_ACCESS)
LocalCoefficientRow
paper_transform_factors_for_test(const LocalCoefficientRow &row, double k0,
                                 double k1, double k2, double k3) {
  return paper_transform_factors(row, k0, k1, k2, k3);
}
#endif
} // namespace detail

} // namespace hundun::immersed
