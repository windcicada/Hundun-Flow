// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/fvm_immersed_operator.hpp"

#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "fvm_immersed_boundary_authority_detail.hpp"
#include "ib_deterministic_qr_detail.hpp"
#include "ib_periodic_surface_window_detail.hpp"
#include "ib_quadratic_reconstruction_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace hundun::finite_volume {
namespace {

using immersed::CellRegion;
using mesh::EntityOwnership;
using mesh::LocalCellId;
using mesh::LocalFaceId;
using runtime::Error;
using runtime::Int3;
using runtime::Real3;

constexpr std::size_t kInvalidFace = std::numeric_limits<std::size_t>::max();
constexpr std::uint64_t kTermHashOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kTermHashPrime = UINT64_C(1099511628211);

thread_local const std::vector<detail::ImmersedWallNormalGradient>
    *scoped_wall_normal_gradients = nullptr;
thread_local bool scoped_force_authority_evaluation = false;

class ScopedWallNormalGradients final {
public:
  explicit ScopedWallNormalGradients(
      const std::vector<detail::ImmersedWallNormalGradient> &gradients) {
    if (scoped_wall_normal_gradients != nullptr)
      throw Error("immersed operator wall authority scope is already active");
    scoped_wall_normal_gradients = &gradients;
  }
  ~ScopedWallNormalGradients() { scoped_wall_normal_gradients = nullptr; }

  ScopedWallNormalGradients(const ScopedWallNormalGradients &) = delete;
  ScopedWallNormalGradients &
  operator=(const ScopedWallNormalGradients &) = delete;
};

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (unsigned byte = 0U; byte < 8U; ++byte) {
    hash ^= (value >> (byte * 8U)) & UINT64_C(0xff);
    hash *= kTermHashPrime;
  }
}

bool same(Int3 lhs, Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 subtract(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 add(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 multiply(double scale, Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Real3 cross(Real3 lhs, Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y,
          lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double norm(Real3 value) {
  const double result = std::sqrt(dot(value, value));
  if (!(result > 0.0) || !std::isfinite(result))
    throw Error("immersed operator geometric distance is invalid");
  return result;
}

struct OperationGuard final {
  explicit OperationGuard(bool &active) : active_(&active) {
    if (active)
      throw Error("immersed operator operation is already active");
    active = true;
  }
  ~OperationGuard() { *active_ = false; }
  bool *active_;
};

std::size_t direction_occurrence(Int3 fluid, Int3 solid) {
  const Int3 delta{solid.x - fluid.x, solid.y - fluid.y, solid.z - fluid.z};
  if (delta.x == -1 && delta.y == 0 && delta.z == 0)
    return 0U;
  if (delta.x == 1 && delta.y == 0 && delta.z == 0)
    return 1U;
  if (delta.x == 0 && delta.y == -1 && delta.z == 0)
    return 2U;
  if (delta.x == 0 && delta.y == 1 && delta.z == 0)
    return 3U;
  if (delta.x == 0 && delta.y == 0 && delta.z == -1)
    return 4U;
  if (delta.x == 0 && delta.y == 0 && delta.z == 1)
    return 5U;
  throw Error("immersed operator link is not a direct logical neighbour");
}

Int3 occurrence_offset(std::size_t occurrence) {
  constexpr std::array<Int3, 6> offsets{Int3{-1, 0, 0}, Int3{1, 0, 0},
                                        Int3{0, -1, 0}, Int3{0, 1, 0},
                                        Int3{0, 0, -1}, Int3{0, 0, 1}};
  if (occurrence >= offsets.size())
    throw Error("immersed operator logical occurrence is out of range");
  return offsets[occurrence];
}

struct CellPairFace final {
  mesh::GlobalCellId first{};
  mesh::GlobalCellId second{};
  LocalFaceId face{kInvalidFace};
};

struct SharedFace final {
  LocalFaceId face{kInvalidFace};
  LocalCellId owner{};
  LocalCellId neighbour{};
  bool owner_owned{};
  bool neighbour_owned{};
  Real3 area{};
  Real3 displacement{};
  Real3 owner_center{};
  Real3 face_center{};
};

struct WallLink final {
  immersed::ImmersedLinkId id{};
  LocalCellId fluid{};
  LocalFaceId face{kInvalidFace};
  std::size_t occurrence{};
  Real3 area_from_fluid{};
  double signed_wall_measure_m2{};
  Real3 wall_intercept_m{};
  double normal_scale{};
  Real3 solid_to_fluid_normal{};
  Real3 pressure_quadrature_m{};
  Real3 surface_measure_m2{};
  Real3 surface_patch_centroid_m{};
  immersed::LocalCoefficientRow background_row;
  immersed::LocalCoefficientRow transformed_row;
  std::array<Real3, 7> local_sample_points_m{};
};

enum class PhysicalTermKind : std::uint8_t {
  convective_direct,
  pressure_direct,
  viscous_orthogonal,
  viscous_deferred_gradient
};

struct PhysicalTerm final {
  std::uint64_t stable_id{};
  immersed::ImmersedLinkId link{};
  PhysicalTermKind kind{PhysicalTermKind::convective_direct};
  std::uint32_t algebraic_occurrence{};
  std::uint32_t output_component{};
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

enum class BoundaryReplacementKind : std::uint8_t {
  pressure_face,
  pressure_diagonal_defect,
  pressure_neighbour_defect,
  viscous_wall,
  viscous_diagonal_defect,
  viscous_neighbour_defect
};

struct BoundaryReplacementTerm final {
  std::uint64_t stable_id{};
  immersed::ImmersedLinkId link{};
  BoundaryReplacementKind kind{BoundaryReplacementKind::pressure_face};
  std::uint32_t occurrence{};
  std::uint32_t component{};
  double scale{};
};

struct BoundaryEvaluationGroup final {
  std::uint64_t stable_id{};
  std::vector<immersed::ImmersedLinkId> links;
  std::vector<std::size_t> physical_term_indices;
  std::vector<std::size_t> replacement_term_indices;
};

enum class AffineInputKind : std::uint8_t { pressure, velocity };

struct AffineDonorKey final {
  AffineInputKind input_kind{AffineInputKind::pressure};
  mesh::GlobalCellId donor{};
  std::uint32_t input_component{};
  std::uint32_t output_component{};

  friend bool operator<(const AffineDonorKey &lhs,
                        const AffineDonorKey &rhs) noexcept {
    return std::tie(lhs.input_kind, lhs.donor, lhs.input_component,
                    lhs.output_component) <
           std::tie(rhs.input_kind, rhs.donor, rhs.input_component,
                    rhs.output_component);
  }
};

struct AffineInputValueKey final {
  AffineInputKind input_kind{AffineInputKind::pressure};
  mesh::GlobalCellId donor{};
  std::uint32_t input_component{};

  friend bool operator<(const AffineInputValueKey &lhs,
                        const AffineInputValueKey &rhs) noexcept {
    return std::tie(lhs.input_kind, lhs.donor, lhs.input_component) <
           std::tie(rhs.input_kind, rhs.donor, rhs.input_component);
  }
};

struct AffineLinkCoefficient final {
  immersed::ImmersedLinkId link{};
  std::size_t wall_link_index{};
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

struct AffineDonorTerm final {
  std::uint64_t stable_id{};
  AffineDonorKey key;
  std::size_t input_snapshot_index{};
  std::vector<AffineLinkCoefficient> link_coefficients;
};

struct AffineWallGradientTerm final {
  std::uint64_t stable_id{};
  immersed::ImmersedLinkId link{};
  std::uint32_t output_component{};
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

struct PressureAffinePlan final {
  std::vector<AffineDonorTerm> donor_terms;
  std::vector<AffineWallGradientTerm> wall_gradient_terms;
};

struct BoundaryRowPlan final {
  LocalCellId cell{};
  std::vector<std::size_t> wall_links;
  std::optional<immersed::QuadraticReconstruction> row_reconstruction;
  std::vector<AffineInputValueKey> affine_input_keys;
  std::vector<PhysicalTerm> physical_terms;
  std::vector<BoundaryReplacementTerm> replacement_terms;
  std::vector<BoundaryEvaluationGroup> evaluation_groups;
  PressureAffinePlan background_pressure_unconstrained;
  PressureAffinePlan background_pressure_constrained;
  std::vector<AffineDonorTerm> background_viscous_terms;
  PressureAffinePlan pressure_unconstrained;
  PressureAffinePlan pressure_constrained;
  PressureAffinePlan pressure_constrained_force;
  std::vector<AffineDonorTerm> viscous_terms;
  std::uint64_t affine_plan_fingerprint{};
  std::uint64_t fingerprint{};
};

struct BoundaryRowEvaluation final {
  ImmersedResidualParts residual;
  ImmersedResidualParts background;
  ImmersedResidualParts removed_background;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  struct EvaluatedAffineDonorTerm final {
    const AffineDonorTerm *descriptor{};
    double effective_coefficient{};
  };
  std::vector<EvaluatedAffineDonorTerm> affine_donor_terms;
  std::vector<EvaluatedAffineDonorTerm> background_affine_donor_terms;
  const std::vector<AffineWallGradientTerm> *affine_wall_gradient_terms{};
  const std::vector<AffineWallGradientTerm>
      *background_affine_wall_gradient_terms{};
#endif
  std::uint64_t affine_plan_fingerprint{};
  std::uint64_t evaluated_group_count{};
  std::uint64_t simultaneous_substitution_count{};
  std::uint64_t canonical_affine_row_evaluation_count{};
  std::uint64_t link_local_runtime_evaluation_count{};
  std::uint64_t immutable_input_snapshot_count{};
  std::uint64_t background_functional_evaluation_count{};
  std::uint64_t background_removal_count{};
  std::uint64_t final_row_write_count{};
};

std::uint64_t physical_term_id(mesh::GlobalCellId cell,
                               immersed::ImmersedLinkId link,
                               PhysicalTermKind kind,
                               std::uint32_t occurrence,
                               std::uint32_t component) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e4454524d31));
  hash_u64(hash, cell);
  hash_u64(hash, link);
  hash_u64(hash, static_cast<std::uint64_t>(kind));
  hash_u64(hash, occurrence);
  hash_u64(hash, component);
  return hash;
}

std::uint64_t replacement_term_id(mesh::GlobalCellId cell,
                                  immersed::ImmersedLinkId link,
                                  BoundaryReplacementKind kind,
                                  std::uint32_t occurrence,
                                  std::uint32_t component) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e4457525031));
  hash_u64(hash, cell);
  hash_u64(hash, link);
  hash_u64(hash, static_cast<std::uint64_t>(kind));
  hash_u64(hash, occurrence);
  hash_u64(hash, component);
  return hash;
}

void hash_double(std::uint64_t &hash, double value) noexcept;

std::uint64_t row_fingerprint(
    mesh::GlobalCellId cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links, const std::vector<PhysicalTerm> &terms,
    const std::vector<BoundaryReplacementTerm> &replacement_terms,
    const std::vector<BoundaryEvaluationGroup> &evaluation_groups,
    std::uint64_t affine_fingerprint) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e44524f5731));
  hash_u64(hash, cell);
  hash_u64(hash, link_indices.size());
  for (const auto index : link_indices)
    hash_u64(hash, links[index].id);
  hash_u64(hash, terms.size());
  for (const auto &term : terms) {
    hash_u64(hash, term.stable_id);
    hash_double(hash, term.coefficient);
    hash_u64(hash, term.source_term_ids.size());
    for (const auto source : term.source_term_ids)
      hash_u64(hash, source);
  }
  hash_u64(hash, replacement_terms.size());
  for (const auto &term : replacement_terms) {
    hash_u64(hash, term.stable_id);
    hash_double(hash, term.scale);
  }
  hash_u64(hash, evaluation_groups.size());
  for (const auto &group : evaluation_groups) {
    hash_u64(hash, group.stable_id);
    hash_u64(hash, group.links.size());
    for (const auto link : group.links)
      hash_u64(hash, link);
    hash_u64(hash, group.physical_term_indices.size());
    for (const auto index : group.physical_term_indices)
      hash_u64(hash, terms[index].stable_id);
    hash_u64(hash, group.replacement_term_indices.size());
    for (const auto index : group.replacement_term_indices)
      hash_u64(hash, replacement_terms[index].stable_id);
  }
  hash_u64(hash, affine_fingerprint);
  return hash;
}

std::uint64_t
boundary_group_id(
    mesh::GlobalCellId cell,
    const std::vector<immersed::ImmersedLinkId> &links,
    const std::vector<PhysicalTerm> &physical_terms,
    const std::vector<std::size_t> &physical_term_indices,
    const std::vector<BoundaryReplacementTerm> &replacement_terms,
    const std::vector<std::size_t> &replacement_term_indices) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e444a4f494e));
  hash_u64(hash, cell);
  hash_u64(hash, links.size());
  for (const auto link : links)
    hash_u64(hash, link);
  hash_u64(hash, physical_term_indices.size());
  for (const auto index : physical_term_indices)
    hash_u64(hash, physical_terms[index].stable_id);
  hash_u64(hash, replacement_term_indices.size());
  for (const auto index : replacement_term_indices)
    hash_u64(hash, replacement_terms[index].stable_id);
  return hash == 0U ? UINT64_C(1) : hash;
}

void validate_boundary_evaluation_groups(
    const std::vector<PhysicalTerm> &physical_terms,
    const std::vector<BoundaryReplacementTerm> &replacement_terms,
    const std::vector<BoundaryEvaluationGroup> &groups) {
  std::vector<std::uint8_t> physical_covered(physical_terms.size(), 0U);
  std::vector<std::uint8_t> replacement_covered(replacement_terms.size(), 0U);
  for (const auto &group : groups) {
    if (group.stable_id == 0U || group.links.empty() ||
        group.physical_term_indices.empty() ||
        group.replacement_term_indices.empty())
      throw Error("immersed operator boundary evaluation group is invalid");
    if (!std::is_sorted(group.links.begin(), group.links.end()) ||
        std::adjacent_find(group.links.begin(), group.links.end()) !=
            group.links.end())
      throw Error("immersed operator boundary group links are invalid");
    for (const auto term_index : group.physical_term_indices) {
      if (term_index >= physical_terms.size() ||
          physical_covered[term_index] != 0U)
        throw Error("immersed operator physical term coverage is invalid");
      const auto &term = physical_terms[term_index];
      if (!std::isfinite(term.coefficient))
        throw Error("immersed operator physical term coefficient is non-finite");
      if (!std::binary_search(group.links.begin(), group.links.end(),
                              term.link))
        throw Error("immersed operator physical term link is inconsistent");
      physical_covered[term_index] = 1U;
    }
    for (const auto term_index : group.replacement_term_indices) {
      if (term_index >= replacement_terms.size() ||
          replacement_covered[term_index] != 0U)
        throw Error("immersed operator boundary term coverage is invalid");
      if (!std::isfinite(replacement_terms[term_index].scale))
        throw Error("immersed operator boundary term scale is non-finite");
      if (!std::binary_search(group.links.begin(), group.links.end(),
                              replacement_terms[term_index].link))
        throw Error("immersed operator boundary group link is inconsistent");
      replacement_covered[term_index] = 1U;
    }
  }
  if (std::find(physical_covered.begin(), physical_covered.end(), 0U) !=
      physical_covered.end())
    throw Error("immersed operator physical term coverage is incomplete");
  if (std::find(replacement_covered.begin(), replacement_covered.end(), 0U) !=
      replacement_covered.end())
    throw Error("immersed operator boundary term coverage is incomplete");
}

void hash_double(std::uint64_t &hash, double value) noexcept {
  std::uint64_t encoded{};
  std::memcpy(&encoded, &value, sizeof(encoded));
  hash_u64(hash, encoded);
}

std::uint64_t affine_donor_term_id(mesh::GlobalCellId cell,
                                   const AffineDonorKey &key) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e4441464431));
  hash_u64(hash, cell);
  hash_u64(hash, static_cast<std::uint64_t>(key.input_kind));
  hash_u64(hash, key.donor);
  hash_u64(hash, key.input_component);
  hash_u64(hash, key.output_component);
  return hash == 0U ? UINT64_C(1) : hash;
}

std::uint64_t affine_wall_gradient_term_id(
    mesh::GlobalCellId cell, immersed::ImmersedLinkId link,
    std::uint32_t output_component) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e4441465731));
  hash_u64(hash, cell);
  hash_u64(hash, link);
  hash_u64(hash, output_component);
  return hash == 0U ? UINT64_C(1) : hash;
}

struct LinkCoefficientAccumulator final {
  std::size_t wall_link_index{};
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

using AffineDonorAccumulator =
    std::map<AffineDonorKey,
             std::map<immersed::ImmersedLinkId, LinkCoefficientAccumulator>>;
struct WallGradientAccumulatorValue final {
  double coefficient{};
  std::vector<std::uint64_t> source_term_ids;
};

using WallGradientAccumulator =
    std::map<std::pair<immersed::ImmersedLinkId, std::uint32_t>,
             WallGradientAccumulatorValue>;

void add_affine_weights(AffineDonorAccumulator &target,
                        AffineInputKind input_kind,
                        std::uint32_t input_component,
                        std::uint32_t output_component,
                        immersed::ImmersedLinkId link,
                        std::size_t wall_link_index, double scale,
                        std::uint64_t source_term_id,
                        const std::vector<immersed::WeightedDonor> &weights) {
  if (!std::isfinite(scale) || source_term_id == 0U)
    throw Error("immersed operator affine functional scale is non-finite");
  for (const auto &donor : weights) {
    const double contribution = scale * donor.weight;
    if (!std::isfinite(contribution))
      throw Error("immersed operator affine donor coefficient is non-finite");
    auto &by_link = target[{input_kind, donor.global_cell, input_component,
                            output_component}];
    auto [position, inserted] = by_link.try_emplace(
        link, LinkCoefficientAccumulator{wall_link_index, 0.0, {}});
    if (!inserted && position->second.wall_link_index != wall_link_index)
      throw Error("immersed operator affine link identity is inconsistent");
    position->second.coefficient += contribution;
    position->second.source_term_ids.push_back(source_term_id);
    if (!std::isfinite(position->second.coefficient))
      throw Error("immersed operator affine donor coefficient is non-finite");
  }
}

void add_wall_gradient_weight(WallGradientAccumulator &target,
                              immersed::ImmersedLinkId link,
                              std::uint32_t output_component, double value,
                              std::uint64_t source_term_id) {
  if (!std::isfinite(value) || source_term_id == 0U)
    throw Error("immersed operator affine wall coefficient is non-finite");
  auto &entry = target[{link, output_component}];
  entry.coefficient += value;
  entry.source_term_ids.push_back(source_term_id);
  if (!std::isfinite(entry.coefficient))
    throw Error("immersed operator affine wall coefficient is non-finite");
}

std::vector<AffineDonorTerm>
finalize_affine_donor_terms(mesh::GlobalCellId cell,
                            const AffineDonorAccumulator &accumulator) {
  std::vector<AffineDonorTerm> result;
  result.reserve(accumulator.size());
  for (const auto &[key, by_link] : accumulator) {
    AffineDonorTerm term{};
    term.stable_id = affine_donor_term_id(cell, key);
    term.key = key;
    for (const auto &[link, value] : by_link) {
      if (value.coefficient == 0.0)
        continue;
      auto source_term_ids = value.source_term_ids;
      std::sort(source_term_ids.begin(), source_term_ids.end());
      source_term_ids.erase(
          std::unique(source_term_ids.begin(), source_term_ids.end()),
          source_term_ids.end());
      term.link_coefficients.push_back({link, value.wall_link_index,
                                        value.coefficient,
                                        std::move(source_term_ids)});
    }
    if (!term.link_coefficients.empty())
      result.push_back(std::move(term));
  }
  return result;
}

std::vector<AffineWallGradientTerm>
finalize_wall_gradient_terms(mesh::GlobalCellId cell,
                             const WallGradientAccumulator &accumulator) {
  std::vector<AffineWallGradientTerm> result;
  result.reserve(accumulator.size());
  for (const auto &[key, value] : accumulator) {
    if (value.coefficient == 0.0)
      continue;
    auto source_term_ids = value.source_term_ids;
    std::sort(source_term_ids.begin(), source_term_ids.end());
    source_term_ids.erase(
        std::unique(source_term_ids.begin(), source_term_ids.end()),
        source_term_ids.end());
    result.push_back({affine_wall_gradient_term_id(cell, key.first, key.second),
                      key.first, key.second, value.coefficient,
                      std::move(source_term_ids)});
  }
  return result;
}

std::uint64_t affine_plan_fingerprint(
    mesh::GlobalCellId cell,
    const PressureAffinePlan &background_unconstrained,
    const PressureAffinePlan &background_constrained,
    const std::vector<AffineDonorTerm> &background_viscous_terms,
    const PressureAffinePlan &unconstrained,
    const PressureAffinePlan &constrained,
    const PressureAffinePlan &constrained_force,
    const std::vector<AffineDonorTerm> &viscous_terms) noexcept {
  std::uint64_t hash = kTermHashOffset;
  hash_u64(hash, UINT64_C(0x48554e4441465031));
  hash_u64(hash, cell);
  const auto hash_donor_terms = [&](const std::vector<AffineDonorTerm> &terms) {
    hash_u64(hash, terms.size());
    for (const auto &term : terms) {
      hash_u64(hash, term.stable_id);
      hash_u64(hash, term.input_snapshot_index);
      hash_u64(hash, term.link_coefficients.size());
      for (const auto &coefficient : term.link_coefficients) {
        hash_u64(hash, coefficient.link);
        hash_double(hash, coefficient.coefficient);
        hash_u64(hash, coefficient.source_term_ids.size());
        for (const auto source : coefficient.source_term_ids)
          hash_u64(hash, source);
      }
    }
  };
  const auto hash_wall_terms =
      [&](const std::vector<AffineWallGradientTerm> &terms) {
        hash_u64(hash, terms.size());
        for (const auto &term : terms) {
          hash_u64(hash, term.stable_id);
          hash_double(hash, term.coefficient);
          hash_u64(hash, term.source_term_ids.size());
          for (const auto source : term.source_term_ids)
            hash_u64(hash, source);
        }
      };
  hash_donor_terms(background_unconstrained.donor_terms);
  hash_wall_terms(background_unconstrained.wall_gradient_terms);
  hash_donor_terms(background_constrained.donor_terms);
  hash_wall_terms(background_constrained.wall_gradient_terms);
  hash_donor_terms(background_viscous_terms);
  hash_donor_terms(unconstrained.donor_terms);
  hash_wall_terms(unconstrained.wall_gradient_terms);
  hash_donor_terms(constrained.donor_terms);
  hash_wall_terms(constrained.wall_gradient_terms);
  hash_donor_terms(constrained_force.donor_terms);
  hash_wall_terms(constrained_force.wall_gradient_terms);
  hash_donor_terms(viscous_terms);
  return hash == 0U ? UINT64_C(1) : hash;
}

bool bits_zero(double value) noexcept {
  std::uint64_t encoded{};
  std::memcpy(&encoded, &value, sizeof(encoded));
  return encoded == 0U;
}

} // namespace

struct ImmersedOperatorAdapter::Impl final {
  const mesh::MeshTopology *topology{};
  const immersed::GhostStencilPlan *ghost_plan{};
  const runtime::MpiContext *mpi{};
  runtime::Box3 owned_box{};
  Int3 local_extent{};
  Int3 global_extent{};
  std::vector<std::uint8_t> local_active;
  std::vector<SharedFace> shared_faces;
  std::vector<WallLink> wall_links;
  std::vector<BoundaryRowPlan> active_rows;
  std::vector<std::size_t> boundary_row_by_wall_link;
  mutable std::vector<double> scratch;
  mutable std::vector<double> affine_input_values;
  mutable std::vector<double> viscosity_by_wall_link;
  mutable std::vector<double> pressure_wall_normal_gradient_by_link;
  std::vector<std::uint8_t> local_pressure_wall_link;
  mutable ImmersedOperatorReport last_report{};
  mutable std::uint64_t last_wall_functional_evaluation_count{};
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  mutable std::vector<detail::ImmersedBoundaryRowEvaluationRecord>
      last_boundary_row_evaluations;
  mutable std::uint64_t last_boundary_authority_lookup_probe_count{};
#endif
  mutable bool active{};

  Int3 field_index(LocalCellId cell, int ghost_width) const {
    return field_index_global(topology->global_cell_id(cell), ghost_width);
  }

  Int3 field_index_global(mesh::GlobalCellId global_id,
                          int ghost_width) const {
    const auto nx = static_cast<std::uint64_t>(global_extent.x);
    const auto ny = static_cast<std::uint64_t>(global_extent.y);
    const auto nz = static_cast<std::uint64_t>(global_extent.z);
    if (nx == 0U || ny == 0U || nz == 0U ||
        nx > std::numeric_limits<std::uint64_t>::max() / ny ||
        nx * ny > std::numeric_limits<std::uint64_t>::max() / nz ||
        global_id >= nx * ny * nz)
      throw Error("immersed operator global donor is out of range");
    const auto plane = nx * ny;
    const Int3 global{static_cast<int>(global_id % nx),
                      static_cast<int>((global_id / nx) % ny),
                      static_cast<int>(global_id / plane)};
    int coordinates[3]{global.x, global.y, global.z};
    const int begin[3]{owned_box.begin.x, owned_box.begin.y, owned_box.begin.z};
    const int extent[3]{local_extent.x, local_extent.y, local_extent.z};
    const int global_size[3]{global_extent.x, global_extent.y, global_extent.z};
    int result[3]{};
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      result[axis] = coordinates[axis] - begin[axis];
      const int lower = result[axis] - global_size[axis];
      const int upper = result[axis] + global_size[axis];
      const auto inside = [&](int value) {
        return value >= -ghost_width && value < extent[axis] + ghost_width;
      };
      if (!inside(result[axis]) && inside(lower))
        result[axis] = lower;
      else if (!inside(result[axis]) && inside(upper))
        result[axis] = upper;
      if (!inside(result[axis]))
        throw Error("immersed operator global donor exceeds field halo");
    }
    return {result[0], result[1], result[2]};
  }

  double read_cell(const runtime::FieldView<const double> &view,
                   LocalCellId cell, int component) const {
    if (cell >= local_active.size() || local_active[cell] == 0U)
      throw Error("immersed operator attempted an inactive cell read");
    const auto index = field_index(cell, view.ghost_width());
    const double value = view(index.x, index.y, index.z, component);
    if (!std::isfinite(value))
      throw Error("immersed operator active field value is non-finite");
    return value;
  }

  double read_global(const runtime::FieldView<const double> &view,
                     mesh::GlobalCellId global_id, int component) const {
    const auto index = field_index_global(global_id, view.ghost_width());
    const double value = view(index.x, index.y, index.z, component);
    if (!std::isfinite(value))
      throw Error("immersed operator affine donor value is non-finite");
    return value;
  }

  void add_owned(LocalCellId cell, int component, double value) const {
    if (cell >= topology->owned_cell_count())
      return;
    const double candidate =
        scratch[cell * 3U + static_cast<std::size_t>(component)] + value;
    if (!std::isfinite(candidate))
      throw Error("immersed operator residual accumulation is non-finite");
    scratch[cell * 3U + static_cast<std::size_t>(component)] = candidate;
  }

  std::array<double, 3> viscous_traction(
      const SharedFace &face, const runtime::FieldView<const double> &velocity,
      const runtime::FieldView<const double> &gradients, double mu) const {
    if (!(mu >= 0.0) || !std::isfinite(mu))
      throw Error("immersed operator viscosity is invalid");
    const Real3 u_p{read_cell(velocity, face.owner, 0),
                    read_cell(velocity, face.owner, 1),
                    read_cell(velocity, face.owner, 2)};
    const Real3 u_n{read_cell(velocity, face.neighbour, 0),
                    read_cell(velocity, face.neighbour, 1),
                    read_cell(velocity, face.neighbour, 2)};
    const double a = norm(subtract(face.face_center, face.owner_center));
    const double distance = norm(face.displacement);
    const double b = distance - a;
    if (!(b > 0.0) || !std::isfinite(b))
      throw Error("immersed operator interpolation distance is invalid");
    std::array<double, 9> gradient{};
    for (int value = 0; value < 9; ++value)
      gradient[static_cast<std::size_t>(value)] =
          (b / distance) * read_cell(gradients, face.owner, value) +
          (a / distance) * read_cell(gradients, face.neighbour, value);
    const double d2 = dot(face.displacement, face.displacement);
    const std::array<double, 3> delta{u_n.x - u_p.x, u_n.y - u_p.y,
                                      u_n.z - u_p.z};
    const std::array<double, 3> d{face.displacement.x, face.displacement.y,
                                  face.displacement.z};
    for (std::size_t component = 0U; component < 3U; ++component) {
      const std::size_t base = component * 3U;
      const double projected = gradient[base] * d[0] +
                               gradient[base + 1U] * d[1] +
                               gradient[base + 2U] * d[2];
      const double correction = (delta[component] - projected) / d2;
      for (std::size_t direction = 0U; direction < 3U; ++direction)
        gradient[base + direction] += correction * d[direction];
    }
    const double divergence = gradient[0] + gradient[4] + gradient[8];
    const std::array<double, 3> area{face.area.x, face.area.y, face.area.z};
    std::array<double, 3> traction{};
    for (std::size_t component = 0U; component < 3U; ++component)
      for (std::size_t direction = 0U; direction < 3U; ++direction) {
        double stress = mu * (gradient[component * 3U + direction] +
                              gradient[direction * 3U + component]);
        if (component == direction)
          stress -= mu * (2.0 / 3.0) * divergence;
        traction[component] += stress * area[direction];
      }
    for (const double value : traction)
      if (!std::isfinite(value))
        throw Error("immersed operator viscous traction is non-finite");
    return traction;
  }

  double scalar_diffusive_flux(
      const SharedFace &face, const runtime::FieldView<const double> &values,
      const runtime::FieldView<const double> &gradients, double gamma) const {
    if (!(gamma >= 0.0) || !std::isfinite(gamma))
      throw Error("immersed operator diffusion coefficient is invalid");
    const double sd = dot(face.area, face.displacement);
    const double factor = dot(face.area, face.area) / sd;
    if (!(sd > 0.0) || !std::isfinite(factor))
      throw Error("immersed operator diffusion projection is invalid");
    const double distance = norm(face.displacement);
    const double a = norm(subtract(face.face_center, face.owner_center));
    const double b = distance - a;
    if (!(b > 0.0) || !std::isfinite(b))
      throw Error("immersed operator diffusion distance is invalid");
    Real3 gradient{};
    const double weights[2]{b / distance, a / distance};
    double *components[3]{&gradient.x, &gradient.y, &gradient.z};
    for (int direction = 0; direction < 3; ++direction)
      *components[direction] =
          weights[0] * read_cell(gradients, face.owner, direction) +
          weights[1] * read_cell(gradients, face.neighbour, direction);
    const Real3 nonorth =
        subtract(face.area, multiply(factor, face.displacement));
    const double flux = gamma * ((read_cell(values, face.neighbour, 0) -
                                  read_cell(values, face.owner, 0)) *
                                     factor +
                                 dot(gradient, nonorth));
    if (!std::isfinite(flux))
      throw Error("immersed operator diffusive flux is non-finite");
    return flux;
  }
};

namespace {

std::string current_exception_message() {
  try {
    throw;
  } catch (const std::exception &error) {
    return error.what();
  } catch (...) {
    return "immersed operator failed with an unknown exception";
  }
}

template <class ImplType, class Function>
void prepare_collectively(const ImplType &impl, Function &&function) {
  bool local_ok = true;
  std::string message;
  try {
    function();
  } catch (...) {
    local_ok = false;
    message = current_exception_message();
  }
  const auto status = runtime::collective_status(*impl.mpi, local_ok, message);
  if (!status.ok)
    throw Error(status.message + " (lowest failing rank " +
                std::to_string(status.failing_rank) + ")");
}

template <class ImplType, class View>
void validate_cell_view(const ImplType &impl, const View &view,
                        std::uint32_t components, bool require_ghost) {
  if (!same(view.interior_extent(), impl.local_extent) ||
      view.components() != components ||
      (require_ghost && view.ghost_width() < 1))
    throw Error("immersed operator cell field layout is invalid");
}

template <class ImplType, class View>
void validate_face_view(const ImplType &impl, const View &view,
                        std::uint32_t components) {
  if (view.face_count() != impl.topology->local_face_count() ||
      view.components() != components)
    throw Error("immersed operator face field layout is invalid");
}

using AffineBoundaryFunctional =
    immersed::detail::QuadraticReconstructionWeights::AffineBoundaryFunctional;

AffineBoundaryFunctional unconstrained_value_functional(
    const immersed::QuadraticReconstruction &reconstruction,
    runtime::Real3 point_m) {
  return {immersed::detail::QuadraticReconstructionWeights::value_weights(
              reconstruction, point_m),
          0.0};
}

AffineBoundaryFunctional normal_gradient_constrained_value_functional(
    const immersed::QuadraticReconstruction &reconstruction,
    runtime::Real3 point_m) {
  return immersed::detail::QuadraticReconstructionWeights::
      origin_normal_gradient_constrained_value_weights(reconstruction,
                                                       point_m);
}

void add_pressure_functional(AffineDonorAccumulator &donors,
                             WallGradientAccumulator &wall_gradients,
                             std::uint32_t output_component,
                             immersed::ImmersedLinkId link,
                             std::size_t wall_link_index, double scale,
                             std::uint64_t source_term_id,
                             const AffineBoundaryFunctional &functional) {
  add_affine_weights(donors, AffineInputKind::pressure, 0U, output_component,
                     link, wall_link_index, scale, source_term_id,
                     functional.donors);
  add_wall_gradient_weight(wall_gradients, link, output_component,
                           scale * functional.boundary_coefficient,
                           source_term_id);
}

PressureAffinePlan build_background_pressure_affine_plan(
    mesh::GlobalCellId cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links,
    const immersed::GhostStencilPlan &ghost_plan) {
  AffineDonorAccumulator donors;
  WallGradientAccumulator wall_gradients;
  for (const auto link_index : link_indices) {
    const auto &link = links[link_index];
    const auto &reconstruction = ghost_plan.reconstruction(link.id);
    const auto face_value = unconstrained_value_functional(
        reconstruction, link.pressure_quadrature_m);
    const double area[3]{link.area_from_fluid.x, link.area_from_fluid.y,
                         link.area_from_fluid.z};
    for (std::uint32_t output = 0U; output < 3U; ++output) {
      const auto source = physical_term_id(
          cell, link.id, PhysicalTermKind::pressure_direct,
          static_cast<std::uint32_t>(link.occurrence), output);
      add_pressure_functional(donors, wall_gradients, output, link.id,
                              link_index, area[output], source, face_value);
    }
  }
  return {finalize_affine_donor_terms(cell, donors),
          finalize_wall_gradient_terms(cell, wall_gradients)};
}

std::vector<AffineDonorTerm> build_background_viscous_affine_plan(
    mesh::GlobalCellId cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links,
    const immersed::GhostStencilPlan &ghost_plan) {
  AffineDonorAccumulator donors;
  constexpr std::array<Real3, 3> axes{Real3{1.0, 0.0, 0.0},
                                      Real3{0.0, 1.0, 0.0},
                                      Real3{0.0, 0.0, 1.0}};
  for (const auto link_index : link_indices) {
    const auto &link = links[link_index];
    const auto &reconstruction = ghost_plan.reconstruction(link.id);
    const auto displacement = subtract(
        link.local_sample_points_m[link.occurrence],
        link.local_sample_points_m[6]);
    const double displacement_squared = dot(displacement, displacement);
    if (!(displacement_squared > 0.0) ||
        !std::isfinite(displacement_squared))
      throw Error("immersed operator background displacement is invalid");
    const double d[3]{displacement.x, displacement.y, displacement.z};
    const double area[3]{link.area_from_fluid.x, link.area_from_fluid.y,
                         link.area_from_fluid.z};
    const auto add_derivative =
        [&](std::uint32_t input, std::uint32_t output,
            std::uint32_t direction, double stress_scale,
            PhysicalTermKind kind, std::uint64_t source) {
          const double orthogonal_scale =
              stress_scale * d[direction] / displacement_squared;
          if (kind == PhysicalTermKind::viscous_orthogonal) {
            add_affine_weights(
                donors, AffineInputKind::velocity, input, output, link.id,
                link_index, orthogonal_scale, source,
                immersed::detail::QuadraticReconstructionWeights::
                    origin_constrained_directional_gradient_weights(
                        reconstruction, link.pressure_quadrature_m,
                        displacement));
            return;
          }
          add_affine_weights(
              donors, AffineInputKind::velocity, input, output, link.id,
              link_index, stress_scale, source,
              immersed::detail::QuadraticReconstructionWeights::
                  origin_constrained_directional_gradient_weights(
                      reconstruction, link.pressure_quadrature_m,
                      axes[direction]));
          add_affine_weights(
              donors, AffineInputKind::velocity, input, output, link.id,
              link_index, -orthogonal_scale, source,
              immersed::detail::QuadraticReconstructionWeights::
                  origin_constrained_directional_gradient_weights(
                      reconstruction, link.pressure_quadrature_m,
                      displacement));
        };
    for (std::uint32_t output = 0U; output < 3U; ++output) {
      for (const auto kind : {PhysicalTermKind::viscous_orthogonal,
                              PhysicalTermKind::viscous_deferred_gradient}) {
        const auto source = physical_term_id(
            cell, link.id, kind,
            static_cast<std::uint32_t>(link.occurrence), output);
        for (std::uint32_t input = 0U; input < 3U; ++input) {
          if (input == output)
            for (std::uint32_t direction = 0U; direction < 3U; ++direction)
              add_derivative(input, output, direction, -area[direction], kind,
                             source);
          add_derivative(input, output, output, -area[input], kind, source);
          add_derivative(input, output, input,
                         (2.0 / 3.0) * area[output], kind, source);
        }
      }
    }
  }
  return finalize_affine_donor_terms(cell, donors);
}

PressureAffinePlan build_pressure_affine_plan(
    mesh::GlobalCellId cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links,
    const immersed::GhostStencilPlan &ghost_plan,
    bool constrain_wall_normal_gradient) {
  AffineDonorAccumulator donors;
  WallGradientAccumulator wall_gradients;
  for (const auto link_index : link_indices) {
    const auto &link = links[link_index];
    const auto &reconstruction = ghost_plan.reconstruction(link.id);
    const auto value_functional = [&](runtime::Real3 point) {
      return constrain_wall_normal_gradient
                 ? normal_gradient_constrained_value_functional(
                       reconstruction, point)
                 : unconstrained_value_functional(reconstruction, point);
    };
    const auto face_value = value_functional(link.pressure_quadrature_m);
    const auto add_remainder = [&](std::uint32_t output, double scale,
                                   std::uint64_t source,
                                   runtime::Real3 point) {
      add_pressure_functional(donors, wall_gradients, output, link.id,
                              link_index, scale, source,
                              value_functional(point));
      add_pressure_functional(donors, wall_gradients, output, link.id,
                              link_index, -scale, source,
                              value_functional(link.wall_intercept_m));
      const auto direction = subtract(point, link.wall_intercept_m);
      const auto gradient_functional =
          constrain_wall_normal_gradient
              ? immersed::detail::QuadraticReconstructionWeights::
                    origin_normal_gradient_constrained_directional_gradient_weights(
                        reconstruction, link.wall_intercept_m, direction)
              : AffineBoundaryFunctional{
                    immersed::detail::QuadraticReconstructionWeights::
                        directional_gradient_weights(
                            reconstruction, link.wall_intercept_m, direction),
                    0.0};
      add_pressure_functional(donors, wall_gradients, output, link.id,
                              link_index, -scale, source,
                              gradient_functional);
    };
    const double area[3]{link.area_from_fluid.x, link.area_from_fluid.y,
                         link.area_from_fluid.z};
    if (link.transformed_row.source != link.background_row.source)
      throw Error("immersed pressure LFP source difference is unsupported");
    for (std::uint32_t output = 0U; output < 3U; ++output) {
      const auto face_term = replacement_term_id(
          cell, link.id, BoundaryReplacementKind::pressure_face, 7U, output);
      const auto diagonal_defect_term = replacement_term_id(
          cell, link.id, BoundaryReplacementKind::pressure_diagonal_defect,
          6U, output);
      add_pressure_functional(donors, wall_gradients, output, link.id,
                              link_index, area[output], face_term, face_value);
      add_remainder(output,
                    area[output] * (link.transformed_row.diagonal -
                                    link.background_row.diagonal),
                    diagonal_defect_term, link.local_sample_points_m[6]);
      for (std::uint32_t occurrence = 0U; occurrence < 6U; ++occurrence) {
        const auto neighbour_defect_term = replacement_term_id(
            cell, link.id, BoundaryReplacementKind::pressure_neighbour_defect,
            occurrence, output);
        add_remainder(output,
                      area[output] *
                          (link.transformed_row.neighbour[occurrence] -
                           link.background_row.neighbour[occurrence]),
                      neighbour_defect_term,
                      link.local_sample_points_m[occurrence]);
      }
    }
  }
  return {finalize_affine_donor_terms(cell, donors),
          finalize_wall_gradient_terms(cell, wall_gradients)};
}

Int3 global_cell_from_id(mesh::GlobalCellId id, Int3 extent) {
  const auto nx = static_cast<std::uint64_t>(extent.x);
  const auto ny = static_cast<std::uint64_t>(extent.y);
  const auto plane = nx * ny;
  return {static_cast<int>(id % nx), static_cast<int>((id / nx) % ny),
          static_cast<int>(id / plane)};
}

std::array<Int3, 4> logical_face_vertices(mesh::LogicalFace face) {
  const auto coordinate = face.coordinate;
  switch (face.axis) {
  case mesh::FaceAxis::x:
    return {Int3{coordinate.x, coordinate.y, coordinate.z},
            Int3{coordinate.x, coordinate.y + 1, coordinate.z},
            Int3{coordinate.x, coordinate.y + 1, coordinate.z + 1},
            Int3{coordinate.x, coordinate.y, coordinate.z + 1}};
  case mesh::FaceAxis::y:
    return {Int3{coordinate.x, coordinate.y, coordinate.z},
            Int3{coordinate.x, coordinate.y, coordinate.z + 1},
            Int3{coordinate.x + 1, coordinate.y, coordinate.z + 1},
            Int3{coordinate.x + 1, coordinate.y, coordinate.z}};
  case mesh::FaceAxis::z:
    return {Int3{coordinate.x, coordinate.y, coordinate.z},
            Int3{coordinate.x + 1, coordinate.y, coordinate.z},
            Int3{coordinate.x + 1, coordinate.y + 1, coordinate.z},
            Int3{coordinate.x, coordinate.y + 1, coordinate.z}};
  }
  throw Error("immersed pressure row face axis is invalid");
}

std::array<std::array<double, immersed::detail::kQuadraticBasisSize>, 3>
single_face_pressure_functional(
    LocalFaceId face, LocalCellId cell,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::detail::QuadraticFrame &frame) {
  constexpr std::array<std::array<std::size_t, 3>, 2> triangles{{
      {{0U, 1U, 2U}},
      {{0U, 2U, 3U}},
  }};
  constexpr std::array<std::array<double, 3>, 3> barycentric{{
      {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}},
  }};
  std::array<std::array<double, immersed::detail::kQuadraticBasisSize>, 3>
      result{};
  const auto neighbour = topology.neighbour(face);
  const bool owner_side = topology.owner(face) == cell;
  const bool neighbour_side = neighbour.has_value() && *neighbour == cell;
  if (!owner_side && !neighbour_side)
    throw Error("immersed pressure face association is invalid");
  const auto logical = logical_face_vertices(topology.logical_face(face));
  std::array<Real3, 4> vertices{};
  for (std::size_t index = 0U; index < vertices.size(); ++index)
    vertices[index] = geometry.vertex_position_m(logical[index]);
  const auto oriented_area = geometry.face_area_vector_m2(
      face, owner_side ? mesh::FaceSide::owner : mesh::FaceSide::neighbour);
  Real3 raw_area{};
  std::array<std::array<double, immersed::detail::kQuadraticBasisSize>, 3>
      raw_integral{};
  for (const auto &triangle : triangles) {
    const auto &a = vertices[triangle[0]];
    const auto &b = vertices[triangle[1]];
    const auto &c = vertices[triangle[2]];
    const auto vector_area = multiply(
        0.5, cross(subtract(b, a), subtract(c, a)));
    raw_area = add(raw_area, vector_area);
    for (const auto &weight : barycentric) {
      const auto point =
          add(add(multiply(weight[0], a), multiply(weight[1], b)),
              multiply(weight[2], c));
      const auto basis = immersed::detail::quadratic_basis_at(point, frame);
      const double area_components[3]{vector_area.x, vector_area.y,
                                      vector_area.z};
      for (std::size_t component = 0U; component < 3U; ++component)
        for (std::size_t mode = 0U; mode < basis.size(); ++mode)
          raw_integral[component][mode] +=
              area_components[component] * basis[mode] / 3.0;
    }
  }
  const double orientation = dot(raw_area, oriented_area) >= 0.0 ? 1.0 : -1.0;
  for (std::size_t component = 0U; component < 3U; ++component)
    for (std::size_t mode = 0U;
         mode < immersed::detail::kQuadraticBasisSize; ++mode)
      result[component][mode] += orientation * raw_integral[component][mode];
  return result;
}

std::array<std::array<double, immersed::detail::kQuadraticBasisSize>, 3>
complete_cell_pressure_functionals(
    LocalCellId cell, const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry,
    const immersed::detail::QuadraticFrame &frame,
    const std::vector<LocalFaceId> &cell_faces) {
  std::array<std::array<double, immersed::detail::kQuadraticBasisSize>, 3>
      result{};
  if (cell_faces.size() != 6U)
    throw Error("immersed pressure row must have six background faces");
  for (const auto face : cell_faces) {
    if (face >= topology.local_face_count())
      throw Error("immersed pressure row face index is invalid");
    const auto face_integral = single_face_pressure_functional(
        face, cell, topology, geometry, frame);
    for (std::size_t component = 0U; component < 3U; ++component)
      for (std::size_t mode = 0U;
           mode < immersed::detail::kQuadraticBasisSize; ++mode)
        result[component][mode] += face_integral[component][mode];
  }
  return result;
}

PressureAffinePlan build_complete_pressure_boundary_row_plan(
    LocalCellId row_cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links,
    const std::vector<SharedFace> &shared_faces,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::QuadraticReconstruction &row_reconstruction,
    const immersed::GhostStencilPlan &ghost_plan,
    bool force_authority,
    const std::vector<BoundaryReplacementTerm> &replacement_terms,
    const std::vector<LocalFaceId> &cell_faces,
    const std::vector<std::size_t> &cell_shared_faces) {
  constexpr std::size_t basis_size = immersed::detail::kQuadraticBasisSize;
  if (link_indices.empty() || link_indices.size() >= basis_size)
    throw Error("immersed pressure row constraint count is invalid");
  const mesh::GlobalCellId row_global = topology.global_cell_id(row_cell);
  immersed::detail::QuadraticFrame frame{
      geometry.cell_center_m(row_cell), Real3{1.0, 0.0, 0.0},
      Real3{0.0, 1.0, 0.0}, Real3{0.0, 0.0, 1.0},
      std::cbrt(geometry.cell_volume_m3(row_cell))};

  const auto &authority_donors =
      immersed::detail::QuadraticReconstructionWeights::donor_global_ids(
          row_reconstruction);
  std::vector<mesh::GlobalCellId> donors(authority_donors.begin(),
                                          authority_donors.end());
  if (donors.size() < basis_size)
    throw Error("immersed pressure row donor union is incomplete");

  std::vector<double> design(donors.size() * basis_size, 0.0);
  for (std::size_t donor = 0U; donor < donors.size(); ++donor) {
    const auto moment = immersed::detail::quadratic_cell_average_basis(
        global_cell_from_id(donors[donor], topology.global_extent()), frame,
        topology, geometry);
    std::copy(moment.begin(), moment.end(),
              design.begin() + static_cast<std::ptrdiff_t>(donor * basis_size));
  }

  const std::size_t constraint_count = link_indices.size();
  std::vector<double> constraints(constraint_count * basis_size, 0.0);
  for (std::size_t constraint = 0U; constraint < constraint_count;
       ++constraint) {
    const auto &link = links[link_indices[constraint]];
    const auto derivative =
        immersed::detail::quadratic_directional_derivative_basis(
            link.wall_intercept_m, link.solid_to_fluid_normal, frame);
    std::copy(derivative.begin(), derivative.end(),
              constraints.begin() +
                  static_cast<std::ptrdiff_t>(constraint * basis_size));
  }

  std::vector<double> constraint_transpose(basis_size * constraint_count,
                                           0.0);
  for (std::size_t row = 0U; row < basis_size; ++row)
    for (std::size_t column = 0U; column < constraint_count; ++column)
      constraint_transpose[row * constraint_count + column] =
          constraints[column * basis_size + row];
  const auto constraint_qr = immersed::detail::factorize_design_matrix(
      constraint_transpose, basis_size, constraint_count);
  std::vector<double> particular(basis_size * constraint_count, 0.0);
  for (std::size_t constraint = 0U; constraint < constraint_count;
       ++constraint) {
    std::vector<double> functional(constraint_count, 0.0);
    functional[constraint] = 1.0;
    const auto column = constraint_qr.functional_weights(functional);
    for (std::size_t coefficient = 0U; coefficient < basis_size;
         ++coefficient)
      particular[coefficient * constraint_count + constraint] =
          column[coefficient];
  }

  const std::size_t free_count = basis_size - constraint_count;
  std::vector<double> nullspace;
  nullspace.reserve(basis_size * free_count);
  std::size_t accepted_columns = 0U;
  for (std::size_t axis = 0U;
       axis < basis_size && accepted_columns < free_count; ++axis) {
    std::vector<double> vector(basis_size, 0.0);
    vector[axis] = 1.0;
    for (std::size_t column = 0U; column < constraint_count; ++column) {
      const double projection = constraint_qr.thin_q[axis * constraint_count + column];
      for (std::size_t row = 0U; row < basis_size; ++row)
        vector[row] -=
            projection * constraint_qr.thin_q[row * constraint_count + column];
    }
    for (int pass = 0; pass < 2; ++pass)
      for (std::size_t column = 0U; column < accepted_columns; ++column) {
        double projection = 0.0;
        for (std::size_t row = 0U; row < basis_size; ++row)
          projection += vector[row] * nullspace[row * free_count + column];
        for (std::size_t row = 0U; row < basis_size; ++row)
          vector[row] -= projection * nullspace[row * free_count + column];
      }
    double magnitude = 0.0;
    for (const double value : vector)
      magnitude = std::hypot(magnitude, value);
    if (magnitude <= 4096.0 * std::numeric_limits<double>::epsilon())
      continue;
    if (nullspace.empty())
      nullspace.resize(basis_size * free_count, 0.0);
    for (std::size_t row = 0U; row < basis_size; ++row)
      nullspace[row * free_count + accepted_columns] = vector[row] / magnitude;
    ++accepted_columns;
  }
  if (accepted_columns != free_count)
    throw Error("immersed pressure row constraint nullspace is rank deficient");

  std::vector<double> reduced_design(donors.size() * free_count, 0.0);
  for (std::size_t donor = 0U; donor < donors.size(); ++donor)
    for (std::size_t free = 0U; free < free_count; ++free)
      for (std::size_t coefficient = 0U; coefficient < basis_size;
           ++coefficient)
        reduced_design[donor * free_count + free] +=
            design[donor * basis_size + coefficient] *
            nullspace[coefficient * free_count + free];
  const auto reduced_qr = immersed::detail::factorize_design_matrix(
      reduced_design, donors.size(), free_count);
  auto complete_functionals =
      complete_cell_pressure_functionals(row_cell, topology, geometry, frame,
                                         cell_faces);

  // Force-authority variant (M2): when `force_authority` is set, the row's
  // wall-face pressure value is re-anchored from the background grid-face
  // value to the shared per-link wall-anchored authority reconstruction at
  // the wall intercept, weighted by the link's true body-surface measure:
  //   flux_c = A_surface_c * p_authority(wall)
  //   (replacing A_face_c * p_face, both expressed over the row donor union
  //   and the wall-gradient constraints)
  // This variant is used only by the final force collection; the solve rows
  // keep the background face-flux structure (see the force-authority
  // semantic amendment). For a constant field the replacement preserves the
  // closed-surface annihilation in the aggregate.
  struct WallValueFunctional final {
    std::vector<double> donor_weights;
    std::vector<double> wall_weights;
  };
  std::vector<WallValueFunctional> face_wall_values;
  std::vector<AffineBoundaryFunctional> authority_wall_values;
  if (force_authority) {
    face_wall_values.resize(constraint_count);
    authority_wall_values.resize(constraint_count);
    for (std::size_t constraint = 0U; constraint < constraint_count;
         ++constraint) {
      const auto &link = links[link_indices[constraint]];
      const auto value_basis = immersed::detail::quadratic_basis_at(
          link.pressure_quadrature_m, frame);
      std::vector<double> free_value(free_count, 0.0);
      for (std::size_t free = 0U; free < free_count; ++free)
        for (std::size_t coefficient = 0U; coefficient < basis_size;
             ++coefficient)
          free_value[free] +=
              value_basis[coefficient] *
              nullspace[coefficient * free_count + free];
      auto &face_value = face_wall_values[constraint];
      face_value.donor_weights = reduced_qr.functional_weights(free_value);
      face_value.wall_weights.assign(constraint_count, 0.0);
      for (std::size_t wall = 0U; wall < constraint_count; ++wall) {
        for (std::size_t coefficient = 0U; coefficient < basis_size;
             ++coefficient)
          face_value.wall_weights[wall] +=
              value_basis[coefficient] *
              particular[coefficient * constraint_count + wall];
        for (std::size_t donor = 0U; donor < donors.size(); ++donor) {
          double design_particular = 0.0;
          for (std::size_t coefficient = 0U; coefficient < basis_size;
               ++coefficient)
            design_particular +=
                design[donor * basis_size + coefficient] *
                particular[coefficient * constraint_count + wall];
          face_value.wall_weights[wall] -=
              face_value.donor_weights[donor] * design_particular;
        }
      }
      {
        // Canonical sharp-interface boundary pressure: the per-link authority
        // reconstruction's extrapolated value at the body-surface patch
        // centroid (the ghost-cell-IBM extrapolate_scalar structure). The
        // constrained wall-gradient value is deliberately not used here: its
        // gradient-datum coupling introduces a coherent O(h) quadrature
        // component (see the force-authority semantic amendment).
        AffineBoundaryFunctional wall_value{};
        for (const auto &donor :
             immersed::detail::QuadraticReconstructionWeights::value_weights(
                 ghost_plan.reconstruction(link.id),
                 link.surface_patch_centroid_m))
          wall_value.donors.push_back(donor);
        authority_wall_values[constraint] = std::move(wall_value);
      }
    }
  }

  PressureAffinePlan result{};
  for (std::uint32_t output = 0U; output < 3U; ++output) {
    std::vector<double> free_functional(free_count, 0.0);
    for (std::size_t free = 0U; free < free_count; ++free)
      for (std::size_t coefficient = 0U; coefficient < basis_size;
           ++coefficient)
        free_functional[free] +=
            nullspace[coefficient * free_count + free] *
            complete_functionals[output][coefficient];
    auto donor_weights = reduced_qr.functional_weights(free_functional);
    std::vector<double> wall_weights(constraint_count, 0.0);
    for (std::size_t constraint = 0U; constraint < constraint_count;
         ++constraint) {
      for (std::size_t coefficient = 0U; coefficient < basis_size;
           ++coefficient)
        wall_weights[constraint] +=
            complete_functionals[output][coefficient] *
            particular[coefficient * constraint_count + constraint];
      for (std::size_t donor = 0U; donor < donors.size(); ++donor) {
        double design_particular = 0.0;
        for (std::size_t coefficient = 0U; coefficient < basis_size;
             ++coefficient)
          design_particular += design[donor * basis_size + coefficient] *
                               particular[coefficient * constraint_count +
                                          constraint];
        wall_weights[constraint] -= donor_weights[donor] * design_particular;
      }
    }

    for (const auto face_index : cell_shared_faces) {
      if (face_index >= shared_faces.size())
        throw Error("immersed pressure row shared-face index is invalid");
      const auto &face = shared_faces[face_index];
      const bool row_owner = face.owner == row_cell;
      const bool row_neighbour = face.neighbour == row_cell;
      if (!row_owner && !row_neighbour)
        throw Error("immersed pressure row shared-face association is invalid");
      const double distance = norm(face.displacement);
      const double owner_to_face =
          norm(subtract(face.face_center, face.owner_center));
      const double owner_weight = (distance - owner_to_face) / distance;
      const double neighbour_weight = owner_to_face / distance;
      const double oriented_area =
          (row_owner ? 1.0 : -1.0) *
          (output == 0U ? face.area.x
                        : output == 1U ? face.area.y : face.area.z);
      const auto subtract_shared = [&](LocalCellId cell, double weight) {
        const auto id = topology.global_cell_id(cell);
        const auto found = std::lower_bound(donors.begin(), donors.end(), id);
        if (found == donors.end() || *found != id)
          throw Error("immersed pressure row shared donor is missing");
        donor_weights[static_cast<std::size_t>(found - donors.begin())] -=
            oriented_area * weight;
      };
      subtract_shared(face.owner, owner_weight);
      subtract_shared(face.neighbour, neighbour_weight);
    }

    if (force_authority) {
      for (std::size_t constraint = 0U; constraint < constraint_count;
           ++constraint) {
        const auto &link = links[link_indices[constraint]];
        const double surface_component =
            output == 0U   ? link.surface_measure_m2.x
            : output == 1U ? link.surface_measure_m2.y
                           : link.surface_measure_m2.z;
        const double face_component =
            output == 0U   ? link.area_from_fluid.x
            : output == 1U ? link.area_from_fluid.y
                           : link.area_from_fluid.z;
        const auto &authority_value = authority_wall_values[constraint];
        const auto &face_value = face_wall_values[constraint];
        for (const auto &donor : authority_value.donors) {
          const auto found =
              std::lower_bound(donors.begin(), donors.end(),
                               donor.global_cell);
          if (found == donors.end() || *found != donor.global_cell)
            throw Error("immersed pressure row authority donor is missing");
          donor_weights[static_cast<std::size_t>(
              found - donors.begin())] +=
              -surface_component * donor.weight;
        }
        for (std::size_t donor = 0U; donor < donors.size(); ++donor)
          donor_weights[donor] -=
              face_component * face_value.donor_weights[donor];
        wall_weights[constraint] +=
            -surface_component * authority_value.boundary_coefficient;
        for (std::size_t wall = 0U; wall < constraint_count; ++wall)
          wall_weights[wall] -=
              face_component * face_value.wall_weights[wall];
      }
    }

    std::vector<std::uint64_t> sources;
    for (const auto &term : replacement_terms)
      if (term.component == output &&
          (term.kind == BoundaryReplacementKind::pressure_face ||
           term.kind == BoundaryReplacementKind::pressure_diagonal_defect ||
           term.kind == BoundaryReplacementKind::pressure_neighbour_defect) &&
          term.scale != 0.0)
        sources.push_back(term.stable_id);
    std::sort(sources.begin(), sources.end());
    sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
    for (std::size_t donor = 0U; donor < donors.size(); ++donor) {
      const double coefficient = donor_weights[donor];
      if (coefficient == 0.0)
        continue;
      AffineDonorTerm term{};
      term.key = {AffineInputKind::pressure, donors[donor], 0U, output};
      term.stable_id = affine_donor_term_id(row_global, term.key);
      double assigned = 0.0;
      for (std::size_t index = 0U; index < link_indices.size(); ++index) {
        const double share =
            index + 1U == link_indices.size()
                ? coefficient - assigned
                : coefficient / static_cast<double>(link_indices.size());
        assigned += share;
        const auto link_index = link_indices[index];
        term.link_coefficients.push_back(
            {links[link_index].id, link_index, share,
             index == 0U ? sources : std::vector<std::uint64_t>{}});
      }
      result.donor_terms.push_back(std::move(term));
    }
    for (std::size_t constraint = 0U; constraint < constraint_count;
         ++constraint) {
      if (wall_weights[constraint] == 0.0)
        continue;
      const auto link = links[link_indices[constraint]].id;
      result.wall_gradient_terms.push_back(
          {affine_wall_gradient_term_id(row_global, link, output), link,
           output, wall_weights[constraint], {}});
    }
  }
  std::sort(result.donor_terms.begin(), result.donor_terms.end(),
            [](const auto &left, const auto &right) {
              return left.key < right.key;
            });
  std::sort(result.wall_gradient_terms.begin(),
            result.wall_gradient_terms.end(), [](const auto &left,
                                                 const auto &right) {
              return std::tie(left.link, left.output_component) <
                     std::tie(right.link, right.output_component);
            });
  return result;
}

immersed::QuadraticReconstruction build_boundary_row_reconstruction(
    LocalCellId row_cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links,
    const std::vector<SharedFace> &shared_faces,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::GhostStencilPlan &ghost_plan,
    const std::vector<std::size_t> &cell_shared_faces) {
  std::vector<mesh::GlobalCellId> donor_ids;
  for (const auto link_index : link_indices) {
    const auto &ids = immersed::detail::QuadraticReconstructionWeights::
        donor_global_ids(ghost_plan.reconstruction(links[link_index].id));
    donor_ids.insert(donor_ids.end(), ids.begin(), ids.end());
  }
  for (const auto face_index : cell_shared_faces) {
    if (face_index >= shared_faces.size())
      throw Error("immersed boundary row shared-face index is invalid");
    const auto &face = shared_faces[face_index];
    if (face.owner != row_cell && face.neighbour != row_cell)
      throw Error("immersed boundary row shared-face association is invalid");
    donor_ids.push_back(topology.global_cell_id(face.owner));
    donor_ids.push_back(topology.global_cell_id(face.neighbour));
  }
  std::sort(donor_ids.begin(), donor_ids.end());
  donor_ids.erase(std::unique(donor_ids.begin(), donor_ids.end()),
                  donor_ids.end());
  std::vector<Int3> donors;
  donors.reserve(donor_ids.size());
  for (const auto donor : donor_ids)
    donors.push_back(global_cell_from_id(donor, topology.global_extent()));
  const auto anchor = topology.global_cell(row_cell);
  immersed::detail::BoundaryAuthorityCoverageScope coverage;
  return immersed::QuadraticReconstruction::create(
      geometry.cell_center_m(row_cell), Real3{1.0, 0.0, 0.0},
      Real3{0.0, 1.0, 0.0}, Real3{0.0, 0.0, 1.0},
      std::cbrt(geometry.cell_volume_m3(row_cell)), anchor, donors, topology,
      geometry);
}

std::vector<AffineDonorTerm> build_viscous_affine_plan(
    mesh::GlobalCellId cell, const std::vector<std::size_t> &link_indices,
    const std::vector<WallLink> &links,
    const immersed::GhostStencilPlan &ghost_plan) {
  AffineDonorAccumulator donors;
  constexpr std::array<Real3, 3> axes{Real3{1.0, 0.0, 0.0},
                                      Real3{0.0, 1.0, 0.0},
                                      Real3{0.0, 0.0, 1.0}};
  for (const auto link_index : link_indices) {
    const auto &link = links[link_index];
    const auto &reconstruction = ghost_plan.reconstruction(link.id);
    const double area[3]{link.area_from_fluid.x, link.area_from_fluid.y,
                         link.area_from_fluid.z};
    const auto add_traction = [&](runtime::Real3 point, double scale,
                                  BoundaryReplacementKind kind,
                                  std::uint32_t occurrence) {
      for (std::uint32_t output = 0U; output < 3U; ++output) {
        const auto source_term_id = replacement_term_id(
            cell, link.id, kind, occurrence, output);
        for (std::uint32_t input = 0U; input < 3U; ++input) {
          if (input == output) {
            for (std::uint32_t direction = 0U; direction < 3U; ++direction)
              add_affine_weights(
                  donors, AffineInputKind::velocity, input, output, link.id,
                  link_index, scale * area[direction], source_term_id,
                  immersed::detail::QuadraticReconstructionWeights::
                      origin_constrained_directional_gradient_weights(
                          reconstruction, point, axes[direction]));
          }
          add_affine_weights(
              donors, AffineInputKind::velocity, input, output, link.id,
              link_index, scale * area[input], source_term_id,
              immersed::detail::QuadraticReconstructionWeights::
                  origin_constrained_directional_gradient_weights(
                      reconstruction, point, axes[output]));
          add_affine_weights(
              donors, AffineInputKind::velocity, input, output, link.id,
              link_index, scale * (-2.0 / 3.0) * area[output], source_term_id,
              immersed::detail::QuadraticReconstructionWeights::
                  origin_constrained_directional_gradient_weights(
                      reconstruction, point, axes[input]));
        }
      }
    };
    add_traction(link.wall_intercept_m, -1.0,
                 BoundaryReplacementKind::viscous_wall, 7U);
    add_traction(link.local_sample_points_m[6],
                 -(link.transformed_row.diagonal -
                   link.background_row.diagonal),
                 BoundaryReplacementKind::viscous_diagonal_defect, 6U);
    for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence)
      add_traction(link.local_sample_points_m[occurrence],
                   -(link.transformed_row.neighbour[occurrence] -
                     link.background_row.neighbour[occurrence]),
                   BoundaryReplacementKind::viscous_neighbour_defect,
                   static_cast<std::uint32_t>(occurrence));
  }
  return finalize_affine_donor_terms(cell, donors);
}

std::vector<std::uint64_t>
affine_source_term_ids(const std::vector<AffineDonorTerm> &donor_terms,
                       const std::vector<AffineWallGradientTerm> &wall_terms =
                           {}) {
  std::vector<std::uint64_t> result;
  for (const auto &term : donor_terms)
    for (const auto &coefficient : term.link_coefficients)
      result.insert(result.end(), coefficient.source_term_ids.begin(),
                    coefficient.source_term_ids.end());
  for (const auto &term : wall_terms)
    result.insert(result.end(), term.source_term_ids.begin(),
                  term.source_term_ids.end());
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

bool pressure_replacement(BoundaryReplacementKind kind) noexcept {
  return kind == BoundaryReplacementKind::pressure_face ||
         kind == BoundaryReplacementKind::pressure_diagonal_defect ||
         kind == BoundaryReplacementKind::pressure_neighbour_defect;
}

void validate_background_source_coverage(
    const std::vector<PhysicalTerm> &physical_terms,
    const PressureAffinePlan &pressure_unconstrained,
    const PressureAffinePlan &pressure_constrained,
    const std::vector<AffineDonorTerm> &viscous_terms) {
  const auto unconstrained_sources = affine_source_term_ids(
      pressure_unconstrained.donor_terms,
      pressure_unconstrained.wall_gradient_terms);
  const auto constrained_sources = affine_source_term_ids(
      pressure_constrained.donor_terms,
      pressure_constrained.wall_gradient_terms);
  const auto viscous_sources = affine_source_term_ids(viscous_terms);
  std::vector<std::uint64_t> declared;
  declared.reserve(physical_terms.size());
  for (const auto &term : physical_terms) {
    declared.push_back(term.stable_id);
    if (term.coefficient == 0.0)
      continue;
    const auto contains = [&](const std::vector<std::uint64_t> &sources) {
      return std::binary_search(sources.begin(), sources.end(),
                                term.stable_id);
    };
    if (term.kind == PhysicalTermKind::pressure_direct) {
      if (!contains(unconstrained_sources) || !contains(constrained_sources))
        throw Error(
            "immersed background pressure source coverage is incomplete");
    } else if ((term.kind == PhysicalTermKind::viscous_orthogonal ||
                term.kind ==
                    PhysicalTermKind::viscous_deferred_gradient) &&
               !contains(viscous_sources)) {
      throw Error(
          "immersed background viscous source coverage is incomplete");
    }
  }
  std::sort(declared.begin(), declared.end());
  const auto require_declared = [&](const std::vector<std::uint64_t> &sources) {
    for (const auto source : sources)
      if (!std::binary_search(declared.begin(), declared.end(), source))
        throw Error("immersed background affine source term is undeclared");
  };
  require_declared(unconstrained_sources);
  require_declared(constrained_sources);
  require_declared(viscous_sources);
}

void validate_affine_source_coverage(
    const std::vector<BoundaryReplacementTerm> &replacement_terms,
    const PressureAffinePlan &pressure_unconstrained,
    const PressureAffinePlan &pressure_constrained,
    const std::vector<AffineDonorTerm> &viscous_terms) {
  const auto pressure_sources = affine_source_term_ids(
      pressure_unconstrained.donor_terms,
      pressure_unconstrained.wall_gradient_terms);
  const auto constrained_sources = affine_source_term_ids(
      pressure_constrained.donor_terms,
      pressure_constrained.wall_gradient_terms);
  const auto viscous_sources = affine_source_term_ids(viscous_terms);
  for (const auto &term : replacement_terms) {
    if (term.scale == 0.0)
      continue;
    const auto contains = [&](const std::vector<std::uint64_t> &sources) {
      return std::binary_search(sources.begin(), sources.end(),
                                term.stable_id);
    };
    if (pressure_replacement(term.kind)) {
      if (!contains(pressure_sources) || !contains(constrained_sources))
        throw Error("immersed pressure affine source coverage is incomplete");
    } else if (!contains(viscous_sources)) {
      throw Error("immersed viscous affine source coverage is incomplete");
    }
  }
  std::vector<std::uint64_t> declared;
  declared.reserve(replacement_terms.size());
  for (const auto &term : replacement_terms)
    declared.push_back(term.stable_id);
  std::sort(declared.begin(), declared.end());
  const auto require_declared = [&](const std::vector<std::uint64_t> &sources) {
    for (const auto source : sources)
      if (!std::binary_search(declared.begin(), declared.end(), source))
        throw Error("immersed affine source term is undeclared");
  };
  require_declared(pressure_sources);
  require_declared(constrained_sources);
  require_declared(viscous_sources);
}

std::vector<AffineInputValueKey> bind_affine_input_snapshot_indices(
    PressureAffinePlan &background_unconstrained,
    PressureAffinePlan &background_constrained,
    std::vector<AffineDonorTerm> &background_viscous,
    PressureAffinePlan &wall_unconstrained,
    PressureAffinePlan &wall_constrained,
    PressureAffinePlan &wall_constrained_force,
    std::vector<AffineDonorTerm> &wall_viscous) {
  const std::array<std::vector<AffineDonorTerm> *, 7> plans{
      &background_unconstrained.donor_terms,
      &background_constrained.donor_terms,
      &background_viscous,
      &wall_unconstrained.donor_terms,
      &wall_constrained.donor_terms,
      &wall_constrained_force.donor_terms,
      &wall_viscous};
  std::vector<AffineInputValueKey> keys;
  for (const auto *terms : plans)
    for (const auto &term : *terms)
      keys.push_back(
          {term.key.input_kind, term.key.donor, term.key.input_component});
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end(),
                         [](const auto &left, const auto &right) {
                           return !(left < right) && !(right < left);
                         }),
             keys.end());
  for (auto *terms : plans)
    for (auto &term : *terms) {
      const AffineInputValueKey key{term.key.input_kind, term.key.donor,
                                    term.key.input_component};
      const auto position = std::lower_bound(keys.begin(), keys.end(), key);
      if (position == keys.end() || key < *position || *position < key)
        throw Error("immersed operator affine input binding is incomplete");
      term.input_snapshot_index =
          static_cast<std::size_t>(position - keys.begin());
    }
  return keys;
}

template <class ImplType>
BoundaryRowEvaluation evaluate_boundary_row_once(
    const ImplType &impl, const BoundaryRowPlan &row,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &pressure,
    const runtime::FaceFieldView<const double> &dynamic_viscosity_by_face,
    const std::vector<double> *wall_normal_gradients) {
  BoundaryRowEvaluation result{};
  const auto &pressure_plan = wall_normal_gradients == nullptr
                                  ? row.pressure_unconstrained
                                  : scoped_force_authority_evaluation
                                        ? row.pressure_constrained_force
                                        : row.pressure_constrained;
  const auto &background_pressure_plan =
      row.background_pressure_unconstrained;
  if (impl.affine_input_values.size() < row.affine_input_keys.size())
    throw Error("immersed operator affine input workspace is incomplete");
  for (std::size_t index = 0U; index < row.affine_input_keys.size(); ++index) {
    const auto &key = row.affine_input_keys[index];
    const auto &field = key.input_kind == AffineInputKind::pressure
                            ? pressure
                            : velocity;
    impl.affine_input_values[index] =
        impl.read_global(field, key.donor,
                         static_cast<int>(key.input_component));
  }
  if (impl.viscosity_by_wall_link.size() < impl.wall_links.size())
    throw Error("immersed operator viscosity workspace is incomplete");
  for (const auto link_index : row.wall_links) {
    const double value = dynamic_viscosity_by_face(
        impl.wall_links[link_index].face, 0);
    if (!(value >= 0.0) || !std::isfinite(value))
      throw Error("immersed operator viscosity is invalid");
    impl.viscosity_by_wall_link[link_index] = value;
  }
  const auto evaluate_donor_terms =
      [&](const std::vector<AffineDonorTerm> &terms,
          bool viscosity_scaled, std::array<double, 3> &target,
          [[maybe_unused]] bool background) {
        for (const auto &term : terms) {
          double effective_coefficient = 0.0;
          for (const auto &link_coefficient : term.link_coefficients) {
            if (link_coefficient.wall_link_index >= impl.wall_links.size() ||
                impl.wall_links[link_coefficient.wall_link_index].id !=
                    link_coefficient.link)
              throw Error(
                  "immersed operator affine link coefficient is invalid");
            double scale = 1.0;
            if (viscosity_scaled) {
              scale =
                  impl.viscosity_by_wall_link[link_coefficient.wall_link_index];
              if (!(scale >= 0.0) || !std::isfinite(scale))
                throw Error("immersed operator viscosity is invalid");
            }
            effective_coefficient += link_coefficient.coefficient * scale;
          }
          if (term.input_snapshot_index >= row.affine_input_keys.size())
            throw Error("immersed operator affine input snapshot is incomplete");
          const AffineInputValueKey expected{term.key.input_kind,
                                             term.key.donor,
                                             term.key.input_component};
          const auto &bound = row.affine_input_keys[term.input_snapshot_index];
          if (expected < bound || bound < expected)
            throw Error("immersed operator affine input binding is invalid");
          const double value =
              impl.affine_input_values[term.input_snapshot_index];
          target[term.key.output_component] += effective_coefficient * value;
          if (!std::isfinite(effective_coefficient) ||
              !std::isfinite(target[term.key.output_component]))
            throw Error("immersed operator affine row is non-finite");
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
          auto &captured = background ? result.background_affine_donor_terms
                                      : result.affine_donor_terms;
          captured.push_back({&term, effective_coefficient});
#endif
        }
      };
  evaluate_donor_terms(background_pressure_plan.donor_terms, false,
                       result.background.pressure, true);
  evaluate_donor_terms(row.background_viscous_terms, true,
                       result.background.viscous, true);
  evaluate_donor_terms(pressure_plan.donor_terms, false,
                       result.residual.pressure, false);
  evaluate_donor_terms(row.viscous_terms, true, result.residual.viscous,
                       false);
  const auto evaluate_wall_gradient_terms =
      [&](const std::vector<AffineWallGradientTerm> &terms,
          std::array<double, 3> &target) {
        if (wall_normal_gradients == nullptr && !terms.empty())
          throw Error("immersed operator affine wall input is missing");
        for (const auto &term : terms) {
          if (wall_normal_gradients == nullptr ||
              term.link >= wall_normal_gradients->size())
            throw Error("immersed operator affine wall link is invalid");
          target[term.output_component] +=
              term.coefficient * (*wall_normal_gradients)[term.link];
          if (!std::isfinite(target[term.output_component]))
            throw Error("immersed operator affine wall row is non-finite");
        }
      };
  if (wall_normal_gradients != nullptr) {
    evaluate_wall_gradient_terms(background_pressure_plan.wall_gradient_terms,
                                 result.background.pressure);
    evaluate_wall_gradient_terms(pressure_plan.wall_gradient_terms,
                                 result.residual.pressure);
  } else {
    if (!background_pressure_plan.wall_gradient_terms.empty() ||
        !pressure_plan.wall_gradient_terms.empty())
      throw Error("immersed operator unconstrained wall row is invalid");
  }
  result.removed_background = result.background;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  result.affine_wall_gradient_terms = &pressure_plan.wall_gradient_terms;
  result.background_affine_wall_gradient_terms =
      &background_pressure_plan.wall_gradient_terms;
#endif
  result.affine_plan_fingerprint = row.affine_plan_fingerprint;
  result.evaluated_group_count = row.evaluation_groups.size();
  result.simultaneous_substitution_count = 1U;
  result.canonical_affine_row_evaluation_count = 1U;
  result.link_local_runtime_evaluation_count = 0U;
  result.immutable_input_snapshot_count = 1U;
  result.background_functional_evaluation_count = 1U;
  result.background_removal_count = 1U;
  result.final_row_write_count = 1U;
  return result;
}

} // namespace

ImmersedOperatorAdapter ImmersedOperatorAdapter::create(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const immersed::ImmersedDomain &domain,
    const immersed::GhostStencilPlan &ghost_plan,
    const immersed::LocalFlowPatternTransform &transform,
    const ImmersedReconstruction &reconstruction) {
  geometry.require_compatible(topology);
  reconstruction.require_immersed_operator_compatible(topology, geometry,
                                                      domain, ghost_plan);
  auto impl = std::make_unique<Impl>();
  impl->topology = &topology;
  impl->ghost_plan = &ghost_plan;
  impl->mpi = &reconstruction.mpi_for_immersed_operator();
  impl->owned_box = topology.owned_global_box();
  impl->local_extent = {impl->owned_box.end.x - impl->owned_box.begin.x,
                        impl->owned_box.end.y - impl->owned_box.begin.y,
                        impl->owned_box.end.z - impl->owned_box.begin.z};
  impl->global_extent = topology.global_extent();
  impl->local_active.resize(topology.local_cell_count(), 0U);
  for (LocalCellId cell = 0U; cell < topology.local_cell_count(); ++cell)
    impl->local_active[cell] =
        domain.region(cell) == CellRegion::fluid ? 1U : 0U;

  std::vector<CellPairFace> pair_faces;
  pair_faces.reserve(topology.local_face_count());
  std::vector<std::vector<LocalFaceId>> topology_faces_by_cell(
      topology.local_cell_count());
  for (LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    if (owner >= topology_faces_by_cell.size())
      throw Error("immersed operator face owner is out of range");
    topology_faces_by_cell[owner].push_back(face);
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value())
      continue;
    if (*neighbour >= topology_faces_by_cell.size())
      throw Error("immersed operator face neighbour is out of range");
    if (*neighbour != owner)
      topology_faces_by_cell[*neighbour].push_back(face);
    const auto owner_id = topology.global_cell_id(owner);
    const auto neighbour_id = topology.global_cell_id(*neighbour);
    pair_faces.push_back({std::min(owner_id, neighbour_id),
                          std::max(owner_id, neighbour_id), face});
  }
  std::sort(pair_faces.begin(), pair_faces.end(),
            [](const auto &lhs, const auto &rhs) {
              return std::array{lhs.first, lhs.second, lhs.face} <
                     std::array{rhs.first, rhs.second, rhs.face};
            });

  std::vector<std::vector<std::size_t>> shared_faces_by_cell(
      topology.local_cell_count());
  for (LocalFaceId face = 0U; face < topology.local_face_count(); ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!neighbour.has_value() ||
        impl->local_active[topology.owner(face)] == 0U ||
        impl->local_active[*neighbour] == 0U)
      continue;
    SharedFace record{};
    record.face = face;
    record.owner = topology.owner(face);
    record.neighbour = *neighbour;
    record.owner_owned =
        topology.cell_ownership(record.owner) == EntityOwnership::owned;
    record.neighbour_owned =
        topology.cell_ownership(record.neighbour) == EntityOwnership::owned &&
        !topology.periodic_pair(face).has_value();
    record.area = geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    record.displacement = geometry.face_displacement_m(face);
    record.owner_center = geometry.cell_center_m(record.owner);
    record.face_center = geometry.face_center_m(face);
    if (!finite(record.area) || !finite(record.displacement) ||
        !finite(record.owner_center) || !finite(record.face_center))
      throw Error("immersed operator shared-face geometry is non-finite");
    const auto shared_index = impl->shared_faces.size();
    impl->shared_faces.push_back(record);
    if (record.owner >= shared_faces_by_cell.size() ||
        record.neighbour >= shared_faces_by_cell.size())
      throw Error("immersed operator shared-face cell is out of range");
    shared_faces_by_cell[record.owner].push_back(shared_index);
    if (record.neighbour != record.owner)
      shared_faces_by_cell[record.neighbour].push_back(shared_index);
  }

  std::vector<std::vector<std::size_t>> row_links(topology.owned_cell_count());
  for (std::size_t link_index = 0U;
       link_index < ghost_plan.immersed_operator_link_count(); ++link_index) {
    const auto &link = ghost_plan.link_for_immersed_operator(
        static_cast<immersed::ImmersedLinkId>(link_index));
    const auto fluid = topology.find_local_cell(link.fluid_cell);
    if (!fluid.has_value() ||
        topology.cell_ownership(*fluid) != EntityOwnership::owned)
      continue;
    if (impl->local_active[*fluid] == 0U)
      throw Error("immersed operator link fluid row is inactive");
    const auto key = std::array{std::min(link.fluid_cell, link.solid_cell),
                                std::max(link.fluid_cell, link.solid_cell)};
    const auto found = std::lower_bound(
        pair_faces.begin(), pair_faces.end(), key,
        [](const CellPairFace &candidate, const auto &target) {
          return std::array{candidate.first, candidate.second} < target;
        });
    if (found == pair_faces.end() ||
        std::array{found->first, found->second} != key)
      throw Error("immersed operator link face is missing");
    const auto nx = static_cast<std::uint64_t>(impl->global_extent.x);
    const auto ny = static_cast<std::uint64_t>(impl->global_extent.y);
    const auto plane = nx * ny;
    const Int3 solid{static_cast<int>(link.solid_cell % nx),
                     static_cast<int>((link.solid_cell / nx) % ny),
                     static_cast<int>(link.solid_cell / plane)};
    const auto occurrence =
        direction_occurrence(topology.global_cell(*fluid), solid);
    Real3 area =
        geometry.face_area_vector_m2(found->face, mesh::FaceSide::owner);
    if (topology.owner(found->face) != *fluid)
      area = multiply(-1.0, area);
    const auto solid_local = topology.find_local_cell(link.solid_cell);
    if (!solid_local.has_value())
      throw Error("immersed operator link solid cell is not local");
    const Real3 fluid_center = geometry.cell_center_m(*fluid);
    const Real3 solid_center = geometry.cell_center_m(*solid_local);
    const Real3 face_center = geometry.face_center_m(found->face);
    const double center_distance = norm(subtract(solid_center, fluid_center));
    const double owner_to_face = norm(subtract(face_center, fluid_center));
    if (!(owner_to_face > 0.0) || !(owner_to_face < center_distance))
      throw Error("immersed operator wall interpolation distance is invalid");
    immersed::LocalCoefficientRow background{};
    background.neighbour[occurrence] = owner_to_face / center_distance;
    background.diagonal = 1.0 - background.neighbour[occurrence];
    const auto transformed = transform.transform_full(
        background, link.fluid_to_wall_fraction, link.solid_to_fluid_normal);
    const double signed_wall_measure = -dot(area, link.solid_to_fluid_normal);
    if (!(signed_wall_measure > 0.0) || !std::isfinite(signed_wall_measure))
      throw Error("immersed operator signed wall measure is invalid");
    std::array<Real3, 7> local_sample_points{};
    const Int3 fluid_global = topology.global_cell(*fluid);
    const immersed::detail::PeriodicCellMapper periodic_cells(
        topology.global_extent(), topology.periodicity(), geometry.length_m());
    for (std::size_t sample = 0U; sample < 6U; ++sample) {
      const auto offset = occurrence_offset(sample);
      const Int3 logical{fluid_global.x + offset.x, fluid_global.y + offset.y,
                         fluid_global.z + offset.z};
      const auto image = periodic_cells.image(logical);
      if (!image.has_value())
        throw Error("immersed operator LFP sample cell is outside the mesh");
      const auto sample_id = topology.global_cell_id(image->canonical);
      const auto sample_local = topology.find_local_cell(sample_id);
      if (!sample_local.has_value())
        throw Error("immersed operator local LFP sample cell is unavailable");
      local_sample_points[sample] =
          add(geometry.cell_center_m(*sample_local), image->shift_m);
    }
    local_sample_points[6] = fluid_center;
    const auto surface_measure =
        ghost_plan.surface_measure_vector_m2_for_immersed_operator(link.id);
    const auto surface_centroid =
        ghost_plan.surface_patch_centroid_m_for_immersed_operator(link.id);
    if (!finite(surface_measure))
      throw Error("immersed operator link surface measure is non-finite");
    if (!finite(surface_centroid))
      throw Error("immersed operator link surface centroid is non-finite");
    const std::size_t index = impl->wall_links.size();
    impl->wall_links.push_back(
        {link.id, *fluid, found->face, occurrence, area, signed_wall_measure,
         link.wall_intercept_m, link.fluid_to_wall_fraction,
         link.solid_to_fluid_normal, face_center, surface_measure,
         surface_centroid, background, transformed, local_sample_points});
    row_links[*fluid].push_back(index);
  }
  impl->last_report.row_fingerprint = UINT64_C(1469598103934665603);
  hash_u64(impl->last_report.row_fingerprint,
           UINT64_C(0x48554e444c465031));
  double replacement_coefficient_square_sum = 0.0;
  for (LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell) {
    if (impl->local_active[cell] == 0U)
      continue;
    auto &links = row_links[cell];
    std::sort(links.begin(), links.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                return impl->wall_links[lhs].id < impl->wall_links[rhs].id;
              });
    std::vector<PhysicalTerm> physical_terms;
    std::vector<BoundaryReplacementTerm> replacement_terms;
    std::vector<BoundaryEvaluationGroup> evaluation_groups;
    std::optional<immersed::QuadraticReconstruction> row_reconstruction;
    PressureAffinePlan background_pressure_unconstrained;
    PressureAffinePlan background_pressure_constrained;
    std::vector<AffineDonorTerm> background_viscous_terms;
    PressureAffinePlan pressure_unconstrained;
    PressureAffinePlan pressure_constrained;
    PressureAffinePlan pressure_constrained_force;
    std::vector<AffineDonorTerm> viscous_terms;
    std::vector<AffineInputValueKey> affine_input_keys;
    std::uint64_t affine_fingerprint = 0U;
    std::uint64_t fingerprint = 0U;
    if (!links.empty()) {
      row_reconstruction.emplace(build_boundary_row_reconstruction(
          cell, links, impl->wall_links, impl->shared_faces, topology,
          geometry, ghost_plan, shared_faces_by_cell[cell]));
      std::vector<immersed::ImmersedLinkId> link_ids;
      link_ids.reserve(links.size());
      for (std::size_t slot = 0U; slot < links.size(); ++slot) {
        const auto &link = impl->wall_links[links[slot]];
        link_ids.push_back(link.id);
        for (std::size_t occurrence = 0U; occurrence < 6U; ++occurrence) {
          const double difference =
              link.transformed_row.neighbour[occurrence] -
              link.background_row.neighbour[occurrence];
          replacement_coefficient_square_sum += difference * difference;
        }
        const double diagonal_difference =
            link.transformed_row.diagonal - link.background_row.diagonal;
        const double source_difference =
            link.transformed_row.source - link.background_row.source;
        replacement_coefficient_square_sum +=
            diagonal_difference * diagonal_difference +
            source_difference * source_difference;
        const double area[3]{link.area_from_fluid.x,
                             link.area_from_fluid.y,
                             link.area_from_fluid.z};
        for (std::uint32_t component = 0U; component < 3U; ++component) {
          const auto append_physical = [&](PhysicalTermKind kind,
                                           double coefficient) {
            const auto stable_id = physical_term_id(
                topology.global_cell_id(cell), link.id, kind,
                static_cast<std::uint32_t>(link.occurrence), component);
            std::vector<std::uint64_t> source_ids;
            if (coefficient != 0.0)
              source_ids.push_back(stable_id);
            physical_terms.push_back(
                {stable_id, link.id, kind,
                 static_cast<std::uint32_t>(link.occurrence), component,
                 coefficient, std::move(source_ids)});
          };
          append_physical(PhysicalTermKind::convective_direct, 0.0);
          append_physical(PhysicalTermKind::pressure_direct, area[component]);
          append_physical(PhysicalTermKind::viscous_orthogonal, 1.0);
          append_physical(PhysicalTermKind::viscous_deferred_gradient, 1.0);
          const auto append_replacement = [&](BoundaryReplacementKind kind,
                                              std::uint32_t occurrence,
                                              double scale) {
            replacement_terms.push_back(
                {replacement_term_id(topology.global_cell_id(cell), link.id,
                                     kind, occurrence, component),
                 link.id, kind, occurrence, component, scale});
          };
          append_replacement(BoundaryReplacementKind::pressure_face, 7U,
                             area[component]);
          append_replacement(BoundaryReplacementKind::pressure_diagonal_defect,
                             6U,
                             area[component] *
                                 (link.transformed_row.diagonal -
                                  link.background_row.diagonal));
          for (std::uint32_t occurrence = 0U; occurrence < 6U; ++occurrence)
            append_replacement(
                BoundaryReplacementKind::pressure_neighbour_defect,
                occurrence,
                area[component] *
                    (link.transformed_row.neighbour[occurrence] -
                     link.background_row.neighbour[occurrence]));
          append_replacement(BoundaryReplacementKind::viscous_wall, 7U,
                             -1.0);
          append_replacement(BoundaryReplacementKind::viscous_diagonal_defect,
                             6U,
                             -(link.transformed_row.diagonal -
                               link.background_row.diagonal));
          for (std::uint32_t occurrence = 0U; occurrence < 6U; ++occurrence) {
            append_replacement(
                BoundaryReplacementKind::viscous_neighbour_defect, occurrence,
                -(link.transformed_row.neighbour[occurrence] -
                  link.background_row.neighbour[occurrence]));
          }
        }
      }
      std::sort(physical_terms.begin(), physical_terms.end(),
                [](const auto &left, const auto &right) {
                  return left.stable_id < right.stable_id;
                });
      if (std::adjacent_find(physical_terms.begin(), physical_terms.end(),
                             [](const auto &left, const auto &right) {
                               return left.stable_id == right.stable_id;
                             }) != physical_terms.end())
        throw Error("immersed operator physical term ID is not unique");
      std::sort(replacement_terms.begin(), replacement_terms.end(),
                [](const auto &left, const auto &right) {
                  return left.stable_id < right.stable_id;
                });
      if (std::adjacent_find(replacement_terms.begin(), replacement_terms.end(),
                             [](const auto &left, const auto &right) {
                               return left.stable_id == right.stable_id;
                             }) != replacement_terms.end())
        throw Error("immersed operator replacement term ID is not unique");
      BoundaryEvaluationGroup group{};
      group.links = link_ids;
      group.physical_term_indices.resize(physical_terms.size());
      std::iota(group.physical_term_indices.begin(),
                group.physical_term_indices.end(), std::size_t{0});
      group.replacement_term_indices.resize(replacement_terms.size());
      std::iota(group.replacement_term_indices.begin(),
                group.replacement_term_indices.end(), std::size_t{0});
      group.stable_id = boundary_group_id(
          topology.global_cell_id(cell), group.links, physical_terms,
          group.physical_term_indices, replacement_terms,
          group.replacement_term_indices);
      evaluation_groups.push_back(std::move(group));
      validate_boundary_evaluation_groups(
          physical_terms, replacement_terms, evaluation_groups);
      background_pressure_unconstrained =
          build_background_pressure_affine_plan(
              topology.global_cell_id(cell), links, impl->wall_links,
              ghost_plan);
      background_pressure_constrained = background_pressure_unconstrained;
      background_viscous_terms = build_background_viscous_affine_plan(
          topology.global_cell_id(cell), links, impl->wall_links, ghost_plan);
      validate_background_source_coverage(
          physical_terms, background_pressure_unconstrained,
          background_pressure_constrained, background_viscous_terms);
      pressure_unconstrained = build_pressure_affine_plan(
          topology.global_cell_id(cell), links, impl->wall_links, ghost_plan,
          false);
      pressure_constrained = build_complete_pressure_boundary_row_plan(
          cell, links, impl->wall_links, impl->shared_faces, topology,
          geometry, *row_reconstruction, ghost_plan, false, replacement_terms,
          topology_faces_by_cell[cell], shared_faces_by_cell[cell]);
      pressure_constrained_force = build_complete_pressure_boundary_row_plan(
          cell, links, impl->wall_links, impl->shared_faces, topology,
          geometry, *row_reconstruction, ghost_plan, true, replacement_terms,
          topology_faces_by_cell[cell], shared_faces_by_cell[cell]);
      viscous_terms = build_viscous_affine_plan(
          topology.global_cell_id(cell), links, impl->wall_links, ghost_plan);
      validate_affine_source_coverage(replacement_terms,
                                      pressure_unconstrained,
                                      pressure_constrained, viscous_terms);
      validate_affine_source_coverage(replacement_terms,
                                      pressure_unconstrained,
                                      pressure_constrained_force, viscous_terms);
      affine_input_keys = bind_affine_input_snapshot_indices(
          background_pressure_unconstrained,
          background_pressure_constrained, background_viscous_terms,
          pressure_unconstrained, pressure_constrained,
          pressure_constrained_force, viscous_terms);
      affine_fingerprint = affine_plan_fingerprint(
          topology.global_cell_id(cell), background_pressure_unconstrained,
          background_pressure_constrained, background_viscous_terms,
          pressure_unconstrained, pressure_constrained,
          pressure_constrained_force, viscous_terms);
      fingerprint = row_fingerprint(topology.global_cell_id(cell), links,
                                    impl->wall_links, physical_terms,
                                    replacement_terms, evaluation_groups,
                                    affine_fingerprint);
      impl->last_report.replacement_group_count +=
          static_cast<std::uint64_t>(evaluation_groups.size());
      impl->last_report.algebraic_occurrence_count +=
          static_cast<std::uint64_t>(replacement_terms.size());
      hash_u64(impl->last_report.row_fingerprint,
               topology.global_cell_id(cell));
      hash_u64(impl->last_report.row_fingerprint, fingerprint);
      ++impl->last_report.simultaneous_substitution_count;
    }
    impl->active_rows.push_back({cell, links, std::move(row_reconstruction),
                                 std::move(affine_input_keys),
                                 std::move(physical_terms),
                                 std::move(replacement_terms),
                                 std::move(evaluation_groups),
                                 std::move(background_pressure_unconstrained),
                                 std::move(background_pressure_constrained),
                                 std::move(background_viscous_terms),
                                 std::move(pressure_unconstrained),
                                 std::move(pressure_constrained),
                                 std::move(pressure_constrained_force),
                                 std::move(viscous_terms), affine_fingerprint,
                                 fingerprint});
  }
  constexpr auto missing_boundary_row =
      std::numeric_limits<std::size_t>::max();
  impl->boundary_row_by_wall_link.assign(
      ghost_plan.immersed_operator_link_count(), missing_boundary_row);
  for (std::size_t row_index = 0U; row_index < impl->active_rows.size();
       ++row_index) {
    for (const auto wall_link_index : impl->active_rows[row_index].wall_links) {
      if (wall_link_index >= impl->wall_links.size())
        throw Error("immersed boundary row link index is invalid");
      const auto link = impl->wall_links[wall_link_index].id;
      if (link >= impl->boundary_row_by_wall_link.size() ||
          impl->boundary_row_by_wall_link[link] != missing_boundary_row)
        throw Error("immersed boundary row link association is invalid");
      impl->boundary_row_by_wall_link[link] = row_index;
    }
  }
  impl->last_report.active_row_count = impl->active_rows.size();
  impl->last_report.replacement_coefficient_l2 =
      std::sqrt(replacement_coefficient_square_sum);
  impl->last_report.limiting_case_status =
      std::isfinite(impl->last_report.replacement_coefficient_l2) &&
              ((impl->last_report.replacement_group_count == 0U &&
                impl->last_report.algebraic_occurrence_count == 0U) ||
               (impl->last_report.replacement_group_count > 0U &&
                impl->last_report.algebraic_occurrence_count > 0U))
          ? 1U
          : 0U;
  if (impl->last_report.row_fingerprint == 0U)
    impl->last_report.row_fingerprint = 1U;
  impl->scratch.resize(topology.owned_cell_count() * 3U, 0.0);
  std::size_t max_affine_input_count = 0U;
  for (const auto &row : impl->active_rows)
    max_affine_input_count =
        std::max(max_affine_input_count, row.affine_input_keys.size());
  impl->affine_input_values.resize(max_affine_input_count);
  impl->viscosity_by_wall_link.resize(impl->wall_links.size());
  impl->pressure_wall_normal_gradient_by_link.resize(
      ghost_plan.immersed_operator_link_count());
  impl->local_pressure_wall_link.assign(
      ghost_plan.immersed_operator_link_count(), 0U);
  for (const auto &link : impl->wall_links) {
    if (link.id >= impl->local_pressure_wall_link.size())
      throw Error("immersed operator wall link ID is outside the plan");
    impl->local_pressure_wall_link[link.id] = 1U;
  }
  return ImmersedOperatorAdapter(std::move(impl));
}

ImmersedOperatorAdapter::ImmersedOperatorAdapter(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
ImmersedOperatorAdapter::~ImmersedOperatorAdapter() noexcept = default;
ImmersedOperatorAdapter::ImmersedOperatorAdapter(
    ImmersedOperatorAdapter &&) noexcept = default;

ImmersedOperatorAdapter::Impl &ImmersedOperatorAdapter::require_impl() const {
  if (!impl_)
    throw Error("immersed operator has been moved from");
  return *impl_;
}

void ImmersedOperatorAdapter::accumulate_momentum(
    const FaceMassFlux &mass_flux,
    const runtime::FaceFieldView<const double> &face_velocity,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &pressure,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FaceFieldView<const double> &dynamic_viscosity_by_face,
    const runtime::FieldView<double> &residual) const {
  const auto *wall_normal_gradients = scoped_wall_normal_gradients;
  auto &impl = require_impl();
  OperationGuard operation(impl.active);
  ImmersedResidualParts local_reaction{};
  std::uint64_t local_wall_functional_evaluation_count = 0U;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  std::vector<detail::ImmersedBoundaryRowEvaluationRecord>
      local_row_evaluations;
#endif
  prepare_collectively(impl, [&] {
    mass_flux.require_immersed_operator_compatible(*impl.topology);
    validate_face_view(impl, face_velocity, 3U);
    validate_face_view(impl, dynamic_viscosity_by_face, 1U);
    validate_cell_view(impl, velocity, 3U, true);
    validate_cell_view(impl, pressure, 1U, true);
    validate_cell_view(impl, velocity_gradient, 9U, true);
    validate_cell_view(impl, residual, 3U, false);
    if (wall_normal_gradients != nullptr) {
      std::fill(impl.pressure_wall_normal_gradient_by_link.begin(),
                impl.pressure_wall_normal_gradient_by_link.end(),
                std::numeric_limits<double>::quiet_NaN());
      for (const auto &condition : *wall_normal_gradients) {
        if (condition.link >=
                impl.pressure_wall_normal_gradient_by_link.size() ||
            impl.local_pressure_wall_link[condition.link] == 0U ||
            !std::isfinite(condition.value) ||
            std::isfinite(
                impl.pressure_wall_normal_gradient_by_link[condition.link]))
          throw Error("immersed operator pressure wall condition is invalid");
        impl.pressure_wall_normal_gradient_by_link[condition.link] =
            condition.value;
      }
      for (const auto &link : impl.wall_links)
        if (!std::isfinite(impl.pressure_wall_normal_gradient_by_link[link.id]))
          throw Error("immersed operator pressure wall condition is missing");
    }

    for (LocalCellId cell = 0U; cell < impl.topology->owned_cell_count();
         ++cell) {
      const auto index = impl.field_index(cell, residual.ghost_width());
      for (int component = 0; component < 3; ++component) {
        const double value = residual(index.x, index.y, index.z, component);
        if (!std::isfinite(value))
          throw Error("immersed operator residual input is non-finite");
        impl.scratch[cell * 3U + static_cast<std::size_t>(component)] = value;
      }
    }
    for (const auto &face : impl.shared_faces) {
      const double flux = mass_flux.value_for_immersed_operator(face.face);
      const double mu = dynamic_viscosity_by_face(face.face, 0);
      const auto traction =
          impl.viscous_traction(face, velocity, velocity_gradient, mu);
      const double distance = norm(face.displacement);
      const double a = norm(subtract(face.face_center, face.owner_center));
      const double b = distance - a;
      if (!(b > 0.0))
        throw Error("immersed operator pressure distance is invalid");
      const double face_pressure =
          (b / distance) * impl.read_cell(pressure, face.owner, 0) +
          (a / distance) * impl.read_cell(pressure, face.neighbour, 0);
      for (int component = 0; component < 3; ++component) {
        const double transported = face_velocity(face.face, component);
        const double convective = flux * transported;
        const double area = component == 0   ? face.area.x
                            : component == 1 ? face.area.y
                                             : face.area.z;
        const double pressure_part = face_pressure * area;
        const double total = convective + pressure_part -
                             traction[static_cast<std::size_t>(component)];
        if (!std::isfinite(transported) || !std::isfinite(total))
          throw Error("immersed operator momentum face value is non-finite");
        if (face.owner_owned)
          impl.add_owned(face.owner, component, total);
        if (face.neighbour_owned)
          impl.add_owned(face.neighbour, component, -total);
      }
    }

    for (const auto &row : impl.active_rows) {
      if (row.wall_links.empty())
        continue;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
      detail::ImmersedBoundaryRowEvaluationRecord evaluation_snapshot{};
      evaluation_snapshot.active_cell = impl.topology->global_cell_id(row.cell);
      evaluation_snapshot.row_fingerprint = row.fingerprint;
      for (std::size_t component = 0U; component < 3U; ++component)
        evaluation_snapshot.residual_before_wall[component] =
            impl.scratch[row.cell * 3U + component];
#endif
      for (const auto link_index : row.wall_links) {
        const auto &link = impl.wall_links[link_index];
        const double wall_flux =
            mass_flux.value_for_immersed_operator(link.face);
        if (!bits_zero(wall_flux))
          throw Error("immersed operator wall mass flux must be positive zero");
      }
      const auto evaluated = evaluate_boundary_row_once(
          impl, row, velocity, pressure, dynamic_viscosity_by_face,
          wall_normal_gradients == nullptr
              ? nullptr
              : &impl.pressure_wall_normal_gradient_by_link);
      ++local_wall_functional_evaluation_count;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
      const auto reaction_before_row = local_reaction;
#endif
      for (std::size_t component = 0U; component < 3U; ++component) {
        const double background =
            evaluated.background.convective[component] +
            evaluated.background.pressure[component] +
            evaluated.background.viscous[component];
        const double removed =
            evaluated.removed_background.convective[component] +
            evaluated.removed_background.pressure[component] +
            evaluated.removed_background.viscous[component];
        const double candidate =
            (impl.scratch[row.cell * 3U + component] + background) - removed +
            evaluated.residual.pressure[component] +
            evaluated.residual.viscous[component];
        if (!std::isfinite(candidate))
          throw Error("immersed operator final boundary row is non-finite");
        impl.scratch[row.cell * 3U + component] = candidate;
        local_reaction.pressure[component] -=
            evaluated.residual.pressure[component];
        local_reaction.viscous[component] -=
            evaluated.residual.viscous[component];
      }
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
      evaluation_snapshot.wall_contribution.convective =
          evaluated.residual.convective;
      evaluation_snapshot.background_contribution.convective =
          evaluated.background.convective;
      evaluation_snapshot.background_contribution.pressure =
          evaluated.background.pressure;
      evaluation_snapshot.background_contribution.viscous =
          evaluated.background.viscous;
      evaluation_snapshot.removed_background_contribution.convective =
          evaluated.removed_background.convective;
      evaluation_snapshot.removed_background_contribution.pressure =
          evaluated.removed_background.pressure;
      evaluation_snapshot.removed_background_contribution.viscous =
          evaluated.removed_background.viscous;
      evaluation_snapshot.wall_contribution.pressure =
          evaluated.residual.pressure;
      evaluation_snapshot.wall_contribution.viscous =
          evaluated.residual.viscous;
      for (std::size_t component = 0U; component < 3U; ++component) {
        evaluation_snapshot.budget_reaction_delta.convective[component] =
            local_reaction.convective[component] -
            reaction_before_row.convective[component];
        evaluation_snapshot.budget_reaction_delta.pressure[component] =
            local_reaction.pressure[component] -
            reaction_before_row.pressure[component];
        evaluation_snapshot.budget_reaction_delta.viscous[component] =
            local_reaction.viscous[component] -
            reaction_before_row.viscous[component];
      }
      const auto authority_affine_kind = [](AffineInputKind kind) {
        return static_cast<std::uint8_t>(kind);
      };
      evaluation_snapshot.affine_plan_fingerprint =
          evaluated.affine_plan_fingerprint;
      const auto append_affine_terms =
              [&](const std::vector<BoundaryRowEvaluation::
                                    EvaluatedAffineDonorTerm> &terms,
              std::vector<detail::ImmersedBoundaryAffineDonorTermRecord>
                  &target) {
            target.reserve(terms.size());
            for (const auto &term : terms) {
              std::vector<std::uint64_t> source_term_ids;
              for (const auto &coefficient : term.descriptor->link_coefficients)
                source_term_ids.insert(source_term_ids.end(),
                                       coefficient.source_term_ids.begin(),
                                       coefficient.source_term_ids.end());
              std::sort(source_term_ids.begin(), source_term_ids.end());
              source_term_ids.erase(
                  std::unique(source_term_ids.begin(), source_term_ids.end()),
                  source_term_ids.end());
              target.push_back(
                  {term.descriptor->stable_id,
                   authority_affine_kind(term.descriptor->key.input_kind),
                   term.descriptor->key.donor,
                   term.descriptor->key.input_component,
                   term.descriptor->key.output_component,
                   term.effective_coefficient,
                   static_cast<std::uint32_t>(
                       term.descriptor->link_coefficients.size()),
                   std::move(source_term_ids)});
            }
          };
      append_affine_terms(evaluated.background_affine_donor_terms,
                          evaluation_snapshot.background_affine_donor_terms);
      std::vector<std::uint64_t> group_by_term(row.replacement_terms.size(),
                                               0U);
      for (const auto &group : row.evaluation_groups)
        for (const auto term_index : group.replacement_term_indices) {
          if (term_index >= group_by_term.size() ||
              group_by_term[term_index] != 0U)
            throw Error("immersed operator evaluated term group is invalid");
          group_by_term[term_index] = group.stable_id;
        }
      evaluation_snapshot.replacement_terms.reserve(
          row.replacement_terms.size());
      for (std::size_t term_index = 0U;
           term_index < row.replacement_terms.size(); ++term_index) {
        const auto &term = row.replacement_terms[term_index];
        const auto authority_kind =
            static_cast<std::uint8_t>(term.kind);
        evaluation_snapshot.replacement_terms.push_back(
            {term.stable_id, term.link, authority_kind, term.occurrence,
             term.component, term.scale, group_by_term[term_index]});
      }
      append_affine_terms(evaluated.affine_donor_terms,
                          evaluation_snapshot.affine_donor_terms);
      if (evaluated.background_affine_wall_gradient_terms == nullptr)
        throw Error("immersed operator background wall plan is missing");
      evaluation_snapshot.background_affine_wall_gradient_terms.reserve(
          evaluated.background_affine_wall_gradient_terms->size());
      for (const auto &term :
           *evaluated.background_affine_wall_gradient_terms)
        evaluation_snapshot.background_affine_wall_gradient_terms.push_back(
            {term.stable_id, term.link, term.output_component,
             term.coefficient, term.source_term_ids});
      if (evaluated.affine_wall_gradient_terms == nullptr)
        throw Error("immersed operator affine wall plan is missing");
      evaluation_snapshot.affine_wall_gradient_terms.reserve(
          evaluated.affine_wall_gradient_terms->size());
      for (const auto &term : *evaluated.affine_wall_gradient_terms)
        evaluation_snapshot.affine_wall_gradient_terms.push_back(
            {term.stable_id, term.link, term.output_component,
             term.coefficient, term.source_term_ids});
      evaluation_snapshot.canonical_affine_row_evaluation_count =
          evaluated.canonical_affine_row_evaluation_count;
      evaluation_snapshot.link_local_runtime_evaluation_count =
          evaluated.link_local_runtime_evaluation_count;
      for (std::size_t component = 0U; component < 3U; ++component)
        evaluation_snapshot.residual_after_wall[component] =
            impl.scratch[row.cell * 3U + component];
      evaluation_snapshot.evaluated_group_count =
          evaluated.evaluated_group_count;
      evaluation_snapshot.simultaneous_substitution_count =
          evaluated.simultaneous_substitution_count;
      evaluation_snapshot.immutable_input_snapshot_count =
          evaluated.immutable_input_snapshot_count;
      evaluation_snapshot.background_functional_evaluation_count =
          evaluated.background_functional_evaluation_count;
      evaluation_snapshot.background_removal_count =
          evaluated.background_removal_count;
      evaluation_snapshot.final_row_write_count =
          evaluated.final_row_write_count;
      local_row_evaluations.push_back(std::move(evaluation_snapshot));
#endif
    }
  });

  for (LocalCellId cell = 0U; cell < impl.topology->owned_cell_count();
       ++cell) {
    const auto index = impl.field_index(cell, residual.ghost_width());
    for (int component = 0; component < 3; ++component)
      residual(index.x, index.y, index.z, component) =
          impl.scratch[cell * 3U + static_cast<std::size_t>(component)];
  }
  auto next = impl.last_report;
  next.budget_reaction_N = local_reaction;
  impl.last_report = next;
  impl.last_wall_functional_evaluation_count =
      local_wall_functional_evaluation_count;
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  impl.last_boundary_row_evaluations = std::move(local_row_evaluations);
#endif
}

void ImmersedOperatorAdapter::accumulate_transport(
    FiniteVolumeQuantity quantity, const FaceMassFlux &mass_flux,
    const runtime::FaceFieldView<const double> &face_values,
    const runtime::FieldView<const double> &values,
    const runtime::FieldView<const double> &gradients,
    const runtime::FaceFieldView<const double> &gamma_by_face,
    const runtime::FieldView<double> &residual) const {
  auto &impl = require_impl();
  OperationGuard operation(impl.active);
  prepare_collectively(impl, [&] {
    switch (quantity.kind) {
    case FiniteVolumeQuantityKind::density:
    case FiniteVolumeQuantityKind::enthalpy:
    case FiniteVolumeQuantityKind::scalar:
      break;
    case FiniteVolumeQuantityKind::velocity:
    case FiniteVolumeQuantityKind::pressure:
    default:
      throw Error("immersed operator transport quantity is invalid");
    }
    if (quantity.kind != FiniteVolumeQuantityKind::scalar &&
        quantity.scalar_index != 0U)
      throw Error("immersed operator transport quantity is invalid");
    mass_flux.require_immersed_operator_compatible(*impl.topology);
    validate_face_view(impl, face_values, 1U);
    validate_face_view(impl, gamma_by_face, 1U);
    validate_cell_view(impl, values, 1U, true);
    validate_cell_view(impl, gradients, 3U, true);
    validate_cell_view(impl, residual, 1U, false);
    for (LocalCellId cell = 0U; cell < impl.topology->owned_cell_count();
         ++cell) {
      const auto index = impl.field_index(cell, residual.ghost_width());
      const double value = residual(index.x, index.y, index.z, 0);
      if (!std::isfinite(value))
        throw Error("immersed operator residual input is non-finite");
      impl.scratch[cell * 3U] = value;
    }
    for (const auto &face : impl.shared_faces) {
      const double flux = mass_flux.value_for_immersed_operator(face.face);
      const double transported = face_values(face.face, 0);
      const double diffusive = impl.scalar_diffusive_flux(
          face, values, gradients, gamma_by_face(face.face, 0));
      const double total = flux * transported - diffusive;
      if (!std::isfinite(transported) || !std::isfinite(total))
        throw Error("immersed operator transport face value is non-finite");
      if (face.owner_owned)
        impl.add_owned(face.owner, 0, total);
      if (face.neighbour_owned)
        impl.add_owned(face.neighbour, 0, -total);
    }
    for (const auto &link : impl.wall_links) {
      if (!bits_zero(mass_flux.value_for_immersed_operator(link.face)))
        throw Error("immersed operator wall mass flux must be positive zero");
      const double gamma = gamma_by_face(link.face, 0);
      if (!(gamma >= 0.0) || !std::isfinite(gamma))
        throw Error("immersed operator wall diffusion is invalid");
    }
  });
  for (LocalCellId cell = 0U; cell < impl.topology->owned_cell_count();
       ++cell) {
    const auto index = impl.field_index(cell, residual.ghost_width());
    residual(index.x, index.y, index.z, 0) = impl.scratch[cell * 3U];
  }
}

ImmersedOperatorReport ImmersedOperatorAdapter::report() const {
  return require_impl().last_report;
}

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
namespace {

using InterfacePressureDonorKey =
    std::pair<mesh::GlobalCellId, std::uint32_t>;
using InterfacePressureWallKey =
    std::pair<immersed::ImmersedLinkId, std::uint32_t>;
using InterfacePressureDonorMap =
    std::map<InterfacePressureDonorKey, double>;
using InterfacePressureWallMap =
    std::map<InterfacePressureWallKey, double>;

struct InterfacePressureMaps final {
  InterfacePressureDonorMap background;
  InterfacePressureDonorMap a22;
  InterfacePressureDonorMap legacy_unconstrained_lfp;
  InterfacePressureDonorMap difference;
  InterfacePressureWallMap background_wall;
  InterfacePressureWallMap a22_wall;
  InterfacePressureWallMap legacy_unconstrained_lfp_wall;
  InterfacePressureWallMap difference_wall;
};

void aggregate_pressure_plan(const PressureAffinePlan &plan,
                             InterfacePressureDonorMap &donors,
                             InterfacePressureWallMap &walls) {
  for (const auto &term : plan.donor_terms) {
    if (term.key.input_kind != AffineInputKind::pressure ||
        term.key.input_component != 0U || term.key.output_component >= 3U)
      throw Error("immersed interface pressure donor kind is invalid");
    double effective = 0.0;
    for (const auto &coefficient : term.link_coefficients)
      effective += coefficient.coefficient;
    auto &value = donors[{term.key.donor, term.key.output_component}];
    value += effective;
    if (!std::isfinite(effective) || !std::isfinite(value))
      throw Error("immersed interface pressure donor is non-finite");
  }
  for (const auto &term : plan.wall_gradient_terms) {
    if (term.output_component >= 3U || !std::isfinite(term.coefficient))
      throw Error("immersed interface pressure wall term is invalid");
    auto &value = walls[{term.link, term.output_component}];
    value += term.coefficient;
    if (!std::isfinite(value))
      throw Error("immersed interface pressure wall term is non-finite");
  }
}

template <class Key>
std::map<Key, double> subtract_pressure_maps(const std::map<Key, double> &left,
                                             const std::map<Key, double> &right) {
  auto result = left;
  for (const auto &[key, value] : right)
    result[key] -= value;
  for (auto entry = result.begin(); entry != result.end();) {
    if (entry->second == 0.0)
      entry = result.erase(entry);
    else
      ++entry;
  }
  return result;
}

InterfacePressureMaps
interface_pressure_maps(const BoundaryRowPlan &row) {
  if (row.wall_links.empty())
    throw Error("immersed interface pressure row has no wall link");
  InterfacePressureMaps result;
  aggregate_pressure_plan(row.background_pressure_constrained,
                          result.background, result.background_wall);
  aggregate_pressure_plan(row.pressure_constrained, result.a22,
                          result.a22_wall);
  aggregate_pressure_plan(row.pressure_unconstrained,
                          result.legacy_unconstrained_lfp,
                          result.legacy_unconstrained_lfp_wall);
  result.difference = subtract_pressure_maps(result.a22, result.background);
  result.difference_wall =
      subtract_pressure_maps(result.a22_wall, result.background_wall);
  return result;
}

InterfacePressureMaps interface_pressure_force_maps(
    const BoundaryRowPlan &row) {
  if (row.wall_links.empty())
    throw Error("immersed interface pressure force row has no wall link");
  InterfacePressureMaps result;
  aggregate_pressure_plan(row.background_pressure_constrained,
                          result.background, result.background_wall);
  aggregate_pressure_plan(row.pressure_constrained_force, result.a22,
                          result.a22_wall);
  aggregate_pressure_plan(row.pressure_unconstrained,
                          result.legacy_unconstrained_lfp,
                          result.legacy_unconstrained_lfp_wall);
  result.difference = subtract_pressure_maps(result.a22, result.background);
  result.difference_wall =
      subtract_pressure_maps(result.a22_wall, result.background_wall);
  return result;
}

} // namespace
#endif

namespace detail {

ForceAuthorityEvaluationScope::ForceAuthorityEvaluationScope() {
  if (scoped_force_authority_evaluation)
    throw Error("immersed operator force authority scope is already active");
  scoped_force_authority_evaluation = true;
}

ForceAuthorityEvaluationScope::~ForceAuthorityEvaluationScope() {
  scoped_force_authority_evaluation = false;
}

const immersed::QuadraticReconstruction &
ImmersedBoundaryAuthorityAccess::row_reconstruction(
    const ImmersedOperatorAdapter &adapter, immersed::ImmersedLinkId link) {
  auto &impl = adapter.require_impl();
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  impl.last_boundary_authority_lookup_probe_count = 0U;
#endif
  constexpr auto missing_boundary_row =
      std::numeric_limits<std::size_t>::max();
  if (link >= impl.boundary_row_by_wall_link.size())
    throw Error("immersed boundary row link is outside the operator plan");
  const auto row_index = impl.boundary_row_by_wall_link[link];
  if (row_index == missing_boundary_row)
    throw Error("immersed boundary row link is outside the operator plan");
  if (row_index >= impl.active_rows.size())
    throw Error("immersed boundary row index is invalid");
#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  impl.last_boundary_authority_lookup_probe_count = 1U;
#endif
  const auto &row = impl.active_rows[row_index];
  if (!row.row_reconstruction.has_value())
    throw Error("immersed boundary row authority is unavailable");
  return *row.row_reconstruction;
}

void accumulate_momentum_with_wall_normal_constraints(
    const ImmersedOperatorAdapter &adapter, const FaceMassFlux &mass_flux,
    const runtime::FaceFieldView<const double> &face_velocity,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &pressure,
    const std::vector<ImmersedWallNormalGradient> &wall_normal_gradients,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FaceFieldView<const double> &dynamic_viscosity_by_face,
    const runtime::FieldView<double> &residual) {
  ScopedWallNormalGradients scope(wall_normal_gradients);
  adapter.accumulate_momentum(mass_flux, face_velocity, velocity, pressure,
                              velocity_gradient, dynamic_viscosity_by_face,
                              residual);
}

} // namespace detail

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
namespace detail {

std::vector<ImmersedBoundaryRowRecord>
ImmersedBoundaryAuthorityAccess::rows(const ImmersedOperatorAdapter &adapter) {
  const auto &impl = adapter.require_impl();
  std::vector<ImmersedBoundaryRowRecord> result;
  result.reserve(impl.active_rows.size());
  for (const auto &row : impl.active_rows) {
    ImmersedBoundaryRowRecord record{};
    record.active_cell = impl.topology->global_cell_id(row.cell);
    record.row_replacement_fingerprint =
        row.wall_links.empty() ? 0U : row.fingerprint;
    record.replacement_group_count =
        static_cast<std::uint64_t>(row.evaluation_groups.size());
    std::vector<std::uint64_t> group_by_physical(row.physical_terms.size(),
                                                 0U);
    for (const auto &group : row.evaluation_groups)
      for (const auto index : group.physical_term_indices) {
        if (index >= group_by_physical.size() || group_by_physical[index] != 0U)
          throw Error("immersed operator physical snapshot group is invalid");
        group_by_physical[index] = group.stable_id;
      }
    for (std::size_t index = 0U; index < row.physical_terms.size(); ++index) {
      const auto &term = row.physical_terms[index];
      record.covered_physical_terms.push_back(
          {term.stable_id, term.link, static_cast<std::uint8_t>(term.kind),
           term.algebraic_occurrence, term.output_component, term.coefficient,
           group_by_physical[index], term.source_term_ids});
    }
    record.links.reserve(row.wall_links.size());
    for (const auto index : row.wall_links) {
      const auto &link = impl.wall_links[index];
      record.links.push_back({link.id, link.occurrence, link.normal_scale,
                              link.solid_to_fluid_normal, row.fingerprint,
                              link.wall_intercept_m, link.area_from_fluid,
                              link.signed_wall_measure_m2,
                              link.pressure_quadrature_m,
                              link.surface_measure_m2,
                              link.surface_patch_centroid_m});
    }
    result.push_back(std::move(record));
  }
  return result;
}

std::uint64_t
ImmersedBoundaryAuthorityAccess::last_wall_functional_evaluation_count(
    const ImmersedOperatorAdapter &adapter) {
  return adapter.require_impl().last_wall_functional_evaluation_count;
}

std::uint64_t
ImmersedBoundaryAuthorityAccess::last_boundary_authority_lookup_probe_count(
    const ImmersedOperatorAdapter &adapter) {
  return adapter.require_impl().last_boundary_authority_lookup_probe_count;
}

std::vector<ImmersedBoundaryRowEvaluationRecord>
ImmersedBoundaryAuthorityAccess::last_boundary_row_evaluations(
    const ImmersedOperatorAdapter &adapter) {
  return adapter.require_impl().last_boundary_row_evaluations;
}

std::vector<ImmersedInterfacePressureRowRecord>
ImmersedBoundaryAuthorityAccess::interface_pressure_rows(
    const ImmersedOperatorAdapter &adapter) {
  const auto &impl = adapter.require_impl();
  std::vector<ImmersedInterfacePressureRowRecord> result;
  const auto append_donors = [](const InterfacePressureDonorMap &source,
                                auto &target) {
    for (const auto &[key, coefficient] : source)
      if (coefficient != 0.0)
        target.push_back({key.first, key.second, coefficient});
  };
  const auto append_walls = [](const InterfacePressureWallMap &source,
                               auto &target) {
    for (const auto &[key, coefficient] : source)
      if (coefficient != 0.0)
        target.push_back({key.first, key.second, coefficient});
  };
  for (const auto &row : impl.active_rows) {
    if (row.wall_links.empty())
      continue;
    const auto maps = interface_pressure_maps(row);
    ImmersedInterfacePressureRowRecord record{};
    record.momentum_cell = impl.topology->global_cell_id(row.cell);
    record.authority_fingerprint = row.fingerprint;
    append_donors(maps.background, record.background_donor_terms);
    append_donors(maps.a22, record.a22_donor_terms);
    append_donors(maps.legacy_unconstrained_lfp,
                  record.legacy_unconstrained_lfp_donor_terms);
    append_donors(maps.difference, record.difference_donor_terms);
    append_walls(maps.background_wall, record.background_wall_terms);
    append_walls(maps.a22_wall, record.a22_wall_terms);
    append_walls(maps.legacy_unconstrained_lfp_wall,
                 record.legacy_unconstrained_lfp_wall_terms);
    append_walls(maps.difference_wall, record.difference_wall_terms);
    result.push_back(std::move(record));
  }
  return result;
}

std::vector<ImmersedInterfacePressureRowRecord>
ImmersedBoundaryAuthorityAccess::interface_pressure_force_rows(
    const ImmersedOperatorAdapter &adapter) {
  const auto &impl = adapter.require_impl();
  std::vector<ImmersedInterfacePressureRowRecord> result;
  const auto append_donors = [](const InterfacePressureDonorMap &source,
                                auto &target) {
    for (const auto &[key, coefficient] : source)
      if (coefficient != 0.0)
        target.push_back({key.first, key.second, coefficient});
  };
  const auto append_walls = [](const InterfacePressureWallMap &source,
                               auto &target) {
    for (const auto &[key, coefficient] : source)
      if (coefficient != 0.0)
        target.push_back({key.first, key.second, coefficient});
  };
  for (const auto &row : impl.active_rows) {
    if (row.wall_links.empty())
      continue;
    const auto maps = interface_pressure_force_maps(row);
    ImmersedInterfacePressureRowRecord record{};
    record.momentum_cell = impl.topology->global_cell_id(row.cell);
    record.authority_fingerprint = row.fingerprint;
    append_donors(maps.background, record.background_donor_terms);
    append_donors(maps.a22, record.a22_donor_terms);
    append_donors(maps.legacy_unconstrained_lfp,
                  record.legacy_unconstrained_lfp_donor_terms);
    append_donors(maps.difference, record.difference_donor_terms);
    append_walls(maps.background_wall, record.background_wall_terms);
    append_walls(maps.a22_wall, record.a22_wall_terms);
    append_walls(maps.legacy_unconstrained_lfp_wall,
                 record.legacy_unconstrained_lfp_wall_terms);
    append_walls(maps.difference_wall, record.difference_wall_terms);
    result.push_back(std::move(record));
  }
  return result;
}

} // namespace detail
#endif

} // namespace hundun::finite_volume
