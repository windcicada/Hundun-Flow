// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once
#include "models_esf_detail.hpp"
namespace hundun::v04::tcr::detail {
enum class Status : std::uint8_t {
  success,
  invalid_input,
  negative_discriminant,
  weak_denominator,
  no_admissible_root,
  branch_required,
  fold_unresolved,
  branch_crossing,
  statistics_unavailable,
  mapping_unavailable,
  stale_revision,
  scientific_evidence_unavailable
};
enum class Mode : std::uint8_t { off, shadow, experimental, validated };
struct AlgebraInput {
  double eta{}, rate_ratio{};
};
struct AlgebraReport {
  Status status{Status::invalid_input};
  double discriminant{}, signed_root{}, control{};
  std::array<double, 2> roots{}; // minus, plus
  std::array<bool, 2> admissible{};
  std::uint32_t admissible_count{};
  int selected_sign{};
};
// TCR control is in [0,1], but it controls mixing, not PaSR reaction
// increments.
AlgebraReport algebra(AlgebraInput, int accepted_sign = 0) noexcept;
struct Statistics {
  Status status{Status::statistics_unavailable};
  std::size_t count{};
  double scalar_mean{}, scalar_variance{}, rate_mean{},
      scalar_rate_covariance{};
  bool minimal_pair_limited{};
};
Statistics statistics(const double *scalar, const double *rates,
                      std::size_t count) noexcept;
struct MappingReport {
  Status status{Status::mapping_unavailable};
  AlgebraInput input{};
  std::uint64_t mapping_identity{};
};
class StatisticsMappingProvider {
public:
  virtual ~StatisticsMappingProvider() = default;
  virtual MappingReport map(const Statistics &) const noexcept = 0;
};
MappingReport map_statistics(const Statistics &,
                             const StatisticsMappingProvider &) noexcept;
struct ReactantMappingInput {
  portable::Revision revision{}, expected_revision{};
  std::uint64_t composition_fingerprint{}, rate_composition_fingerprint{};
  const double *mean_mass_fractions{};
  const double *molecular_weights_kg_per_kmol{};
  std::size_t species_count{};
  const std::size_t *reactant_indices{};
  std::size_t reactant_count{};
  const double
      *field_progress_rates{}; // same units and progress definition as PSR
  std::size_t field_count{};
  double psr_progress_rate{}, weak_rate_absolute_threshold{};
};
// Author theory Eq 2.80,2.88: eta=sum(X_reactants), R=mean(omega)/omega_PSR.
MappingReport
ideal_gas_reactant_mole_fraction_v1(const ReactantMappingInput &) noexcept;
struct History {
  portable::Revision revision{};
  bool initialized{};
  AlgebraInput input{};
  double signed_root{}, control{};
  int branch_sign{};
  std::uint64_t mapping_identity{};
  int initialization_sign{};
  std::uint64_t fold_count{};
  AlgebraInput last_fold{};
  portable::Revision last_fold_base_revision{};
};
// Explicit, revision-bound fold evidence on a piecewise continuation path.
struct FoldEvent {
  bool supplied{};
  portable::Revision base_revision{};
  AlgebraInput at_fold{};
  int departure_sign{};
};
struct TrialRequest {
  portable::Revision expected_revision{};
  MappingReport mapping{};
  Mode mode{Mode::off};
  int initialization_sign{};
  FoldEvent fold{};
};
struct Trial {
  Status status{Status::invalid_input};
  AlgebraReport observed{};
  History candidate{};
  portable::Revision base_revision{};
  Mode mode{Mode::off};
  bool available{}, feedback{};
  double mixer_control{1};
  TrialRequest request{};
};
bool valid_history(const History &) noexcept;
Trial prepare(const History &, const TrialRequest &) noexcept;
// Produces accepted-value candidate; neither function mutates its inputs.
Trial accept(const History &, const Trial &, portable::Revision next) noexcept;
Trial restore(const History &snapshot, portable::Revision expected) noexcept;
struct MixReport {
  Trial tcr{};
  esf::detail::Report esf{};
  bool available{};
};
MixReport mix(const History &, const TrialRequest &, esf::detail::Request,
              esf::detail::Workspace &) noexcept;
} // namespace hundun::v04::tcr::detail
