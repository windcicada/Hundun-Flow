// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_portable.hpp"
#include "models_spray_breakup_detail.hpp"
#include "models_spray_transfer_detail.hpp"
#include <array>

namespace hundun::v04::spray::detail {
inline constexpr std::size_t kParcelEventSegmentCapacity = 256;
inline constexpr std::size_t kParcelEventCapacity = 64;
enum class ParcelPass : std::uint8_t { predictor, corrector };
struct ParcelAuxiliaryState {
  double tab_deformation{};
  double tab_deformation_rate_per_s{};
  std::uint64_t breakup_ordinal{};
};
struct ParcelTabIntervalReport {
  portable::Revision revision{};
  TabBreakupReport tab{};
};
class ParcelTabEvolutionProvider {
public:
  virtual ~ParcelTabEvolutionProvider() = default;
  virtual ParcelTabIntervalReport
  advance(const SprayParcelState &, const ParcelAuxiliaryState &,
          double elapsed_s, double duration_s, ParcelPass,
          portable::Revision) const noexcept = 0;
};
enum class ParcelEventKind : std::uint8_t {
  complete_evaporation,
  physical_outlet,
  wall_collision,
  breakup,
  internal_cell_crossing
};
enum class ParcelEventsStatus : std::uint8_t {
  success,
  invalid_input,
  provider_failure,
  stale_revision,
  capacity_exceeded,
  minimum_step,
  non_finite,
  budget_failure,
  event_failure
};
struct ParcelIntervalReport {
  bool available{};
  portable::Revision revision{};
  SprayParcelState parcel{};
  double elapsed_duration_s{};
  double initial_liquid_absolute_enthalpy_j_per_kg{};
  double liquid_absolute_enthalpy_j_per_kg{};
  // Caloric slope at both endpoint temperatures, supplied by the liquid
  // service. Used solely to bound binary64 endpoint quantization, not LTE.
  double liquid_heat_capacity_bound_j_per_kg_k{};
  bool complete_evaporation{};
  TransferExchangeCandidate exchange{};
};
class ParcelIntervalProvider {
public:
  virtual ~ParcelIntervalProvider() = default;
  // Exact terminal time is relative to this interval. Terminal mass and
  // diameter must be exactly zero; nonterminal outputs advance the full
  // requested dt.
  virtual ParcelIntervalReport advance(const SprayParcelState &,
                                       double elapsed_s, double duration_s,
                                       ParcelPass,
                                       portable::Revision) const noexcept = 0;
};
struct ParcelEvent {
  ParcelEventKind kind{ParcelEventKind::internal_cell_crossing};
  double elapsed_time_s{};
  std::uint64_t identity{};
  std::uint64_t next_global_cell{};
  Vector3 wall_normal{};
  double restitution{1.0};
  bool applied{};
  Vector3 contact_position_m{};
  bool has_contact_position{};
};
struct ParcelEventQueryReport {
  bool available{};
  portable::Revision revision{};
  std::size_t count{};
  std::array<ParcelEvent, 8> events{};
};
class ParcelEventGeometryProvider {
public:
  virtual ~ParcelEventGeometryProvider() = default;
  virtual ParcelEventQueryReport query(const SprayParcelState &begin,
                                       const SprayParcelState &end,
                                       double begin_time, double end_time,
                                       ParcelPass,
                                       portable::Revision) const noexcept = 0;
};
class ParcelEventBreakupProvider {
public:
  virtual ~ParcelEventBreakupProvider() = default;
  virtual BreakupChildReport split(const SprayParcelState &,
                                   const ParcelEvent &,
                                   const ParcelAuxiliaryState &,
                                   const TabBreakupReport *, ParcelPass,
                                   portable::Revision) const noexcept = 0;
};
struct ParcelEventSegment {
  std::uint32_t ordinal{};
  std::uint64_t global_cell{};
  double begin_time_s{}, end_time_s{};
  SprayParcelState begin{}, end{};
  TransferExchangeCandidate exchange{}; // already extensive, multiplicity once
};
struct ParcelEventsInput {
  SprayParcelState accepted_parcel{};
  portable::Revision revision{};
  const ParcelIntervalProvider *interval{};
  const ParcelEventGeometryProvider *geometry{};
  const ParcelEventBreakupProvider *breakup{};
  ParcelAuxiliaryState accepted_auxiliary{};
  const ParcelTabEvolutionProvider
      *tab_evolution{}; // nullptr: explicit off, preserve auxiliary values
  // Provider/event/segment times are measured from the common accepted-step
  // origin. A child continuation starts at parent offset + advanced duration.
  double elapsed_offset_s{};
  double duration_s{};
  double initial_substep_s{};
  double minimum_substep_s{1e-12};
  double relative_tolerance{1e-6};
  double position_absolute_tolerance_m{1e-10};
  double velocity_absolute_tolerance_m_per_s{1e-10};
  double mass_absolute_tolerance_kg{1e-24};
  double temperature_absolute_tolerance_k{1e-8};
  double momentum_absolute_tolerance_kg_m_per_s{1e-20};
  double energy_absolute_tolerance_j{1e-16};
  double event_time_tolerance_s{1e-12};
  std::size_t maximum_segments{kParcelEventSegmentCapacity};
  std::size_t maximum_events{kParcelEventCapacity};
  std::size_t maximum_attempts{2048};
};
struct ParcelEventsReport {
  ParcelEventsStatus status{ParcelEventsStatus::invalid_input};
  std::string_view model_id{"adaptive_predictor_corrector_parcel_events_v1"};
  bool available{};
  portable::Revision revision{};
  std::uint32_t predictor_passes{}, corrector_passes{};
  std::size_t attempts{}, rejected_substeps{}, provider_queries{};
  double failure_time_s{};
  std::size_t failure_segment{};
  ParcelPass failure_pass{ParcelPass::predictor};
  double advanced_duration_s{};
  double children_remaining_duration_s{};
  SprayParcelState parcel{};
  bool parent_removed{}, complete_evaporation{}, physical_outlet{},
      breakup_requested{};
  ParcelAuxiliaryState auxiliary{};
  bool tab_evolved{};
  std::size_t segment_count{}, event_count{};
  std::array<ParcelEventSegment, kParcelEventSegmentCapacity> segments{};
  std::array<ParcelEvent, kParcelEventCapacity> events{};
  TransferExchangeCandidate exchange{};
  portable::ExchangeDelta outlet_inventory{};
  portable::ExchangeDelta wall_exchange{}; // wall gain, not gas phase exchange
  BreakupChildCandidate children{};
  BreakupConservationReport breakup_budget{};
};
[[nodiscard]] ParcelEventsReport
integrate_parcel_events(const ParcelEventsInput &) noexcept;

// Existing A--S + Schiller--Naumann fixed interval adapter. The environment
// provider supplies a freshly sampled film and absolute liquid enthalpy for
// each queried state; it can bind the neutral P4 film/asset module.
struct ParcelTransferEnvironment {
  bool available{};
  portable::Revision revision{};
  GasTransferEnvironment gas{};
  const LiquidPropertyService *liquid{};
  double liquid_absolute_enthalpy_j_per_kg{};
  double far_gas_density_kg_per_m3{}; // never substitute one-third film density
                                      // for TAB forcing
};
class ParcelTransferEnvironmentProvider {
public:
  virtual ~ParcelTransferEnvironmentProvider() = default;
  virtual ParcelTransferEnvironment
  sample(const SprayParcelState &, double, ParcelPass,
         portable::Revision) const noexcept = 0;
};
class FixedAsParcelIntervalProvider final : public ParcelIntervalProvider {
public:
  explicit FixedAsParcelIntervalProvider(
      const ParcelTransferEnvironmentProvider &provider)
      : provider_(provider) {}
  ParcelIntervalReport advance(const SprayParcelState &, double, double,
                               ParcelPass,
                               portable::Revision) const noexcept override;

private:
  const ParcelTransferEnvironmentProvider &provider_;
};
class FixedTabEvolutionProvider final : public ParcelTabEvolutionProvider {
public:
  explicit FixedTabEvolutionProvider(
      const ParcelTransferEnvironmentProvider &provider)
      : provider_(provider) {}
  ParcelTabIntervalReport advance(const SprayParcelState &,
                                  const ParcelAuxiliaryState &, double, double,
                                  ParcelPass,
                                  portable::Revision) const noexcept override;

private:
  const ParcelTransferEnvironmentProvider &provider_;
};
class FixedTabEventBreakupProvider final : public ParcelEventBreakupProvider {
public:
  FixedTabEventBreakupProvider(
      const ParcelTransferEnvironmentProvider &provider,
      std::uint32_t child_parcels)
      : provider_(provider), child_parcels_(child_parcels) {}
  BreakupChildReport split(const SprayParcelState &, const ParcelEvent &,
                           const ParcelAuxiliaryState &,
                           const TabBreakupReport *, ParcelPass,
                           portable::Revision) const noexcept override;

private:
  const ParcelTransferEnvironmentProvider &provider_;
  std::uint32_t child_parcels_{};
};
} // namespace hundun::v04::spray::detail
