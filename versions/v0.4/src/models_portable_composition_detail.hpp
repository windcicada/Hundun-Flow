// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_esf_detail.hpp"
#include "models_exchange_batch_detail.hpp"
#include "models_exchange_routing_detail.hpp"
#include "models_spray_events_detail.hpp"
#include "models_spray_migration_detail.hpp"
#include "models_spray_parcel_detail.hpp"
#include "models_spray_properties_detail.hpp"
#include "models_tcr_detail.hpp"

namespace hundun::v04::portable {
bool valid_vapor_mapping(const spray::SprayParcelState &,
                         const spray::detail::LiquidAsset &,
                         const GasIdentity &, std::size_t) noexcept;
// Prescribed cells, not a mesh/flow solver. All fields use the cell inventory
// density for common-source conversion; EOS densities are reported separately.
struct CompositionCellInput {
  ExchangeCell inventory{};
  double pressure_pa{};
  const double *fields{}; // N x (Ns+1), complete Y then total h
  esf::detail::Request transport{};
  tcr::detail::History tcr_history{};
  tcr::detail::TrialRequest tcr_trial{};
  const std::size_t *reactant_indices{};
  std::size_t reactant_count{};
  const double *progress_species_weights{}; // c=sum(w_i Y_i), source in 1/s
  double weak_progress_rate_per_s{1e-12};
};
struct CompositionInput {
  Revision revision{};
  double start_time_s{}, duration_s{};
  std::size_t field_count{};
  const combustion::ChemistryIdentity *identity{};
  GasQueryProvider *gas{};
  GasAdvanceProvider *chemistry{};
  bool reaction_enabled{true};
  const CompositionCellInput *cells{};
  std::size_t cell_count{};
  const spray::detail::ParcelEventsInput *parcels{};
  const std::size_t *vapor_species_indices{};
  const spray::detail::LiquidAsset *const *parcel_materials{};
  std::size_t parcel_count{};
  const spray::detail::InjectorValueSnapshot *candidate_injectors{};
  std::size_t injector_count{};
  std::uint64_t parcel_rng_seed{};
  spray::detail::ParcelMigrationPlan *migration{};
  const spray::detail::ParcelLocationProvider *location{};
  double thermal_absolute_tolerance_j{1e-14};
  double thermal_relative_tolerance{1e-8};
};
struct CompositionCellCandidate {
  ExchangeCell inventory{};
  const double *fields{};
  std::array<double, 4> eos_densities_kg_per_m3{};
  tcr::detail::Trial tcr{};
  double heat_release_report_j_per_m3{};
};
struct BreakupEnergyLedger {
  double deformation_consumed_j{}, surface_increase_j{},
      bulk_kinetic_increase_j{}, residual_j{};
  std::size_t event_count{};
};
struct CompositionReport {
  Status status{Status::invalid_input};
  std::string_view model_id{"portable_source_transport_reaction_v1"};
  Revision revision{};
  std::uint64_t generation{};
  std::size_t failure_index{};
  std::uint32_t failure_module{}; // 1 preflight,2 parcel,3 exchange,4 gas,5
                                  // TCR,6 ESF,7 migration
  int lowest_failing_rank{-1};
  bool available{};
  const CompositionCellCandidate *cells{};
  std::size_t cell_count{}, field_count{}, species_count{};
  const spray::detail::ParcelMigrationValue *parcels{};
  std::size_t parcel_count{};
  ExchangeDelta wall{}, outlet{};
  BreakupEnergyLedger
      breakup_energy{}; // separate internal TAB reservoir, never gas source
  std::uint64_t chemistry_calls{};
  bool reaction_enabled{};
};
class CompositionWorkspace {
public:
  CompositionWorkspace(std::size_t cells, std::size_t parcels,
                       std::size_t segments, std::size_t species);
  CompositionWorkspace(const CompositionWorkspace &) = delete;
  CompositionWorkspace &operator=(const CompositionWorkspace &) = delete;
  // All ranks participate, even when local work fails. No state is published.
  // Order is source -> transport/IEM/TCR -> two chemistry half-intervals.
  // This first-order composition is NOT labelled Strang or product acceptance.
  CompositionReport prepare(MPI_Comm, const CompositionInput &) noexcept;
  bool current(const CompositionReport &) const noexcept;

private:
  CompositionReport prepare_parcels(const CompositionInput &) noexcept;
  CompositionReport prepare_cells(const CompositionInput &,
                                  CompositionReport) noexcept;
  std::size_t species_{}, parcel_capacity_{}, segment_capacity_{};
  std::uint64_t generation_{};
  ExchangeWorkspace exchange_;
  ExchangeRoutingWorkspace routing_;
  std::size_t segment_count_{};
  std::size_t id_count_{};
  std::vector<spray::ParcelId> attempt_ids_;
  esf::detail::Workspace esf_;
  std::vector<ExchangeCell> cell_inputs_;
  std::vector<ExchangeSegment> segments_;
  std::vector<CompositionCellCandidate> cells_;
  std::vector<double> fields_, scratch_;
  std::vector<spray::detail::ParcelMigrationValue> parcels_;
  std::vector<spray::detail::ParcelEventsInput> jobs_;
  std::vector<std::size_t> vapor_indices_;
};
struct CompositionCellValue {
  ExchangeCell inventory{};
  double pressure_pa{};
  std::vector<double> fields;
  tcr::detail::History tcr_history{};
  esf::detail::CounterAddress random{};
};
struct CompositionSnapshot {
  std::uint32_t version{1};
  Revision accepted_revision{};
  GasIdentity gas_identity;
  std::size_t field_count{};
  double accepted_time_s{};
  bool reaction_enabled{};
  std::uint64_t parcel_rng_seed{};
  std::uint32_t parcel_rng_algorithm_version{1};
  std::vector<CompositionCellValue> cells;
  std::vector<spray::detail::ParcelMigrationValue> parcels;
  std::vector<spray::detail::InjectorValueSnapshot> injectors;
};
struct CompositionValueReport {
  Status status{Status::invalid_input};
  bool available{};
  CompositionSnapshot candidate;
};
// Cold ordinary-value preparation. A complete return is assigned by the caller
// as one value; no products, fields or MPI participants are committed here.
CompositionValueReport prepare_accepted_values(MPI_Comm,
                                               const CompositionWorkspace &,
                                               const CompositionInput &,
                                               const CompositionReport &,
                                               Revision next) noexcept;
CompositionValueReport restore_composition_values(
    const CompositionSnapshot &, const GasIdentity &, Revision expected,
    std::size_t max_cells, std::size_t max_parcels,
    const spray::detail::ParcelLocationProvider *) noexcept;
} // namespace hundun::v04::portable
