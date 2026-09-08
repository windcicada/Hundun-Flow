// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "hundun/v04_portable.hpp"
#include "models_spray_source_detail.hpp"
#include <vector>

namespace hundun::v04::portable {
enum class ExchangeChannel : std::uint8_t {
  interphase,
  wall,
  outlet,
  external
};
struct ExchangeSegment {
  Revision revision{};
  std::uint64_t global_cell{};
  spray::ParcelId parcel_id{};
  std::uint64_t segment_ordinal{};
  double deposition_weight{1.0};
  ExchangeChannel channel{ExchangeChannel::interphase};
  spray::detail::ParcelExchangeBudget
      delta{}; // multiplicity included, weight not yet applied
  std::size_t vapor_species_index{};
};
struct ExchangeCell {
  std::uint64_t global_cell{};
  double volume_m3{};
  double gas_mass_kg{};
  spray::Vector3 gas_momentum_kg_m_per_s{};
};
struct CellExchange {
  std::uint64_t global_cell{};
  double volume_m3{};
  spray::detail::GasCellExchangeCandidate gas{};
  const double *gas_species_mass_delta_kg{};
  std::size_t species_count{};
};
struct BreakupEnergyLedger {
  double deformation_consumed_j{}, surface_increase_j{},
      bulk_kinetic_increase_j{}, residual_j{};
  std::size_t event_count{};
};
struct ExchangeBatchReport {
  Status status{Status::invalid_input};
  std::size_t failure_index{};
  std::uint64_t generation{};
  Revision revision{};
  const CellExchange *cells{};
  std::size_t cell_count{};
  ExchangeDelta wall{}, outlet{}, external{};
  bool available{};
};
// A borrowed, revision-bound certificate for rows whose complete deposition
// stencil was audited before owner routing. Invalidated by the next plan call.
class RoutedExchangeSegments {
public:
  RoutedExchangeSegments() = default;

private:
  friend class OwnerExchangeRoutingPlan;
  friend class ExchangeWorkspace;
  const ExchangeSegment *data_{};
  std::size_t size_{};
  Revision revision_{};
  const std::uint64_t *generation_{};
  std::uint64_t captured_generation_{};
};
class ExchangeWorkspace {
public:
  ExchangeWorkspace(std::size_t cell_capacity, std::size_t segment_capacity,
                    std::size_t species_count = 1);
  ExchangeWorkspace(const ExchangeWorkspace &) = delete;
  ExchangeWorkspace &operator=(const ExchangeWorkspace &) = delete;
  ExchangeBatchReport evaluate(Revision, const ExchangeCell *, std::size_t,
                               const ExchangeSegment *, std::size_t,
                               double thermal_atol_j,
                               double thermal_rtol) noexcept;
  ExchangeBatchReport evaluate_routed(Revision, const ExchangeCell *,
                                      std::size_t, RoutedExchangeSegments,
                                      double thermal_atol_j,
                                      double thermal_rtol) noexcept;
  bool current(const ExchangeBatchReport &report) const noexcept;

private:
  ExchangeBatchReport evaluate_impl(Revision, const ExchangeCell *, std::size_t,
                                    const ExchangeSegment *, std::size_t,
                                    double, double,
                                    bool complete_stencils) noexcept;
  std::vector<CellExchange> cells_;
  std::vector<std::size_t> order_;
  std::size_t species_{};
  std::vector<double> species_mass_;
  std::uint64_t generation_{};
};
} // namespace hundun::v04::portable
