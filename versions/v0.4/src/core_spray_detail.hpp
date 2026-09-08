// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core_reaction_detail.hpp"
#include "core_spray_advance_detail.hpp"
#include "core_spray_history_detail.hpp"

namespace hundun::v04::detail {
// Product-owned cold services, ordinary candidate lane and native restart
// participant. Neither this owner nor the portable models commit gas fields.
class ProductSpray {
public:
  Status configure_local(const ValidatedModel &, const std::filesystem::path &,
                         const ProductReactionSources &,
                         const CartesianGeometryPlan &, MeshPatch, int rank,
                         RestartCellRecordsView tcr);
  Status configure_collective(MPI_Comm, Span<const RemoteDonorFieldSpec>);
  Status configure_source_transport(Span<const RemoteDonorFieldSpec>);
  Status exchange_source_transport(Span<FieldView>, RevisionToken boundary,
                                   Span<ThermophysicalGhostAuthority>) noexcept;
  Status prepare(portable::Revision, double duration, double pressure_reference,
                 ConstFieldView density, FieldView pressure, FieldView enthalpy,
                 FieldView velocity,
                 Span<const FieldView> independent) noexcept;
  Status stage(RestartCellRecordsView tcr) noexcept;
  Status stage_restore(const RestartImage &) noexcept;
  Status validate_restored_positions() noexcept;
  Status preflight_commit() const noexcept {
    return history.preflight_commit();
  }
  void commit() noexcept {
    if (enabled())
      history.commit();
  }
  void discard() noexcept {
    history.discard();
    if (advance_)
      advance_->discard();
  }
  bool enabled() const noexcept { return identity_ != 0; }
  PlanFingerprint fingerprint() const noexcept { return identity_; }
  std::uint64_t owned_bytes() const noexcept { return owned_bytes_; }
  RemoteDonorExchangeStats halo_stats() const noexcept {
    auto result = halo_.stats();
    const auto source = source_halo_.stats();
    result.received_cells += source.received_cells;
    result.supplied_cells += source.supplied_cells;
    result.bytes_per_exchange += source.bytes_per_exchange;
    result.peer_messages += source.peer_messages;
    return result;
  }
  portable::ExchangeBatchReport exchange() const noexcept {
    return advance_ ? advance_->exchange() : portable::ExchangeBatchReport{};
  }
  ProductSprayHistory history;

private:
  Status agree(Status) const noexcept;
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10340}; }
  MPI_Comm comm_{MPI_COMM_NULL};
  int rank_{};
  PlanFingerprint identity_{};
  std::uint64_t owned_bytes_{}, maximum_bytes_{}, local_bytes_{};
  const CartesianGeometryPlan *geometry_{};
  MeshPatch patch_{};
  SpraySpec spec_;
  spray::detail::LiquidAsset asset_;
  ProductParcelGas gas_;
  ProductParcelGeometry events_;
  std::unique_ptr<spray::detail::FilmEnvironmentBridge> film_;
  std::vector<std::unique_ptr<ProductSprayHistory::Injector>> injectors_;
  std::vector<ProductSprayHistory::Injector *> injector_views_;
  std::unique_ptr<ProductParcelAdvance> advance_;
  std::array<bool, 3> periodic_{};
  std::vector<GlobalCellId> halo_cells_;
  std::vector<Int3> halo_indices_;
  std::vector<FieldView> halo_views_;
  std::vector<ConstFieldView> species_views_;
  RemoteDonorExchangePlan halo_;
  RemoteDonorExchangePlan source_halo_;
  bool source_halo_bound_{};
  spray::detail::ParcelMigrationPlan restore_migration_;
  portable::Revision revision_{};
};
} // namespace hundun::v04::detail
