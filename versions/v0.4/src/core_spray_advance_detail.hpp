// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core_spray_events_detail.hpp"
#include "core_spray_gas_detail.hpp"
#include "models_exchange_owner_detail.hpp"
#include "models_spray_parcel_detail.hpp"
#include <memory>

namespace hundun::v04::detail {
// One native gas attempt owns the publication decision. This lane only builds
// ordinary parcel/source candidates and stages the borrowed injectors. Its gas
// sampler must already be bound to accepted MeanState with complete Halo.
class ProductParcelAdvance {
public:
  using Parcel = spray::detail::ParcelMigrationValue;
  using Injector = spray::detail::DeterministicInjector;
  ProductParcelAdvance(ProductParcelGas &gas, ProductParcelGeometry &geometry,
                       const spray::detail::LiquidAsset &asset,
                       spray::detail::FilmEnvironmentBridge &film)
      : gas_(gas), events_(geometry), asset_(asset), film_(film),
        interval_(film), tab_(film), breakup_(film, 2) {}
  Status configure(MPI_Comm, const CartesianGeometryPlan &, MeshPatch,
                   const SpraySpec &, std::uint64_t maximum_bytes) noexcept;
  Status prepare(Span<const Parcel>, Span<Injector *const>, portable::Revision,
                 double duration_s, ConstFieldView density,
                 ConstFieldView velocity) noexcept;
  Span<const Parcel> parcels() const noexcept {
    return available_ ? Span<const Parcel>{current_.data(), current_.size()}
                      : Span<const Parcel>{};
  }
  portable::ExchangeBatchReport exchange() const noexcept {
    return available_ ? exchange_ : portable::ExchangeBatchReport{};
  }
  portable::BreakupEnergyLedger breakup_energy() const noexcept {
    return available_ ? breakup_energy_ : portable::BreakupEnergyLedger{};
  }
  void discard() noexcept {
    available_ = false;
    route_.discard();
    migration_.discard();
  }
  std::uint64_t owned_bytes() const noexcept { return owned_bytes_; }
  std::uint32_t waves() const noexcept { return waves_; }

private:
  struct Job {
    Parcel value;
    double start{}, duration{};
  };
  Status agree(Status) noexcept;
  Status move() noexcept;
  Status wave(portable::Revision, double start, double duration,
              std::uint32_t index) noexcept;
  Status append(const portable::ExchangeSegment &) noexcept;
  Status audit_id(spray::ParcelId, portable::Revision) noexcept;
  Status inventory(const spray::SprayParcelState &,
                   portable::Revision) noexcept;
  Status ledger(const portable::ExchangeDelta &, const Parcel &,
                portable::Revision, portable::ExchangeChannel,
                std::uint64_t ordinal) noexcept;
  ProductParcelGas &gas_;
  ProductParcelGeometry &events_;
  const spray::detail::LiquidAsset &asset_;
  spray::detail::FilmEnvironmentBridge &film_;
  spray::detail::FixedAsParcelIntervalProvider interval_;
  spray::detail::FixedTabEvolutionProvider tab_;
  spray::detail::FixedTabEventBreakupProvider breakup_;
  MPI_Comm comm_{MPI_COMM_NULL};
  int rank_{};
  const CartesianGeometryPlan *geometry_{};
  MeshPatch patch_{};
  std::size_t parcel_capacity_{}, segment_capacity_{};
  double maximum_step_{}, minimum_step_{}, relative_tolerance_{},
      minimum_width_{};
  bool tab_enabled_{}, available_{};
  std::uint64_t owned_bytes_{}, local_anchor_{};
  std::uint32_t waves_{};
  portable::Revision revision_{};
  std::vector<Parcel> current_, next_;
  std::vector<Job> jobs_;
  std::vector<portable::ExchangeSegment> segments_, ids_;
  std::vector<portable::ExchangeCell> cells_;
  spray::detail::ParcelMigrationPlan migration_;
  portable::OwnerExchangeRoutingPlan route_;
  std::unique_ptr<portable::ExchangeWorkspace> batch_;
  portable::ExchangeBatchReport exchange_{};
  portable::BreakupEnergyLedger breakup_energy_{};
};
} // namespace hundun::v04::detail
