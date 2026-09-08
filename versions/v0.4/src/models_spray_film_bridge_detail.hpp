// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_spray_events_detail.hpp"
#include "models_spray_properties_detail.hpp"
namespace hundun::v04::spray::detail {
struct ParcelGasSample {
  portable::Status status{portable::Status::unavailable};
  portable::Revision revision{};
  std::uint64_t composition_fingerprint{};
  double pressure_pa{}, enthalpy_j_per_kg{};
  Vector3 velocity_m_per_s{};
  std::size_t species_count{};
};
// Geometry/field ownership stays external. The sampler supplies full Y in the
// caller-owned buffer at this position and accepted-step-relative time.
class ParcelGasStateProvider {
public:
  virtual ~ParcelGasStateProvider() = default;
  virtual ParcelGasSample sample(const SprayParcelState &,
                                 double elapsed_step_time_s, ParcelPass,
                                 portable::Revision, double *full_y,
                                 std::size_t capacity) const noexcept = 0;
};
struct FilmEnvironmentReport {
  portable::Status status{portable::Status::invalid_input};
  ParcelTransferEnvironment environment{};
};
// Exclusive prepared lane. Asset, neutral gas provider and sampler must outlive
// it. Does not own product fields or publish parcel/gas state. Copying
// forbidden because LiquidPropertyService borrows the exact asset pack.
class FilmEnvironmentBridge final : public ParcelTransferEnvironmentProvider {
public:
  FilmEnvironmentBridge(const LiquidAsset &, portable::GasQueryProvider &,
                        const ParcelGasStateProvider &, portable::Revision);
  FilmEnvironmentBridge(const FilmEnvironmentBridge &) = delete;
  FilmEnvironmentBridge &operator=(const FilmEnvironmentBridge &) = delete;
  // Reuse the cold-reserved exclusive lane for the next native attempt.
  // No sampled state is cached across queries; failed binding preserves the
  // current revision and all owned capacities.
  portable::Status bind_revision(portable::Revision revision) noexcept {
    if (revision.algorithm_version != 1 || !revision.input_revision)
      return portable::Status::invalid_input;
    revision_ = revision;
    return portable::Status::success;
  }
  FilmEnvironmentReport query(const SprayParcelState &, double, ParcelPass,
                              portable::Revision) const noexcept;
  ParcelTransferEnvironment sample(const SprayParcelState &, double, ParcelPass,
                                   portable::Revision) const noexcept override;

private:
  const LiquidAsset &asset_;
  portable::GasQueryProvider &gas_;
  const ParcelGasStateProvider &sampler_;
  portable::Revision revision_;
  LiquidPropertyService liquid_;
  mutable FilmQueryWorkspace film_;
  mutable std::vector<double> y_;
};
} // namespace hundun::v04::spray::detail
