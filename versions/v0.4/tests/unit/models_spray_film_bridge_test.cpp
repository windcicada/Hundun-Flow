// SPDX-License-Identifier: Apache-2.0
#include "models_chemistry_adapter_detail.hpp"
#include "models_spray_film_bridge_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04;
using namespace hundun::v04::spray::detail;
class PrescribedGas final : public ParcelGasStateProvider {
public:
  std::uint64_t fingerprint{};
  mutable unsigned queries{};
  bool fail{};
  ParcelGasSample sample(const spray::SprayParcelState &, double elapsed,
                         ParcelPass, portable::Revision revision, double *y,
                         std::size_t cap) const noexcept override {
    ++queries;
    if (fail || cap < 2)
      return {};
    y[0] = 0.01;
    y[1] = 0.99;
    return {portable::Status::success, revision, fingerprint, 100000, 101000,
            {elapsed, 0, 0},           2};
  }
};
bool check(bool ok, const char *msg) {
  if (!ok)
    std::cerr << "FAIL: " << msg << '\n';
  return ok;
}
int main(int argc, char **argv) {
  if (argc != 2)
    return 1;
  chemistry::detail::AnalyticIsomerBackend gas;
  auto asset = load_liquid_asset(argv[1], UINT64_C(6004043157121730787),
                                 gas.gas_identity());
  if (!asset.available)
    return 1;
  PrescribedGas sampler;
  sampler.fingerprint = gas.composition().fingerprint;
  portable::Revision revision{4, 9, 1};
  FilmEnvironmentBridge bridge(asset.asset, gas, sampler, revision);
  spray::SprayParcelState parcel{
      {1, 2}, {0, 0, 0}, {0, 0, 0}, 4.188790204786391e-10,
      1e-4,   1,         298.15,    asset.asset.pack.material_fingerprint,
      0,      0};
  auto environment = bridge.query(parcel, 0.5, ParcelPass::corrector, revision);
  bool ok = check(environment.environment.available,
                  "full-composition film bridge produces A-S environment");
  ok &= check(environment.environment.gas.gas_velocity_m_per_s[0] == 0.5,
              "accepted-step time coordinate reaches external sampler");
  ok &= check(environment.environment.liquid_absolute_enthalpy_j_per_kg ==
                  -100000,
              "bridge supplies absolute liquid endpoint enthalpy");
  ok &= check(environment.environment.far_gas_density_kg_per_m3 > 0 &&
                  environment.environment.far_gas_density_kg_per_m3 !=
                      environment.environment.gas.film_density_kg_per_m3,
              "TAB consumes far gas density, not one-third film density");
  auto changed = parcel;
  changed.temperature_k += 10;
  auto warmer = bridge.query(changed, 0.6, ParcelPass::predictor, revision);
  ok &=
      check(warmer.environment.available &&
                std::abs(warmer.environment.liquid_absolute_enthalpy_j_per_kg +
                         90000) < 1e-8,
            "changed surface temperature is freshly queried without stale film "
            "cache");
  auto wrong = revision;
  ++wrong.input_revision;
  ok &= check(bridge.query(parcel, 0, ParcelPass::predictor, wrong).status ==
                  portable::Status::stale_revision,
              "bridge bound revision refuses another candidate");
  sampler.fail = true;
  ok &= check(
      !bridge.sample(parcel, 0, ParcelPass::predictor, revision).available,
      "sampler failure does not return old environment");
  sampler.fail = false;
  FixedAsParcelIntervalProvider interval(bridge);
  const auto integrated =
      interval.advance(parcel, 0, 1e-7, ParcelPass::corrector, revision);
  ok &= check(integrated.available && integrated.exchange.available &&
                  integrated.parcel.droplet_mass_kg < parcel.droplet_mass_kg,
              "P4 full-gas film and asset bridge drives real A-S interval");
  ok &= check(integrated.available &&
                  std::abs(integrated.exchange.parcel_liquid_mass_delta_kg +
                           integrated.exchange.gas_mass_delta_kg) < 1e-25,
              "A-S physical transfer counts evaporated mass once");
  ok &= check(sampler.queries >= 5,
              "each interval state obtains a new external gas sample");
  const portable::Revision next{5, 10, 1};
  ok &= check(
      bridge.bind_revision(next) == portable::Status::success &&
          interval
              .advance(integrated.parcel, 0, 1e-7, ParcelPass::corrector, next)
              .available &&
          bridge.query(parcel, 0, ParcelPass::predictor, revision).status ==
              portable::Status::stale_revision,
      "next accepted step reuses the actual film and A-S lane with a new "
      "revision");
  auto invalid_revision = next;
  invalid_revision.algorithm_version = 2;
  ok &= check(bridge.bind_revision(invalid_revision) ==
                      portable::Status::invalid_input &&
                  bridge.query(parcel, 0, ParcelPass::predictor, next)
                      .environment.available,
              "failed revision binding preserves the current prepared lane");
  return ok ? 0 : 1;
}
