// SPDX-License-Identifier: Apache-2.0
#include "models_chemistry_adapter_detail.hpp"
#include "models_spray_properties_detail.hpp"
#include <cmath>
#include <iostream>
#include <limits>
using namespace hundun::v04;
class FailureGas final : public portable::GasQueryProvider {
public:
  chemistry::detail::AnalyticIsomerBackend &gas;
  int calls{}, fail_at{3};
  bool stale{};
  explicit FailureGas(chemistry::detail::AnalyticIsomerBackend &g) : gas(g) {}
  const portable::GasIdentity &gas_identity() const noexcept override {
    return gas.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &o) noexcept override {
    if (++calls == fail_at) {
      o.sample = {};
      return portable::Status::provider_failure;
    }
    auto s = gas.query_gas(q, o);
    if (stale)
      ++o.sample.revision.input_revision;
    return s;
  }
};
bool check(bool ok, const char *name) {
  if (!ok)
    std::cerr << "FAIL: " << name << '\n';
  return ok;
}
int main(int argc, char **argv) {
  spray::detail::LiquidAsset asset;
  asset.pack.material_fingerprint = 1;
  asset.pack.minimum_temperature_k = 280;
  asset.pack.maximum_temperature_k = 400;
  asset.pack.cp_j_per_kg_k = {
      spray::detail::TemperatureCorrelationKind::polynomial_cubic,
      300,
      {500, 2, 0, 0}};
  asset.reference_temperature_k = 300;
  asset.reference_liquid_enthalpy_j_per_kg = -300000;
  auto h = spray::detail::evaluate_liquid_enthalpy(asset, 310);
  // Integral_300^310 [500+2(T-300)]dT=5100 J/kg.
  bool ok =
      check(h.available && std::abs(h.liquid_enthalpy_j_per_kg + 294900) < 1e-9,
            "absolute liquid enthalpy includes analytic cp integral");
  if (!check(argc == 3, "two synthetic asset paths required"))
    return 1;
  chemistry::detail::AnalyticIsomerBackend gas;
  auto loaded =
      spray::detail::load_liquid_asset(argv[1], 1, gas.gas_identity());
  ok &= check(!loaded.available &&
                  loaded.status == portable::Status::identity_mismatch,
              "wrong raw-content identity rejects file before publication");
  loaded = spray::detail::load_liquid_asset(
      argv[1], UINT64_C(6004043157121730787), gas.gas_identity());
  ok &=
      check(loaded.available,
            "constant synthetic asset loads by exact content and gas identity");
  double y[]{0.01, 0.99};
  spray::detail::FilmQueryInput input;
  input.expected_revision = {4, 9, 1};
  input.surface_temperature_k = 298.15;
  input.far_gas = {input.expected_revision,
                   gas.composition().fingerprint,
                   portable::GasStateCoordinates::pressure_enthalpy,
                   100000,
                   101000,
                   0,
                   y,
                   2};
  spray::detail::FilmQueryWorkspace film(2);
  auto sampled = film.query(loaded.asset, gas, input);
  ok &=
      check(sampled.available && std::abs(sampled.film.reference_temperature_k -
                                          331.48333333333335) < 1e-10,
            "one-third temperature uses complete far PH gas query");
  ok &=
      check(sampled.available && std::abs(film.film_mass_fractions()[0] -
                                          0.043333333333333334) < 1e-14,
            "surface equilibrium and carrier composition produce full film Y");
  ok &= check(std::abs(sampled.liquid_enthalpy_j_per_kg + 100000) < 1e-10 &&
                  std::abs(sampled.vapor_enthalpy_j_per_kg - 100000) < 1e-10,
              "absolute liquid and vapor references close latent heat");
  auto changed = input;
  changed.far_gas.revision.input_revision++;
  sampled = film.query(loaded.asset, gas, changed);
  ok &= check(!sampled.available &&
                  sampled.status == portable::Status::stale_revision &&
                  !film.film_mass_fractions(),
              "stale far sample produces no film candidate");
  auto broken = loaded.asset;
  broken.reference_liquid_enthalpy_j_per_kg += 1000;
  sampled = film.query(broken, gas, input);
  ok &= check(!sampled.available &&
                  sampled.status == portable::Status::conservation_failure,
              "inconsistent liquid vapor enthalpy reference is rejected");
  FailureGas failed(gas);
  sampled = film.query(loaded.asset, failed, input);
  ok &= check(!sampled.available &&
                  sampled.status == portable::Status::provider_failure &&
                  !film.film_mass_fractions(),
              "film provider failure discards all three query results");
  failed.calls = 0;
  failed.fail_at = 0;
  failed.stale = true;
  sampled = film.query(loaded.asset, failed, input);
  ok &= check(!sampled.available &&
                  sampled.status == portable::Status::stale_revision,
              "provider cannot return old revision");
  changed = input;
  changed.surface_temperature_k = 401;
  ok &= check(!film.query(loaded.asset, gas, changed).available,
              "liquid range boundary is enforced");
  changed = input;
  changed.far_gas.pressure_pa = 5000;
  ok &= check(film.query(loaded.asset, gas, changed).status ==
                  portable::Status::unavailable,
              "boiling conditions explicitly unsupported");
  chemistry::detail::AnalyticIsomerBackend reordered(3, true, 1200);
  auto beta = spray::detail::load_liquid_asset(
      argv[2], UINT64_C(668675689539421851), reordered.gas_identity());
  ok &= check(beta.available && beta.asset.vapor_species_index == 0,
              "second synthetic package binds reordered B vapor");
  auto beta_h = spray::detail::evaluate_liquid_enthalpy(beta.asset, 308.15);
  ok &= check(beta_h.available &&
                  std::abs(beta_h.liquid_enthalpy_j_per_kg + 294900) < 1e-8,
              "second package has independent nonconstant cp integral");
  changed = input;
  changed.far_gas.composition_fingerprint = reordered.composition().fingerprint;
  changed.far_gas.enthalpy_j_per_kg = 219000;
  sampled = film.query(beta.asset, reordered, changed);
  ok &=
      check(sampled.available && std::abs(film.film_mass_fractions()[0] -
                                          0.056666666666666664) < 1e-13,
            "Clausius reference and reordered vapor map produce distinct film");
  ok &= check(!spray::detail::load_liquid_asset(
                   argv[2], UINT64_C(668675689539421851), gas.gas_identity())
                   .available,
              "asset refuses unapproved gas species ordering");
  return ok ? 0 : 1;
}
