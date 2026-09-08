// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_cantera.hpp"
#include "models_chemistry_adapter_detail.hpp"
#include "models_esf_detail.hpp"
#include "models_spray_properties_detail.hpp"

#include "../support/chemistry_test_support.hpp"

#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {

std::filesystem::path mechanism_path;

hundun::v04::chemistry::CanteraBackendConfig config() {
  hundun::v04::chemistry::CanteraBackendConfig value;
  value.mechanism.file = mechanism_path;
  value.mechanism.sha256 =
      "c518a07cada5f1bddcdb308f0a2f695d92cc6373e173ffd87e96312530b52aee";
  value.mechanism.phase = "synthetic-gas";
  value.species_names = {"A", "B"};
  value.chemistry.relative_tolerance = 1.0e-10;
  value.chemistry.absolute_tolerance = 1.0e-18;
  value.chemistry.maximum_internal_steps = 5000;
  return value;
}

void test_exact_identity_and_independent_lanes() {
  auto runtime =
      std::make_shared<hundun::v04::chemistry::CanteraBackendRuntime>(config());
  HUNDUN_CHECK(runtime->composition().species.size() == 2U);
  HUNDUN_CHECK(runtime->composition().species[0].name == "A");
  HUNDUN_CHECK(runtime->composition().species[1].name == "B");
  HUNDUN_CHECK(runtime->mechanism_sha256() == config().mechanism.sha256);

  hundun::v04::chemistry::CanteraWorkspacePool pool(runtime, 2U);
  HUNDUN_CHECK(pool.workspace_count() == 2U);
  HUNDUN_CHECK(pool.workspaces_are_distinct());
  auto lane0 = hundun::v04::chemistry::make_cantera_backend(config(), pool);
  auto lane1 = hundun::v04::chemistry::make_cantera_backend(config(), pool);
  HUNDUN_CHECK(lane0->lane_index() == 0U);
  HUNDUN_CHECK(lane1->lane_index() == 1U);
  HUNDUN_CHECK(lane0->composition().fingerprint ==
               lane1->composition().fingerprint);

  bool missing = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::make_cantera_backend(config(), pool));
  } catch (const std::runtime_error &) {
    missing = true;
  }
  HUNDUN_CHECK(missing);
}

void test_identity_mismatch_and_runtime_lifetime() {
  auto wrong_hash = config();
  wrong_hash.mechanism.sha256[0] = '0';
  bool rejected = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::CanteraBackendRuntime(wrong_hash));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto runtime =
      std::make_shared<hundun::v04::chemistry::CanteraBackendRuntime>(config());
  hundun::v04::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  auto changed = config();
  changed.species_names = {"B", "A"};
  rejected = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::make_cantera_backend(changed, pool));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  runtime.reset();
  rejected = false;
  try {
    static_cast<void>(
        hundun::v04::chemistry::make_cantera_backend(config(), pool));
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_controls_are_validated_without_product_schema() {
  auto runtime =
      std::make_shared<hundun::v04::chemistry::CanteraBackendRuntime>(config());
  hundun::v04::chemistry::CanteraWorkspacePool pool(runtime, 1U);
  for (int mutation = 0; mutation < 4; ++mutation) {
    auto invalid = config();
    if (mutation == 0)
      invalid.chemistry.relative_tolerance = 0.0;
    if (mutation == 1)
      invalid.chemistry.absolute_tolerance = -1.0;
    if (mutation == 2)
      invalid.chemistry.maximum_internal_steps = 0;
    if (mutation == 3)
      invalid.chemistry.relative_tolerance =
          std::numeric_limits<double>::quiet_NaN();
    bool rejected = false;
    try {
      static_cast<void>(
          hundun::v04::chemistry::make_cantera_backend(invalid, pool));
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  // Rejected controls must not consume a lane.
  auto valid = hundun::v04::chemistry::make_cantera_backend(config(), pool);
  HUNDUN_CHECK(valid->lane_index() == 0U);
}

void test_neutral_gas_query_and_closure_bridge() {
  using namespace hundun::v04;
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(config());
  chemistry::CanteraWorkspacePool pool(runtime, 1);
  auto backend = chemistry::make_cantera_backend(config(), pool);
  double y[]{1, 0}, d[2]{}, h[2]{}, rates[2]{};
  portable::GasQuery query{{1, 2, 1},
                           backend->composition().fingerprint,
                           portable::GasStateCoordinates::pressure_temperature,
                           101325,
                           0,
                           1000,
                           y,
                           2};
  portable::GasQueryOutput output{{}, d, h, rates, 2};
  const auto first_status = backend->query_gas(query, output);
  HUNDUN_CHECK(first_status == portable::Status::success);
  // NASA7 A: cp/R=2.5, h/RT=2.5. H atomic weight=1.008 kg/kmol.
  HUNDUN_CHECK_NEAR(output.sample.enthalpy_j_per_kg, 20621187.04899117, 1e-5);
  HUNDUN_CHECK_NEAR(rates[0] + rates[1], 0, 1e-13);
  HUNDUN_CHECK(rates[0] < 0 && d[0] >= 0 && output.sample.viscosity_pa_s > 0);
  double mixed_y[]{0.5, 0.5};
  auto mixed = query;
  mixed.mass_fractions = mixed_y;
  portable::GasQueryOutput mixed_output{{}, d, h, rates, 2};
  HUNDUN_CHECK(backend->query_gas(mixed, mixed_output) ==
               portable::Status::success);
  HUNDUN_CHECK(d[0] > 0 && d[1] > 0);
  query.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  query.enthalpy_j_per_kg = output.sample.enthalpy_j_per_kg;
  HUNDUN_CHECK(backend->query_gas(query, output) == portable::Status::success);
  HUNDUN_CHECK_NEAR(output.sample.temperature_k, 1000, 1e-8);
  double final_y[2], delta_y[2];
  portable::GasAdvanceOutput bounded{{}, final_y, delta_y, 2};
  HUNDUN_CHECK(backend->advance_gas({query, 0, 1e-6}, bounded) ==
               portable::Status::success);
  HUNDUN_CHECK(final_y[0] < 1 && final_y[1] > 0 &&
               bounded.internal_step_count > 0);
  HUNDUN_CHECK_NEAR(delta_y[0] + delta_y[1], 0, 1e-13);
  HUNDUN_CHECK_NEAR(bounded.final_sample.enthalpy_j_per_kg,
                    query.enthalpy_j_per_kg, 1e-5);
  esf::detail::Workspace persistent(2);
  double fields[]{1, 0, query.enthalpy_j_per_kg, 1, 0, query.enthalpy_j_per_kg,
                  1, 0, query.enthalpy_j_per_kg, 1, 0, query.enthalpy_j_per_kg};
  double pressures[]{101325, 101325, 101325, 101325}, densities[4]{};
  for (double &density : densities)
    density = output.sample.density_kg_per_m3;
  for (std::size_t n : {2U, 4U}) {
    const auto ensemble = persistent.react(
        {{query.revision, backend->closure_identity().fingerprint, n, 2,
          fields},
         query.revision,
         &backend->closure_identity(),
         &backend->gas_identity(),
         pressures,
         densities,
         0,
         1e-6},
        *backend);
    HUNDUN_CHECK(ensemble.status == portable::Status::success);
    HUNDUN_CHECK(ensemble.chemistry_call_count == 2 * n);
    if (ensemble.status == portable::Status::success) {
      HUNDUN_CHECK_NEAR(ensemble.candidate.values[0], final_y[0], 1e-8);
      HUNDUN_CHECK_NEAR(ensemble.candidate.values[0] +
                            ensemble.candidate.values[1],
                        1, 1e-13);
      HUNDUN_CHECK(ensemble.candidate.values[2] == query.enthalpy_j_per_kg);
    }
  }
  chemistry::detail::BackendAdapter adapter(
      *backend, *backend, backend->closure_identity(), query.revision);
  combustion::CombustionClosureRequest r;
  r.composition_fingerprint = adapter.identity().fingerprint;
  r.mean_state = {
      101325, output.sample.density_kg_per_m3, query.enthalpy_j_per_kg, {1, 0}};
  r.duration_s = 1e-6;
  combustion::CombustionClosureConfig cfg;
  const auto finite =
      combustion::evaluate_combustion_closure(r, cfg, adapter, &adapter);
  HUNDUN_CHECK(finite.succeeded());
  HUNDUN_CHECK(finite.chemistry_call_count == 2);
  HUNDUN_CHECK(
      finite.candidate.integrated_thermochemical_enthalpy_delta_j_per_m3 == 0);
  cfg.tci_closure = combustion::TciClosureKind::pasr_algebraic_v1;
  cfg.mixing_time = {0.001, 1e-5, 1e-5, 0.7, 0.5};
  const auto pasr =
      combustion::evaluate_combustion_closure(r, cfg, adapter, &adapter);
  HUNDUN_CHECK(pasr.succeeded());
  HUNDUN_CHECK(pasr.chemistry_call_count == 2);
  cfg.tci_closure = combustion::TciClosureKind::esf_tpdf;
  r.stochastic_fields.assign(4, r.mean_state);
  const auto ensemble =
      combustion::evaluate_combustion_closure(r, cfg, adapter, &adapter);
  HUNDUN_CHECK(ensemble.succeeded());
  HUNDUN_CHECK(ensemble.chemistry_call_count == 8);
  query.temperature_k = 100;
  query.coordinates = portable::GasStateCoordinates::pressure_temperature;
  HUNDUN_CHECK(backend->query_gas(query, output) ==
               portable::Status::invalid_input);
  HUNDUN_CHECK(output.sample.composition_fingerprint == 0);

  // Independent synthetic liquid consistent with NASA7 A: cp=2.5R/M,
  // h_A(298.15)=6148206.918656717 J/kg; choose L=1e6 J/kg.
  spray::detail::LiquidAsset liquid;
  liquid.gas_identity = backend->gas_identity();
  liquid.vapor_species_name = "A";
  liquid.vapor_species_index = 0;
  liquid.vapor_molecular_weight_kg_per_kmol = 1.008;
  liquid.reference_temperature_k = 298.15;
  liquid.reference_liquid_enthalpy_j_per_kg = 5148206.918656717;
  auto &pack = liquid.pack;
  pack.material_fingerprint = 321;
  pack.minimum_temperature_k = 280;
  pack.maximum_temperature_k = 400;
  pack.density_kg_per_m3.c[0] = 800;
  pack.cp_j_per_kg_k.c[0] = 20621.18704899117;
  pack.latent_heat_j_per_kg.c[0] = 1e6;
  pack.surface_tension_n_per_m.c[0] = 0.025;
  pack.viscosity_pa_s.c[0] = 0.001;
  pack.saturation_pressure.pressure_scale_pa = 6000;
  double far_y[]{0.01, 0.99};
  auto far_query = query;
  far_query.temperature_k = 1000;
  far_query.pressure_pa = 100000;
  far_query.mass_fractions = far_y;
  HUNDUN_CHECK(backend->query_gas(far_query, output) ==
               portable::Status::success);
  far_query.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  far_query.enthalpy_j_per_kg = output.sample.enthalpy_j_per_kg;
  spray::detail::FilmQueryWorkspace film(2);
  spray::detail::FilmQueryInput film_input{
      far_query, query.revision, {}, 298.15};
  const auto film_result = film.query(liquid, *backend, film_input);
  HUNDUN_CHECK(film_result.available);
  HUNDUN_CHECK_NEAR(film_result.film.reference_temperature_k, 532.1, 1e-8);
  HUNDUN_CHECK_NEAR(film.film_mass_fractions()[0], 0.043333333333333334, 1e-13);
  HUNDUN_CHECK_NEAR(film_result.latent_heat_consistency_residual_j_per_kg, 0,
                    1e-6);
  HUNDUN_CHECK(film_result.film.vapor_diffusivity_m2_per_s > 0);
}

void test_expired_pool_returns_explicit_failure() {
  using namespace hundun::v04;
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(config());
  std::unique_ptr<chemistry::CanteraBackend> backend;
  {
    chemistry::CanteraWorkspacePool pool(runtime, 1);
    backend = chemistry::make_cantera_backend(config(), pool);
  }
  portable::GasQuery query;
  portable::GasQueryOutput output;
  HUNDUN_CHECK(backend->query_gas(query, output) ==
               portable::Status::unavailable);
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    mechanism_path = argv[1];
    test_exact_identity_and_independent_lanes();
    test_identity_mismatch_and_runtime_lifetime();
    test_controls_are_validated_without_product_schema();
    test_neutral_gas_query_and_closure_bridge();
    test_expired_pool_returns_explicit_failure();
  });
}
