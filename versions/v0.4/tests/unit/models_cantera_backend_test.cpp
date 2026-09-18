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
  for (std::size_t n : {2U, 4U})
    for (auto intervals : {esf::detail::ReactionIntervals::full,
                           esf::detail::ReactionIntervals::two_halves}) {
      const auto ensemble = persistent.react(
          {{query.revision, backend->closure_identity().fingerprint, n, 2,
            fields},
           query.revision,
           &backend->closure_identity(),
           &backend->gas_identity(),
           pressures,
           densities,
           0,
           1e-6,
           intervals},
          *backend);
      HUNDUN_CHECK(ensemble.status == portable::Status::success);
      HUNDUN_CHECK(
          ensemble.chemistry_call_count ==
          n * (intervals == esf::detail::ReactionIntervals::full ? 1 : 2));
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

void test_interval_call_order(bool reference=false) {
  using namespace hundun::v04;
  auto c = config();
  if(reference)c.chemistry={0.,1e-10,5000,true};
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(c);
  chemistry::CanteraWorkspacePool pool(runtime, 1);
  auto gas = chemistry::make_cantera_backend(c, pool);
  double y[]{1., 0.}, other_y[]{.3, .7}, d[2]{}, h[2]{}, w[2]{};
  portable::GasQuery q{{1, 1, 1},
                       gas->composition().fingerprint,
                       portable::GasStateCoordinates::pressure_temperature,
                       101325,
                       0,
                       1200,
                       y,
                       2};
  portable::GasQueryOutput sample{{}, d, h, w, 2};
  HUNDUN_CHECK(gas->query_gas(q, sample) == portable::Status::success);
  q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  q.enthalpy_j_per_kg = sample.sample.enthalpy_j_per_kg;
  double final[2]{}, delta[2]{}, expected[2]{}, temperature{};
  portable::GasAdvanceOutput output{{}, final, delta, 2};
  for (unsigned repeat = 0; repeat < 32; ++repeat) {
    HUNDUN_CHECK(gas->advance_gas({q, 0, 1e-4}, output) ==
                 portable::Status::success);
    if (repeat == 0) {
      expected[0] = final[0];
      expected[1] = final[1];
      temperature = output.final_sample.temperature_k;
    }
    HUNDUN_CHECK(final[0] == expected[0] && final[1] == expected[1]);
    HUNDUN_CHECK(output.final_sample.temperature_k == temperature);
    auto other = q;
    other.coordinates = portable::GasStateCoordinates::pressure_temperature;
    other.temperature_k = 1800;
    other.pressure_pa = 202650;
    other.mass_fractions = other_y;
    HUNDUN_CHECK(gas->query_gas(other, sample) == portable::Status::success);
    other.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
    other.enthalpy_j_per_kg = sample.sample.enthalpy_j_per_kg;
    HUNDUN_CHECK(gas->advance_gas({other, 0, 3e-4}, output) ==
                 portable::Status::success);
  }
}

void test_autonomous_interval_epoch() {
  using namespace hundun::v04;
  const auto c = config();
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(c);
  chemistry::CanteraWorkspacePool pool(runtime, 1);
  auto gas = chemistry::make_cantera_backend(c, pool);
  double y[2]{.8, .2}, d[2], h[2], rates[2], result[2], delta[2];
  portable::GasQuery q{{1, 1, 1},
                       gas->composition().fingerprint,
                       portable::GasStateCoordinates::pressure_temperature,
                       101325,
                       0,
                       1200,
                       y,
                       2};
  portable::GasQueryOutput sample{{}, d, h, rates, 2};
  HUNDUN_CHECK(gas->query_gas(q, sample) == portable::Status::success);
  q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  q.enthalpy_j_per_kg = sample.sample.enthalpy_j_per_kg;
  portable::GasAdvanceOutput out{{}, result, delta, 2};
  std::array<double, 2> expected{}, expected_delta{};
  double temperature = 0;
  unsigned steps = 0;
  for (double epoch : {0., .021806297823786736, 1., 1000000.}) {
    HUNDUN_CHECK(gas->advance_gas({q, epoch, 1e-6}, out) ==
                 portable::Status::success);
    if (epoch == 0) {
      std::copy_n(result, 2, expected.data());
      std::copy_n(delta, 2, expected_delta.data());
      temperature = out.final_sample.temperature_k;
      steps = out.internal_step_count;
    }
    HUNDUN_CHECK(std::equal(expected.begin(), expected.end(), result));
    HUNDUN_CHECK(
        std::equal(expected_delta.begin(), expected_delta.end(), delta));
    HUNDUN_CHECK(out.final_sample.temperature_k == temperature);
    HUNDUN_CHECK(out.internal_step_count == steps &&
                 out.completed_duration_s == 1e-6);
  }
}

void test_interval_budget_preserves_output() {
  using namespace hundun::v04;
  auto c = config();
  c.chemistry.maximum_internal_steps = 1;
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(c);
  chemistry::CanteraWorkspacePool pool(runtime, 1);
  auto gas = chemistry::make_cantera_backend(c, pool);
  double y[2]{.8, .2}, d[2], h[2], rates[2];
  portable::GasQuery q{{1, 1, 1},
                       gas->composition().fingerprint,
                       portable::GasStateCoordinates::pressure_temperature,
                       101325,
                       0,
                       1200,
                       y,
                       2};
  portable::GasQueryOutput sample{{}, d, h, rates, 2};
  HUNDUN_CHECK(gas->query_gas(q, sample) == portable::Status::success);
  q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  q.enthalpy_j_per_kg = sample.sample.enthalpy_j_per_kg;
  double result[2]{17., 19.}, delta[2]{23., 29.};
  portable::GasAdvanceOutput out{{}, result, delta, 2};
  HUNDUN_CHECK(gas->advance_gas({q, .021806297823786736, 1e-4}, out) ==
               portable::Status::provider_failure);
  HUNDUN_CHECK(result[0] == 17. && result[1] == 19. && delta[0] == 23. &&
               delta[1] == 29.);
  HUNDUN_CHECK(out.completed_duration_s == 0 && out.internal_step_count == 0);
  // A rejected interval leaves the shared workspace ready for the next query.
  HUNDUN_CHECK(gas->advance_gas({q, 1., 0.}, out) == portable::Status::success);
  HUNDUN_CHECK(result[0] == y[0] && result[1] == y[1]);
  HUNDUN_CHECK(delta[0] == 0. && delta[1] == 0.);
}

void test_continuous_nasa_and_identity(const std::filesystem::path &path) {
  using namespace hundun::v04;
  chemistry::CanteraBackendConfig source;
  source.mechanism = {
      path, "6e4f08f1a7280178f6925ad70df92f575d1007e6b3a5a189d49ef865bdfab29c",
      "nitrogen"};
  source.species_names = {"N2"};
  source.chemistry = {1e-10, 1e-18, 5000};
  auto continuous = source;
  continuous.continuous_enthalpy = true;
  auto raw_runtime = std::make_shared<chemistry::CanteraBackendRuntime>(source);
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(continuous);
  chemistry::CanteraWorkspacePool raw_pool(raw_runtime, 1), pool(runtime, 2);
  bool mismatched = false;
  try {
    static_cast<void>(chemistry::make_cantera_backend(source, pool));
  } catch (const std::invalid_argument &) {
    mismatched = true;
  }
  HUNDUN_CHECK(mismatched);
  auto raw = chemistry::make_cantera_backend(source, raw_pool);
  auto gas = chemistry::make_cantera_backend(continuous, pool);
  auto second = chemistry::make_cantera_backend(continuous, pool);
  HUNDUN_CHECK(
      !portable::same_gas_identity(raw->gas_identity(), gas->gas_identity()));
  HUNDUN_CHECK(
      portable::same_gas_identity(gas->gas_identity(), second->gas_identity()));
  HUNDUN_CHECK(raw->gas_identity().enthalpy_reference ==
               gas->gas_identity().enthalpy_reference);
  HUNDUN_CHECK(raw->closure_identity().fingerprint ==
               gas->closure_identity().fingerprint);
  double y = 1, d{}, h{}, rate{}, raw_d{}, raw_h{}, raw_rate{};
  portable::GasQuery query{{1, 1, 1},
                           gas->composition().fingerprint,
                           portable::GasStateCoordinates::pressure_temperature,
                           101325,
                           0,
                           300,
                           &y,
                           1};
  portable::GasQueryOutput output{{}, &d, &h, &rate, 1};
  portable::GasQueryOutput original{{}, &raw_d, &raw_h, &raw_rate, 1};
  for (double temperature :
       {300., 999.999999, 1000., 1000.000001, 1200., 2500.}) {
    query.coordinates = portable::GasStateCoordinates::pressure_temperature;
    query.temperature_k = temperature;
    HUNDUN_CHECK(raw->query_gas(query, original) == portable::Status::success);
    HUNDUN_CHECK(gas->query_gas(query, output) == portable::Status::success);
    HUNDUN_CHECK_NEAR(output.sample.cp_j_per_kg_k,
                      original.sample.cp_j_per_kg_k, 1e-10);
    // Direct integration of the published low/high NASA7 polynomials at
    // 1000 K gives the high-interval shift below; cp and low-interval h stay
    // fixed.
    HUNDUN_CHECK_NEAR(output.sample.enthalpy_j_per_kg -
                          original.sample.enthalpy_j_per_kg,
                      temperature > 1000. ? 0.18555729389759495 : 0., 2e-9);
    query.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
    query.enthalpy_j_per_kg = output.sample.enthalpy_j_per_kg;
    const auto ph_status = second->query_gas(query, output);
    HUNDUN_CHECK(ph_status == portable::Status::success);
    HUNDUN_CHECK_NEAR(output.sample.temperature_k, temperature, 1e-7);
  }
}

void test_case_temperature_interval(const std::filesystem::path &path) {
  using namespace hundun::v04;
  chemistry::CanteraBackendConfig source;
  source.mechanism = {
      path, "6e4f08f1a7280178f6925ad70df92f575d1007e6b3a5a189d49ef865bdfab29c",
      "nitrogen"};
  source.species_names = {"N2"};
  source.chemistry = {1e-10, 1e-18, 5000};
  source.continuous_enthalpy = true;
  auto extended = source;
  extended.minimum_temperature = 273.15;
  extended.maximum_temperature = 3500;
  auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(extended);
  chemistry::CanteraWorkspacePool pool(runtime, 1);
  auto gas = chemistry::make_cantera_backend(extended, pool);
  double y = 1, d = 0, h = 0, w = 0;
  portable::GasQuery q{{1, 1, 1},
                       gas->composition().fingerprint,
                       portable::GasStateCoordinates::pressure_temperature,
                       100000,
                       0,
                       295,
                       &y,
                       1};
  portable::GasQueryOutput out{{}, &d, &h, &w, 1};
  HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::success);
  // Original N2 low-temperature NASA polynomial, evaluated at the actual inlet
  // temperature.
  const double t = 295, r = 8314.46261815324 / 28.014;
  const double cp =
      r * (3.298677 +
           t * (.0014082404 +
                t * (-3.963222e-6 + t * (5.641515e-9 + t * (-2.444854e-12)))));
  HUNDUN_CHECK_NEAR(out.sample.cp_j_per_kg_k, cp, 1e-9);
  q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  q.enthalpy_j_per_kg = out.sample.enthalpy_j_per_kg;
  HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::success);
  HUNDUN_CHECK_NEAR(out.sample.temperature_k, 295, 1e-8);
  for (double target : {0., 5., -5., 50., -50.}) {
    q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
    q.enthalpy_j_per_kg = target;
    HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::success);
    HUNDUN_CHECK_NEAR(out.sample.enthalpy_j_per_kg, target,
                      1e-10 * std::max(1., std::abs(target)));
  }
  for (double temperature : {273.15, 3500.}) {
    q.coordinates = portable::GasStateCoordinates::pressure_temperature;
    q.temperature_k = temperature;
    HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::success);
  }
  for (double temperature : {273.14, 3500.01}) {
    q.temperature_k = temperature;
    HUNDUN_CHECK(gas->query_gas(q, out) == portable::Status::invalid_input);
  }
  auto raw_runtime = std::make_shared<chemistry::CanteraBackendRuntime>(source);
  chemistry::CanteraWorkspacePool raw_pool(raw_runtime, 1);
  auto raw = chemistry::make_cantera_backend(source, raw_pool);
  q.temperature_k = 295;
  HUNDUN_CHECK(raw->query_gas(q, out) == portable::Status::invalid_input);
  chemistry::CanteraWorkspacePool mismatch_pool(runtime, 1);
  bool mismatch = false;
  try {
    static_cast<void>(chemistry::make_cantera_backend(source, mismatch_pool));
  } catch (const std::invalid_argument &) {
    mismatch = true;
  }
  HUNDUN_CHECK(mismatch);
  HUNDUN_CHECK(
      !portable::same_gas_identity(raw->gas_identity(), gas->gas_identity()));
  HUNDUN_CHECK(raw->closure_identity().fingerprint ==
               gas->closure_identity().fingerprint);
  for (auto bounds :
       {std::array<double, 2>{273, 0}, std::array<double, 2>{3500, 273},
        std::array<double, 2>{-1, 3500}}) {
    auto bad = extended;
    bad.minimum_temperature = bounds[0];
    bad.maximum_temperature = bounds[1];
    bool rejected = false;
    try {
      chemistry::CanteraBackendRuntime ignored(bad);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
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
    HUNDUN_CHECK(argc == 3);
    mechanism_path = argv[1];
    test_interval_call_order();
    test_interval_call_order(true);
    test_autonomous_interval_epoch();
    test_interval_budget_preserves_output();
    test_continuous_nasa_and_identity(argv[2]);
    test_case_temperature_interval(argv[2]);
    test_exact_identity_and_independent_lanes();
    test_identity_mismatch_and_runtime_lifetime();
    test_controls_are_validated_without_product_schema();
    test_neutral_gas_query_and_closure_bridge();
    test_expired_pool_returns_explicit_failure();
  });
}
