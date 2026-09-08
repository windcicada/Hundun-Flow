// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

// Small, prescribed-state model cases. No product driver or coupled flow.
// One invocation selects exactly one manifest configuration; no --all mode.
#include "models_chemistry_adapter_detail.hpp"
#include "models_esf_detail.hpp"
#include "models_spray_breakup_detail.hpp"
#include "models_spray_mechanics_detail.hpp"
#include "models_spray_transfer_detail.hpp"
#include "models_tcr_detail.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string_view>
#include <vector>

namespace {
using namespace hundun::v04;
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;
constexpr double pi = 3.141592653589793238462643383279502884;
bool near(double actual, double reference, double absolute, double relative) {
  return std::isfinite(actual) && std::isfinite(reference) &&
         std::abs(actual - reference) <=
             absolute + relative * std::abs(reference);
}
bool check(bool passed, const char *observation) {
  if (!passed)
    std::cerr << "FAIL " << observation << '\n';
  return passed;
}
void emit(const char *key, double value) {
  std::cout << key << '=' << std::setprecision(17) << value << '\n';
}
int count(std::string_view text) {
  if (text == "1")
    return 1;
  if (text == "2")
    return 2;
  if (text == "4")
    return 4;
  if (text == "10")
    return 10;
  if (text == "20")
    return 20;
  if (text == "40")
    return 40;
  return 0;
}

bool chemistry_case(bool beta, int steps) {
  chemistry::detail::AnalyticIsomerBackend backend(beta ? 3 : 2, beta,
                                                   beta ? 1200 : 1000);
  std::array<double, 2> y =
      beta ? std::array<double, 2>{0, 1} : std::array<double, 2>{1, 0};
  const auto &gas = backend.gas_identity();
  const auto &closure = backend.closure_identity();
  double diffusion[2]{}, species_h[2]{}, rates[2]{}, next_y[2]{}, delta[2]{};
  portable::GasQuery query{{0, 1, 1},
                           gas.composition_fingerprint,
                           portable::GasStateCoordinates::pressure_temperature,
                           101325,
                           0,
                           300,
                           y.data(),
                           2};
  portable::GasQueryOutput sampled{{}, diffusion, species_h, rates, 2};
  if (!check(backend.query_gas(query, sampled) == portable::Status::success,
             "initial PT query"))
    return false;
  const double initial_h = beta ? 102220 : 101850;
  bool passed =
      check(near(sampled.sample.enthalpy_j_per_kg, initial_h, 1e-8, 1e-13),
            "absolute initial h");
  query.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
  query.enthalpy_j_per_kg = initial_h;
  const double dt = .5 / steps;
  portable::GasAdvanceOutput out{{}, next_y, delta, 2};
  for (int i = 0; i < steps; ++i) {
    if (!check(backend.advance_gas({query, i * dt, dt}, out) ==
                   portable::Status::success,
               "0D interval"))
      return false;
    passed &= check(out.completed_duration_s == dt, "complete interval");
    passed &= check(near(delta[0] + delta[1], 0, 1e-14, 1e-13),
                    "mass and single-element balance");
    const double heat =
        -closure.species[0].formation_enthalpy_j_per_kg * delta[0] -
        closure.species[1].formation_enthalpy_j_per_kg * delta[1];
    passed &=
        check(near(out.integrated_heat_release_j_per_m3, heat, 1e-9, 1e-12),
              "formation heat report");
    passed &=
        check(near(out.final_sample.enthalpy_j_per_kg, initial_h, 1e-8, 1e-13),
              "closed total enthalpy");
    y = {next_y[0], next_y[1]};
  }
  const double expected_y = beta ? .223130160148429828933280470764013
                                 : .367879441171442321595523770161461;
  const double expected_t = beta ? 364.739153320964180922226627436332
                                 : 363.212055882855767840447622983854;
  passed &= check(near(y[beta ? 1 : 0], expected_y, 1e-14, 1e-12),
                  "independent exp(-kt) reference");
  passed &= check(near(out.final_sample.temperature_k, expected_t, 1e-9, 1e-12),
                  "independent constant-cp temperature");
  emit("reactant_mass_fraction", y[beta ? 1 : 0]);
  emit("temperature_k", out.final_sample.temperature_k);
  return passed;
}

bool iem_case(int fields, int steps) {
  combustion::ChemistryIdentity id{1, {{1, {1}, 0}, {1, {1}, 0}}, 42};
  std::array<double, 12> values{.2, .8, 100, .8, .2, 300,
                                .2, .8, 100, .8, .2, 300};
  esf::detail::Workspace workspace(2);
  esf::detail::Request q;
  q.identity = &id;
  q.dt_s = 2.0 / steps;
  q.mixing_time_s = 1;
  esf::detail::Report report;
  for (int i = 0; i < steps; ++i) {
    q.accepted = {{std::uint64_t(i), std::uint64_t(i + 1), 1},
                  42,
                  std::size_t(fields),
                  2,
                  values.data()};
    q.expected_revision = q.accepted.revision;
    q.random.accepted_step = i;
    report = workspace.advance(q);
    if (!check(report.status == portable::Status::success,
               "uniform IEM interval"))
      return false;
    std::copy_n(report.candidate.values, fields * 3, values.data());
  }
  emit("species_variance", report.variances[0]);
  emit("mean_enthalpy_j_per_kg", report.means[2]);
  return check(near(report.variances[0], .012180175491295142270459954547524,
                    1e-14, 1e-12),
               "variance .09 exp(-2)") &&
         check(near(report.means[0], .5, 1e-14, 1e-13) &&
                   near(report.means[2], 200, 1e-11, 1e-13),
               "IEM mean conservation");
}

bool tcr_case() {
  using namespace tcr::detail;
  History h;
  h.revision = {0, 1, 1};
  constexpr double ratio[]{.5, .75, 1};
  constexpr double reference[]{.139620389971936778000184409261214,
                               .225708114822568234916397374393457, 1.0 / 3};
  for (unsigned i = 0; i < 3; ++i) {
    TrialRequest q;
    q.expected_revision = h.revision;
    q.mode = Mode::experimental;
    q.mapping = {tcr::detail::Status::success, {.25, ratio[i]}, 71};
    q.initialization_sign = i == 0 ? -1 : 0;
    const auto trial = prepare(h, q);
    if (!check(trial.available &&
                   trial.status == tcr::detail::Status::success &&
                   trial.candidate.branch_sign == -1,
               "continuous accepted minus branch"))
      return false;
    if (!check(near(trial.mixer_control, reference[i], 1e-14, 1e-12),
               "independent TCR root"))
      return false;
    const auto committed = accept(h, trial, {i + 1, i + 2, 1});
    if (!check(committed.available, "TCR accepted value candidate"))
      return false;
    h = committed.candidate;
    emit("tcr_control", h.control);
  }
  return check(restore(h, h.revision).available,
               "accepted TCR history restore");
}

SprayParcelState droplet(double rho = 750) {
  SprayParcelState p;
  p.id = {17, 29};
  p.droplet_diameter_m = 1e-4;
  p.droplet_mass_kg = rho * pi * 1e-12 / 6;
  p.multiplicity = 4;
  p.temperature_k = 350;
  p.liquid_material_fingerprint = 91;
  return p;
}
TemperatureCorrelation constant(double x) {
  TemperatureCorrelation c;
  c.c[0] = x;
  return c;
}
LiquidPropertyPack heating_pack() {
  LiquidPropertyPack p;
  p.material_fingerprint = 91;
  p.minimum_temperature_k = 250;
  p.maximum_temperature_k = 1200;
  p.density_kg_per_m3 = constant(750);
  p.cp_j_per_kg_k = constant(2200);
  p.latent_heat_j_per_kg = constant(2.5e5);
  p.surface_tension_n_per_m = constant(.025);
  p.viscosity_pa_s = constant(8e-4);
  p.saturation_pressure.antoine_a = 4;
  p.saturation_pressure.pressure_scale_pa = 1;
  return p;
}
bool heating_case(int steps) {
  const auto pack = heating_pack();
  const LiquidPropertyService liquid(&pack, 1);
  FixedExchangeInput input;
  input.committed_parcel = droplet();
  input.liquid_properties = &liquid;
  input.duration_s = .001;
  input.internal_substeps = steps;
  input.maximum_internal_substeps = 40;
  input.predictor_environment = {{0, 0, 0}, 950,  .1,   100000, .8, 2e-5,
                                 .05,       2e-5, 1100, 28,     28, 1.02e6};
  input.corrector_environment = input.predictor_environment;
  const auto r = integrate_fixed_exchange(input);
  if (!check(r.succeeded(), "saturated single-drop fixed predictor/corrector"))
    return false;
  const double absolute = steps == 10 ? 6e-5 : (steps == 20 ? 1.5e-5 : 4e-6);
  emit("temperature_k", r.candidate_parcel.temperature_k);
  emit("absolute_reference_error_k",
       std::abs(r.candidate_parcel.temperature_k -
                371.426252621165355296643781549286));
  return check(
             near(r.candidate_parcel.temperature_k,
                  371.426252621165355296643781549286, absolute, 1e-12),
             "analytic exponential heating and registered refinement bound") &&
         check(r.exchange.gas_mass_delta_kg == 0,
               "saturated zero evaporation") &&
         check(near(r.exchange.parcel_thermochemical_enthalpy_delta_j +
                        r.exchange.thermal_exchange_to_gas_j,
                    0, 1e-14, 1e-12),
               "closed heating budget");
}
bool evaporation_case(int steps) {
  AbramzonSirignanoInput q;
  q.parcel = droplet();
  q.duration_s = .01 / steps;
  q.liquid_properties = {750, 2200, 2.5e5, 20000, .025, 8e-4};
  q.vapor_absolute_thermochemical_enthalpy_j_per_kg = 1.02e6;
  double gas_gain = 0;
  const double initial_mass = q.parcel.droplet_mass_kg * q.parcel.multiplicity;
  for (int i = 0; i < steps; ++i) {
    OneThirdFilmInput film{{0, 0, 0}, q.parcel.temperature_k,
                           350,       .2,
                           .02,       100000,
                           .8,        2e-5,
                           .05,       2e-5,
                           1100};
    q.film = sample_one_third_film(film);
    const auto r = evaluate_abramzon_sirignano(q);
    if (!check(r.succeeded() && r.bt_converged && !r.complete_evaporation,
               "Re0 A-S interval"))
      return false;
    if (!check(near(r.diameter_squared_rate_m2_per_s,
                    -3.463523737543514578914129737142e-8, 1e-20, 1e-12),
               "independent constant-film d2 slope"))
      return false;
    gas_gain += r.exchange.gas_mass_delta_kg;
    q.parcel = r.candidate_parcel;
  }
  const double d2 = q.parcel.droplet_diameter_m * q.parcel.droplet_diameter_m;
  emit("diameter_squared_m2", d2);
  emit("gas_mass_gain_kg", gas_gain);
  return check(near(d2, 9.653647626245648542108587026286e-9, 1e-20, 1e-12),
               "independent A-S d2 endpoint") &&
         check(near(q.parcel.droplet_mass_kg * q.parcel.multiplicity + gas_gain,
                    initial_mass, 1e-22, 1e-12),
               "parcel plus gas mass");
}
ParcelAccelerationSample zero_acceleration(const ParcelKinematicSample &,
                                           const void *) noexcept {
  return {true, {0, 0, 0}};
}
ParcelTrajectoryCandidate path(Vector3 position, Vector3 velocity,
                               double duration) {
  ParcelTrajectoryInput q;
  q.parcel = droplet();
  q.parcel.position_m = position;
  q.parcel.velocity_m_per_s = velocity;
  q.duration_s = duration;
  q.characteristic_length_m = .4;
  q.maximum_particle_cfl = .5;
  q.acceleration = zero_acceleration;
  return make_parcel_trajectory_candidate(q);
}
bool ballistic_case() {
  const auto r = path({1, 2, 3}, {2, -1, .5}, 2);
  if (!check(r.succeeded(), "ballistic interval"))
    return false;
  emit("x_m", r.parcel.position_m[0]);
  return check(near(r.parcel.position_m[0], 5, 1e-12, 1e-12) &&
                   near(r.parcel.position_m[1], 0, 1e-12, 1e-12) &&
                   near(r.parcel.position_m[2], 4, 1e-12, 1e-12),
               "x=x0+ut");
}
bool rebound_case() {
  CartesianMeshSpec spec;
  spec.kind = GeometryKind::uniform;
  spec.lower = {0, 0, 0};
  spec.upper = {4, 4, 4};
  spec.has_exact_cells = true;
  spec.exact_cells = {4, 4, 4};
  spec.minimum_spacing = {1e-9, 1e-9, 1e-9};
  spec.max_growth_ratio = 1;
  spec.limits = {64, 1048576};
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!check(bool(CartesianGeometryCompiler::compile(
                 MPI_COMM_SELF, spec, GeometryBudget{}, geometry, patch)),
             "4 cubed static geometry"))
    return false;
  const std::array<Real3, 8> v{{{1, 1, 1},
                                {1, 1, 2},
                                {1, 2, 1},
                                {1, 2, 2},
                                {2, 1, 1},
                                {2, 1, 2},
                                {2, 2, 1},
                                {2, 2, 2}}};
  constexpr unsigned indices[12][3]{{0, 1, 3}, {0, 3, 2}, {4, 6, 7}, {4, 7, 5},
                                    {0, 4, 5}, {0, 5, 1}, {2, 3, 7}, {2, 7, 6},
                                    {0, 2, 6}, {0, 6, 4}, {1, 5, 7}, {1, 7, 3}};
  std::array<TriangleInput, 12> triangles{};
  for (unsigned i = 0; i < 12; ++i)
    triangles[i] = {v[indices[i][0]], v[indices[i][1]], v[indices[i][2]]};
  StlScanPlan scan;
  ImmersedSurfacePlan surface;
  const StlScanBudget budget{1048576, 1048576, 1024, 1024, 1};
  if (!check(bool(StlScanCompiler::compile_triangles(
                 geometry, patch, {triangles.data(), triangles.size()},
                 CartesianAxis::y, budget, scan)) &&
                 bool(ImmersedSurfaceCompiler::compile(scan, surface)),
             "closed cube surface"))
    return false;
  const auto incoming = path({3, 1.5, 1.5}, {-2, .5, 0}, 1);
  const auto r = resolve_static_ibm_rebounds(
      incoming, surface, {ImmersedFluidSide::outside, .5, .25, 2});
  if (!check(r.succeeded() && r.event_count == 1, "one static IBM reflection"))
    return false;
  emit("collision_time_s", r.events[0].event_time_s);
  emit("final_x_m", r.parcel.position_m[0]);
  return check(near(r.events[0].event_time_s, .5, 1e-12, 1e-12) &&
                   near(r.parcel.position_m[0], 2.5, 1e-12, 1e-12) &&
                   near(r.parcel.position_m[1], 1.8125, 1e-12, 1e-12),
               "analytic normal/tangential restitution path");
}
bool tab_case(bool split) {
  TabBreakupInput q{.2, .3, 0, 1, 1, 0, 4, 1, 10, .4, {1, 1, 1, 1}};
  if (!split) {
    const auto r = evaluate_tab_breakup(q);
    if (!check(r.succeeded(), "unforced TAB oscillator"))
      return false;
    emit("deformation", r.candidate.deformation);
    return check(near(r.candidate.deformation, .246944755504361498428226187916,
                      1e-14, 1e-12) &&
                     near(r.candidate.deformation_rate_per_s,
                          -.077930423555659478374644849740, 1e-14, 1e-12),
                 "independent harmonic oscillator");
  }
  q = {1, 20000, 10, 1, 800, 8e-4, .025, 5e-5, 1, 1e-4, {2.0 / 3, .5, 8, 10}};
  TabRepresentativeSplitInput input;
  input.parent = droplet(800);
  input.parent.multiplicity = 100;
  input.parent.velocity_m_per_s = {20, -3, 4};
  input.tab_trigger = evaluate_tab_breakup(q);
  input.liquid_density_kg_per_m3 = 800;
  input.surface_tension_n_per_m = .025;
  input.liquid_absolute_thermochemical_enthalpy_j_per_kg = -2e5;
  input.child_parcel_count = 4;
  const auto r = generate_tab_representative_children(input);
  if (!check(r.succeeded(), "TAB representative children"))
    return false;
  emit("representative_diameter_m", r.representative_diameter_m);
  emit("energy_residual_j", r.total_energy_residual_j);
  return check(near(r.representative_diameter_m,
                    3.947368421052631578947368421053e-5, 1e-18, 1e-12) &&
                   near(r.transverse_speed_m_per_s, .5, 1e-14, 1e-12),
               "independent TAB diameter and transverse speed") &&
         check(near(r.transverse_kinetic_energy_j,
                    5.235987755982988730771072305466e-9, 1e-21, 1e-12) &&
                   near(r.total_energy_residual_j, 0, 1e-18, 0),
               "explicit TAB energy budget");
}
} // namespace
int main(int argc, char **argv) {
  bool passed = false;
  const std::string_view mode = argc > 1 ? argv[1] : "";
  if (mode == "chemistry" && argc == 4 &&
      (std::string_view(argv[2]) == "alpha" ||
       std::string_view(argv[2]) == "beta") &&
      (count(argv[3]) == 1 || count(argv[3]) == 2 || count(argv[3]) == 4))
    passed =
        chemistry_case(std::string_view(argv[2]) == "beta", count(argv[3]));
  else if (mode == "iem" && argc == 4 &&
           (count(argv[2]) == 2 || count(argv[2]) == 4) &&
           (count(argv[3]) == 1 || count(argv[3]) == 2 || count(argv[3]) == 4))
    passed = iem_case(count(argv[2]), count(argv[3]));
  else if (mode == "tcr" && argc == 2)
    passed = tcr_case();
  else if ((mode == "heating" || mode == "evaporation") && argc == 3 &&
           (count(argv[2]) == 10 || count(argv[2]) == 20 ||
            count(argv[2]) == 40))
    passed = mode == "heating" ? heating_case(count(argv[2]))
                               : evaporation_case(count(argv[2]));
  else if (mode == "ballistic" && argc == 2)
    passed = ballistic_case();
  else if (mode == "rebound" && argc == 2) {
    MPI_Init(&argc, &argv);
    passed = rebound_case();
    MPI_Finalize();
  } else if ((mode == "tab_oscillator" || mode == "tab_split") && argc == 2)
    passed = tab_case(mode == "tab_split");
  else {
    std::cerr << "Select one configuration from portable_v1_manifest.json\n";
    return 2;
  }
  std::cout << (passed ? "V1_CASE_PASS" : "V1_CASE_FAIL") << '\n';
  return passed ? 0 : 1;
}
