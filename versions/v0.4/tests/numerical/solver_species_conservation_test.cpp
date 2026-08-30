// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <= 1.0e-13 * std::max(1.0, std::abs(expected));
}

struct OwnedField {
  std::vector<double> bytes;
  FieldView view{};
};

struct OwnedFaceField {
  std::vector<double> bytes;
  FaceFieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision) {
  OwnedField field;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field.bytes.assign(nx * ny * nz * components, 0.0);
  field.view.base = field.bytes.data() + ghosts + ghosts * nx + ghosts * nx * ny;
  field.view.interior = cells;
  field.view.ghosts = {ghosts, ghosts, ghosts};
  field.view.components = components;
  field.view.stride_y = nx;
  field.view.stride_z = nx * ny;
  field.view.component_stride = nx * ny * nz;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = 2000U + id;
  field.view.revision_domain = 9201U;
  return field;
}

OwnedFaceField make_face_field(CartesianAxis axis, Int3 cells,
                               StorageIdentity identity) {
  OwnedFaceField field;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) {
    ++extents.x;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
  } else {
    ++extents.z;
  }
  field.bytes.assign(static_cast<std::size_t>(extents.x) * extents.y *
                         extents.z,
                     0.0);
  field.view = {field.bytes.data(),
                extents,
                static_cast<std::size_t>(extents.x),
                static_cast<std::size_t>(extents.x) * extents.y,
                axis,
                identity,
                9301U};
  return field;
}

void fill_field(OwnedField& field, double value) {
  std::fill(field.bytes.begin(), field.bytes.end(), value);
}

CartesianMeshSpec mesh_spec(std::int32_t n) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  mesh.minimum_spacing = {1.0 / n, 1.0 / n, 1.0 / n};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(n) * n * n;
  mesh.limits.max_memory_bytes_per_rank = 1U << 28U;
  return mesh;
}

SpeciesThermophysicalSpec thermophysical_species(std::string_view name,
                                                  double molecular_weight) {
  SpeciesThermophysicalSpec species;
  species.stable_name.assign(name.data(), name.size());
  species.molecular_weight = molecular_weight;
  species.temperature_switch = 1000.0;
  species.nasa7_low[0U] = 3.5;
  species.nasa7_high[0U] = 3.5;
  species.viscosity_reference = 1.8e-5;
  species.conductivity = 0.026;
  return species;
}

ThermophysicalSpec thermo_spec() {
  ThermophysicalSpec value;
  value.data_file = "analytic.d";
  value.minimum_temperature = 200.0;
  value.maximum_temperature = 2000.0;
  value.temperature_relative_tolerance = 1.0e-12;
  value.maximum_temperature_iterations = 64U;
  value.closed_mass_relative_tolerance = 1.0e-12;
  value.maximum_closed_mass_iterations = 32U;
  value.maximum_closed_mass_relative_step = 0.2;
  value.species.push_back(thermophysical_species("species_a", 28.0));
  value.species.push_back(thermophysical_species("species_b", 32.0));
  value.species.back().nasa7_low[0U] = 4.25;
  value.species.back().nasa7_high[0U] = 4.25;
  return value;
}

constexpr FieldId kDensity = 0U;
constexpr FieldId kVelocity = 1U;
constexpr FieldId kPressure = 2U;
constexpr FieldId kEnthalpy = 3U;
constexpr FieldId kTemperature = 4U;
constexpr FieldId kSpecies = 5U;
constexpr FieldId kPassive = 6U;
constexpr FieldId kEffectiveViscosity = 7U;
constexpr FieldId kCompressibility = 8U;
constexpr FieldId kVelocityGradient = 9U;

struct ProductionFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  EquationPlanSet equations;
};

bool make_production_fixture(std::int32_t n, ProductionFixture& out) {
  const CartesianMeshSpec mesh = mesh_spec(n);
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14a5ca1aU;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  model.transported_scalars = {
      {"species_a", TransportedScalarRole::species, 0.5, 2.0},
      {"tracer", TransportedScalarRole::passive_scalar, 2.0, 4.0}};
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;

  FieldRegistry registry;
  FieldId density = 0U;
  FieldId velocity = 0U;
  FieldId pressure = 0U;
  FieldId enthalpy = 0U;
  FieldId temperature = 0U;
  if (!registry.require_field("rho", 1U, 2U, density) ||
      density != kDensity ||
      !registry.require_field("U", 3U, 2U, velocity) ||
      velocity != kVelocity ||
      !registry.require_field("pi", 1U, 2U, pressure) ||
      pressure != kPressure ||
      !registry.require_field("h", 1U, 2U, enthalpy) ||
      enthalpy != kEnthalpy ||
      !registry.require_field("T", 1U, 2U, temperature) ||
      temperature != kTemperature ||
      !CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes,
                                 out.time)) {
    return false;
  }

  const ThermophysicalSpec thermophysics = thermo_spec();
  if (!ThermodynamicsPlan::compile(
          thermophysics,
          {model.transported_scalars.data(), model.transported_scalars.size()},
          out.thermodynamics) ||
      !TransportPlan::compile(thermophysics, out.thermodynamics,
                              out.transport)) {
    return false;
  }
  const std::array<FieldId, 10U> declared{
      kDensity,          kVelocity,       kPressure,
      kEnthalpy,         kTemperature,    kSpecies,
      kPassive,          kEffectiveViscosity,
      kCompressibility,  kVelocityGradient};
  if (!out.contributions.configure({declared.data(), declared.size()}) ||
      !out.contributions.freeze()) {
    return false;
  }

  const std::array<ScalarEquationSpec, 2U> scalar_specs{{
      {kSpecies, TransportedScalarRole::species, 0.5, 2.0},
      {kPassive, TransportedScalarRole::passive_scalar, 2.0, 4.0},
  }};
  EquationPlanSpec spec;
  spec.density = kDensity;
  spec.velocity = kVelocity;
  spec.pressure_perturbation = kPressure;
  spec.enthalpy = kEnthalpy;
  spec.temperature = kTemperature;
  spec.effective_viscosity = kEffectiveViscosity;
  spec.pressure_compressibility = kCompressibility;
  spec.velocity_gradient = kVelocityGradient;
  spec.pressure_reference = PressureReferenceKind::closed_mass;
  spec.scalars = {scalar_specs.data(), scalar_specs.size()};
  spec.closed_mass_service_stage = 1U;
  spec.maximum_cells_per_rank = static_cast<std::size_t>(n) * n * n;
  return static_cast<bool>(EquationPlanSet::compile(
      MPI_COMM_SELF, out.schemes, out.geometry, out.patch, out.boundary,
      out.contributions, out.thermodynamics, out.transport, spec,
      out.equations));
}

struct FinalFluxFixture {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxWriter writer;
};

bool make_linear_final_flux(const CartesianKernelPlan& kernels, Int3 cells,
                            FinalFluxFixture& fixture,
                            ConstFaceFluxView& committed) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("flux_dependency", 1U, 0U,
                              fixture.dependency) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{ArenaFieldRequest{
      fixture.dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  if (!ArenaLayout::compile(schema, {requests.data(), requests.size()},
                            layout) ||
      !StateLayers::allocate(layout, fixture.layers) ||
      !AttemptTransaction::create(fixture.layers.field_count(), 1U,
                                  fixture.layers.field_count(),
                                  fixture.transaction) ||
      !FaceFluxStorage::allocate_final(cells, fixture.storage)) {
    return false;
  }
  FinalFaceFluxAuthority authority;
  if (!authority.claim(41U, 0U, fixture.transaction, fixture.writer) ||
      !fixture.transaction.begin(fixture.layers) ||
      !fixture.transaction.revise_trial(fixture.dependency)) {
    return false;
  }
  const RevisionDependency dependency{
      AttemptTransaction::field_revision_source(fixture.dependency),
      fixture.transaction.trial_revision(fixture.dependency)};
  PendingFaceFluxView pending;
  if (!fixture.writer.begin_pending(fixture.transaction, fixture.storage,
                                    pending)) {
    return false;
  }

  OwnedField rho = make_field(80U, cells, 1U, 2U, 501U);
  OwnedField velocity = make_field(81U, cells, 3U, 2U, 502U);
  const double spacing = 1.0 / static_cast<double>(cells.x);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const Int3 cell{i, j, k};
        rho.view.unchecked(cell, 0U) = 1.0;
        velocity.view.unchecked(cell, 0U) =
            (static_cast<double>(i) + 0.5) * spacing;
        velocity.view.unchecked(cell, 1U) = 0.0;
        velocity.view.unchecked(cell, 2U) = 0.0;
      }
    }
  }
  const std::array<ConstFieldView, 2U> reads{as_const(rho.view),
                                             as_const(velocity.view)};
  const KernelInvocation call{{reads.data(), reads.size()}, {},
                              {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U, nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(kernels, call, pending)) &&
         static_cast<bool>(fixture.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(fixture.transaction.collective_finish(MPI_COMM_SELF,
                                                                  Status{})) &&
         static_cast<bool>(fixture.writer.committed(fixture.storage,
                                                     committed));
}

bool test_independent_species_closure() {
  const std::array<double, 3U> independent{0.15, 0.20, 0.25};
  const std::array<double, 3U> fluxes{1.5e-4, -4.0e-5, 7.0e-5};
  SpeciesClosure closure;
  bool passed = expect(static_cast<bool>(close_independent_species(
                           {independent.data(), independent.size()},
                           {fluxes.data(), fluxes.size()}, closure)),
                       "valid N-1 species close without clipping");
  passed &= expect(close(closure.dependent_mass_fraction, 0.40),
                   "dependent species is exactly one minus the independent sum");
  passed &= expect(close(closure.dependent_diffusive_flux, -1.8e-4),
                   "dependent diffusion flux is exactly minus the independent sum");

  std::array<double, 2U> invalid_negative{-1.0e-12, 0.2};
  std::array<double, 2U> invalid_sum{0.8, 0.3};
  const SpeciesClosure sentinel{9.0, 11.0};
  closure = sentinel;
  passed &= expect(close_independent_species(
                       {invalid_negative.data(), invalid_negative.size()}, {},
                       closure).code == StatusCode::numerical_failure &&
                       close(closure.dependent_mass_fraction,
                             sentinel.dependent_mass_fraction),
                   "negative independent species rejects atomically, without clipping");
  closure = sentinel;
  passed &= expect(close_independent_species(
                       {invalid_sum.data(), invalid_sum.size()}, {}, closure)
                           .code == StatusCode::numerical_failure &&
                       close(closure.dependent_mass_fraction,
                             sentinel.dependent_mass_fraction),
                   "independent sum above one rejects atomically, without renormalizing");
  return passed;
}

bool test_composition_dependent_production_eos() {
  ProductionFixture fixture;
  bool passed = expect(make_production_fixture(2, fixture),
                       "composition-dependent production thermo plan compiles");
  if (!passed) {
    return false;
  }

  constexpr double p_ref = 101325.0;
  constexpr double pi = 1750.0;
  constexpr double p_abs = p_ref + pi;
  constexpr double h = 650000.0;
  constexpr double mw_a = 28.0;
  constexpr double mw_b = 32.0;
  constexpr double cp_over_r_a = 3.5;
  constexpr double cp_over_r_b = 4.25;
  const std::array<double, 2U> independent_a_mass_fractions{0.15, 0.75};
  std::array<ThermoState, 2U> states{};
  std::array<double, 2U> oracle_psi{};

  for (std::size_t case_index = 0U;
       case_index < independent_a_mass_fractions.size(); ++case_index) {
    const double y_a = independent_a_mass_fractions[case_index];
    const double y_b = 1.0 - y_a;
    const std::array<double, 1U> independent{y_a};
    const double gas_constant =
        kUniversalGasConstant * (y_a / mw_a + y_b / mw_b);
    const double cp = kUniversalGasConstant *
                      (y_a * cp_over_r_a / mw_a +
                       y_b * cp_over_r_b / mw_b);
    const double temperature = h / cp;
    oracle_psi[case_index] = 1.0 / (gas_constant * temperature);
    const double density = p_abs * oracle_psi[case_index];

    passed &= expect(
        static_cast<bool>(fixture.thermodynamics.evaluate_from_reference_pressure(
            p_ref, pi, h, {independent.data(), independent.size()}, {},
            states[case_index])) &&
            close(states[case_index].rho, density) &&
            close(states[case_index].drho_dp_hY, oracle_psi[case_index]),
        "production EOS matches the independent ideal-gas mixture oracle");
  }

  passed &= expect(
      !close(states[0U].rho, states[1U].rho) &&
          !close(states[0U].drho_dp_hY, states[1U].drho_dp_hY),
      "different valid N-1 compositions change rho and drho/dp through the dependent species closure");
  const double mutated_fixed_r_density = p_abs * oracle_psi[0U];
  passed &= expect(
      !close(mutated_fixed_r_density, states[1U].rho) &&
          !close(oracle_psi[0U], states[1U].drho_dp_hY),
      "a composition-blind fixed-gas-constant mutation is killed numerically");
  return passed;
}

bool test_scalar_catalog_contract() {
  std::array<ScalarEquationSpec, 2U> scalars{{
      {7U, TransportedScalarRole::species, 0.7, 0.9},
      {8U, TransportedScalarRole::passive_scalar, 1.1, 0.8},
  }};
  bool passed = expect(scalars[0].field != scalars[1].field &&
                           scalars[0].role == TransportedScalarRole::species &&
                           scalars[1].role ==
                               TransportedScalarRole::passive_scalar,
                       "species and passive scalar share one typed catalog");
  passed &= expect(scalars[0].molecular_schmidt > 0.0 &&
                       scalars[0].turbulent_schmidt > 0.0,
                   "molecular and turbulent Schmidt numbers are explicit");
  return passed;
}

bool test_scalar_mass_diffusivity_oracle_and_atomicity() {
  const Int3 cells{4, 3, 2};
  OwnedField molecular = make_field(20U, cells, 1U, 0U, 51U);
  OwnedField turbulent = make_field(21U, cells, 1U, 0U, 52U);
  OwnedField diffusivity = make_field(22U, cells, 1U, 0U, 53U);
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        molecular.view.unchecked({i, j, k}, 0U) = 2.0 + i;
        turbulent.view.unchecked({i, j, k}, 0U) = 3.0 + j;
      }
    }
  }
  const ScalarEquationSpec spec{7U, TransportedScalarRole::species, 0.5, 2.0};
  bool passed = expect(static_cast<bool>(form_scalar_mass_diffusivity(
                           spec, as_const(molecular.view),
                           as_const(turbulent.view), {{0, 0, 0}, cells},
                           diffusivity.view)),
                       "scalar diffusivity is formed by the production hot path");
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double expected = (2.0 + i) / 0.5 + (3.0 + j) / 2.0;
        passed &= expect(close(diffusivity.view.unchecked({i, j, k}, 0U),
                               expected),
                         "rhoD=mu/Sc+mu_t/Sc_t exactly");
      }
    }
  }
  std::fill(diffusivity.bytes.begin(), diffusivity.bytes.end(), 91.0);
  turbulent.view.unchecked({2, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  passed &= expect(form_scalar_mass_diffusivity(
                       spec, as_const(molecular.view),
                       as_const(turbulent.view), {{0, 0, 0}, cells},
                       diffusivity.view).code == StatusCode::numerical_failure,
                   "non-finite turbulent viscosity rejects diffusivity formation");
  passed &= expect(std::all_of(diffusivity.bytes.begin(), diffusivity.bytes.end(),
                               [](double value) { return value == 91.0; }),
                   "diffusivity formation failure leaves all output bytes untouched");
  return passed;
}

bool test_outward_flux_sign_contract() {
  double heat_gradient = 0.0;
  double scalar_gradient = 0.0;
  bool passed = expect(static_cast<bool>(resolve_heat_flux_normal_gradient(
                           12.0, 3.0, heat_gradient)) &&
                           close(heat_gradient, -4.0),
                       "positive outward heat flux gives negative outward T gradient");
  passed &= expect(static_cast<bool>(resolve_scalar_flux_normal_gradient(
                       10.0, 2.0, scalar_gradient)) &&
                       close(scalar_gradient, -5.0),
                   "positive outward scalar flux gives negative outward q gradient");
  return passed;
}

bool test_conservative_global_balance_oracle() {
  // A one-cell finite-volume balance with supplied final outward mass fluxes.
  const std::array<double, 6U> outward_mdot{-0.5, 0.8, -0.1,
                                            0.2, -0.4, 0.3};
  constexpr double y = 0.25;
  long double boundary = 0.0L;
  for (const double flux : outward_mdot) {
    boundary += flux * y;
  }
  constexpr double rho_volume_old = 1.2;
  constexpr double dt = 0.1;
  const double rho_y_volume_new = rho_volume_old * y -
                                  dt * static_cast<double>(boundary);
  const double residual =
      (rho_y_volume_new - rho_volume_old * y) / dt +
      static_cast<double>(boundary);
  return expect(close(residual, 0.0),
                "N-1 transport global change equals supplied boundary flux");
}

bool output_is(const OwnedField& diagonal, const OwnedField& rhs,
               const OwnedField& residual, double value) {
  const auto matches = [value](const std::vector<double>& bytes) {
    return std::all_of(bytes.begin(), bytes.end(),
                       [value](double actual) { return actual == value; });
  };
  return matches(diagonal.bytes) && matches(rhs.bytes) &&
         matches(residual.bytes);
}

bool faces_are(const OwnedFaceField& x, const OwnedFaceField& y,
               const OwnedFaceField& z, double value) {
  const auto matches = [value](const std::vector<double>& bytes) {
    return std::all_of(bytes.begin(), bytes.end(),
                       [value](double actual) { return actual == value; });
  };
  return matches(x.bytes) && matches(y.bytes) && matches(z.bytes);
}

bool same_certificate(const EquationAssemblyCertificate& left,
                      const EquationAssemblyCertificate& right) {
  return left.plan == right.plan && left.scope == right.scope &&
         left.time == right.time && left.geometry == right.geometry &&
         left.face_flux == right.face_flux && left.state == right.state &&
         left.dt == right.dt;
}

bool test_production_species_and_passive_assembly() {
  constexpr std::int32_t n = 4;
  ProductionFixture fixture;
  bool passed = expect(make_production_fixture(n, fixture),
                       "species/passive production EquationPlanSet compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  const double dx = 1.0 / static_cast<double>(n);
  const double volume = dx * dx * dx;

  OwnedField rho = make_field(kDensity, cells, 1U, 2U, 601U);
  OwnedField species = make_field(kSpecies, cells, 1U, 2U, 602U);
  OwnedField passive = make_field(kPassive, cells, 1U, 2U, 603U);
  OwnedField molecular = make_field(20U, cells, 1U, 1U, 604U);
  OwnedField turbulent = make_field(21U, cells, 1U, 1U, 605U);
  OwnedField species_diffusivity =
      make_field(22U, cells, 1U, 1U, 606U);
  OwnedField passive_diffusivity =
      make_field(23U, cells, 1U, 1U, 607U);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const Int3 cell{i, j, k};
        const double x = (static_cast<double>(i) + 0.5) * dx;
        rho.view.unchecked(cell, 0U) = 1.0;
        species.view.unchecked(cell, 0U) = 0.20 + 0.10 * x;
        passive.view.unchecked(cell, 0U) = 0.30 - 0.04 * x;
        if (i >= -1 && i < cells.x + 1 && j >= -1 &&
            j < cells.y + 1 && k >= -1 && k < cells.z + 1) {
          molecular.view.unchecked(cell, 0U) = 2.0;
          turbulent.view.unchecked(cell, 0U) = 4.0;
          species_diffusivity.view.unchecked(cell, 0U) = 6.0;
          passive_diffusivity.view.unchecked(cell, 0U) = 2.0;
        }
      }
    }
  }
  fill_field(species_diffusivity, 6.0);
  fill_field(passive_diffusivity, 2.0);
  const ScalarEquationSpec* species_spec = fixture.equations.species().spec(0U);
  const ScalarEquationSpec* passive_spec = fixture.equations.scalars().spec(0U);
  passed &= expect(species_spec != nullptr && passive_spec != nullptr &&
                       static_cast<bool>(form_scalar_mass_diffusivity(
                           *species_spec, as_const(molecular.view),
                           as_const(turbulent.view), {{0, 0, 0}, cells},
                           species_diffusivity.view)) &&
                       static_cast<bool>(form_scalar_mass_diffusivity(
                           *passive_spec, as_const(molecular.view),
                           as_const(turbulent.view), {{0, 0, 0}, cells},
                           passive_diffusivity.view)),
                   "production Schmidt laws form both scalar diffusivities");
  passed &= expect(close(species_diffusivity.view.unchecked({1, 1, 1}, 0U),
                         2.0 / 0.5 + 4.0 / 2.0) &&
                       close(passive_diffusivity.view.unchecked({1, 1, 1}, 0U),
                             2.0 / 2.0 + 4.0 / 4.0),
                   "species and passive operators use their own Sc and Sc_t");

  FinalFluxFixture final_flux;
  ConstFaceFluxView committed;
  passed &= expect(make_linear_final_flux(fixture.equations.kernels(), cells,
                                          final_flux, committed),
                   "FinalFaceFluxWriter commits the production mass flux");
  if (!passed) {
    return false;
  }

  const PrimitiveHistory density_history{as_const(rho.view),
                                          as_const(rho.view),
                                          as_const(rho.view)};
  const PrimitiveHistory species_history{as_const(species.view),
                                          as_const(species.view),
                                          as_const(species.view)};
  const PrimitiveHistory passive_history{as_const(passive.view),
                                          as_const(passive.view),
                                          as_const(passive.view)};
  const std::array independent{species_history};
  const std::array passives{passive_history};
  EquationStateView state;
  state.density = density_history;
  state.independent_species = {independent.data(), independent.size()};
  state.passive_scalars = {passives.data(), passives.size()};
  const std::array diffusivities{as_const(species_diffusivity.view),
                                  as_const(passive_diffusivity.view)};
  EquationMaterialView material;
  material.scalar_mass_diffusivity = {diffusivities.data(),
                                      diffusivities.size()};
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 701U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.contribution_stage = 1U;
  context.face_flux = committed.revision;
  context.face_flux_authority = committed.certificate.authority();
  context.face_flux_storage = committed.certificate.storage();
  context.face_flux_revision_domain = committed.certificate.revision_domain();
  context.scope = EquationAssemblyScope::final_conservative;
  context.mass_flux = committed;
  context.provisional_mass_flux = false;

  OwnedField diagonal = make_field(30U, cells, 1U, 0U, 610U);
  OwnedField rhs = make_field(31U, cells, 1U, 0U, 611U);
  OwnedField residual = make_field(32U, cells, 1U, 0U, 612U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 613U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 614U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 615U);
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};

  EquationAssemblyCertificate species_certificate;
  passed &= expect(static_cast<bool>(assemble_species(
                       fixture.equations.species(), 0U, state, material, {},
                       context, system, species_certificate)),
                   "production N-1 species assembly succeeds");
  long double species_sum = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const double x = (static_cast<double>(i) + 0.5) * dx;
        const double q = 0.20 + 0.10 * x;
        const double expected_residual = (0.20 + 0.20 * x) * volume;
        const double expected_diagonal = 10.0 * volume + 6.0 * 6.0 * dx;
        passed &= expect(close(residual.view.unchecked(cell, 0U),
                               expected_residual) &&
                             close(diagonal.view.unchecked(cell, 0U),
                                   expected_diagonal) &&
                             close(rhs.view.unchecked(cell, 0U),
                                   expected_diagonal * q - expected_residual),
                         "N-1 species cell-integral residual/operator matches oracle");
        species_sum += residual.view.unchecked(cell, 0U);
      }
    }
  }
  long double boundary_species = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      const double q_min = 0.5 *
          (species.view.unchecked({-1, j, k}, 0U) +
           species.view.unchecked({0, j, k}, 0U));
      const double q_max = 0.5 *
          (species.view.unchecked({cells.x - 1, j, k}, 0U) +
           species.view.unchecked({cells.x, j, k}, 0U));
      boundary_species += committed.x.unchecked({cells.x, j, k}) * q_max -
                          committed.x.unchecked({0, j, k}) * q_min;
    }
  }
  passed &= expect(close(static_cast<double>(species_sum),
                         static_cast<double>(boundary_species)),
                   "global N-1 residual equals the committed boundary flux");
  passed &= expect(species_certificate.valid() &&
                       close(ax.view.unchecked({2, 1, 1}), 6.0 * dx),
                   "N-1 assembly publishes certificate and Schmidt coefficient");

  EquationAssemblyCertificate passive_certificate;
  passed &= expect(static_cast<bool>(assemble_scalar(
                       fixture.equations.scalars(), 0U, state, material, {},
                       context, system, passive_certificate)),
                   "production passive-scalar assembly succeeds");
  const Int3 probe{2, 1, 1};
  const double probe_x = (static_cast<double>(probe.x) + 0.5) * dx;
  const double probe_q = 0.30 - 0.04 * probe_x;
  const double passive_residual = (0.30 - 0.08 * probe_x) * volume;
  const double passive_diagonal = 10.0 * volume + 6.0 * 2.0 * dx;
  passed &= expect(close(residual.view.unchecked(probe, 0U),
                         passive_residual) &&
                       close(diagonal.view.unchecked(probe, 0U),
                             passive_diagonal) &&
                       close(rhs.view.unchecked(probe, 0U),
                             passive_diagonal * probe_q - passive_residual) &&
                       close(ax.view.unchecked({2, 1, 1}), 2.0 * dx) &&
                       passive_certificate.valid(),
                   "passive residual/operator uses its independent Schmidt law");

  const auto reset_outputs = [&]() {
    fill_field(diagonal, 91.0);
    fill_field(rhs, 91.0);
    fill_field(residual, 91.0);
    std::fill(ax.bytes.begin(), ax.bytes.end(), 91.0);
    std::fill(ay.bytes.begin(), ay.bytes.end(), 91.0);
    std::fill(az.bytes.begin(), az.bytes.end(), 91.0);
  };

  reset_outputs();
  EquationAssemblyCertificate rejected = passive_certificate;
  PrimitiveHistory aliased_passive = passive_history;
  aliased_passive.accepted = as_const(residual.view);
  aliased_passive.accepted.field = kPassive;
  aliased_passive.accepted.revision = passive_history.accepted.revision;
  const std::array aliased_passives{aliased_passive};
  EquationStateView aliased_state = state;
  aliased_state.passive_scalars = {aliased_passives.data(),
                                   aliased_passives.size()};
  passed &= expect(
      assemble_scalar(fixture.equations.scalars(), 0U, aliased_state,
                      material, {}, context, system, rejected)
                  .code == StatusCode::invalid_plan &&
          output_is(diagonal, rhs, residual, 91.0) &&
          faces_are(ax, ay, az, 91.0) &&
          same_certificate(rejected, passive_certificate),
      "zero-contribution scalar history/output alias rejects atomically");

  reset_outputs();
  rejected = passive_certificate;
  context.box = {{1, 0, 0}, {cells.x - 1, cells.y, cells.z}};
  passed &= expect(
      assemble_scalar(fixture.equations.scalars(), 0U, state, material, {},
                      context, system, rejected)
                  .code == StatusCode::invalid_plan &&
          output_is(diagonal, rhs, residual, 91.0) &&
          faces_are(ax, ay, az, 91.0) &&
          same_certificate(rejected, passive_certificate),
      "zero-contribution partial-box scalar assembly rejects atomically");
  context.box = {};

  reset_outputs();
  rejected = passive_certificate;
  ConstFaceFluxView stale = committed;
  ++stale.revision;
  context.face_flux = stale.revision;
  context.mass_flux = stale;
  passed &= expect(assemble_scalar(fixture.equations.scalars(), 0U, state,
                                   material, {}, context, system, rejected)
                           .code == StatusCode::invalid_plan &&
                       output_is(diagonal, rhs, residual, 91.0) &&
                       faces_are(ax, ay, az, 91.0) &&
                       rejected.plan == passive_certificate.plan,
                   "stale final authority rejects atomically");

  FaceFluxStorage workspace;
  FaceFluxView provisional;
  reset_outputs();
  rejected = passive_certificate;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, workspace)) &&
                       static_cast<bool>(workspace.workspace_view(
                           0U, 777U, provisional)),
                   "provisional flux mutation allocates");
  context.face_flux = 777U;
  context.face_flux_authority = 0U;
  context.face_flux_storage = provisional.x.storage_identity;
  context.face_flux_revision_domain = provisional.x.revision_domain;
  context.mass_flux = as_const(provisional);
  context.provisional_mass_flux = true;
  passed &= expect(assemble_scalar(fixture.equations.scalars(), 0U, state,
                                   material, {}, context, system, rejected)
                           .code == StatusCode::invalid_plan &&
                       output_is(diagonal, rhs, residual, 91.0) &&
                       faces_are(ax, ay, az, 91.0) &&
                       rejected.plan == passive_certificate.plan,
                   "provisional final flux rejects atomically");

  reset_outputs();
  rejected = species_certificate;
  context.face_flux = committed.revision;
  context.face_flux_authority = committed.certificate.authority();
  context.face_flux_storage = committed.certificate.storage();
  context.face_flux_revision_domain = committed.certificate.revision_domain();
  context.mass_flux = committed;
  context.provisional_mass_flux = false;
  species.view.unchecked({1, 1, 1}, 0U) = 1.01;
  passed &= expect(assemble_species(fixture.equations.species(), 0U, state,
                                    material, {}, context, system, rejected)
                           .code == StatusCode::numerical_failure &&
                       output_is(diagonal, rhs, residual, 91.0) &&
                       faces_are(ax, ay, az, 91.0) &&
                       rejected.plan == species_certificate.plan,
                   "invalid N-1 composition rejects before all writes");
  return passed;
}

bool test_reacting_registration_stays_out_of_scope() {
  const std::array<FieldId, 3U> declared{0U, 1U, 2U};
  ContributionRegistry registry;
  bool passed = expect(static_cast<bool>(
                           registry.configure({declared.data(), declared.size()})),
                       "inert contribution registry configures");
  ContributionSpec chemistry;
  chemistry.conserved_quantity = 0U;
  chemistry.stage = 1U;
  chemistry.explicit_source = 1U;
  chemistry.capability = ContributionCapability::chemistry;
  passed &= expect(registry.register_contribution(chemistry).code ==
                       StatusCode::invalid_plan,
                   "v0.4 species equations reject chemistry registration");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_independent_species_closure();
  passed &= test_composition_dependent_production_eos();
  passed &= test_scalar_catalog_contract();
  passed &= test_scalar_mass_diffusivity_oracle_and_atomicity();
  passed &= test_outward_flux_sign_contract();
  passed &= test_conservative_global_balance_oracle();
  passed &= test_production_species_and_passive_assembly();
  passed &= test_reacting_registration_stays_out_of_scope();
  MPI_Finalize();
  return passed ? 0 : 1;
}
