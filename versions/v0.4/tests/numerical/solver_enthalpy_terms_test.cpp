// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr double kPi = 3.141592653589793238462643383279502884;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <= 1.0e-12 * std::max(1.0, std::abs(expected));
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
  field.view.storage_identity = static_cast<StorageIdentity>(1000U + id);
  field.view.revision_domain = 9001U;
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
                     -777.0);
  field.view = {field.bytes.data(),
                extents,
                static_cast<std::size_t>(extents.x),
                static_cast<std::size_t>(extents.x) * extents.y,
                axis,
                identity,
                9002U};
  return field;
}

bool copy_face_flux(ConstFaceFluxView source, FaceFluxView destination) {
  const std::array sources{source.x, source.y, source.z};
  const std::array destinations{destination.x, destination.y, destination.z};
  for (std::size_t axis = 0U; axis < sources.size(); ++axis) {
    if (sources[axis].extents.x != destinations[axis].extents.x ||
        sources[axis].extents.y != destinations[axis].extents.y ||
        sources[axis].extents.z != destinations[axis].extents.z) {
      return false;
    }
    for (std::int32_t z = 0; z < sources[axis].extents.z; ++z) {
      for (std::int32_t y = 0; y < sources[axis].extents.y; ++y) {
        for (std::int32_t x = 0; x < sources[axis].extents.x; ++x) {
          const Int3 face{x, y, z};
          destinations[axis].unchecked(face) = sources[axis].unchecked(face);
        }
      }
    }
  }
  return true;
}

CartesianMeshSpec mesh_spec(std::int32_t n, bool stretched = false) {
  CartesianMeshSpec mesh;
  mesh.kind = stretched ? GeometryKind::tensor_stretched
                        : GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  const double inverse_n = 1.0 / static_cast<double>(n);
  if (stretched) {
    mesh.has_base_spacing = true;
    mesh.base_spacing = {1.08 * inverse_n, 1.08 * inverse_n,
                         1.08 * inverse_n};
    mesh.minimum_spacing = {0.92 * inverse_n, 0.92 * inverse_n,
                            0.92 * inverse_n};
    mesh.max_growth_ratio = 1.0 + 0.8 * inverse_n;
    mesh.focus_regions.push_back(
        {{0.35, 0.35, 0.35}, {0.65, 0.65, 0.65},
         {0.96 * inverse_n, 0.96 * inverse_n, 0.96 * inverse_n}});
  } else {
    mesh.minimum_spacing = {inverse_n, inverse_n, inverse_n};
    mesh.max_growth_ratio = 1.0;
  }
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(n) * n * n;
  mesh.limits.max_memory_bytes_per_rank = 1U << 28U;
  return mesh;
}

SpeciesThermophysicalSpec air() {
  SpeciesThermophysicalSpec species;
  species.stable_name = "air";
  species.molecular_weight = 28.96546;
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
  value.species.push_back(air());
  return value;
}

struct Fixture {
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

bool make_fixture(std::int32_t n, Fixture& out, bool stretched = false) {
  const CartesianMeshSpec mesh = mesh_spec(n, stretched);
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14e17111U;
  model.pressure_reference = PressureReferenceKind::closed_mass;
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
  if (!registry.require_field("rho", 1U, 2U, density) || density != 0U ||
      !registry.require_field("U", 3U, 2U, velocity) || velocity != 1U ||
      !registry.require_field("pi", 1U, 2U, pressure) || pressure != 2U ||
      !registry.require_field("h", 1U, 2U, enthalpy) || enthalpy != 3U ||
      !registry.require_field("T", 1U, 2U, temperature) || temperature != 4U ||
      !CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes,
                                 out.time)) {
    return false;
  }
  const ThermophysicalSpec thermo = thermo_spec();
  if (!ThermodynamicsPlan::compile(thermo, {}, out.thermodynamics) ||
      !TransportPlan::compile(thermo, out.thermodynamics, out.transport)) {
    return false;
  }
  const std::array<FieldId, 8U> declared{0U, 1U, 2U, 3U,
                                         4U, 5U, 6U, 7U};
  if (!out.contributions.configure({declared.data(), declared.size()}) ||
      !out.contributions.freeze()) {
    return false;
  }
  EquationPlanSpec spec;
  spec.density = 0U;
  spec.velocity = 1U;
  spec.pressure_perturbation = 2U;
  spec.enthalpy = 3U;
  spec.temperature = 4U;
  spec.effective_viscosity = 5U;
  spec.pressure_compressibility = 6U;
  spec.velocity_gradient = 7U;
  spec.pressure_reference = PressureReferenceKind::closed_mass;
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

bool make_final_flux(const CartesianKernelPlan& kernels, Int3 cells,
                     double velocity_x, FinalFluxFixture& fixture,
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
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        rho.view.unchecked({i, j, k}, 0U) = 1.0;
        velocity.view.unchecked({i, j, k}, 0U) = velocity_x;
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

double periodic_center(const AxisMetrics& axis, std::int32_t index,
                       std::int32_t cells) {
  const std::int32_t wrapped = (index % cells + cells) % cells;
  const std::int32_t periods = (index - wrapped) / cells;
  return axis.centres().data[static_cast<std::size_t>(wrapped)] +
         static_cast<double>(periods) *
             (axis.faces().data[static_cast<std::size_t>(cells)] -
              axis.faces().data[0U]);
}

double local_volume(const CartesianGeometryPlan& geometry, Int3 cell) {
  return geometry.x().widths().data[static_cast<std::size_t>(cell.x)] *
         geometry.y().widths().data[static_cast<std::size_t>(cell.y)] *
         geometry.z().widths().data[static_cast<std::size_t>(cell.z)];
}

double run_enthalpy_mms(std::int32_t n, bool stretched,
                        bool pressure_only) {
  Fixture fixture;
  if (!make_fixture(n, fixture, stretched)) {
    return HUGE_VAL;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField rho = make_field(0U, cells, 1U, 2U, 2001U);
  OwnedField velocity = make_field(1U, cells, 3U, 2U, 2002U);
  OwnedField pi = make_field(2U, cells, 1U, 1U, 2003U);
  OwnedField h = make_field(3U, cells, 1U, 2U, 2004U);
  OwnedField temperature = make_field(4U, cells, 1U, 1U, 2005U);
  OwnedField mu = make_field(5U, cells, 1U, 1U, 2006U);
  OwnedField lambda = make_field(6U, cells, 1U, 1U, 2007U);
  OwnedField cp = make_field(7U, cells, 1U, 1U, 2008U);
  OwnedField lambda_over_cp = make_field(8U, cells, 1U, 1U, 2009U);
  OwnedField gradients = make_field(7U, cells, 9U, 1U, 2010U);
  OwnedField diagonal = make_field(30U, cells, 1U, 0U, 2011U);
  OwnedField rhs = make_field(31U, cells, 1U, 0U, 2012U);
  OwnedField residual = make_field(32U, cells, 1U, 0U, 2013U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 2014U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 2015U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 2016U);

  constexpr double velocity_x = 0.4;
  constexpr double pressure_amplitude = 0.3;
  constexpr double temperature_amplitude = 5.0;
  constexpr double heat_capacity = 4.0;
  constexpr double conductivity = 2.0;
  constexpr double enthalpy_amplitude =
      heat_capacity * temperature_amplitude;
  constexpr double wave = 2.0 * kPi;
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const double x = periodic_center(fixture.geometry.x(), i, n);
        const Int3 cell{i, j, k};
        rho.view.unchecked(cell, 0U) = 1.0;
        velocity.view.unchecked(cell, 0U) = velocity_x;
        velocity.view.unchecked(cell, 1U) = 0.0;
        velocity.view.unchecked(cell, 2U) = 0.0;
        h.view.unchecked(cell, 0U) =
            300000.0 +
            (pressure_only ? 0.0
                           : enthalpy_amplitude * std::sin(wave * x));
        if (i >= -1 && i < cells.x + 1 && j >= -1 &&
            j < cells.y + 1 && k >= -1 && k < cells.z + 1) {
          pi.view.unchecked(cell, 0U) =
              pressure_amplitude * std::sin(wave * x);
          temperature.view.unchecked(cell, 0U) =
              300.0 +
              (pressure_only
                   ? 0.0
                   : temperature_amplitude * std::sin(wave * x));
          mu.view.unchecked(cell, 0U) = 0.0;
          lambda.view.unchecked(cell, 0U) = conductivity;
          cp.view.unchecked(cell, 0U) = heat_capacity;
          lambda_over_cp.view.unchecked(cell, 0U) =
              conductivity / heat_capacity;
        }
      }
    }
  }
  std::fill(gradients.bytes.begin(), gradients.bytes.end(), 0.0);

  FinalFluxFixture flux_fixture;
  ConstFaceFluxView flux;
  if (!make_final_flux(fixture.equations.kernels(), cells, velocity_x,
                       flux_fixture, flux)) {
    return HUGE_VAL;
  }
  const PrimitiveHistory density{as_const(rho.view), as_const(rho.view),
                                 as_const(rho.view)};
  const PrimitiveHistory vector_velocity{
      as_const(velocity.view), as_const(velocity.view),
      as_const(velocity.view)};
  const PrimitiveHistory pressure{as_const(pi.view), as_const(pi.view),
                                  as_const(pi.view)};
  const PrimitiveHistory enthalpy{as_const(h.view), as_const(h.view),
                                  as_const(h.view)};
  const PrimitiveHistory thermal{as_const(temperature.view),
                                 as_const(temperature.view),
                                 as_const(temperature.view)};
  EquationStateView state{density, vector_velocity, pressure, enthalpy,
                          thermal, 101325.0, {}, {}, 101325.0, 101325.0};
  EquationMaterialView material;
  material.effective_viscosity = as_const(mu.view);
  material.thermal_conductivity = as_const(lambda.view);
  material.heat_capacity = as_const(cp.view);
  material.enthalpy_diffusivity = as_const(lambda_over_cp.view);
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 2020U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = flux.revision;
  context.face_flux_authority = flux.certificate.authority();
  context.face_flux_storage = flux.certificate.storage();
  context.face_flux_revision_domain = flux.certificate.revision_domain();
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::final_conservative;
  context.mass_flux = flux;
  context.provisional_mass_flux = false;
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};
  EquationAssemblyCertificate certificate;
  if (!assemble_enthalpy(fixture.equations.enthalpy(), state, material,
                         as_const(gradients.view), {}, context, system,
                         certificate)) {
    return HUGE_VAL;
  }

  long double error = 0.0L;
  long double measure = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const double x = fixture.geometry.x().centres().data[
            static_cast<std::size_t>(i)];
        const double c = std::cos(wave * x);
        const double s = std::sin(wave * x);
        const double pressure_work =
            velocity_x * pressure_amplitude * wave * c;
        const double expected_density =
            pressure_only
                ? -pressure_work
                : velocity_x * enthalpy_amplitude * wave * c -
                      pressure_work +
                      conductivity * temperature_amplitude * wave * wave * s;
        const double expected = expected_density *
                                local_volume(fixture.geometry, cell);
        const double difference =
            residual.view.unchecked(cell, 0U) - expected;
        error += difference * difference;
        measure += expected * expected;
      }
    }
  }
  return std::sqrt(static_cast<double>(error / measure));
}

double observed_order(double coarse, double fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

bool test_enthalpy_and_pressure_work_mms_orders() {
  bool passed = true;
  for (const bool pressure_only : {false, true}) {
    for (const bool stretched : {false, true}) {
      const double e16 = run_enthalpy_mms(16, stretched, pressure_only);
      const double e32 = run_enthalpy_mms(32, stretched, pressure_only);
      const double e64 = run_enthalpy_mms(64, stretched, pressure_only);
      passed &= expect(std::isfinite(e16) && e16 > e32 && e32 > e64,
                       pressure_only
                           ? "production pressure-work MMS converges"
                           : "production enthalpy MMS converges");
      passed &= expect(observed_order(e16, e32) >= 1.8 &&
                           observed_order(e32, e64) >= 1.8,
                       pressure_only
                           ? "uniform/stretched pressure-work order >= 1.8"
                           : "uniform/stretched enthalpy order >= 1.8");
    }
  }
  return passed;
}

bool test_pressure_material_derivative_oracle() {
  // At the cell centre, p=7+2t+3x, U=(4,0,0): Dp/Dt=2+12=14.
  const PressureWorkPoint point{7.0, 7.0 - 0.1 * 2.0,
                                7.0 - 0.2 * 2.0, {3.0, 0.0, 0.0},
                                {4.0, 0.0, 0.0}};
  const BdfCoefficients bdf{15.0, -20.0, 5.0, 2U};
  double material = 0.0;
  bool passed = expect(static_cast<bool>(evaluate_pressure_material_derivative(
                           bdf, point, material)),
                       "pressure material derivative evaluates");
  passed &= expect(close(material, 14.0), "Dp/Dt uses BDF(p_abs)+U dot grad(p_abs)");

  // div(pU)=U.grad(p)+p div(U). This mutation contributes 35 and must not be
  // accepted when p=7 and div(U)=5.
  PressureWorkPoint divergent = point;
  divergent.velocity_divergence = 5.0;
  double mutation = 0.0;
  passed &= expect(static_cast<bool>(evaluate_pressure_material_derivative(
                           bdf, divergent, mutation)) && close(mutation, 14.0),
                   "pressure work kills the div(pU) mutation");
  passed &= expect(!close(mutation, 49.0),
                   "pressure work excludes p_abs times div(U)");
  return passed;
}

bool test_energy_term_signs_and_units() {
  const std::array<EnthalpyTermPoint, 5U> isolated{{
      {EnthalpyTermKind::unsteady, 2.5, 0.0, 0.0, 0.0},
      {EnthalpyTermKind::advection, -3.0, 0.0, 0.0, 0.0},
      {EnthalpyTermKind::pressure_work, 7.0, 0.0, 0.0, 0.0},
      {EnthalpyTermKind::conduction, -11.0, 0.0, 0.0, 0.0},
      {EnthalpyTermKind::viscous_dissipation, 13.0, 0.0, 0.0, 0.0},
  }};
  // Residual contract:
  // BDF(rho h)+div(mdot h)-Dp/Dt-div(lambda grad T)-Phi-source.
  const std::array<double, 5U> expected{2.5, -3.0, -7.0, 11.0, -13.0};
  bool passed = true;
  for (std::size_t i = 0U; i < isolated.size(); ++i) {
    double residual = 0.0;
    passed &= expect(static_cast<bool>(assemble_isolated_enthalpy_term(
                         isolated[i], 0.125, residual)),
                     "isolated energy term evaluates");
    passed &= expect(close(residual, expected[i] * 0.125),
                     "isolated energy term has cell-integral sign and units");
  }
  EnthalpyTermPoint source{EnthalpyTermKind::source, 17.0, 0.0, 0.0, 0.0};
  double source_residual = 0.0;
  passed &= expect(static_cast<bool>(assemble_isolated_enthalpy_term(
                       source, 0.125, source_residual)) &&
                       close(source_residual, -17.0 * 0.125),
                   "enthalpy source is subtracted from the residual");
  return passed;
}

bool test_viscous_dissipation_uses_complete_tau() {
  VelocityGradient gradient;
  gradient.value = {1.0, 2.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, -2.0};
  double phi = 0.0;
  bool passed = expect(static_cast<bool>(newtonian_viscous_dissipation(
                           gradient, 0.5, phi)),
                       "complete Newtonian dissipation evaluates");
  // Independent tensor oracle: tau=mu*(G+G^T-2/3 tr(G) I), Phi=tau:G.
  long double oracle = 0.0L;
  const long double trace = 1.0L + 4.0L - 2.0L;
  for (std::size_t i = 0U; i < 3U; ++i) {
    for (std::size_t j = 0U; j < 3U; ++j) {
      const long double gij = gradient.value[3U * i + j];
      const long double gji = gradient.value[3U * j + i];
      const long double tau =
          0.5L * (gij + gji - (i == j ? 2.0L * trace / 3.0L : 0.0L));
      oracle += tau * gij;
    }
  }
  passed &= expect(close(phi, static_cast<double>(oracle)),
                   "Phi is tau:grad(U), not a componentwise Laplacian proxy");
  return passed;
}

bool same_certificate(const EquationAssemblyCertificate& left,
                      const EquationAssemblyCertificate& right) {
  return left.plan == right.plan && left.scope == right.scope &&
         left.time == right.time && left.geometry == right.geometry &&
         left.face_flux == right.face_flux && left.state == right.state &&
         left.dt == right.dt;
}

bool test_production_enthalpy_assembly_oracle() {
  constexpr std::int32_t n = 6;
  constexpr double velocity_x = 4.0;
  constexpr double enthalpy_slope = 3.0;
  constexpr double conductivity = 2.0;
  constexpr double heat_capacity = 4.0;
  constexpr double viscosity = 0.5;
  const double spacing = 1.0 / static_cast<double>(n);
  const double volume = spacing * spacing * spacing;
  const Int3 oracle_cell{3, 3, 3};

  Fixture fixture;
  bool passed = expect(make_fixture(n, fixture),
                       "enthalpy production fixture compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField rho_trial = make_field(0U, cells, 1U, 2U, 601U);
  OwnedField rho_accepted = make_field(0U, cells, 1U, 2U, 602U);
  OwnedField rho_previous = make_field(0U, cells, 1U, 2U, 603U);
  OwnedField velocity = make_field(1U, cells, 3U, 2U, 604U);
  OwnedField pi_trial = make_field(2U, cells, 1U, 2U, 605U);
  OwnedField pi_accepted = make_field(2U, cells, 1U, 2U, 606U);
  OwnedField pi_previous = make_field(2U, cells, 1U, 2U, 607U);
  OwnedField h_trial = make_field(3U, cells, 1U, 2U, 608U);
  OwnedField h_accepted = make_field(3U, cells, 1U, 2U, 609U);
  OwnedField h_previous = make_field(3U, cells, 1U, 2U, 610U);
  OwnedField temperature = make_field(4U, cells, 1U, 2U, 611U);
  OwnedField mu = make_field(5U, cells, 1U, 2U, 612U);
  OwnedField lambda = make_field(6U, cells, 1U, 2U, 613U);
  OwnedField cp = make_field(7U, cells, 1U, 2U, 614U);
  OwnedField lambda_over_cp = make_field(8U, cells, 1U, 2U, 615U);
  OwnedField gradients = make_field(7U, cells, 9U, 1U, 616U);
  OwnedField diagonal = make_field(30U, cells, 1U, 0U, 617U);
  OwnedField rhs = make_field(31U, cells, 1U, 0U, 618U);
  OwnedField residual = make_field(32U, cells, 1U, 0U, 619U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 620U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 621U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 622U);

  FinalFluxFixture stationary_flux;
  FinalFluxFixture moving_flux;
  ConstFaceFluxView stationary;
  ConstFaceFluxView moving;
  passed &= expect(make_final_flux(fixture.equations.kernels(), cells, 0.0,
                                   stationary_flux, stationary),
                   "stationary final flux commits");
  passed &= expect(make_final_flux(fixture.equations.kernels(), cells,
                                   velocity_x, moving_flux, moving),
                   "advective final flux commits");
  if (!passed) {
    return false;
  }

  const PrimitiveHistory density{as_const(rho_trial.view),
                                 as_const(rho_accepted.view),
                                 as_const(rho_previous.view)};
  const PrimitiveHistory vector_velocity{as_const(velocity.view),
                                         as_const(velocity.view),
                                         as_const(velocity.view)};
  const PrimitiveHistory pressure{as_const(pi_trial.view),
                                  as_const(pi_accepted.view),
                                  as_const(pi_previous.view)};
  const PrimitiveHistory enthalpy{as_const(h_trial.view),
                                  as_const(h_accepted.view),
                                  as_const(h_previous.view)};
  const PrimitiveHistory thermal{as_const(temperature.view),
                                 as_const(temperature.view),
                                 as_const(temperature.view)};
  EquationStateView state{density, vector_velocity, pressure, enthalpy,
                          thermal, 100.0, {}, {}, 100.0, 100.0};
  EquationMaterialView material;
  material.effective_viscosity = as_const(mu.view);
  material.thermal_conductivity = as_const(lambda.view);
  material.heat_capacity = as_const(cp.view);
  material.enthalpy_diffusivity = as_const(lambda_over_cp.view);
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 700U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::final_conservative;
  context.provisional_mass_flux = false;

  auto select_flux = [&](ConstFaceFluxView flux) {
    context.face_flux = flux.revision;
    context.face_flux_authority = flux.certificate.authority();
    context.face_flux_storage = flux.certificate.storage();
    context.face_flux_revision_domain = flux.certificate.revision_domain();
    context.mass_flux = flux;
  };
  auto reset_fields = [&]() {
    for (std::int32_t k = -2; k < cells.z + 2; ++k) {
      for (std::int32_t j = -2; j < cells.y + 2; ++j) {
        for (std::int32_t i = -2; i < cells.x + 2; ++i) {
          const Int3 cell{i, j, k};
          rho_trial.view.unchecked(cell, 0U) = 1.0;
          rho_accepted.view.unchecked(cell, 0U) = 1.0;
          rho_previous.view.unchecked(cell, 0U) = 1.0;
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            velocity.view.unchecked(cell, component) = 0.0;
          }
          pi_trial.view.unchecked(cell, 0U) = 0.0;
          pi_accepted.view.unchecked(cell, 0U) = 0.0;
          pi_previous.view.unchecked(cell, 0U) = 0.0;
          h_trial.view.unchecked(cell, 0U) = 300000.0;
          h_accepted.view.unchecked(cell, 0U) = 300000.0;
          h_previous.view.unchecked(cell, 0U) = 300000.0;
          temperature.view.unchecked(cell, 0U) = 300.0;
          mu.view.unchecked(cell, 0U) = 0.0;
          lambda.view.unchecked(cell, 0U) = conductivity;
          cp.view.unchecked(cell, 0U) = heat_capacity;
          lambda_over_cp.view.unchecked(cell, 0U) =
              conductivity / heat_capacity;
        }
      }
    }
    std::fill(gradients.bytes.begin(), gradients.bytes.end(), 0.0);
    std::fill(diagonal.bytes.begin(), diagonal.bytes.end(), -11.0);
    std::fill(rhs.bytes.begin(), rhs.bytes.end(), -12.0);
    std::fill(residual.bytes.begin(), residual.bytes.end(), -13.0);
    std::fill(ax.bytes.begin(), ax.bytes.end(), -14.0);
    std::fill(ay.bytes.begin(), ay.bytes.end(), -15.0);
    std::fill(az.bytes.begin(), az.bytes.end(), -16.0);
    state.pressure_reference = 100.0;
    state.accepted_pressure_reference = 100.0;
    state.previous_pressure_reference = 100.0;
    select_flux(stationary);
  };
  auto assemble_and_check = [&](double expected,
                                std::string_view description) {
    EquationAssemblyCertificate certificate;
    const Status status = assemble_enthalpy(
        fixture.equations.enthalpy(), state, material,
        as_const(gradients.view), {}, context, system, certificate);
    passed &= expect(static_cast<bool>(status), description);
    passed &= expect(close(residual.view.unchecked(oracle_cell, 0U), expected),
                     "production term matches independent cell-integral oracle");
    passed &= expect(certificate.valid() &&
                         certificate.scope ==
                             EquationAssemblyScope::final_conservative,
                     "production term publishes final certificate");
    return certificate;
  };

  // BDF(rho*h): rho=1 and h changes by 2 J/kg in 0.1 s, hence
  // +20 W/m^3 before multiplying by the cell volume.
  reset_fields();
  std::fill(h_trial.bytes.begin(), h_trial.bytes.end(), 300002.0);
  EquationAssemblyCertificate certificate = assemble_and_check(
      20.0 * volume, "unsteady enthalpy isolates in production assembler");

  // div(mdot*h): Ux=4 m/s and dh/dx=3 J/(kg m), giving +12 W/m^3.
  reset_fields();
  select_flux(moving);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const Int3 cell{i, j, k};
        const double x = (static_cast<double>(i) + 0.5) * spacing;
        const double value = 300000.0 + enthalpy_slope * x;
        h_trial.view.unchecked(cell, 0U) = value;
        h_accepted.view.unchecked(cell, 0U) = value;
        h_previous.view.unchecked(cell, 0U) = value;
        velocity.view.unchecked(cell, 0U) = velocity_x;
      }
    }
  }
  certificate = assemble_and_check(
      velocity_x * enthalpy_slope * volume,
      "advective enthalpy isolates in production assembler");

  // Dp/Dt: p_abs changes by 0.2 Pa in 0.1 s, and energy subtracts it.
  reset_fields();
  state.accepted_pressure_reference = 99.8;
  state.previous_pressure_reference = 99.8;
  certificate = assemble_and_check(
      -2.0 * volume, "pressure material derivative isolates in production assembler");

  // -div(lambda*grad(T)): T=300+5*x^2 and lambda=2, hence -20 W/m^3.
  reset_fields();
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * spacing;
        temperature.view.unchecked({i, j, k}, 0U) = 300.0 + 5.0 * x * x;
      }
    }
  }
  certificate = assemble_and_check(
      -conductivity * 10.0 * volume,
      "temperature conduction isolates in production assembler");
  const double transmissibility =
      (conductivity / heat_capacity) * spacing;
  passed &= expect(close(ax.view.unchecked({4, 3, 3}), transmissibility) &&
                       close(diagonal.view.unchecked(oracle_cell, 0U),
                             10.0 * volume + 6.0 * transmissibility),
                   "temperature residual and frozen dT/dh operator stay distinct");

  // -Phi with the complete Newtonian tau:grad(U) tensor.
  reset_fields();
  const std::array<double, 9U> gradient_values{
      1.0, 2.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, -2.0};
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        mu.view.unchecked({i, j, k}, 0U) = viscosity;
        for (std::uint8_t component = 0U; component < 9U; ++component) {
          gradients.view.unchecked({i, j, k}, component) =
              gradient_values[component];
        }
      }
    }
  }
  VelocityGradient point_gradient;
  point_gradient.value = gradient_values;
  double dissipation = 0.0;
  passed &= expect(static_cast<bool>(newtonian_viscous_dissipation(
                           point_gradient, viscosity, dissipation)),
                   "viscous oracle evaluates independently");
  certificate = assemble_and_check(
      -dissipation * volume,
      "viscous dissipation isolates in production assembler");

  const std::vector<double> full_diagonal = diagonal.bytes;
  const std::vector<double> full_rhs = rhs.bytes;
  const std::vector<double> full_residual = residual.bytes;
  const std::vector<double> full_ax = ax.bytes;
  const std::vector<double> full_ay = ay.bytes;
  const std::vector<double> full_az = az.bytes;

  // C1/C2 must evaluate the same conservative residual against the current
  // target-layer flux before that flux is eligible for final publication.
  const EquationAssemblyContext final_context = context;
  FaceFluxStorage target_flux_storage;
  FaceFluxView target_flux;
  constexpr RevisionToken target_flux_revision = 1701U;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, target_flux_storage)) &&
          static_cast<bool>(target_flux_storage.workspace_view(
              0U, target_flux_revision, target_flux)) &&
          copy_face_flux(final_context.mass_flux, target_flux),
      "target-coupled flux copies the current conservative face values");
  context.scope = EquationAssemblyScope::target_coupled;
  context.face_flux = target_flux_revision;
  context.face_flux_authority = 0U;
  context.face_flux_storage = 0U;
  context.face_flux_revision_domain = 0U;
  context.mass_flux = as_const(target_flux);
  context.provisional_mass_flux = false;
  std::fill(diagonal.bytes.begin(), diagonal.bytes.end(), -61.0);
  std::fill(rhs.bytes.begin(), rhs.bytes.end(), -62.0);
  std::fill(residual.bytes.begin(), residual.bytes.end(), -63.0);
  std::fill(ax.bytes.begin(), ax.bytes.end(), -64.0);
  std::fill(ay.bytes.begin(), ay.bytes.end(), -65.0);
  std::fill(az.bytes.begin(), az.bytes.end(), -66.0);
  EquationAssemblyCertificate target_certificate;
  const Status target_status = assemble_enthalpy(
      fixture.equations.enthalpy(), state, material,
      as_const(gradients.view), {}, context, system, target_certificate);
  if (!target_status) {
    std::cerr << "target-coupled status code="
              << static_cast<unsigned>(target_status.code)
              << " detail=" << target_status.detail << '\n';
  }
  passed &= expect(static_cast<bool>(target_status),
                   "target-coupled scope accepts an attempt-local flux");
  passed &= expect(diagonal.bytes == full_diagonal && rhs.bytes == full_rhs &&
                       residual.bytes == full_residual && ax.bytes == full_ax &&
                       ay.bytes == full_ay && az.bytes == full_az &&
                       target_certificate.valid() &&
                       target_certificate.scope ==
                           EquationAssemblyScope::target_coupled,
                   "target-coupled assembly equals the final conservative residual");
  context = final_context;
  const KernelBox lower{{0, 0, 0}, {cells.x / 2, cells.y, cells.z}};
  const KernelBox upper{{cells.x / 2, 0, 0},
                        {cells.x - cells.x / 2, cells.y, cells.z}};

  std::fill(diagonal.bytes.begin(), diagonal.bytes.end(), -71.0);
  std::fill(rhs.bytes.begin(), rhs.bytes.end(), -72.0);
  std::fill(residual.bytes.begin(), residual.bytes.end(), -73.0);
  std::fill(ax.bytes.begin(), ax.bytes.end(), -74.0);
  std::fill(ay.bytes.begin(), ay.bytes.end(), -75.0);
  std::fill(az.bytes.begin(), az.bytes.end(), -76.0);
  AssemblyEpoch epoch;
  EquationAssemblyCertificate tiled = certificate;
  passed &= expect(static_cast<bool>(epoch.begin(context, system)) &&
                       static_cast<bool>(assemble_tile(
                           epoch, fixture.equations.enthalpy(), state,
                           material, as_const(gradients.view), {}, upper)) &&
                       epoch.completed_tiles() == 1U &&
                       same_certificate(tiled, certificate) &&
                       static_cast<bool>(assemble_tile(
                           epoch, fixture.equations.enthalpy(), state,
                           material, as_const(gradients.view), {}, lower)) &&
                       static_cast<bool>(epoch.finalize(tiled)),
                   "AssemblyEpoch publishes only after complete tiled coverage");
  passed &= expect(diagonal.bytes == full_diagonal && rhs.bytes == full_rhs &&
                       residual.bytes == full_residual &&
                       ax.bytes == full_ax && ay.bytes == full_ay &&
                       az.bytes == full_az &&
                       same_certificate(tiled, certificate),
                   "reverse-order tiles equal full-domain production assembly");

  std::fill(diagonal.bytes.begin(), diagonal.bytes.end(), -81.0);
  std::fill(rhs.bytes.begin(), rhs.bytes.end(), -82.0);
  std::fill(residual.bytes.begin(), residual.bytes.end(), -83.0);
  std::fill(ax.bytes.begin(), ax.bytes.end(), -84.0);
  std::fill(ay.bytes.begin(), ay.bytes.end(), -85.0);
  std::fill(az.bytes.begin(), az.bytes.end(), -86.0);
  AssemblyEpoch incomplete;
  EquationAssemblyCertificate incomplete_certificate = certificate;
  passed &= expect(static_cast<bool>(incomplete.begin(context, system)) &&
                       static_cast<bool>(assemble_tile(
                           incomplete, fixture.equations.enthalpy(), state,
                           material, as_const(gradients.view), {}, lower)) &&
                       incomplete.finalize(incomplete_certificate).code ==
                           StatusCode::invalid_plan &&
                       same_certificate(incomplete_certificate, certificate),
                   "incomplete tile coverage cannot publish a certificate");

  AssemblyEpoch changed_dependency;
  EquationStateView changed_state = state;
  changed_state.accepted_pressure_reference += 1.0;
  passed &= expect(static_cast<bool>(changed_dependency.begin(context,
                                                              system)) &&
                       static_cast<bool>(assemble_tile(
                           changed_dependency, fixture.equations.enthalpy(),
                           state, material, as_const(gradients.view), {},
                           lower)) &&
                       assemble_tile(changed_dependency,
                                     fixture.equations.enthalpy(),
                                     changed_state, material,
                                     as_const(gradients.view), {}, upper)
                               .code == StatusCode::invalid_plan,
                   "tile dependency mutation poisons the epoch");

  auto seed_outputs = [&]() {
    std::fill(diagonal.bytes.begin(), diagonal.bytes.end(), 41.0);
    std::fill(rhs.bytes.begin(), rhs.bytes.end(), 42.0);
    std::fill(residual.bytes.begin(), residual.bytes.end(), 43.0);
    std::fill(ax.bytes.begin(), ax.bytes.end(), 44.0);
    std::fill(ay.bytes.begin(), ay.bytes.end(), 45.0);
    std::fill(az.bytes.begin(), az.bytes.end(), 46.0);
  };
  auto regular_outputs_unchanged = [&]() {
    return std::all_of(diagonal.bytes.begin(), diagonal.bytes.end(),
                       [](double value) { return value == 41.0; }) &&
           std::all_of(rhs.bytes.begin(), rhs.bytes.end(),
                       [](double value) { return value == 42.0; }) &&
           std::all_of(residual.bytes.begin(), residual.bytes.end(),
                       [](double value) { return value == 43.0; }) &&
           std::all_of(ax.bytes.begin(), ax.bytes.end(),
                       [](double value) { return value == 44.0; }) &&
           std::all_of(ay.bytes.begin(), ay.bytes.end(),
                       [](double value) { return value == 45.0; }) &&
           std::all_of(az.bytes.begin(), az.bytes.end(),
                       [](double value) { return value == 46.0; });
  };

  // The full input allocation, not only the first cell, is protected when an
  // output view aliases trial enthalpy.
  seed_outputs();
  const std::vector<double> enthalpy_snapshot = h_trial.bytes;
  EquationSystemView aliased = system;
  aliased.residual = h_trial.view;
  EquationAssemblyCertificate rejected = certificate;
  const Status alias_status = assemble_enthalpy(
      fixture.equations.enthalpy(), state, material,
      as_const(gradients.view), {}, context, aliased, rejected);
  passed &= expect(alias_status.code == StatusCode::invalid_plan,
                   "enthalpy rejects input/output aliasing");
  passed &= expect(h_trial.bytes == enthalpy_snapshot &&
                       std::all_of(diagonal.bytes.begin(), diagonal.bytes.end(),
                                   [](double value) { return value == 41.0; }) &&
                       std::all_of(rhs.bytes.begin(), rhs.bytes.end(),
                                   [](double value) { return value == 42.0; }) &&
                       std::all_of(ax.bytes.begin(), ax.bytes.end(),
                                   [](double value) { return value == 44.0; }) &&
                       std::all_of(ay.bytes.begin(), ay.bytes.end(),
                                   [](double value) { return value == 45.0; }) &&
                       std::all_of(az.bytes.begin(), az.bytes.end(),
                                   [](double value) { return value == 46.0; }) &&
                       same_certificate(rejected, certificate),
                   "alias rejection preserves every output and certificate");

  seed_outputs();
  context.box = {{1, 0, 0}, {cells.x - 1, cells.y, cells.z}};
  rejected = certificate;
  const Status partial_status = assemble_enthalpy(
      fixture.equations.enthalpy(), state, material,
      as_const(gradients.view), {}, context, system, rejected);
  passed &= expect(partial_status.code == StatusCode::invalid_plan,
                   "enthalpy rejects partial equation boxes");
  passed &= expect(regular_outputs_unchanged() &&
                       same_certificate(rejected, certificate),
                   "partial-box rejection preserves all outputs and certificate");

  // An attempt-local target flux must never cross the final-energy boundary,
  // even though it is not a momentum-predictor/provisional flux.
  context.box = {};
  seed_outputs();
  context.face_flux = target_flux_revision;
  context.face_flux_authority = 0U;
  context.face_flux_storage = 0U;
  context.face_flux_revision_domain = 0U;
  context.mass_flux = as_const(target_flux);
  context.provisional_mass_flux = false;
  rejected = certificate;
  passed &= expect(assemble_enthalpy(
                       fixture.equations.enthalpy(), state, material,
                       as_const(gradients.view), {}, context, system,
                       rejected).code == StatusCode::invalid_plan,
                   "final enthalpy rejects an uncommitted target flux");
  passed &= expect(regular_outputs_unchanged() &&
                       same_certificate(rejected, certificate),
                   "uncommitted target-flux rejection is output-atomic");

  // A momentum-predictor/provisional flux must also never cross the final
  // energy authority boundary.
  context.box = {};
  FaceFluxStorage workspace;
  FaceFluxView provisional;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, workspace)) &&
                       static_cast<bool>(workspace.workspace_view(
                           0U, 777U, provisional)),
                   "provisional mutation allocates");
  seed_outputs();
  context.face_flux = 777U;
  context.mass_flux = as_const(provisional);
  context.provisional_mass_flux = true;
  rejected = certificate;
  passed &= expect(assemble_enthalpy(
                       fixture.equations.enthalpy(), state, material,
                       as_const(gradients.view), {}, context,
                       system, rejected).code == StatusCode::invalid_plan,
                   "enthalpy rejects provisional flux at the final authority boundary");
  passed &= expect(regular_outputs_unchanged() &&
                       same_certificate(rejected, certificate),
                   "final-flux authority failure is output-atomic");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_pressure_material_derivative_oracle();
  passed &= test_enthalpy_and_pressure_work_mms_orders();
  passed &= test_energy_term_signs_and_units();
  passed &= test_viscous_dissipation_uses_complete_tau();
  passed &= test_production_enthalpy_assembly_oracle();
  MPI_Finalize();
  return passed ? 0 : 1;
}
