// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr FieldId kCompressibility = 6U;
constexpr StageId kClosedMassServiceStage = 37U;
constexpr double kSentinel = -913.25;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             1.0e-13 * std::max(1.0, std::abs(expected));
}

struct OwnedField {
  std::vector<double> values;
  FieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, RevisionToken revision,
                      StorageIdentity storage) {
  OwnedField field;
  const std::size_t count = static_cast<std::size_t>(cells.x) *
                            static_cast<std::size_t>(cells.y) *
                            static_cast<std::size_t>(cells.z);
  field.values.assign(count, 0.0);
  field.view.base = field.values.data();
  field.view.interior = cells;
  field.view.components = 1U;
  field.view.stride_y = static_cast<std::size_t>(cells.x);
  field.view.stride_z = static_cast<std::size_t>(cells.x) *
                        static_cast<std::size_t>(cells.y);
  field.view.component_stride = count;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = storage;
  field.view.revision_domain = 81001U;
  return field;
}

CartesianMeshSpec mesh_spec(bool stretched) {
  constexpr std::int32_t n = 16;
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
  mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(n) * n * n;
  mesh.limits.max_memory_bytes_per_rank = 1U << 30U;
  return mesh;
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec face;
  face.flow_kind = BoundaryKind::no_slip_wall;
  face.thermal_kind = BoundaryKind::adiabatic_wall;
  face.mach_limit = 0.95;
  return face;
}

ValidatedModel model_for(CartesianMeshSpec mesh,
                         PressureReferenceKind pressure_reference) {
  ValidatedModel model;
  model.mesh = std::move(mesh);
  model.fingerprint = pressure_reference == PressureReferenceKind::closed_mass
                          ? 0x14a0c501U
                          : 0x14a0c502U;
  model.pressure_reference = pressure_reference;
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  model.schemes.diffusion = DiffusionScheme::central2;
  if (pressure_reference == PressureReferenceKind::closed_mass) {
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::periodic;
      face.thermal_kind = BoundaryKind::none;
      face.mach_limit = 0.95;
    }
    return model;
  }

  for (BoundaryFaceSpec& face : model.boundaries) {
    face = wall();
  }
  BoundaryFaceSpec& inlet = model.boundaries[0U];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.thermal_kind = BoundaryKind::none;
  inlet.velocity = {1.0, 0.0, 0.0};
  inlet.temperature = 300.0;
  BoundaryFaceSpec& outlet = model.boundaries[1U];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.thermal_kind = BoundaryKind::none;
  outlet.pressure = 101325.0;
  return model;
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
  ThermophysicalSpec spec;
  spec.data_file = "analytic.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  spec.species.push_back(air());
  return spec;
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

bool make_fixture(bool stretched, PressureReferenceKind pressure_reference,
                  Fixture& out) {
  const CartesianMeshSpec mesh = mesh_spec(stretched);
  const ValidatedModel model = model_for(mesh, pressure_reference);
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
      !registry.require_field("T", 1U, 2U, temperature) ||
      temperature != 4U ||
      !CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {},
                                          out.geometry, out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry,
                                 out.patch, registry, out.boundary,
                                 out.schemes, out.time)) {
    return false;
  }

  const ThermophysicalSpec thermophysics = thermo_spec();
  if (!ThermodynamicsPlan::compile(thermophysics, {}, out.thermodynamics) ||
      !TransportPlan::compile(thermophysics, out.thermodynamics,
                              out.transport)) {
    return false;
  }
  const std::array<FieldId, 8U> declared{0U, 1U, 2U, 3U,
                                         4U, 5U, 6U, 7U};
  if (!out.contributions.configure({declared.data(), declared.size()}) ||
      !out.contributions.freeze()) {
    return false;
  }

  EquationPlanSpec spec;
  spec.density = density;
  spec.velocity = velocity;
  spec.pressure_perturbation = pressure;
  spec.enthalpy = enthalpy;
  spec.temperature = temperature;
  spec.effective_viscosity = 5U;
  spec.pressure_compressibility = kCompressibility;
  spec.velocity_gradient = 7U;
  spec.pressure_reference = pressure_reference;
  spec.closed_mass_service_stage =
      pressure_reference == PressureReferenceKind::closed_mass
          ? kClosedMassServiceStage
          : 0U;
  spec.maximum_cells_per_rank =
      static_cast<std::size_t>(out.patch.cells.x) *
      static_cast<std::size_t>(out.patch.cells.y) *
      static_cast<std::size_t>(out.patch.cells.z);
  return static_cast<bool>(EquationPlanSet::compile(
      MPI_COMM_SELF, out.schemes, out.geometry, out.patch, out.boundary,
      out.contributions, out.thermodynamics, out.transport, spec,
      out.equations));
}

EquationAssemblyContext context_for(const Fixture& fixture,
                                    StageId contribution_stage) {
  EquationAssemblyContext context;
  context.dt = 0.25;
  context.bdf = {2.75, -3.0, 0.25, 2U};
  context.time = 4101U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = 4102U;
  context.contribution_stage = contribution_stage;
  context.scope = EquationAssemblyScope::final_conservative;
  context.box = {{0, 0, 0}, fixture.patch.cells};
  return context;
}

double volume(const Fixture& fixture, Int3 cell) {
  const std::size_t x =
      static_cast<std::size_t>(fixture.patch.begin.x + cell.x);
  const std::size_t y =
      static_cast<std::size_t>(fixture.patch.begin.y + cell.y);
  const std::size_t z =
      static_cast<std::size_t>(fixture.patch.begin.z + cell.z);
  return fixture.geometry.x().widths().data[x] *
         fixture.geometry.y().widths().data[y] *
         fixture.geometry.z().widths().data[z];
}

EquationAssemblyCertificate marker_certificate() {
  return {991U, EquationAssemblyScope::momentum_predictor, 992U, 993U,
          994U, 995U};
}

PressureReferenceCertificate reference_certificate(
    const Fixture& fixture, const EquationAssemblyContext& context) {
  const PressureReferencePlan& plan = fixture.equations.pressure_reference();
  return {plan.fingerprint(),
          fixture.equations.thermophysical_predictor().fingerprint(),
          fixture.thermodynamics.fingerprint(),
          9981U,
          context.time,
          9982U,
          plan.kind()};
}

bool same_certificate(const EquationAssemblyCertificate& left,
                      const EquationAssemblyCertificate& right) {
  return left.plan == right.plan && left.scope == right.scope &&
         left.time == right.time && left.geometry == right.geometry &&
         left.face_flux == right.face_flux && left.state == right.state;
}

bool all_equal(const std::vector<double>& values, double expected) {
  return std::all_of(values.begin(), values.end(), [expected](double value) {
    return value == expected;
  });
}

void fill_compressibility(OwnedField& psi) {
  const Int3 cells = psi.view.interior;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        psi.view.unchecked({x, y, z}, 0U) =
            1.0e-5 * (1.0 + 0.01 * x + 0.02 * y + 0.03 * z);
      }
    }
  }
}

bool test_oracle(bool stretched) {
  Fixture fixture;
  bool passed = expect(make_fixture(stretched,
                                    PressureReferenceKind::closed_mass,
                                    fixture),
                       stretched ? "stretched closed-mass fixture compiles"
                                 : "uniform closed-mass fixture compiles");
  if (!passed) {
    return false;
  }
  const PressureReferencePlan& plan = fixture.equations.pressure_reference();
  passed &= expect(plan.kind() == PressureReferenceKind::closed_mass &&
                       plan.gauge() ==
                           PressureGauge::compressibility_weighted_zero_mean &&
                       plan.service_stage() == kClosedMassServiceStage &&
                       plan.compressibility_field() == kCompressibility,
                   "closed-mass pressure metadata is compiled");

  EquationAssemblyContext context =
      context_for(fixture, kClosedMassServiceStage);
  OwnedField psi = make_field(kCompressibility, fixture.patch.cells,
                              context.thermo, 82001U);
  OwnedField diagonal =
      make_field(20U, fixture.patch.cells, 81002U, 82002U);
  fill_compressibility(psi);
  std::fill(diagonal.values.begin(), diagonal.values.end(), kSentinel);
  EquationAssemblyCertificate certificate = marker_certificate();
  const PressureReferenceCertificate reference =
      reference_certificate(fixture, context);
  const Status status = assemble_pressure_storage(
      plan, as_const(psi.view), reference, context, diagonal.view,
      certificate);
  passed &= expect(static_cast<bool>(status),
                   stretched ? "stretched pressure storage assembles"
                             : "uniform pressure storage assembles");

  const Int3 cells = fixture.patch.cells;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double expected = context.bdf.a0 *
                                psi.view.unchecked(cell, 0U) *
                                volume(fixture, cell);
        passed &= close(diagonal.view.unchecked(cell, 0U), expected);
      }
    }
  }
  passed &= expect(certificate.plan == plan.fingerprint() &&
                       certificate.scope == context.scope &&
                       certificate.time == context.time &&
                       certificate.geometry == context.geometry &&
                       certificate.face_flux == context.face_flux &&
                       certificate.state != 0U &&
                       certificate.valid(),
                   "pressure storage publishes its dependency certificate");
  PressureReferenceCertificate revised_reference = reference;
  ++revised_reference.pressure_reference;
  EquationAssemblyCertificate revised_certificate;
  passed &= expect(assemble_pressure_storage(
                       plan, as_const(psi.view), revised_reference, context,
                       diagonal.view, revised_certificate) &&
                       revised_certificate.valid() &&
                       revised_certificate.state != certificate.state,
                   "p_ref revision invalidates pressure-storage identity");
  return passed;
}

bool test_failure_atomicity() {
  Fixture fixture;
  bool passed = expect(make_fixture(false, PressureReferenceKind::closed_mass,
                                    fixture),
                       "failure-atomicity fixture compiles");
  if (!passed) {
    return false;
  }
  const PressureReferencePlan& plan = fixture.equations.pressure_reference();
  EquationAssemblyContext context =
      context_for(fixture, kClosedMassServiceStage);
  OwnedField psi = make_field(kCompressibility, fixture.patch.cells,
                              context.thermo, 83001U);
  OwnedField diagonal =
      make_field(21U, fixture.patch.cells, 83002U, 83002U);
  fill_compressibility(psi);
  const EquationAssemblyCertificate marker = marker_certificate();
  const PressureReferenceCertificate reference =
      reference_certificate(fixture, context);

  auto failure_is_atomic = [&](ConstFieldView input,
                               PressureReferenceCertificate authority,
                               const EquationAssemblyContext& bad_context,
                               std::string_view description) {
    std::fill(diagonal.values.begin(), diagonal.values.end(), kSentinel);
    EquationAssemblyCertificate certificate = marker;
    const Status status = assemble_pressure_storage(
        plan, input, authority, bad_context, diagonal.view, certificate);
    return expect(!status && all_equal(diagonal.values, kSentinel) &&
                      same_certificate(certificate, marker),
                  description);
  };

  EquationAssemblyContext bad_stage = context;
  bad_stage.contribution_stage = kClosedMassServiceStage + 1U;
  passed &= failure_is_atomic(as_const(psi.view), reference, bad_stage,
                              "closed-mass wrong service stage is atomic");

  PressureReferenceCertificate wrong_authority = reference;
  ++wrong_authority.plan;
  passed &= failure_is_atomic(as_const(psi.view), wrong_authority, context,
                              "foreign p_ref authority is atomic");

  ConstFieldView wrong_field = as_const(psi.view);
  wrong_field.field = 5U;
  passed &= failure_is_atomic(wrong_field, reference, context,
                              "wrong compressibility field is atomic");
  ConstFieldView stale = as_const(psi.view);
  stale.revision ^= 1U;
  if (stale.revision == 0U) {
    stale.revision = 1U;
  }
  passed &= failure_is_atomic(stale, reference, context,
                              "stale compressibility field is atomic");

  const std::size_t bad_cell = psi.values.size() / 2U;
  const double saved = psi.values[bad_cell];
  psi.values[bad_cell] = std::numeric_limits<double>::quiet_NaN();
  passed &= failure_is_atomic(as_const(psi.view), reference, context,
                              "non-finite compressibility is atomic");
  psi.values[bad_cell] = -saved;
  passed &= failure_is_atomic(as_const(psi.view), reference, context,
                              "non-positive compressibility is atomic");
  psi.values[bad_cell] = saved;

  const std::vector<double> input_snapshot = psi.values;
  EquationAssemblyCertificate alias_certificate = marker;
  const Status alias_status = assemble_pressure_storage(
      plan, as_const(psi.view), reference, context, psi.view,
      alias_certificate);
  passed &= expect(!alias_status && psi.values == input_snapshot &&
                       same_certificate(alias_certificate, marker),
                   "aliased pressure output is rejected atomically");
  return passed;
}

bool test_open_pressure_reference() {
  Fixture fixture;
  bool passed = expect(
      make_fixture(false, PressureReferenceKind::boundary_absolute, fixture),
      "open pressure-reference fixture compiles");
  if (!passed) {
    return false;
  }
  const PressureReferencePlan& plan = fixture.equations.pressure_reference();
  passed &= expect(plan.kind() == PressureReferenceKind::boundary_absolute &&
                       plan.gauge() ==
                           PressureGauge::absolute_boundary_dirichlet &&
                       plan.service_stage() == 0U &&
                       plan.compressibility_field() == kCompressibility,
                   "open pressure-reference metadata has no mass-service stage");

  // The source/contribution stage remains meaningful to other equations, but
  // an absolute-boundary pressure correction must not acquire the closed-mass
  // service-stage dependency.
  EquationAssemblyContext context = context_for(fixture, 911U);
  OwnedField psi = make_field(kCompressibility, fixture.patch.cells,
                              context.thermo, 84001U);
  OwnedField diagonal =
      make_field(22U, fixture.patch.cells, 84002U, 84002U);
  fill_compressibility(psi);
  EquationAssemblyCertificate certificate = marker_certificate();
  const PressureReferenceCertificate reference =
      reference_certificate(fixture, context);
  passed &= expect(assemble_pressure_storage(plan, as_const(psi.view),
                                             reference, context, diagonal.view,
                                             certificate) &&
                       certificate.valid(),
                   "open pressure storage is independent of mass-service stage");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_oracle(false);
  passed &= test_oracle(true);
  passed &= test_failure_atomicity();
  passed &= test_open_pressure_reference();
  const int finalize_status = MPI_Finalize();
  return passed && finalize_status == MPI_SUCCESS ? 0 : 1;
}
