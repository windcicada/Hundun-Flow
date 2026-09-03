// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_flow.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr StageId kFirstStage = 31U;
constexpr StageId kSecondStage = 47U;
constexpr FieldId kEnthalpy = 3U;
constexpr FieldId kFirstSource = 8U;
constexpr FieldId kSecondSource = 9U;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             1.0e-12 * std::max(1.0, std::abs(expected));
}

struct OwnedField {
  std::vector<double> values;
  FieldView view{};
};

struct OwnedFaceField {
  std::vector<double> values;
  FaceFieldView view{};
};

OwnedField make_field(FieldId id, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision) {
  OwnedField field;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  field.values.assign(nx * ny * nz * components, 0.0);
  field.view.base =
      field.values.data() + ghosts + ghosts * nx + ghosts * nx * ny;
  field.view.interior = cells;
  field.view.ghosts = {ghosts, ghosts, ghosts};
  field.view.components = components;
  field.view.stride_y = nx;
  field.view.stride_z = nx * ny;
  field.view.component_stride = nx * ny * nz;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = static_cast<StorageIdentity>(1000U + id);
  field.view.revision_domain = 15001U;
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
  field.values.assign(static_cast<std::size_t>(extents.x) * extents.y *
                          extents.z,
                      0.0);
  field.view = {field.values.data(),
                extents,
                static_cast<std::size_t>(extents.x),
                static_cast<std::size_t>(extents.x) * extents.y,
                axis,
                identity,
                15002U};
  return field;
}

CartesianMeshSpec mesh_spec(std::int32_t cells_per_axis) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {cells_per_axis, cells_per_axis, cells_per_axis};
  mesh.minimum_spacing = {1.0 / cells_per_axis, 1.0 / cells_per_axis,
                          1.0 / cells_per_axis};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells =
      static_cast<std::uint64_t>(cells_per_axis) * cells_per_axis *
      cells_per_axis;
  mesh.limits.max_memory_bytes_per_rank = 1U << 26U;
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

UnitDimension source_units() {
  UnitDimension units;
  units.si_exponents = {1, -1, -3, 0, 0, 0, 0};
  return units;
}

UnitDimension wrong_source_units() {
  UnitDimension units;
  units.si_exponents = {1, 0, -1, 0, 0, 0, 0};
  return units;
}

ContributionSpec contribution(StageId stage, FieldId source,
                              Span<const FieldId> reads) {
  ContributionSpec spec;
  spec.conserved_quantity = kEnthalpy;
  spec.units = source_units();
  spec.stage = stage;
  spec.reads = reads;
  spec.explicit_source = source;
  spec.capability = ContributionCapability::inert_source;
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

bool make_fixture(std::int32_t cells_per_axis, Fixture& out) {
  const CartesianMeshSpec mesh = mesh_spec(cells_per_axis);
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14a057a9U;
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
      !registry.require_field("h", 1U, 2U, enthalpy) || enthalpy != kEnthalpy ||
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

  const std::array<FieldId, 10U> declared{0U, 1U, 2U, 3U, 4U,
                                          5U, 6U, 7U, 8U, 9U};
  const std::array<FieldId, 1U> reads{0U};
  if (!out.contributions.configure({declared.data(), declared.size()}) ||
      !out.contributions.register_contribution(contribution(
          kSecondStage, kSecondSource, {reads.data(), reads.size()})) ||
      !out.contributions.register_contribution(contribution(
          kFirstStage, kFirstSource, {reads.data(), reads.size()})) ||
      !out.contributions.freeze()) {
    return false;
  }

  EquationPlanSpec spec;
  spec.density = 0U;
  spec.velocity = 1U;
  spec.pressure_perturbation = 2U;
  spec.enthalpy = kEnthalpy;
  spec.temperature = 4U;
  spec.effective_viscosity = 5U;
  spec.pressure_compressibility = 6U;
  spec.velocity_gradient = 7U;
  spec.pressure_reference = PressureReferenceKind::closed_mass;
  spec.closed_mass_service_stage = 1U;
  spec.maximum_cells_per_rank =
      static_cast<std::size_t>(cells_per_axis) * cells_per_axis *
      cells_per_axis;
  return static_cast<bool>(EquationPlanSet::compile(
      MPI_COMM_SELF, out.schemes, out.geometry, out.patch, out.boundary,
      out.contributions, out.thermodynamics, out.transport, spec,
      out.equations));
}

bool test_wrong_equation_source_units_reject_at_compile() {
  constexpr std::int32_t kCellsPerAxis = 4;
  const CartesianMeshSpec mesh = mesh_spec(kCellsPerAxis);
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14a057b0U;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  FieldRegistry fields;
  FieldId id = 0U;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  const ThermophysicalSpec thermo = thermo_spec();
  const std::array<FieldId, 10U> declared{0U, 1U, 2U, 3U, 4U,
                                          5U, 6U, 7U, 8U, 9U};
  const std::array<FieldId, 1U> reads{0U};
  ContributionSpec wrong =
      contribution(kFirstStage, kFirstSource,
                   {reads.data(), reads.size()});
  wrong.units = wrong_source_units();
  bool ready = fields.require_field("rho", 1U, 2U, id) && id == 0U &&
               fields.require_field("U", 3U, 2U, id) && id == 1U &&
               fields.require_field("pi", 1U, 2U, id) && id == 2U &&
               fields.require_field("h", 1U, 2U, id) && id == 3U &&
               fields.require_field("T", 1U, 2U, id) && id == 4U &&
               CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {},
                                                   geometry, patch) &&
               BoundaryCompiler::compile(MPI_COMM_SELF, model, geometry,
                                         patch, fields, boundary, schemes,
                                         time) &&
               ThermodynamicsPlan::compile(thermo, {}, thermodynamics) &&
               TransportPlan::compile(thermo, thermodynamics, transport) &&
               contributions.configure({declared.data(), declared.size()}) &&
               contributions.register_contribution(wrong) &&
               contributions.freeze();
  if (!ready) {
    return expect(false, "wrong-unit compile fixture prepares");
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
  spec.maximum_cells_per_rank = 64U;
  EquationPlanSet equations;
  return expect(EquationPlanSet::compile(
                    MPI_COMM_SELF, schemes, geometry, patch, boundary,
                    contributions, thermodynamics, transport, spec,
                    equations).code == StatusCode::invalid_plan,
                "wrong physical source units reject during equation compile");
}

struct FinalFluxFixture {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxWriter writer;
};

bool make_zero_final_flux(const CartesianKernelPlan& kernels, Int3 cells,
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
  if (!authority.claim(61U, 0U, fixture.transaction, fixture.writer) ||
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
  OwnedField density = make_field(80U, cells, 1U, 2U, 501U);
  OwnedField velocity = make_field(81U, cells, 3U, 2U, 502U);
  std::fill(density.values.begin(), density.values.end(), 1.0);
  const std::array<ConstFieldView, 2U> reads{as_const(density.view),
                                             as_const(velocity.view)};
  const KernelInvocation call{{reads.data(), reads.size()},
                              {},
                              {{0, 0, 0}, cells},
                              0U,
                              0U,
                              1U,
                              0U,
                              nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(kernels, call, pending)) &&
         static_cast<bool>(fixture.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(fixture.transaction.collective_finish(MPI_COMM_SELF,
                                                                  Status{})) &&
         static_cast<bool>(fixture.writer.committed(fixture.storage,
                                                     committed));
}

EquationContributionView contribution_view(StageId stage, FieldId source,
                                           ConstFieldView values) {
  EquationContributionView view;
  view.explicit_source_density = values;
  view.has_implicit_sink = false;
  view.conserved_quantity = kEnthalpy;
  view.units = source_units();
  view.stage = stage;
  view.explicit_source_field = source;
  return view;
}

struct OutputSnapshot {
  std::vector<double> diagonal;
  std::vector<double> rhs;
  std::vector<double> residual;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  EquationAssemblyCertificate certificate{};
};

OutputSnapshot snapshot(const OwnedField& diagonal, const OwnedField& rhs,
                        const OwnedField& residual, const OwnedFaceField& x,
                        const OwnedFaceField& y, const OwnedFaceField& z,
                        EquationAssemblyCertificate certificate) {
  return {diagonal.values, rhs.values, residual.values, x.values, y.values,
          z.values, certificate};
}

bool same_certificate(EquationAssemblyCertificate left,
                      EquationAssemblyCertificate right) {
  return left.plan == right.plan && left.scope == right.scope &&
         left.time == right.time && left.geometry == right.geometry &&
         left.face_flux == right.face_flux && left.state == right.state &&
         left.dt == right.dt;
}

bool unchanged(const OutputSnapshot& before, const OwnedField& diagonal,
               const OwnedField& rhs, const OwnedField& residual,
               const OwnedFaceField& x, const OwnedFaceField& y,
               const OwnedFaceField& z,
               EquationAssemblyCertificate certificate) {
  return before.diagonal == diagonal.values && before.rhs == rhs.values &&
         before.residual == residual.values && before.x == x.values &&
         before.y == y.values && before.z == z.values &&
         same_certificate(before.certificate, certificate);
}

void reset_outputs(OwnedField& diagonal, OwnedField& rhs,
                   OwnedField& residual, OwnedFaceField& x,
                   OwnedFaceField& y, OwnedFaceField& z, double sentinel) {
  std::fill(diagonal.values.begin(), diagonal.values.end(), sentinel - 1.0);
  std::fill(rhs.values.begin(), rhs.values.end(), sentinel - 2.0);
  std::fill(residual.values.begin(), residual.values.end(), sentinel - 3.0);
  std::fill(x.values.begin(), x.values.end(), sentinel - 4.0);
  std::fill(y.values.begin(), y.values.end(), sentinel - 5.0);
  std::fill(z.values.begin(), z.values.end(), sentinel - 6.0);
}

bool test_enthalpy_sources_are_selected_by_stage() {
  constexpr std::int32_t kCellsPerAxis = 4;
  Fixture fixture;
  bool passed = expect(make_fixture(kCellsPerAxis, fixture),
                       "two-stage enthalpy equation plan compiles");
  if (!passed) {
    return false;
  }

  const Int3 cells = fixture.patch.cells;
  OwnedField density = make_field(0U, cells, 1U, 0U, 601U);
  OwnedField velocity = make_field(1U, cells, 3U, 1U, 602U);
  OwnedField pressure = make_field(2U, cells, 1U, 1U, 603U);
  OwnedField enthalpy = make_field(kEnthalpy, cells, 1U, 2U, 604U);
  OwnedField temperature = make_field(4U, cells, 1U, 1U, 605U);
  OwnedField viscosity = make_field(5U, cells, 1U, 0U, 606U);
  OwnedField conductivity = make_field(10U, cells, 1U, 1U, 607U);
  OwnedField enthalpy_diffusivity =
      make_field(11U, cells, 1U, 1U, 608U);
  OwnedField velocity_gradient = make_field(7U, cells, 9U, 0U, 609U);
  OwnedField first_source = make_field(kFirstSource, cells, 1U, 0U, 610U);
  OwnedField second_source = make_field(kSecondSource, cells, 1U, 0U, 611U);
  std::fill(density.values.begin(), density.values.end(), 1.0);
  std::fill(enthalpy.values.begin(), enthalpy.values.end(), 3.0e5);
  std::fill(temperature.values.begin(), temperature.values.end(), 300.0);
  std::fill(viscosity.values.begin(), viscosity.values.end(), 1.0);
  std::fill(conductivity.values.begin(), conductivity.values.end(), 1.0);
  std::fill(enthalpy_diffusivity.values.begin(),
            enthalpy_diffusivity.values.end(), 1.0);
  std::fill(first_source.values.begin(), first_source.values.end(), 2.0);
  std::fill(second_source.values.begin(), second_source.values.end(), 5.0);

  FinalFluxFixture final_flux;
  ConstFaceFluxView committed;
  passed &= expect(make_zero_final_flux(fixture.equations.kernels(), cells,
                                        final_flux, committed),
                   "zero final flux is certified");
  if (!passed) {
    return false;
  }

  EquationStateView state;
  state.density = {as_const(density.view), as_const(density.view),
                   as_const(density.view)};
  state.velocity = {as_const(velocity.view), as_const(velocity.view),
                    as_const(velocity.view)};
  state.pressure_perturbation = {
      as_const(pressure.view), as_const(pressure.view),
      as_const(pressure.view)};
  state.enthalpy = {as_const(enthalpy.view), as_const(enthalpy.view),
                    as_const(enthalpy.view)};
  state.temperature = {as_const(temperature.view), as_const(temperature.view),
                       as_const(temperature.view)};
  state.pressure_reference = 101325.0;
  state.accepted_pressure_reference = 101325.0;
  state.previous_pressure_reference = 101325.0;

  EquationMaterialView material;
  material.effective_viscosity = as_const(viscosity.view);
  material.thermal_conductivity = as_const(conductivity.view);
  material.enthalpy_diffusivity = as_const(enthalpy_diffusivity.view);

  OwnedField diagonal = make_field(30U, cells, 1U, 0U, 701U);
  OwnedField rhs = make_field(31U, cells, 1U, 0U, 702U);
  OwnedField residual = make_field(32U, cells, 1U, 0U, 703U);
  OwnedFaceField x = make_face_field(CartesianAxis::x, cells, 704U);
  OwnedFaceField y = make_face_field(CartesianAxis::y, cells, 705U);
  OwnedFaceField z = make_face_field(CartesianAxis::z, cells, 706U);
  EquationSystemView system{diagonal.view, rhs.view, residual.view, x.view,
                            y.view, z.view};

  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 801U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = committed.revision;
  context.face_flux_authority = committed.certificate.authority();
  context.face_flux_storage = committed.certificate.storage();
  context.face_flux_revision_domain =
      committed.certificate.revision_domain();
  context.scope = EquationAssemblyScope::final_conservative;
  context.mass_flux = committed;
  context.provisional_mass_flux = false;

  const EquationContributionView first = contribution_view(
      kFirstStage, kFirstSource, as_const(first_source.view));
  const EquationContributionView second = contribution_view(
      kSecondStage, kSecondSource, as_const(second_source.view));
  const double volume = 1.0 / static_cast<double>(kCellsPerAxis *
                                                  kCellsPerAxis *
                                                  kCellsPerAxis);
  const Int3 probe{1, 1, 1};
  EquationAssemblyCertificate certificate;

  context.contribution_stage = kFirstStage;
  const std::array first_only{first};
  passed &= expect(static_cast<bool>(assemble_enthalpy(
                       fixture.equations.enthalpy(), state, material,
                       as_const(velocity_gradient.view),
                       {first_only.data(), first_only.size()}, context, system,
                       certificate)),
                   "first stage accepts only its compiled contribution");
  passed &= expect(close(residual.view.unchecked(probe, 0U), -2.0 * volume),
                   "first stage assembles only the first source");

  reset_outputs(diagonal, rhs, residual, x, y, z, -1000.0);
  certificate = {};
  context.contribution_stage = kSecondStage;
  const std::array second_only{second};
  passed &= expect(static_cast<bool>(assemble_enthalpy(
                       fixture.equations.enthalpy(), state, material,
                       as_const(velocity_gradient.view),
                       {second_only.data(), second_only.size()}, context,
                       system, certificate)),
                   "second stage accepts only its compiled contribution");
  passed &= expect(close(residual.view.unchecked(probe, 0U), -5.0 * volume),
                   "second stage assembles only the second source");

  reset_outputs(diagonal, rhs, residual, x, y, z, -2000.0);
  certificate = {901U, EquationAssemblyScope::momentum_predictor, 902U, 903U,
                 904U, 905U, 906.0};
  const OutputSnapshot before_wrong_stage =
      snapshot(diagonal, rhs, residual, x, y, z, certificate);
  context.contribution_stage = kFirstStage;
  const Status wrong_stage = assemble_enthalpy(
      fixture.equations.enthalpy(), state, material,
      as_const(velocity_gradient.view),
      {second_only.data(), second_only.size()}, context, system, certificate);
  passed &= expect(wrong_stage.code == StatusCode::invalid_plan,
                   "a contribution from the wrong stage is rejected");
  passed &= expect(unchanged(before_wrong_stage, diagonal, rhs, residual, x, y,
                             z, certificate),
                   "wrong-stage rejection leaves every output atomic");

  reset_outputs(diagonal, rhs, residual, x, y, z, -3000.0);
  certificate = {911U, EquationAssemblyScope::momentum_predictor, 912U, 913U,
                 914U, 915U, 916.0};
  const OutputSnapshot before_wrong_descriptor =
      snapshot(diagonal, rhs, residual, x, y, z, certificate);
  EquationContributionView wrong_descriptor = first;
  ++wrong_descriptor.units.si_exponents[0U];
  const std::array descriptor_mismatch{wrong_descriptor};
  const Status wrong_descriptor_status = assemble_enthalpy(
      fixture.equations.enthalpy(), state, material,
      as_const(velocity_gradient.view),
      {descriptor_mismatch.data(), descriptor_mismatch.size()}, context,
      system, certificate);
  passed &= expect(wrong_descriptor_status.code == StatusCode::invalid_plan,
                   "a view that disagrees with its compiled descriptor is rejected");
  passed &= expect(unchanged(before_wrong_descriptor, diagonal, rhs, residual,
                             x, y, z, certificate),
                   "descriptor rejection leaves every output atomic");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  const bool passed = test_enthalpy_sources_are_selected_by_stage();
  const bool units_passed = test_wrong_equation_source_units_reject_at_compile();
  MPI_Finalize();
  return passed && units_passed ? 0 : 1;
}
