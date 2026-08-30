// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_flow.hpp"

#include <mpi.h>

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
  field.view.revision_domain = 7001U;
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
                8001U};
  return field;
}

CartesianMeshSpec mesh_spec(std::int32_t n, bool stretched) {
  CartesianMeshSpec mesh;
  mesh.kind = stretched ? GeometryKind::tensor_stretched : GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {n, n, n};
  if (stretched) {
    const double inverse_n = 1.0 / static_cast<double>(n);
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
    mesh.minimum_spacing = {1.0 / n, 1.0 / n, 1.0 / n};
    mesh.max_growth_ratio = 1.0;
  }
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(n) * n * n;
  mesh.limits.max_memory_bytes_per_rank = 1U << 30U;
  return mesh;
}

ValidatedModel periodic_model(const CartesianMeshSpec& mesh) {
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x14a0c111U;
  model.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  model.schemes.diffusion = DiffusionScheme::central2;
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

bool make_fixture(std::int32_t n, bool stretched, Fixture& out) {
  const CartesianMeshSpec mesh = mesh_spec(n, stretched);
  const ValidatedModel model = periodic_model(mesh);
  FieldRegistry registry;
  FieldId density = 0U;
  if (!registry.require_field("rho", 1U, 2U, density) || density != 0U) {
    return false;
  }
  if (!CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes, out.time)) {
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

double x_center(const CartesianGeometryPlan& geometry, std::int32_t i) {
  return geometry.x().centres().data[static_cast<std::size_t>(i)];
}

double y_center(const CartesianGeometryPlan& geometry, std::int32_t j) {
  return geometry.y().centres().data[static_cast<std::size_t>(j)];
}

double z_center(const CartesianGeometryPlan& geometry, std::int32_t k) {
  return geometry.z().centres().data[static_cast<std::size_t>(k)];
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

double cell_volume(const CartesianGeometryPlan& geometry, Int3 cell) {
  return geometry.x().widths().data[static_cast<std::size_t>(cell.x)] *
         geometry.y().widths().data[static_cast<std::size_t>(cell.y)] *
         geometry.z().widths().data[static_cast<std::size_t>(cell.z)];
}

double rho_value(double x, double y, double z) {
  return 1.1 + 0.08 * std::sin(2.0 * kPi * x) +
         0.06 * std::cos(2.0 * kPi * y) +
         0.04 * std::sin(2.0 * kPi * z);
}

Real3 velocity_value(double x, double y, double z) {
  return {0.35 + 0.12 * std::sin(2.0 * kPi * y),
          -0.18 + 0.09 * std::cos(2.0 * kPi * z),
          0.11 + 0.07 * std::sin(2.0 * kPi * x)};
}

double continuity_oracle(double x, double y, double z) {
  const double w = 2.0 * kPi;
  const Real3 u = velocity_value(x, y, z);
  const double drdx = 0.08 * w * std::cos(w * x);
  const double drdy = -0.06 * w * std::sin(w * y);
  const double drdz = 0.04 * w * std::cos(w * z);
  return u.x * drdx + u.y * drdy + u.z * drdz;
}


double run_momentum_mms(std::int32_t n, bool stretched) {
  Fixture fixture;
  if (!make_fixture(n, stretched, fixture)) {
    return HUGE_VAL;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField rho = make_field(0U, cells, 1U, 2U, 1101U);
  OwnedField velocity = make_field(1U, cells, 3U, 2U, 1102U);
  OwnedField pi = make_field(2U, cells, 1U, 1U, 1103U);
  OwnedField viscosity = make_field(5U, cells, 1U, 1U, 1104U);
  OwnedField gradients = make_field(7U, cells, 9U, 1U, 1105U);
  OwnedField diagonal = make_field(50U, cells, 3U, 0U, 1106U);
  OwnedField rhs = make_field(51U, cells, 3U, 0U, 1107U);
  OwnedField residual = make_field(52U, cells, 3U, 0U, 1108U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 1109U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 1110U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 1111U);

  constexpr double u0 = 0.37;
  constexpr double amplitude = 0.11;
  constexpr double pressure_amplitude = 0.23;
  constexpr double mu = 0.041;
  const double wave = 2.0 * kPi;
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    const double z = periodic_center(fixture.geometry.z(), k, n);
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      const double y = periodic_center(fixture.geometry.y(), j, n);
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const double x = periodic_center(fixture.geometry.x(), i, n);
        const Int3 cell{i, j, k};
        rho.view.unchecked(cell, 0U) = 1.0;
        velocity.view.unchecked(cell, 0U) =
            u0 + amplitude * std::sin(wave * x);
        velocity.view.unchecked(cell, 1U) = 0.0;
        velocity.view.unchecked(cell, 2U) = 0.0;
        if (i >= -1 && i < cells.x + 1 && j >= -1 &&
            j < cells.y + 1 && k >= -1 && k < cells.z + 1) {
          pi.view.unchecked(cell, 0U) =
              pressure_amplitude * std::cos(wave * x);
          viscosity.view.unchecked(cell, 0U) = mu;
          gradients.view.unchecked(cell, 0U) =
              amplitude * wave * std::cos(wave * x);
        }
        (void)y;
        (void)z;
      }
    }
  }

  FaceFluxStorage flux_storage;
  FaceFluxView flux;
  if (!FaceFluxStorage::allocate_workspace(cells, 1U, flux_storage) ||
      !flux_storage.workspace_view(0U, 1112U, flux)) {
    return HUGE_VAL;
  }
  const std::array<ConstFieldView, 2U> mass_reads{as_const(rho.view),
                                                  as_const(velocity.view)};
  const KernelInvocation mass_call{{mass_reads.data(), mass_reads.size()}, {},
                                   {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U,
                                   nullptr};
  if (!reconstruct_mass_flux(fixture.equations.kernels(), mass_call, flux)) {
    return HUGE_VAL;
  }

  EquationStateView state;
  state.density = {as_const(rho.view), as_const(rho.view), as_const(rho.view)};
  state.velocity = {as_const(velocity.view), as_const(velocity.view),
                    as_const(velocity.view)};
  state.pressure_perturbation = {as_const(pi.view), as_const(pi.view),
                                 as_const(pi.view)};
  EquationMaterialView material;
  material.effective_viscosity = as_const(viscosity.view);
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 1113U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = 1112U;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::momentum_predictor;
  context.mass_flux = as_const(flux);
  context.provisional_mass_flux = true;
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};
  EquationAssemblyCertificate certificate;
  if (!assemble_momentum(fixture.equations.momentum(), state, material,
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
        const double x = x_center(fixture.geometry, i);
        const double s = std::sin(wave * x);
        const double c = std::cos(wave * x);
        const double u = u0 + amplitude * s;
        const double du = amplitude * wave * c;
        const double d2u = -amplitude * wave * wave * s;
        const double dpdx = -pressure_amplitude * wave * s;
        // The production equation is conservative.  Because this manufactured
        // velocity has non-zero divergence, div(U tensor U)=2*u*du/dx rather
        // than the non-conservative U dot grad(U)=u*du/dx.
        const double expected_density =
            2.0 * u * du + dpdx - (4.0 / 3.0) * mu * d2u;
        const double expected = expected_density *
                                cell_volume(fixture.geometry, cell);
        const double difference =
            residual.view.unchecked(cell, 0U) - expected;
        error += difference * difference;
        measure += expected * expected;
        error += residual.view.unchecked(cell, 1U) *
                 residual.view.unchecked(cell, 1U);
        error += residual.view.unchecked(cell, 2U) *
                 residual.view.unchecked(cell, 2U);
      }
    }
  }
  return std::sqrt(static_cast<double>(error / measure));
}

double run_continuity(std::int32_t n, bool stretched) {
  Fixture fixture;
  if (!make_fixture(n, stretched, fixture)) {
    return HUGE_VAL;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField rho = make_field(0U, cells, 1U, 2U, 11U);
  OwnedField velocity = make_field(1U, cells, 3U, 2U, 12U);
  OwnedField pi = make_field(2U, cells, 1U, 2U, 13U);
  OwnedField h = make_field(3U, cells, 1U, 2U, 14U);
  OwnedField temperature = make_field(4U, cells, 1U, 2U, 15U);
  OwnedField residual = make_field(20U, cells, 1U, 0U, 16U);
  FaceFluxStorage flux_storage;
  FaceFluxView flux;
  if (!FaceFluxStorage::allocate_workspace(cells, 1U, flux_storage) ||
      !flux_storage.workspace_view(0U, 101U, flux)) {
    return HUGE_VAL;
  }
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    const std::int32_t gk = (k % n + n) % n;
    const double z = z_center(fixture.geometry, gk);
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      const std::int32_t gj = (j % n + n) % n;
      const double y = y_center(fixture.geometry, gj);
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const std::int32_t gi = (i % n + n) % n;
        const double x = x_center(fixture.geometry, gi);
        const Int3 cell{i, j, k};
        rho.view.unchecked(cell, 0U) = rho_value(x, y, z);
        const Real3 u = velocity_value(x, y, z);
        velocity.view.unchecked(cell, 0U) = u.x;
        velocity.view.unchecked(cell, 1U) = u.y;
        velocity.view.unchecked(cell, 2U) = u.z;
        pi.view.unchecked(cell, 0U) = 300.0 * std::sin(2.0 * kPi * x);
        h.view.unchecked(cell, 0U) = 3.0e5 + 1000.0 * y;
        temperature.view.unchecked(cell, 0U) = 300.0 + 10.0 * y;
      }
    }
  }
  const std::array<ConstFieldView, 2U> mass_reads{as_const(rho.view),
                                                  as_const(velocity.view)};
  const KernelInvocation mass_call{{mass_reads.data(), mass_reads.size()}, {},
                                   {{0, 0, 0}, cells}, 0U, 0U, 1U, 0U,
                                   nullptr};
  if (!reconstruct_mass_flux(fixture.equations.kernels(), mass_call, flux)) {
    return HUGE_VAL;
  }
  const PrimitiveHistory density{as_const(rho.view), as_const(rho.view),
                                 as_const(rho.view)};
  const PrimitiveHistory vector_velocity{as_const(velocity.view),
                                         as_const(velocity.view),
                                         as_const(velocity.view)};
  const PrimitiveHistory pressure{as_const(pi.view), as_const(pi.view),
                                  as_const(pi.view)};
  const PrimitiveHistory enthalpy{as_const(h.view), as_const(h.view),
                                  as_const(h.view)};
  const PrimitiveHistory thermal{as_const(temperature.view),
                                 as_const(temperature.view),
                                 as_const(temperature.view)};
  EquationStateView state{density, vector_velocity, pressure, enthalpy, thermal,
                          101325.0, {}, {}};
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 100U;
  context.geometry = fixture.geometry.topology_revision();
  context.face_flux = 101U;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::momentum_predictor;
  context.mass_flux = as_const(flux);
  context.provisional_mass_flux = true;
  EquationSystemView system;
  system.residual = residual.view;
  EquationAssemblyCertificate certificate;
  if (!assemble_continuity(fixture.equations.continuity(), state, context,
                           system, certificate)) {
    return HUGE_VAL;
  }
  long double error = 0.0L;
  long double measure = 0.0L;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    const double z = z_center(fixture.geometry, k);
    for (std::int32_t j = 0; j < cells.y; ++j) {
      const double y = y_center(fixture.geometry, j);
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double x = x_center(fixture.geometry, i);
        const Int3 cell{i, j, k};
        const double volume = cell_volume(fixture.geometry, cell);
        const double expected = continuity_oracle(x, y, z) * volume;
        const double difference = residual.view.unchecked(cell, 0U) - expected;
        error += difference * difference;
        measure += expected * expected;
      }
    }
  }
  return std::sqrt(static_cast<double>(error / measure));
}

double order(double coarse, double fine) {
  return std::log(coarse / fine) / std::log(2.0);
}

bool test_mms_orders() {
  bool passed = true;
  for (const bool stretched : {false, true}) {
    const double e16 = run_continuity(16, stretched);
    const double e32 = run_continuity(32, stretched);
    const double e64 = run_continuity(64, stretched);
    passed &= expect(std::isfinite(e16) && e16 > e32 && e32 > e64,
                     "variable-density cell-integral residual converges");
    passed &= expect(order(e16, e32) >= 1.8 && order(e32, e64) >= 1.8,
                     "uniform/stretched continuity has two adjacent orders >= 1.8");
    const double m16 = run_momentum_mms(16, stretched);
    const double m32 = run_momentum_mms(32, stretched);
    const double m64 = run_momentum_mms(64, stretched);
    passed &= expect(std::isfinite(m16) && m16 > m32 && m32 > m64,
                     "production full-stress momentum MMS converges");
    passed &= expect(order(m16, m32) >= 1.8 &&
                         order(m32, m64) >= 1.8,
                     "uniform/stretched momentum has two adjacent orders >= 1.8");
  }
  return passed;
}

bool test_full_newtonian_stress_contract() {
  static_assert(static_cast<int>(MomentumStressForm::full_newtonian) !=
                    static_cast<int>(MomentumStressForm::laplacian_only),
                "full Newtonian stress must be a distinct selectable contract");
  EquationAssemblyScope predictor = EquationAssemblyScope::momentum_predictor;
  EquationAssemblyScope final = EquationAssemblyScope::final_conservative;
  return expect(predictor != final,
                "predictor and final conservative authorities stay distinct");
}

bool test_momentum_full_stress_and_operator_oracle() {
  Fixture fixture;
  constexpr std::int32_t n = 8;
  bool passed = expect(make_fixture(n, false, fixture),
                       "momentum oracle fixture compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField rho = make_field(0U, cells, 1U, 2U, 201U);
  OwnedField velocity = make_field(1U, cells, 3U, 2U, 202U);
  OwnedField pi = make_field(2U, cells, 1U, 2U, 203U);
  OwnedField h = make_field(3U, cells, 1U, 2U, 204U);
  OwnedField temperature = make_field(4U, cells, 1U, 2U, 205U);
  OwnedField viscosity = make_field(5U, cells, 1U, 2U, 206U);
  OwnedField diagonal = make_field(30U, cells, 3U, 0U, 207U);
  OwnedField rhs = make_field(31U, cells, 3U, 0U, 208U);
  OwnedField residual = make_field(32U, cells, 3U, 0U, 209U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 210U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 211U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 212U);

  constexpr double mu = 0.75;
  const double spacing = 1.0 / static_cast<double>(n);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const Int3 cell{i, j, k};
        const double x = (static_cast<double>(i) + 0.5) * spacing;
        rho.view.unchecked(cell, 0U) = 1.0;
        velocity.view.unchecked(cell, 0U) = x * x;
        velocity.view.unchecked(cell, 1U) = 0.0;
        velocity.view.unchecked(cell, 2U) = 0.0;
        pi.view.unchecked(cell, 0U) = 0.0;
        h.view.unchecked(cell, 0U) = 3.0e5;
        temperature.view.unchecked(cell, 0U) = 300.0;
        viscosity.view.unchecked(cell, 0U) = mu;
      }
    }
  }
  std::fill(diagonal.bytes.begin(), diagonal.bytes.end(), -901.0);
  std::fill(rhs.bytes.begin(), rhs.bytes.end(), -902.0);
  std::fill(residual.bytes.begin(), residual.bytes.end(), -903.0);

  FaceFluxStorage flux_storage;
  FaceFluxView flux;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, flux_storage)) &&
                       static_cast<bool>(flux_storage.workspace_view(
                           0U, 301U, flux)),
                   "zero provisional mass flux allocates");
  if (!passed) {
    return false;
  }
  // FaceFluxStorage is zero-initialized; this isolates unsteady, pressure and
  // the complete Newtonian stress without relying on a convection oracle.
  const PrimitiveHistory density{as_const(rho.view), as_const(rho.view),
                                 as_const(rho.view)};
  const PrimitiveHistory vector_velocity{as_const(velocity.view),
                                         as_const(velocity.view),
                                         as_const(velocity.view)};
  const PrimitiveHistory pressure{as_const(pi.view), as_const(pi.view),
                                  as_const(pi.view)};
  const PrimitiveHistory enthalpy{as_const(h.view), as_const(h.view),
                                  as_const(h.view)};
  const PrimitiveHistory thermal{as_const(temperature.view),
                                 as_const(temperature.view),
                                 as_const(temperature.view)};
  EquationStateView state{density, vector_velocity, pressure, enthalpy, thermal,
                          101325.0, {}, {}, 101325.0, 101325.0};
  EquationMaterialView material;
  material.effective_viscosity = as_const(viscosity.view);
  OwnedField gradients = make_field(7U, cells, 9U, 1U, 213U);
  for (std::int32_t k = -1; k < cells.z + 1; ++k) {
    for (std::int32_t j = -1; j < cells.y + 1; ++j) {
      for (std::int32_t i = -1; i < cells.x + 1; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * spacing;
        gradients.view.unchecked({i, j, k}, 0U) = 2.0 * x;
      }
    }
  }
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {10.0, -10.0, 0.0, 1U};
  context.time = 300U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = 301U;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::momentum_predictor;
  context.mass_flux = as_const(flux);
  context.provisional_mass_flux = true;
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};
  EquationAssemblyCertificate certificate;
  const Status momentum_status = assemble_momentum(
      fixture.equations.momentum(), state, material, as_const(gradients.view),
      {}, context, system, certificate);
  passed &= expect(static_cast<bool>(momentum_status),
                   "production momentum assembler succeeds");

  const Int3 cell{3, 3, 3};
  const double volume = spacing * spacing * spacing;
  const double transmissibility = mu * spacing;
  const double expected_diagonal = 10.0 * volume + 6.0 * transmissibility;
  const double expected_full_stress = -(8.0 / 3.0) * mu * volume;
  const double laplacian_only_mutation = -2.0 * mu * volume;
  passed &= expect(std::abs(residual.view.unchecked(cell, 0U) -
                            expected_full_stress) < 2.0e-14,
                   "momentum residual uses div of full Newtonian stress");
  passed &= expect(std::abs(residual.view.unchecked(cell, 0U) -
                            laplacian_only_mutation) > 1.0e-5,
                   "momentum oracle kills the Laplacian-only mutation");
  passed &= expect(std::abs(residual.view.unchecked(cell, 1U)) < 2.0e-14 &&
                       std::abs(residual.view.unchecked(cell, 2U)) < 2.0e-14,
                   "full stress does not leak into transverse momentum");
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    passed &= expect(std::abs(diagonal.view.unchecked(cell, component) -
                              expected_diagonal) < 2.0e-14,
                     "momentum diagonal is cell-integral BDF plus six faces");
  }
  passed &= expect(std::abs(ax.view.unchecked({4, 3, 3}) -
                            transmissibility) < 2.0e-14 &&
                       std::abs(ay.view.unchecked({3, 4, 3}) -
                                transmissibility) < 2.0e-14 &&
                       std::abs(az.view.unchecked({3, 3, 4}) -
                                transmissibility) < 2.0e-14,
                   "operator exposes positive harmonic face transmissibilities");
  passed &= expect(certificate.valid() &&
                       certificate.scope ==
                           EquationAssemblyScope::momentum_predictor,
                   "momentum assembly publishes a typed predictor certificate");
  return passed;
}

bool test_momentum_unsteady_advection_oracle() {
  Fixture fixture;
  constexpr std::int32_t n = 8;
  bool passed = expect(make_fixture(n, false, fixture),
                       "unsteady/advection momentum fixture compiles");
  if (!passed) {
    return false;
  }
  const Int3 cells = fixture.patch.cells;
  OwnedField rho = make_field(0U, cells, 1U, 2U, 401U);
  OwnedField velocity = make_field(1U, cells, 3U, 2U, 402U);
  OwnedField accepted_velocity = make_field(1U, cells, 3U, 2U, 403U);
  OwnedField previous_velocity = make_field(1U, cells, 3U, 2U, 404U);
  OwnedField pi = make_field(2U, cells, 1U, 2U, 405U);
  OwnedField h = make_field(3U, cells, 1U, 2U, 406U);
  OwnedField temperature = make_field(4U, cells, 1U, 2U, 407U);
  OwnedField viscosity = make_field(5U, cells, 1U, 1U, 408U);
  OwnedField gradients = make_field(7U, cells, 9U, 1U, 409U);
  OwnedField diagonal = make_field(40U, cells, 3U, 0U, 410U);
  OwnedField rhs = make_field(41U, cells, 3U, 0U, 411U);
  OwnedField residual = make_field(42U, cells, 3U, 0U, 412U);
  OwnedFaceField ax = make_face_field(CartesianAxis::x, cells, 413U);
  OwnedFaceField ay = make_face_field(CartesianAxis::y, cells, 414U);
  OwnedFaceField az = make_face_field(CartesianAxis::z, cells, 415U);

  constexpr std::array<double, 3U> base{0.70, -0.35, 0.21};
  constexpr std::array<double, 3U> slope_x{0.16, -0.08, 0.04};
  constexpr std::array<double, 3U> slope_y{-0.05, 0.11, 0.07};
  constexpr std::array<double, 3U> slope_z{0.09, 0.03, -0.06};
  constexpr std::array<double, 3U> temporal_delta{0.025, -0.014, 0.009};
  constexpr double density = 1.25;
  constexpr double mu = 0.031;
  const double spacing = 1.0 / static_cast<double>(n);
  for (std::int32_t k = -2; k < cells.z + 2; ++k) {
    const double z = (static_cast<double>(k) + 0.5) * spacing;
    for (std::int32_t j = -2; j < cells.y + 2; ++j) {
      const double y = (static_cast<double>(j) + 0.5) * spacing;
      for (std::int32_t i = -2; i < cells.x + 2; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * spacing;
        const Int3 cell{i, j, k};
        rho.view.unchecked(cell, 0U) = density;
        pi.view.unchecked(cell, 0U) = 0.0;
        h.view.unchecked(cell, 0U) = 3.0e5;
        temperature.view.unchecked(cell, 0U) = 300.0;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double q = base[component] + slope_x[component] * x +
                           slope_y[component] * y + slope_z[component] * z;
          velocity.view.unchecked(cell, component) = q;
          accepted_velocity.view.unchecked(cell, component) =
              q - temporal_delta[component];
          previous_velocity.view.unchecked(cell, component) =
              q - 2.0 * temporal_delta[component];
        }
      }
    }
  }
  for (std::int32_t k = -1; k < cells.z + 1; ++k) {
    for (std::int32_t j = -1; j < cells.y + 1; ++j) {
      for (std::int32_t i = -1; i < cells.x + 1; ++i) {
        const Int3 cell{i, j, k};
        viscosity.view.unchecked(cell, 0U) = mu;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          gradients.view.unchecked(cell, 3U * component) = slope_x[component];
          gradients.view.unchecked(cell, 3U * component + 1U) =
              slope_y[component];
          gradients.view.unchecked(cell, 3U * component + 2U) =
              slope_z[component];
        }
      }
    }
  }

  FaceFluxStorage flux_storage;
  FaceFluxView flux;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, flux_storage)) &&
                       static_cast<bool>(flux_storage.workspace_view(
                           0U, 501U, flux)),
                   "nonzero provisional mass flux allocates");
  if (!passed) {
    return false;
  }
  constexpr double mass_x = 0.40;
  constexpr double mass_y = -0.23;
  constexpr double mass_z = 0.17;
  for (std::int32_t k = 0; k < flux.x.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.x.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.x.extents.x; ++i) {
        flux.x.unchecked({i, j, k}) = mass_x;
      }
    }
  }
  for (std::int32_t k = 0; k < flux.y.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.y.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.y.extents.x; ++i) {
        flux.y.unchecked({i, j, k}) = mass_y;
      }
    }
  }
  for (std::int32_t k = 0; k < flux.z.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.z.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.z.extents.x; ++i) {
        flux.z.unchecked({i, j, k}) = mass_z;
      }
    }
  }

  EquationStateView state;
  state.density = {as_const(rho.view), as_const(rho.view),
                   as_const(rho.view)};
  state.velocity = {as_const(velocity.view), as_const(accepted_velocity.view),
                    as_const(previous_velocity.view)};
  state.pressure_perturbation = {as_const(pi.view), as_const(pi.view),
                                 as_const(pi.view)};
  state.enthalpy = {as_const(h.view), as_const(h.view), as_const(h.view)};
  state.temperature = {as_const(temperature.view), as_const(temperature.view),
                       as_const(temperature.view)};
  state.pressure_reference = 101325.0;
  state.accepted_pressure_reference = 101325.0;
  state.previous_pressure_reference = 101325.0;
  EquationMaterialView material;
  material.effective_viscosity = as_const(viscosity.view);
  EquationAssemblyContext context;
  context.dt = 0.1;
  context.bdf = {15.0, -20.0, 5.0, 2U};
  context.time = 500U;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = 501U;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::momentum_predictor;
  context.mass_flux = as_const(flux);
  context.provisional_mass_flux = true;
  EquationSystemView system{diagonal.view, rhs.view, residual.view,
                            ax.view, ay.view, az.view};
  EquationAssemblyCertificate certificate;
  const Status momentum_status = assemble_momentum(
      fixture.equations.momentum(), state, material, as_const(gradients.view),
      {}, context, system, certificate);
  passed &= expect(static_cast<bool>(momentum_status),
                   "production momentum assembles nonzero BDF2 and advection");
  if (!passed) {
    return false;
  }

  const Int3 cell{3, 3, 3};
  const double volume = spacing * spacing * spacing;
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    const double unsteady = 10.0 * density * temporal_delta[component];
    const double advection =
        spacing * (mass_x * slope_x[component] +
                   mass_y * slope_y[component] +
                   mass_z * slope_z[component]);
    const double expected = unsteady * volume + advection;
    const double actual = residual.view.unchecked(cell, component);
    passed &= expect(std::abs(actual - expected) <=
                         4.0e-14 * std::max(1.0, std::abs(expected)),
                     "momentum residual matches exact BDF2 plus face-flux oracle");
    passed &= expect(std::abs(actual - advection) > 1.0e-5,
                     "oracle kills omitted-unsteady mutation");
    passed &= expect(std::abs(actual - unsteady * volume) > 1.0e-5,
                     "oracle kills omitted-advection mutation");
    passed &= expect(std::abs(rhs.view.unchecked(cell, component) -
                              (diagonal.view.unchecked(cell, component) *
                                   velocity.view.unchecked(cell, component) -
                               actual)) <= 4.0e-14,
                     "linear RHS remains algebraically consistent");
  }
  passed &= expect(certificate.valid(),
                   "unsteady/advection assembly publishes a certificate");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_mms_orders();
  passed &= test_full_newtonian_stress_contract();
  passed &= test_momentum_full_stress_and_operator_oracle();
  passed &= test_momentum_unsteady_advection_oracle();
  MPI_Finalize();
  return passed ? 0 : 1;
}
