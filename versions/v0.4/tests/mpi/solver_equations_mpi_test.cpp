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

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool same_collective_status(MPI_Comm communicator, Status status) {
  const std::array<std::uint64_t, 2U> local{
      static_cast<std::uint64_t>(status.code), status.detail};
  std::array<std::uint64_t, 2U> minimum{};
  std::array<std::uint64_t, 2U> maximum{};
  return MPI_Allreduce(local.data(), minimum.data(), 2, MPI_UINT64_T, MPI_MIN,
                       communicator) == MPI_SUCCESS &&
         MPI_Allreduce(local.data(), maximum.data(), 2, MPI_UINT64_T, MPI_MAX,
                       communicator) == MPI_SUCCESS &&
         minimum == maximum;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::tensor_stretched;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {17, 11, 7};
  mesh.has_base_spacing = true;
  mesh.base_spacing = {1.04 / 17.0, 1.04 / 11.0, 1.04 / 7.0};
  mesh.minimum_spacing = {0.96 / 17.0, 0.96 / 11.0, 0.96 / 7.0};
  mesh.max_growth_ratio = 1.0 + 0.8 / 17.0;
  mesh.focus_regions.push_back(
      {{0.35, 0.35, 0.35}, {0.65, 0.65, 0.65},
       {0.98 / 17.0, 0.98 / 11.0, 0.98 / 7.0}});
  mesh.limits.max_global_cells = 17U * 11U * 7U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 28U;
  return mesh;
}

ValidatedModel model(
    const CartesianMeshSpec& mesh,
    ConvectionScheme momentum = ConvectionScheme::central2) {
  ValidatedModel value;
  value.mesh = mesh;
  value.pressure_reference = PressureReferenceKind::closed_mass;
  value.fingerprint = 0x14e90123U;
  for (BoundaryFaceSpec& face : value.boundaries) {
    face.flow_kind = BoundaryKind::no_slip_wall;
    face.thermal_kind = BoundaryKind::adiabatic_wall;
  }
  value.schemes.momentum = momentum;
  value.schemes.enthalpy = ConvectionScheme::central2;
  value.schemes.species = ConvectionScheme::central2;
  value.schemes.passive_scalar = ConvectionScheme::central2;
  return value;
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
  SpeciesThermophysicalSpec air;
  air.stable_name = "air";
  air.molecular_weight = 28.96546;
  air.temperature_switch = 1000.0;
  air.nasa7_low[0U] = 3.5;
  air.nasa7_high[0U] = 3.5;
  air.viscosity_reference = 1.8e-5;
  air.conductivity = 0.026;
  spec.species.push_back(air);
  return spec;
}

struct Dependencies {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
};

bool make_dependencies(
    MPI_Comm communicator, Dependencies& out,
    ConvectionScheme momentum = ConvectionScheme::central2) {
  const CartesianMeshSpec mesh = mesh_spec();
  const ValidatedModel input = model(mesh, momentum);
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
      !registry.require_field("T", 1U, 2U, temperature) || temperature != 4U) {
    return false;
  }
  if (!CartesianGeometryCompiler::compile(communicator, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(communicator, input, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes, out.time)) {
    return false;
  }
  const ThermophysicalSpec thermo = thermo_spec();
  if (!ThermodynamicsPlan::compile(thermo, {}, out.thermodynamics) ||
      !TransportPlan::compile(thermo, out.thermodynamics, out.transport)) {
    return false;
  }
  const std::array<FieldId, 8U> fields{0U, 1U, 2U, 3U,
                                      4U, 5U, 6U, 7U};
  return out.contributions.configure({fields.data(), fields.size()}) &&
         out.contributions.freeze();
}

EquationPlanSpec plan_spec() {
  EquationPlanSpec spec;
  spec.density = 0U;
  spec.velocity = 1U;
  spec.pressure_perturbation = 2U;
  spec.enthalpy = 3U;
  spec.temperature = 4U;
  spec.effective_viscosity = 5U;
  spec.velocity_gradient = 7U;
  spec.pressure_compressibility = 6U;
  spec.pressure_reference = PressureReferenceKind::closed_mass;
  spec.closed_mass_service_stage = 1U;
  spec.maximum_cells_per_rank = 17U * 11U * 7U;
  return spec;
}

Status compile(MPI_Comm communicator, const Dependencies& input,
               const EquationPlanSpec& spec, EquationPlanSet& out,
               EquationCompileDiagnostics* diagnostics = nullptr) {
  return EquationPlanSet::compile(
      communicator, input.schemes, input.geometry, input.patch, input.boundary,
      input.contributions, input.thermodynamics, input.transport, spec, out,
      diagnostics);
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
                      std::uint8_t ghosts, RevisionToken revision,
                      StorageIdentity identity) {
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
  field.view.storage_identity = identity;
  field.view.revision_domain = 14001U;
  return field;
}

OwnedFaceField make_face(CartesianAxis axis, Int3 cells,
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
                14002U};
  return field;
}

void fill(OwnedField& field, double value) {
  std::fill(field.bytes.begin(), field.bytes.end(), value);
}

double rho_at(Int3 global) {
  return 1.0 + 0.01 * global.x + 0.02 * global.y + 0.03 * global.z;
}

double flux_x(Int3 global_face) {
  return 0.11 + 0.013 * global_face.x - 0.007 * global_face.y +
         0.005 * global_face.z;
}

double flux_y(Int3 global_face) {
  return -0.09 + 0.004 * global_face.x + 0.017 * global_face.y -
         0.006 * global_face.z;
}

double flux_z(Int3 global_face) {
  return 0.03 - 0.008 * global_face.x + 0.009 * global_face.y +
         0.019 * global_face.z;
}

bool test_distributed_continuity_global_id_oracle(MPI_Comm world, int rank) {
  Dependencies dependencies;
  EquationPlanSet plan;
  bool passed = expect(make_dependencies(world, dependencies) &&
                           static_cast<bool>(compile(
                               world, dependencies, plan_spec(), plan)),
                       rank, "tensor-stretched equation plan compiles for numerical oracle");
  if (!passed) {
    return false;
  }
  const MeshPatch patch = dependencies.patch;
  const Int3 cells = patch.cells;
  OwnedField rho = make_field(0U, cells, 1U, 0U, 101U,
                              10000U + static_cast<unsigned>(rank));
  OwnedField residual = make_field(30U, cells, 1U, 0U, 102U,
                                   11000U + static_cast<unsigned>(rank));
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        rho.view.unchecked({i, j, k}, 0U) =
            rho_at({patch.begin.x + i, patch.begin.y + j,
                    patch.begin.z + k});
      }
    }
  }
  FaceFluxStorage storage;
  FaceFluxView flux;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, storage)) &&
                       static_cast<bool>(storage.workspace_view(0U, 777U,
                                                                flux)),
                   rank, "rank-local flux workspace allocates");
  if (!passed) {
    return false;
  }
  for (std::int32_t k = 0; k < flux.x.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.x.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.x.extents.x; ++i) {
        flux.x.unchecked({i, j, k}) = flux_x(
            {patch.begin.x + i, patch.begin.y + j, patch.begin.z + k});
      }
    }
  }
  for (std::int32_t k = 0; k < flux.y.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.y.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.y.extents.x; ++i) {
        flux.y.unchecked({i, j, k}) = flux_y(
            {patch.begin.x + i, patch.begin.y + j, patch.begin.z + k});
      }
    }
  }
  for (std::int32_t k = 0; k < flux.z.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.z.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.z.extents.x; ++i) {
        flux.z.unchecked({i, j, k}) = flux_z(
            {patch.begin.x + i, patch.begin.y + j, patch.begin.z + k});
      }
    }
  }
  EquationStateView state;
  state.density = {as_const(rho.view), as_const(rho.view), as_const(rho.view)};
  EquationAssemblyContext context;
  context.dt = 0.2;
  context.bdf = {5.0, -5.0, 0.0, 1U};
  context.time = 103U;
  context.geometry = dependencies.geometry.topology_revision();
  context.face_flux = 777U;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::momentum_predictor;
  context.mass_flux = as_const(flux);
  context.provisional_mass_flux = true;
  EquationSystemView system;
  system.residual = residual.view;
  EquationAssemblyCertificate certificate;
  passed &= expect(static_cast<bool>(assemble_continuity(
                       plan.continuity(), state, context, system,
                       certificate)),
                   rank, "distributed production continuity assembler succeeds");
  if (!passed) {
    return false;
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 local{i, j, k};
        const Int3 global{patch.begin.x + i, patch.begin.y + j,
                          patch.begin.z + k};
        const double expected =
            flux_x({global.x + 1, global.y, global.z}) - flux_x(global) +
            flux_y({global.x, global.y + 1, global.z}) - flux_y(global) +
            flux_z({global.x, global.y, global.z + 1}) - flux_z(global);
        const double actual = residual.view.unchecked(local, 0U);
        passed &= expect(std::isfinite(actual) &&
                             std::abs(actual - expected) <=
                                 2.0e-13 *
                                     std::max(1.0, std::abs(expected)),
                         rank, "local output matches COMM_SELF global-ID face oracle");
      }
    }
  }
  return passed;
}

bool test_distributed_momentum_tensor_stretched_oracle(MPI_Comm world,
                                                              int rank) {
  Dependencies dependencies;
  EquationPlanSet plan;
  bool passed = expect(make_dependencies(world, dependencies,
                                         ConvectionScheme::limited_central2) &&
                           static_cast<bool>(compile(
                               world, dependencies, plan_spec(), plan)),
                       rank,
                       "tensor-stretched equation plan compiles for momentum oracle");
  if (!passed) {
    return false;
  }
  const MeshPatch patch = dependencies.patch;
  const Int3 cells = patch.cells;
  const std::uint8_t reach = plan.kernels().reach();
  OwnedField rho = make_field(0U, cells, 1U, 0U, 201U,
                              12000U + static_cast<unsigned>(rank));
  OwnedField velocity = make_field(1U, cells, 3U, reach, 202U,
                                   13000U + static_cast<unsigned>(rank));
  OwnedField pressure = make_field(2U, cells, 1U, 1U, 203U,
                                   14000U + static_cast<unsigned>(rank));
  OwnedField viscosity = make_field(5U, cells, 1U, 1U, 204U,
                                    15000U + static_cast<unsigned>(rank));
  OwnedField gradient = make_field(7U, cells, 9U, 1U, 205U,
                                   16000U + static_cast<unsigned>(rank));
  OwnedField residual = make_field(31U, cells, 3U, 0U, 206U,
                                   17000U + static_cast<unsigned>(rank));
  constexpr std::array<double, 3U> uniform_velocity{0.4, -0.25, 0.17};
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        rho.view.unchecked({i, j, k}, 0U) = 1.18;
      }
    }
  }
  for (std::int32_t k = -static_cast<std::int32_t>(reach);
       k < cells.z + reach; ++k) {
    for (std::int32_t j = -static_cast<std::int32_t>(reach);
         j < cells.y + reach; ++j) {
      for (std::int32_t i = -static_cast<std::int32_t>(reach);
           i < cells.x + reach; ++i) {
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          velocity.view.unchecked({i, j, k}, component) =
              uniform_velocity[component];
        }
      }
    }
  }
  for (std::int32_t k = -1; k < cells.z + 1; ++k) {
    for (std::int32_t j = -1; j < cells.y + 1; ++j) {
      for (std::int32_t i = -1; i < cells.x + 1; ++i) {
        const Int3 cell{i, j, k};
        pressure.view.unchecked(cell, 0U) = 0.0;
        viscosity.view.unchecked(cell, 0U) = 0.021;
        for (std::uint8_t component = 0U; component < 9U; ++component) {
          gradient.view.unchecked(cell, component) = 0.0;
        }
      }
    }
  }

  FaceFluxStorage storage;
  FaceFluxView flux;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                           cells, 1U, storage)) &&
                       static_cast<bool>(storage.workspace_view(0U, 888U,
                                                                flux)),
                   rank, "momentum oracle flux workspace allocates");
  if (!passed) {
    return false;
  }
  for (std::int32_t k = 0; k < flux.x.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.x.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.x.extents.x; ++i) {
        flux.x.unchecked({i, j, k}) = flux_x(
            {patch.begin.x + i, patch.begin.y + j, patch.begin.z + k});
      }
    }
  }
  for (std::int32_t k = 0; k < flux.y.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.y.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.y.extents.x; ++i) {
        flux.y.unchecked({i, j, k}) = flux_y(
            {patch.begin.x + i, patch.begin.y + j, patch.begin.z + k});
      }
    }
  }
  for (std::int32_t k = 0; k < flux.z.extents.z; ++k) {
    for (std::int32_t j = 0; j < flux.z.extents.y; ++j) {
      for (std::int32_t i = 0; i < flux.z.extents.x; ++i) {
        flux.z.unchecked({i, j, k}) = flux_z(
            {patch.begin.x + i, patch.begin.y + j, patch.begin.z + k});
      }
    }
  }

  EquationStateView state;
  state.density = {as_const(rho.view), as_const(rho.view), as_const(rho.view)};
  state.velocity = {as_const(velocity.view), as_const(velocity.view),
                    as_const(velocity.view)};
  state.pressure_perturbation = {as_const(pressure.view),
                                 as_const(pressure.view),
                                 as_const(pressure.view)};
  EquationMaterialView material;
  material.effective_viscosity = as_const(viscosity.view);
  EquationAssemblyContext context;
  context.dt = 0.2;
  context.bdf = {5.0, -5.0, 0.0, 1U};
  context.time = 207U;
  context.geometry = dependencies.geometry.topology_revision();
  context.boundary = dependencies.boundary.revision();
  context.transport = dependencies.transport.fingerprint();
  context.face_flux = 888U;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::momentum_predictor;
  context.mass_flux = as_const(flux);
  context.provisional_mass_flux = true;
  EquationSystemView system;
  system.residual = residual.view;
  EquationAssemblyCertificate certificate;
  passed &= expect(static_cast<bool>(assemble_momentum(
                       plan.momentum(), state, material,
                       as_const(gradient.view), {}, context, system,
                       certificate)),
                   rank,
                   "distributed tensor-stretched production momentum assembles");
  if (!passed) {
    return false;
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 local{i, j, k};
        const Int3 global{patch.begin.x + i, patch.begin.y + j,
                          patch.begin.z + k};
        const double divergence =
            flux_x({global.x + 1, global.y, global.z}) - flux_x(global) +
            flux_y({global.x, global.y + 1, global.z}) - flux_y(global) +
            flux_z({global.x, global.y, global.z + 1}) - flux_z(global);
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double expected = uniform_velocity[component] * divergence;
          const double actual = residual.view.unchecked(local, component);
          passed &= expect(std::isfinite(actual) &&
                               std::abs(actual - expected) <=
                                   5.0e-13 *
                                       std::max(1.0, std::abs(expected)),
                           rank,
                           "distributed momentum matches global-ID advection oracle");
        }
      }
    }
  }
  passed &= expect(certificate.valid(), rank,
                   "distributed momentum publishes a typed certificate");
  const std::vector<double> residual_only = residual.bytes;
  OwnedField linear_diagonal = make_field(
      32U, cells, 3U, 0U, 207U, 18000U + static_cast<unsigned>(rank));
  OwnedField linear_rhs = make_field(
      33U, cells, 3U, 0U, 208U, 19000U + static_cast<unsigned>(rank));
  OwnedField low_order_rhs_delta = make_field(
      41U, cells, 3U, 0U, 209U, 19500U + static_cast<unsigned>(rank));
  OwnedFaceField linear_x = make_face(CartesianAxis::x, cells,
                                      20000U + static_cast<unsigned>(rank));
  OwnedFaceField linear_y = make_face(CartesianAxis::y, cells,
                                      21000U + static_cast<unsigned>(rank));
  OwnedFaceField linear_z = make_face(CartesianAxis::z, cells,
                                      22000U + static_cast<unsigned>(rank));
  EquationSystemView linear_system;
  linear_system.diagonal = linear_diagonal.view;
  linear_system.rhs = linear_rhs.view;
  linear_system.residual = residual.view;
  linear_system.x_coefficient = linear_x.view;
  linear_system.y_coefficient = linear_y.view;
  linear_system.z_coefficient = linear_z.view;
  EquationAssemblyCertificate linear_certificate;
  const bool linear_status = static_cast<bool>(assemble_momentum_predictor(
      plan.momentum(), state, material, as_const(gradient.view), {}, context,
      linear_system, low_order_rhs_delta.view, linear_certificate));
  bool future_linear_ok = linear_status && linear_certificate.valid() &&
                          residual.bytes == residual_only;
  bool has_high_cfl_cell = false;
  const auto residual_offset = [&](Int3 cell,
                                   std::uint8_t component) noexcept {
    return static_cast<std::size_t>(cell.x) +
           static_cast<std::size_t>(cell.y) * residual.view.stride_y +
           static_cast<std::size_t>(cell.z) * residual.view.stride_z +
           static_cast<std::size_t>(component) * residual.view.component_stride;
  };
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const double a_conv =
            std::max(flux.x.unchecked({i + 1, j, k}), 0.0) +
            std::max(-flux.x.unchecked({i, j, k}), 0.0) +
            std::max(flux.y.unchecked({i, j + 1, k}), 0.0) +
            std::max(-flux.y.unchecked({i, j, k}), 0.0) +
            std::max(flux.z.unchecked({i, j, k + 1}), 0.0) +
            std::max(-flux.z.unchecked({i, j, k}), 0.0);
        const double cell_volume =
            dependencies.geometry.x().widths().data[
                static_cast<std::size_t>(patch.begin.x + i)] *
            dependencies.geometry.y().widths().data[
                static_cast<std::size_t>(patch.begin.y + j)] *
            dependencies.geometry.z().widths().data[
                static_cast<std::size_t>(patch.begin.z + k)];
        const double a_base =
            5.0 * 1.18 * cell_volume +
            linear_x.view.unchecked({i, j, k}) +
            linear_x.view.unchecked({i + 1, j, k}) +
            linear_y.view.unchecked({i, j, k}) +
            linear_y.view.unchecked({i, j + 1, k}) +
            linear_z.view.unchecked({i, j, k}) +
            linear_z.view.unchecked({i, j, k + 1});
        has_high_cfl_cell = has_high_cfl_cell || a_conv > a_base;
        const double expected_diagonal = a_base + a_conv;
        bool cell_linear_ok = linear_status;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          if (!linear_status) {
            continue;
          }
          const double actual_diagonal =
              linear_diagonal.view.unchecked(cell, component);
          const double expected_rhs =
              expected_diagonal * uniform_velocity[component] -
              residual_only[residual_offset(cell, component)];
          const double actual_rhs = linear_rhs.view.unchecked(cell, component);
          cell_linear_ok =
              cell_linear_ok && std::isfinite(actual_diagonal) &&
              std::isfinite(actual_rhs) &&
              std::abs(low_order_rhs_delta.view.unchecked(cell, component)) <=
                  5.0e-13 &&
              std::abs(actual_diagonal - expected_diagonal) <=
                  5.0e-13 * std::max(1.0, std::abs(expected_diagonal)) &&
              std::abs(actual_rhs - expected_rhs) <=
                  5.0e-13 * std::max(1.0, std::abs(expected_rhs));
        }
        future_linear_ok = future_linear_ok && cell_linear_ok;
      }
    }
  }
  passed &= expect(
      future_linear_ok && has_high_cfl_cell, rank,
      "additive outgoing-mass implicit-upwind diagonal oracle");
  const std::vector<double> full_residual = residual.bytes;
  std::fill(residual.bytes.begin(), residual.bytes.end(), -901.0);
  const std::int32_t split = cells.x / 2;
  const KernelBox lower{{0, 0, 0}, {split, cells.y, cells.z}};
  const KernelBox upper{{split, 0, 0},
                        {cells.x - split, cells.y, cells.z}};
  AssemblyEpoch epoch;
  EquationAssemblyCertificate tiled = certificate;
  passed &= expect(
      split > 0 && split < cells.x &&
          static_cast<bool>(epoch.begin(context, system)) &&
          static_cast<bool>(assemble_tile(
              epoch, plan.momentum(), state, material,
              as_const(gradient.view), {}, upper)) &&
          static_cast<bool>(assemble_tile(
              epoch, plan.momentum(), state, material,
              as_const(gradient.view), {}, lower)) &&
          static_cast<bool>(epoch.finalize(tiled)) &&
          tiled.valid() && residual.bytes == full_residual,
      rank, "distributed reverse-order momentum tiles equal full assembly");

  // The assembled diagonal/RHS oracle above is intentionally insufficient:
  // the production predictor must solve the six-neighbour matrix rather than
  // return its one-point Jacobi value.  Exercise that exact public seam, then
  // form A*U-b independently after a velocity halo exchange.
  ReductionEngine momentum_reductions;
  LinearWorkspaceRequirements momentum_requirements;
  passed &= expect(
      static_cast<bool>(make_linear_workspace_requirements(
          LinearAlgorithm::fgmres, cells, 1U, 12U,
          ReductionMode::mpi_allreduce, 910U, momentum_requirements)) &&
          static_cast<bool>(ReductionEngine::compile(
              world, ReductionMode::mpi_allreduce,
              momentum_requirements.reduction_capacity,
              momentum_reductions)),
      rank, "implicit momentum Krylov resources compile");
  OwnedField momentum_vectors = make_field(
      92U, cells, momentum_requirements.vector_slots, 1U, 911U,
      23000U);
  OwnedField momentum_scalars = make_field(
      93U,
      {static_cast<std::int32_t>(momentum_requirements.scalar_doubles), 1,
       1},
      1U, 0U, 912U, 23000U);
  SolverWorkspace momentum_workspace;
  passed &= expect(
      static_cast<bool>(SolverWorkspace::bind(
          momentum_requirements, momentum_vectors.view,
          momentum_scalars.view, momentum_workspace)),
      rank, "implicit momentum Krylov workspace binds");
  const std::array<HaloFieldSpec, 1U> momentum_halo_fields{{
      {momentum_vectors.view.field, 1U, 1U}}};
  HaloEngine momentum_halo;
  passed &= expect(
      static_cast<bool>(momentum_halo.reserve(
          world, patch,
          {momentum_halo_fields.data(), momentum_halo_fields.size()},
          dependencies.boundary.halo_topology())),
      rank, "implicit momentum operator halo reserves");
  if (!passed) return false;

  std::vector<double> diagonal_quotient(
      static_cast<std::size_t>(cells.x) * cells.y * cells.z * 3U);
  const std::vector<double> assembled_momentum_diagonal =
      linear_diagonal.bytes;
  const std::vector<double> assembled_momentum_rhs = linear_rhs.bytes;
  const auto packed_offset = [cells](Int3 cell,
                                     std::uint8_t component) noexcept {
    return static_cast<std::size_t>(cell.x) +
           static_cast<std::size_t>(cells.x) *
               (static_cast<std::size_t>(cell.y) +
                static_cast<std::size_t>(cells.y) *
                    (static_cast<std::size_t>(cell.z) +
                     static_cast<std::size_t>(cells.z) * component));
  };
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          diagonal_quotient[packed_offset(cell, component)] =
              linear_rhs.view.unchecked(cell, component) /
              linear_diagonal.view.unchecked(cell, component);
        }
      }
    }
  }
  MomentumPredictorLimiterReport limiter_report;
  Status solve_status = limit_momentum_predictor_correction(
      as_const(velocity.view), as_const(flux), {}, linear_system,
      as_const(low_order_rhs_delta.view),
      momentum_reductions, limiter_report);
  MomentumPredictorSolveReport solve_report;
  ResourceCounters solve_resources;
  if (solve_status) {
    solve_status = solve_momentum_predictor(
        world, plan.momentum(), dependencies.boundary, patch,
        linear_certificate, as_const(flux), {}, linear_system, velocity.view,
        momentum_halo, momentum_workspace, momentum_reductions,
        &solve_resources, solve_report);
  }
  bool solve_results_valid = solve_report.solve_calls == 3U;
  std::uint64_t reported_iterations = 0U;
  for (const LinearSolveResult& result : solve_report.components) {
    solve_results_valid =
        solve_results_valid && result.status &&
        (result.termination == LinearTermination::converged ||
         result.termination == LinearTermination::zero_rhs) &&
        std::isfinite(result.final_true_residual) &&
        result.final_true_residual <= result.initial_true_residual;
    reported_iterations += result.iterations;
  }
  passed &= expect(static_cast<bool>(solve_status) &&
                       same_collective_status(world, solve_status) &&
                       solve_results_valid &&
                       reported_iterations == solve_resources.linear_iterations,
                   rank, "three implicit momentum components converge and report work");
  if (!passed) {
    std::cerr << "rank " << rank << " momentum solve status="
              << static_cast<unsigned>(solve_status.code) << '/'
              << solve_status.detail << " calls="
              << static_cast<unsigned>(solve_report.solve_calls)
              << " limiter=" << limiter_report.limited << '/'
              << limiter_report.theta << " iterations="
              << reported_iterations << '/'
              << solve_resources.linear_iterations;
    for (std::size_t component = 0U;
         component < solve_report.components.size(); ++component) {
      const LinearSolveResult& result = solve_report.components[component];
      std::cerr << " c" << component << '='
                << static_cast<unsigned>(result.status.code) << '/'
                << result.status.detail << '/'
                << static_cast<unsigned>(result.termination) << '/'
                << result.iterations << '/' << result.initial_true_residual
                << '/' << result.final_true_residual << '/'
                << result.operator_applies << '/' << result.reduction_calls;
    }
    std::cerr << '\n';
  }
  if (!passed) return false;

  const std::array<HaloFieldSpec, 1U> verification_halo_fields{{
      {velocity.view.field, 1U, 3U}}};
  HaloEngine verification_halo;
  passed &= expect(
      static_cast<bool>(verification_halo.reserve(
          world, patch,
          {verification_halo_fields.data(), verification_halo_fields.size()},
          dependencies.boundary.halo_topology())),
      rank, "solved velocity verification halo reserves");
  std::array<FieldView, 1U> velocity_fields{velocity.view};
  HaloTicket velocity_ticket;
  Status velocity_halo_status = verification_halo.begin(
      32U, {velocity_fields.data(), velocity_fields.size()}, velocity_ticket);
  if (velocity_halo_status) {
    velocity_halo_status = verification_halo.finish(
        velocity_ticket, {velocity_fields.data(), velocity_fields.size()});
  }
  passed &= expect(static_cast<bool>(velocity_halo_status) &&
                       same_collective_status(world, velocity_halo_status),
                   rank, "solved velocity halo exchange succeeds");
  velocity.view = velocity_fields[0U];

  double local_residual_squared = 0.0;
  double local_rhs_squared = 0.0;
  double local_jacobi_difference_squared = 0.0;
  const auto coefficient = [&](CartesianAxis axis, Int3 face) noexcept {
    return axis == CartesianAxis::x
               ? linear_x.view.unchecked(face)
               : (axis == CartesianAxis::y
                      ? linear_y.view.unchecked(face)
                      : linear_z.view.unchecked(face));
  };
  const auto face_mass = [&](CartesianAxis axis, Int3 face) noexcept {
    return axis == CartesianAxis::x
               ? flux.x.unchecked(face)
               : (axis == CartesianAxis::y ? flux.y.unchecked(face)
                                             : flux.z.unchecked(face));
  };
  const auto shifted_cell = [](Int3 cell, CartesianAxis axis,
                               int direction) noexcept {
    if (axis == CartesianAxis::x)
      cell.x += direction;
    else if (axis == CartesianAxis::y)
      cell.y += direction;
    else
      cell.z += direction;
    return cell;
  };
  const Int3 global_cells = dependencies.geometry.global_cells();
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double centre = velocity.view.unchecked(cell, component);
          double neighbor_sum = 0.0;
          for (CartesianAxis axis : {CartesianAxis::x, CartesianAxis::y,
                                     CartesianAxis::z}) {
            for (int direction : {-1, 1}) {
              const Int3 neighbor = shifted_cell(cell, axis, direction);
              const Int3 face =
                  direction < 0 ? cell : shifted_cell(cell, axis, 1);
              const double mass = face_mass(axis, face);
              const double incoming =
                  direction < 0 ? std::max(mass, 0.0)
                                : std::max(-mass, 0.0);
              const Int3 global_neighbor{patch.begin.x + neighbor.x,
                                         patch.begin.y + neighbor.y,
                                         patch.begin.z + neighbor.z};
              const std::int32_t global_coordinate =
                  axis == CartesianAxis::x
                      ? global_neighbor.x
                      : (axis == CartesianAxis::y ? global_neighbor.y
                                                   : global_neighbor.z);
              const std::int32_t global_extent =
                  axis == CartesianAxis::x
                      ? global_cells.x
                      : (axis == CartesianAxis::y ? global_cells.y
                                                   : global_cells.z);
              // All six model boundaries in this fixture are no-slip walls,
              // whose homogeneous ghost perturbation is -U_centre.
              const double neighbor_value =
                  global_coordinate < 0 || global_coordinate >= global_extent
                      ? -centre
                      : velocity.view.unchecked(neighbor, component);
              neighbor_sum +=
                  (coefficient(axis, face) + incoming) * neighbor_value;
            }
          }
          const double rhs_value =
              linear_rhs.view.unchecked(cell, component);
          const double residual_value =
              linear_diagonal.view.unchecked(cell, component) * centre -
              neighbor_sum - rhs_value;
          const double jacobi_difference =
              centre - diagonal_quotient[packed_offset(cell, component)];
          local_residual_squared += residual_value * residual_value;
          local_rhs_squared += rhs_value * rhs_value;
          local_jacobi_difference_squared +=
              jacobi_difference * jacobi_difference;
        }
      }
    }
  }
  std::array<double, 3U> local_norms{local_residual_squared,
                                     local_rhs_squared,
                                     local_jacobi_difference_squared};
  std::array<double, 3U> global_norms{};
  const int norm_status = MPI_Allreduce(local_norms.data(),
                                        global_norms.data(), 3, MPI_DOUBLE,
                                        MPI_SUM, world);
  const double relative_residual =
      std::sqrt(global_norms[0U]) /
      std::max(1.0e-30, std::sqrt(global_norms[1U]));
  passed &= expect(norm_status == MPI_SUCCESS &&
                       std::isfinite(relative_residual) &&
                       relative_residual <= 1.05e-4 &&
                       global_norms[2U] > 1.0e-24,
                   rank,
                   "implicit momentum solution satisfies independent A*U-b oracle and differs from diagonal quotient");

  // An inactive IBM cell is an identity/zero row outside the active-fluid
  // operator. Its caller-side initial value must not enter a neighbouring
  // fluid RHS or alter the fluid Krylov trajectory.
  const std::size_t local_cell_count =
      static_cast<std::size_t>(cells.x) * cells.y * cells.z;
  const Int3 inactive_cell{cells.x / 2, cells.y / 2, cells.z / 2};
  const std::size_t inactive_offset = packed_offset(inactive_cell, 0U);
  std::vector<std::uint8_t> active(local_cell_count, 1U);
  active[inactive_offset] = 0U;
  const Int3 x_face_extents{cells.x + 1, cells.y, cells.z};
  const Int3 y_face_extents{cells.x, cells.y + 1, cells.z};
  const Int3 z_face_extents{cells.x, cells.y, cells.z + 1};
  const auto face_count = [](Int3 extents) noexcept {
    return static_cast<std::size_t>(extents.x) * extents.y * extents.z;
  };
  const auto face_offset = [](Int3 extents, Int3 face) noexcept {
    return static_cast<std::size_t>(face.x) +
           static_cast<std::size_t>(extents.x) *
               (static_cast<std::size_t>(face.y) +
                static_cast<std::size_t>(extents.y) *
                    static_cast<std::size_t>(face.z));
  };
  std::vector<std::uint8_t> active_x_faces(face_count(x_face_extents), 1U);
  std::vector<std::uint8_t> active_y_faces(face_count(y_face_extents), 1U);
  std::vector<std::uint8_t> active_z_faces(face_count(z_face_extents), 1U);
  active_x_faces[face_offset(x_face_extents, inactive_cell)] = 0U;
  active_x_faces[face_offset(
      x_face_extents,
      {inactive_cell.x + 1, inactive_cell.y, inactive_cell.z})] = 0U;
  active_y_faces[face_offset(y_face_extents, inactive_cell)] = 0U;
  active_y_faces[face_offset(
      y_face_extents,
      {inactive_cell.x, inactive_cell.y + 1, inactive_cell.z})] = 0U;
  active_z_faces[face_offset(z_face_extents, inactive_cell)] = 0U;
  active_z_faces[face_offset(
      z_face_extents,
      {inactive_cell.x, inactive_cell.y, inactive_cell.z + 1})] = 0U;
  const MgDomainActivityView momentum_activity{
      {active.data(), active.size()},
      {active_x_faces.data(), active_x_faces.size()},
      {active_y_faces.data(), active_y_faces.size()},
      {active_z_faces.data(), active_z_faces.size()},
      991U + static_cast<unsigned>(rank), 992U};
  const auto solve_inactive_probe = [&](double inactive_value,
                                        std::vector<double>& fluid_solution) {
    linear_diagonal.bytes = assembled_momentum_diagonal;
    linear_rhs.bytes = assembled_momentum_rhs;
    for (std::int32_t k = -static_cast<std::int32_t>(reach);
         k < cells.z + reach; ++k) {
      for (std::int32_t j = -static_cast<std::int32_t>(reach);
           j < cells.y + reach; ++j) {
        for (std::int32_t i = -static_cast<std::int32_t>(reach);
             i < cells.x + reach; ++i) {
          for (std::uint8_t component = 0U; component < 3U; ++component)
            velocity.view.unchecked({i, j, k}, component) =
                uniform_velocity[component];
        }
      }
    }
    for (std::uint8_t component = 0U; component < 3U; ++component)
      velocity.view.unchecked(inactive_cell, component) =
          inactive_value * static_cast<double>(component + 1U);
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 cell{i, j, k};
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            linear_rhs.view.unchecked(cell, component) =
                linear_diagonal.view.unchecked(cell, component) *
                uniform_velocity[component];
          }
        }
      }
    }
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      linear_diagonal.view.unchecked(inactive_cell, component) = 1.0;
      linear_rhs.view.unchecked(inactive_cell, component) = 0.0;
    }
    MomentumPredictorSolveReport probe_report;
    ResourceCounters probe_resources;
    const Status probe_status = solve_momentum_predictor(
        world, plan.momentum(), dependencies.boundary, patch,
        linear_certificate, as_const(flux),
        momentum_activity, linear_system, velocity.view,
        momentum_halo, momentum_workspace, momentum_reductions,
        &probe_resources, probe_report);
    fluid_solution.clear();
    fluid_solution.reserve((local_cell_count - 1U) * 3U);
    if (probe_status) {
      for (std::int32_t k = 0; k < cells.z; ++k) {
        for (std::int32_t j = 0; j < cells.y; ++j) {
          for (std::int32_t i = 0; i < cells.x; ++i) {
            const Int3 cell{i, j, k};
            if (packed_offset(cell, 0U) == inactive_offset) continue;
            for (std::uint8_t component = 0U; component < 3U; ++component)
              fluid_solution.push_back(
                  velocity.view.unchecked(cell, component));
          }
        }
      }
    }
    return probe_status;
  };
  std::vector<double> zero_solid_solution;
  std::vector<double> nonzero_solid_solution;
  const Status zero_solid =
      solve_inactive_probe(0.0, zero_solid_solution);
  const Status nonzero_solid =
      solve_inactive_probe(17.0, nonzero_solid_solution);
  double local_active_difference = 0.0;
  if (zero_solid_solution.size() == nonzero_solid_solution.size()) {
    for (std::size_t index = 0U; index < zero_solid_solution.size(); ++index)
      local_active_difference = std::max(
          local_active_difference,
          std::abs(zero_solid_solution[index] -
                   nonzero_solid_solution[index]));
  } else {
    local_active_difference = std::numeric_limits<double>::infinity();
  }
  double global_active_difference = 0.0;
  const int active_difference_status =
      MPI_Allreduce(&local_active_difference, &global_active_difference, 1,
                    MPI_DOUBLE, MPI_MAX, world);
  passed &= expect(
      zero_solid && nonzero_solid && active_difference_status == MPI_SUCCESS &&
          std::isfinite(global_active_difference) &&
          global_active_difference <= 5.0e-12,
      rank,
      "inactive IBM initial velocity is isolated from active momentum solve");
  return passed;
}

bool test_conservative_momentum_predictor_limiter(MPI_Comm world, int rank) {
  Dependencies dependencies;
  bool passed = expect(make_dependencies(world, dependencies,
                                          ConvectionScheme::tvd2),
                       rank, "momentum-limiter dependencies compile");
  if (!passed) return false;

  const Int3 cells = dependencies.patch.cells;
  const StorageIdentity salt = 26000U + 100U * static_cast<unsigned>(rank);
  OwnedField velocity = make_field(1U, cells, 3U, 1U, 401U, salt + 1U);
  OwnedField diagonal = make_field(30U, cells, 3U, 0U, 402U, salt + 2U);
  OwnedField rhs = make_field(31U, cells, 3U, 0U, 403U, salt + 3U);
  OwnedField low_delta = make_field(41U, cells, 3U, 0U, 404U, salt + 4U);
  OwnedField residual = make_field(42U, cells, 3U, 0U, 405U, salt + 5U);
  OwnedFaceField x = make_face(CartesianAxis::x, cells, salt + 5U);
  OwnedFaceField y = make_face(CartesianAxis::y, cells, salt + 5U);
  OwnedFaceField z = make_face(CartesianAxis::z, cells, salt + 5U);
  OwnedFaceField coefficient_x =
      make_face(CartesianAxis::x, cells, salt + 6U);
  OwnedFaceField coefficient_y =
      make_face(CartesianAxis::y, cells, salt + 7U);
  OwnedFaceField coefficient_z =
      make_face(CartesianAxis::z, cells, salt + 8U);
  const EquationSystemView limiter_system{
      diagonal.view, rhs.view, residual.view, coefficient_x.view,
      coefficient_y.view, coefficient_z.view};
  fill(velocity, 0.0);
  fill(diagonal, 1.0);
  fill(rhs, 0.0);
  fill(low_delta, 0.0);
  std::fill(x.bytes.begin(), x.bytes.end(), 1.0);
  std::fill(y.bytes.begin(), y.bytes.end(), 0.0);
  std::fill(z.bytes.begin(), z.bytes.end(), 0.0);

  const Int3 upstream{1, 0, 0};
  const Int3 downstream{2, 0, 0};
  const bool owns_pair = rank == 0 && cells.x > downstream.x;
  if (owns_pair) {
    velocity.view.unchecked(upstream, 0U) = 1.0;
    rhs.view.unchecked(upstream, 0U) = -1.0;
    low_delta.view.unchecked(upstream, 0U) = 2.0;
    rhs.view.unchecked(downstream, 0U) = 2.0;
    low_delta.view.unchecked(downstream, 0U) = -2.0;
  }

  ReductionEngine reductions;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          world, ReductionMode::mpi_allreduce, 1U, reductions)),
      rank, "momentum-limiter reduction compiles");
  if (!passed) return false;

  const double local_before = owns_pair
                                  ? rhs.view.unchecked(upstream, 0U) +
                                        rhs.view.unchecked(downstream, 0U)
                                  : 0.0;
  const ConstFaceFluxView flux{as_const(x.view), as_const(y.view),
                               as_const(z.view), 405U};
  MomentumPredictorLimiterReport report;
  const Status status = limit_momentum_predictor_correction(
      as_const(velocity.view), flux, {}, limiter_system,
      as_const(low_delta.view), reductions, report);
  const double local_after = owns_pair
                                 ? rhs.view.unchecked(upstream, 0U) +
                                       rhs.view.unchecked(downstream, 0U)
                                 : 0.0;
  passed &= expect(static_cast<bool>(status) &&
                       same_collective_status(world, status) &&
                       report.limited && report.activations == 1U &&
                       report.theta == 0.5,
                   rank, "paired correction selects one global theta");
  if (owns_pair) {
    passed &= expect(rhs.view.unchecked(upstream, 0U) == 0.0 &&
                         rhs.view.unchecked(downstream, 0U) == 1.0 &&
                         local_before == local_after,
                     rank,
                     "paired correction is bounded and conservative");
  }
  fill(low_delta, 0.0);
  const std::vector<double> fast_path_rhs = rhs.bytes;
  const LinearReductionCounters counters_before = reductions.counters();
  MomentumPredictorLimiterReport fast_path_report;
  const Status fast_path_status = limit_momentum_predictor_correction(
      as_const(velocity.view), flux, {}, limiter_system,
      as_const(low_delta.view), reductions, fast_path_report);
  const LinearReductionCounters counters_after = reductions.counters();
  passed &= expect(
      static_cast<bool>(fast_path_status) && !fast_path_report.limited &&
          fast_path_report.activations == 0U &&
          fast_path_report.theta == 1.0 && rhs.bytes == fast_path_rhs &&
          counters_after.calls == counters_before.calls + 1U &&
          counters_after.blocking_operations ==
              counters_before.blocking_operations + 1U,
      rank, "theta-one limiter is byte-identical with one global reduction");

  // A bounded point endpoint is not sufficient when the converged implicit
  // row is not an M-matrix. Two strong incoming faces create an exact active
  // neighbour sum of 15 against a unit diagonal.
  fill(velocity, 0.0);
  fill(diagonal, 1.0);
  fill(rhs, 0.0);
  fill(low_delta, 0.0);
  std::fill(x.bytes.begin(), x.bytes.end(), 0.0);
  std::fill(y.bytes.begin(), y.bytes.end(), 0.0);
  std::fill(z.bytes.begin(), z.bytes.end(), 0.0);
  const Int3 convergent_cell{1, 1, 1};
  const bool owns_convergent_cell = rank == 0 &&
                                    cells.x > convergent_cell.x + 1 &&
                                    cells.y > convergent_cell.y &&
                                    cells.z > convergent_cell.z;
  if (owns_convergent_cell) {
    velocity.view.unchecked(convergent_cell, 0U) = 2.0;
    rhs.view.unchecked(convergent_cell, 0U) = 3.0;
    x.view.unchecked(convergent_cell) = 8.0;
    x.view.unchecked(
        {convergent_cell.x + 1, convergent_cell.y, convergent_cell.z}) = -7.0;
  }
  MomentumPredictorLimiterReport majorant_report;
  const Status majorant_status = limit_momentum_predictor_correction(
      as_const(velocity.view), as_const(FaceFluxView{
                                   x.view, y.view, z.view, 405U, {}}),
      {}, limiter_system, as_const(low_delta.view), reductions,
      majorant_report);
  passed &= expect(
      majorant_status && !majorant_report.limited &&
          (!owns_convergent_cell ||
           (diagonal.view.unchecked(convergent_cell, 0U) == 15.0 &&
            rhs.view.unchecked(convergent_cell, 0U) == 31.0 &&
            diagonal.view.unchecked(convergent_cell, 0U) *
                        velocity.view.unchecked(convergent_cell, 0U) -
                    rhs.view.unchecked(convergent_cell, 0U) ==
                -1.0)),
      rank,
      "high-inflow active row receives residual-preserving exact M-matrix majorant");
  return passed;
}

PisoPlanSpec piso_spec() {
  PisoPlanSpec spec;
  spec.pressure_correctors = 2U;
  spec.pressure_stage = 211U;
  spec.final_flux_slot = 0U;
  spec.pressure_solve = {1.0e-15, 1.0e-13, 400U, 4U, 12U};
  spec.eos_tolerance = 1.0e-10;
  spec.continuity_tolerance = 1.0e-10;
  spec.closed_mass_tolerance = 1.0e-10;
  spec.gauge_tolerance = 1.0e-12;
  return spec;
}

std::int32_t axis_value(Int3 value, CartesianAxis axis) {
  return axis == CartesianAxis::x
             ? value.x
             : (axis == CartesianAxis::y ? value.y : value.z);
}

bool test_piso_intermediate_halo_oracle(MPI_Comm world, int rank) {
  Dependencies dependencies;
  EquationPlanSet equations;
  bool passed = expect(
      make_dependencies(world, dependencies) &&
          static_cast<bool>(compile(world, dependencies, plan_spec(),
                                    equations)),
      rank, "distributed PISO dependencies compile");
  if (!passed) {
    return false;
  }
  PisoPlan piso;
  const Status piso_status =
      PisoPlan::compile(world, equations, piso_spec(), piso);
  passed &= expect(static_cast<bool>(piso_status) &&
                       same_collective_status(world, piso_status),
                   rank, "exactly-two PISO plan compiles collectively");
  if (!passed) {
    return false;
  }

  const MeshPatch patch = dependencies.patch;
  const Int3 cells = patch.cells;
  const std::uint8_t reach = equations.kernels().reach();
  const StorageIdentity salt = 30000U + 100U * static_cast<unsigned>(rank);
  OwnedField density = make_field(0U, cells, 1U, reach, 301U, salt + 1U);
  OwnedField boundary_pressure =
      make_field(2U, cells, 1U, reach, 318U, salt + 11U);
  OwnedField boundary_enthalpy =
      make_field(3U, cells, 1U, reach, 319U, salt + 12U);
  OwnedField boundary_temperature =
      make_field(4U, cells, 1U, reach, 320U, salt + 13U);
  std::array<OwnedField, 6U> boundary_thermo_aux{
      make_field(50U, cells, 1U, reach, 321U, salt + 14U),
      make_field(51U, cells, 1U, reach, 322U, salt + 15U),
      make_field(52U, cells, 1U, reach, 323U, salt + 16U),
      make_field(53U, cells, 1U, reach, 324U, salt + 17U),
      make_field(54U, cells, 1U, reach, 325U, salt + 18U),
      make_field(55U, cells, 1U, reach, 326U, salt + 19U)};
  OwnedField velocity = make_field(1U, cells, 3U, 0U, 302U, salt + 2U);
  OwnedField diagonal = make_field(30U, cells, 3U, 0U, 303U, salt + 3U);
  OwnedField rhs = make_field(31U, cells, 3U, 0U, 304U, salt + 4U);
  OwnedField r_au = make_field(40U, cells, 3U, 1U, 305U, salt + 5U);
  OwnedField h_by_a =
      make_field(41U, cells, 3U, reach, 306U, salt + 6U);
  OwnedField pressure_gradient =
      make_field(42U, cells, 3U, 0U, 308U, salt + 10U);
  OwnedFaceField ax = make_face(CartesianAxis::x, cells, salt + 7U);
  OwnedFaceField ay = make_face(CartesianAxis::y, cells, salt + 8U);
  OwnedFaceField az = make_face(CartesianAxis::z, cells, salt + 9U);
  fill(density, 1.0);
  double boundary_enthalpy_value = 0.0;
  double boundary_cp = 0.0;
  double boundary_gas = 0.0;
  passed &= expect(
      static_cast<bool>(dependencies.thermodynamics.mixture_enthalpy(
          300.0, {}, boundary_enthalpy_value, boundary_cp, boundary_gas)),
      rank, "distributed wall thermophysical state evaluates");
  fill(boundary_pressure, 0.0);
  fill(boundary_enthalpy, boundary_enthalpy_value);
  fill(boundary_temperature, 300.0);
  for (OwnedField& field : boundary_thermo_aux) {
    fill(field, 1.0);
  }
  fill(velocity, 0.0);
  constexpr std::array<double, 3U> diagonals{2.0, 4.0, 8.0};
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 local{i, j, k};
        const Int3 global{patch.begin.x + i, patch.begin.y + j,
                          patch.begin.z + k};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          diagonal.view.unchecked(local, component) = diagonals[component];
          const CartesianAxis axis = static_cast<CartesianAxis>(component);
          const double predictor =
              static_cast<double>(axis_value(global, axis) + 1);
          velocity.view.unchecked(local, component) = predictor;
          rhs.view.unchecked(local, component) =
              diagonals[component] * predictor;
        }
      }
    }
  }

  FaceFluxStorage flux_storage;
  FaceFluxView phi;
  FaceFluxStorage trial_flux_storage;
  FaceFluxView trial_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(cells, 1U,
                                                            flux_storage)) &&
          static_cast<bool>(flux_storage.workspace_view(0U, 307U, phi)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              cells, 1U, trial_flux_storage)) &&
          static_cast<bool>(trial_flux_storage.workspace_view(
              0U, 312U, trial_flux)),
      rank, "distributed phiHbyA and trial-flux workspaces allocate");
  FieldRegistry flux_registry;
  FieldSchema flux_schema;
  FieldId flux_dependency = 0U;
  const std::array flux_requests{ArenaFieldRequest{
      0U, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout flux_layout;
  StateLayers flux_layers;
  AttemptTransaction flux_transaction;
  FaceFluxStorage committed_flux_storage;
  FinalFaceFluxAuthority committed_flux_authority;
  FinalFaceFluxWriter committed_flux_writer;
  ConstFaceFluxView committed_flux;
  passed &= expect(
      static_cast<bool>(flux_registry.declare_field(
          "distributed.piso_flux_history", 1U, 0U, flux_dependency)) &&
          flux_dependency == 0U &&
          static_cast<bool>(flux_registry.freeze(flux_schema)) &&
          static_cast<bool>(ArenaLayout::compile(
              flux_schema,
              {flux_requests.data(), flux_requests.size()}, flux_layout)) &&
          static_cast<bool>(StateLayers::allocate(flux_layout, flux_layers)) &&
          static_cast<bool>(AttemptTransaction::create(
              flux_layers.field_count(), 1U, flux_layers.field_count(),
              flux_transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(
              cells, committed_flux_storage)) &&
          static_cast<bool>(committed_flux_authority.claim(
              214U, 0U, flux_transaction, committed_flux_writer)) &&
          static_cast<bool>(committed_flux_writer.initialize_committed(
              committed_flux_storage, as_const(trial_flux))) &&
          static_cast<bool>(committed_flux_writer.committed(
              committed_flux_storage, committed_flux)),
      rank, "distributed committed BE face-flux authority initializes");
  const PisoCouplerWorkspace workspace{
      r_au.view, h_by_a.view, pressure_gradient.view,
      ax.view, ay.view, az.view, phi};
  const std::array<HaloFieldSpec, 3U> halo_fields{{
      {density.view.field, 1U, 1U},
      {r_au.view.field, 1U, 3U}, {h_by_a.view.field, reach, 3U}}};
  HaloEngine halo;
  const Status halo_status = halo.reserve(
      world, patch, {halo_fields.data(), halo_fields.size()},
      dependencies.boundary.halo_topology());
  passed &= expect(static_cast<bool>(halo_status), rank,
                   "distributed PISO halo reserves");
  constexpr FieldId correction_field = 90U;
  const std::array<HaloFieldSpec, 1U> correction_halo_fields{{
      {correction_field, 1U, 1U}}};
  HaloEngine correction_halo;
  const Status correction_halo_status = correction_halo.reserve(
      world, patch,
      {correction_halo_fields.data(), correction_halo_fields.size()},
      dependencies.boundary.halo_topology());
  passed &= expect(static_cast<bool>(correction_halo_status), rank,
                   "distributed pressure-correction halo reserves");
  const PisoCouplerServices services{world,
                                     &dependencies.geometry,
                                     patch,
                                     &dependencies.boundary,
                                     &dependencies.thermodynamics,
                                     &halo,
                                     212U,
                                     density.view.field,
                                     &correction_halo,
                                     213U,
                                     correction_field};
  PressureVelocityCoupler coupler;
  const Status bind_status = PressureVelocityCoupler::bind(
      piso, equations, services, workspace, coupler);
  passed &= expect(static_cast<bool>(bind_status) &&
                       same_collective_status(world, bind_status),
                   rank, "distributed PISO coupler binds collectively");
  if (!passed) {
    return false;
  }

  PisoIntermediateInput input;
  input.momentum = {equations.momentum().fingerprint(),
                    EquationAssemblyScope::momentum_predictor,
                    311U,
                    dependencies.geometry.topology_revision(),
                    312U,
                    313U};
  input.predictor.plan =
      equations.thermophysical_predictor().fingerprint();
  input.predictor.time = 311U;
  input.predictor.geometry = dependencies.geometry.topology_revision();
  input.predictor.accepted_face_flux = committed_flux.revision;
  input.predictor.committed_face_flux_authority =
      committed_flux.certificate.authority();
  input.predictor.committed_face_flux_storage =
      committed_flux.certificate.storage();
  input.predictor.committed_face_flux_revision_domain =
      committed_flux.certificate.revision_domain();
  input.predictor.predicted_density = density.view.revision;
  input.predictor.predicted_density_storage = density.view.storage_identity;
  input.predictor.predicted_density_revision_domain =
      density.view.revision_domain;
  input.predictor.paired_face_flux = trial_flux.revision;
  input.predictor.paired_face_flux_storage = trial_flux.x.storage_identity;
  input.predictor.paired_face_flux_revision_domain =
      trial_flux.x.revision_domain;
  input.predictor.state = 315U;
  input.predictor.order = 1U;
  input.pressure_reference = {
      equations.pressure_reference().fingerprint(),
      equations.thermophysical_predictor().fingerprint(),
      dependencies.thermodynamics.fingerprint(),
      316U,
      311U,
      317U,
      PressureReferenceKind::closed_mass};
  input.density = density.view;
  input.trial_velocity = as_const(velocity.view);
  input.trial_flux = as_const(trial_flux);
  input.temporal_reference = as_const(phi);
  input.committed_face_history.accepted = committed_flux;
  input.momentum_system.diagonal = diagonal.view;
  input.momentum_system.rhs = rhs.view;
  input.bdf = {5.0, -5.0, 0.0, 1U};
  input.numeric_boundary = dependencies.boundary.revision();
  input.corrector = 1U;
  const std::array<BoundaryGhostFieldAuthority, 2U>
      thermophysical_ghost_authorities{{
          make_boundary_ghost_field_authority(
              as_const(boundary_pressure.view)),
          make_boundary_ghost_field_authority(
              as_const(boundary_enthalpy.view)),
      }};
  const BoundaryThermophysicalGhostAuthority thermophysical_ghost_authority{
      327U,
      dependencies.boundary.revision(),
      dependencies.boundary.local_layout_fingerprint(),
      cells,
      {reach, reach, reach},
      reach,
      {thermophysical_ghost_authorities.data(),
       thermophysical_ghost_authorities.size()}};
  const BoundaryThermophysicalGhostContext thermophysical_ghost_context{
      input.momentum.time,
      dependencies.geometry.fingerprint(),
      input.pressure_reference.pressure_reference,
      input.numeric_boundary,
      BoundaryThermophysicalGhostPhase::corrector_one};
  BoundaryThermophysicalGhostCertificate thermophysical_ghost_certificate;
  const Status thermophysical_ghost_status =
      BoundaryThermophysicalFaceClosure::close(
          dependencies.boundary, dependencies.thermodynamics,
          dependencies.transport,
          {101325.0, as_const(boundary_pressure.view),
           as_const(boundary_enthalpy.view), {},
           thermophysical_ghost_authority},
          {density.view, boundary_temperature.view,
           boundary_thermo_aux[0U].view, boundary_thermo_aux[1U].view,
           boundary_thermo_aux[2U].view, boundary_thermo_aux[3U].view,
           boundary_thermo_aux[4U].view, boundary_thermo_aux[5U].view},
          thermophysical_ghost_context, thermophysical_ghost_certificate);
  passed &= expect(static_cast<bool>(thermophysical_ghost_status) &&
                       thermophysical_ghost_certificate.valid() &&
                       same_collective_status(world,
                                              thermophysical_ghost_status),
                   rank,
                   "distributed wall thermophysical ghosts are certified");
  input.thermophysical_boundary = {
      thermophysical_ghost_certificate,
      {101325.0, as_const(boundary_pressure.view),
       as_const(boundary_enthalpy.view), {}, as_const(density.view)}};
  PisoIntermediateCertificate certificate;
  const Status refresh = coupler.refresh(input, certificate);
  if (!refresh || !certificate.valid()) {
    std::cerr << "rank " << rank << " distributed PISO refresh status="
              << static_cast<unsigned>(refresh.code) << '/' << refresh.detail
              << " certificate=" << certificate.valid()
              << " thermo-status="
              << static_cast<unsigned>(thermophysical_ghost_status.code) << '/'
              << thermophysical_ghost_status.detail
              << " thermo-certificate="
              << thermophysical_ghost_certificate.valid()
              << " density-revision=" << density.view.revision
              << " density-ghost-revision="
              << halo.ghost_revision(density.view.field) << '\n';
  }
  passed &= expect(static_cast<bool>(refresh) && certificate.valid() &&
                       same_collective_status(world, refresh),
                   rank, "distributed PISO intermediate refresh succeeds");
  if (!passed) {
    return false;
  }

  const Int3 global_cells = dependencies.geometry.global_cells();
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    const CartesianAxis axis = static_cast<CartesianAxis>(component);
    const AxisMetrics& metric = dependencies.geometry.axis(axis);
    for (const bool high : {false, true}) {
      const std::int32_t local_normal =
          high ? axis_value(cells, axis) : 0;
      const std::int32_t global_face =
          axis_value(patch.begin, axis) + local_normal;
      const std::int32_t global_extent = axis_value(global_cells, axis);
      const bool physical = global_face == 0 || global_face == global_extent;
      const Int3 face = axis == CartesianAxis::x
                            ? Int3{local_normal, 0, 0}
                            : (axis == CartesianAxis::y
                                   ? Int3{0, local_normal, 0}
                                   : Int3{0, 0, local_normal});
      const ConstFaceFieldView flux_face =
          axis == CartesianAxis::x
              ? as_const(phi).x
              : (axis == CartesianAxis::y ? as_const(phi).y
                                           : as_const(phi).z);
      const FaceFieldView coefficient =
          axis == CartesianAxis::x
              ? ax.view
              : (axis == CartesianAxis::y ? ay.view : az.view);
      double area = 1.0;
      if (axis != CartesianAxis::x) {
        area *= dependencies.geometry.x().widths().data[
            static_cast<std::size_t>(patch.begin.x)];
      }
      if (axis != CartesianAxis::y) {
        area *= dependencies.geometry.y().widths().data[
            static_cast<std::size_t>(patch.begin.y)];
      }
      if (axis != CartesianAxis::z) {
        area *= dependencies.geometry.z().widths().data[
            static_cast<std::size_t>(patch.begin.z)];
      }
      double expected_flux = 0.0;
      double expected_coefficient = 0.0;
      if (!physical) {
        const double left = static_cast<double>(global_face);
        const double right = static_cast<double>(global_face + 1);
        const double face_location =
            metric.faces().data[static_cast<std::size_t>(global_face)];
        const double left_centre =
            metric.centres().data[static_cast<std::size_t>(global_face - 1)];
        const double right_centre =
            metric.centres().data[static_cast<std::size_t>(global_face)];
        expected_flux =
            area * ((right_centre - face_location) * left +
                    (face_location - left_centre) * right) /
            (right_centre - left_centre);
        const auto cell_volume_for_normal = [&](std::int32_t normal) {
          double volume = 1.0;
          volume *= dependencies.geometry.x().widths().data[
              axis == CartesianAxis::x
                  ? static_cast<std::size_t>(normal)
                  : static_cast<std::size_t>(patch.begin.x)];
          volume *= dependencies.geometry.y().widths().data[
              axis == CartesianAxis::y
                  ? static_cast<std::size_t>(normal)
                  : static_cast<std::size_t>(patch.begin.y)];
          volume *= dependencies.geometry.z().widths().data[
              axis == CartesianAxis::z
                  ? static_cast<std::size_t>(normal)
                  : static_cast<std::size_t>(patch.begin.z)];
          return volume;
        };
        const double left_coefficient =
            cell_volume_for_normal(global_face - 1) / diagonals[component];
        const double right_coefficient =
            cell_volume_for_normal(global_face) / diagonals[component];
        const double left_distance = face_location - left_centre;
        const double right_distance = right_centre - face_location;
        expected_coefficient =
            area / (left_distance / left_coefficient +
                    right_distance / right_coefficient);
      }
      passed &= expect(
          std::abs(flux_face.unchecked(face) - expected_flux) <=
                  5.0e-13 * std::max(1.0, std::abs(expected_flux)) &&
              std::abs(coefficient.unchecked(face) - expected_coefficient) <=
                  5.0e-13 *
                      std::max(1.0, std::abs(expected_coefficient)),
          rank,
          physical
              ? "physical no-slip face has zero flux/Neumann coefficient"
              : "MPI partition face uses remote HbyA and exact transmissibility");
    }
  }
  return passed;
}

bool test_collective_compile_and_self_reference(MPI_Comm world, int rank,
                                                int size) {
  Dependencies distributed;
  bool passed = expect(make_dependencies(world, distributed), rank,
                       "non-divisible distributed dependencies compile");
  EquationPlanSet plan;
  const Status status = compile(world, distributed, plan_spec(), plan);
  passed &= expect(static_cast<bool>(status) && same_collective_status(world, status),
                   rank, "equation compile succeeds collectively");
  passed &= expect(plan.fingerprint() != 0U &&
                       plan.global_cells().x == 17 &&
                       plan.global_cells().y == 11 &&
                       plan.global_cells().z == 7,
                   rank, "distributed plan retains global identity");

  Dependencies reference;
  passed &= expect(make_dependencies(MPI_COMM_SELF, reference), rank,
                   "each rank builds an in-process full-domain reference");
  EquationPlanSet full;
  passed &= expect(static_cast<bool>(compile(MPI_COMM_SELF, reference,
                                             plan_spec(), full)),
                   rank, "full-domain reference equation plan compiles");
  passed &= expect(plan.semantic_fingerprint() == full.semantic_fingerprint(),
                   rank, "1/2/4 decomposition shares the COMM_SELF semantic plan");
  passed &= expect(plan.local_cells() ==
                       static_cast<std::size_t>(distributed.patch.cells.x) *
                           distributed.patch.cells.y * distributed.patch.cells.z &&
                       plan.local_cells() < 17U * 11U * 7U +
                                                (size == 1 ? 1U : 0U),
                   rank, "plan owns only its non-divisible local partition");
  return passed;
}

bool test_collective_atomic_failure(MPI_Comm world, int rank, int size) {
  Dependencies input;
  if (!make_dependencies(world, input)) {
    return false;
  }
  EquationPlanSet published;
  EquationPlanSpec good = plan_spec();
  if (!compile(world, input, good, published)) {
    return false;
  }
  const PlanFingerprint before = published.fingerprint();
  EquationPlanSpec divergent = good;
  if (size > 1 && rank == 1) {
    divergent.velocity = divergent.density;
  } else if (size == 1) {
    divergent.velocity = divergent.density;
  }
  EquationCompileDiagnostics diagnostics;
  const Status failure = compile(world, input, divergent, published, &diagnostics);
  bool passed = expect(failure.code == StatusCode::invalid_plan &&
                           same_collective_status(world, failure),
                       rank, "lowest-rank invalid spec fails collectively");
  passed &= expect(diagnostics.lowest_failing_rank == (size > 1 ? 1 : 0), rank,
                   "compile diagnoses the lowest failing rank");
  passed &= expect(published.fingerprint() == before, rank,
                   "collective compile failure is atomic");
  passed &= expect(static_cast<bool>(compile(world, input, good, published)), rank,
                   "collective compile retries cleanly after failure");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  bool passed = test_collective_compile_and_self_reference(MPI_COMM_WORLD, rank,
                                                           size);
  passed &= test_collective_atomic_failure(MPI_COMM_WORLD, rank, size);
  passed &= test_distributed_continuity_global_id_oracle(MPI_COMM_WORLD, rank);
  passed &= test_distributed_momentum_tensor_stretched_oracle(MPI_COMM_WORLD,
                                                              rank);
  passed &= test_conservative_momentum_predictor_limiter(MPI_COMM_WORLD,
                                                         rank);
  passed &= test_piso_intermediate_halo_oracle(MPI_COMM_WORLD, rank);
  MPI_Finalize();
  return passed ? 0 : 1;
}
