// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <vector>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void *allocate(std::size_t size) {
  observe();
  void *result = std::malloc(size == 0U ? 1U : size);
  if (result == nullptr) throw std::bad_alloc{};
  return result;
}

void *allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void *result = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&result, alignment, requested) != 0 || result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard &) = delete;
  Guard &operator=(const Guard &) = delete;
};

}  // namespace allocation_observer

void *operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void *operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}
void *operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void *operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void *pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

constexpr Int3 kCells{3, 2, 2};
constexpr std::int32_t kGhosts = 2;
constexpr double kSentinel = -777777.0;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool near(double actual, double expected, double tolerance = 1.0e-12) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             tolerance * std::max(1.0, std::abs(expected));
}

BoundaryThermophysicalGhostContext ghost_context(
    const BoundaryPlan& boundary,
    BoundaryThermophysicalGhostPhase phase =
        BoundaryThermophysicalGhostPhase::corrector_one) {
  return {8101U, boundary.geometry_fingerprint(), 8103U,
          boundary.revision(), phase};
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(Int3 cells, std::int32_t ghosts, FieldId id,
                      RevisionToken revision,
                      StorageIdentity storage_identity) {
  OwnedField owned;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  owned.storage.assign(nx * ny * nz, kSentinel);
  owned.view.base =
      owned.storage.data() + ghosts + ghosts * nx + ghosts * nx * ny;
  owned.view.interior = cells;
  owned.view.ghosts = {ghosts, ghosts, ghosts};
  owned.view.components = 1U;
  owned.view.stride_y = nx;
  owned.view.stride_z = nx * ny;
  owned.view.component_stride = nx * ny * nz;
  owned.view.field = id;
  owned.view.revision = revision;
  owned.view.storage_identity = storage_identity;
  owned.view.revision_domain = 7001U;
  return owned;
}

OwnedField make_field(FieldId id, RevisionToken revision,
                      StorageIdentity storage_identity) {
  return make_field(kCells, kGhosts, id, revision, storage_identity);
}

SpeciesThermophysicalSpec constant_air() {
  SpeciesThermophysicalSpec species;
  species.stable_name = "air";
  species.molecular_weight = 28.96546;
  species.temperature_switch = 1000.0;
  const double gas = kUniversalGasConstant / species.molecular_weight;
  species.nasa7_low[0U] = 1005.0 / gas;
  species.nasa7_high[0U] = 1005.0 / gas;
  species.viscosity_reference = 1.8e-5;
  species.conductivity = 0.026;
  return species;
}

ThermophysicalSpec constant_spec() {
  ThermophysicalSpec spec;
  spec.data_file = "bc-thermo-constant.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  spec.species.push_back(constant_air());
  return spec;
}

SpeciesThermophysicalSpec varying_species(std::string_view name, double mw,
                                          double a1, double a2,
                                          double viscosity,
                                          double conductivity) {
  SpeciesThermophysicalSpec species;
  species.stable_name.assign(name.data(), name.size());
  species.molecular_weight = mw;
  species.temperature_switch = 1000.0;
  species.nasa7_low = {a1, a2, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double high_a2 = 0.5 * a2;
  const double high_a1 = a1 + (a2 - high_a2) * species.temperature_switch;
  const double low_h_switch =
      a1 * species.temperature_switch +
      0.5 * a2 * species.temperature_switch * species.temperature_switch;
  const double high_h_without_offset =
      high_a1 * species.temperature_switch +
      0.5 * high_a2 * species.temperature_switch * species.temperature_switch;
  species.nasa7_high = {high_a1, high_a2, 0.0,
                        0.0,     0.0,     low_h_switch - high_h_without_offset,
                        0.0};
  species.viscosity_reference = viscosity;
  species.conductivity = conductivity;
  return species;
}

ThermophysicalSpec mixture_spec() {
  ThermophysicalSpec spec;
  spec.data_file = "bc-thermo-mixture.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  spec.species.push_back(
      varying_species("A", 24.0, 3.2, 4.0e-4, 1.7e-5, 0.024));
  spec.species.push_back(
      varying_species("B", 40.0, 4.1, 1.0e-4, 2.2e-5, 0.031));
  return spec;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {3.0, 2.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kCells;
  mesh.minimum_spacing = {1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 64U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 20U;
  return mesh;
}

CartesianMeshSpec mpi_mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {2.0, 8.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {2, 8, 2};
  mesh.minimum_spacing = {1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 64U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 20U;
  return mesh;
}

ValidatedModel boundary_model() {
  ValidatedModel model;
  model.fingerprint = UINT64_C(0x4243544845524d4f);
  model.pressure_reference = PressureReferenceKind::boundary_absolute;
  model.mesh = mesh_spec();
  model.schemes.momentum = ConvectionScheme::limited_central2;
  model.schemes.enthalpy = ConvectionScheme::limited_central2;
  model.schemes.species = ConvectionScheme::tvd2;
  model.schemes.passive_scalar = ConvectionScheme::tvd2;
  for (BoundaryFaceSpec &face : model.boundaries) {
    face.flow_kind = BoundaryKind::no_slip_wall;
    face.thermal_kind = BoundaryKind::adiabatic_wall;
    face.temperature = 500.0;
    face.mach_limit = 0.95;
  }
  BoundaryFaceSpec &inlet = model.boundaries[0U];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.thermal_kind = BoundaryKind::none;
  inlet.velocity = {1.0, 0.0, 0.0};
  inlet.temperature = 500.0;
  BoundaryFaceSpec &outlet = model.boundaries[1U];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.thermal_kind = BoundaryKind::none;
  outlet.pressure = 101325.0;
  return model;
}

ValidatedModel mixture_boundary_model() {
  ValidatedModel model = boundary_model();
  model.fingerprint = UINT64_C(0x42434d4958545552);
  model.transported_scalars.push_back(
      {"A", TransportedScalarRole::species, 1.0, 1.0});
  for (std::size_t index = 0U; index < model.boundaries.size(); ++index) {
    BoundaryFaceSpec &face = model.boundaries[index];
    const bool inlet = index == 0U;
    face.scalars.push_back({"A",
                            inlet ? ScalarBoundaryKind::dirichlet
                                  : ScalarBoundaryKind::zero_gradient,
                            inlet ? 0.35 : 0.0,
                            ScalarBoundaryKind::zero_gradient, 0.0});
  }
  return model;
}

ValidatedModel mpi_boundary_model() {
  ValidatedModel model = boundary_model();
  model.fingerprint = UINT64_C(0x42434d504947484f);
  model.mesh = mpi_mesh_spec();
  model.boundaries[4U].flow_kind = BoundaryKind::periodic;
  model.boundaries[4U].thermal_kind = BoundaryKind::none;
  model.boundaries[5U].flow_kind = BoundaryKind::periodic;
  model.boundaries[5U].thermal_kind = BoundaryKind::none;
  return model;
}

struct Fixture {
  BoundaryPlan boundary;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  OwnedField pressure = make_field(0U, 11U, 101U);
  OwnedField enthalpy = make_field(0U, 12U, 102U);
  std::array<OwnedField, 8U> output{
      make_field(20U, 21U, 120U), make_field(21U, 22U, 121U),
      make_field(22U, 23U, 122U), make_field(23U, 24U, 123U),
      make_field(24U, 25U, 124U), make_field(25U, 26U, 125U),
      make_field(26U, 27U, 126U), make_field(27U, 28U, 127U)};
  std::array<BoundaryGhostFieldAuthority, 2U> field_authority{};

  bool initialize() {
    CartesianGeometryPlan geometry;
    MeshPatch patch;
    FieldRegistry registry;
    SchemePlan schemes;
    TimeSchemePlan time;
    const ValidatedModel model = boundary_model();
    if (!CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh_spec(), {},
                                            geometry, patch) ||
        !BoundaryCompiler::compile(MPI_COMM_SELF, model, geometry, patch,
                                   registry, boundary, schemes, time)) {
      return false;
    }
    pressure.view.field = boundary.pressure_field();
    enthalpy.view.field = boundary.enthalpy_field();
    const ThermophysicalSpec thermo = constant_spec();
    if (!ThermodynamicsPlan::compile(thermo, {}, thermodynamics) ||
        !TransportPlan::compile(thermo, thermodynamics, transport)) {
      return false;
    }
    field_authority = {{
        make_boundary_ghost_field_authority(as_const(pressure.view)),
        make_boundary_ghost_field_authority(as_const(enthalpy.view)),
    }};
    return true;
  }

  BoundaryThermophysicalGhostInput input() const {
    return {100000.0,
            as_const(pressure.view),
            as_const(enthalpy.view),
            {},
            {9001U,
             boundary.revision(),
             boundary.local_layout_fingerprint(),
             kCells,
             {kGhosts, kGhosts, kGhosts},
             static_cast<std::uint8_t>(kGhosts),
             {field_authority.data(), field_authority.size()}}};
  }

  BoundaryThermophysicalGhostOutput outputs() const {
    return {output[0U].view, output[1U].view, output[2U].view, output[3U].view,
            output[4U].view, output[5U].view, output[6U].view, output[7U].view};
  }
};

struct MixtureFixture {
  BoundaryPlan boundary;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  OwnedField pressure = make_field(0U, 31U, 201U);
  OwnedField enthalpy = make_field(0U, 32U, 202U);
  OwnedField species = make_field(0U, 33U, 203U);
  std::array<OwnedField, 8U> output{
      make_field(30U, 41U, 220U), make_field(31U, 42U, 221U),
      make_field(32U, 43U, 222U), make_field(33U, 44U, 223U),
      make_field(34U, 45U, 224U), make_field(35U, 46U, 225U),
      make_field(36U, 47U, 226U), make_field(37U, 48U, 227U)};
  std::array<ConstFieldView, 1U> species_views{};
  std::array<BoundaryGhostFieldAuthority, 3U> field_authority{};

  bool initialize() {
    CartesianGeometryPlan geometry;
    MeshPatch patch;
    FieldRegistry registry;
    SchemePlan schemes;
    TimeSchemePlan time;
    const ValidatedModel model = mixture_boundary_model();
    const Status geometry_status = CartesianGeometryCompiler::compile(
        MPI_COMM_SELF, mesh_spec(), {}, geometry, patch);
    const Status boundary_status =
        geometry_status
            ? BoundaryCompiler::compile(MPI_COMM_SELF, model, geometry, patch,
                                        registry, boundary, schemes, time)
            : geometry_status;
    if (!geometry_status || !boundary_status) {
      std::cerr << "mixture geometry/boundary status="
                << static_cast<unsigned>(geometry_status.code) << "/"
                << geometry_status.detail << " "
                << static_cast<unsigned>(boundary_status.code) << "/"
                << boundary_status.detail << '\n';
      return false;
    }
    pressure.view.field = boundary.pressure_field();
    enthalpy.view.field = boundary.enthalpy_field();
    const Span<const BoundaryTransportedField> transported =
        boundary.transported_fields();
    if (transported.size != 1U ||
        transported.data[0U].role != BoundaryScalarRole::species) {
      return false;
    }
    species.view.field = transported.data[0U].field;
    species_views[0U] = as_const(species.view);
    const ThermophysicalSpec thermo = mixture_spec();
    const std::array<TransportedScalarSpec, 1U> catalog{{
        {"A", TransportedScalarRole::species, 1.0, 1.0},
    }};
    const Status thermo_status = ThermodynamicsPlan::compile(
        thermo, {catalog.data(), catalog.size()}, thermodynamics);
    const Status transport_status =
        thermo_status
            ? TransportPlan::compile(thermo, thermodynamics, transport)
            : thermo_status;
    if (!thermo_status || !transport_status) {
      std::cerr << "mixture thermo/transport status="
                << static_cast<unsigned>(thermo_status.code) << "/"
                << thermo_status.detail << " "
                << static_cast<unsigned>(transport_status.code) << "/"
                << transport_status.detail << '\n';
      return false;
    }
    field_authority = {{
        make_boundary_ghost_field_authority(as_const(pressure.view)),
        make_boundary_ghost_field_authority(as_const(enthalpy.view)),
        make_boundary_ghost_field_authority(as_const(species.view)),
    }};
    return true;
  }

  BoundaryThermophysicalGhostInput input() const {
    return {100000.0,
            as_const(pressure.view),
            as_const(enthalpy.view),
            {species_views.data(), species_views.size()},
            {9002U,
             boundary.revision(),
             boundary.local_layout_fingerprint(),
             kCells,
             {kGhosts, kGhosts, kGhosts},
             static_cast<std::uint8_t>(kGhosts),
             {field_authority.data(), field_authority.size()}}};
  }

  BoundaryThermophysicalGhostOutput outputs() const {
    return {output[0U].view, output[1U].view, output[2U].view, output[3U].view,
            output[4U].view, output[5U].view, output[6U].view, output[7U].view};
  }
};

bool test_constant_cp_physical_ghost_closure() {
  Fixture fixture;
  bool passed = expect(fixture.initialize(), "constant-cp fixture compiles");
  if (!passed) return false;

  constexpr double temperature = 500.0;
  constexpr double cp = 1005.0;
  constexpr double pi = 1325.0;
  const double h = cp * temperature;
  std::fill(fixture.pressure.storage.begin(), fixture.pressure.storage.end(),
            pi);
  std::fill(fixture.enthalpy.storage.begin(), fixture.enthalpy.storage.end(),
            h);
  for (OwnedField &field : fixture.output) {
    std::fill(field.storage.begin(), field.storage.end(), kSentinel);
    for (std::int32_t z = 0; z < kCells.z; ++z) {
      for (std::int32_t y = 0; y < kCells.y; ++y) {
        for (std::int32_t x = 0; x < kCells.x; ++x) {
          field.view.unchecked({x, y, z}, 0U) = 1000.0 + field.view.field;
        }
      }
    }
  }

  const std::array<BoundaryGhostFieldAuthority, 2U> legacy_authorities{{
      {fixture.pressure.view.field, fixture.pressure.view.revision,
       fixture.pressure.view.storage_identity,
       fixture.pressure.view.revision_domain},
      {fixture.enthalpy.view.field, fixture.enthalpy.view.revision,
       fixture.enthalpy.view.storage_identity,
       fixture.enthalpy.view.revision_domain},
  }};
  BoundaryThermophysicalGhostInput legacy_input = fixture.input();
  legacy_input.authority.fields = {legacy_authorities.data(),
                                   legacy_authorities.size()};
  Status legacy_status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      legacy_input, fixture.outputs());
  passed &= expect(static_cast<bool>(legacy_status),
                   "legacy four-field authority and close overload remain source compatible");

  Status status;
  BoundaryThermophysicalGhostCertificate certificate;
  std::size_t allocations = 99U;
  {
    allocation_observer::Guard guard;
    status = BoundaryThermophysicalFaceClosure::close(
        fixture.boundary, fixture.thermodynamics, fixture.transport,
        fixture.input(), fixture.outputs(), ghost_context(fixture.boundary),
        certificate);
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(status),
                   "constant-cp physical closure succeeds");
  const BoundaryThermophysicalGhostBinding certified_binding{
      100000.0, as_const(fixture.pressure.view),
      as_const(fixture.enthalpy.view), {},
      as_const(fixture.output[0U].view)};
  passed &= expect(
      certificate.valid() && certificate.matches(
                                 fixture.boundary,
                                 ghost_context(fixture.boundary),
                                 certified_binding),
      "closure publishes a current physical-ghost certificate");
  passed &= expect(allocations == 0U, "closure hot call allocates no memory");

  const Int3 mutated_ghost{-1, 0, 0};
  const double pressure_before =
      fixture.pressure.view.unchecked(mutated_ghost, 0U);
  fixture.pressure.view.unchecked(mutated_ghost, 0U) += 1.0;
  passed &= expect(!certificate.matches(fixture.boundary,
                                         ghost_context(fixture.boundary),
                                         certified_binding),
                   "same-revision physical pressure mutation invalidates");
  fixture.pressure.view.unchecked(mutated_ghost, 0U) = pressure_before;
  const double density_before =
      fixture.output[0U].view.unchecked(mutated_ghost, 0U);
  fixture.output[0U].view.unchecked(mutated_ghost, 0U) *= 1.01;
  passed &= expect(!certificate.matches(fixture.boundary,
                                         ghost_context(fixture.boundary),
                                         certified_binding),
                   "same-revision physical density mutation invalidates");
  fixture.output[0U].view.unchecked(mutated_ghost, 0U) = density_before;
  BoundaryThermophysicalGhostContext wrong_phase =
      ghost_context(fixture.boundary);
  wrong_phase.phase = BoundaryThermophysicalGhostPhase::corrector_two;
  passed &= expect(!certificate.matches(fixture.boundary, wrong_phase,
                                         certified_binding),
                   "wrong C1/C2 phase invalidates");
  BoundaryThermophysicalGhostContext wrong_time =
      ghost_context(fixture.boundary);
  ++wrong_time.target_time;
  passed &= expect(!certificate.matches(fixture.boundary, wrong_time,
                                         certified_binding),
                   "wrong target time invalidates");
  BoundaryThermophysicalGhostCertificate next_time_certificate;
  const Status next_time_status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), fixture.outputs(), wrong_time, next_time_certificate);
  passed &= expect(
      static_cast<bool>(next_time_status) && next_time_certificate.valid() &&
          next_time_certificate.collective_lineage() !=
              certificate.collective_lineage(),
      "target-time change produces a different collective lineage");

  ThermoState expected_thermo;
  MolecularTransportState expected_transport;
  passed &=
      expect(static_cast<bool>(fixture.thermodynamics.evaluate(
                 101325.0, h, {}, Real3{400.0, -50.0, 7.0}, expected_thermo)) &&
                 static_cast<bool>(fixture.transport.evaluate(
                     temperature, {}, expected_transport)),
             "independent thermophysical oracle evaluates");
  const std::array<double, 8U> expected{
      expected_thermo.rho,
      expected_thermo.temperature,
      expected_thermo.cp,
      expected_thermo.drho_dp_hY,
      expected_thermo.drho_dh_pY,
      expected_transport.viscosity,
      expected_transport.conductivity,
      expected_transport.conductivity / expected_thermo.cp};
  const std::array<Int3, 4U> probes{{
      {-1, 0, 0},
      {kCells.x, 1, 1},
      {-2, 1, 0},
      {-2, -2, -2},
  }};
  for (std::size_t field = 0U; field < fixture.output.size(); ++field) {
    for (const Int3 probe : probes) {
      passed &=
          expect(near(fixture.output[field].view.unchecked(probe, 0U),
                      expected[field]),
                 "low/high, multi-layer, and corner ghosts use EOS/transport");
    }
    passed &= expect(fixture.output[field].view.unchecked({1, 1, 1}, 0U) ==
                         1000.0 + fixture.output[field].view.field,
                     "closure leaves interior unchanged");
  }
  return passed;
}

bool test_nasa_mixture_uses_fixed_ph_y() {
  MixtureFixture fixture;
  bool passed = expect(fixture.initialize(), "NASA mixture fixture compiles");
  if (!passed) return false;

  constexpr double temperature = 650.0;
  constexpr double fraction = 0.35;
  constexpr double perturbation = 4321.0;
  const std::array<double, 1U> fractions{fraction};
  double enthalpy = 0.0;
  double cp = 0.0;
  double gas = 0.0;
  passed &= expect(static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
                       temperature, {fractions.data(), fractions.size()},
                       enthalpy, cp, gas)),
                   "NASA mixture target enthalpy evaluates");
  std::fill(fixture.pressure.storage.begin(), fixture.pressure.storage.end(),
            perturbation);
  std::fill(fixture.enthalpy.storage.begin(), fixture.enthalpy.storage.end(),
            enthalpy);
  std::fill(fixture.species.storage.begin(), fixture.species.storage.end(),
            fraction);
  for (OwnedField &field : fixture.output) {
    std::fill(field.storage.begin(), field.storage.end(), kSentinel);
  }

  BoundaryThermophysicalGhostCertificate certificate;
  const Status status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), fixture.outputs(), ghost_context(fixture.boundary),
      certificate);
  passed &= expect(static_cast<bool>(status),
                   "NASA mixture physical closure succeeds");
  const BoundaryThermophysicalGhostBinding certified_binding{
      100000.0, as_const(fixture.pressure.view),
      as_const(fixture.enthalpy.view),
      {fixture.species_views.data(), fixture.species_views.size()},
      as_const(fixture.output[0U].view)};
  const Int3 mutated_species_ghost{-1, 0, 0};
  const double species_before =
      fixture.species.view.unchecked(mutated_species_ghost, 0U);
  fixture.species.view.unchecked(mutated_species_ghost, 0U) += 0.01;
  passed &= expect(
      !certificate.matches(fixture.boundary, ghost_context(fixture.boundary),
                           certified_binding),
      "same-revision physical species mutation invalidates");
  fixture.species.view.unchecked(mutated_species_ghost, 0U) = species_before;

  ThermoState expected_thermo;
  MolecularTransportState expected_transport;
  passed &= expect(
      static_cast<bool>(fixture.thermodynamics.evaluate(
          104321.0, enthalpy, {fractions.data(), fractions.size()},
          Real3{-900.0, 700.0, 50.0}, expected_thermo)) &&
          static_cast<bool>(fixture.transport.evaluate(
              temperature, {fractions.data(), fractions.size()},
              expected_transport)),
      "mixture oracle is independent of arbitrary velocity for closure fields");
  const std::array<double, 8U> expected{
      expected_thermo.rho,
      expected_thermo.temperature,
      expected_thermo.cp,
      expected_thermo.drho_dp_hY,
      expected_thermo.drho_dh_pY,
      expected_transport.viscosity,
      expected_transport.conductivity,
      expected_transport.conductivity / expected_thermo.cp};
  const Int3 probe{kCells.x + 1, kCells.y + 1, kCells.z + 1};
  for (std::size_t field = 0U; field < fixture.output.size(); ++field) {
    passed &= expect(near(fixture.output[field].view.unchecked(probe, 0U),
                          expected[field], 2.0e-11),
                     "mixture corner ghost closes from the same p/h/Y state");
  }
  return passed;
}

bool test_invalid_state_and_contract_are_atomic() {
  MixtureFixture fixture;
  bool passed =
      expect(fixture.initialize(), "atomic-rejection mixture fixture compiles");
  if (!passed) return false;

  constexpr double temperature = 650.0;
  constexpr double fraction = 0.35;
  const std::array<double, 1U> fractions{fraction};
  double enthalpy = 0.0;
  double cp = 0.0;
  double gas = 0.0;
  passed &= expect(static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
                       temperature, {fractions.data(), fractions.size()},
                       enthalpy, cp, gas)),
                   "atomic-rejection target enthalpy evaluates");
  if (!passed) return false;

  const auto reset_valid_inputs = [&]() {
    std::fill(fixture.pressure.storage.begin(), fixture.pressure.storage.end(),
              1325.0);
    std::fill(fixture.enthalpy.storage.begin(), fixture.enthalpy.storage.end(),
              enthalpy);
    std::fill(fixture.species.storage.begin(), fixture.species.storage.end(),
              fraction);
  };
  const auto reset_outputs = [&]() {
    for (std::size_t index = 0U; index < fixture.output.size(); ++index) {
      std::fill(fixture.output[index].storage.begin(),
                fixture.output[index].storage.end(),
                kSentinel - static_cast<double>(index));
    }
  };
  const auto output_snapshot = [&]() {
    std::array<std::vector<double>, 8U> snapshot;
    for (std::size_t index = 0U; index < snapshot.size(); ++index) {
      snapshot[index] = fixture.output[index].storage;
    }
    return snapshot;
  };
  const auto outputs_unchanged =
      [&](const std::array<std::vector<double>, 8U> &snapshot) {
        for (std::size_t index = 0U; index < snapshot.size(); ++index) {
          if (snapshot[index] != fixture.output[index].storage) return false;
        }
        return true;
      };
  constexpr Int3 last_physical_corner{
      kCells.x + kGhosts - 1, kCells.y + kGhosts - 1, kCells.z + kGhosts - 1};

  reset_valid_inputs();
  reset_outputs();
  fixture.pressure.view.unchecked(last_physical_corner, 0U) = -100000.0;
  auto snapshot = output_snapshot();
  BoundaryThermophysicalGhostCertificate certificate;
  Status status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), fixture.outputs(), ghost_context(fixture.boundary),
      certificate);
  passed &= expect(!status && status.code == StatusCode::numerical_failure &&
                       outputs_unchanged(snapshot) && !certificate.valid(),
                   "nonpositive p_abs rejects before any output mutation");

  reset_valid_inputs();
  reset_outputs();
  fixture.enthalpy.view.unchecked(last_physical_corner, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), fixture.outputs(), ghost_context(fixture.boundary),
      certificate);
  passed &= expect(!status && status.code == StatusCode::numerical_failure &&
                       outputs_unchanged(snapshot) && !certificate.valid(),
                   "nonfinite h rejects before any output mutation");

  reset_valid_inputs();
  reset_outputs();
  fixture.species.view.unchecked(last_physical_corner, 0U) = 1.2;
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), fixture.outputs(), ghost_context(fixture.boundary),
      certificate);
  passed &= expect(!status && status.code == StatusCode::numerical_failure &&
                       outputs_unchanged(snapshot) && !certificate.valid(),
                   "out-of-simplex Y rejects before any output mutation");

  reset_valid_inputs();
  reset_outputs();
  std::array<BoundaryGhostFieldAuthority, 3U> stale = fixture.field_authority;
  ++stale[0U].revision;
  BoundaryThermophysicalGhostInput stale_input = fixture.input();
  stale_input.authority.fields = {stale.data(), stale.size()};
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport, stale_input,
      fixture.outputs(), ghost_context(fixture.boundary), certificate);
  passed &= expect(!status && status.code == StatusCode::invalid_plan &&
                       outputs_unchanged(snapshot) && !certificate.valid(),
                   "stale ghost revision rejects atomically");

  reset_outputs();
  BoundaryThermophysicalGhostOutput wrong_shape = fixture.outputs();
  wrong_shape.temperature.ghosts.x = 1;
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), wrong_shape, ghost_context(fixture.boundary),
      certificate);
  passed &= expect(!status && status.code == StatusCode::invalid_plan &&
                       outputs_unchanged(snapshot),
                   "mismatched output shape rejects atomically");

  reset_outputs();
  const std::vector<double> enthalpy_snapshot = fixture.enthalpy.storage;
  BoundaryThermophysicalGhostOutput aliased = fixture.outputs();
  aliased.density = fixture.enthalpy.view;
  aliased.density.field = fixture.output[0U].view.field;
  aliased.density.revision = fixture.output[0U].view.revision;
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      fixture.input(), aliased, ghost_context(fixture.boundary), certificate);
  passed &= expect(!status && status.code == StatusCode::invalid_plan &&
                       outputs_unchanged(snapshot) &&
                       fixture.enthalpy.storage == enthalpy_snapshot,
                   "input/output alias rejects without mutating either side");

  reset_valid_inputs();
  reset_outputs();
  OwnedField foreign_pressure = make_field(
      fixture.pressure.view.field, fixture.pressure.view.revision,
      fixture.pressure.view.storage_identity);
  foreign_pressure.view.replica = fixture.pressure.view.replica;
  std::fill(foreign_pressure.storage.begin(), foreign_pressure.storage.end(),
            1325.0);
  BoundaryThermophysicalGhostInput foreign_input = fixture.input();
  foreign_input.pressure_perturbation = as_const(foreign_pressure.view);
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      foreign_input, fixture.outputs(), ghost_context(fixture.boundary),
      certificate);
  passed &= expect(!status && status.code == StatusCode::invalid_plan &&
                       outputs_unchanged(snapshot),
                   "same-revision foreign base rejects atomically");

  BoundaryThermophysicalGhostInput foreign_replica_input = fixture.input();
  ConstFieldView foreign_replica = as_const(fixture.pressure.view);
  ++foreign_replica.replica;
  foreign_replica_input.pressure_perturbation = foreign_replica;
  snapshot = output_snapshot();
  status = BoundaryThermophysicalFaceClosure::close(
      fixture.boundary, fixture.thermodynamics, fixture.transport,
      foreign_replica_input, fixture.outputs(), ghost_context(fixture.boundary),
      certificate);
  passed &= expect(!status && status.code == StatusCode::invalid_plan &&
                       outputs_unchanged(snapshot),
                   "same-base foreign replica rejects atomically");
  return passed;
}

bool test_mpi_and_periodic_ghosts_remain_halo_owned(int rank, int size) {
  if (size != 2) return true;

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  FieldRegistry registry;
  SchemePlan schemes;
  TimeSchemePlan time;
  const ValidatedModel model = mpi_boundary_model();
  bool passed =
      expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                 MPI_COMM_WORLD, mpi_mesh_spec(), {}, geometry, patch)) &&
                 static_cast<bool>(BoundaryCompiler::compile(
                     MPI_COMM_WORLD, model, geometry, patch, registry, boundary,
                     schemes, time)),
             "2-rank MPI/periodic boundary fixture compiles");
  passed &= expect(patch.process_grid.x == 1 && patch.process_grid.y == 2 &&
                       patch.process_grid.z == 1,
                   "2-rank fixture decomposes along y");
  if (!passed) return false;

  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  const ThermophysicalSpec thermo = constant_spec();
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                       thermo, {}, thermodynamics)) &&
                       static_cast<bool>(TransportPlan::compile(
                           thermo, thermodynamics, transport)),
                   "2-rank constant thermophysics compiles");
  if (!passed) return false;

  constexpr std::int32_t ghosts = 2;
  OwnedField pressure =
      make_field(patch.cells, ghosts, boundary.pressure_field(), 61U,
                 301U + static_cast<unsigned>(rank));
  OwnedField enthalpy =
      make_field(patch.cells, ghosts, boundary.enthalpy_field(), 62U,
                 311U + static_cast<unsigned>(rank));
  std::array<OwnedField, 8U> outputs;
  for (std::size_t index = 0U; index < outputs.size(); ++index) {
    outputs[index] =
        make_field(patch.cells, ghosts, static_cast<FieldId>(40U + index),
                   static_cast<RevisionToken>(70U + index),
                   static_cast<StorageIdentity>(400U + 16U * rank + index));
  }
  constexpr double temperature = 500.0;
  constexpr double cp = 1005.0;
  constexpr double pi = 1325.0;
  const double h = cp * temperature;
  std::fill(pressure.storage.begin(), pressure.storage.end(), pi);
  std::fill(enthalpy.storage.begin(), enthalpy.storage.end(), h);
  for (OwnedField &output : outputs) {
    std::fill(output.storage.begin(), output.storage.end(), kSentinel);
  }
  const std::array<BoundaryGhostFieldAuthority, 2U> field_authority{{
      make_boundary_ghost_field_authority(as_const(pressure.view)),
      make_boundary_ghost_field_authority(as_const(enthalpy.view)),
  }};
  const BoundaryThermophysicalGhostInput input{
      100000.0,
      as_const(pressure.view),
      as_const(enthalpy.view),
      {},
      {static_cast<std::uintptr_t>(9100U + rank),
       boundary.revision(),
       boundary.local_layout_fingerprint(),
       patch.cells,
       {ghosts, ghosts, ghosts},
       static_cast<std::uint8_t>(ghosts),
       {field_authority.data(), field_authority.size()}}};
  const BoundaryThermophysicalGhostOutput output{
      outputs[0U].view, outputs[1U].view, outputs[2U].view, outputs[3U].view,
      outputs[4U].view, outputs[5U].view, outputs[6U].view, outputs[7U].view};
  BoundaryThermophysicalGhostCertificate certificate;
  const Status status = BoundaryThermophysicalFaceClosure::close(
      boundary, thermodynamics, transport, input, output,
      ghost_context(boundary), certificate);
  passed &=
      expect(static_cast<bool>(status) && certificate.valid(),
             "2-rank physical closure succeeds with a certificate");
  const std::array<std::uint64_t, 5U> local_tokens{
      certificate.collective_semantics(), certificate.collective_lineage(),
      certificate.collective_target(), certificate.rank_local_lineage(),
      certificate.rank_local_binding()};
  std::array<std::uint64_t, 5U> minimum_tokens{};
  std::array<std::uint64_t, 5U> maximum_tokens{};
  passed &= expect(
      MPI_Allreduce(local_tokens.data(), minimum_tokens.data(),
                    static_cast<int>(local_tokens.size()), MPI_UINT64_T,
                    MPI_MIN, MPI_COMM_WORLD) == MPI_SUCCESS &&
          MPI_Allreduce(local_tokens.data(), maximum_tokens.data(),
                        static_cast<int>(local_tokens.size()), MPI_UINT64_T,
                        MPI_MAX, MPI_COMM_WORLD) == MPI_SUCCESS &&
          minimum_tokens[0U] == maximum_tokens[0U] &&
          minimum_tokens[1U] == maximum_tokens[1U] &&
          minimum_tokens[2U] == maximum_tokens[2U] &&
          minimum_tokens[3U] != maximum_tokens[3U] &&
          minimum_tokens[4U] != maximum_tokens[4U],
      "collective semantics/lineage/target are rank invariant while stable "
      "storage lineage and concrete binding stay rank local");

  ThermoState expected;
  passed &= expect(
      static_cast<bool>(thermodynamics.evaluate(101325.0, h, {}, {}, expected)),
      "2-rank thermodynamic oracle evaluates");
  const std::int32_t physical_y = rank == 0 ? -1 : patch.cells.y;
  const std::int32_t mpi_y = rank == 0 ? patch.cells.y : -1;
  passed &=
      expect(near(outputs[0U].view.unchecked({-1, 0, 0}, 0U), expected.rho),
             "locally owned x physical ghost is closed");
  passed &= expect(
      near(outputs[0U].view.unchecked({0, physical_y, 0}, 0U), expected.rho),
      "locally owned outer y ghost is closed");
  passed &= expect(outputs[0U].view.unchecked({0, mpi_y, 0}, 0U) == kSentinel,
                   "internal y MPI ghost remains halo-owned");
  passed &= expect(outputs[0U].view.unchecked({0, 0, -1}, 0U) == kSentinel,
                   "periodic z ghost remains halo-owned");
  passed &= expect(outputs[0U].view.unchecked({-1, mpi_y, 0}, 0U) == kSentinel,
                   "x/MPI corner remains halo-owned");
  passed &= expect(
      near(outputs[0U].view.unchecked({-1, physical_y, 0}, 0U), expected.rho),
      "fully physical x/y corner is closed");

  const BoundaryThermophysicalGhostBinding certified_binding{
      100000.0, as_const(pressure.view), as_const(enthalpy.view), {},
      as_const(outputs[0U].view)};
  const Int3 periodic_ghost{0, 0, -1};
  const Int3 mpi_ghost{0, mpi_y, 0};
  pressure.view.unchecked(periodic_ghost, 0U) += 17.0;
  outputs[0U].view.unchecked(periodic_ghost, 0U) += 19.0;
  pressure.view.unchecked(mpi_ghost, 0U) += 23.0;
  outputs[0U].view.unchecked(mpi_ghost, 0U) += 29.0;
  passed &= expect(
      certificate.matches(boundary, ghost_context(boundary),
                          certified_binding),
      "post-issuance MPI and periodic ghost mutation is outside physical-ghost authority");

  const int local = passed ? 1 : 0;
  int global = 0;
  return MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ==
             MPI_SUCCESS &&
         global != 0;
}

}  // namespace

int main(int argc, char **argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
    MPI_Finalize();
    return 2;
  }
  const bool passed =
      test_constant_cp_physical_ghost_closure() &&
      test_nasa_mixture_uses_fixed_ph_y() &&
      test_invalid_state_and_contract_are_atomic() &&
      test_mpi_and_periodic_ghosts_remain_halo_owned(rank, size);
  const int finalized = MPI_Finalize();
  return passed && finalized == MPI_SUCCESS ? 0 : 1;
}
