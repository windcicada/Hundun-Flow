// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_physics.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  void* result = std::malloc(size == 0U ? 1U : size);
  if (result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* result = nullptr;
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

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool close(double actual, double expected, double relative = 1.0e-12,
           double absolute = 1.0e-14) {
  return std::abs(actual - expected) <=
         absolute + relative * std::max(std::abs(actual), std::abs(expected));
}

SpeciesThermophysicalSpec species(std::string_view name, double molecular_weight,
                                  double cp_over_r, double viscosity,
                                  double conductivity) {
  SpeciesThermophysicalSpec result;
  result.stable_name = std::string{name};
  result.molecular_weight = molecular_weight;
  result.temperature_switch = 1000.0;
  result.nasa7_low[0U] = cp_over_r;
  result.nasa7_high[0U] = cp_over_r;
  result.transport_law = TransportLaw::constant;
  result.viscosity_reference = viscosity;
  result.conductivity = conductivity;
  return result;
}

ThermophysicalSpec two_species_spec(bool sutherland) {
  ThermophysicalSpec spec;
  spec.data_file = "thermophysics.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 3000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 50U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 20U;
  spec.maximum_closed_mass_relative_step = 0.2;
  spec.species.push_back(species("light", 20.0, 3.5, 1.0e-5, 0.020));
  spec.species.push_back(species("heavy", 40.0, 4.0, 4.0e-5, 0.080));
  if (sutherland) {
    for (std::size_t index = 0U; index < spec.species.size(); ++index) {
      auto& item = spec.species[index];
      item.transport_law = TransportLaw::sutherland;
      item.transport_reference_temperature = 300.0;
      item.sutherland_temperature = 100.0 + 10.0 * index;
      item.prandtl = 0.70 + 0.05 * index;
      item.conductivity = 0.0;
    }
  }
  return spec;
}

std::array<TransportedScalarSpec, 1U> scalar_catalog() {
  return {TransportedScalarSpec{"light", TransportedScalarRole::species}};
}

bool compile_transport(const ThermophysicalSpec& spec,
                       ThermodynamicsPlan& thermodynamics,
                       TransportPlan& transport) {
  const auto scalars = scalar_catalog();
  return static_cast<bool>(ThermodynamicsPlan::compile(
             spec, {scalars.data(), scalars.size()}, thermodynamics)) &&
         static_cast<bool>(TransportPlan::compile(spec, thermodynamics,
                                                  transport));
}

double wilke_phi(double mu_i, double mu_j, double weight_i,
                 double weight_j) {
  const double numerator =
      1.0 + std::sqrt(mu_i / mu_j) * std::pow(weight_j / weight_i, 0.25);
  return numerator * numerator /
         std::sqrt(8.0 * (1.0 + weight_i / weight_j));
}

double sutherland_viscosity(const SpeciesThermophysicalSpec& item,
                            double temperature) {
  const double ratio =
      temperature / item.transport_reference_temperature;
  return item.viscosity_reference * ratio * std::sqrt(ratio) *
         (item.transport_reference_temperature +
          item.sutherland_temperature) /
         (temperature + item.sutherland_temperature);
}

MolecularTransportState independent_wilke_oracle(
    const ThermophysicalSpec& spec, double temperature,
    double first_mass_fraction) {
  const double second_mass_fraction = 1.0 - first_mass_fraction;
  const double first_moles =
      first_mass_fraction / spec.species[0U].molecular_weight;
  const double second_moles =
      second_mass_fraction / spec.species[1U].molecular_weight;
  const double first_mole_fraction =
      first_moles / (first_moles + second_moles);
  const double second_mole_fraction = 1.0 - first_mole_fraction;

  std::array<double, 2U> viscosity{};
  std::array<double, 2U> conductivity{};
  for (std::size_t index = 0U; index < spec.species.size(); ++index) {
    const auto& item = spec.species[index];
    if (item.transport_law == TransportLaw::constant) {
      viscosity[index] = item.viscosity_reference;
      conductivity[index] = item.conductivity;
    } else {
      viscosity[index] = sutherland_viscosity(item, temperature);
      const double cp = kUniversalGasConstant * item.nasa7_low[0U] /
                        item.molecular_weight;
      conductivity[index] = viscosity[index] * cp / item.prandtl;
    }
  }

  const double phi01 =
      wilke_phi(viscosity[0U], viscosity[1U],
                spec.species[0U].molecular_weight,
                spec.species[1U].molecular_weight);
  const double phi10 =
      wilke_phi(viscosity[1U], viscosity[0U],
                spec.species[1U].molecular_weight,
                spec.species[0U].molecular_weight);
  return {
      first_mole_fraction * viscosity[0U] /
              (first_mole_fraction + second_mole_fraction * phi01) +
          second_mole_fraction * viscosity[1U] /
              (first_mole_fraction * phi10 + second_mole_fraction),
      first_mole_fraction * conductivity[0U] /
              (first_mole_fraction + second_mole_fraction * phi01) +
          second_mole_fraction * conductivity[1U] /
              (first_mole_fraction * phi10 + second_mole_fraction),
  };
}

bool test_constant_and_wilke_transport() {
  bool passed = true;
  ThermophysicalSpec constant_spec = two_species_spec(false);
  ThermodynamicsPlan constant_thermo;
  TransportPlan constant;
  passed &= expect(compile_transport(constant_spec, constant_thermo, constant),
                   "constant transport compiles through thermo authority");
  passed &= expect(constant.kernel() == TransportKernel::constant &&
                       constant.fingerprint() != 0U,
                   "all-constant data select the constant kernel");

  const std::array independent{0.25};
  MolecularTransportState state;
  passed &= expect(static_cast<bool>(constant.evaluate(
                       400.0, {independent.data(), independent.size()}, state)),
                   "constant mixture evaluates");
  const MolecularTransportState expected = independent_wilke_oracle(
      constant_spec, 400.0, independent[0U]);
  passed &= expect(close(state.viscosity, expected.viscosity) &&
                       close(state.conductivity, expected.conductivity),
                   "constant mixture matches an independent Wilke oracle");

  MolecularTransportState repeated;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  Status repeated_status;
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < 1000U; ++iteration) {
      repeated_status = constant.evaluate(
          400.0, {independent.data(), independent.size()}, repeated);
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(repeated_status) &&
                       hot_allocations == 0U &&
                       repeated.viscosity == state.viscosity &&
                       repeated.conductivity == state.conductivity,
                   "repeated constant transport is deterministic and allocation-free");

  const std::array invalid_sum{1.01};
  repeated = {7.25, 8.5};
  const Status invalid = constant.evaluate(
      400.0, {invalid_sum.data(), invalid_sum.size()}, repeated);
  passed &= expect(invalid.code == StatusCode::numerical_failure &&
                       repeated.viscosity == 7.25 &&
                       repeated.conductivity == 8.5,
                   "failed composition leaves the caller output untouched");
  const Status invalid_temperature = constant.evaluate(
      199.0, {independent.data(), independent.size()}, repeated);
  passed &= expect(invalid_temperature.code == StatusCode::numerical_failure &&
                       repeated.viscosity == 7.25 &&
                       repeated.conductivity == 8.5,
                   "failed temperature leaves the caller output untouched");
  return passed;
}

bool test_sutherland_temperature_dependence() {
  bool passed = true;
  ThermophysicalSpec spec = two_species_spec(true);
  auto& constant_species = spec.species[1U];
  constant_species.transport_law = TransportLaw::constant;
  constant_species.transport_reference_temperature = 0.0;
  constant_species.sutherland_temperature = 0.0;
  constant_species.prandtl = 0.0;
  constant_species.conductivity = 0.080;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  passed &= expect(compile_transport(spec, thermodynamics, transport),
                   "Sutherland transport compiles");
  passed &= expect(transport.kernel() == TransportKernel::sutherland_wilke,
                   "mixed Sutherland/constant input selects its static kernel");
  const std::array independent{0.40};
  MolecularTransportState cold;
  MolecularTransportState hot;
  passed &= expect(static_cast<bool>(transport.evaluate(
                       300.0, {independent.data(), independent.size()}, cold)) &&
                       static_cast<bool>(transport.evaluate(
                           900.0, {independent.data(), independent.size()}, hot)),
                   "temperature-dependent mixture evaluates at two temperatures");
  passed &= expect(hot.viscosity > cold.viscosity &&
                       hot.conductivity > cold.conductivity &&
                       std::isfinite(hot.viscosity) &&
                       std::isfinite(hot.conductivity),
                   "Sutherland viscosity and cp/Pr conductivity grow smoothly");
  const MolecularTransportState expected_cold =
      independent_wilke_oracle(spec, 300.0, independent[0U]);
  const MolecularTransportState expected_hot =
      independent_wilke_oracle(spec, 900.0, independent[0U]);
  passed &= expect(close(cold.viscosity, expected_cold.viscosity) &&
                       close(cold.conductivity,
                             expected_cold.conductivity) &&
                       close(hot.viscosity, expected_hot.viscosity) &&
                       close(hot.conductivity, expected_hot.conductivity),
                   "mixed Sutherland transport matches independent temperature oracles");

  MolecularTransportState repeated;
  Status repeated_status;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < 1000U; ++iteration) {
      repeated_status = transport.evaluate(
          675.0, {independent.data(), independent.size()}, repeated);
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const MolecularTransportState expected_repeated =
      independent_wilke_oracle(spec, 675.0, independent[0U]);
  passed &= expect(static_cast<bool>(repeated_status) &&
                       hot_allocations == 0U &&
                       close(repeated.viscosity,
                             expected_repeated.viscosity) &&
                       close(repeated.conductivity,
                             expected_repeated.conductivity),
                   "Sutherland hot path is allocation-free and matches its oracle");
  return passed;
}

bool test_transport_rejects_mismatched_thermodynamics() {
  bool passed = true;
  ThermophysicalSpec thermodynamic_spec = two_species_spec(true);
  const auto scalars = scalar_catalog();
  ThermodynamicsPlan thermodynamics;
  passed &= expect(static_cast<bool>(ThermodynamicsPlan::compile(
                       thermodynamic_spec,
                       {scalars.data(), scalars.size()}, thermodynamics)),
                   "reference thermodynamics compiles for mismatch tests");

  ThermophysicalSpec changed_nasa = thermodynamic_spec;
  changed_nasa.species[0U].nasa7_low[0U] += 0.25;
  changed_nasa.species[0U].nasa7_high[0U] += 0.25;
  TransportPlan rejected;
  passed &= expect(!TransportPlan::compile(changed_nasa, thermodynamics,
                                           rejected),
                   "transport rejects NASA data from another thermo authority");

  ThermophysicalSpec reordered = thermodynamic_spec;
  std::swap(reordered.species[0U], reordered.species[1U]);
  passed &= expect(!TransportPlan::compile(reordered, thermodynamics, rejected),
                   "transport rejects reordered species from another authority");
  return passed;
}

CartesianMeshSpec uniform_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {-1.0, -1.0, -1.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {6, 5, 4};
  mesh.minimum_spacing = {1.0e-6, 1.0e-6, 1.0e-6};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {1000000U, 1U << 28U};
  return mesh;
}

CartesianMeshSpec stretched_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::tensor_stretched;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {4.0, 2.0, 1.0};
  mesh.has_base_spacing = true;
  mesh.base_spacing = {0.25, 0.25, 0.25};
  mesh.minimum_spacing = {0.025, 0.025, 0.025};
  mesh.max_growth_ratio = 1.25;
  mesh.focus_regions.push_back(
      {{1.0, 0.5, 0.25}, {3.0, 1.5, 0.75}, {0.1, 0.1, 0.1}});
  mesh.limits = {100000000U, 1U << 30U};
  return mesh;
}

struct GradientFixture {
  FieldId velocity{};
  FieldId pressure{};
  FieldSchema schema;
  StateLayers layers;
};

bool make_gradient_fixture(Int3 cells, GradientFixture& fixture) {
  FieldRegistry registry;
  if (!registry.declare_field("velocity", 3U, 1U, fixture.velocity) ||
      !registry.declare_field("pressure", 1U, 1U, fixture.pressure) ||
      !registry.freeze(fixture.schema)) {
    return false;
  }
  const std::array requests{
      ArenaFieldRequest{fixture.velocity, cells, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{fixture.pressure, cells, FieldPlacement{0U},
                        FieldLifetime::state_layer},
  };
  ArenaLayout layout;
  return ArenaLayout::compile(
             fixture.schema, {requests.data(), requests.size()}, layout) &&
         StateLayers::allocate(layout, fixture.layers);
}

void fill_linear_velocity(FieldView velocity, const CartesianGeometryPlan& geometry,
                          const MeshPatch& patch) {
  for (std::int32_t z = -1; z <= patch.cells.z; ++z) {
    for (std::int32_t y = -1; y <= patch.cells.y; ++y) {
      for (std::int32_t x = -1; x <= patch.cells.x; ++x) {
        const double xc = geometry.x().uniform_width() *
                          static_cast<double>(patch.begin.x + x);
        const double yc = geometry.y().uniform_width() *
                          static_cast<double>(patch.begin.y + y);
        const double zc = geometry.z().uniform_width() *
                          static_cast<double>(patch.begin.z + z);
        const Int3 index{x, y, z};
        velocity.unchecked(index, 0U) = 1.0 + 2.0 * xc + 3.0 * yc + 5.0 * zc;
        velocity.unchecked(index, 1U) = -2.0 + 7.0 * xc - 11.0 * yc + 13.0 * zc;
        velocity.unchecked(index, 2U) = 4.0 - 17.0 * xc + 19.0 * yc - 23.0 * zc;
      }
    }
  }
}

double extended_centre(const AxisMetrics& axis,
                       std::int32_t global_index) {
  const Span<const double> centres = axis.centres();
  const Span<const double> widths = axis.widths();
  if (global_index < 0) {
    return centres.data[0U] +
           static_cast<double>(global_index) * widths.data[0U];
  }
  if (static_cast<std::size_t>(global_index) >= centres.size) {
    return centres.data[centres.size - 1U] +
           static_cast<double>(global_index -
                               static_cast<std::int32_t>(centres.size - 1U)) *
               widths.data[widths.size - 1U];
  }
  return centres.data[static_cast<std::size_t>(global_index)];
}

void fill_quadratic_velocity(FieldView velocity,
                             const CartesianGeometryPlan& geometry,
                             const MeshPatch& patch) {
  for (std::int32_t z = -1; z <= patch.cells.z; ++z) {
    const double zc = extended_centre(geometry.z(), patch.begin.z + z);
    for (std::int32_t y = -1; y <= patch.cells.y; ++y) {
      const double yc = extended_centre(geometry.y(), patch.begin.y + y);
      for (std::int32_t x = -1; x <= patch.cells.x; ++x) {
        const double xc = extended_centre(geometry.x(), patch.begin.x + x);
        const Int3 index{x, y, z};
        velocity.unchecked(index, 0U) =
            1.0 + 2.0 * xc + 3.0 * xc * xc + 4.0 * yc +
            5.0 * yc * yc + 6.0 * zc + 7.0 * zc * zc;
        velocity.unchecked(index, 1U) =
            -2.0 + 7.0 * xc - 2.0 * xc * xc - 11.0 * yc +
            3.0 * yc * yc + 13.0 * zc - 4.0 * zc * zc;
        velocity.unchecked(index, 2U) =
            4.0 - 17.0 * xc + 5.0 * xc * xc + 19.0 * yc -
            6.0 * yc * yc - 23.0 * zc + 8.0 * zc * zc;
      }
    }
  }
}

bool revise_complete_attempt(AttemptTransaction& transaction,
                             FieldId velocity, FieldId pressure) {
  return transaction.revise_trial(velocity) &&
         transaction.revise_trial(pressure);
}

bool test_derived_gradient_lifecycle() {
  bool passed = true;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, uniform_mesh(), GeometryBudget{}, geometry,
                       patch)),
                   "gradient geometry compiles");
  GradientFixture fixture;
  passed &= expect(make_gradient_fixture(patch.cells, fixture),
                   "gradient field/state fixture compiles");

  const std::array declared{fixture.velocity, fixture.pressure};
  constexpr RevisionSlotId cache_slot = 0U;
  constexpr RevisionSourceId geometry_source = 70001U;
  constexpr RevisionSourceId boundary_source = 70002U;
  constexpr RevisionSourceId turbulence_source = 70003U;
  std::size_t cell_count = static_cast<std::size_t>(patch.cells.x) *
                           static_cast<std::size_t>(patch.cells.y) *
                           static_cast<std::size_t>(patch.cells.z);
  DerivedFieldPlan derived;
  passed &= expect(static_cast<bool>(DerivedFieldPlan::compile(
                       fixture.velocity, {declared.data(), declared.size()},
                       cache_slot, geometry_source, boundary_source,
                       turbulence_source, cell_count, derived)),
                   "derived plan reserves active and pending gradient buffers");

  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       fixture.schema.size(), 1U,
                       fixture.schema.size() + 3U, transaction)),
                   "gradient transaction reserves exact cache dependencies");
  passed &= expect(static_cast<bool>(transaction.begin(fixture.layers)) &&
                       revise_complete_attempt(transaction, fixture.velocity,
                                               fixture.pressure),
                   "gradient attempt begins with complete write set");
  FieldView trial_velocity;
  passed &= expect(static_cast<bool>(fixture.layers.view(
                       StateRole::trial, fixture.velocity, trial_velocity)),
                   "trial velocity view is available");
  fill_linear_velocity(trial_velocity, geometry, patch);
  constexpr PlanFingerprint boundary_plan_identity = 11001U;
  constexpr PlanFingerprint turbulence_plan_identity = 21001U;
  const DerivedRevisionTuple first = make_derived_revision_tuple(
      as_const(trial_velocity), geometry, patch, 11U,
      boundary_plan_identity, 21U, turbulence_plan_identity);
  Span<const VelocityGradient> pending;
  passed &= expect(static_cast<bool>(derived.prepare_velocity_gradient(
                       as_const(trial_velocity), geometry, patch, first,
                       transaction, pending)),
                   "first gradient is prepared as pending");
  passed &= expect(pending.size == cell_count &&
                       derived.gradient_compute_count() == 1U,
                   "one complete local gradient field is computed");
  passed &= expect(
      derived.finalize_velocity_gradient(transaction).code ==
          StatusCode::invalid_plan,
      "cache finalize rejects an active attempt before collective outcome");
  Span<const VelocityGradient> still_pending;
  passed &= expect(static_cast<bool>(derived.prepare_velocity_gradient(
                       as_const(trial_velocity), geometry, patch, first,
                       transaction, still_pending)) &&
                       still_pending.data == pending.data &&
                       derived.gradient_compute_count() == 1U,
                   "early finalize rejection preserves pending cache bytes");
  GradientFixture wrong_fixture;
  AttemptTransaction wrong_transaction;
  passed &= expect(make_gradient_fixture(patch.cells, wrong_fixture) &&
                       static_cast<bool>(AttemptTransaction::create(
                           wrong_fixture.schema.size(), 1U,
                           wrong_fixture.schema.size() + 3U,
                           wrong_transaction)) &&
                       static_cast<bool>(wrong_transaction.begin(
                           wrong_fixture.layers)) &&
                       revise_complete_attempt(wrong_transaction,
                                               wrong_fixture.velocity,
                                               wrong_fixture.pressure),
                   "an unrelated transaction begins independently");
  passed &= expect(
      derived.prepare_velocity_gradient(
          as_const(trial_velocity), geometry, patch, first,
          wrong_transaction, still_pending).code == StatusCode::invalid_plan,
      "pending cache reuse rejects a different active transaction identity");
  passed &= expect(static_cast<bool>(wrong_transaction.collective_finish(
                       MPI_COMM_SELF, Status{})),
                   "the unrelated transaction reaches a committed outcome");
  passed &= expect(
      derived.finalize_velocity_gradient(wrong_transaction).code ==
          StatusCode::invalid_plan,
      "cache finalize rejects a different committed transaction identity");
  passed &= expect(static_cast<bool>(derived.prepare_velocity_gradient(
                       as_const(trial_velocity), geometry, patch, first,
                       transaction, still_pending)) &&
                       still_pending.data == pending.data,
                   "wrong-transaction rejection preserves pending ownership");
  const std::array expected{2.0, 3.0, 5.0, 7.0, -11.0, 13.0,
                            -17.0, 19.0, -23.0};
  for (const double actual : pending.data[cell_count / 2U].value) {
    passed &= expect(std::isfinite(actual),
                     "analytic gradient entries are finite");
  }
  for (std::size_t component = 0U; component < expected.size(); ++component) {
    passed &= expect(close(pending.data[cell_count / 2U].value[component],
                           expected[component]),
                     "linear velocity has exact Cartesian gradient");
  }
  Span<const VelocityGradient> repeated;
  Status repeated_status;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    repeated_status = derived.prepare_velocity_gradient(
        as_const(trial_velocity), geometry, patch, first, transaction,
        repeated);
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(repeated_status) &&
                       repeated.data == pending.data &&
                       derived.gradient_compute_count() == 1U &&
                       hot_allocations == 0U,
                   "matching pending tuple reuses one buffer with zero allocation");
  Span<const VelocityGradient> invisible;
  passed &= expect(derived.velocity_gradient(first, invisible).code ==
                       StatusCode::invalid_plan,
                   "pending gradient is invisible to committed readers");
  passed &= expect(static_cast<bool>(transaction.collective_finish(
                       MPI_COMM_SELF, Status{})) &&
                       static_cast<bool>(
                           derived.finalize_velocity_gradient(transaction)),
                   "successful transaction publishes pending gradient bytes");
  Span<const VelocityGradient> active;
  passed &= expect(static_cast<bool>(derived.velocity_gradient(first, active)) &&
                       active.size == cell_count &&
                       derived.cache_revision() != 0U,
                   "exact committed tuple reads the active cache");

  const auto rejects_stale = [&](DerivedRevisionTuple stale,
                                 std::string_view description) {
    return expect(derived.velocity_gradient(stale, invisible).code ==
                      StatusCode::invalid_plan,
                  description);
  };
  DerivedRevisionTuple stale_velocity = first;
  ++stale_velocity.velocity;
  passed &= rejects_stale(stale_velocity,
                          "a velocity revision change rejects stale gradients");
  DerivedRevisionTuple stale_geometry = first;
  ++stale_geometry.geometry;
  passed &= rejects_stale(stale_geometry,
                          "a geometry revision change rejects stale gradients");
  DerivedRevisionTuple stale_boundary = first;
  ++stale_boundary.boundary;
  passed &= rejects_stale(stale_boundary,
                          "a boundary revision change rejects stale gradients");
  DerivedRevisionTuple stale_turbulence = first;
  ++stale_turbulence.turbulence;
  passed &= rejects_stale(
      stale_turbulence,
      "a turbulence revision change rejects stale gradients");
  DerivedRevisionTuple stale_storage = first;
  ++stale_storage.velocity_storage;
  passed &= rejects_stale(
      stale_storage, "a velocity storage change rejects stale gradients");
  DerivedRevisionTuple stale_domain = first;
  ++stale_domain.velocity_revision_domain;
  passed &= rejects_stale(
      stale_domain, "a revision-domain change rejects stale gradients");
  DerivedRevisionTuple stale_patch = first;
  ++stale_patch.patch_identity;
  passed &= rejects_stale(stale_patch,
                          "a patch identity change rejects stale gradients");
  DerivedRevisionTuple stale_geometry_plan = first;
  ++stale_geometry_plan.geometry_plan_identity;
  passed &= rejects_stale(
      stale_geometry_plan,
      "a geometry-plan change rejects an equal-token stale gradient");
  DerivedRevisionTuple stale_boundary_plan = first;
  ++stale_boundary_plan.boundary_plan_identity;
  passed &= rejects_stale(
      stale_boundary_plan,
      "a boundary-plan change rejects an equal-token stale gradient");
  DerivedRevisionTuple stale_turbulence_plan = first;
  ++stale_turbulence_plan.turbulence_plan_identity;
  passed &= rejects_stale(
      stale_turbulence_plan,
      "a turbulence-plan change rejects an equal-token stale gradient");

  passed &= expect(static_cast<bool>(transaction.begin(fixture.layers)) &&
                       revise_complete_attempt(transaction, fixture.velocity,
                                               fixture.pressure),
                   "a rejected follow-on attempt begins");
  passed &= expect(static_cast<bool>(fixture.layers.view(
                       StateRole::trial, fixture.velocity, trial_velocity)),
                   "rejected attempt gets a fresh velocity revision");
  fill_linear_velocity(trial_velocity, geometry, patch);
  const DerivedRevisionTuple rejected = make_derived_revision_tuple(
      as_const(trial_velocity), geometry, patch, 12U,
      boundary_plan_identity, 22U, turbulence_plan_identity);
  passed &= expect(static_cast<bool>(derived.prepare_velocity_gradient(
                       as_const(trial_velocity), geometry, patch, rejected,
                       transaction, pending)) &&
                       derived.gradient_compute_count() == 2U,
                   "changed exact tuple computes one pending replacement");
  const RevisionToken active_revision = derived.cache_revision();
  passed &= expect(transaction.collective_finish(
                       MPI_COMM_SELF,
                       Status{StatusCode::rejected_step, 9U}).code ==
                       StatusCode::rejected_step &&
                       static_cast<bool>(
                           derived.finalize_velocity_gradient(transaction)),
                   "rollback explicitly discards the pending buffer");
  passed &= expect(derived.cache_revision() == active_revision &&
                       static_cast<bool>(derived.velocity_gradient(first, active)) &&
                       derived.velocity_gradient(rejected, invisible).code ==
                           StatusCode::invalid_plan,
                   "rollback preserves prior active bytes and revision");
  return passed;
}

bool test_quadratic_gradient(const CartesianMeshSpec& mesh,
                             bool expect_stretched) {
  bool passed = true;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  passed &= expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                       MPI_COMM_SELF, mesh, GeometryBudget{}, geometry,
                       patch)),
                   "quadratic-gradient geometry compiles");
  if (!passed) {
    return false;
  }
  passed &= expect(geometry.x().uniform() != expect_stretched,
                   "quadratic-gradient test selects the requested metric path");

  GradientFixture fixture;
  passed &= expect(make_gradient_fixture(patch.cells, fixture),
                   "quadratic-gradient field fixture compiles");
  const std::array declared{fixture.velocity, fixture.pressure};
  const std::size_t cell_count =
      static_cast<std::size_t>(patch.cells.x) *
      static_cast<std::size_t>(patch.cells.y) *
      static_cast<std::size_t>(patch.cells.z);
  DerivedFieldPlan derived;
  passed &= expect(static_cast<bool>(DerivedFieldPlan::compile(
                       fixture.velocity, {declared.data(), declared.size()},
                       0U, 71001U, 71002U, 71003U, cell_count, derived)),
                   "quadratic-gradient cache storage compiles");
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       fixture.schema.size(), 1U,
                       fixture.schema.size() + 3U, transaction)) &&
                       static_cast<bool>(transaction.begin(fixture.layers)) &&
                       revise_complete_attempt(transaction, fixture.velocity,
                                               fixture.pressure),
                   "quadratic-gradient attempt begins");
  FieldView velocity;
  passed &= expect(static_cast<bool>(fixture.layers.view(
                       StateRole::trial, fixture.velocity, velocity)),
                   "quadratic-gradient trial velocity is available");
  if (!passed) {
    return false;
  }
  fill_quadratic_velocity(velocity, geometry, patch);
  const DerivedRevisionTuple revisions = make_derived_revision_tuple(
      as_const(velocity), geometry, patch, 31U, 31001U, 41U, 41001U);
  Span<const VelocityGradient> gradients;
  Status gradient_status;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    gradient_status = derived.prepare_velocity_gradient(
        as_const(velocity), geometry, patch, revisions, transaction,
        gradients);
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(gradient_status) &&
                       gradients.size == cell_count && hot_allocations == 0U,
                   "quadratic-gradient kernel is allocation-free");
  if (!gradient_status) {
    return false;
  }

  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    const double zc = extended_centre(geometry.z(), patch.begin.z + z);
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      const double yc = extended_centre(geometry.y(), patch.begin.y + y);
      for (std::int32_t x = 0; x < patch.cells.x; ++x, ++flat) {
        const double xc = extended_centre(geometry.x(), patch.begin.x + x);
        const std::array expected{
            2.0 + 6.0 * xc, 4.0 + 10.0 * yc, 6.0 + 14.0 * zc,
            7.0 - 4.0 * xc, -11.0 + 6.0 * yc, 13.0 - 8.0 * zc,
            -17.0 + 10.0 * xc, 19.0 - 12.0 * yc,
            -23.0 + 16.0 * zc};
        for (std::size_t entry = 0U; entry < expected.size(); ++entry) {
          passed &= expect(close(gradients.data[flat].value[entry],
                                 expected[entry], 2.0e-11, 2.0e-12),
                           "three-point gradient is exact for a quadratic field");
        }
      }
    }
  }
  passed &= expect(static_cast<bool>(transaction.collective_finish(
                       MPI_COMM_SELF, Status{})) &&
                       static_cast<bool>(
                           derived.finalize_velocity_gradient(transaction)),
                   "quadratic-gradient cache commits normally");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_constant_and_wilke_transport();
  passed &= test_sutherland_temperature_dependence();
  passed &= test_transport_rejects_mismatched_thermodynamics();
  passed &= test_derived_gradient_lifecycle();
  passed &= test_quadratic_gradient(uniform_mesh(), false);
  passed &= test_quadratic_gradient(stretched_mesh(), true);
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 transport and derived-cache tests passed\n";
  return 0;
}
