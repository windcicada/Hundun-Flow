// SPDX-License-Identifier: Apache-2.0

#ifndef HUNDUN_V04_ENABLE_TEST_ACCESS
#define HUNDUN_V04_ENABLE_TEST_ACCESS 1
#endif
#include "hundun/v04_flow.hpp"
#include "../support/candidate_boundary_fixture.hpp"
#include "../support/ibm_force_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <vector>

namespace candidate_boundary_allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed))
    count.fetch_add(1U, std::memory_order_relaxed);
}

void* allocate(std::size_t size) {
  observe();
  if (void* pointer = std::malloc(size == 0U ? 1U : size)) return pointer;
  throw std::bad_alloc{};
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* pointer = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&pointer, alignment, requested) == 0 &&
      pointer != nullptr)
    return pointer;
  throw std::bad_alloc{};
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

}  // namespace candidate_boundary_allocation_observer

void* operator new(std::size_t size) {
  return candidate_boundary_allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return candidate_boundary_allocation_observer::allocate(size);
}
void* operator new(std::size_t size, std::align_val_t alignment) {
  return candidate_boundary_allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t size, std::align_val_t alignment) {
  return candidate_boundary_allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return candidate_boundary_allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return candidate_boundary_allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t,
                     std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool same_double_bits(double left, double right) {
  std::uint64_t left_bits = 0U;
  std::uint64_t right_bits = 0U;
  std::memcpy(&left_bits, &left, sizeof(left));
  std::memcpy(&right_bits, &right, sizeof(right));
  return left_bits == right_bits;
}

struct OwnedField {
  std::vector<double> bytes;
  FieldView view{};
};

struct OwnedFaceField {
  std::vector<double> bytes;
  FaceFieldView view{};
};

class IdentityCertificateOperator final : public LinearOperator {
 public:
  explicit IdentityCertificateOperator(LinearOperatorCertificate certificate)
      : certificate_(certificate) {}

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_;
  }

  Status apply(FieldView input, FieldView output) const noexcept override {
    const Int3 cells = certificate_.local_shape;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          output.unchecked({x, y, z}, 0U) =
              input.unchecked({x, y, z}, 0U);
    return {};
  }

 private:
  LinearOperatorCertificate certificate_{};
};

OwnedField make_field(FieldId field, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision,
                      StorageIdentity storage) {
  OwnedField owned;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  owned.bytes.assign(nx * ny * nz * components, 0.0);
  owned.view.base = owned.bytes.data() + ghosts + ghosts * nx +
                    ghosts * nx * ny;
  owned.view.interior = cells;
  owned.view.ghosts = {ghosts, ghosts, ghosts};
  owned.view.components = components;
  owned.view.stride_y = nx;
  owned.view.stride_z = nx * ny;
  owned.view.component_stride = nx * ny * nz;
  owned.view.field = field;
  owned.view.revision = revision;
  owned.view.storage_identity = storage;
  owned.view.revision_domain = 15001U;
  return owned;
}

OwnedFaceField make_face(CartesianAxis axis, Int3 cells,
                         StorageIdentity storage) {
  OwnedFaceField owned;
  Int3 extents = cells;
  if (axis == CartesianAxis::x) {
    ++extents.x;
  } else if (axis == CartesianAxis::y) {
    ++extents.y;
  } else {
    ++extents.z;
  }
  owned.bytes.assign(static_cast<std::size_t>(extents.x) * extents.y *
                         extents.z,
                     0.0);
  owned.view = {owned.bytes.data(),
                extents,
                static_cast<std::size_t>(extents.x),
                static_cast<std::size_t>(extents.x) * extents.y,
                axis,
                storage,
                15002U};
  return owned;
}

void fill(OwnedField& field, double value) {
  std::fill(field.bytes.begin(), field.bytes.end(), value);
}

Status fill_exact_thermodynamic_candidate(
    const ThermodynamicsPlan& thermodynamics, double pressure_reference,
    ConstFieldView base_pressure, ConstFieldView pressure_correction,
    ConstFieldView candidate_enthalpy, ConstFieldView candidate_velocity,
    Span<const ConstFieldView> independent_species, FieldView density,
    FieldView temperature, FieldView pressure_compressibility) {
  const Int3 cells = density.interior;
  if (independent_species.size !=
          thermodynamics.independent_species_count() ||
      (independent_species.size != 0U &&
       independent_species.data == nullptr))
    return {StatusCode::invalid_plan, 1U};
  std::vector<double> composition(independent_species.size);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        for (std::size_t species = 0U; species < independent_species.size;
             ++species)
          composition[species] =
              independent_species.data[species].unchecked(cell, 0U);
        ThermoState replayed;
        const Status status =
            thermodynamics.evaluate_from_reference_pressure(
                pressure_reference,
                base_pressure.unchecked(cell, 0U) +
                    pressure_correction.unchecked(cell, 0U),
                candidate_enthalpy.unchecked(cell, 0U),
                {composition.data(), composition.size()},
                {candidate_velocity.unchecked(cell, 0U),
                 candidate_velocity.unchecked(cell, 1U),
                 candidate_velocity.unchecked(cell, 2U)},
                replayed);
        if (!status) return status;
        density.unchecked(cell, 0U) = replayed.rho;
        temperature.unchecked(cell, 0U) = replayed.temperature;
        pressure_compressibility.unchecked(cell, 0U) =
            replayed.drho_dp_hY;
      }
  return {};
}

void fill_face_flux(FaceFluxView flux, double value) {
  const std::array<FaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (FaceFieldView face : faces) {
    for (std::int32_t z = 0; z < face.extents.z; ++z) {
      for (std::int32_t y = 0; y < face.extents.y; ++y) {
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          face.unchecked({x, y, z}) = value;
        }
      }
    }
  }
}

bool face_flux_all_equal(ConstFaceFluxView flux, double value) {
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (ConstFaceFieldView face : faces)
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          if (!same_double_bits(face.unchecked({x, y, z}), value))
            return false;
  return true;
}

std::vector<double> face_flux_values(ConstFaceFluxView flux) {
  std::vector<double> values;
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (ConstFaceFieldView face : faces) {
    for (std::int32_t z = 0; z < face.extents.z; ++z) {
      for (std::int32_t y = 0; y < face.extents.y; ++y) {
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          values.push_back(face.unchecked({x, y, z}));
        }
      }
    }
  }
  return values;
}

struct CommittedFluxHistory {
  FieldId dependency{};
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage storage;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
};

bool initialize_flux_history(Int3 cells, CommittedFluxHistory& history) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("piso.authority_flux_history", 1U, 0U,
                              history.dependency) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{ArenaFieldRequest{
      history.dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  return static_cast<bool>(ArenaLayout::compile(
             schema, {requests.data(), requests.size()}, layout)) &&
         static_cast<bool>(StateLayers::allocate(layout, history.layers)) &&
         static_cast<bool>(AttemptTransaction::create(
             history.layers.field_count(), 1U, history.layers.field_count(),
             history.transaction)) &&
         static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                           history.storage)) &&
         static_cast<bool>(history.authority.claim(
             71U, 0U, history.transaction, history.writer));
}

bool commit_zero_flux(const CartesianKernelPlan& kernels, Int3 cells,
                      RevisionToken seed, CommittedFluxHistory& history,
                      ConstFaceFluxView& committed) {
  if (!history.transaction.begin(history.layers) ||
      !history.transaction.revise_trial(history.dependency)) {
    return false;
  }
  const RevisionDependency dependency{
      AttemptTransaction::field_revision_source(history.dependency),
      history.transaction.trial_revision(history.dependency)};
  PendingFaceFluxView pending;
  if (!history.writer.begin_pending(history.transaction, history.storage,
                                    pending)) {
    return false;
  }
  OwnedField density = make_field(120U, cells, 1U, 2U, seed, seed + 1000U);
  OwnedField velocity =
      make_field(121U, cells, 3U, 2U, seed + 1U, seed + 1001U);
  fill(density, 1.0);
  fill(velocity, 0.0);
  const std::array<ConstFieldView, 2U> reads{
      as_const(density.view), as_const(velocity.view)};
  const KernelInvocation invocation{
      {reads.data(), reads.size()}, {}, {{0, 0, 0}, cells},
      0U, 0U, 1U, 0U, nullptr};
  const std::array dependencies{dependency};
  return static_cast<bool>(reconstruct_mass_flux(kernels, invocation,
                                                 pending)) &&
         static_cast<bool>(history.writer.publish_pending(
             {dependencies.data(), dependencies.size()}, pending)) &&
         static_cast<bool>(history.transaction.collective_finish(
             MPI_COMM_SELF, Status{})) &&
         static_cast<bool>(history.writer.committed(history.storage,
                                                     committed));
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  // Keep this fixture deliberately non-unit and anisotropic.  The PISO
  // workspace stores V/A_U, so a unit-cube oracle would not distinguish the
  // integral diagonal from the legacy unit-volume inverse.
  mesh.upper = {2.0, 3.0, 4.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {4, 4, 4};
  mesh.minimum_spacing = {0.5, 0.75, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 64U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 24U;
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

ThermophysicalSpec thermophysical_spec() {
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
  CartesianKernelPlan kernels;
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  EquationPlanSet equations;
};

double cell_width(const Fixture& fixture, CartesianAxis axis,
                  std::int32_t index) {
  const Span<const double> widths = fixture.geometry.axis(axis).widths();
  return widths.data[static_cast<std::size_t>(index)];
}

double cell_volume(const Fixture& fixture, Int3 cell) {
  return cell_width(fixture, CartesianAxis::x, cell.x) *
         cell_width(fixture, CartesianAxis::y, cell.y) *
         cell_width(fixture, CartesianAxis::z, cell.z);
}

double time_step_for_bdf(BdfCoefficients bdf) {
  if (bdf.order == 1U) return 1.0 / bdf.a0;
  const double ratio = 1.0 / (std::sqrt(-bdf.a1 / bdf.a2) - 1.0);
  return (1.0 + ratio) / -bdf.a1;
}

double face_area(const Fixture& fixture, CartesianAxis axis, Int3 cell) {
  if (axis == CartesianAxis::x) {
    return cell_width(fixture, CartesianAxis::y, cell.y) *
           cell_width(fixture, CartesianAxis::z, cell.z);
  }
  if (axis == CartesianAxis::y) {
    return cell_width(fixture, CartesianAxis::x, cell.x) *
           cell_width(fixture, CartesianAxis::z, cell.z);
  }
  return cell_width(fixture, CartesianAxis::x, cell.x) *
         cell_width(fixture, CartesianAxis::y, cell.y);
}

bool make_fixture(Fixture& out, bool open_boundary = false,
                  bool require_two_layer_boundary = false) {
  const CartesianMeshSpec mesh = mesh_spec();
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = 0x15010001U;
  model.pressure_reference =
      open_boundary ? PressureReferenceKind::boundary_absolute
                    : PressureReferenceKind::closed_mass;
  if (open_boundary) {
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::slip;
      face.thermal_kind = BoundaryKind::none;
    }
    model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
    model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
    model.boundaries[0U].temperature = 300.0;
    model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
    model.boundaries[1U].pressure = 101325.0;
  } else {
    for (BoundaryFaceSpec& face : model.boundaries) {
      face.flow_kind = BoundaryKind::periodic;
      face.thermal_kind = BoundaryKind::none;
    }
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  // E_h currently requests a two-layer homogeneous variation boundary even
  // when its own reconstruction is central2.  A test that executes E_h can
  // opt into that plan contract through an otherwise-unused scalar scheme.
  model.schemes.passive_scalar = require_two_layer_boundary
                                     ? ConvectionScheme::tvd2
                                     : ConvectionScheme::central2;

  FieldRegistry registry;
  FieldId id = 0U;
  if (!registry.require_field("rho", 1U, 2U, id) || id != 0U ||
      !registry.require_field("U", 3U, 2U, id) || id != 1U ||
      !registry.require_field("pi", 1U, 2U, id) || id != 2U ||
      !registry.require_field("h", 1U, 2U, id) || id != 3U ||
      !registry.require_field("T", 1U, 2U, id) || id != 4U ||
      !CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {},
                                          out.geometry, out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry,
                                 out.patch, registry, out.boundary,
                                 out.schemes, out.time) ||
      !CartesianKernelPlan::compile(out.schemes, out.geometry, out.patch,
                                    out.boundary, out.kernels)) {
    return false;
  }
  const ThermophysicalSpec thermophysics = thermophysical_spec();
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
  EquationPlanSpec equations;
  equations.density = 0U;
  equations.velocity = 1U;
  equations.pressure_perturbation = 2U;
  equations.enthalpy = 3U;
  equations.temperature = 4U;
  equations.effective_viscosity = 5U;
  equations.pressure_compressibility = 6U;
  equations.velocity_gradient = 7U;
  equations.pressure_reference = model.pressure_reference;
  equations.closed_mass_service_stage = open_boundary ? 0U : 1U;
  equations.maximum_cells_per_rank = 64U;
  return static_cast<bool>(EquationPlanSet::compile(
      MPI_COMM_SELF, out.schemes, out.geometry, out.patch, out.boundary,
      out.contributions, out.thermodynamics, out.transport, equations,
      out.equations));
}

PisoPlanSpec valid_spec() {
  PisoPlanSpec spec;
  spec.pressure_correctors = 2U;
  spec.pressure_stage = 21U;
  spec.final_flux_slot = 0U;
  spec.pressure_solve = {1.0e-15, 1.0e-13, 400U, 4U, 12U};
  spec.eos_tolerance = 1.0e-10;
  spec.continuity_tolerance = 1.0e-10;
  spec.energy_tolerance = 1.0e-10;
  spec.closed_mass_tolerance = 1.0e-10;
  spec.gauge_tolerance = 1.0e-12;
  return spec;
}

bool test_exactly_two_compile_authority() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture), "PISO dependency fixture compiles");
  if (!passed) {
    return false;
  }
  PisoPlan plan;
  const PisoPlanSpec accepted = valid_spec();
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, accepted, plan)) &&
                       plan.fingerprint() != 0U &&
                       plan.pressure_correctors() == 2U &&
                       plan.pressure_algorithm() == LinearAlgorithm::fgmres &&
                       plan.pressure_stage() == accepted.pressure_stage &&
                       plan.final_flux_slot() == accepted.final_flux_slot,
                   "PisoPlan freezes exactly two correctors and one writer slot");

  for (const std::uint8_t count : {0U, 1U, 3U, 255U}) {
    PisoPlanSpec mutation = accepted;
    mutation.pressure_correctors = count;
    PisoPlan rejected;
    passed &= expect(PisoPlan::compile(MPI_COMM_SELF, fixture.equations,
                                      mutation, rejected)
                             .code == StatusCode::invalid_plan &&
                         rejected.fingerprint() == 0U,
                     "every non-two corrector count rejects cold");
  }

  PisoPlanSpec mutation = accepted;
  mutation.pressure_stage = 0U;
  PisoPlan rejected;
  passed &= expect(PisoPlan::compile(MPI_COMM_SELF, fixture.equations,
                                    mutation, rejected)
                           .code == StatusCode::invalid_plan,
                   "stage zero cannot own pressure correction");
  mutation = accepted;
  mutation.pressure_solve.restart = 0U;
  passed &= expect(PisoPlan::compile(MPI_COMM_SELF, fixture.equations,
                                    mutation, rejected)
                           .code == StatusCode::invalid_plan,
                   "fixed FGMRES pressure solve rejects a zero restart");
  mutation = accepted;
  mutation.energy_tolerance =
      std::numeric_limits<double>::quiet_NaN();
  passed &= expect(PisoPlan::compile(MPI_COMM_SELF, fixture.equations,
                                    mutation, rejected)
                           .code == StatusCode::invalid_plan,
                   "non-finite terminal energy gate rejects cold");
  mutation = accepted;
  mutation.continuity_tolerance =
      std::numeric_limits<double>::quiet_NaN();
  passed &= expect(PisoPlan::compile(MPI_COMM_SELF, fixture.equations,
                                    mutation, rejected)
                           .code == StatusCode::invalid_plan,
                   "non-finite terminal gate rejects cold");

  PisoAttemptReport report;
  static_assert(report.pressure.size() == 2U,
                "PISO report has exactly two pressure solve records");
  passed &= expect(report.pressure_solve_calls == 0U &&
                       report.energy_residual == 0.0 &&
                       report.final_flux_revision == 0U,
                   "attempt report starts unpublished");
  PisoPlanSpec legacy = accepted;
  legacy.energy_tolerance = 0.0;
  PisoPlan legacy_plan;
  passed &= expect(
      static_cast<bool>(PisoPlan::compile(MPI_COMM_SELF, fixture.equations,
                                          legacy, legacy_plan)),
      "zero energy tolerance preserves the legacy terminal contract");
  return passed;
}

bool same_intermediate(const PisoIntermediateCertificate& left,
                       const PisoIntermediateCertificate& right) {
  return left.plan == right.plan && left.r_au == right.r_au &&
         left.h_by_a == right.h_by_a &&
         left.pressure_face_coefficient ==
             right.pressure_face_coefficient &&
         left.phi_h_by_a == right.phi_h_by_a &&
         left.trial_face_flux == right.trial_face_flux &&
         left.temporal_face_flux == right.temporal_face_flux &&
         left.committed_face_history == right.committed_face_history &&
         left.dependency == right.dependency &&
         left.corrector == right.corrector &&
         left.thermophysical_boundary_semantics ==
             right.thermophysical_boundary_semantics &&
         left.thermophysical_boundary_target ==
             right.thermophysical_boundary_target &&
         left.thermophysical_boundary_rank_local_binding ==
             right.thermophysical_boundary_rank_local_binding &&
         left.thermophysical_boundary_collective_lineage ==
             right.thermophysical_boundary_collective_lineage &&
         left.thermophysical_boundary_rank_local_lineage ==
             right.thermophysical_boundary_rank_local_lineage &&
         left.pressure_energy_refinement ==
             right.pressure_energy_refinement &&
         left.pressure_energy_refinement_collective_lineage ==
             right.pressure_energy_refinement_collective_lineage &&
         left.pressure_energy_refinement_lineage ==
             right.pressure_energy_refinement_lineage;
}

bool same_state_correction(const PisoStateCorrectionCertificate& left,
                           const PisoStateCorrectionCertificate& right) {
  const auto same_reference = [](const PressureReferenceCertificate& lhs,
                                 const PressureReferenceCertificate& rhs) {
    return lhs.plan == rhs.plan && lhs.predictor == rhs.predictor &&
           lhs.thermodynamics == rhs.thermodynamics &&
           lhs.closure == rhs.closure && lhs.time == rhs.time &&
           lhs.pressure_reference == rhs.pressure_reference &&
           lhs.kind == rhs.kind;
  };
  return left.plan == right.plan &&
         left.pressure_system == right.pressure_system &&
         left.correction == right.correction &&
         left.enthalpy_correction == right.enthalpy_correction &&
         left.velocity == right.velocity && left.pressure == right.pressure &&
         left.enthalpy == right.enthalpy && left.density == right.density &&
         left.temperature == right.temperature &&
         left.face_flux == right.face_flux &&
         left.exact_eos_closure == right.exact_eos_closure &&
         left.state == right.state && left.corrector == right.corrector &&
         left.closure == right.closure &&
         left.thermophysical_boundary_semantics ==
             right.thermophysical_boundary_semantics &&
         left.thermophysical_boundary_target ==
             right.thermophysical_boundary_target &&
         left.thermophysical_boundary_rank_local_binding ==
             right.thermophysical_boundary_rank_local_binding &&
         left.thermophysical_boundary_collective_lineage ==
             right.thermophysical_boundary_collective_lineage &&
         left.thermophysical_boundary_rank_local_lineage ==
             right.thermophysical_boundary_rank_local_lineage &&
         same_reference(left.input_pressure_reference,
                        right.input_pressure_reference) &&
         same_reference(left.output_pressure_reference,
                        right.output_pressure_reference) &&
         left.closed_gauge_collective_transaction ==
             right.closed_gauge_collective_transaction &&
         left.closed_gauge_rank_local_transaction ==
             right.closed_gauge_rank_local_transaction;
}

Status prepare_closed_gauge(
    const PressureReferencePlan& reference_plan,
    const PressureCorrectionCertificate& pressure,
    const PressureReferenceCertificate& predecessor,
    double pressure_reference, ConstFieldView pressure_perturbation,
    ConstFieldView raw_pressure_correction,
    ConstFieldView candidate_pressure_compressibility,
    PlanFingerprint exact_eos_closure, ReductionEngine& reductions,
    ClosedGaugeCorrectionCertificate& certificate) {
  ClosedGaugeCorrectionPrepareInput input;
  input.predecessor = predecessor;
  input.pressure_reference = pressure_reference;
  input.corrector = pressure.corrector;
  input.time = pressure.time;
  input.geometry = pressure.geometry;
  input.pressure_correction_authority = pressure.state;
  input.target_thermodynamic_closure = exact_eos_closure;
  input.pressure_perturbation = pressure_perturbation;
  input.raw_pressure_correction = raw_pressure_correction;
  input.candidate_pressure_compressibility =
      candidate_pressure_compressibility;
  return reference_plan.prepare_closed_gauge_correction(
      input, reductions, certificate);
}

bool test_coupler_lifecycle_and_mutations() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture), "coupler fixture compiles");
  if (!passed) {
    return false;
  }
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "coupler PisoPlan compiles");
  const Int3 cells = fixture.patch.cells;
  const std::uint8_t reach = fixture.equations.kernels().reach();
  OwnedField density = make_field(0U, cells, 1U, reach, 1601U, 2601U);
  OwnedField velocity = make_field(1U, cells, 3U, 0U, 1602U, 2602U);
  OwnedField target_pressure =
      make_field(2U, cells, 1U, reach, 1608U, 2611U);
  OwnedField diagonal = make_field(30U, cells, 3U, 0U, 1603U, 2603U);
  OwnedField rhs = make_field(31U, cells, 3U, 0U, 1604U, 2604U);
  OwnedField r_au = make_field(40U, cells, 3U, 1U, 1605U, 2605U);
  OwnedField h_by_a = make_field(41U, cells, 3U, reach, 1606U, 2606U);
  OwnedField pressure_gradient =
      make_field(42U, cells, 3U, 0U, 1607U, 2610U);
  OwnedFaceField ax = make_face(CartesianAxis::x, cells, 2607U);
  OwnedFaceField ay = make_face(CartesianAxis::y, cells, 2608U);
  OwnedFaceField az = make_face(CartesianAxis::z, cells, 2609U);
  FaceFluxStorage flux_storage;
  FaceFluxView phi_h_by_a;
  FaceFluxView trial_flux_authority;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_workspace(
                       cells, 2U, flux_storage)) &&
                       static_cast<bool>(flux_storage.workspace_view(
                           0U, 1610U, phi_h_by_a)) &&
                       static_cast<bool>(flux_storage.workspace_view(
                           1U, 1702U, trial_flux_authority)),
                   "phiHbyA workspace allocates");
  fill(density, 1.0);
  fill(velocity, 0.0);
  fill(target_pressure, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        target_pressure.view.unchecked({x, y, z}, 0U) =
            -26.0 * static_cast<double>(x) / 9.0;
      }
    }
  }
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        diagonal.view.unchecked(cell, 0U) = 2.0;
        diagonal.view.unchecked(cell, 1U) = 4.0;
        diagonal.view.unchecked(cell, 2U) = 8.0;
        rhs.view.unchecked(cell, 0U) = 2.0;
        rhs.view.unchecked(cell, 1U) = 8.0;
        rhs.view.unchecked(cell, 2U) = 24.0;
      }
    }
  }
  const PisoCouplerWorkspace workspace{
      r_au.view, h_by_a.view, pressure_gradient.view,
      ax.view, ay.view, az.view, phi_h_by_a};
  const std::array<HaloFieldSpec, 3U> halo_fields{{
      {density.view.field, 1U, 1U},
      {r_au.view.field, 1U, 3U},
      {h_by_a.view.field, reach, 3U}}};
  HaloEngine halo;
  passed &= expect(static_cast<bool>(halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {halo_fields.data(), halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "PISO intermediate halo reserves cold");
  constexpr FieldId correction_field = 90U;
  const std::array<HaloFieldSpec, 1U> correction_halo_fields{{
      {correction_field, 1U, 1U}}};
  HaloEngine correction_halo;
  passed &= expect(static_cast<bool>(correction_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {correction_halo_fields.data(),
                        correction_halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "PISO correction halo reserves cold");
  const PisoCouplerServices services{MPI_COMM_SELF,
                                     &fixture.geometry,
                                     fixture.patch,
                                     &fixture.boundary,
                                     &fixture.thermodynamics,
                                     &halo,
                                     160U,
                                     density.view.field,
                                     &correction_halo,
                                     162U,
                                     correction_field};
  PressureVelocityCoupler coupler;
  passed &= expect(static_cast<bool>(PressureVelocityCoupler::bind(
                       piso, fixture.equations, services, workspace,
                       coupler)) &&
                       coupler.fingerprint() != 0U &&
                       coupler.workspace_storage_address() ==
                           reinterpret_cast<std::uintptr_t>(r_au.view.base),
                   "coupler binds disjoint preallocated workspace");

  CommittedFluxHistory flux_history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(
      initialize_flux_history(cells, flux_history) &&
          commit_zero_flux(fixture.equations.kernels(), cells, 1650U,
                           flux_history, previous_flux) &&
          commit_zero_flux(fixture.equations.kernels(), cells, 1660U,
                           flux_history, accepted_flux),
      "two certified committed face-flux history layers publish");

  PisoIntermediateInput input;
  input.momentum = {fixture.equations.momentum().fingerprint(),
                    EquationAssemblyScope::momentum_predictor,
                    1701U,
                    fixture.geometry.topology_revision(),
                    1702U,
                    1703U,
                    time_step_for_bdf({10.0, -15.0, 5.0, 2U})};
  input.predictor.plan =
      fixture.equations.thermophysical_predictor().fingerprint();
  input.predictor.time = 1701U;
  input.predictor.geometry = fixture.geometry.topology_revision();
  input.predictor.accepted_face_flux = accepted_flux.revision;
  input.predictor.previous_face_flux = previous_flux.revision;
  input.predictor.committed_face_flux_authority =
      accepted_flux.certificate.authority();
  input.predictor.committed_face_flux_storage =
      accepted_flux.certificate.storage();
  input.predictor.committed_face_flux_revision_domain =
      accepted_flux.certificate.revision_domain();
  input.predictor.predicted_density = density.view.revision;
  input.predictor.predicted_density_storage = density.view.storage_identity;
  input.predictor.predicted_density_revision_domain =
      density.view.revision_domain;
  input.predictor.paired_face_flux = trial_flux_authority.revision;
  input.predictor.paired_face_flux_storage =
      trial_flux_authority.x.storage_identity;
  input.predictor.paired_face_flux_revision_domain =
      trial_flux_authority.x.revision_domain;
  input.predictor.state = 1705U;
  input.predictor.order = 2U;
  input.pressure_reference = {
      fixture.equations.pressure_reference().fingerprint(),
      fixture.equations.thermophysical_predictor().fingerprint(),
      fixture.thermodynamics.fingerprint(),
      1706U,
      1701U,
      1707U,
      PressureReferenceKind::closed_mass};
  input.density = density.view;
  input.trial_velocity = as_const(velocity.view);
  input.trial_flux = as_const(trial_flux_authority);
  input.temporal_reference = as_const(phi_h_by_a);
  input.committed_face_history = {accepted_flux, previous_flux};
  input.momentum_system.diagonal = diagonal.view;
  input.momentum_system.rhs = rhs.view;
  input.bdf = {10.0, -15.0, 5.0, 2U};
  input.numeric_boundary = fixture.boundary.revision();
  input.thermophysical_boundary.binding.pressure_reference = 101325.0;
  input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(target_pressure.view);

  FaceFluxView cross_axis_output = phi_h_by_a;
  cross_axis_output.x.base =
      const_cast<double*>(accepted_flux.y.base);
  const PisoCouplerWorkspace cross_axis_workspace{
      r_au.view, h_by_a.view, pressure_gradient.view,
      ax.view, ay.view, az.view, cross_axis_output};
  PressureVelocityCoupler cross_axis_coupler;
  passed &= expect(
      static_cast<bool>(PressureVelocityCoupler::bind(
          piso, fixture.equations, services, cross_axis_workspace,
          cross_axis_coupler)),
      "cross-axis history-alias fixture binds before attempt inputs exist");
  PisoIntermediateInput cross_axis_input = input;
  cross_axis_input.corrector = 1U;
  cross_axis_input.temporal_reference = as_const(cross_axis_output);
  PisoIntermediateCertificate cross_axis_certificate;
  passed &= expect(
      cross_axis_coupler.refresh(cross_axis_input,
                                 cross_axis_certificate).code ==
              StatusCode::invalid_plan &&
          !cross_axis_certificate.valid(),
      "committed history rejects cross-axis output alias atomically");

  FieldView history_aliased_r_au = r_au.view;
  const double* const history_storage_begin =
      std::min(accepted_flux.x.base, previous_flux.x.base);
  const std::size_t r_au_prefix =
      static_cast<std::size_t>(history_aliased_r_au.ghosts.x) +
      static_cast<std::size_t>(history_aliased_r_au.ghosts.y) *
          history_aliased_r_au.stride_y +
      static_cast<std::size_t>(history_aliased_r_au.ghosts.z) *
          history_aliased_r_au.stride_z;
  history_aliased_r_au.base =
      const_cast<double*>(history_storage_begin) + r_au_prefix;
  const PisoCouplerWorkspace history_cell_alias_workspace{
      history_aliased_r_au, h_by_a.view, pressure_gradient.view,
      ax.view, ay.view, az.view, phi_h_by_a};
  PressureVelocityCoupler history_cell_alias_coupler;
  passed &= expect(
      static_cast<bool>(PressureVelocityCoupler::bind(
          piso, fixture.equations, services, history_cell_alias_workspace,
          history_cell_alias_coupler)),
      "cell-workspace history-alias fixture binds before attempt inputs exist");
  const std::vector<double> accepted_before =
      face_flux_values(accepted_flux);
  const std::vector<double> previous_before =
      face_flux_values(previous_flux);
  PisoIntermediateInput history_cell_alias_input = input;
  history_cell_alias_input.corrector = 1U;
  const PisoIntermediateCertificate preserved_certificate{
      11U, 12U, 13U, 14U, 15U, 16U, 17U, 18U, 19U, 1U};
  PisoIntermediateCertificate history_cell_alias_certificate =
      preserved_certificate;
  passed &= expect(
      history_cell_alias_coupler
                  .refresh(history_cell_alias_input,
                           history_cell_alias_certificate)
                  .code == StatusCode::invalid_plan &&
          same_intermediate(history_cell_alias_certificate,
                            preserved_certificate) &&
          face_flux_values(accepted_flux) == accepted_before &&
          face_flux_values(previous_flux) == previous_before,
      "committed history rejects writable cell-workspace alias before mutation");

  input.corrector = 2U;
  PisoIntermediateCertificate certificate;
  passed &= expect(coupler.refresh(input, certificate).code ==
                       StatusCode::invalid_plan &&
                       !certificate.valid(),
                   "corrector two cannot bypass corrector one");

  input.corrector = 1U;
  passed &= expect(static_cast<bool>(coupler.refresh(input, certificate)) &&
                       certificate.valid() && certificate.corrector == 1U,
                   "corrector one certifies current intermediates");
  ConstFaceFluxView inspected_intermediate_flux;
  passed &= expect(
      static_cast<bool>(coupler.inspect_intermediate_flux(
          certificate, inspected_intermediate_flux)) &&
          inspected_intermediate_flux.x.base == phi_h_by_a.x.base &&
          inspected_intermediate_flux.y.base == phi_h_by_a.y.base &&
          inspected_intermediate_flux.z.base == phi_h_by_a.z.base &&
          inspected_intermediate_flux.revision != 0U &&
          !inspected_intermediate_flux.certificate.valid(),
      "current intermediate certificate exposes attempt-local phiHbyA read-only");
  PisoIntermediateCertificate stale_intermediate_flux = certificate;
  ++stale_intermediate_flux.dependency;
  const ConstFaceFluxView preserved_intermediate_flux =
      inspected_intermediate_flux;
  passed &= expect(
      coupler.inspect_intermediate_flux(stale_intermediate_flux,
                                        inspected_intermediate_flux)
              .code == StatusCode::invalid_plan &&
          inspected_intermediate_flux.x.base ==
              preserved_intermediate_flux.x.base &&
          inspected_intermediate_flux.y.base ==
              preserved_intermediate_flux.y.base &&
          inspected_intermediate_flux.z.base ==
              preserved_intermediate_flux.z.base &&
          inspected_intermediate_flux.revision ==
              preserved_intermediate_flux.revision,
      "stale intermediate certificate cannot replace inspected phiHbyA view");
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double volume = cell_volume(fixture, cell);
        passed &= expect(
            std::abs(r_au.view.unchecked(cell, 0U) - volume / 2.0) <
                    1.0e-14 &&
                std::abs(r_au.view.unchecked(cell, 1U) - volume / 4.0) <
                    1.0e-14 &&
                std::abs(r_au.view.unchecked(cell, 2U) - volume / 8.0) <
                    1.0e-14 &&
                std::abs(h_by_a.view.unchecked(cell, 0U)) < 1.0e-14 &&
                std::abs(h_by_a.view.unchecked(cell, 1U)) < 1.0e-14 &&
                std::abs(h_by_a.view.unchecked(cell, 2U)) < 1.0e-14,
            "rAU uses the diagonal and HbyA uses the converged predictor");
      }
    }
  }
  const PisoIntermediateCertificate first = certificate;

  auto refresh_first = [&](PisoIntermediateInput selected,
                           PisoIntermediateCertificate& out) {
    selected.corrector = 1U;
    selected.prior_corrector = 0U;
    return coupler.refresh(selected, out);
  };
  PisoIntermediateInput mutation = input;
  ++mutation.momentum_system.rhs.revision;
  PisoIntermediateCertificate rhs_mutated;
  passed &= expect(static_cast<bool>(refresh_first(mutation, rhs_mutated)) &&
                       rhs_mutated.r_au == first.r_au &&
                       rhs_mutated.h_by_a != first.h_by_a &&
                       rhs_mutated.pressure_face_coefficient ==
                           first.pressure_face_coefficient &&
                       rhs_mutated.phi_h_by_a != first.phi_h_by_a,
                   "momentum RHS invalidates HbyA/phiHbyA but reuses rAU");

  mutation = input;
  ++mutation.trial_velocity.revision;
  PisoIntermediateCertificate velocity_mutated;
  passed &= expect(
      static_cast<bool>(refresh_first(mutation, velocity_mutated)) &&
          velocity_mutated.r_au == first.r_au &&
          velocity_mutated.h_by_a != first.h_by_a &&
          velocity_mutated.pressure_face_coefficient ==
              first.pressure_face_coefficient &&
          velocity_mutated.phi_h_by_a != first.phi_h_by_a,
      "trial velocity invalidates trial-dependent intermediates only");

  mutation = input;
  ++mutation.numeric_boundary;
  PisoIntermediateCertificate boundary_mutated;
  passed &= expect(
      refresh_first(mutation, boundary_mutated).code ==
              StatusCode::invalid_plan &&
          !boundary_mutated.valid(),
      "foreign numeric boundary fails closed before intermediates");

  mutation = input;
  ++mutation.pressure_reference.pressure_reference;
  PisoIntermediateCertificate pressure_reference_mutated;
  passed &= expect(
      static_cast<bool>(refresh_first(mutation, pressure_reference_mutated)) &&
          pressure_reference_mutated.r_au == first.r_au &&
          pressure_reference_mutated.h_by_a == first.h_by_a &&
          pressure_reference_mutated.pressure_face_coefficient !=
              first.pressure_face_coefficient &&
          pressure_reference_mutated.phi_h_by_a != first.phi_h_by_a,
      "p_ref authority invalidates pressure-face/phiHbyA only");

  mutation = input;
  ++mutation.density.revision;
  mutation.predictor.predicted_density = mutation.density.revision;
  PisoIntermediateCertificate density_mutated;
  passed &= expect(static_cast<bool>(refresh_first(mutation, density_mutated)) &&
                       density_mutated.r_au != first.r_au &&
                       density_mutated.h_by_a != first.h_by_a &&
                       density_mutated.pressure_face_coefficient !=
                           first.pressure_face_coefficient &&
                       density_mutated.phi_h_by_a != first.phi_h_by_a,
                   "density authority invalidates the complete dependent chain");

  mutation = input;
  mutation.bdf.a1 = -14.0;
  mutation.bdf.a2 = 4.0;
  PisoIntermediateCertificate history_bdf_mutated;
  passed &= expect(
      static_cast<bool>(refresh_first(mutation, history_bdf_mutated)) &&
          history_bdf_mutated.r_au == first.r_au &&
          history_bdf_mutated.h_by_a == first.h_by_a &&
          history_bdf_mutated.pressure_face_coefficient ==
              first.pressure_face_coefficient &&
          history_bdf_mutated.phi_h_by_a != first.phi_h_by_a,
      "history-only BDF mutation invalidates phiHbyA without rebuilding rAU");

  mutation = input;
  ++mutation.committed_face_history.accepted.revision;
  PisoIntermediateCertificate stale_committed_history = first;
  passed &= expect(
      refresh_first(mutation, stale_committed_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(stale_committed_history, first),
      "stale committed face-history revision rejects atomically");

  mutation = input;
  ++mutation.predictor.committed_face_flux_authority;
  PisoIntermediateCertificate authority_mismatched_history = first;
  passed &= expect(
      refresh_first(mutation, authority_mismatched_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(authority_mismatched_history, first),
      "committed face-history writer authority rejects independently");

  mutation = input;
  ++mutation.predictor.committed_face_flux_storage;
  PisoIntermediateCertificate storage_mismatched_history = first;
  passed &= expect(
      refresh_first(mutation, storage_mismatched_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(storage_mismatched_history, first),
      "committed face-history storage identity rejects independently");

  mutation = input;
  ++mutation.predictor.committed_face_flux_revision_domain;
  PisoIntermediateCertificate domain_mismatched_history = first;
  passed &= expect(
      refresh_first(mutation, domain_mismatched_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(domain_mismatched_history, first),
      "committed face-history revision domain rejects independently");

  mutation = input;
  mutation.committed_face_history.previous = {};
  PisoIntermediateCertificate missing_committed_history = first;
  passed &= expect(
      refresh_first(mutation, missing_committed_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(missing_committed_history, first),
      "BDF2 requires the certified previous face-history layer");

  mutation = input;
  mutation.committed_face_history.previous = accepted_flux;
  mutation.predictor.previous_face_flux = accepted_flux.revision;
  PisoIntermediateCertificate duplicate_committed_history = first;
  passed &= expect(
      refresh_first(mutation, duplicate_committed_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(duplicate_committed_history, first),
      "BDF2 rejects one committed layer presented as both history levels");

  CommittedFluxHistory foreign_flux_history;
  ConstFaceFluxView foreign_previous_flux;
  ConstFaceFluxView foreign_accepted_flux;
  passed &= expect(
      initialize_flux_history(cells, foreign_flux_history) &&
          commit_zero_flux(fixture.equations.kernels(), cells, 1670U,
                           foreign_flux_history, foreign_previous_flux) &&
          commit_zero_flux(fixture.equations.kernels(), cells, 1680U,
                           foreign_flux_history, foreign_accepted_flux),
      "foreign committed face-history authority publishes matching revisions");
  mutation = input;
  mutation.committed_face_history = {foreign_accepted_flux,
                                     foreign_previous_flux};
  PisoIntermediateCertificate foreign_committed_history = first;
  passed &= expect(
      refresh_first(mutation, foreign_committed_history).code ==
              StatusCode::invalid_plan &&
          same_intermediate(foreign_committed_history, first),
      "same-numbered history from a foreign writer authority rejects");

  mutation = input;
  ++mutation.trial_flux.revision;
  mutation.momentum.face_flux = mutation.trial_flux.revision;
  mutation.predictor.paired_face_flux = mutation.trial_flux.revision;
  PisoIntermediateCertificate trial_flux_mutated;
  passed &= expect(
      static_cast<bool>(refresh_first(mutation, trial_flux_mutated)) &&
          trial_flux_mutated.r_au == first.r_au &&
          trial_flux_mutated.h_by_a == first.h_by_a &&
          trial_flux_mutated.pressure_face_coefficient ==
              first.pressure_face_coefficient &&
          trial_flux_mutated.phi_h_by_a != first.phi_h_by_a &&
          trial_flux_mutated.trial_face_flux != first.trial_face_flux,
      "trial face-flux revision invalidates phiHbyA only");

  passed &= expect(static_cast<bool>(refresh_first(input, certificate)),
                   "base corrector-one authority can be rebuilt after mutations");
  const PisoIntermediateCertificate authorized = certificate;
  input.corrector = 2U;
  input.prior_corrector = authorized.dependency;
  ConstFieldView revised_velocity = as_const(velocity.view);
  ++revised_velocity.revision;
  input.trial_velocity = revised_velocity;
  ConstFaceFluxView revised_trial_flux = as_const(trial_flux_authority);
  ++revised_trial_flux.revision;
  input.trial_flux = revised_trial_flux;
  input.temporal_reference = {};
  input.committed_face_history = {};
  PisoIntermediateCertificate second;
  passed &= expect(static_cast<bool>(coupler.refresh(input, second)) &&
                       second.valid() && second.corrector == 2U &&
                       second.dependency != first.dependency,
                   "corrector two rebuilds after trial-U revision");

  OwnedField density_accepted =
      make_field(50U, cells, 1U, 0U, 1801U, 2801U);
  OwnedField density_previous =
      make_field(51U, cells, 1U, 0U, 1802U, 2802U);
  OwnedField drho_dp = make_field(6U, cells, 1U, 0U, 1803U, 2803U);
  OwnedField pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 1804U, 2804U);
  OwnedField pressure_rhs =
      make_field(54U, cells, 1U, 0U, 1805U, 2805U);
  fill(density_accepted, 0.95);
  fill(density_previous, 0.7);
  fill(drho_dp, 0.02);
  fill(pressure_diagonal, -17.0);
  fill(pressure_rhs, -19.0);
  PressureCorrectionInput pressure_input;
  pressure_input.intermediate = second;
  pressure_input.pressure_reference = input.pressure_reference;
  pressure_input.density_trial = as_const(density.view);
  pressure_input.density_accepted = as_const(density_accepted.view);
  pressure_input.density_previous = as_const(density_previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = input.bdf;
  pressure_input.time = input.momentum.time;
  pressure_input.geometry = input.momentum.geometry;
  pressure_input.numeric_boundary = input.numeric_boundary;
  const PressureCorrectionSystemView pressure_system{
      pressure_diagonal.view, pressure_rhs.view};
  PressureCorrectionCertificate pressure_certificate;
  passed &= expect(
      static_cast<bool>(coupler.assemble_pressure_system(
          pressure_input, pressure_system, pressure_certificate)) &&
          pressure_certificate.valid() &&
          pressure_certificate.corrector == 2U,
      "pressure system consumes the current corrector-two intermediates");

  // Full pressure-work RED.  At cell (1,1,1), the independent central-
  // difference oracle has V=3/8, rAU_x=3/16, G_x(dp)=2 and
  // G_x(pi)=-52/9.  With U=0 the missing derivative is therefore
  // V*rAU_x*G_x(dp)*G_x(pi)=-0.8125.  The legacy flux-only E_p action
  // returns zero because both its temporal and frozen-h face terms vanish.
  ConstFaceFluxView pressure_energy_target_flux;
  passed &= expect(static_cast<bool>(coupler.inspect_intermediate_flux(
                       second, pressure_energy_target_flux)),
                   "full pressure-work RED inspects the current C2 flux");
  OwnedField pressure_energy_temporal =
      make_field(91U, cells, 1U, 0U, 1812U, 2812U);
  OwnedField pressure_energy_direction =
      make_field(correction_field, cells, 1U, 1U, 1813U, 2813U);
  OwnedField pressure_energy_response =
      make_field(92U, cells, 1U, 0U, 1814U, 2814U);
  OwnedFaceField frozen_h_x =
      make_face(CartesianAxis::x, cells, 2815U);
  OwnedFaceField frozen_h_y =
      make_face(CartesianAxis::y, cells, 2816U);
  OwnedFaceField frozen_h_z =
      make_face(CartesianAxis::z, cells, 2817U);
  fill(pressure_energy_temporal, 0.0);
  fill(pressure_energy_direction, 0.0);
  fill(pressure_energy_response, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        pressure_energy_direction.view.unchecked({x, y, z}, 0U) =
            static_cast<double>(x);
      }
    }
  }
  PressureEnergyFrozenFaceEnthalpy frozen_h{
      as_const(frozen_h_x.view), as_const(frozen_h_y.view),
      as_const(frozen_h_z.view), 1815U, 1816U, 0U};
  frozen_h.local_binding =
      pressure_energy_frozen_face_enthalpy_local_binding(frozen_h);
  PressureEnergyPressureFluxBinding pressure_energy_binding;
  pressure_energy_binding.geometry = &fixture.geometry;
  pressure_energy_binding.boundary = &fixture.boundary;
  pressure_energy_binding.patch = fixture.patch;
  pressure_energy_binding.services = {MPI_COMM_SELF, &correction_halo, 163U,
                                      correction_field, 1U};
  pressure_energy_binding.intermediate = second;
  pressure_energy_binding.pressure = pressure_certificate;
  pressure_energy_binding.temporal_diagonal =
      as_const(pressure_energy_temporal.view);
  pressure_energy_binding.x_pressure_coefficient = as_const(ax.view);
  pressure_energy_binding.y_pressure_coefficient = as_const(ay.view);
  pressure_energy_binding.z_pressure_coefficient = as_const(az.view);
  pressure_energy_binding.target_flux = pressure_energy_target_flux;
  pressure_energy_binding.frozen_face_enthalpy = frozen_h;
  PisoCartesianPressureWorkLinearization pressure_work_linearization;
  passed &= expect(
      static_cast<bool>(coupler.inspect_cartesian_pressure_work_linearization(
          second, pressure_certificate, as_const(target_pressure.view),
          input.trial_velocity, pressure_work_linearization)) &&
          pressure_work_linearization.valid(),
      "coupler privately signs the current Cartesian pressure-work target");
  PisoCartesianPressureWorkLinearization rejected_pressure_work =
      pressure_work_linearization;
  PisoIntermediateCertificate stale_pressure_work_intermediate = second;
  ++stale_pressure_work_intermediate.dependency;
  passed &= expect(
      coupler
                  .inspect_cartesian_pressure_work_linearization(
                      stale_pressure_work_intermediate, pressure_certificate,
                      as_const(target_pressure.view), input.trial_velocity,
                      rejected_pressure_work)
                  .code == StatusCode::invalid_plan &&
          rejected_pressure_work.valid(),
      "stale intermediate cannot replace an issued pressure-work capability");
  PressureCorrectionCertificate stale_pressure_work_system =
      pressure_certificate;
  ++stale_pressure_work_system.state;
  passed &= expect(
      coupler
                  .inspect_cartesian_pressure_work_linearization(
                      second, stale_pressure_work_system,
                      as_const(target_pressure.view), input.trial_velocity,
                      rejected_pressure_work)
                  .code == StatusCode::invalid_plan &&
          rejected_pressure_work.valid(),
      "foreign pressure system cannot replace an issued pressure-work capability");
  const double certified_target_pressure =
      target_pressure.view.unchecked({1, 1, 1}, 0U);
  target_pressure.view.unchecked({1, 1, 1}, 0U) += 1.0;
  passed &= expect(
      coupler
                  .inspect_cartesian_pressure_work_linearization(
                      second, pressure_certificate,
                      as_const(target_pressure.view), input.trial_velocity,
                      rejected_pressure_work)
                  .code == StatusCode::invalid_plan &&
          rejected_pressure_work.valid(),
      "same-revision target-p mutation fails closed against refresh");
  target_pressure.view.unchecked({1, 1, 1}, 0U) =
      certified_target_pressure;
  velocity.view.unchecked({1, 1, 1}, 0U) = 1.0;
  passed &= expect(
      coupler
                  .inspect_cartesian_pressure_work_linearization(
                      second, pressure_certificate,
                      as_const(target_pressure.view), input.trial_velocity,
                      rejected_pressure_work)
                  .code == StatusCode::invalid_plan &&
          rejected_pressure_work.valid(),
      "same-revision target-U mutation fails closed against HbyA");
  velocity.view.unchecked({1, 1, 1}, 0U) = 0.0;
  pressure_energy_binding.pressure_work = pressure_work_linearization;
  pressure_energy_binding.identity = {1911U, 1912U, 1913U, 1914U, 1915U};
  PressureEnergyPressureFluxOperator pressure_energy_operator;
  PressureEnergyPressureFluxCertificate pressure_energy_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          pressure_energy_binding, pressure_energy_operator,
          pressure_energy_certificate)) &&
          static_cast<bool>(pressure_energy_operator.apply(
              pressure_energy_direction.view,
              pressure_energy_response.view)),
      "full pressure-work RED binds and applies the current E_p action");
  passed &= expect(
      pressure_energy_certificate.pressure_work_scope ==
              PressureEnergyPressureWorkScope::exact_cartesian &&
          pressure_energy_certificate.full_cartesian_pressure_work &&
          !pressure_energy_certificate.flux_only_quasi_newton &&
          pressure_energy_certificate.pressure_work_linearization != 0U,
      "exact E_p certificate exposes a typed Cartesian pressure-work scope");
  passed &= expect(
      std::abs(pressure_energy_response.view.unchecked({1, 1, 1}, 0U) +
               0.8125) < 1.0e-13,
      "full E_p matches the independent -0.8125 pressure-work FD oracle");
  OwnedField gauge_shifted_direction =
      make_field(correction_field, cells, 1U, 1U, 1817U, 2818U);
  OwnedField gauge_shifted_response =
      make_field(104U, cells, 1U, 0U, 1818U, 2819U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        gauge_shifted_direction.view.unchecked({x, y, z}, 0U) =
            pressure_energy_direction.view.unchecked({x, y, z}, 0U) + 17.0;
  passed &= expect(
      static_cast<bool>(pressure_energy_operator.apply(
          gauge_shifted_direction.view, gauge_shifted_response.view)),
      "exact spatial E_p applies a constant-gauge-shifted direction");
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        passed &= expect(
            std::abs(gauge_shifted_response.view.unchecked({x, y, z}, 0U) -
                     pressure_energy_response.view.unchecked({x, y, z}, 0U)) <
                1.0e-13,
            "exact spatial E_p is invariant to a constant pressure gauge");
  const Int3 pressure_work_cell{1, 1, 1};
  const double dx = cell_width(fixture, CartesianAxis::x, 1);
  const double base_gradient =
      (target_pressure.view.unchecked({2, 1, 1}, 0U) -
       target_pressure.view.unchecked({0, 1, 1}, 0U)) /
      (2.0 * dx);
  const double direction_gradient =
      (pressure_energy_direction.view.unchecked({2, 1, 1}, 0U) -
       pressure_energy_direction.view.unchecked({0, 1, 1}, 0U)) /
      (2.0 * dx);
  constexpr double pressure_work_epsilon = 1.0e-7;
  const double volume = cell_volume(fixture, pressure_work_cell);
  const double reciprocal = r_au.view.unchecked(pressure_work_cell, 0U);
  const double base_velocity = h_by_a.view.unchecked(pressure_work_cell, 0U);
  const double perturbed_velocity =
      base_velocity - pressure_work_epsilon * reciprocal *
                          direction_gradient;
  const double base_pressure_work =
      -volume * base_velocity * base_gradient;
  const double perturbed_pressure_work =
      -volume * perturbed_velocity *
      (base_gradient + pressure_work_epsilon * direction_gradient);
  const double pressure_work_fd =
      (perturbed_pressure_work - base_pressure_work) /
      pressure_work_epsilon;
  passed &= expect(
      std::abs(pressure_energy_response.view.unchecked(pressure_work_cell,
                                                       0U) -
               pressure_work_fd) < 3.0e-8,
      "exact E_p action matches an independent nonlinear pressure-work FD");
  passed &= expect(
      pressure_energy_operator
                  .apply(pressure_energy_direction.view,
                         target_pressure.view)
                  .code == StatusCode::invalid_plan,
      "exact E_p rejects output aliasing its certified target pressure");

  constexpr FieldId enthalpy_variation_field = 93U;
  constexpr FieldId temperature_variation_field = 94U;
  const std::array<HaloFieldSpec, 2U> energy_halo_fields{{
      {enthalpy_variation_field, 2U, 1U},
      {temperature_variation_field, 1U, 1U},
  }};
  HaloEngine energy_halo;
  passed &= expect(static_cast<bool>(energy_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {energy_halo_fields.data(), energy_halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "exact Cartesian E_h halo reserves both variations");
  OwnedField assembled_enthalpy =
      make_field(95U, cells, 1U, 0U, 1920U, 2820U);
  OwnedField target_enthalpy =
      make_field(3U, cells, 1U, 2U, 1921U, 2821U);
  OwnedField density_enthalpy_derivative =
      make_field(96U, cells, 1U, 0U, 1922U, 2822U);
  OwnedField heat_capacity =
      make_field(97U, cells, 1U, 1U, 1923U, 2823U);
  OwnedField thermal_conductivity =
      make_field(98U, cells, 1U, 1U, 1924U, 2824U);
  OwnedField enthalpy_diffusivity =
      make_field(99U, cells, 1U, 1U, 1925U, 2825U);
  OwnedField delta_temperature = make_field(
      temperature_variation_field, cells, 1U, 1U, 1926U, 2826U);
  fill(assembled_enthalpy, 20.0);
  fill(target_enthalpy, 300000.0);
  fill(density_enthalpy_derivative, -1.0e-6);
  fill(heat_capacity, 1000.0);
  fill(thermal_conductivity, 2.0);
  fill(enthalpy_diffusivity, 0.002);
  fill(delta_temperature, 0.0);
  OwnedFaceField exact_frozen_h_x =
      make_face(CartesianAxis::x, cells, 2827U);
  OwnedFaceField exact_frozen_h_y =
      make_face(CartesianAxis::y, cells, 2827U);
  OwnedFaceField exact_frozen_h_z =
      make_face(CartesianAxis::z, cells, 2827U);
  OwnedFaceField directional_h_x =
      make_face(CartesianAxis::x, cells, 2830U);
  OwnedFaceField directional_h_y =
      make_face(CartesianAxis::y, cells, 2830U);
  OwnedFaceField directional_h_z =
      make_face(CartesianAxis::z, cells, 2830U);
  const FrozenConvectionContext exact_convection_context{
      1930U, fixture.boundary.revision()};
  FrozenConvectionFaceField exact_frozen_h;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          fixture.kernels, fixture.schemes.enthalpy(),
          pressure_energy_target_flux, as_const(target_enthalpy.view), 0U,
          exact_convection_context,
          {exact_frozen_h_x.view, exact_frozen_h_y.view,
           exact_frozen_h_z.view},
          exact_frozen_h)),
      "exact Cartesian E_h freezes the same target flux and enthalpy layer");
  PressureEnergyEnthalpyBinding exact_enthalpy_binding;
  exact_enthalpy_binding.geometry = &fixture.geometry;
  exact_enthalpy_binding.kernels = &fixture.kernels;
  exact_enthalpy_binding.boundary = &fixture.boundary;
  exact_enthalpy_binding.patch = fixture.patch;
  exact_enthalpy_binding.convection = fixture.schemes.enthalpy();
  exact_enthalpy_binding.services = {
      MPI_COMM_SELF, &energy_halo, 1931U, enthalpy_variation_field,
      temperature_variation_field};
  exact_enthalpy_binding.authority = {
      input.bdf,
      input.momentum.time,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision(),
      1932U,
      1933U,
      exact_convection_context.collective_semantics,
      1934U,
      1935U};
  exact_enthalpy_binding.assembled_diagonal =
      as_const(assembled_enthalpy.view);
  exact_enthalpy_binding.target_enthalpy = as_const(target_enthalpy.view);
  exact_enthalpy_binding.density_enthalpy_derivative =
      as_const(density_enthalpy_derivative.view);
  exact_enthalpy_binding.heat_capacity = as_const(heat_capacity.view);
  exact_enthalpy_binding.thermal_conductivity =
      as_const(thermal_conductivity.view);
  exact_enthalpy_binding.enthalpy_diffusivity =
      as_const(enthalpy_diffusivity.view);
  exact_enthalpy_binding.target_flux = pressure_energy_target_flux;
  exact_enthalpy_binding.convection_context = exact_convection_context;
  exact_enthalpy_binding.frozen_face_enthalpy = exact_frozen_h;
  exact_enthalpy_binding.workspace = {
      delta_temperature.view,
      {directional_h_x.view, directional_h_y.view, directional_h_z.view}};
  exact_enthalpy_binding.identity = pressure_energy_binding.identity;
  exact_enthalpy_binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope;
  PressureEnergyEnthalpyOperator exact_enthalpy_operator;
  PressureEnergyEnthalpyCertificate exact_enthalpy_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyEnthalpyOperator::bind(
          exact_enthalpy_binding, exact_enthalpy_operator,
          exact_enthalpy_certificate)) &&
          exact_enthalpy_certificate.valid(),
      "exact Cartesian E_h binds from a real frozen-spatial operator");

  PressureEnergySchurBlockAuthority exact_block_authority;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurBlockAuthority::exact_cartesian(
          pressure_energy_operator, exact_enthalpy_operator,
          exact_block_authority)) &&
          exact_block_authority.valid() &&
          exact_block_authority.scope() == PressureEnergySchurBlockScope::
                                               exact_cartesian_frozen_spatial,
      "real complete E_p and E_h operators issue exact Cartesian authority");

  PressureEnergyPressureFluxBinding flux_only_binding =
      pressure_energy_binding;
  flux_only_binding.pressure_work = {};
  PressureEnergyPressureFluxOperator flux_only_operator;
  PressureEnergyPressureFluxCertificate flux_only_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          flux_only_binding, flux_only_operator, flux_only_certificate)) &&
          flux_only_certificate.valid() &&
          flux_only_certificate.pressure_work_scope ==
              PressureEnergyPressureWorkScope::flux_only_quasi_newton &&
          PressureEnergySchurBlockAuthority::exact_cartesian(
              flux_only_operator, exact_enthalpy_operator,
              exact_block_authority)
                  .code == StatusCode::invalid_plan &&
          exact_block_authority.valid(),
      "a real legacy flux-only E_p operator cannot promote exact scope");

  LinearOperatorCertificate continuity_pressure_certificate =
      pressure_energy_operator.certificate();
  continuity_pressure_certificate.collective_fingerprint ^=
      UINT64_C(0x3333333333333333);
  if (continuity_pressure_certificate.collective_fingerprint == 0U)
    continuity_pressure_certificate.collective_fingerprint = 1U;
  continuity_pressure_certificate.operator_class = LinearOperatorClass::spd;
  IdentityCertificateOperator continuity_pressure_operator(
      continuity_pressure_certificate);
  OwnedField continuity_enthalpy_diagonal =
      make_field(100U, cells, 1U, 0U, 1936U, 2836U);
  OwnedField continuity_enthalpy_scale =
      make_field(101U, cells, 1U, 0U, 1937U, 2837U);
  OwnedField continuity_workspace =
      make_field(102U, cells, 1U, 0U, 1938U, 2838U);
  OwnedField eliminated_workspace = make_field(
      enthalpy_variation_field, cells, 1U, 2U, 1939U, 2839U);
  OwnedField energy_workspace =
      make_field(103U, cells, 1U, 0U, 1940U, 2840U);
  fill(continuity_enthalpy_diagonal, -2.0);
  fill(continuity_enthalpy_scale, 2.0);
  PressureEnergySchurBinding exact_schur_binding;
  exact_schur_binding.continuity_pressure = &continuity_pressure_operator;
  exact_schur_binding.energy_pressure = &pressure_energy_operator;
  exact_schur_binding.energy_enthalpy = &exact_enthalpy_operator;
  exact_schur_binding.continuity_enthalpy_diagonal =
      as_const(continuity_enthalpy_diagonal.view);
  exact_schur_binding.continuity_enthalpy_row_scale =
      as_const(continuity_enthalpy_scale.view);
  exact_schur_binding.workspace = {
      continuity_workspace.view, eliminated_workspace.view,
      energy_workspace.view};
  exact_schur_binding.scaled_pivot_floor = 1.0e-12;
  PressureEnergySchurOperator missing_scope_schur;
  PressureEnergyJacobianCertificate missing_scope_certificate;
  passed &= expect(
      PressureEnergySchurOperator::bind(
          exact_schur_binding, missing_scope_schur,
          missing_scope_certificate)
                  .code == StatusCode::invalid_plan &&
          !missing_scope_certificate.valid(),
      "Cartesian Schur fails closed without typed block authority");
  exact_schur_binding.block_authority = exact_block_authority;
  PressureEnergySchurOperator exact_schur;
  PressureEnergyJacobianCertificate exact_schur_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurOperator::bind(
          exact_schur_binding, exact_schur, exact_schur_certificate)) &&
          exact_schur_certificate.valid() &&
          exact_schur_certificate.jacobian_scope ==
              PressureEnergyJacobianScope::exact_cartesian_frozen_spatial &&
          !exact_schur_certificate.full_nonlinear_jacobian,
      "Schur derives exact Cartesian scope only from real typed operators");
  PressureEnergySchurBinding foreign_component_binding = exact_schur_binding;
  foreign_component_binding.energy_pressure = &flux_only_operator;
  PressureEnergySchurOperator foreign_component_schur;
  PressureEnergyJacobianCertificate foreign_component_certificate =
      exact_schur_certificate;
  passed &= expect(
      PressureEnergySchurOperator::bind(
          foreign_component_binding, foreign_component_schur,
          foreign_component_certificate)
                  .code == StatusCode::invalid_plan &&
          !foreign_component_certificate.valid(),
      "typed Cartesian authority rejects a foreign energy-pressure block");
  if (!passed)
    return false;

  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const double volume = cell_volume(fixture, cell);
        const double expected_rhs = 0.75 * volume;
        const double face_sum =
            ax.view.unchecked(cell) +
            ax.view.unchecked({xc + 1, yc, zc}) +
            ay.view.unchecked(cell) +
            ay.view.unchecked({xc, yc + 1, zc}) +
            az.view.unchecked(cell) +
            az.view.unchecked({xc, yc, zc + 1});
        passed &= expect(
            std::abs(pressure_diagonal.view.unchecked(cell, 0U) -
                     (face_sum + 10.0 * 0.02 * volume)) < 1.0e-14 &&
                std::abs(pressure_rhs.view.unchecked(cell, 0U) -
                         expected_rhs) < 1.0e-14,
            "pressure matrix contains full BDF defect and a0*V*drho_dp once");
      }
    }
  }
  OwnedField pressure_solution =
      make_field(60U, cells, 1U, 1U, 1810U, 2810U);
  OwnedField pressure_applied =
      make_field(61U, cells, 1U, 0U, 1811U, 2811U);
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        pressure_solution.view.unchecked({xc, yc, zc}, 0U) =
            static_cast<double>(xc + 1);
      }
    }
  }
  const std::array<HaloFieldSpec, 1U> pressure_halo_fields{{
      {pressure_solution.view.field, 1U, 1U}}};
  HaloEngine pressure_halo;
  passed &= expect(static_cast<bool>(pressure_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {pressure_halo_fields.data(),
                        pressure_halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "pressure Krylov halo reserves cold");
  PressureLinearOperator pressure_operator;
  const PressureOperatorServices pressure_services{
      MPI_COMM_SELF, &pressure_halo, 161U, pressure_solution.view.field};
  passed &= expect(
      static_cast<bool>(coupler.bind_pressure_operator(
          pressure_services, pressure_system, pressure_operator)) &&
          pressure_operator.fingerprint() != 0U &&
          pressure_operator.coefficient_storage_address() ==
              reinterpret_cast<std::uintptr_t>(pressure_diagonal.view.base),
      "pressure operator persistently binds current coefficient storage");
  const LinearIdentity pressure_identity{1901U, 1902U, 1903U, 1904U,
                                         1905U};
  passed &= expect(
      static_cast<bool>(pressure_operator.refresh(
          {pressure_certificate, pressure_identity, 1906U})) &&
          pressure_operator.certificate().operator_class ==
              LinearOperatorClass::spd,
      "pressure operator refreshes numeric identity without rebinding");
  passed &= expect(static_cast<bool>(pressure_operator.apply(
                       pressure_solution.view, pressure_applied.view)),
                   "pressure operator applies through persistent halo");
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const double centre = static_cast<double>(xc + 1);
        const double left = static_cast<double>((xc + cells.x - 1) % cells.x +
                                                1);
        const double right =
            static_cast<double>((xc + 1) % cells.x + 1);
        const double expected =
            pressure_diagonal.view.unchecked(cell, 0U) * centre -
            ax.view.unchecked(cell) * left -
            ax.view.unchecked({xc + 1, yc, zc}) * right -
            (ay.view.unchecked(cell) +
             ay.view.unchecked({xc, yc + 1, zc}) +
             az.view.unchecked(cell) +
             az.view.unchecked({xc, yc, zc + 1})) *
                centre;
        passed &= expect(
            std::abs(pressure_applied.view.unchecked(cell, 0U) - expected) <
                1.0e-14,
            "pressure operator matches independent periodic stencil oracle");
      }
    }
  }
  const LinearOperatorCertificate operator_marker =
      pressure_operator.certificate();
  PressureOperatorRevision stale_operator{
      pressure_certificate, pressure_identity, 1906U};
  stale_operator.pressure.plan = 0U;
  passed &= expect(
      pressure_operator.refresh(stale_operator).code ==
              StatusCode::invalid_plan &&
          pressure_operator.certificate().identity.fingerprint ==
              operator_marker.identity.fingerprint,
      "stale pressure certificate cannot replace published operator identity");

  PressureCorrectionCertificate pressure_mutated;
  ++pressure_input.density_accepted.revision;
  passed &= expect(
      static_cast<bool>(coupler.assemble_pressure_system(
          pressure_input, pressure_system, pressure_mutated)) &&
          pressure_mutated.state != pressure_certificate.state,
      "accepted-density revision invalidates pressure numeric identity");
  --pressure_input.density_accepted.revision;

  const double saved_drho = drho_dp.view.unchecked({1, 1, 1}, 0U);
  drho_dp.view.unchecked({1, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  fill(pressure_diagonal, -23.0);
  fill(pressure_rhs, -29.0);
  const PressureCorrectionCertificate pressure_marker = pressure_certificate;
  pressure_mutated = pressure_marker;
  passed &= expect(
      coupler
                  .assemble_pressure_system(pressure_input, pressure_system,
                                            pressure_mutated)
                  .code == StatusCode::numerical_failure &&
          pressure_diagonal.view.unchecked({1, 1, 1}, 0U) == -23.0 &&
          pressure_rhs.view.unchecked({1, 1, 1}, 0U) == -29.0 &&
          pressure_mutated.state == pressure_marker.state,
      "non-finite drho/dp rejects atomically before pressure publication");
  drho_dp.view.unchecked({1, 1, 1}, 0U) = saved_drho;

  PressureCorrectionSystemView pressure_aliased = pressure_system;
  pressure_aliased.diagonal = density_accepted.view;
  passed &= expect(
      coupler
                  .assemble_pressure_system(pressure_input, pressure_aliased,
                                            pressure_mutated)
                  .code == StatusCode::invalid_plan,
      "pressure assembly rejects input/output aliasing");

  passed &= expect(coupler.refresh(input, certificate).code ==
                       StatusCode::invalid_plan,
                   "consumed corrector-one token cannot authorize a third solve");

  input.trial_flux = as_const(trial_flux_authority);
  input.temporal_reference = as_const(phi_h_by_a);
  input.committed_face_history = {accepted_flux, previous_flux};
  passed &= expect(static_cast<bool>(refresh_first(input, certificate)),
                   "new corrector-one token is available for failure poisoning");
  input.corrector = 2U;
  input.temporal_reference = {};
  input.committed_face_history = {};
  input.prior_corrector = certificate.dependency;
  revised_trial_flux = as_const(trial_flux_authority);
  revised_trial_flux.revision = certificate.trial_face_flux;
  input.trial_flux = revised_trial_flux;
  PisoIntermediateCertificate stale_flux = second;
  passed &= expect(coupler.refresh(input, stale_flux).code ==
                       StatusCode::invalid_plan &&
                       same_intermediate(stale_flux, second),
                   "corrector two rejects the stale corrector-one trial flux");

  input.trial_flux = as_const(trial_flux_authority);
  input.temporal_reference = as_const(phi_h_by_a);
  input.committed_face_history = {accepted_flux, previous_flux};
  passed &= expect(static_cast<bool>(refresh_first(input, certificate)),
                   "fresh corrector-one token is available after stale flux rejection");
  input.corrector = 2U;
  input.temporal_reference = {};
  input.committed_face_history = {};
  input.prior_corrector = certificate.dependency;
  revised_trial_flux.revision = certificate.trial_face_flux + 1U;
  input.trial_flux = revised_trial_flux;
  const double corrector_two_saved =
      diagonal.view.unchecked({1, 1, 1}, 0U);
  diagonal.view.unchecked({1, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  PisoIntermediateCertificate failed_second = second;
  passed &= expect(coupler.refresh(input, failed_second).code ==
                       StatusCode::numerical_failure &&
                       same_intermediate(failed_second, second),
                   "failed corrector two publishes no intermediate certificate");
  diagonal.view.unchecked({1, 1, 1}, 0U) = corrector_two_saved;
  passed &= expect(coupler.refresh(input, failed_second).code ==
                       StatusCode::invalid_plan,
                   "failed corrector two consumes its one-shot authorization");

  PisoCouplerWorkspace aliased = workspace;
  aliased.h_by_a = r_au.view;
  PressureVelocityCoupler rejected;
  passed &= expect(PressureVelocityCoupler::bind(
                       piso, fixture.equations, services, aliased, rejected)
                           .code == StatusCode::invalid_plan &&
                       rejected.fingerprint() == 0U,
                   "coupler rejects cross-workspace aliasing");

  input.corrector = 1U;
  input.prior_corrector = 0U;
  input.trial_velocity = as_const(velocity.view);
  input.trial_flux = as_const(trial_flux_authority);
  input.temporal_reference = as_const(phi_h_by_a);
  input.committed_face_history = {accepted_flux, previous_flux};
  const double saved = diagonal.view.unchecked({1, 1, 1}, 0U);
  diagonal.view.unchecked({1, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  const PisoIntermediateCertificate marker = second;
  certificate = marker;
  passed &= expect(coupler.refresh(input, certificate).code ==
                       StatusCode::numerical_failure &&
                       same_intermediate(certificate, marker),
                   "non-finite momentum diagonal rejects before certificate");
  diagonal.view.unchecked({1, 1, 1}, 0U) = saved;

  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        velocity.view.unchecked({xc, yc, zc}, 0U) =
            static_cast<double>(xc + 1);
      }
    }
  }
  input.corrector = 1U;
  input.prior_corrector = 0U;
  input.trial_velocity = as_const(velocity.view);
  fill_face_flux(phi_h_by_a, 0.0);
  fill_face_flux(trial_flux_authority, 0.0);
  passed &= expect(static_cast<bool>(coupler.refresh(input, certificate)),
                   "periodic halo rebuild succeeds for nonuniform HbyA");
  const double x_face_area = face_area(fixture, CartesianAxis::x, {0, 1, 1});
  const double periodic_face_flux = 2.5 * x_face_area;
  passed &= expect(
      std::abs(phi_h_by_a.x.unchecked({0, 1, 1}) -
                   periodic_face_flux) < 1.0e-14 &&
          std::abs(phi_h_by_a.x.unchecked({cells.x, 1, 1}) -
                   periodic_face_flux) < 1.0e-14 &&
          std::abs(phi_h_by_a.x.unchecked({1, 1, 1}) - 1.5 * x_face_area) <
              1.0e-14,
      "periodic face flux interpolates opposite ghosts rather than clamping");
  const double x_width = cell_width(fixture, CartesianAxis::x, 0);
  const double left_volume = cell_volume(fixture, {cells.x - 1, 1, 1});
  const double right_volume = cell_volume(fixture, {0, 1, 1});
  const double expected_periodic_x_coefficient =
      x_face_area /
      (0.5 * x_width / (left_volume / 2.0) +
       0.5 * x_width / (right_volume / 2.0));
  passed &= expect(
      std::abs(ax.view.unchecked({0, 1, 1}) -
               expected_periodic_x_coefficient) < 1.0e-14 &&
          std::abs(ax.view.unchecked({cells.x, 1, 1}) -
                   expected_periodic_x_coefficient) < 1.0e-14,
      "periodic pressure coefficients use harmonic rho*rAU and wrapped centre distance");

  fill(density_accepted, 0.9);
  ++density_accepted.view.revision;
  pressure_input.intermediate = certificate;
  pressure_input.density_trial = as_const(density.view);
  pressure_input.density_accepted = as_const(density_accepted.view);
  pressure_input.density_previous = as_const(density_previous.view);
  pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  pressure_input.bdf = input.bdf;
  pressure_input.time = input.momentum.time;
  pressure_input.geometry = input.momentum.geometry;
  pressure_input.numeric_boundary = input.numeric_boundary;
  PressureCorrectionCertificate solve_pressure_one;
  passed &= expect(static_cast<bool>(coupler.assemble_pressure_system(
                       pressure_input, pressure_system,
                       solve_pressure_one)) &&
                       solve_pressure_one.corrector == 1U,
                   "corrector-one pressure system prepares the solve epoch");

  LinearWorkspaceRequirements krylov_requirements;
  passed &= expect(static_cast<bool>(make_linear_workspace_requirements(
                       LinearAlgorithm::fgmres, cells, 1U,
                       piso.pressure_solve().restart,
                       ReductionMode::mpi_allreduce, 2001U,
                       krylov_requirements)),
                   "pressure FGMRES workspace requirements compile");
  constexpr StorageIdentity krylov_storage = 2901U;
  OwnedField krylov_vectors = make_field(
      70U, cells, krylov_requirements.vector_slots, 1U, 2002U,
      krylov_storage);
  OwnedField krylov_scalars = make_field(
      71U,
      {static_cast<std::int32_t>(krylov_requirements.scalar_doubles), 1, 1},
      1U, 0U, 2003U, krylov_storage);
  SolverWorkspace krylov_workspace;
  passed &= expect(static_cast<bool>(SolverWorkspace::bind(
                       krylov_requirements, krylov_vectors.view,
                       krylov_scalars.view, krylov_workspace)),
                   "pressure FGMRES workspace binds once");
  constexpr StorageIdentity coupled_krylov_storage = 2930U;
  OwnedField coupled_krylov_vectors = make_field(
      70U, cells, krylov_requirements.vector_slots, 1U, 2015U,
      coupled_krylov_storage);
  OwnedField coupled_krylov_scalars = make_field(
      71U,
      {static_cast<std::int32_t>(krylov_requirements.scalar_doubles), 1, 1},
      1U, 0U, 2016U, coupled_krylov_storage);
  SolverWorkspace coupled_krylov_workspace;
  passed &= expect(static_cast<bool>(SolverWorkspace::bind(
                       krylov_requirements, coupled_krylov_vectors.view,
                       coupled_krylov_scalars.view,
                       coupled_krylov_workspace)),
                   "coupled pressure-energy FGMRES workspace binds independently");
  ReductionEngine reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_SELF, ReductionMode::mpi_allreduce,
                       krylov_requirements.reduction_capacity, reductions)),
                   "pressure reduction engine compiles once");

  // The coupled pressure/enthalpy block owns positivity as one joint
  // candidate.  Its pressure correction therefore must not be silently
  // rescaled by the legacy pressure-only 0.9 density-depletion limiter.
  OwnedField probe_r_au =
      make_field(40U, cells, 3U, 1U, 2004U, 2911U);
  OwnedField probe_h_by_a =
      make_field(41U, cells, 3U, reach, 2005U, 2912U);
  OwnedField probe_gradient =
      make_field(42U, cells, 3U, 0U, 2006U, 2913U);
  OwnedFaceField probe_ax = make_face(CartesianAxis::x, cells, 2914U);
  OwnedFaceField probe_ay = make_face(CartesianAxis::y, cells, 2915U);
  OwnedFaceField probe_az = make_face(CartesianAxis::z, cells, 2916U);
  FaceFluxStorage probe_phi_storage;
  FaceFluxView probe_phi;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, probe_phi_storage)) &&
          static_cast<bool>(probe_phi_storage.workspace_view(
              0U, 2007U, probe_phi)),
      "coupled correction scale probe workspace allocates");
  const PisoCouplerWorkspace probe_workspace{
      probe_r_au.view, probe_h_by_a.view, probe_gradient.view,
      probe_ax.view, probe_ay.view, probe_az.view, probe_phi};
  PressureVelocityCoupler coupled_scale_probe;
  passed &= expect(static_cast<bool>(PressureVelocityCoupler::bind(
                       piso, fixture.equations, services, probe_workspace,
                       coupled_scale_probe)),
                   "coupled correction scale probe binds independently");
  OwnedField probe_density =
      make_field(0U, cells, 1U, reach, 2008U, 2917U);
  OwnedField probe_velocity =
      make_field(1U, cells, 3U, 0U, 2009U, 2918U);
  fill(probe_density, 1.0);
  fill(probe_velocity, 0.0);
  PisoIntermediateInput probe_input = input;
  constexpr double probe_pressure_reference_value = 101325.0;
  probe_input.thermophysical_boundary.binding.pressure_reference =
      probe_pressure_reference_value;
  probe_input.density = probe_density.view;
  probe_input.trial_velocity = as_const(probe_velocity.view);
  probe_input.temporal_reference = as_const(probe_phi);
  probe_input.predictor.predicted_density = probe_density.view.revision;
  probe_input.predictor.predicted_density_storage =
      probe_density.view.storage_identity;
  probe_input.predictor.predicted_density_revision_domain =
      probe_density.view.revision_domain;
  PisoIntermediateCertificate probe_intermediate;
  passed &= expect(static_cast<bool>(coupled_scale_probe.refresh(
                       probe_input, probe_intermediate)),
                   "coupled correction scale probe refreshes C1");
  OwnedField probe_pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 2010U, 2919U);
  OwnedField probe_pressure_rhs =
      make_field(54U, cells, 1U, 0U, 2011U, 2920U);
  PressureCorrectionInput probe_pressure_input = pressure_input;
  probe_pressure_input.intermediate = probe_intermediate;
  probe_pressure_input.density_trial = as_const(probe_density.view);
  const PressureCorrectionSystemView probe_pressure_system{
      probe_pressure_diagonal.view, probe_pressure_rhs.view};
  PressureCorrectionCertificate probe_pressure;
  passed &= expect(static_cast<bool>(coupled_scale_probe.assemble_pressure_system(
                       probe_pressure_input, probe_pressure_system,
                       probe_pressure)),
                   "coupled correction scale probe assembles C1 pressure");
  ConstFaceFluxView probe_pressure_energy_target_flux;
  passed &= expect(
      static_cast<bool>(coupled_scale_probe.inspect_intermediate_flux(
          probe_intermediate, probe_pressure_energy_target_flux)),
      "correction-entry stale probe inspects its current target flux");
  PisoCartesianPressureWorkLinearization probe_pressure_work;
  passed &= expect(
      static_cast<bool>(
          coupled_scale_probe.inspect_cartesian_pressure_work_linearization(
              probe_intermediate, probe_pressure,
              as_const(target_pressure.view),
              as_const(probe_velocity.view), probe_pressure_work)) &&
          probe_pressure_work.valid(),
      "correction-entry stale probe issues a live pressure-work capability");
  PressureEnergyPressureFluxBinding probe_energy_pressure_binding;
  probe_energy_pressure_binding.geometry = &fixture.geometry;
  probe_energy_pressure_binding.boundary = &fixture.boundary;
  probe_energy_pressure_binding.patch = fixture.patch;
  probe_energy_pressure_binding.services = {
      MPI_COMM_SELF, &correction_halo, 163U, correction_field, 1U};
  probe_energy_pressure_binding.intermediate = probe_intermediate;
  probe_energy_pressure_binding.pressure = probe_pressure;
  probe_energy_pressure_binding.temporal_diagonal =
      as_const(pressure_energy_temporal.view);
  probe_energy_pressure_binding.x_pressure_coefficient =
      as_const(probe_ax.view);
  probe_energy_pressure_binding.y_pressure_coefficient =
      as_const(probe_ay.view);
  probe_energy_pressure_binding.z_pressure_coefficient =
      as_const(probe_az.view);
  probe_energy_pressure_binding.target_flux =
      probe_pressure_energy_target_flux;
  probe_energy_pressure_binding.frozen_face_enthalpy = frozen_h;
  probe_energy_pressure_binding.identity = pressure_energy_binding.identity;
  probe_energy_pressure_binding.pressure_work = probe_pressure_work;
  PressureEnergyPressureFluxOperator probe_energy_pressure_operator;
  PressureEnergyPressureFluxCertificate probe_energy_pressure_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          probe_energy_pressure_binding, probe_energy_pressure_operator,
          probe_energy_pressure_certificate)),
      "correction-entry stale probe binds a real complete E_p operator");
  PressureEnergySchurBlockAuthority probe_scope_authority;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurBlockAuthority::exact_cartesian(
          probe_energy_pressure_operator, exact_enthalpy_operator,
          probe_scope_authority)) &&
          probe_scope_authority.valid(),
      "live correction-entry E_p can sign exact scope before state mutation");
  OwnedField probe_correction =
      make_field(correction_field, cells, 1U, 1U, 2012U, 2921U);
  OwnedField probe_pressure_perturbation =
      make_field(2U, cells, 1U, 1U, 2013U, 2922U);
  OwnedField probe_enthalpy =
      make_field(3U, cells, 1U, 1U, 2017U, 2923U);
  OwnedField probe_temperature =
      make_field(4U, cells, 1U, 1U, 2018U, 2924U);
  OwnedField probe_enthalpy_correction =
      make_field(63U, cells, 1U, 0U, 2019U, 2925U);
  OwnedField probe_candidate_enthalpy =
      make_field(64U, cells, 1U, 0U, 2020U, 2926U);
  OwnedField probe_candidate_density =
      make_field(65U, cells, 1U, 0U, 2021U, 2927U);
  OwnedField probe_candidate_temperature =
      make_field(66U, cells, 1U, 0U, 2022U, 2928U);
  OwnedField probe_candidate_compressibility =
      make_field(67U, cells, 1U, 0U, 2023U, 2929U);
  // The affine pressure-only update is 1 + 0.02*(-75) = -0.5.  The
  // caller-certified EOS candidate remains positive and is the only density
  // authority the coupled transaction may publish.
  double probe_base_enthalpy = 0.0;
  double probe_candidate_enthalpy_value = 0.0;
  double probe_cp = 0.0;
  double probe_gas_constant = 0.0;
  passed &= expect(
      static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
          300.0, {}, probe_base_enthalpy, probe_cp,
          probe_gas_constant)) &&
          static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
              310.0, {}, probe_candidate_enthalpy_value, probe_cp,
              probe_gas_constant)),
      "coupled correction probe constructs physical base/candidate enthalpy");
  fill(probe_correction, -75.0);
  fill(probe_pressure_perturbation, 0.0);
  fill(probe_enthalpy, probe_base_enthalpy);
  fill(probe_temperature, 300.0);
  fill(probe_enthalpy_correction,
       probe_candidate_enthalpy_value - probe_base_enthalpy);
  fill(probe_candidate_enthalpy, probe_candidate_enthalpy_value);
  fill(probe_candidate_density, 0.0);
  fill(probe_candidate_temperature, 0.0);
  fill(probe_candidate_compressibility, 0.0);
  FaceFluxStorage probe_flux_storage;
  FaceFluxView probe_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, probe_flux_storage)) &&
          static_cast<bool>(probe_flux_storage.workspace_view(
              0U, 2014U, probe_flux)),
      "coupled correction scale probe flux allocates");
  const PisoCoupledStateView probe_state{
      probe_velocity.view, probe_pressure_perturbation.view,
      probe_enthalpy.view, probe_density.view, probe_temperature.view};
  PisoExactEosClosureIdentity probe_closure;
  probe_closure.thermodynamics =
      probe_pressure_input.pressure_reference.thermodynamics;
  probe_closure.pressure_reference =
      probe_pressure_input.pressure_reference.pressure_reference;
  probe_closure.composition = exact_composition_identity_for_test(
      probe_closure.thermodynamics, {}, cells);
  probe_closure.pressure_state = make_piso_field_revision_identity(
      as_const(probe_pressure_perturbation.view));
  probe_closure.pressure_correction = make_piso_field_revision_identity(
      as_const(probe_correction.view));
  probe_closure.enthalpy_state =
      make_piso_field_revision_identity(as_const(probe_enthalpy.view));
  probe_closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(probe_enthalpy_correction.view));
  probe_closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(probe_candidate_enthalpy.view));
  probe_closure.candidate_density = make_piso_field_revision_identity(
      as_const(probe_candidate_density.view));
  probe_closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(probe_candidate_temperature.view));
  probe_closure.closure = 2930U;
  PisoExactThermodynamicCandidateView probe_candidate;
  probe_candidate.enthalpy = as_const(probe_candidate_enthalpy.view);
  probe_candidate.density = as_const(probe_candidate_density.view);
  probe_candidate.temperature = as_const(probe_candidate_temperature.view);
  probe_candidate.closure = probe_closure;
  probe_candidate.pressure_compressibility =
      as_const(probe_candidate_compressibility.view);
  probe_candidate.pressure_compressibility.field = drho_dp.view.field;
  passed &= expect(static_cast<bool>(fill_exact_thermodynamic_candidate(
                       fixture.thermodynamics, probe_pressure_reference_value,
                       as_const(probe_pressure_perturbation.view),
                       as_const(probe_correction.view),
                       as_const(probe_candidate_enthalpy.view),
                       as_const(probe_velocity.view), {},
                       probe_candidate_density.view,
                       probe_candidate_temperature.view,
                       probe_candidate_compressibility.view)),
                   "coupled correction probe replays its exact EOS state");
  const double probe_expected_density =
      probe_candidate_density.view.unchecked({1, 1, 1}, 0U);
  const double probe_expected_temperature =
      probe_candidate_temperature.view.unchecked({1, 1, 1}, 0U);
  passed &= expect(
      static_cast<bool>(prepare_closed_gauge(
          fixture.equations.pressure_reference(), probe_pressure,
          probe_pressure_input.pressure_reference,
          probe_pressure_reference_value,
          as_const(probe_pressure_perturbation.view),
          as_const(probe_correction.view),
          probe_candidate.pressure_compressibility,
          probe_candidate.closure.closure, reductions,
          probe_candidate.closed_gauge)),
      "coupled correction gauge transaction certifies before publication");
  const std::vector<double> probe_base_flux =
      face_flux_values(as_const(probe_phi));
  PisoStateCorrectionCertificate probe_correction_certificate;
  passed &= expect(
      static_cast<bool>(coupled_scale_probe.correct_coupled_trial_state(
          probe_pressure, probe_correction.view,
          as_const(probe_enthalpy_correction.view), probe_state,
          probe_candidate, probe_flux, reductions,
          probe_correction_certificate)) &&
          std::abs(probe_correction.view.unchecked({1, 1, 1}, 0U) + 75.0) <
              1.0e-14 &&
          std::abs(probe_candidate.closed_gauge.shift + 75.0) < 1.0e-14 &&
          std::abs(probe_candidate.closed_gauge.next_pressure_reference -
                   101250.0) < 1.0e-14 &&
          std::abs(probe_pressure_perturbation.view.unchecked({1, 1, 1},
                                                               0U)) <
              1.0e-14 &&
          std::abs(probe_candidate.closed_gauge.next_pressure_reference +
                       probe_pressure_perturbation.view.unchecked(
                           {1, 1, 1}, 0U) -
                   101250.0) < 1.0e-14 &&
          std::abs(probe_enthalpy.view.unchecked({1, 1, 1}, 0U) -
                   probe_candidate_enthalpy_value) <
              1.0e-14 &&
          std::abs(probe_density.view.unchecked({1, 1, 1}, 0U) -
                   probe_expected_density) <
              1.0e-14 &&
          std::abs(probe_temperature.view.unchecked({1, 1, 1}, 0U) -
                   probe_expected_temperature) <
              1.0e-14 &&
          std::abs(probe_velocity.view.unchecked({1, 1, 1}, 0U)) <
              1.0e-14 &&
          face_flux_values(as_const(probe_flux)) == probe_base_flux &&
          probe_correction_certificate.closure ==
              PisoStateClosure::exact_eos &&
          probe_correction_certificate.input_pressure_reference
                  .pressure_reference ==
              probe_pressure_input.pressure_reference.pressure_reference &&
          probe_correction_certificate.output_pressure_reference
                  .pressure_reference ==
              probe_candidate.closed_gauge.rank_local_transaction &&
          probe_correction_certificate.closed_gauge_collective_transaction ==
              probe_candidate.closed_gauge.collective_transaction &&
          probe_correction_certificate.enthalpy ==
              probe_enthalpy.view.revision &&
          probe_correction_certificate.temperature ==
              probe_temperature.view.revision &&
          probe_correction_certificate.exact_eos_closure != 0U,
      "constant raw dp moves only the closed pressure reference while p-absolute/EOS/U/flux remain invariant");
  passed &= expect(
      PressureEnergySchurBlockAuthority::exact_cartesian(
          probe_energy_pressure_operator, exact_enthalpy_operator,
          probe_scope_authority)
                  .code == StatusCode::invalid_plan &&
          probe_scope_authority.valid(),
      "correction entry invalidates the old E_p capability before any new "
      "exact Schur scope can be signed");
  PressureEnergySchurBinding stale_signed_scope_binding = exact_schur_binding;
  stale_signed_scope_binding.energy_pressure =
      &probe_energy_pressure_operator;
  stale_signed_scope_binding.block_authority = probe_scope_authority;
  PressureEnergySchurOperator stale_signed_scope_schur;
  PressureEnergyJacobianCertificate stale_signed_scope_certificate =
      exact_schur_certificate;
  passed &= expect(
      PressureEnergySchurOperator::bind(
          stale_signed_scope_binding, stale_signed_scope_schur,
          stale_signed_scope_certificate)
                  .code == StatusCode::invalid_plan &&
          !stale_signed_scope_certificate.valid(),
      "a previously signed exact authority cannot bind Schur after its E_p "
      "issuer is stale");

  auto reset_probe_state = [&]() {
    probe_input.thermophysical_boundary.binding.pressure_reference =
        probe_pressure_reference_value;
    fill(probe_correction, -75.0);
    fill(probe_pressure_perturbation, 0.0);
    fill(probe_enthalpy, probe_base_enthalpy);
    fill(probe_density, 1.0);
    fill(probe_temperature, 300.0);
    fill(probe_enthalpy_correction,
         probe_candidate_enthalpy_value - probe_base_enthalpy);
    fill(probe_candidate_enthalpy, probe_candidate_enthalpy_value);
    fill(probe_candidate_density, 0.0);
    fill(probe_candidate_temperature, 0.0);
    fill(probe_candidate_compressibility, 0.0);
    fill(probe_velocity, 0.0);
    passed &= expect(static_cast<bool>(fill_exact_thermodynamic_candidate(
                         fixture.thermodynamics,
                         probe_pressure_reference_value,
                         as_const(probe_pressure_perturbation.view),
                         as_const(probe_correction.view),
                         as_const(probe_candidate_enthalpy.view),
                         as_const(probe_velocity.view), {},
                         probe_candidate_density.view,
                         probe_candidate_temperature.view,
                         probe_candidate_compressibility.view)),
                     "probe reset restores a physical exact EOS candidate");
    fill_face_flux(probe_flux, 7.0);
  };
  auto rebuild_probe_c1 = [&]() {
    PisoIntermediateCertificate rebuilt;
    probe_input.corrector = 1U;
    probe_input.prior_corrector = 0U;
    Status rebuilt_status =
        coupled_scale_probe.refresh(probe_input, rebuilt);
    if (rebuilt_status) {
      probe_pressure_input.intermediate = rebuilt;
      rebuilt_status = coupled_scale_probe.assemble_pressure_system(
          probe_pressure_input, probe_pressure_system, probe_pressure);
    }
    return rebuilt_status;
  };
  auto prepare_probe_gauge =
      [&](PisoExactThermodynamicCandidateView& selected) {
        selected = probe_candidate;
        selected.closed_gauge = {};
        return prepare_closed_gauge(
            fixture.equations.pressure_reference(), probe_pressure,
            probe_pressure_input.pressure_reference,
            probe_pressure_reference_value,
            as_const(probe_pressure_perturbation.view),
            as_const(probe_correction.view),
            selected.pressure_compressibility,
            selected.closure.closure, reductions, selected.closed_gauge);
      };
  auto atomically_unchanged = [&](const std::vector<double>& velocity_before,
                                  const std::vector<double>& pressure_before,
                                  const std::vector<double>& enthalpy_before,
                                  const std::vector<double>& density_before,
                                  const std::vector<double>& temperature_before,
                                  const std::vector<double>& flux_before) {
    return probe_velocity.bytes == velocity_before &&
           probe_pressure_perturbation.bytes == pressure_before &&
           probe_enthalpy.bytes == enthalpy_before &&
           probe_density.bytes == density_before &&
           probe_temperature.bytes == temperature_before &&
           face_flux_values(as_const(probe_flux)) == flux_before;
  };
  auto reject_probe_candidate =
      [&](PisoExactThermodynamicCandidateView selected,
          StatusCode expected_code) {
        const std::vector<double> velocity_before = probe_velocity.bytes;
        const std::vector<double> pressure_before =
            probe_pressure_perturbation.bytes;
        const std::vector<double> enthalpy_before = probe_enthalpy.bytes;
        const std::vector<double> density_before = probe_density.bytes;
        const std::vector<double> temperature_before =
            probe_temperature.bytes;
        const std::vector<double> flux_before =
            face_flux_values(as_const(probe_flux));
        const PisoStateCorrectionCertificate certificate_before =
            probe_correction_certificate;
        const Status rejected =
            coupled_scale_probe.correct_coupled_trial_state(
                probe_pressure, probe_correction.view,
                as_const(probe_enthalpy_correction.view), probe_state,
                selected, probe_flux, reductions,
                probe_correction_certificate);
        return rejected.code == expected_code &&
               same_state_correction(probe_correction_certificate,
                                     certificate_before) &&
               atomically_unchanged(
                   velocity_before, pressure_before, enthalpy_before,
                   density_before, temperature_before, flux_before);
      };

  reset_probe_state();
  fill(probe_candidate_density, -0.1);
  passed &= expect(
      static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(probe_candidate,
                                 StatusCode::numerical_failure),
      "non-positive exact-EOS density rejects collectively without a partial state or flux write");

  reset_probe_state();
  fill(probe_candidate_temperature,
       std::numeric_limits<double>::quiet_NaN());
  passed &= expect(
      static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(probe_candidate,
                                 StatusCode::numerical_failure),
      "non-finite exact-EOS temperature rejects collectively without a partial state or flux write");

  reset_probe_state();
  const Status corrected_velocity_rebuilt = rebuild_probe_c1();
  PisoExactThermodynamicCandidateView corrected_velocity_candidate;
  const Status corrected_velocity_gauge =
      corrected_velocity_rebuilt
          ? prepare_probe_gauge(corrected_velocity_candidate)
          : corrected_velocity_rebuilt;
  // The base-layer U is finite and the pressure direction is nonzero, but
  // this staged same-target U has a finite component whose square overflows.
  // An EOS replay that incorrectly used the old U would accept it.
  fill(probe_h_by_a, 1.0e200);
  passed &= expect(
      static_cast<bool>(corrected_velocity_gauge) &&
          reject_probe_candidate(corrected_velocity_candidate,
                                 StatusCode::invalid_plan),
      "legacy exact EOS replays the corrected target-layer velocity and rejects a non-finite Mach without partial publication");

  reset_probe_state();
  fill(probe_pressure_perturbation, -probe_pressure_reference_value);
  fill(probe_correction, 0.0);
  passed &= expect(
          static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(probe_candidate,
                                 StatusCode::invalid_plan),
      "zero exact-EOS absolute pressure invalidates the prepared gauge transaction without a partial write");

  reset_probe_state();
  probe_input.thermophysical_boundary.binding.pressure_reference =
      std::numeric_limits<double>::max();
  fill(probe_pressure_perturbation,
       std::numeric_limits<double>::max());
  fill(probe_correction, 0.0);
  passed &= expect(
          static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(probe_candidate,
                                 StatusCode::invalid_plan),
      "infinite exact-EOS absolute pressure invalidates the prepared gauge transaction without a partial write");

  reset_probe_state();
  PisoExactThermodynamicCandidateView empty_closure_candidate =
      probe_candidate;
  empty_closure_candidate.closure = {};
  passed &= expect(
      static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(empty_closure_candidate,
                                 StatusCode::invalid_plan),
      "empty exact-EOS closure authority rejects before mutation");

  reset_probe_state();
  PisoExactThermodynamicCandidateView stale_revision_candidate =
      probe_candidate;
  ++stale_revision_candidate.temperature.revision;
  passed &= expect(
      static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(stale_revision_candidate,
                                 StatusCode::invalid_plan),
      "stale candidate revision cannot consume a fresh exact-EOS closure");

  reset_probe_state();
  OwnedField foreign_candidate_density = make_field(
      probe_candidate_density.view.field, cells, 1U, 0U,
      probe_candidate_density.view.revision, 3931U);
  fill(foreign_candidate_density, 0.4);
  PisoExactThermodynamicCandidateView foreign_candidate = probe_candidate;
  foreign_candidate.density = as_const(foreign_candidate_density.view);
  passed &= expect(
      static_cast<bool>(rebuild_probe_c1()) &&
          reject_probe_candidate(foreign_candidate,
                                 StatusCode::invalid_plan),
      "same-numbered candidate revision from foreign storage cannot consume exact-EOS closure authority");

  const auto reject_gauge_mutation = [&](auto mutate) {
    reset_probe_state();
    PisoExactThermodynamicCandidateView selected;
    const Status rebuilt = rebuild_probe_c1();
    const Status prepared = rebuilt ? prepare_probe_gauge(selected) : rebuilt;
    if (!prepared) return false;
    mutate(selected);
    return reject_probe_candidate(selected, StatusCode::invalid_plan);
  };
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.closed_gauge.predecessor_pressure_reference;
      }),
      "wrong closed-gauge predecessor rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        selected.closed_gauge.corrector = 2U;
      }),
      "wrong closed-gauge corrector rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.closed_gauge.time;
      }),
      "wrong closed-gauge target time rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.closed_gauge.geometry;
      }),
      "wrong closed-gauge geometry rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.closed_gauge.activity_local_fingerprint;
      }),
      "wrong closed-gauge IBM activity rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.pressure_compressibility.revision;
      }),
      "stale pressure compressibility rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        selected.closed_gauge.shift =
            std::nextafter(selected.closed_gauge.shift,
                           std::numeric_limits<double>::infinity());
      }),
      "mutated closed-gauge shift rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.closed_gauge.output_pressure_reference.pressure_reference;
      }),
      "mutated closed-gauge output-reference token rejects before any state/flux write");
  passed &= expect(
      reject_gauge_mutation([](PisoExactThermodynamicCandidateView& selected) {
        ++selected.closed_gauge.collective_transaction;
      }),
      "mutated closed-gauge collective token rejects before any state/flux write");

  reset_probe_state();
  PisoExactThermodynamicCandidateView alias_candidate;
  const Status alias_rebuilt = rebuild_probe_c1();
  const Status alias_prepared =
      alias_rebuilt ? prepare_probe_gauge(alias_candidate) : alias_rebuilt;
  FaceFluxView internal_cell_aliased_flux = probe_flux;
  internal_cell_aliased_flux.x.base = probe_gradient.view.base;
  internal_cell_aliased_flux.x.storage_identity =
      internal_cell_aliased_flux.y.storage_identity;
  internal_cell_aliased_flux.x.revision_domain =
      internal_cell_aliased_flux.y.revision_domain;
  const std::vector<double> alias_velocity_before = probe_velocity.bytes;
  const std::vector<double> alias_pressure_before =
      probe_pressure_perturbation.bytes;
  const std::vector<double> alias_enthalpy_before = probe_enthalpy.bytes;
  const std::vector<double> alias_density_before = probe_density.bytes;
  const std::vector<double> alias_temperature_before = probe_temperature.bytes;
  const std::vector<double> alias_flux_before =
      face_flux_values(as_const(internal_cell_aliased_flux));
  const std::vector<double> alias_gradient_before = probe_gradient.bytes;
  const PisoStateCorrectionCertificate alias_certificate_before =
      probe_correction_certificate;
  const Status alias_rejected =
      alias_prepared
          ? coupled_scale_probe.correct_coupled_trial_state(
                probe_pressure, probe_correction.view,
                as_const(probe_enthalpy_correction.view), probe_state,
                alias_candidate, internal_cell_aliased_flux, reductions,
                probe_correction_certificate)
          : alias_prepared;
  passed &= expect(
      alias_rejected.code == StatusCode::invalid_plan &&
          same_state_correction(probe_correction_certificate,
                                alias_certificate_before) &&
          probe_velocity.bytes == alias_velocity_before &&
          probe_pressure_perturbation.bytes == alias_pressure_before &&
          probe_enthalpy.bytes == alias_enthalpy_before &&
          probe_density.bytes == alias_density_before &&
          probe_temperature.bytes == alias_temperature_before &&
          face_flux_values(as_const(internal_cell_aliased_flux)) ==
              alias_flux_before &&
          probe_gradient.bytes == alias_gradient_before,
      "output flux aliasing internal pressure-gradient storage rejects before gradient or external writes");

  const std::array<HaloFieldSpec, 1U> krylov_halo_fields{{
      {krylov_vectors.view.field, 1U, 1U}}};
  HaloEngine krylov_halo;
  passed &= expect(static_cast<bool>(krylov_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {krylov_halo_fields.data(),
                        krylov_halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "pressure operator Krylov halo reserves once");
  PressureLinearOperator solve_operator;
  passed &= expect(static_cast<bool>(coupler.bind_pressure_operator(
                       {MPI_COMM_SELF, &krylov_halo, 201U,
                        krylov_vectors.view.field},
                       pressure_system, solve_operator)),
                   "pressure solve operator binds Krylov field once");

  NativeCartesianMgSpec mg_spec;
  mg_spec.communicator = MPI_COMM_SELF;
  mg_spec.geometry = &fixture.geometry;
  mg_spec.patch = fixture.patch;
  mg_spec.boundaries = {MgBoundaryKind::periodic,
                        MgBoundaryKind::periodic,
                        MgBoundaryKind::periodic,
                        MgBoundaryKind::periodic,
                        MgBoundaryKind::periodic,
                        MgBoundaryKind::periodic};
  mg_spec.null_space = MgNullSpace::none;
  const LinearIdentity solve_identity_one{
      2101U, 2102U, 2103U, krylov_workspace.fingerprint(), 2105U};
  const LinearIdentity coupled_solve_identity_one{
      solve_identity_one.symbolic, solve_identity_one.numeric,
      solve_identity_one.hierarchy, coupled_krylov_workspace.fingerprint(),
      2125U};
  mg_spec.identity = solve_identity_one;
  mg_spec.coefficients = {2106U, solve_pressure_one.state, 0.0};
  NativeCartesianMgSpec produced_mg_spec;
  passed &= expect(
      static_cast<bool>(coupler.make_native_pressure_mg_spec(
          MPI_COMM_SELF, solve_identity_one, mg_spec.coefficients,
          produced_mg_spec)) &&
          produced_mg_spec.policy.pre_sweeps == 1U &&
          produced_mg_spec.policy.post_sweeps == 2U &&
          produced_mg_spec.policy.point_smoother ==
              MgPointSmootherKind::chebyshev_jacobi &&
          produced_mg_spec.policy.cycle == MgCycleKind::f_cycle &&
          produced_mg_spec.policy.chebyshev_lower_spectrum_fraction == 0.3,
      "pressure producer publishes certified Chebyshev 1/2 F-cycle policy");
  NativeCartesianMgSpec coupled_mg_spec = produced_mg_spec;
  coupled_mg_spec.identity = coupled_solve_identity_one;
  NativeCartesianMgSpec predecessor_mg_spec = produced_mg_spec;
  predecessor_mg_spec.policy.pre_sweeps = 3U;
  predecessor_mg_spec.policy.post_sweeps = 3U;
  MgWorkspaceRequirements mg_requirements;
  passed &= expect(static_cast<bool>(make_mg_workspace_requirements(
                       MPI_COMM_SELF, fixture.geometry, fixture.patch,
                       produced_mg_spec.policy, 2107U, mg_requirements)),
                   "pressure Native MG workspace requirements compile");
  MgWorkspaceRequirements predecessor_mg_requirements;
  passed &= expect(
      static_cast<bool>(make_mg_workspace_requirements(
          MPI_COMM_SELF, fixture.geometry, fixture.patch,
          predecessor_mg_spec.policy, 2107U, predecessor_mg_requirements)) &&
          predecessor_mg_requirements.fingerprint ==
              mg_requirements.fingerprint &&
          predecessor_mg_requirements.collective_fingerprint ==
              mg_requirements.collective_fingerprint &&
          predecessor_mg_requirements.total_doubles ==
              mg_requirements.total_doubles &&
          predecessor_mg_requirements.arena_shape.x ==
              mg_requirements.arena_shape.x &&
          predecessor_mg_requirements.arena_shape.y ==
              mg_requirements.arena_shape.y &&
          predecessor_mg_requirements.arena_shape.z ==
              mg_requirements.arena_shape.z,
      "pressure F/1/2 and F/3/3 retain the same exact layout workspace capacity");
  constexpr StorageIdentity mg_storage = 2902U;
  OwnedField mg_vectors = make_field(80U, mg_requirements.arena_shape, 1U,
                                     1U, 2108U, mg_storage);
  MgWorkspace mg_workspace;
  passed &= expect(static_cast<bool>(MgWorkspace::bind(
                       mg_requirements, mg_vectors.view, mg_workspace)),
                   "pressure Native MG workspace binds once");
  const std::array<HaloFieldSpec, 1U> mg_halo_fields{{
      {mg_vectors.view.field, 1U, 1U}}};
  HaloEngine mg_halo;
  passed &= expect(static_cast<bool>(mg_halo.reserve(
                       MPI_COMM_SELF, mg_requirements.levels[0U].patch,
                       {mg_halo_fields.data(), mg_halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "pressure finest MG halo reserves once");
  std::vector<HaloEngine> coarse_halos(mg_requirements.level_count - 1U);
  std::vector<HaloEngine*> coarse_halo_pointers(coarse_halos.size());
  for (std::size_t level = 1U; level < mg_requirements.level_count; ++level) {
    passed &= expect(static_cast<bool>(coarse_halos[level - 1U].reserve(
                         MPI_COMM_SELF,
                         mg_requirements.levels[level].patch,
                         {mg_halo_fields.data(), mg_halo_fields.size()},
                         fixture.boundary.halo_topology())),
                     "pressure coarse MG halo reserves once");
    coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
  }
  const MgRuntimeServices mg_services{
      &mg_halo, &reductions, &mg_workspace,
      {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
  NativeCartesianMgPlan pressure_mg;
  MgPlanCounters mg_counters;
  passed &= expect(static_cast<bool>(coupler.compile_native_pressure_mg(
                       solve_pressure_one, produced_mg_spec, mg_services,
                       pressure_system, pressure_mg, &mg_counters)) &&
                       pressure_mg.certificate().preconditioner_class ==
                           LinearPreconditionerClass::flexible &&
                       pressure_mg.certificate().status_scope ==
                           LinearPreconditionerStatusScope::collective &&
                       pressure_mg.certificate().apply_lifecycle ==
                           LinearPreconditionerApplyLifecycle::prepared_batch,
                   "pressure Native MG declares collective status after internal consensus");
  NativeCartesianMgPlan coupled_pressure_mg;
  MgPlanCounters coupled_mg_counters;
  passed &= expect(static_cast<bool>(coupler.compile_native_pressure_mg(
                       solve_pressure_one, coupled_mg_spec, mg_services,
                       pressure_system, coupled_pressure_mg,
                       &coupled_mg_counters)),
                   "coupled pressure-energy MG lifecycle compiles independently");
  NativeCartesianMgPlan predecessor_pressure_mg;
  MgPlanCounters predecessor_mg_counters;
  const Status predecessor_compiled = coupler.compile_native_pressure_mg(
      solve_pressure_one, predecessor_mg_spec, mg_services, pressure_system,
      predecessor_pressure_mg, &predecessor_mg_counters);
  passed &= expect(
      static_cast<bool>(predecessor_compiled) &&
          pressure_mg.symbolic_fingerprint() !=
              predecessor_pressure_mg.symbolic_fingerprint() &&
          pressure_mg.certificate().collective_fingerprint !=
              predecessor_pressure_mg.certificate().collective_fingerprint &&
          pressure_mg.workspace_storage_address() ==
              predecessor_pressure_mg.workspace_storage_address() &&
          predecessor_mg_counters.symbolic_builds == 1U &&
          predecessor_mg_counters.hierarchy_rebuilds == 1U,
      "pressure F/1/2 changes structural and prepared identity from F/3/3 without changing bound workspace");
  const std::uintptr_t krylov_address =
      krylov_workspace.vector_storage_address();
  const std::uintptr_t mg_address = pressure_mg.workspace_storage_address();

  OwnedField pressure_correction =
      make_field(correction_field, cells, 1U, 1U, 2109U, 2903U);
  OwnedField coupled_pressure_correction =
      make_field(correction_field, cells, 1U, 1U, 2127U, 2931U);
  fill(pressure_correction, 0.0);
  fill(coupled_pressure_correction, 0.0);
  {
    PisoPressureSolveEpoch abandoned_epoch;
    Status abandoned = abandoned_epoch.begin(piso);
    if (abandoned) {
      abandoned = abandoned_epoch.solve(
          piso, 1U, solve_pressure_one, solve_identity_one,
          mg_spec.coefficients, coupler, solve_operator, pressure_mg,
          pressure_system, pressure_correction.view, krylov_workspace,
          reductions, nullptr, &mg_counters);
    }
    passed &= expect(
        static_cast<bool>(abandoned) &&
            abandoned_epoch.solve_calls() == 1U &&
            krylov_workspace.recycle_correction_count_for_test() > 0U,
        "a scoped C1 epoch owns a live capture session before abandonment");
  }
  fill(pressure_correction, 0.0);
  PisoPressureSolveEpoch solve_epoch;
  passed &= expect(static_cast<bool>(solve_epoch.begin(piso)),
                   "pressure solve epoch begins once");
  passed &= expect(static_cast<bool>(solve_epoch.prepare_linear_lifecycle(
                       piso, 1U, solve_pressure_one, solve_identity_one,
                       mg_spec.coefficients, coupler, solve_operator,
                       pressure_mg, pressure_system,
                       pressure_correction.view, krylov_workspace,
                       &mg_counters)) &&
                       static_cast<bool>(solve_epoch.solve_prepared(
                           solve_operator,
                           PisoPressureSolveContract::pressure_continuity,
                           reductions)) &&
                       solve_epoch.solve_calls() == 1U,
                   "prepared corrector one performs exactly one FGMRES solve");
  passed &= expect(krylov_workspace.recycle_correction_count_for_test() > 0U,
                   "corrector one leaves captured corrections before mutation");

  const auto seed_capture = [&]() {
    fill(pressure_correction, 0.0);
    Status seeded = krylov_workspace.recycle_begin_capture_for_test(
        cells, solve_identity_one.fingerprint);
    if (seeded)
      seeded = krylov_workspace.recycle_capture_cycle_start_for_test(
          as_const(pressure_correction.view), reductions);
    if (seeded) {
      fill(pressure_correction, 1.0);
      seeded = krylov_workspace.recycle_capture_cycle_publish_for_test(
          as_const(pressure_correction.view), reductions);
    }
    return static_cast<bool>(seeded) &&
           krylov_workspace.recycle_correction_count_for_test() == 1U;
  };

  PisoPressureSolveEpoch failed_epoch;
  passed &= expect(static_cast<bool>(failed_epoch.begin(piso)),
                   "a failed pressure retry probe begins independently");
  const Status invalid_retry = failed_epoch.solve(
      piso, 2U, solve_pressure_one, solve_identity_one,
      mg_spec.coefficients, coupler, solve_operator, pressure_mg,
      pressure_system, pressure_correction.view, krylov_workspace, reductions,
      nullptr, &mg_counters);
  const Status stale_projection = krylov_workspace
                                      .recycle_begin_projection_for_test(
                                          cells, solve_identity_one.fingerprint);
  passed &= expect(
      invalid_retry.code == StatusCode::invalid_plan &&
          krylov_workspace.recycle_correction_count_for_test() == 0U &&
          stale_projection.code == StatusCode::invalid_plan,
      "invalid pressure retry clears captured C1 session before any projection");

  passed &= expect(seed_capture(),
                   "view-failure probe seeds one captured correction");
  PisoPressureSolveEpoch view_failure_epoch;
  passed &= expect(static_cast<bool>(view_failure_epoch.begin(piso)),
                   "view-failure epoch begins independently");
  FieldView invalid_correction = pressure_correction.view;
  invalid_correction.interior = {cells.x + 1, cells.y, cells.z};
  const Status view_failure = view_failure_epoch.solve(
      piso, 1U, solve_pressure_one, solve_identity_one,
      mg_spec.coefficients, coupler, solve_operator, pressure_mg,
      pressure_system, invalid_correction, krylov_workspace, reductions,
      nullptr, &mg_counters);
  const Status view_projection = krylov_workspace.recycle_begin_projection_for_test(
      cells, solve_identity_one.fingerprint);
  passed &= expect(
      view_failure.code == StatusCode::invalid_plan &&
          krylov_workspace.recycle_correction_count_for_test() == 0U &&
          view_projection.code == StatusCode::invalid_plan,
      "invalid pressure correction view clears captured C1 session");

  passed &= expect(seed_capture(),
                   "plan-failure probe seeds one captured correction");
  PisoPressureSolveEpoch plan_failure_epoch;
  passed &= expect(static_cast<bool>(plan_failure_epoch.begin(piso)),
                   "plan-failure epoch begins independently");
  const Status plan_failure = plan_failure_epoch.solve(
      PisoPlan{}, 1U, solve_pressure_one, solve_identity_one,
      mg_spec.coefficients, coupler, solve_operator, pressure_mg,
      pressure_system, pressure_correction.view, krylov_workspace, reductions,
      nullptr, &mg_counters);
  const Status plan_projection = krylov_workspace.recycle_begin_projection_for_test(
      cells, solve_identity_one.fingerprint);
  passed &= expect(
      plan_failure.code == StatusCode::invalid_plan &&
          krylov_workspace.recycle_correction_count_for_test() == 0U &&
          plan_projection.code == StatusCode::invalid_plan,
      "invalid pressure plan clears captured C1 session");

  passed &= expect(seed_capture(),
                   "refresh-failure probe seeds one captured correction");
  PisoPressureSolveEpoch refresh_failure_epoch;
  passed &= expect(static_cast<bool>(refresh_failure_epoch.begin(piso)),
                   "refresh-failure epoch begins independently");
  const MgCoefficientIdentity invalid_coefficients{2117U,
                                                   solve_pressure_one.state + 1U,
                                                   0.0};
  const Status refresh_failure = refresh_failure_epoch.solve(
      piso, 1U, solve_pressure_one, solve_identity_one,
      invalid_coefficients, coupler, solve_operator, pressure_mg,
      pressure_system, pressure_correction.view, krylov_workspace, reductions,
      nullptr, &mg_counters);
  const Status refresh_projection = krylov_workspace
                                        .recycle_begin_projection_for_test(
                                            cells,
                                            solve_identity_one.fingerprint);
  passed &= expect(
      refresh_failure.code == StatusCode::invalid_plan &&
          krylov_workspace.recycle_correction_count_for_test() == 0U &&
          refresh_projection.code == StatusCode::invalid_plan,
      "pressure lifecycle refresh failure clears captured C1 session");

  passed &= expect(seed_capture(),
                   "exact-identity failure probe seeds one captured correction");
  PressureLinearOperator invalid_exact_operator;
  PisoPressureSolveEpoch exact_failure_epoch;
  passed &= expect(static_cast<bool>(exact_failure_epoch.begin(piso)),
                   "exact-identity failure epoch begins independently");
  const Status exact_failure = exact_failure_epoch.solve(
      piso, 1U, solve_pressure_one, solve_identity_one,
      mg_spec.coefficients, coupler, solve_operator, invalid_exact_operator,
      pressure_mg, pressure_system, pressure_correction.view, krylov_workspace,
      reductions, nullptr, &mg_counters);
  const Status exact_projection = krylov_workspace
                                     .recycle_begin_projection_for_test(
                                         cells, solve_identity_one.fingerprint);
  passed &= expect(
      exact_failure.code == StatusCode::invalid_plan &&
          krylov_workspace.recycle_correction_count_for_test() == 0U &&
          exact_projection.code == StatusCode::invalid_plan,
      "exact-operator identity failure clears captured C1 session");

  fill(pressure_correction, 0.0);
  PisoPressureSolveEpoch retry_epoch;
  passed &= expect(static_cast<bool>(retry_epoch.begin(piso)) &&
                       static_cast<bool>(retry_epoch.solve(
                           piso, 1U, solve_pressure_one, solve_identity_one,
                           mg_spec.coefficients, coupler, solve_operator,
                           pressure_mg, pressure_system,
                           pressure_correction.view, krylov_workspace,
                           reductions, nullptr, &mg_counters)) &&
                       retry_epoch.solve_calls() == 1U &&
                       krylov_workspace.recycle_correction_count_for_test() <=
                           kLinearRecycleMaximumDirections,
                   "C1 can re-begin and complete after the failed retry");
  PisoPressureSolveEpoch coupled_epoch;
  const Status coupled_begin = coupled_epoch.begin(piso);
  const Status coupled_prepare_one =
      coupled_begin
          ? coupled_epoch.prepare_linear_lifecycle(
                piso, 1U, solve_pressure_one, coupled_solve_identity_one,
                mg_spec.coefficients, coupler, solve_operator,
                coupled_pressure_mg, pressure_system,
                coupled_pressure_correction.view, coupled_krylov_workspace,
                &coupled_mg_counters)
          : coupled_begin;
  const Status coupled_solve_one =
      coupled_prepare_one
          ? coupled_epoch.solve_prepared(
                solve_operator,
                PisoPressureSolveContract::continuity_energy_coupled,
                reductions)
          : coupled_prepare_one;
  passed &= expect(static_cast<bool>(coupled_solve_one) &&
                       coupled_epoch.solve_calls() == 1U,
                   "coupled pressure-energy C1 solves after caller-visible preparation");

  OwnedField trial_pressure =
      make_field(2U, cells, 1U, 1U, 2110U, 2904U);
  OwnedField trial_enthalpy =
      make_field(3U, cells, 1U, 1U, 2112U, 2932U);
  OwnedField trial_temperature =
      make_field(4U, cells, 1U, 1U, 2113U, 2933U);
  OwnedField c1_enthalpy_correction =
      make_field(91U, cells, 1U, 0U, 2114U, 2934U);
  OwnedField c1_candidate_enthalpy =
      make_field(92U, cells, 1U, 0U, 2115U, 2935U);
  OwnedField c1_candidate_density =
      make_field(93U, cells, 1U, 0U, 2116U, 2936U);
  OwnedField c1_candidate_temperature =
      make_field(94U, cells, 1U, 0U, 2117U, 2937U);
  OwnedField c1_candidate_compressibility =
      make_field(99U, cells, 1U, 0U, 2118U, 2938U);
  double trial_enthalpy_value = 0.0;
  double c1_candidate_enthalpy_value = 0.0;
  double c2_candidate_enthalpy_value = 0.0;
  double exact_cp = 0.0;
  double exact_gas_constant = 0.0;
  passed &= expect(
      static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
          300.0, {}, trial_enthalpy_value, exact_cp,
          exact_gas_constant)) &&
          static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
              301.0, {}, c1_candidate_enthalpy_value, exact_cp,
              exact_gas_constant)) &&
          static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
              302.0, {}, c2_candidate_enthalpy_value, exact_cp,
              exact_gas_constant)),
      "coupled C1/C2 fixture constructs physical enthalpy layers");
  fill(trial_pressure, 0.0);
  fill(trial_enthalpy, trial_enthalpy_value);
  fill(trial_temperature, 300.0);
  fill(c1_enthalpy_correction,
       c1_candidate_enthalpy_value - trial_enthalpy_value);
  fill(c1_candidate_enthalpy, c1_candidate_enthalpy_value);
  fill(c1_candidate_temperature, 0.0);
  fill(c1_candidate_compressibility, 0.0);
  FaceFluxStorage trial_flux_storage;
  FaceFluxView trial_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, trial_flux_storage)) &&
          static_cast<bool>(trial_flux_storage.workspace_view(
              0U, 2111U, trial_flux)),
      "corrector-one trial flux workspace allocates cold");
  FieldView corrected_velocity = velocity.view;
  FieldView corrected_density = density.view;
  FieldView corrected_enthalpy = trial_enthalpy.view;
  FieldView corrected_temperature = trial_temperature.view;
  ++corrected_velocity.revision;
  ++corrected_density.revision;
  ++corrected_enthalpy.revision;
  ++corrected_temperature.revision;
  const PisoCoupledStateView corrected_state{
      corrected_velocity, trial_pressure.view, corrected_enthalpy,
      corrected_density, corrected_temperature};
  Status c1_profile_status;
  for (std::int32_t zc = 0; zc < cells.z && c1_profile_status; ++zc)
    for (std::int32_t yc = 0; yc < cells.y && c1_profile_status; ++yc)
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const double target_density =
            corrected_density.unchecked(cell, 0U) +
            drho_dp.view.unchecked(cell, 0U) *
                pressure_correction.view.unchecked(cell, 0U);
        const double target_absolute_pressure =
            101325.0 + trial_pressure.view.unchecked(cell, 0U) +
            pressure_correction.view.unchecked(cell, 0U);
        if (!(target_density > 0.0) ||
            !(target_absolute_pressure > 0.0)) {
          c1_profile_status = {StatusCode::numerical_failure, 1U};
          break;
        }
        const double target_temperature =
            target_absolute_pressure /
            (exact_gas_constant * target_density);
        double candidate_h = 0.0;
        double candidate_cp = 0.0;
        double candidate_gas_constant = 0.0;
        c1_profile_status = fixture.thermodynamics.mixture_enthalpy(
            target_temperature, {}, candidate_h, candidate_cp,
            candidate_gas_constant);
        if (!c1_profile_status) break;
        c1_candidate_enthalpy.view.unchecked(cell, 0U) = candidate_h;
        c1_enthalpy_correction.view.unchecked(cell, 0U) =
            candidate_h - trial_enthalpy.view.unchecked(cell, 0U);
      }
  passed &= expect(static_cast<bool>(c1_profile_status),
                   "C1 exact profile preserves its assembled density target");
  PisoExactEosClosureIdentity c1_closure;
  c1_closure.thermodynamics =
      pressure_input.pressure_reference.thermodynamics;
  c1_closure.pressure_reference =
      pressure_input.pressure_reference.pressure_reference;
  c1_closure.composition = exact_composition_identity_for_test(
      c1_closure.thermodynamics, {}, cells);
  c1_closure.pressure_state = make_piso_field_revision_identity(
      as_const(trial_pressure.view));
  c1_closure.pressure_correction = make_piso_field_revision_identity(
      as_const(pressure_correction.view));
  c1_closure.enthalpy_state = make_piso_field_revision_identity(
      as_const(corrected_enthalpy));
  c1_closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(c1_enthalpy_correction.view));
  c1_closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(c1_candidate_enthalpy.view));
  c1_closure.candidate_density = make_piso_field_revision_identity(
      as_const(c1_candidate_density.view));
  c1_closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(c1_candidate_temperature.view));
  c1_closure.closure = 2939U;
  PisoExactThermodynamicCandidateView c1_candidate;
  c1_candidate.enthalpy = as_const(c1_candidate_enthalpy.view);
  c1_candidate.density = as_const(c1_candidate_density.view);
  c1_candidate.temperature = as_const(c1_candidate_temperature.view);
  c1_candidate.closure = c1_closure;
  c1_candidate.pressure_compressibility =
      as_const(c1_candidate_compressibility.view);
  c1_candidate.pressure_compressibility.field = drho_dp.view.field;
  passed &= expect(static_cast<bool>(fill_exact_thermodynamic_candidate(
                       fixture.thermodynamics, 101325.0,
                       as_const(trial_pressure.view),
                       as_const(pressure_correction.view),
                       as_const(c1_candidate_enthalpy.view),
                       as_const(corrected_velocity), {},
                       c1_candidate_density.view,
                       c1_candidate_temperature.view,
                       c1_candidate_compressibility.view)),
                   "C1 exact candidate is replayed from p/h/Y/U");
  passed &= expect(
      static_cast<bool>(prepare_closed_gauge(
          fixture.equations.pressure_reference(), solve_pressure_one,
          pressure_input.pressure_reference, 101325.0,
          as_const(trial_pressure.view),
          as_const(pressure_correction.view),
          c1_candidate.pressure_compressibility,
          c1_candidate.closure.closure, reductions,
          c1_candidate.closed_gauge)),
      "C1 exact state prepares one closed-gauge transition");
  PisoStateCorrectionCertificate corrected_one;
  passed &= expect(static_cast<bool>(coupler.correct_coupled_trial_state(
                       solve_pressure_one, pressure_correction.view,
                       as_const(c1_enthalpy_correction.view), corrected_state,
                       c1_candidate, trial_flux, reductions,
                       corrected_one)) &&
                       corrected_one.valid() &&
                       corrected_one.corrector == 1U &&
                       corrected_one.closure == PisoStateClosure::exact_eos &&
                       corrected_one.enthalpy == corrected_enthalpy.revision &&
                       corrected_one.temperature ==
                           corrected_temperature.revision,
                   "coupled corrector one atomically publishes U/pi/h/rho/T and trial flux under exact-EOS authority");
  const auto cell_offset = [cells](Int3 cell) {
    return static_cast<std::size_t>(cell.x) +
           static_cast<std::size_t>(cells.x) *
               (static_cast<std::size_t>(cell.y) +
                static_cast<std::size_t>(cells.y) *
                    static_cast<std::size_t>(cell.z));
  };
  const auto face_offset = [](Int3 extents, Int3 face) {
    return static_cast<std::size_t>(face.x) +
           static_cast<std::size_t>(extents.x) *
               (static_cast<std::size_t>(face.y) +
                static_cast<std::size_t>(extents.y) *
                    static_cast<std::size_t>(face.z));
  };
  const std::size_t cell_count = static_cast<std::size_t>(cells.x) *
                                 static_cast<std::size_t>(cells.y) *
                                 static_cast<std::size_t>(cells.z);
  std::vector<double> c1_pressure(cell_count);
  std::vector<double> c1_velocity(3U * cell_count);
  double maximum_c1_velocity_correction = 0.0;
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const std::size_t offset = cell_offset(cell);
        c1_pressure[offset] = trial_pressure.view.unchecked(cell, 0U);
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double value = corrected_velocity.unchecked(cell, component);
          c1_velocity[3U * offset + component] = value;
          maximum_c1_velocity_correction = std::max(
              maximum_c1_velocity_correction,
              std::abs(value - h_by_a.view.unchecked(cell, component)));
        }
      }
    }
  }
  const ConstFaceFluxView c1_flux_view = as_const(trial_flux);
  const std::array<ConstFaceFieldView, 3U> c1_flux_faces{
      c1_flux_view.x, c1_flux_view.y, c1_flux_view.z};
  std::array<std::vector<double>, 3U> c1_flux;
  double maximum_c1_flux_correction = 0.0;
  const ConstFaceFluxView c1_base_flux = as_const(phi_h_by_a);
  const std::array<ConstFaceFieldView, 3U> c1_base_faces{
      c1_base_flux.x, c1_base_flux.y, c1_base_flux.z};
  for (std::size_t axis = 0U; axis < c1_flux_faces.size(); ++axis) {
    const Int3 extents = c1_flux_faces[axis].extents;
    c1_flux[axis].resize(static_cast<std::size_t>(extents.x) *
                         static_cast<std::size_t>(extents.y) *
                         static_cast<std::size_t>(extents.z));
    for (std::int32_t zf = 0; zf < extents.z; ++zf) {
      for (std::int32_t yf = 0; yf < extents.y; ++yf) {
        for (std::int32_t xf = 0; xf < extents.x; ++xf) {
          const Int3 face{xf, yf, zf};
          const double value = c1_flux_faces[axis].unchecked(face);
          c1_flux[axis][face_offset(extents, face)] = value;
          maximum_c1_flux_correction = std::max(
              maximum_c1_flux_correction,
              std::abs(value - c1_base_faces[axis].unchecked(face)));
        }
      }
    }
  }
  passed &= expect(
      std::max(maximum_c1_velocity_correction,
               maximum_c1_flux_correction) > 1.0e-12,
      "two-corrector fixture exercises a nonzero C1 state correction");
  double maximum_continuity = 0.0;
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const double density_defect =
            (input.bdf.a0 * corrected_density.unchecked(cell, 0U) +
             input.bdf.a1 * density_accepted.view.unchecked(cell, 0U) +
             input.bdf.a2 * density_previous.view.unchecked(cell, 0U)) *
            cell_volume(fixture, cell);
        const double flux_divergence =
            trial_flux.x.unchecked({xc + 1, yc, zc}) -
            trial_flux.x.unchecked(cell) +
            trial_flux.y.unchecked({xc, yc + 1, zc}) -
            trial_flux.y.unchecked(cell) +
            trial_flux.z.unchecked({xc, yc, zc + 1}) -
            trial_flux.z.unchecked(cell);
        maximum_continuity = std::max(
            maximum_continuity,
            std::abs(density_defect + flux_divergence));
        passed &= expect(
            std::abs(trial_pressure.view.unchecked(cell, 0U) -
                     (pressure_correction.view.unchecked(cell, 0U) -
                      c1_candidate.closed_gauge.shift)) <
                    1.0e-12 &&
                std::abs(corrected_density.unchecked(cell, 0U) -
                         (1.0 + 0.02 *
                                    pressure_correction.view.unchecked(
                                        cell, 0U))) <
                    1.0e-12 &&
                std::abs(corrected_velocity.unchecked(cell, 0U) -
                         (h_by_a.view.unchecked(cell, 0U) -
                          r_au.view.unchecked(cell, 0U) *
                              pressure_gradient.view.unchecked(cell, 0U))) <
                    1.0e-12,
            "U/pi/rho use the same pressure-correction gradient path");
      }
    }
  }
  passed &= expect(maximum_continuity < 2.0e-8,
                   "corrector-one flux satisfies the assembled continuity equation");

  input.corrector = 2U;
  input.temporal_reference = {};
  input.committed_face_history = {};
  input.prior_corrector = corrected_one.state;
  input.pressure_reference = corrected_one.output_pressure_reference;
  input.thermophysical_boundary.binding.pressure_reference =
      c1_candidate.closed_gauge.next_pressure_reference;
  input.trial_velocity = as_const(corrected_velocity);
  input.density = corrected_density;
  ConstFaceFluxView corrected_trial_flux = as_const(trial_flux);
  corrected_trial_flux.revision = corrected_one.state;
  input.trial_flux = corrected_trial_flux;
  ++input.momentum.state;
  PisoIntermediateCertificate solve_intermediate_two;
  passed &= expect(static_cast<bool>(coupler.refresh(
                       input, solve_intermediate_two)) &&
                       solve_intermediate_two.corrector == 2U,
                   "corrector two rebuilds intermediates before solve");
  double maximum_c2_cell_base_error = 0.0;
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const std::size_t offset = cell_offset(cell);
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          maximum_c2_cell_base_error = std::max(
              maximum_c2_cell_base_error,
              std::abs(h_by_a.view.unchecked(cell, component) -
                       c1_velocity[3U * offset + component]));
        }
      }
    }
  }
  passed &= expect(
      maximum_c2_cell_base_error < 1.0e-12,
      "corrector-two cell base is the C1 corrected velocity");
  double maximum_c2_face_base_error = 0.0;
  const ConstFaceFluxView c2_base_flux = as_const(phi_h_by_a);
  const std::array<ConstFaceFieldView, 3U> c2_base_faces{
      c2_base_flux.x, c2_base_flux.y, c2_base_flux.z};
  for (std::size_t axis = 0U; axis < c2_base_faces.size(); ++axis) {
    const Int3 extents = c2_base_faces[axis].extents;
    for (std::int32_t zf = 0; zf < extents.z; ++zf) {
      for (std::int32_t yf = 0; yf < extents.y; ++yf) {
        for (std::int32_t xf = 0; xf < extents.x; ++xf) {
          const Int3 face{xf, yf, zf};
          maximum_c2_face_base_error = std::max(
              maximum_c2_face_base_error,
              std::abs(c2_base_faces[axis].unchecked(face) -
                       c1_flux[axis][face_offset(extents, face)]));
        }
      }
    }
  }
  passed &= expect(
      maximum_c2_face_base_error < 1.0e-12,
      "corrector-two face base is the C1 corrected trial flux");
  pressure_input.intermediate = solve_intermediate_two;
  pressure_input.pressure_reference = input.pressure_reference;
  pressure_input.density_trial = as_const(corrected_density);
  PressureCorrectionCertificate solve_pressure_two;
  passed &= expect(static_cast<bool>(coupler.assemble_pressure_system(
                       pressure_input, pressure_system,
                       solve_pressure_two)) &&
                       solve_pressure_two.corrector == 2U,
                   "corrector-two pressure system refreshes exact numerics");
  LinearIdentity solve_identity_two = solve_identity_one;
  solve_identity_two.numeric = 2112U;
  solve_identity_two.hierarchy = 2113U;
  solve_identity_two.fingerprint = 2115U;
  const MgCoefficientIdentity coefficients_two{
      2116U, solve_pressure_two.state, 0.0};
  fill(pressure_correction, 0.0);
  const Status corrector_two_solve = retry_epoch.solve(
      piso, 2U, solve_pressure_two, solve_identity_two, coefficients_two,
      coupler, solve_operator, pressure_mg, pressure_system,
      pressure_correction.view, krylov_workspace, reductions, nullptr,
      &mg_counters);
  passed &= expect(corrector_two_solve.code == StatusCode::rejected_step &&
                       retry_epoch.solve_calls() == 2U,
                   "pressure-only corrector two rejects a zero linear residual whose real-EOS continuity audit is nonzero");
  PisoAttemptReport observed_pressure;
  passed &= expect(
      static_cast<bool>(retry_epoch.observe(observed_pressure)) &&
          observed_pressure.pressure_solve_calls == 2U &&
          observed_pressure.pressure[0U].convergence_audits == 0U &&
          observed_pressure.pressure[1U].status.code ==
              StatusCode::rejected_step &&
          observed_pressure.pressure[1U].termination ==
              LinearTermination::convergence_audit_failure &&
          observed_pressure.pressure[1U].iterations == 0U &&
          observed_pressure.pressure[1U].initial_true_residual == 0.0 &&
          observed_pressure.pressure[1U].convergence_audits == 1U &&
          observed_pressure.pressure[1U].convergence_rejections == 1U &&
          observed_pressure.pressure[0U].recycle_offered_directions == 0U &&
          !observed_pressure.pressure[0U].recycle_projection_attempted &&
          !observed_pressure.pressure[0U].recycle_projection_accepted &&
          observed_pressure.pressure[0U].recycle_capture_cycle_attempts >=
              observed_pressure.pressure[0U].recycle_cycle_corrections &&
          observed_pressure.pressure[0U].recycle_capture_vector_passes ==
              2U * observed_pressure.pressure[0U]
                       .recycle_capture_cycle_attempts &&
          observed_pressure.pressure[0U].recycle_capture_reduction_calls ==
              observed_pressure.pressure[0U].recycle_capture_cycle_attempts &&
          observed_pressure.pressure[0U]
                  .recycle_capture_blocking_operations ==
              2U * observed_pressure.pressure[0U]
                       .recycle_capture_cycle_attempts &&
          observed_pressure.pressure[1U].recycle_cycle_corrections == 0U &&
          observed_pressure.pressure[1U].recycle_capture_vector_passes == 0U &&
          observed_pressure.pressure[1U].recycle_capture_cycle_attempts == 0U &&
          observed_pressure.pressure[1U].recycle_capture_reduction_calls ==
              0U &&
          observed_pressure.pressure[1U]
                  .recycle_capture_blocking_operations == 0U &&
          observed_pressure.pressure[1U].recycle_offered_directions == 0U &&
          !observed_pressure.pressure[1U].recycle_projection_attempted &&
          !observed_pressure.pressure[1U].recycle_projection_accepted &&
          observed_pressure.pressure[1U].final_convergence_metric >
              observed_pressure.pressure[1U].convergence_limit,
      "the real-EOS audit exposes and rejects the pressure-only C2 split before publication");
  LinearIdentity coupled_solve_identity_two = solve_identity_two;
  coupled_solve_identity_two.workspace = coupled_krylov_workspace.fingerprint();
  coupled_solve_identity_two.fingerprint = 2126U;
  fill(coupled_pressure_correction, 0.0);
  passed &= expect(
      static_cast<bool>(coupled_epoch.prepare_linear_lifecycle(
          piso, 2U, solve_pressure_two, coupled_solve_identity_two,
          coefficients_two, coupler, solve_operator, coupled_pressure_mg,
          pressure_system, coupled_pressure_correction.view,
          coupled_krylov_workspace, &coupled_mg_counters)) &&
          static_cast<bool>(coupled_epoch.solve_prepared(
              solve_operator,
              PisoPressureSolveContract::continuity_energy_coupled,
              reductions)) &&
          coupled_epoch.solve_calls() == 2U,
      "coupled pressure-energy C2 solves without pressure-only sealing");
  PisoAttemptReport observed_coupled_pressure;
  passed &= expect(
      static_cast<bool>(coupled_epoch.observe(observed_coupled_pressure)) &&
          observed_coupled_pressure.pressure_solve_calls == 2U &&
          observed_coupled_pressure.pressure[0U].convergence_audits == 0U &&
          observed_coupled_pressure.pressure[1U].convergence_audits == 0U &&
          observed_coupled_pressure.pressure[1U].convergence_rejections == 0U,
      "coupled solve contract bypasses the legacy C2 continuity audit");

  FieldRegistry transaction_registry;
  FieldSchema transaction_schema;
  FieldId transaction_dependency = 0U;
  passed &= expect(
      static_cast<bool>(transaction_registry.declare_field(
          "piso.final_state", 1U, 0U, transaction_dependency)) &&
          static_cast<bool>(transaction_registry.freeze(transaction_schema)),
      "PISO final-state transaction schema freezes");
  const std::array transaction_requests{ArenaFieldRequest{
      transaction_dependency, {1, 1, 1}, {0U},
      FieldLifetime::state_layer}};
  ArenaLayout transaction_layout;
  StateLayers transaction_layers;
  AttemptTransaction transaction;
  passed &= expect(
      static_cast<bool>(ArenaLayout::compile(
          transaction_schema,
          {transaction_requests.data(), transaction_requests.size()},
          transaction_layout)) &&
          static_cast<bool>(StateLayers::allocate(transaction_layout,
                                                  transaction_layers)) &&
          static_cast<bool>(AttemptTransaction::create(
              transaction_layers.field_count(), 1U,
              transaction_layers.field_count(), transaction)),
      "PISO final-state transaction resources allocate cold");
  FaceFluxStorage final_flux_storage;
  FinalFaceFluxAuthority final_flux_authority;
  FinalFaceFluxWriter final_flux_writer;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_final(
          cells, final_flux_storage)) &&
          static_cast<bool>(final_flux_authority.claim(
              piso.pressure_stage(), piso.final_flux_slot(), transaction,
              final_flux_writer)) &&
          static_cast<bool>(transaction.begin(transaction_layers)) &&
          static_cast<bool>(transaction.revise_trial(
              transaction_dependency)),
      "only the second pressure stage acquires final-flux authority");
  const RevisionDependency final_dependency{
      AttemptTransaction::field_revision_source(transaction_dependency),
      transaction.trial_revision(transaction_dependency)};
  PendingFaceFluxView pending_final_flux;
  passed &= expect(static_cast<bool>(final_flux_writer.begin_pending(
                       transaction, final_flux_storage,
                       pending_final_flux)),
                   "corrector two acquires one pending final-flux replica");
  FieldView final_velocity = corrected_velocity;
  FieldView final_density = corrected_density;
  FieldView final_pressure = trial_pressure.view;
  FieldView final_enthalpy = corrected_enthalpy;
  FieldView final_temperature = corrected_temperature;
  ++final_velocity.revision;
  ++final_density.revision;
  ++final_pressure.revision;
  ++final_enthalpy.revision;
  ++final_temperature.revision;
  OwnedField c2_enthalpy_correction =
      make_field(95U, cells, 1U, 0U, 2211U, 2940U);
  OwnedField c2_candidate_enthalpy =
      make_field(96U, cells, 1U, 0U, 2212U, 2941U);
  OwnedField c2_candidate_density =
      make_field(97U, cells, 1U, 0U, 2213U, 2942U);
  OwnedField c2_candidate_temperature =
      make_field(98U, cells, 1U, 0U, 2214U, 2943U);
  OwnedField c2_candidate_compressibility =
      make_field(100U, cells, 1U, 0U, 2215U, 2945U);
  fill(c2_enthalpy_correction,
       c2_candidate_enthalpy_value - c1_candidate_enthalpy_value);
  fill(c2_candidate_enthalpy, c2_candidate_enthalpy_value);
  fill(c2_candidate_temperature, 0.0);
  fill(c2_candidate_compressibility, 0.0);
  const PisoCoupledStateView final_state{
      final_velocity, final_pressure, final_enthalpy, final_density,
      final_temperature};
  Status c2_profile_status;
  for (std::int32_t zc = 0; zc < cells.z && c2_profile_status; ++zc)
    for (std::int32_t yc = 0; yc < cells.y && c2_profile_status; ++yc)
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const double target_density =
            final_density.unchecked(cell, 0U) +
            drho_dp.view.unchecked(cell, 0U) *
                coupled_pressure_correction.view.unchecked(cell, 0U);
        const double target_absolute_pressure =
            c1_candidate.closed_gauge.next_pressure_reference +
            final_pressure.unchecked(cell, 0U) +
            coupled_pressure_correction.view.unchecked(cell, 0U);
        if (!(target_density > 0.0) ||
            !(target_absolute_pressure > 0.0)) {
          c2_profile_status = {StatusCode::numerical_failure, 1U};
          break;
        }
        const double target_temperature =
            target_absolute_pressure /
            (exact_gas_constant * target_density);
        double candidate_h = 0.0;
        double candidate_cp = 0.0;
        double candidate_gas_constant = 0.0;
        c2_profile_status = fixture.thermodynamics.mixture_enthalpy(
            target_temperature, {}, candidate_h, candidate_cp,
            candidate_gas_constant);
        if (!c2_profile_status) break;
        c2_candidate_enthalpy.view.unchecked(cell, 0U) = candidate_h;
        c2_enthalpy_correction.view.unchecked(cell, 0U) =
            candidate_h - final_enthalpy.unchecked(cell, 0U);
      }
  passed &= expect(static_cast<bool>(c2_profile_status),
                   "C2 exact profile preserves its assembled density target");
  PisoExactEosClosureIdentity c2_closure;
  c2_closure.thermodynamics =
      pressure_input.pressure_reference.thermodynamics;
  c2_closure.pressure_reference =
      pressure_input.pressure_reference.pressure_reference;
  c2_closure.composition = c1_closure.composition;
  c2_closure.pressure_state =
      make_piso_field_revision_identity(as_const(final_pressure));
  c2_closure.pressure_correction = make_piso_field_revision_identity(
      as_const(coupled_pressure_correction.view));
  c2_closure.enthalpy_state =
      make_piso_field_revision_identity(as_const(final_enthalpy));
  c2_closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(c2_enthalpy_correction.view));
  c2_closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(c2_candidate_enthalpy.view));
  c2_closure.candidate_density = make_piso_field_revision_identity(
      as_const(c2_candidate_density.view));
  c2_closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(c2_candidate_temperature.view));
  c2_closure.closure = 2944U;
  PisoExactThermodynamicCandidateView c2_candidate;
  c2_candidate.enthalpy = as_const(c2_candidate_enthalpy.view);
  c2_candidate.density = as_const(c2_candidate_density.view);
  c2_candidate.temperature = as_const(c2_candidate_temperature.view);
  c2_candidate.closure = c2_closure;
  c2_candidate.pressure_compressibility =
      as_const(c2_candidate_compressibility.view);
  c2_candidate.pressure_compressibility.field = drho_dp.view.field;
  passed &= expect(static_cast<bool>(fill_exact_thermodynamic_candidate(
                       fixture.thermodynamics,
                       c1_candidate.closed_gauge.next_pressure_reference,
                       as_const(final_pressure),
                       as_const(coupled_pressure_correction.view),
                       as_const(c2_candidate_enthalpy.view),
                       as_const(final_velocity), {},
                       c2_candidate_density.view,
                       c2_candidate_temperature.view,
                       c2_candidate_compressibility.view)),
                   "C2 exact candidate is replayed from p/h/Y/U");
  passed &= expect(
      static_cast<bool>(prepare_closed_gauge(
          fixture.equations.pressure_reference(), solve_pressure_two,
          pressure_input.pressure_reference,
          c1_candidate.closed_gauge.next_pressure_reference,
          as_const(final_pressure),
          as_const(coupled_pressure_correction.view),
          c2_candidate.pressure_compressibility,
          c2_candidate.closure.closure, reductions,
          c2_candidate.closed_gauge)),
      "C2 exact state prepares the successor closed-gauge transition");
  PisoStateCorrectionCertificate corrected_two;
  passed &= expect(
      static_cast<bool>(coupler.correct_coupled_pending_state(
          solve_pressure_two, coupled_pressure_correction.view,
          as_const(c2_enthalpy_correction.view), final_state, c2_candidate,
          pending_final_flux, reductions, corrected_two)) &&
          corrected_two.valid() && corrected_two.corrector == 2U &&
          corrected_two.face_flux == pending_final_flux.revision() &&
          corrected_two.closure == PisoStateClosure::exact_eos &&
          corrected_two.enthalpy == final_enthalpy.revision &&
          corrected_two.temperature == final_temperature.revision &&
          corrected_two.exact_eos_closure !=
              corrected_one.exact_eos_closure,
      "coupled corrector two atomically writes U/pi/h/rho/T and pressure-equation pending flux under a fresh exact-EOS closure");
  ConstFaceFluxView corrected_pending_read;
  passed &= expect(
      static_cast<bool>(coupler.inspect_corrected_pending(
          corrected_two, pending_final_flux, corrected_pending_read)) &&
          corrected_pending_read.revision == corrected_two.face_flux &&
          corrected_pending_read.certificate.matches(corrected_pending_read),
      "C2 corrected pending flux is readable before terminal audit");
  PendingFaceFluxView foreign_pending_flux;
  passed &= expect(
      static_cast<bool>(foreign_flux_history.transaction.begin(
          foreign_flux_history.layers)) &&
          static_cast<bool>(foreign_flux_history.transaction.revise_trial(
              foreign_flux_history.dependency)) &&
          static_cast<bool>(foreign_flux_history.writer.begin_pending(
              foreign_flux_history.transaction, foreign_flux_history.storage,
              foreign_pending_flux)),
      "foreign pending flux acquires independent writer/storage authority");
  ConstFaceFluxView rejected_corrected_read;
  passed &= expect(
      coupler.inspect_corrected_pending(corrected_two, foreign_pending_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan &&
          rejected_corrected_read.revision == 0U &&
          foreign_flux_history.transaction.collective_finish(
              MPI_COMM_SELF, {StatusCode::rejected_step, 2203U})
                  .code == StatusCode::rejected_step,
      "foreign pending writer/storage cannot use the C2 correction authority");
  PisoStateCorrectionCertificate stale_correction = corrected_two;
  ++stale_correction.state;
  passed &= expect(
      coupler.inspect_corrected_pending(stale_correction,
                                        pending_final_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan &&
          rejected_corrected_read.revision == 0U,
      "stale C2 correction state cannot inspect pending flux");
  stale_correction = corrected_two;
  ++stale_correction.face_flux;
  passed &= expect(
      coupler.inspect_corrected_pending(stale_correction,
                                        pending_final_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan,
      "foreign C2 flux revision cannot inspect pending flux");
  stale_correction = corrected_two;
  ++stale_correction.enthalpy;
  passed &= expect(
      coupler.inspect_corrected_pending(stale_correction,
                                        pending_final_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan,
      "stale C2 enthalpy revision cannot inspect pending flux");
  stale_correction = corrected_two;
  ++stale_correction.exact_eos_closure;
  passed &= expect(
      coupler.inspect_corrected_pending(stale_correction,
                                        pending_final_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan,
      "foreign C2 exact-EOS closure cannot inspect pending flux");
  stale_correction = corrected_two;
  stale_correction.corrector = 1U;
  passed &= expect(
      coupler.inspect_corrected_pending(stale_correction,
                                        pending_final_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan,
      "non-C2 correction cannot inspect pending flux");
  double maximum_cumulative_pressure_error = 0.0;
  double maximum_incremental_velocity_error = 0.0;
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const std::size_t offset = cell_offset(cell);
        maximum_cumulative_pressure_error = std::max(
            maximum_cumulative_pressure_error,
            std::abs(final_pressure.unchecked(cell, 0U) -
                     (c1_pressure[offset] +
                      coupled_pressure_correction.view.unchecked(cell, 0U) -
                      c2_candidate.closed_gauge.shift)));
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double expected =
              c1_velocity[3U * offset + component] -
              r_au.view.unchecked(cell, component) *
                  pressure_gradient.view.unchecked(cell, component);
          maximum_incremental_velocity_error = std::max(
              maximum_incremental_velocity_error,
              std::abs(final_velocity.unchecked(cell, component) - expected));
        }
      }
    }
  }
  passed &= expect(
      maximum_cumulative_pressure_error < 1.0e-12,
      "corrector two accumulates delta-p onto the C1 pressure");
  passed &= expect(
      maximum_incremental_velocity_error < 1.0e-12,
      "corrector two applies delta-U to the C1 corrected velocity");
  const std::array final_dependencies{final_dependency};
  PisoTerminalAuditInput audit_input;
  audit_input.correction = corrected_two;
  audit_input.pressure_reference = corrected_two.output_pressure_reference;
  audit_input.density = as_const(final_density);
  audit_input.eos_density = as_const(final_density);
  audit_input.density_accepted = as_const(density_accepted.view);
  audit_input.density_previous = as_const(density_previous.view);
  audit_input.pressure_perturbation = as_const(final_pressure);
  audit_input.drho_dp_h_y = c2_candidate.pressure_compressibility;
  audit_input.bdf = input.bdf;
  audit_input.step_dt = time_step_for_bdf(input.bdf);
  audit_input.convective_cfl_limit = 1.0e6;
  audit_input.closed_mass_target = 0.0;
  for (std::int32_t zc = 0; zc < cells.z; ++zc)
    for (std::int32_t yc = 0; yc < cells.y; ++yc)
      for (std::int32_t xc = 0; xc < cells.x; ++xc)
        audit_input.closed_mass_target +=
            cell_volume(fixture, {xc, yc, zc}) *
            density.view.unchecked({xc, yc, zc}, 0U);
  audit_input.boundary_closure_residual = 0.0;

  OwnedField bad_eos =
      make_field(92U, cells, 1U, 0U, 2201U, 2910U);
  for (std::int32_t zc = 0; zc < cells.z; ++zc)
    for (std::int32_t yc = 0; yc < cells.y; ++yc)
      for (std::int32_t xc = 0; xc < cells.x; ++xc)
        bad_eos.view.unchecked({xc, yc, zc}, 0U) =
            final_density.unchecked({xc, yc, zc}, 0U) + 1.0e-5;
  PisoAttemptReport audit_report;
  PisoTerminalCertificate terminal;
  PisoTerminalAuditInput mutation_audit = audit_input;
  mutation_audit.correction.output_pressure_reference =
      corrected_two.input_pressure_reference;
  mutation_audit.correction.closed_gauge_rank_local_transaction =
      corrected_two.input_pressure_reference.pressure_reference;
  mutation_audit.pressure_reference =
      corrected_two.input_pressure_reference;
  passed &= expect(
      mutation_audit.correction.valid() &&
          coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                      reductions, audit_report, terminal)
                  .code == StatusCode::invalid_plan &&
          !terminal.valid() && pending_final_flux.valid(),
      "terminal audit rejects a self-consistent forged correction/audit pair carrying the old pressure reference");

  mutation_audit = audit_input;
  mutation_audit.pressure_reference =
      corrected_two.input_pressure_reference;
  terminal = {};
  passed &= expect(
      corrected_two.input_pressure_reference.pressure_reference !=
              corrected_two.output_pressure_reference.pressure_reference &&
          coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                      reductions, audit_report, terminal)
                  .code == StatusCode::invalid_plan &&
          !terminal.valid() && pending_final_flux.valid(),
      "terminal audit rejects the pre-C2 pressure reference without consuming pending flux");

  OwnedField foreign_terminal_chi = make_field(
      c2_candidate.pressure_compressibility.field, cells, 1U, 0U,
      c2_candidate.pressure_compressibility.revision, 3945U);
  for (std::int32_t zc = 0; zc < cells.z; ++zc)
    for (std::int32_t yc = 0; yc < cells.y; ++yc)
      for (std::int32_t xc = 0; xc < cells.x; ++xc)
        foreign_terminal_chi.view.unchecked({xc, yc, zc}, 0U) =
            c2_candidate.pressure_compressibility.unchecked(
                {xc, yc, zc}, 0U);
  mutation_audit = audit_input;
  mutation_audit.drho_dp_h_y = as_const(foreign_terminal_chi.view);
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::invalid_plan &&
          !terminal.valid() && pending_final_flux.valid(),
      "terminal audit rejects same-revision compressibility from foreign storage");

  std::vector<std::uint8_t> foreign_terminal_activity(cell_count, 1U);
  mutation_audit = audit_input;
  mutation_audit.active = {foreign_terminal_activity.data(),
                           foreign_terminal_activity.size()};
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::invalid_plan &&
          !terminal.valid() && pending_final_flux.valid(),
      "terminal audit rejects an activity mask not bound by the C2 gauge transaction");

  mutation_audit = audit_input;
  mutation_audit.eos_density = as_const(bad_eos.view);
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::rejected_step &&
          audit_report.eos_residual > valid_spec().eos_tolerance &&
          !terminal.valid() && pending_final_flux.valid(),
      "EOS terminal mutation rejects without consuming pending flux");

  OwnedField bad_accepted =
      make_field(93U, cells, 1U, 0U, 2202U, 2911U);
  fill(bad_accepted, 0.9001);
  mutation_audit = audit_input;
  mutation_audit.density_accepted = as_const(bad_accepted.view);
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::rejected_step &&
          audit_report.continuity_residual >
              valid_spec().continuity_tolerance &&
          !terminal.valid() && pending_final_flux.valid(),
      "continuity terminal mutation rejects without consuming pending flux");

  mutation_audit = audit_input;
  mutation_audit.energy_residual = 1.0e-5;
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::rejected_step &&
          audit_report.energy_residual > valid_spec().energy_tolerance &&
          !terminal.valid() && pending_final_flux.valid(),
      "energy terminal mutation rejects without consuming pending flux");

  mutation_audit = audit_input;
  mutation_audit.closed_mass_target = audit_input.closed_mass_target * 1.01;
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::rejected_step &&
          audit_report.closed_mass_residual >
              valid_spec().closed_mass_tolerance &&
          !terminal.valid() && pending_final_flux.valid(),
      "closed-mass terminal mutation rejects without consuming pending flux");

  OwnedField bad_pressure = make_field(
      94U, cells, 1U, 0U, final_pressure.revision, 2912U);
  fill(bad_pressure, 1.0e-4);
  mutation_audit = audit_input;
  mutation_audit.pressure_perturbation = as_const(bad_pressure.view);
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(mutation_audit, pending_final_flux,
                                  reductions, audit_report, terminal)
                  .code == StatusCode::rejected_step &&
          audit_report.gauge_residual > valid_spec().gauge_tolerance &&
          !terminal.valid() && pending_final_flux.valid(),
      "gauge terminal mutation rejects without consuming pending flux");

  passed &= expect(
      coupler.publish_pending_final(
                  PisoTerminalCertificate{},
                  {final_dependencies.data(), final_dependencies.size()},
                  {}, reductions, final_flux_writer, pending_final_flux)
                  .code == StatusCode::invalid_plan &&
          pending_final_flux.valid(),
      "pending final flux cannot publish without a terminal certificate");
  audit_report = {};
  terminal = {};
  const Status accepted_audit = coupler.audit_pending_final(
      audit_input, pending_final_flux, reductions, audit_report, terminal);
  if (!accepted_audit) {
    std::cerr << "terminal audit status="
              << static_cast<unsigned>(accepted_audit.code)
              << " eos=" << audit_report.eos_residual
              << " continuity=" << audit_report.continuity_residual
              << " mass=" << audit_report.closed_mass_residual
              << " gauge=" << audit_report.gauge_residual << '\n';
  }
  passed &= expect(
      static_cast<bool>(accepted_audit) &&
          terminal.valid() &&
          terminal.final_flux == pending_final_flux.revision() &&
          terminal.pressure_reference.plan ==
              corrected_two.output_pressure_reference.plan &&
          terminal.pressure_reference.predictor ==
              corrected_two.output_pressure_reference.predictor &&
          terminal.pressure_reference.thermodynamics ==
              corrected_two.output_pressure_reference.thermodynamics &&
          terminal.pressure_reference.closure ==
              corrected_two.output_pressure_reference.closure &&
          terminal.pressure_reference.time ==
              corrected_two.output_pressure_reference.time &&
          terminal.pressure_reference.pressure_reference ==
              corrected_two.output_pressure_reference.pressure_reference &&
          terminal.pressure_reference.kind ==
              corrected_two.output_pressure_reference.kind &&
          terminal.closed_gauge_collective_transaction ==
              corrected_two.closed_gauge_collective_transaction &&
          audit_report.eos_residual <= valid_spec().eos_tolerance &&
          audit_report.continuity_residual <=
              valid_spec().continuity_tolerance &&
          audit_report.energy_residual <= valid_spec().energy_tolerance &&
          audit_report.closed_mass_residual <=
              valid_spec().closed_mass_tolerance &&
          audit_report.gauge_residual <= valid_spec().gauge_tolerance,
      "all five terminal gates certify the pending final state");
  PisoTerminalCertificate stale_terminal = terminal;
  stale_terminal.pressure_reference =
      corrected_two.input_pressure_reference;
  ConstFaceFluxView rejected_pending;
  passed &= expect(
      coupler.inspect_pending_final(stale_terminal, pending_final_flux,
                                    rejected_pending)
                  .code == StatusCode::invalid_plan &&
          rejected_pending.revision == 0U && pending_final_flux.valid(),
      "terminal inspection rejects the old pressure-reference authority");
  stale_terminal = terminal;
  ++stale_terminal.closed_gauge_collective_transaction;
  passed &= expect(
      coupler.publish_pending_final(
                  stale_terminal,
                  {final_dependencies.data(), final_dependencies.size()},
                  {}, reductions, final_flux_writer, pending_final_flux)
                  .code == StatusCode::invalid_plan &&
          pending_final_flux.valid(),
      "terminal publication rejects a mutated collective gauge token");
  terminal = {};
  passed &= expect(
      static_cast<bool>(coupler.audit_pending_final(
          audit_input, pending_final_flux, reductions, audit_report,
          terminal)) &&
          terminal.valid(),
      "terminal authority is recertified after a rejected forged publication");
  ConstFaceFluxView inspected_pending;
  passed &= expect(
      static_cast<bool>(coupler.inspect_pending_final(
          terminal, pending_final_flux, inspected_pending)) &&
          inspected_pending.certificate.matches(inspected_pending) &&
          inspected_pending.revision == terminal.final_flux,
      "terminal audit exposes a certified read-only pending flux view");
  double maximum_incremental_flux_error =
      std::numeric_limits<double>::infinity();
  if (inspected_pending.revision != 0U &&
      inspected_pending.x.base != nullptr &&
      inspected_pending.y.base != nullptr &&
      inspected_pending.z.base != nullptr) {
    maximum_incremental_flux_error = 0.0;
    const std::array<CartesianAxis, 3U> axes{
        CartesianAxis::x, CartesianAxis::y, CartesianAxis::z};
    const std::array<ConstFaceFieldView, 3U> final_faces{
        inspected_pending.x, inspected_pending.y, inspected_pending.z};
    const std::array<ConstFaceFieldView, 3U> pressure_coefficients{
        as_const(ax.view), as_const(ay.view), as_const(az.view)};
    for (std::size_t axis_index = 0U; axis_index < axes.size();
         ++axis_index) {
      const Int3 extents = final_faces[axis_index].extents;
      for (std::int32_t zf = 0; zf < extents.z; ++zf) {
        for (std::int32_t yf = 0; yf < extents.y; ++yf) {
          for (std::int32_t xf = 0; xf < extents.x; ++xf) {
            const Int3 face{xf, yf, zf};
            Int3 left = face;
            if (axes[axis_index] == CartesianAxis::x) {
              --left.x;
            } else if (axes[axis_index] == CartesianAxis::y) {
              --left.y;
            } else {
              --left.z;
            }
            const double correction_jump =
                coupled_pressure_correction.view.unchecked(face, 0U) -
                coupled_pressure_correction.view.unchecked(left, 0U);
            const double expected =
                c1_flux[axis_index][face_offset(extents, face)] -
                pressure_coefficients[axis_index].unchecked(face) *
                    correction_jump;
            maximum_incremental_flux_error = std::max(
                maximum_incremental_flux_error,
                std::abs(final_faces[axis_index].unchecked(face) - expected));
          }
        }
      }
    }
  }
  passed &= expect(
      maximum_incremental_flux_error < 1.0e-12,
      "corrector two applies delta-phi to the C1 corrected trial flux");
  ConstFaceFluxView rejected_inspection;
  passed &= expect(
      coupler.inspect_pending_final(PisoTerminalCertificate{},
                                    pending_final_flux,
                                    rejected_inspection)
                  .code == StatusCode::invalid_plan &&
          rejected_inspection.revision == 0U,
      "pending flux inspection requires the exact terminal certificate");
  passed &= expect(
      static_cast<bool>(coupler.publish_pending_final(
          terminal,
          {final_dependencies.data(), final_dependencies.size()},
          {}, reductions, final_flux_writer, pending_final_flux)) &&
          static_cast<bool>(transaction.collective_finish(MPI_COMM_SELF,
                                                           Status{})),
      "pending pressure-equation flux publishes only at collective commit");
  rejected_corrected_read = {};
  passed &= expect(
      coupler.inspect_corrected_pending(corrected_two, pending_final_flux,
                                        rejected_corrected_read)
                  .code == StatusCode::invalid_plan &&
          rejected_corrected_read.revision == 0U,
      "consumed pending flux cannot be inspected through C2 authority");
  ConstFaceFluxView committed_final_flux;
  const bool committed_flux_ok =
      static_cast<bool>(final_flux_writer.committed(final_flux_storage,
                                                    committed_final_flux)) &&
          committed_final_flux.revision == corrected_two.face_flux &&
          committed_final_flux.certificate.valid();
  passed &= expect(committed_flux_ok,
      "committed final flux names the exact corrector-two revision");
  if (!committed_flux_ok) {
    return false;
  }
  double maximum_final_continuity = 0.0;
  for (std::int32_t zc = 0; zc < cells.z; ++zc) {
    for (std::int32_t yc = 0; yc < cells.y; ++yc) {
      for (std::int32_t xc = 0; xc < cells.x; ++xc) {
        const Int3 cell{xc, yc, zc};
        const double density_defect =
            (input.bdf.a0 * final_density.unchecked(cell, 0U) +
             input.bdf.a1 * density_accepted.view.unchecked(cell, 0U) +
             input.bdf.a2 * density_previous.view.unchecked(cell, 0U)) *
            cell_volume(fixture, cell);
        const double flux_divergence =
            committed_final_flux.x.unchecked({xc + 1, yc, zc}) -
            committed_final_flux.x.unchecked(cell) +
            committed_final_flux.y.unchecked({xc, yc + 1, zc}) -
            committed_final_flux.y.unchecked(cell) +
            committed_final_flux.z.unchecked({xc, yc, zc + 1}) -
            committed_final_flux.z.unchecked(cell);
        maximum_final_continuity = std::max(
            maximum_final_continuity,
            std::abs(density_defect + flux_divergence));
      }
    }
  }
  passed &= expect(
      maximum_final_continuity < 2.0e-8,
      "committed final flux, not reconstructed final U, closes continuity");
  const Status finalized = retry_epoch.finalize(audit_report);
  passed &= expect(finalized.code == StatusCode::invalid_plan &&
                       audit_report.pressure_solve_calls == 0U &&
                       audit_report.continuity_residual <=
                           valid_spec().continuity_tolerance &&
                       audit_report.final_flux_revision ==
                           committed_final_flux.revision &&
                       pressure_mg.counters().numeric_refreshes == 1U &&
                       krylov_workspace.vector_storage_address() ==
                           krylov_address &&
                       pressure_mg.workspace_storage_address() == mg_address,
                   "a rejected pressure-only C2 epoch cannot overwrite the independently certified coupled terminal report or workspaces");
  passed &= expect(
      retry_epoch
                  .solve(piso, 2U, solve_pressure_two, solve_identity_two,
                         coefficients_two, coupler, solve_operator,
                         pressure_mg, pressure_system,
                         pressure_correction.view, krylov_workspace,
                         reductions)
                  .code == StatusCode::invalid_plan,
      "a third pressure solve is impossible after epoch finalization");

  PisoIntermediateInput missing_ghost_authority = input;
  missing_ghost_authority.density.field = 99U;
  passed &= expect(coupler.refresh(missing_ghost_authority, certificate).code ==
                       StatusCode::invalid_plan,
                   "density halo requires the cold-bound field authority");
  return passed;
}

bool test_open_boundary_terminal_authority() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, true),
                       "open-boundary PISO fixture compiles");
  if (!passed) {
    return false;
  }
  CartesianMeshSpec foreign_mesh = mesh_spec();
  foreign_mesh.upper.x = 4.0;
  foreign_mesh.minimum_spacing.x = 1.0;
  CartesianGeometryPlan foreign_geometry;
  MeshPatch foreign_patch;
  PressureCorrectionBoundaryPlan foreign_pressure_boundary;
  passed &= expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, foreign_mesh, {}, foreign_geometry, foreign_patch)) &&
          foreign_geometry.fingerprint() != fixture.geometry.fingerprint(),
      "same-size foreign metric geometry compiles distinctly");
  passed &= expect(
      PressureCorrectionBoundaryPlan::compile(
          foreign_geometry, foreign_patch, fixture.boundary,
          foreign_pressure_boundary)
              .code == StatusCode::invalid_plan &&
          !foreign_pressure_boundary.certificate().valid(),
      "pressure boundary rejects a same-size foreign source geometry");
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "open-boundary PISO plan compiles");
  const Int3 cells = fixture.patch.cells;
  const std::uint8_t reach = fixture.equations.kernels().reach();
  OwnedField density = make_field(0U, cells, 1U, reach, 3001U, 4001U);
  OwnedField r_au = make_field(40U, cells, 3U, 1U, 3002U, 4002U);
  OwnedField h_by_a = make_field(41U, cells, 3U, reach, 3003U, 4003U);
  OwnedField gradient = make_field(42U, cells, 3U, 0U, 3004U, 4004U);
  OwnedFaceField ax = make_face(CartesianAxis::x, cells, 4005U);
  OwnedFaceField ay = make_face(CartesianAxis::y, cells, 4006U);
  OwnedFaceField az = make_face(CartesianAxis::z, cells, 4007U);
  FaceFluxStorage phi_storage;
  FaceFluxView phi;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(cells, 1U,
                                                            phi_storage)) &&
          static_cast<bool>(phi_storage.workspace_view(0U, 3005U, phi)),
      "open-boundary phiHbyA workspace allocates");
  const PisoCouplerWorkspace workspace{r_au.view, h_by_a.view, gradient.view,
                                       ax.view, ay.view, az.view, phi};
  const std::array<HaloFieldSpec, 3U> halo_fields{{
      {density.view.field, 1U, 1U}, {r_au.view.field, 1U, 3U},
      {h_by_a.view.field, reach, 3U}}};
  HaloEngine halo;
  constexpr FieldId correction_field = 90U;
  const std::array<HaloFieldSpec, 1U> correction_fields{{
      {correction_field, 1U, 1U}}};
  HaloEngine correction_halo;
  passed &= expect(
      static_cast<bool>(halo.reserve(
          MPI_COMM_SELF, fixture.patch,
          {halo_fields.data(), halo_fields.size()},
          fixture.boundary.halo_topology())) &&
          static_cast<bool>(correction_halo.reserve(
              MPI_COMM_SELF, fixture.patch,
              {correction_fields.data(), correction_fields.size()},
              fixture.boundary.halo_topology())),
      "open-boundary coupler halos reserve");
  const PisoCouplerServices services{
      MPI_COMM_SELF, &fixture.geometry, fixture.patch, &fixture.boundary,
      &fixture.thermodynamics, &halo, 301U, density.view.field,
      &correction_halo, 302U,
      correction_field};
  IbmForceFixture donor_fixture;
  CandidateBoundaryFixture donor_interface_fixture;
  CandidateBoundaryFixtureSpec donor_interface_spec;
  donor_interface_spec.cells_per_axis = 16;
  donor_interface_spec.immersed = true;
  RemoteDonorExchangePlan live_pressure_donors;
  RemoteDonorExchangePlan candidate_pressure_donors;
  const std::array<RemoteDonorFieldSpec, 1U> live_pressure_fields{{
      {correction_field, 1U},
  }};
  constexpr FieldId candidate_pressure_field = 200U;
  const std::array<RemoteDonorFieldSpec, 1U> candidate_pressure_fields{{
      {candidate_pressure_field, 1U},
  }};
  constexpr StageId live_pressure_donor_stage = 303U;
  constexpr StageId candidate_pressure_donor_stage = 304U;
  passed &= expect(
      donor_fixture.initialize(MPI_COMM_SELF, 16) &&
          donor_interface_fixture.initialize(MPI_COMM_SELF,
                                              donor_interface_spec) &&
          static_cast<bool>(RemoteDonorExchangePlan::compile(
              MPI_COMM_SELF, donor_fixture.geometry.global_cells(),
              donor_fixture.patch, donor_fixture.boundary.reconstruction(),
              {live_pressure_fields.data(), live_pressure_fields.size()},
              live_pressure_donor_stage, live_pressure_donors)) &&
          static_cast<bool>(RemoteDonorExchangePlan::compile(
              MPI_COMM_SELF, donor_fixture.geometry.global_cells(),
              donor_fixture.patch, donor_fixture.boundary.reconstruction(),
              {candidate_pressure_fields.data(),
               candidate_pressure_fields.size()},
              candidate_pressure_donor_stage,
              candidate_pressure_donors)),
      "independent live/candidate IBM pressure donor plans compile");
  PisoCouplerServices donor_services = services;
  donor_services.pressure_correction_donors = &live_pressure_donors;
  donor_services.pressure_correction_donor_stage =
      live_pressure_donor_stage;
  donor_services.candidate_pressure_correction_donors =
      &candidate_pressure_donors;
  donor_services.candidate_pressure_correction_donor_stage =
      candidate_pressure_donor_stage;
  donor_services.candidate_pressure_correction_field =
      candidate_pressure_field;
  donor_services.candidate_pressure_correction_donor_reach =
      candidate_pressure_donors.reach();
  donor_services.candidate_pressure_correction_donor_fingerprint =
      candidate_pressure_donors.fingerprint();
  donor_services.immersed_interface =
      &donor_interface_fixture.immersed_interface;
  PressureVelocityCoupler donor_contract_coupler;
  passed &= expect(
      static_cast<bool>(PressureVelocityCoupler::bind(
          piso, fixture.equations, donor_services, workspace,
          donor_contract_coupler)),
      "candidate IBM donor contract binds distinct field/stage/reach/fingerprint");
  const auto rejects_candidate_donor_contract = [&](
      PisoCouplerServices mutation, std::string_view description) {
    PressureVelocityCoupler rejected;
    passed &= expect(
        PressureVelocityCoupler::bind(
            piso, fixture.equations, mutation, workspace, rejected)
                .code == StatusCode::invalid_plan,
        description);
  };
  PisoCouplerServices bad_candidate_donor = donor_services;
  ++bad_candidate_donor.candidate_pressure_correction_donor_fingerprint;
  rejects_candidate_donor_contract(
      bad_candidate_donor,
      "candidate IBM donor rejects a foreign plan fingerprint");
  bad_candidate_donor = donor_services;
  bad_candidate_donor.candidate_pressure_correction_donors =
      &live_pressure_donors;
  rejects_candidate_donor_contract(
      bad_candidate_donor,
      "candidate IBM donor rejects the live donor plan instance");
  bad_candidate_donor = donor_services;
  bad_candidate_donor.candidate_pressure_correction_field =
      correction_field;
  rejects_candidate_donor_contract(
      bad_candidate_donor,
      "candidate IBM donor rejects the live correction FieldId");
  bad_candidate_donor = donor_services;
  bad_candidate_donor.candidate_pressure_correction_donor_stage =
      live_pressure_donor_stage;
  rejects_candidate_donor_contract(
      bad_candidate_donor,
      "candidate IBM donor rejects the live correction StageId");
  bad_candidate_donor = donor_services;
  ++bad_candidate_donor.candidate_pressure_correction_donor_reach;
  rejects_candidate_donor_contract(
      bad_candidate_donor,
      "candidate IBM donor rejects a foreign nonzero reach");
  PisoCouplerServices foreign_geometry_services = services;
  foreign_geometry_services.geometry = &foreign_geometry;
  foreign_geometry_services.patch = foreign_patch;
  PressureVelocityCoupler foreign_geometry_coupler;
  passed &= expect(
      PressureVelocityCoupler::bind(
          piso, fixture.equations, foreign_geometry_services, workspace,
          foreign_geometry_coupler)
                  .code == StatusCode::invalid_plan,
      "PISO bind rejects a same-size foreign metric geometry");
  PressureVelocityCoupler coupler;
  passed &= expect(static_cast<bool>(PressureVelocityCoupler::bind(
                       piso, fixture.equations, services, workspace,
                       coupler)),
                   "open-boundary coupler binds");
  PressureEnergyCandidateBoundaryFinalizerBinding finalizer_binding;
  finalizer_binding.communicator = MPI_COMM_SELF;
  finalizer_binding.geometry = &fixture.geometry;
  finalizer_binding.patch = fixture.patch;
  finalizer_binding.boundary = &fixture.boundary;
  finalizer_binding.kernels = &fixture.equations.kernels();
  finalizer_binding.thermodynamics = &fixture.thermodynamics;
  finalizer_binding.transport = &fixture.transport;
  finalizer_binding.coupler = &coupler;
  PressureEnergyCandidateBoundaryFinalizer open_finalizer;
  passed &= expect(
      static_cast<bool>(PressureEnergyCandidateBoundaryFinalizer::bind(
          finalizer_binding, open_finalizer)),
      "open candidate boundary finalizer cold-binds independently of its coupler issuer");
  Fixture foreign_finalizer_fixture;
  PressureEnergyCandidateBoundaryFinalizer foreign_finalizer;
  PressureEnergyCandidateBoundaryFinalizerBinding foreign_finalizer_binding =
      finalizer_binding;
  passed &= expect(
      make_fixture(foreign_finalizer_fixture, true) &&
          (foreign_finalizer_binding.geometry =
               &foreign_finalizer_fixture.geometry,
           foreign_finalizer_binding.patch = foreign_finalizer_fixture.patch,
           foreign_finalizer_binding.boundary =
               &foreign_finalizer_fixture.boundary,
           foreign_finalizer_binding.kernels =
               &foreign_finalizer_fixture.equations.kernels(),
           foreign_finalizer_binding.thermodynamics =
               &foreign_finalizer_fixture.thermodynamics,
           foreign_finalizer_binding.transport =
               &foreign_finalizer_fixture.transport,
           true) &&
          PressureEnergyCandidateBoundaryFinalizer::bind(
              foreign_finalizer_binding, foreign_finalizer)
                  .code == StatusCode::invalid_plan &&
          !foreign_finalizer.ready(),
      "candidate boundary finalizer rejects same-shape foreign geometry/boundary/thermo identities for a claimed coupler");

  OwnedField boundary_pressure = make_field(
      fixture.boundary.pressure_field(), cells, 1U, reach, 3008U, 4090U);
  OwnedField boundary_enthalpy = make_field(
      fixture.boundary.enthalpy_field(), cells, 1U, reach, 3009U, 4091U);
  std::array<OwnedField, 7U> boundary_thermo_aux{
      make_field(50U, cells, 1U, reach, 3050U, 4092U),
      make_field(51U, cells, 1U, reach, 3051U, 4093U),
      make_field(52U, cells, 1U, reach, 3052U, 4094U),
      make_field(53U, cells, 1U, reach, 3053U, 4095U),
      make_field(54U, cells, 1U, reach, 3054U, 4096U),
      make_field(55U, cells, 1U, reach, 3055U, 4097U),
      make_field(56U, cells, 1U, reach, 3056U, 4098U)};
  double boundary_enthalpy_value = 0.0;
  double boundary_cp = 0.0;
  double boundary_gas = 0.0;
  passed &= expect(
      static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
          500.0, {}, boundary_enthalpy_value, boundary_cp, boundary_gas)),
      "open-boundary thermophysical ghost state evaluates");
  fill(boundary_pressure, 0.0);
  fill(boundary_enthalpy, boundary_enthalpy_value);
  fill(boundary_thermo_aux[0U], 500.0);
  const ThermodynamicsPlan* issuing_thermodynamics =
      &fixture.thermodynamics;
  const TransportPlan* issuing_transport = &fixture.transport;
  const auto certify_boundary_thermophysics = [&](
      BoundaryThermophysicalGhostContext context,
      BoundaryThermophysicalGhostCertificate& certificate) {
    const std::array<BoundaryGhostFieldAuthority, 2U> authorities{{
        make_boundary_ghost_field_authority(
            as_const(boundary_pressure.view)),
        make_boundary_ghost_field_authority(
            as_const(boundary_enthalpy.view)),
    }};
    const BoundaryThermophysicalGhostAuthority authority{
        static_cast<std::uintptr_t>(3100U +
                                    static_cast<std::uint8_t>(context.phase)),
        fixture.boundary.revision(),
        fixture.boundary.local_layout_fingerprint(), cells,
        {reach, reach, reach}, reach,
        {authorities.data(), authorities.size()}};
    return BoundaryThermophysicalFaceClosure::close(
        fixture.boundary, *issuing_thermodynamics, *issuing_transport,
        {101325.0, as_const(boundary_pressure.view),
         as_const(boundary_enthalpy.view), {}, authority},
        {density.view, boundary_thermo_aux[0U].view,
         boundary_thermo_aux[1U].view, boundary_thermo_aux[2U].view,
         boundary_thermo_aux[3U].view, boundary_thermo_aux[4U].view,
         boundary_thermo_aux[5U].view, boundary_thermo_aux[6U].view},
        context, certificate);
  };

  OwnedField velocity = make_field(1U, cells, 3U, reach, 3006U, 4008U);
  OwnedField diagonal = make_field(30U, cells, 3U, 0U, 3016U, 4016U);
  OwnedField rhs = make_field(31U, cells, 3U, 0U, 3017U, 4017U);
  FaceFluxStorage paired_storage;
  FaceFluxView paired_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, paired_storage)) &&
          static_cast<bool>(paired_storage.workspace_view(0U, 3018U,
                                                          paired_flux)),
      "open-boundary paired-flux authority allocates");
  fill(density, 1.0);
  // The thermophysical boundary closure may produce a non-owner EOS density
  // at a pressure/enthalpy Dirichlet face.  Preserve a deliberately distinct
  // high-x ghost so the pressure coefficient cannot silently fall back to an
  // owner zero-gradient closure.
  for (std::int32_t zc = 0; zc < cells.z; ++zc)
    for (std::int32_t yc = 0; yc < cells.y; ++yc)
      density.view.unchecked({cells.x, yc, zc}, 0U) = 4.0;
  fill(velocity, 0.0);
  const BdfCoefficients boundary_bdf{1.5, -2.0, 0.5, 2U};
  fill(diagonal, boundary_bdf.a0 * cell_volume(fixture, {0, 0, 0}));
  fill(rhs, 0.0);
  fill_face_flux(phi, 0.0);
  fill_face_flux(paired_flux, 9.0);
  CommittedFluxHistory boundary_history;
  ConstFaceFluxView boundary_previous;
  ConstFaceFluxView boundary_accepted;
  passed &= expect(
      initialize_flux_history(cells, boundary_history) &&
          commit_zero_flux(fixture.equations.kernels(), cells, 3020U,
                           boundary_history, boundary_previous) &&
          commit_zero_flux(fixture.equations.kernels(), cells, 3030U,
                           boundary_history, boundary_accepted),
      "open-boundary committed face-history layers publish");
  PisoIntermediateInput boundary_input;
  boundary_input.momentum = {
      fixture.equations.momentum().fingerprint(),
      EquationAssemblyScope::momentum_predictor,
      3040U,
      fixture.geometry.topology_revision(),
      paired_flux.revision,
      3041U,
      time_step_for_bdf(boundary_bdf)};
  boundary_input.predictor.plan =
      fixture.equations.thermophysical_predictor().fingerprint();
  boundary_input.predictor.time = boundary_input.momentum.time;
  boundary_input.predictor.geometry = boundary_input.momentum.geometry;
  boundary_input.predictor.accepted_face_flux = boundary_accepted.revision;
  boundary_input.predictor.previous_face_flux = boundary_previous.revision;
  boundary_input.predictor.committed_face_flux_authority =
      boundary_accepted.certificate.authority();
  boundary_input.predictor.committed_face_flux_storage =
      boundary_accepted.certificate.storage();
  boundary_input.predictor.committed_face_flux_revision_domain =
      boundary_accepted.certificate.revision_domain();
  boundary_input.predictor.predicted_density = density.view.revision + 1U;
  boundary_input.predictor.predicted_density_storage =
      density.view.storage_identity;
  boundary_input.predictor.predicted_density_revision_domain =
      density.view.revision_domain;
  boundary_input.predictor.paired_face_flux = paired_flux.revision;
  boundary_input.predictor.paired_face_flux_storage =
      paired_flux.x.storage_identity;
  boundary_input.predictor.paired_face_flux_revision_domain =
      paired_flux.x.revision_domain;
  boundary_input.predictor.state = 3042U;
  boundary_input.predictor.order = boundary_bdf.order;
  boundary_input.pressure_reference = {
      fixture.equations.pressure_reference().fingerprint(),
      fixture.equations.thermophysical_predictor().fingerprint(),
      fixture.thermodynamics.fingerprint(),
      3043U,
      boundary_input.momentum.time,
      3044U,
      PressureReferenceKind::boundary_absolute};
  boundary_input.density = density.view;
  boundary_input.trial_velocity = as_const(velocity.view);
  boundary_input.trial_flux = as_const(paired_flux);
  boundary_input.momentum_system = {diagonal.view, rhs.view};
  boundary_input.bdf = boundary_bdf;
  boundary_input.numeric_boundary = fixture.boundary.revision();
  boundary_input.corrector = 1U;
  boundary_input.temporal_reference = as_const(phi);
  boundary_input.committed_face_history = {boundary_accepted,
                                           boundary_previous};
  BoundaryThermophysicalGhostCertificate boundary_c1_thermophysics;
  const BoundaryThermophysicalGhostContext boundary_c1_context{
      boundary_input.momentum.time, fixture.geometry.fingerprint(),
      boundary_input.pressure_reference.pressure_reference,
      boundary_input.numeric_boundary,
      BoundaryThermophysicalGhostPhase::corrector_one};
  passed &= expect(
      static_cast<bool>(certify_boundary_thermophysics(
          boundary_c1_context, boundary_c1_thermophysics)),
      "open-boundary C1 physical ghosts are certified");
  boundary_input.thermophysical_boundary = {
      boundary_c1_thermophysics,
      {101325.0, as_const(boundary_pressure.view),
       as_const(boundary_enthalpy.view), {}, as_const(density.view)}};
  PisoIntermediateCertificate boundary_certificate;
  const Status boundary_refresh =
      coupler.refresh(boundary_input, boundary_certificate);
  const bool boundary_authority_ok =
      static_cast<bool>(boundary_refresh) &&
          boundary_certificate.valid() &&
          std::abs(phi.x.unchecked({0, 0, 0}) - 9.0) < 1.0e-12 &&
          std::abs(phi.x.unchecked({1, 0, 0})) < 1.0e-12 &&
          std::abs(phi.x.unchecked({cells.x, 0, 0})) < 1.0e-12;
  if (!boundary_authority_ok) {
    std::cerr << "open C1 status="
              << static_cast<unsigned>(boundary_refresh.code) << '/'
              << boundary_refresh.detail << " phi="
              << phi.x.unchecked({0, 0, 0}) << ','
              << phi.x.unchecked({1, 0, 0}) << ','
              << phi.x.unchecked({cells.x, 0, 0}) << '\n';
  }
  passed &= expect(
      boundary_authority_ok,
      "C1 keeps paired flux only at fixed physical boundaries");
  const Int3 high_owner{cells.x - 1, 0, 0};
  const double high_distance =
      0.5 * cell_width(fixture, CartesianAxis::x, high_owner.x);
  const double owner_mobility = 1.0 / boundary_bdf.a0;
  const double ghost_mobility =
      density.view.unchecked({cells.x, 0, 0}, 0U) / boundary_bdf.a0;
  const double expected_high_coefficient =
      2.0 * face_area(fixture, CartesianAxis::x, high_owner) /
      (high_distance / owner_mobility + high_distance / ghost_mobility);
  const double owner_only_coefficient =
      face_area(fixture, CartesianAxis::x, high_owner) * owner_mobility /
      high_distance;
  passed &= expect(
      std::abs(ax.view.unchecked({cells.x, 0, 0}) - expected_high_coefficient) <
              1.0e-13 &&
          std::abs(ax.view.unchecked({cells.x, 0, 0}) -
                   owner_only_coefficient) > 1.0e-3,
      "pressure Dirichlet coefficient consumes the EOS-closed rho ghost");

  PisoIntermediateInput boundary_c2_input = boundary_input;
  boundary_c2_input.corrector = 2U;
  boundary_c2_input.prior_corrector = boundary_certificate.dependency;
  boundary_c2_input.temporal_reference = {};
  boundary_c2_input.committed_face_history = {};
  ConstFaceFluxView boundary_c2_flux = as_const(paired_flux);
  ++boundary_c2_flux.revision;
  boundary_c2_input.trial_flux = boundary_c2_flux;
  BoundaryThermophysicalGhostCertificate boundary_c2_thermophysics;
  BoundaryThermophysicalGhostContext boundary_c2_context =
      boundary_c1_context;
  boundary_c2_context.phase =
      BoundaryThermophysicalGhostPhase::corrector_two;
  passed &= expect(
      static_cast<bool>(certify_boundary_thermophysics(
          boundary_c2_context, boundary_c2_thermophysics)),
      "open-boundary C2 physical ghosts are independently certified");
  boundary_c2_input.thermophysical_boundary = {
      boundary_c2_thermophysics,
      {101325.0, as_const(boundary_pressure.view),
       as_const(boundary_enthalpy.view), {}, as_const(density.view)}};
  fill_face_flux(paired_flux, 0.0);
  PisoIntermediateCertificate boundary_c2_certificate;
  passed &= expect(
      static_cast<bool>(coupler.refresh(boundary_c2_input,
                                        boundary_c2_certificate)) &&
          boundary_c2_certificate.valid(),
      "open-boundary C1 to C2 phase chain advances");

  PressureVelocityCoupler mutation_coupler;
  passed &= expect(static_cast<bool>(PressureVelocityCoupler::bind(
                       piso, fixture.equations, services, workspace,
                       mutation_coupler)),
                   "independent open-boundary mutation coupler binds");

  ThermophysicalSpec foreign_thermophysical_spec = thermophysical_spec();
  foreign_thermophysical_spec.data_file = "foreign-issuer.d";
  foreign_thermophysical_spec.species[0U].conductivity *= 1.25;
  ThermodynamicsPlan foreign_thermodynamics;
  TransportPlan foreign_transport;
  passed &= expect(
      static_cast<bool>(ThermodynamicsPlan::compile(
          foreign_thermophysical_spec, {}, foreign_thermodynamics)) &&
          static_cast<bool>(TransportPlan::compile(
              foreign_thermophysical_spec, foreign_thermodynamics,
              foreign_transport)),
      "foreign thermodynamics and transport issuer compiles");
  issuing_thermodynamics = &foreign_thermodynamics;
  issuing_transport = &foreign_transport;
  BoundaryThermophysicalGhostCertificate foreign_issuer_certificate;
  passed &= expect(
      static_cast<bool>(certify_boundary_thermophysics(
          boundary_c1_context, foreign_issuer_certificate)) &&
          foreign_issuer_certificate.collective_semantics() !=
              boundary_c1_thermophysics.collective_semantics() &&
          foreign_issuer_certificate.collective_lineage() !=
              boundary_c1_thermophysics.collective_lineage(),
      "foreign thermophysical issuer changes collective semantics/lineage");
  boundary_input.thermophysical_boundary.certificate =
      foreign_issuer_certificate;
  PisoIntermediateCertificate foreign_issuer_rejected;
  passed &= expect(
      mutation_coupler.refresh(boundary_input, foreign_issuer_rejected).code ==
              StatusCode::invalid_plan &&
          !foreign_issuer_rejected.valid(),
      "PISO rejects a foreign thermodynamics/transport ghost issuer");
  issuing_thermodynamics = &fixture.thermodynamics;
  issuing_transport = &fixture.transport;
  boundary_input.thermophysical_boundary.certificate =
      boundary_c1_thermophysics;

  const Int3 certified_physical_ghost{cells.x, 0, 0};
  const auto rejects_boundary_mutation = [&](
      std::string_view description) {
    PisoIntermediateCertificate rejected;
    passed &= expect(
        mutation_coupler.refresh(boundary_input, rejected).code ==
                StatusCode::invalid_plan &&
            !rejected.valid(),
        description);
  };
  const double certified_pressure =
      boundary_pressure.view.unchecked(certified_physical_ghost, 0U);
  boundary_pressure.view.unchecked(certified_physical_ghost, 0U) += 1.0;
  rejects_boundary_mutation(
      "same-revision physical p ghost mutation fails closed in PISO");
  boundary_pressure.view.unchecked(certified_physical_ghost, 0U) =
      certified_pressure;
  const double certified_enthalpy =
      boundary_enthalpy.view.unchecked(certified_physical_ghost, 0U);
  boundary_enthalpy.view.unchecked(certified_physical_ghost, 0U) += 1.0;
  rejects_boundary_mutation(
      "same-revision physical h ghost mutation fails closed in PISO");
  boundary_enthalpy.view.unchecked(certified_physical_ghost, 0U) =
      certified_enthalpy;
  const double certified_density =
      density.view.unchecked(certified_physical_ghost, 0U);
  density.view.unchecked(certified_physical_ghost, 0U) *= 1.01;
  rejects_boundary_mutation(
      "same-revision physical rho ghost mutation fails closed in PISO");
  density.view.unchecked(certified_physical_ghost, 0U) = certified_density;
  OwnedField foreign_boundary_pressure = make_field(
      boundary_pressure.view.field, cells, 1U, reach,
      boundary_pressure.view.revision,
      boundary_pressure.view.storage_identity);
  foreign_boundary_pressure.view.revision_domain =
      boundary_pressure.view.revision_domain;
  foreign_boundary_pressure.view.replica =
      boundary_pressure.view.replica + 1U;
  fill(foreign_boundary_pressure, 0.0);
  boundary_input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(foreign_boundary_pressure.view);
  rejects_boundary_mutation(
      "foreign physical p base/replica fails closed in PISO");
  boundary_input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(boundary_pressure.view);
  ++boundary_input.thermophysical_boundary.binding.pressure_perturbation
        .storage_identity;
  rejects_boundary_mutation(
      "foreign physical p storage fails closed in PISO");
  boundary_input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(boundary_pressure.view);
  ++boundary_input.thermophysical_boundary.binding.pressure_perturbation
        .revision_domain;
  rejects_boundary_mutation(
      "foreign physical p revision domain fails closed in PISO");
  boundary_input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(boundary_pressure.view);
  const auto rejects_foreign_context = [&](
      BoundaryThermophysicalGhostContext context,
      std::string_view description) {
    BoundaryThermophysicalGhostCertificate foreign_context;
    passed &= expect(
        static_cast<bool>(certify_boundary_thermophysics(context,
                                                         foreign_context)),
        "foreign thermophysical context can be independently issued");
    boundary_input.thermophysical_boundary.certificate = foreign_context;
    rejects_boundary_mutation(description);
    boundary_input.thermophysical_boundary.certificate =
        boundary_c1_thermophysics;
  };
  BoundaryThermophysicalGhostContext foreign_context = boundary_c1_context;
  foreign_context.phase = BoundaryThermophysicalGhostPhase::corrector_two;
  rejects_foreign_context(foreign_context,
                          "C2 physical-ghost certificate cannot enter C1");
  foreign_context = boundary_c1_context;
  ++foreign_context.target_time;
  rejects_foreign_context(
      foreign_context,
      "foreign target-time physical-ghost certificate cannot enter C1");
  foreign_context = boundary_c1_context;
  foreign_context.geometry = foreign_geometry.fingerprint();
  BoundaryThermophysicalGhostCertificate foreign_geometry_certificate;
  passed &= expect(
      certify_boundary_thermophysics(foreign_context,
                                     foreign_geometry_certificate)
                  .code == StatusCode::invalid_plan &&
          !foreign_geometry_certificate.valid(),
      "boundary source rejects a foreign geometry before certificate issue");
  foreign_context = boundary_c1_context;
  ++foreign_context.pressure_reference;
  rejects_foreign_context(
      foreign_context,
      "foreign p-ref physical-ghost certificate cannot enter C1");

  FieldRegistry registry;
  FieldSchema schema;
  FieldId dependency = 0U;
  passed &= expect(
      static_cast<bool>(registry.declare_field("open.final", 1U, 0U,
                                               dependency)) &&
          static_cast<bool>(registry.freeze(schema)),
      "open-boundary transaction schema freezes");
  const std::array requests{ArenaFieldRequest{
      dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage final_storage;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
  passed &= expect(
      static_cast<bool>(ArenaLayout::compile(
          schema, {requests.data(), requests.size()}, layout)) &&
          static_cast<bool>(StateLayers::allocate(layout, layers)) &&
          static_cast<bool>(AttemptTransaction::create(
              layers.field_count(), 1U, layers.field_count(), transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(authority.claim(piso.pressure_stage(),
                                            piso.final_flux_slot(),
                                            transaction, writer)) &&
          static_cast<bool>(transaction.begin(layers)) &&
          static_cast<bool>(transaction.revise_trial(dependency)),
      "open-boundary pending transaction begins");
  PendingFaceFluxView pending;
  passed &= expect(static_cast<bool>(writer.begin_pending(
                       transaction, final_storage, pending)),
                   "open-boundary pending final flux acquires");
  fill(density, 1.0);
  fill(velocity, 0.0);
  const std::array<ConstFieldView, 2U> flux_reads{
      as_const(density.view), as_const(velocity.view)};
  const KernelInvocation flux_invocation{
      {flux_reads.data(), flux_reads.size()}, {}, {{0, 0, 0}, cells},
      0U, 0U, 1U, 0U, nullptr};
  passed &= expect(static_cast<bool>(reconstruct_mass_flux(
                       fixture.equations.kernels(), flux_invocation,
                       pending)),
                   "open-boundary zero pending flux reconstructs for audit");

  OwnedField accepted = make_field(91U, cells, 1U, 0U, 3007U, 4009U);
  OwnedField previous = make_field(92U, cells, 1U, 0U, 3009U, 4010U);
  OwnedField drho_dp = make_field(6U, cells, 1U, 0U, 3010U, 4011U);
  OwnedField pressure_diagonal =
      make_field(53U, cells, 1U, 0U, 3011U, 4012U);
  OwnedField pressure_rhs =
      make_field(54U, cells, 1U, 0U, 3012U, 4013U);
  OwnedField pressure_correction =
      make_field(90U, cells, 1U, 1U, 3013U, 4014U);
  fill(accepted, 1.0);
  fill(previous, 1.0);
  fill(drho_dp, 0.02);
  fill(pressure_correction, 0.0);
  boundary_pressure.view.revision = 3008U;
  fill(boundary_pressure, 9.0);
  PressureCorrectionInput terminal_pressure_input;
  terminal_pressure_input.intermediate = boundary_c2_certificate;
  terminal_pressure_input.pressure_reference =
      boundary_c2_input.pressure_reference;
  terminal_pressure_input.density_trial = as_const(density.view);
  terminal_pressure_input.density_accepted = as_const(accepted.view);
  terminal_pressure_input.density_previous = as_const(previous.view);
  terminal_pressure_input.drho_dp_h_y = as_const(drho_dp.view);
  terminal_pressure_input.bdf = boundary_bdf;
  terminal_pressure_input.time = boundary_c2_input.momentum.time;
  terminal_pressure_input.geometry = boundary_c2_input.momentum.geometry;
  terminal_pressure_input.numeric_boundary =
      boundary_c2_input.numeric_boundary;
  const PressureCorrectionSystemView terminal_pressure_system{
      pressure_diagonal.view, pressure_rhs.view};
  PressureCorrectionCertificate terminal_pressure;
  passed &= expect(
      static_cast<bool>(coupler.assemble_pressure_system(
          terminal_pressure_input, terminal_pressure_system,
          terminal_pressure)) &&
          terminal_pressure.corrector == 2U,
      "open-boundary C2 pressure system assembles for real affine correction");
  PisoFrozenMomentumStageAuthority open_frozen;
  passed &= expect(
      static_cast<bool>(coupler.make_frozen_momentum_stage_authority(
          boundary_c2_certificate, terminal_pressure, open_frozen)) &&
          open_frozen.valid() &&
          open_frozen.scope() ==
              PisoFrozenMomentumStageScope::cartesian_open_boundary_ibm,
      "boundary-absolute corrector issues open/IBM frozen-momentum authority");
  OwnedField open_raw_dp = make_field(
      correction_field, cells, 1U, 1U, 3060U, 4060U);
  OwnedField open_scaled_dp = make_field(
      200U, cells, 1U, 1U, 3061U, 4061U);
  OwnedField open_candidate_velocity = make_field(
      201U, cells, 3U, 1U, 3062U, 4062U);
  fill(open_raw_dp, 2.0);
  fill(open_scaled_dp, 77.0);
  const std::array<HaloFieldSpec, 1U> open_candidate_halo_fields{{
      {open_scaled_dp.view.field, 1U, 1U},
  }};
  HaloEngine open_candidate_halo;
  passed &= expect(
      static_cast<bool>(open_candidate_halo.reserve(
          MPI_COMM_SELF, fixture.patch,
          {open_candidate_halo_fields.data(),
           open_candidate_halo_fields.size()},
          fixture.boundary.halo_topology())),
      "open candidate pressure-correction halo reserves independently");
  PisoFrozenMomentumPressureStageCertificate open_pressure_stage;
  passed &= expect(
      static_cast<bool>(coupler.form_frozen_momentum_scaled_pressure(
          open_frozen, as_const(open_raw_dp.view), open_candidate_halo, 1.0,
          open_scaled_dp.view, open_pressure_stage)),
      "open candidate scales the owned pressure direction");
  std::array<FieldView, 1U> open_candidate_halo_views{open_scaled_dp.view};
  HaloTicket open_candidate_ticket;
  Status open_candidate_exchange = open_candidate_halo.begin(
      3063U,
      {open_candidate_halo_views.data(), open_candidate_halo_views.size()},
      open_candidate_ticket);
  if (open_candidate_exchange)
    open_candidate_exchange = open_candidate_halo.finish(
        open_candidate_ticket,
        {open_candidate_halo_views.data(), open_candidate_halo_views.size()});
  open_scaled_dp.view = open_candidate_halo_views[0U];
  PisoFrozenMomentumVelocityStageCertificate open_velocity_stage;
  const Status open_velocity_status =
      open_candidate_exchange
          ? coupler.stage_frozen_momentum_velocity(
                open_frozen, open_pressure_stage, open_candidate_halo,
                open_scaled_dp.view, open_candidate_velocity.view,
                open_velocity_stage)
          : open_candidate_exchange;
  passed &= expect(
      static_cast<bool>(open_velocity_status) && open_velocity_stage.valid() &&
          open_scaled_dp.view.unchecked({-1, 0, 0}, 0U) == 2.0 &&
          open_scaled_dp.view.unchecked({cells.x, 0, 0}, 0U) == -2.0,
      "open velocity stage owns homogeneous correction ghosts before its gradient");
  OwnedField open_candidate_density = make_field(
      202U, cells, 1U, 1U, 3064U, 4064U);
  fill(open_candidate_density, 1.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        ThermoState candidate_thermo;
        const Real3 candidate_u{
            open_candidate_velocity.view.unchecked(cell, 0U),
            open_candidate_velocity.view.unchecked(cell, 1U),
            open_candidate_velocity.view.unchecked(cell, 2U)};
        passed &= expect(
            static_cast<bool>(fixture.thermodynamics.evaluate(
                101325.0 + boundary_pressure.view.unchecked(cell, 0U),
                boundary_enthalpy.view.unchecked(cell, 0U), {},
                candidate_u, candidate_thermo)),
            "open candidate owned EOS state evaluates");
        open_candidate_density.view.unchecked(cell, 0U) =
            candidate_thermo.rho;
        boundary_thermo_aux[0U].view.unchecked(cell, 0U) =
            candidate_thermo.temperature;
      }
  const std::array<BoundaryGhostFieldAuthority, 2U>
      open_candidate_boundary_authorities{{
          make_boundary_ghost_field_authority(
              as_const(boundary_pressure.view)),
          make_boundary_ghost_field_authority(
              as_const(boundary_enthalpy.view)),
      }};
  const BoundaryThermophysicalGhostAuthority
      open_candidate_boundary_authority{
          3164U,
          fixture.boundary.revision(),
          fixture.boundary.local_layout_fingerprint(),
          cells,
          {reach, reach, reach},
          reach,
          {open_candidate_boundary_authorities.data(),
           open_candidate_boundary_authorities.size()}};
  const BoundaryThermophysicalGhostContext open_candidate_boundary_context{
      terminal_pressure.time,
      fixture.geometry.fingerprint(),
      terminal_pressure_input.pressure_reference.pressure_reference,
      fixture.boundary.revision(),
      BoundaryThermophysicalGhostPhase::corrector_two};
  BoundaryThermophysicalGhostCertificate open_candidate_boundary_thermo;
  passed &= expect(
      static_cast<bool>(BoundaryThermophysicalFaceClosure::close(
          fixture.boundary, fixture.thermodynamics, fixture.transport,
          {101325.0, as_const(boundary_pressure.view),
           as_const(boundary_enthalpy.view), {},
           open_candidate_boundary_authority},
          {open_candidate_density.view, boundary_thermo_aux[0U].view,
           boundary_thermo_aux[1U].view, boundary_thermo_aux[2U].view,
           boundary_thermo_aux[3U].view, boundary_thermo_aux[4U].view,
           boundary_thermo_aux[5U].view, boundary_thermo_aux[6U].view},
          open_candidate_boundary_context,
          open_candidate_boundary_thermo)),
      "open candidate physical p/h state receives same-target EOS closure");
  const std::array<HaloFieldSpec, 5U> open_state_halo_fields{{
      {boundary_pressure.view.field, 1U, 1U},
      {boundary_enthalpy.view.field, 1U, 1U},
      {open_candidate_density.view.field, 1U, 1U},
      {boundary_thermo_aux[0U].view.field, 1U, 1U},
      {open_candidate_velocity.view.field, 1U, 3U},
  }};
  HaloEngine open_state_halo;
  std::array<FieldView, 5U> open_state_halo_views{{
      boundary_pressure.view, boundary_enthalpy.view,
      open_candidate_density.view, boundary_thermo_aux[0U].view,
      open_candidate_velocity.view}};
  HaloTicket open_state_halo_ticket;
  Status open_state_halo_status = open_state_halo.reserve(
      MPI_COMM_SELF, fixture.patch,
      {open_state_halo_fields.data(), open_state_halo_fields.size()},
      fixture.boundary.halo_topology());
  if (open_state_halo_status)
    open_state_halo_status = open_state_halo.begin(
        3165U,
        {open_state_halo_views.data(), open_state_halo_views.size()},
        open_state_halo_ticket);
  if (open_state_halo_status)
    open_state_halo_status = open_state_halo.finish(
        open_state_halo_ticket,
        {open_state_halo_views.data(), open_state_halo_views.size()});
  boundary_pressure.view = open_state_halo_views[0U];
  boundary_enthalpy.view = open_state_halo_views[1U];
  open_candidate_density.view = open_state_halo_views[2U];
  boundary_thermo_aux[0U].view = open_state_halo_views[3U];
  open_candidate_velocity.view = open_state_halo_views[4U];
  passed &= expect(static_cast<bool>(open_state_halo_status),
                   "open candidate state halo completes before finalization");
  FaceFluxStorage open_mechanical_flux_storage;
  FaceFluxView open_mechanical_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, open_mechanical_flux_storage)) &&
          static_cast<bool>(open_mechanical_flux_storage.workspace_view(
              0U, 3065U, open_mechanical_flux)),
      "open candidate mechanical-flux scratch allocates");
  fill_face_flux(open_mechanical_flux, 123.0);
  PisoFrozenMomentumFluxStageCertificate open_flux_stage;
  const Status open_flux_status =
      open_velocity_status
          ? coupler.stage_frozen_momentum_flux(
                open_frozen, open_velocity_stage,
                as_const(open_candidate_density.view), open_mechanical_flux,
                open_flux_stage)
          : open_velocity_status;
  const Int3 outlet_face{cells.x, 0, 0};
  const double candidate_owner_mobility =
      open_candidate_density.view.unchecked(high_owner, 0U) *
      r_au.view.unchecked(high_owner, 0U);
  const double candidate_ghost_mobility =
      open_candidate_density.view.unchecked(outlet_face, 0U) *
      r_au.view.unchecked(outlet_face, 0U);
  const double candidate_outlet_coefficient =
      2.0 * face_area(fixture, CartesianAxis::x, high_owner) /
      (high_distance / candidate_owner_mobility +
       high_distance / candidate_ghost_mobility);
  const double outlet_pressure_response =
      pressure_correction_mass_flux_response(
          as_const(open_scaled_dp.view), fixture.geometry, fixture.patch,
          fixture.boundary, CartesianAxis::x, outlet_face,
          candidate_outlet_coefficient);
  passed &= expect(
      static_cast<bool>(open_flux_status) && open_flux_stage.valid() &&
          open_flux_stage.nonphysical_flux_provenance() != 0U &&
          open_flux_stage.pressure_outlet_provisional_provenance() != 0U &&
          std::abs(open_mechanical_flux.x.unchecked(outlet_face) -
                   (phi.x.unchecked(outlet_face) +
                    outlet_pressure_response)) < 1.0e-12 &&
          open_mechanical_flux.x.unchecked({0, 0, 0}) == 0.0,
      "open mechanical flux separates pressure-outlet provisional authority from finalizer-owned faces");
  FaceFluxStorage open_final_flux_storage;
  FaceFluxView open_final_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, open_final_flux_storage)) &&
          static_cast<bool>(open_final_flux_storage.workspace_view(
              0U, 3166U, open_final_flux)),
      "open final boundary-flux scratch allocates independently");
  fill_face_flux(open_final_flux, -777.0);
  ReductionEngine open_finalizer_reductions;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_SELF, ReductionMode::mpi_allreduce, 8U,
          open_finalizer_reductions)),
      "open boundary finalizer reductions compile");
  PressureEnergyCandidateBoundaryFinalizeInput open_finalize_input;
  open_finalize_input.authority = open_frozen;
  open_finalize_input.pressure_stage = open_pressure_stage;
  open_finalize_input.velocity_stage = open_velocity_stage;
  open_finalize_input.flux_stage = open_flux_stage;
  open_finalize_input.pressure_reference =
      terminal_pressure_input.pressure_reference;
  open_finalize_input.absolute_pressure_reference = 101325.0;
  open_finalize_input.pressure_perturbation =
      as_const(boundary_pressure.view);
  open_finalize_input.enthalpy = as_const(boundary_enthalpy.view);
  open_finalize_input.density = as_const(open_candidate_density.view);
  open_finalize_input.temperature =
      as_const(boundary_thermo_aux[0U].view);
  open_finalize_input.velocity = as_const(open_candidate_velocity.view);
  open_finalize_input.composition_identity = 3167U;
  open_finalize_input.thermophysical_boundary = {
      open_candidate_boundary_thermo,
      {101325.0, as_const(boundary_pressure.view),
       as_const(boundary_enthalpy.view), {},
       as_const(open_candidate_density.view)}};
  open_finalize_input.state_halo = &open_state_halo;
  open_finalize_input.mechanical_flux = as_const(open_mechanical_flux);
  open_finalize_input.final_flux = open_final_flux;
  FinalBoundaryFluxCertificate open_final_boundary;
  const Status open_finalize_status = open_finalizer.finalize(
      open_finalize_input, open_finalizer_reductions,
      open_final_boundary);
  double inlet_enthalpy = 0.0;
  double inlet_cp = 0.0;
  double inlet_gas = 0.0;
  ThermoState inlet_thermo;
  const double inlet_absolute_pressure =
      101325.0 + boundary_pressure.view.unchecked({0, 0, 0}, 0U);
  passed &= expect(
      static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
          300.0, {}, inlet_enthalpy, inlet_cp, inlet_gas)) &&
          static_cast<bool>(fixture.thermodynamics.evaluate(
              inlet_absolute_pressure, inlet_enthalpy, {},
              {1.0, 0.0, 0.0}, inlet_thermo)),
      "independent velocity-inlet EOS oracle evaluates");
  const double expected_inlet_flux =
      inlet_thermo.rho * face_area(fixture, CartesianAxis::x, {0, 0, 0});
  passed &= expect(
      static_cast<bool>(open_finalize_status) &&
          open_final_boundary.valid() &&
          std::abs(open_final_flux.x.unchecked({0, 0, 0}) -
                   expected_inlet_flux) < 1.0e-12 &&
          same_double_bits(open_final_flux.x.unchecked(outlet_face),
                           open_mechanical_flux.x.unchecked(outlet_face)) &&
          same_double_bits(open_final_flux.x.unchecked({1, 0, 0}),
                           open_mechanical_flux.x.unchecked({1, 0, 0})) &&
          same_double_bits(open_final_flux.y.unchecked({0, 0, 0}), 0.0) &&
          same_double_bits(open_final_flux.y.unchecked({0, cells.y, 0}),
                           0.0) &&
          same_double_bits(open_final_flux.z.unchecked({0, 0, 0}), 0.0) &&
          same_double_bits(open_final_flux.z.unchecked({0, 0, cells.z}),
                           0.0),
      "finalizer rebuilds EOS velocity-inlet flux, zeros wall/symmetry faces, and bit-preserves provisional/internal faces");

  const auto rejects_boundary_finalize_without_write = [&] (
      PressureEnergyCandidateBoundaryFinalizeInput mutation,
      std::string_view description) {
    constexpr double sentinel = -991.0;
    fill_face_flux(open_final_flux, sentinel);
    FinalBoundaryFluxCertificate stale = open_final_boundary;
    const Status rejected = open_finalizer.finalize(
        mutation, open_finalizer_reductions, stale);
    passed &= expect(
        rejected.code == StatusCode::invalid_plan && !stale.valid() &&
            face_flux_all_equal(as_const(open_final_flux), sentinel),
        description);
  };
  PressureEnergyCandidateBoundaryFinalizeInput reference_mutation =
      open_finalize_input;
  ++reference_mutation.pressure_reference.pressure_reference;
  rejects_boundary_finalize_without_write(
      reference_mutation,
      "finalizer clears a reused certificate and zero-writes on a foreign pressure reference");

  const double saved_candidate_u =
      open_candidate_velocity.view.unchecked({0, 0, 0}, 0U);
  open_candidate_velocity.view.unchecked({0, 0, 0}, 0U) += 0.25;
  rejects_boundary_finalize_without_write(
      open_finalize_input,
      "finalizer rejects a no-revision candidate velocity mutation with zero write");
  open_candidate_velocity.view.unchecked({0, 0, 0}, 0U) =
      saved_candidate_u;

  const double saved_candidate_density =
      open_candidate_density.view.unchecked({0, 0, 0}, 0U);
  open_candidate_density.view.unchecked({0, 0, 0}, 0U) *= 1.01;
  rejects_boundary_finalize_without_write(
      open_finalize_input,
      "finalizer rejects a no-revision candidate density mutation with zero write");
  open_candidate_density.view.unchecked({0, 0, 0}, 0U) =
      saved_candidate_density;

  const double saved_candidate_temperature =
      boundary_thermo_aux[0U].view.unchecked({0, 0, 0}, 0U);
  boundary_thermo_aux[0U].view.unchecked({0, 0, 0}, 0U) += 10.0;
  rejects_boundary_finalize_without_write(
      open_finalize_input,
      "finalizer rejects a no-revision candidate temperature mutation with zero write");
  boundary_thermo_aux[0U].view.unchecked({0, 0, 0}, 0U) =
      saved_candidate_temperature;

  const Status open_finalize_replay_status = open_finalizer.finalize(
      open_finalize_input, open_finalizer_reductions,
      open_final_boundary);
  passed &= expect(
      static_cast<bool>(open_finalize_replay_status) &&
          open_final_boundary.valid(),
      "failed finalizer probes leave the same candidate replayable");
  ReductionEngine final_boundary_gate_reductions;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_SELF, ReductionMode::mpi_allreduce, 1U,
          final_boundary_gate_reductions)),
      "open final-boundary certificate gate reductions compile");
  PisoFrozenMomentumExactCandidateInput missing_final_boundary_input;
  PisoFrozenMomentumExactCandidateCertificate missing_final_boundary_output;
  passed &= expect(
      coupler.certify_frozen_momentum_exact_baseline(
                 open_frozen, open_pressure_stage, open_velocity_stage,
                 open_flux_stage, missing_final_boundary_input,
                 final_boundary_gate_reductions,
                 missing_final_boundary_output)
                  .code == StatusCode::invalid_plan &&
          !missing_final_boundary_output.valid(),
      "open exact baseline rejects before scratch access when final-boundary authority is missing");
  FinalBoundaryFluxCertificate foreign_final_boundary =
      FinalBoundaryFluxCertificateTestAccess::make_foreign(
          &mutation_coupler, as_const(open_mechanical_flux));
  PisoFrozenMomentumExactCandidateInput foreign_final_boundary_input;
  foreign_final_boundary_input.final_boundary_flux = foreign_final_boundary;
  PisoFrozenMomentumExactCandidateCertificate foreign_final_boundary_output;
  passed &= expect(
      foreign_final_boundary.valid() &&
          coupler.certify_frozen_momentum_exact_baseline(
                     open_frozen, open_pressure_stage, open_velocity_stage,
                     open_flux_stage, foreign_final_boundary_input,
                     final_boundary_gate_reductions,
                     foreign_final_boundary_output)
                  .code == StatusCode::invalid_plan &&
          !foreign_final_boundary_output.valid(),
      "open exact baseline collectively rejects a foreign finalizer certificate");

  fill(open_raw_dp, -1000.0);
  ++open_scaled_dp.view.revision;
  ++open_candidate_velocity.view.revision;
  ++open_mechanical_flux.revision;
  PisoFrozenMomentumPressureStageCertificate backflow_pressure_stage;
  Status backflow_status = coupler.form_frozen_momentum_scaled_pressure(
      open_frozen, as_const(open_raw_dp.view), open_candidate_halo, 1.0,
      open_scaled_dp.view, backflow_pressure_stage);
  std::array<FieldView, 1U> backflow_halo_views{open_scaled_dp.view};
  HaloTicket backflow_halo_ticket;
  if (backflow_status)
    backflow_status = open_candidate_halo.begin(
        3168U,
        {backflow_halo_views.data(), backflow_halo_views.size()},
        backflow_halo_ticket);
  if (backflow_status)
    backflow_status = open_candidate_halo.finish(
        backflow_halo_ticket,
        {backflow_halo_views.data(), backflow_halo_views.size()});
  open_scaled_dp.view = backflow_halo_views[0U];
  PisoFrozenMomentumVelocityStageCertificate backflow_velocity_stage;
  if (backflow_status)
    backflow_status = coupler.stage_frozen_momentum_velocity(
        open_frozen, backflow_pressure_stage, open_candidate_halo,
        open_scaled_dp.view, open_candidate_velocity.view,
        backflow_velocity_stage);
  PisoFrozenMomentumFluxStageCertificate backflow_flux_stage;
  if (backflow_status)
    backflow_status = coupler.stage_frozen_momentum_flux(
        open_frozen, backflow_velocity_stage,
        as_const(open_candidate_density.view), open_mechanical_flux,
        backflow_flux_stage);
  std::array<FieldView, 5U> backflow_state_views{{
      boundary_pressure.view, boundary_enthalpy.view,
      open_candidate_density.view, boundary_thermo_aux[0U].view,
      open_candidate_velocity.view}};
  HaloTicket backflow_state_ticket;
  if (backflow_status)
    backflow_status = open_state_halo.begin(
        3169U, {backflow_state_views.data(), backflow_state_views.size()},
        backflow_state_ticket);
  if (backflow_status)
    backflow_status = open_state_halo.finish(
        backflow_state_ticket,
        {backflow_state_views.data(), backflow_state_views.size()});
  boundary_pressure.view = backflow_state_views[0U];
  boundary_enthalpy.view = backflow_state_views[1U];
  open_candidate_density.view = backflow_state_views[2U];
  boundary_thermo_aux[0U].view = backflow_state_views[3U];
  open_candidate_velocity.view = backflow_state_views[4U];
  constexpr double rejected_backflow_sentinel = -992.0;
  fill_face_flux(open_final_flux, rejected_backflow_sentinel);
  PressureEnergyCandidateBoundaryFinalizeInput rejected_backflow_input =
      open_finalize_input;
  rejected_backflow_input.pressure_stage = backflow_pressure_stage;
  rejected_backflow_input.velocity_stage = backflow_velocity_stage;
  rejected_backflow_input.flux_stage = backflow_flux_stage;
  rejected_backflow_input.velocity =
      as_const(open_candidate_velocity.view);
  rejected_backflow_input.mechanical_flux =
      as_const(open_mechanical_flux);
  FinalBoundaryFluxCertificate rejected_backflow_certificate =
      open_final_boundary;
  const Status rejected_backflow_status =
      backflow_status
          ? open_finalizer.finalize(
                rejected_backflow_input, open_finalizer_reductions,
                rejected_backflow_certificate)
          : backflow_status;
  passed &= expect(
      static_cast<bool>(backflow_status) &&
          open_mechanical_flux.x.unchecked(outlet_face) < 0.0 &&
          rejected_backflow_status.code == StatusCode::rejected_step &&
          !rejected_backflow_certificate.valid() &&
          face_flux_all_equal(as_const(open_final_flux),
                              rejected_backflow_sentinel),
      "first inward pressure-outlet face rejects collectively when backflow is disabled and zero-writes final flux");
  fill(open_raw_dp, 2.0);
  PisoFrozenMomentumStageAuthority post_backflow_replay;
  passed &= expect(
      static_cast<bool>(coupler.make_frozen_momentum_stage_authority(
          boundary_c2_certificate, terminal_pressure,
          post_backflow_replay)) &&
          post_backflow_replay.valid(),
      "rejected pressure-outlet backflow consumes no current frozen authority");
  ReductionEngine correction_reductions;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_SELF, ReductionMode::mpi_allreduce, 1U,
          correction_reductions)),
      "open-boundary affine correction reductions compile");
  PisoStateCorrectionCertificate correction;
  passed &= expect(
      static_cast<bool>(coupler.correct_pending_state(
          terminal_pressure, pressure_correction.view,
          {velocity.view, boundary_pressure.view, density.view,
           as_const(drho_dp.view)},
          pending, correction_reductions, correction)) &&
          correction.valid() && correction.corrector == 2U &&
          correction.closure == PisoStateClosure::pressure_affine &&
          correction.input_pressure_reference.pressure_reference ==
              correction.output_pressure_reference.pressure_reference &&
          correction.closed_gauge_collective_transaction == 0U &&
          correction.closed_gauge_rank_local_transaction == 0U,
      "open-boundary affine C2 publishes input-equals-output reference and zero gauge authority");
  PressureReferenceCertificate pressure_reference =
      boundary_input.pressure_reference;
  BoundaryThermophysicalGhostCertificate boundary_terminal_thermophysics;
  const BoundaryThermophysicalGhostContext boundary_terminal_context{
      pressure_reference.time, fixture.geometry.fingerprint(),
      pressure_reference.pressure_reference, fixture.boundary.revision(),
      BoundaryThermophysicalGhostPhase::terminal};
  passed &= expect(
      static_cast<bool>(certify_boundary_thermophysics(
          boundary_terminal_context, boundary_terminal_thermophysics)),
      "open-boundary terminal physical ghosts are recertified");
  PisoTerminalAuditInput audit;
  audit.correction = correction;
  audit.pressure_reference = pressure_reference;
  audit.density = as_const(density.view);
  audit.eos_density = as_const(density.view);
  audit.density_accepted = as_const(accepted.view);
  audit.density_previous = as_const(previous.view);
  audit.pressure_perturbation = as_const(boundary_pressure.view);
  audit.bdf = boundary_bdf;
  audit.step_dt = time_step_for_bdf(audit.bdf);
  audit.convective_cfl_limit = 1.0e6;
  audit.closed_mass_target = 0.0;
  audit.boundary_closure_residual = 0.0;
  audit.boundary_closure_samples = 1U;
  audit.thermophysical_boundary = {
      boundary_terminal_thermophysics,
      {101325.0, as_const(boundary_pressure.view),
       as_const(boundary_enthalpy.view), {}, as_const(density.view)}};
  ReductionEngine reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_SELF, ReductionMode::mpi_allreduce, 7U,
                       reductions)),
                   "open-boundary audit reductions compile");
  PisoTerminalAuditInput arbitrary_c2_authority = audit;
  arbitrary_c2_authority.correction.thermophysical_boundary_target = 3016U;
  arbitrary_c2_authority.correction
      .thermophysical_boundary_rank_local_binding = 3017U;
  PisoAttemptReport arbitrary_c2_report;
  PisoTerminalCertificate arbitrary_c2_terminal;
  passed &= expect(
      coupler.audit_pending_final(arbitrary_c2_authority, pending, reductions,
                                  arbitrary_c2_report,
                                  arbitrary_c2_terminal)
                  .code == StatusCode::invalid_plan &&
          !arbitrary_c2_terminal.valid(),
      "terminal rejects arbitrary C2 phase target and local binding");
  BoundaryThermophysicalGhostContext wrong_terminal_context =
      boundary_terminal_context;
  wrong_terminal_context.phase =
      BoundaryThermophysicalGhostPhase::corrector_two;
  BoundaryThermophysicalGhostCertificate wrong_terminal_thermophysics;
  passed &= expect(
      static_cast<bool>(certify_boundary_thermophysics(
          wrong_terminal_context, wrong_terminal_thermophysics)),
      "C2 physical ghosts can be recertified for terminal phase test");
  audit.thermophysical_boundary.certificate =
      wrong_terminal_thermophysics;
  PisoAttemptReport wrong_phase_report;
  PisoTerminalCertificate wrong_phase_terminal;
  passed &= expect(
      coupler.audit_pending_final(audit, pending, reductions,
                                  wrong_phase_report,
                                  wrong_phase_terminal)
                  .code == StatusCode::invalid_plan &&
          !wrong_phase_terminal.valid(),
      "C2 physical-ghost certificate cannot enter terminal audit");
  audit.thermophysical_boundary.certificate =
      boundary_terminal_thermophysics;
  PisoTerminalAuditInput invalid_cfl_audit = audit;
  invalid_cfl_audit.step_dt *= 2.0;
  PisoAttemptReport invalid_cfl_report;
  PisoTerminalCertificate invalid_cfl_terminal;
  passed &= expect(
      coupler.audit_pending_final(invalid_cfl_audit, pending, reductions,
                                  invalid_cfl_report,
                                  invalid_cfl_terminal)
                  .code == StatusCode::invalid_plan &&
          invalid_cfl_report.final_flux_revision == 0U &&
          !invalid_cfl_terminal.valid(),
      "terminal CFL audit rejects a positive time step inconsistent with "
      "the frozen BDF target and remains unavailable");
  PisoAttemptReport report;
  PisoTerminalCertificate terminal;
  passed &= expect(
      static_cast<bool>(coupler.audit_pending_final(
          audit, pending, reductions, report, terminal)) &&
          terminal.valid() && report.eos_residual == 0.0 &&
          report.continuity_residual == 0.0 &&
          report.energy_residual == 0.0 &&
          report.closed_mass_residual == 0.0 &&
          report.gauge_residual == 0.0 &&
          report.committed_convective_cfl_out_max == 0.0 &&
          report.committed_convective_cfl_abs_max == 0.0 &&
          report.committed_convective_cfl_limit == 1.0e6 &&
          !report.continuity_witness.valid,
      "open-boundary audit uses boundary closure, ignores mass target, and "
      "does not build a success-path continuity witness");
  const Int3 witness_cell{0, 0, 0};
  const double accepted_density = accepted.view.unchecked(witness_cell, 0U);
  accepted.view.unchecked(witness_cell, 0U) = accepted_density + 0.25;
  PisoAttemptReport continuity_failure_report;
  PisoTerminalCertificate continuity_failure_terminal;
  const Status continuity_failure = coupler.audit_pending_final(
      audit, pending, reductions, continuity_failure_report,
      continuity_failure_terminal);
  passed &= expect(
      continuity_failure.code == StatusCode::rejected_step &&
          !continuity_failure_terminal.valid() &&
          continuity_failure_report.continuity_residual > 0.0 &&
          continuity_failure_report.continuity_witness.valid &&
          continuity_failure_report.continuity_witness.rank == 0 &&
          continuity_failure_report.continuity_witness.global_cell == 0U &&
          std::isfinite(
              continuity_failure_report.continuity_witness.raw_balance) &&
          std::isfinite(continuity_failure_report.continuity_witness.scale) &&
          continuity_failure_report.continuity_witness.scale > 0.0,
      "continuity failure reconstructs a valid detailed witness only after "
      "the global gate fails");
  accepted.view.unchecked(witness_cell, 0U) = accepted_density;
  audit.boundary_closure_residual = 1.0e-4;
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(audit, pending, reductions, report,
                                  terminal)
                  .code == StatusCode::rejected_step &&
          report.gauge_residual == audit.boundary_closure_residual &&
          !terminal.valid(),
      "open-boundary closure mutation rejects terminal publication");
  audit.boundary_closure_residual = 0.0;
  audit.boundary_closure_samples = 0U;
  terminal = {};
  passed &= expect(
      coupler.audit_pending_final(audit, pending, reductions, report,
                                  terminal)
                  .code == StatusCode::invalid_plan &&
          !terminal.valid(),
      "open-boundary audit requires a global outlet-face sample");
  audit.boundary_closure_samples = 1U;
  audit.pressure_reference.kind = PressureReferenceKind::closed_mass;
  passed &= expect(
      coupler.audit_pending_final(audit, pending, reductions, report,
                                  terminal)
                  .code == StatusCode::invalid_plan,
      "open/closed pressure-reference kind cannot be forged");
  const Status rolled_back = transaction.collective_finish(
      MPI_COMM_SELF, {StatusCode::rejected_step, 3016U});
  passed &= expect(rolled_back.code == StatusCode::rejected_step,
                   "rejected open-boundary audit rolls back pending flux");
  return passed;
}

constexpr double kFullRefreshPressureReference = 101325.0;
constexpr BdfCoefficients kFullRefreshBdf{15.0, -20.0, 5.0, 2U};
constexpr RevisionToken kFullRefreshTime = 9101U;
constexpr RevisionToken kFullRefreshPressureReferenceRevision = 9102U;
constexpr FieldId kFullRefreshCorrectionField = 90U;
constexpr FieldId kFullRefreshEnthalpyDirectionField = 93U;
constexpr FieldId kFullRefreshTemperatureDirectionField = 94U;
constexpr FieldId kFrozenCandidateCorrectionField = 200U;
constexpr FieldId kFrozenCandidateVelocityField = 201U;
constexpr FieldId kFrozenCandidateDensityField = 202U;

double full_refresh_phase(std::int32_t index, std::int32_t extent, bool face) {
  constexpr double two_pi = 6.283185307179586476925286766559;
  const double offset = face ? 0.0 : 0.5;
  return two_pi * (static_cast<double>(index) + offset) /
         static_cast<double>(extent);
}

double full_refresh_temperature(Int3 cell, Int3 cells) {
  return 500.0 + 18.0 * std::sin(full_refresh_phase(cell.x, cells.x, false)) +
         11.0 * std::cos(full_refresh_phase(cell.y, cells.y, false)) +
         7.0 * std::sin(full_refresh_phase(cell.z, cells.z, false));
}

double full_refresh_h_by_a(Int3 cell, Int3 cells, std::uint8_t component) {
  const double x = full_refresh_phase(cell.x, cells.x, false);
  const double y = full_refresh_phase(cell.y, cells.y, false);
  const double z = full_refresh_phase(cell.z, cells.z, false);
  if (component == 0U)
    return 0.62 + 0.09 * std::sin(x) + 0.04 * std::cos(y);
  if (component == 1U)
    return -0.31 + 0.07 * std::cos(y) - 0.03 * std::sin(z);
  return 0.18 + 0.05 * std::sin(z) + 0.02 * std::cos(x);
}

double full_refresh_temporal_face(CartesianAxis axis, Int3 face, Int3 cells) {
  const double x = full_refresh_phase(face.x, cells.x, true);
  const double y = full_refresh_phase(face.y, cells.y, true);
  const double z = full_refresh_phase(face.z, cells.z, true);
  if (axis == CartesianAxis::x)
    return 0.23 + 0.031 * std::sin(x) + 0.017 * std::cos(y);
  if (axis == CartesianAxis::y)
    return -0.14 + 0.029 * std::cos(y) - 0.013 * std::sin(z);
  return 0.09 + 0.023 * std::sin(z) + 0.011 * std::cos(x);
}

double full_refresh_c2_face(CartesianAxis axis, Int3 face, Int3 cells) {
  const double phase = axis == CartesianAxis::x
                           ? full_refresh_phase(face.y, cells.y, true)
                           : (axis == CartesianAxis::y
                                  ? full_refresh_phase(face.z, cells.z, true)
                                  : full_refresh_phase(face.x, cells.x, true));
  return 0.37 + 0.05 * static_cast<double>(static_cast<std::uint8_t>(axis)) +
         0.007 * std::sin(phase);
}

ConstFaceFieldView full_refresh_face(ConstFaceFluxView flux,
                                     CartesianAxis axis) {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

FaceFieldView full_refresh_face(FaceFluxView flux, CartesianAxis axis) {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

ConstFaceFieldView full_refresh_face(const FrozenConvectionFaceField& flux,
                                     CartesianAxis axis) {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

ConstFaceFieldView full_refresh_face(ConstFaceFieldView x, ConstFaceFieldView y,
                                     ConstFaceFieldView z, CartesianAxis axis) {
  return axis == CartesianAxis::x ? x : (axis == CartesianAxis::y ? y : z);
}

Int3 full_refresh_wrap(Int3 cell, Int3 cells) {
  const auto wrap = [](std::int32_t value, std::int32_t extent) {
    const std::int32_t remainder = value % extent;
    return remainder < 0 ? remainder + extent : remainder;
  };
  return {wrap(cell.x, cells.x), wrap(cell.y, cells.y), wrap(cell.z, cells.z)};
}

std::size_t full_refresh_offset(Int3 cell, Int3 cells) {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(cell.z));
}

struct FullRefreshSample {
  OwnedField density;
  OwnedField predictor_velocity;
  OwnedField velocity;
  OwnedField pressure;
  OwnedField enthalpy;
  OwnedField temperature;
  OwnedField viscosity;
  OwnedField conductivity;
  OwnedField heat_capacity;
  OwnedField enthalpy_diffusivity;
  OwnedField pressure_compressibility;
  OwnedField enthalpy_compressibility;
  OwnedField velocity_gradient;
  OwnedField momentum_diagonal;
  OwnedField momentum_rhs;
  OwnedField r_au;
  OwnedField h_by_a;
  OwnedField pressure_gradient;
  OwnedFaceField pressure_x;
  OwnedFaceField pressure_y;
  OwnedFaceField pressure_z;
  FaceFluxStorage piso_flux_storage;
  FaceFluxView piso_flux{};
  FaceFluxView c1_trial_flux{};
  FaceFluxView c2_trial_flux{};
  FaceFluxStorage total_flux_storage;
  FaceFluxView total_flux{};
  OwnedField density_accepted;
  OwnedField density_previous;
  OwnedField pressure_diagonal;
  OwnedField pressure_rhs;
  OwnedField history_density;
  OwnedField history_velocity;
  OwnedField history_pressure;
  OwnedField history_enthalpy;
  OwnedField history_temperature;
  OwnedField energy_diagonal;
  OwnedField energy_rhs;
  OwnedField energy_residual;
  OwnedFaceField energy_x;
  OwnedFaceField energy_y;
  OwnedFaceField energy_z;
  HaloEngine state_halo;
  HaloEngine piso_halo;
  HaloEngine correction_halo;
  PressureVelocityCoupler coupler;
  PisoIntermediateCertificate intermediate{};
  PressureCorrectionCertificate pressure_certificate{};
  ConstFaceFluxView intermediate_flux{};
  PisoIntermediateInput refresh_input{};
};

bool full_refresh_exchange_state(const Fixture& fixture,
                                 FullRefreshSample& sample) {
  const std::uint8_t reach = fixture.equations.kernels().reach();
  const std::array<HaloFieldSpec, 11U> specs{{
      {sample.density.view.field, reach, 1U},
      {sample.velocity.view.field, reach, 3U},
      {sample.pressure.view.field, reach, 1U},
      {sample.enthalpy.view.field, reach, 1U},
      {sample.temperature.view.field, reach, 1U},
      {sample.viscosity.view.field, reach, 1U},
      {sample.pressure_compressibility.view.field, reach, 1U},
      {sample.conductivity.view.field, reach, 1U},
      {sample.heat_capacity.view.field, reach, 1U},
      {sample.enthalpy_diffusivity.view.field, reach, 1U},
      {sample.enthalpy_compressibility.view.field, reach, 1U},
  }};
  if (!sample.state_halo.reserve(MPI_COMM_SELF, fixture.patch,
                                 {specs.data(), specs.size()},
                                 fixture.boundary.halo_topology())) {
    std::cerr << "full-refresh state halo reserve failed\n";
    return false;
  }
  std::array<FieldView, 11U> fields{{
      sample.density.view,
      sample.velocity.view,
      sample.pressure.view,
      sample.enthalpy.view,
      sample.temperature.view,
      sample.viscosity.view,
      sample.pressure_compressibility.view,
      sample.conductivity.view,
      sample.heat_capacity.view,
      sample.enthalpy_diffusivity.view,
      sample.enthalpy_compressibility.view,
  }};
  HaloTicket ticket;
  Status status =
      sample.state_halo.begin(311U, {fields.data(), fields.size()}, ticket);
  const Status begin_status = status;
  if (status)
    status = sample.state_halo.finish(ticket, {fields.data(), fields.size()});
  if (!status) {
    std::cerr << "full-refresh state halo begin/status="
              << static_cast<unsigned>(begin_status.code) << '/'
              << begin_status.detail << " -> "
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << '\n';
    return false;
  }
  sample.density.view = fields[0U];
  sample.velocity.view = fields[1U];
  sample.pressure.view = fields[2U];
  sample.enthalpy.view = fields[3U];
  sample.temperature.view = fields[4U];
  sample.viscosity.view = fields[5U];
  sample.pressure_compressibility.view = fields[6U];
  sample.conductivity.view = fields[7U];
  sample.heat_capacity.view = fields[8U];
  sample.enthalpy_diffusivity.view = fields[9U];
  sample.enthalpy_compressibility.view = fields[10U];

  FieldView pressure = sample.pressure.view;
  FieldView enthalpy = sample.enthalpy.view;
  FieldView velocity = sample.velocity.view;
  const Status p_boundary = apply_boundary_ghosts(
      BoundaryStage::pressure, fixture.boundary, {&pressure, 1U}, {});
  const Status h_boundary = apply_boundary_ghosts(
      BoundaryStage::enthalpy, fixture.boundary, {&enthalpy, 1U}, {});
  const Status u_boundary = apply_boundary_ghosts(
      BoundaryStage::momentum, fixture.boundary, {&velocity, 1U}, {});
  if (!p_boundary || !h_boundary || !u_boundary) {
    std::cerr << "full-refresh boundary status="
              << static_cast<unsigned>(p_boundary.code) << '/'
              << p_boundary.detail << ','
              << static_cast<unsigned>(h_boundary.code) << '/'
              << h_boundary.detail << ','
              << static_cast<unsigned>(u_boundary.code) << '/'
              << u_boundary.detail << '\n';
    return false;
  }
  return true;
}

bool build_full_refresh_sample(const Fixture& fixture, const PisoPlan& piso,
                               const CommittedFluxHistory& history,
                               ConstFaceFluxView accepted_flux,
                               ConstFaceFluxView previous_flux,
                               ConstFieldView dp, ConstFieldView dh,
                               double pressure_scale, double enthalpy_scale,
                               RevisionToken seed, bool use_c2_byte_copy,
                               FullRefreshSample& sample) {
  // Keep the owner explicit: accepted_flux/previous_flux borrow this storage.
  (void)history;
  const auto fail = [&](std::string_view stage) {
    std::cerr << "full-refresh sample seed=" << seed << " failed at " << stage
              << '\n';
    return false;
  };
  const Int3 cells = fixture.patch.cells;
  const std::uint8_t reach = fixture.equations.kernels().reach();
  const auto field = [&](FieldId id, std::uint8_t components,
                         std::uint8_t ghosts, RevisionToken delta) {
    return make_field(id, cells, components, ghosts, seed + delta,
                      seed + 1000U + delta);
  };
  sample.density = field(0U, 1U, reach, 1U);
  sample.predictor_velocity = field(1U, 3U, 0U, 60U);
  sample.velocity = field(1U, 3U, reach, 2U);
  sample.pressure = field(2U, 1U, reach, 3U);
  // The coupled E_h operator intentionally requires a two-layer target-h
  // view even for central2, while the state refresh itself exchanges reach 1.
  sample.enthalpy = field(3U, 1U, 2U, 4U);
  sample.temperature = field(4U, 1U, reach, 5U);
  sample.viscosity = field(5U, 1U, reach, 6U);
  sample.pressure_compressibility = field(6U, 1U, reach, 7U);
  sample.velocity_gradient = field(7U, 9U, 1U, 8U);
  sample.conductivity = field(120U, 1U, reach, 9U);
  sample.heat_capacity = field(121U, 1U, reach, 10U);
  sample.enthalpy_diffusivity = field(122U, 1U, reach, 11U);
  sample.enthalpy_compressibility = field(123U, 1U, reach, 12U);
  sample.momentum_diagonal = field(30U, 3U, 0U, 13U);
  sample.momentum_rhs = field(31U, 3U, 0U, 14U);
  sample.r_au = field(40U, 3U, 1U, 15U);
  sample.h_by_a = field(41U, 3U, reach, 16U);
  sample.pressure_gradient = field(42U, 3U, 0U, 17U);
  sample.pressure_x = make_face(CartesianAxis::x, cells, seed + 1018U);
  sample.pressure_y = make_face(CartesianAxis::y, cells, seed + 1019U);
  sample.pressure_z = make_face(CartesianAxis::z, cells, seed + 1020U);
  sample.density_accepted = field(50U, 1U, 0U, 21U);
  sample.density_previous = field(51U, 1U, 0U, 22U);
  sample.pressure_diagonal = field(53U, 1U, 0U, 23U);
  sample.pressure_rhs = field(54U, 1U, 0U, 24U);
  sample.history_density = field(0U, 1U, reach, 25U);
  sample.history_velocity = field(1U, 3U, reach, 26U);
  sample.history_pressure = field(2U, 1U, reach, 27U);
  sample.history_enthalpy = field(3U, 1U, reach, 28U);
  sample.history_temperature = field(4U, 1U, reach, 29U);
  sample.energy_diagonal = field(130U, 1U, 0U, 30U);
  sample.energy_rhs = field(131U, 1U, 0U, 31U);
  sample.energy_residual = field(132U, 1U, 0U, 32U);
  sample.energy_x = make_face(CartesianAxis::x, cells, seed + 1033U);
  sample.energy_y = make_face(CartesianAxis::y, cells, seed + 1034U);
  sample.energy_z = make_face(CartesianAxis::z, cells, seed + 1035U);
  fill(sample.velocity_gradient, 0.0);

  double reference_h = 0.0;
  double reference_cp = 0.0;
  double reference_gas = 0.0;
  if (!fixture.thermodynamics.mixture_enthalpy(500.0, {}, reference_h,
                                               reference_cp, reference_gas)) {
    return fail("reference enthalpy");
  }
  ThermoState reference_thermo;
  if (!fixture.thermodynamics.evaluate(kFullRefreshPressureReference,
                                       reference_h, {}, {}, reference_thermo)) {
    return fail("reference EOS");
  }
  fill(sample.history_density, reference_thermo.rho);
  fill(sample.history_velocity, 0.0);
  fill(sample.history_pressure, 0.0);
  fill(sample.history_enthalpy, reference_h);
  fill(sample.history_temperature, 500.0);
  fill(sample.density_accepted, reference_thermo.rho);
  fill(sample.density_previous, reference_thermo.rho);

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double pi = pressure_scale * dp.unchecked(cell, 0U);
        double base_h = 0.0;
        double base_cp = 0.0;
        double base_gas = 0.0;
        if (!fixture.thermodynamics.mixture_enthalpy(
                full_refresh_temperature(cell, cells), {}, base_h, base_cp,
                base_gas)) {
          return fail("target enthalpy");
        }
        const double enthalpy =
            base_h + enthalpy_scale * dh.unchecked(cell, 0U);
        sample.pressure.view.unchecked(cell, 0U) = pi;
        sample.enthalpy.view.unchecked(cell, 0U) = enthalpy;
        const double volume = cell_volume(fixture, cell);
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double diagonal = kFullRefreshBdf.a0 * volume *
                                  (1.0 + 0.2 * static_cast<double>(component));
          const double h_by_a = full_refresh_h_by_a(cell, cells, component);
          sample.predictor_velocity.view.unchecked(cell, component) = h_by_a;
          sample.momentum_diagonal.view.unchecked(cell, component) = diagonal;
          sample.momentum_rhs.view.unchecked(cell, component) =
              diagonal * h_by_a;
        }
      }
    }
  }

  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        Real3 gradient{};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          Int3 minus = cell;
          Int3 plus = cell;
          if (component == 0U) {
            --minus.x;
            ++plus.x;
          } else if (component == 1U) {
            --minus.y;
            ++plus.y;
          } else {
            --minus.z;
            ++plus.z;
          }
          minus = full_refresh_wrap(minus, cells);
          plus = full_refresh_wrap(plus, cells);
          const CartesianAxis axis = static_cast<CartesianAxis>(component);
          const double spacing = cell_width(
              fixture, axis,
              component == 0U ? cell.x : (component == 1U ? cell.y : cell.z));
          const double value = (sample.pressure.view.unchecked(plus, 0U) -
                                sample.pressure.view.unchecked(minus, 0U)) /
                               (2.0 * spacing);
          if (component == 0U)
            gradient.x = value;
          else if (component == 1U)
            gradient.y = value;
          else
            gradient.z = value;
        }
        const double volume = cell_volume(fixture, cell);
        Real3 velocity{};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double diagonal =
              sample.momentum_diagonal.view.unchecked(cell, component);
          const double reciprocal = volume / diagonal;
          const double grad = component == 0U
                                  ? gradient.x
                                  : (component == 1U ? gradient.y : gradient.z);
          const double value =
              full_refresh_h_by_a(cell, cells, component) - reciprocal * grad;
          sample.velocity.view.unchecked(cell, component) = value;
          if (component == 0U)
            velocity.x = value;
          else if (component == 1U)
            velocity.y = value;
          else
            velocity.z = value;
        }
        ThermoState thermo;
        if (!fixture.thermodynamics.evaluate_from_reference_pressure(
                kFullRefreshPressureReference,
                sample.pressure.view.unchecked(cell, 0U),
                sample.enthalpy.view.unchecked(cell, 0U), {}, velocity,
                thermo)) {
          return fail("target EOS");
        }
        MolecularTransportState transport;
        if (!fixture.transport.evaluate(thermo.temperature, {}, transport)) {
          return fail("target transport");
        }
        sample.density.view.unchecked(cell, 0U) = thermo.rho;
        sample.temperature.view.unchecked(cell, 0U) = thermo.temperature;
        sample.viscosity.view.unchecked(cell, 0U) = transport.viscosity;
        sample.conductivity.view.unchecked(cell, 0U) = transport.conductivity;
        sample.heat_capacity.view.unchecked(cell, 0U) = thermo.cp;
        sample.enthalpy_diffusivity.view.unchecked(cell, 0U) =
            transport.conductivity / thermo.cp;
        sample.pressure_compressibility.view.unchecked(cell, 0U) =
            thermo.drho_dp_hY;
        sample.enthalpy_compressibility.view.unchecked(cell, 0U) =
            thermo.drho_dh_pY;
      }
    }
  }
  if (!full_refresh_exchange_state(fixture, sample))
    return fail("state halo/boundary");

  if (!FaceFluxStorage::allocate_workspace(cells, 3U,
                                           sample.piso_flux_storage) ||
      !sample.piso_flux_storage.workspace_view(0U, seed + 40U,
                                               sample.piso_flux) ||
      !sample.piso_flux_storage.workspace_view(1U, seed + 41U,
                                               sample.c1_trial_flux) ||
      !sample.piso_flux_storage.workspace_view(2U, seed + 42U,
                                               sample.c2_trial_flux) ||
      !FaceFluxStorage::allocate_workspace(cells, 1U,
                                           sample.total_flux_storage) ||
      !sample.total_flux_storage.workspace_view(0U, seed + 43U,
                                                sample.total_flux)) {
    return fail("face workspaces");
  }
  fill_face_flux(sample.c1_trial_flux, 0.0);
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    FaceFieldView temporal = full_refresh_face(sample.piso_flux, axis);
    FaceFieldView c2 = full_refresh_face(sample.c2_trial_flux, axis);
    for (std::int32_t z = 0; z < temporal.extents.z; ++z)
      for (std::int32_t y = 0; y < temporal.extents.y; ++y)
        for (std::int32_t x = 0; x < temporal.extents.x; ++x) {
          const Int3 face{x, y, z};
          temporal.unchecked(face) =
              full_refresh_temporal_face(axis, face, cells);
          c2.unchecked(face) = full_refresh_c2_face(axis, face, cells);
        }
  }

  const PisoCouplerWorkspace workspace{sample.r_au.view,
                                       sample.h_by_a.view,
                                       sample.pressure_gradient.view,
                                       sample.pressure_x.view,
                                       sample.pressure_y.view,
                                       sample.pressure_z.view,
                                       sample.piso_flux};
  const std::array<HaloFieldSpec, 3U> piso_specs{{
      {sample.density.view.field, 1U, 1U},
      {sample.r_au.view.field, 1U, 3U},
      {sample.h_by_a.view.field, reach, 3U},
  }};
  const std::array<HaloFieldSpec, 1U> correction_specs{{
      {kFullRefreshCorrectionField, 1U, 1U},
  }};
  if (!sample.piso_halo.reserve(MPI_COMM_SELF, fixture.patch,
                                {piso_specs.data(), piso_specs.size()},
                                fixture.boundary.halo_topology()) ||
      !sample.correction_halo.reserve(
          MPI_COMM_SELF, fixture.patch,
          {correction_specs.data(), correction_specs.size()},
          fixture.boundary.halo_topology())) {
    return fail("PISO halos");
  }
  const PisoCouplerServices services{MPI_COMM_SELF,
                                     &fixture.geometry,
                                     fixture.patch,
                                     &fixture.boundary,
                                     &fixture.thermodynamics,
                                     &sample.piso_halo,
                                     312U,
                                     sample.density.view.field,
                                     &sample.correction_halo,
                                     313U,
                                     kFullRefreshCorrectionField};
  if (!PressureVelocityCoupler::bind(piso, fixture.equations, services,
                                     workspace, sample.coupler)) {
    return fail("PISO bind");
  }

  PisoIntermediateInput input;
  input.momentum = {fixture.equations.momentum().fingerprint(),
                    EquationAssemblyScope::momentum_predictor,
                    kFullRefreshTime,
                    fixture.geometry.topology_revision(),
                    sample.c1_trial_flux.revision,
                    seed + 44U,
                    time_step_for_bdf(kFullRefreshBdf)};
  input.predictor.plan =
      fixture.equations.thermophysical_predictor().fingerprint();
  input.predictor.time = kFullRefreshTime;
  input.predictor.geometry = fixture.geometry.topology_revision();
  input.predictor.accepted_face_flux = accepted_flux.revision;
  input.predictor.previous_face_flux = previous_flux.revision;
  input.predictor.committed_face_flux_authority =
      accepted_flux.certificate.authority();
  input.predictor.committed_face_flux_storage =
      accepted_flux.certificate.storage();
  input.predictor.committed_face_flux_revision_domain =
      accepted_flux.certificate.revision_domain();
  input.predictor.predicted_density = sample.density.view.revision;
  input.predictor.predicted_density_storage =
      sample.density.view.storage_identity;
  input.predictor.predicted_density_revision_domain =
      sample.density.view.revision_domain;
  input.predictor.paired_face_flux = sample.c1_trial_flux.revision;
  input.predictor.paired_face_flux_storage =
      sample.c1_trial_flux.x.storage_identity;
  input.predictor.paired_face_flux_revision_domain =
      sample.c1_trial_flux.x.revision_domain;
  input.predictor.state = seed + 45U;
  input.predictor.order = kFullRefreshBdf.order;
  input.pressure_reference = {
      fixture.equations.pressure_reference().fingerprint(),
      fixture.equations.thermophysical_predictor().fingerprint(),
      fixture.thermodynamics.fingerprint(),
      seed + 46U,
      kFullRefreshTime,
      kFullRefreshPressureReferenceRevision,
      PressureReferenceKind::closed_mass};
  input.density = sample.density.view;
  // C1 must see the frozen momentum predictor.  The distinct candidate
  // velocity carries -rAU*grad(pi) only into EOS and the energy residual.
  input.trial_velocity = as_const(sample.predictor_velocity.view);
  input.trial_flux = as_const(sample.c1_trial_flux);
  input.momentum_system = {sample.momentum_diagonal.view,
                           sample.momentum_rhs.view};
  input.bdf = kFullRefreshBdf;
  input.numeric_boundary = fixture.boundary.revision();
  input.corrector = 1U;
  input.temporal_reference = as_const(sample.piso_flux);
  input.committed_face_history = {accepted_flux, previous_flux};
  input.thermophysical_boundary.binding.pressure_reference =
      kFullRefreshPressureReference;
  input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(sample.pressure.view);
  input.thermophysical_boundary.binding.enthalpy =
      as_const(sample.enthalpy.view);
  input.thermophysical_boundary.binding.density = as_const(sample.density.view);
  sample.refresh_input = input;
  if (!sample.coupler.refresh(input, sample.intermediate) ||
      !sample.coupler.inspect_intermediate_flux(sample.intermediate,
                                                sample.intermediate_flux)) {
    return fail("C1 refresh/inspect");
  }

  if (use_c2_byte_copy) {
    PisoIntermediateInput c2 = input;
    c2.corrector = 2U;
    c2.prior_corrector = sample.intermediate.dependency;
    c2.trial_flux = as_const(sample.c2_trial_flux);
    c2.temporal_reference = {};
    c2.committed_face_history = {};
    if (!sample.coupler.refresh(c2, sample.intermediate) ||
        !sample.coupler.inspect_intermediate_flux(sample.intermediate,
                                                  sample.intermediate_flux)) {
      return fail("C2 refresh/inspect");
    }
    return true;
  }

  const PressureCorrectionInput pressure_input{
      sample.intermediate,
      input.pressure_reference,
      as_const(sample.density.view),
      as_const(sample.density_accepted.view),
      as_const(sample.density_previous.view),
      as_const(sample.pressure_compressibility.view),
      kFullRefreshBdf,
      kFullRefreshTime,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision()};
  const PressureCorrectionSystemView pressure_system{
      sample.pressure_diagonal.view, sample.pressure_rhs.view};
  if (!sample.coupler.assemble_pressure_system(pressure_input, pressure_system,
                                               sample.pressure_certificate)) {
    return fail("pressure assembly");
  }

  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const ConstFaceFieldView source =
        full_refresh_face(sample.intermediate_flux, axis);
    const ConstFaceFieldView coefficient = full_refresh_face(
        as_const(sample.pressure_x.view), as_const(sample.pressure_y.view),
        as_const(sample.pressure_z.view), axis);
    FaceFieldView total = full_refresh_face(sample.total_flux, axis);
    for (std::int32_t z = 0; z < total.extents.z; ++z)
      for (std::int32_t y = 0; y < total.extents.y; ++y)
        for (std::int32_t x = 0; x < total.extents.x; ++x) {
          const Int3 face{x, y, z};
          total.unchecked(face) =
              source.unchecked(face) + pressure_correction_mass_flux_response(
                                           as_const(sample.pressure.view),
                                           fixture.geometry, fixture.patch,
                                           fixture.boundary, axis, face,
                                           coefficient.unchecked(face));
        }
  }

  const PrimitiveHistory density_history{as_const(sample.density.view),
                                         as_const(sample.history_density.view),
                                         as_const(sample.history_density.view)};
  const PrimitiveHistory velocity_history{
      as_const(sample.velocity.view), as_const(sample.history_velocity.view),
      as_const(sample.history_velocity.view)};
  const PrimitiveHistory pressure_history{
      as_const(sample.pressure.view), as_const(sample.history_pressure.view),
      as_const(sample.history_pressure.view)};
  const PrimitiveHistory enthalpy_history{
      as_const(sample.enthalpy.view), as_const(sample.history_enthalpy.view),
      as_const(sample.history_enthalpy.view)};
  const PrimitiveHistory temperature_history{
      as_const(sample.temperature.view),
      as_const(sample.history_temperature.view),
      as_const(sample.history_temperature.view)};
  const EquationStateView state{density_history,
                                velocity_history,
                                pressure_history,
                                enthalpy_history,
                                temperature_history,
                                kFullRefreshPressureReference,
                                {},
                                {},
                                kFullRefreshPressureReference,
                                kFullRefreshPressureReference};
  EquationMaterialView material;
  material.molecular_viscosity = as_const(sample.viscosity.view);
  material.effective_viscosity = as_const(sample.viscosity.view);
  material.thermal_conductivity = as_const(sample.conductivity.view);
  material.heat_capacity = as_const(sample.heat_capacity.view);
  material.enthalpy_diffusivity = as_const(sample.enthalpy_diffusivity.view);
  material.pressure_compressibility =
      as_const(sample.pressure_compressibility.view);
  EquationAssemblyContext context;
  context.dt = time_step_for_bdf(kFullRefreshBdf);
  context.bdf = kFullRefreshBdf;
  context.time = kFullRefreshTime;
  context.geometry = fixture.geometry.topology_revision();
  context.boundary = fixture.boundary.revision();
  context.thermo = fixture.thermodynamics.fingerprint();
  context.transport = fixture.transport.fingerprint();
  context.face_flux = sample.total_flux.revision;
  context.contribution_stage = 1U;
  context.scope = EquationAssemblyScope::target_coupled;
  context.mass_flux = as_const(sample.total_flux);
  context.provisional_mass_flux = false;
  const EquationSystemView energy_system{
      sample.energy_diagonal.view, sample.energy_rhs.view,
      sample.energy_residual.view, sample.energy_x.view,
      sample.energy_y.view,        sample.energy_z.view};
  EquationAssemblyCertificate energy_certificate;
  const Status energy_status =
      assemble_enthalpy(fixture.equations.enthalpy(), state, material,
                        as_const(sample.velocity_gradient.view), {}, context,
                        energy_system, energy_certificate);
  if (!energy_status) {
    std::cerr << "full-refresh energy status="
              << static_cast<unsigned>(energy_status.code) << '/'
              << energy_status.detail << '\n';
    return fail("energy assembly");
  }
  return true;
}

struct FrozenExactScratch {
  OwnedField scaled_pressure;
  OwnedField scaled_enthalpy;
  OwnedField pressure;
  OwnedField enthalpy;
  OwnedField density;
  OwnedField temperature;
  OwnedField compressibility;
  OwnedField velocity;
  FaceFluxStorage flux_storage;
  FaceFluxView flux{};
  HaloEngine correction_halo;
  PisoFrozenMomentumPressureStageCertificate pressure_stage{};
  PisoFrozenMomentumVelocityStageCertificate velocity_stage{};
  PisoFrozenMomentumFluxStageCertificate flux_stage{};
  PisoExactThermodynamicCandidateView thermodynamic{};
};

PressureReferenceCertificate full_refresh_pressure_reference(
    const Fixture& fixture, RevisionToken seed) {
  return {fixture.equations.pressure_reference().fingerprint(),
          fixture.equations.thermophysical_predictor().fingerprint(),
          fixture.thermodynamics.fingerprint(),
          seed + 46U,
          kFullRefreshTime,
          kFullRefreshPressureReferenceRevision,
          PressureReferenceKind::closed_mass};
}

bool stage_frozen_exact_scratch(
    const Fixture& fixture, FullRefreshSample& sample,
    const PisoFrozenMomentumStageAuthority& authority,
    ConstFieldView raw_pressure_direction,
    ConstFieldView raw_enthalpy_direction, double alpha, RevisionToken seed,
    FrozenExactScratch& scratch) {
  const Int3 cells = fixture.patch.cells;
  scratch.scaled_pressure = make_field(
      kFrozenCandidateCorrectionField, cells, 1U, 1U, seed + 1U,
      seed + 101U);
  scratch.scaled_enthalpy =
      make_field(207U, cells, 1U, 0U, seed + 2U, seed + 102U);
  scratch.pressure =
      make_field(208U, cells, 1U, 0U, seed + 3U, seed + 103U);
  scratch.enthalpy =
      make_field(203U, cells, 1U, 0U, seed + 4U, seed + 104U);
  scratch.density = make_field(kFrozenCandidateDensityField, cells, 1U, 1U,
                               seed + 5U, seed + 105U);
  scratch.temperature =
      make_field(204U, cells, 1U, 0U, seed + 6U, seed + 106U);
  scratch.compressibility = make_field(
      sample.pressure_compressibility.view.field, cells, 1U, 0U,
      seed + 10U, seed + 110U);
  scratch.velocity = make_field(kFrozenCandidateVelocityField, cells, 3U, 1U,
                                seed + 7U, seed + 107U);
  const std::array<HaloFieldSpec, 1U> correction_spec{{
      {scratch.scaled_pressure.view.field, 1U, 1U}}};
  if (!scratch.correction_halo.reserve(
          MPI_COMM_SELF, fixture.patch,
          {correction_spec.data(), correction_spec.size()},
          fixture.boundary.halo_topology()) ||
      !FaceFluxStorage::allocate_workspace(cells, 1U,
                                           scratch.flux_storage) ||
      !scratch.flux_storage.workspace_view(0U, seed + 8U, scratch.flux)) {
    return false;
  }
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const Int3 cell{x, y, z};
        const unsigned outside =
            static_cast<unsigned>(x < 0 || x >= cells.x) +
            static_cast<unsigned>(y < 0 || y >= cells.y) +
            static_cast<unsigned>(z < 0 || z >= cells.z);
        if (outside <= 1U)
          scratch.density.view.unchecked(cell, 0U) =
              sample.density.view.unchecked(cell, 0U);
      }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double scaled_h =
            alpha == 0.0
                ? 0.0
                : alpha * raw_enthalpy_direction.unchecked(cell, 0U);
        scratch.scaled_enthalpy.view.unchecked(cell, 0U) = scaled_h;
        scratch.enthalpy.view.unchecked(cell, 0U) =
            sample.enthalpy.view.unchecked(cell, 0U) + scaled_h;
        scratch.temperature.view.unchecked(cell, 0U) =
            sample.temperature.view.unchecked(cell, 0U);
      }
  if (!sample.coupler.form_frozen_momentum_scaled_pressure(
          authority, raw_pressure_direction, scratch.correction_halo, alpha,
          scratch.scaled_pressure.view, scratch.pressure_stage)) {
    return false;
  }
  std::array<FieldView, 1U> halo_fields{scratch.scaled_pressure.view};
  HaloTicket ticket;
  Status status = scratch.correction_halo.begin(
      seed + 9U, {halo_fields.data(), halo_fields.size()}, ticket);
  if (status)
    status = scratch.correction_halo.finish(
        ticket, {halo_fields.data(), halo_fields.size()});
  scratch.scaled_pressure.view = halo_fields[0U];
  return status && sample.coupler.stage_frozen_momentum_velocity(
                       authority, scratch.pressure_stage,
                       scratch.correction_halo,
                       as_const(scratch.scaled_pressure.view),
                       scratch.velocity.view, scratch.velocity_stage);
}

bool prepare_frozen_exact_thermodynamics(
    const Fixture& fixture, FullRefreshSample& sample,
    const PisoFrozenMomentumStageAuthority& authority,
    const PressureReferenceCertificate& predecessor,
    double absolute_pressure_reference, PlanFingerprint closure_token,
    FrozenExactScratch& scratch, ReductionEngine& reductions) {
  PisoExactEosClosureIdentity closure;
  closure.thermodynamics = predecessor.thermodynamics;
  closure.pressure_reference = predecessor.pressure_reference;
  closure.composition = exact_composition_identity_for_test(
      predecessor.thermodynamics, {}, fixture.patch.cells);
  closure.pressure_state = make_piso_field_revision_identity(
      as_const(sample.pressure.view));
  closure.pressure_correction = make_piso_field_revision_identity(
      as_const(scratch.scaled_pressure.view));
  closure.enthalpy_state = make_piso_field_revision_identity(
      as_const(sample.enthalpy.view));
  closure.enthalpy_correction = make_piso_field_revision_identity(
      as_const(scratch.scaled_enthalpy.view));
  closure.candidate_enthalpy = make_piso_field_revision_identity(
      as_const(scratch.enthalpy.view));
  closure.candidate_density = make_piso_field_revision_identity(
      as_const(scratch.density.view));
  closure.candidate_temperature = make_piso_field_revision_identity(
      as_const(scratch.temperature.view));
  closure.closure = closure_token;
  scratch.thermodynamic.enthalpy = as_const(scratch.enthalpy.view);
  scratch.thermodynamic.density = as_const(scratch.density.view);
  scratch.thermodynamic.temperature = as_const(scratch.temperature.view);
  scratch.thermodynamic.closure = closure;
  scratch.thermodynamic.pressure_compressibility =
      as_const(scratch.compressibility.view);
  scratch.thermodynamic.independent_species = {};
  const Int3 cells = fixture.patch.cells;
  const bool alpha_zero = scratch.pressure_stage.alpha() == 0.0;
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const unsigned outside =
            static_cast<unsigned>(x < 0 || x >= cells.x) +
            static_cast<unsigned>(y < 0 || y >= cells.y) +
            static_cast<unsigned>(z < 0 || z >= cells.z);
        if (outside > 1U) continue;
        const Int3 wrapped{(x % cells.x + cells.x) % cells.x,
                           (y % cells.y + cells.y) % cells.y,
                           (z % cells.z + cells.z) % cells.z};
        const double candidate_pressure =
            sample.pressure.view.unchecked(wrapped, 0U) +
            scratch.scaled_pressure.view.unchecked({x, y, z}, 0U);
        const double candidate_enthalpy =
            scratch.enthalpy.view.unchecked(wrapped, 0U);
        const Real3 candidate_velocity{
            scratch.velocity.view.unchecked(wrapped, 0U),
            scratch.velocity.view.unchecked(wrapped, 1U),
            scratch.velocity.view.unchecked(wrapped, 2U)};
        ThermoState replayed;
        if (!fixture.thermodynamics.evaluate_from_reference_pressure(
                absolute_pressure_reference, candidate_pressure,
                candidate_enthalpy, {}, candidate_velocity, replayed,
                sample.temperature.view.unchecked(wrapped, 0U)))
          return false;
        if (!alpha_zero)
          scratch.density.view.unchecked({x, y, z}, 0U) = replayed.rho;
        if (outside == 0U) {
          if (!alpha_zero)
            scratch.temperature.view.unchecked(wrapped, 0U) =
                replayed.temperature;
          scratch.compressibility.view.unchecked(wrapped, 0U) =
              replayed.drho_dp_hY;
          scratch.pressure.view.unchecked(wrapped, 0U) = candidate_pressure;
        }
      }
  if (!prepare_closed_gauge(
          fixture.equations.pressure_reference(), sample.pressure_certificate,
          predecessor, absolute_pressure_reference,
          as_const(sample.pressure.view),
          as_const(scratch.scaled_pressure.view),
          scratch.thermodynamic.pressure_compressibility, closure_token,
          reductions, scratch.thermodynamic.closed_gauge)) {
    return false;
  }
  const double shift = scratch.thermodynamic.closed_gauge.shift;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        scratch.pressure.view.unchecked(cell, 0U) =
            sample.pressure.view.unchecked(cell, 0U) +
            scratch.scaled_pressure.view.unchecked(cell, 0U) - shift;
      }
  return static_cast<bool>(sample.coupler.stage_frozen_momentum_flux(
      authority, scratch.velocity_stage,
      as_const(scratch.density.view), scratch.flux, scratch.flux_stage));
}

PisoFrozenMomentumExactCandidateInput frozen_exact_input(
    FullRefreshSample& sample, ConstFieldView raw_enthalpy_direction,
    const FrozenExactScratch& scratch) {
  return {raw_enthalpy_direction,
          as_const(scratch.scaled_pressure.view),
          as_const(scratch.scaled_enthalpy.view),
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          as_const(scratch.pressure.view),
          scratch.thermodynamic,
          as_const(scratch.velocity.view),
          as_const(scratch.flux)};
}

struct RefinementAuthoritySample {
  OwnedField raw_pressure_direction;
  OwnedField raw_enthalpy_direction;
  FullRefreshSample state;
  FrozenExactScratch baseline;
  FrozenExactScratch selected;
  FaceFluxStorage provisional_flux_storage;
  FaceFluxView provisional_flux{};
  PisoStateCorrectionCertificate correction{};
  PisoPressureEnergyRefinementStateCertificate authority{};
  PisoIntermediateInput refresh_input{};
};

bool select_refinement_exact_candidate(
    const PisoFrozenMomentumExactCandidateCertificate& exact_baseline,
    const PisoFrozenMomentumExactCandidateCertificate& exact_selected,
    RevisionToken seed,
    PressureEnergyGlobalizationSelectionCertificate& selection) {
  PressureEnergyGlobalizationSample baseline_sample;
  baseline_sample.alpha = 0.0;
  baseline_sample.global_normalized_continuity = 1.0;
  baseline_sample.global_normalized_energy = 1.0;
  baseline_sample.thermodynamically_admissible = true;
  baseline_sample.state_and_flux_finite = true;
  baseline_sample.corrector = exact_baseline.corrector();
  baseline_sample.target_time = exact_baseline.target_time();
  baseline_sample.correction_direction = exact_baseline.correction_direction();
  baseline_sample.state_provenance =
      exact_baseline.candidate_state_provenance();
  baseline_sample.mass_flux_provenance =
      exact_baseline.candidate_mass_flux_provenance();
  std::array<PressureEnergyGlobalizationSample,
             kPressureEnergyGlobalizationCandidateCount>
      candidates{};
  for (std::size_t index = 0U; index < candidates.size(); ++index) {
    candidates[index] = baseline_sample;
    candidates[index].alpha = std::ldexp(1.0, -static_cast<int>(index));
    candidates[index].global_normalized_continuity = 2.0;
    candidates[index].global_normalized_energy = 2.0;
    candidates[index].state_provenance = seed + 100U + index;
    candidates[index].mass_flux_provenance = seed + 200U + index;
  }
  candidates[1U].global_normalized_continuity = 0.25;
  candidates[1U].global_normalized_energy = 0.25;
  candidates[1U].state_provenance =
      exact_selected.candidate_state_provenance();
  candidates[1U].mass_flux_provenance =
      exact_selected.candidate_mass_flux_provenance();
  return select_pressure_energy_globalization(
             baseline_sample, {candidates.data(), candidates.size()},
             selection) &&
         selection.valid() && selection.alpha == 0.5;
}

bool build_refinement_authority_sample(
    const Fixture& fixture, const PisoPlan& piso,
    const CommittedFluxHistory& history, ConstFaceFluxView accepted_flux,
    ConstFaceFluxView previous_flux, RevisionToken seed,
    ReductionEngine& reductions, RefinementAuthoritySample& sample) {
  const auto fail = [&](std::string_view stage) {
    std::cerr << "refinement authority sample seed=" << seed << " failed at "
              << stage << '\n';
    return false;
  };
  const Int3 cells = fixture.patch.cells;
  sample.raw_pressure_direction =
      make_field(kFullRefreshCorrectionField, cells, 1U, 1U, seed + 1U,
                 seed + 101U);
  sample.raw_enthalpy_direction =
      make_field(kFullRefreshEnthalpyDirectionField, cells, 1U, 2U,
                 seed + 2U, seed + 102U);
  fill(sample.raw_pressure_direction, 0.0);
  fill(sample.raw_enthalpy_direction, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        sample.raw_pressure_direction.view.unchecked(cell, 0U) =
            180.0 *
            (std::sin(full_refresh_phase(x, cells.x, false)) -
             0.2 * std::cos(full_refresh_phase(y, cells.y, false)));
        sample.raw_enthalpy_direction.view.unchecked(cell, 0U) =
            2.5 * std::cos(full_refresh_phase(z, cells.z, false)) +
            0.3 * std::sin(full_refresh_phase(x, cells.x, false));
      }

  const RevisionToken state_seed = seed + 200U;
  if (!build_full_refresh_sample(
          fixture, piso, history, accepted_flux, previous_flux,
          as_const(sample.raw_pressure_direction.view),
          as_const(sample.raw_enthalpy_direction.view), 0.0, 0.0, state_seed,
          false, sample.state)) {
    return fail("ordinary C1 refresh");
  }
  const PressureReferenceCertificate initial_reference =
      full_refresh_pressure_reference(fixture, state_seed);
  PisoFrozenMomentumStageAuthority c1_stage;
  if (!sample.state.coupler.make_frozen_momentum_stage_authority(
          sample.state.intermediate, sample.state.pressure_certificate,
          c1_stage) ||
      c1_stage.corrector() != 1U) {
    return fail("C1 frozen-momentum authority");
  }
  FrozenExactScratch c1_baseline;
  FrozenExactScratch c1_selected;
  if (!stage_frozen_exact_scratch(
          fixture, sample.state, c1_stage,
          as_const(sample.raw_pressure_direction.view),
          as_const(sample.raw_enthalpy_direction.view), 0.0, seed + 300U,
          c1_baseline) ||
      !prepare_frozen_exact_thermodynamics(
          fixture, sample.state, c1_stage, initial_reference,
          kFullRefreshPressureReference, seed + 400U, c1_baseline,
          reductions)) {
    return fail("exact C1 baseline staging");
  }
  PisoFrozenMomentumExactCandidateCertificate c1_exact_baseline;
  if (!sample.state.coupler.certify_frozen_momentum_exact_baseline(
          c1_stage, c1_baseline.pressure_stage, c1_baseline.velocity_stage,
          c1_baseline.flux_stage,
          frozen_exact_input(sample.state,
                             as_const(sample.raw_enthalpy_direction.view),
                             c1_baseline),
          reductions, c1_exact_baseline)) {
    return fail("exact C1 baseline certification");
  }
  if (!stage_frozen_exact_scratch(
          fixture, sample.state, c1_stage,
          as_const(sample.raw_pressure_direction.view),
          as_const(sample.raw_enthalpy_direction.view), 0.5, seed + 500U,
          c1_selected) ||
      !prepare_frozen_exact_thermodynamics(
          fixture, sample.state, c1_stage, initial_reference,
          kFullRefreshPressureReference, seed + 600U, c1_selected,
          reductions)) {
    return fail("selected exact C1 staging");
  }
  PisoFrozenMomentumExactCandidateCertificate c1_exact_selected;
  if (!sample.state.coupler.certify_frozen_momentum_exact_candidate(
          c1_stage, c1_exact_baseline, c1_selected.pressure_stage,
          c1_selected.velocity_stage, c1_selected.flux_stage,
          frozen_exact_input(sample.state,
                             as_const(sample.raw_enthalpy_direction.view),
                             c1_selected),
          reductions, c1_exact_selected)) {
    return fail("selected exact C1 certification");
  }
  PressureEnergyGlobalizationSelectionCertificate c1_selection;
  if (!select_refinement_exact_candidate(c1_exact_baseline, c1_exact_selected,
                                         seed + 700U, c1_selection) ||
      !sample.state.piso_flux_storage.workspace_view(
          1U, seed + 1000U, sample.state.c1_trial_flux)) {
    return fail("selected exact C1 publication preflight");
  }
  PisoStateCorrectionCertificate c1_correction;
  if (!sample.state.coupler.commit_frozen_momentum_coupled_trial_state(
          c1_stage, c1_exact_selected, c1_selection,
          {sample.state.velocity.view, sample.state.pressure.view,
           sample.state.enthalpy.view, sample.state.density.view,
           sample.state.temperature.view},
          sample.state.c1_trial_flux, reductions, c1_correction)) {
    return fail("selected exact C1 commit");
  }

  PisoIntermediateInput c2_input = sample.state.refresh_input;
  c2_input.corrector = 2U;
  c2_input.prior_corrector = c1_correction.state;
  c2_input.temporal_reference = {};
  c2_input.committed_face_history = {};
  c2_input.pressure_reference = c1_correction.output_pressure_reference;
  c2_input.density = sample.state.density.view;
  c2_input.trial_velocity = as_const(sample.state.velocity.view);
  c2_input.trial_flux = as_const(sample.state.c1_trial_flux);
  c2_input.thermophysical_boundary = {};
  c2_input.thermophysical_boundary.binding = {
      c1_selected.thermodynamic.closed_gauge.next_pressure_reference,
      as_const(sample.state.pressure.view),
      as_const(sample.state.enthalpy.view),
      {},
      as_const(sample.state.density.view)};
  if (!sample.state.coupler.refresh(c2_input, sample.state.intermediate) ||
      !sample.state.coupler.inspect_intermediate_flux(
          sample.state.intermediate, sample.state.intermediate_flux)) {
    return fail("ordinary C2 refresh");
  }

  const PressureCorrectionInput pressure_input{
      sample.state.intermediate,
      c1_correction.output_pressure_reference,
      as_const(sample.state.density.view),
      as_const(sample.state.density_accepted.view),
      as_const(sample.state.density_previous.view),
      as_const(sample.state.pressure_compressibility.view),
      kFullRefreshBdf,
      kFullRefreshTime,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision()};
  if (!sample.state.coupler.assemble_pressure_system(
          pressure_input,
          {sample.state.pressure_diagonal.view, sample.state.pressure_rhs.view},
          sample.state.pressure_certificate)) {
    return fail("C2 pressure assembly");
  }
  PisoFrozenMomentumStageAuthority stage;
  if (!sample.state.coupler.make_frozen_momentum_stage_authority(
          sample.state.intermediate, sample.state.pressure_certificate,
          stage) ||
      stage.corrector() != 2U) {
    return fail("C2 frozen-momentum authority");
  }

  if (!stage_frozen_exact_scratch(
          fixture, sample.state, stage,
          as_const(sample.raw_pressure_direction.view),
          as_const(sample.raw_enthalpy_direction.view), 0.0, seed + 400U,
          sample.baseline) ||
      !prepare_frozen_exact_thermodynamics(
          fixture, sample.state, stage,
          c1_correction.output_pressure_reference,
          c1_selected.thermodynamic.closed_gauge.next_pressure_reference,
          seed + 1200U, sample.baseline, reductions)) {
    return fail("exact C2 baseline staging");
  }
  PisoFrozenMomentumExactCandidateCertificate exact_baseline;
  if (!sample.state.coupler.certify_frozen_momentum_exact_baseline(
          stage, sample.baseline.pressure_stage,
          sample.baseline.velocity_stage, sample.baseline.flux_stage,
          frozen_exact_input(sample.state,
                             as_const(sample.raw_enthalpy_direction.view),
                             sample.baseline),
          reductions, exact_baseline)) {
    return fail("exact C2 baseline certification");
  }

  if (!stage_frozen_exact_scratch(
          fixture, sample.state, stage,
          as_const(sample.raw_pressure_direction.view),
          as_const(sample.raw_enthalpy_direction.view), 0.5, seed + 600U,
          sample.selected) ||
      !prepare_frozen_exact_thermodynamics(
          fixture, sample.state, stage,
          c1_correction.output_pressure_reference,
          c1_selected.thermodynamic.closed_gauge.next_pressure_reference,
          seed + 1400U, sample.selected, reductions)) {
    return fail("selected exact C2 staging");
  }
  PisoFrozenMomentumExactCandidateCertificate exact_selected;
  if (!sample.state.coupler.certify_frozen_momentum_exact_candidate(
          stage, exact_baseline, sample.selected.pressure_stage,
          sample.selected.velocity_stage, sample.selected.flux_stage,
          frozen_exact_input(sample.state,
                             as_const(sample.raw_enthalpy_direction.view),
                             sample.selected),
          reductions, exact_selected)) {
    return fail("selected exact C2 certification");
  }

  PressureEnergyGlobalizationSelectionCertificate selection;
  if (!select_refinement_exact_candidate(exact_baseline, exact_selected,
                                         seed + 1500U, selection)) {
    return fail("exact C2 selection");
  }

  if (!FaceFluxStorage::allocate_workspace(
      cells, 1U, sample.provisional_flux_storage) ||
      !sample.provisional_flux_storage.workspace_view(
          0U, seed + 1600U, sample.provisional_flux) ||
      !sample.state.coupler
           .commit_frozen_momentum_coupled_refinement_trial_state(
               stage, exact_selected, selection,
               {sample.state.velocity.view, sample.state.pressure.view,
                sample.state.enthalpy.view, sample.state.density.view,
                sample.state.temperature.view},
               sample.provisional_flux, 1U, reductions, sample.correction,
               sample.authority)) {
    return fail("typed provisional C2 commit");
  }

  sample.refresh_input = sample.state.refresh_input;
  sample.refresh_input.corrector = 2U;
  sample.refresh_input.prior_corrector = sample.correction.state;
  sample.refresh_input.temporal_reference = {};
  sample.refresh_input.committed_face_history = {};
  sample.refresh_input.pressure_reference =
      sample.correction.output_pressure_reference;
  sample.refresh_input.density = sample.state.density.view;
  sample.refresh_input.trial_velocity = as_const(sample.state.velocity.view);
  sample.refresh_input.trial_flux = as_const(sample.provisional_flux);
  sample.refresh_input.thermophysical_boundary = {};
  sample.refresh_input.thermophysical_boundary.binding = {
      sample.selected.thermodynamic.closed_gauge.next_pressure_reference,
      as_const(sample.state.pressure.view),
      as_const(sample.state.enthalpy.view),
      {},
      as_const(sample.state.density.view)};
  return sample.correction.valid() && sample.authority.valid();
}

bool test_pressure_energy_refinement_authority() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, false, true),
                       "refinement authority periodic fixture compiles");
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "refinement authority PISO plan compiles");
  const Int3 cells = fixture.patch.cells;
  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(initialize_flux_history(cells, history) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        36001U, history, previous_flux) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        36002U, history, accepted_flux),
                   "refinement authority BDF2 history initializes");
  ReductionEngine reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_SELF, ReductionMode::mpi_allreduce, 2U,
                       reductions)),
                   "refinement authority reductions compile");
  if (!passed) return false;

  RefinementAuthoritySample legal;
  passed &= expect(build_refinement_authority_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       36100U, reductions, legal),
                   "selected exact C2 publishes typed provisional refinement authority");
  if (!passed) return false;

  legal.state.velocity.view.unchecked({-1, 0, 0}, 0U) += 1.0;
  legal.state.pressure.view.unchecked({-1, 0, 0}, 0U) += 2.0;
  legal.state.enthalpy.view.unchecked({-1, 0, 0}, 0U) += 3.0;
  legal.state.density.view.unchecked({-1, 0, 0}, 0U) += 4.0;
  legal.state.temperature.view.unchecked({-1, 0, 0}, 0U) += 5.0;
  PisoIntermediateCertificate refined;
  const Status refined_status =
      legal.state.coupler.refresh_pressure_energy_refinement(
          legal.authority, legal.refresh_input,
          {legal.state.velocity.view, legal.state.pressure.view,
           legal.state.enthalpy.view, legal.state.density.view,
           legal.state.temperature.view},
          refined);
  passed &= expect(
      static_cast<bool>(refined_status) && refined.valid() &&
          refined.corrector == 2U && refined.pressure_energy_refinement == 1U &&
          refined.pressure_energy_refinement_collective_lineage ==
              legal.authority.collective_lineage() &&
          refined.pressure_energy_refinement_lineage == legal.authority.lineage(),
      "ghost-only state mutation preserves typed refinement consumption");

  const PressureCorrectionInput refined_pressure_input{
      refined,
      legal.correction.output_pressure_reference,
      as_const(legal.state.density.view),
      as_const(legal.state.density_accepted.view),
      as_const(legal.state.density_previous.view),
      as_const(legal.state.pressure_compressibility.view),
      kFullRefreshBdf,
      kFullRefreshTime,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision()};
  PressureCorrectionCertificate refined_pressure;
  passed &= expect(
      static_cast<bool>(legal.state.coupler.assemble_pressure_system(
          refined_pressure_input,
          {legal.state.pressure_diagonal.view, legal.state.pressure_rhs.view},
          refined_pressure)) &&
          refined_pressure.valid() &&
          refined_pressure.pressure_energy_refinement == 1U &&
          refined_pressure.pressure_energy_refinement_collective_lineage ==
              legal.authority.collective_lineage() &&
          refined_pressure.pressure_energy_refinement_lineage ==
              legal.authority.lineage(),
      "rank-invariant and rank-local refinement lineages reach the pressure certificate");

  PisoAttemptReport legal_report;
  legal_report.pressure_energy_refinement_solve_calls = 1U;
  legal_report.pressure_energy_refinement_termination =
      PressureEnergyRefinementTermination::component_residuals_converged;
  PisoPressureEnergyRefinementSolveReport& legal_entry =
      legal_report.pressure_energy_refinement[0U];
  legal_entry.target_generation = refined_pressure.time;
  legal_entry.collective_lineage =
      refined_pressure.pressure_energy_refinement_collective_lineage;
  legal_entry.pressure_state = refined_pressure.state;
  legal_entry.linear_identity = {37001U, 37002U, 37003U, 37004U, 37005U};
  legal_entry.ordinal = 1U;
  passed &= expect(
      PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(legal_report),
      "typed refinement report accepts a legal same-target active prefix");

  PisoAttemptReport stale_suffix = legal_report;
  stale_suffix.pressure_energy_refinement[1U] = legal_entry;
  stale_suffix.pressure_energy_refinement[1U].ordinal = 2U;
  ++stale_suffix.pressure_energy_refinement[1U].collective_lineage;
  passed &= expect(
      !PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(stale_suffix),
      "typed refinement report rejects a stale active suffix beyond count");

  PisoAttemptReport two_entry = stale_suffix;
  two_entry.pressure_energy_refinement_solve_calls = 2U;
  passed &= expect(
      PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(two_entry),
      "typed refinement report accepts ordered unique same-target entries");
  PisoAttemptReport wrong_order = two_entry;
  wrong_order.pressure_energy_refinement[1U].ordinal = 3U;
  passed &= expect(
      !PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(wrong_order),
      "typed refinement report rejects a non-consecutive ordinal");
  PisoAttemptReport wrong_target = two_entry;
  ++wrong_target.pressure_energy_refinement[1U].target_generation;
  passed &= expect(
      !PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(wrong_target),
      "typed refinement report rejects mixed target generations");
  PisoAttemptReport duplicate_lineage = two_entry;
  duplicate_lineage.pressure_energy_refinement[1U].collective_lineage =
      duplicate_lineage.pressure_energy_refinement[0U].collective_lineage;
  passed &= expect(
      !PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(
              duplicate_lineage),
      "typed refinement report rejects duplicate collective lineage");
  PisoAttemptReport rejected_termination = legal_report;
  rejected_termination.pressure_energy_refinement_termination =
      PressureEnergyRefinementTermination::rejected_candidate;
  passed &= expect(
      !PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(
              rejected_termination),
      "rejected-candidate termination cannot form a successful final report");
  PisoAttemptReport wrong_capacity = legal_report;
  wrong_capacity.pressure_energy_refinement_termination =
      PressureEnergyRefinementTermination::iteration_capacity_exhausted;
  passed &= expect(
      !PisoPressureSolveEpoch::
          validate_pressure_energy_refinement_report_for_test(wrong_capacity),
      "capacity-exhausted termination requires the full active prefix");

  PisoIntermediateCertificate replayed_refinement;
  passed &= expect(
      legal.state.coupler
                  .refresh_pressure_energy_refinement(
                      legal.authority, legal.refresh_input,
                      {legal.state.velocity.view, legal.state.pressure.view,
                       legal.state.enthalpy.view, legal.state.density.view,
                       legal.state.temperature.view},
                      replayed_refinement)
                  .code == StatusCode::invalid_plan &&
          !replayed_refinement.valid(),
      "typed refinement authority is single-use and rejects replay");

  RefinementAuthoritySample ordinary_repeat;
  passed &= expect(
      build_refinement_authority_sample(
          fixture, piso, history, accepted_flux, previous_flux, 38100U,
          reductions, ordinary_repeat),
      "ordinary-repeat fixture publishes typed provisional refinement authority");
  if (!passed) return false;
  PisoIntermediateCertificate ordinary_repeated_c2;
  passed &= expect(
      ordinary_repeat.state.coupler
                  .refresh(ordinary_repeat.refresh_input, ordinary_repeated_c2)
                  .code == StatusCode::invalid_plan &&
          !ordinary_repeated_c2.valid(),
      "ordinary repeated C2 cannot bypass typed refinement authority");

  RefinementAuthoritySample interior_mutation;
  passed &= expect(
      build_refinement_authority_sample(
          fixture, piso, history, accepted_flux, previous_flux, 40100U,
          reductions, interior_mutation),
      "interior-mutation fixture publishes typed provisional refinement authority");
  if (!passed) return false;
  const Int3 probe{0, 0, 0};
  const double pressure_before =
      interior_mutation.state.pressure.view.unchecked(probe, 0U);
  interior_mutation.state.pressure.view.unchecked(probe, 0U) =
      std::nextafter(pressure_before,
                     std::numeric_limits<double>::infinity());
  PisoIntermediateCertificate interior_rejected;
  passed &= expect(
      interior_mutation.state.coupler
                  .refresh_pressure_energy_refinement(
                      interior_mutation.authority,
                      interior_mutation.refresh_input,
                      {interior_mutation.state.velocity.view,
                       interior_mutation.state.pressure.view,
                       interior_mutation.state.enthalpy.view,
                       interior_mutation.state.density.view,
                       interior_mutation.state.temperature.view},
                      interior_rejected)
                  .code == StatusCode::invalid_plan &&
          !interior_rejected.valid(),
      "unrevisioned interior state-bit mutation invalidates typed refinement authority");

  RefinementAuthoritySample flux_mutation;
  passed &= expect(
      build_refinement_authority_sample(
          fixture, piso, history, accepted_flux, previous_flux, 42100U,
          reductions, flux_mutation),
      "flux-mutation fixture publishes typed provisional refinement authority");
  if (!passed) return false;
  const double flux_before =
      flux_mutation.provisional_flux.x.unchecked(probe);
  flux_mutation.provisional_flux.x.unchecked(probe) =
      std::nextafter(flux_before, std::numeric_limits<double>::infinity());
  PisoIntermediateCertificate flux_rejected;
  passed &= expect(
      flux_mutation.state.coupler
                  .refresh_pressure_energy_refinement(
                      flux_mutation.authority, flux_mutation.refresh_input,
                      {flux_mutation.state.velocity.view,
                       flux_mutation.state.pressure.view,
                       flux_mutation.state.enthalpy.view,
                       flux_mutation.state.density.view,
                       flux_mutation.state.temperature.view},
                      flux_rejected)
                  .code == StatusCode::invalid_plan &&
          !flux_rejected.valid(),
      "unrevisioned provisional-flux bit mutation invalidates typed refinement authority");

  RefinementAuthoritySample state_revision;
  passed &= expect(
      build_refinement_authority_sample(
          fixture, piso, history, accepted_flux, previous_flux, 44100U,
          reductions, state_revision),
      "state-revision fixture publishes typed provisional refinement authority");
  if (!passed) return false;
  PisoCoupledStateView replaced_state{
      state_revision.state.velocity.view, state_revision.state.pressure.view,
      state_revision.state.enthalpy.view, state_revision.state.density.view,
      state_revision.state.temperature.view};
  ++replaced_state.pressure_perturbation.revision;
  PisoIntermediateInput replaced_state_input = state_revision.refresh_input;
  replaced_state_input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(replaced_state.pressure_perturbation);
  PisoIntermediateCertificate state_revision_rejected;
  passed &= expect(
      state_revision.state.coupler
                  .refresh_pressure_energy_refinement(
                      state_revision.authority, replaced_state_input,
                      replaced_state, state_revision_rejected)
                  .code == StatusCode::invalid_plan &&
          !state_revision_rejected.valid(),
      "same-storage state revision replacement invalidates typed refinement authority");

  RefinementAuthoritySample flux_revision;
  passed &= expect(
      build_refinement_authority_sample(
          fixture, piso, history, accepted_flux, previous_flux, 46100U,
          reductions, flux_revision),
      "flux-revision fixture publishes typed provisional refinement authority");
  if (!passed) return false;
  FaceFluxView replaced_flux;
  passed &= expect(
      static_cast<bool>(flux_revision.provisional_flux_storage.workspace_view(
          0U, flux_revision.provisional_flux.revision + 1U, replaced_flux)),
      "same-storage replacement provisional flux view is available");
  if (!passed) return false;
  PisoIntermediateInput replaced_flux_input = flux_revision.refresh_input;
  replaced_flux_input.trial_flux = as_const(replaced_flux);
  PisoIntermediateCertificate flux_revision_rejected;
  passed &= expect(
      flux_revision.state.coupler
                  .refresh_pressure_energy_refinement(
                      flux_revision.authority, replaced_flux_input,
                      {flux_revision.state.velocity.view,
                       flux_revision.state.pressure.view,
                       flux_revision.state.enthalpy.view,
                       flux_revision.state.density.view,
                       flux_revision.state.temperature.view},
                      flux_revision_rejected)
                  .code == StatusCode::invalid_plan &&
          !flux_revision_rejected.valid(),
      "same-storage flux revision replacement invalidates typed refinement authority");
  return passed;
}

bool test_frozen_momentum_candidate_alpha_zero_exact() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, false, true),
                       "candidate-stage periodic fixture compiles");
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "candidate-stage PISO plan compiles");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(initialize_flux_history(cells, history) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        11201U, history, previous_flux) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        11202U, history, accepted_flux),
                   "candidate-stage BDF2 history initializes");
  if (!passed) return false;

  OwnedField raw_dp = make_field(kFullRefreshCorrectionField, cells, 1U, 1U,
                                 11203U, 12203U);
  OwnedField dh = make_field(kFullRefreshEnthalpyDirectionField, cells, 1U, 2U,
                             11204U, 12204U);
  fill(raw_dp, 0.0);
  fill(dh, 0.0);
  FullRefreshSample baseline;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(raw_dp.view), as_const(dh.view), 0.0, 0.0,
                       11300U, false, baseline),
                   "candidate-stage baseline refreshes and assembles pressure");
  if (!passed) return false;

  OwnedField scaled_dp = make_field(kFrozenCandidateCorrectionField, cells, 1U, 1U,
                                    11350U, 12350U);
  OwnedField candidate_velocity =
      make_field(kFrozenCandidateVelocityField, cells, 3U, 0U, 11351U,
                 12351U);
  OwnedField candidate_density =
      make_field(kFrozenCandidateDensityField, cells, 1U, 1U, 11353U,
                 12353U);
  fill(scaled_dp, -17.0);
  fill(candidate_velocity, -19.0);
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const unsigned outside =
            static_cast<unsigned>(x < 0 || x >= cells.x) +
            static_cast<unsigned>(y < 0 || y >= cells.y) +
            static_cast<unsigned>(z < 0 || z >= cells.z);
        if (outside <= 1U)
          candidate_density.view.unchecked({x, y, z}, 0U) =
              baseline.density.view.unchecked({x, y, z}, 0U);
      }
  HaloEngine candidate_correction_halo;
  const std::array<HaloFieldSpec, 1U> candidate_correction_specs{{
      {scaled_dp.view.field, 1U, 1U}}};
  passed &= expect(
      static_cast<bool>(candidate_correction_halo.reserve(
          MPI_COMM_SELF, fixture.patch,
          {candidate_correction_specs.data(), candidate_correction_specs.size()},
          fixture.boundary.halo_topology())) &&
          candidate_correction_halo.instance_identity() !=
              baseline.correction_halo.instance_identity(),
      "candidate correction halo is independently reserved");
  FaceFluxStorage candidate_flux_storage;
  FaceFluxView candidate_flux;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, candidate_flux_storage)) &&
          static_cast<bool>(candidate_flux_storage.workspace_view(
              0U, 11352U, candidate_flux)),
      "candidate-stage total-flux scratch allocates");
  if (!passed) return false;
  fill_face_flux(candidate_flux, -23.0);

  PisoFrozenMomentumStageAuthority stage;
  PisoFrozenMomentumPressureStageCertificate pressure_stage;
  PisoFrozenMomentumVelocityStageCertificate velocity_stage;
  PisoFrozenMomentumFluxStageCertificate flux_stage;
  passed &= expect(
      static_cast<bool>(baseline.coupler.make_frozen_momentum_stage_authority(
          baseline.intermediate, baseline.pressure_certificate, stage)) &&
          stage.valid() &&
          stage.scope() ==
              PisoFrozenMomentumStageScope::cartesian_periodic &&
          stage.corrector() == 1U,
      "current C1 intermediate and pressure issue Cartesian frozen-momentum authority");
  passed &= expect(
      scaled_dp.view.revision == 11350U &&
      static_cast<bool>(baseline.coupler.form_frozen_momentum_scaled_pressure(
          stage, as_const(raw_dp.view), candidate_correction_halo, 0.0,
          scaled_dp.view,
          pressure_stage)) &&
          scaled_dp.view.revision == 11350U &&
          pressure_stage.valid() && pressure_stage.alpha() == 0.0 &&
          pressure_stage.canonical_lineage() != 0U &&
          pressure_stage.scratch_binding() != 0U,
      "alpha-zero pressure stage writes owned scaled correction only");
  PisoFrozenMomentumVelocityStageCertificate stale_ghost_rejected;
  passed &= expect(
      baseline.coupler
              .stage_frozen_momentum_velocity(
                  stage, pressure_stage, candidate_correction_halo,
                  as_const(scaled_dp.view),
                  candidate_velocity.view, stale_ghost_rejected)
              .code == StatusCode::invalid_plan &&
          !stale_ghost_rejected.valid(),
      "velocity stage rejects scaled correction before its candidate halo");
  std::array<FieldView, 1U> correction_halo_fields{scaled_dp.view};
  HaloTicket correction_ticket;
  Status correction_halo_status = candidate_correction_halo.begin(
      314U, {correction_halo_fields.data(), correction_halo_fields.size()},
      correction_ticket);
  if (correction_halo_status)
    correction_halo_status = candidate_correction_halo.finish(
        correction_ticket,
        {correction_halo_fields.data(), correction_halo_fields.size()});
  scaled_dp.view = correction_halo_fields[0U];
  passed &= expect(static_cast<bool>(correction_halo_status),
                   "candidate scaled correction halo completes");
  PisoFrozenMomentumPressureStageCertificate reused_revision_rejected;
  passed &= expect(
      baseline.coupler
                  .form_frozen_momentum_scaled_pressure(
                      stage, as_const(raw_dp.view), candidate_correction_halo,
                      0.0, scaled_dp.view, reused_revision_rejected)
                  .code == StatusCode::invalid_plan &&
          !reused_revision_rejected.valid(),
      "candidate pressure stage preserves caller revision and rejects an already-exchanged write revision");
  passed &= expect(
      static_cast<bool>(baseline.coupler.stage_frozen_momentum_velocity(
          stage, pressure_stage, candidate_correction_halo,
          as_const(scaled_dp.view),
          candidate_velocity.view, velocity_stage)) &&
          velocity_stage.valid() && velocity_stage.alpha() == 0.0 &&
          velocity_stage.canonical_lineage() != 0U &&
          velocity_stage.scratch_binding() != 0U,
      "halo-complete alpha-zero velocity stage issues bound canonical provenance");
  const Status flux_stage_status =
      baseline.coupler.stage_frozen_momentum_flux(
          stage, velocity_stage, as_const(candidate_density.view),
          candidate_flux, flux_stage);
  if (!flux_stage_status)
    std::cerr << "candidate alpha-zero flux status="
              << static_cast<unsigned>(flux_stage_status.code) << '/'
              << flux_stage_status.detail << '\n';
  passed &= expect(
      static_cast<bool>(flux_stage_status) &&
          flux_stage.valid() && flux_stage.alpha() == 0.0 &&
          flux_stage.canonical_lineage() != 0U &&
          flux_stage.scratch_binding() != 0U,
      "alpha-zero density/flux stage issues bound canonical provenance");

  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        passed &= expect(
            scaled_dp.view.unchecked(cell, 0U) == 0.0,
            "alpha-zero staged correction is exactly zero");
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double actual =
              candidate_velocity.view.unchecked(cell, component);
          const double expected =
              baseline.h_by_a.view.unchecked(cell, component);
          passed &= expect(std::memcmp(&actual, &expected, sizeof(double)) == 0,
                           "alpha-zero candidate U is byte-equal HbyA");
        }
      }
  const std::vector<double> actual_flux =
      face_flux_values(as_const(candidate_flux));
  const std::vector<double> expected_flux =
      face_flux_values(baseline.intermediate_flux);
  passed &= expect(actual_flux.size() == expected_flux.size(),
                   "alpha-zero candidate/base flux shapes match");
  for (std::size_t index = 0U;
       index < std::min(actual_flux.size(), expected_flux.size()); ++index) {
    passed &= expect(std::memcmp(&actual_flux[index], &expected_flux[index],
                                 sizeof(double)) == 0,
                     "alpha-zero candidate flux is byte-equal base flux");
  }

  ConstFaceFluxView still_current;
  PisoFrozenMomentumStageAuthority replay;
  passed &= expect(
      static_cast<bool>(baseline.coupler.inspect_intermediate_flux(
          baseline.intermediate, still_current)) &&
          static_cast<bool>(baseline.coupler.make_frozen_momentum_stage_authority(
              baseline.intermediate, baseline.pressure_certificate, replay)) &&
          replay.valid(),
      "candidate staging leaves current intermediate and pressure authority reusable");

  return passed;
}

bool test_frozen_momentum_candidate_c1_c2_mapping() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, false, true),
                       "candidate mapping periodic fixture compiles");
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "candidate mapping PISO plan compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;

  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(initialize_flux_history(cells, history) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        11401U, history, previous_flux) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        11402U, history, accepted_flux),
                   "candidate mapping BDF2 history initializes");
  if (!passed) return false;

  OwnedField raw_dp = make_field(kFullRefreshCorrectionField, cells, 1U, 1U,
                                 11403U, 12403U);
  OwnedField dh = make_field(kFullRefreshEnthalpyDirectionField, cells, 1U, 2U,
                             11404U, 12404U);
  fill(raw_dp, 0.0);
  fill(dh, 0.0);
  FullRefreshSample c1;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(raw_dp.view), as_const(dh.view), 0.0, 0.0,
                       11500U, false, c1),
                   "candidate mapping C1 baseline closes");
  if (!passed) return false;

  PisoFrozenMomentumStageAuthority c1_stage;
  passed &= expect(
      static_cast<bool>(c1.coupler.make_frozen_momentum_stage_authority(
          c1.intermediate, c1.pressure_certificate, c1_stage)),
      "candidate mapping C1 authority issues");
  OwnedField scaled_dp = make_field(kFrozenCandidateCorrectionField, cells,
                                    1U, 1U,
                                    11540U, 12540U);
  OwnedField candidate_velocity =
      make_field(kFrozenCandidateVelocityField, cells, 3U, 1U, 11541U,
                 12541U);
  OwnedField candidate_density =
      make_field(kFrozenCandidateDensityField, cells, 1U, 1U, 11542U,
                 12542U);
  FaceFluxStorage candidate_flux_storage;
  FaceFluxView candidate_flux;
  HaloEngine candidate_halo;
  HaloEngine foreign_halo;
  const std::array<HaloFieldSpec, 1U> correction_spec{{
      {scaled_dp.view.field, 1U, 1U}}};
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          cells, 1U, candidate_flux_storage)) &&
          static_cast<bool>(candidate_flux_storage.workspace_view(
              0U, 11543U, candidate_flux)) &&
          static_cast<bool>(candidate_halo.reserve(
              MPI_COMM_SELF, fixture.patch,
              {correction_spec.data(), correction_spec.size()},
              fixture.boundary.halo_topology())) &&
          static_cast<bool>(foreign_halo.reserve(
              MPI_COMM_SELF, fixture.patch,
              {correction_spec.data(), correction_spec.size()},
              fixture.boundary.halo_topology())),
      "candidate mapping scratch and two independent halos allocate");
  if (!passed) return false;

  const std::vector<double> frozen_density_bytes = c1.density.bytes;
  const std::vector<double> frozen_h_by_a_bytes = c1.h_by_a.bytes;
  const std::vector<double> frozen_base_flux =
      face_flux_values(c1.intermediate_flux);
  std::array<PlanFingerprint, 3U> state_lineages{};
  std::array<PlanFingerprint, 3U> flux_lineages{};

  const auto component_index = [](CartesianAxis axis) {
    return static_cast<std::uint8_t>(axis);
  };
  const auto normal_index = [](Int3 value, CartesianAxis axis) {
    return axis == CartesianAxis::x
               ? value.x
               : (axis == CartesianAxis::y ? value.y : value.z);
  };
  const auto shift = [](Int3 value, CartesianAxis axis, int amount) {
    if (axis == CartesianAxis::x)
      value.x += amount;
    else if (axis == CartesianAxis::y)
      value.y += amount;
    else
      value.z += amount;
    return value;
  };
  const auto density_factor = [&](Int3 cell) {
    const Int3 wrapped = full_refresh_wrap(cell, cells);
    return 1.0 + 0.025 *
                     (std::sin(full_refresh_phase(wrapped.x, cells.x, false)) +
                      0.4 * std::cos(
                                full_refresh_phase(wrapped.y, cells.y, false)) -
                      0.2 * std::sin(
                                full_refresh_phase(wrapped.z, cells.z, false)));
  };
  const auto fill_direction = [&](bool enabled) {
    std::fill(raw_dp.bytes.begin(), raw_dp.bytes.end(),
              std::numeric_limits<double>::quiet_NaN());
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          raw_dp.view.unchecked(cell, 0U) =
              enabled
                  ? 3200.0 *
                        (std::sin(full_refresh_phase(x, cells.x, false)) -
                         0.35 * std::cos(
                                    full_refresh_phase(y, cells.y, false)) +
                         0.15 * std::sin(
                                    full_refresh_phase(z, cells.z, false)))
                  : 0.0;
        }
  };
  const auto fill_candidate_density = [&](bool enabled) {
    for (std::int32_t z = -1; z <= cells.z; ++z)
      for (std::int32_t y = -1; y <= cells.y; ++y)
        for (std::int32_t x = -1; x <= cells.x; ++x) {
          const Int3 cell{x, y, z};
          const unsigned outside =
              static_cast<unsigned>(x < 0 || x >= cells.x) +
              static_cast<unsigned>(y < 0 || y >= cells.y) +
              static_cast<unsigned>(z < 0 || z >= cells.z);
          if (outside > 1U) {
            candidate_density.view.unchecked(cell, 0U) =
                std::numeric_limits<double>::quiet_NaN();
            continue;
          }
          const Int3 wrapped = full_refresh_wrap(cell, cells);
          const double base = c1.density.view.unchecked(wrapped, 0U);
          candidate_density.view.unchecked(cell, 0U) =
              enabled ? base * density_factor(wrapped) : base;
        }
  };
  const auto exchange_scaled = [&](HaloEngine& halo, StageId stage) {
    std::array<FieldView, 1U> views{scaled_dp.view};
    HaloTicket ticket;
    Status status = halo.begin(stage, {views.data(), views.size()}, ticket);
    if (status) status = halo.finish(ticket, {views.data(), views.size()});
    scaled_dp.view = views[0U];
    return status;
  };
  const auto expected_face = [&](CartesianAxis axis, Int3 face,
                                 bool density_changed) {
    const std::uint8_t component = component_index(axis);
    const Int3 left = full_refresh_wrap(shift(face, axis, -1), cells);
    const Int3 right = full_refresh_wrap(face, cells);
    const double base_left_rho = c1.density.view.unchecked(left, 0U);
    const double base_right_rho = c1.density.view.unchecked(right, 0U);
    const double candidate_left_rho =
        density_changed ? base_left_rho * density_factor(left) : base_left_rho;
    const double candidate_right_rho = density_changed
                                           ? base_right_rho * density_factor(right)
                                           : base_right_rho;
    const double left_h = c1.h_by_a.view.unchecked(left, component);
    const double right_h = c1.h_by_a.view.unchecked(right, component);
    const double left_r = c1.r_au.view.unchecked(left, component);
    const double right_r = c1.r_au.view.unchecked(right, component);
    const double area = face_area(fixture, axis, right);
    const double theta = -full_refresh_temporal_face(axis, face, cells);
    const auto q = [&](double rho_left, double rho_right) {
      return area * 0.5 * (rho_left * left_h + rho_right * right_h) +
             kFullRefreshBdf.a0 * 0.5 *
                 (rho_left * left_r + rho_right * right_r) * theta;
    };
    const std::int32_t left_normal = normal_index(left, axis);
    const std::int32_t right_normal = normal_index(right, axis);
    const double left_width = cell_width(fixture, axis, left_normal);
    const double right_width = cell_width(fixture, axis, right_normal);
    const double coefficient =
        area /
        (0.5 * left_width / (candidate_left_rho * left_r) +
         0.5 * right_width / (candidate_right_rho * right_r));
    const double jump = scaled_dp.view.unchecked(right, 0U) -
                        scaled_dp.view.unchecked(left, 0U);
    return full_refresh_face(c1.intermediate_flux, axis).unchecked(face) +
           q(candidate_left_rho, candidate_right_rho) -
           q(base_left_rho, base_right_rho) - coefficient * jump;
  };

  struct MappingCase {
    bool pressure;
    bool density;
    double alpha;
  };
  constexpr std::array<MappingCase, 3U> cases{{
      {true, false, 0.5}, {false, true, 1.0}, {true, true, 0.25}}};
  for (std::size_t case_index = 0U; case_index < cases.size(); ++case_index) {
    const MappingCase selected = cases[case_index];
    if (case_index != 0U) ++scaled_dp.view.revision;
    fill_direction(selected.pressure);
    fill_candidate_density(selected.density);
    fill(scaled_dp, -777.0);
    fill(candidate_velocity, -888.0);
    fill_face_flux(candidate_flux, -999.0);
    PisoFrozenMomentumPressureStageCertificate pressure_stage;
    PisoFrozenMomentumVelocityStageCertificate velocity_stage;
    PisoFrozenMomentumFluxStageCertificate flux_stage;
    passed &= expect(
        static_cast<bool>(c1.coupler.form_frozen_momentum_scaled_pressure(
            c1_stage, as_const(raw_dp.view), candidate_halo, selected.alpha,
            scaled_dp.view, pressure_stage)),
        "C1 candidate forms deterministic owned scaled pressure");
    if (case_index == 0U) {
      PisoFrozenMomentumVelocityStageCertificate skipped;
      passed &= expect(
          c1.coupler
                  .stage_frozen_momentum_velocity(
                      c1_stage, pressure_stage, candidate_halo,
                      as_const(scaled_dp.view), candidate_velocity.view, skipped)
                  .code == StatusCode::invalid_plan &&
              !skipped.valid(),
          "C1 velocity rejects skipped candidate correction halo");
      passed &= expect(static_cast<bool>(exchange_scaled(foreign_halo, 317U)),
                       "foreign candidate halo exchanges the scaled field");
      PisoFrozenMomentumVelocityStageCertificate foreign;
      passed &= expect(
          c1.coupler
                  .stage_frozen_momentum_velocity(
                      c1_stage, pressure_stage, foreign_halo,
                      as_const(scaled_dp.view), candidate_velocity.view, foreign)
                  .code == StatusCode::invalid_plan &&
              !foreign.valid(),
          "pressure-stage authority rejects a foreign compatible halo");
    }
    passed &= expect(static_cast<bool>(exchange_scaled(candidate_halo, 316U)),
                     "bound C1 candidate correction halo exchanges");
    passed &= expect(
        static_cast<bool>(c1.coupler.stage_frozen_momentum_velocity(
            c1_stage, pressure_stage, candidate_halo,
            as_const(scaled_dp.view), candidate_velocity.view,
            velocity_stage)) &&
            static_cast<bool>(c1.coupler.stage_frozen_momentum_flux(
                c1_stage, velocity_stage, as_const(candidate_density.view),
                candidate_flux, flux_stage)),
        "C1 pressure/density candidate stages velocity and total flux");
    if (!velocity_stage.valid() || !flux_stage.valid()) continue;
    state_lineages[case_index] = velocity_stage.canonical_lineage();
    flux_lineages[case_index] = flux_stage.canonical_lineage();

    double maximum_velocity_error = 0.0;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const auto axis = static_cast<CartesianAxis>(component);
            const Int3 minus = full_refresh_wrap(shift(cell, axis, -1), cells);
            const Int3 plus = full_refresh_wrap(shift(cell, axis, 1), cells);
            const double gradient =
                (scaled_dp.view.unchecked(plus, 0U) -
                 scaled_dp.view.unchecked(minus, 0U)) /
                (2.0 * cell_width(fixture, axis, normal_index(cell, axis)));
            const double expected =
                c1.h_by_a.view.unchecked(cell, component) -
                c1.r_au.view.unchecked(cell, component) * gradient;
            maximum_velocity_error = std::max(
                maximum_velocity_error,
                std::abs(candidate_velocity.view.unchecked(cell, component) -
                         expected));
          }
        }
    double maximum_flux_error = 0.0;
    for (CartesianAxis axis : {CartesianAxis::x, CartesianAxis::y,
                               CartesianAxis::z}) {
      const ConstFaceFieldView actual =
          full_refresh_face(as_const(candidate_flux), axis);
      for (std::int32_t z = 0; z < actual.extents.z; ++z)
        for (std::int32_t y = 0; y < actual.extents.y; ++y)
          for (std::int32_t x = 0; x < actual.extents.x; ++x) {
            const Int3 face{x, y, z};
            maximum_flux_error = std::max(
                maximum_flux_error,
                std::abs(actual.unchecked(face) -
                         expected_face(axis, face, selected.density)));
          }
    }
    passed &= expect(maximum_velocity_error < 2.0e-12 &&
                         maximum_flux_error < 2.0e-12,
                     "C1 dp-only/rho-only/mixed mapping matches independent oracle");

    if (case_index == 0U) {
      const PisoFrozenMomentumVelocityStageCertificate stale_velocity =
          velocity_stage;
      ++scaled_dp.view.revision;
      passed &= expect(
          static_cast<bool>(exchange_scaled(candidate_halo, 318U)),
          "candidate correction halo can advance to a newer ghost revision");
      PisoFrozenMomentumFluxStageCertificate stale_halo_rejected;
      passed &= expect(
          c1.coupler
                  .stage_frozen_momentum_flux(
                      c1_stage, stale_velocity,
                      as_const(candidate_density.view), candidate_flux,
                      stale_halo_rejected)
                  .code == StatusCode::invalid_plan &&
              !stale_halo_rejected.valid(),
          "flux stage rejects a velocity certificate after its bound correction halo advances");
    }

    if (case_index == 2U) {
      const double corner = candidate_density.view.unchecked({-1, -1, -1}, 0U);
      candidate_density.view.unchecked({-1, -1, -1}, 0U) = 12345.0;
      PisoFrozenMomentumFluxStageCertificate corner_replay;
      passed &= expect(
          static_cast<bool>(c1.coupler.stage_frozen_momentum_flux(
              c1_stage, velocity_stage, as_const(candidate_density.view),
              candidate_flux, corner_replay)) &&
              corner_replay.canonical_lineage() ==
                  flux_stage.canonical_lineage(),
          "undefined edge/corner density poison does not enter candidate lineage");
      candidate_density.view.unchecked({-1, -1, -1}, 0U) = corner;
      const double pressure_corner =
          scaled_dp.view.unchecked({-1, -1, -1}, 0U);
      scaled_dp.view.unchecked({-1, -1, -1}, 0U) = -54321.0;
      PisoFrozenMomentumFluxStageCertificate pressure_corner_replay;
      passed &= expect(
          static_cast<bool>(c1.coupler.stage_frozen_momentum_flux(
              c1_stage, velocity_stage, as_const(candidate_density.view),
              candidate_flux, pressure_corner_replay)) &&
              pressure_corner_replay.canonical_lineage() ==
                  flux_stage.canonical_lineage(),
          "undefined edge/corner correction poison does not enter replay lineage");
      scaled_dp.view.unchecked({-1, -1, -1}, 0U) = pressure_corner;
      const double face_ghost = scaled_dp.view.unchecked({-1, 0, 0}, 0U);
      scaled_dp.view.unchecked({-1, 0, 0}, 0U) = face_ghost + 1.0;
      PisoFrozenMomentumFluxStageCertificate face_ghost_rejected;
      passed &= expect(
          c1.coupler
                  .stage_frozen_momentum_flux(
                      c1_stage, velocity_stage,
                      as_const(candidate_density.view), candidate_flux,
                      face_ghost_rejected)
                  .code == StatusCode::invalid_plan &&
              !face_ghost_rejected.valid(),
          "defined face-ghost mutation invalidates velocity/flux replay");
      scaled_dp.view.unchecked({-1, 0, 0}, 0U) = face_ghost;
    }
  }
  passed &= expect(state_lineages[0U] != state_lineages[1U] &&
                       state_lineages[0U] != state_lineages[2U] &&
                       state_lineages[1U] != state_lineages[2U] &&
                       flux_lineages[0U] != flux_lineages[1U] &&
                       flux_lineages[0U] != flux_lineages[2U] &&
                       flux_lineages[1U] != flux_lineages[2U],
                   "each C1 alpha/state candidate has unique canonical provenance");
  ConstFaceFluxView still_current;
  PisoFrozenMomentumStageAuthority still_current_stage;
  passed &= expect(c1.density.bytes == frozen_density_bytes &&
                       c1.h_by_a.bytes == frozen_h_by_a_bytes &&
                       face_flux_values(c1.intermediate_flux) ==
                           frozen_base_flux &&
                       static_cast<bool>(c1.coupler.inspect_intermediate_flux(
                           c1.intermediate, still_current)) &&
                       static_cast<bool>(
                           c1.coupler.make_frozen_momentum_stage_authority(
                               c1.intermediate, c1.pressure_certificate,
                               still_current_stage)),
                   "repeated C1 staging consumes no live state/intermediate/pressure authority");

  PisoFrozenMomentumPressureStageCertificate live_halo_rejected;
  passed &= expect(
      c1.coupler
              .form_frozen_momentum_scaled_pressure(
                  c1_stage, as_const(raw_dp.view), c1.correction_halo, 0.5,
                  scaled_dp.view, live_halo_rejected)
              .code == StatusCode::invalid_plan &&
          !live_halo_rejected.valid(),
      "candidate pressure stage rejects the coupler live correction halo");

  FullRefreshSample c2;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(raw_dp.view), as_const(dh.view), 0.0, 0.0,
                       11600U, true, c2),
                   "candidate mapping C2 byte-copy baseline refreshes");
  PressureReferenceCertificate c2_reference{
      fixture.equations.pressure_reference().fingerprint(),
      fixture.equations.thermophysical_predictor().fingerprint(),
      fixture.thermodynamics.fingerprint(),
      11646U,
      kFullRefreshTime,
      kFullRefreshPressureReferenceRevision,
      PressureReferenceKind::closed_mass};
  const PressureCorrectionInput c2_pressure_input{
      c2.intermediate,
      c2_reference,
      as_const(c2.density.view),
      as_const(c2.density_accepted.view),
      as_const(c2.density_previous.view),
      as_const(c2.pressure_compressibility.view),
      kFullRefreshBdf,
      kFullRefreshTime,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision()};
  passed &= expect(static_cast<bool>(c2.coupler.assemble_pressure_system(
                       c2_pressure_input,
                       {c2.pressure_diagonal.view, c2.pressure_rhs.view},
                       c2.pressure_certificate)),
                   "candidate mapping C2 pressure authority assembles");
  PisoFrozenMomentumStageAuthority c2_stage;
  passed &= expect(
      static_cast<bool>(c2.coupler.make_frozen_momentum_stage_authority(
          c2.intermediate, c2.pressure_certificate, c2_stage)) &&
          c2_stage.corrector() == 2U,
      "C2 issues its own frozen-momentum authority");

  OwnedField c2_scaled = make_field(kFrozenCandidateCorrectionField, cells,
                                    1U, 1U,
                                    11650U, 12650U);
  OwnedField c2_velocity = make_field(kFrozenCandidateVelocityField, cells,
                                      3U, 1U, 11651U, 12651U);
  OwnedField c2_density = make_field(kFrozenCandidateDensityField, cells, 1U,
                                     1U, 11652U, 12652U);
  FaceFluxStorage c2_flux_storage;
  FaceFluxView c2_flux;
  HaloEngine c2_halo;
  passed &= expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(cells, 1U,
                                                             c2_flux_storage)) &&
          static_cast<bool>(c2_flux_storage.workspace_view(0U, 11653U,
                                                           c2_flux)) &&
          static_cast<bool>(c2_halo.reserve(
              MPI_COMM_SELF, fixture.patch,
              {correction_spec.data(), correction_spec.size()},
              fixture.boundary.halo_topology())),
      "C2 candidate scratch allocates");
  fill(raw_dp, 0.0);
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const Int3 cell{x, y, z};
        const unsigned outside =
            static_cast<unsigned>(x < 0 || x >= cells.x) +
            static_cast<unsigned>(y < 0 || y >= cells.y) +
            static_cast<unsigned>(z < 0 || z >= cells.z);
        if (outside > 1U) continue;
        const Int3 wrapped = full_refresh_wrap(cell, cells);
        c2_density.view.unchecked(cell, 0U) =
            c2.density.view.unchecked(wrapped, 0U) * density_factor(wrapped);
      }
  PisoFrozenMomentumPressureStageCertificate c2_pressure_stage;
  PisoFrozenMomentumVelocityStageCertificate c2_velocity_stage;
  PisoFrozenMomentumFluxStageCertificate c2_flux_stage;
  passed &= expect(
      static_cast<bool>(c2.coupler.form_frozen_momentum_scaled_pressure(
          c2_stage, as_const(raw_dp.view), c2_halo, 1.0,
          c2_scaled.view, c2_pressure_stage)),
      "C2 forms zero pressure direction with nonzero density direction");
  std::array<FieldView, 1U> c2_halo_views{c2_scaled.view};
  HaloTicket c2_ticket;
  Status c2_halo_status =
      c2_halo.begin(318U, {c2_halo_views.data(), c2_halo_views.size()},
                    c2_ticket);
  if (c2_halo_status)
    c2_halo_status = c2_halo.finish(
        c2_ticket, {c2_halo_views.data(), c2_halo_views.size()});
  c2_scaled.view = c2_halo_views[0U];
  passed &= expect(
      static_cast<bool>(c2_halo_status) &&
          static_cast<bool>(c2.coupler.stage_frozen_momentum_velocity(
              c2_stage, c2_pressure_stage, c2_halo,
              as_const(c2_scaled.view), c2_velocity.view,
              c2_velocity_stage)) &&
          static_cast<bool>(c2.coupler.stage_frozen_momentum_flux(
              c2_stage, c2_velocity_stage, as_const(c2_density.view), c2_flux,
              c2_flux_stage)),
      "C2 stages density-refreshed total flux from explicit offset");
  double c2_density_response = 0.0;
  double c2_oracle_error = 0.0;
  for (CartesianAxis axis : {CartesianAxis::x, CartesianAxis::y,
                             CartesianAxis::z}) {
    const ConstFaceFieldView base = c2.intermediate_flux.x.axis == axis
                                        ? c2.intermediate_flux.x
                                        : (axis == CartesianAxis::y
                                               ? c2.intermediate_flux.y
                                               : c2.intermediate_flux.z);
    const ConstFaceFieldView actual = full_refresh_face(as_const(c2_flux), axis);
    const std::uint8_t component = component_index(axis);
    for (std::int32_t z = 0; z < actual.extents.z; ++z)
      for (std::int32_t y = 0; y < actual.extents.y; ++y)
        for (std::int32_t x = 0; x < actual.extents.x; ++x) {
          const Int3 face{x, y, z};
          const Int3 left = full_refresh_wrap(shift(face, axis, -1), cells);
          const Int3 right = full_refresh_wrap(face, cells);
          const double left_delta =
              c2.density.view.unchecked(left, 0U) *
              (density_factor(left) - 1.0);
          const double right_delta =
              c2.density.view.unchecked(right, 0U) *
              (density_factor(right) - 1.0);
          const double expected =
              base.unchecked(face) +
              face_area(fixture, axis, right) * 0.5 *
                  (left_delta * c2.h_by_a.view.unchecked(left, component) +
                   right_delta * c2.h_by_a.view.unchecked(right, component));
          c2_density_response = std::max(
              c2_density_response,
              std::abs(actual.unchecked(face) - base.unchecked(face)));
          c2_oracle_error = std::max(
              c2_oracle_error, std::abs(actual.unchecked(face) - expected));
        }
  }
  passed &= expect(c2_density_response > 1.0e-8 && c2_oracle_error < 2.0e-12,
                   "C2 explicit Phi0-Q(rho0) offset preserves nonzero exact density sensitivity");

  PisoFrozenMomentumPressureStageCertificate foreign_authority_rejected;
  passed &= expect(
      c1.coupler
              .form_frozen_momentum_scaled_pressure(
                  c2_stage, as_const(raw_dp.view), candidate_halo, 0.5,
                  scaled_dp.view, foreign_authority_rejected)
              .code == StatusCode::invalid_plan &&
          !foreign_authority_rejected.valid(),
      "C1 coupler rejects a foreign C2 stage authority");
  PressureCorrectionCertificate discarded;
  OwnedField poison_diagonal =
      make_field(180U, cells, 1U, 0U, 11660U, 12660U);
  OwnedField poison_rhs =
      make_field(181U, cells, 1U, 0U, 11661U, 12661U);
  (void)c1.coupler.assemble_pressure_system(
      {}, {poison_diagonal.view, poison_rhs.view}, discarded);
  PisoFrozenMomentumPressureStageCertificate stale_rejected;
  passed &= expect(
      c1.coupler
              .form_frozen_momentum_scaled_pressure(
                  c1_stage, as_const(raw_dp.view), candidate_halo, 0.5,
                  scaled_dp.view, stale_rejected)
              .code == StatusCode::invalid_plan &&
          !stale_rejected.valid(),
      "cleared pressure lifecycle makes the former candidate authority stale");
  return passed;
}

bool test_frozen_momentum_exact_publication_authority() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, false, true),
                       "frozen exact publication periodic fixture compiles");
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "frozen exact publication PISO plan compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(initialize_flux_history(cells, history) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        11701U, history, previous_flux) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        11702U, history, accepted_flux),
                   "frozen exact publication BDF2 history initializes");
  OwnedField raw_dp = make_field(kFullRefreshCorrectionField, cells, 1U, 1U,
                                 11703U, 12703U);
  OwnedField raw_dh = make_field(kFullRefreshEnthalpyDirectionField, cells,
                                 1U, 2U, 11704U, 12704U);
  for (std::int32_t z = -1; z <= cells.z; ++z)
    for (std::int32_t y = -1; y <= cells.y; ++y)
      for (std::int32_t x = -1; x <= cells.x; ++x) {
        const Int3 cell{x, y, z};
        const Int3 wrapped = full_refresh_wrap(cell, cells);
        raw_dp.view.unchecked(cell, 0U) =
            250.0 *
            (std::sin(full_refresh_phase(wrapped.x, cells.x, false)) -
             0.25 * std::cos(
                        full_refresh_phase(wrapped.y, cells.y, false)));
        raw_dh.view.unchecked(cell, 0U) =
            3.0 * std::cos(
                      full_refresh_phase(wrapped.z, cells.z, false)) +
            0.4 * std::sin(full_refresh_phase(wrapped.x, cells.x, false));
      }
  constexpr RevisionToken sample_seed = 11750U;
  FullRefreshSample sample;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(raw_dp.view), as_const(raw_dh.view), 0.0, 0.0,
                       sample_seed, false, sample),
                   "frozen exact publication base refreshes");
  if (!passed) return false;
  PisoFrozenMomentumStageAuthority authority;
  passed &= expect(
      static_cast<bool>(sample.coupler.make_frozen_momentum_stage_authority(
          sample.intermediate, sample.pressure_certificate, authority)),
      "frozen exact publication stage authority issues");
  ReductionEngine reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_SELF, ReductionMode::mpi_allreduce, 2U,
                       reductions)),
                   "frozen exact publication reductions compile");
  if (!passed) return false;
  const PressureReferenceCertificate predecessor =
      full_refresh_pressure_reference(fixture, sample_seed);

  FrozenExactScratch baseline_scratch;
  passed &= expect(
      stage_frozen_exact_scratch(
          fixture, sample, authority, as_const(raw_dp.view),
          as_const(raw_dh.view), 0.0, 11800U, baseline_scratch) &&
          prepare_frozen_exact_thermodynamics(
              fixture, sample, authority, predecessor,
              kFullRefreshPressureReference,
              11850U, baseline_scratch, reductions),
      "alpha-zero full exact state stages and closes its gauge");
  PisoFrozenMomentumExactCandidateCertificate exact_baseline;
  passed &= expect(
      static_cast<bool>(sample.coupler.certify_frozen_momentum_exact_baseline(
          authority, baseline_scratch.pressure_stage,
          baseline_scratch.velocity_stage, baseline_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), baseline_scratch),
          reductions, exact_baseline)) &&
          exact_baseline.valid() && exact_baseline.alpha() == 0.0 &&
          exact_baseline.baseline_state_provenance() ==
              exact_baseline.candidate_state_provenance() &&
          exact_baseline.baseline_mass_flux_provenance() ==
              exact_baseline.candidate_mass_flux_provenance() &&
          exact_baseline.base_state_provenance() != 0U,
      "alpha-zero exact certificate is the canonical selector baseline");
  PisoFrozenMomentumExactCandidateInput closed_without_compressibility =
      frozen_exact_input(sample, as_const(raw_dh.view), baseline_scratch);
  closed_without_compressibility.thermodynamic.pressure_compressibility = {};
  PisoFrozenMomentumExactCandidateCertificate cleared_closed_baseline =
      exact_baseline;
  passed &= expect(
      sample.coupler
                  .certify_frozen_momentum_exact_baseline(
                      authority, baseline_scratch.pressure_stage,
                      baseline_scratch.velocity_stage,
                      baseline_scratch.flux_stage,
                      closed_without_compressibility, reductions,
                      cleared_closed_baseline)
                  .code == StatusCode::invalid_plan &&
          !cleared_closed_baseline.valid(),
      "closed exact baseline rejects missing pressure-compressibility authority and clears a reused certificate");
  const double baseline_temperature =
      baseline_scratch.temperature.view.unchecked({0, 0, 0}, 0U);
  baseline_scratch.temperature.view.unchecked({0, 0, 0}, 0U) =
      baseline_temperature + 1.0;
  PisoFrozenMomentumExactCandidateCertificate reused_exact_baseline =
      exact_baseline;
  passed &= expect(
      !sample.coupler.certify_frozen_momentum_exact_baseline(
          authority, baseline_scratch.pressure_stage,
          baseline_scratch.velocity_stage, baseline_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), baseline_scratch),
          reductions, reused_exact_baseline) &&
          !reused_exact_baseline.valid(),
      "exact baseline rejects a positive but non-base temperature and clears a reused output certificate");
  baseline_scratch.temperature.view.unchecked({0, 0, 0}, 0U) =
      baseline_temperature;
  const Int3 stationary_probe{0, 0, 0};
  const double stationary_velocity =
      baseline_scratch.velocity.view.unchecked(stationary_probe, 0U);
  baseline_scratch.velocity.view.unchecked(stationary_probe, 0U) =
      stationary_velocity + 1.0;
  PressureEnergyStationaryCertificate nonbaseline_stationary;
  passed &= expect(
      sample.coupler
                  .certify_frozen_momentum_stationary(
                      authority, exact_baseline, 0.0, 0.0, reductions,
                      nonbaseline_stationary)
                  .code == StatusCode::invalid_plan &&
          !nonbaseline_stationary.valid(),
      "stationary authority rejects an alpha-zero state that is not byte-equal to the live baseline");
  baseline_scratch.velocity.view.unchecked(stationary_probe, 0U) =
      stationary_velocity;
  PressureEnergyStationaryCertificate stationary;
  const Status stationary_status =
      sample.coupler.certify_frozen_momentum_stationary(
          authority, exact_baseline, 0.0, 0.0, reductions, stationary);
  passed &= expect(
      static_cast<bool>(stationary_status) &&
          stationary.valid(),
      "zero direction and terminal residuals issue typed stationary authority");
  PressureEnergyStationaryCertificate reused_stationary = stationary;
  passed &= expect(
      !sample.coupler.certify_frozen_momentum_stationary(
          authority, exact_baseline, -1.0, 0.0, reductions,
          reused_stationary) &&
          !reused_stationary.valid(),
      "stationary failure clears a reused terminal certificate without consuming exact authority");
  const std::vector<double> pressure_before = sample.pressure.bytes;
  const std::vector<double> enthalpy_before = sample.enthalpy.bytes;
  const std::vector<double> density_before = sample.density.bytes;
  const std::vector<double> temperature_before = sample.temperature.bytes;
  const std::vector<double> velocity_before = sample.velocity.bytes;
  const std::vector<double> output_before =
      face_flux_values(as_const(sample.c1_trial_flux));
  PisoStateCorrectionCertificate rejected_alpha_zero;
  passed &= expect(
      sample.coupler
                  .commit_frozen_momentum_coupled_trial_state(
                      authority, exact_baseline,
                      PressureEnergyGlobalizationSelectionCertificate{},
                      {sample.velocity.view, sample.pressure.view,
                       sample.enthalpy.view, sample.density.view,
                       sample.temperature.view},
                      sample.c1_trial_flux, reductions, rejected_alpha_zero)
                  .code == StatusCode::invalid_plan &&
          !rejected_alpha_zero.valid() &&
          sample.pressure.bytes == pressure_before &&
          sample.enthalpy.bytes == enthalpy_before &&
          sample.density.bytes == density_before &&
          sample.temperature.bytes == temperature_before &&
          sample.velocity.bytes == velocity_before &&
          face_flux_values(as_const(sample.c1_trial_flux)) == output_before,
      "alpha-zero exact stage cannot publish without stationary authority");

  FrozenExactScratch selected_scratch;
  passed &= expect(
      stage_frozen_exact_scratch(
          fixture, sample, authority, as_const(raw_dp.view),
          as_const(raw_dh.view), 0.5, 11900U, selected_scratch) &&
          prepare_frozen_exact_thermodynamics(
              fixture, sample, authority, predecessor,
              kFullRefreshPressureReference,
              11950U, selected_scratch, reductions),
      "half-step full exact state stages and closes its gauge");
  PisoFrozenMomentumExactCandidateCertificate exact_selected;
  passed &= expect(
      static_cast<bool>(sample.coupler.certify_frozen_momentum_exact_candidate(
          authority, exact_baseline, selected_scratch.pressure_stage,
          selected_scratch.velocity_stage, selected_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), selected_scratch),
          reductions, exact_selected)) &&
          exact_selected.valid() && exact_selected.alpha() == 0.5 &&
          exact_selected.correction_direction() ==
              exact_baseline.correction_direction() &&
          exact_selected.baseline_state_provenance() ==
              exact_baseline.candidate_state_provenance() &&
          exact_selected.baseline_mass_flux_provenance() ==
              exact_baseline.candidate_mass_flux_provenance(),
      "selected exact state binds joint raw dp/dh direction and exact baseline");
  OwnedField foreign_base_velocity =
      make_field(sample.velocity.view.field, cells, 3U, 0U,
                 sample.velocity.view.revision, 11998U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component)
          foreign_base_velocity.view.unchecked({x, y, z}, component) =
              sample.velocity.view.unchecked({x, y, z}, component);
  PisoFrozenMomentumExactCandidateInput foreign_base_input =
      frozen_exact_input(sample, as_const(raw_dh.view), selected_scratch);
  foreign_base_input.base_state.velocity = foreign_base_velocity.view;
  PisoFrozenMomentumExactCandidateCertificate foreign_base_rejected =
      exact_selected;
  passed &= expect(
      !sample.coupler.certify_frozen_momentum_exact_candidate(
          authority, exact_baseline, selected_scratch.pressure_stage,
          selected_scratch.velocity_stage, selected_scratch.flux_stage,
          foreign_base_input, reductions, foreign_base_rejected) &&
          !foreign_base_rejected.valid(),
      "selected exact candidate rejects a numeric-equal foreign base identity and clears reused output");
  PisoFrozenMomentumExactCandidateCertificate reused_exact_candidate =
      exact_selected;
  passed &= expect(
      !sample.coupler.certify_frozen_momentum_exact_candidate(
          authority, PisoFrozenMomentumExactCandidateCertificate{},
          selected_scratch.pressure_stage, selected_scratch.velocity_stage,
          selected_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), selected_scratch),
          reductions, reused_exact_candidate) &&
          !reused_exact_candidate.valid(),
      "selected exact failure clears a reused candidate certificate and preserves the prior current replay");
  PressureEnergyGlobalizationSample baseline_sample;
  baseline_sample.alpha = 0.0;
  baseline_sample.global_normalized_continuity = 1.0;
  baseline_sample.global_normalized_energy = 1.0;
  baseline_sample.thermodynamically_admissible = true;
  baseline_sample.state_and_flux_finite = true;
  baseline_sample.corrector = exact_baseline.corrector();
  baseline_sample.target_time = exact_baseline.target_time();
  baseline_sample.correction_direction =
      exact_baseline.correction_direction();
  baseline_sample.state_provenance =
      exact_baseline.candidate_state_provenance();
  baseline_sample.mass_flux_provenance =
      exact_baseline.candidate_mass_flux_provenance();
  std::array<PressureEnergyGlobalizationSample,
             kPressureEnergyGlobalizationCandidateCount>
      candidates{};
  for (std::size_t index = 0U; index < candidates.size(); ++index) {
    candidates[index] = baseline_sample;
    candidates[index].alpha = std::ldexp(1.0, -static_cast<int>(index));
    candidates[index].global_normalized_continuity = 2.0;
    candidates[index].global_normalized_energy = 2.0;
    candidates[index].state_provenance = UINT64_C(0x71000000) + index;
    candidates[index].mass_flux_provenance = UINT64_C(0x72000000) + index;
  }
  candidates[1U].global_normalized_continuity = 0.2;
  candidates[1U].global_normalized_energy = 0.2;
  candidates[1U].state_provenance =
      exact_selected.candidate_state_provenance();
  candidates[1U].mass_flux_provenance =
      exact_selected.candidate_mass_flux_provenance();
  PressureEnergyGlobalizationSelectionCertificate selection;
  passed &= expect(
      static_cast<bool>(select_pressure_energy_globalization(
          baseline_sample, {candidates.data(), candidates.size()},
          selection)) &&
          selection.valid() && selection.alpha == 0.5,
      "selector signs only the matching half-step exact state and flux");

  const Int3 probe{0, 0, 0};
  const double raw_h_before = raw_dh.view.unchecked(probe, 0U);
  raw_dh.view.unchecked(probe, 0U) = raw_h_before + 1.0;
  PisoStateCorrectionCertificate raw_h_rejected;
  const Status raw_h_status =
      sample.coupler.commit_frozen_momentum_coupled_trial_state(
          authority, exact_selected, selection,
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          sample.c1_trial_flux, reductions, raw_h_rejected);
  passed &= expect(
      !raw_h_status &&
          !raw_h_rejected.valid() && sample.pressure.bytes == pressure_before &&
          face_flux_values(as_const(sample.c1_trial_flux)) == output_before,
      "raw enthalpy direction mutation rejects with zero publication");
  raw_dh.view.unchecked(probe, 0U) = raw_h_before;
  const double scaled_p_before =
      selected_scratch.scaled_pressure.view.unchecked(probe, 0U);
  selected_scratch.scaled_pressure.view.unchecked(probe, 0U) =
      raw_dp.view.unchecked(probe, 0U);
  PisoStateCorrectionCertificate raw_as_scaled_rejected;
  passed &= expect(
      sample.coupler
                  .commit_frozen_momentum_coupled_trial_state(
                      authority, exact_selected, selection,
                      {sample.velocity.view, sample.pressure.view,
                       sample.enthalpy.view, sample.density.view,
                       sample.temperature.view},
                      sample.c1_trial_flux, reductions,
                      raw_as_scaled_rejected)
                  .code == StatusCode::invalid_plan &&
          !raw_as_scaled_rejected.valid() &&
          sample.pressure.bytes == pressure_before &&
          face_flux_values(as_const(sample.c1_trial_flux)) == output_before,
      "raw dp cannot be substituted for the selected scaled dp");
  selected_scratch.scaled_pressure.view.unchecked(probe, 0U) =
      scaled_p_before;

  const PisoFrozenMomentumExactCandidateCertificate old_exact = exact_selected;
  passed &= expect(
      prepare_frozen_exact_thermodynamics(
          fixture, sample, authority, predecessor,
          kFullRefreshPressureReference,
          11951U, selected_scratch, reductions),
      "same numeric candidate refreshes a distinct EOS scratch closure");
  PisoFrozenMomentumExactCandidateCertificate replay_exact;
  passed &= expect(
      static_cast<bool>(sample.coupler.certify_frozen_momentum_exact_candidate(
          authority, exact_baseline, selected_scratch.pressure_stage,
          selected_scratch.velocity_stage, selected_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), selected_scratch),
          reductions, replay_exact)) &&
          replay_exact.correction_direction() ==
              old_exact.correction_direction() &&
          replay_exact.candidate_state_provenance() ==
              old_exact.candidate_state_provenance() &&
          replay_exact.candidate_mass_flux_provenance() ==
              old_exact.candidate_mass_flux_provenance() &&
          replay_exact.canonical_lineage() == old_exact.canonical_lineage() &&
          replay_exact.scratch_binding() != old_exact.scratch_binding(),
      "same values under a new closure revision preserve canonical provenance but change scratch binding");
  PisoStateCorrectionCertificate stale_exact_rejected;
  passed &= expect(
      sample.coupler
                  .commit_frozen_momentum_coupled_trial_state(
                      authority, old_exact, selection,
                      {sample.velocity.view, sample.pressure.view,
                       sample.enthalpy.view, sample.density.view,
                       sample.temperature.view},
                      sample.c1_trial_flux, reductions, stale_exact_rejected)
                  .code == StatusCode::invalid_plan &&
          !stale_exact_rejected.valid() &&
          sample.pressure.bytes == pressure_before &&
          face_flux_values(as_const(sample.c1_trial_flux)) == output_before,
      "stale scratch-bound exact replay cannot publish canonical-equal values");

  const double material_before = sample.viscosity.view.unchecked(probe, 0U);
  sample.viscosity.view.unchecked(probe, 0U) = material_before * 1.125;
  PisoFrozenMomentumExactCandidateCertificate material_replay;
  passed &= expect(
      static_cast<bool>(sample.coupler.certify_frozen_momentum_exact_candidate(
          authority, exact_baseline, selected_scratch.pressure_stage,
          selected_scratch.velocity_stage, selected_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), selected_scratch),
          reductions, material_replay)) &&
          material_replay.candidate_state_provenance() ==
              replay_exact.candidate_state_provenance() &&
          material_replay.candidate_mass_flux_provenance() ==
              replay_exact.candidate_mass_flux_provenance(),
      "material-only replay is outside five-field selector provenance");
  sample.viscosity.view.unchecked(probe, 0U) = material_before;

  const auto rejects_unrevisioned_base_mutation =
      [&](OwnedField& field, std::uint8_t component,
          std::string_view description) {
        const double original = field.view.unchecked(probe, component);
        field.view.unchecked(probe, component) = original + 2.0;
        const std::vector<double> pressure_at_entry = sample.pressure.bytes;
        const std::vector<double> enthalpy_at_entry = sample.enthalpy.bytes;
        const std::vector<double> density_at_entry = sample.density.bytes;
        const std::vector<double> temperature_at_entry =
            sample.temperature.bytes;
        const std::vector<double> velocity_at_entry = sample.velocity.bytes;
        const std::vector<double> flux_at_entry =
            face_flux_values(as_const(sample.c1_trial_flux));
        PisoStateCorrectionCertificate rejected;
        const Status rejected_status =
            sample.coupler.commit_frozen_momentum_coupled_trial_state(
                authority, material_replay, selection,
                {sample.velocity.view, sample.pressure.view,
                 sample.enthalpy.view, sample.density.view,
                 sample.temperature.view},
                sample.c1_trial_flux, reductions, rejected);
        const bool unchanged =
            !rejected_status && !rejected.valid() &&
            sample.pressure.bytes == pressure_at_entry &&
            sample.enthalpy.bytes == enthalpy_at_entry &&
            sample.density.bytes == density_at_entry &&
            sample.temperature.bytes == temperature_at_entry &&
            sample.velocity.bytes == velocity_at_entry &&
            face_flux_values(as_const(sample.c1_trial_flux)) == flux_at_entry;
        field.view.unchecked(probe, component) = original;
        passed &= expect(unchanged, description);
      };
  rejects_unrevisioned_base_mutation(
      sample.velocity, 0U,
      "post-certificate unrevisioned base U mutation is collectively rejected with zero writes");
  rejects_unrevisioned_base_mutation(
      sample.density, 0U,
      "post-certificate unrevisioned base rho mutation is collectively rejected with zero writes");
  rejects_unrevisioned_base_mutation(
      sample.temperature, 0U,
      "post-certificate unrevisioned base T mutation is collectively rejected with zero writes");
  rejects_unrevisioned_base_mutation(
      sample.pressure, 0U,
      "post-certificate unrevisioned base p mutation is collectively rejected with zero writes");
  rejects_unrevisioned_base_mutation(
      sample.enthalpy, 0U,
      "post-certificate unrevisioned base h mutation is collectively rejected with zero writes");

  PisoStateCorrectionCertificate committed;
  const auto same_owned_cells = [&](ConstFieldView left, ConstFieldView right,
                                    std::uint8_t components) {
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          for (std::uint8_t component = 0U; component < components;
               ++component)
            if (std::memcmp(
                    &left.unchecked({x, y, z}, component),
                    &right.unchecked({x, y, z}, component), sizeof(double)) !=
                0)
              return false;
    return true;
  };
  passed &= expect(
      static_cast<bool>(sample.coupler.commit_frozen_momentum_coupled_trial_state(
          authority, material_replay, selection,
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          sample.c1_trial_flux, reductions, committed)) &&
          committed.valid() && committed.corrector == 1U &&
          committed.correction ==
              selected_scratch.scaled_pressure.view.revision &&
          committed.correction != raw_dp.view.revision &&
          same_owned_cells(as_const(sample.pressure.view),
                           as_const(selected_scratch.pressure.view), 1U) &&
          same_owned_cells(as_const(sample.enthalpy.view),
                           as_const(selected_scratch.enthalpy.view), 1U) &&
          same_owned_cells(as_const(sample.density.view),
                           as_const(selected_scratch.density.view), 1U) &&
          same_owned_cells(as_const(sample.temperature.view),
                           as_const(selected_scratch.temperature.view), 1U) &&
          same_owned_cells(as_const(sample.velocity.view),
                           as_const(selected_scratch.velocity.view), 3U) &&
          face_flux_values(as_const(sample.c1_trial_flux)) ==
              face_flux_values(as_const(selected_scratch.flux)),
      "selected replay atomically publishes p/h/rho/T/U/flux and records scaled-dp revision");
  PisoStateCorrectionCertificate reused_commit = committed;
  passed &= expect(
      !sample.coupler.commit_frozen_momentum_coupled_trial_state(
          authority, material_replay, selection,
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          sample.c1_trial_flux, reductions, reused_commit) &&
          !reused_commit.valid(),
      "failed commit after authority consumption clears a reused correction certificate");
  PisoFrozenMomentumStageAuthority consumed;
  passed &= expect(
      sample.coupler
                  .make_frozen_momentum_stage_authority(
                      sample.intermediate, sample.pressure_certificate,
                      consumed)
                  .code == StatusCode::invalid_plan &&
          !consumed.valid(),
      "successful exact publication alone consumes current stage authority");
  return passed;
}

bool test_frozen_momentum_stationary_pending_lineage() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, false, true),
                       "stationary pending periodic fixture compiles");
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "stationary pending PISO plan compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(initialize_flux_history(cells, history) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        12001U, history, previous_flux) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        12002U, history, accepted_flux),
                   "stationary pending BDF2 history initializes");
  OwnedField raw_dp = make_field(kFullRefreshCorrectionField, cells, 1U, 1U,
                                 12003U, 13003U);
  OwnedField raw_dh = make_field(kFullRefreshEnthalpyDirectionField, cells,
                                 1U, 2U, 12004U, 13004U);
  fill(raw_dp, 0.0);
  fill(raw_dh, 0.0);
  constexpr RevisionToken sample_seed = 12050U;
  FullRefreshSample sample;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(raw_dp.view), as_const(raw_dh.view), 0.0, 0.0,
                       sample_seed, false, sample),
                   "stationary pending C1 base refreshes");
  ReductionEngine reductions;
  passed &= expect(static_cast<bool>(ReductionEngine::compile(
                       MPI_COMM_SELF, ReductionMode::mpi_allreduce, 2U,
                       reductions)),
                   "stationary pending reductions compile");
  if (!passed) return false;

  PisoFrozenMomentumStageAuthority c1_authority;
  passed &= expect(
      static_cast<bool>(sample.coupler.make_frozen_momentum_stage_authority(
          sample.intermediate, sample.pressure_certificate, c1_authority)),
      "stationary pending C1 stage authority issues");
  FrozenExactScratch c1_scratch;
  const PressureReferenceCertificate c1_predecessor =
      full_refresh_pressure_reference(fixture, sample_seed);
  passed &= expect(
      stage_frozen_exact_scratch(
          fixture, sample, c1_authority, as_const(raw_dp.view),
          as_const(raw_dh.view), 0.0, 12100U, c1_scratch) &&
          prepare_frozen_exact_thermodynamics(
              fixture, sample, c1_authority, c1_predecessor,
              kFullRefreshPressureReference,
              12150U, c1_scratch, reductions),
      "stationary pending C1 exact baseline stages");
  PisoFrozenMomentumExactCandidateCertificate c1_exact;
  PressureEnergyStationaryCertificate c1_stationary;
  passed &= expect(
      static_cast<bool>(sample.coupler.certify_frozen_momentum_exact_baseline(
          c1_authority, c1_scratch.pressure_stage,
          c1_scratch.velocity_stage, c1_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), c1_scratch),
          reductions, c1_exact)) &&
          static_cast<bool>(sample.coupler.certify_frozen_momentum_stationary(
              c1_authority, c1_exact, 0.0, 0.0, reductions,
              c1_stationary)),
      "stationary pending C1 exact baseline reaches typed terminal gate");
  PisoStateCorrectionCertificate corrected_one;
  passed &= expect(
      static_cast<bool>(sample.coupler.commit_frozen_momentum_stationary_trial_state(
          c1_authority, c1_exact, c1_stationary,
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          sample.c1_trial_flux, reductions, corrected_one)) &&
          corrected_one.valid() && corrected_one.corrector == 1U,
      "stationary C1 commits only through its typed terminal authority");
  if (!passed) return false;

  PisoIntermediateInput c2_input = sample.refresh_input;
  c2_input.corrector = 2U;
  c2_input.temporal_reference = {};
  c2_input.committed_face_history = {};
  c2_input.prior_corrector = corrected_one.state;
  c2_input.pressure_reference = corrected_one.output_pressure_reference;
  c2_input.thermophysical_boundary.binding.pressure_reference =
      c1_scratch.thermodynamic.closed_gauge.next_pressure_reference;
  c2_input.trial_velocity = as_const(sample.velocity.view);
  c2_input.density = sample.density.view;
  ConstFaceFluxView c2_trial = as_const(sample.c1_trial_flux);
  c2_trial.revision = corrected_one.state;
  c2_input.trial_flux = c2_trial;
  ++c2_input.momentum.state;
  PisoIntermediateCertificate c2_intermediate;
  passed &= expect(static_cast<bool>(sample.coupler.refresh(
                       c2_input, c2_intermediate)) &&
                       c2_intermediate.corrector == 2U,
                   "stationary pending C2 refresh consumes corrected C1 lineage");
  const PressureCorrectionInput c2_pressure_input{
      c2_intermediate,
      c2_input.pressure_reference,
      as_const(sample.density.view),
      as_const(sample.density_accepted.view),
      as_const(sample.density_previous.view),
      as_const(sample.pressure_compressibility.view),
      kFullRefreshBdf,
      kFullRefreshTime,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision()};
  sample.intermediate = c2_intermediate;
  passed &= expect(static_cast<bool>(sample.coupler.assemble_pressure_system(
                       c2_pressure_input,
                       {sample.pressure_diagonal.view, sample.pressure_rhs.view},
                       sample.pressure_certificate)) &&
                       sample.pressure_certificate.corrector == 2U,
                   "stationary pending C2 pressure authority assembles");
  PisoFrozenMomentumStageAuthority c2_authority;
  passed &= expect(
      static_cast<bool>(sample.coupler.make_frozen_momentum_stage_authority(
          sample.intermediate, sample.pressure_certificate, c2_authority)) &&
          c2_authority.corrector() == 2U,
      "stationary pending C2 stage authority issues");
  FrozenExactScratch c2_scratch;
  passed &= expect(
      stage_frozen_exact_scratch(
          fixture, sample, c2_authority, as_const(raw_dp.view),
          as_const(raw_dh.view), 0.0, 12200U, c2_scratch) &&
          prepare_frozen_exact_thermodynamics(
              fixture, sample, c2_authority,
              corrected_one.output_pressure_reference,
              c1_scratch.thermodynamic.closed_gauge.next_pressure_reference,
              12250U, c2_scratch, reductions),
      "stationary pending C2 exact baseline stages from corrected C1 flux");
  PisoFrozenMomentumExactCandidateCertificate c2_exact;
  PressureEnergyStationaryCertificate c2_stationary;
  passed &= expect(
      static_cast<bool>(sample.coupler.certify_frozen_momentum_exact_baseline(
          c2_authority, c2_scratch.pressure_stage,
          c2_scratch.velocity_stage, c2_scratch.flux_stage,
          frozen_exact_input(sample, as_const(raw_dh.view), c2_scratch),
          reductions, c2_exact)) &&
          static_cast<bool>(sample.coupler.certify_frozen_momentum_stationary(
              c2_authority, c2_exact, 0.0, 0.0, reductions,
              c2_stationary)),
      "stationary pending C2 exact baseline reaches typed terminal gate");

  FieldRegistry registry;
  FieldSchema schema;
  FieldId dependency = 0U;
  passed &= expect(static_cast<bool>(registry.declare_field(
                       "frozen.pending", 1U, 0U, dependency)) &&
                       static_cast<bool>(registry.freeze(schema)),
                   "stationary pending transaction schema freezes");
  const std::array requests{ArenaFieldRequest{
      dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  StateLayers layers;
  AttemptTransaction transaction;
  FaceFluxStorage final_storage;
  FinalFaceFluxAuthority final_authority;
  FinalFaceFluxWriter final_writer;
  PendingFaceFluxView pending;
  passed &= expect(
      static_cast<bool>(ArenaLayout::compile(
          schema, {requests.data(), requests.size()}, layout)) &&
          static_cast<bool>(StateLayers::allocate(layout, layers)) &&
          static_cast<bool>(AttemptTransaction::create(
              layers.field_count(), 1U, layers.field_count(), transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(final_authority.claim(
              piso.pressure_stage(), piso.final_flux_slot(), transaction,
              final_writer)) &&
          static_cast<bool>(transaction.begin(layers)) &&
          static_cast<bool>(transaction.revise_trial(dependency)) &&
          static_cast<bool>(final_writer.begin_pending(
              transaction, final_storage, pending)),
      "stationary pending C2 acquires opaque final-flux payload");
  PisoStateCorrectionCertificate corrected_two = corrected_one;
  passed &= expect(
      !sample.coupler.commit_frozen_momentum_stationary_pending_state(
          c2_authority, c2_exact, PressureEnergyStationaryCertificate{},
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          pending, reductions, corrected_two) &&
          !corrected_two.valid(),
      "failed C2 pending entry clears a reused correction while preserving stage authority");
  passed &= expect(
      static_cast<bool>(sample.coupler.commit_frozen_momentum_stationary_pending_state(
          c2_authority, c2_exact, c2_stationary,
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          pending, reductions, corrected_two)) &&
          corrected_two.valid() && corrected_two.corrector == 2U &&
          corrected_two.face_flux == pending.revision(),
      "stationary C2 writes selected flux directly into opaque pending payload");
  ConstFaceFluxView inspected;
  passed &= expect(
      static_cast<bool>(sample.coupler.inspect_corrected_pending(
          corrected_two, pending, inspected)) &&
          inspected.revision == corrected_two.face_flux &&
          face_flux_values(inspected) ==
              face_flux_values(as_const(c2_scratch.flux)),
      "C2 frozen commit establishes corrected-pending terminal lineage without recomputing flux");
  const PisoStateCorrectionCertificate published_two = corrected_two;
  corrected_two = published_two;
  passed &= expect(
      !sample.coupler.commit_frozen_momentum_stationary_pending_state(
          c2_authority, c2_exact, PressureEnergyStationaryCertificate{},
          {sample.velocity.view, sample.pressure.view, sample.enthalpy.view,
           sample.density.view, sample.temperature.view},
          pending, reductions, corrected_two) &&
          !corrected_two.valid(),
      "failed C2 pending re-entry clears caller and internal corrected-pending authority");
  ConstFaceFluxView cleared_pending;
  passed &= expect(
      sample.coupler
                  .inspect_corrected_pending(published_two, pending,
                                             cleared_pending)
                  .code == StatusCode::invalid_plan,
      "cleared corrected-pending state cannot be inspected through an old certificate");
  return passed;
}

std::vector<double> full_refresh_cells(ConstFieldView field) {
  const Int3 cells = field.interior;
  std::vector<double> result(static_cast<std::size_t>(cells.x) * cells.y *
                             cells.z);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        result[full_refresh_offset(cell, cells)] = field.unchecked(cell, 0U);
      }
  return result;
}

std::vector<double> full_refresh_continuity(const Fixture& fixture,
                                            const FullRefreshSample& sample) {
  const Int3 cells = fixture.patch.cells;
  std::vector<double> result(static_cast<std::size_t>(cells.x) * cells.y *
                             cells.z);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        double value = kFullRefreshBdf.a0 * cell_volume(fixture, cell) *
                       sample.density.view.unchecked(cell, 0U);
        for (CartesianAxis axis :
             {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
          Int3 plus = cell;
          if (axis == CartesianAxis::x)
            ++plus.x;
          else if (axis == CartesianAxis::y)
            ++plus.y;
          else
            ++plus.z;
          const ConstFaceFieldView face =
              full_refresh_face(as_const(sample.total_flux), axis);
          value += face.unchecked(plus) - face.unchecked(cell);
        }
        result[full_refresh_offset(cell, cells)] = value;
      }
  return result;
}

std::vector<double>
full_refresh_central_difference(const std::vector<double>& plus,
                                const std::vector<double>& minus,
                                double epsilon) {
  std::vector<double> result(plus.size());
  for (std::size_t index = 0U; index < result.size(); ++index)
    result[index] = (plus[index] - minus[index]) / (2.0 * epsilon);
  return result;
}

std::vector<double> full_refresh_subtract(const std::vector<double>& left,
                                          const std::vector<double>& right) {
  std::vector<double> result(left.size());
  for (std::size_t index = 0U; index < result.size(); ++index)
    result[index] = left[index] - right[index];
  return result;
}

std::vector<double> full_refresh_add(const std::vector<double>& left,
                                     const std::vector<double>& right) {
  std::vector<double> result(left.size());
  for (std::size_t index = 0U; index < result.size(); ++index)
    result[index] = left[index] + right[index];
  return result;
}

double full_refresh_norm(const std::vector<double>& values) {
  long double sum = 0.0L;
  for (double value : values)
    sum += value * value;
  return std::sqrt(static_cast<double>(sum));
}

double full_refresh_relative_error(const std::vector<double>& actual,
                                   const std::vector<double>& expected) {
  return full_refresh_norm(full_refresh_subtract(actual, expected)) /
         std::max(full_refresh_norm(expected), 1.0e-30);
}

struct FullRefreshMissingResponse {
  OwnedFaceField x;
  OwnedFaceField y;
  OwnedFaceField z;
  std::vector<double> continuity;
  std::vector<double> energy;
};

FullRefreshMissingResponse full_refresh_missing_response(
    const Fixture& fixture, const FullRefreshSample& baseline,
    ConstFieldView state_direction, bool pressure_direction,
    const FrozenConvectionFaceField& frozen_enthalpy,
    StorageIdentity storage_seed) {
  const Int3 cells = fixture.patch.cells;
  FullRefreshMissingResponse response{
      make_face(CartesianAxis::x, cells, storage_seed),
      make_face(CartesianAxis::y, cells, storage_seed + 1U),
      make_face(CartesianAxis::z, cells, storage_seed + 2U),
      std::vector<double>(static_cast<std::size_t>(cells.x) * cells.y *
                          cells.z),
      std::vector<double>(static_cast<std::size_t>(cells.x) * cells.y *
                          cells.z)};
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    FaceFieldView output =
        axis == CartesianAxis::x
            ? response.x.view
            : (axis == CartesianAxis::y ? response.y.view : response.z.view);
    const std::uint8_t component = static_cast<std::uint8_t>(axis);
    for (std::int32_t z = 0; z < output.extents.z; ++z)
      for (std::int32_t y = 0; y < output.extents.y; ++y)
        for (std::int32_t x = 0; x < output.extents.x; ++x) {
          const Int3 face{x, y, z};
          Int3 left = face;
          if (axis == CartesianAxis::x)
            --left.x;
          else if (axis == CartesianAxis::y)
            --left.y;
          else
            --left.z;
          const Int3 wrapped_left = full_refresh_wrap(left, cells);
          const Int3 wrapped_right = full_refresh_wrap(face, cells);
          const auto density_variation = [&](Int3 cell) {
            const double derivative =
                pressure_direction
                    ? baseline.pressure_compressibility.view.unchecked(cell, 0U)
                    : baseline.enthalpy_compressibility.view.unchecked(cell,
                                                                       0U);
            return derivative * state_direction.unchecked(cell, 0U);
          };
          const double left_drho = density_variation(wrapped_left);
          const double right_drho = density_variation(wrapped_right);
          const double left_h =
              baseline.h_by_a.view.unchecked(wrapped_left, component);
          const double right_h =
              baseline.h_by_a.view.unchecked(wrapped_right, component);
          const double left_r =
              baseline.r_au.view.unchecked(wrapped_left, component);
          const double right_r =
              baseline.r_au.view.unchecked(wrapped_right, component);
          const double area = face_area(fixture, axis, face);
          const double b = -full_refresh_temporal_face(axis, face, cells);
          output.unchecked(face) =
              area * 0.5 * (left_drho * left_h + right_drho * right_h) +
              kFullRefreshBdf.a0 * b * 0.5 *
                  (left_drho * left_r + right_drho * right_r);
        }
  }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        double continuity = 0.0;
        double energy = 0.0;
        for (CartesianAxis axis :
             {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
          Int3 plus = cell;
          if (axis == CartesianAxis::x)
            ++plus.x;
          else if (axis == CartesianAxis::y)
            ++plus.y;
          else
            ++plus.z;
          const ConstFaceFieldView missing =
              axis == CartesianAxis::x
                  ? as_const(response.x.view)
                  : (axis == CartesianAxis::y ? as_const(response.y.view)
                                              : as_const(response.z.view));
          const ConstFaceFieldView h_face =
              full_refresh_face(frozen_enthalpy, axis);
          continuity += missing.unchecked(plus) - missing.unchecked(cell);
          energy += h_face.unchecked(plus) * missing.unchecked(plus) -
                    h_face.unchecked(cell) * missing.unchecked(cell);
        }
        const std::size_t offset = full_refresh_offset(cell, cells);
        response.continuity[offset] = continuity;
        response.energy[offset] = energy;
      }
  return response;
}

bool test_c1_full_refresh_four_block_sensitivity_red() {
  Fixture fixture;
  bool passed = expect(make_fixture(fixture, false, true),
                       "full-refresh periodic fixture compiles");
  passed &=
      expect(fixture.schemes.momentum() == ConvectionScheme::central2 &&
                 fixture.schemes.enthalpy() == ConvectionScheme::central2 &&
                 fixture.schemes.passive_scalar() == ConvectionScheme::tvd2 &&
                 fixture.boundary.transported_fields().size == 0U,
             "full-refresh physics stays central2 while an unregistered "
             "passive scalar only lifts the E_h boundary reach");
  if (!passed)
    return false;
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, fixture.equations, valid_spec(), piso)),
                   "full-refresh PISO plan compiles");
  if (!passed)
    return false;
  const Int3 cells = fixture.patch.cells;

  CommittedFluxHistory history;
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(initialize_flux_history(cells, history) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        9201U, history, previous_flux) &&
                       commit_zero_flux(fixture.equations.kernels(), cells,
                                        9202U, history, accepted_flux),
                   "full-refresh BDF2 history publishes two zero-flux layers");
  if (!passed)
    return false;

  OwnedField dp =
      make_field(kFullRefreshCorrectionField, cells, 1U, 2U, 9203U, 10203U);
  OwnedField dh = make_field(kFullRefreshEnthalpyDirectionField, cells, 1U, 2U,
                             9204U, 10204U);
  fill(dp, 0.0);
  fill(dh, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double xp = full_refresh_phase(x, cells.x, false);
        const double yp = full_refresh_phase(y, cells.y, false);
        const double zp = full_refresh_phase(z, cells.z, false);
        dp.view.unchecked(cell, 0U) =
            1.0e4 * (std::sin(xp) + 0.35 * std::cos(yp) - 0.20 * std::sin(zp));
        dh.view.unchecked(cell, 0U) =
            1.0e5 *
            (0.60 * std::cos(xp) - 0.40 * std::sin(yp) + 0.30 * std::cos(zp));
      }

  FullRefreshSample baseline;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(dp.view), as_const(dh.view), 0.0, 0.0, 9300U,
                       false, baseline),
                   "baseline EOS/transport/halo/boundary/C1 refresh closes");
  if (!passed)
    return false;

  long double gauge_numerator = 0.0L;
  long double gauge_weight = 0.0L;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double weight =
            cell_volume(fixture, cell) *
            baseline.pressure_compressibility.view.unchecked(cell, 0U);
        gauge_numerator += weight * dp.view.unchecked(cell, 0U);
        gauge_weight += weight;
      }
  const double gauge_shift =
      static_cast<double>(gauge_numerator / gauge_weight);
  long double post_gauge = 0.0L;
  long double post_gauge_scale = 0.0L;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        dp.view.unchecked(cell, 0U) -= gauge_shift;
        const double weighted =
            cell_volume(fixture, cell) *
            baseline.pressure_compressibility.view.unchecked(cell, 0U) *
            dp.view.unchecked(cell, 0U);
        post_gauge += weighted;
        post_gauge_scale += std::abs(weighted);
      }
  passed &= expect(
      std::abs(static_cast<double>(post_gauge)) <=
          2.0e-15 * static_cast<double>(post_gauge_scale),
      "pressure direction satisfies the closed-gauge weighted zero mean");

  OwnedFaceField frozen_x = make_face(CartesianAxis::x, cells, 10310U);
  OwnedFaceField frozen_y = make_face(CartesianAxis::y, cells, 10310U);
  OwnedFaceField frozen_z = make_face(CartesianAxis::z, cells, 10310U);
  const FrozenConvectionContext convection_context{
      fixture.equations.semantic_fingerprint(), fixture.boundary.revision()};
  FrozenConvectionFaceField frozen_enthalpy;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          fixture.kernels, fixture.schemes.enthalpy(),
          baseline.intermediate_flux, as_const(baseline.enthalpy.view), 0U,
          convection_context, {frozen_x.view, frozen_y.view, frozen_z.view},
          frozen_enthalpy)) &&
          frozen_enthalpy.valid(),
      "baseline target enthalpy faces freeze from the real C1 flux");
  if (!passed)
    return false;

  OwnedField cp_response = make_field(140U, cells, 1U, 0U, 9310U, 10320U);
  OwnedField ch_response = make_field(141U, cells, 1U, 0U, 9311U, 10321U);
  OwnedField ep_response = make_field(142U, cells, 1U, 0U, 9312U, 10322U);
  OwnedField eh_response = make_field(143U, cells, 1U, 0U, 9313U, 10323U);
  OwnedField energy_temporal = make_field(144U, cells, 1U, 0U, 9314U, 10324U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double volume = cell_volume(fixture, cell);
        ch_response.view.unchecked(cell, 0U) =
            kFullRefreshBdf.a0 * volume *
            baseline.enthalpy_compressibility.view.unchecked(cell, 0U) *
            dh.view.unchecked(cell, 0U);
        energy_temporal.view.unchecked(cell, 0U) =
            kFullRefreshBdf.a0 * volume *
            (baseline.enthalpy.view.unchecked(cell, 0U) *
                 baseline.pressure_compressibility.view.unchecked(cell, 0U) -
             1.0);
      }

  const PressureCorrectionSystemView pressure_system{
      baseline.pressure_diagonal.view, baseline.pressure_rhs.view};
  PressureLinearOperator continuity_pressure;
  const LinearIdentity identity{9315U, 9316U, 9317U, 9318U, 9319U};
  passed &=
      expect(static_cast<bool>(baseline.coupler.bind_pressure_operator(
                 {MPI_COMM_SELF, &baseline.correction_halo, 313U,
                  kFullRefreshCorrectionField},
                 pressure_system, continuity_pressure)) &&
                 static_cast<bool>(continuity_pressure.refresh(
                     {baseline.pressure_certificate, identity, 9320U})) &&
                 static_cast<bool>(
                     continuity_pressure.apply(dp.view, cp_response.view)),
             "current frozen C_p applies through the real pressure operator");

  PisoCartesianPressureWorkLinearization pressure_work;
  passed &= expect(
      static_cast<bool>(
          baseline.coupler.inspect_cartesian_pressure_work_linearization(
              baseline.intermediate, baseline.pressure_certificate,
              as_const(baseline.pressure.view),
              as_const(baseline.predictor_velocity.view), pressure_work)) &&
          pressure_work.valid(),
      "baseline C1 issues exact Cartesian pressure-work authority");
  PressureEnergyFrozenFaceEnthalpy pressure_energy_frozen{
      frozen_enthalpy.x,
      frozen_enthalpy.y,
      frozen_enthalpy.z,
      frozen_enthalpy.revision,
      frozen_enthalpy.reconstruction,
      0U};
  pressure_energy_frozen.local_binding =
      pressure_energy_frozen_face_enthalpy_local_binding(
          pressure_energy_frozen);
  PressureEnergyPressureFluxBinding ep_binding;
  ep_binding.geometry = &fixture.geometry;
  ep_binding.boundary = &fixture.boundary;
  ep_binding.patch = fixture.patch;
  ep_binding.services = {MPI_COMM_SELF, &baseline.correction_halo, 313U,
                         kFullRefreshCorrectionField, 1U};
  ep_binding.intermediate = baseline.intermediate;
  ep_binding.pressure = baseline.pressure_certificate;
  ep_binding.temporal_diagonal = as_const(energy_temporal.view);
  ep_binding.x_pressure_coefficient = as_const(baseline.pressure_x.view);
  ep_binding.y_pressure_coefficient = as_const(baseline.pressure_y.view);
  ep_binding.z_pressure_coefficient = as_const(baseline.pressure_z.view);
  ep_binding.target_flux = baseline.intermediate_flux;
  ep_binding.frozen_face_enthalpy = pressure_energy_frozen;
  ep_binding.identity = identity;
  ep_binding.pressure_work = pressure_work;
  PressureEnergyPressureFluxOperator ep_operator;
  PressureEnergyPressureFluxCertificate ep_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          ep_binding, ep_operator, ep_certificate)) &&
          static_cast<bool>(ep_operator.apply(dp.view, ep_response.view)),
      "current frozen E_p applies through the complete pressure-work block");

  OwnedField delta_temperature = make_field(
      kFullRefreshTemperatureDirectionField, cells, 1U, 1U, 9321U, 10331U);
  OwnedFaceField directional_x = make_face(CartesianAxis::x, cells, 10332U);
  OwnedFaceField directional_y = make_face(CartesianAxis::y, cells, 10332U);
  OwnedFaceField directional_z = make_face(CartesianAxis::z, cells, 10332U);
  const std::array<HaloFieldSpec, 2U> energy_specs{{
      {kFullRefreshEnthalpyDirectionField, 2U, 1U},
      {kFullRefreshTemperatureDirectionField, 1U, 1U},
  }};
  HaloEngine energy_halo;
  passed &= expect(static_cast<bool>(energy_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {energy_specs.data(), energy_specs.size()},
                       fixture.boundary.halo_topology())),
                   "current frozen E_h variation halo reserves");
  PressureEnergyEnthalpyBinding eh_binding;
  eh_binding.geometry = &fixture.geometry;
  eh_binding.kernels = &fixture.kernels;
  eh_binding.boundary = &fixture.boundary;
  eh_binding.patch = fixture.patch;
  eh_binding.convection = fixture.schemes.enthalpy();
  eh_binding.services = {MPI_COMM_SELF, &energy_halo, 314U,
                         kFullRefreshEnthalpyDirectionField,
                         kFullRefreshTemperatureDirectionField};
  eh_binding.authority = {kFullRefreshBdf,
                          kFullRefreshTime,
                          fixture.geometry.topology_revision(),
                          fixture.boundary.revision(),
                          fixture.thermodynamics.fingerprint(),
                          fixture.transport.fingerprint(),
                          fixture.equations.semantic_fingerprint(),
                          fixture.thermodynamics.fingerprint(),
                          fixture.transport.fingerprint()};
  eh_binding.assembled_diagonal = as_const(baseline.energy_diagonal.view);
  eh_binding.target_enthalpy = as_const(baseline.enthalpy.view);
  eh_binding.density_enthalpy_derivative =
      as_const(baseline.enthalpy_compressibility.view);
  eh_binding.heat_capacity = as_const(baseline.heat_capacity.view);
  eh_binding.thermal_conductivity = as_const(baseline.conductivity.view);
  eh_binding.enthalpy_diffusivity =
      as_const(baseline.enthalpy_diffusivity.view);
  eh_binding.target_flux = baseline.intermediate_flux;
  eh_binding.convection_context = convection_context;
  eh_binding.frozen_face_enthalpy = frozen_enthalpy;
  eh_binding.workspace = {
      delta_temperature.view,
      {directional_x.view, directional_y.view, directional_z.view}};
  eh_binding.identity = identity;
  eh_binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope;
  PressureEnergyEnthalpyOperator eh_operator;
  PressureEnergyEnthalpyCertificate eh_certificate;
  const Status eh_bind = PressureEnergyEnthalpyOperator::bind(
      eh_binding, eh_operator, eh_certificate);
  const Status eh_apply =
      eh_bind ? eh_operator.apply(dh.view, eh_response.view) : eh_bind;
  if (!eh_bind || !eh_apply) {
    std::cerr << "full-refresh E_h bind/apply status="
              << static_cast<unsigned>(eh_bind.code) << '/' << eh_bind.detail
              << " -> " << static_cast<unsigned>(eh_apply.code) << '/'
              << eh_apply.detail << '\n';
  }
  passed &= expect(static_cast<bool>(eh_bind) && static_cast<bool>(eh_apply),
                   "current frozen E_h applies through the real spatial block");
  if (!passed)
    return false;

  const std::vector<double> cp = full_refresh_cells(as_const(cp_response.view));
  const std::vector<double> ch = full_refresh_cells(as_const(ch_response.view));
  const std::vector<double> ep = full_refresh_cells(as_const(ep_response.view));
  const std::vector<double> eh = full_refresh_cells(as_const(eh_response.view));
  const FullRefreshMissingResponse missing_p = full_refresh_missing_response(
      fixture, baseline, as_const(dp.view), true, frozen_enthalpy, 10401U);
  const FullRefreshMissingResponse missing_h = full_refresh_missing_response(
      fixture, baseline, as_const(dh.view), false, frozen_enthalpy, 10411U);
  std::cout << "full-refresh missing norms C_p/C_h/E_p/E_h="
            << full_refresh_norm(missing_p.continuity) << '/'
            << full_refresh_norm(missing_h.continuity) << '/'
            << full_refresh_norm(missing_p.energy) << '/'
            << full_refresh_norm(missing_h.energy) << '\n';
  passed &=
      expect(full_refresh_norm(missing_p.continuity) > 1.0e-6 &&
                 full_refresh_norm(missing_h.continuity) > 1.0e-6 &&
                 full_refresh_norm(missing_p.energy) > 1.0e-3 &&
                 full_refresh_norm(missing_h.energy) > 1.0e-3,
             "independent D(dphi*) and D(hhat*dphi*) oracles are nonzero");

  struct LadderRow {
    double epsilon{};
    double cp_error{};
    double ch_error{};
    double ep_error{};
    double eh_error{};
    std::vector<double> c_p_fd;
    std::vector<double> c_h_fd;
    std::vector<double> e_p_fd;
    std::vector<double> e_h_fd;
  };
  constexpr std::array<double, 6U> epsilons{8.0e-4, 4.0e-4, 2.0e-4,
                                            1.0e-4, 5.0e-5, 2.5e-5};
  std::array<LadderRow, epsilons.size()> ladder{};
  for (std::size_t level = 0U; level < epsilons.size(); ++level) {
    const double epsilon = epsilons[level];
    FullRefreshSample p_plus;
    FullRefreshSample p_minus;
    FullRefreshSample h_plus;
    FullRefreshSample h_minus;
    const RevisionToken level_seed = 11000U + 1000U * level;
    const bool built =
        build_full_refresh_sample(fixture, piso, history, accepted_flux,
                                  previous_flux, as_const(dp.view),
                                  as_const(dh.view), epsilon, 0.0,
                                  level_seed + 100U, false, p_plus) &&
        build_full_refresh_sample(fixture, piso, history, accepted_flux,
                                  previous_flux, as_const(dp.view),
                                  as_const(dh.view), -epsilon, 0.0,
                                  level_seed + 200U, false, p_minus) &&
        build_full_refresh_sample(fixture, piso, history, accepted_flux,
                                  previous_flux, as_const(dp.view),
                                  as_const(dh.view), 0.0, epsilon,
                                  level_seed + 300U, false, h_plus) &&
        build_full_refresh_sample(fixture, piso, history, accepted_flux,
                                  previous_flux, as_const(dp.view),
                                  as_const(dh.view), 0.0, -epsilon,
                                  level_seed + 400U, false, h_minus);
    passed &= expect(built, "epsilon-ladder EOS/transport/C1 samples rebuild");
    if (!built)
      return false;
    LadderRow& row = ladder[level];
    row.epsilon = epsilon;
    row.c_p_fd = full_refresh_central_difference(
        full_refresh_continuity(fixture, p_plus),
        full_refresh_continuity(fixture, p_minus), epsilon);
    row.c_h_fd = full_refresh_central_difference(
        full_refresh_continuity(fixture, h_plus),
        full_refresh_continuity(fixture, h_minus), epsilon);
    row.e_p_fd = full_refresh_central_difference(
        full_refresh_cells(as_const(p_plus.energy_residual.view)),
        full_refresh_cells(as_const(p_minus.energy_residual.view)), epsilon);
    row.e_h_fd = full_refresh_central_difference(
        full_refresh_cells(as_const(h_plus.energy_residual.view)),
        full_refresh_cells(as_const(h_minus.energy_residual.view)), epsilon);
    row.cp_error = full_refresh_relative_error(
        full_refresh_subtract(row.c_p_fd, cp), missing_p.continuity);
    row.ch_error = full_refresh_relative_error(
        full_refresh_subtract(row.c_h_fd, ch), missing_h.continuity);
    row.ep_error = full_refresh_relative_error(
        full_refresh_subtract(row.e_p_fd, ep), missing_p.energy);
    row.eh_error = full_refresh_relative_error(
        full_refresh_subtract(row.e_h_fd, eh), missing_h.energy);
    std::cout << "full-refresh epsilon=" << epsilon
              << " mismatch-relative-error Cp/Ch/Ep/Eh=" << row.cp_error << '/'
              << row.ch_error << '/' << row.ep_error << '/' << row.eh_error
              << '\n';
    if (level + 1U == epsilons.size()) {
      std::cout << "full-refresh finest norms fd Cp/Ch/Ep/Eh="
                << full_refresh_norm(row.c_p_fd) << '/'
                << full_refresh_norm(row.c_h_fd) << '/'
                << full_refresh_norm(row.e_p_fd) << '/'
                << full_refresh_norm(row.e_h_fd)
                << " frozen=" << full_refresh_norm(cp) << '/'
                << full_refresh_norm(ch) << '/' << full_refresh_norm(ep) << '/'
                << full_refresh_norm(eh) << " mismatch="
                << full_refresh_norm(full_refresh_subtract(row.c_p_fd, cp))
                << '/'
                << full_refresh_norm(full_refresh_subtract(row.c_h_fd, ch))
                << '/'
                << full_refresh_norm(full_refresh_subtract(row.e_p_fd, ep))
                << '/'
                << full_refresh_norm(full_refresh_subtract(row.e_h_fd, eh))
                << '\n';
    }
  }

  const LadderRow& finest = ladder.back();
  passed &= expect(
      finest.cp_error < 1.0e-6 && finest.ch_error < 1.0e-6 &&
          finest.ep_error < 1.0e-6 && finest.eh_error < 1.0e-6,
      "C1 four-block mismatch equals the independent omitted-flux response");
  const auto second_order_before_roundoff = [&](auto member) {
    bool observed = false;
    for (std::size_t level = 0U; level + 1U < ladder.size(); ++level) {
      const double coarse = ladder[level].*member;
      const double fine = ladder[level + 1U].*member;
      if (coarse > 5.0e-10) {
        const double ratio = fine / coarse;
        observed = observed || (ratio > 0.12 && ratio < 0.42);
      }
    }
    return observed;
  };
  passed &= expect(second_order_before_roundoff(&LadderRow::cp_error) &&
                       second_order_before_roundoff(&LadderRow::ch_error) &&
                       second_order_before_roundoff(&LadderRow::ep_error) &&
                       second_order_before_roundoff(&LadderRow::eh_error),
                   "central full-refresh FD remainder quarters on epsilon "
                   "halving before roundoff");

  const double mixed_epsilon = epsilons.back();
  FullRefreshSample mixed_plus;
  FullRefreshSample mixed_minus;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(dp.view), as_const(dh.view), mixed_epsilon,
                       mixed_epsilon, 16000U, false, mixed_plus) &&
                       build_full_refresh_sample(
                           fixture, piso, history, accepted_flux, previous_flux,
                           as_const(dp.view), as_const(dh.view), -mixed_epsilon,
                           -mixed_epsilon, 17000U, false, mixed_minus),
                   "mixed pressure-enthalpy full-refresh samples rebuild");
  const std::vector<double> mixed_c = full_refresh_central_difference(
      full_refresh_continuity(fixture, mixed_plus),
      full_refresh_continuity(fixture, mixed_minus), mixed_epsilon);
  const std::vector<double> mixed_e = full_refresh_central_difference(
      full_refresh_cells(as_const(mixed_plus.energy_residual.view)),
      full_refresh_cells(as_const(mixed_minus.energy_residual.view)),
      mixed_epsilon);
  const double mixed_c_error = full_refresh_relative_error(
      full_refresh_subtract(mixed_c, full_refresh_add(cp, ch)),
      full_refresh_add(missing_p.continuity, missing_h.continuity));
  const double mixed_e_error = full_refresh_relative_error(
      full_refresh_subtract(mixed_e, full_refresh_add(ep, eh)),
      full_refresh_add(missing_p.energy, missing_h.energy));
  std::cout << "full-refresh mixed mismatch-relative-error C/E="
            << mixed_c_error << '/' << mixed_e_error << '\n';
  passed &= expect(
      mixed_c_error < 2.0e-6 && mixed_e_error < 2.0e-6,
      "mixed direction closes as the sum of the four certified block gaps");

  FullRefreshSample c2_plus;
  FullRefreshSample c2_minus;
  passed &= expect(build_full_refresh_sample(
                       fixture, piso, history, accepted_flux, previous_flux,
                       as_const(dp.view), as_const(dh.view), mixed_epsilon,
                       mixed_epsilon, 18000U, true, c2_plus) &&
                       build_full_refresh_sample(
                           fixture, piso, history, accepted_flux, previous_flux,
                           as_const(dp.view), as_const(dh.view), -mixed_epsilon,
                           -mixed_epsilon, 19000U, true, c2_minus),
                   "C2 byte-copy negative-control samples refresh");
  const std::vector<double> c2_plus_flux =
      face_flux_values(c2_plus.intermediate_flux);
  const std::vector<double> c2_minus_flux =
      face_flux_values(c2_minus.intermediate_flux);
  const double c2_fd_norm = full_refresh_norm(full_refresh_central_difference(
      c2_plus_flux, c2_minus_flux, mixed_epsilon));
  std::cout << "full-refresh C2 byte-copy equal="
            << (c2_plus_flux == c2_minus_flux) << " face-FD-norm=" << c2_fd_norm
            << '\n';
  passed &= expect(c2_plus_flux == c2_minus_flux &&
                       c2_plus_flux ==
                           face_flux_values(as_const(c2_plus.c2_trial_flux)) &&
                       c2_fd_norm == 0.0 &&
                       full_refresh_norm(missing_p.continuity) +
                               full_refresh_norm(missing_h.continuity) >
                           1.0e-6,
                   "C2 byte-copy is an explicit false-negative control for the "
                   "nonzero C1 flux sensitivity");
  return passed;
}

bool test_candidate_boundary_compiler_fixture() {
  bool passed = true;
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.inlet = CandidateBoundaryInlet::velocity;
  spec.inlet_velocity = 1.25;
  passed &= expect(fixture.initialize(MPI_COMM_SELF, spec),
                   "real open-boundary compiler fixture initializes");
  if (fixture.diagnostic_step != 0U)
    std::cerr << "candidate-boundary initialize step="
              << fixture.diagnostic_step << " status="
              << static_cast<unsigned>(fixture.diagnostic_status.code) << '/'
              << fixture.diagnostic_status.detail << '\n';
  CandidateBoundaryScratch baseline;
  const bool staged = fixture.stage(0.0, 16.0, 4.0, 24000U, baseline);
  passed &= expect(
      staged && baseline.final_boundary.valid(),
      "real compiler chain finalizes an alpha-zero velocity-inlet baseline");

  const Int3 cells = fixture.patch.cells;
  if (staged) {
    double inlet_h = 0.0;
    double inlet_cp = 0.0;
    double inlet_gas = 0.0;
    ThermoState inlet_state;
    const Int3 owner{0, 0, 0};
    passed &= expect(
        static_cast<bool>(fixture.thermodynamics.mixture_enthalpy(
            spec.inlet_temperature, {}, inlet_h, inlet_cp, inlet_gas)) &&
            static_cast<bool>(fixture.thermodynamics.evaluate(
                CandidateBoundaryFixture::absolute_pressure_reference +
                    baseline.pressure.view.unchecked(owner, 0U),
                inlet_h, {}, {spec.inlet_velocity, 0.0, 0.0},
                inlet_state, spec.inlet_temperature)),
        "velocity-inlet independent EOS oracle evaluates");
    const double expected_inlet =
        inlet_state.rho * spec.inlet_velocity *
        fixture.face_area(CartesianAxis::x, {0, 0, 0});
    passed &= expect(
        std::abs(baseline.final_flux.x.unchecked({0, 0, 0}) -
                 expected_inlet) < 1.0e-12 &&
            same_double_bits(
                baseline.final_flux.x.unchecked({cells.x, 0, 0}),
                baseline.mechanical_flux.x.unchecked({cells.x, 0, 0})) &&
            same_double_bits(
                baseline.final_flux.x.unchecked({1, 0, 0}),
                baseline.mechanical_flux.x.unchecked({1, 0, 0})) &&
            same_double_bits(
                baseline.final_flux.y.unchecked({0, 0, 0}), 0.0) &&
            same_double_bits(
                baseline.final_flux.z.unchecked({0, 0, 0}), 0.0),
        "velocity inlet uses same-target EOS rho while outlet/internal remain mechanical and symmetry is +0");

    PressureEnergyCandidateBoundaryFinalizeInput bad_alias =
        fixture.finalizer_input(baseline);
    ++bad_alias.thermophysical_boundary.binding.pressure_perturbation.field;
    const std::vector<double> before_alias_rejection =
        face_flux_values(as_const(baseline.final_flux));
    FinalBoundaryFluxCertificate cleared = baseline.final_boundary;
    passed &= expect(
        fixture.finalizer
                    .finalize(bad_alias, fixture.reductions, cleared)
                    .code == StatusCode::invalid_plan &&
            !cleared.valid() &&
            face_flux_values(as_const(baseline.final_flux)) ==
                before_alias_rejection,
        "foreign candidate thermophysical semantic alias clears a reused certificate and zero-writes flux");

    Status hot_status;
    {
      candidate_boundary_allocation_observer::Guard guard;
      for (std::uint32_t repetition = 0U; repetition < 16U;
           ++repetition) {
        FinalBoundaryFluxCertificate replayed;
        const Status status = fixture.finalizer.finalize(
            fixture.finalizer_input(baseline), fixture.reductions,
            replayed);
        if (hot_status && (!status || !replayed.valid()))
          hot_status = status ? Status{StatusCode::invalid_plan, 0U}
                              : status;
      }
    }
    passed &= expect(
        static_cast<bool>(hot_status) &&
            candidate_boundary_allocation_observer::count.load(
                std::memory_order_relaxed) == 0U,
        "prebound final-boundary prepare/replay/commit path performs zero heap allocations");

    PressureEnergyCandidateBoundaryFinalizeInput revised_output =
        fixture.finalizer_input(baseline);
    ++revised_output.final_flux.revision;
    FinalBoundaryFluxCertificate replayed_boundary;
    passed &= expect(
        static_cast<bool>(fixture.finalizer.finalize(
            revised_output, fixture.reductions, replayed_boundary)) &&
            replayed_boundary.valid() &&
            replayed_boundary.canonical_lineage() ==
                baseline.final_boundary.canonical_lineage() &&
            replayed_boundary.scratch_binding() !=
                baseline.final_boundary.scratch_binding(),
        "same final physical flux under a new output revision preserves canonical provenance but changes scratch binding");
  }

  CandidateBoundaryFixture multispecies_fixture;
  CandidateBoundaryFixtureSpec multispecies_spec;
  multispecies_spec.multispecies = true;
  multispecies_spec.allow_backflow = true;
  CandidateBoundaryScratch multispecies_baseline;
  CandidateBoundaryScratch multispecies_positive;
  CandidateBoundaryScratch multispecies_positive_replay;
  CandidateBoundaryScratch multispecies_backflow;
  PisoFrozenMomentumExactCandidateCertificate multispecies_exact_baseline;
  PisoFrozenMomentumExactCandidateCertificate multispecies_exact_positive;
  PisoFrozenMomentumExactCandidateCertificate
      multispecies_exact_positive_replay;
  const bool multispecies_initialized =
      multispecies_fixture.initialize(MPI_COMM_SELF, multispecies_spec);
  const bool multispecies_baseline_staged =
      multispecies_initialized &&
      multispecies_fixture.stage(0.0, 16.0, 0.0, 24400U,
                                 multispecies_baseline);
  const auto aliased_flux = [&](test::OwnedField& target, RevisionToken revision,
                                FaceFluxStorage& storage,
                                FaceFluxView& output) {
    Status status = FaceFluxStorage::allocate_workspace(
        multispecies_fixture.patch.cells, 1U, storage);
    if (status) status = storage.workspace_view(0U, revision, output);
    if (!status) return status;
    fill_face_flux(output, -993.0);
    output.x.base = target.view.base;
    output.x.storage_identity = output.y.storage_identity;
    output.x.revision_domain = output.y.revision_domain;
    return Status{};
  };
  const auto rejects_hot_alias = [&](test::OwnedField& target,
                                     RevisionToken revision) {
    FaceFluxStorage storage;
    FaceFluxView output;
    Status status = aliased_flux(target, revision, storage, output);
    const std::vector<double> target_before = target.storage;
    const std::vector<double> output_before =
        status ? face_flux_values(as_const(output)) : std::vector<double>{};
    FinalBoundaryFluxCertificate certificate =
        multispecies_baseline.final_boundary;
    if (status) {
      PressureEnergyCandidateBoundaryFinalizeInput input =
          multispecies_fixture.finalizer_input(multispecies_baseline);
      input.final_flux = output;
      status = multispecies_fixture.finalizer.finalize(
          input, multispecies_fixture.reductions, certificate);
    }
    const bool rejected =
        status.code == StatusCode::invalid_plan && !certificate.valid() &&
        target.storage == target_before &&
        face_flux_values(as_const(output)) == output_before;
    std::copy(target_before.begin(), target_before.end(), target.storage.begin());
    return rejected;
  };
  const auto rejects_fresh_alias = [&](test::OwnedField& target,
                                       RevisionToken revision) {
    FaceFluxStorage storage;
    FaceFluxView output;
    Status status = aliased_flux(target, revision, storage, output);
    const std::vector<double> target_before = target.storage;
    const std::vector<double> output_before =
        status ? face_flux_values(as_const(output)) : std::vector<double>{};
    if (status)
      status = multispecies_fixture.finalizer.close_fresh_physical_flux(
          {CandidateBoundaryFixture::absolute_pressure_reference,
           multispecies_baseline.thermophysical_pressure_alias,
           as_const(multispecies_baseline.velocity.view),
           {multispecies_baseline.thermophysical_species_aliases.data(),
            multispecies_baseline.thermophysical_species_aliases.size()},
           as_const(multispecies_baseline.mechanical_flux), output},
          multispecies_fixture.reductions);
    const bool rejected =
        status.code == StatusCode::invalid_plan &&
        target.storage == target_before &&
        face_flux_values(as_const(output)) == output_before;
    std::copy(target_before.begin(), target_before.end(), target.storage.begin());
    return rejected;
  };
  if (multispecies_baseline_staged) {
    std::array<test::OwnedField*, 6U> hot_inputs{{
        &multispecies_baseline.pressure, &multispecies_baseline.enthalpy,
        &multispecies_baseline.density, &multispecies_baseline.temperature,
        &multispecies_baseline.velocity,
        &multispecies_baseline.independent_species[0U]}};
    bool hot_aliases_rejected = true;
    for (std::size_t input = 0U; input < hot_inputs.size(); ++input)
      hot_aliases_rejected =
          rejects_hot_alias(*hot_inputs[input], 24410U + input) &&
          hot_aliases_rejected;
    passed &= expect(
        hot_aliases_rejected,
        "hot final-boundary closure rejects final flux aliasing p/h/rho/T/U/Y without writes");

    std::array<test::OwnedField*, 3U> fresh_inputs{{
        &multispecies_baseline.pressure, &multispecies_baseline.velocity,
        &multispecies_baseline.independent_species[0U]}};
    bool fresh_aliases_rejected = true;
    for (std::size_t input = 0U; input < fresh_inputs.size(); ++input)
      fresh_aliases_rejected =
          rejects_fresh_alias(*fresh_inputs[input], 24420U + input) &&
          fresh_aliases_rejected;
    passed &= expect(
        fresh_aliases_rejected,
        "Fresh final-boundary closure rejects final flux aliasing p/U/Y without writes");
  }
  const Status multispecies_baseline_status =
      multispecies_baseline_staged
          ? multispecies_fixture.coupler
                .certify_frozen_momentum_exact_baseline(
                    multispecies_fixture.authority,
                    multispecies_baseline.pressure_stage,
                    multispecies_baseline.velocity_stage,
                    multispecies_baseline.flux_stage,
                    multispecies_fixture.exact_input(multispecies_baseline),
                    multispecies_fixture.reductions,
                    multispecies_exact_baseline)
          : Status{StatusCode::invalid_plan, 0U};
  const bool multispecies_positive_staged =
      multispecies_baseline_status &&
      multispecies_fixture.stage(0.5, 16.0, 0.0, 24500U,
                                 multispecies_positive);
  const Status multispecies_positive_status =
      multispecies_positive_staged
          ? multispecies_fixture.coupler
                .certify_frozen_momentum_exact_candidate(
                    multispecies_fixture.authority,
                    multispecies_exact_baseline,
                    multispecies_positive.pressure_stage,
                    multispecies_positive.velocity_stage,
                    multispecies_positive.flux_stage,
                    multispecies_fixture.exact_input(multispecies_positive),
                    multispecies_fixture.reductions,
                    multispecies_exact_positive)
          : Status{StatusCode::invalid_plan, 0U};
  const bool multispecies_positive_replay_staged =
      multispecies_positive_status &&
      multispecies_fixture.stage(0.5, 16.0, 0.0, 24550U,
                                 multispecies_positive_replay);
  const Status multispecies_positive_replay_status =
      multispecies_positive_replay_staged
          ? multispecies_fixture.coupler
                .certify_frozen_momentum_exact_candidate(
                    multispecies_fixture.authority,
                    multispecies_exact_baseline,
                    multispecies_positive_replay.pressure_stage,
                    multispecies_positive_replay.velocity_stage,
                    multispecies_positive_replay.flux_stage,
                    multispecies_fixture.exact_input(
                        multispecies_positive_replay),
                    multispecies_fixture.reductions,
                    multispecies_exact_positive_replay)
          : Status{StatusCode::invalid_plan, 0U};
  bool multispecies_no_revision_mutation_rejected = false;
  if (multispecies_positive_replay_status &&
      !multispecies_positive_replay.independent_species.empty()) {
    double& poisoned =
        multispecies_positive_replay.independent_species[0U]
            .view.unchecked({0, 0, 0}, 0U);
    const double original = poisoned;
    poisoned = std::nextafter(original,
                              std::numeric_limits<double>::infinity());
    PisoFrozenMomentumExactCandidateCertificate reused_after_poison =
        multispecies_exact_positive_replay;
    const Status poisoned_status =
        multispecies_fixture.coupler
            .certify_frozen_momentum_exact_candidate(
                multispecies_fixture.authority,
                multispecies_exact_baseline,
                multispecies_positive_replay.pressure_stage,
                multispecies_positive_replay.velocity_stage,
                multispecies_positive_replay.flux_stage,
                multispecies_fixture.exact_input(
                    multispecies_positive_replay),
                multispecies_fixture.reductions, reused_after_poison);
    multispecies_no_revision_mutation_rejected =
        poisoned_status.code == StatusCode::invalid_plan &&
        !reused_after_poison.valid();
    poisoned = original;
  }
  const bool multispecies_backflow_staged =
      multispecies_initialized &&
      multispecies_fixture.stage(1.0, -1000.0, 0.0, 24600U,
                                 multispecies_backflow);
  if (!(multispecies_baseline_status &&
        multispecies_exact_baseline.valid() &&
        multispecies_positive_status && multispecies_exact_positive.valid() &&
        multispecies_positive_replay_status &&
        multispecies_exact_positive_replay.valid() &&
        multispecies_backflow_staged))
    std::cerr << "candidate-boundary multispecies init/stages="
              << multispecies_initialized << '/'
              << multispecies_baseline_staged << '/'
              << static_cast<unsigned>(multispecies_baseline_status.code)
              << ':' << multispecies_baseline_status.detail << '/'
              << multispecies_positive_staged << '/'
              << static_cast<unsigned>(multispecies_positive_status.code)
              << ':' << multispecies_positive_status.detail << '/'
              << multispecies_backflow_staged << " fixture="
              << multispecies_fixture.diagnostic_step << '/'
              << static_cast<unsigned>(
                     multispecies_fixture.diagnostic_status.code)
              << ':' << multispecies_fixture.diagnostic_status.detail
              << '\n';
  const std::array<double, 2U> inlet_composition{{0.2, 0.3}};
  const std::array<double, 2U> backflow_composition{{0.25, 0.35}};
  double inlet_h = 0.0;
  double inlet_cp = 0.0;
  double inlet_gas = 0.0;
  double backflow_h = 0.0;
  double backflow_cp = 0.0;
  double backflow_gas = 0.0;
  ThermoState inlet_thermo;
  ThermoState backflow_thermo;
  ThermoState wrong_backflow_thermo;
  Status multispecies_oracle_status =
      multispecies_initialized
          ? multispecies_fixture.thermodynamics.mixture_enthalpy(
                multispecies_spec.inlet_temperature,
                {inlet_composition.data(), inlet_composition.size()},
                inlet_h, inlet_cp, inlet_gas)
          : Status{StatusCode::invalid_plan, 0U};
  if (multispecies_oracle_status && multispecies_positive_staged) {
    const Int3 owner{0, 0, 0};
    multispecies_oracle_status = multispecies_fixture.thermodynamics.evaluate(
        CandidateBoundaryFixture::absolute_pressure_reference +
            multispecies_positive.pressure.view.unchecked(owner, 0U),
        inlet_h, {inlet_composition.data(), inlet_composition.size()},
        {multispecies_spec.inlet_velocity, 0.0, 0.0}, inlet_thermo,
        multispecies_spec.inlet_temperature);
  }
  if (multispecies_oracle_status) {
    multispecies_oracle_status =
        multispecies_fixture.thermodynamics.mixture_enthalpy(
            multispecies_spec.backflow_temperature,
            {backflow_composition.data(), backflow_composition.size()},
            backflow_h, backflow_cp, backflow_gas);
  }
  if (multispecies_oracle_status) {
    multispecies_oracle_status = multispecies_fixture.thermodynamics.evaluate(
        multispecies_spec.outlet_pressure, backflow_h,
        {backflow_composition.data(), backflow_composition.size()},
        {multispecies_spec.backflow_velocity, 0.0, 0.0}, backflow_thermo,
        multispecies_spec.backflow_temperature);
  }
  if (multispecies_oracle_status) {
    double wrong_h = 0.0;
    double wrong_cp = 0.0;
    double wrong_gas = 0.0;
    multispecies_oracle_status =
        multispecies_fixture.thermodynamics.mixture_enthalpy(
            multispecies_spec.backflow_temperature,
            {inlet_composition.data(), inlet_composition.size()}, wrong_h,
            wrong_cp, wrong_gas);
    if (multispecies_oracle_status)
      multispecies_oracle_status =
          multispecies_fixture.thermodynamics.evaluate(
              multispecies_spec.outlet_pressure, wrong_h,
              {inlet_composition.data(), inlet_composition.size()},
              {multispecies_spec.backflow_velocity, 0.0, 0.0},
              wrong_backflow_thermo,
              multispecies_spec.backflow_temperature);
  }
  const double inlet_area =
      multispecies_initialized
          ? multispecies_fixture.face_area(CartesianAxis::x, {0, 0, 0})
          : 0.0;
  const Int3 outlet_face{multispecies_fixture.patch.cells.x, 0, 0};
  const Int3 outlet_owner{multispecies_fixture.patch.cells.x - 1, 0, 0};
  const double outlet_area =
      multispecies_initialized
          ? multispecies_fixture.face_area(CartesianAxis::x, outlet_face)
          : 0.0;
  const double expected_inlet_flux =
      inlet_thermo.rho * multispecies_spec.inlet_velocity * inlet_area;
  const double expected_backflow_flux =
      backflow_thermo.rho * multispecies_spec.backflow_velocity * outlet_area;
  const double wrong_backflow_flux = wrong_backflow_thermo.rho *
                                     multispecies_spec.backflow_velocity *
                                     outlet_area;
  const auto close_flux = [](double actual, double expected) noexcept {
    return std::abs(actual - expected) <=
           64.0 * std::numeric_limits<double>::epsilon() *
               std::max({1.0, std::abs(actual), std::abs(expected)});
  };
  const bool inlet_ghost_targets =
      multispecies_positive_staged &&
      multispecies_positive.thermophysical_species_aliases.size() == 2U &&
      multispecies_positive.thermophysical_species_aliases[0U].unchecked(
          {-1, 0, 0}, 0U) == inlet_composition[0U] &&
      multispecies_positive.thermophysical_species_aliases[1U].unchecked(
          {-1, 0, 0}, 0U) == inlet_composition[1U];
  const bool backflow_face_targets =
      multispecies_backflow_staged &&
      multispecies_backflow.thermophysical_species_aliases.size() == 2U &&
      0.5 * (multispecies_backflow.independent_species[0U].view.unchecked(
                 outlet_owner, 0U) +
             multispecies_backflow.thermophysical_species_aliases[0U]
                 .unchecked(outlet_face, 0U)) == backflow_composition[0U] &&
      0.5 * (multispecies_backflow.independent_species[1U].view.unchecked(
                 outlet_owner, 0U) +
             multispecies_backflow.thermophysical_species_aliases[1U]
                 .unchecked(outlet_face, 0U)) == backflow_composition[1U];
  passed &= expect(
      multispecies_baseline_status && multispecies_exact_baseline.valid() &&
          multispecies_positive_status &&
          multispecies_exact_positive.valid() &&
          multispecies_positive_replay_status &&
          multispecies_exact_positive_replay.valid() &&
          multispecies_no_revision_mutation_rejected &&
          multispecies_positive.pressure_stage.canonical_lineage() ==
              multispecies_positive_replay.pressure_stage
                  .canonical_lineage() &&
          multispecies_positive.pressure_stage.scratch_binding() !=
              multispecies_positive_replay.pressure_stage.scratch_binding() &&
          multispecies_positive.velocity_stage.canonical_lineage() ==
              multispecies_positive_replay.velocity_stage
                  .canonical_lineage() &&
          multispecies_positive.velocity_stage.scratch_binding() !=
              multispecies_positive_replay.velocity_stage.scratch_binding() &&
          multispecies_positive.flux_stage.canonical_lineage() ==
              multispecies_positive_replay.flux_stage.canonical_lineage() &&
          multispecies_positive.flux_stage.scratch_binding() !=
              multispecies_positive_replay.flux_stage.scratch_binding() &&
          multispecies_positive.final_boundary.canonical_lineage() ==
              multispecies_positive_replay.final_boundary
                  .canonical_lineage() &&
          multispecies_positive.final_boundary.scratch_binding() !=
              multispecies_positive_replay.final_boundary.scratch_binding() &&
          multispecies_exact_positive.candidate_state_provenance() ==
              multispecies_exact_positive_replay
                  .candidate_state_provenance() &&
          multispecies_exact_positive.candidate_mass_flux_provenance() ==
              multispecies_exact_positive_replay
                  .candidate_mass_flux_provenance() &&
          multispecies_backflow_staged &&
          multispecies_backflow.final_boundary.valid() &&
          multispecies_positive.independent_species.size() == 2U &&
          multispecies_positive.thermophysical_species_aliases.size() == 2U &&
          multispecies_positive.independent_species[0U].view.field == 230U &&
          multispecies_positive.independent_species[1U].view.field == 231U &&
          multispecies_positive.thermophysical_species_aliases[0U].field ==
              8U &&
          multispecies_positive.thermophysical_species_aliases[1U].field ==
              9U &&
          multispecies_positive.independent_species[0U].view.field !=
              multispecies_positive.thermophysical_species_aliases[0U].field &&
          multispecies_positive.independent_species[1U].view.field !=
              multispecies_positive.thermophysical_species_aliases[1U].field &&
          multispecies_oracle_status && inlet_ghost_targets &&
          backflow_face_targets &&
          close_flux(multispecies_positive.final_flux.x.unchecked({0, 0, 0}),
                     expected_inlet_flux) &&
          close_flux(multispecies_backflow.final_flux.x.unchecked(outlet_face),
                     expected_backflow_flux) &&
          !close_flux(multispecies_backflow.final_flux.x.unchecked(outlet_face),
                      wrong_backflow_flux),
      "three-species/two-independent real chain binds semantic live scalar IDs to distinct candidate payloads and certifies inlet/backflow EOS flux");

  CandidateBoundaryFixture mass_fixture;
  CandidateBoundaryFixtureSpec mass_spec;
  mass_spec.inlet = CandidateBoundaryInlet::mass_flow;
  mass_spec.mass_flow_rate = 0.375;
  CandidateBoundaryScratch mass_candidate;
  const bool mass_initialized = mass_fixture.initialize(MPI_COMM_SELF, mass_spec);
  const bool mass_staged =
      mass_initialized &&
      mass_fixture.stage(0.0, 8.0, 0.0, 25000U, mass_candidate);
  if (!mass_staged)
    std::cerr << "candidate-boundary mass=" << mass_initialized << " step="
              << mass_fixture.diagnostic_step << " status="
              << static_cast<unsigned>(mass_fixture.diagnostic_status.code)
              << '/' << mass_fixture.diagnostic_status.detail << '\n';
  passed &= expect(
      mass_staged,
      "real compiler chain finalizes a mass-flow inlet baseline");
  if (mass_candidate.final_boundary.valid()) {
    double achieved = 0.0;
    for (std::int32_t z = 0; z < mass_fixture.patch.cells.z; ++z)
      for (std::int32_t y = 0; y < mass_fixture.patch.cells.y; ++y)
        achieved += mass_candidate.final_flux.x.unchecked({0, y, z});
    passed &= expect(
        std::abs(achieved - mass_spec.mass_flow_rate) <=
            64.0 * std::numeric_limits<double>::epsilon(),
        "mass-flow inlet globally normalizes capacity to the configured target");
  }

  CandidateBoundaryFixture backflow_fixture;
  CandidateBoundaryFixtureSpec backflow_spec;
  backflow_spec.allow_backflow = true;
  backflow_spec.backflow_velocity = -0.4;
  backflow_spec.backflow_temperature = 360.0;
  CandidateBoundaryScratch inward;
  CandidateBoundaryScratch recovered;
  const bool backflow_initialized =
      backflow_fixture.initialize(MPI_COMM_SELF, backflow_spec);
  const bool inward_staged =
      backflow_initialized &&
      backflow_fixture.stage(1.0, -1000.0, 0.0, 26000U, inward);
  if (!inward_staged)
    std::cerr << "candidate-boundary inward=" << backflow_initialized << " step="
              << backflow_fixture.diagnostic_step << " status="
              << static_cast<unsigned>(backflow_fixture.diagnostic_status.code)
              << '/' << backflow_fixture.diagnostic_status.detail << '\n';
  passed &= expect(
      inward_staged &&
          inward.mechanical_flux.x.unchecked(
              {backflow_fixture.patch.cells.x, 0, 0}) < 0.0 &&
          inward.final_flux.x.unchecked(
              {backflow_fixture.patch.cells.x, 0, 0}) < 0.0 &&
          inward.final_boundary.outlet_fixed_point_iterations() == 1U,
      "pressure outlet classifies first inward provisional flux and closes configured backflow");
  const bool recovery_staged = backflow_initialized &&
      backflow_fixture.stage(1.0, 1000.0, 0.0, 27000U, recovered);
  if (!recovery_staged)
    std::cerr << "candidate-boundary recovery step="
              << backflow_fixture.diagnostic_step << " status="
              << static_cast<unsigned>(backflow_fixture.diagnostic_status.code)
              << '/' << backflow_fixture.diagnostic_status.detail << '\n';
  passed &= expect(
      recovery_staged &&
          recovered.mechanical_flux.x.unchecked(
              {backflow_fixture.patch.cells.x, 0, 0}) > 0.0 &&
          same_double_bits(
              recovered.final_flux.x.unchecked(
                  {backflow_fixture.patch.cells.x, 0, 0}),
              recovered.mechanical_flux.x.unchecked(
                  {backflow_fixture.patch.cells.x, 0, 0})) &&
          recovered.final_boundary.outlet_fixed_point_iterations() == 0U,
      "pressure outlet recovery returns bitwise to provisional outflow");

  CandidateBoundaryFixture disabled_backflow_fixture;
  CandidateBoundaryFixtureSpec disabled_backflow_spec;
  disabled_backflow_spec.allow_backflow = false;
  CandidateBoundaryScratch rejected_inward;
  CandidateBoundaryScratch rejected_replay;
  const bool disabled_initialized = disabled_backflow_fixture.initialize(
      MPI_COMM_SELF, disabled_backflow_spec);
  const bool disabled_inward_staged =
      disabled_initialized && disabled_backflow_fixture.stage(
                                  1.0, -1000.0, 0.0, 27500U,
                                  rejected_inward);
  passed &= expect(
      disabled_initialized && !disabled_inward_staged &&
          disabled_backflow_fixture.diagnostic_status.code ==
              StatusCode::rejected_step &&
          !rejected_inward.final_boundary.valid() &&
          face_flux_all_equal(as_const(rejected_inward.final_flux), -991.0) &&
          disabled_backflow_fixture.stage(
              1.0, 1000.0, 0.0, 27600U, rejected_replay) &&
          rejected_replay.final_boundary.valid(),
      "real compiler chain collectively rejects disabled first backflow with zero output and preserves replay authority");

  CandidateBoundaryFixture ibm_fixture;
  CandidateBoundaryFixtureSpec ibm_spec;
  ibm_spec.immersed = true;
  ibm_spec.cells_per_axis = 16;
  CandidateBoundaryScratch ibm_baseline;
  const bool ibm_staged =
      ibm_fixture.initialize(MPI_COMM_SELF, ibm_spec) &&
      ibm_fixture.stage(0.0, 16.0, 0.0, 27700U, ibm_baseline);
  passed &= expect(
      ibm_staged && ibm_fixture.immersed_link_count() > 0U &&
          static_cast<bool>(ibm_fixture.validate_zero_interface_flux(
              as_const(ibm_baseline.final_flux))) &&
          ibm_baseline.final_boundary.ibm_donor_lineage() != 0U &&
          ibm_baseline.final_boundary.ibm_geometry_lineage() != 0U &&
          ibm_baseline.final_boundary.ibm_zero_interface_lineage() != 0U,
      "real compiler chain zeros every IBM control face and seals its geometry/donor lineage");
  if (ibm_staged) {
    PressureEnergyCandidateBoundaryFinalizer foreign;
    PressureEnergyCandidateBoundaryFinalizerBinding foreign_binding =
        ibm_fixture.finalizer_binding();
    foreign_binding.candidate_pressure_correction_donors =
        &ibm_fixture.live_pressure_donors;
    passed &= expect(
        PressureEnergyCandidateBoundaryFinalizer::bind(
            foreign_binding, foreign)
                .code == StatusCode::invalid_plan &&
            !foreign.ready(),
        "real IBM finalizer rejects a foreign live-donor authority even when geometry is shape-compatible");

    PressureEnergyCandidateBoundaryFinalizer foreign_reach;
    PressureEnergyCandidateBoundaryFinalizerBinding foreign_reach_binding =
        ibm_fixture.finalizer_binding();
    ++foreign_reach_binding.candidate_pressure_correction_donor_reach;
    passed &= expect(
        PressureEnergyCandidateBoundaryFinalizer::bind(
            foreign_reach_binding, foreign_reach)
                .code == StatusCode::invalid_plan &&
            !foreign_reach.ready(),
        "real IBM finalizer rejects a foreign nonzero donor reach");

    CandidateBoundaryFixture foreign_geometry_fixture;
    const bool foreign_initialized =
        foreign_geometry_fixture.initialize(MPI_COMM_SELF, ibm_spec);
    PressureEnergyCandidateBoundaryFinalizer foreign_geometry_finalizer;
    PressureEnergyCandidateBoundaryFinalizerBinding
        foreign_geometry_binding = ibm_fixture.finalizer_binding();
    foreign_geometry_binding.immersed_interface =
        &foreign_geometry_fixture.immersed_interface;
    passed &= expect(
        foreign_initialized &&
            PressureEnergyCandidateBoundaryFinalizer::bind(
                foreign_geometry_binding, foreign_geometry_finalizer)
                    .code == StatusCode::invalid_plan &&
            !foreign_geometry_finalizer.ready(),
        "real IBM finalizer rejects a foreign same-shape interface geometry authority");

    PisoFrozenMomentumExactCandidateCertificate ibm_exact_baseline;
    const Status ibm_baseline_status =
        ibm_fixture.coupler.certify_frozen_momentum_exact_baseline(
            ibm_fixture.authority, ibm_baseline.pressure_stage,
            ibm_baseline.velocity_stage, ibm_baseline.flux_stage,
            ibm_fixture.exact_input(ibm_baseline), ibm_fixture.reductions,
            ibm_exact_baseline);
    CandidateBoundaryScratch ibm_selected;
    PisoFrozenMomentumExactCandidateCertificate ibm_exact_selected;
    const bool ibm_selected_staged =
        ibm_baseline_status &&
        ibm_fixture.stage(0.5, 16.0, 0.0, 27800U, ibm_selected);
    const Status ibm_selected_status =
        ibm_selected_staged
            ? ibm_fixture.coupler.certify_frozen_momentum_exact_candidate(
                  ibm_fixture.authority, ibm_exact_baseline,
                  ibm_selected.pressure_stage, ibm_selected.velocity_stage,
                  ibm_selected.flux_stage,
                  ibm_fixture.exact_input(ibm_selected),
                  ibm_fixture.reductions, ibm_exact_selected)
            : Status{StatusCode::invalid_plan, 0U};
    passed &= expect(
        ibm_baseline_status && ibm_exact_baseline.valid() &&
            ibm_selected_status && ibm_exact_selected.valid() &&
            static_cast<bool>(ibm_fixture.validate_zero_interface_flux(
                as_const(ibm_selected.final_flux))),
        "real IBM alpha-zero and positive exact candidates retain the final-boundary issuer chain");

    if (ibm_exact_selected.valid()) {
      PressureEnergyGlobalizationSample ibm_baseline_sample;
      ibm_baseline_sample.alpha = 0.0;
      ibm_baseline_sample.global_normalized_continuity = 1.0;
      ibm_baseline_sample.global_normalized_energy = 1.0;
      ibm_baseline_sample.thermodynamically_admissible = true;
      ibm_baseline_sample.state_and_flux_finite = true;
      ibm_baseline_sample.corrector = ibm_exact_baseline.corrector();
      ibm_baseline_sample.target_time = ibm_exact_baseline.target_time();
      ibm_baseline_sample.correction_direction =
          ibm_exact_baseline.correction_direction();
      ibm_baseline_sample.state_provenance =
          ibm_exact_baseline.candidate_state_provenance();
      ibm_baseline_sample.mass_flux_provenance =
          ibm_exact_baseline.candidate_mass_flux_provenance();
      std::array<PressureEnergyGlobalizationSample,
                 kPressureEnergyGlobalizationCandidateCount>
          ibm_samples{};
      for (std::size_t index = 0U; index < ibm_samples.size(); ++index) {
        ibm_samples[index] = ibm_baseline_sample;
        ibm_samples[index].alpha =
            std::ldexp(1.0, -static_cast<int>(index));
        ibm_samples[index].global_normalized_continuity = 2.0;
        ibm_samples[index].global_normalized_energy = 2.0;
        ibm_samples[index].state_provenance = 32000U + index;
        ibm_samples[index].mass_flux_provenance = 33000U + index;
      }
      ibm_samples[1U].global_normalized_continuity = 0.2;
      ibm_samples[1U].global_normalized_energy = 0.2;
      ibm_samples[1U].state_provenance =
          ibm_exact_selected.candidate_state_provenance();
      ibm_samples[1U].mass_flux_provenance =
          ibm_exact_selected.candidate_mass_flux_provenance();
      PressureEnergyGlobalizationSelectionCertificate ibm_selection;
      FaceFluxStorage ibm_commit_storage;
      FaceFluxView ibm_committed_flux;
      PisoStateCorrectionCertificate ibm_committed;
      passed &= expect(
          static_cast<bool>(select_pressure_energy_globalization(
              ibm_baseline_sample,
              {ibm_samples.data(), ibm_samples.size()}, ibm_selection)) &&
              ibm_selection.valid() &&
              static_cast<bool>(FaceFluxStorage::allocate_workspace(
                  ibm_fixture.patch.cells, 1U, ibm_commit_storage)) &&
              static_cast<bool>(ibm_commit_storage.workspace_view(
                  0U, 34000U, ibm_committed_flux)) &&
              static_cast<bool>(
                  ibm_fixture.coupler
                      .commit_frozen_momentum_coupled_trial_state(
                          ibm_fixture.authority, ibm_exact_selected,
                          ibm_selection,
                          {ibm_fixture.velocity.view,
                           ibm_fixture.pressure.view,
                           ibm_fixture.enthalpy.view,
                           ibm_fixture.density.view,
                           ibm_fixture.temperature.view},
                          ibm_committed_flux, ibm_fixture.reductions,
                          ibm_committed)) &&
              ibm_committed.valid() &&
              static_cast<bool>(ibm_fixture.validate_zero_interface_flux(
                  as_const(ibm_committed_flux))),
          "real IBM selector and atomic commit publish the exact zero-interface final flux");
    }
  }

  if (staged) {
    PisoFrozenMomentumExactCandidateCertificate exact_baseline;
    const Status exact_baseline_status =
        fixture.coupler.certify_frozen_momentum_exact_baseline(
            fixture.authority, baseline.pressure_stage,
            baseline.velocity_stage, baseline.flux_stage,
            fixture.exact_input(baseline), fixture.reductions,
            exact_baseline);
    if (!exact_baseline_status)
      std::cerr << "candidate-boundary exact-base="
                << static_cast<unsigned>(exact_baseline_status.code) << '/'
                << exact_baseline_status.detail << '\n';
    passed &= expect(
        static_cast<bool>(exact_baseline_status) &&
            exact_baseline.valid() && exact_baseline.alpha() == 0.0,
        "open exact baseline consumes the independently issued final-boundary certificate");

    PisoFrozenMomentumExactCandidateInput open_with_compressibility =
        fixture.exact_input(baseline);
    open_with_compressibility.thermodynamic.pressure_compressibility =
        as_const(baseline.pressure_compressibility.view);
    PisoFrozenMomentumExactCandidateCertificate cleared_open_baseline =
        exact_baseline;
    passed &= expect(
        fixture.coupler
                    .certify_frozen_momentum_exact_baseline(
                        fixture.authority, baseline.pressure_stage,
                        baseline.velocity_stage, baseline.flux_stage,
                        open_with_compressibility, fixture.reductions,
                        cleared_open_baseline)
                    .code == StatusCode::invalid_plan &&
            !cleared_open_baseline.valid(),
        "boundary-absolute exact baseline rejects foreign closed-gauge pressure-compressibility authority and clears a reused certificate");

    CandidateBoundaryScratch selected;
    PisoFrozenMomentumExactCandidateCertificate exact_selected;
    const bool selected_staged =
        fixture.stage(0.5, 16.0, 4.0, 28000U, selected);
    Status exact_selected_status{StatusCode::invalid_plan, 0U};
    if (selected_staged)
      exact_selected_status =
          fixture.coupler.certify_frozen_momentum_exact_candidate(
              fixture.authority, exact_baseline, selected.pressure_stage,
              selected.velocity_stage, selected.flux_stage,
              fixture.exact_input(selected), fixture.reductions,
              exact_selected);
    if (!selected_staged || !exact_selected_status)
      std::cerr << "candidate-boundary selected-stage=" << selected_staged
                << " step=" << fixture.diagnostic_step << " stage-status="
                << static_cast<unsigned>(fixture.diagnostic_status.code) << '/'
                << fixture.diagnostic_status.detail << " exact="
                << static_cast<unsigned>(exact_selected_status.code) << '/'
                << exact_selected_status.detail << '\n';
    passed &= expect(
        selected_staged && static_cast<bool>(exact_selected_status) &&
            exact_selected.valid() && exact_selected.alpha() == 0.5 &&
            exact_selected.correction_direction() ==
                exact_baseline.correction_direction(),
        "open positive candidate replays the same raw dp/dh direction and final physical flux");

    PressureEnergyGlobalizationSample baseline_sample;
    baseline_sample.alpha = 0.0;
    baseline_sample.global_normalized_continuity = 1.0;
    baseline_sample.global_normalized_energy = 1.0;
    baseline_sample.thermodynamically_admissible = true;
    baseline_sample.state_and_flux_finite = true;
    baseline_sample.corrector = exact_baseline.corrector();
    baseline_sample.target_time = exact_baseline.target_time();
    baseline_sample.correction_direction =
        exact_baseline.correction_direction();
    baseline_sample.state_provenance =
        exact_baseline.candidate_state_provenance();
    baseline_sample.mass_flux_provenance =
        exact_baseline.candidate_mass_flux_provenance();
    std::array<PressureEnergyGlobalizationSample,
               kPressureEnergyGlobalizationCandidateCount>
        samples{};
    for (std::size_t index = 0U; index < samples.size(); ++index) {
      samples[index] = baseline_sample;
      samples[index].alpha = std::ldexp(1.0, -static_cast<int>(index));
      samples[index].global_normalized_continuity = 2.0;
      samples[index].global_normalized_energy = 2.0;
      samples[index].state_provenance = 29000U + index;
      samples[index].mass_flux_provenance = 30000U + index;
    }
    samples[1U].global_normalized_continuity = 0.2;
    samples[1U].global_normalized_energy = 0.2;
    samples[1U].state_provenance =
        exact_selected.candidate_state_provenance();
    samples[1U].mass_flux_provenance =
        exact_selected.candidate_mass_flux_provenance();
    PressureEnergyGlobalizationSelectionCertificate selection;
    passed &= expect(
        static_cast<bool>(select_pressure_energy_globalization(
            baseline_sample, {samples.data(), samples.size()}, selection)) &&
            selection.valid() && selection.alpha == 0.5,
        "open selector binds exact candidate state and final boundary flux provenance");

    FaceFluxStorage commit_storage;
    FaceFluxView committed_flux;
    PisoStateCorrectionCertificate committed;
    passed &= expect(
        static_cast<bool>(FaceFluxStorage::allocate_workspace(
            cells, 1U, commit_storage)) &&
            static_cast<bool>(commit_storage.workspace_view(
                0U, 31000U, committed_flux)) &&
            static_cast<bool>(
                fixture.coupler.commit_frozen_momentum_coupled_trial_state(
                    fixture.authority, exact_selected, selection,
                    {fixture.velocity.view, fixture.pressure.view,
                     fixture.enthalpy.view, fixture.density.view,
                     fixture.temperature.view},
                    committed_flux, fixture.reductions, committed)) &&
            committed.valid() &&
            face_flux_values(as_const(committed_flux)) ==
                face_flux_values(as_const(selected.final_flux)),
        "open selected candidate atomically publishes its final physical flux");
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_exactly_two_compile_authority();
  passed &= test_coupler_lifecycle_and_mutations();
  passed &= test_frozen_momentum_candidate_alpha_zero_exact();
  passed &= test_frozen_momentum_candidate_c1_c2_mapping();
  passed &= test_frozen_momentum_exact_publication_authority();
  passed &= test_pressure_energy_refinement_authority();
  passed &= test_frozen_momentum_stationary_pending_lineage();
  passed &= test_c1_full_refresh_four_block_sensitivity_red();
  passed &= test_open_boundary_terminal_authority();
  passed &= test_candidate_boundary_compiler_fixture();
  MPI_Finalize();
  return passed ? 0 : 1;
}
