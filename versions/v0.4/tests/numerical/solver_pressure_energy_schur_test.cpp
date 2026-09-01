// SPDX-License-Identifier: Apache-2.0

#ifndef HUNDUN_V04_ENABLE_TEST_ACCESS
#define HUNDUN_V04_ENABLE_TEST_ACCESS 1
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <type_traits>
#include <vector>

#include "hundun/v04_flow.hpp"

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void* allocate(std::size_t bytes) {
  if (enabled.load(std::memory_order_relaxed))
    count.fetch_add(1U, std::memory_order_relaxed);
  void* const result = std::malloc(bytes == 0U ? 1U : bytes);
  if (result == nullptr)
    throw std::bad_alloc{};
  return result;
}

void* allocate_aligned(std::size_t bytes, std::size_t alignment) {
  if (enabled.load(std::memory_order_relaxed))
    count.fetch_add(1U, std::memory_order_relaxed);
  void* result = nullptr;
  if (posix_memalign(&result, alignment, bytes == 0U ? alignment : bytes) !=
          0 ||
      result == nullptr)
    throw std::bad_alloc{};
  return result;
}

class Guard {
public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }
};

} // namespace allocation_observer

void* operator new(std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new[](std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      bytes, static_cast<std::size_t>(alignment));
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      bytes, static_cast<std::size_t>(alignment));
}
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(bytes);
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t bytes, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(bytes);
  } catch (...) {
    return nullptr;
  }
}
void* operator new(std::size_t bytes, std::align_val_t alignment,
                   const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate_aligned(
        bytes, static_cast<std::size_t>(alignment));
  } catch (...) {
    return nullptr;
  }
}
void* operator new[](std::size_t bytes, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate_aligned(
        bytes, static_cast<std::size_t>(alignment));
  } catch (...) {
    return nullptr;
  }
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
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t,
                     const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t,
                       const std::nothrow_t&) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

static_assert(
    std::is_base_of_v<LinearOperator, PressureEnergyPressureFluxOperator>,
    "the conservative pressure-to-energy face response is a Schur block");
static_assert(std::is_base_of_v<LinearOperator, PressureEnergyEnthalpyOperator>,
              "the frozen-spatial enthalpy response is a Schur block");

constexpr std::size_t kSize = 3U;
using Vector = std::array<double, kSize>;
using Matrix = std::array<double, kSize * kSize>;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

bool close(double actual, double expected, double tolerance = 2.0e-13) {
  return std::isfinite(actual) &&
         std::abs(actual - expected) <=
             tolerance * std::max(1.0, std::abs(expected));
}

Vector multiply(const Matrix& matrix, const Vector& value) {
  Vector result{};
  for (std::size_t row = 0U; row < kSize; ++row)
    for (std::size_t column = 0U; column < kSize; ++column)
      result[row] += matrix[row * kSize + column] * value[column];
  return result;
}

Matrix exact_schur(const Matrix& continuity_pressure,
                   const Vector& continuity_enthalpy,
                   const Matrix& energy_pressure,
                   const Matrix& energy_enthalpy) {
  Matrix result = energy_pressure;
  for (std::size_t row = 0U; row < kSize; ++row)
    for (std::size_t column = 0U; column < kSize; ++column)
      for (std::size_t eliminated = 0U; eliminated < kSize; ++eliminated)
        result[row * kSize + column] -=
            energy_enthalpy[row * kSize + eliminated] /
            continuity_enthalpy[eliminated] *
            continuity_pressure[eliminated * kSize + column];
  return result;
}

bool solve_dense(Matrix matrix, Vector rhs, Vector& solution) {
  for (std::size_t pivot = 0U; pivot < kSize; ++pivot) {
    std::size_t selected = pivot;
    for (std::size_t row = pivot + 1U; row < kSize; ++row)
      if (std::abs(matrix[row * kSize + pivot]) >
          std::abs(matrix[selected * kSize + pivot]))
        selected = row;
    if (!(std::abs(matrix[selected * kSize + pivot]) > 1.0e-14)) return false;
    if (selected != pivot) {
      for (std::size_t column = 0U; column < kSize; ++column)
        std::swap(matrix[pivot * kSize + column],
                  matrix[selected * kSize + column]);
      std::swap(rhs[pivot], rhs[selected]);
    }
    for (std::size_t row = pivot + 1U; row < kSize; ++row) {
      const double factor = matrix[row * kSize + pivot] /
                            matrix[pivot * kSize + pivot];
      for (std::size_t column = pivot; column < kSize; ++column)
        matrix[row * kSize + column] -=
            factor * matrix[pivot * kSize + column];
      rhs[row] -= factor * rhs[pivot];
    }
  }
  for (std::size_t reverse = 0U; reverse < kSize; ++reverse) {
    const std::size_t row = kSize - 1U - reverse;
    double value = rhs[row];
    for (std::size_t column = row + 1U; column < kSize; ++column)
      value -= matrix[row * kSize + column] * solution[column];
    solution[row] = value / matrix[row * kSize + row];
  }
  return true;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

struct OwnedFaces {
  std::vector<double> storage;
  FaceFieldView x{};
  FaceFieldView y{};
  FaceFieldView z{};
};
OwnedField field(FieldId id, RevisionToken revision,
                 StorageIdentity storage, Vector values = {}) {
  OwnedField result;
  result.storage.assign(values.begin(), values.end());
  result.view = {result.storage.data(),
                 {static_cast<std::int32_t>(kSize), 1, 1},
                 {},
                 1U,
                 kSize,
                 kSize,
                 kSize,
                 0U,
                 id,
                 revision,
                 storage,
                 7001U};
  return result;
}

OwnedField shaped_field(FieldId id, RevisionToken revision,
                        StorageIdentity storage, Int3 shape,
                        double value) {
  const std::size_t count = static_cast<std::size_t>(shape.x) *
                            static_cast<std::size_t>(shape.y) *
                            static_cast<std::size_t>(shape.z);
  OwnedField result;
  result.storage.assign(count, value);
  result.view = {result.storage.data(),
                 shape,
                 {},
                 1U,
                 static_cast<std::size_t>(shape.x),
                 static_cast<std::size_t>(shape.x) *
                     static_cast<std::size_t>(shape.y),
                 count,
                 0U,
                 id,
                 revision,
                 storage,
                 7001U};
  return result;
}

OwnedField ghosted_field(FieldId id, RevisionToken revision,
                         StorageIdentity storage, Int3 shape,
                         std::uint8_t ghosts, double value = 0.0) {
  const std::size_t nx = static_cast<std::size_t>(shape.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(shape.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(shape.z + 2 * ghosts);
  OwnedField result;
  result.storage.assign(nx * ny * nz, value);
  result.view = {
      result.storage.data() + ghosts + ghosts * nx + ghosts * nx * ny,
      shape,
      {ghosts, ghosts, ghosts},
      1U,
      nx,
      nx * ny,
      nx * ny * nz,
      0U,
      id,
      revision,
      storage,
      7001U};
  return result;
}

OwnedField ghosted_components_field(FieldId id, RevisionToken revision,
                                     StorageIdentity storage, Int3 shape,
                                     std::uint8_t ghosts,
                                     std::uint8_t components,
                                     double value = 0.0) {
  const std::size_t nx = static_cast<std::size_t>(shape.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(shape.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(shape.z + 2 * ghosts);
  OwnedField result;
  result.storage.assign(nx * ny * nz * components, value);
  result.view = {
      result.storage.data() + ghosts + ghosts * nx + ghosts * nx * ny,
      shape,
      {ghosts, ghosts, ghosts},
      components,
      nx,
      nx * ny,
      nx * ny * nz,
      0U,
      id,
      revision,
      storage,
      7001U};
  return result;
}

OwnedFaces face_bundle(Int3 cells, StorageIdentity storage,
                       RevisionDomainIdentity domain, double value = 0.0) {
  const Int3 x_shape{cells.x + 1, cells.y, cells.z};
  const Int3 y_shape{cells.x, cells.y + 1, cells.z};
  const Int3 z_shape{cells.x, cells.y, cells.z + 1};
  const auto count = [](Int3 shape) {
    return static_cast<std::size_t>(shape.x) *
           static_cast<std::size_t>(shape.y) *
           static_cast<std::size_t>(shape.z);
  };
  const std::size_t x_count = count(x_shape);
  const std::size_t y_count = count(y_shape);
  const std::size_t z_count = count(z_shape);
  OwnedFaces result;
  result.storage.assign(x_count + y_count + z_count, value);
  result.x = {result.storage.data(),
              x_shape,
              static_cast<std::size_t>(x_shape.x),
              static_cast<std::size_t>(x_shape.x) * x_shape.y,
              CartesianAxis::x,
              storage,
              domain};
  result.y = {result.storage.data() + x_count,
              y_shape,
              static_cast<std::size_t>(y_shape.x),
              static_cast<std::size_t>(y_shape.x) * y_shape.y,
              CartesianAxis::y,
              storage,
              domain};
  result.z = {result.storage.data() + x_count + y_count,
              z_shape,
              static_cast<std::size_t>(z_shape.x),
              static_cast<std::size_t>(z_shape.x) * z_shape.y,
              CartesianAxis::z,
              storage,
              domain};
  return result;
}

ConstFaceFluxView flux_view(const OwnedFaces& faces, RevisionToken revision) {
  return {
      as_const(faces.x), as_const(faces.y), as_const(faces.z), revision, {}};
}

enum class PressureFluxBoundary { neumann_dirichlet, periodic };

struct PressureFluxFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  HaloEngine halo;
};

struct EnthalpySpatialFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  CartesianKernelPlan kernels;
  HaloEngine halo;
};

bool make_enthalpy_spatial_fixture(EnthalpySpatialFixture& out,
                                   MPI_Comm communicator = MPI_COMM_SELF,
                                   bool periodic = false,
                                   bool limited_enthalpy = false) {
  int communicator_size = 0;
  if (MPI_Comm_size(communicator, &communicator_size) != MPI_SUCCESS ||
      communicator_size <= 0) {
    return false;
  }
  const std::int32_t edge = 4 * communicator_size;
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {static_cast<double>(edge), static_cast<double>(edge),
                static_cast<double>(edge)};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {edge, edge, edge};
  mesh.minimum_spacing = {1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = static_cast<std::uint64_t>(edge) * edge * edge;
  mesh.limits.max_memory_bytes_per_rank = 1U << 26U;
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = periodic ? 8700U : 8701U;
  model.turbulence = TurbulenceKind::none;
  model.pressure_reference = periodic
                                 ? PressureReferenceKind::closed_mass
                                 : PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = periodic ? BoundaryKind::periodic : BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
  }
  if (!periodic) {
    model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
    model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
    model.boundaries[0U].temperature = 300.0;
    model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
    model.boundaries[1U].pressure = 101325.0;
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = limited_enthalpy ? ConvectionScheme::limited_central2
                                            : ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::tvd2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  FieldRegistry registry;
  FieldId id = 0U;
  const bool fields = registry.require_field("rho", 1U, 1U, id) &&
                      registry.require_field("U", 3U, 1U, id) &&
                      registry.require_field("pi", 1U, 1U, id) &&
                      registry.require_field("h", 1U, 2U, id) && id == 3U &&
                      registry.require_field("T", 1U, 1U, id) &&
                      registry.require_field("dh", 1U, 2U, id) && id == 5U &&
                      registry.require_field("deltaT", 1U, 1U, id) && id == 6U;
  if (!fields ||
      !CartesianGeometryCompiler::compile(communicator, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(communicator, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes,
                                 out.time) ||
      !CartesianKernelPlan::compile(out.schemes, out.geometry, out.patch,
                                    out.boundary, out.kernels)) {
    return false;
  }
  const std::array<HaloFieldSpec, 2U> fields_to_exchange{{
      {5U, 2U, 1U},
      {6U, 1U, 1U},
  }};
  return static_cast<bool>(
      out.halo.reserve(communicator, out.patch,
                       {fields_to_exchange.data(), fields_to_exchange.size()},
                       out.boundary.halo_topology()));
}

bool prepare_frozen_enthalpy(const EnthalpySpatialFixture& fixture,
                             ConstFaceFluxView target_flux,
                             ConstFieldView target,
                             FrozenConvectionFaceOutput output,
                             FrozenConvectionContext context,
                             FrozenConvectionFaceField& frozen) {
  return static_cast<bool>(freeze_cartesian_target_convection_faces(
      fixture.kernels, fixture.schemes.enthalpy(), target_flux, target, 0U,
      context, output, frozen));
}

bool make_pressure_flux_fixture(PressureFluxBoundary kind,
                                PressureFluxFixture& out) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {3.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = {3, 1, 1};
  mesh.minimum_spacing = {1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 3U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 20U;
  ValidatedModel model;
  model.mesh = mesh;
  model.fingerprint = kind == PressureFluxBoundary::periodic ? 8101U : 8102U;
  model.turbulence = TurbulenceKind::none;
  model.pressure_reference = kind == PressureFluxBoundary::periodic
                                 ? PressureReferenceKind::closed_mass
                                 : PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : model.boundaries) {
    face.flow_kind = kind == PressureFluxBoundary::periodic
                         ? BoundaryKind::periodic
                         : BoundaryKind::symmetry;
    face.thermal_kind = BoundaryKind::none;
  }
  if (kind == PressureFluxBoundary::neumann_dirichlet) {
    model.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
    model.boundaries[0U].velocity = {1.0, 0.0, 0.0};
    model.boundaries[0U].temperature = 300.0;
    model.boundaries[1U].flow_kind = BoundaryKind::pressure_outlet;
    model.boundaries[1U].pressure = 101325.0;
  }
  model.schemes.momentum = ConvectionScheme::central2;
  model.schemes.enthalpy = ConvectionScheme::central2;
  model.schemes.species = ConvectionScheme::central2;
  model.schemes.passive_scalar = ConvectionScheme::central2;
  FieldRegistry registry;
  FieldId id = 0U;
  const bool compiled_fields =
      registry.require_field("rho", 1U, 1U, id) && id == 0U &&
      registry.require_field("U", 3U, 1U, id) && id == 1U &&
      registry.require_field("pi", 1U, 1U, id) && id == 2U &&
      registry.require_field("h", 1U, 1U, id) && id == 3U &&
      registry.require_field("T", 1U, 1U, id) && id == 4U;
  if (!compiled_fields ||
      !CartesianGeometryCompiler::compile(MPI_COMM_SELF, mesh, {}, out.geometry,
                                          out.patch) ||
      !BoundaryCompiler::compile(MPI_COMM_SELF, model, out.geometry, out.patch,
                                 registry, out.boundary, out.schemes,
                                 out.time)) {
    return false;
  }
  const std::array<HaloFieldSpec, 1U> fields{{{2U, 1U, 1U}}};
  return static_cast<bool>(out.halo.reserve(MPI_COMM_SELF, out.patch,
                                            {fields.data(), fields.size()},
                                            out.boundary.halo_topology()));
}

PressureEnergyPressureFluxBinding pressure_flux_binding(
    PressureFluxFixture& fixture, ConstFieldView temporal,
    const OwnedFaces& coefficient, const OwnedFaces& target,
    const OwnedFaces& enthalpy, PressureContinuityActivityView activity = {},
    LinearIdentity identity = {8201U, 8202U, 8203U, 8204U, 8205U}) {
  PisoIntermediateCertificate intermediate;
  intermediate.plan = 8211U;
  intermediate.r_au = 8212U;
  intermediate.h_by_a = 8213U;
  intermediate.pressure_face_coefficient = 8214U;
  intermediate.phi_h_by_a = 8215U;
  intermediate.trial_face_flux = 8216U;
  intermediate.dependency = 8217U;
  intermediate.corrector = 2U;
  intermediate.thermophysical_boundary_semantics = 8221U;
  intermediate.thermophysical_boundary_target = 8222U;
  intermediate.thermophysical_boundary_rank_local_binding = 8223U;
  intermediate.thermophysical_boundary_collective_lineage = 8224U;
  intermediate.thermophysical_boundary_rank_local_lineage = 8225U;
  PressureCorrectionCertificate pressure;
  pressure.plan = intermediate.plan;
  pressure.intermediate = intermediate.dependency;
  pressure.time = 8218U;
  pressure.geometry = fixture.geometry.topology_revision();
  pressure.numeric_boundary = fixture.boundary.revision();
  pressure.state = 8219U;
  pressure.corrector = 2U;
  pressure.thermophysical_boundary_semantics =
      intermediate.thermophysical_boundary_semantics;
  pressure.thermophysical_boundary_target =
      intermediate.thermophysical_boundary_target;
  pressure.thermophysical_boundary_rank_local_binding =
      intermediate.thermophysical_boundary_rank_local_binding;
  pressure.thermophysical_boundary_collective_lineage =
      intermediate.thermophysical_boundary_collective_lineage;
  pressure.thermophysical_boundary_rank_local_lineage =
      intermediate.thermophysical_boundary_rank_local_lineage;
  PressureEnergyPressureFluxBinding binding;
  binding.geometry = &fixture.geometry;
  binding.boundary = &fixture.boundary;
  binding.patch = fixture.patch;
  binding.services = {MPI_COMM_SELF, &fixture.halo, 8220U, 2U, 1U};
  binding.intermediate = intermediate;
  binding.pressure = pressure;
  binding.temporal_diagonal = temporal;
  binding.x_pressure_coefficient = as_const(coefficient.x);
  binding.y_pressure_coefficient = as_const(coefficient.y);
  binding.z_pressure_coefficient = as_const(coefficient.z);
  binding.target_flux = flux_view(target, 8221U);
  binding.frozen_face_enthalpy = {as_const(enthalpy.x),
                                  as_const(enthalpy.y),
                                  as_const(enthalpy.z),
                                  8222U,
                                  8223U,
                                  0U};
  binding.frozen_face_enthalpy.local_binding =
      pressure_energy_frozen_face_enthalpy_local_binding(
          binding.frozen_face_enthalpy);
  binding.activity = activity;
  binding.identity = identity;
  return binding;
}

Vector values(ConstFieldView field_view) {
  Vector result{};
  for (std::size_t index = 0U; index < kSize; ++index)
    result[index] = field_view.unchecked(
        {static_cast<std::int32_t>(index), 0, 0}, 0U);
  return result;
}

class DenseOperator final : public LinearOperator {
 public:
  DenseOperator(Matrix matrix, LinearIdentity identity,
                PlanFingerprint collective_fingerprint,
                LinearOperatorClass operator_class)
      : matrix_(matrix) {
    certificate_.identity = identity;
    certificate_.collective_fingerprint = collective_fingerprint;
    certificate_.local_shape = {static_cast<std::int32_t>(kSize), 1, 1};
    certificate_.operator_class = operator_class;
  }

  LinearOperatorCertificate certificate() const noexcept override {
    return certificate_;
  }

  Status apply(FieldView input, FieldView output) const noexcept override {
    const Vector result = multiply(matrix_, values(as_const(input)));
    for (std::size_t index = 0U; index < kSize; ++index)
      output.unchecked({static_cast<std::int32_t>(index), 0, 0}, 0U) =
          result[index];
    return {};
  }

  void replace_certificate(LinearOperatorCertificate certificate) noexcept {
    certificate_ = certificate;
  }

 private:
  Matrix matrix_{};
  LinearOperatorCertificate certificate_{};
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

bool authorize_generic_schur(PressureEnergySchurBinding& binding) {
  return binding.energy_pressure != nullptr &&
         binding.energy_enthalpy != nullptr &&
         static_cast<bool>(PressureEnergySchurBlockAuthority::
                               generic_algebraic_quasi_newton_for_test(
                                   binding.energy_pressure->certificate(),
                                   binding.energy_enthalpy->certificate(),
                                   binding.activity,
                                   binding.block_authority));
}

bool test_exact_schur_and_recovery() {
  const Matrix continuity_pressure{{
      2.0, -1.0, 0.0,
      0.5, 3.0, -0.2,
      0.0, 1.0, 1.5,
  }};
  const Vector continuity_enthalpy{{-2.0, -3.0, -4.0}};
  const Matrix energy_pressure{{
      1.0, 0.4, 0.0,
      -0.3, 2.0, 0.6,
      0.2, 0.0, 1.2,
  }};
  const Matrix energy_enthalpy{{
      0.5, -0.1, 0.2,
      0.3, 0.8, 0.0,
      -0.4, 0.2, 1.1,
  }};
  const LinearIdentity identity{401U, 402U, 403U, 404U, 405U};
  DenseOperator cp(continuity_pressure, identity, 101U,
                   LinearOperatorClass::spd);
  DenseOperator ep(energy_pressure, identity, 201U,
                   LinearOperatorClass::nonsymmetric);
  DenseOperator eh(energy_enthalpy, identity, 301U,
                   LinearOperatorClass::nonsymmetric);
  OwnedField ch = field(1U, 11U, 1011U, continuity_enthalpy);
  OwnedField ch_scale = field(2U, 12U, 1012U, {{2.0, 3.0, 4.0}});
  OwnedField continuity_workspace = field(3U, 13U, 1013U);
  OwnedField eliminated_workspace = field(4U, 14U, 1014U);
  OwnedField energy_workspace = field(5U, 15U, 1015U);

  PressureEnergySchurBinding binding;
  binding.continuity_pressure = &cp;
  binding.energy_pressure = &ep;
  binding.energy_enthalpy = &eh;
  binding.continuity_enthalpy_diagonal = as_const(ch.view);
  binding.continuity_enthalpy_row_scale = as_const(ch_scale.view);
  binding.workspace = {continuity_workspace.view, eliminated_workspace.view,
                       energy_workspace.view};
  binding.scaled_pivot_floor = 1.0e-12;
  bool passed = expect(authorize_generic_schur(binding),
                       "dense Schur receives an explicit generic-quasi authority");

  PressureEnergySchurOperator schur;
  PressureEnergyJacobianCertificate certificate;
  passed &= expect(static_cast<bool>(PressureEnergySchurOperator::bind(
                       binding, schur, certificate)),
                   "exact pressure-energy Schur binds");
  passed &= expect(certificate.valid() &&
                       certificate.schur.operator_class ==
                           LinearOperatorClass::nonsymmetric &&
                       certificate.sign_class ==
                           PressureEnergySchurSignClass::general &&
                       certificate.exact_block_equivalent &&
                       certificate.cell_local_continuity_enthalpy &&
                       certificate.native_mg_preconditioner_only &&
                       certificate.jacobian_scope ==
                           PressureEnergyJacobianScope::
                               generic_algebraic_quasi_newton &&
                       certificate.exact_algebraic_schur &&
                       !certificate.full_nonlinear_jacobian &&
                       close(certificate.minimum_scaled_abs_c_h, 1.0),
                   "Schur certificate is exact and never claims SPD/M-matrix");

  const Matrix expected_schur = exact_schur(
      continuity_pressure, continuity_enthalpy, energy_pressure,
      energy_enthalpy);
  passed &= expect(!close(expected_schur[1U], expected_schur[3U], 1.0e-14),
                   "oracle Schur is genuinely nonsymmetric");

  const Vector probe{{0.7, -1.2, 0.4}};
  OwnedField pressure = field(6U, 16U, 1016U, probe);
  OwnedField applied = field(7U, 17U, 1017U);
  passed &= expect(static_cast<bool>(schur.apply(pressure.view, applied.view)),
                   "exact Schur applies");
  const Vector expected_applied = multiply(expected_schur, probe);
  const Vector actual_applied = values(as_const(applied.view));
  for (std::size_t index = 0U; index < kSize; ++index)
    passed &= expect(close(actual_applied[index], expected_applied[index]),
                     "Schur action equals dense block elimination");

  const Vector continuity_residual{{0.3, -0.8, 0.5}};
  const Vector energy_residual{{-1.1, 0.2, 0.9}};
  OwnedField rc = field(8U, 18U, 1018U, continuity_residual);
  OwnedField re = field(9U, 19U, 1019U, energy_residual);
  OwnedField rhs = field(10U, 20U, 1020U);
  passed &= expect(static_cast<bool>(schur.form_pressure_rhs(
                       as_const(rc.view), as_const(re.view), rhs.view)),
                   "Schur right-hand side forms");
  Vector inv_ch_rc{};
  for (std::size_t index = 0U; index < kSize; ++index)
    inv_ch_rc[index] = continuity_residual[index] /
                       continuity_enthalpy[index];
  const Vector eh_inv_ch_rc = multiply(energy_enthalpy, inv_ch_rc);
  Vector expected_rhs{};
  for (std::size_t index = 0U; index < kSize; ++index)
    expected_rhs[index] = -energy_residual[index] + eh_inv_ch_rc[index];
  const Vector actual_rhs = values(as_const(rhs.view));
  for (std::size_t index = 0U; index < kSize; ++index)
    passed &= expect(close(actual_rhs[index], expected_rhs[index]),
                     "Schur right-hand side matches block oracle");

  Vector delta_pressure{};
  passed &= expect(solve_dense(expected_schur, expected_rhs, delta_pressure),
                   "dense Schur oracle solves");
  OwnedField dp = field(11U, 21U, 1021U, delta_pressure);
  OwnedField dh = field(12U, 22U, 1022U);
  passed &= expect(static_cast<bool>(schur.recover_enthalpy(
                       as_const(rc.view), dp.view, dh.view)),
                   "enthalpy correction recovers from the same block");
  const Vector delta_enthalpy = values(as_const(dh.view));
  const Vector cp_dp = multiply(continuity_pressure, delta_pressure);
  const Vector ep_dp = multiply(energy_pressure, delta_pressure);
  const Vector eh_dh = multiply(energy_enthalpy, delta_enthalpy);
  for (std::size_t index = 0U; index < kSize; ++index) {
    passed &= expect(close(cp_dp[index] +
                               continuity_enthalpy[index] *
                                   delta_enthalpy[index],
                           -continuity_residual[index]),
                     "recovered correction satisfies continuity block row");
    passed &= expect(close(ep_dp[index] + eh_dh[index],
                           -energy_residual[index]),
                     "recovered correction satisfies energy block row");
  }

  ch.view.unchecked({1, 0, 0}, 0U) = 0.0;
  PressureEnergySchurOperator rejected;
  PressureEnergyJacobianCertificate rejected_certificate = certificate;
  const Status pivot = PressureEnergySchurOperator::bind(
      binding, rejected, rejected_certificate);
  passed &= expect(pivot.code == StatusCode::rejected_step &&
                       !rejected_certificate.valid(),
                   "zero C_h pivot rejects without a false certificate");
  return passed;
}

bool test_enthalpy_spatial_binding_rejects_an_empty_contract() {
  PressureEnergyEnthalpyOperator operation;
  PressureEnergyEnthalpyCertificate certificate;
  const Status status =
      PressureEnergyEnthalpyOperator::bind({}, operation, certificate);
  bool passed =
      expect(status.code == StatusCode::invalid_plan && !certificate.valid(),
             "frozen-spatial E_h fails closed without target "
             "authorities");
  passed &= expect(operation.apply({}, {}).code == StatusCode::invalid_plan,
                   "unbound frozen-spatial E_h rejects apply without "
                   "dereferencing an absent halo service");
  return passed;
}

bool test_enthalpy_spatial_target_contract_binds() {
  EnthalpySpatialFixture fixture;
  bool passed = expect(make_enthalpy_spatial_fixture(fixture),
                       "frozen-spatial E_h fixture compiles");
  if (!passed)
    return false;
  const Int3 cells = fixture.patch.cells;
  OwnedField assembled = shaped_field(20U, 8702U, 8703U, cells, 20.0);
  OwnedField target = ghosted_field(3U, 8704U, 8705U, cells, 2U);
  OwnedField rho_h = shaped_field(21U, 8706U, 8707U, cells, -1.0e-6);
  OwnedField cp = ghosted_field(22U, 8708U, 8709U, cells, 1U, 1000.0);
  OwnedField lambda = ghosted_field(23U, 8710U, 8711U, cells, 1U, 2.0);
  OwnedField lambda_over_cp =
      ghosted_field(24U, 8712U, 8713U, cells, 1U, 0.002);
  OwnedField delta_temperature = ghosted_field(6U, 8714U, 8715U, cells, 1U);
  bool has_enthalpy_dirichlet = false;
  bool has_enthalpy_neumann = false;
  const Span<const BoundaryIndexSpan> boundary_spans = fixture.boundary.spans();
  for (std::size_t index = 0U; index < boundary_spans.size; ++index) {
    const BoundaryIndexSpan& span = boundary_spans.data[index];
    if (span.stage != BoundaryStage::enthalpy || span.field != 3U)
      continue;
    has_enthalpy_dirichlet =
        has_enthalpy_dirichlet || span.relation == BoundaryRelation::dirichlet;
    has_enthalpy_neumann = has_enthalpy_neumann ||
                           span.relation == BoundaryRelation::zero_gradient ||
                           span.relation == BoundaryRelation::normal_gradient;
  }
  passed &= expect(has_enthalpy_dirichlet && has_enthalpy_neumann,
                   "E_h fixture carries both Dirichlet and Neumann "
                   "homogeneous enthalpy closures");
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        target.view.unchecked({x, y, z}, 0U) =
            300000.0 + 31.0 * x + 7.0 * y + 3.0 * z;
        const double local_cp = 910.0 + 37.0 * x + 11.0 * y + 5.0 * z;
        const double local_lambda = 1.4 + 0.31 * x + 0.17 * y + 0.09 * z;
        cp.view.unchecked({x, y, z}, 0U) = local_cp;
        lambda.view.unchecked({x, y, z}, 0U) = local_lambda;
        lambda_over_cp.view.unchecked({x, y, z}, 0U) = local_lambda / local_cp;
      }
  passed &= expect(
      static_cast<bool>(apply_homogeneous_scalar_boundary_ghosts(
          BoundaryStage::enthalpy, fixture.boundary, 3U, target.view, 2U)),
      "target enthalpy owns a complete test boundary closure");
  std::array<FieldView, 3U> material_fields{cp.view, lambda.view,
                                            lambda_over_cp.view};
  passed &= expect(
      static_cast<bool>(apply_physical_zero_gradient(
          fixture.boundary, {material_fields.data(), material_fields.size()})),
      "nonuniform cp/lambda fields own finite coefficient ghosts");

  OwnedFaces target_flux = face_bundle(cells, 8716U, 8717U);
  for (std::int32_t z = 0; z < target_flux.x.extents.z; ++z)
    for (std::int32_t y = 0; y < target_flux.x.extents.y; ++y)
      for (std::int32_t x = 0; x < target_flux.x.extents.x; ++x)
        target_flux.x.unchecked({x, y, z}) = (x & 1) == 0 ? 0.25 : -0.15;
  for (std::int32_t z = 0; z < target_flux.y.extents.z; ++z)
    for (std::int32_t y = 0; y < target_flux.y.extents.y; ++y)
      for (std::int32_t x = 0; x < target_flux.y.extents.x; ++x)
        target_flux.y.unchecked({x, y, z}) = (y & 1) == 0 ? -0.08 : 0.12;
  for (std::int32_t z = 0; z < target_flux.z.extents.z; ++z)
    for (std::int32_t y = 0; y < target_flux.z.extents.y; ++y)
      for (std::int32_t x = 0; x < target_flux.z.extents.x; ++x)
        target_flux.z.unchecked({x, y, z}) = (z & 1) == 0 ? 0.04 : -0.06;
  const ConstFaceFluxView flux = flux_view(target_flux, 8718U);
  OwnedFaces frozen_storage = face_bundle(cells, 8719U, 8720U);
  OwnedFaces directional_storage = face_bundle(cells, 8721U, 8722U);
  const FrozenConvectionContext convection_context{8723U,
                                                   fixture.boundary.revision()};
  FrozenConvectionFaceField frozen;
  passed &= expect(prepare_frozen_enthalpy(
                       fixture, flux, as_const(target.view),
                       {frozen_storage.x, frozen_storage.y, frozen_storage.z},
                       convection_context, frozen),
                   "target enthalpy reconstruction freezes exactly");
  if (!passed)
    return false;

  PressureEnergyEnthalpyBinding binding;
  binding.geometry = &fixture.geometry;
  binding.kernels = &fixture.kernels;
  binding.boundary = &fixture.boundary;
  binding.patch = fixture.patch;
  binding.convection = fixture.schemes.enthalpy();
  binding.services = {MPI_COMM_SELF, &fixture.halo, 8724U, 5U, 6U};
  binding.authority = {{2.0, -2.0, 0.0, 1U},
                       8725U,
                       fixture.geometry.topology_revision(),
                       fixture.boundary.revision(),
                       8726U,
                       8727U,
                       convection_context.collective_semantics,
                       8728U,
                       8729U};
  binding.assembled_diagonal = as_const(assembled.view);
  binding.target_enthalpy = as_const(target.view);
  binding.density_enthalpy_derivative = as_const(rho_h.view);
  binding.heat_capacity = as_const(cp.view);
  binding.thermal_conductivity = as_const(lambda.view);
  binding.enthalpy_diffusivity = as_const(lambda_over_cp.view);
  binding.target_flux = flux;
  binding.convection_context = convection_context;
  binding.frozen_face_enthalpy = frozen;
  binding.workspace = {
      delta_temperature.view,
      {directional_storage.x, directional_storage.y, directional_storage.z}};
  binding.identity = {8730U, 8731U, 8732U, 8733U, 8734U};
  binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope;

  PressureEnergyEnthalpyOperator operation;
  PressureEnergyEnthalpyCertificate certificate;
  const Status status =
      PressureEnergyEnthalpyOperator::bind(binding, operation, certificate);
  passed &= expect(
      static_cast<bool>(status) && certificate.valid() &&
          certificate.exact_cartesian_spatial_response &&
          certificate.exact_temperature_space_conduction &&
          !certificate.ibm_spatial_derivative &&
          certificate.inactive_rows_identity &&
          certificate.inactive_interfaces_zero &&
          certificate.allocation_free_apply &&
          certificate.generalized_face_count == 0U &&
          certificate.geometry_fingerprint == fixture.geometry.fingerprint() &&
          certificate.halo_stage == 8724U &&
          certificate.enthalpy_variation_field == 5U &&
          certificate.temperature_variation_field == 6U &&
          certificate.delta_temperature_revision ==
              delta_temperature.view.revision &&
          certificate.delta_temperature_storage ==
              delta_temperature.view.storage_identity &&
          certificate.delta_temperature_revision_domain ==
              delta_temperature.view.revision_domain &&
          certificate.directional_enthalpy_storage ==
              directional_storage.x.storage_identity &&
          certificate.directional_enthalpy_revision_domain ==
              directional_storage.x.revision_domain &&
          certificate.halo_instance == fixture.halo.instance_identity() &&
          certificate.linearization_policy ==
              FrozenConvectionLinearizationPolicy::
                  semismooth_generalized_zero_slope,
      "frozen-spatial E_h binds target, material, branch and "
      "boundary authorities");
  if (!passed)
    return false;

  PressureEnergyEnthalpyBinding foreign_geometry = binding;
  ++foreign_geometry.authority.geometry;
  PressureEnergyEnthalpyOperator rejected_geometry;
  PressureEnergyEnthalpyCertificate rejected_geometry_certificate;
  passed &= expect(
      PressureEnergyEnthalpyOperator::bind(foreign_geometry, rejected_geometry,
                                           rejected_geometry_certificate)
                  .code == StatusCode::invalid_plan &&
          !rejected_geometry_certificate.valid(),
      "E_h refuses a geometry authority foreign to the bound metric plan");
  if (!passed)
    return false;

  HaloEngine reverse_id_halo;
  const std::array<HaloFieldSpec, 2U> reverse_id_specs{{
      {6U, 2U, 1U},
      {5U, 1U, 1U},
  }};
  passed &= expect(static_cast<bool>(reverse_id_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {reverse_id_specs.data(), reverse_id_specs.size()},
                       fixture.boundary.halo_topology())),
                   "reverse-ID E_h halo canonicalizes its field contract");
  OwnedField reverse_delta_temperature =
      ghosted_field(5U, 8950U, 8951U, cells, 1U);
  PressureEnergyEnthalpyBinding reverse_id_binding = binding;
  reverse_id_binding.services = {MPI_COMM_SELF, &reverse_id_halo, 8952U, 6U,
                                 5U};
  reverse_id_binding.workspace.delta_temperature =
      reverse_delta_temperature.view;
  reverse_id_binding.identity = {8953U, 8954U, 8955U, 8956U, 8957U};
  PressureEnergyEnthalpyOperator reverse_id_operation;
  PressureEnergyEnthalpyCertificate reverse_id_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyEnthalpyOperator::bind(
          reverse_id_binding, reverse_id_operation, reverse_id_certificate)) &&
          reverse_id_certificate.enthalpy_variation_field == 6U &&
          reverse_id_certificate.temperature_variation_field == 5U,
      "E_h binds when canonical Halo FieldId order reverses dh/deltaT "
      "semantic order");
  OwnedField reverse_direction =
      ghosted_field(6U, 8958U, 8959U, cells, 2U, 0.0);
  OwnedField reverse_output = shaped_field(27U, 8960U, 8961U, cells, -1.0);
  passed &= expect(static_cast<bool>(reverse_id_operation.apply(
                       reverse_direction.view, reverse_output.view)),
                   "reverse-ID E_h applies one dual-field exchange");
  if (!passed)
    return false;

  OwnedField variation = ghosted_field(5U, 8735U, 8736U, cells, 2U);
  OwnedField applied = shaped_field(25U, 8737U, 8738U, cells, -99.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        variation.view.unchecked({x, y, z}, 0U) =
            0.8 + 0.21 * x - 0.13 * y + 0.09 * z;
  const Status apply_status = operation.apply(variation.view, applied.view);
  passed &= expect(static_cast<bool>(apply_status),
                   "frozen-spatial E_h applies one complete target response");
  if (!passed)
    return false;
  OwnedFaces zero_flux_storage = face_bundle(cells, 8970U, 8971U, 0.0);
  const ConstFaceFluxView zero_flux = flux_view(zero_flux_storage, 8972U);
  OwnedFaces zero_flux_frozen_storage = face_bundle(cells, 8973U, 8974U);
  OwnedFaces zero_flux_directional_storage = face_bundle(cells, 8975U, 8976U);
  FrozenConvectionFaceField zero_flux_frozen;
  passed &= expect(prepare_frozen_enthalpy(
                       fixture, zero_flux, as_const(target.view),
                       {zero_flux_frozen_storage.x, zero_flux_frozen_storage.y,
                        zero_flux_frozen_storage.z},
                       convection_context, zero_flux_frozen),
                   "conduction-only oracle freezes a zero target flux");
  PressureEnergyEnthalpyBinding conduction_binding = binding;
  conduction_binding.target_flux = zero_flux;
  conduction_binding.frozen_face_enthalpy = zero_flux_frozen;
  conduction_binding.workspace.directional_enthalpy = {
      zero_flux_directional_storage.x, zero_flux_directional_storage.y,
      zero_flux_directional_storage.z};
  conduction_binding.identity = {8977U, 8978U, 8979U, 8980U, 8981U};
  PressureEnergyEnthalpyOperator conduction_operation;
  PressureEnergyEnthalpyCertificate conduction_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyEnthalpyOperator::bind(
          conduction_binding, conduction_operation, conduction_certificate)),
      "zero-flux E_h binds an explicit conduction-only spatial response");
  OwnedField conduction_applied = shaped_field(28U, 8982U, 8983U, cells, -98.0);
  passed &= expect(static_cast<bool>(conduction_operation.apply(
                       variation.view, conduction_applied.view)),
                   "zero-flux E_h applies the conduction-only response");
  if (!passed)
    return false;

  // Independent residual finite difference: reconstruct the nonlinear
  // target faces afresh at h+/-eps*dh and form conduction directly from the
  // uniform finite-volume metric.  It deliberately does not call the
  // operator's directional-face or diffusion helpers.
  constexpr double epsilon = 1.0e-3;
  OwnedField plus_h = ghosted_field(3U, 8739U, 8740U, cells, 2U);
  OwnedField minus_h = ghosted_field(3U, 8741U, 8742U, cells, 2U);
  OwnedField plus_t = ghosted_field(4U, 8743U, 8744U, cells, 1U);
  OwnedField minus_t = ghosted_field(4U, 8745U, 8746U, cells, 1U);
  for (std::int32_t z = -2; z < cells.z + 2; ++z) {
    for (std::int32_t y = -2; y < cells.y + 2; ++y) {
      for (std::int32_t x = -2; x < cells.x + 2; ++x) {
        const Int3 cell{x, y, z};
        const double base = target.view.unchecked(cell, 0U);
        const double direction = variation.view.unchecked(cell, 0U);
        plus_h.view.unchecked(cell, 0U) = base + epsilon * direction;
        minus_h.view.unchecked(cell, 0U) = base - epsilon * direction;
      }
    }
  }
  for (std::int32_t z = -1; z < cells.z + 1; ++z) {
    for (std::int32_t y = -1; y < cells.y + 1; ++y) {
      for (std::int32_t x = -1; x < cells.x + 1; ++x) {
        const Int3 cell{x, y, z};
        const double base =
            target.view.unchecked(cell, 0U) / cp.view.unchecked(cell, 0U);
        const double direction = delta_temperature.view.unchecked(cell, 0U);
        plus_t.view.unchecked(cell, 0U) = base + epsilon * direction;
        minus_t.view.unchecked(cell, 0U) = base - epsilon * direction;
      }
    }
  }
  OwnedFaces plus_faces = face_bundle(cells, 8747U, 8748U);
  OwnedFaces minus_faces = face_bundle(cells, 8749U, 8750U);
  FrozenConvectionFaceField plus_frozen;
  FrozenConvectionFaceField minus_frozen;
  passed &= expect(
      prepare_frozen_enthalpy(fixture, flux, as_const(plus_h.view),
                              {plus_faces.x, plus_faces.y, plus_faces.z},
                              convection_context, plus_frozen) &&
          prepare_frozen_enthalpy(fixture, flux, as_const(minus_h.view),
                                  {minus_faces.x, minus_faces.y, minus_faces.z},
                                  convection_context, minus_frozen),
      "finite-difference oracle reconstructs both nonlinear residuals");
  if (!passed)
    return false;

  const double dx = fixture.geometry.x().uniform_width();
  const double dy = fixture.geometry.y().uniform_width();
  const double dz = fixture.geometry.z().uniform_width();
  const double volume = dx * dy * dz;
  const auto shifted = [](Int3 cell, CartesianAxis axis, int amount) noexcept {
    if (axis == CartesianAxis::x)
      cell.x += amount;
    else if (axis == CartesianAxis::y)
      cell.y += amount;
    else
      cell.z += amount;
    return cell;
  };
  const auto material_transmissibility = [&](ConstFieldView material,
                                             CartesianAxis axis, Int3 face) {
    const Int3 left = shifted(face, axis, -1);
    const double left_value = material.unchecked(left, 0U);
    const double right_value = material.unchecked(face, 0U);
    const double spacing =
        axis == CartesianAxis::x ? dx : (axis == CartesianAxis::y ? dy : dz);
    const double area = axis == CartesianAxis::x
                            ? dy * dz
                            : (axis == CartesianAxis::y ? dx * dz : dx * dy);
    return area * 2.0 * left_value * right_value /
           (spacing * (left_value + right_value));
  };
  const auto selected_face = [](ConstFaceFieldView x, ConstFaceFieldView y,
                                ConstFaceFieldView z,
                                CartesianAxis axis) noexcept {
    return axis == CartesianAxis::x ? x : (axis == CartesianAxis::y ? y : z);
  };
  double largest_new_error = 0.0;
  double largest_conduction_operator_error = 0.0;
  double largest_convection_operator_error = 0.0;
  double largest_diagonal_red = 0.0;
  double largest_convection_only = 0.0;
  double largest_conduction_only = 0.0;
  OwnedField old_diagonal = shaped_field(26U, 8751U, 8752U, cells, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        double proxy_diagonal = 0.0;
        double plus_convection = 0.0;
        double minus_convection = 0.0;
        double plus_negative_diffusion = 0.0;
        double minus_negative_diffusion = 0.0;
        for (CartesianAxis axis :
             {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
          const Int3 plus_face = shifted(cell, axis, 1);
          const ConstFaceFieldView rate =
              selected_face(flux.x, flux.y, flux.z, axis);
          const ConstFaceFieldView plus_value =
              selected_face(plus_frozen.x, plus_frozen.y, plus_frozen.z, axis);
          const ConstFaceFieldView minus_value = selected_face(
              minus_frozen.x, minus_frozen.y, minus_frozen.z, axis);
          plus_convection +=
              rate.unchecked(plus_face) * plus_value.unchecked(plus_face) -
              rate.unchecked(cell) * plus_value.unchecked(cell);
          minus_convection +=
              rate.unchecked(plus_face) * minus_value.unchecked(plus_face) -
              rate.unchecked(cell) * minus_value.unchecked(cell);
          for (int direction : {-1, 1}) {
            const Int3 face = direction < 0 ? cell : shifted(cell, axis, 1);
            const Int3 neighbour = shifted(cell, axis, direction);
            proxy_diagonal += material_transmissibility(
                as_const(lambda_over_cp.view), axis, face);
            const double conductance =
                material_transmissibility(as_const(lambda.view), axis, face);
            plus_negative_diffusion +=
                conductance * (plus_t.view.unchecked(cell, 0U) -
                               plus_t.view.unchecked(neighbour, 0U));
            minus_negative_diffusion +=
                conductance * (minus_t.view.unchecked(cell, 0U) -
                               minus_t.view.unchecked(neighbour, 0U));
          }
        }
        const double temporal_density = binding.authority.bdf.a0 * volume *
                                        target.view.unchecked(cell, 0U) *
                                        rho_h.view.unchecked(cell, 0U);
        const double local = assembled.view.unchecked(cell, 0U) -
                             proxy_diagonal + temporal_density;
        const double residual_plus = local * plus_h.view.unchecked(cell, 0U) +
                                     plus_convection + plus_negative_diffusion;
        const double residual_minus = local * minus_h.view.unchecked(cell, 0U) +
                                      minus_convection +
                                      minus_negative_diffusion;
        const double expected =
            (residual_plus - residual_minus) / (2.0 * epsilon);
        const double convection_only =
            (plus_convection - minus_convection) / (2.0 * epsilon);
        const double conduction_only =
            (plus_negative_diffusion - minus_negative_diffusion) /
            (2.0 * epsilon);
        largest_convection_only =
            std::max(largest_convection_only, std::abs(convection_only));
        largest_conduction_only =
            std::max(largest_conduction_only, std::abs(conduction_only));
        const double actual = applied.view.unchecked(cell, 0U);
        const double actual_conduction =
            conduction_applied.view.unchecked(cell, 0U);
        const double expected_conduction =
            local * variation.view.unchecked(cell, 0U) + conduction_only;
        largest_new_error =
            std::max(largest_new_error, std::abs(actual - expected) /
                                            std::max(1.0, std::abs(expected)));
        largest_conduction_operator_error =
            std::max(largest_conduction_operator_error,
                     std::abs(actual_conduction - expected_conduction) /
                         std::max(1.0, std::abs(expected_conduction)));
        largest_convection_operator_error =
            std::max(largest_convection_operator_error,
                     std::abs((actual - actual_conduction) - convection_only) /
                         std::max(1.0, std::abs(convection_only)));
        old_diagonal.view.unchecked(cell, 0U) =
            assembled.view.unchecked(cell, 0U) + temporal_density;
        const double old_response = old_diagonal.view.unchecked(cell, 0U) *
                                    variation.view.unchecked(cell, 0U);
        largest_diagonal_red =
            std::max(largest_diagonal_red, std::abs(old_response - expected));
      }
    }
  }
  passed &= expect(largest_new_error < 2.0e-7,
                   "frozen-spatial E_h matches an independent nonlinear "
                   "residual finite difference");
  passed &=
      expect(largest_conduction_operator_error < 2.0e-7 &&
                 largest_convection_operator_error < 2.0e-7,
             "zero-flux conduction-only action and combined-minus-conduction "
             "convection-only action match independent residual derivatives");
  passed &= expect(largest_convection_only > 1.0e-3 &&
                       largest_conduction_only > 1.0e-6,
                   "positive/negative target flux and nonuniform lambda/cp "
                   "exercise convection-only, conduction-only and combined "
                   "responses");
  passed &= expect(largest_diagonal_red > 1.0e-3,
                   "the legacy diagonal E_h misses the certified spatial "
                   "response RED");

  OwnedField hot_output = shaped_field(27U, 8753U, 8754U, cells, -77.0);
  Status hot_status;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    hot_status = operation.apply(variation.view, hot_output.view);
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(hot_status) && hot_allocations == 0U,
                   "frozen-spatial E_h hot apply performs zero allocation");

  const Status alias_status = operation.apply(variation.view, variation.view);
  passed &= expect(alias_status.code == StatusCode::invalid_plan,
                   "frozen-spatial E_h rejects input/output aliasing");

  PressureEnergyEnthalpyCertificate& mutable_certificate =
      const_cast<PressureEnergyEnthalpyCertificate&>(
          operation.enthalpy_certificate());
  const StorageIdentity saved_delta_temperature_storage =
      mutable_certificate.delta_temperature_storage;
  ++mutable_certificate.delta_temperature_storage;
  const Status foreign_scratch =
      operation.apply(variation.view, hot_output.view);
  mutable_certificate.delta_temperature_storage =
      saved_delta_temperature_storage;
  passed &= expect(foreign_scratch.code == StatusCode::invalid_plan,
                   "foreign-but-shape-valid scratch identity invalidates "
                   "the bound E_h certificate");

  std::fill(hot_output.storage.begin(), hot_output.storage.end(), -1234.0);
  const std::vector<double> rejected_snapshot = hot_output.storage;
  const double saved_direction = variation.view.unchecked({0, 0, 0}, 0U);
  variation.view.unchecked({0, 0, 0}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  const Status nonfinite_direction =
      operation.apply(variation.view, hot_output.view);
  variation.view.unchecked({0, 0, 0}, 0U) = saved_direction;
  passed &= expect(nonfinite_direction.code == StatusCode::numerical_failure &&
                       hot_output.storage == rejected_snapshot,
                   "nonfinite dh fails collectively before output commit");

  const double saved_target = target.view.unchecked({1, 1, 1}, 0U);
  target.view.unchecked({1, 1, 1}, 0U) = saved_target + 0.25;
  const Status stale_target = operation.apply(variation.view, hot_output.view);
  target.view.unchecked({1, 1, 1}, 0U) = saved_target;
  passed &= expect(stale_target.code == StatusCode::invalid_plan &&
                       hot_output.storage == rejected_snapshot,
                   "raw target mutation invalidates the frozen branch before "
                   "output commit");

  const double saved_lambda = lambda.view.unchecked({1, 1, 1}, 0U);
  lambda.view.unchecked({1, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  const Status nonfinite_transport =
      operation.apply(variation.view, hot_output.view);
  lambda.view.unchecked({1, 1, 1}, 0U) = saved_lambda;
  passed &= expect(nonfinite_transport.code == StatusCode::numerical_failure &&
                       hot_output.storage == rejected_snapshot,
                   "nonfinite lambda fails atomically in exact temperature "
                   "diffusion");
  return passed;
}

bool test_enthalpy_spatial_periodic_mpi_and_inactive_interfaces() {
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS || size <= 0)
    return false;
  EnthalpySpatialFixture fixture;
  bool passed =
      expect(make_enthalpy_spatial_fixture(fixture, MPI_COMM_WORLD, true),
             "periodic MPI frozen-spatial E_h fixture compiles");
  if (!passed)
    return false;
  const Int3 cells = fixture.patch.cells;
  const Int3 global = fixture.geometry.global_cells();
  const auto count = [](Int3 shape) {
    return static_cast<std::size_t>(shape.x) *
           static_cast<std::size_t>(shape.y) *
           static_cast<std::size_t>(shape.z);
  };
  const auto offset = [](Int3 shape, Int3 index) {
    return static_cast<std::size_t>(index.x) +
           static_cast<std::size_t>(shape.x) *
               (static_cast<std::size_t>(index.y) +
                static_cast<std::size_t>(shape.y) *
                    static_cast<std::size_t>(index.z));
  };
  constexpr double pi = 3.141592653589793238462643383279502884;
  const auto phase = [pi](std::int32_t global_index, std::int32_t extent) {
    return 2.0 * pi * (static_cast<double>(global_index) + 0.5) /
           static_cast<double>(extent);
  };

  OwnedField assembled = shaped_field(30U, 8801U, 8802U, cells, 30.0);
  OwnedField target = ghosted_field(3U, 8803U, 8804U, cells, 2U);
  OwnedField rho_h = shaped_field(31U, 8805U, 8806U, cells, -1.2e-6);
  OwnedField cp = ghosted_field(32U, 8807U, 8808U, cells, 1U);
  OwnedField lambda = ghosted_field(33U, 8809U, 8810U, cells, 1U);
  OwnedField lambda_over_cp = ghosted_field(34U, 8811U, 8812U, cells, 1U);
  OwnedField delta_temperature = ghosted_field(6U, 8813U, 8814U, cells, 1U);
  OwnedField variation = ghosted_field(5U, 8815U, 8816U, cells, 2U);
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::int32_t gx = fixture.patch.begin.x + x;
        const std::int32_t gy = fixture.patch.begin.y + y;
        const std::int32_t gz = fixture.patch.begin.z + z;
        const double sx = std::sin(phase(gx, global.x));
        const double sy = std::sin(phase(gy, global.y));
        const double sz = std::sin(phase(gz, global.z));
        const double local_cp = 980.0 + 19.0 * sx + 13.0 * sy + 7.0 * sz;
        const double local_lambda = 1.8 + 0.21 * sx - 0.12 * sy + 0.08 * sz;
        target.view.unchecked(cell, 0U) =
            300000.0 + 41.0 * sx + 23.0 * sy - 17.0 * sz;
        cp.view.unchecked(cell, 0U) = local_cp;
        lambda.view.unchecked(cell, 0U) = local_lambda;
        lambda_over_cp.view.unchecked(cell, 0U) = local_lambda / local_cp;
        variation.view.unchecked(cell, 0U) =
            0.7 + 0.17 * sx - 0.11 * sy + 0.09 * sz;
      }
    }
  }

  HaloEngine target_halo;
  const std::array<HaloFieldSpec, 4U> target_specs{
      {{3U, 2U, 1U}, {32U, 1U, 1U}, {33U, 1U, 1U}, {34U, 1U, 1U}}};
  passed &= expect(static_cast<bool>(target_halo.reserve(
                       MPI_COMM_WORLD, fixture.patch,
                       {target_specs.data(), target_specs.size()},
                       fixture.boundary.halo_topology())),
                   "periodic target/material halo reserves");
  std::array<FieldView, 4U> target_fields{target.view, cp.view, lambda.view,
                                          lambda_over_cp.view};
  HaloTicket target_ticket;
  Status target_exchange = target_halo.begin(
      8817U, {target_fields.data(), target_fields.size()}, target_ticket);
  if (target_exchange)
    target_exchange = target_halo.finish(
        target_ticket, {target_fields.data(), target_fields.size()});
  passed &= expect(static_cast<bool>(target_exchange),
                   "periodic target/material ghosts exchange");
  if (!passed)
    return false;

  OwnedFaces target_flux = face_bundle(cells, 8818U, 8819U);
  const auto fill_axis_flux = [&](FaceFieldView face, CartesianAxis axis) {
    for (std::int32_t z = 0; z < face.extents.z; ++z) {
      for (std::int32_t y = 0; y < face.extents.y; ++y) {
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          const std::int32_t gx = fixture.patch.begin.x + x;
          const std::int32_t gy = fixture.patch.begin.y + y;
          const std::int32_t gz = fixture.patch.begin.z + z;
          const std::int32_t coordinate =
              axis == CartesianAxis::x ? gx
                                       : (axis == CartesianAxis::y ? gy : gz);
          const std::int32_t extent =
              axis == CartesianAxis::x
                  ? global.x
                  : (axis == CartesianAxis::y ? global.y : global.z);
          const double amplitude =
              axis == CartesianAxis::x
                  ? 0.23
                  : (axis == CartesianAxis::y ? -0.14 : 0.09);
          face.unchecked({x, y, z}) =
              amplitude * std::cos(2.0 * pi * coordinate / extent);
        }
      }
    }
  };
  fill_axis_flux(target_flux.x, CartesianAxis::x);
  fill_axis_flux(target_flux.y, CartesianAxis::y);
  fill_axis_flux(target_flux.z, CartesianAxis::z);
  const ConstFaceFluxView flux = flux_view(target_flux, 8820U);
  OwnedFaces frozen_storage = face_bundle(cells, 8821U, 8822U);
  OwnedFaces directional_storage = face_bundle(cells, 8823U, 8824U);
  const FrozenConvectionContext convection_context{8825U,
                                                   fixture.boundary.revision()};
  FrozenConvectionFaceField frozen;
  passed &= expect(prepare_frozen_enthalpy(
                       fixture, flux, as_const(target.view),
                       {frozen_storage.x, frozen_storage.y, frozen_storage.z},
                       convection_context, frozen),
                   "periodic target reconstruction freezes");
  if (!passed)
    return false;

  const Int3 x_shape{cells.x + 1, cells.y, cells.z};
  const Int3 y_shape{cells.x, cells.y + 1, cells.z};
  const Int3 z_shape{cells.x, cells.y, cells.z + 1};
  std::vector<std::uint8_t> active_cells(count(cells), 1U);
  std::vector<std::uint8_t> active_x(count(x_shape), 1U);
  std::vector<std::uint8_t> active_y(count(y_shape), 1U);
  std::vector<std::uint8_t> active_z(count(z_shape), 1U);
  const Int3 solid{cells.x / 2, cells.y / 2, cells.z / 2};
  active_cells[offset(cells, solid)] = 0U;
  active_x[offset(x_shape, solid)] = 0U;
  active_x[offset(x_shape, {solid.x + 1, solid.y, solid.z})] = 0U;
  active_y[offset(y_shape, solid)] = 0U;
  active_y[offset(y_shape, {solid.x, solid.y + 1, solid.z})] = 0U;
  active_z[offset(z_shape, solid)] = 0U;
  active_z[offset(z_shape, {solid.x, solid.y, solid.z + 1})] = 0U;
  PressureContinuityActivityView activity{
      {active_cells.data(), active_cells.size()},
      {active_x.data(), active_x.size()},
      {active_y.data(), active_y.size()},
      {active_z.data(), active_z.size()},
      static_cast<PlanFingerprint>(8830U + rank),
      8839U};

  PressureEnergyEnthalpyBinding binding;
  binding.geometry = &fixture.geometry;
  binding.kernels = &fixture.kernels;
  binding.boundary = &fixture.boundary;
  binding.patch = fixture.patch;
  binding.convection = fixture.schemes.enthalpy();
  binding.services = {MPI_COMM_WORLD, &fixture.halo, 8826U, 5U, 6U};
  binding.authority = {{2.0, -2.0, 0.0, 1U},
                       8827U,
                       fixture.geometry.topology_revision(),
                       fixture.boundary.revision(),
                       8828U,
                       8829U,
                       convection_context.collective_semantics,
                       8831U,
                       8832U};
  binding.assembled_diagonal = as_const(assembled.view);
  binding.target_enthalpy = as_const(target.view);
  binding.density_enthalpy_derivative = as_const(rho_h.view);
  binding.heat_capacity = as_const(cp.view);
  binding.thermal_conductivity = as_const(lambda.view);
  binding.enthalpy_diffusivity = as_const(lambda_over_cp.view);
  binding.target_flux = flux;
  binding.convection_context = convection_context;
  binding.frozen_face_enthalpy = frozen;
  binding.workspace = {
      delta_temperature.view,
      {directional_storage.x, directional_storage.y, directional_storage.z}};
  binding.activity = activity;
  binding.identity = {8840U, static_cast<RevisionToken>(8841U + rank), 8842U,
                      static_cast<RevisionToken>(8843U + rank),
                      static_cast<PlanFingerprint>(8844U + rank)};
  binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope;
  PressureEnergyEnthalpyOperator operation;
  PressureEnergyEnthalpyCertificate certificate;
  const Status bind_status =
      PressureEnergyEnthalpyOperator::bind(binding, operation, certificate);
  passed &= expect(static_cast<bool>(bind_status) && certificate.valid() &&
                       certificate.inactive_cells == 1U &&
                       !certificate.ibm_spatial_derivative,
                   "periodic MPI E_h binds inactive cell/face authority");
  unsigned long long local_collective = static_cast<unsigned long long>(
      certificate.linear.collective_fingerprint);
  unsigned long long minimum_collective = 0U;
  unsigned long long maximum_collective = 0U;
  MPI_Allreduce(&local_collective, &minimum_collective, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_collective, &maximum_collective, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
  passed &= expect(minimum_collective == maximum_collective &&
                       minimum_collective != 0U,
                   "E_h collective identity excludes rank-local partitions");
  if (!passed)
    return false;

  // Production IBM uses the same masked Cartesian spatial responses as this
  // fixture.  E_p deliberately omits pressure-work: the resulting typed
  // authority is a spatial quasi-Newton scope, never an exact Cartesian or
  // full nonlinear Jacobian claim.
  HaloEngine pressure_halo;
  const std::array<HaloFieldSpec, 1U> pressure_halo_specs{{{2U, 1U, 1U}}};
  passed &= expect(static_cast<bool>(pressure_halo.reserve(
                       MPI_COMM_WORLD, fixture.patch,
                       {pressure_halo_specs.data(), pressure_halo_specs.size()},
                       fixture.boundary.halo_topology())),
                   "IBM spatial quasi-Newton pressure halo reserves");
  OwnedField energy_pressure_temporal =
      shaped_field(37U, 8851U, 8852U, cells, 0.125);
  OwnedFaces pressure_coefficients = face_bundle(cells, 8853U, 8854U, 0.5);
  PisoIntermediateCertificate pressure_intermediate;
  pressure_intermediate.plan = 8855U;
  pressure_intermediate.r_au = 8856U;
  pressure_intermediate.h_by_a = 8857U;
  pressure_intermediate.pressure_face_coefficient = 8858U;
  pressure_intermediate.phi_h_by_a = 8859U;
  pressure_intermediate.trial_face_flux = 8860U;
  pressure_intermediate.dependency = 8861U;
  pressure_intermediate.corrector = 2U;
  pressure_intermediate.thermophysical_boundary_semantics = 8862U;
  pressure_intermediate.thermophysical_boundary_target = 8863U;
  pressure_intermediate.thermophysical_boundary_rank_local_binding = 8864U;
  pressure_intermediate.thermophysical_boundary_collective_lineage = 8865U;
  pressure_intermediate.thermophysical_boundary_rank_local_lineage = 8866U;
  PressureCorrectionCertificate pressure_certificate;
  pressure_certificate.plan = pressure_intermediate.plan;
  pressure_certificate.intermediate = pressure_intermediate.dependency;
  pressure_certificate.time = 8867U;
  pressure_certificate.geometry = fixture.geometry.topology_revision();
  pressure_certificate.numeric_boundary = fixture.boundary.revision();
  pressure_certificate.state = 8868U;
  pressure_certificate.corrector = 2U;
  pressure_certificate.thermophysical_boundary_semantics =
      pressure_intermediate.thermophysical_boundary_semantics;
  pressure_certificate.thermophysical_boundary_target =
      pressure_intermediate.thermophysical_boundary_target;
  pressure_certificate.thermophysical_boundary_rank_local_binding =
      pressure_intermediate.thermophysical_boundary_rank_local_binding;
  pressure_certificate.thermophysical_boundary_collective_lineage =
      pressure_intermediate.thermophysical_boundary_collective_lineage;
  pressure_certificate.thermophysical_boundary_rank_local_lineage =
      pressure_intermediate.thermophysical_boundary_rank_local_lineage;
  PressureEnergyFrozenFaceEnthalpy pressure_frozen_enthalpy{
      frozen.x, frozen.y, frozen.z, frozen.revision, frozen.reconstruction, 0U};
  pressure_frozen_enthalpy.local_binding =
      pressure_energy_frozen_face_enthalpy_local_binding(
          pressure_frozen_enthalpy);
  PressureEnergyPressureFluxBinding energy_pressure_binding;
  energy_pressure_binding.geometry = &fixture.geometry;
  energy_pressure_binding.boundary = &fixture.boundary;
  energy_pressure_binding.patch = fixture.patch;
  energy_pressure_binding.services =
      {MPI_COMM_WORLD, &pressure_halo, 8869U, 2U, 1U};
  energy_pressure_binding.intermediate = pressure_intermediate;
  energy_pressure_binding.pressure = pressure_certificate;
  energy_pressure_binding.temporal_diagonal =
      as_const(energy_pressure_temporal.view);
  energy_pressure_binding.x_pressure_coefficient =
      as_const(pressure_coefficients.x);
  energy_pressure_binding.y_pressure_coefficient =
      as_const(pressure_coefficients.y);
  energy_pressure_binding.z_pressure_coefficient =
      as_const(pressure_coefficients.z);
  energy_pressure_binding.target_flux = flux;
  energy_pressure_binding.frozen_face_enthalpy = pressure_frozen_enthalpy;
  energy_pressure_binding.activity = activity;
  energy_pressure_binding.identity = binding.identity;
  PressureEnergyPressureFluxOperator energy_pressure;
  PressureEnergyPressureFluxCertificate energy_pressure_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          energy_pressure_binding, energy_pressure,
          energy_pressure_certificate)) &&
          energy_pressure_certificate.valid() &&
          energy_pressure_certificate.pressure_work_scope ==
              PressureEnergyPressureWorkScope::flux_only_quasi_newton &&
          energy_pressure_certificate.flux_only_quasi_newton &&
          !energy_pressure_certificate.full_cartesian_pressure_work &&
          energy_pressure_certificate.activity_local_fingerprint ==
              certificate.activity_local_fingerprint &&
          energy_pressure_certificate.activity_collective_fingerprint ==
              certificate.activity_collective_fingerprint,
      "IBM E_p binds an activity-aware flux-only spatial response");

  PressureEnergySchurBlockAuthority spatial_authority;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurBlockAuthority::
                            ibm_cartesian_spatial_quasi_newton(
                                energy_pressure, operation,
                                spatial_authority)) &&
          spatial_authority.valid() &&
          spatial_authority.scope() == PressureEnergySchurBlockScope::
                                           ibm_cartesian_spatial_quasi_newton,
      "typed authority accepts only the masked Cartesian IBM spatial pair");

  LinearOperatorCertificate continuity_pressure_certificate =
      energy_pressure_certificate.linear;
  continuity_pressure_certificate.collective_fingerprint ^=
      UINT64_C(0x3131313131313131);
  if (continuity_pressure_certificate.collective_fingerprint == 0U)
    continuity_pressure_certificate.collective_fingerprint = 1U;
  continuity_pressure_certificate.operator_class = LinearOperatorClass::spd;
  IdentityCertificateOperator continuity_pressure(
      continuity_pressure_certificate);
  OwnedField continuity_enthalpy =
      shaped_field(38U, 8870U, 8871U, cells, -2.0);
  OwnedField continuity_enthalpy_scale =
      shaped_field(39U, 8872U, 8873U, cells, 2.0);
  OwnedField continuity_workspace =
      shaped_field(40U, 8874U, 8875U, cells, 0.0);
  OwnedField eliminated_workspace =
      shaped_field(41U, 8876U, 8877U, cells, 0.0);
  OwnedField energy_workspace =
      shaped_field(42U, 8878U, 8879U, cells, 0.0);
  PressureEnergySchurBinding schur_binding;
  schur_binding.continuity_pressure = &continuity_pressure;
  schur_binding.energy_pressure = &energy_pressure;
  schur_binding.energy_enthalpy = &operation;
  schur_binding.continuity_enthalpy_diagonal =
      as_const(continuity_enthalpy.view);
  schur_binding.continuity_enthalpy_row_scale =
      as_const(continuity_enthalpy_scale.view);
  schur_binding.workspace = {continuity_workspace.view,
                             eliminated_workspace.view,
                             energy_workspace.view};
  schur_binding.activity = {activity.cells, activity.local_fingerprint,
                            activity.collective_fingerprint};
  schur_binding.scaled_pivot_floor = 1.0e-12;
  schur_binding.block_authority = spatial_authority;
  PressureEnergySchurOperator spatial_schur;
  PressureEnergyJacobianCertificate spatial_schur_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurOperator::bind(
          schur_binding, spatial_schur, spatial_schur_certificate)) &&
          spatial_schur_certificate.valid() &&
          spatial_schur_certificate.jacobian_scope ==
              PressureEnergyJacobianScope::
                  ibm_cartesian_spatial_quasi_newton &&
          spatial_schur_certificate.exact_algebraic_schur &&
          !spatial_schur_certificate.full_nonlinear_jacobian &&
          !certificate.ibm_spatial_derivative,
      "IBM Schur certifies the algebraic spatial quasi-Newton block without "
      "claiming full nonlinear or IBM donor derivatives");
  if (!passed)
    return false;

  OwnedField output = shaped_field(35U, 8845U, 8846U, cells, -9.0);
  const Status apply_status = operation.apply(variation.view, output.view);
  passed &= expect(static_cast<bool>(apply_status),
                   "periodic MPI E_h applies with one dual-field exchange");
  if (!passed)
    return false;
  const HaloRuntimeCounters halo_counters = fixture.halo.runtime_counters();
  passed &= expect(halo_counters.begin_calls == 1U &&
                       halo_counters.finish_calls == 1U,
                   "dh(reach2) and deltaT(reach1) share one halo exchange");

  const double dx = fixture.geometry.x().uniform_width();
  const double dy = fixture.geometry.y().uniform_width();
  const double dz = fixture.geometry.z().uniform_width();
  const double volume = dx * dy * dz;
  const auto shift = [](Int3 cell, CartesianAxis axis, int side) {
    if (axis == CartesianAxis::x)
      cell.x += side;
    else if (axis == CartesianAxis::y)
      cell.y += side;
    else
      cell.z += side;
    return cell;
  };
  const auto face_view = [](ConstFaceFieldView x, ConstFaceFieldView y,
                            ConstFaceFieldView z, CartesianAxis axis) {
    return axis == CartesianAxis::x ? x : (axis == CartesianAxis::y ? y : z);
  };
  const auto active_face = [&](CartesianAxis axis, Int3 face) {
    if (axis == CartesianAxis::x)
      return active_x[offset(x_shape, face)] != 0U;
    if (axis == CartesianAxis::y)
      return active_y[offset(y_shape, face)] != 0U;
    return active_z[offset(z_shape, face)] != 0U;
  };
  const auto transmissibility = [&](ConstFieldView material, CartesianAxis axis,
                                    Int3 face) {
    const double left = material.unchecked(shift(face, axis, -1), 0U);
    const double right = material.unchecked(face, 0U);
    const double spacing =
        axis == CartesianAxis::x ? dx : (axis == CartesianAxis::y ? dy : dz);
    const double area = axis == CartesianAxis::x
                            ? dy * dz
                            : (axis == CartesianAxis::y ? dx * dz : dx * dy);
    return area * 2.0 * left * right / (spacing * (left + right));
  };
  double maximum_error = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (active_cells[offset(cells, cell)] == 0U) {
          maximum_error = std::max(
              maximum_error, std::abs(output.view.unchecked(cell, 0U) -
                                      variation.view.unchecked(cell, 0U)));
          continue;
        }
        double proxy_diagonal = 0.0;
        double convection = 0.0;
        double conduction = 0.0;
        for (CartesianAxis axis :
             {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
          const ConstFaceFieldView rates =
              face_view(flux.x, flux.y, flux.z, axis);
          for (int side : {-1, 1}) {
            const Int3 face = side < 0 ? cell : shift(cell, axis, 1);
            proxy_diagonal +=
                transmissibility(as_const(lambda_over_cp.view), axis, face);
            if (!active_face(axis, face))
              continue;
            const Int3 neighbour = shift(cell, axis, side);
            const double direction_face =
                0.5 * (variation.view.unchecked(cell, 0U) +
                       variation.view.unchecked(neighbour, 0U));
            convection += side * rates.unchecked(face) * direction_face;
            conduction += transmissibility(as_const(lambda.view), axis, face) *
                          (delta_temperature.view.unchecked(cell, 0U) -
                           delta_temperature.view.unchecked(neighbour, 0U));
          }
        }
        const double local = assembled.view.unchecked(cell, 0U) -
                             proxy_diagonal +
                             binding.authority.bdf.a0 * volume *
                                 target.view.unchecked(cell, 0U) *
                                 rho_h.view.unchecked(cell, 0U);
        const double expected = local * variation.view.unchecked(cell, 0U) +
                                convection + conduction;
        maximum_error =
            std::max(maximum_error,
                     std::abs(output.view.unchecked(cell, 0U) - expected) /
                         std::max(1.0, std::abs(expected)));
      }
    }
  }
  passed &= expect(maximum_error < 5.0e-13,
                   "periodic MPI response matches independent conservative "
                   "face and temperature-space oracle");

  OwnedField variation_from_solid = ghosted_field(5U, 8847U, 8848U, cells, 2U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        variation_from_solid.view.unchecked({x, y, z}, 0U) =
            variation.view.unchecked({x, y, z}, 0U);
  variation_from_solid.view.unchecked(solid, 0U) += 1000.0;
  OwnedField output_from_solid = shaped_field(36U, 8849U, 8850U, cells, -8.0);
  passed &= expect(static_cast<bool>(operation.apply(variation_from_solid.view,
                                                     output_from_solid.view)),
                   "inactive-cell perturbation applies");
  double active_change = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (active_cells[offset(cells, cell)] != 0U)
          active_change =
              std::max(active_change,
                       std::abs(output_from_solid.view.unchecked(cell, 0U) -
                                output.view.unchecked(cell, 0U)));
      }
  passed &= expect(active_change < 5.0e-13 &&
                       close(output_from_solid.view.unchecked(solid, 0U),
                             variation_from_solid.view.unchecked(solid, 0U)),
                   "inactive row is identity and all incident interfaces are "
                   "exactly zero");

  const std::vector<double> collective_rejected_snapshot =
      output_from_solid.storage;
  const double saved_rank_direction = variation.view.unchecked({0, 0, 0}, 0U);
  if (rank == 0)
    variation.view.unchecked({0, 0, 0}, 0U) =
        std::numeric_limits<double>::quiet_NaN();
  const Status collective_nonfinite =
      operation.apply(variation.view, output_from_solid.view);
  variation.view.unchecked({0, 0, 0}, 0U) = saved_rank_direction;
  passed &= expect(
      collective_nonfinite.code == StatusCode::numerical_failure &&
          output_from_solid.storage == collective_rejected_snapshot,
      "one-rank nonfinite dh rejects collectively without a partial output "
      "commit");

  int local_pass = passed ? 1 : 0;
  int global_pass = 0;
  MPI_Allreduce(&local_pass, &global_pass, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  return global_pass == 1;
}

bool test_enthalpy_semismooth_limiter_certificate() {
  EnthalpySpatialFixture fixture;
  bool passed =
      expect(make_enthalpy_spatial_fixture(fixture, MPI_COMM_SELF, false, true),
             "limited-convection E_h fixture compiles");
  if (!passed)
    return false;
  const Int3 cells = fixture.patch.cells;
  OwnedField assembled = shaped_field(40U, 8901U, 8902U, cells, 25.0);
  OwnedField target = ghosted_field(3U, 8903U, 8904U, cells, 2U, 0.0);
  OwnedField rho_h = shaped_field(41U, 8905U, 8906U, cells, -1.0e-6);
  OwnedField cp = ghosted_field(42U, 8907U, 8908U, cells, 1U, 1000.0);
  OwnedField lambda = ghosted_field(43U, 8909U, 8910U, cells, 1U, 2.0);
  OwnedField diffusivity = ghosted_field(44U, 8911U, 8912U, cells, 1U, 0.002);
  OwnedField delta_temperature = ghosted_field(6U, 8913U, 8914U, cells, 1U);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        target.view.unchecked({x, y, z}, 0U) = 300000.0;
  passed &= expect(
      static_cast<bool>(apply_homogeneous_scalar_boundary_ghosts(
          BoundaryStage::enthalpy, fixture.boundary, 3U, target.view, 2U)),
      "constant limiter target closes physical ghosts");
  std::array<FieldView, 3U> material{cp.view, lambda.view, diffusivity.view};
  passed &= expect(static_cast<bool>(apply_physical_zero_gradient(
                       fixture.boundary, {material.data(), material.size()})),
                   "constant limiter material ghosts close");
  OwnedFaces flux_storage = face_bundle(cells, 8915U, 8916U, 0.125);
  const ConstFaceFluxView flux = flux_view(flux_storage, 8917U);
  OwnedFaces frozen_storage = face_bundle(cells, 8918U, 8919U);
  OwnedFaces directional_storage = face_bundle(cells, 8920U, 8921U);
  OwnedFaces compiled_conductance = face_bundle(cells, 8938U, 8939U);
  OwnedField compiled_local = shaped_field(46U, 8940U, 8941U, cells, 0.0);
  OwnedField compiled_response = shaped_field(47U, 8942U, 8943U, cells, 0.0);
  const auto face_values = [](ConstFaceFieldView face) {
    return static_cast<std::size_t>(face.extents.x) *
           static_cast<std::size_t>(face.extents.y) *
           static_cast<std::size_t>(face.extents.z);
  };
  const std::size_t branch_count =
      face_values(as_const(flux_storage.x)) +
      face_values(as_const(flux_storage.y)) +
      face_values(as_const(flux_storage.z));
  std::vector<std::uint16_t> compiled_branches(branch_count, UINT16_C(0));
  const FrozenConvectionContext context{8922U, fixture.boundary.revision()};
  FrozenConvectionFaceField frozen;
  passed &= expect(prepare_frozen_enthalpy(
                       fixture, flux, as_const(target.view),
                       {frozen_storage.x, frozen_storage.y, frozen_storage.z},
                       context, frozen),
                   "constant target freezes at limiter kinks");
  if (!passed)
    return false;

  PressureEnergyEnthalpyBinding binding;
  binding.geometry = &fixture.geometry;
  binding.kernels = &fixture.kernels;
  binding.boundary = &fixture.boundary;
  binding.patch = fixture.patch;
  binding.convection = fixture.schemes.enthalpy();
  binding.services = {MPI_COMM_SELF, &fixture.halo, 8923U, 5U, 6U};
  binding.authority = {
      {2.0, -2.0, 0.0, 1U},         8924U, fixture.geometry.topology_revision(),
      fixture.boundary.revision(),  8925U, 8926U,
      context.collective_semantics, 8927U, 8928U};
  binding.assembled_diagonal = as_const(assembled.view);
  binding.target_enthalpy = as_const(target.view);
  binding.density_enthalpy_derivative = as_const(rho_h.view);
  binding.heat_capacity = as_const(cp.view);
  binding.thermal_conductivity = as_const(lambda.view);
  binding.enthalpy_diffusivity = as_const(diffusivity.view);
  binding.target_flux = flux;
  binding.convection_context = context;
  binding.frozen_face_enthalpy = frozen;
  binding.workspace = {
      delta_temperature.view,
      {directional_storage.x, directional_storage.y, directional_storage.z},
      {compiled_local.view,
       {compiled_conductance.x, compiled_conductance.y,
        compiled_conductance.z},
       {{compiled_branches.data(), compiled_branches.size()}},
       compiled_response.view}};
  binding.identity = {8929U, 8930U, 8931U, 8932U, 8933U};
  binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope;
  PressureEnergyEnthalpyOperator semismooth;
  PressureEnergyEnthalpyCertificate semismooth_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyEnthalpyOperator::bind(
          binding, semismooth, semismooth_certificate)) &&
          semismooth_certificate.valid() &&
          semismooth_certificate.compiled_factored_apply &&
          semismooth_certificate.compiled_numeric_revision != 0U &&
          semismooth_certificate.compiled_local_binding != 0U &&
          semismooth_certificate.generalized_face_count > 0U &&
          semismooth_certificate.linearization_policy ==
              FrozenConvectionLinearizationPolicy::
                  semismooth_generalized_zero_slope,
      "semismooth zero-slope policy certifies every encountered limiter kink");

  PressureEnergyEnthalpyBinding classical_binding = binding;
  classical_binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::classical_active_branch;
  PressureEnergyEnthalpyOperator classical;
  PressureEnergyEnthalpyCertificate classical_certificate;
  const Status classical_status = PressureEnergyEnthalpyOperator::bind(
      classical_binding, classical, classical_certificate);
  passed &= expect(classical_status.code == StatusCode::numerical_failure &&
                       !classical_certificate.valid(),
                   "classical policy refuses nondifferentiable limiter kinks");
  if (!passed)
    return false;

  OwnedField direction = ghosted_field(5U, 8934U, 8935U, cells, 2U, 0.0);
  OwnedField output = shaped_field(45U, 8936U, 8937U, cells, -1.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        direction.view.unchecked({x, y, z}, 0U) =
            1.0 + 0.1 * x - 0.05 * y + 0.03 * z;
  passed &=
      expect(static_cast<bool>(semismooth.apply(direction.view, output.view)),
             "certified semismooth limiter action remains finite");
  PressureEnergyEnthalpyPreparedEpoch prepared_epoch;
  OwnedField prepared_output = shaped_field(48U, 8944U, 8945U, cells, -2.0);
  const Status prepare_status =
      semismooth.prepare_repeated_apply(prepared_epoch);
  const Status prepared_status = semismooth.apply_prepared(
      direction.view, prepared_output.view, prepared_epoch);
  const Status close_status = semismooth.close_repeated_apply(prepared_epoch);
  bool prepared_equal = true;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        prepared_equal =
            prepared_equal &&
            close(prepared_output.view.unchecked({x, y, z}, 0U),
                  output.view.unchecked({x, y, z}, 0U));
  passed &= expect(static_cast<bool>(prepare_status) &&
                       static_cast<bool>(prepared_status) &&
                       static_cast<bool>(close_status) && prepared_equal &&
                       !prepared_epoch.valid(),
                   "prepared limited E_h epoch preserves the exact action");

  std::fill(prepared_output.storage.begin(), prepared_output.storage.end(),
            -3.0);
  const std::vector<double> rejected_snapshot = prepared_output.storage;
  compiled_branches[0U] ^= UINT16_C(1);
  const Status stale_branch =
      semismooth.apply(direction.view, prepared_output.view);
  compiled_branches[0U] ^= UINT16_C(1);
  passed &= expect(stale_branch.code == StatusCode::invalid_plan &&
                       prepared_output.storage == rejected_snapshot,
                   "raw compiled branch mutation fails atomically");

  const double saved_target = target.view.unchecked({1, 1, 1}, 0U);
  target.view.unchecked({1, 1, 1}, 0U) = saved_target + 0.25;
  const Status stale_target =
      semismooth.apply(direction.view, prepared_output.view);
  target.view.unchecked({1, 1, 1}, 0U) = saved_target;
  passed &= expect(stale_target.code == StatusCode::invalid_plan &&
                       prepared_output.storage == rejected_snapshot,
                   "compiled limited E_h detects raw target mutation before "
                   "output commit");

  const double saved_lambda = lambda.view.unchecked({1, 1, 1}, 0U);
  lambda.view.unchecked({1, 1, 1}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  const Status stale_transport =
      semismooth.apply(direction.view, prepared_output.view);
  lambda.view.unchecked({1, 1, 1}, 0U) = saved_lambda;
  passed &= expect(stale_transport.code == StatusCode::numerical_failure &&
                       prepared_output.storage == rejected_snapshot,
                   "compiled limited E_h detects nonfinite raw transport "
                   "before output commit");
  return passed;
}

bool test_diagonal_operator_activity_and_identity() {
  const LinearIdentity identity{501U, 502U, 503U, 504U, 505U};
  OwnedField diagonal = field(31U, 41U, 1041U,
                              {{2.0, std::numeric_limits<double>::quiet_NaN(),
                                -4.0}});
  const std::array<std::uint8_t, kSize> active{{1U, 0U, 1U}};
  PressureEnergyCellActivity activity{
      {active.data(), active.size()}, 601U, 602U};
  PressureEnergyDiagonalBinding binding;
  binding.diagonal = as_const(diagonal.view);
  binding.activity = activity;
  binding.identity = identity;
  binding.inactive_diagonal = 0.0;
  PressureEnergyDiagonalOperator op;
  PressureEnergyDiagonalCertificate certificate;
  bool passed = expect(static_cast<bool>(PressureEnergyDiagonalOperator::bind(
                           binding, op, certificate)),
                       "cell-local Ep diagonal binds with an IBM activity row");
  passed &= expect(certificate.valid() &&
                       certificate.linear.identity.fingerprint ==
                           identity.fingerprint &&
                       certificate.diagonal_revision == diagonal.view.revision &&
                       certificate.active_cells == 2U &&
                       certificate.inactive_cells == 1U,
                   "diagonal certificate binds identity, revision, and mask");
  OwnedField input = field(32U, 42U, 1042U, {{3.0, 7.0, -2.0}});
  OwnedField output = field(33U, 43U, 1043U);
  passed &= expect(static_cast<bool>(op.apply(input.view, output.view)),
                   "masked cell-local diagonal applies");
  const Vector result = values(as_const(output.view));
  passed &= expect(close(result[0U], 6.0) && close(result[1U], 0.0) &&
                       close(result[2U], 8.0),
                   "inactive Ep row is explicit zero and ignores a NaN payload");

  binding.inactive_diagonal = 1.0;
  PressureEnergyDiagonalOperator eh;
  PressureEnergyDiagonalCertificate eh_certificate;
  passed &= expect(static_cast<bool>(PressureEnergyDiagonalOperator::bind(
                           binding, eh, eh_certificate)) &&
                       static_cast<bool>(eh.apply(input.view, output.view)) &&
                       close(values(as_const(output.view))[1U], 7.0),
                   "inactive Eh row can be an explicit identity row");
  return passed;
}

bool test_ibm_double_diagonal_typed_schur_authority() {
  const LinearIdentity identity{9101U, 9102U, 9103U, 9104U, 9105U};
  const std::array<std::uint8_t, kSize> active{{1U, 0U, 1U}};
  const PressureEnergyCellActivity activity{
      {active.data(), active.size()}, 9106U, 9107U};
  OwnedField energy_pressure_diagonal =
      field(50U, 9108U, 9109U,
            {{2.0, std::numeric_limits<double>::quiet_NaN(), -4.0}});
  OwnedField energy_enthalpy_diagonal =
      field(51U, 9110U, 9111U,
            {{3.0, std::numeric_limits<double>::quiet_NaN(), 5.0}});

  PressureEnergyDiagonalBinding energy_pressure_binding;
  energy_pressure_binding.diagonal = as_const(energy_pressure_diagonal.view);
  energy_pressure_binding.activity = activity;
  energy_pressure_binding.identity = identity;
  energy_pressure_binding.inactive_diagonal = 0.0;
  PressureEnergyDiagonalOperator energy_pressure;
  PressureEnergyDiagonalCertificate energy_pressure_certificate;
  PressureEnergyDiagonalBinding energy_enthalpy_binding;
  energy_enthalpy_binding.diagonal = as_const(energy_enthalpy_diagonal.view);
  energy_enthalpy_binding.activity = activity;
  energy_enthalpy_binding.identity = identity;
  energy_enthalpy_binding.inactive_diagonal = 1.0;
  PressureEnergyDiagonalOperator energy_enthalpy;
  PressureEnergyDiagonalCertificate energy_enthalpy_certificate;
  bool passed = expect(
      static_cast<bool>(PressureEnergyDiagonalOperator::bind(
          energy_pressure_binding, energy_pressure,
          energy_pressure_certificate)) &&
          static_cast<bool>(PressureEnergyDiagonalOperator::bind(
              energy_enthalpy_binding, energy_enthalpy,
              energy_enthalpy_certificate)) &&
          energy_pressure_certificate.inactive_diagonal == 0.0 &&
          energy_enthalpy_certificate.inactive_diagonal == 1.0,
      "IBM fallback binds typed Ep=zero and Eh=identity inactive rows");
  if (!passed)
    return false;

  PressureEnergySchurBlockAuthority ibm_authority;
  passed &= expect(
      static_cast<bool>(
          PressureEnergySchurBlockAuthority::ibm_double_diagonal(
              energy_pressure, energy_enthalpy, ibm_authority)) &&
          ibm_authority.valid() &&
          ibm_authority.scope() == PressureEnergySchurBlockScope::
                                       ibm_double_diagonal_quasi_newton,
      "the Ep/Eh double diagonal pair issues IBM quasi-Newton authority");

  PressureEnergySchurBlockAuthority preserved_authority = ibm_authority;
  passed &= expect(
      PressureEnergySchurBlockAuthority::ibm_double_diagonal(
          energy_enthalpy, energy_pressure, preserved_authority)
                  .code == StatusCode::invalid_plan &&
          preserved_authority.valid() &&
          preserved_authority.scope() == PressureEnergySchurBlockScope::
                                             ibm_double_diagonal_quasi_newton,
      "swapped inactive diagonals fail closed and preserve prior authority");
  PressureEnergyDiagonalBinding foreign_activity_diagonal_binding =
      energy_enthalpy_binding;
  ++foreign_activity_diagonal_binding.activity.local_fingerprint;
  PressureEnergyDiagonalOperator foreign_activity_enthalpy;
  PressureEnergyDiagonalCertificate foreign_activity_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyDiagonalOperator::bind(
          foreign_activity_diagonal_binding, foreign_activity_enthalpy,
          foreign_activity_certificate)) &&
          PressureEnergySchurBlockAuthority::ibm_double_diagonal(
              energy_pressure, foreign_activity_enthalpy,
              preserved_authority)
                  .code == StatusCode::invalid_plan,
      "IBM double-diagonal authority rejects mismatched activity lineage");

  LinearOperatorCertificate continuity_pressure_certificate =
      energy_pressure_certificate.linear;
  continuity_pressure_certificate.collective_fingerprint ^=
      UINT64_C(0x7171717171717171);
  if (continuity_pressure_certificate.collective_fingerprint == 0U)
    continuity_pressure_certificate.collective_fingerprint = 1U;
  continuity_pressure_certificate.operator_class = LinearOperatorClass::spd;
  IdentityCertificateOperator continuity_pressure(
      continuity_pressure_certificate);
  OwnedField continuity_enthalpy =
      field(52U, 9112U, 9113U,
            {{-2.0, std::numeric_limits<double>::quiet_NaN(), -4.0}});
  OwnedField continuity_enthalpy_scale =
      field(53U, 9114U, 9115U,
            {{2.0, std::numeric_limits<double>::quiet_NaN(), 4.0}});
  OwnedField continuity_workspace = field(54U, 9116U, 9117U);
  OwnedField eliminated_workspace = field(55U, 9118U, 9119U);
  OwnedField energy_workspace = field(56U, 9120U, 9121U);
  PressureEnergySchurBinding binding;
  binding.continuity_pressure = &continuity_pressure;
  binding.energy_pressure = &energy_pressure;
  binding.energy_enthalpy = &energy_enthalpy;
  binding.continuity_enthalpy_diagonal =
      as_const(continuity_enthalpy.view);
  binding.continuity_enthalpy_row_scale =
      as_const(continuity_enthalpy_scale.view);
  binding.workspace = {continuity_workspace.view,
                       eliminated_workspace.view, energy_workspace.view};
  binding.activity = activity;
  binding.scaled_pivot_floor = 1.0e-12;

  PressureEnergySchurOperator missing_authority_operator;
  PressureEnergyJacobianCertificate missing_authority_certificate;
  passed &= expect(
      PressureEnergySchurOperator::bind(
          binding, missing_authority_operator, missing_authority_certificate)
                  .code == StatusCode::invalid_plan &&
          !missing_authority_certificate.valid(),
      "IBM Schur rejects a bare binding without typed diagonal authority");

  binding.block_authority = ibm_authority;
  PressureEnergySchurOperator schur;
  PressureEnergyJacobianCertificate schur_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurOperator::bind(
          binding, schur, schur_certificate)) &&
          schur_certificate.valid() &&
          schur_certificate.jacobian_scope ==
              PressureEnergyJacobianScope::
                  ibm_double_diagonal_quasi_newton &&
          schur_certificate.active_cells == 2U &&
          schur_certificate.inactive_cells == 1U &&
          !schur_certificate.full_nonlinear_jacobian,
      "Schur derives IBM double-diagonal scope without a nonlinear-Jacobian "
      "claim");

  OwnedField direction = field(57U, 9122U, 9123U, {{0.5, 7.0, -0.25}});
  OwnedField output = field(58U, 9124U, 9125U, {{-9.0, -9.0, -9.0}});
  passed &= expect(static_cast<bool>(schur.apply(direction.view, output.view)) &&
                       close(output.view.unchecked({1, 0, 0}, 0U), 7.0),
                   "IBM Schur preserves the inactive pressure identity row");

  PressureEnergySchurBinding foreign_activity_binding = binding;
  ++foreign_activity_binding.activity.local_fingerprint;
  PressureEnergySchurOperator foreign_activity_operator;
  PressureEnergyJacobianCertificate foreign_activity_schur_certificate =
      schur_certificate;
  passed &= expect(
      PressureEnergySchurOperator::bind(
          foreign_activity_binding, foreign_activity_operator,
          foreign_activity_schur_certificate)
                  .code == StatusCode::invalid_plan &&
          !foreign_activity_schur_certificate.valid(),
      "typed IBM authority cannot be reused with a foreign activity map");
  return passed;
}

bool test_mass_flow_three_cell_pressure_flux_red() {
  PressureFluxFixture fixture;
  bool passed = expect(make_pressure_flux_fixture(
                           PressureFluxBoundary::neumann_dirichlet, fixture),
                       "three-cell pressure-flux RED fixture compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  OwnedField temporal = shaped_field(80U, 8301U, 8302U, cells, 0.0);
  OwnedFaces coefficient = face_bundle(cells, 8303U, 8304U, 0.0);
  coefficient.x.unchecked({1, 0, 0}) = 1.0;
  coefficient.x.unchecked({2, 0, 0}) = 1.0;
  OwnedFaces target = face_bundle(cells, 8305U, 8306U, 0.0);
  OwnedFaces enthalpy = face_bundle(cells, 8307U, 8308U, 316470.0);
  constexpr double delta_phi_1 = -3.23719569791e-8;
  constexpr double delta_phi_2 = +1.61855552962e-8;
  constexpr double terminal_0 = -0.01024416458;
  constexpr double terminal_1 = +0.01536632105;
  constexpr double terminal_2 = -0.00512216550;
  enthalpy.x.unchecked({1, 0, 0}) = -terminal_0 / -delta_phi_1;
  enthalpy.x.unchecked({2, 0, 0}) = -terminal_2 / delta_phi_2;
  OwnedField dp = ghosted_field(2U, 8309U, 8310U, cells, 1U);
  dp.view.unchecked({0, 0, 0}, 0U) = 0.0;
  dp.view.unchecked({1, 0, 0}, 0U) = -delta_phi_1;
  dp.view.unchecked({2, 0, 0}, 0U) = -delta_phi_1 - delta_phi_2;
  OwnedField output = shaped_field(81U, 8311U, 8312U, cells, 0.0);
  const PressureEnergyPressureFluxBinding binding = pressure_flux_binding(
      fixture, as_const(temporal.view), coefficient, target, enthalpy);
  PressureEnergyPressureFluxOperator op;
  PressureEnergyPressureFluxCertificate certificate;
  passed &= expect(
      static_cast<bool>(
          PressureEnergyPressureFluxOperator::bind(binding, op, certificate)) &&
          certificate.valid() && certificate.conservative_face_response &&
          certificate.exact_piso_flux_jump &&
          certificate.target_flux_revision == binding.target_flux.revision &&
          certificate.pressure_face_coefficients ==
              binding.intermediate.pressure_face_coefficient,
      "mass-flow E_p binds the pressure faces and frozen target h_f");
  passed &= expect(static_cast<bool>(op.apply(dp.view, output.view)),
                   "mass-flow E_p applies");
  const Vector response = values(as_const(output.view));
  const double predicted_1 = -terminal_0 - terminal_2;
  passed &=
      expect(close(response[0U], terminal_0, 2.0e-13) &&
                 close(response[1U], predicted_1, 2.0e-13) &&
                 close(response[2U], terminal_2, 2.0e-13),
             "D(h_f dphi) reproduces the conservative three-cell pattern");
  passed &= expect(
      std::abs(response[0U] - terminal_0) < 1.0e-12 &&
          std::abs(response[1U] - terminal_1) < 1.0e-8 &&
          std::abs(response[2U] - terminal_2) < 1.0e-12,
      "pressure-flux response explains over 99.99 percent of terminal RED");
  const double face_one = pressure_correction_mass_flux_response(
      as_const(dp.view), fixture.geometry, fixture.patch, fixture.boundary,
      CartesianAxis::x, {1, 0, 0}, 1.0);
  const double face_two = pressure_correction_mass_flux_response(
      as_const(dp.view), fixture.geometry, fixture.patch, fixture.boundary,
      CartesianAxis::x, {2, 0, 0}, 1.0);
  passed &= expect(close(face_one, delta_phi_1) && close(face_two, delta_phi_2),
                   "RED uses the measured internal pressure-flux increments");
  return passed;
}

bool test_boundary_constant_h_and_directional_derivative() {
  PressureFluxFixture fixture;
  bool passed = expect(make_pressure_flux_fixture(
                           PressureFluxBoundary::neumann_dirichlet, fixture),
                       "physical-boundary E_p fixture compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  OwnedField temporal = shaped_field(82U, 8401U, 8402U, cells, 0.0);
  OwnedFaces coefficient = face_bundle(cells, 8403U, 8404U, 0.0);
  coefficient.x.unchecked({0, 0, 0}) = 2.0;
  coefficient.x.unchecked({1, 0, 0}) = 1.0;
  coefficient.x.unchecked({2, 0, 0}) = 1.0;
  coefficient.x.unchecked({3, 0, 0}) = 2.0;
  OwnedFaces target = face_bundle(cells, 8405U, 8406U, 0.0);
  target.x.unchecked({0, 0, 0}) = 0.3;
  target.x.unchecked({1, 0, 0}) = -0.2;
  target.x.unchecked({2, 0, 0}) = 0.4;
  target.x.unchecked({3, 0, 0}) = 0.1;
  OwnedFaces unit_h = face_bundle(cells, 8407U, 8408U, 1.0);
  OwnedFaces ten_h = face_bundle(cells, 8409U, 8410U, 10.0);
  OwnedField dp = ghosted_field(2U, 8411U, 8412U, cells, 1U);
  dp.view.unchecked({0, 0, 0}, 0U) = 1.0;
  dp.view.unchecked({1, 0, 0}, 0U) = 2.0;
  dp.view.unchecked({2, 0, 0}, 0U) = 4.0;
  OwnedField unit_output = shaped_field(83U, 8413U, 8414U, cells, 0.0);
  OwnedField ten_output = shaped_field(84U, 8415U, 8416U, cells, 0.0);
  PressureEnergyPressureFluxOperator unit_op;
  PressureEnergyPressureFluxOperator ten_op;
  PressureEnergyPressureFluxCertificate unit_certificate;
  PressureEnergyPressureFluxCertificate ten_certificate;
  const auto unit_binding = pressure_flux_binding(
      fixture, as_const(temporal.view), coefficient, target, unit_h);
  const auto ten_binding = pressure_flux_binding(
      fixture, as_const(temporal.view), coefficient, target, ten_h);
  passed &=
      expect(static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
                 unit_binding, unit_op, unit_certificate)) &&
                 static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
                     ten_binding, ten_op, ten_certificate)) &&
                 static_cast<bool>(unit_op.apply(dp.view, unit_output.view)) &&
                 static_cast<bool>(ten_op.apply(dp.view, ten_output.view)),
             "constant-h face response operators bind and apply");
  const Vector unit = values(as_const(unit_output.view));
  const Vector ten = values(as_const(ten_output.view));
  passed &= expect(
      close(unit[0U], -1.0) && close(unit[1U], -1.0) && close(unit[2U], 10.0),
      "physical Neumann and Dirichlet jumps have exact signs");
  for (std::size_t cell = 0U; cell < kSize; ++cell) {
    passed &= expect(close(ten[cell], 10.0 * unit[cell]),
                     "D(h dphi)=h D(dphi) for constant frozen h");
  }
  passed &= expect(
      close(pressure_correction_jump(as_const(dp.view), fixture.geometry,
                                     fixture.patch, fixture.boundary,
                                     CartesianAxis::x, {0, 0, 0}),
            0.0) &&
          close(pressure_correction_jump(as_const(dp.view), fixture.geometry,
                                         fixture.patch, fixture.boundary,
                                         CartesianAxis::x, {3, 0, 0}),
                -4.0),
      "single jump authority selects physical Neumann and Dirichlet closure");

  temporal.view.unchecked({0, 0, 0}, 0U) = 2.0;
  temporal.view.unchecked({1, 0, 0}, 0U) = 3.0;
  temporal.view.unchecked({2, 0, 0}, 0U) = 4.0;
  ++temporal.view.revision;
  dp.view.unchecked({0, 0, 0}, 0U) = 0.7;
  dp.view.unchecked({1, 0, 0}, 0U) = -0.3;
  dp.view.unchecked({2, 0, 0}, 0U) = 0.4;
  OwnedField derivative = shaped_field(85U, 8417U, 8418U, cells, 0.0);
  PressureEnergyPressureFluxOperator derivative_op;
  PressureEnergyPressureFluxCertificate derivative_certificate;
  const auto derivative_binding = pressure_flux_binding(
      fixture, as_const(temporal.view), coefficient, target, ten_h);
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          derivative_binding, derivative_op, derivative_certificate)) &&
          static_cast<bool>(derivative_op.apply(dp.view, derivative.view)),
      "directional-derivative E_p binds and applies");
  auto advanced_boundary_target = derivative_binding;
  ++advanced_boundary_target.pressure.numeric_boundary;
  PressureEnergyPressureFluxOperator advanced_boundary_op;
  PressureEnergyPressureFluxCertificate advanced_boundary_certificate;
  passed &= expect(
      !PressureEnergyPressureFluxOperator::bind(
          advanced_boundary_target, advanced_boundary_op,
          advanced_boundary_certificate) &&
          !advanced_boundary_certificate.valid(),
      "foreign numeric boundary fails closed against the compiled plan");
  OwnedField base = ghosted_field(2U, 8419U, 8420U, cells, 1U);
  OwnedField shifted = ghosted_field(2U, 8421U, 8422U, cells, 1U);
  constexpr double epsilon = 1.0e-6;
  for (std::int32_t x = 0; x < cells.x; ++x) {
    const Int3 cell{x, 0, 0};
    const double value = 1.2 + 0.2 * x;
    base.view.unchecked(cell, 0U) = value;
    shifted.view.unchecked(cell, 0U) =
        value + epsilon * dp.view.unchecked(cell, 0U);
  }
  const auto residual = [&](ConstFieldView pressure, Int3 cell) {
    const Int3 plus{cell.x + 1, 0, 0};
    const double minus_flux =
        target.x.unchecked(cell) +
        pressure_correction_mass_flux_response(
            pressure, fixture.geometry, fixture.patch, fixture.boundary,
            CartesianAxis::x, cell, coefficient.x.unchecked(cell));
    const double plus_flux =
        target.x.unchecked(plus) +
        pressure_correction_mass_flux_response(
            pressure, fixture.geometry, fixture.patch, fixture.boundary,
            CartesianAxis::x, plus, coefficient.x.unchecked(plus));
    return temporal.view.unchecked(cell, 0U) * pressure.unchecked(cell, 0U) +
           ten_h.x.unchecked(plus) * plus_flux -
           ten_h.x.unchecked(cell) * minus_flux;
  };
  const Vector applied = values(as_const(derivative.view));
  for (std::int32_t x = 0; x < cells.x; ++x) {
    const Int3 cell{x, 0, 0};
    const double finite_difference = (residual(as_const(shifted.view), cell) -
                                      residual(as_const(base.view), cell)) /
                                     epsilon;
    passed &= expect(
        close(applied[static_cast<std::size_t>(x)], finite_difference, 3.0e-8),
        "E_p matches the target residual directional derivative");
  }

  FieldView foreign = dp.view;
  foreign.field = 3U;
  passed &=
      expect(derivative_op.apply(foreign, derivative.view).code ==
                     StatusCode::invalid_plan &&
                 derivative_op.apply(dp.view, dp.view).code ==
                     StatusCode::invalid_plan,
             "pressure-flux E_p rejects foreign and aliased Krylov views");
  auto stale = derivative_binding;
  ++stale.pressure.geometry;
  auto foreign_frozen = derivative_binding;
  foreign_frozen.frozen_face_enthalpy.x = as_const(unit_h.x);
  auto foreign_thermophysical_local = derivative_binding;
  ++foreign_thermophysical_local.pressure
        .thermophysical_boundary_rank_local_binding;
  auto foreign_thermophysical_target = derivative_binding;
  ++foreign_thermophysical_target.pressure
        .thermophysical_boundary_target;
  PressureEnergyPressureFluxOperator rejected;
  PressureEnergyPressureFluxCertificate rejected_certificate;
  passed &=
      expect(PressureEnergyPressureFluxOperator::bind(stale, rejected,
                                                      rejected_certificate)
                         .code == StatusCode::invalid_plan &&
                 PressureEnergyPressureFluxOperator::bind(
                     foreign_frozen, rejected, rejected_certificate)
                         .code == StatusCode::invalid_plan &&
                 PressureEnergyPressureFluxOperator::bind(
                     foreign_thermophysical_local, rejected,
                     rejected_certificate)
                         .code == StatusCode::invalid_plan &&
                 PressureEnergyPressureFluxOperator::bind(
                     foreign_thermophysical_target, rejected,
                     rejected_certificate)
                         .code == StatusCode::invalid_plan &&
                 !rejected_certificate.valid(),
             "stale geometry, frozen-h, and thermophysical authorities "
             "cannot bind");
  return passed;
}

// Analytic oracle for one frozen Cartesian target layer.  This deliberately
// does not refresh EOS, transport, boundary coefficients, HbyA, rAU, or phi;
// it must never be cited as a full-refresh/nonlinear Jacobian certificate.
bool test_analytic_frozen_target_cartesian_four_block_fd_certificate() {
  EnthalpySpatialFixture fixture;
  bool passed = expect(
      make_enthalpy_spatial_fixture(fixture, MPI_COMM_SELF, true, false),
      "analytic frozen-target fixture compiles a periodic Cartesian layer");
  if (!passed) return false;

  SpeciesThermophysicalSpec air;
  air.stable_name = "air";
  air.molecular_weight = 28.96546;
  air.temperature_switch = 1000.0;
  air.nasa7_low[0U] = 3.5;
  air.nasa7_high[0U] = 3.5;
  air.viscosity_reference = 1.8e-5;
  air.conductivity = 0.026;
  ThermophysicalSpec thermophysical;
  thermophysical.data_file = "analytic-frozen-target-four-block.d";
  thermophysical.minimum_temperature = 200.0;
  thermophysical.maximum_temperature = 2000.0;
  thermophysical.temperature_relative_tolerance = 1.0e-12;
  thermophysical.maximum_temperature_iterations = 64U;
  thermophysical.closed_mass_relative_tolerance = 1.0e-12;
  thermophysical.maximum_closed_mass_iterations = 32U;
  thermophysical.maximum_closed_mass_relative_step = 0.2;
  thermophysical.species.push_back(air);
  ThermodynamicsPlan thermodynamics;
  TransportPlan transport;
  ContributionRegistry contributions;
  EquationPlanSet equations;
  const std::array<FieldId, 10U> declared{{0U, 1U, 2U, 3U, 4U,
                                           5U, 6U, 7U, 8U, 9U}};
  EquationPlanSpec equation_spec;
  equation_spec.density = 0U;
  equation_spec.velocity = 1U;
  equation_spec.pressure_perturbation = 2U;
  equation_spec.enthalpy = 3U;
  equation_spec.temperature = 4U;
  equation_spec.effective_viscosity = 7U;
  equation_spec.pressure_compressibility = 8U;
  equation_spec.velocity_gradient = 9U;
  equation_spec.pressure_reference = PressureReferenceKind::closed_mass;
  equation_spec.closed_mass_service_stage = 1U;
  equation_spec.maximum_cells_per_rank = 64U;
  passed &= expect(
      static_cast<bool>(ThermodynamicsPlan::compile(thermophysical, {},
                                                    thermodynamics)) &&
          static_cast<bool>(TransportPlan::compile(
              thermophysical, thermodynamics, transport)) &&
          static_cast<bool>(contributions.configure(
              {declared.data(), declared.size()})) &&
          static_cast<bool>(contributions.freeze()) &&
          static_cast<bool>(EquationPlanSet::compile(
              MPI_COMM_SELF, fixture.schemes, fixture.geometry,
              fixture.patch, fixture.boundary, contributions, thermodynamics,
              transport, equation_spec, equations)),
      "analytic frozen-target fixture compiles real pressure/enthalpy plans");
  if (!passed) return false;

  PisoPlanSpec piso_spec;
  piso_spec.pressure_correctors = 2U;
  piso_spec.pressure_stage = 910U;
  piso_spec.final_flux_slot = 0U;
  piso_spec.pressure_solve = {1.0e-15, 1.0e-13, 400U, 4U, 12U};
  piso_spec.eos_tolerance = 1.0e-10;
  piso_spec.continuity_tolerance = 1.0e-10;
  piso_spec.energy_tolerance = 1.0e-10;
  piso_spec.closed_mass_tolerance = 1.0e-10;
  piso_spec.gauge_tolerance = 1.0e-12;
  PisoPlan piso;
  passed &= expect(static_cast<bool>(PisoPlan::compile(
                       MPI_COMM_SELF, equations, piso_spec, piso)),
                   "analytic frozen-target fixture compiles the real PISO plan");
  if (!passed) return false;

  const Int3 cells = fixture.patch.cells;
  const std::uint8_t reach = equations.kernels().reach();
  constexpr double pressure_reference = 101325.0;
  constexpr double ideal_density_factor = 3.5;
  const BdfCoefficients bdf{2.0, -3.0, 1.0, 2U};
  const double step_dt =
      (1.0 + 1.0 / (std::sqrt(-bdf.a1 / bdf.a2) - 1.0)) / -bdf.a1;
  auto wrap = [](std::int32_t value, std::int32_t extent) noexcept {
    const std::int32_t remainder = value % extent;
    return remainder < 0 ? remainder + extent : remainder;
  };
  auto wrap_cell = [&](Int3 cell) noexcept {
    cell.x = wrap(cell.x, cells.x);
    cell.y = wrap(cell.y, cells.y);
    cell.z = wrap(cell.z, cells.z);
    return cell;
  };
  const auto fill_periodic_scalar = [&](FieldView view,
                                        std::uint8_t ghosts) {
    for (std::int32_t z = -static_cast<std::int32_t>(ghosts);
         z < cells.z + static_cast<std::int32_t>(ghosts); ++z) {
      for (std::int32_t y = -static_cast<std::int32_t>(ghosts);
           y < cells.y + static_cast<std::int32_t>(ghosts); ++y) {
        for (std::int32_t x = -static_cast<std::int32_t>(ghosts);
             x < cells.x + static_cast<std::int32_t>(ghosts); ++x) {
          const Int3 destination{x, y, z};
          view.unchecked(destination, 0U) =
              view.unchecked(wrap_cell(destination), 0U);
        }
      }
    }
  };
  const auto cell_index = [&](Int3 cell) {
    return static_cast<std::size_t>(cell.x) +
           static_cast<std::size_t>(cells.x) *
               (static_cast<std::size_t>(cell.y) +
                static_cast<std::size_t>(cells.y) *
                    static_cast<std::size_t>(cell.z));
  };
  const std::size_t local_cells = static_cast<std::size_t>(cells.x) *
                                  static_cast<std::size_t>(cells.y) *
                                  static_cast<std::size_t>(cells.z);

  OwnedField density = ghosted_field(0U, 9001U, 10001U, cells, 2U);
  OwnedField velocity = ghosted_components_field(
      1U, 9002U, 10002U, cells, reach, 3U);
  OwnedField target_pressure =
      ghosted_field(2U, 9003U, 10003U, cells, 2U);
  OwnedField target_enthalpy =
      ghosted_field(3U, 9004U, 10004U, cells, 2U);
  OwnedField target_temperature =
      ghosted_field(4U, 9005U, 10005U, cells, 1U);
  OwnedField momentum_diagonal = ghosted_components_field(
      30U, 9006U, 10006U, cells, 0U, 3U);
  OwnedField momentum_rhs = ghosted_components_field(
      31U, 9007U, 10007U, cells, 0U, 3U);
  constexpr std::array<double, 3U> target_velocity{{0.25, -0.125,
                                                    0.0625}};
  constexpr std::array<double, 3U> momentum_diagonal_value{{2.0, 4.0,
                                                            8.0}};
  for (std::int32_t z = 0; z < cells.z; ++z) {
    for (std::int32_t y = 0; y < cells.y; ++y) {
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double pressure =
            37.0 * x - 19.0 * y + 13.0 * z + 3.0 * x * y;
        const double enthalpy =
            300000.0 + 170.0 * x - 90.0 * y + 55.0 * z;
        target_pressure.view.unchecked(cell, 0U) = pressure;
        target_enthalpy.view.unchecked(cell, 0U) = enthalpy;
        target_temperature.view.unchecked(cell, 0U) = enthalpy / 1000.0;
        density.view.unchecked(cell, 0U) =
            ideal_density_factor * (pressure_reference + pressure) / enthalpy;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          velocity.view.unchecked(cell, component) =
              target_velocity[component];
          momentum_diagonal.view.unchecked(cell, component) =
              momentum_diagonal_value[component];
          momentum_rhs.view.unchecked(cell, component) =
              momentum_diagonal_value[component] *
              target_velocity[component];
        }
      }
    }
  }
  fill_periodic_scalar(target_pressure.view, 2U);
  fill_periodic_scalar(target_enthalpy.view, 2U);
  fill_periodic_scalar(target_temperature.view, 1U);
  fill_periodic_scalar(density.view, 2U);
  std::fill(velocity.storage.begin(), velocity.storage.end(), 0.0);
  for (std::int32_t z = -static_cast<std::int32_t>(reach);
       z < cells.z + static_cast<std::int32_t>(reach); ++z)
    for (std::int32_t y = -static_cast<std::int32_t>(reach);
         y < cells.y + static_cast<std::int32_t>(reach); ++y)
      for (std::int32_t x = -static_cast<std::int32_t>(reach);
           x < cells.x + static_cast<std::int32_t>(reach); ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component)
          velocity.view.unchecked({x, y, z}, component) =
              target_velocity[component];

  OwnedFaces trial_flux_storage = face_bundle(cells, 10008U, 11008U);
  FaceFluxView trial_flux{trial_flux_storage.x, trial_flux_storage.y,
                          trial_flux_storage.z, 9008U, {}};
  const std::array<ConstFieldView, 2U> flux_reads{
      as_const(density.view), as_const(velocity.view)};
  const KernelInvocation flux_invocation{
      {flux_reads.data(), flux_reads.size()}, {}, {{0, 0, 0}, cells},
      0U, 0U, 1U, 0U, nullptr};
  passed &= expect(static_cast<bool>(reconstruct_mass_flux(
                       equations.kernels(), flux_invocation, trial_flux)),
                   "analytic frozen-target trial flux reconstructs");

  OwnedField r_au = ghosted_components_field(
      40U, 9009U, 10009U, cells, 1U, 3U);
  OwnedField h_by_a = ghosted_components_field(
      41U, 9010U, 10010U, cells, reach, 3U);
  OwnedField pressure_gradient = ghosted_components_field(
      42U, 9011U, 10011U, cells, 1U, 3U);
  OwnedFaces pressure_coefficient = face_bundle(cells, 10012U, 11012U);
  OwnedFaces phi_h_by_a_storage = face_bundle(cells, 10013U, 11013U);
  FaceFluxView phi_h_by_a{phi_h_by_a_storage.x, phi_h_by_a_storage.y,
                          phi_h_by_a_storage.z, 9012U, {}};
  const PisoCouplerWorkspace coupler_workspace{
      r_au.view,
      h_by_a.view,
      pressure_gradient.view,
      pressure_coefficient.x,
      pressure_coefficient.y,
      pressure_coefficient.z,
      phi_h_by_a};
  const std::array<HaloFieldSpec, 4U> coupler_halo_fields{{
      {density.view.field, 1U, 1U},
      {r_au.view.field, 1U, 3U},
      {h_by_a.view.field, reach, 3U},
      {pressure_gradient.view.field, 1U, 3U},
  }};
  HaloEngine coupler_halo;
  constexpr FieldId pressure_direction_field = 90U;
  const std::array<HaloFieldSpec, 1U> pressure_halo_fields{{
      {pressure_direction_field, 1U, 1U},
  }};
  HaloEngine pressure_halo;
  passed &= expect(
      static_cast<bool>(coupler_halo.reserve(
          MPI_COMM_SELF, fixture.patch,
          {coupler_halo_fields.data(), coupler_halo_fields.size()},
          fixture.boundary.halo_topology())) &&
          static_cast<bool>(pressure_halo.reserve(
              MPI_COMM_SELF, fixture.patch,
              {pressure_halo_fields.data(), pressure_halo_fields.size()},
              fixture.boundary.halo_topology())),
      "analytic frozen-target fixture reserves PISO and pressure halos");
  const PisoCouplerServices coupler_services{
      MPI_COMM_SELF, &fixture.geometry, fixture.patch, &fixture.boundary,
      &thermodynamics, &coupler_halo, 911U, density.view.field,
      &pressure_halo, 912U,
      pressure_direction_field};
  PressureVelocityCoupler coupler;
  passed &= expect(static_cast<bool>(PressureVelocityCoupler::bind(
                       piso, equations, coupler_services, coupler_workspace,
                       coupler)),
                   "analytic frozen-target binds the pressure-velocity coupler");
  if (!passed) return false;

  FieldRegistry history_registry;
  FieldSchema history_schema;
  FieldId history_dependency = 0U;
  StateLayers history_layers;
  AttemptTransaction history_transaction;
  FaceFluxStorage history_storage;
  FinalFaceFluxAuthority history_authority;
  FinalFaceFluxWriter history_writer;
  const std::array<ArenaFieldRequest, 1U> history_requests{{
      {0U, {1, 1, 1}, {0U}, FieldLifetime::state_layer},
  }};
  ArenaLayout history_layout;
  passed &= expect(
      static_cast<bool>(history_registry.declare_field(
          "full_block_fd_flux_history", 1U, 0U, history_dependency)) &&
          history_dependency == 0U &&
          static_cast<bool>(history_registry.freeze(history_schema)) &&
          static_cast<bool>(ArenaLayout::compile(
              history_schema,
              {history_requests.data(), history_requests.size()},
              history_layout)) &&
          static_cast<bool>(StateLayers::allocate(history_layout,
                                                  history_layers)) &&
          static_cast<bool>(AttemptTransaction::create(
              history_layers.field_count(), 1U,
              history_layers.field_count(), history_transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(
              cells, history_storage)) &&
          static_cast<bool>(history_authority.claim(
              913U, 0U, history_transaction, history_writer)),
      "analytic frozen-target allocates certified BDF2 face history");
  auto commit_history = [&](ConstFaceFluxView& committed) {
    if (!history_transaction.begin(history_layers) ||
        !history_transaction.revise_trial(history_dependency))
      return false;
    const RevisionDependency dependency{
        AttemptTransaction::field_revision_source(history_dependency),
        history_transaction.trial_revision(history_dependency)};
    PendingFaceFluxView pending;
    if (!history_writer.begin_pending(history_transaction, history_storage,
                                      pending))
      return false;
    const std::array dependencies{dependency};
    return static_cast<bool>(reconstruct_mass_flux(
               equations.kernels(), flux_invocation, pending)) &&
           static_cast<bool>(history_writer.publish_pending(
               {dependencies.data(), dependencies.size()}, pending)) &&
           static_cast<bool>(history_transaction.collective_finish(
               MPI_COMM_SELF, Status{})) &&
           static_cast<bool>(history_writer.committed(history_storage,
                                                      committed));
  };
  ConstFaceFluxView previous_flux;
  ConstFaceFluxView accepted_flux;
  passed &= expect(commit_history(previous_flux) &&
                       commit_history(accepted_flux),
                   "analytic frozen-target publishes committed flux history");
  if (!passed) return false;

  PisoIntermediateInput intermediate_input;
  intermediate_input.momentum = {
      equations.momentum().fingerprint(), EquationAssemblyScope::momentum_predictor,
      9201U, fixture.geometry.topology_revision(), trial_flux.revision, 9202U,
      step_dt};
  intermediate_input.predictor.plan =
      equations.thermophysical_predictor().fingerprint();
  intermediate_input.predictor.time = 9201U;
  intermediate_input.predictor.geometry =
      fixture.geometry.topology_revision();
  intermediate_input.predictor.accepted_face_flux = accepted_flux.revision;
  intermediate_input.predictor.previous_face_flux = previous_flux.revision;
  intermediate_input.predictor.committed_face_flux_authority =
      accepted_flux.certificate.authority();
  intermediate_input.predictor.committed_face_flux_storage =
      accepted_flux.certificate.storage();
  intermediate_input.predictor.committed_face_flux_revision_domain =
      accepted_flux.certificate.revision_domain();
  intermediate_input.predictor.predicted_density = density.view.revision;
  intermediate_input.predictor.predicted_density_storage =
      density.view.storage_identity;
  intermediate_input.predictor.predicted_density_revision_domain =
      density.view.revision_domain;
  intermediate_input.predictor.paired_face_flux = trial_flux.revision;
  intermediate_input.predictor.paired_face_flux_storage =
      trial_flux.x.storage_identity;
  intermediate_input.predictor.paired_face_flux_revision_domain =
      trial_flux.x.revision_domain;
  intermediate_input.predictor.state = 9203U;
  intermediate_input.predictor.order = 2U;
  intermediate_input.pressure_reference = {
      equations.pressure_reference().fingerprint(),
      equations.thermophysical_predictor().fingerprint(),
      thermodynamics.fingerprint(), 9204U, 9201U, 9205U,
      PressureReferenceKind::closed_mass};
  intermediate_input.density = density.view;
  intermediate_input.trial_velocity = as_const(velocity.view);
  intermediate_input.trial_flux = as_const(trial_flux);
  intermediate_input.momentum_system = {momentum_diagonal.view,
                                        momentum_rhs.view};
  intermediate_input.bdf = bdf;
  intermediate_input.numeric_boundary = fixture.boundary.revision();
  intermediate_input.corrector = 1U;
  intermediate_input.temporal_reference = as_const(phi_h_by_a);
  intermediate_input.committed_face_history = {accepted_flux, previous_flux};
  intermediate_input.thermophysical_boundary.binding.pressure_reference =
      pressure_reference;
  intermediate_input.thermophysical_boundary.binding.pressure_perturbation =
      as_const(target_pressure.view);
  PisoIntermediateCertificate intermediate;
  passed &= expect(static_cast<bool>(coupler.refresh(intermediate_input,
                                                     intermediate)) &&
                       intermediate.valid() && intermediate.corrector == 1U,
                   "analytic frozen-target refreshes the real C1 layer");
  if (!passed) return false;

  OwnedField density_accepted =
      shaped_field(50U, 9206U, 10206U, cells, 0.0);
  OwnedField density_previous =
      shaped_field(51U, 9207U, 10207U, cells, 0.0);
  OwnedField density_pressure_derivative =
      shaped_field(8U, 9208U, 10208U, cells, 0.0);
  OwnedField density_enthalpy_derivative =
      shaped_field(52U, 9209U, 10209U, cells, 0.0);
  OwnedField pressure_diagonal =
      shaped_field(53U, 9210U, 10210U, cells, 0.0);
  OwnedField pressure_rhs =
      shaped_field(54U, 9211U, 10211U, cells, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double rho = density.view.unchecked(cell, 0U);
        const double h = target_enthalpy.view.unchecked(cell, 0U);
        density_accepted.view.unchecked(cell, 0U) = 0.99 * rho;
        density_previous.view.unchecked(cell, 0U) = 1.01 * rho;
        density_pressure_derivative.view.unchecked(cell, 0U) =
            ideal_density_factor / h;
        density_enthalpy_derivative.view.unchecked(cell, 0U) = -rho / h;
      }
  const PressureCorrectionInput pressure_input{
      intermediate,
      intermediate_input.pressure_reference,
      as_const(density.view),
      as_const(density_accepted.view),
      as_const(density_previous.view),
      as_const(density_pressure_derivative.view),
      bdf,
      intermediate_input.momentum.time,
      intermediate_input.momentum.geometry,
      intermediate_input.numeric_boundary};
  const PressureCorrectionSystemView pressure_system{pressure_diagonal.view,
                                                       pressure_rhs.view};
  PressureCorrectionCertificate pressure_certificate;
  passed &= expect(static_cast<bool>(coupler.assemble_pressure_system(
                       pressure_input, pressure_system,
                       pressure_certificate)) &&
                       pressure_certificate.valid(),
                   "analytic frozen-target assembles the real C_p system");

  PressureLinearOperator continuity_pressure;
  passed &= expect(
      static_cast<bool>(coupler.bind_pressure_operator(
          {MPI_COMM_SELF, &pressure_halo, 914U, pressure_direction_field, 1U},
          pressure_system, continuity_pressure)) &&
          static_cast<bool>(continuity_pressure.refresh(
              {pressure_certificate,
               {9401U, 9402U, 9403U, 9404U, 9405U}, 9306U})) &&
          continuity_pressure.certificate().operator_class ==
              LinearOperatorClass::spd,
      "analytic frozen-target binds the actual Cartesian C_p operator");
  ConstFaceFluxView target_flux;
  passed &= expect(static_cast<bool>(coupler.inspect_intermediate_flux(
                       intermediate, target_flux)),
                   "analytic frozen-target inspects the same frozen mass flux");
  if (!passed) return false;

  OwnedFaces frozen_enthalpy_storage =
      face_bundle(cells, 10212U, 11212U);
  OwnedFaces directional_enthalpy_storage =
      face_bundle(cells, 10213U, 11213U);
  const FrozenConvectionContext convection_context{
      equations.enthalpy().fingerprint(), fixture.boundary.revision()};
  FrozenConvectionFaceField frozen_enthalpy;
  passed &= expect(static_cast<bool>(freeze_cartesian_target_convection_faces(
                       equations.kernels(), fixture.schemes.enthalpy(),
                       target_flux, as_const(target_enthalpy.view), 0U,
                       convection_context,
                       {frozen_enthalpy_storage.x,
                        frozen_enthalpy_storage.y,
                        frozen_enthalpy_storage.z},
                       frozen_enthalpy)),
                   "analytic frozen-target freezes enthalpy reconstruction");
  if (!passed) return false;

  const auto select_const_face = [](ConstFaceFieldView x,
                                    ConstFaceFieldView y,
                                    ConstFaceFieldView z,
                                    CartesianAxis axis) noexcept {
    return axis == CartesianAxis::x ? x
           : axis == CartesianAxis::y ? y
                                      : z;
  };
  const auto face_cells = [&](CartesianAxis axis, Int3 face,
                              Int3& left, Int3& right) {
    right = face;
    left = face;
    if (axis == CartesianAxis::x) {
      left.x = face.x - 1;
      right.x = face.x;
    } else if (axis == CartesianAxis::y) {
      left.y = face.y - 1;
      right.y = face.y;
    } else {
      left.z = face.z - 1;
      right.z = face.z;
    }
    left = wrap_cell(left);
    right = wrap_cell(right);
  };
  double largest_frozen_face_error = 0.0;
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const ConstFaceFieldView frozen = select_const_face(
        frozen_enthalpy.x, frozen_enthalpy.y, frozen_enthalpy.z, axis);
    for (std::int32_t z = 0; z < frozen.extents.z; ++z)
      for (std::int32_t y = 0; y < frozen.extents.y; ++y)
        for (std::int32_t x = 0; x < frozen.extents.x; ++x) {
          Int3 left;
          Int3 right;
          face_cells(axis, {x, y, z}, left, right);
          const double independent_central =
              0.5 * (target_enthalpy.view.unchecked(left, 0U) +
                     target_enthalpy.view.unchecked(right, 0U));
          largest_frozen_face_error =
              std::max(largest_frozen_face_error,
                       std::abs(frozen.unchecked({x, y, z}) -
                                independent_central));
        }
  }
  passed &= expect(largest_frozen_face_error < 1.0e-12,
                   "analytic frozen-target oracle verifies central target "
                   "face enthalpy");

  OwnedField continuity_enthalpy =
      shaped_field(55U, 9212U, 10214U, cells, 0.0);
  OwnedField continuity_enthalpy_scale =
      shaped_field(56U, 9213U, 10215U, cells, 0.0);
  OwnedField energy_pressure_temporal =
      shaped_field(57U, 9214U, 10216U, cells, 0.0);
  OwnedField assembled_enthalpy =
      shaped_field(58U, 9215U, 10217U, cells, 0.0);
  OwnedField heat_capacity =
      ghosted_field(59U, 9216U, 10218U, cells, 1U, 1000.0);
  OwnedField thermal_conductivity =
      ghosted_field(60U, 9217U, 10219U, cells, 1U, 2.0);
  OwnedField enthalpy_diffusivity =
      ghosted_field(61U, 9218U, 10220U, cells, 1U, 0.002);
  OwnedField delta_temperature =
      ghosted_field(92U, 9219U, 10221U, cells, 1U, 0.0);
  const double diffusion_proxy = 6.0 * 0.002;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double volume = 1.0;
        const double rho = density.view.unchecked(cell, 0U);
        const double h = target_enthalpy.view.unchecked(cell, 0U);
        const double rho_p =
            density_pressure_derivative.view.unchecked(cell, 0U);
        const double rho_h =
            density_enthalpy_derivative.view.unchecked(cell, 0U);
        continuity_enthalpy.view.unchecked(cell, 0U) =
            bdf.a0 * volume * rho_h;
        continuity_enthalpy_scale.view.unchecked(cell, 0U) =
            std::abs(bdf.a0 * volume * rho_h);
        energy_pressure_temporal.view.unchecked(cell, 0U) =
            bdf.a0 * volume * (h * rho_p - 1.0);
        assembled_enthalpy.view.unchecked(cell, 0U) =
            bdf.a0 * volume * rho + diffusion_proxy;
      }

  PressureEnergyFrozenFaceEnthalpy pressure_frozen_enthalpy{
      frozen_enthalpy.x, frozen_enthalpy.y, frozen_enthalpy.z,
      frozen_enthalpy.revision, frozen_enthalpy.reconstruction, 0U};
  pressure_frozen_enthalpy.local_binding =
      pressure_energy_frozen_face_enthalpy_local_binding(
          pressure_frozen_enthalpy);
  PisoCartesianPressureWorkLinearization pressure_work;
  passed &= expect(
      static_cast<bool>(coupler.inspect_cartesian_pressure_work_linearization(
          intermediate, pressure_certificate, as_const(target_pressure.view),
          as_const(velocity.view), pressure_work)) &&
          pressure_work.valid(),
      "analytic frozen-target obtains live Cartesian pressure-work "
      "authority");
  const LinearIdentity block_identity{9401U, 9402U, 9403U, 9404U, 9405U};
  PressureEnergyPressureFluxBinding energy_pressure_binding;
  energy_pressure_binding.geometry = &fixture.geometry;
  energy_pressure_binding.boundary = &fixture.boundary;
  energy_pressure_binding.patch = fixture.patch;
  energy_pressure_binding.services = {
      MPI_COMM_SELF, &pressure_halo, 914U, pressure_direction_field, 1U};
  energy_pressure_binding.intermediate = intermediate;
  energy_pressure_binding.pressure = pressure_certificate;
  energy_pressure_binding.temporal_diagonal =
      as_const(energy_pressure_temporal.view);
  energy_pressure_binding.x_pressure_coefficient =
      as_const(pressure_coefficient.x);
  energy_pressure_binding.y_pressure_coefficient =
      as_const(pressure_coefficient.y);
  energy_pressure_binding.z_pressure_coefficient =
      as_const(pressure_coefficient.z);
  energy_pressure_binding.target_flux = target_flux;
  energy_pressure_binding.frozen_face_enthalpy = pressure_frozen_enthalpy;
  energy_pressure_binding.identity = block_identity;
  energy_pressure_binding.pressure_work = pressure_work;
  PressureEnergyPressureFluxOperator energy_pressure;
  PressureEnergyPressureFluxCertificate energy_pressure_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          energy_pressure_binding, energy_pressure,
          energy_pressure_certificate)) &&
          energy_pressure_certificate.full_cartesian_pressure_work &&
          !energy_pressure_certificate.flux_only_quasi_newton,
      "analytic frozen-target binds the complete Cartesian E_p operator");

  constexpr FieldId enthalpy_direction_field = 91U;
  const std::array<HaloFieldSpec, 2U> energy_halo_fields{{
      {enthalpy_direction_field, 2U, 1U},
      {delta_temperature.view.field, 1U, 1U},
  }};
  HaloEngine energy_halo;
  passed &= expect(static_cast<bool>(energy_halo.reserve(
                       MPI_COMM_SELF, fixture.patch,
                       {energy_halo_fields.data(), energy_halo_fields.size()},
                       fixture.boundary.halo_topology())),
                   "analytic frozen-target reserves the exact E_h halo");
  PressureEnergyEnthalpyBinding energy_enthalpy_binding;
  energy_enthalpy_binding.geometry = &fixture.geometry;
  energy_enthalpy_binding.kernels = &equations.kernels();
  energy_enthalpy_binding.boundary = &fixture.boundary;
  energy_enthalpy_binding.patch = fixture.patch;
  energy_enthalpy_binding.convection = fixture.schemes.enthalpy();
  energy_enthalpy_binding.services = {
      MPI_COMM_SELF, &energy_halo, 916U, enthalpy_direction_field,
      delta_temperature.view.field};
  energy_enthalpy_binding.authority = {
      bdf,
      intermediate_input.momentum.time,
      fixture.geometry.topology_revision(),
      fixture.boundary.revision(),
      9406U,
      9407U,
      convection_context.collective_semantics,
      thermodynamics.fingerprint(),
      transport.fingerprint()};
  energy_enthalpy_binding.assembled_diagonal =
      as_const(assembled_enthalpy.view);
  energy_enthalpy_binding.target_enthalpy = as_const(target_enthalpy.view);
  energy_enthalpy_binding.density_enthalpy_derivative =
      as_const(density_enthalpy_derivative.view);
  energy_enthalpy_binding.heat_capacity = as_const(heat_capacity.view);
  energy_enthalpy_binding.thermal_conductivity =
      as_const(thermal_conductivity.view);
  energy_enthalpy_binding.enthalpy_diffusivity =
      as_const(enthalpy_diffusivity.view);
  energy_enthalpy_binding.target_flux = target_flux;
  energy_enthalpy_binding.convection_context = convection_context;
  energy_enthalpy_binding.frozen_face_enthalpy = frozen_enthalpy;
  energy_enthalpy_binding.workspace = {
      delta_temperature.view,
      {directional_enthalpy_storage.x, directional_enthalpy_storage.y,
       directional_enthalpy_storage.z}};
  energy_enthalpy_binding.identity = block_identity;
  energy_enthalpy_binding.linearization_policy =
      FrozenConvectionLinearizationPolicy::semismooth_generalized_zero_slope;
  PressureEnergyEnthalpyOperator energy_enthalpy;
  PressureEnergyEnthalpyCertificate energy_enthalpy_certificate;
  passed &= expect(
      static_cast<bool>(PressureEnergyEnthalpyOperator::bind(
          energy_enthalpy_binding, energy_enthalpy,
          energy_enthalpy_certificate)) &&
          energy_enthalpy_certificate.exact_cartesian_spatial_response &&
          energy_enthalpy_certificate.exact_temperature_space_conduction,
      "analytic frozen-target binds the frozen-spatial E_h operator");

  PressureEnergyDiagonalOperator continuity_enthalpy_operator;
  PressureEnergyDiagonalCertificate continuity_enthalpy_certificate;
  passed &= expect(static_cast<bool>(PressureEnergyDiagonalOperator::bind(
                       {as_const(continuity_enthalpy.view), {},
                        {9411U, 9412U, 9413U, 9414U, 9415U}, 0.0},
                       continuity_enthalpy_operator,
                       continuity_enthalpy_certificate)),
                   "analytic frozen-target binds the production C_h diagonal");

  PressureEnergySchurBlockAuthority block_authority;
  OwnedField schur_continuity_workspace =
      shaped_field(62U, 9220U, 10222U, cells, 0.0);
  OwnedField schur_enthalpy_workspace = ghosted_field(
      enthalpy_direction_field, 9221U, 10223U, cells, 2U, 0.0);
  OwnedField schur_energy_workspace =
      shaped_field(63U, 9222U, 10224U, cells, 0.0);
  PressureEnergySchurBinding schur_binding;
  schur_binding.continuity_pressure = &continuity_pressure;
  schur_binding.energy_pressure = &energy_pressure;
  schur_binding.energy_enthalpy = &energy_enthalpy;
  schur_binding.continuity_enthalpy_diagonal =
      as_const(continuity_enthalpy.view);
  schur_binding.continuity_enthalpy_row_scale =
      as_const(continuity_enthalpy_scale.view);
  schur_binding.workspace = {schur_continuity_workspace.view,
                             schur_enthalpy_workspace.view,
                             schur_energy_workspace.view};
  schur_binding.scaled_pivot_floor = 1.0e-12;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurBlockAuthority::exact_cartesian(
          energy_pressure, energy_enthalpy, block_authority)) &&
          block_authority.scope() == PressureEnergySchurBlockScope::
                                           exact_cartesian_frozen_spatial,
      "analytic frozen-target derives typed Cartesian block authority");
  schur_binding.block_authority = block_authority;
  PressureEnergySchurOperator schur;
  PressureEnergyJacobianCertificate jacobian;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurOperator::bind(
          schur_binding, schur, jacobian)) &&
          jacobian.valid() &&
          jacobian.jacobian_scope == PressureEnergyJacobianScope::
                                           exact_cartesian_frozen_spatial &&
          jacobian.exact_algebraic_schur &&
          !jacobian.full_nonlinear_jacobian &&
          jacobian.shared_pressure.valid() &&
          jacobian.shared_pressure.scope() ==
              PressureEnergySharedPressureScope::cartesian,
      "analytic frozen-target binds exact Schur and typed shared-pressure halo");
  if (!passed) return false;

  auto mismatched_pressure_binding = energy_pressure_binding;
  mismatched_pressure_binding.services.halo_stage = 915U;
  PressureEnergyPressureFluxOperator mismatched_energy_pressure;
  PressureEnergyPressureFluxCertificate mismatched_energy_certificate;
  PressureEnergySchurBlockAuthority mismatched_block_authority;
  PressureEnergySchurBinding mismatched_schur_binding = schur_binding;
  PressureEnergySchurOperator mismatched_schur;
  PressureEnergyJacobianCertificate mismatched_jacobian;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          mismatched_pressure_binding, mismatched_energy_pressure,
          mismatched_energy_certificate)) &&
          static_cast<bool>(PressureEnergySchurBlockAuthority::exact_cartesian(
              mismatched_energy_pressure, energy_enthalpy,
              mismatched_block_authority)),
      "stage-mutation fixture keeps the same mathematical E_p/E_h block");
  mismatched_schur_binding.energy_pressure = &mismatched_energy_pressure;
  mismatched_schur_binding.block_authority = mismatched_block_authority;
  passed &= expect(
      static_cast<bool>(PressureEnergySchurOperator::bind(
          mismatched_schur_binding, mismatched_schur,
          mismatched_jacobian)) &&
          mismatched_jacobian.valid() &&
          !mismatched_jacobian.shared_pressure.available(),
      "a mismatched halo stage keeps generic Schur composition and cannot "
      "reuse a ghost revision");
  if (!passed) return false;

  OwnedField pressure_direction = ghosted_field(
      pressure_direction_field, 9223U, 10225U, cells, 1U, 0.0);
  OwnedField enthalpy_direction = ghosted_field(
      enthalpy_direction_field, 9224U, 10226U, cells, 2U, 0.0);
  OwnedField continuity_pressure_action =
      shaped_field(64U, 9225U, 10227U, cells, 0.0);
  OwnedField continuity_enthalpy_action =
      shaped_field(65U, 9226U, 10228U, cells, 0.0);
  OwnedField energy_pressure_action =
      shaped_field(66U, 9227U, 10229U, cells, 0.0);
  OwnedField energy_enthalpy_action =
      shaped_field(67U, 9228U, 10230U, cells, 0.0);
  std::vector<double> base_pressure(local_cells, 0.0);
  std::vector<double> base_enthalpy(local_cells, 0.0);
  std::vector<double> pressure_probe(local_cells, 0.0);
  std::vector<double> enthalpy_probe(local_cells, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::size_t index = cell_index(cell);
        base_pressure[index] = target_pressure.view.unchecked(cell, 0U);
        base_enthalpy[index] = target_enthalpy.view.unchecked(cell, 0U);
        pressure_probe[index] =
            4.0 + 0.7 * x - 0.45 * y + 0.31 * z + 0.08 * x * z;
        enthalpy_probe[index] =
            90.0 - 8.0 * x + 5.5 * y - 3.25 * z + 0.6 * y * z;
      }

  const auto face_view = [&](const OwnedFaces& faces,
                             CartesianAxis axis) noexcept {
    return axis == CartesianAxis::x ? as_const(faces.x)
           : axis == CartesianAxis::y ? as_const(faces.y)
                                      : as_const(faces.z);
  };
  const auto target_flux_face = [&](CartesianAxis axis) noexcept {
    return axis == CartesianAxis::x ? target_flux.x
           : axis == CartesianAxis::y ? target_flux.y
                                      : target_flux.z;
  };
  const auto shifted = [&](Int3 cell, CartesianAxis axis,
                           int amount) noexcept {
    if (axis == CartesianAxis::x)
      cell.x += amount;
    else if (axis == CartesianAxis::y)
      cell.y += amount;
    else
      cell.z += amount;
    return wrap_cell(cell);
  };

  struct ResidualPair {
    std::vector<double> continuity;
    std::vector<double> energy;
  };
  const auto residual = [&](double pressure_factor,
                            double enthalpy_factor) {
    ResidualPair result{{}, {}};
    result.continuity.assign(local_cells, 0.0);
    result.energy.assign(local_cells, 0.0);
    const auto p = [&](Int3 cell) {
      const std::size_t index = cell_index(wrap_cell(cell));
      return base_pressure[index] +
             pressure_factor * pressure_probe[index];
    };
    const auto h = [&](Int3 cell) {
      const std::size_t index = cell_index(wrap_cell(cell));
      return base_enthalpy[index] +
             enthalpy_factor * enthalpy_probe[index];
    };
    const auto delta_p = [&](Int3 cell) {
      const std::size_t index = cell_index(wrap_cell(cell));
      return pressure_factor * pressure_probe[index];
    };
    const auto face_flux_delta = [&](CartesianAxis axis, Int3 face) {
      Int3 left;
      Int3 right;
      face_cells(axis, face, left, right);
      const double coefficient =
          face_view(pressure_coefficient, axis).unchecked(face);
      return -coefficient * (delta_p(right) - delta_p(left));
    };
    const auto target_face_h = [&](CartesianAxis axis, Int3 face) {
      Int3 left;
      Int3 right;
      face_cells(axis, face, left, right);
      return 0.5 * (base_enthalpy[cell_index(left)] +
                    base_enthalpy[cell_index(right)]);
    };
    const auto face_h_delta = [&](CartesianAxis axis, Int3 face) {
      Int3 left;
      Int3 right;
      face_cells(axis, face, left, right);
      return 0.5 * enthalpy_factor *
             (enthalpy_probe[cell_index(left)] +
              enthalpy_probe[cell_index(right)]);
    };
    const auto energy_face_flux_delta = [&](CartesianAxis axis, Int3 face) {
      const double delta_flux = face_flux_delta(axis, face);
      const double delta_enthalpy = face_h_delta(axis, face);
      return target_flux_face(axis).unchecked(face) * delta_enthalpy +
             delta_flux * target_face_h(axis, face) +
             delta_flux * delta_enthalpy;
    };
    const auto gradient = [&](auto&& value, Int3 cell,
                              CartesianAxis axis) {
      return 0.5 * (value(shifted(cell, axis, 1)) -
                    value(shifted(cell, axis, -1)));
    };
    for (std::int32_t z = 0; z < cells.z; ++z) {
      for (std::int32_t y = 0; y < cells.y; ++y) {
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const std::size_t index = cell_index(cell);
          const double pressure_absolute = pressure_reference + p(cell);
          const double rho =
              ideal_density_factor * pressure_absolute / h(cell);
          const double target_rho = density.view.unchecked(cell, 0U);
          // For this calorically-perfect ideal-gas oracle, rho*h-p reduces
          // analytically to (cp/R-1)*p.  Evaluate that closed form so the FD
          // measures the block, not cancellation of two O(1e5) products.
          const double q =
              (ideal_density_factor - 1.0) * pressure_absolute;
          const double target_pressure_absolute =
              pressure_reference + base_pressure[index];
          const double target_q = (ideal_density_factor - 1.0) *
                                  target_pressure_absolute;
          double continuity_flux_change = 0.0;
          double energy_flux_change = 0.0;
          double conduction_change = 0.0;
          double pressure_work = 0.0;
          double target_pressure_work = 0.0;
          for (std::size_t component = 0U; component < 3U; ++component) {
            const CartesianAxis axis = static_cast<CartesianAxis>(component);
            Int3 plus_face = cell;
            if (axis == CartesianAxis::x)
              ++plus_face.x;
            else if (axis == CartesianAxis::y)
              ++plus_face.y;
            else
              ++plus_face.z;
            continuity_flux_change += face_flux_delta(axis, plus_face) -
                                      face_flux_delta(axis, cell);
            energy_flux_change +=
                energy_face_flux_delta(axis, plus_face) -
                energy_face_flux_delta(axis, cell);
            for (int direction : {-1, 1}) {
              const Int3 neighbour = shifted(cell, axis, direction);
              const double temperature = h(cell) / 1000.0;
              const double neighbour_temperature = h(neighbour) / 1000.0;
              const double target_temperature_value =
                  base_enthalpy[index] / 1000.0;
              const double target_neighbour_temperature =
                  base_enthalpy[cell_index(neighbour)] / 1000.0;
              conduction_change +=
                  2.0 * ((temperature - neighbour_temperature) -
                         (target_temperature_value -
                          target_neighbour_temperature));
            }
            const double target_gradient = gradient(
                [&](Int3 selected) {
                  return base_pressure[cell_index(selected)];
                },
                cell, axis);
            const double pressure_gradient_value =
                gradient(p, cell, axis);
            const double correction_gradient =
                gradient(delta_p, cell, axis);
            const double reciprocal =
                r_au.view.unchecked(cell,
                                    static_cast<std::uint8_t>(component));
            const double predictor_velocity =
                h_by_a.view.unchecked(cell,
                                      static_cast<std::uint8_t>(component));
            pressure_work +=
                -(predictor_velocity - reciprocal * correction_gradient) *
                pressure_gradient_value;
            target_pressure_work +=
                -predictor_velocity * target_gradient;
          }
          result.continuity[index] =
              bdf.a0 * (rho - target_rho) + continuity_flux_change;
          result.energy[index] =
              bdf.a0 * (q - target_q) + energy_flux_change +
              conduction_change + pressure_work - target_pressure_work;
        }
      }
    }
    return result;
  };

  struct DirectionCase {
    double pressure{};
    double enthalpy{};
    std::string_view name{};
  };
  constexpr std::array<DirectionCase, 3U> directions{{
      {1.0, 0.0, "p-only"},
      {0.0, 1.0, "h-only"},
      {1.0, 1.0, "mixed"},
  }};
  constexpr std::array<double, 3U> epsilons{{1.0e-2, 1.0e-3, 1.0e-4}};
  std::array<double, 3U> continuity_errors{};
  std::array<double, 3U> energy_errors{};
  std::array<double, 3U> continuity_platform_errors{};
  std::array<double, 3U> energy_platform_errors{};
  std::array<double, 3U> continuity_signal{};
  std::array<double, 3U> energy_signal{};
  for (std::size_t direction_index = 0U;
       direction_index < directions.size(); ++direction_index) {
    const DirectionCase selected = directions[direction_index];
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const std::size_t index = cell_index(cell);
          pressure_direction.view.unchecked(cell, 0U) =
              selected.pressure * pressure_probe[index];
          enthalpy_direction.view.unchecked(cell, 0U) =
              selected.enthalpy * enthalpy_probe[index];
        }
    passed &= expect(
        static_cast<bool>(continuity_pressure.apply(
            pressure_direction.view, continuity_pressure_action.view)) &&
            static_cast<bool>(continuity_enthalpy_operator.apply(
                enthalpy_direction.view,
                continuity_enthalpy_action.view)) &&
            static_cast<bool>(energy_pressure.apply(
                pressure_direction.view, energy_pressure_action.view)) &&
            static_cast<bool>(energy_enthalpy.apply(
                enthalpy_direction.view, energy_enthalpy_action.view)),
        "analytic frozen-target FD applies all four real Jacobian blocks");
    if (!passed) return false;
    std::vector<double> previous_continuity_fd(local_cells, 0.0);
    std::vector<double> previous_energy_fd(local_cells, 0.0);
    bool has_previous_epsilon = false;
    for (double epsilon : epsilons) {
      const ResidualPair plus = residual(epsilon * selected.pressure,
                                         epsilon * selected.enthalpy);
      const ResidualPair minus = residual(-epsilon * selected.pressure,
                                          -epsilon * selected.enthalpy);
      double level_continuity_error = 0.0;
      double level_energy_error = 0.0;
      double worst_energy_action = 0.0;
      double worst_energy_fd = 0.0;
      Int3 worst_energy_cell{};
      for (std::int32_t z = 0; z < cells.z; ++z)
        for (std::int32_t y = 0; y < cells.y; ++y)
          for (std::int32_t x = 0; x < cells.x; ++x) {
            const Int3 cell{x, y, z};
            const std::size_t index = cell_index(cell);
            const double continuity_fd =
                (plus.continuity[index] - minus.continuity[index]) /
                (2.0 * epsilon);
            const double energy_fd =
                (plus.energy[index] - minus.energy[index]) /
                (2.0 * epsilon);
            const double continuity_action =
                continuity_pressure_action.view.unchecked(cell, 0U) +
                continuity_enthalpy_action.view.unchecked(cell, 0U);
            const double energy_action =
                energy_pressure_action.view.unchecked(cell, 0U) +
                energy_enthalpy_action.view.unchecked(cell, 0U);
            continuity_signal[direction_index] =
                std::max(continuity_signal[direction_index],
                         std::abs(continuity_fd));
            energy_signal[direction_index] =
                std::max(energy_signal[direction_index], std::abs(energy_fd));
            const double continuity_error =
                std::abs(continuity_action - continuity_fd) /
                std::max({1.0, std::abs(continuity_action),
                          std::abs(continuity_fd)});
            const double energy_error =
                std::abs(energy_action - energy_fd) /
                std::max({1.0, std::abs(energy_action), std::abs(energy_fd)});
            level_continuity_error =
                std::max(level_continuity_error, continuity_error);
            if (energy_error > level_energy_error) {
              level_energy_error = energy_error;
              worst_energy_action = energy_action;
              worst_energy_fd = energy_fd;
              worst_energy_cell = cell;
            }
            if (has_previous_epsilon) {
              continuity_platform_errors[direction_index] = std::max(
                  continuity_platform_errors[direction_index],
                  std::abs(continuity_fd - previous_continuity_fd[index]) /
                      std::max({1.0, std::abs(continuity_fd),
                                std::abs(previous_continuity_fd[index])}));
              energy_platform_errors[direction_index] = std::max(
                  energy_platform_errors[direction_index],
                  std::abs(energy_fd - previous_energy_fd[index]) /
                      std::max({1.0, std::abs(energy_fd),
                                std::abs(previous_energy_fd[index])}));
            }
            previous_continuity_fd[index] = continuity_fd;
            previous_energy_fd[index] = energy_fd;
          }
      continuity_errors[direction_index] =
          std::max(continuity_errors[direction_index],
                   level_continuity_error);
      energy_errors[direction_index] =
          std::max(energy_errors[direction_index], level_energy_error);
      std::cout << "analytic-frozen-target-four-block-fd-direction "
                << selected.name << " epsilon=" << epsilon
                << " continuity=" << level_continuity_error
                << " energy=" << level_energy_error
                << " energy-action=" << worst_energy_action
                << " energy-fd=" << worst_energy_fd << " cell=("
                << worst_energy_cell.x << ',' << worst_energy_cell.y << ','
                << worst_energy_cell.z << ")\n";
      passed &= expect(
          level_continuity_error < 2.0e-7 && level_energy_error < 2.0e-7,
          "analytic frozen-target epsilon level matches all four blocks");
      has_previous_epsilon = true;
    }
    passed &= expect(
        continuity_platform_errors[direction_index] < 5.0e-8 &&
            energy_platform_errors[direction_index] < 5.0e-8,
        "analytic frozen-target FD has an epsilon-stable derivative platform");
  }
  passed &= expect(
      continuity_signal[0U] > 1.0e-5 && energy_signal[0U] > 1.0e-3 &&
          continuity_signal[1U] > 1.0e-7 && energy_signal[1U] > 1.0e-3 &&
          continuity_signal[2U] > 1.0e-5 && energy_signal[2U] > 1.0e-3,
      "analytic frozen-target directions exercise C_p/C_h and E_p/E_h "
      "without a zero-signal pass");

  OwnedField schur_action =
      shaped_field(68U, 9229U, 10231U, cells, 0.0);
  OwnedField recovery_continuity_residual =
      shaped_field(69U, 9230U, 10232U, cells, 0.0);
  OwnedField recovered_enthalpy = ghosted_field(
      enthalpy_direction_field, 9231U, 10233U, cells, 2U, 0.0);
  OwnedField recovered_energy_action =
      shaped_field(70U, 9232U, 10234U, cells, 0.0);
  OwnedField residual_elimination_energy_action =
      shaped_field(71U, 9233U, 10235U, cells, 0.0);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::size_t index = cell_index(cell);
        pressure_direction.view.unchecked(cell, 0U) =
            0.125 * pressure_probe[index];
      }
  passed &= expect(
      static_cast<bool>(continuity_pressure.apply(
          pressure_direction.view, continuity_pressure_action.view)) &&
          static_cast<bool>(energy_pressure.apply(
              pressure_direction.view, energy_pressure_action.view)),
      "analytic frozen-target applies the independent real C_p and E_p");
  const HaloRuntimeCounters pressure_halo_before_schur =
      pressure_halo.runtime_counters();
  const HaloRuntimeCounters energy_halo_before_schur =
      energy_halo.runtime_counters();
  passed &= expect(
      static_cast<bool>(
          schur.apply(pressure_direction.view, schur_action.view)),
      "analytic frozen-target applies the exact Schur through its shared halo");
  const HaloRuntimeCounters pressure_halo_after_schur =
      pressure_halo.runtime_counters();
  const HaloRuntimeCounters energy_halo_after_schur =
      energy_halo.runtime_counters();
  passed &= expect(
      pressure_halo_after_schur.begin_calls -
                  pressure_halo_before_schur.begin_calls ==
              1U &&
          pressure_halo_after_schur.finish_calls -
                  pressure_halo_before_schur.finish_calls ==
              1U &&
          energy_halo_after_schur.begin_calls -
                  energy_halo_before_schur.begin_calls ==
              1U &&
          energy_halo_after_schur.finish_calls -
                  energy_halo_before_schur.finish_calls ==
              1U,
      "one Schur apply uses one shared C_p/E_p halo plus one E_h halo");
  if (!passed) return false;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        enthalpy_direction.view.unchecked(cell, 0U) =
            continuity_pressure_action.view.unchecked(cell, 0U) /
            continuity_enthalpy.view.unchecked(cell, 0U);
      }
  passed &= expect(
      static_cast<bool>(energy_enthalpy.apply(
          enthalpy_direction.view, energy_enthalpy_action.view)),
      "analytic frozen-target applies E_h to the explicit C_h elimination");
  if (!passed) return false;
  double schur_identity_error = 0.0;
  double schur_identity_signal = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double expected =
            energy_pressure_action.view.unchecked(cell, 0U) -
            energy_enthalpy_action.view.unchecked(cell, 0U);
        const double actual = schur_action.view.unchecked(cell, 0U);
        schur_identity_signal =
            std::max(schur_identity_signal, std::abs(expected));
        schur_identity_error = std::max(
            schur_identity_error,
            std::abs(actual - expected) /
                std::max({1.0, std::abs(actual), std::abs(expected)}));
      }
  passed &= expect(
      schur_identity_signal > 1.0e-3 && schur_identity_error < 5.0e-12,
      "real exact Schur satisfies Sdp=E_pdp-E_h(C_pdp/C_h)");

  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        recovery_continuity_residual.view.unchecked(cell, 0U) =
            0.02 + 0.003 * x - 0.002 * y + 0.001 * z;
      }
  passed &= expect(
      static_cast<bool>(schur.recover_enthalpy(
          as_const(recovery_continuity_residual.view),
          pressure_direction.view, recovered_enthalpy.view)) &&
          static_cast<bool>(continuity_pressure.apply(
              pressure_direction.view, continuity_pressure_action.view)) &&
          static_cast<bool>(energy_pressure.apply(
              pressure_direction.view, energy_pressure_action.view)) &&
          static_cast<bool>(energy_enthalpy.apply(
              recovered_enthalpy.view, recovered_energy_action.view)),
      "real exact Schur recovers h and reapplies both block rows");
  if (!passed) return false;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        enthalpy_direction.view.unchecked(cell, 0U) =
            recovery_continuity_residual.view.unchecked(cell, 0U) /
            continuity_enthalpy.view.unchecked(cell, 0U);
      }
  passed &= expect(
      static_cast<bool>(energy_enthalpy.apply(
          enthalpy_direction.view,
          residual_elimination_energy_action.view)),
      "analytic frozen-target applies E_h to the residual elimination");
  if (!passed) return false;
  double recovered_continuity_row_error = 0.0;
  double recovered_energy_row_error = 0.0;
  double recovered_energy_row_signal = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double continuity_row =
            recovery_continuity_residual.view.unchecked(cell, 0U) +
            continuity_pressure_action.view.unchecked(cell, 0U) +
            continuity_enthalpy.view.unchecked(cell, 0U) *
                recovered_enthalpy.view.unchecked(cell, 0U);
        const double continuity_scale = std::max(
            {1.0,
             std::abs(recovery_continuity_residual.view.unchecked(cell, 0U)),
             std::abs(continuity_pressure_action.view.unchecked(cell, 0U)),
             std::abs(continuity_enthalpy.view.unchecked(cell, 0U) *
                      recovered_enthalpy.view.unchecked(cell, 0U))});
        recovered_continuity_row_error =
            std::max(recovered_continuity_row_error,
                     std::abs(continuity_row) / continuity_scale);
        const double recovered_energy_row =
            energy_pressure_action.view.unchecked(cell, 0U) +
            recovered_energy_action.view.unchecked(cell, 0U);
        const double expected_energy_row =
            schur_action.view.unchecked(cell, 0U) -
            residual_elimination_energy_action.view.unchecked(cell, 0U);
        recovered_energy_row_signal =
            std::max(recovered_energy_row_signal,
                     std::abs(expected_energy_row));
        recovered_energy_row_error = std::max(
            recovered_energy_row_error,
            std::abs(recovered_energy_row - expected_energy_row) /
                std::max({1.0, std::abs(recovered_energy_row),
                          std::abs(expected_energy_row)}));
      }
  passed &= expect(
      recovered_continuity_row_error < 5.0e-12 &&
          recovered_energy_row_signal > 1.0e-3 &&
          recovered_energy_row_error < 5.0e-12,
      "real recover_enthalpy closes continuity and the corresponding energy row");

  std::cout << "analytic-frozen-target-four-block-fd-not-full-refresh"
            << " cp=" << continuity_errors[0U]
            << " ep=" << energy_errors[0U]
            << " ch=" << continuity_errors[1U]
            << " eh=" << energy_errors[1U]
            << " mixed-c=" << continuity_errors[2U]
            << " mixed-e=" << energy_errors[2U]
            << " cp-platform=" << continuity_platform_errors[0U]
            << " ep-platform=" << energy_platform_errors[0U]
            << " ch-platform=" << continuity_platform_errors[1U]
            << " eh-platform=" << energy_platform_errors[1U]
            << " mixed-c-platform=" << continuity_platform_errors[2U]
            << " mixed-e-platform=" << energy_platform_errors[2U]
            << " schur-identity=" << schur_identity_error
            << " recovery-continuity=" << recovered_continuity_row_error
            << " recovery-energy=" << recovered_energy_row_error << '\n';
  return passed;
}

bool test_periodic_and_ibm_pressure_flux_semantics() {
  PressureFluxFixture periodic_fixture;
  bool passed = expect(make_pressure_flux_fixture(
                           PressureFluxBoundary::periodic, periodic_fixture),
                       "periodic E_p fixture compiles");
  if (!passed) return false;
  const Int3 cells = periodic_fixture.patch.cells;
  OwnedField temporal = shaped_field(86U, 8501U, 8502U, cells, 0.0);
  OwnedFaces coefficient = face_bundle(cells, 8503U, 8504U, 0.0);
  for (std::int32_t x = 0; x <= cells.x; ++x)
    coefficient.x.unchecked({x, 0, 0}) = 1.0;
  OwnedFaces target = face_bundle(cells, 8505U, 8506U, 0.0);
  OwnedFaces enthalpy = face_bundle(cells, 8507U, 8508U, 10.0);
  OwnedField dp = ghosted_field(2U, 8509U, 8510U, cells, 1U);
  dp.view.unchecked({0, 0, 0}, 0U) = 1.0;
  dp.view.unchecked({1, 0, 0}, 0U) = 2.0;
  dp.view.unchecked({2, 0, 0}, 0U) = 4.0;
  OwnedField output = shaped_field(87U, 8511U, 8512U, cells, 0.0);
  PressureEnergyPressureFluxOperator periodic_op;
  PressureEnergyPressureFluxCertificate periodic_certificate;
  const auto periodic_binding = pressure_flux_binding(
      periodic_fixture, as_const(temporal.view), coefficient, target, enthalpy);
  passed &=
      expect(static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
                 periodic_binding, periodic_op, periodic_certificate)) &&
                 static_cast<bool>(periodic_op.apply(dp.view, output.view)),
             "periodic pressure-flux E_p binds and exchanges its halo");
  const Vector periodic = values(as_const(output.view));
  passed &= expect(close(periodic[0U], -40.0) && close(periodic[1U], -10.0) &&
                       close(periodic[2U], 50.0) &&
                       close(periodic[0U] + periodic[1U] + periodic[2U], 0.0),
                   "periodic face response is pairwise conservative");
  passed &=
      expect(close(pressure_correction_mass_flux_response(
                       as_const(dp.view), periodic_fixture.geometry,
                       periodic_fixture.patch, periodic_fixture.boundary,
                       CartesianAxis::x, {0, 0, 0}, 1.0),
                   3.0) &&
                 close(pressure_correction_mass_flux_response(
                           as_const(dp.view), periodic_fixture.geometry,
                           periodic_fixture.patch, periodic_fixture.boundary,
                           CartesianAxis::x, {3, 0, 0}, 1.0),
                       3.0),
             "periodic domain-end faces share the same jump authority");

  PressureFluxFixture ibm_fixture;
  passed &= expect(make_pressure_flux_fixture(
                       PressureFluxBoundary::neumann_dirichlet, ibm_fixture),
                   "IBM-masked E_p fixture compiles");
  if (!passed) return false;
  OwnedField ibm_temporal = shaped_field(88U, 8513U, 8514U, cells, 0.0);
  ibm_temporal.view.unchecked({2, 0, 0}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  OwnedFaces ibm_coefficient = face_bundle(cells, 8515U, 8516U, 0.0);
  ibm_coefficient.x.unchecked({1, 0, 0}) = 1.0;
  ibm_coefficient.x.unchecked({2, 0, 0}) =
      std::numeric_limits<double>::quiet_NaN();
  ibm_coefficient.x.unchecked({3, 0, 0}) =
      std::numeric_limits<double>::quiet_NaN();
  OwnedFaces ibm_target = face_bundle(cells, 8517U, 8518U, 0.0);
  OwnedFaces ibm_enthalpy = face_bundle(cells, 8519U, 8520U, 10.0);
  ibm_enthalpy.x.unchecked({2, 0, 0}) =
      std::numeric_limits<double>::quiet_NaN();
  ibm_enthalpy.x.unchecked({3, 0, 0}) =
      std::numeric_limits<double>::quiet_NaN();
  const std::array<std::uint8_t, 3U> active_cells{{1U, 1U, 0U}};
  const std::array<std::uint8_t, 4U> active_x{{1U, 1U, 0U, 0U}};
  const std::array<std::uint8_t, 6U> active_y{{1U, 1U, 0U, 1U, 1U, 0U}};
  const std::array<std::uint8_t, 6U> active_z{{1U, 1U, 0U, 1U, 1U, 0U}};
  const PressureContinuityActivityView activity{
      {active_cells.data(), active_cells.size()},
      {active_x.data(), active_x.size()},
      {active_y.data(), active_y.size()},
      {active_z.data(), active_z.size()},
      8521U,
      8522U};
  OwnedField ibm_dp = ghosted_field(2U, 8523U, 8524U, cells, 1U);
  ibm_dp.view.unchecked({0, 0, 0}, 0U) = 1.0;
  ibm_dp.view.unchecked({1, 0, 0}, 0U) = 2.0;
  ibm_dp.view.unchecked({2, 0, 0}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  OwnedField ibm_output = shaped_field(89U, 8525U, 8526U, cells, 99.0);
  PressureEnergyPressureFluxOperator ibm_op;
  PressureEnergyPressureFluxCertificate ibm_certificate;
  const auto ibm_binding = pressure_flux_binding(
      ibm_fixture, as_const(ibm_temporal.view), ibm_coefficient, ibm_target,
      ibm_enthalpy, activity);
  passed &=
      expect(static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
                 ibm_binding, ibm_op, ibm_certificate)) &&
                 ibm_certificate.active_cells == 2U &&
                 ibm_certificate.inactive_cells == 1U &&
                 static_cast<bool>(ibm_op.apply(ibm_dp.view, ibm_output.view)),
             "IBM E_p ignores inactive NaN rows and interface payloads");
  const Vector ibm = values(as_const(ibm_output.view));
  passed &= expect(
      close(ibm[0U], -10.0) && close(ibm[1U], 10.0) && close(ibm[2U], 0.0),
      "IBM interface response is zero and solid row is closed");
  auto invalid_activity = activity;
  std::array<std::uint8_t, 4U> leaking_x = active_x;
  leaking_x[2U] = 1U;
  invalid_activity.x_faces = {leaking_x.data(), leaking_x.size()};
  auto leaking_binding = ibm_binding;
  leaking_binding.activity = invalid_activity;
  PressureEnergyPressureFluxOperator rejected;
  PressureEnergyPressureFluxCertificate rejected_certificate;
  passed &=
      expect(PressureEnergyPressureFluxOperator::bind(leaking_binding, rejected,
                                                      rejected_certificate)
                     .code == StatusCode::invalid_plan,
             "IBM solid cell cannot advertise a live pressure-energy face");
  return passed;
}

bool test_masked_solid_rows_rhs_sign_and_stale_components() {
  const Matrix identity_matrix{{
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
  }};
  const Matrix energy_pressure{{
      2.0, 0.0, 0.0,
      0.0, 3.0, 0.0,
      0.0, 0.0, 4.0,
  }};
  const Matrix energy_enthalpy{{
      5.0, 0.0, 0.0,
      0.0, 6.0, 0.0,
      0.0, 0.0, 7.0,
  }};
  const LinearIdentity identity{701U, 702U, 703U, 704U, 705U};
  DenseOperator cp(identity_matrix, identity, 711U, LinearOperatorClass::spd);
  DenseOperator ep(energy_pressure, identity, 712U,
                   LinearOperatorClass::nonsymmetric);
  DenseOperator eh(energy_enthalpy, identity, 713U,
                   LinearOperatorClass::nonsymmetric);
  OwnedField ch = field(41U, 51U, 1051U,
                        {{-2.0, 0.0, -4.0}});
  OwnedField scale = field(42U, 52U, 1052U,
                           {{2.0, 0.0, 4.0}});
  OwnedField cwork = field(43U, 53U, 1053U);
  OwnedField hwork = field(44U, 54U, 1054U);
  OwnedField ework = field(45U, 55U, 1055U);
  const std::array<std::uint8_t, kSize> active{{1U, 0U, 1U}};
  PressureEnergySchurBinding binding;
  binding.continuity_pressure = &cp;
  binding.energy_pressure = &ep;
  binding.energy_enthalpy = &eh;
  binding.continuity_enthalpy_diagonal = as_const(ch.view);
  binding.continuity_enthalpy_row_scale = as_const(scale.view);
  binding.workspace = {cwork.view, hwork.view, ework.view};
  binding.activity = {{active.data(), active.size()}, 721U, 722U};
  binding.scaled_pivot_floor = 1.0e-12;
  bool passed = expect(authorize_generic_schur(binding),
                       "masked dense Schur receives typed generic authority");
  PressureEnergySchurOperator schur;
  PressureEnergyJacobianCertificate certificate;
  passed &= expect(static_cast<bool>(PressureEnergySchurOperator::bind(
                       binding, schur, certificate)) &&
                           certificate.active_cells == 2U &&
                           certificate.inactive_cells == 1U,
                       "solid C_h=0 is excluded from the physical pivot gate");
  PressureEnergySchurBinding revised_binding = binding;
  ConstFieldView revised_ch = revised_binding.continuity_enthalpy_diagonal;
  ++revised_ch.revision;
  revised_binding.continuity_enthalpy_diagonal = revised_ch;
  PressureEnergySchurOperator revised_schur;
  PressureEnergyJacobianCertificate revised_certificate;
  passed &= expect(static_cast<bool>(PressureEnergySchurOperator::bind(
                           revised_binding, revised_schur,
                           revised_certificate)) &&
                       revised_certificate.schur.collective_fingerprint ==
                           certificate.schur.collective_fingerprint &&
                       revised_certificate.block_jacobian !=
                           certificate.block_jacobian,
                   "Schur keeps local C_h revision out of collective identity");

  OwnedField pressure = field(46U, 56U, 1056U, {{1.0, 9.0, -2.0}});
  OwnedField applied = field(47U, 57U, 1057U);
  passed &= expect(static_cast<bool>(schur.apply(pressure.view, applied.view)) &&
                       close(values(as_const(applied.view))[1U], 9.0),
                   "solid Schur row is dp=0 identity");

  const Vector pressure_system_rhs{{0.3, 88.0, -0.5}};  // equals -R_C
  const Vector energy_residual{{-1.0, 77.0, 2.0}};
  OwnedField bc = field(48U, 58U, 1058U, pressure_system_rhs);
  OwnedField re = field(49U, 59U, 1059U, energy_residual);
  OwnedField rhs = field(50U, 60U, 1060U);
  passed &= expect(static_cast<bool>(
                       schur.form_pressure_rhs_from_continuity_system_rhs(
                           as_const(bc.view), as_const(re.view), rhs.view)),
                   "Schur RHS consumes pressure_system.rhs=-R_C directly");
  const Vector formed = values(as_const(rhs.view));
  passed &= expect(close(formed[0U], 1.0 - 5.0 * 0.3 / -2.0) &&
                       close(formed[1U], 0.0) &&
                       close(formed[2U], -2.0 - 7.0 * -0.5 / -4.0),
                   "system-RHS convenience has the certified sign and zero solid row");

  OwnedField dp = field(51U, 61U, 1061U, {{0.2, 0.0, -0.1}});
  OwnedField dh = field(52U, 62U, 1062U);
  passed &= expect(static_cast<bool>(
                       schur.recover_enthalpy_from_continuity_system_rhs(
                           as_const(bc.view), dp.view, dh.view)),
                   "enthalpy recovery consumes the same pressure RHS");
  const Vector recovered = values(as_const(dh.view));
  passed &= expect(close(recovered[0U], (0.3 - 0.2) / -2.0) &&
                       close(recovered[1U], 0.0) &&
                       close(recovered[2U], (-0.5 - -0.1) / -4.0),
                   "system-RHS recovery preserves -R_C and closes solid dh=0");

  const std::array<std::uint8_t, kSize> all_solid{{0U, 0U, 0U}};
  PressureEnergySchurBinding solid_binding = binding;
  solid_binding.activity =
      {{all_solid.data(), all_solid.size()}, 731U, 732U};
  passed &= expect(authorize_generic_schur(solid_binding),
                   "all-solid dense Schur refreshes its typed activity authority");
  PressureEnergySchurOperator solid_schur;
  PressureEnergyJacobianCertificate solid_certificate;
  passed &= expect(static_cast<bool>(PressureEnergySchurOperator::bind(
                           solid_binding, solid_schur, solid_certificate)) &&
                       solid_certificate.valid() &&
                       solid_certificate.active_cells == 0U &&
                       close(solid_certificate.minimum_scaled_abs_c_h, 0.0) &&
                       static_cast<bool>(solid_schur.apply(pressure.view,
                                                           applied.view)) &&
                       values(as_const(applied.view)) ==
                           values(as_const(pressure.view)),
                   "an all-solid rank has a certified identity Schur block");

  LinearOperatorCertificate stale = ep.certificate();
  ++stale.collective_fingerprint;
  ep.replace_certificate(stale);
  const Status stale_apply = schur.apply(pressure.view, applied.view);
  passed &= expect(stale_apply.code == StatusCode::invalid_plan,
                   "a component revision refreshed after bind invalidates the Schur certificate");

  DenseOperator mixed(energy_pressure,
                      {identity.symbolic, identity.numeric + 1U,
                       identity.hierarchy, identity.workspace,
                       identity.fingerprint + 1U},
                      714U, LinearOperatorClass::nonsymmetric);
  binding.energy_pressure = &mixed;
  PressureEnergySchurOperator rejected;
  PressureEnergyJacobianCertificate rejected_certificate;
  const Status mixed_status = PressureEnergySchurOperator::bind(
      binding, rejected, rejected_certificate);
  passed &= expect(mixed_status.code == StatusCode::invalid_plan &&
                       !rejected_certificate.valid(),
                   "mixed complete LinearIdentity components cannot bind");
  return passed;
}

bool test_thermodynamic_tangent_and_frozen_mutations() {
  constexpr double pressure = 101325.0;
  constexpr double temperature = 300.0;
  constexpr double gas_constant = 287.05;
  constexpr double cp = 1004.5;
  constexpr double enthalpy = cp * temperature;
  ThermoState state;
  state.rho = pressure / (gas_constant * temperature);
  state.temperature = temperature;
  state.cp = cp;
  state.gas_constant = gas_constant;
  state.drho_dp_hY = 1.0 / (gas_constant * temperature);
  state.drho_dh_pY = -state.rho / (temperature * cp);
  PressureEnergyThermoJacobian jacobian;
  bool passed = expect(static_cast<bool>(form_pressure_energy_thermo_jacobian(
                           pressure, enthalpy, state, jacobian)),
                       "pressure-energy EOS tangent forms");
  passed &= expect(close(jacobian.drho_dp_hY, state.drho_dp_hY) &&
                       close(jacobian.drho_dh_pY, state.drho_dh_pY) &&
                       close(jacobian.dq_dp_hY,
                             enthalpy * state.drho_dp_hY - 1.0) &&
                       close(jacobian.dq_dh_pY,
                             state.rho + enthalpy * state.drho_dh_pY),
                   "q=rho*h-p Jacobian uses the certified EOS tangent");

  constexpr double delta_pressure = 250.0;
  constexpr double a0 = 250.0;
  const double frozen_h_energy_defect =
      a0 * jacobian.dq_dp_hY * delta_pressure;
  passed &= expect(std::abs(frozen_h_energy_defect) > 1.0e4,
                   "frozen-h pressure correction leaves target energy RED");

  const double rho_star = state.rho;
  const double artificial_delta_h = delta_pressure / rho_star;
  const double p_next = pressure + delta_pressure;
  const double h_next = enthalpy + artificial_delta_h;
  const double rho_next = p_next * cp / (gas_constant * h_next);
  const double q_before = state.rho * enthalpy - pressure;
  const double q_after = rho_next * h_next - p_next;
  passed &= expect(close(q_after - q_before,
                         (cp / gas_constant - 1.0) * delta_pressure,
                         5.0e-13) &&
                       std::abs(q_after - q_before) > 1.0,
                   "frozen-rho delta_h=delta_p/rho* fails constant-cp q");
  return passed;
}

bool test_pressure_rate_modal_oracle() {
  const auto parasitic_radius = [](double gain) {
    const std::complex<double> discriminant{gain * (gain - 1.0), 0.0};
    const std::complex<double> root = std::sqrt(discriminant);
    return std::max(std::abs(std::complex<double>{gain, 0.0} + root),
                    std::abs(std::complex<double>{gain, 0.0} - root));
  };
  bool passed = expect(parasitic_radius(0.285) < 1.0,
                       "uniform ideal-gas gain has damped EX2 pressure roots");
  passed &= expect(close(parasitic_radius(-1.0 / 3.0), 1.0, 2.0e-14),
                   "old BDF(Dp/Dt)-EX2 loop reaches its stability edge");
  passed &= expect(parasitic_radius(-0.5) > 1.0 && std::abs(-0.5) < 1.0,
                   "old loop is unstable where lagged BE remains stable");
  // The direct target-time block has only the physical BDF2 factor
  // 3*r^2-4*r+1=(r-1)(3*r-1); it has no gain-dependent quadratic.
  passed &= expect(close(1.0, 1.0) && close(1.0 / 3.0, 1.0 / 3.0),
                   "target-time block removes the gain-dependent roots");
  return passed;
}

ThermodynamicsPlan constant_cp_thermodynamics(double cp,
                                              double gas_constant) {
  ThermophysicalSpec spec;
  spec.data_file = "pressure-energy-analytic.d";
  spec.minimum_temperature = 200.0;
  spec.maximum_temperature = 2000.0;
  spec.temperature_relative_tolerance = 1.0e-12;
  spec.maximum_temperature_iterations = 64U;
  spec.closed_mass_relative_tolerance = 1.0e-12;
  spec.maximum_closed_mass_iterations = 32U;
  spec.maximum_closed_mass_relative_step = 0.2;
  SpeciesThermophysicalSpec species;
  species.stable_name = "gas";
  species.molecular_weight = kUniversalGasConstant / gas_constant;
  species.temperature_switch = 1000.0;
  species.nasa7_low[0U] = cp / gas_constant;
  species.nasa7_high[0U] = cp / gas_constant;
  species.viscosity_reference = 1.8e-5;
  species.conductivity = 0.026;
  spec.species.push_back(species);
  ThermodynamicsPlan plan;
  const Status status = ThermodynamicsPlan::compile(spec, {}, plan);
  if (!status) std::cerr << "FAIL: constant-cp temporal plan compiles\n";
  return plan;
}

BdfCoefficients bdf(double dt, double ratio) {
  if (ratio == 0.0) return {1.0 / dt, -1.0 / dt, 0.0, 1U};
  return {(1.0 + 2.0 * ratio) / ((1.0 + ratio) * dt),
          -(1.0 + ratio) / dt,
          ratio * ratio / ((1.0 + ratio) * dt), 2U};
}

bool test_target_time_temporal_closed_loop() {
  constexpr double cp = 1004.5;
  constexpr double gas = 287.05;
  ThermodynamicsPlan thermodynamics = constant_cp_thermodynamics(cp, gas);
  bool passed = thermodynamics.fingerprint() != 0U;
  for (const double ratio : {0.0, 0.5, 1.0, 2.0}) {
    const BdfCoefficients coefficients = bdf(0.02, ratio);
    const double p_n = 101325.0;
    const double h_n = cp * 300.0;
    const double p_nm1 = 100900.0;
    const double h_nm1 = cp * 299.7;
    ThermoState accepted;
    ThermoState previous;
    passed &= expect(static_cast<bool>(thermodynamics.evaluate(
                         p_n, h_n, {}, {}, accepted)) &&
                         (coefficients.order == 1U ||
                          static_cast<bool>(thermodynamics.evaluate(
                              p_nm1, h_nm1, {}, {}, previous))),
                     "temporal history closes through the production EOS");
    const double rho_target =
        -(coefficients.a1 * accepted.rho +
          coefficients.a2 *
              (coefficients.order == 2U ? previous.rho : 0.0)) /
        coefficients.a0;
    const double p_target =
        -(coefficients.a1 * p_n +
          coefficients.a2 *
              (coefficients.order == 2U ? p_nm1 : 0.0)) /
        coefficients.a0;
    const double h_target = p_target * cp / (gas * rho_target);
    double pressure = p_target * 1.0001;
    double enthalpy = h_target * 0.9998;
    double first_norm = 0.0;
    for (std::uint8_t corrector = 1U; corrector <= 2U; ++corrector) {
      ThermoState state;
      PressureEnergyThermoJacobian thermo;
      passed &= expect(static_cast<bool>(thermodynamics.evaluate(
                           pressure, enthalpy, {}, {}, state)) &&
                           static_cast<bool>(
                               form_pressure_energy_thermo_jacobian(
                                   pressure, enthalpy, state, thermo)),
                       "corrector recomputes the target EOS tangent");
      PressureEnergyTemporalPoint point;
      point.bdf = coefficients;
      point.cell_volume = 0.125;
      point.pressure_absolute = pressure;
      point.density = state.rho;
      point.enthalpy = enthalpy;
      point.accepted_pressure_absolute = p_n;
      point.accepted_density = accepted.rho;
      point.accepted_enthalpy = h_n;
      point.previous_pressure_absolute =
          coefficients.order == 2U
              ? p_nm1
              : std::numeric_limits<double>::quiet_NaN();
      point.previous_density =
          coefficients.order == 2U
              ? previous.rho
              : std::numeric_limits<double>::quiet_NaN();
      point.previous_enthalpy =
          coefficients.order == 2U
              ? h_nm1
              : std::numeric_limits<double>::quiet_NaN();
      point.target_thermo = thermo;
      PressureEnergyTemporalLinearization temporal;
      passed &= expect(static_cast<bool>(
                           linearize_pressure_energy_temporal(point,
                                                              temporal)),
                       "target-time BDF(rho*h-p) linearizes");
      const double norm =
          std::max(std::abs(temporal.continuity_residual),
                   std::abs(temporal.energy_residual));
      if (corrector == 1U) first_norm = norm;
      const double determinant =
          temporal.continuity_pressure * temporal.energy_enthalpy -
          temporal.continuity_enthalpy * temporal.energy_pressure;
      passed &= expect(std::isfinite(determinant) &&
                           std::abs(determinant) > 1.0e-12,
                       "temporal p-h block is nonsingular");
      const double delta_pressure =
          (-temporal.continuity_residual * temporal.energy_enthalpy +
           temporal.continuity_enthalpy * temporal.energy_residual) /
          determinant;
      const double delta_enthalpy =
          (-temporal.continuity_pressure * temporal.energy_residual +
           temporal.continuity_residual * temporal.energy_pressure) /
          determinant;
      pressure += delta_pressure;
      enthalpy += delta_enthalpy;
    }
    ThermoState final_state;
    PressureEnergyThermoJacobian final_thermo;
    passed &= expect(static_cast<bool>(thermodynamics.evaluate(
                         pressure, enthalpy, {}, {}, final_state)) &&
                         static_cast<bool>(form_pressure_energy_thermo_jacobian(
                             pressure, enthalpy, final_state, final_thermo)),
                     "C2 state recloses through EOS");
    PressureEnergyTemporalPoint final_point;
    final_point.bdf = coefficients;
    final_point.cell_volume = 0.125;
    final_point.pressure_absolute = pressure;
    final_point.density = final_state.rho;
    final_point.enthalpy = enthalpy;
    final_point.accepted_pressure_absolute = p_n;
    final_point.accepted_density = accepted.rho;
    final_point.accepted_enthalpy = h_n;
    final_point.previous_pressure_absolute =
        coefficients.order == 2U
            ? p_nm1
            : std::numeric_limits<double>::quiet_NaN();
    final_point.previous_density =
        coefficients.order == 2U
            ? previous.rho
            : std::numeric_limits<double>::quiet_NaN();
    final_point.previous_enthalpy =
        coefficients.order == 2U
            ? h_nm1
            : std::numeric_limits<double>::quiet_NaN();
    final_point.target_thermo = final_thermo;
    PressureEnergyTemporalLinearization final_temporal;
    passed &= expect(static_cast<bool>(linearize_pressure_energy_temporal(
                         final_point, final_temporal)),
                     "terminal target-time residual evaluates");
    const double final_norm =
        std::max(std::abs(final_temporal.continuity_residual),
                 std::abs(final_temporal.energy_residual));
    passed &= expect(final_norm < first_norm * 1.0e-6 &&
                         close(pressure, p_target, 2.0e-12) &&
                         close(enthalpy, h_target, 2.0e-8),
                     "two target-time corrections close pressure-energy-EOS-continuity");
  }
  return passed;
}

bool test_collective_identity_is_decomposition_independent() {
  const LinearIdentity wide_identity{901U, 902U, 903U, 904U, 905U};
  const LinearIdentity narrow_identity{901U, 912U, 903U, 914U, 915U};
  OwnedField wide = shaped_field(71U, 72U, 73U, {5, 3, 2}, 2.0);
  OwnedField narrow = shaped_field(71U, 72U, 74U, {4, 3, 2}, 2.0);
  narrow.view.revision = 75U;
  narrow.view.replica = 1U;
  PressureEnergyDiagonalOperator wide_operator;
  PressureEnergyDiagonalOperator narrow_operator;
  PressureEnergyDiagonalCertificate wide_certificate;
  PressureEnergyDiagonalCertificate narrow_certificate;
  const Status wide_status = PressureEnergyDiagonalOperator::bind(
      {as_const(wide.view), {}, wide_identity, 0.0}, wide_operator,
      wide_certificate);
  const Status narrow_status = PressureEnergyDiagonalOperator::bind(
      {as_const(narrow.view), {}, narrow_identity, 0.0}, narrow_operator,
      narrow_certificate);
  if (!wide_status || !narrow_status)
    std::cerr << "decomposition identity bind status="
              << static_cast<unsigned>(wide_status.code) << ':'
              << wide_status.detail << ','
              << static_cast<unsigned>(narrow_status.code) << ':'
              << narrow_status.detail << '\n';
  bool passed = expect(
      static_cast<bool>(wide_status) && static_cast<bool>(narrow_status),
      "unequal local pressure-energy partitions bind independently");
  passed &= expect(
      wide_certificate.linear.collective_fingerprint ==
              narrow_certificate.linear.collective_fingerprint &&
          wide_certificate.linear.identity.numeric !=
              narrow_certificate.linear.identity.numeric &&
          wide_certificate.linear.identity.workspace !=
              narrow_certificate.linear.identity.workspace &&
          wide_certificate.linear.local_shape.x == 5 &&
          narrow_certificate.linear.local_shape.x == 4,
      "collective pressure-energy identity excludes decomposition-local shape");

  PressureFluxFixture fixture;
  passed &= expect(make_pressure_flux_fixture(
                       PressureFluxBoundary::neumann_dirichlet, fixture),
                   "pressure-flux collective identity fixture compiles");
  if (!passed) return false;
  const Int3 cells = fixture.patch.cells;
  OwnedField temporal_a = shaped_field(90U, 73U, 8601U, cells, 2.0);
  OwnedField temporal_b = shaped_field(90U, 73U, 8602U, cells, 2.0);
  temporal_b.view.revision = 8615U;
  temporal_b.view.replica = 1U;
  OwnedFaces coefficient_a = face_bundle(cells, 8603U, 8604U, 1.0);
  OwnedFaces coefficient_b = face_bundle(cells, 8605U, 8606U, 1.0);
  OwnedFaces target_a = face_bundle(cells, 8607U, 8608U, 0.0);
  OwnedFaces target_b = face_bundle(cells, 8609U, 8610U, 0.0);
  OwnedFaces enthalpy_a = face_bundle(cells, 8611U, 8612U, 300000.0);
  OwnedFaces enthalpy_b = face_bundle(cells, 8613U, 8614U, 300000.0);
  const LinearIdentity identity_a{911U, 912U, 913U, 914U, 915U};
  const LinearIdentity identity_b{911U, 922U, 913U, 924U, 925U};
  const auto binding_a =
      pressure_flux_binding(fixture, as_const(temporal_a.view), coefficient_a,
                            target_a, enthalpy_a, {}, identity_a);
  auto binding_b =
      pressure_flux_binding(fixture, as_const(temporal_b.view), coefficient_b,
                            target_b, enthalpy_b, {}, identity_b);
  // Target flux and frozen-face revisions are rank-local producer tokens in
  // production.  Their exact values belong to binding_revision/local_binding,
  // while collective identity may only certify their non-zero presence.
  binding_b.target_flux.revision = 8616U;
  binding_b.frozen_face_enthalpy.revision = 8617U;
  binding_b.intermediate.plan = 8618U;
  binding_b.pressure.plan = binding_b.intermediate.plan;
  ++binding_b.intermediate
        .thermophysical_boundary_rank_local_binding;
  binding_b.pressure.thermophysical_boundary_rank_local_binding =
      binding_b.intermediate
          .thermophysical_boundary_rank_local_binding;
  ++binding_b.intermediate.thermophysical_boundary_rank_local_lineage;
  binding_b.pressure.thermophysical_boundary_rank_local_lineage =
      binding_b.intermediate.thermophysical_boundary_rank_local_lineage;
  binding_b.frozen_face_enthalpy.local_binding =
      pressure_energy_frozen_face_enthalpy_local_binding(
          binding_b.frozen_face_enthalpy);
  PressureEnergyPressureFluxOperator operator_a;
  PressureEnergyPressureFluxOperator operator_b;
  PressureEnergyPressureFluxCertificate certificate_a;
  PressureEnergyPressureFluxCertificate certificate_b;
  passed &= expect(
      static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
          binding_a, operator_a, certificate_a)) &&
          static_cast<bool>(PressureEnergyPressureFluxOperator::bind(
              binding_b, operator_b, certificate_b)) &&
          certificate_a.linear.collective_fingerprint ==
              certificate_b.linear.collective_fingerprint &&
          certificate_a.binding_revision != certificate_b.binding_revision &&
          certificate_a.frozen_face_enthalpy_local_binding !=
              certificate_b.frozen_face_enthalpy_local_binding &&
          certificate_a.linear.identity.numeric !=
              certificate_b.linear.identity.numeric &&
          certificate_a.linear.identity.workspace !=
              certificate_b.linear.identity.workspace,
      "pressure-flux collective identity excludes rank-local storage and "
      "replicas");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
  bool passed = test_exact_schur_and_recovery();
  passed &= test_enthalpy_spatial_binding_rejects_an_empty_contract();
  passed &= test_enthalpy_spatial_target_contract_binds();
  passed &= test_enthalpy_spatial_periodic_mpi_and_inactive_interfaces();
  passed &= test_enthalpy_semismooth_limiter_certificate();
  passed &= test_diagonal_operator_activity_and_identity();
  passed &= test_ibm_double_diagonal_typed_schur_authority();
  passed &= test_mass_flow_three_cell_pressure_flux_red();
  passed &= test_boundary_constant_h_and_directional_derivative();
  passed &= test_periodic_and_ibm_pressure_flux_semantics();
  passed &=
      test_analytic_frozen_target_cartesian_four_block_fd_certificate();
  passed &= test_masked_solid_rows_rhs_sign_and_stale_components();
  passed &= test_thermodynamic_tangent_and_frozen_mutations();
  passed &= test_pressure_rate_modal_oracle();
  passed &= test_target_time_temporal_closed_loop();
  passed &= test_collective_identity_is_decomposition_independent();
  MPI_Finalize();
  return passed ? 0 : 1;
}
