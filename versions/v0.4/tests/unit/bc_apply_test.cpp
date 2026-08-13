// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include <mpi.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <utility>
#include <vector>

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

constexpr Int3 kInterior{4, 3, 2};
constexpr Int3 kGhosts{2, 2, 2};
constexpr double kSentinel = -987654.25;
constexpr std::size_t kCanaryDoubles = 16U;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

struct OwnedView {
  std::vector<double> storage;
  FieldView view;

  OwnedView(FieldId field, std::uint8_t components,
            Int3 interior = kInterior, Int3 ghosts = kGhosts)
      : storage(), view{} {
    const std::size_t extent_x = static_cast<std::size_t>(
        interior.x + 2 * ghosts.x);
    const std::size_t extent_y = static_cast<std::size_t>(
        interior.y + 2 * ghosts.y);
    const std::size_t extent_z = static_cast<std::size_t>(
        interior.z + 2 * ghosts.z);
    const std::size_t stride_y = extent_x;
    const std::size_t stride_z = stride_y * extent_y;
    const std::size_t component_stride = stride_z * extent_z;
    storage.assign(kCanaryDoubles +
                       component_stride * static_cast<std::size_t>(components) +
                       kCanaryDoubles,
                   kSentinel);
    view.base = storage.data() + kCanaryDoubles +
                static_cast<std::size_t>(ghosts.x) +
                static_cast<std::size_t>(ghosts.y) * stride_y +
                static_cast<std::size_t>(ghosts.z) * stride_z;
    view.interior = interior;
    view.ghosts = ghosts;
    view.components = components;
    view.stride_y = stride_y;
    view.stride_z = stride_z;
    view.component_stride = component_stride;
    view.field = field;
  }
};

BoundaryFaceSpec wall() {
  BoundaryFaceSpec value;
  value.flow_kind = BoundaryKind::no_slip_wall;
  value.thermal_kind = BoundaryKind::adiabatic_wall;
  value.mach_limit = 0.95;
  return value;
}

ValidatedModel model() {
  ValidatedModel value;
  value.fingerprint = 0xA991U;
  value.pressure_reference = PressureReferenceKind::boundary_absolute;
  value.transported_scalars.push_back(
      TransportedScalarSpec{"mixture_fraction",
                            TransportedScalarRole::passive_scalar});
  value.transported_scalars.push_back(
      TransportedScalarSpec{"progress",
                            TransportedScalarRole::passive_scalar});
  for (BoundaryFaceSpec& face : value.boundaries) {
    face = wall();
    face.scalars.push_back(
        ScalarBoundarySpec{"mixture_fraction",
                           ScalarBoundaryKind::zero_gradient, 0.0});
    face.scalars.push_back(
        ScalarBoundarySpec{"progress", ScalarBoundaryKind::zero_gradient,
                           0.0});
  }
  BoundaryFaceSpec& inlet = value.boundaries[0];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.thermal_kind = BoundaryKind::none;
  inlet.velocity = Real3{2.0, 3.0, 4.0};
  inlet.temperature = 300.0;
  inlet.scalars[0] =
      ScalarBoundarySpec{"mixture_fraction", ScalarBoundaryKind::dirichlet,
                         0.25};
  inlet.scalars[1] =
      ScalarBoundarySpec{"progress", ScalarBoundaryKind::normal_flux, 0.5};

  BoundaryFaceSpec& outlet = value.boundaries[1];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.thermal_kind = BoundaryKind::none;
  outlet.pressure = 101325.0;
  outlet.scalars[0].kind = ScalarBoundaryKind::convective;
  BoundaryFaceSpec& isothermal_wall = value.boundaries[2];
  isothermal_wall.thermal_kind = BoundaryKind::isothermal_wall;
  isothermal_wall.temperature = 425.0;
  BoundaryFaceSpec& mass_inlet = value.boundaries[4];
  mass_inlet.flow_kind = BoundaryKind::mass_flow_inlet;
  mass_inlet.thermal_kind = BoundaryKind::none;
  mass_inlet.mass_flow_rate = 1.0;
  mass_inlet.temperature = 310.0;
  mass_inlet.direction = Real3{0.0, 0.0, 1.0};
  value.schemes = SchemeSpec{};
  value.time = TimeControlSpec{};
  return value;
}

bool compile_plan(BoundaryPlan& plan) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = Real3{0.0, 0.0, 0.0};
  mesh.upper = Real3{4.0, 3.0, 2.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = kInterior;
  mesh.minimum_spacing = Real3{1.0, 1.0, 1.0};
  mesh.max_growth_ratio = 1.0;
  mesh.limits.max_global_cells = 1024U;
  mesh.limits.max_memory_bytes_per_rank = 1U << 20U;
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch)) {
    return false;
  }
  FieldRegistry registry;
  SchemePlan schemes;
  TimeSchemePlan time;
  return static_cast<bool>(BoundaryCompiler::compile(
      MPI_COMM_SELF, model(), geometry, patch, registry, plan, schemes, time));
}

void initialise_interior(OwnedView& field) {
  for (std::uint8_t component = 0U; component < field.view.components;
       ++component) {
    for (std::int32_t z = 0; z < field.view.interior.z; ++z) {
      for (std::int32_t y = 0; y < field.view.interior.y; ++y) {
        for (std::int32_t x = 0; x < field.view.interior.x; ++x) {
          field.view.unchecked(Int3{x, y, z}, component) =
              1000.0 * component + 100.0 * z + 10.0 * y + x + 1.0;
        }
      }
    }
  }
}

bool unchanged(const OwnedView& field,
               const std::vector<double>& before) {
  return field.storage == before;
}

bool canaries_intact(const OwnedView& field) {
  return std::all_of(field.storage.begin(),
                     field.storage.begin() +
                         static_cast<std::ptrdiff_t>(kCanaryDoubles),
                     [](double value) { return value == kSentinel; }) &&
         std::all_of(field.storage.end() -
                         static_cast<std::ptrdiff_t>(kCanaryDoubles),
                     field.storage.end(),
                     [](double value) { return value == kSentinel; });
}

const BoundaryIndexSpan* find_span(const BoundaryPlan& plan,
                                   BoundaryStage stage,
                                   CartesianFace face, FieldId field) {
  for (std::size_t index = 0U; index < plan.spans().size; ++index) {
    const BoundaryIndexSpan& span = plan.spans().data[index];
    if (span.stage == stage && span.face == face && span.field == field) {
      return &span;
    }
  }
  return nullptr;
}

FieldId scalar_field(const BoundaryPlan& plan, CartesianFace face,
                     ScalarBoundaryKind kind) {
  for (std::size_t index = 0U; index < plan.spans().size; ++index) {
    const BoundaryIndexSpan& span = plan.spans().data[index];
    if (span.stage == BoundaryStage::scalar && span.face == face &&
        span.parameter < plan.scalar_kinds().size &&
        plan.scalar_kinds().data[span.parameter] == kind) {
      return span.field;
    }
  }
  return std::numeric_limits<FieldId>::max();
}

BoundaryResolvedValues resolved_values(
    const BoundaryPlan& plan, std::vector<double>& scalar,
    std::vector<Real3>& vector, std::vector<double>& gradient) {
  scalar.assign(plan.resolved_scalar_count(), 0.0);
  vector.assign(plan.resolved_vector_count(), Real3{});
  gradient.assign(plan.resolved_normal_gradient_count(), 0.0);
  return BoundaryResolvedValues{
      Span<const double>{scalar.data(), scalar.size()},
      Span<const Real3>{vector.data(), vector.size()},
      Span<const double>{gradient.data(), gradient.size()}};
}

std::size_t resolved_index(const BoundaryIndexSpan& span,
                           std::uint32_t inner,
                           std::uint32_t outer) {
  return static_cast<std::size_t>(span.resolved_begin) +
         static_cast<std::size_t>(outer) * span.tangent_inner_count + inner;
}

template <class Value>
void fill_resolved_slice(std::vector<Value>& values,
                         const BoundaryIndexSpan& span,
                         const Value& value) {
  const std::size_t begin = span.resolved_begin;
  const std::size_t end = begin + span.resolved_stride;
  if (end <= values.size()) {
    std::fill(values.begin() + static_cast<std::ptrdiff_t>(begin),
              values.begin() + static_cast<std::ptrdiff_t>(end), value);
  }
}

std::vector<FieldId> scalar_fields(const BoundaryPlan& plan) {
  std::vector<FieldId> fields;
  for (std::size_t index = 0U; index < plan.spans().size; ++index) {
    const BoundaryIndexSpan& span = plan.spans().data[index];
    if (span.stage == BoundaryStage::scalar &&
        std::find(fields.begin(), fields.end(), span.field) == fields.end()) {
      fields.push_back(span.field);
    }
  }
  return fields;
}

bool test_empty_plan() {
  BoundaryPlan plan;
  OwnedView field{7U, 1U};
  const std::vector<double> before = field.storage;
  const Status status = apply_boundary_ghosts(
      BoundaryStage::momentum, plan, Span<FieldView>{&field.view, 1U}, {});
  return expect(static_cast<bool>(status),
                "empty compiled stage is an exact no-op") &&
         expect(unchanged(field, before),
                "empty stage leaves every storage byte untouched");
}

bool test_stage_atomic_when_second_field_is_missing(const BoundaryPlan& plan) {
  const std::vector<FieldId> ids = scalar_fields(plan);
  bool passed = expect(ids.size() == 2U,
                       "fixture compiles two distinct scalar fields");
  if (ids.size() != 2U) {
    return false;
  }
  OwnedView first{ids[0], 1U};
  initialise_interior(first);
  const std::vector<double> before = first.storage;
  std::vector<double> scalar;
  std::vector<Real3> vector;
  std::vector<double> gradient;
  const BoundaryResolvedValues resolved =
      resolved_values(plan, scalar, vector, gradient);
  const Status status = apply_boundary_ghosts(
      BoundaryStage::scalar, plan, Span<FieldView>{&first.view, 1U},
      resolved);
  passed &= expect(status.code == StatusCode::invalid_plan,
                   "missing second scalar field rejects the whole stage");
  passed &= expect(unchanged(first, before),
                   "late span failure cannot partially mutate an earlier span");
  passed &= expect(canaries_intact(first),
                   "late span failure preserves storage canaries");
  return passed;
}

bool test_view_preflight_is_atomic(const BoundaryPlan& plan) {
  bool passed = true;
  OwnedView velocity{plan.velocity_field(), 3U};
  initialise_interior(velocity);
  std::vector<double> scalar;
  std::vector<Real3> vector;
  std::vector<double> gradient;
  const BoundaryResolvedValues resolved =
      resolved_values(plan, scalar, vector, gradient);

  const auto rejects_without_writes = [&](FieldView malformed,
                                          std::string_view reason) {
    const std::vector<double> before = velocity.storage;
    const Status status = apply_boundary_ghosts(
        BoundaryStage::momentum, plan, Span<FieldView>{&malformed, 1U},
        resolved);
    passed &= expect(status.code == StatusCode::invalid_plan, reason);
    passed &= expect(unchanged(velocity, before),
                     "rejected view leaves all field bytes unchanged");
    passed &= expect(canaries_intact(velocity),
                     "rejected view preserves front and rear canaries");
  };

  FieldView wrong_shape = velocity.view;
  --wrong_shape.interior.y;
  rejects_without_writes(wrong_shape,
                         "tangential shape mismatch is rejected pre-write");

  FieldView shallow_normal = velocity.view;
  shallow_normal.interior.x = 1;
  rejects_without_writes(shallow_normal,
                         "normal depth smaller than stencil is rejected");

  FieldView wrong_stride = velocity.view;
  wrong_stride.stride_y = 1U;
  rejects_without_writes(wrong_stride,
                         "overlapping row stride is rejected pre-write");

  FieldView wrong_plane_stride = velocity.view;
  wrong_plane_stride.stride_z = wrong_plane_stride.stride_y;
  rejects_without_writes(wrong_plane_stride,
                         "overlapping plane stride is rejected pre-write");

  FieldView wrong_component_stride = velocity.view;
  wrong_component_stride.component_stride = wrong_component_stride.stride_z;
  rejects_without_writes(wrong_component_stride,
                         "overlapping component stride is rejected pre-write");

  const std::size_t ptrdiff_max = static_cast<std::size_t>(
      std::numeric_limits<std::ptrdiff_t>::max());
  FieldView unrepresentable_component_stride = velocity.view;
  unrepresentable_component_stride.component_stride = ptrdiff_max + 1U;
  rejects_without_writes(
      unrepresentable_component_stride,
      "component stride outside ptrdiff_t is rejected pre-write");

  FieldView overflowing_affine_offset = velocity.view;
  overflowing_affine_offset.component_stride = ptrdiff_max / 2U + 1U;
  rejects_without_writes(
      overflowing_affine_offset,
      "combined ghosted component affine offset is rejected pre-write");
  return passed;
}

bool test_resolved_shape_preflight_is_atomic(const BoundaryPlan& plan) {
  OwnedView velocity{plan.velocity_field(), 3U};
  initialise_interior(velocity);
  const std::vector<double> before = velocity.storage;
  std::vector<Real3> vector(plan.resolved_vector_count(), Real3{});
  BoundaryResolvedValues missing;
  Status status = apply_boundary_ghosts(
      BoundaryStage::momentum, plan, Span<FieldView>{&velocity.view, 1U},
      missing);
  bool passed = expect(status.code == StatusCode::invalid_plan,
                       "missing current-stage resolved buffer is rejected") &&
                expect(unchanged(velocity, before),
                       "missing resolved buffer rejects before mutation");
  BoundaryResolvedValues short_values;
  short_values.vector = Span<const Real3>{
      vector.data(), vector.empty() ? 0U : vector.size() - 1U};
  status = apply_boundary_ghosts(
      BoundaryStage::momentum, plan, Span<FieldView>{&velocity.view, 1U},
      short_values);
  passed &= expect(status.code == StatusCode::invalid_plan,
                   "short current-stage resolved arrays are rejected");
  passed &= expect(unchanged(velocity, before),
                   "resolved parameter shape rejects before field mutation");
  passed &= expect(canaries_intact(velocity),
                   "resolved parameter failure preserves canaries");

  OwnedView pressure{plan.pressure_field(), 1U};
  initialise_interior(pressure);
  const std::vector<double> before_pressure = pressure.storage;
  status = apply_boundary_ghosts(
      BoundaryStage::pressure, plan, Span<FieldView>{&pressure.view, 1U}, {});
  passed &= expect(status.code == StatusCode::invalid_plan,
                   "missing current-stage scalar buffer is rejected");
  passed &= expect(unchanged(pressure, before_pressure),
                   "missing scalar buffer rejection is atomic");
  const BoundaryIndexSpan* pressure_span = find_span(
      plan, BoundaryStage::pressure, CartesianFace::x_max,
      plan.pressure_field());
  passed &= expect(pressure_span != nullptr,
                   "short-scalar fixture has a resolved pressure span");
  if (pressure_span != nullptr) {
    const std::size_t required =
        static_cast<std::size_t>(pressure_span->resolved_begin) +
        pressure_span->resolved_stride;
    std::vector<double> short_scalar(required - 1U, 0.0);
    BoundaryResolvedValues short_scalar_values;
    short_scalar_values.scalar =
        Span<const double>{short_scalar.data(), short_scalar.size()};
    status = apply_boundary_ghosts(
        BoundaryStage::pressure, plan,
        Span<FieldView>{&pressure.view, 1U}, short_scalar_values);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "short current-stage scalar buffer is rejected");
    passed &= expect(unchanged(pressure, before_pressure),
                     "short scalar buffer rejection is atomic");
  }

  const std::vector<FieldId> ids = scalar_fields(plan);
  OwnedView first{ids[0], 1U};
  OwnedView second{ids[1], 1U};
  initialise_interior(first);
  initialise_interior(second);
  const std::vector<double> before_first = first.storage;
  const std::vector<double> before_second = second.storage;
  FieldView fields[]{first.view, second.view};
  std::vector<double> scalar(plan.resolved_scalar_count(), 0.0);
  BoundaryResolvedValues missing_gradient;
  missing_gradient.scalar =
      Span<const double>{scalar.data(), scalar.size()};
  status = apply_boundary_ghosts(BoundaryStage::scalar, plan,
                                 Span<FieldView>{fields, 2U},
                                 missing_gradient);
  passed &= expect(status.code == StatusCode::invalid_plan,
                   "missing current-stage gradient buffer is rejected");
  passed &= expect(unchanged(first, before_first) &&
                       unchanged(second, before_second),
                   "missing gradient rejects the complete stage atomically");
  const FieldId flux_field = scalar_field(
      plan, CartesianFace::x_min, ScalarBoundaryKind::normal_flux);
  const BoundaryIndexSpan* flux_span = find_span(
      plan, BoundaryStage::scalar, CartesianFace::x_min, flux_field);
  passed &= expect(flux_span != nullptr,
                   "short-gradient fixture has a resolved flux span");
  if (flux_span != nullptr) {
    const std::size_t required =
        static_cast<std::size_t>(flux_span->resolved_begin) +
        flux_span->resolved_stride;
    std::vector<double> short_gradient(required - 1U, 0.0);
    BoundaryResolvedValues short_gradient_values;
    short_gradient_values.scalar =
        Span<const double>{scalar.data(), scalar.size()};
    short_gradient_values.normal_gradient = Span<const double>{
        short_gradient.data(), short_gradient.size()};
    status = apply_boundary_ghosts(BoundaryStage::scalar, plan,
                                   Span<FieldView>{fields, 2U},
                                   short_gradient_values);
    passed &= expect(status.code == StatusCode::invalid_plan,
                     "short current-stage gradient buffer is rejected");
    passed &= expect(unchanged(first, before_first) &&
                         unchanged(second, before_second),
                     "short gradient rejects the complete stage atomically");
  }
  return passed;
}

bool test_resolved_field_space_relations(const BoundaryPlan& plan) {
  std::vector<double> resolved_scalar;
  std::vector<Real3> resolved_vector;
  std::vector<double> resolved_gradient;
  BoundaryResolvedValues resolved = resolved_values(
      plan, resolved_scalar, resolved_vector, resolved_gradient);
  bool passed = true;

  const BoundaryIndexSpan* inlet_enthalpy_span = find_span(
      plan, BoundaryStage::enthalpy, CartesianFace::x_min,
      plan.enthalpy_field());
  passed &= expect(inlet_enthalpy_span != nullptr &&
                       inlet_enthalpy_span->relation ==
                           BoundaryRelation::dirichlet &&
                       inlet_enthalpy_span->value_source ==
                           BoundaryValueSource::resolved_scalar,
                   "velocity inlet compiles a field-space enthalpy slot");
  if (inlet_enthalpy_span != nullptr &&
      inlet_enthalpy_span->resolved_begin < resolved_scalar.size()) {
    constexpr double kInletEnthalpyTarget = 321000.0;
    fill_resolved_slice(resolved_scalar, *inlet_enthalpy_span,
                        kInletEnthalpyTarget);
    OwnedView enthalpy{plan.enthalpy_field(), 1U};
    initialise_interior(enthalpy);
    passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                         BoundaryStage::enthalpy, plan,
                         Span<FieldView>{&enthalpy.view, 1U}, resolved)),
                     "resolved velocity-inlet enthalpy stage applies");
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior =
          enthalpy.view.unchecked(Int3{layer - 1, 1, 1}, 0U);
      const double ghost =
          enthalpy.view.unchecked(Int3{-layer, 1, 1}, 0U);
      passed &= expect(ghost == 2.0 * kInletEnthalpyTarget - interior,
                       "velocity inlet uses resolved h from its T/Y state");
    }
    passed &= expect(canaries_intact(enthalpy),
                     "resolved inlet enthalpy preserves canaries");
  }

  const BoundaryIndexSpan* pressure_span = find_span(
      plan, BoundaryStage::pressure, CartesianFace::x_max,
      plan.pressure_field());
  passed &= expect(pressure_span != nullptr &&
                       pressure_span->value_source ==
                           BoundaryValueSource::resolved_scalar,
                   "pressure outlet compiles a field-space scalar slot");
  if (pressure_span != nullptr &&
      pressure_span->resolved_begin < resolved_scalar.size()) {
    constexpr double kPiTarget = 17.25;
    fill_resolved_slice(resolved_scalar, *pressure_span, kPiTarget);
    OwnedView pressure{plan.pressure_field(), 1U};
    initialise_interior(pressure);
    passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                         BoundaryStage::pressure, plan,
                         Span<FieldView>{&pressure.view, 1U}, resolved)),
                     "resolved pressure stage applies");
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior = pressure.view.unchecked(
          Int3{kInterior.x - layer, 1, 1}, 0U);
      const double ghost = pressure.view.unchecked(
          Int3{kInterior.x - 1 + layer, 1, 1}, 0U);
      passed &= expect(ghost == 2.0 * kPiTarget - interior,
                       "pressure uses resolved pi, not absolute input pressure");
    }
    passed &= expect(canaries_intact(pressure),
                     "resolved pressure preserves canaries");
  }

  const BoundaryIndexSpan* thermal_span = find_span(
      plan, BoundaryStage::enthalpy, CartesianFace::y_min,
      plan.enthalpy_field());
  passed &= expect(thermal_span != nullptr &&
                       thermal_span->value_source ==
                           BoundaryValueSource::resolved_scalar,
                   "isothermal wall compiles a resolved enthalpy slot");
  if (thermal_span != nullptr &&
      thermal_span->resolved_begin < resolved_scalar.size()) {
    constexpr double kEnthalpyTarget = 612345.0;
    fill_resolved_slice(resolved_scalar, *thermal_span, kEnthalpyTarget);
    OwnedView enthalpy{plan.enthalpy_field(), 1U};
    initialise_interior(enthalpy);
    passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                         BoundaryStage::enthalpy, plan,
                         Span<FieldView>{&enthalpy.view, 1U}, resolved)),
                     "resolved enthalpy stage applies");
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior =
          enthalpy.view.unchecked(Int3{1, layer - 1, 1}, 0U);
      const double ghost =
          enthalpy.view.unchecked(Int3{1, -layer, 1}, 0U);
      passed &= expect(ghost == 2.0 * kEnthalpyTarget - interior,
                       "thermal wall uses resolved h, not temperature");
    }
    passed &= expect(canaries_intact(enthalpy),
                     "resolved enthalpy preserves canaries");
  }

  const FieldId gradient_field = scalar_field(
      plan, CartesianFace::x_min, ScalarBoundaryKind::normal_flux);
  const BoundaryIndexSpan* gradient_span = find_span(
      plan, BoundaryStage::scalar, CartesianFace::x_min, gradient_field);
  passed &= expect(gradient_span != nullptr &&
                       gradient_span->value_source ==
                           BoundaryValueSource::resolved_normal_gradient,
                   "normal-flux scalar compiles a resolved gradient slot");
  if (gradient_span != nullptr &&
      gradient_span->resolved_begin < resolved_gradient.size()) {
    constexpr double kFieldGradient = 2.5;
    fill_resolved_slice(resolved_gradient, *gradient_span, kFieldGradient);
    const std::vector<FieldId> ids = scalar_fields(plan);
    OwnedView first{ids[0], 1U};
    OwnedView second{ids[1], 1U};
    initialise_interior(first);
    initialise_interior(second);
    FieldView fields[]{first.view, second.view};
    OwnedView& scalar = ids[0] == gradient_field ? first : second;
    passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                         BoundaryStage::scalar, plan,
                         Span<FieldView>{fields, 2U}, resolved)),
                     "resolved scalar gradient stage applies");
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior =
          scalar.view.unchecked(Int3{layer - 1, 1, 1}, 0U);
      const double distance =
          layer == 1 ? plan.normal_distance_1().data[gradient_span->parameter]
                     : plan.normal_distance_2().data[gradient_span->parameter];
      const double ghost =
          scalar.view.unchecked(Int3{-layer, 1, 1}, 0U);
      passed &= expect(ghost == interior + kFieldGradient * distance,
                       "normal flux uses the resolved field-space gradient");
    }
    passed &= expect(canaries_intact(scalar),
                     "resolved scalar gradient preserves canaries");
  }

  const FieldId convective_field = scalar_field(
      plan, CartesianFace::x_max, ScalarBoundaryKind::convective);
  const BoundaryIndexSpan* convective_span = find_span(
      plan, BoundaryStage::scalar, CartesianFace::x_max, convective_field);
  passed &= expect(convective_span != nullptr &&
                       convective_span->value_source ==
                           BoundaryValueSource::resolved_scalar,
                   "convective scalar compiles a resolved field target");
  if (convective_span != nullptr &&
      convective_span->resolved_begin < resolved_scalar.size()) {
    constexpr double kConvectiveTarget = 0.875;
    fill_resolved_slice(resolved_scalar, *convective_span, kConvectiveTarget);
    const std::vector<FieldId> ids = scalar_fields(plan);
    OwnedView first{ids[0], 1U};
    OwnedView second{ids[1], 1U};
    initialise_interior(first);
    initialise_interior(second);
    FieldView fields[]{first.view, second.view};
    OwnedView& scalar = ids[0] == convective_field ? first : second;
    passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                         BoundaryStage::scalar, plan,
                         Span<FieldView>{fields, 2U}, resolved)),
                     "resolved convective scalar stage applies");
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior = scalar.view.unchecked(
          Int3{kInterior.x - layer, 1, 1}, 0U);
      const double ghost = scalar.view.unchecked(
          Int3{kInterior.x - 1 + layer, 1, 1}, 0U);
      passed &= expect(ghost == 2.0 * kConvectiveTarget - interior,
                       "convective fill uses the resolved field-space target");
    }
    passed &= expect(canaries_intact(scalar),
                     "resolved convective fill preserves canaries");
  }
  return passed;
}

bool test_resolved_values_are_face_cell_local(const BoundaryPlan& plan) {
  std::vector<double> scalar;
  std::vector<Real3> vector;
  std::vector<double> gradient;
  BoundaryResolvedValues values = resolved_values(plan, scalar, vector,
                                                  gradient);
  bool passed = true;

  const BoundaryIndexSpan* pressure = find_span(
      plan, BoundaryStage::pressure, CartesianFace::x_max,
      plan.pressure_field());
  const BoundaryIndexSpan* momentum = find_span(
      plan, BoundaryStage::momentum, CartesianFace::z_min,
      plan.velocity_field());
  const FieldId flux_field = scalar_field(
      plan, CartesianFace::x_min, ScalarBoundaryKind::normal_flux);
  const BoundaryIndexSpan* flux = find_span(
      plan, BoundaryStage::scalar, CartesianFace::x_min, flux_field);
  passed &= expect(pressure != nullptr && pressure->resolved_stride >= 2U,
                   "pressure resolved slice covers each tangential face cell");
  passed &= expect(momentum != nullptr && momentum->resolved_stride >= 2U,
                   "vector resolved slice covers each tangential face cell");
  passed &= expect(flux != nullptr && flux->resolved_stride >= 2U,
                   "gradient resolved slice covers each tangential face cell");
  if (pressure == nullptr || momentum == nullptr || flux == nullptr ||
      pressure->resolved_stride < 2U || momentum->resolved_stride < 2U ||
      flux->resolved_stride < 2U) {
    return false;
  }

  const std::size_t pressure0 = resolved_index(*pressure, 0U, 0U);
  const std::size_t pressure1 = resolved_index(*pressure, 1U, 0U);
  scalar[pressure0] = 11.0;
  scalar[pressure1] = 29.0;
  OwnedView pi{plan.pressure_field(), 1U};
  initialise_interior(pi);
  passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                       BoundaryStage::pressure, plan,
                       Span<FieldView>{&pi.view, 1U}, values)),
                   "face-local pressure values apply");
  const double pi_interior0 =
      pi.view.unchecked(Int3{kInterior.x - 1, 0, 0}, 0U);
  const double pi_interior1 =
      pi.view.unchecked(Int3{kInterior.x - 1, 1, 0}, 0U);
  passed &= expect(pi.view.unchecked(Int3{kInterior.x, 0, 0}, 0U) ==
                       2.0 * scalar[pressure0] - pi_interior0,
                   "first pressure face cell uses its own resolved value");
  passed &= expect(pi.view.unchecked(Int3{kInterior.x, 1, 0}, 0U) ==
                       2.0 * scalar[pressure1] - pi_interior1,
                   "second pressure face cell uses its own resolved value");

  const std::size_t vector0 = resolved_index(*momentum, 0U, 0U);
  const std::size_t vector1 = resolved_index(*momentum, 1U, 0U);
  vector[vector0] = Real3{1.0, 2.0, 3.0};
  vector[vector1] = Real3{4.0, 5.0, 6.0};
  OwnedView velocity{plan.velocity_field(), 3U};
  initialise_interior(velocity);
  passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                       BoundaryStage::momentum, plan,
                       Span<FieldView>{&velocity.view, 1U}, values)),
                   "face-local vector values apply");
  passed &= expect(velocity.view.unchecked(Int3{0, 0, -1}, 2U) ==
                       2.0 * vector[vector0].z -
                           velocity.view.unchecked(Int3{0, 0, 0}, 2U),
                   "first vector face cell uses its own resolved value");
  passed &= expect(velocity.view.unchecked(Int3{1, 0, -1}, 2U) ==
                       2.0 * vector[vector1].z -
                           velocity.view.unchecked(Int3{1, 0, 0}, 2U),
                   "second vector face cell uses its own resolved value");

  const std::size_t gradient0 = resolved_index(*flux, 0U, 0U);
  const std::size_t gradient1 = resolved_index(*flux, 1U, 0U);
  gradient[gradient0] = 2.0;
  gradient[gradient1] = 7.0;
  const std::vector<FieldId> ids = scalar_fields(plan);
  OwnedView first{ids[0], 1U};
  OwnedView second{ids[1], 1U};
  initialise_interior(first);
  initialise_interior(second);
  FieldView fields[]{first.view, second.view};
  OwnedView& transported = ids[0] == flux_field ? first : second;
  passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                       BoundaryStage::scalar, plan,
                       Span<FieldView>{fields, 2U}, values)),
                   "face-local normal gradients apply");
  const double distance = plan.normal_distance_1().data[flux->parameter];
  passed &= expect(transported.view.unchecked(Int3{-1, 0, 0}, 0U) ==
                       transported.view.unchecked(Int3{0, 0, 0}, 0U) +
                           gradient[gradient0] * distance,
                   "first flux face cell uses its own resolved gradient");
  passed &= expect(transported.view.unchecked(Int3{-1, 1, 0}, 0U) ==
                       transported.view.unchecked(Int3{0, 1, 0}, 0U) +
                           gradient[gradient1] * distance,
                   "second flux face cell uses its own resolved gradient");
  return passed;
}

bool test_stage_local_finite_preflight_is_atomic(const BoundaryPlan& plan) {
  std::vector<double> scalar;
  std::vector<Real3> vector;
  std::vector<double> gradient;
  BoundaryResolvedValues values = resolved_values(plan, scalar, vector,
                                                  gradient);
  const BoundaryIndexSpan* pressure = find_span(
      plan, BoundaryStage::pressure, CartesianFace::x_max,
      plan.pressure_field());
  const BoundaryIndexSpan* momentum = find_span(
      plan, BoundaryStage::momentum, CartesianFace::z_min,
      plan.velocity_field());
  const FieldId flux_field = scalar_field(
      plan, CartesianFace::x_min, ScalarBoundaryKind::normal_flux);
  const BoundaryIndexSpan* flux = find_span(
      plan, BoundaryStage::scalar, CartesianFace::x_min, flux_field);
  bool passed = expect(pressure != nullptr && momentum != nullptr &&
                           flux != nullptr,
                       "finite-preflight fixture has all resolved sources");
  if (pressure == nullptr || momentum == nullptr || flux == nullptr) {
    return false;
  }

  std::fill(scalar.begin(), scalar.end(),
            std::numeric_limits<double>::quiet_NaN());
  std::fill(vector.begin(), vector.end(),
            Real3{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0});
  std::fill(gradient.begin(), gradient.end(),
            std::numeric_limits<double>::infinity());
  fill_resolved_slice(scalar, *pressure, 13.0);
  OwnedView pi{plan.pressure_field(), 1U};
  initialise_interior(pi);
  passed &= expect(static_cast<bool>(apply_boundary_ghosts(
                       BoundaryStage::pressure, plan,
                       Span<FieldView>{&pi.view, 1U}, values)),
                   "unreferenced-stage NaN and Inf do not block pressure");

  scalar[pressure->resolved_begin] =
      std::numeric_limits<double>::quiet_NaN();
  OwnedView rejected_pi{plan.pressure_field(), 1U};
  initialise_interior(rejected_pi);
  const std::vector<double> before_pi = rejected_pi.storage;
  passed &= expect(apply_boundary_ghosts(
                       BoundaryStage::pressure, plan,
                       Span<FieldView>{&rejected_pi.view, 1U}, values)
                           .code == StatusCode::invalid_plan,
                   "current-stage NaN scalar is rejected");
  passed &= expect(unchanged(rejected_pi, before_pi),
                   "NaN scalar rejection is stage atomic");

  fill_resolved_slice(vector, *momentum, Real3{1.0, 2.0, 3.0});
  vector[momentum->resolved_begin].y =
      std::numeric_limits<double>::infinity();
  OwnedView rejected_velocity{plan.velocity_field(), 3U};
  initialise_interior(rejected_velocity);
  const std::vector<double> before_velocity = rejected_velocity.storage;
  passed &= expect(apply_boundary_ghosts(
                       BoundaryStage::momentum, plan,
                       Span<FieldView>{&rejected_velocity.view, 1U}, values)
                           .code == StatusCode::invalid_plan,
                   "current-stage Inf vector is rejected");
  passed &= expect(unchanged(rejected_velocity, before_velocity),
                   "Inf vector rejection is stage atomic");

  std::fill(scalar.begin(), scalar.end(), 0.0);
  std::fill(gradient.begin(), gradient.end(), 0.0);
  gradient[flux->resolved_begin] =
      std::numeric_limits<double>::quiet_NaN();
  const std::vector<FieldId> ids = scalar_fields(plan);
  OwnedView first{ids[0], 1U};
  OwnedView second{ids[1], 1U};
  initialise_interior(first);
  initialise_interior(second);
  const std::vector<double> before_first = first.storage;
  const std::vector<double> before_second = second.storage;
  FieldView fields[]{first.view, second.view};
  passed &= expect(apply_boundary_ghosts(
                       BoundaryStage::scalar, plan,
                       Span<FieldView>{fields, 2U}, values)
                           .code == StatusCode::invalid_plan,
                   "current-stage NaN gradient is rejected");
  passed &= expect(unchanged(first, before_first) &&
                       unchanged(second, before_second),
                   "NaN gradient rejects the complete scalar stage");
  return passed;
}

bool test_two_layer_relations_and_zero_allocations(const BoundaryPlan& plan) {
  OwnedView velocity{plan.velocity_field(), 3U};
  initialise_interior(velocity);
  std::vector<double> resolved_scalar;
  std::vector<Real3> resolved_vector;
  std::vector<double> resolved_gradient;
  const BoundaryResolvedValues resolved = resolved_values(
      plan, resolved_scalar, resolved_vector, resolved_gradient);
  bool passed = expect(static_cast<bool>(apply_boundary_ghosts(
                           BoundaryStage::momentum, plan,
                           Span<FieldView>{&velocity.view, 1U}, resolved)),
                       "valid non-empty momentum stage applies");
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    const double target = component == 0U ? 2.0 : (component == 1U ? 3.0 : 4.0);
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior =
          velocity.view.unchecked(Int3{layer - 1, 1, 1}, component);
      const double ghost =
          velocity.view.unchecked(Int3{-layer, 1, 1}, component);
      passed &= expect(ghost == 2.0 * target - interior,
                       "Dirichlet ghost mirrors the matching interior layer");
    }
  }
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    for (std::int32_t layer = 1; layer <= 2; ++layer) {
      const double interior = velocity.view.unchecked(
          Int3{kInterior.x - layer, 1, 1}, component);
      const double ghost = velocity.view.unchecked(
          Int3{kInterior.x - 1 + layer, 1, 1}, component);
      passed &= expect(ghost == interior,
                       "zero-gradient ghost copies the matching interior layer");
    }
  }
  passed &= expect(canaries_intact(velocity),
                   "valid two-layer application preserves canaries");

  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard observe;
    for (std::size_t iteration = 0U; iteration < 128U; ++iteration) {
      const Status status = apply_boundary_ghosts(
          BoundaryStage::momentum, plan,
          Span<FieldView>{&velocity.view, 1U}, resolved);
      if (!status) {
        passed = false;
        break;
      }
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(hot_allocations == 0U,
                   "preflight plus repeated apply performs zero C++ allocations");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  BoundaryPlan plan;
  bool passed = test_empty_plan();
  passed &= expect(compile_plan(plan), "non-empty boundary fixture compiles");
  if (plan.semantic_fingerprint() != 0U) {
    passed &= test_stage_atomic_when_second_field_is_missing(plan);
    passed &= test_view_preflight_is_atomic(plan);
    passed &= test_resolved_shape_preflight_is_atomic(plan);
    passed &= test_resolved_field_space_relations(plan);
    passed &= test_resolved_values_are_face_cell_local(plan);
    passed &= test_stage_local_finite_preflight_is_atomic(plan);
    passed &= test_two_layer_relations_and_zero_allocations(plan);
  }
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 boundary apply tests passed\n";
  return 0;
}
