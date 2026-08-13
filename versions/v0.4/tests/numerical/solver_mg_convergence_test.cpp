// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <algorithm>
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
#include <vector>

namespace allocation_observer {
std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};
void* allocate(std::size_t bytes) {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
  void* const result = std::malloc(bytes == 0U ? 1U : bytes);
  if (result == nullptr) {
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
};
}  // namespace allocation_observer

void* operator new(std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new[](std::size_t bytes) {
  return allocation_observer::allocate(bytes);
}
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  void* result = nullptr;
  if (allocation_observer::enabled.load(std::memory_order_relaxed)) {
    allocation_observer::count.fetch_add(1U, std::memory_order_relaxed);
  }
  if (posix_memalign(&result, static_cast<std::size_t>(alignment),
                     bytes == 0U ? static_cast<std::size_t>(alignment)
                                 : bytes) != 0 ||
      result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}
void* operator new[](std::size_t bytes, std::align_val_t alignment) {
  return ::operator new(bytes, alignment);
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
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
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

constexpr double kPi = 3.141592653589793238462643383279502884;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

CartesianMeshSpec uniform_mesh(Int3 cells) {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::uniform;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {1.0, 1.0, 1.0};
  mesh.has_exact_cells = true;
  mesh.exact_cells = cells;
  mesh.minimum_spacing = {1.0e-8, 1.0e-8, 1.0e-8};
  mesh.max_growth_ratio = 1.0;
  mesh.limits = {10000000U, 1U << 30U};
  return mesh;
}

CartesianMeshSpec stretched_variable_mesh() {
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::tensor_stretched;
  mesh.lower = {0.0, 0.0, 0.0};
  mesh.upper = {4.0, 4.0, 1.0};
  mesh.has_base_spacing = true;
  mesh.base_spacing = {0.25, 0.25, 0.125};
  mesh.minimum_spacing = {0.025, 0.025, 0.125};
  mesh.max_growth_ratio = 1.25;
  mesh.focus_regions.push_back(
      {{1.0, 1.0, 0.0}, {3.0, 3.0, 1.0}, {0.05, 0.05, 0.125}});
  mesh.limits = {10000000U, 1U << 30U};
  return mesh;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField field(Int3 cells, std::uint8_t components, std::uint8_t ghosts,
                 FieldId id, StorageIdentity storage) {
  OwnedField result;
  const std::size_t nx = static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t ny = static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t nz = static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.storage.assign(nx * ny * nz * components, 0.0);
  result.view.base = result.storage.data() + ghosts + ghosts * nx +
                     ghosts * nx * ny;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.component_stride = nx * ny * nz;
  result.view.field = id;
  result.view.revision = 1U;
  result.view.storage_identity = storage;
  result.view.revision_domain = 901U;
  return result;
}

struct OwnedFace {
  std::vector<double> storage;
  FaceFieldView view{};
};

OwnedFace face(Int3 extents, CartesianAxis axis, StorageIdentity storage) {
  OwnedFace result;
  const auto nx = static_cast<std::size_t>(extents.x);
  const auto ny = static_cast<std::size_t>(extents.y);
  const auto nz = static_cast<std::size_t>(extents.z);
  result.storage.assign(nx * ny * nz, 0.0);
  result.view.base = result.storage.data();
  result.view.extents = extents;
  result.view.stride_y = nx;
  result.view.stride_z = nx * ny;
  result.view.axis = axis;
  result.view.storage_identity = storage;
  result.view.revision_domain = 902U;
  return result;
}

ConstFaceFieldView constant(FaceFieldView view) noexcept {
  return {view.base,          view.extents, view.stride_y,
          view.stride_z,     view.axis,    view.storage_identity,
          view.revision_domain};
}

MgBoundaryKind minimum_boundary(const MgBoundarySet& boundaries,
                                CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? boundaries.x_min
             : (axis == CartesianAxis::y ? boundaries.y_min
                                         : boundaries.z_min);
}

MgBoundaryKind maximum_boundary(const MgBoundarySet& boundaries,
                                CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? boundaries.x_max
             : (axis == CartesianAxis::y ? boundaries.y_max
                                         : boundaries.z_max);
}

double coordinate(const AxisMetrics& axis, std::int32_t index) noexcept {
  return axis.centres().data[static_cast<std::size_t>(index)];
}

double width(const AxisMetrics& axis, std::int32_t index) noexcept {
  return axis.widths().data[static_cast<std::size_t>(index)];
}

double gamma(double x, double y, double z, bool variable) noexcept {
  (void)z;
  // The variable case is periodic in z, so its coefficient must also be
  // periodic there.  x/y variation still exercises metric restriction and
  // variable coarse coefficients without making the oracle nonsymmetric.
  return variable ? 1.0 + 0.30 * x + 0.20 * y : 1.0;
}

struct NumericalFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch{};
  MgBoundarySet boundaries{};
  bool variable{};
  OwnedField diagonal;
  OwnedFace x_faces;
  OwnedFace y_faces;
  OwnedFace z_faces;
  OwnedField vectors;
  OwnedField exact;
  OwnedField rhs;
  OwnedField solution;
  OwnedField residual;
  OwnedField correction;
  OwnedField applied;
  MgWorkspaceRequirements workspace_requirements{};
  MgWorkspace workspace;
  ReductionEngine reductions;
  HaloEngine halo;
  std::vector<HaloEngine> coarse_halos;
  std::vector<HaloEngine*> coarse_halo_pointers;
  NativeCartesianMgPlan plan;

  bool create(const CartesianMeshSpec& mesh, MgBoundarySet selected,
              MgNullSpace null_space, bool variable_coefficient) {
    boundaries = selected;
    variable = variable_coefficient;
    if (!CartesianGeometryCompiler::compile(
            MPI_COMM_SELF, mesh, GeometryBudget{}, geometry, patch)) {
      return false;
    }
    const Int3 cells = patch.cells;
    diagonal = field(cells, 1U, 0U, 1U, 700U);
    x_faces = face({cells.x + 1, cells.y, cells.z}, CartesianAxis::x,
                   710U);
    y_faces = face({cells.x, cells.y + 1, cells.z}, CartesianAxis::y,
                   711U);
    z_faces = face({cells.x, cells.y, cells.z + 1}, CartesianAxis::z,
                   712U);
    exact = field(cells, 1U, 0U, 20U, 720U);
    rhs = field(cells, 1U, 0U, 21U, 721U);
    solution = field(cells, 1U, 0U, 22U, 722U);
    residual = field(cells, 1U, 0U, 23U, 723U);
    correction = field(cells, 1U, 0U, 24U, 724U);
    applied = field(cells, 1U, 0U, 25U, 725U);
    fill_coefficients();

    NativeCartesianMgSpec spec;
    spec.communicator = MPI_COMM_SELF;
    spec.geometry = &geometry;
    spec.patch = patch;
    spec.boundaries = boundaries;
    spec.null_space = null_space;
    spec.policy.anisotropy_threshold = 4.0;
    spec.policy.coefficient_change_rebuild_ratio = 0.25;
    spec.policy.pre_sweeps = 2U;
    spec.policy.post_sweeps = 2U;
    spec.policy.maximum_levels = 20U;
    spec.policy.coarse_sweeps = 32U;
    spec.policy.minimum_coarse_extent = 3U;
    spec.policy.line_relaxation_maximum_extent = 4096U;
    spec.identity = {301U, 302U, 303U, 304U, 305U};
    spec.coefficients = {1U, variable ? 402U : 401U, 0.0};
    if (!make_mg_workspace_requirements(MPI_COMM_SELF, geometry, patch,
                                        spec.policy, 81U,
                                        workspace_requirements)) {
      return false;
    }
    vectors = field(workspace_requirements.arena_shape, 1U, 0U, 10U, 500U);
    if (!MgWorkspace::bind(workspace_requirements, vectors.view, workspace) ||
        !ReductionEngine::compile(MPI_COMM_SELF,
                                  ReductionMode::mpi_allreduce, 4U,
                                  reductions)) {
      return false;
    }
    const std::array<HaloFieldSpec, 1U> halo_fields{{{10U, 1U, 1U}}};
    const HaloTopology topology{
        boundaries.x_min == MgBoundaryKind::periodic,
        boundaries.y_min == MgBoundaryKind::periodic,
        boundaries.z_min == MgBoundaryKind::periodic};
    if (!halo.reserve(MPI_COMM_SELF,
                      workspace_requirements.levels[0U].patch,
                      {halo_fields.data(), halo_fields.size()}, topology)) {
      return false;
    }
    coarse_halos.resize(workspace_requirements.level_count - 1U);
    coarse_halo_pointers.resize(coarse_halos.size());
    for (std::size_t level = 1U;
         level < workspace_requirements.level_count; ++level) {
      if (!coarse_halos[level - 1U].reserve(
              MPI_COMM_SELF, workspace_requirements.levels[level].patch,
              {halo_fields.data(), halo_fields.size()}, topology)) {
        return false;
      }
      coarse_halo_pointers[level - 1U] = &coarse_halos[level - 1U];
    }
    const MgRuntimeServices services{
        &halo, &reductions, &workspace,
        {coarse_halo_pointers.data(), coarse_halo_pointers.size()}};
    return static_cast<bool>(NativeCartesianMgPlan::compile(
        spec, services, coefficient_views(), plan));
  }

  MgCoefficientViews coefficient_views() const noexcept {
    return {as_const(diagonal.view), constant(x_faces.view),
            constant(y_faces.view), constant(z_faces.view)};
  }

  double face_coefficient(CartesianAxis axis, Int3 face_index) const noexcept {
    const AxisMetrics& mx = geometry.x();
    const AxisMetrics& my = geometry.y();
    const AxisMetrics& mz = geometry.z();
    const auto& metric = geometry.axis(axis);
    const std::int32_t normal = axis == CartesianAxis::x
                                    ? face_index.x
                                    : (axis == CartesianAxis::y
                                           ? face_index.y
                                           : face_index.z);
    const std::int32_t count = axis == CartesianAxis::x
                                   ? patch.cells.x
                                   : (axis == CartesianAxis::y ? patch.cells.y
                                                               : patch.cells.z);
    const MgBoundaryKind boundary =
        normal == 0 ? minimum_boundary(boundaries, axis)
                    : maximum_boundary(boundaries, axis);
    if ((normal == 0 || normal == count) &&
        boundary == MgBoundaryKind::neumann) {
      return 0.0;
    }
    const double x = axis == CartesianAxis::x
                         ? metric.faces().data[static_cast<std::size_t>(normal)]
                         : coordinate(mx, face_index.x);
    const double y = axis == CartesianAxis::y
                         ? metric.faces().data[static_cast<std::size_t>(normal)]
                         : coordinate(my, face_index.y);
    const double z = axis == CartesianAxis::z
                         ? metric.faces().data[static_cast<std::size_t>(normal)]
                         : coordinate(mz, face_index.z);
    const double area = axis == CartesianAxis::x
                            ? width(my, face_index.y) * width(mz, face_index.z)
                            : (axis == CartesianAxis::y
                                   ? width(mx, face_index.x) *
                                         width(mz, face_index.z)
                                   : width(mx, face_index.x) *
                                         width(my, face_index.y));
    double distance = 0.0;
    if (normal > 0 && normal < count) {
      distance = metric.centres().data[static_cast<std::size_t>(normal)] -
                 metric.centres().data[static_cast<std::size_t>(normal - 1)];
    } else if (boundary == MgBoundaryKind::periodic) {
      distance = 0.5 * (metric.widths().data[0U] +
                        metric.widths().data[metric.widths().size - 1U]);
    } else {
      distance = 0.5 * metric.widths().data[
                           normal == 0 ? 0U : metric.widths().size - 1U];
    }
    return gamma(x, y, z, variable) * area / distance;
  }

  void fill_coefficients() noexcept {
    const Int3 cells = patch.cells;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i <= cells.x; ++i) {
          x_faces.view.unchecked({i, j, k}) =
              face_coefficient(CartesianAxis::x, {i, j, k});
        }
      }
    }
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j <= cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          y_faces.view.unchecked({i, j, k}) =
              face_coefficient(CartesianAxis::y, {i, j, k});
        }
      }
    }
    for (std::int32_t k = 0; k <= cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          z_faces.view.unchecked({i, j, k}) =
              face_coefficient(CartesianAxis::z, {i, j, k});
        }
      }
    }
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          diagonal.view.unchecked({i, j, k}, 0U) =
              x_faces.view.unchecked({i, j, k}) +
              x_faces.view.unchecked({i + 1, j, k}) +
              y_faces.view.unchecked({i, j, k}) +
              y_faces.view.unchecked({i, j + 1, k}) +
              z_faces.view.unchecked({i, j, k}) +
              z_faces.view.unchecked({i, j, k + 1});
        }
      }
    }
  }

  double neighbor(ConstFieldView input, Int3 cell, CartesianAxis axis,
                  int direction) const noexcept {
    Int3 selected = cell;
    std::int32_t* coordinate_index = axis == CartesianAxis::x
                                         ? &selected.x
                                         : (axis == CartesianAxis::y
                                                ? &selected.y
                                                : &selected.z);
    const std::int32_t count = axis == CartesianAxis::x
                                   ? patch.cells.x
                                   : (axis == CartesianAxis::y ? patch.cells.y
                                                               : patch.cells.z);
    *coordinate_index += direction;
    if (*coordinate_index >= 0 && *coordinate_index < count) {
      return input.unchecked(selected, 0U);
    }
    const MgBoundaryKind boundary =
        direction < 0 ? minimum_boundary(boundaries, axis)
                      : maximum_boundary(boundaries, axis);
    if (boundary == MgBoundaryKind::periodic) {
      *coordinate_index = direction < 0 ? count - 1 : 0;
      return input.unchecked(selected, 0U);
    }
    return 0.0;
  }

  void apply_operator(ConstFieldView input, FieldView output) const noexcept {
    const Int3 cells = patch.cells;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 cell{i, j, k};
          double value = diagonal.view.unchecked(cell, 0U) *
                         input.unchecked(cell, 0U);
          value -= x_faces.view.unchecked({i, j, k}) *
                   neighbor(input, cell, CartesianAxis::x, -1);
          value -= x_faces.view.unchecked({i + 1, j, k}) *
                   neighbor(input, cell, CartesianAxis::x, 1);
          value -= y_faces.view.unchecked({i, j, k}) *
                   neighbor(input, cell, CartesianAxis::y, -1);
          value -= y_faces.view.unchecked({i, j + 1, k}) *
                   neighbor(input, cell, CartesianAxis::y, 1);
          value -= z_faces.view.unchecked({i, j, k}) *
                   neighbor(input, cell, CartesianAxis::z, -1);
          value -= z_faces.view.unchecked({i, j, k + 1}) *
                   neighbor(input, cell, CartesianAxis::z, 1);
          output.unchecked(cell, 0U) = value;
        }
      }
    }
  }

  void make_manufactured_solution(bool null_space) noexcept {
    const Int3 cells = patch.cells;
    double sum = 0.0;
    std::size_t count = 0U;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const double x = coordinate(geometry.x(), i);
          const double y = coordinate(geometry.y(), j);
          const double z = coordinate(geometry.z(), k);
          const double value = null_space
                                   ? std::cos(2.0 * kPi * x) *
                                         std::cos(2.0 * kPi * y) *
                                         std::cos(2.0 * kPi * z)
                                   : std::sin(kPi * x) *
                                         std::cos(kPi * y) *
                                         std::sin(2.0 * kPi * z);
          exact.view.unchecked({i, j, k}, 0U) = value;
          sum += value;
          ++count;
        }
      }
    }
    if (null_space) {
      const double mean = sum / static_cast<double>(count);
      for (std::int32_t k = 0; k < cells.z; ++k) {
        for (std::int32_t j = 0; j < cells.y; ++j) {
          for (std::int32_t i = 0; i < cells.x; ++i) {
            exact.view.unchecked({i, j, k}, 0U) -= mean;
          }
        }
      }
    }
    apply_operator(as_const(exact.view), rhs.view);
  }

  double residual_norm() noexcept {
    apply_operator(as_const(solution.view), applied.view);
    double sum = 0.0;
    const Int3 cells = patch.cells;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const double value = rhs.view.unchecked({i, j, k}, 0U) -
                               applied.view.unchecked({i, j, k}, 0U);
          residual.view.unchecked({i, j, k}, 0U) = value;
          sum += value * value;
        }
      }
    }
    return std::sqrt(sum);
  }

  bool correction_is_finite() const noexcept {
    const Int3 cells = patch.cells;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          if (!std::isfinite(correction.view.unchecked({i, j, k}, 0U))) {
            return false;
          }
        }
      }
    }
    return true;
  }

  double correction_mean() const noexcept {
    double sum = 0.0;
    std::size_t count = 0U;
    const Int3 cells = patch.cells;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          sum += correction.view.unchecked({i, j, k}, 0U);
          ++count;
        }
      }
    }
    return sum / static_cast<double>(count);
  }

  bool run_cycles(std::uint32_t cycles, double maximum_ratio,
                  bool require_projection) noexcept {
    double previous = residual_norm();
    const double initial = previous;
    bool passed = std::isfinite(initial) && initial > 0.0;
    {
      allocation_observer::Guard guard;
      for (std::uint32_t cycle = 0U; cycle < cycles; ++cycle) {
        std::fill(correction.storage.begin(), correction.storage.end(), 91.0);
        const Status applied_status = plan.apply(
            as_const(residual.view), correction.view, cycle);
        passed = static_cast<bool>(applied_status) &&
                 correction_is_finite() && passed;
        if (require_projection) {
          passed = std::abs(correction_mean()) <= 1.0e-12 && passed;
        }
        const Int3 cells = patch.cells;
        for (std::int32_t k = 0; k < cells.z; ++k) {
          for (std::int32_t j = 0; j < cells.y; ++j) {
            for (std::int32_t i = 0; i < cells.x; ++i) {
              solution.view.unchecked({i, j, k}, 0U) +=
                  correction.view.unchecked({i, j, k}, 0U);
            }
          }
        }
        const double current = residual_norm();
        passed = std::isfinite(current) &&
                 current <= previous * (1.0 + 1.0e-12) && passed;
        previous = current;
      }
    }
    return passed && previous <= maximum_ratio * initial &&
           allocation_observer::count.load(std::memory_order_relaxed) == 0U;
  }
};

bool test_uniform_mixed_boundaries() {
  const MgBoundarySet mixed{
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet,
      MgBoundaryKind::neumann,   MgBoundaryKind::neumann,
      MgBoundaryKind::periodic,  MgBoundaryKind::periodic};
  NumericalFixture fixture;
  bool passed = expect(
      fixture.create(uniform_mesh({16, 12, 8}), mixed, MgNullSpace::none,
                     false),
      "uniform mixed-boundary hierarchy compiles");
  if (!passed) {
    return false;
  }
  fixture.make_manufactured_solution(false);
  passed &= expect(fixture.run_cycles(10U, 0.10, false),
                   "uniform mixed-boundary V-cycles decrease true residual monotonically within a fixed bound and allocate nothing");
  return passed;
}

bool test_stretched_variable_coefficients() {
  const MgBoundarySet mixed{
      MgBoundaryKind::dirichlet, MgBoundaryKind::dirichlet,
      MgBoundaryKind::neumann,   MgBoundaryKind::neumann,
      MgBoundaryKind::periodic,  MgBoundaryKind::periodic};
  NumericalFixture fixture;
  bool passed = expect(
      fixture.create(stretched_variable_mesh(), mixed, MgNullSpace::none,
                     true),
      "tensor-stretched variable-coefficient hierarchy compiles");
  if (!passed) {
    return false;
  }
  fixture.make_manufactured_solution(false);
  passed &= expect(fixture.run_cycles(12U, 0.20, false),
                   "stretched variable-coefficient V-cycles remain finite, monotone, and bounded");
  return passed;
}

bool test_constant_null_space_projection() {
  const MgBoundarySet singular{
      MgBoundaryKind::neumann,  MgBoundaryKind::neumann,
      MgBoundaryKind::neumann,  MgBoundaryKind::neumann,
      MgBoundaryKind::periodic, MgBoundaryKind::periodic};
  NumericalFixture fixture;
  bool passed = expect(
      fixture.create(uniform_mesh({16, 12, 8}), singular,
                     MgNullSpace::constant, false),
      "singular Neumann-periodic hierarchy compiles with constant null space");
  if (!passed) {
    return false;
  }
  fixture.make_manufactured_solution(true);
  passed &= expect(fixture.run_cycles(10U, 0.15, true),
                   "singular V-cycles project the constant mode and reduce the independent true residual");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  bool passed = true;
  passed &= test_uniform_mixed_boundaries();
  passed &= test_stretched_variable_coefficients();
  passed &= test_constant_null_space_projection();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 native Cartesian MG numerical tests passed\n";
  return 0;
}
