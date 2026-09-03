// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_boundary.hpp"
#include "hundun/v04_execution.hpp"

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
#include <memory>
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

void* allocate(std::size_t size) {
  observe();
  void* const result = std::malloc(size == 0U ? 1U : size);
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

constexpr Int3 kCells{6, 5, 4};

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

std::size_t cell_count(Int3 cells) {
  return static_cast<std::size_t>(cells.x) *
         static_cast<std::size_t>(cells.y) *
         static_cast<std::size_t>(cells.z);
}

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool close(double actual, double expected, double scale = 1.0) {
  const double tolerance =
      128.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(scale), std::abs(expected)});
  return std::isfinite(actual) && std::abs(actual - expected) <= tolerance;
}

bool direction_close(double actual, double expected) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <= 5.0e-8 * scale;
}

CartesianMeshSpec mesh_spec() {
  CartesianMeshSpec value;
  value.kind = GeometryKind::uniform;
  value.lower = {0.0, 0.0, 0.0};
  value.upper = {static_cast<double>(kCells.x),
                 static_cast<double>(kCells.y),
                 static_cast<double>(kCells.z)};
  value.has_exact_cells = true;
  value.exact_cells = kCells;
  value.minimum_spacing = {1.0, 1.0, 1.0};
  value.max_growth_ratio = 1.0;
  value.limits.max_global_cells = cell_count(kCells);
  value.limits.max_memory_bytes_per_rank = 1U << 24U;
  return value;
}

ValidatedModel periodic_model(double limiter = 1.0) {
  ValidatedModel value;
  value.fingerprint = 0x92a4c731U;
  value.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& face : value.boundaries) {
    face.flow_kind = BoundaryKind::periodic;
    face.thermal_kind = BoundaryKind::none;
    face.mach_limit = 0.95;
  }
  value.schemes.momentum = ConvectionScheme::central2;
  value.schemes.enthalpy = ConvectionScheme::central2;
  value.schemes.species = ConvectionScheme::central2;
  value.schemes.passive_scalar = ConvectionScheme::central2;
  value.schemes.diffusion = DiffusionScheme::central2;
  value.schemes.limiter = limiter;
  value.time = TimeControlSpec{};
  return value;
}

struct KernelFixture {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  BoundaryPlan boundary;
  SchemePlan schemes;
  CartesianKernelPlan kernels;
};

bool make_kernel_fixture(KernelFixture& fixture,
                         const CartesianMeshSpec& mesh,
                         double limiter) {
  FieldRegistry registry;
  TimeSchemePlan time;
  const Status geometry_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh, GeometryBudget{}, fixture.geometry,
      fixture.patch);
  if (!geometry_status) {
    std::cerr << "geometry compile detail=" << geometry_status.detail << '\n';
  }
  bool passed = expect(static_cast<bool>(geometry_status),
                       "Cartesian geometry compiles");
  passed &= expect(
      static_cast<bool>(BoundaryCompiler::compile(
          MPI_COMM_SELF, periodic_model(limiter), fixture.geometry,
          fixture.patch, registry, fixture.boundary, fixture.schemes, time)),
      "periodic boundary and central schemes compile");
  passed &= expect(
      static_cast<bool>(CartesianKernelPlan::compile(
          fixture.schemes, fixture.geometry, fixture.patch, fixture.boundary,
          fixture.kernels)),
      "Cartesian kernel plan compiles once outside the hot path");
  return passed;
}

bool make_kernel_fixture(KernelFixture& fixture) {
  return make_kernel_fixture(fixture, mesh_spec(), 1.0);
}

CartesianMeshSpec stretched_mesh_spec() {
  CartesianMeshSpec value;
  value.kind = GeometryKind::tensor_stretched;
  value.lower = {0.0, 0.0, 0.0};
  value.upper = {1.0, 1.0, 1.0};
  value.has_exact_cells = true;
  value.exact_cells = kCells;
  value.has_base_spacing = true;
  value.base_spacing = {1.08 / static_cast<double>(kCells.x),
                        1.08 / static_cast<double>(kCells.y),
                        1.08 / static_cast<double>(kCells.z)};
  value.minimum_spacing = {0.92 / static_cast<double>(kCells.x),
                           0.92 / static_cast<double>(kCells.y),
                           0.92 / static_cast<double>(kCells.z)};
  value.max_growth_ratio = 1.2;
  value.focus_regions.push_back(
      {{0.3, 0.3, 0.3}, {0.7, 0.7, 0.7},
       {1.0 / static_cast<double>(kCells.x),
        1.0 / static_cast<double>(kCells.y),
        1.0 / static_cast<double>(kCells.z)}});
  value.limits.max_global_cells = cell_count(kCells);
  value.limits.max_memory_bytes_per_rank = 1U << 24U;
  return value;
}

struct OwnedField {
  std::vector<double> allocation;
  FieldView view{};
};

OwnedField make_field(FieldId field, Int3 cells, std::uint8_t components,
                      std::uint8_t ghosts, RevisionToken revision) {
  OwnedField result;
  const std::size_t width =
      static_cast<std::size_t>(cells.x + 2 * ghosts);
  const std::size_t height =
      static_cast<std::size_t>(cells.y + 2 * ghosts);
  const std::size_t depth =
      static_cast<std::size_t>(cells.z + 2 * ghosts);
  result.allocation.assign(width * height * depth * components, 0.0);
  result.view.base =
      result.allocation.data() + ghosts +
      static_cast<std::size_t>(ghosts) * width +
      static_cast<std::size_t>(ghosts) * width * height;
  result.view.interior = cells;
  result.view.ghosts = {ghosts, ghosts, ghosts};
  result.view.components = components;
  result.view.stride_y = width;
  result.view.stride_z = width * height;
  result.view.component_stride = width * height * depth;
  result.view.field = field;
  result.view.revision = revision;
  result.view.storage_identity = static_cast<StorageIdentity>(field) + 1001U;
  result.view.revision_domain = 9001U;
  return result;
}

void fill_faces(FaceFluxView flux, double value) {
  for (std::int32_t z = 0; z < flux.x.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.x.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.x.extents.x; ++x) {
        flux.x.unchecked({x, y, z}) = value;
      }
    }
  }
  for (std::int32_t z = 0; z < flux.y.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.y.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.y.extents.x; ++x) {
        flux.y.unchecked({x, y, z}) = value;
      }
    }
  }
  for (std::int32_t z = 0; z < flux.z.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.z.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.z.extents.x; ++x) {
        flux.z.unchecked({x, y, z}) = value;
      }
    }
  }
}

template <class FaceView>
std::uint64_t face_checksum(FaceView flux) {
  std::uint64_t hash = 14695981039346656037ULL;
  const auto mix = [&hash](double value) {
    std::uint64_t word = bits(value);
    for (std::size_t byte = 0U; byte < sizeof(word); ++byte) {
      hash ^= (word >> (8U * byte)) & 0xffU;
      hash *= 1099511628211ULL;
    }
  };
  for (std::int32_t z = 0; z < flux.x.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.x.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.x.extents.x; ++x) {
        mix(flux.x.unchecked({x, y, z}));
      }
    }
  }
  for (std::int32_t z = 0; z < flux.y.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.y.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.y.extents.x; ++x) {
        mix(flux.y.unchecked({x, y, z}));
      }
    }
  }
  for (std::int32_t z = 0; z < flux.z.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.z.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.z.extents.x; ++x) {
        mix(flux.z.unchecked({x, y, z}));
      }
    }
  }
  return hash;
}

std::uint64_t allocation_checksum(const OwnedField& field) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (double value : field.allocation) {
    const std::uint64_t word = bits(value);
    for (std::size_t byte = 0U; byte < sizeof(word); ++byte) {
      hash ^= (word >> (8U * byte)) & 0xffU;
      hash *= 1099511628211ULL;
    }
  }
  return hash;
}

bool same_counters(const KernelCounters& left,
                   const KernelCounters& right) {
  return left.invocations == right.invocations &&
         left.cells == right.cells && left.faces == right.faces &&
         left.logical_bytes_read == right.logical_bytes_read &&
         left.logical_bytes_written == right.logical_bytes_written;
}

std::uint64_t reconstructed_face_count(KernelBox box) {
  const auto x = static_cast<std::uint64_t>(box.cells.x);
  const auto y = static_cast<std::uint64_t>(box.cells.y);
  const auto z = static_cast<std::uint64_t>(box.cells.z);
  return (x + (box.begin.x == 0 ? 1U : 0U)) * y * z +
         x * (y + (box.begin.y == 0 ? 1U : 0U)) * z +
         x * y * (z + (box.begin.z == 0 ? 1U : 0U));
}

void fill_cell_field(OwnedField& field, double offset) {
  const Int3 ghosts{static_cast<std::int32_t>(field.view.ghosts.x),
                    static_cast<std::int32_t>(field.view.ghosts.y),
                    static_cast<std::int32_t>(field.view.ghosts.z)};
  for (std::int32_t z = -ghosts.z;
       z < field.view.interior.z + ghosts.z; ++z) {
    for (std::int32_t y = -ghosts.y;
         y < field.view.interior.y + ghosts.y; ++y) {
      for (std::int32_t x = -ghosts.x;
           x < field.view.interior.x + ghosts.x; ++x) {
        for (std::uint8_t component = 0U;
             component < field.view.components; ++component) {
          field.view.unchecked({x, y, z}, component) =
              offset + 0.125 * x - 0.0625 * y + 0.03125 * z +
              0.5 * static_cast<double>(component);
        }
      }
    }
  }
}

bool test_face_layout(FaceFluxStorage& storage, FaceFluxView& first) {
  Status allocation_status;
  std::size_t observed_allocations =
      std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    allocation_status =
        FaceFluxStorage::allocate_workspace(kCells, 2U, storage);
    observed_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  bool passed = expect(static_cast<bool>(allocation_status) &&
                           observed_allocations == 1U,
                       "active/pending face storage performs one allocation");
  FaceFluxView second;
  passed &= expect(static_cast<bool>(storage.workspace_view(0U, 17U, first)) &&
                       static_cast<bool>(storage.workspace_view(1U, 18U, second)),
                   "both preallocated face replicas expose borrowed views");
  passed &= expect(first.x.extents.x == kCells.x + 1 &&
                       first.x.extents.y == kCells.y &&
                       first.x.extents.z == kCells.z &&
                       first.y.extents.x == kCells.x &&
                       first.y.extents.y == kCells.y + 1 &&
                       first.y.extents.z == kCells.z &&
                       first.z.extents.x == kCells.x &&
                       first.z.extents.y == kCells.y &&
                       first.z.extents.z == kCells.z + 1,
                   "x/y/z arrays have exact staggered face extents");
  const std::array<const double*, 6U> bases{
      first.x.base, first.y.base, first.z.base,
      second.x.base, second.y.base, second.z.base};
  for (const double* base : bases) {
    passed &= expect(reinterpret_cast<std::uintptr_t>(base) % 64U == 0U,
                     "every directional replica begins on a 64-byte boundary");
  }
  passed &= expect(first.x.stride_y >=
                           static_cast<std::size_t>(first.x.extents.x) &&
                       first.y.stride_y >=
                           static_cast<std::size_t>(first.y.extents.x) &&
                       first.z.stride_y >=
                           static_cast<std::size_t>(first.z.extents.x) &&
                       first.x.stride_y % 8U == 0U &&
                       first.y.stride_y % 8U == 0U &&
                       first.z.stride_y % 8U == 0U,
                   "face rows are padded while x remains unit stride");
  passed &= expect(&first.x.unchecked({1, 0, 0}) -
                           &first.x.unchecked({0, 0, 0}) ==
                       1 &&
                       &first.y.unchecked({1, 0, 0}) -
                               &first.y.unchecked({0, 0, 0}) ==
                           1 &&
                       &first.z.unchecked({1, 0, 0}) -
                               &first.z.unchecked({0, 0, 0}) ==
                           1,
                   "all directional arrays are x-contiguous");
  const FaceFluxStorageCounters counters = storage.counters();
  const std::uint64_t logical_faces =
      static_cast<std::uint64_t>(kCells.x + 1) * kCells.y * kCells.z +
      static_cast<std::uint64_t>(kCells.x) * (kCells.y + 1) * kCells.z +
      static_cast<std::uint64_t>(kCells.x) * kCells.y * (kCells.z + 1);
  passed &= expect(counters.aligned_payload_allocations == 1U &&
                       counters.aligned_payload_bytes >=
                           logical_faces * 2U * sizeof(double) &&
                       counters.replicas == 2U &&
                       counters.directional_blocks == 6U,
                   "one allocation owns both replicas and six directional blocks");
  passed &= expect(first.x.storage_identity == first.y.storage_identity &&
                       first.x.storage_identity == first.z.storage_identity &&
                       first.x.storage_identity == second.x.storage_identity &&
                       first.x.revision_domain == first.y.revision_domain &&
                       first.x.revision_domain == first.z.revision_domain,
                   "one storage/revision domain identifies the entire face arena");
  return passed;
}

bool run_divergence(const CartesianKernelPlan& plan, ConstFaceFluxView flux,
                    FieldView output, KernelCounters* counters = nullptr) {
  const std::array<FieldView, 1U> writes{output};
  const KernelInvocation invocation{
      {}, Span<const FieldView>{writes.data(), writes.size()},
      {{0, 0, 0}, output.interior}, 0U, 0U, 1U, flux.revision, counters};
  return static_cast<bool>(
      cartesian_provisional_face_divergence(plan, flux, invocation));
}

bool test_shared_face_and_global_conservation(const CartesianKernelPlan& plan,
                                              FaceFluxView flux) {
  constexpr FieldId kDivergence = 50U;
  OwnedField divergence = make_field(kDivergence, kCells, 1U, 0U, 29U);
  bool passed = true;

  fill_faces(flux, 0.0);
  constexpr double shared_value = 0x1.9p+3;
  flux.x.unchecked({2, 2, 1}) = shared_value;
  passed &= expect(run_divergence(plan, as_const(flux), divergence.view),
                   "divergence consumes the one shared-face array");
  const double left = divergence.view.unchecked({1, 2, 1}, 0U);
  const double right = divergence.view.unchecked({2, 2, 1}, 0U);
  constexpr std::uint64_t sign = UINT64_C(1) << 63U;
  passed &= expect(bits(left) == bits(shared_value) &&
                       bits(right) == (bits(shared_value) ^ sign),
                   "adjacent cells register one identical face bit pattern with opposite signs");

  for (std::int32_t z = 0; z < flux.x.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.x.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.x.extents.x; ++x) {
        flux.x.unchecked({x, y, z}) =
            x == kCells.x ? 0.125 * (y + 2 * z)
                          : 0.125 * (x + y + 2 * z);
      }
    }
  }
  for (std::int32_t z = 0; z < flux.y.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.y.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.y.extents.x; ++x) {
        flux.y.unchecked({x, y, z}) =
            y == kCells.y ? -0.0625 * (2 * x + z)
                          : -0.0625 * (2 * x + y + z);
      }
    }
  }
  for (std::int32_t z = 0; z < flux.z.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.z.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.z.extents.x; ++x) {
        flux.z.unchecked({x, y, z}) =
            z == kCells.z ? 0.03125 * (x - y)
                          : 0.03125 * (x - y + 3 * z);
      }
    }
  }
  passed &= expect(run_divergence(plan, as_const(flux), divergence.view),
                   "periodic face divergence evaluates without a gather temporary");
  long double global_raw_divergence = 0.0L;
  double absolute_sum = 0.0;
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t y = 0; y < kCells.y; ++y) {
      for (std::int32_t x = 0; x < kCells.x; ++x) {
        const double value = divergence.view.unchecked({x, y, z}, 0U);
        global_raw_divergence += static_cast<long double>(value);
        absolute_sum += std::abs(value);
      }
    }
  }
  passed &= expect(close(static_cast<double>(global_raw_divergence), 0.0,
                         absolute_sum),
                   "periodic global raw divergence closes at roundoff");

  // Break periodic equality deliberately. The volume integral must then be
  // exactly the independently accumulated outward boundary mass flux.
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t y = 0; y < kCells.y; ++y) {
      flux.x.unchecked({0, y, z}) = -0.25;
      flux.x.unchecked({kCells.x, y, z}) = 0.75;
    }
  }
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t x = 0; x < kCells.x; ++x) {
      flux.y.unchecked({x, 0, z}) = 0.125;
      flux.y.unchecked({x, kCells.y, z}) = -0.375;
    }
  }
  for (std::int32_t y = 0; y < kCells.y; ++y) {
    for (std::int32_t x = 0; x < kCells.x; ++x) {
      flux.z.unchecked({x, y, 0}) = -0.5;
      flux.z.unchecked({x, y, kCells.z}) = -0.125;
    }
  }
  passed &= expect(run_divergence(plan, as_const(flux), divergence.view),
                   "open-boundary flux fixture evaluates");
  long double volume_change = 0.0L;
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t y = 0; y < kCells.y; ++y) {
      for (std::int32_t x = 0; x < kCells.x; ++x) {
        volume_change += static_cast<long double>(
            divergence.view.unchecked({x, y, z}, 0U));
      }
    }
  }
  long double boundary_flux = 0.0L;
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t y = 0; y < kCells.y; ++y) {
      boundary_flux += flux.x.unchecked({kCells.x, y, z}) -
                       flux.x.unchecked({0, y, z});
    }
  }
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t x = 0; x < kCells.x; ++x) {
      boundary_flux += flux.y.unchecked({x, kCells.y, z}) -
                       flux.y.unchecked({x, 0, z});
    }
  }
  for (std::int32_t y = 0; y < kCells.y; ++y) {
    for (std::int32_t x = 0; x < kCells.x; ++x) {
      boundary_flux += flux.z.unchecked({x, y, kCells.z}) -
                       flux.z.unchecked({x, y, 0});
    }
  }
  passed &= expect(close(static_cast<double>(volume_change),
                         static_cast<double>(boundary_flux),
                         static_cast<double>(std::abs(boundary_flux))),
                   "boundary mass flux equals the finite-volume volume change");
  return passed;
}

struct TransactionFixture {
  FieldId state{};
  StateLayers layers;
  AttemptTransaction transaction;
};

bool make_transaction_fixture(TransactionFixture& fixture) {
  FieldRegistry registry;
  FieldSchema schema;
  bool passed = expect(
      static_cast<bool>(registry.declare_field("state", 1U, 0U,
                                               fixture.state)),
      "transaction state field declares");
  passed &= expect(static_cast<bool>(registry.freeze(schema)),
                   "transaction field schema freezes");
  const std::array requests{
      ArenaFieldRequest{fixture.state, {2, 1, 1}, {0U},
                        FieldLifetime::state_layer}};
  ArenaLayout layout;
  passed &= expect(static_cast<bool>(ArenaLayout::compile(
                       schema,
                       Span<const ArenaFieldRequest>{requests.data(),
                                                     requests.size()},
                       layout)),
                   "transaction arena compiles");
  passed &= expect(static_cast<bool>(StateLayers::allocate(layout,
                                                           fixture.layers)),
                   "transaction state layers allocate");
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       fixture.layers.field_count(), 1U,
                       fixture.layers.field_count(), fixture.transaction)),
                   "transaction reserves one final-flux cache slot");
  return passed;
}

bool begin_complete_state_write(TransactionFixture& fixture,
                                RevisionDependency& dependency) {
  if (!fixture.transaction.begin(fixture.layers) ||
      !fixture.transaction.revise_trial(fixture.state)) {
    return false;
  }
  dependency = {
      AttemptTransaction::field_revision_source(fixture.state),
      fixture.transaction.trial_revision(fixture.state)};
  return dependency.revision != 0U;
}

bool reconstruct_pending(const CartesianKernelPlan& plan,
                         PendingFaceFluxView& pending, double velocity_value) {
  OwnedField rho = make_field(70U, kCells, 1U, 2U, 170U);
  OwnedField velocity = make_field(71U, kCells, 3U, 2U, 171U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        const Int3 cell{x, y, z};
        rho.view.unchecked(cell, 0U) = 1.0;
        velocity.view.unchecked(cell, 0U) = velocity_value;
        velocity.view.unchecked(cell, 1U) = velocity_value;
        velocity.view.unchecked(cell, 2U) = velocity_value;
      }
    }
  }
  const std::array<ConstFieldView, 2U> reads{as_const(rho.view),
                                             as_const(velocity.view)};
  const KernelInvocation invocation{
      Span<const ConstFieldView>{reads.data(), reads.size()}, {},
      {{0, 0, 0}, kCells}, 0U, 0U, 1U, 0U, nullptr};
  return static_cast<bool>(reconstruct_mass_flux(plan, invocation, pending));
}

bool reconstruct_diagnostic_pending(const CartesianKernelPlan& plan,
                                    PendingFaceFluxView& pending) {
  OwnedField density = make_field(75U, kCells, 1U, 2U, 175U);
  OwnedField velocity = make_field(76U, kCells, 3U, 2U, 176U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        const Int3 cell{x, y, z};
        density.view.unchecked(cell, 0U) =
            1.0 + 0.01 * x - 0.006 * y + 0.004 * z;
        velocity.view.unchecked(cell, 0U) =
            -0.4 + 0.13 * x + 0.07 * y - 0.05 * z;
        velocity.view.unchecked(cell, 1U) =
            0.3 - 0.08 * x + 0.12 * y + 0.04 * z;
        velocity.view.unchecked(cell, 2U) =
            -0.2 + 0.05 * x - 0.1 * y + 0.15 * z;
      }
    }
  }
  const std::array<ConstFieldView, 2U> reads{as_const(density.view),
                                             as_const(velocity.view)};
  const KernelInvocation invocation{
      Span<const ConstFieldView>{reads.data(), reads.size()}, {},
      {{0, 0, 0}, kCells}, 0U, 0U, 1U, 0U, nullptr};
  return static_cast<bool>(reconstruct_mass_flux(plan, invocation, pending));
}

bool publish_attempt(FinalFaceFluxWriter& writer,
                     TransactionFixture& fixture, FaceFluxStorage& storage,
                     const CartesianKernelPlan& plan, double value,
                     Status outcome, PendingFaceFluxView& pending,
                     RevisionToken& revision) {
  RevisionDependency dependency;
  if (!begin_complete_state_write(fixture, dependency) ||
      !writer.begin_pending(fixture.transaction, storage, pending)) {
    return false;
  }
  revision = pending.revision();
  if (!reconstruct_pending(plan, pending, value)) {
    return false;
  }
  const std::array dependencies{dependency};
  if (!writer.publish_pending(
          Span<const RevisionDependency>{dependencies.data(),
                                         dependencies.size()}, pending) ||
      pending.valid()) {
    return false;
  }
  const Status finished =
      fixture.transaction.collective_finish(MPI_COMM_SELF, outcome);
  if (outcome && !finished) {
    return false;
  }
  if (!outcome && (finished.code != outcome.code ||
                   finished.detail != outcome.detail)) {
    return false;
  }
  return true;
}

bool test_final_flux_authority(const CartesianKernelPlan& plan) {
  TransactionFixture fixture;
  bool passed = make_transaction_fixture(fixture);
  FaceFluxStorage storage;
  passed &= expect(static_cast<bool>(FaceFluxStorage::allocate_final(
                       kCells, storage)),
                   "final face-flux storage allocates exactly two replicas");
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
  FinalFaceFluxWriter duplicate;
  constexpr StageId kFinalPressureStage = 41U;
  constexpr RevisionSlotId kFinalFluxSlot = 0U;
  passed &= expect(authority.claim(0U, kFinalFluxSlot, fixture.transaction,
                                   writer).code == StatusCode::invalid_plan,
                   "stage zero cannot claim the final-flux authority");
  passed &= expect(static_cast<bool>(authority.claim(
                       kFinalPressureStage, kFinalFluxSlot,
                       fixture.transaction, writer)),
                   "the final pressure stage claims the only final-flux writer");
  passed &= expect(authority.claim(kFinalPressureStage, kFinalFluxSlot,
                                   fixture.transaction, duplicate).code ==
                       StatusCode::invalid_plan,
                   "a second final-flux writer is rejected explicitly");
  FinalFaceFluxAuthority competing_authority;
  passed &= expect(competing_authority.claim(
                       kFinalPressureStage, kFinalFluxSlot,
                       fixture.transaction, duplicate).code ==
                       StatusCode::invalid_plan,
                   "the transaction rejects a competing authority object");

  RevisionDependency missing_begin_dependency;
  passed &= expect(begin_complete_state_write(fixture,
                                              missing_begin_dependency),
                   "missing-begin mutation starts a complete state write");
  passed &= expect(fixture.transaction.collective_finish(
                       MPI_COMM_SELF, Status{}).code ==
                       StatusCode::invalid_plan &&
                       !fixture.transaction.active(),
                   "collective success without begin_pending rolls back");

  RevisionDependency missing_publish_dependency;
  PendingFaceFluxView unpublished;
  passed &= expect(begin_complete_state_write(fixture,
                                              missing_publish_dependency) &&
                       static_cast<bool>(writer.begin_pending(
                           fixture.transaction, storage, unpublished)) &&
                       reconstruct_pending(plan, unpublished, -4.0),
                   "missing-publish mutation constructs pending bytes");
  passed &= expect(fixture.transaction.collective_finish(
                       MPI_COMM_SELF, Status{}).code ==
                       StatusCode::invalid_plan &&
                       !fixture.transaction.active(),
                   "collective success without publish_pending rolls back");
  passed &= expect(!reconstruct_pending(plan, unpublished, -5.0),
                   "a pending handle cannot write after its attempt rolls back");

  RevisionDependency lease_retry_dependency;
  PendingFaceFluxView lease_retry;
  passed &= expect(begin_complete_state_write(fixture,
                                              lease_retry_dependency) &&
                       static_cast<bool>(writer.begin_pending(
                           fixture.transaction, storage, lease_retry)),
                   "a new attempt acquires a fresh pending lease");
  passed &= expect(!reconstruct_pending(plan, unpublished, -6.0),
                   "an old pending handle stays invalid during the next attempt");
  passed &= expect(fixture.transaction.collective_finish(
                       MPI_COMM_SELF,
                       Status{StatusCode::rejected_step, 906U}).code ==
                       StatusCode::rejected_step,
                   "lease-isolation mutation rolls back explicitly");

  RevisionDependency unauthorized_dependency;
  passed &= expect(begin_complete_state_write(fixture,
                                              unauthorized_dependency),
                   "unauthorized publication mutation starts a complete attempt");
  const std::array unauthorized_dependencies{unauthorized_dependency};
  passed &= expect(
      fixture.transaction.publish_pending_cache(
          kFinalFluxSlot,
          Span<const RevisionDependency>{unauthorized_dependencies.data(),
                                         unauthorized_dependencies.size()},
          PendingCacheStamp{99U}).code == StatusCode::invalid_plan,
      "generic cache publication cannot forge the claimed final-flux slot");
  passed &= expect(fixture.transaction.collective_finish(
                       MPI_COMM_SELF, Status{}).code ==
                       StatusCode::invalid_plan,
                   "a forged final-flux publication forces rollback");

  ConstFaceFluxView unavailable;
  passed &= expect(writer.committed(storage, unavailable).code ==
                       StatusCode::invalid_plan,
                   "no final flux is visible before the first committed attempt");
  passed &= expect(storage.view(0U, 1U, unavailable).code ==
                       StatusCode::invalid_plan,
                   "final storage cannot expose arbitrary replica bytes");

  RevisionDependency first_dependency;
  PendingFaceFluxView first_pending;
  passed &= expect(begin_complete_state_write(fixture, first_dependency) &&
                       static_cast<bool>(writer.begin_pending(
                           fixture.transaction, storage, first_pending)),
                   "first attempt borrows only the pending face replica");
  const RevisionToken first_revision = first_pending.revision();
  passed &= expect(reconstruct_pending(plan, first_pending, 3.25),
                   "only an authorized conservative kernel writes pending flux");
  const std::array first_dependencies{first_dependency};
  passed &= expect(static_cast<bool>(writer.publish_pending(
                       Span<const RevisionDependency>{first_dependencies.data(),
                                                      first_dependencies.size()},
                       first_pending)) &&
                       !first_pending.valid(),
                   "pending final flux publishes into the transaction log");
  passed &= expect(writer.committed(storage, unavailable).code ==
                       StatusCode::invalid_plan,
                   "pending bytes remain invisible before collective success");
  passed &= expect(static_cast<bool>(fixture.transaction.collective_finish(
                       MPI_COMM_SELF, Status{})),
                   "the first final-flux attempt commits collectively");
  ConstFaceFluxView committed;
  passed &= expect(static_cast<bool>(writer.committed(storage, committed)) &&
                       committed.revision == first_revision &&
                       committed.revision != 0U &&
                       committed.certificate.valid(),
                   "committed view carries the writer authority and uses handle rotation");

  OwnedField final_divergence = make_field(72U, kCells, 1U, 0U, 172U);
  OwnedField final_scalar = make_field(73U, kCells, 1U, 1U, 173U);
  OwnedField final_convection = make_field(74U, kCells, 1U, 0U, 174U);
  fill_cell_field(final_scalar, 2.0);
  const std::array<FieldView, 1U> final_divergence_writes{
      final_divergence.view};
  const std::array<ConstFieldView, 1U> final_convection_reads{
      as_const(final_scalar.view)};
  const std::array<FieldView, 1U> final_convection_writes{
      final_convection.view};
  const KernelInvocation final_divergence_call{
      {}, {final_divergence_writes.data(), final_divergence_writes.size()},
      {{0, 0, 0}, kCells}, 0U, 0U, 1U, committed.revision, nullptr};
  const KernelInvocation final_convection_call{
      {final_convection_reads.data(), final_convection_reads.size()},
      {final_convection_writes.data(), final_convection_writes.size()},
      {{0, 0, 0}, kCells}, 0U, 0U, 1U, committed.revision, nullptr};
  passed &= expect(static_cast<bool>(cartesian_face_divergence(
                       plan, committed, final_divergence_call)) &&
                       static_cast<bool>(cartesian_convection(
                           plan, ConvectionScheme::central2, committed,
                           final_convection_call)),
                   "formal conservative kernels accept only a committed flux certificate");

  std::array<FaceFluxConsumer, 6U> consumers;
  for (FaceFluxConsumer& consumer : consumers) {
    passed &= expect(static_cast<bool>(consumer.bind(committed.certificate)),
                     "every conservative consumer binds one final revision");
  }


  const auto rejects_tamper = [&committed](ConstFaceFluxView tampered) {
    FaceFluxConsumer consumer;
    return static_cast<bool>(consumer.bind(committed.certificate)) &&
           consumer.consume(tampered).code == StatusCode::invalid_plan;
  };
  ConstFaceFluxView tampered_y_base = committed;
  ++tampered_y_base.y.base;
  ConstFaceFluxView tampered_z_base = committed;
  ++tampered_z_base.z.base;
  ConstFaceFluxView tampered_y_axis = committed;
  tampered_y_axis.y.axis = CartesianAxis::x;
  ConstFaceFluxView tampered_z_axis = committed;
  tampered_z_axis.z.axis = CartesianAxis::y;
  ConstFaceFluxView tampered_y_stride = committed;
  ++tampered_y_stride.y.stride_y;
  ConstFaceFluxView tampered_z_stride = committed;
  ++tampered_z_stride.z.stride_z;
  passed &= expect(rejects_tamper(tampered_y_base) &&
                       rejects_tamper(tampered_z_base) &&
                       rejects_tamper(tampered_y_axis) &&
                       rejects_tamper(tampered_z_axis) &&
                       rejects_tamper(tampered_y_stride) &&
                       rejects_tamper(tampered_z_stride),
                   "certificate rejects y/z base, axis, and stride tampering");

  TransactionFixture foreign_fixture;
  FinalFaceFluxAuthority foreign_authority;
  FinalFaceFluxWriter foreign_writer;
  passed &= make_transaction_fixture(foreign_fixture);
  passed &= expect(static_cast<bool>(foreign_authority.claim(
                       kFinalPressureStage, kFinalFluxSlot,
                       foreign_fixture.transaction, foreign_writer)),
                   "a separate transaction can claim its own writer identity");
  RevisionDependency foreign_dependency;
  PendingFaceFluxView foreign_pending;
  passed &= expect(begin_complete_state_write(foreign_fixture,
                                              foreign_dependency),
                   "competing transaction begins normally");
  passed &= expect(foreign_writer.begin_pending(
                       foreign_fixture.transaction, storage,
                       foreign_pending).code == StatusCode::invalid_plan,
                   "a competing transaction cannot write the bound final storage");
  passed &= expect(foreign_fixture.transaction.collective_finish(
                       MPI_COMM_SELF,
                       Status{StatusCode::rejected_step, 907U}).code ==
                       StatusCode::rejected_step,
                   "competing transaction rolls back without touching final flux");
  ConstFaceFluxView stale = committed;
  ++stale.revision;
  for (FaceFluxConsumer& consumer : consumers) {
    passed &= expect(consumer.consume(stale).code == StatusCode::invalid_plan,
                     "every conservative consumer rejects a stale flux revision");
    passed &= expect(static_cast<bool>(consumer.consume(committed)) &&
                         consumer.consumed_revision() == committed.revision,
                     "density/momentum/h/species/passive/diagnostics consume one revision");
  }

  const std::uint64_t accepted_checksum = face_checksum(committed);
  const std::array<const double*, 3U> accepted_bases{
      committed.x.base, committed.y.base, committed.z.base};
  PendingFaceFluxView rejected_pending;
  RevisionToken rejected_revision{};
  passed &= expect(publish_attempt(
                       writer, fixture, storage, plan, -91.0,
                       Status{StatusCode::rejected_step, 901U},
                       rejected_pending, rejected_revision),
                   "a failed attempt is finalized as rollback");
  passed &= expect(rejected_revision != 0U,
                   "rollback mutation produced a distinct pending revision");
  ConstFaceFluxView after_rejection;
  passed &= expect(static_cast<bool>(writer.committed(storage,
                                                      after_rejection)) &&
                       after_rejection.revision == committed.revision &&
                       face_checksum(after_rejection) == accepted_checksum &&
                       after_rejection.x.base == accepted_bases[0U] &&
                       after_rejection.y.base == accepted_bases[1U] &&
                       after_rejection.z.base == accepted_bases[2U],
                   "rollback leaves active bytes, pointers, and revision untouched");

  RevisionDependency next_dependency;
  PendingFaceFluxView next_pending;
  passed &= expect(begin_complete_state_write(fixture, next_dependency) &&
                       static_cast<bool>(writer.begin_pending(
                           fixture.transaction, storage, next_pending)),
                   "a retry gets the inactive face replica");
  const RevisionToken next_revision = next_pending.revision();
  passed &= expect(reconstruct_pending(plan, next_pending, 7.5),
                   "replacement pending flux is reconstructed completely");
  const std::array next_dependencies{next_dependency};
  passed &= expect(static_cast<bool>(writer.publish_pending(
                       Span<const RevisionDependency>{next_dependencies.data(),
                                                      next_dependencies.size()},
                       next_pending)) &&
                       !next_pending.valid(),
                   "replacement final flux enters the pending log");
  ConstFaceFluxView before_switch;
  passed &= expect(static_cast<bool>(writer.committed(storage,
                                                      before_switch)) &&
                       face_checksum(before_switch) == accepted_checksum,
                   "active readers see the old complete revision during a retry");
  passed &= expect(static_cast<bool>(fixture.transaction.collective_finish(
                       MPI_COMM_SELF, Status{})),
                   "replacement attempt reaches collective success");
  ConstFaceFluxView before_finalize;
  passed &= expect(static_cast<bool>(writer.committed(storage,
                                                      before_finalize)) &&
                       face_checksum(before_finalize) != accepted_checksum,
                   "collective commit atomically rotates state and final flux");
  ConstFaceFluxView replaced;
  passed &= expect(static_cast<bool>(writer.committed(storage, replaced)) &&
                       replaced.revision == next_revision &&
                       replaced.revision != committed.revision &&
                       replaced.x.base != accepted_bases[0U] &&
                       face_checksum(replaced) != accepted_checksum,
                   "one atomic handle switch exposes the complete new revision");
  return passed;
}

enum class LifetimeObject : std::uint8_t {
  transaction,
  writer,
  storage,
};

constexpr std::array<std::array<LifetimeObject, 3U>, 6U>
    kLifetimeDestructionOrders{{
        {LifetimeObject::transaction, LifetimeObject::writer,
         LifetimeObject::storage},
        {LifetimeObject::transaction, LifetimeObject::storage,
         LifetimeObject::writer},
        {LifetimeObject::writer, LifetimeObject::transaction,
         LifetimeObject::storage},
        {LifetimeObject::writer, LifetimeObject::storage,
         LifetimeObject::transaction},
        {LifetimeObject::storage, LifetimeObject::transaction,
         LifetimeObject::writer},
        {LifetimeObject::storage, LifetimeObject::writer,
         LifetimeObject::transaction},
    }};

struct LifetimeFixture {
  FieldId state{};
  StateLayers layers;
  std::unique_ptr<AttemptTransaction> transaction;
  std::unique_ptr<FinalFaceFluxWriter> writer;
  std::unique_ptr<FaceFluxStorage> storage;
};

bool make_lifetime_fixture(LifetimeFixture& fixture) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("lifetime_state", 1U, 0U, fixture.state) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{
      ArenaFieldRequest{fixture.state, {2, 1, 1}, {0U},
                        FieldLifetime::state_layer}};
  ArenaLayout layout;
  if (!ArenaLayout::compile(
          schema,
          Span<const ArenaFieldRequest>{requests.data(), requests.size()},
          layout) ||
      !StateLayers::allocate(layout, fixture.layers)) {
    return false;
  }

  fixture.transaction = std::make_unique<AttemptTransaction>();
  fixture.writer = std::make_unique<FinalFaceFluxWriter>();
  fixture.storage = std::make_unique<FaceFluxStorage>();
  if (!AttemptTransaction::create(fixture.layers.field_count(), 1U,
                                  fixture.layers.field_count(),
                                  *fixture.transaction) ||
      !FaceFluxStorage::allocate_final(kCells, *fixture.storage)) {
    return false;
  }
  FinalFaceFluxAuthority authority;
  return static_cast<bool>(authority.claim(73U, 0U, *fixture.transaction,
                                           *fixture.writer));
}

std::array<std::size_t, 3U> state_handles(const LifetimeFixture& fixture) {
  return {fixture.layers.handle(StateRole::accepted_n),
          fixture.layers.handle(StateRole::accepted_n_minus_one),
          fixture.layers.handle(StateRole::trial)};
}

bool begin_lifetime_attempt(LifetimeFixture& fixture,
                            const CartesianKernelPlan& plan,
                            PendingFaceFluxView& pending,
                            bool publish) {
  if (!fixture.transaction->begin(fixture.layers) ||
      !fixture.transaction->revise_trial(fixture.state)) {
    return false;
  }
  const RevisionDependency dependency{
      AttemptTransaction::field_revision_source(fixture.state),
      fixture.transaction->trial_revision(fixture.state)};
  if (dependency.revision == 0U ||
      !fixture.writer->begin_pending(*fixture.transaction, *fixture.storage,
                                     pending) ||
      !reconstruct_pending(plan, pending, publish ? 6.0 : -6.0)) {
    return false;
  }
  if (!publish) {
    return true;
  }
  const std::array dependencies{dependency};
  return static_cast<bool>(fixture.writer->publish_pending(
             Span<const RevisionDependency>{dependencies.data(),
                                             dependencies.size()},
             pending)) &&
         !pending.valid();
}

void destroy_lifetime_object(LifetimeFixture& fixture,
                             LifetimeObject object) {
  switch (object) {
    case LifetimeObject::transaction:
      fixture.transaction.reset();
      break;
    case LifetimeObject::writer:
      fixture.writer.reset();
      break;
    case LifetimeObject::storage:
      fixture.storage.reset();
      break;
  }
}

bool run_lifetime_destruction_order(
    const CartesianKernelPlan& plan,
    const std::array<LifetimeObject, 3U>& order, bool publish) {
  LifetimeFixture fixture;
  PendingFaceFluxView pending;
  if (!make_lifetime_fixture(fixture) ||
      !begin_lifetime_attempt(fixture, plan, pending, publish)) {
    return false;
  }
  const std::array<std::size_t, 3U> handles_before =
      state_handles(fixture);

  destroy_lifetime_object(fixture, order[0U]);
  bool passed = state_handles(fixture) == handles_before;
  if (order[0U] != LifetimeObject::transaction) {
    const Status finished =
        fixture.transaction->collective_finish(MPI_COMM_SELF, Status{});
    passed = passed && finished.code == StatusCode::invalid_plan &&
             !fixture.transaction->active() &&
             state_handles(fixture) == handles_before;
  }
  destroy_lifetime_object(fixture, order[1U]);
  passed = passed && state_handles(fixture) == handles_before;
  destroy_lifetime_object(fixture, order[2U]);
  passed = passed && state_handles(fixture) == handles_before &&
           !pending.valid();
  return passed;
}

bool test_final_flux_lifetime_fail_closed(
    const CartesianKernelPlan& plan) {
  bool unpublished_passed = true;
  bool published_passed = true;
  for (const auto& order : kLifetimeDestructionOrders) {
    unpublished_passed =
        run_lifetime_destruction_order(plan, order, false) &&
        unpublished_passed;
    published_passed = run_lifetime_destruction_order(plan, order, true) &&
                       published_passed;
  }
  bool passed = expect(
      unpublished_passed,
      "all 3! active-unpublished transaction/writer/storage destruction orders fail closed");
  passed &= expect(
      published_passed,
      "all 3! published-pending transaction/writer/storage destruction orders fail closed");

  LifetimeFixture fixture;
  passed &= expect(make_lifetime_fixture(fixture),
                   "pending-view lifetime fixture initializes");
  if (!fixture.transaction || !fixture.writer || !fixture.storage) {
    return false;
  }
  const std::array<std::size_t, 3U> handles_before =
      state_handles(fixture);
  auto pending = std::make_unique<PendingFaceFluxView>();
  passed &= expect(begin_lifetime_attempt(fixture, plan, *pending, false),
                   "unpublished pending view acquires an active lease");
  pending.reset();
  const Status finished =
      fixture.transaction->collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(finished.code == StatusCode::invalid_plan &&
                       !fixture.transaction->active() &&
                       state_handles(fixture) == handles_before,
                   "destroying an unpublished pending view forces rollback without handle rotation");
  return passed;
}

bool test_point_convection_diagnostic(const CartesianKernelPlan& plan) {
  LifetimeFixture fixture;
  bool passed = expect(make_lifetime_fixture(fixture),
                       "point-diagnostic final-flux fixture initializes");
  if (!passed || !fixture.transaction || !fixture.writer ||
      !fixture.storage) {
    return false;
  }

  PendingFaceFluxView pending;
  passed &= expect(fixture.transaction->begin(fixture.layers) &&
                       fixture.transaction->revise_trial(fixture.state),
                   "point-diagnostic attempt begins");
  const RevisionDependency dependency{
      AttemptTransaction::field_revision_source(fixture.state),
      fixture.transaction->trial_revision(fixture.state)};
  passed &= expect(dependency.revision != 0U &&
                       static_cast<bool>(fixture.writer->begin_pending(
                           *fixture.transaction, *fixture.storage, pending)) &&
                       reconstruct_diagnostic_pending(plan, pending),
                   "point-diagnostic fixture reconstructs varied mass flux");
  const std::array dependencies{dependency};
  passed &= expect(static_cast<bool>(fixture.writer->publish_pending(
                       Span<const RevisionDependency>{dependencies.data(),
                                                      dependencies.size()},
                       pending)) &&
                       !pending.valid(),
                   "point-diagnostic fixture publishes final mass flux");
  passed &= expect(static_cast<bool>(fixture.transaction->collective_finish(
                       MPI_COMM_SELF, Status{})),
                   "point-diagnostic fixture commits final mass flux");

  ConstFaceFluxView committed;
  passed &= expect(static_cast<bool>(fixture.writer->committed(
                       *fixture.storage, committed)) &&
                       committed.certificate.valid(),
                   "point-diagnostic fixture exposes a certified flux");
  if (!passed) {
    return false;
  }

  constexpr std::uint8_t kComponents = 2U;
  constexpr Int3 kProbeCell{2, 2, 2};
  OwnedField transported =
      make_field(77U, kCells, kComponents, 2U, 177U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        const double xx = static_cast<double>(x);
        const double yy = static_cast<double>(y);
        const double zz = static_cast<double>(z);
        transported.view.unchecked({x, y, z}, 0U) =
            2.0 + 0.07 * xx - 0.04 * yy + 0.03 * zz +
            0.006 * xx * yy - 0.003 * yy * zz;
        transported.view.unchecked({x, y, z}, 1U) = 1.0;
      }
    }
  }
  OwnedField production = make_field(78U, kCells, kComponents, 0U, 178U);
  const std::array<ConstFieldView, 1U> reads{as_const(transported.view)};
  const std::array<FieldView, 1U> writes{production.view};
  const KernelInvocation invocation{
      Span<const ConstFieldView>{reads.data(), reads.size()},
      Span<const FieldView>{writes.data(), writes.size()},
      {{0, 0, 0}, kCells}, 0U, 0U, kComponents, committed.revision, nullptr};

  for (const ConvectionScheme scheme :
       {ConvectionScheme::central2, ConvectionScheme::tvd2}) {
    std::fill(production.allocation.begin(), production.allocation.end(), 0.0);
    passed &= expect(static_cast<bool>(cartesian_convection(
                         plan, scheme, committed, invocation)),
                     scheme == ConvectionScheme::tvd2
                         ? "TVD2 production convection evaluates"
                         : "donor-cell production convection evaluates");

    ConvectionPointDiagnostic scalar_diagnostic;
    ConvectionPointDiagnostic donor_diagnostic;
    passed &= expect(static_cast<bool>(diagnose_cartesian_convection_point(
                         plan, scheme, committed, as_const(transported.view),
                         0U, kProbeCell, scalar_diagnostic)) &&
                         static_cast<bool>(diagnose_cartesian_convection_point(
                             plan, scheme, committed,
                             as_const(transported.view), 1U, kProbeCell,
                             donor_diagnostic)),
                     scheme == ConvectionScheme::tvd2
                         ? "TVD2 point diagnostics evaluate"
                         : "donor-cell point diagnostics evaluate");
    const double production_scalar =
        production.view.unchecked(kProbeCell, 0U);
    const double production_mass =
        production.view.unchecked(kProbeCell, 1U);
    passed &= expect(close(scalar_diagnostic.divergence, production_scalar,
                           std::abs(production_scalar)) &&
                         close(donor_diagnostic.divergence, production_mass,
                               std::abs(production_mass)) &&
                         close(scalar_diagnostic.mass_divergence,
                               donor_diagnostic.mass_divergence,
                               std::abs(donor_diagnostic.mass_divergence)) &&
                         close(scalar_diagnostic.mass_divergence,
                               production_mass,
                               std::abs(production_mass)),
                     scheme == ConvectionScheme::tvd2
                         ? "TVD2 point divergence and mass divergence match the production field"
                         : "donor-cell point divergence and mass divergence match the production field");

    if (scheme == ConvectionScheme::tvd2) {
      const double scale = std::max(
          {1.0, std::abs(scalar_diagnostic.selected_face_value),
           std::abs(scalar_diagnostic.selected_donor_minimum),
           std::abs(scalar_diagnostic.selected_donor_maximum)});
      passed &= expect(scalar_diagnostic.face_envelope_checked &&
                           scalar_diagnostic.face_envelope_valid,
                       "TVD2 point diagnostic checks a valid donor envelope");
      passed &= expect(
          scalar_diagnostic.maximum_face_envelope_violation <=
              128.0 * std::numeric_limits<double>::epsilon() * scale,
          "TVD2 face-envelope violation remains at roundoff");
      passed &= expect(
          scalar_diagnostic.selected_face_value >=
                  scalar_diagnostic.selected_donor_minimum -
                      128.0 * std::numeric_limits<double>::epsilon() * scale &&
              scalar_diagnostic.selected_face_value <=
                  scalar_diagnostic.selected_donor_maximum +
                      128.0 * std::numeric_limits<double>::epsilon() * scale,
          "TVD2 selected face remains inside its donor-cell envelope");
    }
  }
  return passed;
}

bool test_exact_kernel_counters_and_allocations(
    const CartesianKernelPlan& plan, FaceFluxView flux) {
  constexpr std::uint8_t kComponents = 2U;
  const KernelBox full_box{{0, 0, 0}, kCells};
  const KernelBox subtile{{2, 1, 1}, {3, 2, 2}};
  const std::uint64_t cells = cell_count(kCells);
  const std::uint64_t tile_cells = cell_count(subtile.cells);
  const std::uint64_t full_faces = reconstructed_face_count(full_box);
  const std::uint64_t tile_faces = reconstructed_face_count(subtile);

  OwnedField density = make_field(80U, kCells, 1U, 2U, 180U);
  OwnedField velocity = make_field(81U, kCells, 3U, 2U, 181U);
  OwnedField transported = make_field(82U, kCells, kComponents, 2U, 182U);
  OwnedField gradient = make_field(83U, kCells, 6U, 0U, 183U);
  OwnedField convection = make_field(84U, kCells, kComponents, 0U, 184U);
  OwnedField diffusion = make_field(85U, kCells, kComponents, 0U, 185U);
  OwnedField diffusivity = make_field(86U, kCells, 1U, 1U, 186U);
  fill_cell_field(density, 2.0);
  fill_cell_field(velocity, 0.75);
  fill_cell_field(transported, -1.0);
  fill_cell_field(diffusivity, 1.5);

  const std::array<ConstFieldView, 1U> transported_reads{
      as_const(transported.view)};
  const std::array<FieldView, 1U> gradient_writes{gradient.view};
  const std::array<FieldView, 1U> convection_writes{convection.view};
  const std::array<FieldView, 1U> diffusion_writes{diffusion.view};
  const std::array<ConstFieldView, 2U> mass_reads{as_const(density.view),
                                                 as_const(velocity.view)};

  KernelCounters gradient_count;
  KernelInvocation gradient_call{
      Span<const ConstFieldView>{transported_reads.data(),
                                 transported_reads.size()},
      Span<const FieldView>{gradient_writes.data(), gradient_writes.size()},
      full_box, 0U, 0U, kComponents, 0U, &gradient_count};
  KernelCounters full_mass_count;
  KernelInvocation full_mass_call{
      Span<const ConstFieldView>{mass_reads.data(), mass_reads.size()}, {},
      full_box, 0U, 0U, 1U, 0U, &full_mass_count};
  KernelCounters tile_mass_count;
  KernelInvocation tile_mass_call = full_mass_call;
  tile_mass_call.box = subtile;
  tile_mass_call.counters = &tile_mass_count;
  KernelCounters diffusion_count;
  KernelInvocation diffusion_call{
      Span<const ConstFieldView>{transported_reads.data(),
                                 transported_reads.size()},
      Span<const FieldView>{diffusion_writes.data(), diffusion_writes.size()},
      full_box, 0U, 0U, kComponents, 0U, &diffusion_count};
  KernelCounters limited_count;
  KernelInvocation limited_call{
      Span<const ConstFieldView>{transported_reads.data(),
                                 transported_reads.size()},
      Span<const FieldView>{convection_writes.data(),
                            convection_writes.size()},
      full_box, 0U, 0U, kComponents, flux.revision, &limited_count};
  KernelCounters tvd_count;
  KernelInvocation tvd_call = limited_call;
  tvd_call.counters = &tvd_count;

  bool passed = expect(static_cast<bool>(cartesian_gradient(plan,
                                                             gradient_call)),
                       "gradient evaluates for exact counter accounting");
  passed &= expect(static_cast<bool>(reconstruct_mass_flux(
                       plan, full_mass_call, flux)),
                   "full-box mass flux reconstruction evaluates");
  passed &= expect(static_cast<bool>(reconstruct_mass_flux(
                       plan, tile_mass_call, flux)),
                   "subtile mass flux reconstruction evaluates");
  passed &= expect(static_cast<bool>(cartesian_diffusion(
                       plan, as_const(diffusivity.view), diffusion_call)),
                   "diffusion evaluates for exact counter accounting");
  passed &= expect(static_cast<bool>(cartesian_provisional_convection(
                       plan, ConvectionScheme::limited_central2,
                       as_const(flux), limited_call)),
                   "limited convection evaluates for exact counter accounting");
  passed &= expect(static_cast<bool>(cartesian_provisional_convection(
                       plan, ConvectionScheme::tvd2, as_const(flux),
                       tvd_call)),
                   "TVD convection evaluates for exact counter accounting");

  passed &= expect(
      gradient_count.invocations == 1U && gradient_count.cells == cells &&
          gradient_count.faces == 0U &&
          gradient_count.logical_bytes_read ==
              7U * kComponents * cells * sizeof(double) &&
          gradient_count.logical_bytes_written ==
              3U * kComponents * cells * sizeof(double),
      "gradient counters match exact logical scalar traffic");
  passed &= expect(
      full_mass_count.invocations == 1U &&
          full_mass_count.cells == cells &&
          full_mass_count.faces == full_faces &&
          full_mass_count.logical_bytes_read ==
              4U * full_faces * sizeof(double) &&
          full_mass_count.logical_bytes_written ==
              full_faces * sizeof(double),
      "full-box mass counters count each uniquely owned face once");
  passed &= expect(
      tile_mass_count.invocations == 1U &&
          tile_mass_count.cells == tile_cells &&
          tile_mass_count.faces == tile_faces &&
          tile_mass_count.logical_bytes_read ==
              4U * tile_faces * sizeof(double) &&
          tile_mass_count.logical_bytes_written ==
              tile_faces * sizeof(double),
      "subtile mass counters exclude faces owned by a neighboring tile");
  passed &= expect(
      diffusion_count.invocations == 1U && diffusion_count.cells == cells &&
          diffusion_count.faces == 6U * cells &&
          diffusion_count.logical_bytes_read ==
              (12U + 12U * kComponents) * cells * sizeof(double) &&
          diffusion_count.logical_bytes_written ==
              kComponents * cells * sizeof(double),
      "diffusion counters include coefficient and transported traffic");
  passed &= expect(
      limited_count.invocations == 1U && limited_count.cells == cells &&
          limited_count.faces == 6U * cells &&
          limited_count.logical_bytes_read ==
              (6U + 24U * kComponents) * cells * sizeof(double) &&
          limited_count.logical_bytes_written ==
              kComponents * cells * sizeof(double),
      "limited convection counters include both reconstructed face states");
  passed &= expect(
      tvd_count.invocations == 1U && tvd_count.cells == cells &&
          tvd_count.faces == 6U * cells &&
          tvd_count.logical_bytes_read ==
              (6U + 24U * kComponents) * cells * sizeof(double) &&
          tvd_count.logical_bytes_written ==
              kComponents * cells * sizeof(double),
      "TVD counters include both reconstructed face states");

  KernelCounters repeated;
  gradient_call.counters = &repeated;
  full_mass_call.counters = &repeated;
  tile_mass_call.counters = &repeated;
  diffusion_call.counters = &repeated;
  limited_call.counters = &repeated;
  tvd_call.counters = &repeated;
  bool hot_ok = true;
  std::size_t allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t repetition = 0U; repetition < 100U; ++repetition) {
      hot_ok = hot_ok && static_cast<bool>(cartesian_gradient(
                             plan, gradient_call));
      hot_ok = hot_ok && static_cast<bool>(reconstruct_mass_flux(
                             plan, full_mass_call, flux));
      hot_ok = hot_ok && static_cast<bool>(reconstruct_mass_flux(
                             plan, tile_mass_call, flux));
      hot_ok = hot_ok && static_cast<bool>(cartesian_diffusion(
                             plan, as_const(diffusivity.view),
                             diffusion_call));
      hot_ok = hot_ok && static_cast<bool>(cartesian_provisional_convection(
                             plan, ConvectionScheme::limited_central2,
                             as_const(flux), limited_call));
      hot_ok = hot_ok && static_cast<bool>(cartesian_provisional_convection(
                             plan, ConvectionScheme::tvd2, as_const(flux),
                             tvd_call));
    }
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const std::uint64_t expected_cells =
      100U * (5U * cells + tile_cells);
  const std::uint64_t expected_faces =
      100U * (full_faces + tile_faces + 18U * cells);
  const std::uint64_t expected_read_doubles =
      100U * (7U * kComponents * cells + 4U * full_faces +
              4U * tile_faces +
              (12U + 12U * kComponents) * cells +
              2U * (6U + 24U * kComponents) * cells);
  const std::uint64_t expected_written_doubles =
      100U * (3U * kComponents * cells + full_faces + tile_faces +
              3U * kComponents * cells);
  passed &= expect(hot_ok && allocations == 0U,
                   "100 repeated gradient/flux/transport kernels allocate zero bytes");
  passed &= expect(
      repeated.invocations == 600U && repeated.cells == expected_cells &&
          repeated.faces == expected_faces &&
          repeated.logical_bytes_read ==
              expected_read_doubles * sizeof(double) &&
          repeated.logical_bytes_written ==
              expected_written_doubles * sizeof(double),
      "repeated mixed kernels retain exact additive work counters");
  return passed;
}

bool test_alias_mutations(const CartesianKernelPlan& plan,
                          FaceFluxView flux) {
  constexpr std::uint8_t kComponents = 2U;
  const KernelBox full_box{{0, 0, 0}, kCells};
  OwnedField diffusivity = make_field(90U, kCells, 1U, 1U, 190U);
  fill_cell_field(diffusivity, 1.25);
  bool passed = true;

  OwnedField gradient_alias = make_field(91U, kCells, 6U, 2U, 191U);
  fill_cell_field(gradient_alias, -2.0);
  FieldView shifted_gradient_output = gradient_alias.view;
  ++shifted_gradient_output.base;
  const std::array<ConstFieldView, 1U> gradient_reads{
      as_const(gradient_alias.view)};
  const std::array<FieldView, 1U> gradient_writes{shifted_gradient_output};
  KernelCounters gradient_counters{3U, 5U, 7U, 11U, 13U};
  const KernelCounters gradient_before = gradient_counters;
  const std::uint64_t gradient_checksum = allocation_checksum(gradient_alias);
  const KernelInvocation gradient_call{
      Span<const ConstFieldView>{gradient_reads.data(), gradient_reads.size()},
      Span<const FieldView>{gradient_writes.data(), gradient_writes.size()},
      full_box, 0U, 0U, kComponents, 0U, &gradient_counters};
  passed &= expect(cartesian_gradient(plan, gradient_call).code ==
                           StatusCode::invalid_plan &&
                       allocation_checksum(gradient_alias) ==
                           gradient_checksum &&
                       same_counters(gradient_counters, gradient_before),
                   "gradient rejects overlapping shifted output without mutation");

  OwnedField diffusion_alias = make_field(92U, kCells, kComponents, 2U, 192U);
  fill_cell_field(diffusion_alias, 0.5);
  FieldView shifted_diffusion_output = diffusion_alias.view;
  ++shifted_diffusion_output.base;
  const std::array<ConstFieldView, 1U> diffusion_reads{
      as_const(diffusion_alias.view)};
  const std::array<FieldView, 1U> diffusion_writes{shifted_diffusion_output};
  KernelCounters diffusion_counters{17U, 19U, 23U, 29U, 31U};
  const KernelCounters diffusion_before = diffusion_counters;
  const std::uint64_t diffusion_checksum =
      allocation_checksum(diffusion_alias);
  const KernelInvocation diffusion_call{
      Span<const ConstFieldView>{diffusion_reads.data(), diffusion_reads.size()},
      Span<const FieldView>{diffusion_writes.data(), diffusion_writes.size()},
      full_box, 0U, 0U, kComponents, 0U, &diffusion_counters};
  passed &= expect(
      cartesian_diffusion(plan, as_const(diffusivity.view),
                          diffusion_call).code == StatusCode::invalid_plan &&
          allocation_checksum(diffusion_alias) == diffusion_checksum &&
          same_counters(diffusion_counters, diffusion_before),
      "diffusion rejects overlapping shifted output without mutation");

  OwnedField convection_alias = make_field(93U, kCells, kComponents, 2U, 193U);
  fill_cell_field(convection_alias, 1.75);
  FieldView shifted_convection_output = convection_alias.view;
  ++shifted_convection_output.base;
  const std::array<ConstFieldView, 1U> convection_reads{
      as_const(convection_alias.view)};
  const std::array<FieldView, 1U> convection_writes{
      shifted_convection_output};
  KernelCounters convection_counters{37U, 41U, 43U, 47U, 53U};
  const KernelCounters convection_before = convection_counters;
  const std::uint64_t convection_checksum =
      allocation_checksum(convection_alias);
  const KernelInvocation convection_call{
      Span<const ConstFieldView>{convection_reads.data(), convection_reads.size()},
      Span<const FieldView>{convection_writes.data(), convection_writes.size()},
      full_box, 0U, 0U, kComponents, flux.revision, &convection_counters};
  passed &= expect(
      cartesian_provisional_convection(
          plan, ConvectionScheme::central2, as_const(flux),
                           convection_call).code ==
              StatusCode::invalid_plan &&
          allocation_checksum(convection_alias) == convection_checksum &&
          same_counters(convection_counters, convection_before),
      "convection rejects overlapping shifted output without mutation");
  return passed;
}

bool test_convection_revision_and_hot_counters(
    const CartesianKernelPlan& plan, FaceFluxView flux) {
  constexpr FieldId kTransported = 61U;
  constexpr FieldId kDivergence = 62U;
  constexpr FieldId kConvection = 63U;
  OwnedField transported = make_field(kTransported, kCells, 2U, 1U, 71U);
  OwnedField divergence = make_field(kDivergence, kCells, 1U, 0U, 72U);
  OwnedField convection = make_field(kConvection, kCells, 2U, 0U, 73U);
  for (std::int32_t z = -1; z <= kCells.z; ++z) {
    for (std::int32_t y = -1; y <= kCells.y; ++y) {
      for (std::int32_t x = -1; x <= kCells.x; ++x) {
        transported.view.unchecked({x, y, z}, 0U) = 2.0;
        transported.view.unchecked({x, y, z}, 1U) = -3.5;
      }
    }
  }
  for (std::int32_t z = 0; z < flux.x.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.x.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.x.extents.x; ++x) {
        flux.x.unchecked({x, y, z}) = 0.125 * (x + y + z);
      }
    }
  }
  for (std::int32_t z = 0; z < flux.y.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.y.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.y.extents.x; ++x) {
        flux.y.unchecked({x, y, z}) = -0.0625 * (2 * x + y + z);
      }
    }
  }
  for (std::int32_t z = 0; z < flux.z.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.z.extents.y; ++y) {
      for (std::int32_t x = 0; x < flux.z.extents.x; ++x) {
        flux.z.unchecked({x, y, z}) = 0.03125 * (x - y + z);
      }
    }
  }

  std::array<FieldView, 1U> divergence_writes{divergence.view};
  std::array<ConstFieldView, 1U> convection_reads{
      as_const(transported.view)};
  std::array<FieldView, 1U> convection_writes{convection.view};
  KernelCounters one_divergence;
  KernelCounters one_convection;
  KernelInvocation divergence_call{
      {}, Span<const FieldView>{divergence_writes.data(),
                                divergence_writes.size()},
      {{0, 0, 0}, kCells}, 0U, 0U, 1U, flux.revision, &one_divergence};
  KernelInvocation convection_call{
      Span<const ConstFieldView>{convection_reads.data(),
                                 convection_reads.size()},
      Span<const FieldView>{convection_writes.data(),
                            convection_writes.size()},
      {{0, 0, 0}, kCells}, 0U, 0U, 2U, flux.revision, &one_convection};
  bool passed = expect(static_cast<bool>(cartesian_provisional_face_divergence(
                           plan, as_const(flux), divergence_call)),
                       "face divergence accepts the exact required revision");
  passed &= expect(static_cast<bool>(cartesian_provisional_convection(
                       plan, ConvectionScheme::central2, as_const(flux),
                       convection_call)),
                   "central convection consumes the same face revision");
  passed &= expect(cartesian_face_divergence(
                       plan, as_const(flux), divergence_call).code ==
                       StatusCode::invalid_plan &&
                       cartesian_convection(
                           plan, ConvectionScheme::central2, as_const(flux),
                           convection_call).code == StatusCode::invalid_plan,
                   "formal conservative kernels reject provisional workspace flux");
  for (std::int32_t z = 0; z < kCells.z; ++z) {
    for (std::int32_t y = 0; y < kCells.y; ++y) {
      for (std::int32_t x = 0; x < kCells.x; ++x) {
        const Int3 index{x, y, z};
        const double div = divergence.view.unchecked(index, 0U);
        passed &= expect(close(convection.view.unchecked(index, 0U),
                               2.0 * div, std::abs(div)),
                         "constant component convection equals q times divergence");
        passed &= expect(close(convection.view.unchecked(index, 1U),
                               -3.5 * div, std::abs(div)),
                         "all transported components reuse the same face flux");
      }
    }
  }
  KernelInvocation stale_call = convection_call;
  ++stale_call.required_face_flux_revision;
  passed &= expect(cartesian_provisional_convection(
                       plan, ConvectionScheme::central2, as_const(flux),
                       stale_call).code == StatusCode::invalid_plan,
                   "convection rejects a stale required final-flux revision");

  const std::uint64_t cells = cell_count(kCells);
  passed &= expect(one_divergence.invocations == 1U &&
                       one_divergence.cells == cells &&
                       one_divergence.faces == 6U * cells &&
                       one_divergence.logical_bytes_read ==
                           6U * cells * sizeof(double) &&
                       one_divergence.logical_bytes_written ==
                           cells * sizeof(double),
                   "divergence counters report the exact minimum logical traffic");
  passed &= expect(one_convection.invocations == 1U &&
                       one_convection.cells == cells &&
                       one_convection.faces == 6U * cells &&
                       one_convection.logical_bytes_read ==
                           (6U + 12U * 2U) * cells * sizeof(double) &&
                       one_convection.logical_bytes_written ==
                           2U * cells * sizeof(double),
                   "convection counters report the exact minimum logical traffic");

  KernelCounters repeated;
  divergence_call.counters = &repeated;
  convection_call.counters = &repeated;
  bool hot_ok = true;
  std::size_t allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::size_t repetition = 0U; repetition < 100U; ++repetition) {
      hot_ok = hot_ok && static_cast<bool>(cartesian_provisional_face_divergence(
                             plan, as_const(flux), divergence_call));
      hot_ok = hot_ok && static_cast<bool>(cartesian_provisional_convection(
                             plan, ConvectionScheme::central2,
                             as_const(flux), convection_call));
    }
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(hot_ok && allocations == 0U,
                   "100 repeated conservative kernels allocate zero bytes");
  passed &= expect(
      repeated.invocations == 200U && repeated.cells == 200U * cells &&
          repeated.faces == 1200U * cells &&
          repeated.logical_bytes_read == 100U * (6U + 6U + 24U) *
                                                 cells * sizeof(double) &&
          repeated.logical_bytes_written ==
              100U * (1U + 2U) * cells * sizeof(double),
      "100 repeated kernels preserve strict linear byte and work counters");
  return passed;
}

bool test_frozen_target_convection_faces(const CartesianKernelPlan& plan,
                                         FaceFluxStorage& storage,
                                         FaceFluxView target_flux) {
  FaceFluxView frozen_storage;
  bool passed = expect(
      static_cast<bool>(storage.workspace_view(1U, 18U, frozen_storage)),
      "a disjoint preallocated face replica is available for frozen states");
  OwnedField enthalpy = make_field(95U, kCells, 1U, 2U, 195U);
  OwnedField convection = make_field(96U, kCells, 1U, 0U, 196U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        enthalpy.view.unchecked({x, y, z}, 0U) =
            3.0 + 0.2 * x * x - 0.15 * y + 0.07 * z * z +
            0.03 * x * y;
      }
    }
  }
  const auto fill_directional_flux = [](FaceFieldView faces,
                                        double scale) noexcept {
    for (std::int32_t z = 0; z < faces.extents.z; ++z) {
      for (std::int32_t y = 0; y < faces.extents.y; ++y) {
        for (std::int32_t x = 0; x < faces.extents.x; ++x) {
          const int parity = (x + 2 * y + 3 * z) & 1;
          faces.unchecked({x, y, z}) =
              (parity == 0 ? 1.0 : -1.0) *
              scale * (1.0 + 0.1 * x + 0.05 * y + 0.025 * z);
        }
      }
    }
  };
  fill_directional_flux(target_flux.x, 0.19);
  fill_directional_flux(target_flux.y, 0.13);
  fill_directional_flux(target_flux.z, 0.11);

  const std::array<ConvectionScheme, 3U> schemes{
      ConvectionScheme::central2, ConvectionScheme::limited_central2,
      ConvectionScheme::tvd2};
  for (ConvectionScheme scheme : schemes) {
    fill_faces(frozen_storage, -991.0);
    FrozenConvectionFaceField frozen;
    std::size_t allocations = std::numeric_limits<std::size_t>::max();
    Status frozen_status;
    {
      allocation_observer::Guard guard;
      frozen_status = freeze_cartesian_target_convection_faces(
          plan, scheme, as_const(target_flux), as_const(enthalpy.view), 0U,
          {UINT64_C(0x8f86d9f314a2c571), 219U},
          {frozen_storage.x, frozen_storage.y, frozen_storage.z}, frozen);
      allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    passed &= expect(static_cast<bool>(frozen_status) && frozen.valid() &&
                         allocations == 0U,
                     "target convection freezes one exact allocation-free face state");

    const std::array<ConstFieldView, 1U> reads{as_const(enthalpy.view)};
    const std::array<FieldView, 1U> writes{convection.view};
    const KernelInvocation invocation{
        {reads.data(), reads.size()}, {writes.data(), writes.size()},
        {{0, 0, 0}, kCells}, 0U, 0U, 1U, target_flux.revision, nullptr};
    passed &= expect(static_cast<bool>(cartesian_target_convection(
                         plan, scheme, as_const(target_flux), invocation)),
                     "target convection residual evaluates against the frozen flux");
    for (std::int32_t z = 0; z < kCells.z; ++z) {
      for (std::int32_t y = 0; y < kCells.y; ++y) {
        for (std::int32_t x = 0; x < kCells.x; ++x) {
          const Int3 cell{x, y, z};
          const double reconstructed =
              target_flux.x.unchecked({x + 1, y, z}) *
                  frozen.x.unchecked({x + 1, y, z}) -
              target_flux.x.unchecked({x, y, z}) *
                  frozen.x.unchecked({x, y, z}) +
              target_flux.y.unchecked({x, y + 1, z}) *
                  frozen.y.unchecked({x, y + 1, z}) -
              target_flux.y.unchecked({x, y, z}) *
                  frozen.y.unchecked({x, y, z}) +
              target_flux.z.unchecked({x, y, z + 1}) *
                  frozen.z.unchecked({x, y, z + 1}) -
              target_flux.z.unchecked({x, y, z}) *
                  frozen.z.unchecked({x, y, z});
          passed &= expect(
              close(reconstructed,
                    convection.view.unchecked(cell, 0U), reconstructed),
              "frozen face products reproduce the production convection residual");
        }
      }
    }
  }

  fill_faces(frozen_storage, 817.0);
  const std::uint64_t before = face_checksum(as_const(frozen_storage));
  const double saved = enthalpy.view.unchecked({0, 0, 0}, 0U);
  enthalpy.view.unchecked({0, 0, 0}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  FrozenConvectionFaceField rejected;
  const Status rejected_status = freeze_cartesian_target_convection_faces(
      plan, ConvectionScheme::limited_central2, as_const(target_flux),
      as_const(enthalpy.view), 0U,
      {UINT64_C(0x8f86d9f314a2c571), 219U},
      {frozen_storage.x, frozen_storage.y, frozen_storage.z}, rejected);
  passed &= expect(rejected_status.code == StatusCode::numerical_failure &&
                       !rejected.valid() &&
                       face_checksum(as_const(frozen_storage)) == before,
                   "non-finite frozen reconstruction rejects atomically");
  enthalpy.view.unchecked({0, 0, 0}, 0U) = saved;
  return passed;
}

bool test_frozen_target_convection_directional_derivative(
    const CartesianKernelPlan& plan, ConstFaceFluxView target_flux) {
  FaceFluxStorage frozen_owner;
  FaceFluxStorage derivative_owner;
  FaceFluxStorage plus_owner;
  FaceFluxStorage minus_owner;
  FaceFluxView frozen_faces;
  FaceFluxView derivative_faces;
  FaceFluxView plus_faces;
  FaceFluxView minus_faces;
  bool passed = expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          kCells, 1U, frozen_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, derivative_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, plus_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, minus_owner)) &&
          static_cast<bool>(frozen_owner.workspace_view(0U, 301U,
                                                        frozen_faces)) &&
          static_cast<bool>(derivative_owner.workspace_view(
              0U, 302U, derivative_faces)) &&
          static_cast<bool>(plus_owner.workspace_view(0U, 303U,
                                                      plus_faces)) &&
          static_cast<bool>(minus_owner.workspace_view(0U, 304U,
                                                       minus_faces)),
      "directional-derivative test owns four disjoint face workspaces");
  if (!passed) return false;

  bool saw_positive_flux = false;
  bool saw_negative_flux = false;
  const auto inspect_flux_signs = [&](ConstFaceFieldView faces) {
    for (std::int32_t z = 0; z < faces.extents.z; ++z) {
      for (std::int32_t y = 0; y < faces.extents.y; ++y) {
        for (std::int32_t x = 0; x < faces.extents.x; ++x) {
          const double value = faces.unchecked({x, y, z});
          saw_positive_flux = saw_positive_flux || value > 0.0;
          saw_negative_flux = saw_negative_flux || value < 0.0;
        }
      }
    }
  };
  inspect_flux_signs(target_flux.x);
  inspect_flux_signs(target_flux.y);
  inspect_flux_signs(target_flux.z);
  passed &= expect(saw_positive_flux && saw_negative_flux,
                   "directional derivative exercises both donor directions");

  OwnedField target = make_field(97U, kCells, 1U, 2U, 401U);
  OwnedField variation = make_field(98U, kCells, 1U, 2U, 402U);
  OwnedField plus = make_field(99U, kCells, 1U, 2U, 403U);
  OwnedField minus = make_field(100U, kCells, 1U, 2U, 404U);
  constexpr double epsilon = 2.0e-6;
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        const double xd = static_cast<double>(x);
        const double yd = static_cast<double>(y);
        const double zd = static_cast<double>(z);
        const double target_value =
            7.0 + 0.173 * xd + 0.0191 * xd * xd + 0.00231 * xd * xd * xd +
            0.127 * yd + 0.0137 * yd * yd + 0.00173 * yd * yd * yd +
            0.091 * zd + 0.0109 * zd * zd + 0.00137 * zd * zd * zd +
            0.0047 * xd * yd - 0.0031 * yd * zd + 0.0023 * xd * zd;
        const double direction =
            -0.41 + 0.071 * xd * xd - 0.053 * yd + 0.029 * zd * zd +
            0.017 * xd * yd - 0.011 * xd * zd + 0.007 * yd * zd;
        target.view.unchecked({x, y, z}, 0U) = target_value;
        variation.view.unchecked({x, y, z}, 0U) = direction;
        plus.view.unchecked({x, y, z}, 0U) =
            target_value + epsilon * direction;
        minus.view.unchecked({x, y, z}, 0U) =
            target_value - epsilon * direction;
      }
    }
  }

  const FrozenConvectionContext context{
      UINT64_C(0x8f86d9f314a2c571), 501U};
  const std::array<ConvectionScheme, 3U> schemes{
      ConvectionScheme::central2, ConvectionScheme::limited_central2,
      ConvectionScheme::tvd2};
  for (ConvectionScheme scheme : schemes) {
    FrozenConvectionFaceField frozen;
    passed &= expect(
        static_cast<bool>(freeze_cartesian_target_convection_faces(
            plan, scheme, target_flux, as_const(target.view), 0U, context,
            {frozen_faces.x, frozen_faces.y, frozen_faces.z}, frozen)),
        "smooth target reconstruction freezes before differentiation");

    FrozenConvectionFaceDirectionalDerivative derivative;
    std::size_t allocations = std::numeric_limits<std::size_t>::max();
    Status differentiated;
    {
      allocation_observer::Guard guard;
      differentiated =
          differentiate_frozen_cartesian_target_convection_faces(
              plan, scheme, target_flux, as_const(target.view), 0U, context,
              FrozenConvectionLinearizationPolicy::classical_active_branch,
              frozen, as_const(variation.view), 0U,
              {derivative_faces.x, derivative_faces.y, derivative_faces.z},
              derivative);
      allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    passed &= expect(static_cast<bool>(differentiated) && derivative.valid() &&
                         derivative.reconstruction != frozen.reconstruction &&
                         derivative.branch_authority != 0U &&
                         derivative.policy ==
                             FrozenConvectionLinearizationPolicy::
                                 classical_active_branch &&
                         derivative.classical_everywhere &&
                         derivative.generalized_face_count == 0U &&
                         allocations == 0U,
                     "frozen target branch differentiates allocation-free");

    FrozenConvectionFaceField plus_frozen;
    FrozenConvectionFaceField minus_frozen;
    passed &= expect(
        static_cast<bool>(freeze_cartesian_target_convection_faces(
            plan, scheme, target_flux, as_const(plus.view), 0U,
            {context.collective_semantics, 502U},
            {plus_faces.x, plus_faces.y, plus_faces.z}, plus_frozen)) &&
            static_cast<bool>(freeze_cartesian_target_convection_faces(
                plan, scheme, target_flux, as_const(minus.view), 0U,
                {context.collective_semantics, 503U},
                {minus_faces.x, minus_faces.y, minus_faces.z}, minus_frozen)),
        "centered-FD reference evaluates exact target reconstruction twice");
    passed &= expect(plus_frozen.reconstruction == frozen.reconstruction &&
                         minus_frozen.reconstruction ==
                             frozen.reconstruction &&
                         plus_frozen.revision != frozen.revision &&
                         plus_frozen.local_binding != frozen.local_binding,
                     "collective reconstruction identity excludes exact "
                     "rank-local target and storage revisions");

    const auto compare_axis = [&](ConstFaceFieldView actual,
                                  ConstFaceFieldView positive,
                                  ConstFaceFieldView negative) {
      bool axis_passed = true;
      for (std::int32_t z = 0; z < actual.extents.z; ++z) {
        for (std::int32_t y = 0; y < actual.extents.y; ++y) {
          for (std::int32_t x = 0; x < actual.extents.x; ++x) {
            const Int3 face{x, y, z};
            const double expected =
                (positive.unchecked(face) - negative.unchecked(face)) /
                (2.0 * epsilon);
            axis_passed &= direction_close(actual.unchecked(face), expected);
          }
        }
      }
      return axis_passed;
    };
    passed &= expect(
        compare_axis(derivative.x, plus_frozen.x, minus_frozen.x) &&
            compare_axis(derivative.y, plus_frozen.y, minus_frozen.y) &&
            compare_axis(derivative.z, plus_frozen.z, minus_frozen.z),
        "central2/limited-central2/TVD2 frozen branch derivative matches "
        "centered FD for positive and negative face flux");
    FrozenConvectionFaceDirectionalDerivative plus_derivative;
    passed &= expect(
        static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
            plan, scheme, target_flux, as_const(plus.view), 0U,
            {context.collective_semantics, 502U},
            FrozenConvectionLinearizationPolicy::classical_active_branch,
            plus_frozen, as_const(variation.view), 0U,
            {derivative_faces.x, derivative_faces.y, derivative_faces.z},
            plus_derivative)) &&
            plus_derivative.reconstruction == derivative.reconstruction &&
            plus_derivative.branch_authority !=
                derivative.branch_authority &&
            plus_derivative.local_binding != derivative.local_binding,
        "exact target/closure/view identity stays rank-local to branch and "
        "storage authority");
    if (scheme != ConvectionScheme::central2) {
      FrozenConvectionFaceField nonlinear_direction;
      passed &= expect(
          static_cast<bool>(freeze_cartesian_target_convection_faces(
              plan, scheme, target_flux, as_const(variation.view), 0U,
              {context.collective_semantics, 504U},
              {plus_faces.x, plus_faces.y, plus_faces.z},
              nonlinear_direction)),
          "nonlinear variation reconstruction is available as a negative "
          "control");
      bool differs_from_nonlinear_variation = false;
      const auto find_difference = [&](ConstFaceFieldView actual,
                                       ConstFaceFieldView nonlinear) {
        for (std::int32_t z = 0; z < actual.extents.z; ++z) {
          for (std::int32_t y = 0; y < actual.extents.y; ++y) {
            for (std::int32_t x = 0; x < actual.extents.x; ++x) {
              const Int3 face{x, y, z};
              differs_from_nonlinear_variation =
                  differs_from_nonlinear_variation ||
                  !direction_close(actual.unchecked(face),
                                   nonlinear.unchecked(face));
            }
          }
        }
      };
      find_difference(derivative.x, nonlinear_direction.x);
      find_difference(derivative.y, nonlinear_direction.y);
      find_difference(derivative.z, nonlinear_direction.z);
      passed &= expect(
          differs_from_nonlinear_variation,
          "variation is not passed through a fresh nonlinear limiter");
    }
  }
  return passed;
}

bool test_frozen_direction_analytic_branches(
    const CartesianKernelPlan& plan) {
  FaceFluxStorage flux_owner;
  FaceFluxStorage frozen_owner;
  FaceFluxStorage derivative_owner;
  FaceFluxView flux;
  FaceFluxView frozen_faces;
  FaceFluxView derivative_faces;
  bool passed = expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          kCells, 1U, flux_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, frozen_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, derivative_owner)) &&
          static_cast<bool>(flux_owner.workspace_view(0U, 551U, flux)) &&
          static_cast<bool>(frozen_owner.workspace_view(0U, 552U,
                                                        frozen_faces)) &&
          static_cast<bool>(derivative_owner.workspace_view(
              0U, 553U, derivative_faces)),
      "analytic branch oracle owns disjoint face storage");
  if (!passed) return false;
  fill_faces(flux, 1.0);
  for (std::int32_t z = 0; z < flux.x.extents.z; ++z) {
    for (std::int32_t y = 0; y < flux.x.extents.y; ++y) {
      flux.x.unchecked({2, y, z}) = -1.0;
    }
  }

  const std::array<double, 9U> x_increments{
      1.0, -1.0, -1.7, 4.0, 1.0, 1.4, 0.7, 2.0, 0.9};
  const std::array<double, 8U> y_increments{
      1.0, 2.0, 0.7, 4.0, 1.1, 0.45, 2.2, 0.8};
  const std::array<double, 7U> z_increments{
      1.0, 4.0, 1.2, 0.5, 2.4, 0.9, 1.7};
  std::array<double, 10U> x_values{};
  std::array<double, 9U> y_values{};
  std::array<double, 8U> z_values{};
  for (std::size_t i = 0U; i < x_increments.size(); ++i)
    x_values[i + 1U] = x_values[i] + x_increments[i];
  for (std::size_t i = 0U; i < y_increments.size(); ++i)
    y_values[i + 1U] = y_values[i] + y_increments[i];
  for (std::size_t i = 0U; i < z_increments.size(); ++i)
    z_values[i + 1U] = z_values[i] + z_increments[i];

  OwnedField target = make_field(103U, kCells, 1U, 2U, 751U);
  OwnedField variation = make_field(104U, kCells, 1U, 2U, 752U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        target.view.unchecked({x, y, z}, 0U) =
            x_values[static_cast<std::size_t>(x + 2)] +
            y_values[static_cast<std::size_t>(y + 2)] +
            z_values[static_cast<std::size_t>(z + 2)];
        const double xd = static_cast<double>(x);
        const double yd = static_cast<double>(y);
        const double zd = static_cast<double>(z);
        variation.view.unchecked({x, y, z}, 0U) =
            10.0 + 0.31 * xd * xd - 0.23 * yd * yd +
            0.41 * zd * zd + 0.07 * xd * yd - 0.05 * yd * zd;
      }
    }
  }
  const FrozenConvectionContext context{
      UINT64_C(0x9b715cf2a80634de), 753U};
  FrozenConvectionFaceField frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, frozen)),
      "analytic target branch field freezes exactly");
  FrozenConvectionFaceDirectionalDerivative derivative;
  passed &= expect(
      static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch,
          frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          derivative)) &&
          derivative.valid() && derivative.classical_everywhere,
      "all analytic target branches have classical derivatives");

  const Int3 x_zero_face{0, 1, 1};
  const double x_zero_oracle =
      variation.view.unchecked({-1, 1, 1}, 0U);
  const Int3 y_centred_face{1, 0, 1};
  const double y_left = variation.view.unchecked({1, -1, 1}, 0U);
  const double y_delta_left =
      y_left - variation.view.unchecked({1, -2, 1}, 0U);
  const double y_delta_right =
      variation.view.unchecked({1, 0, 1}, 0U) - y_left;
  const double y_centred_oracle =
      y_left + 0.25 * (y_delta_left + y_delta_right);
  const Int3 z_left_delta_face{1, 1, 0};
  const double z_left = variation.view.unchecked({1, 1, -1}, 0U);
  const double z_left_delta_oracle =
      z_left +
      (z_left - variation.view.unchecked({1, 1, -2}, 0U));
  const Int3 x_right_delta_face{2, 1, 1};
  const double x_right = variation.view.unchecked({2, 1, 1}, 0U);
  const double x_right_delta_oracle =
      x_right -
      (variation.view.unchecked({3, 1, 1}, 0U) - x_right);
  passed &= expect(
      close(derivative.x.unchecked(x_zero_face), x_zero_oracle,
            x_zero_oracle) &&
          close(derivative.y.unchecked(y_centred_face), y_centred_oracle,
                y_centred_oracle) &&
          close(derivative.z.unchecked(z_left_delta_face),
                z_left_delta_oracle, z_left_delta_oracle) &&
          close(derivative.x.unchecked(x_right_delta_face),
                x_right_delta_oracle, x_right_delta_oracle),
      "independent oracle covers zero/centred/left-delta/right-delta across "
      "three axes, both donors, and boundary ghosts");

  FrozenConvectionFaceDirectionalDerivative semismooth_on_smooth;
  passed &= expect(
      static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::
              semismooth_generalized_zero_slope,
          frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          semismooth_on_smooth)) &&
          semismooth_on_smooth.classical_everywhere &&
          semismooth_on_smooth.generalized_face_count == 0U &&
          semismooth_on_smooth.reconstruction != derivative.reconstruction,
      "linearization policy changes collective semantics even when all "
      "faces are classical");

  FrozenConvectionFaceField limited_frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, limited_frozen)),
      "mixed limited-central target branches freeze");
  FrozenConvectionFaceDirectionalDerivative limited_derivative;
  passed &= expect(
      static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch,
          limited_frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          limited_derivative)),
      "mixed limited-central target branches differentiate");
  const double mixed_left = variation.view.unchecked({-1, 1, 1}, 0U);
  const double mixed_right = variation.view.unchecked({0, 1, 1}, 0U);
  const double mixed_right_delta_left = mixed_right - mixed_left;
  const double mixed_right_delta_right =
      variation.view.unchecked({1, 1, 1}, 0U) - mixed_right;
  const double mixed_right_centred =
      0.5 * (mixed_right_delta_left + mixed_right_delta_right);
  const double mixed_limited_oracle =
      0.5 * (mixed_left + mixed_right - 0.5 * mixed_right_centred);
  passed &= expect(
      close(limited_derivative.x.unchecked(x_zero_face),
            mixed_limited_oracle, mixed_limited_oracle),
      "limited-central independently combines zero and centred branches");
  return passed;
}

bool test_constant_semismooth_and_zero_limiter(
    const CartesianKernelPlan& classical_plan,
    const CartesianKernelPlan& zero_limiter_plan) {
  FaceFluxStorage flux_owner;
  FaceFluxStorage frozen_owner;
  FaceFluxStorage derivative_owner;
  FaceFluxView flux;
  FaceFluxView frozen_faces;
  FaceFluxView derivative_faces;
  bool passed = expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          kCells, 1U, flux_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, frozen_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, derivative_owner)) &&
          static_cast<bool>(flux_owner.workspace_view(0U, 561U, flux)) &&
          static_cast<bool>(frozen_owner.workspace_view(0U, 562U,
                                                        frozen_faces)) &&
          static_cast<bool>(derivative_owner.workspace_view(
              0U, 563U, derivative_faces)),
      "constant-policy tests own disjoint face storage");
  if (!passed) return false;
  const auto fill_signed_flux = [](FaceFieldView faces, int axis) {
    for (std::int32_t z = 0; z < faces.extents.z; ++z) {
      for (std::int32_t y = 0; y < faces.extents.y; ++y) {
        for (std::int32_t x = 0; x < faces.extents.x; ++x) {
          faces.unchecked({x, y, z}) =
              ((x + 2 * y + 3 * z + axis) & 1) == 0 ? 0.7 : -0.9;
        }
      }
    }
  };
  fill_signed_flux(flux.x, 0);
  fill_signed_flux(flux.y, 1);
  fill_signed_flux(flux.z, 2);
  OwnedField target = make_field(105U, kCells, 1U, 2U, 761U);
  OwnedField variation = make_field(106U, kCells, 1U, 2U, 762U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        target.view.unchecked({x, y, z}, 0U) = 7.25;
        variation.view.unchecked({x, y, z}, 0U) =
            -0.4 + 0.11 * x - 0.07 * y + 0.13 * z + 0.017 * x * y;
      }
    }
  }
  const FrozenConvectionContext context{
      UINT64_C(0x36a19df807c254be), 763U};
  const std::uint64_t face_count =
      static_cast<std::uint64_t>(kCells.x + 1) * kCells.y * kCells.z +
      static_cast<std::uint64_t>(kCells.x) * (kCells.y + 1) * kCells.z +
      static_cast<std::uint64_t>(kCells.x) * kCells.y * (kCells.z + 1);
  const auto analytic_faces = [&](ConvectionScheme scheme,
                                  const FrozenConvectionFaceDirectionalDerivative&
                                      derivative) {
    bool correct = true;
    const auto check_axis = [&](ConstFaceFieldView actual,
                                ConstFaceFieldView rates,
                                CartesianAxis axis) {
      bool axis_correct = true;
      for (std::int32_t z = 0; z < actual.extents.z; ++z) {
        for (std::int32_t y = 0; y < actual.extents.y; ++y) {
          for (std::int32_t x = 0; x < actual.extents.x; ++x) {
            const Int3 right{x, y, z};
            Int3 left = right;
            if (axis == CartesianAxis::x)
              --left.x;
            else if (axis == CartesianAxis::y)
              --left.y;
            else
              --left.z;
            const double left_value =
                variation.view.unchecked(left, 0U);
            const double right_value =
                variation.view.unchecked(right, 0U);
            const double expected =
                scheme == ConvectionScheme::limited_central2
                    ? 0.5 * (left_value + right_value)
                    : (rates.unchecked(right) >= 0.0 ? left_value
                                                     : right_value);
            axis_correct &=
                close(actual.unchecked(right), expected, expected);
          }
        }
      }
      return axis_correct;
    };
    correct &= check_axis(derivative.x, as_const(flux.x), CartesianAxis::x);
    correct &= check_axis(derivative.y, as_const(flux.y), CartesianAxis::y);
    correct &= check_axis(derivative.z, as_const(flux.z), CartesianAxis::z);
    return correct;
  };

  const std::array<ConvectionScheme, 2U> schemes{
      ConvectionScheme::limited_central2, ConvectionScheme::tvd2};
  for (ConvectionScheme scheme : schemes) {
    FrozenConvectionFaceField frozen;
    passed &= expect(
        static_cast<bool>(freeze_cartesian_target_convection_faces(
            classical_plan, scheme, as_const(flux), as_const(target.view), 0U,
            context, {frozen_faces.x, frozen_faces.y, frozen_faces.z},
            frozen)),
        "constant target freezes under nonzero limiter");
    fill_faces(derivative_faces, 919.0);
    const std::uint64_t before = face_checksum(as_const(derivative_faces));
    FrozenConvectionFaceDirectionalDerivative classical_rejected;
    const Status classical_status =
        differentiate_frozen_cartesian_target_convection_faces(
            classical_plan, scheme, as_const(flux), as_const(target.view), 0U,
            context,
            FrozenConvectionLinearizationPolicy::classical_active_branch,
            frozen, as_const(variation.view), 0U,
            {derivative_faces.x, derivative_faces.y, derivative_faces.z},
            classical_rejected);
    passed &= expect(classical_status.code == StatusCode::numerical_failure &&
                         !classical_rejected.valid() &&
                         face_checksum(as_const(derivative_faces)) == before,
                     "classical constant limiter kink fails before writes");

    FrozenConvectionFaceDirectionalDerivative generalized;
    passed &= expect(
        static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
            classical_plan, scheme, as_const(flux), as_const(target.view), 0U,
            context,
            FrozenConvectionLinearizationPolicy::
                semismooth_generalized_zero_slope,
            frozen, as_const(variation.view), 0U,
            {derivative_faces.x, derivative_faces.y, derivative_faces.z},
            generalized)) &&
            generalized.valid() && !generalized.classical_everywhere &&
            generalized.generalized_face_count == face_count &&
            analytic_faces(scheme, generalized),
        "semismooth constant field uses zero slope on every face");

    FrozenConvectionFaceField zero_frozen;
    passed &= expect(
        static_cast<bool>(freeze_cartesian_target_convection_faces(
            zero_limiter_plan, scheme, as_const(flux), as_const(target.view),
            0U, context,
            {frozen_faces.x, frozen_faces.y, frozen_faces.z}, zero_frozen)),
        "constant target freezes under limiter zero");
    FrozenConvectionFaceDirectionalDerivative zero_derivative;
    passed &= expect(
        static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
            zero_limiter_plan, scheme, as_const(flux), as_const(target.view),
            0U, context,
            FrozenConvectionLinearizationPolicy::classical_active_branch,
            zero_frozen, as_const(variation.view), 0U,
            {derivative_faces.x, derivative_faces.y, derivative_faces.z},
            zero_derivative)) &&
            zero_derivative.valid() && zero_derivative.classical_everywhere &&
            zero_derivative.generalized_face_count == 0U &&
            analytic_faces(scheme, zero_derivative),
        "limiter zero is globally classical with zero slope derivative");
  }
  return passed;
}

bool test_subnormal_same_sign_limiter(const CartesianKernelPlan& plan) {
  FaceFluxStorage flux_owner;
  FaceFluxStorage frozen_owner;
  FaceFluxStorage derivative_owner;
  FaceFluxView flux;
  FaceFluxView frozen_faces;
  FaceFluxView derivative_faces;
  bool passed = expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          kCells, 1U, flux_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, frozen_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, derivative_owner)) &&
          static_cast<bool>(flux_owner.workspace_view(0U, 571U, flux)) &&
          static_cast<bool>(frozen_owner.workspace_view(0U, 572U,
                                                        frozen_faces)) &&
          static_cast<bool>(derivative_owner.workspace_view(
              0U, 573U, derivative_faces)),
      "subnormal limiter test owns disjoint face storage");
  if (!passed) return false;
  fill_faces(flux, 1.0);
  const double denormal = std::numeric_limits<double>::denorm_min();
  passed &= expect(denormal > 0.0,
                   "platform exposes positive binary64 subnormals");
  OwnedField target = make_field(107U, kCells, 1U, 2U, 771U);
  OwnedField variation = make_field(108U, kCells, 1U, 2U, 772U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        const double coefficient =
            static_cast<double>(40 + 2 * x + 4 * y + 8 * z);
        target.view.unchecked({x, y, z}, 0U) = denormal * coefficient;
        variation.view.unchecked({x, y, z}, 0U) =
            3.0 + 0.21 * x * x - 0.17 * y + 0.09 * z * z;
      }
    }
  }
  const FrozenConvectionContext context{
      UINT64_C(0x7ca15e2904b863df), 773U};
  FrozenConvectionFaceField frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, frozen)),
      "same-sign subnormal TVD target freezes");
  const Int3 selected_face{1, 1, 1};
  const double target_left = target.view.unchecked({0, 1, 1}, 0U);
  const double target_right = target.view.unchecked({1, 1, 1}, 0U);
  const double target_oracle = 0.5 * (target_left + target_right);
  passed &= expect(bits(frozen.x.unchecked(selected_face)) ==
                       bits(target_oracle) &&
                       bits(target_oracle) != bits(target_left),
                   "production minmod preserves same-sign subnormal slope");

  FrozenConvectionFaceDirectionalDerivative derivative;
  const Status derivative_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch,
          frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          derivative);
  const bool derivative_available =
      static_cast<bool>(derivative_status) && derivative.valid() &&
      derivative.classical_everywhere;
  passed &= expect(
      derivative_available,
      "same-sign subnormal target selects a classical centred branch");
  const double variation_left =
      variation.view.unchecked({0, 1, 1}, 0U);
  const double variation_delta_left =
      variation_left - variation.view.unchecked({-1, 1, 1}, 0U);
  const double variation_delta_right =
      variation.view.unchecked({1, 1, 1}, 0U) - variation_left;
  const double direction_oracle =
      variation_left +
      0.25 * (variation_delta_left + variation_delta_right);
  passed &= expect(derivative_available &&
                       close(derivative.x.unchecked(selected_face),
                             direction_oracle, direction_oracle),
                   "subnormal production value and selected derivative "
                   "branch share one sign-safe minmod decision");
  return passed;
}

bool test_stretched_direction_oracle(const KernelFixture& fixture) {
  FaceFluxStorage flux_owner;
  FaceFluxStorage frozen_owner;
  FaceFluxStorage derivative_owner;
  FaceFluxView flux;
  FaceFluxView frozen_faces;
  FaceFluxView derivative_faces;
  bool passed = expect(
      fixture.geometry.kind() == GeometryKind::tensor_stretched &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, flux_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, frozen_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, derivative_owner)) &&
          static_cast<bool>(flux_owner.workspace_view(0U, 581U, flux)) &&
          static_cast<bool>(frozen_owner.workspace_view(0U, 582U,
                                                        frozen_faces)) &&
          static_cast<bool>(derivative_owner.workspace_view(
              0U, 583U, derivative_faces)),
      "stretched oracle owns a tensor metric and disjoint face storage");
  if (!passed) return false;
  const auto nonuniform = [](Span<const double> widths) {
    for (std::size_t i = 1U; i < widths.size; ++i) {
      if (bits(widths.data[i]) != bits(widths.data[0U])) return true;
    }
    return false;
  };
  passed &= expect(nonuniform(fixture.geometry.x().widths()) ||
                       nonuniform(fixture.geometry.y().widths()) ||
                       nonuniform(fixture.geometry.z().widths()),
                   "tensor metric contains genuinely stretched widths");
  fill_faces(flux, 0.6);
  OwnedField target = make_field(109U, kCells, 1U, 2U, 781U);
  OwnedField variation = make_field(110U, kCells, 1U, 2U, 782U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        target.view.unchecked({x, y, z}, 0U) =
            2.0 + 0.13 * x - 0.09 * y + 0.07 * z + 0.011 * x * z;
        variation.view.unchecked({x, y, z}, 0U) =
            -1.0 + 0.17 * x * x + 0.08 * y - 0.12 * z * z +
            0.019 * x * y;
      }
    }
  }
  const FrozenConvectionContext context{
      UINT64_C(0x51e29ab8047c63df), 783U};
  FrozenConvectionFaceField frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          fixture.kernels, ConvectionScheme::central2, as_const(flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, frozen)),
      "stretched central target freezes");
  FrozenConvectionFaceDirectionalDerivative derivative;
  passed &= expect(
      static_cast<bool>(differentiate_frozen_cartesian_target_convection_faces(
          fixture.kernels, ConvectionScheme::central2, as_const(flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch,
          frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          derivative)) &&
          derivative.classical_everywhere,
      "stretched central target differentiates classically");

  const auto interpolate = [](Span<const double> faces,
                              Span<const double> centres,
                              std::int32_t face, double left,
                              double right) {
    const double face_coordinate = faces.data[face];
    const double left_distance = face_coordinate - centres.data[face - 1];
    const double right_distance = centres.data[face] - face_coordinate;
    return (right_distance * left + left_distance * right) /
           (left_distance + right_distance);
  };
  const Int3 x_face{2, 1, 1};
  const Int3 y_face{1, 2, 1};
  const Int3 z_face{1, 1, 2};
  const double x_oracle = interpolate(
      fixture.geometry.x().faces(), fixture.geometry.x().centres(), 2,
      variation.view.unchecked({1, 1, 1}, 0U),
      variation.view.unchecked({2, 1, 1}, 0U));
  const double y_oracle = interpolate(
      fixture.geometry.y().faces(), fixture.geometry.y().centres(), 2,
      variation.view.unchecked({1, 1, 1}, 0U),
      variation.view.unchecked({1, 2, 1}, 0U));
  const double z_oracle = interpolate(
      fixture.geometry.z().faces(), fixture.geometry.z().centres(), 2,
      variation.view.unchecked({1, 1, 1}, 0U),
      variation.view.unchecked({1, 1, 2}, 0U));
  passed &= expect(
      close(derivative.x.unchecked(x_face), x_oracle, x_oracle) &&
          close(derivative.y.unchecked(y_face), y_oracle, y_oracle) &&
          close(derivative.z.unchecked(z_face), z_oracle, z_oracle),
      "independent metric-distance oracle covers all stretched axes");
  return passed;
}

bool test_frozen_target_direction_fail_closed(
    const CartesianKernelPlan& plan, FaceFluxView target_flux) {
  FaceFluxStorage frozen_owner;
  FaceFluxStorage derivative_owner;
  FaceFluxView frozen_faces;
  FaceFluxView derivative_faces;
  bool passed = expect(
      static_cast<bool>(FaceFluxStorage::allocate_workspace(
          kCells, 1U, frozen_owner)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              kCells, 1U, derivative_owner)) &&
          static_cast<bool>(frozen_owner.workspace_view(0U, 601U,
                                                        frozen_faces)) &&
          static_cast<bool>(derivative_owner.workspace_view(
              0U, 602U, derivative_faces)),
      "failure tests own disjoint frozen and derivative face storage");
  if (!passed) return false;

  OwnedField target = make_field(101U, kCells, 1U, 2U, 701U);
  OwnedField variation = make_field(102U, kCells, 1U, 2U, 702U);
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        const double xd = static_cast<double>(x);
        const double yd = static_cast<double>(y);
        const double zd = static_cast<double>(z);
        target.view.unchecked({x, y, z}, 0U) =
            4.0 + 0.181 * xd + 0.0217 * xd * xd +
            0.00213 * xd * xd * xd + 0.119 * yd + 0.0113 * yd * yd +
            0.00191 * yd * yd * yd + 0.083 * zd + 0.0149 * zd * zd +
            0.00119 * zd * zd * zd + 0.0037 * xd * yd;
        variation.view.unchecked({x, y, z}, 0U) =
            -0.2 + 0.037 * xd * xd - 0.029 * yd + 0.023 * zd * zd +
            0.013 * xd * zd;
      }
    }
  }
  const FrozenConvectionContext context{
      UINT64_C(0xd204e5a6713cb98f), 703U};
  FrozenConvectionFaceField frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(target_flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, frozen)),
      "smooth limited target freezes for failure-path tests");

  fill_faces(derivative_faces, 811.0);
  const std::uint64_t untouched = face_checksum(as_const(derivative_faces));
  const auto rejected_without_write =
      [&](const FrozenConvectionFaceField& authority,
          ConstFieldView direction,
          FrozenConvectionFaceOutput output) {
        FrozenConvectionFaceDirectionalDerivative rejected;
        const Status status =
            differentiate_frozen_cartesian_target_convection_faces(
                plan, ConvectionScheme::limited_central2,
                as_const(target_flux), as_const(target.view), 0U, context,
                FrozenConvectionLinearizationPolicy::classical_active_branch,
                authority, direction, 0U, output, rejected);
        return status.code == StatusCode::invalid_plan && !rejected.valid() &&
               face_checksum(as_const(derivative_faces)) == untouched;
      };
  FrozenConvectionFaceField stale_revision = frozen;
  ++stale_revision.revision;
  FrozenConvectionFaceField stale_local = frozen;
  ++stale_local.local_binding;
  passed &= expect(
      rejected_without_write(
          stale_revision, as_const(variation.view),
          {derivative_faces.x, derivative_faces.y, derivative_faces.z}) &&
          rejected_without_write(
              stale_local, as_const(variation.view),
              {derivative_faces.x, derivative_faces.y,
               derivative_faces.z}),
      "stale frozen revision/local binding fail closed");

  const double saved_frozen_byte = frozen_faces.z.unchecked(
      {frozen_faces.z.extents.x - 1, frozen_faces.z.extents.y - 1,
       frozen_faces.z.extents.z - 1});
  frozen_faces.z.unchecked(
      {frozen_faces.z.extents.x - 1, frozen_faces.z.extents.y - 1,
       frozen_faces.z.extents.z - 1}) = saved_frozen_byte + 0.25;
  passed &= expect(
      rejected_without_write(
          frozen, as_const(variation.view),
          {derivative_faces.x, derivative_faces.y, derivative_faces.z}),
      "raw frozen face-byte mutation fails numeric sealing");
  frozen_faces.z.unchecked(
      {frozen_faces.z.extents.x - 1, frozen_faces.z.extents.y - 1,
       frozen_faces.z.extents.z - 1}) = saved_frozen_byte;

  const std::vector<double> saved_target_bytes = target.allocation;
  for (double& value : target.allocation) value += 0.375;
  passed &= expect(
      rejected_without_write(
          frozen, as_const(variation.view),
          {derivative_faces.x, derivative_faces.y, derivative_faces.z}),
      "same-branch raw target shift fails frozen numeric comparison");
  std::copy(saved_target_bytes.begin(), saved_target_bytes.end(),
            target.allocation.begin());

  const std::uint64_t flux_before = face_checksum(as_const(target_flux));
  FrozenConvectionFaceDirectionalDerivative alias_rejected;
  const Status alias_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(target_flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch, frozen,
          as_const(variation.view), 0U,
          {target_flux.x, target_flux.y, target_flux.z}, alias_rejected);
  passed &= expect(alias_status.code == StatusCode::invalid_plan &&
                       !alias_rejected.valid() &&
                       face_checksum(as_const(target_flux)) == flux_before,
                   "direction output may not alias target flux storage");

  const double saved_nan = variation.view.unchecked({0, 0, -2}, 0U);
  variation.view.unchecked({0, 0, -2}, 0U) =
      std::numeric_limits<double>::quiet_NaN();
  FrozenConvectionFaceDirectionalDerivative nonfinite_rejected;
  const Status nonfinite_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(target_flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch, frozen,
          as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          nonfinite_rejected);
  passed &= expect(nonfinite_status.code == StatusCode::numerical_failure &&
                       !nonfinite_rejected.valid() &&
                       face_checksum(as_const(derivative_faces)) == untouched,
                   "z-only non-finite direction rejects before x/y/z writes");
  variation.view.unchecked({0, 0, -2}, 0U) = saved_nan;

  const std::array<double, 10U> tied_x{
      0.0, 1.0, 4.0, 8.75, 14.2, 20.9, 28.8, 37.95, 48.4, 60.2};
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        target.view.unchecked({x, y, z}, 0U) =
            tied_x[static_cast<std::size_t>(x + 2)] + 0.125 * y +
            0.0625 * z;
      }
    }
  }
  ++target.view.revision;
  FrozenConvectionFaceField tied_frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(target_flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, tied_frozen)),
      "finite limiter-tie target still has an exact nonlinear reconstruction");
  FrozenConvectionFaceDirectionalDerivative tie_rejected;
  const Status tie_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(target_flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch,
          tied_frozen,
          as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          tie_rejected);
  passed &= expect(tie_status.code == StatusCode::numerical_failure &&
                       !tie_rejected.valid() &&
                       face_checksum(as_const(derivative_faces)) == untouched,
                   "limiter active-branch tie has no fabricated Jacobian");
  FrozenConvectionFaceDirectionalDerivative tie_generalized;
  const Status tie_generalized_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::limited_central2, as_const(target_flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::
              semismooth_generalized_zero_slope,
          tied_frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          tie_generalized);
  const double tie_oracle =
      0.5 * (variation.view.unchecked({-1, 0, 0}, 0U) +
             variation.view.unchecked({0, 0, 0}, 0U));
  passed &= expect(
      static_cast<bool>(tie_generalized_status) &&
          tie_generalized.valid() && !tie_generalized.classical_everywhere &&
          tie_generalized.generalized_face_count > 0U &&
          tie_generalized.policy ==
              FrozenConvectionLinearizationPolicy::
                  semismooth_generalized_zero_slope &&
          close(derivative_faces.x.unchecked({0, 0, 0}), tie_oracle,
                tie_oracle),
      "explicit semismooth policy gives limiter tie zero slope derivative");
  fill_faces(derivative_faces, 811.0);

  const std::array<double, 10U> kinked_x{
      0.0, 0.0, 1.0, 3.5, 7.2, 12.1, 18.4, 26.0, 34.9, 45.3};
  for (std::int32_t z = -2; z < kCells.z + 2; ++z) {
    for (std::int32_t y = -2; y < kCells.y + 2; ++y) {
      for (std::int32_t x = -2; x < kCells.x + 2; ++x) {
        target.view.unchecked({x, y, z}, 0U) =
            kinked_x[static_cast<std::size_t>(x + 2)] + 0.125 * y +
            0.0625 * z;
      }
    }
  }
  ++target.view.revision;
  FrozenConvectionFaceField kinked_frozen;
  passed &= expect(
      static_cast<bool>(freeze_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(target_flux),
          as_const(target.view), 0U, context,
          {frozen_faces.x, frozen_faces.y, frozen_faces.z}, kinked_frozen)),
      "finite zero-slope kink still has an exact nonlinear reconstruction");
  FrozenConvectionFaceDirectionalDerivative kink_rejected;
  const Status kink_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(target_flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::classical_active_branch,
          kinked_frozen,
          as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          kink_rejected);
  passed &= expect(kink_status.code == StatusCode::numerical_failure &&
                       !kink_rejected.valid() &&
                       face_checksum(as_const(derivative_faces)) == untouched,
                   "zero-slope minmod kink fails closed before writes");
  FrozenConvectionFaceDirectionalDerivative kink_generalized;
  const Status kink_generalized_status =
      differentiate_frozen_cartesian_target_convection_faces(
          plan, ConvectionScheme::tvd2, as_const(target_flux),
          as_const(target.view), 0U, context,
          FrozenConvectionLinearizationPolicy::
              semismooth_generalized_zero_slope,
          kinked_frozen, as_const(variation.view), 0U,
          {derivative_faces.x, derivative_faces.y, derivative_faces.z},
          kink_generalized);
  const double kink_oracle =
      target_flux.x.unchecked({0, 0, 0}) >= 0.0
          ? variation.view.unchecked({-1, 0, 0}, 0U)
          : variation.view.unchecked({0, 0, 0}, 0U);
  passed &= expect(
      static_cast<bool>(kink_generalized_status) &&
          kink_generalized.valid() &&
          kink_generalized.generalized_face_count > 0U &&
          close(derivative_faces.x.unchecked({0, 0, 0}), kink_oracle,
                kink_oracle),
      "explicit semismooth TVD kink keeps only the frozen donor value");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  KernelFixture fixture;
  KernelFixture zero_limiter_fixture;
  KernelFixture stretched_fixture;
  FaceFluxStorage storage;
  FaceFluxView raw_flux;
  bool passed = make_kernel_fixture(fixture);
  passed &= make_kernel_fixture(zero_limiter_fixture, mesh_spec(), 0.0);
  passed &= make_kernel_fixture(stretched_fixture, stretched_mesh_spec(), 1.0);
  passed &= test_face_layout(storage, raw_flux);
  if (passed) {
    passed &= test_shared_face_and_global_conservation(fixture.kernels,
                                                       raw_flux);
    passed &= test_exact_kernel_counters_and_allocations(fixture.kernels,
                                                         raw_flux);
    passed &= test_alias_mutations(fixture.kernels, raw_flux);
    passed &= test_convection_revision_and_hot_counters(fixture.kernels,
                                                        raw_flux);
    passed &= test_frozen_target_convection_faces(fixture.kernels, storage,
                                                   raw_flux);
    passed &= test_frozen_target_convection_directional_derivative(
        fixture.kernels, as_const(raw_flux));
    passed &= test_frozen_direction_analytic_branches(fixture.kernels);
    passed &= test_constant_semismooth_and_zero_limiter(
        fixture.kernels, zero_limiter_fixture.kernels);
    passed &= test_subnormal_same_sign_limiter(fixture.kernels);
    passed &= test_stretched_direction_oracle(stretched_fixture);
    passed &= test_frozen_target_direction_fail_closed(fixture.kernels,
                                                       raw_flux);
    passed &= test_final_flux_authority(fixture.kernels);
    passed &= test_final_flux_lifetime_fail_closed(fixture.kernels);
    passed &= test_point_convection_diagnostic(fixture.kernels);
  }
  if (passed) {
    std::cout << "v0.4 conservative face-flux tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
