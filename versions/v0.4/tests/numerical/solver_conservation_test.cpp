// SPDX-License-Identifier: Apache-2.0

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

ValidatedModel periodic_model() {
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
  value.schemes.limiter = 1.0;
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

bool make_kernel_fixture(KernelFixture& fixture) {
  FieldRegistry registry;
  TimeSchemePlan time;
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh_spec(), GeometryBudget{}, fixture.geometry,
          fixture.patch)),
      "uniform Cartesian geometry compiles");
  passed &= expect(
      static_cast<bool>(BoundaryCompiler::compile(
          MPI_COMM_SELF, periodic_model(), fixture.geometry, fixture.patch,
          registry, fixture.boundary, fixture.schemes, time)),
      "periodic boundary and central schemes compile");
  passed &= expect(
      static_cast<bool>(CartesianKernelPlan::compile(
          fixture.schemes, fixture.geometry, fixture.patch, fixture.boundary,
          fixture.kernels)),
      "Cartesian kernel plan compiles once outside the hot path");
  return passed;
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

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  KernelFixture fixture;
  FaceFluxStorage storage;
  FaceFluxView raw_flux;
  bool passed = make_kernel_fixture(fixture);
  passed &= test_face_layout(storage, raw_flux);
  if (passed) {
    passed &= test_shared_face_and_global_conservation(fixture.kernels,
                                                       raw_flux);
    passed &= test_exact_kernel_counters_and_allocations(fixture.kernels,
                                                         raw_flux);
    passed &= test_alias_mutations(fixture.kernels, raw_flux);
    passed &= test_convection_revision_and_hot_counters(fixture.kernels,
                                                        raw_flux);
    passed &= test_final_flux_authority(fixture.kernels);
    passed &= test_final_flux_lifetime_fail_closed(fixture.kernels);
  }
  if (passed) {
    std::cout << "v0.4 conservative face-flux tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}
