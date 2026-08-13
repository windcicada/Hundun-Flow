// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_parallel.hpp"
#include "parallel_halo_detail.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
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

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  if (void* pointer = std::malloc(size == 0U ? 1U : size)) {
    return pointer;
  }
  throw std::bad_alloc{};
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* pointer = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&pointer, alignment, requested) == 0 &&
      pointer != nullptr) {
    return pointer;
  }
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

void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}

void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}

void* operator new(std::size_t size, std::align_val_t alignment,
                   const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate_aligned(
        size, static_cast<std::size_t>(alignment));
  } catch (...) {
    return nullptr;
  }
}

void* operator new[](std::size_t size, std::align_val_t alignment,
                     const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate_aligned(
        size, static_cast<std::size_t>(alignment));
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

constexpr Int3 kInterior{4, 3, 2};
constexpr Int3 kGhosts{1, 1, 1};
constexpr std::uint8_t kVelocityComponents = 2U;
constexpr FieldId kVelocity = 0U;
constexpr FieldId kPressure = 1U;

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local, MPI_Comm communicator = MPI_COMM_WORLD) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN, communicator) ==
             MPI_SUCCESS &&
         output != 0;
}

std::size_t extent(std::int32_t interior, std::int32_t ghosts) {
  return static_cast<std::size_t>(interior + 2 * ghosts);
}

OwnedField make_field(FieldId id, std::uint8_t components,
                      RevisionToken revision) {
  OwnedField result;
  const std::size_t stride_y = extent(kInterior.x, kGhosts.x);
  const std::size_t stride_z = stride_y * extent(kInterior.y, kGhosts.y);
  const std::size_t component_stride =
      stride_z * extent(kInterior.z, kGhosts.z);
  result.storage.assign(component_stride * components, -777777.0);
  result.view.base = result.storage.data() +
                     static_cast<std::ptrdiff_t>(kGhosts.x) +
                     static_cast<std::ptrdiff_t>(kGhosts.y) *
                         static_cast<std::ptrdiff_t>(stride_y) +
                     static_cast<std::ptrdiff_t>(kGhosts.z) *
                         static_cast<std::ptrdiff_t>(stride_z);
  result.view.interior = kInterior;
  result.view.ghosts = kGhosts;
  result.view.components = components;
  result.view.stride_y = stride_y;
  result.view.stride_z = stride_z;
  result.view.component_stride = component_stride;
  result.view.field = id;
  result.view.revision = revision;
  result.view.storage_identity = static_cast<StorageIdentity>(100U + id);
  result.view.revision_domain = 9001U;
  return result;
}

double value_for(int rank, FieldId field, std::uint8_t component,
                 Int3 index) {
  return static_cast<double>(rank * 1000000 + static_cast<int>(field) * 100000 +
                             static_cast<int>(component) * 10000 +
                             index.z * 100 + index.y * 10 + index.x);
}

void fill_interior(OwnedField& field, int rank) {
  for (std::uint8_t component = 0U; component < field.view.components;
       ++component) {
    for (std::int32_t z = 0; z < kInterior.z; ++z) {
      for (std::int32_t y = 0; y < kInterior.y; ++y) {
        for (std::int32_t x = 0; x < kInterior.x; ++x) {
          const Int3 index{x, y, z};
          field.view.unchecked(index, component) =
              value_for(rank, field.view.field, component, index);
        }
      }
    }
  }
}

void reset_ghosts(OwnedField& field, double sentinel) {
  for (std::uint8_t component = 0U; component < field.view.components;
       ++component) {
    for (std::int32_t z = -1; z <= kInterior.z; ++z) {
      for (std::int32_t y = -1; y <= kInterior.y; ++y) {
        for (std::int32_t x = -1; x <= kInterior.x; ++x) {
          const bool interior = x >= 0 && x < kInterior.x && y >= 0 &&
                                y < kInterior.y && z >= 0 && z < kInterior.z;
          if (!interior) {
            field.view.unchecked(Int3{x, y, z}, component) = sentinel;
          }
        }
      }
    }
  }
}

bool verify_face(const OwnedField& field, CartesianAxis axis, bool high,
                 int source_rank, int rank, std::string_view description) {
  bool passed = true;
  for (std::uint8_t component = 0U; component < field.view.components;
       ++component) {
    const std::int32_t x0 = axis == CartesianAxis::x
                                ? (high ? kInterior.x : -1)
                                : 0;
    const std::int32_t x1 = axis == CartesianAxis::x ? x0 + 1 : kInterior.x;
    const std::int32_t y0 = axis == CartesianAxis::y
                                ? (high ? kInterior.y : -1)
                                : 0;
    const std::int32_t y1 = axis == CartesianAxis::y ? y0 + 1 : kInterior.y;
    const std::int32_t z0 = axis == CartesianAxis::z
                                ? (high ? kInterior.z : -1)
                                : 0;
    const std::int32_t z1 = axis == CartesianAxis::z ? z0 + 1 : kInterior.z;
    for (std::int32_t z = z0; z < z1; ++z) {
      for (std::int32_t y = y0; y < y1; ++y) {
        for (std::int32_t x = x0; x < x1; ++x) {
          Int3 source{x, y, z};
          if (axis == CartesianAxis::x) {
            source.x = high ? 0 : kInterior.x - 1;
          } else if (axis == CartesianAxis::y) {
            source.y = high ? 0 : kInterior.y - 1;
          } else {
            source.z = high ? 0 : kInterior.z - 1;
          }
          passed &= expect(field.view.unchecked(Int3{x, y, z}, component) ==
                               value_for(source_rank, field.view.field,
                                         component, source),
                           rank, description);
        }
      }
    }
  }
  return passed;
}

bool verify_untouched_face(const OwnedField& field, CartesianAxis axis,
                           bool high, double sentinel, int rank,
                           std::string_view description) {
  bool passed = true;
  for (std::uint8_t component = 0U; component < field.view.components;
       ++component) {
    const std::int32_t x0 = axis == CartesianAxis::x
                                ? (high ? kInterior.x : -1)
                                : 0;
    const std::int32_t x1 = axis == CartesianAxis::x ? x0 + 1 : kInterior.x;
    const std::int32_t y0 = axis == CartesianAxis::y
                                ? (high ? kInterior.y : -1)
                                : 0;
    const std::int32_t y1 = axis == CartesianAxis::y ? y0 + 1 : kInterior.y;
    const std::int32_t z0 = axis == CartesianAxis::z
                                ? (high ? kInterior.z : -1)
                                : 0;
    const std::int32_t z1 = axis == CartesianAxis::z ? z0 + 1 : kInterior.z;
    for (std::int32_t z = z0; z < z1; ++z) {
      for (std::int32_t y = y0; y < y1; ++y) {
        for (std::int32_t x = x0; x < x1; ++x) {
          passed &= expect(field.view.unchecked(Int3{x, y, z}, component) ==
                               sentinel,
                           rank, description);
        }
      }
    }
  }
  return passed;
}

MeshPatch patch_for(int rank, int size) {
  return MeshPatch{Int3{rank * kInterior.x, 0, 0}, kInterior,
                   Int3{size, 1, 1}, Int3{rank, 0, 0}};
}

std::array<HaloFieldSpec, 2U> field_specs() {
  return {HaloFieldSpec{kVelocity, 1U, kVelocityComponents},
          HaloFieldSpec{kPressure, 1U, 1U}};
}

std::array<FieldView, 2U> field_views(OwnedField& velocity,
                                     OwnedField& pressure) {
  return {pressure.view, velocity.view};  // Deliberately not specification order.
}

bool exchange(HaloEngine& engine, StageId stage,
              std::array<FieldView, 2U>& views, int rank,
              bool check_publication_boundary) {
  HaloTicket ticket;
  bool passed = expect(static_cast<bool>(engine.begin(
                           stage, Span<const FieldView>{views.data(), views.size()},
                           ticket)),
                       rank, "halo begin succeeds");
  passed &= expect(ticket.active() && ticket.stage() == stage, rank,
                   "begin publishes one active ticket for the requested stage");
  if (check_publication_boundary) {
    passed &= expect(engine.ghost_revision(kVelocity) == 0U &&
                         engine.ghost_revision(kPressure) == 0U,
                     rank, "begin does not publish a ghost certificate");
  }
  passed &= expect(static_cast<bool>(engine.finish(
                       ticket, Span<FieldView>{views.data(), views.size()})),
                   rank, "halo finish succeeds");
  passed &= expect(!ticket.active() &&
                       engine.ghost_revision(kVelocity) == 17U &&
                       engine.ghost_revision(kPressure) == 23U,
                   rank, "finish alone publishes exact field revisions");
  return passed;
}

bool test_begin_atomicity(int rank, int size) {
  OwnedField velocity = make_field(kVelocity, kVelocityComponents, 17U);
  OwnedField pressure = make_field(kPressure, 1U, 23U);
  fill_interior(velocity, rank);
  fill_interior(pressure, rank);
  reset_ghosts(velocity, -737373.0);
  reset_ghosts(pressure, -737373.0);
  const auto specs = field_specs();
  auto views = field_views(velocity, pressure);
  HaloEngine engine;
  bool passed = expect(
      static_cast<bool>(engine.reserve(
          MPI_COMM_WORLD, patch_for(rank, size),
          Span<const HaloFieldSpec>{specs.data(), specs.size()},
          HaloTopology{true, false, false})),
      rank, "atomicity halo plan reserves");

  HaloTicket first;
  passed &= expect(
      static_cast<bool>(engine.begin(
          61U, Span<const FieldView>{views.data(), views.size()}, first)) &&
          first.active(),
      rank, "first exchange begins before re-entrant attempt");
  HaloTicket second;
  const Status reentrant = engine.begin(
      62U, Span<const FieldView>{views.data(), views.size()}, second);
  passed &= expect(reentrant.code == StatusCode::invalid_plan &&
                       !second.active() && first.active(),
                   rank, "re-entrant begin is collectively rejected");
  passed &= expect(
      static_cast<bool>(engine.finish(
          first, Span<FieldView>{views.data(), views.size()})) &&
          engine.ghost_revision(kVelocity) == 17U &&
          engine.ghost_revision(kPressure) == 23U,
      rank, "rejected begin preserves first exchange revision metadata");
  const int low_rank = rank == 0 ? size - 1 : rank - 1;
  passed &= verify_face(velocity, CartesianAxis::x, false, low_rank, rank,
                        "first exchange remains exact after rejected begin");

  std::array<FieldView, 2U> malformed = views;
  malformed[1].components = 1U;
  HaloTicket malformed_ticket;
  const Status malformed_result = engine.begin(
      63U, Span<const FieldView>{malformed.data(), malformed.size()},
      malformed_ticket);
  passed &= expect(malformed_result.code == StatusCode::invalid_plan &&
                       !malformed_ticket.active() &&
                       engine.ghost_revision(kVelocity) == 0U &&
                       engine.ghost_revision(kPressure) == 0U,
                   rank,
                   "malformed field set is rejected without a certificate");
  passed &= exchange(engine, 64U, views, rank, true);
  passed &= verify_face(velocity, CartesianAxis::x, false, low_rank, rank,
                        "valid exchange recovers after partial validation");
  return all_true(passed);
}

bool test_forced_chunking(int rank, int size) {
  OwnedField velocity = make_field(kVelocity, kVelocityComponents, 17U);
  OwnedField pressure = make_field(kPressure, 1U, 23U);
  fill_interior(velocity, rank);
  fill_interior(pressure, rank);
  reset_ghosts(velocity, -626262.0);
  reset_ghosts(pressure, -626262.0);
  const auto specs = field_specs();
  auto views = field_views(velocity, pressure);
  detail::set_halo_maximum_chunk_doubles_for_test(5U);
  HaloEngine engine;
  const Status reserved = engine.reserve(
      MPI_COMM_WORLD, patch_for(rank, size),
      Span<const HaloFieldSpec>{specs.data(), specs.size()},
      HaloTopology{true, true, true});
  detail::clear_halo_maximum_chunk_doubles_for_test();
  bool passed = expect(static_cast<bool>(reserved), rank,
                       "forced multi-chunk halo plan reserves");
  const HaloPlanStats stats = engine.plan_stats();
  if (size > 1) {
    passed &= expect(stats.maximum_messages_per_exchange >
                         stats.transport_peer_count,
                     rank,
                     "small test chunks produce multiple tagged messages per peer");
  }
  passed &= exchange(engine, 65U, views, rank, true);
  const int low_rank = rank == 0 ? size - 1 : rank - 1;
  const int high_rank = rank + 1 == size ? 0 : rank + 1;
  passed &= verify_face(velocity, CartesianAxis::x, false, low_rank, rank,
                        "forced chunks preserve x-minus payload ordering");
  passed &= verify_face(pressure, CartesianAxis::x, true, high_rank, rank,
                        "forced chunks preserve x-plus payload ordering");
  return all_true(passed);
}

bool test_tag_upper_bound_normalization(int rank, int size) {
  const auto specs = field_specs();
  const auto reserve_with_bound = [&](int injected_bound,
                                      int expected_bound,
                                      std::string_view description) {
    detail::set_halo_tag_upper_bound_for_test(injected_bound);
    HaloEngine engine;
    const Status reserved = engine.reserve(
        MPI_COMM_WORLD, patch_for(rank, size),
        Span<const HaloFieldSpec>{specs.data(), specs.size()},
        HaloTopology{true, false, false});
    detail::clear_halo_tag_upper_bound_for_test();
    return expect(static_cast<bool>(reserved) && engine.ready() &&
                      engine.plan_stats().maximum_tag == expected_bound,
                  rank, description);
  };

  bool passed = reserve_with_bound(
      -1, 32767,
      "missing MPI_TAG_UB attribute falls back to the guaranteed bound");
  const int different_bound =
      size == 1 ? 41001 : (rank == size - 1 ? 35001 : 41001);
  passed &= reserve_with_bound(
      different_bound, size == 1 ? 41001 : 35001,
      "rank-local MPI_TAG_UB values use one collective minimum");
  const int anomalous_bound = rank == size - 1 ? 17 : 41001;
  passed &= reserve_with_bound(
      anomalous_bound, 17,
      "anomalously low MPI_TAG_UB values use one conservative collective minimum");
  return all_true(passed);
}

bool test_contract_query_is_exact_local_and_read_only(int rank, int size) {
  OwnedField velocity = make_field(kVelocity, kVelocityComponents, 17U);
  OwnedField pressure = make_field(kPressure, 1U, 23U);
  fill_interior(velocity, rank);
  fill_interior(pressure, rank);
  reset_ghosts(velocity, -525252.0);
  reset_ghosts(pressure, -525252.0);
  const auto specs = field_specs();
  constexpr HaloTopology topology{true, false, false};
  HaloEngine engine;
  bool passed = expect(
      static_cast<bool>(engine.reserve(
          MPI_COMM_WORLD, patch_for(rank, size),
          Span<const HaloFieldSpec>{specs.data(), specs.size()},
          topology)),
      rank, "contract-query halo plan reserves");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  MPI_Comm congruent = MPI_COMM_NULL;
  MPI_Comm reversed = MPI_COMM_NULL;
  passed &= expect(MPI_Comm_dup(MPI_COMM_WORLD, &congruent) == MPI_SUCCESS &&
                       MPI_Comm_split(MPI_COMM_WORLD, 0, size - 1 - rank,
                                      &reversed) == MPI_SUCCESS,
                   rank, "halo comparison communicators compile");
  if (!all_true(passed)) {
    if (congruent != MPI_COMM_NULL) (void)MPI_Comm_free(&congruent);
    if (reversed != MPI_COMM_NULL) (void)MPI_Comm_free(&reversed);
    return false;
  }

  const HaloPlanStats plan_before = engine.plan_stats();
  const HaloRuntimeCounters counters_before = engine.runtime_counters();
  MeshPatch wrong_patch = patch_for(rank, size);
  ++wrong_patch.begin.x;
  auto wrong_fields = specs;
  ++wrong_fields[0].components;
  const std::array<HaloFieldSpec, 2U> reversed_fields{
      specs[1], specs[0]};
  const Status identical = engine.validate_contract(
      MPI_COMM_WORLD, patch_for(rank, size),
      Span<const HaloFieldSpec>{specs.data(), specs.size()}, topology);
  const Status compatible = engine.validate_contract(
      congruent, patch_for(rank, size),
      Span<const HaloFieldSpec>{specs.data(), specs.size()}, topology);
  const Status similar = engine.validate_contract(
      reversed, patch_for(rank, size),
      Span<const HaloFieldSpec>{specs.data(), specs.size()}, topology);
  const Status patch_mismatch = engine.validate_contract(
      MPI_COMM_WORLD, wrong_patch,
      Span<const HaloFieldSpec>{specs.data(), specs.size()}, topology);
  const Status field_mismatch = engine.validate_contract(
      MPI_COMM_WORLD, patch_for(rank, size),
      Span<const HaloFieldSpec>{wrong_fields.data(), wrong_fields.size()},
      topology);
  const Status order_mismatch = engine.validate_contract(
      MPI_COMM_WORLD, patch_for(rank, size),
      Span<const HaloFieldSpec>{reversed_fields.data(),
                                reversed_fields.size()},
      topology);
  const Status null = engine.validate_contract(
      MPI_COMM_NULL, patch_for(rank, size),
      Span<const HaloFieldSpec>{specs.data(), specs.size()}, topology);
  const HaloPlanStats plan_after = engine.plan_stats();
  const HaloRuntimeCounters counters_after = engine.runtime_counters();
  passed &= expect(static_cast<bool>(identical) &&
                       static_cast<bool>(compatible) &&
                       (size == 1 ? static_cast<bool>(similar)
                                  : similar.code == StatusCode::invalid_plan) &&
                       patch_mismatch.code == StatusCode::invalid_plan &&
                       field_mismatch.code == StatusCode::invalid_plan &&
                       order_mismatch.code == StatusCode::invalid_plan &&
                       null.code == StatusCode::invalid_plan &&
                       plan_after.request_storage_address ==
                           plan_before.request_storage_address &&
                       plan_after.send_storage_address ==
                           plan_before.send_storage_address &&
                       plan_after.receive_storage_address ==
                           plan_before.receive_storage_address &&
                       counters_after.begin_calls == counters_before.begin_calls &&
                       counters_after.finish_calls == counters_before.finish_calls &&
                       !engine.active(),
                   rank,
                   "canonical fields and CONGRUENT pass; SIMILAR/wrong/order contracts reject without mutation");

  auto views = field_views(velocity, pressure);
  passed &= exchange(engine, 66U, views, rank, true);
  (void)MPI_Comm_free(&reversed);
  (void)MPI_Comm_free(&congruent);
  return all_true(passed);
}

bool test_contract_query_rejects_bidirectional_topology_mutation(int rank,
                                                                 int size) {
  const auto specs = field_specs();
  const MeshPatch patch = patch_for(rank, size);
  constexpr HaloTopology forward{true, false, true};
  constexpr HaloTopology reverse{false, true, false};
  HaloEngine forward_engine;
  HaloEngine reverse_engine;
  bool passed = expect(
      static_cast<bool>(forward_engine.reserve(
          MPI_COMM_WORLD, patch,
          Span<const HaloFieldSpec>{specs.data(), specs.size()}, forward)) &&
          static_cast<bool>(reverse_engine.reserve(
              MPI_COMM_WORLD, patch,
              Span<const HaloFieldSpec>{specs.data(), specs.size()}, reverse)),
      rank, "opposite halo topology plans reserve");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const auto validate = [&](const HaloEngine& engine,
                            HaloTopology topology) noexcept {
    return engine.validate_contract(
        MPI_COMM_WORLD, patch,
        Span<const HaloFieldSpec>{specs.data(), specs.size()}, topology);
  };
  const Status forward_exact = validate(forward_engine, forward);
  const Status reverse_exact = validate(reverse_engine, reverse);
  const Status forward_x_mutation =
      validate(forward_engine, HaloTopology{false, false, true});
  const Status forward_y_mutation =
      validate(forward_engine, HaloTopology{true, true, true});
  const Status forward_z_mutation =
      validate(forward_engine, HaloTopology{true, false, false});
  const Status reverse_x_mutation =
      validate(reverse_engine, HaloTopology{true, true, false});
  const Status reverse_y_mutation =
      validate(reverse_engine, HaloTopology{false, false, false});
  const Status reverse_z_mutation =
      validate(reverse_engine, HaloTopology{false, true, true});

  const auto rejects = [](Status status) noexcept {
    return status.code == StatusCode::invalid_plan;
  };
  passed &= expect(
      static_cast<bool>(forward_exact) && static_cast<bool>(reverse_exact) &&
          rejects(forward_x_mutation) && rejects(forward_y_mutation) &&
          rejects(forward_z_mutation) && rejects(reverse_x_mutation) &&
          rejects(reverse_y_mutation) && rejects(reverse_z_mutation),
      rank,
      "halo contract compares every periodic axis in both mutation directions");
  return all_true(passed);
}

bool test_nonperiodic(int rank, int size) {
  OwnedField velocity = make_field(kVelocity, kVelocityComponents, 17U);
  OwnedField pressure = make_field(kPressure, 1U, 23U);
  fill_interior(velocity, rank);
  fill_interior(pressure, rank);
  constexpr double sentinel = -818181.0;
  reset_ghosts(velocity, sentinel);
  reset_ghosts(pressure, sentinel);

  const auto specs = field_specs();
  HaloEngine engine;
  bool passed = expect(static_cast<bool>(engine.reserve(
                           MPI_COMM_WORLD, patch_for(rank, size),
                           Span<const HaloFieldSpec>{specs.data(), specs.size()})),
                       rank, "nonperiodic halo plan reserves");
  const HaloPlanStats before = engine.plan_stats();
  passed &= expect(before.local_peer_count == 0U, rank,
                   "nonperiodic decomposition has no self-periodic peer");
  const std::size_t expected_peers =
      size == 1 ? 0U : static_cast<std::size_t>((rank > 0 ? 1 : 0) +
                                                (rank + 1 < size ? 1 : 0));
  passed &= expect(before.transport_peer_count == expected_peers, rank,
                   "nonperiodic plan contains only direct face ranks");
  passed &= expect(before.persistent_request_count ==
                           2U * before.maximum_messages_per_exchange,
                   rank, "each outgoing persistent chunk has one matched receive");

  auto views = field_views(velocity, pressure);
  passed &= exchange(engine, 31U, views, rank, true);
  if (rank > 0) {
    passed &= verify_face(velocity, CartesianAxis::x, false, rank - 1, rank,
                          "x-minus vector ghost equals neighbor interior");
    passed &= verify_face(pressure, CartesianAxis::x, false, rank - 1, rank,
                          "x-minus scalar ghost equals neighbor interior");
  } else {
    passed &= verify_untouched_face(velocity, CartesianAxis::x, false, sentinel,
                                    rank, "physical x-minus face stays untouched");
  }
  if (rank + 1 < size) {
    passed &= verify_face(velocity, CartesianAxis::x, true, rank + 1, rank,
                          "x-plus vector ghost equals neighbor interior");
    passed &= verify_face(pressure, CartesianAxis::x, true, rank + 1, rank,
                          "x-plus scalar ghost equals neighbor interior");
  } else {
    passed &= verify_untouched_face(pressure, CartesianAxis::x, true, sentinel,
                                    rank, "physical x-plus face stays untouched");
  }
  passed &= verify_untouched_face(velocity, CartesianAxis::y, false, sentinel,
                                  rank, "nonperiodic y-minus stays untouched");
  passed &= verify_untouched_face(pressure, CartesianAxis::z, true, sentinel,
                                  rank, "nonperiodic z-plus stays untouched");

  const HaloPlanStats after = engine.plan_stats();
  const HaloRuntimeCounters counters = engine.runtime_counters();
  passed &= expect(before.request_storage_address == after.request_storage_address &&
                       before.send_storage_address == after.send_storage_address &&
                       before.receive_storage_address ==
                           after.receive_storage_address,
                   rank, "persistent request and buffer addresses remain stable");
  passed &= expect(counters.begin_calls == 1U && counters.finish_calls == 1U &&
                       counters.bytes_packed ==
                           before.send_capacity_doubles * sizeof(double) &&
                       counters.bytes_unpacked ==
                           before.receive_capacity_doubles * sizeof(double) &&
                       counters.messages_started ==
                           before.maximum_messages_per_exchange &&
                       before.maximum_bytes_per_exchange ==
                           before.send_capacity_doubles * sizeof(double),
                   rank, "one exchange exactly matches compiled byte/message bounds");
  return all_true(passed);
}

bool test_periodic_and_hot_path(int rank, int size) {
  OwnedField velocity = make_field(kVelocity, kVelocityComponents, 17U);
  OwnedField pressure = make_field(kPressure, 1U, 23U);
  fill_interior(velocity, rank);
  fill_interior(pressure, rank);
  reset_ghosts(velocity, -919191.0);
  reset_ghosts(pressure, -919191.0);

  const auto specs = field_specs();
  HaloEngine engine;
  const HaloTopology periodic{true, true, true};
  bool passed = expect(static_cast<bool>(engine.reserve(
                           MPI_COMM_WORLD, patch_for(rank, size),
                           Span<const HaloFieldSpec>{specs.data(), specs.size()},
                           periodic)),
                       rank, "fully periodic halo plan reserves");
  const HaloPlanStats compiled = engine.plan_stats();
  passed &= expect(compiled.persistent_request_count ==
                           2U * compiled.maximum_messages_per_exchange &&
                       compiled.maximum_bytes_per_exchange ==
                           compiled.send_capacity_doubles * sizeof(double),
                   rank, "periodic persistent requests and byte capacity agree");
  const std::size_t expected_transport_peers =
      size == 1 ? 0U : (size == 2 ? 1U : 2U);
  passed &= expect(compiled.transport_peer_count == expected_transport_peers &&
                       compiled.local_peer_count == 1U &&
                       compiled.maximum_messages_per_exchange ==
                           expected_transport_peers,
                   rank,
                   "periodic directions and fields merge into one small message per rank peer");

  auto views = field_views(velocity, pressure);
  passed &= exchange(engine, 41U, views, rank, true);
  const int low_rank = rank == 0 ? size - 1 : rank - 1;
  const int high_rank = rank + 1 == size ? 0 : rank + 1;
  passed &= verify_face(velocity, CartesianAxis::x, false, low_rank, rank,
                        "periodic x-minus wraps exactly");
  passed &= verify_face(velocity, CartesianAxis::x, true, high_rank, rank,
                        "periodic x-plus wraps exactly");
  passed &= verify_face(pressure, CartesianAxis::y, false, rank, rank,
                        "self-periodic y-minus copies opposite interior");
  passed &= verify_face(pressure, CartesianAxis::y, true, rank, rank,
                        "self-periodic y-plus copies opposite interior");
  passed &= verify_face(velocity, CartesianAxis::z, false, rank, rank,
                        "self-periodic z-minus copies opposite interior");
  passed &= verify_face(velocity, CartesianAxis::z, true, rank, rank,
                        "self-periodic z-plus copies opposite interior");

  constexpr std::size_t repetitions = 32U;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  const HaloRuntimeCounters counters_before = engine.runtime_counters();
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0U; iteration < repetitions; ++iteration) {
      HaloTicket ticket;
      const Status begin = engine.begin(
          static_cast<StageId>(100U + iteration),
          Span<const FieldView>{views.data(), views.size()}, ticket);
      const Status finish = engine.finish(
          ticket, Span<FieldView>{views.data(), views.size()});
      passed &= static_cast<bool>(begin) && static_cast<bool>(finish);
    }
    hot_allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  const HaloRuntimeCounters counters_after = engine.runtime_counters();
  const HaloPlanStats after = engine.plan_stats();
  passed &= expect(hot_allocations == 0U, rank,
                   "repeated begin/finish performs zero C++ heap allocations");
  passed &= expect(after.request_storage_address ==
                           compiled.request_storage_address &&
                       after.send_storage_address == compiled.send_storage_address &&
                       after.receive_storage_address ==
                           compiled.receive_storage_address,
                   rank, "hot repetition preserves all persistent addresses");
  passed &= expect(counters_after.begin_calls - counters_before.begin_calls ==
                           repetitions &&
                       counters_after.finish_calls - counters_before.finish_calls ==
                           repetitions,
                   rank, "hot repetition increments begin/finish exactly once");
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  const auto uninitialized_specs = field_specs();
  const Status uninitialized = HaloEngine{}.validate_contract(
      MPI_COMM_WORLD, MeshPatch{},
      Span<const HaloFieldSpec>{uninitialized_specs.data(),
                                uninitialized_specs.size()},
      HaloTopology{});
  if (uninitialized.code != StatusCode::invalid_plan) {
    return 3;
  }
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  bool passed = test_nonperiodic(rank, size);
  passed &= test_periodic_and_hot_path(rank, size);
  passed &= test_begin_atomicity(rank, size);
  passed &= test_forced_chunking(rank, size);
  passed &= test_tag_upper_bound_normalization(rank, size);
  passed &= test_contract_query_is_exact_local_and_read_only(rank, size);
  passed &=
      test_contract_query_rejects_bidirectional_topology_mutation(rank, size);

  const auto finalized_specs = field_specs();
  const MeshPatch finalized_patch = patch_for(rank, size);
  constexpr HaloTopology finalized_topology{true, false, false};
  HaloEngine finalized_engine;
  passed &= static_cast<bool>(finalized_engine.reserve(
      MPI_COMM_WORLD, finalized_patch,
      Span<const HaloFieldSpec>{finalized_specs.data(), finalized_specs.size()},
      finalized_topology));
  passed = all_true(passed);
  if (rank == 0 && !passed) {
    std::cerr << "v0.4 halo MPI test failed for " << size << " ranks\n";
  }
  MPI_Finalize();
  const Status finalized = finalized_engine.validate_contract(
      MPI_COMM_WORLD, finalized_patch,
      Span<const HaloFieldSpec>{finalized_specs.data(), finalized_specs.size()},
      finalized_topology);
  passed &= finalized.code == StatusCode::invalid_plan;
  if (rank == 0 && finalized.code != StatusCode::invalid_plan) {
    std::cerr << "v0.4 halo finalized contract query did not reject\n";
  }
  return passed ? 0 : 1;
}
