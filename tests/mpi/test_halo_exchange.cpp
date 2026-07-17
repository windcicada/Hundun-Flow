// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_descriptor.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "runtime/src/halo_detail.hpp"
#include "runtime/src/halo_test_access.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

namespace {

using hundun::runtime::Box3;
using hundun::runtime::Error;
using hundun::runtime::ExchangePlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::HaloExchange;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;

static_assert(std::is_copy_constructible_v<ExchangePlan>);
static_assert(!std::is_copy_constructible_v<HaloExchange>);
static_assert(!std::is_copy_assignable_v<HaloExchange>);
static_assert(std::is_nothrow_move_constructible_v<HaloExchange>);
static_assert(!std::is_move_assignable_v<HaloExchange>);

constexpr Int3 kGlobal{12, 10, 8};
constexpr std::array<bool, 3> kPeriodic{true, true, true};

bool same(Int3 left, Int3 right) {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(Box3 left, Box3 right) {
  return same(left.begin, right.begin) && same(left.end, right.end);
}

std::size_t box_cells(Box3 box) {
  const auto nx = static_cast<std::size_t>(box.end.x - box.begin.x);
  const auto ny = static_cast<std::size_t>(box.end.y - box.begin.y);
  const auto nz = static_cast<std::size_t>(box.end.z - box.begin.z);
  return nx * ny * nz;
}

FieldDescriptor descriptor(std::string name, ScalarType scalar_type,
                           std::uint32_t components, int ghost_width) {
  return FieldDescriptor{std::move(name),
                         "1",
                         "halo_test",
                         FunctionSpace::cell_average,
                         scalar_type,
                         components,
                         ghost_width,
                         false,
                         RestartPolicy::transient,
                         OutputPolicy::never};
}

template <class Function>
std::string expect_collective_error(const MpiContext& context,
                                    Function&& function) {
  bool threw = false;
  std::string message;
  try {
    function();
  } catch (const Error& error) {
    threw = true;
    message = error.what();
  }
  const int local_threw = threw ? 1 : 0;
  int total_threw = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_threw, &total_threw, 1, MPI_INT, MPI_SUM,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(total_threw == context.size());
  HUNDUN_CHECK(!message.empty());
  HUNDUN_CHECK(message.size() < 256U);
  std::array<char, 256> reference{};
  if (context.rank() == 0) {
    std::copy(message.begin(), message.end(), reference.begin());
  }
  HUNDUN_CHECK(MPI_Bcast(reference.data(), static_cast<int>(reference.size()),
                         MPI_CHAR, 0, context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(message == std::string(reference.data()));
  return message;
}

template <class Function>
void expect_local_error(Function&& function) {
  bool threw = false;
  try {
    function();
  } catch (const Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()).empty() == false);
  }
  HUNDUN_CHECK(threw);
}

enum class ExpectedHook { post, wait, cleanup_wait, plan_collective,
                          wire_collective };

unsigned long long collective_count(const MpiContext& context,
                                    std::size_t local_count) {
  const auto local = static_cast<unsigned long long>(local_count);
  unsigned long long total = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local, &total, 1, MPI_UNSIGNED_LONG_LONG,
                             MPI_SUM, context.comm()) == MPI_SUCCESS);
  return total;
}

void check_only_hook_fired(
    const MpiContext& context,
    const hundun::runtime::detail::HaloTestSnapshot& snapshot,
    ExpectedHook expected) {
  HUNDUN_CHECK(collective_count(context, snapshot.post_errors_injected) ==
               (expected == ExpectedHook::post ? 1U : 0U));
  HUNDUN_CHECK(collective_count(context, snapshot.wait_errors_injected) ==
               (expected == ExpectedHook::wait ? 1U : 0U));
  HUNDUN_CHECK(
      collective_count(context, snapshot.cleanup_wait_errors_injected) ==
      (expected == ExpectedHook::cleanup_wait ? 1U : 0U));
  HUNDUN_CHECK(
      collective_count(
          context, snapshot.plan_width_first_collective_errors_injected) ==
      (expected == ExpectedHook::plan_collective ? 1U : 0U));
  HUNDUN_CHECK(
      collective_count(context,
                       snapshot.wire_first_collective_errors_injected) ==
      (expected == ExpectedHook::wire_collective ? 1U : 0U));
}

int wrapped(int value, int extent) {
  int result = value % extent;
  return result < 0 ? result + extent : result;
}

bool oracle_global_cell(const StructuredDecomposition& decomposition,
                        Int3 local_cell, std::array<bool, 3> periodic,
                        Int3& global_cell) {
  const Box3 owned = decomposition.owned_box();
  const Int3 global = decomposition.global_extent();
  global_cell = Int3{owned.begin.x + local_cell.x,
                     owned.begin.y + local_cell.y,
                     owned.begin.z + local_cell.z};
  int* coordinates[3]{&global_cell.x, &global_cell.y, &global_cell.z};
  const int extents[3]{global.x, global.y, global.z};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (*coordinates[axis] < 0 || *coordinates[axis] >= extents[axis]) {
      if (!periodic[axis]) {
        return false;
      }
      *coordinates[axis] = wrapped(*coordinates[axis], extents[axis]);
    }
  }
  return true;
}

std::uint64_t oracle_global_id(Int3 global_cell, Int3 global_extent) {
  return ((static_cast<std::uint64_t>(global_cell.z) *
           static_cast<std::uint64_t>(global_extent.y)) +
          static_cast<std::uint64_t>(global_cell.y)) *
             static_cast<std::uint64_t>(global_extent.x) +
         static_cast<std::uint64_t>(global_cell.x);
}

template <class T>
T sample_value(std::uint64_t global_id, std::uint32_t component,
               int pattern);

template <>
double sample_value<double>(std::uint64_t global_id,
                            std::uint32_t component, int pattern) {
  return static_cast<double>(global_id) +
         static_cast<double>(component) * 2000.0 +
         static_cast<double>(pattern) * 10000.0;
}

template <>
std::int32_t sample_value<std::int32_t>(std::uint64_t global_id,
                                       std::uint32_t component, int pattern) {
  return static_cast<std::int32_t>(global_id) +
         static_cast<std::int32_t>(component) * 2000 + pattern * 10000;
}

template <>
std::uint8_t sample_value<std::uint8_t>(std::uint64_t global_id,
                                       std::uint32_t component, int pattern) {
  return static_cast<std::uint8_t>(
      (global_id + static_cast<std::uint64_t>(component) * 53U +
       static_cast<std::uint64_t>(pattern) * 29U) % 251U);
}

template <class T>
T sentinel_value();

template <>
double sentinel_value<double>() {
  return -9001.0;
}

template <>
std::int32_t sentinel_value<std::int32_t>() {
  return -9001;
}

template <>
std::uint8_t sentinel_value<std::uint8_t>() {
  return static_cast<std::uint8_t>(255);
}

bool owned_coordinate(Int3 cell, Int3 extent) {
  return cell.x >= 0 && cell.x < extent.x && cell.y >= 0 &&
         cell.y < extent.y && cell.z >= 0 && cell.z < extent.z;
}

bool within_exchange_padding(Int3 cell, Int3 extent, int width) {
  return cell.x >= -width && cell.x < extent.x + width &&
         cell.y >= -width && cell.y < extent.y + width &&
         cell.z >= -width && cell.z < extent.z + width;
}

template <class T>
void initialize_field(FieldStorage& storage, FieldId id,
                      const StructuredDecomposition& decomposition,
                      int pattern) {
  auto view = storage.view<T>(id);
  const Int3 extent = view.interior_extent();
  const int ghost = view.ghost_width();
  for (int k = -ghost; k < extent.z + ghost; ++k) {
    for (int j = -ghost; j < extent.y + ghost; ++j) {
      for (int i = -ghost; i < extent.x + ghost; ++i) {
        for (std::uint32_t component = 0; component < view.components();
             ++component) {
          if (owned_coordinate(Int3{i, j, k}, extent)) {
            Int3 global_cell{};
            HUNDUN_CHECK(oracle_global_cell(decomposition, Int3{i, j, k},
                                            kPeriodic, global_cell));
            view(i, j, k, static_cast<int>(component)) = sample_value<T>(
                oracle_global_id(global_cell, decomposition.global_extent()),
                component, pattern);
          } else {
            view(i, j, k, static_cast<int>(component)) = sentinel_value<T>();
          }
        }
      }
    }
  }
}

template <class T>
void set_owned_pattern(FieldStorage& storage, FieldId id,
                       const StructuredDecomposition& decomposition,
                       int pattern) {
  auto view = storage.view<T>(id);
  const Int3 extent = view.interior_extent();
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        Int3 global_cell{};
        HUNDUN_CHECK(oracle_global_cell(decomposition, Int3{i, j, k},
                                        kPeriodic, global_cell));
        for (std::uint32_t component = 0; component < view.components();
             ++component) {
          view(i, j, k, static_cast<int>(component)) = sample_value<T>(
              oracle_global_id(global_cell, decomposition.global_extent()),
              component, pattern);
        }
      }
    }
  }
}

template <class T>
void check_field(const FieldStorage& storage, FieldId id,
                 const StructuredDecomposition& decomposition,
                 std::array<bool, 3> periodic, int exchange_width,
                 int owned_pattern, int ghost_pattern) {
  const auto view = storage.view<T>(id);
  const Int3 extent = view.interior_extent();
  const int ghost = view.ghost_width();
  for (int k = -ghost; k < extent.z + ghost; ++k) {
    for (int j = -ghost; j < extent.y + ghost; ++j) {
      for (int i = -ghost; i < extent.x + ghost; ++i) {
        const Int3 local{i, j, k};
        const bool owned = owned_coordinate(local, extent);
        Int3 global_cell{};
        const bool has_source = oracle_global_cell(
            decomposition, local, periodic, global_cell);
        const bool exchanged = !owned && has_source &&
                               within_exchange_padding(local, extent,
                                                       exchange_width);
        for (std::uint32_t component = 0; component < view.components();
             ++component) {
          T expected = sentinel_value<T>();
          if (owned || exchanged) {
            expected = sample_value<T>(
                oracle_global_id(global_cell, decomposition.global_extent()),
                component, owned ? owned_pattern : ghost_pattern);
          }
          HUNDUN_CHECK(view(i, j, k, static_cast<int>(component)) == expected);
        }
      }
    }
  }
}

Box3 expected_axis_boxes(Int3 extent, int ghost, Int3 offset,
                         bool receive) {
  const auto axis = [ghost, receive](int n, int o) {
    if (o < 0) {
      return receive ? std::array<int, 2>{-ghost, 0}
                     : std::array<int, 2>{0, ghost};
    }
    if (o == 0) {
      return std::array<int, 2>{0, n};
    }
    return receive ? std::array<int, 2>{n, n + ghost}
                   : std::array<int, 2>{n - ghost, n};
  };
  const auto x = axis(extent.x, offset.x);
  const auto y = axis(extent.y, offset.y);
  const auto z = axis(extent.z, offset.z);
  return Box3{Int3{x[0], y[0], z[0]}, Int3{x[1], y[1], z[1]}};
}

void test_plan_geometry(const MpiContext& context) {
  auto periodic = StructuredDecomposition::create(context, kGlobal, kPeriodic);
  const Int3 extent = periodic.local_extent();
  for (int ghost : {0, 1, 2}) {
    const auto plan = ExchangePlan::create(periodic, extent, ghost);
    HUNDUN_CHECK(plan.ghost_width() == ghost);
    HUNDUN_CHECK(plan.regions().size() == 26U);
    int previous_code = -1;
    std::set<int> codes;
    std::set<std::tuple<int, int, int>> received_cells;
    for (const auto& region : plan.regions()) {
      const int code = hundun::runtime::detail::halo_offset_code(region.offset);
      HUNDUN_CHECK(code > previous_code);
      previous_code = code;
      HUNDUN_CHECK(codes.insert(code).second);
      HUNDUN_CHECK(region.neighbor_rank ==
                   periodic.neighbor_rank(region.offset));
      HUNDUN_CHECK(same(region.send_box,
                        expected_axis_boxes(extent, ghost, region.offset,
                                            false)));
      HUNDUN_CHECK(same(region.receive_box,
                        expected_axis_boxes(extent, ghost, region.offset,
                                            true)));
      HUNDUN_CHECK(box_cells(region.send_box) ==
                   box_cells(region.receive_box));
      for (int k = region.receive_box.begin.z;
           k < region.receive_box.end.z; ++k) {
        for (int j = region.receive_box.begin.y;
             j < region.receive_box.end.y; ++j) {
          for (int i = region.receive_box.begin.x;
               i < region.receive_box.end.x; ++i) {
            HUNDUN_CHECK(!owned_coordinate(Int3{i, j, k}, extent));
            HUNDUN_CHECK(received_cells.insert({i, j, k}).second);
          }
        }
      }
    }
    const auto padded_x = static_cast<std::size_t>(extent.x + 2 * ghost);
    const auto padded_y = static_cast<std::size_t>(extent.y + 2 * ghost);
    const auto padded_z = static_cast<std::size_t>(extent.z + 2 * ghost);
    const auto interior = static_cast<std::size_t>(extent.x) *
                          static_cast<std::size_t>(extent.y) *
                          static_cast<std::size_t>(extent.z);
    HUNDUN_CHECK(received_cells.size() ==
                 padded_x * padded_y * padded_z - interior);
  }

  expect_local_error([&] {
    static_cast<void>(ExchangePlan::create(
        periodic, Int3{extent.x + 1, extent.y, extent.z}, 1));
  });
  expect_local_error([&] {
    static_cast<void>(ExchangePlan::create(periodic, extent, -1));
  });
  expect_local_error([&] {
    static_cast<void>(ExchangePlan::create(
        periodic, extent, std::min({extent.x, extent.y, extent.z}) + 1));
  });

  auto nonperiodic = StructuredDecomposition::create(
      context, kGlobal, std::array<bool, 3>{false, false, false});
  const auto boundary_plan = ExchangePlan::create(
      nonperiodic, nonperiodic.local_extent(), 1);
  bool saw_null = false;
  for (const auto& region : boundary_plan.regions()) {
    if (region.neighbor_rank == MPI_PROC_NULL) {
      saw_null = true;
    }
  }
  int local_saw_null = saw_null ? 1 : 0;
  int any_saw_null = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_saw_null, &any_saw_null, 1, MPI_INT,
                             MPI_MAX, context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(any_saw_null == 1);

  // Independently exchange inverse payload sizes on a q/r-irregular local
  // partition. This is test evidence, not a production peer handshake.
  auto irregular = StructuredDecomposition::create(
      context, Int3{13, 11, 9}, kPeriodic);
  const auto irregular_plan = ExchangePlan::create(
      irregular, irregular.local_extent(), 1);
  for (const auto& region : irregular_plan.regions()) {
    const Int3 inverse_offset{-region.offset.x, -region.offset.y,
                              -region.offset.z};
    const auto inverse = std::find_if(
        irregular_plan.regions().begin(), irregular_plan.regions().end(),
        [inverse_offset](const hundun::runtime::ExchangeRegion& candidate) {
          return same(candidate.offset, inverse_offset);
        });
    HUNDUN_CHECK(inverse != irregular_plan.regions().end());
    const std::array<std::uint64_t, 2> local_inverse{
        box_cells(inverse->receive_box),
        static_cast<std::uint64_t>(inverse->neighbor_rank)};
    std::array<std::uint64_t, 2> peer_inverse{};
    HUNDUN_CHECK(MPI_Sendrecv(
                     local_inverse.data(),
                     static_cast<int>(local_inverse.size()), MPI_UINT64_T,
                     inverse->neighbor_rank,
                     hundun::runtime::detail::halo_offset_code(region.offset),
                     peer_inverse.data(),
                     static_cast<int>(peer_inverse.size()), MPI_UINT64_T,
                     region.neighbor_rank,
                     hundun::runtime::detail::halo_offset_code(region.offset),
                     irregular.comm(), MPI_STATUS_IGNORE) == MPI_SUCCESS);
    if (region.neighbor_rank != MPI_PROC_NULL) {
      int local_rank = -1;
      HUNDUN_CHECK(MPI_Comm_rank(irregular.comm(), &local_rank) == MPI_SUCCESS);
      HUNDUN_CHECK(box_cells(region.send_box) == peer_inverse[0]);
      HUNDUN_CHECK(peer_inverse[1] ==
                   static_cast<std::uint64_t>(local_rank));
    }
  }

  if (context.size() == 1) {
    auto exact_limit = StructuredDecomposition::create(
        context, Int3{INT_MAX - 1, 1, 1},
        std::array<bool, 3>{false, false, false});
    const auto exact_limit_plan = ExchangePlan::create(
        exact_limit, exact_limit.local_extent(), 1);
    const auto positive_x = std::find_if(
        exact_limit_plan.regions().begin(), exact_limit_plan.regions().end(),
        [](const hundun::runtime::ExchangeRegion& region) {
          return same(region.offset, Int3{1, 0, 0});
        });
    HUNDUN_CHECK(positive_x != exact_limit_plan.regions().end());
    HUNDUN_CHECK(positive_x->receive_box.end.x == INT_MAX);

    auto huge = StructuredDecomposition::create(
        context, Int3{INT_MAX, 1, 1},
        std::array<bool, 3>{false, false, false});
    expect_local_error([&] {
      static_cast<void>(ExchangePlan::create(
          huge, huge.local_extent(), 1));
    });
  }
}

void test_irregular_partition_exchange(const MpiContext& context) {
  auto decomposition = StructuredDecomposition::create(
      context, Int3{13, 11, 9}, kPeriodic);
  FieldRegistry registry;
  const FieldId id = registry.declare_field(
      descriptor("irregular", ScalarType::float64, 3U, 2));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_field<double>(storage, id, decomposition, 3);
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  halo.exchange(storage, id);
  check_field<double>(storage, id, decomposition, kPeriodic, 2, 3, 3);
}

template <class T>
void test_typed_periodic_exchange(const MpiContext& context,
                                  ScalarType scalar_type,
                                  std::uint32_t components, int width) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  FieldRegistry registry;
  const FieldId id = registry.declare_field(
      descriptor("typed", scalar_type, components, width));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_field<T>(storage, id, decomposition, 1);
  auto plan = ExchangePlan::create(decomposition,
                                   decomposition.local_extent(), width);
  auto halo = HaloExchange::create(decomposition, std::move(plan));
  halo.exchange(storage, id);
  check_field<T>(storage, id, decomposition, kPeriodic, width, 1, 1);
}

void test_all_scalar_types(const MpiContext& context) {
  for (int width : {1, 2}) {
    test_typed_periodic_exchange<double>(context, ScalarType::float64, 1U,
                                         width);
    test_typed_periodic_exchange<double>(context, ScalarType::float64, 3U,
                                         width);
    test_typed_periodic_exchange<std::int32_t>(
        context, ScalarType::int32, 1U, width);
    test_typed_periodic_exchange<std::int32_t>(
        context, ScalarType::int32, 3U, width);
    test_typed_periodic_exchange<std::uint8_t>(
        context, ScalarType::uint8, 1U, width);
    test_typed_periodic_exchange<std::uint8_t>(
        context, ScalarType::uint8, 3U, width);
  }
}

void test_async_snapshot_and_outer_ghost(const MpiContext& context) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  FieldRegistry registry;
  const FieldId async_id = registry.declare_field(
      descriptor("async", ScalarType::float64, 3U, 2));
  const FieldId sync_id = registry.declare_field(
      descriptor("sync", ScalarType::float64, 3U, 2));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_field<double>(storage, async_id, decomposition, 1);
  initialize_field<double>(storage, sync_id, decomposition, 1);

  auto async_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  auto sync_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  async_halo.begin(storage, async_id);
  context.barrier();
  check_field<double>(storage, async_id, decomposition, kPeriodic, 0, 1, 1);
  set_owned_pattern<double>(storage, async_id, decomposition, 2);
  async_halo.wait(storage, async_id);
  sync_halo.exchange(storage, sync_id);

  check_field<double>(storage, async_id, decomposition, kPeriodic, 1, 2, 1);
  check_field<double>(storage, sync_id, decomposition, kPeriodic, 1, 1, 1);
  const auto async_view = storage.view<double>(async_id);
  const auto sync_view = storage.view<double>(sync_id);
  const Int3 extent = decomposition.local_extent();
  for (int k = -2; k < extent.z + 2; ++k) {
    for (int j = -2; j < extent.y + 2; ++j) {
      for (int i = -2; i < extent.x + 2; ++i) {
        if (!owned_coordinate(Int3{i, j, k}, extent)) {
          for (int component = 0; component < 3; ++component) {
            HUNDUN_CHECK(async_view(i, j, k, component) ==
                         sync_view(i, j, k, component));
          }
        }
      }
    }
  }
}

void test_nonperiodic_and_zero_width(const MpiContext& context) {
  const std::array<bool, 3> nonperiodic{false, false, false};
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, nonperiodic);
  FieldRegistry registry;
  const FieldId id = registry.declare_field(
      descriptor("nonperiodic", ScalarType::float64, 1U, 1));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_field<double>(storage, id, decomposition, 1);

  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  halo.exchange(storage, id);
  check_field<double>(storage, id, decomposition, nonperiodic, 1, 1, 1);

  FieldRegistry zero_registry;
  const FieldId zero_id = zero_registry.declare_field(
      descriptor("zero_width", ScalarType::float64, 1U, 0));
  zero_registry.freeze();
  FieldStorage zero_storage(zero_registry, decomposition.local_extent());
  initialize_field<double>(zero_storage, zero_id, decomposition, 1);
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions zero_options;
  zero_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(zero_options);
  auto zero = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 0));
  zero.begin(zero_storage, zero_id);
  static_cast<void>(expect_collective_error(
      context, [&] { zero.begin(zero_storage, zero_id); }));
  zero.wait(zero_storage, zero_id);
  const auto snapshot = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(snapshot.receive_posts == 0U);
  HUNDUN_CHECK(snapshot.send_posts == 0U);
  check_field<double>(zero_storage, zero_id, decomposition, nonperiodic, 0, 1,
                      1);
  hundun::runtime::detail::set_halo_test_options({});
}

void test_state_machine_and_moves(const MpiContext& context) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  FieldRegistry registry;
  const FieldId first = registry.declare_field(
      descriptor("first", ScalarType::float64, 1U, 1));
  const FieldId second = registry.declare_field(
      descriptor("second", ScalarType::float64, 1U, 1));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  FieldStorage other(registry, decomposition.local_extent());
  FieldRegistry wrong_layout_registry;
  const FieldId wrong_layout_id = wrong_layout_registry.declare_field(
      descriptor("wrong_layout", ScalarType::float64, 3U, 1));
  wrong_layout_registry.freeze();
  FieldStorage wrong_layout(wrong_layout_registry,
                            decomposition.local_extent());
  initialize_field<double>(storage, first, decomposition, 1);
  initialize_field<double>(storage, second, decomposition, 1);
  initialize_field<double>(other, first, decomposition, 1);

  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  static_cast<void>(expect_collective_error(
      context, [&] { halo.wait(storage, first); }));

  halo.begin(storage, first);
  check_field<double>(storage, first, decomposition, kPeriodic, 0, 1, 1);
  static_cast<void>(expect_collective_error(
      context, [&] { halo.begin(storage, first); }));
  static_cast<void>(expect_collective_error(
      context, [&] { halo.exchange(storage, first); }));
  static_cast<void>(expect_collective_error(
      context, [&] { halo.wait(other, first); }));
  static_cast<void>(expect_collective_error(
      context, [&] { halo.wait(wrong_layout, wrong_layout_id); }));
  static_cast<void>(expect_collective_error(
      context, [&] { halo.wait(storage, second); }));
  halo.wait(storage, first);
  check_field<double>(storage, first, decomposition, kPeriodic, 1, 1, 1);

  // The same object is reusable after both state errors and success.
  initialize_field<double>(storage, first, decomposition, 2);
  halo.exchange(storage, first);
  check_field<double>(storage, first, decomposition, kPeriodic, 1, 2, 2);

  // FieldStorage's unique entry array is the pending identity. Moving the
  // owner preserves that identity; the moved-from wrapper must be rejected.
  initialize_field<double>(storage, first, decomposition, 3);
  halo.begin(storage, first);
  FieldStorage moved_storage(std::move(storage));
  static_cast<void>(expect_collective_error(
      context, [&] { halo.wait(storage, first); }));
  halo.wait(moved_storage, first);
  check_field<double>(moved_storage, first, decomposition, kPeriodic, 1, 3, 3);

  initialize_field<double>(moved_storage, first, decomposition, 4);
  halo.begin(moved_storage, first);
  HaloExchange moved_halo(std::move(halo));
  expect_local_error([&] { halo.begin(moved_storage, first); });
  moved_halo.wait(moved_storage, first);
  check_field<double>(moved_storage, first, decomposition, kPeriodic, 1, 4, 4);
}

void test_context_isolation(const MpiContext& context) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  FieldRegistry registry;
  const FieldId id = registry.declare_field(
      descriptor("isolated", ScalarType::float64, 3U, 1));
  registry.freeze();
  FieldStorage first(registry, decomposition.local_extent());
  FieldStorage second(registry, decomposition.local_extent());
  initialize_field<double>(first, id, decomposition, 5);
  initialize_field<double>(second, id, decomposition, 6);

  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions isolation_options;
  isolation_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(isolation_options);
  auto first_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  auto snapshot = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(snapshot.communicator_is_distinct_congruent);
  HUNDUN_CHECK(snapshot.communicator_uses_errors_return);
  auto second_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  snapshot = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(snapshot.communicator_is_distinct_congruent);
  HUNDUN_CHECK(snapshot.communicator_uses_errors_return);

  first_halo.begin(first, id);
  second_halo.begin(second, id);
  second_halo.wait(second, id);
  first_halo.wait(first, id);
  check_field<double>(first, id, decomposition, kPeriodic, 1, 5, 5);
  check_field<double>(second, id, decomposition, kPeriodic, 1, 6, 6);
  hundun::runtime::detail::set_halo_test_options({});
}

void test_external_lifetimes_and_active_destructor(const MpiContext& context) {
  FieldRegistry registry;
  const FieldId id = registry.declare_field(
      descriptor("lifetime", ScalarType::float64, 1U, 1));
  registry.freeze();

  std::optional<HaloExchange> detached_halo;
  std::optional<FieldStorage> detached_storage;
  {
    auto decomposition =
        StructuredDecomposition::create(context, kGlobal, kPeriodic);
    detached_storage.emplace(registry, decomposition.local_extent());
    initialize_field<double>(*detached_storage, id, decomposition, 7);
    auto external_plan = ExchangePlan::create(
        decomposition, decomposition.local_extent(), 1);
    detached_halo.emplace(
        HaloExchange::create(decomposition, std::move(external_plan)));
    detached_halo->begin(*detached_storage, id);
  }
  // Both the decomposition and caller-owned plan are gone. The Halo owns all
  // communication state required to complete.
  detached_halo->wait(*detached_storage, id);
  {
    auto verifier = StructuredDecomposition::create(context, kGlobal, kPeriodic);
    check_field<double>(*detached_storage, id, verifier, kPeriodic, 1, 7, 7);
  }
  detached_halo.reset();
  detached_storage.reset();

  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions lifetime_options;
  lifetime_options.observe = true;
  lifetime_options.inject_cleanup_wait_error_rank =
      context.size() > 1 ? 1 : 0;
  hundun::runtime::detail::set_halo_test_options(lifetime_options);
  std::optional<HaloExchange> active;
  active.emplace(HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1)));
  {
    FieldStorage temporary(registry, decomposition.local_extent());
    initialize_field<double>(temporary, id, decomposition, 8);
    active->begin(temporary, id);
  }
  // Storage is already gone. Destruction drains and discards without touching
  // it, then later communication on the decomposition remains usable.
  active.reset();
  context.barrier();
  HUNDUN_CHECK(
      hundun::runtime::detail::halo_test_snapshot().destructor_drains == 1U);
  const auto cleanup_injections = static_cast<unsigned long long>(
      hundun::runtime::detail::halo_test_snapshot()
          .cleanup_wait_errors_injected);
  unsigned long long total_cleanup_injections = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&cleanup_injections, &total_cleanup_injections, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(total_cleanup_injections == 1U);
  check_only_hook_fired(
      context, hundun::runtime::detail::halo_test_snapshot(),
      ExpectedHook::cleanup_wait);

  hundun::runtime::detail::set_halo_test_options({});
  FieldStorage fresh(registry, decomposition.local_extent());
  initialize_field<double>(fresh, id, decomposition, 9);
  auto fresh_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  fresh_halo.exchange(fresh, id);
  check_field<double>(fresh, id, decomposition, kPeriodic, 1, 9, 9);
}

void test_collective_preflight_mismatches(const MpiContext& context) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);

  const int injection_rank = context.size() > 1 ? 1 : 0;
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions plan_collective_options;
  plan_collective_options.inject_plan_width_first_collective_error_rank =
      injection_rank;
  plan_collective_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(plan_collective_options);
  auto matching_plan = ExchangePlan::create(
      decomposition, decomposition.local_extent(), 1);
  const std::string plan_collective_message = expect_collective_error(
      context, [&] {
        static_cast<void>(
            HaloExchange::create(decomposition, std::move(matching_plan)));
      });
  HUNDUN_CHECK(plan_collective_message.find("agreement") !=
               std::string::npos);
  const auto plan_collective_snapshot =
      hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(plan_collective_snapshot.plan_width_second_collective_entries ==
               1U);
  const auto local_plan_injections = static_cast<unsigned long long>(
      plan_collective_snapshot.plan_width_first_collective_errors_injected);
  unsigned long long total_plan_injections = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local_plan_injections, &total_plan_injections, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(total_plan_injections == 1U);
  check_only_hook_fired(context, plan_collective_snapshot,
                        ExpectedHook::plan_collective);
  hundun::runtime::detail::set_halo_test_options({});

  auto differing_plan = ExchangePlan::create(
      decomposition, decomposition.local_extent(),
      context.rank() == 0 ? 1 : 2);
  static_cast<void>(expect_collective_error(context, [&] {
    static_cast<void>(
        HaloExchange::create(decomposition, std::move(differing_plan)));
  }));

  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions mismatch_options;
  mismatch_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(mismatch_options);

  FieldRegistry normal_registry;
  const FieldId normal_id = normal_registry.declare_field(
      descriptor("normal", ScalarType::float64, 1U, 1));
  normal_registry.freeze();
  FieldStorage normal(normal_registry, decomposition.local_extent());
  initialize_field<double>(normal, normal_id, decomposition, 1);

  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions wire_collective_options;
  wire_collective_options.inject_wire_first_collective_error_rank =
      injection_rank;
  wire_collective_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(wire_collective_options);
  const std::string wire_collective_message = expect_collective_error(
      context, [&] { halo.begin(normal, normal_id); });
  HUNDUN_CHECK(wire_collective_message.find("agreement") !=
               std::string::npos);
  const auto wire_collective_snapshot =
      hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(wire_collective_snapshot.wire_second_collective_entries ==
               1U);
  HUNDUN_CHECK(wire_collective_snapshot.receive_posts == 0U);
  HUNDUN_CHECK(wire_collective_snapshot.send_posts == 0U);
  const auto local_wire_injections = static_cast<unsigned long long>(
      wire_collective_snapshot.wire_first_collective_errors_injected);
  unsigned long long total_wire_injections = 0U;
  HUNDUN_CHECK(MPI_Allreduce(&local_wire_injections, &total_wire_injections, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             context.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(total_wire_injections == 1U);
  check_only_hook_fired(context, wire_collective_snapshot,
                        ExpectedHook::wire_collective);
  hundun::runtime::detail::set_halo_test_options({});
  halo.exchange(normal, normal_id);
  check_field<double>(normal, normal_id, decomposition, kPeriodic, 1, 1, 1);
  initialize_field<double>(normal, normal_id, decomposition, 1);

  FieldRegistry ids_registry;
  const FieldId id_zero = ids_registry.declare_field(
      descriptor("id_zero", ScalarType::float64, 1U, 1));
  const FieldId id_one = ids_registry.declare_field(
      descriptor("id_one", ScalarType::float64, 1U, 1));
  ids_registry.freeze();
  FieldStorage ids(ids_registry, decomposition.local_extent());
  initialize_field<double>(ids, id_zero, decomposition, 1);
  initialize_field<double>(ids, id_one, decomposition, 1);
  static_cast<void>(expect_collective_error(context, [&] {
    halo.begin(ids, context.rank() == 0 ? id_zero : id_one);
  }));

  FieldRegistry scalar_registry;
  const ScalarType rank_scalar =
      context.rank() == 0 ? ScalarType::float64 : ScalarType::int32;
  const FieldId scalar_id = scalar_registry.declare_field(
      descriptor("scalar", rank_scalar, 1U, 1));
  scalar_registry.freeze();
  FieldStorage scalar_storage(scalar_registry, decomposition.local_extent());
  static_cast<void>(expect_collective_error(
      context, [&] { halo.begin(scalar_storage, scalar_id); }));

  FieldRegistry component_registry;
  const FieldId component_id = component_registry.declare_field(descriptor(
      "components", ScalarType::float64,
      context.rank() == 0 ? 1U : 3U, 1));
  component_registry.freeze();
  FieldStorage component_storage(component_registry,
                                 decomposition.local_extent());
  static_cast<void>(expect_collective_error(
      context, [&] { halo.begin(component_storage, component_id); }));

  FieldRegistry ghost_registry;
  const FieldId ghost_id = ghost_registry.declare_field(descriptor(
      "ghost", ScalarType::float64, 1U, context.rank() == 0 ? 1 : 2));
  ghost_registry.freeze();
  FieldStorage ghost_storage(ghost_registry, decomposition.local_extent());
  static_cast<void>(expect_collective_error(
      context, [&] { halo.begin(ghost_storage, ghost_id); }));

  const FieldId invalid_id =
      context.rank() == 0 ? std::numeric_limits<FieldId>::max() : normal_id;
  static_cast<void>(expect_collective_error(
      context, [&] { halo.begin(normal, invalid_id); }));

  Int3 rank_extent = decomposition.local_extent();
  if (context.rank() == 0) {
    ++rank_extent.x;
  }
  FieldStorage wrong_extent(normal_registry, rank_extent);
  static_cast<void>(expect_collective_error(
      context, [&] { halo.begin(wrong_extent, normal_id); }));

  const auto snapshot = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(snapshot.receive_posts == 0U);
  HUNDUN_CHECK(snapshot.send_posts == 0U);
  check_field<double>(normal, normal_id, decomposition, kPeriodic, 0, 1, 1);

  hundun::runtime::detail::set_halo_test_options({});
  halo.exchange(normal, normal_id);
  check_field<double>(normal, normal_id, decomposition, kPeriodic, 1, 1, 1);
}

void test_small_chunks_and_reuse(const MpiContext& context) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  FieldRegistry registry;
  const FieldId wide = registry.declare_field(
      descriptor("wide", ScalarType::float64, 3U, 2));
  const FieldId narrow = registry.declare_field(
      descriptor("narrow", ScalarType::float64, 1U, 2));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  initialize_field<double>(storage, wide, decomposition, 1);
  initialize_field<double>(storage, narrow, decomposition, 1);

  auto observed_plan = ExchangePlan::create(
      decomposition, decomposition.local_extent(), 2);
  std::size_t expected_row_copies = 0U;
  for (const auto& region : observed_plan.regions()) {
    if (region.neighbor_rank != MPI_PROC_NULL &&
        region.send_box.begin.x < region.send_box.end.x) {
      expected_row_copies +=
          static_cast<std::size_t>(region.send_box.end.y -
                                   region.send_box.begin.y) *
          static_cast<std::size_t>(region.send_box.end.z -
                                   region.send_box.begin.z);
    }
  }
  auto halo = HaloExchange::create(decomposition, std::move(observed_plan));
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions small_chunk_options;
  small_chunk_options.chunk_limit = 7U;
  small_chunk_options.waitall_limit = 5U;
  small_chunk_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(small_chunk_options);
  halo.exchange(storage, wide);
  check_field<double>(storage, wide, decomposition, kPeriodic, 2, 1, 1);
  const auto first = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(first.receive_posts > 26U);
  HUNDUN_CHECK(first.receive_posts == first.send_posts);
  HUNDUN_CHECK(first.all_receives_preceded_sends);
  HUNDUN_CHECK(first.first_send_sequence == first.receive_posts);
  HUNDUN_CHECK(first.chunk_offsets_ordered);
  HUNDUN_CHECK(first.pack_row_copy_events == expected_row_copies);
  HUNDUN_CHECK(first.unpack_row_copy_events == expected_row_copies);

  initialize_field<double>(storage, wide, decomposition, 2);
  halo.exchange(storage, wide);
  check_field<double>(storage, wide, decomposition, kPeriodic, 2, 2, 2);
  const auto second = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(second.send_buffer_capacity == first.send_buffer_capacity);
  HUNDUN_CHECK(second.receive_buffer_capacity == first.receive_buffer_capacity);
  HUNDUN_CHECK(second.request_capacity == first.request_capacity);
  HUNDUN_CHECK(second.chunk_metadata_capacity ==
               first.chunk_metadata_capacity);
  HUNDUN_CHECK(second.wait_batch_capacity == first.wait_batch_capacity);

  // A different field may resize the reusable working set only while idle.
  halo.exchange(storage, narrow);
  check_field<double>(storage, narrow, decomposition, kPeriodic, 2, 1, 1);
  hundun::runtime::detail::set_halo_test_options({});
}

void test_recoverable_failures(const MpiContext& context) {
  auto decomposition =
      StructuredDecomposition::create(context, kGlobal, kPeriodic);
  FieldRegistry registry;
  const FieldId id = registry.declare_field(
      descriptor("failure", ScalarType::float64, 3U, 1));
  registry.freeze();
  FieldStorage storage(registry, decomposition.local_extent());
  const int injection_rank = context.size() > 1 ? 1 : 0;

  initialize_field<double>(storage, id, decomposition, 1);
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions post_options;
  post_options.inject_post_error_rank = injection_rank;
  post_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(post_options);
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  HUNDUN_CHECK(hundun::runtime::detail::halo_test_snapshot()
                   .context_generation == 1U);
  const std::string post_message = expect_collective_error(
      context, [&] { halo.begin(storage, id); });
  HUNDUN_CHECK(post_message.find("post") != std::string::npos);
  HUNDUN_CHECK(post_message.find("rank=" + std::to_string(injection_rank)) !=
               std::string::npos);
  HUNDUN_CHECK(post_message.find("operation=MPI_Irecv") !=
               std::string::npos);
  HUNDUN_CHECK(post_message.find("result=" + std::to_string(MPI_ERR_OTHER)) !=
               std::string::npos);
  HUNDUN_CHECK(post_message.find("region=0") != std::string::npos);
  HUNDUN_CHECK(post_message.find("chunk_offset=0") != std::string::npos);
  HUNDUN_CHECK(post_message.find("tag=26") != std::string::npos);
  check_field<double>(storage, id, decomposition, kPeriodic, 0, 1, 1);
  const auto post_snapshot = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(post_snapshot.context_generation == 2U);
  HUNDUN_CHECK(post_snapshot.context_replacements == 1U);
  HUNDUN_CHECK(post_snapshot.last_context_replacement_distinct_congruent);
  HUNDUN_CHECK(post_snapshot.cancel_calls > 0U);
  HUNDUN_CHECK(post_snapshot.cancellation_status_checks ==
               post_snapshot.cancel_calls);
  HUNDUN_CHECK(post_snapshot.cancelled_requests +
                   post_snapshot.completed_requests ==
               post_snapshot.cancellation_status_checks);
  check_only_hook_fired(context, post_snapshot, ExpectedHook::post);

  hundun::runtime::detail::set_halo_test_options({});
  halo.exchange(storage, id);
  check_field<double>(storage, id, decomposition, kPeriodic, 1, 1, 1);

  initialize_field<double>(storage, id, decomposition, 2);
  hundun::runtime::detail::reset_halo_test_observation();
  hundun::runtime::detail::HaloTestOptions wait_options;
  wait_options.inject_wait_error_rank = injection_rank;
  wait_options.observe = true;
  hundun::runtime::detail::set_halo_test_options(wait_options);
  halo.begin(storage, id);
  HUNDUN_CHECK(hundun::runtime::detail::halo_test_snapshot()
                   .context_generation == 2U);
  const std::string wait_message = expect_collective_error(
      context, [&] { halo.wait(storage, id); });
  HUNDUN_CHECK(wait_message.find("completion") != std::string::npos);
  HUNDUN_CHECK(wait_message.find("rank=" + std::to_string(injection_rank)) !=
               std::string::npos);
  HUNDUN_CHECK(wait_message.find("operation=MPI_Waitall") !=
               std::string::npos);
  HUNDUN_CHECK(wait_message.find("result=" + std::to_string(MPI_ERR_OTHER)) !=
               std::string::npos);
  HUNDUN_CHECK(wait_message.find("chunk_offset=0") != std::string::npos);
  check_field<double>(storage, id, decomposition, kPeriodic, 0, 2, 2);
  const auto wait_snapshot = hundun::runtime::detail::halo_test_snapshot();
  HUNDUN_CHECK(wait_snapshot.context_generation == 3U);
  HUNDUN_CHECK(wait_snapshot.context_replacements == 1U);
  HUNDUN_CHECK(wait_snapshot.last_context_replacement_distinct_congruent);
  check_only_hook_fired(context, wait_snapshot, ExpectedHook::wait);

  hundun::runtime::detail::set_halo_test_options({});
  halo.exchange(storage, id);
  check_field<double>(storage, id, decomposition, kPeriodic, 1, 2, 2);
}

void run_full(const MpiContext& context) {
  HUNDUN_CHECK(context.size() == 1 || context.size() == 2 ||
               context.size() == 4);
  test_plan_geometry(context);
  test_all_scalar_types(context);
  test_irregular_partition_exchange(context);
  test_async_snapshot_and_outer_ghost(context);
  test_nonperiodic_and_zero_width(context);
}

void run_state(const MpiContext& context) {
  test_state_machine_and_moves(context);
  test_context_isolation(context);
}

void run_lifetime(const MpiContext& context) {
  test_external_lifetimes_and_active_destructor(context);
}

void run_mismatch(const MpiContext& context) {
  test_collective_preflight_mismatches(context);
}

void run_small_chunks(const MpiContext& context) {
  test_small_chunks_and_reuse(context);
}

void run_failure(const MpiContext& context) {
  test_recoverable_failures(context);
}

int run_finalized_idle(int argc, char** argv) {
  std::optional<HaloExchange> halo;
  int active_result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    auto context = MpiContext::duplicate(MPI_COMM_WORLD);
    active_result = hundun::test::run([&] {
      auto decomposition =
          StructuredDecomposition::create(context, kGlobal, kPeriodic);
      halo.emplace(HaloExchange::create(
          decomposition,
          ExchangePlan::create(decomposition, decomposition.local_extent(),
                               0)));
    });
  }
  if (active_result != EXIT_SUCCESS) {
    return active_result;
  }
  return hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);
    halo.reset();
  });
}

}  // namespace

int main(int argc, char** argv) {
  const std::string_view requested_mode = argc > 1 ? argv[1] : "full";
  if (requested_mode == "finalized_idle") {
    return run_finalized_idle(argc, argv);
  }
  MpiEnvironment environment(argc, argv);
  auto context = MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    if (requested_mode == "full") {
      run_full(context);
      return;
    }
    if (requested_mode == "state") {
      run_state(context);
      return;
    }
    if (requested_mode == "lifetime") {
      run_lifetime(context);
      return;
    }
    if (requested_mode == "mismatch") {
      run_mismatch(context);
      return;
    }
    if (requested_mode == "small_chunks") {
      run_small_chunks(context);
      return;
    }
    if (requested_mode == "failure") {
      run_failure(context);
      return;
    }
    throw Error("unknown halo MPI test mode");
  });
}
