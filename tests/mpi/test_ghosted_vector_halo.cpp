// SPDX-License-Identifier: Apache-2.0

#include "hundun/execution/execution.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/linear/ghosted_vector_halo.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "linear/src/ghosted_vector_halo_detail.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::execution::BackendIdentity;
using hundun::execution::CpuReferenceContext;
using hundun::execution::ExecutionCapability;
using hundun::execution::ExecutionContext;
using hundun::execution::ExecutionSpace;
using hundun::linear::BufferHaloPath;
using hundun::linear::GhostedVector;
using hundun::linear::GhostedVectorHalo;
using hundun::linear::VectorLayout;
using hundun::mesh::MeshTopology;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Error;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kExtent{7, 5, 3};

static_assert(!std::is_copy_constructible_v<GhostedVectorHalo>);
static_assert(!std::is_copy_assignable_v<GhostedVectorHalo>);
static_assert(std::is_nothrow_move_constructible_v<GhostedVectorHalo>);
static_assert(!std::is_move_assignable_v<GhostedVectorHalo>);

Int3 process_grid_for(int ranks) {
  switch (ranks) {
    case 1:
      return {1, 1, 1};
    case 2:
      return {2, 1, 1};
    case 4:
      return {2, 2, 1};
    default:
      throw Error("unsupported GhostedVector Halo test rank count");
  }
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
  HUNDUN_CHECK(message.size() < 512U);
  std::array<char, 512> reference{};
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
    HUNDUN_CHECK(!std::string(error.what()).empty());
  }
  HUNDUN_CHECK(threw);
}

double sample(hundun::mesh::GlobalCellId id, int pattern) {
  return static_cast<double>(id) + 10000.0 * static_cast<double>(pattern);
}

void fill_vector(GhostedVector& vector, const MeshTopology& topology,
                 int pattern, double ghost_sentinel) {
  auto view = vector.local_view();
  for (std::size_t local = 0; local < topology.owned_cell_count(); ++local) {
    view[local] = sample(topology.global_cell_id(local), pattern);
  }
  for (std::size_t local = topology.owned_cell_count();
       local < topology.local_cell_count(); ++local) {
    view[local] = ghost_sentinel;
  }
}

void check_vector(const GhostedVector& vector, const MeshTopology& topology,
                  int pattern, std::optional<std::size_t> changed_owned = {},
                  double changed_delta = 0.0) {
  const auto view = vector.local_view();
  HUNDUN_CHECK(view.size() == topology.local_cell_count());
  for (std::size_t local = 0; local < topology.local_cell_count(); ++local) {
    double expected = sample(topology.global_cell_id(local), pattern);
    if (changed_owned.has_value() && local == *changed_owned) {
      expected += changed_delta;
    }
    HUNDUN_CHECK_NEAR(view[local], expected, 0.0);
  }
}

void check_ghost_sentinel(const GhostedVector& vector,
                          const MeshTopology& topology, double sentinel) {
  const auto ghosts = vector.ghost_view();
  HUNDUN_CHECK(ghosts.size() == topology.ghost_cell_count());
  for (std::size_t index = 0; index < ghosts.size(); ++index) {
    HUNDUN_CHECK_NEAR(ghosts[index], sentinel, 0.0);
  }
}

class AlternateHostContext final : public ExecutionContext {
 public:
  explicit AlternateHostContext(BackendIdentity identity) noexcept
      : identity_(identity) {}
  std::string_view backend_name() const noexcept override {
    return "alternate_host_test";
  }
  BackendIdentity backend_identity() const noexcept override {
    return identity_;
  }
  ExecutionSpace space() const noexcept override {
    return ExecutionSpace::host;
  }
  bool ordered() const noexcept override { return true; }
  bool supports(ExecutionCapability capability) const noexcept override {
    return capability == ExecutionCapability::buffer_allocation ||
           capability == ExecutionCapability::host_access ||
           capability == ExecutionCapability::transfer;
  }

 private:
  BackendIdentity identity_;
};

void run_exchange_case(const MpiContext& world,
                       std::array<bool, 3> periodic) {
  const auto decomposition = StructuredDecomposition::create(
      world, kExtent, periodic,
      DecompositionOptions{process_grid_for(world.size())});
  const MeshTopology topology(decomposition);
  CpuReferenceContext execution;
  const VectorLayout layout = VectorLayout::from_topology(topology);
  HUNDUN_CHECK(layout.owned_count() == topology.owned_cell_count());
  HUNDUN_CHECK(layout.ghost_count() == topology.ghost_cell_count());
  for (std::size_t local = 0; local < layout.local_count(); ++local) {
    HUNDUN_CHECK(layout.global_ids()[local] == topology.global_cell_id(local));
  }

  GhostedVector vector(execution, layout);
  GhostedVectorHalo halo =
      GhostedVectorHalo::create(decomposition, topology, execution);
  HUNDUN_CHECK(halo.path() == BufferHaloPath::host_direct);
  HUNDUN_CHECK(halo.owned_count() == topology.owned_cell_count());
  HUNDUN_CHECK(halo.ghost_count() == topology.ghost_cell_count());
  HUNDUN_CHECK(halo.receive_value_count() == topology.ghost_cell_count());
  const auto identity = vector.allocation_identity();
  const auto epoch = vector.epoch();

  fill_vector(vector, topology, 1, -1.0);
  halo.exchange(vector);
  check_vector(vector, topology, 1);
  HUNDUN_CHECK(vector.allocation_identity() == identity);
  HUNDUN_CHECK(vector.epoch() == epoch);

  hundun::linear::detail::set_vector_halo_test_options(
      hundun::linear::detail::VectorHaloTestOptions{1U, true});
  hundun::linear::detail::reset_vector_halo_test_observation();
  GhostedVectorHalo small_chunk_halo =
      GhostedVectorHalo::create(decomposition, topology, execution);
  fill_vector(vector, topology, 2, -2.0);
  small_chunk_halo.exchange(vector);
  check_vector(vector, topology, 2);
  const auto first_snapshot =
      hundun::linear::detail::vector_halo_test_snapshot();
  HUNDUN_CHECK(first_snapshot.receives_preceded_sends);
  HUNDUN_CHECK(first_snapshot.chunk_offsets_ordered);
  HUNDUN_CHECK(first_snapshot.receive_posts >=
               (topology.ghost_cell_count() == 0U ? 0U : 1U));
  HUNDUN_CHECK(first_snapshot.send_wire_identities.size() ==
               first_snapshot.receive_wire_identities.size());
  std::size_t received_values = 0U;
  std::size_t sent_values = 0U;
  bool send_seen = false;
  std::vector<std::size_t> next_receive_offset(
      static_cast<std::size_t>(world.size()), 0U);
  std::vector<std::size_t> next_send_offset(
      static_cast<std::size_t>(world.size()), 0U);
  for (const auto& event : first_snapshot.post_events) {
    HUNDUN_CHECK(event.peer >= 0 && event.peer < world.size());
    HUNDUN_CHECK(event.count == 1);
    HUNDUN_CHECK(event.tag >= 0);
    const std::size_t peer = static_cast<std::size_t>(event.peer);
    if (event.receive) {
      HUNDUN_CHECK(!send_seen);
      HUNDUN_CHECK(event.offset == next_receive_offset[peer]);
      next_receive_offset[peer] += static_cast<std::size_t>(event.count);
      received_values += static_cast<std::size_t>(event.count);
    } else {
      send_seen = true;
      HUNDUN_CHECK(event.offset == next_send_offset[peer]);
      next_send_offset[peer] += static_cast<std::size_t>(event.count);
      sent_values += static_cast<std::size_t>(event.count);
    }
  }
  HUNDUN_CHECK(received_values == small_chunk_halo.receive_value_count());
  HUNDUN_CHECK(sent_values == small_chunk_halo.send_value_count());
  HUNDUN_CHECK(first_snapshot.post_events.size() ==
               first_snapshot.receive_posts + first_snapshot.send_posts);

  fill_vector(vector, topology, 3, -3.0);
  small_chunk_halo.begin(vector);
  check_ghost_sentinel(vector, topology, -3.0);
  std::optional<std::size_t> changed;
  if (vector.owned_count() != 0U) {
    changed = 0U;
    vector.owned_view()[0] += 777.0;
  }
  small_chunk_halo.wait(vector);
  check_vector(vector, topology, 3, changed, changed ? 777.0 : 0.0);
  const auto second_snapshot =
      hundun::linear::detail::vector_halo_test_snapshot();
  HUNDUN_CHECK(second_snapshot.send_wire_identities ==
               first_snapshot.send_wire_identities);
  HUNDUN_CHECK(second_snapshot.receive_wire_identities ==
               first_snapshot.receive_wire_identities);
  HUNDUN_CHECK(second_snapshot.request_capacity ==
               first_snapshot.request_capacity);

  GhostedVectorHalo other =
      GhostedVectorHalo::create(decomposition, topology, execution);
  fill_vector(vector, topology, 4, -4.0);
  small_chunk_halo.begin(vector);
  other.begin(vector);
  other.wait(vector);
  small_chunk_halo.wait(vector);
  check_vector(vector, topology, 4);

  {
    GhostedVectorHalo draining =
        GhostedVectorHalo::create(decomposition, topology, execution);
    fill_vector(vector, topology, 5, -5.0);
    draining.begin(vector);
  }
  fill_vector(vector, topology, 6, -6.0);
  small_chunk_halo.exchange(vector);
  check_vector(vector, topology, 6);
  hundun::linear::detail::set_vector_halo_test_options({});
}

void run_full(const MpiContext& world) {
  run_exchange_case(world, {true, true, true});
  run_exchange_case(world, {false, false, false});
}

void run_mismatch(const MpiContext& world) {
  const auto decomposition = StructuredDecomposition::create(
      world, kExtent, {true, true, true},
      DecompositionOptions{process_grid_for(world.size())});
  const MeshTopology topology(decomposition);
  CpuReferenceContext execution;
  const VectorLayout layout = VectorLayout::from_topology(topology);
  GhostedVector vector(execution, layout);
  GhostedVectorHalo halo =
      GhostedVectorHalo::create(decomposition, topology, execution);

  expect_collective_error(world, [&] { halo.wait(vector); });

  fill_vector(vector, topology, 1, -1.0);
  const std::string operation_message = expect_collective_error(world, [&] {
    if (world.rank() == 0) {
      halo.begin(vector);
    } else {
      halo.exchange(vector);
    }
  });
  HUNDUN_CHECK(operation_message.find("operation") != std::string::npos);

  std::vector<hundun::mesh::GlobalCellId> wrong_ids = layout.global_ids();
  if (world.rank() == 1 && !wrong_ids.empty()) {
    wrong_ids[0] += topology.global_cell_count() + 100U;
  }
  GhostedVector wrong_layout_vector(
      execution, VectorLayout(layout.owned_count(), std::move(wrong_ids)));
  const std::string layout_message = expect_collective_error(
      world, [&] { halo.exchange(world.rank() == 1 ? wrong_layout_vector
                                                   : vector); });
  HUNDUN_CHECK(layout_message.find("layout") != std::string::npos);

  std::size_t wrong_owned = layout.owned_count();
  if (world.rank() == 1 && layout.local_count() != 0U) {
    wrong_owned = wrong_owned == 0U ? 1U : wrong_owned - 1U;
  }
  GhostedVector wrong_count_vector(
      execution, VectorLayout(wrong_owned, layout.global_ids()));
  const std::string count_message = expect_collective_error(
      world, [&] { halo.exchange(world.rank() == 1 ? wrong_count_vector
                                                   : vector); });
  HUNDUN_CHECK(count_message.find("layout") != std::string::npos);

  AlternateHostContext alternate(execution.backend_identity() + 99U);
  GhostedVector alternate_vector(alternate, layout);
  const std::string backend_message = expect_collective_error(
      world, [&] { halo.exchange(world.rank() == 1 ? alternate_vector
                                                   : vector); });
  HUNDUN_CHECK(backend_message.find("backend") != std::string::npos);

  fill_vector(vector, topology, 2, -2.0);
  halo.begin(vector);
  expect_collective_error(world, [&] { halo.begin(vector); });

  GhostedVector target(execution, layout);
  const std::string target_message = expect_collective_error(
      world,
      [&] { halo.wait(world.rank() == 1 ? target : vector); });
  HUNDUN_CHECK(target_message.find("target") != std::string::npos);
  halo.wait(vector);
  check_vector(vector, topology, 2);

  GhostedVector moved_target(std::move(vector));
  fill_vector(moved_target, topology, 3, -3.0);
  halo.begin(moved_target);
  GhostedVector moved_again(std::move(moved_target));
  halo.wait(moved_again);
  check_vector(moved_again, topology, 3);
}

void run_failure(const MpiContext& world) {
  const auto decomposition = StructuredDecomposition::create(
      world, kExtent, {true, true, true},
      DecompositionOptions{process_grid_for(world.size())});
  const MeshTopology topology(decomposition);
  CpuReferenceContext execution;
  const VectorLayout layout = VectorLayout::from_topology(topology);

  hundun::linear::detail::VectorHaloTestOptions options{};
  options.inject_request_id_mismatch_rank = 1;
  hundun::linear::detail::set_vector_halo_test_options(options);
  const std::string request_message = expect_collective_error(world, [&] {
    static_cast<void>(
        GhostedVectorHalo::create(decomposition, topology, execution));
  });
  HUNDUN_CHECK(request_message.find("request") != std::string::npos);

  options = {};
  options.observe = true;
  hundun::linear::detail::set_vector_halo_test_options(options);
  GhostedVectorHalo halo =
      GhostedVectorHalo::create(decomposition, topology, execution);
  GhostedVector vector(execution, layout);

  options.inject_post_failure_rank = 1;
  options.observe = true;
  hundun::linear::detail::set_vector_halo_test_options(options);
  fill_vector(vector, topology, 1, -1.0);
  const std::string post_message =
      expect_collective_error(world, [&] { halo.begin(vector); });
  HUNDUN_CHECK(post_message.find("post") != std::string::npos);
  HUNDUN_CHECK(post_message.find("rank=1") != std::string::npos);
  check_ghost_sentinel(vector, topology, -1.0);
  HUNDUN_CHECK(hundun::linear::detail::vector_halo_test_snapshot()
                   .context_replacements == 1U);

  options = {};
  options.inject_completion_failure_rank = 2;
  options.observe = true;
  hundun::linear::detail::set_vector_halo_test_options(options);
  fill_vector(vector, topology, 2, -2.0);
  halo.begin(vector);
  const std::string completion_message =
      expect_collective_error(world, [&] { halo.wait(vector); });
  HUNDUN_CHECK(completion_message.find("completion") != std::string::npos);
  HUNDUN_CHECK(completion_message.find("rank=2") != std::string::npos);
  check_ghost_sentinel(vector, topology, -2.0);
  HUNDUN_CHECK(hundun::linear::detail::vector_halo_test_snapshot()
                   .context_replacements == 2U);

  hundun::linear::detail::set_vector_halo_test_options({});
  fill_vector(vector, topology, 3, -3.0);
  halo.exchange(vector);
  check_vector(vector, topology, 3);
}

int run_finalized_idle(int argc, char** argv) {
  std::optional<GhostedVectorHalo> halo;
  std::optional<StructuredDecomposition> decomposition;
  std::optional<MeshTopology> topology;
  std::optional<MpiContext> context;
  std::optional<GhostedVector> vector;
  CpuReferenceContext execution;
  int active_result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    context.emplace(MpiContext::duplicate(MPI_COMM_WORLD));
    active_result = hundun::test::run([&] {
      decomposition.emplace(StructuredDecomposition::create(
          *context, kExtent, {false, false, false},
          DecompositionOptions{Int3{1, 1, 1}}));
      topology.emplace(*decomposition);
      vector.emplace(execution, VectorLayout::from_topology(*topology));
      halo.emplace(
          GhostedVectorHalo::create(*decomposition, *topology, execution));
      GhostedVectorHalo moved(std::move(*halo));
      expect_local_error([&] { static_cast<void>(halo->path()); });
      halo.emplace(std::move(moved));
      HUNDUN_CHECK(halo->path() == BufferHaloPath::host_direct);
    });
  }
  if (active_result != EXIT_SUCCESS) {
    return active_result;
  }
  return hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);
    HUNDUN_CHECK(halo->owned_count() == topology->owned_cell_count());
    halo.reset();
    vector.reset();
    topology.reset();
    decomposition.reset();
    context.reset();
  });
}

}  // namespace

int main(int argc, char** argv) {
  const std::string_view mode = argc > 1 ? argv[1] : "full";
  if (mode == "finalized_idle") {
    return run_finalized_idle(argc, argv);
  }
  MpiEnvironment environment(argc, argv);
  auto world = MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    if (mode == "full") {
      run_full(world);
    } else if (mode == "mismatch") {
      run_mismatch(world);
    } else if (mode == "failure") {
      run_failure(world);
    } else {
      throw Error("unknown GhostedVector Halo test mode");
    }
  });
}
