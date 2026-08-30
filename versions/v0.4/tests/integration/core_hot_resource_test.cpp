// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_app.hpp"

#include "../support/product_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <new>

namespace allocation_observer {

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
  if (posix_memalign(&pointer, alignment, size == 0U ? alignment : size) == 0 &&
      pointer != nullptr)
    return pointer;
  throw std::bad_alloc{};
}

}  // namespace allocation_observer

void* operator new(std::size_t size) { return allocation_observer::allocate(size); }
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
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

bool expect(bool condition, const char* description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

ResourceCounters observe_contract(const FrozenExecutionGraph& graph) noexcept {
  ResourceCounters counters;
  counters.peak_workspace_bytes = graph.resources().max_live_workspace_bytes;
  const Span<const FrozenStage> stages = graph.stages();
  for (std::size_t index = 0U; index < stages.size; ++index) {
    const FrozenStage& stage = stages.data[index];
    ResourceCounters increment;
    increment.merged_halo_messages = stage.resources.merged_halo_messages;
    increment.merged_halo_bytes = stage.resources.merged_halo_bytes;
    increment.numeric_refills = stage.resources.numeric_refills;
    increment.hierarchy_rebuilds = stage.resources.hierarchy_rebuilds;
    increment.cache_publishes = stage.resources.cache_publishes;
    increment.linear_iterations = stage.resources.linear_iterations;
    increment.stage_wall_nanoseconds = stage.resources.stage_wall_nanoseconds;
    if (!add_resource_counters(counters, increment)) {
      counters.allocations = UINT64_MAX;
      return counters;
    }
  }
  return counters;
}

bool run() {
  CompiledCasePlan plan;
  const Status status = ProductCompiler::compile(
      MPI_COMM_SELF, test::product_model({18, 12, 8}), {}, plan);
  if (!expect(static_cast<bool>(status), "resource fixture freezes"))
    return false;
  const FrozenExecutionGraph* graph = plan.execution_graph();
  if (!expect(graph != nullptr, "sealed graph is available")) return false;
  const std::uintptr_t state = plan.state_storage_address();
  const std::uintptr_t krylov = plan.krylov_storage_address();
  const std::uintptr_t mg = plan.mg_storage_address();
  const FrozenStage* first_corrector = graph->stage(40U);
  const FrozenStage* final_corrector = graph->stage(50U);
  bool passed = true;
  passed &= expect(first_corrector != nullptr && final_corrector != nullptr &&
                       first_corrector->resources.cache_publishes == 0U &&
                       final_corrector->resources.cache_publishes == 1U,
                   "only corrector two publishes current-attempt final flux");
  passed &= expect(!plan.summary().exact_numeric_certified &&
                       !plan.summary().preconditioner_setup_certified,
                   "never-filled numeric state cannot masquerade as current");

  allocation_observer::count.store(0U, std::memory_order_relaxed);
  allocation_observer::enabled.store(true, std::memory_order_release);
  for (std::size_t repetition = 0U; repetition < 100U; ++repetition) {
    const ResourceCounters counters = observe_contract(*graph);
    passed &= static_cast<bool>(
        validate_resource_counters(graph->resources(), counters));
    passed &= plan.state_storage_address() == state;
    passed &= plan.krylov_storage_address() == krylov;
    passed &= plan.mg_storage_address() == mg;
    passed &= plan.summary().pressure_correctors == 2U;
  }
  allocation_observer::enabled.store(false, std::memory_order_release);
  passed &= expect(allocation_observer::count.load(std::memory_order_relaxed) ==
                       0U,
                   "100 hot contract traversals allocate nothing");
  passed &= expect(plan.state_storage_address() == state &&
                       plan.krylov_storage_address() == krylov &&
                       plan.mg_storage_address() == mg,
                   "all sealed storage addresses remain stable");

  ResourceCounters late = observe_contract(*graph);
  late.allocations = graph->resources().allocation_allowance + 1U;
  passed &= expect(validate_resource_counters(graph->resources(), late).code ==
                       StatusCode::invalid_plan,
                   "late capacity use is rejected");

  ProductDriver driver;
  Status runtime_status =
      ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (runtime_status) runtime_status = driver.initialize(initial);
  DriverStepReport warmup;
  if (runtime_status)
    runtime_status =
        driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, warmup);
  passed &= expect(runtime_status && warmup.accepted &&
                       warmup.piso.pressure_solve_calls == 2U,
                   "cold numeric setup completes before measured hot steps");
  CommittedOutputSnapshot baseline_output;
  RestartSnapshot baseline_restart;
  if (runtime_status)
    runtime_status = driver.committed_output_snapshot(baseline_output);
  if (runtime_status)
    runtime_status = driver.committed_restart_snapshot(baseline_restart);
  std::array<StorageIdentity, 3U> state_storage{};
  for (std::size_t index = 0U;
       index < baseline_output.fields.size && index < state_storage.size();
       ++index)
    state_storage[index] =
        baseline_output.fields.data[index].values.storage_identity;
  const StorageIdentity flux_storage =
      baseline_restart.final_mass_flux.x.storage_identity;

  allocation_observer::count.store(0U, std::memory_order_relaxed);
  allocation_observer::enabled.store(true, std::memory_order_release);
  DriverStepReport last_step;
  for (std::size_t repetition = 0U; repetition < 100U && runtime_status;
       ++repetition) {
    DriverStepReport step;
    runtime_status =
        driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
    last_step = step;
    if (!runtime_status || !step.accepted || step.attempts != 1U ||
        step.piso.pressure_solve_calls != 2U)
      break;
    CommittedOutputSnapshot output;
    RestartSnapshot restart;
    runtime_status = driver.committed_output_snapshot(output);
    if (runtime_status)
      runtime_status = driver.committed_restart_snapshot(restart);
    for (std::size_t index = 0U;
         index < output.fields.size && index < state_storage.size(); ++index)
      if (output.fields.data[index].values.storage_identity !=
          state_storage[index])
        runtime_status = {StatusCode::invalid_plan, 1U};
    if (restart.final_mass_flux.x.storage_identity != flux_storage ||
        restart.final_mass_flux.y.storage_identity != flux_storage ||
        restart.final_mass_flux.z.storage_identity != flux_storage)
      runtime_status = {StatusCode::invalid_plan, 2U};
  }
  allocation_observer::enabled.store(false, std::memory_order_release);
  if (!runtime_status)
    std::cerr << "hot status=" << static_cast<unsigned>(runtime_status.code)
              << "/" << runtime_status.detail
              << " stage=" << last_step.failed_stage
              << " attempts=" << last_step.attempts << '\n';
  passed &= expect(static_cast<bool>(runtime_status),
                   "100 real hot steps complete");
  passed &= expect(allocation_observer::count.load(std::memory_order_relaxed) ==
                       0U,
                   "100 real hot steps allocate no C++ heap storage");
  return passed;
}

bool run_immersed_hot_loop() {
  ValidatedModel model = test::product_model({16, 16, 16});
  model.mesh.lower = {-2.0, -2.0, -2.0};
  model.mesh.upper = {2.0, 2.0, 2.0};
  model.mesh.minimum_spacing = {0.25, 0.25, 0.25};
  model.time.initial_dt = 2.5e-4;
  model.immersed_boundary = ImmersedBoundarySpec{
      "cylinder_ascii.stl", ImmersedFluidSide::outside};
  const std::filesystem::path data_root =
      std::filesystem::path{HUNDUN_V04_SOURCE_ROOT} / "tests" / "data";
  CompiledCasePlan plan;
  Status status =
      ProductCompiler::compile(MPI_COMM_SELF, model, data_root, plan);
  ProductDriver driver;
  if (status)
    status = ProductDriver::create(MPI_COMM_SELF, std::move(plan), driver);
  DriverInitialState initial;
  initial.pressure_reference = 98000.0;
  initial.temperature = 315.0;
  if (status) status = driver.initialize(initial);
  DriverStepReport warmup;
  if (status)
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, warmup);
  bool passed = expect(status && warmup.accepted,
                       "immersed numeric and wall-law warmup completes");
  CommittedOutputSnapshot baseline;
  RestartSnapshot restart;
  if (status) status = driver.committed_output_snapshot(baseline);
  if (status) status = driver.committed_restart_snapshot(restart);
  const StorageIdentity velocity_storage =
      status ? baseline.fields.data[0U].values.storage_identity : 0U;
  const StorageIdentity flux_storage =
      status ? restart.final_mass_flux.x.storage_identity : 0U;

  allocation_observer::count.store(0U, std::memory_order_relaxed);
  allocation_observer::enabled.store(true, std::memory_order_release);
  DriverStepReport last;
  for (std::size_t repetition = 0U; repetition < 100U && status;
       ++repetition) {
    status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, last);
    if (!status || !last.accepted || last.attempts != 1U ||
        last.piso.pressure_solve_calls != 2U ||
        last.resources.ibm_exchanges == 0U ||
        last.resources.structured_exchanges == 0U ||
        last.resources.reduction_collectives == 0U)
      break;
    CommittedOutputSnapshot output;
    RestartSnapshot current;
    status = driver.committed_output_snapshot(output);
    if (status) status = driver.committed_restart_snapshot(current);
    if (status &&
        (output.fields.data[0U].values.storage_identity != velocity_storage ||
         current.final_mass_flux.x.storage_identity != flux_storage))
      status = {StatusCode::invalid_plan, 3U};
  }
  allocation_observer::enabled.store(false, std::memory_order_release);
  if (!status)
    std::cerr << "immersed hot status="
              << static_cast<unsigned>(status.code) << '/' << status.detail
              << " stage=" << last.failed_stage
              << " attempts=" << last.attempts << '\n';
  passed &= expect(static_cast<bool>(status),
                   "100 immersed production hot steps complete");
  passed &= expect(allocation_observer::count.load(std::memory_order_relaxed) ==
                       0U,
                   "100 immersed hot steps allocate no C++ heap storage");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const bool passed = run() && run_immersed_hot_loop();
  MPI_Finalize();
  return passed ? 0 : 1;
}
