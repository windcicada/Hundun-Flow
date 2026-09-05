// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include <mpi.h>

#include <array>
#include <cstdlib>
#include <iostream>
#include <new>

#include "../support/product_fixture.hpp"
#include "hundun/v04_app.hpp"
#include "hundun/v04_linear.hpp"

namespace {
struct Allocation {
  void* base{};
  std::size_t bytes{};
};
std::array<Allocation, 128U> allocations{};
std::size_t allocation_count = 0U;
bool observing = false;
}  // namespace
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  void* pointer = nullptr;
  if (::posix_memalign(&pointer, static_cast<std::size_t>(alignment),
                       bytes == 0U ? 1U : bytes) != 0)
    throw std::bad_alloc{};
  if (observing && allocation_count < allocations.size())
    allocations[allocation_count++] = {pointer, bytes};
  return pointer;
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}

int main(int argc, char** argv) {
  using namespace hundun::v04;
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  bool passed = true;
  {
    ValidatedModel model = test::product_model({16, 16, 16});
    CompiledCasePlan plan;
    observing = true;
    const Status compiled =
        ProductCompiler::compile(MPI_COMM_WORLD, model, {}, plan);
    observing = false;
    passed = static_cast<bool>(compiled);
    const ArenaFieldLayout* mg = nullptr;
    if (passed)
      for (const FieldDescriptor& field : *plan.field_schema())
        if (field.stable_name == "mg_arena")
          mg = plan.arena_layout()->field(field.id);
    passed &= mg != nullptr;
    if (mg != nullptr) {
      const auto* owner = plan.arena_layout()->owner(mg->owner_index);
      const std::uintptr_t owner_base =
          plan.mg_storage_address() - mg->offset_doubles * sizeof(double);
      std::size_t actual_bytes = 0U;
      for (std::size_t i = 0U; i < allocation_count; ++i)
        if (reinterpret_cast<std::uintptr_t>(allocations[i].base) == owner_base)
          actual_bytes = allocations[i].bytes;
      const auto required = static_cast<std::size_t>(mg->interior.x);
      // This field is a raw linear bundle. Per-level ghost/alignment is
      // already in its declared length; only final 64-byte rounding is allowed.
      passed &= mg->interior.y == 1 && mg->interior.z == 1 &&
                mg->span_doubles >= required &&
                mg->span_doubles - required < 8U && owner != nullptr &&
                actual_bytes == owner->span_doubles * sizeof(double);
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::cout << "MG capacity rank=" << rank << " bundle_doubles=" << required
                << " outer_doubles=" << mg->span_doubles
                << " owner_allocated_bytes=" << actual_bytes << '\n';
      CartesianGeometryPlan geometry;
      MeshPatch patch;
      Status bounds = CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, model.mesh, {}, geometry, patch);
      MgWorkspaceRequirements requirements;
      MgHierarchyPolicy policy;
      policy.cycle = MgCycleKind::f_cycle;
      if (bounds)
        bounds = make_mg_workspace_requirements(MPI_COMM_WORLD, geometry, patch,
                                                policy, 1U, requirements);
      FieldView bundle;
      bundle.base = reinterpret_cast<double*>(plan.mg_storage_address());
      bundle.interior = mg->interior;
      bundle.components = 1U;
      bundle.stride_y = bundle.stride_z = bundle.component_stride = required;
      bundle.revision = bundle.storage_identity = bundle.revision_domain = 1U;
      MgWorkspace workspace;
      if (bounds) bounds = MgWorkspace::bind(requirements, bundle, workspace);
      passed &= bounds && requirements.total_doubles == required;
      const auto begin = plan.mg_storage_address();
      const auto end = begin + required * sizeof(double);
      std::uintptr_t previous_end = begin;
      for (std::size_t level = 0; bounds && level < workspace.level_count();
           ++level) {
        for (auto slot :
             {MgWorkspaceSlot::solution, MgWorkspaceSlot::rhs,
              MgWorkspaceSlot::residual, MgWorkspaceSlot::temporary}) {
          const FieldView view = workspace.level(level, slot);
          const auto raw =
              reinterpret_cast<std::uintptr_t>(view.base) -
              (1U + view.stride_y + view.stride_z) * sizeof(double);
          const auto last = reinterpret_cast<std::uintptr_t>(view.base) +
                            (view.interior.x + view.stride_y * view.interior.y +
                             view.stride_z * view.interior.z + 1U) *
                                sizeof(double);
          passed &= view.ghosts.x == 1 && view.ghosts.y == 1 &&
                    view.ghosts.z == 1 && raw >= previous_end && last <= end &&
                    last > raw;
          previous_end = last;
        }
      }
      std::cout << "MG bounded level views rank=" << rank
                << " levels=" << workspace.level_count()
                << " slots=" << requirements.slots_per_level
                << " valid=" << passed << '\n';
    }
    // Exercise the real bound MG levels, not just their storage description.
    ProductDriver driver;
    Status status = compiled;
    if (status)
      status = ProductDriver::create(MPI_COMM_WORLD, std::move(plan), driver);
    if (status) status = driver.initialize({});
    DriverStepReport step;
    if (status) status = driver.advance({1.0, 1.0, 1.0, 1.0, 1.0}, step);
    passed &= status && step.accepted;
  }
  int okay = passed ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!okay)
    std::cerr << "FAIL MG bundle storage exceeds its declared raw capacity\n";
  MPI_Finalize();
  return okay ? 0 : 1;
}
