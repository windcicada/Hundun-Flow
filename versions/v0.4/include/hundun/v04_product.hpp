// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_io.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>

namespace hundun::v04 {

enum class ProductFreezePhase : std::uint8_t {
  geometry_and_decomposition,
  capability_registration,
  logical_analysis,
  schema_and_allocation,
  plan_instantiation,
  numeric_capacity,
  communication_binding,
  view_and_graph_binding,
  validation,
  sealed
};

struct PlanSummary {
  Int3 global_cells{};
  Int3 local_cells{};
  std::size_t field_count{};
  std::size_t arena_doubles{};
  std::size_t graph_stage_count{};
  std::size_t graph_node_count{};
  std::size_t maximum_workspace_bytes{};
  std::size_t service_staging_bytes{};
  std::uint8_t pressure_correctors{};
  double pressure_absolute_tolerance{};
  double pressure_relative_tolerance{};
  std::uint32_t pressure_maximum_iterations{};
  std::uint32_t pressure_true_residual_interval{};
  std::uint32_t pressure_krylov_restart{};
  double terminal_eos_tolerance{};
  double terminal_continuity_tolerance{};
  double terminal_closed_mass_tolerance{};
  double terminal_gauge_tolerance{};
  bool immersed{};
  bool exact_numeric_certified{};
  bool preconditioner_setup_certified{};
  bool sealed{};
};

class CompiledCasePlan {
 public:
  CompiledCasePlan() noexcept = default;
  ~CompiledCasePlan() noexcept;
  CompiledCasePlan(const CompiledCasePlan&) = delete;
  CompiledCasePlan& operator=(const CompiledCasePlan&) = delete;
  CompiledCasePlan(CompiledCasePlan&&) noexcept;
  CompiledCasePlan& operator=(CompiledCasePlan&&) noexcept;

  PlanFingerprint fingerprint() const noexcept;
  PlanSummary summary() const noexcept;
  Span<const ProductFreezePhase> freeze_order() const noexcept;
  const FieldSchema* field_schema() const noexcept;
  const ArenaLayout* arena_layout() const noexcept;
  const FrozenExecutionGraph* execution_graph() const noexcept;
  const IoServicePlan* io_services() const noexcept;
  const ContributionPlan* contribution_plan() const noexcept;
  const CpuExecutionPlan* cpu_execution_plan() const noexcept;
  PlanFingerprint cpu_plan_fingerprint() const noexcept;
  PlanFingerprint stl_fingerprint() const noexcept;
  std::uintptr_t state_storage_address() const noexcept;
  std::uintptr_t krylov_storage_address() const noexcept;
  std::uintptr_t mg_storage_address() const noexcept;

 private:
  friend class ProductCompiler;
  friend class ProductDriver;
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

class ProductCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const ValidatedModel& model,
                        const std::filesystem::path& case_root,
                        CompiledCasePlan& out) noexcept;
};

}  // namespace hundun::v04
