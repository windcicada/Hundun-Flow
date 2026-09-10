// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_combustion.hpp"
#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_io.hpp"
#include "hundun/v04_portable.hpp"

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
  CouplingKind coupling{CouplingKind::piso};
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
  bool unity_lewis_enthalpy{};
  bool immersed{};
  IbmReconstructionAudit ibm_boundary_reconstruction{};
  IbmReconstructionAudit ibm_surface_reconstruction{};
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

// An explicit source-control continuation witness. This interface transports
// evidence; it never infers a fold or a departure branch from endpoint roots.
struct ProductTcrFoldQuery {
  std::uint64_t global_cell{};
  portable::Revision history_revision{}, input_revision{};
  bool initialized{};
  double accepted_eta{}, accepted_rate_ratio{}, eta{}, rate_ratio{};
  int accepted_branch_sign{};
};
struct ProductTcrFoldEvidence {
  bool supplied{};
  PlanFingerprint source_identity{};
  std::uint64_t global_cell{};
  portable::Revision history_revision{}, input_revision{};
  double eta{}, rate_ratio{};
  int departure_sign{};
};
class ProductTcrFoldProvider {
public:
  virtual ~ProductTcrFoldProvider() = default;
  // Stable content/algorithm identity, included in the product/restart
  // contract.
  virtual PlanFingerprint fingerprint() const noexcept = 0;
  // Pure, bounded, allocation-free query of immutable evidence. Providers must
  // not advance an independent clock/history. The native transaction owns all
  // accepted branch/fold history; repeated attempts query the same evidence.
  virtual Status query(const ProductTcrFoldQuery &,
                       ProductTcrFoldEvidence &) const noexcept = 0;
};

// Borrowed, exclusive-lane providers. Their owners must outlive the compiled
// plan and its driver. Case identity and species/thermo binding are validated
// collectively before a provider can supply a product source.
struct ProductCouplingBindings {
  portable::GasQueryProvider* gas_query{};
  portable::GasAdvanceProvider *gas_advance{};
  const combustion::ChemistryIdentity *chemistry_identity{};
  const ProductTcrFoldProvider *tcr_fold{};
};

class ProductCompiler {
 public:
  static Status compile(MPI_Comm communicator,
                        const ValidatedModel& model,
                        const std::filesystem::path& case_root,
                        CompiledCasePlan& out,
                        ProductCouplingBindings coupling = {}) noexcept;
  // Bind one validated source plan for explicit transport-history recovery.
  static Status compile_transport_restart(MPI_Comm communicator,
      const ValidatedModel& source, const std::filesystem::path& source_root,
      const ValidatedModel& target, const std::filesystem::path& target_root,
      CompiledCasePlan& out) noexcept;
};

}  // namespace hundun::v04
