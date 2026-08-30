// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_execution.hpp"
#include "hundun/v04_mesh.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04 {

enum class CpuKernelVariant : std::uint8_t { scalar, avx2, avx512 };

struct CpuExecutionRequest {
  std::uint32_t threads_per_rank{};
  bool pure_mpi{};
  bool bind_threads{true};
  Span<const std::int32_t> core_ids{};
};

struct CpuExecutionDiagnostics {
  std::uint32_t allowed_core_count{};
  std::uint32_t numa_node_count{};
  std::uint32_t local_rank_count{};
  std::uint32_t recommended_local_rank_count{};
  bool duplicate_core_warning{};
  bool affinity_constrained{};
};

class CpuExecutionPlan {
 public:
  CpuExecutionPlan() noexcept = default;
  CpuExecutionPlan(const CpuExecutionPlan&) = delete;
  CpuExecutionPlan& operator=(const CpuExecutionPlan&) = delete;
  CpuExecutionPlan(CpuExecutionPlan&&) noexcept = default;
  CpuExecutionPlan& operator=(CpuExecutionPlan&&) noexcept = default;

  static Status compile(MPI_Comm communicator, CpuExecutionRequest request,
                        CpuExecutionPlan& out) noexcept;

  std::size_t worker_count() const noexcept { return core_ids_.size(); }
  bool pure_mpi() const noexcept { return pure_mpi_; }
  bool binds_threads() const noexcept { return bind_threads_; }
  CpuKernelVariant kernel_variant() const noexcept { return kernel_variant_; }
  Int3 tile_shape() const noexcept { return tile_shape_; }
  CpuExecutionDiagnostics diagnostics() const noexcept { return diagnostics_; }
  Span<const std::int32_t> core_ids() const noexcept {
    return {core_ids_.data(), core_ids_.size()};
  }
  PlanFingerprint semantic_fingerprint() const noexcept {
    if (core_ids_.empty()) return 0U;
    std::uint64_t hash = UINT64_C(1469598103934665603);
    const auto mix = [&](std::uint64_t value) {
      hash ^= value;
      hash *= UINT64_C(1099511628211);
    };
    // Physical core ids are launcher placement evidence, not portable plan
    // semantics.  The frozen execution identity covers the worker topology
    // and selected kernel/tile policy.
    mix(static_cast<std::uint64_t>(core_ids_.size()));
    mix(pure_mpi_ ? 1U : 0U);
    mix(bind_threads_ ? 1U : 0U);
    mix(static_cast<std::uint8_t>(kernel_variant_));
    mix(static_cast<std::uint32_t>(tile_shape_.x));
    mix(static_cast<std::uint32_t>(tile_shape_.y));
    mix(static_cast<std::uint32_t>(tile_shape_.z));
    return hash == 0U ? 1U : hash;
  }

 private:
  friend class CpuThreadTeam;
  std::vector<std::int32_t> core_ids_;
  CpuKernelVariant kernel_variant_{CpuKernelVariant::scalar};
  Int3 tile_shape_{32, 4, 4};
  CpuExecutionDiagnostics diagnostics_{};
  bool pure_mpi_{true};
  bool bind_threads_{};
};

using CpuTask = void (*)(std::size_t worker, void* context) noexcept;

class CpuThreadTeam {
 public:
  CpuThreadTeam() noexcept = default;
  ~CpuThreadTeam() noexcept;
  CpuThreadTeam(const CpuThreadTeam&) = delete;
  CpuThreadTeam& operator=(const CpuThreadTeam&) = delete;
  CpuThreadTeam(CpuThreadTeam&& other) noexcept;
  CpuThreadTeam& operator=(CpuThreadTeam&& other) noexcept;

  static Status create(const CpuExecutionPlan& plan,
                       CpuThreadTeam& out) noexcept;
  Status run(CpuTask task, void* context) noexcept;
  std::size_t worker_count() const noexcept;

 private:
  struct Impl;
  explicit CpuThreadTeam(Impl* implementation) noexcept
      : implementation_(implementation) {}
  void release() noexcept;
  Impl* implementation_{};
};

struct HaloFieldSpec {
  FieldId field{};
  std::uint8_t width{};
  std::uint8_t components{};
};

struct HaloTopology {
  bool periodic_x{};
  bool periodic_y{};
  bool periodic_z{};
};

class HaloTicket {
 public:
  bool active() const noexcept { return active_; }
  StageId stage() const noexcept { return stage_; }

 private:
  friend class HaloEngine;
  std::uint64_t engine_identity_{};
  std::uint64_t generation_{};
  StageId stage_{};
  bool active_{};
};

struct HaloPlanStats {
  std::size_t transport_peer_count{};
  std::size_t local_peer_count{};
  std::size_t persistent_request_count{};
  std::uint64_t send_capacity_doubles{};
  std::uint64_t receive_capacity_doubles{};
  std::uint64_t maximum_messages_per_exchange{};
  std::uint64_t maximum_bytes_per_exchange{};
  std::uintptr_t request_storage_address{};
  std::uintptr_t send_storage_address{};
  std::uintptr_t receive_storage_address{};
  int maximum_tag{};
};

struct HaloRuntimeCounters {
  std::uint64_t begin_calls{};
  std::uint64_t finish_calls{};
  std::uint64_t messages_started{};
  std::uint64_t bytes_packed{};
  std::uint64_t bytes_unpacked{};
};

class HaloEngine {
 public:
  HaloEngine() noexcept = default;
  ~HaloEngine() noexcept;
  HaloEngine(const HaloEngine&) = delete;
  HaloEngine& operator=(const HaloEngine&) = delete;
  HaloEngine(HaloEngine&& other) noexcept;
  HaloEngine& operator=(HaloEngine&& other) noexcept;

  Status reserve(MPI_Comm communicator, const MeshPatch& patch,
                 Span<const HaloFieldSpec> fields) noexcept;
  Status reserve(MPI_Comm communicator, const MeshPatch& patch,
                 Span<const HaloFieldSpec> fields,
                 HaloTopology topology) noexcept;
  Status begin(StageId stage, Span<const FieldView> fields,
               HaloTicket& ticket) noexcept;
  Status begin(StageId stage, Span<const FieldView> fields,
               Status prerequisite, HaloTicket& ticket) noexcept;
  Status finish(HaloTicket& ticket, Span<FieldView> fields) noexcept;
  Status validate_contract(MPI_Comm communicator, const MeshPatch& patch,
                           Span<const HaloFieldSpec> fields,
                           HaloTopology topology) const noexcept;

  RevisionToken ghost_revision(FieldId field) const noexcept;
  int lowest_failing_rank() const noexcept;
  std::uintptr_t instance_identity() const noexcept {
    return reinterpret_cast<std::uintptr_t>(implementation_);
  }
  bool ready() const noexcept;
  bool active() const noexcept;
  HaloPlanStats plan_stats() const noexcept;
  HaloRuntimeCounters runtime_counters() const noexcept;

 private:
  struct Impl;
  void release() noexcept;
  Impl* implementation_{};
};

}  // namespace hundun::v04
