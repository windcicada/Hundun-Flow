// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

#include "solver_krylov_test_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
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

namespace {

using namespace hundun::v04;

constexpr std::int32_t kGlobalCells = 67;
constexpr RevisionToken kExecutionRevision = 7001U;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local, MPI_Comm communicator) {
  const int value = local ? 1 : 0;
  int result = 0;
  const int reduced =
      MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, communicator);
  return reduced == MPI_SUCCESS && result != 0;
}

struct Partition {
  std::int32_t begin{};
  std::int32_t cells{};
  int rank{};
  int size{};
};

Partition partition(int rank, int size) noexcept {
  const std::int32_t base = kGlobalCells / size;
  const std::int32_t remainder = kGlobalCells % size;
  return {rank * base + std::min(rank, static_cast<int>(remainder)),
          base + (rank < remainder ? 1 : 0), rank, size};
}

struct SolutionSnapshot {
  std::array<double, static_cast<std::size_t>(kGlobalCells)> values{};
  std::int32_t cells{};
};

SolutionSnapshot snapshot_solution(const FieldView& solution) noexcept {
  SolutionSnapshot result;
  result.cells = solution.interior.x;
  for (std::int32_t cell = 0; cell < solution.interior.x; ++cell) {
    result.values[static_cast<std::size_t>(cell)] =
        solution.unchecked({cell, 0, 0}, 0U);
  }
  return result;
}

bool same_solution(const FieldView& solution,
                   const SolutionSnapshot& snapshot) noexcept {
  if (solution.interior.x != snapshot.cells) {
    return false;
  }
  for (std::int32_t cell = 0; cell < solution.interior.x; ++cell) {
    if (solution.unchecked({cell, 0, 0}, 0U) !=
        snapshot.values[static_cast<std::size_t>(cell)]) {
      return false;
    }
  }
  return true;
}

bool same_field_values(const FieldView& left, const FieldView& right) noexcept {
  if (left.interior.x != right.interior.x ||
      left.interior.y != right.interior.y ||
      left.interior.z != right.interior.z) {
    return false;
  }
  for (std::int32_t z = 0; z < left.interior.z; ++z) {
    for (std::int32_t y = 0; y < left.interior.y; ++y) {
      for (std::int32_t x = 0; x < left.interior.x; ++x) {
        if (left.unchecked({x, y, z}, 0U) !=
            right.unchecked({x, y, z}, 0U)) {
          return false;
        }
      }
    }
  }
  return true;
}

bool same_reduction_work(const LinearReductionCounters& left,
                         const LinearReductionCounters& right) noexcept {
  return left.calls == right.calls && left.scalars == right.scalars &&
         left.logical_bytes == right.logical_bytes &&
         left.tree_messages == right.tree_messages &&
         left.blocking_operations == right.blocking_operations;
}

LinearReductionCounters reduction_delta(
    const LinearReductionCounters& after,
    const LinearReductionCounters& before) noexcept {
  LinearReductionCounters result;
  result.calls = after.calls - before.calls;
  result.scalars = after.scalars - before.scalars;
  result.logical_bytes = after.logical_bytes - before.logical_bytes;
  result.tree_messages = after.tree_messages - before.tree_messages;
  result.blocking_operations =
      after.blocking_operations - before.blocking_operations;
  result.wall_nanoseconds = after.wall_nanoseconds - before.wall_nanoseconds;
  return result;
}

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

OwnedField make_field(FieldId id, std::int32_t local_cells,
                      std::uint8_t ghosts, std::uint8_t components,
                      RevisionToken revision, StorageIdentity storage_id,
                      RevisionDomainIdentity domain) {
  OwnedField field;
  const std::size_t padded =
      static_cast<std::size_t>(local_cells + 2 * ghosts);
  const std::size_t component_stride = padded;
  field.storage.assign(component_stride * components, 0.0);
  field.view.base = field.storage.data() + ghosts;
  field.view.interior = {local_cells, 1, 1};
  field.view.ghosts = {ghosts, 0, 0};
  field.view.components = components;
  field.view.stride_y = padded;
  field.view.stride_z = padded;
  field.view.component_stride = component_stride;
  field.view.field = id;
  field.view.revision = revision;
  field.view.storage_identity = storage_id;
  field.view.revision_domain = domain;
  return field;
}

LinearIdentity identity(std::uint64_t salt) noexcept {
  return {11U + salt, 12U + salt, 13U + salt, 14U + salt,
          15U + salt};
}

struct FlexibleApplicationTrace {
  std::array<double, static_cast<std::size_t>(kGlobalCells)> expected{};
  std::uint32_t observed_operator_inputs{};
  std::uint32_t mismatched_operator_inputs{};
  std::uint32_t overwritten_outputs{};
  double previous_temporal_factor{};
  bool pending{};
  bool saw_spatial_variation{};
  bool saw_temporal_variation{};
};

class TridiagonalOperator final : public LinearOperator {
 public:
  TridiagonalOperator(MPI_Comm communicator, Partition local,
                      LinearIdentity expected, double lower, double diagonal,
                      double upper, bool spd, int failing_rank = -1,
                      std::uint32_t failing_call = 0U,
                      FlexibleApplicationTrace* flexible_trace = nullptr) noexcept
      : communicator_(communicator),
        local_(local),
        expected_(expected),
        lower_(lower),
        diagonal_(diagonal),
        upper_(upper),
        spd_(spd),
        failing_rank_(failing_rank),
        failing_call_(failing_call),
        flexible_trace_(flexible_trace) {}

  LinearOperatorCertificate certificate() const noexcept override {
    LinearOperatorCertificate result;
    result.identity = expected_;
    result.collective_fingerprint = 0x12345678U;
    result.local_shape = {local_.cells, 1, 1};
    result.operator_class =
        spd_ ? LinearOperatorClass::spd : LinearOperatorClass::nonsymmetric;
    return result;
  }

  Status apply(FieldView x, FieldView y) const noexcept override {
    ++calls_;
    if (flexible_trace_ != nullptr && flexible_trace_->pending) {
      bool matches = true;
      for (std::int32_t cell = 0; cell < local_.cells; ++cell) {
        matches =
            x.unchecked({cell, 0, 0}, 0U) ==
                flexible_trace_->expected[static_cast<std::size_t>(cell)] &&
            matches;
      }
      if (matches) {
        ++flexible_trace_->observed_operator_inputs;
      } else {
        ++flexible_trace_->mismatched_operator_inputs;
      }
      flexible_trace_->pending = false;
    }
    double left = 0.0;
    double right = 0.0;
    MPI_Request requests[4]{MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                            MPI_REQUEST_NULL, MPI_REQUEST_NULL};
    int request_count = 0;
    bool communication_ok = true;
    if (local_.rank > 0) {
      int mpi_status = MPI_Irecv(&left, 1, MPI_DOUBLE, local_.rank - 1, 1501,
                                 communicator_, &requests[request_count]);
      communication_ok = mpi_status == MPI_SUCCESS && communication_ok;
      request_count += mpi_status == MPI_SUCCESS ? 1 : 0;
      mpi_status = MPI_Isend(&x.unchecked({0, 0, 0}, 0U), 1, MPI_DOUBLE,
                             local_.rank - 1, 1502, communicator_,
                             &requests[request_count]);
      communication_ok = mpi_status == MPI_SUCCESS && communication_ok;
      request_count += mpi_status == MPI_SUCCESS ? 1 : 0;
    }
    if (local_.rank + 1 < local_.size) {
      int mpi_status = MPI_Irecv(&right, 1, MPI_DOUBLE, local_.rank + 1, 1502,
                                 communicator_, &requests[request_count]);
      communication_ok = mpi_status == MPI_SUCCESS && communication_ok;
      request_count += mpi_status == MPI_SUCCESS ? 1 : 0;
      mpi_status = MPI_Isend(
          &x.unchecked({local_.cells - 1, 0, 0}, 0U), 1, MPI_DOUBLE,
          local_.rank + 1, 1501, communicator_, &requests[request_count]);
      communication_ok = mpi_status == MPI_SUCCESS && communication_ok;
      request_count += mpi_status == MPI_SUCCESS ? 1 : 0;
    }
    const int wait_status =
        request_count == 0
            ? MPI_SUCCESS
            : MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
    if (!communication_ok || wait_status != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, 801U};
    }
    for (std::int32_t cell = 0; cell < local_.cells; ++cell) {
      const std::int32_t global = local_.begin + cell;
      const double west = global == 0
                              ? 0.0
                              : (cell == 0
                                     ? left
                                     : x.unchecked({cell - 1, 0, 0}, 0U));
      const double east = global + 1 == kGlobalCells
                              ? 0.0
                              : (cell + 1 == local_.cells
                                     ? right
                                     : x.unchecked({cell + 1, 0, 0}, 0U));
      y.unchecked({cell, 0, 0}, 0U) =
          lower_ * west + diagonal_ * x.unchecked({cell, 0, 0}, 0U) +
          upper_ * east;
    }
    if (local_.rank == failing_rank_ && calls_ == failing_call_) {
      return {StatusCode::numerical_failure, 802U};
    }
    return {};
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  MPI_Comm communicator_{MPI_COMM_NULL};
  Partition local_{};
  LinearIdentity expected_{};
  double lower_{};
  double diagonal_{};
  double upper_{};
  bool spd_{};
  int failing_rank_{-1};
  std::uint32_t failing_call_{};
  FlexibleApplicationTrace* flexible_trace_{};
  mutable std::uint32_t calls_{};
};

// Exact identity on the solvable subspace and an exact null row at global
// cell zero.  A correction supported only on cell zero is therefore a finite,
// nonzero solution-space direction whose current A2 image is bitwise zero.
class NullspaceIdentityOperator final : public LinearOperator {
 public:
  NullspaceIdentityOperator(LinearIdentity expected, Partition local) noexcept
      : expected_(expected), local_(local) {}

  LinearOperatorCertificate certificate() const noexcept override {
    return {expected_, 0x12345678U, {local_.cells, 1, 1},
            LinearOperatorClass::nonsymmetric};
  }

  Status apply(FieldView input, FieldView output) const noexcept override {
    ++calls_;
    for (std::int32_t cell = 0; cell < local_.cells; ++cell) {
      output.unchecked({cell, 0, 0}, 0U) =
          local_.begin + cell == 0
              ? 0.0
              : input.unchecked({cell, 0, 0}, 0U);
    }
    return {};
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  LinearIdentity expected_{};
  Partition local_{};
  mutable std::uint32_t calls_{};
};

class ScalingPreconditioner final : public LinearPreconditioner {
 public:
  ScalingPreconditioner(LinearIdentity expected, double scale, bool spd,
                        bool variable = false,
                        std::int32_t global_begin = 0, int local_rank = 0,
                        int failing_rank = -1,
                        std::uint32_t failing_call = 0U,
                        FlexibleApplicationTrace* flexible_trace = nullptr,
                        LinearPreconditionerStatusScope status_scope =
                            LinearPreconditionerStatusScope::rank_local) noexcept
      : expected_(expected),
        scale_(scale),
        global_begin_(global_begin),
        local_rank_(local_rank),
        failing_rank_(failing_rank),
        failing_call_(failing_call),
        flexible_trace_(flexible_trace),
        spd_(spd),
        variable_(variable),
        status_scope_(status_scope) {}

  LinearPreconditionerCertificate certificate() const noexcept override {
    return {expected_, 0x87654321U,
            spd_ ? LinearPreconditionerClass::fixed_spd
                 : (variable_ ? LinearPreconditionerClass::flexible
                              : LinearPreconditionerClass::fixed_general),
            status_scope_};
  }

  Status apply(ConstFieldView residual, FieldView output,
               std::uint32_t iteration) noexcept override {
    ++calls_;
    const double temporal_factor =
        variable_ ? 0.82 + 0.09 * static_cast<double>(iteration % 4U) : 1.0;
    if (flexible_trace_ != nullptr) {
      if (flexible_trace_->pending) {
        ++flexible_trace_->overwritten_outputs;
      }
      if (calls_ > 1U &&
          temporal_factor != flexible_trace_->previous_temporal_factor) {
        flexible_trace_->saw_temporal_variation = true;
      }
      flexible_trace_->previous_temporal_factor = temporal_factor;
      flexible_trace_->pending = true;
    }
    for (std::int32_t cell = 0; cell < residual.interior.x; ++cell) {
      const std::int32_t global = global_begin_ + cell;
      const double spatial_factor =
          variable_
              ? 0.76 + 0.04 * static_cast<double>((global * 5) % 13)
              : 1.0;
      const double value = scale_ * temporal_factor * spatial_factor *
                           residual.unchecked({cell, 0, 0}, 0U);
      output.unchecked({cell, 0, 0}, 0U) = value;
      if (flexible_trace_ != nullptr) {
        flexible_trace_->expected[static_cast<std::size_t>(cell)] = value;
        if (cell > 0) {
          const std::int32_t previous_global = global - 1;
          const double previous_spatial =
              0.76 +
              0.04 * static_cast<double>((previous_global * 5) % 13);
          flexible_trace_->saw_spatial_variation =
              spatial_factor != previous_spatial ||
              flexible_trace_->saw_spatial_variation;
        }
      }
    }
    if (local_rank_ == failing_rank_ && calls_ == failing_call_) {
      return {StatusCode::numerical_failure, 803U};
    }
    return {};
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  LinearIdentity expected_{};
  double scale_{};
  std::int32_t global_begin_{};
  int local_rank_{};
  int failing_rank_{-1};
  std::uint32_t failing_call_{};
  FlexibleApplicationTrace* flexible_trace_{};
  bool spd_{};
  bool variable_{};
  LinearPreconditionerStatusScope status_scope_{
      LinearPreconditionerStatusScope::rank_local};
  std::uint32_t calls_{};
};

// This fixture deliberately performs its own status agreement outside the
// solver reduction engine on success.  A failing callback additionally uses
// the engine consensus so the production lowest-failing-rank authority is
// preserved for the failure-path test.  The direct MPI agreement models the
// stronger callback contract while leaving the solver counter able to expose
// removal of the redundant outer consensus on successful columns.
class CollectiveStatusPreconditioner final : public LinearPreconditioner {
 public:
  CollectiveStatusPreconditioner(
      MPI_Comm communicator, ReductionEngine& reductions,
      LinearIdentity expected, double scale, Partition local,
      LinearPreconditionerStatusScope status_scope,
      int failing_rank = -1, std::uint32_t failing_call = 0U) noexcept
      : communicator_(communicator),
        reductions_(&reductions),
        expected_(expected),
        scale_(scale),
        local_(local),
        status_scope_(status_scope),
        failing_rank_(failing_rank),
        failing_call_(failing_call) {}

  LinearPreconditionerCertificate certificate() const noexcept override {
    return {expected_, 0x87654321U,
            LinearPreconditionerClass::fixed_general,
            status_scope_};
  }

  Status apply(ConstFieldView residual, FieldView output,
               std::uint32_t) noexcept override {
    ++calls_;
    for (std::int32_t cell = 0; cell < residual.interior.x; ++cell) {
      output.unchecked({cell, 0, 0}, 0U) =
          scale_ * residual.unchecked({cell, 0, 0}, 0U);
    }
    const int local_failure =
        local_.rank == failing_rank_ && calls_ == failing_call_ ? 1 : 0;
    int global_failure = 0;
    if (MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                      communicator_) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, 804U};
    }
    if (global_failure != 0) {
      return reductions_->consensus(
          local_failure != 0
              ? Status{StatusCode::numerical_failure, 805U}
              : Status{});
    }
    return {};
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  MPI_Comm communicator_{MPI_COMM_NULL};
  ReductionEngine* reductions_{};
  LinearIdentity expected_{};
  double scale_{};
  Partition local_{};
  LinearPreconditionerStatusScope status_scope_{
      LinearPreconditionerStatusScope::rank_local};
  int failing_rank_{-1};
  std::uint32_t failing_call_{};
  std::uint32_t calls_{};
};

class PreparedBatchPreconditioner final : public LinearPreconditioner {
 public:
  PreparedBatchPreconditioner(
      ReductionEngine& reductions, LinearIdentity expected, double scale,
      int rank,
      LinearPreconditionerStatusScope status_scope =
                     LinearPreconditionerStatusScope::collective,
      LinearPreconditionerApplyLifecycle lifecycle =
          LinearPreconditionerApplyLifecycle::prepared_batch,
      int failing_prepare_rank = -1, int failing_apply_rank = -1,
      std::uint32_t failing_apply_call = 0U) noexcept
      : reductions_(&reductions),
        expected_(expected),
        scale_(scale),
        rank_(rank),
        status_scope_(status_scope),
        lifecycle_(lifecycle),
        failing_prepare_rank_(failing_prepare_rank),
        failing_apply_rank_(failing_apply_rank),
        failing_apply_call_(failing_apply_call) {}

  LinearPreconditionerCertificate certificate() const noexcept override {
    return {expected_, 0x87654321U, LinearPreconditionerClass::flexible,
            status_scope_, lifecycle_};
  }

  Status prepare_batch(const LinearPreconditionerBatchDescriptor& descriptor,
                       LinearPreconditionerBatchTicket& ticket) noexcept
      override {
    ++prepare_calls_;
    descriptor_workspace_ = descriptor.workspace;
    descriptor_shape_ = descriptor.shape;
    descriptor_slot_count_ = descriptor.slot_count;
    descriptor_maximum_applications_ = descriptor.maximum_applications;
    if (descriptor.workspace == nullptr || descriptor.slot_count == 0U ||
        descriptor.maximum_applications == 0U ||
        (rank_ == failing_prepare_rank_)) {
      return {StatusCode::invalid_plan, 880U};
    }
    issue_batch_ticket(ticket, this, descriptor);
    return {};
  }

  Status apply(ConstFieldView input, FieldView output,
               std::uint32_t iteration) noexcept override {
    ++direct_calls_;
    return apply_common(input, output, iteration);
  }

  Status apply_prepared(ConstFieldView input, FieldView output,
                        std::uint32_t iteration,
                        const LinearPreconditionerBatchTicket&) noexcept
      override {
    ++prepared_calls_;
    return apply_common(input, output, iteration);
  }

  std::uint32_t prepare_calls() const noexcept { return prepare_calls_; }
  std::uint32_t direct_calls() const noexcept { return direct_calls_; }
  std::uint32_t prepared_calls() const noexcept { return prepared_calls_; }
  const SolverWorkspace* descriptor_workspace() const noexcept {
    return descriptor_workspace_;
  }
  Int3 descriptor_shape() const noexcept { return descriptor_shape_; }
  std::uint8_t descriptor_slot_count() const noexcept {
    return descriptor_slot_count_;
  }
  std::uint32_t descriptor_maximum_applications() const noexcept {
    return descriptor_maximum_applications_;
  }

 private:
  Status apply_common(ConstFieldView input, FieldView output,
                      std::uint32_t) noexcept {
    ++application_calls_;
    for (std::int32_t cell = 0; cell < input.interior.x; ++cell) {
      output.unchecked({cell, 0, 0}, 0U) =
          scale_ * input.unchecked({cell, 0, 0}, 0U);
    }
    const Status local = rank_ == failing_apply_rank_ &&
                                 application_calls_ == failing_apply_call_
                             ? Status{StatusCode::numerical_failure, 882U}
                             : Status{};
    return reductions_->consensus(local);
  }

  ReductionEngine* reductions_{};
  LinearIdentity expected_{};
  double scale_{};
  int rank_{};
  LinearPreconditionerStatusScope status_scope_{
      LinearPreconditionerStatusScope::collective};
  LinearPreconditionerApplyLifecycle lifecycle_{
      LinearPreconditionerApplyLifecycle::prepared_batch};
  int failing_prepare_rank_{-1};
  int failing_apply_rank_{-1};
  std::uint32_t failing_apply_call_{};
  const SolverWorkspace* descriptor_workspace_{};
  Int3 descriptor_shape_{};
  std::uint8_t descriptor_slot_count_{};
  std::uint32_t descriptor_maximum_applications_{};
  std::uint32_t application_calls_{};
  std::uint32_t prepare_calls_{};
  std::uint32_t direct_calls_{};
  std::uint32_t prepared_calls_{};
};

class RankVariantControlOperator final : public LinearOperator {
 public:
  RankVariantControlOperator(LinearIdentity expected, Partition local) noexcept
      : expected_(expected), local_(local) {}

  LinearOperatorCertificate certificate() const noexcept override {
    return {expected_, 0x12345678U, {local_.cells, 1, 1},
            LinearOperatorClass::spd};
  }

  Status apply(FieldView, FieldView) const noexcept override {
    ++calls_;
    return {};
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  LinearIdentity expected_{};
  Partition local_{};
  mutable std::uint32_t calls_{};
};

class MismatchedProvenanceOperator final : public LinearOperator {
 public:
  MismatchedProvenanceOperator(LinearIdentity expected, Partition local)
      noexcept
      : expected_(expected), local_(local) {}

  LinearOperatorCertificate certificate() const noexcept override {
    return {expected_, 0x12345678U, {local_.cells, 1, 1},
            LinearOperatorClass::nonsymmetric};
  }

  Status apply(FieldView, FieldView) const noexcept override {
    ++calls_;
    return {StatusCode::numerical_failure, 806U};
  }

  LinearOperatorFailureProvenance failure_provenance() const
      noexcept override {
    return {{StatusCode::numerical_failure, 807U},
            LinearOperatorStatusScope::collective, local_.size - 1};
  }

  std::uint32_t calls() const noexcept { return calls_; }

 private:
  LinearIdentity expected_{};
  Partition local_{};
  mutable std::uint32_t calls_{};
};

struct SolveFixture {
  Partition local{};
  OwnedField rhs;
  OwnedField solution;
  OwnedField vectors;
  OwnedField scalars;
  LinearWorkspaceRequirements requirements{};
  SolverWorkspace workspace;
  ReductionEngine reductions;
  LinearIdentity expected{};
};

bool initialize_fixture(MPI_Comm communicator, LinearAlgorithm algorithm,
                        std::uint32_t restart, SolveFixture& fixture) {
  int rank = 0;
  int size = 0;
  const bool communicator_ready =
      MPI_Comm_rank(communicator, &rank) == MPI_SUCCESS &&
      MPI_Comm_size(communicator, &size) == MPI_SUCCESS && size > 0;
  if (!communicator_ready) {
    return false;
  }
  bool local_ready = true;
  try {
    fixture.local = partition(rank, size);
    fixture.expected = identity(static_cast<std::uint64_t>(algorithm));
    local_ready = static_cast<bool>(make_linear_workspace_requirements(
        algorithm, {fixture.local.cells, 1, 1}, 0U, restart,
        ReductionMode::mpi_allreduce, kExecutionRevision,
        fixture.requirements));
    if (local_ready) {
      const StorageIdentity storage =
          9001U + static_cast<std::uint64_t>(algorithm);
      const RevisionDomainIdentity domain = 9101U;
      fixture.rhs = make_field(1U, fixture.local.cells, 0U, 1U, 1U,
                               storage + 100U, domain);
      fixture.solution = make_field(2U, fixture.local.cells, 0U, 1U, 2U,
                                    storage + 200U, domain);
      fixture.vectors = make_field(
          3U, fixture.local.cells, 0U, fixture.requirements.vector_slots, 3U,
          storage, domain);
      fixture.scalars = make_field(
          4U, static_cast<std::int32_t>(fixture.requirements.scalar_doubles),
          0U, 1U, 4U, storage, domain);
    }
  } catch (...) {
    local_ready = false;
  }
  LinearLifecycleCounters counters{};
  if (local_ready) {
    local_ready = static_cast<bool>(SolverWorkspace::bind(
        fixture.requirements, fixture.vectors.view, fixture.scalars.view,
        fixture.workspace, &counters));
    if (local_ready) {
      fixture.expected.workspace = fixture.workspace.fingerprint();
    }
  }
  const bool every_rank_ready = all_true(local_ready, communicator);
  if (!every_rank_ready) {
    return false;
  }
  return static_cast<bool>(ReductionEngine::compile(
      communicator, ReductionMode::mpi_allreduce,
      fixture.requirements.reduction_capacity, fixture.reductions));
}

double exact_solution(std::int32_t global) noexcept {
  constexpr double pi = 3.141592653589793238462643383279502884;
  return std::sin(pi * static_cast<double>(global + 1) /
                  static_cast<double>(kGlobalCells + 1));
}

void fill_system(SolveFixture& fixture, double lower, double diagonal,
                 double upper, double initial = 0.0) noexcept {
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    const std::int32_t global = fixture.local.begin + cell;
    const double west = global == 0 ? 0.0 : exact_solution(global - 1);
    const double east = global + 1 == kGlobalCells
                            ? 0.0
                            : exact_solution(global + 1);
    fixture.rhs.view.unchecked({cell, 0, 0}, 0U) =
        lower * west + diagonal * exact_solution(global) + upper * east;
    fixture.solution.view.unchecked({cell, 0, 0}, 0U) = initial;
  }
}

struct ResidualOracle {
  double residual{std::numeric_limits<double>::infinity()};
  double rhs_norm{std::numeric_limits<double>::infinity()};
  bool communication_ok{};
};

double scaled_l2(
    const std::array<double, static_cast<std::size_t>(kGlobalCells)>& values)
    noexcept {
  double scale = 0.0;
  double sum_squares = 1.0;
  for (double value : values) {
    const double magnitude = std::abs(value);
    if (magnitude == 0.0) {
      continue;
    }
    if (scale < magnitude) {
      const double ratio = scale / magnitude;
      sum_squares = 1.0 + sum_squares * ratio * ratio;
      scale = magnitude;
    } else {
      const double ratio = magnitude / scale;
      sum_squares += ratio * ratio;
    }
  }
  return scale == 0.0 ? 0.0 : scale * std::sqrt(sum_squares);
}

ResidualOracle independent_true_residual(const SolveFixture& fixture,
                                         MPI_Comm communicator, double lower,
                                         double diagonal,
                                         double upper) noexcept {
  ResidualOracle result;
  if (fixture.local.size <= 0 || fixture.local.size > 4) {
    return result;
  }
  std::array<int, 4U> counts{};
  std::array<int, 4U> displacements{};
  for (int other = 0; other < fixture.local.size; ++other) {
    const Partition other_partition = partition(other, fixture.local.size);
    counts[static_cast<std::size_t>(other)] = other_partition.cells;
    displacements[static_cast<std::size_t>(other)] = other_partition.begin;
  }
  std::array<double, static_cast<std::size_t>(kGlobalCells)> solution{};
  std::array<double, static_cast<std::size_t>(kGlobalCells)> rhs{};
  const int solution_gather = MPI_Allgatherv(
      fixture.solution.view.base, fixture.local.cells, MPI_DOUBLE,
      solution.data(), counts.data(), displacements.data(), MPI_DOUBLE,
      communicator);
  const int rhs_gather = MPI_Allgatherv(
      fixture.rhs.view.base, fixture.local.cells, MPI_DOUBLE, rhs.data(),
      counts.data(), displacements.data(), MPI_DOUBLE, communicator);
  result.communication_ok =
      solution_gather == MPI_SUCCESS && rhs_gather == MPI_SUCCESS;
  if (!result.communication_ok) {
    return result;
  }
  std::array<double, static_cast<std::size_t>(kGlobalCells)> residual{};
  for (std::int32_t global = 0; global < kGlobalCells; ++global) {
    const double west = global == 0
                            ? 0.0
                            : solution[static_cast<std::size_t>(global - 1)];
    const double east =
        global + 1 == kGlobalCells
            ? 0.0
            : solution[static_cast<std::size_t>(global + 1)];
    const double applied =
        lower * west + diagonal * solution[static_cast<std::size_t>(global)] +
        upper * east;
    residual[static_cast<std::size_t>(global)] =
        rhs[static_cast<std::size_t>(global)] - applied;
  }
  result.residual = scaled_l2(residual);
  result.rhs_norm = scaled_l2(rhs);
  return result;
}

bool residual_report_matches(double reported,
                             const ResidualOracle& oracle) noexcept {
  const double tolerance =
      512.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, oracle.rhs_norm);
  return oracle.communication_ok && std::isfinite(reported) &&
         std::isfinite(oracle.residual) &&
         std::abs(reported - oracle.residual) <= tolerance;
}

bool residual_is_accepted(const ResidualOracle& oracle,
                          LinearSolveControl selected) noexcept {
  return oracle.communication_ok && std::isfinite(oracle.residual) &&
         oracle.residual <=
             std::max(selected.absolute_tolerance,
                      selected.relative_tolerance * oracle.rhs_norm);
}

bool same_call_count(std::uint32_t local, MPI_Comm communicator) noexcept {
  std::uint32_t minimum = 0U;
  std::uint32_t maximum = 0U;
  const int minimum_status = MPI_Allreduce(
      &local, &minimum, 1, MPI_UINT32_T, MPI_MIN, communicator);
  const int maximum_status = MPI_Allreduce(
      &local, &maximum, 1, MPI_UINT32_T, MPI_MAX, communicator);
  return minimum_status == MPI_SUCCESS && maximum_status == MPI_SUCCESS &&
         minimum == maximum;
}

double global_error(const SolveFixture& fixture, MPI_Comm communicator) {
  double local = 0.0;
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    const double difference =
        fixture.solution.view.unchecked({cell, 0, 0}, 0U) -
        exact_solution(fixture.local.begin + cell);
    local = std::max(local, std::abs(difference));
  }
  double result = 0.0;
  return MPI_Allreduce(&local, &result, 1, MPI_DOUBLE, MPI_MAX,
                       communicator) == MPI_SUCCESS
             ? result
             : std::numeric_limits<double>::infinity();
}

bool finite_solution(const SolveFixture& fixture) noexcept {
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    if (!std::isfinite(
            fixture.solution.view.unchecked({cell, 0, 0}, 0U))) {
      return false;
    }
  }
  return true;
}

LinearSolveControl control(std::uint32_t maximum = 300U,
                           std::uint32_t restart = 0U) noexcept {
  LinearSolveControl result;
  result.absolute_tolerance = 1.0e-12;
  result.relative_tolerance = 1.0e-10;
  result.maximum_iterations = maximum;
  result.true_residual_interval = 8U;
  result.restart = restart;
  return result;
}

LinearSolveInvocation invocation(const SolveFixture& fixture,
                                 LinearSolveControl selected) noexcept {
  return {as_const(fixture.rhs.view), fixture.solution.view, fixture.expected,
          selected};
}

class RejectOnceConvergenceAudit final : public LinearConvergenceAudit {
 public:
  LinearConvergenceAuditCertificate certificate() const noexcept override {
    return {UINT64_C(0x7630346175646974)};
  }

  Status evaluate(ConstFieldView, ConstFieldView,
                  ReductionEngine& reductions,
                  LinearConvergenceAuditResult& result) noexcept override {
    const double local_metric = evaluations_++ == 0U ? 2.0 : 0.5;
    double global_metric = 0.0;
    const Status status = reductions.checked_max(
        {&local_metric, 1U}, {&global_metric, 1U});
    if (!status) return status;
    result = {global_metric <= 1.0, global_metric, 1.0};
    return {};
  }

 private:
  std::uint32_t evaluations_{};
};

class AcceptingConvergenceAudit final : public LinearConvergenceAudit {
 public:
  LinearConvergenceAuditCertificate certificate() const noexcept override {
    return {UINT64_C(0x7630346163636570)};
  }

  Status evaluate(ConstFieldView, ConstFieldView residual,
                  ReductionEngine& reductions,
                  LinearConvergenceAuditResult& result) noexcept override {
    ++evaluations_;
    double local_maximum = 0.0;
    Status local_status{};
    for (std::int32_t cell = 0; cell < residual.interior.x; ++cell) {
      const double value = residual.unchecked({cell, 0, 0}, 0U);
      if (!std::isfinite(value)) {
        local_status = {StatusCode::numerical_failure, 0x61636365U};
        break;
      }
      local_maximum = std::max(local_maximum, std::abs(value));
    }
    double global_maximum = 0.0;
    const Status reduced = reductions.checked_max(
        {&local_maximum, 1U}, {&global_maximum, 1U}, local_status);
    if (!reduced) return reduced;
    result = {true, global_maximum, 1.0};
    return {};
  }

  std::uint32_t evaluations() const noexcept { return evaluations_; }

 private:
  std::uint32_t evaluations_{};
};

bool test_pcg_spd_and_zero_rhs(MPI_Comm communicator, int rank) {
  SolveFixture fixture;
  bool passed = expect(initialize_fixture(communicator, LinearAlgorithm::pcg,
                                          0U, fixture),
                       rank, "PCG fixture initializes");
  if (!all_true(passed, communicator)) {
    return false;
  }
  fill_system(fixture, -1.0, 2.0, -1.0);
  TridiagonalOperator op(communicator, fixture.local, fixture.expected, -1.0,
                         2.0, -1.0, true);
  ScalingPreconditioner preconditioner(fixture.expected, 0.5, true);
  const LinearSolveControl selected = control();
  const LinearSolveResult result =
      solve_pcg(op, preconditioner, invocation(fixture, selected),
                fixture.workspace, fixture.reductions);
  const ResidualOracle oracle = independent_true_residual(
      fixture, communicator, -1.0, 2.0, -1.0);
  const double error = global_error(fixture, communicator);
  passed &= expect(result.status.code == StatusCode::ok &&
                       result.termination == LinearTermination::converged &&
                       result.iterations > 0U &&
                       residual_report_matches(result.final_true_residual,
                                               oracle) &&
                       residual_is_accepted(oracle, selected) &&
                       result.operator_applies >= result.iterations + 2U &&
                       finite_solution(fixture) && error < 1.0e-8,
                   rank,
                   "PCG solves certified SPD system with explicit true residual");

  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    fixture.rhs.view.unchecked({cell, 0, 0}, 0U) = 0.0;
    fixture.solution.view.unchecked({cell, 0, 0}, 0U) = 9.0;
  }
  const LinearSolveResult zero =
      solve_pcg(op, preconditioner, invocation(fixture, selected),
                fixture.workspace, fixture.reductions);
  const ResidualOracle zero_oracle = independent_true_residual(
      fixture, communicator, -1.0, 2.0, -1.0);
  bool zero_solution = true;
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    zero_solution = fixture.solution.view.unchecked({cell, 0, 0}, 0U) == 0.0 &&
                    zero_solution;
  }
  passed &= expect(zero.status.code == StatusCode::ok &&
                       zero.termination == LinearTermination::zero_rhs &&
                       zero.iterations == 0U && zero_solution &&
                       residual_report_matches(zero.final_true_residual,
                                               zero_oracle) &&
                       zero_oracle.residual == 0.0,
                   rank, "zero RHS returns exact zero solution without iteration");
  return all_true(passed, communicator);
}

bool test_nonsymmetric_solvers(MPI_Comm communicator, int rank) {
  bool passed = true;
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 12U, fixture);
    passed &= expect(initialized, rank, "FGMRES fixture initializes");
    if (!all_true(initialized, communicator)) {
      return false;
    }
    fill_system(fixture, -1.2, 3.0, -0.7);
    FlexibleApplicationTrace trace;
    TridiagonalOperator op(communicator, fixture.local, fixture.expected, -1.2,
                           3.0, -0.7, false, -1, 0U, &trace);
    ScalingPreconditioner flexible(
        fixture.expected, 1.0 / 3.0, false, true, fixture.local.begin,
        fixture.local.rank, -1, 0U, &trace);
    LinearSolveControl selected = control(250U, 12U);
    selected.true_residual_interval = 1U;
    const LinearSolveResult result = solve_fgmres(
        op, flexible, invocation(fixture, selected),
        fixture.workspace, fixture.reductions);
    const ResidualOracle oracle = independent_true_residual(
        fixture, communicator, -1.2, 3.0, -0.7);
    const double error = global_error(fixture, communicator);
    const bool calls_match = same_call_count(flexible.calls(), communicator);
    passed &= expect(result.status.code == StatusCode::ok &&
                         result.termination == LinearTermination::converged &&
                         result.iterations > 1U &&
                         result.iterations < 80U &&
                         residual_report_matches(result.final_true_residual,
                                                 oracle) &&
                         residual_is_accepted(oracle, selected) &&
                         finite_solution(fixture) && error < 1.0e-8 &&
                         flexible.calls() == result.iterations &&
                         flexible.calls() == result.preconditioner_applies &&
                         op.calls() == result.operator_applies &&
                         result.reduction_calls ==
                             4U + 3U * result.iterations &&
                         result.operator_applies ==
                             1U + 2U * result.iterations &&
                         result.norm_breakdown_restarts == 0U &&
                         trace.observed_operator_inputs == flexible.calls() &&
                         trace.mismatched_operator_inputs == 0U &&
                         trace.overwritten_outputs == 0U && !trace.pending &&
                         trace.saw_spatial_variation &&
                         trace.saw_temporal_variation && calls_match,
                     rank,
                     "true-residual audits preserve the flexible FGMRES Arnoldi cycle on a nonnormal operator");
  }
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 12U, fixture);
    passed &= expect(initialized, rank,
                     "ill-scaled FGMRES fixture initializes");
    if (!all_true(initialized, communicator)) {
      return false;
    }
    constexpr double scale = 1.0e8;
    fill_system(fixture, -1.2 * scale, 3.0 * scale, -0.7 * scale);
    TridiagonalOperator op(communicator, fixture.local, fixture.expected,
                           -1.2 * scale, 3.0 * scale, -0.7 * scale, false);
    ScalingPreconditioner fixed(fixture.expected, 1.0 / (3.0 * scale),
                                false);
    const LinearSolveControl selected = control(250U, 12U);
    const LinearSolveResult result = solve_fgmres(
        op, fixed, invocation(fixture, selected), fixture.workspace,
        fixture.reductions);
    const ResidualOracle oracle = independent_true_residual(
        fixture, communicator, -1.2 * scale, 3.0 * scale, -0.7 * scale);
    const double error = global_error(fixture, communicator);
    passed &= expect(result.status.code == StatusCode::ok &&
                         result.termination == LinearTermination::converged &&
                         result.iterations > 0U &&
                         residual_report_matches(result.final_true_residual,
                                                 oracle) &&
                         residual_is_accepted(oracle, selected) &&
                         result.norm_breakdown_restarts == 0U &&
                         finite_solution(fixture) && error < 1.0e-8,
                     rank,
                     "single-reduction FGMRES solves a 1e8-scaled nonnormal system");
  }
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::bicgstab, 0U, fixture);
    passed &= expect(initialized, rank, "BiCGSTAB fixture initializes");
    if (!all_true(initialized, communicator)) {
      return false;
    }
    fill_system(fixture, -1.2, 3.0, -0.7);
    TridiagonalOperator op(communicator, fixture.local, fixture.expected, -1.2,
                           3.0, -0.7, false);
    ScalingPreconditioner fixed(fixture.expected, 1.0 / 3.0, false);
    const LinearSolveControl selected = control(250U);
    const LinearSolveResult result = solve_bicgstab(
        op, fixed, invocation(fixture, selected), fixture.workspace,
        fixture.reductions);
    const ResidualOracle oracle = independent_true_residual(
        fixture, communicator, -1.2, 3.0, -0.7);
    const double error = global_error(fixture, communicator);
    passed &= expect(result.status.code == StatusCode::ok &&
                         result.termination == LinearTermination::converged &&
                         result.iterations > 0U &&
                         residual_report_matches(result.final_true_residual,
                                                 oracle) &&
                         residual_is_accepted(oracle, selected) &&
                         result.operator_applies >= result.iterations + 2U &&
                         fixed.calls() >= result.iterations &&
                         finite_solution(fixture) && error < 1.0e-8,
                     rank,
                     "BiCGSTAB solves a nonsymmetric advection-diffusion system");
  }
  return all_true(passed, communicator);
}

bool test_rejections_are_transactional(MPI_Comm communicator, int rank,
                                       int size) {
  SolveFixture fixture;
  const bool initialized = initialize_fixture(communicator, LinearAlgorithm::pcg,
                                              0U, fixture);
  bool passed = expect(initialized, rank,
                       "transactional failure fixture initializes");
  if (!all_true(initialized, communicator)) {
    return false;
  }
  fill_system(fixture, -1.2, 3.0, -0.7, 4.25);
  const SolutionSnapshot before = snapshot_solution(fixture.solution.view);
  TridiagonalOperator nonsymmetric(communicator, fixture.local,
                                   fixture.expected, -1.2, 3.0, -0.7, false);
  ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0, true);
  const LinearSolveResult invalid = solve_pcg(
      nonsymmetric, preconditioner, invocation(fixture, control()),
      fixture.workspace, fixture.reductions);
  passed &= expect(invalid.status.code == StatusCode::invalid_plan &&
                       nonsymmetric.calls() == 0U &&
                       preconditioner.calls() == 0U &&
                       same_solution(fixture.solution.view, before),
                   rank,
                   "PCG rejects a non-SPD operator before callbacks or publishing x");

  TridiagonalOperator alias_operator(communicator, fixture.local,
                                     fixture.expected, -1.0, 2.0, -1.0,
                                     true);
  ScalingPreconditioner alias_preconditioner(fixture.expected, 0.5, true);
  LinearSolveInvocation shifted_alias = invocation(fixture, control());
  shifted_alias.rhs = as_const(fixture.solution.view);
  ++shifted_alias.rhs.base;
  shifted_alias.rhs.storage_identity =
      fixture.solution.view.storage_identity + 100000U;
  const LinearSolveResult alias_result = solve_pcg(
      alias_operator, alias_preconditioner, shifted_alias, fixture.workspace,
      fixture.reductions);
  passed &= expect(
      alias_result.status.code == StatusCode::invalid_plan &&
          alias_operator.calls() == 0U && alias_preconditioner.calls() == 0U &&
          same_solution(fixture.solution.view, before),
      rank,
      "shifted overlapping caller views are rejected by address, independent of identity");

  LinearIdentity rank_local_identity = fixture.expected;
  if (rank == size - 1) {
    ++rank_local_identity.fingerprint;
  }
  TridiagonalOperator certificate_mismatch(
      communicator, fixture.local, rank_local_identity, -1.0, 2.0, -1.0,
      true);
  ScalingPreconditioner certificate_preconditioner(fixture.expected, 0.5,
                                                    true);
  const LinearSolveResult mismatched = solve_pcg(
      certificate_mismatch, certificate_preconditioner,
      invocation(fixture, control()), fixture.workspace, fixture.reductions);
  passed &= expect(mismatched.status.code == StatusCode::invalid_plan &&
                       certificate_mismatch.calls() == 0U &&
                       certificate_preconditioner.calls() == 0U &&
                       same_solution(fixture.solution.view, before),
                   rank,
                   "rank-local certificate mismatch is collective before callbacks");

  RankVariantControlOperator control_operator(fixture.expected,
                                               fixture.local);
  ScalingPreconditioner control_preconditioner(fixture.expected, 0.5, true);
  LinearSolveControl rank_variant_control = control();
  if (size > 1 && rank == size - 1) {
    ++rank_variant_control.maximum_iterations;
  }
  const LinearSolveResult control_mismatch = solve_pcg(
      control_operator, control_preconditioner,
      invocation(fixture, rank_variant_control), fixture.workspace,
      fixture.reductions);
  if (size > 1) {
    passed &= expect(
        control_mismatch.status.code == StatusCode::invalid_plan &&
            control_operator.calls() == 0U &&
            control_preconditioner.calls() == 0U &&
            same_solution(fixture.solution.view, before),
        rank,
        "rank-local valid control mismatch is rejected before callbacks");
  } else {
    passed &= expect(control_mismatch.status.code != StatusCode::invalid_plan,
                     rank,
                     "single-rank matching control passes collective contract validation");
  }

  TridiagonalOperator failing(communicator, fixture.local, fixture.expected,
                              -1.0, 2.0, -1.0, true, size - 1, 2U);
  const LinearSolveResult failure =
      solve_pcg(failing, preconditioner, invocation(fixture, control()),
                fixture.workspace, fixture.reductions);
  passed &= expect(failure.status.code == StatusCode::numerical_failure &&
                       failure.lowest_failing_rank == size - 1 &&
                       same_solution(fixture.solution.view, before) &&
                       finite_solution(fixture),
                   rank,
                   "rank-local operator failure is collective and leaves caller x unchanged");

  ScalingPreconditioner failing_preconditioner(
      fixture.expected, 0.5, true, false, fixture.local.begin,
      fixture.local.rank, size - 1, 2U);
  TridiagonalOperator healthy(communicator, fixture.local, fixture.expected,
                              -1.0, 2.0, -1.0, true);
  const LinearSolveResult preconditioner_failure = solve_pcg(
      healthy, failing_preconditioner, invocation(fixture, control()),
      fixture.workspace, fixture.reductions);
  passed &= expect(
      preconditioner_failure.status.code == StatusCode::numerical_failure &&
          preconditioner_failure.lowest_failing_rank == size - 1 &&
          same_solution(fixture.solution.view, before) &&
          finite_solution(fixture),
      rank,
      "rank-local preconditioner failure preserves collective progress and caller x");

  fill_system(fixture, -1.0, 2.0, -1.0);
  const LinearSolveControl selected = control();
  const LinearSolveResult recovered = solve_pcg(
      healthy, preconditioner, invocation(fixture, selected),
      fixture.workspace, fixture.reductions);
  const ResidualOracle recovered_oracle = independent_true_residual(
      fixture, communicator, -1.0, 2.0, -1.0);
  const double recovered_error = global_error(fixture, communicator);
  passed &= expect(recovered.status.code == StatusCode::ok &&
                       recovered.termination == LinearTermination::converged &&
                       residual_report_matches(recovered.final_true_residual,
                                               recovered_oracle) &&
                       residual_is_accepted(recovered_oracle, selected) &&
                       recovered_error < 1.0e-8,
                   rank,
                   "same reduction engine and workspace recover after local callback failures");
  return all_true(passed, communicator);
}

bool test_breakdown_maxiter_and_stale_identity(MPI_Comm communicator,
                                                int rank) {
  SolveFixture fixture;
  const bool initialized = initialize_fixture(
      communicator, LinearAlgorithm::bicgstab, 0U, fixture);
  bool passed = expect(initialized, rank, "breakdown fixture initializes");
  if (!all_true(initialized, communicator)) {
    return false;
  }
  fill_system(fixture, 0.0, 0.0, 0.0, -3.0);
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    fixture.rhs.view.unchecked({cell, 0, 0}, 0U) = 1.0;
  }
  TridiagonalOperator zero(communicator, fixture.local, fixture.expected, 0.0,
                           0.0, 0.0, false);
  ScalingPreconditioner fixed(fixture.expected, 1.0, false);
  const SolutionSnapshot before_breakdown =
      snapshot_solution(fixture.solution.view);
  const LinearSolveResult breakdown = solve_bicgstab(
      zero, fixed, invocation(fixture, control()), fixture.workspace,
      fixture.reductions);
  passed &= expect(breakdown.status.code == StatusCode::numerical_failure &&
                       breakdown.termination == LinearTermination::breakdown &&
                       same_solution(fixture.solution.view, before_breakdown) &&
                       finite_solution(fixture),
                   rank, "zero operator produces deterministic finite breakdown");

  fill_system(fixture, -1.2, 3.0, -0.7, 1.25);
  const SolutionSnapshot before_maximum =
      snapshot_solution(fixture.solution.view);
  TridiagonalOperator regular(communicator, fixture.local, fixture.expected,
                              -1.2, 3.0, -0.7, false);
  const LinearSolveResult maximum = solve_bicgstab(
      regular, fixed, invocation(fixture, control(1U)), fixture.workspace,
      fixture.reductions);
  passed &= expect(maximum.status.code == StatusCode::rejected_step &&
                       maximum.termination ==
                           LinearTermination::maximum_iterations &&
                       same_solution(fixture.solution.view, before_maximum) &&
                       finite_solution(fixture),
                   rank, "iteration cap cannot masquerade as convergence");

  const LinearIdentity live = fixture.expected;
  fixture.expected.fingerprint += 1U;
  const std::uint32_t calls_before_stale = regular.calls();
  const std::uint32_t preconditioner_calls_before_stale = fixed.calls();
  const SolutionSnapshot before_stale = snapshot_solution(fixture.solution.view);
  const LinearSolveResult stale = solve_bicgstab(
      regular, fixed, invocation(fixture, control()), fixture.workspace,
      fixture.reductions);
  passed &= expect(stale.status.code == StatusCode::invalid_plan &&
                       regular.calls() == calls_before_stale &&
                       fixed.calls() == preconditioner_calls_before_stale &&
                       same_solution(fixture.solution.view, before_stale),
                   rank, "stale invocation identity is rejected collectively");
  fixture.expected = live;
  return all_true(passed, communicator);
}

bool test_fgmres_norm_breakdown_lifecycle(MPI_Comm communicator, int rank) {
  bool passed = true;
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, fixture);
    passed &= expect(initialized, rank,
                     "first-column FGMRES breakdown fixture initializes");
    if (!all_true(initialized, communicator)) return false;
    fill_system(fixture, -1.2, 3.0, -0.7, 2.5);
    const SolutionSnapshot before = snapshot_solution(fixture.solution.view);
    TridiagonalOperator op(communicator, fixture.local, fixture.expected,
                           -1.2, 3.0, -0.7, false);
    ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0,
                                         false);
    LinearSolveControl selected = control(300U, 4U);
    selected.true_residual_interval = 1U;
    detail::force_single_reduction_fgmres_breakdown_for_test(1U);
    const LinearSolveResult result = solve_fgmres(
        op, preconditioner, invocation(fixture, selected),
        fixture.workspace, fixture.reductions);
    detail::force_single_reduction_fgmres_breakdown_for_test(0U);
    const ResidualOracle oracle = independent_true_residual(
        fixture, communicator, -1.2, 3.0, -0.7);
    const double error = global_error(fixture, communicator);
    passed &= expect(
        result.status.code == StatusCode::ok &&
            result.termination == LinearTermination::converged &&
            result.iterations > 1U && result.iterations < 80U &&
            result.norm_breakdown_restarts == 0U &&
            result.preconditioner_applies == result.iterations &&
            result.operator_applies == 1U + 2U * result.iterations &&
            result.reduction_calls ==
                4U + 3U * result.iterations +
                    3U * ((result.iterations + selected.restart - 1U) /
                          selected.restart) &&
            !same_solution(fixture.solution.view, before) &&
            residual_report_matches(result.final_true_residual, oracle) &&
            residual_is_accepted(oracle, selected) &&
            finite_solution(fixture) && error < 1.0e-8,
        rank,
        "first-column unsafe norm is explicitly recovered without a restart");
  }
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, fixture);
    passed &= expect(initialized, rank,
                     "bounded FGMRES breakdown-restart fixture initializes");
    if (!all_true(initialized, communicator)) return false;
    fill_system(fixture, -1.2, 3.0, -0.7, 2.5);
    const SolutionSnapshot before = snapshot_solution(fixture.solution.view);
    TridiagonalOperator op(communicator, fixture.local, fixture.expected,
                           -1.2, 3.0, -0.7, false);
    ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0,
                                         false);
    detail::force_single_reduction_fgmres_breakdown_for_test(2U);
    const LinearSolveResult result = solve_fgmres(
        op, preconditioner, invocation(fixture, control(4U, 4U)),
        fixture.workspace, fixture.reductions);
    detail::force_single_reduction_fgmres_breakdown_for_test(0U);
    passed &= expect(
            result.status.code == StatusCode::rejected_step &&
            result.termination == LinearTermination::maximum_iterations &&
            result.iterations == 4U &&
            result.norm_breakdown_restarts == 1U &&
            result.preconditioner_applies == 4U &&
            result.operator_applies == 7U &&
            result.reduction_calls == 12U &&
            same_solution(fixture.solution.view, before) &&
            finite_solution(fixture),
        rank,
        "later unsafe norm replaces the true residual and remains bounded");
  }
  return all_true(passed, communicator);
}

bool test_fgmres_rank_selective_failures(MPI_Comm communicator, int rank,
                                         int size) {
  bool passed = true;
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, fixture);
    passed &= expect(initialized, rank,
                     "FGMRES rank-local operator failure fixture initializes");
    if (!all_true(initialized, communicator)) return false;
    fill_system(fixture, -1.2, 3.0, -0.7, 3.25);
    const SolutionSnapshot before = snapshot_solution(fixture.solution.view);
    TridiagonalOperator failing(communicator, fixture.local, fixture.expected,
                                -1.2, 3.0, -0.7, false, size - 1, 2U);
    ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0,
                                         false);
    const LinearSolveResult result = solve_fgmres(
        failing, preconditioner, invocation(fixture, control(8U, 4U)),
        fixture.workspace, fixture.reductions);
    passed &= expect(result.status.code == StatusCode::numerical_failure &&
                         result.lowest_failing_rank == size - 1 &&
                         same_solution(fixture.solution.view, before) &&
                         finite_solution(fixture),
                     rank,
                     "rank-local FGMRES operator failure is collective and transactional");
  }
  {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, fixture);
    passed &= expect(
        initialized, rank,
        "FGMRES rank-local preconditioner failure fixture initializes");
    if (!all_true(initialized, communicator)) return false;
    fill_system(fixture, -1.2, 3.0, -0.7, 3.25);
    const SolutionSnapshot before = snapshot_solution(fixture.solution.view);
    TridiagonalOperator healthy(communicator, fixture.local, fixture.expected,
                                -1.2, 3.0, -0.7, false);
    ScalingPreconditioner failing(
        fixture.expected, 1.0 / 3.0, false, false, fixture.local.begin,
        fixture.local.rank, size - 1, 1U);
    const LinearSolveResult result = solve_fgmres(
        healthy, failing, invocation(fixture, control(8U, 4U)),
        fixture.workspace, fixture.reductions);
    passed &= expect(
        result.status.code == StatusCode::numerical_failure &&
            result.lowest_failing_rank == size - 1 &&
            same_solution(fixture.solution.view, before) &&
            finite_solution(fixture),
        rank,
        "rank-local FGMRES preconditioner failure is collective and transactional");
  }
  return all_true(passed, communicator);
}

bool test_operator_provenance_matching(MPI_Comm communicator, int rank) {
  SolveFixture fixture;
  const bool initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, fixture);
  bool passed = expect(initialized, rank,
                       "operator provenance matching fixture initializes");
  if (!all_true(initialized, communicator)) return false;
  fill_system(fixture, -1.2, 3.0, -0.7, 3.25);
  const SolutionSnapshot before = snapshot_solution(fixture.solution.view);
  MismatchedProvenanceOperator mismatched(fixture.expected, fixture.local);
  ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0, false);
  const LinearSolveResult result = solve_fgmres(
      mismatched, preconditioner, invocation(fixture, control(8U, 4U)),
      fixture.workspace, fixture.reductions);
  passed &= expect(
      result.status.code == StatusCode::numerical_failure &&
          result.status.detail == 806U &&
          result.termination == LinearTermination::operator_failure &&
          result.lowest_failing_rank == 0 && mismatched.calls() == 1U &&
          same_solution(fixture.solution.view, before),
      rank,
      "non-matching collective metadata falls back to rank-local consensus");
  return all_true(passed, communicator);
}

bool test_fgmres_collective_status_scope(MPI_Comm communicator, int rank,
                                          int size) {
  bool passed = true;
  const LinearSolveControl selected = control(250U, 12U);

  SolveFixture rank_local_fixture;
  const bool rank_local_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 12U, rank_local_fixture);
  passed &= expect(rank_local_initialized, rank,
                   "rank-local status-scope fixture initializes");
  if (!all_true(rank_local_initialized, communicator)) return false;
  fill_system(rank_local_fixture, -1.2, 3.0, -0.7);
  TridiagonalOperator rank_local_operator(
      communicator, rank_local_fixture.local, rank_local_fixture.expected,
      -1.2, 3.0, -0.7, false);
  CollectiveStatusPreconditioner rank_local_preconditioner(
      communicator, rank_local_fixture.reductions,
      rank_local_fixture.expected, 1.0 / 3.0, rank_local_fixture.local,
      LinearPreconditionerStatusScope::rank_local);
  passed &= expect(
      rank_local_preconditioner.certificate().status_scope ==
          LinearPreconditionerStatusScope::rank_local,
      rank, "rank-local preconditioner status scope remains the default");
  const LinearSolveResult rank_local_result = solve_fgmres(
      rank_local_operator, rank_local_preconditioner,
      invocation(rank_local_fixture, selected), rank_local_fixture.workspace,
      rank_local_fixture.reductions);

  SolveFixture collective_fixture;
  const bool collective_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 12U, collective_fixture);
  passed &= expect(collective_initialized, rank,
                   "collective status-scope fixture initializes");
  if (!all_true(collective_initialized, communicator)) return false;
  fill_system(collective_fixture, -1.2, 3.0, -0.7);
  TridiagonalOperator collective_operator(
      communicator, collective_fixture.local, collective_fixture.expected,
      -1.2, 3.0, -0.7, false);
  CollectiveStatusPreconditioner collective_preconditioner(
      communicator, collective_fixture.reductions, collective_fixture.expected,
      1.0 / 3.0, collective_fixture.local,
      LinearPreconditionerStatusScope::collective);
  passed &= expect(
      collective_preconditioner.certificate().status_scope ==
          LinearPreconditionerStatusScope::collective,
      rank, "collective preconditioner explicitly declares collective status");
  const LinearSolveResult collective_result = solve_fgmres(
      collective_operator, collective_preconditioner,
      invocation(collective_fixture, selected), collective_fixture.workspace,
      collective_fixture.reductions);
  const LinearReductionCounters rank_local_reduction_counters =
      rank_local_fixture.reductions.counters();
  const LinearReductionCounters collective_reduction_counters =
      collective_fixture.reductions.counters();
  const ResidualOracle collective_oracle = independent_true_residual(
      collective_fixture, communicator, -1.2, 3.0, -0.7);
  passed &= expect(
      rank_local_result.status.code == StatusCode::ok &&
          rank_local_result.termination == LinearTermination::converged &&
          collective_result.status.code == StatusCode::ok &&
          collective_result.termination == LinearTermination::converged &&
          collective_result.iterations == rank_local_result.iterations &&
          collective_preconditioner.calls() == collective_result.iterations &&
          collective_reduction_counters.blocking_operations +
                  collective_result.iterations ==
              rank_local_reduction_counters.blocking_operations &&
          residual_report_matches(collective_result.final_true_residual,
                                  collective_oracle) &&
          residual_is_accepted(collective_oracle, selected) &&
          finite_solution(collective_fixture),
      rank,
      "collective FGMRES status removes exactly one counted outer reduction per successful column");

  SolveFixture rank_local_failure_fixture;
  const bool rank_local_failure_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, rank_local_failure_fixture);
  passed &= expect(rank_local_failure_initialized, rank,
                   "rank-local collective-status failure fixture initializes");
  if (!all_true(rank_local_failure_initialized, communicator)) return false;
  fill_system(rank_local_failure_fixture, -1.2, 3.0, -0.7, 2.5);
  const SolutionSnapshot rank_local_before =
      snapshot_solution(rank_local_failure_fixture.solution.view);
  TridiagonalOperator rank_local_failure_operator(
      communicator, rank_local_failure_fixture.local,
      rank_local_failure_fixture.expected, -1.2, 3.0, -0.7, false);
  CollectiveStatusPreconditioner rank_local_failure_preconditioner(
      communicator, rank_local_failure_fixture.reductions,
      rank_local_failure_fixture.expected, 1.0 / 3.0,
      rank_local_failure_fixture.local,
      LinearPreconditionerStatusScope::rank_local, size - 1, 1U);
  const LinearSolveResult rank_local_failure = solve_fgmres(
      rank_local_failure_operator, rank_local_failure_preconditioner,
      invocation(rank_local_failure_fixture, control(8U, 4U)),
      rank_local_failure_fixture.workspace,
      rank_local_failure_fixture.reductions);

  SolveFixture collective_failure_fixture;
  const bool collective_failure_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, collective_failure_fixture);
  passed &= expect(collective_failure_initialized, rank,
                   "collective status failure fixture initializes");
  if (!all_true(collective_failure_initialized, communicator)) return false;
  fill_system(collective_failure_fixture, -1.2, 3.0, -0.7, 2.5);
  const SolutionSnapshot collective_before =
      snapshot_solution(collective_failure_fixture.solution.view);
  TridiagonalOperator collective_failure_operator(
      communicator, collective_failure_fixture.local,
      collective_failure_fixture.expected, -1.2, 3.0, -0.7, false);
  CollectiveStatusPreconditioner collective_failure_preconditioner(
      communicator, collective_failure_fixture.reductions,
      collective_failure_fixture.expected, 1.0 / 3.0,
      collective_failure_fixture.local,
      LinearPreconditionerStatusScope::collective, size - 1, 1U);
  const LinearSolveResult collective_failure = solve_fgmres(
      collective_failure_operator, collective_failure_preconditioner,
      invocation(collective_failure_fixture, control(8U, 4U)),
      collective_failure_fixture.workspace,
      collective_failure_fixture.reductions);
  const LinearReductionCounters rank_local_failure_reductions =
      rank_local_failure_fixture.reductions.counters();
  const LinearReductionCounters collective_failure_reductions =
      collective_failure_fixture.reductions.counters();
  passed &= expect(
      rank_local_failure.status.code == StatusCode::numerical_failure &&
          collective_failure.status.code == StatusCode::numerical_failure &&
          rank_local_failure.termination ==
              LinearTermination::preconditioner_failure &&
          collective_failure.termination ==
              LinearTermination::preconditioner_failure &&
          rank_local_failure.lowest_failing_rank == 0 &&
          collective_failure.lowest_failing_rank == size - 1 &&
          rank_local_failure_reductions.blocking_operations ==
              collective_failure_reductions.blocking_operations + 1U &&
          same_solution(rank_local_failure_fixture.solution.view,
                        rank_local_before) &&
          same_solution(collective_failure_fixture.solution.view,
                        collective_before) &&
          finite_solution(collective_failure_fixture),
      rank,
      "collective preconditioner failure preserves its internal lowest rank, avoids one redundant reduction, and leaves caller x unchanged");

  SolveFixture invalid_fixture;
  const bool invalid_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, invalid_fixture);
  passed &= expect(invalid_initialized, rank,
                   "invalid status-scope fixture initializes");
  if (!all_true(invalid_initialized, communicator)) return false;
  fill_system(invalid_fixture, -1.2, 3.0, -0.7, 1.5);
  const SolutionSnapshot invalid_before =
      snapshot_solution(invalid_fixture.solution.view);
  TridiagonalOperator invalid_operator(
      communicator, invalid_fixture.local, invalid_fixture.expected, -1.2,
      3.0, -0.7, false);
  ScalingPreconditioner invalid_preconditioner(
      invalid_fixture.expected, 1.0 / 3.0, false, false,
      invalid_fixture.local.begin, invalid_fixture.local.rank, -1, 0U,
      nullptr, static_cast<LinearPreconditionerStatusScope>(255U));
  const LinearSolveResult invalid = solve_fgmres(
      invalid_operator, invalid_preconditioner,
      invocation(invalid_fixture, control(8U, 4U)), invalid_fixture.workspace,
      invalid_fixture.reductions);
  passed &= expect(
      invalid.status.code == StatusCode::invalid_plan &&
          invalid_operator.calls() == 0U &&
          invalid_preconditioner.calls() == 0U &&
          same_solution(invalid_fixture.solution.view, invalid_before),
      rank, "invalid preconditioner status scope is rejected before callbacks");

  if (size > 1) {
    SolveFixture mismatch_fixture;
    const bool mismatch_initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, mismatch_fixture);
    passed &= expect(mismatch_initialized, rank,
                     "rank-mismatched status-scope fixture initializes");
    if (!all_true(mismatch_initialized, communicator)) return false;
    fill_system(mismatch_fixture, -1.2, 3.0, -0.7, 1.5);
    const SolutionSnapshot mismatch_before =
        snapshot_solution(mismatch_fixture.solution.view);
    TridiagonalOperator mismatch_operator(
        communicator, mismatch_fixture.local, mismatch_fixture.expected, -1.2,
        3.0, -0.7, false);
    const LinearPreconditionerStatusScope mismatch_scope =
        rank == size - 1 ? LinearPreconditionerStatusScope::collective
                         : LinearPreconditionerStatusScope::rank_local;
    ScalingPreconditioner mismatch_preconditioner(
        mismatch_fixture.expected, 1.0 / 3.0, false, false,
        mismatch_fixture.local.begin, mismatch_fixture.local.rank, -1, 0U,
        nullptr, mismatch_scope);
    const LinearSolveResult mismatch = solve_fgmres(
        mismatch_operator, mismatch_preconditioner,
        invocation(mismatch_fixture, control(8U, 4U)), mismatch_fixture.workspace,
        mismatch_fixture.reductions);
    passed &= expect(
        mismatch.status.code == StatusCode::invalid_plan &&
            mismatch_operator.calls() == 0U &&
            mismatch_preconditioner.calls() == 0U &&
            same_solution(mismatch_fixture.solution.view, mismatch_before),
        rank,
        "rank-mismatched preconditioner status scope rejects before callbacks");
  }
  return all_true(passed, communicator);
}

bool test_fgmres_prepared_batch_lifecycle(MPI_Comm communicator, int rank,
                                           int size) {
  bool passed = true;
  const LinearSolveControl selected = control(250U, 12U);

  SolveFixture baseline_fixture;
  const bool baseline_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, selected.restart,
      baseline_fixture);
  passed &= expect(baseline_initialized, rank,
                   "prepared-batch baseline fixture initializes");
  if (!all_true(baseline_initialized, communicator)) return false;
  fill_system(baseline_fixture, -1.2, 3.0, -0.7);
  TridiagonalOperator baseline_operator(
      communicator, baseline_fixture.local, baseline_fixture.expected, -1.2,
      3.0, -0.7, false);
  PreparedBatchPreconditioner baseline_preconditioner(
      baseline_fixture.reductions, baseline_fixture.expected, 1.0 / 3.0,
      rank, LinearPreconditionerStatusScope::collective,
      LinearPreconditionerApplyLifecycle::per_call_checked);
  const LinearSolveResult baseline = solve_fgmres(
      baseline_operator, baseline_preconditioner,
      invocation(baseline_fixture, selected), baseline_fixture.workspace,
      baseline_fixture.reductions);

  SolveFixture prepared_fixture;
  const bool prepared_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, selected.restart,
      prepared_fixture);
  passed &= expect(prepared_initialized, rank,
                   "prepared-batch fixture initializes");
  if (!all_true(prepared_initialized, communicator)) return false;
  fill_system(prepared_fixture, -1.2, 3.0, -0.7);
  TridiagonalOperator prepared_operator(
      communicator, prepared_fixture.local, prepared_fixture.expected, -1.2,
      3.0, -0.7, false);
  PreparedBatchPreconditioner prepared_preconditioner(
      prepared_fixture.reductions, prepared_fixture.expected, 1.0 / 3.0,
      rank);
  const LinearSolveResult prepared = solve_fgmres(
      prepared_operator, prepared_preconditioner,
      invocation(prepared_fixture, selected), prepared_fixture.workspace,
      prepared_fixture.reductions);
  const double prepared_error = global_error(prepared_fixture, communicator);
  const double baseline_error = global_error(baseline_fixture, communicator);
  const LinearReductionCounters baseline_counters =
      baseline_fixture.reductions.counters();
  const LinearReductionCounters prepared_counters =
      prepared_fixture.reductions.counters();
  passed &= expect(
      baseline.status.code == StatusCode::ok &&
      prepared.status.code == StatusCode::ok &&
          baseline.termination == LinearTermination::converged &&
          prepared.termination == LinearTermination::converged &&
          baseline_preconditioner.certificate().status_scope ==
              LinearPreconditionerStatusScope::collective &&
          baseline_preconditioner.certificate().apply_lifecycle ==
              LinearPreconditionerApplyLifecycle::per_call_checked &&
          prepared_preconditioner.certificate().status_scope ==
              LinearPreconditionerStatusScope::collective &&
          prepared_preconditioner.certificate().apply_lifecycle ==
              LinearPreconditionerApplyLifecycle::prepared_batch &&
          prepared.iterations == baseline.iterations &&
          prepared.operator_applies == baseline.operator_applies &&
          prepared.preconditioner_applies == baseline.preconditioner_applies &&
          prepared.reduction_calls == baseline.reduction_calls &&
          baseline_preconditioner.prepare_calls() == 0U &&
          baseline_preconditioner.direct_calls() == baseline.iterations &&
          baseline_preconditioner.prepared_calls() == 0U &&
          prepared_preconditioner.prepare_calls() == 1U &&
          prepared_preconditioner.direct_calls() == 0U &&
          prepared_preconditioner.prepared_calls() == prepared.iterations &&
          prepared_preconditioner.descriptor_workspace() ==
              &prepared_fixture.workspace &&
          prepared_preconditioner.descriptor_shape().x ==
              prepared_fixture.local.cells &&
          prepared_preconditioner.descriptor_shape().y == 1 &&
          prepared_preconditioner.descriptor_shape().z == 1 &&
          prepared_preconditioner.descriptor_slot_count() == selected.restart &&
          prepared_preconditioner.descriptor_maximum_applications() ==
              selected.maximum_iterations &&
          same_reduction_work(baseline_counters, prepared_counters) &&
          same_field_values(baseline_fixture.solution.view,
                            prepared_fixture.solution.view) &&
          std::abs(prepared.final_true_residual -
                   baseline.final_true_residual) < 1.0e-14 &&
          prepared_error < 1.0e-8 && baseline_error < 1.0e-8,
      rank,
      "per-call and prepared FGMRES have identical internal consensus work and numerical result; only dispatch differs");

  SolveFixture hot_failure_fixture;
  const bool hot_failure_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, hot_failure_fixture);
  passed &= expect(hot_failure_initialized, rank,
                   "prepared hot-failure fixture initializes");
  if (!all_true(hot_failure_initialized, communicator)) return false;
  fill_system(hot_failure_fixture, -1.2, 3.0, -0.7, 2.25);
  const SolutionSnapshot hot_failure_before =
      snapshot_solution(hot_failure_fixture.solution.view);
  TridiagonalOperator hot_failure_operator(
      communicator, hot_failure_fixture.local, hot_failure_fixture.expected,
      -1.2, 3.0, -0.7, false);
  PreparedBatchPreconditioner hot_failure_preconditioner(
      hot_failure_fixture.reductions, hot_failure_fixture.expected, 1.0 / 3.0,
      rank, LinearPreconditionerStatusScope::collective,
      LinearPreconditionerApplyLifecycle::prepared_batch, -1, size - 1, 1U);
  const LinearSolveResult hot_failure = solve_fgmres(
      hot_failure_operator, hot_failure_preconditioner,
      invocation(hot_failure_fixture, control(8U, 4U)),
      hot_failure_fixture.workspace, hot_failure_fixture.reductions);
  passed &= expect(
      hot_failure.status.code == StatusCode::numerical_failure &&
          hot_failure.termination == LinearTermination::preconditioner_failure &&
          hot_failure.lowest_failing_rank == size - 1 &&
          hot_failure.preconditioner_applies == 1U &&
          hot_failure.operator_applies == 1U &&
          hot_failure_operator.calls() == 1U &&
          hot_failure_preconditioner.prepare_calls() == 1U &&
          hot_failure_preconditioner.direct_calls() == 0U &&
          hot_failure_preconditioner.prepared_calls() == 1U &&
          same_solution(hot_failure_fixture.solution.view, hot_failure_before) &&
          finite_solution(hot_failure_fixture),
      rank,
      "rank-selective prepared hot failure agrees internally, reports the lowest rank, enters no operator after the initial residual, and preserves caller x");

  SolveFixture rank_local_fixture;
  const bool rank_local_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, rank_local_fixture);
  passed &= expect(rank_local_initialized, rank,
                   "prepared rank-local rejection fixture initializes");
  if (!all_true(rank_local_initialized, communicator)) return false;
  fill_system(rank_local_fixture, -1.2, 3.0, -0.7, 2.25);
  const SolutionSnapshot rank_local_before =
      snapshot_solution(rank_local_fixture.solution.view);
  TridiagonalOperator rank_local_operator(
      communicator, rank_local_fixture.local, rank_local_fixture.expected, -1.2,
      3.0, -0.7, false);
  PreparedBatchPreconditioner rank_local_preconditioner(
      rank_local_fixture.reductions, rank_local_fixture.expected, 1.0 / 3.0,
      rank,
      LinearPreconditionerStatusScope::rank_local);
  const LinearSolveResult rank_local = solve_fgmres(
      rank_local_operator, rank_local_preconditioner,
      invocation(rank_local_fixture, control(8U, 4U)),
      rank_local_fixture.workspace, rank_local_fixture.reductions);
  passed &= expect(
      rank_local.status.code == StatusCode::invalid_plan &&
          rank_local_operator.calls() == 0U &&
          rank_local_preconditioner.prepare_calls() == 0U &&
          rank_local_preconditioner.direct_calls() == 0U &&
          rank_local_preconditioner.prepared_calls() == 0U &&
          same_solution(rank_local_fixture.solution.view, rank_local_before),
      rank, "prepared lifecycle rejects rank-local status before callbacks");

  SolveFixture invalid_fixture;
  const bool invalid_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, invalid_fixture);
  passed &= expect(invalid_initialized, rank,
                   "invalid prepared lifecycle fixture initializes");
  if (!all_true(invalid_initialized, communicator)) return false;
  fill_system(invalid_fixture, -1.2, 3.0, -0.7, 2.25);
  TridiagonalOperator invalid_operator(
      communicator, invalid_fixture.local, invalid_fixture.expected, -1.2,
      3.0, -0.7, false);
  PreparedBatchPreconditioner invalid_preconditioner(
      invalid_fixture.reductions, invalid_fixture.expected, 1.0 / 3.0, rank,
      LinearPreconditionerStatusScope::collective,
      static_cast<LinearPreconditionerApplyLifecycle>(255U));
  const LinearSolveResult invalid = solve_fgmres(
      invalid_operator, invalid_preconditioner,
      invocation(invalid_fixture, control(8U, 4U)), invalid_fixture.workspace,
      invalid_fixture.reductions);
  passed &= expect(
      invalid.status.code == StatusCode::invalid_plan &&
          invalid_operator.calls() == 0U &&
          invalid_preconditioner.prepare_calls() == 0U &&
          invalid_preconditioner.direct_calls() == 0U &&
          invalid_preconditioner.prepared_calls() == 0U,
      rank, "invalid prepared lifecycle rejects before callbacks");

  if (size > 1) {
    SolveFixture mismatch_fixture;
    const bool mismatch_initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, mismatch_fixture);
    passed &= expect(mismatch_initialized, rank,
                     "prepared lifecycle mismatch fixture initializes");
    if (!all_true(mismatch_initialized, communicator)) return false;
    fill_system(mismatch_fixture, -1.2, 3.0, -0.7, 2.25);
    TridiagonalOperator mismatch_operator(
        communicator, mismatch_fixture.local, mismatch_fixture.expected, -1.2,
        3.0, -0.7, false);
    const LinearPreconditionerApplyLifecycle mismatch_lifecycle =
        rank == size - 1
            ? LinearPreconditionerApplyLifecycle::prepared_batch
            : LinearPreconditionerApplyLifecycle::per_call_checked;
    PreparedBatchPreconditioner mismatch_preconditioner(
        mismatch_fixture.reductions, mismatch_fixture.expected, 1.0 / 3.0,
        rank,
        LinearPreconditionerStatusScope::collective, mismatch_lifecycle);
    const LinearSolveResult mismatch = solve_fgmres(
        mismatch_operator, mismatch_preconditioner,
        invocation(mismatch_fixture, control(8U, 4U)),
        mismatch_fixture.workspace, mismatch_fixture.reductions);
    passed &= expect(
        mismatch.status.code == StatusCode::invalid_plan &&
            mismatch_operator.calls() == 0U &&
            mismatch_preconditioner.prepare_calls() == 0U &&
            mismatch_preconditioner.direct_calls() == 0U &&
            mismatch_preconditioner.prepared_calls() == 0U,
        rank, "prepared lifecycle mismatch rejects before every callback");
  }

  SolveFixture failure_fixture;
  const bool failure_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, failure_fixture);
  passed &= expect(failure_initialized, rank,
                   "prepared cold-failure fixture initializes");
  if (!all_true(failure_initialized, communicator)) return false;
  fill_system(failure_fixture, -1.2, 3.0, -0.7, 2.25);
  const SolutionSnapshot failure_before =
      snapshot_solution(failure_fixture.solution.view);
  TridiagonalOperator failure_operator(
      communicator, failure_fixture.local, failure_fixture.expected, -1.2,
      3.0, -0.7, false);
  PreparedBatchPreconditioner failure_preconditioner(
      failure_fixture.reductions, failure_fixture.expected, 1.0 / 3.0, rank,
      LinearPreconditionerStatusScope::collective,
      LinearPreconditionerApplyLifecycle::prepared_batch, size - 1);
  const LinearSolveResult failure = solve_fgmres(
      failure_operator, failure_preconditioner,
      invocation(failure_fixture, control(8U, 4U)), failure_fixture.workspace,
      failure_fixture.reductions);
  passed &= expect(
      failure.status.code == StatusCode::invalid_plan &&
          failure.lowest_failing_rank == size - 1 &&
          failure_operator.calls() == 0U &&
          failure_preconditioner.prepare_calls() == 1U &&
          failure_preconditioner.direct_calls() == 0U &&
          failure_preconditioner.prepared_calls() == 0U &&
          same_solution(failure_fixture.solution.view, failure_before),
      rank,
      "rank-selective cold prepare failure is collective and transactional");
  return all_true(passed, communicator);
}

bool test_fgmres_reduction_capacity(MPI_Comm communicator, int rank) {
  bool passed = true;
  for (const std::uint32_t restart : {2U, 64U}) {
    SolveFixture fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, restart, fixture);
    passed &= expect(initialized, rank,
                     "FGMRES reduction-capacity fixture initializes");
    if (!all_true(initialized, communicator)) return false;
    const std::size_t expected_capacity =
        restart + 1U < kLinearRecycleMaximumDirections
            ? kLinearRecycleMaximumDirections
            : restart + 1U;
    passed &= expect(
        fixture.requirements.reduction_capacity == expected_capacity, rank,
        "FGMRES combined dot/norm capacity also covers the recycle QR rank");
    fill_system(fixture, -1.2, 3.0, -0.7);
    TridiagonalOperator op(communicator, fixture.local, fixture.expected,
                           -1.2, 3.0, -0.7, false);
    ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0,
                                         false);
    LinearSolveControl selected = control(250U, restart);
    selected.true_residual_interval = 4U;
    const LinearSolveResult result = solve_fgmres(
        op, preconditioner, invocation(fixture, selected), fixture.workspace,
        fixture.reductions);
    const ResidualOracle oracle = independent_true_residual(
        fixture, communicator, -1.2, 3.0, -0.7);
    passed &= expect(result.status.code == StatusCode::ok &&
                         result.termination == LinearTermination::converged &&
                         residual_report_matches(result.final_true_residual,
                                                 oracle) &&
                         residual_is_accepted(oracle, selected),
                     rank,
                     "restart 2/64 combined reduction reaches true residual");
  }
  return all_true(passed, communicator);
}

bool test_fgmres_supplemental_convergence_audit(MPI_Comm communicator,
                                                int rank) {
  SolveFixture fixture;
  bool passed = expect(initialize_fixture(communicator,
                                          LinearAlgorithm::fgmres, 12U,
                                          fixture),
                       rank, "supplemental audit fixture initializes");
  if (!all_true(passed, communicator)) return false;
  TridiagonalOperator op(communicator, fixture.local, fixture.expected, -1.2,
                         3.0, -0.7, false);
  ScalingPreconditioner preconditioner(fixture.expected, 1.0 / 3.0, false);
  const LinearSolveControl selected = control(250U, 12U);

  fill_system(fixture, -1.2, 3.0, -0.7);
  const LinearSolveResult baseline = solve_fgmres(
      op, preconditioner, invocation(fixture, selected), fixture.workspace,
      fixture.reductions);

  fill_system(fixture, -1.2, 3.0, -0.7);
  RejectOnceConvergenceAudit audit;
  LinearSolveInvocation audited_invocation = invocation(fixture, selected);
  audited_invocation.convergence_audit = &audit;
  const LinearSolveResult audited =
      solve_fgmres(op, preconditioner, audited_invocation, fixture.workspace,
                   fixture.reductions);
  const ResidualOracle oracle = independent_true_residual(
      fixture, communicator, -1.2, 3.0, -0.7);
  passed &= expect(
      baseline.status.code == StatusCode::ok &&
          audited.status.code == StatusCode::ok &&
          audited.termination == LinearTermination::converged &&
          audited.iterations > baseline.iterations &&
          audited.convergence_audits == 2U &&
          audited.convergence_rejections == 1U &&
          audited.final_convergence_metric == 0.5 &&
          audited.convergence_limit == 1.0 &&
          residual_is_accepted(oracle, selected),
      rank,
      "FGMRES keeps the canonical true-residual gate and continues after one native audit rejection");
  return all_true(passed, communicator);
}

void set_uniform(FieldView field, double value) noexcept {
  for (std::int32_t cell = 0; cell < field.interior.x; ++cell) {
    field.unchecked({cell, 0, 0}, 0U) = value;
  }
  ++field.revision;
}

bool field_is_uniform(ConstFieldView field, double value,
                      double tolerance = 0.0) noexcept {
  for (std::int32_t cell = 0; cell < field.interior.x; ++cell) {
    if (std::abs(field.unchecked({cell, 0, 0}, 0U) - value) > tolerance) {
      return false;
    }
  }
  return true;
}

using DenseVector =
    std::array<double, static_cast<std::size_t>(kGlobalCells)>;

bool gather_dense_field(ConstFieldView field, const Partition& local,
                        MPI_Comm communicator, DenseVector& values) noexcept {
  if (local.size <= 0 || local.size > 4 || field.base == nullptr ||
      field.interior.x != local.cells) {
    return false;
  }
  std::array<int, 4U> counts{};
  std::array<int, 4U> displacements{};
  for (int other = 0; other < local.size; ++other) {
    const Partition other_partition = partition(other, local.size);
    counts[static_cast<std::size_t>(other)] = other_partition.cells;
    displacements[static_cast<std::size_t>(other)] = other_partition.begin;
  }
  return MPI_Allgatherv(field.base, local.cells, MPI_DOUBLE, values.data(),
                        counts.data(), displacements.data(), MPI_DOUBLE,
                        communicator) == MPI_SUCCESS;
}

DenseVector dense_difference(const DenseVector& after,
                             const DenseVector& before) noexcept {
  DenseVector result{};
  for (std::size_t index = 0U; index < result.size(); ++index) {
    result[index] = after[index] - before[index];
  }
  return result;
}

bool same_dense(const DenseVector& left, const DenseVector& right) noexcept {
  for (std::size_t index = 0U; index < left.size(); ++index) {
    if (left[index] != right[index]) return false;
  }
  return true;
}

struct CycleReconstructionProbe {
  LinearSolveResult result{};
  DenseVector internal_solution{};
  std::array<DenseVector, kLinearRecycleMaximumDirections> corrections{};
  std::array<DenseVector, kLinearRecycleMaximumDirections> raw_z{};
  std::size_t correction_count{};
  bool caller_unchanged{};
  bool variable_temporal{};
  bool variable_spatial{};
};

bool run_cycle_reconstruction_probe(
    MPI_Comm communicator, std::uint32_t restart,
    std::uint32_t maximum_iterations, bool variable_preconditioner,
    double caller_seed, CycleReconstructionProbe& output) {
  SolveFixture fixture;
  if (!initialize_fixture(communicator, LinearAlgorithm::fgmres, restart,
                          fixture)) {
    return false;
  }
  fill_system(fixture, -1.2, 3.0, -0.7, caller_seed);
  const SolutionSnapshot caller_before = snapshot_solution(fixture.solution.view);
  const Int3 shape{fixture.local.cells, 1, 1};
  if (!fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) {
    return false;
  }
  FlexibleApplicationTrace trace;
  TridiagonalOperator linear_operator(
      communicator, fixture.local, fixture.expected, -1.2, 3.0, -0.7, false,
      -1, 0U, variable_preconditioner ? &trace : nullptr);
  ScalingPreconditioner preconditioner(
      fixture.expected, 1.0 / 3.0, false, variable_preconditioner,
      fixture.local.begin, fixture.local.rank, -1, 0U,
      variable_preconditioner ? &trace : nullptr);
  LinearSolveControl selected = control(maximum_iterations, restart);
  selected.absolute_tolerance = 1.0e-300;
  selected.relative_tolerance = 0.0;
  output.result = solve_fgmres(
      linear_operator, preconditioner, invocation(fixture, selected),
      fixture.workspace, fixture.reductions);
  output.caller_unchanged = same_solution(fixture.solution.view, caller_before);
  output.variable_temporal = trace.saw_temporal_variation;
  output.variable_spatial = trace.saw_spatial_variation;
  bool gathered = gather_dense_field(
      as_const(fixture.workspace.vector(0U, shape)), fixture.local,
      communicator, output.internal_solution);

  const std::uint64_t local_count =
      fixture.workspace.recycle_correction_count_for_test();
  std::uint64_t minimum_count = 0U;
  std::uint64_t maximum_count = 0U;
  const bool count_agreed =
      MPI_Allreduce(&local_count, &minimum_count, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) == MPI_SUCCESS &&
      MPI_Allreduce(&local_count, &maximum_count, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) == MPI_SUCCESS &&
      minimum_count == maximum_count &&
      maximum_count <= kLinearRecycleMaximumDirections;
  if (!count_agreed) return false;
  output.correction_count = static_cast<std::size_t>(maximum_count);
  for (std::size_t index = 0U; index < output.correction_count; ++index) {
    gathered = gather_dense_field(
                   fixture.workspace.recycle_correction_for_test(index, shape),
                   fixture.local, communicator, output.corrections[index]) &&
               gathered;
  }
  if (variable_preconditioner) {
    const std::uint8_t z_begin = static_cast<std::uint8_t>(restart + 2U);
    for (std::size_t index = 0U; index < output.raw_z.size(); ++index) {
      gathered = gather_dense_field(
                     as_const(fixture.workspace.vector(
                         static_cast<std::uint8_t>(z_begin + index), shape)),
                     fixture.local, communicator, output.raw_z[index]) &&
                 gathered;
    }
  }
  return gathered;
}

void set_partition_values(FieldView field, const Partition& local,
                          const DenseVector& values) noexcept {
  for (std::int32_t cell = 0; cell < field.interior.x; ++cell) {
    field.unchecked({cell, 0, 0}, 0U) =
        values[static_cast<std::size_t>(local.begin + cell)];
  }
  ++field.revision;
}

struct GuardSelectionProbe {
  LinearSolveResult result{};
  LinearReductionCounters reduction_work{};
  DenseVector solution{};
  DenseVector correction{};
  std::uint32_t operator_calls{};
  std::uint32_t preconditioner_calls{};
  std::size_t correction_count{};
  std::size_t allocations{};
};

bool run_guard_selection_probe(MPI_Comm communicator, double seed_scale,
                               double absolute_tolerance, bool capture,
                               GuardSelectionProbe& output) {
  SolveFixture fixture;
  if (!initialize_fixture(communicator, LinearAlgorithm::fgmres, 4U,
                          fixture)) {
    return false;
  }
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    const double rhs = exact_solution(fixture.local.begin + cell);
    fixture.rhs.view.unchecked({cell, 0, 0}, 0U) = rhs;
    fixture.solution.view.unchecked({cell, 0, 0}, 0U) = seed_scale * rhs;
  }
  const Int3 shape{fixture.local.cells, 1, 1};
  if (capture &&
      !fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) {
    return false;
  }
  TridiagonalOperator linear_operator(
      communicator, fixture.local, fixture.expected, 0.0, 1.0, 0.0, false);
  ScalingPreconditioner preconditioner(fixture.expected, 1.0, false);
  LinearSolveControl selected = control(8U, 4U);
  selected.absolute_tolerance = absolute_tolerance;
  selected.relative_tolerance = 0.0;
  const LinearReductionCounters before = fixture.reductions.counters();
  {
    allocation_observer::Guard guard;
    output.result = solve_fgmres(
        linear_operator, preconditioner, invocation(fixture, selected),
        fixture.workspace, fixture.reductions);
    output.allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  output.reduction_work =
      reduction_delta(fixture.reductions.counters(), before);
  output.operator_calls = linear_operator.calls();
  output.preconditioner_calls = preconditioner.calls();
  if (!gather_dense_field(as_const(fixture.solution.view), fixture.local,
                          communicator, output.solution)) {
    return false;
  }
  output.correction_count =
      fixture.workspace.recycle_correction_count_for_test();
  if (output.correction_count > 0U &&
      !gather_dense_field(
          fixture.workspace.recycle_correction_for_test(0U, shape),
          fixture.local, communicator, output.correction)) {
    return false;
  }
  return output.correction_count <= kLinearRecycleMaximumDirections;
}

bool same_ordinary_result(const LinearSolveResult& left,
                          const LinearSolveResult& right) noexcept {
  return left.status.code == right.status.code &&
         left.status.detail == right.status.detail &&
         left.termination == right.termination &&
         left.iterations == right.iterations &&
         left.initial_true_residual == right.initial_true_residual &&
         left.final_true_residual == right.final_true_residual &&
         left.recursive_residual == right.recursive_residual &&
         left.reduction_calls == right.reduction_calls &&
         left.operator_applies == right.operator_applies &&
         left.preconditioner_applies == right.preconditioner_applies &&
         left.norm_breakdown_restarts == right.norm_breakdown_restarts &&
         left.convergence_audits == right.convergence_audits &&
         left.convergence_rejections == right.convergence_rejections &&
         left.final_convergence_metric == right.final_convergence_metric &&
         left.convergence_limit == right.convergence_limit &&
         left.lowest_failing_rank == right.lowest_failing_rank;
}

bool same_complete_result(const LinearSolveResult& left,
                          const LinearSolveResult& right) noexcept {
  return same_ordinary_result(left, right) &&
         left.recycle_offered_directions ==
             right.recycle_offered_directions &&
         left.recycle_retained_directions ==
             right.recycle_retained_directions &&
         left.recycle_operator_applies == right.recycle_operator_applies &&
         left.recycle_reduction_calls == right.recycle_reduction_calls &&
         left.recycle_projection_attempted ==
             right.recycle_projection_attempted &&
         left.recycle_projection_accepted ==
             right.recycle_projection_accepted &&
         left.recycle_projected_true_residual ==
             right.recycle_projected_true_residual &&
         left.recycle_cycle_corrections ==
             right.recycle_cycle_corrections &&
         left.recycle_capture_vector_passes ==
             right.recycle_capture_vector_passes &&
         left.recycle_capture_cycle_attempts ==
             right.recycle_capture_cycle_attempts &&
         left.recycle_capture_reduction_calls ==
             right.recycle_capture_reduction_calls &&
         left.recycle_capture_blocking_operations ==
             right.recycle_capture_blocking_operations;
}

bool capture_work_is_one_cycle(const GuardSelectionProbe& capture,
                               const GuardSelectionProbe& ordinary) noexcept {
  return capture.result.recycle_cycle_corrections == 1U &&
         capture.result.recycle_capture_vector_passes == 2U &&
         capture.result.recycle_capture_cycle_attempts == 1U &&
         capture.result.recycle_capture_reduction_calls == 1U &&
         capture.result.recycle_capture_blocking_operations == 2U &&
         capture.reduction_work.calls == ordinary.reduction_work.calls + 1U &&
         capture.reduction_work.scalars ==
             ordinary.reduction_work.scalars + 1U &&
         capture.reduction_work.logical_bytes ==
             ordinary.reduction_work.logical_bytes + sizeof(double) &&
         capture.reduction_work.tree_messages ==
             ordinary.reduction_work.tree_messages &&
         capture.reduction_work.blocking_operations ==
             ordinary.reduction_work.blocking_operations + 2U;
}

DenseVector dense_apply(const DenseVector& input, double lower,
                        double diagonal, double upper) noexcept {
  DenseVector output{};
  for (std::int32_t global = 0; global < kGlobalCells; ++global) {
    const double west = global == 0
                            ? 0.0
                            : input[static_cast<std::size_t>(global - 1)];
    const double east = global + 1 == kGlobalCells
                            ? 0.0
                            : input[static_cast<std::size_t>(global + 1)];
    output[static_cast<std::size_t>(global)] =
        lower * west + diagonal * input[static_cast<std::size_t>(global)] +
        upper * east;
  }
  return output;
}

double dense_dot(const DenseVector& left, const DenseVector& right) noexcept {
  double result = 0.0;
  for (std::size_t index = 0U; index < left.size(); ++index) {
    result += left[index] * right[index];
  }
  return result;
}

double dense_norm(const DenseVector& value) noexcept {
  return std::sqrt(std::max(0.0, dense_dot(value, value)));
}

struct DenseProjectionOracle {
  DenseVector candidate{};
  DenseVector residual{};
  std::array<DenseVector, kLinearRecycleMaximumDirections> images{};
  std::array<DenseVector, kLinearRecycleMaximumDirections> directions{};
  std::array<std::size_t, kLinearRecycleMaximumDirections> retained_indices{};
  std::size_t retained{};
  double initial_residual{};
  double projected_residual{};
};

DenseProjectionOracle dense_projection_oracle(
    const std::array<DenseVector, kLinearRecycleMaximumDirections>& input,
    const DenseVector& rhs, const DenseVector& initial, double lower,
    double diagonal, double upper) noexcept {
  DenseProjectionOracle result;
  result.candidate = initial;
  const DenseVector initial_operator =
      dense_apply(initial, lower, diagonal, upper);
  for (std::size_t index = 0U; index < rhs.size(); ++index) {
    result.residual[index] = rhs[index] - initial_operator[index];
  }
  result.initial_residual = dense_norm(result.residual);
  double maximum_norm = 0.0;
  for (std::size_t index = 0U; index < input.size(); ++index) {
    result.directions[index] = input[index];
    result.images[index] = dense_apply(input[index], lower, diagonal, upper);
    maximum_norm = std::max(maximum_norm, dense_norm(result.images[index]));
  }
  const double rank_threshold =
      64.0 * std::numeric_limits<double>::epsilon() * maximum_norm;
  for (std::size_t index = 0U; index < input.size(); ++index) {
    if (!(dense_norm(result.images[index]) > rank_threshold)) continue;
    for (std::size_t pass = 0U; pass < 2U; ++pass) {
      for (std::size_t prior = 0U; prior < result.retained; ++prior) {
        const std::size_t prior_index = result.retained_indices[prior];
        const double coefficient =
            dense_dot(result.images[prior_index], result.images[index]);
        for (std::size_t cell = 0U; cell < result.images[index].size();
             ++cell) {
          result.images[index][cell] -=
              coefficient * result.images[prior_index][cell];
          result.directions[index][cell] -=
              coefficient * result.directions[prior_index][cell];
        }
      }
    }
    const double image_norm = dense_norm(result.images[index]);
    if (!(image_norm > rank_threshold)) continue;
    const double inverse_norm = 1.0 / image_norm;
    for (double& value : result.images[index]) value *= inverse_norm;
    for (double& value : result.directions[index]) value *= inverse_norm;
    result.retained_indices[result.retained] = index;
    ++result.retained;
  }
  for (std::size_t index = 0U; index < result.retained; ++index) {
    const std::size_t retained_index = result.retained_indices[index];
    const double alpha =
        dense_dot(result.images[retained_index], result.residual);
    for (std::size_t cell = 0U; cell < result.candidate.size(); ++cell) {
      result.candidate[cell] +=
          alpha * result.directions[retained_index][cell];
    }
  }
  const DenseVector candidate_operator =
      dense_apply(result.candidate, lower, diagonal, upper);
  for (std::size_t index = 0U; index < result.residual.size(); ++index) {
    result.residual[index] = rhs[index] - candidate_operator[index];
  }
  result.projected_residual = dense_norm(result.residual);
  return result;
}

bool test_c1_residual_guarded_warm_start(MPI_Comm communicator, int rank) {
  bool passed = true;
  DenseVector rhs{};
  DenseVector zero_seed{};
  DenseVector better_seed{};
  DenseVector tie_seed{};
  for (std::int32_t global = 0; global < kGlobalCells; ++global) {
    const std::size_t index = static_cast<std::size_t>(global);
    rhs[index] = exact_solution(global);
    better_seed[index] = 0.5 * rhs[index];
    tie_seed[index] = 2.0 * rhs[index];
  }
  const double rhs_norm = scaled_l2(rhs);
  const double norm_tolerance =
      512.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, rhs_norm);

  GuardSelectionProbe cold;
  GuardSelectionProbe worse;
  const bool cold_ran =
      run_guard_selection_probe(communicator, 0.0, 1.0e-13, true, cold);
  const bool worse_ran =
      run_guard_selection_probe(communicator, 3.0, 1.0e-13, true, worse);
  const bool worse_matches_cold =
      cold_ran && worse_ran && cold.result.status.code == StatusCode::ok &&
      worse.result.status.code == StatusCode::ok &&
      cold.result.termination == LinearTermination::converged &&
      worse.result.termination == LinearTermination::converged &&
      same_complete_result(cold.result, worse.result) &&
      same_reduction_work(cold.reduction_work, worse.reduction_work) &&
      cold.operator_calls == worse.operator_calls &&
      cold.preconditioner_calls == worse.preconditioner_calls &&
      cold.correction_count == 1U && worse.correction_count == 1U &&
      same_dense(cold.solution, worse.solution) &&
      same_dense(cold.correction, worse.correction) &&
      std::abs(worse.result.initial_true_residual - rhs_norm) <=
          norm_tolerance &&
      cold.allocations == 0U && worse.allocations == 0U;
  if (!worse_matches_cold && rank == 0) {
    std::cerr << "guard cold/worse init/iter/op/pc/red/capture/alloc="
              << cold.result.initial_true_residual << '/'
              << worse.result.initial_true_residual << ' '
              << cold.result.iterations << '/' << worse.result.iterations
              << ' ' << cold.operator_calls << '/' << worse.operator_calls
              << ' ' << cold.preconditioner_calls << '/'
              << worse.preconditioner_calls << ' '
              << cold.reduction_work.blocking_operations << '/'
              << worse.reduction_work.blocking_operations << ' '
              << cold.correction_count << '/' << worse.correction_count
              << ' ' << cold.allocations << '/' << worse.allocations << '\n';
  }
  passed &= expect(
      worse_matches_cold, rank,
      "strictly worse C1 seed falls back to the exact zero-seed result, work, chronological correction and zero-allocation path");

  const auto retained_probe = [&](double seed_scale,
                                  const DenseVector& seed,
                                  double expected_fraction,
                                  std::string_view description) {
    GuardSelectionProbe capture;
    GuardSelectionProbe ordinary;
    const bool capture_ran = run_guard_selection_probe(
        communicator, seed_scale, 1.0e-13, true, capture);
    const bool ordinary_ran = run_guard_selection_probe(
        communicator, seed_scale, 1.0e-13, false, ordinary);
    const DenseVector expected_correction =
        dense_difference(capture.solution, seed);
    const bool retained =
        capture_ran && ordinary_ran &&
        capture.result.status.code == StatusCode::ok &&
        ordinary.result.status.code == StatusCode::ok &&
        same_ordinary_result(capture.result, ordinary.result) &&
        capture.operator_calls == ordinary.operator_calls &&
        capture.preconditioner_calls == ordinary.preconditioner_calls &&
        same_dense(capture.solution, ordinary.solution) &&
        capture.correction_count == 1U &&
        same_dense(capture.correction, expected_correction) &&
        std::abs(capture.result.initial_true_residual -
                 expected_fraction * rhs_norm) <= norm_tolerance &&
        capture_work_is_one_cycle(capture, ordinary) &&
        capture.allocations == 0U && ordinary.allocations == 0U;
    if (!retained && rank == 0) {
      std::cerr << "guard retained scale/init/iter/op/pc/work=" << seed_scale
                << ' ' << capture.result.initial_true_residual << '/'
                << ordinary.result.initial_true_residual << ' '
                << capture.result.iterations << '/'
                << ordinary.result.iterations << ' '
                << capture.operator_calls << '/' << ordinary.operator_calls
                << ' ' << capture.preconditioner_calls << '/'
                << ordinary.preconditioner_calls << ' '
                << capture.reduction_work.calls << '/'
                << ordinary.reduction_work.calls << ' '
                << capture.reduction_work.blocking_operations << '/'
                << ordinary.reduction_work.blocking_operations << '\n';
    }
    return expect(retained, rank, description);
  };
  passed &= retained_probe(
      0.5, better_seed, 0.5,
      "strictly better C1 seed is retained with its independent true residual and only the frozen one-cycle capture work");
  passed &= retained_probe(
      2.0, tie_seed, 1.0,
      "exact residual-norm tie retains the seed and its signed accepted correction rather than resetting to zero");

  GuardSelectionProbe cold_early;
  GuardSelectionProbe worse_early;
  const bool cold_early_ran = run_guard_selection_probe(
      communicator, 0.0, 1.0e6, true, cold_early);
  const bool worse_early_ran = run_guard_selection_probe(
      communicator, 3.0, 1.0e6, true, worse_early);
  const bool early_exception =
      cold_early_ran && worse_early_ran &&
      cold_early.result.status.code == StatusCode::ok &&
      worse_early.result.status.code == StatusCode::ok &&
      same_complete_result(cold_early.result, worse_early.result) &&
      same_dense(cold_early.solution, zero_seed) &&
      same_dense(worse_early.solution, zero_seed) &&
      cold_early.operator_calls == 1U && worse_early.operator_calls == 1U &&
      cold_early.preconditioner_calls == 0U &&
      worse_early.preconditioner_calls == 0U &&
      cold_early.correction_count == 0U &&
      worse_early.correction_count == 0U &&
      cold_early.reduction_work.calls == worse_early.reduction_work.calls &&
      cold_early.reduction_work.scalars ==
          worse_early.reduction_work.scalars &&
      cold_early.reduction_work.logical_bytes ==
          worse_early.reduction_work.logical_bytes &&
      cold_early.reduction_work.tree_messages ==
          worse_early.reduction_work.tree_messages &&
      worse_early.reduction_work.blocking_operations ==
          cold_early.reduction_work.blocking_operations + 1U &&
      cold_early.allocations == 0U && worse_early.allocations == 0U;
  if (!early_exception && rank == 0) {
    std::cerr << "guard early init/iter/op/pc/work/alloc="
              << cold_early.result.initial_true_residual << '/'
              << worse_early.result.initial_true_residual << ' '
              << cold_early.result.iterations << '/'
              << worse_early.result.iterations << ' '
              << cold_early.operator_calls << '/'
              << worse_early.operator_calls << ' '
              << cold_early.preconditioner_calls << '/'
              << worse_early.preconditioner_calls << ' '
              << cold_early.reduction_work.blocking_operations << '/'
              << worse_early.reduction_work.blocking_operations << ' '
              << cold_early.allocations << '/' << worse_early.allocations
              << '\n';
  }
  passed &= expect(
      early_exception, rank,
      "worse-seed selected-zero early convergence adds exactly one consensus and no operator, preconditioner, checked-reduction or allocation work");

  SolveFixture prerequisite_fixture;
  const bool prerequisite_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, prerequisite_fixture);
  passed &= expect(prerequisite_initialized, rank,
                   "guard prerequisite failure fixture initializes");
  if (!all_true(prerequisite_initialized, communicator)) return false;
  fill_system(prerequisite_fixture, 0.0, 1.0, 0.0, 3.0);
  const SolutionSnapshot prerequisite_caller =
      snapshot_solution(prerequisite_fixture.solution.view);
  const Int3 prerequisite_shape{prerequisite_fixture.local.cells, 1, 1};
  const int failing_rank = prerequisite_fixture.local.size - 1;
  const bool prerequisite_capture = static_cast<bool>(
      prerequisite_fixture.workspace.recycle_begin_capture_for_test(
          prerequisite_shape, prerequisite_fixture.expected.fingerprint));
  const Status local_prerequisite =
      prerequisite_fixture.local.rank == failing_rank
          ? Status{StatusCode::numerical_failure, 890U}
          : Status{};
  const Status prerequisite_failure =
      prerequisite_fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(prerequisite_fixture.solution.view),
          prerequisite_fixture.reductions, local_prerequisite);
  const bool prerequisite_failed_cleanly =
      prerequisite_capture &&
      prerequisite_failure.code == StatusCode::numerical_failure &&
      prerequisite_failure.detail == 890U &&
      prerequisite_fixture.reductions.lowest_failing_rank() == failing_rank &&
      same_solution(prerequisite_fixture.solution.view,
                    prerequisite_caller) &&
      prerequisite_fixture.workspace.recycle_correction_count_for_test() ==
          0U;
  const bool prerequisite_retry =
      static_cast<bool>(
          prerequisite_fixture.workspace.recycle_begin_capture_for_test(
              prerequisite_shape,
              prerequisite_fixture.expected.fingerprint)) &&
      static_cast<bool>(
          prerequisite_fixture.workspace.recycle_capture_cycle_start_for_test(
              as_const(prerequisite_fixture.solution.view),
              prerequisite_fixture.reductions));
  set_uniform(prerequisite_fixture.solution.view, 4.0);
  const bool prerequisite_retry_published = static_cast<bool>(
      prerequisite_fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(prerequisite_fixture.solution.view),
          prerequisite_fixture.reductions));
  passed &= expect(
      prerequisite_failed_cleanly && prerequisite_retry &&
          prerequisite_retry_published &&
          prerequisite_fixture.workspace.recycle_correction_count_for_test() ==
              1U,
      rank,
      "rank-local fallback revision prerequisite preserves the true lowest rank, clears capture once and permits a clean publish retry");
  prerequisite_fixture.workspace.recycle_clear_for_test();
  return all_true(passed, communicator);
}

bool test_recycle_capture_projection(MPI_Comm communicator, int rank) {
  SolveFixture fixture;
  bool passed = expect(initialize_fixture(communicator, LinearAlgorithm::fgmres,
                                           4U, fixture),
                       rank, "cycle-correction fixture initializes");
  if (!all_true(passed, communicator)) return false;
  const Int3 shape{fixture.local.cells, 1, 1};
  passed &= expect(
      static_cast<bool>(fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)),
      rank, "cycle-correction capture begins with the source identity");
  if (!all_true(passed, communicator)) return false;

  double base = 0.0;
  const std::uint8_t initial_snapshot_slot =
      fixture.workspace.recycle_snapshot_slot_for_test();
  std::array<std::uint8_t, 6U> snapshot_after_publish{};
  for (std::size_t cycle = 0U; cycle < 6U; ++cycle) {
    set_uniform(fixture.solution.view, base);
    const Status started =
        fixture.workspace.recycle_capture_cycle_start_for_test(
            as_const(fixture.solution.view), fixture.reductions);
    passed &= expect(static_cast<bool>(started), rank,
                     "cycle snapshot captures the accepted-cycle start");
    const double delta = static_cast<double>(cycle + 1U);
    base += delta;
    set_uniform(fixture.solution.view, base);
    const Status published =
        fixture.workspace.recycle_capture_cycle_publish_for_test(
            as_const(fixture.solution.view), fixture.reductions);
    passed &= expect(static_cast<bool>(published), rank,
                     "cycle correction publishes after the accepted cycle");
    snapshot_after_publish[cycle] =
        fixture.workspace.recycle_snapshot_slot_for_test();
  }
  // A completed cycle whose accepted solution is unchanged still performs
  // both local vector passes and one global nonzero check, but must not enter
  // the offered ring.
  set_uniform(fixture.solution.view, base);
  const Status zero_started =
      fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(fixture.solution.view), fixture.reductions);
  const Status zero_published =
      fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(fixture.solution.view), fixture.reductions);
  passed &= expect(static_cast<bool>(zero_started) &&
                       static_cast<bool>(zero_published),
                   rank, "zero correction completes without offering a direction");
  const std::array<std::uint8_t, 6U> expected_snapshots{
      static_cast<std::uint8_t>(initial_snapshot_slot + 1U),
      static_cast<std::uint8_t>(initial_snapshot_slot + 2U),
      static_cast<std::uint8_t>(initial_snapshot_slot + 3U),
      static_cast<std::uint8_t>(initial_snapshot_slot + 4U),
      initial_snapshot_slot,
      static_cast<std::uint8_t>(initial_snapshot_slot + 1U)};
  bool physical_order =
      snapshot_after_publish == expected_snapshots &&
      fixture.workspace.recycle_snapshot_slot_for_test() ==
          expected_snapshots.back();
  const std::array<std::uint8_t, 4U> expected_physical_order{
      static_cast<std::uint8_t>(initial_snapshot_slot + 2U),
      static_cast<std::uint8_t>(initial_snapshot_slot + 3U),
      static_cast<std::uint8_t>(initial_snapshot_slot + 4U),
      initial_snapshot_slot};
  for (std::size_t index = 0U; index < expected_physical_order.size();
       ++index) {
    physical_order =
        fixture.workspace.recycle_correction_physical_slot_for_test(index) ==
            expected_physical_order[index] &&
        physical_order;
  }
  passed &= expect(
      fixture.workspace.recycle_correction_count_for_test() == 4U &&
          fixture.workspace.recycle_capture_cycle_corrections_for_test() ==
              6U &&
          fixture.workspace.recycle_capture_cycle_attempts_for_test() == 7U &&
          fixture.workspace.recycle_capture_vector_passes_for_test() == 14U &&
          fixture.workspace.recycle_capture_reduction_calls_for_test() == 7U &&
          fixture.workspace.recycle_capture_blocking_operations_for_test() ==
              14U && physical_order,
      rank,
      "capture reports zero-cycle work separately from the cap-four ring");
  const std::array<double, 4U> expected_directions{3.0, 4.0, 5.0, 6.0};
  for (std::size_t index = 0U; index < expected_directions.size(); ++index) {
    const ConstFieldView correction =
        fixture.workspace.recycle_correction_for_test(index, shape);
    if (correction.base != nullptr &&
        !field_is_uniform(correction, expected_directions[index])) {
      std::cerr << "rank " << rank << " ring index " << index << " value="
                << correction.unchecked({0, 0, 0}, 0U) << " expected="
                << expected_directions[index] << '\n';
    }
    passed &= expect(
        correction.base != nullptr &&
            field_is_uniform(correction, expected_directions[index]),
        rank, "retained corrections follow chronological last-four order");
  }
  if (!all_true(passed, communicator)) return false;

  fixture.workspace.recycle_clear_for_test();
  set_uniform(fixture.solution.view, 0.0);
  const bool identity_capture =
      static_cast<bool>(fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) &&
      static_cast<bool>(fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  set_uniform(fixture.solution.view, 1.0);
  const bool identity_publish = static_cast<bool>(
      fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  set_uniform(fixture.solution.view, 0.0);
  TridiagonalOperator identity_operator(communicator, fixture.local,
                                        fixture.expected, 0.0, 1.0, 0.0,
                                        false);
  ScalingPreconditioner identity_preconditioner(fixture.expected, 1.0, false);
  const LinearSolveControl selected = control(50U, 4U);
  set_uniform(fixture.rhs.view, 5.0);
  passed &= expect(identity_capture && identity_publish, rank,
                   "identity projection captures one finite correction");
  passed &= expect(static_cast<bool>(
                       fixture.workspace.recycle_begin_projection_for_test(
                           shape, fixture.expected.fingerprint)),
                   rank, "projection begins on the refreshed current identity");
  if (!all_true(passed, communicator)) return false;
  const std::size_t correction_calls_before = identity_operator.calls();
  const LinearReductionCounters projection_counters_before =
      fixture.reductions.counters();
  LinearSolveResult projected;
  std::size_t projected_allocations = 0U;
  {
    allocation_observer::Guard guard;
    projected = solve_fgmres(identity_operator, identity_preconditioner,
                             invocation(fixture, selected), fixture.workspace,
                             fixture.reductions);
    projected_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const std::size_t correction_calls_after = identity_operator.calls();
  const LinearReductionCounters projection_counters_after =
      fixture.reductions.counters();
  passed &= expect(
      projected.status.code == StatusCode::ok &&
          projected.termination == LinearTermination::converged &&
          projected.iterations == 0U &&
          projected.recycle_offered_directions == 1U &&
          projected.recycle_retained_directions == 1U &&
          projected.recycle_projection_attempted &&
          projected.recycle_projection_accepted &&
          projected.recycle_operator_applies == 2U &&
          projected.recycle_reduction_calls == 5U &&
          projected.recycle_projected_true_residual <= 1.0e-12 &&
          projected.operator_applies == 1U &&
          correction_calls_after - correction_calls_before == 3U &&
          projection_counters_after.calls - projection_counters_before.calls ==
              9U &&
          projection_counters_after.scalars -
                  projection_counters_before.scalars ==
              9U &&
          projection_counters_after.logical_bytes -
                  projection_counters_before.logical_bytes ==
              9U * sizeof(double) &&
          projection_counters_after.blocking_operations -
                  projection_counters_before.blocking_operations ==
              19U &&
          projection_counters_after.tree_messages -
                  projection_counters_before.tree_messages ==
              0U &&
          projected_allocations == 0U &&
          field_is_uniform(as_const(fixture.solution.view), 5.0, 1.0e-12),
      rank,
      "current-operator QR projection admits the dense identity oracle with ordinary work separated");
  fixture.workspace.recycle_clear_for_test();

  // Independent dense oracle: use four globally distinct solution-space
  // corrections and a nonnormal tridiagonal current operator.  The oracle
  // below forms A*d, performs two-pass CGS, computes the least-squares
  // coefficients, and reconstructs both the candidate and its exact defect
  // without reading production projection storage or counters.
  SolveFixture dense_fixture;
  const bool dense_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, dense_fixture);
  passed &= expect(dense_initialized, rank,
                   "dense cycle-correction oracle fixture initializes");
  if (!all_true(passed, communicator)) return false;
  constexpr double dense_lower = -1.2;
  constexpr double dense_diagonal = 3.0;
  constexpr double dense_upper = -0.7;
  const Int3 dense_shape{dense_fixture.local.cells, 1, 1};
  const std::array<std::int32_t, kLinearRecycleMaximumDirections> impulses{
      7, 22, 39, 58};
  std::array<DenseVector, kLinearRecycleMaximumDirections> dense_corrections{};
  DenseVector accepted_solution{};
  passed &= expect(static_cast<bool>(
                       dense_fixture.workspace.recycle_begin_capture_for_test(
                           dense_shape, dense_fixture.expected.fingerprint)),
                   rank, "dense oracle capture begins");
  for (std::size_t cycle = 0U;
       cycle < kLinearRecycleMaximumDirections; ++cycle) {
    DenseVector direction{};
    if (cycle < 3U) {
      direction[static_cast<std::size_t>(impulses[cycle])] =
          1.0 + static_cast<double>(cycle) * 0.25;
    } else {
      // A deliberately dependent fourth correction exercises deterministic
      // post-CGS deflation while retaining a nontrivial A2 image.
      for (std::size_t cell = 0U; cell < direction.size(); ++cell) {
        direction[cell] = 0.4 * dense_corrections[0U][cell] +
                          0.8 * dense_corrections[1U][cell];
      }
    }
    dense_corrections[cycle] = direction;
    set_partition_values(dense_fixture.solution.view, dense_fixture.local,
                         accepted_solution);
    const Status started =
        dense_fixture.workspace.recycle_capture_cycle_start_for_test(
            as_const(dense_fixture.solution.view), dense_fixture.reductions);
    for (std::size_t cell = 0U; cell < accepted_solution.size(); ++cell) {
      accepted_solution[cell] += direction[cell];
    }
    set_partition_values(dense_fixture.solution.view, dense_fixture.local,
                         accepted_solution);
    const Status published =
        dense_fixture.workspace.recycle_capture_cycle_publish_for_test(
            as_const(dense_fixture.solution.view), dense_fixture.reductions);
    passed &= expect(static_cast<bool>(started) && static_cast<bool>(published),
                     rank, "dense oracle publishes independent correction");
  }
  DenseVector target{};
  const std::array<double, kLinearRecycleMaximumDirections> target_coefficients{
      0.75, -1.25, 0.5, 1.1};
  for (std::size_t direction = 0U;
       direction < kLinearRecycleMaximumDirections; ++direction) {
    for (std::size_t cell = 0U; cell < target.size(); ++cell) {
      target[cell] += target_coefficients[direction] *
                      dense_corrections[direction][cell];
    }
  }
  const DenseVector dense_rhs =
      dense_apply(target, dense_lower, dense_diagonal, dense_upper);
  set_partition_values(dense_fixture.rhs.view, dense_fixture.local, dense_rhs);
  DenseVector dense_zero{};
  set_partition_values(dense_fixture.solution.view, dense_fixture.local,
                       dense_zero);
  passed &= expect(static_cast<bool>(
                       dense_fixture.workspace.recycle_begin_projection_for_test(
                           dense_shape, dense_fixture.expected.fingerprint)),
                   rank, "dense oracle projection begins on current operator");
  TridiagonalOperator dense_operator(
      communicator, dense_fixture.local, dense_fixture.expected, dense_lower,
      dense_diagonal, dense_upper, false);
  ScalingPreconditioner dense_preconditioner(dense_fixture.expected, 1.0 / 3.0,
                                             false);
  LinearSolveResult dense_result;
  std::size_t dense_allocations = 0U;
  {
    allocation_observer::Guard guard;
    dense_result = solve_fgmres(
        dense_operator, dense_preconditioner,
        invocation(dense_fixture, selected), dense_fixture.workspace,
        dense_fixture.reductions);
    dense_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const DenseProjectionOracle dense_oracle = dense_projection_oracle(
      dense_corrections, dense_rhs, dense_zero, dense_lower, dense_diagonal,
      dense_upper);
  std::array<int, 4U> dense_counts{};
  std::array<int, 4U> dense_displacements{};
  for (int other = 0; other < dense_fixture.local.size; ++other) {
    const Partition other_partition =
        partition(other, dense_fixture.local.size);
    dense_counts[static_cast<std::size_t>(other)] = other_partition.cells;
    dense_displacements[static_cast<std::size_t>(other)] =
        other_partition.begin;
  }
  DenseVector dense_solution{};
  const int dense_gather_status = MPI_Allgatherv(
      dense_fixture.solution.view.base, dense_fixture.local.cells, MPI_DOUBLE,
      dense_solution.data(), dense_counts.data(), dense_displacements.data(),
      MPI_DOUBLE, communicator);
  bool dense_candidate_matches = dense_gather_status == MPI_SUCCESS;
  bool dense_orthogonal = true;
  bool dense_a2u_relation = true;
  for (std::size_t left = 0U; left < dense_oracle.retained; ++left) {
    const std::size_t left_index = dense_oracle.retained_indices[left];
    const DenseVector transformed_direction = dense_apply(
        dense_oracle.directions[left_index], dense_lower, dense_diagonal,
        dense_upper);
    for (std::size_t cell = 0U; cell < transformed_direction.size(); ++cell) {
      dense_a2u_relation =
          std::abs(transformed_direction[cell] -
                   dense_oracle.images[left_index][cell]) <= 1.0e-9 &&
          dense_a2u_relation;
    }
    dense_orthogonal =
        std::abs(dense_dot(dense_oracle.images[left_index],
                           dense_oracle.images[left_index]) -
                 1.0) <= 1.0e-9 &&
        dense_orthogonal;
    dense_orthogonal =
        std::abs(dense_dot(dense_oracle.images[left_index],
                           dense_oracle.residual)) <= 1.0e-9 &&
        dense_orthogonal;
    for (std::size_t right = left + 1U; right < dense_oracle.retained;
         ++right) {
      const std::size_t right_index = dense_oracle.retained_indices[right];
      dense_orthogonal =
          std::abs(dense_dot(dense_oracle.images[left_index],
                             dense_oracle.images[right_index])) <= 1.0e-9 &&
          dense_orthogonal;
    }
  }
  for (std::size_t cell = 0U; cell < dense_solution.size(); ++cell) {
    dense_candidate_matches =
        std::abs(dense_solution[cell] - dense_oracle.candidate[cell]) <=
            1.0e-9 &&
        dense_candidate_matches;
  }
  passed &= expect(
      dense_result.status.code == StatusCode::ok &&
          dense_result.termination == LinearTermination::converged &&
          dense_result.iterations == 0U &&
          dense_result.recycle_offered_directions == 4U &&
          dense_result.recycle_retained_directions == dense_oracle.retained &&
          dense_oracle.retained == 3U && dense_result.recycle_projection_attempted &&
          dense_result.recycle_projection_accepted &&
          dense_result.recycle_operator_applies == 5U &&
          dense_result.recycle_reduction_calls == 23U &&
          dense_result.operator_applies == 1U && dense_candidate_matches &&
          dense_orthogonal && dense_a2u_relation && dense_allocations == 0U &&
          std::abs(dense_result.initial_true_residual -
                   dense_oracle.initial_residual) <= 1.0e-9 &&
          std::abs(dense_result.recycle_projected_true_residual -
                   dense_oracle.projected_residual) <= 1.0e-9 &&
          dense_result.recycle_projected_true_residual <
              dense_result.initial_true_residual,
      rank,
      "independent dense A2U/CGS2 oracle verifies orthogonality and minimizer admission");
  dense_fixture.workspace.recycle_clear_for_test();

  // The smallest legal FGMRES restart still has to support the full
  // last-four recycle ring.  At maximum_restart=1 the fourth naive image
  // scratch slot is also the first recycle-pool slot, so this exact identity
  // oracle detects any direction/image storage alias before ordinary FGMRES
  // can hide it by converging from a damaged projection.
  {
    SolveFixture minimum_restart_fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 1U,
        minimum_restart_fixture);
    passed &= expect(initialized, rank,
                     "restart-one four-direction fixture initializes");
    if (!all_true(passed, communicator)) return false;
    const Int3 minimum_shape{minimum_restart_fixture.local.cells, 1, 1};
    const std::array<std::int32_t, kLinearRecycleMaximumDirections>
        minimum_impulses{7, 22, 39, 58};
    const std::array<double, kLinearRecycleMaximumDirections>
        minimum_coefficients{1.0, 2.0, 3.0, 4.0};
    DenseVector accepted{};
    DenseVector target{};
    bool captured = static_cast<bool>(
        minimum_restart_fixture.workspace.recycle_begin_capture_for_test(
            minimum_shape,
            minimum_restart_fixture.expected.fingerprint));
    for (std::size_t direction = 0U;
         direction < kLinearRecycleMaximumDirections && captured;
         ++direction) {
      set_partition_values(minimum_restart_fixture.solution.view,
                           minimum_restart_fixture.local, accepted);
      captured = static_cast<bool>(
          minimum_restart_fixture.workspace
              .recycle_capture_cycle_start_for_test(
                  as_const(minimum_restart_fixture.solution.view),
                  minimum_restart_fixture.reductions));
      const std::size_t impulse =
          static_cast<std::size_t>(minimum_impulses[direction]);
      accepted[impulse] += 1.0;
      target[impulse] = minimum_coefficients[direction];
      set_partition_values(minimum_restart_fixture.solution.view,
                           minimum_restart_fixture.local, accepted);
      captured = static_cast<bool>(
          minimum_restart_fixture.workspace
              .recycle_capture_cycle_publish_for_test(
                  as_const(minimum_restart_fixture.solution.view),
                  minimum_restart_fixture.reductions)) &&
                 captured;
    }
    DenseVector minimum_zero{};
    set_partition_values(minimum_restart_fixture.solution.view,
                         minimum_restart_fixture.local, minimum_zero);
    set_partition_values(minimum_restart_fixture.rhs.view,
                         minimum_restart_fixture.local, target);
    const bool projection_begun = static_cast<bool>(
        minimum_restart_fixture.workspace.recycle_begin_projection_for_test(
            minimum_shape,
            minimum_restart_fixture.expected.fingerprint));
    TridiagonalOperator minimum_operator(
        communicator, minimum_restart_fixture.local,
        minimum_restart_fixture.expected, 0.0, 1.0, 0.0, false);
    ScalingPreconditioner minimum_preconditioner(
        minimum_restart_fixture.expected, 1.0, false);
    LinearSolveResult minimum_result;
    std::size_t minimum_allocations = 0U;
    {
      allocation_observer::Guard guard;
      minimum_result = solve_fgmres(
          minimum_operator, minimum_preconditioner,
          invocation(minimum_restart_fixture, control(12U, 1U)),
          minimum_restart_fixture.workspace,
          minimum_restart_fixture.reductions);
      minimum_allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    DenseVector minimum_solution{};
    const bool gathered = gather_dense_field(
        as_const(minimum_restart_fixture.solution.view),
        minimum_restart_fixture.local, communicator, minimum_solution);
    const bool minimum_matrix =
        captured && projection_begun && gathered &&
            minimum_result.status.code == StatusCode::ok &&
            minimum_result.termination == LinearTermination::converged &&
            minimum_result.iterations == 0U &&
            minimum_result.recycle_offered_directions == 4U &&
            minimum_result.recycle_retained_directions == 4U &&
            minimum_result.recycle_operator_applies == 5U &&
            minimum_result.recycle_reduction_calls == 22U &&
            minimum_result.recycle_projection_attempted &&
            minimum_result.recycle_projection_accepted &&
            minimum_result.recycle_projected_true_residual == 0.0 &&
            same_dense(minimum_solution, target) &&
            minimum_allocations == 0U;
    if (!minimum_matrix) {
      std::size_t first_mismatch = minimum_solution.size();
      for (std::size_t index = 0U; index < minimum_solution.size(); ++index) {
        if (minimum_solution[index] != target[index]) {
          first_mismatch = index;
          break;
        }
      }
      std::cerr << "rank " << rank << " restart-one status/detail/term="
                << static_cast<unsigned>(minimum_result.status.code) << '/'
                << minimum_result.status.detail << '/'
                << static_cast<unsigned>(minimum_result.termination)
                << " iter/offered/retained/op/red="
                << minimum_result.iterations << '/'
                << minimum_result.recycle_offered_directions << '/'
                << minimum_result.recycle_retained_directions << '/'
                << minimum_result.recycle_operator_applies << '/'
                << minimum_result.recycle_reduction_calls
                << " attempted/accepted/projected="
                << minimum_result.recycle_projection_attempted << '/'
                << minimum_result.recycle_projection_accepted << '/'
                << minimum_result.recycle_projected_true_residual
                << " alloc/mismatch=" << minimum_allocations << '/'
                << first_mismatch << '\n';
    }
    passed &= expect(
        minimum_matrix,
        rank,
        "restart-one full recycle ring keeps image scratch disjoint and reaches the exact projection without ordinary iterations");
  }

  auto cycle_boundary_probe = [&](std::uint32_t maximum_iterations,
                                  std::uint32_t expected_cycles,
                                  std::string_view description,
                                  double absolute_tolerance,
                                  bool expect_success) noexcept {
    SolveFixture cycle_fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 30U, cycle_fixture);
    bool probe_passed = expect(initialized, rank,
                               "restart-30 cycle-boundary fixture initializes");
    if (!all_true(probe_passed, communicator)) return false;
    fill_system(cycle_fixture, -1.2, 3.0, -0.7);
    TridiagonalOperator cycle_operator(
        communicator, cycle_fixture.local, cycle_fixture.expected, -1.2, 3.0,
        -0.7, false);
    ScalingPreconditioner cycle_preconditioner(cycle_fixture.expected, 1.0 / 3.0,
                                               false);
    LinearSolveControl boundary_control =
        control(maximum_iterations, 30U);
    boundary_control.absolute_tolerance = absolute_tolerance;
    boundary_control.relative_tolerance = 0.0;
    const Int3 cycle_shape{cycle_fixture.local.cells, 1, 1};
    probe_passed &= expect(static_cast<bool>(
                               cycle_fixture.workspace
                                   .recycle_begin_capture_for_test(
                                       cycle_shape,
                                       cycle_fixture.expected.fingerprint)),
                           rank, "restart-30 capture begins before cycle probe");
    const LinearSolveResult cycle_result = solve_fgmres(
        cycle_operator, cycle_preconditioner,
        invocation(cycle_fixture, boundary_control), cycle_fixture.workspace,
        cycle_fixture.reductions);
    const bool terminal =
        expect_success
            ? cycle_result.status.code == StatusCode::ok &&
                  cycle_result.termination == LinearTermination::converged
            : cycle_result.status.code == StatusCode::rejected_step &&
                  cycle_result.termination ==
                      LinearTermination::maximum_iterations;
    probe_passed &= expect(
        terminal &&
            cycle_result.iterations == maximum_iterations &&
            cycle_result.recycle_cycle_corrections == expected_cycles &&
            cycle_result.recycle_capture_cycle_attempts == expected_cycles &&
            cycle_result.recycle_capture_vector_passes ==
                2U * expected_cycles &&
            cycle_result.recycle_capture_reduction_calls == expected_cycles &&
            cycle_result.recycle_capture_blocking_operations ==
                2U * expected_cycles &&
            cycle_fixture.workspace.recycle_correction_count_for_test() ==
                expected_cycles,
        rank, description);
    // A direct workspace probe deliberately leaves the capture session to its
    // owning PISO epoch.  Production failure cleanup is asserted separately;
    // clearing here before the assertion would not prove that lifecycle.
    return all_true(probe_passed, communicator);
  };
  passed &= cycle_boundary_probe(
      65U, 3U, "restart-30 bounded 30/30/5 publishes capture",
      1.0e-300, false);
  passed &= cycle_boundary_probe(
      54U, 2U, "restart-30 bounded 30/24 publishes capture",
      1.0e-300, false);
  passed &= cycle_boundary_probe(
      30U, 1U, "restart-30 successful terminal cycle publishes capture",
      5.0e-12, true);

  CycleReconstructionProbe cycle_30;
  CycleReconstructionProbe cycle_60;
  CycleReconstructionProbe cycle_65;
  CycleReconstructionProbe cycle_30_for_54;
  CycleReconstructionProbe cycle_54;
  const bool reconstructed =
      run_cycle_reconstruction_probe(communicator, 30U, 30U, false, 0.0,
                                     cycle_30) &&
      run_cycle_reconstruction_probe(communicator, 30U, 60U, false, 0.0,
                                     cycle_60) &&
      run_cycle_reconstruction_probe(communicator, 30U, 65U, false, 0.0,
                                     cycle_65) &&
      run_cycle_reconstruction_probe(communicator, 30U, 30U, false, 0.0,
                                     cycle_30_for_54) &&
      run_cycle_reconstruction_probe(communicator, 30U, 54U, false, 0.0,
                                     cycle_54);
  DenseVector zero_solution{};
  const std::array<DenseVector, 3U> expected_65{
      dense_difference(cycle_30.internal_solution, zero_solution),
      dense_difference(cycle_60.internal_solution,
                       cycle_30.internal_solution),
      dense_difference(cycle_65.internal_solution,
                       cycle_60.internal_solution)};
  const std::array<DenseVector, 2U> expected_54{
      dense_difference(cycle_30_for_54.internal_solution, zero_solution),
      dense_difference(cycle_54.internal_solution,
                       cycle_30_for_54.internal_solution)};
  bool chronological_65 = cycle_65.correction_count == expected_65.size();
  for (std::size_t index = 0U; index < expected_65.size(); ++index) {
    chronological_65 =
        same_dense(cycle_65.corrections[index], expected_65[index]) &&
        chronological_65;
  }
  bool chronological_54 = cycle_54.correction_count == expected_54.size();
  for (std::size_t index = 0U; index < expected_54.size(); ++index) {
    chronological_54 =
        same_dense(cycle_54.corrections[index], expected_54[index]) &&
        chronological_54;
  }
  const auto failed_at = [](const CycleReconstructionProbe& probe,
                            std::uint32_t iterations,
                            std::size_t corrections) noexcept {
    return probe.result.status.code == StatusCode::rejected_step &&
           probe.result.termination == LinearTermination::maximum_iterations &&
           probe.result.iterations == iterations &&
           probe.correction_count == corrections && probe.caller_unchanged;
  };
  passed &= expect(
      reconstructed && failed_at(cycle_30, 30U, 1U) &&
          failed_at(cycle_60, 60U, 2U) &&
          failed_at(cycle_65, 65U, 3U) && chronological_65,
      rank,
      "independent x30/x60/x65 reconstruction proves chronological 30/30/5 corrections and unchanged failed-solve caller bytes");
  passed &= expect(
      reconstructed && failed_at(cycle_30_for_54, 30U, 1U) &&
          failed_at(cycle_54, 54U, 2U) && chronological_54,
      rank,
      "independent x30/x54 reconstruction proves chronological 30/24 corrections and unchanged failed-solve caller bytes");

  {
    SolveFixture variable_fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, variable_fixture);
    passed &= expect(initialized, rank,
                     "variable-Z ancestry fixture initializes");
    if (!all_true(passed, communicator)) return false;
    fill_system(variable_fixture, -1.2, 3.0, -0.7);
    FlexibleApplicationTrace variable_trace;
    TridiagonalOperator variable_operator(
        communicator, variable_fixture.local, variable_fixture.expected, -1.2,
        3.0, -0.7, false, -1, 0U, &variable_trace);
    ScalingPreconditioner variable_preconditioner(
        variable_fixture.expected, 1.0 / 3.0, false, true,
        variable_fixture.local.begin, variable_fixture.local.rank, -1, 0U,
        &variable_trace);
    LinearSolveControl variable_control = control(4U, 4U);
    variable_control.absolute_tolerance = 1.0e-300;
    variable_control.relative_tolerance = 0.0;
    const Int3 variable_shape{variable_fixture.local.cells, 1, 1};
    passed &= expect(static_cast<bool>(
                         variable_fixture.workspace
                             .recycle_begin_capture_for_test(
                                 variable_shape,
                                 variable_fixture.expected.fingerprint)),
                     rank, "variable-Z capture begins from solution state");
    const LinearSolveResult variable_result = solve_fgmres(
        variable_operator, variable_preconditioner,
        invocation(variable_fixture, variable_control),
        variable_fixture.workspace, variable_fixture.reductions);
    const ConstFieldView variable_correction =
        variable_fixture.workspace.recycle_correction_for_test(0U,
                                                                 variable_shape);
    const FieldView variable_workspace_solution =
        variable_fixture.workspace.vector(0U, variable_shape);
    bool solution_space_correction = variable_correction.base != nullptr &&
                                     variable_workspace_solution.base != nullptr;
    for (std::int32_t cell = 0; cell < variable_fixture.local.cells; ++cell) {
      solution_space_correction =
          variable_correction.unchecked({cell, 0, 0}, 0U) ==
              variable_workspace_solution.unchecked({cell, 0, 0}, 0U) &&
          solution_space_correction;
    }
    passed &= expect(
        variable_result.status.code == StatusCode::rejected_step &&
            variable_result.iterations == 4U &&
            variable_result.recycle_cycle_corrections == 1U &&
            variable_result.recycle_capture_cycle_attempts == 1U &&
            variable_result.recycle_capture_vector_passes == 2U &&
            variable_result.recycle_capture_reduction_calls == 1U &&
            solution_space_correction && variable_trace.saw_temporal_variation &&
            variable_trace.saw_spatial_variation,
        rank,
        "cycle correction follows accepted x delta rather than variable-Z ancestry");
    // Test-only cleanup after the mathematical assertion; no lifecycle claim
    // is derived from this explicit clear.
    variable_fixture.workspace.recycle_clear_for_test();
  }

  CycleReconstructionProbe variable_probe;
  const bool variable_reconstructed = run_cycle_reconstruction_probe(
      communicator, 4U, 4U, true, 0.25, variable_probe);
  DenseVector variable_seed{};
  variable_seed.fill(0.25);
  const DenseVector variable_delta = dense_difference(
      variable_probe.internal_solution, variable_seed);
  bool not_any_raw_z = true;
  for (const DenseVector& raw_z : variable_probe.raw_z) {
    not_any_raw_z = !same_dense(variable_delta, raw_z) && not_any_raw_z;
  }
  passed &= expect(
      variable_reconstructed && variable_probe.caller_unchanged &&
          variable_probe.result.status.code == StatusCode::rejected_step &&
          variable_probe.result.iterations == 4U &&
          variable_probe.correction_count == 1U &&
          same_dense(variable_probe.corrections[0U], variable_delta) &&
          not_any_raw_z && variable_probe.variable_temporal &&
          variable_probe.variable_spatial,
      rank,
      "independent variable-preconditioner x4-x0 reconstruction matches the captured correction and differs from every raw Z vector");
  if (!all_true(passed, communicator)) return false;

  set_uniform(fixture.solution.view, 0.0);
  passed &= expect(static_cast<bool>(
                       fixture.workspace.recycle_begin_capture_for_test(
                           shape, fixture.expected.fingerprint)) &&
                       static_cast<bool>(
                           fixture.workspace.recycle_capture_cycle_start_for_test(
                               as_const(fixture.solution.view),
                               fixture.reductions)),
                   rank, "single correction capture is repeatable");
  set_uniform(fixture.solution.view, 1.0);
  passed &= expect(static_cast<bool>(
                       fixture.workspace.recycle_capture_cycle_publish_for_test(
                           as_const(fixture.solution.view), fixture.reductions)),
                   rank, "single correction publishes for the seed audit");
  set_uniform(fixture.solution.view, 1.0);
  set_uniform(fixture.rhs.view, 1.0);
  passed &= expect(static_cast<bool>(
                       fixture.workspace.recycle_begin_projection_for_test(
                           shape, fixture.expected.fingerprint)),
                   rank, "already-converged seed enters the supplemental audit");
  const LinearReductionCounters skipped_before = fixture.reductions.counters();
  LinearSolveResult skipped;
  std::size_t skipped_allocations = 0U;
  {
    allocation_observer::Guard guard;
    skipped = solve_fgmres(identity_operator, identity_preconditioner,
                           invocation(fixture, selected), fixture.workspace,
                           fixture.reductions);
    skipped_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const LinearReductionCounters skipped_work = reduction_delta(
      fixture.reductions.counters(), skipped_before);

  SolveFixture skipped_baseline_fixture;
  const bool skipped_baseline_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, skipped_baseline_fixture);
  passed &= expect(skipped_baseline_initialized, rank,
                   "already-converged plain baseline initializes");
  if (!all_true(passed, communicator)) return false;
  set_uniform(skipped_baseline_fixture.solution.view, 1.0);
  set_uniform(skipped_baseline_fixture.rhs.view, 1.0);
  TridiagonalOperator skipped_baseline_operator(
      communicator, skipped_baseline_fixture.local,
      skipped_baseline_fixture.expected, 0.0, 1.0, 0.0, false);
  ScalingPreconditioner skipped_baseline_preconditioner(
      skipped_baseline_fixture.expected, 1.0, false);
  const LinearReductionCounters skipped_baseline_before =
      skipped_baseline_fixture.reductions.counters();
  LinearSolveResult skipped_baseline;
  std::size_t skipped_baseline_allocations = 0U;
  {
    allocation_observer::Guard guard;
    skipped_baseline = solve_fgmres(
        skipped_baseline_operator, skipped_baseline_preconditioner,
        invocation(skipped_baseline_fixture, selected),
        skipped_baseline_fixture.workspace,
        skipped_baseline_fixture.reductions);
    skipped_baseline_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const LinearReductionCounters skipped_baseline_work = reduction_delta(
      skipped_baseline_fixture.reductions.counters(), skipped_baseline_before);
  passed &= expect(
      skipped.status.code == StatusCode::ok &&
          skipped.recycle_offered_directions == 1U &&
          !skipped.recycle_projection_attempted &&
          skipped.recycle_retained_directions == 0U &&
          skipped.recycle_operator_applies == 0U &&
          skipped.recycle_reduction_calls == 0U &&
          !skipped.recycle_projection_accepted &&
          skipped.iterations == 0U &&
          skipped.status.code == skipped_baseline.status.code &&
          skipped.termination == skipped_baseline.termination &&
          skipped.initial_true_residual ==
              skipped_baseline.initial_true_residual &&
          skipped.final_true_residual == skipped_baseline.final_true_residual &&
          skipped.operator_applies == skipped_baseline.operator_applies &&
          skipped.preconditioner_applies ==
              skipped_baseline.preconditioner_applies &&
          skipped.reduction_calls == skipped_baseline.reduction_calls &&
          skipped_work.calls == skipped_baseline_work.calls &&
          skipped_work.scalars == skipped_baseline_work.scalars &&
          skipped_work.logical_bytes ==
              skipped_baseline_work.logical_bytes &&
          skipped_work.tree_messages ==
              skipped_baseline_work.tree_messages &&
          skipped_work.blocking_operations ==
              skipped_baseline_work.blocking_operations + 1U &&
          same_field_values(fixture.solution.view,
                            skipped_baseline_fixture.solution.view) &&
          skipped_allocations == 0U && skipped_baseline_allocations == 0U,
      rank, "already-converged C2 seed audits without consuming corrections");
  fixture.workspace.recycle_clear_for_test();

  set_uniform(fixture.solution.view, 0.0);
  set_uniform(fixture.rhs.view, 1.0);
  const bool no_cycle_capture = static_cast<bool>(
      fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint));
  const bool no_cycle_projection = static_cast<bool>(
      fixture.workspace.recycle_begin_projection_for_test(
          shape, fixture.expected.fingerprint));
  const LinearSolveResult no_cycle = solve_fgmres(
      identity_operator, identity_preconditioner,
      invocation(fixture, selected), fixture.workspace, fixture.reductions);
  passed &= expect(
      no_cycle_capture && no_cycle_projection && no_cycle.status.code ==
          StatusCode::ok && no_cycle.recycle_offered_directions == 0U &&
          !no_cycle.recycle_projection_attempted &&
          no_cycle.recycle_projected_true_residual == 0.0 &&
          no_cycle.iterations == 1U,
      rank, "empty C1 session records a zero-work C2 projection fallback");
  fixture.workspace.recycle_clear_for_test();

  {
    SolveFixture deflated_fixture;
    SolveFixture deflated_baseline_fixture;
    const bool initialized =
        initialize_fixture(communicator, LinearAlgorithm::fgmres, 4U,
                           deflated_fixture) &&
        initialize_fixture(communicator, LinearAlgorithm::fgmres, 4U,
                           deflated_baseline_fixture);
    passed &= expect(initialized, rank,
                     "exact all-deflated and plain fixtures initialize");
    if (!all_true(passed, communicator)) return false;
    const Int3 deflated_shape{deflated_fixture.local.cells, 1, 1};
    DenseVector zero{};
    DenseVector null_direction{};
    DenseVector solvable_rhs{};
    null_direction[0U] = 1.0;
    solvable_rhs[1U] = 2.0;
    set_partition_values(deflated_fixture.solution.view,
                         deflated_fixture.local, zero);
    const bool captured =
        static_cast<bool>(
            deflated_fixture.workspace.recycle_begin_capture_for_test(
                deflated_shape, deflated_fixture.expected.fingerprint)) &&
        static_cast<bool>(
            deflated_fixture.workspace.recycle_capture_cycle_start_for_test(
                as_const(deflated_fixture.solution.view),
                deflated_fixture.reductions));
    set_partition_values(deflated_fixture.solution.view,
                         deflated_fixture.local, null_direction);
    const bool published = static_cast<bool>(
        deflated_fixture.workspace.recycle_capture_cycle_publish_for_test(
            as_const(deflated_fixture.solution.view),
            deflated_fixture.reductions));
    set_partition_values(deflated_fixture.solution.view,
                         deflated_fixture.local, zero);
    set_partition_values(deflated_fixture.rhs.view, deflated_fixture.local,
                         solvable_rhs);
    const bool projection_begun = static_cast<bool>(
        deflated_fixture.workspace.recycle_begin_projection_for_test(
            deflated_shape, deflated_fixture.expected.fingerprint));
    set_partition_values(deflated_baseline_fixture.solution.view,
                         deflated_baseline_fixture.local, zero);
    set_partition_values(deflated_baseline_fixture.rhs.view,
                         deflated_baseline_fixture.local, solvable_rhs);
    const SolutionSnapshot deflated_seed =
        snapshot_solution(deflated_fixture.solution.view);
    const SolutionSnapshot baseline_seed =
        snapshot_solution(deflated_baseline_fixture.solution.view);
    NullspaceIdentityOperator deflated_operator(deflated_fixture.expected,
                                                deflated_fixture.local);
    NullspaceIdentityOperator deflated_baseline_operator(
        deflated_baseline_fixture.expected,
        deflated_baseline_fixture.local);
    ScalingPreconditioner deflated_preconditioner(
        deflated_fixture.expected, 1.0, false);
    ScalingPreconditioner deflated_baseline_preconditioner(
        deflated_baseline_fixture.expected, 1.0, false);
    const LinearReductionCounters deflated_before =
        deflated_fixture.reductions.counters();
    const LinearReductionCounters deflated_baseline_before =
        deflated_baseline_fixture.reductions.counters();
    LinearSolveResult deflated;
    LinearSolveResult deflated_baseline;
    std::size_t deflated_allocations = 0U;
    std::size_t deflated_baseline_allocations = 0U;
    {
      allocation_observer::Guard guard;
      deflated = solve_fgmres(
          deflated_operator, deflated_preconditioner,
          invocation(deflated_fixture, selected), deflated_fixture.workspace,
          deflated_fixture.reductions);
      deflated_allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    {
      allocation_observer::Guard guard;
      deflated_baseline = solve_fgmres(
          deflated_baseline_operator, deflated_baseline_preconditioner,
          invocation(deflated_baseline_fixture, selected),
          deflated_baseline_fixture.workspace,
          deflated_baseline_fixture.reductions);
      deflated_baseline_allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    const LinearReductionCounters deflated_work = reduction_delta(
        deflated_fixture.reductions.counters(), deflated_before);
    const LinearReductionCounters deflated_baseline_work = reduction_delta(
        deflated_baseline_fixture.reductions.counters(),
        deflated_baseline_before);
    const bool deflated_matrix =
        captured && published && projection_begun &&
            deflated_seed.cells == baseline_seed.cells &&
            deflated_seed.values == baseline_seed.values &&
            deflated.status.code == deflated_baseline.status.code &&
            deflated.termination == deflated_baseline.termination &&
            deflated.iterations == deflated_baseline.iterations &&
            deflated.operator_applies ==
                deflated_baseline.operator_applies &&
            deflated.preconditioner_applies ==
                deflated_baseline.preconditioner_applies &&
            deflated.reduction_calls ==
                deflated_baseline.reduction_calls &&
            deflated.initial_true_residual ==
                deflated_baseline.initial_true_residual &&
            deflated.final_true_residual ==
                deflated_baseline.final_true_residual &&
            same_field_values(deflated_fixture.solution.view,
                              deflated_baseline_fixture.solution.view) &&
            deflated.recycle_offered_directions == 1U &&
            deflated.recycle_retained_directions == 0U &&
            deflated.recycle_projection_attempted &&
            !deflated.recycle_projection_accepted &&
            deflated.recycle_operator_applies == 1U &&
            deflated.recycle_reduction_calls == 1U &&
            deflated.recycle_projected_true_residual == 0.0 &&
            deflated_work.calls == deflated_baseline_work.calls + 1U &&
            deflated_work.scalars ==
                deflated_baseline_work.scalars + 1U &&
            deflated_work.logical_bytes ==
                deflated_baseline_work.logical_bytes + sizeof(double) &&
            deflated_work.tree_messages ==
                deflated_baseline_work.tree_messages &&
            deflated_work.blocking_operations ==
                deflated_baseline_work.blocking_operations + 2U &&
            deflated_allocations == 0U &&
            deflated_baseline_allocations == 0U;
    if (!deflated_matrix) {
      std::cerr
          << "rank " << rank << " all-deflated status="
          << static_cast<unsigned>(deflated.status.code) << '/'
          << static_cast<unsigned>(deflated_baseline.status.code)
          << " term=" << static_cast<unsigned>(deflated.termination) << '/'
          << static_cast<unsigned>(deflated_baseline.termination)
          << " iter=" << deflated.iterations << '/'
          << deflated_baseline.iterations << " offered/retained/attempted="
          << deflated.recycle_offered_directions << '/'
          << deflated.recycle_retained_directions << '/'
          << deflated.recycle_projection_attempted
          << " recycle op/red/proj=" << deflated.recycle_operator_applies
          << '/' << deflated.recycle_reduction_calls << '/'
          << deflated.recycle_projected_true_residual
          << " ordinary op/pc/red=" << deflated.operator_applies << '/'
          << deflated.preconditioner_applies << '/'
          << deflated.reduction_calls << " baseline="
          << deflated_baseline.operator_applies << '/'
          << deflated_baseline.preconditioner_applies << '/'
          << deflated_baseline.reduction_calls << " aggregate delta="
          << (deflated_work.calls - deflated_baseline_work.calls) << '/'
          << (deflated_work.scalars - deflated_baseline_work.scalars) << '/'
          << (deflated_work.logical_bytes -
              deflated_baseline_work.logical_bytes)
          << '/'
          << (deflated_work.blocking_operations -
              deflated_baseline_work.blocking_operations)
          << " alloc=" << deflated_allocations << '/'
          << deflated_baseline_allocations << '\n';
    }
    passed &= expect(
        deflated_matrix, rank,
        "exact current-A2 null images take the all-deflated branch with O=1 work and byte-identical ordinary solve state");
  }

  {
    SolveFixture ordinary_fixture;
    SolveFixture ordinary_baseline_fixture;
    const bool initialized =
        initialize_fixture(communicator, LinearAlgorithm::fgmres, 4U,
                           ordinary_fixture) &&
        initialize_fixture(communicator, LinearAlgorithm::fgmres, 4U,
                           ordinary_baseline_fixture);
    passed &= expect(initialized, rank,
                     "accepted-then-ordinary and plain fixtures initialize");
    if (!all_true(passed, communicator)) return false;
    const Int3 ordinary_shape{ordinary_fixture.local.cells, 1, 1};
    DenseVector zero{};
    DenseVector span_direction{};
    DenseVector mixed_rhs{};
    span_direction[0U] = 1.0;
    mixed_rhs[0U] = 1.0;
    mixed_rhs[1U] = 1.0;
    set_partition_values(ordinary_fixture.solution.view,
                         ordinary_fixture.local, zero);
    const bool captured =
        static_cast<bool>(
            ordinary_fixture.workspace.recycle_begin_capture_for_test(
                ordinary_shape, ordinary_fixture.expected.fingerprint)) &&
        static_cast<bool>(
            ordinary_fixture.workspace.recycle_capture_cycle_start_for_test(
                as_const(ordinary_fixture.solution.view),
                ordinary_fixture.reductions));
    set_partition_values(ordinary_fixture.solution.view,
                         ordinary_fixture.local, span_direction);
    const bool published = static_cast<bool>(
        ordinary_fixture.workspace.recycle_capture_cycle_publish_for_test(
            as_const(ordinary_fixture.solution.view),
            ordinary_fixture.reductions));
    set_partition_values(ordinary_fixture.solution.view,
                         ordinary_fixture.local, zero);
    set_partition_values(ordinary_fixture.rhs.view, ordinary_fixture.local,
                         mixed_rhs);
    const bool projection_begun = static_cast<bool>(
        ordinary_fixture.workspace.recycle_begin_projection_for_test(
            ordinary_shape, ordinary_fixture.expected.fingerprint));
    set_partition_values(ordinary_baseline_fixture.solution.view,
                         ordinary_baseline_fixture.local, zero);
    set_partition_values(ordinary_baseline_fixture.rhs.view,
                         ordinary_baseline_fixture.local, mixed_rhs);
    TridiagonalOperator ordinary_operator(
        communicator, ordinary_fixture.local, ordinary_fixture.expected, 0.0,
        1.0, 0.0, false);
    TridiagonalOperator ordinary_baseline_operator(
        communicator, ordinary_baseline_fixture.local,
        ordinary_baseline_fixture.expected, 0.0, 1.0, 0.0, false);
    ScalingPreconditioner ordinary_preconditioner(
        ordinary_fixture.expected, 1.0, false);
    ScalingPreconditioner ordinary_baseline_preconditioner(
        ordinary_baseline_fixture.expected, 1.0, false);
    AcceptingConvergenceAudit ordinary_audit;
    AcceptingConvergenceAudit ordinary_baseline_audit;
    LinearSolveInvocation ordinary_invocation =
        invocation(ordinary_fixture, selected);
    LinearSolveInvocation ordinary_baseline_invocation =
        invocation(ordinary_baseline_fixture, selected);
    ordinary_invocation.convergence_audit = &ordinary_audit;
    ordinary_baseline_invocation.convergence_audit =
        &ordinary_baseline_audit;
    const LinearReductionCounters ordinary_before =
        ordinary_fixture.reductions.counters();
    const LinearReductionCounters ordinary_baseline_before =
        ordinary_baseline_fixture.reductions.counters();
    LinearSolveResult ordinary;
    LinearSolveResult ordinary_baseline;
    std::size_t ordinary_allocations = 0U;
    std::size_t ordinary_baseline_allocations = 0U;
    {
      allocation_observer::Guard guard;
      ordinary = solve_fgmres(
          ordinary_operator, ordinary_preconditioner, ordinary_invocation,
          ordinary_fixture.workspace, ordinary_fixture.reductions);
      ordinary_allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    {
      allocation_observer::Guard guard;
      ordinary_baseline = solve_fgmres(
          ordinary_baseline_operator, ordinary_baseline_preconditioner,
          ordinary_baseline_invocation, ordinary_baseline_fixture.workspace,
          ordinary_baseline_fixture.reductions);
      ordinary_baseline_allocations =
          allocation_observer::count.load(std::memory_order_relaxed);
    }
    const LinearReductionCounters ordinary_work = reduction_delta(
        ordinary_fixture.reductions.counters(), ordinary_before);
    const LinearReductionCounters ordinary_baseline_work = reduction_delta(
        ordinary_baseline_fixture.reductions.counters(),
        ordinary_baseline_before);
    const bool ordinary_matrix =
        captured && published && projection_begun &&
            ordinary.status.code == StatusCode::ok &&
            ordinary.termination == LinearTermination::converged &&
            ordinary.recycle_offered_directions == 1U &&
            ordinary.recycle_retained_directions == 1U &&
            ordinary.recycle_projection_attempted &&
            ordinary.recycle_projection_accepted &&
            ordinary.recycle_operator_applies == 2U &&
            ordinary.recycle_reduction_calls == 5U &&
            ordinary.iterations > 0U &&
            ordinary.recycle_projected_true_residual <
                ordinary.initial_true_residual &&
            ordinary.final_true_residual !=
                ordinary.recycle_projected_true_residual &&
            ordinary.convergence_audits == 1U &&
            ordinary.convergence_rejections == 0U &&
            ordinary_audit.evaluations() == 1U &&
            ordinary.status.code == ordinary_baseline.status.code &&
            ordinary.termination == ordinary_baseline.termination &&
            ordinary.iterations == ordinary_baseline.iterations &&
            ordinary.operator_applies == ordinary_baseline.operator_applies &&
            ordinary.preconditioner_applies ==
                ordinary_baseline.preconditioner_applies &&
            ordinary.reduction_calls ==
                ordinary_baseline.reduction_calls &&
            same_field_values(ordinary_fixture.solution.view,
                              ordinary_baseline_fixture.solution.view) &&
            ordinary_work.calls == ordinary_baseline_work.calls + 5U &&
            ordinary_work.scalars == ordinary_baseline_work.scalars + 5U &&
            ordinary_work.logical_bytes ==
                ordinary_baseline_work.logical_bytes +
                    5U * sizeof(double) &&
            ordinary_work.tree_messages ==
                ordinary_baseline_work.tree_messages &&
            ordinary_work.blocking_operations ==
                ordinary_baseline_work.blocking_operations + 9U &&
            ordinary_allocations == 0U &&
            ordinary_baseline_allocations == 0U;
    if (!ordinary_matrix) {
      std::cerr
          << "rank " << rank << " accepted-ordinary status="
          << static_cast<unsigned>(ordinary.status.code) << '/'
          << static_cast<unsigned>(ordinary_baseline.status.code)
          << " term=" << static_cast<unsigned>(ordinary.termination) << '/'
          << static_cast<unsigned>(ordinary_baseline.termination)
          << " iter=" << ordinary.iterations << '/'
          << ordinary_baseline.iterations << " init/proj/final="
          << ordinary.initial_true_residual << '/'
          << ordinary.recycle_projected_true_residual << '/'
          << ordinary.final_true_residual << " offered/retained/flags="
          << ordinary.recycle_offered_directions << '/'
          << ordinary.recycle_retained_directions << '/'
          << ordinary.recycle_projection_attempted << '/'
          << ordinary.recycle_projection_accepted << " recycle op/red="
          << ordinary.recycle_operator_applies << '/'
          << ordinary.recycle_reduction_calls << " ordinary op/pc/red="
          << ordinary.operator_applies << '/'
          << ordinary.preconditioner_applies << '/'
          << ordinary.reduction_calls << " baseline="
          << ordinary_baseline.operator_applies << '/'
          << ordinary_baseline.preconditioner_applies << '/'
          << ordinary_baseline.reduction_calls << " audits="
          << ordinary.convergence_audits << '/'
          << ordinary_audit.evaluations() << " aggregate delta="
          << (ordinary_work.calls - ordinary_baseline_work.calls) << '/'
          << (ordinary_work.scalars - ordinary_baseline_work.scalars) << '/'
          << (ordinary_work.logical_bytes -
              ordinary_baseline_work.logical_bytes)
          << '/'
          << (ordinary_work.blocking_operations -
              ordinary_baseline_work.blocking_operations)
          << " alloc=" << ordinary_allocations << '/'
          << ordinary_baseline_allocations << '\n';
    }
    passed &= expect(
        ordinary_matrix, rank,
        "accepted projection removes only the recycle-span component, then canonical audited ordinary FGMRES converges with exact 5/5/40/9 reduction blocking work");
  }

  set_uniform(fixture.solution.view, 0.0);
  const bool fallback_capture =
      static_cast<bool>(fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) &&
      static_cast<bool>(fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  set_uniform(fixture.solution.view, 1.0);
  const bool fallback_publish = static_cast<bool>(
      fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  passed &= expect(fallback_capture && fallback_publish, rank,
                   "non-improving fallback correction captures");
  for (std::int32_t cell = 0; cell < fixture.local.cells; ++cell) {
    fixture.rhs.view.unchecked({cell, 0, 0}, 0U) =
        fixture.local.begin + cell == 0 ? -66.0 : 1.0;
  }
  ++fixture.rhs.view.revision;
  set_uniform(fixture.solution.view, 0.0);
  passed &= expect(static_cast<bool>(
                       fixture.workspace.recycle_begin_projection_for_test(
                           shape, fixture.expected.fingerprint)),
                   rank, "non-improving seed enters projection");
  ScalingPreconditioner normal_preconditioner(fixture.expected, 1.0, false);
  SolveFixture baseline_fixture;
  const bool baseline_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 4U, baseline_fixture);
  passed &= expect(baseline_initialized, rank,
                   "fallback baseline fixture initializes");
  if (!all_true(passed, communicator)) return false;
  for (std::int32_t cell = 0; cell < baseline_fixture.local.cells; ++cell) {
    baseline_fixture.rhs.view.unchecked({cell, 0, 0}, 0U) =
        baseline_fixture.local.begin + cell == 0 ? -66.0 : 1.0;
    baseline_fixture.solution.view.unchecked({cell, 0, 0}, 0U) = 0.0;
  }
  ++baseline_fixture.rhs.view.revision;
  ++baseline_fixture.solution.view.revision;
  TridiagonalOperator baseline_identity_operator(
      communicator, baseline_fixture.local, baseline_fixture.expected, 0.0,
      1.0, 0.0, false);
  ScalingPreconditioner baseline_preconditioner(
      baseline_fixture.expected, 1.0, false);
  const LinearReductionCounters baseline_before =
      baseline_fixture.reductions.counters();
  const LinearReductionCounters fallback_before = fixture.reductions.counters();
  LinearSolveResult baseline;
  LinearSolveResult fallback;
  std::size_t baseline_allocations = 0U;
  std::size_t fallback_allocations = 0U;
  {
    allocation_observer::Guard guard;
    baseline = solve_fgmres(
        baseline_identity_operator, baseline_preconditioner,
        invocation(baseline_fixture, selected), baseline_fixture.workspace,
        baseline_fixture.reductions);
    baseline_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  {
    allocation_observer::Guard guard;
    fallback = solve_fgmres(
        identity_operator, normal_preconditioner,
        invocation(fixture, selected), fixture.workspace, fixture.reductions);
    fallback_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const LinearReductionCounters baseline_work = reduction_delta(
      baseline_fixture.reductions.counters(), baseline_before);
  const LinearReductionCounters fallback_work =
      reduction_delta(fixture.reductions.counters(), fallback_before);
  const bool fallback_matrix =
      baseline.status.code == fallback.status.code &&
          baseline.termination == fallback.termination &&
          baseline.iterations == fallback.iterations &&
          baseline.operator_applies == fallback.operator_applies &&
          baseline.preconditioner_applies == fallback.preconditioner_applies &&
          baseline.reduction_calls == fallback.reduction_calls &&
          baseline.norm_breakdown_restarts == fallback.norm_breakdown_restarts &&
          baseline.initial_true_residual == fallback.initial_true_residual &&
          baseline.final_true_residual == fallback.final_true_residual &&
          same_field_values(fixture.solution.view,
                            baseline_fixture.solution.view) &&
          fallback.recycle_projection_attempted &&
          !fallback.recycle_projection_accepted &&
          fallback.recycle_retained_directions == 1U &&
          fallback.recycle_operator_applies == 2U &&
          fallback.recycle_reduction_calls == 5U &&
          fallback.recycle_projected_true_residual ==
              fallback.initial_true_residual &&
          fallback_work.calls == baseline_work.calls + 5U &&
          fallback_work.scalars == baseline_work.scalars + 5U &&
          fallback_work.logical_bytes ==
              baseline_work.logical_bytes + 5U * sizeof(double) &&
          fallback_work.tree_messages == baseline_work.tree_messages &&
          fallback_work.blocking_operations ==
              baseline_work.blocking_operations + 8U &&
          fallback_work.blocking_operations -
                      baseline_work.blocking_operations +
                  fallback.recycle_operator_applies ==
              10U &&
          baseline_allocations == 0U && fallback_allocations == 0U;
  if (!fallback_matrix) {
    std::cerr << "rank " << rank << " fallback status="
              << static_cast<unsigned>(fallback.status.code) << '/'
              << static_cast<unsigned>(baseline.status.code)
              << " term=" << static_cast<unsigned>(fallback.termination)
              << '/' << static_cast<unsigned>(baseline.termination)
              << " iter=" << fallback.iterations << '/'
              << baseline.iterations << " init/proj/final="
              << fallback.initial_true_residual << '/'
              << fallback.recycle_projected_true_residual << '/'
              << fallback.final_true_residual << " retained/op/red="
              << fallback.recycle_retained_directions << '/'
              << fallback.recycle_operator_applies << '/'
              << fallback.recycle_reduction_calls << " aggregate delta="
              << (fallback_work.calls - baseline_work.calls) << '/'
              << (fallback_work.scalars - baseline_work.scalars) << '/'
              << (fallback_work.logical_bytes - baseline_work.logical_bytes)
              << '/'
              << (fallback_work.blocking_operations -
                  baseline_work.blocking_operations)
              << " alloc=" << fallback_allocations << '/'
              << baseline_allocations << '\n';
  }
  passed &= expect(
      fallback_matrix, rank,
      "finite non-improving projection restores the byte-identical seed with exact 5/5/40 and aggregate 10 blocking operations");
  fixture.workspace.recycle_clear_for_test();

  {
    SolveFixture fault_fixture;
    const bool initialized = initialize_fixture(
        communicator, LinearAlgorithm::fgmres, 4U, fault_fixture);
    passed &= expect(initialized, rank,
                     "QR reduction fault-seam fixture initializes");
    if (!all_true(passed, communicator)) return false;
    const Int3 fault_shape{fault_fixture.local.cells, 1, 1};
    DenseVector fault_zero{};
    DenseVector fault_first{};
    DenseVector fault_second{};
    DenseVector fault_rhs{};
    fault_rhs.fill(1.0);
    fault_first[7U] = 1.0;
    fault_second[22U] = 1.25;
    auto capture_fault_directions = [&]() noexcept {
      DenseVector accepted{};
      bool captured = static_cast<bool>(
          fault_fixture.workspace.recycle_begin_capture_for_test(
              fault_shape, fault_fixture.expected.fingerprint));
      set_partition_values(fault_fixture.solution.view, fault_fixture.local,
                           accepted);
      captured = static_cast<bool>(
                     fault_fixture.workspace
                         .recycle_capture_cycle_start_for_test(
                             as_const(fault_fixture.solution.view),
                             fault_fixture.reductions)) &&
                 captured;
      for (std::size_t cell = 0U; cell < accepted.size(); ++cell) {
        accepted[cell] += fault_first[cell];
      }
      set_partition_values(fault_fixture.solution.view, fault_fixture.local,
                           accepted);
      captured = static_cast<bool>(
                     fault_fixture.workspace
                         .recycle_capture_cycle_publish_for_test(
                             as_const(fault_fixture.solution.view),
                             fault_fixture.reductions)) &&
                 captured;
      set_partition_values(fault_fixture.solution.view, fault_fixture.local,
                           accepted);
      captured = static_cast<bool>(
                     fault_fixture.workspace
                         .recycle_capture_cycle_start_for_test(
                             as_const(fault_fixture.solution.view),
                             fault_fixture.reductions)) &&
                 captured;
      for (std::size_t cell = 0U; cell < accepted.size(); ++cell) {
        accepted[cell] += fault_second[cell];
      }
      set_partition_values(fault_fixture.solution.view, fault_fixture.local,
                           accepted);
      captured = static_cast<bool>(
                     fault_fixture.workspace
                         .recycle_capture_cycle_publish_for_test(
                             as_const(fault_fixture.solution.view),
                             fault_fixture.reductions)) &&
                 captured;
      set_partition_values(fault_fixture.solution.view, fault_fixture.local,
                           fault_zero);
      set_partition_values(fault_fixture.rhs.view, fault_fixture.local,
                           fault_rhs);
      captured = static_cast<bool>(
                     fault_fixture.workspace.recycle_begin_projection_for_test(
                         fault_shape, fault_fixture.expected.fingerprint)) &&
                 captured;
      return captured;
    };
    const bool fault_captured = capture_fault_directions();
    const SolutionSnapshot caller_before =
        snapshot_solution(fault_fixture.solution.view);
    const int failing_rank = fault_fixture.local.size - 1;
    const Status armed =
        fault_fixture.reductions.arm_checked_sum_fault_for_test(5U,
                                                                 failing_rank);
    TridiagonalOperator fault_operator(
        communicator, fault_fixture.local, fault_fixture.expected, 0.0, 1.0,
        0.0, false);
    ScalingPreconditioner fault_preconditioner(fault_fixture.expected, 1.0,
                                               false);
    const LinearSolveControl fault_control = control(20U, 4U);
    const LinearSolveResult fault_result = solve_fgmres(
        fault_operator, fault_preconditioner,
        invocation(fault_fixture, fault_control), fault_fixture.workspace,
        fault_fixture.reductions);
    passed &= expect(
        fault_captured && static_cast<bool>(armed) &&
            fault_result.status.code == StatusCode::numerical_failure &&
            fault_result.lowest_failing_rank == failing_rank &&
            fault_result.recycle_projection_attempted &&
            !fault_result.recycle_projection_accepted &&
            same_solution(fault_fixture.solution.view, caller_before),
        rank,
        "one-shot QR checked-sum fault fails together at the true lowest rank without publishing x");

    // The failed projection must consume its private session itself.  Rebuild
    // immediately, without a test clear or re-arm, and prove the same
    // workspace can complete an ordinary projection/solve.
    const bool retry_captured = capture_fault_directions();
    const LinearSolveResult retry_result = solve_fgmres(
        fault_operator, fault_preconditioner,
        invocation(fault_fixture, fault_control), fault_fixture.workspace,
        fault_fixture.reductions);
    passed &= expect(
        retry_captured && retry_result.status.code == StatusCode::ok &&
            retry_result.termination == LinearTermination::converged &&
            retry_result.recycle_projection_attempted &&
            retry_result.recycle_projection_accepted,
        rank, "consumed QR fault arm does not poison the next private session");
    fault_fixture.reductions.clear_checked_sum_fault_for_test();

    const auto projection_operator_failure_probe =
        [&](std::string_view description) noexcept {
          const bool captured = capture_fault_directions();
          const SolutionSnapshot before =
              snapshot_solution(fault_fixture.solution.view);
          TridiagonalOperator failing_operator(
              communicator, fault_fixture.local, fault_fixture.expected, 0.0,
              1.0, 0.0, false, failing_rank, 2U);
          ScalingPreconditioner healthy_preconditioner(
              fault_fixture.expected, 1.0, false);
          const LinearSolveResult failed = solve_fgmres(
              failing_operator, healthy_preconditioner,
              invocation(fault_fixture, fault_control),
              fault_fixture.workspace, fault_fixture.reductions);
          const Status stale_projection =
              fault_fixture.workspace.recycle_begin_projection_for_test(
                  fault_shape, fault_fixture.expected.fingerprint);
          const bool failed_solution_unchanged =
              same_solution(fault_fixture.solution.view, before);
          const bool retry_captured = capture_fault_directions();
          TridiagonalOperator retry_operator(
              communicator, fault_fixture.local, fault_fixture.expected, 0.0,
              1.0, 0.0, false);
          ScalingPreconditioner retry_preconditioner(
              fault_fixture.expected, 1.0, false);
          const LinearSolveResult retry = solve_fgmres(
              retry_operator, retry_preconditioner,
              invocation(fault_fixture, fault_control),
              fault_fixture.workspace, fault_fixture.reductions);
          const bool probe_passed =
              captured &&
              failed.status.code == StatusCode::numerical_failure &&
              failed.termination == LinearTermination::operator_failure &&
              failed.lowest_failing_rank == failing_rank &&
              failed.recycle_offered_directions == 2U &&
              failed.recycle_projection_attempted &&
              !failed.recycle_projection_accepted &&
              failed.recycle_operator_applies == 1U &&
              failed_solution_unchanged &&
              stale_projection.code == StatusCode::invalid_plan &&
              retry_captured && retry.status.code == StatusCode::ok &&
              retry.termination == LinearTermination::converged &&
              retry.recycle_projection_attempted &&
              retry.recycle_projection_accepted;
          return expect(probe_passed, rank, description);
        };
    passed &= projection_operator_failure_probe(
        "rank-selective current-A2 failure is collective, transactional, single-use and immediately retryable");
  }

  set_uniform(fixture.solution.view, 0.0);
  const bool nonfinite_capture =
      static_cast<bool>(fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) &&
      static_cast<bool>(fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  set_uniform(fixture.solution.view, std::numeric_limits<double>::quiet_NaN());
  const Status nonfinite_publish =
      fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(fixture.solution.view), fixture.reductions);
  passed &= expect(nonfinite_capture &&
                       nonfinite_publish.detail != 0U &&
                       nonfinite_publish.code != StatusCode::ok,
                   rank, "nonfinite capture fails collectively before offering a direction");

  // No explicit clear: the failed publish must leave the workspace ready for
  // a new capture on the same resources.
  set_uniform(fixture.solution.view, 0.0);
  const bool rank_selective_capture =
      static_cast<bool>(fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) &&
      static_cast<bool>(fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  if (fixture.local.rank == fixture.local.size - 1) {
    set_uniform(fixture.solution.view, std::numeric_limits<double>::quiet_NaN());
  }
  const Status rank_selective_publish =
      fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(fixture.solution.view), fixture.reductions);
  passed &= expect(
      rank_selective_capture && rank_selective_publish.code != StatusCode::ok &&
          rank_selective_publish.detail != 0U &&
          fixture.workspace.recycle_correction_count_for_test() == 0U,
      rank, "rank-selective nonfinite capture fails together without offering work");

  set_uniform(fixture.solution.view, 0.0);
  const bool capture_retry_begun =
      static_cast<bool>(fixture.workspace.recycle_begin_capture_for_test(
          shape, fixture.expected.fingerprint)) &&
      static_cast<bool>(fixture.workspace.recycle_capture_cycle_start_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  set_uniform(fixture.solution.view, 1.0);
  const bool capture_retry_published = static_cast<bool>(
      fixture.workspace.recycle_capture_cycle_publish_for_test(
          as_const(fixture.solution.view), fixture.reductions));
  set_uniform(fixture.solution.view, 0.0);
  set_uniform(fixture.rhs.view, 1.0);
  const bool capture_retry_projection = static_cast<bool>(
      fixture.workspace.recycle_begin_projection_for_test(
          shape, fixture.expected.fingerprint));
  const LinearSolveResult capture_retry = solve_fgmres(
      identity_operator, identity_preconditioner,
      invocation(fixture, selected), fixture.workspace, fixture.reductions);
  passed &= expect(
      capture_retry_begun && capture_retry_published &&
          capture_retry_projection && capture_retry.status.code == StatusCode::ok &&
          capture_retry.termination == LinearTermination::converged &&
          capture_retry.recycle_projection_attempted &&
          capture_retry.recycle_projection_accepted,
      rank,
      "rank-selective capture failure is self-consumed and the divergent revision path immediately retries on the same workspace");
  return all_true(passed, communicator);
}

bool test_repeated_hot_solves_allocate_zero(MPI_Comm communicator, int rank) {
  SolveFixture fixture;
  bool passed = expect(initialize_fixture(communicator, LinearAlgorithm::pcg,
                                          0U, fixture),
                       rank, "repeated solve fixture initializes");
  if (!all_true(passed, communicator)) {
    return false;
  }
  TridiagonalOperator op(communicator, fixture.local, fixture.expected, -1.0,
                         2.0, -1.0, true);
  ScalingPreconditioner preconditioner(fixture.expected, 0.5, true);
  const std::uintptr_t vector_address =
      fixture.workspace.vector_storage_address();
  const std::uintptr_t scalar_address =
      fixture.workspace.scalar_storage_address();
  bool hot_ok = true;
  std::size_t allocations = 0U;
  {
    allocation_observer::Guard guard;
    for (std::size_t repetition = 0U; repetition < 100U; ++repetition) {
      fill_system(fixture, -1.0, 2.0, -1.0);
      const LinearSolveResult result =
          solve_pcg(op, preconditioner, invocation(fixture, control()),
                    fixture.workspace, fixture.reductions);
      hot_ok = result.status.code == StatusCode::ok &&
               fixture.workspace.vector_storage_address() == vector_address &&
               fixture.workspace.scalar_storage_address() == scalar_address &&
               finite_solution(fixture) && hot_ok;
    }
    allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(hot_ok && allocations == 0U, rank,
                   "100 repeated PCG solves preserve workspace and allocate zero heap bytes");
  const ResidualOracle oracle = independent_true_residual(
      fixture, communicator, -1.0, 2.0, -1.0);
  const double error = global_error(fixture, communicator);
  passed &= expect(residual_is_accepted(oracle, control()) && error < 1.0e-8,
                   rank,
                   "hot-loop allocation guard is followed by an independent residual oracle");

  SolveFixture fgmres_fixture;
  const bool fgmres_initialized = initialize_fixture(
      communicator, LinearAlgorithm::fgmres, 12U, fgmres_fixture);
  passed &= expect(fgmres_initialized, rank,
                   "repeated FGMRES solve fixture initializes");
  if (!all_true(fgmres_initialized, communicator)) return false;
  TridiagonalOperator fgmres_op(communicator, fgmres_fixture.local,
                                fgmres_fixture.expected, -1.2, 3.0, -0.7,
                                false);
  ScalingPreconditioner fgmres_preconditioner(
      fgmres_fixture.expected, 1.0 / 3.0, false);
  const std::uintptr_t fgmres_vector_address =
      fgmres_fixture.workspace.vector_storage_address();
  const std::uintptr_t fgmres_scalar_address =
      fgmres_fixture.workspace.scalar_storage_address();
  bool fgmres_hot_ok = true;
  std::size_t fgmres_allocations = 0U;
  {
    allocation_observer::Guard guard;
    for (std::size_t repetition = 0U; repetition < 32U; ++repetition) {
      fill_system(fgmres_fixture, -1.2, 3.0, -0.7);
      const LinearSolveResult result = solve_fgmres(
          fgmres_op, fgmres_preconditioner,
          invocation(fgmres_fixture, control(250U, 12U)),
          fgmres_fixture.workspace, fgmres_fixture.reductions);
      fgmres_hot_ok =
          result.status.code == StatusCode::ok &&
          result.termination == LinearTermination::converged &&
          result.norm_breakdown_restarts == 0U &&
          fgmres_fixture.workspace.vector_storage_address() ==
              fgmres_vector_address &&
          fgmres_fixture.workspace.scalar_storage_address() ==
              fgmres_scalar_address &&
          finite_solution(fgmres_fixture) && fgmres_hot_ok;
    }
    fgmres_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(
      fgmres_hot_ok && fgmres_allocations == 0U, rank,
      "32 repeated FGMRES solves preserve workspace and allocate zero heap bytes");
  const ResidualOracle fgmres_oracle = independent_true_residual(
      fgmres_fixture, communicator, -1.2, 3.0, -0.7);
  const double fgmres_error = global_error(fgmres_fixture, communicator);
  passed &= expect(
      residual_is_accepted(fgmres_oracle, control(250U, 12U)) &&
          fgmres_error < 1.0e-8,
      rank,
      "FGMRES hot-loop allocation guard is followed by an independent residual oracle");
  return all_true(passed, communicator);
}

bool test_fgmres_basis_workspace_nonoverlap(MPI_Comm communicator, int rank) {
  constexpr std::uint32_t restart = 30U;
  SolveFixture fixture;
  bool passed = expect(initialize_fixture(communicator,
                                          LinearAlgorithm::fgmres, restart,
                                          fixture),
                       rank, "restart-30 basis-layout fixture initializes");
  if (!all_true(passed, communicator)) return false;

  const Int3 shape{fixture.local.cells, 1, 1};
  const std::uint8_t ax_slot = static_cast<std::uint8_t>(2U * restart + 2U);
  FieldView destination = fixture.workspace.vector(ax_slot, shape);
  std::array<ConstFieldView, restart + 1U> sources{};
  for (std::size_t row = 0U; row < sources.size(); ++row) {
    sources[row] = as_const(fixture.workspace.vector(
        static_cast<std::uint8_t>(1U + row), shape));
  }
  passed &= expect(
      detail::krylov_basis_update_inputs_disjoint_for_test(
          destination, sources.data(), sources.size()),
      rank,
      "restart-30 Arnoldi basis slots are address-disjoint from ax storage");
  sources[0U] = as_const(destination);
  passed &= expect(
      !detail::krylov_basis_update_inputs_disjoint_for_test(
          destination, sources.data(), sources.size()),
      rank, "basis-layout authority detects an overlapping destination slot");
  return all_true(passed, communicator);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 1;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  bool passed = size == 1 || size == 2 || size == 4;
  passed &= test_pcg_spd_and_zero_rhs(MPI_COMM_WORLD, rank);
  passed &= test_nonsymmetric_solvers(MPI_COMM_WORLD, rank);
  passed &= test_rejections_are_transactional(MPI_COMM_WORLD, rank, size);
  passed &= test_breakdown_maxiter_and_stale_identity(MPI_COMM_WORLD, rank);
  passed &= test_fgmres_norm_breakdown_lifecycle(MPI_COMM_WORLD, rank);
  passed &= test_fgmres_rank_selective_failures(MPI_COMM_WORLD, rank, size);
  passed &= test_operator_provenance_matching(MPI_COMM_WORLD, rank);
  passed &= test_fgmres_collective_status_scope(MPI_COMM_WORLD, rank, size);
  passed &= test_fgmres_prepared_batch_lifecycle(MPI_COMM_WORLD, rank, size);
  passed &= test_fgmres_reduction_capacity(MPI_COMM_WORLD, rank);
  passed &=
      test_fgmres_supplemental_convergence_audit(MPI_COMM_WORLD, rank);
  passed &= test_c1_residual_guarded_warm_start(MPI_COMM_WORLD, rank);
  passed &= test_recycle_capture_projection(MPI_COMM_WORLD, rank);
  passed &= test_fgmres_basis_workspace_nonoverlap(MPI_COMM_WORLD, rank);
  passed &= test_repeated_hot_solves_allocate_zero(MPI_COMM_WORLD, rank);
  const bool global = all_true(passed, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
