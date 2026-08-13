// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"

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

class ScalingPreconditioner final : public LinearPreconditioner {
 public:
  ScalingPreconditioner(LinearIdentity expected, double scale, bool spd,
                        bool variable = false,
                        std::int32_t global_begin = 0, int local_rank = 0,
                        int failing_rank = -1,
                        std::uint32_t failing_call = 0U,
                        FlexibleApplicationTrace* flexible_trace = nullptr) noexcept
      : expected_(expected),
        scale_(scale),
        global_begin_(global_begin),
        local_rank_(local_rank),
        failing_rank_(failing_rank),
        failing_call_(failing_call),
        flexible_trace_(flexible_trace),
        spd_(spd),
        variable_(variable) {}

  LinearPreconditionerCertificate certificate() const noexcept override {
    return {expected_, 0x87654321U,
            spd_ ? LinearPreconditionerClass::fixed_spd
                 : (variable_ ? LinearPreconditionerClass::flexible
                              : LinearPreconditionerClass::fixed_general)};
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
  std::uint32_t calls_{};
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
    const LinearSolveControl selected = control(250U, 12U);
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
                         residual_report_matches(result.final_true_residual,
                                                 oracle) &&
                         residual_is_accepted(oracle, selected) &&
                         finite_solution(fixture) && error < 1.0e-8 &&
                         flexible.calls() == result.iterations &&
                         flexible.calls() == result.preconditioner_applies &&
                         op.calls() == result.operator_applies &&
                         trace.observed_operator_inputs == flexible.calls() &&
                         trace.mismatched_operator_inputs == 0U &&
                         trace.overwritten_outputs == 0U && !trace.pending &&
                         trace.saw_spatial_variation &&
                         trace.saw_temporal_variation && calls_match,
                     rank,
                     "right-preconditioned FGMRES supports a flexible preconditioner");
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
  passed &= test_repeated_hot_solves_allocate_zero(MPI_COMM_WORLD, rank);
  const bool global = all_true(passed, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
