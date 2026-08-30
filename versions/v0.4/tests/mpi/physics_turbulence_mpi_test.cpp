// SPDX-License-Identifier: Apache-2.0

#include "../support/turbulence_fixture.hpp"

#include <mpi.h>

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <new>
#include <string_view>

namespace allocation_observer {
std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};
void* allocate(std::size_t size) {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
  if (void* result = std::malloc(size == 0U ? 1U : size)) return result;
  throw std::bad_alloc{};
}
class Guard {
 public:
  Guard() {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }
};
}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}
void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
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
void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local) {
  const int value = local ? 1 : 0;
  int global = 0;
  return MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         global != 0;
}

TurbulenceOwnedField gradient_field(const TurbulenceFixture& fixture,
                                    FieldId field, RevisionToken revision,
                                    StorageIdentity storage) {
  TurbulenceOwnedField result =
      make_turbulence_field(field, fixture.patch.cells, revision, storage);
  const std::size_t scalar_size = result.storage.size();
  result.storage.resize(9U * scalar_size, 0.0);
  result.view.base = result.storage.data();
  result.view.components = 9U;
  result.view.component_stride = scalar_size;
  for (std::size_t cell = 0U; cell < fixture.gradients.size(); ++cell) {
    const std::int32_t x = static_cast<std::int32_t>(
        cell % static_cast<std::size_t>(fixture.patch.cells.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(fixture.patch.cells.x);
    const std::int32_t y = static_cast<std::int32_t>(
        yz % static_cast<std::size_t>(fixture.patch.cells.y));
    const std::int32_t z = static_cast<std::int32_t>(
        yz / static_cast<std::size_t>(fixture.patch.cells.y));
    for (std::uint8_t component = 0U; component < 9U; ++component) {
      result.view.unchecked({x, y, z}, component) =
          fixture.gradients[cell].value[component];
    }
  }
  return result;
}

bool test_collective_binding_and_hot_update(int rank, int size) {
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(static_cast<bool>(CartesianGeometryCompiler::compile(
                           MPI_COMM_WORLD,
                           turbulence_mesh(GeometryKind::uniform, 5), {},
                           geometry, patch)),
                       rank, "non-divisible turbulence geometry compiles");
  const std::array<FieldId, 3U> fields{0U, 1U, 2U};
  ContributionRegistry divergent_registry;
  passed &= expect(static_cast<bool>(divergent_registry.configure(
                       {fields.data(), fields.size()})),
                   rank, "divergent registry configures");
  TurbulencePlanSpec divergent_spec;
  if (size > 1 && rank == 1) {
    divergent_spec.kind = TurbulenceKind::wale;
  }
  TurbulencePlan rejected;
  const Status divergent = TurbulencePlan::compile(
      MPI_COMM_WORLD, divergent_spec, geometry, patch, 2U, 41U,
      divergent_registry, rejected);
  if (size > 1) {
    passed &= expect(divergent.code == StatusCode::invalid_plan &&
                         rejected.fingerprint() == 0U,
                     rank, "rank-divergent model selection rejects collectively");
  }

  TurbulenceFixture fixture;
  passed &= expect(fixture.initialize(TurbulencePlanSpec{},
                                      GeometryKind::uniform, MPI_COMM_WORLD,
                                      5),
                   rank, "default distributed Vreman plan compiles");
  if (!all_true(passed)) return false;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  const std::uint64_t fingerprint = fixture.plan.fingerprint();
  MPI_Allreduce(&fingerprint, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&fingerprint, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(minimum == maximum && minimum != 0U, rank,
                   "1/2/4 ranks share one turbulence semantic plan");

  for (VelocityGradient& gradient : fixture.gradients) {
    gradient.value[0U] = 0.7;
    gradient.value[4U] = -0.4;
    gradient.value[8U] = -0.3;
    gradient.value[1U] = 0.2;
  }
  TurbulenceCertificate certificate;
  passed &= expect(static_cast<bool>(fixture.plan.update(
                       fixture.input(), fixture.effective.view, certificate)),
                   rank, "distributed Vreman field updates");
  Status hot;
  {
    allocation_observer::Guard guard;
    for (RevisionToken revision = 200U; revision < 264U; ++revision) {
      const Status status = fixture.plan.update(
          fixture.input(revision), fixture.effective.view, certificate);
      if (!status && hot) hot = status;
    }
  }
  passed &= expect(static_cast<bool>(hot) &&
                       allocation_observer::count.load(
                           std::memory_order_relaxed) == 0U &&
                       fixture.plan.update_count() == 65U,
                   rank, "repeated turbulence updates are allocation-free");
  for (double value : fixture.effective.storage) {
    passed &= expect(std::isfinite(value) && value > 1.8e-5, rank,
                     "every rank publishes finite positive mu_eff");
  }

  const std::vector<double> live_before = fixture.effective.storage;
  const std::uint64_t count_before = fixture.plan.update_count();
  TurbulenceOwnedField candidate_density = make_turbulence_field(
      10U, fixture.patch.cells, 301U,
      static_cast<StorageIdentity>(401U + rank));
  TurbulenceOwnedField candidate_molecular = make_turbulence_field(
      11U, fixture.patch.cells, 302U,
      static_cast<StorageIdentity>(501U + rank));
  std::copy(fixture.density.storage.begin(), fixture.density.storage.end(),
            candidate_density.storage.begin());
  std::copy(fixture.molecular.storage.begin(), fixture.molecular.storage.end(),
            candidate_molecular.storage.begin());
  TurbulenceOwnedField candidate_gradient = gradient_field(
      fixture, 12U, 303U, static_cast<StorageIdentity>(601U + rank));
  TurbulenceOwnedField candidate_effective = make_turbulence_field(
      13U, fixture.patch.cells, 304U,
      static_cast<StorageIdentity>(701U + rank));
  const TurbulenceCandidateInput candidate_input{
      as_const(candidate_density.view), as_const(candidate_molecular.view),
      as_const(candidate_gradient.view), candidate_gradient.view.revision};
  TurbulenceCandidateCertificate candidate_certificate;
  Status candidate_status;
  std::size_t candidate_allocations = 0U;
  {
    allocation_observer::Guard guard;
    candidate_status = fixture.plan.evaluate_candidate_effective_viscosity(
        candidate_input, candidate_effective.view, candidate_certificate);
    candidate_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(candidate_status && candidate_certificate.valid() &&
                       candidate_allocations == 0U &&
                       candidate_effective.storage == live_before &&
                       fixture.plan.update_count() == count_before,
                   rank,
                   "distributed candidate evaluation is exact, stateless, "
                   "and allocation-free");

  std::fill(candidate_effective.storage.begin(),
            candidate_effective.storage.end(), -19.0);
  const std::vector<double> poison_sentinel = candidate_effective.storage;
  const int poison_rank = size - 1;
  const double saved_density = candidate_density.storage.back();
  if (rank == poison_rank) {
    candidate_density.storage.back() =
        std::numeric_limits<double>::quiet_NaN();
  }
  candidate_status = fixture.plan.evaluate_candidate_effective_viscosity(
      candidate_input, candidate_effective.view, candidate_certificate);
  const int local_failure = candidate_status ? 0 : 1;
  int any_failure = 0;
  const bool consensus =
      MPI_Allreduce(&local_failure, &any_failure, 1, MPI_INT, MPI_MAX,
                    MPI_COMM_WORLD) == MPI_SUCCESS;
  passed &= expect(
      consensus && any_failure == 1 &&
          (rank != poison_rank ||
           (!candidate_status && !candidate_certificate.valid() &&
            candidate_effective.storage == poison_sentinel)) &&
          fixture.effective.storage == live_before &&
          fixture.plan.update_count() == count_before,
      rank,
      "single-rank candidate poison is collectively skippable without live "
      "writes");
  if (rank == poison_rank) {
    candidate_density.storage.back() = saved_density;
  }
  passed &= expect(
      fixture.plan.evaluate_candidate_effective_viscosity(
          candidate_input, candidate_effective.view, candidate_certificate) &&
          candidate_effective.storage == live_before &&
          fixture.effective.storage == live_before &&
          fixture.plan.update_count() == count_before,
      rank, "clean candidate evaluation recovers after distributed poison");
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const bool passed = test_collective_binding_and_hot_update(rank, size);
  MPI_Finalize();
  return passed ? 0 : 1;
}
