// SPDX-License-Identifier: Apache-2.0

#include "../support/ibm_force_fixture.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>

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
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

bool close(double left, double right, double tolerance = 2.0e-10) {
  return std::abs(left - right) <= tolerance;
}

bool test_decomposition_and_cache(int rank, int size) {
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(MPI_COMM_WORLD), rank,
                       "distributed final-force fixture compiles");
  if (!all_true(passed)) return false;
  fixture.fill_analytic();
  const FinalSurfaceState state = fixture.state(fixture.committed_flux.flux());
  SurfaceForce force;
  const Status evaluated = evaluate_surface_force(
      MPI_COMM_WORLD, fixture.quadrature, state, force);
  passed &= expect(
      static_cast<bool>(evaluated) && close(force.pressure.x, -0.7) &&
          close(force.pressure.y, 0.4) && close(force.pressure.z, -0.2) &&
          close(force.viscous.x, 0.0) && close(force.viscous.y, 0.0) &&
          close(force.viscous.z, 0.0) && close(force.moment.x, 0.0) &&
          close(force.moment.y, 0.0) && close(force.moment.z, 0.0),
      rank, "1/2/4-rank quadratic cube force matches divergence-theorem oracle");
  std::uint64_t minimum = force.revision;
  std::uint64_t maximum = force.revision;
  MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(minimum == maximum && minimum != 0U, rank,
                   "all ranks publish one final-force revision");

  FieldRegistry registry;
  FieldSchema schema;
  FieldId dependency = 0U;
  const std::array requests{ArenaFieldRequest{
      dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  StateLayers layers;
  AttemptTransaction transaction;
  FinalForceCache cache;
  passed &= expect(registry.declare_field("force.mpi", 1U, 0U, dependency) &&
                       registry.freeze(schema) &&
                       ArenaLayout::compile(
                           schema, {requests.data(), requests.size()}, layout) &&
                       StateLayers::allocate(layout, layers) &&
                       AttemptTransaction::create(
                           layers.field_count(), 1U, layers.field_count(),
                           transaction) &&
                       FinalForceCache::bind(
                           MPI_COMM_WORLD, fixture.quadrature,
                           fixture.geometry.topology_revision(), 91U, 0U,
                           cache),
                   rank, "distributed final-force cache binds");
  if (!all_true(passed)) return false;

  const auto attempt = [&](bool fail, SurfaceForce& committed) {
    FinalForceCertificate prepared;
    Status status = transaction.begin(layers);
    if (status) status = transaction.revise_trial(dependency);
    const RevisionDependency stamp{
        AttemptTransaction::field_revision_source(dependency),
        transaction.trial_revision(dependency)};
    const std::array dependencies{stamp};
    if (status) {
      status = cache.prepare(state, {dependencies.data(), dependencies.size()},
                             transaction, prepared);
    }
    if (status) {
      const Status outcome =
          fail && ((size == 1 && rank == 0) || (size > 1 && rank == 1))
              ? Status{StatusCode::rejected_step, 92U}
              : Status{};
      status = transaction.collective_finish(MPI_COMM_WORLD, outcome);
    }
    if (transaction.finished()) {
      const Status finalized = cache.finalize(transaction);
      if (status && !finalized) status = finalized;
    }
    FinalForceCertificate certificate;
    const Status visible = cache.committed(committed, certificate);
    return std::array<Status, 2U>{status, visible};
  };

  SurfaceForce committed;
  const auto rejected = attempt(true, committed);
  passed &= expect(rejected[0U].code == StatusCode::rejected_step &&
                       rejected[1U].code == StatusCode::invalid_plan,
                   rank, "one-rank failed attempt leaves force cache invisible everywhere");
  const auto accepted = attempt(false, committed);
  passed &= expect(static_cast<bool>(accepted[0U]) &&
                       static_cast<bool>(accepted[1U]) &&
                       close(committed.pressure.x, -0.7) &&
                       close(committed.pressure.y, 0.4) &&
                       close(committed.pressure.z, -0.2),
                   rank, "collective success publishes the exact final force once");
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const bool passed = test_decomposition_and_cache(rank, size);
  MPI_Finalize();
  return passed ? 0 : 1;
}
