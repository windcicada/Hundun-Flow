// SPDX-License-Identifier: Apache-2.0

#include "../support/ibm_force_fixture.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string_view>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) std::cerr << "FAIL: " << description << '\n';
  return condition;
}

double difference(const SurfaceForce& left, const SurfaceForce& right) {
  return std::abs(left.pressure.x - right.pressure.x) +
         std::abs(left.pressure.y - right.pressure.y) +
         std::abs(left.pressure.z - right.pressure.z) +
         std::abs(left.viscous.x - right.viscous.x) +
         std::abs(left.viscous.y - right.viscous.y) +
         std::abs(left.viscous.z - right.viscous.z) +
         std::abs(left.moment.x - right.moment.x) +
         std::abs(left.moment.y - right.moment.y) +
         std::abs(left.moment.z - right.moment.z);
}

bool test_revision_and_flux_mutations() {
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(), "force mutation fixture compiles");
  if (!passed) return false;
  fixture.fill_analytic();
  const FinalSurfaceState baseline_state =
      fixture.state(fixture.committed_flux.flux());
  SurfaceForce baseline;
  passed &= expect(static_cast<bool>(evaluate_surface_force(
                       MPI_COMM_SELF, fixture.quadrature, baseline_state,
                       baseline)),
                   "baseline final-state force evaluates");
  const SurfaceForce marker = baseline;

  FinalSurfaceState mutation = baseline_state;
  ++mutation.velocity_gradient.revision;
  SurfaceForce rejected = marker;
  passed &= expect(evaluate_surface_force(MPI_COMM_SELF, fixture.quadrature,
                                          mutation, rejected)
                           .code == StatusCode::invalid_plan &&
                       difference(rejected, marker) == 0.0,
                   "corrector-1/stale gradient mutation is RED and atomic");

  mutation = baseline_state;
  ++mutation.final_velocity.revision;
  rejected = marker;
  passed &= expect(evaluate_surface_force(MPI_COMM_SELF, fixture.quadrature,
                                          mutation, rejected)
                           .code == StatusCode::invalid_plan &&
                       difference(rejected, marker) == 0.0,
                   "stale final-U mutation is RED and atomic");

  mutation = baseline_state;
  ++mutation.effective_viscosity.revision;
  rejected = marker;
  passed &= expect(evaluate_surface_force(MPI_COMM_SELF, fixture.quadrature,
                                          mutation, rejected)
                           .code == StatusCode::invalid_plan &&
                       difference(rejected, marker) == 0.0,
                   "stale mu_eff mutation is RED and atomic");

  FaceFluxStorage provisional_storage;
  FaceFluxView provisional;
  passed &= expect(FaceFluxStorage::allocate_workspace(
                       fixture.patch.cells, 1U, provisional_storage) &&
                       provisional_storage.workspace_view(
                           0U, baseline_state.final_flux, provisional),
                   "provisional flux fixture allocates");
  mutation = baseline_state;
  mutation.face_flux = as_const(provisional);
  rejected = marker;
  passed &= expect(evaluate_surface_force(MPI_COMM_SELF, fixture.quadrature,
                                          mutation, rejected)
                           .code == StatusCode::invalid_plan &&
                       difference(rejected, marker) == 0.0,
                   "provisional-flux substitution is RED and atomic");
  return passed;
}

bool test_sign_area_normal_mutations() {
  SurfaceTractionPoint point;
  point.position = {0.7, -0.2, 0.4};
  point.solid_to_fluid_normal = {1.0, 0.0, 0.0};
  point.weight = 2.5;
  point.absolute_pressure = 3.0;
  point.effective_viscosity = 0.02;
  point.velocity_gradient.value =
      {0.2, 0.3, 0.0, -0.1, -0.05, 0.2, 0.0, 0.1, -0.15};
  const std::array baseline_points{point};
  SurfaceForce baseline;
  bool passed = expect(static_cast<bool>(integrate_surface_traction(
                           {baseline_points.data(), baseline_points.size()},
                           {}, baseline)),
                       "single-point traction baseline evaluates");

  auto mutated = baseline_points;
  mutated[0U].solid_to_fluid_normal.x = -1.0;
  SurfaceForce result;
  passed &= expect(integrate_surface_traction(
                       {mutated.data(), mutated.size()}, {}, result) &&
                       difference(result, baseline) > 1.0,
                   "reversed-normal mutation is RED");
  mutated = baseline_points;
  mutated[0U].weight = 1.0;
  passed &= expect(integrate_surface_traction(
                       {mutated.data(), mutated.size()}, {}, result) &&
                       difference(result, baseline) > 1.0,
                   "omitted-area mutation is RED");
  mutated = baseline_points;
  mutated[0U].absolute_pressure *= -1.0;
  passed &= expect(integrate_surface_traction(
                       {mutated.data(), mutated.size()}, {}, result) &&
                       difference(result, baseline) > 1.0,
                   "wrong-pressure-sign mutation is RED");
  return passed;
}

bool test_force_cache_rollback_and_commit() {
  IbmForceFixture fixture;
  bool passed = expect(fixture.initialize(), "force cache fixture compiles");
  if (!passed) return false;
  fixture.fill_analytic();
  const FinalSurfaceState state = fixture.state(fixture.committed_flux.flux());

  FieldRegistry registry;
  FieldSchema schema;
  FieldId dependency = 0U;
  passed &= expect(registry.declare_field("force.dependency", 1U, 0U,
                                         dependency) &&
                       registry.freeze(schema),
                   "force cache transaction schema freezes");
  const std::array requests{ArenaFieldRequest{
      dependency, {1, 1, 1}, {0U}, FieldLifetime::state_layer}};
  ArenaLayout layout;
  StateLayers layers;
  AttemptTransaction transaction;
  FinalForceCache cache;
  passed &= expect(ArenaLayout::compile(schema,
                                        {requests.data(), requests.size()},
                                        layout) &&
                       StateLayers::allocate(layout, layers) &&
                       AttemptTransaction::create(layers.field_count(), 1U,
                                                  layers.field_count(),
                                                  transaction) &&
                       FinalForceCache::bind(
                           MPI_COMM_SELF, fixture.quadrature,
                           fixture.geometry.topology_revision(), 88U, 0U,
                           cache),
                   "force cache binds cold transaction resources");

  const auto run_attempt = [&](Status outcome,
                               FinalForceCertificate& certificate) {
    if (!transaction.begin(layers) || !transaction.revise_trial(dependency)) {
      return Status{StatusCode::invalid_plan, 1U};
    }
    const RevisionDependency dependency_stamp{
        AttemptTransaction::field_revision_source(dependency),
        transaction.trial_revision(dependency)};
    const std::array dependencies{dependency_stamp};
    Status status = cache.prepare(
        state, {dependencies.data(), dependencies.size()}, transaction,
        certificate);
    if (status) {
      status = transaction.collective_finish(MPI_COMM_SELF, outcome);
    }
    if (transaction.finished()) {
      const Status finalized = cache.finalize(transaction);
      if (status && !finalized) status = finalized;
    }
    return status;
  };

  FinalForceCertificate pending;
  const Status rejected = run_attempt(
      {StatusCode::rejected_step, 77U}, pending);
  SurfaceForce committed_force;
  FinalForceCertificate committed_certificate;
  passed &= expect(rejected.code == StatusCode::rejected_step &&
                       cache.committed(committed_force, committed_certificate)
                               .code == StatusCode::invalid_plan,
                   "failed attempt publishes no final-force cache");
  const Status accepted = run_attempt(Status{}, pending);
  passed &= expect(static_cast<bool>(accepted) &&
                       cache.committed(committed_force,
                                       committed_certificate) &&
                       committed_certificate.valid() &&
                       committed_force.revision == committed_certificate.force,
                   "successful commit publishes exactly one final-force cache");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  bool passed = test_revision_and_flux_mutations();
  passed &= test_sign_area_normal_mutations();
  passed &= test_force_cache_rollback_and_commit();
  MPI_Finalize();
  return passed ? 0 : 1;
}
