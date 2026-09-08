// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/candidate_boundary_fixture.hpp"
#include <mpi.h>
#include <cmath>
#include <iostream>

using namespace hundun::v04;
using namespace hundun::v04::test;

int run() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.inlet = CandidateBoundaryInlet::mass_flow;
  spec.mass_flow_rate = 0.375;
  spec.outlet_kind = BoundaryKind::zero_gradient_mass_outlet;
  spec.predictor_velocity = 0.25;
  const bool initialized = fixture.initialize(MPI_COMM_WORLD, spec);
  int okay = initialized ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!okay) {
    std::cerr << "FAIL: mass-balanced outlet real-chain initialization rank="
              << rank << " stage=" << fixture.diagnostic_step << '\n';
    return 1;
  }
  CandidateBoundaryScratch candidate;
  const bool staged = fixture.stage(0.5, 8.0, 0.0, 71000U, candidate);
  okay = staged ? 1 : 0;
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!okay)
    std::cerr << "FAIL: candidate staging rank=" << rank
              << " stage=" << fixture.diagnostic_step
              << " status=" << unsigned(fixture.diagnostic_status.code)
              << '/' << fixture.diagnostic_status.detail << '\n';
  if (okay) {
    double q = 0.0, storage = 0.0;
    const auto cells = fixture.patch.cells;
    for (int z = 0; z < cells.z; ++z)
      for (int y = 0; y < cells.y; ++y) {
        if (fixture.local_face_owner(CartesianFace::x_max))
          q += candidate.final_flux.x.unchecked({cells.x, y, z});
        for (int x = 0; x < cells.x; ++x) {
          const Int3 c{x,y,z};
          storage += fixture.cell_volume(c) *
              (fixture.bdf.a0 * candidate.density.view.unchecked(c,0U) +
               fixture.bdf.a1 * fixture.accepted_density.view.unchecked(c,0U));
        }
      }
    const double local[2]{q, storage};
    double global[2]{};
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    okay = std::abs(global[0] + global[1] - spec.mass_flow_rate) < 1e-13 &&
        std::abs(global[1]) > 1e-8;
    if (fixture.local_face_owner(CartesianFace::x_max)) {
      okay &= candidate.pressure.view.unchecked({cells.x,0,0},0U) ==
          candidate.pressure.view.unchecked({cells.x-1,0,0},0U);
      okay &= candidate.scaled_pressure.view.unchecked({cells.x,0,0},0U) == 4.0;
    }
    if (!okay) std::cerr << "FAIL: Neumann/BDF storage closure out=" << global[0]
                        << " storage=" << global[1] << '\n';
  }
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (okay) {
    const auto cells = fixture.patch.cells;
    // The fresh physical closure shares the production finalizer engine.
    // Supply a nonuniform outlet with one reversed face; preserve that shape.
    if (fixture.local_face_owner(CartesianFace::x_max))
      for (int z = 0; z < cells.z; ++z)
        for (int y = 0; y < cells.y; ++y)
          candidate.mechanical_flux.x.unchecked({cells.x, y, z}) =
              y == 0 && z == 0 ? -0.01 : 0.02;
    const FreshPhysicalBoundaryFluxClosureInput input{
        spec.outlet_pressure, candidate.thermophysical_pressure_alias,
        as_const(candidate.velocity.view),
        {candidate.semantic_independent_species_views.data(),
         candidate.semantic_independent_species_views.size()},
        as_const(candidate.mechanical_flux), candidate.final_flux};
    const Status status = fixture.finalizer.close_fresh_physical_flux(
        input, fixture.reductions);
    double local[3]{}, global[3]{};
    if (fixture.local_face_owner(CartesianFace::x_max)) {
      for (int z = 0; z < cells.z; ++z)
        for (int y = 0; y < cells.y; ++y) {
          local[0] += candidate.final_flux.x.unchecked({cells.x, y, z});
          local[2] += candidate.mechanical_flux.x.unchecked({cells.x, y, z});
        }
      local[1] = candidate.final_flux.x.unchecked({cells.x, 0, 0});
    }
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    okay = status && std::abs(global[0] - spec.mass_flow_rate) < 1e-13 &&
           global[1] < 0.0;
    if (!okay)
      std::cerr << "FAIL: signed outlet shape/global mass closure status="
                << unsigned(status.code) << '/' << status.detail
                << " out=" << global[0] << " reversed=" << global[1] << '\n';
    if (fixture.local_face_owner(CartesianFace::x_max))
      for (int z = 0; z < cells.z; ++z)
        for (int y = 0; y < cells.y; ++y) {
          const double expected = (y == 0 && z == 0 ? -0.01 : 0.02) *
              spec.mass_flow_rate / global[2];
          okay &= std::abs(candidate.final_flux.x.unchecked({cells.x,y,z}) - expected) < 1e-13;
        }
    if (!okay) std::cerr << "FAIL: uniform signed shape rank=" << rank << '\n';
    // A zero or negative net outlet pool cannot be repaired by inventing an
    // area-weighted flow. Failure must leave all published face bytes intact.
    for (const double bad : {0.0, -0.02}) {
      std::vector<double> before;
      for (const auto face : {candidate.final_flux.x, candidate.final_flux.y, candidate.final_flux.z})
        for (int z = 0; z < face.extents.z; ++z)
          for (int y = 0; y < face.extents.y; ++y)
            for (int x = 0; x < face.extents.x; ++x)
              before.push_back(face.unchecked({x,y,z}));
      if (fixture.local_face_owner(CartesianFace::x_max))
        for (int z = 0; z < cells.z; ++z)
          for (int y = 0; y < cells.y; ++y)
            candidate.mechanical_flux.x.unchecked({cells.x,y,z}) = bad;
      const Status rejected = fixture.finalizer.close_fresh_physical_flux(input, fixture.reductions);
      okay &= rejected.code == StatusCode::rejected_step;
      std::size_t i = 0;
      for (const auto face : {candidate.final_flux.x, candidate.final_flux.y, candidate.final_flux.z})
        for (int z = 0; z < face.extents.z; ++z)
          for (int y = 0; y < face.extents.y; ++y)
            for (int x = 0; x < face.extents.x; ++x) {
              const double after = face.unchecked({x,y,z});
              if (before[i] != after) {
                std::cerr << "FAIL: atomic face rank=" << rank << " index=" << i
                    << " before=" << before[i] << " after=" << after << '\n';
                okay = 0;
              }
              ++i;
            }
      if (!okay) std::cerr << "FAIL: shape/atomic rejection rank=" << rank
          << " pool=" << bad << " status=" << unsigned(rejected.code)
          << '/' << rejected.detail << '\n';
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &okay, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0 && okay) std::cout << "PASS: mass-balanced outlet closure\n";
  return okay ? 0 : 1;
}

int patch_test() {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.inlet = CandidateBoundaryInlet::mass_flow;
  spec.mass_flow_rate = 0.375;
  if (!fixture.initialize(MPI_COMM_WORLD, spec)) return 1;
  std::array<std::vector<Int3>, 2U> support;
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z = 0; z < fixture.patch.cells.z; ++z)
      for (int y = 0; y < fixture.patch.cells.y; ++y)
        support[y + fixture.patch.begin.y < fixture.geometry.global_cells().y / 2
            ? 0 : 1].push_back({0,y,z});
  const std::array<PhysicalMassFlowPatch, 2U> patches{{
      {CartesianFace::x_min, {1,0,0}, 0.125, {support[0].data(),support[0].size()}},
      {CartesianFace::x_min, {1,0,0}, 0.250, {support[1].data(),support[1].size()}}}};
  auto binding = fixture.finalizer_binding();
  binding.mass_flow_patches = {patches.data(), patches.size()};
  if (!PressureEnergyCandidateBoundaryFinalizer::bind(binding, fixture.finalizer)) return 1;
  const auto fingerprint = fixture.finalizer.fingerprint();
  auto duplicate = patches;
  duplicate[1].local_faces = duplicate[0].local_faces;
  binding.mass_flow_patches = {duplicate.data(), duplicate.size()};
  const Status rejected = PressureEnergyCandidateBoundaryFinalizer::bind(binding, fixture.finalizer);
  if (rejected.code != StatusCode::invalid_plan ||
      fixture.finalizer.fingerprint() != fingerprint) return 1;
  CandidateBoundaryScratch candidate;
  if (!fixture.stage(0.0, 8.0, 0.0, 72000U, candidate)) return 1;
  double local[2]{}, global[2]{};
  for (std::size_t p = 0; p < 2; ++p)
    for (const Int3 face : support[p])
      local[p] += candidate.final_flux.x.unchecked(face);
  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const bool okay = std::abs(global[0]-0.125)<1e-13 && std::abs(global[1]-0.250)<1e-13;
  if (!okay) std::cerr << "FAIL: independent patch targets " << global[0] << ',' << global[1] << '\n';
  return okay ? 0 : 1;
}

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 2;
  const int result = run() | patch_test();
  MPI_Finalize();
  return result;
}
