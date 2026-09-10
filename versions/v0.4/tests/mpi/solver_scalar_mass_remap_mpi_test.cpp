// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/candidate_boundary_fixture.hpp"
#include "../../src/solver_scalar_mass_remap_detail.hpp"

#include <mpi.h>
#include <array>
#include <iostream>
#include <limits>

using namespace hundun::v04;
using namespace hundun::v04::test;

bool coupling_forcing_replay() {
  // Captured GTMC trajectory: the first two approximate flow solutions
  // already resolve the composition change, while the final sweep must
  // evaluate species on the fully converged, certified flux.
  detail::SpeciesCouplingForcing forcing;
  constexpr double c = 16.0 * std::numeric_limits<double>::epsilon();
  constexpr double e = 1.0e-10;
  forcing.begin(1U);
  bool passed = forcing.can_update_composition(1.8484e-9, 1.5598e-7, c, e);
  forcing.observe(1.6819e-5);
  forcing.begin(2U);
  passed &= forcing.can_update_composition(1.8489e-9, 1.5596e-7, c, e);
  forcing.observe(9.4519e-7);
  forcing.begin(3U);
  passed &= !forcing.can_update_composition(1.8489e-9, 1.5596e-7, c, e);
  passed &= forcing.can_update_composition(7.9334e-14, 6.1633e-11, c, e);
  // A stalled outer update requests a complete inner solve even if its
  // current C/E sample would otherwise qualify for an early guess.
  forcing.observe(9.0e-7);
  forcing.begin(4U);
  passed &= !forcing.can_update_composition(0.0, 0.0, c, e);
  forcing.observe(2.0e-12);
  forcing.begin(5U);
  passed &= !forcing.can_update_composition(0.0, 0.0, c, e);
  // A new proposal cannot inherit the old proposal's forcing or failure.
  forcing.begin(1U);
  passed &= forcing.can_update_composition(1.8484e-9, 1.5598e-7, c, e);
  forcing.observe(std::numeric_limits<double>::quiet_NaN());
  forcing.begin(2U);
  passed &= !forcing.can_update_composition(0.0, 0.0, c, e);
  forcing.begin(4U);
  passed &= forcing.can_update_composition(1.8484e-9, 1.5598e-7, c, e);
  if (!passed) std::cerr << "FAIL: composition forcing trajectory and lifecycle\n";
  return passed;
}

bool history_contract() {
  detail::SpeciesCouplingHistory history;
  ReductionEngine reductions;
  Status status = ReductionEngine::compile(MPI_COMM_WORLD,
      ReductionMode::mpi_allreduce,
      detail::SpeciesCouplingHistory::reduction_capacity, reductions);
  if (status) status = history.reserve(2U);
  if (!status) return false;
  const auto bytes = history.owned_payload_bytes();
  std::array<double,2U> input{0.1,0.1}, solved{};
  std::array<FieldView,2U> in{}, out{};
  const std::array<std::size_t,2U> species{0U,1U};
  for (std::size_t s=0; s<2U; ++s) {
    in[s].base=&input[s];out[s].base=&solved[s];
    in[s].interior=out[s].interior={1,1,1};
  }
  const auto update = [&](unsigned sweep, double residual, bool& used) {
    history.begin(sweep);
    return history.update({1,1,1},{in.data(),2U},{out.data(),2U},
        {species.data(),2U},{},residual,reductions,used);
  };
  bool passed=true, accelerated=false;
  for (unsigned sweep=1U; sweep<=4U; ++sweep) {
    solved={0.5*input[0]+0.125, 0.25*input[1]+0.15};
    const auto actual=solved;
    const double residual=std::max(std::abs(solved[0]-input[0]),
                                   std::abs(solved[1]-input[1]));
    bool used=false;
    status=update(sweep,residual,used);
    passed &= status && history.valid() && solved==actual;
    accelerated |= used;
    if (!status) return false;
    input={history.guess(0U),history.guess(1U)};
  }
  passed &= accelerated && std::abs(input[0]-0.25)<1e-12 &&
      std::abs(input[1]-0.2)<1e-12 && history.owned_payload_bytes()==bytes;
  history.begin(1U);
  passed &= !history.valid();
  // A negative extrapolated fraction falls back to the entire solved cell.
  input={0.9,0.0};solved={0.1,0.0};bool used=false;
  status=update(1U,0.8,used);
  input={0.1,0.0};solved={0.01,0.0};
  if (status) status=update(2U,0.09,used);
  passed &= status && used && history.guess(0U)==solved[0] &&
      history.guess(1U)==solved[1];
  // A stalled map and a gap both discard the previous fitting directions.
  if (status) status=update(3U,0.1,used);
  passed &= status && !used;
  history.begin(5U);
  passed &= !history.valid() && history.owned_payload_bytes()==bytes;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  status=reductions.arm_checked_sum_fault_for_test(1U,size-1);
  if (status) status=update(6U,0.01,used);
  passed &= !status && !history.valid() &&
      reductions.lowest_failing_rank()==size-1 && solved[0]==0.01;
  reductions.clear_checked_sum_fault_for_test();
#endif
  // Unequal direction magnitudes retain the same least-squares solution.
  std::array<double,21U> gram{};
  gram[0]=1e-20;gram[5]=1.0;gram[16]=2e-21;gram[17]=0.3;
  gram[20]=0.09+4e-22;
  std::array<double,4U> gamma{};
  passed &= detail::SpeciesCouplingHistory::coefficients(gram,2U,gamma) &&
      std::abs(gamma[0]-0.2)<1e-12 && std::abs(gamma[1]-0.3)<1e-12;
  gram[5]=0.0;
  passed &= !detail::SpeciesCouplingHistory::coefficients(gram,2U,gamma);
  if (!passed) std::cerr << "FAIL: species history map, bounds, lifecycle or failure\n";
  return passed;
}

bool run() {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.multispecies = true;
  if (!fixture.initialize(MPI_COMM_WORLD, spec)) return false;
  std::array<FieldId, 2> fields{};
  std::array<FieldView, 2> scalars{};
  const std::array<TransportedScalarRole, 2> roles{
      TransportedScalarRole::species, TransportedScalarRole::species};
  for (std::size_t s=0; s<2; ++s) {
    fields[s] = fixture.independent_species[s].view.field;
    scalars[s] = fixture.independent_species[s].view;
    fill(fixture.independent_species[s], 0.0);
  }
  const auto cells = fixture.patch.cells;
  for (int z=0; z<cells.z; ++z)
    for (int y=0; y<cells.y; ++y)
      for (int x=0; x<cells.x; ++x)
        fixture.density.view.unchecked({x,y,z},0U) =
            1.0 / fixture.cell_volume({x,y,z});
  detail::ScalarMassRemap remap;
  Status status = remap.allocate(fixture.patch, {fields.data(),2U},
                                 {roles.data(),2U}, 2U);
  if (status) status = remap.bind(MPI_COMM_WORLD, fixture.boundary);
  FaceFluxStorage storage;
  FaceFluxView base, final;
  if (status) status = FaceFluxStorage::allocate_workspace(cells,2U,storage);
  if (status) status = storage.workspace_view(0U,81U,base);
  if (status) status = storage.workspace_view(1U,82U,final);
  if (!status) return false;
  for (auto flux : {base,final})
    for (auto face : {flux.x,flux.y,flux.z})
      for (int z=0; z<face.extents.z; ++z)
        for (int y=0; y<face.extents.y; ++y)
          for (int x=0; x<face.extents.x; ++x) face.unchecked({x,y,z})=0.0;
  status = remap.capture(fixture.kernels, as_const(fixture.density.view),
                         {scalars.data(),2U}, as_const(base), 1.0);
  // Fixed physical donor is 0.2/0.3. The extended-precision incoming
  // quantity is nonzero, but its correctly rounded binary64 q is zero.
  // A relative residual of one must not demand an unrepresentable update.
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y)
        final.x.unchecked({0,y,z}) = std::numeric_limits<double>::denorm_min();
  const BoundaryResolvedValues values{
      {fixture.boundary_scalar_values.data(),fixture.boundary_scalar_values.size()},
      {fixture.boundary_vector_values.data(),fixture.boundary_vector_values.size()},
      {fixture.boundary_normal_gradient_values.data(),fixture.boundary_normal_gradient_values.size()}};
  detail::ScalarMassRemap::Report report;
  if (status) status = remap.solve(fixture.kernels, as_const(fixture.density.view),
      as_const(fixture.velocity.view), as_const(final), {scalars.data(),2U}, {},
      values, fixture.reductions, report);
  if (!status) std::cerr << "FAIL: representable trace remap "
      << unsigned(status.code) << '/' << status.detail << " residual="
      << report.residual << " iterations=" << report.iterations << '\n';
  bool passed = status && report.residual == 1.0 &&
      report.convergence_residual == 0.0 && report.iterations == 0U;
  remap.copy_solution({scalars.data(),2U}, TransportedScalarRole::species);
  for (auto scalar : scalars)
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y)
        for (int x=0; x<cells.x; ++x)
          passed &= scalar.unchecked({x,y,z},0U) == 0.0;
  // Ordinary representable transport still requires a real update; the
  // quantization gate must not mistake an initially zero field for closure.
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y) final.x.unchecked({0,y,z}) = 1.0e-3;
  status = remap.solve(fixture.kernels, as_const(fixture.density.view),
      as_const(fixture.velocity.view), as_const(final), {scalars.data(),2U}, {},
      values, fixture.reductions, report);
  passed &= status && report.iterations > 0U &&
      report.residual <= detail::ScalarMassRemap::tolerance;
  remap.copy_solution({scalars.data(),2U}, TransportedScalarRole::species);
  if (fixture.local_face_owner(CartesianFace::x_min))
    for (int z=0; z<cells.z; ++z)
      for (int y=0; y<cells.y; ++y) {
        passed &= std::abs(scalars[0].unchecked({0,y,z},0U)-2.0e-4) < 1e-18;
        passed &= std::abs(scalars[1].unchecked({0,y,z},0U)-3.0e-4) < 1e-18;
      }
  if (!passed) std::cerr << "FAIL: quantization and ordinary update contracts\n";
  return passed;
}

int main(int argc, char** argv) {
  if (MPI_Init(&argc,&argv) != MPI_SUCCESS) return 2;
  int local=(run() && coupling_forcing_replay() && history_contract()) ? 1 : 0, global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
