// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "../support/candidate_boundary_fixture.hpp"
#include "../../src/solver_scalar_mass_remap_detail.hpp"
#include "../support/dense_solve.hpp"

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

bool frozen_transport_pairing() {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.multispecies=true; spec.cells_per_axis=9;
  if (!fixture.initialize(MPI_COMM_WORLD,spec)) return false;
  const auto cells=fixture.patch.cells;
  const int nx=fixture.geometry.global_cells().x;
  constexpr double a0=2.3;
  const double pi=std::acos(-1.);
  const auto face=[&](int x,double amplitude) {
    return x==0 || x==nx ? 0. : amplitude*std::sin(2*pi*x/nx);
  };
  const auto star_mass=[&](int x) {return 1-(face(x+1,.12)-face(x,.12))/a0;};
  const auto final_mass=[&](int x) {return 1-(face(x+1,.07)-face(x,.07))/a0;};
  const auto value=[&](unsigned field,unsigned species,int x) {
    const double angle=2*pi*(x+.5)/nx;
    const double signal=std::sin(angle)+.4*std::cos(angle);
    return species==0 ? .25+(.02+.01*field)*signal : .2-(.01+.005*field)*signal;
  };
  std::array<FieldId,2> ids;
  std::array<FieldView,2> scalars;
  const std::array<TransportedScalarRole,2> roles{
      TransportedScalarRole::species,TransportedScalarRole::species};
  for(unsigned s=0;s<2;++s) {scalars[s]=fixture.independent_species[s].view;ids[s]=scalars[s].field;}
  detail::ScalarMassRemap remap;
  auto status=remap.allocate(fixture.patch,{ids.data(),2},{roles.data(),2},2);
  if(status)status=remap.bind(MPI_COMM_WORLD,fixture.boundary);
  FaceFluxStorage storage; FaceFluxView base,final;
  if(status)status=FaceFluxStorage::allocate_workspace(cells,2,storage);
  if(status)status=storage.workspace_view(0,81,base);
  if(status)status=storage.workspace_view(1,82,final);
  if(!status)return false;
  for(auto flux:{base,final})for(auto f:{flux.x,flux.y,flux.z})
    for(int z=0;z<f.extents.z;++z)for(int y=0;y<f.extents.y;++y)for(int x=0;x<f.extents.x;++x)
      f.unchecked({x,y,z})=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<=cells.x;++x) {
    base.x.unchecked({x,y,z})=face(x+fixture.patch.begin.x,.12);
    final.x.unchecked({x,y,z})=face(x+fixture.patch.begin.x,.07);
  }
  const BoundaryResolvedValues boundary_values{
      {fixture.boundary_scalar_values.data(),fixture.boundary_scalar_values.size()},
      {fixture.boundary_vector_values.data(),fixture.boundary_vector_values.size()},
      {fixture.boundary_normal_gradient_values.data(),fixture.boundary_normal_gradient_values.size()}};
  bool passed=true;
  double local[4]{};
  for(unsigned field=0;field<4;++field) {
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
      const Int3 c{x,y,z};
      fixture.density.view.unchecked(c,0)=1/fixture.cell_volume(c);
      for(unsigned s=0;s<2;++s)scalars[s].unchecked(c,0)=value(field,s,x+fixture.patch.begin.x);
    }
    // Each stochastic realization keeps its own conservative transported
    // inventory while sharing the same pressure-driven mass correction.
    status=remap.capture_frozen_density(fixture.kernels,as_const(fixture.density.view),
        {scalars.data(),2},as_const(base),a0);
    if(!status)return false;
    const auto bytes=remap.owned_payload_bytes();
    const Int3 last{cells.x-1,cells.y-1,cells.z-1};
    const double saved=scalars[1].unchecked(last,0);
    scalars[1].unchecked(last,0)=std::numeric_limits<double>::quiet_NaN();
    const auto invalid=remap.capture_frozen_density(fixture.kernels,as_const(fixture.density.view),
        {scalars.data(),2},as_const(base),a0*2);
    scalars[1].unchecked(last,0)=saved;
    passed &= invalid.code==StatusCode::rejected_step && remap.owned_payload_bytes()==bytes;
    const double rho_before=fixture.density.view.unchecked(last,0);
    fixture.density.view.unchecked(last,0)=-1;
    const auto negative=remap.capture_frozen_density(fixture.kernels,as_const(fixture.density.view),
        {scalars.data(),2},as_const(base),a0*2);
    fixture.density.view.unchecked(last,0)=rho_before;
    passed &= negative.code==StatusCode::rejected_step;
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
      fixture.density.view.unchecked({x,y,z},0)=final_mass(x+fixture.patch.begin.x)/fixture.cell_volume({x,y,z});
    detail::ScalarMassRemap::Report report;
    status=remap.solve(fixture.kernels,as_const(fixture.density.view),as_const(fixture.velocity.view),
        as_const(final),{scalars.data(),2},{},boundary_values,fixture.reductions,report);
    if(!status)return false;
    remap.copy_solution({scalars.data(),2},TransportedScalarRole::species);
    local[1]=std::max(local[1],report.mass_pairing_residual);
    for(unsigned s=0;s<2;++s) {
      std::vector<std::vector<double>> matrix(nx,std::vector<double>(nx));
      std::vector<double> rhs(nx);
      for(int x=0;x<nx;++x) {
        const double west=-(face(x,.07)-face(x,.12))/a0;
        const double east=(face(x+1,.07)-face(x+1,.12))/a0;
        matrix[x][x]=final_mass(x)+std::max(west,0.)+std::max(east,0.);
        if(x>0)matrix[x][x-1]=std::min(west,0.);
        if(x+1<nx)matrix[x][x+1]=std::min(east,0.);
        rhs[x]=star_mass(x)*value(field,s,x);
      }
      const auto expected=dense_solve(matrix,rhs);
      long double inventories[2]{};
      for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
        const int gx=x+fixture.patch.begin.x;
        const double q=scalars[s].unchecked({x,y,z},0);
        local[0]=std::max(local[0],std::abs(q-expected[gx]));
        inventories[0]+=star_mass(gx)*value(field,s,gx);
        inventories[1]+=final_mass(gx)*q;
        passed &= q>=.1 && q<=.4;
      }
      MPI_Allreduce(MPI_IN_PLACE,inventories,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      local[2]=std::max(local[2],double(std::abs(inventories[1]-inventories[0])/inventories[0]));
    }
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
      passed &= scalars[0].unchecked({x,y,z},0)+scalars[1].unchecked({x,y,z},0)<=1;
  }
  double global[4]{};
  MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(rank==0)std::cout<<"frozen_transport dense_error="<<global[0]<<" mass_pairing="<<global[1]<<" species_inventory="<<global[2]<<'\n';
  return passed && global[0]<1e-12 && global[1]<1e-14 && global[2]<1e-12;
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
  int local=(run() && coupling_forcing_replay() && history_contract() && frozen_transport_pairing()) ? 1 : 0, global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}
