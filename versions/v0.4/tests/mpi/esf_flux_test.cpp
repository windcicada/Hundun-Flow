// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#include "../support/candidate_boundary_fixture.hpp"
#include "../support/dense_solve.hpp"
#include "../../src/solver_esf_flux_detail.hpp"
#include <iostream>
#include <limits>

using namespace hundun::v04;
using namespace hundun::v04::test;

bool run(int mode, bool immersed=false) {
  const bool periodic=std::abs(mode)==2;
  const int direction=mode==0 ? 0 : mode>0 ? 1 : -1;
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;
  spec.multispecies=true; spec.cells_per_axis=9;
  PeriodicPisoFixture periodic_fixture;
  if (periodic ? !periodic_fixture.initialize(9,MPI_COMM_WORLD,true) :
                 !fixture.initialize(MPI_COMM_WORLD,spec)) return false;
  const auto patch=periodic ? periodic_fixture.patch : fixture.patch;
  const auto cells=patch.cells;
  const auto& geometry=periodic ? periodic_fixture.geometry : fixture.geometry;
  const auto& kernels=periodic ? periodic_fixture.equations.kernels() : fixture.kernels;
  const auto& boundary=periodic ? periodic_fixture.boundary : fixture.boundary;
  const auto density=periodic ? periodic_fixture.density.view : fixture.density.view;
  const auto velocity=periodic ? periodic_fixture.velocity.view : fixture.velocity.view;
  ReductionEngine reductions;
  if(!ReductionEngine::compile(MPI_COMM_WORLD,ReductionMode::mpi_allreduce,4,reductions))return false;
  const int nx=geometry.global_cells().x;
  int rank{},size{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  const std::size_t count=std::size_t(cells.x)*cells.y*cells.z;
  constexpr double dt=.43;
  const double pi=std::acos(-1.);
  const auto face=[&](int x,double amplitude,double offset) {
    return immersed && (x==nx/2 || x==nx/2+1) ? 0. :
        offset+(x==0 || x==nx ? 0. : amplitude*std::sin(2*pi*x/nx));
  };
  const auto fs=[&](int x){return face(x,.12,.02*direction);};
  const auto ff=[&](int x){return face(x,.07,.06*direction);};
  const auto mass=[&](int x,bool final){return 1-dt*(final ? ff(x+1)-ff(x) : fs(x+1)-fs(x));};
  const auto value=[&](unsigned f,unsigned s,int x,bool shifted) {
    if(immersed && x==nx/2)return s<3 ? 0. : -19.;
    const double phase=2*pi*(x+.5)/nx;
    const double a=f==5 ? -.0003+.001*std::sin(phase) :
        .25+(.02+.01*(f==4 ? 1.5 : f))*std::sin(phase);
    const double b=.2-.03*std::cos(phase);
    const double q[4]{a,b,1-a-b,150000+30000*std::cos(phase)-300000*a+200000*b};
    return q[s]+(s==3 && shifted ? 2000000*a-1000000*b+500000*(1-a-b) : 0);
  };
  auto seed=make_field(200,cells,4,0,1,2000);
  auto iterate=make_field(201,cells,4,1,1,2001);
  auto next=make_field(202,cells,4,0,1,2002);
  HaloEngine halo;
  const HaloFieldSpec halo_spec{201,1,4};
  auto status=halo.reserve(MPI_COMM_WORLD,patch,{&halo_spec,1},boundary.halo_topology());
  FaceFluxStorage storage; FaceFluxView frozen,final;
  if(status) status=FaceFluxStorage::allocate_workspace(cells,2,storage);
  if(status) status=storage.workspace_view(0,81,frozen);
  if(status) status=storage.workspace_view(1,82,final);
  if(!status)return false;
  for(auto flux:{frozen,final})for(auto f:{flux.x,flux.y,flux.z})
    for(int z=0;z<f.extents.z;++z)for(int y=0;y<f.extents.y;++y)for(int x=0;x<f.extents.x;++x)
      f.unchecked({x,y,z})=0;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<=cells.x;++x) {
    frozen.x.unchecked({x,y,z})=fs(x+patch.begin.x);
    final.x.unchecked({x,y,z})=ff(x+patch.begin.x);
  }
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
    density.unchecked({x,y,z},0)=1/detail::cell_volume(kernels,{x,y,z});
  std::vector<std::uint8_t> active(count,1);
  std::size_t offset{};
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++offset)
    if(immersed && x+patch.begin.x==nx/2) active[offset]=0;
  std::vector<double> solutions(6*4*count);
  double maximum[5]{};
  bool passed=true;
  for(unsigned field=0;field<6;++field)for(bool shifted:{false,true}) {
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
      for(unsigned s=0;s<4;++s)seed.view.unchecked({x,y,z},s)=value(field,s,x+patch.begin.x,shifted);
    const auto immutable=seed.storage;
    // An open boundary has a prescribed thermochemical tuple; h follows its
    // own coordinate and affine reference, while Y includes the dependent species.
    const auto close=[&](FieldView v) {
      if(periodic)return Status{};
      for(unsigned f=0;f<6;++f) {
        const auto face_id=static_cast<CartesianFace>(f);
        if(!fixture.local_face_owner(face_id))continue;
        const int axis=f/2, n=axis==0 ? cells.x : axis==1 ? cells.y : cells.z;
        for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
          Int3 c{x,y,z},ghost=c;
          int& coordinate=axis==0 ? ghost.x : axis==1 ? ghost.y : ghost.z;
          if(coordinate!=(f%2 ? n-1 : 0))continue;
          coordinate=f%2 ? n : -1;
          for(unsigned s=0;s<4;++s) {
            const double q=axis==0 ? value(field,s,f%2 ? nx : -1,shifted) : v.unchecked(c,s);
            v.unchecked(ghost,s)=2*q-v.unchecked(c,s);
          }
        }
      }
      return Status{};
    };
    detail::StatisticalFluxReport report;
    const auto correct=[&](ConstFieldView initial=ConstFieldView{}) {
      return detail::correct_statistical_flux(kernels,as_const(density),
          as_const(seed.view),as_const(frozen),as_const(final),dt,{active.data(),active.size()},field==5,
          boundary,as_const(velocity),iterate.view,next.view,halo,178,
          reductions,close,report,initial,false);
    };
    status=correct();
    if(!status) {
      std::cerr<<"flux status="<<unsigned(status.code)<<'/'<<status.detail<<" residual="<<report.residual<<'\n';
      return false;
    }
    passed &= seed.storage==immutable;
    maximum[1]=std::max(maximum[1],report.mass_pairing_residual);
    maximum[2]=std::max(maximum[2],report.composition_closure);
    for(unsigned s=0;s<4;++s) {
      std::vector<std::vector<double>> matrix(nx,std::vector<double>(nx));
      std::vector<double> rhs(nx);
      for(int x=0;x<nx;++x) {
        const double west=-dt*(ff(x)-fs(x)),east=dt*(ff(x+1)-fs(x+1));
        matrix[x][x]=mass(x,true)+std::max(west,0.)+std::max(east,0.);
        rhs[x]=mass(x,false)*value(field,s,x,shifted);
        if(x>0 || periodic)matrix[x][(x+nx-1)%nx]=std::min(west,0.);
        else rhs[x]-=std::min(west,0.)*value(field,s,-1,shifted);
        if(x+1<nx || periodic)matrix[x][(x+1)%nx]=std::min(east,0.);
        else rhs[x]-=std::min(east,0.)*value(field,s,nx,shifted);
      }
      const auto expected=dense_solve(matrix,rhs);
      long double inventory[2]{};
      std::size_t i{};
      for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++i) {
        const Int3 c{x,y,z};const int gx=x+patch.begin.x;
        const double q=iterate.view.unchecked(c,s);
        maximum[0]=std::max(maximum[0],std::abs(q-expected[gx])/std::max(1.,std::abs(expected[gx])));
        inventory[0]+=mass(gx,false)*value(field,s,gx,shifted);
        inventory[1]+=mass(gx,true)*q;
        const auto slot=(field*4+s)*count+i;
        if(!shifted)solutions[slot]=q;
        else if(s<3)passed &= q==solutions[slot];
        else {
          const double h=solutions[slot]+2000000*solutions[(field*4)*count+i]
              -1000000*solutions[(field*4+1)*count+i]+500000*solutions[(field*4+2)*count+i];
          maximum[3]=std::max(maximum[3],std::abs(q-h)/std::max(1.,std::abs(h)));
        }
        if(field<5 && s<3)passed &= q>=0 && q<=1;
        if(field==4 && !shifted) {
          double mean{};for(unsigned f=0;f<4;++f)mean+=solutions[(f*4+s)*count+i]/4;
          maximum[4]=std::max(maximum[4],std::abs(q-mean)/std::max(1.,std::abs(mean)));
        }
      }
      MPI_Allreduce(MPI_IN_PLACE,inventory,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      if(direction==0 || periodic)passed &= std::abs(inventory[1]-inventory[0])<1e-12L*std::max(1.L,std::abs(inventory[0]));
    }
    const auto solved=iterate.storage;
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
      for(unsigned s=0;s<4;++s)iterate.view.unchecked({x,y,z},s)=-91.;
    status=correct();
    passed &= status && iterate.storage==solved;
    // Audit evaluates the submitted tuple without solving away its defect.
    auto candidate=make_field(203,cells,4,0,1,2003);
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
      for(unsigned s=0;s<4;++s)candidate.view.unchecked({x,y,z},s)=iterate.view.unchecked({x,y,z},s);
    const auto audit=[&](ConstFieldView submitted) {
      return detail::correct_statistical_flux(kernels,as_const(density),
          as_const(seed.view),as_const(frozen),as_const(final),dt,{active.data(),active.size()},field==5,
          boundary,as_const(velocity),iterate.view,next.view,halo,178,
          reductions,close,report,submitted,true);
    };
    const auto submitted=candidate.storage;
    // A converged initial guess keeps the frozen physical seed and equation.
    status=correct(as_const(candidate.view));
    passed &= status && report.iterations==0 && iterate.storage==solved &&
        candidate.storage==submitted && seed.storage==immutable;
    status=audit(as_const(candidate.view));
    passed &= status && report.convergence_residual<1e-12 &&
        candidate.storage==submitted && seed.storage==immutable;
    std::size_t ci{};
    for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x,++ci)
      if(active[ci]) {
        candidate.view.unchecked({x,y,z},0)+=.00001;
        candidate.view.unchecked({x,y,z},1)-=.00001;
      }
    const auto perturbed=candidate.storage;
    status=audit(as_const(candidate.view));
    passed &= status && report.convergence_residual>1e-8 &&
        report.composition_closure<2e-12 && candidate.storage==perturbed && seed.storage==immutable;
    passed &= !audit(as_const(iterate.view)) && !audit({});
    // A displaced initial guess converges to the same dense-reference solution.
    status=correct(as_const(candidate.view));
    passed &= status && candidate.storage==perturbed && seed.storage==immutable;
    for(std::size_t i=0;i<solved.size();++i)
      passed &= std::abs(iterate.storage[i]-solved[i])
          <=6e-14*std::max(1.,std::abs(solved[i]));
    status=correct();
    passed &= status && iterate.storage==solved;
    // Late validation failures preserve the previous candidate on all ranks.
    const Int3 last{cells.x-1,cells.y-1,cells.z-1};
    const double q=seed.view.unchecked(last,3);
    if(rank==size-1)seed.view.unchecked(last,3)=std::numeric_limits<double>::quiet_NaN();
    status=correct();
    passed &= !status && iterate.storage==solved;
    seed.view.unchecked(last,3)=q;
    if(immersed) {
      const int sx=nx/2-patch.begin.x;
      if(sx>=0 && sx<cells.x)final.x.unchecked({sx,0,0})=.1;
      status=correct();
      passed &= !status && iterate.storage==solved;
      if(sx>=0 && sx<cells.x)final.x.unchecked({sx,0,0})=0;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,maximum,5,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  if(rank==0)std::cout<<"tuple_flux direction="<<direction<<" periodic="<<periodic<<" solid="<<immersed<<" dense="<<maximum[0]
      <<" mass_pair="<<maximum[1]<<" sumY="<<maximum[2]<<" href="<<maximum[3]<<" mean="<<maximum[4]<<'\n';
  return passed && maximum[0]<1e-11 && maximum[1]<1e-14 && maximum[2]<2e-12 &&
      maximum[3]<1e-11 && maximum[4]<1e-11;
}
bool native_storage() {
  CandidateBoundaryFixture fixture;
  CandidateBoundaryFixtureSpec spec;spec.cells_per_axis=9;
  if(!fixture.initialize(MPI_COMM_WORLD,spec))return false;
  FieldRegistry registry;FieldSchema schema;
  std::array<FieldId,3> ids{};
  for(unsigned i=0;i<3;++i)
    if(!registry.declare_field("tuple"+std::to_string(i),4,1,ids[i]))return false;
  if(!registry.freeze(schema))return false;
  const auto cells=fixture.patch.cells;
  std::array<ArenaFieldRequest,3> requests;
  for(unsigned i=0;i<3;++i)requests[i]={ids[i],cells,{0},FieldLifetime::state_layer};
  ArenaLayout layout;StateLayers layers;
  if(!ArenaLayout::compile(schema,{requests.data(),3},layout) || !StateLayers::allocate(layout,layers))return false;
  FieldView seed,iterate,next;
  if(!layers.view(StateRole::accepted_n,ids[0],seed) ||
     !layers.view(StateRole::trial,ids[1],iterate) ||
     !layers.view(StateRole::trial,ids[2],next))return false;
  bool passed=seed.storage_identity==iterate.storage_identity && iterate.storage_identity==next.storage_identity;
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
    const Int3 c{x,y,z};seed.unchecked(c,0)=.2;seed.unchecked(c,1)=.3;
    seed.unchecked(c,2)=.5;seed.unchecked(c,3)=12000;
  }
  FaceFluxStorage storage;FaceFluxView flux;
  if(!FaceFluxStorage::allocate_workspace(cells,1,storage) || !storage.workspace_view(0,81,flux))return false;
  for(auto f:{flux.x,flux.y,flux.z})
    for(int z=0;z<f.extents.z;++z)for(int y=0;y<f.extents.y;++y)for(int x=0;x<f.extents.x;++x)f.unchecked({x,y,z})=0;
  HaloEngine halo;const HaloFieldSpec hs{iterate.field,1,4};
  if(!halo.reserve(MPI_COMM_WORLD,fixture.patch,{&hs,1},fixture.boundary.halo_topology()))return false;
  detail::StatisticalFluxReport report;
  const auto solve=[&](FieldView scratch) {
    return detail::correct_statistical_flux(fixture.kernels,as_const(fixture.density.view),as_const(seed),
        as_const(flux),as_const(flux),.01,{},false,fixture.boundary,as_const(fixture.velocity.view),
        iterate,scratch,halo,178,fixture.reductions,[](FieldView){return Status{};},report);
  };
  const auto status=solve(next);
  if(!status)std::cerr<<"native_tuple shared_arena status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';
  passed &= bool(status);
  if(status)for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x)
    for(unsigned c=0;c<4;++c)passed &= iterate.unchecked({x,y,z},c)==seed.unchecked({x,y,z},c);
  auto alias=iterate;alias.field=next.field;
  passed &= !solve(alias);
  return passed;
}
int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  bool passed=true;
  for(int direction:{0,1,-1,2,-2})passed=run(direction) && passed;
  passed=run(0,true) && passed;passed=run(2,true) && passed;
  passed=native_storage() && passed;
  int local=passed,global{};
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Finalize();return global ? 0 : 1;
}
