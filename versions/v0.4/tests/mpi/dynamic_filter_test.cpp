// SPDX-License-Identifier: Apache-2.0
#include "../../src/core_tcr_dynamic_detail.hpp"
#include <iostream>
using namespace hundun::v04;
namespace {
struct Field {
  std::vector<double> storage;
  FieldView view;
  Field(Int3 n,unsigned components):storage(std::size_t(n.x)*n.y*n.z*components) {
    view.base=storage.data();view.interior=n;view.components=components;
    view.stride_y=n.x;view.stride_z=n.x*n.y;view.component_stride=std::size_t(n.x)*n.y*n.z;
    view.field=1;view.revision=1;view.storage_identity=reinterpret_cast<StorageIdentity>(storage.data());
    view.revision_domain=17;
  }
};
bool compute(MPI_Comm comm,const CartesianGeometryPlan &geometry,MeshPatch patch,
    int mode,double scale,std::vector<double> &out,bool dyn711) {
  const auto n=patch.cells;
  Field a(n,4),b(n,4),aux(n,4),rho(n,1);
  const auto count=std::size_t(n.x)*n.y*n.z;
  std::vector<std::uint8_t> active(count,1);
  std::size_t i{};
  for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x,++i) {
    const int gx=x+patch.begin.x;
    const double u=.02*std::sin(1.7*gx+2.1*y+.9*z),v=.015*std::cos(.8*gx-1.3*y+2.3*z);
    for(auto f:{a.view,b.view,aux.view}) {
      const double factor=f.base==aux.view.base ? 6 : 1;
      const double q=.2+factor*u,r=.3+factor*v;
      f.unchecked({x,y,z},0)=q;f.unchecked({x,y,z},1)=r;
      f.unchecked({x,y,z},2)=1-q-r;f.unchecked({x,y,z},3)=0;
    }
    rho.view.unchecked({x,y,z},0)=scale*(mode>=2 ? .5+.2*((gx+2*y+z)%5) : 1.);
    if(mode==3 && gx==3 && y==1)active[i]=0;
  }
  detail::DynamicTcrHistory history;
  history.configure(73,count,3,2);
  detail::DynamicTcrPlan plan;
  const std::array<bool,3> periodic{mode!=0,true,true};
  auto status=plan.configure(comm,geometry,patch,periodic,3,{0,1,2},1u<<28);
  if(dyn711) {
    detail::Dyn711History windows;
    windows.configure(73,count,3,2);
    const auto component=[](FieldView view,unsigned q) {
      view.base+=q*view.component_stride;view.components=1;return as_const(view);
    };
    const std::array<ConstFieldView,3> means{component(a.view,0),component(a.view,1),component(a.view,2)};
    const std::array<ConstFieldView,3> gradients{component(aux.view,0),component(aux.view,1),component(aux.view,2)};
    std::vector<double> updated(count);
    for(unsigned step=0;step<3 && status;++step) {
      std::vector<std::uint8_t> proposal;
      // Candidate recomputation starts from exactly the same accepted rates,
      // coefficient and clocks; a rejected attempt publishes no new history.
      for(unsigned attempt=0;attempt<2 && status;++attempt) {
        status=windows.begin(step);
        for(std::size_t c=0;c<count && status;++c) {
          if(!active[c]) {status=windows.stage_inactive(c);continue;}
          for(unsigned q=0;q<3 && status;++q)
            status=windows.stage_rate(c,q,.001,1,2,.3,1,1,1e-30);
        }
        if(status)status=plan.finish(windows,means,gradients,as_const(rho.view),{active.data(),active.size()});
        if(!status)break;
        const auto prepared=windows.prepared_snapshot().values;
        if(attempt==0) {
          proposal.assign(prepared.data,prepared.data+prepared.size);
          windows.discard();
        } else if(!std::equal(proposal.begin(),proposal.end(),prepared.data))return false;
      }
      if(!status)break;
      windows.commit();
      if(windows.statistics_calls()!=step+1)return false;
      for(std::size_t c=0;c<count;++c) {
        if(step==1)updated[c]=windows.cphi(c);
        if(step==2 && windows.cphi(c)!=updated[c])return false;
        if(!active[c] && windows.cphi(c)!=2)return false;
      }
    }
    if(!status){std::cerr<<"dyn711 filter status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
    out.resize(3*count);
    for(i=0;i<count;++i)for(unsigned g=0;g<3;++g)out[3*i+g]=windows.cphi(i);
    int rank{},ranks{};MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    if(ranks>1 && mode==0 && scale==1) {
      if(rank==0)for(unsigned step=3;step<8;++step) {
        if(!windows.begin(step))return false;
        for(std::size_t c=0;c<count;++c)if(!windows.stage_inactive(c))return false;
        if(!windows.seal())return false;
        windows.commit();
      }
      // A mismatched history clock must reject collectively before one rank
      // enters halo exchange while another rank takes the held-Cphi path.
      const auto rejected=plan.finish(windows,means,gradients,as_const(rho.view),{active.data(),active.size()});
      if(rejected.code!=StatusCode::invalid_plan)return false;
    }
  } else {
  if(status)status=history.begin(0);
  const std::array<FieldView,2> fields{a.view,b.view};
  if(status)status=plan.finish(history,{fields.data(),fields.size()},as_const(aux.view),
      as_const(rho.view),{active.data(),active.size()},0);
  if(!status){std::cerr<<"dynamic filter status="<<unsigned(status.code)<<'/'<<status.detail<<'\n';return false;}
  history.commit();
  out.resize(3*count);
  for(i=0;i<count;++i)for(unsigned g=0;g<3;++g)out[3*i+g]=history.accepted(i)[15+g];
  }
  return true;
}
}
int main(int argc,char **argv) {
  MPI_Init(&argc,&argv);
  int rank{},ranks{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  bool okay=true;
  {
    CartesianMeshSpec mesh;
    mesh.kind=GeometryKind::uniform;mesh.lower={0,0,0};mesh.upper={7,5,5};
    mesh.has_exact_cells=true;mesh.exact_cells={7,5,5};mesh.minimum_spacing={1e-9,1e-9,1e-9};
    mesh.max_growth_ratio=1;mesh.limits={10000,1u<<28};
    CartesianGeometryPlan geometry;MeshPatch global;
    auto status=CartesianGeometryCompiler::compile(MPI_COMM_SELF,mesh,{1024,128},geometry,global);
    if(!status){MPI_Abort(MPI_COMM_WORLD,1);}
    auto local=global;local.begin.x=rank*(7/ranks)+std::min(rank,7%ranks);
    local.cells.x=7/ranks+(rank<7%ranks ? 1:0);
    local.process_grid={ranks,1,1};local.process_coord={rank,0,0};
    for(bool dyn711:{false,true}) {
    double invariance{},partition{};unsigned interior_coefficients{};
    for(int mode=0;mode<4;++mode) {
      std::vector<double> serial,distributed,scaled;
      okay=compute(MPI_COMM_SELF,geometry,global,mode,1,serial,dyn711)&&okay;
      okay=compute(MPI_COMM_WORLD,geometry,local,mode,1,distributed,dyn711)&&okay;
      okay=compute(MPI_COMM_WORLD,geometry,local,mode,2,scaled,dyn711)&&okay;
      if(!okay)MPI_Abort(MPI_COMM_WORLD,2);
      std::size_t i{};
      for(int z=0;z<local.cells.z;++z)for(int y=0;y<local.cells.y;++y)
        for(int x=0;x<local.cells.x;++x,++i)for(unsigned g=0;g<3;++g) {
          const auto j=std::size_t(x+local.begin.x+7*y+35*z);
          partition=std::max(partition,std::abs(distributed[3*i+g]-serial[3*j+g]));
          invariance=std::max(invariance,std::abs(distributed[3*i+g]-scaled[3*i+g]));
          if(distributed[3*i+g]>1.01 && distributed[3*i+g]<15.99)++interior_coefficients;
        }
    }
    double errors[]{partition,invariance};MPI_Allreduce(MPI_IN_PLACE,errors,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&interior_coefficients,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
    okay=okay && errors[0]<1e-11 && errors[1]<1e-11 && interior_coefficients>0;
    if(rank==0)std::cout<<(dyn711 ? "dyn711" : "624cf")<<" dynamic filter partition="<<errors[0]<<" uniform-density-scale="<<errors[1]
        <<" interior-coefficients="<<interior_coefficients<<'\n';
    }
  }
  int result=okay?0:1;MPI_Allreduce(MPI_IN_PLACE,&result,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Finalize();return result;
}
