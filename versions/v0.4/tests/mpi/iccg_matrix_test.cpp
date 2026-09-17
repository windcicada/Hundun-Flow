// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_iccg.hpp"
#include <iostream>
using namespace hundun::v04;
namespace {
struct Field {
  std::vector<double> data;
  FieldView view;
  explicit Field(Int3 n):data(std::size_t(n.x+2)*(n.y+2)*(n.z+2)) {
    view.interior=n;view.components=1;view.ghosts={1,1,1};
    view.stride_y=n.x+2;view.stride_z=std::size_t(n.x+2)*(n.y+2);
    view.component_stride=data.size();view.base=data.data()+1+view.stride_y+view.stride_z;
    view.field=1;view.revision=1;view.storage_identity=100;view.revision_domain=99;
  }
};
bool run(int rank,int ranks,int width,unsigned mode) {
  const Int3 n{width,2,1},global{width*ranks,2,1};
  const MeshPatch patch{{rank*width,0,0},n};
  const std::array<bool,3> periodic{mode==1,false,true};
  const bool solid=mode==2 || mode==3;
  const auto is_solid=[&](int x) {return solid && x==global.x/2;};
  const auto flat=[&](Int3 c){return std::size_t(c.x)+n.x*(c.y+std::size_t(n.y)*c.z);};
  const auto count=std::size_t(n.x)*n.y*n.z;
  std::vector<detail::ColdPressureRow> rows(count);
  for(int z=0;z<n.z;++z)for(int y=0;y<n.y;++y)for(int x=0;x<n.x;++x) {
    const Int3 c{x,y,z};auto& row=rows[flat(c)];const int gx=x+patch.begin.x;
    if(is_solid(gx)){row.diagonal=1;continue;}
    for(unsigned f=0;f<6;++f) {
      Int3 nb{gx,y,z};auto& at=f/2==0 ? nb.x : f/2==1 ? nb.y : nb.z;
      const int extent=f/2==0 ? global.x : f/2==1 ? global.y : global.z;
      at+=f%2 ? 1:-1;
      if(periodic[f/2])at=(at+extent)%extent;
      if(at<0 || at>=extent || is_solid(nb.x))continue;
      row.neighbour[f]=1;
    }
    if(mode==4 && gx==global.x-1 && y==0 && gx>0)row.neighbour[0]=2;
    for(auto a:row.neighbour)row.diagonal+=a;
    if((mode!=5 && gx==0 && y==0) || (mode==2 && gx==global.x-1 && y==0) || mode==4)
      row.diagonal+=1;
  }
  if(mode==6 && rank==ranks-1)rows.back().diagonal=std::numeric_limits<double>::quiet_NaN();
  Field scratch(n);
  std::vector<std::size_t> parent(count);
  std::vector<std::uint8_t> grounded(count);
  ReductionEngine reductions;
  auto status=ReductionEngine::compile(MPI_COMM_WORLD,ReductionMode::mpi_allreduce,4,reductions);
  if(!status)return false;
  const auto exchange=[&](FieldView field) -> Status {
    double left[2],right[2],from_left[2]{},from_right[2]{};
    for(int y=0;y<2;++y) {left[y]=field.unchecked({0,y,0},0);right[y]=field.unchecked({n.x-1,y,0},0);}
    const int lower=rank>0 ? rank-1 : periodic[0] ? ranks-1 : MPI_PROC_NULL;
    const int upper=rank+1<ranks ? rank+1 : periodic[0] ? 0 : MPI_PROC_NULL;
    int code=MPI_Sendrecv(right,2,MPI_DOUBLE,upper,90,from_left,2,MPI_DOUBLE,lower,90,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    code|=MPI_Sendrecv(left,2,MPI_DOUBLE,lower,91,from_right,2,MPI_DOUBLE,upper,91,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    if(code!=MPI_SUCCESS)return {StatusCode::mpi_failure,1};
    for(int y=0;y<2;++y) {
      field.unchecked({-1,y,0},0)=from_left[y];field.unchecked({n.x,y,0},0)=from_right[y];
      for(int x=0;x<n.x;++x) {
        field.unchecked({x,y,-1},0)=field.unchecked({x,y,0},0);
        field.unchecked({x,y,1},0)=field.unchecked({x,y,0},0);
      }
    }
    return {};
  };
  detail::IccgMatrixAdmissionReport report;
  status=detail::admit_iccg_matrix(rows,patch,global,periodic,scratch.view,
      {parent.data(),parent.size()},{grounded.data(),grounded.size()},exchange,reductions,report);
  const bool expected=mode<=2;
  const unsigned detail=mode==4 ? 17871:mode==6 ? 17869:17872;
  bool okay=bool(status)==expected && report.spd==expected && (expected || status.detail==detail);
  int all=okay;MPI_Allreduce(MPI_IN_PLACE,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(rank==0)std::cout<<"ICCG matrix width="<<width<<" mode="<<mode<<" status="<<unsigned(status.code)<<'/'<<status.detail
      <<" grounding_exchanges="<<report.grounding_exchanges<<" pass="<<all<<'\n';
  return all;
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank{},ranks{};
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  bool okay=true;
  for(int width:{1,3})for(unsigned mode=0;mode<7;++mode) {
    // A separated ungrounded fluid component needs at least three x cells.
    if((mode==2 || mode==3 || mode==4) && width*ranks<3)continue;
    okay=run(rank,ranks,width,mode)&&okay;
  }
  MPI_Finalize();return okay ? 0:1;
}
