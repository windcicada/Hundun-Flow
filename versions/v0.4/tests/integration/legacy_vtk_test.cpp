// SPDX-License-Identifier: Apache-2.0
#include "../support/piso_fixture.hpp"
#include "hundun/v04_io.hpp"
#include <iostream>
using namespace hundun::v04;
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank{};MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if(argc!=2){MPI_Finalize();return 2;}
  CartesianMeshSpec spec;spec.kind=GeometryKind::tensor_stretched;
  spec.lower={0,0,0};spec.upper={1,1,1};spec.has_exact_cells=true;spec.exact_cells={17,11,7};
  spec.has_base_spacing=true;spec.base_spacing={1.04/17,1.04/11,1.04/7};
  spec.minimum_spacing={.96/17,.96/11,.96/7};spec.max_growth_ratio=1+.8/17;
  spec.focus_regions.push_back({{.35,.35,.35},{.65,.65,.65},{.98/17,.98/11,.98/7}});
  spec.limits.max_global_cells=17*11*7;spec.limits.max_memory_bytes_per_rank=134217728;
  CartesianGeometryPlan geometry;MeshPatch patch;
  Status status=CartesianGeometryCompiler::compile(MPI_COMM_WORLD,spec,{},geometry,patch);
  if(!status){MPI_Abort(MPI_COMM_WORLD,2);}
  auto u=test::make_field(0,patch.cells,3,0,11,21);
  auto p=test::make_field(1,patch.cells,1,0,12,22);
  auto t=test::make_field(2,patch.cells,1,0,13,23);
  for(int z=0;z<patch.cells.z;++z)for(int y=0;y<patch.cells.y;++y)for(int x=0;x<patch.cells.x;++x) {
    const double X=geometry.x().centres().data[patch.begin.x+x],Y=geometry.y().centres().data[patch.begin.y+y],Z=geometry.z().centres().data[patch.begin.z+z];
    u.view.unchecked({x,y,z},0)=1+X;
    p.view.unchecked({x,y,z},0)=X+2*Y+4*Z;
    t.view.unchecked({x,y,z},0)=300+2*X+3*Y+5*Z;
  }
  const std::array<SnapshotFieldSpec,3> specs{{{0,3},{1,1},{2,1}}};
  const auto n=static_cast<std::size_t>(patch.cells.x)*patch.cells.y*patch.cells.z;
  std::array<RuntimeServiceCapacity,5> capacities{};
  for(unsigned i=0;i<5;++i)capacities[i]={static_cast<RuntimeServiceKind>(i),static_cast<StageId>(200+i),n*5*8,4194304,32};
  IoServicePlan services;status=IoServicePlan::compile({specs.data(),specs.size()},{capacities.data(),capacities.size()},n,services);
  const std::array<SnapshotFieldView,3> fields{{{"U",as_const(u.view),11},{"pi",as_const(p.view),12},{"T",as_const(t.view),13}}};
  CommittedOutputSnapshot snapshot{&geometry,patch,42,43,.125,1,{fields.data(),fields.size()},true,{},101325};
  if(status)status=VisitWriter::write_legacy(MPI_COMM_WORLD,argv[1],services,snapshot);
  snapshot.step=2;snapshot.time=.3;
  if(status)status=VisitWriter::write_legacy(MPI_COMM_WORLD,argv[1],services,snapshot);
  // Already published frames reject replay before writing rank files.
  if(status && VisitWriter::write_legacy(MPI_COMM_WORLD,argv[1],services,snapshot))status={StatusCode::invalid_plan,1};
  if(rank==0)std::cout<<"legacy_vtk="<<bool(status)<<" grid=17x11x7 frames=2 replay=protected\n";
  MPI_Finalize();return status ? 0 : 1;
}
