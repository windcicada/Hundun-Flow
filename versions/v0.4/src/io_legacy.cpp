// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "hundun/v04_io.hpp"
#include "io_output_detail.hpp"
#include <algorithm>
#include <array>
#include <climits>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace hundun::v04 {
namespace {
using namespace detail;
using Box=std::array<int,6>;
std::string filename(std::uint64_t step,int rank) {
  return "s"+std::to_string(step)+"r"+std::to_string(rank)+".vtk";
}
void ascii(std::vector<std::uint8_t>& out,const std::string& text) {
  out.insert(out.end(),text.begin(),text.end());
}
void real(std::vector<std::uint8_t>& out,double value) {
  std::uint64_t bits;std::memcpy(&bits,&value,8);
  for(int i=56;i>=0;i-=8)out.push_back(static_cast<std::uint8_t>(bits>>i));
}
std::string field_name(std::string_view name) {
  if(name=="U")return "Velocity";
  if(name=="pi")return "PressureGauge";
  if(name=="T")return "Temperature";
  if(name=="rho")return "Density";
  if(name=="h")return "Enthalpy";
  std::string out(name);
  for(auto& c:out) if(!std::isalnum(static_cast<unsigned char>(c)) && c!='_')c='_';
  return out;
}
}
Status VisitWriter::write_legacy(MPI_Comm comm,const std::filesystem::path& directory,
    const IoServicePlan& services,const CommittedOutputSnapshot& snapshot,
    IoFailureContext* failure) noexcept {
  IoFailureCapture capture(failure);
  int rank{},ranks{};
  if(comm==MPI_COMM_NULL || MPI_Comm_rank(comm,&rank)!=MPI_SUCCESS ||
     MPI_Comm_size(comm,&ranks)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,kOutputInput};
  auto status=output_collective_stage(comm,[&]() -> Status {
    if(directory.empty())return {StatusCode::invalid_plan,kOutputInput};
    return validate_output_snapshot(services,RuntimeServiceKind::visit,snapshot);
  },&capture.context);
  if(!status)return status;
  const auto capacity=output_service(services,RuntimeServiceKind::visit)->maximum_staging_bytes_per_rank;
  std::string index;
  status=output_collective_stage(comm,[&]() -> Status {
    if(rank!=0)return {};
    const auto index_path=directory/"solution.visit";
    const std::string first="!NBLOCKS "+std::to_string(ranks)+"\n";
    index.clear();
    if(std::filesystem::exists(index_path)) {
      std::ifstream old(index_path,std::ios::binary|std::ios::ate);
      const auto length=old.tellg();
      if(!old || length<0 || static_cast<std::uint64_t>(length)>capacity/8)
        return {StatusCode::io_failure,kOutputCapacity};
      index.resize(static_cast<std::size_t>(length));old.seekg(0);
      if(!old.read(index.data(),length) || index.rfind(first,0)!=0)
        return {StatusCode::invalid_plan,kOutputInput};
      std::istringstream prior(index);std::string line;std::getline(prior,line);
      double previous=-1.;std::uint64_t previous_step=0;
      while(std::getline(prior,line)) {
        if(line.rfind("!TIME ",0)!=0)return {StatusCode::invalid_plan,kOutputInput};
        std::istringstream time_text(line.substr(6));double time;std::string extra;
        if(!(time_text>>time) || (time_text>>extra) || !std::isfinite(time) || time<=previous)
          return {StatusCode::invalid_plan,kOutputInput};
        previous=time;
        for(int r=0;r<ranks;++r) {
          if(!std::getline(prior,line) || line.empty() || line[0]!='s')
            return {StatusCode::invalid_plan,kOutputInput};
          const auto separator=line.find('r');
          if(separator==std::string::npos)return {StatusCode::invalid_plan,kOutputInput};
          std::uint64_t step{};
          const auto parsed=std::from_chars(line.data()+1,line.data()+separator,step);
          if(parsed.ec!=std::errc{} || parsed.ptr!=line.data()+separator ||
             line!=filename(step,r) || (r==0 ? step<=previous_step : step!=previous_step))
            return {StatusCode::invalid_plan,kOutputInput};
          previous_step=step;
        }
      }
      if(snapshot.time<=previous || snapshot.step<=previous_step)
        return {StatusCode::invalid_plan,kOutputInput};
    } else index=first;
    std::ostringstream frame;frame<<std::setprecision(17)<<"!TIME "<<snapshot.time<<'\n';
    for(int r=0;r<ranks;++r)frame<<filename(snapshot.step,r)<<'\n';
    index+=frame.str();
    if(index.size()>capacity/8)return {StatusCode::invalid_plan,kOutputCapacity};
    return {};
  },&capture.context);
  if(!status)return status;
  const auto begin=snapshot.patch.begin,cells=snapshot.patch.cells,global=snapshot.geometry->global_cells();
  const Box own{begin.x,begin.y,begin.z,cells.x,cells.y,cells.z};
  const std::array<int,3> full{global.x,global.y,global.z};
  std::array<int,3> dims{};
  std::vector<Box> boxes;
  std::vector<double> values,send,receive;
  std::vector<std::uint8_t> bytes;
  std::size_t components{},points{};
  status=output_collective_stage(comm,[&]() -> Status {
    boxes.resize(ranks);
    for(std::size_t f=0;f<snapshot.fields.size;++f)components+=snapshot.fields.data[f].values.components;
    points=1;
    for(int a=0;a<3;++a) {
      dims[a]=own[a+3]+(own[a]+own[a+3]<full[a] ? 1 : 0);
      if(static_cast<std::size_t>(dims[a])>SIZE_MAX/points)return {StatusCode::invalid_plan,kOutputCapacity};
      points*=dims[a];
    }
    // Bounds the field buffer, both face buffers, binary payload and headers.
    if(components>SIZE_MAX-6 || points>capacity/(5*(components+6)*sizeof(double)))
      return {StatusCode::invalid_plan,kOutputCapacity};
    values.resize(points*components);
    send.resize(points*components);receive.resize(points*components);
    bytes.reserve(points*(components+6)*sizeof(double)+16384);
    return {};
  },&capture.context);
  if(!status)return status;
  if(MPI_Allgather(own.data(),6,MPI_INT,boxes.data(),6,MPI_INT,comm)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,kOutputCollective};
  const auto cell=[&](int x,int y,int z){return static_cast<std::size_t>(x)+dims[0]*(static_cast<std::size_t>(y)+dims[1]*static_cast<std::size_t>(z));};
  for(int z=0;z<cells.z;++z)for(int y=0;y<cells.y;++y)for(int x=0;x<cells.x;++x) {
    std::size_t c=0;
    for(std::size_t f=0;f<snapshot.fields.size;++f)for(unsigned d=0;d<snapshot.fields.data[f].values.components;++d)
      values[cell(x,y,z)*components+c++]=snapshot.fields.data[f].values.unchecked({x,y,z},d);
  }
  for(int axis=0;axis<3;++axis) {
    int low=MPI_PROC_NULL,high=MPI_PROC_NULL;
    status=output_collective_stage(comm,[&]() -> Status {
      for(int r=0;r<ranks;++r) {
        const auto& b=boxes[r];bool aligned=true;
        for(int a=0;a<3;++a)if(a!=axis)aligned&=b[a]==own[a] && b[a+3]==own[a+3];
        if(!aligned)continue;
        if(b[axis]+b[axis+3]==own[axis])low=r;
        if(own[axis]+own[axis+3]==b[axis])high=r;
      }
      if((own[axis]>0 && low==MPI_PROC_NULL) ||
         (own[axis]+own[axis+3]<full[axis] && high==MPI_PROC_NULL))
        return {StatusCode::invalid_plan,kOutputInput};
      return {};
    },&capture.context);
    if(!status)return status;
    std::size_t count=0;
    for(int z=0;z<dims[2];++z)for(int y=0;y<dims[1];++y)for(int x=0;x<dims[0];++x) {
      const int coordinate[]{x,y,z};if(coordinate[axis]!=0)continue;
      for(std::size_t c=0;c<components;++c)send[count++]=values[cell(x,y,z)*components+c];
    }
    // Plan capacities make the allocation cold; chunk payloads for MPI int counts.
    std::uint64_t maximum_count=count;
    if(MPI_Allreduce(MPI_IN_PLACE,&maximum_count,1,MPI_UINT64_T,MPI_MAX,comm)!=MPI_SUCCESS)
      return {StatusCode::mpi_failure,kOutputCollective};
    for(std::size_t offset=0;offset<maximum_count;offset+=INT_MAX) {
      const auto start=std::min(offset,count);
      const int n=static_cast<int>(std::min<std::size_t>(count-start,INT_MAX));
      const int result=MPI_Sendrecv(send.data()+start,n,MPI_DOUBLE,low,180+axis,
          receive.data()+start,n,MPI_DOUBLE,high,180+axis,comm,MPI_STATUS_IGNORE);
      status=output_collective_status(comm,result==MPI_SUCCESS ? Status{} : Status{StatusCode::mpi_failure,kOutputCollective});
      if(!status)return status;
    }
    if(high!=MPI_PROC_NULL) {
      std::size_t i=0;
      for(int z=0;z<dims[2];++z)for(int y=0;y<dims[1];++y)for(int x=0;x<dims[0];++x) {
        const int coordinate[]{x,y,z};if(coordinate[axis]!=dims[axis]-1)continue;
        for(std::size_t c=0;c<components;++c)values[cell(x,y,z)*components+c]=receive[i++];
      }
    }
  }
  status=output_create_directory(comm,rank,directory,&capture.context);
  if(!status)return status;
  status=output_collective_stage(comm,[&]() -> Status {
    std::ostringstream header;header<<std::setprecision(17);
    header<<"# vtk DataFile Version 3.0\nHundun accepted cell centres; SI; p_ref="<<snapshot.pressure_reference
          <<"\nBINARY\nDATASET STRUCTURED_GRID\nDIMENSIONS "<<dims[0]<<' '<<dims[1]<<' '<<dims[2]
          <<"\nPOINTS "<<points<<" double\n";
    ascii(bytes,header.str());
    const Span<const double> coordinates[]{snapshot.geometry->x().centres(),snapshot.geometry->y().centres(),snapshot.geometry->z().centres()};
    for(int z=0;z<dims[2];++z)for(int y=0;y<dims[1];++y)for(int x=0;x<dims[0];++x) {
      real(bytes,coordinates[0].data[begin.x+x]);real(bytes,coordinates[1].data[begin.y+y]);real(bytes,coordinates[2].data[begin.z+z]);
    }
    ascii(bytes,"\nPOINT_DATA "+std::to_string(points)+"\n");
    std::size_t component=0;
    for(std::size_t field=0;field<snapshot.fields.size;++field) {
      const auto& f=snapshot.fields.data[field];
      const auto name=field_name(f.stable_name);const auto width=f.values.components;
      if(width==3)ascii(bytes,"VECTORS "+name+" double\n");
      else if(width<=4)ascii(bytes,"SCALARS "+name+" double "+std::to_string(width)+"\nLOOKUP_TABLE default\n");
      else ascii(bytes,"FIELD FieldData 1\n"+name+" "+std::to_string(width)+" "+std::to_string(points)+" double\n");
      for(std::size_t p=0;p<points;++p)for(unsigned c=0;c<width;++c)real(bytes,values[p*components+component+c]);
      ascii(bytes,"\n");
      if(f.stable_name=="pi") {
        ascii(bytes,"SCALARS Pressure double 1\nLOOKUP_TABLE default\n");
        for(std::size_t p=0;p<points;++p)real(bytes,snapshot.pressure_reference+values[p*components+component]);
        ascii(bytes,"\n");
      }
      component+=width;
    }
    ascii(bytes,"SCALARS vtkGhostType unsigned_char 1\nLOOKUP_TABLE default\n");
    for(int z=0;z<dims[2];++z)for(int y=0;y<dims[1];++y)for(int x=0;x<dims[0];++x)
      bytes.push_back(x>=cells.x || y>=cells.y || z>=cells.z ? 1 : 0);
    ascii(bytes,"\n");
    if(bytes.size()>capacity)return {StatusCode::invalid_plan,kOutputCapacity};
    const auto name=filename(snapshot.step,rank);
    if(!output_write_file(directory/(name+".tmp"),bytes,false,&capture.context))
      return {StatusCode::io_failure,kOutputFile};
    std::filesystem::rename(directory/(name+".tmp"),directory/name);
    return {};
  },&capture.context);
  if(!status)return status;
  return output_collective_stage(comm,[&]() -> Status {
    if(rank!=0)return {};
    if(!output_write_file(directory/"solution.tmp",index,false,&capture.context))
      return {StatusCode::io_failure,kOutputFile};
    std::filesystem::rename(directory/"solution.tmp",directory/"solution.visit");
    return output_sync_directory(directory,&capture.context) ? Status{} : Status{StatusCode::io_failure,kOutputFile};
  },&capture.context);
}
} // namespace hundun::v04
