// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "hundun/v04_app.hpp"
#include "io_output_detail.hpp"
#include "yyjson.h"

#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>

namespace hundun::v04 {
namespace {
constexpr std::uint32_t kRunInput=10510U;
constexpr std::size_t kRunLimit=65536U;
struct Document {
  yyjson_doc* value{};
  ~Document(){yyjson_doc_free(value);}
};
Status run_tokens(const std::filesystem::path& file, const std::string& case_root,
                  std::vector<std::string>& tokens) {
  std::ifstream input(file,std::ios::binary|std::ios::ate);
  if(!input) return {StatusCode::io_failure,kRunInput};
  const auto size=input.tellg();
  if(size<=0 || size>static_cast<std::streamoff>(kRunLimit))
    return {StatusCode::invalid_case,kRunInput};
  std::string text(static_cast<std::size_t>(size),'\0');
  input.seekg(0);
  if(!input.read(text.data(),size))return {StatusCode::io_failure,kRunInput};
  Document doc{yyjson_read(text.data(),text.size(),0)};
  auto* root=doc.value ? yyjson_doc_get_root(doc.value) : nullptr;
  if(!yyjson_is_obj(root))return {StatusCode::invalid_case,kRunInput};
  std::set<std::string> keys;
  std::size_t idx,max;
  yyjson_val *key,*value;
  const std::set<std::string> allowed{"mode","steps","end_time","output","restart","output_interval",
      "restart_interval","monitor_interval","diagnostics_interval","max_dt","initial_state","visit_format"};
  yyjson_obj_foreach(root,idx,max,key,value) {
    const std::string name(yyjson_get_str(key),yyjson_get_len(key));
    if(!allowed.count(name) || !keys.insert(name).second)
      return {StatusCode::invalid_case,kRunInput};
  }
  const auto string_value=[](yyjson_val* v,std::string& out) {
    if(!yyjson_is_str(v) || yyjson_get_len(v)==0)return false;
    out.assign(yyjson_get_str(v),yyjson_get_len(v));
    return out.find('\0')==std::string::npos;
  };
  std::string mode,output,restart;
  if(!string_value(yyjson_obj_get(root,"mode"),mode) ||
     (mode!="new" && mode!="restart") ||
     !string_value(yyjson_obj_get(root,"output"),output))
    return {StatusCode::invalid_case,kRunInput};
  if(mode=="restart") {
    if(!string_value(yyjson_obj_get(root,"restart"),restart) || keys.count("initial_state"))
      return {StatusCode::invalid_case,kRunInput};
  } else if(keys.count("restart")) return {StatusCode::invalid_case,kRunInput};
  auto* steps=yyjson_obj_get(root,"steps");
  auto* end_time=yyjson_obj_get(root,"end_time");
  if((steps && (!yyjson_is_uint(steps) || yyjson_get_uint(steps)==0)) ||
     (!steps && !end_time) || (end_time && (!yyjson_is_num(end_time) ||
       !std::isfinite(yyjson_get_num(end_time)) || yyjson_get_num(end_time)<=0)))
    return {StatusCode::invalid_case,kRunInput};
  const auto path=[&](const std::string& name) {
    const std::filesystem::path p{name};
    return (p.is_absolute() ? p : file.parent_path()/p).lexically_normal().string();
  };
  tokens={"run",case_root,"--output",path(output),"--steps",std::to_string(steps ? yyjson_get_uint(steps) : INT64_MAX)};
  if(end_time) {
    std::ostringstream out;out<<std::setprecision(17)<<yyjson_get_num(end_time);
    tokens.push_back("--until");tokens.push_back(out.str());
  }
  if(mode=="restart") { tokens.push_back("--restart");tokens.push_back(path(restart)); }
  for(auto name:{"output_interval","restart_interval","monitor_interval","diagnostics_interval"}) {
    if(auto* v=yyjson_obj_get(root,name)) {
      if(!yyjson_is_uint(v))return {StatusCode::invalid_case,kRunInput};
      std::string flag=name;
      std::replace(flag.begin(),flag.end(),'_','-');
      tokens.push_back("--"+flag);tokens.push_back(std::to_string(yyjson_get_uint(v)));
    }
  }
  if(auto* v=yyjson_obj_get(root,"visit_format")) {
    std::string format;
    if(!string_value(v,format) || (format!="legacy" && format!="xml"))
      return {StatusCode::invalid_case,kRunInput};
    tokens.push_back("--visit-format");tokens.push_back(format);
  }
  if(auto* v=yyjson_obj_get(root,"max_dt")) {
    const double dt=yyjson_get_num(v);
    if(!yyjson_is_num(v) || !std::isfinite(dt) || dt<=0)
      return {StatusCode::invalid_case,kRunInput};
    std::ostringstream out;out<<std::setprecision(17)<<dt;
    tokens.push_back("--max-dt");tokens.push_back(out.str());
  }
  if(auto* v=yyjson_obj_get(root,"initial_state")) {
    if(!yyjson_is_arr(v) || yyjson_arr_size(v)<5 || yyjson_arr_size(v)>69)
      return {StatusCode::invalid_case,kRunInput};
    std::ostringstream out;out<<std::setprecision(17);
    std::size_t i,count; yyjson_val* component;
    yyjson_arr_foreach(v,i,count,component) {
      if(!yyjson_is_num(component) || !std::isfinite(yyjson_get_num(component)))
        return {StatusCode::invalid_case,kRunInput};
      if(i)out<<',';
      out<<yyjson_get_num(component);
    }
    tokens.push_back("--initial-state");tokens.push_back(out.str());
  }
  return {};
}
} // namespace

Status expand_run_arguments(MPI_Comm comm,int argc,char** argv,
                            std::vector<std::string>& arguments) noexcept {
  int rank{};
  if(comm==MPI_COMM_NULL || MPI_Comm_rank(comm,&rank)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,kRunInput};
  std::string wire;
  auto status=detail::output_collective_stage(comm,[&]() -> Status {
    if(rank!=0)return {};
    std::vector<std::string> tokens;
    bool configured=argc==1;
    std::filesystem::path file=std::filesystem::current_path()/"run.json";
    std::string case_root=std::filesystem::current_path().string();
    std::vector<std::string> overrides;
    if(argc>=3 && std::string_view(argv[1])=="run") {
      case_root=argv[2];
      for(int i=3;i<argc;++i) {
        if(std::string_view(argv[i])=="--config") {
          if(configured || i+1==argc)return {StatusCode::invalid_case,kRunInput};
          file=std::filesystem::absolute(argv[++i]);configured=true;
        } else overrides.emplace_back(argv[i]);
      }
    }
    if(configured) {
      auto result=run_tokens(file,case_root,tokens);
      if(!result)return result;
      // A CLI key replaces its matching run.json key. Cross-key conflicts
      // (for example fresh initial state plus Restart) remain explicit errors.
      std::set<std::string> flags;
      for(const auto& word:overrides) if(word.rfind("--",0)==0)flags.insert(word);
      std::vector<std::string> merged{tokens[0],tokens[1]};
      for(std::size_t i=2;i<tokens.size();i+=2) if(!flags.count(tokens[i])) {
        merged.push_back(tokens[i]);merged.push_back(tokens[i+1]);
      }
      merged.insert(merged.end(),overrides.begin(),overrides.end());
      tokens=std::move(merged);
    } else for(int i=1;i<argc;++i)tokens.emplace_back(argv[i]);
    tokens.insert(tokens.begin(),argv[0]);
    for(const auto& token:tokens) {wire+=token;wire+='\0';}
    if(wire.size()>kRunLimit)return {StatusCode::invalid_case,kRunInput};
    return {};
  });
  if(!status)return status;
  std::uint64_t length=wire.size();
  if(MPI_Bcast(&length,1,MPI_UINT64_T,0,comm)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,kRunInput};
  status=detail::output_collective_stage(comm,[&]() -> Status {wire.resize(length);return {};});
  if(!status)return status;
  if(MPI_Bcast(wire.data(),static_cast<int>(length),MPI_CHAR,0,comm)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,kRunInput};
  return detail::output_collective_stage(comm,[&]() -> Status {
    arguments.clear();
    for(std::size_t start=0;start<wire.size();) {
      const auto end=wire.find('\0',start);
      if(end==std::string::npos)return {StatusCode::invalid_case,kRunInput};
      arguments.push_back(wire.substr(start,end-start));start=end+1;
    }
    return {};
  });
}
} // namespace hundun::v04
