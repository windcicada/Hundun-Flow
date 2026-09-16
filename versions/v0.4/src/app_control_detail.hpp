// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once

#include "io_output_detail.hpp"
#include <chrono>
#include <iomanip>
#include <sstream>

namespace hundun::v04::detail {

// All filesystem work is root-owned and converged before the next collective.
// Pending requests survive failed output; a new flag arriving during writing
// remains a separate request for the next accepted step.
inline Status application_requests(MPI_Comm comm, int rank,
    const std::filesystem::path& directory, std::array<int, 2>& requests) {
  auto status = output_collective_stage(comm, [&]() -> Status {
    if (rank != 0) return {};
    const char* names[]{"stop", "output"};
    for (unsigned i=0; i<2; ++i) {
      const auto flag=directory/names[i], pending=directory/(std::string(names[i])+".pending");
      if (std::filesystem::exists(pending)) {
        if (!std::filesystem::is_regular_file(pending))
          return {StatusCode::io_failure, kOutputFile};
        requests[i]=1;
      } else if (std::filesystem::exists(flag)) {
        if (!std::filesystem::is_regular_file(flag))
          return {StatusCode::io_failure, kOutputFile};
        std::filesystem::rename(flag,pending);
        requests[i]=1;
      }
    }
    return {};
  });
  if (status && MPI_Bcast(requests.data(),2,MPI_INT,0,comm)!=MPI_SUCCESS)
    return {StatusCode::mpi_failure,kOutputFile};
  return status;
}

inline Status application_status(MPI_Comm comm, int rank,
    const std::filesystem::path& directory, const char* phase,
    std::uint64_t step, double time) {
  return output_collective_stage(comm,[&]() -> Status {
    if (rank != 0) return {};
    std::filesystem::create_directories(directory);
    std::ostringstream json;
    const auto now=std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    json << std::setprecision(17) << "{\"phase\":\"" << phase
         << "\",\"step\":" << step << ",\"time\":" << time
         << ",\"updated_unix_ms\":" << now << "}\n";
    const auto temporary=directory/"status.tmp";
    if (!output_write_file(temporary,json.str()))
      return {StatusCode::io_failure,kOutputFile};
    std::filesystem::rename(temporary,directory/"status.json");
    return {};
  });
}

inline Status application_acknowledge(MPI_Comm comm,int rank,
    const std::filesystem::path& directory,const std::array<int,2>& requests,
    std::uint64_t step,double time) {
  return output_collective_stage(comm,[&]() -> Status {
    if (rank!=0) return {};
    const char* names[]{"stop","output"};
    for(unsigned i=0;i<2;++i) if(requests[i]) {
      std::ostringstream receipt;
      receipt << std::setprecision(17) << "{\"request\":\"" << names[i]
              << "\",\"step\":" << step << ",\"time\":" << time
              << ",\"published\":\"" << (i==0 ? "Restart/current" : "Visit")
              << "\"}\n";
      if(!output_write_file(directory/"control.jsonl",receipt.str(),true))
        return {StatusCode::io_failure,kOutputFile};
      std::filesystem::remove(directory/(std::string(names[i])+".pending"));
    }
    return output_sync_directory(directory) ? Status{} : Status{StatusCode::io_failure,kOutputFile};
  });
}
} // namespace hundun::v04::detail
