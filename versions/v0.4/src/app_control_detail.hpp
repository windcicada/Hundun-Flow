// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once

#include "io_output_detail.hpp"
#include "hundun/v04_app.hpp"
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

// A local writer also serves in-step observers. Their failures are converged
// by the application after advance; filesystem errors leave physics untouched.
inline Status application_status_local(
    const std::filesystem::path& directory, const char* phase,
    std::uint64_t step, double time,
    const DriverAttemptProgress* progress = nullptr) noexcept {
  try {
    std::filesystem::create_directories(directory);
    std::ostringstream json;
    json.imbue(std::locale::classic());
    const auto now=std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    json << std::setprecision(17) << "{\"phase\":\"" << phase
         << "\",\"step\":" << step << ",\"time\":" << time
         << ",\"updated_unix_ms\":" << now;
    if (progress) {
      const auto& failure=progress->previous_failure;
      json << ",\"target_step\":" << progress->proposal.accepted_step+1U
           << ",\"attempt\":" << progress->attempt
           << ",\"coupling_sweep\":" << progress->coupling_sweep
           << ",\"dt\":" << progress->proposal.dt
           << ",\"retry_kind\":\"" << (progress->coupling_sweep>1 ? "scalar_coupling" :
               progress->attempt>1 ? "time_step" : "none")
           << "\",\"previous_failure\":{\"code\":" << unsigned(failure.failure.code)
           << ",\"detail\":" << failure.failure.detail
           << ",\"stage\":" << failure.stage
           << ",\"attempt\":" << failure.attempt
           << ",\"dt\":" << failure.dt << '}';
    }
    json << "}\n";
    const auto temporary=directory/"status.tmp";
    if (!output_write_file(temporary,json.str()))
      return {StatusCode::io_failure,kOutputFile};
    std::filesystem::rename(temporary,directory/"status.json");
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure,kOutputCapacity};
  } catch (...) {
    return {StatusCode::io_failure,kOutputFile};
  }
}

inline Status application_status(MPI_Comm comm, int rank,
    const std::filesystem::path& directory, const char* phase,
    std::uint64_t step, double time) {
  return output_collective_status(comm,rank==0
      ? application_status_local(directory,phase,step,time) : Status{});
}

struct ApplicationAttemptObserver {
  const std::filesystem::path* directory{};
  int rank{};
  Status status{};
  static void observe(void* context, const DriverAttemptProgress& progress) noexcept {
    auto& self=*static_cast<ApplicationAttemptObserver*>(context);
    if (self.rank!=0 || !self.status) return;
    self.status=application_status_local(*self.directory,
        progress.previous_failure.attempt ? "retrying" : "solving",
        progress.proposal.accepted_step,progress.proposal.time,&progress);
  }
};

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
