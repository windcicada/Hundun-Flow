// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_io.hpp"

#include "io_output_detail.hpp"

#include <mpi.h>

#include <iomanip>
#include <locale>
#include <new>
#include <sstream>
#include <string>

namespace hundun::v04 {

Status MonitorWriter::append(MPI_Comm communicator,
                             const std::filesystem::path& monitor_file,
                             const IoServicePlan& services,
                             const CommittedOutputSnapshot& snapshot,
                             std::string_view json_payload) noexcept {
  int rank = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS)
    return {StatusCode::invalid_plan, detail::kOutputInput};
  std::filesystem::path parent;
  Status status =
      detail::output_collective_stage(communicator, [&]() -> Status {
        parent = monitor_file.parent_path();
        if (monitor_file.empty() || parent.empty() ||
            json_payload.size() < 2U || json_payload.front() != '{' ||
            json_payload.back() != '}' ||
            json_payload.find('\n') != std::string_view::npos ||
            json_payload.find('\r') != std::string_view::npos)
          return {StatusCode::invalid_plan, detail::kOutputInput};
        return detail::validate_output_snapshot(
            services, RuntimeServiceKind::monitor, snapshot);
      });
  if (!status) return status;
  status = detail::output_create_directory(communicator, rank, parent);
  if (!status) return status;
  return detail::output_collective_stage(communicator, [&]() -> Status {
    if (rank != 0) return {};
    std::ostringstream line;
    line.exceptions(std::ios::badbit | std::ios::failbit);
    line.imbue(std::locale::classic());
    line << std::setprecision(17) << "{\"step\":" << snapshot.step
         << ",\"time\":" << snapshot.time << ",\"plan\":"
         << snapshot.plan << ",\"payload\":" << json_payload << "}\n";
    const std::string text = line.str();
    const RuntimeServiceCapacity* capacity =
        detail::output_service(services, RuntimeServiceKind::monitor);
    return text.size() <= capacity->maximum_staging_bytes_per_rank &&
                   detail::output_write_file(monitor_file, text, true)
               ? Status{}
               : Status{StatusCode::io_failure, detail::kOutputFile};
  });
}

}  // namespace hundun::v04
