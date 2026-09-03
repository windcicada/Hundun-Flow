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
                             std::string_view json_payload) noexcept try {
  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      monitor_file.empty() || monitor_file.parent_path().empty() ||
      json_payload.size() < 2U || json_payload.front() != '{' ||
      json_payload.back() != '}' ||
      json_payload.find('\n') != std::string_view::npos ||
      json_payload.find('\r') != std::string_view::npos)
    return {StatusCode::invalid_plan, detail::kOutputInput};
  Status status = detail::validate_output_snapshot(
      services, RuntimeServiceKind::monitor, snapshot);
  status = detail::output_collective_status(communicator, status);
  if (!status) return status;
  status = detail::output_create_directory(communicator, rank,
                                           monitor_file.parent_path());
  if (!status) return status;
  if (rank == 0) {
    std::ostringstream line;
    line.imbue(std::locale::classic());
    line << std::setprecision(17) << "{\"step\":" << snapshot.step
         << ",\"time\":" << snapshot.time << ",\"plan\":"
         << snapshot.plan << ",\"payload\":" << json_payload << "}\n";
    const std::string text = line.str();
    const RuntimeServiceCapacity* capacity =
        detail::output_service(services, RuntimeServiceKind::monitor);
    status = text.size() <= capacity->maximum_staging_bytes_per_rank &&
                     detail::output_write_file(monitor_file, text, true)
                 ? Status{}
                 : Status{StatusCode::io_failure, detail::kOutputFile};
  }
  return detail::output_collective_status(communicator, status);
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, detail::kOutputCapacity};
} catch (...) {
  return {StatusCode::io_failure, detail::kOutputFile};
}

}  // namespace hundun::v04
