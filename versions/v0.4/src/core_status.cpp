// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_status.hpp"

namespace hundun::v04 {

std::string_view status_message(Status status) noexcept {
  switch (status.code) {
    case StatusCode::ok:
      return "ok";
    case StatusCode::invalid_case:
      return "invalid case";
    case StatusCode::invalid_plan:
      return "invalid plan";
    case StatusCode::allocation_failure:
      return "allocation failure";
    case StatusCode::mpi_failure:
      return "MPI failure";
    case StatusCode::numerical_failure:
      return "numerical failure";
    case StatusCode::rejected_step:
      return "rejected step";
    case StatusCode::io_failure:
      return "I/O failure";
  }
  return "unknown status";
}

}  // namespace hundun::v04
