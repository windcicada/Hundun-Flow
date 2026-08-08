// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_collective_status.hpp"

#include "hundun/rt_error.hpp"
#include "rt_mpi_error_detail.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <utility>

namespace hundun::runtime {
CollectiveStatus collective_status(const MpiContext& context, bool local_ok,
                                   std::string_view local_message) {
  detail::require_mpi_active("collect MPI status");
  const MPI_Comm comm = context.comm();
  if (comm == MPI_COMM_NULL) {
    throw Error("collective_status requires a valid MPI context");
  }

  const int rank = context.rank();
  const int size = context.size();

  const int local_failing_rank = local_ok ? size : rank;
  int failing_rank = size;
  detail::check_mpi(MPI_Allreduce(&local_failing_rank, &failing_rank, 1,
                                  MPI_INT, MPI_MIN, comm),
                    "MPI_Allreduce");
  if (failing_rank == size) {
    return {true, -1, ""};
  }

  const auto local_length = local_message.size();
  const bool fits_u64 =
      local_length <= std::numeric_limits<std::uint64_t>::max();
  const bool local_oversized =
      rank == failing_rank &&
      (!fits_u64 ||
       local_length >
           static_cast<std::string_view::size_type>(
               std::numeric_limits<int>::max()));
  const int local_oversized_value = local_oversized ? 1 : 0;
  int any_oversized = 0;
  detail::check_mpi(MPI_Allreduce(&local_oversized_value, &any_oversized, 1,
                                  MPI_INT, MPI_MAX, comm),
                    "MPI_Allreduce");
  if (any_oversized != 0) {
    throw Error("collective status message exceeds MPI count limit");
  }

  std::uint64_t message_length = 0;
  if (rank == failing_rank) {
    message_length = static_cast<std::uint64_t>(local_length);
  }
  detail::check_mpi(
      MPI_Bcast(&message_length, 1, MPI_UINT64_T, failing_rank, comm),
      "MPI_Bcast");

  const int local_size_oversized =
      message_length > std::numeric_limits<std::size_t>::max() ? 1 : 0;
  int any_size_oversized = 0;
  detail::check_mpi(MPI_Allreduce(&local_size_oversized, &any_size_oversized,
                                  1, MPI_INT, MPI_MAX, comm),
                    "MPI_Allreduce");
  if (any_size_oversized != 0) {
    throw Error("collective status message exceeds local size limit");
  }
  const auto local_size = static_cast<std::size_t>(message_length);

  std::string message;
  bool allocation_ok = true;
  try {
    message.resize(local_size);
  } catch (...) {
    allocation_ok = false;
  }
  const int local_allocation_ok = allocation_ok ? 1 : 0;
  int every_allocation_ok = 0;
  detail::check_mpi(MPI_Allreduce(&local_allocation_ok,
                                  &every_allocation_ok, 1, MPI_INT, MPI_MIN,
                                  comm),
                    "MPI_Allreduce");
  if (every_allocation_ok == 0) {
    throw Error("unable to allocate collective status message");
  }

  if (rank == failing_rank) {
    std::copy(local_message.begin(), local_message.end(), message.begin());
  }
  const int byte_count = static_cast<int>(message_length);
  detail::check_mpi(
      MPI_Bcast(message.data(), byte_count, MPI_BYTE, failing_rank, comm),
      "MPI_Bcast");
  return {false, failing_rank, std::move(message)};
}

void require_expected_ranks(const MpiContext& context,
                            std::optional<int> expected_ranks) {
  const int size = context.size();

  const bool local_ok = !expected_ranks || *expected_ranks == size;
  std::array<char, 96> buffer{};
  std::string_view local_message;
  if (!local_ok) {
    const int count = std::snprintf(buffer.data(), buffer.size(),
                                    "expected MPI rank count %d, got %d",
                                    *expected_ranks, size);
    if (count < 0 || static_cast<std::size_t>(count) >= buffer.size()) {
      local_message = "MPI rank count mismatch";
    } else {
      local_message =
          std::string_view{buffer.data(), static_cast<std::size_t>(count)};
    }
  }

  const CollectiveStatus status =
      collective_status(context, local_ok, local_message);
  if (!status.ok) {
    throw Error(status.message);
  }
}

}  // namespace hundun::runtime
