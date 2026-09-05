// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include <fcntl.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "field_view_interval_detail.hpp"
#include "hundun/v04_io.hpp"

namespace hundun::v04::detail {

inline constexpr std::uint32_t kOutputInput = 10401U;
inline constexpr std::uint32_t kOutputCapacity = 10402U;
inline constexpr std::uint32_t kOutputCollective = 10403U;
inline constexpr std::uint32_t kOutputFile = 10404U;

struct IoFailureCapture {
  IoFailureContext context{};
  IoFailureContext* destination;
  explicit IoFailureCapture(IoFailureContext* out) noexcept : destination(out) {
    if (out != nullptr) *out = {};
  }
  ~IoFailureCapture() noexcept {
    if (destination != nullptr) *destination = context;
  }
};

inline void output_record_failure(IoFailureContext* out,
                                  IoFailureOperation operation, int error,
                                  const std::filesystem::path& path) noexcept {
  if (out == nullptr || out->valid)
    return;  // Preserve the first failed syscall.
  out->valid = true;
  out->operation = operation;
  out->system_error = error;
  const std::size_t length = std::strlen(path.c_str());
  const std::size_t copied = std::min(length, out->path.size() - 1U);
  std::memcpy(out->path.data(), path.c_str(), copied);
  out->path[copied] = '\0';
  out->path_truncated = copied != length;
}

inline const RuntimeServiceCapacity* output_service(
    const IoServicePlan& plan, RuntimeServiceKind kind) noexcept {
  for (std::size_t index = 0U; index < plan.services().size; ++index) {
    if (plan.services().data[index].kind == kind)
      return &plan.services().data[index];
  }
  return nullptr;
}

inline bool same_output_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

inline Status validate_output_snapshot(
    const IoServicePlan& plan, RuntimeServiceKind kind,
    const CommittedOutputSnapshot& snapshot) noexcept {
  const RuntimeServiceCapacity* service = output_service(plan, kind);
  if (service == nullptr || plan.fingerprint() == 0U ||
      snapshot.geometry == nullptr || snapshot.plan == 0U ||
      snapshot.schema == 0U || !snapshot.committed ||
      snapshot.fields.data == nullptr ||
      snapshot.fields.size != plan.snapshot_fields().size ||
      snapshot.patch.cells.x <= 0 || snapshot.patch.cells.y <= 0 ||
      snapshot.patch.cells.z <= 0)
    return {StatusCode::invalid_plan, kOutputInput};
  std::size_t cells = static_cast<std::size_t>(snapshot.patch.cells.x);
  if (cells > std::numeric_limits<std::size_t>::max() /
                  static_cast<std::size_t>(snapshot.patch.cells.y))
    return {StatusCode::invalid_plan, kOutputInput};
  cells *= static_cast<std::size_t>(snapshot.patch.cells.y);
  if (cells > std::numeric_limits<std::size_t>::max() /
                  static_cast<std::size_t>(snapshot.patch.cells.z))
    return {StatusCode::invalid_plan, kOutputInput};
  cells *= static_cast<std::size_t>(snapshot.patch.cells.z);
  std::size_t components = 0U;
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const SnapshotFieldSpec sealed = plan.snapshot_fields().data[index];
    const SnapshotFieldView field = snapshot.fields.data[index];
    FieldStorageInterval interval{};
    if (field.stable_name.empty() || field.values.field != sealed.field ||
        field.values.components != sealed.components ||
        field.accepted_revision == 0U ||
        field.accepted_revision != field.values.revision ||
        !same_output_shape(field.values.interior, snapshot.patch.cells) ||
        !field_storage_interval(field.values, interval) ||
        sealed.components >
            std::numeric_limits<std::size_t>::max() - components)
      return {StatusCode::invalid_plan, kOutputInput};
    components += sealed.components;
  }
  if (components != 0U &&
      cells > std::numeric_limits<std::size_t>::max() / components)
    return {StatusCode::invalid_plan, kOutputInput};
  const std::size_t doubles = cells * components;
  if (doubles > std::numeric_limits<std::size_t>::max() / sizeof(double) ||
      doubles * sizeof(double) > service->maximum_snapshot_bytes_per_rank)
    return {StatusCode::invalid_plan, kOutputCapacity};
  return {};
}

inline Status output_collective_status(
    MPI_Comm communicator, Status local,
    IoFailureContext* failure = nullptr) noexcept {
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0)
    return {StatusCode::mpi_failure, kOutputCollective};
  const int candidate = local ? size : rank;
  int failing = size;
  if (MPI_Allreduce(&candidate, &failing, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kOutputCollective};
  if (failing == size) return {};
  std::uint64_t wire = 0U;
  if (rank == failing)
    wire = (static_cast<std::uint64_t>(local.code) << 32U) | local.detail;
  if (MPI_Bcast(&wire, 1, MPI_UINT64_T, failing, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kOutputCollective};
  if (failure != nullptr) {
    std::array<int, 4U> metadata{
        {failure->valid ? 1 : 0, static_cast<int>(failure->operation),
         failure->system_error, failure->path_truncated ? 1 : 0}};
    if (MPI_Bcast(metadata.data(), static_cast<int>(metadata.size()), MPI_INT,
                  failing, communicator) != MPI_SUCCESS ||
        MPI_Bcast(failure->path.data(), static_cast<int>(failure->path.size()),
                  MPI_CHAR, failing, communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kOutputCollective};
    failure->valid = metadata[0U] != 0;
    failure->operation = static_cast<IoFailureOperation>(metadata[1U]);
    failure->system_error = metadata[2U];
    failure->path_truncated = metadata[3U] != 0;
    failure->rank = failing;
  }
  return {static_cast<StatusCode>(wire >> 32U),
          static_cast<std::uint32_t>(wire)};
}

// The callback must contain local work only. Catch before the next collective
// so a failed allocation cannot make one rank leave the collective sequence.
template <class LocalWork>
inline Status output_collective_stage(
    MPI_Comm communicator, LocalWork&& work,
    IoFailureContext* failure = nullptr) noexcept {
  Status status;
  try {
    status = work();
  } catch (const std::bad_alloc&) {
    status = {StatusCode::allocation_failure, kOutputCapacity};
  } catch (...) {
    status = {StatusCode::io_failure, kOutputFile};
  }
  return output_collective_status(communicator, status, failure);
}

inline bool output_write_file(const std::filesystem::path& path,
                              const std::uint8_t* data, std::size_t bytes,
                              bool append,
                              IoFailureContext* failure = nullptr) noexcept {
  const int flags = append ? O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC
                           : O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC;
  int descriptor = -1;
  do {
    descriptor = ::open(path.c_str(), flags, 0644);
  } while (descriptor < 0 && errno == EINTR);
  if (descriptor < 0) {
    output_record_failure(failure, IoFailureOperation::open, errno, path);
    return false;
  }
  std::size_t cursor = 0U;
  bool okay = true;
  while (cursor < bytes) {
    const ssize_t count = ::write(descriptor, data + cursor, bytes - cursor);
    if (count < 0 && errno == EINTR) continue;
    if (count <= 0) {
      output_record_failure(failure, IoFailureOperation::write,
                            count == 0 ? EIO : errno, path);
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(count);
  }
  if (okay) {
    int synced;
    do {
      synced = ::fsync(descriptor);
    } while (synced != 0 && errno == EINTR);
    if (synced != 0) {
      output_record_failure(failure, IoFailureOperation::sync, errno, path);
      okay = false;
    }
  }
  // Never retry close: the descriptor may already have been released/reused.
  if (::close(descriptor) != 0) {
    output_record_failure(failure, IoFailureOperation::close, errno, path);
    okay = false;
  }
  return okay;
}

inline bool output_write_file(const std::filesystem::path& path,
                              const std::vector<std::uint8_t>& bytes,
                              bool append = false,
                              IoFailureContext* failure = nullptr) noexcept {
  return output_write_file(path, bytes.data(), bytes.size(), append, failure);
}

inline bool output_write_file(const std::filesystem::path& path,
                              std::string_view text, bool append = false,
                              IoFailureContext* failure = nullptr) noexcept {
  return output_write_file(path,
                           reinterpret_cast<const std::uint8_t*>(text.data()),
                           text.size(), append, failure);
}

inline bool output_sync_directory(
    const std::filesystem::path& path,
    IoFailureContext* failure = nullptr) noexcept {
  int descriptor;
  do {
    descriptor = ::open(path.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
  } while (descriptor < 0 && errno == EINTR);
  if (descriptor < 0) {
    output_record_failure(failure, IoFailureOperation::open, errno, path);
    return false;
  }
  int synced;
  do {
    synced = ::fsync(descriptor);
  } while (synced != 0 && errno == EINTR);
  if (synced != 0)
    output_record_failure(failure, IoFailureOperation::sync, errno, path);
  const int closed = ::close(descriptor);
  if (closed != 0)
    output_record_failure(failure, IoFailureOperation::close, errno, path);
  return synced == 0 && closed == 0;
}

inline Status output_create_directory(
    MPI_Comm communicator, int rank, const std::filesystem::path& path,
    IoFailureContext* failure = nullptr) noexcept {
  return output_collective_stage(
      communicator,
      [&]() -> Status {
        if (rank == 0) {
          std::error_code error;
          std::filesystem::create_directories(path, error);
          if (error)
            output_record_failure(failure, IoFailureOperation::create_directory,
                                  error.value(), path);
          if (error || !output_sync_directory(path, failure))
            return {StatusCode::io_failure, kOutputFile};
        }
        return {};
      },
      failure);
}

inline std::string output_json_escape(std::string_view input) {
  std::string out;
  out.reserve(input.size() + 8U);
  for (const char character : input) {
    switch (character) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (static_cast<unsigned char>(character) < 0x20U)
          throw std::invalid_argument("control character");
        out += character;
    }
  }
  return out;
}

}  // namespace hundun::v04::detail
