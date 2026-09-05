// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include <fcntl.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

inline Status output_collective_status(MPI_Comm communicator,
                                       Status local) noexcept {
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
  return {static_cast<StatusCode>(wire >> 32U),
          static_cast<std::uint32_t>(wire)};
}

// The callback must contain local work only. Catch before the next collective
// so a failed allocation cannot make one rank leave the collective sequence.
template <class LocalWork>
inline Status output_collective_stage(MPI_Comm communicator,
                                      LocalWork&& work) noexcept {
  Status status;
  try {
    status = work();
  } catch (const std::bad_alloc&) {
    status = {StatusCode::allocation_failure, kOutputCapacity};
  } catch (...) {
    status = {StatusCode::io_failure, kOutputFile};
  }
  return output_collective_status(communicator, status);
}

inline bool output_write_file(const std::filesystem::path& path,
                              const std::uint8_t* data, std::size_t bytes,
                              bool append) noexcept {
  const int flags = append ? O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC
                           : O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC;
  const int descriptor = ::open(path.c_str(), flags, 0644);
  if (descriptor < 0) return false;
  std::size_t cursor = 0U;
  bool okay = true;
  while (cursor < bytes) {
    const ssize_t count = ::write(descriptor, data + cursor, bytes - cursor);
    if (count <= 0) {
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(count);
  }
  if (okay && ::fsync(descriptor) != 0) okay = false;
  if (::close(descriptor) != 0) okay = false;
  return okay;
}

inline bool output_write_file(const std::filesystem::path& path,
                              const std::vector<std::uint8_t>& bytes,
                              bool append = false) noexcept {
  return output_write_file(path, bytes.data(), bytes.size(), append);
}

inline bool output_write_file(const std::filesystem::path& path,
                              std::string_view text,
                              bool append = false) noexcept {
  return output_write_file(
      path, reinterpret_cast<const std::uint8_t*>(text.data()), text.size(),
      append);
}

inline bool output_sync_directory(const std::filesystem::path& path) noexcept {
  const int descriptor = ::open(path.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
  if (descriptor < 0) return false;
  const bool okay = ::fsync(descriptor) == 0;
  return ::close(descriptor) == 0 && okay;
}

inline Status output_create_directory(MPI_Comm communicator, int rank,
                                      const std::filesystem::path& path) noexcept {
  return output_collective_stage(communicator, [&]() -> Status {
    if (rank == 0) {
      std::error_code error;
      std::filesystem::create_directories(path, error);
      if (error || !output_sync_directory(path))
        return {StatusCode::io_failure, kOutputFile};
    }
    return {};
  });
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
