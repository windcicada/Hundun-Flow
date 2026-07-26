// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/diagnostic_session.hpp"

#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

#include <mpi.h>

#include <array>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <system_error>

namespace hundun::diagnostics {
namespace {

std::filesystem::path output_path(const std::filesystem::path& directory,
                                  int rank, std::uint64_t step) {
  std::ostringstream name;
  name.imbue(std::locale::classic());
  name << "diagnostics.v1.rank-" << std::setw(6) << std::setfill('0') << rank
       << ".step-" << std::setw(20) << step << ".jsonl";
  return directory / name.str();
}

std::uint64_t stable_hash(std::string_view value) noexcept {
  std::uint64_t result = 1469598103934665603ULL;
  for (char raw : value) {
    const auto byte = static_cast<unsigned char>(raw);
    result ^= byte;
    result *= 1099511628211ULL;
  }
  return result;
}

void require_path_agreement(const runtime::MpiContext& mpi,
                            const std::filesystem::path& path) {
  const auto text = path.lexically_normal().generic_string();
  const std::array<std::uint64_t, 2> local{
      static_cast<std::uint64_t>(text.size()), stable_hash(text)};
  std::array<std::uint64_t, 2> minimum = local;
  std::array<std::uint64_t, 2> maximum = local;
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), minimum.data(), 2, MPI_UINT64_T, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce diagnostic path minimum");
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), maximum.data(), 2, MPI_UINT64_T, MPI_MAX,
                    mpi.comm()),
      "MPI_Allreduce diagnostic path maximum");
  if (minimum != maximum)
    throw runtime::Error("diagnostic output path differs across ranks");
}

void require_schedule_agreement(const runtime::MpiContext& mpi,
                                std::uint64_t step, int write_interval) {
  const std::array<std::uint64_t, 2> local{
      step, static_cast<std::uint64_t>(write_interval)};
  std::array<std::uint64_t, 2> minimum = local;
  std::array<std::uint64_t, 2> maximum = local;
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), minimum.data(), 2, MPI_UINT64_T, MPI_MIN,
                    mpi.comm()),
      "MPI_Allreduce diagnostic schedule minimum");
  runtime::check_mpi_result(
      MPI_Allreduce(local.data(), maximum.data(), 2, MPI_UINT64_T, MPI_MAX,
                    mpi.comm()),
      "MPI_Allreduce diagnostic schedule maximum");
  if (minimum != maximum)
    throw runtime::Error("diagnostic output schedule differs across ranks");
}

}  // namespace

std::size_t DiagnosticBatch::size() const noexcept { return records_.size(); }
bool DiagnosticBatch::empty() const noexcept { return records_.empty(); }

std::string DiagnosticBatch::canonical_json_lines() const {
  std::string result;
  for (const auto& record : records_) {
    result += to_canonical_json(record);
    result.push_back('\n');
  }
  return result;
}

DiagnosticBatchSink::DiagnosticBatchSink(DiagnosticDescriptor descriptor,
                                         DiagnosticRequest request,
                                         DiagnosticBatch& batch)
    : descriptor_(descriptor), request_(std::move(request)), batch_(&batch) {
  validate(descriptor_);
  validate(request_, descriptor_);
}

void DiagnosticBatchSink::submit(const DiagnosticRecord& record) {
  try {
    validate(record, descriptor_, request_);
    batch_->records_.push_back(record);
  } catch (const DiagnosticCollectionError&) {
    throw;
  } catch (const std::exception& error) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", -1,
                                    error.what());
  }
}

DiagnosticSession::DiagnosticSession(std::filesystem::path directory,
                                     int write_interval, int rank)
    : directory_(std::move(directory)),
      write_interval_(write_interval),
      rank_(rank) {
  if (directory_.empty() || write_interval_ < 1 || rank_ < 0)
    throw runtime::Error("invalid diagnostic session configuration");
}

bool DiagnosticSession::due(std::uint64_t step) const noexcept {
  return step != 0U &&
         step % static_cast<std::uint64_t>(write_interval_) == 0U;
}

const std::filesystem::path& DiagnosticSession::directory() const noexcept {
  return directory_;
}
int DiagnosticSession::write_interval() const noexcept {
  return write_interval_;
}
int DiagnosticSession::rank() const noexcept { return rank_; }

void DiagnosticSession::publish(const runtime::MpiContext& mpi,
                                std::uint64_t step,
                                const DiagnosticBatch& batch) const {
  if (rank_ != mpi.rank())
    throw runtime::Error("diagnostic session rank mismatch");
  require_schedule_agreement(mpi, step, write_interval_);
  const auto path = output_path(directory_, rank_, step);
  const auto temporary = std::filesystem::path(path.string() + ".tmp");
  require_path_agreement(mpi, directory_);
  const auto bytes = batch.canonical_json_lines();
  bool local_ok = true;
  std::string message;
  try {
    std::filesystem::create_directories(directory_);
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream)
      throw runtime::Error("unable to open diagnostic staging file");
    stream.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    stream.flush();
    if (!stream)
      throw runtime::Error("unable to stage diagnostic records");
    stream.close();
    if (!stream)
      throw runtime::Error("unable to close diagnostic staging file");
  } catch (const std::exception& error) {
    local_ok = false;
    message = error.what();
  }
  auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    throw runtime::Error(status.message);
  }
  std::error_code error;
  std::filesystem::rename(temporary, path, error);
  status = runtime::collective_status(
      mpi, !error, error ? "unable to publish diagnostic records" : "");
  if (!status.ok) {
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    throw runtime::Error(status.message);
  }
}

}  // namespace hundun::diagnostics
