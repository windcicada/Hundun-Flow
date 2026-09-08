// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_io.hpp"

#include "app_identity_detail.hpp"
#include "field_view_interval_detail.hpp"
#include "io_restart_detail.hpp"
#include "io_output_detail.hpp"

#include <mpi.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

namespace fs = std::filesystem;

constexpr std::uint32_t kRestartInput = 10301U;
constexpr std::uint32_t kRestartCollective = 10302U;
constexpr std::uint32_t kRestartDirectory = 10303U;
constexpr std::uint32_t kRestartRankFile = 10304U;
constexpr std::uint32_t kRestartManifest = 10305U;
constexpr std::uint32_t kRestartIntegrity = 10306U;
constexpr std::uint32_t kRestartMismatch = 10307U;
constexpr std::uint32_t kRestartCoverage = 10308U;
constexpr std::uint32_t kRestartPublication = 10309U;
constexpr std::uint32_t kRestartReadBudget = 10310U;
constexpr std::size_t kCurrentMaximumBytes = 256U;
constexpr std::uint32_t kLegacyFormatVersion = 1U;
constexpr std::uint32_t kExactHistoryFormatVersion = 2U;
constexpr std::uint32_t kSignedHistoryFormatVersion = 3U;
constexpr std::uint32_t kCellRecordFormatVersion = 4U;
constexpr std::uint32_t kVariableCellRecordFormatVersion = 5U;
constexpr std::array<char, 8U> kRankMagic{{'H', '4', 'R', 'A', 'N', 'K', '0', '1'}};
constexpr std::array<char, 8U> kManifestMagic{{'H', '4', 'M', 'A', 'N', 'I', '0', '1'}};
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
std::atomic<int> g_restart_failure_point{
    static_cast<int>(detail::RestartFailurePoint::none)};
std::atomic<int> g_restart_failure_rank{-1};

bool injected(detail::RestartFailurePoint point, int rank) noexcept {
  return g_restart_failure_point.load(std::memory_order_relaxed) ==
             static_cast<int>(point) &&
         g_restart_failure_rank.load(std::memory_order_relaxed) == rank;
}
#endif

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool valid_global_patch(Int3 global, const MeshPatch& patch) noexcept {
  if (global.x <= 0 || global.y <= 0 || global.z <= 0 ||
      patch.begin.x < 0 || patch.begin.y < 0 || patch.begin.z < 0 ||
      patch.cells.x <= 0 || patch.cells.y <= 0 || patch.cells.z <= 0) {
    return false;
  }
  return patch.begin.x <= global.x - patch.cells.x &&
         patch.begin.y <= global.y - patch.cells.y &&
         patch.begin.z <= global.z - patch.cells.z;
}

bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  out = left * right;
  return true;
}

bool cell_count(Int3 cells, std::size_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) return false;
  std::size_t plane = 0U;
  return checked_multiply(static_cast<std::size_t>(cells.x),
                          static_cast<std::size_t>(cells.y), plane) &&
         checked_multiply(plane, static_cast<std::size_t>(cells.z), out);
}

std::uint64_t hash_bytes(const std::uint8_t* data, std::size_t size) noexcept {
  std::uint64_t hash = kFnvOffset;
  for (std::size_t index = 0U; index < size; ++index) {
    hash ^= data[index];
    hash *= kFnvPrime;
  }
  return hash == 0U ? 1U : hash;
}

class Encoder {
 public:
  explicit Encoder(std::size_t capacity = 0U) { data_.reserve(capacity); }
  void bytes(const void* data, std::size_t size) {
    const auto* begin = static_cast<const std::uint8_t*>(data);
    data_.insert(data_.end(), begin, begin + size);
  }
  void u8(std::uint8_t value) { data_.push_back(value); }
  void u16(std::uint16_t value) {
    for (unsigned shift = 0U; shift < 16U; shift += 8U)
      data_.push_back(static_cast<std::uint8_t>(value >> shift));
  }
  void u32(std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      data_.push_back(static_cast<std::uint8_t>(value >> shift));
  }
  void i32(std::int32_t value) { u32(static_cast<std::uint32_t>(value)); }
  void u64(std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      data_.push_back(static_cast<std::uint8_t>(value >> shift));
  }
  void real(double value) {
    std::uint64_t bits = 0U;
    static_assert(sizeof(bits) == sizeof(value), "binary64 restart format");
    std::memcpy(&bits, &value, sizeof(bits));
    u64(bits);
  }
  void int3(Int3 value) {
    i32(value.x);
    i32(value.y);
    i32(value.z);
  }
  void append_integrity() { u64(hash_bytes(data_.data(), data_.size())); }
  const std::vector<std::uint8_t>& data() const noexcept { return data_; }
  std::vector<std::uint8_t> take() noexcept { return std::move(data_); }

 private:
  std::vector<std::uint8_t> data_;
};

class Decoder {
 public:
  explicit Decoder(const std::vector<std::uint8_t>& data) noexcept
      : data_(data.data()), size_(data.size()) {}
  bool bytes(void* out, std::size_t count) noexcept {
    if (count > size_ - cursor_) return false;
    std::memcpy(out, data_ + cursor_, count);
    cursor_ += count;
    return true;
  }
  bool u8(std::uint8_t& value) noexcept {
    if (cursor_ == size_) return false;
    value = data_[cursor_++];
    return true;
  }
  bool u16(std::uint16_t& value) noexcept {
    std::uint64_t wide = 0U;
    if (!integer(2U, wide)) return false;
    value = static_cast<std::uint16_t>(wide);
    return true;
  }
  bool u32(std::uint32_t& value) noexcept {
    std::uint64_t wide = 0U;
    if (!integer(4U, wide)) return false;
    value = static_cast<std::uint32_t>(wide);
    return true;
  }
  bool i32(std::int32_t& value) noexcept {
    std::uint32_t bits = 0U;
    if (!u32(bits)) return false;
    value = static_cast<std::int32_t>(bits);
    return true;
  }
  bool u64(std::uint64_t& value) noexcept { return integer(8U, value); }
  bool real(double& value) noexcept {
    std::uint64_t bits = 0U;
    if (!u64(bits)) return false;
    std::memcpy(&value, &bits, sizeof(value));
    return true;
  }
  bool int3(Int3& value) noexcept {
    return i32(value.x) && i32(value.y) && i32(value.z);
  }
  std::size_t remaining() const noexcept { return size_ - cursor_; }

 private:
  bool integer(std::size_t bytes, std::uint64_t& value) noexcept {
    if (bytes > size_ - cursor_) return false;
    value = 0U;
    for (std::size_t index = 0U; index < bytes; ++index)
      value |= static_cast<std::uint64_t>(data_[cursor_++]) << (8U * index);
    return true;
  }
  const std::uint8_t* data_{};
  std::size_t size_{};
  std::size_t cursor_{};
};

struct FieldMeta {
  RestartFieldRole role{RestartFieldRole::velocity};
  FieldId field{};
  std::uint8_t components{};
};
constexpr std::size_t kCommonMetadataBound = 128U * sizeof(FieldMeta);

struct RankRecord {
  Int3 begin{};
  Int3 cells{};
  std::uint64_t bytes{};
  std::uint64_t hash{};
};

struct Manifest {
  std::uint32_t format_version{kLegacyFormatVersion};
  std::uint32_t rank_count{};
  Int3 global_cells{};
  PlanFingerprint plan{};
  PlanFingerprint schema{};
  PlanFingerprint geometry{};
  double time{};
  double dt{};
  double pressure_reference{};
  std::uint64_t step{};
  std::uint64_t controller_state{};
  double previous_pressure_reference{};
  double closed_mass_target{};
  RevisionToken final_mass_flux_revision{};
  RevisionToken previous_mass_flux_revision{};
  PlanFingerprint method_history_signature{};
  std::vector<FieldMeta> fields;
  std::vector<FieldMeta> rate_fields;
  std::vector<RankRecord> ranks;
  PlanFingerprint cell_record_identity{};
  std::uint32_t cell_record_bytes{};
};

struct RankBlock {
  std::uint32_t rank_count{};
  std::uint32_t rank{};
  Manifest common;
  MeshPatch patch{};
  std::vector<RestartImageField> fields;
  std::vector<RestartImageField> previous_fields;
  std::vector<RestartImageField> accepted_rate_fields;
  std::vector<RestartImageField> previous_rate_fields;
  std::array<std::vector<double>, 3U> flux;
  std::array<std::vector<double>, 3U> previous_flux;
  std::vector<std::uint8_t> cell_records;
  std::vector<std::uint32_t> cell_record_lengths;
  std::vector<std::size_t> cell_record_offsets;
};

struct BulkBytes {
  std::size_t value{};
  bool add(std::size_t count, std::size_t width = 1U) noexcept {
    std::size_t bytes = 0U;
    if (!checked_multiply(count, width, bytes) || bytes > SIZE_MAX - value) return false;
    value += bytes;
    return true;
  }
};

bool retained_image_bytes(const RestartImage& image, std::size_t& out) noexcept {
  BulkBytes bytes;
  if (!bytes.add(image.cell_records.capacity()) ||
      !bytes.add(image.cell_record_lengths.capacity(), sizeof(std::uint32_t)))
    return false;
  for (const auto* catalog : {&image.fields, &image.previous_fields,
                              &image.accepted_rate_fields, &image.previous_rate_fields}) {
    if (!bytes.add(catalog->capacity(), sizeof(RestartImageField))) return false;
    for (const auto& field : *catalog)
      if (!bytes.add(field.values.capacity(), sizeof(double))) return false;
  }
  for (const auto* flux : {&image.final_mass_flux, &image.previous_mass_flux})
    for (const auto& axis : *flux)
      if (!bytes.add(axis.capacity(), sizeof(double))) return false;
  out = bytes.value;
  return true;
}

// A source block owns only authoritative faces; a target image also stores its
// lower/upper partition interface faces. Count physical allocations, not views.
struct ReadLayout {
  std::size_t owned{}, encoded{}, coverage{};
};

bool read_layout(const Manifest& manifest, const MeshPatch& patch,
                 bool target, ReadLayout& out) noexcept {
  std::size_t cells = 0U;
  if (!cell_count(patch.cells, cells)) return false;
  const bool exact = manifest.format_version >= kExactHistoryFormatVersion;
  const std::size_t layers = exact ? 2U : 1U;
  BulkBytes values, descriptors, vectors, coverage;
  if (!coverage.add(cells)) return false;
  const auto fields = [&](const std::vector<FieldMeta>& catalog, std::size_t copies) {
    for (const auto& field : catalog) {
      std::size_t count = 0U;
      if (!checked_multiply(cells, field.components, count) ||
          !values.add(count, sizeof(double) * copies)) return false;
    }
    return descriptors.add(catalog.size(), sizeof(RestartImageField) * copies) &&
           vectors.add(catalog.size(), copies);
  };
  if (!fields(manifest.fields, layers) ||
      (exact && !fields(manifest.rate_fields, 2U))) return false;
  for (int axis = 0; axis < 3; ++axis) {
    Int3 extents = patch.cells;
    auto& extent = axis == 0 ? extents.x : (axis == 1 ? extents.y : extents.z);
    const int begin = axis == 0 ? patch.begin.x : (axis == 1 ? patch.begin.y : patch.begin.z);
    const int global = axis == 0 ? manifest.global_cells.x :
                       (axis == 1 ? manifest.global_cells.y : manifest.global_cells.z);
    if (target || begin + extent == global) {
      if (extent == INT_MAX) return false;
      ++extent;
    }
    std::size_t faces = 0U;
    if (!cell_count(extents, faces) || !values.add(faces, sizeof(double) * layers) ||
        !vectors.add(layers) || !coverage.add(faces)) return false;
  }
  if (manifest.format_version >= kCellRecordFormatVersion &&
      (!values.add(cells, manifest.cell_record_bytes) || !vectors.add(1U)))
    return false;
  const bool variable =
      manifest.format_version == kVariableCellRecordFormatVersion;
  if (variable && !values.add(cells, sizeof(std::uint32_t)))
    return false;
  BulkBytes owned{values.value};
  if (!owned.add(descriptors.value)) return false;
  if (variable && !(target ? coverage.add(cells + 1, sizeof(std::size_t))
                           : owned.add(cells + 1, sizeof(std::size_t))))
    return false;
  // Common header: 80 + four bytes per field, optional 36 + rates, optional signature.
  BulkBytes encoded{52U + 80U + manifest.fields.size() * 4U};
  if ((exact &&
       (!encoded.add(36U) || !encoded.add(manifest.rate_fields.size(), 4U))) ||
      (manifest.format_version >= kSignedHistoryFormatVersion &&
       !encoded.add(8U)) ||
      (manifest.format_version >= kCellRecordFormatVersion &&
       !encoded.add(12U)) ||
      !encoded.add(vectors.value, 8U) || !encoded.add(values.value))
    return false;
  out = {owned.value, encoded.value, target ? coverage.value : 0U};
  return true;
}

Status plan_read_bulk(const Manifest &manifest, const RestartExpected &expected,
                      std::size_t retained, RestartReadLimits limits,
                      RestartReadReport &report,
                      std::size_t variable_target_bytes = 0U) noexcept {
  ReadLayout target;
  if (!read_layout(manifest, expected.target_patch, true, target))
    return {StatusCode::allocation_failure, kRestartReadBudget};
  if (variable_target_bytes > SIZE_MAX - target.owned)
    return {StatusCode::allocation_failure, kRestartReadBudget};
  target.owned += variable_target_bytes;
  std::size_t largest_block = 0U;
  for (const auto& record : manifest.ranks) {
    ReadLayout source;
    if (!read_layout(manifest, MeshPatch{record.begin, record.cells, {}, {}},
                     false, source))
      return {StatusCode::io_failure, kRestartManifest};
    if (manifest.format_version == kVariableCellRecordFormatVersion) {
      if (record.bytes < source.encoded || record.bytes > SIZE_MAX ||
          record.bytes - source.encoded > SIZE_MAX - source.owned)
        return {StatusCode::io_failure, kRestartManifest};
      source.owned += std::size_t(record.bytes) - source.encoded;
      source.encoded = std::size_t(record.bytes);
    } else if (record.bytes != source.encoded)
      return {StatusCode::io_failure, kRestartManifest};
    BulkBytes block{source.encoded};
    // Includes bytes + decoded values + catalogs at the same time. A malformed
    // common header may use all 64+64 metadata slots before it is compared.
    if (!block.add(source.owned) || !block.add(kCommonMetadataBound))
      return {StatusCode::allocation_failure, kRestartReadBudget};
    largest_block = std::max(largest_block, block.value);
  }
  BulkBytes peak{retained};
  if (!peak.add(manifest.fields.capacity(), sizeof(FieldMeta)) ||
      !peak.add(manifest.rate_fields.capacity(), sizeof(FieldMeta)) ||
      !peak.add(manifest.ranks.capacity(), sizeof(RankRecord)) ||
      !peak.add(target.owned) || !peak.add(target.coverage) || !peak.add(largest_block))
    return {StatusCode::allocation_failure, kRestartReadBudget};
  report.new_image_bytes = target.owned;
  report.peak_bulk_bytes = std::max(report.peak_bulk_bytes, peak.value);
  return report.peak_bulk_bytes <= limits.maximum_bulk_bytes ? Status{} :
      Status{StatusCode::allocation_failure, kRestartReadBudget};
}

Status collective_status(MPI_Comm communicator, Status local) noexcept {
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0)
    return {StatusCode::mpi_failure, kRestartCollective};
  const int candidate = local ? size : rank;
  int failing = size;
  if (MPI_Allreduce(&candidate, &failing, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRestartCollective};
  if (failing == size) return {};
  std::uint64_t wire = 0U;
  if (rank == failing)
    wire = (static_cast<std::uint64_t>(local.code) << 32U) | local.detail;
  if (MPI_Bcast(&wire, 1, MPI_UINT64_T, failing, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRestartCollective};
  return {static_cast<StatusCode>(wire >> 32U),
          static_cast<std::uint32_t>(wire)};
}

Status consensus_u64(MPI_Comm communicator, std::uint64_t value) noexcept {
  std::uint64_t minimum = value;
  std::uint64_t maximum = value;
  if (MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRestartCollective};
  return minimum == maximum
             ? Status{}
             : Status{StatusCode::invalid_plan, kRestartMismatch};
}

template <class LocalWork>
Status restart_local_stage(MPI_Comm communicator, LocalWork&& work,
                           IoFailureContext* failure = nullptr) noexcept {
  Status status;
  try {
    status = work();
  } catch (const std::bad_alloc&) {
    status = {StatusCode::allocation_failure, kRestartInput};
  } catch (...) {
    status = {StatusCode::io_failure, kRestartInput};
  }
  return failure != nullptr
             ? detail::output_collective_status(communicator, status, failure)
             : collective_status(communicator, status);
}

bool write_file_sync(const fs::path& path,
                     const std::vector<std::uint8_t>& bytes,
                     IoFailureContext* failure = nullptr) noexcept {
  return detail::output_write_file(path, bytes.data(), bytes.size(),
                                   detail::OutputFileMode::exclusive, failure);
}

bool read_file(const fs::path& path, std::vector<std::uint8_t>& bytes,
               std::size_t maximum_bytes, std::uint64_t expected_bytes,
               IoFailureContext* failure = nullptr) {
  int descriptor;
  do {
    // O_NONBLOCK lets fstat reject a FIFO without waiting for an external writer.
    descriptor = ::open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NONBLOCK);
  } while (descriptor < 0 && errno == EINTR);
  if (descriptor < 0) {
    detail::output_record_failure(failure, IoFailureOperation::open, errno, path);
    return false;
  }
  struct stat info {};
  int stated;
  do {
    stated = ::fstat(descriptor, &info);
  } while (stated != 0 && errno == EINTR);
  const bool representable = stated == 0 && info.st_size >= 0 &&
      static_cast<std::uintmax_t>(info.st_size) <= SIZE_MAX;
  bool okay = representable && S_ISREG(info.st_mode) &&
      static_cast<std::uintmax_t>(info.st_size) <= maximum_bytes &&
      (expected_bytes == 0U || static_cast<std::uintmax_t>(info.st_size) == expected_bytes);
  if (!okay)
    detail::output_record_failure(failure, IoFailureOperation::stat,
        stated != 0 ? errno : (!representable ? EOVERFLOW :
        (!S_ISREG(info.st_mode) ? EINVAL :
         (static_cast<std::uintmax_t>(info.st_size) > maximum_bytes ? EFBIG : EINVAL))), path);
  if (okay) {
    try {
      bytes.resize(static_cast<std::size_t>(info.st_size));
    } catch (...) {
      ::close(descriptor);
      throw;  // The enclosing local stage distinguishes allocation from I/O.
    }
  }
  std::size_t cursor = 0U;
  while (okay && cursor < bytes.size()) {
    const ssize_t count =
        ::read(descriptor, bytes.data() + cursor, bytes.size() - cursor);
    if (count < 0 && errno == EINTR) continue;
    if (count <= 0) {
      detail::output_record_failure(failure, IoFailureOperation::read,
                                    count == 0 ? EIO : errno, path);
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(count);
  }
  if (okay) {
    std::uint8_t extra = 0U;
    ssize_t count;
    do { count = ::read(descriptor, &extra, 1U); } while (count < 0 && errno == EINTR);
    if (count != 0) {
      detail::output_record_failure(failure, IoFailureOperation::read,
                                    count < 0 ? errno : EIO, path);
      okay = false;  // Detect growth after the preallocation fstat as well.
    }
  }
  if (::close(descriptor) != 0) {
    detail::output_record_failure(failure, IoFailureOperation::close, errno, path);
    okay = false;
  }
  return okay;
}

bool sync_directory(const fs::path& path,
                     IoFailureContext* failure = nullptr) noexcept {
  return detail::output_sync_directory(path, failure);
}

std::string rank_name(std::uint32_t rank) {
  std::string number = std::to_string(rank);
  return "rank-" + std::string(8U - std::min<std::size_t>(8U, number.size()), '0') +
         number + ".bin";
}

bool verified_integrity(const std::vector<std::uint8_t>& bytes) noexcept {
  if (bytes.size() < sizeof(std::uint64_t)) return false;
  std::uint64_t stored = 0U;
  for (std::size_t index = 0U; index < sizeof(stored); ++index)
    stored |= static_cast<std::uint64_t>(
                  bytes[bytes.size() - sizeof(stored) + index])
              << (8U * index);
  return stored == hash_bytes(bytes.data(), bytes.size() - sizeof(stored));
}

bool expected_face_extents(ConstFaceFluxView flux, Int3 cells) noexcept {
  return flux.revision != 0U && flux.certificate.valid() &&
         flux.certificate.matches(flux) &&
         same(flux.x.extents, {cells.x + 1, cells.y, cells.z}) &&
         same(flux.y.extents, {cells.x, cells.y + 1, cells.z}) &&
         same(flux.z.extents, {cells.x, cells.y, cells.z + 1}) &&
         flux.x.axis == CartesianAxis::x &&
         flux.y.axis == CartesianAxis::y &&
         flux.z.axis == CartesianAxis::z;
}

bool valid_field_catalog(Span<const RestartFieldView> fields,
                         Int3 cells, bool rates) noexcept {
  if (fields.data == nullptr || fields.size == 0U || fields.size > 64U)
    return false;
  bool velocity = false;
  bool pressure = false;
  bool enthalpy = false;
  bool enthalpy_rate = false;
  for (std::size_t index = 0U; index < fields.size; ++index) {
    const RestartFieldView& field = fields.data[index];
    detail::FieldStorageInterval interval{};
    if (!detail::field_storage_interval(field.values, interval) ||
        !same(field.values.interior, cells) ||
        field.values.revision == 0U) {
      return false;
    }
    for (std::size_t prior = 0U; prior < index; ++prior) {
      if (fields.data[prior].values.field == field.values.field)
        return false;
    }
    switch (field.role) {
      case RestartFieldRole::velocity:
        if (rates || velocity || field.values.components != 3U) return false;
        velocity = true;
        break;
      case RestartFieldRole::pressure_perturbation:
      case RestartFieldRole::pressure_absolute:
        if (rates || pressure || field.values.components != 1U) return false;
        pressure = true;
        break;
      case RestartFieldRole::enthalpy:
        if (rates || enthalpy || field.values.components != 1U) return false;
        enthalpy = true;
        break;
      case RestartFieldRole::independent_species:
      case RestartFieldRole::transported_scalar:
        if (rates || field.values.components != 1U) return false;
        break;
      case RestartFieldRole::stochastic_field:
        if (rates || field.values.components < 2U)
        return false;
        break;
      case RestartFieldRole::stochastic_transport:
        if (rates || field.values.components != 2U)
        return false;
        break;
      case RestartFieldRole::enthalpy_nonadvective_rate:
        if (!rates || enthalpy_rate || field.values.components != 1U)
          return false;
        enthalpy_rate = true;
        break;
      case RestartFieldRole::scalar_nonadvective_rate:
        if (!rates || field.values.components != 1U) return false;
        break;
      default:
        return false;
    }
  }
  return rates ? enthalpy_rate : (velocity && pressure && enthalpy);
}

bool same_field_layout(Span<const RestartFieldView> left,
                       Span<const RestartFieldView> right) noexcept {
  if (left.size != right.size || left.data == nullptr || right.data == nullptr)
    return false;
  for (std::size_t index = 0U; index < left.size; ++index) {
    if (left.data[index].role != right.data[index].role ||
        left.data[index].values.field != right.data[index].values.field ||
        left.data[index].values.components !=
            right.data[index].values.components)
      return false;
  }
  return true;
}

bool has_exact_history(const RestartSnapshot& snapshot) noexcept {
  return snapshot.previous_fields.size != 0U ||
         snapshot.accepted_rate_fields.size != 0U ||
         snapshot.previous_rate_fields.size != 0U ||
         snapshot.previous_mass_flux.revision != 0U ||
         snapshot.previous_pressure_reference != 0.0 ||
         snapshot.closed_mass_target != 0.0 ||
         snapshot.method_history_signature != 0U;
}

std::uint32_t snapshot_format(const RestartSnapshot& snapshot) noexcept {
  if (snapshot.cell_records.identity != 0U)
    return snapshot.cell_records.record_bytes == 0U
               ? kVariableCellRecordFormatVersion
               : kCellRecordFormatVersion;
  return !has_exact_history(snapshot) ? kLegacyFormatVersion
      : snapshot.method_history_signature != 0U ? kSignedHistoryFormatVersion
                                               : kExactHistoryFormatVersion;
}

bool valid_snapshot(const RestartSnapshot& snapshot) noexcept {
  if (!valid_global_patch(snapshot.global_cells, snapshot.patch) ||
      snapshot.plan == 0U || snapshot.schema == 0U ||
      snapshot.geometry == 0U || !std::isfinite(snapshot.time) ||
      !std::isfinite(snapshot.dt) || snapshot.dt <= 0.0 ||
      !std::isfinite(snapshot.pressure_reference) || snapshot.step == 0U ||
      !valid_field_catalog(snapshot.fields, snapshot.patch.cells, false) ||
      !expected_face_extents(snapshot.final_mass_flux, snapshot.patch.cells)) {
    return false;
  }
  const auto &records = snapshot.cell_records;
  if (records.identity != 0U) {
    std::size_t cells = 0U, bytes = 0U;
    if (!has_exact_history(snapshot) ||
        snapshot.method_history_signature == 0U ||
        !cell_count(snapshot.patch.cells, cells))
      return false;
    if (records.record_bytes == 0U) {
      if (records.variable_cell_bytes.size != cells ||
          !records.variable_cell_bytes.data)
          return false;
      for (std::size_t i = 0; i < cells; ++i) {
          if (records.variable_cell_bytes.data[i] > SIZE_MAX - bytes)
            return false;
          bytes += records.variable_cell_bytes.data[i];
      }
    } else if (records.variable_cell_bytes.size != 0U ||
               records.variable_cell_bytes.data != nullptr ||
               !checked_multiply(cells, records.record_bytes, bytes))
      return false;
    if (records.values.size != bytes || (bytes != 0U && !records.values.data))
      return false;
  } else if (records.record_bytes != 0U || records.values.size != 0U ||
             records.values.data != nullptr ||
             records.variable_cell_bytes.size != 0U ||
             records.variable_cell_bytes.data != nullptr)
    return false;
  if (!has_exact_history(snapshot)) return true;
  return snapshot.controller_state != 0U &&
         valid_field_catalog(snapshot.previous_fields,
                             snapshot.patch.cells, false) &&
         same_field_layout(snapshot.fields, snapshot.previous_fields) &&
         valid_field_catalog(snapshot.accepted_rate_fields,
                             snapshot.patch.cells, true) &&
         valid_field_catalog(snapshot.previous_rate_fields,
                             snapshot.patch.cells, true) &&
         same_field_layout(snapshot.accepted_rate_fields,
                           snapshot.previous_rate_fields) &&
         expected_face_extents(snapshot.previous_mass_flux,
                               snapshot.patch.cells) &&
         snapshot.previous_mass_flux.revision <
             snapshot.final_mass_flux.revision &&
         snapshot.final_mass_flux.revision !=
             std::numeric_limits<RevisionToken>::max() &&
         std::isfinite(snapshot.previous_pressure_reference) &&
         std::isfinite(snapshot.closed_mass_target) &&
         snapshot.closed_mass_target > 0.0;
}

std::uint64_t snapshot_signature(const RestartSnapshot& snapshot) {
  Encoder encoder;
  const std::uint32_t version = snapshot_format(snapshot);
  encoder.u32(version);
  encoder.int3(snapshot.global_cells);
  encoder.u64(snapshot.plan);
  encoder.u64(snapshot.schema);
  encoder.u64(snapshot.geometry);
  encoder.real(snapshot.time);
  encoder.real(snapshot.dt);
  encoder.real(snapshot.pressure_reference);
  encoder.u64(snapshot.step);
  encoder.u64(snapshot.controller_state);
  encoder.u32(static_cast<std::uint32_t>(snapshot.fields.size));
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView field = snapshot.fields.data[index];
    encoder.u8(static_cast<std::uint8_t>(field.role));
    encoder.u16(field.values.field);
    encoder.u8(field.values.components);
  }
  if (version >= kExactHistoryFormatVersion) {
    encoder.real(snapshot.previous_pressure_reference);
    encoder.real(snapshot.closed_mass_target);
    encoder.u64(snapshot.final_mass_flux.revision);
    encoder.u64(snapshot.previous_mass_flux.revision);
    encoder.u32(
        static_cast<std::uint32_t>(snapshot.accepted_rate_fields.size));
    for (std::size_t index = 0U;
         index < snapshot.accepted_rate_fields.size; ++index) {
      const RestartFieldView field =
          snapshot.accepted_rate_fields.data[index];
      encoder.u8(static_cast<std::uint8_t>(field.role));
      encoder.u16(field.values.field);
      encoder.u8(field.values.components);
    }
  }
  if (version >= kSignedHistoryFormatVersion)
    encoder.u64(snapshot.method_history_signature);
  if (version >= kCellRecordFormatVersion) {
    encoder.u64(snapshot.cell_records.identity);
    encoder.u32(snapshot.cell_records.record_bytes);
  }
  return hash_bytes(encoder.data().data(), encoder.data().size());
}

void encode_common(Encoder& encoder, const RestartSnapshot& snapshot,
                   std::uint32_t version) {
  encoder.int3(snapshot.global_cells);
  encoder.u64(snapshot.plan);
  encoder.u64(snapshot.schema);
  encoder.u64(snapshot.geometry);
  encoder.real(snapshot.time);
  encoder.real(snapshot.dt);
  encoder.real(snapshot.pressure_reference);
  encoder.u64(snapshot.step);
  encoder.u64(snapshot.controller_state);
  encoder.u32(static_cast<std::uint32_t>(snapshot.fields.size));
  for (std::size_t index = 0U; index < snapshot.fields.size; ++index) {
    const RestartFieldView field = snapshot.fields.data[index];
    encoder.u8(static_cast<std::uint8_t>(field.role));
    encoder.u16(field.values.field);
    encoder.u8(field.values.components);
  }
  if (version >= kExactHistoryFormatVersion) {
    encoder.real(snapshot.previous_pressure_reference);
    encoder.real(snapshot.closed_mass_target);
    encoder.u64(snapshot.final_mass_flux.revision);
    encoder.u64(snapshot.previous_mass_flux.revision);
    encoder.u32(
        static_cast<std::uint32_t>(snapshot.accepted_rate_fields.size));
    for (std::size_t index = 0U;
         index < snapshot.accepted_rate_fields.size; ++index) {
      const RestartFieldView field =
          snapshot.accepted_rate_fields.data[index];
      encoder.u8(static_cast<std::uint8_t>(field.role));
      encoder.u16(field.values.field);
      encoder.u8(field.values.components);
    }
  }
  if (version >= kSignedHistoryFormatVersion)
    encoder.u64(snapshot.method_history_signature);
  if (version >= kCellRecordFormatVersion) {
    encoder.u64(snapshot.cell_records.identity);
    encoder.u32(snapshot.cell_records.record_bytes);
  }
}

bool decode_common(Decoder& decoder, std::uint32_t version,
                   Manifest& manifest) {
  std::uint32_t field_count = 0U;
  manifest.format_version = version;
  if (!decoder.int3(manifest.global_cells) || !decoder.u64(manifest.plan) ||
      !decoder.u64(manifest.schema) || !decoder.u64(manifest.geometry) ||
      !decoder.real(manifest.time) || !decoder.real(manifest.dt) ||
      !decoder.real(manifest.pressure_reference) ||
      !decoder.u64(manifest.step) ||
      !decoder.u64(manifest.controller_state) ||
      !decoder.u32(field_count) || field_count < 3U || field_count > 64U ||
      manifest.plan == 0U || manifest.schema == 0U ||
      manifest.geometry == 0U || !std::isfinite(manifest.time) ||
      !std::isfinite(manifest.dt) || manifest.dt <= 0.0 ||
      !std::isfinite(manifest.pressure_reference) || manifest.step == 0U) {
    return false;
  }
  manifest.fields.resize(field_count);
  for (FieldMeta& field : manifest.fields) {
    std::uint8_t role = 0U;
    if (!decoder.u8(role) || !decoder.u16(field.field) ||
        !decoder.u8(field.components) ||
        (role >
             static_cast<std::uint8_t>(RestartFieldRole::transported_scalar) &&
         role !=
             static_cast<std::uint8_t>(RestartFieldRole::stochastic_field) &&
         role != static_cast<std::uint8_t>(
                     RestartFieldRole::stochastic_transport)) ||
        field.components == 0U) {
      return false;
    }
    field.role = static_cast<RestartFieldRole>(role);
  }
  if (version >= kExactHistoryFormatVersion) {
    std::uint32_t rate_count = 0U;
    if (!decoder.real(manifest.previous_pressure_reference) ||
        !decoder.real(manifest.closed_mass_target) ||
        !decoder.u64(manifest.final_mass_flux_revision) ||
        !decoder.u64(manifest.previous_mass_flux_revision) ||
        !decoder.u32(rate_count) || rate_count == 0U || rate_count > 64U ||
        manifest.controller_state == 0U ||
        !std::isfinite(manifest.previous_pressure_reference) ||
        !std::isfinite(manifest.closed_mass_target) ||
        !(manifest.closed_mass_target > 0.0) ||
        manifest.previous_mass_flux_revision == 0U ||
        manifest.previous_mass_flux_revision >=
            manifest.final_mass_flux_revision ||
        manifest.final_mass_flux_revision ==
            std::numeric_limits<RevisionToken>::max()) {
      return false;
    }
    manifest.rate_fields.resize(rate_count);
    for (FieldMeta& field : manifest.rate_fields) {
      std::uint8_t role = 0U;
      if (!decoder.u8(role) || !decoder.u16(field.field) ||
          !decoder.u8(field.components) || field.components != 1U ||
          role < static_cast<std::uint8_t>(
                     RestartFieldRole::enthalpy_nonadvective_rate) ||
          role > static_cast<std::uint8_t>(
                     RestartFieldRole::scalar_nonadvective_rate)) {
        return false;
      }
      field.role = static_cast<RestartFieldRole>(role);
    }
  } else if (version != kLegacyFormatVersion) {
    return false;
  }
  if (version >= kSignedHistoryFormatVersion &&
      (!decoder.u64(manifest.method_history_signature) ||
       manifest.method_history_signature == 0U))
    return false;
  if (version >= kCellRecordFormatVersion &&
      (!decoder.u64(manifest.cell_record_identity) ||
       !decoder.u32(manifest.cell_record_bytes) ||
       manifest.cell_record_identity == 0U ||
       (version == kVariableCellRecordFormatVersion
            ? manifest.cell_record_bytes != 0U
            : manifest.cell_record_bytes == 0U)))
    return false;
  return true;
}

bool rank_block_size(const RestartSnapshot& snapshot, std::size_t& bytes) {
  Encoder header;
  encode_common(header, snapshot, snapshot_format(snapshot));
  // Magic, version/size/rank, patch, and trailing integrity.
  bytes = header.data().size() + 8U + 12U + 24U + 8U;
  const auto add_values = [&](std::size_t count) {
    std::size_t payload = 0U;
    if (!checked_multiply(count, sizeof(double), payload) ||
        bytes > SIZE_MAX - 8U || payload > SIZE_MAX - bytes - 8U) return false;
    bytes += 8U + payload;
    return true;
  };
  std::size_t cells = 0U;
  if (!cell_count(snapshot.patch.cells, cells)) return false;
  const auto fields = [&](Span<const RestartFieldView> views) {
    for (std::size_t i = 0U; i < views.size; ++i) {
      const auto& view = views.data[i];
      std::size_t values = 0U;
      if (!checked_multiply(cells, view.values.components, values) ||
          !add_values(values)) return false;
    }
    return true;
  };
  if (!fields(snapshot.fields)) return false;
  const bool exact = has_exact_history(snapshot);
  if (exact && (!fields(snapshot.previous_fields) ||
                !fields(snapshot.accepted_rate_fields) ||
                !fields(snapshot.previous_rate_fields))) return false;
  for (int history = 0; history < (exact ? 2 : 1); ++history) {
    for (int axis = 0; axis < 3; ++axis) {
      Int3 owned = snapshot.patch.cells;
      if (axis == 0 && snapshot.patch.begin.x + owned.x == snapshot.global_cells.x) ++owned.x;
      if (axis == 1 && snapshot.patch.begin.y + owned.y == snapshot.global_cells.y) ++owned.y;
      if (axis == 2 && snapshot.patch.begin.z + owned.z == snapshot.global_cells.z) ++owned.z;
      std::size_t count = 0U;
      if (!cell_count(owned, count) || !add_values(count)) return false;
    }
  }
  if (snapshot.cell_records.identity != 0U) {
    if (bytes > SIZE_MAX - 8U ||
        snapshot.cell_records.values.size > SIZE_MAX - bytes - 8U)
      return false;
    bytes += 8U + snapshot.cell_records.values.size;
    if (snapshot_format(snapshot) == kVariableCellRecordFormatVersion) {
      std::size_t lengths = 0;
      if (!checked_multiply(cells, sizeof(std::uint32_t), lengths) ||
          lengths > SIZE_MAX - bytes)
        return false;
      bytes += lengths;
    }
  }
  return true;
}

Status encode_rank_block(const RestartSnapshot& snapshot, int size, int rank,
                         std::vector<std::uint8_t>& out,
                         std::size_t maximum_bytes) {
  std::size_t bytes = 0U;
  if (!rank_block_size(snapshot, bytes))
    return {StatusCode::invalid_plan, kRestartInput};
  if (maximum_bytes != 0U && bytes > maximum_bytes)
    return {StatusCode::allocation_failure, kRestartRankFile};
  Encoder encoder(bytes);
  const std::uint32_t version = snapshot_format(snapshot);
  encoder.bytes(kRankMagic.data(), kRankMagic.size());
  encoder.u32(version);
  encoder.u32(static_cast<std::uint32_t>(size));
  encoder.u32(static_cast<std::uint32_t>(rank));
  encode_common(encoder, snapshot, version);
  encoder.int3(snapshot.patch.begin);
  encoder.int3(snapshot.patch.cells);
  std::size_t cells = 0U;
  if (!cell_count(snapshot.patch.cells, cells))
    return {StatusCode::invalid_plan, kRestartInput};
  const auto encode_fields = [&](Span<const RestartFieldView> fields) {
    for (std::size_t field_index = 0U; field_index < fields.size;
         ++field_index) {
      const ConstFieldView view = fields.data[field_index].values;
      std::size_t values = 0U;
      if (!checked_multiply(cells, view.components, values))
        return Status{StatusCode::invalid_plan, kRestartInput};
      encoder.u64(values);
      for (std::int32_t z = 0; z < view.interior.z; ++z)
        for (std::int32_t y = 0; y < view.interior.y; ++y)
          for (std::int32_t x = 0; x < view.interior.x; ++x)
            for (std::uint8_t component = 0U;
                 component < view.components; ++component) {
              const double value = view.unchecked({x, y, z}, component);
              if (!std::isfinite(value))
                return Status{StatusCode::numerical_failure, kRestartInput};
              encoder.real(value);
            }
    }
    return Status{};
  };
  const auto encode_flux = [&](ConstFaceFluxView flux) {
    const std::array<ConstFaceFieldView, 3U> faces{
        flux.x, flux.y, flux.z};
    for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
      Int3 owned = snapshot.patch.cells;
      const std::int32_t patch_end =
          axis == 0U ? snapshot.patch.begin.x + snapshot.patch.cells.x
                     : (axis == 1U
                            ? snapshot.patch.begin.y + snapshot.patch.cells.y
                            : snapshot.patch.begin.z + snapshot.patch.cells.z);
      const std::int32_t global_end =
          axis == 0U ? snapshot.global_cells.x
                     : (axis == 1U ? snapshot.global_cells.y
                                   : snapshot.global_cells.z);
      if (patch_end == global_end) {
        if (axis == 0U)
          ++owned.x;
        else if (axis == 1U)
          ++owned.y;
        else
          ++owned.z;
      }
      std::size_t values = 0U;
      if (!cell_count(owned, values))
        return Status{StatusCode::invalid_plan, kRestartInput};
      encoder.u64(values);
      for (std::int32_t z = 0; z < owned.z; ++z)
        for (std::int32_t y = 0; y < owned.y; ++y)
          for (std::int32_t x = 0; x < owned.x; ++x) {
            const double value = faces[axis].unchecked({x, y, z});
            if (!std::isfinite(value))
              return Status{StatusCode::numerical_failure, kRestartInput};
            encoder.real(value);
          }
    }
    return Status{};
  };
  Status status = encode_fields(snapshot.fields);
  if (status && version >= kExactHistoryFormatVersion)
    status = encode_fields(snapshot.previous_fields);
  if (status && version >= kExactHistoryFormatVersion)
    status = encode_fields(snapshot.accepted_rate_fields);
  if (status && version >= kExactHistoryFormatVersion)
    status = encode_fields(snapshot.previous_rate_fields);
  if (status) status = encode_flux(snapshot.final_mass_flux);
  if (status && version >= kExactHistoryFormatVersion)
    status = encode_flux(snapshot.previous_mass_flux);
  if (status && version >= kCellRecordFormatVersion) {
    encoder.u64(snapshot.cell_records.values.size);
    if (snapshot_format(snapshot) == kVariableCellRecordFormatVersion)
      for (std::size_t i = 0;
           i < snapshot.cell_records.variable_cell_bytes.size; ++i)
        encoder.u32(snapshot.cell_records.variable_cell_bytes.data[i]);
    if (snapshot.cell_records.values.size != 0U)
      encoder.bytes(snapshot.cell_records.values.data,
                    snapshot.cell_records.values.size);
  }
  if (!status) return status;
  encoder.append_integrity();
  if (encoder.data().size() != bytes)
    return {StatusCode::invalid_plan, kRestartRankFile};
  out = encoder.take();
  return {};
}

Status encode_manifest(const RestartSnapshot& snapshot, int size,
                       const std::vector<RankRecord>& records,
                       std::vector<std::uint8_t>& out) {
  if (records.size() != static_cast<std::size_t>(size))
    return {StatusCode::invalid_plan, kRestartManifest};
  const std::uint32_t version = snapshot_format(snapshot);
  Encoder common;
  encode_common(common, snapshot, version);
  if (records.size() > (SIZE_MAX - common.data().size() - 24U) / 40U)
    return {StatusCode::invalid_plan, kRestartManifest};
  Encoder encoder(common.data().size() + 24U + records.size() * 40U);
  encoder.bytes(kManifestMagic.data(), kManifestMagic.size());
  encoder.u32(version);
  encoder.u32(static_cast<std::uint32_t>(size));
  encode_common(encoder, snapshot, version);
  for (const RankRecord& record : records) {
    encoder.int3(record.begin);
    encoder.int3(record.cells);
    encoder.u64(record.bytes);
    encoder.u64(record.hash);
  }
  encoder.append_integrity();
  out = encoder.take();
  return {};
}

Status parse_manifest(const std::vector<std::uint8_t>& bytes,
                      Manifest& out, std::uint32_t maximum_ranks) noexcept {
  if (!verified_integrity(bytes))
    return {StatusCode::io_failure, kRestartIntegrity};
  try {
    Decoder decoder(bytes);
    std::array<char, 8U> magic{};
    std::uint32_t version = 0U;
    Manifest candidate;
    if (!decoder.bytes(magic.data(), magic.size()) || magic != kManifestMagic ||
        !decoder.u32(version) ||
        (version != kLegacyFormatVersion &&
         version != kExactHistoryFormatVersion &&
         version != kSignedHistoryFormatVersion &&
         version != kCellRecordFormatVersion &&
         version != kVariableCellRecordFormatVersion) ||
        !decoder.u32(candidate.rank_count) || candidate.rank_count == 0U ||
        candidate.rank_count > maximum_ranks ||
        !decode_common(decoder, version, candidate) ||
        !valid_global_patch(
            candidate.global_cells,
            MeshPatch{{0, 0, 0}, candidate.global_cells, {}, {}})) {
      return {StatusCode::io_failure, kRestartManifest};
    }
    // Each record is exactly 40 encoded bytes, followed by one integrity word.
    // Reject malicious counts before allocating the rank directory.
    if (decoder.remaining() < 8U ||
        candidate.rank_count != (decoder.remaining() - 8U) / 40U ||
        (decoder.remaining() - 8U) % 40U != 0U)
      return {StatusCode::io_failure, kRestartManifest};
    candidate.ranks.resize(candidate.rank_count);
    for (RankRecord& record : candidate.ranks) {
      if (!decoder.int3(record.begin) || !decoder.int3(record.cells) ||
          !decoder.u64(record.bytes) || !decoder.u64(record.hash) ||
          record.bytes < 16U || record.hash == 0U ||
          !valid_global_patch(candidate.global_cells,
                              MeshPatch{record.begin, record.cells, {}, {}})) {
        return {StatusCode::io_failure, kRestartManifest};
      }
    }
    std::uint64_t integrity = 0U;
    if (!decoder.u64(integrity) || decoder.remaining() != 0U) {
      return {StatusCode::io_failure, kRestartManifest};
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kRestartManifest};
  } catch (...) {
    return {StatusCode::io_failure, kRestartManifest};
  }
}

bool same_common(const Manifest& left, const Manifest& right) noexcept;

Status parse_rank_block(const std::vector<std::uint8_t>& bytes,
                        RankBlock& out, const Manifest& expected,
                        std::uint32_t source) noexcept {
  if (!verified_integrity(bytes))
    return {StatusCode::io_failure, kRestartIntegrity};
  try {
    Decoder decoder(bytes);
    std::array<char, 8U> magic{};
    std::uint32_t version = 0U;
    RankBlock candidate;
    if (!decoder.bytes(magic.data(), magic.size()) || magic != kRankMagic ||
        !decoder.u32(version) ||
        (version != kLegacyFormatVersion &&
         version != kExactHistoryFormatVersion &&
         version != kSignedHistoryFormatVersion &&
         version != kCellRecordFormatVersion &&
         version != kVariableCellRecordFormatVersion) ||
        !decoder.u32(candidate.rank_count) || candidate.rank_count == 0U ||
        !decoder.u32(candidate.rank) ||
        candidate.rank >= candidate.rank_count ||
        !decode_common(decoder, version, candidate.common) ||
        !decoder.int3(candidate.patch.begin) ||
        !decoder.int3(candidate.patch.cells) ||
        !valid_global_patch(candidate.common.global_cells, candidate.patch)) {
      return {StatusCode::io_failure, kRestartRankFile};
    }
    const auto& record = expected.ranks[source];
    if (candidate.rank != source || candidate.rank_count != expected.rank_count ||
        !same(candidate.patch.begin, record.begin) || !same(candidate.patch.cells, record.cells) ||
        !same_common(candidate.common, expected))
      return {StatusCode::io_failure, kRestartMismatch};
    std::size_t cells = 0U;
    if (!cell_count(candidate.patch.cells, cells))
      return {StatusCode::io_failure, kRestartRankFile};
    const auto decode_fields = [&](const std::vector<FieldMeta>& metadata,
                                   std::vector<RestartImageField>& fields) {
      fields.resize(metadata.size());
      for (std::size_t field_index = 0U; field_index < fields.size();
           ++field_index) {
        const FieldMeta meta = metadata[field_index];
        RestartImageField& field = fields[field_index];
        field.role = meta.role;
        field.field = meta.field;
        field.components = meta.components;
        std::size_t expected_values = 0U;
        std::uint64_t stored_values = 0U;
        if (!checked_multiply(cells, meta.components, expected_values) ||
            !decoder.u64(stored_values) || stored_values != expected_values ||
            expected_values > decoder.remaining() / sizeof(double))
          return false;
        field.values.resize(expected_values);
        for (double& value : field.values)
          if (!decoder.real(value) || !std::isfinite(value)) return false;
      }
      return true;
    };
    const auto decode_flux = [&](std::array<std::vector<double>, 3U>& flux) {
      for (std::size_t axis = 0U; axis < flux.size(); ++axis) {
        Int3 owned = candidate.patch.cells;
        const std::int32_t patch_end =
            axis == 0U
                ? candidate.patch.begin.x + candidate.patch.cells.x
                : (axis == 1U
                       ? candidate.patch.begin.y + candidate.patch.cells.y
                       : candidate.patch.begin.z + candidate.patch.cells.z);
        const std::int32_t global_end =
            axis == 0U ? candidate.common.global_cells.x
                       : (axis == 1U ? candidate.common.global_cells.y
                                     : candidate.common.global_cells.z);
        if (patch_end == global_end) {
          if (axis == 0U)
            ++owned.x;
          else if (axis == 1U)
            ++owned.y;
          else
            ++owned.z;
        }
        std::size_t expected_values = 0U;
        std::uint64_t stored_values = 0U;
        if (!cell_count(owned, expected_values) ||
            !decoder.u64(stored_values) || stored_values != expected_values ||
            expected_values > decoder.remaining() / sizeof(double))
          return false;
        flux[axis].resize(expected_values);
        for (double& value : flux[axis])
          if (!decoder.real(value) || !std::isfinite(value)) return false;
      }
      return true;
    };
    if (!decode_fields(candidate.common.fields, candidate.fields))
      return {StatusCode::io_failure, kRestartRankFile};
    if (version >= kExactHistoryFormatVersion &&
        (!decode_fields(candidate.common.fields, candidate.previous_fields) ||
         !decode_fields(candidate.common.rate_fields,
                        candidate.accepted_rate_fields) ||
         !decode_fields(candidate.common.rate_fields,
                        candidate.previous_rate_fields)))
      return {StatusCode::io_failure, kRestartRankFile};
    if (!decode_flux(candidate.flux) ||
        (version >= kExactHistoryFormatVersion &&
         !decode_flux(candidate.previous_flux)))
      return {StatusCode::io_failure, kRestartRankFile};
    if (version >= kCellRecordFormatVersion) {
      std::size_t expected_bytes = 0U;
      std::uint64_t encoded_bytes = 0U;
      if (!decoder.u64(encoded_bytes))
        return {StatusCode::io_failure, kRestartRankFile};
      if (version == kVariableCellRecordFormatVersion) {
        if (cells > decoder.remaining() / sizeof(std::uint32_t))
          return {StatusCode::io_failure, kRestartRankFile};
        candidate.cell_record_lengths.resize(cells);
        candidate.cell_record_offsets.resize(cells + 1);
        for (std::size_t i = 0; i < cells; ++i) {
          auto &n = candidate.cell_record_lengths[i];
          if (!decoder.u32(n) || n > SIZE_MAX - expected_bytes)
            return {StatusCode::io_failure, kRestartRankFile};
          candidate.cell_record_offsets[i] = expected_bytes;
          expected_bytes += n;
        }
        candidate.cell_record_offsets[cells] = expected_bytes;
      } else if (!checked_multiply(cells, candidate.common.cell_record_bytes,
                                   expected_bytes))
        return {StatusCode::io_failure, kRestartRankFile};
      if (encoded_bytes != expected_bytes ||
          expected_bytes > decoder.remaining())
        return {StatusCode::io_failure, kRestartRankFile};
      candidate.cell_records.resize(expected_bytes);
      if (expected_bytes != 0U &&
          !decoder.bytes(candidate.cell_records.data(), expected_bytes))
        return {StatusCode::io_failure, kRestartRankFile};
    }

    std::uint64_t integrity = 0U;
    if (!decoder.u64(integrity) || decoder.remaining() != 0U) {
      return {StatusCode::io_failure, kRestartRankFile};
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kRestartRankFile};
  } catch (...) {
    return {StatusCode::io_failure, kRestartRankFile};
  }
}

bool same_common(const Manifest& left, const Manifest& right) noexcept {
  if (left.format_version != right.format_version ||
      !same(left.global_cells, right.global_cells) || left.plan != right.plan ||
      left.schema != right.schema || left.geometry != right.geometry ||
      left.time != right.time || left.dt != right.dt ||
      left.pressure_reference != right.pressure_reference ||
      left.step != right.step ||
      left.controller_state != right.controller_state ||
      left.previous_pressure_reference != right.previous_pressure_reference ||
      left.closed_mass_target != right.closed_mass_target ||
      left.method_history_signature != right.method_history_signature ||
      left.cell_record_identity != right.cell_record_identity ||
      left.cell_record_bytes != right.cell_record_bytes ||
      left.final_mass_flux_revision != right.final_mass_flux_revision ||
      left.previous_mass_flux_revision != right.previous_mass_flux_revision ||
      left.fields.size() != right.fields.size() ||
      left.rate_fields.size() != right.rate_fields.size()) {
    return false;
  }
  const auto same_metadata = [](const std::vector<FieldMeta>& a,
                                const std::vector<FieldMeta>& b) noexcept {
    if (a.size() != b.size()) return false;
    for (std::size_t index = 0U; index < a.size(); ++index)
      if (a[index].role != b[index].role ||
          a[index].field != b[index].field ||
          a[index].components != b[index].components)
        return false;
    return true;
  };
  return same_metadata(left.fields, right.fields) &&
         same_metadata(left.rate_fields, right.rate_fields);
}

Status broadcast_bytes(MPI_Comm communicator, int rank,
                       std::vector<std::uint8_t>& bytes,
                       std::size_t maximum_bytes = SIZE_MAX) noexcept {
  std::uint64_t size = rank == 0 ? bytes.size() : 0U;
  if (MPI_Bcast(&size, 1, MPI_UINT64_T, 0, communicator) != MPI_SUCCESS ||
      size == 0U || size > static_cast<std::uint64_t>(INT_MAX))
    return {StatusCode::mpi_failure, kRestartCollective};
  Status status;
  status = collective_status(communicator, size <= maximum_bytes ? Status{} :
      Status{StatusCode::allocation_failure, kRestartReadBudget});
  if (!status) return status;
  try {
    if (rank != 0) bytes.resize(static_cast<std::size_t>(size));
  } catch (...) {
    status = {StatusCode::allocation_failure, kRestartCollective};
  }
  status = collective_status(communicator, status);
  if (!status) return status;
  return MPI_Bcast(bytes.data(), static_cast<int>(bytes.size()), MPI_BYTE, 0,
                   communicator) == MPI_SUCCESS
             ? Status{}
             : Status{StatusCode::mpi_failure, kRestartCollective};
}

Status broadcast_string(MPI_Comm communicator, int rank,
                        std::string& value) noexcept {
  std::vector<std::uint8_t> bytes;
  Status status;
  try {
    if (rank == 0) bytes.assign(value.begin(), value.end());
  } catch (...) {
    status = {StatusCode::allocation_failure, kRestartCollective};
  }
  status = collective_status(communicator, status);
  if (!status) return status;
  status = broadcast_bytes(communicator, rank, bytes, kCurrentMaximumBytes);
  if (status && rank != 0) {
    try {
      value.assign(bytes.begin(), bytes.end());
    } catch (...) {
      status = {StatusCode::allocation_failure, kRestartCollective};
    }
  }
  return collective_status(communicator, status);
}

Status read_current_name(const fs::path& directory, std::string& out,
                         IoFailureContext* failure) noexcept
    try {
  std::vector<std::uint8_t> bytes;
  if (!read_file(directory / "current", bytes, kCurrentMaximumBytes, 0U, failure) || bytes.empty())
    return {StatusCode::io_failure, kRestartDirectory};
  if (bytes.back() == '\n') bytes.pop_back();
  if (bytes.empty()) return {StatusCode::io_failure, kRestartDirectory};
  for (std::uint8_t value : bytes) {
    const char character = static_cast<char>(value);
    if (!(character == '-' || (character >= '0' && character <= '9') ||
          (character >= 'a' && character <= 'z'))) {
      return {StatusCode::io_failure, kRestartDirectory};
    }
  }
  try {
    out.assign(bytes.begin(), bytes.end());
  } catch (...) {
    return {StatusCode::allocation_failure, kRestartDirectory};
  }
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kRestartDirectory};
} catch (...) {
  return {StatusCode::io_failure, kRestartDirectory};
}

Status validate_expected(const RestartExpected& expected,
                         const Manifest& manifest) noexcept {
  const bool current_identity =
      expected.plan == manifest.plan && expected.schema == manifest.schema;
  const bool compatible_identity =
      manifest.format_version >= kExactHistoryFormatVersion &&
      expected.compatible_storage_plan != 0U &&
      expected.compatible_storage_schema != 0U &&
      expected.compatible_storage_plan == manifest.plan &&
      expected.compatible_storage_schema == manifest.schema;
  const bool compatible_method =
      manifest.format_version >= kExactHistoryFormatVersion &&
      expected.compatible_method_plan != 0U &&
      expected.compatible_method_plan == manifest.plan &&
      expected.schema == manifest.schema;
  if (!valid_global_patch(expected.global_cells, expected.target_patch) ||
      !same(expected.global_cells, manifest.global_cells) ||
      expected.plan == 0U || expected.schema == 0U || expected.geometry == 0U ||
      (!current_identity && !compatible_identity && !compatible_method) ||
      expected.geometry != manifest.geometry ||
      expected.cell_record_identity != manifest.cell_record_identity ||
      expected.cell_record_bytes != manifest.cell_record_bytes ||
      expected.fields.data == nullptr ||
      expected.fields.size != manifest.fields.size()) {
    return {StatusCode::invalid_plan, kRestartMismatch};
  }
  for (std::size_t index = 0U; index < expected.fields.size; ++index) {
    const RestartExpectedField field = expected.fields.data[index];
    const FieldMeta stored = manifest.fields[index];
    if (field.role != stored.role || field.field != stored.field ||
        field.components != stored.components || field.components == 0U) {
      return {StatusCode::invalid_plan, kRestartMismatch};
    }
  }
  if (manifest.format_version >= kExactHistoryFormatVersion) {
    if (expected.rate_fields.data == nullptr ||
        expected.rate_fields.size != manifest.rate_fields.size())
      return {StatusCode::invalid_plan, kRestartMismatch};
    for (std::size_t index = 0U; index < expected.rate_fields.size; ++index) {
      const RestartExpectedField field = expected.rate_fields.data[index];
      const FieldMeta stored = manifest.rate_fields[index];
      if (field.role != stored.role || field.field != stored.field ||
          field.components != stored.components || field.components != 1U)
        return {StatusCode::invalid_plan, kRestartMismatch};
    }
  }
  return {};
}

struct GenerationKey {
  std::uint64_t step{};
  std::uint64_t tick{};
};

bool generation_key(std::string_view name, GenerationKey& key) noexcept {
  constexpr std::string_view prefix = "generation-";
  if (name.substr(0U, prefix.size()) != prefix) return false;
  name.remove_prefix(prefix.size());
  const auto dash = name.find('-');
  if (dash == 0U || dash == std::string_view::npos || dash + 1U == name.size())
    return false;
  const auto step = std::from_chars(name.data(), name.data() + dash, key.step);
  const auto tick = std::from_chars(name.data() + dash + 1U,
                                   name.data() + name.size(), key.tick);
  return step.ec == std::errc{} && step.ptr == name.data() + dash &&
         tick.ec == std::errc{} && tick.ptr == name.data() + name.size();
}

Status prune_generations(const fs::path& directory, std::uint32_t keep_last,
                          const std::string& current, IoFailureContext* failure) {
  std::error_code error;
  struct Generation {
    fs::path path;
    GenerationKey key;
  };
  std::vector<Generation> generations;
  for (fs::directory_iterator iterator(directory, error), end;
       !error && iterator != end; iterator.increment(error)) {
    const std::string name = iterator->path().filename().string();
    GenerationKey key;
    if (fs::is_directory(iterator->symlink_status(error)) && !error &&
        generation_key(name, key)) {
      generations.push_back({iterator->path(), key});
    }
  }
  if (error) {
    detail::output_record_failure(failure, IoFailureOperation::read, error.value(), directory);
    return {StatusCode::io_failure, kRestartPublication};
  }
  // Keep current even after a rollback restart; the remaining generations
  // are ranked by accepted step, with a deterministic numeric tick tie-break.
  std::sort(generations.begin(), generations.end(),
            [](const Generation& a, const Generation& b) {
              return a.key.step < b.key.step ||
                     (a.key.step == b.key.step && a.key.tick < b.key.tick);
            });
  while (generations.size() > keep_last) {
    auto selected = generations.begin();
    while (selected != generations.end() &&
           selected->path.filename().string() == current)
      ++selected;
    if (selected == generations.end()) break;
    fs::remove_all(selected->path, error);
    if (error) {
      detail::output_record_failure(failure, IoFailureOperation::remove,
                                    error.value(), selected->path);
      return {StatusCode::io_failure, kRestartPublication};
    }
    generations.erase(selected);
  }
  return sync_directory(directory, failure)
             ? Status{} : Status{StatusCode::io_failure, kRestartPublication};
}

bool remove_stale_pending(const fs::path& directory, IoFailureContext* failure) {
  std::error_code error;
  for (fs::directory_iterator iterator(directory, error), end;
       !error && iterator != end; iterator.increment(error)) {
    const std::string name = iterator->path().filename().string();
    if (name.size() >= 8U &&
        name.compare(name.size() - 8U, 8U, "-pending") == 0) {
      std::string_view owned(name.data(), name.size() - 8U);
      const bool pointer = owned.substr(0U, 8U) == "current-";
      if (pointer) owned.remove_prefix(8U);
      GenerationKey key;
      if (!generation_key(owned, key)) continue;
      const auto type = iterator->symlink_status(error);
      if (error) return false;
      if (pointer ? !fs::is_regular_file(type) : !fs::is_directory(type))
        continue;
      fs::remove_all(iterator->path(), error);
      if (error) {
        detail::output_record_failure(failure, IoFailureOperation::remove,
                                      error.value(), iterator->path());
        return false;
      }
    }
  }
  if (error)
    detail::output_record_failure(failure, IoFailureOperation::read, error.value(), directory);
  return !error && sync_directory(directory, failure);
}

Status publish_generation(const fs::path& directory,
                          const fs::path& pending,
                          const std::string& generation,
                          RestartWriteReport& report) {
  std::error_code error;
  const fs::path final = directory / generation;
  fs::rename(pending, final, error);
  if (error)
    detail::output_record_failure(&report.failure, IoFailureOperation::rename,
                                  error.value(), final);
  if (error || !sync_directory(directory, &report.failure))
    return {StatusCode::io_failure, kRestartPublication};
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (injected(detail::RestartFailurePoint::after_generation_rename, 0))
    return {StatusCode::io_failure, kRestartPublication};
#endif
  const std::string pointer_text = generation + "\n";
  std::vector<std::uint8_t> pointer(pointer_text.begin(), pointer_text.end());
  const fs::path pointer_pending = directory / ("current-" + generation + "-pending");
  if (!write_file_sync(pointer_pending, pointer, &report.failure))
    return {StatusCode::io_failure, kRestartPublication};
  const auto current = directory / "current";
  fs::rename(pointer_pending, current, error);
  if (error) {
    detail::output_record_failure(&report.failure, IoFailureOperation::rename,
                                  error.value(), current);
    return {StatusCode::io_failure, kRestartPublication};
  }
  report.publication = RestartPublicationState::visible_not_durable;
  if (!sync_directory(directory, &report.failure))
    return {StatusCode::io_failure, kRestartPublication};
  report.publication = RestartPublicationState::durable;
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (injected(detail::RestartFailurePoint::after_current_switch, 0))
    return {StatusCode::io_failure, kRestartPublication};
#endif
  return {};
}

}  // namespace

#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
namespace detail {

void set_restart_failure_for_test(RestartFailurePoint point,
                                  int rank) noexcept {
  g_restart_failure_rank.store(rank, std::memory_order_relaxed);
  g_restart_failure_point.store(static_cast<int>(point),
                                std::memory_order_release);
}

void clear_restart_failure_for_test() noexcept {
  g_restart_failure_point.store(static_cast<int>(RestartFailurePoint::none),
                                std::memory_order_release);
  g_restart_failure_rank.store(-1, std::memory_order_relaxed);
}

}  // namespace detail
#endif

void RestartImage::clear() noexcept {
  storage_layout_migrated = false;
  global_cells = {};
  patch = {};
  plan = 0U;
  schema = 0U;
  geometry = 0U;
  time = 0.0;
  dt = 0.0;
  pressure_reference = 0.0;
  step = 0U;
  controller_state = 0U;
  source_manifest_sha256 = {};
  fields.clear();
  for (std::vector<double>& axis : final_mass_flux) axis.clear();
  backward_euler_recovery = true;
  previous_fields.clear();
  accepted_rate_fields.clear();
  previous_rate_fields.clear();
  for (std::vector<double>& axis : previous_mass_flux) axis.clear();
  previous_pressure_reference = 0.0;
  closed_mass_target = 0.0;
  final_mass_flux_revision = 0U;
  previous_mass_flux_revision = 0U;
  storage_layout_migrated = false;
  source_format_version = kLegacyFormatVersion;
  method_history_signature = 0U;
  cell_record_identity = 0U;
  cell_record_bytes = 0U;
  cell_records.clear();
  cell_record_lengths.clear();
}

Status RestartWriter::write(MPI_Comm communicator,
                            const std::filesystem::path& restart_directory,
                            const RestartSnapshot& snapshot,
                            RestartWriteOptions options) noexcept {
  struct Capture {
    RestartWriteReport value{};
    RestartWriteReport* destination;
    ~Capture() noexcept { if (destination != nullptr) *destination = value; }
  } capture{{}, options.report};
  auto& report = capture.value;
  const auto local_stage = [&](auto&& work) noexcept {
    return restart_local_stage(communicator, std::forward<decltype(work)>(work),
                               &report.failure);
  };
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::invalid_plan, kRestartInput};
  }
  Status status = !restart_directory.empty() && options.keep_last != 0U &&
                          valid_snapshot(snapshot)
                      ? Status{}
                      : Status{StatusCode::invalid_plan, kRestartInput};
  status = collective_status(communicator, status);
  if (!status) return status;
  std::uint64_t signature = 0U;
  status = local_stage( [&]() -> Status {
    signature = snapshot_signature(snapshot);
    return {};
  });
  if (!status) return status;
  status = consensus_u64(communicator, signature);
  if (!status) return status;

  std::size_t local_bytes = 0U;
  status = local_stage([&]() -> Status {
    return rank_block_size(snapshot, local_bytes) ? Status{}
        : Status{StatusCode::invalid_plan, kRestartInput};
  });
  if (!status) return status;
  std::uint64_t maximum_bytes = local_bytes;
  if (MPI_Allreduce(MPI_IN_PLACE, &maximum_bytes, 1, MPI_UINT64_T,
                     MPI_MAX, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRestartCollective};
  status = local_stage([&]() -> Status {
    std::size_t required = local_bytes;
    if (rank == 0) {
      Encoder common;
      encode_common(common, snapshot, snapshot_format(snapshot));
      std::size_t metadata = 0U;
      if (!checked_multiply(static_cast<std::size_t>(size),
                             64U + sizeof(RankRecord) + 40U, metadata) ||
          metadata > SIZE_MAX - common.data().size() - 24U ||
          maximum_bytes > SIZE_MAX - metadata - common.data().size() - 24U)
        return {StatusCode::invalid_plan, kRestartInput};
      required = static_cast<std::size_t>(maximum_bytes) + metadata + common.data().size() + 24U;
    }
    return options.maximum_bulk_staging_bytes != 0U &&
                   required > options.maximum_bulk_staging_bytes
        ? Status{StatusCode::allocation_failure, kRestartRankFile} : Status{};
  });
  if (!status) return status;

  std::string generation;
  status = local_stage( [&]() -> Status {
    if (rank == 0) {
      const auto tick = static_cast<std::uint64_t>(
          std::chrono::steady_clock::now().time_since_epoch().count());
      generation = "generation-" + std::to_string(snapshot.step) + "-" +
                   std::to_string(tick);
    }
    return {};
  });
  if (!status) return status;
  status = broadcast_string(communicator, rank, generation);
  if (!status) return status;
  fs::path pending;
  status = local_stage( [&]() -> Status {
    pending = restart_directory / (generation + "-pending");
    return {};
  });
  if (!status) return status;
  status = local_stage( [&]() -> Status {
    if (rank != 0) return {};
    std::error_code error;
    fs::create_directories(restart_directory, error);
    if (!error && !remove_stale_pending(restart_directory, &report.failure))
      error = std::make_error_code(std::errc::io_error);
    if (!error) fs::create_directory(pending, error);
    if (error) detail::output_record_failure(&report.failure,
        IoFailureOperation::create_directory, error.value(), pending);
    return !error && sync_directory(restart_directory, &report.failure)
                 ? Status{}
                 : Status{StatusCode::io_failure, kRestartDirectory};
  });
  if (!status) return status;
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (injected(detail::RestartFailurePoint::after_directory, rank))
    status = {StatusCode::io_failure, kRestartDirectory};
  status = collective_status(communicator, status);
  if (!status) return status;
#endif

  std::vector<std::uint8_t> rank_bytes;
  status = local_stage( [&]() -> Status {
    Status local = encode_rank_block(snapshot, size, rank, rank_bytes,
                                     options.maximum_bulk_staging_bytes);
    report.rank_payload_bytes = rank_bytes.size();
    report.peak_bulk_staging_bytes = rank_bytes.capacity();
    if (local &&
        !write_file_sync(pending / rank_name(static_cast<std::uint32_t>(rank)),
                       rank_bytes, &report.failure)) {
      local = {StatusCode::io_failure, kRestartRankFile};
    }
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
    if (local && injected(detail::RestartFailurePoint::after_rank_file, rank))
      local = {StatusCode::io_failure, kRestartRankFile};
#endif
    return local;
  });
  if (!status) return status;

  const std::array<std::uint64_t, 8U> local_record{{
      static_cast<std::uint64_t>(snapshot.patch.begin.x),
      static_cast<std::uint64_t>(snapshot.patch.begin.y),
      static_cast<std::uint64_t>(snapshot.patch.begin.z),
      static_cast<std::uint64_t>(snapshot.patch.cells.x),
      static_cast<std::uint64_t>(snapshot.patch.cells.y),
      static_cast<std::uint64_t>(snapshot.patch.cells.z),
      static_cast<std::uint64_t>(rank_bytes.size()),
      hash_bytes(rank_bytes.data(), rank_bytes.size())}};
  // The durable rank file and its fixed-size record now own the information.
  // Do not keep our full payload alive while root verifies every rank file.
  std::vector<std::uint8_t>().swap(rank_bytes);
  std::vector<std::uint64_t> gathered;
  status = local_stage( [&]() -> Status {
    if (rank == 0) gathered.resize(static_cast<std::size_t>(size) * 8U);
    return {};
  });
  if (!status) return status;
  if (MPI_Gather(local_record.data(), static_cast<int>(local_record.size()),
                 MPI_UINT64_T, rank == 0 ? gathered.data() : nullptr,
                 static_cast<int>(local_record.size()), MPI_UINT64_T, 0,
                 communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kRestartCollective};
  }
  status = local_stage( [&]() -> Status {
    if (rank != 0) return {};
    std::vector<RankRecord> records(static_cast<std::size_t>(size));
    for (int source = 0; source < size; ++source) {
      const std::size_t base = static_cast<std::size_t>(source) * 8U;
      records[static_cast<std::size_t>(source)] = {
          {static_cast<std::int32_t>(gathered[base]),
           static_cast<std::int32_t>(gathered[base + 1U]),
           static_cast<std::int32_t>(gathered[base + 2U])},
          {static_cast<std::int32_t>(gathered[base + 3U]),
           static_cast<std::int32_t>(gathered[base + 4U]),
           static_cast<std::int32_t>(gathered[base + 5U])},
          gathered[base + 6U], gathered[base + 7U]};
    }
    std::vector<std::uint8_t> manifest_bytes;
    Status local = encode_manifest(snapshot, size, records, manifest_bytes);
    if (local && !write_file_sync(pending / "manifest.bin", manifest_bytes, &report.failure))
      local = {StatusCode::io_failure, kRestartManifest};
    std::size_t maximum_rank_bytes = 0U;
    for (const auto& record : records) {
      if (record.bytes > SIZE_MAX) return {StatusCode::invalid_plan, kRestartInput};
      maximum_rank_bytes = std::max(maximum_rank_bytes, static_cast<std::size_t>(record.bytes));
    }
    const std::size_t metadata_bytes = gathered.capacity() * sizeof(std::uint64_t) +
        records.capacity() * sizeof(RankRecord) + manifest_bytes.capacity();
    if (maximum_rank_bytes > SIZE_MAX - metadata_bytes)
      return {StatusCode::invalid_plan, kRestartInput};
    const std::size_t peak = maximum_rank_bytes + metadata_bytes;
    if (options.maximum_bulk_staging_bytes != 0U &&
        peak > options.maximum_bulk_staging_bytes)
      return {StatusCode::allocation_failure, kRestartIntegrity};
    std::vector<std::uint8_t> verify;
    verify.reserve(maximum_rank_bytes);
    report.peak_bulk_staging_bytes = std::max(report.peak_bulk_staging_bytes, peak);
    for (int source = 0; source < size && local; ++source) {
      const RankRecord& record = records[static_cast<std::size_t>(source)];
      if (!read_file(pending / rank_name(static_cast<std::uint32_t>(source)),
                     verify, maximum_rank_bytes, record.bytes, &report.failure) ||
          verify.size() != record.bytes ||
          hash_bytes(verify.data(), verify.size()) != record.hash ||
          !verified_integrity(verify)) {
        local = {StatusCode::io_failure, kRestartIntegrity};
      }
    }
    if (local && !sync_directory(pending, &report.failure))
      local = {StatusCode::io_failure, kRestartDirectory};
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
    if (local &&
        injected(detail::RestartFailurePoint::after_manifest, rank))
      local = {StatusCode::io_failure, kRestartManifest};
#endif
    return local;
  });
  if (!status) return status;
  status = local_stage([&]() -> Status {
    return rank == 0
               ? publish_generation(restart_directory, pending, generation, report)
               : Status{};
  });
  auto publication = static_cast<std::uint8_t>(report.publication);
  if (MPI_Bcast(&publication, 1, MPI_UINT8_T, 0, communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kRestartCollective};
  report.publication = static_cast<RestartPublicationState>(publication);
  if (!status) return status;
  report.cleanup_status = restart_local_stage(communicator, [&]() -> Status {
    return rank == 0
               ? prune_generations(restart_directory, options.keep_last,
                                   generation, &report.cleanup_failure)
               : Status{};
  }, &report.cleanup_failure);
  // current is durable. A cleanup warning must not be reported as a failed
  // publication or trigger an attempt to overwrite/re-publish this checkpoint.
  return report.cleanup_status.code == StatusCode::mpi_failure
             ? report.cleanup_status : Status{};
}

Status RestartReader::load(MPI_Comm communicator,
                           const std::filesystem::path& restart_directory,
                           const RestartExpected& expected,
                           RestartImage& out,
                           RestartReadReport* report,
                           RestartReadLimits limits) noexcept try {
  detail::IoFailureCapture failure_capture(report ? &report->failure : nullptr);
  if (report) *report = {};
  const auto local_stage = [&](auto&& work) noexcept {
    return restart_local_stage(communicator, std::forward<decltype(work)>(work),
                                &failure_capture.context);
  };
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::invalid_plan, kRestartInput};
  }
  RestartReadReport budget;
  Status status = local_stage([&]() -> Status {
    if (limits.maximum_bulk_bytes == 0U || limits.maximum_manifest_bytes == 0U ||
        limits.maximum_source_ranks == 0U)
      return {StatusCode::invalid_plan, kRestartReadBudget};
    BulkBytes initial;
    if (!retained_image_bytes(out, budget.retained_image_bytes) ||
        !initial.add(budget.retained_image_bytes) || !initial.add(kCurrentMaximumBytes))
      return {StatusCode::allocation_failure, kRestartReadBudget};
    budget.peak_bulk_bytes = initial.value;
    return initial.value <= limits.maximum_bulk_bytes ? Status{} :
        Status{StatusCode::allocation_failure, kRestartReadBudget};
  });
  const auto publish_budget = [&]() noexcept {
    if (!report) return;
    report->retained_image_bytes = budget.retained_image_bytes;
    report->new_image_bytes = budget.new_image_bytes;
    report->peak_bulk_bytes = budget.peak_bulk_bytes;
  };
  publish_budget();
  if (!status) return status;
  const auto manifest_limit = std::min({limits.maximum_manifest_bytes,
      limits.maximum_bulk_bytes - budget.retained_image_bytes - kCurrentMaximumBytes,
      static_cast<std::size_t>(INT_MAX)});
  std::string generation;
  std::vector<std::uint8_t> manifest_bytes;
  status = local_stage([&]() -> Status {
    if (restart_directory.empty()) return {StatusCode::invalid_plan, kRestartInput};
    Status local;
    if (rank == 0) {
      local = read_current_name(restart_directory, generation, &failure_capture.context);
      if (local && !read_file(restart_directory / generation / "manifest.bin",
                              manifest_bytes, manifest_limit, 0U, &failure_capture.context))
        local = {StatusCode::io_failure, kRestartManifest};
    }
    return local;
  });
  if (!status) return status;
  status = broadcast_string(communicator, rank, generation);
  if (!status) return status;
  status = broadcast_bytes(communicator, rank, manifest_bytes, manifest_limit);
  status = collective_status(communicator, status);
  if (!status) return status;
  Manifest manifest;
  RuntimeSha256Digest source_digest{};
  status = local_stage([&]() -> Status {
    BulkBytes parsing{budget.retained_image_bytes};
    // The encoded file itself bounds the rank directory count. Check this
    // conservative metadata allocation bound before invoking the parser.
    if (!parsing.add(manifest_bytes.capacity()) || !parsing.add(kCommonMetadataBound) ||
        !parsing.add(manifest_bytes.size() / 40U, sizeof(RankRecord)))
      return {StatusCode::allocation_failure, kRestartReadBudget};
    budget.peak_bulk_bytes = std::max(budget.peak_bulk_bytes, parsing.value);
    if (parsing.value > limits.maximum_bulk_bytes)
      return {StatusCode::allocation_failure, kRestartReadBudget};
    Status local = parse_manifest(manifest_bytes, manifest, limits.maximum_source_ranks);
    if (local) local = validate_expected(expected, manifest);
    if (local && !detail::runtime_sha256_bytes(
          {manifest_bytes.data(), manifest_bytes.size()}, source_digest))
      local = {StatusCode::io_failure, kRestartIntegrity};
    if (!local) return local;
    // The digest and parsed manifest suffice for all following integrity checks.
    // Do not keep another raw manifest allocation alive beside the new image.
    std::vector<std::uint8_t>().swap(manifest_bytes);
    return plan_read_bulk(manifest, expected, budget.retained_image_bytes, limits, budget);
  });
  publish_budget();
  if (!status) return status;

  fs::path generation_directory;
  status = local_stage([&]() -> Status {
    generation_directory = restart_directory / generation;
    for (std::uint32_t source = static_cast<std::uint32_t>(rank);
         source < manifest.rank_count && status;
         source += static_cast<std::uint32_t>(size)) {
      std::vector<std::uint8_t> bytes;
      RankBlock block;
      const RankRecord record = manifest.ranks[source];
      if (!read_file(generation_directory / rank_name(source), bytes,
                      static_cast<std::size_t>(record.bytes), record.bytes, &failure_capture.context) ||
          bytes.size() != record.bytes ||
          hash_bytes(bytes.data(), bytes.size()) != record.hash) {
        status = {StatusCode::io_failure, kRestartIntegrity};
        break;
      }
      status = parse_rank_block(bytes, block, manifest, source);
      if (report) {
        ++report->integrity_blocks;
        report->rank_file_bytes_read += bytes.size();
      }
      if (status &&
          (block.rank_count != manifest.rank_count || block.rank != source ||
           !same(block.patch.begin, record.begin) ||
           !same(block.patch.cells, record.cells) ||
           !same_common(block.common, manifest))) {
        status = {StatusCode::io_failure, kRestartMismatch};
      }
    }
    return status;
  });
  if (!status) return status;

  RestartImage candidate;
  std::vector<std::size_t> target_record_offsets;
  if (manifest.format_version == kVariableCellRecordFormatVersion) {
    // First collect only the lengths owned by the target patch. The initial
    // budget already includes these arrays and the largest decoded block.
    // Allocate the exact target payload only after this bounded census.
    status = local_stage([&]() -> Status {
      std::size_t cells = 0U;
      if (!cell_count(expected.target_patch.cells, cells))
        return {StatusCode::invalid_plan, kRestartInput};
      candidate.cell_record_lengths.resize(cells);
      target_record_offsets.resize(cells + 1);
      const auto index = [](Int3 p, Int3 n) {
        return (std::size_t(p.z) * n.y + p.y) * n.x + p.x;
      };
      for (std::uint32_t source = 0; source < manifest.rank_count; ++source) {
        const auto &record = manifest.ranks[source];
        const auto &target = expected.target_patch;
        const Int3 begin{std::max(record.begin.x, target.begin.x),
                         std::max(record.begin.y, target.begin.y),
                         std::max(record.begin.z, target.begin.z)};
        const Int3 end{std::min(record.begin.x + record.cells.x,
                                target.begin.x + target.cells.x),
                       std::min(record.begin.y + record.cells.y,
                                target.begin.y + target.cells.y),
                       std::min(record.begin.z + record.cells.z,
                                target.begin.z + target.cells.z)};
        if (begin.x >= end.x || begin.y >= end.y || begin.z >= end.z)
          continue;
        std::vector<std::uint8_t> bytes;
        RankBlock block;
        if (!read_file(generation_directory / rank_name(source), bytes,
                       std::size_t(record.bytes), record.bytes,
                       &failure_capture.context) ||
            hash_bytes(bytes.data(), bytes.size()) != record.hash)
          return {StatusCode::io_failure, kRestartRankFile};
        const auto parsed = parse_rank_block(bytes, block, manifest, source);
        if (report)
          report->rank_file_bytes_read += bytes.size();
        if (!parsed)
          return parsed;
        for (int z = begin.z; z < end.z; ++z)
          for (int y = begin.y; y < end.y; ++y)
            for (int x = begin.x; x < end.x; ++x) {
              const auto old = index(
                  {x - record.begin.x, y - record.begin.y, z - record.begin.z},
                  record.cells);
              const auto current = index(
                  {x - target.begin.x, y - target.begin.y, z - target.begin.z},
                  target.cells);
              candidate.cell_record_lengths[current] =
                  block.cell_record_lengths[old];
            }
      }
      for (std::size_t i = 0; i < cells; ++i) {
        if (candidate.cell_record_lengths[i] >
            SIZE_MAX - target_record_offsets[i])
          return {StatusCode::allocation_failure, kRestartReadBudget};
        target_record_offsets[i + 1] =
            target_record_offsets[i] + candidate.cell_record_lengths[i];
      }
      return plan_read_bulk(manifest, expected, budget.retained_image_bytes,
                            limits, budget, target_record_offsets.back());
    });
    publish_budget();
    if (!status)
      return status;
  }
  status = local_stage([&]() -> Status {
    candidate.global_cells = manifest.global_cells;
    candidate.patch = expected.target_patch;
    candidate.plan = manifest.plan;
    candidate.schema = manifest.schema;
    candidate.storage_layout_migrated =
        (manifest.plan != expected.plan || manifest.schema != expected.schema) &&
        manifest.plan == expected.compatible_storage_plan &&
        manifest.schema == expected.compatible_storage_schema;
    candidate.geometry = manifest.geometry;
    candidate.time = manifest.time;
    candidate.dt = manifest.dt;
    candidate.pressure_reference = manifest.pressure_reference;
    candidate.step = manifest.step;
    candidate.controller_state = manifest.controller_state;
    candidate.source_format_version = manifest.format_version;
    candidate.method_history_signature = manifest.method_history_signature;
    candidate.cell_record_identity = manifest.cell_record_identity;
    candidate.cell_record_bytes = manifest.cell_record_bytes;
    candidate.source_manifest_sha256 = source_digest;
    const bool exact_history =
        manifest.format_version >= kExactHistoryFormatVersion;
    candidate.backward_euler_recovery = !exact_history;
    candidate.previous_pressure_reference =
        manifest.previous_pressure_reference;
    candidate.closed_mass_target = manifest.closed_mass_target;
    candidate.final_mass_flux_revision = manifest.final_mass_flux_revision;
    candidate.previous_mass_flux_revision =
        manifest.previous_mass_flux_revision;
    std::size_t target_cells = 0U;
    if (!cell_count(expected.target_patch.cells, target_cells))
      status = {StatusCode::invalid_plan, kRestartInput};
    if (status && manifest.format_version >= kCellRecordFormatVersion) {
      std::size_t bytes = 0U;
      if (manifest.format_version == kVariableCellRecordFormatVersion)
        bytes = target_record_offsets.back();
      else if (!checked_multiply(target_cells, manifest.cell_record_bytes,
                                 bytes))
        return {StatusCode::allocation_failure, kRestartReadBudget};
      candidate.cell_records.resize(bytes);
    }
    const auto allocate_fields = [&](const std::vector<FieldMeta>& metadata,
                                     std::vector<RestartImageField>& fields) {
      fields.resize(metadata.size());
      for (std::size_t index = 0U; index < fields.size(); ++index) {
        const FieldMeta meta = metadata[index];
        RestartImageField& field = fields[index];
        field.role = meta.role;
        field.field = meta.field;
        field.components = meta.components;
        std::size_t values = 0U;
        if (!checked_multiply(target_cells, meta.components, values))
          return false;
        field.values.resize(values);
      }
      return true;
    };
    if (status && !allocate_fields(manifest.fields, candidate.fields))
      status = {StatusCode::invalid_plan, kRestartInput};
    if (status && exact_history &&
        (!allocate_fields(manifest.fields, candidate.previous_fields) ||
         !allocate_fields(manifest.rate_fields,
                          candidate.accepted_rate_fields) ||
         !allocate_fields(manifest.rate_fields,
                          candidate.previous_rate_fields)))
      status = {StatusCode::invalid_plan, kRestartInput};
    std::vector<std::uint8_t> cell_coverage;
    std::array<std::vector<std::uint8_t>, 3U> face_coverage;
    std::array<Int3, 3U> target_face_extents{
        {{expected.target_patch.cells.x + 1, expected.target_patch.cells.y,
          expected.target_patch.cells.z},
         {expected.target_patch.cells.x, expected.target_patch.cells.y + 1,
          expected.target_patch.cells.z},
         {expected.target_patch.cells.x, expected.target_patch.cells.y,
          expected.target_patch.cells.z + 1}}};
    if (status) {
      cell_coverage.assign(target_cells, 0U);
      for (std::size_t axis = 0U; axis < 3U; ++axis) {
        std::size_t count = 0U;
        if (!cell_count(target_face_extents[axis], count)) {
          status = {StatusCode::invalid_plan, kRestartInput};
          break;
        }
        candidate.final_mass_flux[axis].assign(count, 0.0);
        if (exact_history)
          candidate.previous_mass_flux[axis].assign(count, 0.0);
        face_coverage[axis].assign(count, 0U);
      }
    }

    const auto dense_index = [](Int3 local, Int3 cells) noexcept {
      return (static_cast<std::size_t>(local.z) *
                  static_cast<std::size_t>(cells.y) +
              static_cast<std::size_t>(local.y)) *
                 static_cast<std::size_t>(cells.x) +
             static_cast<std::size_t>(local.x);
    };
    for (std::uint32_t source = 0U; source < manifest.rank_count && status;
         ++source) {
      // A cell-disjoint block can still own the target's upper interface face.
      // Filter using manifest-owned face extents as well as cell extents. The
      // preceding distributed integrity scan still checks EVERY source file.
      const auto& record = manifest.ranks[source];
      const auto intersects = [](Int3 a, Int3 na, Int3 b, Int3 nb) {
        return a.x < b.x + nb.x && b.x < a.x + na.x &&
               a.y < b.y + nb.y && b.y < a.y + na.y &&
               a.z < b.z + nb.z && b.z < a.z + na.z;
      };
      bool needed = intersects(record.begin, record.cells,
                                expected.target_patch.begin, expected.target_patch.cells);
      for (std::size_t axis = 0U; axis < 3U && !needed; ++axis) {
        Int3 owned = record.cells;
        if (axis == 0U && record.begin.x + owned.x == manifest.global_cells.x) ++owned.x;
        if (axis == 1U && record.begin.y + owned.y == manifest.global_cells.y) ++owned.y;
        if (axis == 2U && record.begin.z + owned.z == manifest.global_cells.z) ++owned.z;
        needed = intersects(record.begin, owned, expected.target_patch.begin,
                            target_face_extents[axis]);
      }
      if (!needed) continue;
      std::vector<std::uint8_t> bytes;
      RankBlock block;
      if (!read_file(generation_directory / rank_name(source), bytes,
                      static_cast<std::size_t>(record.bytes), record.bytes, &failure_capture.context) ||
          hash_bytes(bytes.data(), bytes.size()) != record.hash) {
        status = {StatusCode::io_failure, kRestartRankFile};
        break;
      }
      status = parse_rank_block(bytes, block, manifest, source);
      if (report) {
        ++report->restoration_blocks;
        report->rank_file_bytes_read += bytes.size();
      }
      if (!status) break;
      const Int3 begin{
          std::max(block.patch.begin.x, expected.target_patch.begin.x),
          std::max(block.patch.begin.y, expected.target_patch.begin.y),
          std::max(block.patch.begin.z, expected.target_patch.begin.z)};
      const Int3 end{std::min(block.patch.begin.x + block.patch.cells.x,
                              expected.target_patch.begin.x +
                                  expected.target_patch.cells.x),
                     std::min(block.patch.begin.y + block.patch.cells.y,
                              expected.target_patch.begin.y +
                                  expected.target_patch.cells.y),
                     std::min(block.patch.begin.z + block.patch.cells.z,
                              expected.target_patch.begin.z +
                                  expected.target_patch.cells.z)};
      for (std::int32_t z = begin.z; z < end.z; ++z)
        for (std::int32_t y = begin.y; y < end.y; ++y)
          for (std::int32_t x = begin.x; x < end.x; ++x) {
            const Int3 old_local{x - block.patch.begin.x,
                                 y - block.patch.begin.y,
                                 z - block.patch.begin.z};
            const Int3 new_local{x - expected.target_patch.begin.x,
                                 y - expected.target_patch.begin.y,
                                 z - expected.target_patch.begin.z};
            const std::size_t old_cell =
                dense_index(old_local, block.patch.cells);
            const std::size_t new_cell =
                dense_index(new_local, expected.target_patch.cells);
            if (cell_coverage[new_cell] != 0U) {
              status = {StatusCode::io_failure, kRestartCoverage};
              break;
            }
            cell_coverage[new_cell] = 1U;
            if (manifest.format_version >= kCellRecordFormatVersion) {
              std::size_t width = manifest.cell_record_bytes;
              std::size_t source_offset = old_cell * width,
                          target_offset = new_cell * width;
              if (manifest.format_version == kVariableCellRecordFormatVersion) {
                width = block.cell_record_lengths[old_cell];
                if (width != candidate.cell_record_lengths[new_cell]) {
                  status = {StatusCode::io_failure, kRestartMismatch};
                  break;
                }
                source_offset = block.cell_record_offsets[old_cell];
                target_offset = target_record_offsets[new_cell];
              }
              if (width != 0U)
                std::copy_n(block.cell_records.data() + source_offset, width,
                            candidate.cell_records.data() + target_offset);
            }
            const auto copy_fields =
                [&](std::vector<RestartImageField>& target,
                    const std::vector<RestartImageField>& source_fields) {
                  for (std::size_t field_index = 0U;
                       field_index < target.size(); ++field_index) {
                    const std::size_t components =
                        target[field_index].components;
                    for (std::size_t component = 0U; component < components;
                         ++component)
                      target[field_index]
                          .values[new_cell * components + component] =
                          source_fields[field_index]
                              .values[old_cell * components + component];
                  }
                };
            copy_fields(candidate.fields, block.fields);
            if (exact_history) {
              copy_fields(candidate.previous_fields, block.previous_fields);
              copy_fields(candidate.accepted_rate_fields,
                          block.accepted_rate_fields);
              copy_fields(candidate.previous_rate_fields,
                          block.previous_rate_fields);
            }
          }
      if (!status) break;

      for (std::size_t axis = 0U; axis < 3U && status; ++axis) {
        Int3 old_owned = block.patch.cells;
        const std::int32_t old_end =
            axis == 0U
                ? block.patch.begin.x + block.patch.cells.x
                : (axis == 1U ? block.patch.begin.y + block.patch.cells.y
                              : block.patch.begin.z + block.patch.cells.z);
        const std::int32_t global_end =
            axis == 0U ? manifest.global_cells.x
                       : (axis == 1U ? manifest.global_cells.y
                                     : manifest.global_cells.z);
        if (old_end == global_end) {
          if (axis == 0U)
            ++old_owned.x;
          else if (axis == 1U)
            ++old_owned.y;
          else
            ++old_owned.z;
        }
        for (std::int32_t z = 0; z < old_owned.z; ++z)
          for (std::int32_t y = 0; y < old_owned.y; ++y)
            for (std::int32_t x = 0; x < old_owned.x; ++x) {
              const Int3 global{block.patch.begin.x + x,
                                block.patch.begin.y + y,
                                block.patch.begin.z + z};
              const Int3 target_local{global.x - expected.target_patch.begin.x,
                                      global.y - expected.target_patch.begin.y,
                                      global.z - expected.target_patch.begin.z};
              const Int3 extent = target_face_extents[axis];
              if (target_local.x < 0 || target_local.y < 0 ||
                  target_local.z < 0 || target_local.x >= extent.x ||
                  target_local.y >= extent.y || target_local.z >= extent.z)
                continue;
              const std::size_t old_index = dense_index({x, y, z}, old_owned);
              const std::size_t target_index =
                  dense_index(target_local, extent);
              if (face_coverage[axis][target_index] != 0U) {
                status = {StatusCode::io_failure, kRestartCoverage};
                break;
              }
              face_coverage[axis][target_index] = 1U;
              candidate.final_mass_flux[axis][target_index] =
                  block.flux[axis][old_index];
              if (exact_history)
                candidate.previous_mass_flux[axis][target_index] =
                    block.previous_flux[axis][old_index];
            }
      }
    }
    if (status && std::any_of(cell_coverage.begin(), cell_coverage.end(),
                              [](std::uint8_t value) { return value != 1U; }))
      status = {StatusCode::io_failure, kRestartCoverage};
    for (std::size_t axis = 0U; axis < 3U && status; ++axis) {
      if (std::any_of(face_coverage[axis].begin(), face_coverage[axis].end(),
                      [](std::uint8_t value) { return value != 1U; }))
        status = {StatusCode::io_failure, kRestartCoverage};
    }
    return status;
  });
  if (!status) return status;
  out = std::move(candidate);
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kRestartInput};
} catch (...) {
  return {StatusCode::io_failure, kRestartInput};
}

}  // namespace hundun::v04
