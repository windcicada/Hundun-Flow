// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_io.hpp"

#include "app_identity_detail.hpp"
#include "field_view_interval_detail.hpp"
#include "io_restart_detail.hpp"

#include <mpi.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
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
constexpr std::uint32_t kLegacyFormatVersion = 1U;
constexpr std::uint32_t kExactHistoryFormatVersion = 2U;
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
  std::vector<FieldMeta> fields;
  std::vector<FieldMeta> rate_fields;
  std::vector<RankRecord> ranks;
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
};

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

bool write_file_sync(const fs::path& path,
                     const std::vector<std::uint8_t>& bytes) noexcept {
  const int descriptor =
      ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0644);
  if (descriptor < 0) return false;
  std::size_t cursor = 0U;
  bool okay = true;
  while (cursor < bytes.size()) {
    const ssize_t written =
        ::write(descriptor, bytes.data() + cursor, bytes.size() - cursor);
    if (written <= 0) {
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(written);
  }
  if (okay && ::fsync(descriptor) != 0) okay = false;
  if (::close(descriptor) != 0) okay = false;
  return okay;
}

bool read_file(const fs::path& path,
               std::vector<std::uint8_t>& bytes) noexcept {
  const int descriptor = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
  if (descriptor < 0) return false;
  struct stat info {};
  bool okay = ::fstat(descriptor, &info) == 0 && info.st_size >= 0;
  if (okay) {
    try {
      bytes.resize(static_cast<std::size_t>(info.st_size));
    } catch (...) {
      okay = false;
    }
  }
  std::size_t cursor = 0U;
  while (okay && cursor < bytes.size()) {
    const ssize_t count =
        ::read(descriptor, bytes.data() + cursor, bytes.size() - cursor);
    if (count <= 0) {
      okay = false;
      break;
    }
    cursor += static_cast<std::size_t>(count);
  }
  if (::close(descriptor) != 0) okay = false;
  return okay;
}

bool sync_directory(const fs::path& path) noexcept {
  const int descriptor = ::open(path.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
  if (descriptor < 0) return false;
  const bool okay = ::fsync(descriptor) == 0;
  return ::close(descriptor) == 0 && okay;
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
         snapshot.closed_mass_target != 0.0;
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

std::uint64_t snapshot_signature(const RestartSnapshot& snapshot) noexcept {
  Encoder encoder;
  const std::uint32_t version = has_exact_history(snapshot)
                                    ? kExactHistoryFormatVersion
                                    : kLegacyFormatVersion;
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
  if (version == kExactHistoryFormatVersion) {
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
  if (version == kExactHistoryFormatVersion) {
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
}

bool decode_common(Decoder& decoder, std::uint32_t version,
                   Manifest& manifest) noexcept {
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
  try {
    manifest.fields.resize(field_count);
  } catch (...) {
    return false;
  }
  for (FieldMeta& field : manifest.fields) {
    std::uint8_t role = 0U;
    if (!decoder.u8(role) || !decoder.u16(field.field) ||
        !decoder.u8(field.components) ||
        role > static_cast<std::uint8_t>(
                   RestartFieldRole::transported_scalar) ||
        field.components == 0U) {
      return false;
    }
    field.role = static_cast<RestartFieldRole>(role);
  }
  if (version == kExactHistoryFormatVersion) {
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
    try {
      manifest.rate_fields.resize(rate_count);
    } catch (...) {
      return false;
    }
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
  return true;
}

Status encode_rank_block(const RestartSnapshot& snapshot, int size, int rank,
                         std::vector<std::uint8_t>& out) {
  Encoder encoder;
  const std::uint32_t version = has_exact_history(snapshot)
                                    ? kExactHistoryFormatVersion
                                    : kLegacyFormatVersion;
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
  if (status && version == kExactHistoryFormatVersion)
    status = encode_fields(snapshot.previous_fields);
  if (status && version == kExactHistoryFormatVersion)
    status = encode_fields(snapshot.accepted_rate_fields);
  if (status && version == kExactHistoryFormatVersion)
    status = encode_fields(snapshot.previous_rate_fields);
  if (status) status = encode_flux(snapshot.final_mass_flux);
  if (status && version == kExactHistoryFormatVersion)
    status = encode_flux(snapshot.previous_mass_flux);
  if (!status) return status;
  encoder.append_integrity();
  out = encoder.data();
  return {};
}

Status encode_manifest(const RestartSnapshot& snapshot, int size,
                       const std::vector<RankRecord>& records,
                       std::vector<std::uint8_t>& out) {
  if (records.size() != static_cast<std::size_t>(size))
    return {StatusCode::invalid_plan, kRestartManifest};
  Encoder encoder;
  const std::uint32_t version = has_exact_history(snapshot)
                                    ? kExactHistoryFormatVersion
                                    : kLegacyFormatVersion;
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
  out = encoder.data();
  return {};
}

Status parse_manifest(const std::vector<std::uint8_t>& bytes,
                      Manifest& out) noexcept {
  if (!verified_integrity(bytes))
    return {StatusCode::io_failure, kRestartIntegrity};
  try {
    Decoder decoder(bytes);
    std::array<char, 8U> magic{};
    std::uint32_t version = 0U;
    Manifest candidate;
    if (!decoder.bytes(magic.data(), magic.size()) ||
        magic != kManifestMagic || !decoder.u32(version) ||
        (version != kLegacyFormatVersion &&
         version != kExactHistoryFormatVersion) ||
        !decoder.u32(candidate.rank_count) ||
        candidate.rank_count == 0U ||
        !decode_common(decoder, version, candidate) ||
        !valid_global_patch(candidate.global_cells,
                            MeshPatch{{0, 0, 0}, candidate.global_cells, {}, {}})) {
      return {StatusCode::io_failure, kRestartManifest};
    }
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

Status parse_rank_block(const std::vector<std::uint8_t>& bytes,
                        RankBlock& out) noexcept {
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
         version != kExactHistoryFormatVersion) ||
        !decoder.u32(candidate.rank_count) || candidate.rank_count == 0U ||
        !decoder.u32(candidate.rank) || candidate.rank >= candidate.rank_count ||
        !decode_common(decoder, version, candidate.common) ||
        !decoder.int3(candidate.patch.begin) ||
        !decoder.int3(candidate.patch.cells) ||
        !valid_global_patch(candidate.common.global_cells, candidate.patch)) {
      return {StatusCode::io_failure, kRestartRankFile};
    }
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
            !decoder.u64(stored_values) || stored_values != expected_values)
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
            !decoder.u64(stored_values) || stored_values != expected_values)
          return false;
        flux[axis].resize(expected_values);
        for (double& value : flux[axis])
          if (!decoder.real(value) || !std::isfinite(value)) return false;
      }
      return true;
    };
    if (!decode_fields(candidate.common.fields, candidate.fields))
      return {StatusCode::io_failure, kRestartRankFile};
    if (version == kExactHistoryFormatVersion &&
        (!decode_fields(candidate.common.fields, candidate.previous_fields) ||
         !decode_fields(candidate.common.rate_fields,
                        candidate.accepted_rate_fields) ||
         !decode_fields(candidate.common.rate_fields,
                        candidate.previous_rate_fields)))
      return {StatusCode::io_failure, kRestartRankFile};
    if (!decode_flux(candidate.flux) ||
        (version == kExactHistoryFormatVersion &&
         !decode_flux(candidate.previous_flux)))
      return {StatusCode::io_failure, kRestartRankFile};
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
      !same(left.global_cells, right.global_cells) ||
      left.plan != right.plan || left.schema != right.schema ||
      left.geometry != right.geometry || left.time != right.time ||
      left.dt != right.dt ||
      left.pressure_reference != right.pressure_reference ||
      left.step != right.step ||
      left.controller_state != right.controller_state ||
      left.previous_pressure_reference != right.previous_pressure_reference ||
      left.closed_mass_target != right.closed_mass_target ||
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
                       std::vector<std::uint8_t>& bytes) noexcept {
  std::uint64_t size = rank == 0 ? bytes.size() : 0U;
  if (MPI_Bcast(&size, 1, MPI_UINT64_T, 0, communicator) != MPI_SUCCESS ||
      size == 0U || size > static_cast<std::uint64_t>(INT_MAX))
    return {StatusCode::mpi_failure, kRestartCollective};
  try {
    if (rank != 0) bytes.resize(static_cast<std::size_t>(size));
  } catch (...) {
    return {StatusCode::allocation_failure, kRestartCollective};
  }
  return MPI_Bcast(bytes.data(), static_cast<int>(bytes.size()), MPI_BYTE, 0,
                   communicator) == MPI_SUCCESS
             ? Status{}
             : Status{StatusCode::mpi_failure, kRestartCollective};
}

Status broadcast_string(MPI_Comm communicator, int rank,
                        std::string& value) noexcept {
  std::vector<std::uint8_t> bytes;
  try {
    if (rank == 0) bytes.assign(value.begin(), value.end());
  } catch (...) {
    return {StatusCode::allocation_failure, kRestartCollective};
  }
  Status status = broadcast_bytes(communicator, rank, bytes);
  if (status && rank != 0) {
    try {
      value.assign(bytes.begin(), bytes.end());
    } catch (...) {
      status = {StatusCode::allocation_failure, kRestartCollective};
    }
  }
  return collective_status(communicator, status);
}

Status read_current_name(const fs::path& directory,
                         std::string& out) noexcept {
  std::vector<std::uint8_t> bytes;
  if (!read_file(directory / "current", bytes) || bytes.empty() ||
      bytes.size() > 256U)
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
}

Status validate_expected(const RestartExpected& expected,
                         const Manifest& manifest) noexcept {
  if (!valid_global_patch(expected.global_cells, expected.target_patch) ||
      !same(expected.global_cells, manifest.global_cells) ||
      expected.plan == 0U || expected.schema == 0U ||
      expected.geometry == 0U || expected.plan != manifest.plan ||
      expected.schema != manifest.schema ||
      expected.geometry != manifest.geometry || expected.fields.data == nullptr ||
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
  if (manifest.format_version == kExactHistoryFormatVersion) {
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

void prune_generations(const fs::path& directory, std::uint32_t keep_last,
                       const std::string& current) noexcept {
  std::error_code error;
  std::vector<fs::path> generations;
  for (fs::directory_iterator iterator(directory, error), end;
       !error && iterator != end; iterator.increment(error)) {
    const std::string name = iterator->path().filename().string();
    if (iterator->is_directory(error) && !error &&
        name.rfind("generation-", 0U) == 0U &&
        name.find("-pending") == std::string::npos) {
      generations.push_back(iterator->path());
    }
  }
  if (error) return;
  std::sort(generations.begin(), generations.end());
  while (generations.size() > keep_last) {
    auto selected = generations.begin();
    while (selected != generations.end() &&
           selected->filename().string() == current)
      ++selected;
    if (selected == generations.end()) break;
    fs::remove_all(*selected, error);
    if (error) return;
    generations.erase(selected);
  }
  (void)sync_directory(directory);
}

bool remove_stale_pending(const fs::path& directory) noexcept {
  std::error_code error;
  for (fs::directory_iterator iterator(directory, error), end;
       !error && iterator != end; iterator.increment(error)) {
    const std::string name = iterator->path().filename().string();
    if (name.size() >= 8U &&
        name.compare(name.size() - 8U, 8U, "-pending") == 0) {
      fs::remove_all(iterator->path(), error);
      if (error) return false;
    }
  }
  return !error && sync_directory(directory);
}

Status publish_generation(const fs::path& directory,
                          const fs::path& pending,
                          const std::string& generation,
                          std::uint32_t keep_last) noexcept {
  std::error_code error;
  const fs::path final = directory / generation;
  fs::rename(pending, final, error);
  if (error || !sync_directory(directory))
    return {StatusCode::io_failure, kRestartPublication};
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (injected(detail::RestartFailurePoint::after_generation_rename, 0))
    return {StatusCode::io_failure, kRestartPublication};
#endif
  const std::string pointer_text = generation + "\n";
  std::vector<std::uint8_t> pointer(pointer_text.begin(), pointer_text.end());
  const fs::path pointer_pending = directory / ("current-" + generation + "-pending");
  if (!write_file_sync(pointer_pending, pointer))
    return {StatusCode::io_failure, kRestartPublication};
  fs::rename(pointer_pending, directory / "current", error);
  if (error || !sync_directory(directory))
    return {StatusCode::io_failure, kRestartPublication};
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (injected(detail::RestartFailurePoint::after_current_switch, 0))
    return {StatusCode::io_failure, kRestartPublication};
#endif
  prune_generations(directory, keep_last, generation);
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
}

Status RestartWriter::write(MPI_Comm communicator,
                            const std::filesystem::path& restart_directory,
                            const RestartSnapshot& snapshot,
                            RestartWriteOptions options) noexcept try {
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0 ||
      restart_directory.empty() || options.keep_last == 0U) {
    return {StatusCode::invalid_plan, kRestartInput};
  }
  Status status = valid_snapshot(snapshot)
                      ? Status{}
                      : Status{StatusCode::invalid_plan, kRestartInput};
  status = collective_status(communicator, status);
  if (!status) return status;
  status = consensus_u64(communicator, snapshot_signature(snapshot));
  if (!status) return status;

  std::string generation;
  if (rank == 0) {
    const auto tick = static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    generation = "generation-" + std::to_string(snapshot.step) + "-" +
                 std::to_string(tick);
  }
  status = broadcast_string(communicator, rank, generation);
  if (!status) return status;
  const fs::path pending = restart_directory / (generation + "-pending");
  if (rank == 0) {
    std::error_code error;
    fs::create_directories(restart_directory, error);
    if (!error && !remove_stale_pending(restart_directory))
      error = std::make_error_code(std::errc::io_error);
    if (!error) fs::create_directory(pending, error);
    status = !error && sync_directory(restart_directory)
                 ? Status{}
                 : Status{StatusCode::io_failure, kRestartDirectory};
  }
  status = collective_status(communicator, status);
  if (!status) return status;
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (injected(detail::RestartFailurePoint::after_directory, rank))
    status = {StatusCode::io_failure, kRestartDirectory};
  status = collective_status(communicator, status);
  if (!status) return status;
#endif

  std::vector<std::uint8_t> rank_bytes;
  status = encode_rank_block(snapshot, size, rank, rank_bytes);
  if (status &&
      !write_file_sync(pending / rank_name(static_cast<std::uint32_t>(rank)),
                       rank_bytes)) {
    status = {StatusCode::io_failure, kRestartRankFile};
  }
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
  if (status && injected(detail::RestartFailurePoint::after_rank_file, rank))
    status = {StatusCode::io_failure, kRestartRankFile};
#endif
  status = collective_status(communicator, status);
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
  std::vector<std::uint64_t> gathered;
  if (rank == 0) gathered.resize(static_cast<std::size_t>(size) * 8U);
  if (MPI_Gather(local_record.data(), static_cast<int>(local_record.size()),
                 MPI_UINT64_T, rank == 0 ? gathered.data() : nullptr,
                 static_cast<int>(local_record.size()), MPI_UINT64_T, 0,
                 communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kRestartCollective};
  }
  if (rank == 0) {
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
    status = encode_manifest(snapshot, size, records, manifest_bytes);
    if (status && !write_file_sync(pending / "manifest.bin", manifest_bytes))
      status = {StatusCode::io_failure, kRestartManifest};
    for (int source = 0; source < size && status; ++source) {
      std::vector<std::uint8_t> verify;
      const RankRecord& record = records[static_cast<std::size_t>(source)];
      if (!read_file(pending / rank_name(static_cast<std::uint32_t>(source)),
                     verify) ||
          verify.size() != record.bytes ||
          hash_bytes(verify.data(), verify.size()) != record.hash ||
          !verified_integrity(verify)) {
        status = {StatusCode::io_failure, kRestartIntegrity};
      }
    }
    if (status && !sync_directory(pending))
      status = {StatusCode::io_failure, kRestartDirectory};
#ifdef HUNDUN_V04_ENABLE_TEST_ACCESS
    if (status &&
        injected(detail::RestartFailurePoint::after_manifest, rank))
      status = {StatusCode::io_failure, kRestartManifest};
#endif
  }
  status = collective_status(communicator, status);
  if (!status) return status;
  if (rank == 0)
    status = publish_generation(restart_directory, pending, generation,
                                options.keep_last);
  return collective_status(communicator, status);
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kRestartInput};
} catch (...) {
  return {StatusCode::io_failure, kRestartInput};
}

Status RestartReader::load(MPI_Comm communicator,
                           const std::filesystem::path& restart_directory,
                           const RestartExpected& expected,
                           RestartImage& out) noexcept try {
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0 ||
      restart_directory.empty()) {
    return {StatusCode::invalid_plan, kRestartInput};
  }
  std::string generation;
  std::vector<std::uint8_t> manifest_bytes;
  Status status;
  if (rank == 0) {
    status = read_current_name(restart_directory, generation);
    if (status &&
        !read_file(restart_directory / generation / "manifest.bin",
                   manifest_bytes))
      status = {StatusCode::io_failure, kRestartManifest};
  }
  status = collective_status(communicator, status);
  if (!status) return status;
  status = broadcast_string(communicator, rank, generation);
  if (!status) return status;
  status = broadcast_bytes(communicator, rank, manifest_bytes);
  status = collective_status(communicator, status);
  if (!status) return status;
  Manifest manifest;
  status = parse_manifest(manifest_bytes, manifest);
  if (status) status = validate_expected(expected, manifest);
  status = collective_status(communicator, status);
  if (!status) return status;

  const fs::path generation_directory = restart_directory / generation;
  for (std::uint32_t source = static_cast<std::uint32_t>(rank);
       source < manifest.rank_count && status;
       source += static_cast<std::uint32_t>(size)) {
    std::vector<std::uint8_t> bytes;
    RankBlock block;
    const RankRecord record = manifest.ranks[source];
    if (!read_file(generation_directory / rank_name(source), bytes) ||
        bytes.size() != record.bytes ||
        hash_bytes(bytes.data(), bytes.size()) != record.hash) {
      status = {StatusCode::io_failure, kRestartIntegrity};
      break;
    }
    status = parse_rank_block(bytes, block);
    if (status &&
        (block.rank_count != manifest.rank_count || block.rank != source ||
         !same(block.patch.begin, record.begin) ||
         !same(block.patch.cells, record.cells) ||
         !same_common(block.common, manifest))) {
      status = {StatusCode::io_failure, kRestartMismatch};
    }
  }
  status = collective_status(communicator, status);
  if (!status) return status;

  RestartImage candidate;
  candidate.global_cells = manifest.global_cells;
  candidate.patch = expected.target_patch;
  candidate.plan = manifest.plan;
  candidate.schema = manifest.schema;
  candidate.geometry = manifest.geometry;
  candidate.time = manifest.time;
  candidate.dt = manifest.dt;
  candidate.pressure_reference = manifest.pressure_reference;
  candidate.step = manifest.step;
  candidate.controller_state = manifest.controller_state;
  if (!detail::runtime_sha256_bytes(
          {manifest_bytes.data(), manifest_bytes.size()},
          candidate.source_manifest_sha256))
    status = {StatusCode::io_failure, kRestartIntegrity};
  const bool exact_history =
      manifest.format_version == kExactHistoryFormatVersion;
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
  std::array<Int3, 3U> target_face_extents{{
      {expected.target_patch.cells.x + 1, expected.target_patch.cells.y,
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
  for (std::uint32_t source = 0U;
       source < manifest.rank_count && status; ++source) {
    std::vector<std::uint8_t> bytes;
    RankBlock block;
    if (!read_file(generation_directory / rank_name(source), bytes)) {
      status = {StatusCode::io_failure, kRestartRankFile};
      break;
    }
    status = parse_rank_block(bytes, block);
    if (!status) break;
    const Int3 begin{
        std::max(block.patch.begin.x, expected.target_patch.begin.x),
        std::max(block.patch.begin.y, expected.target_patch.begin.y),
        std::max(block.patch.begin.z, expected.target_patch.begin.z)};
    const Int3 end{
        std::min(block.patch.begin.x + block.patch.cells.x,
                 expected.target_patch.begin.x + expected.target_patch.cells.x),
        std::min(block.patch.begin.y + block.patch.cells.y,
                 expected.target_patch.begin.y + expected.target_patch.cells.y),
        std::min(block.patch.begin.z + block.patch.cells.z,
                 expected.target_patch.begin.z + expected.target_patch.cells.z)};
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
          const auto copy_fields = [&](std::vector<RestartImageField>& target,
                                       const std::vector<RestartImageField>&
                                           source_fields) {
            for (std::size_t field_index = 0U;
                 field_index < target.size(); ++field_index) {
              const std::size_t components = target[field_index].components;
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
          axis == 0U ? block.patch.begin.x + block.patch.cells.x
                     : (axis == 1U
                            ? block.patch.begin.y + block.patch.cells.y
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
            const Int3 target_local{
                global.x - expected.target_patch.begin.x,
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
  if (status &&
      std::any_of(cell_coverage.begin(), cell_coverage.end(),
                  [](std::uint8_t value) { return value != 1U; }))
    status = {StatusCode::io_failure, kRestartCoverage};
  for (std::size_t axis = 0U; axis < 3U && status; ++axis) {
    if (std::any_of(face_coverage[axis].begin(), face_coverage[axis].end(),
                    [](std::uint8_t value) { return value != 1U; }))
      status = {StatusCode::io_failure, kRestartCoverage};
  }
  status = collective_status(communicator, status);
  if (!status) return status;
  out = std::move(candidate);
  return {};
} catch (const std::bad_alloc&) {
  return {StatusCode::allocation_failure, kRestartInput};
} catch (...) {
  return {StatusCode::io_failure, kRestartInput};
}

}  // namespace hundun::v04
