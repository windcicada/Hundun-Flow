// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/rt_restart_binary.hpp"

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "rt_mpi_error_detail.hpp"
#include "rt_restart_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <new>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::runtime {
namespace {

constexpr std::array<char, 8> kRankMagic{'H', 'U', 'N', 'D',
                                         'U', 'N', 'R', '1'};
constexpr std::array<char, 8> kManifestMagic{'H', 'U', 'N', 'D',
                                             'M', 'A', 'N', '1'};
constexpr std::array<char, 8> kCompletedMagic{'H', 'U', 'N', 'D',
                                              'C', 'M', 'P', '1'};
constexpr std::uint32_t kVersion = 1U;
constexpr std::uint32_t kEndian = UINT32_C(0x01020304);

static_assert(sizeof(std::int32_t) == 4U);
static_assert(sizeof(std::uint32_t) == 4U);
static_assert(sizeof(std::int64_t) == 8U);
static_assert(sizeof(std::uint64_t) == 8U);
static_assert(sizeof(double) == 8U);
static_assert(std::numeric_limits<double>::is_iec559);

bool same(Int3 lhs, Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(Box3 lhs, Box3 rhs) noexcept {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

Int3 box_extent(Box3 box) noexcept {
  return Int3{box.end.x - box.begin.x, box.end.y - box.begin.y,
              box.end.z - box.begin.z};
}

std::uint64_t checked_add(std::uint64_t lhs, std::uint64_t rhs,
                          std::string_view operation) {
  if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
    throw Error(std::string(operation) + " overflows uint64");
  }
  return lhs + rhs;
}

std::uint64_t checked_multiply(std::uint64_t lhs, std::uint64_t rhs,
                               std::string_view operation) {
  if (rhs != 0U && lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
    throw Error(std::string(operation) + " overflows uint64");
  }
  return lhs * rhs;
}

std::size_t checked_size(std::uint64_t value, std::string_view operation) {
  if (value > std::numeric_limits<std::size_t>::max()) {
    throw Error(std::string(operation) + " exceeds local size range");
  }
  return static_cast<std::size_t>(value);
}

std::int32_t checked_i32(int value, std::string_view operation) {
  if (value < std::numeric_limits<std::int32_t>::min() ||
      value > std::numeric_limits<std::int32_t>::max()) {
    throw Error(std::string(operation) + " exceeds int32 range");
  }
  return static_cast<std::int32_t>(value);
}

std::uint32_t checked_u32(std::size_t value, std::string_view operation) {
  if (value > std::numeric_limits<std::uint32_t>::max()) {
    throw Error(std::string(operation) + " exceeds uint32 range");
  }
  return static_cast<std::uint32_t>(value);
}

class Encoder final {
public:
  template <class T> void append(T value) {
    static_assert(std::is_trivially_copyable_v<T>);
    append_bytes(&value, sizeof(T));
  }

  void append_bytes(const void *data, std::size_t size) {
    if (size > bytes_.max_size() - bytes_.size()) {
      throw Error("restart byte stream exceeds vector capacity");
    }
    const auto *begin = static_cast<const std::byte *>(data);
    bytes_.insert(bytes_.end(), begin, begin + size);
  }

  void append_magic(const std::array<char, 8> &magic) {
    append_bytes(magic.data(), magic.size());
  }

  std::vector<std::byte> finish() && { return std::move(bytes_); }

private:
  std::vector<std::byte> bytes_;
};

class Decoder final {
public:
  explicit Decoder(const std::vector<std::byte> &bytes) noexcept
      : bytes_(bytes) {}

  template <class T> T read(std::string_view label) {
    static_assert(std::is_trivially_copyable_v<T>);
    require(sizeof(T), label);
    T value{};
    std::memcpy(&value, bytes_.data() + cursor_, sizeof(T));
    cursor_ += sizeof(T);
    return value;
  }

  void require_magic(const std::array<char, 8> &expected,
                     std::string_view label) {
    require(expected.size(), label);
    if (std::memcmp(bytes_.data() + cursor_, expected.data(),
                    expected.size()) != 0) {
      throw Error(std::string(label) + " magic does not match");
    }
    cursor_ += expected.size();
  }

  std::string read_string(std::string_view label) {
    const std::uint32_t length = read<std::uint32_t>(label);
    if (length == 0U) {
      throw Error(std::string(label) + " must not be empty");
    }
    require(length, label);
    const auto *begin = reinterpret_cast<const char *>(bytes_.data() + cursor_);
    std::string result(begin, static_cast<std::size_t>(length));
    cursor_ += length;
    return result;
  }

  std::vector<std::byte> read_bytes(std::size_t size, std::string_view label) {
    require(size, label);
    const auto begin = bytes_.begin() + static_cast<std::ptrdiff_t>(cursor_);
    std::vector<std::byte> result(begin,
                                  begin + static_cast<std::ptrdiff_t>(size));
    cursor_ += size;
    return result;
  }

  void require_end(std::string_view label) const {
    if (cursor_ != bytes_.size()) {
      throw Error(std::string(label) + " has trailing bytes");
    }
  }

private:
  void require(std::size_t size, std::string_view label) const {
    if (cursor_ > bytes_.size() || size > bytes_.size() - cursor_) {
      throw Error(std::string(label) + " is truncated");
    }
  }

  const std::vector<std::byte> &bytes_;
  std::size_t cursor_{};
};

bool valid_utf8(std::string_view text) noexcept {
  std::size_t index = 0;
  while (index < text.size()) {
    const auto lead = static_cast<unsigned char>(text[index]);
    std::size_t continuation_count = 0;
    std::uint32_t code_point = 0;
    if (lead <= 0x7FU) {
      ++index;
      continue;
    }
    if (lead >= 0xC2U && lead <= 0xDFU) {
      continuation_count = 1;
      code_point = lead & 0x1FU;
    } else if (lead >= 0xE0U && lead <= 0xEFU) {
      continuation_count = 2;
      code_point = lead & 0x0FU;
    } else if (lead >= 0xF0U && lead <= 0xF4U) {
      continuation_count = 3;
      code_point = lead & 0x07U;
    } else {
      return false;
    }
    if (continuation_count > text.size() - index - 1U) {
      return false;
    }
    for (std::size_t continuation = 0; continuation < continuation_count;
         ++continuation) {
      const auto byte =
          static_cast<unsigned char>(text[index + 1U + continuation]);
      if ((byte & 0xC0U) != 0x80U) {
        return false;
      }
      code_point = (code_point << 6U) | (byte & 0x3FU);
    }
    if ((continuation_count == 2U && code_point < 0x800U) ||
        (continuation_count == 3U && code_point < 0x10000U) ||
        (code_point >= 0xD800U && code_point <= 0xDFFFU) ||
        code_point > 0x10FFFFU) {
      return false;
    }
    index += continuation_count + 1U;
  }
  return true;
}

std::uint32_t scalar_code(ScalarType type) {
  switch (type) {
  case ScalarType::float64:
    return 1U;
  case ScalarType::int32:
    return 2U;
  case ScalarType::uint8:
    break;
  }
  throw Error("Restart v1 supports persistent Float64 and Int32 fields only");
}

std::size_t scalar_size(ScalarType type) {
  switch (type) {
  case ScalarType::float64:
    return sizeof(double);
  case ScalarType::int32:
    return sizeof(std::int32_t);
  case ScalarType::uint8:
    break;
  }
  throw Error("Restart v1 supports persistent Float64 and Int32 fields only");
}

std::vector<FieldId> persistent_fields(const FieldRegistry &registry) {
  if (!registry.frozen()) {
    throw Error("Restart v1 requires a frozen field registry");
  }
  std::vector<FieldId> fields;
  fields.reserve(registry.size());
  for (std::size_t index = 0; index < registry.size(); ++index) {
    const auto field_id = static_cast<FieldId>(index);
    const FieldDescriptor &descriptor = registry.descriptor(field_id);
    if (descriptor.restart != RestartPolicy::persistent) {
      continue;
    }
    if (descriptor.space != FunctionSpace::cell_average) {
      throw Error("Restart v1 supports persistent cell_average fields only");
    }
    if (descriptor.name.empty() || !valid_utf8(descriptor.name)) {
      throw Error("Restart v1 persistent field names must be nonempty UTF-8");
    }
    static_cast<void>(checked_u32(descriptor.name.size(), "field-name size"));
    static_cast<void>(scalar_code(descriptor.scalar_type));
    if (descriptor.components == 0U) {
      throw Error("Restart v1 field component count must be positive");
    }
    fields.push_back(field_id);
  }
  static_cast<void>(checked_u32(fields.size(), "persistent-field count"));
  return fields;
}

void validate_extent(Int3 extent, std::string_view label) {
  if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0) {
    throw Error(std::string(label) + " must be positive");
  }
}

void validate_box(Box3 box, Int3 global_extent, std::string_view label) {
  validate_extent(global_extent, "global extent");
  const bool valid = box.begin.x >= 0 && box.begin.y >= 0 && box.begin.z >= 0 &&
                     box.end.x > box.begin.x && box.end.y > box.begin.y &&
                     box.end.z > box.begin.z && box.end.x <= global_extent.x &&
                     box.end.y <= global_extent.y &&
                     box.end.z <= global_extent.z;
  if (!valid) {
    throw Error(std::string(label) + " is invalid");
  }
}

void validate_rank_metadata(const detail::RestartRankMetadata &metadata,
                            const FieldStorage &storage) {
  if (metadata.ranks <= 0 || metadata.rank < 0 ||
      metadata.rank >= metadata.ranks) {
    throw Error("Restart v1 rank metadata is invalid");
  }
  if (metadata.step < 0) {
    throw Error("Restart v1 step must be nonnegative");
  }
  if (!std::isfinite(metadata.time_s)) {
    throw Error("Restart v1 time must be finite");
  }
  validate_box(metadata.owned_box, metadata.global_extent,
               "Restart v1 owned box");
  if (!same(storage.interior_extent(), box_extent(metadata.owned_box))) {
    throw Error("Restart v1 storage extent does not match the owned box");
  }
}

template <class T>
void validate_writer_view_layout(const FieldView<const T> &view,
                                 Int3 expected_extent,
                                 const FieldDescriptor &descriptor) {
  if (!same(view.interior_extent(), expected_extent) ||
      view.components() != descriptor.components ||
      view.ghost_width() != descriptor.ghost_width) {
    throw Error("Restart v1 field storage layout does not match the registry");
  }
}

std::string format_step_leaf(std::int64_t step) {
  if (step < 0) {
    throw Error("Restart v1 step must be nonnegative");
  }
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "step" << std::setfill('0') << std::setw(8) << step;
  if (!output.good()) {
    throw Error("unable to format Restart v1 step directory");
  }
  return output.str();
}

std::int64_t parse_step_leaf(const std::filesystem::path &path) {
  const std::string leaf = path.filename().string();
  if (leaf.size() < 12U || leaf.compare(0U, 4U, "step") != 0) {
    throw Error("Restart v1 directory leaf is not a canonical step name");
  }
  std::int64_t step = 0;
  const char *begin = leaf.data() + 4U;
  const char *end = leaf.data() + leaf.size();
  const auto result = std::from_chars(begin, end, step);
  if (result.ec != std::errc{} || result.ptr != end || step < 0 ||
      format_step_leaf(step) != leaf) {
    throw Error("Restart v1 directory leaf is not a canonical step name");
  }
  return step;
}

std::string format_rank_filename(int rank) {
  if (rank < 0) {
    throw Error("Restart v1 rank must be nonnegative");
  }
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << "restart.rank" << std::setfill('0') << std::setw(6) << rank
         << ".bin";
  if (!output.good()) {
    throw Error("unable to format Restart v1 rank filename");
  }
  return output.str();
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

template <class Function>
void converge_phase(const MpiContext &context, std::string_view phase,
                    Function &&function) {
  bool local_ok = true;
  std::string local_message_storage;
  std::string_view local_message;
  const auto set_rich_message = [&](std::string_view detail) noexcept {
    try {
      local_message_storage.reserve(phase.size() + 2U + detail.size());
      local_message_storage.append(phase);
      local_message_storage.append(": ");
      local_message_storage.append(detail);
      local_message = local_message_storage;
    } catch (...) {
      local_message = "Restart v1 failure diagnostic unavailable";
    }
  };
  try {
    function();
  } catch (const std::bad_alloc &) {
    local_ok = false;
    local_message = "Restart v1 allocation failure";
  } catch (const std::exception &error) {
    local_ok = false;
    set_rich_message(error.what());
  } catch (...) {
    local_ok = false;
    set_rich_message("unknown failure");
  }
  const CollectiveStatus status =
      collective_status(context, local_ok, local_message);
  if (!status.ok) {
    throw Error(status.message);
  }
}

template <class T>
std::size_t checked_vector_count(int ranks, std::size_t values_per_rank,
                                 std::string_view operation) {
  if (ranks <= 0) {
    throw Error(std::string(operation) + " requires a positive rank count");
  }
  if (values_per_rank > std::numeric_limits<std::uint64_t>::max()) {
    throw Error(std::string(operation) + " exceeds uint64 range");
  }
  const std::uint64_t count =
      checked_multiply(static_cast<std::uint64_t>(ranks),
                       static_cast<std::uint64_t>(values_per_rank), operation);
  const std::size_t local_count = checked_size(count, operation);
  const std::vector<T> empty;
  if (local_count > empty.max_size()) {
    throw Error(std::string(operation) + " exceeds vector capacity");
  }
  return local_count;
}

void validate_communicators(const MpiContext &context,
                            const StructuredDecomposition &decomposition) {
  if (context.comm() == MPI_COMM_NULL ||
      decomposition.comm() == MPI_COMM_NULL) {
    throw Error("Restart v1 requires live MPI communicators");
  }
  int comparison = MPI_UNEQUAL;
  detail::check_mpi(
      MPI_Comm_compare(context.comm(), decomposition.comm(), &comparison),
      "MPI_Comm_compare");
  if (comparison != MPI_IDENT && comparison != MPI_CONGRUENT) {
    throw Error(
        "Restart v1 decomposition communicator must preserve context ranks");
  }
}

void require_path_agreement(const MpiContext &context,
                            const std::string &local_path,
                            detail::RestartFailureInjection injection = {}) {
  const bool local_size_ok =
      local_path.size() <=
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  const CollectiveStatus size_status = collective_status(
      context, local_size_ok,
      local_size_ok ? "" : "Restart v1 path exceeds MPI count range");
  if (!size_status.ok) {
    throw Error(size_status.message);
  }

  std::uint64_t root_length =
      context.rank() == 0 ? static_cast<std::uint64_t>(local_path.size()) : 0U;
  detail::check_mpi(MPI_Bcast(&root_length, 1, MPI_UINT64_T, 0, context.comm()),
                    "MPI_Bcast");
  std::string root_path;
  converge_phase(context, "Restart v1 path allocation", [&] {
    root_path.resize(checked_size(root_length, "Restart v1 path length"));
    if (context.rank() == 0) {
      root_path = local_path;
    }
    if (injection.phase == detail::RestartFailurePhase::path_preparation &&
        injection.rank == context.rank()) {
      throw std::bad_alloc{};
    }
  });
  detail::check_mpi(MPI_Bcast(root_path.data(), static_cast<int>(root_length),
                              MPI_BYTE, 0, context.comm()),
                    "MPI_Bcast");
  const bool matches = local_path == root_path;
  const CollectiveStatus status = collective_status(
      context, matches,
      matches ? "" : "Restart v1 step-directory paths differ across ranks");
  if (!status.ok) {
    throw Error(status.message);
  }
}

template <std::size_t Size>
void require_int_agreement(const MpiContext &context,
                           const std::array<int, Size> &local,
                           std::string_view message,
                           detail::RestartFailureInjection injection = {}) {
  std::array<int, Size> minimum{};
  std::array<int, Size> maximum{};
  detail::check_mpi(MPI_Allreduce(local.data(), minimum.data(),
                                  static_cast<int>(Size), MPI_INT, MPI_MIN,
                                  context.comm()),
                    "MPI_Allreduce");
  detail::check_mpi(MPI_Allreduce(local.data(), maximum.data(),
                                  static_cast<int>(Size), MPI_INT, MPI_MAX,
                                  context.comm()),
                    "MPI_Allreduce");
  if (injection.phase == detail::RestartFailurePhase::agreement_preparation) {
    converge_phase(context, "Restart v1 agreement preparation", [&] {
      if (injection.rank == context.rank()) {
        throw std::bad_alloc{};
      }
    });
  }
  const bool matches = minimum == maximum;
  const CollectiveStatus status = collective_status(
      context, matches, matches ? std::string_view{} : message);
  if (!status.ok) {
    throw Error(status.message);
  }
}

void require_u64_agreement(const MpiContext &context, std::uint64_t local,
                           std::string_view message,
                           detail::RestartFailureInjection injection = {}) {
  std::uint64_t minimum = 0;
  std::uint64_t maximum = 0;
  detail::check_mpi(
      MPI_Allreduce(&local, &minimum, 1, MPI_UINT64_T, MPI_MIN, context.comm()),
      "MPI_Allreduce");
  detail::check_mpi(
      MPI_Allreduce(&local, &maximum, 1, MPI_UINT64_T, MPI_MAX, context.comm()),
      "MPI_Allreduce");
  if (injection.phase == detail::RestartFailurePhase::agreement_preparation) {
    converge_phase(context, "Restart v1 agreement preparation", [&] {
      if (injection.rank == context.rank()) {
        throw std::bad_alloc{};
      }
    });
  }
  const bool matches = minimum == maximum;
  const CollectiveStatus status = collective_status(
      context, matches, matches ? std::string_view{} : message);
  if (!status.ok) {
    throw Error(status.message);
  }
}

std::vector<Box3> gather_boxes(const MpiContext &context, Box3 local_box,
                               detail::RestartFailureInjection injection = {}) {
  const std::array<int, 6> local{local_box.begin.x, local_box.begin.y,
                                 local_box.begin.z, local_box.end.x,
                                 local_box.end.y,   local_box.end.z};
  std::vector<int> gathered;
  converge_phase(context, "Restart v1 owned-box gather preparation", [&] {
    if (injection.phase == detail::RestartFailurePhase::owned_box_preparation &&
        injection.rank == context.rank()) {
      throw std::bad_alloc{};
    }
    gathered.resize(checked_vector_count<int>(
        context.size(), local.size(), "Restart v1 owned-box gather count"));
  });
  detail::check_mpi(MPI_Allgather(local.data(), static_cast<int>(local.size()),
                                  MPI_INT, gathered.data(),
                                  static_cast<int>(local.size()), MPI_INT,
                                  context.comm()),
                    "MPI_Allgather");
  std::vector<Box3> boxes;
  converge_phase(context, "Restart v1 owned-box result preparation", [&] {
    boxes.reserve(checked_vector_count<Box3>(
        context.size(), 1U, "Restart v1 owned-box result count"));
    for (int process = 0; process < context.size(); ++process) {
      const std::size_t offset =
          static_cast<std::size_t>(process) * local.size();
      boxes.push_back(Box3{
          Int3{gathered[offset], gathered[offset + 1U], gathered[offset + 2U]},
          Int3{gathered[offset + 3U], gathered[offset + 4U],
               gathered[offset + 5U]}});
    }
  });
  return boxes;
}

bool boxes_overlap(Box3 lhs, Box3 rhs) noexcept {
  return lhs.begin.x < rhs.end.x && rhs.begin.x < lhs.end.x &&
         lhs.begin.y < rhs.end.y && rhs.begin.y < lhs.end.y &&
         lhs.begin.z < rhs.end.z && rhs.begin.z < lhs.end.z;
}

void validate_complete_boxes(const std::vector<Box3> &boxes,
                             Int3 global_extent) {
  std::uint64_t total = 0;
  for (std::size_t index = 0; index < boxes.size(); ++index) {
    validate_box(boxes[index], global_extent, "Restart v1 rank owned box");
    const Int3 extent = box_extent(boxes[index]);
    total = checked_add(
        total,
        checked_multiply(checked_multiply(static_cast<std::uint64_t>(extent.x),
                                          static_cast<std::uint64_t>(extent.y),
                                          "owned-box volume"),
                         static_cast<std::uint64_t>(extent.z),
                         "owned-box volume"),
        "owned-box volume sum");
    for (std::size_t other = 0; other < index; ++other) {
      if (boxes_overlap(boxes[index], boxes[other])) {
        throw Error("Restart v1 rank owned boxes overlap");
      }
    }
  }
  const std::uint64_t global_volume = checked_multiply(
      checked_multiply(static_cast<std::uint64_t>(global_extent.x),
                       static_cast<std::uint64_t>(global_extent.y),
                       "global volume"),
      static_cast<std::uint64_t>(global_extent.z), "global volume");
  if (total != global_volume) {
    throw Error("Restart v1 rank owned boxes do not cover the global extent");
  }
}

void require_metadata_agreement(
    const MpiContext &context, const StructuredDecomposition &decomposition,
    std::int64_t step, double time_s, std::uint64_t schema_fingerprint,
    detail::RestartFailureInjection injection = {}) {
  const Int3 global = decomposition.global_extent();
  const Int3 grid = decomposition.process_grid();
  const auto periodic = decomposition.periodic();
  const std::array<int, 10> local{context.size(),
                                  global.x,
                                  global.y,
                                  global.z,
                                  periodic[0] ? 1 : 0,
                                  periodic[1] ? 1 : 0,
                                  periodic[2] ? 1 : 0,
                                  grid.x,
                                  grid.y,
                                  grid.z};
  require_int_agreement(
      context, local, "Restart v1 decomposition metadata differs across ranks",
      injection);
  std::uint64_t step_bits = 0;
  std::memcpy(&step_bits, &step, sizeof(step));
  require_u64_agreement(context, step_bits,
                        "Restart v1 step differs across ranks", injection);
  require_u64_agreement(context, double_bits(time_s),
                        "Restart v1 time bits differ across ranks", injection);
  require_u64_agreement(context, schema_fingerprint,
                        "Restart v1 persistent schemas differ across ranks",
                        injection);
}

std::vector<std::byte> read_file(const std::filesystem::path &path) {
  const std::uintmax_t file_size = std::filesystem::file_size(path);
  if (file_size > std::numeric_limits<std::size_t>::max() ||
      file_size > static_cast<std::uintmax_t>(
                      std::numeric_limits<std::streamsize>::max())) {
    throw Error("Restart v1 file is too large to read");
  }
  std::vector<std::byte> bytes(static_cast<std::size_t>(file_size));
  std::ifstream input(path, std::ios::binary);
  if (!input.is_open()) {
    throw Error("unable to open Restart v1 file: " + path.string());
  }
  if (!bytes.empty()) {
    input.read(reinterpret_cast<char *>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
    if (input.gcount() != static_cast<std::streamsize>(bytes.size())) {
      throw Error("unable to read complete Restart v1 file: " + path.string());
    }
  }
  if (input.bad()) {
    throw Error("Restart v1 file read failed: " + path.string());
  }
  return bytes;
}

void write_atomic_file(const std::filesystem::path &final_path,
                       const std::vector<std::byte> &bytes) {
  std::filesystem::path temporary = final_path;
  temporary += ".tmp";
  if (std::filesystem::exists(final_path) ||
      std::filesystem::exists(temporary)) {
    throw Error("Restart v1 refuses to overwrite an existing file");
  }
  if (bytes.size() >
      static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())) {
    throw Error("Restart v1 output exceeds stream size range");
  }
  std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
  if (!output.is_open()) {
    throw Error("unable to open Restart v1 temporary file: " +
                temporary.string());
  }
  if (!bytes.empty()) {
    output.write(reinterpret_cast<const char *>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
  }
  output.flush();
  if (!output.good()) {
    throw Error("unable to write Restart v1 temporary file: " +
                temporary.string());
  }
  output.close();
  if (output.fail()) {
    throw Error("unable to close Restart v1 temporary file: " +
                temporary.string());
  }
  std::filesystem::rename(temporary, final_path);
}

struct RankRecord final {
  int rank;
  Box3 owned_box;
  std::string filename;
  std::uint64_t byte_size;
  std::uint64_t crc64;
};

struct Manifest final {
  std::int64_t step;
  double time_s;
  int ranks;
  Int3 global_extent;
  std::array<bool, 3> periodic;
  Int3 process_grid;
  std::uint32_t persistent_field_count;
  std::uint64_t schema_fingerprint;
  std::vector<RankRecord> records;
};

std::vector<std::byte> encode_manifest(const Manifest &manifest) {
  if (manifest.ranks <= 0 ||
      manifest.records.size() != static_cast<std::size_t>(manifest.ranks)) {
    throw Error("Restart v1 manifest rank records are incomplete");
  }
  Encoder encoder;
  encoder.append_magic(kManifestMagic);
  encoder.append(kVersion);
  encoder.append(kEndian);
  encoder.append(manifest.step);
  encoder.append(manifest.time_s);
  encoder.append(checked_i32(manifest.ranks, "rank count"));
  encoder.append(checked_i32(manifest.global_extent.x, "global x extent"));
  encoder.append(checked_i32(manifest.global_extent.y, "global y extent"));
  encoder.append(checked_i32(manifest.global_extent.z, "global z extent"));
  for (const bool value : manifest.periodic) {
    encoder.append(static_cast<std::uint8_t>(value ? 1U : 0U));
  }
  encoder.append(checked_i32(manifest.process_grid.x, "process-grid x"));
  encoder.append(checked_i32(manifest.process_grid.y, "process-grid y"));
  encoder.append(checked_i32(manifest.process_grid.z, "process-grid z"));
  encoder.append(manifest.persistent_field_count);
  encoder.append(manifest.schema_fingerprint);
  encoder.append(checked_u32(manifest.records.size(), "rank-record count"));
  for (std::size_t index = 0; index < manifest.records.size(); ++index) {
    const RankRecord &record = manifest.records[index];
    if (record.rank != static_cast<int>(index) ||
        record.filename != format_rank_filename(record.rank)) {
      throw Error("Restart v1 manifest rank records are not canonical");
    }
    validate_box(record.owned_box, manifest.global_extent,
                 "Restart v1 manifest owned box");
    encoder.append(checked_i32(record.rank, "rank ID"));
    encoder.append(checked_i32(record.owned_box.begin.x, "owned begin x"));
    encoder.append(checked_i32(record.owned_box.begin.y, "owned begin y"));
    encoder.append(checked_i32(record.owned_box.begin.z, "owned begin z"));
    encoder.append(checked_i32(record.owned_box.end.x, "owned end x"));
    encoder.append(checked_i32(record.owned_box.end.y, "owned end y"));
    encoder.append(checked_i32(record.owned_box.end.z, "owned end z"));
    encoder.append(checked_u32(record.filename.size(), "rank filename size"));
    encoder.append_bytes(record.filename.data(), record.filename.size());
    encoder.append(record.byte_size);
    encoder.append(record.crc64);
  }
  return std::move(encoder).finish();
}

Manifest decode_manifest(const std::vector<std::byte> &bytes) {
  Decoder decoder(bytes);
  decoder.require_magic(kManifestMagic, "Restart v1 manifest");
  if (decoder.read<std::uint32_t>("manifest version") != kVersion) {
    throw Error("Restart v1 manifest version does not match");
  }
  if (decoder.read<std::uint32_t>("manifest endian") != kEndian) {
    throw Error("Restart v1 manifest endian marker does not match");
  }
  Manifest manifest{};
  manifest.step = decoder.read<std::int64_t>("manifest step");
  manifest.time_s = decoder.read<double>("manifest time");
  manifest.ranks = decoder.read<std::int32_t>("manifest rank count");
  manifest.global_extent =
      Int3{decoder.read<std::int32_t>("manifest global x"),
           decoder.read<std::int32_t>("manifest global y"),
           decoder.read<std::int32_t>("manifest global z")};
  for (std::size_t axis = 0; axis < manifest.periodic.size(); ++axis) {
    const std::uint8_t value = decoder.read<std::uint8_t>("periodic flag");
    if (value > 1U) {
      throw Error("Restart v1 manifest periodic flag is invalid");
    }
    manifest.periodic[axis] = value != 0U;
  }
  manifest.process_grid =
      Int3{decoder.read<std::int32_t>("manifest process x"),
           decoder.read<std::int32_t>("manifest process y"),
           decoder.read<std::int32_t>("manifest process z")};
  manifest.persistent_field_count =
      decoder.read<std::uint32_t>("manifest persistent-field count");
  manifest.schema_fingerprint =
      decoder.read<std::uint64_t>("manifest schema fingerprint");
  const std::uint32_t record_count =
      decoder.read<std::uint32_t>("manifest rank-record count");
  if (manifest.ranks <= 0 ||
      record_count != static_cast<std::uint32_t>(manifest.ranks)) {
    throw Error("Restart v1 manifest rank-record count does not match");
  }
  validate_extent(manifest.global_extent, "manifest global extent");
  validate_extent(manifest.process_grid, "manifest process grid");
  if (manifest.step < 0 || !std::isfinite(manifest.time_s)) {
    throw Error("Restart v1 manifest step/time is invalid");
  }
  manifest.records.reserve(record_count);
  for (std::uint32_t index = 0; index < record_count; ++index) {
    RankRecord record{};
    record.rank = decoder.read<std::int32_t>("manifest rank ID");
    record.owned_box =
        Box3{Int3{decoder.read<std::int32_t>("manifest owned begin x"),
                  decoder.read<std::int32_t>("manifest owned begin y"),
                  decoder.read<std::int32_t>("manifest owned begin z")},
             Int3{decoder.read<std::int32_t>("manifest owned end x"),
                  decoder.read<std::int32_t>("manifest owned end y"),
                  decoder.read<std::int32_t>("manifest owned end z")}};
    record.filename = decoder.read_string("manifest rank filename");
    record.byte_size = decoder.read<std::uint64_t>("manifest rank byte size");
    record.crc64 = decoder.read<std::uint64_t>("manifest rank CRC");
    if (record.rank != static_cast<int>(index) ||
        record.filename != format_rank_filename(record.rank)) {
      throw Error(
          "Restart v1 manifest rank IDs or filenames are not canonical");
    }
    validate_box(record.owned_box, manifest.global_extent,
                 "Restart v1 manifest owned box");
    manifest.records.push_back(std::move(record));
  }
  decoder.require_end("Restart v1 manifest");
  validate_complete_boxes(
      [&] {
        std::vector<Box3> boxes;
        boxes.reserve(manifest.records.size());
        for (const RankRecord &record : manifest.records) {
          boxes.push_back(record.owned_box);
        }
        return boxes;
      }(),
      manifest.global_extent);
  return manifest;
}

std::vector<std::byte> encode_marker(std::uint64_t manifest_size,
                                     std::uint64_t manifest_crc) {
  Encoder encoder;
  encoder.append_magic(kCompletedMagic);
  encoder.append(kVersion);
  encoder.append(kEndian);
  encoder.append(manifest_size);
  encoder.append(manifest_crc);
  return std::move(encoder).finish();
}

std::pair<std::uint64_t, std::uint64_t>
decode_marker(const std::vector<std::byte> &bytes) {
  Decoder decoder(bytes);
  decoder.require_magic(kCompletedMagic, "Restart v1 completion marker");
  if (decoder.read<std::uint32_t>("completion version") != kVersion) {
    throw Error("Restart v1 completion marker version does not match");
  }
  if (decoder.read<std::uint32_t>("completion endian") != kEndian) {
    throw Error("Restart v1 completion marker endian does not match");
  }
  const std::uint64_t size =
      decoder.read<std::uint64_t>("completion manifest size");
  const std::uint64_t crc =
      decoder.read<std::uint64_t>("completion manifest CRC");
  decoder.require_end("Restart v1 completion marker");
  return {size, crc};
}

void validate_checkpoint_preconditions(
    const MpiContext &context, const StructuredDecomposition &decomposition,
    const FieldRegistry &registry, const FieldStorage &storage,
    const std::filesystem::path &step_directory, std::int64_t step,
    double time_s, detail::RestartFailureInjection injection,
    std::string &path_text, std::uint64_t &schema_fingerprint,
    std::vector<std::byte> &rank_bytes) {
  validate_communicators(context, decomposition);
  if (step < 0 || !std::isfinite(time_s)) {
    throw Error("Restart v1 step/time is invalid");
  }
  if (step_directory.filename().string() != format_step_leaf(step)) {
    throw Error("Restart v1 directory leaf does not match the step");
  }
  path_text = step_directory.string();
  if (path_text.empty()) {
    throw Error("Restart v1 step directory must not be empty");
  }
  schema_fingerprint = detail::restart_schema_fingerprint(registry);
  if (!same(storage.interior_extent(), decomposition.local_extent())) {
    throw Error("Restart v1 storage extent does not match the decomposition");
  }
  if (injection.phase != detail::RestartFailurePhase::none &&
      (injection.rank < 0 || injection.rank >= context.size())) {
    throw Error("Restart v1 test failure injection rank is invalid");
  }
  const detail::RestartRankMetadata metadata{context.rank(),
                                             context.size(),
                                             decomposition.global_extent(),
                                             decomposition.owned_box(),
                                             step,
                                             time_s};
  rank_bytes = detail::encode_restart_rank(metadata, registry, storage);
}

} // namespace

namespace detail {

std::uint64_t crc64_ecma(const std::byte *data, std::size_t size) noexcept {
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  std::uint64_t crc = 0;
  for (std::size_t index = 0; index < size; ++index) {
    crc ^=
        static_cast<std::uint64_t>(std::to_integer<unsigned char>(data[index]))
        << 56U;
    for (int bit = 0; bit < 8; ++bit) {
      crc = (crc & (UINT64_C(1) << 63U)) != 0U ? (crc << 1U) ^ polynomial
                                               : crc << 1U;
    }
  }
  return crc;
}

std::uint64_t restart_schema_fingerprint(const FieldRegistry &registry) {
  const auto fields = persistent_fields(registry);
  Encoder encoder;
  encoder.append(checked_u32(fields.size(), "persistent-field count"));
  for (const FieldId field_id : fields) {
    const FieldDescriptor &descriptor = registry.descriptor(field_id);
    encoder.append(checked_u32(descriptor.name.size(), "field-name size"));
    encoder.append_bytes(descriptor.name.data(), descriptor.name.size());
    encoder.append(scalar_code(descriptor.scalar_type));
    encoder.append(descriptor.components);
  }
  const auto bytes = std::move(encoder).finish();
  return crc64_ecma(bytes.data(), bytes.size());
}

std::uint64_t checked_owned_value_count(Int3 extent, std::uint32_t components) {
  validate_extent(extent, "Restart v1 owned extent");
  if (components == 0U) {
    throw Error("Restart v1 component count must be positive");
  }
  return checked_multiply(
      checked_multiply(checked_multiply(static_cast<std::uint64_t>(extent.x),
                                        static_cast<std::uint64_t>(extent.y),
                                        "Restart v1 value count"),
                       static_cast<std::uint64_t>(extent.z),
                       "Restart v1 value count"),
      components, "Restart v1 value count");
}

std::vector<std::byte> encode_restart_rank(RestartRankMetadata metadata,
                                           const FieldRegistry &registry,
                                           const FieldStorage &storage) {
  validate_rank_metadata(metadata, storage);
  const auto fields = persistent_fields(registry);
  const Int3 local_extent = box_extent(metadata.owned_box);
  Encoder encoder;
  encoder.append_magic(kRankMagic);
  encoder.append(kVersion);
  encoder.append(kEndian);
  encoder.append(checked_i32(metadata.rank, "rank"));
  encoder.append(checked_i32(metadata.ranks, "rank count"));
  encoder.append(checked_i32(metadata.global_extent.x, "global x extent"));
  encoder.append(checked_i32(metadata.global_extent.y, "global y extent"));
  encoder.append(checked_i32(metadata.global_extent.z, "global z extent"));
  encoder.append(checked_i32(metadata.owned_box.begin.x, "owned begin x"));
  encoder.append(checked_i32(metadata.owned_box.begin.y, "owned begin y"));
  encoder.append(checked_i32(metadata.owned_box.begin.z, "owned begin z"));
  encoder.append(checked_i32(local_extent.x, "local x extent"));
  encoder.append(checked_i32(local_extent.y, "local y extent"));
  encoder.append(checked_i32(local_extent.z, "local z extent"));
  encoder.append(metadata.step);
  encoder.append(metadata.time_s);
  encoder.append(checked_u32(fields.size(), "persistent-field count"));
  for (const FieldId field_id : fields) {
    const FieldDescriptor &descriptor = registry.descriptor(field_id);
    const std::uint64_t value_count =
        checked_owned_value_count(local_extent, descriptor.components);
    encoder.append(checked_u32(descriptor.name.size(), "field-name size"));
    encoder.append_bytes(descriptor.name.data(), descriptor.name.size());
    encoder.append(scalar_code(descriptor.scalar_type));
    encoder.append(descriptor.components);
    encoder.append(value_count);
    if (descriptor.scalar_type == ScalarType::float64) {
      const auto view = storage.view<double>(field_id);
      validate_writer_view_layout(view, local_extent, descriptor);
      for (int k = 0; k < local_extent.z; ++k) {
        for (int j = 0; j < local_extent.y; ++j) {
          for (int i = 0; i < local_extent.x; ++i) {
            for (std::uint32_t component = 0; component < descriptor.components;
                 ++component) {
              encoder.append(view(i, j, k, static_cast<int>(component)));
            }
          }
        }
      }
    } else {
      const auto view = storage.view<std::int32_t>(field_id);
      validate_writer_view_layout(view, local_extent, descriptor);
      for (int k = 0; k < local_extent.z; ++k) {
        for (int j = 0; j < local_extent.y; ++j) {
          for (int i = 0; i < local_extent.x; ++i) {
            for (std::uint32_t component = 0; component < descriptor.components;
                 ++component) {
              encoder.append(view(i, j, k, static_cast<int>(component)));
            }
          }
        }
      }
    }
  }
  return std::move(encoder).finish();
}

StagedRestartRank decode_restart_rank(const std::vector<std::byte> &bytes,
                                      RestartRankMetadata expected,
                                      const FieldRegistry &registry) {
  const auto fields = persistent_fields(registry);
  Decoder decoder(bytes);
  decoder.require_magic(kRankMagic, "Restart v1 rank file");
  if (decoder.read<std::uint32_t>("rank-file version") != kVersion) {
    throw Error("Restart v1 rank-file version does not match");
  }
  if (decoder.read<std::uint32_t>("rank-file endian") != kEndian) {
    throw Error("Restart v1 rank-file endian marker does not match");
  }
  const int rank = decoder.read<std::int32_t>("rank-file rank");
  const int ranks = decoder.read<std::int32_t>("rank-file rank count");
  const Int3 global_extent{decoder.read<std::int32_t>("rank-file global x"),
                           decoder.read<std::int32_t>("rank-file global y"),
                           decoder.read<std::int32_t>("rank-file global z")};
  const Int3 begin{decoder.read<std::int32_t>("rank-file begin x"),
                   decoder.read<std::int32_t>("rank-file begin y"),
                   decoder.read<std::int32_t>("rank-file begin z")};
  const Int3 local_extent{decoder.read<std::int32_t>("rank-file local x"),
                          decoder.read<std::int32_t>("rank-file local y"),
                          decoder.read<std::int32_t>("rank-file local z")};
  const std::int64_t step = decoder.read<std::int64_t>("rank-file step");
  const double time_s = decoder.read<double>("rank-file time");
  const std::uint32_t field_count =
      decoder.read<std::uint32_t>("rank-file field count");
  const Int3 expected_local_extent = box_extent(expected.owned_box);
  if (rank != expected.rank || ranks != expected.ranks ||
      !same(global_extent, expected.global_extent) ||
      !same(begin, expected.owned_box.begin) ||
      !same(local_extent, expected_local_extent) || step != expected.step ||
      double_bits(time_s) != double_bits(expected.time_s)) {
    throw Error("Restart v1 rank-file metadata does not match");
  }
  if (field_count != fields.size()) {
    throw Error("Restart v1 rank-file field count does not match");
  }

  StagedRestartRank staged{expected, {}};
  staged.fields.reserve(fields.size());
  for (const FieldId field_id : fields) {
    const FieldDescriptor &descriptor = registry.descriptor(field_id);
    const std::string name = decoder.read_string("rank-file field name");
    const std::uint32_t type_code =
        decoder.read<std::uint32_t>("rank-file scalar type");
    const std::uint32_t components =
        decoder.read<std::uint32_t>("rank-file components");
    const std::uint64_t value_count =
        decoder.read<std::uint64_t>("rank-file value count");
    const std::uint64_t expected_count =
        checked_owned_value_count(local_extent, descriptor.components);
    if (name != descriptor.name || !valid_utf8(name) ||
        type_code != scalar_code(descriptor.scalar_type) ||
        components != descriptor.components || value_count != expected_count) {
      throw Error("Restart v1 rank-file field schema does not match");
    }
    const std::uint64_t byte_count =
        checked_multiply(value_count, scalar_size(descriptor.scalar_type),
                         "Restart v1 field payload size");
    staged.fields.push_back(StagedRestartField{
        field_id, descriptor.scalar_type, descriptor.components,
        decoder.read_bytes(checked_size(byte_count, "field payload size"),
                           "rank-file field payload")});
  }
  decoder.require_end("Restart v1 rank file");
  return staged;
}

RestartCommitPlan make_restart_commit_plan(const FieldRegistry &registry,
                                           FieldStorage &storage,
                                           Int3 expected_extent) {
  validate_extent(expected_extent, "Restart v1 destination extent");
  if (!same(storage.interior_extent(), expected_extent)) {
    throw Error("Restart v1 storage extent does not match the decomposition");
  }

  const auto fields = persistent_fields(registry);
  RestartCommitPlan plan{expected_extent, {}};
  plan.fields.reserve(fields.size());
  for (const FieldId field_id : fields) {
    const FieldDescriptor &descriptor = registry.descriptor(field_id);
    if (descriptor.scalar_type == ScalarType::float64) {
      auto view = storage.view<double>(field_id);
      if (!same(view.interior_extent(), expected_extent) ||
          view.components() != descriptor.components ||
          view.ghost_width() != descriptor.ghost_width) {
        throw Error(
            "Restart v1 destination storage does not match the registry");
      }
      plan.fields.push_back(RestartCommitField{field_id, descriptor.scalar_type,
                                               descriptor.components,
                                               RestartCommitView{view}});
    } else {
      auto view = storage.view<std::int32_t>(field_id);
      if (!same(view.interior_extent(), expected_extent) ||
          view.components() != descriptor.components ||
          view.ghost_width() != descriptor.ghost_width) {
        throw Error(
            "Restart v1 destination storage does not match the registry");
      }
      plan.fields.push_back(RestartCommitField{field_id, descriptor.scalar_type,
                                               descriptor.components,
                                               RestartCommitView{view}});
    }
  }
  return plan;
}

void commit_restart_rank(const StagedRestartRank &staged,
                         const RestartCommitPlan &plan) noexcept {
  const Int3 extent = plan.extent;
  for (std::size_t field_index = 0; field_index < staged.fields.size();
       ++field_index) {
    const StagedRestartField &staged_field = staged.fields[field_index];
    const RestartCommitField &commit_field = plan.fields[field_index];
    std::size_t cursor = 0;
    if (commit_field.scalar_type == ScalarType::float64) {
      auto *view = std::get_if<FieldView<double>>(&commit_field.view);
      for (int k = 0; k < extent.z; ++k) {
        for (int j = 0; j < extent.y; ++j) {
          for (int i = 0; i < extent.x; ++i) {
            for (std::uint32_t component = 0;
                 component < commit_field.components; ++component) {
              double value = 0.0;
              std::memcpy(&value, staged_field.values.data() + cursor,
                          sizeof(value));
              cursor += sizeof(value);
              (*view)(i, j, k, static_cast<int>(component)) = value;
            }
          }
        }
      }
    } else {
      auto *view = std::get_if<FieldView<std::int32_t>>(&commit_field.view);
      for (int k = 0; k < extent.z; ++k) {
        for (int j = 0; j < extent.y; ++j) {
          for (int i = 0; i < extent.x; ++i) {
            for (std::uint32_t component = 0;
                 component < commit_field.components; ++component) {
              std::int32_t value = 0;
              std::memcpy(&value, staged_field.values.data() + cursor,
                          sizeof(value));
              cursor += sizeof(value);
              (*view)(i, j, k, static_cast<int>(component)) = value;
            }
          }
        }
      }
    }
  }
}

void write_restart_checkpoint_with_failure(
    const MpiContext &context, const StructuredDecomposition &decomposition,
    const FieldRegistry &registry, const FieldStorage &storage,
    const std::filesystem::path &step_directory, std::int64_t step,
    double time_s, RestartFailureInjection injection) {
  require_mpi_active("write Restart v1 checkpoint");
  std::string path_text;
  std::uint64_t schema_fingerprint = 0;
  std::vector<std::byte> rank_bytes;
  converge_phase(context, "Restart v1 preflight", [&] {
    validate_checkpoint_preconditions(
        context, decomposition, registry, storage, step_directory, step, time_s,
        injection, path_text, schema_fingerprint, rank_bytes);
  });
  require_path_agreement(context, path_text, injection);
  require_metadata_agreement(context, decomposition, step, time_s,
                             schema_fingerprint, injection);
  const auto boxes =
      gather_boxes(context, decomposition.owned_box(), injection);
  converge_phase(context, "Restart v1 owned-box validation", [&] {
    validate_complete_boxes(boxes, decomposition.global_extent());
  });
  converge_phase(context, "Restart v1 existing-directory check", [&] {
    if (context.rank() == 0 && std::filesystem::exists(step_directory)) {
      throw Error("Restart v1 step directory already exists");
    }
  });

  converge_phase(context, "Restart v1 directory creation", [&] {
    if (context.rank() == 0 &&
        !std::filesystem::create_directories(step_directory)) {
      throw Error("unable to create Restart v1 step directory");
    }
  });

  std::string local_filename;
  converge_phase(context, "Restart v1 rank-file write", [&] {
    local_filename = format_rank_filename(context.rank());
    if (injection.phase == RestartFailurePhase::filename_preparation &&
        injection.rank == context.rank()) {
      throw std::bad_alloc{};
    }
    if (injection.phase == RestartFailurePhase::rank_file &&
        injection.rank == context.rank()) {
      throw Error("injected Restart v1 rank-file failure");
    }
    write_atomic_file(step_directory / local_filename, rank_bytes);
  });

  const std::array<int, 7> local_record{context.rank(),
                                        decomposition.owned_box().begin.x,
                                        decomposition.owned_box().begin.y,
                                        decomposition.owned_box().begin.z,
                                        decomposition.owned_box().end.x,
                                        decomposition.owned_box().end.y,
                                        decomposition.owned_box().end.z};
  const std::array<std::uint64_t, 2> local_integrity{
      static_cast<std::uint64_t>(rank_bytes.size()),
      crc64_ecma(rank_bytes.data(), rank_bytes.size())};
  std::vector<int> gathered_records;
  std::vector<std::uint64_t> gathered_integrity;
  converge_phase(context, "Restart v1 record-gather preparation", [&] {
    if (injection.phase == RestartFailurePhase::record_gather_preparation &&
        injection.rank == context.rank()) {
      throw std::bad_alloc{};
    }
    const std::size_t record_count = checked_vector_count<int>(
        context.size(), local_record.size(), "Restart v1 record-gather count");
    const std::size_t integrity_count = checked_vector_count<std::uint64_t>(
        context.size(), local_integrity.size(),
        "Restart v1 integrity-gather count");
    if (context.rank() == 0) {
      gathered_records.resize(record_count);
      gathered_integrity.resize(integrity_count);
    }
  });
  detail::check_mpi(
      MPI_Gather(
          local_record.data(), static_cast<int>(local_record.size()), MPI_INT,
          context.rank() == 0 ? gathered_records.data() : nullptr,
          static_cast<int>(local_record.size()), MPI_INT, 0, context.comm()),
      "MPI_Gather");
  detail::check_mpi(
      MPI_Gather(local_integrity.data(),
                 static_cast<int>(local_integrity.size()), MPI_UINT64_T,
                 context.rank() == 0 ? gathered_integrity.data() : nullptr,
                 static_cast<int>(local_integrity.size()), MPI_UINT64_T, 0,
                 context.comm()),
      "MPI_Gather");

  std::vector<std::byte> manifest_bytes;
  converge_phase(context, "Restart v1 manifest write", [&] {
    if (context.rank() != 0) {
      return;
    }
    if (injection.phase == RestartFailurePhase::manifest &&
        injection.rank == 0) {
      throw Error("injected Restart v1 manifest failure");
    }
    std::vector<RankRecord> records;
    records.reserve(static_cast<std::size_t>(context.size()));
    for (int process = 0; process < context.size(); ++process) {
      const std::size_t record_offset =
          static_cast<std::size_t>(process) * local_record.size();
      const std::size_t integrity_offset =
          static_cast<std::size_t>(process) * local_integrity.size();
      records.push_back(RankRecord{
          gathered_records[record_offset],
          Box3{Int3{gathered_records[record_offset + 1U],
                    gathered_records[record_offset + 2U],
                    gathered_records[record_offset + 3U]},
               Int3{gathered_records[record_offset + 4U],
                    gathered_records[record_offset + 5U],
                    gathered_records[record_offset + 6U]}},
          format_rank_filename(process), gathered_integrity[integrity_offset],
          gathered_integrity[integrity_offset + 1U]});
    }
    const Manifest manifest{step,
                            time_s,
                            context.size(),
                            decomposition.global_extent(),
                            decomposition.periodic(),
                            decomposition.process_grid(),
                            checked_u32(persistent_fields(registry).size(),
                                        "persistent-field count"),
                            schema_fingerprint,
                            std::move(records)};
    manifest_bytes = encode_manifest(manifest);
    write_atomic_file(step_directory / "manifest.v1.bin", manifest_bytes);
  });

  converge_phase(context, "Restart v1 completion marker write", [&] {
    if (context.rank() != 0) {
      return;
    }
    if (injection.phase == RestartFailurePhase::marker && injection.rank == 0) {
      throw Error("injected Restart v1 marker failure");
    }
    const auto marker =
        encode_marker(static_cast<std::uint64_t>(manifest_bytes.size()),
                      crc64_ecma(manifest_bytes.data(), manifest_bytes.size()));
    write_atomic_file(step_directory / "COMPLETED", marker);
  });
  context.barrier();
}

} // namespace detail

void write_restart_checkpoint(const MpiContext &context,
                              const StructuredDecomposition &decomposition,
                              const FieldRegistry &registry,
                              const FieldStorage &storage,
                              const std::filesystem::path &step_directory,
                              std::int64_t step, double time_s) {
  detail::write_restart_checkpoint_with_failure(
      context, decomposition, registry, storage, step_directory, step, time_s,
      {});
}

RestartMetadata
read_restart_checkpoint(const MpiContext &context,
                        const StructuredDecomposition &decomposition,
                        const FieldRegistry &registry, FieldStorage &storage,
                        const std::filesystem::path &step_directory) {
  detail::require_mpi_active("read Restart v1 checkpoint");
  std::string path_text;
  std::int64_t path_step = 0;
  std::uint64_t schema_fingerprint = 0;
  detail::RestartCommitPlan commit_plan{};
  converge_phase(context, "Restart v1 read preflight", [&] {
    validate_communicators(context, decomposition);
    path_step = parse_step_leaf(step_directory);
    path_text = step_directory.string();
    if (path_text.empty()) {
      throw Error("Restart v1 step directory must not be empty");
    }
    schema_fingerprint = detail::restart_schema_fingerprint(registry);
    commit_plan = detail::make_restart_commit_plan(
        registry, storage, decomposition.local_extent());
  });
  require_path_agreement(context, path_text);
  require_metadata_agreement(context, decomposition, path_step, 0.0,
                             schema_fingerprint);
  const auto boxes = gather_boxes(context, decomposition.owned_box());
  converge_phase(context, "Restart v1 read owned-box validation", [&] {
    validate_complete_boxes(boxes, decomposition.global_extent());
  });

  std::vector<std::byte> marker_bytes;
  std::uint64_t expected_manifest_size = 0;
  std::uint64_t expected_manifest_crc = 0;
  converge_phase(context, "Restart v1 completion marker read", [&] {
    marker_bytes = read_file(step_directory / "COMPLETED");
    const auto marker = decode_marker(marker_bytes);
    expected_manifest_size = marker.first;
    expected_manifest_crc = marker.second;
  });

  std::vector<std::byte> manifest_bytes;
  Manifest manifest{};
  converge_phase(context, "Restart v1 manifest read", [&] {
    manifest_bytes = read_file(step_directory / "manifest.v1.bin");
    if (manifest_bytes.size() != expected_manifest_size ||
        detail::crc64_ecma(manifest_bytes.data(), manifest_bytes.size()) !=
            expected_manifest_crc) {
      throw Error("Restart v1 manifest size or CRC does not match marker");
    }
    manifest = decode_manifest(manifest_bytes);
    if (manifest.step != path_step || manifest.ranks != context.size() ||
        !same(manifest.global_extent, decomposition.global_extent()) ||
        manifest.periodic != decomposition.periodic() ||
        !same(manifest.process_grid, decomposition.process_grid()) ||
        manifest.schema_fingerprint != schema_fingerprint ||
        manifest.persistent_field_count != persistent_fields(registry).size() ||
        manifest.records.size() != boxes.size()) {
      throw Error("Restart v1 manifest does not match the current runtime");
    }
    for (std::size_t index = 0; index < boxes.size(); ++index) {
      if (!same(manifest.records[index].owned_box, boxes[index])) {
        throw Error("Restart v1 manifest owned partition does not match");
      }
    }
  });
  require_u64_agreement(context,
                        static_cast<std::uint64_t>(manifest_bytes.size()),
                        "Restart v1 manifest sizes differ across ranks");
  require_u64_agreement(context, expected_manifest_crc,
                        "Restart v1 manifest CRCs differ across ranks");
  require_u64_agreement(context, double_bits(manifest.time_s),
                        "Restart v1 manifest time differs across ranks");

  detail::StagedRestartRank staged{};
  converge_phase(context, "Restart v1 rank-file read", [&] {
    const RankRecord &record =
        manifest.records[static_cast<std::size_t>(context.rank())];
    const auto rank_bytes = read_file(step_directory / record.filename);
    if (rank_bytes.size() != record.byte_size ||
        detail::crc64_ecma(rank_bytes.data(), rank_bytes.size()) !=
            record.crc64) {
      throw Error("Restart v1 rank-file size or CRC does not match manifest");
    }
    const detail::RestartRankMetadata expected{context.rank(),
                                               context.size(),
                                               decomposition.global_extent(),
                                               decomposition.owned_box(),
                                               manifest.step,
                                               manifest.time_s};
    staged = detail::decode_restart_rank(rank_bytes, expected, registry);
  });

  detail::commit_restart_rank(staged, commit_plan);
  return RestartMetadata{manifest.step, manifest.time_s};
}

} // namespace hundun::runtime
