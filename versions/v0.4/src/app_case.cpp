// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_case.hpp"

#include "app_case_detail.hpp"
#include "yyjson.h"

#include <mpi.h>

#include <cerrno>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fcntl.h>
#include <initializer_list>
#include <limits>
#include <new>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <sys/types.h>
#include <type_traits>
#include <unordered_set>
#include <unistd.h>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

namespace fs = std::filesystem;

constexpr std::uint8_t kWireVersion = 1U;
constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr std::size_t kMaxJsonDepth = 32U;
constexpr std::string_view kSemanticContract =
    "HUNDUN-FLOW-v0.4-case-schema-v1|units=SI|"
    "flow=single_phase_low_mach_compressible|"
    "pressure_closure=local_absolute_pressure_drho_dp|reacting=false|"
    "coupling=PISO|pressure_correctors=2";

static_assert(std::is_nothrow_move_assignable_v<ValidatedModel>,
              "ValidatedModel publication must not allocate or throw");

enum Detail : std::uint32_t {
  detail_none = 0,
  detail_case_root = 1,
  detail_case_json = 2,
  detail_json_too_large = 3,
  detail_json_syntax = 4,
  detail_json_duplicate_key = 5,
  detail_json_schema = 6,
  detail_json_value = 7,
  detail_reference_count = 8,
  detail_reference_path = 9,
  detail_reference_missing = 10,
  detail_reference_too_large = 11,
  detail_wire = 12,
  detail_collective = 13,
};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<detail::FileOpenObserver> g_file_open_observer{nullptr};
std::atomic<int> g_last_lowest_failing_rank{-1};
#endif

Status invalid_case(Detail detail_code) noexcept {
  return {StatusCode::invalid_case, static_cast<std::uint32_t>(detail_code)};
}

void notify_file_open(int rank) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const auto observer = g_file_open_observer.load(std::memory_order_relaxed);
  if (observer != nullptr) {
    observer(rank);
  }
#else
  static_cast<void>(rank);
#endif
}

void record_lowest_failing_rank(int rank) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  g_last_lowest_failing_rank.store(rank, std::memory_order_relaxed);
#else
  static_cast<void>(rank);
#endif
}

class Hash64 {
 public:
  void bytes(const void* data, std::size_t size) noexcept {
    const auto* input = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
      value_ ^= static_cast<std::uint64_t>(input[index]);
      value_ *= kFnvPrime;
    }
  }

  void text(std::string_view value) noexcept {
    const std::uint64_t size = static_cast<std::uint64_t>(value.size());
    integer(size);
    bytes(value.data(), value.size());
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t byte = 0; byte < sizeof(Unsigned); ++byte) {
      const unsigned char part =
          static_cast<unsigned char>((bits >> (byte * 8U)) & 0xffU);
      bytes(&part, 1U);
    }
  }

  PlanFingerprint finish() const noexcept {
    return value_ == 0U ? 1U : value_;
  }

 private:
  std::uint64_t value_{kFnvOffset};
};

class Document {
 public:
  explicit Document(yyjson_doc* document) noexcept : document_(document) {}
  ~Document() { yyjson_doc_free(document_); }
  Document(const Document&) = delete;
  Document& operator=(const Document&) = delete;

  yyjson_doc* get() const noexcept { return document_; }

 private:
  yyjson_doc* document_{};
};

class UniqueFd {
 public:
  UniqueFd() noexcept = default;
  explicit UniqueFd(int descriptor) noexcept : descriptor_(descriptor) {}
  ~UniqueFd() {
    if (descriptor_ >= 0) {
      static_cast<void>(::close(descriptor_));
    }
  }
  UniqueFd(const UniqueFd&) = delete;
  UniqueFd& operator=(const UniqueFd&) = delete;
  UniqueFd(UniqueFd&& other) noexcept
      : descriptor_(std::exchange(other.descriptor_, -1)) {}
  UniqueFd& operator=(UniqueFd&& other) noexcept {
    if (this != &other) {
      if (descriptor_ >= 0) {
        static_cast<void>(::close(descriptor_));
      }
      descriptor_ = std::exchange(other.descriptor_, -1);
    }
    return *this;
  }

  explicit operator bool() const noexcept { return descriptor_ >= 0; }
  int get() const noexcept { return descriptor_; }

 private:
  int descriptor_{-1};
};

Status open_regular_at(int root_descriptor, const fs::path& relative,
                       int rank, Detail open_failure,
                       Detail non_regular_failure, UniqueFd& out,
                       struct stat& metadata) {
  notify_file_open(rank);
  const int descriptor = ::openat(
      root_descriptor, relative.c_str(),
      O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK | O_NOCTTY);
  if (descriptor < 0) {
    return invalid_case(errno == ELOOP ? non_regular_failure : open_failure);
  }
  UniqueFd candidate(descriptor);
  if (::fstat(candidate.get(), &metadata) != 0) {
    return invalid_case(open_failure);
  }
  if (!S_ISREG(metadata.st_mode)) {
    return invalid_case(non_regular_failure);
  }
  out = std::move(candidate);
  return {};
}

Status read_bounded_text(const UniqueFd& file, const struct stat& metadata,
                         std::uint64_t byte_limit, Detail read_failure,
                         Detail too_large, std::string& out) {
  if (metadata.st_size < 0 ||
      static_cast<std::uint64_t>(metadata.st_size) > byte_limit) {
    return invalid_case(too_large);
  }

  std::array<char, 64U * 1024U> buffer{};
  out.clear();
  out.reserve(static_cast<std::size_t>(metadata.st_size));
  std::uint64_t total = 0U;
  for (;;) {
    const ssize_t count = ::read(file.get(), buffer.data(), buffer.size());
    if (count < 0) {
      if (errno == EINTR) {
        continue;
      }
      return invalid_case(read_failure);
    }
    if (count == 0) {
      return {};
    }
    const auto unsigned_count = static_cast<std::uint64_t>(count);
    if (unsigned_count > byte_limit - total) {
      return invalid_case(too_large);
    }
    total += unsigned_count;
    out.append(buffer.data(), static_cast<std::size_t>(count));
  }
}

Status hash_bounded_file(const UniqueFd& file, const struct stat& metadata,
                         Hash64& hash) {
  if (metadata.st_size < 0 ||
      static_cast<std::uint64_t>(metadata.st_size) >
          detail::kMaxReferencedFileBytes) {
    return invalid_case(detail_reference_too_large);
  }

  std::array<char, 64U * 1024U> buffer{};
  std::uint64_t total = 0U;
  for (;;) {
    const ssize_t count = ::read(file.get(), buffer.data(), buffer.size());
    if (count < 0) {
      if (errno == EINTR) {
        continue;
      }
      return invalid_case(detail_reference_missing);
    }
    if (count == 0) {
      hash.integer(total);
      return {};
    }
    const auto unsigned_count = static_cast<std::uint64_t>(count);
    if (unsigned_count > detail::kMaxReferencedFileBytes - total) {
      return invalid_case(detail_reference_too_large);
    }
    total += unsigned_count;
    hash.bytes(buffer.data(), static_cast<std::size_t>(count));
  }
}

bool has_duplicate_keys(yyjson_val* value, std::size_t depth) {
  if (value == nullptr || depth > kMaxJsonDepth) {
    return true;
  }
  if (yyjson_is_obj(value)) {
    std::unordered_set<std::string_view> keys;
    const std::size_t object_size = yyjson_obj_size(value);
    keys.reserve(object_size);
    std::size_t index = 0U;
    std::size_t maximum = 0U;
    yyjson_val* key = nullptr;
    yyjson_val* child = nullptr;
    yyjson_obj_foreach(value, index, maximum, key, child) {
      const char* key_text = yyjson_get_str(key);
      if (key_text == nullptr) {
        return true;
      }
      const std::string_view view{key_text, yyjson_get_len(key)};
      if (!keys.insert(view).second || has_duplicate_keys(child, depth + 1U)) {
        return true;
      }
    }
    return false;
  }
  if (yyjson_is_arr(value)) {
    std::size_t index = 0U;
    std::size_t maximum = 0U;
    yyjson_val* child = nullptr;
    yyjson_arr_foreach(value, index, maximum, child) {
      if (has_duplicate_keys(child, depth + 1U)) {
        return true;
      }
    }
  }
  return false;
}

bool object_has_exact_keys(yyjson_val* object,
                           std::initializer_list<std::string_view> required) {
  if (!yyjson_is_obj(object) || yyjson_obj_size(object) != required.size()) {
    return false;
  }
  for (const std::string_view key : required) {
    if (yyjson_obj_getn(object, key.data(), key.size()) == nullptr) {
      return false;
    }
  }
  return true;
}

bool root_has_case_keys(yyjson_val* root) {
  if (!yyjson_is_obj(root)) {
    return false;
  }
  constexpr std::array<std::string_view, 6U> required{
      "schema_version", "units", "mesh", "flow", "solver", "time"};
  for (const std::string_view key : required) {
    if (yyjson_obj_getn(root, key.data(), key.size()) == nullptr) {
      return false;
    }
  }
  const bool has_turbulence = yyjson_obj_get(root, "turbulence") != nullptr;
  return yyjson_obj_size(root) == required.size() +
                                         (has_turbulence ? 1U : 0U);
}

std::optional<std::string_view> string_value(yyjson_val* object,
                                             std::string_view key) {
  yyjson_val* value = yyjson_obj_getn(object, key.data(), key.size());
  if (!yyjson_is_str(value)) {
    return std::nullopt;
  }
  return std::string_view{yyjson_get_str(value), yyjson_get_len(value)};
}

bool equals(std::optional<std::string_view> value,
            std::string_view expected) noexcept {
  return value.has_value() && *value == expected;
}

bool parse_geometry(std::string_view value, GeometryKind& out) noexcept {
  if (value == "uniform") {
    out = GeometryKind::uniform;
    return true;
  }
  if (value == "tensor_stretched") {
    out = GeometryKind::tensor_stretched;
    return true;
  }
  return false;
}

bool parse_turbulence(std::string_view value, TurbulenceKind& out) noexcept {
  if (value == "none") {
    out = TurbulenceKind::none;
    return true;
  }
  if (value == "wale") {
    out = TurbulenceKind::wale;
    return true;
  }
  if (value == "vreman_wall_function") {
    out = TurbulenceKind::vreman_wall_function;
    return true;
  }
  return false;
}

bool parse_time_control(std::string_view value,
                        TimeControlKind& out) noexcept {
  if (value == "fixed") {
    out = TimeControlKind::fixed;
    return true;
  }
  if (value == "adaptive_flow") {
    out = TimeControlKind::adaptive_flow;
    return true;
  }
  if (value == "adaptive_acoustic") {
    out = TimeControlKind::adaptive_acoustic;
    return true;
  }
  return false;
}

bool valid_direct_name(std::string_view text, std::string_view extension,
                       fs::path& out) {
  if (text.empty() || text.size() > detail::kMaxRelativePathBytes ||
      text.find('\\') != std::string_view::npos ||
      text.find('\0') != std::string_view::npos) {
    return false;
  }
  fs::path candidate{std::string(text)};
  if (candidate.is_absolute() || candidate.has_root_path() ||
      !candidate.parent_path().empty() || candidate.filename() != candidate ||
      candidate == "." || candidate == ".." ||
      candidate.extension().generic_string() != extension) {
    return false;
  }
  out = std::move(candidate);
  return true;
}

Status open_direct_file(int root_descriptor, std::string_view name,
                        std::string_view extension, int rank,
                        fs::path& relative, UniqueFd& descriptor,
                        struct stat& metadata) {
  if (!valid_direct_name(name, extension, relative)) {
    return invalid_case(detail_reference_path);
  }
  return open_regular_at(root_descriptor, relative, rank,
                         detail_reference_missing, detail_reference_path,
                         descriptor, metadata);
}

class WireWriter {
 public:
  void byte(std::uint8_t value) { bytes_.push_back(value); }

  void u16(std::uint16_t value) {
    byte(static_cast<std::uint8_t>(value & 0xffU));
    byte(static_cast<std::uint8_t>((value >> 8U) & 0xffU));
  }

  void u64(std::uint64_t value) {
    for (unsigned int shift = 0U; shift < 64U; shift += 8U) {
      byte(static_cast<std::uint8_t>((value >> shift) & 0xffU));
    }
  }

  bool text(std::string_view value) {
    if (value.size() > std::numeric_limits<std::uint16_t>::max()) {
      return false;
    }
    u16(static_cast<std::uint16_t>(value.size()));
    bytes_.insert(bytes_.end(), value.begin(), value.end());
    return bytes_.size() <= detail::kMaxWireBytes;
  }

  std::vector<std::uint8_t> take() && { return std::move(bytes_); }

 private:
  std::vector<std::uint8_t> bytes_;
};

class WireReader {
 public:
  explicit WireReader(const std::vector<std::uint8_t>& bytes) noexcept
      : bytes_(bytes) {}

  bool byte(std::uint8_t& out) noexcept {
    if (position_ >= bytes_.size()) {
      return false;
    }
    out = bytes_[position_++];
    return true;
  }

  bool u16(std::uint16_t& out) noexcept {
    std::uint8_t low = 0U;
    std::uint8_t high = 0U;
    if (!byte(low) || !byte(high)) {
      return false;
    }
    out = static_cast<std::uint16_t>(low) |
          static_cast<std::uint16_t>(static_cast<std::uint16_t>(high) << 8U);
    return true;
  }

  bool u64(std::uint64_t& out) noexcept {
    out = 0U;
    for (unsigned int shift = 0U; shift < 64U; shift += 8U) {
      std::uint8_t part = 0U;
      if (!byte(part)) {
        return false;
      }
      out |= static_cast<std::uint64_t>(part) << shift;
    }
    return true;
  }

  bool text(std::string& out) {
    std::uint16_t size = 0U;
    if (!u16(size) || size > detail::kMaxRelativePathBytes ||
        position_ + size > bytes_.size()) {
      return false;
    }
    out.assign(reinterpret_cast<const char*>(bytes_.data() + position_), size);
    position_ += size;
    return true;
  }

  bool finished() const noexcept { return position_ == bytes_.size(); }

 private:
  const std::vector<std::uint8_t>& bytes_;
  std::size_t position_{};
};

Status serialize_model(const ValidatedModel& model,
                       std::vector<std::uint8_t>& out) {
  try {
    if (model.data_files.size() > detail::kMaxReferencedFiles ||
        model.data_files.size() + (model.stl_file.has_value() ? 1U : 0U) >
            detail::kMaxReferencedFiles ||
        model.data_files.size() > std::numeric_limits<std::uint16_t>::max()) {
      return invalid_case(detail_wire);
    }
    WireWriter writer;
    writer.byte(kWireVersion);
    writer.byte(static_cast<std::uint8_t>(model.geometry));
    writer.byte(static_cast<std::uint8_t>(model.turbulence));
    writer.byte(static_cast<std::uint8_t>(model.time_control));
    writer.u16(static_cast<std::uint16_t>(model.data_files.size()));
    writer.byte(model.stl_file.has_value() ? 1U : 0U);
    for (const fs::path& path : model.data_files) {
      if (!writer.text(path.generic_string())) {
        return invalid_case(detail_wire);
      }
    }
    if (model.stl_file.has_value() &&
        !writer.text(model.stl_file->generic_string())) {
      return invalid_case(detail_wire);
    }
    writer.u64(model.fingerprint);
    out = std::move(writer).take();
    if (out.empty() || out.size() > detail::kMaxWireBytes) {
      return invalid_case(detail_wire);
    }
    return {};
  } catch (...) {
    return {StatusCode::allocation_failure, detail_wire};
  }
}

Status deserialize_model(const std::vector<std::uint8_t>& bytes,
                         ValidatedModel& out) {
  try {
    WireReader reader(bytes);
    std::uint8_t version = 0U;
    std::uint8_t geometry = 0U;
    std::uint8_t turbulence = 0U;
    std::uint8_t time_control = 0U;
    std::uint16_t data_count = 0U;
    std::uint8_t has_stl = 0U;
    if (!reader.byte(version) || version != kWireVersion ||
        !reader.byte(geometry) ||
        geometry > static_cast<std::uint8_t>(GeometryKind::tensor_stretched) ||
        !reader.byte(turbulence) ||
        turbulence >
            static_cast<std::uint8_t>(TurbulenceKind::vreman_wall_function) ||
        !reader.byte(time_control) ||
        time_control >
            static_cast<std::uint8_t>(TimeControlKind::adaptive_acoustic) ||
        !reader.u16(data_count) || data_count > detail::kMaxReferencedFiles ||
        !reader.byte(has_stl) || has_stl > 1U ||
        static_cast<std::size_t>(data_count) + (has_stl != 0U ? 1U : 0U) >
            detail::kMaxReferencedFiles) {
      return invalid_case(detail_wire);
    }

    ValidatedModel model;
    model.geometry = static_cast<GeometryKind>(geometry);
    model.turbulence = static_cast<TurbulenceKind>(turbulence);
    model.time_control = static_cast<TimeControlKind>(time_control);
    model.data_files.reserve(data_count);
    for (std::uint16_t index = 0U; index < data_count; ++index) {
      std::string path;
      if (!reader.text(path)) {
        return invalid_case(detail_wire);
      }
      fs::path parsed;
      if (!valid_direct_name(path, ".d", parsed)) {
        return invalid_case(detail_wire);
      }
      model.data_files.push_back(std::move(parsed));
    }
    if (has_stl != 0U) {
      std::string path;
      fs::path parsed;
      if (!reader.text(path) || !valid_direct_name(path, ".stl", parsed)) {
        return invalid_case(detail_wire);
      }
      model.stl_file = std::move(parsed);
    }
    if (!reader.u64(model.fingerprint) || model.fingerprint == 0U ||
        !reader.finished()) {
      return invalid_case(detail_wire);
    }
    out = std::move(model);
    return {};
  } catch (...) {
    return {StatusCode::allocation_failure, detail_wire};
  }
}

Status compile_on_root(const fs::path& case_root, int rank,
                       ValidatedModel& model,
                       std::vector<std::uint8_t>& payload) {
  try {
    std::error_code error;
    const fs::path canonical_root = fs::canonical(case_root, error);
    if (error) {
      return invalid_case(detail_case_root);
    }
    UniqueFd root_descriptor(::open(canonical_root.c_str(),
                                    O_RDONLY | O_CLOEXEC | O_DIRECTORY |
                                        O_NOFOLLOW | O_NONBLOCK | O_NOCTTY));
    struct stat root_metadata {};
    if (!root_descriptor ||
        ::fstat(root_descriptor.get(), &root_metadata) != 0 ||
        !S_ISDIR(root_metadata.st_mode)) {
      return invalid_case(detail_case_root);
    }

    UniqueFd json_descriptor;
    struct stat json_metadata {};
    const Status json_open =
        open_regular_at(root_descriptor.get(), fs::path{"case.json"}, rank,
                        detail_case_json, detail_case_json, json_descriptor,
                        json_metadata);
    if (!json_open) {
      return json_open;
    }
    std::string json;
    const Status json_read =
        read_bounded_text(json_descriptor, json_metadata,
                          detail::kMaxJsonBytes, detail_case_json,
                          detail_json_too_large, json);
    if (!json_read) {
      return json_read;
    }

    yyjson_read_err read_error{};
    Document document(yyjson_read_opts(json.data(), json.size(), 0U, nullptr,
                                       &read_error));
    if (document.get() == nullptr) {
      if (read_error.code == YYJSON_READ_ERROR_MEMORY_ALLOCATION) {
        return {StatusCode::allocation_failure, detail_json_syntax};
      }
      return invalid_case(detail_json_syntax);
    }
    yyjson_val* root = yyjson_doc_get_root(document.get());
    if (has_duplicate_keys(root, 0U)) {
      return invalid_case(detail_json_duplicate_key);
    }
    if (!root_has_case_keys(root)) {
      return invalid_case(detail_json_schema);
    }

    yyjson_val* schema_version = yyjson_obj_get(root, "schema_version");
    if (!yyjson_is_uint(schema_version) ||
        yyjson_get_uint(schema_version) != 1U ||
        !equals(string_value(root, "units"), "SI")) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* mesh = yyjson_obj_get(root, "mesh");
    yyjson_val* flow = yyjson_obj_get(root, "flow");
    yyjson_val* solver = yyjson_obj_get(root, "solver");
    yyjson_val* turbulence = yyjson_obj_get(root, "turbulence");
    yyjson_val* time = yyjson_obj_get(root, "time");
    if (!object_has_exact_keys(mesh, {"kind", "data_files", "stl_file"}) ||
        !object_has_exact_keys(flow,
                               {"model", "pressure_closure", "reacting"}) ||
        !object_has_exact_keys(solver,
                               {"coupling", "pressure_correctors"}) ||
        (turbulence != nullptr &&
         !object_has_exact_keys(turbulence, {"model"})) ||
        !object_has_exact_keys(time, {"control"})) {
      return invalid_case(detail_json_schema);
    }

    const auto geometry_text = string_value(mesh, "kind");
    const auto turbulence_text =
        turbulence == nullptr
            ? std::optional<std::string_view>{"vreman_wall_function"}
            : string_value(turbulence, "model");
    const auto time_text = string_value(time, "control");
    if (!geometry_text || !turbulence_text || !time_text ||
        !parse_geometry(*geometry_text, model.geometry) ||
        !parse_turbulence(*turbulence_text, model.turbulence) ||
        !parse_time_control(*time_text, model.time_control) ||
        !equals(string_value(flow, "model"),
                "single_phase_low_mach_compressible") ||
        !equals(string_value(flow, "pressure_closure"),
                "local_absolute_pressure_drho_dp") ||
        !equals(string_value(solver, "coupling"), "PISO")) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* reacting = yyjson_obj_get(flow, "reacting");
    yyjson_val* pressure_correctors =
        yyjson_obj_get(solver, "pressure_correctors");
    if (!yyjson_is_bool(reacting) || yyjson_get_bool(reacting) ||
        !yyjson_is_uint(pressure_correctors) ||
        yyjson_get_uint(pressure_correctors) != 2U) {
      return invalid_case(detail_json_value);
    }

    yyjson_val* data_files = yyjson_obj_get(mesh, "data_files");
    yyjson_val* stl_file = yyjson_obj_get(mesh, "stl_file");
    if (!yyjson_is_arr(data_files) ||
        (!yyjson_is_null(stl_file) && !yyjson_is_str(stl_file))) {
      return invalid_case(detail_reference_count);
    }
    const std::size_t data_file_count = yyjson_arr_size(data_files);
    const std::size_t total_reference_count =
        data_file_count + (yyjson_is_str(stl_file) ? 1U : 0U);
    if (data_file_count > detail::kMaxReferencedFiles ||
        total_reference_count > detail::kMaxReferencedFiles) {
      return invalid_case(detail_reference_count);
    }

    Hash64 hash;
    hash.text(kSemanticContract);
    hash.integer(static_cast<std::uint8_t>(model.geometry));
    hash.integer(static_cast<std::uint8_t>(model.turbulence));
    hash.integer(static_cast<std::uint8_t>(model.time_control));
    hash.integer(static_cast<std::uint16_t>(data_file_count));

    model.data_files.reserve(data_file_count);
    std::set<std::pair<dev_t, ino_t>> referenced_targets;
    std::size_t index = 0U;
    std::size_t maximum = 0U;
    yyjson_val* file_value = nullptr;
    yyjson_arr_foreach(data_files, index, maximum, file_value) {
      if (!yyjson_is_str(file_value)) {
        return invalid_case(detail_reference_path);
      }
      const std::string_view name{yyjson_get_str(file_value),
                                  yyjson_get_len(file_value)};
      fs::path relative;
      UniqueFd descriptor;
      struct stat metadata {};
      const Status opened = open_direct_file(
          root_descriptor.get(), name, ".d", rank, relative, descriptor,
          metadata);
      if (!opened) {
        return opened;
      }
      if (!referenced_targets
               .insert(std::make_pair(metadata.st_dev, metadata.st_ino))
               .second) {
        return invalid_case(detail_reference_path);
      }
      hash.text("data");
      hash.text(relative.generic_string());
      const Status hashed = hash_bounded_file(descriptor, metadata, hash);
      if (!hashed) {
        return hashed;
      }
      model.data_files.push_back(std::move(relative));
    }

    if (yyjson_is_str(stl_file)) {
      const std::string_view name{yyjson_get_str(stl_file),
                                  yyjson_get_len(stl_file)};
      fs::path relative;
      UniqueFd descriptor;
      struct stat metadata {};
      const Status opened = open_direct_file(
          root_descriptor.get(), name, ".stl", rank, relative, descriptor,
          metadata);
      if (!opened) {
        return opened;
      }
      if (!referenced_targets
               .insert(std::make_pair(metadata.st_dev, metadata.st_ino))
               .second) {
        return invalid_case(detail_reference_path);
      }
      hash.text("stl");
      hash.text(relative.generic_string());
      const Status hashed = hash_bounded_file(descriptor, metadata, hash);
      if (!hashed) {
        return hashed;
      }
      model.stl_file = std::move(relative);
    } else {
      hash.text("no-stl");
    }

    model.fingerprint = hash.finish();
    return serialize_model(model, payload);
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, detail_none};
  } catch (...) {
    return invalid_case(detail_case_root);
  }
}

Status bcast_header(std::array<std::uint64_t, 4U>& header, MPI_Comm communicator) {
  if (MPI_Bcast(header.data(), static_cast<int>(header.size()), MPI_UINT64_T, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }
  return {};
}

}  // namespace

namespace detail {

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void set_file_open_observer_for_test(FileOpenObserver observer) noexcept {
  g_file_open_observer.store(observer, std::memory_order_relaxed);
}

int last_lowest_failing_rank_for_test() noexcept {
  return g_last_lowest_failing_rank.load(std::memory_order_relaxed);
}
#endif

}  // namespace detail

Status CaseCompiler::load_and_compile(MPI_Comm communicator,
                                      const fs::path& case_root,
                                      ValidatedModel& out) {
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, detail_collective};
  }

  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, detail_collective};
  }

  Status root_status{};
  ValidatedModel root_model;
  std::vector<std::uint8_t> payload;
  if (rank == 0) {
    root_status = compile_on_root(case_root, rank, root_model, payload);
  }

  std::array<std::uint64_t, 4U> header{};
  if (rank == 0) {
    header[0] = static_cast<std::uint64_t>(root_status.code);
    header[1] = root_status.detail;
    header[2] = root_status ? static_cast<std::uint64_t>(payload.size()) : 0U;
    header[3] = root_status ? std::numeric_limits<std::uint64_t>::max() : 0U;
  }
  const Status header_status = bcast_header(header, communicator);
  if (!header_status) {
    return header_status;
  }
  if (header[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      header[1] > std::numeric_limits<std::uint32_t>::max()) {
    return {StatusCode::invalid_plan, detail_wire};
  }
  const Status collective_status{static_cast<StatusCode>(header[0]),
                                 static_cast<std::uint32_t>(header[1])};
  if (!collective_status) {
    record_lowest_failing_rank(
        header[3] <= static_cast<std::uint64_t>(std::numeric_limits<int>::max())
            ? static_cast<int>(header[3])
            : -1);
    return collective_status;
  }
  record_lowest_failing_rank(-1);
  if (header[2] == 0U || header[2] > detail::kMaxWireBytes ||
      header[2] > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    return {StatusCode::invalid_plan, detail_wire};
  }

  int allocation_ok = 1;
  if (rank != 0) {
    try {
      payload.resize(static_cast<std::size_t>(header[2]));
    } catch (...) {
      allocation_ok = 0;
    }
  }
  int all_allocations_ok = 0;
  if (MPI_Allreduce(&allocation_ok, &all_allocations_ok, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }
  if (all_allocations_ok == 0) {
    const int candidate = allocation_ok != 0 ? size : rank;
    int lowest = size;
    if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
        MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail_collective};
    }
    record_lowest_failing_rank(lowest);
    return {StatusCode::allocation_failure,
            static_cast<std::uint32_t>(lowest)};
  }

  if (MPI_Bcast(payload.data(), static_cast<int>(payload.size()), MPI_BYTE, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }

  ValidatedModel decoded;
  const Status decode_status = deserialize_model(payload, decoded);
  int decode_ok = decode_status ? 1 : 0;
  int all_decode_ok = 0;
  if (MPI_Allreduce(&decode_ok, &all_decode_ok, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, detail_collective};
  }
  if (all_decode_ok == 0) {
    const int local_allocation_failure =
        decode_status.code == StatusCode::allocation_failure ? 1 : 0;
    int any_allocation_failure = 0;
    if (MPI_Allreduce(&local_allocation_failure, &any_allocation_failure, 1,
                      MPI_INT, MPI_MAX, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail_collective};
    }
    const bool selected_failure =
        any_allocation_failure != 0 ? local_allocation_failure != 0
                                    : decode_ok == 0;
    const int candidate = selected_failure ? rank : size;
    int lowest = size;
    if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
        MPI_SUCCESS) {
      return {StatusCode::mpi_failure, detail_collective};
    }
    record_lowest_failing_rank(lowest);
    return any_allocation_failure != 0
               ? Status{StatusCode::allocation_failure,
                        static_cast<std::uint32_t>(lowest)}
               : Status{StatusCode::invalid_plan, detail_wire};
  }

  out = std::move(decoded);
  return {};
}

}  // namespace hundun::v04
